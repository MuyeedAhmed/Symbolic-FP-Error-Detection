import ast
import os
import textwrap
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple


class Node:
    """Node"""


@dataclass(frozen=True)
class Variable(Node):
    name: str
    width: str
    shape: Tuple[str, ...] = ()


@dataclass(frozen=True)
class Constant(Node):
    value: float
    width: str


@dataclass(frozen=True)
class Operation(Node):
    operation: str
    width: str
    strategy: str
    arguments: Tuple[Node, ...]
    detail: str = ""


def show(node: Node) -> str:
    if isinstance(node, Variable):
        return node.name
    if isinstance(node, Constant):
        return repr(node.value)
    head = node.operation + node.width + (f"^{node.strategy}" if node.strategy else "")
    if node.operation == "cmp":
        head = f"cmp{node.width}[{node.detail}]"
    return f"{head}({', '.join(show(argument) for argument in node.arguments)})"


def width_of(node: Node) -> str:
    return node.width


def leaves(node: Node, out: Optional[List[Node]] = None) -> List[Node]:
    out = [] if out is None else out
    if isinstance(node, (Variable, Constant)):
        if node not in out:
            out.append(node)
    else:
        for argument in node.arguments:
            leaves(argument, out)
    return out

@dataclass(frozen=True)
class Type:
    kind: str
    width: str = ""
    shape: Tuple[str, ...] = ()

    def __str__(self):
        if self.kind != "real":
            return self.kind
        return f"float{self.width}" + (f"[{','.join(self.shape)}]" if self.shape else "")


@dataclass
class Value:
    type: Type
    tree: Node
    
def float32(*shape) -> Type:
    return Type("real", "32", tuple(shape))


def float64(*shape) -> Type:
    return Type("real", "64", tuple(shape))


LIBRARY = {
    "np.sum": ("sum", "left"), "gpu_sum": ("sum", "tree"),
    "np.dot": ("dot", "left"), "np.matmul": ("dot", "left"), "gpu_gemm": ("dot", "tree"),
    "np.sqrt": ("sqrt", ""), "math.sqrt": ("sqrt", ""),
    "np.float32": ("cast", "32"), "np.float64": ("cast", "64"),
}
CASTS = {"np.float32": "32", "np.float64": "64"}
BINARY = {ast.Add: "+", ast.Sub: "-", ast.Mult: "x", ast.Div: "/", ast.MatMult: "dot"}
COMPARE = {ast.Lt: "<", ast.LtE: "<=", ast.Gt: ">", ast.GtE: ">=", ast.Eq: "==", ast.NotEq: "!="}


class NotInLanguage(Exception):
    pass


def _unify(left: Value, right: Value, where: str) -> str:
    lw, rw = left.type.width, right.type.width
    if lw == "any":
        return rw
    if rw == "any":
        return lw
    if lw == rw:
        return lw
    raise TypeError(f"mixed widths float{lw} and float{rw} in `{where}`: write the cast explicitly")


def _with_width(node: Node, width: str) -> Node:
    if isinstance(node, Constant) and node.width == "any":
        return Constant(node.value, width)
    return node


def _broadcast(left: Type, right: Type) -> Tuple[str, ...]:
    return left.shape if len(left.shape) >= len(right.shape) else right.shape


def type_expression(expression: ast.AST, environment: Dict[str, Value], line: int) -> Value:
    where = f"{ast.unparse(expression)} at L{line}"
    if isinstance(expression, ast.Constant):
        if isinstance(expression.value, bool) or not isinstance(expression.value, (int, float)):
            raise NotInLanguage(f"literal `{expression.value!r}` at L{line}")
        return Value(Type("real", "any"), Constant(float(expression.value), "any"))
    if isinstance(expression, ast.Name):
        if expression.id not in environment:
            raise TypeError(f"`{expression.id}` at L{line} is not an input and was not assigned")
        return environment[expression.id]
    if isinstance(expression, ast.Attribute) and expression.attr == "T":
        base = type_expression(expression.value, environment, line)
        transposed = Type("real", base.type.width, tuple(reversed(base.type.shape)))
        return Value(transposed, Operation("transpose", base.type.width, "", (base.tree,)))
    if isinstance(expression, ast.UnaryOp) and isinstance(expression.op, ast.USub):
        operand = type_expression(expression.operand, environment, line)
        if isinstance(operand.tree, Constant):
            return Value(operand.type, Constant(-operand.tree.value, operand.tree.width))
        zero = Constant(0.0, operand.type.width)
        return Value(operand.type, Operation("-", operand.type.width, "", (zero, operand.tree)))
    if isinstance(expression, ast.BinOp):
        if type(expression.op) not in BINARY:
            raise NotInLanguage(f"operator `{ast.unparse(expression)}` at L{line}")
        left = type_expression(expression.left, environment, line)
        right = type_expression(expression.right, environment, line)
        width = _unify(left, right, where)
        operation = BINARY[type(expression.op)]
        if operation == "dot":
            return _dot(left, right, width, "left")
        shape = _broadcast(left.type, right.type)
        node = Operation(operation, width, "", (_with_width(left.tree, width), _with_width(right.tree, width)))
        return Value(Type("real", width, shape), node)
    if isinstance(expression, ast.Compare):
        if len(expression.ops) != 1:
            raise NotInLanguage(f"chained comparison at L{line}")
        left = type_expression(expression.left, environment, line)
        right = type_expression(expression.comparators[0], environment, line)
        width = _unify(left, right, where)
        operator = COMPARE[type(expression.ops[0])]
        node = Operation("cmp", width, "", (_with_width(left.tree, width), _with_width(right.tree, width)), operator)
        return Value(Type("bool"), node)
    if isinstance(expression, ast.Call):
        name = ast.unparse(expression.func)
        if isinstance(expression.func, ast.Attribute) and expression.func.attr == "astype":
            base = type_expression(expression.func.value, environment, line)
            width = CASTS[ast.unparse(expression.args[0])]
            return Value(Type("real", width, base.type.shape), Operation("cast", width, "", (base.tree,)))
        if name not in LIBRARY:
            raise NotInLanguage(f"call `{name}` at L{line}")
        operation, strategy = LIBRARY[name]
        arguments = [type_expression(argument, environment, line) for argument in expression.args]
        if operation == "cast":
            base = arguments[0]
            return Value(Type("real", strategy, base.type.shape), Operation("cast", strategy, "", (base.tree,)))
        if operation == "sum":
            base = arguments[0]
            return Value(Type("real", base.type.width), Operation("sum", base.type.width, strategy, (base.tree,)))
        if operation == "sqrt":
            base = arguments[0]
            return Value(base.type, Operation("sqrt", base.type.width, "", (base.tree,)))
        if operation == "dot":
            left, right = arguments
            return _dot(left, right, _unify(left, right, where), strategy)
    raise NotInLanguage(f"`{ast.unparse(expression)}` at L{line}")


def _dot(left: Value, right: Value, width: str, strategy: str) -> Value:
    ls, rs = left.type.shape, right.type.shape
    shape = (ls[:-1] + rs[1:]) if len(ls) >= 2 and len(rs) >= 2 else (ls[:1] if len(ls) == 2 else ())
    return Value(Type("real", width, shape), Operation("dot", width, strategy, (left.tree, right.tree)))


@dataclass
class Statement:
    line: int
    name: str
    source: str
    value: Value


@dataclass
class Analysis:
    label: str
    inputs: Dict[str, Type]
    statements: List[Statement]
    output_name: str
    output: Value


def _key(node: Node):
    if isinstance(node, Variable):
        return ("var", node.name)
    if isinstance(node, Constant):
        return ("const", node.value, node.width)
    return ("op", node.operation, node.width, node.strategy, node.detail, len(node.arguments))


def first_difference(a: Node, b: Node) -> Optional[Tuple[Node, Node]]:
    if _key(a) != _key(b):
        return a, b
    if isinstance(a, Operation):
        for argument_a, argument_b in zip(a.arguments, b.arguments):
            found = first_difference(argument_a, argument_b)
            if found:
                return found
    return None


def analyse(source: str, inputs: Dict[str, Type], label: str) -> Analysis:
    source = textwrap.dedent(source).strip("\n")
    environment: Dict[str, Value] = {
        name: Value(input_type, Variable(name, input_type.width, input_type.shape)) for name, input_type in inputs.items()}
    versions: Dict[str, int] = {}
    statements: List[Statement] = []
    for statement in ast.parse(source).body:
        if not isinstance(statement, ast.Assign) or not isinstance(statement.targets[0], ast.Name):
            raise NotInLanguage(f"statement `{ast.unparse(statement).splitlines()[0]}` at L{statement.lineno}")
        name = statement.targets[0].id
        value = type_expression(statement.value, environment, statement.lineno)
        environment[name] = value
        versions[name] = versions.get(name, 0) + 1
        shown_name = name if versions[name] == 1 else f"{name}#{versions[name]}"
        statements.append(Statement(statement.lineno, shown_name, ast.unparse(statement.value), value))
    last = statements[-1]
    return Analysis(label, inputs, statements, last.name.split("#")[0], last.value)

def print_analysis(analysis: Analysis):
    print(f"  {analysis.label}")
    print(f"    inputs: " + ", ".join(f"{name} : {input_type}" for name, input_type in analysis.inputs.items()))
    for statement in analysis.statements:
        print(f"    L{statement.line}  {statement.name} = {show(statement.value.tree)}")
    print(f"    output {analysis.output_name} : {analysis.output.type}")


def print_z3_sketch(a: Analysis, b: Analysis, difference: Optional[Tuple[Node, Node]]):
    print("Z3 Encoding:")
    declared = []
    for leaf in leaves(a.output.tree) + leaves(b.output.tree):
        if isinstance(leaf, Variable) and leaf.name not in declared:
            declared.append(leaf.name)
            print(f"    declare  : {leaf.name} : {a.inputs[leaf.name]}")
        elif isinstance(leaf, Constant) and repr(leaf) not in declared:
            declared.append(repr(leaf))
            print(f"    constant : {leaf.value!r} as float{leaf.width}")
    print(f"    A        : {a.output_name}_A = {show(a.output.tree)}")
    print(f"    B        : {b.output_name}_B = {show(b.output.tree)}")
    print(f"    assert   : {a.output_name}_A != {b.output_name}_B")

    used = {node.operation for node in _all_nodes(a.output.tree) + _all_nodes(b.output.tree) if isinstance(node, Operation)}



def _all_nodes(node: Node, out: Optional[List[Node]] = None) -> List[Node]:
    out = [] if out is None else out
    if node not in out:
        out.append(node)
        if isinstance(node, Operation):
            for argument in node.arguments:
                _all_nodes(argument, out)
    return out


def to_dot(a: Analysis, b: Analysis, difference: Optional[Tuple[Node, Node]], title: str) -> str:
    lines = ["digraph effects {", f'  label="{title}"; labelloc=t; fontsize=14;',
             "  rankdir=BT; node [shape=box, fontname=Helvetica, fontsize=11];"]
    def cluster(analysis: Analysis, tag: str, red: Optional[Node]):
        ids: Dict[Node, str] = {}
        lines.append(f'  subgraph cluster_{tag} {{ label="{analysis.label}"; style=rounded; color=gray50;')
        for index, node in enumerate(_all_nodes(analysis.output.tree)):
            ids[node] = f"{tag}{index}"
            if isinstance(node, Variable):
                text, style = f"{node.name} : {analysis.inputs[node.name]}", "shape=ellipse"
            elif isinstance(node, Constant):
                text, style = f"{node.value!r} (float{node.width})", "shape=plaintext"
            else:
                head = node.operation + node.width + (f"^{node.strategy}" if node.strategy else "")
                text = f"cmp{node.width} {node.detail}" if node.operation == "cmp" else head
                style = "style=filled, fillcolor=gray92"
            if red is not None and node == red:
                style += ", color=red, penwidth=2.5, fontcolor=red"
            if node is analysis.output.tree:
                style += ", peripheries=2"
            lines.append(f'    {ids[node]} [label="{text}", {style}];')
        for node in ids:
            if isinstance(node, Operation):
                for position, argument in enumerate(node.arguments):
                    edge_label = f' [label="{position}"]' if len(node.arguments) > 1 and node.operation in ("-", "/", "cmp", "dot") else ""
                    lines.append(f"    {ids[argument]} -> {ids[node]}{edge_label};")
        lines.append("  }")

    cluster(a, "A", difference[0] if difference else None)
    cluster(b, "B", difference[1] if difference else None)
    lines.append("}")
    return "\n".join(lines) + "\n"


def compare(title: str, file_name: str, source_a: str, source_b: str, inputs: Dict[str, Type],
            label_a: str = "A (sklearn-style)", label_b: str = "B (cuML-style)", dot_directory: str = "dot"):
    print("=" * 100)
    print(f"  {title}")
    print("=" * 100)
    a, b = analyse(source_a, inputs, label_a), analyse(source_b, inputs, label_b)
    print_analysis(a)
    print_analysis(b)
    difference = first_difference(a.output.tree, b.output.tree)
    print("  Comparison of the two output trees:")
    if difference is None:
        print("    Identical")
    else:
        print("    Different")
    print_z3_sketch(a, b, difference)
    os.makedirs(dot_directory, exist_ok=True)
    path = os.path.join(dot_directory, f"{file_name}.dot")
    with open(path, "w") as handle:
        handle.write(to_dot(a, b, difference, title))
    
    return a, b, difference

input = dict(title="(A @ B) @ C vs A @ (B @ C)", file_name="Graph",
         inputs={"A": float32("m", "k"), "B": float32("k", "l"), "C": float32("l", "n")},
         source_a="Z = (A @ B) @ C",
         source_b="Z = A @ (B @ C)")


if __name__ == "__main__":
    here = os.path.dirname(os.path.abspath(__file__))
    results = compare(**input, dot_directory=os.path.join(here, "dot"))