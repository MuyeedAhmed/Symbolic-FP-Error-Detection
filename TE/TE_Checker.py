import ast
import copy
import textwrap
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Set, Tuple


class Behaviour:
    """Base"""

@dataclass(frozen=True)
class Empty(Behaviour):
    pass

@dataclass(frozen=True)
class Variable(Behaviour):
    name: str

@dataclass(frozen=True)
class Atom(Behaviour):
    operation: str
    precision: str = ""
    strategy: str = ""
    operands: str = ""
    source_position: str = field(default="", compare=False, hash=False)


@dataclass(frozen=True)
class Sequence(Behaviour):
    parts: Tuple[Behaviour, ...]


@dataclass(frozen=True)
class Choice(Behaviour):
    parts: Tuple[Behaviour, ...]


@dataclass(frozen=True)
class Recursion(Behaviour):
    variable: str
    body: Behaviour


EMPTY = Empty()
MARKERS = {"in"}


def show_atom(atom: Atom, labels: bool = True) -> str:
    text = atom.operation + atom.precision
    if atom.strategy:
        text += f"^{atom.strategy}"
    if labels and atom.operands:
        text += f"[{atom.operands}]"
    return text


def show(behaviour: Behaviour, labels: bool = True) -> str:
    if isinstance(behaviour, Empty):
        return "empty"
    if isinstance(behaviour, Variable):
        return behaviour.name
    if isinstance(behaviour, Atom):
        return show_atom(behaviour, labels)
    if isinstance(behaviour, Sequence):
        return " ; ".join(f"({show(part, labels)})" if isinstance(part, Choice) else show(part, labels)
                          for part in behaviour.parts)
    if isinstance(behaviour, Choice):
        return " | ".join(f"({show(part, labels)})" if isinstance(part, Sequence) else show(part, labels)
                          for part in behaviour.parts)
    if isinstance(behaviour, Recursion):
        return f"rec {behaviour.variable}.({show(behaviour.body, labels)})"
    raise TypeError(behaviour)


def mentions(behaviour: Behaviour, variable: str) -> bool:
    if isinstance(behaviour, Variable):
        return behaviour.name == variable
    if isinstance(behaviour, (Sequence, Choice)):
        return any(mentions(part, variable) for part in behaviour.parts)
    if isinstance(behaviour, Recursion):
        return behaviour.variable != variable and mentions(behaviour.body, variable)
    return False


def substitute(behaviour: Behaviour, variable: str, replacement: Behaviour) -> Behaviour:
    if isinstance(behaviour, Variable):
        return replacement if behaviour.name == variable else behaviour
    if isinstance(behaviour, Sequence):
        return Sequence(tuple(substitute(part, variable, replacement) for part in behaviour.parts))
    if isinstance(behaviour, Choice):
        return Choice(tuple(substitute(part, variable, replacement) for part in behaviour.parts))
    if isinstance(behaviour, Recursion):
        if behaviour.variable == variable:
            return behaviour
        return Recursion(behaviour.variable, substitute(behaviour.body, variable, replacement))
    return behaviour


def _fold_accumulator(variable: str, body: Behaviour) -> Optional[Behaviour]:
    arms = list(body.parts) if isinstance(body, Choice) else [body]
    increment, other_arms = None, []
    for arm in arms:
        last = arm.parts[-1] if isinstance(arm, Sequence) else None
        if (increment is None and isinstance(arm, Sequence) and len(arm.parts) >= 2
                and arm.parts[0] == Variable(variable)
                and isinstance(last, Atom) and last.operation == "+" and not last.strategy
                and not any(mentions(part, variable) for part in arm.parts[1:])):
            increment = arm
        else:
            other_arms.append(arm)
    if increment is None or any(mentions(arm, variable) for arm in other_arms):
        return None
    plus = increment.parts[-1]
    left_fold = Atom("+", plus.precision, "left", plus.operands, plus.source_position)
    return sequence(choice(*other_arms) if other_arms else EMPTY, *increment.parts[1:-1], left_fold)


def normalise(behaviour: Behaviour) -> Behaviour:
    if isinstance(behaviour, Sequence):
        parts: List[Behaviour] = []
        for part in behaviour.parts:
            normalised = normalise(part)
            if isinstance(normalised, Empty):
                continue
            parts.extend(normalised.parts if isinstance(normalised, Sequence) else [normalised])
        if not parts:
            return EMPTY
        return parts[0] if len(parts) == 1 else Sequence(tuple(parts))
    if isinstance(behaviour, Choice):
        flat: List[Behaviour] = []
        for part in behaviour.parts:
            normalised = normalise(part)
            flat.extend(normalised.parts if isinstance(normalised, Choice) else [normalised])
        unique: List[Behaviour] = []
        for part in flat:
            if part not in unique:
                unique.append(part)
        unique.sort(key=lambda part: (not isinstance(part, Empty), show(part)))
        return unique[0] if len(unique) == 1 else Choice(tuple(unique))
    if isinstance(behaviour, Recursion):
        variable = behaviour.variable
        body = normalise(behaviour.body)
        if isinstance(body, Choice) and Variable(variable) in body.parts:
            body = normalise(Choice(tuple(arm for arm in body.parts if arm != Variable(variable))))
        if body == Variable(variable):
            return EMPTY
        if not mentions(body, variable):
            return body
        folded = _fold_accumulator(variable, body)
        return folded if folded is not None else Recursion(variable, body)
    return behaviour


def sequence(*parts: Behaviour) -> Behaviour:
    return normalise(Sequence(tuple(parts)))


def choice(*parts: Behaviour) -> Behaviour:
    return normalise(Choice(tuple(parts)))


def recursion(variable: str, body: Behaviour) -> Behaviour:
    return normalise(Recursion(variable, body))

def atoms_in_order(behaviour: Behaviour) -> List[Atom]:
    atoms: List[Atom] = []

    def walk(term: Behaviour):
        if isinstance(term, Atom):
            atoms.append(term)
        elif isinstance(term, (Sequence, Choice)):
            for part in term.parts:
                walk(part)
        elif isinstance(term, Recursion):
            walk(term.body)
    walk(behaviour)
    return atoms

def shown(atoms: List[Atom]) -> List[str]:
    texts: List[str] = []
    for atom in atoms:
        if atom.operation in MARKERS:
            continue
        text = show_atom(atom, False)
        if text not in texts:
            texts.append(text)
    return texts

def compact(behaviour: Behaviour) -> str:
    return "{" + " ".join(shown(atoms_in_order(behaviour))) + "}"


def is_float_atom(atom: Atom) -> bool:
    return atom.precision in ("32", "64")


@dataclass(frozen=True)
class Type:
    precision: str = ""
    kind: str = "real"
    shape: Tuple = ()

    def __str__(self):
        if self.kind != "real":
            return self.kind
        precision = "?" if self.precision == "any" else self.precision
        return f"float{precision}" + (f"[{','.join(map(str, self.shape))}]" if self.shape else "")


INT, BOOL, STR, FUNCTION = Type("", "int"), Type("", "bool"), Type("", "str"), Type("", "fn")


def real32(*shape) -> Type:
    return Type("32", "real", shape)


def real64(*shape) -> Type:
    return Type("64", "real", shape)


@dataclass
class Value:
    type: Type
    provenance: Behaviour = EMPTY
    elements: Tuple["Value", ...] = ()
    function: Optional[ast.FunctionDef] = None


Environment = Dict[str, Value]


def unify_precision(left: str, right: str, context: str) -> str:
    if left == "any":
        return right
    if right == "any":
        return left
    if left == right:
        return left
    if "|" in left or "|" in right:
        return "32|64"
    raise TypeError(f"both operands must have one width, got float{left} and float{right} in `{context}`: "
                    f"write the cast explicitly (np.float32/np.float64/.astype) so it appears as an effect")


def join_type(left: Type, right: Type) -> Type:
    if left == right or right.kind != "real" and left.kind == "real" and left.precision != "any":
        return left
    if left.kind != "real" or left.precision == "any":
        return right
    if right.kind != "real" or right.precision == "any":
        return left
    return Type("32|64", "real", left.shape or right.shape)


def join_value(left: Optional[Value], right: Optional[Value]) -> Value:
    if left is None:
        return right
    if right is None:
        return left
    return Value(join_type(left.type, right.type), choice(left.provenance, right.provenance),
                 left.elements or right.elements, left.function or right.function)


def merge_environments(first: Environment, second: Environment) -> Environment:
    names = list(first) + [name for name in second if name not in first]
    return {name: join_value(first.get(name), second.get(name)) for name in names}



LIBRARY = {
    "np.dot": ("matmul", "left"), "np.matmul": ("matmul", "left"), ".dot": ("matmul", "left"),
    "np.einsum": ("matmul", "left"), "np.sum": ("reduce", "left"),
    "np.sqrt": ("sqrt", ""), "math.sqrt": ("sqrt", ""), "sqrt": ("sqrt", ""),
    "np.argmin": ("argmin", ""), "np.argmax": ("argmax", ""),
    "np.argsort": ("argsort", ""), "np.argpartition": ("argsort", ""),
    "np.maximum": ("same", ""), "np.minimum": ("same", ""), "np.abs": ("same", ""), "abs": ("same", ""),
    "np.zeros": ("fresh_real", ""), "np.zeros_like": ("like", ""), "np.full": ("fresh_int", ""),
    "np.arange": ("fresh_int", ""), "np.bincount": ("count", ""),
    "np.float64": ("cast", "64"), "np.float32": ("cast", "32"),
    "len": ("shape", ""), "range": ("shape", ""), "int": ("int", ""), "min": ("same", ""), "max": ("same", ""),
    "gpu_rowsum": ("reduce", "tree"), "gpu_sum": ("reduce", "tree"),
    "gpu_gemm": ("matmul", "tree"), "gpu_reduce_rows_by_key": ("reduce", "tree"),
    "gpu_argmin_rows": ("argmin", ""), "gpu_select_k": ("argsort", ""), "gpu_bincount": ("count", ""),
}
DTYPES = {"np.float64": "64", "np.float32": "32", "float64": "64", "float32": "32"}
BINARY_OPERATORS = {ast.Mult: "x", ast.MatMult: "x", ast.Add: "+", ast.Sub: "-", ast.Div: "/", ast.Pow: "x"}


class _StripLabel(ast.NodeTransformer):
    """Operand labels ignore indexing and casts: x[i][k], np.float64(x), x.astype(..) all name x."""

    def visit_Subscript(self, node):
        return self.visit(node.value)

    def visit_Call(self, node):
        function_name = ast.unparse(node.func)
        if function_name in DTYPES or function_name.endswith(".astype"):
            base = node.func.value if function_name.endswith(".astype") else node.args[0]
            return self.visit(base)
        return self.generic_visit(node)


def label(expression: ast.AST) -> str:
    return ast.unparse(_StripLabel().visit(copy.deepcopy(expression))).replace(" ", "")


@dataclass
class Decision:
    kind: str
    atom: Atom
    slice: Behaviour
    source_position: str
    source: str

    @property
    def is_float(self) -> bool:
        return is_float_atom(self.atom)

    def inputs(self) -> List[str]:
        return sorted({atom.operands for atom in atoms_in_order(self.slice) if atom.operation == "in"})


@dataclass
class TypedStatement:
    kind: str
    source: str
    behaviour: Behaviour
    line: int
    source_position: str
    effect_note: str = ""
    then_branch: List["TypedStatement"] = field(default_factory=list)
    else_branch: List["TypedStatement"] = field(default_factory=list)
    body: List["TypedStatement"] = field(default_factory=list)
    inlined: List[Tuple[str, List["TypedStatement"]]] = field(default_factory=list)
    is_decision: bool = False


class Context:
    def __init__(self):
        self.decisions: List[Decision] = []
        self.control: List[Behaviour] = []
        self.function_stack: List[str] = []
        self.returns: List[List[Value]] = []
        self.inlined: List[List[Tuple[str, List[TypedStatement]]]] = [[]]
        self.depth = 0

    def source_position(self, node) -> str:
        position = f"L{getattr(node, 'lineno', '?')}"
        return position + (f" in {self.function_stack[-1]}" if self.function_stack else "")


def _target_names(target: ast.AST) -> Set[str]:
    if isinstance(target, ast.Name):
        return {target.id}
    if isinstance(target, ast.Subscript):
        return _target_names(target.value)
    if isinstance(target, ast.Tuple):
        return set().union(*(_target_names(element) for element in target.elts))
    return set()


def _assigned_names(statements: List[ast.stmt]) -> Set[str]:
    names: Set[str] = set()
    for statement in statements:
        for node in ast.walk(statement):
            if isinstance(node, ast.Assign):
                for target in node.targets:
                    names |= _target_names(target)
            elif isinstance(node, ast.AugAssign):
                names |= _target_names(node.target)
            elif (isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute)
                  and node.func.attr in ("append", "pop", "extend")):
                names |= _target_names(node.func.value)
    return names


def _element(value: Value) -> Value:
    return Value(Type(value.type.precision, value.type.kind, value.type.shape[1:]), value.provenance)


def _real_precision(value: Value) -> str:
    return value.type.precision if value.type.kind == "real" else "any"


def type_expression(expression: ast.AST, environment: Environment, context: Context) -> Tuple[Value, Behaviour]:
    if isinstance(expression, ast.Constant):
        if isinstance(expression.value, bool):
            return Value(BOOL), EMPTY
        if isinstance(expression.value, int):
            return Value(INT), EMPTY
        if isinstance(expression.value, float):
            return Value(Type("any")), EMPTY
        if isinstance(expression.value, str):
            return Value(STR), EMPTY
        return Value(INT), EMPTY
    if isinstance(expression, ast.Name):
        if expression.id not in environment:
            raise TypeError(f"free variable `{expression.id}` at {context.source_position(expression)}: declare it in `inputs`")
        return environment[expression.id], EMPTY
    if isinstance(expression, ast.Attribute):
        dotted = ast.unparse(expression)
        if dotted in environment:
            return environment[dotted], EMPTY
        if dotted in DTYPES:
            return Value(STR), EMPTY
        base, behaviour = type_expression(expression.value, environment, context)
        if expression.attr == "T":
            transposed = Type(base.type.precision, base.type.kind, tuple(reversed(base.type.shape)))
            return Value(transposed, base.provenance), behaviour
        if expression.attr in ("shape", "size", "ndim"):
            return Value(INT), behaviour
        raise TypeError(f"unsupported attribute `{dotted}` at {context.source_position(expression)}")
    if isinstance(expression, ast.Subscript):
        base, behaviour = type_expression(expression.value, environment, context)
        index_parts = expression.slice.elts if isinstance(expression.slice, ast.Tuple) else [expression.slice]
        index_provenances, dropped_dimensions = [], 0
        for index in index_parts:
            if isinstance(index, ast.Slice) or (isinstance(index, ast.Constant) and index.value is None):
                continue
            index_value, index_behaviour = type_expression(index, environment, context)
            behaviour = sequence(behaviour, index_behaviour)
            index_provenances.append(index_value.provenance)
            dropped_dimensions += 1
        result_type = Type(base.type.precision, base.type.kind, base.type.shape[dropped_dimensions:])
        return Value(result_type, sequence(base.provenance, *index_provenances)), behaviour
    if isinstance(expression, ast.UnaryOp):
        value, behaviour = type_expression(expression.operand, environment, context)
        return (Value(BOOL, value.provenance) if isinstance(expression.op, ast.Not) else value), behaviour
    if isinstance(expression, ast.BoolOp):
        values, behaviours = zip(*(type_expression(operand, environment, context) for operand in expression.values))
        return Value(BOOL, sequence(*(value.provenance for value in values))), sequence(*behaviours)
    if isinstance(expression, ast.IfExp):
        test, test_behaviour = type_expression(expression.test, environment, context)
        then, then_behaviour = type_expression(expression.body, environment, context)
        other, other_behaviour = type_expression(expression.orelse, environment, context)
        return (Value(join_type(then.type, other.type), sequence(test.provenance, choice(then.provenance, other.provenance))),
                sequence(test_behaviour, choice(then_behaviour, other_behaviour)))
    if isinstance(expression, (ast.Tuple, ast.List)):
        values, behaviours = [], []
        for element in expression.elts:
            value, behaviour = type_expression(element, environment, context)
            values.append(value)
            behaviours.append(behaviour)
        kind = "tuple" if isinstance(expression, ast.Tuple) else "int"
        return Value(Type("", kind), sequence(*(value.provenance for value in values)), tuple(values)), sequence(*behaviours)
    if isinstance(expression, ast.BinOp):
        left, left_behaviour = type_expression(expression.left, environment, context)
        right, right_behaviour = type_expression(expression.right, environment, context)
        operands = label(expression)
        position = context.source_position(expression)
        if left.type.kind != "real" and right.type.kind != "real":            # exact integer arithmetic
            return Value(INT, sequence(left.provenance, right.provenance)), sequence(left_behaviour, right_behaviour)
        precision = unify_precision(_real_precision(left), _real_precision(right), f"{operands} at {position}")
        shape = left.type.shape if len(left.type.shape) >= len(right.type.shape) else right.type.shape
        if isinstance(expression.op, ast.MatMult):
            shape = (left.type.shape[:1] + right.type.shape[-1:]) if left.type.shape and right.type.shape else ()
            atoms = [Atom("x", precision, "", operands, position), Atom("+", precision, "left", operands, position)]
        else:
            atoms = [Atom(BINARY_OPERATORS[type(expression.op)], precision, "", operands, position)]
        return (Value(Type(precision, "real", shape), sequence(left.provenance, right.provenance, *atoms)),
                sequence(left_behaviour, right_behaviour, *atoms))
    if isinstance(expression, ast.Compare):
        left, behaviour = type_expression(expression.left, environment, context)
        provenance = left.provenance
        position = context.source_position(expression)
        for right_expression in expression.comparators:
            right, right_behaviour = type_expression(right_expression, environment, context)
            behaviour = sequence(behaviour, right_behaviour)
            if left.type.kind == "real" or right.type.kind == "real":
                precision = unify_precision(_real_precision(left), _real_precision(right), f"{label(expression)} at {position}")
            else:
                precision = "int"
            atom = Atom("cmp", precision, "", label(expression), position)
            provenance = sequence(provenance, right.provenance, atom)
            behaviour = sequence(behaviour, atom)
            context.decisions.append(Decision("cmp", atom, sequence(*context.control, provenance), position, label(expression)))
            left = right
        return Value(BOOL, provenance), behaviour
    if isinstance(expression, ast.Call):
        return type_call(expression, environment, context)
    raise TypeError(f"unsupported expression `{ast.unparse(expression)}` at {context.source_position(expression)}")


def type_call(call: ast.Call, environment: Environment, context: Context) -> Tuple[Value, Behaviour]:
    function_name = ast.unparse(call.func)
    position = context.source_position(call)
    if isinstance(call.func, ast.Attribute) and function_name not in LIBRARY and function_name not in environment:
        receiver_name = ast.unparse(call.func.value)                                    # method call
        if call.func.attr == "astype":
            receiver, behaviour = type_expression(call.func.value, environment, context)
            precision = DTYPES[ast.unparse(call.args[0]).replace("'", "").replace('"', "")]
            atom = Atom("cast", precision, "", label(call.func.value), position)
            return Value(Type(precision, "real", receiver.type.shape), sequence(receiver.provenance, atom)), sequence(behaviour, atom)
        if call.func.attr in ("append", "extend"):
            argument, behaviour = type_expression(call.args[0], environment, context)
            current = environment[receiver_name]
            environment[receiver_name] = Value(current.type, choice(current.provenance, sequence(*context.control, argument.provenance)))
            return Value(INT), behaviour
        if call.func.attr == "pop":
            receiver, behaviour = type_expression(call.func.value, environment, context)
            return _element(receiver), behaviour
        if call.func.attr == "dot":
            receiver, behaviour = type_expression(call.func.value, environment, context)
            argument, argument_behaviour = type_expression(call.args[0], environment, context)
            return _apply_library_primitive("matmul", "left", [receiver, argument], call, context, sequence(behaviour, argument_behaviour))
        raise TypeError(f"unsupported method `{function_name}` at {position}")
    if function_name in environment and environment[function_name].type.kind == "fn":     # user function: inline it
        return inline_call(environment[function_name].function, call, environment, context)
    if function_name not in LIBRARY:
        raise TypeError(f"unknown primitive `{function_name}` at {position}; add it to LIBRARY with its type-and-effect signature")
    kind, strategy = LIBRARY[function_name]
    arguments, behaviour = [], EMPTY
    for argument_expression in call.args:
        argument, argument_behaviour = type_expression(argument_expression, environment, context)
        arguments.append(argument)
        behaviour = sequence(behaviour, argument_behaviour)
    for keyword in call.keywords:
        _, keyword_behaviour = type_expression(keyword.value, environment, context)
        behaviour = sequence(behaviour, keyword_behaviour)
    return _apply_library_primitive(kind, strategy, arguments, call, context, behaviour)


def _apply_library_primitive(kind: str, strategy: str, arguments: List[Value], call: ast.Call,
                             context: Context, behaviour: Behaviour) -> Tuple[Value, Behaviour]:
    position = context.source_position(call)
    reals = [argument for argument in arguments if argument.type.kind == "real"]
    operands = ",".join(label(argument) for argument in call.args
                        if not (isinstance(argument, ast.Constant) and isinstance(argument.value, str)))
    provenance_in = sequence(*(argument.provenance for argument in arguments))
    if kind == "cast":
        value = arguments[0]
        atom = Atom("cast", strategy, "", operands, position)
        return Value(Type(strategy, "real", value.type.shape), sequence(value.provenance, atom)), sequence(behaviour, atom)
    if kind == "matmul":
        precision = "any"
        for real in reals:
            precision = unify_precision(precision, real.type.precision, f"{operands} at {position}")
        atoms = [Atom("x", precision, "", operands, position), Atom("+", precision, strategy, operands, position)]
        shape = ()
        if len(reals) >= 2 and all(real.type.shape for real in reals):
            shape = reals[0].type.shape[:1] + reals[-1].type.shape[-1:]
        same_operand_twice = (len(reals) >= 2 and reals[0] is reals[-1]) or \
            (len(call.args) >= 2 and ast.unparse(call.args[-1]) == ast.unparse(call.args[-2]))
        if same_operand_twice:
            shape = reals[0].type.shape[:1]                                   # einsum('ij,ij->i', X, X): row norms
        return Value(Type(precision, "real", shape), sequence(provenance_in, *atoms)), sequence(behaviour, *atoms)
    if kind == "reduce":
        value = arguments[0]
        if value.type.kind != "real":                                         # counting booleans / ints: exact
            return Value(INT, provenance_in), behaviour
        atom = Atom("+", value.type.precision, strategy, operands, position)
        shape = value.type.shape[1:] if len(value.type.shape) > 1 and "row" in ast.unparse(call.func) else ()
        return Value(Type(value.type.precision, "real", shape), sequence(provenance_in, atom)), sequence(behaviour, atom)
    if kind == "sqrt":
        value = arguments[0]
        atom = Atom("sqrt", value.type.precision, "", operands, position)
        return Value(Type(value.type.precision, "real", value.type.shape), sequence(value.provenance, atom)), sequence(behaviour, atom)
    if kind in ("argmin", "argsort", "argmax"):
        value = arguments[0]
        precision = value.type.precision if value.type.kind == "real" else "int"
        atom = Atom(kind, precision, "", operands, position)
        context.decisions.append(Decision(kind, atom, sequence(*context.control, provenance_in, atom), position, operands))
        return Value(Type("", "int", value.type.shape[:1]), sequence(provenance_in, atom)), sequence(behaviour, atom)
    if kind == "same":
        value = reals[0] if reals else arguments[0]
        return Value(value.type, provenance_in), behaviour
    if kind == "like":
        return Value(arguments[0].type, EMPTY), behaviour
    if kind == "fresh_real":
        return Value(Type("any", "real", ("n",))), behaviour
    if kind == "shape":                                                       # len / range
        value = arguments[0]
        return Value(Type("", "int", ("n",)), value.provenance if value.type.kind == "int" else EMPTY), behaviour
    if kind in ("fresh_int", "count", "int"):
        return Value(Type("", "int", ("n",)), provenance_in), behaviour
    raise TypeError(kind)


def inline_call(function: ast.FunctionDef, call: ast.Call, environment: Environment, context: Context) -> Tuple[Value, Behaviour]:
    behaviour = EMPTY
    local: Environment = {name: value for name, value in environment.items()
                          if value.type.kind == "fn" or name.startswith("np") or name.startswith("math")}
    for parameter, argument_expression in zip(function.args.args, call.args):
        argument, argument_behaviour = type_expression(argument_expression, environment, context)
        behaviour = sequence(behaviour, argument_behaviour)
        local[parameter.arg] = argument
    context.function_stack.append(function.name)
    context.returns.append([])
    context.inlined.append([])
    body = type_block(function.body, local, context)
    context.inlined.pop()
    returned = context.returns.pop()
    context.function_stack.pop()
    context.inlined[-1].append((f"{function.name}({', '.join(ast.unparse(argument) for argument in call.args)})", body))
    result = Value(INT)
    for returned_value in returned:
        result = join_value(result, returned_value) if result is not None and result.type != INT else returned_value
    return result, sequence(behaviour, *(statement.behaviour for statement in body))


def _bind(target: ast.AST, value: Value, environment: Environment, context: Context):
    provenance = sequence(*context.control, value.provenance)
    if isinstance(target, ast.Name):
        environment[target.id] = Value(value.type, provenance, value.elements, value.function)
    elif isinstance(target, ast.Subscript):                                   # partial update: join with the old value
        base = ast.unparse(target.value)
        if base not in environment:
            raise TypeError(f"store into unknown `{base}`")
        current = environment[base]
        environment[base] = Value(join_type(current.type, Type(value.type.precision, value.type.kind, current.type.shape)),
                                  choice(current.provenance, provenance))
    elif isinstance(target, ast.Tuple):
        for element_target, element_value in zip(target.elts, value.elements):
            _bind(element_target, element_value, environment, context)
    else:
        raise TypeError(f"unsupported assignment target `{ast.unparse(target)}`")


def type_block(statements: List[ast.stmt], environment: Environment, context: Context) -> List[TypedStatement]:
    typed: List[TypedStatement] = []
    for statement in statements:
        position = context.source_position(statement)
        line = statement.lineno
        decisions_before = len(context.decisions)
        context.inlined.append([])
        if isinstance(statement, ast.FunctionDef):
            environment[statement.name] = Value(FUNCTION, function=statement)
            typed_statement = TypedStatement("def", f"def {statement.name}(...)", EMPTY, line, position)
        elif isinstance(statement, (ast.Pass, ast.Break, ast.Continue)):
            typed_statement = TypedStatement("pass", ast.unparse(statement), EMPTY, line, position)
        elif isinstance(statement, ast.Return):
            value, behaviour = type_expression(statement.value, environment, context)
            if context.returns:
                context.returns[-1].append(Value(value.type, sequence(*context.control, value.provenance), value.elements))
            typed_statement = TypedStatement("return", ast.unparse(statement), behaviour, line, position)
        elif isinstance(statement, ast.Expr):
            _, behaviour = type_expression(statement.value, environment, context)
            typed_statement = TypedStatement("expr", ast.unparse(statement), behaviour, line, position)
        elif isinstance(statement, (ast.Assign, ast.AugAssign)):
            if isinstance(statement, ast.AugAssign):
                value_expression = ast.BinOp(left=copy.deepcopy(statement.target), op=statement.op, right=statement.value)
                ast.copy_location(value_expression, statement)
                ast.fix_missing_locations(value_expression)
                target = statement.target
            else:
                value_expression, target = statement.value, statement.targets[0]
            value, behaviour = type_expression(value_expression, environment, context)
            _bind(target, value, environment, context)
            typed_statement = TypedStatement("assign", ast.unparse(statement), behaviour, line, position)
        elif isinstance(statement, ast.If):                                   # if: test ; (then | else)
            test, test_behaviour = type_expression(statement.test, environment, context)
            then_environment, else_environment = dict(environment), dict(environment)
            context.control.append(test.provenance)
            context.depth += 1
            then_branch = type_block(statement.body, then_environment, context)
            else_branch = type_block(statement.orelse, else_environment, context)
            context.depth -= 1
            context.control.pop()
            environment.clear()
            environment.update(merge_environments(then_environment, else_environment))
            then_behaviour = sequence(*(branch_statement.behaviour for branch_statement in then_branch))
            else_behaviour = sequence(*(branch_statement.behaviour for branch_statement in else_branch))
            typed_statement = TypedStatement(
                "if", "if " + ast.unparse(statement.test) + ":", sequence(test_behaviour, choice(then_behaviour, else_behaviour)),
                line, position, then_branch=then_branch, else_branch=else_branch,
                effect_note="" if then_behaviour != else_behaviour else "(if==else)")
        elif isinstance(statement, (ast.For, ast.While)):
            typed_statement = type_loop(statement, environment, context)
        else:
            raise TypeError(f"unsupported statement `{ast.unparse(statement).splitlines()[0]}` at {position}")
        typed_statement.inlined = context.inlined.pop()
        typed_statement.is_decision = any(decision.source_position == position
                                          for decision in context.decisions[decisions_before:])
        typed.append(typed_statement)
    return typed


def type_loop(loop, environment: Environment, context: Context) -> TypedStatement:
    position, line = context.source_position(loop), loop.lineno
    loop_name = "loop" if context.depth == 0 else f"loop{context.depth}"
    carried = sorted(name for name in _assigned_names(loop.body) if name in environment)
    before = {name: environment[name] for name in carried}
    decisions_before = len(context.decisions)
    for name in carried:
        environment[name] = Value(environment[name].type, Variable(f"loop_{name}"), environment[name].elements)
    is_for = isinstance(loop, ast.For)
    if is_for:
        iterable, _ = type_expression(loop.iter, environment, context)
        if not isinstance(loop.target, ast.Name):
            raise TypeError(f"unsupported loop target at {position}")
        loop_variable = loop.target.id
    body: List[TypedStatement] = []
    body_environment = dict(environment)
    condition_behaviour = EMPTY
    for _ in range(4):
        del context.decisions[decisions_before:]
        trial = dict(body_environment)
        if is_for:
            trial[loop_variable] = _element(iterable)
            control = []
        else:
            condition, condition_behaviour = type_expression(loop.test, trial, context)
            control = [condition.provenance]
        context.control.extend(control)
        context.depth += 1
        body = type_block(loop.body, trial, context)
        context.depth -= 1
        del context.control[len(context.control) - len(control):]
        merged = dict(body_environment)
        for name in carried:
            merged[name] = Value(join_type(body_environment[name].type, trial[name].type), Variable(f"loop_{name}"), trial[name].elements)
        stable = all(merged[name].type == body_environment[name].type for name in carried)
        body_environment = merged
        if stable:
            break

    def close(provenance: Behaviour, skip: str = "") -> Behaviour:
        for name in carried:
            if name != skip:
                provenance = substitute(provenance, f"loop_{name}", after[name])
        return normalise(provenance)

    after: Dict[str, Behaviour] = {name: recursion(f"loop_{name}", choice(before[name].provenance, trial[name].provenance))
                                   for name in carried}
    for name in carried:
        after[name] = close(after[name], skip=name)
    for name in carried:
        environment[name] = Value(join_type(before[name].type, trial[name].type), after[name], trial[name].elements)
    for name, value in trial.items():                   # names first assigned inside the loop
        if name not in environment and (not is_for or name != loop_variable):
            environment[name] = Value(value.type, close(value.provenance), value.elements, value.function)
    for decision in context.decisions[decisions_before:]:
        decision.slice = close(decision.slice)
    body_behaviour = sequence(*(statement.behaviour for statement in body))
    if is_for:
        loop_behaviour = recursion(loop_name, choice(EMPTY, sequence(body_behaviour, Variable(loop_name))))
        source = f"for {ast.unparse(loop.target)} in {ast.unparse(loop.iter)}:"
    else:
        loop_behaviour = recursion(loop_name, sequence(condition_behaviour, choice(EMPTY, sequence(body_behaviour, Variable(loop_name)))))
        source = f"while {ast.unparse(loop.test)}:"
    accumulators = [name for name in carried if isinstance(after[name], Sequence)
                    and any(atom.strategy == "left" for atom in atoms_in_order(after[name]))
                    and not any(atom.strategy == "left" for atom in atoms_in_order(before[name].provenance))]
    note = ""
    if accumulators:
        note = f"(accumulator{'s' if len(accumulators) > 1 else ''} {', '.join(accumulators)}: a sequential sum, written +^left)"
    return TypedStatement("for" if is_for else "while", source, loop_behaviour, line, position, effect_note=note, body=body)


@dataclass
class Trace:
    label: str
    source: str
    block: List[TypedStatement]
    decisions: List[Decision]
    behaviour: Behaviour
    inputs: Dict[str, Type]


def analyse(source: str, inputs: Dict[str, Type], trace_label: str) -> Trace:
    source = textwrap.dedent(source).strip("\n")
    tree = ast.parse(source)
    context = Context()
    environment: Environment = {name: Value(input_type, Atom("in", "", "", name)) for name, input_type in inputs.items()}
    block = type_block(tree.body, environment, context)
    for statement in reversed(tree.body):
        if isinstance(statement, ast.Assign):
            position = context.source_position(statement)
            for name in sorted(_target_names(statement.targets[0])):
                atom = Atom("out", "", "", name, position)
                context.decisions.append(Decision("out", atom, sequence(environment[name].provenance, atom), position, name))
            break
    return Trace(trace_label, source, block, context.decisions, sequence(*(statement.behaviour for statement in block)), inputs)


# ------ Print ------

def _statements_in_order(block: List[TypedStatement], depth: int = 0, function: str = ""):
    for statement in block:
        yield statement, depth, function
        for callee_name, callee_body in statement.inlined:
            yield from _statements_in_order(callee_body, depth + 1, callee_name)
        yield from _statements_in_order(statement.then_branch, depth + 1, function)
        yield from _statements_in_order(statement.else_branch, depth + 1, function)
        yield from _statements_in_order(statement.body, depth + 1, function)


def print_analysis(trace: Trace):
    print(f"-- {trace.label}")
    print(f"inputs: " + ", ".join(f"{name} : {input_type}" for name, input_type in trace.inputs.items()))
    print()
    for statement, depth, function in _statements_in_order(trace.block):
        indent = "  " * depth
        source = statement.source
        function_tag = f"  [in {function}]" if function else ""
        decision_tag = ""
        function_tag = ""
        print(f"L{statement.line:<3} {indent}{source}{function_tag}{decision_tag}")
        if isinstance(statement.behaviour, Recursion):
            effect = f"rec ... {compact(statement.behaviour)}"
        else:
            effect = show(statement.behaviour)
        if effect != "empty" or statement.effect_note:
            print(f"effect: {effect} {statement.effect_note}")
            
    for index, decision in enumerate(trace.decisions):
        slice_text = " ; ".join(shown(atoms_in_order(decision.slice))) or "empty"
        print(f"- D{index + 1:<2} {show_atom(decision.atom):<30} {decision.source_position:<24}"
              f"steps: {slice_text} ")
    print()
    print("-" * 100)
    print()
    print()


def analyse_pair(title: str, source_a: str, source_b: str, inputs: Dict[str, Type], name_a: str, name_b: str):
    print("=" * 104)
    print(f"  {title}")
    print("=" * 104)
    trace_a = analyse(source_a, inputs, name_a)
    trace_b = analyse(source_b, inputs, name_b)
    print_analysis(trace_a)
    print_analysis(trace_b)
    print()
    return trace_a, trace_b


DBSCAN_SKLEARN = """
def rdist(a, b):
    d = 0.0
    for k in range(len(a)):
        tmp = a[k] - b[k]
        d += tmp * tmp
    return d

def query_radius(tree_data, i, r2):
    neighbors = []
    for j in range(len(tree_data)):
        if rdist(tree_data[i], tree_data[j]) <= r2:
            neighbors.append(j)
    return neighbors

def dbscan_fit(X, eps, min_samples):
    tree_data = X.astype(np.float64)
    r2 = eps * eps
    n = len(X)
    neighborhoods = []
    for i in range(n):
        neighborhoods.append(query_radius(tree_data, i, r2))
    core = np.zeros(n)
    for i in range(n):
        core[i] = len(neighborhoods[i]) >= min_samples
    labels = np.full(n, -1)
    cluster = 0
    for i in range(n):
        if labels[i] != -1 or not core[i]:
            continue
        stack = [i]
        while stack:
            p = stack.pop()
            if labels[p] == -1:
                labels[p] = cluster
                if core[p]:
                    for q in neighborhoods[p]:
                        if labels[q] == -1:
                            stack.append(q)
        cluster += 1
    return labels

labels = dbscan_fit(X, eps, min_samples)
"""

DBSCAN_CUML = """
def dbscan_fit(X, eps, min_samples):
    n = X.shape[0]
    eps2 = np.float32(eps * eps)
    # vertexdeg/algo.cuh: raft::distance L2 expanded, fused with the GEMM
    norms = gpu_rowsum(X * X)       
    dots = gpu_gemm(X, X.T)
    D2 = norms[:, None] + norms[None, :] - 2.0 * dots
    adj = D2 <= eps2
    deg = gpu_rowsum(adj)
    core = deg >= min_samples
    labels = np.full(n, -1)
    cluster = 0
    for i in range(n):
        if labels[i] != -1 or not core[i]:
            continue
        frontier = [i]
        while frontier:
            p = frontier.pop()
            if labels[p] == -1:
                labels[p] = cluster
                if core[p]:
                    for q in range(n):
                        if adj[p][q] and labels[q] == -1:
                            frontier.append(q)
        cluster += 1
    return labels

labels = dbscan_fit(X, eps, min_samples)
"""

if __name__ == "__main__":
    analyse_pair(
        "DBSCAN", DBSCAN_SKLEARN, DBSCAN_CUML,
        inputs={"X": real32("n", "d"), "eps": real64(), "min_samples": INT},
        name_a="sklearn", name_b="cuML",
    )
