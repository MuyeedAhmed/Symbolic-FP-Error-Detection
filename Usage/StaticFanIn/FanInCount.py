import ast
import os
import argparse
import pandas as pd
from collections import defaultdict
from pathlib import Path

class CallFinder(ast.NodeVisitor):
    def __init__(self):
        self.definitions = []
        self.calls = []
        self.imports = {}
        self.current_function = None
        self.current_class = None

    def visit_Import(self, node):
        for alias in node.names:
            self.imports[alias.asname or alias.name] = alias.name
        self.generic_visit(node)

    def visit_ImportFrom(self, node):
        module = node.module
        for alias in node.names:
            full_name = f"{module}.{alias.name}" if module else alias.name
            self.imports[alias.asname or alias.name] = full_name
        self.generic_visit(node)

    def visit_ClassDef(self, node):
        old_class = self.current_class
        self.current_class = node.name
        self.generic_visit(node)
        self.current_class = old_class

    def visit_FunctionDef(self, node):
        self.definitions.append({
            'name': node.name,
            'class': self.current_class,
            'line': node.lineno,
            'end_line': getattr(node, 'end_lineno', node.lineno)
        })
        old_func = self.current_function
        self.current_function = node.name
        self.generic_visit(node)
        self.current_function = old_func

    def visit_AsyncFunctionDef(self, node):
        self.visit_FunctionDef(node)

    def visit_Call(self, node):
        func_name = self._get_func_name(node.func)
        if func_name:
            self.calls.append({
                'called_name': func_name,
                'caller_func': self.current_function,
                'caller_class': self.current_class,
                'line': node.lineno
            })
        self.generic_visit(node)

    def visit_BinOp(self, node):
        op_map = {
            ast.MatMult: '@',
            ast.Add: '+',
            ast.Sub: '-',
            ast.Mult: '*',
            ast.Div: '/',
        }
        op_type = type(node.op)
        if op_type in op_map:
            op_sym = op_map[op_type]
            self.calls.append({
                'called_name': op_sym,
                'caller_func': self.current_function,
                'caller_class': self.current_class,
                'line': node.lineno
            })
        self.generic_visit(node)

    def _get_func_name(self, node):
        if isinstance(node, ast.Name):
            return node.id
        elif isinstance(node, ast.Attribute):
            value = self._get_func_name(node.value)
            if value:
                return f"{value}.{node.attr}"
        return None

def get_module_path(file_path, root_dir):
    try:
        rel_path = os.path.relpath(file_path, root_dir)
        parts = list(Path(rel_path).parts)
        if parts[-1] == '__init__.py':
            parts.pop()
        else:
            parts[-1] = os.path.splitext(parts[-1])[0]
        
        pkg_name = os.path.basename(root_dir.rstrip(os.sep))
        return f"{pkg_name}." + ".".join(parts)
    except:
        return "Unknown"

def analyze_directory(directory, targets):
    abs_root = os.path.abspath(directory)
    codebase_definitions = {}
    codebase_calls = []

    for root, _, files in os.walk(directory):
        for file in files:
            if file.endswith('.py'):
                path = os.path.join(root, file)
                try:
                    with open(path, 'r', encoding='utf-8') as f:
                        tree = ast.parse(f.read(), filename=path)
                    
                    module_path = get_module_path(os.path.abspath(path), abs_root)
                    finder = CallFinder()
                    finder.visit(tree)
                    
                    for defn in finder.definitions:
                        key = (module_path, defn['class'], defn['name'])
                        codebase_definitions[key] = {
                            'file': path, 
                            'line': defn['line'],
                            'size': defn['end_line'] - defn['line'] + 1
                        }
                    
                    for call in finder.calls:
                        codebase_calls.append({
                            'module': module_path,
                            'class': call['caller_class'],
                            'function': call['caller_func'],
                            'called_name': call['called_name'],
                            'line': call['line'],
                            'imports': finder.imports
                        })
                except Exception as e:
                    print(f"Error parsing {path}: {e}")

    call_graph = defaultdict(list)
    def_by_name = defaultdict(list)
    for key in codebase_definitions:
        def_by_name[key[2]].append(key)

    fan_in_counter = defaultdict(int)

    for call in codebase_calls:
        caller_key = (call['module'], call['class'], call['function'] or "Global scope")
        called_name = call['called_name']
        resolved = None
        
        for target in targets:
            if called_name == target or called_name.endswith('.' + target):
                resolved = {'type': 'target', 'value': target}
                break
        
        if not resolved:
            key_in_module = (call['module'], call['class'], called_name)
            if key_in_module in codebase_definitions:
                resolved = {'type': 'func', 'value': key_in_module}
            else:
                key_top_level = (call['module'], None, called_name)
                if key_top_level in codebase_definitions:
                    resolved = {'type': 'func', 'value': key_top_level}

            if not resolved and called_name.startswith('self.'):
                method_name = called_name[len('self.'):]
                key_in_same_class = (call['module'], call['class'], method_name)
                if key_in_same_class in codebase_definitions:
                    resolved = {'type': 'func', 'value': key_in_same_class}

            if not resolved:
                parts = called_name.split('.')
                prefix = parts[0]
                if prefix in call['imports']:
                    resolved_prefix = call['imports'][prefix]
                    full_resolved_name = ".".join([resolved_prefix] + parts[1:])
                    for key in codebase_definitions:
                        fqn = f"{key[0]}.{key[1]}.{key[2]}" if key[1] else f"{key[0]}.{key[2]}"
                        if fqn == full_resolved_name:
                            resolved = {'type': 'func', 'value': key}
                            break
            
            if not resolved:
                if called_name in def_by_name:
                    resolved = {'type': 'func', 'value': def_by_name[called_name][0]}

            if not resolved:
                parts = called_name.split('.')
                if len(parts) > 1:
                    suffix = parts[-1]
                    if suffix in def_by_name:
                        resolved = {'type': 'func', 'value': def_by_name[suffix][0]}

        if resolved:
            call_graph[caller_key].append(resolved)
            if resolved['type'] == 'func':
                fan_in_counter[resolved['value']] += 1

    reachable = defaultdict(set)
    for caller_key, resolved_calls in call_graph.items():
        for res in resolved_calls:
            if res['type'] == 'target':
                reachable[caller_key].add(res['value'])

    changed = True
    while changed:
        changed = False
        for caller_key, resolved_calls in call_graph.items():
            for res in resolved_calls:
                if res['type'] == 'func':
                    callee_key = res['value']
                    for target in reachable[callee_key]:
                        if target not in reachable[caller_key]:
                            reachable[caller_key].add(target)
                            changed = True

    final_data = []
    for caller_key, targets_reached in reachable.items():
        mod, cls, func = caller_key
        method_path = f"{mod}.{cls}.{func}" if cls else f"{mod}.{func}"
        
        size = 0
        file_path = "Unknown"
        if caller_key in codebase_definitions:
            size = codebase_definitions[caller_key]['size']
            file_path = codebase_definitions[caller_key]['file']
        
        fan_in = fan_in_counter[caller_key]
        
        for t in targets_reached:
            final_data.append({
                'Method Path': method_path,
                'Target': t,
                'Fan-In Count': fan_in,
                'Function Size': size,
                'File Path': file_path
            })
            
    return final_data

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("directory")
    parser.add_argument("--targets", nargs="+", default=['matmul', 'inv', 'det', '@'])
    parser.add_argument("--output", default="Detailed_FanIn.xlsx")
    args = parser.parse_args()

    results = analyze_directory(args.directory, args.targets)
    if results:
        df = pd.DataFrame(results)
        df.to_excel(args.output, index=False)
        print(f"Results saved to {args.output}")
    else:
        print("No paths found.")

if __name__ == "__main__":
    main()
