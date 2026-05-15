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
            'line': node.lineno
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

def is_public(name):
    if not name or name == "Global scope" or name == "N/A":
        return True
    return not name.startswith('_')

def analyze_directory(directory, targets):
    abs_root = os.path.abspath(directory)
    codebase_definitions = {}
    codebase_calls = []
    file_contents = {}

    print("Step 1")
    for root, _, files in os.walk(directory):
        for file in files:
            if file.endswith('.py'):
                path = os.path.join(root, file)
                try:
                    with open(path, 'r', encoding='utf-8') as f:
                        lines = f.readlines()
                    file_contents[path] = lines
                    tree = ast.parse("".join(lines), filename=path)
                    
                    module_path = get_module_path(os.path.abspath(path), abs_root)
                    finder = CallFinder()
                    finder.visit(tree)
                    
                    for defn in finder.definitions:
                        key = (module_path, defn['class'], defn['name'])
                        codebase_definitions[key] = {'file': path, 'line': defn['line']}
                    
                    for call in finder.calls:
                        codebase_calls.append({
                            'module': module_path,
                            'class': call['caller_class'],
                            'function': call['caller_func'],
                            'called_name': call['called_name'],
                            'line': call['line'],
                            'file': path,
                            'imports': finder.imports
                        })
                except Exception as e:
                    print(f"Error parsing {path}: {e}")

    print(f"Found {len(codebase_definitions)} definitions and {len(codebase_calls)} calls.")

    print("Step 2")
    call_graph = defaultdict(list)
    def_by_name = defaultdict(list)
    for key in codebase_definitions:
        def_by_name[key[2]].append(key)

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
            call_graph[caller_key].append({
                'resolved': resolved,
                'line': call['line'],
                'file': call['file'],
                'full_call': called_name
            })

    print("Step 3")
    reachable = defaultdict(set)
    reachable_info = defaultdict(lambda: defaultdict(list))

    for caller_key, calls in call_graph.items():
        for call in calls:
            if call['resolved']['type'] == 'target':
                target = call['resolved']['value']
                reachable[caller_key].add(target)
                if call not in reachable_info[caller_key][target]:
                    reachable_info[caller_key][target].append(call)

    changed = True
    iteration = 0
    while changed and iteration < 20:
        changed = False
        iteration += 1
        for caller_key, calls in call_graph.items():
            for call in calls:
                if call['resolved']['type'] == 'func':
                    callee_key = call['resolved']['value']
                    for target in reachable[callee_key]:
                        if target not in reachable[caller_key]:
                            reachable[caller_key].add(target)
                            reachable_info[caller_key][target].append(call)
                            changed = True

    print(f"Transitive closure completed in {iteration} iterations.")

    print("Step 4")
    all_results = []
    for caller_key, targets_found in reachable.items():
        module_path, class_name, func_name = caller_key
        
        module_parts = module_path.split('.')
        module_public = all(not p.startswith('_') or p == '__init__' for p in module_parts)
        func_public = is_public(func_name)
        class_public = is_public(class_name) or class_name is None
        visibility = "Public" if (module_public and func_public and class_public) else "Internal"

        for target in targets_found:
            for info in reachable_info[caller_key][target]:
                line_idx = info['line'] - 1
                lines = file_contents[info['file']]
                source_line = lines[line_idx].strip() if 0 <= line_idx < len(lines) else "N/A"
                
                all_results.append({
                    'Target': target,
                    'Visibility': visibility,
                    'Module Path': module_path,
                    'Function/Method': f"{class_name}.{func_name}" if class_name else func_name,
                    'Source Code': source_line,
                    'File': info['file'],
                    'Line': info['line'],
                    'Full Call': info['full_call']
                })
    
    return all_results

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("directory", help="Directory to search")
    parser.add_argument("--targets", nargs="+", default=['matmul', 'dot', 'inv', 'det', 'trace', '@'],
                        help="Default values: ['matmul', 'dot', 'inv', 'det', 'trace', '@']")
    parser.add_argument("--output", default="FanIn_Results.xlsx",
                        help="Output Excel file name.")
    
    args = parser.parse_args()
    args.directory = args.directory.rstrip(os.sep)
    
    print(f"Targets: {args.targets}")
    print(f"Directory: {args.directory}")

    results = analyze_directory(args.directory, args.targets)
    
    if not results:
        print("No calls found.")
        return

    df = pd.DataFrame(results)
    cols = ['Target', 'Visibility', 'Module Path', 'Function/Method', 'Source Code', 'File', 'Line', 'Full Call']
    df = df[cols].drop_duplicates()
    
    try:
        df.to_excel(args.output, index=False)
        print(f"\nResults successfully saved to {args.output}")
    except Exception as e:
        print(f"\nError saving to Excel: {e}")

    print("\n" + "="*20 +"Summary"+ "="*20)
    summary = df['Target'].value_counts().sort_index()
    for target, count in summary.items():
        print(f"{target:.<20} {count}")
    
    print("\nVisibility Summary:")
    vis_summary = df['Visibility'].value_counts()
    for vis, count in vis_summary.items():
        print(f"{vis:.<20} {count}")

if __name__ == "__main__":
    main()
