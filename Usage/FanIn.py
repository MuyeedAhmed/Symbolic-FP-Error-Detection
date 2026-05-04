import ast
import os
import argparse
import pandas as pd
from collections import defaultdict
from pathlib import Path

class CallFinder(ast.NodeVisitor):
    def __init__(self, target_names):
        self.target_names = target_names
        self.found_calls = []
        self.current_function = None
        self.current_class = None
        self.imports = {}

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
        old_func = self.current_function
        self.current_function = node.name
        self.generic_visit(node)
        self.current_function = old_func

    def visit_AsyncFunctionDef(self, node):
        self.visit_FunctionDef(node)

    def visit_Call(self, node):
        func_name = self._get_func_name(node.func)
        if func_name:
            for target in self.target_names:
                if func_name == target or func_name.endswith('.' + target):
                    self.found_calls.append({
                        'function': self.current_function,
                        'class': self.current_class,
                        'line': node.lineno,
                        'target': target,
                        'full_call': func_name
                    })
                    break
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
            if op_sym in self.target_names:
                self.found_calls.append({
                    'function': self.current_function,
                    'class': self.current_class,
                    'line': node.lineno,
                    'target': op_sym,
                    'full_call': op_sym
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
    all_results = []
    abs_root = os.path.abspath(directory)
    
    for root, _, files in os.walk(directory):
        for file in files:
            if file.endswith('.py'):
                path = os.path.join(root, file)
                try:
                    with open(path, 'r', encoding='utf-8') as f:
                        lines = f.readlines()
                        tree = ast.parse("".join(lines), filename=path)
                    
                    finder = CallFinder(targets)
                    finder.visit(tree)
                    
                    module_path = get_module_path(os.path.abspath(path), abs_root)
                    
                    for call in finder.found_calls:
                        line_idx = call['line'] - 1
                        source_line = lines[line_idx].strip() if 0 <= line_idx < len(lines) else "N/A"
                        
                        caller_func = call['function'] or "Global scope"
                        caller_class = call['class'] or "N/A"
                        
                        module_parts = module_path.split('.')
                        module_public = all(not p.startswith('_') or p == '__init__' for p in module_parts)
                        func_public = is_public(caller_func)
                        class_public = is_public(caller_class)
                        
                        visibility = "Public" if (module_public and func_public and class_public) else "Internal"
                        
                        all_results.append({
                            'Target': call['target'],
                            'Visibility': visibility,
                            'Module Path': module_path,
                            'Function/Method': f"{caller_class}.{caller_func}" if caller_class != "N/A" else caller_func,
                            'Source Code': source_line,
                            'File': path,
                            'Line': call['line'],
                            'Full Call': call['full_call']
                        })
                except Exception as e:
                    print(f"Error parsing {path}: {e}")
    
    return all_results

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("directory", help="Directory to search (e.g., path to numpy/scipy)")
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
    df = df[cols]
    
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
