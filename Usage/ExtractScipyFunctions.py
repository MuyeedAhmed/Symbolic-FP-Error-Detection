import pandas as pd
import ast
import os

def get_function_source(file_path, method_path):
    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            source = f.read()
        tree = ast.parse(source)
    except Exception as e:
        return f"# Error reading or parsing {file_path}: {e}"

    parts = method_path.split('.')
    func_name = parts[-1]
    
    results = []

    class Visitor(ast.NodeVisitor):
        def __init__(self):
            self.current_class = None
            
        def visit_ClassDef(self, node):
            old_class = self.current_class
            self.current_class = node.name
            self.generic_visit(node)
            self.current_class = old_class
            
        def visit_FunctionDef(self, node):
            if node.name == func_name:
                results.append((self.current_class, node))
            self.generic_visit(node)
            
        def visit_AsyncFunctionDef(self, node):
            self.visit_FunctionDef(node)

    visitor = Visitor()
    visitor.visit(tree)
    
    if not results:
        return f"# Function {func_name} not found in {file_path}"
    
    if len(results) == 1:
        return ast.get_source_segment(source, results[0][1])
    
    if len(parts) >= 2:
        potential_class = parts[-2]
        for cls_name, node in results:
            if cls_name == potential_class:
                return ast.get_source_segment(source, node)
                
    return ast.get_source_segment(source, results[0][1])

def main():
    excel_path = 'Usage/DynamicFanInResults/Scipy_Dynamic_FanIn_Merged_ExcludingEmpty.xlsx'
    output_dir = 'Usage/ExtractedFunctions'
    output_file = os.path.join(output_dir, 'scipy_functions.py')
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        
    df = pd.read_excel(excel_path)
    
    with open(output_file, 'w', encoding='utf-8') as out:
        out.write(f"# Extracted Scipy Functions\n")
        out.write(f"# Source: {excel_path}\n\n")
        
        for index, row in df.iterrows():
            method_path = str(row['Method Path'])
            file_path = str(row['File Path'])
            target = str(row['Target'])
            call_count = str(row['Call Count'])
            total_exec = str(row['Total Executions'])
            func_size = str(row['Function Size'])
            
            divider = f"\n{'#'*80}\n"
            divider += f"# Method Path: {method_path}\n"
            divider += f"# File Path:   {file_path}\n"
            divider += f"# Target:      {target}\n"
            divider += f"# Call Count:  {call_count}\n"
            divider += f"# Total Exec:  {total_exec}\n"
            divider += f"# Func Size:   {func_size}\n"
            divider += f"{'#'*80}\n\n"
            
            print(f"[{index+1}/{len(df)}] Extracting {method_path}...")
            
            func_source = get_function_source(file_path, method_path)
            
            out.write(divider)
            out.write(func_source)
            out.write("\n")


if __name__ == "__main__":
    main()
