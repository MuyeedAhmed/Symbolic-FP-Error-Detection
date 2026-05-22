import requests
import os
import zipfile
import io
import ast
import time
import pandas as pd

GITHUB_API_URL = "https://api.github.com/search/repositories"
TOKEN = "TOKEN_PLACEHOLDER"
GITHUB_TOKEN = os.environ.get("GITHUB_TOKEN", TOKEN)

HEADERS = {"Authorization": f"token {GITHUB_TOKEN}"} if GITHUB_TOKEN != "YOUR_GITHUB_TOKEN_HERE" else {}

EXTRACTED_BASE_DIR = "ExtractedFunctions"

def get_top_python_repos(limit=100):
    repos = []
    page = 1
    while len(repos) < limit:
        params = {
            "q": "language:python",
            "sort": "stars",
            "order": "desc",
            "per_page": min(limit - len(repos), 100),
            "page": page
        }
        response = requests.get(GITHUB_API_URL, params=params, headers=HEADERS)
        if response.status_code != 200:
            print(f"Error: {response.status_code} {response.text}")
            break
        data = response.json()
        repos.extend(data.get("items", []))
        if len(data.get("items", [])) == 0:
            break
        page += 1
        if not HEADERS: time.sleep(2) 
    return repos[:limit]

def search_linalg_in_file(content):
    class LinearAlgebraVisitor(ast.NodeVisitor):
        def __init__(self, source_lines):
            self.source_lines = source_lines
            self.current_function_node = None
            self.found_functions = {}
            self.counts = {'dot': 0, 'inv': 0, 'matmul': 0}
            self.details = []
            self.numpy_aliases = set()
            self.linalg_aliases = set()
            self.dot_imported = False
            self.inv_imported = False

        def visit_Import(self, node):
            for alias in node.names:
                if alias.name == 'numpy': self.numpy_aliases.add(alias.asname or 'numpy')
                if alias.name == 'numpy.linalg': self.linalg_aliases.add(alias.asname or 'linalg')
            self.generic_visit(node)

        def visit_ImportFrom(self, node):
            if node.module == 'numpy':
                for alias in node.names:
                    if alias.name == 'dot': self.dot_imported = True
                    if alias.name == 'linalg': self.linalg_aliases.add(alias.asname or 'linalg')
            if node.module == 'numpy.linalg':
                for alias in node.names:
                    if alias.name == 'inv': self.inv_imported = True
            self.generic_visit(node)

        def visit_FunctionDef(self, node):
            old = self.current_function_node
            self.current_function_node = node
            self.generic_visit(node)
            self.current_function_node = old

        def visit_AsyncFunctionDef(self, node):
            old = self.current_function_node
            self.current_function_node = node
            self.generic_visit(node)
            self.current_function_node = old

        def visit_BinOp(self, node):
            if isinstance(node.op, ast.MatMult):
                self.counts['matmul'] += 1
                if self.current_function_node:
                    self._extract(self.current_function_node, 'matmul')
            self.generic_visit(node)

        def visit_Call(self, node):
            type_found = None
            if isinstance(node.func, ast.Attribute):
                if node.func.attr == 'dot': type_found = 'dot'
                elif node.func.attr == 'inv': type_found = 'inv'
            elif isinstance(node.func, ast.Name):
                if node.func.id == 'dot' and self.dot_imported: type_found = 'dot'
                elif node.func.id == 'inv' and self.inv_imported: type_found = 'inv'

            if type_found:
                self.counts[type_found] += 1
                if self.current_function_node:
                    self._extract(self.current_function_node, type_found)
            self.generic_visit(node)

        def _extract(self, node, op_type):
            start, end = node.lineno - 1, getattr(node, 'end_lineno', node.lineno)
            self.found_functions[node.name] = "\n".join(self.source_lines[start:end])
            self.details.append({'function': node.name, 'type': op_type})

    try:
        source_lines = content.splitlines()
        tree = ast.parse(content)
        v = LinearAlgebraVisitor(source_lines)
        v.visit(tree)
        return v.counts, v.details, v.found_functions
    except: return {'dot': 0, 'inv': 0, 'matmul': 0}, [], {}

def main():
    if not os.path.exists(EXTRACTED_BASE_DIR): os.makedirs(EXTRACTED_BASE_DIR)
    repos = get_top_python_repos(100)
    summary, detailed = [], []
    for repo in repos:
        name = repo['full_name']
        repo_summary = {'repo': name, 'dot': 0, 'inv': 0, 'matmul': 0}
        url = f"https://github.com/{name}/archive/refs/heads/{repo['default_branch']}.zip"
        try:
            res = requests.get(url, timeout=60)
            if res.status_code == 200:
                with zipfile.ZipFile(io.BytesIO(res.content)) as z:
                    for info in z.infolist():
                        if info.filename.endswith(".py"):
                            with z.open(info) as f:
                                try:
                                    content = f.read().decode('utf-8', errors='ignore')
                                    counts, details, funcs = search_linalg_in_file(content)
                                    for k in repo_summary:
                                        if k != 'repo': repo_summary[k] += counts[k]
                                    for d in details:
                                        detailed.append({'repo': name, 'file': info.filename, 'function': d['function'], 'type': d['type']})
                                    if funcs:
                                        repo_dir = os.path.join(EXTRACTED_BASE_DIR, name.replace("/", "_"))
                                        if not os.path.exists(repo_dir): os.makedirs(repo_dir)
                                        for fn, src in funcs.items():
                                            with open(os.path.join(repo_dir, f"{info.filename.replace('/', '_')}_{fn}.py"), "w", encoding="utf-8") as out:
                                                out.write(src)
                                except: continue
        except: pass
        summary.append(repo_summary)
        print(f"{name}: dot={repo_summary['dot']}, inv={repo_summary['inv']}, matmul={repo_summary['matmul']}")
        time.sleep(1)

    os.makedirs("Results", exist_ok=True)
    pd.DataFrame(summary).to_excel("Results/Summary_Top_100.xlsx", index=False)
    pd.DataFrame(detailed).to_excel("Results/Detailed_Top_100.xlsx", index=False)

if __name__ == "__main__":
    main()
