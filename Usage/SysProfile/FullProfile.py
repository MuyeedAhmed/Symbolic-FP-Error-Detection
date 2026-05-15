import sys
import os
import argparse
import pandas as pd
import inspect
from collections import defaultdict
import runpy
import traceback
from pathlib import Path


# python Usage/FullProfile.py -m pytest --pyargs scipy.linalg --ignore=/Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/linalg/tests/_cython_examples

class FullAnalyzer:
    def __init__(self):
        self.calls = defaultdict(int)
        self.executions = defaultdict(int)
        self.metadata = {}

    def _get_key(self, code):
        return (code.co_filename, code.co_name, code.co_firstlineno)

    def profile_func(self, frame, event, arg):
        if event == 'call':
            code = frame.f_code
            key = self._get_key(code)
            self.executions[key] += 1
            
            if key not in self.metadata:
                self.metadata[key] = self._get_info(code, frame)
            
            parent = frame.f_back
            if parent:
                p_key = self._get_key(parent.f_code)
                if p_key not in self.metadata:
                    self.metadata[p_key] = self._get_info(parent.f_code, parent)
                self.calls[(p_key, key)] += 1

        elif event == 'c_call':
            p_key = self._get_key(frame.f_code)
            if p_key not in self.metadata:
                self.metadata[p_key] = self._get_info(frame.f_code, frame)
            
            c_name = arg.__name__ if hasattr(arg, '__name__') else str(arg)
            if ' ' in c_name:
                c_name = c_name.split(' ')[-1].rstrip('>')
            
            c_key = ("<C-function>", c_name, 0)
            if c_key not in self.metadata:
                self.metadata[c_key] = {
                    'file': '<C-function>',
                    'name': c_name,
                    'class': None,
                    'module': 'builtin',
                    'size': 0
                }
            
            self.executions[c_key] += 1
            self.calls[(p_key, c_key)] += 1
        
        return self.profile_func

    def _get_info(self, code, frame):
        cls_name = None
        try:
            if 'self' in frame.f_locals:
                cls_name = frame.f_locals['self'].__class__.__name__
            elif 'cls' in frame.f_locals:
                cls_name = frame.f_locals['cls'].__name__
        except:
            pass

        size = 0
        try:
            lines, _ = inspect.getsourcelines(code)
            size = len(lines)
        except:
            pass

        return {
            'file': code.co_filename,
            'name': code.co_name,
            'class': cls_name,
            'module': frame.f_globals.get('__name__', 'unknown'),
            'size': size
        }

def main():
    parser = argparse.ArgumentParser(description="Full Dynamic Profiler")
    parser.add_argument("-m", "--module", action="store_true", help="Run library module as a script")
    parser.add_argument("--output", default="Full_Profile.xlsx", help="Output Excel filename")
    parser.add_argument("--exclude-libs", action="store_true", help="Exclude calls where both caller and callee are in libs")
    parser.add_argument("script", help="Script path to run or module name if -m is used")
    parser.add_argument("script_args", nargs=argparse.REMAINDER, help="Arguments for the script")
    
    args = parser.parse_args()
    if not args.script:
        parser.print_help()
        return

    analyzer = FullAnalyzer()
    
    orig_argv = sys.argv
    orig_path = sys.path[:]
    
    if args.module:
        target_mod = args.script
        sys.argv = [args.script] + args.script_args
    else:
        sys.argv = [args.script] + args.script_args
        target_mod = None
        script_dir = os.path.dirname(os.path.abspath(args.script))
        sys.path.insert(0, script_dir)

    sys.setprofile(analyzer.profile_func)
    try:
        if target_mod:
            runpy.run_module(target_mod, run_name="__main__", alter_sys=True)
        else:
            runpy.run_path(args.script, run_name="__main__")
    except SystemExit:
        pass
    except Exception:
        print("\n--- Target Script Error Exp ---")
        traceback.print_exc()
        print("-------------------------------\n")
    finally:
        sys.setprofile(None)
        sys.argv = orig_argv
        sys.path = orig_path

    results = []
    def is_lib(path):
        return 'site-packages' in path or 'lib/python' in path or path == '<C-function>'

    for (p_key, c_key), call_count in analyzer.calls.items():
        p_meta = analyzer.metadata.get(p_key, {})
        c_meta = analyzer.metadata.get(c_key, {})
        
        p_file = p_meta.get('file', 'unknown')
        c_file = c_meta.get('file', 'unknown')
        
        if args.exclude_libs:
            if is_lib(p_file) and is_lib(c_file):
                continue

        p_mod = p_meta.get('module', 'unknown')
        p_cls = p_meta.get('class')
        p_func = p_meta.get('name', 'unknown')
        p_path = f"{p_mod}.{p_cls}.{p_func}" if p_cls else f"{p_mod}.{p_func}"

        c_mod = c_meta.get('module', 'unknown')
        c_cls = c_meta.get('class')
        c_func = c_meta.get('name', 'unknown')
        c_path = f"{c_mod}.{c_cls}.{c_func}" if c_cls else f"{c_mod}.{c_func}"
        
        results.append({
            'Caller Path': p_path,
            'Caller File': p_file,
            'Callee Path': c_path,
            'Callee File': c_file,
            'Call Count': call_count,
            'Callee Total Executions': analyzer.executions.get(c_key, 0),
            'Callee Size': c_meta.get('size', 0)
        })

    if results:
        df = pd.DataFrame(results)
        df = df.sort_values(by='Call Count', ascending=False)
        try:
            df.to_excel(args.output, index=False)
            print(f"Results saved to {args.output}")
        except Exception as e:
            csv_output = args.output.replace('.xlsx', '.csv')
            df.to_csv(csv_output, index=False)
            print(f"Failed to save to Excel (possibly too many rows), saved to {csv_output} instead. Error: {e}")
    else:
        print("No profiling data collected.")

if __name__ == "__main__":
    main()
