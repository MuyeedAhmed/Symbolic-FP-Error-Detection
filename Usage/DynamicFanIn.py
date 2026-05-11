import sys
import os
import argparse
import pandas as pd
import inspect
from collections import defaultdict
import runpy
import traceback
from pathlib import Path


TARGET_MAP = {
    'inv': 'inv',
    'det': 'det',
    'trace': 'trace',
    'dot': 'dot',
    'matmul': 'matmul',
    '@': 'matmul'
}

class DynamicAnalyzer:
    def __init__(self, targets):
        self.targets = targets
        self.hits = defaultdict(int)
        self.executions = defaultdict(int)
        self.metadata = {}
        
        self.profile_to_target = {}
        for t in targets:
            profile_name = TARGET_MAP.get(t, t)
            self.profile_to_target[profile_name] = t

    def _get_key(self, code):
        return (code.co_filename, code.co_name, code.co_firstlineno)

    def profile_func(self, frame, event, arg):
        if event == 'call':
            code = frame.f_code
            key = self._get_key(code)
            self.executions[key] += 1
            
            if key not in self.metadata:
                self.metadata[key] = self._get_info(code, frame)
            
            name = code.co_name
            if name in self.profile_to_target:
                target = self.profile_to_target[name]
                parent = frame.f_back
                if parent:
                    p_key = self._get_key(parent.f_code)
                    self.hits[(p_key, target)] += 1

        elif event == 'c_call':
            c_name = arg.__name__ if hasattr(arg, '__name__') else str(arg)
            if ' ' in c_name:
                c_name = c_name.split(' ')[-1].rstrip('>')
            
            if c_name in self.profile_to_target:
                target = self.profile_to_target[c_name]
                p_key = self._get_key(frame.f_code)
                self.hits[(p_key, target)] += 1
        
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
    parser = argparse.ArgumentParser(description="Dynamic Fan-In Analyzer")
    parser.add_argument("-m", "--module", action="store_true", help="Run library module as a script")
    parser.add_argument("script", help="Script path to run or module name if -m is used")
    parser.add_argument("script_args", nargs=argparse.REMAINDER, help="Arguments for the script")
    parser.add_argument("--targets", nargs="+", default=['matmul', 'dot', 'inv', 'det', 'trace', '@'],
                        help="Targets to track (default: matmul dot inv det trace @)")
    parser.add_argument("--output", default="Dynamic_FanIn.xlsx", help="Output Excel filename")
    parser.add_argument("--exclude-libs", action="store_true", help="Exclude calls from site-packages/stdlib")
    
    args = parser.parse_args()
    if not args.script:
        parser.print_help()
        return

    analyzer = DynamicAnalyzer(args.targets)
    
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
    for (key, target), call_count in analyzer.hits.items():
        meta = analyzer.metadata.get(key, {})
        file_path = meta.get('file', 'unknown')
        
        if args.exclude_libs:
            if 'site-packages' in file_path or 'lib/python' in file_path:
                continue

        mod = meta.get('module', 'unknown')
        cls = meta.get('class')
        func = meta.get('name', 'unknown')
        
        method_path = f"{mod}.{cls}.{func}" if cls else f"{mod}.{func}"
        
        results.append({
            'Method Path': method_path,
            'Target': target,
            'Call Count': call_count,
            'Total Executions': analyzer.executions.get(key, 0),
            'Function Size': meta.get('size', 0),
            'File Path': file_path
        })

    if results:
        df = pd.DataFrame(results)
        df = df.sort_values(by='Call Count', ascending=False)
        df.to_excel(args.output, index=False)


if __name__ == "__main__":
    main()
