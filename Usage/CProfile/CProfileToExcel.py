import argparse
import cProfile
import pstats
import pandas as pd
import sys
import os
import pytest


def main():
    parser = argparse.ArgumentParser(description="Run cProfile on a pytest suite and export results to Excel")
    parser.add_argument("--output", default="cprofile_result.xlsx", help="Output Excel filename")
    parser.add_argument("--package", nargs='+', help="The package to test")
    parser.add_argument("remaining", nargs=argparse.REMAINDER, help="Additional pytest arguments")
    args = parser.parse_args()

    pytest_commands = ["--pyargs"] + args.package + args.remaining
    
    xlsx_output = args.output
    
    profiler = cProfile.Profile()
    profiler.enable()
    try:
        exit_code = pytest.main(pytest_commands)
        print(f"Pytest finished with exit code {exit_code}")
    except Exception as e:
        print(f"An error occurred while running pytest: {e}")
    finally:
        profiler.disable()
    
    ps = pstats.Stats(profiler)
    
    data = []
    for func, stat in ps.stats.items():
        primitive_calls = stat[0]
        total_calls = stat[1]
        total_time = stat[2]
        cumulative_time = stat[3]
        
        res = {
            "File Path": func[0],
            "Line": func[1],
            "Function": func[2],
            "Total Calls": total_calls,
            "Primitive Calls": primitive_calls,
            "Total Time (s)": total_time,
            "Time per Call (s)": total_time / total_calls if total_calls > 0 else 0,
            "Cumulative Time (s)": cumulative_time,
            "Cumulative Time per Call (s)": cumulative_time / total_calls if total_calls > 0 else 0
        }
        data.append(res)
    
    df = pd.DataFrame(data)
    
    df = df.sort_values(by="Total Calls", ascending=False)
    df.to_excel(xlsx_output, index=False)
    

if __name__ == "__main__":
    main()
