import pandas as pd
import os


def ReplaceTextInExcelColumn(file_path, old_text, new_text):
    df = pd.read_excel(file_path)
    df[column] = df[column].astype(str).str.replace(old_text, new_text, regex=False)
    df.to_excel(file_path, index=False)

if __name__ == "__main__":
    file_path = "Detailed_FanIn.xlsx"
    # file_path = "Dynamic_FanIn.xlsx"
    column = "File Path"
    old_text = "/Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/"
    new_text = "https://github.com/MuyeedAhmed/scipy/tree/main/"
    
    ReplaceTextInExcelColumn(file_path, old_text, new_text)