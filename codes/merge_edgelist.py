import argparse
import os
import pandas as pd
import sys

def main(*args):
    if len(args) < 3:
        raise ValueError("[merge edgelist] matrix, project and file list arguments needed.")
    matrix = args[0]
    project = args[1]
    file_list = args[2:]
    list_df = []
    for file_dir in file_list:
        if not os.path.isabs(file_dir):
            if os.path.exists(os.path.join(os.getcwd(), file_list)):
                file_dir = os.path.exists(os.path.join(os.getcwd(), file_list))
                print (f"FILE directory is NOT absolute, try to find file in {dir} <--- CONFIRM PLZ")
                pass
            elif os.path.exists(os.path.join(os.path.dirname(matrix), file_list)):
                file_dir = os.path.join(os.path.dirname(matrix), file_list)
                print (f"FILE directory is NOT absolute, try to find file in {dir} <--- CONFIRM PLZ")
                pass 
            elif os.path.exists(os.path.join(os.path.dirname(matrix), project, file_list)):
                file_dir = os.path.join(os.path.dirname(matrix), project, file_list)
                print (f"FILE directory is NOT absolute, try to find file in {dir} <--- CONFIRM PLZ")
                pass
            else:
                raise FileNotFoundError()
        df = pd.read_csv(file_dir, header=None, sep='\t')
        list_df.append(df)
    print ("Merge edgelist succeed. Len=[{}]".format(','.join([str(df.shape[0]) for df in list_df])))
    output_dir = os.path.join(os.path.dirname(matrix), project + "_merged.edgelist")
    print (f"Output file is {output_dir}")
    pd.concat(list_df).reset_index(drop=True).to_csv(output_dir, sep='\t', header=None, index=None)
    
if __name__ == "__main__":
    arguments = sys.argv[1:]

    main(*arguments)