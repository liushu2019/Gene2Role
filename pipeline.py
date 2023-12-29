import subprocess
import argparse
import os
import pandas as pd
from logging import getLogger, INFO, DEBUG
import logging
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing
import copy

def run_script(script_name, script_args):
    script_args = [str(x) for x in script_args]
    subprocess.run(["python", script_name] + script_args, check=True)

def exec_signeds2v(args):
    if args.input == '':
        raise ValueError("[SignedS2V] input argument needed.")
    current_dir = os.getcwd()   
    # Change the working directory to SignedS2V directory
    new_dir = f"{current_dir}/tools/SignedS2V/"
    os.chdir(new_dir)
    print ("Run SignedS2V...")
    logger.info("Run SignedS2V...")
    # if args.output is None:
    args.output = os.path.join(os.path.dirname(args.input), args.project+'.emb')
    signeds2v_args = ['--input', args.input, '--output', args.output, '--dimensions', args.dimensions, '--walk-length', args.walk_length, '--num-walks', args.num_walks, '--window-size', args.window_size, '--until-layer', args.until_layer, '--iter', args.iter, '--workers', args.workers]
    if args.OPT1:
        signeds2v_args = signeds2v_args + ['--OPT1']
    if args.OPT2:
        signeds2v_args = signeds2v_args + ['--OPT2']
    if args.OPT3:
        signeds2v_args = signeds2v_args + ['--OPT3']
    if args.scalefree:
        signeds2v_args = signeds2v_args + ['--scalefree']
    run_script("src/main.py", signeds2v_args)
    os.chdir(current_dir)
    print ("Finished SignedS2V.")
    logger.info("Finished SignedS2V.")
    return args.output

def exec_split_cells(args):
    if (args.matrix is None) or (args.cell_metadata is None):
        raise ValueError("[split cell] matrix and cell_metadata arguments needed.")
    print ("Run split cell...")
    logger.info("Run split cell...")
    split_cell_args = [args.matrix, args.cell_metadata]
    run_script("codes/split_cells.py", split_cell_args)
    print ("Finished split cell.")
    logger.info("Finished split cell.")
    
def exec_spearman(args):
    if (args.matrix is None) or (args.project is None):
        raise ValueError("[spearman] matrix and project arguments needed.")
    print ("Run spearman...")
    logger.info("Run spearman...")
    sub_args = [args.matrix, args.project, '-correlation_threshold', args.correlation_threshold]
    if args.reindex:
        sub_args = sub_args + ['--reindex']
    run_script("codes/spearman.py", sub_args)
    print ("Finished spearman.")
    logger.info("Finished spearman.")
    
def exec_eeisp(args):
    if (args.matrix is None) or (args.project is None):
        raise ValueError("[eeisp] matrix and project arguments needed.")
    print ("Run eeisp...")
    logger.info("Run eeisp...")
    sub_args = [args.matrix, args.project, args.mode, '--threCDI', args.threCDI, '--threEEI', args.threEEI]
    if args.reindex:
        sub_args = sub_args + ['--reindex']
    run_script("codes/eeisp.py", sub_args)
    print ("Finished eeisp.")
    logger.info("Finished eeisp.")
    
def exec_merge_edgelist(args, file_list):
    if (args.matrix is None) or (args.project is None):
        raise ValueError("[merge edgelist] matrix and project arguments needed.")
    print ("Run merge edgelist...")
    logger.info("Run merge edgelist...")
    sub_args = [args.matrix, args.project] + file_list
    run_script("codes/merge_edgelist.py", sub_args)
    print ("Finished merge edgelist.")
    logger.info("Finished merge edgelist.")
    
def exec_(args):
    if (args.matrix is None) or (args.cell_metadata is None):
        raise ValueError("[split cell] matrix and cell_metadata arguments needed.")
    print ("Run ...")
    logger.info("Run ...")
    sub_args = [args.matrix, args.cell_metadata]
    run_script("codes/.py", sub_args)
    print ("Finished .")
    logger.info("Finished .")
    
def main():
    parser = argparse.ArgumentParser(description="Gene2Role pipeline script")

    parser.add_argument("TaskMode", type=int, default=1, help="Task mode. 1: run SignedS2V for an edgelist file. 2: run spearman and SignedS2V from gene X cell count matrix. 3: run eeisp and SignedS2V from gene X cell count matrix.")
    parser.add_argument("CellType", type=int, default=1, help="Cell type. 1: single cell-type. 2: multiple cell-type in one matrix. 3: multiple gene X cell count matrixes.")
    parser.add_argument("EmbeddingMode", type=int, default=1, help="Embedding mode. 1: single network embedding. 2: multiple network embedding, only work if the previous argument is 2.")
    parser.add_argument('input', nargs='+', default='[ALL] Input file. For CellType = 1, an edgelist; for CellType = 2, a gene X cell matrix; for CellType = 3, multiple gene X cell matrixes separated by space.')
    parser.add_argument("--project", type=str, default='sample', help="Project name which will be used as the folder name.")
    parser.add_argument('--reindex', action='store_true', 
                        help='Flag for reindexing genes. Needed if gene name contained in input file.')
    # args for SignedS2V
    parser.add_argument('--output', nargs='?', default=None,
                        help='[SignedS2V] Output emb path, if Not given, follow input file name')
    parser.add_argument('--dimensions', type=int, default=128,
                        help='[SignedS2V] Number of dimensions. Default is 128.')
    parser.add_argument('--walk-length', type=int, default=80,
                        help='[SignedS2V] Length of walk per source. Default is 80.')
    parser.add_argument('--num-walks', type=int, default=100,
                        help='[SignedS2V] Number of walks per source. Default is 100.')
    parser.add_argument('--window-size', type=int, default=10,
                        help='[SignedS2V] Context size for optimization. Default is 10.')
    parser.add_argument('--until-layer', type=int, default=3,
                        help='[SignedS2V] Calculation until the layer. Default is 3.')
    parser.add_argument('--iter', default=5, type=int,
                        help='[SignedS2V] Number of epochs in SGD')
    parser.add_argument('--workers', type=int, default=8,
                        help='Number of parallel workers. Default is 8.')
    parser.add_argument('--OPT1', action='store_true', 
                        help='optimization 1')
    parser.add_argument('--OPT2', action='store_true',
                        help='optimization 2')
    parser.add_argument('--OPT3', action='store_true',
                        help='optimization 3')
    parser.add_argument('--scalefree', action='store_true',
                        help='scale free flag')
    # spearman
    # parser.add_argument("--file_path", 
                        # type=str, 
                        # help="[spearman] Path to the CSV file containing gene expression data.")
    # parser.add_argument("--out_file_name",  #folder / file name / project name
    #                     type=str, 
    #                     help="[spearman] out put file name")
    parser.add_argument("--correlation_threshold", 
                        type=float, 
                        help="[spearman] Threshold for filtering high correlations. (default:0.4)",
                        default=0.4)
    # essip
    parser.add_argument("--matrix", help="[ESSIP/spearman/split cell] CSV file for gene-cell matrix.", type=str)
    # parser.add_argument("--filename", help="[ESSIP] File name", type=str)
    parser.add_argument("--mode", help="[ESSIP] mode for (1) percentage threshold or (2) solid threshold for CDI and EEI (default: 1)", type=int, default=1)
    parser.add_argument("--threCDI", help="[ESSIP] Threshold for CDI (default: 0.5). When mode = 1, filter out top value*100\%; when mode = 1, filter out those >= value.", type=float, default=0.5)
    parser.add_argument("--threEEI", help="[ESSIP] Threshold for EEI (default: 0.5). When mode = 1, filter out top value*100\%; when mode = 1, filter out those >= value.", type=float, default=0.5)
    # split cell
    # parser.add_argument("--count_matrix", type=str, help="[split cell] CSV file for gene-cell matrix.")
    parser.add_argument("--cell_metadata", type=str, help="[split cell] CSV file with cell type information.", default=None)        
    # downstream TODO    
    args = parser.parse_args()
    print (args)
    if args.CellType == 3 and len(args.input) <= 0:
        raise ValueError("Received one input file while CellType = 3.")
    if args.TaskMode == 1 and args.EmbeddingMode == 2 and len(args.input) <= 0:
        raise ValueError("Received one input file while TaskMode = 1 and EmbeddingMode = 2.")
    if args.CellType == 2 and args.cell_metadata is None:
        raise ValueError("No cell metadata received while CellType = 2 or 3.")
    if args.CellType == 3 or (args.TaskMode == 1 and args.EmbeddingMode == 2): # multiple input
        absolute_directory = os.path.abspath(args.input[0])
        args.matrix = args.input = [os.path.abspath(x) for x in args.input]
    else: # single input
        args.input = args.input[0]
        absolute_directory = os.path.abspath(args.input)
        args.matrix = args.input = absolute_directory
        
    output_directory = os.path.dirname(absolute_directory)
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    file_handler = logging.FileHandler(f'{output_directory}/Gene2Role_{args.project}.log', mode="w")
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
    logger.setLevel(DEBUG)
    logger.info(args)
    logger.info(f"Absolute directory: {absolute_directory}")
    output_file_name_list = []
    
    if args.TaskMode == 1: #run SignedS2V
        if args.EmbeddingMode == 1:
            output_file_name_list.append(exec_signeds2v(args))
        else: # multiple edgelist
            args.matrix = args.matrix[0]
            exec_merge_edgelist(args, args.input)
            args.input = os.path.join(os.path.dirname(args.matrix), args.project + "_merged.edgelist")
            output_file_name_list.append(exec_signeds2v(args)) 
    else:
        if args.CellType == 2: #multi cell
            exec_split_cells(args)
            output_directory = os.path.dirname(args.matrix)
            output_directory = os.path.join(output_directory, 'splitMatrix')
            index_tracker = pd.read_csv(os.path.join(output_directory, 'index_tracker.tsv'), sep='\t', index_col=0)
            mapping_file_dir = os.path.join(output_directory, 'index_tracker.tsv')
            list_input_file = []
            list_project = []
            project_org = args.project
            with ProcessPoolExecutor(max_workers=args.workers) as executor:
                for cell_type in index_tracker.columns:
                    args_tmp = copy.deepcopy(args)
                    output_file = os.path.join(output_directory, f'{cell_type}.csv')
                    args_tmp.matrix = output_file
                    args_tmp.project = project_org + '_' + cell_type
                    list_project.append(project_org + '_' + cell_type)
                    if args.TaskMode == 2: #run spearman and SignedS2V from gene X cell count matrix
                        job = executor.submit(exec_spearman, args_tmp)
                        list_input_file.append(os.path.join(os.path.dirname(args_tmp.matrix), args_tmp.project + "_spearman.edgelist"))
                    elif args.TaskMode == 3: #run eeisp and SignedS2V from gene X cell count matrix.
                        job = executor.submit(exec_eeisp, args_tmp)
                        list_input_file.append(os.path.join(os.path.dirname(args_tmp.matrix), args_tmp.project + "_eeisp.edgelist"))
                        
            if args.EmbeddingMode == 2: #multi embedding
                args.project = project_org
                exec_merge_edgelist(args, list_input_file)
                args.input = os.path.join(os.path.dirname(args.matrix), args.project + "_merged.edgelist")
                output_file_name_list.append(exec_signeds2v(args))
            else: # embed separatly
                for input_dir in list(zip(list_input_file, list_project)):
                    args.input = input_dir[0]
                    args.project = input_dir[1]
                    output_file_name_list.append(exec_signeds2v(args))
        elif args.CellType == 3: #multiple matrixes
            input_list_org = args.input
            project_org = args.project
            matrix_org = args.matrix
            list_input_file = []
            list_project = []
            with ProcessPoolExecutor(max_workers=args.workers) as executor:
                for ix, input_dir in enumerate(input_list_org):
                    args_tmp = copy.deepcopy(args)
                    args_tmp.input = input_dir
                    args_tmp.matrix = input_dir
                    args_tmp.project = project_org + '_' + str(ix)
                    list_project.append(project_org + '_' + str(ix))
                    if args.TaskMode == 2: #run spearman and SignedS2V from gene X cell count matrix
                        job = executor.submit(exec_spearman, args_tmp)
                        list_input_file.append(os.path.join(os.path.dirname(args_tmp.matrix), args_tmp.project + "_spearman.edgelist"))
                    elif args.TaskMode == 3: #run eeisp and SignedS2V from gene X cell count matrix.
                        job = executor.submit(exec_eeisp, args_tmp)
                        list_input_file.append(os.path.join(os.path.dirname(args_tmp.matrix), args_tmp.project + "_eeisp.edgelist"))
                        
            if args.EmbeddingMode == 2: #multi embedding
                args.project = project_org
                args.matrix = args.matrix[0]
                exec_merge_edgelist(args, list_input_file)
                args.input = os.path.join(os.path.dirname(input_list_org[0]), args.project + "_merged.edgelist")
                args.matrix = matrix_org
                output_file_name_list.append(exec_signeds2v(args))
            else: # embed separatly
                for input_dir in list(zip(list_input_file, list_project)):
                    args.input = input_dir[0]
                    args.project = input_dir[1]
                    output_file_name_list.append(exec_signeds2v(args))
                    
        elif args.CellType == 1: #single cell
            if args.TaskMode == 2: #run spearman and SignedS2V from gene X cell count matrix
                exec_spearman(args)
                args.input = os.path.join(os.path.dirname(args.matrix), args.project + "_spearman.edgelist")
                output_file_name_list.append(exec_signeds2v(args))
                mapping_file_dir =  os.path.join(os.path.dirname(args.matrix), args.project + "_spearman_nodeID_mapping.tsv")
            elif args.TaskMode == 3: #run eeisp and SignedS2V from gene X cell count matrix.
                exec_eeisp(args)
                args.input = os.path.join(os.path.dirname(args.matrix), args.project + "_eeisp.edgelist")
                output_file_name_list.append(exec_signeds2v(args))
                mapping_file_dir =  os.path.join(os.path.dirname(args.matrix), args.project + "_number_nonzero_exp.txt")
            
    print("Pipeline completed.")

    print ("-------------------OUTPUT EMBEDDING-------------------")
    for output_file_name in output_file_name_list:
        print (f"{output_file_name}")
    if args.TaskMode == 2 or args.reindex:
        print ("------------------- GENE NAME INFO -------------------")
        print (mapping_file_dir)
    print ("------------------------------------------------------")
    
if __name__ == "__main__":
    logger = getLogger("Gen2Role_pipeline")
    main()