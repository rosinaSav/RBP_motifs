import argparse
import itertools as it
import multiprocessing as mp
import os
import re
import subprocess
import sys

def flatten(structured_list):
    '''
    Flatten a structured list.
    '''
    flat_list = list(it.chain(*structured_list))
    return(flat_list)

def get_extension(file_name, extension_length, valid_list = None):
    '''
    Determine the extension at the end of a file name.
    '''
    extension = file_name[-extension_length:]
    if valid_list:
        if extension not in valid_list:
            print("File format must be included in {0}!".format(valid_list))
            sys.exit()
    return(extension)

def line_count(file):
    '''
    Count the number of lines in a file.
    '''
    #not using wc -l because I want the number of lines, not the number of newlines.
    output = run_process(["grep", "-c", "^", file])
    return(int(output))

def list_to_dict(input_list, index1, index2, as_list = False, uniquify = False, floatify = False):
    '''
    Convert the input_list into a dictionary, with the index1th element of each sublist as the key and the index2th element as the value.
    '''
    if as_list and floatify:
        print("_as_list_ and _floatify_ can't both be True!")
        raise Exception
    output_dict = {}
    for i in input_list:
        if not as_list:
            if floatify:
                #convert the value into a float
                output_dict[i[index1]] = float(i[index2])
            else:
                output_dict[i[index1]] = i[index2]
        else:
            #if several sublists can have the same value in index1 and you want all their index2th values as a list
            if i[index1] not in output_dict:
                output_dict[i[index1]] = []
            output_dict[i[index1]].append(i[index2])
    if as_list and uniquify:
        #if the values are lists and you don't want duplicates in the value lists
        output_dict = {i: sorted(list(set(output_dict[i]))) for i in output_dict}
    return(output_dict)

def make_dir(dir_name):
    '''
    Check whether a directory exists and if not, create it.
    '''
    if not os.path.exists(dir_name):
        os.mkdir(dir_name)   

def overlap(a, b):
    '''
    Given two lists, determine which elements appear in both.
    '''
    set_a = set(a)
    intersection = set_a.intersection(b)
    return(list(intersection))

def parse_arguments(description, arguments, floats = None, flags = None, ints = None):
    '''
    Use argparse to parse a set of input arguments from the command line.
    '''
    if not floats:
        floats = []
    if not flags:
        flags = []
    if not ints:
        ints = []
    parser = argparse.ArgumentParser(description = description)
    for pos, argument in enumerate(arguments):
        if pos in flags:
            parser.add_argument("--{0}".format(argument), action = "store_true", help = argument)
        else:
            if pos in floats:
                curr_type = float
            elif pos in ints:
                curr_type = int
            else:
                curr_type = str
            parser.add_argument(argument, type = curr_type, help = argument)
    args = parser.parse_args()
    return(args)
    
def print_elements(input_list):
    '''
    Take a list and print out the elements separated by carriage returns.
    '''
    for i in input_list:
        print(i)
    print("\n")

def remove_file(file_name):
    '''
    Remove a file, if it exists.
    '''
    try:
        os.remove(file_name)
    except FileNotFoundError:
        pass

def run_in_parallel(input_list, args, func, kwargs_dict = None, workers = None, onebyone = False):
    '''
    Take an input list, divide into chunks and then apply a function to each of the chunks in parallel.
    '''
    if not workers:
        #divide by two to get the number of physical cores
        #subtract one to leave one core free
        workers = int(os.cpu_count()/2 - 1)
    elif workers == "all":
        workers = os.cpu_count()
    arg_to_parallelize = args.index("foo")
    if not onebyone:
        chunk_list = [input_list[i::workers] for i in range(workers)]
    else:
        #each element in the input list will constitute a chunk of its own.
        chunk_list = input_list
    pool = mp.Pool(workers)
    results = []
    for i in chunk_list:
        current_args = args.copy()
        current_args[arg_to_parallelize] = i
        if kwargs_dict:
            process = pool.apply_async(func, tuple(current_args), kwargs_dict)
        else:
            process = pool.apply_async(func, tuple(current_args))            
        results.append(process)
    pool.close()
    pool.join()
    return(results)

def run_process(arguments, return_string = True, input_to_pipe = None, file_for_input = None, file_for_output = None, univ_nl = True, shell = False):
    '''
    Run a command on the command line.
    '''
    if file_for_input:
        input_file = open(file_for_input)
        stdin_src = input_file
    else:
        stdin_src = subprocess.PIPE
    if file_for_output:
        output_file = open(file_for_output, "w")
        stdout_dest = output_file
    else:
        stdout_dest = subprocess.PIPE
    arguments = [str(i) for i in arguments]
    if shell:
        arguments = " ".join(arguments)
    process = subprocess.Popen(arguments, shell = shell, stdout = stdout_dest, stderr = subprocess.PIPE,
                               stdin = stdin_src, universal_newlines = univ_nl)
    if input_to_pipe:
        stdout, stderr = process.communicate(input_to_pipe)
    else:
        stdout, stderr = process.communicate()
    return_code = process.poll()
    if return_code != 0:
        print("Process failed!")
        print(arguments)
        print(stderr)
        if file_for_input:
            input_file.close()
        if file_for_output:
            output_file.close()
        return("error")
    if file_for_input:
        input_file.close()
    if file_for_output:
        output_file.close()
    #if the process returns bytes but you want to get a string back.
    if return_string and type(stdout) == bytes:
        stdout = stdout.decode("utf-8")
    return(stdout)
    
def write_to_fasta(identifiers,sequences,file_name):
    '''
    Write a set of sequences to fasta.
    '''
    check_equal_lengths([identifiers,sequences])
    fasta_file = open(file_name,"w")
    for i,j in enumerate(identifiers):
        fasta_file.write(">{0}\n".format(j))
        fasta_file.write("{0}\n".format(sequences[i]))
    fasta_file.close()
