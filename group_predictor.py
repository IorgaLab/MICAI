import pandas as pd
import argparse
import sys
import os
from analysis import analyze
import multiprocessing
from functools import partial


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument( "--path", default = None, type=str, required=True, help="Path of the folder")
    parser.add_argument( "--species", default = None, type=str, required=True, choices=["Kp","Ec","Pa"], help="Espèce (Kp,Ec,Pa)")
    parser.add_argument( "--mode", default = None, type=str, required=True, help="Mode : binary or regression", choices=["binary","regression"])
    parser.add_argument( "--output", default = None, type=str, required=True, help="Name of the csv file to store the data in")
    parser.add_argument( "--size", default = None, type=int, help="Size of n-grams chosen")
    parser.add_argument( "--encoding", default = None, type=str, required=True, help="Encoding : nucl or prot", choices=["nucl","prot"])    
    parser.add_argument( "--format", default="fasta", type=str, help="Format of the given file, fasta by default")
    args = parser.parse_args()


    tmp=os.listdir(args.path)
    list_path=[args.path+"/"+t for t in tmp]

    list_path.sort()
    list_names=[]

    for path in list_path:
        # Remove extension, take only the last part -> name
        list_names.append(path.split(".")[0].split("/")[-1])
    
    # Multiprocessing analysis
    partial_process=partial(analyze,species=args.species,mode=args.mode,encoding=args.encoding,size=args.size,format=args.format)
    with multiprocessing.Pool() as pool:
        results=pool.map(partial_process,list_path)

    names_antibio=results[0].keys()

    # Fill the results into a csv file
    dataframe=pd.DataFrame(columns=names_antibio,index=list_names)

    for k in range(len(results)):
        for antibio in names_antibio:
            dataframe.loc[list_names[k],antibio]=results[k][antibio]

    dataframe.to_csv(args.output)

if __name__ == "__main__":
    main()