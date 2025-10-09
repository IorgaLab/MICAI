import pandas as pd
import argparse
import sys
import os
from analysis import analyze
import multiprocessing
from functools import partial


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument( "--type", default = None, type=str, required=True, choices=["file","folder"], help="file OR folder")
    parser.add_argument( "--path", default = None, type=str, required=True, help="Name of the file/folder")
    parser.add_argument( "--species", default = None, type=str, required=True, choices=["Kp","Ec","Pa"], help="Espèce (Kp,Ec,Pa)")
    parser.add_argument( "--mode", default = None, type=str, required=True, help="Mode : binary or regression", choices=["binary","regression"])
    parser.add_argument( "--output", default = None, type=str, required=True, help="Name of the csv file to store the data in")
    parser.add_argument( "--size", default=None,type=int,help="Size of n-gram/k-mer analysis, default =8 for prot and =24 for nucl",choices=[8,11,14,24,33,42])
    parser.add_argument( "--encoding", default = None, type=str, required=True, help="Encoding : nucl or prot", choices=["nucl","prot"])    
    args = parser.parse_args()

    if args.type=="file":
        #TODO
        sys.exit()
    elif args.type=="folder":
        tmp=os.listdir(args.path)
        list_path=[args.path+"/"+t for t in tmp]
    else:
        print("Error")
    list_path.sort()
    list_names=[]

    for path in list_path:
        list_names.append(path.split(".")[0].split("/")[-1])
    

    partial_process=partial(analyze,species=args.species,mode=args.mode,encoding=args.encoding,size=args.size,format="fasta")
    with multiprocessing.Pool() as pool:
        results=pool.map(partial_process,list_path)

    names_antibio=results[0].keys()

    dataframe=pd.DataFrame(columns=names_antibio,index=list_names)

    for k in range(len(results)):
        for antibio in names_antibio:
            dataframe.loc[list_names[k],antibio]=results[k][antibio]

    dataframe.to_csv(args.output)

if __name__ == "__main__":
    main()