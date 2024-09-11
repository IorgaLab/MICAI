import argparse
from analysis import analyze

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument( "--file", default = None, type=str, required=True, help="Path (relative or absolute) of the file to analyze")
    parser.add_argument( "--species", default = None, type=str, required=True, help="Species (Kp,Ec,Pa)",choices=["Kp","Ec","Pa"])
    parser.add_argument( "--mode", default = None, type=str, required=True, help="Mode : binary or regression", choices=["binary","regression"])
    parser.add_argument( "--encoding", default = None, type=str, required=True, help="Encoding : nucl or prot", choices=["nucl","prot"])
    parser.add_argument( "--size", default=None,type=int,help="Size of n-gram/k-mer analysis, not mandatory",choices=[6,8,11,14,17,18,24,33,42,51])
    parser.add_argument( "--format", default="fasta", type=str, help="Format of the given file, fasta by default")
    args = parser.parse_args()

    results=analyze(args.file,args.species,args.mode,args.encoding,args.size,args.format)

    for k in results.keys():
        print(k,":",results[k])