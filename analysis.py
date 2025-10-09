import numpy as np
import pickle
import pandas as pd
import xgboost as xgb
from Bio import SeqIO
import json
import os

def remove_point(x):
    if x>=1:
        return(str(int(x)))
    else:
        return(str(x))
    
def load_json(file):
    with open(file,"r") as f:
        dico=json.load(f)
    return(dico)

def analyze(file,species,mode,encoding,size,format):

    if size==None:
        if encoding=="nucl":
            N=24
        if encoding=="prot":
            N=8
    else:
        N=size

    if encoding=="prot":
        encoding="grams"
    if encoding=="nucl":
        encoding="mers"

    rep = os.path.abspath(
        os.path.join(os.path.dirname(__file__), "data", f"{species}_{N}-{encoding}_{mode}")
    )

    records = list(SeqIO.parse(file, format))
    s=set()
    for seq_line in records:
        seq=str(seq_line.seq)
        n=len(seq)
        for mot in range(n-N+1):
            s.add(seq[mot:mot+N])

    d_len_pos=load_json(rep+"/len_pos.json")
    name_to_new_pos=load_json(rep+"/name_to_new_pos.json")

    tab=np.zeros(d_len_pos["len_pos"])
    for elem in s:
        val=name_to_new_pos.get(elem)
        if val!=None:
            tab[val]=1
 
    if np.mean(tab)<0.01:
        print(np.mean(tab))
        print("The feature matrix is empty (or almost), are you sure it is the right encoding format? Nucleotides? Proteins?")

    dico_names=load_json(rep+"/dico_names.json")
    dico_alpha=load_json(rep+"/dico_alpha.json")

    if mode=="binary":
        dico_treshold=load_json(rep+"/dico_threshold.json")

    dico_results={}

    for antibio in dico_names.keys():

        equiv_features=load_json(rep+"/equiv_features_"+antibio+".json")

        array=np.zeros((1,len(equiv_features)),dtype=np.uint8)
        for j, cols in enumerate(equiv_features):
            if len(cols) > 1:
                new_cols=[int(k) for k in cols]
                non_zero_counts = np.count_nonzero(tab[new_cols])
                if non_zero_counts / len(new_cols) >= float(dico_alpha[antibio]):
                    array[0, j] = 1
                else:
                    array[0, j] = 0
            elif len(cols)==1:
                array[0, j] = tab[int(cols[0])]
            else:
                pass


        model=xgb.Booster()
        model.load_model(rep+"/model_"+antibio+".json")

        if mode=="binary":
            y_pred=model.inplace_predict(pd.DataFrame(array))[0]>float(dico_treshold[antibio])
            if y_pred:
                dico_results[dico_names[antibio]]="Resistant"
            else:
                dico_results[dico_names[antibio]]="Sensible"
        
        if mode=="regression":  
            with open(rep + f"/dillution_range_{antibio}.json", "r") as f:
                dillution_range = json.load(f)
                
            y_pred=np.clip(round((model.inplace_predict(pd.DataFrame(array)))[0]), np.log2(dillution_range["minimum"]) , np.log2(dillution_range["maximum"])+1)
            y_pred = 2.**y_pred
            
            if y_pred == dillution_range["minimum"]:
                sign = "<="
            elif y_pred == dillution_range["maximum"]*2:
                sign = ">"
                y_pred = dillution_range["maximum"]
            else:
                sign = ""
            
            prediction = remove_point(y_pred)
            dico_results[dico_names[antibio]]=sign+prediction

    return(dico_results)