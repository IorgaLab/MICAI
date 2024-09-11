import numpy as np
import pickle
import pandas as pd
import xgboost as xgb
from Bio import SeqIO
import json


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

    rep="data/"+species+"_"+str(N)+"-"+encoding+"_"+mode+"/"


    records = list(SeqIO.parse(file, format))
    s=set()
    for seq_line in records:
        seq=str(seq_line.seq)
        n=len(seq)
        for mot in range(n-N+1):
            s.add(seq[mot:mot+N])

    d_len_pos=load_json(rep+"len_pos.json")
    name_to_new_pos=load_json(rep+"name_to_new_pos.json")

    tab=np.zeros(d_len_pos["len_pos"])
    for elem in s:
        val=name_to_new_pos.get(elem)
        if val!=None:
            tab[val]=1
 
    if np.mean(tab)<0.01:
       print(file,np.mean(tab))
       raise Exception("The matrix is empty (or almost), are you sure it is the right format?")

    dico_names=load_json(rep+"dico_names.json")
    dico_alpha=load_json(rep+"dico_alpha.json")

    if mode=="regression":
        group=load_json(rep+"groupage.json")

    if mode=="binary":
        dico_treshold=load_json(rep+"dico_threshold.json")

    dico_results={}

    for antibio in dico_names.keys():

        equiv_features=load_json(rep+"equiv_features_"+antibio+".json")

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
        model.load_model(rep+"model_"+antibio+".json")

        if mode=="binary":
            y_pred=model.inplace_predict(pd.DataFrame(array))[0]>float(dico_treshold[antibio])
            if y_pred:
                dico_results[dico_names[antibio]]="Resistant"
            else:
                dico_results[dico_names[antibio]]="Sensible"
        
        if mode=="regression":
            with open(rep+"label_encoder_"+antibio,"rb") as f:
                le=pickle.load(f)
            y_pred=np.clip(round((model.inplace_predict(pd.DataFrame(array)))[0]),0,len(le.classes_)-1)
            y_le=le.inverse_transform([y_pred])[0]

            y_final=False
            for k in group[antibio]:
                if y_le in k:
                    dico_results[dico_names[antibio]]="From "+remove_point(2.**min(k))+" to "+remove_point(2.**max(k))
                    y_final=True
            if y_final==False:
                dico_results[dico_names[antibio]]=remove_point(2.**y_le)

    return(dico_results)