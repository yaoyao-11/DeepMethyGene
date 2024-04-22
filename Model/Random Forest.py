from sklearn.svm import SVR
from sklearn.metrics import mean_squared_error
import pandas as pd
from sklearn.model_selection import KFold
from scipy.stats import pearsonr
import pandas as pd
import math 
from tqdm import tqdm
import numpy as np
from sklearn.model_selection import train_test_split
import os,logging,pickle,random,torch
from sklearn.linear_model import ElasticNet
from sklearn.model_selection import GridSearchCV
from sklearn.metrics import r2_score
# import matplotlib.pyplot as plt
from torch.autograd import Variable
import torch.nn.functional as F
import torch.nn as nn
from sklearn.model_selection import KFold
from collections import  Counter
import warnings
import itertools
import gc
import os
from torch.utils.data import DataLoader
from sklearn.ensemble import RandomForestRegressor
# ...

def cross_validation_random_forest(x, y, n_splits=5):
    kf = KFold(n_splits=n_splits, shuffle=True, random_state=42)
    r2_scores = []
   
    for train_idx, val_idx in kf.split(x):
        x_train, x_val = x[train_idx], x[val_idx]
        y_train, y_val = y[train_idx], y[val_idx]

        # 使用随机森林回归模型
        model = RandomForestRegressor(n_estimators=100, random_state=42)
        model.fit(x_train, y_train.ravel())

        # 验证模型
        y_pred = model.predict(x_val)
        r2 = r2_score(y_val, y_pred)
        r2_scores.append(r2)

    return np.mean(r2_scores)

# ...
# Load data
BRCA_exp=pd.read_csv("../need_gene_exp_13988(13982).csv")
BRCA_meth=pd.read_csv("../BRCA_data_meth.csv")
BRCA_meth_range=pd.read_csv("../BRCA_data_meth_range.csv")
promoter_range=pd.read_csv("../hg19_promoter.txt",sep="\t")
# gene_list=pd.read_csv("../gene_list.csv")

# gene_list=pd.read_csv("../remain_gene_onefold.csv")

BRCA_meth_index=BRCA_meth["Unnamed: 0"]
BRCA_meth=BRCA_meth.drop("Unnamed: 0",axis=1)
BRCA_meth.index=list(BRCA_meth_index)
BRCA_exp_index=list(BRCA_exp["Unnamed: 0"])
BRCA_exp=BRCA_exp.drop(["Unnamed: 0"],axis=1)
BRCA_exp.index=BRCA_exp_index
need_BRCA_exp_t=BRCA_exp.T
gene_list=pd.read_csv("../gene_list.csv")

for i in tqdm(list(gene_list["x"])[6767:7000]):
    gene_name=promoter_range.loc[promoter_range["gene name"]==i]
    need_BRCA_meth_range=BRCA_meth_range.loc[BRCA_meth_range["seqnames"]==list(gene_name["chrID"])[0]]
    need_BRCA_meth_range_index=need_BRCA_meth_range["Unnamed: 0"]
    need_BRCA_meth_range=need_BRCA_meth_range.drop("Unnamed: 0",axis=1)
    need_BRCA_meth_range.index=list(need_BRCA_meth_range_index)
    need_cpg={}
    for ii in range(len(need_BRCA_meth_range)):
        if int(need_BRCA_meth_range.iloc[ii,1])<=int(gene_name["end"]+1e7) and int(need_BRCA_meth_range.iloc[ii,1])>=int(gene_name["start"]-1e7):
            need_cpg[need_BRCA_meth_range.iloc[ii,5]]=need_BRCA_meth_range.iloc[ii,1]
    need_BRCA_meth=BRCA_meth.loc[list(need_cpg.keys())]
    need_BRCA_meth_T=need_BRCA_meth.T
    a=pd.DataFrame(need_BRCA_exp_t[i])
    # ...
    if len(a.columns) == 1:
        y = need_BRCA_exp_t[i]
        df = pd.concat([need_BRCA_meth_T, y], axis=1)
        data = df.sample(frac=1, random_state=1)
        train_x = data.drop([i], axis="columns").values
        train_y = data[i].values
        mean_r2 = cross_validation_random_forest(train_x, train_y)
        File = open("../result_r2/RandomForestnormal.txt", 'a')   
        File.write(i + "\t" + str(mean_r2) + "\n")
        File.flush()
    else:
        mean_=[]
        for len_ in range(len(a.columns)):
            y=a.iloc[:,len_]
            df = pd.concat([need_BRCA_meth_T,y],axis=1)
            data=df.sample(frac=1,random_state=1)
            train_x = data.drop([i], axis="columns").values
            train_y = data[i].values
            mean_r2 = cross_validation_random_forest(train_x, train_y)
            mean_.append(mean_r2)
        File = open("../result_r2/RandomForestnormal.txt", 'a')   
        File.write(i+"\t"+str(max(mean_))+"\n")
        File.flush()
        print(i,mean_)
        print(i+": ",max(mean_))

