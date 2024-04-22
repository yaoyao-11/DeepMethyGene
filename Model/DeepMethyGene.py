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

warnings.filterwarnings("ignore")
SEED = 42
random.seed(SEED)
np.random.seed(SEED)
torch.manual_seed(SEED)
torch.cuda.manual_seed_all(42)

# Load data
BRCA_exp=pd.read_csv("../need_gene_exp_13988(13982).csv")
BRCA_meth=pd.read_csv("../BRCA_data_meth.csv")
BRCA_meth_range=pd.read_csv("../BRCA_data_meth_range.csv")
promoter_range=pd.read_csv("../hg19_promoter.txt",sep="\t")
gene_list=pd.read_csv("../gene_list.csv")
# gene_list=pd.read_csv("./remain_data.csv")
# gene_list=pd.read_csv("../remain_gene_onefold.csv")

BRCA_meth_index=BRCA_meth["Unnamed: 0"]
BRCA_meth=BRCA_meth.drop("Unnamed: 0",axis=1)
BRCA_meth.index=list(BRCA_meth_index)
BRCA_exp_index=list(BRCA_exp["Unnamed: 0"])
BRCA_exp=BRCA_exp.drop(["Unnamed: 0"],axis=1)
BRCA_exp.index=BRCA_exp_index
need_BRCA_exp_t=BRCA_exp.T



torch.cuda.set_device(1)  # 设置当前使用的GPU索引为0
# torch.cuda.set_per_process_memory_fraction(0.6, 0)
device = torch.device('cuda:1' if torch.cuda.is_available() else 'cpu')
class ResidualBlock(nn.Module):
    def __init__(self, channels):
        super(ResidualBlock, self).__init__()
        self.conv1 = nn.Conv1d(channels, channels, kernel_size=3, padding=1)
        self.conv2 = nn.Conv1d(channels, channels, kernel_size=3, padding=1)
        self.leaky_relu = nn.LeakyReLU(0.01)

    def forward(self, x):
        residual = x
        out = self.leaky_relu(self.conv1(x))
        out = self.conv2(out)
        out += residual  # 添加残差
        out = self.leaky_relu(out)
        return out
class AdaptiveRegressionCNN(nn.Module):
    def __init__(self, input_size):
        super(AdaptiveRegressionCNN, self).__init__()

        # out_channels_conv1 = min(64, input_size // 10)
        # out_channels_conv2 = min(32, input_size // 20)
        out_channels_conv1=64
        out_channels_conv2=32

        self.conv1 = nn.Conv1d(in_channels=1, out_channels=out_channels_conv1, kernel_size=3, padding=1)
        self.resblock1 = ResidualBlock(out_channels_conv1)  # 第一个残差块

        self.conv2 = nn.Conv1d(in_channels=out_channels_conv1, out_channels=out_channels_conv2, kernel_size=3, padding=1)
        self.resblock2 = ResidualBlock(out_channels_conv2)  # 第二个残差块

        self.leaky_relu = nn.LeakyReLU(0.01)

        self._to_linear = None
        self._calculate_to_linear(input_size)

        self.fc1 = nn.Linear(self._to_linear, 512)
        self.fc2 = nn.Linear(512, 1)

    def _calculate_to_linear(self, L):
        x = torch.randn(1, 1, L)
        x = self.leaky_relu(self.conv1(x))
        x = self.resblock1(x)
        x = self.leaky_relu(self.conv2(x))
        x = self.resblock2(x)
        self._to_linear = x.numel() // x.size(0)

    def forward(self, x):
        x = self.leaky_relu(self.conv1(x))
        x = self.resblock1(x)
        x = self.leaky_relu(self.conv2(x))
        x = self.resblock2(x)
        x = x.view(x.size(0), -1)
        x = self.leaky_relu(self.fc1(x))
        x = self.fc2(x)
        return x

def cross_validation(Dmodel, x, y, patience, epochs, gene_name,device,n_splits=5):
    kf = KFold(n_splits=n_splits, shuffle=True, random_state=42)
    r2_scores = []
    best_weights_per_fold = []
    r2_scores_per_fold = []     
   
    for train_idx, val_idx in kf.split(x):
        x_train, x_val = x[train_idx], x[val_idx]
        y_train, y_val = y[train_idx], y[val_idx]
        x_train = torch.tensor(x_train, dtype=torch.float32)
        y_train = torch.tensor(y_train, dtype=torch.float32)
        x_val = torch.tensor(x_val, dtype=torch.float32)
        y_val = torch.tensor(y_val, dtype=torch.float32)
        device = device
        model = Dmodel
        model.to(device)
        optimizer = torch.optim.Adam(model.parameters(), lr=0.001)
        criterion = nn.MSELoss()

        best_val_loss = float('inf')
        count = 0

        for epoch in range(epochs):
            model.train()
            optimizer.zero_grad()
            y_pred = model(x_train)
            loss = criterion(y_pred, y_train)
            loss.backward()
            optimizer.step()

            model.eval()
            with torch.no_grad():
                y_pred_val = model(x_val)
                val_loss = criterion(y_pred_val, y_val)
                if val_loss.item() < best_val_loss:
                    best_val_loss = val_loss.item()
                    count = 0
                    best_model_state = model.state_dict() 
                else:
                    count += 1
                    if count >= patience:
                        model.load_state_dict(best_model_state)  # Load the best model weights
                        break

        model.eval()
        with torch.no_grad():
            y_pred_val = model(x_val).to(device)
            last_pre = [i[0] for i in y_pred_val.cpu().numpy().tolist()]
            last_val = [i[0] for i in y_val.cpu().numpy().tolist()]
            r = pearsonr(last_pre, last_val)
            r2 = math.pow(r[0], 2)
            r2_scores.append(r2)
            r2_scores_per_fold.append(r2)
            # best_weights_per_fold.append(deepcopy(model.state_dict()))
        del model, optimizer, criterion, x_train, x_val, y_train, y_val, y_pred_val
        torch.cuda.empty_cache()
        gc.collect()
    # 找到最佳模型的索引
    # best_fold_index = np.argmax(r2_scores_per_fold)
    # 加载最佳模型权重
    # best_model = Dmodel
    # best_model.load_state_dict(best_weights_per_fold[best_fold_index])
    # torch.save(best_model.state_dict(), './CPG_model/{}_best_model.pth'.format(gene_name))
    return np.mean(r2_scores)

not_singal_gene=[]
reslut_dic={}
# for i in tqdm(list(gene_list["x"])):
for i in tqdm(list(gene_list["x"])):
    gene_name=promoter_range.loc[promoter_range["gene name"]==i]
    if len(gene_name)==1:
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
        if len(a.columns)==1:
            y=need_BRCA_exp_t[i]
            df = pd.concat([need_BRCA_meth_T,y],axis=1)
            data=df.sample(frac=1,random_state=1)
            train_x=data.drop([i],axis="columns")
            train_y=data[i]
            x = torch.tensor(train_x.values).float().unsqueeze(1) 
            y=torch.unsqueeze(torch.tensor(train_y.values),dim=1).float()
            x=x.to(device)
            y=y.to(device)
            input_size = x.shape[2]
            # output_size = y.shape[1]
            torch.manual_seed(42)
            if torch.cuda.is_available():
                torch.cuda.manual_seed_all(42)
    
            Dmodel=AdaptiveRegressionCNN(input_size=input_size).to(device)
            mean_r2 = cross_validation(Dmodel,x, y, patience=90,epochs=700,gene_name=i,device=device)
            File = open("../result_r2/CNNnormal.txt",'a')  
            File.write(i+"\t"+str(mean_r2)+"\n")
            File.flush()
        
        else:
            mean_=[]
            for len_ in range(len(a.columns)):
                y=a.iloc[:,len_]
                df = pd.concat([need_BRCA_meth_T,y],axis=1)
                data=df.sample(frac=1,random_state=1)
                train_x = data.drop([i], axis="columns").values
                train_y = data[i].values
                mean_r2 = cross_validation(Dmodel,x, y, patience=90,epochs=700,gene_name=i,device=device)
                mean_.append(mean_r2)
            File = open("../result_r2/CNNnormal.txt",'a')  
            File.write(i+"\t"+str(max(mean_))+"\n")
            File.flush()
            print(i,mean_)
            print(i+": ",max(mean_))
