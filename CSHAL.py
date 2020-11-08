
from Bio import Entrez, SeqIO
import numpy as np
import pandas as pd
import os
import pickle
import torch
import torch.optim as optim
import torch.nn as nn
import torch.utils.data as Data
import torch.nn.functional as F
import time
from sklearn.metrics import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pickle
import sys
from pyfiglet import Figlet
f = Figlet(font='slant')
print f.renderText('C-SHAL')


def one_hot_encode(seq):
    seq = seq.lower()
    mat_list = [np.eye(4)[i] for i in range(4)]
    mapping = dict(zip(['a','c','g','t'],mat_list))
    seq_one_hot =  np.stack([mapping[i] for i in seq]).T   
    return(mapping,seq_one_hot)


# Define the DeepSEA CNN model
class DeepSEA(nn.Module):
    def __init__(self, ):
        super(DeepSEA, self).__init__()
        self.Conv1 = nn.Conv1d(in_channels=4, out_channels=320, kernel_size=8)
        self.Conv2 = nn.Conv1d(in_channels=320, out_channels=480, kernel_size=8)
        self.Conv3 = nn.Conv1d(in_channels=480, out_channels=960, kernel_size=8)
        self.Maxpool = nn.MaxPool1d(kernel_size=4, stride=4)
        self.Drop1 = nn.Dropout(p=0.2)
        self.Drop2 = nn.Dropout(p=0.5)
        self.Linear1 = nn.Linear(53*960, 925)
        self.Linear2 = nn.Linear(925, 919)
    def forward(self, input):
        x = self.Conv1(input)
        x = F.relu(x)
        x = self.Maxpool(x)
        x = self.Drop1(x)
        x = self.Conv2(x)
        x = F.relu(x)
        x = self.Maxpool(x)
        x = self.Drop1(x)
        x = self.Conv3(x)
        x = F.relu(x)
        x = self.Drop2(x)
        x = x.view(-1, 53*960)
        x = self.Linear1(x)
        x = F.relu(x)
        x = self.Linear2(x)
        return x

def log_change(P_ref,P_alt):
  term1 = np.log(P_ref/(1-P_ref))
  term2 = np.log(P_alt/(1-P_alt))
  return abs(term1-term2)

def Run_Deepsea(seq,a1,a2):
    seq_ref = seq[0:500] + a1 + seq[501:]
    mapping,seq_one_hot_ref = one_hot_encode(seq_ref)
    x_ref = torch.tensor(seq_one_hot_ref.reshape(1,4,1000),dtype=torch.float).to(device)  #reshaping, converting to tensor and putting on GPU. We nee shape (1,4,1000) since 1st shape represents no of data points
    seq_alt = seq[0:500] + a2 +seq[501:]
    mapping, seq_one_hot_alt = one_hot_encode(seq_alt)
    x_alt =  torch.tensor(seq_one_hot_alt.reshape(1,4,1000),dtype=torch.float).to(device)
    CNN.eval()
    with torch.no_grad(): 
        y_pred_ref =  CNN(x_ref)
        y_pred_alt = CNN(x_alt)
    P_ref =  torch.sigmoid(y_pred_ref.data)   
    P_alt = torch.sigmoid(y_pred_alt.data)
    P_ref = P_ref.detach().cpu().reshape(-1).numpy()  
    P_alt = P_alt.detach().cpu().reshape(-1).numpy()  
    LC = log_change(P_ref,P_alt)
    sorted_LC = np.sort(LC)[::-1]
    return sorted_LC,P_ref,P_alt

def Get_seq(start_pos,end_pos,chr):
    
    try:
        if int(chr) < 10:
            id_chr = "".join(["NC_00000",chr])
        else:
            id_chr = "".join(["NC_0000",chr])

        handle = Entrez.efetch(db="nucleotide",
                        id = id_chr,
                        rettype = "fasta",
                        strand = 1,
                        seq_start = start_pos,
                        seq_stop  = end_pos)
        record = SeqIO.read(handle,"fasta")
        
        return str(record.seq)
    except:
        print("No proper chromosome found ...")




if 
