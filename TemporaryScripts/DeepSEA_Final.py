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

#We will define some additional imports to use GPU
try:
  if(torch.cuda.is_available()):
      print("GPU successfully detected - ")
      print(torch.cuda.get_device_name(0))
      device = torch.device("cuda:0")
except Exception as e:
  print("GPU not detected. Change the settings as mentioned earlier and run session again")
  device = torch.device("cpu")

#Setting random seed for reproducibility
SEED = 1234
torch.manual_seed(SEED)
np.random.seed(SEED)
torch.backends.cudnn.deterministic = True


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




# Create a DeepSEA network and send it to GPU
CNN = DeepSEA().to(device)

# Read the best parameters shared by the authors of the paper
best_model = torch.load("/N/slate/ppugale/DEEPSEA_FULLDATA/Complete_data/deepsea_bestmodel.pkl")

# load the best parameters to newly created DeepSEA network
CNN.load_state_dict(best_model)
print("best model loaded")

# ------------------------------------------------------------------------------
# Add code here
# Define cost function and send it to GPU.
cost_function = nn.BCEWithLogitsLoss().to(device)

seq_list = pickle.load(open("Sequences.pkl","rb"))
final_results = {}
detailed_final_results = {}
for rsid in seq_list.keys():
    print(f"Running ... {rsid}")
    seq = seq_list[rsid][0]
    a1 = seq_list[rsid][1]
    a2 = seq_list[rsid][2]
    if seq != None:
        log_changes,P_ref,P_alt = Run_Deepsea(seq,a1,a2)
        detailed_final_results[rsid] = [P_ref,P_alt]
        final_results[rsid] = log_changes

print(final_results[rsid])  #Printing the last rsid
output = open('Final_Output.pkl','wb')
pickle.dump(final_results,output)
output.close()

output2 = open('Detailed_Final_Output.pkl','wb')
pickle.dump(detailed_final_results,output2)
output2.close()