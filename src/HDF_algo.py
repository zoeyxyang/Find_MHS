import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from collections import Counter
import time
from read_data import get_matrix
import os
from functions import MatchIdName





'''

Input: the dataframe (matrix) to be run for, the data directory, number of time to be run

Output: a list of the length of the hitting sets

'''
def HDF(df, dataDir, n):
    a = time.time()
    HS_sizes = []
    MHS_NSCs = []
    for i in range(n):
        HSet_NSC = HDF_once(df,dataDir)
        MHS_NSCs.append(HSet_NSC)
        HS_sizes.append(len(HSet_NSC))
        
    print("runtime for the algorithm is:", time.time()-a)   
        
    print("List of Minimal Hitting Sets: ")    
    
    for i in range(n):
        HSet_NSC = MHS_NSCs[i]
        HSet_Name = MatchIdName(HSet_NSC, dataDir)
        print("Minimal Hitting Set {}: {}".format(i+1, HSet_NSC))
        print("Minimal Hitting Set {}: {}".format(i+1, HSet_Name))
        
    return HS_sizes


'''

Input: the dataframe (matrix) to be run for, the data directory

Output: the hitting set in NSC 

'''

def HDF_once(df, dataDir):
    HSets = [] #drugs in NSC index
    currDF = CalculateDegree(df)
    
    while currDF.shape[1] > 1: #stop when no strain is left
        
        #identify largest connection:
        largestNSCs = currDF[currDF["Degree"] == currDF["Degree"].max()].index.tolist()
        #randomly pick one drug with largest connection:
        chosenDrugIdx = np.random.randint(len(largestNSCs),size = 1)
        chosenDrugIdx = int(chosenDrugIdx[0])
        chosenDrugNSC = largestNSCs[chosenDrugIdx]
        #cover that drug
        HSets.append(chosenDrugNSC)
        #update the dataframe
        update = UpdateDF(currDF, chosenDrugNSC)
        currDF = update
 
    return HSets



def CalculateDegree(df):
    arr = df.to_numpy()
    df["Degree"] = np.sum(arr, axis=1)
    return df

def UpdateDF(df,drugNSC):
    matchStrains = [col for col in df.columns if df.loc[drugNSC,[col]].eq(1).any()]
    #print("matchStrains:", matchStrains)
    if 'Degree' in matchStrains :
        matchStrains.remove('Degree')
    updateDF = df.drop(matchStrains,axis=1)
    updateDF = updateDF.drop(drugNSC, axis=0)
    return updateDF
