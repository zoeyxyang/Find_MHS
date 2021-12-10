# columns are cell types, rows are drugs
import pandas as pd
import time
from functions import MatchIdName


def optimal_subset(remaining_set):
    number_of_drugs_targeted = remaining_set.sum(axis=1)
    number_of_drugs_targeted = pd.to_numeric(number_of_drugs_targeted)
    drug = number_of_drugs_targeted.idxmax()
    drug_targets = remaining_set.loc[drug, remaining_set.loc[drug] == 1].keys()
    return drug, drug_targets


def greedy_set_cover(df):
    remaining_df = df
    min_subsets = set()
    while len(remaining_df.columns) != 0:
        drug, drug_targets = optimal_subset(remaining_df)
        remaining_df = remaining_df.drop(index=drug, columns=drug_targets)
        min_subsets.add(drug)
    return min_subsets


def download_as_pd(file):
    df = pd.read_csv(file, index_col=0)
    return df


def GD(bin_df, dataDir):
    
    start = time.time()
    HSet_NSC = list(greedy_set_cover(bin_df))
    end = time.time()
    HSet_Name = MatchIdName(HSet_NSC, dataDir)
    
    print("runtime for the algorithm is:", end - start) 
    print("Minimal Hitting Set {}: {}".format(1, HSet_NSC))
    print("Minimal Hitting Set {}: {}".format(1, HSet_Name))
    
    
    return HSet_Name




if __name__ == "__main__":
    
    
    start = time.time()
    df = download_as_pd("all_zero_one_matrix.csv")
    min_set = greedy_set_cover()

    print("minimum set solution:", min_set)
    end = time.time()
    print("Program Runtime:", end-start, "seconds")
    # time returns seconds since it was run on unix

