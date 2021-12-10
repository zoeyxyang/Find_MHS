import numpy as np
import pandas as pd
import random
import multiprocessing
import time
from functions import MatchIdName


def provide_params(i, df):

    global df_current, strain_to_drug


    subset = random.sample(list(df.index), i)
    df_current = df.loc[subset]  ##Get df corresponding to just the randomly selected drugs. We will do SA on this df
    strain_to_drug = create_graph(df_current)
    min_HS = simulated_annealing(0, i * 10, 0.1, 20)  ##Do SA
    return min_HS

def create_graph(df):
    """Create mapping from strain to drugs
    Input: Hitting matrix with 0s and 1s
    Returns: Dictionary mapping cell-lines to the set of drugs they respond to"""
    strain_to_drug = {}
    for i in range(0, len(df.columns)):
        strain = df.columns[i]
        strain_to_drug[strain] = set([drug for drug in df.index if df.loc[drug, strain] == 1])

    return strain_to_drug

def optimization_fn(drug_set):
    """Function to optimize"""
    drugs = [v for k, v in strain_to_drug.items() if k in drug_set]
    return len(set().union(*drugs)) ##Number of covered drugs - we want to minimize this

def check_drug_removal(drug, HS):
    """Check if drug can be removed from current HS
    Input: Drug being tested, Hitting Set
    Returns: Boolean - True means the drug needs to be retained in HS
                        False means drug can be removed from HS"""
    for k, v in strain_to_drug.items():
        if drug in v: #if drug hits that strain
            common_drugs = set(v).intersection(HS)
            if common_drugs == {drug} and len(common_drugs) == 1:
                return True
    return False

def simulated_annealing(beta_0, t_eq, delta_beta, beta_max):
    """Performs simulated annealing
    Input: Parameters for SA:   beta_0 (initial temperature)
                                t_eq (number of equilibration steps
                                delta_beta (step size for increasing temperature)
                                beta_max (maximum temperature)
    Output: List of Minimal Hitting Sets (MHS) and their sizes.
            MHSs comprise the drug indices in the input matrix"""
    temp = beta_0
    random_HS = []
    all_drugs = df_current.index
    ###Generating a random Hitting Set (HS) for initialization
    for k, v in strain_to_drug.items():
        if len(v) != 0:
            random_HS.append(random.choice(list(v)))

    random_HS = list(set(random_HS))

    ####Perform T_eq equilibration steps
    while temp < beta_max:
        for i in range(0, t_eq):
            selected_drug = random.choice(all_drugs)
            if selected_drug in random_HS: ##if drug is covered
                if check_drug_removal(selected_drug, random_HS) == False: ##If no set becomes uncovered
                    random_HS.remove(selected_drug) ##Drug is redundant; update HS by removing drug
            else: ##if drug is uncovered
                metropolis_criterion = np.exp(-temp)
                if random.uniform(0, 1) < metropolis_criterion:
                    random_HS.append(selected_drug) ##Accept new drug into HS with probability proportional to MC
        temp += delta_beta
    return random_HS


def SA(bin_df, dataDir):
    a = time.time()
    num_drugs = random.choices(range(20, 50), k=300)  ##Pick a random number of drugs to test for MHS
    df = bin_df
    pool = multiprocessing.Pool(8)
    args = [[i, df] for i in num_drugs]
    HSets = pool.starmap(provide_params, args)
    pool.close()

    print("runtime for the algorithm is:", time.time() - a)
    
    lengths = [len(x) for x in HSets]
    print("Length of Minimal Hitting Set(s): {}".format(min(lengths)))
    print("List of Minimal Hitting Sets: ")
    MHS_NSCs = [x for x in HSets if len(x) == min(lengths)] 
    HS_sizes = []
    for i in range(1, len(MHS_NSCs)+1):
        MHS_NSC = MHS_NSCs[i-1]
        HS_sizes.append(len(MHS_NSC))
        print("Minimal Hitting Set {}: {}".format(i, MHS_NSC))
        print("Minimal Hitting Set {}: {}".format(i, MatchIdName(MHS_NSC, dataDir)))
        
    return HS_sizes
    
if __name__ == "__main__":
  
    a = time.time()
    num_drugs = random.choices(range(20, 50), k=300)  ##Pick a random number of drugs to test for MHS
    df = pd.read_csv("all_zero_one_matrix.txt", sep="\t", index_col=0)
    pool = multiprocessing.Pool(8)
    args = [[i, df] for i in num_drugs]
    HSets = pool.starmap(provide_params, args)
    pool.close()

    print(time.time() - a)

    lengths = [len(x) for x in HSets]
    print("Length of Minimal Hitting Set(s): {}".format(min(lengths)))
    print("List of Minimal Hitting Sets: ")
    MHS = [x for x in HSets if len(x) == min(lengths)]
    for i in range(1, len(MHS)+1):
        print("Minimal Hitting Set {}: {}".format(i, MHS[i-1]))
