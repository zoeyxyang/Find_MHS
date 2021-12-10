# Find_MHS
Minimal Drug Combinations Targeting the National Cancer Institute 60-Cell Line Database





## CMU Course Project - 02-712 Computational Methods for Biological Modeling and Simulation
Ben Oppenheimer, Kavya Prasad, Michael Aksu, Tomas Matteson, Zoey Yang

## Short Introduction of the Software
We wrote software in Python to find the minimal drug combinations for the NCI-60 cancer cell line database using multiple different algorithms.

Users can set the parameters for (1) tissue of origin of cancer cell lines they want to analyze (2) a threshold for the normalized drug activities (the program assumes that a drug is active against a specific cell line if the activity is greater than the threshold. Otherwise, it is considered inactive) (3) the algorithm to run for the minimal hitting set problem. 


## Usage


Parameters:


-t: tissue of origin 
(users would be able to choose any combinations of the following organs: breast, cns, colon, leukemia, melanoma, non-small cell lung, ovarian, prostate, renal, and all)
                 
-th: threshold (a float value)  

-a: algorithm 
(users would be able to choose any of the following algorithms to run: SA, ILP, HDF, GD)
(SA: Simulated Annealing,  ILP: Integer Linear Programming, HDF: Highest Degree First, GD:Greedy Algorithm)





Example for finding the whole dataset with a threshold of 0 and Highest Degree First Algortihm


```bash

cd src

python3 main.py -t all -th 0 -a HDF

```




## Introduction
Cancer is a deadly disease, and one scientists are still fighting with. There are a number of factors that make it hard to cure cancer - rapidly mutating cell-lines, fitness advantage of cancer cells, metastasis etc. It can also be highly specific, varying in molecular phenotype from patient to patient. Although there is greater success in curing cancer today as opposed to a few decades ago, we still face some fundamental challenges in treatment. 
One of the fundamental challenges in fighting cancer is drug resistance in cancer. This challenge occurs when some cancerous cells are resistant to a cancer drug. Therefore, the drug is ineffective in preventing the cancer from spreading, because these cells will fill the place of other cancer cells. Consequently, it is important to create a combination of drugs that targets all cancer cell types. It is likely that creating the minimum combination of drugs would be more affordable to the patient and would reduce side-effects. 
A more theoretical abstraction of this problem is known as the minimum set cover problem.  Using this abstraction, Alexei Vasquez developed an approach to finding the minimum combination of drugs [1]. Specifically, given a collection of cancerous cell types and drugs, can we find the smallest set of drugs that targets all of the cell types? 
Besides its clinical relevance, this problem also is computationally challenging. The minimum set cover problem is NP-complete, which means that no one has been able to devise an algorithm that allows a computer to find the optimal solution in polynomial time when the dataset gets large. However, computer scientists have developed algorithms that produce solutions that are at least close to the optimal solution and can be executed in polynomial time. 
To tackle the minimum drug combination problem, we implemented and compared a few heuristic algorithms. We used the National Cancer Institute’s NCI-60 Human Tumor Cell Lines screen. This is where we enter the domain of optimization and simulation algorithms!


## References
[1] Vazquez, Alexei. "Optimal drug combinations and minimal hitting sets." BMC systems biology 3.1 (2009): 1-6.

[2] “NCI-60 Human Tumor Cell Lines Screen.” Developmental Therapeutics Program (DTP), https://dtp.cancer.gov/discovery_development/nci-60/.  

[3] Reinhold, William C et al. “CellMiner: a web-based suite of genomic and pharmacologic tools to explore transcript and drug patterns in the NCI-60 cell line set.” Cancer research vol. 72,14 (2012): 3499-511. doi:10.1158/0008-5472.CAN-12-1370











