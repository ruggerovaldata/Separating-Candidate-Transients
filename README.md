# Separating Candidate Transients

Step 1: 

Inserting the bananas credentials in dblogin.py

Running the ExtractCsv.py file to obtain the .csv files contaninig the eta and V parameters for the datasets dersired. It is possible to select the Filter radius, out of which excluding the sources. 

Step 2: 

Running the FindingCandidateTransient.py. Python version: 3.9. Dependecies: matplotlib, scipy, pandas, tabulate. 

Inserting the .csv file names with the parameters that want to be analysed, without including the .csv extension. 

Indicating the p value that want to be used, i.e. the percentage of certainty with which a sources will be classified as an inlier. For example: if p=0.99 the sources that will be classified as inlier are 99% certain of being inlier. 
