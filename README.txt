Andor Budai (2020); 
Eötvös University, Institute of Physics, 1117 Budapest, Hungary; 
email: arandras@caesar.elte.hu

Data: folder containing the light curves
VarErr: folder containing the randomised light curve variabilities

If you want to analyse the sample -> run analyseError.m

WHAT WE DID:
We have ran anlyseError(10000, 10) to create resultTableErr.csv, table containing the angles and variabilities.
We have ran corrError(10000, 10) to create RP16Sper.csv and RP19Sper.csv, tables containing the correlation test results for the sample containing 16 GRBs and the 8 samples containing 19 GRBs respectively.

FUNCTION TREE:
	1. analyseError
		varMeasureMod
	
	2. corrError.m