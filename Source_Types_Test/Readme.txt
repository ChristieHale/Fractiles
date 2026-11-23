Source types 1, 6 and 7 in one 45.4 fault file: a type 1 fault, a type 6
fault (the geometry of PEER Set 3 Test 3 with coarser steps) and a type 7
source of three ruptures (Hazard/Input/Ruptures_Type7.txt, read by HAZ only).

Each source has a single branch and the run has one GMM, so every Monte Carlo
sample is the same curve: every fractile and the mean must equal HAZ's total
hazard (out3). The fractile program reads the fault file to find each
source's branches, so this checks that it reads types 6 and 7 as HAZ does.
