# qPCR

Analysis of qPCR using the 2∆∆Ct method, based on the article "Analysis of Relative Gene Expression Data Using Real-Time Quantitative PCR and the 2∆∆CtMethod" by Livak, Kenneth J. and Thomas D. Schmittgen (2001).

In order to do the the analysis, we assume that the amplification efficiencies of both the target and reference must be approximately equal. A way of testing that assumption is to see how ∆Ct varies with template dilution. If the assumption is correct, the absolute value of the slope will be 0.  This analysis can be seen in the folder "concentration". It uses a txt file with the data organized in three columns: concentration, Target and CT.

For the qPCR analysis, two primer pairs were used, and compared it a third, used as a reference. PP1 refers to the first primer pair, PP2 to the second, and ref to the third. The strains used were the wild type (WT) and mutant, in combination with four different plasmids. The txt file with the data is organized in 4 columns, with the replicate number (in my case, 3 replicates), the strain name, the target or primer pair, and the CT values obtained.  
