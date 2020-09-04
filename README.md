# qPCR

Analysis of qPCR using the 2∆∆Ct method, based on the article "Analysis of Relative GeneExpression Data Using Real-Time Quantitative PCR and the 2∆∆CtMethod" by Livak, Kenneth J. and Thomas D. Schmittgen (2001).

In order to do the the analysis, we assume that the amplification efficiencies of both the target and reference must be approximately equal. A way of testing that assumption is to see how ∆Ct varies with template dilution. If the assumption is correct, the absolute value of the slope will be 0.  This analysis can be seen in the folder "concentration". It uses a txt file with the data organized in three columns: concentration, Target and CT.
