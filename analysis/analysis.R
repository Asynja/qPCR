library(tidyverse)

#get raw data (with replicate number)
bigqPCRraw<-read.table("qPCR_3replicates.txt", header = T)

#join replicate and target, to be able to spread later
bigqPCRraw %>% 
  transmute(Strain, Target_replica = paste(bigqPCRraw$Target, bigqPCRraw$Replicate, sep = "_"), CT) %>% 
  spread(Target_replica, CT)-> big_spread

#normalize (PP1-ref and PP2-ref):
big_spread %>% 
  transmute(Strain, NormalisedPP1_1 = PP1_1 - ref_1, NormalisedPP1_2 = PP1_2 - ref_2, 
            NormalisedPP1_3 = PP1_3 - ref_3,NormalisedPP1_4 = PP1_4 - ref_4, 
            NormalisedPP1_5 = PP1_5 - ref_5,NormalisedPP1_6 = PP1_6 - ref_6,
            NormalisedPP1_7 = PP1_7 - ref_7, NormalisedPP1_8 = PP1_8 - ref_8, 
            NormalisedPP1_9 = PP1_9 - ref_9,
            NormalisedPP2_1 = PP2_1 - ref_1, NormalisedPP2_2 = PP2_2 - ref_2, 
            NormalisedPP2_3 = PP2_3 - ref_3,NormalisedPP2_4 = PP2_4 - ref_4, 
            NormalisedPP2_5 = PP2_5 - ref_5,NormalisedPP2_6 = PP2_6 - ref_6,
            NormalisedPP2_7 = PP2_7 - ref_7, NormalisedPP2_8 = PP2_8 - ref_8, 
            NormalisedPP2_9 = PP2_9 - ref_9
            ) -> normalised_table

#"flip" the values in order to do the strain subtraction:
normalised_table %>% 
  group_by(Strain) %>% 
  gather(target, value, -Strain) %>% 
  spread(Strain, value) -> pre_subtraction

#Get the average of each target in WT, 1081:
average_refWT = apply(big_spread[1,2:10], 1, mean, na.rm = TRUE)
average_PP1WT = apply(big_spread[1,11:19], 1, mean, na.rm = TRUE)
average_PP2WT = apply(big_spread[1,20:28], 1, mean, na.rm = TRUE)

average_norm_PP1WT = average_PP1WT - average_refWT
average_norm_PP2WT = average_PP2WT - average_refWT

#Now the appropriate subtraction between strains (mutant-WT, with every combination of plasmid)
#In this case, I will separate PP1 and PP2 in two diferent tables, to be able to subtract the average of the WT in each case.
pre_subtraction %>% 
  filter( target == "NormalisedPP1_1"| target == "NormalisedPP1_2"| 
            target == "NormalisedPP1_3" | target == "NormalisedPP1_4" |
            target == "NormalisedPP1_5" | target == "NormalisedPP1_6" |
            target == "NormalisedPP1_7" | target == "NormalisedPP1_8" |
            target == "NormalisedPP1_9" ) -> PP1_table

pre_subtraction %>% 
  filter(  target == "NormalisedPP2_1"| target == "NormalisedPP2_2"| 
             target == "NormalisedPP2_3" | target == "NormalisedPP2_4" |
             target == "NormalisedPP2_5" | target == "NormalisedPP2_6" |
             target == "NormalisedPP2_7" | target == "NormalisedPP2_8" |
             target == "NormalisedPP2_9") -> PP2_table


#Subtract in each case
PP1_table %>% 
  mutate(subtracted_WT = WT - average_norm_PP1WT, 
         subtracted_mutant = mutant - average_norm_PP1WT,
         subtracted_WT_plasmid_1 = WT_plasmid_1 - average_norm_PP1WT,
         subtracted_mutant_plasmid_1 = mutant_plasmid_1 - average_norm_PP1WT,
         subtracted_WT_plasmid_2 = WT_plasmid_2 - average_norm_PP1WT,
         subtracted_mutant_plasmid_2 = mutant_plasmid_2 - average_norm_PP1WT,
         subtracted_mutant_plasmid_3 = mutant_plasmid_3 - average_norm_PP1WT,
         subtracted_WT_plasmid_3 = WT_plasmid_3 - average_norm_PP1WT,
         subtracted_WT_plasmid_4 = WT_plasmid_4 - average_norm_PP1WT,
         subtracted_mutant_plasmid_4 = mutant_plasmid_4 - average_norm_PP1WT
  ) -> Subtracted_PP1

PP2_table %>% 
  mutate(subtracted_WT = WT - average_norm_PP2WT, 
         subtracted_mutant = mutant - average_norm_PP2WT,
         subtracted_WT_plasmid_1 = WT_plasmid_1 - average_norm_PP2WT,
         subtracted_mutant_plasmid_1 = mutant_plasmid_1 - average_norm_PP2WT,
         subtracted_WT_plasmid_2 = WT_plasmid_2 - average_norm_PP2WT,
         subtracted_mutant_plasmid_2 = mutant_plasmid_2 - average_norm_PP2WT,
         subtracted_mutant_plasmid_3 = mutant_plasmid_3 - average_norm_PP2WT,
         subtracted_WT_plasmid_3 = WT_plasmid_3 - average_norm_PP2WT,
         subtracted_WT_plasmid_4 = WT_plasmid_4 - average_norm_PP2WT,
         subtracted_mutant_plasmid_4 = mutant_plasmid_4 - average_norm_PP2WT
  ) -> Subtracted_PP2

#2**
Subtracted_PP1 %>% 
  transmute(target, NoPlasmid_WT = 2**(-subtracted_WT), NoPlasmid_mutant = 2**(-subtracted_mutant), 
            plasmid_1_WT = 2** (-subtracted_WT_plasmid_1),plasmid_1_mutant = 2**(-subtracted_mutant_plasmid_1),
            plasmid_2_WT = 2**(-subtracted_WT_plasmid_2), plasmid_2_mutant = 2**(-subtracted_mutant_plasmid_2), 
            plasmid_3_WT = 2**(-subtracted_WT_plasmid_3), plasmid_3_mutant = 2**(-subtracted_mutant_plasmid_3),
            plasmid_4_WT = 2**(-subtracted_WT_plasmid_4), plasmid_4_mutant = 2**(-subtracted_mutant_plasmid_4)) -> TwoDelta_strains_PP1

Subtracted_PP2 %>% 
  transmute(target, NoPlasmid_WT = 2**(-subtracted_WT), NoPlasmid_mutant = 2**(-subtracted_mutant), 
            plasmid_1_WT = 2** (-subtracted_WT_plasmid_1),plasmid_1_mutant = 2**(-subtracted_mutant_plasmid_1),
            plasmid_2_WT = 2**(-subtracted_WT_plasmid_2), plasmid_2_mutant = 2**(-subtracted_mutant_plasmid_2), 
            plasmid_3_WT = 2**(-subtracted_WT_plasmid_3), plasmid_3_mutant = 2**(-subtracted_mutant_plasmid_3),
            plasmid_4_WT = 2**(-subtracted_WT_plasmid_4), plasmid_4_mutant = 2**(-subtracted_mutant_plasmid_4)) -> TwoDelta_strains_PP2

#gather and spread in order to do the average in both PP1 and PP2
TwoDelta_strains_PP1 %>% 
  group_by(target) %>% 
  gather(Strain, value, -target) %>% 
  spread(target, value) -> table_pre_average_PP1

TwoDelta_strains_PP2 %>% 
  group_by(target) %>% 
  gather(Strain, value, -target) %>% 
  spread(target, value) -> table_pre_average_PP2

#Average in PP1 and PP2
table_pre_average_PP1 %>% 
  mutate(PP1 = apply(table_pre_average_PP1[,2:10], 1, mean, na.rm = TRUE),
         PP1_sd = apply(table_pre_average_PP1[,2:10], 1, sd, na.rm = TRUE))-> average_PP1

table_pre_average_PP2 %>% 
  mutate(PP2 = apply(table_pre_average_PP2[,2:10], 1, mean, na.rm = TRUE),
         PP2_sd = apply(table_pre_average_PP2[,2:10], 1, sd, na.rm = TRUE))-> average_PP2

#PLOT PP1
average_PP1 %>% 
  ggplot(aes(x = Strain, y =PP1, ylim (3))) + geom_col(width = 0.5, fill = "48") + 
  geom_errorbar(aes(ymin = PP1 - PP1_sd, ymax = PP1 + PP1_sd), width = 0.25) + 
  xlab ("Strain") + ylab("PP1") + theme_bw() + 
  coord_cartesian(xlim = NULL, ylim = c(0,100), expand = TRUE) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#PLOT PP2
average_PP2 %>% 
  ggplot(aes(x = Strain, y =PP2)) + geom_col(width = 0.5, fill = "48") + 
  geom_errorbar(aes(ymin = PP2 - PP2_sd, ymax = PP2 + PP2_sd), width = 0.25) + 
  xlab ("Strain") + ylab("PP2") + theme_bw() + 
  coord_cartesian(xlim = NULL, ylim = c(0,0.2), expand = TRUE) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
