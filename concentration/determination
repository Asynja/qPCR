library(tidyverse)
#Remember to set session to qPCR folde

#Plot it as is:

ct_data <- read.table('concentration_determination.txt', header = TRUE)

ct_data %>% 
  filter(concentration > 2) -> filtered_ct


filtered_ct %>% 
  group_by(Target) %>% 
  ggplot() + geom_line(aes(x= concentration, y = CT, color = Target))
  
#To check if the delta delta CT method is actually applicable in this case
  
filtered_ct %>% 
spread(Target, CT) -> filtered_spreaded_ct
  


#Logarithm of the concentration

filtered_spreaded_ct %>% 
  transmute(logarithmic = log10(concentration), fimI = (fimI_1 + fimI_2)/2, OFF = (OFF_1 + OFF_2)/2, ON = (ON_1 + ON_2)/2) -> validation_ct

validation_ct %>% 
  mutate(OFF_relative = OFF - fimI, ON_relative = ON - fimI) -> validation_ct

#Plotting everything:

validation_ct %>% 
  ggplot() + geom_point(aes(x = logarithmic, y = ON_relative)) + 
  geom_smooth(aes(x = logarithmic, y = ON_relative), method = "lm", fill = "blue") +
  geom_point(aes(x = logarithmic, y = OFF_relative)) +
  geom_smooth(aes(x = logarithmic, y = OFF_relative), method = "lm", fill = "red") +
  xlab("cDNA dilution") + ylab(expression(paste(Delta, "Ct "))) + ylim(0,6) + 
  theme_bw()
