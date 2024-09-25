#######################################
#-------------------------------------#
#   correcting for tag leakage        #                          
#-------------------------------------#
#######################################

# first it takes the 0.1 % and then it substracts this 
depr_df <- data_set %>% group_by(Species) %>%
  mutate(countSpec = ((0.001*sum(c_across(X120A:ZVLC))))) 

#roundging of to 0 digits
depr_df$countSpec <- round(depr_df$countSpec, digits = 0)

#substract from the columns of the samples that contain abundanca data
depr_minus <- depr_df %>%
  mutate(across(11:23, ~ . - countSpec))

#Make sure minus value become 0
depr_minus <- as.data.frame(depr_minus)
depr_minus[depr_minus < 0] <- 0