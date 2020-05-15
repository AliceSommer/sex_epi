
head(cgs)
head(gap_value)
head(cgs_interest_melt)

# set the number of randomizations
nrep <- 10^4

# create a matrix where the t_rand will be saved
t_arrays <- matrix(NA, ncol=length(cgs), nrow=nrep)

for(i in 1:length(cgs)){
  print(i)
  dat_cg = cgs_interest_melt[cgs_interest_melt$cg_name == cgs[i],]
  for(j in 1:nrep){
    W_rep = sample(dat_cg$sex)
    # retrieve the values of the imputed POs for the units assigned to treatment
    newy_F = dat_cg$value[W_rep == "F"]
    # retrieve the values of the imputed POs for the units assigned to control
    newy_M = dat_cg$value[W_rep == "M"]
    # calculate the test statistic
    new_gap_stat = stat(newy_F,newy_M)
    # fill t_arrays 
    t_arrays[j,i] = new_gap_stat 
  }
}

## calculate p_value
p_values <- apply(t_arrays, 2, function(x) mean(x >= 75))

## melt t_arrays
t_arrays_data_frame <- data.frame(t_arrays)
colnames(t_arrays_data_frame) <- cgs

t_array_melt <- melt(t_arrays_data_frame)

g_rand <- ggplot(t_array_melt,aes(x = value)) + 
  facet_wrap(~variable) + 
  geom_histogram(binwidth = 0.8) + 
  theme(axis.text.x = element_text(angle = -90, vjust = 0.5, hjust = -1)) 

ggsave(file = 'null_distributions.jpeg', g_rand,
       dpi=300,
       width = 150,
       height = 100,
       units = "mm")

