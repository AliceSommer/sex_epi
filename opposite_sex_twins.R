
phenotype[phenotype$ZYGOSITY == "DZ_OppositeSex",]

DZ_ID <- phenotype$SAMPLE[phenotype$PID %in% c(159,11)]
DZ_ID_11 <- phenotype$SAMPLE[phenotype$PID == 11]
DZ_ID_159 <- phenotype$SAMPLE[phenotype$PID == 159]

# subset data to samples of DZ opposite sex twins
opposite_sex_melt <- cgs_interest_melt[cgs_interest_melt$variable %in% DZ_ID,]
opposite_sex_melt$pair_nb <- "Pair 1"
opposite_sex_melt$pair_nb[opposite_sex_melt$variable %in% DZ_ID_11] <- "Pair 2" # DZ_ID_11 == 2

g_twin <- ggplot(opposite_sex_melt) +
  geom_point(aes(x = value, y = cg_name, colour = sex)) +
  facet_grid(scenario ~ pair_nb, scales = "free", space = "free") + theme_minimal() +
  ggtitle("Opposite sex twins") + xlab('x 100 (%)')
  
ggsave(file = 'twin_opposite_sex.jpeg', g_twin,
         dpi=300,
         width = 150,
         height = 160,
         units = "mm")

