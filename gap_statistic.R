library(reshape2)
library(ggplot2)
library(gridExtra)

# set working directory
setwd('/Users/alicesommer/Desktop/Bureau/data_sex_epi')

###################
# Bind's GAP STAT #
###################
stat <- function(a,b){
  i<-1
  
  d1<-quantile(a,probs=(100:0)/100,na.rm=TRUE)<quantile(b,probs=(0:100)/100,na.rm=TRUE) #female on the right
  d2<-quantile(a,probs=(0:100)/100,na.rm=TRUE)>quantile(b,probs=(100:0)/100,na.rm=TRUE)
  
  if(table(d1)[1]==101) return(100)
  if(table(d1)[1]<=table(d1)[2]) #female on the right
  { while(1*(d1[i])==0){i<-i+1}
    return(((100:0)/100)[i]*100)
  }
  if(table(d1)[1]>table(d1)[2]) #female on the left
  {  while(1*(d2[i])==0){i<-i+1}
    return(((100:0)/100)[i]*100)
  }
}

## twin methylation data
load('betas_preprocessed.RData')
data <- betas.preprocessed
# metadata
phenotype <- read.csv('phenotype.csv')
head(phenotype)

## bob methylation data
load("clean.RData")
data_b <- clean4

################
################
##
## TWIN STUDY 
##
################
################

#####################
#### FORMAT DATA ####
#####################

# order the data by ID and age
phenotype_order <- phenotype[order(phenotype$IID, phenotype$AGE),]
head(phenotype_order,10)

# within IID keep only the youngest age
# because longitudinal (several time points, we want one)
phenotype_order_single <- phenotype_order[!duplicated(phenotype_order$IID),]
head(phenotype_order_single,10)

# subset the methylation data
data_meth <- data[,colnames(data) %in% phenotype_order_single$SAMPLE]
# sex variables
table(phenotype_order_single$SEX)
# get female IDs
female_ID <- phenotype_order_single$SAMPLE[phenotype_order_single$SEX == 'female']
female_ID_num <- colnames(data_meth) %in% female_ID
# get male IDs
male_ID <- phenotype_order_single$SAMPLE[phenotype_order_single$SEX == 'male']
male_ID_num <- colnames(data_meth) %in% male_ID

sum(female_ID %in% phenotype_order_single$SAMPLE)
dim(data_meth[,colnames(data_meth) %in% female_ID])
dim(data_meth[,colnames(data_meth) %in% male_ID])

################################################################################################################################

# try the GAP STAT function
stat(data_meth["cg13869341",colnames(data_meth) %in% female_ID],data_meth["cg13869341",colnames(data_meth) %in% male_ID])

####################################
### APPLY GAP STAT TO TWIN STUDY ###
####################################

# gap_stat <- apply(data_meth, 1, function(x) stat(x[female_ID_num], x[male_ID_num]))
# scenario <- apply(data_meth, 1, function(x) as.numeric(mean(x[female_ID_num]) > mean(x[male_ID_num])))
# save(gap_stat, file="gap_stat.RData")
# save(scenario, file="scenario.RData")
load('gap_stat.RData') 
load('scenario.RData') 

hist(gap_stat, breaks=100)
table(gap_stat)
table(scenario)

# ## check difference in means 
# ## to compare with gap stat
# diff = NULL
# iter = 0
# for (i in c(70:93,95,96,98)) {
#   iter = iter+1
#   
#   num_gap = i
#   mean_girl = mean(data_meth[which(gap_stat == i)[1],female_ID_num]) # the first cg where gap_stat == 1
#   mean_boy = mean(data_meth[which(gap_stat == i)[1],male_ID_num]) # the first cg where gap_stat == 1
#   diff[iter] = abs(mean_girl - mean_boy)
# 
# }
# 
# # example for gap_stat = 98
# hist(data_meth[which(gap_stat == i)[1],female_ID_num],col="pink",breaks=100,main="",xlab="")
# hist(data_meth[which(gap_stat == i)[1],male_ID_num],col="blue",breaks=100,add=TRUE)
# 
# # check a range of difference in means for a range of gap_stats
# c(70:93,95,96,98)
# diff

###################################
### PLOT GAP STAT OF TWIN STUDY ###
###################################
gap <- 80

gap_stats_interest <- gap_stat[gap_stat >= gap]
gap_stats_ordered <- gap_stats_interest[order(gap_stats_interest)]

# store the cgs of interest
cgs <- names(gap_stats_ordered)
head(cgs)

# subset the data to have only the cgs of interest
cgs_interest = data.frame(data_meth[cgs,])
dim(cgs_interest)
head(cgs_interest)

cgs_interest$cg_name <- factor(rownames(cgs_interest), levels = rownames(cgs_interest))

# format data for ggplot
cgs_interest_melt <- melt(cgs_interest, id.vars = "cg_name", measure.vars = 1:444)
head(cgs_interest_melt,20)

# add sex variable 
cgs_interest_melt$sex <- "M"
cgs_interest_melt$sex[cgs_interest_melt$variable %in% female_ID] <- "F"

# add scenario variable
cgs_interest_melt <- merge(cgs_interest_melt, data.frame(cg_names = names(scenario),scenario),
                           by.x = 1, by.y = 1, all.x = TRUE)

# add chromosome, site, distance variables (from bob data) 
data_b$cg_name <- factor(rownames(data_b), levels = rownames(data_b))
cgs_interest_melt <- merge(cgs_interest_melt, data_b[,c("chr","site","distance","cg_name")],
                           by.x = "cg_name", by.y = "cg_name", all.x = TRUE)
cgs_interest_melt$chr <- factor(as.numeric(cgs_interest_melt$chr), levels = 1:22)

gap_value <- data.frame(stat = as.character(gap_stats_ordered), cg_name = names(gap_stats_ordered))
# gap_value <- merge(gap_value, data_b[,c("chr","cg_name")], by.x = "cg_name", by.y = "cg_name", all.x = TRUE)
# gap_value_1 <- gap_value[gap_value$cg_name %in% unique(cgs_interest_melt$cg_name[cgs_interest_melt$scenario == 1]),]
# gap_value_0 <- gap_value[gap_value$cg_name %in% unique(cgs_interest_melt$cg_name[cgs_interest_melt$scenario == 0]),]

g_1 <- ggplot(cgs_interest_melt[cgs_interest_melt$scenario == 1,], aes(x=value)) +
  geom_histogram(alpha=0.5, aes(fill=sex), bins = 100, position="identity") +
  facet_wrap(. ~ chr + cg_name, nrow = 7, scales="free_y") + xlab('x 100 (%)') +
  # geom_text(data = gap_value_1, mapping = aes(x = .5, y = 40, label = stat)) +
  ggtitle('Twin study - Female higher scenario') + theme_minimal() 

g_0 <- ggplot(cgs_interest_melt[cgs_interest_melt$scenario == 0,], aes(x=value)) +
  geom_histogram(alpha=0.5, aes(fill=sex), bins = 100, position="identity") +
  facet_wrap(. ~ chr + cg_name, nrow = 4, scales="free_y") + xlab('x 100 (%)') +
  # geom_text(data = gap_value_0, mapping = aes(x = .5, y = 40, label = stat)) +
  ggtitle('Twin study - Male higher scenario') + theme_minimal() 

ggsave(file = 'twin_scenario1.jpeg', g_1,
       dpi=300,
       width = 170,
       height = 260,
       units = "mm")

ggsave(file = 'twin_scenario0.jpeg', g_0,
       dpi=300,
       width = 170,
       height = 180,
       units = "mm")

   
################
################
##
## BOB DATA
##
################
################


###################################
### APPLY GAP STAT TO BOB DATA ###
##################################

# gap_stat_bob <- apply(data_b[,1:17], 1, function(x) stat(x[c(1,3)], x[c(2,4:17)])) # female vs male 
# save(gap_stat_bob, file="gap_stat_bob.RData")
load('gap_stat_bob.RData')

hist(gap_stat_bob,breaks=100)
table(gap_stat_bob)

num_gap = 100
hist(as.numeric(data_b[which(gap_stat_bob == num_gap)[7],c(1,3)]),col="pink",main="",xlab="",xlim=c(0,1))
hist(as.numeric(data_b[which(gap_stat_bob == num_gap)[7],c(2,4:17)]),col="blue",add=TRUE)

#############################################################
### PLOT GAP STAT OF BOB DATA THAT MATCH CGS OF INTEREST ###
############################################################

# use the same cgs as in twin study
bob_cgs_interest = data_b[cgs,]
dim(bob_cgs_interest)
head(bob_cgs_interest)

# add chromosome
bob_melt <- melt(bob_cgs_interest, id.vars = c("cg_name","chr"), measure.vars = 1:17)
bob_melt$chr <- factor(as.numeric(bob_melt$chr), levels = 1:22)
head(bob_melt,20)

# add sex variable 
bob_melt$sex <- "M"
bob_melt$sex[bob_melt$variable %in% c("X44","X34")] <- "F"

# add scenario variable
bob_melt <- merge(bob_melt, data.frame(cg_names = names(scenario),scenario),
                           by.x = 1, by.y = 1, all.x = TRUE)

gap_value_bob <- data.frame(stat = as.character(gap_stat_bob[cgs]), cg_name = cgs)

ggplot(bob_melt, aes(x=value)) +
  geom_histogram(alpha=0.5, aes(fill=sex), bins = 100, position="identity") +
  facet_wrap(. ~ chr + cg_name, nrow = 7) +
  # geom_text(data = gap_value_bob, mapping = aes(x = .5, y = 6, label = stat)) +
  ggtitle('Bob') 

g_1_bob <- ggplot(bob_melt[bob_melt$scenario == 1,], aes(x=value)) +
  geom_histogram(alpha=0.5, aes(fill=sex), bins = 100, position="identity") +
  facet_wrap(. ~ chr + cg_name, nrow = 7, scales="free_y") + xlab('x 100 (%)') +
  # geom_text(data = gap_value_1, mapping = aes(x = .5, y = 40, label = stat)) +
  ggtitle('Bob data - Female higher scenario') + theme_minimal() 

g_0_bob <- ggplot(bob_melt[bob_melt$scenario == 0,], aes(x=value)) +
  geom_histogram(alpha=0.5, aes(fill=sex), bins = 100, position="identity") +
  facet_wrap(. ~ chr + cg_name, nrow = 4, scales="free_y") + xlab('x 100 (%)') +
  # geom_text(data = gap_value_0, mapping = aes(x = .5, y = 40, label = stat)) +
  ggtitle('Bob data  - Male higher scenario') + theme_minimal() 

ggsave(file = 'bob_scenario1.jpeg', g_1_bob,
       dpi=300,
       width = 170,
       height = 260,
       units = "mm")

ggsave(file = 'bob_scenario0.jpeg', g_0_bob,
       dpi=300,
       width = 170,
       height = 180,
       units = "mm")





