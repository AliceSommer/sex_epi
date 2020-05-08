library(reshape2)
library(ggplot2)
library(Biobase)

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

# methylation data
load('betas_preprocessed.RData')
data <- betas.preprocessed
# metadata
phenotype <- read.csv('phenotype.csv')
head(phenotype)

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
# save(gap_stat, file="gap_stat.RData")
load('gap_stat.RData') 

hist(gap_stat, breaks=100)
table(gap_stat)

## check some difference in 
diff = NULL
iter = 0
for (i in c(70:93,95,96,98)) {
  iter = iter+1
  
  num_gap = i
  mean_girl = mean(data_meth[which(gap_stat == i)[1],female_ID_num]) # the first cg where gap_stat == 1
  mean_boy = mean(data_meth[which(gap_stat == i)[1],male_ID_num]) # the first cg where gap_stat == 1
  diff[iter] = abs(mean_girl - mean_boy)

}

# example for gap_stat = 98
hist(data_meth[which(gap_stat == i)[1],female_ID_num],col="pink",breaks=100,main="",xlab="")
hist(data_meth[which(gap_stat == i)[1],male_ID_num],col="blue",breaks=100,add=TRUE)

# check a gange of difference in means for a range of gap_stats
c(70:93,95,96,98)
diff

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
scenario_1 = data.frame(data_meth[cgs,])
dim(scenario_1)
head(scenario_1)

scenario_1$cg_name <- factor(rownames(scenario_1), levels = rownames(scenario_1))

scenario_1_melt <- melt(scenario_1, id.vars = "cg_name", measure.vars = 1:444)
head(scenario_1_melt,20)

scenario_1_melt$sex <- "M"
scenario_1_melt$sex[scenario_1_melt$variable %in% female_ID] <- "F"

gap_value <- data.frame(stat = as.character(gap_stats_ordered), cg_name = names(gap_stats_ordered))

ggplot(scenario_1_melt, aes(x=value)) +
  geom_histogram(alpha=0.5, aes(fill=sex), bins = 100, position="identity") +
  facet_wrap(. ~ cg_name, nrow = 7) +
  geom_text(data = gap_value, mapping = aes(x = .5, y = 40, label = stat)) +
  ggtitle('Twin') 
   

################
################
##
## BOB DATA
##
################
################

#####################
#### FORMAT DATA #### 
#####################
# doesn't have to be formated because Twin data formated to match Bob's

load("BOB_data/epigen_ozone/clean.RData")

###################################
### APPLY GAP STAT TO BOB DATA ###
##################################

# gap_stat_bob <- apply(clean4[,1:17], 1, function(x) stat(x[c(1,3)], x[c(2,4:17)])) # female vs male 
# save(gap_stat_bob, file="gap_stat_bob.RData")
load('gap_stat_bob.RData')

hist(gap_stat_bob,breaks=100)
table(gap_stat_bob)

num_gap = 100
hist(as.numeric(clean4[which(gap_stat_bob == num_gap)[7],c(1,3)]),col="pink",main="",xlab="",xlim=c(0,1))
hist(as.numeric(clean4[which(gap_stat_bob == num_gap)[7],c(2,4:17)]),col="blue",add=TRUE)

#############################################################
### PLOT GAP STAT OF BOB DATA THAT MATCH CGS OF INTEREST ###
############################################################

# use the same cgs as in twin study
scenario_1x = clean4[cgs,]
dim(scenario_1x)
head(scenario_1x)

scenario_1x$cg_name <- factor(rownames(scenario_1x), levels = rownames(scenario_1))

scenario_1x_melt <- melt(scenario_1x, id.vars = c("cg_name","chr"), measure.vars = 1:17)
head(scenario_1x_melt,20)

scenario_1x_melt$sex <- "M"
scenario_1x_melt$sex[scenario_1x_melt$variable %in% c("X44","X34")] <- "F"

gap_value_bob <- data.frame(stat = as.character(gap_stat_bob[cgs]), cg_name = cgs)

ggplot(scenario_1x_melt, aes(x=value)) +
  geom_histogram(alpha=0.5, aes(fill=sex), bins = 100, position="identity") +
  facet_wrap(. ~ cg_name, nrow = 7) +
  geom_text(data = gap_value_bob, mapping = aes(x = .5, y = 6, label = stat)) +
  ggtitle('Bob') 

tail(colnames(clean4))

table(scenario_1x[,"chr"])



