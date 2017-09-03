library(afex)
library(data.table)
library(plyr)
library(BayesFactor)

####################################
#Orginal study (Heyman et al., 2015)
####################################

#Read in the data
## setwd("C:/Users/u0072526/Dropbox/Artikels/In preparation/Priming load replication/data original")

## If you run this script in non-English OS.
Sys.setlocale(category = "LC_ALL", locale = "English")

setwd("D:/core/Research/projects/OpenSci/checking_papers/WM_semantic_priming")
## original<-as.data.frame(rbindlist(lapply(list.files(), read.table, header=TRUE)))
unzip("Data original.zip")
original<-as.data.frame(rbindlist(lapply(list.files(pattern="*.txt"), read.table, header=TRUE)))
#read in all the original data
extract_participantnumber<-function(s) strsplit(s, "_")[[1]][1] #function to get participant number from the file name; the first number of the file name = the participant number
participants_original<-as.character(sapply(list.files(pattern = "*.txt"), extract_participantnumber)) #vector with all participant numbers
original$participant<-rep(participants_original,each=320) #add the participant numbers to the data (original); there are 320 observations per participant (i.e., the 320 stimuli), thus the number is repeated 320 times
original$block<-ifelse(original$groupNb<33,1,2) #add block numbers, 1 or 2; the first 32 groups (also called cycles in the manuscript)  are in block 1, the last 32 groups are in block 2
original$correct_response<-ifelse(original$type=="PSW","NONWORD","WORD") #add the correct response to every target stimulus; if the target is of the type PSW (pseudoword), the correct response is NONWORD, else it is WORD
original$acc<-ifelse(original$response==original$correct_response,1,0) #add the accuracy (i.e., whether the response given by the participant matches the correct response), 1 for correct, 0 for incorrect; if a response wasn't given before the 3000 ms deadline, the trial automatically terminated; this was also coded as an incorrect response

 
#Participant characteristics
extract_age<-function(s) strsplit(s, "_")[[1]][2] #function to extract the age of every participant; the second number in the file name = the age
age <- sapply(list.files(pattern = "*.txt"), extract_age)
mean(as.numeric(age)) #average age
extract_gender <- function(s) strsplit(s, "_")[[1]][3] #function to extract the gender of every participant; the third element in the file name = the gender (m for male, f for female)
gender <- sapply(list.files(pattern = "*.txt"), extract_gender)
table(gender) #distribution of gender

## Delete raw files after import the raw data
file.remove(list.files(pattern="*.txt"))

#Dot memory task
original_loads<-original[original$indexInGroup==1,] #retains only the data for the first prime-target pair per cycle; it results in a subset of the data with 64 rows per participant (1 row for each of the 64 cycles)
original_ttests<-ddply(original_loads,.(participant,load),summarize,mean=mean(hits),p_value=try(t.test(hits,mu=1)$p.value)) #calculates the average number of hits per participant in both load conditions and performs by-participant and by-load one sample t-tests on the hits to check whether participants are above chance (with chance being 1 hit); in some cases participants always gave the four correct locations, resulting in a perfect score, no variability, and thus some error messages when performing the t-tests; note that this is a slightly different procedure than in the original analyses: here, we take all 32 patterns into account, originally only cycles which featured at least one critical item were taken into account (N<=32); the results are very similar of course 
all(as.numeric(original_ttests$p_value) < .05, na.rm=T) #all p-values are smaller than .05
t.test(mean~load,paired=T,data=original_ttests) #difference between high and low load
tapply(original_ttests$mean,original_ttests$load,mean) #average number of hits per condition (collapsed across participants)

#Removing outliers and error responses
sort(tapply(original$acc,original$participant,mean)) #error responses per participant (most is 86%)
original<-original[original$type!="PSW" & original$type!="Filler_SYM", ] #remove fillers and pseudowords
original_temp<-original[original$acc==1 & original$RT>250,] #temporary datafile with only correct responses and RTs above 250 ms 
original$cutoff<-rep(as.vector(tapply(original_temp$RT,original_temp$participant,mean)+3*tapply(original_temp$RT,original_temp$participant,sd)),each=120) #calculate a cutoff value per participant (i.e., by-participant mean + 3* by-participant standard deviation, based on the temporary dataset which contains no errors and unusually fast responses) and replicate it 120 times (=number of crucial items and thus number of observations per participant in the updated "original" dataset)
original_final<-original[original$acc==1 & original$RT>250 & original$RT<=original$cutoff,]  #final dataframe used in the analyses (errors and outliers have been removed)
1-nrow(original_final)/nrow(original) #proportion of data removed (not counting fillers and pseudowords)

#Analyses of the response times
#Participant analysis with Greenhouse Geisser correction of degrees of freedom
original_anovaF1<-aov_car(RT~type*load*target_prime_relatedness*soa+Error(participant/type*load*target_prime_relatedness*soa),data=original_final, anova_table=list(correction = "GG", es= "pes"))
original_anovaF1

#Item analysis
original_anovaF2<-aov_car(RT~type*load*target_prime_relatedness*soa+Error(target/load*target_prime_relatedness*soa),data=original_final, anova_table=list(correction = "GG", es= "pes"))
original_anovaF2

#By-participant priming effects (Table 2 in Heyman et al., 2015)
original_final$conditions<-as.factor(paste(original_final$load,original_final$type,original_final$soa,original_final$target_prime_relatedness,sep=".")) #creates a condition variable by combining load, type, soa, and relatedness into one value (mainly for cosmetic purposes)
original_participant_means<-as.data.frame(tapply(original_final$RT,list(original_final$participant,original_final$conditions),mean)) #calculate by-participant and by-condition means
apply(original_participant_means,2,mean) #then collapse across participants

#t-tests for all priming effects per condition; they are grouped per row (FA first, then BA, then SYM pairs), from left to right (see Table 2 in Heyman et al., 2015)
t.test(original_participant_means$HIGH_LOAD.FA.LOW_SOA.TARGET_PRIME_UNRELATED,original_participant_means$HIGH_LOAD.FA.LOW_SOA.TARGET_PRIME_RELATED, paired=T)
t.test(original_participant_means$LOW_LOAD.FA.LOW_SOA.TARGET_PRIME_UNRELATED,original_participant_means$LOW_LOAD.FA.LOW_SOA.TARGET_PRIME_RELATED, paired=T)
t.test(original_participant_means$HIGH_LOAD.FA.HIGH_SOA.TARGET_PRIME_UNRELATED,original_participant_means$HIGH_LOAD.FA.HIGH_SOA.TARGET_PRIME_RELATED, paired=T)
t.test(original_participant_means$LOW_LOAD.FA.HIGH_SOA.TARGET_PRIME_UNRELATED,original_participant_means$LOW_LOAD.FA.HIGH_SOA.TARGET_PRIME_RELATED, paired=T)

t.test(original_participant_means$HIGH_LOAD.BA.LOW_SOA.TARGET_PRIME_UNRELATED,original_participant_means$HIGH_LOAD.BA.LOW_SOA.TARGET_PRIME_RELATED, paired=T)
t.test(original_participant_means$LOW_LOAD.BA.LOW_SOA.TARGET_PRIME_UNRELATED,original_participant_means$LOW_LOAD.BA.LOW_SOA.TARGET_PRIME_RELATED, paired=T)
t.test(original_participant_means$HIGH_LOAD.BA.HIGH_SOA.TARGET_PRIME_UNRELATED,original_participant_means$HIGH_LOAD.BA.HIGH_SOA.TARGET_PRIME_RELATED, paired=T)
t.test(original_participant_means$LOW_LOAD.BA.HIGH_SOA.TARGET_PRIME_UNRELATED,original_participant_means$LOW_LOAD.BA.HIGH_SOA.TARGET_PRIME_RELATED, paired=T)

t.test(original_participant_means$HIGH_LOAD.SYM.LOW_SOA.TARGET_PRIME_UNRELATED,original_participant_means$HIGH_LOAD.SYM.LOW_SOA.TARGET_PRIME_RELATED, paired=T)
t.test(original_participant_means$LOW_LOAD.SYM.LOW_SOA.TARGET_PRIME_UNRELATED,original_participant_means$LOW_LOAD.SYM.LOW_SOA.TARGET_PRIME_RELATED, paired=T)
t.test(original_participant_means$HIGH_LOAD.SYM.HIGH_SOA.TARGET_PRIME_UNRELATED,original_participant_means$HIGH_LOAD.SYM.HIGH_SOA.TARGET_PRIME_RELATED, paired=T)
t.test(original_participant_means$LOW_LOAD.SYM.HIGH_SOA.TARGET_PRIME_UNRELATED,original_participant_means$LOW_LOAD.SYM.HIGH_SOA.TARGET_PRIME_RELATED, paired=T)


#######################
#Replication in English
#######################
#Analyses are analogous, so no comments unless clarification is needed

## setwd("C:/Users/u0072526/Dropbox/Artikels/In preparation/Priming load replication/data replication English")
unzip("Experiment 1.zip")
replication_EN<-as.data.frame(rbindlist(lapply(list.files(pattern = "*.txt"), read.table, header=TRUE)))
participants_replication_EN<-as.character(sapply(list.files(pattern = "*.txt"), extract_participantnumber))
replication_EN$participant<-rep(participants_replication_EN,each=320)
replication_EN$block<-ifelse(replication_EN$groupNb<33,1,2)
replication_EN$correct_response<-ifelse(replication_EN$type=="PSW","NONWORD","WORD")
replication_EN$acc<-ifelse(replication_EN$response==replication_EN$correct_response,1,0)

#Participant characteristics
age <- sapply(list.files(pattern = "*.txt"), extract_age)
mean(as.numeric(age))
gender <- sapply(list.files(), extract_gender)
table(gender)

#Dot memory task
replication_EN_loads<-replication_EN[replication_EN$indexInGroup==1,]
replication_EN_ttests<-ddply(replication_EN_loads,.(participant,load),summarize,mean=mean(hits),p_value=try(t.test(hits,mu=1)$p.value))
all(as.numeric(replication_EN_ttests$p_value) < .05, na.rm=T) #not all p-values are smaller than .05
sort(as.numeric(replication_EN_ttests$p_value),decreasing = T) #one participant did not perform significantly above chance: participant 16 in the high load condition, see replication_EN_ttests; that participant will be removed from the analyses 
t.test(replication_EN_loads$hits[replication_EN_loads$participant==16 & replication_EN_loads$load=="HIGH_LOAD"],mu=1)
t.test(mean~load,paired=T,data=replication_EN_ttests) 
tapply(replication_EN_ttests$mean,replication_EN_ttests$load,mean)

## Remove raw text file after import raw data
file.remove(list.files(pattern="*.txt"))

#Removing outliers and error responses
sort(tapply(replication_EN$acc,replication_EN$participant,mean)) #participants 27, 66, 56, 74, 14, 5, 10, and 13 will be removed from the analyses because of high error rates (more than 15%); two participants, 41 and 30 had very high error rates (97%) suggesting that they reversed the response options => their responses will be recoded; participant 27 also had a high error rate (100%), but was removed from the analyses because s/he never responded within the 3000 ms deadline
remove_participants_replication_EN<-c(16,27,66,56,74,14,5,10,13) #participant numbers that will be removed from the analyses
replication_EN$acc[replication_EN$participant==41|replication_EN$participant==30]<-ifelse(replication_EN$acc[replication_EN$participant==41|replication_EN$participant==30]==0,1,0) #accuracies for participants 41 and 30 are reverse scored
replication_EN<-replication_EN[replication_EN$type!="PSW" & replication_EN$type!="Filler_SYM" & !replication_EN$participant %in% remove_participants_replication_EN, ] #remove fillers and pseudowords + participants with many errors
replication_EN_temp<-replication_EN[replication_EN$acc==1 & replication_EN$RT>250,]
replication_EN$cutoff<-rep(as.vector(tapply(replication_EN_temp$RT,replication_EN_temp$participant,mean)+3*tapply(replication_EN_temp$RT,replication_EN_temp$participant,sd)),each=120)
replication_EN_final<-replication_EN[replication_EN$acc==1 & replication_EN$RT>250 & replication_EN$RT<=replication_EN$cutoff,]  
1-nrow(replication_EN_final)/nrow(replication_EN)

#Analyses of the response times
#Participant analysis
#Analysis gave a warning message, because there were missing values for participant 68 => also removed from the analyses
replication_EN_final<-replication_EN_final[replication_EN_final$participant!=68,]
replication_EN_anovaF1<-aov_car(RT~type*load*target_prime_relatedness*soa+Error(participant/type*load*target_prime_relatedness*soa),data=replication_EN_final, anova_table=list(correction = "GG", es= "pes"))
replication_EN_anovaF1

#Item analysis
replication_EN_anovaF2<-aov_car(RT~type*load*target_prime_relatedness*soa+Error(target/load*target_prime_relatedness*soa),data=replication_EN_final, anova_table=list(correction = "GG", es= "pes"))
replication_EN_anovaF2

#By-participant priming effects (Table 2 in Heyman et al., 2015)
replication_EN_final$conditions<-as.factor(paste(replication_EN_final$load,replication_EN_final$type,replication_EN_final$soa,replication_EN_final$target_prime_relatedness,sep="."))
replication_EN_participant_means<-as.data.frame(tapply(replication_EN_final$RT,list(replication_EN_final$participant,replication_EN_final$conditions),mean))
apply(replication_EN_participant_means,2,mean)

t.test(replication_EN_participant_means$HIGH_LOAD.FA.LOW_SOA.TARGET_PRIME_UNRELATED,replication_EN_participant_means$HIGH_LOAD.FA.LOW_SOA.TARGET_PRIME_RELATED, paired=T)
t.test(replication_EN_participant_means$LOW_LOAD.FA.LOW_SOA.TARGET_PRIME_UNRELATED,replication_EN_participant_means$LOW_LOAD.FA.LOW_SOA.TARGET_PRIME_RELATED, paired=T)
t.test(replication_EN_participant_means$HIGH_LOAD.FA.HIGH_SOA.TARGET_PRIME_UNRELATED,replication_EN_participant_means$HIGH_LOAD.FA.HIGH_SOA.TARGET_PRIME_RELATED, paired=T)
t.test(replication_EN_participant_means$LOW_LOAD.FA.HIGH_SOA.TARGET_PRIME_UNRELATED,replication_EN_participant_means$LOW_LOAD.FA.HIGH_SOA.TARGET_PRIME_RELATED, paired=T)

t.test(replication_EN_participant_means$HIGH_LOAD.BA.LOW_SOA.TARGET_PRIME_UNRELATED,replication_EN_participant_means$HIGH_LOAD.BA.LOW_SOA.TARGET_PRIME_RELATED, paired=T)
t.test(replication_EN_participant_means$LOW_LOAD.BA.LOW_SOA.TARGET_PRIME_UNRELATED,replication_EN_participant_means$LOW_LOAD.BA.LOW_SOA.TARGET_PRIME_RELATED, paired=T)
t.test(replication_EN_participant_means$HIGH_LOAD.BA.HIGH_SOA.TARGET_PRIME_UNRELATED,replication_EN_participant_means$HIGH_LOAD.BA.HIGH_SOA.TARGET_PRIME_RELATED, paired=T)
t.test(replication_EN_participant_means$LOW_LOAD.BA.HIGH_SOA.TARGET_PRIME_UNRELATED,replication_EN_participant_means$LOW_LOAD.BA.HIGH_SOA.TARGET_PRIME_RELATED, paired=T)

t.test(replication_EN_participant_means$HIGH_LOAD.SYM.LOW_SOA.TARGET_PRIME_UNRELATED,replication_EN_participant_means$HIGH_LOAD.SYM.LOW_SOA.TARGET_PRIME_RELATED, paired=T)
t.test(replication_EN_participant_means$LOW_LOAD.SYM.LOW_SOA.TARGET_PRIME_UNRELATED,replication_EN_participant_means$LOW_LOAD.SYM.LOW_SOA.TARGET_PRIME_RELATED, paired=T)
t.test(replication_EN_participant_means$HIGH_LOAD.SYM.HIGH_SOA.TARGET_PRIME_UNRELATED,replication_EN_participant_means$HIGH_LOAD.SYM.HIGH_SOA.TARGET_PRIME_RELATED, paired=T)
t.test(replication_EN_participant_means$LOW_LOAD.SYM.HIGH_SOA.TARGET_PRIME_UNRELATED,replication_EN_participant_means$LOW_LOAD.SYM.HIGH_SOA.TARGET_PRIME_RELATED, paired=T)


#####################
#Replication in Dutch
#####################

#Read in the data
##setwd("C:/Users/u0072526/Dropbox/Artikels/In preparation/Priming load replication/data replication Dutch")
unzip("Experiment 2.zip")
replication_DU<-as.data.frame(rbindlist(lapply(list.files(pattern = "*.txt"), read.table, header=TRUE)))
participants_replication_DU<-as.character(sapply(list.files(pattern = "*.txt"), extract_participantnumber))
replication_DU$participant<-rep(participants_replication_DU,each=320)
replication_DU$block<-ifelse(replication_DU$groupNb<33,1,2)
replication_DU$correct_response<-ifelse(replication_DU$type=="PSW","NONWORD","WORD")
replication_DU$acc<-ifelse(replication_DU$response==replication_DU$correct_response,1,0)

#Participant characteristics
age <- sapply(list.files(pattern = "*.txt"), extract_age)
mean(as.numeric(age))
gender <- sapply(list.files(pattern = "*.txt"), extract_gender)
table(gender)

#Dot memory task
replication_DU_loads<-replication_DU[replication_DU$indexInGroup==1,]
replication_DU_ttests<-ddply(replication_DU_loads,.(participant,load),summarize,mean=mean(hits),p_value=try(t.test(hits,mu=1)$p.value))
all(as.numeric(replication_DU_ttests$p_value) < .05, na.rm=T) 
t.test(mean~load,paired=T,data=replication_DU_ttests) 
tapply(replication_DU_ttests$mean,replication_DU_ttests$load,mean)

#Removing outliers and error responses
sort(tapply(replication_DU$acc,replication_DU$participant,mean)) #Participants 13, 162, 36, and 43 will be removed from the analyses because of high error rates (more than 15%).
remove_participants_replication_DU<-c(13,162,36,43)
replication_DU<-replication_DU[replication_DU$type!="PSW" & replication_DU$type!="Filler_SYM" & !replication_DU$participant %in% remove_participants_replication_DU, ] #remove fillers and pseudowords + participants with many errors
replication_DU_temp<-replication_DU[replication_DU$acc==1 & replication_DU$RT>250,]
replication_DU$cutoff<-rep(as.vector(tapply(replication_DU_temp$RT,replication_DU_temp$participant,mean)+3*tapply(replication_DU_temp$RT,replication_DU_temp$participant,sd)),each=120)
replication_DU_final<-replication_DU[replication_DU$acc==1 & replication_DU$RT>250 & replication_DU$RT<=replication_DU$cutoff,]  
1-nrow(replication_DU_final)/nrow(replication_DU)


## Remove raw text file after import raw data
file.remove(list.files(pattern="*.txt"))


#Analyses of the response times
#Participant analysis
replication_DU_anovaF1<-aov_car(RT~type*load*target_prime_relatedness*soa+Error(participant/type*load*target_prime_relatedness*soa),data=replication_DU_final, anova_table=list(correction = "GG", es= "pes"))
replication_DU_anovaF1

#Item analysis
replication_DU_anovaF2<-aov_car(RT~type*load*target_prime_relatedness*soa+Error(target/load*target_prime_relatedness*soa),data=replication_DU_final, anova_table=list(correction = "GG", es= "pes"))
replication_DU_anovaF2

#By-participant priming effects (Table 2 in Heyman et al., 2015)
replication_DU_final$conditions<-as.factor(paste(replication_DU_final$load,replication_DU_final$type,replication_DU_final$soa,replication_DU_final$target_prime_relatedness,sep="."))
replication_DU_participant_means<-as.data.frame(tapply(replication_DU_final$RT,list(replication_DU_final$participant,replication_DU_final$conditions),mean))
apply(replication_DU_participant_means,2,mean)

t.test(replication_DU_participant_means$HIGH_LOAD.FA.LOW_SOA.TARGET_PRIME_UNRELATED,replication_DU_participant_means$HIGH_LOAD.FA.LOW_SOA.TARGET_PRIME_RELATED, paired=T)
t.test(replication_DU_participant_means$LOW_LOAD.FA.LOW_SOA.TARGET_PRIME_UNRELATED,replication_DU_participant_means$LOW_LOAD.FA.LOW_SOA.TARGET_PRIME_RELATED, paired=T)
t.test(replication_DU_participant_means$HIGH_LOAD.FA.HIGH_SOA.TARGET_PRIME_UNRELATED,replication_DU_participant_means$HIGH_LOAD.FA.HIGH_SOA.TARGET_PRIME_RELATED, paired=T)
t.test(replication_DU_participant_means$LOW_LOAD.FA.HIGH_SOA.TARGET_PRIME_UNRELATED,replication_DU_participant_means$LOW_LOAD.FA.HIGH_SOA.TARGET_PRIME_RELATED, paired=T)

t.test(replication_DU_participant_means$HIGH_LOAD.BA.LOW_SOA.TARGET_PRIME_UNRELATED,replication_DU_participant_means$HIGH_LOAD.BA.LOW_SOA.TARGET_PRIME_RELATED, paired=T)
t.test(replication_DU_participant_means$LOW_LOAD.BA.LOW_SOA.TARGET_PRIME_UNRELATED,replication_DU_participant_means$LOW_LOAD.BA.LOW_SOA.TARGET_PRIME_RELATED, paired=T)
t.test(replication_DU_participant_means$HIGH_LOAD.BA.HIGH_SOA.TARGET_PRIME_UNRELATED,replication_DU_participant_means$HIGH_LOAD.BA.HIGH_SOA.TARGET_PRIME_RELATED, paired=T)
t.test(replication_DU_participant_means$LOW_LOAD.BA.HIGH_SOA.TARGET_PRIME_UNRELATED,replication_DU_participant_means$LOW_LOAD.BA.HIGH_SOA.TARGET_PRIME_RELATED, paired=T)

t.test(replication_DU_participant_means$HIGH_LOAD.SYM.LOW_SOA.TARGET_PRIME_UNRELATED,replication_DU_participant_means$HIGH_LOAD.SYM.LOW_SOA.TARGET_PRIME_RELATED, paired=T)
t.test(replication_DU_participant_means$LOW_LOAD.SYM.LOW_SOA.TARGET_PRIME_UNRELATED,replication_DU_participant_means$LOW_LOAD.SYM.LOW_SOA.TARGET_PRIME_RELATED, paired=T)
t.test(replication_DU_participant_means$HIGH_LOAD.SYM.HIGH_SOA.TARGET_PRIME_UNRELATED,replication_DU_participant_means$HIGH_LOAD.SYM.HIGH_SOA.TARGET_PRIME_RELATED, paired=T)
t.test(replication_DU_participant_means$LOW_LOAD.SYM.HIGH_SOA.TARGET_PRIME_UNRELATED,replication_DU_participant_means$LOW_LOAD.SYM.HIGH_SOA.TARGET_PRIME_RELATED, paired=T)

####################
#Joint analyses
####################

#Analysis combining the original data with those from the Dutch (exact) replication
replication_DU_final$participant<-as.numeric(replication_DU_final$participant)+80 #to make unique participant numbers
replication_DU_final$experiment<-"rep"
original_final$experiment<-"orig"
joint<-rbind(original_final,replication_DU_final)
joint$participant<-as.factor(joint$participant)
joint$experiment<-as.factor(joint$experiment)

#Participant analysis
joint_anovaF1<-aov_car(RT~type*load*target_prime_relatedness*soa*experiment+Error(participant/type*load*target_prime_relatedness*soa),data=joint, anova_table=list(correction = "GG", es= "pes"))
joint_anovaF1

#Item analysis
joint_anovaF2<-aov_car(RT~type*load*target_prime_relatedness*soa*experiment+Error(target/load*target_prime_relatedness*soa*experiment),data=joint, anova_table=list(correction = "GG", es= "pes"))
joint_anovaF2


#Calculate by-participant priming effects per condition based on the original data
original_PE_high_load_FA_low_SOA<-original_participant_means$HIGH_LOAD.FA.LOW_SOA.TARGET_PRIME_UNRELATED-original_participant_means$HIGH_LOAD.FA.LOW_SOA.TARGET_PRIME_RELATED
original_PE_low_load_FA_low_SOA<-original_participant_means$LOW_LOAD.FA.LOW_SOA.TARGET_PRIME_UNRELATED-original_participant_means$LOW_LOAD.FA.LOW_SOA.TARGET_PRIME_RELATED
original_PE_high_load_FA_high_SOA<-original_participant_means$HIGH_LOAD.FA.HIGH_SOA.TARGET_PRIME_UNRELATED-original_participant_means$HIGH_LOAD.FA.HIGH_SOA.TARGET_PRIME_RELATED
original_PE_low_load_FA_high_SOA<-original_participant_means$LOW_LOAD.FA.HIGH_SOA.TARGET_PRIME_UNRELATED-original_participant_means$LOW_LOAD.FA.HIGH_SOA.TARGET_PRIME_RELATED

original_PE_high_load_BA_low_SOA<-original_participant_means$HIGH_LOAD.BA.LOW_SOA.TARGET_PRIME_UNRELATED-original_participant_means$HIGH_LOAD.BA.LOW_SOA.TARGET_PRIME_RELATED
original_PE_low_load_BA_low_SOA<-original_participant_means$LOW_LOAD.BA.LOW_SOA.TARGET_PRIME_UNRELATED-original_participant_means$LOW_LOAD.BA.LOW_SOA.TARGET_PRIME_RELATED
original_PE_high_load_BA_high_SOA<-original_participant_means$HIGH_LOAD.BA.HIGH_SOA.TARGET_PRIME_UNRELATED-original_participant_means$HIGH_LOAD.BA.HIGH_SOA.TARGET_PRIME_RELATED
original_PE_low_load_BA_high_SOA<-original_participant_means$LOW_LOAD.BA.HIGH_SOA.TARGET_PRIME_UNRELATED-original_participant_means$LOW_LOAD.BA.HIGH_SOA.TARGET_PRIME_RELATED

original_PE_high_load_SYM_low_SOA<-original_participant_means$HIGH_LOAD.SYM.LOW_SOA.TARGET_PRIME_UNRELATED-original_participant_means$HIGH_LOAD.SYM.LOW_SOA.TARGET_PRIME_RELATED
original_PE_low_load_SYM_low_SOA<-original_participant_means$LOW_LOAD.SYM.LOW_SOA.TARGET_PRIME_UNRELATED-original_participant_means$LOW_LOAD.SYM.LOW_SOA.TARGET_PRIME_RELATED
original_PE_high_load_SYM_high_SOA<-original_participant_means$HIGH_LOAD.SYM.HIGH_SOA.TARGET_PRIME_UNRELATED-original_participant_means$HIGH_LOAD.SYM.HIGH_SOA.TARGET_PRIME_RELATED
original_PE_low_load_SYM_high_SOA<-original_participant_means$LOW_LOAD.SYM.HIGH_SOA.TARGET_PRIME_UNRELATED-original_participant_means$LOW_LOAD.SYM.HIGH_SOA.TARGET_PRIME_RELATED

#Bayesian one-sample t-tests on the resulting priming effects
ttestBF(original_PE_high_load_FA_low_SOA)
ttestBF(original_PE_low_load_FA_low_SOA)
ttestBF(original_PE_high_load_FA_high_SOA)
ttestBF(original_PE_low_load_FA_high_SOA)

ttestBF(original_PE_high_load_BA_low_SOA)
ttestBF(original_PE_low_load_BA_low_SOA)
ttestBF(original_PE_high_load_BA_high_SOA)
ttestBF(original_PE_low_load_BA_high_SOA)

ttestBF(original_PE_high_load_SYM_low_SOA)
ttestBF(original_PE_low_load_SYM_low_SOA)
ttestBF(original_PE_high_load_SYM_high_SOA)
ttestBF(original_PE_low_load_SYM_high_SOA)


#Calculate by-participant priming effects per condition based on all the data
all_PE_high_load_FA_low_SOA<-c(original_participant_means$HIGH_LOAD.FA.LOW_SOA.TARGET_PRIME_UNRELATED,replication_EN_participant_means$HIGH_LOAD.FA.LOW_SOA.TARGET_PRIME_UNRELATED,replication_DU_participant_means$HIGH_LOAD.FA.LOW_SOA.TARGET_PRIME_UNRELATED)-c(original_participant_means$HIGH_LOAD.FA.LOW_SOA.TARGET_PRIME_RELATED,replication_EN_participant_means$HIGH_LOAD.FA.LOW_SOA.TARGET_PRIME_RELATED,replication_DU_participant_means$HIGH_LOAD.FA.LOW_SOA.TARGET_PRIME_RELATED)
all_PE_low_load_FA_low_SOA<-c(original_participant_means$LOW_LOAD.FA.LOW_SOA.TARGET_PRIME_UNRELATED,replication_EN_participant_means$LOW_LOAD.FA.LOW_SOA.TARGET_PRIME_UNRELATED,replication_DU_participant_means$LOW_LOAD.FA.LOW_SOA.TARGET_PRIME_UNRELATED)-c(original_participant_means$LOW_LOAD.FA.LOW_SOA.TARGET_PRIME_RELATED,replication_EN_participant_means$LOW_LOAD.FA.LOW_SOA.TARGET_PRIME_RELATED,replication_DU_participant_means$LOW_LOAD.FA.LOW_SOA.TARGET_PRIME_RELATED)
all_PE_high_load_FA_high_SOA<-c(original_participant_means$HIGH_LOAD.FA.HIGH_SOA.TARGET_PRIME_UNRELATED,replication_EN_participant_means$HIGH_LOAD.FA.HIGH_SOA.TARGET_PRIME_UNRELATED,replication_DU_participant_means$HIGH_LOAD.FA.HIGH_SOA.TARGET_PRIME_UNRELATED)-c(original_participant_means$HIGH_LOAD.FA.HIGH_SOA.TARGET_PRIME_RELATED,replication_EN_participant_means$HIGH_LOAD.FA.HIGH_SOA.TARGET_PRIME_RELATED,replication_DU_participant_means$HIGH_LOAD.FA.HIGH_SOA.TARGET_PRIME_RELATED)
all_PE_low_load_FA_high_SOA<-c(original_participant_means$LOW_LOAD.FA.HIGH_SOA.TARGET_PRIME_UNRELATED,replication_EN_participant_means$LOW_LOAD.FA.HIGH_SOA.TARGET_PRIME_UNRELATED,replication_DU_participant_means$LOW_LOAD.FA.HIGH_SOA.TARGET_PRIME_UNRELATED)-c(original_participant_means$LOW_LOAD.FA.HIGH_SOA.TARGET_PRIME_RELATED,replication_EN_participant_means$LOW_LOAD.FA.HIGH_SOA.TARGET_PRIME_RELATED,replication_DU_participant_means$LOW_LOAD.FA.HIGH_SOA.TARGET_PRIME_RELATED)

all_PE_high_load_BA_low_SOA<-c(original_participant_means$HIGH_LOAD.BA.LOW_SOA.TARGET_PRIME_UNRELATED,replication_EN_participant_means$HIGH_LOAD.BA.LOW_SOA.TARGET_PRIME_UNRELATED,replication_DU_participant_means$HIGH_LOAD.BA.LOW_SOA.TARGET_PRIME_UNRELATED)-c(original_participant_means$HIGH_LOAD.BA.LOW_SOA.TARGET_PRIME_RELATED,replication_EN_participant_means$HIGH_LOAD.BA.LOW_SOA.TARGET_PRIME_RELATED,replication_DU_participant_means$HIGH_LOAD.BA.LOW_SOA.TARGET_PRIME_RELATED)
all_PE_low_load_BA_low_SOA<-c(original_participant_means$LOW_LOAD.BA.LOW_SOA.TARGET_PRIME_UNRELATED,replication_EN_participant_means$LOW_LOAD.BA.LOW_SOA.TARGET_PRIME_UNRELATED,replication_DU_participant_means$LOW_LOAD.BA.LOW_SOA.TARGET_PRIME_UNRELATED)-c(original_participant_means$LOW_LOAD.BA.LOW_SOA.TARGET_PRIME_RELATED,replication_EN_participant_means$LOW_LOAD.BA.LOW_SOA.TARGET_PRIME_RELATED,replication_DU_participant_means$LOW_LOAD.BA.LOW_SOA.TARGET_PRIME_RELATED)
all_PE_high_load_BA_high_SOA<-c(original_participant_means$HIGH_LOAD.BA.HIGH_SOA.TARGET_PRIME_UNRELATED,replication_EN_participant_means$HIGH_LOAD.BA.HIGH_SOA.TARGET_PRIME_UNRELATED,replication_DU_participant_means$HIGH_LOAD.BA.HIGH_SOA.TARGET_PRIME_UNRELATED)-c(original_participant_means$HIGH_LOAD.BA.HIGH_SOA.TARGET_PRIME_RELATED,replication_EN_participant_means$HIGH_LOAD.BA.HIGH_SOA.TARGET_PRIME_RELATED,replication_DU_participant_means$HIGH_LOAD.BA.HIGH_SOA.TARGET_PRIME_RELATED)
all_PE_low_load_BA_high_SOA<-c(original_participant_means$LOW_LOAD.BA.HIGH_SOA.TARGET_PRIME_UNRELATED,replication_EN_participant_means$LOW_LOAD.BA.HIGH_SOA.TARGET_PRIME_UNRELATED,replication_DU_participant_means$LOW_LOAD.BA.HIGH_SOA.TARGET_PRIME_UNRELATED)-c(original_participant_means$LOW_LOAD.BA.HIGH_SOA.TARGET_PRIME_RELATED,replication_EN_participant_means$LOW_LOAD.BA.HIGH_SOA.TARGET_PRIME_RELATED,replication_DU_participant_means$LOW_LOAD.BA.HIGH_SOA.TARGET_PRIME_RELATED)

all_PE_high_load_SYM_low_SOA<-c(original_participant_means$HIGH_LOAD.SYM.LOW_SOA.TARGET_PRIME_UNRELATED,replication_EN_participant_means$HIGH_LOAD.SYM.LOW_SOA.TARGET_PRIME_UNRELATED,replication_DU_participant_means$HIGH_LOAD.SYM.LOW_SOA.TARGET_PRIME_UNRELATED)-c(original_participant_means$HIGH_LOAD.SYM.LOW_SOA.TARGET_PRIME_RELATED,replication_EN_participant_means$HIGH_LOAD.SYM.LOW_SOA.TARGET_PRIME_RELATED,replication_DU_participant_means$HIGH_LOAD.SYM.LOW_SOA.TARGET_PRIME_RELATED)
all_PE_low_load_SYM_low_SOA<-c(original_participant_means$LOW_LOAD.SYM.LOW_SOA.TARGET_PRIME_UNRELATED,replication_EN_participant_means$LOW_LOAD.SYM.LOW_SOA.TARGET_PRIME_UNRELATED,replication_DU_participant_means$LOW_LOAD.SYM.LOW_SOA.TARGET_PRIME_UNRELATED)-c(original_participant_means$LOW_LOAD.SYM.LOW_SOA.TARGET_PRIME_RELATED,replication_EN_participant_means$LOW_LOAD.SYM.LOW_SOA.TARGET_PRIME_RELATED,replication_DU_participant_means$LOW_LOAD.SYM.LOW_SOA.TARGET_PRIME_RELATED)
all_PE_high_load_SYM_high_SOA<-c(original_participant_means$HIGH_LOAD.SYM.HIGH_SOA.TARGET_PRIME_UNRELATED,replication_EN_participant_means$HIGH_LOAD.SYM.HIGH_SOA.TARGET_PRIME_UNRELATED,replication_DU_participant_means$HIGH_LOAD.SYM.HIGH_SOA.TARGET_PRIME_UNRELATED)-c(original_participant_means$HIGH_LOAD.SYM.HIGH_SOA.TARGET_PRIME_RELATED,replication_EN_participant_means$HIGH_LOAD.SYM.HIGH_SOA.TARGET_PRIME_RELATED,replication_DU_participant_means$HIGH_LOAD.SYM.HIGH_SOA.TARGET_PRIME_RELATED)
all_PE_low_load_SYM_high_SOA<-c(original_participant_means$LOW_LOAD.SYM.HIGH_SOA.TARGET_PRIME_UNRELATED,replication_EN_participant_means$LOW_LOAD.SYM.HIGH_SOA.TARGET_PRIME_UNRELATED,replication_DU_participant_means$LOW_LOAD.SYM.HIGH_SOA.TARGET_PRIME_UNRELATED)-c(original_participant_means$LOW_LOAD.SYM.HIGH_SOA.TARGET_PRIME_RELATED,replication_EN_participant_means$LOW_LOAD.SYM.HIGH_SOA.TARGET_PRIME_RELATED,replication_DU_participant_means$LOW_LOAD.SYM.HIGH_SOA.TARGET_PRIME_RELATED)

#Bayesian and frequentist one-sample t-tests on the resulting priming effects (all priming effects are considered to come from one big experiment)
t.test(all_PE_high_load_FA_low_SOA)
ttestBF(all_PE_high_load_FA_low_SOA)
t.test(all_PE_low_load_FA_low_SOA)
ttestBF(all_PE_low_load_FA_low_SOA)
t.test(all_PE_high_load_FA_high_SOA)
ttestBF(all_PE_high_load_FA_high_SOA)
t.test(all_PE_low_load_FA_high_SOA)
ttestBF(all_PE_low_load_FA_high_SOA)

t.test(all_PE_high_load_BA_low_SOA)
ttestBF(all_PE_high_load_BA_low_SOA)
t.test(all_PE_low_load_BA_low_SOA)
ttestBF(all_PE_low_load_BA_low_SOA)
t.test(all_PE_high_load_BA_high_SOA)
ttestBF(all_PE_high_load_BA_high_SOA)
t.test(all_PE_low_load_BA_high_SOA)
ttestBF(all_PE_low_load_BA_high_SOA)

t.test(all_PE_high_load_SYM_low_SOA)
ttestBF(all_PE_high_load_SYM_low_SOA)
t.test(all_PE_low_load_SYM_low_SOA)
ttestBF(all_PE_low_load_SYM_low_SOA)
t.test(all_PE_high_load_SYM_high_SOA)
ttestBF(all_PE_high_load_SYM_high_SOA)
t.test(all_PE_low_load_SYM_high_SOA)
ttestBF(all_PE_low_load_SYM_high_SOA)

#Bayesian paired samples t-tests on the original data to evaluate the (potential) impact of the load manipulation on the conditional priming effects 
ttestBF(original_PE_low_load_FA_low_SOA,original_PE_high_load_FA_low_SOA,paired=T)
ttestBF(original_PE_low_load_FA_high_SOA,original_PE_high_load_FA_high_SOA,paired=T)

ttestBF(original_PE_low_load_BA_low_SOA,original_PE_high_load_BA_low_SOA,paired=T)
ttestBF(original_PE_low_load_BA_high_SOA,original_PE_high_load_BA_high_SOA,paired=T)

ttestBF(original_PE_low_load_SYM_low_SOA,original_PE_high_load_SYM_low_SOA,paired=T)
ttestBF(original_PE_low_load_SYM_high_SOA,original_PE_high_load_SYM_high_SOA,paired=T)

#Bayesian and frequentist paired samples t-tests on all the data to evaluate the (potential) impact of the load manipulation on the conditional priming effects 
ttestBF(all_PE_low_load_FA_low_SOA,all_PE_high_load_FA_low_SOA,paired=T)
t.test(all_PE_low_load_FA_low_SOA,all_PE_high_load_FA_low_SOA,paired=T)
ttestBF(all_PE_low_load_FA_high_SOA,all_PE_high_load_FA_high_SOA,paired=T)
t.test(all_PE_low_load_FA_high_SOA,all_PE_high_load_FA_high_SOA,paired=T)

ttestBF(all_PE_low_load_BA_low_SOA,all_PE_high_load_BA_low_SOA,paired=T)
t.test(all_PE_low_load_BA_low_SOA,all_PE_high_load_BA_low_SOA,paired=T)
ttestBF(all_PE_low_load_BA_high_SOA,all_PE_high_load_BA_high_SOA,paired=T)
t.test(all_PE_low_load_BA_high_SOA,all_PE_high_load_BA_high_SOA,paired=T)

ttestBF(all_PE_low_load_SYM_low_SOA,all_PE_high_load_SYM_low_SOA,paired=T)
t.test(all_PE_low_load_SYM_low_SOA,all_PE_high_load_SYM_low_SOA,paired=T)
ttestBF(all_PE_low_load_SYM_high_SOA,all_PE_high_load_SYM_high_SOA,paired=T)
t.test(all_PE_low_load_SYM_high_SOA,all_PE_high_load_SYM_high_SOA,paired=T)

