library(MOTE)
library(TOSTER)

## Import the priming efffects data
orig_DU_comparisons <- read.csv("orig_DU.csv")
repli_DU_comparisons <- read.csv("repli_DU.csv")
repli_EN_comparisons <- read.csv("repli_EN.csv")

## original_study
## mean RT
apply(orig_DU_comparisons[,c(9,8,21,20,7,6,19,18,5,4,17,16,3,2,15,14,13,12,25,24,11,10,23,22)], 2, mean)
## t test and CI of each original priming effect
Orig_DU_effects <- list(
### FA
### SOA 200, High Load
FA_LSOA_HLOAD = round(c( as.numeric( unlist(t.test(orig_DU_comparisons[,9], orig_DU_comparisons[,8],paired = TRUE))[c(1,2,6,4,5,3)]), unlist(d.dep.t.avg(m1 = mean(orig_DU_comparisons[,9]), m2 = mean(orig_DU_comparisons[,8]), sd1 = sd(orig_DU_comparisons[,9]), sd2 = sd(orig_DU_comparisons[,8]), n = length(orig_DU_comparisons[,9]) ))[1:3]), digits = 3),
### SOA 200, Low Load
FA_LSOA_LLOAD = round(c( as.numeric( unlist(t.test(orig_DU_comparisons[,21], orig_DU_comparisons[,20],paired = TRUE))[c(1,2,6,4,5,3)]), unlist(d.dep.t.avg(m1 = mean(orig_DU_comparisons[,21]), m2 = mean(orig_DU_comparisons[,20]), sd1 = sd(orig_DU_comparisons[,21]), sd2 = sd(orig_DU_comparisons[,20]), n = length(orig_DU_comparisons[,21]) ))[1:3]), digits = 3),
### SOA 2000, High Load
FA_HSOA_HLOAD = round(c( as.numeric( unlist(t.test(orig_DU_comparisons[,7], orig_DU_comparisons[,6],paired = TRUE))[c(1,2,6,4,5,3)]), unlist(d.dep.t.avg(m1 = mean(orig_DU_comparisons[,7]), m2 = mean(orig_DU_comparisons[,6]), sd1 = sd(orig_DU_comparisons[,7]), sd2 = sd(orig_DU_comparisons[,6]), n = length(orig_DU_comparisons[,7]) ))[1:3]), digits = 3),
### SOA 2000, Low Load
FA_HSOA_LLOAD = round(c( as.numeric( unlist(t.test(orig_DU_comparisons[,19], orig_DU_comparisons[,18],paired = TRUE))[c(1,2,6,4,5,3)]), unlist(d.dep.t.avg(m1 = mean(orig_DU_comparisons[,19]), m2 = mean(orig_DU_comparisons[,18]), sd1 = sd(orig_DU_comparisons[,19]), sd2 = sd(orig_DU_comparisons[,18]), n = length(orig_DU_comparisons[,19]) ))[1:3]), digits = 3),
### BA
### SOA 200, High Load
BA_LSOA_HLOAD = round(c( as.numeric( unlist(t.test(orig_DU_comparisons[,5], orig_DU_comparisons[,4],paired = TRUE))[c(1,2,6,4,5,3)]), unlist(d.dep.t.avg(m1 = mean(orig_DU_comparisons[,5]), m2 = mean(orig_DU_comparisons[,4]), sd1 = sd(orig_DU_comparisons[,5]), sd2 = sd(orig_DU_comparisons[,4]), n = length(orig_DU_comparisons[,5]) ))[1:3]), digits = 3),
### SOA 200, Low Load
BA_LSOA_LLOAD = round(c( as.numeric( unlist(t.test(orig_DU_comparisons[,17], orig_DU_comparisons[,16],paired = TRUE))[c(1,2,6,4,5,3)]), unlist(d.dep.t.avg(m1 = mean(orig_DU_comparisons[,17]), m2 = mean(orig_DU_comparisons[,16]), sd1 = sd(orig_DU_comparisons[,17]), sd2 = sd(orig_DU_comparisons[,16]), n = length(orig_DU_comparisons[,17]) ))[1:3]), digits = 3),
### SOA 2000, High Load
BA_HSOA_HLOAD = round(c( as.numeric( unlist(t.test(orig_DU_comparisons[,3], orig_DU_comparisons[,2],paired = TRUE))[c(1,2,6,4,5,3)]), unlist(d.dep.t.avg(m1 = mean(orig_DU_comparisons[,3]), m2 = mean(orig_DU_comparisons[,2]), sd1 = sd(orig_DU_comparisons[,3]), sd2 = sd(orig_DU_comparisons[,2]), n = length(orig_DU_comparisons[,3]) ))[1:3]), digits = 3),
### SOA 2000, Low Load
BA_HSOA_LLOAD = round(c( as.numeric( unlist(t.test(orig_DU_comparisons[,15], orig_DU_comparisons[,14],paired = TRUE))[c(1,2,6,4,5,3)]), unlist(d.dep.t.avg(m1 = mean(orig_DU_comparisons[,15]), m2 = mean(orig_DU_comparisons[,14]), sd1 = sd(orig_DU_comparisons[,15]), sd2 = sd(orig_DU_comparisons[,14]), n = length(orig_DU_comparisons[,15]) ))[1:3]), digits = 3),
### SYM
### SOA 200, High Load
SYM_LSOA_HLOAD = round(c( as.numeric( unlist(t.test(orig_DU_comparisons[,13], orig_DU_comparisons[,12],paired = TRUE))[c(1,2,6,4,5,3)]), unlist(d.dep.t.avg(m1 = mean(orig_DU_comparisons[,13]), m2 = mean(orig_DU_comparisons[,12]), sd1 = sd(orig_DU_comparisons[,13]), sd2 = sd(orig_DU_comparisons[,12]), n = length(orig_DU_comparisons[,13]) ))[1:3]), digits = 3),
### SOA 200, Low Load
SYM_LSOA_LLOAD = round(c( as.numeric( unlist(t.test(orig_DU_comparisons[,25], orig_DU_comparisons[,24],paired = TRUE))[c(1,2,6,4,5,3)]), unlist(d.dep.t.avg(m1 = mean(orig_DU_comparisons[,25]), m2 = mean(orig_DU_comparisons[,24]), sd1 = sd(orig_DU_comparisons[,25]), sd2 = sd(orig_DU_comparisons[,24]), n = length(orig_DU_comparisons[,25]) ))[1:3]), digits = 3),
### SOA 2000, High Load
SYM_HSOA_HLOAD = round(c( as.numeric( unlist(t.test(orig_DU_comparisons[,11], orig_DU_comparisons[,10],paired = TRUE))[c(1,2,6,4,5,3)]), unlist(d.dep.t.avg(m1 = mean(orig_DU_comparisons[,11]), m2 = mean(orig_DU_comparisons[,10]), sd1 = sd(orig_DU_comparisons[,11]), sd2 = sd(orig_DU_comparisons[,10]), n = length(orig_DU_comparisons[,11]) ))[1:3]), digits = 3),
### SOA 2000, Low Load
SYM_HSOA_LLOAD = round(c( as.numeric( unlist(t.test(orig_DU_comparisons[,23], orig_DU_comparisons[,22],paired = TRUE))[c(1,2,6,4,5,3)]), unlist(d.dep.t.avg(m1 = mean(orig_DU_comparisons[,23]), m2 = mean(orig_DU_comparisons[,22]), sd1 = sd(orig_DU_comparisons[,23]), sd2 = sd(orig_DU_comparisons[,22]), n = length(orig_DU_comparisons[,23]) ))[1:3]), digits = 3)
)

## Replication study experiment 2, Dutch
Repli_DU_effects <- list(
### FA
### SOA 200, High Load
FA_LSOA_HLOAD = round(c( as.numeric( unlist(t.test(repli_DU_comparisons[,9], repli_DU_comparisons[,8],paired = TRUE))[c(1,2,6,4,5,3)]), unlist(d.dep.t.avg(m1 = mean(repli_DU_comparisons[,9]), m2 = mean(repli_DU_comparisons[,8]), sd1 = sd(repli_DU_comparisons[,9]), sd2 = sd(repli_DU_comparisons[,8]), n = length(repli_DU_comparisons[,9]) ))[1:3]), digits = 3),
### SOA 200, Low Load
FA_LSOA_LLOAD = round(c( as.numeric( unlist(t.test(repli_DU_comparisons[,21], repli_DU_comparisons[,20],paired = TRUE))[c(1,2,6,4,5,3)]), unlist(d.dep.t.avg(m1 = mean(repli_DU_comparisons[,21]), m2 = mean(repli_DU_comparisons[,20]), sd1 = sd(repli_DU_comparisons[,21]), sd2 = sd(repli_DU_comparisons[,20]), n = length(repli_DU_comparisons[,21]) ))[1:3]), digits = 3),
### SOA 2000, High Load
FA_HSOA_HLOAD = round(c( as.numeric( unlist(t.test(repli_DU_comparisons[,7], repli_DU_comparisons[,6],paired = TRUE))[c(1,2,6,4,5,3)]), unlist(d.dep.t.avg(m1 = mean(repli_DU_comparisons[,7]), m2 = mean(repli_DU_comparisons[,6]), sd1 = sd(repli_DU_comparisons[,7]), sd2 = sd(repli_DU_comparisons[,6]), n = length(repli_DU_comparisons[,7]) ))[1:3]), digits = 3),
### SOA 2000, Low Load
FA_HSOA_LLOAD = round(c( as.numeric( unlist(t.test(repli_DU_comparisons[,19], repli_DU_comparisons[,18],paired = TRUE))[c(1,2,6,4,5,3)]), unlist(d.dep.t.avg(m1 = mean(repli_DU_comparisons[,19]), m2 = mean(repli_DU_comparisons[,18]), sd1 = sd(repli_DU_comparisons[,19]), sd2 = sd(repli_DU_comparisons[,18]), n = length(repli_DU_comparisons[,19]) ))[1:3]), digits = 3),
### BA
### SOA 200, High Load
BA_LSOA_HLOAD = round(c( as.numeric( unlist(t.test(repli_DU_comparisons[,5], repli_DU_comparisons[,4],paired = TRUE))[c(1,2,6,4,5,3)]), unlist(d.dep.t.avg(m1 = mean(repli_DU_comparisons[,5]), m2 = mean(repli_DU_comparisons[,4]), sd1 = sd(repli_DU_comparisons[,5]), sd2 = sd(repli_DU_comparisons[,4]), n = length(repli_DU_comparisons[,5]) ))[1:3]), digits = 3),
### SOA 200, Low Load
BA_LSOA_LLOAD = round(c( as.numeric( unlist(t.test(repli_DU_comparisons[,17], repli_DU_comparisons[,16],paired = TRUE))[c(1,2,6,4,5,3)]), unlist(d.dep.t.avg(m1 = mean(repli_DU_comparisons[,17]), m2 = mean(repli_DU_comparisons[,16]), sd1 = sd(repli_DU_comparisons[,17]), sd2 = sd(repli_DU_comparisons[,16]), n = length(repli_DU_comparisons[,17]) ))[1:3]), digits = 3),
### SOA 2000, High Load
BA_HSOA_HLOAD = round(c( as.numeric( unlist(t.test(repli_DU_comparisons[,3], repli_DU_comparisons[,2],paired = TRUE))[c(1,2,6,4,5,3)]), unlist(d.dep.t.avg(m1 = mean(repli_DU_comparisons[,3]), m2 = mean(repli_DU_comparisons[,2]), sd1 = sd(repli_DU_comparisons[,3]), sd2 = sd(repli_DU_comparisons[,2]), n = length(repli_DU_comparisons[,3]) ))[1:3]), digits = 3),
### SOA 2000, Low Load
BA_HSOA_LLOAD = round(c( as.numeric( unlist(t.test(repli_DU_comparisons[,15], repli_DU_comparisons[,14],paired = TRUE))[c(1,2,6,4,5,3)]), unlist(d.dep.t.avg(m1 = mean(repli_DU_comparisons[,15]), m2 = mean(repli_DU_comparisons[,14]), sd1 = sd(repli_DU_comparisons[,15]), sd2 = sd(repli_DU_comparisons[,14]), n = length(repli_DU_comparisons[,15]) ))[1:3]), digits = 3),
### SYM
### SOA 200, High Load
SYM_LSOA_HLOAD = round(c( as.numeric( unlist(t.test(repli_DU_comparisons[,13], repli_DU_comparisons[,12],paired = TRUE))[c(1,2,6,4,5,3)]), unlist(d.dep.t.avg(m1 = mean(repli_DU_comparisons[,13]), m2 = mean(repli_DU_comparisons[,12]), sd1 = sd(repli_DU_comparisons[,13]), sd2 = sd(repli_DU_comparisons[,12]), n = length(repli_DU_comparisons[,13]) ))[1:3]), digits = 3),
### SOA 200, Low Load
SYM_LSOA_LLOAD = round(c( as.numeric( unlist(t.test(repli_DU_comparisons[,25], repli_DU_comparisons[,24],paired = TRUE))[c(1,2,6,4,5,3)]), unlist(d.dep.t.avg(m1 = mean(repli_DU_comparisons[,25]), m2 = mean(repli_DU_comparisons[,24]), sd1 = sd(repli_DU_comparisons[,25]), sd2 = sd(repli_DU_comparisons[,24]), n = length(repli_DU_comparisons[,25]) ))[1:3]), digits = 3),
### SOA 2000, High Load
SYM_HSOA_HLOAD = round(c( as.numeric( unlist(t.test(repli_DU_comparisons[,11], repli_DU_comparisons[,10],paired = TRUE))[c(1,2,6,4,5,3)]), unlist(d.dep.t.avg(m1 = mean(repli_DU_comparisons[,11]), m2 = mean(repli_DU_comparisons[,10]), sd1 = sd(repli_DU_comparisons[,11]), sd2 = sd(repli_DU_comparisons[,10]), n = length(repli_DU_comparisons[,11]) ))[1:3]), digits = 3),
### SOA 2000, Low Load
SYM_HSOA_LLOAD = round(c( as.numeric( unlist(t.test(repli_DU_comparisons[,23], repli_DU_comparisons[,22],paired = TRUE))[c(1,2,6,4,5,3)]), unlist(d.dep.t.avg(m1 = mean(repli_DU_comparisons[,23]), m2 = mean(repli_DU_comparisons[,22]), sd1 = sd(repli_DU_comparisons[,23]), sd2 = sd(repli_DU_comparisons[,22]), n = length(repli_DU_comparisons[,23]) ))[1:3]), digits = 3)
)

## Replication Study, Experiment 1 English
Repli_EN_effects <- list(
### FA
### SOA 200, High Load
FA_LSOA_HLOAD = round(c( as.numeric( unlist(t.test(repli_EN_comparisons[,9], repli_EN_comparisons[,8],paired = TRUE))[c(1,2,6,4,5,3)]), unlist(d.dep.t.avg(m1 = mean(repli_EN_comparisons[,9]), m2 = mean(repli_EN_comparisons[,8]), sd1 = sd(repli_EN_comparisons[,9]), sd2 = sd(repli_EN_comparisons[,8]), n = length(repli_EN_comparisons[,9]) ))[1:3]), digits = 3),
### SOA 200, Low Load
FA_LSOA_LLOAD = round(c( as.numeric( unlist(t.test(repli_EN_comparisons[,21], repli_EN_comparisons[,20],paired = TRUE))[c(1,2,6,4,5,3)]), unlist(d.dep.t.avg(m1 = mean(repli_EN_comparisons[,21]), m2 = mean(repli_EN_comparisons[,20]), sd1 = sd(repli_EN_comparisons[,21]), sd2 = sd(repli_EN_comparisons[,20]), n = length(repli_EN_comparisons[,21]) ))[1:3]), digits = 3),
### SOA 2000, High Load
FA_HSOA_HLOAD = round(c( as.numeric( unlist(t.test(repli_EN_comparisons[,7], repli_EN_comparisons[,6],paired = TRUE))[c(1,2,6,4,5,3)]), unlist(d.dep.t.avg(m1 = mean(repli_EN_comparisons[,7]), m2 = mean(repli_EN_comparisons[,6]), sd1 = sd(repli_EN_comparisons[,7]), sd2 = sd(repli_EN_comparisons[,6]), n = length(repli_EN_comparisons[,7]) ))[1:3]), digits = 3),
### SOA 2000, Low Load
FA_HSOA_LLOAD = round(c( as.numeric( unlist(t.test(repli_EN_comparisons[,19], repli_EN_comparisons[,18],paired = TRUE))[c(1,2,6,4,5,3)]), unlist(d.dep.t.avg(m1 = mean(repli_EN_comparisons[,19]), m2 = mean(repli_EN_comparisons[,18]), sd1 = sd(repli_EN_comparisons[,19]), sd2 = sd(repli_EN_comparisons[,18]), n = length(repli_EN_comparisons[,19]) ))[1:3]), digits = 3),
### BA
### SOA 200, High Load
BA_LSOA_HLOAD = round(c( as.numeric( unlist(t.test(repli_EN_comparisons[,5], repli_EN_comparisons[,4],paired = TRUE))[c(1,2,6,4,5,3)]), unlist(d.dep.t.avg(m1 = mean(repli_EN_comparisons[,5]), m2 = mean(repli_EN_comparisons[,4]), sd1 = sd(repli_EN_comparisons[,5]), sd2 = sd(repli_EN_comparisons[,4]), n = length(repli_EN_comparisons[,5]) ))[1:3]), digits = 3),
### SOA 200, Low Load
BA_LSOA_LLOAD = round(c( as.numeric( unlist(t.test(repli_EN_comparisons[,17], repli_EN_comparisons[,16],paired = TRUE))[c(1,2,6,4,5,3)]), unlist(d.dep.t.avg(m1 = mean(repli_EN_comparisons[,17]), m2 = mean(repli_EN_comparisons[,16]), sd1 = sd(repli_EN_comparisons[,17]), sd2 = sd(repli_EN_comparisons[,16]), n = length(repli_EN_comparisons[,17]) ))[1:3]), digits = 3),
### SOA 2000, High Load
BA_HSOA_HLOAD = round(c( as.numeric( unlist(t.test(repli_EN_comparisons[,3], repli_EN_comparisons[,2],paired = TRUE))[c(1,2,6,4,5,3)]), unlist(d.dep.t.avg(m1 = mean(repli_EN_comparisons[,3]), m2 = mean(repli_EN_comparisons[,2]), sd1 = sd(repli_EN_comparisons[,3]), sd2 = sd(repli_EN_comparisons[,2]), n = length(repli_EN_comparisons[,3]) ))[1:3]), digits = 3),
### SOA 2000, Low Load
BA_HSOA_LLOAD = round(c( as.numeric( unlist(t.test(repli_EN_comparisons[,15], repli_EN_comparisons[,14],paired = TRUE))[c(1,2,6,4,5,3)]), unlist(d.dep.t.avg(m1 = mean(repli_EN_comparisons[,15]), m2 = mean(repli_EN_comparisons[,14]), sd1 = sd(repli_EN_comparisons[,15]), sd2 = sd(repli_EN_comparisons[,14]), n = length(repli_EN_comparisons[,15]) ))[1:3]), digits = 3),
### SYM
### SOA 200, High Load
SYM_LSOA_HLOAD = round(c( as.numeric( unlist(t.test(repli_EN_comparisons[,13], repli_EN_comparisons[,12],paired = TRUE))[c(1,2,6,4,5,3)]), unlist(d.dep.t.avg(m1 = mean(repli_EN_comparisons[,13]), m2 = mean(repli_EN_comparisons[,12]), sd1 = sd(repli_EN_comparisons[,13]), sd2 = sd(repli_EN_comparisons[,12]), n = length(repli_EN_comparisons[,13]) ))[1:3]), digits = 3),
### SOA 200, Low Load
SYM_LSOA_LLOAD = round(c( as.numeric( unlist(t.test(repli_EN_comparisons[,25], repli_EN_comparisons[,24],paired = TRUE))[c(1,2,6,4,5,3)]), unlist(d.dep.t.avg(m1 = mean(repli_EN_comparisons[,25]), m2 = mean(repli_EN_comparisons[,24]), sd1 = sd(repli_EN_comparisons[,25]), sd2 = sd(repli_EN_comparisons[,24]), n = length(repli_EN_comparisons[,25]) ))[1:3]), digits = 3),
### SOA 2000, High Load
SYM_HSOA_HLOAD = round(c( as.numeric( unlist(t.test(repli_EN_comparisons[,11], repli_EN_comparisons[,10],paired = TRUE))[c(1,2,6,4,5,3)]), unlist(d.dep.t.avg(m1 = mean(repli_EN_comparisons[,11]), m2 = mean(repli_EN_comparisons[,10]), sd1 = sd(repli_EN_comparisons[,11]), sd2 = sd(repli_EN_comparisons[,10]), n = length(repli_EN_comparisons[,11]) ))[1:3]), digits = 3),
### SOA 2000, Low Load
SYM_HSOA_LLOAD = round(c( as.numeric( unlist(t.test(repli_EN_comparisons[,23], repli_EN_comparisons[,22],paired = TRUE))[c(1,2,6,4,5,3)]), unlist(d.dep.t.avg(m1 = mean(repli_EN_comparisons[,23]), m2 = mean(repli_EN_comparisons[,22]), sd1 = sd(repli_EN_comparisons[,23]), sd2 = sd(repli_EN_comparisons[,22]), n = length(repli_EN_comparisons[,23]) ))[1:3]), digits = 3)
)

## TOST for every priming effect of the replication study
## Set the alpha level for every test. I'm thinking a lower alpha level would be better.
alpha_level = .05

## The bounds for every TOST refers to the effect size of original priming effects
## TOST decision function
TOST_decision <- function(TOST_output){
M = TOST_output["M"]
LCI = TOST_output["LL_CI_TOST"]
UCI = TOST_output["UL_CI_TOST"] 
low_bound = TOST_output["low_eqbound.dlow"]
high_bound = TOST_output["high_eqbound.dhigh"]
if(LCI >= low_bound & UCI <= high_bound){
return("Equivalent")
}
else if(LCI <= low_bound & UCI >= high_bound){
return("Inconclusive")
}
else{
if(M > high_bound){
return("Superior")
}
else if(M < low_bound){
return("Inferior")
}else{
return("Non-inferior")
}
}
}

#################################
repli_EN_TOST <- list(
## Replication Experiment 1, FA pairs, SOA = 200 ms, High Load
repli_EN_FA_LSOA_HLOAD_TOST = c(M = mean(repli_EN_comparisons[,9]) - mean(repli_EN_comparisons[,8]) ,unlist(TOSTpaired(n = length(repli_EN_comparisons[,9]), m1 = mean(repli_EN_comparisons[,9]), m2 = mean(repli_EN_comparisons[,8]), sd1 = sd(repli_EN_comparisons[,9]), sd2 = sd(repli_EN_comparisons[,8]), r12 = cor(repli_EN_comparisons[,9], repli_EN_comparisons[,8]), low_eqbound_dz = Orig_DU_effects$FA_LSOA_HLOAD["dlow"], high_eqbound_dz = Orig_DU_effects$FA_LSOA_HLOAD["dhigh"], alpha = alpha_level))),

## Replication Experiment 1, FA pairs, SOA = 200 ms, Low Load
repli_EN_FA_LSOA_LLOAD_TOST = c(M = mean(repli_EN_comparisons[,21]) - mean(repli_EN_comparisons[,20]) ,unlist(TOSTpaired(n = length(repli_EN_comparisons[,21]), m1 = mean(repli_EN_comparisons[,21]), m2 = mean(repli_EN_comparisons[,20]), sd1 = sd(repli_EN_comparisons[,21]), sd2 = sd(repli_EN_comparisons[,20]), r12 = cor(repli_EN_comparisons[,21], repli_EN_comparisons[,20]), low_eqbound_dz = Orig_DU_effects$FA_LSOA_LLOAD["dlow"], high_eqbound_dz = Orig_DU_effects$FA_LSOA_LLOAD["dhigh"], alpha = alpha_level))),

## Replication Experiment 1, FA pairs, SOA = 2000 ms, High Load
repli_EN_FA_HSOA_HLOAD_TOST = c(M = mean(repli_EN_comparisons[,7]) - mean(repli_EN_comparisons[,6]) ,unlist(TOSTpaired(n = length(repli_EN_comparisons[,7]), m1 = mean(repli_EN_comparisons[,7]), m2 = mean(repli_EN_comparisons[,6]), sd1 = sd(repli_EN_comparisons[,7]), sd2 = sd(repli_EN_comparisons[,6]), r12 = cor(repli_EN_comparisons[,7], repli_EN_comparisons[,6]), low_eqbound_dz = Orig_DU_effects$FA_HSOA_HLOAD["dlow"], high_eqbound_dz = Orig_DU_effects$FA_HSOA_HLOAD["dhigh"], alpha = alpha_level))),

## Replication Experiment 1, FA pairs, SOA = 2000 ms, Low Load
repli_EN_FA_HSOA_LLOAD_TOST = c(M = mean(repli_EN_comparisons[,19]) - mean(repli_EN_comparisons[,18]) ,unlist(TOSTpaired(n = length(repli_EN_comparisons[,19]), m1 = mean(repli_EN_comparisons[,19]), m2 = mean(repli_EN_comparisons[,18]), sd1 = sd(repli_EN_comparisons[,19]), sd2 = sd(repli_EN_comparisons[,18]), r12 = cor(repli_EN_comparisons[,19], repli_EN_comparisons[,18]), low_eqbound_dz = Orig_DU_effects$FA_HSOA_LLOAD["dlow"], high_eqbound_dz = Orig_DU_effects$FA_HSOA_LLOAD["dhigh"], alpha = alpha_level))),

## Replication Experiment 1, BA pairs, SOA = 200 ms, High Load
repli_EN_BA_LSOA_HLOAD_TOST = c(M = mean(repli_EN_comparisons[,5]) - mean(repli_EN_comparisons[,4]) ,unlist(TOSTpaired(n = length(repli_EN_comparisons[,5]), m1 = mean(repli_EN_comparisons[,5]), m2 = mean(repli_EN_comparisons[,4]), sd1 = sd(repli_EN_comparisons[,5]), sd2 = sd(repli_EN_comparisons[,4]), r12 = cor(repli_EN_comparisons[,5], repli_EN_comparisons[,4]), low_eqbound_dz = Orig_DU_effects$BA_LSOA_HLOAD["dlow"], high_eqbound_dz = Orig_DU_effects$BA_LSOA_HLOAD["dhigh"], alpha = alpha_level))),

## Replication Experiment 1, BA pairs, SOA = 200 ms, Low Load
repli_EN_BA_LSOA_LLOAD_TOST = c(M = mean(repli_EN_comparisons[,17]) - mean(repli_EN_comparisons[,16]) ,unlist(TOSTpaired(n = length(repli_EN_comparisons[,17]), m1 = mean(repli_EN_comparisons[,17]), m2 = mean(repli_EN_comparisons[,16]), sd1 = sd(repli_EN_comparisons[,17]), sd2 = sd(repli_EN_comparisons[,16]), r12 = cor(repli_EN_comparisons[,17], repli_EN_comparisons[,16]), low_eqbound_dz = Orig_DU_effects$BA_LSOA_LLOAD["dlow"], high_eqbound_dz = Orig_DU_effects$BA_LSOA_LLOAD["dhigh"], alpha = alpha_level))),

## Replication Experiment 1, BA pairs, SOA = 2000 ms, High Load
repli_EN_BA_HSOA_HLOAD_TOST = c(M = mean(repli_EN_comparisons[,3]) - mean(repli_EN_comparisons[,2]), unlist(TOSTpaired(n = length(repli_EN_comparisons[,3]), m1 = mean(repli_EN_comparisons[,3]), m2 = mean(repli_EN_comparisons[,2]), sd1 = sd(repli_EN_comparisons[,3]), sd2 = sd(repli_EN_comparisons[,2]), r12 = cor(repli_EN_comparisons[,3], repli_EN_comparisons[,2]), low_eqbound_dz = Orig_DU_effects$BA_HSOA_HLOAD["dlow"], high_eqbound_dz = Orig_DU_effects$BA_HSOA_HLOAD["dhigh"], alpha = alpha_level))),

## Replication Experiment 1, BA pairs, SOA = 2000 ms, Low Load
repli_EN_BA_HSOA_LLOAD_TOST = c(M = mean(repli_EN_comparisons[,15]) - mean(repli_EN_comparisons[,14]), unlist(TOSTpaired(n = length(repli_EN_comparisons[,15]), m1 = mean(repli_EN_comparisons[,15]), m2 = mean(repli_EN_comparisons[,14]), sd1 = sd(repli_EN_comparisons[,15]), sd2 = sd(repli_EN_comparisons[,14]), r12 = cor(repli_EN_comparisons[,15], repli_EN_comparisons[,14]), low_eqbound_dz = Orig_DU_effects$BA_HSOA_LLOAD["dlow"], high_eqbound_dz = Orig_DU_effects$BA_HSOA_LLOAD["dhigh"], alpha = alpha_level))),

## Replication Experiment 1, SYM pairs, SOA = 200 ms, High Load
repli_EN_SYM_LSOA_HLOAD_TOST = c(M = mean(repli_EN_comparisons[,13]) - mean(repli_EN_comparisons[,12]), unlist(TOSTpaired(n = length(repli_EN_comparisons[,13]), m1 = mean(repli_EN_comparisons[,13]), m2 = mean(repli_EN_comparisons[,12]), sd1 = sd(repli_EN_comparisons[,13]), sd2 = sd(repli_EN_comparisons[,12]), r12 = cor(repli_EN_comparisons[,13], repli_EN_comparisons[,12]), low_eqbound_dz = Orig_DU_effects$SYM_LSOA_HLOAD["dlow"], high_eqbound_dz = Orig_DU_effects$SYM_LSOA_HLOAD["dhigh"], alpha = alpha_level))),

## Replication Experiment 1, SYM pairs, SOA = 200 ms, Low Load
repli_EN_SYM_LSOA_LLOAD_TOST = c(M = mean(repli_EN_comparisons[,25]) - mean(repli_EN_comparisons[,24]), unlist(TOSTpaired(n = length(repli_EN_comparisons[,25]), m1 = mean(repli_EN_comparisons[,25]), m2 = mean(repli_EN_comparisons[,24]), sd1 = sd(repli_EN_comparisons[,25]), sd2 = sd(repli_EN_comparisons[,24]), r12 = cor(repli_EN_comparisons[,25], repli_EN_comparisons[,24]), low_eqbound_dz = Orig_DU_effects$SYM_LSOA_LLOAD["dlow"], high_eqbound_dz = Orig_DU_effects$SYM_LSOA_LLOAD["dhigh"], alpha = alpha_level))),

## Replication Experiment 1, SYM pairs, SOA = 2000 ms, High Load
repli_EN_SYM_HSOA_HLOAD_TOST = c(M = mean(repli_EN_comparisons[,11]) - mean(repli_EN_comparisons[,10]), unlist(TOSTpaired(n = length(repli_EN_comparisons[,11]), m1 = mean(repli_EN_comparisons[,11]), m2 = mean(repli_EN_comparisons[,10]), sd1 = sd(repli_EN_comparisons[,11]), sd2 = sd(repli_EN_comparisons[,10]), r12 = cor(repli_EN_comparisons[,11], repli_EN_comparisons[,10]), low_eqbound_dz = Orig_DU_effects$SYM_HSOA_HLOAD["dlow"], high_eqbound_dz = Orig_DU_effects$SYM_HSOA_HLOAD["dhigh"], alpha = alpha_level))),

## Replication Experiment 1, SYM pairs, SOA = 2000 ms, Low Load
repli_EN_SYM_HSOA_LLOAD_TOST = c(M = mean(repli_EN_comparisons[,23]) - mean(repli_EN_comparisons[,22]), unlist(TOSTpaired(n = length(repli_EN_comparisons[,23]), m1 = mean(repli_EN_comparisons[,23]), m2 = mean(repli_EN_comparisons[,22]), sd1 = sd(repli_EN_comparisons[,23]), sd2 = sd(repli_EN_comparisons[,22]), r12 = cor(repli_EN_comparisons[,23], repli_EN_comparisons[,22]), low_eqbound_dz = Orig_DU_effects$SYM_HSOA_LLOAD["dlow"], high_eqbound_dz = Orig_DU_effects$SYM_HSOA_LLOAD["dhigh"], alpha = alpha_level)))
)

## Check the conclusion of each TOST
repli_EN_TOST_results <- lapply(repli_EN_TOST, TOST_decision)

#################################
repli_DU_TOST <- list(
## Replication Experiment 1, FA pairs, SOA = 200 ms, High Load
repli_DU_FA_LSOA_HLOAD_TOST = c(M = mean(repli_DU_comparisons[,9]) - mean(repli_DU_comparisons[,8]) ,unlist(TOSTpaired(n = length(repli_DU_comparisons[,9]), m1 = mean(repli_DU_comparisons[,9]), m2 = mean(repli_DU_comparisons[,8]), sd1 = sd(repli_DU_comparisons[,9]), sd2 = sd(repli_DU_comparisons[,8]), r12 = cor(repli_DU_comparisons[,9], repli_DU_comparisons[,8]), low_eqbound_dz = Orig_DU_effects$FA_LSOA_HLOAD["dlow"], high_eqbound_dz = Orig_DU_effects$FA_LSOA_HLOAD["dhigh"], alpha = alpha_level))),

## Replication Experiment 1, FA pairs, SOA = 200 ms, Low Load
repli_DU_FA_LSOA_LLOAD_TOST = c(M = mean(repli_DU_comparisons[,21]) - mean(repli_DU_comparisons[,20]) ,unlist(TOSTpaired(n = length(repli_DU_comparisons[,21]), m1 = mean(repli_DU_comparisons[,21]), m2 = mean(repli_DU_comparisons[,20]), sd1 = sd(repli_DU_comparisons[,21]), sd2 = sd(repli_DU_comparisons[,20]), r12 = cor(repli_DU_comparisons[,21], repli_DU_comparisons[,20]), low_eqbound_dz = Orig_DU_effects$FA_LSOA_LLOAD["dlow"], high_eqbound_dz = Orig_DU_effects$FA_LSOA_LLOAD["dhigh"], alpha = alpha_level))),

## Replication Experiment 1, FA pairs, SOA = 2000 ms, High Load
repli_DU_FA_HSOA_HLOAD_TOST = c(M = mean(repli_DU_comparisons[,7]) - mean(repli_DU_comparisons[,6]) ,unlist(TOSTpaired(n = length(repli_DU_comparisons[,7]), m1 = mean(repli_DU_comparisons[,7]), m2 = mean(repli_DU_comparisons[,6]), sd1 = sd(repli_DU_comparisons[,7]), sd2 = sd(repli_DU_comparisons[,6]), r12 = cor(repli_DU_comparisons[,7], repli_DU_comparisons[,6]), low_eqbound_dz = Orig_DU_effects$FA_HSOA_HLOAD["dlow"], high_eqbound_dz = Orig_DU_effects$FA_HSOA_HLOAD["dhigh"], alpha = alpha_level))),

## Replication Experiment 1, FA pairs, SOA = 2000 ms, Low Load
repli_DU_FA_HSOA_LLOAD_TOST = c(M = mean(repli_DU_comparisons[,19]) - mean(repli_DU_comparisons[,18]) ,unlist(TOSTpaired(n = length(repli_DU_comparisons[,19]), m1 = mean(repli_DU_comparisons[,19]), m2 = mean(repli_DU_comparisons[,18]), sd1 = sd(repli_DU_comparisons[,19]), sd2 = sd(repli_DU_comparisons[,18]), r12 = cor(repli_DU_comparisons[,19], repli_DU_comparisons[,18]), low_eqbound_dz = Orig_DU_effects$FA_HSOA_LLOAD["dlow"], high_eqbound_dz = Orig_DU_effects$FA_HSOA_LLOAD["dhigh"], alpha = alpha_level))),

## Replication Experiment 1, BA pairs, SOA = 200 ms, High Load
repli_DU_BA_LSOA_HLOAD_TOST = c(M = mean(repli_DU_comparisons[,5]) - mean(repli_DU_comparisons[,4]) ,unlist(TOSTpaired(n = length(repli_DU_comparisons[,5]), m1 = mean(repli_DU_comparisons[,5]), m2 = mean(repli_DU_comparisons[,4]), sd1 = sd(repli_DU_comparisons[,5]), sd2 = sd(repli_DU_comparisons[,4]), r12 = cor(repli_DU_comparisons[,5], repli_DU_comparisons[,4]), low_eqbound_dz = Orig_DU_effects$BA_LSOA_HLOAD["dlow"], high_eqbound_dz = Orig_DU_effects$BA_LSOA_HLOAD["dhigh"], alpha = alpha_level))),

## Replication Experiment 1, BA pairs, SOA = 200 ms, Low Load
repli_DU_BA_LSOA_LLOAD_TOST = c(M = mean(repli_DU_comparisons[,17]) - mean(repli_DU_comparisons[,16]) ,unlist(TOSTpaired(n = length(repli_DU_comparisons[,17]), m1 = mean(repli_DU_comparisons[,17]), m2 = mean(repli_DU_comparisons[,16]), sd1 = sd(repli_DU_comparisons[,17]), sd2 = sd(repli_DU_comparisons[,16]), r12 = cor(repli_DU_comparisons[,17], repli_DU_comparisons[,16]), low_eqbound_dz = Orig_DU_effects$BA_LSOA_LLOAD["dlow"], high_eqbound_dz = Orig_DU_effects$BA_LSOA_LLOAD["dhigh"], alpha = alpha_level))),

## Replication Experiment 1, BA pairs, SOA = 2000 ms, High Load
repli_DU_BA_HSOA_HLOAD_TOST = c(M = mean(repli_DU_comparisons[,3]) - mean(repli_DU_comparisons[,2]), unlist(TOSTpaired(n = length(repli_DU_comparisons[,3]), m1 = mean(repli_DU_comparisons[,3]), m2 = mean(repli_DU_comparisons[,2]), sd1 = sd(repli_DU_comparisons[,3]), sd2 = sd(repli_DU_comparisons[,2]), r12 = cor(repli_DU_comparisons[,3], repli_DU_comparisons[,2]), low_eqbound_dz = Orig_DU_effects$BA_HSOA_HLOAD["dlow"], high_eqbound_dz = Orig_DU_effects$BA_HSOA_HLOAD["dhigh"], alpha = alpha_level))),

## Replication Experiment 1, BA pairs, SOA = 2000 ms, Low Load
repli_DU_BA_HSOA_LLOAD_TOST = c(M = mean(repli_DU_comparisons[,15]) - mean(repli_DU_comparisons[,14]), unlist(TOSTpaired(n = length(repli_DU_comparisons[,15]), m1 = mean(repli_DU_comparisons[,15]), m2 = mean(repli_DU_comparisons[,14]), sd1 = sd(repli_DU_comparisons[,15]), sd2 = sd(repli_DU_comparisons[,14]), r12 = cor(repli_DU_comparisons[,15], repli_DU_comparisons[,14]), low_eqbound_dz = Orig_DU_effects$BA_HSOA_LLOAD["dlow"], high_eqbound_dz = Orig_DU_effects$BA_HSOA_LLOAD["dhigh"], alpha = alpha_level))),

## Replication Experiment 1, SYM pairs, SOA = 200 ms, High Load
repli_DU_SYM_LSOA_HLOAD_TOST = c(M = mean(repli_DU_comparisons[,13]) - mean(repli_DU_comparisons[,12]), unlist(TOSTpaired(n = length(repli_DU_comparisons[,13]), m1 = mean(repli_DU_comparisons[,13]), m2 = mean(repli_DU_comparisons[,12]), sd1 = sd(repli_DU_comparisons[,13]), sd2 = sd(repli_DU_comparisons[,12]), r12 = cor(repli_DU_comparisons[,13], repli_DU_comparisons[,12]), low_eqbound_dz = Orig_DU_effects$SYM_LSOA_HLOAD["dlow"], high_eqbound_dz = Orig_DU_effects$SYM_LSOA_HLOAD["dhigh"], alpha = alpha_level))),

## Replication Experiment 1, SYM pairs, SOA = 200 ms, Low Load
repli_DU_SYM_LSOA_LLOAD_TOST = c(M = mean(repli_DU_comparisons[,25]) - mean(repli_DU_comparisons[,24]), unlist(TOSTpaired(n = length(repli_DU_comparisons[,25]), m1 = mean(repli_DU_comparisons[,25]), m2 = mean(repli_DU_comparisons[,24]), sd1 = sd(repli_DU_comparisons[,25]), sd2 = sd(repli_DU_comparisons[,24]), r12 = cor(repli_DU_comparisons[,25], repli_DU_comparisons[,24]), low_eqbound_dz = Orig_DU_effects$SYM_LSOA_LLOAD["dlow"], high_eqbound_dz = Orig_DU_effects$SYM_LSOA_LLOAD["dhigh"], alpha = alpha_level))),

## Replication Experiment 1, SYM pairs, SOA = 2000 ms, High Load
repli_DU_SYM_HSOA_HLOAD_TOST = c(M = mean(repli_DU_comparisons[,11]) - mean(repli_DU_comparisons[,10]), unlist(TOSTpaired(n = length(repli_DU_comparisons[,11]), m1 = mean(repli_DU_comparisons[,11]), m2 = mean(repli_DU_comparisons[,10]), sd1 = sd(repli_DU_comparisons[,11]), sd2 = sd(repli_DU_comparisons[,10]), r12 = cor(repli_DU_comparisons[,11], repli_DU_comparisons[,10]), low_eqbound_dz = Orig_DU_effects$SYM_HSOA_HLOAD["dlow"], high_eqbound_dz = Orig_DU_effects$SYM_HSOA_HLOAD["dhigh"], alpha = alpha_level))),

## Replication Experiment 1, SYM pairs, SOA = 2000 ms, Low Load
repli_DU_SYM_HSOA_LLOAD_TOST = c(M = mean(repli_DU_comparisons[,23]) - mean(repli_DU_comparisons[,22]), unlist(TOSTpaired(n = length(repli_DU_comparisons[,23]), m1 = mean(repli_DU_comparisons[,23]), m2 = mean(repli_DU_comparisons[,22]), sd1 = sd(repli_DU_comparisons[,23]), sd2 = sd(repli_DU_comparisons[,22]), r12 = cor(repli_DU_comparisons[,23], repli_DU_comparisons[,22]), low_eqbound_dz = Orig_DU_effects$SYM_HSOA_LLOAD["dlow"], high_eqbound_dz = Orig_DU_effects$SYM_HSOA_LLOAD["dhigh"], alpha = alpha_level)))
)

## Check the conclusion of each TOST
repli_DU_TOST_results <- lapply(repli_DU_TOST, TOST_decision)

## Save the objects for 
save.image('WM_TOST.RData')
