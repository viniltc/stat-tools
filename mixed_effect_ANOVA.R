

## ANOVA - Vinil
## Data must be organized in this form
## DV--IV(Group)--IV(Exposure)--Subj


#getwd() 
#setwd("C:/Users/Asus/Documents/Data")
#list.files()

# read files 
ownvia<-read.table("ownvia_fourcolumn.txt",header=TRUE)

parvia<-read.table("parvia_fourcolumn.txt")
catcht <-read.table("catch_fourcolumn.txt")
intforce<-read.table("intforce_epoch.txt")
intforce1<-read.table("intforce_epoch1.txt")
mindist1<-read.table("mindist_epoch1.txt")

mdfit<-read.table("md_fitting_coeffs.txt")
deltafit<-read.table("deltatc_fitting_coeffs.txt")
#md_sub1<-read.table("md12_fourcolumn.txt")
#md_sub2<-read.table("md21_fourcolumn.txt")
md_sub1<-read.table("md12_three_group_fourcolumn.txt") # three group(H VH PV) , three time (pre mid post) 
md_sub2<-read.table("md21_three_group_fourcolumn.txt")
#score1<-read.table("score_1_fourcolumn.txt")
#score2<-read.table("score_2_fourcolumn.txt")
score1<-read.table("score_1_three_group_fourcolumn.txt")
score2<-read.table("score_2_three_group_fourcolumn.txt")
#intforce<-read.table("intf_fourcolumn.txt")
intforce<-read.table("intf_three_group_fourcolumn.txt")

LIndex_vp<-read.table("lead_vp.txt") # rows 1:10 H, 11:20 VH


lead11<-read.table("lead11_anova.txt")
lead12<-read.table("lead12_anova.txt")
lead21<-read.table("lead21_anova.txt")
lead22<-read.table("lead22_anova.txt")

lead11_late<-lead11[lead11[, "V3"] == 3,]
lead12_late<-lead12[lead12[, "V3"] == 3,]
lead22_late<-lead22[lead22[, "V3"] == 3,]
lead21_late<-lead21[lead21[, "V3"] == 3,]


#load library plyr and nlme
library(plyr)
library(nlme)

# Repeated mixed effect ANOVA using 'lme'

# step 1 : Grouping data

ownvia<-read.table("ownvia_fourcolumn.txt",header=TRUE)
ownvia <- within(ownvia, {
  subject <- factor(subject)
  group <- factor(group)
  #time <- factor(time)

})
ownd<- groupedData(MD~group|subject, data = ownvia)
#attach(ownd)
# step 2 : perform ANOVA using lme
ownd.mod<-lme(MD~group*time, random =~1|subject, data = ownd)
ownd.anova<-anova(ownd.mod)
ownd.anova
# to save ANOVA table as doc file!
capture.output(ownd.anova, file = "ownvia_anova.doc")
#ttest (planned comparisons)
t.test(V1[time==2 & group ==1],V1[time==2 & group==2])  #this is significant
t.test(V1[time==1 & group ==1],V1[time==1 & group==2])  #this is NOT significant
t.test(V1[time==3 & group ==1],V1[time==3 & group==2])  #this is NOT significant


#same steps for other indicators!

# partner via point distance
parvia <- within(parvia, {
  subject <- factor(V4)
  group <- factor(3-V2)
  time <- factor(V3)
})
pard<- groupedData(V1~group|subject, data = parvia)
attach(pard)
pard.mod<-lme(V1~group*time, random =~1|subject, data = pard)
pard.anova<-anova(pard.mod)
pard.anova
capture.output(pard.anova, file = "parvia_anova.doc")

#ttest (planned comparisons)
t.test(V1[time==2 & group ==1],V1[time==2 & group==2])  #this is significant
t.test(V1[time==1 & group ==1],V1[time==1 & group==2])  #this is NOT significant
t.test(V1[time==3 & group ==1],V1[time==3 & group==2])  #this is NOT significant


# Subject 1 , MD
md_sub1 <- within(md_sub1, {
  subject <- factor(V4)
  group <- factor(V2)
  time <- factor(V3)
})
sub1<- groupedData(V1~group|subject, data = md_sub1)
attach(sub1)
sub1.mod<-lme(V1~group*time, random =~1|subject, data = sub1)
sub1.anova<-anova(sub1.mod)
sub1.anova
#capture.output(sub1.anova, file = "md12_anova.doc")
capture.output(sub1.anova, file = "md12_three_group_anova.doc")

#ttest (planned comparisons)
t1<-t.test(V1[time==2 & group ==1],V1[time==2 & group==2])  #this is significant
t2<-t.test(V1[time==1 & group ==1],V1[time==1 & group==2])  #this is NOT significant
t3<-t.test(V1[time==3 & group ==1],V1[time==3 & group==2])  #this is NOT significant
capture.output(pard.anova, file = "md12_PC.doc")

# steps to perform posthoc analysis using Tuckey's procedure

sub1_t1<-subset(md_sub1, md_sub1$time==1)
sub1_t2<-subset(md_sub1, md_sub1$time==2) # create subset of 3 levels of time
sub1_t3<-subset(md_sub1, md_sub1$time==3)

# Subject 2 , MD
aov_sub1_t1<-aov(sub1_t1$V1~sub1_t1$group)
aov_sub1_t2<-aov(sub1_t2$V1~sub1_t2$group) # perform the one-way anova for above subsets
aov_sub1_t3<-aov(sub1_t3$V1~sub1_t3$group)

summary(aov_sub1_t1)
summary(aov_sub1_t2)  # summary anova table..
summary(aov_sub1_t3)

tu1<-TukeyHSD(aov_sub1_t1)
tu2<-TukeyHSD(aov_sub1_t2) # finally pairwise comparison for the simple effect using Tukey
tu3<-TukeyHSD(aov_sub1_t3)
capture.output(tu1, file = "md12_tuckey1.doc")
capture.output(tu2, file = "md12_tuckey2.doc")
capture.output(tu3, file = "md12_tuckey3.doc")

# To explore interaction effect with the help of "simple effects" analysis
#With an interaction effect, 
#Conduct a simple effects analysis of the variable 'group' 
#on the outcome variable errors at each level of 'time'

summary(aov_sub1_t1)
summary(aov_sub1_t2)  # to check if simple effects are signficant!!
summary(aov_sub1_t3)

#The definition of an interaction effect states that -
#the effect of one variable changes across levels of the other variable. -
#For example, we might expect the effect of 'group'(sensory information) to be greater in the late 'time' level than in the early time!.

#In order to really understand the different effect sizes, you should make use of the etaSquare

library(lsr) # for etaSquared() funtion

etaSquared(aov_sub1_t1, anova= TRUE)
etaSquared(aov_sub1_t2, anova= TRUE)
etaSquared(aov_sub1_t3, anova= TRUE)



md_sub2 <- within(md_sub2, {
  subject <- factor(V4)
  group <- factor(V2)
  time <- factor(V3)
})
sub2<- groupedData(V1~group|subject, data = md_sub2)
attach(sub2)
sub2.mod<-lme(V1~group*time, random =~1|subject, data = sub2)
sub2.anova<-anova(sub2.mod)
sub2.anova
#capture.output(sub2.anova, file = "md21_anova.doc")
capture.output(sub2.anova, file = "md21_three_group_anova.doc")

#ttest (planned comparisons)
t1<-t.test(V1[time==2 & group ==1],V1[time==2 & group==2])  #this is significant
t2<-t.test(V1[time==1 & group ==1],V1[time==1 & group==2])  #this is NOT significant
t3<-t.test(V1[time==3 & group ==1],V1[time==3 & group==2])  #this is NOT significant
capture.output(t3, file = "md21_PC.doc")

sub2_t1<-subset(md_sub2, md_sub2$time==1)
sub2_t2<-subset(md_sub2, md_sub2$time==2) # create subset of 3 levels of time
sub2_t3<-subset(md_sub2, md_sub2$time==3)

# Subject 2 , MD
aov_sub2_t1<-aov(sub2_t1$V1~sub2_t1$group)
aov_sub2_t2<-aov(sub2_t2$V1~sub2_t2$group) # perform the one-way anova for above subsets
aov_sub2_t3<-aov(sub2_t3$V1~sub2_t3$group)

summary(aov_sub2_t1)
summary(aov_sub2_t2)  # summary anova table..
summary(aov_sub2_t3)

tu1<-TukeyHSD(aov_sub2_t1)
tu2<-TukeyHSD(aov_sub2_t2) # finally pairwise comparison for the simple effect using Tukey
tu3<-TukeyHSD(aov_sub2_t3)
capture.output(tu1, file = "md21_tuckey1.doc")
capture.output(tu2, file = "md21_tuckey2.doc")
capture.output(tu3, file = "md21_tuckey3.doc")


summary(aov_sub2_t1)
summary(aov_sub2_t2)  # to check if simple effects are signficant!!
summary(aov_sub2_t3)

library(lsr) # for etaSquared() funtion

etaSquared(aov_sub2_t1, anova= TRUE)
etaSquared(aov_sub2_t2, anova= TRUE)
etaSquared(aov_sub2_t3, anova= TRUE)


#  catch trial
catcht <- within(catcht, {
  subject <- factor(V4)
  group <- factor(3-V2)
  time <- factor(V3)
})
catd<- groupedData(V1~group*time|subject, data = catcht)
attach(catd)
catd.mod<-lme(V1~group*time, random =~1|subject, data = catd)
catd.anova<-anova(catd.mod)
catd.anova
capture.output(catd.anova, file = "catd_anova.doc")

#planned comparison
t.test(V1[time==2 & group ==1],V1[time==2 & group==2])  #this is NOT significant
t.test(V1[time==1 & group ==1],V1[time==1 & group==2])  #this is NOT significant
t.test(V1[time==3 & group ==1],V1[time==3 & group==2])  #this is NOT significant




# interaction force
intforce <- within(intforce, {
  subject <- factor(V4)
  group <- factor(V2)
  time <- factor(V3)
})
intd<- groupedData(V1~group*time|subject, data = intforce)
attach(intd)
intd.mod<-lme(V1~group*time, random =~1|subject, data = intd)
intd.anova<-anova(intd.mod)
intd.anova
capture.output(intd.anova, file = "intforce_three_group_anova.doc")

#planned comparison
t1<-t.test(V1[time==2 & group ==1],V1[time==2 & group==2])  #this is significant
t2<-t.test(V1[time==1 & group ==1],V1[time==1 & group==2])  #this is NOT significant
t3<-t.test(V1[time==3 & group ==1],V1[time==3 & group==2])  #this is NOT significant
capture.output(t3, file = "intforce_PC.doc")


intf_t1<-subset(intforce, intforce$time==1)
intf_t2<-subset(intforce, intforce$time==2) # create subset of 3 levels of time
intf_t3<-subset(intforce, intforce$time==3)

# Subject 2 , MD
aov_intf_t1<-aov(intf_t1$V1~intf_t1$group)
aov_intf_t2<-aov(intf_t2$V1~intf_t2$group) # perform the one-way anova for above subsets
aov_intf_t3<-aov(intf_t3$V1~intf_t3$group)

summary(aov_intf_t1)
summary(aov_intf_t2)  # summary anova table..
summary(aov_intf_t3)

tu1<-TukeyHSD(aov_intf_t1)
tu2<-TukeyHSD(aov_intf_t2) # finally pairwise comparison for the simple effect using Tukey
tu3<-TukeyHSD(aov_intf_t3)
capture.output(tu1, file = "intforce_tuckey1.doc")
capture.output(tu2, file = "intforce_tuckey2.doc")
capture.output(tu3, file = "intforce_tuckey3.doc")

summary(aov_intf_t1)
summary(aov_intf_t2)  # to check if simple effects are signficant!!
summary(aov_intf_t3)

library(lsr) # for etaSquared() funtion

etaSquared(aov_intf_t1, anova= TRUE)
etaSquared(aov_intf_t2, anova= TRUE)
etaSquared(aov_intf_t3, anova= TRUE)



# delta crossing times ~ via point1
dct1 <- within(dct1, {
  subject <- factor(V4)
  group <- factor(V2)
  time <- factor(V3)
})
dct1d<- groupedData(V1~group*time|subject, data = dct1)
attach(dct1d)
dct1d.mod<-lme(V1~group*time, random =~1|subject, data = dct1d)
dct1d.anova<-anova(dct1d.mod)
dct1d.anova
capture.output(dct1d.anova, file = "dct1_anova.doc")

#planned comparison
t1<-t.test(V1[time==2 & group ==1],V1[time==2 & group==2])  #this is NOT significant
t2<-t.test(V1[time==1 & group ==1],V1[time==1 & group==2])  #this is NOT significant
t3<-t.test(V1[time==3 & group ==1],V1[time==3 & group==2])  #this is  significant
capture.output(t3, file = "dct1_PC.doc")





# delta crossing times ~ via point2
dct2 <- within(dct2, {
  subject <- factor(V4)
  group <- factor(V2)
  time <- factor(V3)
})
dct2d<- groupedData(V1~group*time|subject, data = dct2)
attach(dct2d)
dct2d.mod<-lme(V1~group*time, random =~1|subject, data = dct2d)
dct2d.anova<-anova(dct2d.mod)
dct2d.anova
capture.output(dct2d.anova, file = "dct2_anova.doc")

#planned comparison
t1<-t.test(V1[time==2 & group ==1],V1[time==2 & group==2])  #this is NOT significant
t2<-t.test(V1[time==1 & group ==1],V1[time==1 & group==2])  #this is NOT significant
t3<-t.test(V1[time==3 & group ==1],V1[time==3 & group==2])  #this is NOT significant
capture.output(dct2d.anova, file = "dct2_PC.doc")

#planned comparison
t1<-t.test(V1[time==2 & group ==1],V1[time==2 & group==2])  #this is NOT significant
t2<-t.test(V1[time==1 & group ==1],V1[time==1 & group==2])  #this is NOT significant
t3<-t.test(V1[time==3 & group ==1],V1[time==3 & group==2])  #this is  significant
capture.output(t3, file = "dct1_PC.doc")

dct2_t1<-subset(dct2, dct2$time==1)
dct2_t2<-subset(dct2, dct2$time==2) # create subset of 3 levels of time
dct2_t3<-subset(dct2, dct2$time==3)

# Subject 2 , DCT
aov_dct2_t1<-aov(dct2_t1$V1~dct2_t1$group)
aov_dct2_t2<-aov(dct2_t2$V1~dct2_t2$group) # perform the one-way anova for above subsets
aov_dct2_t3<-aov(dct2_t3$V1~dct2_t3$group)

summary(aov_dct2_t1)
summary(aov_dct2_t2)  # summary anova table..
summary(aov_dct2_t3)

tu1<-TukeyHSD(aov_dct2_t1)
tu2<-TukeyHSD(aov_dct2_t2) # finally pairwise comparison for the simple effect using Tukey
tu3<-TukeyHSD(aov_dct2_t3)



# lead1
lead11 <- within(lead11, {
  subject <- factor(V4)
  group <- factor(V2)
  time <- factor(V3)
})
lead11d<- groupedData(V1~group*time|subject, data = lead11)
attach(lead11d)
lead11d.mod<-lme(V1~group*time, random =~1|subject, data = lead11d)
lead11d.anova<-anova(lead11d.mod)
lead11d.anova

#planned comparison
t.test(V1[time==3 & group ==1],V1[time==3 & group==2])  #this is NOT significant
t.test(V1[time==3 & group ==1],V1[time==3 & group==3])  #this is NOT significant
t.test(V1[time==3 & group ==2],V1[time==3 & group==3])  #this is NOT significant

# lead2
lead12 <- within(lead12, {
  subject <- factor(V4)
  group <- factor(V2)
  time <- factor(V3)
})
lead12d<- groupedData(V1~group*time|subject, data = lead12)
attach(lead12d)
lead12d.mod<-lme(V1~group*time, random =~1|subject, data = lead12d)
lead12d.anova<-anova(lead12d.mod)
lead12d.anova

#planned comparison
t.test(V1[time==3 & group ==1],V1[time==3 & group==2])  #this is NOT significant
t.test(V1[time==3 & group ==1],V1[time==3 & group==3])  #this is NOT significant
t.test(V1[time==3 & group ==2],V1[time==3 & group==3])  #this is NOT significant

lead21 <- within(lead21, {
  subject <- factor(V4)
  group <- factor(V2)
  time <- factor(V3)
})
lead21d<- groupedData(V1~group*time|subject, data = lead21)
attach(lead21d)
lead21d.mod<-lme(V1~group*time, random =~1|subject, data = lead21d)
lead21d.anova<-anova(lead21d.mod)
lead21d.anova

#planned comparison
t.test(V1[time==3 & group ==1],V1[time==3 & group==2])  #this is NOT significant
t.test(V1[time==3 & group ==1],V1[time==3 & group==3])  #this is NOT significant
t.test(V1[time==3 & group ==2],V1[time==3 & group==3])  #this is NOT significant

# lead2
lead22 <- within(lead22, {
  subject <- factor(V4)
  group <- factor(V2)
  time <- factor(V3)
})
lead22d<- groupedData(V1~group*time|subject, data = lead22)
attach(lead22d)
lead22d.mod<-lme(V1~group*time, random =~1|subject, data = lead22d)
lead22d.anova<-anova(lead22d.mod)
lead22d.anova

#planned comparison
t.test(V1[time==3 & group ==1],V1[time==3 & group==2])  #this is NOT significant
t.test(V1[time==3 & group ==1],V1[time==3 & group==3])  #this is NOT significant
t.test(V1[time==3 & group ==2],V1[time==3 & group==3])  #this is NOT significant

# post-hoc for lead index
#lead 11
lead11_t1<-subset(lead11, lead11$time==1)
lead11_t2<-subset(lead11, lead11$time==2) # create subset of 3 levels of time
lead11_t3<-subset(lead11, lead11$time==3)

aov_lead11_t1<-aov(lead11_t1$V1~lead11_t1$group)
aov_lead11_t2<-aov(lead11_t2$V1~lead11_t2$group) # perform the one-way anova for above subsets
aov_lead11_t3<-aov(lead11_t3$V1~lead11_t3$group)

summary(aov_lead11_t1)
summary(aov_lead11_t2)  # summary anova table..
summary(aov_lead11_t3)

tu1<-TukeyHSD(aov_lead11_t1)
tu2<-TukeyHSD(aov_lead11_t2) # finally pairwise comparison for the simple effect using Tukey
tu3<-TukeyHSD(aov_lead11_t3)

#lead 12
lead12_t1<-subset(lead12, lead12$time==1)
lead12_t2<-subset(lead12, lead12$time==2) # create subset of 3 levels of time
lead12_t3<-subset(lead12, lead12$time==3)

aov_lead12_t1<-aov(lead12_t1$V1~lead12_t1$group)
aov_lead12_t2<-aov(lead12_t2$V1~lead12_t2$group) # perform the one-way anova for above subsets
aov_lead12_t3<-aov(lead12_t3$V1~lead12_t3$group)

summary(aov_lead12_t1)
summary(aov_lead12_t2)  # summary anova table..

summary(aov_lead12_t3)

tu1<-TukeyHSD(aov_lead12_t1)
tu2<-TukeyHSD(aov_lead12_t2) # finally pairwise comparison for the simple effect using Tukey
tu3<-TukeyHSD(aov_lead12_t3)

#lead 21
lead21_t1<-subset(lead21, lead21$time==1)
lead21_t2<-subset(lead21, lead21$time==2) # create subset of 3 levels of time
lead21_t3<-subset(lead21, lead21$time==3)

aov_lead21_t1<-aov(lead21_t1$V1~lead21_t1$group)
aov_lead21_t2<-aov(lead21_t2$V1~lead21_t2$group) # perform the one-way anova for above subsets
aov_lead21_t3<-aov(lead21_t3$V1~lead21_t3$group)

summary(aov_lead21_t1)
summary(aov_lead21_t2)  # summary anova table..
summary(aov_lead21_t3)

tu1<-TukeyHSD(aov_lead21_t1)
tu2<-TukeyHSD(aov_lead21_t2) # finally pairwise comparison for the simple effect using Tukey
tu3<-TukeyHSD(aov_lead21_t3)

#lead 22
lead22_t1<-subset(lead22, lead22$time==1)
lead22_t2<-subset(lead22, lead22$time==2) # create subset of 3 levels of time
lead22_t3<-subset(lead22, lead22$time==3)

aov_lead22_t1<-aov(lead22_t1$V1~lead22_t1$group)
aov_lead22_t2<-aov(lead22_t2$V1~lead22_t2$group) # perform the one-way anova for above subsets
aov_lead22_t3<-aov(lead22_t3$V1~lead22_t3$group)

summary(aov_lead22_t1)
summary(aov_lead22_t2)  # summary anova table..
summary(aov_lead22_t3)

tu1<-TukeyHSD(aov_lead22_t1)
tu2<-TukeyHSD(aov_lead22_t2) # finally pairwise comparison for the simple effect using Tukey
tu3<-TukeyHSD(aov_lead22_t3)



# Score 1
score1 <- within(score1, {
  subject <- factor(V4)
  group <- factor(V2)
  time <- factor(V3)
})
score1d<- groupedData(V1~group*time|subject, data = score1)
attach(score1d)
score1d.mod<-lme(V1~group*time, random =~1|subject, data = score1d)
score1d.anova<-anova(score1d.mod)
score1d.anova
capture.output(score1d.anova, file = "score1_three_group_anova.doc")

#planned comparison
t1<-t.test(V1[time==2 & group ==1],V1[time==2 & group==2])  #this is NOT significant
t2<-t.test(V1[time==1 & group ==1],V1[time==1 & group==2])  #this is NOT significant
t3<-t.test(V1[time==3 & group ==1],V1[time==3 & group==2])  #this is NOT significant
capture.output(t3, file = "score1_PC.doc")

score1_t1<-subset(score1, score1$time==1)
score1_t2<-subset(score1, score1$time==2) # create subset of 3 levels of time
score1_t3<-subset(score1, score1$time==3)

# Subject 2 , MD
aov_score1_t1<-aov(score1_t1$V1~score1_t1$group)
aov_score1_t2<-aov(score1_t2$V1~score1_t2$group) # perform the one-way anova for above subsets
aov_score1_t3<-aov(score1_t3$V1~score1_t3$group)

summary(aov_score1_t1)
summary(aov_score1_t2)  # summary anova table..
summary(aov_score1_t3)

tu1<-TukeyHSD(aov_score1_t1)
tu2<-TukeyHSD(aov_score1_t2) # finally pairwise comparison for the simple effect using Tukey
tu3<-TukeyHSD(aov_score1_t3)
capture.output(tu1, file = "score1_tuckey1.doc")
capture.output(tu2, file = "score1_tuckey2.doc")
capture.output(tu3, file = "score1_tuckey3.doc")


summary(aov_score1_t1)
summary(aov_score1_t2)  #to check if simple effects are signficant!!
summary(aov_score1_t3)

library(lsr) # for etaSquared() funtion

etaSquared(aov_score1_t1, anova= TRUE)
etaSquared(aov_score1_t2, anova= TRUE)
etaSquared(aov_score1_t3, anova= TRUE)


# Score 2
score2 <- within(score2, {
  subject <- factor(V4)
  group <- factor(V2)
  time <- factor(V3)
})
score2d<- groupedData(V1~group*time|subject, data = score2)
attach(score2d)
score2d.mod<-lme(V1~group*time, random =~1|subject, data = score2d)
score2d.anova<-anova(score2d.mod)
score2d.anova
capture.output(score2d.anova, file = "score2_anova.doc")

#planned comparison
t1<-t.test(V1[time==2 & group ==1],V1[time==2 & group==2])  #this is NOT significant
t2<-t.test(V1[time==1 & group ==1],V1[time==1 & group==2])  #this is NOT significant
t3<-t.test(V1[time==3 & group ==1],V1[time==3 & group==2])  #this is NOT significant
capture.output(t3, file = "score2_PC.doc")

score2_t1<-subset(score2, score2$time==1)
score2_t2<-subset(score2, score2$time==2) # create subset of 3 levels of time
score2_t3<-subset(score2, score2$time==3)

# Subject 2 , MD
aov_score2_t1<-aov(score2_t1$V1~score2_t1$group)
aov_score2_t2<-aov(score2_t2$V1~score2_t2$group) # perform the one-way anova for above subsets
aov_score2_t3<-aov(score2_t3$V1~score2_t3$group)

summary(aov_score2_t1)
summary(aov_score2_t2)  # summary anova table..
summary(aov_score2_t3)

TukeyHSD(aov_score2_t1)
TukeyHSD(aov_score2_t2) # finally pairwise comparison for the simple effect using Tukey
TukeyHSD(aov_score2_t3)

summary(aov_score2_t1)
summary(aov_score2_t2)  # to check if simple effects are signficant!!
summary(aov_score2_t3)

library(lsr) # for etaSquared() funtion

etaSquared(aov_score2_t1, anova= TRUE)
etaSquared(aov_score2_t2, anova= TRUE)
etaSquared(aov_score2_t3, anova= TRUE)



# visualize summary statistics using error bar ~ various ways of visualising 
#Using GGPLOT - "Grammer of graphics"!!

library(ggplot2)
library(ggsignif)

#define position
posn.d <- position_dodge(width = 0.1)
posn.jd <- position_jitterdodge(jitter.width = 0.1, dodge.width = 0.2)
posn.j <- position_jitter(width = 0.2)

parvia.graph<-ggplot(parvia, aes(x= factor(V3), y= V1, col = factor(V2), group = factor(V2)))
ownvia.graph<-ggplot(ownvia, aes(x= factor(V3), y= V1, col = factor(V2), group = factor(V2)))
intforce.graph<-ggplot(intforce, aes(x= factor(V2), y= V1, col = factor(V3), group = factor(V3)))


#Interaction force
# Plot 1: Mean and SD - with T-tipped error bars
intforce.graph +stat_summary(geom = "point", fun.y = mean, position = posn.d) +
  stat_summary(geom = "errorbar", fun.data = mean_sdl,position = posn.d, fun.args = list(mult = 1), width = 0.1)+
  labs(x="Exposure", y="Interaction force (N)")+
  scale_color_manual(labels = c("Early", "Middle", "Late"), values = c("blue", "red", "green"))+
  scale_x_discrete("Training", labels = c("H", "VH", "PV"))+
  ggtitle("title goes here") + theme_bw()+
  guides(color=guide_legend("Group"))

# Plot 2: Jittered, dodged scatter plot with transparent points
intforce.graph +geom_point(position = posn.jd, alpha = 0.6)+
  labs(x="Exposure", y="Interaction force (N)")+
  scale_color_manual(labels = c("Early", "Middle", "Late"), values = c("blue", "red", "green"))+
  scale_x_discrete("Training", labels = c("H", "VH", "PV"))+
  ggtitle("title goes here") + theme_bw()+
  guides(color=guide_legend("Group"))

# Plot 3: Mean and 95% CI - the easy way
intforce.graph +stat_summary(fun.data = mean_cl_normal,  position = posn.d)+
  labs(x="Exposure", y="Distance (m)")+
  scale_color_manual(labels =  c("Early", "Middle", "Late"), values = c("blue", "red", "green"))+
  scale_x_discrete("Training", labels = c("H", "VH", "PV"))+
  ggtitle("title goes here") + theme_bw()+
  guides(color=guide_legend("Group"))

# Plot 4: Mean and SD - the easy way
intforce.graph +stat_summary(fun.data = mean_sdl, fun.args = list(mult =1), position = posn.d)+
  labs(x="Exposure", y="Distance (m)")+
  scale_color_manual(labels = c("Early", "Middle", "Late"), values = c("blue", "red", "green"))+
  scale_x_discrete("Training", labels = c("H", "VH", "PV"))+
  ggtitle("title goes here") + theme_bw()+
  guides(color=guide_legend("Group"))





#Own via-point distance
df1 <- data.frame(a = c(1,1,2,2), b = c(0.008, 0.0082, 0.0082, 0.008))
ownvia.graph +stat_summary(geom = "point", fun.y = mean, position = posn.d) +
stat_summary(geom = "errorbar", fun.data = mean_sdl,position = posn.d, fun.args = list(mult = 1), width = 0.1)+
labs(x="Exposure", y="Distance (m)")+
scale_color_manual(labels = c("H", "VH"), values = c("blue", "red"))+
scale_x_discrete("Training", labels = c("Pre", "Post"))+
ggtitle("title goes here") + theme_bw()+
#geom_line(aes(x = c(1,1,2,2), y = c(0.008, 0.0082, 0.0082, 0.008)))+
annotate("text", x = 1, y = .008, label = "*", size = 8)+
guides(color=guide_legend("Group"))




# to save image as 
figure1<-last_plot()
ggsave(figure1,file="figure1.eps",width=4,height=4,dpi = 400)



#  welch two sample t-test for each fitting coeffs

#force
t.test(V2~V1, data=intf_fit)

#deltatc
t.test(V2~V1, data=deltafit)
t.test(V3~V1, data=deltafit)

#md
t.test(V2~V1, data=mdfit)
t.test(V3~V1, data=mdfit)

t.test(V2~V1, data=mindist1)

boxplot(V2~V1, data = deltafit)

# Wilcoxin-Mann-Whitney test
wilcox.test(V2 ~ V1, data=deltafit)


#significance testing on model parameters

#one column
model.intf.lqg<-read.table("H_model_intforce_LQG.txt")
model.intf.nash<-read.table("H_model_intforce_Nash.txt")
model.intp.lqg<-read.table("H_model_intpower_LQG.txt")
model.intp.nash<-read.table("H_model_intpower_Nash.txt")

#4 column
# 11, 22, 12, 21
model.md.lqg<-read.table("H_model_md_all_LQG.txt")
model.md.nash<-read.table("H_model_md_all_Nash.txt")
model.delta.lqg<-read.table("H_model_delta_all_LQG.txt")
model.delta.nash<-read.table("H_model_delta_all_Nash.txt")

#leadindex
model.lead.nash_1<-read.table("H_model_intpowerIndex1_1.txt")
model.lead.nash_2<-read.table("H_model_intpowerIndex2_1.txt")
model.lead.lqg_1<-read.table("H_model_intpowerIndex1_7.txt")
model.lead.lqg_2<-read.table("H_model_intpowerIndex2_7.txt")


#t-Test to compare the means of two groups 
#under the assumption that both samples are random, 
#independent, and come from normally 
#distributed population with unknow but equal variances

#Before proceeding with the t-test, 
#it is necessary to evaluate the sample variances 
#of the two groups, using a Fisher's F-test 
#to verify the homoskedasticity (homogeneity of variances).

var.test(model.intf.lqg$V1,model.intf.nash$V1)
var.test(model.intp.lqg$V1,model.intp.nash$V1)

# md12, delta 1
var.test(model.md.lqg$V3,model.md.nash$V4)
var.test(model.delta.lqg$V3,model.delta.nash$V3)

#md21, delta 2
var.test(model.md.lqg$V4,model.md.nash$V3)
var.test(model.delta.lqg$V4,model.delta.nash$V4)

# if you get  p-value greater than 0.05, 
#then we can assume that the two variances 
#are homogeneous. 
t.test(model.intf.lqg$V1,model.intf.nash$V1, var.equal=TRUE, paired=FALSE)
t.test(model.intp.lqg$V1,model.intp.nash$V1, var.equal=TRUE, paired=FALSE)
t.test(model.md.lqg$V3,model.md.nash$V3, var.equal=TRUE, paired=FALSE)
t.test(model.md.lqg$V4,model.md.nash$V4, var.equal=TRUE, paired=FALSE)


# delta crossing times

t.test((model.delta.lqg$V4-model.delta.lqg$V1),(model.delta.nash$V4-model.delta.nash$V1),var.equal = TRUE, paired = FALSE)
t.test((model.delta.lqg$V3-model.delta.lqg$V2),(model.delta.nash$V3-model.delta.nash$V2),var.equal = TRUE, paired = FALSE)



# visulaize using box plot

boxplot(model.md.lqg$V3,model.md.nash$V3,xlab="model", ylab="minimum distnace", names= c("LQG", "Nash"))

# comparing mid training effect in indicators

t.test(V3~V1,paired=FALSE,var.equal=TRUE, data=intforce1)
t.test(V9~V1,paired=FALSE,var.equal=TRUE, data=mindist1)


#model _  lead index

t.test(model.lead.lqg_1$V1,model.lead.nash_1$V1, var.equal=TRUE, paired=FALSE)
t.test(model.lead.lqg_1$V2,model.lead.nash_1$V2, var.equal=TRUE, paired=FALSE)
t.test(model.lead.lqg_2$V1,model.lead.nash_2$V1, var.equal=TRUE, paired=FALSE)
t.test(model.lead.lqg_2$V2,model.lead.nash_2$V2, var.equal=TRUE, paired=FALSE)


# steps to perform posthoc analysis




