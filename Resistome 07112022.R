#The Gut Resistome during Hematopoietic Stem Cell Transplantation in a Pediatric Cohort
#Heston SM, Young RR, Jenkins K, Martin PL, Stokhuyzen A, Ward DV, Bhattarai SK, Bucci V, Arshad M, Chao NJ, Seed PC, Kelly MS.
#7/11/2022

library(tidyverse)
library(ggplot2)
library(dplyr)
library(phyloseq)
library(BiocGenerics)
library(lme4)
library(lmerTest)
library(MASS)
library(glmmADMB)
library(RColorBrewer)
library(rmcorr)
library(scales)
library(reshape2)
library(ggpubr)
library(ICC)

remove(list=ls())
set.seed(1234)

#####DATA PREPARATION#####
#RESISTOME
phy.card <- readRDS("T:/Projects/PID/Kelly/Intestinal Microbiome of  - Pro00064365/Resistome Analysis/Raw Data/phy.card.04042022_pruned_500k_paired_reads.rds")
phy.card <- prune_taxa(taxa_sums(phy.card)>0, phy.card)
otu.table <- data.frame(otu_table(phy.card))
metadata_card <- data.frame(sample_data(phy.card))
num <- metadata_card %>% distinct(study_id, .keep_all = TRUE) #80
tax.table<- data.frame(tax_table(phy.card)) ## 350 unique ARG

#Table 1
summary(num$age)
table(num$sex)
table(num$diagnosis)
table(num$type_hsct)
table(num$source_hsct)
filter(num, source_hsct == "Other")
table(num$hsct_prep_cat)

sample.count <- metadata_card %>% count(study_id)
summary(sample.count$n) #median number of samples per child 8 (6, 12) Range 1-21

sum(metadata_card$R1_metag_paired) #5624857872 paired sequencing reads
summary(metadata_card$R1_metag_paired) #Paired median 7452427   (3412479  , 12023062 )

#creating tax table with arg antibiotic class
arg.raw <- read.csv("T:/Projects/PID/Kelly/Intestinal Microbiome of  - Pro00064365/Resistome Analysis/Raw Data/aro_index_03282022.csv") #from CARD
arg_card <- data.frame(otu_table(phy.card))
arg_card <- tibble::rownames_to_column(arg_card, var = "concated_column")
arg_card$concated_column[arg_card$concated_column == "gb|AEQ93536_1|ARO_3004468|CrpP"] <- "gb|AEQ93536_1|ARO_3004467|CrpP"
arg_gene <- tidyr::separate(arg_card, "concated_column", c("A", "B", "C", "D"), sep = "\\|", extra = "merge", remove = FALSE)

arg.raw$ARO.Accession <- stringr::str_replace_all(arg.raw$ARO.Accession, ":", "_")
arg.raw$ARO.Name <- stringr::str_replace_all(arg.raw$ARO.Name, "'", "")
arg.table <- arg.raw %>% dplyr::select(Resistance.Mechanism, Drug.Class, AMR.Gene.Family, ARO.Accession)
arg.table <- arg.table[,c(4,3,2,1)]

merge <- left_join(arg_gene, arg.table, by =c("C"="ARO.Accession"), copy = FALSE)
merge <- merge %>% dplyr::select(concated_column, C, Resistance.Mechanism, Drug.Class, AMR.Gene.Family)
merge$Resistance.Mechanism[merge$C == "ARO_3004467"] <- "antibiotic inactivation"
merge$Drug.Class[merge$C == "ARO_3004467"] <- "fluoroquinolone antibiotic"

rownames(merge) <- merge$concated_column
merge <- merge %>% dplyr::select(C, Resistance.Mechanism, Drug.Class, AMR.Gene.Family)
merge <- merge %>% dplyr::rename("ARO.Accession" = "C")
merge <- merge[, c(2, 3, 4, 1)]
remove(arg_card, arg_gene, arg.table, arg.raw)

#Adding clinically relevant antibiotic class indicator variable
merge <- merge %>% separate("Drug.Class", c("class1", "class2", "class3", "class4", "class5", "class6", "class7", "class8", "class9", "class10", "class11", "class12", "class13", "class14", "class15", "class16"), sep = ";", extra = "warn", remove = FALSE)
merge[sapply(merge, is.character)] <- lapply(merge[sapply(merge, is.character)], as.factor)
#Make list of all ARG classes
classes1 <- data.frame(levels(merge$class1))
classes1 <- dplyr::rename(classes1, c("class" = "levels.merge.class1."))
classes2 <- data.frame(levels(merge$class2))
classes2 <- dplyr::rename(classes2, c("class" = "levels.merge.class2."))
classes3 <- data.frame(levels(merge$class3))
classes3 <- dplyr::rename(classes3, c("class" = "levels.merge.class3."))
classes4 <- data.frame(levels(merge$class4))
classes4 <- dplyr::rename(classes4, c("class" = "levels.merge.class4."))
classes5 <- data.frame(levels(merge$class5))
classes5 <- dplyr::rename(classes5, c("class" = "levels.merge.class5."))
classes6 <- data.frame(levels(merge$class6))
classes6 <- dplyr::rename(classes6, c("class" = "levels.merge.class6."))
classes7 <- data.frame(levels(merge$class7))
classes7 <- dplyr::rename(classes7, c("class" = "levels.merge.class7."))
classes8 <- data.frame(levels(merge$class8))
classes8 <- dplyr::rename(classes8, c("class" = "levels.merge.class8."))
classes9 <- data.frame(levels(merge$class9))
classes9 <- dplyr::rename(classes9, c("class" = "levels.merge.class9."))
classes10 <- data.frame(levels(merge$class10))
classes10 <- dplyr::rename(classes10, c("class" = "levels.merge.class10."))
classes11 <- data.frame(levels(merge$class11))
classes11 <- dplyr::rename(classes11, c("class" = "levels.merge.class11."))
classes12 <- data.frame(levels(merge$class12))
classes12 <- dplyr::rename(classes12, c("class" = "levels.merge.class12."))
classes13 <- data.frame(levels(merge$class13))
classes13 <- dplyr::rename(classes13, c("class" = "levels.merge.class13."))
classes14 <- data.frame(levels(merge$class14))
classes14 <- dplyr::rename(classes14, c("class" = "levels.merge.class14."))
classes15 <- data.frame(levels(merge$class15))
classes15 <- dplyr::rename(classes15, c("class" = "levels.merge.class15."))
classes16 <- data.frame(levels(merge$class16))
classes16 <- dplyr::rename(classes16, c("class" = "levels.merge.class16."))
classes <- rbind(classes1, classes2, classes3, classes4, classes5, classes6, classes7, classes8, classes9, classes10, classes11, classes12, classes13, classes14, classes15, classes16)
classes <- distinct(classes)
remove(classes1, classes2, classes3, classes4, classes5, classes6, classes7, classes8, classes9, classes10, classes11, classes12, classes13, classes14, classes15, classes16)
#34 unique ARG classes

#Assigning ARG class to each ARG
merge <- merge %>% rownames_to_column %>% mutate(clinical = case_when(str_detect(AMR.Gene.Family, "van")~ "van", 
                                                                      str_detect(AMR.Gene.Family, "KPC|IMI|SME|OXA|NDM|VIM|IMP|GIM|SPM")~ "carbapenemase",
                                                                      str_detect(AMR.Gene.Family, "SHV|TEM|CTX|PER|GES|VEB|BES|CME|SFO|BEL|TLA")~ "esbl"))
merge <- merge %>% mutate(amino = if_else(str_detect(Drug.Class, "aminoglycoside"), 1, 0)) %>% mutate(carbapenem = if_else(str_detect(Drug.Class, "carbapenem|penem"), 1, 0)) %>% mutate(betalactam = if_else(str_detect(Drug.Class, "cephalosporin|cephamycin|penam|monobactam"), 1, 0)) %>% 
  mutate(trimeth = if_else(str_detect(Drug.Class, "diaminopyrimidine"), 1, 0)) %>% mutate(fluoro = if_else(str_detect(Drug.Class, "fluoroquinolone"), 1, 0)) %>% mutate(macrolide = if_else(str_detect(Drug.Class, "macrolide"), 1, 0)) %>% mutate(phenicol = if_else(str_detect(Drug.Class, "phenicol"), 1, 0)) %>% 
  mutate(sulfa = if_else(str_detect(Drug.Class, "sulfonamide"), 1, 0)) %>% mutate(tetra = if_else(str_detect(Drug.Class, "tetracycline|glycylcycline"), 1, 0)) %>% mutate(clinda = if_else(str_detect(Drug.Class, "lincosamide"), 1, 0)) %>% mutate(vanc = if_else(str_detect(Drug.Class, "glycopeptide"), 1, 0)) %>% 
  mutate(linezolid = if_else(str_detect(Drug.Class, "oxazolidinone"), 1, 0), flagyl = if_else(str_detect(Drug.Class, "nitroimidazole "), 1, 0))
merge <- merge %>% dplyr::select(rowname, Resistance.Mechanism, Drug.Class, AMR.Gene.Family, ARO.Accession, clinical, amino, carbapenem, betalactam,trimeth, fluoro, macrolide, phenicol, sulfa, tetra, clinda, vanc, linezolid, flagyl) %>% column_to_rownames(var = "rowname")

#adding updated ARG classifications to phyloseq
taxa <- as.matrix(merge, rownames.force = NA)
nrow(taxa)
rownames(taxa) <- rownames(otu_table(phy.card))
phy.card<- `tax_table<-`(phy.card, taxa)

#Adding in abundances for each ARG class to metadata
amino.abun <- merge %>% dplyr::select(6)
amino.otu <- cbind(otu.table, amino.abun)
amino.otu <- amino.otu %>% filter(amino==1)
amino.genes <- data.frame(colSums(amino.otu))
amino.genes <- amino.genes %>% dplyr::rename(c("amino.abund" = "colSums.amino.otu.")) %>% rownames_to_column(var="SampleID") %>% separate(SampleID, c(NA, "SampleID"), sep = 1, remove = TRUE)
metad <- data.frame(sample_data(phy.card))
metad <- left_join(metad, amino.genes, by = "SampleID")

carb.abun <- merge %>% dplyr::select(7)
carb.otu <- cbind(otu.table, carb.abun)
carb.otu <- carb.otu %>% filter(carbapenem==1)
carb.genes <- data.frame(colSums(carb.otu))
carb.genes <- carb.genes %>% dplyr::rename(c("carb.abund" = "colSums.carb.otu.")) %>% rownames_to_column(var="SampleID") %>% separate(SampleID, c(NA, "SampleID"), sep = 1, remove = TRUE)
metad <- left_join(metad, carb.genes, by = "SampleID")

beta.abun <- merge %>% dplyr::select(8)
beta.otu <- cbind(otu.table, beta.abun)
beta.otu <- beta.otu %>% filter(betalactam==1)
beta.genes <- data.frame(colSums(beta.otu))
beta.genes <- beta.genes %>% dplyr::rename(c("beta.abund" = "colSums.beta.otu.")) %>% rownames_to_column(var="SampleID") %>% separate(SampleID, c(NA, "SampleID"), sep = 1, remove = TRUE)
metad <- left_join(metad, beta.genes, by = "SampleID")

trimeth.abun <- merge %>% dplyr::select(9)
trimeth.otu <- cbind(otu.table, trimeth.abun)
trimeth.otu <- trimeth.otu %>% filter(trimeth==1)
trimeth.genes <- data.frame(colSums(trimeth.otu))
trimeth.genes <- trimeth.genes %>% dplyr::rename(c("trimeth.abund" = "colSums.trimeth.otu.")) %>% rownames_to_column(var="SampleID") %>% separate(SampleID, c(NA, "SampleID"), sep = 1, remove = TRUE)
metad <- left_join(metad, trimeth.genes, by = "SampleID")

fluoro.abun <- merge %>% dplyr::select(10)
fluoro.otu <- cbind(otu.table, fluoro.abun)
fluoro.otu <- fluoro.otu %>% filter(fluoro==1)
fluoro.genes <- data.frame(colSums(fluoro.otu))
fluoro.genes <- fluoro.genes %>% dplyr::rename(c("fluoro.abund" = "colSums.fluoro.otu.")) %>% rownames_to_column(var="SampleID") %>% separate(SampleID, c(NA, "SampleID"), sep = 1, remove = TRUE)
metad <- left_join(metad, fluoro.genes, by = "SampleID")

macrolide.abun <- merge %>% dplyr::select(11)
macrolide.otu <- cbind(otu.table, macrolide.abun)
macrolide.otu <- macrolide.otu %>% filter(macrolide==1)
macrolide.genes <- data.frame(colSums(macrolide.otu))
macrolide.genes <- macrolide.genes %>% dplyr::rename(c("macrolide.abund" = "colSums.macrolide.otu.")) %>% rownames_to_column(var="SampleID") %>% separate(SampleID, c(NA, "SampleID"), sep = 1, remove = TRUE)
metad <- left_join(metad, macrolide.genes, by = "SampleID")

phenicol.abun <- merge %>% dplyr::select(12)
phenicol.otu <- cbind(otu.table, phenicol.abun)
phenicol.otu <- phenicol.otu %>% filter(phenicol==1)
phenicol.genes <- data.frame(colSums(phenicol.otu))
phenicol.genes <- phenicol.genes %>% dplyr::rename(c("phenicol.abund" = "colSums.phenicol.otu.")) %>% rownames_to_column(var="SampleID") %>% separate(SampleID, c(NA, "SampleID"), sep = 1, remove = TRUE)
metad <- left_join(metad, phenicol.genes, by = "SampleID")

sulfa.abun <- merge %>% dplyr::select(13)
sulfa.otu <- cbind(otu.table, sulfa.abun)
sulfa.otu <- sulfa.otu %>% filter(sulfa==1)
sulfa.genes <- data.frame(colSums(sulfa.otu))
sulfa.genes <- sulfa.genes %>% dplyr::rename(c("sulfa.abund" = "colSums.sulfa.otu.")) %>% rownames_to_column(var="SampleID") %>% separate(SampleID, c(NA, "SampleID"), sep = 1, remove = TRUE)
metad <- left_join(metad, sulfa.genes, by = "SampleID")

tetra.abun <- merge %>% dplyr::select(14)
tetra.otu <- cbind(otu.table, tetra.abun)
tetra.otu <- tetra.otu %>% filter(tetra==1)
tetra.genes <- data.frame(colSums(tetra.otu))
tetra.genes <- tetra.genes %>% dplyr::rename(c("tetra.abund" = "colSums.tetra.otu.")) %>% rownames_to_column(var="SampleID") %>% separate(SampleID, c(NA, "SampleID"), sep = 1, remove = TRUE)
metad <- left_join(metad, tetra.genes, by = "SampleID")

clinda.abun <- merge %>% dplyr::select(15)
clinda.otu <- cbind(otu.table, clinda.abun)
clinda.otu <- clinda.otu %>% filter(clinda==1)
clinda.genes <- data.frame(colSums(clinda.otu))
clinda.genes <- clinda.genes %>% dplyr::rename(c("clinda.abund" = "colSums.clinda.otu.")) %>% rownames_to_column(var="SampleID") %>% separate(SampleID, c(NA, "SampleID"), sep = 1, remove = TRUE)
metad <- left_join(metad, clinda.genes, by = "SampleID")

vanc.abun <- merge %>% dplyr::select(16)
vanc.otu <- cbind(otu.table, vanc.abun)
vanc.otu <- vanc.otu %>% filter(vanc==1)
vanc.genes <- data.frame(colSums(vanc.otu))
vanc.genes <- vanc.genes %>% dplyr::rename(c("vanc.abund" = "colSums.vanc.otu.")) %>% rownames_to_column(var="SampleID") %>% separate(SampleID, c(NA, "SampleID"), sep = 1, remove = TRUE)
metad <- left_join(metad, vanc.genes, by = "SampleID")

linezolid.abun <- merge %>% dplyr::select(17)
linezolid.otu <- cbind(otu.table, linezolid.abun)
linezolid.otu <- linezolid.otu %>% filter(linezolid==1)
linezolid.genes <- data.frame(colSums(linezolid.otu))
linezolid.genes <- linezolid.genes %>% dplyr::rename(c("linezolid.abund" = "colSums.linezolid.otu.")) %>% rownames_to_column(var="SampleID") %>% separate(SampleID, c(NA, "SampleID"), sep = 1, remove = TRUE)
metad <- left_join(metad, linezolid.genes, by = "SampleID")


flagyl.abun <- merge %>% dplyr::select(18)
flagyl.otu <- cbind(otu.table, flagyl.abun)
flagyl.otu <- flagyl.otu %>% filter(flagyl==1)
flagyl.genes <- data.frame(colSums(flagyl.otu))
flagyl.genes <- flagyl.genes %>% dplyr::rename(c("flagyl.abund" = "colSums.flagyl.otu.")) %>% rownames_to_column(var="SampleID") %>% separate(SampleID, c(NA, "SampleID"), sep = 1, remove = TRUE)
metad <- left_join(metad, flagyl.genes, by = "SampleID")
remove(amino.abun, amino.genes, amino.otu, beta.abun, beta.genes, beta.otu, carb.abun, carb.genes, carb.otu, clinda.abun, clinda.genes, clinda.otu, flagyl.abun, flagyl.genes, flagyl.otu, fluoro.abun, fluoro.genes, fluoro.otu,
       linezolid.abun, linezolid.genes, linezolid.otu, macrolide.abun, macrolide.genes, macrolide.otu, phenicol.abun, phenicol.genes, phenicol.otu, sulfa.abun, sulfa.genes, sulfa.otu, tetra.abun, tetra.genes, tetra.otu, trimeth.abun,
       trimeth.genes, trimeth.otu, vanc.abun, vanc.genes, vanc.otu)

#adding total ARG abundance to metadata
total.genes <- data.frame(colSums(otu.table))
total.genes <- total.genes %>% dplyr::rename(c("total.abund" = "colSums.otu.table.")) %>% rownames_to_column(var="SampleID") %>% separate(SampleID, c(NA, "SampleID"), sep = 1, remove = TRUE)
metad <- left_join(metad, total.genes, by = "SampleID")
remove(total.genes)
#determining most common ARGs
most.abund <- metad %>% dplyr::select(47:59)
most.abund <- data.frame(colSums(most.abund))
most.abund <- most.abund %>% dplyr::rename("total.abund" = "colSums.most.abund.") %>% dplyr::arrange(desc(total.abund)) %>% mutate(r.abund = round(total.abund))
head(most.abund, 13) #most abundant ARG classes: tetra, betalactam, fluoro, amino, macrolide, phenicol, clinda, carb, vanc, trimeth, sulfa, flagyl, linezolid
remove(most.abund)
#Calculate number of ARGs per sample
newdata <- data.frame(otu_table(phy.card))
for(i in 1: ncol(newdata)){newdata[,i] <- ifelse(newdata[,i]  >0, 1, 0)}

#number of ARGs per sample
tot.num <- data.frame(colSums(newdata))
tot.num <- tot.num %>% dplyr::rename(num.ARG = colSums.newdata.) %>% rownames_to_column(var = "sample") %>% separate("sample", into = c("X", "SampleID"), sep = 1, extra = "merge") %>% dplyr::select(SampleID, num.ARG)
summary(tot.num$num.ARG) #Median (IQR) ARGs per sample: 28 (15, 56)

#Number of ARGs per sample
common <- data.frame(rowSums(newdata))
common <- common %>% dplyr::arrange(desc(rowSums.newdata.)) 
hist(common$rowSums.newdata.)
common$percent = (common$rowSums.newdata./693)*100 #tetW, tetO, ErmB, dfrF most common >60% of samples 
remove(common)
# Number of ARGs detected in each class.
class.num <- merge %>% dplyr::select(.,6:18)
class.num <- data.frame(colSums(class.num))
remove(class.num)

#number of newly acquired ARGs since the previous sample
new <- data.frame(t(newdata)) %>% rownames_to_column(var  = "sammple") %>% separate("sammple", into = c("X", "SampleID"), sep = 1, extra = "merge") %>%
  separate("SampleID", into = c("study.id", "ST", "sample.num"), sep = c(2, 4), remove = FALSE, extra = "merge") %>% dplyr::select(-c(X, ST))
com.id <- metad %>% dplyr::select(SampleID, day)
new1 <- left_join(new, com.id, by = "SampleID")
new1 <- new1 %>% arrange(study.id, day) %>% group_by(study.id)
new2 <- new1
new4 <- new2 %>% arrange(study.id, day)
new6 <- new2 %>% arrange(study.id, day)
groups(new2)
for(i in 4: 353){new2[,i] <- ifelse((new2[,i] == 1) & (lag(new2[,i]) == 0) & (new2$study.id == lag(new2$study.id)), 1, 0)}
new2$gb.AAD01868_1.ARO_3002867.dfrF[new2$SampleID == "01ST001"] <- 0
new2$new.ARG <- rowSums(new2[,4:353])
na <- filter(new2, is.na(new.ARG))
new3 <- new2 %>% ungroup() %>% dplyr::select(SampleID, new.ARG)

#instability since last sample = jaccard distance
for(i in 4: 353){new4[,i] <- ifelse((new4[,i] == 0) & (lag(new4[,i]) == 1) & (new4$study.id == lag(new4$study.id)), 1, 0)}
new4[1, 4:353] <- 0
new4$lost.ARG <- rowSums(new4[,4:353])
new5 <- new4 %>% ungroup() %>% dplyr::select(SampleID, lost.ARG)
for(i in 4: 353){new6[,i] <- ifelse((new6[,i] == 1) & (lag(new6[,i]) == 1) & (new6$study.id == lag(new6$study.id)), 1, 0)}
new6[1, 4:353] <- 0
new6$stable.ARG <- rowSums(new6[4:353])
new7 <- new6 %>% ungroup() %>% dplyr::select(SampleID, stable.ARG, day, study.id)

instability <- left_join(new3, new5, by = "SampleID")
instability <- left_join(instability, new7, by = "SampleID")
instability <- left_join(instability, tot.num, by = "SampleID")
instability <- instability %>% arrange(study.id, day) %>% group_by(study.id)
instability <- instability %>% mutate(inst.score = (1-((stable.ARG)/(stable.ARG + new.ARG + lost.ARG))), week = if_else(day <= -3, -1,
                                                                                                                if_else(day >-3 & day <=3, 0,
                                                                                                                if_else(day >3 & day <= 10, 1,
                                                                                                                if_else(day >10 & day <=17, 2,
                                                                                                                if_else(day >17 & day <=24, 3,
                                                                                                                if_else(day >24 & day <= 33, 4,
                                                                                                                if_else(day >33 & day <= 40, 5,
                                                                                                                if_else(day >40 & day <= 47, 6, if_else(day >47 & day <= 54, 7, 8))))))))))
instability <- instability %>% ungroup() %>% dplyr::select(SampleID, new.ARG, lost.ARG, stable.ARG, num.ARG, inst.score)
remove(new1, new2, new3, new4, new5, new6, new7)
#Adding in daily antibiotic exposure data
daily <- read.csv("T:/Projects/PID/Kelly/Intestinal Microbiome of  - Pro00064365/Resistome Analysis/Raw Data/pbmt_metadata_daily.csv")
sample.id <- metad %>% mutate(seq =1) %>% dplyr::select(!study_id) %>% dplyr::select(!day)
daily <- daily %>% mutate(study.id = sprintf("%02d", study_id)) %>% mutate(stool.id = str_remove(st, "_")) %>% mutate(SampleID = case_when(stool.id != ""~ paste0(study.id,stool.id)))
# 19ST016, 24ST015, and 26ST060 were the second stool sample collected in same day. Manually writing them in
daily$SampleID[daily$study_id=="19" & daily$day =="9"]<- "19ST016"
daily$SampleID[daily$study_id=="24" & daily$day =="6"]<- "24ST015"
daily$SampleID[daily$study_id=="26" & daily$day =="0"]<- "26ST060"
daily.seq <- left_join(daily, sample.id, by = "SampleID") %>% filter(study_id %in% metad$study_id)

test <- daily.seq %>% dplyr::select(-c(day_study, location, engrafted, st, st2, bl, bl2, or, cu, tpn, bsi_onset, bcx_org, cdiff.x, stool.id, total_contam, pct_metagenomic, primary_dx, donor_hsct, day_engraftment, hla, hsct_prep, gvhd_prophy, prior_probx, mortality_100d, mortality_1y, mortality_2y, agvhd_gutliver, stage_gut, stage_liver, bsi1_day, bsi1, bsi2_day, bsi2, cdiff.y, day_cdiff))
test$seq[is.na(test$seq)] <- 0
test2 <- test %>% mutate(other.predict = tetracycline + other_beta_lactam + clinda + ag, beta.cef.predict = cefepime + other_beta_lactam) %>% mutate(other = if_else(other.predict >0, 1, 0), beta.cef = if_else(beta.cef.predict >0, 1, 0))
test2 <- test2 %>%  dplyr::arrange(study_id, day) %>% group_by(study_id) %>% mutate(seq.int = cumsum(seq)) %>% group_by(study_id, seq.int) 

#Determining antibiotic exposure for each antibiotic
test2.a <- test2  %>%  dplyr::arrange(study_id, day)  %>% group_by(study_id, seq.int) %>% mutate(cefepime.cum = cumsum(cefepime)) %>% ungroup() %>% group_by(study_id) %>% mutate(cefepime.lag= lag(cefepime.cum)) %>% mutate(cefepime.exp= if_else(seq == 1 & cefepime.lag > 0, 1, 0))
check.cefepime <- test2.a %>% filter(day == min(day) & cefepime ==1 & seq ==1) #none were on cefepime on first day with first sequenced sample
test2.a$cefepime.exp[is.na(test2.a$cefepime.exp)] <- 0

test2.b <- test2.a %>%  dplyr::arrange(study_id, day)  %>% group_by(study_id, seq.int) %>% mutate(piptazo.cum = cumsum(piptazo)) %>% ungroup() %>% group_by(study_id) %>% mutate(piptazo.lag= lag(piptazo.cum)) %>% mutate(piptazo.exp= if_else(seq == 1 & piptazo.lag > 0, 1, 0))
check.piptazo <- test2.b %>% filter(day == min(day) & piptazo ==1 & seq ==1) #none were on piptazo on first day with first sequenced sample
test2.b$piptazo.exp[is.na(test2.b$piptazo.exp)] <- 0

test2.c <- test2.b %>%  dplyr::arrange(study_id, day)  %>% group_by(study_id, seq.int) %>% mutate(carbapenem.cum = cumsum(carbapenem)) %>% ungroup() %>% group_by(study_id) %>% mutate(carbapenem.lag= lag(carbapenem.cum)) %>% mutate(carbapenem.exp= if_else(seq == 1 & carbapenem.lag > 0, 1, 0))
check.carbapenem <- test2.c %>% filter(day == min(day) & carbapenem ==1 & seq ==1) #none were on carbapenem on first day with first sequenced sample
test2.c$carbapenem.exp[is.na(test2.c$carbapenem.exp)] <- 0

test2.d <- test2.c  %>%  dplyr::arrange(study_id, day)  %>% group_by(study_id, seq.int) %>% mutate(flagyl.cum = cumsum(flagyl)) %>% ungroup() %>% group_by(study_id) %>% mutate(flagyl.lag= lag(flagyl.cum)) %>% mutate(flagyl.exp= if_else(seq == 1 & flagyl.lag > 0, 1, 0))
check.flagyl <- test2.d %>% filter(day == min(day) & flagyl ==1 & seq ==1) #none were on flagyl on first day with first sequenced sample
test2.d$flagyl.exp[is.na(test2.d$flagyl.exp)] <- 0

test2.d1 <- test2.d  %>%  dplyr::arrange(study_id, day)  %>% group_by(study_id, seq.int) %>% mutate(vanc.cum = cumsum(vanc)) %>% ungroup() %>% group_by(study_id) %>% mutate(vanc.lag= lag(vanc.cum)) %>% mutate(vanc.exp= if_else(seq == 1 & vanc.lag > 0, 1, 0))
check.vanc <- test2.d1 %>% filter(day == min(day) & vanc ==1 & seq ==1) #none were on vanc on first day with first sequenced sample
test2.d1$vanc.exp[is.na(test2.d1$vanc.exp)] <- 0

test2.e <- test2.d1 %>%  dplyr::arrange(study_id, day)   %>% group_by(study_id, seq.int) %>% mutate(quinolone.cum = cumsum(quinolone)) %>% ungroup() %>% group_by(study_id) %>% mutate(quinolone.lag= lag(quinolone.cum)) %>% mutate(quinolone.exp= if_else(seq == 1 & quinolone.lag > 0, 1, 0))
check.quinolone <- test2.e %>% filter(day == min(day) & quinolone ==1 & seq ==1) #none were on quinolone on first day with first sequenced sample
test2.e$quinolone.exp[is.na(test2.e$quinolone.exp)] <- 0

test2.f <- test2.e  %>%  dplyr::arrange(study_id, day) %>% group_by(study_id, seq.int) %>% mutate(tetracycline.cum = cumsum(tetracycline)) %>% ungroup() %>% group_by(study_id) %>% mutate(tetracycline.lag= lag(tetracycline.cum)) %>% mutate(tetracycline.exp= if_else(seq == 1 & tetracycline.lag > 0, 1, 0))
check.tetracycline <- test2.f %>% filter(day == min(day) & tetracycline ==1 & seq ==1) #none were on tetracycline on first day with first sequenced sample
test2.f$tetracycline.exp[is.na(test2.f$tetracycline.exp)] <- 0

test2.g <- test2.f %>%  dplyr::arrange(study_id, day)  %>% group_by(study_id, seq.int) %>% mutate(other_beta_lactam.cum = cumsum(other_beta_lactam)) %>% ungroup() %>% group_by(study_id) %>% mutate(other_beta_lactam.lag= lag(other_beta_lactam.cum)) %>% mutate(other_beta_lactam.exp= if_else(seq == 1 & other_beta_lactam.lag > 0, 1, 0))
check.other_beta_lactam <- test2.g %>% filter(day == min(day) & other_beta_lactam ==1 & seq ==1) #3 were on other_beta_lactam on first day with first sequenced sample
test2.g$other_beta_lactam.exp[is.na(test2.g$other_beta_lactam.exp)] <- 0
test2.g$other_beta_lactam.exp[test2.g$study_id == 36 & test2.g$day == -27 & test2.g$other_beta_lactam == 1] <- 1
test2.g$other_beta_lactam.exp[test2.g$study_id == 47 & test2.g$day == -21 & test2.g$other_beta_lactam == 1] <- 1
test2.g$other_beta_lactam.exp[test2.g$study_id == 50 & test2.g$day == -14 & test2.g$other_beta_lactam == 1] <- 1

test2.h <- test2.g  %>%  dplyr::arrange(study_id, day) %>% group_by(study_id, seq.int) %>% mutate(clinda.cum = cumsum(clinda)) %>% ungroup() %>% group_by(study_id) %>% mutate(clinda.lag= lag(clinda.cum)) %>% mutate(clinda.exp= if_else(seq == 1 & clinda.lag > 0, 1, 0))
check.clinda <- test2.h %>% filter(day == min(day) & clinda ==1 & seq ==1) #none were on clinda on first day with first sequenced sample
test2.h$clinda.exp[is.na(test2.h$clinda.exp)] <- 0

test2.i <- test2.h  %>%  dplyr::arrange(study_id, day) %>% group_by(study_id, seq.int) %>% mutate(ag.cum = cumsum(ag)) %>% ungroup() %>% group_by(study_id) %>% mutate(ag.lag= lag(ag.cum)) %>% mutate(ag.exp= if_else(seq == 1 & ag.lag > 0, 1, 0))
check.ag <- test2.i %>% filter(day == min(day) & ag ==1 & seq ==1) #none were on ag on first day with first sequenced sample
test2.i$ag.exp[is.na(test2.i$ag.exp)] <- 0

test2.j <- test2.i  %>%  dplyr::arrange(study_id, day) %>% group_by(study_id, seq.int) %>% mutate(macrolide.cum = cumsum(macrolide)) %>% ungroup() %>% group_by(study_id) %>% mutate(macrolide.lag= lag(macrolide.cum)) %>% mutate(macrolide.exp= if_else(seq == 1 & macrolide.lag > 0, 1, 0))
check.macrolide <- test2.j %>% filter(day == min(day) & macrolide ==1 & seq ==1) #none were on macrolide on first day with first sequenced sample
test2.j$macrolide.exp[is.na(test2.j$macrolide.exp)] <- 0

test2.k <- test2.j  %>%  dplyr::arrange(study_id, day) %>% group_by(study_id, seq.int) %>% mutate(tmp_smx.cum = cumsum(tmp_smx)) %>% ungroup() %>% group_by(study_id) %>% mutate(tmp_smx.lag= lag(tmp_smx.cum)) %>% mutate(tmp_smx.exp= if_else(seq == 1 & tmp_smx.lag > 0, 1, 0))
check.tmp_smx <- test2.k %>% filter(day == min(day) & tmp_smx ==1 & seq ==1) #5 were on tmp_smx on first day with first sequenced sample
test2.k$tmp_smx.exp[is.na(test2.k$tmp_smx.exp)] <- 0
test2.k$tmp_smx.exp[test2.k$study_id == 3 & test2.k$day == -9 & test2.k$tmp_smx == 1] <- 1
test2.k$tmp_smx.exp[test2.k$study_id == 15 & test2.k$day == -23 & test2.k$tmp_smx == 1] <- 1
test2.k$tmp_smx.exp[test2.k$study_id == 17 & test2.k$day == -7 & test2.k$tmp_smx == 1] <- 1
test2.k$tmp_smx.exp[test2.k$study_id == 54 & test2.k$day == -21 & test2.k$tmp_smx == 1] <- 1
test2.k$tmp_smx.exp[test2.k$study_id == 71 & test2.k$day == -7 & test2.k$tmp_smx == 1] <- 1

test2.l <- test2.k  %>%  dplyr::arrange(study_id, day) %>% group_by(study_id, seq.int) %>% mutate(other.cum = cumsum(other)) %>% ungroup() %>% group_by(study_id) %>% mutate(other.lag= lag(other.cum)) %>% mutate(other.exp= if_else(seq == 1 & other.lag > 0, 1, 0))
check.other <- test2.l %>% filter(day == min(day) & other ==1 & seq ==1) 
test2.l$other.exp[is.na(test2.l$other.exp)] <- 0
test2.l$other.exp[test2.l$study_id == 36 & test2.l$day == -27 & test2.l$other == 1] <- 1
test2.l$other.exp[test2.l$study_id == 47 & test2.l$day == -21 & test2.l$other == 1] <- 1
test2.l$other.exp[test2.l$study_id == 50 & test2.l$day == -14 & test2.l$other == 1] <- 1

test2.m <- test2.l  %>%  dplyr::arrange(study_id, day) %>% group_by(study_id, seq.int) %>% mutate(beta.cef.cum = cumsum(beta.cef)) %>% ungroup() %>% group_by(study_id) %>% mutate(beta.cef.lag= lag(beta.cef.cum)) %>% mutate(beta.cef.exp= if_else(seq == 1 & beta.cef.lag > 0, 1, 0))
check.beta.cef <- test2.m %>% filter(day == min(day) & beta.cef ==1 & seq ==1) #none were on beta.cef on first day with first sequenced sample
test2.m$beta.cef.exp[is.na(test2.m$beta.cef.exp)] <- 0
test2.m$beta.cef.exp[test2.m$study_id == 36 & test2.m$day == -27 & test2.m$beta.cef == 1] <- 1
test2.m$beta.cef.exp[test2.m$study_id == 47 & test2.m$day == -21 & test2.m$beta.cef == 1] <- 1
test2.m$beta.cef.exp[test2.m$study_id == 50 & test2.m$day == -14 & test2.m$beta.cef == 1] <- 1

test3 <- test2.m %>% dplyr::select(-c(contains(".cum"), contains(".lag"), contains("predict"), cefepime, piptazo, carbapenem, vanc, flagyl,
                                      clinda, ag, macrolide, tmp_smx, other_beta_lactam, beta.cef, other, quinolone, tetracycline)) %>% ungroup()
remove(check.cefepime, check.piptazo, check.carbapenem, check.vanc, check.flagyl, check.clinda, check.ag,
       check.macrolide, check.tmp_smx, check.other_beta_lactam, check.beta.cef, check.other, check.quinolone, check.tetracycline,
       test2, test2.a, test2.b, test2.c, test2.d, test2.d1, test2.e, test2.f, test2.g, test2.h, test2.i, test2.j, test2.k, test2.l, test2.m)

test4 <- left_join(test3, instability, by = "SampleID")

#long by patient day with sequencing data
test4 <- test4 %>% filter(seq ==1)

#antibiotic exposures
abx.exp <- test4 %>% summarize(cefepime = sum(cefepime.exp), vanc = sum(vanc.exp), carbapenem = sum(carbapenem.exp), quinolone = sum(quinolone.exp), ag = sum(ag.exp), clinda = sum(clinda.exp), macrolide= sum(macrolide.exp), tmp_smx = sum(tmp_smx.exp), tetracycline = sum(tetracycline.exp), flagyl = sum(flagyl.exp), piptazo = sum(piptazo.exp), other_beta_lactam = sum(other_beta_lactam.exp))
abx.exp1 <- data.frame(t(abx.exp)) %>% dplyr::rename(freq = t.abx.exp.) %>% mutate(pct.sample = round(freq/693 *100, digits = 0)) %>% dplyr::arrange(desc(freq))
#percent of samples exposed to each antibiotic
abx.exp2 <- test4 %>% dplyr::arrange(study_id) %>% group_by(study_id) %>% summarize(cefepime = sum(cefepime.exp), vanc = sum(vanc.exp), carb = sum(carbapenem.exp), quin = sum(quinolone.exp), amino = sum(ag.exp), clinda = sum(clinda.exp), macro= sum(macrolide.exp), tmp = sum(tmp_smx.exp), tetra = sum(tetracycline.exp), flagyl = sum(flagyl.exp), piptazo = sum(piptazo.exp), other_beta = sum(other_beta_lactam.exp)) %>%
  mutate(cef.pt = if_else(cefepime >0, 1, 0), vanc.pt = if_else(vanc >0, 1, 0), carb.pt = if_else(carb >0, 1, 0), quin.pt = if_else(quin >0, 1, 0), amino.pt = if_else(amino >0, 1, 0), clinda.pt = if_else(clinda >0, 1, 0), macro.pt = if_else(macro >0, 1, 0), tmp.pt = if_else(tmp >0, 1, 0), tetra.pt = if_else(tetra >0, 1, 0), flagyl.pt = if_else(flagyl >0, 1, 0), piptazo.pt = if_else(piptazo >0, 1, 0), other_beta.pt = if_else(other_beta >0, 1, 0))
abx.pt3 <- abx.exp2 %>% summarize(cefepime = sum(cef.pt), vanc = sum(vanc.pt), carb = sum(carb.pt), quin = sum(quin.pt), amino = sum(amino.pt), clinda = sum(clinda.pt), macro= sum(macro.pt), tmp = sum(tmp.pt), tetra = sum(tetra.pt), flagyl = sum(flagyl.pt), piptazo = sum(piptazo.pt), other_beta = sum(other_beta.pt))
abx.exp4 <- data.frame(t(abx.pt3)) %>% dplyr::rename(freq = t.abx.pt3.) %>% dplyr::arrange(desc(freq)) %>% mutate(pct.patient = round(freq/80, digits = 2) *100)
#percent of patients exposed to each antibiotic
remove(abx.exp, abx.exp1, abx.exp2, abx.pt3, abx.exp4)

#preparing variables for multivariable analyses
test6 <- test4 %>% mutate(prep = if_else(hsct_prep_cat == "myelo", 1, 2), log.depth = log(total_metagenomic))
cols <- c("study_id", "sex", "race", "diagnosis", "type_hsct", "prep")
test6[cols] <- lapply(test6[cols], factor)
levels(test6$diagnosis)
test6$diagnosis <- relevel(test6$diagnosis, ref = "Heme_malignancy")
table(test6$sex)
levels(test6$type_hsct)
test6$type_hsct <- relevel(test6$type_hsct, ref = "Autologous")

#Antibiotic exposure on resistome
check <- test6 %>% dplyr::select(SampleID, contains(".exp")) 
check <- check %>% column_to_rownames(var = "SampleID")
check$num.abx = rowSums(check[,1:12])
check <- check %>% mutate(any.exp = if_else(num.abx >0, 1, 0))
check <- check %>% rownames_to_column(var = "SampleID") %>% dplyr::select(SampleID, num.abx, any.exp)
test7 <- left_join(test6, check, by = "SampleID")

col1 <- c("cefepime.exp", "piptazo.exp", "carbapenem.exp", "flagyl.exp", "vanc.exp", "quinolone.exp", "tetracycline.exp", "other_beta_lactam.exp", "clinda.exp", "ag.exp", "macrolide.exp", "tmp_smx.exp", "other.exp", "beta.cef.exp", "any.exp")
test7[col1] <- lapply(test7[col1], factor)

#Aerobic vs anaerobic classification
test10 <- test7 %>% mutate(aerobic.exp = if_else(cefepime.exp == 1, 1, if_else(vanc.exp ==1, 1, if_else(quinolone.exp ==1, 1, if_else(ag.exp == 1, 1, if_else(macrolide.exp == 1, 1, if_else(tmp_smx.exp ==1, 1, 0)))))), anaerobic.exp = if_else(piptazo.exp == 1, 1, if_else(carbapenem.exp == 1, 1, if_else(flagyl.exp == 1, 1, if_else(clinda.exp == 1, 1, 0)))))
col2 <- c("aerobic.exp", "anaerobic.exp")
test10[col2] <- lapply(test10[col2], factor)

#MICROBIOME
phy.bmt <- readRDS("T:/Projects/PID/Kelly/Intestinal Microbiome of  - Pro00064365/Resistome Analysis/Raw Data/phy.bmt.004042022_pruned_25K_paired_reads.rds")
p.bmt <- prune_samples(phy.bmt@sam_data$R1_metag_paired > 500000, phy.bmt)
nsamples(p.bmt) #693
ntaxa(p.bmt) #599
otu.bmt <- data.frame(otu_table(p.bmt)) 
tax.bmt <- data.frame(tax_table(p.bmt))
meta.bmt <- data.frame(sample_data(p.bmt))

#cleaning Taxa table and rebuilding phyloseq
tax.bmt <- data.frame(tax_table(p.bmt))
relative_bmt <- psmelt(p.bmt)
species <- aggregate(relative_bmt$Abundance, by=list(Species=relative_bmt$Species), FUN=mean)
species <- arrange(species, desc(x))
topspecies <- unique(species$Species[1:200])
top.200 <- data.frame(topspecies) %>% filter(str_detect(topspecies, "[:digit:]"))
tax.bmt <- rownames_to_column(tax.bmt, var = "OTU") %>% arrange(OTU)

tax.bmt <- tax.bmt %>% filter(OTU != "Bacteroides_faecis_CAG_32")
otu.bmt <- aggregate(otu.bmt, list(Group=replace(rownames(otu.bmt), rownames(otu.bmt) %in% c("Bacteroides_faecis_CAG_32", "Bacteroides_faecis"), "Bacteroides_faecis")), sum) 

tax.bmt$Species[tax.bmt$OTU == "Blautia_sp_CAG_257"] <- "Blautia_caecimuris"
tax.bmt$OTU[tax.bmt$OTU == "Blautia_sp_CAG_257"] <- "Blautia_caecimuris"
otu.bmt$Group[otu.bmt$Group == "Blautia_sp_CAG_257"] <- "Blautia_caecimuris"

tax.bmt$Species[tax.bmt$OTU == "Clostridium_sp_7_2_43FAA"] <- "Clostridium_tertium"
tax.bmt$OTU[tax.bmt$OTU == "Clostridium_sp_7_2_43FAA"] <- "Clostridium_tertium"
otu.bmt$Group[otu.bmt$Group == "Clostridium_sp_7_2_43FAA"] <- "Clostridium_tertium"

tax.bmt$Species[tax.bmt$OTU == "Firmicutes_bacterium_CAG_83"] <- "Vescimonas_coprocola"
tax.bmt$OTU[tax.bmt$OTU == "Firmicutes_bacterium_CAG_83"] <- "Vescimonas_coprocola"
otu.bmt$Group[otu.bmt$Group == "Firmicutes_bacterium_CAG_83"] <- "Vescimonas_coprocola"

tax.bmt$Species[tax.bmt$OTU == "Firmicutes_bacterium_CAG_646"] <- "Faecalibacterium_gallinarum"
tax.bmt$OTU[tax.bmt$OTU == "Firmicutes_bacterium_CAG_646"] <- "Faecalibacterium_gallinarum"
otu.bmt$Group[otu.bmt$Group == "Firmicutes_bacterium_CAG_646"] <- "Faecalibacterium_gallinarum"

tax.bmt$Species[tax.bmt$OTU == "Abiotrophia_sp_HMSC24B09"] <- "Kingella_denitrificans"
tax.bmt$OTU[tax.bmt$OTU == "Abiotrophia_sp_HMSC24B09"] <- "Kingella_denitrificans"
otu.bmt$Group[otu.bmt$Group == "Abiotrophia_sp_HMSC24B09"] <- "Kingella_denitrificans"

tax.bmt <- tax.bmt %>% filter(OTU != "Clostridium_bolteae_CAG_59")
otu.bmt <- column_to_rownames(otu.bmt, var = "Group") 
otu.bmt <- aggregate(otu.bmt, list(Group=replace(rownames(otu.bmt), rownames(otu.bmt) %in% c("Clostridium_bolteae_CAG_59", "Clostridium_bolteae"), "Clostridium_bolteae")), sum) %>% arrange(Group) %>% column_to_rownames(var = "Group")

otu.bmt <- as.matrix(otu.bmt, rownames.force = NA)
nrow(otu.bmt)
OTU=otu_table(otu.bmt, taxa_are_rows = TRUE)
nrow(OTU)
tax.bmt <- tax.bmt %>% arrange(OTU) %>% column_to_rownames(var = "OTU")
tax.bmt <- as.matrix(tax.bmt, rownames.force = NA)
TAX = tax_table(tax.bmt)
nrow(TAX)
meta <- meta.bmt %>% mutate(SampleID=paste0("X", SampleID)) 
rownames(meta) <- NULL
meta <- meta %>% column_to_rownames(var="SampleID")
META = sample_data(meta)
p.bmt <- phyloseq(OTU, TAX, META)

#number of species per sample
num.species <- data.frame(otu_table(p.bmt)) 
for(i in 1: ncol(num.species)){num.species[,i] <- ifelse(num.species[,i] >0, 1, 0)}

spec.num <- data.frame(colSums(num.species))
spec.num <- spec.num %>% rownames_to_column(var = "SampleID") %>% dplyr::rename("num.species" = "colSums.num.species.")
summary(spec.num$num.species) #median (IQR) number of species per sample: 28 (15, 50)
meta.bmt <- meta %>% rownames_to_column(var = "SampleID")
meta_bmt <- left_join(meta.bmt, spec.num, by = "SampleID")

#Newly acquired species since the previous sample
new.sp <- data.frame(t(num.species)) %>% rownames_to_column(var  = "sammple") %>% separate("sammple", into = c("X", "SampleID"), sep = 1, extra = "merge") %>%
  separate("SampleID", into = c("study.id", "ST", "sample.num"), sep = c(2, 4), remove = FALSE, extra = "merge") %>% dplyr::select(-c(X, ST))
com.id <- metad %>% dplyr::select(SampleID, day)
new.sp1 <- left_join(new.sp, com.id, by = "SampleID")
new.sp1 <- new.sp1 %>% arrange(study.id, day) %>% group_by(study.id)
new.sp2 <- new.sp1
new.sp4 <- new.sp1 %>% arrange(study.id, day)
new.sp6 <- new.sp1 %>% arrange(study.id, day)
groups(new.sp2)
for(i in 4: 600){new.sp2[,i] <- ifelse((new.sp2[,i] == 1) & (lag(new.sp2[,i]) == 0) & (new.sp2$study.id == lag(new.sp2$study.id)), 1, 0)}
new.sp2$Ruminococcus_gnavus[new.sp2$SampleID == "01ST001"] <- 0
new.sp2$Terrisporobacter_othiniensis[new.sp2$SampleID == "01ST001"] <- 0
new.sp2$Clostridium_innocuum[new.sp2$SampleID == "01ST001"] <- 0
new.sp2$Erysipelatoclostridium_ramosum[new.sp2$SampleID == "01ST001"] <- 0
new.sp2$new.sp <- rowSums(new.sp2[,4:600])
na <- filter(new.sp2, is.na(new.sp))
new.sp3 <- new.sp2 %>% ungroup() %>% dplyr::select(SampleID, new.sp)
summary(new.sp3$new.sp)
#Instability of bacterial species
for(i in 4:600){new.sp4[,i] <- ifelse((new.sp4[,i] == 0) & (lag(new.sp4[,i]) == 1) & (new.sp4$study.id == lag(new.sp4$study.id)), 1, 0)}
new.sp4[1, 4:600] <- 0
new.sp4$lost.sp <- rowSums(new.sp4[,4:600])
new.sp5 <- new.sp4 %>% ungroup() %>% dplyr::select(SampleID, lost.sp)

for(i in 4:600){new.sp6[,i] <- ifelse((new.sp6[,i] == 1) & (lag(new.sp6[,i]) == 1) & (new.sp6$study.id == lag(new.sp6$study.id)), 1, 0)}
new.sp6[1, 4:600] <- 0
new.sp6$stable.sp <- rowSums(new.sp6[4:600])
new.7sp <- new.sp6 %>% ungroup() %>% dplyr::select(SampleID, stable.sp, day, study.id)

instability.sp <- left_join(new.sp3, new.sp5, by = "SampleID")
instability.sp <- left_join(instability.sp, new.7sp, by = "SampleID")
instability.sp <- left_join(instability.sp, spec.num, by = "SampleID")
instability.sp <- instability.sp %>% arrange(study.id, day) %>% group_by(study.id)
instability.sp <- instability.sp %>% mutate(inst.score = (1-((stable.sp)/(stable.sp + new.sp + lost.sp))), week = if_else(day <= -3, -1,
                                                                                                                  if_else(day >-3 & day <=3, 0,
                                                                                                                  if_else(day >3 & day <= 10, 1,
                                                                                                                  if_else(day >10 & day <=17, 2,
                                                                                                                  if_else(day >17 & day <=24, 3,
                                                                                                                  if_else(day >24 & day <= 33, 4,
                                                                                                                  if_else(day >33 & day <= 40, 5,
                                                                                                                  if_else(day >40 & day <= 47, 6,
                                                                                                                  if_else(day >47 & day <= 54, 7, 8))))))))))

new.sp8 <- instability.sp %>% ungroup() %>% mutate(SampleID = paste0("X", SampleID)) %>% dplyr::select(SampleID, new.sp, inst.score)
meta_bmt <- left_join(meta_bmt, new.sp8, by = "SampleID")
remove(new.sp, new.sp1, new.sp2, new.sp3, new.sp4, new.sp5, new.sp6)
#Preparing variables for analyses
levels(meta_bmt$diagnosis)
meta_bmt$diagnosis <- relevel(meta_bmt$diagnosis, ref = "Heme_malignancy")
levels(meta_bmt$type_hsct)
meta_bmt$type_hsct <- relevel(meta_bmt$type_hsct, ref = "Autologous")
table(meta_bmt$race)
meta_bmt <- meta_bmt %>% mutate(prep = if_else(hsct_prep_cat == "myelo", 1, 2))
cols2 <- c("study_id", "sex", "race", "prep")
meta_bmt[cols2] <- lapply(meta_bmt[cols2], factor)
#adding in abx exposures
abx.exps <- test10 %>% dplyr::select(SampleID, cefepime.exp, piptazo.exp, carbapenem.exp, flagyl.exp, vanc.exp, quinolone.exp, macrolide.exp, tmp_smx.exp, other.exp)
meta_bmt <- left_join(meta_bmt, abx.exps, by = "SampleID") 
meta_bmt <- meta_bmt %>% column_to_rownames(var = "SampleID")

relative_bmt <- psmelt(p.bmt)
species <- aggregate(relative_bmt$Abundance, by=list(Phylum=relative_bmt$Phylum, Genus = relative_bmt$Genus, Species=relative_bmt$Species), FUN=mean)
species <- arrange(species, desc(x))
topspecies <- unique(species$Species[1:75])
p.bmt.t <- prune_taxa(topspecies, p.bmt) #phyloseq with only 75 most abundant species
otu.t <- data.frame(otu_table(p.bmt.t))
tax.t <- data.frame(tax_table(p.bmt.t))
bmt.otu <- data.frame(otu_table(p.bmt))
meta.t <- data.frame(sample_data(p.bmt.t))

#Most abundant specieS
relative.bmt <- relative_bmt 
genera <- aggregate(relative.bmt$Abundance, by=list(Genus = relative.bmt$Genus), FUN = mean)
topgenera <- genera %>% filter(!str_detect(Genus, "unclassified")) %>% arrange(desc(x)) %>% slice_max(n=19, order_by = x)
sum(topgenera$x)
relative.bmt <- relative.bmt %>% mutate(plot.Genera = if_else(Genus %in% topgenera$Genus, Genus, "Other"))
table(relative.bmt$plot.Genera)

relative.bmt$plot.Genera <- factor(relative.bmt$plot.Genera, levels = c("Akkermansia", "Alistipes", "Bacteroides", "Bifidobacterium", "Blautia", "Clostridioides", "Eggerthella", "Enterobacter", "Enterococcus", "Erysipelatoclostridium", "Escherichia", "Faecalibacterium", "Flavonifractor", "Hungatella", "Intestinibacter", "Klebsiella", "Parabacteroides", "Ruthenibacterium", "Sellimonas", "Other"))
relative.bmt$week <- factor(relative.bmt$week, levels = c(-2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14))
table(relative.bmt$week)
relative.bmt$plot.day[relative.bmt$day >= -34 & relative.bmt$day <= -25] <- -30
relative.bmt$plot.day[relative.bmt$day >= -24 & relative.bmt$day <= -15] <- -20
relative.bmt$plot.day[relative.bmt$day >= -14 & relative.bmt$day <= -5] <- -10
relative.bmt$plot.day[relative.bmt$day >= -4 & relative.bmt$day <= 5] <- 0
relative.bmt$plot.day[relative.bmt$day >= 6 & relative.bmt$day <= 15] <- 10
relative.bmt$plot.day[relative.bmt$day >= 16 & relative.bmt$day <= 25] <- 20
relative.bmt$plot.day[relative.bmt$day >= 26 & relative.bmt$day <= 35] <- 30
relative.bmt$plot.day[relative.bmt$day >= 36 & relative.bmt$day <= 45] <- 40
relative.bmt$plot.day[relative.bmt$day >= 46 & relative.bmt$day <= 55] <- 50
relative.bmt$plot.day[relative.bmt$day >= 56 & relative.bmt$day <= 65] <- 60
relative.bmt$plot.day[relative.bmt$day >= 66 & relative.bmt$day <= 75] <- 70
relative.bmt$plot.day[relative.bmt$day >= 76 & relative.bmt$day <= 85] <- 80
relative.bmt$plot.day[relative.bmt$day >= 86 & relative.bmt$day <= 95] <- 90
relative.bmt$plot.day[relative.bmt$day >= 96 & relative.bmt$day <= 105] <- 100

#Antibiotic effect on microbiome number of species
meta_bmt <- rownames_to_column(meta_bmt, var="SampleID")
micro.num <- meta_bmt %>% dplyr::select(SampleID, num.species, new.sp, inst.score) %>% dplyr::rename("inst.score.sp" = "inst.score")
test10 <- test10 %>% mutate(SampleID = paste0("X", SampleID))
test8 <- left_join(test10, micro.num, by = "SampleID")

#Data prep for individual antibiotic effects on individual ARG classes
#beta lactam
phy.betalactam <- subset_taxa(phy.card, betalactam==1)
taxa.betalactam <- data.frame(tax_table(phy.betalactam))
nrow(taxa.betalactam)
summary(sample_sums(phy.betalactam))
otu.betalactam <- data.frame(otu_table(phy.betalactam))
for(i in 1: ncol(otu.betalactam)){otu.betalactam[,i] <- ifelse(otu.betalactam[,i] >0, 1, 0)}
tot.num.betalactam <- data.frame(colSums(otu.betalactam))
tot.num.betalactam <- tot.num.betalactam %>% dplyr::rename("num.betalactam.ARG" = "colSums.otu.betalactam.") %>% rownames_to_column(var = "sample") %>% separate("sample", into = c("X", "SampleID"), sep = 1, extra = "merge") %>% dplyr::select(SampleID, num.betalactam.ARG)
summary(tot.num.betalactam$num.betalactam.ARG)

#vancomycin
phy.vanc <- subset_taxa(phy.card, vanc==1)
taxa.vanc <- data.frame(tax_table(phy.vanc))
nrow(taxa.vanc)
summary(sample_sums(phy.vanc))
otu.vanc <- data.frame(otu_table(phy.vanc))
for(i in 1: ncol(otu.vanc)){otu.vanc[,i] <- ifelse(otu.vanc[,i] >0, 1, 0)}
tot.num.vanc <- data.frame(colSums(otu.vanc))
tot.num.vanc <- tot.num.vanc %>% dplyr::rename("num.vanc.ARG" = "colSums.otu.vanc.") %>% rownames_to_column(var = "sample") %>% separate("sample", into = c("X", "SampleID"), sep = 1, extra = "merge") %>% dplyr::select(SampleID, num.vanc.ARG)
summary(tot.num.vanc$num.vanc.ARG)

#carbapenem
phy.carbapenem <- subset_taxa(phy.card, carbapenem==1)
taxa.carbapenem <- data.frame(tax_table(phy.carbapenem))
nrow(taxa.carbapenem)
summary(sample_sums(phy.carbapenem))
otu.carbapenem <- data.frame(otu_table(phy.carbapenem))
for(i in 1: ncol(otu.carbapenem)){otu.carbapenem[,i] <- ifelse(otu.carbapenem[,i] >0, 1, 0)}
tot.num.carbapenem <- data.frame(colSums(otu.carbapenem))
tot.num.carbapenem <- tot.num.carbapenem %>% dplyr::rename("num.carbapenem.ARG" = "colSums.otu.carbapenem.") %>% rownames_to_column(var = "sample") %>% separate("sample", into = c("X", "SampleID"), sep = 1, extra = "merge") %>% dplyr::select(SampleID, num.carbapenem.ARG)
summary(tot.num.carbapenem$num.carbapenem.ARG)

#quinolone
phy.quinolone <- subset_taxa(phy.card, fluoro==1)
taxa.quinolone <- data.frame(tax_table(phy.quinolone))
nrow(taxa.quinolone)
summary(sample_sums(phy.quinolone))
otu.quinolone <- data.frame(otu_table(phy.quinolone))
for(i in 1: ncol(otu.quinolone)){otu.quinolone[,i] <- ifelse(otu.quinolone[,i] >0, 1, 0)}
tot.num.quinolone <- data.frame(colSums(otu.quinolone))
tot.num.quinolone <- tot.num.quinolone %>% dplyr::rename("num.quinolone.ARG" = "colSums.otu.quinolone.") %>% rownames_to_column(var = "sample") %>% separate("sample", into = c("X", "SampleID"), sep = 1, extra = "merge") %>% dplyr::select(SampleID, num.quinolone.ARG)
summary(tot.num.quinolone$num.quinolone.ARG)

#amino
phy.amino <- subset_taxa(phy.card, amino==1)
taxa.amino <- data.frame(tax_table(phy.amino))
nrow(taxa.amino)
summary(sample_sums(phy.amino))
otu.amino <- data.frame(otu_table(phy.amino))
for(i in 1: ncol(otu.amino)){otu.amino[,i] <- ifelse(otu.amino[,i] >0, 1, 0)}
tot.num.amino <- data.frame(colSums(otu.amino))
tot.num.amino <- tot.num.amino %>% dplyr::rename("num.amino.ARG" = "colSums.otu.amino.") %>% rownames_to_column(var = "sample") %>% separate("sample", into = c("X", "SampleID"), sep = 1, extra = "merge") %>% dplyr::select(SampleID, num.amino.ARG)
summary(tot.num.amino$num.amino.ARG)

#clinda
phy.clinda <- subset_taxa(phy.card, clinda==1)
taxa.clinda <- data.frame(tax_table(phy.clinda))
nrow(taxa.clinda)
summary(sample_sums(phy.clinda))
otu.clinda <- data.frame(otu_table(phy.clinda))
for(i in 1: ncol(otu.clinda)){otu.clinda[,i] <- ifelse(otu.clinda[,i] >0, 1, 0)}
tot.num.clinda <- data.frame(colSums(otu.clinda))
tot.num.clinda <- tot.num.clinda %>% dplyr::rename("num.clinda.ARG" = "colSums.otu.clinda.") %>% rownames_to_column(var = "sample") %>% separate("sample", into = c("X", "SampleID"), sep = 1, extra = "merge") %>% dplyr::select(SampleID, num.clinda.ARG)
summary(tot.num.clinda$num.clinda.ARG)

#macrolide
phy.macrolide <- subset_taxa(phy.card, macrolide==1)
taxa.macrolide <- data.frame(tax_table(phy.macrolide))
nrow(taxa.macrolide)
summary(sample_sums(phy.macrolide))
otu.macrolide <- data.frame(otu_table(phy.macrolide))
for(i in 1: ncol(otu.macrolide)){otu.macrolide[,i] <- ifelse(otu.macrolide[,i] >0, 1, 0)}
tot.num.macrolide <- data.frame(colSums(otu.macrolide))
tot.num.macrolide <- tot.num.macrolide %>% dplyr::rename("num.macrolide.ARG" = "colSums.otu.macrolide.") %>% rownames_to_column(var = "sample") %>% separate("sample", into = c("X", "SampleID"), sep = 1, extra = "merge") %>% dplyr::select(SampleID, num.macrolide.ARG)
summary(tot.num.macrolide$num.macrolide.ARG)

#trimeth
phy.trimeth <- subset_taxa(phy.card, trimeth==1)
taxa.trimeth <- data.frame(tax_table(phy.trimeth))
nrow(taxa.trimeth)
summary(sample_sums(phy.trimeth))
otu.trimeth <- data.frame(otu_table(phy.trimeth))
for(i in 1: ncol(otu.trimeth)){otu.trimeth[,i] <- ifelse(otu.trimeth[,i] >0, 1, 0)}
tot.num.trimeth <- data.frame(colSums(otu.trimeth))
tot.num.trimeth <- tot.num.trimeth %>% dplyr::rename("num.trimeth.ARG" = "colSums.otu.trimeth.") %>% rownames_to_column(var = "sample") %>% separate("sample", into = c("X", "SampleID"), sep = 1, extra = "merge") %>% dplyr::select(SampleID, num.trimeth.ARG)
summary(tot.num.trimeth$num.trimeth.ARG)

#tetra
phy.tetra <- subset_taxa(phy.card, tetra==1)
taxa.tetra <- data.frame(tax_table(phy.tetra))
nrow(taxa.tetra)
summary(sample_sums(phy.tetra))
otu.tetra <- data.frame(otu_table(phy.tetra))
for(i in 1: ncol(otu.tetra)){otu.tetra[,i] <- ifelse(otu.tetra[,i] >0, 1, 0)}
tot.num.tetra <- data.frame(colSums(otu.tetra))
tot.num.tetra <- tot.num.tetra %>% dplyr::rename("num.tetra.ARG" = "colSums.otu.tetra.") %>% rownames_to_column(var = "sample") %>% separate("sample", into = c("X", "SampleID"), sep = 1, extra = "merge") %>% dplyr::select(SampleID, num.tetra.ARG)
summary(tot.num.tetra$num.tetra.ARG)

#flagyl
phy.flagyl <- subset_taxa(phy.card, flagyl==1)
taxa.flagyl <- data.frame(tax_table(phy.flagyl))
nrow(taxa.flagyl)
summary(sample_sums(phy.flagyl))
otu.flagyl <- data.frame(otu_table(phy.flagyl))
for(i in 1: ncol(otu.flagyl)){otu.flagyl[,i] <- ifelse(otu.flagyl[,i] >0, 1, 0)}
tot.num.flagyl <- data.frame(colSums(otu.flagyl))
tot.num.flagyl <- tot.num.flagyl %>% dplyr::rename("num.flagyl.ARG" = "colSums.otu.flagyl.") %>% rownames_to_column(var = "sample") %>% separate("sample", into = c("X", "SampleID"), sep = 1, extra = "merge") %>% dplyr::select(SampleID, num.flagyl.ARG)
summary(tot.num.flagyl$num.flagyl.ARG)

#phenicol
phy.phenicol <- subset_taxa(phy.card, phenicol==1)
taxa.phenicol <- data.frame(tax_table(phy.phenicol))
nrow(taxa.phenicol)
summary(sample_sums(phy.phenicol))
otu.phenicol <- data.frame(otu_table(phy.phenicol))
for(i in 1: ncol(otu.phenicol)){otu.phenicol[,i] <- ifelse(otu.phenicol[,i] >0, 1, 0)}
tot.num.phenicol <- data.frame(colSums(otu.phenicol))
tot.num.phenicol <- tot.num.phenicol %>% dplyr::rename("num.phenicol.ARG" = "colSums.otu.phenicol.") %>% rownames_to_column(var = "sample") %>% separate("sample", into = c("X", "SampleID"), sep = 1, extra = "merge") %>% dplyr::select(SampleID, num.phenicol.ARG)
summary(tot.num.phenicol$num.phenicol.ARG)

#linezolid
phy.linezolid <- subset_taxa(phy.card, linezolid==1)
taxa.linezolid <- data.frame(tax_table(phy.linezolid))
nrow(taxa.linezolid)
summary(sample_sums(phy.linezolid))
otu.linezolid <- data.frame(otu_table(phy.linezolid))
for(i in 1: ncol(otu.linezolid)){otu.linezolid[,i] <- ifelse(otu.linezolid[,i] >0, 1, 0)}
tot.num.linezolid <- data.frame(colSums(otu.linezolid))
tot.num.linezolid <- tot.num.linezolid %>% dplyr::rename("num.linezolid.ARG" = "colSums.otu.linezolid.") %>% rownames_to_column(var = "sample") %>% separate("sample", into = c("X", "SampleID"), sep = 1, extra = "merge") %>% dplyr::select(SampleID, num.linezolid.ARG)
summary(tot.num.linezolid$num.linezolid.ARG)

#sulfa
phy.sulfa <- subset_taxa(phy.card, sulfa==1)
taxa.sulfa <- data.frame(tax_table(phy.sulfa))
nrow(taxa.sulfa)
summary(sample_sums(phy.sulfa))
otu.sulfa <- data.frame(otu_table(phy.sulfa))
for(i in 1: ncol(otu.sulfa)){otu.sulfa[,i] <- ifelse(otu.sulfa[,i] >0, 1, 0)}
tot.num.sulfa <- data.frame(colSums(otu.sulfa))
tot.num.sulfa <- tot.num.sulfa %>% dplyr::rename("num.sulfa.ARG" = "colSums.otu.sulfa.") %>% rownames_to_column(var = "sample") %>% separate("sample", into = c("X", "SampleID"), sep = 1, extra = "merge") %>% dplyr::select(SampleID, num.sulfa.ARG)
summary(tot.num.sulfa$num.sulfa.ARG)

class.num <- left_join(tot.num.betalactam, tot.num.vanc, by = "SampleID")
class.num <- left_join(class.num, tot.num.carbapenem, by = "SampleID")
class.num <- left_join(class.num, tot.num.quinolone, by = "SampleID")
class.num <- left_join(class.num, tot.num.amino, by = "SampleID")
class.num <- left_join(class.num, tot.num.clinda, by = "SampleID")
class.num <- left_join(class.num, tot.num.macrolide, by = "SampleID")
class.num <- left_join(class.num, tot.num.trimeth, by = "SampleID")
class.num <- left_join(class.num, tot.num.tetra, by = "SampleID")
class.num <- left_join(class.num, tot.num.flagyl, by = "SampleID")
test9 <- left_join(test8, class.num, by = "SampleID")
test9 <- test9 %>% mutate(other.abx.exp = if_else(tetracycline.exp == 1, 1, if_else(clinda.exp == 1, 1, if_else(ag.exp == 1, 1, if_else(other_beta_lactam.exp == 1, 1, 0)))))
sum(test9$other.abx.exp)


#####ANALYSES AND FIGURES#####
#MICROBIOME
#Figure 1
colors_20  <-  c("gold1", "#E31A1C", "#6A3D9A", "dodgerblue2", "gray50",
                 "palegreen2", "green4", "#FF7F00", "skyblue2", "#CAB2D6", 
                 "#FDBF6F", # lt orange
                 "maroon","deeppink1","blue1","steelblue4", 
                 "darkturquoise", "darkorange4","brown", "#FB9A99", "black")

plot_days <- ggplot(arrange(relative.bmt, plot.Genera), aes(x=plot.day, y=Abundance, fill=plot.Genera)) +
  geom_bar(stat="identity", position="fill") +
  guides(fill = guide_legend(ncol=1, byrow=TRUE, title="Genera")) + 
  scale_fill_manual(values=colors_20) + 
  scale_y_continuous(limits = c(0, 1), expand = expansion(add =c(0,0.02)))+
  scale_x_continuous(limits = c(-35, 105), breaks = c(-30, -20, -10, 0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100), expand = c(0,0))+
  theme_classic()+
  xlab("Day Relative to HCT") + ylab("Relative Abundance")+
  theme(legend.text=element_text(size=24, face = "italic"),
        legend.title=element_text(size=24), 
        axis.title.y = element_text(size=24), 
        axis.text = element_text(size=24, colour="black"), 
        axis.title.x = element_text(size=24), 
        plot.title = element_text(size=24, hjust=1.4), 
        strip.text.x = element_text(size = 24), strip.background=element_rect(fill="white")) 
plot_days

#MaAsLin
library(Maaslin2)
fit_data = Maaslin2(
  input_data = otu_table(p.bmt.t),
  input_metadata = meta_bmt,
  output = "All",
  fixed_effects = c("cefepime.exp", "piptazo.exp", "carbapenem.exp", "vanc.exp", "quinolone.exp", "tmp_smx.exp", "macrolide.exp", "flagyl.exp", "other.exp", "age", "sex", "diagnosis", "prep", "type_hsct", "day", "log.depth"),
  random_effects = c("study_id"),
  normalization = "NONE",
  min_abundance = 0.01,
  min_prevalence = 0.05,
  max_significance = 0.10)

#Figure 2
significant.species <- read.csv("T:/Projects/PID/Kelly/Intestinal Microbiome of  - Pro00064365/Resistome Analysis/Analytic Data/Microbiome species and antibiotics maaslin significant_results 05182022.csv") 
significant.species <- significant.species %>% mutate(abs.effect = abs(coef)) %>% arrange(desc(abs.effect)) 
table(significant.species$metadata)
table(significant.species$value)
table(significant.species$feature)
just.sig <- significant.species %>% mutate(Effect = coef, Bug = feature, Covariate = if_else(metadata == "age", "Age", if_else(metadata=="cefepime.exp", "Cefepime", if_else(metadata == "day", "Day Relative to HCT", if_else(metadata == "diagnosis" & value == "Congenital_immunodef", "Immunodeficiency",
                                                                                                                                                                                                                               if_else(metadata =="diagnosis" & value == "Heme_non-malignancy", "Non-malignant Heme", if_else(metadata=="diagnosis" & value=="Metabolic_disease", "Metabolic Disease", if_else(metadata=="diagnosis" & value == "Solid_tumor", "Solid Tumor", 
                                                                                                                                                                                                                                                                                                                                                                                                               if_else(metadata == "flagyl.exp", "Metronidazole", if_else(metadata == "log.depth", "Sequencing Depth", if_else(metadata == "macrolide.exp", "Macrolide",
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               if_else(metadata == "other.exp", "Other Antibiotics", if_else(metadata== "piptazo.exp", "Pip-Tazo", if_else(metadata=="quinolone.exp", "Fluoroquinolone", if_else(metadata == "tmp_smx.exp", "TMP-SMX", if_else(metadata=="type_hsct", "Allogeneic", if_else(metadata=="vanc.exp", "Vancomycin", "check")))))))))))))))))
just.sig$Bug[just.sig$feature == "Alistipes_finegoldii"] <- "Alistipes finegoldii"
just.sig$Bug[just.sig$feature == "Alistipes_putredinis"] <- "Alistipes putredinis"
just.sig$Bug[just.sig$feature == "Bacteroides_dorei"] <- "Bacteroides dorei"
just.sig$Bug[just.sig$feature == "Bacteroides_fragilis"] <- "Bacteroides fragilis"
just.sig$Bug[just.sig$feature == "Bacteroides_ovatus"] <- "Bacteroides ovatus"
just.sig$Bug[just.sig$feature == "Bacteroides_thetaiotaomicron"] <- "Bacteroides thetaiotaomicron"
just.sig$Bug[just.sig$feature == "Bacteroides_uniformis"] <- "Bacteroides uniformis"
just.sig$Bug[just.sig$feature == "Bacteroides_vulgatus"] <- "Bacteroides vulgatus"
just.sig$Bug[just.sig$feature == "Bacteroides_xylanisolvens"] <- "Bacteroides xylanisolvens"
just.sig$Bug[just.sig$feature == "Bifidobacterium_breve"] <- "Bifidobacterium breve"
just.sig$Bug[just.sig$feature == "Bifidobacterium_longum"] <- "Bifidobacterium longum"
just.sig$Bug[just.sig$feature == "Blautia_coccoides"] <- "Blautia coccoides"
just.sig$Bug[just.sig$feature == "Blautia_wexlerae"] <- "Blautia wexlerae"
just.sig$Bug[just.sig$feature == "Clostridioides_difficile"] <- "Clostridioides difficile"
just.sig$Bug[just.sig$feature == "Clostridium_clostridioforme"] <- "Clostridium clostridioforme"
just.sig$Bug[just.sig$feature == "Clostridium_innocuum"] <- "Clostridium innocuum"
just.sig$Bug[just.sig$feature == "Clostridium_spiroforme"] <- "Clostridium spiroforme"
just.sig$Bug[just.sig$feature == "Clostridium_symbiosum"] <- "Clostridium symbiosum"
just.sig$Bug[just.sig$feature == "Eggerthella_lenta"] <- "Eggerthella lenta"
just.sig$Bug[just.sig$feature == "Enterococcus_faecalis"] <- "Enterococcus faecalis"
just.sig$Bug[just.sig$feature == "Enterococcus_faecium"] <- "Enterococcus faecium"
just.sig$Bug[just.sig$feature == "Enterococcus_gallinarum"] <- "Enterococcus gallinarum"
just.sig$Bug[just.sig$feature == "Erysipelatoclostridium_ramosum"] <- "Erysipelatoclostridium ramosum"
just.sig$Bug[just.sig$feature == "Escherichia_coli"] <- "Escherichia coli"
just.sig$Bug[just.sig$feature == "Eubacterium_rectale"] <- "Eubacterium rectale"
just.sig$Bug[just.sig$feature == "Eubacterium_sp_CAG_180"] <- "Eubacterium sp CAG 180"
just.sig$Bug[just.sig$feature == "Faecalibacterium_prausnitzii"] <- "Faecalibacterium prausnitzii"
just.sig$Bug[just.sig$feature == "Flavonifractor_plautii"] <- "Flavonifractor plautii"
just.sig$Bug[just.sig$feature == "Fusicatenibacter_saccharivorans"] <- "Fusicatenibacter saccharivorans"
just.sig$Bug[just.sig$feature == "Hungatella_hathewayi"] <- "Hungatella hathewayi"
just.sig$Bug[just.sig$feature == "Intestinibacter_bartlettii"] <- "Intestinibacter bartlettii"
just.sig$Bug[just.sig$feature == "Klebsiella_pneumoniae"] <- "Klebsiella pneumoniae"
just.sig$Bug[just.sig$feature == "Klebsiella_variicola"] <- "Klebsiella variicola"
just.sig$Bug[just.sig$feature == "Lactobacillus_rhamnosus"] <- "Lactobacillus rhamnosus"
just.sig$Bug[just.sig$feature == "Parabacteroides_distasonis"] <- "Parabacteroides distasonis"
just.sig$Bug[just.sig$feature == "Parabacteroides_merdae"] <- "Parabacteroides merdae"
just.sig$Bug[just.sig$feature == "Roseburia_faecis"] <- "Roseburia faecis"
just.sig$Bug[just.sig$feature == "Roseburia_intestinalis"] <- "Roseburia intestinalis"
just.sig$Bug[just.sig$feature == "Ruminococcus_gnavus"] <- "Ruminococcus gnavus"
just.sig$Bug[just.sig$feature == "Sellimonas_intestinalis"] <- "Sellimonas intestinalis"
just.sig$Bug[just.sig$feature == "Streptococcus_parasanguinis"] <- "Streptococcus parasanguinis"
just.sig$Bug[just.sig$feature == "Streptococcus_salivarius"] <- "Streptococcus salivarius"
just.sig$Bug[just.sig$feature == "Streptococcus_thermophilus"] <- "Streptococcus thermophilus"
just.sig$Bug[just.sig$feature == "Veillonella_parvula"] <- "Veillonella parvula"

supp.table <- just.sig %>% mutate(Species = Bug, Standard_Error = stderr, P_Value = pval, Q_Value = qval) %>%
  dplyr::select(Species, Covariate, Effect, Standard_Error, P_Value, Q_Value)
table(supp.table$Covariate)
supp.table$Covariate[supp.table$Covariate=="Pip-Tazo"] <- "Piperacillin-Tazobactam"
supp.table$Covariate[supp.table$Covariate=="Non-malignant Heme"] <- "Hematologic Non-Malignant Disease"
supp.table$Covariate[supp.table$Covariate=="TMP-SMX"] <- "Trimethoprim-Sulfamethoxazole"
supp.table <- supp.table %>% arrange(Covariate) %>% arrange(Species)
write.csv(supp.table, "T:/Projects/PID/Kelly/Intestinal Microbiome of  - Pro00064365/Resistome Analysis/Analytic Data/sig maaslin output for supplemental table 05182022.csv")

just.sig <- just.sig %>% dplyr::select(Bug, Covariate, Effect)

#rescaling IRR so that colors are more different
rescale.plot.neg <- just.sig %>% filter(Effect < 0)
#rescale.plot.neg$Effect[rescale.plot.neg$Effect < -1.5] <- -1.5
normalize <- function(x, na.rm = TRUE) {
  return((x- min(x)) /(max(x)-min(x)))
}
rescale.plot.neg$rescale.Effect = normalize(rescale.plot.neg$Effect)
rescale.plot.neg <- rescale.plot.neg %>% mutate(plot.Effect = rescale.Effect-1) 
summary(rescale.plot.neg$rescale.Effect)
filter(rescale.plot.neg, rescale.Effect == median(rescale.plot.neg$rescale.Effect))
rescale.plot.neg <- rescale.plot.neg %>% dplyr::select(-rescale.Effect) 
rescale.plot.pos <- just.sig %>% filter(Effect >= 0) 
#rescale.plot.pos$Effect[rescale.plot.pos$Effect >1.5] <- 1.5
hist(rescale.plot.pos$Effect)
rescale.plot.pos$plot.Effect = normalize(rescale.plot.pos$Effect)
summary(rescale.plot.pos$Effect)
filter(rescale.plot.pos, plot.Effect == median(rescale.plot.pos$plot.Effect))

rescale.plot <- rbind(rescale.plot.neg, rescale.plot.pos)
rescale.plot <- rescale.plot %>% mutate(sig = 1, common = paste0(Bug, Covariate))
significant <- rescale.plot %>% dplyr::select(common, sig)

all.result <- read.csv("T:/Projects/PID/Kelly/Intestinal Microbiome of  - Pro00064365/Resistome Analysis/Analytic Data/Microbiome species and antibiotics maaslin all results 05182022.csv") #Supplementary Table 1
all.result <- all.result %>% mutate(Effect = coef)
table(all.result$feature)
all.result$Bug[all.result$feature == "Alistipes_finegoldii"] <- "Alistipes finegoldii"
all.result$Bug[all.result$feature == "Alistipes_putredinis"] <- "Alistipes putredinis"
all.result$Bug[all.result$feature == "Bacteroides_caccae"] <- "Bacteroides caccae"
all.result$Bug[all.result$feature == "Bacteroides_dorei"] <- "Bacteroides dorei"
all.result$Bug[all.result$feature == "Bacteroides_fragilis"] <- "Bacteroides fragilis"
all.result$Bug[all.result$feature == "Bacteroides_ovatus"] <- "Bacteroides ovatus"
all.result$Bug[all.result$feature == "Bacteroides_stercoris"] <- "Bacteroides stercoris"
all.result$Bug[all.result$feature == "Bacteroides_thetaiotaomicron"] <- "Bacteroides thetaiotaomicron"
all.result$Bug[all.result$feature == "Bacteroides_uniformis"] <- "Bacteroides uniformis"
all.result$Bug[all.result$feature == "Bacteroides_vulgatus"] <- "Bacteroides vulgatus"
all.result$Bug[all.result$feature == "Bacteroides_xylanisolvens"] <- "Bacteroides xylanisolvens"
all.result$Bug[all.result$feature == "Bifidobacterium_breve"] <- "Bifidobacterium breve"
all.result$Bug[all.result$feature == "Bifidobacterium_longum"] <- "Bifidobacterium longum"
all.result$Bug[all.result$feature == "Blautia_coccoides"] <- "Blautia coccoides"
all.result$Bug[all.result$feature == "Blautia_wexlerae"] <- "Blautia wexlerae"
all.result$Bug[all.result$feature == "Clostridioides_difficile"] <- "Clostridioides difficile"
all.result$Bug[all.result$feature == "Clostridium_bolteae"] <- "Clostridium bolteae"
all.result$Bug[all.result$feature == "Clostridium_clostridioforme"] <- "Clostridium clostridioforme"
all.result$Bug[all.result$feature == "Clostridium_innocuum"] <- "Clostridium innocuum"
all.result$Bug[all.result$feature == "Clostridium_spiroforme"] <- "Clostridium spiroforme"
all.result$Bug[all.result$feature == "Clostridium_symbiosum"] <- "Clostridium symbiosum"
all.result$Bug[all.result$feature == "Eggerthella_lenta"] <- "Eggerthella lenta"
all.result$Bug[all.result$feature == "Enterococcus_faecalis"] <- "Enterococcus faecalis"
all.result$Bug[all.result$feature == "Enterococcus_faecium"] <- "Enterococcus faecium"
all.result$Bug[all.result$feature == "Enterococcus_gallinarum"] <- "Enterococcus gallinarum"
all.result$Bug[all.result$feature == "Erysipelatoclostridium_ramosum"] <- "Erysipelatoclostridium ramosum"
all.result$Bug[all.result$feature == "Escherichia_coli"] <- "Escherichia coli"
all.result$Bug[all.result$feature == "Eubacterium_rectale"] <- "Eubacterium rectale"
all.result$Bug[all.result$feature == "Eubacterium_sp_CAG_180"] <- "Eubacterium sp CAG 180"
all.result$Bug[all.result$feature == "Faecalibacterium_prausnitzii"] <- "Faecalibacterium prausnitzii"
all.result$Bug[all.result$feature == "Flavonifractor_plautii"] <- "Flavonifractor plautii"
all.result$Bug[all.result$feature == "Fusicatenibacter_saccharivorans"] <- "Fusicatenibacter saccharivorans"
all.result$Bug[all.result$feature == "Hungatella_hathewayi"] <- "Hungatella hathewayi"
all.result$Bug[all.result$feature == "Intestinibacter_bartlettii"] <- "Intestinibacter bartlettii"
all.result$Bug[all.result$feature == "Klebsiella_pneumoniae"] <- "Klebsiella pneumoniae"
all.result$Bug[all.result$feature == "Klebsiella_variicola"] <- "Klebsiella variicola"
all.result$Bug[all.result$feature == "Lactobacillus_rhamnosus"] <- "Lactobacillus rhamnosus"
all.result$Bug[all.result$feature == "Parabacteroides_distasonis"] <- "Parabacteroides distasonis"
all.result$Bug[all.result$feature == "Parabacteroides_merdae"] <- "Parabacteroides merdae"
all.result$Bug[all.result$feature == "Roseburia_faecis"] <- "Roseburia faecis"
all.result$Bug[all.result$feature == "Roseburia_intestinalis"] <- "Roseburia intestinalis"
all.result$Bug[all.result$feature == "Ruminococcus_gnavus"] <- "Ruminococcus gnavus"
all.result$Bug[all.result$feature == "Ruthenibacterium_lactatiformans"] <- "Ruthenibacterium lactatiformans"
all.result$Bug[all.result$feature == "Sellimonas_intestinalis"] <- "Sellimonas intestinalis"
all.result$Bug[all.result$feature == "Streptococcus_parasanguinis"] <- "Streptococcus parasanguinis"
all.result$Bug[all.result$feature == "Streptococcus_salivarius"] <- "Streptococcus salivarius"
all.result$Bug[all.result$feature == "Streptococcus_thermophilus"] <- "Streptococcus thermophilus"
all.result$Bug[all.result$feature == "Veillonella_parvula"] <- "Veillonella parvula"

all.result$Covariate[all.result$metadata=="age"] <- "Age"
all.result$Covariate[all.result$metadata=="carbapenem.exp"] <- "Carbapenem"
all.result$Covariate[all.result$metadata=="cefepime.exp"] <- "Cefepime"
all.result$Covariate[all.result$metadata=="day"] <- "Day Relative to HCT"
all.result$Covariate[all.result$metadata=="diagnosis" & all.result$value == "Congenital_immunodef"] <- "Immunodeficiency"
all.result$Covariate[all.result$metadata=="flagyl.exp"] <- "Metronidazole"
all.result$Covariate[all.result$metadata=="log.depth"] <- "Sequencing Depth"
all.result$Covariate[all.result$metadata=="macrolide.exp"] <- "Macrolide"
all.result$Covariate[all.result$metadata=="other.exp"] <- "Other Antibiotics"
all.result$Covariate[all.result$metadata=="piptazo.exp"] <- "Pip-Tazo"
all.result$Covariate[all.result$metadata=="prep"] <- "Non-myelo/RIC"
all.result$Covariate[all.result$metadata=="quinolone.exp"] <- "Fluoroquinolone"
all.result$Covariate[all.result$metadata=="sex"] <- "Male"
all.result$Covariate[all.result$metadata=="tmp_smx.exp"] <- "TMP-SMX"
all.result$Covariate[all.result$metadata=="type_hsct"] <- "Allogeneic"
all.result$Covariate[all.result$metadata=="vanc.exp"] <- "Vancomycin"
all.result$Covariate[all.result$metadata=="diagnosis" & all.result$value =="Heme_non-malignancy"] <- "Non-malignant Heme"
all.result$Covariate[all.result$metadata=="diagnosis" & all.result$value =="Metabolic_disease"] <- "Metabolic Disease"
all.result$Covariate[all.result$metadata=="diagnosis" & all.result$value =="Solid_tumor"] <- "Solid Tumor"

all.rescale.plot.neg <- all.result %>% filter(Effect < 0)
all.rescale.plot.neg$rescale.Effect = normalize(all.rescale.plot.neg$Effect)
all.rescale.plot.neg <- all.rescale.plot.neg %>% mutate(plot.Effect = rescale.Effect-1) 
summary(all.rescale.plot.neg$rescale.Effect)
filter(all.rescale.plot.neg, rescale.Effect == median(all.rescale.plot.neg$rescale.Effect))
all.rescale.plot.neg <- all.rescale.plot.neg %>% dplyr::select(-rescale.Effect) 
all.rescale.plot.pos <- all.result %>% filter(Effect >= 0) 
all.rescale.plot.pos$plot.Effect = normalize(all.rescale.plot.pos$Effect)
summary(all.rescale.plot.pos$plot.Effect)

all.rescale.plot <- rbind(all.rescale.plot.neg, all.rescale.plot.pos)
all.rescale.plot$common = paste0(all.rescale.plot$Bug, all.rescale.plot$Covariate)
all.rescale.plot <- left_join(all.rescale.plot, significant, by = "common")
table(all.rescale.plot$sig)
all.rescale.plot$sig <- as.factor(all.rescale.plot$sig)
#Figure 2
heat.species.all <- ggplot(all.rescale.plot, aes(x =Covariate, y= Bug, fill = plot.Effect, color = ""))+
  geom_tile(color = "black", size = 0.5) +
  geom_point(aes(size = sig),
             shape = 8, na.rm = TRUE, stroke= 0.5,
             show.legend = FALSE,
             color="black")+
  scale_size_manual(values=c(2))+
  scale_fill_gradient2(low = "navy",  high = "orangered2", midpoint = 0,
                       breaks = c(-1, 0, 1),
                       labels = c(expression(""<=-1.5), "1.0", expression("">=1.5)),
                       name = expression(beta),
                       na.value = "gray 70")+
  labs(y = "", x = "")+
  scale_x_discrete(
    limits = c("Cefepime", "Pip-Tazo", "Carbapenem", "Vancomycin", "Metronidazole", "Fluoroquinolone", "TMP-SMX", "Macrolide","Other Antibiotics", "Age", "Male", "Non-malignant Heme", "Metabolic Disease", "Immunodeficiency", "Solid Tumor", "Non-myelo/RIC", "Allogeneic", "Day Relative to HCT", "Sequencing Depth"),
    position = "top")+
  scale_y_discrete(
    limits = rev,
    position = "left")+ 
  theme_minimal()+
  theme(axis.text.x.top = element_text(angle = 90, hjust = 0, vjust = 0.5, color = "black"),# size = 26),
        axis.title.x.top = element_text(vjust = 2, color = "black"),#, size = 30),
        axis.title.y = element_text(vjust = 3),# size = 30),
        axis.text.y = element_text(color = "black", face = "italic"),# size = 26),
        axis.ticks = element_blank(),
        plot.margin = margin(t = 2.5,
                             b=0,
                             l=0, r=0, unit = "cm"))+
  coord_fixed()
heat.species.all

#Antibiotics effects on microbiome composition (Supplementary Table 2)
#Number of Species
numspec.anaerobe.adj <- glmmadmb(num.species ~ aerobic.exp + anaerobic.exp + age + sex + diagnosis + prep + type_hsct + day + log.depth + (1 |study_id), data  = test8, family = "nbinom1")
summary(numspec.anaerobe.adj)
numspec.anaerobe.adj.model <- cbind(Estimate = coef(numspec.anaerobe.adj), confint(numspec.anaerobe.adj))
exp(numspec.anaerobe.adj.model)

#new species since previous sample
newspec.anaerobe.adj <- glmmadmb(new.sp ~ aerobic.exp + anaerobic.exp + age + sex + diagnosis + prep + type_hsct + day + log.depth + (1 |study_id), data  = test8, family = "nbinom1")
summary(newspec.anaerobe.adj)
newspec.anaerobe.adj.model <- cbind(Estimate = coef(newspec.anaerobe.adj), confint(newspec.anaerobe.adj))
exp(newspec.anaerobe.adj.model)

#Species instability
instspec.anaerobe.adj <- lmer(inst.score.sp ~ aerobic.exp + anaerobic.exp + age + sex + diagnosis + prep + type_hsct + day + log.depth + (1 |study_id), data  = test8)
summary(instspec.anaerobe.adj)
confint(instspec.anaerobe.adj)

#RESISTOME
#Antibiotics effects on resistome (Supplementary Table 2)
#Abundance
rpkm.anaerobe.adj <- glmmadmb(total.abund ~ aerobic.exp + anaerobic.exp + age + sex + diagnosis + prep + type_hsct + day + log.depth + (1 |study_id), data  = test10, family = "nbinom2")
summary(rpkm.anaerobe.adj)
RPKM.anaerobe.adj.model <- cbind(Estimate = coef(rpkm.anaerobe.adj), confint(rpkm.anaerobe.adj))
anaerobe.rpkm <- data.frame(exp(RPKM.anaerobe.adj.model))
exp(RPKM.anaerobe.adj.model)
anaerobe.rpkm <- anaerobe.rpkm %>% rownames_to_column(var = "row") %>% filter(row == "aerobic.exp1" | row == "anaerobic.exp1") %>% mutate(effect = round(Estimate, 2), lower= round(X2.5.., 2), upper = round(X97.5.., 2)) %>%
  mutate(drug = if_else(row == "aerobic.exp1", "Aerobic Antibiotics", "Anaerobic Antibiotics")) %>% dplyr::select(drug, effect, lower, upper)

#Number of ARGs
numARG.anaerobe.adj <- glmmadmb(num.ARG ~ aerobic.exp + anaerobic.exp + age + sex + diagnosis + prep + type_hsct + day + log.depth + (1 |study_id), data  = test10, family = "nbinom1")
summary(numARG.anaerobe.adj)
numARG.anaerobe.adj.model <- cbind(Estimate = coef(numARG.anaerobe.adj), confint(numARG.anaerobe.adj))
anaerobe.num <- data.frame(exp(numARG.anaerobe.adj.model))
exp(numARG.anaerobe.adj.model)
anaerobe.num <- anaerobe.num %>% rownames_to_column(var = "row") %>% filter(row == "aerobic.exp1" | row == "anaerobic.exp1") %>% mutate(effect = round(Estimate, 2), lower= round(X2.5.., 2), upper = round(X97.5.., 2)) %>%
  mutate(drug = if_else(row == "aerobic.exp1", "Aerobic Antibiotics", "Anaerobic Antibiotics")) %>% dplyr::select(drug, effect, lower, upper)

#Number of new ARGs
newARG.anaerobe.adj <- glmmadmb(new.ARG ~ aerobic.exp + anaerobic.exp + age + sex + diagnosis + prep + type_hsct + day + log.depth + (1 |study_id), data  = test10, family = "nbinom1")
summary(newARG.anaerobe.adj)
newARG.anaerobe.adj.model <- cbind(Estimate = coef(newARG.anaerobe.adj), confint(newARG.anaerobe.adj))
anaerobe.new <- data.frame(exp(newARG.anaerobe.adj.model))
exp(newARG.anaerobe.adj.model)
anaerobe.new <- anaerobe.new %>% rownames_to_column(var = "row") %>% filter(row == "aerobic.exp1" | row == "anaerobic.exp1") %>% mutate(effect = round(Estimate, 2), lower= round(X2.5.., 2), upper = round(X97.5.., 2)) %>%
  mutate(drug = if_else(row == "aerobic.exp1", "Aerobic Antibiotics", "Anaerobic Antibiotics")) %>% dplyr::select(drug, effect, lower, upper)

#instability
inst.anaerobe.adj <- lmer(inst.score ~ aerobic.exp + anaerobic.exp + age + sex + diagnosis + prep + type_hsct + day + log.depth + (1 |study_id), data  = test10)
summary(inst.anaerobe.adj)
confint(inst.anaerobe.adj)
AIC(inst.anaerobe.adj)

inst.anaerobe.exp.adj.model <- data.frame(cbind(coef(summary(inst.anaerobe.adj))), confint(inst.anaerobe.adj, parm = c("(Intercept)", "aerobic.exp1", "anaerobic.exp1", "age" , "sex2" , "diagnosisCongenital_immunodef", "diagnosisHeme_non-malignancy", "diagnosisMetabolic_disease", "diagnosisSolid_tumor", "prep2" , "type_hsctAllogeneic" , "day" , "log.depth" ))) %>%
  rownames_to_column(var = "row") %>% filter(row == "aerobic.exp1" | row == "anaerobic.exp1")  %>%
  mutate(effect = round(Estimate, 3), lower = round(X2.5.., 3), upper = round(X97.5.., 3), drug = if_else(row == "aerobic.exp1", "Aerobic Antibiotics", "Anaerobic Antibiotics")) %>%
  dplyr::select(drug, effect, lower, upper)

#Correlation between microbiome and resistome measures
resistome <- test10 %>% dplyr::select(study_id, SampleID, day, num.ARG, new.ARG, inst.score) %>% dplyr::rename("inst.score.ARG" = "inst.score") 
meta_bmt <- meta_bmt %>% rownames_to_column(var="SampleID")
composition <- meta_bmt %>% dplyr::select(study_id, SampleID, day, num.species, new.sp, inst.score) %>% dplyr::rename("inst.score.sp" = "inst.score")
res.comp <- left_join(composition, resistome, by= "SampleID") %>% dplyr::rename("day" = "day.x", "study_id" = "study_id.x") %>% dplyr::select(-c(day.y, study_id.y))
plot.df.c <- composition %>% dplyr::select(-c(SampleID, study_id)) %>% mutate(Source = "Composition") %>% dplyr::rename("New" = "new.sp", "Num" = "num.species", "inst.score" = "inst.score.sp")
plot.df.r <- resistome %>% dplyr::select(-c(SampleID, study_id)) %>% mutate(Source = "Resistome") %>% dplyr::rename("New" = "new.ARG", "Num" = "num.ARG", "inst.score" = "inst.score.ARG")
plot.df <- rbind(plot.df.c, plot.df.r)
plot.df$Source <- as.factor(plot.df$Source)

#Figure 3
colors <- c("Resistome" = "blue3", "Composition" = "red3")

plot1.2 <- ggplot(data = plot.df, aes(x = day, y = Num, colour = Source))+
  geom_point(aes(x=day, y=Num),
             color = "indianred2",
             data = filter(plot.df, Source == "Composition")) +
  geom_point(aes(x = day, y=Num),
             color = "turquoise",
             data = filter(plot.df, Source == "Resistome"))+
  geom_smooth(aes(x=day, y=Num),
              color = "indianred4",
              fill = "rosybrown3",
              data = filter(plot.df, Source == "Composition"))+
  geom_smooth(aes(x = day, y = Num),
              color = "turquoise4",
              fill = "paleturquoise3",              
              data = filter(plot.df, Source == "Resistome"))+
  theme_classic()+
  theme(text = element_text(size = 18),
        axis.text = element_text(color = "black"),
        plot.margin = unit(c(1,0.5,0.5,1), "cm"))+
  ylab("Number of Species/ARGs")+
  xlab("Day Relative to HCT")+
  scale_y_continuous(limits = c(0, 130), breaks = c(0, 20, 40, 60, 80, 100, 130), expand = c(0,0))+
  scale_x_continuous(limits = c(-30, 101), breaks = c(-30, 0, 30, 60, 100), expand = c(0,0))
plot1.2

summary(plot.df$New)
plot1.4 <- ggplot(data = plot.df, aes(x = day, y = New, colour = Source))+
  geom_point(aes(x=day, y=New),
             color = "indianred2",
             data = filter(plot.df, Source == "Composition")) +
  geom_point(aes(x = day, y=New),
             color = "turquoise",
             data = filter(plot.df, Source == "Resistome"))+
  geom_smooth(aes(x=day, y=New),
              color = "indianred4",
              fill = "rosybrown3",
              data = filter(plot.df, Source == "Composition"))+
  geom_smooth(aes(x = day, y = New),
              color = "turquoise4",
              fill = "paleturquoise3",              
              data = filter(plot.df, Source == "Resistome"))+
  theme_classic()+
  theme(text = element_text(size = 18),
        axis.text = element_text(color = "black"),
        plot.margin = unit(c(0.5,0.5,0.5,1), "cm"))+
  ylab("Number of New Species/ARGs")+
  xlab("Day Relative to HCT")+
  scale_y_continuous(limits = c(0, 100), breaks = c(0, 20, 40, 60, 80, 100), expand = c(0,0))+
  scale_x_continuous(limits = c(-30, 101), breaks = c(-30, 0, 30, 60, 100), expand = c(0,0))
plot1.4

summary(plot.df$inst.score)
plot2.2 <- ggplot(data = plot.df, aes(x = day, y = inst.score, colour = Source))+
  geom_point(aes(x=day, y=inst.score),
             color = "indianred2",
             data = filter(plot.df, Source == "Composition")) +
  geom_point(aes(x = day, y=inst.score),
             color = "turquoise",
             data = filter(plot.df, Source == "Resistome"))+
  geom_smooth(aes(x=day, y=inst.score),
              color = "indianred4",
              fill = "rosybrown3",
              data = filter(plot.df, Source == "Composition"))+
  geom_smooth(aes(x = day, y = inst.score),
              color = "turquoise4",
              fill = "paleturquoise3",              
              data = filter(plot.df, Source == "Resistome"))+
  theme_classic()+
  theme(text = element_text(size = 20),
        axis.text = element_text(color = "black"),
        plot.margin = unit(c(0.5,0.5,0.5,1), "cm"))+
  ylab("Jaccard Distance")+
  xlab("Day Relative to HCT")+
  scale_y_continuous(limits = c(0, 1.03), breaks = c(0, 0.25, 0.5, 0.75, 1), labels =c("0", "0.25", "0.5", "0.75", "1.0"), expand = c(0,0))+
  scale_x_continuous(limits = c(-30, 101), breaks = c(-30, 0, 30, 60, 100), expand = c(0,0))
plot2.2

#Repeated measures correlations
#Num ARG
rmc.out.num <- rmcorr(study_id,
                      num.ARG, num.species, res.comp, CI.level = 0.95, CIs =
                        c("analytic", "bootstrap"), nreps = 100,
                      bstrap.out = F)
rmc.out.num #r = 0.5241923 (95% CI 0.46-0.58), p<0.0001

plot(rmc.out.num,
     xlab = "Resistome",
     ylab = "Microbiome",
     overall = TRUE,
     overall.col = "gray60",
     overall.lwd = 3,
     overall.lty=2)

res.comp.num <- ggplot(res.comp, aes(x = num.ARG, y = num.species)) +
  geom_point()+
  theme_classic()+
  labs(x = "Number of ARGs", y = "Number of Species")+
  annotate(geom = "text", x=25, y= 115, label = expression('r'[rm]*'=0.52 (0.46, 0.58)'), size = 6)+
  scale_x_continuous(limits = c(0, 120), expand = c(0,0), breaks= c(0, 20, 40, 60, 80, 100, 120))+
  scale_y_continuous(limits = c(0,130), breaks = c(0, 20, 40, 60, 80, 100, 130), expand = c(0,0))+
  theme(text = element_text(size = 20, color = "black"), plot.title = element_text(size = 20, hjust = 0.5),
        axis.text = element_text(color = "black"),
        plot.margin = unit(c(1,1,0.5,0.5), "cm")) +
  geom_smooth(method = "lm", se = FALSE, size = 2)
res.comp.num

#number of new ARGs
rmc.out.new <- rmcorr(study_id,
                      new.ARG, new.sp, res.comp, CI.level = 0.95, CIs =
                        c("analytic", "bootstrap"), nreps = 100,
                      bstrap.out = F)
rmc.out.new #r = 0.527644 (95% CI 0.47-0.58), p<0.0001

res.comp.new <- ggplot(res.comp, aes(x = new.ARG, y = new.sp)) +
  geom_point()+
  theme_classic()+
  geom_smooth(method = "lm", se = FALSE, size = 2)+
  scale_y_continuous(limits =c(0, 100), breaks = c(0, 20, 40, 60 ,80 ,100), expand = c(0,0))+
  scale_x_continuous(limits =c(0, 100), breaks = c(0, 20, 40, 60 ,80 ,100), expand = c(0,0))+
  labs(x = "Number of New ARGs", y = "Number of New Species")+
  annotate(geom = "text", x=21, y= 90, label = expression('r'[rm]*'=0.53 (0.47, 0.58)'), size = 6)+
  theme(text = element_text(size = 20, color = "black"), plot.title = element_text(size = 20, hjust = 0.5),
        axis.text = element_text(color = "black"),
        plot.margin = unit(c(0.5,1,0.5,0.5), "cm"))
res.comp.new

#instability
rmc.out.inst.score <- rmcorr(study_id,
                             inst.score.ARG, inst.score.sp, res.comp, CI.level = 0.95, CIs =
                               c("analytic", "bootstrap"), nreps = 100,
                             bstrap.out = F)
rmc.out.inst.score #r = 0.6111856 (95% CI 0.56-0.66), p<0.0001

res.comp.inst.score <- ggplot(res.comp, aes(x = inst.score.ARG, y = inst.score.sp)) +
  geom_point()+
  theme_classic()+
  geom_smooth(method = "lm", se = FALSE, size = 2) +
  scale_y_continuous(limits = c(0, 1.05), breaks = c(0, 0.25, 0.5, 0.75, 1.0), labels = c("0", "0.25", "0.5", "0.75", "1.0"), expand = c(0,0))+
  scale_x_continuous(limits = c(0, 1.05), breaks = c(0, 0.25, 0.5, 0.75, 1.0), labels = c("0", "0.25", "0.5", "0.75", "1.0"), expand = c(0,0))+
  labs(x = "Jaccard Distance of ARGs", y = "Jaccard Distance of Species")+
  annotate(geom = "text", x=0.22, y= 1.02, label = expression('r'[rm]*'=0.61 (0.56, 0.66)'), size = 6)+
  theme(text = element_text(size = 20, color = "black"), plot.title = element_text(size = 20, hjust = 0.1),
        axis.text = element_text(color = "black"),
        plot.margin = unit(c(0.5, 1, 0.5, 0.5), "cm"))
res.comp.inst.score

#arranging Figure 3
ggarrange(plot1.2, res.comp.num, plot1.4, res.comp.new, plot2.2, res.comp.inst.score,
          ncol = 2, nrow = 3,
          labels = c("A","", "B", "", "C", ""), 
          #hjust = c(-0.1, 0, -0.2, 0, 0),
          vjust = 1,
          font.label = list(size = 24, color = "black"))+
  theme(plot.margin = margin())


#sensitivity testing of individual antibiotics (Supplementary Table 3)
#Number of ARGs
numARG.cefepimeexp.adj <- glmmadmb(num.ARG ~ cefepime.exp + age + sex + diagnosis + prep + type_hsct + day + log.depth + (1 |study_id), data  = test10, family = "nbinom1")
summary(numARG.cefepimeexp.adj)
numARG.cefepimeexp.adj.model <- cbind(Estimate = coef(numARG.cefepimeexp.adj), confint(numARG.cefepimeexp.adj))
cefepime.num <- data.frame(exp(numARG.cefepimeexp.adj.model))
exp(numARG.cefepimeexp.adj.model)
cefepime.num <- cefepime.num %>% rownames_to_column(var = "row") %>% filter(row == "cefepime.exp1") %>% mutate(effect = round(Estimate, 2), lower= round(X2.5.., 2), upper = round(X97.5.., 2)) %>%
  mutate(drug = if_else(row == "cefepime.exp1", "Cefepime", "one")) %>% dplyr::select(drug, effect, lower, upper)

numARG.vancexp.adj <- glmmadmb(num.ARG ~ vanc.exp + age + sex + diagnosis + prep + type_hsct + day + log.depth + (1 |study_id), data  = test10, family = "nbinom1")
summary(numARG.vancexp.adj)
numARG.vancexp.adj.model <- cbind(Estimate = coef(numARG.vancexp.adj), confint(numARG.vancexp.adj))
vanco.num <- data.frame(exp(numARG.vancexp.adj.model))
exp(numARG.vancexp.adj.model)
vanco.num <- vanco.num %>% rownames_to_column(var = "row") %>% filter(row == "vanc.exp1") %>% mutate(effect = round(Estimate, 2), lower= round(X2.5.., 2), upper = round(X97.5.., 2)) %>%
  mutate(drug = if_else(row == "vanc.exp1", "Vancomycin", "one")) %>% dplyr::select(drug, effect, lower, upper)

numARG.quinoloneexp.adj <- glmmadmb(num.ARG ~ quinolone.exp + age + sex + diagnosis + prep + type_hsct + day + log.depth + (1 |study_id), data  = test10, family = "nbinom1")
summary(numARG.quinoloneexp.adj)
numARG.quinoloneexp.adj.model <- cbind(Estimate = coef(numARG.quinoloneexp.adj), confint(numARG.quinoloneexp.adj))
quinolone.num <- data.frame(exp(numARG.quinoloneexp.adj.model))
exp(numARG.quinoloneexp.adj.model)
quinolone.num <- quinolone.num %>% rownames_to_column(var = "row") %>% filter(row == "quinolone.exp1") %>% mutate(effect = round(Estimate, 2), lower= round(X2.5.., 2), upper = round(X97.5.., 2)) %>%
  mutate(drug = if_else(row == "quinolone.exp1", "Fluoroquinolone", "one")) %>% dplyr::select(drug, effect, lower, upper)

numARG.agexp.adj <- glmmadmb(num.ARG ~ ag.exp + age + sex + diagnosis + prep + type_hsct + day + log.depth + (1 |study_id), data  = test10, family = "nbinom1")
summary(numARG.agexp.adj)
numARG.agexp.adj.model <- cbind(Estimate = coef(numARG.agexp.adj), confint(numARG.agexp.adj))
ag.num <- data.frame(exp(numARG.agexp.adj.model))
exp(numARG.agexp.adj.model)
ag.num <- ag.num %>% rownames_to_column(var = "row") %>% filter(row == "ag.exp1") %>% mutate(effect = round(Estimate, 2), lower= round(X2.5.., 2), upper = round(X97.5.., 2)) %>%
  mutate(drug = if_else(row == "ag.exp1", "Aminoglycoside", "one")) %>% dplyr::select(drug, effect, lower, upper)

numARG.macrolideexp.adj <- glmmadmb(num.ARG ~ macrolide.exp + age + sex + diagnosis + prep + type_hsct + day + log.depth + (1 |study_id), data  = test10, family = "nbinom1")
summary(numARG.macrolideexp.adj)
numARG.macrolideexp.adj.model <- cbind(Estimate = coef(numARG.macrolideexp.adj), confint(numARG.macrolideexp.adj))
macrolide.num <- data.frame(exp(numARG.macrolideexp.adj.model))
exp(numARG.macrolideexp.adj.model)
macrolide.num <- macrolide.num %>% rownames_to_column(var = "row") %>% filter(row == "macrolide.exp1") %>% mutate(effect = round(Estimate, 2), lower= round(X2.5.., 2), upper = round(X97.5.., 2)) %>%
  mutate(drug = if_else(row == "macrolide.exp1", "Macrolide", "one")) %>% dplyr::select(drug, effect, lower, upper)

numARG.piptazoexp.adj <- glmmadmb(num.ARG ~ piptazo.exp + age + sex + diagnosis + prep + type_hsct + day + log.depth + (1 |study_id), data  = test10, family = "nbinom1")
summary(numARG.piptazoexp.adj)
numARG.piptazoexp.adj.model <- cbind(Estimate = coef(numARG.piptazoexp.adj), confint(numARG.piptazoexp.adj))
piptazo.num <- data.frame(exp(numARG.piptazoexp.adj.model))
exp(numARG.piptazoexp.adj.model)
piptazo.num <- piptazo.num %>% rownames_to_column(var = "row") %>% filter(row == "piptazo.exp1") %>% mutate(effect = round(Estimate, 2), lower= round(X2.5.., 2), upper = round(X97.5.., 2)) %>%
  mutate(drug = if_else(row == "piptazo.exp1", "Pip-Tazo", "one")) %>% dplyr::select(drug, effect, lower, upper)

numARG.carbapenemexp.adj <- glmmadmb(num.ARG ~ carbapenem.exp + age + sex + diagnosis + prep + type_hsct + day + log.depth + (1 |study_id), data  = test10, family = "nbinom1")
summary(numARG.carbapenemexp.adj)
numARG.carbapenemexp.adj.model <- cbind(Estimate = coef(numARG.carbapenemexp.adj), confint(numARG.carbapenemexp.adj))
carbapenem.num <- data.frame(exp(numARG.carbapenemexp.adj.model))
exp(numARG.carbapenemexp.adj.model)
carbapenem.num <- carbapenem.num %>% rownames_to_column(var = "row") %>% filter(row == "carbapenem.exp1") %>% mutate(effect = round(Estimate, 2), lower= round(X2.5.., 2), upper = round(X97.5.., 2)) %>%
  mutate(drug = if_else(row == "carbapenem.exp1", "Carbapenem", "one")) %>% dplyr::select(drug, effect, lower, upper)

numARG.flagylexp.adj <- glmmadmb(num.ARG ~ flagyl.exp + age + sex + diagnosis + prep + type_hsct + day + log.depth + (1 |study_id), data  = test10, family = "nbinom1")
summary(numARG.flagylexp.adj)
numARG.flagylexp.adj.model <- cbind(Estimate = coef(numARG.flagylexp.adj), confint(numARG.flagylexp.adj))
flagyl.num <- data.frame(exp(numARG.flagylexp.adj.model))
exp(numARG.flagylexp.adj.model)
flagyl.num <- flagyl.num %>% rownames_to_column(var = "row") %>% filter(row == "flagyl.exp1") %>% mutate(effect = round(Estimate, 2), lower= round(X2.5.., 2), upper = round(X97.5.., 2)) %>%
  mutate(drug = if_else(row == "flagyl.exp1", "Metronidazole", "one")) %>% dplyr::select(drug, effect, lower, upper)

numARG.clindaexp.adj <- glmmadmb(num.ARG ~ clinda.exp + age + sex + diagnosis + prep + type_hsct + day + log.depth + (1 |study_id), data  = test10, family = "nbinom1")
summary(numARG.clindaexp.adj)
numARG.clindaexp.adj.model <- cbind(Estimate = coef(numARG.clindaexp.adj), confint(numARG.clindaexp.adj))
clinda.num <- data.frame(exp(numARG.clindaexp.adj.model))
exp(numARG.clindaexp.adj.model)
clinda.num <- clinda.num %>% rownames_to_column(var = "row") %>% filter(row == "clinda.exp1") %>% mutate(effect = round(Estimate, 2), lower= round(X2.5.., 2), upper = round(X97.5.., 2)) %>%
  mutate(drug = if_else(row == "clinda.exp1", "Clindamycin", "one")) %>% dplyr::select(drug, effect, lower, upper)

numARG.tmp_smxexp.adj <- glmmadmb(num.ARG ~ tmp_smx.exp + age + sex + diagnosis + prep + type_hsct + day + log.depth + (1 |study_id), data  = test10, family = "nbinom1")
summary(numARG.tmp_smxexp.adj)
numARG.tmp_smxexp.adj.model <- cbind(Estimate = coef(numARG.tmp_smxexp.adj), confint(numARG.tmp_smxexp.adj))
tmp_smx.num <- data.frame(exp(numARG.tmp_smxexp.adj.model))
exp(numARG.tmp_smxexp.adj.model)
tmp_smx.num <- tmp_smx.num %>% rownames_to_column(var = "row") %>% filter(row == "tmp_smx.exp1") %>% mutate(effect = round(Estimate, 2), lower= round(X2.5.., 2), upper = round(X97.5.., 2)) %>%
  mutate(drug = if_else(row == "tmp_smx.exp1", "TMP-SMX", "one")) %>% dplyr::select(drug, effect, lower, upper)


#New ARGs
newARG.cefepimeexp.adj <- glmmadmb(new.ARG ~ cefepime.exp + age + sex + diagnosis + prep + type_hsct + day + log.depth + (1 |study_id), data  = test10, family = "nbinom1")
summary(newARG.cefepimeexp.adj)
newARG.cefepimeexp.adj.model <- cbind(Estimate = coef(newARG.cefepimeexp.adj), confint(newARG.cefepimeexp.adj))
cefepime.new <- data.frame(exp(newARG.cefepimeexp.adj.model))
exp(newARG.cefepimeexp.adj.model)
cefepime.new <- cefepime.new %>% rownames_to_column(var = "row") %>% filter(row == "cefepime.exp1") %>% mutate(effect = round(Estimate, 2), lower= round(X2.5.., 2), upper = round(X97.5.., 2)) %>%
  mutate(drug = if_else(row == "cefepime.exp1", "Cefepime", "one")) %>% dplyr::select(drug, effect, lower, upper)

newARG.vancexp.adj <- glmmadmb(new.ARG ~ vanc.exp + age + sex + diagnosis + prep + type_hsct + day + log.depth + (1 |study_id), data  = test10, family = "nbinom1")
summary(newARG.vancexp.adj)
newARG.vancexp.adj.model <- cbind(Estimate = coef(newARG.vancexp.adj), confint(newARG.vancexp.adj))
vanc.new <- data.frame(exp(newARG.vancexp.adj.model))
exp(newARG.vancexp.adj.model)
vanc.new <- vanc.new %>% rownames_to_column(var = "row") %>% filter(row == "vanc.exp1") %>% mutate(effect = round(Estimate, 2), lower= round(X2.5.., 2), upper = round(X97.5.., 2)) %>%
  mutate(drug = if_else(row == "vanc.exp1", "Vancomycin", "one")) %>% dplyr::select(drug, effect, lower, upper)

newARG.quinoloneexp.adj <- glmmadmb(new.ARG ~ quinolone.exp + age + sex + diagnosis + prep + type_hsct + day + log.depth + (1 |study_id), data  = test10, family = "nbinom1")
summary(newARG.quinoloneexp.adj)
newARG.quinoloneexp.adj.model <- cbind(Estimate = coef(newARG.quinoloneexp.adj), confint(newARG.quinoloneexp.adj))
quinolone.new <- data.frame(exp(newARG.quinoloneexp.adj.model))
exp(newARG.quinoloneexp.adj.model)
quinolone.new <- quinolone.new %>% rownames_to_column(var = "row") %>% filter(row == "quinolone.exp1") %>% mutate(effect = round(Estimate, 2), lower= round(X2.5.., 2), upper = round(X97.5.., 2)) %>%
  mutate(drug = if_else(row == "quinolone.exp1", "Fluoroquinolone", "one")) %>% dplyr::select(drug, effect, lower, upper)

newARG.agexp.adj <- glmmadmb(new.ARG ~ ag.exp + age + sex + diagnosis + prep + type_hsct + day + log.depth + (1 |study_id), data  = test10, family = "nbinom1")
summary(newARG.agexp.adj)
newARG.agexp.adj.model <- cbind(Estimate = coef(newARG.agexp.adj), confint(newARG.agexp.adj))
ag.new <- data.frame(exp(newARG.agexp.adj.model))
exp(newARG.agexp.adj.model)
ag.new <- ag.new %>% rownames_to_column(var = "row") %>% filter(row == "ag.exp1") %>% mutate(effect = round(Estimate, 2), lower= round(X2.5.., 2), upper = round(X97.5.., 2)) %>%
  mutate(drug = if_else(row == "ag.exp1", "Aminoglycoside", "one")) %>% dplyr::select(drug, effect, lower, upper)

newARG.macrolideexp.adj <- glmmadmb(new.ARG ~ macrolide.exp + age + sex + diagnosis + prep + type_hsct + day + log.depth + (1 |study_id), data  = test10, family = "nbinom1")
summary(newARG.macrolideexp.adj)
newARG.macrolideexp.adj.model <- cbind(Estimate = coef(newARG.macrolideexp.adj), confint(newARG.macrolideexp.adj))
macrolide.new <- data.frame(exp(newARG.macrolideexp.adj.model))
exp(newARG.macrolideexp.adj.model)
macrolide.new <- macrolide.new %>% rownames_to_column(var = "row") %>% filter(row == "macrolide.exp1") %>% mutate(effect = round(Estimate, 2), lower= round(X2.5.., 2), upper = round(X97.5.., 2)) %>%
  mutate(drug = if_else(row == "macrolide.exp1", "Macrolide", "one")) %>% dplyr::select(drug, effect, lower, upper)

newARG.tmp_smxexp.adj <- glmmadmb(new.ARG ~ tmp_smx.exp + age + sex + diagnosis + prep + type_hsct + day + log.depth + (1 |study_id), data  = test10, family = "nbinom1")
summary(newARG.tmp_smxexp.adj)
newARG.tmp_smxexp.adj.model <- cbind(Estimate = coef(newARG.tmp_smxexp.adj), confint(newARG.tmp_smxexp.adj))
tmp_smx.new <- data.frame(exp(newARG.tmp_smxexp.adj.model))
exp(newARG.tmp_smxexp.adj.model)
tmp_smx.new <- tmp_smx.new %>% rownames_to_column(var = "row") %>% filter(row == "tmp_smx.exp1") %>% mutate(effect = round(Estimate, 2), lower= round(X2.5.., 2), upper = round(X97.5.., 2)) %>%
  mutate(drug = if_else(row == "tmp_smx.exp1", "TMP-SMX", "one")) %>% dplyr::select(drug, effect, lower, upper)

newARG.piptazoexp.adj <- glmmadmb(new.ARG ~ piptazo.exp + age + sex + diagnosis + prep + type_hsct + day + log.depth + (1 |study_id), data  = test10, family = "nbinom1")
summary(newARG.piptazoexp.adj)
newARG.piptazoexp.adj.model <- cbind(Estimate = coef(newARG.piptazoexp.adj), confint(newARG.piptazoexp.adj))
piptazo.new <- data.frame(exp(newARG.piptazoexp.adj.model))
exp(newARG.piptazoexp.adj.model)
piptazo.new <- piptazo.new %>% rownames_to_column(var = "row") %>% filter(row == "piptazo.exp1") %>% mutate(effect = round(Estimate, 2), lower= round(X2.5.., 2), upper = round(X97.5.., 2)) %>%
  mutate(drug = if_else(row == "piptazo.exp1", "Pip-Tazo", "one")) %>% dplyr::select(drug, effect, lower, upper)

newARG.carbapenemexp.adj <- glmmadmb(new.ARG ~ carbapenem.exp + age + sex + diagnosis + prep + type_hsct + day + log.depth + (1 |study_id), data  = test10, family = "nbinom1")
summary(newARG.carbapenemexp.adj)
newARG.carbapenemexp.adj.model <- cbind(Estimate = coef(newARG.carbapenemexp.adj), confint(newARG.carbapenemexp.adj))
carbapenem.new <- data.frame(exp(newARG.carbapenemexp.adj.model))
exp(newARG.carbapenemexp.adj.model)
carbapenem.new <- carbapenem.new %>% rownames_to_column(var = "row") %>% filter(row == "carbapenem.exp1") %>% mutate(effect = round(Estimate, 2), lower= round(X2.5.., 2), upper = round(X97.5.., 2)) %>%
  mutate(drug = if_else(row == "carbapenem.exp1", "Carbapenem", "one")) %>% dplyr::select(drug, effect, lower, upper)

newARG.flagylexp.adj <- glmmadmb(new.ARG ~ flagyl.exp + age + sex + diagnosis + prep + type_hsct + day + log.depth + (1 |study_id), data  = test10, family = "nbinom1")
summary(newARG.flagylexp.adj)
newARG.flagylexp.adj.model <- cbind(Estimate = coef(newARG.flagylexp.adj), confint(newARG.flagylexp.adj))
flagyl.new <- data.frame(exp(newARG.flagylexp.adj.model))
exp(newARG.flagylexp.adj.model)
flagyl.new <- flagyl.new %>% rownames_to_column(var = "row") %>% filter(row == "flagyl.exp1") %>% mutate(effect = round(Estimate, 2), lower= round(X2.5.., 2), upper = round(X97.5.., 2)) %>%
  mutate(drug = if_else(row == "flagyl.exp1", "Metronidazole", "one")) %>% dplyr::select(drug, effect, lower, upper)

newARG.clindaexp.adj <- glmmadmb(new.ARG ~ clinda.exp + age + sex + diagnosis + prep + type_hsct + day + log.depth + (1 |study_id), data  = test10, family = "nbinom1")
summary(newARG.clindaexp.adj)
newARG.clindaexp.adj.model <- cbind(Estimate = coef(newARG.clindaexp.adj), confint(newARG.clindaexp.adj))
clinda.new <- data.frame(exp(newARG.clindaexp.adj.model))
exp(newARG.clindaexp.adj.model)
clinda.new <- clinda.new %>% rownames_to_column(var = "row") %>% filter(row == "clinda.exp1") %>% mutate(effect = round(Estimate, 2), lower= round(X2.5.., 2), upper = round(X97.5.., 2)) %>%
  mutate(drug = if_else(row == "clinda.exp1", "Clindamycin", "one")) %>% dplyr::select(drug, effect, lower, upper)

#Instability
inst.cefepimeexp.adj <- lmer(inst.score ~ cefepime.exp + age + sex + diagnosis + prep + type_hsct + day + log.depth + (1 |study_id), data  = test10)
summary(inst.cefepimeexp.adj)
confint(inst.cefepimeexp.adj)
inst.cefepime.exp.adj.model <- data.frame(cbind(coef(summary(inst.cefepimeexp.adj))), confint(inst.cefepimeexp.adj, parm = c("(Intercept)", "cefepime.exp1", "age" , "sex2" , "diagnosisCongenital_immunodef", "diagnosisHeme_non-malignancy", "diagnosisMetabolic_disease", "diagnosisSolid_tumor", "prep2" , "type_hsctAllogeneic" , "day" , "log.depth" ))) %>%
  rownames_to_column(var = "row") %>% filter(row == "cefepime.exp1") %>%
  mutate(effect = round(Estimate, 3), lower = round(X2.5.., 3), upper = round(X97.5.., 3), drug = if_else(row == "cefepime.exp1", "Cefepime", "")) %>%
  dplyr::select(drug, effect, lower, upper)

inst.vancexp.adj <- lmer(inst.score ~ vanc.exp + age + sex + diagnosis + prep + type_hsct + day + log.depth + (1 |study_id), data  = test10)
summary(inst.vancexp.adj)
confint(inst.vancexp.adj)
inst.vanc.exp.adj.model <- data.frame(cbind(coef(summary(inst.vancexp.adj))), confint(inst.vancexp.adj, parm = c("(Intercept)", "vanc.exp1", "age" , "sex2" , "diagnosisCongenital_immunodef", "diagnosisHeme_non-malignancy", "diagnosisMetabolic_disease", "diagnosisSolid_tumor", "prep2" , "type_hsctAllogeneic" , "day" , "log.depth" ))) %>%
  rownames_to_column(var = "row") %>% filter(row == "vanc.exp1") %>%
  mutate(effect = round(Estimate, 3), lower = round(X2.5.., 3), upper = round(X97.5.., 3), drug = if_else(row == "vanc.exp1", "Vancomycin", "")) %>%
  dplyr::select(drug, effect, lower, upper)

inst.quinoloneexp.adj <- lmer(inst.score ~ quinolone.exp + age + sex + diagnosis + prep + type_hsct + day + log.depth + (1 |study_id), data  = test10)
summary(inst.quinoloneexp.adj)
confint(inst.quinoloneexp.adj)
inst.quinolone.exp.adj.model <- data.frame(cbind(coef(summary(inst.quinoloneexp.adj))), confint(inst.quinoloneexp.adj, parm = c("(Intercept)", "quinolone.exp1", "age" , "sex2" , "diagnosisCongenital_immunodef", "diagnosisHeme_non-malignancy", "diagnosisMetabolic_disease", "diagnosisSolid_tumor", "prep2" , "type_hsctAllogeneic" , "day" , "log.depth" ))) %>%
  rownames_to_column(var = "row") %>% filter(row == "quinolone.exp1") %>%
  mutate(effect = round(Estimate, 3), lower = round(X2.5.., 3), upper = round(X97.5.., 3), drug = if_else(row == "quinolone.exp1", "Fluoroquinolone", "")) %>%
  dplyr::select(drug, effect, lower, upper)

inst.agexp.adj <- lmer(inst.score ~ ag.exp + age + sex + diagnosis + prep + type_hsct + day + log.depth + (1 |study_id), data  = test10)
summary(inst.agexp.adj)
confint(inst.agexp.adj)
inst.ag.exp.adj.model <- data.frame(cbind(coef(summary(inst.agexp.adj))), confint(inst.agexp.adj, parm = c("(Intercept)", "ag.exp1", "age" , "sex2" , "diagnosisCongenital_immunodef", "diagnosisHeme_non-malignancy", "diagnosisMetabolic_disease", "diagnosisSolid_tumor", "prep2" , "type_hsctAllogeneic" , "day" , "log.depth" ))) %>%
  rownames_to_column(var = "row") %>% filter(row == "ag.exp1") %>%
  mutate(effect = round(Estimate, 3), lower = round(X2.5.., 3), upper = round(X97.5.., 3), drug = if_else(row == "ag.exp1", "Aminoglycoside", "")) %>%
  dplyr::select(drug, effect, lower, upper)

inst.macrolideexp.adj <- lmer(inst.score ~ macrolide.exp + age + sex + diagnosis + prep + type_hsct + day + log.depth + (1 |study_id), data  = test10)
summary(inst.macrolideexp.adj)
confint(inst.macrolideexp.adj)
inst.macrolide.exp.adj.model <- data.frame(cbind(coef(summary(inst.macrolideexp.adj))), confint(inst.macrolideexp.adj, parm = c("(Intercept)", "macrolide.exp1", "age" , "sex2" , "diagnosisCongenital_immunodef", "diagnosisHeme_non-malignancy", "diagnosisMetabolic_disease", "diagnosisSolid_tumor", "prep2" , "type_hsctAllogeneic" , "day" , "log.depth" ))) %>%
  rownames_to_column(var = "row") %>% filter(row == "macrolide.exp1") %>%
  mutate(effect = round(Estimate, 3), lower = round(X2.5.., 3), upper = round(X97.5.., 3), drug = if_else(row == "macrolide.exp1", "Macrolide", "")) %>%
  dplyr::select(drug, effect, lower, upper)

inst.tmp_smxexp.adj <- lmer(inst.score ~ tmp_smx.exp + age + sex + diagnosis + prep + type_hsct + day + log.depth + (1 |study_id), data  = test10)
summary(inst.tmp_smxexp.adj)
confint(inst.tmp_smxexp.adj)
inst.tmp_smx.exp.adj.model <- data.frame(cbind(coef(summary(inst.tmp_smxexp.adj))), confint(inst.tmp_smxexp.adj, parm = c("(Intercept)", "tmp_smx.exp1", "age" , "sex2" , "diagnosisCongenital_immunodef", "diagnosisHeme_non-malignancy", "diagnosisMetabolic_disease", "diagnosisSolid_tumor", "prep2" , "type_hsctAllogeneic" , "day" , "log.depth" ))) %>%
  rownames_to_column(var = "row") %>% filter(row == "tmp_smx.exp1") %>%
  mutate(effect = round(Estimate, 3), lower = round(X2.5.., 3), upper = round(X97.5.., 3), drug = if_else(row == "tmp_smx.exp1", "TMP-SMX", "")) %>%
  dplyr::select(drug, effect, lower, upper)

inst.piptazoexp.adj <- lmer(inst.score ~ piptazo.exp + age + sex + diagnosis + prep + type_hsct + day + log.depth + (1 |study_id), data  = test10)
summary(inst.piptazoexp.adj)
confint(inst.piptazoexp.adj)
inst.piptazo.exp.adj.model <- data.frame(cbind(coef(summary(inst.piptazoexp.adj))), confint(inst.piptazoexp.adj, parm = c("(Intercept)", "piptazo.exp1", "age" , "sex2" , "diagnosisCongenital_immunodef", "diagnosisHeme_non-malignancy", "diagnosisMetabolic_disease", "diagnosisSolid_tumor", "prep2" , "type_hsctAllogeneic" , "day" , "log.depth" ))) %>%
  rownames_to_column(var = "row") %>% filter(row == "piptazo.exp1") %>%
  mutate(effect = round(Estimate, 3), lower = round(X2.5.., 3), upper = round(X97.5.., 3), drug = if_else(row == "piptazo.exp1", "Pip-Tazo", "")) %>%
  dplyr::select(drug, effect, lower, upper)

inst.carbapenemexp.adj <- lmer(inst.score ~ carbapenem.exp + age + sex + diagnosis + prep + type_hsct + day + log.depth + (1 |study_id), data  = test10)
summary(inst.carbapenemexp.adj)
confint(inst.carbapenemexp.adj)
inst.carbapenem.exp.adj.model <- data.frame(cbind(coef(summary(inst.carbapenemexp.adj))), confint(inst.carbapenemexp.adj, parm = c("(Intercept)", "carbapenem.exp1", "age" , "sex2" , "diagnosisCongenital_immunodef", "diagnosisHeme_non-malignancy", "diagnosisMetabolic_disease", "diagnosisSolid_tumor", "prep2" , "type_hsctAllogeneic" , "day" , "log.depth" ))) %>%
  rownames_to_column(var = "row") %>% filter(row == "carbapenem.exp1") %>%
  mutate(effect = round(Estimate, 3), lower = round(X2.5.., 3), upper = round(X97.5.., 3), drug = if_else(row == "carbapenem.exp1", "Carbapenem", "")) %>%
  dplyr::select(drug, effect, lower, upper)

inst.flagylexp.adj <- lmer(inst.score ~ flagyl.exp + age + sex + diagnosis + prep + type_hsct + day + log.depth + (1 |study_id), data  = test10)
summary(inst.flagylexp.adj)
confint(inst.flagylexp.adj)
inst.flagyl.exp.adj.model <- data.frame(cbind(coef(summary(inst.flagylexp.adj))), confint(inst.flagylexp.adj, parm = c("(Intercept)", "flagyl.exp1", "age" , "sex2" , "diagnosisCongenital_immunodef", "diagnosisHeme_non-malignancy", "diagnosisMetabolic_disease", "diagnosisSolid_tumor", "prep2" , "type_hsctAllogeneic" , "day" , "log.depth" ))) %>%
  rownames_to_column(var = "row") %>% filter(row == "flagyl.exp1") %>%
  mutate(effect = round(Estimate, 3), lower = round(X2.5.., 3), upper = round(X97.5.., 3), drug = if_else(row == "flagyl.exp1", "Metronidazole", "")) %>%
  dplyr::select(drug, effect, lower, upper)

inst.clindaexp.adj <- lmer(inst.score ~ clinda.exp + age + sex + diagnosis + prep + type_hsct + day + log.depth + (1 |study_id), data  = test10)
summary(inst.clindaexp.adj)
confint(inst.clindaexp.adj)
inst.clinda.exp.adj.model <- data.frame(cbind(coef(summary(inst.clindaexp.adj))), confint(inst.clindaexp.adj, parm = c("(Intercept)", "clinda.exp1", "age" , "sex2" , "diagnosisCongenital_immunodef", "diagnosisHeme_non-malignancy", "diagnosisMetabolic_disease", "diagnosisSolid_tumor", "prep2" , "type_hsctAllogeneic" , "day" , "log.depth" ))) %>%
  rownames_to_column(var = "row") %>% filter(row == "clinda.exp1") %>%
  mutate(effect = round(Estimate, 3), lower = round(X2.5.., 3), upper = round(X97.5.., 3), drug = if_else(row == "clinda.exp1", "Clindamycin", "")) %>%
  dplyr::select(drug, effect, lower, upper)

#ARG Abundance
abundARG.cefepimeexp.adj <- glmmadmb(total.abund ~ cefepime.exp + age + sex + diagnosis + prep + type_hsct + day + log.depth + (1 |study_id), data  = test10, family = "nbinom1")
summary(abundARG.cefepimeexp.adj)
abundARG.cefepimeexp.adj.model <- cbind(Estimate = coef(abundARG.cefepimeexp.adj), confint(abundARG.cefepimeexp.adj))
cefepime.rpkm <- data.frame(exp(abundARG.cefepimeexp.adj.model))
exp(abundARG.cefepimeexp.adj.model)
cefepime.rpkm <- cefepime.rpkm %>% rownames_to_column(var = "row") %>% filter(row == "cefepime.exp1") %>% mutate(effect = round(Estimate, 2), lower= round(X2.5.., 2), upper = round(X97.5.., 2)) %>%
  mutate(drug = if_else(row == "cefepime.exp1", "Cefepime", "one")) %>% dplyr::select(drug, effect, lower, upper)

abundARG.vancexp.adj <- glmmadmb(total.abund ~ vanc.exp + age + sex + diagnosis + prep + type_hsct + day + log.depth + (1 |study_id), data  = test10, family = "nbinom1")
summary(abundARG.vancexp.adj)
abundARG.vancexp.adj.model <- cbind(Estimate = coef(abundARG.vancexp.adj), confint(abundARG.vancexp.adj))
vanc.rpkm <- data.frame(exp(abundARG.vancexp.adj.model))
exp(abundARG.vancexp.adj.model)
vanc.rpkm <- vanc.rpkm %>% rownames_to_column(var = "row") %>% filter(row == "vanc.exp1") %>% mutate(effect = round(Estimate, 2), lower= round(X2.5.., 2), upper = round(X97.5.., 2)) %>%
  mutate(drug = if_else(row == "vanc.exp1", "Vancomycin", "one")) %>% dplyr::select(drug, effect, lower, upper)

abundARG.quinoloneexp.adj <- glmmadmb(total.abund ~ quinolone.exp + age + sex + diagnosis + prep + type_hsct + day + log.depth + (1 |study_id), data  = test10, family = "nbinom2")
summary(abundARG.quinoloneexp.adj)
abundARG.quinoloneexp.adj.model <- cbind(Estimate = coef(abundARG.quinoloneexp.adj), confint(abundARG.quinoloneexp.adj))
quinolone.rpkm <- data.frame(exp(abundARG.quinoloneexp.adj.model))
exp(abundARG.quinoloneexp.adj.model)
quinolone.rpkm <- quinolone.rpkm %>% rownames_to_column(var = "row") %>% filter(row == "quinolone.exp1") %>% mutate(effect = round(Estimate, 2), lower= round(X2.5.., 2), upper = round(X97.5.., 2)) %>%
  mutate(drug = if_else(row == "quinolone.exp1", "Fluoroquinolone", "one")) %>% dplyr::select(drug, effect, lower, upper)

abundARG.agexp.adj <- glmmadmb(total.abund ~ ag.exp + age + sex + diagnosis + prep + type_hsct + day + log.depth + (1 |study_id), data  = test10, family = "nbinom2")
summary(abundARG.agexp.adj)
abundARG.agexp.adj.model <- cbind(Estimate = coef(abundARG.agexp.adj), confint(abundARG.agexp.adj))
ag.rpkm <- data.frame(exp(abundARG.agexp.adj.model))
exp(abundARG.agexp.adj.model)
ag.rpkm <- ag.rpkm %>% rownames_to_column(var = "row") %>% filter(row == "ag.exp1") %>% mutate(effect = round(Estimate, 2), lower= round(X2.5.., 2), upper = round(X97.5.., 2)) %>%
  mutate(drug = if_else(row == "ag.exp1", "Aminoglycoside", "one")) %>% dplyr::select(drug, effect, lower, upper)

abundARG.macrolideexp.adj <- glmmadmb(total.abund ~ macrolide.exp + age + sex + diagnosis + prep + type_hsct + day + log.depth + (1 |study_id), data  = test10, family = "nbinom1")
summary(abundARG.macrolideexp.adj)
abundARG.macrolideexp.adj.model <- cbind(Estimate = coef(abundARG.macrolideexp.adj), confint(abundARG.macrolideexp.adj))
macrolide.rpkm <- data.frame(exp(abundARG.macrolideexp.adj.model))
exp(abundARG.macrolideexp.adj.model)
macrolide.rpkm <- macrolide.rpkm %>% rownames_to_column(var = "row") %>% filter(row == "macrolide.exp1") %>% mutate(effect = round(Estimate, 2), lower= round(X2.5.., 2), upper = round(X97.5.., 2)) %>%
  mutate(drug = if_else(row == "macrolide.exp1", "Macrolide", "one")) %>% dplyr::select(drug, effect, lower, upper)

abundARG.tmp_smxexp.adj <- glmmadmb(total.abund ~ tmp_smx.exp + age + sex + diagnosis + prep + type_hsct + day + log.depth + (1 |study_id), data  = test10, family = "nbinom1")
summary(abundARG.tmp_smxexp.adj)
abundARG.tmp_smxexp.adj.model <- cbind(Estimate = coef(abundARG.tmp_smxexp.adj), confint(abundARG.tmp_smxexp.adj))
tmp_smx.rpkm <- data.frame(exp(abundARG.tmp_smxexp.adj.model))
exp(abundARG.tmp_smxexp.adj.model)
tmp_smx.rpkm <- tmp_smx.rpkm %>% rownames_to_column(var = "row") %>% filter(row == "tmp_smx.exp1") %>% mutate(effect = round(Estimate, 2), lower= round(X2.5.., 2), upper = round(X97.5.., 2)) %>%
  mutate(drug = if_else(row == "tmp_smx.exp1", "TMP-SMX", "one")) %>% dplyr::select(drug, effect, lower, upper)

abundARG.piptazoexp.adj <- glmmadmb(total.abund ~ piptazo.exp + age + sex + diagnosis + prep + type_hsct + day + log.depth + (1 |study_id), data  = test10, family = "nbinom1")
summary(abundARG.piptazoexp.adj)
abundARG.piptazoexp.adj.model <- cbind(Estimate = coef(abundARG.piptazoexp.adj), confint(abundARG.piptazoexp.adj))
piptazo.rpkm <- data.frame(exp(abundARG.piptazoexp.adj.model))
exp(abundARG.piptazoexp.adj.model)
piptazo.rpkm <- piptazo.rpkm %>% rownames_to_column(var = "row") %>% filter(row == "piptazo.exp1") %>% mutate(effect = round(Estimate, 2), lower= round(X2.5.., 2), upper = round(X97.5.., 2)) %>%
  mutate(drug = if_else(row == "piptazo.exp1", "Pip-Tazo", "one")) %>% dplyr::select(drug, effect, lower, upper)

abundARG.carbapenemexp.adj <- glmmadmb(total.abund ~ carbapenem.exp + age + sex + diagnosis + prep + type_hsct + day + log.depth + (1 |study_id), data  = test10, family = "nbinom1")
summary(abundARG.carbapenemexp.adj)
abundARG.carbapenemexp.adj.model <- cbind(Estimate = coef(abundARG.carbapenemexp.adj), confint(abundARG.carbapenemexp.adj))
carbapenem.rpkm <- data.frame(exp(abundARG.carbapenemexp.adj.model))
exp(abundARG.carbapenemexp.adj.model)
carbapenem.rpkm <- carbapenem.rpkm %>% rownames_to_column(var = "row") %>% filter(row == "carbapenem.exp1") %>% mutate(effect = round(Estimate, 2), lower= round(X2.5.., 2), upper = round(X97.5.., 2)) %>%
  mutate(drug = if_else(row == "carbapenem.exp1", "Carbapenem", "one")) %>% dplyr::select(drug, effect, lower, upper)

abundARG.flagylexp.adj <- glmmadmb(total.abund ~ flagyl.exp + age + sex + diagnosis + prep + type_hsct + day + log.depth + (1 |study_id), data  = test10, family = "nbinom1")
summary(abundARG.flagylexp.adj)
abundARG.flagylexp.adj.model <- cbind(Estimate = coef(abundARG.flagylexp.adj), confint(abundARG.flagylexp.adj))
flagyl.rpkm <- data.frame(exp(abundARG.flagylexp.adj.model))
exp(abundARG.flagylexp.adj.model)
flagyl.rpkm <- flagyl.rpkm %>% rownames_to_column(var = "row") %>% filter(row == "flagyl.exp1") %>% mutate(effect = round(Estimate, 2), lower= round(X2.5.., 2), upper = round(X97.5.., 2)) %>%
  mutate(drug = if_else(row == "flagyl.exp1", "Metronidazole", "one")) %>% dplyr::select(drug, effect, lower, upper)

abundARG.clindaexp.adj <- glmmadmb(total.abund ~ clinda.exp + age + sex + diagnosis + prep + type_hsct + day + log.depth + (1 |study_id), data  = test10, family = "nbinom2")
summary(abundARG.clindaexp.adj)
abundARG.clindaexp.adj.model <- cbind(Estimate = coef(abundARG.clindaexp.adj), confint(abundARG.clindaexp.adj))
clinda.rpkm <- data.frame(exp(abundARG.clindaexp.adj.model))
exp(abundARG.clindaexp.adj.model)
clinda.rpkm <- clinda.rpkm %>% rownames_to_column(var = "row") %>% filter(row == "clinda.exp1") %>% mutate(effect = round(Estimate, 2), lower= round(X2.5.., 2), upper = round(X97.5.., 2)) %>%
  mutate(drug = if_else(row == "clinda.exp1", "Clindamycin", "one")) %>% dplyr::select(drug, effect, lower, upper)

#Figure 4
fp.num <- rbind(anaerobe.num, piptazo.num, carbapenem.num, flagyl.num, clinda.num, cefepime.num, vanco.num, quinolone.num, ag.num, macrolide.num, tmp_smx.num) %>%
  mutate(index = if_else(drug == "Aerobic Antibiotics", 12, if_else(drug == 'Cefepime', 11, if_else(drug== "Vancomycin", 10, if_else(drug == "Fluoroquinolone", 9, if_else(drug== "Aminoglycoside", 8, if_else(drug=="Macrolide", 7, if_else(drug=="TMP-SMX", 6, if_else(drug=="Anaerobic Antibiotics", 5, if_else(drug=="Pip-Tazo", 4, if_else(drug=="Carbapenem", 3, if_else(drug== "Metronidazole", 2, 1)))))))))))) %>%
  arrange(index) %>% mutate(color = if_else(drug == "Aerobic Antibiotics", "Main", if_else(drug == 'Cefepime', "Sub", if_else(drug== "Vancomycin", "Sub", if_else(drug == "Fluoroquinolone", "Sub", if_else(drug== "Aminoglycoside", "Sub", if_else(drug=="Macrolide", "Sub", if_else(drug=="TMP-SMX", "Sub", if_else(drug=="Anaerobic Antibiotics", "Main", if_else(drug=="Pip-Tazo", "Sub1", if_else(drug=="Carbapenem", "Sub1", if_else(drug== "Metronidazole", "Sub1", "Sub1"))))))))))))
a <- ifelse(fp.num$color == "Main", "black", ifelse(fp.num$color== "Sub", "dodgerblue3", "forestgreen"))
b <- ifelse(fp.num$color == "Main", "bold", "plain")

num.h <- ggplot(data = fp.num, aes(y=index, x=effect, xmin=lower, xmax = upper)) +
  geom_point(color = a, size = 4) +
  geom_errorbarh(height = 0.2, color = a, size = 1.25) +
  scale_y_continuous(breaks = 1:12, labels = fp.num$drug, expand=c(0, 0.5)) +
  scale_x_continuous(limits = c(0.5, 1.5), breaks = c(0.5, 0.75, 1.0, 1.25, 1.5), labels = c("0.5", "0.75", "1.0", "1.25", "1.5")) +
  labs(title = "Number of ARGs", x = expression(beta)) +
  geom_vline(xintercept = 1, color = "black", linetype = "longdash", alpha = 0.5) +
  theme(panel.background = element_blank(),
        plot.title = element_text(size = 20, hjust = 0.5),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_blank(),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(size = 15, color = a, face = b),
        axis.line.x = element_line(linetype = "solid"),
        panel.grid = element_blank(),
        axis.ticks.y = element_blank())
num.h

fp.new <- rbind(anaerobe.new, piptazo.new, carbapenem.new, flagyl.new, clinda.new, cefepime.new, vanc.new, quinolone.new, ag.new, macrolide.new, tmp_smx.new) %>%
  mutate(index = if_else(drug == "Aerobic Antibiotics", 12, if_else(drug == 'Cefepime', 11, if_else(drug== "Vancomycin", 10, if_else(drug == "Fluoroquinolone", 9, if_else(drug== "Aminoglycoside", 8, if_else(drug=="Macrolide", 7, if_else(drug=="TMP-SMX", 6, if_else(drug=="Anaerobic Antibiotics", 5, if_else(drug=="Pip-Tazo", 4, if_else(drug=="Carbapenem", 3, if_else(drug== "Metronidazole", 2, 1)))))))))))) %>%
  arrange(index) %>% mutate(color = if_else(drug == "Aerobic Antibiotics", "Main", if_else(drug == 'Cefepime', "Sub", if_else(drug== "Vancomycin", "Sub", if_else(drug == "Fluoroquinolone", "Sub", if_else(drug== "Aminoglycoside", "Sub", if_else(drug=="Macrolide", "Sub", if_else(drug=="TMP-SMX", "Sub", if_else(drug=="Anaerobic Antibiotics", "Main", if_else(drug=="Pip-Tazo", "Sub1", if_else(drug=="Carbapenem", "Sub1", if_else(drug== "Metronidazole", "Sub1", "Sub1"))))))))))))

new.h <- ggplot(data = fp.new, aes(y=index, x=effect, xmin=lower, xmax = upper)) +
  geom_point(color = a, size = 4) +
  geom_errorbarh(height = 0.2, color = a, size = 1.25) +
  scale_y_continuous(breaks = 1:12, labels = fp.new$drug, expand=c(0, 0.5)) +
  scale_x_continuous(limits = c(0.5, 2.0), breaks = c(0.5, 0.7, 1.0, 1.3, 1.6, 2.0)) +
  labs(title = "Number of New ARGs", x = expression(beta)) +
  geom_vline(xintercept = 1, color = "black", linetype = "longdash", alpha = 0.5) +
  theme(panel.background = element_blank(),
        plot.title = element_text(size = 20, hjust = 0.5),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_blank(),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_blank(),
        axis.line.x = element_line(linetype = "solid"),
        panel.grid = element_blank(),
        axis.ticks.y = element_blank())
new.h

fp.inst <- rbind(inst.anaerobe.exp.adj.model, inst.cefepime.exp.adj.model, inst.vanc.exp.adj.model, inst.quinolone.exp.adj.model, inst.ag.exp.adj.model, inst.macrolide.exp.adj.model, inst.tmp_smx.exp.adj.model, inst.piptazo.exp.adj.model, inst.carbapenem.exp.adj.model, inst.flagyl.exp.adj.model, inst.clinda.exp.adj.model) %>%
  mutate(index = if_else(drug == "Aerobic Antibiotics", 12, if_else(drug == 'Cefepime', 11, if_else(drug== "Vancomycin", 10, if_else(drug == "Fluoroquinolone", 9, if_else(drug== "Aminoglycoside", 8, if_else(drug=="Macrolide", 7, if_else(drug=="TMP-SMX", 6, if_else(drug=="Anaerobic Antibiotics", 5, if_else(drug=="Pip-Tazo", 4, if_else(drug=="Carbapenem", 3, if_else(drug== "Metronidazole", 2, 1)))))))))))) %>%
  arrange(index) %>% mutate(color = if_else(drug == "Aerobic Antibiotics", "Main", if_else(drug == 'Cefepime', "Sub", if_else(drug== "Vancomycin", "Sub", if_else(drug == "Fluoroquinolone", "Sub", if_else(drug== "Aminoglycoside", "Sub", if_else(drug=="Macrolide", "Sub", if_else(drug=="TMP-SMX", "Sub", if_else(drug=="Anaerobic Antibiotics", "Main", if_else(drug=="Pip-Tazo", "Sub1", if_else(drug=="Carbapenem", "Sub1", if_else(drug== "Metronidazole", "Sub1", "Sub1"))))))))))))

inst.h <- ggplot(data = fp.inst, aes(y=index, x=effect, xmin=lower, xmax = upper)) +
  geom_point(color = a, size = 4) +
  geom_errorbarh(height = 0.2, color = a, size = 1.25) +
  scale_y_continuous(breaks = 1:12, labels = fp.inst$drug, expand=c(0, 0.5)) +
  scale_x_continuous(limits = c(-0.25, 0.25), breaks = c(-0.25, -0.1, 0, 0.1, 0.25), labels = c("-0.25", "-0.1", "0.0", "0.1", "0.25")) +
  labs(title = "Instability of ARGs", x = "Change in Jaccard Distance") +
  geom_vline(xintercept = 0, color = "black", linetype = "longdash", alpha = 0.5) +
  theme(panel.background = element_blank(),
        plot.title = element_text(size = 20, hjust = 0.5),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_blank(),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_blank(),
        axis.line.x = element_line(linetype = "solid"),
        panel.grid = element_blank(),
        axis.ticks.y = element_blank())
inst.h

fp.rpkm <- rbind(anaerobe.rpkm, piptazo.rpkm, carbapenem.rpkm, flagyl.rpkm, clinda.rpkm, cefepime.rpkm, vanc.rpkm, quinolone.rpkm, ag.rpkm, macrolide.rpkm, tmp_smx.rpkm) %>%
  mutate(index = if_else(drug == "Aerobic Antibiotics", 12, if_else(drug == 'Cefepime', 11, if_else(drug== "Vancomycin", 10, if_else(drug == "Fluoroquinolone", 9, if_else(drug== "Aminoglycoside", 8, if_else(drug=="Macrolide", 7, if_else(drug=="TMP-SMX", 6, if_else(drug=="Anaerobic Antibiotics", 5, if_else(drug=="Pip-Tazo", 4, if_else(drug=="Carbapenem", 3, if_else(drug== "Metronidazole", 2, 1)))))))))))) %>%
  arrange(index) %>% mutate(color = if_else(drug == "Aerobic Antibiotics", "Main", if_else(drug == 'Cefepime', "Sub", if_else(drug== "Vancomycin", "Sub", if_else(drug == "Fluoroquinolone", "Sub", if_else(drug== "Aminoglycoside", "Sub", if_else(drug=="Macrolide", "Sub", if_else(drug=="TMP-SMX", "Sub", if_else(drug=="Anaerobic Antibiotics", "Main", if_else(drug=="Pip-Tazo", "Sub1", if_else(drug=="Carbapenem", "Sub1", if_else(drug== "Metronidazole", "Sub1", "Sub1"))))))))))))

rpkm.h <- ggplot(data = fp.rpkm, aes(y=index, x=effect, xmin=lower, xmax = upper)) +
  geom_point(color = a, size = 4) +
  geom_errorbarh(height = 0.2, color = a, size = 1.25) +
  scale_y_continuous(breaks = 1:12, labels = fp.rpkm$drug, expand=c(0, 0.5)) +
  scale_x_continuous(limits = c(0.5, 3.6), breaks = c(0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.6), labels = c("0.5", "1.0", "1.5", "2.0", "2.5", "3.0", "3.6")) +
  labs(title = "Abundance of ARGs", x = expression(beta)) +
  geom_vline(xintercept = 1, color = "black", linetype = "longdash", alpha = 0.5) +
  theme(panel.background = element_blank(),
        plot.title = element_text(size = 20, hjust = 0.5),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_blank(),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_blank(),
        axis.line.x = element_line(linetype = "solid"),
        panel.grid = element_blank(),
        axis.ticks.y = element_blank())
rpkm.h

ggarrange(num.h, new.h, rpkm.h, inst.h, nrow =1, widths = c(1.7,1,1,1), align = "h")

#Individual antibiotics on individual ARG class abundance and numbers (Supplementary Tables 4 and 5)
#cefepime exposure
rpkm.cefepime.betalactam <- glmmadmb(beta.abund ~ cefepime.exp + (1|study_id), data = test9, family = "nbinom1")
summary(rpkm.cefepime.betalactam)
rpkm.cefepime.betalactam.model <- cbind(Estimate = coef(rpkm.cefepime.betalactam), confint(rpkm.cefepime.betalactam))
exp(rpkm.cefepime.betalactam.model)

num.cefepime.betalactam <- glmmadmb(num.betalactam.ARG ~ cefepime.exp + (1|study_id), data = test9, family = "nbinom1")
summary(num.cefepime.betalactam)
num.cefepime.betalactam.model <- cbind(Estimate = coef(num.cefepime.betalactam), confint(num.cefepime.betalactam))
exp(num.cefepime.betalactam.model)

rpkm.cefepime.vanc <- glmmadmb(vanc.abund ~ cefepime.exp + (1|study_id), data = test9, family = "nbinom1")
summary(rpkm.cefepime.vanc)
rpkm.cefepime.vanc.model <- cbind(Estimate = coef(rpkm.cefepime.vanc), confint(rpkm.cefepime.vanc))
exp(rpkm.cefepime.vanc.model)

num.cefepime.vanc <- glmmadmb(num.vanc.ARG ~ cefepime.exp + (1|study_id), data = test9, family = "nbinom1")
summary(num.cefepime.vanc)
num.cefepime.vanc.model <- cbind(Estimate = coef(num.cefepime.vanc), confint(num.cefepime.vanc))
exp(num.cefepime.vanc.model)

rpkm.cefepime.carbapenem <- glmmadmb(carb.abund ~ cefepime.exp + (1|study_id), data = test9, family = "nbinom1")
summary(rpkm.cefepime.carbapenem)
rpkm.cefepime.carbapenem.model <- cbind(Estimate = coef(rpkm.cefepime.carbapenem), confint(rpkm.cefepime.carbapenem))
exp(rpkm.cefepime.carbapenem.model)

num.cefepime.carbapenem <- glmmadmb(num.carbapenem.ARG ~ cefepime.exp + (1|study_id), data = test9, family = "nbinom1")
summary(num.cefepime.carbapenem)
num.cefepime.carbapenem.model <- cbind(Estimate = coef(num.cefepime.carbapenem), confint(num.cefepime.carbapenem))
exp(num.cefepime.carbapenem.model)

rpkm.cefepime.quinolone <- glmmadmb(fluoro.abund ~ cefepime.exp + (1|study_id), data = test9, family = "nbinom1")
summary(rpkm.cefepime.quinolone)
rpkm.cefepime.quinolone.model <- cbind(Estimate = coef(rpkm.cefepime.quinolone), confint(rpkm.cefepime.quinolone))
exp(rpkm.cefepime.quinolone.model)

num.cefepime.quinolone <- glmmadmb(num.quinolone.ARG ~ cefepime.exp + (1|study_id), data = test9, family = "nbinom1")
summary(num.cefepime.quinolone)
num.cefepime.quinolone.model <- cbind(Estimate = coef(num.cefepime.quinolone), confint(num.cefepime.quinolone))
exp(num.cefepime.quinolone.model)

rpkm.cefepime.amino <- glmmadmb(amino.abund ~ cefepime.exp + (1|study_id), data = test9, family = "nbinom1")
summary(rpkm.cefepime.amino)
rpkm.cefepime.amino.model <- cbind(Estimate = coef(rpkm.cefepime.amino), confint(rpkm.cefepime.amino))
exp(rpkm.cefepime.amino.model)

num.cefepime.amino <- glmmadmb(num.amino.ARG ~ cefepime.exp + (1|study_id), data = test9, family = "nbinom1")
summary(num.cefepime.amino)
num.cefepime.amino.model <- cbind(Estimate = coef(num.cefepime.amino), confint(num.cefepime.amino))
exp(num.cefepime.amino.model)

rpkm.cefepime.clinda <- glmmadmb(clinda.abund ~ cefepime.exp + (1|study_id), data = test9, family = "nbinom1")
summary(rpkm.cefepime.clinda)
rpkm.cefepime.clinda.model <- cbind(Estimate = coef(rpkm.cefepime.clinda), confint(rpkm.cefepime.clinda))
exp(rpkm.cefepime.clinda.model)

num.cefepime.clinda <- glmmadmb(num.clinda.ARG ~ cefepime.exp + (1|study_id), data = test9, family = "nbinom1")
summary(num.cefepime.clinda)
num.cefepime.clinda.model <- cbind(Estimate = coef(num.cefepime.clinda), confint(num.cefepime.clinda))
exp(num.cefepime.clinda.model)

rpkm.cefepime.macrolide <- glmmadmb(macrolide.abund ~ cefepime.exp + (1|study_id), data = test9, family = "nbinom2")
summary(rpkm.cefepime.macrolide)
rpkm.cefepime.macrolide.model <- cbind(Estimate = coef(rpkm.cefepime.macrolide), confint(rpkm.cefepime.macrolide))
exp(rpkm.cefepime.macrolide.model)

num.cefepime.macrolide <- glmmadmb(num.macrolide.ARG ~ cefepime.exp + (1|study_id), data = test9, family = "nbinom1")
summary(num.cefepime.macrolide)
num.cefepime.macrolide.model <- cbind(Estimate = coef(num.cefepime.macrolide), confint(num.cefepime.macrolide))
exp(num.cefepime.macrolide.model)

rpkm.cefepime.trimeth <- glmmadmb(trimeth.abund ~ cefepime.exp + (1|study_id), data = test9, family = "nbinom1")
summary(rpkm.cefepime.trimeth)
rpkm.cefepime.trimeth.model <- cbind(Estimate = coef(rpkm.cefepime.trimeth), confint(rpkm.cefepime.trimeth))
exp(rpkm.cefepime.trimeth.model)

num.cefepime.trimeth <- glmmadmb(num.trimeth.ARG ~ cefepime.exp + (1|study_id), data = test9, family = "nbinom2")
summary(num.cefepime.trimeth)
num.cefepime.trimeth.model <- cbind(Estimate = coef(num.cefepime.trimeth), confint(num.cefepime.trimeth))
exp(num.cefepime.trimeth.model)

rpkm.cefepime.tetra <- glmmadmb(tetra.abund ~ cefepime.exp + (1|study_id), data = test9, family = "nbinom1")
summary(rpkm.cefepime.tetra)
rpkm.cefepime.tetra.model <- cbind(Estimate = coef(rpkm.cefepime.tetra), confint(rpkm.cefepime.tetra))
exp(rpkm.cefepime.tetra.model)

num.cefepime.tetra <- glmmadmb(num.tetra.ARG ~ cefepime.exp + (1|study_id), data = test9, family = "nbinom1")
summary(num.cefepime.tetra)
num.cefepime.tetra.model <- cbind(Estimate = coef(num.cefepime.tetra), confint(num.cefepime.tetra))
exp(num.cefepime.tetra.model)

rpkm.cefepime.flagyl <- glmmadmb(flagyl.abund ~ cefepime.exp + (1|study_id), data = test9, family = "nbinom1")
summary(rpkm.cefepime.flagyl)
rpkm.cefepime.flagyl.model <- cbind(Estimate = coef(rpkm.cefepime.flagyl), confint(rpkm.cefepime.flagyl))
exp(rpkm.cefepime.flagyl.model)
#did not converge
#num.cefepime.flagyl <- glmmadmb(num.flagyl.ARG ~ cefepime.exp + (1|study_id), data = test9, family = "nbinom2")
#summary(num.cefepime.flagyl)
#num.cefepime.flagyl.model <- cbind(Estimate = coef(num.cefepime.flagyl), confint(num.cefepime.flagyl))
#exp(num.cefepime.flagyl.model)

#adjusted P value for multiple comparisons
pvalues <- c(  0.25,   0.32,   1.3e-06,   0.039,   0.64,   0.75,   0.0026,   0.023,   0.26,   7.1e-09,   1.1e-09,   0.84,   3.1e-07,   4.8e-09,   4.5e-08,   0.00014,   2.4e-05,   0.0026,   1.1e-11,   NA)
p.adjust(pvalues, method = "BH")
arg.class <- c("BetalactamAbund", "VancomycinAbund", "CarbapenemAbund", "QuinoloneAbund", "AminoglycosidAbund", "ClindamycinAbund", "MacrolideAbund", "TrimethoprimAbund", "TetracyclineAbund", "MetronidazoleAbund",
               "Betalactam.num", "Vancomycin.num", "Carbapenem.num", "Quinolone.num", "Aminoglycosid.num", "Clindamycin.num", "Macrolide.num", "Trimethoprim.num", "Tetracycline.num", "Metronidazole.num")
adj.p.cefepime<- c(
  3.293333e-01,3.800000e-01 ,3.528571e-06 ,5.700000e-02, 7.152941e-01, 7.916667e-01, 4.490909e-03, 3.641667e-02,
  3.293333e-01,3.372500e-08 ,1.045000e-08 ,8.400000e-01, 9.816667e-07, 3.040000e-08, 1.710000e-07, 2.955556e-04,
  5.700000e-05, 4.490909e-03, 2.090000e-10,           NA)
cefepime.bh <- data.frame(adj.p.cefepime, rownames = arg.class) %>% column_to_rownames(var = "rownames")
cefepime.bh <- cefepime.bh %>% format(adj.p.cefepime, scientific = FALSE) 
cefepime.bh$adj.p.cefepime <- as.numeric(cefepime.bh$adj.p.cefepime)
cefepime.bh$adj.p.cefepime <- signif(cefepime.bh$adj.p.cefepime, 2)


#piptazo exposure
rpkm.piptazo.betalactam <- glmmadmb(beta.abund ~ piptazo.exp + (1|study_id), data = test9, family = "nbinom1")
summary(rpkm.piptazo.betalactam)
rpkm.piptazo.betalactam.model <- cbind(Estimate = coef(rpkm.piptazo.betalactam), confint(rpkm.piptazo.betalactam))
exp(rpkm.piptazo.betalactam.model)

num.piptazo.betalactam <- glmmadmb(num.betalactam.ARG ~ piptazo.exp + (1|study_id), data = test9, family = "nbinom1")
summary(num.piptazo.betalactam)
num.piptazo.betalactam.model <- cbind(Estimate = coef(num.piptazo.betalactam), confint(num.piptazo.betalactam))
exp(num.piptazo.betalactam.model)

rpkm.piptazo.vanc <- glmmadmb(vanc.abund ~ piptazo.exp + (1|study_id), data = test9, family = "nbinom1")
summary(rpkm.piptazo.vanc)
rpkm.piptazo.vanc.model <- cbind(Estimate = coef(rpkm.piptazo.vanc), confint(rpkm.piptazo.vanc))
exp(rpkm.piptazo.vanc.model)

num.piptazo.vanc <- glmmadmb(num.vanc.ARG ~ piptazo.exp + (1|study_id), data = test9, family = "nbinom1")
summary(num.piptazo.vanc)
num.piptazo.vanc.model <- cbind(Estimate = coef(num.piptazo.vanc), confint(num.piptazo.vanc))
exp(num.piptazo.vanc.model)

rpkm.piptazo.carbapenem <- glmmadmb(carb.abund ~ piptazo.exp + (1|study_id), data = test9, family = "nbinom1")
summary(rpkm.piptazo.carbapenem)
rpkm.piptazo.carbapenem.model <- cbind(Estimate = coef(rpkm.piptazo.carbapenem), confint(rpkm.piptazo.carbapenem))
exp(rpkm.piptazo.carbapenem.model)

num.piptazo.carbapenem <- glmmadmb(num.carbapenem.ARG ~ piptazo.exp + (1|study_id), data = test9, family = "nbinom1")
summary(num.piptazo.carbapenem)
num.piptazo.carbapenem.model <- cbind(Estimate = coef(num.piptazo.carbapenem), confint(num.piptazo.carbapenem))
exp(num.piptazo.carbapenem.model)

rpkm.piptazo.quinolone <- glmmadmb(fluoro.abund ~ piptazo.exp + (1|study_id), data = test9, family = "nbinom1")
summary(rpkm.piptazo.quinolone)
rpkm.piptazo.quinolone.model <- cbind(Estimate = coef(rpkm.piptazo.quinolone), confint(rpkm.piptazo.quinolone))
exp(rpkm.piptazo.quinolone.model)

num.piptazo.quinolone <- glmmadmb(num.quinolone.ARG ~ piptazo.exp + (1|study_id), data = test9, family = "nbinom1")
summary(num.piptazo.quinolone)
num.piptazo.quinolone.model <- cbind(Estimate = coef(num.piptazo.quinolone), confint(num.piptazo.quinolone))
exp(num.piptazo.quinolone.model)

rpkm.piptazo.amino <- glmmadmb(amino.abund ~ piptazo.exp + (1|study_id), data = test9, family = "nbinom1")
summary(rpkm.piptazo.amino)
rpkm.piptazo.amino.model <- cbind(Estimate = coef(rpkm.piptazo.amino), confint(rpkm.piptazo.amino))
exp(rpkm.piptazo.amino.model)

num.piptazo.amino <- glmmadmb(num.amino.ARG ~ piptazo.exp + (1|study_id), data = test9, family = "nbinom1")
summary(num.piptazo.amino)
num.piptazo.amino.model <- cbind(Estimate = coef(num.piptazo.amino), confint(num.piptazo.amino))
exp(num.piptazo.amino.model)

rpkm.piptazo.clinda <- glmmadmb(clinda.abund ~ piptazo.exp + (1|study_id), data = test9, family = "nbinom1")
summary(rpkm.piptazo.clinda)
rpkm.piptazo.clinda.model <- cbind(Estimate = coef(rpkm.piptazo.clinda), confint(rpkm.piptazo.clinda))
exp(rpkm.piptazo.clinda.model)

num.piptazo.clinda <- glmmadmb(num.clinda.ARG ~ piptazo.exp + (1|study_id), data = test9, family = "nbinom2")
summary(num.piptazo.clinda)
num.piptazo.clinda.model <- cbind(Estimate = coef(num.piptazo.clinda), confint(num.piptazo.clinda))
exp(num.piptazo.clinda.model)

rpkm.piptazo.macrolide <- glmmadmb(macrolide.abund ~ piptazo.exp + (1|study_id), data = test9, family = "nbinom1")
summary(rpkm.piptazo.macrolide)
rpkm.piptazo.macrolide.model <- cbind(Estimate = coef(rpkm.piptazo.macrolide), confint(rpkm.piptazo.macrolide))
exp(rpkm.piptazo.macrolide.model)

num.piptazo.macrolide <- glmmadmb(num.macrolide.ARG ~ piptazo.exp + (1|study_id), data = test9, family = "nbinom1")
summary(num.piptazo.macrolide)
num.piptazo.macrolide.model <- cbind(Estimate = coef(num.piptazo.macrolide), confint(num.piptazo.macrolide))
exp(num.piptazo.macrolide.model)

rpkm.piptazo.trimeth <- glmmadmb(trimeth.abund ~ piptazo.exp + (1|study_id), data = test9, family = "nbinom1")
summary(rpkm.piptazo.trimeth)
rpkm.piptazo.trimeth.model <- cbind(Estimate = coef(rpkm.piptazo.trimeth), confint(rpkm.piptazo.trimeth))
exp(rpkm.piptazo.trimeth.model)

num.piptazo.trimeth <- glmmadmb(num.trimeth.ARG ~ piptazo.exp + (1|study_id), data = test9, family = "nbinom2")
summary(num.piptazo.trimeth)
num.piptazo.trimeth.model <- cbind(Estimate = coef(num.piptazo.trimeth), confint(num.piptazo.trimeth))
exp(num.piptazo.trimeth.model)

rpkm.piptazo.tetra <- glmmadmb(tetra.abund ~ piptazo.exp + (1|study_id), data = test9, family = "nbinom1")
rpkm.piptazo.tetra.model <- cbind(Estimate = coef(rpkm.piptazo.tetra), confint(rpkm.piptazo.tetra))
exp(rpkm.piptazo.tetra.model)
summary(rpkm.piptazo.tetra)

num.piptazo.tetra <- glmmadmb(num.tetra.ARG ~ piptazo.exp + (1|study_id), data = test9, family = "nbinom1")
num.piptazo.tetra.model <- cbind(Estimate = coef(num.piptazo.tetra), confint(num.piptazo.tetra))
exp(num.piptazo.tetra.model)
summary(num.piptazo.tetra)

rpkm.piptazo.flagyl <- glmmadmb(flagyl.abund ~ piptazo.exp + (1|study_id), data = test9, family = "nbinom1")
rpkm.piptazo.flagyl.model <- cbind(Estimate = coef(rpkm.piptazo.flagyl), confint(rpkm.piptazo.flagyl))
exp(rpkm.piptazo.flagyl.model)
summary(rpkm.piptazo.flagyl)

num.piptazo.flagyl <- glmmadmb(num.flagyl.ARG ~ piptazo.exp + (1|study_id), data = test9, family = "nbinom1")
num.piptazo.flagyl.model <- cbind(Estimate = coef(num.piptazo.flagyl), confint(num.piptazo.flagyl))
exp(num.piptazo.flagyl.model)
summary(num.piptazo.flagyl)

#adjusted P value for multiple comparisons
pvalues <- c(0.0024, 2e-16, 0.65, 0.00084, 2e-16, 0.0093, 0.0061, 0.25, 2e-16, 0.036, 6.5e-06, 0.30, 0.08, 0.0037, 0.0039, 0.054, 0.0033, 0.5, 1.5e-05,   0.099)

p.adjust(pvalues, method = "BH")
arg.class <- c("BetalactamAbund", "VancomycinAbund", "CarbapenemAbund", "QuinoloneAbund", "AminoglycosidAbund", "ClindamycinAbund", "MacrolideAbund", "TrimethoprimAbund", "TetracyclineAbund", "MetronidazoleAbund",
               "Betalactam.num", "Vancomycin.num", "Carbapenem.num", "Quinolone.num", "Aminoglycosid.num", "Clindamycin.num", "Macrolide.num", "Trimethoprim.num", "Tetracycline.num", "Metronidazole.num")
adj.p.piptazo<- c(6.857143e-03, 1.333333e-15, 6.500000e-01, 2.800000e-03, 1.333333e-15, 1.550000e-02, 1.109091e-02, 2.941176e-01,
                  1.333333e-15, 5.538462e-02, 3.250000e-05, 3.333333e-01, 1.066667e-01, 7.800000e-03, 7.800000e-03, 7.714286e-02,
                  7.800000e-03, 5.263158e-01, 6.000000e-05, 1.237500e-01)

piptazo.bh <- data.frame(adj.p.piptazo, rownames = arg.class) %>% column_to_rownames(var = "rownames")
piptazo.bh <- piptazo.bh %>% format(adj.p.piptazo, scientific = FALSE) 
piptazo.bh$adj.p.piptazo <- as.numeric(piptazo.bh$adj.p.piptazo)
piptazo.bh$adj.p.piptazo <- signif(piptazo.bh$adj.p.piptazo, 2)


#carbapenem exposure
rpkm.carbapenem.betalactam <- glmmadmb(beta.abund ~ carbapenem.exp + (1|study_id), data = test9, family = "nbinom1")
rpkm.carbapenem.betalactam.model <- cbind(Estimate = coef(rpkm.carbapenem.betalactam), confint(rpkm.carbapenem.betalactam))
exp(rpkm.carbapenem.betalactam.model)
summary(rpkm.carbapenem.betalactam)

num.carbapenem.betalactam <- glmmadmb(num.betalactam.ARG ~ carbapenem.exp + (1|study_id), data = test9, family = "nbinom1")
num.carbapenem.betalactam.model <- cbind(Estimate = coef(num.carbapenem.betalactam), confint(num.carbapenem.betalactam))
exp(num.carbapenem.betalactam.model)
summary(num.carbapenem.betalactam)

rpkm.carbapenem.vanc <- glmmadmb(vanc.abund ~ carbapenem.exp + (1|study_id), data = test9, family = "nbinom1")
rpkm.carbapenem.vanc.model <- cbind(Estimate = coef(rpkm.carbapenem.vanc), confint(rpkm.carbapenem.vanc))
exp(rpkm.carbapenem.vanc.model)
summary(rpkm.carbapenem.vanc)

num.carbapenem.vanc <- glmmadmb(num.vanc.ARG ~ carbapenem.exp + (1|study_id), data = test9, family = "nbinom1")
num.carbapenem.vanc.model <- cbind(Estimate = coef(num.carbapenem.vanc), confint(num.carbapenem.vanc))
exp(num.carbapenem.vanc.model)
summary(num.carbapenem.vanc)

rpkm.carbapenem.carbapenem <- glmmadmb(carb.abund ~ carbapenem.exp + (1|study_id), data = test9, family = "nbinom1")
rpkm.carbapenem.carbapenem.model <- cbind(Estimate = coef(rpkm.carbapenem.carbapenem), confint(rpkm.carbapenem.carbapenem))
exp(rpkm.carbapenem.carbapenem.model)
summary(rpkm.carbapenem.carbapenem)

num.carbapenem.carbapenem <- glmmadmb(num.carbapenem.ARG ~ carbapenem.exp + (1|study_id), data = test9, family = "nbinom1")
num.carbapenem.carbapenem.model <- cbind(Estimate = coef(num.carbapenem.carbapenem), confint(num.carbapenem.carbapenem))
exp(num.carbapenem.carbapenem.model)
summary(num.carbapenem.carbapenem)

rpkm.carbapenem.quinolone <- glmmadmb(fluoro.abund ~ carbapenem.exp + (1|study_id), data = test9, family = "nbinom1")
rpkm.carbapenem.quinolone.model <- cbind(Estimate = coef(rpkm.carbapenem.quinolone), confint(rpkm.carbapenem.quinolone))
exp(rpkm.carbapenem.quinolone.model)
summary(rpkm.carbapenem.quinolone)

num.carbapenem.quinolone <- glmmadmb(num.quinolone.ARG ~ carbapenem.exp + (1|study_id), data = test9, family = "nbinom1")
num.carbapenem.quinolone.model <- cbind(Estimate = coef(num.carbapenem.quinolone), confint(num.carbapenem.quinolone))
exp(num.carbapenem.quinolone.model)
summary(num.carbapenem.quinolone)

rpkm.carbapenem.amino <- glmmadmb(amino.abund ~ carbapenem.exp + (1|study_id), data = test9, family = "nbinom1")
rpkm.carbapenem.amino.model <- cbind(Estimate = coef(rpkm.carbapenem.amino), confint(rpkm.carbapenem.amino))
exp(rpkm.carbapenem.amino.model)
summary(rpkm.carbapenem.amino)

num.carbapenem.amino <- glmmadmb(num.amino.ARG ~ carbapenem.exp + (1|study_id), data = test9, family = "nbinom1")
num.carbapenem.amino.model <- cbind(Estimate = coef(num.carbapenem.amino), confint(num.carbapenem.amino))
exp(num.carbapenem.amino.model)
summary(num.carbapenem.amino)

rpkm.carbapenem.clinda <- glmmadmb(clinda.abund ~ carbapenem.exp + (1|study_id), data = test9, family = "nbinom1")
rpkm.carbapenem.clinda.model <- cbind(Estimate = coef(rpkm.carbapenem.clinda), confint(rpkm.carbapenem.clinda))
exp(rpkm.carbapenem.clinda.model)
summary(rpkm.carbapenem.clinda)

num.carbapenem.clinda <- glmmadmb(num.clinda.ARG ~ carbapenem.exp + (1|study_id), data = test9, family = "nbinom1")
num.carbapenem.clinda.model <- cbind(Estimate = coef(num.carbapenem.clinda), confint(num.carbapenem.clinda))
exp(num.carbapenem.clinda.model)
summary(num.carbapenem.clinda)

rpkm.carbapenem.macrolide <- glmmadmb(macrolide.abund ~ carbapenem.exp + (1|study_id), data = test9, family = "nbinom1")
rpkm.carbapenem.macrolide.model <- cbind(Estimate = coef(rpkm.carbapenem.macrolide), confint(rpkm.carbapenem.macrolide))
exp(rpkm.carbapenem.macrolide.model)
summary(rpkm.carbapenem.macrolide)

num.carbapenem.macrolide <- glmmadmb(num.macrolide.ARG ~ carbapenem.exp + (1|study_id), data = test9, family = "nbinom1")
num.carbapenem.macrolide.model <- cbind(Estimate = coef(num.carbapenem.macrolide), confint(num.carbapenem.macrolide))
exp(num.carbapenem.macrolide.model)
summary(num.carbapenem.macrolide)

rpkm.carbapenem.trimeth <- glmmadmb(trimeth.abund ~ carbapenem.exp + (1|study_id), data = test9, family = "nbinom1")
rpkm.carbapenem.trimeth.model <- cbind(Estimate = coef(rpkm.carbapenem.trimeth), confint(rpkm.carbapenem.trimeth))
exp(rpkm.carbapenem.trimeth.model)
summary(rpkm.carbapenem.trimeth)

num.carbapenem.trimeth <- glmmadmb(num.trimeth.ARG ~ carbapenem.exp + (1|study_id), data = test9, family = "nbinom2")
num.carbapenem.trimeth.model <- cbind(Estimate = coef(num.carbapenem.trimeth), confint(num.carbapenem.trimeth))
exp(num.carbapenem.trimeth.model)
summary(num.carbapenem.trimeth)

rpkm.carbapenem.tetra <- glmmadmb(tetra.abund ~ carbapenem.exp + (1|study_id), data = test9, family = "nbinom1")
rpkm.carbapenem.tetra.model <- cbind(Estimate = coef(rpkm.carbapenem.tetra), confint(rpkm.carbapenem.tetra))
exp(rpkm.carbapenem.tetra.model)
summary(rpkm.carbapenem.tetra)

num.carbapenem.tetra <- glmmadmb(num.tetra.ARG ~ carbapenem.exp + (1|study_id), data = test9, family = "nbinom1")
num.carbapenem.tetra.model <- cbind(Estimate = coef(num.carbapenem.tetra), confint(num.carbapenem.tetra))
exp(num.carbapenem.tetra.model)
summary(num.carbapenem.tetra)

rpkm.carbapenem.flagyl <- glmmadmb(flagyl.abund ~ carbapenem.exp + (1|study_id), data = test9, family = "nbinom1")
rpkm.carbapenem.flagyl.model <- cbind(Estimate = coef(rpkm.carbapenem.flagyl), confint(rpkm.carbapenem.flagyl))
exp(rpkm.carbapenem.flagyl.model)
summary(rpkm.carbapenem.flagyl)
#did not converge
#num.carbapenem.flagyl <- glmmadmb(num.flagyl.ARG ~ carbapenem.exp + (1|study_id), data = test9, family = "nbinom")
#num.carbapenem.flagyl.model <- cbind(Estimate = coef(num.carbapenem.flagyl), confint(num.carbapenem.flagyl))
#exp(num.carbapenem.flagyl.model)
#summary(num.carbapenem.flagyl)

#adjusted P value for multiple comparisons
pvalues <- c(0.03, 0.77, 0.6, 0.0064, 0.33, 0.16, 0.016, 0.1, 0.92, 0.74, 0.27, 0.25, 0.85, 0.76, 0.96, 0.8, 0.91, 0.074, 0.17, NA)

p.adjust(pvalues, method = "BH", n=20)
arg.class <- c("BetalactamAbund", "VancomycinAbund", "CarbapenemAbund", "QuinoloneAbund", "AminoglycosidAbund", "ClindamycinAbund", "MacrolideAbund", "TrimethoprimAbund", "TetracyclineAbund", "MetronidazoleAbund",
               "Betalactam.num", "Vancomycin.num", "Carbapenem.num", "Quinolone.num", "Aminoglycosid.num", "Clindamycin.num", "Macrolide.num", "Trimethoprim.num", "Tetracycline.num", "Metronidazole.num")
adj.p.carbapenem<- c(
  0.2000000, 1.0000000, 1.0000000, 0.1280000, 0.6600000, 0.4857143, 0.1600000, 0.4000000, 1.0000000, 1.0000000, 
  0.6000000, 0.6000000, 1.0000000, 1.0000000, 1.0000000, 1.0000000, 1.0000000, 0.3700000, 0.4857143,        NA)

carbapenem.bh <- data.frame(adj.p.carbapenem, rownames = arg.class) %>% column_to_rownames(var = "rownames")
carbapenem.bh <- carbapenem.bh %>% format(adj.p.carbapenem, scientific = FALSE) 
carbapenem.bh$adj.p.carbapenem <- as.numeric(carbapenem.bh$adj.p.carbapenem)
carbapenem.bh$adj.p.carbapenem <- signif(carbapenem.bh$adj.p.carbapenem, 2)

#vancomycin exposure
rpkm.vancomycin.betalactam <- glmmadmb(beta.abund ~ vanc.exp + (1|study_id), data = test9, family = "nbinom1")
rpkm.vancomycin.betalactam.model <- cbind(Estimate = coef(rpkm.vancomycin.betalactam), confint(rpkm.vancomycin.betalactam))
exp(rpkm.vancomycin.betalactam.model)
summary(rpkm.vancomycin.betalactam)

num.vancomycin.betalactam <- glmmadmb(num.betalactam.ARG ~ vanc.exp + (1|study_id), data = test9, family = "nbinom1")
num.vancomycin.betalactam.model <- cbind(Estimate = coef(num.vancomycin.betalactam), confint(num.vancomycin.betalactam))
exp(num.vancomycin.betalactam.model)
summary(num.vancomycin.betalactam)

rpkm.vancomycin.vanc <- glmmadmb(vanc.abund ~ vanc.exp + (1|study_id), data = test9, family = "nbinom1")
rpkm.vancomycin.vanc.model <- cbind(Estimate = coef(rpkm.vancomycin.vanc), confint(rpkm.vancomycin.vanc))
exp(rpkm.vancomycin.vanc.model)
summary(rpkm.vancomycin.vanc)

num.vancomycin.vanc <- glmmadmb(num.vanc.ARG ~ vanc.exp + (1|study_id), data = test9, family = "nbinom1")
num.vancomycin.vanc.model <- cbind(Estimate = coef(num.vancomycin.vanc), confint(num.vancomycin.vanc))
exp(num.vancomycin.vanc.model)
summary(num.vancomycin.vanc)

rpkm.vancomycin.vancomycin <- glmmadmb(carb.abund ~ vanc.exp + (1|study_id), data = test9, family = "nbinom1")
rpkm.vancomycin.vancomycin.model <- cbind(Estimate = coef(rpkm.vancomycin.vancomycin), confint(rpkm.vancomycin.vancomycin))
exp(rpkm.vancomycin.vancomycin.model)
summary(rpkm.vancomycin.vancomycin)

num.vancomycin.vancomycin <- glmmadmb(num.carbapenem.ARG ~ vanc.exp + (1|study_id), data = test9, family = "nbinom1")
num.vancomycin.vancomycin.model <- cbind(Estimate = coef(num.vancomycin.vancomycin), confint(num.vancomycin.vancomycin))
exp(num.vancomycin.vancomycin.model)
summary(num.vancomycin.vancomycin)

rpkm.vancomycin.quinolone <- glmmadmb(fluoro.abund ~ vanc.exp + (1|study_id), data = test9, family = "nbinom2")
rpkm.vancomycin.quinolone.model <- cbind(Estimate = coef(rpkm.vancomycin.quinolone), confint(rpkm.vancomycin.quinolone))
exp(rpkm.vancomycin.quinolone.model)
summary(rpkm.vancomycin.quinolone)

num.vancomycin.quinolone <- glmmadmb(num.quinolone.ARG ~ vanc.exp + (1|study_id), data = test9, family = "nbinom1")
num.vancomycin.quinolone.model <- cbind(Estimate = coef(num.vancomycin.quinolone), confint(num.vancomycin.quinolone))
exp(num.vancomycin.quinolone.model)
summary(num.vancomycin.quinolone)

rpkm.vancomycin.amino <- glmmadmb(amino.abund ~ vanc.exp + (1|study_id), data = test9, family = "nbinom1")
rpkm.vancomycin.amino.model <- cbind(Estimate = coef(rpkm.vancomycin.amino), confint(rpkm.vancomycin.amino))
exp(rpkm.vancomycin.amino.model)
summary(rpkm.vancomycin.amino)

num.vancomycin.amino <- glmmadmb(num.amino.ARG ~ vanc.exp + (1|study_id), data = test9, family = "nbinom1")
num.vancomycin.amino.model <- cbind(Estimate = coef(num.vancomycin.amino), confint(num.vancomycin.amino))
exp(num.vancomycin.amino.model)
summary(num.vancomycin.amino)

rpkm.vancomycin.clinda <- glmmadmb(clinda.abund ~ vanc.exp + (1|study_id), data = test9, family = "nbinom1")
rpkm.vancomycin.clinda.model <- cbind(Estimate = coef(rpkm.vancomycin.clinda), confint(rpkm.vancomycin.clinda))
exp(rpkm.vancomycin.clinda.model)
summary(rpkm.vancomycin.clinda)

num.vancomycin.clinda <- glmmadmb(num.clinda.ARG ~ vanc.exp + (1|study_id), data = test9, family = "nbinom1")
num.vancomycin.clinda.model <- cbind(Estimate = coef(num.vancomycin.clinda), confint(num.vancomycin.clinda))
exp(num.vancomycin.clinda.model)
summary(num.vancomycin.clinda)

rpkm.vancomycin.macrolide <- glmmadmb(macrolide.abund ~ vanc.exp + (1|study_id), data = test9, family = "nbinom1")
rpkm.vancomycin.macrolide.model <- cbind(Estimate = coef(rpkm.vancomycin.macrolide), confint(rpkm.vancomycin.macrolide))
exp(rpkm.vancomycin.macrolide.model)
summary(rpkm.vancomycin.macrolide)

num.vancomycin.macrolide <- glmmadmb(num.macrolide.ARG ~ vanc.exp + (1|study_id), data = test9, family = "nbinom1")
num.vancomycin.macrolide.model <- cbind(Estimate = coef(num.vancomycin.macrolide), confint(num.vancomycin.macrolide))
exp(num.vancomycin.macrolide.model)
summary(num.vancomycin.macrolide)

rpkm.vancomycin.trimeth <- glmmadmb(trimeth.abund ~ vanc.exp + (1|study_id), data = test9, family = "nbinom1")
rpkm.vancomycin.trimeth.model <- cbind(Estimate = coef(rpkm.vancomycin.trimeth), confint(rpkm.vancomycin.trimeth))
exp(rpkm.vancomycin.trimeth.model)
summary(rpkm.vancomycin.trimeth)

num.vancomycin.trimeth <- glmmadmb(num.trimeth.ARG ~ vanc.exp + (1|study_id), data = test9, family = "nbinom2")
num.vancomycin.trimeth.model <- cbind(Estimate = coef(num.vancomycin.trimeth), confint(num.vancomycin.trimeth))
exp(num.vancomycin.trimeth.model)
summary(num.vancomycin.trimeth)

rpkm.vancomycin.tetra <- glmmadmb(tetra.abund ~ vanc.exp + (1|study_id), data = test9, family = "nbinom1")
rpkm.vancomycin.tetra.model <- cbind(Estimate = coef(rpkm.vancomycin.tetra), confint(rpkm.vancomycin.tetra))
exp(rpkm.vancomycin.tetra.model)
summary(rpkm.vancomycin.tetra)

num.vancomycin.tetra <- glmmadmb(num.tetra.ARG ~ vanc.exp + (1|study_id), data = test9, family = "nbinom1")
num.vancomycin.tetra.model <- cbind(Estimate = coef(num.vancomycin.tetra), confint(num.vancomycin.tetra))
exp(num.vancomycin.tetra.model)
summary(num.vancomycin.tetra)

rpkm.vancomycin.flagyl <- glmmadmb(flagyl.abund ~ vanc.exp + (1|study_id), data = test9, family = "nbinom1")
rpkm.vancomycin.flagyl.model <- cbind(Estimate = coef(rpkm.vancomycin.flagyl), confint(rpkm.vancomycin.flagyl))
exp(rpkm.vancomycin.flagyl.model)
summary(rpkm.vancomycin.flagyl)

num.vancomycin.flagyl <- glmmadmb(num.flagyl.ARG ~ vanc.exp + (1|study_id), data = test9, family = "nbinom2")
num.vancomycin.flagyl.model <- cbind(Estimate = coef(num.vancomycin.flagyl), confint(num.vancomycin.flagyl))
exp(num.vancomycin.flagyl.model)
summary(num.vancomycin.flagyl)

#adjusted P value for multiple comparisons
pvalues <- c(0.035, 0.57, 0.17, 0.041, 7.7e-07, 0.37, 0.0085, 0.015, 8.3e-05, 0.00057, 7.4e-06, 0.57, 0.0023, 2.4e-05, 1.3e-06, 1.1e-05, 3.3e-05, 0.0097, 1.3e-09, 0.0025)

p.adjust(pvalues, method = "BH", n=20)
arg.class <- c("BetalactamAbund", "VancomycinAbund", "vancomycinAbund", "QuinoloneAbund", "AminoglycosidAbund", "ClindamycinAbund", "MacrolideAbund", "TrimethoprimAbund", "TetracyclineAbund", "MetronidazoleAbund",
               "Betalactam.num", "Vancomycin.num", "vancomycin.num", "Quinolone.num", "Aminoglycosid.num", "Clindamycin.num", "Macrolide.num", "Trimethoprim.num", "Tetracycline.num", "Metronidazole.num")
adj.p.vancomycin<- c(
  4.666667e-02, 5.700000e-01, 2.000000e-01, 5.125000e-02, 7.700000e-06, 4.111111e-01, 1.416667e-02, 2.142857e-02,
  2.075000e-04, 1.266667e-03, 3.700000e-05, 5.700000e-01, 4.545455e-03, 8.000000e-05, 8.666667e-06, 4.400000e-05,
  9.428571e-05, 1.492308e-02, 2.600000e-08, 4.545455e-03)

vancomycin.bh <- data.frame(adj.p.vancomycin, rownames = arg.class) %>% column_to_rownames(var = "rownames")
vancomycin.bh$adj.p.vancomycin <- as.numeric(vancomycin.bh$adj.p.vancomycin)
vancomycin.bh$adj.p.vancomycin <- signif(vancomycin.bh$adj.p.vancomycin, 2) %>% format(adj.p.vancomycin, scientific = FALSE) 

#flagyl exposure
rpkm.flagyl.betalactam <- glmmadmb(beta.abund ~ flagyl.exp + (1|study_id), data = test9, family = "nbinom1")
rpkm.flagyl.betalactam.model <- cbind(Estimate = coef(rpkm.flagyl.betalactam), confint(rpkm.flagyl.betalactam))
exp(rpkm.flagyl.betalactam.model)
summary(rpkm.flagyl.betalactam)

num.flagyl.betalactam <- glmmadmb(num.betalactam.ARG ~ flagyl.exp + (1|study_id), data = test9, family = "nbinom1")
num.flagyl.betalactam.model <- cbind(Estimate = coef(num.flagyl.betalactam), confint(num.flagyl.betalactam))
exp(num.flagyl.betalactam.model)
summary(num.flagyl.betalactam)

rpkm.flagyl.flagyl <- glmmadmb(vanc.abund ~ flagyl.exp + (1|study_id), data = test9, family = "nbinom1")
rpkm.flagyl.flagyl.model <- cbind(Estimate = coef(rpkm.flagyl.flagyl), confint(rpkm.flagyl.flagyl))
exp(rpkm.flagyl.flagyl.model)
summary(rpkm.flagyl.flagyl)

num.flagyl.flagyl <- glmmadmb(num.vanc.ARG ~ flagyl.exp + (1|study_id), data = test9, family = "nbinom1")
num.flagyl.flagyl.model <- cbind(Estimate = coef(num.flagyl.flagyl), confint(num.flagyl.flagyl))
exp(num.flagyl.flagyl.model)
summary(num.flagyl.flagyl)

rpkm.flagyl.flagyl <- glmmadmb(carb.abund ~ flagyl.exp + (1|study_id), data = test9, family = "nbinom1")
rpkm.flagyl.flagyl.model <- cbind(Estimate = coef(rpkm.flagyl.flagyl), confint(rpkm.flagyl.flagyl))
exp(rpkm.flagyl.flagyl.model)
summary(rpkm.flagyl.flagyl)

num.flagyl.flagyl <- glmmadmb(num.carbapenem.ARG ~ flagyl.exp + (1|study_id), data = test9, family = "nbinom1")
num.flagyl.flagyl.model <- cbind(Estimate = coef(num.flagyl.flagyl), confint(num.flagyl.flagyl))
exp(num.flagyl.flagyl.model)
summary(num.flagyl.flagyl)

rpkm.flagyl.quinolone <- glmmadmb(fluoro.abund ~ flagyl.exp + (1|study_id), data = test9, family = "nbinom1")
rpkm.flagyl.quinolone.model <- cbind(Estimate = coef(rpkm.flagyl.quinolone), confint(rpkm.flagyl.quinolone))
exp(rpkm.flagyl.quinolone.model)
summary(rpkm.flagyl.quinolone)

num.flagyl.quinolone <- glmmadmb(num.quinolone.ARG ~ flagyl.exp + (1|study_id), data = test9, family = "nbinom1")
num.flagyl.quinolone.model <- cbind(Estimate = coef(num.flagyl.quinolone), confint(num.flagyl.quinolone))
exp(num.flagyl.quinolone.model)
summary(num.flagyl.quinolone)

rpkm.flagyl.amino <- glmmadmb(amino.abund ~ flagyl.exp + (1|study_id), data = test9, family = "nbinom1")
rpkm.flagyl.amino.model <- cbind(Estimate = coef(rpkm.flagyl.amino), confint(rpkm.flagyl.amino))
exp(rpkm.flagyl.amino.model)
summary(rpkm.flagyl.amino)

num.flagyl.amino <- glmmadmb(num.amino.ARG ~ flagyl.exp + (1|study_id), data = test9, family = "nbinom1")
num.flagyl.amino.model <- cbind(Estimate = coef(num.flagyl.amino), confint(num.flagyl.amino))
exp(num.flagyl.amino.model)
summary(num.flagyl.amino)

rpkm.flagyl.clinda <- glmmadmb(clinda.abund ~ flagyl.exp + (1|study_id), data = test9, family = "nbinom1")
rpkm.flagyl.clinda.model <- cbind(Estimate = coef(rpkm.flagyl.clinda), confint(rpkm.flagyl.clinda))
exp(rpkm.flagyl.clinda.model)
summary(rpkm.other.abx.clinda)

num.flagyl.clinda <- glmmadmb(num.clinda.ARG ~ flagyl.exp + (1|study_id), data = test9, family = "nbinom1")
num.flagyl.clinda.model <- cbind(Estimate = coef(num.flagyl.clinda), confint(num.flagyl.clinda))
exp(num.flagyl.clinda.model)
summary(num.flagyl.clinda)

rpkm.flagyl.macrolide <- glmmadmb(macrolide.abund ~ flagyl.exp + (1|study_id), data = test9, family = "nbinom1")
rpkm.flagyl.macrolide.model <- cbind(Estimate = coef(rpkm.flagyl.macrolide), confint(rpkm.flagyl.macrolide))
exp(rpkm.flagyl.macrolide.model)
summary(rpkm.flagyl.macrolide)

num.flagyl.macrolide <- glmmadmb(num.macrolide.ARG ~ flagyl.exp + (1|study_id), data = test9, family = "nbinom1")
num.flagyl.macrolide.model <- cbind(Estimate = coef(num.flagyl.macrolide), confint(num.flagyl.macrolide))
exp(num.flagyl.macrolide.model)
summary(num.flagyl.macrolide)

rpkm.flagyl.trimeth <- glmmadmb(trimeth.abund ~ flagyl.exp + (1|study_id), data = test9, family = "nbinom1")
rpkm.flagyl.trimeth.model <- cbind(Estimate = coef(rpkm.flagyl.trimeth), confint(rpkm.flagyl.trimeth))
exp(rpkm.flagyl.trimeth.model)
summary(rpkm.flagyl.trimeth)

num.flagyl.trimeth <- glmmadmb(num.trimeth.ARG ~ flagyl.exp + (1|study_id), data = test9, family = "nbinom2")
num.flagyl.trimeth.model <- cbind(Estimate = coef(num.flagyl.trimeth), confint(num.flagyl.trimeth))
exp(num.flagyl.trimeth.model)
summary(num.flagyl.trimeth)

rpkm.flagyl.tetra <- glmmadmb(tetra.abund ~ flagyl.exp + (1|study_id), data = test9, family = "nbinom1")
rpkm.flagyl.tetra.model <- cbind(Estimate = coef(rpkm.flagyl.tetra), confint(rpkm.flagyl.tetra))
exp(rpkm.flagyl.tetra.model)
summary(rpkm.flagyl.tetra)

num.flagyl.tetra <- glmmadmb(num.tetra.ARG ~ flagyl.exp + (1|study_id), data = test9, family = "nbinom1")
num.flagyl.tetra.model <- cbind(Estimate = coef(num.flagyl.tetra), confint(num.flagyl.tetra))
exp(num.flagyl.tetra.model)
summary(num.flagyl.tetra)

rpkm.flagyl.flagyl <- glmmadmb(flagyl.abund ~ flagyl.exp + (1|study_id), data = test9, family = "nbinom1")
rpkm.flagyl.flagyl.model <- cbind(Estimate = coef(rpkm.flagyl.flagyl), confint(rpkm.flagyl.flagyl))
exp(rpkm.flagyl.flagyl.model)
summary(rpkm.flagyl.flagyl)

num.flagyl.flagyl <- glmmadmb(num.flagyl.ARG ~ flagyl.exp + (1|study_id), data = test9, family = "nbinom", zeroInflation = TRUE)
num.flagyl.flagyl.model <- cbind(Estimate = coef(num.flagyl.flagyl), confint(num.flagyl.flagyl))
exp(num.flagyl.flagyl.model)
summary(num.flagyl.flagyl)

#adjusted P value for multiple comparisons
pvalues <- c(  0.34, 0.0062, 0.33, 0.066, 4e-04, 0.088, 0.0021, 0.1, 0.12, 0.2, 0.00044, 0.08, 0.16, 0.14, 0.027, 0.001, 0.011, 0.69, 3.4e-05, 0.15)
p.adjust(pvalues, method = "BH", n=20)
arg.class <- c("BetalactamAbund", "vancAbund", "flagylAbund", "QuinoloneAbund", "AminoglycosidAbund", "ClindamycinAbund", "MacrolideAbund", "TrimethoprimAbund", "TetracyclineAbund", "MetronidazoleAbund",
               "Betalactam.num", "vanc.num", "flagyl.num", "Quinolone.num", "Aminoglycosid.num", "Clindamycin.num", "Macrolide.num", "Trimethoprim.num", "Tetracycline.num", "Metronidazole.num")
adj.p.flagyl<- c( 
  0.357894737, 0.020666667, 0.357894737, 0.146666667,
  0.002933333, 0.160000000, 0.008400000, 0.166666667,
  0.184615385, 0.235294118, 0.002933333, 0.160000000,
  0.200000000, 0.200000000, 0.067500000, 0.005000000,
  0.031428571, 0.690000000, 0.000680000, 0.200000000)
flagyl.bh <- data.frame(adj.p.flagyl, rownames = arg.class) %>% column_to_rownames(var = "rownames")
flagyl.bh <- flagyl.bh %>% format(adj.p.flagyl, scientific = FALSE) 
flagyl.bh$adj.p.flagyl <- as.numeric(flagyl.bh$adj.p.flagyl)
flagyl.bh$adj.p.flagyl <- signif(flagyl.bh$adj.p.flagyl, 2)

#quinolone exposure
rpkm.quinolone.betalactam <- glmmadmb(beta.abund ~ quinolone.exp + (1|study_id), data = test9, family = "nbinom1")
rpkm.quinolone.betalactam.model <- cbind(Estimate = coef(rpkm.quinolone.betalactam), confint(rpkm.quinolone.betalactam))
exp(rpkm.quinolone.betalactam.model)
summary(rpkm.quinolone.betalactam)

num.quinolone.betalactam <- glmmadmb(num.betalactam.ARG ~ quinolone.exp + (1|study_id), data = test9, family = "nbinom1")
num.quinolone.betalactam.model <- cbind(Estimate = coef(num.quinolone.betalactam), confint(num.quinolone.betalactam))
exp(num.quinolone.betalactam.model)
summary(num.quinolone.betalactam)

rpkm.quinolone.vanc <- glmmadmb(vanc.abund ~ quinolone.exp + (1|study_id), data = test9, family = "nbinom1")
rpkm.quinolone.vanc.model <- cbind(Estimate = coef(rpkm.quinolone.vanc), confint(rpkm.quinolone.vanc))
exp(rpkm.quinolone.vanc.model)
summary(rpkm.quinolone.vanc)

num.quinolone.vanc <- glmmadmb(num.vanc.ARG ~ quinolone.exp + (1|study_id), data = test9, family = "nbinom1")
num.quinolone.vanc.model <- cbind(Estimate = coef(num.quinolone.vanc), confint(num.quinolone.vanc))
exp(num.quinolone.vanc.model)
summary(num.quinolone.vanc)

rpkm.quinolone.carb <- glmmadmb(carb.abund ~ quinolone.exp + (1|study_id), data = test9, family = "nbinom2")
rpkm.quinolone.carb.model <- cbind(Estimate = coef(rpkm.quinolone.carb), confint(rpkm.quinolone.carb))
exp(rpkm.quinolone.carb.model)
summary(rpkm.quinolone.carb)

num.quinolone.carb <- glmmadmb(num.carbapenem.ARG ~ quinolone.exp + (1|study_id), data = test9, family = "nbinom1")
num.quinolone.carb.model <- cbind(Estimate = coef(num.quinolone.carb), confint(num.quinolone.carb))
exp(num.quinolone.carb.model)
summary(num.quinolone.carb)

rpkm.quinolone.quinolone <- glmmadmb(fluoro.abund ~ quinolone.exp + (1|study_id), data = test9, family = "nbinom1")
rpkm.quinolone.quinolone.model <- cbind(Estimate = coef(rpkm.quinolone.quinolone), confint(rpkm.quinolone.quinolone))
exp(rpkm.quinolone.quinolone.model)
summary(rpkm.quinolone.quinolone)

num.quinolone.quinolone <- glmmadmb(num.quinolone.ARG ~ quinolone.exp + (1|study_id), data = test9, family = "nbinom1")
num.quinolone.quinolone.model <- cbind(Estimate = coef(num.quinolone.quinolone), confint(num.quinolone.quinolone))
exp(num.quinolone.quinolone.model)
summary(num.quinolone.quinolone)

rpkm.quinolone.amino <- glmmadmb(amino.abund ~ quinolone.exp + (1|study_id), data = test9, family = "nbinom1")
rpkm.quinolone.amino.model <- cbind(Estimate = coef(rpkm.quinolone.amino), confint(rpkm.quinolone.amino))
exp(rpkm.quinolone.amino.model)
summary(rpkm.quinolone.amino)

num.quinolone.amino <- glmmadmb(num.amino.ARG ~ quinolone.exp + (1|study_id), data = test9, family = "nbinom1")
num.quinolone.amino.model <- cbind(Estimate = coef(num.quinolone.amino), confint(num.quinolone.amino))
exp(num.quinolone.amino.model)
summary(num.quinolone.amino)

rpkm.quinolone.clinda <- glmmadmb(clinda.abund ~ quinolone.exp + (1|study_id), data = test9, family = "nbinom1")
rpkm.quinolone.clinda.model <- cbind(Estimate = coef(rpkm.quinolone.clinda), confint(rpkm.quinolone.clinda))
exp(rpkm.quinolone.clinda.model)
summary(rpkm.quinolone.clinda)

num.quinolone.clinda <- glmmadmb(num.clinda.ARG ~ quinolone.exp + (1|study_id), data = test9, family = "nbinom1")
num.quinolone.clinda.model <- cbind(Estimate = coef(num.quinolone.clinda), confint(num.quinolone.clinda))
exp(num.quinolone.clinda.model)
summary(num.quinolone.clinda)

rpkm.quinolone.macrolide <- glmmadmb(macrolide.abund ~ quinolone.exp + (1|study_id), data = test9, family = "nbinom1")
rpkm.quinolone.macrolide.model <- cbind(Estimate = coef(rpkm.quinolone.macrolide), confint(rpkm.quinolone.macrolide))
exp(rpkm.quinolone.macrolide.model)
summary(rpkm.quinolone.macrolide)

num.quinolone.macrolide <- glmmadmb(num.macrolide.ARG ~ quinolone.exp + (1|study_id), data = test9, family = "nbinom1")
num.quinolone.macrolide.model <- cbind(Estimate = coef(num.quinolone.macrolide), confint(num.quinolone.macrolide))
exp(num.quinolone.macrolide.model)
summary(num.quinolone.macrolide)

rpkm.quinolone.trimeth <- glmmadmb(trimeth.abund ~ quinolone.exp + (1|study_id), data = test9, family = "nbinom1")
rpkm.quinolone.trimeth.model <- cbind(Estimate = coef(rpkm.quinolone.trimeth), confint(rpkm.quinolone.trimeth))
exp(rpkm.quinolone.trimeth.model)
summary(rpkm.quinolone.trimeth)

num.quinolone.trimeth <- glmmadmb(num.trimeth.ARG ~ quinolone.exp + (1|study_id), data = test9, family = "nbinom1")
num.quinolone.trimeth.model <- cbind(Estimate = coef(num.quinolone.trimeth), confint(num.quinolone.trimeth))
exp(num.quinolone.trimeth.model)
summary(num.quinolone.trimeth)

rpkm.quinolone.tetra <- glmmadmb(tetra.abund ~ quinolone.exp + (1|study_id), data = test9, family = "nbinom1")
rpkm.quinolone.tetra.model <- cbind(Estimate = coef(rpkm.quinolone.tetra), confint(rpkm.quinolone.tetra))
exp(rpkm.quinolone.tetra.model)
summary(rpkm.quinolone.tetra)

num.quinolone.tetra <- glmmadmb(num.tetra.ARG ~ quinolone.exp + (1|study_id), data = test9, family = "nbinom1")
num.quinolone.tetra.model <- cbind(Estimate = coef(num.quinolone.tetra), confint(num.quinolone.tetra))
exp(num.quinolone.tetra.model)
summary(num.quinolone.tetra)

rpkm.quinolone.flagyl <- glmmadmb(flagyl.abund ~ quinolone.exp + (1|study_id), data = test9, family = "nbinom2")
rpkm.quinolone.flagyl.model <- cbind(Estimate = coef(rpkm.quinolone.flagyl), confint(rpkm.quinolone.flagyl))
exp(rpkm.quinolone.flagyl.model)
summary(rpkm.quinolone.flagyl)

num.quinolone.flagyl <- glmmadmb(num.flagyl.ARG ~ quinolone.exp + (1|study_id), data = test9, family = "nbinom2", zeroInflation = TRUE)
num.quinolone.flagyl.model <- cbind(Estimate = coef(num.quinolone.flagyl), confint(num.quinolone.flagyl))
exp(num.quinolone.flagyl.model)
summary(num.quinolone.flagyl)

#adjusted P value for multiple comparisons
pvalues <- c(0.077, 0.75, 0.6, 0.23, 0.13, 0.017, 0.036, 0.15, 0.67, 0.74, 0.16, 0.45, 0.41, 0.068, 0.17, 0.12, 0.095, 0.064, 0.066, 0.85)
p.adjust(pvalues, method = "BH", n=20)
arg.class <- c("BetalactamAbund", "vancAbund", "quinoloneAbund", "QuinoloneAbund", "AminoglycosidAbund", "ClindamycinAbund", "MacrolideAbund", "TrimethoprimAbund", "TetracyclineAbund", "MetronidazoleAbund",
               "Betalactam.num", "vanc.num", "quinolone.num", "Quinolone.num", "Aminoglycosid.num", "Clindamycin.num", "Macrolide.num", "Trimethoprim.num", "Tetracycline.num", "Metronidazole.num")
adj.p.quinolone<- c(
  0.2566667, 0.7894737, 0.7500000, 0.3538462, 0.2833333,
  0.2566667, 0.2566667, 0.2833333, 0.7882353, 0.7894737,
  0.2833333, 0.6000000, 0.5857143, 0.2566667, 0.2833333,
  0.2833333, 0.2714286, 0.2566667, 0.2566667, 0.8500000)
quinolone.bh <- data.frame(adj.p.quinolone, rownames = arg.class) %>% column_to_rownames(var = "rownames")
quinolone.bh <- quinolone.bh %>% format(adj.p.quinolone, scientific = FALSE) 
quinolone.bh$adj.p.quinolone <- as.numeric(quinolone.bh$adj.p.quinolone)
quinolone.bh$adj.p.quinolone <- signif(quinolone.bh$adj.p.quinolone, 2)

#trimethoprim exposure
rpkm.trimeth.betalactam <- glmmadmb(beta.abund ~ tmp_smx.exp + (1|study_id), data = test9, family = "nbinom1")
rpkm.trimeth.betalactam.model <- cbind(Estimate = coef(rpkm.trimeth.betalactam), confint(rpkm.trimeth.betalactam))
exp(rpkm.trimeth.betalactam.model)
summary(rpkm.trimeth.betalactam)

num.trimeth.betalactam <- glmmadmb(num.betalactam.ARG ~ tmp_smx.exp + (1|study_id), data = test9, family = "nbinom1")
num.trimeth.betalactam.model <- cbind(Estimate = coef(num.trimeth.betalactam), confint(num.trimeth.betalactam))
exp(num.trimeth.betalactam.model)
summary(num.trimeth.betalactam)

rpkm.trimeth.vanc <- glmmadmb(vanc.abund ~ tmp_smx.exp + (1|study_id), data = test9, family = "nbinom1")
rpkm.trimeth.vanc.model <- cbind(Estimate = coef(rpkm.trimeth.vanc), confint(rpkm.trimeth.vanc))
exp(rpkm.trimeth.vanc.model)
summary(rpkm.trimeth.vanc)

num.trimeth.vanc <- glmmadmb(num.vanc.ARG ~ tmp_smx.exp + (1|study_id), data = test9, family = "nbinom1")
num.trimeth.vanc.model <- cbind(Estimate = coef(num.trimeth.vanc), confint(num.trimeth.vanc))
exp(num.trimeth.vanc.model)
summary(num.trimeth.vanc)

rpkm.trimeth.carb <- glmmadmb(carb.abund ~ tmp_smx.exp + (1|study_id), data = test9, family = "nbinom1")
rpkm.trimeth.carb.model <- cbind(Estimate = coef(rpkm.trimeth.carb), confint(rpkm.trimeth.carb))
exp(rpkm.trimeth.carb.model)
summary(rpkm.trimeth.carb)

num.trimeth.carb <- glmmadmb(num.carbapenem.ARG ~ tmp_smx.exp + (1|study_id), data = test9, family = "nbinom1")
num.trimeth.carb.model <- cbind(Estimate = coef(num.trimeth.carb), confint(num.trimeth.carb))
exp(num.trimeth.carb.model)
summary(num.trimeth.carb)

rpkm.trimeth.trimeth <- glmmadmb(fluoro.abund ~ tmp_smx.exp + (1|study_id), data = test9, family = "nbinom2")
rpkm.trimeth.trimeth.model <- cbind(Estimate = coef(rpkm.trimeth.trimeth), confint(rpkm.trimeth.trimeth))
exp(rpkm.trimeth.trimeth.model)
summary(rpkm.trimeth.trimeth)

num.trimeth.trimeth <- glmmadmb(num.quinolone.ARG ~ tmp_smx.exp + (1|study_id), data = test9, family = "nbinom1")
num.trimeth.trimeth.model <- cbind(Estimate = coef(num.trimeth.trimeth), confint(num.trimeth.trimeth))
exp(num.trimeth.trimeth.model)
summary(num.trimeth.trimeth)

rpkm.trimeth.amino <- glmmadmb(amino.abund ~ tmp_smx.exp + (1|study_id), data = test9, family = "nbinom1")
rpkm.trimeth.amino.model <- cbind(Estimate = coef(rpkm.trimeth.amino), confint(rpkm.trimeth.amino))
exp(rpkm.trimeth.amino.model)
summary(rpkm.trimeth.amino)

num.trimeth.amino <- glmmadmb(num.amino.ARG ~ tmp_smx.exp + (1|study_id), data = test9, family = "nbinom1")
num.trimeth.amino.model <- cbind(Estimate = coef(num.trimeth.amino), confint(num.trimeth.amino))
exp(num.trimeth.amino.model)
summary(num.trimeth.amino)

rpkm.trimeth.clinda <- glmmadmb(clinda.abund ~ tmp_smx.exp + (1|study_id), data = test9, family = "nbinom1")
rpkm.trimeth.clinda.model <- cbind(Estimate = coef(rpkm.trimeth.clinda), confint(rpkm.trimeth.clinda))
exp(rpkm.trimeth.clinda.model)
summary(rpkm.trimeth.clinda)

num.trimeth.clinda <- glmmadmb(num.clinda.ARG ~ tmp_smx.exp + (1|study_id), data = test9, family = "nbinom1")
num.trimeth.clinda.model <- cbind(Estimate = coef(num.trimeth.clinda), confint(num.trimeth.clinda))
exp(num.trimeth.clinda.model)
summary(num.trimeth.clinda)

rpkm.trimeth.macrolide <- glmmadmb(macrolide.abund ~ tmp_smx.exp + (1|study_id), data = test9, family = "nbinom1")
rpkm.trimeth.macrolide.model <- cbind(Estimate = coef(rpkm.trimeth.macrolide), confint(rpkm.trimeth.macrolide))
exp(rpkm.trimeth.macrolide.model)
summary(rpkm.trimeth.macrolide)

num.trimeth.macrolide <- glmmadmb(num.macrolide.ARG ~ tmp_smx.exp + (1|study_id), data = test9, family = "nbinom1")
num.trimeth.macrolide.model <- cbind(Estimate = coef(num.trimeth.macrolide), confint(num.trimeth.macrolide))
exp(num.trimeth.macrolide.model)
summary(num.trimeth.macrolide)

rpkm.trimeth.trimeth <- glmmadmb(trimeth.abund ~ tmp_smx.exp + (1|study_id), data = test9, family = "nbinom1")
rpkm.trimeth.trimeth.model <- cbind(Estimate = coef(rpkm.trimeth.trimeth), confint(rpkm.trimeth.trimeth))
exp(rpkm.trimeth.trimeth.model)
summary(rpkm.trimeth.trimeth)

num.trimeth.trimeth <- glmmadmb(num.trimeth.ARG ~ tmp_smx.exp + (1|study_id), data = test9, family = "nbinom1")
num.trimeth.trimeth.model <- cbind(Estimate = coef(num.trimeth.trimeth), confint(num.trimeth.trimeth))
exp(num.trimeth.trimeth.model)
summary(num.trimeth.trimeth)

rpkm.trimeth.tetra <- glmmadmb(tetra.abund ~ tmp_smx.exp + (1|study_id), data = test9, family = "nbinom1")
rpkm.trimeth.tetra.model <- cbind(Estimate = coef(rpkm.trimeth.tetra), confint(rpkm.trimeth.tetra))
exp(rpkm.trimeth.tetra.model)
summary(rpkm.trimeth.tetra)

num.trimeth.tetra <- glmmadmb(num.tetra.ARG ~ tmp_smx.exp + (1|study_id), data = test9, family = "nbinom1")
num.trimeth.tetra.model <- cbind(Estimate = coef(num.trimeth.tetra), confint(num.trimeth.tetra))
exp(num.trimeth.tetra.model)
summary(num.trimeth.tetra)

rpkm.trimeth.flagyl <- glmmadmb(flagyl.abund ~ tmp_smx.exp + (1|study_id), data = test9, family = "nbinom2")
rpkm.trimeth.flagyl.model <- cbind(Estimate = coef(rpkm.trimeth.flagyl), confint(rpkm.trimeth.flagyl))
exp(rpkm.trimeth.flagyl.model)
summary(rpkm.trimeth.flagyl)
#did not converge
#num.trimeth.flagyl <- glmmadmb(num.flagyl.ARG ~ tmp_smx.exp + (1|study_id), data = test9, family = "nbinom", zeroInflation = TRUE)
#num.trimeth.flagyl.model <- cbind(Estimate = coef(num.trimeth.flagyl), confint(num.trimeth.flagyl))
#exp(num.trimeth.flagyl.model)
#summary(num.trimeth.flagyl)

#adjusted P value for multiple comparisons
pvalues <- c(0.22, 0.0041, 0.22, 2.1e-11, 4e-06, 0.037, 9e-04, 0.21, 0.0078, 0.007, 0.31, 7.9e-05, 0.077, 0.0018, 0.059, 0.001, 0.0063, 0.024, 0.01, NA)
p.adjust(pvalues, method = "BH", n=20)
arg.class <- c("BetalactamAbund", "vancAbund", "carbAbund", "quinoloneAbund", "AminoglycosidAbund", "ClindamycinAbund", "MacrolideAbund", "TrimethoprimAbund", "TetracyclineAbund", "MetronidazoleAbund",
               "Betalactam.num", "vanc.num", "carb.num", "quinolone.num", "Aminoglycosid.num", "Clindamycin.num", "Macrolide.num", "Trimethoprim.num", "Tetracycline.num", "Metronidazole.num")
adj.p.trimeth<- c(
  2.444444e-01, 1.171429e-02, 2.444444e-01, 4.200000e-10,
  4.000000e-05, 5.692308e-02, 4.000000e-03, 2.444444e-01,
  1.560000e-02, 1.555556e-02, 3.263158e-01, 5.266667e-04,
  1.026667e-01, 6.000000e-03, 8.428571e-02, 4.000000e-03,
  1.555556e-02, 4.000000e-02, 1.818182e-02,           NA)
trimeth.bh <- data.frame(adj.p.trimeth, rownames = arg.class) %>% column_to_rownames(var = "rownames")
trimeth.bh <- trimeth.bh %>% format(adj.p.trimeth, scientific = FALSE) 
trimeth.bh$adj.p.trimeth <- as.numeric(trimeth.bh$adj.p.trimeth)
trimeth.bh$adj.p.trimeth <- signif(trimeth.bh$adj.p.trimeth, 2)

#macrolide exposure
rpkm.macrolide.betalactam <- glmmadmb(beta.abund ~ macrolide.exp + (1|study_id), data = test9, family = "nbinom1")
rpkm.macrolide.betalactam.model <- cbind(Estimate = coef(rpkm.macrolide.betalactam), confint(rpkm.macrolide.betalactam))
exp(rpkm.macrolide.betalactam.model)
summary(rpkm.macrolide.betalactam)

num.macrolide.betalactam <- glmmadmb(num.betalactam.ARG ~ macrolide.exp + (1|study_id), data = test9, family = "nbinom1")
num.macrolide.betalactam.model <- cbind(Estimate = coef(num.macrolide.betalactam), confint(num.macrolide.betalactam))
exp(num.macrolide.betalactam.model)
summary(num.macrolide.betalactam)

rpkm.macrolide.vanc <- glmmadmb(vanc.abund ~ macrolide.exp + (1|study_id), data = test9, family = "nbinom1")
rpkm.macrolide.vanc.model <- cbind(Estimate = coef(rpkm.macrolide.vanc), confint(rpkm.macrolide.vanc))
exp(rpkm.macrolide.vanc.model)
summary(rpkm.macrolide.vanc)

num.macrolide.vanc <- glmmadmb(num.vanc.ARG ~ macrolide.exp + (1|study_id), data = test9, family = "nbinom1")
num.macrolide.vanc.model <- cbind(Estimate = coef(num.macrolide.vanc), confint(num.macrolide.vanc))
exp(num.macrolide.vanc.model)
summary(num.macrolide.vanc)

rpkm.macrolide.carb <- glmmadmb(carb.abund ~ macrolide.exp + (1|study_id), data = test9, family = "nbinom1")
rpkm.macrolide.carb.model <- cbind(Estimate = coef(rpkm.macrolide.carb), confint(rpkm.macrolide.carb))
exp(rpkm.macrolide.carb.model)
summary(rpkm.macrolide.carb)

num.macrolide.carb <- glmmadmb(num.carbapenem.ARG ~ macrolide.exp + (1|study_id), data = test9, family = "nbinom1")
num.macrolide.carb.model <- cbind(Estimate = coef(num.macrolide.carb), confint(num.macrolide.carb))
exp(num.macrolide.carb.model)
summary(num.macrolide.carb)

rpkm.macrolide.macrolide <- glmmadmb(fluoro.abund ~ macrolide.exp + (1|study_id), data = test9, family = "nbinom2")
rpkm.macrolide.macrolide.model <- cbind(Estimate = coef(rpkm.macrolide.macrolide), confint(rpkm.macrolide.macrolide))
exp(rpkm.macrolide.macrolide.model)
summary(rpkm.macrolide.macrolide)

num.macrolide.macrolide <- glmmadmb(num.quinolone.ARG ~ macrolide.exp + (1|study_id), data = test9, family = "nbinom1")
num.macrolide.macrolide.model <- cbind(Estimate = coef(num.macrolide.macrolide), confint(num.macrolide.macrolide))
exp(num.macrolide.macrolide.model)
summary(num.macrolide.macrolide)

rpkm.macrolide.amino <- glmmadmb(amino.abund ~ macrolide.exp + (1|study_id), data = test9, family = "nbinom1")
rpkm.macrolide.amino.model <- cbind(Estimate = coef(rpkm.macrolide.amino), confint(rpkm.macrolide.amino))
exp(rpkm.macrolide.amino.model)
summary(rpkm.macrolide.amino)

num.macrolide.amino <- glmmadmb(num.amino.ARG ~ macrolide.exp + (1|study_id), data = test9, family = "nbinom1")
num.macrolide.amino.model <- cbind(Estimate = coef(num.macrolide.amino), confint(num.macrolide.amino))
exp(num.macrolide.amino.model)
summary(num.macrolide.amino)

rpkm.macrolide.clinda <- glmmadmb(clinda.abund ~ macrolide.exp + (1|study_id), data = test9, family = "nbinom1")
rpkm.macrolide.clinda.model <- cbind(Estimate = coef(rpkm.macrolide.clinda), confint(rpkm.macrolide.clinda))
exp(rpkm.macrolide.clinda.model)
summary(rpkm.macrolide.clinda)

num.macrolide.clinda <- glmmadmb(num.clinda.ARG ~ macrolide.exp + (1|study_id), data = test9, family = "nbinom1")
num.macrolide.clinda.model <- cbind(Estimate = coef(num.macrolide.clinda), confint(num.macrolide.clinda))
exp(num.macrolide.clinda.model)
summary(num.macrolide.clinda)

rpkm.macrolide.macrolide <- glmmadmb(macrolide.abund ~ macrolide.exp + (1|study_id), data = test9, family = "nbinom1")
rpkm.macrolide.macrolide.model <- cbind(Estimate = coef(rpkm.macrolide.macrolide), confint(rpkm.macrolide.macrolide))
exp(rpkm.macrolide.macrolide.model)
summary(rpkm.macrolide.macrolide)

num.macrolide.macrolide <- glmmadmb(num.macrolide.ARG ~ macrolide.exp + (1|study_id), data = test9, family = "nbinom1")
num.macrolide.macrolide.model <- cbind(Estimate = coef(num.macrolide.macrolide), confint(num.macrolide.macrolide))
exp(num.macrolide.macrolide.model)
summary(num.macrolide.macrolide)

rpkm.macrolide.trimeth <- glmmadmb(trimeth.abund ~ macrolide.exp + (1|study_id), data = test9, family = "nbinom1")
rpkm.macrolide.trimeth.model <- cbind(Estimate = coef(rpkm.macrolide.trimeth), confint(rpkm.macrolide.trimeth))
exp(rpkm.macrolide.trimeth.model)
summary(rpkm.macrolide.trimeth)

num.macrolide.trimeth <- glmmadmb(num.trimeth.ARG ~ macrolide.exp + (1|study_id), data = test9, family = "nbinom2")
num.macrolide.trimeth.model <- cbind(Estimate = coef(num.macrolide.trimeth), confint(num.macrolide.trimeth))
exp(num.macrolide.trimeth.model)
summary(num.macrolide.trimeth)

rpkm.macrolide.tetra <- glmmadmb(tetra.abund ~ macrolide.exp + (1|study_id), data = test9, family = "nbinom1")
rpkm.macrolide.tetra.model <- cbind(Estimate = coef(rpkm.macrolide.tetra), confint(rpkm.macrolide.tetra))
exp(rpkm.macrolide.tetra.model)
summary(rpkm.macrolide.tetra)

num.macrolide.tetra <- glmmadmb(num.tetra.ARG ~ macrolide.exp + (1|study_id), data = test9, family = "nbinom1")
num.macrolide.tetra.model <- cbind(Estimate = coef(num.macrolide.tetra), confint(num.macrolide.tetra))
exp(num.macrolide.tetra.model)
summary(num.macrolide.tetra)

rpkm.macrolide.flagyl <- glmmadmb(flagyl.abund ~ macrolide.exp + (1|study_id), data = test9, family = "nbinom2")
rpkm.macrolide.flagyl.model <- cbind(Estimate = coef(rpkm.macrolide.flagyl), confint(rpkm.macrolide.flagyl))
exp(rpkm.macrolide.flagyl.model)
summary(rpkm.macrolide.flagyl)

num.macrolide.flagyl <- glmmadmb(num.flagyl.ARG ~ macrolide.exp + (1|study_id), data = test9, family = "nbinom", zeroInflation = TRUE)
num.macrolide.flagyl.model <- cbind(Estimate = coef(num.macrolide.flagyl), confint(num.macrolide.flagyl))
exp(num.macrolide.flagyl.model)
summary(num.macrolide.flagyl)

#adjusted P value for multiple comparisons
pvalues <- c(0.67, 1.2e-06, 0.85, 0.11, 0.015, 0.049, 0.15, 0.22, 0.066, 0.82, 0.65, 0.00085, 0.82, 0.7, 0.21, 0.0036, 0.32, 0.34, 0.78, 0.99)
p.adjust(pvalues, method = "BH", n=20)
arg.class <- c("BetalactamAbund", "vancAbund", "carbAbund", "quinoloneAbund", "AminoglycosidAbund", "ClindamycinAbund", "MacrolideAbund", "macrolideoprimAbund", "TetracyclineAbund", "MetronidazoleAbund",
               "Betalactam.num", "vanc.num", "carb.num", "quinolone.num", "Aminoglycosid.num", "Clindamycin.num", "Macrolide.num", "macrolideoprim.num", "Tetracycline.num", "Metronidazole.num")
adj.p.macrolide<- c(
  0.8947368, 0.0000240, 0.8947368, 0.3142857, 0.0750000,
  0.1960000, 0.3750000, 0.4400000, 0.2200000, 0.8947368,
  0.8947368, 0.0085000, 0.8947368, 0.8947368, 0.4400000,
  0.0240000, 0.5666667, 0.5666667, 0.8947368, 0.9900000)
macrolide.bh <- data.frame(adj.p.macrolide, rownames = arg.class) %>% column_to_rownames(var = "rownames")
macrolide.bh <- macrolide.bh %>% format(adj.p.macrolide, scientific = FALSE) 
macrolide.bh$adj.p.macrolide <- as.numeric(macrolide.bh$adj.p.macrolide)
macrolide.bh$adj.p.macrolide <- signif(macrolide.bh$adj.p.macrolide, 2)

#other.abx exposure
rpkm.other.abx.betalactam <- glmmadmb(beta.abund ~ other.abx.exp + (1|study_id), data = test9, family = "nbinom1")
rpkm.other.abx.betalactam.model <- cbind(Estimate = coef(rpkm.other.abx.betalactam), confint(rpkm.other.abx.betalactam))
exp(rpkm.other.abx.betalactam.model)
summary(rpkm.other.abx.betalactam)

num.other.abx.betalactam <- glmmadmb(num.betalactam.ARG ~ other.abx.exp + (1|study_id), data = test9, family = "nbinom1")
num.other.abx.betalactam.model <- cbind(Estimate = coef(num.other.abx.betalactam), confint(num.other.abx.betalactam))
exp(num.other.abx.betalactam.model)
summary(num.other.abx.betalactam)

rpkm.other.abx.vanc <- glmmadmb(vanc.abund ~ other.abx.exp + (1|study_id), data = test9, family = "nbinom1")
rpkm.other.abx.vanc.model <- cbind(Estimate = coef(rpkm.other.abx.vanc), confint(rpkm.other.abx.vanc))
exp(rpkm.other.abx.vanc.model)
summary(rpkm.other.abx.vanc)

num.other.abx.vanc <- glmmadmb(num.vanc.ARG ~ other.abx.exp + (1|study_id), data = test9, family = "nbinom1")
num.other.abx.vanc.model <- cbind(Estimate = coef(num.other.abx.vanc), confint(num.other.abx.vanc))
exp(num.other.abx.vanc.model)
summary(num.other.abx.vanc)

rpkm.other.abx.carb <- glmmadmb(carb.abund ~ other.abx.exp + (1|study_id), data = test9, family = "nbinom1")
rpkm.other.abx.carb.model <- cbind(Estimate = coef(rpkm.other.abx.carb), confint(rpkm.other.abx.carb))
exp(rpkm.other.abx.carb.model)
summary(rpkm.other.abx.carb)

num.other.abx.carb <- glmmadmb(num.carbapenem.ARG ~ other.abx.exp + (1|study_id), data = test9, family = "nbinom1")
num.other.abx.carb.model <- cbind(Estimate = coef(num.other.abx.carb), confint(num.other.abx.carb))
exp(num.other.abx.carb.model)
summary(num.other.abx.carb)

rpkm.other.abx.other.abx <- glmmadmb(fluoro.abund ~ other.abx.exp + (1|study_id), data = test9, family = "nbinom2")
rpkm.other.abx.other.abx.model <- cbind(Estimate = coef(rpkm.other.abx.other.abx), confint(rpkm.other.abx.other.abx))
exp(rpkm.other.abx.other.abx.model)
summary(rpkm.other.abx.other.abx)

num.other.abx.other.abx <- glmmadmb(num.quinolone.ARG ~ other.abx.exp + (1|study_id), data = test9, family = "nbinom1")
num.other.abx.other.abx.model <- cbind(Estimate = coef(num.other.abx.other.abx), confint(num.other.abx.other.abx))
exp(num.other.abx.other.abx.model)
summary(num.other.abx.other.abx)

rpkm.other.abx.amino <- glmmadmb(amino.abund ~ other.abx.exp + (1|study_id), data = test9, family = "nbinom1")
rpkm.other.abx.amino.model <- cbind(Estimate = coef(rpkm.other.abx.amino), confint(rpkm.other.abx.amino))
exp(rpkm.other.abx.amino.model)
summary(rpkm.other.abx.amino)

num.other.abx.amino <- glmmadmb(num.amino.ARG ~ other.abx.exp + (1|study_id), data = test9, family = "nbinom1")
num.other.abx.amino.model <- cbind(Estimate = coef(num.other.abx.amino), confint(num.other.abx.amino))
exp(num.other.abx.amino.model)
summary(num.other.abx.amino)

rpkm.other.abx.clinda <- glmmadmb(clinda.abund ~ other.abx.exp + (1|study_id), data = test9, family = "nbinom1")
rpkm.other.abx.clinda.model <- cbind(Estimate = coef(rpkm.other.abx.clinda), confint(rpkm.other.abx.clinda))
exp(rpkm.other.abx.clinda.model)
summary(rpkm.other.abx.clinda)

num.other.abx.clinda <- glmmadmb(num.clinda.ARG ~ other.abx.exp + (1|study_id), data = test9, family = "nbinom1")
num.other.abx.clinda.model <- cbind(Estimate = coef(num.other.abx.clinda), confint(num.other.abx.clinda))
exp(num.other.abx.clinda.model)
summary(num.other.abx.clinda)

rpkm.other.abx.macrolide <- glmmadmb(macrolide.abund ~ other.abx.exp + (1|study_id), data = test9, family = "nbinom1")
rpkm.other.abx.macrolide.model <- cbind(Estimate = coef(rpkm.other.abx.macrolide), confint(rpkm.other.abx.macrolide))
exp(rpkm.other.abx.macrolide.model)
summary(rpkm.other.abx.macrolide)

num.other.abx.other.abx <- glmmadmb(num.macrolide.ARG ~ other.abx.exp + (1|study_id), data = test9, family = "nbinom1")
num.other.abx.other.abx.model <- cbind(Estimate = coef(num.other.abx.other.abx), confint(num.other.abx.other.abx))
exp(num.other.abx.other.abx.model)
summary(num.other.abx.other.abx)

rpkm.other.abx.trimeth <- glmmadmb(trimeth.abund ~ other.abx.exp + (1|study_id), data = test9, family = "nbinom1")
rpkm.other.abx.trimeth.model <- cbind(Estimate = coef(rpkm.other.abx.trimeth), confint(rpkm.other.abx.trimeth))
exp(rpkm.other.abx.trimeth.model)
summary(rpkm.other.abx.trimeth)

num.other.abx.trimeth <- glmmadmb(num.trimeth.ARG ~ other.abx.exp + (1|study_id), data = test9, family = "nbinom2")
num.other.abx.trimeth.model <- cbind(Estimate = coef(num.other.abx.trimeth), confint(num.other.abx.trimeth))
exp(num.other.abx.trimeth.model)
summary(num.other.abx.trimeth)

rpkm.other.abx.tetra <- glmmadmb(tetra.abund ~ other.abx.exp + (1|study_id), data = test9, family = "nbinom1")
rpkm.other.abx.tetra.model <- cbind(Estimate = coef(rpkm.other.abx.tetra), confint(rpkm.other.abx.tetra))
exp(rpkm.other.abx.tetra.model)
summary(rpkm.other.abx.tetra)

num.other.abx.tetra <- glmmadmb(num.tetra.ARG ~ other.abx.exp + (1|study_id), data = test9, family = "nbinom1")
num.other.abx.tetra.model <- cbind(Estimate = coef(num.other.abx.tetra), confint(num.other.abx.tetra))
exp(num.other.abx.tetra.model)
summary(num.other.abx.tetra)

rpkm.other.abx.flagyl <- glmmadmb(flagyl.abund ~ other.abx.exp + (1|study_id), data = test9, family = "nbinom2")
rpkm.other.abx.flagyl.model <- cbind(Estimate = coef(rpkm.other.abx.flagyl), confint(rpkm.other.abx.flagyl))
exp(rpkm.other.abx.flagyl.model)
summary(rpkm.other.abx.flagyl)

num.other.abx.flagyl <- glmmadmb(num.flagyl.ARG ~ other.abx.exp + (1|study_id), data = test9, family = "nbinom", zeroInflation = TRUE)
num.other.abx.flagyl.model <- cbind(Estimate = coef(num.other.abx.flagyl), confint(num.other.abx.flagyl))
exp(num.other.abx.flagyl.model)
summary(num.other.abx.flagyl)

#adjusted P value for multiple comparisons
pvalues <- c(6.5e-08, 0.57, 0.13, 0.11, 0.00011, 0.28, 0.28, 0.058, 0.68, 0.44, 0.2, 0.31, 0.061, 0.45, 0.042, 0.87, 0.63, 0.71, 0.36, 0.24)
p.adjust(pvalues, method = "BH", n=20)
arg.class <- c("BetalactamAbund", "vancAbund", "carbAbund", "quinoloneAbund", "AminoglycosidAbund", "ClindamycinAbund", "other.abxAbund", "other.abxoprimAbund", "TetracyclineAbund", "MetronidazoleAbund",
               "Betalactam.num", "vanc.num", "carb.num", "quinolone.num", "Aminoglycosid.num", "Clindamycin.num", "other.abx.num", "other.abxoprim.num", "Tetracycline.num")
adj.p.other.abx<- c(
  0.0000013, 0.7125000, 0.3714286, 0.3666667, 0.0011000,
  0.5090909, 0.5090909, 0.2440000, 0.7473684, 0.6000000,
  0.5000000, 0.5166667, 0.2440000, 0.6000000, 0.2440000,
  0.8700000, 0.7411765, 0.7473684, 0.5538462, 0.5090909)
other.abx.bh <- data.frame(adj.p.other.abx, rownames = arg.class) %>% column_to_rownames(var = "rownames")
other.abx.bh <- other.abx.bh %>% format(adj.p.other.abx, scientific = FALSE) 
other.abx.bh$adj.p.other.abx <- as.numeric(other.abx.bh$adj.p.other.abx)
other.abx.bh$adj.p.other.abx <- signif(other.abx.bh$adj.p.other.abx, 2)

#Figure 5
#Abundance heatmap
cef.exp<- c(0.91, 1.13, 0.58, 0.84, 1.04, 1.14, 1.36, 0.82, 0.93, 0.44)
pip.tazo.exp <- c(0.67, 11.80,  0.92,  1.43,  2.70,  1.38,  1.32,  1.19,  3.37,  0.62)
carb.exp <- c(1.57, 1.12, 1.18, 1.69, 1.20, 1.30, 1.50, 1.41, 1.02, 0.87)
vanc.ecp <- c(1.17, 1.08, 0.85, 1.47, 1.43, 1.07, 1.20, 0.79, 1.28, 0.59)
flagyl.exp <- c(0.90, 1.67, 0.86, 1.26, 1.46, 1.19, 1.34, 1.21, 1.16, 1.68)
quin.exp <- c(1.44, 0.87, 1.17, 1.25, 0.78, 1.53, 1.42, 1.36, 1.10, 1.12)
trimeth.exp <- c(0.88, 0.61, 0.86, 0.51, 0.64, 0.84, 0.77, 0.89, 0.81, 0.65)
mac.exp <- c(1.08, 3.50, 1.05, 1.31, 1.48, 1.36, 1.23, 1.25, 1.30, 1.07)
other.beta.exp <-c(0.48, 1.16, 1.40, 0.75, 0.50, 0.85, 0.85, 1.35, 0.94, 1.53)
other.exp <- c(0.58, 1.09, 1.43, 1.22, 0.82, 1.18, 1.17, 1.24, 1.06, 1.11)
abund.h <- rbind(cef.exp, pip.tazo.exp, carb.exp, other.beta.exp, vanc.ecp, flagyl.exp, quin.exp, trimeth.exp, mac.exp, other.exp)
colnames(abund.h) <- c("ABeta.Lactam", "CGlycopeptide", "BCarbapenem", "EFluoroquinolone", "IAminoglycoside", "JLincosamide", "GMacrolide", "FDiaminopyrimidine", "HTetracycline", "DNitroimidazole")
rownames(abund.h) <- c("Acefepime exposure", "Bpip.tazo.exp", "Dcarbapenem.exp", "Cother.beta.lactam.exp", "Evancomycin.exp", "Fmetronidazole.exp", "Gquinolone.exp", "Htrimethoprim.exp", "Imacrolide.exp",  "Jother.exp")
col.order <- c("ABeta.Lactam", "BCarbapenem", "CGlycopeptide", "DNitroimidazole", "EFluoroquinolone", "FDiaminopyrimidine", "GMacrolide", "HTetracycline", "IAminoglycoside", "JLincosamide")
row.order <- c("Acefepime exposure", "Bpip.tazo.exp", "Cother.beta.lactam.exp", "Dcarbapenem.exp", "EVancomycin.exp", "Fmetronidazole.exp", "Gquinolone.exp", "Htrimethoprim.exp", "Imacrolide.exp",  "Jother.exp")

abund.df <- data.frame(abund.h) %>% rownames_to_column(var = "Exposure")
abund.df <- melt(abund.df, id.id.vars = Exposure) %>% rename(Gene = variable, IRR = value)
abund.df <- abund.df %>%  mutate(Exposure = factor(Exposure))
abund.df$sig[abund.df$Exposure== "Acefepime exposure" & abund.df$Gene=="BCarbapenem" ] <- 1
abund.df$sig[abund.df$Exposure== "Acefepime exposure" & abund.df$Gene=="GMacrolide" ] <- 1
abund.df$sig[abund.df$Exposure== "Acefepime exposure" & abund.df$Gene== "FDiaminopyrimidine"] <- 1
abund.df$sig[abund.df$Exposure== "Acefepime exposure" & abund.df$Gene== "DNitroimidazole"] <- 1

abund.df$sig[abund.df$Exposure== "Bpip.tazo.exp" & abund.df$Gene=="ABeta.Lactam" ] <- 1
abund.df$sig[abund.df$Exposure== "Bpip.tazo.exp" & abund.df$Gene=="EFluoroquinolone" ] <- 1
abund.df$sig[abund.df$Exposure== "Bpip.tazo.exp" & abund.df$Gene=="JLincosamide" ] <- 1
abund.df$sig[abund.df$Exposure== "Bpip.tazo.exp" & abund.df$Gene== "CGlycopeptide"] <- 1
abund.df$sig[abund.df$Exposure== "Bpip.tazo.exp" & abund.df$Gene== "IAminoglycoside"] <- 1
abund.df$sig[abund.df$Exposure== "Bpip.tazo.exp" & abund.df$Gene== "HTetracycline"] <- 1
abund.df$sig[abund.df$Exposure== "Bpip.tazo.exp" & abund.df$Gene== "GMacrolide" ] <- 1

abund.df$sig[abund.df$Exposure== "Evancomycin.exp" & abund.df$Gene== "ABeta.Lactam" ] <- 1
abund.df$sig[abund.df$Exposure== "Evancomycin.exp" & abund.df$Gene== "FDiaminopyrimidine" ] <- 1
abund.df$sig[abund.df$Exposure== "Evancomycin.exp" & abund.df$Gene== "IAminoglycoside" ] <- 1
abund.df$sig[abund.df$Exposure== "Evancomycin.exp" & abund.df$Gene== "DNitroimidazole" ] <- 1
abund.df$sig[abund.df$Exposure== "Evancomycin.exp" & abund.df$Gene== "GMacrolide" ] <- 1
abund.df$sig[abund.df$Exposure== "Evancomycin.exp" & abund.df$Gene== "HTetracycline" ] <- 1

abund.df$sig[abund.df$Exposure== "Fmetronidazole.exp" & abund.df$Gene== "CGlycopeptide" ] <- 1
abund.df$sig[abund.df$Exposure== "Fmetronidazole.exp" & abund.df$Gene== "IAminoglycoside" ] <- 1
abund.df$sig[abund.df$Exposure== "Fmetronidazole.exp" & abund.df$Gene== "GMacrolide" ] <- 1

abund.df$sig[abund.df$Exposure== "Htrimethoprim.exp" & abund.df$Gene== "GMacrolide" ] <- 1
abund.df$sig[abund.df$Exposure== "Htrimethoprim.exp" & abund.df$Gene== "CGlycopeptide" ] <- 1
abund.df$sig[abund.df$Exposure== "Htrimethoprim.exp" & abund.df$Gene== "IAminoglycoside" ] <- 1
abund.df$sig[abund.df$Exposure== "Htrimethoprim.exp" & abund.df$Gene== "HTetracycline" ] <- 1
abund.df$sig[abund.df$Exposure== "Htrimethoprim.exp" & abund.df$Gene== "DNitroimidazole" ] <- 1
abund.df$sig[abund.df$Exposure== "Htrimethoprim.exp" & abund.df$Gene== "EFluoroquinolone" ] <- 1

abund.df$sig[abund.df$Exposure== "Imacrolide.exp" & abund.df$Gene== "CGlycopeptide" ] <- 1

abund.df$sig[abund.df$Exposure== "Jother.exp" & abund.df$Gene== "ABeta.Lactam" ] <- 1


library(scales)
median(abund.df$IRR)
#rescaling IRR so that colors are more different
rescale.plot.neg <- abund.df %>% filter(IRR < 1)
rescale.plot.neg$IRR[rescale.plot.neg$IRR < 0.5] <- 0.5
normalize <- function(x, na.rm = TRUE) {
  return((x- min(x)) /(max(x)-min(x)))
}
rescale.plot.neg$rescale.IRR = normalize(rescale.plot.neg$IRR)
rescale.plot.neg <- rescale.plot.neg %>% mutate(plot.IRR = rescale.IRR-1) 
summary(rescale.plot.neg$rescale.IRR)
filter(rescale.plot.neg, rescale.IRR == median(rescale.plot.neg$rescale.IRR))
rescale.plot.neg <- rescale.plot.neg %>% dplyr::select(-rescale.IRR) 
rescale.plot.pos <- abund.df %>% filter(IRR >= 1) 
rescale.plot.pos$IRR[rescale.plot.pos$IRR >1.5] <- 1.5
hist(rescale.plot.pos$IRR)
rescale.plot.pos$plot.IRR = normalize(rescale.plot.pos$IRR)
summary(rescale.plot.pos$IRR)
filter(rescale.plot.pos, plot.IRR == median(rescale.plot.pos$plot.IRR))

rescale.plot <- rbind(rescale.plot.neg, rescale.plot.pos)


heat.abund <- ggplot(rescale.plot, aes(x =Gene, y=Exposure, fill = plot.IRR))+ 
  geom_tile(color = "black", size = 0.5) +
  geom_point(aes(x =Gene, y = Exposure, size = sig),
             shape = 8, na.rm = TRUE, stroke = 1.5,
             show.legend = FALSE)+
  scale_fill_gradient2(low = "navy",  high = "orangered2", midpoint = 0,
                       breaks = c(-1, 0, 1),
                       labels = c(expression(""<=0.5), "1.0",expression("">="1.5")),
                       name = expression(beta),
                       na.value = "gray70")+
  labs(y = "Antibiotic Exposure", x = "Antibiotic Resistance Gene Class")+
  scale_x_discrete(
    limits = c("ABeta.Lactam", "BCarbapenem", "CGlycopeptide", "DNitroimidazole", "EFluoroquinolone", "FDiaminopyrimidine", "GMacrolide", "HTetracycline","IAminoglycoside", "JLincosamide"),
    labels = c("Beta Lactam", "Carbapenem", "Glycopeptide", "Nitroimidazole", "Fluoroquinolone", "Diaminopyrimidine", "Macrolide", "Tetracycline", "Aminoglycoside", "Lincosamide"),
    position = "top")+
  scale_y_discrete(
    limits = rev(c("Acefepime exposure", "Bpip.tazo.exp", "Dcarbapenem.exp", "Evancomycin.exp", "Fmetronidazole.exp", "Gquinolone.exp", "Htrimethoprim.exp", "Imacrolide.exp")),
    labels = rev(c("CEF", "TZP", "CBP", "VAN", "MTZ", "FLQ", "SXT", "MAC")),
    position = "left")+
  theme_minimal()+
  theme(axis.text.x.top = element_text(angle = 90, hjust = 0, vjust = 0.5, size = 26, color = "black"),
        axis.title.x.top = element_text(vjust = 2, size = 30),
        axis.title.y = element_text(vjust = 3, size = 26),
        axis.text.y = element_text(size = 26, color = "black"),
        axis.title = element_text(size = 30),
        axis.ticks = element_blank(),
        legend.title = element_text(size = 24),
        legend.text = element_text(size = 24),
        plot.margin = margin(t = 2.5,
                             b=0,
                             l=0, r=0, unit = "cm"))+
  coord_fixed()

heat.abund

#Number heatmap
cef.exp<- c(0.63, 1.02, 0.62, 0.65, 0.71, 0.83, 0.81, 0.83, 0.72, NA)
pip.tazo.exp <- c(0.56, 1.23, 0.76, 0.71, 0.74, 0.85, 0.78, 0.93, 0.69, 0.69)
carb.exp <- c(0.77, 1.43, 0.95, 0.94, 0.99, 0.97, 1.02, 1.31, 0.81, NA)
vanc.ecp <- c(0.70, 0.93, 0.74, 0.72, 0.72, 0.79, 0.80, 0.84, 0.72, 0.63)
flagyl.exp <- c(0.67, 1.31, 0.83, 0.86, 0.82, 0.79, 0.83, 0.97, 0.73, 0.77)
quin.exp <- c(1.29, 1.33, 1.21, 1.40, 1.24, 1.21, 1.24, 1.34, 1.26, 0.94)
trimeth.exp <- c(0.92, 0.56, 0.83, 0.76, 0.87, 0.83, 0.85, 0.85, 0.87, NA)
mac.exp <- c(1.08, 2.07, 0.95, 1.07, 1.19, 1.33, 1.12, 1.14, 1.03, 1.00)
other.beta.exp <-c(1.21, 1.28, 1.39, 1.12, 1.30, 0.98, 1.05, 1.05, 1.10, 1.30)
other.exp <- c(1.07, 1.16, 1.27, 1.03, 1.15, 0.94, 1.02, 1.11, 0.99, 1.09)
num.h <- rbind(cef.exp, pip.tazo.exp, carb.exp, other.beta.exp, vanc.ecp, flagyl.exp, quin.exp, trimeth.exp, mac.exp, other.exp)
colnames(num.h) <- c("ABeta.Lactam", "CGlycopeptide", "BCarbapenem", "EFluoroquinolone", "IAminoglycoside", "JLincosamide", "GMacrolide", "FDiaminopyrimidine", "HTetracycline", "DNitroimidazole")
rownames(num.h) <- c("Acefepime exposure", "Bpip.tazo.exp", "Dcarbapenem.exp", "Cother.beta.lactam.exp", "Evancomycin.exp", "Fmetronidazole.exp", "Gquinolone.exp", "Htrimethoprim.exp", "Imacrolide.exp",  "Jother.exp")
col.order <- c("ABeta.Lactam", "BCarbapenem", "CGlycopeptide", "DNitroimidazole", "EFluoroquinolone", "FDiaminopyrimidine", "GMacrolide", "HTetracycline", "IAminoglycoside", "JLincosamide")
row.order <- c("Acefepime exposure", "Bpip.tazo.exp", "Cother.beta.lactam.exp", "Dcarbapenem.exp", "EVancomycin.exp", "Fmetronidazole.exp", "Gquinolone.exp", "Htrimethoprim.exp", "Imacrolide.exp",  "Jother.exp")

num.df <- data.frame(num.h) %>% rownames_to_column(var = "Exposure")
num.df <- melt(num.df, id.id.vars = Exposure) %>% rename(Gene = variable, IRR = value)
num.df <- num.df %>%  mutate(Exposure = factor(Exposure))

num.df$sig[num.df$Exposure== "Acefepime exposure" & num.df$Gene== "ABeta.Lactam"] <- 1
num.df$sig[num.df$Exposure== "Acefepime exposure" & num.df$Gene== "BCarbapenem"] <- 1
num.df$sig[num.df$Exposure== "Acefepime exposure" & num.df$Gene== "EFluoroquinolone"] <- 1
num.df$sig[num.df$Exposure== "Acefepime exposure" & num.df$Gene== "IAminoglycoside"] <- 1
num.df$sig[num.df$Exposure== "Acefepime exposure" & num.df$Gene== "JLincosamide"] <- 1
num.df$sig[num.df$Exposure== "Acefepime exposure" & num.df$Gene== "GMacrolide"] <- 1
num.df$sig[num.df$Exposure== "Acefepime exposure" & num.df$Gene== "HTetracycline"] <- 1
num.df$sig[num.df$Exposure== "Acefepime exposure" & num.df$Gene== "FDiaminopyrimidine"] <- 1

num.df$sig[num.df$Exposure== "Bpip.tazo.exp" & num.df$Gene== "ABeta.Lactam"] <- 1
num.df$sig[num.df$Exposure== "Bpip.tazo.exp" & num.df$Gene== "EFluoroquinolone"] <- 1
num.df$sig[num.df$Exposure== "Bpip.tazo.exp" & num.df$Gene== "IAminoglycoside"] <- 1
num.df$sig[num.df$Exposure== "Bpip.tazo.exp" & num.df$Gene== "GMacrolide"] <- 1
num.df$sig[num.df$Exposure== "Bpip.tazo.exp" & num.df$Gene== "HTetracycline"] <- 1

num.df$sig[num.df$Exposure== "Evancomycin.exp" & num.df$Gene== "ABeta.Lactam"] <- 1
num.df$sig[num.df$Exposure== "Evancomycin.exp" & num.df$Gene== "BCarbapenem"] <- 1
num.df$sig[num.df$Exposure== "Evancomycin.exp" & num.df$Gene== "EFluoroquinolone"] <- 1
num.df$sig[num.df$Exposure== "Evancomycin.exp" & num.df$Gene== "IAminoglycoside"] <- 1
num.df$sig[num.df$Exposure== "Evancomycin.exp" & num.df$Gene== "JLincosamide"] <- 1
num.df$sig[num.df$Exposure== "Evancomycin.exp" & num.df$Gene== "GMacrolide"] <- 1
num.df$sig[num.df$Exposure== "Evancomycin.exp" & num.df$Gene== "FDiaminopyrimidine"] <- 1
num.df$sig[num.df$Exposure== "Evancomycin.exp" & num.df$Gene== "HTetracycline"] <- 1
num.df$sig[num.df$Exposure== "Evancomycin.exp" & num.df$Gene== "DNitroimidazole"] <- 1

num.df$sig[num.df$Exposure== "Fmetronidazole.exp" & num.df$Gene== "ABeta.Lactam"] <- 1
num.df$sig[num.df$Exposure== "Fmetronidazole.exp" & num.df$Gene== "JLincosamide"] <- 1
num.df$sig[num.df$Exposure== "Fmetronidazole.exp" & num.df$Gene== "GMacrolide"] <- 1
num.df$sig[num.df$Exposure== "Fmetronidazole.exp" & num.df$Gene== "HTetracycline"] <- 1

num.df$sig[num.df$Exposure== "Htrimethoprim.exp" & num.df$Gene== "CGlycopeptide"] <- 1
num.df$sig[num.df$Exposure== "Htrimethoprim.exp" & num.df$Gene== "JLincosamide"] <- 1
num.df$sig[num.df$Exposure== "Htrimethoprim.exp" & num.df$Gene== "HTetracycline"] <- 1
num.df$sig[num.df$Exposure== "Htrimethoprim.exp" & num.df$Gene== "EFluoroquinolone"] <- 1
num.df$sig[num.df$Exposure== "Htrimethoprim.exp" & num.df$Gene== "GMacrolide"] <- 1
num.df$sig[num.df$Exposure== "Htrimethoprim.exp" & num.df$Gene== "FDiaminopyrimidine"] <- 1

num.df$sig[num.df$Exposure== "Imacrolide.exp" & num.df$Gene== "CGlycopeptide" ] <- 1
num.df$sig[num.df$Exposure== "Imacrolide.exp" & num.df$Gene== "JLincosamide" ] <- 1


#rescaling IRR so that colors are more different
rescale.plot.negnum <- num.df %>% filter(IRR < 1)
hist(rescale.plot.negnum$IRR)
rescale.plot.negnum$IRR[rescale.plot.negnum$IRR <0.5] <- 0.5
rescale.plot.numna <- num.df %>% filter(is.na(IRR)) %>% mutate(plot.IRR = NA)
normalize <- function(x, na.rm = TRUE) {
  return((x- min(x)) /(max(x)-min(x)))
}
rescale.plot.negnum$rescale.IRR = normalize(rescale.plot.negnum$IRR)
rescale.plot.negnum <- rescale.plot.negnum %>% mutate(plot.IRR = rescale.IRR-1) 
summary(rescale.plot.negnum$rescale.IRR)
filter(rescale.plot.negnum, rescale.IRR == median(rescale.plot.negnum$rescale.IRR))
rescale.plot.negnum <- rescale.plot.negnum %>% dplyr::select(-rescale.IRR) 
rescale.plot.posnum <- num.df %>% filter(IRR >= 1)
summary(rescale.plot.posnum$IRR)
hist(rescale.plot.posnum$IRR)
rescale.plot.posnum$IRR[rescale.plot.posnum$IRR >1.5]<- 1.5
rescale.plot.posnum$plot.IRR = normalize(rescale.plot.posnum$IRR)
filter(rescale.plot.posnum, plot.IRR == median(rescale.plot.posnum$plot.IRR))
rescale.plotnum <- rbind(rescale.plot.negnum, rescale.plot.posnum, rescale.plot.numna)

heat.num <- ggplot(rescale.plotnum, aes(x =Gene, y=Exposure, fill = plot.IRR, color = ""))+ 
  geom_tile(color = "black", size = 0.5) +
  geom_point(aes(x =Gene, y=Exposure, size = sig),
             shape = 8, na.rm = TRUE, stroke = 1.5,
             show.legend = FALSE,
             color="black")+
  scale_fill_gradient2(low = "navy",  high = "orangered2", midpoint = 0,
                       breaks = c(-1, 0, 1),
                       labels = c(expression(""<=0.5), "1.0", expression("">=1.5)),
                       name = expression(beta),
                       na.value = "gray 70")+
  labs(y = "Antibiotic Exposure", x = "Antibiotic Resistance Gene Class")+
  scale_x_discrete(
    limits = c("ABeta.Lactam", "BCarbapenem", "CGlycopeptide", "DNitroimidazole", "EFluoroquinolone", "FDiaminopyrimidine", "GMacrolide", "HTetracycline","IAminoglycoside", "JLincosamide"),
    labels = c("Beta Lactam", "Carbapenem", "Glycopeptide", "Nitroimidazole", "Fluoroquinolone", "Diaminopyrimidine", "Macrolide", "Tetracycline", "Aminoglycoside", "Lincosamide"),
    position = "top")+
  scale_y_discrete(
    limits = rev(c("Acefepime exposure", "Bpip.tazo.exp", "Dcarbapenem.exp", "Evancomycin.exp", "Fmetronidazole.exp", "Gquinolone.exp", "Htrimethoprim.exp", "Imacrolide.exp")),
    labels = rev(c("CEF", "TZP", "CBP", "VAN", "MTZ", "FLQ", "SXT", "MAC")),
    position = "left")+ 
  theme_minimal()+
  theme(axis.text.x.top = element_text(angle = 90, hjust = 0, vjust = 0.5, size = 26, color = "black"),
        axis.title.x.top = element_text(vjust = 2, color = "black", size = 30),
        axis.title.y = element_text(vjust = 3, size = 30),
        axis.text.y = element_text(size = 26, color = "black"),
        axis.title = element_text(size = 30),
        axis.ticks = element_blank(),
        legend.position = "none",
        plot.margin = margin(t = 2.5,
                             b=0,
                             l=0, r=0, unit = "cm"))+
  coord_fixed()

heat.num

heat <- ggarrange(heat.num , NULL, heat.abund + rremove("ylab"),  
                  labels = c("A. Number of ARGs", "", "B. ARG Abundance"), 
                  vjust = 1,
                  hjust = -0.5,
                  font.label = list(size = 32),
                  nrow = 1, widths = c(1, -0.1, 1.1))
heat