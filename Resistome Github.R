#The Effects of Antibiotic Exposures on the Gut Resistome During Hematopoietic Cell Transplantation in Children
#Heston SM, Young RR, Jenkins K, Martin PL, Stokhuyzen A, Ward DV, Bhattarai SK, Bucci V, Arshad M, Chao NJ, Seed PC, Kelly MS.

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
phy.card <- readRDS("T:/Projects/PID/Kelly/Intestinal Microbiome of  - Pro00064365/Data Preprocessing/phy.card.09222023_pruned_500k_paired_reads.rds")
phy.card <- prune_taxa(taxa_sums(phy.card)>0, phy.card)
otu.table <- data.frame(otu_table(phy.card))
metadata_card <- data.frame(sample_data(phy.card)) #691 samples
num <- metadata_card %>% distinct(study_id, .keep_all = TRUE) #80
tax.table<- data.frame(tax_table(phy.card)) 
nrow(tax.table)## 372 unique ARG

#creating tax table with arg antibiotic class
arg.raw <- read.csv("T:/Projects/PID/Kelly/Intestinal Microbiome of  - Pro00064365/Resistome Analysis/Raw Data/aro_index_09222023.csv") #from CARD
colnames(arg.raw)[1] = "ARO.Accession"
arg_card <- data.frame(otu_table(phy.card))
arg_card <- tibble::rownames_to_column(arg_card, var = "concated_column")

arg.raw$ARO.Accession <- stringr::str_replace_all(arg.raw$ARO.Accession, ":", "_")
arg.table <- arg.raw %>% dplyr::select(Resistance.Mechanism, Drug.Class, ARO.Name, ARO.Accession)
arg.table <- arg.table[,c(4,3,2,1)]

merge <- arg.table %>% filter(ARO.Accession %in% arg_card$concated_column)

remove(arg_card, arg.raw)

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

#Assigning ARG class to each ARG
merge <- merge %>% mutate(clinical = case_when(str_detect(ARO.Name, "van")~ "van", 
                                                                      str_detect(ARO.Name, "KPC|IMI|SME|OXA|NDM|VIM|IMP|GIM|SPM")~ "carbapenemase",
                                                                      str_detect(ARO.Name, "SHV|TEM|CTX|PER|GES|VEB|BES|CME|SFO|BEL|TLA")~ "esbl"))
merge <- merge %>% mutate(amino = if_else(str_detect(Drug.Class, "aminoglycoside"), 1, 0)) %>% mutate(carbapenem = if_else(str_detect(Drug.Class, "carbapenem|penem"), 1, 0)) %>% mutate(betalactam = if_else(str_detect(Drug.Class, "cephalosporin|cephamycin|penam|monobactam"), 1, 0)) %>% 
  mutate(trimeth = if_else(str_detect(Drug.Class, "diaminopyrimidine"), 1, 0)) %>% mutate(fluoro = if_else(str_detect(Drug.Class, "fluoroquinolone"), 1, 0)) %>% mutate(macrolide = if_else(str_detect(Drug.Class, "macrolide"), 1, 0)) %>% mutate(phenicol = if_else(str_detect(Drug.Class, "phenicol"), 1, 0)) %>% 
  mutate(sulfa = if_else(str_detect(Drug.Class, "sulfonamide"), 1, 0)) %>% mutate(tetra = if_else(str_detect(Drug.Class, "tetracycline|glycylcycline"), 1, 0)) %>% mutate(clinda = if_else(str_detect(Drug.Class, "lincosamide"), 1, 0)) %>% mutate(vanc = if_else(str_detect(Drug.Class, "glycopeptide"), 1, 0)) %>% 
  mutate(linezolid = if_else(str_detect(Drug.Class, "oxazolidinone"), 1, 0), flagyl = if_else(str_detect(Drug.Class, "nitroimidazole "), 1, 0))
merge <- merge %>% dplyr::select(Resistance.Mechanism, Drug.Class, ARO.Name, ARO.Accession, clinical, amino, carbapenem, betalactam,trimeth, fluoro, macrolide, phenicol, sulfa, tetra, clinda, vanc, linezolid, flagyl) %>% column_to_rownames(var = "ARO.Accession")
merge <- merge %>% rownames_to_column(var="OTU")

#adding updated ARG classifications to phyloseq
otu.card <- as.matrix(otu.table, rownames.force = NA)
nrow(otu.card)
OTU=otu_table(otu.card, taxa_are_rows = TRUE)
nrow(OTU)
tax.card <- merge %>% arrange(OTU) %>% column_to_rownames(var = "OTU")
tax.card <- as.matrix(tax.card, rownames.force = NA)
TAX = tax_table(tax.card)
nrow(TAX)
meta <- metadata_card %>% mutate(SampleID=paste0("X", SampleID)) 
rownames(meta) <- NULL
meta <- meta %>% column_to_rownames(var="SampleID")
META = sample_data(meta)
phy.card <- phyloseq(OTU, TAX, META)

#Adding in abundances for each ARG class to metadata
colnames(merge)
merge <- merge %>% column_to_rownames(var="OTU")
arg_otu <- otu.table %>% rownames_to_column(var="ARG")
amino.abun <- merge %>% dplyr::select(5) %>% rownames_to_column(var="ARG")
amino.otu <- left_join(arg_otu, amino.abun, by="ARG")
amino.otu <- amino.otu %>% filter(amino==1) 
nrow(amino.otu) 
amino.otu <- column_to_rownames(amino.otu, var="ARG")
amino.genes <- data.frame(colSums(amino.otu))
amino.genes <- amino.genes %>% dplyr::rename(c("amino.abund" = "colSums.amino.otu.")) %>% rownames_to_column(var="SampleID") %>% separate(SampleID, c(NA, "SampleID"), sep = 1, remove = TRUE)
metad <- data.frame(sample_data(phy.card))
metad <- rownames_to_column(metad, var="SampleID")
metad$SampleID <- str_remove_all(metad$SampleID, "X")
metad <- left_join(metad, amino.genes, by = "SampleID")

carb.abun <- merge %>% dplyr::select(6) %>% rownames_to_column(var="ARG")
carb.otu <- left_join(arg_otu, carb.abun, by="ARG")
carb.otu <- carb.otu %>% filter(carbapenem==1)
carb.otu <- column_to_rownames(carb.otu, var="ARG")
carb.genes <- data.frame(colSums(carb.otu))
carb.genes <- carb.genes %>% dplyr::rename(c("carb.abund" = "colSums.carb.otu.")) %>% rownames_to_column(var="SampleID") %>% separate(SampleID, c(NA, "SampleID"), sep = 1, remove = TRUE)
metad <- left_join(metad, carb.genes, by = "SampleID")

beta.abun <- merge %>% dplyr::select(7)%>% rownames_to_column(var="ARG")
beta.otu <- left_join(arg_otu, beta.abun, by="ARG")
beta.otu <- beta.otu %>% filter(betalactam==1)
beta.otu <- column_to_rownames(beta.otu, var="ARG")
beta.genes <- data.frame(colSums(beta.otu))
beta.genes <- beta.genes %>% dplyr::rename(c("beta.abund" = "colSums.beta.otu.")) %>% rownames_to_column(var="SampleID") %>% separate(SampleID, c(NA, "SampleID"), sep = 1, remove = TRUE)
metad <- left_join(metad, beta.genes, by = "SampleID")

trimeth.abun <- merge %>% dplyr::select(8)%>% rownames_to_column(var="ARG")
trimeth.otu <- left_join(arg_otu, trimeth.abun, by="ARG")
trimeth.otu <- trimeth.otu %>% filter(trimeth==1)
trimeth.otu <- column_to_rownames(trimeth.otu, var="ARG")
trimeth.genes <- data.frame(colSums(trimeth.otu))
trimeth.genes <- trimeth.genes %>% dplyr::rename(c("trimeth.abund" = "colSums.trimeth.otu.")) %>% rownames_to_column(var="SampleID") %>% separate(SampleID, c(NA, "SampleID"), sep = 1, remove = TRUE)
metad <- left_join(metad, trimeth.genes, by = "SampleID")

fluoro.abun <- merge %>% dplyr::select(9)%>% rownames_to_column(var="ARG")
fluoro.otu <- left_join(arg_otu, fluoro.abun, by="ARG")
fluoro.otu <- fluoro.otu %>% filter(fluoro==1)
fluoro.otu <- column_to_rownames(fluoro.otu, var="ARG")
fluoro.genes <- data.frame(colSums(fluoro.otu))
fluoro.genes <- fluoro.genes %>% dplyr::rename(c("fluoro.abund" = "colSums.fluoro.otu.")) %>% rownames_to_column(var="SampleID") %>% separate(SampleID, c(NA, "SampleID"), sep = 1, remove = TRUE)
metad <- left_join(metad, fluoro.genes, by = "SampleID")

macrolide.abun <- merge %>% dplyr::select(10)%>% rownames_to_column(var="ARG")
macrolide.otu <- left_join(arg_otu, macrolide.abun, by="ARG")
macrolide.otu <- macrolide.otu %>% filter(macrolide==1)
macrolide.otu <- column_to_rownames(macrolide.otu, var="ARG")
macrolide.genes <- data.frame(colSums(macrolide.otu))
macrolide.genes <- macrolide.genes %>% dplyr::rename(c("macrolide.abund" = "colSums.macrolide.otu.")) %>% rownames_to_column(var="SampleID") %>% separate(SampleID, c(NA, "SampleID"), sep = 1, remove = TRUE)
metad <- left_join(metad, macrolide.genes, by = "SampleID")

phenicol.abun <- merge %>% dplyr::select(11)%>% rownames_to_column(var="ARG")
phenicol.otu <- left_join(arg_otu, phenicol.abun, by="ARG")
phenicol.otu <- phenicol.otu %>% filter(phenicol==1)
phenicol.otu <- column_to_rownames(phenicol.otu, var="ARG")
phenicol.genes <- data.frame(colSums(phenicol.otu))
phenicol.genes <- phenicol.genes %>% dplyr::rename(c("phenicol.abund" = "colSums.phenicol.otu.")) %>% rownames_to_column(var="SampleID") %>% separate(SampleID, c(NA, "SampleID"), sep = 1, remove = TRUE)
metad <- left_join(metad, phenicol.genes, by = "SampleID")

sulfa.abun <- merge %>% dplyr::select(12)%>% rownames_to_column(var="ARG")
sulfa.otu <- left_join(arg_otu, sulfa.abun, by="ARG")
sulfa.otu <- sulfa.otu %>% filter(sulfa==1)
sulfa.otu <- column_to_rownames(sulfa.otu, var="ARG")
sulfa.genes <- data.frame(colSums(sulfa.otu))
sulfa.genes <- sulfa.genes %>% dplyr::rename(c("sulfa.abund" = "colSums.sulfa.otu.")) %>% rownames_to_column(var="SampleID") %>% separate(SampleID, c(NA, "SampleID"), sep = 1, remove = TRUE)
metad <- left_join(metad, sulfa.genes, by = "SampleID")

tetra.abun <- merge %>% dplyr::select(13)%>% rownames_to_column(var="ARG")
tetra.otu <- left_join(arg_otu, tetra.abun, by="ARG")
tetra.otu <- tetra.otu %>% filter(tetra==1)
tetra.otu <- column_to_rownames(tetra.otu, var="ARG")
tetra.genes <- data.frame(colSums(tetra.otu))
tetra.genes <- tetra.genes %>% dplyr::rename(c("tetra.abund" = "colSums.tetra.otu.")) %>% rownames_to_column(var="SampleID") %>% separate(SampleID, c(NA, "SampleID"), sep = 1, remove = TRUE)
metad <- left_join(metad, tetra.genes, by = "SampleID")

clinda.abun <- merge %>% dplyr::select(14)%>% rownames_to_column(var="ARG")
clinda.otu <- left_join(arg_otu, clinda.abun, by="ARG")
clinda.otu <- clinda.otu %>% filter(clinda==1)
clinda.otu <- column_to_rownames(clinda.otu, var="ARG")
clinda.genes <- data.frame(colSums(clinda.otu))
clinda.genes <- clinda.genes %>% dplyr::rename(c("clinda.abund" = "colSums.clinda.otu.")) %>% rownames_to_column(var="SampleID") %>% separate(SampleID, c(NA, "SampleID"), sep = 1, remove = TRUE)
metad <- left_join(metad, clinda.genes, by = "SampleID")

vanc.abun <- merge %>% dplyr::select(15)%>% rownames_to_column(var="ARG")
vanc.otu <- left_join(arg_otu, vanc.abun, by="ARG")
vanc.otu <- vanc.otu %>% filter(vanc==1)
vanc.otu <- column_to_rownames(vanc.otu, var="ARG")
vanc.genes <- data.frame(colSums(vanc.otu))
vanc.genes <- vanc.genes %>% dplyr::rename(c("vanc.abund" = "colSums.vanc.otu.")) %>% rownames_to_column(var="SampleID") %>% separate(SampleID, c(NA, "SampleID"), sep = 1, remove = TRUE)
table(vanc.abun$vanc)
metad <- left_join(metad, vanc.genes, by = "SampleID")

linezolid.abun <- merge %>% dplyr::select(16)%>% rownames_to_column(var="ARG")
linezolid.otu <- left_join(arg_otu, linezolid.abun, by="ARG")
linezolid.otu <- linezolid.otu %>% filter(linezolid==1)
linezolid.otu <- column_to_rownames(linezolid.otu, var="ARG")
linezolid.genes <- data.frame(colSums(linezolid.otu))
linezolid.genes <- linezolid.genes %>% dplyr::rename(c("linezolid.abund" = "colSums.linezolid.otu.")) %>% rownames_to_column(var="SampleID") %>% separate(SampleID, c(NA, "SampleID"), sep = 1, remove = TRUE)
metad <- left_join(metad, linezolid.genes, by = "SampleID")


flagyl.abun <- merge %>% dplyr::select(17)%>% rownames_to_column(var="ARG")
flagyl.otu <- left_join(arg_otu, flagyl.abun, by="ARG")
flagyl.otu <- flagyl.otu %>% filter(flagyl==1)
flagyl.otu <- column_to_rownames(flagyl.otu, var="ARG")
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

#determining most abundant ARG class
most.abund <- metad %>% dplyr::select(48:60)
most.abund <- data.frame(colSums(most.abund))
most.abund <- most.abund %>% dplyr::rename("total.abund" = "colSums.most.abund.") %>% dplyr::arrange(desc(total.abund)) %>% mutate(r.abund = round(total.abund))
head(most.abund, 13) 
remove(most.abund)

#Calculate number of ARGs per sample
newdata <- data.frame(otu_table(phy.card))
for(i in 1: ncol(newdata)){newdata[,i] <- ifelse(newdata[,i]  >0, 1, 0)}

#number of ARGs per sample
tot.num <- data.frame(colSums(newdata))
tot.num <- tot.num %>% dplyr::rename(num.ARG = colSums.newdata.) %>% rownames_to_column(var = "sample") %>% separate("sample", into = c("X", "SampleID"), sep = 1, extra = "merge") %>% dplyr::select(SampleID, num.ARG)
summary(tot.num$num.ARG) 

#Number of ARGs per sample
common <- data.frame(rowSums(newdata))
common <- common %>% dplyr::arrange(desc(rowSums.newdata.)) 
hist(common$rowSums.newdata.)
common <- rownames_to_column(common, var="ARO")
gene.name <- rownames_to_column(merge, var="ARO") %>%
  dplyr::select(ARO, ARO.Name)
common <- left_join(common, gene.name, by = "ARO")
common$percent = (common$rowSums.newdata./691)*100
remove(common)

# Number of ARGs detected in each class.
class.num <- merge %>% dplyr::select(.,5:17)
class.num <- data.frame(colSums(class.num)) %>% arrange(desc(colSums.class.num.))
remove(class.num)

#most abundant ARGs
arg.abund <- data.frame(otu_table(phy.card))
arg.abund <- arg.abund %>% mutate(total.abund = rowSums(arg.abund))
arg.abund <- rownames_to_column(arg.abund, var="ARG.Name")
arg.abund <- dplyr::select(arg.abund, c(ARG.Name, total.abund))
arg.abund <- left_join(arg.abund, arg.table, by=c("ARG.Name" = "ARO.Accession"))
top20.ARG.abund <- arg.abund %>% arrange(desc(total.abund)) %>% slice_max(order_by=total.abund, n=20)

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
for(i in 4: 375){new2[,i] <- ifelse((new2[,i] == 1) & (lag(new2[,i]) == 0) & (new2$study.id == lag(new2$study.id)), 1, 0)}
any(is.na(new2))
which(is.na(new2),arr.ind = TRUE)
new2[,190]
new1[1:2,190]
new2[,365]
new1[1:2,365]
new2$ARO_3002867[new2$SampleID == "01ST001"] <- 0
new2$ARO_3007047[new2$SampleID == "01ST001"] <- 0
new2$new.ARG <- rowSums(new2[,4:375])
na <- filter(new2, is.na(new.ARG))
new3 <- new2 %>% ungroup() %>% dplyr::select(SampleID, new.ARG)

#instability since last sample = jaccard distance
for(i in 4: 375){new4[,i] <- ifelse((new4[,i] == 0) & (lag(new4[,i]) == 1) & (new4$study.id == lag(new4$study.id)), 1, 0)}
new4[1, 4:375] <- 0
new4$lost.ARG <- rowSums(new4[,4:375])
new5 <- new4 %>% ungroup() %>% dplyr::select(SampleID, lost.ARG)
for(i in 4: 375){new6[,i] <- ifelse((new6[,i] == 1) & (lag(new6[,i]) == 1) & (new6$study.id == lag(new6$study.id)), 1, 0)}
new6[1, 4:375] <- 0
new6$stable.ARG <- rowSums(new6[4:375])
new7 <- new6 %>% ungroup() %>% dplyr::select(SampleID, stable.ARG, day, study.id)
instability <- left_join(new3, new5, by = "SampleID")
instability <- left_join(instability, new7, by = "SampleID")
instability <- left_join(instability, tot.num, by = "SampleID")
instability <- instability %>% arrange(study.id, day) %>% group_by(study.id)
hist(instability$day)
instability <- instability %>% mutate(inst.score = (1-((stable.ARG)/(stable.ARG + new.ARG + lost.ARG))), week = if_else(day <= -3, -1,
                                                                                                                if_else(day >-3 & day <=3, 0,
                                                                                                                if_else(day >3 & day <= 10, 1,
                                                                                                                if_else(day >10 & day <=17, 2,
                                                                                                                if_else(day >17 & day <=24, 3,
                                                                                                                if_else(day >24 & day <= 33, 4,
                                                                                                                if_else(day >33 & day <= 40, 5,
                                                                                                                if_else(day >40 & day <= 47, 6, 
                                                                                                                if_else(day >47 & day <= 54, 7, 
                                                                                                                if_else(day >54 & day <= 61, 8,
                                                                                                                if_else(day >61 & day <= 68, 9,
                                                                                                                if_else(day >68 & day <= 75, 10,
                                                                                                                if_else(day >75 & day <= 82, 11,
                                                                                                                if_else(day >82 & day <= 89, 12,13)))))))))))))))

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

test <- daily.seq %>% dplyr::select(-c(day_study, location, engrafted, st, st2, bl, bl2, or, cu, tpn, bsi_onset, bcx_org, cdiff.x, stool.id, pct_metagenomic, primary_dx, donor_hsct, day_engraftment, hla, hsct_prep, gvhd_prophy, prior_probx, mortality_100d, mortality_1y, mortality_2y, agvhd_gutliver, stage_gut, stage_liver, bsi1_day, bsi1, bsi2_day, bsi2, cdiff.y, day_cdiff))
test$seq[is.na(test$seq)] <- 0
test2 <- test %>% mutate(other.predict = tetracycline + other_beta_lactam + clinda + ag, beta.cef.predict = cefepime + other_beta_lactam) %>% mutate(other = if_else(other.predict >0, 1, 0), beta.cef = if_else(beta.cef.predict >0, 1, 0))
test2 <- test2 %>%  dplyr::arrange(study_id, day) %>% group_by(study_id) %>% mutate(seq.int = cumsum(seq)) %>% group_by(study_id, seq.int) 

#Determining antibiotic exposure for each antibiotic
test2.a <- test2  %>%  dplyr::arrange(study_id, day)  %>% group_by(study_id, seq.int) %>% mutate(cefepime.cum = cumsum(cefepime)) %>% ungroup() %>% group_by(study_id) %>% mutate(cefepime.lag= lag(cefepime.cum)) %>% mutate(cefepime.exp= if_else(seq == 1 & cefepime.lag > 0, 1, 0))
check.cefepime <- test2.a %>% filter(day == min(day) & cefepime ==1 & seq ==1) 
test2.a$cefepime.exp[is.na(test2.a$cefepime.exp)] <- 0

test2.b <- test2.a %>%  dplyr::arrange(study_id, day)  %>% group_by(study_id, seq.int) %>% mutate(piptazo.cum = cumsum(piptazo)) %>% ungroup() %>% group_by(study_id) %>% mutate(piptazo.lag= lag(piptazo.cum)) %>% mutate(piptazo.exp= if_else(seq == 1 & piptazo.lag > 0, 1, 0))
check.piptazo <- test2.b %>% filter(day == min(day) & piptazo ==1 & seq ==1) 
test2.b$piptazo.exp[is.na(test2.b$piptazo.exp)] <- 0

test2.c <- test2.b %>%  dplyr::arrange(study_id, day)  %>% group_by(study_id, seq.int) %>% mutate(carbapenem.cum = cumsum(carbapenem)) %>% ungroup() %>% group_by(study_id) %>% mutate(carbapenem.lag= lag(carbapenem.cum)) %>% mutate(carbapenem.exp= if_else(seq == 1 & carbapenem.lag > 0, 1, 0))
check.carbapenem <- test2.c %>% filter(day == min(day) & carbapenem ==1 & seq ==1) 
test2.c$carbapenem.exp[is.na(test2.c$carbapenem.exp)] <- 0

test2.d <- test2.c  %>%  dplyr::arrange(study_id, day)  %>% group_by(study_id, seq.int) %>% mutate(flagyl.cum = cumsum(flagyl)) %>% ungroup() %>% group_by(study_id) %>% mutate(flagyl.lag= lag(flagyl.cum)) %>% mutate(flagyl.exp= if_else(seq == 1 & flagyl.lag > 0, 1, 0))
check.flagyl <- test2.d %>% filter(day == min(day) & flagyl ==1 & seq ==1) 
test2.d$flagyl.exp[is.na(test2.d$flagyl.exp)] <- 0

test2.d1 <- test2.d  %>%  dplyr::arrange(study_id, day)  %>% group_by(study_id, seq.int) %>% mutate(vanc.cum = cumsum(vanc)) %>% ungroup() %>% group_by(study_id) %>% mutate(vanc.lag= lag(vanc.cum)) %>% mutate(vanc.exp= if_else(seq == 1 & vanc.lag > 0, 1, 0))
check.vanc <- test2.d1 %>% filter(day == min(day) & vanc ==1 & seq ==1) 
test2.d1$vanc.exp[is.na(test2.d1$vanc.exp)] <- 0

test2.e <- test2.d1 %>%  dplyr::arrange(study_id, day)   %>% group_by(study_id, seq.int) %>% mutate(quinolone.cum = cumsum(quinolone)) %>% ungroup() %>% group_by(study_id) %>% mutate(quinolone.lag= lag(quinolone.cum)) %>% mutate(quinolone.exp= if_else(seq == 1 & quinolone.lag > 0, 1, 0))
check.quinolone <- test2.e %>% filter(day == min(day) & quinolone ==1 & seq ==1) 
test2.e$quinolone.exp[is.na(test2.e$quinolone.exp)] <- 0

test2.f <- test2.e  %>%  dplyr::arrange(study_id, day) %>% group_by(study_id, seq.int) %>% mutate(tetracycline.cum = cumsum(tetracycline)) %>% ungroup() %>% group_by(study_id) %>% mutate(tetracycline.lag= lag(tetracycline.cum)) %>% mutate(tetracycline.exp= if_else(seq == 1 & tetracycline.lag > 0, 1, 0))
check.tetracycline <- test2.f %>% filter(day == min(day) & tetracycline ==1 & seq ==1) 
test2.f$tetracycline.exp[is.na(test2.f$tetracycline.exp)] <- 0

test2.g <- test2.f %>%  dplyr::arrange(study_id, day)  %>% group_by(study_id, seq.int) %>% mutate(other_beta_lactam.cum = cumsum(other_beta_lactam)) %>% ungroup() %>% group_by(study_id) %>% mutate(other_beta_lactam.lag= lag(other_beta_lactam.cum)) %>% mutate(other_beta_lactam.exp= if_else(seq == 1 & other_beta_lactam.lag > 0, 1, 0))
check.other_beta_lactam <- test2.g %>% filter(day == min(day) & other_beta_lactam ==1 & seq ==1) 
test2.g$other_beta_lactam.exp[is.na(test2.g$other_beta_lactam.exp)] <- 0
test2.g$other_beta_lactam.exp[test2.g$study_id == 36 & test2.g$day == -27 & test2.g$other_beta_lactam == 1] <- 1
test2.g$other_beta_lactam.exp[test2.g$study_id == 47 & test2.g$day == -21 & test2.g$other_beta_lactam == 1] <- 1
test2.g$other_beta_lactam.exp[test2.g$study_id == 50 & test2.g$day == -14 & test2.g$other_beta_lactam == 1] <- 1

test2.h <- test2.g  %>%  dplyr::arrange(study_id, day) %>% group_by(study_id, seq.int) %>% mutate(clinda.cum = cumsum(clinda)) %>% ungroup() %>% group_by(study_id) %>% mutate(clinda.lag= lag(clinda.cum)) %>% mutate(clinda.exp= if_else(seq == 1 & clinda.lag > 0, 1, 0))
check.clinda <- test2.h %>% filter(day == min(day) & clinda ==1 & seq ==1) 
test2.h$clinda.exp[is.na(test2.h$clinda.exp)] <- 0

test2.i <- test2.h  %>%  dplyr::arrange(study_id, day) %>% group_by(study_id, seq.int) %>% mutate(ag.cum = cumsum(ag)) %>% ungroup() %>% group_by(study_id) %>% mutate(ag.lag= lag(ag.cum)) %>% mutate(ag.exp= if_else(seq == 1 & ag.lag > 0, 1, 0))
check.ag <- test2.i %>% filter(day == min(day) & ag ==1 & seq ==1) 
test2.i$ag.exp[is.na(test2.i$ag.exp)] <- 0

test2.j <- test2.i  %>%  dplyr::arrange(study_id, day) %>% group_by(study_id, seq.int) %>% mutate(macrolide.cum = cumsum(macrolide)) %>% ungroup() %>% group_by(study_id) %>% mutate(macrolide.lag= lag(macrolide.cum)) %>% mutate(macrolide.exp= if_else(seq == 1 & macrolide.lag > 0, 1, 0))
check.macrolide <- test2.j %>% filter(day == min(day) & macrolide ==1 & seq ==1) 
test2.j$macrolide.exp[is.na(test2.j$macrolide.exp)] <- 0

test2.k <- test2.j  %>%  dplyr::arrange(study_id, day) %>% group_by(study_id, seq.int) %>% mutate(tmp_smx.cum = cumsum(tmp_smx)) %>% ungroup() %>% group_by(study_id) %>% mutate(tmp_smx.lag= lag(tmp_smx.cum)) %>% mutate(tmp_smx.exp= if_else(seq == 1 & tmp_smx.lag > 0, 1, 0))
check.tmp_smx <- test2.k %>% filter(day == min(day) & tmp_smx ==1 & seq ==1) 
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
check.beta.cef <- test2.m %>% filter(day == min(day) & beta.cef ==1 & seq ==1) 
test2.m$beta.cef.exp[is.na(test2.m$beta.cef.exp)] <- 0
test2.m$beta.cef.exp[test2.m$study_id == 36 & test2.m$day == -27 & test2.m$beta.cef == 1] <- 1
test2.m$beta.cef.exp[test2.m$study_id == 47 & test2.m$day == -21 & test2.m$beta.cef == 1] <- 1
test2.m$beta.cef.exp[test2.m$study_id == 50 & test2.m$day == -14 & test2.m$beta.cef == 1] <- 1

test3 <- test2.m %>% dplyr::select(-c(contains(".cum"), contains(".lag"), contains("predict"), cefepime, piptazo, carbapenem, vanc, flagyl,
                                      clinda, ag, macrolide, tmp_smx, other_beta_lactam, beta.cef, other, quinolone, tetracycline)) %>% ungroup()
remove(check.cefepime, check.piptazo, check.carbapenem, check.vanc, check.flagyl, check.clinda, check.ag,
       check.macrolide, check.tmp_smx, check.other_beta_lactam, check.beta.cef, check.other, check.quinolone, check.tetracycline,
       test2, test2.a, test2.b, test2.c, test2.d, test2.d1, test2.e, test2.f, test2.g, test2.h, test2.i, test2.j, test2.k, test2.l)

test4 <- left_join(test3, instability, by = "SampleID")

#long by patient day with sequencing data
test4 <- test4 %>% filter(seq ==1)

#antibiotic exposures
#percent of samples exposed to each antibiotic
abx.exp <- test4 %>% summarize(cefepime = sum(cefepime.exp), vanc = sum(vanc.exp), carbapenem = sum(carbapenem.exp), quinolone = sum(quinolone.exp), ag = sum(ag.exp), clinda = sum(clinda.exp), macrolide= sum(macrolide.exp), tmp_smx = sum(tmp_smx.exp), tetracycline = sum(tetracycline.exp), flagyl = sum(flagyl.exp), piptazo = sum(piptazo.exp), other_beta_lactam = sum(other_beta_lactam.exp))
abx.exp1 <- data.frame(t(abx.exp)) %>% dplyr::rename(freq = t.abx.exp.) %>% mutate(pct.sample = round(freq/691 *100, digits = 0)) %>% dplyr::arrange(desc(freq))

#percent of children exposed to each antibiotic
abx.exp2 <- test4 %>% dplyr::arrange(study_id) %>% group_by(study_id) %>% summarise(cefepime = sum(cefepime.exp), vanc = sum(vanc.exp), carb = sum(carbapenem.exp), quin = sum(quinolone.exp), amino = sum(ag.exp), clinda = sum(clinda.exp), macro= sum(macrolide.exp), tmp = sum(tmp_smx.exp), tetra = sum(tetracycline.exp), flagyl = sum(flagyl.exp), piptazo = sum(piptazo.exp), other_beta = sum(other_beta_lactam.exp)) %>%
  mutate(cef.pt = if_else(cefepime >0, 1, 0), vanc.pt = if_else(vanc >0, 1, 0), carb.pt = if_else(carb >0, 1, 0), quin.pt = if_else(quin >0, 1, 0), amino.pt = if_else(amino >0, 1, 0), clinda.pt = if_else(clinda >0, 1, 0), macro.pt = if_else(macro >0, 1, 0), tmp.pt = if_else(tmp >0, 1, 0), tetra.pt = if_else(tetra >0, 1, 0), flagyl.pt = if_else(flagyl >0, 1, 0), piptazo.pt = if_else(piptazo >0, 1, 0), other_beta.pt = if_else(other_beta >0, 1, 0))
groups(abx.exp2)
abx.pt3 <- abx.exp2 %>% summarize(cefepime = sum(cef.pt), vanc = sum(vanc.pt), carb = sum(carb.pt), quin = sum(quin.pt), amino = sum(amino.pt), clinda = sum(clinda.pt), macro= sum(macro.pt), tmp = sum(tmp.pt), tetra = sum(tetra.pt), flagyl = sum(flagyl.pt), piptazo = sum(piptazo.pt), other_beta = sum(other_beta.pt))
abx.exp4 <- data.frame(t(abx.pt3)) %>% dplyr::rename(freq = t.abx.pt3.) %>% dplyr::arrange(desc(freq)) %>% mutate(pct.patient = round(freq/80, digits = 2) *100)

#anaerobic exposures
abx.exp5 <- test4 %>% dplyr::arrange(study_id) %>% group_by(study_id) %>% mutate(anaerobe=if_else(sum(piptazo.exp) >0 | sum(carbapenem.exp) >0 | sum(flagyl.exp) >0 | sum(clinda.exp) >0, 1, 0))
abx.exp5 <- abx.exp5[c(1:4, 45:64)]
abx.exp5 <- abx.exp5 %>% distinct(study_id, .keep_all = TRUE)
table(abx.exp5$anaerobe) #58 received anaerobic
remove(abx.exp, abx.exp1, abx.exp2, abx.pt3, abx.exp4)

#preparing variables for multivariable analyses
test6 <- test4 %>% mutate(prep = if_else(hsct_prep_cat == "myelo", 1, 2), log.depth = log(final.pair1))
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
phy.bmt <- readRDS("T:/Projects/PID/Kelly/Intestinal Microbiome of  - Pro00064365/Data Preprocessing/phy.bmt.09222023_pruned_25K_paired_reads.rds")
p.bmt <- prune_samples(phy.bmt@sam_data$SampleID %in% metad$SampleID, phy.bmt)
nsamples(p.bmt) 
p.bmt <- prune_taxa(taxa_sums(p.bmt)>0, p.bmt)
ntaxa(p.bmt)
otu.bmt <- data.frame(otu_table(p.bmt)) 
tax.bmt <- data.frame(tax_table(p.bmt))
meta.bmt <- data.frame(sample_data(p.bmt))
relative_bmt <- psmelt(p.bmt)
species <- aggregate(relative_bmt$Abundance, by=list(Species=relative_bmt$Species), FUN=mean)
species <- arrange(species, desc(x))
topspecies <- unique(species$Species[1:200])
top.200 <- data.frame(topspecies) %>% filter(str_detect(topspecies, "[:digit:]"))
meta <- meta.bmt %>% mutate(SampleID=paste0("X", SampleID)) 

#number of species per sample
num.species <- data.frame(otu_table(p.bmt)) 
for(i in 1: ncol(num.species)){num.species[,i] <- ifelse(num.species[,i] >0, 1, 0)}

spec.num <- data.frame(colSums(num.species))
spec.num <- spec.num %>% rownames_to_column(var = "SampleID") %>% dplyr::rename("num.species" = "colSums.num.species.") 
summary(spec.num$num.species) 
meta_bmt <- left_join(meta, spec.num, by = "SampleID")

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
for(i in 4: 1098){new.sp2[,i] <- ifelse((new.sp2[,i] == 1) & (lag(new.sp2[,i]) == 0) & (new.sp2$study.id == lag(new.sp2$study.id)), 1, 0)}
which(is.na(new.sp2),arr.ind = TRUE)
new.sp2$Ruminococcus_gnavus[new.sp2$SampleID == "01ST001"] <- 0
new.sp2$Terrisporobacter_othiniensis[new.sp2$SampleID == "01ST001"] <- 0
new.sp2$Clostridium_innocuum[new.sp2$SampleID == "01ST001"] <- 0
new.sp2$Erysipelatoclostridium_ramosum[new.sp2$SampleID == "01ST001"] <- 0
new.sp2$new.sp <- rowSums(new.sp2[,4:1098])
na <- filter(new.sp2, is.na(new.sp))
new.sp3 <- new.sp2 %>% ungroup() %>% dplyr::select(SampleID, new.sp)
summary(new.sp3$new.sp)
#Instability of bacterial species
for(i in 4:1098){new.sp4[,i] <- ifelse((new.sp4[,i] == 0) & (lag(new.sp4[,i]) == 1) & (new.sp4$study.id == lag(new.sp4$study.id)), 1, 0)}
new.sp4[1, 4:1098] <- 0
new.sp4$lost.sp <- rowSums(new.sp4[,4:1098])
new.sp5 <- new.sp4 %>% ungroup() %>% dplyr::select(SampleID, lost.sp)

for(i in 4:1098){new.sp6[,i] <- ifelse((new.sp6[,i] == 1) & (lag(new.sp6[,i]) == 1) & (new.sp6$study.id == lag(new.sp6$study.id)), 1, 0)}
new.sp6[1, 4:1098] <- 0
new.sp6$stable.sp <- rowSums(new.sp6[4:1098])
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

new.sp8 <- instability.sp %>% ungroup() %>% 
  mutate(SampleID = paste0("X", SampleID)) %>% 
  dplyr::select(SampleID, new.sp, inst.score)
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
abx.exps <- abx.exps %>% mutate(SampleID = paste0("X", SampleID))
meta_bmt <- left_join(meta_bmt, abx.exps, by = "SampleID") 
meta_bmt <- meta_bmt %>% column_to_rownames(var = "SampleID")
relative_bmt <- psmelt(p.bmt)
species <- aggregate(relative_bmt$Abundance, by=list(Phylum=relative_bmt$Phylum, Genus = relative_bmt$Genus, Species=relative_bmt$Species), FUN=mean)
species <- arrange(species, desc(x))
topspecies <- unique(species$Species[1:75])

#Most abundant specieS
relative.bmt <- relative_bmt 
genera <- aggregate(relative.bmt$Abundance, by=list(Genus = relative.bmt$Genus), FUN = mean)
topgenera <- genera %>% filter(!str_detect(Genus, "unclassified")) %>% arrange(desc(x)) %>% slice_max(n=19, order_by = x)
sum(topgenera$x)
relative.bmt <- relative.bmt %>% mutate(plot.Genera = if_else(Genus %in% topgenera$Genus, Genus, "Other"))
table(relative.bmt$plot.Genera)

top.species <- aggregate(relative.bmt$Abundance, by=list(Species = relative.bmt$Species), FUN=mean)
top.species2 <- top.species %>% filter(!str_detect(Species, "unclassified")) %>% arrange(desc(x)) %>% slice_max(n=19, order_by=x) %>%
  mutate(pct = x*100)
sum(top.species2$x)
relative.bmt <- relative.bmt %>% mutate(plot.Species= if_else(Species %in% top.species2$Species, Species, "Other"))
table(relative.bmt$plot.Species)

relative.bmt$plot.Genera <- factor(relative.bmt$plot.Genera, levels = c("Akkermansia", "Alistipes", "Bacteroides", "Bifidobacterium", "Blautia", "Clostridioides", "Eggerthella", "Enterobacter", "Enterococcus", "Erysipelatoclostridium", "Escherichia", "Faecalibacterium", "Flavonifractor", "Hungatella", "Intestinibacter", "Klebsiella", "Parabacteroides", "Ruthenibacterium", "Sellimonas", "Other"))
relative.bmt$plot.Species <- str_replace_all(relative.bmt$plot.Species, "_", " ")

relative.bmt$plot.Species <- factor(relative.bmt$plot.Species, levels = c("Escherichia coli","Enterococcus faecalis", "Phocaeicola vulgatus", "Enterococcus faecium", "Bacteroides stercoris",
 "Bacteroides fragilis", "Ruminococcus gnavus", "Bifidobacterium breve", "Parabacteroides distasonis", "Alistipes onderdonkii", "Enterocloster bolteae", "Bacteroides ovatus", "Bifidobacterium longum", 
 "Faecalibacterium prausnitzii", "Bacteroides uniformis", "Bacteroides thetaiotaomicron", "Klebsiella pneumoniae", "Blautia wexlerae", "Ruthenibacterium lactatiformans", "Other"))

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
tax.card <- data.frame(tax_table(phy.card))
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
class.num$SampleID <- paste0("X", class.num$SampleID)
test9 <- left_join(test8, class.num, by = "SampleID")
test9 <- test9 %>% mutate(other.abx.exp = if_else(tetracycline.exp == 1, 1, if_else(clinda.exp == 1, 1, if_else(ag.exp == 1, 1, if_else(other_beta_lactam.exp == 1, 1, 0)))))

#One antibiotic exposure variable
one.abx.test8 <- test8 %>% mutate(antibiotic.exposure =
                                    if_else(aerobic.exp==1 & anaerobic.exp==0, 1,
                                            if_else(aerobic.exp==0 & anaerobic.exp==1, 2,
                                                    if_else(aerobic.exp==0& anaerobic.exp==0, 0, 3))))


#####ANALYSES AND FIGURES#####
#MICROBIOME
#Figure 1
colors_20  <-  c("gold1", "#E31A1C", "#6A3D9A", "dodgerblue2", "gray50",
                 "palegreen2", "green4", "#FF7F00", "skyblue2", "#CAB2D6", 
                 "#FDBF6F", # lt orange
                 "maroon","deeppink1","blue1","steelblue4", 
                 "darkturquoise", "darkorange4","brown", "#FB9A99", "black")
#relative.bmt$week <- as.numeric(relative.bmt$week) 
class(relative.bmt$plot.Species)
levels(relative.bmt$plot.Species)
relative.bmt$plot.Species<- factor(relative.bmt$plot.Species, levels=rev(c("Other",
  "Escherichia coli",                "Enterococcus faecalis"  ,        
  "Phocaeicola vulgatus",            "Enterococcus faecium"   ,        
  "Bacteroides stercoris",           "Bacteroides fragilis"   ,        
  "Ruminococcus gnavus",             "Bifidobacterium breve"  ,        
  "Parabacteroides distasonis",      "Alistipes onderdonkii"  ,        
  "Enterocloster bolteae",           "Bacteroides ovatus"             ,
  "Bifidobacterium longum",          "Faecalibacterium prausnitzii"   ,
  "Bacteroides uniformis",           "Bacteroides thetaiotaomicron"   ,
  "Klebsiella pneumoniae",           "Blautia wexlerae"               ,
  "Ruthenibacterium lactatiformans"
)))
plot_weeks <- ggplot(arrange(relative.bmt, Abundance), aes(x=week, y=Abundance, fill=plot.Species)) +
  geom_bar(stat="identity", position="fill") +
  guides(fill = guide_legend(ncol=1, title="Species")) +#byrow=TRUE 
  scale_fill_manual(values=colors_20, labels=rev(c(
    "Other",
    expression(italic("Escherichia coli",))               , expression(italic("Enterococcus faecalis"))  ,        
    expression(italic("Phocaeicola vulgatus",))           , expression(italic("Enterococcus faecium"))   ,        
    expression(italic("Bacteroides stercoris",))          , expression(italic("Bacteroides fragilis"))   ,        
    expression(italic("Ruminococcus gnavus",))            , expression(italic("Bifidobacterium breve"))  ,        
    expression(italic("Parabacteroides distasonis",))     , expression(italic("Alistipes onderdonkii"))  ,        
    expression(italic("Enterocloster bolteae",))          , expression(italic("Bacteroides ovatus"))             ,
    expression(italic("Bifidobacterium longum",))         , expression(italic("Faecalibacterium prausnitzii"))   ,
    expression(italic("Bacteroides uniformis",))          , expression(italic("Bacteroides thetaiotaomicron"))   ,
    expression(italic("Klebsiella pneumoniae",))          , expression(italic("Blautia wexlerae"))               ,
                      expression(italic("Ruthenibacterium lactatiformans")))))+ 
  scale_y_continuous(limits = c(0, 1), expand = expansion(add =c(0,0.02)))+
  scale_x_discrete(expand = c(0,0))+ #limits = c(-2.5, 14.5), #breaks = c(-30, -20, -10, 0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100), 
  theme_classic()+
  xlab("Week Relative to HCT") + ylab("Relative Abundance")+
  theme(legend.text=element_text(size=24),
        legend.text.align = 0,
        legend.title=element_text(size=24), 
        axis.title.y = element_text(size=24), 
        axis.text = element_text(size=24, colour="black"), 
        axis.title.x = element_text(size=24), 
        plot.title = element_text(size=24, hjust=1.4), 
        strip.text.x = element_text(size = 24), strip.background=element_rect(fill="white")) 
plot_weeks

#Bray-Curtis distances
bray.plot <- subset_samples(p.bmt, week== -2 |week==0|week==3|week==7 )

ord_theme   <-  theme(panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA),
                      axis.text.y = element_text(size=18, family="Arial Black", color="grey20"),
                      axis.title.y = element_text(size=20, family="Arial Black", color="grey20", margin = margin(t = 0, r = 10, b = 0, l = 0)),
                      axis.text.x  = element_text(size=18, angle=90, hjust=1, vjust=0.4, family="Arial Black", color="grey20"),
                      axis.title.x = element_text(size=20, family="Arial Black", color="grey20", margin = margin(t = 10, r = 0, b = 0, l = 0)),
                      strip.text.x=element_text(size=18,angle=0, family="Arial Black", color="grey20"),
                      strip.background=element_rect(fill="black"), legend.position = "right",
                      legend.text = element_text(size=20, family="Arial Black", color="grey20"),
                      legend.title = element_text(size=20))

ord <- ordinate(bray.plot, method = "PCoA", distance="bray")

bray.time <- plot_ordination(bray.plot, ord, color = "week")+
  geom_point(size = 5) + ord_theme +
  stat_ellipse(aes(color=week), geom="polygon", alpha=0, type="t", level=0.8, size=1) +
  xlab("PCo 1 (9.9%)") + ylab("PCo 2 (7%)")+
  scale_color_discrete(breaks=c(-2, 0, 3, 7), name="Week Relative to HCT",
                     labels=c(" -2 Pre-HCT", "  0 HCT", "+3 Peri-engraftment", "+7 Post-engraftment"))
  
ggarrange(bray.time, plot_weeks,
          nrow = 2,
          labels = c("A","B"), 
          widths = c(2, 1),
          vjust = 1,
          font.label = list(size = 24, color = "black"))
 
#MaAsLin
library(Maaslin2)
meta_bmt_mas <- meta_bmt %>% mutate(log.depth = log(final.pair1))
meta_bmt_mas <- column_to_rownames(meta_bmt_mas, var="SampleID")
otu_bmt_mas <- data.frame(otu_table(p.bmt))

levels(meta_bmt_mas$diagnosis)
colnames(meta_bmt_mas)
class(meta_bmt_mas$age)

fit_data = Maaslin2(
  input_data = otu_bmt_mas,
  input_metadata = meta_bmt_mas,
  output = "All.mas_Feb24",
  fixed_effects = c("cefepime.exp", "piptazo.exp", "carbapenem.exp", "vanc.exp", "quinolone.exp", "tmp_smx.exp", "macrolide.exp", "flagyl.exp", "other.exp", "age", "sex", "diagnosis", "prep", "type_hsct", "day", "log.depth"),
  reference=c("diagnosis, Heme_malignancy"),
    random_effects = c("study_id"),
  normalization = "NONE",
  min_abundance = 0.01,
  min_prevalence = 0.05,
  max_significance = 0.10)

#Figure 2
significant.species <- read.csv("T:/Projects/PID/Kelly/Intestinal Microbiome of  - Pro00064365/Resistome Analysis/Analytic Data/Microbiome species and antibiotics maaslin significant_results 02192024.csv") 
significant.species <- significant.species %>% mutate(abs.effect = abs(coef)) %>% arrange(desc(abs.effect)) 
table(significant.species$metadata)
table(significant.species$value)
table(significant.species$feature)
just.sig <- significant.species %>% mutate(Effect = coef, Bug = feature, Covariate = if_else(metadata == "age", "Age", if_else(metadata=="cefepime.exp", "Cefepime", if_else(metadata == "day", "Day Relative to HCT", if_else(metadata == "diagnosis" & value == "Congenital_immunodef", "Immunodeficiency",
                                                                                                                                                                                                                               if_else(metadata =="diagnosis" & value == "Heme_non-malignancy", "Non-malignant Heme", if_else(metadata=="diagnosis" & value=="Metabolic_disease", "Metabolic Disease", if_else(metadata=="diagnosis" & value == "Solid_tumor", "Solid Tumor", 
                                                                                                                                                                                                                                                                                                                                                                                                               if_else(metadata == "flagyl.exp", "Metronidazole", if_else(metadata == "log.depth", "Sequencing Depth", if_else(metadata == "macrolide.exp", "Macrolide",
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               if_else(metadata == "other.exp", "Other Antibiotics", if_else(metadata== "piptazo.exp", "Pip-Tazo", if_else(metadata=="quinolone.exp", "Fluoroquinolone", if_else(metadata == "tmp_smx.exp", "TMP-SMX", if_else(metadata=="type_hsct", "Allogeneic", if_else(metadata=="vanc.exp", "Vancomycin", "check")))))))))))))))))
just.sig$Bug <- str_replace_all(just.sig$Bug, "_", " ")
supp.table <- just.sig %>% mutate(Species = Bug, Standard_Error = stderr, P_Value = pval, Q_Value = qval) %>%
  dplyr::select(Species, Covariate, Effect, Standard_Error, P_Value, Q_Value)
table(supp.table$Covariate)
supp.table$Covariate[supp.table$Covariate=="Pip-Tazo"] <- "Piperacillin-Tazobactam"
supp.table$Covariate[supp.table$Covariate=="Non-malignant Heme"] <- "Hematologic Non-Malignant Disease"
supp.table$Covariate[supp.table$Covariate=="TMP-SMX"] <- "Trimethoprim-Sulfamethoxazole"
supp.table <- supp.table %>% arrange(Covariate) %>% arrange(Species)
supp.table$Effect <- round(supp.table$Effect, 4)
supp.table$Standard_Error <- round(supp.table$Standard_Error, 4)
supp.table$P_Value <- round(supp.table$P_Value, 4)
supp.table$Q_Value <- round(supp.table$Q_Value, 4)

#all results in supp table
all.species <- read.csv("T:/Projects/PID/Kelly/Intestinal Microbiome of  - Pro00064365/Resistome Analysis/Analytic Data/Microbiome species and antibiotics maaslin all_results 02192024.csv") 
all.species <- all.species %>% mutate(abs.effect = abs(coef)) %>% arrange(desc(abs.effect)) 
table(all.species$metadata)
table(all.species$value)
table(all.species$feature)
just.all <- all.species %>% mutate(Effect = coef, Bug = feature, Covariate = if_else(metadata == "age", "Age", if_else(metadata=="cefepime.exp", "Cefepime", if_else(metadata == "day", "Day Relative to HCT", if_else(metadata == "diagnosis" & value == "Congenital_immunodef", "Immunodeficiency",
                                                                                                                                                                                                                               if_else(metadata =="diagnosis" & value == "Heme_non-malignancy", "Non-malignant Heme", if_else(metadata=="diagnosis" & value=="Metabolic_disease", "Metabolic Disease", if_else(metadata=="diagnosis" & value == "Solid_tumor", "Solid Tumor", 
                                                                                                                                                                                                                                                                                                                                                                                                               if_else(metadata == "flagyl.exp", "Metronidazole", if_else(metadata == "log.depth", "Sequencing Depth", if_else(metadata == "macrolide.exp", "Macrolide",
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               if_else(metadata == "other.exp", "Other Antibiotics", if_else(metadata== "piptazo.exp", "Pip-Tazo", if_else(metadata=="quinolone.exp", "Fluoroquinolone", if_else(metadata == "tmp_smx.exp", "TMP-SMX", if_else(metadata=="type_hsct", "Allogeneic", if_else(metadata=="vanc.exp", "Vancomycin", "check")))))))))))))))))
just.all$Bug <- str_replace_all(just.all$Bug, "_", " ")
supp.table <- just.all %>% mutate(Species = Bug, Standard_Error = stderr, P_Value = pval, Q_Value = qval) %>%
  dplyr::select(Species, Covariate, Effect, Standard_Error, P_Value, Q_Value)
table(supp.table$Covariate)
supp.table$Covariate[supp.table$Covariate=="Pip-Tazo"] <- "Piperacillin-Tazobactam"
supp.table$Covariate[supp.table$Covariate=="Non-malignant Heme"] <- "Hematologic Non-Malignant Disease"
supp.table$Covariate[supp.table$Covariate=="TMP-SMX"] <- "Trimethoprim-Sulfamethoxazole"
supp.table <- supp.table %>% arrange(Covariate) %>% arrange(Species)
just.sig <- just.sig %>% dplyr::select(Bug, Covariate, Effect)

#rescaling IRR so that colors are more different
rescale.plot.neg <- just.sig %>% filter(Effect < 0)
normalize <- function(x, na.rm = TRUE) {
  return((x- min(x)) /(max(x)-min(x)))
}
rescale.plot.neg$rescale.Effect = normalize(rescale.plot.neg$Effect)
rescale.plot.neg <- rescale.plot.neg %>% mutate(plot.Effect = rescale.Effect-1) 
summary(rescale.plot.neg$rescale.Effect)
filter(rescale.plot.neg, rescale.Effect == median(rescale.plot.neg$rescale.Effect))
rescale.plot.neg <- rescale.plot.neg %>% dplyr::select(-rescale.Effect) 
rescale.plot.pos <- just.sig %>% filter(Effect >= 0) 
hist(rescale.plot.pos$Effect)
rescale.plot.pos$plot.Effect = normalize(rescale.plot.pos$Effect)
summary(rescale.plot.pos$Effect)
filter(rescale.plot.pos, plot.Effect == median(rescale.plot.pos$plot.Effect))

rescale.plot <- rbind(rescale.plot.neg, rescale.plot.pos)
rescale.plot <- rescale.plot %>% mutate(sig = 1, common = paste0(Bug, Covariate))
significant <- rescale.plot %>% dplyr::select(common, sig)

all.result <- read.csv("T:/Projects/PID/Kelly/Intestinal Microbiome of  - Pro00064365/Resistome Analysis/Analytic Data/Microbiome species and antibiotics maaslin all_results 02192024.csv") #Supplementary Table 1
all.result <- all.result %>% mutate(Effect = coef)
table(all.result$feature)
all.result$Bug <- str_replace_all(all.result$feature, "_", " ")
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

heat.species.all <- ggplot(all.rescale.plot, aes(x =Covariate, y= Bug, fill = plot.Effect, color = ""))+
  geom_tile(color = "black", size = 0.5) +
  geom_point(aes(size = sig),
             shape = 8, na.rm = TRUE, stroke= 0.5,
             show.legend = FALSE,
             color="black")+
  scale_size_manual(values=c(2))+
  scale_fill_gradient2(low = "navy",  high = "orangered2", midpoint = 0,
                       breaks = c(-1, 0, 1),
                       labels = c("\u2264 1.5", "1.0", "\u2265 1.5"),
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

#RESISTOME
#Antibiotics effects on resistome (Supplementary Table 2)
#Abundance
one.abx.test8$antibiotic.exposure <- factor(one.abx.test8$antibiotic.exposure)
levels(one.abx.test8$antibiotic.exposure)
table(one.abx.test8$antibiotic.exposure) 

summary(one.abx.test8$total.abund)
gamma <- one.abx.test8 %>% mutate(nonzero.abund = total.abund + 1)
class(gamma$nonzero.abund)
summary(gamma$nonzero.abund)

#Number of Species
numspec.anaerobe.adj <- glmmadmb(num.species ~ antibiotic.exposure + age + sex + diagnosis + prep + type_hsct + day + log.depth + (1 |study_id), data  = one.abx.test8, family = "nbinom1")
summary(numspec.anaerobe.adj)
numspec.anaerobe.adj.model <- cbind(Estimate = coef(numspec.anaerobe.adj), confint(numspec.anaerobe.adj))
numspec.results <- data.frame(round(exp(numspec.anaerobe.adj.model),2))

#new species since previous sample
newspec.anaerobe.adj <- glmmadmb(new.sp ~ antibiotic.exposure + age + sex + diagnosis + prep + type_hsct + day + log.depth + (1 |study_id), data  = one.abx.test8, family = "nbinom1")
summary(newspec.anaerobe.adj)
newspec.anaerobe.adj.model <- cbind(Estimate = coef(newspec.anaerobe.adj), confint(newspec.anaerobe.adj))
newspec.results <- data.frame(round(exp(newspec.anaerobe.adj.model),2))

#Species instability
instspec.anaerobe.adj <- lmer(inst.score.sp ~ antibiotic.exposure + age + sex + diagnosis + prep + type_hsct + day + log.depth + (1 |study_id), data  = one.abx.test8)
summary(instspec.anaerobe.adj)
confint(instspec.anaerobe.adj)

#RESISTOME 
#Abundance
gamma.anaerobe <- glmer(nonzero.abund ~ antibiotic.exposure + age + sex + diagnosis + prep+ type_hsct + day + log.depth + (1 |study_id), data  = gamma, family=Gamma(link="log"))#prep  + diagnosislog.depth
summary(gamma.anaerobe)
gamma.an <- summary(gamma.anaerobe)
gamma.an <- as.data.frame(gamma.an$coefficients)
gamma.com <- confint(gamma.anaerobe, method="Wald")
gamma.com <- gamma.com[-c(1, 2),]
gamma.an1 <- cbind(gamma.an, gamma.com)
gamma.an1 <- gamma.an1[-1,]
gamma.an1 <- gamma.an1 %>% mutate(effect = round(exp(Estimate), digits =2), lower=round(exp(`2.5 %`), digits =2), upper=round(exp(`97.5 %`), digits =2)) %>%
  mutate(day.est = round(exp(Estimate), digits =3), day.lower = round(exp(`2.5 %`), digits =3), day.upper = round(exp(`97.5 %`), digits =4)) %>%
  dplyr::select(effect, lower, upper, day.est, day.lower, day.upper, `Pr(>|z|)`)%>%
  rownames_to_column(var = "row") 
gamma.an2 <- gamma.an1 %>% filter(row == "antibiotic.exposure1" | row == "antibiotic.exposure2"|row == "antibiotic.exposure3") %>% 
  mutate(drug = if_else(row == "antibiotic.exposure1", "Aerobic Antibiotics", if_else(row=="antibiotic.exposure2", "Anaerobic Antibiotics", if_else(row=="antibiotic.exposure3", "Aerobic and Anaerobic Antibiotics", "check")))) %>% 
  dplyr::select(drug, effect, lower, upper)


#Number of ARGs
numARG.anaerobe.adj <- glmmadmb(num.ARG ~ antibiotic.exposure + age + sex + diagnosis + prep + type_hsct + day + log.depth + (1 |study_id), data  = one.abx.test8, family = "nbinom1")
summary(numARG.anaerobe.adj)
numARG.anaerobe.adj.model <- cbind(Estimate = coef(numARG.anaerobe.adj), confint(numARG.anaerobe.adj))
anaerobe.num <- data.frame(exp(numARG.anaerobe.adj.model))
exp(numARG.anaerobe.adj.model)
anaerobe.num <- anaerobe.num %>% rownames_to_column(var = "row") %>% filter(row == "antibiotic.exposure1" | row == "antibiotic.exposure2"|row == "antibiotic.exposure3") %>% mutate(effect = round(Estimate, 2), lower= round(X2.5.., 2), upper = round(X97.5.., 2)) %>%
  mutate(drug = if_else(row == "antibiotic.exposure1", "Aerobic Antibiotics", if_else(row=="antibiotic.exposure2", "Anaerobic Antibiotics", if_else(row=="antibiotic.exposure3", "Aerobic and Anaerobic Antibiotics", "check")))) %>% dplyr::select(drug, effect, lower, upper)

#Number of new ARGs
newARG.anaerobe.adj <- glmmadmb(new.ARG ~ antibiotic.exposure + age + sex + diagnosis + prep + type_hsct + day + log.depth + (1 |study_id), data  = one.abx.test8, family = "nbinom1")
summary(newARG.anaerobe.adj)
newARG.anaerobe.adj.model <- cbind(Estimate = coef(newARG.anaerobe.adj), confint(newARG.anaerobe.adj))
anaerobe.new <- data.frame(exp(newARG.anaerobe.adj.model))
exp(newARG.anaerobe.adj.model)
anaerobe.new <- anaerobe.new %>% rownames_to_column(var = "row") %>% filter(row == "antibiotic.exposure1" | row == "antibiotic.exposure2"|row == "antibiotic.exposure3") %>% mutate(effect = round(Estimate, 2), lower= round(X2.5.., 2), upper = round(X97.5.., 2)) %>%
  mutate(drug = if_else(row == "antibiotic.exposure1", "Aerobic Antibiotics", if_else(row=="antibiotic.exposure2", "Anaerobic Antibiotics", if_else(row=="antibiotic.exposure3", "Aerobic and Anaerobic Antibiotics", "check")))) %>% dplyr::select(drug, effect, lower, upper)

#instability
inst.anaerobe.adj <- lmer(inst.score ~ antibiotic.exposure + age + sex + diagnosis + prep + type_hsct + day + log.depth + (1 |study_id), data  = one.abx.test8)
summary(inst.anaerobe.adj)
confint(inst.anaerobe.adj)
inst.anaerobe.exp.adj.model <- data.frame(cbind(coef(summary(inst.anaerobe.adj))), confint(inst.anaerobe.adj, parm = c("(Intercept)", "antibiotic.exposure1", "antibiotic.exposure2", "antibiotic.exposure3", "age" , "sexM" , "diagnosisCongenital_immunodef", "diagnosisHeme_non-malignancy", "diagnosisMetabolic_disease", "diagnosisSolid_tumor", "prep2" , "type_hsctAllogeneic" , "day" , "log.depth" ))) %>%
  rownames_to_column(var = "row") %>% filter(row == "antibiotic.exposure1" | row == "antibiotic.exposure2"|row == "antibiotic.exposure3") %>% mutate(effect = round(Estimate, 2), lower= round(X2.5.., 2), upper = round(X97.5.., 2)) %>%
  mutate(drug = if_else(row == "antibiotic.exposure1", "Aerobic Antibiotics", if_else(row=="antibiotic.exposure2", "Anaerobic Antibiotics", if_else(row=="antibiotic.exposure3", "Aerobic and Anaerobic Antibiotics", "check")))) %>% dplyr::select(drug, effect, lower, upper) %>%
dplyr::select(drug, effect, lower, upper)

#Correlation between microbiome and resistome measures
resistome <- test8 %>% dplyr::select(study_id, SampleID, day, num.ARG, new.ARG, inst.score) %>% dplyr::rename("inst.score.ARG" = "inst.score") 
composition <- meta_bmt %>% dplyr::select(study_id, SampleID, day, num.species, new.sp, inst.score) %>% dplyr::rename("inst.score.sp" = "inst.score")
res.comp <- left_join(composition, resistome, by= "SampleID") %>% dplyr::rename("day" = "day.x", "study_id" = "study_id.x") %>% dplyr::select(-c(day.y, study_id.y))
plot.df.c <- composition %>% dplyr::select(-c(SampleID, study_id)) %>% mutate(Source = "Composition") %>% dplyr::rename("New" = "new.sp", "Num" = "num.species", "inst.score" = "inst.score.sp")
plot.df.r <- resistome %>% dplyr::select(-c(SampleID, study_id)) %>% mutate(Source = "Resistome") %>% dplyr::rename("New" = "new.ARG", "Num" = "num.ARG", "inst.score" = "inst.score.ARG")
plot.df <- rbind(plot.df.c, plot.df.r)
plot.df$Source <- as.factor(plot.df$Source)

#Figure 3
colors <- c("Resistome" = "blue3", "Composition" = "red3")

summary(plot.df$Num)
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
  scale_y_continuous(limits = c(0, 250), breaks = c(0, 50, 100, 150, 200, 250), expand = c(0,0))+
  scale_x_continuous(limits = c(-30, 101), breaks = c(-30, 0, 30, 60, 100), expand = c(0,0))
plot1.2

summary(plot.df$New)
table(plot.df$Source)
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
  scale_y_continuous(limits = c(0, 180), breaks = c(0, 30, 60, 90, 120, 150, 180), expand = c(0,0))+
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
rmc.out.num 

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
  annotate(geom = "text", x=25, y= 225, label = expression('r'[rm]*'=0.47 (0.40, 0.53)'), size = 6)+
  scale_x_continuous(limits = c(0, 125), expand = c(0,0), breaks= c(0, 25, 50, 75, 100, 125))+
  scale_y_continuous(limits = c(0,250), breaks = c(0, 50, 100, 150, 200, 250), expand = c(0,0))+
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
rmc.out.new 

res.comp.new <- ggplot(res.comp, aes(x = new.ARG, y = new.sp)) +
  geom_point()+
  theme_classic()+
  geom_smooth(method = "lm", se = FALSE, size = 2)+
  scale_y_continuous(limits =c(0, 180), breaks = c(0, 30, 60, 90 ,120 , 150, 180), expand = c(0,0))+
  scale_x_continuous(limits =c(0, 100), breaks = c(0, 20, 40, 60 ,80 ,100), expand = c(0,0))+
  labs(x = "Number of New ARGs", y = "Number of New Species")+
  annotate(geom = "text", x=20, y= 150, label = expression('r'[rm]*'=0.50 (0.44, 0.56)'), size = 6)+
  theme(text = element_text(size = 20, color = "black"), plot.title = element_text(size = 20, hjust = 0.5),
        axis.text = element_text(color = "black"),
        plot.margin = unit(c(0.5,1,0.5,0.5), "cm"))
res.comp.new

#instability
rmc.out.inst.score <- rmcorr(study_id,
                             inst.score.ARG, inst.score.sp, res.comp, CI.level = 0.95, CIs =
                               c("analytic", "bootstrap"), nreps = 100,
                             bstrap.out = F)
rmc.out.inst.score 

res.comp.inst.score <- ggplot(res.comp, aes(x = inst.score.ARG, y = inst.score.sp)) +
  geom_point()+
  theme_classic()+
  geom_smooth(method = "lm", se = FALSE, size = 2) +
  scale_y_continuous(limits = c(0, 1.05), breaks = c(0, 0.25, 0.5, 0.75, 1.0), labels = c("0", "0.25", "0.5", "0.75", "1.0"), expand = c(0,0))+
  scale_x_continuous(limits = c(0, 1.05), breaks = c(0, 0.25, 0.5, 0.75, 1.0), labels = c("0", "0.25", "0.5", "0.75", "1.0"), expand = c(0,0))+
  labs(x = "Jaccard Distance of ARGs", y = "Jaccard Distance of Species")+
  annotate(geom = "text", x=0.22, y= 1.02, label = expression('r'[rm]*'=0.56 (0.50, 0.62)'), size = 6)+
  theme(text = element_text(size = 20, color = "black"), plot.title = element_text(size = 20, hjust = 0.1),
        axis.text = element_text(color = "black"),
        plot.margin = unit(c(0.5, 1, 0.5, 0.5), "cm"))
res.comp.inst.score

#arranging Figure 3
ggarrange(plot1.2, res.comp.num, plot1.4, res.comp.new, plot2.2, res.comp.inst.score,
          ncol = 2, nrow = 3,
          labels = c("A","", "B", "", "C", ""), 
          vjust = 1,
          font.label = list(size = 24, color = "black"))+
  theme(plot.margin = margin())

#####CLINICAL OUTCOMES FIGURES######
outcomes <- meta.bmt %>% dplyr::select(study_id,  agvhd_gutliver, mortality_2y) %>%
  distinct(study_id, .keep_all=TRUE) %>%
  dplyr::rename(gvhd=agvhd_gutliver, mortality=mortality_2y)
outcomes$study_id <- factor(outcomes$study_id)
test8.outcomes<- left_join(test8, outcomes, by="study_id")
test8.outcomes <- test8.outcomes %>% dplyr::select(study_id, SampleID, day, num.ARG, new.ARG, inst.score, total.abund, gvhd, bsi, mortality)
colnames(test8.outcomes)

test8.outcomes$gvhd<- as.factor(test8.outcomes$gvhd)
test8.outcomes$bsi<- as.factor(test8.outcomes$bsi)
test8.outcomes$mortality<- as.factor(test8.outcomes$mortality)
class(test8.outcomes$mortality)

library(ggbreak)
plot.gvhd.abund <- ggplot(data = test8.outcomes, aes(x = day, y = total.abund, colour = gvhd))+
  geom_point(aes(x=day, y=total.abund),
             color = "blue",
             data = filter(test8.outcomes, gvhd == 1)) +
  geom_point(aes(x = day, y=total.abund),
             color = "darkorange",
             data = filter(test8.outcomes, gvhd == 0))+
  geom_smooth(aes(x=day, y=total.abund),
              color = "blue4",
              fill = "steelblue1",
              data = filter(test8.outcomes, gvhd == 1))+
  geom_smooth(aes(x = day, y = total.abund),
              color = "darkorange3",
              fill = "peachpuff2",              
              data = filter(test8.outcomes, gvhd == 0))+
  theme_classic()+
  labs(y="ARG Abundance", x="Day Relative to HCT", color = "GVHD Status")+
  scale_y_continuous(limits = c(0, 100100), breaks = c(0, 10000, 20000, 30000, 40000, 50000, 90000, 100000), labels = c(0, 10000, 20000, 30000, 40000, 50000, 90000, "100000"), expand = c(0,0))+
  scale_y_break(c(50000, 89100))+
  scale_x_continuous(limits = c(-30, 101), breaks = c(-30, 0, 30, 60, 100), expand = c(0,0))+
  theme(text = element_text(size = 20),
        axis.text = element_text(color = "black"),
        plot.margin = unit(c(0.5,0.5,0.5,1), "cm"),
        axis.text.y.right = element_blank(),
        axis.line.y.right = element_blank(),
        axis.ticks.y.right = element_blank())
plot.gvhd.abund
#blue = yes gvhd; orange = no gvhd

plot.bsi.abund <- ggplot(data = test8.outcomes, aes(x = day, y = total.abund, colour = bsi))+
  geom_point(aes(x=day, y=total.abund),
             color = "blue",
             data = filter(test8.outcomes, bsi ==1 | bsi ==2)) +
  geom_point(aes(x = day, y=total.abund),
             color = "darkorange",
             data = filter(test8.outcomes, bsi == 0))+
  geom_smooth(aes(x=day, y=total.abund),
              color = "blue4",
              fill = "steelblue1",
              data = filter(test8.outcomes, bsi ==1 | bsi ==2))+
  geom_smooth(aes(x = day, y = total.abund),
              color = "darkorange3",
              fill = "peachpuff2",              
              data = filter(test8.outcomes, bsi == 0))+
  theme_classic()+
  labs(y="ARG Abundance", x="Day Relative to HCT", color = "BSI Status")+
  scale_y_continuous(limits = c(0, 100100), breaks = c(0, 10000, 20000, 30000, 40000, 50000, 90000, 100000), labels = c(0, 10000, 20000, 30000, 40000, 50000, 90000, "100000"), expand = c(0,0))+
  scale_y_break(c(50000, 89100))+
  scale_x_continuous(limits = c(-30, 101), breaks = c(-30, 0, 30, 60, 100), expand = c(0,0))+
  theme(text = element_text(size = 20),
        axis.text = element_text(color = "black"),
        plot.margin = unit(c(0.5,0.5,0.5,1), "cm"),
        axis.text.y.right = element_blank(),
        axis.line.y.right = element_blank(),
        axis.ticks.y.right = element_blank())
plot.bsi.abund
#blue = yes bsi; orange = no bsi

plot.mort.abund <- ggplot(data = test8.outcomes, aes(x = day, y = total.abund, colour = mortality))+
  geom_point(aes(x=day, y=total.abund),
             color = "blue",
             data = filter(test8.outcomes, mortality == 1)) +
  geom_point(aes(x = day, y=total.abund),
             color = "darkorange",
             data = filter(test8.outcomes, mortality == 0))+
  geom_smooth(aes(x=day, y=total.abund),
              color = "blue4",
              fill = "steelblue1",
              data = filter(test8.outcomes, mortality == 1))+
  geom_smooth(aes(x = day, y = total.abund),
              color = "darkorange3",
              fill = "peachpuff2",              
              data = filter(test8.outcomes, mortality == 0))+
  theme_classic()+
  labs(y="ARG Abundance", x="Day Relative to HCT", color = "GVHD Status")+
  scale_y_continuous(limits = c(0, 100100), breaks = c(0, 10000, 20000, 30000, 40000, 50000, 90000, 100000), labels = c(0, 10000, 20000, 30000, 40000, 50000, 90000, "100000"), expand = c(0,0))+
  scale_y_break(c(50000, 89100))+
  scale_x_continuous(limits = c(-30, 101), breaks = c(-30, 0, 30, 60, 100), expand = c(0,0))+
  theme(text = element_text(size = 20),
        axis.text = element_text(color = "black"),
        plot.margin = unit(c(0.5,0.5,0.5,1), "cm"),
        axis.text.y.right = element_blank(),
        axis.line.y.right = element_blank(),
        axis.ticks.y.right = element_blank())
plot.mort.abund
#blue = yes mortality; orange = no mortality

ggarrange(print(plot.gvhd.abund), print(plot.bsi.abund), print(plot.mort.abund),
          ncol=1,
          labels = c("A","B","C"), 
          vjust = 1,
          font.label = list(size = 24, color = "black"))+
  theme(plot.margin = margin())
          
#####

#sensitivity testing of individual antibiotics (Supplementary Table 3)
#Number of ARGs
numARG.cefepimeexp.adj <- glmmadmb(num.ARG ~ cefepime.exp + age + sex + diagnosis + prep + type_hsct + day + log.depth + (1 |study_id), data  = test8, family = "nbinom1")
summary(numARG.cefepimeexp.adj)
numARG.cefepimeexp.adj.model <- cbind(Estimate = coef(numARG.cefepimeexp.adj), confint(numARG.cefepimeexp.adj))
cefepime.num <- data.frame(exp(numARG.cefepimeexp.adj.model))
exp(numARG.cefepimeexp.adj.model)
cefepime.num <- cefepime.num %>% rownames_to_column(var = "row") %>% filter(row == "cefepime.exp1") %>% mutate(effect = round(Estimate, 2), lower= round(X2.5.., 2), upper = round(X97.5.., 2)) %>%
  mutate(drug = if_else(row == "cefepime.exp1", "Cefepime", "one")) %>% dplyr::select(drug, effect, lower, upper)

numARG.vancexp.adj <- glmmadmb(num.ARG ~ vanc.exp + age + sex + diagnosis + prep + type_hsct + day + log.depth + (1 |study_id), data  = test8, family = "nbinom1")
summary(numARG.vancexp.adj)
numARG.vancexp.adj.model <- cbind(Estimate = coef(numARG.vancexp.adj), confint(numARG.vancexp.adj))
vanco.num <- data.frame(exp(numARG.vancexp.adj.model))
exp(numARG.vancexp.adj.model)
vanco.num <- vanco.num %>% rownames_to_column(var = "row") %>% filter(row == "vanc.exp1") %>% mutate(effect = round(Estimate, 2), lower= round(X2.5.., 2), upper = round(X97.5.., 2)) %>%
  mutate(drug = if_else(row == "vanc.exp1", "Vancomycin", "one")) %>% dplyr::select(drug, effect, lower, upper)

numARG.quinoloneexp.adj <- glmmadmb(num.ARG ~ quinolone.exp + age + sex + diagnosis + prep + type_hsct + day + log.depth + (1 |study_id), data  = test8, family = "nbinom1")
summary(numARG.quinoloneexp.adj)
numARG.quinoloneexp.adj.model <- cbind(Estimate = coef(numARG.quinoloneexp.adj), confint(numARG.quinoloneexp.adj))
quinolone.num <- data.frame(exp(numARG.quinoloneexp.adj.model))
exp(numARG.quinoloneexp.adj.model)
quinolone.num <- quinolone.num %>% rownames_to_column(var = "row") %>% filter(row == "quinolone.exp1") %>% mutate(effect = round(Estimate, 2), lower= round(X2.5.., 2), upper = round(X97.5.., 2)) %>%
  mutate(drug = if_else(row == "quinolone.exp1", "Fluoroquinolone", "one")) %>% dplyr::select(drug, effect, lower, upper)

numARG.agexp.adj <- glmmadmb(num.ARG ~ ag.exp + age + sex + diagnosis + prep + type_hsct + day + log.depth + (1 |study_id), data  = test8, family = "nbinom1")
summary(numARG.agexp.adj)
numARG.agexp.adj.model <- cbind(Estimate = coef(numARG.agexp.adj), confint(numARG.agexp.adj))
ag.num <- data.frame(exp(numARG.agexp.adj.model))
exp(numARG.agexp.adj.model)
ag.num <- ag.num %>% rownames_to_column(var = "row") %>% filter(row == "ag.exp1") %>% mutate(effect = round(Estimate, 2), lower= round(X2.5.., 2), upper = round(X97.5.., 2)) %>%
  mutate(drug = if_else(row == "ag.exp1", "Aminoglycoside", "one")) %>% dplyr::select(drug, effect, lower, upper)

numARG.macrolideexp.adj <- glmmadmb(num.ARG ~ macrolide.exp + age + sex + diagnosis + prep + type_hsct + day + log.depth + (1 |study_id), data  = test8, family = "nbinom1")
summary(numARG.macrolideexp.adj)
numARG.macrolideexp.adj.model <- cbind(Estimate = coef(numARG.macrolideexp.adj), confint(numARG.macrolideexp.adj))
macrolide.num <- data.frame(exp(numARG.macrolideexp.adj.model))
exp(numARG.macrolideexp.adj.model)
macrolide.num <- macrolide.num %>% rownames_to_column(var = "row") %>% filter(row == "macrolide.exp1") %>% mutate(effect = round(Estimate, 2), lower= round(X2.5.., 2), upper = round(X97.5.., 2)) %>%
  mutate(drug = if_else(row == "macrolide.exp1", "Macrolide", "one")) %>% dplyr::select(drug, effect, lower, upper)

numARG.piptazoexp.adj <- glmmadmb(num.ARG ~ piptazo.exp + age + sex + diagnosis + prep + type_hsct + day + log.depth + (1 |study_id), data  = test8, family = "nbinom1")
summary(numARG.piptazoexp.adj)
numARG.piptazoexp.adj.model <- cbind(Estimate = coef(numARG.piptazoexp.adj), confint(numARG.piptazoexp.adj))
piptazo.num <- data.frame(exp(numARG.piptazoexp.adj.model))
exp(numARG.piptazoexp.adj.model)
piptazo.num <- piptazo.num %>% rownames_to_column(var = "row") %>% filter(row == "piptazo.exp1") %>% mutate(effect = round(Estimate, 2), lower= round(X2.5.., 2), upper = round(X97.5.., 2)) %>%
  mutate(drug = if_else(row == "piptazo.exp1", "Pip-Tazo", "one")) %>% dplyr::select(drug, effect, lower, upper)

numARG.carbapenemexp.adj <- glmmadmb(num.ARG ~ carbapenem.exp + age + sex + diagnosis + prep + type_hsct + day + log.depth + (1 |study_id), data  = test8, family = "nbinom1")
summary(numARG.carbapenemexp.adj)
numARG.carbapenemexp.adj.model <- cbind(Estimate = coef(numARG.carbapenemexp.adj), confint(numARG.carbapenemexp.adj))
carbapenem.num <- data.frame(exp(numARG.carbapenemexp.adj.model))
exp(numARG.carbapenemexp.adj.model)
carbapenem.num <- carbapenem.num %>% rownames_to_column(var = "row") %>% filter(row == "carbapenem.exp1") %>% mutate(effect = round(Estimate, 2), lower= round(X2.5.., 2), upper = round(X97.5.., 2)) %>%
  mutate(drug = if_else(row == "carbapenem.exp1", "Carbapenem", "one")) %>% dplyr::select(drug, effect, lower, upper)

numARG.flagylexp.adj <- glmmadmb(num.ARG ~ flagyl.exp + age + sex + diagnosis + prep + type_hsct + day + log.depth + (1 |study_id), data  = test8, family = "nbinom1")
summary(numARG.flagylexp.adj)
numARG.flagylexp.adj.model <- cbind(Estimate = coef(numARG.flagylexp.adj), confint(numARG.flagylexp.adj))
flagyl.num <- data.frame(exp(numARG.flagylexp.adj.model))
exp(numARG.flagylexp.adj.model)
flagyl.num <- flagyl.num %>% rownames_to_column(var = "row") %>% filter(row == "flagyl.exp1") %>% mutate(effect = round(Estimate, 2), lower= round(X2.5.., 2), upper = round(X97.5.., 2)) %>%
  mutate(drug = if_else(row == "flagyl.exp1", "Metronidazole", "one")) %>% dplyr::select(drug, effect, lower, upper)

numARG.clindaexp.adj <- glmmadmb(num.ARG ~ clinda.exp + age + sex + diagnosis + prep + type_hsct + day + log.depth + (1 |study_id), data  = test8, family = "nbinom1")
summary(numARG.clindaexp.adj)
numARG.clindaexp.adj.model <- cbind(Estimate = coef(numARG.clindaexp.adj), confint(numARG.clindaexp.adj))
clinda.num <- data.frame(exp(numARG.clindaexp.adj.model))
exp(numARG.clindaexp.adj.model)
clinda.num <- clinda.num %>% rownames_to_column(var = "row") %>% filter(row == "clinda.exp1") %>% mutate(effect = round(Estimate, 2), lower= round(X2.5.., 2), upper = round(X97.5.., 2)) %>%
  mutate(drug = if_else(row == "clinda.exp1", "Clindamycin", "one")) %>% dplyr::select(drug, effect, lower, upper)

numARG.tmp_smxexp.adj <- glmmadmb(num.ARG ~ tmp_smx.exp + age + sex + diagnosis + prep + type_hsct + day + log.depth + (1 |study_id), data  = test8, family = "nbinom1")
summary(numARG.tmp_smxexp.adj)
numARG.tmp_smxexp.adj.model <- cbind(Estimate = coef(numARG.tmp_smxexp.adj), confint(numARG.tmp_smxexp.adj))
tmp_smx.num <- data.frame(exp(numARG.tmp_smxexp.adj.model))
exp(numARG.tmp_smxexp.adj.model)
tmp_smx.num <- tmp_smx.num %>% rownames_to_column(var = "row") %>% filter(row == "tmp_smx.exp1") %>% mutate(effect = round(Estimate, 2), lower= round(X2.5.., 2), upper = round(X97.5.., 2)) %>%
  mutate(drug = if_else(row == "tmp_smx.exp1", "TMP-SMX", "one")) %>% dplyr::select(drug, effect, lower, upper)


#New ARGs
newARG.cefepimeexp.adj <- glmmadmb(new.ARG ~ cefepime.exp + age + sex + diagnosis + prep + type_hsct + day + log.depth + (1 |study_id), data  = test8, family = "nbinom1")
summary(newARG.cefepimeexp.adj)
newARG.cefepimeexp.adj.model <- cbind(Estimate = coef(newARG.cefepimeexp.adj), confint(newARG.cefepimeexp.adj))
cefepime.new <- data.frame(exp(newARG.cefepimeexp.adj.model))
exp(newARG.cefepimeexp.adj.model)
cefepime.new <- cefepime.new %>% rownames_to_column(var = "row") %>% filter(row == "cefepime.exp1") %>% mutate(effect = round(Estimate, 2), lower= round(X2.5.., 2), upper = round(X97.5.., 2)) %>%
  mutate(drug = if_else(row == "cefepime.exp1", "Cefepime", "one")) %>% dplyr::select(drug, effect, lower, upper)

newARG.vancexp.adj <- glmmadmb(new.ARG ~ vanc.exp + age + sex + diagnosis + prep + type_hsct + day + log.depth + (1 |study_id), data  = test8, family = "nbinom1")
summary(newARG.vancexp.adj)
newARG.vancexp.adj.model <- cbind(Estimate = coef(newARG.vancexp.adj), confint(newARG.vancexp.adj))
vanc.new <- data.frame(exp(newARG.vancexp.adj.model))
exp(newARG.vancexp.adj.model)
vanc.new <- vanc.new %>% rownames_to_column(var = "row") %>% filter(row == "vanc.exp1") %>% mutate(effect = round(Estimate, 2), lower= round(X2.5.., 2), upper = round(X97.5.., 2)) %>%
  mutate(drug = if_else(row == "vanc.exp1", "Vancomycin", "one")) %>% dplyr::select(drug, effect, lower, upper)

newARG.quinoloneexp.adj <- glmmadmb(new.ARG ~ quinolone.exp + age + sex + diagnosis + prep + type_hsct + day + log.depth + (1 |study_id), data  = test8, family = "nbinom1")
summary(newARG.quinoloneexp.adj)
newARG.quinoloneexp.adj.model <- cbind(Estimate = coef(newARG.quinoloneexp.adj), confint(newARG.quinoloneexp.adj))
quinolone.new <- data.frame(exp(newARG.quinoloneexp.adj.model))
exp(newARG.quinoloneexp.adj.model)
quinolone.new <- quinolone.new %>% rownames_to_column(var = "row") %>% filter(row == "quinolone.exp1") %>% mutate(effect = round(Estimate, 2), lower= round(X2.5.., 2), upper = round(X97.5.., 2)) %>%
  mutate(drug = if_else(row == "quinolone.exp1", "Fluoroquinolone", "one")) %>% dplyr::select(drug, effect, lower, upper)

newARG.agexp.adj <- glmmadmb(new.ARG ~ ag.exp + age + sex + diagnosis + prep + type_hsct + day + log.depth + (1 |study_id), data  = test8, family = "nbinom1")
summary(newARG.agexp.adj)
newARG.agexp.adj.model <- cbind(Estimate = coef(newARG.agexp.adj), confint(newARG.agexp.adj))
ag.new <- data.frame(exp(newARG.agexp.adj.model))
exp(newARG.agexp.adj.model)
ag.new <- ag.new %>% rownames_to_column(var = "row") %>% filter(row == "ag.exp1") %>% mutate(effect = round(Estimate, 2), lower= round(X2.5.., 2), upper = round(X97.5.., 2)) %>%
  mutate(drug = if_else(row == "ag.exp1", "Aminoglycoside", "one")) %>% dplyr::select(drug, effect, lower, upper)

newARG.macrolideexp.adj <- glmmadmb(new.ARG ~ macrolide.exp + age + sex + diagnosis + prep + type_hsct + day + log.depth + (1 |study_id), data  = test8, family = "nbinom1")
summary(newARG.macrolideexp.adj)
newARG.macrolideexp.adj.model <- cbind(Estimate = coef(newARG.macrolideexp.adj), confint(newARG.macrolideexp.adj))
macrolide.new <- data.frame(exp(newARG.macrolideexp.adj.model))
exp(newARG.macrolideexp.adj.model)
macrolide.new <- macrolide.new %>% rownames_to_column(var = "row") %>% filter(row == "macrolide.exp1") %>% mutate(effect = round(Estimate, 2), lower= round(X2.5.., 2), upper = round(X97.5.., 2)) %>%
  mutate(drug = if_else(row == "macrolide.exp1", "Macrolide", "one")) %>% dplyr::select(drug, effect, lower, upper)

newARG.tmp_smxexp.adj <- glmmadmb(new.ARG ~ tmp_smx.exp + age + sex + diagnosis + prep + type_hsct + day + log.depth + (1 |study_id), data  = test8, family = "nbinom1")
summary(newARG.tmp_smxexp.adj)
newARG.tmp_smxexp.adj.model <- cbind(Estimate = coef(newARG.tmp_smxexp.adj), confint(newARG.tmp_smxexp.adj))
tmp_smx.new <- data.frame(exp(newARG.tmp_smxexp.adj.model))
exp(newARG.tmp_smxexp.adj.model)
tmp_smx.new <- tmp_smx.new %>% rownames_to_column(var = "row") %>% filter(row == "tmp_smx.exp1") %>% mutate(effect = round(Estimate, 2), lower= round(X2.5.., 2), upper = round(X97.5.., 2)) %>%
  mutate(drug = if_else(row == "tmp_smx.exp1", "TMP-SMX", "one")) %>% dplyr::select(drug, effect, lower, upper)

newARG.piptazoexp.adj <- glmmadmb(new.ARG ~ piptazo.exp + age + sex + diagnosis + prep + type_hsct + day + log.depth + (1 |study_id), data  = test8, family = "nbinom1")
summary(newARG.piptazoexp.adj)
newARG.piptazoexp.adj.model <- cbind(Estimate = coef(newARG.piptazoexp.adj), confint(newARG.piptazoexp.adj))
piptazo.new <- data.frame(exp(newARG.piptazoexp.adj.model))
exp(newARG.piptazoexp.adj.model)
piptazo.new <- piptazo.new %>% rownames_to_column(var = "row") %>% filter(row == "piptazo.exp1") %>% mutate(effect = round(Estimate, 2), lower= round(X2.5.., 2), upper = round(X97.5.., 2)) %>%
  mutate(drug = if_else(row == "piptazo.exp1", "Pip-Tazo", "one")) %>% dplyr::select(drug, effect, lower, upper)

newARG.carbapenemexp.adj <- glmmadmb(new.ARG ~ carbapenem.exp + age + sex + diagnosis + prep + type_hsct + day + log.depth + (1 |study_id), data  = test8, family = "nbinom1")
summary(newARG.carbapenemexp.adj)
newARG.carbapenemexp.adj.model <- cbind(Estimate = coef(newARG.carbapenemexp.adj), confint(newARG.carbapenemexp.adj))
carbapenem.new <- data.frame(exp(newARG.carbapenemexp.adj.model))
exp(newARG.carbapenemexp.adj.model)
carbapenem.new <- carbapenem.new %>% rownames_to_column(var = "row") %>% filter(row == "carbapenem.exp1") %>% mutate(effect = round(Estimate, 2), lower= round(X2.5.., 2), upper = round(X97.5.., 2)) %>%
  mutate(drug = if_else(row == "carbapenem.exp1", "Carbapenem", "one")) %>% dplyr::select(drug, effect, lower, upper)

newARG.flagylexp.adj <- glmmadmb(new.ARG ~ flagyl.exp + age + sex + diagnosis + prep + type_hsct + day + log.depth + (1 |study_id), data  = test8, family = "nbinom1")
summary(newARG.flagylexp.adj)
newARG.flagylexp.adj.model <- cbind(Estimate = coef(newARG.flagylexp.adj), confint(newARG.flagylexp.adj))
flagyl.new <- data.frame(exp(newARG.flagylexp.adj.model))
exp(newARG.flagylexp.adj.model)
flagyl.new <- flagyl.new %>% rownames_to_column(var = "row") %>% filter(row == "flagyl.exp1") %>% mutate(effect = round(Estimate, 2), lower= round(X2.5.., 2), upper = round(X97.5.., 2)) %>%
  mutate(drug = if_else(row == "flagyl.exp1", "Metronidazole", "one")) %>% dplyr::select(drug, effect, lower, upper)

newARG.clindaexp.adj <- glmmadmb(new.ARG ~ clinda.exp + age + sex + diagnosis + prep + type_hsct + day + log.depth + (1 |study_id), data  = test8, family = "nbinom1")
summary(newARG.clindaexp.adj)
newARG.clindaexp.adj.model <- cbind(Estimate = coef(newARG.clindaexp.adj), confint(newARG.clindaexp.adj))
clinda.new <- data.frame(exp(newARG.clindaexp.adj.model))
exp(newARG.clindaexp.adj.model)
clinda.new <- clinda.new %>% rownames_to_column(var = "row") %>% filter(row == "clinda.exp1") %>% mutate(effect = round(Estimate, 2), lower= round(X2.5.., 2), upper = round(X97.5.., 2)) %>%
  mutate(drug = if_else(row == "clinda.exp1", "Clindamycin", "one")) %>% dplyr::select(drug, effect, lower, upper)

#Instability
inst.cefepimeexp.adj <- lmer(inst.score ~ cefepime.exp + age + sex + diagnosis + prep + type_hsct + day + log.depth + (1 |study_id), data  = test8)
summary(inst.cefepimeexp.adj)
confint(inst.cefepimeexp.adj)
inst.cefepime.exp.adj.model <- data.frame(cbind(coef(summary(inst.cefepimeexp.adj))), confint(inst.cefepimeexp.adj, parm = c("(Intercept)", "cefepime.exp1", "age" , "sexM" , "diagnosisCongenital_immunodef", "diagnosisHeme_non-malignancy", "diagnosisMetabolic_disease", "diagnosisSolid_tumor", "prep2" , "type_hsctAllogeneic" , "day" , "log.depth" ))) %>%
  rownames_to_column(var = "row") %>% filter(row == "cefepime.exp1") %>%
  mutate(effect = round(Estimate, 3), lower = round(X2.5.., 3), upper = round(X97.5.., 3), drug = if_else(row == "cefepime.exp1", "Cefepime", "")) %>%
  dplyr::select(drug, effect, lower, upper)

inst.vancexp.adj <- lmer(inst.score ~ vanc.exp + age + sex + diagnosis + prep + type_hsct + day + log.depth + (1 |study_id), data  = test8)
summary(inst.vancexp.adj)
confint(inst.vancexp.adj)
inst.vanc.exp.adj.model <- data.frame(cbind(coef(summary(inst.vancexp.adj))), confint(inst.vancexp.adj, parm = c("(Intercept)", "vanc.exp1", "age" , "sexM" , "diagnosisCongenital_immunodef", "diagnosisHeme_non-malignancy", "diagnosisMetabolic_disease", "diagnosisSolid_tumor", "prep2" , "type_hsctAllogeneic" , "day" , "log.depth" ))) %>%
  rownames_to_column(var = "row") %>% filter(row == "vanc.exp1") %>%
  mutate(effect = round(Estimate, 3), lower = round(X2.5.., 3), upper = round(X97.5.., 3), drug = if_else(row == "vanc.exp1", "Vancomycin", "")) %>%
  dplyr::select(drug, effect, lower, upper)

inst.quinoloneexp.adj <- lmer(inst.score ~ quinolone.exp + age + sex + diagnosis + prep + type_hsct + day + log.depth + (1 |study_id), data  = test8)
summary(inst.quinoloneexp.adj)
confint(inst.quinoloneexp.adj)
inst.quinolone.exp.adj.model <- data.frame(cbind(coef(summary(inst.quinoloneexp.adj))), confint(inst.quinoloneexp.adj, parm = c("(Intercept)", "quinolone.exp1", "age" , "sexM" , "diagnosisCongenital_immunodef", "diagnosisHeme_non-malignancy", "diagnosisMetabolic_disease", "diagnosisSolid_tumor", "prep2" , "type_hsctAllogeneic" , "day" , "log.depth" ))) %>%
  rownames_to_column(var = "row") %>% filter(row == "quinolone.exp1") %>%
  mutate(effect = round(Estimate, 3), lower = round(X2.5.., 3), upper = round(X97.5.., 3), drug = if_else(row == "quinolone.exp1", "Fluoroquinolone", "")) %>%
  dplyr::select(drug, effect, lower, upper)

inst.agexp.adj <- lmer(inst.score ~ ag.exp + age + sex + diagnosis + prep + type_hsct + day + log.depth + (1 |study_id), data  = test8)
summary(inst.agexp.adj)
confint(inst.agexp.adj)
inst.ag.exp.adj.model <- data.frame(cbind(coef(summary(inst.agexp.adj))), confint(inst.agexp.adj, parm = c("(Intercept)", "ag.exp1", "age" , "sexM" , "diagnosisCongenital_immunodef", "diagnosisHeme_non-malignancy", "diagnosisMetabolic_disease", "diagnosisSolid_tumor", "prep2" , "type_hsctAllogeneic" , "day" , "log.depth" ))) %>%
  rownames_to_column(var = "row") %>% filter(row == "ag.exp1") %>%
  mutate(effect = round(Estimate, 3), lower = round(X2.5.., 3), upper = round(X97.5.., 3), drug = if_else(row == "ag.exp1", "Aminoglycoside", "")) %>%
  dplyr::select(drug, effect, lower, upper)

inst.macrolideexp.adj <- lmer(inst.score ~ macrolide.exp + age + sex + diagnosis + prep + type_hsct + day + log.depth + (1 |study_id), data  = test8)
summary(inst.macrolideexp.adj)
confint(inst.macrolideexp.adj)
inst.macrolide.exp.adj.model <- data.frame(cbind(coef(summary(inst.macrolideexp.adj))), confint(inst.macrolideexp.adj, parm = c("(Intercept)", "macrolide.exp1", "age" , "sexM" , "diagnosisCongenital_immunodef", "diagnosisHeme_non-malignancy", "diagnosisMetabolic_disease", "diagnosisSolid_tumor", "prep2" , "type_hsctAllogeneic" , "day" , "log.depth" ))) %>%
  rownames_to_column(var = "row") %>% filter(row == "macrolide.exp1") %>%
  mutate(effect = round(Estimate, 3), lower = round(X2.5.., 3), upper = round(X97.5.., 3), drug = if_else(row == "macrolide.exp1", "Macrolide", "")) %>%
  dplyr::select(drug, effect, lower, upper)

inst.tmp_smxexp.adj <- lmer(inst.score ~ tmp_smx.exp + age + sex + diagnosis + prep + type_hsct + day + log.depth + (1 |study_id), data  = test8)
summary(inst.tmp_smxexp.adj)
confint(inst.tmp_smxexp.adj)
inst.tmp_smx.exp.adj.model <- data.frame(cbind(coef(summary(inst.tmp_smxexp.adj))), confint(inst.tmp_smxexp.adj, parm = c("(Intercept)", "tmp_smx.exp1", "age" , "sexM" , "diagnosisCongenital_immunodef", "diagnosisHeme_non-malignancy", "diagnosisMetabolic_disease", "diagnosisSolid_tumor", "prep2" , "type_hsctAllogeneic" , "day" , "log.depth" ))) %>%
  rownames_to_column(var = "row") %>% filter(row == "tmp_smx.exp1") %>%
  mutate(effect = round(Estimate, 3), lower = round(X2.5.., 3), upper = round(X97.5.., 3), drug = if_else(row == "tmp_smx.exp1", "TMP-SMX", "")) %>%
  dplyr::select(drug, effect, lower, upper)

inst.piptazoexp.adj <- lmer(inst.score ~ piptazo.exp + age + sex + diagnosis + prep + type_hsct + day + log.depth + (1 |study_id), data  = test8)
summary(inst.piptazoexp.adj)
confint(inst.piptazoexp.adj)
inst.piptazo.exp.adj.model <- data.frame(cbind(coef(summary(inst.piptazoexp.adj))), confint(inst.piptazoexp.adj, parm = c("(Intercept)", "piptazo.exp1", "age" , "sexM" , "diagnosisCongenital_immunodef", "diagnosisHeme_non-malignancy", "diagnosisMetabolic_disease", "diagnosisSolid_tumor", "prep2" , "type_hsctAllogeneic" , "day" , "log.depth" ))) %>%
  rownames_to_column(var = "row") %>% filter(row == "piptazo.exp1") %>%
  mutate(effect = round(Estimate, 3), lower = round(X2.5.., 3), upper = round(X97.5.., 3), drug = if_else(row == "piptazo.exp1", "Pip-Tazo", "")) %>%
  dplyr::select(drug, effect, lower, upper)

inst.carbapenemexp.adj <- lmer(inst.score ~ carbapenem.exp + age + sex + diagnosis + prep + type_hsct + day + log.depth + (1 |study_id), data  = test8)
summary(inst.carbapenemexp.adj)
confint(inst.carbapenemexp.adj)
inst.carbapenem.exp.adj.model <- data.frame(cbind(coef(summary(inst.carbapenemexp.adj))), confint(inst.carbapenemexp.adj, parm = c("(Intercept)", "carbapenem.exp1", "age" , "sexM" , "diagnosisCongenital_immunodef", "diagnosisHeme_non-malignancy", "diagnosisMetabolic_disease", "diagnosisSolid_tumor", "prep2" , "type_hsctAllogeneic" , "day" , "log.depth" ))) %>%
  rownames_to_column(var = "row") %>% filter(row == "carbapenem.exp1") %>%
  mutate(effect = round(Estimate, 3), lower = round(X2.5.., 3), upper = round(X97.5.., 3), drug = if_else(row == "carbapenem.exp1", "Carbapenem", "")) %>%
  dplyr::select(drug, effect, lower, upper)

inst.flagylexp.adj <- lmer(inst.score ~ flagyl.exp + age + sex + diagnosis + prep + type_hsct + day + log.depth + (1 |study_id), data  = test8)
summary(inst.flagylexp.adj)
confint(inst.flagylexp.adj)
inst.flagyl.exp.adj.model <- data.frame(cbind(coef(summary(inst.flagylexp.adj))), confint(inst.flagylexp.adj, parm = c("(Intercept)", "flagyl.exp1", "age" , "sexM" , "diagnosisCongenital_immunodef", "diagnosisHeme_non-malignancy", "diagnosisMetabolic_disease", "diagnosisSolid_tumor", "prep2" , "type_hsctAllogeneic" , "day" , "log.depth" ))) %>%
  rownames_to_column(var = "row") %>% filter(row == "flagyl.exp1") %>%
  mutate(effect = round(Estimate, 3), lower = round(X2.5.., 3), upper = round(X97.5.., 3), drug = if_else(row == "flagyl.exp1", "Metronidazole", "")) %>%
  dplyr::select(drug, effect, lower, upper)

inst.clindaexp.adj <- lmer(inst.score ~ clinda.exp + age + sex + diagnosis + prep + type_hsct + day + log.depth + (1 |study_id), data  = test8)
summary(inst.clindaexp.adj)
confint(inst.clindaexp.adj)
inst.clinda.exp.adj.model <- data.frame(cbind(coef(summary(inst.clindaexp.adj))), confint(inst.clindaexp.adj, parm = c("(Intercept)", "clinda.exp1", "age" , "sexM" , "diagnosisCongenital_immunodef", "diagnosisHeme_non-malignancy", "diagnosisMetabolic_disease", "diagnosisSolid_tumor", "prep2" , "type_hsctAllogeneic" , "day" , "log.depth" ))) %>%
  rownames_to_column(var = "row") %>% filter(row == "clinda.exp1") %>%
  mutate(effect = round(Estimate, 3), lower = round(X2.5.., 3), upper = round(X97.5.., 3), drug = if_else(row == "clinda.exp1", "Clindamycin", "")) %>%
  dplyr::select(drug, effect, lower, upper)

#ARG Abundance
gamma.cefepime <- glmer(nonzero.abund ~ cefepime.exp + age + sex + diagnosis + prep+ type_hsct + day + log.depth + (1 |study_id), data  = gamma, family=Gamma(link="log"))
summary(gamma.cefepime)
gamma.cefepime2 <- summary(gamma.cefepime)
gamma.cefepime2 <- as.data.frame(gamma.cefepime2$coefficients)
gamma.cef.com <- confint(gamma.cefepime, method="Wald")
gamma.cef.com <- gamma.cef.com[-c(1, 2),]
gamma.cefepime1 <- cbind(gamma.cefepime2, gamma.cef.com)
gamma.cefepime1 <- gamma.cefepime1[-1,]
gamma.cefepime1 <- gamma.cefepime1 %>% mutate(effect = round(exp(Estimate), digits =2), lower=round(exp(`2.5 %`), digits =2), upper=round(exp(`97.5 %`), digits =2)) %>%
  mutate(day.est = round(exp(Estimate), digits =3), day.lower = round(exp(`2.5 %`), digits =3), day.upper = round(exp(`97.5 %`), digits =4)) %>%
  dplyr::select(effect, lower, upper, day.est, day.lower, day.upper, `Pr(>|z|)`) %>%
  mutate(drug="Cefepime") %>%
  dplyr::select(drug, effect, lower, upper)
gamma.cefepime1 <- gamma.cefepime1[1,]

gamma.vanc <- glmer(nonzero.abund ~ vanc.exp + age + sex + diagnosis + prep+ type_hsct + day + log.depth + (1 |study_id), data  = gamma, family=Gamma(link="log"))
summary(gamma.vanc)
gamma.vanc2 <- summary(gamma.vanc)
gamma.vanc2 <- as.data.frame(gamma.vanc2$coefficients)
gamma.cef.com <- confint(gamma.vanc, method="Wald")
gamma.cef.com <- gamma.cef.com[-c(1, 2),]
gamma.vanc1 <- cbind(gamma.vanc2, gamma.cef.com)
gamma.vanc1 <- gamma.vanc1[-1,]
gamma.vanc1 <- gamma.vanc1 %>% mutate(effect = round(exp(Estimate), digits =2), lower=round(exp(`2.5 %`), digits =2), upper=round(exp(`97.5 %`), digits =2)) %>%
  mutate(day.est = round(exp(Estimate), digits =3), day.lower = round(exp(`2.5 %`), digits =3), day.upper = round(exp(`97.5 %`), digits =4)) %>%
  dplyr::select(effect, lower, upper, day.est, day.lower, day.upper, `Pr(>|z|)`)%>%
  mutate(drug="Vancomycin") %>%
  dplyr::select(drug, effect, lower, upper)
gamma.vanc1 <- gamma.vanc1[1,]

gamma.quinolone <- glmer(nonzero.abund ~ quinolone.exp + age + sex + diagnosis + prep+ type_hsct + day + log.depth + (1 |study_id), data  = gamma, family=Gamma(link="log"))
summary(gamma.quinolone)
gamma.quinolone2 <- summary(gamma.quinolone)
gamma.quinolone2 <- as.data.frame(gamma.quinolone2$coefficients)
gamma.cef.com <- confint(gamma.quinolone, method="Wald")
gamma.cef.com <- gamma.cef.com[-c(1, 2),]
gamma.quinolone1 <- cbind(gamma.quinolone2, gamma.cef.com)
gamma.quinolone1 <- gamma.quinolone1[-1,]
gamma.quinolone1 <- gamma.quinolone1 %>% mutate(effect = round(exp(Estimate), digits =2), lower=round(exp(`2.5 %`), digits =2), upper=round(exp(`97.5 %`), digits =2)) %>%
  mutate(day.est = round(exp(Estimate), digits =3), day.lower = round(exp(`2.5 %`), digits =3), day.upper = round(exp(`97.5 %`), digits =4)) %>%
  dplyr::select(effect, lower, upper, day.est, day.lower, day.upper, `Pr(>|z|)`)%>%
  mutate(drug="Fluoroquinolone") %>%
  dplyr::select(drug, effect, lower, upper)
gamma.quinolone1 <- gamma.quinolone1[1,]

gamma.ag <- glmer(nonzero.abund ~ ag.exp + age + sex + diagnosis + prep+ type_hsct + day + log.depth + (1 |study_id), data  = gamma, family=Gamma(link="log"))
summary(gamma.ag)
gamma.ag2 <- summary(gamma.ag)
gamma.ag2 <- as.data.frame(gamma.ag2$coefficients)
gamma.cef.com <- confint(gamma.ag, method="Wald")
gamma.cef.com <- gamma.cef.com[-c(1, 2),]
gamma.ag1 <- cbind(gamma.ag2, gamma.cef.com)
gamma.ag1 <- gamma.ag1[-1,]
gamma.ag1 <- gamma.ag1 %>% mutate(effect = round(exp(Estimate), digits =2), lower=round(exp(`2.5 %`), digits =2), upper=round(exp(`97.5 %`), digits =2)) %>%
  mutate(day.est = round(exp(Estimate), digits =3), day.lower = round(exp(`2.5 %`), digits =3), day.upper = round(exp(`97.5 %`), digits =4)) %>%
  dplyr::select(effect, lower, upper, day.est, day.lower, day.upper, `Pr(>|z|)`)%>%
  mutate(drug="Aminoglycoside") %>%
  dplyr::select(drug, effect, lower, upper)
gamma.ag1 <- gamma.ag1[1,]

gamma.macrolide <- glmer(nonzero.abund ~ macrolide.exp + age + sex + diagnosis + prep+ type_hsct + day + log.depth + (1 |study_id), data  = gamma, family=Gamma(link="log"))
summary(gamma.macrolide)
gamma.macrolide2 <- summary(gamma.macrolide)
gamma.macrolide2 <- as.data.frame(gamma.macrolide2$coefficients)
gamma.cef.com <- confint(gamma.macrolide, method="Wald")
gamma.cef.com <- gamma.cef.com[-c(1, 2),]
gamma.macrolide1 <- cbind(gamma.macrolide2, gamma.cef.com)
gamma.macrolide1 <- gamma.macrolide1[-1,]
gamma.macrolide1 <- gamma.macrolide1 %>% mutate(effect = round(exp(Estimate), digits =2), lower=round(exp(`2.5 %`), digits =2), upper=round(exp(`97.5 %`), digits =2)) %>%
  mutate(day.est = round(exp(Estimate), digits =3), day.lower = round(exp(`2.5 %`), digits =3), day.upper = round(exp(`97.5 %`), digits =4)) %>%
  dplyr::select(effect, lower, upper, day.est, day.lower, day.upper, `Pr(>|z|)`)%>%
  mutate(drug="Macrolide") %>%
  dplyr::select(drug, effect, lower, upper)
gamma.macrolide1 <- gamma.macrolide1[1,]

gamma.tmp_smx <- glmer(nonzero.abund ~ tmp_smx.exp + age + sex + diagnosis + prep+ type_hsct + day + log.depth + (1 |study_id), data  = gamma, family=Gamma(link="log"))
summary(gamma.tmp_smx)
gamma.tmp_smx2 <- summary(gamma.tmp_smx)
gamma.tmp_smx2 <- as.data.frame(gamma.tmp_smx2$coefficients)
gamma.cef.com <- confint(gamma.tmp_smx, method="Wald")
gamma.cef.com <- gamma.cef.com[-c(1, 2),]
gamma.tmp_smx1 <- cbind(gamma.tmp_smx2, gamma.cef.com)
gamma.tmp_smx1 <- gamma.tmp_smx1[-1,]
gamma.tmp_smx1 <- gamma.tmp_smx1 %>% mutate(effect = round(exp(Estimate), digits =2), lower=round(exp(`2.5 %`), digits =2), upper=round(exp(`97.5 %`), digits =2)) %>%
  mutate(day.est = round(exp(Estimate), digits =3), day.lower = round(exp(`2.5 %`), digits =3), day.upper = round(exp(`97.5 %`), digits =4)) %>%
  dplyr::select(effect, lower, upper, day.est, day.lower, day.upper, `Pr(>|z|)`)%>%
  mutate(drug="TMP-SMX") %>%
  dplyr::select(drug, effect, lower, upper)
gamma.tmp_smx1 <- gamma.tmp_smx1[1,]

gamma.piptazo <- glmer(nonzero.abund ~ piptazo.exp + age + sex + diagnosis + prep+ type_hsct + day + log.depth + (1 |study_id), data  = gamma, family=Gamma(link="log"))
summary(gamma.piptazo)
gamma.piptazo2 <- summary(gamma.piptazo)
gamma.piptazo2 <- as.data.frame(gamma.piptazo2$coefficients)
gamma.cef.com <- confint(gamma.piptazo, method="Wald")
gamma.cef.com <- gamma.cef.com[-c(1, 2),]
gamma.piptazo1 <- cbind(gamma.piptazo2, gamma.cef.com)
gamma.piptazo1 <- gamma.piptazo1[-1,]
gamma.piptazo1 <- gamma.piptazo1 %>% mutate(effect = round(exp(Estimate), digits =2), lower=round(exp(`2.5 %`), digits =2), upper=round(exp(`97.5 %`), digits =2)) %>%
  mutate(day.est = round(exp(Estimate), digits =3), day.lower = round(exp(`2.5 %`), digits =3), day.upper = round(exp(`97.5 %`), digits =4)) %>%
  dplyr::select(effect, lower, upper, day.est, day.lower, day.upper, `Pr(>|z|)`)%>%
  mutate(drug="Pip-Tazo") %>%
  dplyr::select(drug, effect, lower, upper)
gamma.piptazo1 <- gamma.piptazo1[1,]

gamma.carbapenem <- glmer(nonzero.abund ~ carbapenem.exp + age + sex + diagnosis + prep+ type_hsct + day + log.depth + (1 |study_id), data  = gamma, family=Gamma(link="log"))
summary(gamma.carbapenem)
gamma.carbapenem2 <- summary(gamma.carbapenem)
gamma.carbapenem2 <- as.data.frame(gamma.carbapenem2$coefficients)
gamma.cef.com <- confint(gamma.carbapenem, method="Wald")
gamma.cef.com <- gamma.cef.com[-c(1, 2),]
gamma.carbapenem1 <- cbind(gamma.carbapenem2, gamma.cef.com)
gamma.carbapenem1 <- gamma.carbapenem1[-1,]
gamma.carbapenem1 <- gamma.carbapenem1 %>% mutate(effect = round(exp(Estimate), digits =2), lower=round(exp(`2.5 %`), digits =2), upper=round(exp(`97.5 %`), digits =2)) %>%
  mutate(day.est = round(exp(Estimate), digits =3), day.lower = round(exp(`2.5 %`), digits =3), day.upper = round(exp(`97.5 %`), digits =4)) %>%
  dplyr::select(effect, lower, upper, day.est, day.lower, day.upper, `Pr(>|z|)`)%>%
  mutate(drug="Carbapenem") %>%
  dplyr::select(drug, effect, lower, upper)
gamma.carbapenem1 <- gamma.carbapenem1[1,]

gamma.flagyl <- glmer(nonzero.abund ~ flagyl.exp + age + sex + diagnosis + prep+ type_hsct + day + log.depth + (1 |study_id), data  = gamma, family=Gamma(link="log"))
summary(gamma.flagyl)
gamma.flagyl2 <- summary(gamma.flagyl)
gamma.flagyl2 <- as.data.frame(gamma.flagyl2$coefficients)
gamma.cef.com <- confint(gamma.flagyl, method="Wald")
gamma.cef.com <- gamma.cef.com[-c(1, 2),]
gamma.flagyl1 <- cbind(gamma.flagyl2, gamma.cef.com)
gamma.flagyl1 <- gamma.flagyl1[-1,]
gamma.flagyl1 <- gamma.flagyl1 %>% mutate(effect = round(exp(Estimate), digits =2), lower=round(exp(`2.5 %`), digits =2), upper=round(exp(`97.5 %`), digits =2)) %>%
  mutate(day.est = round(exp(Estimate), digits =3), day.lower = round(exp(`2.5 %`), digits =3), day.upper = round(exp(`97.5 %`), digits =4)) %>%
  dplyr::select(effect, lower, upper, day.est, day.lower, day.upper, `Pr(>|z|)`)%>%
  mutate(drug="Metronidazole") %>%
  dplyr::select(drug, effect, lower, upper)
gamma.flagyl1 <- gamma.flagyl1[1,]

gamma.clinda <- glmer(nonzero.abund ~ clinda.exp + age + sex + diagnosis + prep+ type_hsct + day + log.depth + (1 |study_id), data  = gamma, family=Gamma(link="log"))
summary(gamma.clinda)
gamma.clinda2 <- summary(gamma.clinda)
gamma.clinda2 <- as.data.frame(gamma.clinda2$coefficients)
gamma.cef.com <- confint(gamma.clinda, method="Wald")
gamma.cef.com <- gamma.cef.com[-c(1, 2),]
gamma.clinda1 <- cbind(gamma.clinda2, gamma.cef.com)
gamma.clinda1 <- gamma.clinda1[-1,]
gamma.clinda1 <- gamma.clinda1 %>% mutate(effect = round(exp(Estimate), digits =2), lower=round(exp(`2.5 %`), digits =2), upper=round(exp(`97.5 %`), digits =2)) %>%
  mutate(day.est = round(exp(Estimate), digits =3), day.lower = round(exp(`2.5 %`), digits =3), day.upper = round(exp(`97.5 %`), digits =4)) %>%
  dplyr::select(effect, lower, upper, day.est, day.lower, day.upper, `Pr(>|z|)`)%>%
  mutate(drug="Clindamycin") %>%
  dplyr::select(drug, effect, lower, upper)
gamma.clinda1 <- gamma.clinda1[1,]

#Figure 4
fp.num <- rbind(anaerobe.num, piptazo.num, carbapenem.num, flagyl.num, clinda.num, cefepime.num, vanco.num, quinolone.num, ag.num, macrolide.num, tmp_smx.num) %>%
  mutate(index = if_else(drug == "Aerobic Antibiotics", 12, if_else(drug == 'Cefepime', 11, if_else(drug== "Vancomycin", 10, if_else(drug == "Fluoroquinolone", 9, if_else(drug== "Aminoglycoside", 8, if_else(drug=="Macrolide", 7, if_else(drug=="TMP-SMX", 6, if_else(drug=="Anaerobic Antibiotics", 5, if_else(drug=="Pip-Tazo", 4, if_else(drug=="Carbapenem", 3, if_else(drug== "Metronidazole", 2, if_else(drug=="Clindamycin", 1, 0))))))))))))) %>%
  arrange(index) %>% mutate(color = if_else(drug == "Aerobic Antibiotics", "Main", if_else(drug == 'Cefepime', "Sub", if_else(drug== "Vancomycin", "Sub", if_else(drug == "Fluoroquinolone", "Sub", if_else(drug== "Aminoglycoside", "Sub", if_else(drug=="Macrolide", "Sub", if_else(drug=="TMP-SMX", "Sub", if_else(drug=="Anaerobic Antibiotics", "Main", if_else(drug=="Pip-Tazo", "Sub1", if_else(drug=="Carbapenem", "Sub1", if_else(drug== "Metronidazole", "Sub1", if_else(drug=="Aerobic and Anaerobic Antibiotics", "Main", "Sub1")))))))))))))
fp.num <- fp.num %>% mutate(index = index +1)
a <- ifelse(fp.num$color == "Main", "black", ifelse(fp.num$color== "Sub", "dodgerblue3", "forestgreen"))
b <- ifelse(fp.num$color == "Main", "bold", "plain")

num.h <- ggplot(data = fp.num, aes(y=index, x=effect, xmin=lower, xmax = upper)) +
  geom_point(color = a, size = 3.5) +
  geom_errorbarh(height = 0.2, color = a, size = 1) +
  scale_y_continuous(breaks = 1:13, labels = fp.num$drug, expand=c(0, 0.5)) +
  scale_x_continuous(limits = c(0.5, 1.5), breaks = c(0.5, 0.75, 1.0, 1.25, 1.5), labels = c("0.5", "0.75", "1.0", "1.25", "1.5")) +
  labs(title = "Number of ARGs", x = expression(beta)) +
  geom_vline(xintercept = 1, color = "black", linetype = "longdash", alpha = 0.5) +
  theme(panel.background = element_blank(),
        plot.title = element_text(size = 18, hjust = 0.5),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_blank(),
        axis.text.x = element_text(color = "black", size=12),
        axis.text.y = element_text(size = 15, color = a, face = b),
        axis.line.x = element_line(linetype = "solid"),
        panel.grid = element_blank(),
        axis.ticks.y = element_blank())
num.h

fp.new <- rbind(anaerobe.new, piptazo.new, carbapenem.new, flagyl.new, clinda.new, cefepime.new, vanc.new, quinolone.new, ag.new, macrolide.new, tmp_smx.new) %>%
  mutate(index = if_else(drug == "Aerobic Antibiotics", 12, if_else(drug == 'Cefepime', 11, if_else(drug== "Vancomycin", 10, if_else(drug == "Fluoroquinolone", 9, if_else(drug== "Aminoglycoside", 8, if_else(drug=="Macrolide", 7, if_else(drug=="TMP-SMX", 6, if_else(drug=="Anaerobic Antibiotics", 5, if_else(drug=="Pip-Tazo", 4, if_else(drug=="Carbapenem", 3, if_else(drug== "Metronidazole", 2, if_else(drug=="Clindamycin", 1, 0))))))))))))) %>%
  arrange(index) %>% mutate(color = if_else(drug == "Aerobic Antibiotics", "Main", if_else(drug == 'Cefepime', "Sub", if_else(drug== "Vancomycin", "Sub", if_else(drug == "Fluoroquinolone", "Sub", if_else(drug== "Aminoglycoside", "Sub", if_else(drug=="Macrolide", "Sub", if_else(drug=="TMP-SMX", "Sub", if_else(drug=="Anaerobic Antibiotics", "Main", if_else(drug=="Pip-Tazo", "Sub1", if_else(drug=="Carbapenem", "Sub1", if_else(drug== "Metronidazole", "Sub1", if_else(drug=="Aerobic and Anaerobic Antibiotics", "Main", "Sub1")))))))))))))
fp.new <- fp.new %>% mutate(index = index +1)

new.h <- ggplot(data = fp.new, aes(y=index, x=effect, xmin=lower, xmax = upper)) +
  geom_point(color = a, size = 3.5) +
  geom_errorbarh(height = 0.2, color = a, size = 1) +
  scale_y_continuous(breaks = 1:13, labels = fp.new$drug, expand=c(0, 0.5)) +
  scale_x_continuous(limits = c(0.5, 2.05), breaks = c(0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.05)) +
  labs(title = "Number of New ARGs", x = expression(beta)) +
  geom_vline(xintercept = 1, color = "black", linetype = "longdash", alpha = 0.5) +
  theme(panel.background = element_blank(),
        plot.title = element_text(size = 18, hjust = 0.5),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_blank(),
        axis.text.x = element_text(color = "black", size=12),
        axis.text.y = element_blank(),
        axis.line.x = element_line(linetype = "solid"),
        panel.grid = element_blank(),
        axis.ticks.y = element_blank())
new.h

fp.inst <- rbind(inst.anaerobe.exp.adj.model, inst.piptazo.exp.adj.model, inst.carbapenem.exp.adj.model, inst.flagyl.exp.adj.model, inst.clinda.exp.adj.model, inst.cefepime.exp.adj.model, inst.vanc.exp.adj.model, inst.quinolone.exp.adj.model, inst.ag.exp.adj.model, inst.macrolide.exp.adj.model, inst.tmp_smx.exp.adj.model) %>%
  mutate(index = if_else(drug == "Aerobic Antibiotics", 12, if_else(drug == 'Cefepime', 11, if_else(drug== "Vancomycin", 10, if_else(drug == "Fluoroquinolone", 9, if_else(drug== "Aminoglycoside", 8, if_else(drug=="Macrolide", 7, if_else(drug=="TMP-SMX", 6, if_else(drug=="Anaerobic Antibiotics", 5, if_else(drug=="Pip-Tazo", 4, if_else(drug=="Carbapenem", 3, if_else(drug== "Metronidazole", 2, if_else(drug=="Clindamycin", 1, 0))))))))))))) %>%
  arrange(index) %>% mutate(color = if_else(drug == "Aerobic Antibiotics", "Main", if_else(drug == 'Cefepime', "Sub", if_else(drug== "Vancomycin", "Sub", if_else(drug == "Fluoroquinolone", "Sub", if_else(drug== "Aminoglycoside", "Sub", if_else(drug=="Macrolide", "Sub", if_else(drug=="TMP-SMX", "Sub", if_else(drug=="Anaerobic Antibiotics", "Main", if_else(drug=="Pip-Tazo", "Sub1", if_else(drug=="Carbapenem", "Sub1", if_else(drug== "Metronidazole", "Sub1", if_else(drug=="Aerobic and Anaerobic Antibiotics", "Main", "Sub1")))))))))))))
fp.inst <- fp.inst %>% mutate(index = index +1)

inst.h <- ggplot(data = fp.inst, aes(y=index, x=effect, xmin=lower, xmax = upper)) +
  geom_point(color = a, size = 3.5) +
  geom_errorbarh(height = 0.2, color = a, size = 1) +
  scale_y_continuous(breaks = 1:13, labels = fp.inst$drug, expand=c(0, 0.5)) +
  scale_x_continuous(limits = c(-0.25, 0.3), breaks = c(-0.25, -0.15, 0, 0.15, 0.3), labels = c("-0.25", "-0.15", "0.0", "0.15", "0.3")) +
  labs(title = "Fluctuation of ARGs", x = "Change in Jaccard Distance") +
  geom_vline(xintercept = 0, color = "black", linetype = "longdash", alpha = 0.5) +
  theme(panel.background = element_blank(),
        plot.title = element_text(size = 18, hjust = 0.5),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_blank(),
        axis.text.x = element_text(color = "black", size=12),
        axis.text.y = element_blank(),
        axis.line.x = element_line(linetype = "solid"),
        panel.grid = element_blank(),
        axis.ticks.y = element_blank())
inst.h

fp.rpkm <- rbind(gamma.an2, gamma.piptazo1, gamma.carbapenem1, gamma.flagyl1, gamma.clinda1, gamma.cefepime1, gamma.vanc1, gamma.quinolone1, gamma.ag1, gamma.macrolide1, gamma.tmp_smx1) %>%
  mutate(index = if_else(drug == "Aerobic Antibiotics", 12, if_else(drug == 'Cefepime', 11, if_else(drug== "Vancomycin", 10, if_else(drug == "Fluoroquinolone", 9, if_else(drug== "Aminoglycoside", 8, if_else(drug=="Macrolide", 7, if_else(drug=="TMP-SMX", 6, if_else(drug=="Anaerobic Antibiotics", 5, if_else(drug=="Pip-Tazo", 4, if_else(drug=="Carbapenem", 3, if_else(drug== "Metronidazole", 2, if_else(drug=="Clindamycin", 1, 0))))))))))))) %>%
  arrange(index) %>% mutate(color = if_else(drug == "Aerobic Antibiotics", "Main", if_else(drug == 'Cefepime', "Sub", if_else(drug== "Vancomycin", "Sub", if_else(drug == "Fluoroquinolone", "Sub", if_else(drug== "Aminoglycoside", "Sub", if_else(drug=="Macrolide", "Sub", if_else(drug=="TMP-SMX", "Sub", if_else(drug=="Anaerobic Antibiotics", "Main", if_else(drug=="Pip-Tazo", "Sub1", if_else(drug=="Carbapenem", "Sub1", if_else(drug== "Metronidazole", "Sub1", if_else(drug=="Aerobic and Anaerobic Antibiotics", "Main", "Sub1")))))))))))))
fp.rpkm <- fp.rpkm %>% mutate(index = index +1)

rpkm.h <- ggplot(data = fp.rpkm, aes(y=index, x=effect, xmin=lower, xmax = upper)) +
  geom_point(color = a, size = 3.5) +
  geom_errorbarh(height = 0.2, color = a, size = 1) +
  scale_y_continuous(breaks = 1:13, labels = fp.rpkm$drug, expand=c(0, 0.5)) +
  scale_x_continuous(limits = c(0.4, 2.8), breaks = c(0.5, 1.0, 1.5, 2.0, 2.5, 2.8), labels = c("0.5", "1.0", "1.5", "2.0", "2.5", "2.8")) +
  labs(title = "Abundance of ARGs", x = expression(beta)) +
  geom_vline(xintercept = 1, color = "black", linetype = "longdash", alpha = 0.5) +
  theme(panel.background = element_blank(),
        plot.title = element_text(size = 18, hjust = 0.5),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_blank(),
        axis.text.x = element_text(color = "black", size=12),
        axis.text.y = element_blank(),
        axis.line.x = element_line(linetype = "solid"),
        panel.grid = element_blank(),
        axis.ticks.y = element_blank())
rpkm.h

ggarrange(num.h, new.h, rpkm.h, inst.h, nrow =1, widths = c(1.7,1,1,1), align = "h")

#Individual antibiotics on individual ARG class abundance and numbers (Supplementary Tables 4 and 5)
summary(test9$beta.abund)
cols <- seq.int(29, 41)
test9[cols] <- test9[cols] + 1
summary(test9$beta.abund)
summary(test9$num.betalactam.ARG)
library(broom.mixed)  

#cefepime exposure
rpkm.cefepime.beta <- glmer(beta.abund ~ cefepime.exp + (1 |study_id), data  = test9, family=Gamma(link="log"))
summary(rpkm.cefepime.beta)
rpkm.cefepime.beta2 <- summary(rpkm.cefepime.beta)
rpkm.cefepime.beta2 <- as.data.frame(rpkm.cefepime.beta2$coefficients)
rpkm.cef.com <- confint(rpkm.cefepime.beta, method="Wald")
rpkm.cef.com <- rpkm.cef.com[-c(1, 2),]
rpkm.cefepime1 <- cbind(rpkm.cefepime.beta2, rpkm.cef.com)
rpkm.cefepime1 <- rpkm.cefepime1[-1,]
rpkm.cefepime.beta <- rpkm.cefepime1 %>% mutate(effect = round(exp(Estimate), digits =2), lower=round(exp(`2.5 %`), digits =2), upper=round(exp(`97.5 %`), digits =2)) %>%
  mutate(exposure="Cefepime", ARG_class = "Beta-Lactam", model= "Abundance", p.value = `Pr(>|z|)`) %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)

num.cefepime.betalactam <- glmmadmb(num.betalactam.ARG ~ cefepime.exp + (1|study_id), data = test9, family = "nbinom1")
summary(num.cefepime.betalactam)
num.cefepime.beta <- tidy(num.cefepime.betalactam,
     effects = "fixed",
     component = "cond",
     scales=NULL,
     ran_prefix = NULL,
     conf.int = TRUE,
     conf.level = 0.95,
     conf.method = "Wald")
num.cefepime.beta <- num.cefepime.beta %>%
  mutate(exposure = "Cefepime", ARG_class = "Beta-Lactam", model = "Number", effect=round(exp(estimate),2),
         lower=round(exp(conf.low), 2), upper = round(exp(conf.high), 2)) %>%
  filter(term != "(Intercept)") %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)

rpkm.cefepime.vanc <- glmer(vanc.abund ~ cefepime.exp + (1 |study_id), data  = test9, family=Gamma(link="log"))
summary(rpkm.cefepime.vanc)
rpkm.cefepime.vanc2 <- summary(rpkm.cefepime.vanc)
rpkm.cefepime.vanc2 <- as.data.frame(rpkm.cefepime.vanc2$coefficients)
rpkm.cef.com <- confint(rpkm.cefepime.vanc, method="Wald")
rpkm.cef.com <- rpkm.cef.com[-c(1, 2),]
rpkm.cefepime1 <- cbind(rpkm.cefepime.vanc2, rpkm.cef.com)
rpkm.cefepime1 <- rpkm.cefepime1[-1,]
rpkm.cefepime.vanc <- rpkm.cefepime1 %>% mutate(effect = round(exp(Estimate), digits =2), lower=round(exp(`2.5 %`), digits =2), upper=round(exp(`97.5 %`), digits =2)) %>%
  mutate(exposure="Cefepime", ARG_class = "Glycopeptide", model= "Abundance", p.value = `Pr(>|z|)`) %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)

num.cefepime.vanc <- glmmadmb(num.vanc.ARG ~ cefepime.exp + (1|study_id), data = test9, family = "nbinom1")
summary(num.cefepime.vanc)
num.cefepime.vanc <- tidy(num.cefepime.vanc,
                          effects = "fixed",
                          component = "cond",
                          scales=NULL,
                          ran_prefix = NULL,
                          conf.int = TRUE,
                          conf.level = 0.95,
                          conf.method = "Wald")
num.cefepime.vanc <- num.cefepime.vanc %>%
  mutate(exposure = "Cefepime", ARG_class = "Glycopeptide", model = "Number", effect=round(exp(estimate),2),
         lower=round(exp(conf.low), 2), upper = round(exp(conf.high), 2)) %>%
  filter(term != "(Intercept)") %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)

rpkm.cefepime.carb <- glmer(carb.abund ~ cefepime.exp + (1 |study_id), data  = test9, family=Gamma(link="log"))
summary(rpkm.cefepime.carb)
rpkm.cefepime.carb2 <- summary(rpkm.cefepime.carb)
rpkm.cefepime.carb2 <- as.data.frame(rpkm.cefepime.carb2$coefficients)
rpkm.cef.com <- confint(rpkm.cefepime.carb, method="Wald")
rpkm.cef.com <- rpkm.cef.com[-c(1, 2),]
rpkm.cefepime1 <- cbind(rpkm.cefepime.carb2, rpkm.cef.com)
rpkm.cefepime1 <- rpkm.cefepime1[-1,]
rpkm.cefepime.carb <- rpkm.cefepime1 %>% mutate(effect = round(exp(Estimate), digits =2), lower=round(exp(`2.5 %`), digits =2), upper=round(exp(`97.5 %`), digits =2)) %>%
  mutate(exposure="Cefepime", ARG_class = "Carbapenem", model= "Abundance", p.value = `Pr(>|z|)`) %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)

num.cefepime.carb <- glmmadmb(num.carbapenem.ARG ~ cefepime.exp + (1|study_id), data = test9, family = "nbinom1")
summary(num.cefepime.carb)
num.cefepime.carb <- tidy(num.cefepime.carb,
                          effects = "fixed",
                          component = "cond",
                          scales=NULL,
                          ran_prefix = NULL,
                          conf.int = TRUE,
                          conf.level = 0.95,
                          conf.method = "Wald")
num.cefepime.carb <- num.cefepime.carb %>%
  mutate(exposure = "Cefepime", ARG_class = "Carbapenem", model = "Number", effect=round(exp(estimate),2),
         lower=round(exp(conf.low), 2), upper = round(exp(conf.high), 2)) %>%
  filter(term != "(Intercept)") %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)

rpkm.cefepime.fluoro <- glmer(fluoro.abund ~ cefepime.exp + (1 |study_id), data  = test9, family=Gamma(link="log"))
summary(rpkm.cefepime.fluoro)
rpkm.cefepime.fluoro2 <- summary(rpkm.cefepime.fluoro)
rpkm.cefepime.fluoro2 <- as.data.frame(rpkm.cefepime.fluoro2$coefficients)
rpkm.cef.com <- confint(rpkm.cefepime.fluoro, method="Wald")
rpkm.cef.com <- rpkm.cef.com[-c(1, 2),]
rpkm.cefepime1 <- cbind(rpkm.cefepime.fluoro2, rpkm.cef.com)
rpkm.cefepime1 <- rpkm.cefepime1[-1,]
rpkm.cefepime.fluoro <- rpkm.cefepime1 %>% mutate(effect = round(exp(Estimate), digits =2), lower=round(exp(`2.5 %`), digits =2), upper=round(exp(`97.5 %`), digits =2)) %>%
  mutate(exposure="Cefepime", ARG_class = "Fluoroquinolone", model= "Abundance", p.value = `Pr(>|z|)`) %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)

num.cefepime.fluoro <- glmmadmb(num.quinolone.ARG ~ cefepime.exp + (1|study_id), data = test9, family = "nbinom1")
summary(num.cefepime.fluoro)
num.cefepime.fluoro <- tidy(num.cefepime.fluoro,
                          effects = "fixed",
                          component = "cond",
                          scales=NULL,
                          ran_prefix = NULL,
                          conf.int = TRUE,
                          conf.level = 0.95,
                          conf.method = "Wald")
num.cefepime.fluoro <- num.cefepime.fluoro %>%
  mutate(exposure = "Cefepime", ARG_class = "Fluoroquinolone", model = "Number", effect=round(exp(estimate),2),
         lower=round(exp(conf.low), 2), upper = round(exp(conf.high), 2)) %>%
  filter(term != "(Intercept)") %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)

rpkm.cefepime.amino <- glmer(amino.abund ~ cefepime.exp + (1 |study_id), data  = test9, family=Gamma(link="log"))
summary(rpkm.cefepime.amino)
rpkm.cefepime.amino2 <- summary(rpkm.cefepime.amino)
rpkm.cefepime.amino2 <- as.data.frame(rpkm.cefepime.amino2$coefficients)
rpkm.cef.com <- confint(rpkm.cefepime.amino, method="Wald")
rpkm.cef.com <- rpkm.cef.com[-c(1, 2),]
rpkm.cefepime1 <- cbind(rpkm.cefepime.amino2, rpkm.cef.com)
rpkm.cefepime1 <- rpkm.cefepime1[-1,]
rpkm.cefepime.amino <- rpkm.cefepime1 %>% mutate(effect = round(exp(Estimate), digits =2), lower=round(exp(`2.5 %`), digits =2), upper=round(exp(`97.5 %`), digits =2)) %>%
  mutate(exposure="Cefepime", ARG_class = "Aminoglycoside", model= "Abundance", p.value = `Pr(>|z|)`) %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)

num.cefepime.amino <- glmmadmb(num.amino.ARG ~ cefepime.exp + (1|study_id), data = test9, family = "nbinom1")
summary(num.cefepime.amino)
num.cefepime.amino <- tidy(num.cefepime.amino,
                          effects = "fixed",
                          component = "cond",
                          scales=NULL,
                          ran_prefix = NULL,
                          conf.int = TRUE,
                          conf.level = 0.95,
                          conf.method = "Wald")
num.cefepime.amino <- num.cefepime.amino %>%
  mutate(exposure = "Cefepime", ARG_class = "Aminoglycoside", model = "Number", effect=round(exp(estimate),2),
         lower=round(exp(conf.low), 2), upper = round(exp(conf.high), 2)) %>%
  filter(term != "(Intercept)") %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)

rpkm.cefepime.clinda <- glmer(clinda.abund ~ cefepime.exp + (1 |study_id), data  = test9, family=Gamma(link="log"))
summary(rpkm.cefepime.clinda)
rpkm.cefepime.clinda2 <- summary(rpkm.cefepime.clinda)
rpkm.cefepime.clinda2 <- as.data.frame(rpkm.cefepime.clinda2$coefficients)
rpkm.cef.com <- confint(rpkm.cefepime.clinda, method="Wald")
rpkm.cef.com <- rpkm.cef.com[-c(1, 2),]
rpkm.cefepime1 <- cbind(rpkm.cefepime.clinda2, rpkm.cef.com)
rpkm.cefepime1 <- rpkm.cefepime1[-1,]
rpkm.cefepime.clinda <- rpkm.cefepime1 %>% mutate(effect = round(exp(Estimate), digits =2), lower=round(exp(`2.5 %`), digits =2), upper=round(exp(`97.5 %`), digits =2)) %>%
  mutate(exposure="Cefepime", ARG_class = "Lincosamide", model= "Abundance", p.value = `Pr(>|z|)`) %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)

num.cefepime.clinda <- glmmadmb(num.clinda.ARG ~ cefepime.exp + (1|study_id), data = test9, family = "nbinom1")
summary(num.cefepime.clinda)
num.cefepime.clinda <- tidy(num.cefepime.clinda,
                           effects = "fixed",
                           component = "cond",
                           scales=NULL,
                           ran_prefix = NULL,
                           conf.int = TRUE,
                           conf.level = 0.95,
                           conf.method = "Wald")
num.cefepime.clinda <- num.cefepime.clinda %>%
  mutate(exposure = "Cefepime", ARG_class = "Lincosamide", model = "Number", effect=round(exp(estimate),2),
         lower=round(exp(conf.low), 2), upper = round(exp(conf.high), 2)) %>%
  filter(term != "(Intercept)") %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)

rpkm.cefepime.macrolide <- glmer(macrolide.abund ~ cefepime.exp + (1 |study_id), data  = test9, family=Gamma(link="log"))
summary(rpkm.cefepime.macrolide)
rpkm.cefepime.macrolide2 <- summary(rpkm.cefepime.macrolide)
rpkm.cefepime.macrolide2 <- as.data.frame(rpkm.cefepime.macrolide2$coefficients)
rpkm.cef.com <- confint(rpkm.cefepime.macrolide, method="Wald")
rpkm.cef.com <- rpkm.cef.com[-c(1, 2),]
rpkm.cefepime1 <- cbind(rpkm.cefepime.macrolide2, rpkm.cef.com)
rpkm.cefepime1 <- rpkm.cefepime1[-1,]
rpkm.cefepime.macrolide <- rpkm.cefepime1 %>% mutate(effect = round(exp(Estimate), digits =2), lower=round(exp(`2.5 %`), digits =2), upper=round(exp(`97.5 %`), digits =2)) %>%
  mutate(exposure="Cefepime", ARG_class = "Macrolide", model= "Abundance", p.value = `Pr(>|z|)`) %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)

num.cefepime.macrolide <- glmmadmb(num.macrolide.ARG ~ cefepime.exp + (1|study_id), data = test9, family = "nbinom1")
summary(num.cefepime.macrolide)
num.cefepime.macrolide <- tidy(num.cefepime.macrolide,
                           effects = "fixed",
                           component = "cond",
                           scales=NULL,
                           ran_prefix = NULL,
                           conf.int = TRUE,
                           conf.level = 0.95,
                           conf.method = "Wald")
num.cefepime.macrolide <- num.cefepime.macrolide %>%
  mutate(exposure = "Cefepime", ARG_class = "Macrolide", model = "Number", effect=round(exp(estimate),2),
         lower=round(exp(conf.low), 2), upper = round(exp(conf.high), 2)) %>%
  filter(term != "(Intercept)") %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)

rpkm.cefepime.trimeth <- glmer(trimeth.abund ~ cefepime.exp + (1 |study_id), data  = test9, family=Gamma(link="log"))
summary(rpkm.cefepime.trimeth)
rpkm.cefepime.trimeth2 <- summary(rpkm.cefepime.trimeth)
rpkm.cefepime.trimeth2 <- as.data.frame(rpkm.cefepime.trimeth2$coefficients)
rpkm.cef.com <- confint(rpkm.cefepime.trimeth, method="Wald")
rpkm.cef.com <- rpkm.cef.com[-c(1, 2),]
rpkm.cefepime1 <- cbind(rpkm.cefepime.trimeth2, rpkm.cef.com)
rpkm.cefepime1 <- rpkm.cefepime1[-1,]
rpkm.cefepime.trimeth <- rpkm.cefepime1 %>% mutate(effect = round(exp(Estimate), digits =2), lower=round(exp(`2.5 %`), digits =2), upper=round(exp(`97.5 %`), digits =2)) %>%
  mutate(exposure="Cefepime", ARG_class = "Diaminopyrimidine", model= "Abundance", p.value = `Pr(>|z|)`) %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)

num.cefepime.trimeth <- glmmadmb(num.trimeth.ARG ~ cefepime.exp + (1|study_id), data = test9, family = "nbinom1")
summary(num.cefepime.trimeth)
num.cefepime.trimeth <- tidy(num.cefepime.trimeth,
                           effects = "fixed",
                           component = "cond",
                           scales=NULL,
                           ran_prefix = NULL,
                           conf.int = TRUE,
                           conf.level = 0.95,
                           conf.method = "Wald")
num.cefepime.trimeth <- num.cefepime.trimeth %>%
  mutate(exposure = "Cefepime", ARG_class = "Diaminopyrimidine", model = "Number", effect=round(exp(estimate),2),
         lower=round(exp(conf.low), 2), upper = round(exp(conf.high), 2)) %>%
  filter(term != "(Intercept)") %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)

rpkm.cefepime.tetra <- glmer(tetra.abund ~ cefepime.exp + (1 |study_id), data  = test9, family=Gamma(link="log"))
summary(rpkm.cefepime.tetra)
rpkm.cefepime.tetra2 <- summary(rpkm.cefepime.tetra)
rpkm.cefepime.tetra2 <- as.data.frame(rpkm.cefepime.tetra2$coefficients)
rpkm.cef.com <- confint(rpkm.cefepime.tetra, method="Wald")
rpkm.cef.com <- rpkm.cef.com[-c(1, 2),]
rpkm.cefepime1 <- cbind(rpkm.cefepime.tetra2, rpkm.cef.com)
rpkm.cefepime1 <- rpkm.cefepime1[-1,]
rpkm.cefepime.tetra <- rpkm.cefepime1 %>% mutate(effect = round(exp(Estimate), digits =2), lower=round(exp(`2.5 %`), digits =2), upper=round(exp(`97.5 %`), digits =2)) %>%
  mutate(exposure="Cefepime", ARG_class = "Tetracycline", model= "Abundance", p.value = `Pr(>|z|)`) %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)

num.cefepime.tetra <- glmmadmb(num.tetra.ARG ~ cefepime.exp + (1|study_id), data = test9, family = "nbinom1")
summary(num.cefepime.tetra)
num.cefepime.tetra <- tidy(num.cefepime.tetra,
                           effects = "fixed",
                           component = "cond",
                           scales=NULL,
                           ran_prefix = NULL,
                           conf.int = TRUE,
                           conf.level = 0.95,
                           conf.method = "Wald")
num.cefepime.tetra <- num.cefepime.tetra %>%
  mutate(exposure = "Cefepime", ARG_class = "Tetracycline", model = "Number", effect=round(exp(estimate),2),
         lower=round(exp(conf.low), 2), upper = round(exp(conf.high), 2)) %>%
  filter(term != "(Intercept)") %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)

rpkm.cefepime.flagyl <- glmer(flagyl.abund ~ cefepime.exp + (1 |study_id), data  = test9, family=Gamma(link="log"))
summary(rpkm.cefepime.flagyl)
rpkm.cefepime.flagyl2 <- summary(rpkm.cefepime.flagyl)
rpkm.cefepime.flagyl2 <- as.data.frame(rpkm.cefepime.flagyl2$coefficients)
rpkm.cef.com <- confint(rpkm.cefepime.flagyl, method="Wald")
rpkm.cef.com <- rpkm.cef.com[-c(1, 2),]
rpkm.cefepime1 <- cbind(rpkm.cefepime.flagyl2, rpkm.cef.com)
rpkm.cefepime1 <- rpkm.cefepime1[-1,]
rpkm.cefepime.flagyl <- rpkm.cefepime1 %>% mutate(effect = round(exp(Estimate), digits =2), lower=round(exp(`2.5 %`), digits =2), upper=round(exp(`97.5 %`), digits =2)) %>%
  mutate(exposure="Cefepime", ARG_class = "Nitroimidazole", model= "Abundance", p.value = `Pr(>|z|)`) %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)

num.cefepime.flagyl <- glmmadmb(num.flagyl.ARG ~ cefepime.exp + (1|study_id), data = test9, family = "nbinom1")
summary(num.cefepime.flagyl)
num.cefepime.flagyl <- tidy(num.cefepime.flagyl,
                           effects = "fixed",
                           component = "cond",
                           scales=NULL,
                           ran_prefix = NULL,
                           conf.int = TRUE,
                           conf.level = 0.95,
                           conf.method = "Wald")
num.cefepime.flagyl <- num.cefepime.flagyl %>%
  mutate(exposure = "Cefepime", ARG_class = "Nitroimidazole", model = "Number", effect=round(exp(estimate),2),
         lower=round(exp(conf.low), 2), upper = round(exp(conf.high), 2)) %>%
  filter(term != "(Intercept)") %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)

cefepime <- rbind(rpkm.cefepime.beta,
                  num.cefepime.beta,
                  rpkm.cefepime.carb,
                  num.cefepime.carb,
                  rpkm.cefepime.vanc,
                  num.cefepime.vanc,
                  rpkm.cefepime.flagyl,
                  num.cefepime.flagyl,
                  rpkm.cefepime.fluoro,
                  num.cefepime.fluoro,
                  rpkm.cefepime.trimeth,
                  num.cefepime.trimeth,
                  rpkm.cefepime.macrolide,
                  num.cefepime.macrolide,
                  rpkm.cefepime.tetra,
                  num.cefepime.tetra,
                  rpkm.cefepime.amino,
                  num.cefepime.amino,
                  rpkm.cefepime.clinda,
                  num.cefepime.clinda)
rownames(cefepime) <- NULL

cefepime <- cefepime %>%
  mutate(p.adj = p.adjust(p.value, method = "BH"))
cefepime$p.adj <- signif(cefepime$p.adj, 2)
cefepime$p.adj <- format(cefepime$p.adj, scientific = FALSE)

#piptazo exposure
rpkm.piptazo.beta <- glmer(beta.abund ~ piptazo.exp + (1 |study_id), data  = test9, family=Gamma(link="log"))
summary(rpkm.piptazo.beta)
rpkm.piptazo.beta2 <- summary(rpkm.piptazo.beta)
rpkm.piptazo.beta2 <- as.data.frame(rpkm.piptazo.beta2$coefficients)
rpkm.cef.com <- confint(rpkm.piptazo.beta, method="Wald")
rpkm.cef.com <- rpkm.cef.com[-c(1, 2),]
rpkm.piptazo1 <- cbind(rpkm.piptazo.beta2, rpkm.cef.com)
rpkm.piptazo1 <- rpkm.piptazo1[-1,]
rpkm.piptazo.beta <- rpkm.piptazo1 %>% mutate(effect = round(exp(Estimate), digits =2), lower=round(exp(`2.5 %`), digits =2), upper=round(exp(`97.5 %`), digits =2)) %>%
  mutate(exposure="piptazo", ARG_class = "Beta-Lactam", model= "Abundance", p.value = `Pr(>|z|)`) %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)

num.piptazo.betalactam <- glmmadmb(num.betalactam.ARG ~ piptazo.exp + (1|study_id), data = test9, family = "nbinom1")
summary(num.piptazo.betalactam)
num.piptazo.beta <- tidy(num.piptazo.betalactam,
                          effects = "fixed",
                          component = "cond",
                          scales=NULL,
                          ran_prefix = NULL,
                          conf.int = TRUE,
                          conf.level = 0.95,
                          conf.method = "Wald")
num.piptazo.beta <- num.piptazo.beta %>%
  mutate(exposure = "piptazo", ARG_class = "Beta-Lactam", model = "Number", effect=round(exp(estimate),2),
         lower=round(exp(conf.low), 2), upper = round(exp(conf.high), 2)) %>%
  filter(term != "(Intercept)") %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)

rpkm.piptazo.vanc <- glmer(vanc.abund ~ piptazo.exp + (1 |study_id), data  = test9, family=Gamma(link="log"))
summary(rpkm.piptazo.vanc)
rpkm.piptazo.vanc2 <- summary(rpkm.piptazo.vanc)
rpkm.piptazo.vanc2 <- as.data.frame(rpkm.piptazo.vanc2$coefficients)
rpkm.cef.com <- confint(rpkm.piptazo.vanc, method="Wald")
rpkm.cef.com <- rpkm.cef.com[-c(1, 2),]
rpkm.piptazo1 <- cbind(rpkm.piptazo.vanc2, rpkm.cef.com)
rpkm.piptazo1 <- rpkm.piptazo1[-1,]
rpkm.piptazo.vanc <- rpkm.piptazo1 %>% mutate(effect = round(exp(Estimate), digits =2), lower=round(exp(`2.5 %`), digits =2), upper=round(exp(`97.5 %`), digits =2)) %>%
  mutate(exposure="piptazo", ARG_class = "Glycopeptide", model= "Abundance", p.value = `Pr(>|z|)`) %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)

num.piptazo.vanc <- glmmadmb(num.vanc.ARG ~ piptazo.exp + (1|study_id), data = test9, family = "nbinom1")
summary(num.piptazo.vanc)
num.piptazo.vanc <- tidy(num.piptazo.vanc,
                          effects = "fixed",
                          component = "cond",
                          scales=NULL,
                          ran_prefix = NULL,
                          conf.int = TRUE,
                          conf.level = 0.95,
                          conf.method = "Wald")
num.piptazo.vanc <- num.piptazo.vanc %>%
  mutate(exposure = "piptazo", ARG_class = "Glycopeptide", model = "Number", effect=round(exp(estimate),2),
         lower=round(exp(conf.low), 2), upper = round(exp(conf.high), 2)) %>%
  filter(term != "(Intercept)") %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)


rpkm.piptazo.carb <- glmer(carb.abund ~ piptazo.exp + (1 |study_id), data  = test9, family=Gamma(link="log"))
summary(rpkm.piptazo.carb)
rpkm.piptazo.carb2 <- summary(rpkm.piptazo.carb)
rpkm.piptazo.carb2 <- as.data.frame(rpkm.piptazo.carb2$coefficients)
rpkm.cef.com <- confint(rpkm.piptazo.carb, method="Wald")
rpkm.cef.com <- rpkm.cef.com[-c(1, 2),]
rpkm.piptazo1 <- cbind(rpkm.piptazo.carb2, rpkm.cef.com)
rpkm.piptazo1 <- rpkm.piptazo1[-1,]
rpkm.piptazo.carb <- rpkm.piptazo1 %>% mutate(effect = round(exp(Estimate), digits =2), lower=round(exp(`2.5 %`), digits =2), upper=round(exp(`97.5 %`), digits =2)) %>%
  mutate(exposure="piptazo", ARG_class = "Carbapenem", model= "Abundance", p.value = `Pr(>|z|)`) %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)

num.piptazo.carb <- glmmadmb(num.carbapenem.ARG ~ piptazo.exp + (1|study_id), data = test9, family = "nbinom1")
summary(num.piptazo.carb)
num.piptazo.carb <- tidy(num.piptazo.carb,
                          effects = "fixed",
                          component = "cond",
                          scales=NULL,
                          ran_prefix = NULL,
                          conf.int = TRUE,
                          conf.level = 0.95,
                          conf.method = "Wald")
num.piptazo.carb <- num.piptazo.carb %>%
  mutate(exposure = "piptazo", ARG_class = "Carbapenem", model = "Number", effect=round(exp(estimate),2),
         lower=round(exp(conf.low), 2), upper = round(exp(conf.high), 2)) %>%
  filter(term != "(Intercept)") %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)

rpkm.piptazo.fluoro <- glmer(fluoro.abund ~ piptazo.exp + (1 |study_id), data  = test9, family=Gamma(link="log"))
summary(rpkm.piptazo.fluoro)
rpkm.piptazo.fluoro2 <- summary(rpkm.piptazo.fluoro)
rpkm.piptazo.fluoro2 <- as.data.frame(rpkm.piptazo.fluoro2$coefficients)
rpkm.cef.com <- confint(rpkm.piptazo.fluoro, method="Wald")
rpkm.cef.com <- rpkm.cef.com[-c(1, 2),]
rpkm.piptazo1 <- cbind(rpkm.piptazo.fluoro2, rpkm.cef.com)
rpkm.piptazo1 <- rpkm.piptazo1[-1,]
rpkm.piptazo.fluoro <- rpkm.piptazo1 %>% mutate(effect = round(exp(Estimate), digits =2), lower=round(exp(`2.5 %`), digits =2), upper=round(exp(`97.5 %`), digits =2)) %>%
  mutate(exposure="piptazo", ARG_class = "Fluoroquinolone", model= "Abundance", p.value = `Pr(>|z|)`) %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)

num.piptazo.fluoro <- glmmadmb(num.quinolone.ARG ~ piptazo.exp + (1|study_id), data = test9, family = "nbinom1")
summary(num.piptazo.fluoro)
num.piptazo.fluoro <- tidy(num.piptazo.fluoro,
                            effects = "fixed",
                            component = "cond",
                            scales=NULL,
                            ran_prefix = NULL,
                            conf.int = TRUE,
                            conf.level = 0.95,
                            conf.method = "Wald")
num.piptazo.fluoro <- num.piptazo.fluoro %>%
  mutate(exposure = "piptazo", ARG_class = "Fluoroquinolone", model = "Number", effect=round(exp(estimate),2),
         lower=round(exp(conf.low), 2), upper = round(exp(conf.high), 2)) %>%
  filter(term != "(Intercept)") %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)


rpkm.piptazo.amino <- glmer(amino.abund ~ piptazo.exp + (1 |study_id), data  = test9, family=Gamma(link="log"))
summary(rpkm.piptazo.amino)
rpkm.piptazo.amino2 <- summary(rpkm.piptazo.amino)
rpkm.piptazo.amino2 <- as.data.frame(rpkm.piptazo.amino2$coefficients)
rpkm.cef.com <- confint(rpkm.piptazo.amino, method="Wald")
rpkm.cef.com <- rpkm.cef.com[-c(1, 2),]
rpkm.piptazo1 <- cbind(rpkm.piptazo.amino2, rpkm.cef.com)
rpkm.piptazo1 <- rpkm.piptazo1[-1,]
rpkm.piptazo.amino <- rpkm.piptazo1 %>% mutate(effect = round(exp(Estimate), digits =2), lower=round(exp(`2.5 %`), digits =2), upper=round(exp(`97.5 %`), digits =2)) %>%
  mutate(exposure="piptazo", ARG_class = "Aminoglycoside", model= "Abundance", p.value = `Pr(>|z|)`) %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)

num.piptazo.amino <- glmmadmb(num.amino.ARG ~ piptazo.exp + (1|study_id), data = test9, family = "nbinom1")
summary(num.piptazo.amino)
num.piptazo.amino <- tidy(num.piptazo.amino,
                           effects = "fixed",
                           component = "cond",
                           scales=NULL,
                           ran_prefix = NULL,
                           conf.int = TRUE,
                           conf.level = 0.95,
                           conf.method = "Wald")
num.piptazo.amino <- num.piptazo.amino %>%
  mutate(exposure = "piptazo", ARG_class = "Aminoglycoside", model = "Number", effect=round(exp(estimate),2),
         lower=round(exp(conf.low), 2), upper = round(exp(conf.high), 2)) %>%
  filter(term != "(Intercept)") %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)


rpkm.piptazo.clinda <- glmer(clinda.abund ~ piptazo.exp + (1 |study_id), data  = test9, family=Gamma(link="log"))
summary(rpkm.piptazo.clinda)
rpkm.piptazo.clinda2 <- summary(rpkm.piptazo.clinda)
rpkm.piptazo.clinda2 <- as.data.frame(rpkm.piptazo.clinda2$coefficients)
rpkm.cef.com <- confint(rpkm.piptazo.clinda, method="Wald")
rpkm.cef.com <- rpkm.cef.com[-c(1, 2),]
rpkm.piptazo1 <- cbind(rpkm.piptazo.clinda2, rpkm.cef.com)
rpkm.piptazo1 <- rpkm.piptazo1[-1,]
rpkm.piptazo.clinda <- rpkm.piptazo1 %>% mutate(effect = round(exp(Estimate), digits =2), lower=round(exp(`2.5 %`), digits =2), upper=round(exp(`97.5 %`), digits =2)) %>%
  mutate(exposure="piptazo", ARG_class = "Lincosamide", model= "Abundance", p.value = `Pr(>|z|)`) %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)

num.piptazo.clinda <- glmmadmb(num.clinda.ARG ~ piptazo.exp + (1|study_id), data = test9, family = "nbinom1")
summary(num.piptazo.clinda)
num.piptazo.clinda <- tidy(num.piptazo.clinda,
                            effects = "fixed",
                            component = "cond",
                            scales=NULL,
                            ran_prefix = NULL,
                            conf.int = TRUE,
                            conf.level = 0.95,
                            conf.method = "Wald")
num.piptazo.clinda <- num.piptazo.clinda %>%
  mutate(exposure = "piptazo", ARG_class = "Lincosamide", model = "Number", effect=round(exp(estimate),2),
         lower=round(exp(conf.low), 2), upper = round(exp(conf.high), 2)) %>%
  filter(term != "(Intercept)") %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)



rpkm.piptazo.macrolide <- glmer(macrolide.abund ~ piptazo.exp + (1 |study_id), data  = test9, family=Gamma(link="log"))
summary(rpkm.piptazo.macrolide)
rpkm.piptazo.macrolide2 <- summary(rpkm.piptazo.macrolide)
rpkm.piptazo.macrolide2 <- as.data.frame(rpkm.piptazo.macrolide2$coefficients)
rpkm.cef.com <- confint(rpkm.piptazo.macrolide, method="Wald")
rpkm.cef.com <- rpkm.cef.com[-c(1, 2),]
rpkm.piptazo1 <- cbind(rpkm.piptazo.macrolide2, rpkm.cef.com)
rpkm.piptazo1 <- rpkm.piptazo1[-1,]
rpkm.piptazo.macrolide <- rpkm.piptazo1 %>% mutate(effect = round(exp(Estimate), digits =2), lower=round(exp(`2.5 %`), digits =2), upper=round(exp(`97.5 %`), digits =2)) %>%
  mutate(exposure="piptazo", ARG_class = "Macrolide", model= "Abundance", p.value = `Pr(>|z|)`) %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)

num.piptazo.macrolide <- glmmadmb(num.macrolide.ARG ~ piptazo.exp + (1|study_id), data = test9, family = "nbinom1")
summary(num.piptazo.macrolide)
num.piptazo.macrolide <- tidy(num.piptazo.macrolide,
                               effects = "fixed",
                               component = "cond",
                               scales=NULL,
                               ran_prefix = NULL,
                               conf.int = TRUE,
                               conf.level = 0.95,
                               conf.method = "Wald")
num.piptazo.macrolide <- num.piptazo.macrolide %>%
  mutate(exposure = "piptazo", ARG_class = "Macrolide", model = "Number", effect=round(exp(estimate),2),
         lower=round(exp(conf.low), 2), upper = round(exp(conf.high), 2)) %>%
  filter(term != "(Intercept)") %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)


rpkm.piptazo.trimeth <- glmer(trimeth.abund ~ piptazo.exp + (1 |study_id), data  = test9, family=Gamma(link="log"))
summary(rpkm.piptazo.trimeth)
rpkm.piptazo.trimeth2 <- summary(rpkm.piptazo.trimeth)
rpkm.piptazo.trimeth2 <- as.data.frame(rpkm.piptazo.trimeth2$coefficients)
rpkm.cef.com <- confint(rpkm.piptazo.trimeth, method="Wald")
rpkm.cef.com <- rpkm.cef.com[-c(1, 2),]
rpkm.piptazo1 <- cbind(rpkm.piptazo.trimeth2, rpkm.cef.com)
rpkm.piptazo1 <- rpkm.piptazo1[-1,]
rpkm.piptazo.trimeth <- rpkm.piptazo1 %>% mutate(effect = round(exp(Estimate), digits =2), lower=round(exp(`2.5 %`), digits =2), upper=round(exp(`97.5 %`), digits =2)) %>%
  mutate(exposure="piptazo", ARG_class = "Diaminopyrimidine", model= "Abundance", p.value = `Pr(>|z|)`) %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)

num.piptazo.trimeth <- glmmadmb(num.trimeth.ARG ~ piptazo.exp + (1|study_id), data = test9, family = "nbinom1")
summary(num.piptazo.trimeth)
num.piptazo.trimeth <- tidy(num.piptazo.trimeth,
                             effects = "fixed",
                             component = "cond",
                             scales=NULL,
                             ran_prefix = NULL,
                             conf.int = TRUE,
                             conf.level = 0.95,
                             conf.method = "Wald")
num.piptazo.trimeth <- num.piptazo.trimeth %>%
  mutate(exposure = "piptazo", ARG_class = "Diaminopyrimidine", model = "Number", effect=round(exp(estimate),2),
         lower=round(exp(conf.low), 2), upper = round(exp(conf.high), 2)) %>%
  filter(term != "(Intercept)") %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)


rpkm.piptazo.tetra <- glmer(tetra.abund ~ piptazo.exp + (1 |study_id), data  = test9, family=Gamma(link="log"))
summary(rpkm.piptazo.tetra)
rpkm.piptazo.tetra2 <- summary(rpkm.piptazo.tetra)
rpkm.piptazo.tetra2 <- as.data.frame(rpkm.piptazo.tetra2$coefficients)
rpkm.cef.com <- confint(rpkm.piptazo.tetra, method="Wald")
rpkm.cef.com <- rpkm.cef.com[-c(1, 2),]
rpkm.piptazo1 <- cbind(rpkm.piptazo.tetra2, rpkm.cef.com)
rpkm.piptazo1 <- rpkm.piptazo1[-1,]
rpkm.piptazo.tetra <- rpkm.piptazo1 %>% mutate(effect = round(exp(Estimate), digits =2), lower=round(exp(`2.5 %`), digits =2), upper=round(exp(`97.5 %`), digits =2)) %>%
  mutate(exposure="piptazo", ARG_class = "Tetracycline", model= "Abundance", p.value = `Pr(>|z|)`) %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)

num.piptazo.tetra <- glmmadmb(num.tetra.ARG ~ piptazo.exp + (1|study_id), data = test9, family = "nbinom1")
summary(num.piptazo.tetra)
num.piptazo.tetra <- tidy(num.piptazo.tetra,
                           effects = "fixed",
                           component = "cond",
                           scales=NULL,
                           ran_prefix = NULL,
                           conf.int = TRUE,
                           conf.level = 0.95,
                           conf.method = "Wald")
num.piptazo.tetra <- num.piptazo.tetra %>%
  mutate(exposure = "piptazo", ARG_class = "Tetracycline", model = "Number", effect=round(exp(estimate),2),
         lower=round(exp(conf.low), 2), upper = round(exp(conf.high), 2)) %>%
  filter(term != "(Intercept)") %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)


rpkm.piptazo.flagyl <- glmer(flagyl.abund ~ piptazo.exp + (1 |study_id), data  = test9, family=Gamma(link="log"))
summary(rpkm.piptazo.flagyl)
rpkm.piptazo.flagyl2 <- summary(rpkm.piptazo.flagyl)
rpkm.piptazo.flagyl2 <- as.data.frame(rpkm.piptazo.flagyl2$coefficients)
rpkm.cef.com <- confint(rpkm.piptazo.flagyl, method="Wald")
rpkm.cef.com <- rpkm.cef.com[-c(1, 2),]
rpkm.piptazo1 <- cbind(rpkm.piptazo.flagyl2, rpkm.cef.com)
rpkm.piptazo1 <- rpkm.piptazo1[-1,]
rpkm.piptazo.flagyl <- rpkm.piptazo1 %>% mutate(effect = round(exp(Estimate), digits =2), lower=round(exp(`2.5 %`), digits =2), upper=round(exp(`97.5 %`), digits =2)) %>%
  mutate(exposure="piptazo", ARG_class = "Nitroimidazole", model= "Abundance", p.value = `Pr(>|z|)`) %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)

num.piptazo.flagyl <- glmmadmb(num.flagyl.ARG ~ piptazo.exp + (1|study_id), data = test9, family = "nbinom1")
summary(num.piptazo.flagyl)
num.piptazo.flagyl <- tidy(num.piptazo.flagyl,
                            effects = "fixed",
                            component = "cond",
                            scales=NULL,
                            ran_prefix = NULL,
                            conf.int = TRUE,
                            conf.level = 0.95,
                            conf.method = "Wald")
num.piptazo.flagyl <- num.piptazo.flagyl %>%
  mutate(exposure = "piptazo", ARG_class = "Nitroimidazole", model = "Number", effect=round(exp(estimate),2),
         lower=round(exp(conf.low), 2), upper = round(exp(conf.high), 2)) %>%
  filter(term != "(Intercept)") %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)

piptazo <- rbind(rpkm.piptazo.beta,
                  num.piptazo.beta,
                  rpkm.piptazo.carb,
                  num.piptazo.carb,
                  rpkm.piptazo.vanc,
                  num.piptazo.vanc,
                  rpkm.piptazo.flagyl,
                  num.piptazo.flagyl,
                  rpkm.piptazo.fluoro,
                  num.piptazo.fluoro,
                  rpkm.piptazo.trimeth,
                  num.piptazo.trimeth,
                  rpkm.piptazo.macrolide,
                  num.piptazo.macrolide,
                  rpkm.piptazo.tetra,
                  num.piptazo.tetra,
                  rpkm.piptazo.amino,
                  num.piptazo.amino,
                  rpkm.piptazo.clinda,
                  num.piptazo.clinda)
rownames(piptazo) <- NULL

piptazo <- piptazo %>%
  mutate(p.adj = p.adjust(p.value, method = "BH"))
piptazo$p.adj <- signif(piptazo$p.adj, 2)
piptazo$p.adj <- format(piptazo$p.adj, scientific = FALSE)

#carbapenem exposure
rpkm.carbapenem.beta <- glmer(beta.abund ~ carbapenem.exp + (1 |study_id), data  = test9, family=Gamma(link="log"))
summary(rpkm.carbapenem.beta)
rpkm.carbapenem.beta2 <- summary(rpkm.carbapenem.beta)
rpkm.carbapenem.beta2 <- as.data.frame(rpkm.carbapenem.beta2$coefficients)
rpkm.cef.com <- confint(rpkm.carbapenem.beta, method="Wald")
rpkm.cef.com <- rpkm.cef.com[-c(1, 2),]
rpkm.carbapenem1 <- cbind(rpkm.carbapenem.beta2, rpkm.cef.com)
rpkm.carbapenem1 <- rpkm.carbapenem1[-1,]
rpkm.carbapenem.beta <- rpkm.carbapenem1 %>% mutate(effect = round(exp(Estimate), digits =2), lower=round(exp(`2.5 %`), digits =2), upper=round(exp(`97.5 %`), digits =2)) %>%
  mutate(exposure="carbapenem", ARG_class = "Beta-Lactam", model= "Abundance", p.value = `Pr(>|z|)`) %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)

num.carbapenem.betalactam <- glmmadmb(num.betalactam.ARG ~ carbapenem.exp + (1|study_id), data = test9, family = "nbinom1")
summary(num.carbapenem.betalactam)
num.carbapenem.beta <- tidy(num.carbapenem.betalactam,
                         effects = "fixed",
                         component = "cond",
                         scales=NULL,
                         ran_prefix = NULL,
                         conf.int = TRUE,
                         conf.level = 0.95,
                         conf.method = "Wald")
num.carbapenem.beta <- num.carbapenem.beta %>%
  mutate(exposure = "carbapenem", ARG_class = "Beta-Lactam", model = "Number", effect=round(exp(estimate),2),
         lower=round(exp(conf.low), 2), upper = round(exp(conf.high), 2)) %>%
  filter(term != "(Intercept)") %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)

rpkm.carbapenem.vanc <- glmer(vanc.abund ~ carbapenem.exp + (1 |study_id), data  = test9, family=Gamma(link="log"))
summary(rpkm.carbapenem.vanc)
rpkm.carbapenem.vanc2 <- summary(rpkm.carbapenem.vanc)
rpkm.carbapenem.vanc2 <- as.data.frame(rpkm.carbapenem.vanc2$coefficients)
rpkm.cef.com <- confint(rpkm.carbapenem.vanc, method="Wald")
rpkm.cef.com <- rpkm.cef.com[-c(1, 2),]
rpkm.carbapenem1 <- cbind(rpkm.carbapenem.vanc2, rpkm.cef.com)
rpkm.carbapenem1 <- rpkm.carbapenem1[-1,]
rpkm.carbapenem.vanc <- rpkm.carbapenem1 %>% mutate(effect = round(exp(Estimate), digits =2), lower=round(exp(`2.5 %`), digits =2), upper=round(exp(`97.5 %`), digits =2)) %>%
  mutate(exposure="carbapenem", ARG_class = "Glycopeptide", model= "Abundance", p.value = `Pr(>|z|)`) %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)

num.carbapenem.vanc <- glmmadmb(num.vanc.ARG ~ carbapenem.exp + (1|study_id), data = test9, family = "nbinom1")
summary(num.carbapenem.vanc)
num.carbapenem.vanc <- tidy(num.carbapenem.vanc,
                         effects = "fixed",
                         component = "cond",
                         scales=NULL,
                         ran_prefix = NULL,
                         conf.int = TRUE,
                         conf.level = 0.95,
                         conf.method = "Wald")
num.carbapenem.vanc <- num.carbapenem.vanc %>%
  mutate(exposure = "carbapenem", ARG_class = "Glycopeptide", model = "Number", effect=round(exp(estimate),2),
         lower=round(exp(conf.low), 2), upper = round(exp(conf.high), 2)) %>%
  filter(term != "(Intercept)") %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)


rpkm.carbapenem.carb <- glmer(carb.abund ~ carbapenem.exp + (1 |study_id), data  = test9, family=Gamma(link="log"))
summary(rpkm.carbapenem.carb)
rpkm.carbapenem.carb2 <- summary(rpkm.carbapenem.carb)
rpkm.carbapenem.carb2 <- as.data.frame(rpkm.carbapenem.carb2$coefficients)
rpkm.cef.com <- confint(rpkm.carbapenem.carb, method="Wald")
rpkm.cef.com <- rpkm.cef.com[-c(1, 2),]
rpkm.carbapenem1 <- cbind(rpkm.carbapenem.carb2, rpkm.cef.com)
rpkm.carbapenem1 <- rpkm.carbapenem1[-1,]
rpkm.carbapenem.carb <- rpkm.carbapenem1 %>% mutate(effect = round(exp(Estimate), digits =2), lower=round(exp(`2.5 %`), digits =2), upper=round(exp(`97.5 %`), digits =2)) %>%
  mutate(exposure="carbapenem", ARG_class = "Carbapenem", model= "Abundance", p.value = `Pr(>|z|)`) %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)

num.carbapenem.carb <- glmmadmb(num.carbapenem.ARG ~ carbapenem.exp + (1|study_id), data = test9, family = "nbinom1")
summary(num.carbapenem.carb)
num.carbapenem.carb <- tidy(num.carbapenem.carb,
                         effects = "fixed",
                         component = "cond",
                         scales=NULL,
                         ran_prefix = NULL,
                         conf.int = TRUE,
                         conf.level = 0.95,
                         conf.method = "Wald")
num.carbapenem.carb <- num.carbapenem.carb %>%
  mutate(exposure = "carbapenem", ARG_class = "Carbapenem", model = "Number", effect=round(exp(estimate),2),
         lower=round(exp(conf.low), 2), upper = round(exp(conf.high), 2)) %>%
  filter(term != "(Intercept)") %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)

rpkm.carbapenem.fluoro <- glmer(fluoro.abund ~ carbapenem.exp + (1 |study_id), data  = test9, family=Gamma(link="log"))
summary(rpkm.carbapenem.fluoro)
rpkm.carbapenem.fluoro2 <- summary(rpkm.carbapenem.fluoro)
rpkm.carbapenem.fluoro2 <- as.data.frame(rpkm.carbapenem.fluoro2$coefficients)
rpkm.cef.com <- confint(rpkm.carbapenem.fluoro, method="Wald")
rpkm.cef.com <- rpkm.cef.com[-c(1, 2),]
rpkm.carbapenem1 <- cbind(rpkm.carbapenem.fluoro2, rpkm.cef.com)
rpkm.carbapenem1 <- rpkm.carbapenem1[-1,]
rpkm.carbapenem.fluoro <- rpkm.carbapenem1 %>% mutate(effect = round(exp(Estimate), digits =2), lower=round(exp(`2.5 %`), digits =2), upper=round(exp(`97.5 %`), digits =2)) %>%
  mutate(exposure="carbapenem", ARG_class = "Fluoroquinolone", model= "Abundance", p.value = `Pr(>|z|)`) %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)

num.carbapenem.fluoro <- glmmadmb(num.quinolone.ARG ~ carbapenem.exp + (1|study_id), data = test9, family = "nbinom1")
summary(num.carbapenem.fluoro)
num.carbapenem.fluoro <- tidy(num.carbapenem.fluoro,
                           effects = "fixed",
                           component = "cond",
                           scales=NULL,
                           ran_prefix = NULL,
                           conf.int = TRUE,
                           conf.level = 0.95,
                           conf.method = "Wald")
num.carbapenem.fluoro <- num.carbapenem.fluoro %>%
  mutate(exposure = "carbapenem", ARG_class = "Fluoroquinolone", model = "Number", effect=round(exp(estimate),2),
         lower=round(exp(conf.low), 2), upper = round(exp(conf.high), 2)) %>%
  filter(term != "(Intercept)") %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)


rpkm.carbapenem.amino <- glmer(amino.abund ~ carbapenem.exp + (1 |study_id), data  = test9, family=Gamma(link="log"))
summary(rpkm.carbapenem.amino)
rpkm.carbapenem.amino2 <- summary(rpkm.carbapenem.amino)
rpkm.carbapenem.amino2 <- as.data.frame(rpkm.carbapenem.amino2$coefficients)
rpkm.cef.com <- confint(rpkm.carbapenem.amino, method="Wald")
rpkm.cef.com <- rpkm.cef.com[-c(1, 2),]
rpkm.carbapenem1 <- cbind(rpkm.carbapenem.amino2, rpkm.cef.com)
rpkm.carbapenem1 <- rpkm.carbapenem1[-1,]
rpkm.carbapenem.amino <- rpkm.carbapenem1 %>% mutate(effect = round(exp(Estimate), digits =2), lower=round(exp(`2.5 %`), digits =2), upper=round(exp(`97.5 %`), digits =2)) %>%
  mutate(exposure="carbapenem", ARG_class = "Aminoglycoside", model= "Abundance", p.value = `Pr(>|z|)`) %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)

num.carbapenem.amino <- glmmadmb(num.amino.ARG ~ carbapenem.exp + (1|study_id), data = test9, family = "nbinom1")
summary(num.carbapenem.amino)
num.carbapenem.amino <- tidy(num.carbapenem.amino,
                          effects = "fixed",
                          component = "cond",
                          scales=NULL,
                          ran_prefix = NULL,
                          conf.int = TRUE,
                          conf.level = 0.95,
                          conf.method = "Wald")
num.carbapenem.amino <- num.carbapenem.amino %>%
  mutate(exposure = "carbapenem", ARG_class = "Aminoglycoside", model = "Number", effect=round(exp(estimate),2),
         lower=round(exp(conf.low), 2), upper = round(exp(conf.high), 2)) %>%
  filter(term != "(Intercept)") %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)


rpkm.carbapenem.clinda <- glmer(clinda.abund ~ carbapenem.exp + (1 |study_id), data  = test9, family=Gamma(link="log"))
summary(rpkm.carbapenem.clinda)
rpkm.carbapenem.clinda2 <- summary(rpkm.carbapenem.clinda)
rpkm.carbapenem.clinda2 <- as.data.frame(rpkm.carbapenem.clinda2$coefficients)
rpkm.cef.com <- confint(rpkm.carbapenem.clinda, method="Wald")
rpkm.cef.com <- rpkm.cef.com[-c(1, 2),]
rpkm.carbapenem1 <- cbind(rpkm.carbapenem.clinda2, rpkm.cef.com)
rpkm.carbapenem1 <- rpkm.carbapenem1[-1,]
rpkm.carbapenem.clinda <- rpkm.carbapenem1 %>% mutate(effect = round(exp(Estimate), digits =2), lower=round(exp(`2.5 %`), digits =2), upper=round(exp(`97.5 %`), digits =2)) %>%
  mutate(exposure="carbapenem", ARG_class = "Lincosamide", model= "Abundance", p.value = `Pr(>|z|)`) %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)

num.carbapenem.clinda <- glmmadmb(num.clinda.ARG ~ carbapenem.exp + (1|study_id), data = test9, family = "nbinom1")
summary(num.carbapenem.clinda)
num.carbapenem.clinda <- tidy(num.carbapenem.clinda,
                           effects = "fixed",
                           component = "cond",
                           scales=NULL,
                           ran_prefix = NULL,
                           conf.int = TRUE,
                           conf.level = 0.95,
                           conf.method = "Wald")
num.carbapenem.clinda <- num.carbapenem.clinda %>%
  mutate(exposure = "carbapenem", ARG_class = "Lincosamide", model = "Number", effect=round(exp(estimate),2),
         lower=round(exp(conf.low), 2), upper = round(exp(conf.high), 2)) %>%
  filter(term != "(Intercept)") %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)



rpkm.carbapenem.macrolide <- glmer(macrolide.abund ~ carbapenem.exp + (1 |study_id), data  = test9, family=Gamma(link="log"))
summary(rpkm.carbapenem.macrolide)
rpkm.carbapenem.macrolide2 <- summary(rpkm.carbapenem.macrolide)
rpkm.carbapenem.macrolide2 <- as.data.frame(rpkm.carbapenem.macrolide2$coefficients)
rpkm.cef.com <- confint(rpkm.carbapenem.macrolide, method="Wald")
rpkm.cef.com <- rpkm.cef.com[-c(1, 2),]
rpkm.carbapenem1 <- cbind(rpkm.carbapenem.macrolide2, rpkm.cef.com)
rpkm.carbapenem1 <- rpkm.carbapenem1[-1,]
rpkm.carbapenem.macrolide <- rpkm.carbapenem1 %>% mutate(effect = round(exp(Estimate), digits =2), lower=round(exp(`2.5 %`), digits =2), upper=round(exp(`97.5 %`), digits =2)) %>%
  mutate(exposure="carbapenem", ARG_class = "Macrolide", model= "Abundance", p.value = `Pr(>|z|)`) %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)

num.carbapenem.macrolide <- glmmadmb(num.macrolide.ARG ~ carbapenem.exp + (1|study_id), data = test9, family = "nbinom1")
summary(num.carbapenem.macrolide)
num.carbapenem.macrolide <- tidy(num.carbapenem.macrolide,
                              effects = "fixed",
                              component = "cond",
                              scales=NULL,
                              ran_prefix = NULL,
                              conf.int = TRUE,
                              conf.level = 0.95,
                              conf.method = "Wald")
num.carbapenem.macrolide <- num.carbapenem.macrolide %>%
  mutate(exposure = "carbapenem", ARG_class = "Macrolide", model = "Number", effect=round(exp(estimate),2),
         lower=round(exp(conf.low), 2), upper = round(exp(conf.high), 2)) %>%
  filter(term != "(Intercept)") %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)


rpkm.carbapenem.trimeth <- glmer(trimeth.abund ~ carbapenem.exp + (1 |study_id), data  = test9, family=Gamma(link="log"))
summary(rpkm.carbapenem.trimeth)
rpkm.carbapenem.trimeth2 <- summary(rpkm.carbapenem.trimeth)
rpkm.carbapenem.trimeth2 <- as.data.frame(rpkm.carbapenem.trimeth2$coefficients)
rpkm.cef.com <- confint(rpkm.carbapenem.trimeth, method="Wald")
rpkm.cef.com <- rpkm.cef.com[-c(1, 2),]
rpkm.carbapenem1 <- cbind(rpkm.carbapenem.trimeth2, rpkm.cef.com)
rpkm.carbapenem1 <- rpkm.carbapenem1[-1,]
rpkm.carbapenem.trimeth <- rpkm.carbapenem1 %>% mutate(effect = round(exp(Estimate), digits =2), lower=round(exp(`2.5 %`), digits =2), upper=round(exp(`97.5 %`), digits =2)) %>%
  mutate(exposure="carbapenem", ARG_class = "Diaminopyrimidine", model= "Abundance", p.value = `Pr(>|z|)`) %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)

num.carbapenem.trimeth <- glmmadmb(num.trimeth.ARG ~ carbapenem.exp + (1|study_id), data = test9, family = "nbinom1")
summary(num.carbapenem.trimeth)
num.carbapenem.trimeth <- tidy(num.carbapenem.trimeth,
                            effects = "fixed",
                            component = "cond",
                            scales=NULL,
                            ran_prefix = NULL,
                            conf.int = TRUE,
                            conf.level = 0.95,
                            conf.method = "Wald")
num.carbapenem.trimeth <- num.carbapenem.trimeth %>%
  mutate(exposure = "carbapenem", ARG_class = "Diaminopyrimidine", model = "Number", effect=round(exp(estimate),2),
         lower=round(exp(conf.low), 2), upper = round(exp(conf.high), 2)) %>%
  filter(term != "(Intercept)") %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)


rpkm.carbapenem.tetra <- glmer(tetra.abund ~ carbapenem.exp + (1 |study_id), data  = test9, family=Gamma(link="log"))
summary(rpkm.carbapenem.tetra)
rpkm.carbapenem.tetra2 <- summary(rpkm.carbapenem.tetra)
rpkm.carbapenem.tetra2 <- as.data.frame(rpkm.carbapenem.tetra2$coefficients)
rpkm.cef.com <- confint(rpkm.carbapenem.tetra, method="Wald")
rpkm.cef.com <- rpkm.cef.com[-c(1, 2),]
rpkm.carbapenem1 <- cbind(rpkm.carbapenem.tetra2, rpkm.cef.com)
rpkm.carbapenem1 <- rpkm.carbapenem1[-1,]
rpkm.carbapenem.tetra <- rpkm.carbapenem1 %>% mutate(effect = round(exp(Estimate), digits =2), lower=round(exp(`2.5 %`), digits =2), upper=round(exp(`97.5 %`), digits =2)) %>%
  mutate(exposure="carbapenem", ARG_class = "Tetracycline", model= "Abundance", p.value = `Pr(>|z|)`) %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)

num.carbapenem.tetra <- glmmadmb(num.tetra.ARG ~ carbapenem.exp + (1|study_id), data = test9, family = "nbinom1")
summary(num.carbapenem.tetra)
num.carbapenem.tetra <- tidy(num.carbapenem.tetra,
                          effects = "fixed",
                          component = "cond",
                          scales=NULL,
                          ran_prefix = NULL,
                          conf.int = TRUE,
                          conf.level = 0.95,
                          conf.method = "Wald")
num.carbapenem.tetra <- num.carbapenem.tetra %>%
  mutate(exposure = "carbapenem", ARG_class = "Tetracycline", model = "Number", effect=round(exp(estimate),2),
         lower=round(exp(conf.low), 2), upper = round(exp(conf.high), 2)) %>%
  filter(term != "(Intercept)") %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)


rpkm.carbapenem.flagyl <- glmer(flagyl.abund ~ carbapenem.exp + (1 |study_id), data  = test9, family=Gamma(link="log"))
summary(rpkm.carbapenem.flagyl)
rpkm.carbapenem.flagyl2 <- summary(rpkm.carbapenem.flagyl)
rpkm.carbapenem.flagyl2 <- as.data.frame(rpkm.carbapenem.flagyl2$coefficients)
rpkm.cef.com <- confint(rpkm.carbapenem.flagyl, method="Wald")
rpkm.cef.com <- rpkm.cef.com[-c(1, 2),]
rpkm.carbapenem1 <- cbind(rpkm.carbapenem.flagyl2, rpkm.cef.com)
rpkm.carbapenem1 <- rpkm.carbapenem1[-1,]
rpkm.carbapenem.flagyl <- rpkm.carbapenem1 %>% mutate(effect = round(exp(Estimate), digits =2), lower=round(exp(`2.5 %`), digits =2), upper=round(exp(`97.5 %`), digits =2)) %>%
  mutate(exposure="carbapenem", ARG_class = "Nitroimidazole", model= "Abundance", p.value = `Pr(>|z|)`) %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)

num.carbapenem.flagyl <- glmmadmb(num.flagyl.ARG ~ carbapenem.exp + (1|study_id), data = test9, family = "nbinom1")
summary(num.carbapenem.flagyl)
num.carbapenem.flagyl <- tidy(num.carbapenem.flagyl,
                           effects = "fixed",
                           component = "cond",
                           scales=NULL,
                           ran_prefix = NULL,
                           conf.int = TRUE,
                           conf.level = 0.95,
                           conf.method = "Wald")
num.carbapenem.flagyl <- num.carbapenem.flagyl %>%
  mutate(exposure = "carbapenem", ARG_class = "Nitroimidazole", model = "Number", effect=round(exp(estimate),2),
         lower=round(exp(conf.low), 2), upper = round(exp(conf.high), 2)) %>%
  filter(term != "(Intercept)") %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)

carbapenem <- rbind(rpkm.carbapenem.beta,
                 num.carbapenem.beta,
                 rpkm.carbapenem.carb,
                 num.carbapenem.carb,
                 rpkm.carbapenem.vanc,
                 num.carbapenem.vanc,
                 rpkm.carbapenem.flagyl,
                 num.carbapenem.flagyl,
                 rpkm.carbapenem.fluoro,
                 num.carbapenem.fluoro,
                 rpkm.carbapenem.trimeth,
                 num.carbapenem.trimeth,
                 rpkm.carbapenem.macrolide,
                 num.carbapenem.macrolide,
                 rpkm.carbapenem.tetra,
                 num.carbapenem.tetra,
                 rpkm.carbapenem.amino,
                 num.carbapenem.amino,
                 rpkm.carbapenem.clinda,
                 num.carbapenem.clinda)
rownames(carbapenem) <- NULL

carbapenem <- carbapenem %>%
  mutate(p.adj = p.adjust(p.value, method = "BH"))
carbapenem$p.adj <- signif(carbapenem$p.adj, 2)
carbapenem$p.adj <- format(carbapenem$p.adj, scientific = FALSE)


#vancomycin exposure
rpkm.vanc.beta <- glmer(beta.abund ~ vanc.exp + (1 |study_id), data  = test9, family=Gamma(link="log"))
summary(rpkm.vanc.beta)
rpkm.vanc.beta2 <- summary(rpkm.vanc.beta)
rpkm.vanc.beta2 <- as.data.frame(rpkm.vanc.beta2$coefficients)
rpkm.cef.com <- confint(rpkm.vanc.beta, method="Wald")
rpkm.cef.com <- rpkm.cef.com[-c(1, 2),]
rpkm.vanc1 <- cbind(rpkm.vanc.beta2, rpkm.cef.com)
rpkm.vanc1 <- rpkm.vanc1[-1,]
rpkm.vanc.beta <- rpkm.vanc1 %>% mutate(effect = round(exp(Estimate), digits =2), lower=round(exp(`2.5 %`), digits =2), upper=round(exp(`97.5 %`), digits =2)) %>%
  mutate(exposure="vancomycin", ARG_class = "Beta-Lactam", model= "Abundance", p.value = `Pr(>|z|)`) %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)

num.vanc.betalactam <- glmmadmb(num.betalactam.ARG ~ vanc.exp + (1|study_id), data = test9, family = "nbinom1")
summary(num.vanc.betalactam)
num.vanc.beta <- tidy(num.vanc.betalactam,
                          effects = "fixed",
                          component = "cond",
                          scales=NULL,
                          ran_prefix = NULL,
                          conf.int = TRUE,
                          conf.level = 0.95,
                          conf.method = "Wald")
num.vanc.beta <- num.vanc.beta %>%
  mutate(exposure = "vancomycin", ARG_class = "Beta-Lactam", model = "Number", effect=round(exp(estimate),2),
         lower=round(exp(conf.low), 2), upper = round(exp(conf.high), 2)) %>%
  filter(term != "(Intercept)") %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)

rpkm.vanc.vanc <- glmer(vanc.abund ~ vanc.exp + (1 |study_id), data  = test9, family=Gamma(link="log"))
summary(rpkm.vanc.vanc)
rpkm.vanc.vanc2 <- summary(rpkm.vanc.vanc)
rpkm.vanc.vanc2 <- as.data.frame(rpkm.vanc.vanc2$coefficients)
rpkm.cef.com <- confint(rpkm.vanc.vanc, method="Wald")
rpkm.cef.com <- rpkm.cef.com[-c(1, 2),]
rpkm.vanc1 <- cbind(rpkm.vanc.vanc2, rpkm.cef.com)
rpkm.vanc1 <- rpkm.vanc1[-1,]
rpkm.vanc.vanc <- rpkm.vanc1 %>% mutate(effect = round(exp(Estimate), digits =2), lower=round(exp(`2.5 %`), digits =2), upper=round(exp(`97.5 %`), digits =2)) %>%
  mutate(exposure="vancomycin", ARG_class = "Glycopeptide", model= "Abundance", p.value = `Pr(>|z|)`) %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)

num.vanc.vanc <- glmmadmb(num.vanc.ARG ~ vanc.exp + (1|study_id), data = test9, family = "nbinom1")
summary(num.vanc.vanc)
num.vanc.vanc <- tidy(num.vanc.vanc,
                          effects = "fixed",
                          component = "cond",
                          scales=NULL,
                          ran_prefix = NULL,
                          conf.int = TRUE,
                          conf.level = 0.95,
                          conf.method = "Wald")
num.vanc.vanc <- num.vanc.vanc %>%
  mutate(exposure = "vancomycin", ARG_class = "Glycopeptide", model = "Number", effect=round(exp(estimate),2),
         lower=round(exp(conf.low), 2), upper = round(exp(conf.high), 2)) %>%
  filter(term != "(Intercept)") %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)


rpkm.vanc.carb <- glmer(carb.abund ~ vanc.exp + (1 |study_id), data  = test9, family=Gamma(link="log"))
summary(rpkm.vanc.carb)
rpkm.vanc.carb2 <- summary(rpkm.vanc.carb)
rpkm.vanc.carb2 <- as.data.frame(rpkm.vanc.carb2$coefficients)
rpkm.cef.com <- confint(rpkm.vanc.carb, method="Wald")
rpkm.cef.com <- rpkm.cef.com[-c(1, 2),]
rpkm.vanc1 <- cbind(rpkm.vanc.carb2, rpkm.cef.com)
rpkm.vanc1 <- rpkm.vanc1[-1,]
rpkm.vanc.carb <- rpkm.vanc1 %>% mutate(effect = round(exp(Estimate), digits =2), lower=round(exp(`2.5 %`), digits =2), upper=round(exp(`97.5 %`), digits =2)) %>%
  mutate(exposure="vancomycin", ARG_class = "Carbapenem", model= "Abundance", p.value = `Pr(>|z|)`) %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)

num.vanc.carb <- glmmadmb(num.carbapenem.ARG ~ vanc.exp + (1|study_id), data = test9, family = "nbinom1")
summary(num.vanc.carb)
num.vanc.carb <- tidy(num.vanc.carb,
                          effects = "fixed",
                          component = "cond",
                          scales=NULL,
                          ran_prefix = NULL,
                          conf.int = TRUE,
                          conf.level = 0.95,
                          conf.method = "Wald")
num.vanc.carb <- num.vanc.carb %>%
  mutate(exposure = "vancomycin", ARG_class = "Carbapenem", model = "Number", effect=round(exp(estimate),2),
         lower=round(exp(conf.low), 2), upper = round(exp(conf.high), 2)) %>%
  filter(term != "(Intercept)") %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)

rpkm.vanc.fluoro <- glmer(fluoro.abund ~ vanc.exp + (1 |study_id), data  = test9, family=Gamma(link="log"))
summary(rpkm.vanc.fluoro)
rpkm.vanc.fluoro2 <- summary(rpkm.vanc.fluoro)
rpkm.vanc.fluoro2 <- as.data.frame(rpkm.vanc.fluoro2$coefficients)
rpkm.cef.com <- confint(rpkm.vanc.fluoro, method="Wald")
rpkm.cef.com <- rpkm.cef.com[-c(1, 2),]
rpkm.vanc1 <- cbind(rpkm.vanc.fluoro2, rpkm.cef.com)
rpkm.vanc1 <- rpkm.vanc1[-1,]
rpkm.vanc.fluoro <- rpkm.vanc1 %>% mutate(effect = round(exp(Estimate), digits =2), lower=round(exp(`2.5 %`), digits =2), upper=round(exp(`97.5 %`), digits =2)) %>%
  mutate(exposure="vancomycin", ARG_class = "Fluoroquinolone", model= "Abundance", p.value = `Pr(>|z|)`) %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)

num.vanc.fluoro <- glmmadmb(num.quinolone.ARG ~ vanc.exp + (1|study_id), data = test9, family = "nbinom1")
summary(num.vanc.fluoro)
num.vanc.fluoro <- tidy(num.vanc.fluoro,
                            effects = "fixed",
                            component = "cond",
                            scales=NULL,
                            ran_prefix = NULL,
                            conf.int = TRUE,
                            conf.level = 0.95,
                            conf.method = "Wald")
num.vanc.fluoro <- num.vanc.fluoro %>%
  mutate(exposure = "vancomycin", ARG_class = "Fluoroquinolone", model = "Number", effect=round(exp(estimate),2),
         lower=round(exp(conf.low), 2), upper = round(exp(conf.high), 2)) %>%
  filter(term != "(Intercept)") %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)


rpkm.vanc.amino <- glmer(amino.abund ~ vanc.exp + (1 |study_id), data  = test9, family=Gamma(link="log"))
summary(rpkm.vanc.amino)
rpkm.vanc.amino2 <- summary(rpkm.vanc.amino)
rpkm.vanc.amino2 <- as.data.frame(rpkm.vanc.amino2$coefficients)
rpkm.cef.com <- confint(rpkm.vanc.amino, method="Wald")
rpkm.cef.com <- rpkm.cef.com[-c(1, 2),]
rpkm.vanc1 <- cbind(rpkm.vanc.amino2, rpkm.cef.com)
rpkm.vanc1 <- rpkm.vanc1[-1,]
rpkm.vanc.amino <- rpkm.vanc1 %>% mutate(effect = round(exp(Estimate), digits =2), lower=round(exp(`2.5 %`), digits =2), upper=round(exp(`97.5 %`), digits =2)) %>%
  mutate(exposure="vancomycin", ARG_class = "Aminoglycoside", model= "Abundance", p.value = `Pr(>|z|)`) %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)

num.vanc.amino <- glmmadmb(num.amino.ARG ~ vanc.exp + (1|study_id), data = test9, family = "nbinom1")
summary(num.vanc.amino)
num.vanc.amino <- tidy(num.vanc.amino,
                           effects = "fixed",
                           component = "cond",
                           scales=NULL,
                           ran_prefix = NULL,
                           conf.int = TRUE,
                           conf.level = 0.95,
                           conf.method = "Wald")
num.vanc.amino <- num.vanc.amino %>%
  mutate(exposure = "vancomycin", ARG_class = "Aminoglycoside", model = "Number", effect=round(exp(estimate),2),
         lower=round(exp(conf.low), 2), upper = round(exp(conf.high), 2)) %>%
  filter(term != "(Intercept)") %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)


rpkm.vanc.clinda <- glmer(clinda.abund ~ vanc.exp + (1 |study_id), data  = test9, family=Gamma(link="log"))
summary(rpkm.vanc.clinda)
rpkm.vanc.clinda2 <- summary(rpkm.vanc.clinda)
rpkm.vanc.clinda2 <- as.data.frame(rpkm.vanc.clinda2$coefficients)
rpkm.cef.com <- confint(rpkm.vanc.clinda, method="Wald")
rpkm.cef.com <- rpkm.cef.com[-c(1, 2),]
rpkm.vanc1 <- cbind(rpkm.vanc.clinda2, rpkm.cef.com)
rpkm.vanc1 <- rpkm.vanc1[-1,]
rpkm.vanc.clinda <- rpkm.vanc1 %>% mutate(effect = round(exp(Estimate), digits =2), lower=round(exp(`2.5 %`), digits =2), upper=round(exp(`97.5 %`), digits =2)) %>%
  mutate(exposure="vancomycin", ARG_class = "Lincosamide", model= "Abundance", p.value = `Pr(>|z|)`) %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)

num.vanc.clinda <- glmmadmb(num.clinda.ARG ~ vanc.exp + (1|study_id), data = test9, family = "nbinom1")
summary(num.vanc.clinda)
num.vanc.clinda <- tidy(num.vanc.clinda,
                            effects = "fixed",
                            component = "cond",
                            scales=NULL,
                            ran_prefix = NULL,
                            conf.int = TRUE,
                            conf.level = 0.95,
                            conf.method = "Wald")
num.vanc.clinda <- num.vanc.clinda %>%
  mutate(exposure = "vancomycin", ARG_class = "Lincosamide", model = "Number", effect=round(exp(estimate),2),
         lower=round(exp(conf.low), 2), upper = round(exp(conf.high), 2)) %>%
  filter(term != "(Intercept)") %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)



rpkm.vanc.macrolide <- glmer(macrolide.abund ~ vanc.exp + (1 |study_id), data  = test9, family=Gamma(link="log"))
summary(rpkm.vanc.macrolide)
rpkm.vanc.macrolide2 <- summary(rpkm.vanc.macrolide)
rpkm.vanc.macrolide2 <- as.data.frame(rpkm.vanc.macrolide2$coefficients)
rpkm.cef.com <- confint(rpkm.vanc.macrolide, method="Wald")
rpkm.cef.com <- rpkm.cef.com[-c(1, 2),]
rpkm.vanc1 <- cbind(rpkm.vanc.macrolide2, rpkm.cef.com)
rpkm.vanc1 <- rpkm.vanc1[-1,]
rpkm.vanc.macrolide <- rpkm.vanc1 %>% mutate(effect = round(exp(Estimate), digits =2), lower=round(exp(`2.5 %`), digits =2), upper=round(exp(`97.5 %`), digits =2)) %>%
  mutate(exposure="vancomycin", ARG_class = "Macrolide", model= "Abundance", p.value = `Pr(>|z|)`) %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)

num.vanc.macrolide <- glmmadmb(num.macrolide.ARG ~ vanc.exp + (1|study_id), data = test9, family = "nbinom1")
summary(num.vanc.macrolide)
num.vanc.macrolide <- tidy(num.vanc.macrolide,
                               effects = "fixed",
                               component = "cond",
                               scales=NULL,
                               ran_prefix = NULL,
                               conf.int = TRUE,
                               conf.level = 0.95,
                               conf.method = "Wald")
num.vanc.macrolide <- num.vanc.macrolide %>%
  mutate(exposure = "vancomycin", ARG_class = "Macrolide", model = "Number", effect=round(exp(estimate),2),
         lower=round(exp(conf.low), 2), upper = round(exp(conf.high), 2)) %>%
  filter(term != "(Intercept)") %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)


rpkm.vanc.trimeth <- glmer(trimeth.abund ~ vanc.exp + (1 |study_id), data  = test9, family=Gamma(link="log"))
summary(rpkm.vanc.trimeth)
rpkm.vanc.trimeth2 <- summary(rpkm.vanc.trimeth)
rpkm.vanc.trimeth2 <- as.data.frame(rpkm.vanc.trimeth2$coefficients)
rpkm.cef.com <- confint(rpkm.vanc.trimeth, method="Wald")
rpkm.cef.com <- rpkm.cef.com[-c(1, 2),]
rpkm.vanc1 <- cbind(rpkm.vanc.trimeth2, rpkm.cef.com)
rpkm.vanc1 <- rpkm.vanc1[-1,]
rpkm.vanc.trimeth <- rpkm.vanc1 %>% mutate(effect = round(exp(Estimate), digits =2), lower=round(exp(`2.5 %`), digits =2), upper=round(exp(`97.5 %`), digits =2)) %>%
  mutate(exposure="vancomycin", ARG_class = "Diaminopyrimidine", model= "Abundance", p.value = `Pr(>|z|)`) %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)

num.vanc.trimeth <- glmmadmb(num.trimeth.ARG ~ vanc.exp + (1|study_id), data = test9, family = "nbinom1")
summary(num.vanc.trimeth)
num.vanc.trimeth <- tidy(num.vanc.trimeth,
                             effects = "fixed",
                             component = "cond",
                             scales=NULL,
                             ran_prefix = NULL,
                             conf.int = TRUE,
                             conf.level = 0.95,
                             conf.method = "Wald")
num.vanc.trimeth <- num.vanc.trimeth %>%
  mutate(exposure = "vancomycin", ARG_class = "Diaminopyrimidine", model = "Number", effect=round(exp(estimate),2),
         lower=round(exp(conf.low), 2), upper = round(exp(conf.high), 2)) %>%
  filter(term != "(Intercept)") %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)


rpkm.vanc.tetra <- glmer(tetra.abund ~ vanc.exp + (1 |study_id), data  = test9, family=Gamma(link="log"))
summary(rpkm.vanc.tetra)
rpkm.vanc.tetra2 <- summary(rpkm.vanc.tetra)
rpkm.vanc.tetra2 <- as.data.frame(rpkm.vanc.tetra2$coefficients)
rpkm.cef.com <- confint(rpkm.vanc.tetra, method="Wald")
rpkm.cef.com <- rpkm.cef.com[-c(1, 2),]
rpkm.vanc1 <- cbind(rpkm.vanc.tetra2, rpkm.cef.com)
rpkm.vanc1 <- rpkm.vanc1[-1,]
rpkm.vanc.tetra <- rpkm.vanc1 %>% mutate(effect = round(exp(Estimate), digits =2), lower=round(exp(`2.5 %`), digits =2), upper=round(exp(`97.5 %`), digits =2)) %>%
  mutate(exposure="vancomycin", ARG_class = "Tetracycline", model= "Abundance", p.value = `Pr(>|z|)`) %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)

num.vanc.tetra <- glmmadmb(num.tetra.ARG ~ vanc.exp + (1|study_id), data = test9, family = "nbinom1")
summary(num.vanc.tetra)
num.vanc.tetra <- tidy(num.vanc.tetra,
                           effects = "fixed",
                           component = "cond",
                           scales=NULL,
                           ran_prefix = NULL,
                           conf.int = TRUE,
                           conf.level = 0.95,
                           conf.method = "Wald")
num.vanc.tetra <- num.vanc.tetra %>%
  mutate(exposure = "vancomycin", ARG_class = "Tetracycline", model = "Number", effect=round(exp(estimate),2),
         lower=round(exp(conf.low), 2), upper = round(exp(conf.high), 2)) %>%
  filter(term != "(Intercept)") %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)


rpkm.vanc.flagyl <- glmer(flagyl.abund ~ vanc.exp + (1 |study_id), data  = test9, family=Gamma(link="log"))
summary(rpkm.vanc.flagyl)
rpkm.vanc.flagyl2 <- summary(rpkm.vanc.flagyl)
rpkm.vanc.flagyl2 <- as.data.frame(rpkm.vanc.flagyl2$coefficients)
rpkm.cef.com <- confint(rpkm.vanc.flagyl, method="Wald")
rpkm.cef.com <- rpkm.cef.com[-c(1, 2),]
rpkm.vanc1 <- cbind(rpkm.vanc.flagyl2, rpkm.cef.com)
rpkm.vanc1 <- rpkm.vanc1[-1,]
rpkm.vanc.flagyl <- rpkm.vanc1 %>% mutate(effect = round(exp(Estimate), digits =2), lower=round(exp(`2.5 %`), digits =2), upper=round(exp(`97.5 %`), digits =2)) %>%
  mutate(exposure="vancomycin", ARG_class = "Nitroimidazole", model= "Abundance", p.value = `Pr(>|z|)`) %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)

num.vanc.flagyl <- glmmadmb(num.flagyl.ARG ~ vanc.exp + (1|study_id), data = test9, family = "nbinom1")
summary(num.vanc.flagyl)
num.vanc.flagyl <- tidy(num.vanc.flagyl,
                            effects = "fixed",
                            component = "cond",
                            scales=NULL,
                            ran_prefix = NULL,
                            conf.int = TRUE,
                            conf.level = 0.95,
                            conf.method = "Wald")
num.vanc.flagyl <- num.vanc.flagyl %>%
  mutate(exposure = "vancomycin", ARG_class = "Nitroimidazole", model = "Number", effect=round(exp(estimate),2),
         lower=round(exp(conf.low), 2), upper = round(exp(conf.high), 2)) %>%
  filter(term != "(Intercept)") %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)

vanc <- rbind(rpkm.vanc.beta,
                  num.vanc.beta,
                  rpkm.vanc.carb,
                  num.vanc.carb,
                  rpkm.vanc.vanc,
                  num.vanc.vanc,
                  rpkm.vanc.flagyl,
                  num.vanc.flagyl,
                  rpkm.vanc.fluoro,
                  num.vanc.fluoro,
                  rpkm.vanc.trimeth,
                  num.vanc.trimeth,
                  rpkm.vanc.macrolide,
                  num.vanc.macrolide,
                  rpkm.vanc.tetra,
                  num.vanc.tetra,
                  rpkm.vanc.amino,
                  num.vanc.amino,
                  rpkm.vanc.clinda,
                  num.vanc.clinda)
rownames(vanc) <- NULL

vanc <- vanc %>%
  mutate(p.adj = p.adjust(p.value, method = "BH"))
vanc$p.adj <- signif(vanc$p.adj, 2)
vanc$p.adj <- format(vanc$p.adj, scientific = FALSE)

#flagyl exposure
rpkm.flagyl.beta <- glmer(beta.abund ~ flagyl.exp + (1 |study_id), data  = test9, family=Gamma(link="log"))
summary(rpkm.flagyl.beta)
rpkm.flagyl.beta2 <- summary(rpkm.flagyl.beta)
rpkm.flagyl.beta2 <- as.data.frame(rpkm.flagyl.beta2$coefficients)
rpkm.cef.com <- confint(rpkm.flagyl.beta, method="Wald")
rpkm.cef.com <- rpkm.cef.com[-c(1, 2),]
rpkm.flagyl1 <- cbind(rpkm.flagyl.beta2, rpkm.cef.com)
rpkm.flagyl1 <- rpkm.flagyl1[-1,]
rpkm.flagyl.beta <- rpkm.flagyl1 %>% mutate(effect = round(exp(Estimate), digits =2), lower=round(exp(`2.5 %`), digits =2), upper=round(exp(`97.5 %`), digits =2)) %>%
  mutate(exposure="flagyl", ARG_class = "Beta-Lactam", model= "Abundance", p.value = `Pr(>|z|)`) %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)

num.flagyl.betalactam <- glmmadmb(num.betalactam.ARG ~ flagyl.exp + (1|study_id), data = test9, family = "nbinom1")
summary(num.flagyl.betalactam)
num.flagyl.beta <- tidy(num.flagyl.betalactam,
                          effects = "fixed",
                          component = "cond",
                          scales=NULL,
                          ran_prefix = NULL,
                          conf.int = TRUE,
                          conf.level = 0.95,
                          conf.method = "Wald")
num.flagyl.beta <- num.flagyl.beta %>%
  mutate(exposure = "flagyl", ARG_class = "Beta-Lactam", model = "Number", effect=round(exp(estimate),2),
         lower=round(exp(conf.low), 2), upper = round(exp(conf.high), 2)) %>%
  filter(term != "(Intercept)") %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)

rpkm.flagyl.vanc <- glmer(vanc.abund ~ flagyl.exp + (1 |study_id), data  = test9, family=Gamma(link="log"))
summary(rpkm.flagyl.vanc)
rpkm.flagyl.vanc2 <- summary(rpkm.flagyl.vanc)
rpkm.flagyl.vanc2 <- as.data.frame(rpkm.flagyl.vanc2$coefficients)
rpkm.cef.com <- confint(rpkm.flagyl.vanc, method="Wald")
rpkm.cef.com <- rpkm.cef.com[-c(1, 2),]
rpkm.flagyl1 <- cbind(rpkm.flagyl.vanc2, rpkm.cef.com)
rpkm.flagyl1 <- rpkm.flagyl1[-1,]
rpkm.flagyl.vanc <- rpkm.flagyl1 %>% mutate(effect = round(exp(Estimate), digits =2), lower=round(exp(`2.5 %`), digits =2), upper=round(exp(`97.5 %`), digits =2)) %>%
  mutate(exposure="flagyl", ARG_class = "Glycopeptide", model= "Abundance", p.value = `Pr(>|z|)`) %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)

num.flagyl.vanc <- glmmadmb(num.vanc.ARG ~ flagyl.exp + (1|study_id), data = test9, family = "nbinom1")
summary(num.flagyl.vanc)
num.flagyl.vanc <- tidy(num.flagyl.vanc,
                          effects = "fixed",
                          component = "cond",
                          scales=NULL,
                          ran_prefix = NULL,
                          conf.int = TRUE,
                          conf.level = 0.95,
                          conf.method = "Wald")
num.flagyl.vanc <- num.flagyl.vanc %>%
  mutate(exposure = "flagyl", ARG_class = "Glycopeptide", model = "Number", effect=round(exp(estimate),2),
         lower=round(exp(conf.low), 2), upper = round(exp(conf.high), 2)) %>%
  filter(term != "(Intercept)") %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)


rpkm.flagyl.carb <- glmer(carb.abund ~ flagyl.exp + (1 |study_id), data  = test9, family=Gamma(link="log"))
summary(rpkm.flagyl.carb)
rpkm.flagyl.carb2 <- summary(rpkm.flagyl.carb)
rpkm.flagyl.carb2 <- as.data.frame(rpkm.flagyl.carb2$coefficients)
rpkm.cef.com <- confint(rpkm.flagyl.carb, method="Wald")
rpkm.cef.com <- rpkm.cef.com[-c(1, 2),]
rpkm.flagyl1 <- cbind(rpkm.flagyl.carb2, rpkm.cef.com)
rpkm.flagyl1 <- rpkm.flagyl1[-1,]
rpkm.flagyl.carb <- rpkm.flagyl1 %>% mutate(effect = round(exp(Estimate), digits =2), lower=round(exp(`2.5 %`), digits =2), upper=round(exp(`97.5 %`), digits =2)) %>%
  mutate(exposure="flagyl", ARG_class = "Carbapenem", model= "Abundance", p.value = `Pr(>|z|)`) %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)

num.flagyl.carb <- glmmadmb(num.carbapenem.ARG ~ flagyl.exp + (1|study_id), data = test9, family = "nbinom1")
summary(num.flagyl.carb)
num.flagyl.carb <- tidy(num.flagyl.carb,
                          effects = "fixed",
                          component = "cond",
                          scales=NULL,
                          ran_prefix = NULL,
                          conf.int = TRUE,
                          conf.level = 0.95,
                          conf.method = "Wald")
num.flagyl.carb <- num.flagyl.carb %>%
  mutate(exposure = "flagyl", ARG_class = "Carbapenem", model = "Number", effect=round(exp(estimate),2),
         lower=round(exp(conf.low), 2), upper = round(exp(conf.high), 2)) %>%
  filter(term != "(Intercept)") %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)

rpkm.flagyl.fluoro <- glmer(fluoro.abund ~ flagyl.exp + (1 |study_id), data  = test9, family=Gamma(link="log"))
summary(rpkm.flagyl.fluoro)
rpkm.flagyl.fluoro2 <- summary(rpkm.flagyl.fluoro)
rpkm.flagyl.fluoro2 <- as.data.frame(rpkm.flagyl.fluoro2$coefficients)
rpkm.cef.com <- confint(rpkm.flagyl.fluoro, method="Wald")
rpkm.cef.com <- rpkm.cef.com[-c(1, 2),]
rpkm.flagyl1 <- cbind(rpkm.flagyl.fluoro2, rpkm.cef.com)
rpkm.flagyl1 <- rpkm.flagyl1[-1,]
rpkm.flagyl.fluoro <- rpkm.flagyl1 %>% mutate(effect = round(exp(Estimate), digits =2), lower=round(exp(`2.5 %`), digits =2), upper=round(exp(`97.5 %`), digits =2)) %>%
  mutate(exposure="flagyl", ARG_class = "Fluoroquinolone", model= "Abundance", p.value = `Pr(>|z|)`) %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)

num.flagyl.fluoro <- glmmadmb(num.quinolone.ARG ~ flagyl.exp + (1|study_id), data = test9, family = "nbinom1")
summary(num.flagyl.fluoro)
num.flagyl.fluoro <- tidy(num.flagyl.fluoro,
                            effects = "fixed",
                            component = "cond",
                            scales=NULL,
                            ran_prefix = NULL,
                            conf.int = TRUE,
                            conf.level = 0.95,
                            conf.method = "Wald")
num.flagyl.fluoro <- num.flagyl.fluoro %>%
  mutate(exposure = "flagyl", ARG_class = "Fluoroquinolone", model = "Number", effect=round(exp(estimate),2),
         lower=round(exp(conf.low), 2), upper = round(exp(conf.high), 2)) %>%
  filter(term != "(Intercept)") %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)


rpkm.flagyl.amino <- glmer(amino.abund ~ flagyl.exp + (1 |study_id), data  = test9, family=Gamma(link="log"))
summary(rpkm.flagyl.amino)
rpkm.flagyl.amino2 <- summary(rpkm.flagyl.amino)
rpkm.flagyl.amino2 <- as.data.frame(rpkm.flagyl.amino2$coefficients)
rpkm.cef.com <- confint(rpkm.flagyl.amino, method="Wald")
rpkm.cef.com <- rpkm.cef.com[-c(1, 2),]
rpkm.flagyl1 <- cbind(rpkm.flagyl.amino2, rpkm.cef.com)
rpkm.flagyl1 <- rpkm.flagyl1[-1,]
rpkm.flagyl.amino <- rpkm.flagyl1 %>% mutate(effect = round(exp(Estimate), digits =2), lower=round(exp(`2.5 %`), digits =2), upper=round(exp(`97.5 %`), digits =2)) %>%
  mutate(exposure="flagyl", ARG_class = "Aminoglycoside", model= "Abundance", p.value = `Pr(>|z|)`) %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)

num.flagyl.amino <- glmmadmb(num.amino.ARG ~ flagyl.exp + (1|study_id), data = test9, family = "nbinom1")
summary(num.flagyl.amino)
num.flagyl.amino <- tidy(num.flagyl.amino,
                           effects = "fixed",
                           component = "cond",
                           scales=NULL,
                           ran_prefix = NULL,
                           conf.int = TRUE,
                           conf.level = 0.95,
                           conf.method = "Wald")
num.flagyl.amino <- num.flagyl.amino %>%
  mutate(exposure = "flagyl", ARG_class = "Aminoglycoside", model = "Number", effect=round(exp(estimate),2),
         lower=round(exp(conf.low), 2), upper = round(exp(conf.high), 2)) %>%
  filter(term != "(Intercept)") %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)


rpkm.flagyl.clinda <- glmer(clinda.abund ~ flagyl.exp + (1 |study_id), data  = test9, family=Gamma(link="log"))
summary(rpkm.flagyl.clinda)
rpkm.flagyl.clinda2 <- summary(rpkm.flagyl.clinda)
rpkm.flagyl.clinda2 <- as.data.frame(rpkm.flagyl.clinda2$coefficients)
rpkm.cef.com <- confint(rpkm.flagyl.clinda, method="Wald")
rpkm.cef.com <- rpkm.cef.com[-c(1, 2),]
rpkm.flagyl1 <- cbind(rpkm.flagyl.clinda2, rpkm.cef.com)
rpkm.flagyl1 <- rpkm.flagyl1[-1,]
rpkm.flagyl.clinda <- rpkm.flagyl1 %>% mutate(effect = round(exp(Estimate), digits =2), lower=round(exp(`2.5 %`), digits =2), upper=round(exp(`97.5 %`), digits =2)) %>%
  mutate(exposure="flagyl", ARG_class = "Lincosamide", model= "Abundance", p.value = `Pr(>|z|)`) %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)

num.flagyl.clinda <- glmmadmb(num.clinda.ARG ~ flagyl.exp + (1|study_id), data = test9, family = "nbinom1")
summary(num.flagyl.clinda)
num.flagyl.clinda <- tidy(num.flagyl.clinda,
                            effects = "fixed",
                            component = "cond",
                            scales=NULL,
                            ran_prefix = NULL,
                            conf.int = TRUE,
                            conf.level = 0.95,
                            conf.method = "Wald")
num.flagyl.clinda <- num.flagyl.clinda %>%
  mutate(exposure = "flagyl", ARG_class = "Lincosamide", model = "Number", effect=round(exp(estimate),2),
         lower=round(exp(conf.low), 2), upper = round(exp(conf.high), 2)) %>%
  filter(term != "(Intercept)") %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)



rpkm.flagyl.macrolide <- glmer(macrolide.abund ~ flagyl.exp + (1 |study_id), data  = test9, family=Gamma(link="log"))
summary(rpkm.flagyl.macrolide)
rpkm.flagyl.macrolide2 <- summary(rpkm.flagyl.macrolide)
rpkm.flagyl.macrolide2 <- as.data.frame(rpkm.flagyl.macrolide2$coefficients)
rpkm.cef.com <- confint(rpkm.flagyl.macrolide, method="Wald")
rpkm.cef.com <- rpkm.cef.com[-c(1, 2),]
rpkm.flagyl1 <- cbind(rpkm.flagyl.macrolide2, rpkm.cef.com)
rpkm.flagyl1 <- rpkm.flagyl1[-1,]
rpkm.flagyl.macrolide <- rpkm.flagyl1 %>% mutate(effect = round(exp(Estimate), digits =2), lower=round(exp(`2.5 %`), digits =2), upper=round(exp(`97.5 %`), digits =2)) %>%
  mutate(exposure="flagyl", ARG_class = "Macrolide", model= "Abundance", p.value = `Pr(>|z|)`) %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)

num.flagyl.macrolide <- glmmadmb(num.macrolide.ARG ~ flagyl.exp + (1|study_id), data = test9, family = "nbinom1")
summary(num.flagyl.macrolide)
num.flagyl.macrolide <- tidy(num.flagyl.macrolide,
                               effects = "fixed",
                               component = "cond",
                               scales=NULL,
                               ran_prefix = NULL,
                               conf.int = TRUE,
                               conf.level = 0.95,
                               conf.method = "Wald")
num.flagyl.macrolide <- num.flagyl.macrolide %>%
  mutate(exposure = "flagyl", ARG_class = "Macrolide", model = "Number", effect=round(exp(estimate),2),
         lower=round(exp(conf.low), 2), upper = round(exp(conf.high), 2)) %>%
  filter(term != "(Intercept)") %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)


rpkm.flagyl.trimeth <- glmer(trimeth.abund ~ flagyl.exp + (1 |study_id), data  = test9, family=Gamma(link="log"))
summary(rpkm.flagyl.trimeth)
rpkm.flagyl.trimeth2 <- summary(rpkm.flagyl.trimeth)
rpkm.flagyl.trimeth2 <- as.data.frame(rpkm.flagyl.trimeth2$coefficients)
rpkm.cef.com <- confint(rpkm.flagyl.trimeth, method="Wald")
rpkm.cef.com <- rpkm.cef.com[-c(1, 2),]
rpkm.flagyl1 <- cbind(rpkm.flagyl.trimeth2, rpkm.cef.com)
rpkm.flagyl1 <- rpkm.flagyl1[-1,]
rpkm.flagyl.trimeth <- rpkm.flagyl1 %>% mutate(effect = round(exp(Estimate), digits =2), lower=round(exp(`2.5 %`), digits =2), upper=round(exp(`97.5 %`), digits =2)) %>%
  mutate(exposure="flagyl", ARG_class = "Diaminopyrimidine", model= "Abundance", p.value = `Pr(>|z|)`) %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)

num.flagyl.trimeth <- glmmadmb(num.trimeth.ARG ~ flagyl.exp + (1|study_id), data = test9, family = "nbinom2")
summary(num.flagyl.trimeth)
num.flagyl.trimeth <- tidy(num.flagyl.trimeth,
                             effects = "fixed",
                             component = "cond",
                             scales=NULL,
                             ran_prefix = NULL,
                             conf.int = TRUE,
                             conf.level = 0.95,
                             conf.method = "Wald")
num.flagyl.trimeth <- num.flagyl.trimeth %>%
  mutate(exposure = "flagyl", ARG_class = "Diaminopyrimidine", model = "Number", effect=round(exp(estimate),2),
         lower=round(exp(conf.low), 2), upper = round(exp(conf.high), 2)) %>%
  filter(term != "(Intercept)") %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)


rpkm.flagyl.tetra <- glmer(tetra.abund ~ flagyl.exp + (1 |study_id), data  = test9, family=Gamma(link="log"))
summary(rpkm.flagyl.tetra)
rpkm.flagyl.tetra2 <- summary(rpkm.flagyl.tetra)
rpkm.flagyl.tetra2 <- as.data.frame(rpkm.flagyl.tetra2$coefficients)
rpkm.cef.com <- confint(rpkm.flagyl.tetra, method="Wald")
rpkm.cef.com <- rpkm.cef.com[-c(1, 2),]
rpkm.flagyl1 <- cbind(rpkm.flagyl.tetra2, rpkm.cef.com)
rpkm.flagyl1 <- rpkm.flagyl1[-1,]
rpkm.flagyl.tetra <- rpkm.flagyl1 %>% mutate(effect = round(exp(Estimate), digits =2), lower=round(exp(`2.5 %`), digits =2), upper=round(exp(`97.5 %`), digits =2)) %>%
  mutate(exposure="flagyl", ARG_class = "Tetracycline", model= "Abundance", p.value = `Pr(>|z|)`) %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)

num.flagyl.tetra <- glmmadmb(num.tetra.ARG ~ flagyl.exp + (1|study_id), data = test9, family = "nbinom1")
summary(num.flagyl.tetra)
num.flagyl.tetra <- tidy(num.flagyl.tetra,
                           effects = "fixed",
                           component = "cond",
                           scales=NULL,
                           ran_prefix = NULL,
                           conf.int = TRUE,
                           conf.level = 0.95,
                           conf.method = "Wald")
num.flagyl.tetra <- num.flagyl.tetra %>%
  mutate(exposure = "flagyl", ARG_class = "Tetracycline", model = "Number", effect=round(exp(estimate),2),
         lower=round(exp(conf.low), 2), upper = round(exp(conf.high), 2)) %>%
  filter(term != "(Intercept)") %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)


rpkm.flagyl.flagyl <- glmer(flagyl.abund ~ flagyl.exp + (1 |study_id), data  = test9, family=Gamma(link="log"))
summary(rpkm.flagyl.flagyl)
rpkm.flagyl.flagyl2 <- summary(rpkm.flagyl.flagyl)
rpkm.flagyl.flagyl2 <- as.data.frame(rpkm.flagyl.flagyl2$coefficients)
rpkm.cef.com <- confint(rpkm.flagyl.flagyl, method="Wald")
rpkm.cef.com <- rpkm.cef.com[-c(1, 2),]
rpkm.flagyl1 <- cbind(rpkm.flagyl.flagyl2, rpkm.cef.com)
rpkm.flagyl1 <- rpkm.flagyl1[-1,]
rpkm.flagyl.flagyl <- rpkm.flagyl1 %>% mutate(effect = round(exp(Estimate), digits =2), lower=round(exp(`2.5 %`), digits =2), upper=round(exp(`97.5 %`), digits =2)) %>%
  mutate(exposure="flagyl", ARG_class = "Nitroimidazole", model= "Abundance", p.value = `Pr(>|z|)`) %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)

num.flagyl.flagyl <- glmmadmb(num.flagyl.ARG ~ flagyl.exp + (1|study_id), data = test9, family = "nbinom1")
summary(num.flagyl.flagyl)
num.flagyl.flagyl <- tidy(num.flagyl.flagyl,
                            effects = "fixed",
                            component = "cond",
                            scales=NULL,
                            ran_prefix = NULL,
                            conf.int = TRUE,
                            conf.level = 0.95,
                            conf.method = "Wald")
num.flagyl.flagyl <- num.flagyl.flagyl %>%
  mutate(exposure = "flagyl", ARG_class = "Nitroimidazole", model = "Number", effect=round(exp(estimate),2),
         lower=round(exp(conf.low), 2), upper = round(exp(conf.high), 2)) %>%
  filter(term != "(Intercept)") %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)

flagyl <- rbind(rpkm.flagyl.beta,
                  num.flagyl.beta,
                  rpkm.flagyl.carb,
                  num.flagyl.carb,
                  rpkm.flagyl.vanc,
                  num.flagyl.vanc,
                  rpkm.flagyl.flagyl,
                  num.flagyl.flagyl,
                  rpkm.flagyl.fluoro,
                  num.flagyl.fluoro,
                  rpkm.flagyl.trimeth,
                  num.flagyl.trimeth,
                  rpkm.flagyl.macrolide,
                  num.flagyl.macrolide,
                  rpkm.flagyl.tetra,
                  num.flagyl.tetra,
                  rpkm.flagyl.amino,
                  num.flagyl.amino,
                  rpkm.flagyl.clinda,
                  num.flagyl.clinda)
rownames(flagyl) <- NULL

flagyl <- flagyl %>%
  mutate(p.adj = p.adjust(p.value, method = "BH"))
flagyl$p.adj <- signif(flagyl$p.adj, 2)
flagyl$p.adj <- format(flagyl$p.adj, scientific = FALSE)

#quinolone exposure
rpkm.quinolone.beta <- glmer(beta.abund ~ quinolone.exp + (1 |study_id), data  = test9, family=Gamma(link="log"))
summary(rpkm.quinolone.beta)
rpkm.quinolone.beta2 <- summary(rpkm.quinolone.beta)
rpkm.quinolone.beta2 <- as.data.frame(rpkm.quinolone.beta2$coefficients)
rpkm.cef.com <- confint(rpkm.quinolone.beta, method="Wald")
rpkm.cef.com <- rpkm.cef.com[-c(1, 2),]
rpkm.quinolone1 <- cbind(rpkm.quinolone.beta2, rpkm.cef.com)
rpkm.quinolone1 <- rpkm.quinolone1[-1,]
rpkm.quinolone.beta <- rpkm.quinolone1 %>% mutate(effect = round(exp(Estimate), digits =2), lower=round(exp(`2.5 %`), digits =2), upper=round(exp(`97.5 %`), digits =2)) %>%
  mutate(exposure="quinolone", ARG_class = "Beta-Lactam", model= "Abundance", p.value = `Pr(>|z|)`) %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)

num.quinolone.betalactam <- glmmadmb(num.betalactam.ARG ~ quinolone.exp + (1|study_id), data = test9, family = "nbinom1")
summary(num.quinolone.betalactam)
num.quinolone.beta <- tidy(num.quinolone.betalactam,
                          effects = "fixed",
                          component = "cond",
                          scales=NULL,
                          ran_prefix = NULL,
                          conf.int = TRUE,
                          conf.level = 0.95,
                          conf.method = "Wald")
num.quinolone.beta <- num.quinolone.beta %>%
  mutate(exposure = "quinolone", ARG_class = "Beta-Lactam", model = "Number", effect=round(exp(estimate),2),
         lower=round(exp(conf.low), 2), upper = round(exp(conf.high), 2)) %>%
  filter(term != "(Intercept)") %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)

rpkm.quinolone.vanc <- glmer(vanc.abund ~ quinolone.exp + (1 |study_id), data  = test9, family=Gamma(link="log"))
summary(rpkm.quinolone.vanc)
rpkm.quinolone.vanc2 <- summary(rpkm.quinolone.vanc)
rpkm.quinolone.vanc2 <- as.data.frame(rpkm.quinolone.vanc2$coefficients)
rpkm.cef.com <- confint(rpkm.quinolone.vanc, method="Wald")
rpkm.cef.com <- rpkm.cef.com[-c(1, 2),]
rpkm.quinolone1 <- cbind(rpkm.quinolone.vanc2, rpkm.cef.com)
rpkm.quinolone1 <- rpkm.quinolone1[-1,]
rpkm.quinolone.vanc <- rpkm.quinolone1 %>% mutate(effect = round(exp(Estimate), digits =2), lower=round(exp(`2.5 %`), digits =2), upper=round(exp(`97.5 %`), digits =2)) %>%
  mutate(exposure="quinolone", ARG_class = "Glycopeptide", model= "Abundance", p.value = `Pr(>|z|)`) %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)

num.quinolone.vanc <- glmmadmb(num.vanc.ARG ~ quinolone.exp + (1|study_id), data = test9, family = "nbinom1")
summary(num.quinolone.vanc)
num.quinolone.vanc <- tidy(num.quinolone.vanc,
                          effects = "fixed",
                          component = "cond",
                          scales=NULL,
                          ran_prefix = NULL,
                          conf.int = TRUE,
                          conf.level = 0.95,
                          conf.method = "Wald")
num.quinolone.vanc <- num.quinolone.vanc %>%
  mutate(exposure = "quinolone", ARG_class = "Glycopeptide", model = "Number", effect=round(exp(estimate),2),
         lower=round(exp(conf.low), 2), upper = round(exp(conf.high), 2)) %>%
  filter(term != "(Intercept)") %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)


rpkm.quinolone.carb <- glmer(carb.abund ~ quinolone.exp + (1 |study_id), data  = test9, family=Gamma(link="log"))
summary(rpkm.quinolone.carb)
rpkm.quinolone.carb2 <- summary(rpkm.quinolone.carb)
rpkm.quinolone.carb2 <- as.data.frame(rpkm.quinolone.carb2$coefficients)
rpkm.cef.com <- confint(rpkm.quinolone.carb, method="Wald")
rpkm.cef.com <- rpkm.cef.com[-c(1, 2),]
rpkm.quinolone1 <- cbind(rpkm.quinolone.carb2, rpkm.cef.com)
rpkm.quinolone1 <- rpkm.quinolone1[-1,]
rpkm.quinolone.carb <- rpkm.quinolone1 %>% mutate(effect = round(exp(Estimate), digits =2), lower=round(exp(`2.5 %`), digits =2), upper=round(exp(`97.5 %`), digits =2)) %>%
  mutate(exposure="quinolone", ARG_class = "Carbapenem", model= "Abundance", p.value = `Pr(>|z|)`) %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)

num.quinolone.carb <- glmmadmb(num.carbapenem.ARG ~ quinolone.exp + (1|study_id), data = test9, family = "nbinom1")
summary(num.quinolone.carb)
num.quinolone.carb <- tidy(num.quinolone.carb,
                          effects = "fixed",
                          component = "cond",
                          scales=NULL,
                          ran_prefix = NULL,
                          conf.int = TRUE,
                          conf.level = 0.95,
                          conf.method = "Wald")
num.quinolone.carb <- num.quinolone.carb %>%
  mutate(exposure = "quinolone", ARG_class = "Carbapenem", model = "Number", effect=round(exp(estimate),2),
         lower=round(exp(conf.low), 2), upper = round(exp(conf.high), 2)) %>%
  filter(term != "(Intercept)") %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)

rpkm.quinolone.fluoro <- glmer(fluoro.abund ~ quinolone.exp + (1 |study_id), data  = test9, family=Gamma(link="log"))
summary(rpkm.quinolone.fluoro)
rpkm.quinolone.fluoro2 <- summary(rpkm.quinolone.fluoro)
rpkm.quinolone.fluoro2 <- as.data.frame(rpkm.quinolone.fluoro2$coefficients)
rpkm.cef.com <- confint(rpkm.quinolone.fluoro, method="Wald")
rpkm.cef.com <- rpkm.cef.com[-c(1, 2),]
rpkm.quinolone1 <- cbind(rpkm.quinolone.fluoro2, rpkm.cef.com)
rpkm.quinolone1 <- rpkm.quinolone1[-1,]
rpkm.quinolone.fluoro <- rpkm.quinolone1 %>% mutate(effect = round(exp(Estimate), digits =2), lower=round(exp(`2.5 %`), digits =2), upper=round(exp(`97.5 %`), digits =2)) %>%
  mutate(exposure="quinolone", ARG_class = "Fluoroquinolone", model= "Abundance", p.value = `Pr(>|z|)`) %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)

num.quinolone.fluoro <- glmmadmb(num.quinolone.ARG ~ quinolone.exp + (1|study_id), data = test9, family = "nbinom1")
summary(num.quinolone.fluoro)
num.quinolone.fluoro <- tidy(num.quinolone.fluoro,
                            effects = "fixed",
                            component = "cond",
                            scales=NULL,
                            ran_prefix = NULL,
                            conf.int = TRUE,
                            conf.level = 0.95,
                            conf.method = "Wald")
num.quinolone.fluoro <- num.quinolone.fluoro %>%
  mutate(exposure = "quinolone", ARG_class = "Fluoroquinolone", model = "Number", effect=round(exp(estimate),2),
         lower=round(exp(conf.low), 2), upper = round(exp(conf.high), 2)) %>%
  filter(term != "(Intercept)") %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)


rpkm.quinolone.amino <- glmer(amino.abund ~ quinolone.exp + (1 |study_id), data  = test9, family=Gamma(link="log"))
summary(rpkm.quinolone.amino)
rpkm.quinolone.amino2 <- summary(rpkm.quinolone.amino)
rpkm.quinolone.amino2 <- as.data.frame(rpkm.quinolone.amino2$coefficients)
rpkm.cef.com <- confint(rpkm.quinolone.amino, method="Wald")
rpkm.cef.com <- rpkm.cef.com[-c(1, 2),]
rpkm.quinolone1 <- cbind(rpkm.quinolone.amino2, rpkm.cef.com)
rpkm.quinolone1 <- rpkm.quinolone1[-1,]
rpkm.quinolone.amino <- rpkm.quinolone1 %>% mutate(effect = round(exp(Estimate), digits =2), lower=round(exp(`2.5 %`), digits =2), upper=round(exp(`97.5 %`), digits =2)) %>%
  mutate(exposure="quinolone", ARG_class = "Aminoglycoside", model= "Abundance", p.value = `Pr(>|z|)`) %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)

num.quinolone.amino <- glmmadmb(num.amino.ARG ~ quinolone.exp + (1|study_id), data = test9, family = "nbinom1")
summary(num.quinolone.amino)
num.quinolone.amino <- tidy(num.quinolone.amino,
                           effects = "fixed",
                           component = "cond",
                           scales=NULL,
                           ran_prefix = NULL,
                           conf.int = TRUE,
                           conf.level = 0.95,
                           conf.method = "Wald")
num.quinolone.amino <- num.quinolone.amino %>%
  mutate(exposure = "quinolone", ARG_class = "Aminoglycoside", model = "Number", effect=round(exp(estimate),2),
         lower=round(exp(conf.low), 2), upper = round(exp(conf.high), 2)) %>%
  filter(term != "(Intercept)") %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)


rpkm.quinolone.clinda <- glmer(clinda.abund ~ quinolone.exp + (1 |study_id), data  = test9, family=Gamma(link="log"))
summary(rpkm.quinolone.clinda)
rpkm.quinolone.clinda2 <- summary(rpkm.quinolone.clinda)
rpkm.quinolone.clinda2 <- as.data.frame(rpkm.quinolone.clinda2$coefficients)
rpkm.cef.com <- confint(rpkm.quinolone.clinda, method="Wald")
rpkm.cef.com <- rpkm.cef.com[-c(1, 2),]
rpkm.quinolone1 <- cbind(rpkm.quinolone.clinda2, rpkm.cef.com)
rpkm.quinolone1 <- rpkm.quinolone1[-1,]
rpkm.quinolone.clinda <- rpkm.quinolone1 %>% mutate(effect = round(exp(Estimate), digits =2), lower=round(exp(`2.5 %`), digits =2), upper=round(exp(`97.5 %`), digits =2)) %>%
  mutate(exposure="quinolone", ARG_class = "Lincosamide", model= "Abundance", p.value = `Pr(>|z|)`) %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)

num.quinolone.clinda <- glmmadmb(num.clinda.ARG ~ quinolone.exp + (1|study_id), data = test9, family = "nbinom1")
summary(num.quinolone.clinda)
num.quinolone.clinda <- tidy(num.quinolone.clinda,
                            effects = "fixed",
                            component = "cond",
                            scales=NULL,
                            ran_prefix = NULL,
                            conf.int = TRUE,
                            conf.level = 0.95,
                            conf.method = "Wald")
num.quinolone.clinda <- num.quinolone.clinda %>%
  mutate(exposure = "quinolone", ARG_class = "Lincosamide", model = "Number", effect=round(exp(estimate),2),
         lower=round(exp(conf.low), 2), upper = round(exp(conf.high), 2)) %>%
  filter(term != "(Intercept)") %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)



rpkm.quinolone.macrolide <- glmer(macrolide.abund ~ quinolone.exp + (1 |study_id), data  = test9, family=Gamma(link="log"))
summary(rpkm.quinolone.macrolide)
rpkm.quinolone.macrolide2 <- summary(rpkm.quinolone.macrolide)
rpkm.quinolone.macrolide2 <- as.data.frame(rpkm.quinolone.macrolide2$coefficients)
rpkm.cef.com <- confint(rpkm.quinolone.macrolide, method="Wald")
rpkm.cef.com <- rpkm.cef.com[-c(1, 2),]
rpkm.quinolone1 <- cbind(rpkm.quinolone.macrolide2, rpkm.cef.com)
rpkm.quinolone1 <- rpkm.quinolone1[-1,]
rpkm.quinolone.macrolide <- rpkm.quinolone1 %>% mutate(effect = round(exp(Estimate), digits =2), lower=round(exp(`2.5 %`), digits =2), upper=round(exp(`97.5 %`), digits =2)) %>%
  mutate(exposure="quinolone", ARG_class = "Macrolide", model= "Abundance", p.value = `Pr(>|z|)`) %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)

num.quinolone.macrolide <- glmmadmb(num.macrolide.ARG ~ quinolone.exp + (1|study_id), data = test9, family = "nbinom1")
summary(num.quinolone.macrolide)
num.quinolone.macrolide <- tidy(num.quinolone.macrolide,
                               effects = "fixed",
                               component = "cond",
                               scales=NULL,
                               ran_prefix = NULL,
                               conf.int = TRUE,
                               conf.level = 0.95,
                               conf.method = "Wald")
num.quinolone.macrolide <- num.quinolone.macrolide %>%
  mutate(exposure = "quinolone", ARG_class = "Macrolide", model = "Number", effect=round(exp(estimate),2),
         lower=round(exp(conf.low), 2), upper = round(exp(conf.high), 2)) %>%
  filter(term != "(Intercept)") %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)


rpkm.quinolone.trimeth <- glmer(trimeth.abund ~ quinolone.exp + (1 |study_id), data  = test9, family=Gamma(link="log"))
summary(rpkm.quinolone.trimeth)
rpkm.quinolone.trimeth2 <- summary(rpkm.quinolone.trimeth)
rpkm.quinolone.trimeth2 <- as.data.frame(rpkm.quinolone.trimeth2$coefficients)
rpkm.cef.com <- confint(rpkm.quinolone.trimeth, method="Wald")
rpkm.cef.com <- rpkm.cef.com[-c(1, 2),]
rpkm.quinolone1 <- cbind(rpkm.quinolone.trimeth2, rpkm.cef.com)
rpkm.quinolone1 <- rpkm.quinolone1[-1,]
rpkm.quinolone.trimeth <- rpkm.quinolone1 %>% mutate(effect = round(exp(Estimate), digits =2), lower=round(exp(`2.5 %`), digits =2), upper=round(exp(`97.5 %`), digits =2)) %>%
  mutate(exposure="quinolone", ARG_class = "Diaminopyrimidine", model= "Abundance", p.value = `Pr(>|z|)`) %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)

num.quinolone.trimeth <- glmmadmb(num.trimeth.ARG ~ quinolone.exp + (1|study_id), data = test9, family = "nbinom2")
summary(num.quinolone.trimeth)
num.quinolone.trimeth <- tidy(num.quinolone.trimeth,
                             effects = "fixed",
                             component = "cond",
                             scales=NULL,
                             ran_prefix = NULL,
                             conf.int = TRUE,
                             conf.level = 0.95,
                             conf.method = "Wald")
num.quinolone.trimeth <- num.quinolone.trimeth %>%
  mutate(exposure = "quinolone", ARG_class = "Diaminopyrimidine", model = "Number", effect=round(exp(estimate),2),
         lower=round(exp(conf.low), 2), upper = round(exp(conf.high), 2)) %>%
  filter(term != "(Intercept)") %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)


rpkm.quinolone.tetra <- glmer(tetra.abund ~ quinolone.exp + (1 |study_id), data  = test9, family=Gamma(link="log"))
summary(rpkm.quinolone.tetra)
rpkm.quinolone.tetra2 <- summary(rpkm.quinolone.tetra)
rpkm.quinolone.tetra2 <- as.data.frame(rpkm.quinolone.tetra2$coefficients)
rpkm.cef.com <- confint(rpkm.quinolone.tetra, method="Wald")
rpkm.cef.com <- rpkm.cef.com[-c(1, 2),]
rpkm.quinolone1 <- cbind(rpkm.quinolone.tetra2, rpkm.cef.com)
rpkm.quinolone1 <- rpkm.quinolone1[-1,]
rpkm.quinolone.tetra <- rpkm.quinolone1 %>% mutate(effect = round(exp(Estimate), digits =2), lower=round(exp(`2.5 %`), digits =2), upper=round(exp(`97.5 %`), digits =2)) %>%
  mutate(exposure="quinolone", ARG_class = "Tetracycline", model= "Abundance", p.value = `Pr(>|z|)`) %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)

num.quinolone.tetra <- glmmadmb(num.tetra.ARG ~ quinolone.exp + (1|study_id), data = test9, family = "nbinom1")
summary(num.quinolone.tetra)
num.quinolone.tetra <- tidy(num.quinolone.tetra,
                           effects = "fixed",
                           component = "cond",
                           scales=NULL,
                           ran_prefix = NULL,
                           conf.int = TRUE,
                           conf.level = 0.95,
                           conf.method = "Wald")
num.quinolone.tetra <- num.quinolone.tetra %>%
  mutate(exposure = "quinolone", ARG_class = "Tetracycline", model = "Number", effect=round(exp(estimate),2),
         lower=round(exp(conf.low), 2), upper = round(exp(conf.high), 2)) %>%
  filter(term != "(Intercept)") %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)


rpkm.quinolone.flagyl <- glmer(flagyl.abund ~ quinolone.exp + (1 |study_id), data  = test9, family=Gamma(link="log"))
summary(rpkm.quinolone.flagyl)
rpkm.quinolone.flagyl2 <- summary(rpkm.quinolone.flagyl)
rpkm.quinolone.flagyl2 <- as.data.frame(rpkm.quinolone.flagyl2$coefficients)
rpkm.cef.com <- confint(rpkm.quinolone.flagyl, method="Wald")
rpkm.cef.com <- rpkm.cef.com[-c(1, 2),]
rpkm.quinolone1 <- cbind(rpkm.quinolone.flagyl2, rpkm.cef.com)
rpkm.quinolone1 <- rpkm.quinolone1[-1,]
rpkm.quinolone.flagyl <- rpkm.quinolone1 %>% mutate(effect = round(exp(Estimate), digits =2), lower=round(exp(`2.5 %`), digits =2), upper=round(exp(`97.5 %`), digits =2)) %>%
  mutate(exposure="quinolone", ARG_class = "Nitroimidazole", model= "Abundance", p.value = `Pr(>|z|)`) %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)

num.quinolone.flagyl <- glmmadmb(num.flagyl.ARG ~ quinolone.exp + (1|study_id), data = test9, family = "nbinom1")
summary(num.quinolone.flagyl)
num.quinolone.flagyl <- tidy(num.quinolone.flagyl,
                            effects = "fixed",
                            component = "cond",
                            scales=NULL,
                            ran_prefix = NULL,
                            conf.int = TRUE,
                            conf.level = 0.95,
                            conf.method = "Wald")
num.quinolone.flagyl <- num.quinolone.flagyl %>%
  mutate(exposure = "quinolone", ARG_class = "Nitroimidazole", model = "Number", effect=round(exp(estimate),2),
         lower=round(exp(conf.low), 2), upper = round(exp(conf.high), 2)) %>%
  filter(term != "(Intercept)") %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)

quinolone <- rbind(rpkm.quinolone.beta,
                  num.quinolone.beta,
                  rpkm.quinolone.carb,
                  num.quinolone.carb,
                  rpkm.quinolone.vanc,
                  num.quinolone.vanc,
                  rpkm.quinolone.flagyl,
                  num.quinolone.flagyl,
                  rpkm.quinolone.fluoro,
                  num.quinolone.fluoro,
                  rpkm.quinolone.trimeth,
                  num.quinolone.trimeth,
                  rpkm.quinolone.macrolide,
                  num.quinolone.macrolide,
                  rpkm.quinolone.tetra,
                  num.quinolone.tetra,
                  rpkm.quinolone.amino,
                  num.quinolone.amino,
                  rpkm.quinolone.clinda,
                  num.quinolone.clinda)
rownames(quinolone) <- NULL

quinolone <- quinolone %>%
  mutate(p.adj = p.adjust(p.value, method = "BH"))
quinolone$p.adj <- signif(quinolone$p.adj, 2)
quinolone$p.adj <- format(quinolone$p.adj, scientific = FALSE)

#trimethoprim exposure
rpkm.tmp_smx.beta <- glmer(beta.abund ~ tmp_smx.exp + (1 |study_id), data  = test9, family=Gamma(link="log"))
summary(rpkm.tmp_smx.beta)
rpkm.tmp_smx.beta2 <- summary(rpkm.tmp_smx.beta)
rpkm.tmp_smx.beta2 <- as.data.frame(rpkm.tmp_smx.beta2$coefficients)
rpkm.cef.com <- confint(rpkm.tmp_smx.beta, method="Wald")
rpkm.cef.com <- rpkm.cef.com[-c(1, 2),]
rpkm.tmp_smx1 <- cbind(rpkm.tmp_smx.beta2, rpkm.cef.com)
rpkm.tmp_smx1 <- rpkm.tmp_smx1[-1,]
rpkm.tmp_smx.beta <- rpkm.tmp_smx1 %>% mutate(effect = round(exp(Estimate), digits =2), lower=round(exp(`2.5 %`), digits =2), upper=round(exp(`97.5 %`), digits =2)) %>%
  mutate(exposure="TMP-SMX", ARG_class = "Beta-Lactam", model= "Abundance", p.value = `Pr(>|z|)`) %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)

num.tmp_smx.betalactam <- glmmadmb(num.betalactam.ARG ~ tmp_smx.exp + (1|study_id), data = test9, family = "nbinom1")
summary(num.tmp_smx.betalactam)
num.tmp_smx.beta <- tidy(num.tmp_smx.betalactam,
                          effects = "fixed",
                          component = "cond",
                          scales=NULL,
                          ran_prefix = NULL,
                          conf.int = TRUE,
                          conf.level = 0.95,
                          conf.method = "Wald")
num.tmp_smx.beta <- num.tmp_smx.beta %>%
  mutate(exposure = "TMP-SMX", ARG_class = "Beta-Lactam", model = "Number", effect=round(exp(estimate),2),
         lower=round(exp(conf.low), 2), upper = round(exp(conf.high), 2)) %>%
  filter(term != "(Intercept)") %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)

rpkm.tmp_smx.vanc <- glmer(vanc.abund ~ tmp_smx.exp + (1 |study_id), data  = test9, family=Gamma(link="log"))
summary(rpkm.tmp_smx.vanc)
rpkm.tmp_smx.vanc2 <- summary(rpkm.tmp_smx.vanc)
rpkm.tmp_smx.vanc2 <- as.data.frame(rpkm.tmp_smx.vanc2$coefficients)
rpkm.cef.com <- confint(rpkm.tmp_smx.vanc, method="Wald")
rpkm.cef.com <- rpkm.cef.com[-c(1, 2),]
rpkm.tmp_smx1 <- cbind(rpkm.tmp_smx.vanc2, rpkm.cef.com)
rpkm.tmp_smx1 <- rpkm.tmp_smx1[-1,]
rpkm.tmp_smx.vanc <- rpkm.tmp_smx1 %>% mutate(effect = round(exp(Estimate), digits =2), lower=round(exp(`2.5 %`), digits =2), upper=round(exp(`97.5 %`), digits =2)) %>%
  mutate(exposure="TMP-SMX", ARG_class = "Glycopeptide", model= "Abundance", p.value = `Pr(>|z|)`) %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)

num.tmp_smx.vanc <- glmmadmb(num.vanc.ARG ~ tmp_smx.exp + (1|study_id), data = test9, family = "nbinom1")
summary(num.tmp_smx.vanc)
num.tmp_smx.vanc <- tidy(num.tmp_smx.vanc,
                          effects = "fixed",
                          component = "cond",
                          scales=NULL,
                          ran_prefix = NULL,
                          conf.int = TRUE,
                          conf.level = 0.95,
                          conf.method = "Wald")
num.tmp_smx.vanc <- num.tmp_smx.vanc %>%
  mutate(exposure = "TMP-SMX", ARG_class = "Glycopeptide", model = "Number", effect=round(exp(estimate),2),
         lower=round(exp(conf.low), 2), upper = round(exp(conf.high), 2)) %>%
  filter(term != "(Intercept)") %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)


rpkm.tmp_smx.carb <- glmer(carb.abund ~ tmp_smx.exp + (1 |study_id), data  = test9, family=Gamma(link="log"))
summary(rpkm.tmp_smx.carb)
rpkm.tmp_smx.carb2 <- summary(rpkm.tmp_smx.carb)
rpkm.tmp_smx.carb2 <- as.data.frame(rpkm.tmp_smx.carb2$coefficients)
rpkm.cef.com <- confint(rpkm.tmp_smx.carb, method="Wald")
rpkm.cef.com <- rpkm.cef.com[-c(1, 2),]
rpkm.tmp_smx1 <- cbind(rpkm.tmp_smx.carb2, rpkm.cef.com)
rpkm.tmp_smx1 <- rpkm.tmp_smx1[-1,]
rpkm.tmp_smx.carb <- rpkm.tmp_smx1 %>% mutate(effect = round(exp(Estimate), digits =2), lower=round(exp(`2.5 %`), digits =2), upper=round(exp(`97.5 %`), digits =2)) %>%
  mutate(exposure="TMP-SMX", ARG_class = "Carbapenem", model= "Abundance", p.value = `Pr(>|z|)`) %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)

num.tmp_smx.carb <- glmmadmb(num.carbapenem.ARG ~ tmp_smx.exp + (1|study_id), data = test9, family = "nbinom1")
summary(num.tmp_smx.carb)
num.tmp_smx.carb <- tidy(num.tmp_smx.carb,
                          effects = "fixed",
                          component = "cond",
                          scales=NULL,
                          ran_prefix = NULL,
                          conf.int = TRUE,
                          conf.level = 0.95,
                          conf.method = "Wald")
num.tmp_smx.carb <- num.tmp_smx.carb %>%
  mutate(exposure = "TMP-SMX", ARG_class = "Carbapenem", model = "Number", effect=round(exp(estimate),2),
         lower=round(exp(conf.low), 2), upper = round(exp(conf.high), 2)) %>%
  filter(term != "(Intercept)") %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)

rpkm.tmp_smx.fluoro <- glmer(fluoro.abund ~ tmp_smx.exp + (1 |study_id), data  = test9, family=Gamma(link="log"))
summary(rpkm.tmp_smx.fluoro)
rpkm.tmp_smx.fluoro2 <- summary(rpkm.tmp_smx.fluoro)
rpkm.tmp_smx.fluoro2 <- as.data.frame(rpkm.tmp_smx.fluoro2$coefficients)
rpkm.cef.com <- confint(rpkm.tmp_smx.fluoro, method="Wald")
rpkm.cef.com <- rpkm.cef.com[-c(1, 2),]
rpkm.tmp_smx1 <- cbind(rpkm.tmp_smx.fluoro2, rpkm.cef.com)
rpkm.tmp_smx1 <- rpkm.tmp_smx1[-1,]
rpkm.tmp_smx.fluoro <- rpkm.tmp_smx1 %>% mutate(effect = round(exp(Estimate), digits =2), lower=round(exp(`2.5 %`), digits =2), upper=round(exp(`97.5 %`), digits =2)) %>%
  mutate(exposure="TMP-SMX", ARG_class = "Fluoroquinolone", model= "Abundance", p.value = `Pr(>|z|)`) %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)

num.tmp_smx.fluoro <- glmmadmb(num.quinolone.ARG ~ tmp_smx.exp + (1|study_id), data = test9, family = "nbinom1")
summary(num.tmp_smx.fluoro)
num.tmp_smx.fluoro <- tidy(num.tmp_smx.fluoro,
                            effects = "fixed",
                            component = "cond",
                            scales=NULL,
                            ran_prefix = NULL,
                            conf.int = TRUE,
                            conf.level = 0.95,
                            conf.method = "Wald")
num.tmp_smx.fluoro <- num.tmp_smx.fluoro %>%
  mutate(exposure = "TMP-SMX", ARG_class = "Fluoroquinolone", model = "Number", effect=round(exp(estimate),2),
         lower=round(exp(conf.low), 2), upper = round(exp(conf.high), 2)) %>%
  filter(term != "(Intercept)") %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)


rpkm.tmp_smx.amino <- glmer(amino.abund ~ tmp_smx.exp + (1 |study_id), data  = test9, family=Gamma(link="log"))
summary(rpkm.tmp_smx.amino)
rpkm.tmp_smx.amino2 <- summary(rpkm.tmp_smx.amino)
rpkm.tmp_smx.amino2 <- as.data.frame(rpkm.tmp_smx.amino2$coefficients)
rpkm.cef.com <- confint(rpkm.tmp_smx.amino, method="Wald")
rpkm.cef.com <- rpkm.cef.com[-c(1, 2),]
rpkm.tmp_smx1 <- cbind(rpkm.tmp_smx.amino2, rpkm.cef.com)
rpkm.tmp_smx1 <- rpkm.tmp_smx1[-1,]
rpkm.tmp_smx.amino <- rpkm.tmp_smx1 %>% mutate(effect = round(exp(Estimate), digits =2), lower=round(exp(`2.5 %`), digits =2), upper=round(exp(`97.5 %`), digits =2)) %>%
  mutate(exposure="TMP-SMX", ARG_class = "Aminoglycoside", model= "Abundance", p.value = `Pr(>|z|)`) %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)

num.tmp_smx.amino <- glmmadmb(num.amino.ARG ~ tmp_smx.exp + (1|study_id), data = test9, family = "nbinom1")
summary(num.tmp_smx.amino)
num.tmp_smx.amino <- tidy(num.tmp_smx.amino,
                           effects = "fixed",
                           component = "cond",
                           scales=NULL,
                           ran_prefix = NULL,
                           conf.int = TRUE,
                           conf.level = 0.95,
                           conf.method = "Wald")
num.tmp_smx.amino <- num.tmp_smx.amino %>%
  mutate(exposure = "TMP-SMX", ARG_class = "Aminoglycoside", model = "Number", effect=round(exp(estimate),2),
         lower=round(exp(conf.low), 2), upper = round(exp(conf.high), 2)) %>%
  filter(term != "(Intercept)") %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)


rpkm.tmp_smx.clinda <- glmer(clinda.abund ~ tmp_smx.exp + (1 |study_id), data  = test9, family=Gamma(link="log"))
summary(rpkm.tmp_smx.clinda)
rpkm.tmp_smx.clinda2 <- summary(rpkm.tmp_smx.clinda)
rpkm.tmp_smx.clinda2 <- as.data.frame(rpkm.tmp_smx.clinda2$coefficients)
rpkm.cef.com <- confint(rpkm.tmp_smx.clinda, method="Wald")
rpkm.cef.com <- rpkm.cef.com[-c(1, 2),]
rpkm.tmp_smx1 <- cbind(rpkm.tmp_smx.clinda2, rpkm.cef.com)
rpkm.tmp_smx1 <- rpkm.tmp_smx1[-1,]
rpkm.tmp_smx.clinda <- rpkm.tmp_smx1 %>% mutate(effect = round(exp(Estimate), digits =2), lower=round(exp(`2.5 %`), digits =2), upper=round(exp(`97.5 %`), digits =2)) %>%
  mutate(exposure="TMP-SMX", ARG_class = "Lincosamide", model= "Abundance", p.value = `Pr(>|z|)`) %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)

num.tmp_smx.clinda <- glmmadmb(num.clinda.ARG ~ tmp_smx.exp + (1|study_id), data = test9, family = "nbinom1")
summary(num.tmp_smx.clinda)
num.tmp_smx.clinda <- tidy(num.tmp_smx.clinda,
                            effects = "fixed",
                            component = "cond",
                            scales=NULL,
                            ran_prefix = NULL,
                            conf.int = TRUE,
                            conf.level = 0.95,
                            conf.method = "Wald")
num.tmp_smx.clinda <- num.tmp_smx.clinda %>%
  mutate(exposure = "TMP-SMX", ARG_class = "Lincosamide", model = "Number", effect=round(exp(estimate),2),
         lower=round(exp(conf.low), 2), upper = round(exp(conf.high), 2)) %>%
  filter(term != "(Intercept)") %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)



rpkm.tmp_smx.macrolide <- glmer(macrolide.abund ~ tmp_smx.exp + (1 |study_id), data  = test9, family=Gamma(link="log"))
summary(rpkm.tmp_smx.macrolide)
rpkm.tmp_smx.macrolide2 <- summary(rpkm.tmp_smx.macrolide)
rpkm.tmp_smx.macrolide2 <- as.data.frame(rpkm.tmp_smx.macrolide2$coefficients)
rpkm.cef.com <- confint(rpkm.tmp_smx.macrolide, method="Wald")
rpkm.cef.com <- rpkm.cef.com[-c(1, 2),]
rpkm.tmp_smx1 <- cbind(rpkm.tmp_smx.macrolide2, rpkm.cef.com)
rpkm.tmp_smx1 <- rpkm.tmp_smx1[-1,]
rpkm.tmp_smx.macrolide <- rpkm.tmp_smx1 %>% mutate(effect = round(exp(Estimate), digits =2), lower=round(exp(`2.5 %`), digits =2), upper=round(exp(`97.5 %`), digits =2)) %>%
  mutate(exposure="TMP-SMX", ARG_class = "Macrolide", model= "Abundance", p.value = `Pr(>|z|)`) %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)

num.tmp_smx.macrolide <- glmmadmb(num.macrolide.ARG ~ tmp_smx.exp + (1|study_id), data = test9, family = "nbinom1")
summary(num.tmp_smx.macrolide)
num.tmp_smx.macrolide <- tidy(num.tmp_smx.macrolide,
                               effects = "fixed",
                               component = "cond",
                               scales=NULL,
                               ran_prefix = NULL,
                               conf.int = TRUE,
                               conf.level = 0.95,
                               conf.method = "Wald")
num.tmp_smx.macrolide <- num.tmp_smx.macrolide %>%
  mutate(exposure = "TMP-SMX", ARG_class = "Macrolide", model = "Number", effect=round(exp(estimate),2),
         lower=round(exp(conf.low), 2), upper = round(exp(conf.high), 2)) %>%
  filter(term != "(Intercept)") %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)


rpkm.tmp_smx.trimeth <- glmer(trimeth.abund ~ tmp_smx.exp + (1 |study_id), data  = test9, family=Gamma(link="log"))
summary(rpkm.tmp_smx.trimeth)
rpkm.tmp_smx.trimeth2 <- summary(rpkm.tmp_smx.trimeth)
rpkm.tmp_smx.trimeth2 <- as.data.frame(rpkm.tmp_smx.trimeth2$coefficients)
rpkm.cef.com <- confint(rpkm.tmp_smx.trimeth, method="Wald")
rpkm.cef.com <- rpkm.cef.com[-c(1, 2),]
rpkm.tmp_smx1 <- cbind(rpkm.tmp_smx.trimeth2, rpkm.cef.com)
rpkm.tmp_smx1 <- rpkm.tmp_smx1[-1,]
rpkm.tmp_smx.trimeth <- rpkm.tmp_smx1 %>% mutate(effect = round(exp(Estimate), digits =2), lower=round(exp(`2.5 %`), digits =2), upper=round(exp(`97.5 %`), digits =2)) %>%
  mutate(exposure="TMP-SMX", ARG_class = "Diaminopyrimidine", model= "Abundance", p.value = `Pr(>|z|)`) %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)

num.tmp_smx.trimeth <- glmmadmb(num.trimeth.ARG ~ tmp_smx.exp + (1|study_id), data = test9, family = "nbinom2")
summary(num.tmp_smx.trimeth)
num.tmp_smx.trimeth <- tidy(num.tmp_smx.trimeth,
                             effects = "fixed",
                             component = "cond",
                             scales=NULL,
                             ran_prefix = NULL,
                             conf.int = TRUE,
                             conf.level = 0.95,
                             conf.method = "Wald")
num.tmp_smx.trimeth <- num.tmp_smx.trimeth %>%
  mutate(exposure = "TMP-SMX", ARG_class = "Diaminopyrimidine", model = "Number", effect=round(exp(estimate),2),
         lower=round(exp(conf.low), 2), upper = round(exp(conf.high), 2)) %>%
  filter(term != "(Intercept)") %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)


rpkm.tmp_smx.tetra <- glmer(tetra.abund ~ tmp_smx.exp + (1 |study_id), data  = test9, family=Gamma(link="log"))
summary(rpkm.tmp_smx.tetra)
rpkm.tmp_smx.tetra2 <- summary(rpkm.tmp_smx.tetra)
rpkm.tmp_smx.tetra2 <- as.data.frame(rpkm.tmp_smx.tetra2$coefficients)
rpkm.cef.com <- confint(rpkm.tmp_smx.tetra, method="Wald")
rpkm.cef.com <- rpkm.cef.com[-c(1, 2),]
rpkm.tmp_smx1 <- cbind(rpkm.tmp_smx.tetra2, rpkm.cef.com)
rpkm.tmp_smx1 <- rpkm.tmp_smx1[-1,]
rpkm.tmp_smx.tetra <- rpkm.tmp_smx1 %>% mutate(effect = round(exp(Estimate), digits =2), lower=round(exp(`2.5 %`), digits =2), upper=round(exp(`97.5 %`), digits =2)) %>%
  mutate(exposure="TMP-SMX", ARG_class = "Tetracycline", model= "Abundance", p.value = `Pr(>|z|)`) %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)

num.tmp_smx.tetra <- glmmadmb(num.tetra.ARG ~ tmp_smx.exp + (1|study_id), data = test9, family = "nbinom1")
summary(num.tmp_smx.tetra)
num.tmp_smx.tetra <- tidy(num.tmp_smx.tetra,
                           effects = "fixed",
                           component = "cond",
                           scales=NULL,
                           ran_prefix = NULL,
                           conf.int = TRUE,
                           conf.level = 0.95,
                           conf.method = "Wald")
num.tmp_smx.tetra <- num.tmp_smx.tetra %>%
  mutate(exposure = "TMP-SMX", ARG_class = "Tetracycline", model = "Number", effect=round(exp(estimate),2),
         lower=round(exp(conf.low), 2), upper = round(exp(conf.high), 2)) %>%
  filter(term != "(Intercept)") %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)


rpkm.tmp_smx.flagyl <- glmer(flagyl.abund ~ tmp_smx.exp + (1 |study_id), data  = test9, family=Gamma(link="log"))
summary(rpkm.tmp_smx.flagyl)
rpkm.tmp_smx.flagyl2 <- summary(rpkm.tmp_smx.flagyl)
rpkm.tmp_smx.flagyl2 <- as.data.frame(rpkm.tmp_smx.flagyl2$coefficients)
rpkm.cef.com <- confint(rpkm.tmp_smx.flagyl, method="Wald")
rpkm.cef.com <- rpkm.cef.com[-c(1, 2),]
rpkm.tmp_smx1 <- cbind(rpkm.tmp_smx.flagyl2, rpkm.cef.com)
rpkm.tmp_smx1 <- rpkm.tmp_smx1[-1,]
rpkm.tmp_smx.flagyl <- rpkm.tmp_smx1 %>% mutate(effect = round(exp(Estimate), digits =2), lower=round(exp(`2.5 %`), digits =2), upper=round(exp(`97.5 %`), digits =2)) %>%
  mutate(exposure="TMP-SMX", ARG_class = "Nitroimidazole", model= "Abundance", p.value = `Pr(>|z|)`) %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)

num.tmp_smx.flagyl <- glmmadmb(num.flagyl.ARG ~ tmp_smx.exp + (1|study_id), data = test9, family = "nbinom1")
summary(num.tmp_smx.flagyl)
num.tmp_smx.flagyl <- tidy(num.tmp_smx.flagyl,
                            effects = "fixed",
                            component = "cond",
                            scales=NULL,
                            ran_prefix = NULL,
                            conf.int = TRUE,
                            conf.level = 0.95,
                            conf.method = "Wald")
num.tmp_smx.flagyl <- num.tmp_smx.flagyl %>%
  mutate(exposure = "TMP-SMX", ARG_class = "Nitroimidazole", model = "Number", effect=round(exp(estimate),2),
         lower=round(exp(conf.low), 2), upper = round(exp(conf.high), 2)) %>%
  filter(term != "(Intercept)") %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)

tmp_smx <- rbind(rpkm.tmp_smx.beta,
                  num.tmp_smx.beta,
                  rpkm.tmp_smx.carb,
                  num.tmp_smx.carb,
                  rpkm.tmp_smx.vanc,
                  num.tmp_smx.vanc,
                  rpkm.tmp_smx.flagyl,
                  num.tmp_smx.flagyl,
                  rpkm.tmp_smx.fluoro,
                  num.tmp_smx.fluoro,
                  rpkm.tmp_smx.trimeth,
                  num.tmp_smx.trimeth,
                  rpkm.tmp_smx.macrolide,
                  num.tmp_smx.macrolide,
                  rpkm.tmp_smx.tetra,
                  num.tmp_smx.tetra,
                  rpkm.tmp_smx.amino,
                  num.tmp_smx.amino,
                  rpkm.tmp_smx.clinda,
                  num.tmp_smx.clinda)
rownames(tmp_smx) <- NULL

tmp_smx <- tmp_smx %>%
  mutate(p.adj = p.adjust(p.value, method = "BH"))
tmp_smx$p.adj <- signif(tmp_smx$p.adj, 2)
tmp_smx$p.adj <- format(tmp_smx$p.adj, scientific = FALSE)

#macrolide exposure
rpkm.macrolide.beta <- glmer(beta.abund ~ macrolide.exp + (1 |study_id), data  = test9, family=Gamma(link="log"))
summary(rpkm.macrolide.beta)
rpkm.macrolide.beta2 <- summary(rpkm.macrolide.beta)
rpkm.macrolide.beta2 <- as.data.frame(rpkm.macrolide.beta2$coefficients)
rpkm.cef.com <- confint(rpkm.macrolide.beta, method="Wald")
rpkm.cef.com <- rpkm.cef.com[-c(1, 2),]
rpkm.macrolide1 <- cbind(rpkm.macrolide.beta2, rpkm.cef.com)
rpkm.macrolide1 <- rpkm.macrolide1[-1,]
rpkm.macrolide.beta <- rpkm.macrolide1 %>% mutate(effect = round(exp(Estimate), digits =2), lower=round(exp(`2.5 %`), digits =2), upper=round(exp(`97.5 %`), digits =2)) %>%
  mutate(exposure="macrolide", ARG_class = "Beta-Lactam", model= "Abundance", p.value = `Pr(>|z|)`) %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)

num.macrolide.betalactam <- glmmadmb(num.betalactam.ARG ~ macrolide.exp + (1|study_id), data = test9, family = "nbinom1")
summary(num.macrolide.betalactam)
num.macrolide.beta <- tidy(num.macrolide.betalactam,
                          effects = "fixed",
                          component = "cond",
                          scales=NULL,
                          ran_prefix = NULL,
                          conf.int = TRUE,
                          conf.level = 0.95,
                          conf.method = "Wald")
num.macrolide.beta <- num.macrolide.beta %>%
  mutate(exposure = "macrolide", ARG_class = "Beta-Lactam", model = "Number", effect=round(exp(estimate),2),
         lower=round(exp(conf.low), 2), upper = round(exp(conf.high), 2)) %>%
  filter(term != "(Intercept)") %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)

rpkm.macrolide.vanc <- glmer(vanc.abund ~ macrolide.exp + (1 |study_id), data  = test9, family=Gamma(link="log"))
summary(rpkm.macrolide.vanc)
rpkm.macrolide.vanc2 <- summary(rpkm.macrolide.vanc)
rpkm.macrolide.vanc2 <- as.data.frame(rpkm.macrolide.vanc2$coefficients)
rpkm.cef.com <- confint(rpkm.macrolide.vanc, method="Wald")
rpkm.cef.com <- rpkm.cef.com[-c(1, 2),]
rpkm.macrolide1 <- cbind(rpkm.macrolide.vanc2, rpkm.cef.com)
rpkm.macrolide1 <- rpkm.macrolide1[-1,]
rpkm.macrolide.vanc <- rpkm.macrolide1 %>% mutate(effect = round(exp(Estimate), digits =2), lower=round(exp(`2.5 %`), digits =2), upper=round(exp(`97.5 %`), digits =2)) %>%
  mutate(exposure="macrolide", ARG_class = "Glycopeptide", model= "Abundance", p.value = `Pr(>|z|)`) %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)

num.macrolide.vanc <- glmmadmb(num.vanc.ARG ~ macrolide.exp + (1|study_id), data = test9, family = "nbinom1")
summary(num.macrolide.vanc)
num.macrolide.vanc <- tidy(num.macrolide.vanc,
                          effects = "fixed",
                          component = "cond",
                          scales=NULL,
                          ran_prefix = NULL,
                          conf.int = TRUE,
                          conf.level = 0.95,
                          conf.method = "Wald")
num.macrolide.vanc <- num.macrolide.vanc %>%
  mutate(exposure = "macrolide", ARG_class = "Glycopeptide", model = "Number", effect=round(exp(estimate),2),
         lower=round(exp(conf.low), 2), upper = round(exp(conf.high), 2)) %>%
  filter(term != "(Intercept)") %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)


rpkm.macrolide.carb <- glmer(carb.abund ~ macrolide.exp + (1 |study_id), data  = test9, family=Gamma(link="log"))
summary(rpkm.macrolide.carb)
rpkm.macrolide.carb2 <- summary(rpkm.macrolide.carb)
rpkm.macrolide.carb2 <- as.data.frame(rpkm.macrolide.carb2$coefficients)
rpkm.cef.com <- confint(rpkm.macrolide.carb, method="Wald")
rpkm.cef.com <- rpkm.cef.com[-c(1, 2),]
rpkm.macrolide1 <- cbind(rpkm.macrolide.carb2, rpkm.cef.com)
rpkm.macrolide1 <- rpkm.macrolide1[-1,]
rpkm.macrolide.carb <- rpkm.macrolide1 %>% mutate(effect = round(exp(Estimate), digits =2), lower=round(exp(`2.5 %`), digits =2), upper=round(exp(`97.5 %`), digits =2)) %>%
  mutate(exposure="macrolide", ARG_class = "Carbapenem", model= "Abundance", p.value = `Pr(>|z|)`) %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)

num.macrolide.carb <- glmmadmb(num.carbapenem.ARG ~ macrolide.exp + (1|study_id), data = test9, family = "nbinom1")
summary(num.macrolide.carb)
num.macrolide.carb <- tidy(num.macrolide.carb,
                          effects = "fixed",
                          component = "cond",
                          scales=NULL,
                          ran_prefix = NULL,
                          conf.int = TRUE,
                          conf.level = 0.95,
                          conf.method = "Wald")
num.macrolide.carb <- num.macrolide.carb %>%
  mutate(exposure = "macrolide", ARG_class = "Carbapenem", model = "Number", effect=round(exp(estimate),2),
         lower=round(exp(conf.low), 2), upper = round(exp(conf.high), 2)) %>%
  filter(term != "(Intercept)") %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)

rpkm.macrolide.fluoro <- glmer(fluoro.abund ~ macrolide.exp + (1 |study_id), data  = test9, family=Gamma(link="log"))
summary(rpkm.macrolide.fluoro)
rpkm.macrolide.fluoro2 <- summary(rpkm.macrolide.fluoro)
rpkm.macrolide.fluoro2 <- as.data.frame(rpkm.macrolide.fluoro2$coefficients)
rpkm.cef.com <- confint(rpkm.macrolide.fluoro, method="Wald")
rpkm.cef.com <- rpkm.cef.com[-c(1, 2),]
rpkm.macrolide1 <- cbind(rpkm.macrolide.fluoro2, rpkm.cef.com)
rpkm.macrolide1 <- rpkm.macrolide1[-1,]
rpkm.macrolide.fluoro <- rpkm.macrolide1 %>% mutate(effect = round(exp(Estimate), digits =2), lower=round(exp(`2.5 %`), digits =2), upper=round(exp(`97.5 %`), digits =2)) %>%
  mutate(exposure="macrolide", ARG_class = "Fluoroquinolone", model= "Abundance", p.value = `Pr(>|z|)`) %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)

num.macrolide.fluoro <- glmmadmb(num.quinolone.ARG ~ macrolide.exp + (1|study_id), data = test9, family = "nbinom1")
summary(num.macrolide.fluoro)
num.macrolide.fluoro <- tidy(num.macrolide.fluoro,
                            effects = "fixed",
                            component = "cond",
                            scales=NULL,
                            ran_prefix = NULL,
                            conf.int = TRUE,
                            conf.level = 0.95,
                            conf.method = "Wald")
num.macrolide.fluoro <- num.macrolide.fluoro %>%
  mutate(exposure = "macrolide", ARG_class = "Fluoroquinolone", model = "Number", effect=round(exp(estimate),2),
         lower=round(exp(conf.low), 2), upper = round(exp(conf.high), 2)) %>%
  filter(term != "(Intercept)") %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)


rpkm.macrolide.amino <- glmer(amino.abund ~ macrolide.exp + (1 |study_id), data  = test9, family=Gamma(link="log"))
summary(rpkm.macrolide.amino)
rpkm.macrolide.amino2 <- summary(rpkm.macrolide.amino)
rpkm.macrolide.amino2 <- as.data.frame(rpkm.macrolide.amino2$coefficients)
rpkm.cef.com <- confint(rpkm.macrolide.amino, method="Wald")
rpkm.cef.com <- rpkm.cef.com[-c(1, 2),]
rpkm.macrolide1 <- cbind(rpkm.macrolide.amino2, rpkm.cef.com)
rpkm.macrolide1 <- rpkm.macrolide1[-1,]
rpkm.macrolide.amino <- rpkm.macrolide1 %>% mutate(effect = round(exp(Estimate), digits =2), lower=round(exp(`2.5 %`), digits =2), upper=round(exp(`97.5 %`), digits =2)) %>%
  mutate(exposure="macrolide", ARG_class = "Aminoglycoside", model= "Abundance", p.value = `Pr(>|z|)`) %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)

num.macrolide.amino <- glmmadmb(num.amino.ARG ~ macrolide.exp + (1|study_id), data = test9, family = "nbinom1")
summary(num.macrolide.amino)
num.macrolide.amino <- tidy(num.macrolide.amino,
                           effects = "fixed",
                           component = "cond",
                           scales=NULL,
                           ran_prefix = NULL,
                           conf.int = TRUE,
                           conf.level = 0.95,
                           conf.method = "Wald")
num.macrolide.amino <- num.macrolide.amino %>%
  mutate(exposure = "macrolide", ARG_class = "Aminoglycoside", model = "Number", effect=round(exp(estimate),2),
         lower=round(exp(conf.low), 2), upper = round(exp(conf.high), 2)) %>%
  filter(term != "(Intercept)") %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)


rpkm.macrolide.clinda <- glmer(clinda.abund ~ macrolide.exp + (1 |study_id), data  = test9, family=Gamma(link="log"))
summary(rpkm.macrolide.clinda)
rpkm.macrolide.clinda2 <- summary(rpkm.macrolide.clinda)
rpkm.macrolide.clinda2 <- as.data.frame(rpkm.macrolide.clinda2$coefficients)
rpkm.cef.com <- confint(rpkm.macrolide.clinda, method="Wald")
rpkm.cef.com <- rpkm.cef.com[-c(1, 2),]
rpkm.macrolide1 <- cbind(rpkm.macrolide.clinda2, rpkm.cef.com)
rpkm.macrolide1 <- rpkm.macrolide1[-1,]
rpkm.macrolide.clinda <- rpkm.macrolide1 %>% mutate(effect = round(exp(Estimate), digits =2), lower=round(exp(`2.5 %`), digits =2), upper=round(exp(`97.5 %`), digits =2)) %>%
  mutate(exposure="macrolide", ARG_class = "Lincosamide", model= "Abundance", p.value = `Pr(>|z|)`) %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)

num.macrolide.clinda <- glmmadmb(num.clinda.ARG ~ macrolide.exp + (1|study_id), data = test9, family = "nbinom1")
summary(num.macrolide.clinda)
num.macrolide.clinda <- tidy(num.macrolide.clinda,
                            effects = "fixed",
                            component = "cond",
                            scales=NULL,
                            ran_prefix = NULL,
                            conf.int = TRUE,
                            conf.level = 0.95,
                            conf.method = "Wald")
num.macrolide.clinda <- num.macrolide.clinda %>%
  mutate(exposure = "macrolide", ARG_class = "Lincosamide", model = "Number", effect=round(exp(estimate),2),
         lower=round(exp(conf.low), 2), upper = round(exp(conf.high), 2)) %>%
  filter(term != "(Intercept)") %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)



rpkm.macrolide.macrolide <- glmer(macrolide.abund ~ macrolide.exp + (1 |study_id), data  = test9, family=Gamma(link="log"))
summary(rpkm.macrolide.macrolide)
rpkm.macrolide.macrolide2 <- summary(rpkm.macrolide.macrolide)
rpkm.macrolide.macrolide2 <- as.data.frame(rpkm.macrolide.macrolide2$coefficients)
rpkm.cef.com <- confint(rpkm.macrolide.macrolide, method="Wald")
rpkm.cef.com <- rpkm.cef.com[-c(1, 2),]
rpkm.macrolide1 <- cbind(rpkm.macrolide.macrolide2, rpkm.cef.com)
rpkm.macrolide1 <- rpkm.macrolide1[-1,]
rpkm.macrolide.macrolide <- rpkm.macrolide1 %>% mutate(effect = round(exp(Estimate), digits =2), lower=round(exp(`2.5 %`), digits =2), upper=round(exp(`97.5 %`), digits =2)) %>%
  mutate(exposure="macrolide", ARG_class = "Macrolide", model= "Abundance", p.value = `Pr(>|z|)`) %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)

num.macrolide.macrolide <- glmmadmb(num.macrolide.ARG ~ macrolide.exp + (1|study_id), data = test9, family = "nbinom1")
summary(num.macrolide.macrolide)
num.macrolide.macrolide <- tidy(num.macrolide.macrolide,
                               effects = "fixed",
                               component = "cond",
                               scales=NULL,
                               ran_prefix = NULL,
                               conf.int = TRUE,
                               conf.level = 0.95,
                               conf.method = "Wald")
num.macrolide.macrolide <- num.macrolide.macrolide %>%
  mutate(exposure = "macrolide", ARG_class = "Macrolide", model = "Number", effect=round(exp(estimate),2),
         lower=round(exp(conf.low), 2), upper = round(exp(conf.high), 2)) %>%
  filter(term != "(Intercept)") %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)


rpkm.macrolide.trimeth <- glmer(trimeth.abund ~ macrolide.exp + (1 |study_id), data  = test9, family=Gamma(link="log"))
summary(rpkm.macrolide.trimeth)
rpkm.macrolide.trimeth2 <- summary(rpkm.macrolide.trimeth)
rpkm.macrolide.trimeth2 <- as.data.frame(rpkm.macrolide.trimeth2$coefficients)
rpkm.cef.com <- confint(rpkm.macrolide.trimeth, method="Wald")
rpkm.cef.com <- rpkm.cef.com[-c(1, 2),]
rpkm.macrolide1 <- cbind(rpkm.macrolide.trimeth2, rpkm.cef.com)
rpkm.macrolide1 <- rpkm.macrolide1[-1,]
rpkm.macrolide.trimeth <- rpkm.macrolide1 %>% mutate(effect = round(exp(Estimate), digits =2), lower=round(exp(`2.5 %`), digits =2), upper=round(exp(`97.5 %`), digits =2)) %>%
  mutate(exposure="macrolide", ARG_class = "Diaminopyrimidine", model= "Abundance", p.value = `Pr(>|z|)`) %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)

num.macrolide.trimeth <- glmmadmb(num.trimeth.ARG ~ macrolide.exp + (1|study_id), data = test9, family = "nbinom1")
summary(num.macrolide.trimeth)
num.macrolide.trimeth <- tidy(num.macrolide.trimeth,
                             effects = "fixed",
                             component = "cond",
                             scales=NULL,
                             ran_prefix = NULL,
                             conf.int = TRUE,
                             conf.level = 0.95,
                             conf.method = "Wald")
num.macrolide.trimeth <- num.macrolide.trimeth %>%
  mutate(exposure = "macrolide", ARG_class = "Diaminopyrimidine", model = "Number", effect=round(exp(estimate),2),
         lower=round(exp(conf.low), 2), upper = round(exp(conf.high), 2)) %>%
  filter(term != "(Intercept)") %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)


rpkm.macrolide.tetra <- glmer(tetra.abund ~ macrolide.exp + (1 |study_id), data  = test9, family=Gamma(link="log"))
summary(rpkm.macrolide.tetra)
rpkm.macrolide.tetra2 <- summary(rpkm.macrolide.tetra)
rpkm.macrolide.tetra2 <- as.data.frame(rpkm.macrolide.tetra2$coefficients)
rpkm.cef.com <- confint(rpkm.macrolide.tetra, method="Wald")
rpkm.cef.com <- rpkm.cef.com[-c(1, 2),]
rpkm.macrolide1 <- cbind(rpkm.macrolide.tetra2, rpkm.cef.com)
rpkm.macrolide1 <- rpkm.macrolide1[-1,]
rpkm.macrolide.tetra <- rpkm.macrolide1 %>% mutate(effect = round(exp(Estimate), digits =2), lower=round(exp(`2.5 %`), digits =2), upper=round(exp(`97.5 %`), digits =2)) %>%
  mutate(exposure="macrolide", ARG_class = "Tetracycline", model= "Abundance", p.value = `Pr(>|z|)`) %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)

num.macrolide.tetra <- glmmadmb(num.tetra.ARG ~ macrolide.exp + (1|study_id), data = test9, family = "nbinom1")
summary(num.macrolide.tetra)
num.macrolide.tetra <- tidy(num.macrolide.tetra,
                           effects = "fixed",
                           component = "cond",
                           scales=NULL,
                           ran_prefix = NULL,
                           conf.int = TRUE,
                           conf.level = 0.95,
                           conf.method = "Wald")
num.macrolide.tetra <- num.macrolide.tetra %>%
  mutate(exposure = "macrolide", ARG_class = "Tetracycline", model = "Number", effect=round(exp(estimate),2),
         lower=round(exp(conf.low), 2), upper = round(exp(conf.high), 2)) %>%
  filter(term != "(Intercept)") %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)


rpkm.macrolide.flagyl <- glmer(flagyl.abund ~ macrolide.exp + (1 |study_id), data  = test9, family=Gamma(link="log"))
summary(rpkm.macrolide.flagyl)
rpkm.macrolide.flagyl2 <- summary(rpkm.macrolide.flagyl)
rpkm.macrolide.flagyl2 <- as.data.frame(rpkm.macrolide.flagyl2$coefficients)
rpkm.cef.com <- confint(rpkm.macrolide.flagyl, method="Wald")
rpkm.cef.com <- rpkm.cef.com[-c(1, 2),]
rpkm.macrolide1 <- cbind(rpkm.macrolide.flagyl2, rpkm.cef.com)
rpkm.macrolide1 <- rpkm.macrolide1[-1,]
rpkm.macrolide.flagyl <- rpkm.macrolide1 %>% mutate(effect = round(exp(Estimate), digits =2), lower=round(exp(`2.5 %`), digits =2), upper=round(exp(`97.5 %`), digits =2)) %>%
  mutate(exposure="macrolide", ARG_class = "Nitroimidazole", model= "Abundance", p.value = `Pr(>|z|)`) %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)

num.macrolide.flagyl <- glmmadmb(num.flagyl.ARG ~ macrolide.exp + (1|study_id), data = test9, family = "nbinom1")
summary(num.macrolide.flagyl)
num.macrolide.flagyl <- tidy(num.macrolide.flagyl,
                            effects = "fixed",
                            component = "cond",
                            scales=NULL,
                            ran_prefix = NULL,
                            conf.int = TRUE,
                            conf.level = 0.95,
                            conf.method = "Wald")
num.macrolide.flagyl <- num.macrolide.flagyl %>%
  mutate(exposure = "macrolide", ARG_class = "Nitroimidazole", model = "Number", effect=round(exp(estimate),2),
         lower=round(exp(conf.low), 2), upper = round(exp(conf.high), 2)) %>%
  filter(term != "(Intercept)") %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)

macrolide <- rbind(rpkm.macrolide.beta,
                  num.macrolide.beta,
                  rpkm.macrolide.carb,
                  num.macrolide.carb,
                  rpkm.macrolide.vanc,
                  num.macrolide.vanc,
                  rpkm.macrolide.flagyl,
                  num.macrolide.flagyl,
                  rpkm.macrolide.fluoro,
                  num.macrolide.fluoro,
                  rpkm.macrolide.trimeth,
                  num.macrolide.trimeth,
                  rpkm.macrolide.macrolide,
                  num.macrolide.macrolide,
                  rpkm.macrolide.tetra,
                  num.macrolide.tetra,
                  rpkm.macrolide.amino,
                  num.macrolide.amino,
                  rpkm.macrolide.clinda,
                  num.macrolide.clinda)
rownames(macrolide) <- NULL

macrolide <- macrolide %>%
  mutate(p.adj = p.adjust(p.value, method = "BH"))
macrolide$p.adj <- signif(macrolide$p.adj, 2)
macrolide$p.adj <- format(macrolide$p.adj, scientific = FALSE)

#other.abx exposure
rpkm.other.beta <- glmer(beta.abund ~ other.exp + (1 |study_id), data  = test9, family=Gamma(link="log"))
summary(rpkm.other.beta)
rpkm.other.beta2 <- summary(rpkm.other.beta)
rpkm.other.beta2 <- as.data.frame(rpkm.other.beta2$coefficients)
rpkm.cef.com <- confint(rpkm.other.beta, method="Wald")
rpkm.cef.com <- rpkm.cef.com[-c(1, 2),]
rpkm.other1 <- cbind(rpkm.other.beta2, rpkm.cef.com)
rpkm.other1 <- rpkm.other1[-1,]
rpkm.other.beta <- rpkm.other1 %>% mutate(effect = round(exp(Estimate), digits =2), lower=round(exp(`2.5 %`), digits =2), upper=round(exp(`97.5 %`), digits =2)) %>%
  mutate(exposure="other", ARG_class = "Beta-Lactam", model= "Abundance", p.value = `Pr(>|z|)`) %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)

num.other.betalactam <- glmmadmb(num.betalactam.ARG ~ other.exp + (1|study_id), data = test9, family = "nbinom1")
summary(num.other.betalactam)
num.other.beta <- tidy(num.other.betalactam,
                          effects = "fixed",
                          component = "cond",
                          scales=NULL,
                          ran_prefix = NULL,
                          conf.int = TRUE,
                          conf.level = 0.95,
                          conf.method = "Wald")
num.other.beta <- num.other.beta %>%
  mutate(exposure = "other", ARG_class = "Beta-Lactam", model = "Number", effect=round(exp(estimate),2),
         lower=round(exp(conf.low), 2), upper = round(exp(conf.high), 2)) %>%
  filter(term != "(Intercept)") %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)

rpkm.other.vanc <- glmer(vanc.abund ~ other.exp + (1 |study_id), data  = test9, family=Gamma(link="log"))
summary(rpkm.other.vanc)
rpkm.other.vanc2 <- summary(rpkm.other.vanc)
rpkm.other.vanc2 <- as.data.frame(rpkm.other.vanc2$coefficients)
rpkm.cef.com <- confint(rpkm.other.vanc, method="Wald")
rpkm.cef.com <- rpkm.cef.com[-c(1, 2),]
rpkm.other1 <- cbind(rpkm.other.vanc2, rpkm.cef.com)
rpkm.other1 <- rpkm.other1[-1,]
rpkm.other.vanc <- rpkm.other1 %>% mutate(effect = round(exp(Estimate), digits =2), lower=round(exp(`2.5 %`), digits =2), upper=round(exp(`97.5 %`), digits =2)) %>%
  mutate(exposure="other", ARG_class = "Glycopeptide", model= "Abundance", p.value = `Pr(>|z|)`) %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)

num.other.vanc <- glmmadmb(num.vanc.ARG ~ other.exp + (1|study_id), data = test9, family = "nbinom1")
summary(num.other.vanc)
num.other.vanc <- tidy(num.other.vanc,
                          effects = "fixed",
                          component = "cond",
                          scales=NULL,
                          ran_prefix = NULL,
                          conf.int = TRUE,
                          conf.level = 0.95,
                          conf.method = "Wald")
num.other.vanc <- num.other.vanc %>%
  mutate(exposure = "other", ARG_class = "Glycopeptide", model = "Number", effect=round(exp(estimate),2),
         lower=round(exp(conf.low), 2), upper = round(exp(conf.high), 2)) %>%
  filter(term != "(Intercept)") %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)


rpkm.other.carb <- glmer(carb.abund ~ other.exp + (1 |study_id), data  = test9, family=Gamma(link="log"))
summary(rpkm.other.carb)
rpkm.other.carb2 <- summary(rpkm.other.carb)
rpkm.other.carb2 <- as.data.frame(rpkm.other.carb2$coefficients)
rpkm.cef.com <- confint(rpkm.other.carb, method="Wald")
rpkm.cef.com <- rpkm.cef.com[-c(1, 2),]
rpkm.other1 <- cbind(rpkm.other.carb2, rpkm.cef.com)
rpkm.other1 <- rpkm.other1[-1,]
rpkm.other.carb <- rpkm.other1 %>% mutate(effect = round(exp(Estimate), digits =2), lower=round(exp(`2.5 %`), digits =2), upper=round(exp(`97.5 %`), digits =2)) %>%
  mutate(exposure="other", ARG_class = "Carbapenem", model= "Abundance", p.value = `Pr(>|z|)`) %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)

num.other.carb <- glmmadmb(num.carbapenem.ARG ~ other.exp + (1|study_id), data = test9, family = "nbinom1")
summary(num.other.carb)
num.other.carb <- tidy(num.other.carb,
                          effects = "fixed",
                          component = "cond",
                          scales=NULL,
                          ran_prefix = NULL,
                          conf.int = TRUE,
                          conf.level = 0.95,
                          conf.method = "Wald")
num.other.carb <- num.other.carb %>%
  mutate(exposure = "other", ARG_class = "Carbapenem", model = "Number", effect=round(exp(estimate),2),
         lower=round(exp(conf.low), 2), upper = round(exp(conf.high), 2)) %>%
  filter(term != "(Intercept)") %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)

rpkm.other.fluoro <- glmer(fluoro.abund ~ other.exp + (1 |study_id), data  = test9, family=Gamma(link="log"))
summary(rpkm.other.fluoro)
rpkm.other.fluoro2 <- summary(rpkm.other.fluoro)
rpkm.other.fluoro2 <- as.data.frame(rpkm.other.fluoro2$coefficients)
rpkm.cef.com <- confint(rpkm.other.fluoro, method="Wald")
rpkm.cef.com <- rpkm.cef.com[-c(1, 2),]
rpkm.other1 <- cbind(rpkm.other.fluoro2, rpkm.cef.com)
rpkm.other1 <- rpkm.other1[-1,]
rpkm.other.fluoro <- rpkm.other1 %>% mutate(effect = round(exp(Estimate), digits =2), lower=round(exp(`2.5 %`), digits =2), upper=round(exp(`97.5 %`), digits =2)) %>%
  mutate(exposure="other", ARG_class = "Fluoroquinolone", model= "Abundance", p.value = `Pr(>|z|)`) %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)

num.other.fluoro <- glmmadmb(num.quinolone.ARG ~ other.exp + (1|study_id), data = test9, family = "nbinom1")
summary(num.other.fluoro)
num.other.fluoro <- tidy(num.other.fluoro,
                            effects = "fixed",
                            component = "cond",
                            scales=NULL,
                            ran_prefix = NULL,
                            conf.int = TRUE,
                            conf.level = 0.95,
                            conf.method = "Wald")
num.other.fluoro <- num.other.fluoro %>%
  mutate(exposure = "other", ARG_class = "Fluoroquinolone", model = "Number", effect=round(exp(estimate),2),
         lower=round(exp(conf.low), 2), upper = round(exp(conf.high), 2)) %>%
  filter(term != "(Intercept)") %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)


rpkm.other.amino <- glmer(amino.abund ~ other.exp + (1 |study_id), data  = test9, family=Gamma(link="log"))
summary(rpkm.other.amino)
rpkm.other.amino2 <- summary(rpkm.other.amino)
rpkm.other.amino2 <- as.data.frame(rpkm.other.amino2$coefficients)
rpkm.cef.com <- confint(rpkm.other.amino, method="Wald")
rpkm.cef.com <- rpkm.cef.com[-c(1, 2),]
rpkm.other1 <- cbind(rpkm.other.amino2, rpkm.cef.com)
rpkm.other1 <- rpkm.other1[-1,]
rpkm.other.amino <- rpkm.other1 %>% mutate(effect = round(exp(Estimate), digits =2), lower=round(exp(`2.5 %`), digits =2), upper=round(exp(`97.5 %`), digits =2)) %>%
  mutate(exposure="other", ARG_class = "Aminoglycoside", model= "Abundance", p.value = `Pr(>|z|)`) %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)

num.other.amino <- glmmadmb(num.amino.ARG ~ other.exp + (1|study_id), data = test9, family = "nbinom1")
summary(num.other.amino)
num.other.amino <- tidy(num.other.amino,
                           effects = "fixed",
                           component = "cond",
                           scales=NULL,
                           ran_prefix = NULL,
                           conf.int = TRUE,
                           conf.level = 0.95,
                           conf.method = "Wald")
num.other.amino <- num.other.amino %>%
  mutate(exposure = "other", ARG_class = "Aminoglycoside", model = "Number", effect=round(exp(estimate),2),
         lower=round(exp(conf.low), 2), upper = round(exp(conf.high), 2)) %>%
  filter(term != "(Intercept)") %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)


rpkm.other.clinda <- glmer(clinda.abund ~ other.exp + (1 |study_id), data  = test9, family=Gamma(link="log"))
summary(rpkm.other.clinda)
rpkm.other.clinda2 <- summary(rpkm.other.clinda)
rpkm.other.clinda2 <- as.data.frame(rpkm.other.clinda2$coefficients)
rpkm.cef.com <- confint(rpkm.other.clinda, method="Wald")
rpkm.cef.com <- rpkm.cef.com[-c(1, 2),]
rpkm.other1 <- cbind(rpkm.other.clinda2, rpkm.cef.com)
rpkm.other1 <- rpkm.other1[-1,]
rpkm.other.clinda <- rpkm.other1 %>% mutate(effect = round(exp(Estimate), digits =2), lower=round(exp(`2.5 %`), digits =2), upper=round(exp(`97.5 %`), digits =2)) %>%
  mutate(exposure="other", ARG_class = "Lincosamide", model= "Abundance", p.value = `Pr(>|z|)`) %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)

num.other.clinda <- glmmadmb(num.clinda.ARG ~ other.exp + (1|study_id), data = test9, family = "nbinom1")
summary(num.other.clinda)
num.other.clinda <- tidy(num.other.clinda,
                            effects = "fixed",
                            component = "cond",
                            scales=NULL,
                            ran_prefix = NULL,
                            conf.int = TRUE,
                            conf.level = 0.95,
                            conf.method = "Wald")
num.other.clinda <- num.other.clinda %>%
  mutate(exposure = "other", ARG_class = "Lincosamide", model = "Number", effect=round(exp(estimate),2),
         lower=round(exp(conf.low), 2), upper = round(exp(conf.high), 2)) %>%
  filter(term != "(Intercept)") %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)



rpkm.other.macrolide <- glmer(macrolide.abund ~ other.exp + (1 |study_id), data  = test9, family=Gamma(link="log"))
summary(rpkm.other.macrolide)
rpkm.other.macrolide2 <- summary(rpkm.other.macrolide)
rpkm.other.macrolide2 <- as.data.frame(rpkm.other.macrolide2$coefficients)
rpkm.cef.com <- confint(rpkm.other.macrolide, method="Wald")
rpkm.cef.com <- rpkm.cef.com[-c(1, 2),]
rpkm.other1 <- cbind(rpkm.other.macrolide2, rpkm.cef.com)
rpkm.other1 <- rpkm.other1[-1,]
rpkm.other.macrolide <- rpkm.other1 %>% mutate(effect = round(exp(Estimate), digits =2), lower=round(exp(`2.5 %`), digits =2), upper=round(exp(`97.5 %`), digits =2)) %>%
  mutate(exposure="other", ARG_class = "Macrolide", model= "Abundance", p.value = `Pr(>|z|)`) %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)

num.other.macrolide <- glmmadmb(num.macrolide.ARG ~ other.exp + (1|study_id), data = test9, family = "nbinom1")
summary(num.other.macrolide)
num.other.macrolide <- tidy(num.other.macrolide,
                               effects = "fixed",
                               component = "cond",
                               scales=NULL,
                               ran_prefix = NULL,
                               conf.int = TRUE,
                               conf.level = 0.95,
                               conf.method = "Wald")
num.other.macrolide <- num.other.macrolide %>%
  mutate(exposure = "other", ARG_class = "Macrolide", model = "Number", effect=round(exp(estimate),2),
         lower=round(exp(conf.low), 2), upper = round(exp(conf.high), 2)) %>%
  filter(term != "(Intercept)") %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)


rpkm.other.trimeth <- glmer(trimeth.abund ~ other.exp + (1 |study_id), data  = test9, family=Gamma(link="log"))
summary(rpkm.other.trimeth)
rpkm.other.trimeth2 <- summary(rpkm.other.trimeth)
rpkm.other.trimeth2 <- as.data.frame(rpkm.other.trimeth2$coefficients)
rpkm.cef.com <- confint(rpkm.other.trimeth, method="Wald")
rpkm.cef.com <- rpkm.cef.com[-c(1, 2),]
rpkm.other1 <- cbind(rpkm.other.trimeth2, rpkm.cef.com)
rpkm.other1 <- rpkm.other1[-1,]
rpkm.other.trimeth <- rpkm.other1 %>% mutate(effect = round(exp(Estimate), digits =2), lower=round(exp(`2.5 %`), digits =2), upper=round(exp(`97.5 %`), digits =2)) %>%
  mutate(exposure="other", ARG_class = "Diaminopyrimidine", model= "Abundance", p.value = `Pr(>|z|)`) %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)

num.other.trimeth <- glmmadmb(num.trimeth.ARG ~ other.exp + (1|study_id), data = test9, family = "nbinom2")
summary(num.other.trimeth)
num.other.trimeth <- tidy(num.other.trimeth,
                             effects = "fixed",
                             component = "cond",
                             scales=NULL,
                             ran_prefix = NULL,
                             conf.int = TRUE,
                             conf.level = 0.95,
                             conf.method = "Wald")
num.other.trimeth <- num.other.trimeth %>%
  mutate(exposure = "other", ARG_class = "Diaminopyrimidine", model = "Number", effect=round(exp(estimate),2),
         lower=round(exp(conf.low), 2), upper = round(exp(conf.high), 2)) %>%
  filter(term != "(Intercept)") %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)


rpkm.other.tetra <- glmer(tetra.abund ~ other.exp + (1 |study_id), data  = test9, family=Gamma(link="log"))
summary(rpkm.other.tetra)
rpkm.other.tetra2 <- summary(rpkm.other.tetra)
rpkm.other.tetra2 <- as.data.frame(rpkm.other.tetra2$coefficients)
rpkm.cef.com <- confint(rpkm.other.tetra, method="Wald")
rpkm.cef.com <- rpkm.cef.com[-c(1, 2),]
rpkm.other1 <- cbind(rpkm.other.tetra2, rpkm.cef.com)
rpkm.other1 <- rpkm.other1[-1,]
rpkm.other.tetra <- rpkm.other1 %>% mutate(effect = round(exp(Estimate), digits =2), lower=round(exp(`2.5 %`), digits =2), upper=round(exp(`97.5 %`), digits =2)) %>%
  mutate(exposure="other", ARG_class = "Tetracycline", model= "Abundance", p.value = `Pr(>|z|)`) %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)

num.other.tetra <- glmmadmb(num.tetra.ARG ~ other.exp + (1|study_id), data = test9, family = "nbinom1")
summary(num.other.tetra)
num.other.tetra <- tidy(num.other.tetra,
                           effects = "fixed",
                           component = "cond",
                           scales=NULL,
                           ran_prefix = NULL,
                           conf.int = TRUE,
                           conf.level = 0.95,
                           conf.method = "Wald")
num.other.tetra <- num.other.tetra %>%
  mutate(exposure = "other", ARG_class = "Tetracycline", model = "Number", effect=round(exp(estimate),2),
         lower=round(exp(conf.low), 2), upper = round(exp(conf.high), 2)) %>%
  filter(term != "(Intercept)") %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)


rpkm.other.flagyl <- glmer(flagyl.abund ~ other.exp + (1 |study_id), data  = test9, family=Gamma(link="log"))
summary(rpkm.other.flagyl)
rpkm.other.flagyl2 <- summary(rpkm.other.flagyl)
rpkm.other.flagyl2 <- as.data.frame(rpkm.other.flagyl2$coefficients)
rpkm.cef.com <- confint(rpkm.other.flagyl, method="Wald")
rpkm.cef.com <- rpkm.cef.com[-c(1, 2),]
rpkm.other1 <- cbind(rpkm.other.flagyl2, rpkm.cef.com)
rpkm.other1 <- rpkm.other1[-1,]
rpkm.other.flagyl <- rpkm.other1 %>% mutate(effect = round(exp(Estimate), digits =2), lower=round(exp(`2.5 %`), digits =2), upper=round(exp(`97.5 %`), digits =2)) %>%
  mutate(exposure="other", ARG_class = "Nitroimidazole", model= "Abundance", p.value = `Pr(>|z|)`) %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)

num.other.flagyl <- glmmadmb(num.flagyl.ARG ~ other.exp + (1|study_id), data = test9, family = "nbinom1")
summary(num.other.flagyl)
num.other.flagyl <- tidy(num.other.flagyl,
                            effects = "fixed",
                            component = "cond",
                            scales=NULL,
                            ran_prefix = NULL,
                            conf.int = TRUE,
                            conf.level = 0.95,
                            conf.method = "Wald")
num.other.flagyl <- num.other.flagyl %>%
  mutate(exposure = "other", ARG_class = "Nitroimidazole", model = "Number", effect=round(exp(estimate),2),
         lower=round(exp(conf.low), 2), upper = round(exp(conf.high), 2)) %>%
  filter(term != "(Intercept)") %>%
  dplyr::select(model, exposure, ARG_class, effect, lower, upper, p.value)

other <- rbind(rpkm.other.beta,
                  num.other.beta,
                  rpkm.other.carb,
                  num.other.carb,
                  rpkm.other.vanc,
                  num.other.vanc,
                  rpkm.other.flagyl,
                  num.other.flagyl,
                  rpkm.other.fluoro,
                  num.other.fluoro,
                  rpkm.other.trimeth,
                  num.other.trimeth,
                  rpkm.other.macrolide,
                  num.other.macrolide,
                  rpkm.other.tetra,
                  num.other.tetra,
                  rpkm.other.amino,
                  num.other.amino,
                  rpkm.other.clinda,
                  num.other.clinda)
rownames(other) <- NULL

other <- other %>%
  mutate(p.adj = p.adjust(p.value, method = "BH"))
other$p.adj <- signif(other$p.adj, 2)
other$p.adj <- format(other$p.adj, scientific = FALSE)

#Figure 5
all.abx.arg <- rbind(cefepime, piptazo, carbapenem, vanc, flagyl, quinolone, tmp_smx, macrolide, other)

#Abundance heatmap
abund.h <- all.abx.arg %>% filter(model=="Abundance") %>%
  mutate(sig=if_else(p.adj<0.05, 8, NA)) %>%
  dplyr::select(-p.value) %>%
  dplyr::rename("IRR"="effect") %>%
  filter(exposure != "other")
  
library(scales)
median(abund.h$IRR)
summary(abund.h$IRR)
hist(abund.h$IRR)
#rescaling IRR so that colors are more different
rescale.plot.neg <- abund.h %>% filter(IRR < 1)
rescale.plot.neg$IRR[rescale.plot.neg$IRR < 0.5] <- 0.5
normalize <- function(x, na.rm = TRUE) {
  return((x- min(x)) /(max(x)-min(x)))
}
rescale.plot.neg$rescale.IRR = normalize(rescale.plot.neg$IRR)
rescale.plot.neg <- rescale.plot.neg %>% mutate(plot.IRR = rescale.IRR-1) 
summary(rescale.plot.neg$rescale.IRR)
filter(rescale.plot.neg, rescale.IRR == median(rescale.plot.neg$rescale.IRR))
rescale.plot.neg <- rescale.plot.neg %>% dplyr::select(-rescale.IRR) 
rescale.plot.pos <- abund.h %>% filter(IRR >= 1) 
rescale.plot.pos$IRR[rescale.plot.pos$IRR >2.5] <- 2.5
hist(rescale.plot.pos$IRR)
rescale.plot.pos$plot.IRR = normalize(rescale.plot.pos$IRR)
summary(rescale.plot.pos$IRR)
filter(rescale.plot.pos, plot.IRR == median(rescale.plot.pos$plot.IRR))

rescale.plot <- rbind(rescale.plot.neg, rescale.plot.pos)


heat.abund <- ggplot(rescale.plot, aes(x =ARG_class, y=exposure, fill = plot.IRR))+ 
  geom_tile(color = "black", size = 0.5) +
  geom_point(aes(x =ARG_class, y = exposure, size = sig),
             shape = 8, na.rm = TRUE, stroke = 2,
             show.legend = FALSE)+
  scale_fill_gradient2(low = "navy",  high = "orangered2", midpoint = 0,
                       breaks = c(-1, 0, 1),
                       labels = c("\u2264 0.5", "1.0", "\u2265 1.5"),
                       name = expression(beta),
                       na.value = "gray70")+
  labs(y = "Antibiotic Exposure", x = "Antibiotic Resistance Gene Class")+
  scale_x_discrete(
    limits = c("Beta-Lactam", "Carbapenem", "Glycopeptide", "Nitroimidazole", "Fluoroquinolone", "Diaminopyrimidine", "Macrolide", "Tetracycline","Aminoglycoside", "Lincosamide"),
    position = "top")+
  scale_y_discrete(
    limits = rev(c("Cefepime", "piptazo", "carbapenem", "vancomycin", "flagyl", "quinolone", "TMP-SMX", "macrolide")),
    labels = rev(c("CEF", "TZP", "CBP", "VAN", "MTZ", "FLQ", "SXT", "MAC")),
    position = "left")+
  scale_size_identity()+
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
num.h <- all.abx.arg %>% filter(model=="Number") %>%
  mutate(sig=if_else(p.adj<0.05, 8, NA)) %>%
  dplyr::select(-p.value) %>%
  dplyr::rename("IRR"="effect") %>%
  filter(exposure != "other")

library(scales)
median(num.h$IRR)
summary(num.h$IRR)
hist(num.h$IRR)

#rescaling IRR so that colors are more different
rescale.plot.neg <- num.h %>% filter(IRR < 1)
rescale.plot.neg$IRR[rescale.plot.neg$IRR < 0.5] <- 0.5
normalize <- function(x, na.rm = TRUE) {
  return((x- min(x)) /(max(x)-min(x)))
}
rescale.plot.neg$rescale.IRR = normalize(rescale.plot.neg$IRR)
rescale.plot.neg <- rescale.plot.neg %>% mutate(plot.IRR = rescale.IRR-1) 
summary(rescale.plot.neg$rescale.IRR)
filter(rescale.plot.neg, rescale.IRR == median(rescale.plot.neg$rescale.IRR))
rescale.plot.neg <- rescale.plot.neg %>% dplyr::select(-rescale.IRR) 
rescale.plot.pos <- num.h %>% filter(IRR >= 1) 
rescale.plot.pos$IRR[rescale.plot.pos$IRR >2.5] <- 2.5
hist(rescale.plot.pos$IRR)
rescale.plot.pos$plot.IRR = normalize(rescale.plot.pos$IRR)
summary(rescale.plot.pos$IRR)
filter(rescale.plot.pos, plot.IRR == median(rescale.plot.pos$plot.IRR))

rescale.plot <- rbind(rescale.plot.neg, rescale.plot.pos)


heat.num <- ggplot(rescale.plot, aes(x =ARG_class, y=exposure, fill = plot.IRR))+ 
  geom_tile(color = "black", size = 0.5) +
  geom_point(aes(x =ARG_class, y = exposure, size = sig),
             shape = 8, na.rm = TRUE, stroke = 2,
             show.legend = FALSE)+
  scale_fill_gradient2(low = "navy",  high = "orangered2", midpoint = 0,
                       breaks = c(-1, 0, 1),
                       labels = c("\u2264 0.5", "1.0", "\u2265 1.5"),
                       name = expression(beta),
                       na.value = "gray70")+
  labs(y = "Antibiotic Exposure", x = "Antibiotic Resistance Gene Class")+
  scale_x_discrete(
    limits = c("Beta-Lactam", "Carbapenem", "Glycopeptide", "Nitroimidazole", "Fluoroquinolone", "Diaminopyrimidine", "Macrolide", "Tetracycline","Aminoglycoside", "Lincosamide"),
    position = "top")+
  scale_y_discrete(
    limits = rev(c("Cefepime", "piptazo", "carbapenem", "vancomycin", "flagyl", "quinolone", "TMP-SMX", "macrolide")),
    labels = rev(c("CEF", "TZP", "CBP", "VAN", "MTZ", "FLQ", "SXT", "MAC")),
    position = "left")+
  scale_size_identity()+
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

heat <- ggarrange(heat.num , NULL, 
                  heat.abund + rremove("ylab"),  
                  labels = c("A. Number of ARGs", "",
                             "B. ARG Abundance"), 
                  vjust = 1,
                  hjust = -0.5,
                  font.label = list(size = 30),
                  nrow = 1, widths = c(1, -0.1,
                                       1.175))
heat