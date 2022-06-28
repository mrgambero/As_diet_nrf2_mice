###data preparation
#Mice analysis Gabri
#load libraries
library(vegan) #for multivariate statistics
library(ggplot2) #for visualization and plotting
library(performance) #for models
library(see) #for aesthetics
library(ggforce) #for hulls
library(dplyr) #for summary
library(ggpubr) ###for interaction plots  
library(metagenomeSeq) ###for beta diversity analyses
library(corrplot)
library(concaveman)

#set working directory
setwd("/Users/gabri/OneDrive - University of Arizona/Arsenic Project/First_part_amp_sec/preliminary study/data/Mice")
##########################Data Preparation#######################################################################################
#read tables
bac.otu<-read.table('All_Seqs_otu_dada.txt', sep='\t', header=T, row.names=1)
bac.tax.table<-read.table('All_Seqs_silva.txt', sep='\t', header=T, row.names=1)
map<-read.table("mice_mapping1.txt", sep="\t", header=T, row.names=1, check.names=F)


names(map)[6:7] <- c("Diet", "Arsenic") #change names of column 6 and 7 

#remove chloroplasts, mitochondria and Eukaryota, lab contamination
bac.otu<-bac.otu[,-grep("Chloroplast", bac.tax.table$Order)]
bac.tax.table<-bac.tax.table[-grep("Chloroplast", bac.tax.table$Order),]
bac.otu<-bac.otu[,-grep("Mitochondria", bac.tax.table$Family)]
bac.tax.table<-bac.tax.table[-grep("Mitochondria", bac.tax.table$Family),]
bac.otu<-bac.otu[,-grep("Archaea", bac.tax.table$Kingdom)]
bac.tax.table<-bac.tax.table[-grep("Archaea", bac.tax.table$Kingdom),]

#Blanks appear to be ok, so do not remove
#controls<-bac.otu[grep("blank", rownames(map)),] #isolate the blanks from the OTU table (does not work, ise the two lines below)
blanks = rownames(bac.otu)[1:7]
controls<-bac.otu[blanks,]
#contamination<-apply(controls, 1, function(x){which(x>0.0001*sum(x))}) #find high abundance contaminants (higher than 0.0001%)
#bac.otu<-bac.otu[,-unique(as.vector(unlist(contamination)))]
#bac.tax.table<-bac.tax.table[-unique(as.vector(unlist(contamination))),]
#now we can remove control samples
bac.otu <- bac.otu[ ! rownames(bac.otu) %in% blanks, ]
#write.table(bac.otu.tax, file="otu_table_wTax_dada2.txt", sep="\t", quote=F, row.names=T, col.names=NA)

check number of sequences per sample
hist(rowSums(bac.otu))
sort(rowSums(bac.otu))
summary(rowSums(bac.otu))

#remove NRF.As.10 (empty sample)
bac.otu<-bac.otu[-which(rownames(bac.otu)=="NRF.As.10"),]
map<-map[rownames(bac.otu),] #just check again they are all the same

#remove empty OTU
bac.otu<-subset(bac.otu, select=colSums(bac.otu)!=0)
bac.tax.table<-bac.tax.table[colnames(bac.otu),]
#create table with OTU and taxonomic information
bac.otu.tax<-cbind(t(bac.otu), bac.tax.table)

#match map with OTU table
map<-map[rownames(bac.otu),]
map$Genotype<-factor(map$Genotype, levels=c("WT","KO"))
map$Diet<-factor(map$Diet, levels=c("Control","High_fat"))
map$Arsenic<-factor(map$Arsenic, levels=c("Control","Arsenic"))
#saveRDS(map, "mapfile.RDS")
rm(controls, contamination, blanks)

###write_otu and map 

###Substitute sequences with ASV1..ASV2 etc
ASVnames=sprintf("ASV%s",seq(1:ncol(bac.otu))) ##### create a vector with new names
sequences = colnames(bac.otu) #### extract the sequences
sequencesnames=t(rbind(ASVnames, sequences)) #### create a data.frame with all this
#write.csv(sequencesnames, "/Users/gabri/OneDrive - University of Arizona/Arsenic Project/preliminary study/data/sequences.csv") ### write a CSV file with all that inside
colnames(bac.otu)= ASVnames #### substitute the names in the ASV table
rownames(bac.otu.tax) = ASVnames #### substitute the names in the ASV table
rownames(bac.tax.table) = ASVnames #### substitute the names in the ASV table




saveRDS(map, "mapfile.RDS")
saveRDS(bac.otu, "ASV_table.RDS")
saveRDS(bac.tax.table, "TAX_table.RDS")
write.table(sequencesnames, "sequences.txt")





