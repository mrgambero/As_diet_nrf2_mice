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

#check number of sequences per sample
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






#saveRDS(map, "mapfile.RDS")
#saveRDS(bac.otu, "ASV_table.RDS")
#saveRDS(bac.tax.table, "TAX_table.RDS")




#saveRDS(map,"mapfile.rds")
#######physiological data analysis######
map1<-read.table("metadata_mice_imp.txt", sep="\t", header=T, row.names=1, check.names=F) #read the physiological data
#rownames(map1) == rownames(map) #same rownames
#calculate delta weight in percent 
map1$weight_gain =  ((map1$Weight_w20 - map1$Weight_w0) /  map1$Weight_w0 ) * 100 # calculate weight gains
##remove rows that we do not need
colnames(map1)
colnames(map1)[23] = "blood (ml)" # modfy the name (otherwise too long)
map1 = map1[,-c(1,2,3,4,5,6,7,8,9,10,11,12)] # select columns to use for further analysis
#ggsave("plots/cor_mat.pdf", plot = replayPlot(plot1))
colnames(map1)
map1 = map1[,c(1,7,8,9,10,11,12)]# select columns to use for further analysis
colnames(map1)
map1 = decostand(map1, method = "standardize",MARGIN = 2) # standardize to unit of variance
map1dist = vegdist(map1,  method = "euclidean")  # calculate euclidean distances
set.seed(12)
permanova = adonis2(map1dist ~ Genotype*Diet*Arsenic, data= map) ###anyway in the function adonis I can use the method bray too
permanova
####for boxplots of physiological data see R-script "boxplots_phys.R"


#Call:
#  adonis(formula = map1dist ~ Genotype * Diet * Arsenic, data = map) 

#Permutation: free
#Number of permutations: 999

#Terms added sequentially (first to last)

#                       Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#Genotype               1    15.385  15.385   7.948 0.05784  0.003 ** 
#Diet                   1   176.780 176.780  91.326 0.66459  0.001 ***
#Arsenic                1     4.778   4.778   2.468 0.01796  0.110    
#Genotype:Diet          1     1.271   1.271   0.656 0.00478  0.504    
#Genotype:Arsenic       1     0.787   0.787   0.407 0.00296  0.702    
#Diet:Arsenic           1     6.886   6.886   3.558 0.02589  0.036 *  
#Genotype:Diet:Arsenic  1     0.105   0.105   0.054 0.00040  0.995    
#Residuals             31    60.007   1.936         0.22559           
#Total                 38   266.000                 1.00000           
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

####PCA ordination
df_pca <- prcomp(map1dist)
plot(df_pca$x[,1], df_pca$x[,2])

map$Axis01<-df_pca$x[,1]
map$Axis02<-df_pca$x[,2]


pca.summary <- map %>%
  group_by(Diet, Genotype, Arsenic) %>%
  summarise(
    axis01.mean = mean(Axis01),
    axis02.mean = mean(Axis02),
    axis01.sd = sd (Axis01, na.rm=T),
    axis02.sd = sd (Axis02, na.rm=T)
  )
pca.summary
percentage <- round(df_pca$sdev / sum(df_pca$sdev) * 100, 2)
percentage <- paste( c("PC1","PC2"), "(", paste( as.character(percentage), "%", ")", sep="") )


mice_pca_error_p<-ggplot(pca.summary, aes(x=axis01.mean, y=axis02.mean, shape=Genotype))+
  geom_errorbar(aes(color=Arsenic, ymin = axis02.mean-axis02.sd,ymax=axis02.mean+axis02.sd)) + 
  geom_errorbarh(aes(color=Arsenic, xmin = axis01.mean-axis01.sd,xmax = axis01.mean+axis01.sd)) +
  geom_point(aes(color=Arsenic), size=4)+
  geom_mark_hull(aes(group=Diet, label = Diet), concavity=10) +
  scale_color_manual(values = c("#377EB8", "#E41A1C"))+   #doesnt work
  theme_classic(base_size = 15)+
  xlab(percentage[1]) + ylab(percentage[2]) + ggtitle("Physiology")
mice_pca_error_p



####################NMDS_ordination
mic.nmds<-metaMDS(map1dist, k=2, try=100)
mic.nmds$stress #0.08
stress= paste("stress = ", round(mic.nmds$stress, digits = 2))

map$Axis01<-mic.nmds$points[,1]
map$Axis02<-mic.nmds$points[,2]



mice_nmds<-ggplot(map, aes(x=Axis01, y=Axis02))+
  geom_point(aes(color=Arsenic, shape=Genotype), size=4)+
  geom_mark_hull(aes(group=Diet, label = Diet), concavity=10) +
  scale_color_social()+
  theme_classic(base_size = 15)+
  annotate("text",x=min(map$Axis01)+ 0.87,y=max(map$Axis02),hjust=1,label= stress,  size = 3.5)
mice_nmds


############summarize the results of the nmds%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nmds.summary <- map %>%
  group_by(Diet, Genotype, Arsenic) %>%
  summarise(
    axis01.mean = mean(Axis01),
    axis02.mean = mean(Axis02),
    axis01.sd = sd (Axis01, na.rm=T),
    axis02.sd = sd (Axis02, na.rm=T)
  )
nmds.summary



mice_nmds_error_p<-ggplot(nmds.summary, aes(x=axis01.mean, y=axis02.mean, shape=Genotype))+
  geom_errorbar(aes(color=Arsenic, ymin = axis02.mean-axis02.sd,ymax=axis02.mean+axis02.sd)) + 
  geom_errorbarh(aes(color=Arsenic, xmin = axis01.mean-axis01.sd,xmax = axis01.mean+axis01.sd)) +
  geom_point(aes(color=Arsenic), size=4)+
  geom_mark_hull(aes(group=Diet, label = Diet), concavity=10) +
  scale_color_manual(values = c("#377EB8", "#E41A1C"))+   #doesnt work
  theme_classic(base_size = 15)+
  annotate("text",x=min(map$Axis01)+ 6,y=max(map$Axis02),hjust=1,label= stress,  size = 5)+
  ylab("NMDS2") +
  xlab("NMDS1") + ggtitle("Microbiome")
mice_nmds_error_p








####beta diversity patterns####################
metaSeqObject = newMRexperiment(t(bac.otu)) #create a metagenomeseq experiment
metaSeqObject_CSS  = cumNorm( metaSeqObject , p=cumNormStatFast(metaSeqObject) ) #CSS transformation
OTU_read_count_CSS = data.frame(t(MRcounts(metaSeqObject_CSS, norm=TRUE, log=TRUE))) # retransofmr into data frame
bac.dist<-vegdist(OTU_read_count_CSS, method="bray")
set.seed(1)
permanova = adonis2(bac.dist ~ Genotype*Diet*Arsenic, data= map) ###anyway in the function adonis I can use the method bray too
permanova

#Call:
#  adonis(formula = bac.dist ~ Genotype * Diet * Arsenic, data = map) 

#Permutation: free
#Number of permutations: 999

#Terms added sequentially (first to last)

#Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#Genotype               1    0.6864 0.68640  24.472 0.14835  0.001 ***
#  Diet                   1    1.6339 1.63389  58.253 0.35314  0.001 ***
#  Arsenic                1    0.3330 0.33301  11.873 0.07197  0.001 ***
#  Genotype:Diet          1    0.3467 0.34674  12.362 0.07494  0.001 ***
#  Genotype:Arsenic       1    0.2473 0.24726   8.816 0.05344  0.001 ***
#  Diet:Arsenic           1    0.2643 0.26431   9.423 0.05713  0.001 ***
#  Genotype:Diet:Arsenic  1    0.2457 0.24569   8.760 0.05310  0.001 ***
#  Residuals             31    0.8695 0.02805         0.18792           
# Total                 38    4.6268                 1.00000    






####################NMDS_ordination
bac.nmds<-metaMDS(bac.dist, k=2, try=100)
bac.nmds$stress #0.07
stress= paste("stress = ", round(bac.nmds$stress, digits = 2))

map$Axis01<-bac.nmds$points[,1]
map$Axis02<-bac.nmds$points[,2]



mice_nmds<-ggplot(map, aes(x=Axis01, y=Axis02))+
  geom_point(aes(color=Arsenic, shape=Genotype), size=4)+
  geom_mark_hull(aes(group=Diet, label = Diet), concavity=10) +
  scale_color_social()+
  theme_classic(base_size = 15)+
  annotate("text",x=min(map$Axis01)+ 0.17,y=max(map$Axis02),hjust=1,label= stress,  size = 5)
mice_nmds


############summarize the results of the NMDS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


nmds.summary <- map %>%
  group_by(Diet, Genotype, Arsenic) %>%
  summarise(
    axis01.mean = mean(Axis01),
    axis02.mean = mean(Axis02),
    axis01.sd = sd (Axis01, na.rm=T),
    axis02.sd = sd (Axis02, na.rm=T)
  )
nmds.summary



mice_nmds_error_m<-ggplot(nmds.summary, aes(x=axis01.mean, y=axis02.mean, shape=Genotype))+
  geom_errorbar(aes(color=Arsenic, ymin = axis02.mean-axis02.sd,ymax=axis02.mean+axis02.sd)) + 
  geom_errorbarh(aes(color=Arsenic, xmin = axis01.mean-axis01.sd,xmax = axis01.mean+axis01.sd)) +
  geom_point(aes(color=Arsenic), size=4)+
  geom_mark_hull(aes(group=Diet, label = Diet), concavity=10) +
  scale_color_manual(values = c("#377EB8", "#E41A1C"))+ 
  theme_classic(base_size = 15)+
  annotate("text",x=max(map$Axis01),y=max(map$Axis02),hjust=1,label= stress,  size = 5)+
  ylab("NMDS2") +
  xlab("NMDS1") + ggtitle("Microbiome")
mice_nmds_error_m

####
ggarrange(mice_pca_error_p,mice_nmds_error_m, common.legend = TRUE, labels = "AUTO")

###clean the workspace
rm(mic.nmds, metaSeqObject_CSS, metaSeqObject, mice_nmds, nmds.summary, 
   OTU_read_count_CSS, permanova, bac.dist, stress, mice_nmds_error_p,mice_nmds_error_m, cormatrix)





######Richness and shannon################################################
summary(rowSums(bac.otu))
sort(rowSums(bac.otu))
map$Richness.rar<-specnumber(rrarefy(bac.otu, 81000))
map$Shannon.rar<-diversity(rrarefy(bac.otu, 81000), index="shannon")

mean(map$Richness.rar) #274
sd(map$Richness.rar)
mean(map$Shannon.rar) #4.38


####Phylum statistics###
bac.otu_p = data.frame(t(bac.otu))
bac.otu_p$phylum =bac.otu.tax$Phylum
bac.otu_p = aggregate(. ~ phylum ,data = bac.otu_p,  FUN = sum)
rownames(bac.otu_p)  = bac.otu_p[,1]
bac.otu_p = bac.otu_p[,-1]
bac.otu_p = rowSums(bac.otu_p)
bac.otu_p = bac.otu_p/sum(bac.otu_p)
sort(-bac.otu_p)

#####richness graphs
summary(aov(Richness.rar~Genotype*Diet*Arsenic, data=map)) #same results if we invert orders of predictors
par(mfrow=c(2,2))
plot(lm(Richness.rar~Genotype*Diet*Arsenic, data=map))####good!
summary(aov(Richness.rar~Genotype*Diet*Arsenic, data=map))

richness_diet<-ggplot(map, aes(Diet, Richness.rar, fill=Diet))+
  geom_jitter(alpha=0.8)+
  geom_boxplot(aes(fill=Diet), alpha=0.3, outlier.shape = NA)+
  theme_classic()+
  #facet_wrap(~Genotype+Diet, nrow=1)+
  theme(legend.position="none")+
  ylab("Richness") +ggtitle("Diet") + xlab(NULL) +
  #scale_color_manual(values = c("#FEFE62", "#D35FB7"))+
  scale_fill_manual(values = c("#1AFF1A", "#4B0092"))+ theme(axis.text.x = element_blank()) 
richness_diet

map$Genotype = factor(map$Genotype, levels = c("WT","KO"))
richness_Genotype<-ggplot(map, aes(Genotype, Richness.rar, fill=Genotype))+
  geom_jitter(alpha=0.8)+
  geom_boxplot(aes(fill=Genotype), alpha=0.3, outlier.shape = NA)+
  theme_classic()+
  #facet_wrap(~Genotype+Diet, nrow=1)+
  theme(legend.position="none")+
  ylab(NULL) +ggtitle("Genotype")  + xlab(NULL) +
  #scale_color_manual(values = c("#FEFE62", "#D35FB7"))+
  scale_fill_manual(values = c("#00AFBB", "#E7B800")) + theme(axis.text.x = element_blank(),axis.text.y = element_blank()) 
richness_Genotype

richness_As<-ggplot(map, aes(Arsenic, Richness.rar, fill=Arsenic))+
  geom_jitter(alpha=0.8)+
  geom_boxplot(aes(fill=Arsenic), alpha=0.3, outlier.shape = NA)+
  theme_classic()+
  #facet_wrap(~Genotype+Diet, nrow=1)+
  theme(legend.position="none")+
  ylab(NULL) +ggtitle("As intake") +xlab(NULL) + theme(axis.text.x = element_blank(),axis.text.y = element_blank()) +
  scale_color_social() + scale_fill_social()
richness_As



map$center_Ars = scale(map$Arsenic, scale = FALSE, center = TRUE)
Anova(lm(Shannon.rar~Diet*Genotype*Arsenic, data=map))#same results if we invert orders of predictors
#check_model(lm(Richness.rar~Genotype*Diet*Arsenic, data=map)) does not work
par(mfrow=c(2,2))
plot(lm(Shannon.rar~Genotype*Diet*Arsenic, data=map))####good

Shannon_diet<-ggplot(map, aes(Diet, Shannon.rar, fill=Diet))+
  geom_jitter(alpha=0.8)+
  geom_boxplot(aes(fill=Diet), alpha=0.3, outlier.shape = NA)+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  #s cale_color_manual(values = c("#FEFE62", "#D35FB7"))+
  scale_fill_manual(values = c("#1AFF1A", "#4B0092")) +
  #facet_wrap(~Genotype+Diet, nrow=1)+
  theme(legend.position="none")  +
  ylab("Shannon H'") + xlab(NULL) +
  theme()
Shannon_diet

Shannon_Genotype<-ggplot(map, aes(Genotype, Shannon.rar, fill=Genotype))+
  geom_jitter(alpha=0.8)+
  geom_boxplot(aes(fill=Genotype), alpha=0.3, outlier.shape = NA)+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  #facet_wrap(~Genotype+Diet, nrow=1)+
  theme(legend.position="none")+
  ylab(NULL)  + xlab(NULL) +
  #scale_color_manual(values = c("#FEFE62", "#D35FB7"))+
  scale_fill_manual(values = c("#00AFBB", "#E7B800")) +
  theme(axis.text.y = element_blank())
Shannon_Genotype

Shannon_As<-ggplot(map, aes(Arsenic, Shannon.rar, fill=Arsenic))+
  geom_jitter(alpha=0.8)+
  geom_boxplot(aes(fill=Arsenic), alpha=0.3, outlier.shape = NA)+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  #facet_wrap(~Genotype+Diet, nrow=1)+
  theme(legend.position="none")+ xlab(NULL) +
  ylab(NULL)  +
  theme(axis.text.y = element_blank()) +
  scale_color_social() + scale_fill_social()
Shannon_As

a = ggarrange(richness_diet, richness_Genotype, richness_As, Shannon_diet, Shannon_Genotype, Shannon_As, align = "hv")
a
#ggsave("/Users/gabri/OneDrive - University of Arizona/Arsenic Project/First_part_amp_sec/plots/fig1_A.pdf",a,device = "pdf", width = 5, height = 5)

# diet_As = ggline(map, x = "Arsenic", y = "Shannon.rar", 
#                 add = c("mean_se", "jitter"),
#                 color = "Diet", facet.by = "Genotype", palette = c("#1AFF1A", "#4B0092"), ylab = "Shannon H'", title = "Shannon H', Diet:As intake:genotype", xlab = "")
# ggsave("/Users/gabri/OneDrive - University of Arizona/Arsenic Project/First_part_amp_sec/plots/fig1_B.pdf",diet_As,device = "pdf", width = 4, height = 3.5)







#diet_As = ggline(map, x = "Arsenic", y = "Shannon.rar", 
#                 add = c("mean_se", "jitter"),
#                 color = "Diet", palette = c("#1AFF1A", "#4B0092"), ylab = "Shannon H'", title = "Shannon H', Diet:As intake", xlab = "As intake")
#diet_As

library(dplyr)
df.summary2 <- map %>%
  group_by(Arsenic, Diet) %>%
  summarise(
    sd = sd(Shannon.rar),
    Shannon.rar = mean(Shannon.rar)
  )
df.summary2


diet_As1 = ggplot(map, aes( x = Arsenic, Shannon.rar)) +
  geom_jitter(
    aes(color = Diet),
    position = position_jitter(0.2),
    size = 2
  ) + 
  geom_line(
    aes(group = Diet, color = Diet),
    data = df.summary2,
    size = 0.7
  )+
  geom_errorbar(
    aes(ymin = Shannon.rar-sd, ymax = Shannon.rar+sd, color = Diet),
    data = df.summary2, width = 0.5, position = position_dodge(0.05),
    size = 0.7
  )+
  scale_color_manual(values = c("#1AFF1A", "#4B0092")) +
  theme_classic() +
  theme(legend.position="top", legend.key.size = unit(2,"point")) +
  ggtitle("Shannon H',  Diet:As intake")  + xlab ("As intake")

diet_As1


genotype_As = ggline(map, x = "Arsenic", y = "Shannon.rar", 
                     add = c("mean_se", "jitter"),
                     color = "Genotype",shape="Genotype", 
                     palette =c("#40B0A6", "#E1BE6A"), 
                     ylab = "Shannon H'", title = "Shannon H', 
                     Genotype:As intake", xlab = "As intake")
genotype_As

library(dplyr)
df.summary2 <- map %>%
  group_by(Arsenic, Genotype) %>%
  summarise(
    sd = sd(Shannon.rar),
    Shannon.rar = mean(Shannon.rar)
  )
df.summary2


genotype_As1 = ggplot(map, aes( x = Arsenic, Shannon.rar)) +
  geom_jitter(
    aes(color = Genotype, shape = Genotype),
    position = position_jitter(0.2),
    size = 2
  ) + 
  geom_line(
    aes(group = Genotype, color = Genotype),
    data = df.summary2,
    size = 0.7
  )+
  geom_errorbar(
    aes(ymin = Shannon.rar-sd, ymax = Shannon.rar+sd, color = Genotype),
    data = df.summary2, width = 0.5, position = position_dodge(0.05),
    size = 0.7
  )+
  scale_color_manual(values = c("#00AFBB", "#E7B800")) +
  theme_classic() +
  theme(legend.position="top", legend.key.size = unit(2,"point")) +
  guides(fill = guide_legend(override.aes = list(linetype = 0)),
         color = guide_legend(override.aes = list(linetype = 0))) +
  ggtitle("Shannon H',  Genotype:As intake") + xlab ("As intake")

genotype_As1



ggarrange(mice_richness, mice_shannon, diet_As1, genotype_As1, labels = "AUTO")

#clean the workspace
rm(mice_richness, mice_shannon, diet_As, genotype_As, bac.otu_p)



####Differentially abundant ASVs, DEseq2####

#Libraries we need
library(DESeq2)#The package we use
library(pheatmap)#Heat maps
library(reshape2)#Melt and data organization
library(see)#Plots
library(viridis)#Plots

###Substitute sequences with ASV1..ASV2 etc
ASVnames=sprintf("ASV%s",seq(1:ncol(bac.otu))) ##### create a vector with new names
sequences = colnames(bac.otu) #### extract the sequences
sequencesnames=t(rbind(ASVnames, sequences)) #### create a data.frame with all this
#write.csv(sequencesnames, "/Users/gabri/OneDrive - University of Arizona/Arsenic Project/preliminary study/data/sequences.csv") ### write a CSV file with all that inside
colnames(bac.otu)= ASVnames #### substitute the names in the ASV table
rownames(bac.otu.tax) = ASVnames #### substitute the names in the ASV table
rownames(bac.tax.table) = ASVnames #### substitute the names in the ASV table
rm(ASVnames,sequences,sequencesnames) ### remove useless files (keep )

#create a deseq2 object
dds <- DESeqDataSetFromMatrix(countData = t(bac.otu),  
                              colData = map,           
                              design= ~ Genotype*Arsenic*Diet)

diagdds = DESeq(dds, test="Wald", fitType="parametric")
resultsNames(diagdds) ###this are all the conditions it tested
res = results(diagdds,  name=c("GenotypeKO.ArsenicArsenic.DietHigh_fat"), pAdjustMethod = "BH", alpha =  0.05) #builds the result table
summary(res) #100 over-abundant in high fat , 304 less abundant
#DESeq2::plotMA(res) #scatter plot of log2 fold changes
#ggplot(as(res, "data.frame"), aes(x = pvalue)) +   #plot a p-value plot
  #geom_histogram(binwidth = 0.01, fill = "Royalblue", boundary = 0) 
#resLFC <- lfcShrink(diagdds, coef="Diet_High_fat_vs_Control", type="apeglm")
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
#high fat vs normal (icreased in high fat)
sigtab_Diet = cbind(as(sigtab, "data.frame"), as(bac.tax.table[rownames(sigtab), ], "matrix"))
sum(sigtab_Diet$log2FoldChange > 0) #91
###summary graph
###order the table by abundance
sigtab_Diet = sigtab_Diet[order(-sigtab_Diet$baseMean),]



####select only the 100 most abundant ones, otherwise we cannot see anything
sigtab20 = sigtab_Diet[1:20,]
sigtab20$color = ifelse(sigtab20$log2FoldChange > 0, "Increase" , "Decrease")
##### create a column to use as name, that combines ASV and Genus + species when present, otherwise just add sp.
sigtab20$name = paste(rownames(sigtab20), sigtab20$Genus, sep=" ")
species = ifelse(is.na(sigtab20$Species) , "sp." , sigtab20$Species)
sigtab20$name = factor(paste(sigtab20$name, species, sep=" "), levels=paste(sigtab20$name, species, sep=" "))
sigtab20 = sigtab20[order(sigtab20$baseMean),]
###plot it
diet = ggplot(sigtab20, 
       aes(log2FoldChange, y = name, fill = color)) + 
  geom_col(border = NA) + theme_classic() + geom_vline(xintercept =  0, size = 0.3)+
  ylab(element_blank()) +scale_y_discrete(limits = rev) + theme(legend.position = "none") +
  ggtitle("ASV, diet") + scale_fill_manual(values = c("#1AFF1A", "#4B0092")) 
diet

#ggsave( "plots/DEseq2_diet.pdf" , device ="pdf" , height = 25, width = 12)


res = results(diagdds,  name=c("Arsenic_Arsenic_vs_Control"), pAdjustMethod = "BH", alpha =  0.05) #builds the result table
summary(res) #67 over-abundant in AS , 74 less abundant
resLFC <- lfcShrink(diagdds, coef="Arsenic_Arsenic_vs_Control", type="apeglm")
alpha = 0.05
sigtab = resLFC[which(resLFC$padj < alpha), ]
sigtab_AS = cbind(as(sigtab, "data.frame"), as(bac.tax.table[rownames(sigtab), ], "matrix"))
sum(sigtab_AS$log2FoldChange > 0)
sigtab20 = sigtab_AS[1:20,]
sigtab20$color = ifelse(sigtab20$log2FoldChange > 0, "Increase" , "Decrease")
##### create a column to use as name, that combines ASV and Genus + species when present, otherwise just add sp.
sigtab20$name = paste(rownames(sigtab20), sigtab20$Genus, sep=" ")
species = ifelse(is.na(sigtab20$Species) , "sp." , sigtab20$Species)
sigtab20$name = factor(paste(sigtab20$name, species, sep=" "), levels=paste(sigtab20$name, species, sep=" "))
sigtab20 = sigtab20[order(sigtab20$baseMean),]
###plot it
AS = ggplot(sigtab20, 
              aes(log2FoldChange, y = name, fill = color)) + 
  geom_col(border = NA) + theme_classic() + geom_vline(xintercept =  0, size = 0.3)+
  ylab(element_blank()) +scale_y_discrete(limits = rev) + theme(legend.position = "none") +ggtitle("ASV, As intake")  +
  scale_color_social() + scale_fill_social()
AS


res = results(diagdds,  name=c("Genotype_KO_vs_WT"), pAdjustMethod = "BH", alpha =  0.05) #builds the result table
summary(res) #109 over-abundant in AS , 98 less abundant
resLFC <- lfcShrink(diagdds, coef="Genotype_KO_vs_WT", type="apeglm")
alpha = 0.05
sigtab = resLFC[which(resLFC$padj < alpha), ]
sigtab_GT = cbind(as(sigtab, "data.frame"), as(bac.tax.table[rownames(sigtab), ], "matrix"))
sum(sigtab_GT$log2FoldChange > 0)
sigtab20 = sigtab_GT[1:20,]
sigtab20$color = ifelse(sigtab20$log2FoldChange > 0, "Increasee" , "Decrease")
##### create a column to use GT name, that combines GTV and Genus + species when present, otherwise just add sp.
sigtab20$name = paste(rownames(sigtab20), sigtab20$Genus, sep=" ")
species = ifelse(is.na(sigtab20$Species) , "sp." , sigtab20$Species)
sigtab20$name = factor(paste(sigtab20$name, species, sep=" "), levels=paste(sigtab20$name, species, sep=" "))
sigtab20 = sigtab20[order(sigtab20$baseMean),]
###plot it
GT = ggplot(sigtab20, 
            aes(log2FoldChange, y = name, fill = color)) + 
  geom_col(border = NA) + theme_classic() + geom_vline(xintercept =  0, size = 0.3)+
  ylab(element_blank()) +scale_y_discrete(limits = rev) + theme(legend.position = "none")+ggtitle("ASV, genotype") +
  scale_fill_manual(values = c("#40B0A6", "#E1BE6A"))
GT

ggarrange(diet,AS,GT, nrow = 1)
#ggsave("plots/resultsDEseq2.pdf", height = 4, width = 12)

####find common ASVs in the list
significant_diet = rownames(sigtab_Diet)
significant_As = rownames(sigtab_AS)
significant_GT = rownames(sigtab_GT)
list = list(Diet = significant_diet,AS_intake = significant_As, Genotype= significant_GT)
library(ggvenn)
a = ggvenn(list)
a





#########DEseq2 genus level#########
bac.otu_g = bac.otu.tax[,-c(40:44,46)]
bac.otu_g$Genus = ifelse(is.na(bac.otu_g$Genus), "Unassigned", bac.otu_g$Genus)  
bac.otu_g = aggregate(. ~ Genus, FUN = "sum", data = bac.otu_g)
rownames(bac.otu_g) = bac.otu_g[,1]
bac.otu_g = bac.otu_g[,-1]

dds <- DESeqDataSetFromMatrix(countData = bac.otu_g,  
                              colData = map,           
                              design= ~ Genotype + Arsenic + Diet)

diagdds = DESeq(dds, test="Wald", fitType="parametric")
resultsNames(diagdds) ###this are all the conditions it tested
res = results(diagdds,  name=c("Diet_High_fat_vs_Control"), pAdjustMethod = "BH", alpha =  0.05) #builds the result table
summary(res) #92 over-abundant in high fat , 274 less abundant
#DESeq2::plotMA(res) #scatter plot of log2 fold changes
#ggplot(as(res, "data.frame"), aes(x = pvalue)) +   #plot a p-value plot
#geom_histogram(binwidth = 0.01, fill = "Royalblue", boundary = 0) 
resLFC <- lfcShrink(diagdds, coef="Diet_High_fat_vs_Control", type="apeglm")
alpha = 0.05
sigtab = resLFC[which(resLFC$padj < alpha), ]
#high fat vs normal (icreased in high fat)
sigtab_Diet = cbind(as(sigtab, "data.frame"), as(bac.tax.table[rownames(sigtab), ], "matrix"))
sum(sigtab_Diet$log2FoldChange > 0) #13
###summary graph
###order the table by abundance
sigtab_Diet = sigtab_Diet[order(-sigtab_Diet$baseMean),]
####select only the 100 most abundant ones, otherwise we cannot see anything
sigtab20 = sigtab_Diet[1:20,]
sigtab20$color = ifelse(sigtab20$log2FoldChange > 0, "Increase" , "Decrease")
##### create a column to use as name, that combines ASV and Genus + species when present, otherwise just add sp.
sigtab20$name = rownames(sigtab20)
sigtab20$name = factor(sigtab20$name, levels=sigtab20$name)
sigtab20 = sigtab20[order(sigtab20$baseMean),]
###plot it
diet_g = ggplot(sigtab20, 
              aes(log2FoldChange, y = name, fill = color)) + 
  geom_col(border = NA) + theme_classic() + geom_vline(xintercept =  0, size = 0.3)+
  ylab(element_blank()) +scale_y_discrete(limits = rev) + theme(legend.position = "none") +ggtitle("Genus, diet") +
  scale_fill_manual(values = c("#1AFF1A", "#4B0092")) 
  
diet_g

#ggsave( "plots/DEseq2_diet.pdf" , device ="pdf" , height = 25, width = 12)


res = results(diagdds,  name=c("Arsenic_Arsenic_vs_Control"), pAdjustMethod = "BH", alpha =  0.05) #builds the result table
summary(res) #5 over-abundant in AS , 13 less abundant
resLFC <- lfcShrink(diagdds, coef="Arsenic_Arsenic_vs_Control", type="apeglm")
alpha = 0.05
sigtab = resLFC[which(resLFC$padj < alpha), ]
sigtab_AS = cbind(as(sigtab, "data.frame"), as(bac.tax.table[rownames(sigtab), ], "matrix"))
sum(sigtab_AS$log2FoldChange > 0)
sigtab_AS = sigtab_AS[order(-sigtab_AS$baseMean),]
####select only the 100 most abundant ones, otherwise we cannot see anything
sigtab20 = sigtab_AS[1:17,]
sigtab20$color = ifelse(sigtab20$log2FoldChange > 0, "Increase" , "Decrease")
##### create a column to use as name, that combines ASV and Genus + species when present, otherwise just add sp.
sigtab20$name = rownames(sigtab20)
sigtab20$name = factor(sigtab20$name, levels=sigtab20$name)
sigtab20 = sigtab20[order(sigtab20$baseMean),]
###plot it
AS_g = ggplot(sigtab20, 
            aes(log2FoldChange, y = name, fill = color)) + 
  geom_col(border = NA) + theme_classic() + geom_vline(xintercept =  0, size = 0.3)+
  ylab(element_blank()) +scale_y_discrete(limits = rev) + theme(legend.position = "none") +ggtitle("Genus, As intake") +
   scale_color_social() + scale_fill_social()
AS_g


res = results(diagdds,  name=c("Genotype_KO_vs_WT"), pAdjustMethod = "BH", alpha =  0.05) #builds the result table
summary(res) #13 over-abundant in AS , 14 less abundant
resLFC <- lfcShrink(diagdds, coef="Genotype_KO_vs_WT", type="apeglm")
alpha = 0.05
sigtab = resLFC[which(resLFC$padj < alpha), ]
sigtab_GT = cbind(as(sigtab, "data.frame"), as(bac.tax.table[rownames(sigtab), ], "matrix"))
sum(sigtab_GT$log2FoldChange > 0)
sigtab_GT = sigtab_GT[order(-sigtab_GT$baseMean),]
####select only the 100 most abundant ones, otherwise we cannot see anything
sigtab20 = sigtab_GT[1:20,]
sigtab20$color = ifelse(sigtab20$log2FoldChange > 0, "Increase" , "Decrease")
##### create a column to use as name, that combines ASV and Genus + species when present, otherwise just add sp.
sigtab20$name = rownames(sigtab20)
sigtab20$name = factor(sigtab20$name, levels=sigtab20$name)
sigtab20 = sigtab20[order(sigtab20$baseMean),]
###plot it
GT_g = ggplot(sigtab20, 
            aes(log2FoldChange, y = name, fill = color)) + 
  geom_col(border = NA) + theme_classic() + geom_vline(xintercept =  0, size = 0.3)+
  ylab(element_blank()) +scale_y_discrete(limits = rev) + theme(legend.position = "none")+ggtitle("Genus, genotype") +
  scale_fill_manual(values = c("#40B0A6", "#E1BE6A"))
GT_g

ggarrange(diet, GT, AS, diet_g,GT_g, AS_g, nrow = 2, ncol = 3, labels = "AUTO")
ggsave("plots/resultsDEseq2.pdf", height = 8, width = 12)

####find common ASVs in the list
significant_diet = rownames(sigtab_Diet)
significant_As = rownames(sigtab_AS)
significant_GT = rownames(sigtab_GT)
list = list(Diet = significant_diet,AS_intake = significant_As, Genotype= significant_GT)
b = ggvenn(list)
ggarrange(a,b, labels = "AUTO")
#ggsave("plots/venn.pdf", height = 6, width = 10)














####PIcrust2##########
# I have to create a .biom file with the ASVtable and a FSA table with the sequences
# fsa file creation
#re-read the sequences
#sequences = read.csv("/Users/gabri/OneDrive - University of Arizona/Arsenic Project/preliminary study/data/sequences.csv") ### write a CSV file with all that inside
###create FSA (fasta) file
library(seqinr)
#write the fasta file
#write.fasta(sequences = as.list(sequences$sequences), names =  sequences$ASVnames, "/Users/gabri/OneDrive - University of Arizona/Arsenic Project/preliminary study/data/Mice/Picrust2/sequences.fsa", open = "w", nbchar = 60, as.string = FALSE)
#write biom file
library(biomformat)
#biom = make_biom(t(bac.otu), sample_metadata = NULL, observation_metadata = NULL, id = NULL, matrix_element_type = "int")
#write_biom(biom,"/Users/gabri/OneDrive - University of Arizona/Arsenic Project/preliminary study/data/Mice/Picrust2/otutable.biom")

######Run the script like this:
### in the terminal, after activating PICRUSt 2: picrust2_pipeline.py -s sequences.fsa -i otutable.biom -o picrust2_out_pipeline -p 1
### It will produce various files.
### I an interested first of all in all the Arsenic related genes
### See this paper:
### doi: 10.1007/s00248-020-01673-9
### load the file that I need:
Ko_pathways = read.delim("/Users/gabri/OneDrive - University of Arizona/Arsenic Project/First_part_amp_sec/preliminary study/data/Mice/Picrust2/picrust2_out_pipeline/KO_metagenome_out/pred_metagenome_unstrat_descrip.tsv", header = TRUE, row.names = 1,sep="\t")
head(Ko_pathways)
### we are interested in the followng genes ("ArsA", ArsB", "ArsC", "ArsR", "acr3", arsC1, Arsc2, "ArsH")
### they are present at the following codes: K01551, K03893, K18701, K03892, K03325, K00537, K03741, K11811)
Ko_pathways_rel=decostand(Ko_pathways[,-1], method = "total", MARGIN=2) # transform into relative abundances
codes = c("K01551", "K03893", "K18701", "K03892", "K03325", "K00537", "K03741", "K11811")
ko_as = Ko_pathways_rel[rownames(Ko_pathways_rel) %in% codes,]
ko_as_t = data.frame(t(ko_as))#transpose






####KWilcox_test
Wilcox_test_fat = data.frame(matrix(ncol = 2, nrow = ncol(ko_as_t)))
colnames(Wilcox_test_fat) = c("pvalue","chi")
rownames(Wilcox_test_fat) = colnames( ko_as_t )
#results_kruskal = apply(level3KEGGS_t, 2,  function(x){ kruskal.test(x, map$Diet ) })
for(i in 1:ncol(ko_as_t))
{
  a = wilcox.test(ko_as_t[,i] ~  map$Diet)
  Wilcox_test_fat[i,1] = a$p.value
  Wilcox_test_fat[i,2] = a$statistic
}
Wilcox_test_fat$test = "Diet"

Wilcox_test_As = data.frame(matrix(ncol = 2, nrow = ncol(ko_as_t)))
colnames(Wilcox_test_As) = c("pvalue","chi")
rownames(Wilcox_test_As) = colnames( ko_as_t )
#results_kruskal = apply(level3KEGGS_t, 2,  function(x){ kruskal.test(x, map$Diet ) })
for(i in 1:ncol(ko_as_t))
{
  a = wilcox.test(ko_as_t[,i] ~  map$Arsenic)
  Wilcox_test_As[i,1] = a$p.value
  Wilcox_test_As[i,2] = a$statistic
}
Wilcox_test_As$test = "As"




Wilcox_test_Genotype = data.frame(matrix(ncol = 2, nrow = ncol(ko_as_t)))
colnames(Wilcox_test_Genotype) = c("pvalue","chi")
rownames(Wilcox_test_Genotype) = colnames( ko_as_t )
#results_kruskal = apply(level3KEGGS_t, 2,  function(x){ kruskal.test(x, map$Diet ) })
for(i in 1:ncol(ko_as_t))
{
  a = kruskal.test(ko_as_t[,i]  ~ map$Genotype)
  Wilcox_test_Genotype[i,1] = a$p.value
  Wilcox_test_Genotype[i,2] = a$statistic
}
Wilcox_test_Genotype$test = "Genotype"


wilcox_level_3 = rbind(Wilcox_test_fat, Wilcox_test_As, Wilcox_test_Genotype)
wilcox_level_3$adj_p = p.adjust(wilcox_level_3$pvalue, method = "BH", n = length(wilcox_level_3$pvalue))
###########select significant ones R
wilcox_sig = data.frame(wilcox_level_3[which(wilcox_level_3$adj_p<0.05), ,drop = FALSE])
wilcox_sig


#               pvalue     chi        test   adj-P
#K01551  2.880056e-06 357.000000     Diet 2.304045e-05 #arsA, ASNA1, GET3; arsenite/tail-anchored protein-transporting ATPase [EC:3.6.3.16 3.6.3.-]
#K03325  2.232338e-06 345.000000     Diet 2.304045e-05 #ACR3, arsB; arsenite transporter
#K03892  2.232338e-06  35.000000     Diet 2.304045e-05 #ArsR family transcriptional regulator, arsenate/arsenite/antimonite-responsive transcriptional repressor
#K11811  5.472251e-06 352.000000     Diet 3.283351e-05 #arsH; arsenical resistance protein ArsH
#K037411 1.075087e-02 100.000000       As 3.225262e-02 #ARSC2, arsC; arsenate reductase [EC:1.20.4.1]
#K038931 1.157841e-04 324.500000       As 5.557635e-04 #arsB; arsenical pump membrane protein
#K005372 9.738554e-03   6.682105 Genotype 3.225262e-02 #ARSC1, arsC; arsenate reductase [EC:1.20.4.1]
#K038922 1.567520e-02   5.838947 Genotype 4.180054e-02 #arsR; ArsR family transcriptional regulator, arsenate/arsenite/antimonite-responsive transcriptional repressor
#K038932 1.003219e-02   6.629172 Genotype 3.225262e-02 #arsB; arsenical pump membrane protein


ko_as_t_g = ko_as_t 
ko_as_t_d = ko_as_t 
ko_as_t_a = ko_as_t 


ko_as_t_g$var = map$Genotype
ko_as_t_g$ct <- ifelse(ko_as_t_g$var == "KO", "Treatment", "Control")
ko_as_t_g$var = "Genotype"
ko_as_t_d$var = map$Diet
ko_as_t_d$ct <- ifelse(ko_as_t_d$var == "High_fat", "Treatment", "Control")
ko_as_t_d$var = "Diet"
ko_as_t_a$var = map$Arsenic
ko_as_t_a$ct <- ifelse(ko_as_t_a$var == "Arsenic", "Treatment", "Control")
ko_as_t_a$var = "Arsenic"
tot = rbind(ko_as_t_g, ko_as_t_d, ko_as_t_a)
tot_m = tidyr::pivot_longer(tot, col = c(1:8))

gene_names <- c(
  "K01551"="arsA",
  "K00537"="ARSC1",
  "K03325"="acr3",
  "K03741"="ARSC2",
  "K03893"="arsB",
  "K03892"="arsR",
  "K11811"="arsH",
  "K18701"="arsC"
)




p <- ggplot(tot_m, aes(x=var, y=value, fill = ct))  + facet_wrap(~ name, ncol = 4, scales = "free", labeller = as_labeller(gene_names)) +
  geom_boxplot(alpha = 0.6,outlier.shape=NA) +theme_classic() + theme(legend.position="bottom") + ylab("Relative Abundace") +xlab(NULL) + theme(legend.title = element_blank()) +
  geom_point(position = position_jitterdodge(),aes(fill = ct), alpha = 0.45, size = 0.8, color = "black", stroke = 0.9, shape = 21)+
  scale_color_grey() + scale_fill_grey(start = 0.2,
                                       end = 0.9,)
p


ggsave( "/Users/gabri/OneDrive - University of Arizona/Arsenic Project/preliminary study/data/Mice/plots/genes_boxplots.pdf" ,device = "pdf", width = 8, height = 4)




