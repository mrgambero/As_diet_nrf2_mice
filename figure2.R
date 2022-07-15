library(vegan)
library(ggdist)
library(DESeq2)#The package we use
library(pheatmap)#Heat maps
library(reshape2)#Melt and data organization
library(see)#Plots
library(viridis)#Plot
library(ggvenn)
library(ggpubr)
library(gplots)
library(ggtree)
library(ggplot2)
library(phyloseq)
library("doParallel")
library("foreach")
library("DECIPHER")
library("phangorn")
library(ggtree)


####ASV_level_analyses####
setwd("/Users/gabri/OneDrive - University of Arizona/Arsenic Project/First_part_amp_sec/preliminary study/data/Mice")
map = readRDS("mapfile.RDS")
bac.otu = readRDS("ASV_table.RDS")
bac.otu.tax = readRDS("TAX_table.RDS")




###Change seq to ASVs
sequencesnames = read.table("/Users/gabri/OneDrive - University of Arizona/Arsenic Project/First_part_amp_sec/preliminary study/data/Mice/sequences.txt") ### write a CSV file with all that inside

#create a deseq2 object
dds <- DESeqDataSetFromMatrix(countData = t(bac.otu),  
                              colData = map,           
                              design= ~ Genotype*Arsenic*Diet)
diagdds = DESeq(dds, test="Wald", fitType="parametric")
resultsNames(diagdds) ###this are all the conditions it tested
#[1] "Intercept"                              "Genotype_KO_vs_WT"                      "Arsenic_Arsenic_vs_Control"             "Diet_High_fat_vs_Control"              
#[5] "GenotypeKO.ArsenicArsenic"              "GenotypeKO.DietHigh_fat"                "ArsenicArsenic.DietHigh_fat"            "GenotypeKO.ArsenicArsenic.DietHigh_fat"

#Diet
res = results(diagdds,  name=c("Diet_High_fat_vs_Control"), pAdjustMethod = "BH", alpha =  0.05) #builds the result table
summary(res) #87 incresing in High_fat, 174 decreasing
res_diet = data.frame(res) 
res_diet = res_diet[which(res_diet$padj<0.05),]
write.table(res, "tab_diet.txt")



#genotype
res = results(diagdds,  name=c("Genotype_KO_vs_WT"), pAdjustMethod = "BH", alpha =  0.05) #builds the result table
summary(res) #53 incresing in KO, 82 decreasing
res_genotype = data.frame(res) 
res_genotype = res_genotype[which(res_genotype$padj<0.05),]
write.table(res, "tab_genotype.txt")

#Arsenic
res = results(diagdds,  name=c("Arsenic_Arsenic_vs_Control"), pAdjustMethod = "BH", alpha =  0.05) #builds the result table
summary(res) #33 incresing in KO, 51 decreasing
res_arsenic = data.frame(res) 
res_arsenic = res_arsenic[which(res_arsenic$padj<0.05),]
write.table(res, "tab_Arsenic.txt")


#Interaction Genotype Arsenic
res = results(diagdds,  name=c("GenotypeKO.ArsenicArsenic"), pAdjustMethod = "BH", alpha =  0.05) #builds the result table
summary(res) #72 incresing in KO, 70 decreasing
res_gen_arsenic = data.frame(res) 
res_gen_arsenic = res_gen_arsenic[which(res_gen_arsenic$padj<0.05),]
write.table(res, "tab_Gen_Arsenic.txt")


#Interaction Genotype Diet
res = results(diagdds,  name=c("GenotypeKO.DietHigh_fat"), pAdjustMethod = "BH", alpha =  0.05) #builds the result table
summary(res) #62 incresing in KO, 87 decreasing
res_gen_diet = data.frame(res) 
res_gen_diet = res_gen_diet[which(res_gen_diet$padj<0.05),]
write.table(res, "tab_Gen_diet.txt")


#Interaction As Diet
res = results(diagdds,  name=c("ArsenicArsenic.DietHigh_fat"), pAdjustMethod = "BH", alpha =  0.05) #builds the result table
summary(res) #42 incresing in KO, 59 decreasing
res_as_diet = data.frame(res) 
res_as_diet = res_as_diet[which(res_as_diet$padj<0.05),]
write.table(res, "tab_ars_diet.txt")

#Interaction All
res = results(diagdds,  name=c("GenotypeKO.ArsenicArsenic.DietHigh_fat"), pAdjustMethod = "BH", alpha =  0.05) #builds the result table
summary(res) #62 incresing in KO, 53 decreasing
res_all = data.frame(res) 
res_all = res_all[which(res_all$padj<0.05),]
write.table(res, "tab_ars_diet_gen.txt")


###ALL_increasing or decreasing
altogheter = unique(c(rownames(res_diet),rownames(res_genotype),rownames(res_arsenic),
                      rownames(res_gen_arsenic),rownames(res_gen_diet),rownames(res_as_diet),rownames(res_all)))



#######calculate AVG relative abundance
bac.otu = data.frame(t(bac.otu))
bac.otu_rel_ab = decostand(bac.otu, 2, method = "total")
rel_abundances = data.frame(rownames(bac.otu), means = rowMeans(bac.otu_rel_ab),abs = rowSums(bac.otu_rel_ab)  )
rel_abundances <- rel_abundances[order(rel_abundances$means, decreasing = TRUE),] 
rel_abundances_sel = rel_abundances[which(rownames(rel_abundances)%in%altogheter),]
top75 = rownames(rel_abundances_sel)[1:50]


####Now do a tree
sequencesnames75 = data.frame(sequencesnames[sequencesnames[,1] %in% top75,])
sequencesnames75$phylum = bac.otu.tax$Phylum[match(sequencesnames75$ASVnames,rownames(bac.otu.tax))]
sequencesnames75$class = bac.otu.tax$Class[match(sequencesnames75$ASVnames,rownames(bac.otu.tax))]
sequencesnames75$order = bac.otu.tax$Order[match(sequencesnames75$ASVnames,rownames(bac.otu.tax))]
sequencesnames75$family = bac.otu.tax$Family[match(sequencesnames75$ASVnames,rownames(bac.otu.tax))]
sequencesnames75$genus = bac.otu.tax$Genus[match(sequencesnames75$ASVnames,rownames(bac.otu.tax))]
sequencesnames75$species = bac.otu.tax$Species[match(sequencesnames75$ASVnames,rownames(bac.otu.tax))]
sequencesnames75$name = ifelse(!is.na(sequencesnames75$species),paste(sequencesnames75$genus,sequencesnames75$species),sequencesnames75$genus )
sequencesnames75$name = ifelse(is.na(sequencesnames75$name),sequencesnames75$family,sequencesnames75$name)
sequencesnames75$name = ifelse(is.na(sequencesnames75$name),sequencesnames75$order,sequencesnames75$name)
#sequencesnames75$name = sequencesnames75$phylum
unique(sequencesnames75$name)
sequencesnames75$name[which(sequencesnames75$name == "Rikenellaceae_RC9_gut_group")] = "Rikenellaceae_RC9"
sequencesnames75$name[which(sequencesnames75$name == "Lachnospiraceae_NK4A136_group")] = "Lachnospiraceae_NK4A136"
sequencesnames75$name[which(sequencesnames75$name == "Lachnospiraceae_NK4A136_group bacterium")] = "Lachnospiraceae_NK4A136"
sequencesnames75$name[which(sequencesnames75$name == "Clostridiales_vadinBB60_group")] = "Clostridiales_vadinBB60"
sequencesnames75$name = paste(sequencesnames75$ASVnames,sequencesnames75$name , sep = " ")
rownames(sequencesnames75) = sequencesnames75$ASVnames

#saveRDS(sequencesnames75, "Seq_names.RDS")

  seqs_1 = sequencesnames75$sequences
names(seqs_1) =rownames(sequencesnames75)
alignment <- AlignSeqs(DNAStringSet(seqs_1), anchor=NA,verbose=FALSE) 
phangAlign <- phyDat(as(alignment, "matrix"), type="DNA")

seqs_1
dm <- dist.ml(phangAlign) #about 5 mins 
treeNJ1 <- NJ(dm) # Note, tip order != sequence order
#fit = pml(treeNJ, data=phangAlign)
#fitGTR <- update(fit, k=4, inv=0.2)
#fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                   # rearrangement = "stochastic", control = pml.control(trace = 0))

#treeNJ = fitGTR$tree
#treeNJ$tip.label
#phangAlign[[1]]
#labels = names[match(treeNJ$tip.label,sequencesnames75$ASVnames),]
#labels$complete = paste(labels$asv,labels$names)
#table(labels$names) 
#labels$names = make.unique(labels$names, sep = '.')
#table(labels$names) 

treeNJ1$tip.label = sequencesnames75$name


p <- ggtree(treeNJ1, branch.length='none', layout='ellipse')  + 
  geom_tiplab(align= T, linetype=NA,  yscale = 0.2,
              size=2.5, offset=0) +xlim(0, 70)

p


m <- matrix(0, ncol =  7, nrow = nrow(sequencesnames75))
colnames(m) = c("Diet_High_fat_vs_Control","Genotype_KO_vs_WT","Arsenic_Arsenic_vs_Control",
                "GenotypeKO.ArsenicArsenic", "GenotypeKO.DietHigh_fat","ArsenicArsenic.DietHigh_fat",
                "GenotypeKO.ArsenicArsenic.DietHigh_fat")
rownames(m) = sequencesnames75$ASVnames
for (i in 1:nrow(sequencesnames75))
{
  m[i,1] = res_diet$log2FoldChange[match(rownames(m)[i],rownames(res_diet))]
  m[i,2] = res_genotype$log2FoldChange[match(rownames(m)[i], rownames(res_genotype))]  
  m[i,3] = res_arsenic$log2FoldChange[match(rownames(m)[i], rownames(res_arsenic))]  
  m[i,4] = res_gen_arsenic$log2FoldChange[match(rownames(m)[i], rownames(res_gen_arsenic))]
  m[i,5] = res_gen_diet$log2FoldChange[match(rownames(m)[i], rownames(res_gen_diet))]
  m[i,6] = res_as_diet$log2FoldChange[match(rownames(m)[i], rownames(res_as_diet))]
  m[i,7] = res_all$log2FoldChange[match(rownames(m)[i], rownames(res_all))]
}
rownames(m) = sequencesnames75$name
contrasts = m[,c(1,2,3)]
contrasts[contrasts > 2] = 3
contrasts[contrasts < -2] = -3

a = gheatmap(p, contrasts, offset = 20, color="white", 
         colnames_position="top", 
         colnames_angle=90, colnames_offset_y = 0, colnames_offset_x = 0,
         width = 0.2, 
         hjust=0, font.size=2.5)  + scale_fill_gradient2(limits = c(-3, 3), low = "#075AFF",
                                                         mid = "#FFFFCC",
                                                         high = "#FF0000", na.value="white")
a
#ggsave("/Users/gabri/OneDrive - University of Arizona/Arsenic Project/First_part_amp_sec/plots/fig4_contrast.pdf")

interactions = m[,c(4:7)]
interactions[interactions > 2] = 3
interactions[interactions < -2] = -3
gheatmap(p, interactions, offset = 20, color="white", 
         colnames_position="top", 
         colnames_angle=90, colnames_offset_y = 0, colnames_offset_x = 0,
         width = 0.2, 
         hjust=0, font.size=2.5) + scale_fill_gradient2(limits = c(-3, 3), low = "#075AFF",
                                                        mid = "#FFFFCC",
                                                        high = "#FF0000", na.value="white")
#ggsave("/Users/gabri/OneDrive - University of Arizona/Arsenic Project/First_part_amp_sec/plots/fig4_interactions.pdf")




###boxplots
###ASV1 Lactobacillus
###
rownames(map) == names(bac.otu_rel_ab["ASV1",])
df = data.frame(ASV1 = as.numeric(bac.otu_rel_ab["ASV1",]), Diet = map$Diet, Arsenic = map$Arsenic, Genotype = map$Genotype)
ASV1_diet<-ggplot(df, aes(Diet, ASV1, fill=Diet))+
  geom_jitter(alpha=0.8)+
  geom_boxplot(aes(fill=Diet), alpha=0.3, outlier.shape = NA)+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  #scale_color_manual(values = c("#FEFE62", "#D35FB7"))+
  scale_fill_manual(values = c("#1AFF1A", "#4B0092")) +
  #facet_wrap(~Genotype+Diet, nrow=1)+
  theme(legend.position="none")  +
  ylab("Relative abundance") + xlab(NULL) +
  theme(axis.text.x = element_blank())
ASV1_diet

ASV1_Genotype<-ggplot(df, aes(Genotype, ASV1, fill=Genotype))+
  geom_jitter(alpha=0.8)+
  geom_boxplot(aes(fill=Genotype), alpha=0.3, outlier.shape = NA)+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  #facet_wrap(~Genotype+Diet, nrow=1)+
  theme(legend.position="none")+
  ylab(NULL)  + xlab(NULL) +
  #scale_color_manual(values = c("#FEFE62", "#D35FB7"))+
  scale_fill_manual(values = c("#00AFBB", "#E7B800")) +
  theme(axis.text.y = element_blank(), axis.text.x = element_blank())
ASV1_Genotype

ASV1_As<-ggplot(df, aes(Arsenic, ASV1, fill=Arsenic))+
  geom_jitter(alpha=0.8)+
  geom_boxplot(aes(fill=Arsenic), alpha=0.3, outlier.shape = NA)+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  #facet_wrap(~Genotype+Diet, nrow=1)+
  theme(legend.position="none")+ xlab(NULL) +
  ylab(NULL)  +
  theme(axis.text.y = element_blank(),axis.text.x = element_blank()) +
  scale_color_social() + scale_fill_social() 

a = ggarrange(ASV1_diet, ASV1_Genotype, ASV1_As, align = "hv", ncol = 3)
a
#ggsave("/Users/gabri/OneDrive - University of Arizona/Arsenic Project/First_part_amp_sec/plots/fig4_ASV1.pdf", 
 #      width = 5, height = 2)




##ASV7
df = data.frame(ASV7 = as.numeric(bac.otu_rel_ab["ASV7",]), Diet = map$Diet, Arsenic = map$Arsenic, Genotype = map$Genotype)
ASV7_diet<-ggplot(df, aes(Diet, ASV7, fill=Diet))+
  geom_jitter(alpha=0.8)+
  geom_boxplot(aes(fill=Diet), alpha=0.3, outlier.shape = NA)+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  #scale_color_manual(values = c("#FEFE62", "#D35FB7"))+
  scale_fill_manual(values = c("#1AFF1A", "#4B0092")) +
  #facet_wrap(~Genotype+Diet, nrow=1)+
  theme(legend.position="none")  +
  ylab("Relative abundance") + xlab(NULL) +
  theme(axis.text.x = element_blank()) 
ASV7_diet

ASV7_Genotype<-ggplot(df, aes(Genotype, ASV7, fill=Genotype))+
  geom_jitter(alpha=0.8)+
  geom_boxplot(aes(fill=Genotype), alpha=0.3, outlier.shape = NA)+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  #facet_wrap(~Genotype+Diet, nrow=1)+
  theme(legend.position="none")+
  ylab(NULL)  + xlab(NULL) +
  #scale_color_manual(values = c("#FEFE62", "#D35FB7"))+
  scale_fill_manual(values = c("#00AFBB", "#E7B800")) +
  theme(axis.text.y = element_blank(),axis.text.x = element_blank())
ASV7_Genotype

ASV7_As<-ggplot(df, aes(Arsenic, ASV7, fill=Arsenic))+
  geom_jitter(alpha=0.8)+
  geom_boxplot(aes(fill=Arsenic), alpha=0.3, outlier.shape = NA)+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  #facet_wrap(~Genotype+Diet, nrow=1)+
  theme(legend.position="none")+ xlab(NULL) +
  ylab(NULL)  +
  theme(axis.text.y = element_blank(),axis.text.x = element_blank()) +
  scale_color_social() + scale_fill_social()
ASV7_As

a = ggarrange(ASV7_diet, ASV7_Genotype, ASV7_As, align = "hv", ncol = 3)
a








##ASV19
df = data.frame(ASV19 = as.numeric(bac.otu_rel_ab["ASV19",]), Diet = map$Diet, Arsenic = map$Arsenic, Genotype = map$Genotype)
ASV19_diet<-ggplot(df, aes(Diet, ASV19, fill=Diet))+
  geom_jitter(alpha=0.8)+
  geom_boxplot(aes(fill=Diet), alpha=0.3, outlier.shape = NA)+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  #scale_color_manual(values = c("#FEFE62", "#D35FB7"))+
  scale_fill_manual(values = c("#1AFF1A", "#4B0092")) +
  #facet_wrap(~Genotype+Diet, nrow=1)+
  theme(legend.position="none")  +
  ylab("Relative abundance") + xlab(NULL) 
ASV19_diet

ASV19_Genotype<-ggplot(df, aes(Genotype, ASV19, fill=Genotype))+
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
ASV19_Genotype

ASV19_As<-ggplot(df, aes(Arsenic, ASV19, fill=Arsenic))+
  geom_jitter(alpha=0.8)+
  geom_boxplot(aes(fill=Arsenic), alpha=0.3, outlier.shape = NA)+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  #facet_wrap(~Genotype+Diet, nrow=1)+
  theme(legend.position="none")+ xlab(NULL) +
  ylab(NULL)  +
  theme(axis.text.y = element_blank()) +
  scale_color_social() + scale_fill_social() 
ASV19_As

a = ggarrange(ASV1_diet, ASV1_Genotype, ASV1_As,ASV7_diet, ASV7_Genotype, ASV7_As, 
              ASV19_diet, ASV19_Genotype, ASV19_As,align = "hv", ncol = 3, nrow = 3)
a






#ggsave("/Users/gabri/OneDrive - University of Arizona/Arsenic Project/First_part_amp_sec/plots/fig4_ASV7.pdf", 
#      width = 6.5, height = 5.5)




###########Interactions########### all
df = data.frame(ASV3 = as.numeric(bac.otu_rel_ab["ASV3",]), Diet = map$Diet, Arsenic = map$Arsenic, Genotype = map$Genotype)
ASV3 = ggline(df, x = "Arsenic", y = "ASV3", 
                     add = c("mean_se", "jitter"),
                     color = "Genotype",shape="Genotype", 
                     palette =c("#40B0A6", "#E1BE6A"), 
                     ylab = "Relative Abundance", title = "ASV3", xlab = "As intake", facet.by = "Diet")
ASV3

df = data.frame(ASV5 = as.numeric(bac.otu_rel_ab["ASV5",]), Diet = map$Diet, Arsenic = map$Arsenic, Genotype = map$Genotype)
ASV5 = ggline(df, x = "Arsenic", y = "ASV5", 
              add = c("mean_se", "jitter"),
              color = "Genotype",shape="Genotype", 
              palette =c("#40B0A6", "#E1BE6A"), 
              ylab = "Relative Abundance", title = "ASV5", xlab = "As intake", facet.by = "Diet")
ASV5

ggarrange(ASV3, ASV5)
ggsave("/Users/gabri/OneDrive - University of Arizona/Arsenic Project/First_part_amp_sec/plots/fig4_inter2.pdf", 
       width = 7, height = 4.0)
 df = data.frame(ASV16 = as.numeric(bac.otu_rel_ab["ASV16",]), Diet = map$Diet, Arsenic = map$Arsenic, Genotype = map$Genotype)
ASV16 = ggline(df, x = "Diet", y = "ASV16", 
              add = c("mean_se", "jitter"),
              color = "Arsenic",shape="Arsenic", 
              palette =c("#0077b4", "#ce201f"), 
              ylab = "Relative Abundance", title = "ASV16", xlab = "Diet")
ASV16

df = data.frame(ASV10 = as.numeric(bac.otu_rel_ab["ASV10",]), Diet = map$Diet, Arsenic = map$Arsenic, Genotype = map$Genotype)
ASV10 = ggline(df, x = "Arsenic", y = "ASV10", 
              add = c("mean_se", "jitter"),
              color = "Genotype",shape="Genotype", 
              palette =c("#40B0A6", "#E1BE6A"), 
              ylab = "Relative Abundance", title = "ASV10", xlab = "As intake")
ASV10
df = data.frame(ASV13 = as.numeric(bac.otu_rel_ab["ASV13",]), Diet = map$Diet, Arsenic = map$Arsenic, Genotype = map$Genotype)

ASV13 = ggline(df, x = "Genotype", y = "ASV13", 
               add = c("mean_se", "jitter"),
               color = "Diet",shape="Diet", 
               palette =c("#1AFF1A", "#4B0092"), 
               ylab = "Relative Abundance", title = "ASV13", xlab = "Genotype")
ASV13

ggarrange(ASV10,ASV13,ASV16,ncol = 3)
#ggsave("/Users/gabri/OneDrive - University of Arizona/Arsenic Project/First_part_amp_sec/plots/fig4_inter1.pdf", 
#       width = 7, height = 3.0)
 