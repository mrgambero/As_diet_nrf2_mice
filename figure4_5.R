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
library(HH)
library(energy)
library(ecodist)
library(PResiduals)

setwd("/Users/gabri/OneDrive - University of Arizona/Arsenic Project/First_part_amp_sec/preliminary study/data/Mice")
map = readRDS("mapfile.RDS")
bac.otu = readRDS("ASV_table.RDS")
bac.otu.tax = readRDS("TAX_table.RDS")

metab_tab = read.csv("/Users/gabri/OneDrive - University of Arizona/Arsenic Project/First_part_amp_sec/preliminary study/data/metab_and_RNA/metabol_table.csv", row.names = 1)
metab_tab = metab_tab[,colnames(metab_tab) %in% rownames(map)]
metab_tab = metab_tab[, match(rownames(map),colnames(metab_tab)) ]
metab_tab_d = metab_tab[,-which(colnames(metab_tab) == "WT.C.5")]
map = map[-which(rownames(map) == "WT.C.5"),]
bac.otu = bac.otu[-which(rownames(bac.otu) == "WT.C.5"),]
metab_tab_d = log2(metab_tab_d*10000)

met.dist<-vegdist(t(metab_tab_d), method="euclidean") #maybe scale

#save_table_for_mimosa = data.frame(t(bac.otu))
#write.table( save_table_for_mimosa,"/Users/gabri/OneDrive - University of Arizona/Arsenic Project/First_part_amp_sec/preliminary study/mim_table.txt")
metab_metadata = read.csv(file = "/Users/gabri/OneDrive - University of Arizona/Arsenic Project/First_part_amp_sec/preliminary study/data/metab_and_RNA/metab_data.csv")
#metab_tab_d$KEGG = metab_metadata$KEGG[match(metab_metadata$COMP_ID,rownames(metab_tab_d))]
#metab_tab_d[metab_tab_d==""] <- NA
#metab_tab_d = metab_tab_d[!is.na(metab_tab_d$KEGG),]
#write.table( metab_tab_d,"/Users/gabri/OneDrive - University of Arizona/Arsenic Project/First_part_amp_sec/preliminary study/kegg.txt")






permanova = adonis2(met.dist ~ Genotype*Diet*Arsenic, data= map) ###anyway in the function adonis I can use the method bray too
permanova








#######Ordination#####
####PCA ordination
df_pca <- prcomp(t(metab_tab_d))
summary(df_pca)
plot(df_pca$x[,1], df_pca$x[,2])

map$Axis01<-df_pca$x[,1]
map$Axis02<-df_pca$x[,2]


mice_pca<-ggplot(map, aes(x=Axis01, y=Axis02,label = rownames(map)))+
  geom_point(aes(color=Arsenic, shape=Genotype), size=4)+
  geom_mark_hull(aes(group=Diet, label = Diet), concavity=10) +
  scale_color_social()+
  theme_classic(base_size = 15) + xlab("PC1[41.8%]") +ylab("PC2[14.4%]") + ggtitle("PCA, Liver metabolites")+
  theme(legend.position="bottom") #+ 
  #geom_text(aes(fontface=2))
mice_pca
###remove WT.C.5


#ggsave("/Users/gabri/OneDrive - University of Arizona/Arsenic Project/First_part_amp_sec/plots/figure4_new/panelA.pdf", height = 4.5, width = 5.7)

 #####First, let's do ASVs
####Select the most abundant 200 ASVs

bac.otu.rel = decostand(bac.otu, MARGIN = 1, method = "total")
rel_ab = sort(colSums(bac.otu.rel), decreasing = TRUE)[1:100] ##selecting 200 most abudnant
bac.otu.rel_200 = data.frame(t(bac.otu.rel[,colnames(bac.otu.rel) %in% names(rel_ab)]))


##Simple RDA
metab_metadata = read.csv(file = "/Users/gabri/OneDrive - University of Arizona/Arsenic Project/First_part_amp_sec/preliminary study/data/metab_and_RNA/metab_data.csv")
metaSeqObject = newMRexperiment(t(bac.otu)) #create a metagenomeseq experiment
metaSeqObject_CSS  = cumNorm( metaSeqObject , p=cumNormStatFast(metaSeqObject) ) #CSS transformation
OTU_read_count_CSS = data.frame(t(MRcounts(metaSeqObject_CSS, norm=TRUE, log=TRUE))) # retransofmr into data frame
bac.dist<-vegdist(OTU_read_count_CSS, method="bray")

###check stress for NMDS
n = 10
stress <- vector(length = n)
for (i in 1:n) {
  stress[i] <- metaMDS(bac.dist, k = i)$stress
}
names(stress) <- paste0(1:n, "Dim")
# x11(width = 10/2.54, height = 7/2.54)
par(mar = c(3.5,3.5,1,1), mgp = c(2, 0.6, 0), cex = 0.8, las = 2)
barplot(stress, ylab = "stress")

mic.nmds<-metaMDS(bac.dist, k=10)
bacteria = data.frame(mic.nmds$points)
bac.dist1<-vegdist(bacteria, method="euclidean")
###Rcor package
dcorT.test(bac.dist,bac.dist1) 
dcorT.test(met.dist,bac.dist) 

df = data.frame(bac_dist = as.vector(bac.dist), met.dist = as.vector(met.dist))
ggplot(df,aes(x=bac.dist,y=met.dist)) +
  geom_point(alpha = 0.3) + theme_classic() + 
  geom_smooth(method=lm,se=FALSE) + ylab("Metabolic profile distance") + xlab("Bacterial community distance")

#ggsave("/Users/gabri/OneDrive - University of Arizona/Arsenic Project/First_part_amp_sec/plots/figure4_new/panelB.pdf", height = 3, width = 2.5)
  #bias corrected dcor=0.7478206 p-value < 2.2e-16
dcorT.test(met.dist,bac.dist1) 



plot(met.dist~bac.dist )
plot(bac.dist~bac.dist1 )

######FOrward_selection
rda_0<-dbrda( met.dist ~ 1, bacteria)  # Model with intercept only
rda_1<-dbrda( met.dist ~ ., bacteria)  # Model with all explanatory variables
set.seed(12)
ordistep(rda_0, scope=formula(rda_1), direction="forward") #forward selection
sig_bac = bacteria[,c(1,2,7,4,6)]

final_bac = dbrda(formula = met.dist ~ MDS1 + MDS2 + MDS4 + MDS7 + MDS6, data = sig_bac)
anova(final_bac, by = "margin")

mapped = map[,c(5,6,7)]
final_map = dbrda(formula = met.dist ~ ., data = mapped)
anova(final_map, by = "margin")

vp1<-varpart(met.dist,sig_bac,mapped)
plot(vp1)
plot(final_bac)









###Correlations###
meta_tab_t = data.frame(t(metab_tab_d))
meta_tab_t <- data.frame(apply(meta_tab_t, 2, function(x) as.numeric(as.character(x))))
rel_ab = sort(colSums(bac.otu.rel), decreasing = TRUE)[1:100] ##selecting 200 most abudnant
bac.otu.rel_200 = data.frame(t(bac.otu.rel[,colnames(bac.otu.rel) %in% names(rel_ab)]))

bac.otu.rel.200.t = data.frame(t(bac.otu.rel_200))

colnames(meta_tab_t)
corr_mat = cor(bac.otu.rel.200.t,meta_tab_t, method = "spearman", use = "pairwise.complete.obs")
# significance test (fuction definition)
cor2.mtest <- function(mat1,mat2, conf.level = 0.95){
  mat1 <- as.matrix(mat1)
  mat2 <- as.matrix(mat2)
  n <-ncol(mat1)
  l <-ncol (mat2)
  p.mat  <- matrix(NA, n, l)
  
  for(i in 1:n){
    for(j in 1:l){
      tmp <- cor.test(mat1[,i], mat2[,j], conf.level = conf.level)
      p.mat[i,j]<- tmp$p.value
    }
  }
  return(list(p.mat,nrow=n,ncol=l))
}
res1 <- cor2.mtest(bac.otu.rel.200.t,meta_tab_t,0.95)  ######make significance test
pAdj <- p.adjust(c(res1[[1]]), method = "BH") #adjust p-values
resAdj <- matrix(pAdj, ncol = dim(res1[[1]])[1]) #make a matrix of adjusted p-valiues

colnames(corr_mat) = sort(unique(metab_metadata$BIOCHEMICAL))
superpath =metab_metadata$SUPER_PATHWAY[match(colnames(corr_mat), metab_metadata$BIOCHEMICAL)]
annotation <- data.frame(labels =superpath)
names <- sprintf("met[%s]",seq(1:ncol(corr_mat)))

rownames(annotation) <- names
colnames(meta_tab_t) <- names
colnames(corr_mat) = names
library(pheatmap)
dev.off()
a = pheatmap(corr_mat,
             show_rownames= TRUE,
             show_colnames= FALSE,
             #annotation = annotation,
             cutree_rows = 3,
             cutree_cols = 3,
)
a
#pdf(file = "/Users/gabri/OneDrive - University of Arizona/Arsenic Project/First_part_amp_sec/plots/figure4_new/figure4_heatmap.pdf",   # The directory you want to save the file in
#    width = 8, # The width of the plot in inches
#    height = 4) 

a
#dev.off()

row.clusters = as.hclust(a$tree_row)
col.clusters = as.hclust(a$tree_col)
plot(row.clusters)
library(factoextra)
class(row.clusters)
fviz_nbclust(corr_mat, FUN = hcut, method = "silhouette")
fviz_nbclust(t(corr_mat), FUN = hcut, method = "silhouette")
sub_grp <- cutree(row.clusters, k = 3)
sub_grp1 <- cutree(col.clusters, k = 3)
sub_grp
sub_grp1







#Bacteria
bac.otu.rel.200.t
bac.otu.rel.200.t_std = decostand(bac.otu.rel.200.t,MARGIN = 2, method = "standardize")
bac.otu.rel.200.t_std

bac.otu.rel.200.t_std_fat =  bac.otu.rel.200.t_std
bac.otu.rel.200.t_std_fat$diet = map$Diet
bac.otu.rel.200.t_std_fat = aggregate(.~diet, bac.otu.rel.200.t_std_fat, FUN = "mean")
bac.otu.rel.200.t_std_fat = data.frame(t(bac.otu.rel.200.t_std_fat))
colnames(bac.otu.rel.200.t_std_fat) = bac.otu.rel.200.t_std_fat[1,]
bac.otu.rel.200.t_std_fat =  bac.otu.rel.200.t_std_fat[-1,]
bac.otu.rel.200.t_std_fat$group = sub_grp[match(rownames(bac.otu.rel.200.t_std_fat), names(sub_grp))]
melted_fat = data.table::melt(bac.otu.rel.200.t_std_fat, id = "group")
melted_fat$value = as.numeric(melted_fat$value  )
melted_fat$group = factor(melted_fat$group)

levels(melted_fat$group)
sort(table(melted_fat$group))
melted_fat$group <- recode_factor(melted_fat$group , "1" = "B", "2"  = "A", "3" = "C")
melted_fat$group  <- factor(melted_fat$group , levels = c("A", "B", "C"))


p<-ggplot(melted_fat, aes(x=group, y=value, fill=variable)) +
  theme_classic() + ggtitle("ASVs") +ylab("Average z-score") + xlab("cluster")+ 
  scale_fill_manual(values = c("#1AFF1A", "#4B0092")) +
  geom_point(position=position_jitterdodge(),alpha=0.05, size = 1)+
  geom_boxplot( alpha=0.3, outlier.shape = NA) 
p


########Metabolites
metatab_t_std = decostand(meta_tab_t ,  MARGIN = 2, method = "standardize")
metatab_t_std_fat =  metatab_t_std
metatab_t_std_fat$diet = map$Diet
metatab_t_std_fat = aggregate(.~diet, metatab_t_std_fat, FUN = "mean")
metatab_t_std_fat = data.frame(t(metatab_t_std_fat))
colnames(metatab_t_std_fat) = metatab_t_std_fat[1,]
metatab_t_std_fat =  metatab_t_std_fat[-1,]
sub_grp1

plot(col.clusters, cex = 0.6)
rect.hclust(col.clusters, k = 3, border = 2:5)


rownames(metatab_t_std_fat)
metatab_t_std_fat$group = sub_grp1[match(rownames(metatab_t_std_fat), names(sub_grp1))]
melted_fat = data.table::melt(metatab_t_std_fat, id = "group")
melted_fat$value = as.numeric(melted_fat$value  )
levels(melted_fat$group)
melted_fat$group <- recode_factor(melted_fat$group , "1" = "C", "2"  = "A", "3" = "B")
melted_fat$group  <- factor(melted_fat$group , levels = c("A", "B", "C"))


l<-ggplot(melted_fat, aes(x=group, y=value, fill=variable)) +
  theme_classic() + ggtitle("Metabolites") +ylab("Average z-score") + xlab("Cluster")+ 
  scale_fill_manual(values = c("#1AFF1A", "#4B0092")) +
  geom_point(position=position_jitterdodge(),alpha=0.03, size = 1)+
  geom_boxplot( alpha=0.3, outlier.shape = NA) 
l


a = ggarrange(p,l, ncol = 1, common.legend = TRUE)
a


















###Work now with Indicator
###Overview here: https://aasldpubs.onlinelibrary.wiley.com/doi/pdf/10.1002/hep4.1284
####A list of metabolites possibly derived by Microbes
#1 Trimethylamine N-oxide https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7468957/
#2 Succinate https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8771265/
#3 Indoxyl sulfate https://www.nature.com/articles/s41598-021-99845-1
#4-5 Tauroursodeoxycholic acid/tauroursodeoxycholate - ursodeoxycholic acid --- https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4030606/
#6 Taurochenodeoxycholate #taurocholate
#7 Imidazole propionate #https://www.nature.com/articles/s41467-020-19589-w
#8 (4-hydroxyphenyl)lactate https://www.nature.com/articles/s41598-020-71592-9  https://aasldpubs.onlinelibrary.wiley.com/doi/full/10.1002/hep.29892
#10 4-cholesten-3-one https://www.sciencedirect.com/science/article/pii/S0140673671926341 #https://www.mdpi.com/2076-2607/9/9/1881/htm
#11 5-aminopentanoic acid https://www.nature.com/articles/s41586-021-03707-9
table_mic_dev = metab_metadata[metab_metadata$BIOCHEMICAL == "trimethylamine N-oxide" |  
                               metab_metadata$BIOCHEMICAL == "3-indoxyl sulfate" |
                               metab_metadata$BIOCHEMICAL == "deoxycholate" |
                               metab_metadata$BIOCHEMICAL == "imidazole propionate" |
                              metab_metadata$BIOCHEMICAL == "N,N,N-trimethyl-5-aminovalerate",
                               #https://pubmed.ncbi.nlm.nih.gov/32783890/
                               ]
met.dist_mic = metab_tab_d[rownames(metab_tab_d)%in% table_mic_dev$COMP_ID, ]
met.dist_mic1<-vegdist(t(met.dist_mic), method="euclidean") #maybe scale
rda_0<-dbrda( met.dist_mic1 ~ 1, bacteria)  # Model with intercept only
rda_1<-dbrda( met.dist_mic1 ~ ., bacteria)  # Model with all explanatory variables
set.seed(12)
ordistep(rda_0, scope=formula(rda_1), direction="forward") #forward selection
sig_bac = bacteria[,c(1,4,7,2,5,6)]

final_bac = dbrda(formula = met.dist ~ MDS1 +  MDS7 +  MDS5 , data = sig_bac)
anova(final_bac, by = "margin")

mapped = map[,c(5,6,7)]
final_map = dbrda(formula = met.dist ~ ., data = mapped)
anova(final_map, by = "margin")

vp1<-varpart(met.dist_mic1,sig_bac,mapped)
plot(vp1)
plot(final_bac)





####Ok, now lts see how these things vary across the experiments (insert into a loop an make a graph)
models = data.frame()
variation_parts = data.frame()
mapped1 = mapped

plot(as.numeric(met.dist_mic[1,]) ~ mapped$Diet)
met.dist_mic1 = decostand(met.dist_mic, method = "standardize", MARGIN = 1)
plot(as.numeric(met.dist_mic[1,]) ~ mapped$Genotype)

mapped1$chemical = as.numeric(met.dist_mic[1,])
model = lm(chemical ~ Diet*Genotype*Arsenic, data = mapped1)
#model$coefficients



for (i in 1:nrow(met.dist_mic))
{
  mapped1$chemical = as.numeric(met.dist_mic[i,])
  model = lm(chemical ~ Diet*Genotype*Arsenic, data = mapped1)
  a = data.frame(anova(model))
  a$coefficients = model$coefficients[c(2,3,4,5,6,7,8,1)]
  a$std_err = summary(model)$coefficients[c(2,3,4,5,6,7,8,1),2]
  a$chemical = as.character(table_mic_dev$BIOCHEMICAL[i])
  a$treatment = rownames(a)
  models = rbind(models,a)

}



models = models[!models$treatment == "Residuals",]
models$q_value = p.adjust(models$Pr..F., method = "BH")
unique(models$treatment)
models$treatment =  factor(models$treatment , levels = c("Diet", "Genotype", "Arsenic",
                                                         "Diet:Genotype","Diet:Arsenic",
                                                         "Genotype:Arsenic", "Diet:Genotype:Arsenic"))
###ggplot 
f <- ggplot(models, 
            aes(x = treatment, y = coefficients, ymin = coefficients - std_err , ymax = coefficients + std_err)) + 
  facet_wrap(~ chemical, nrow = 1) + geom_pointrange() + theme_classic() + geom_hline(yintercept = 0) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
f
#ggsave("/Users/gabri/OneDrive - University of Arizona/Arsenic Project/First_part_amp_sec/plots/figure5_cleveland.pdf", height = 3, width = 8)








######correlations with the 50 most abundant ASVs
###Correlations###
sequencesnames75 = readRDS("Seq_names.RDS")
bac.otu.rel_200 = data.frame(t(bac.otu.rel[,colnames(bac.otu.rel) %in% sequencesnames75$ASVnames]))
bac.otu.rel.200.t = data.frame(t(bac.otu.rel_200))
firm_bac_ratio  = readRDS("Firm_bac.RDS")[,8]
alpha_div = readRDS("alpha_div.RDS")[,c(10,11)]

  
bac.otu.rel.200.t$firm_bac_ratio = firm_bac_ratio[-13]
bac.otu.rel.200.t$richness = alpha_div$Richness.rar[-13]
bac.otu.rel.200.t$diversity = alpha_div$Shannon.rar[-13]

colnames(met.dist_mic)
corr_mat = cor(bac.otu.rel.200.t,t(met.dist_mic), method = "spearman", use = "pairwise.complete.obs")

# significance test (fuction definition)
resAdj = psych::corr.test(bac.otu.rel.200.t,t(met.dist_mic), method = "spearman", use = "pairwise", adjust = "fdr")
corr_mat = t(resAdj$r)
resAdj = t(resAdj$p.adj) 


cor.test(bac.otu.rel.200.t$ASV5,as.numeric(met.dist_mic[1,]), method = "spearman")

rownames(corr_mat) = metab_metadata[match(rownames(corr_mat),metab_metadata$COMP_ID),3]
colnames(corr_mat) = sequencesnames75[match(colnames(corr_mat),sequencesnames75$ASVnames),9]
colnames(corr_mat)[51] = "Firm/Bac ratio"
colnames(corr_mat)[52] = "Richness"
colnames(corr_mat)[53] = "Diversity"
which(corr_mat == max(corr_mat[2,]),arr.ind=TRUE)




rownames(resAdj) = rownames(corr_mat)
colnames(resAdj) = colnames(corr_mat)
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))

pdf(file="/Users/gabri/OneDrive - University of Arizona/Arsenic Project/First_part_amp_sec/plots/fig_5_corr_matrix.pdf")
plot2 = corrplot(corr_mat,method="ellipse", col=col(200), type = "full",# Add coefficient of correlation
                p.mat = resAdj,  
                 sig.level = 0.05, insig = "blank",
                 #addCoef.col = "black" , number.col="black",
                tl.cex	 = 0.6 , number.cex		= 0.5, 	) ####plot everything
plot2

dev.off()

correlations = data.frame() # 3-indoxyyl, use more bacs 
for(i in 1:500)
  {
      test = cor.test(bac.otu.rel[,i], as.numeric(met.dist_mic[1,]) , method = "spearman")
      cor_part = c(test$statistic,test$parameter,test$p.value, colnames(bac.otu.rel)[i])
      correlations = rbind(correlations,cor_part )
  }
colnames(correlations) = c("parameter","p.value", "name")
correlations$q_value = p.adjust(correlations$p.value, method = "fdr")



#######Look for reported correaltions####
###For 3-(4-hydroxyphenyl)lactate --- good reference https://www.frontiersin.org/articles/10.3389/fphar.2018.00653/full
Firmicutes = rownames(bac.otu.tax)[bac.otu.tax$Phylum =="Firmicutes" ]
Firmicutes = rowSums(bac.otu.rel[,colnames(bac.otu.rel)%in%Firmicutes])

pt0 = data.frame(chm = as.numeric(metab_tab_d["32197",]), 
                 ASV =  Firmicutes      ,#as.numeric(bac.otu.rel.200.t[,"ASV1"]),
                 chem = "3-(4-hydroxyphenyl)lactate", Diet = map$Diet, Genotype = map$Genotype, Arsenic = map$Arsenic) #3-(4-hydroxyphenyl)lactate

a = ggplot(pt0, aes(ASV, chm)) + geom_smooth(method = "lm",se = FALSE, formula = y ~ log(x+0.0001)) + 
  geom_point(aes(color = Arsenic, shape = Genotype, size = Diet)) + theme_classic() + ggtitle("3-(4-hydroxyphenyl)lactate) ~ Firmicutes")+
  scale_color_social() + scale_size_manual(values = c(1.5,3))
a
cor.test(Firmicutes,as.numeric(metab_tab_d["32197",]), method = "spearman")
colnames(metab_tab) == rownames(bac.otu.rel) ##seems correct




###For 3-indoxyl sulfate 
Clostridiales = rownames(bac.otu.tax)[bac.otu.tax$Order =="Clostridiales" ]
Clostridiales = rowSums(bac.otu.rel[,colnames(bac.otu.rel)%in%Clostridiales])

pt = data.frame(chm = as.numeric(metab_tab_d["27672",]), 
                ASV = Clostridiales,
                chem = "3-indoxyl sulfate", Diet = map$Diet, Genotype = map$Genotype, Arsenic = map$Arsenic) #4-cholesten-3-one


b = ggplot(pt, aes(ASV, chm)) + #geom_smooth(method = "lm",se = FALSE ) + 
  geom_point(aes(color = Arsenic, shape = Genotype, size = Diet)) + theme_classic() +
  ggtitle("3-indoxyl sulfate ~ Clostridiales")+
  scale_color_social() + scale_size_manual(values = c(1.5,3))
b #https://pubmed.ncbi.nlm.nih.gov/31815740/
cor.test(Firmicutes,as.numeric(metab_tab_d["27672",]), method = "spearman")





####Now look for correlations that 4-cholesten-3-one ~ ASV103-Lachnospiraceae https://lipidworld.biomedcentral.com/articles/10.1186/s12944-019-1103-7
Bacteroidetes = rownames(bac.otu.tax)[bac.otu.tax$Phylum =="Bacteroidetes" ]
Bacteroidetes = rowSums(bac.otu.rel[,colnames(bac.otu.rel)%in%Bacteroidetes])


pt = data.frame(chm = as.numeric(metab_tab_d["38125",]), 
                ASV = Bacteroidetes,
                chem = "4-cholesten-3-one", Diet = map$Diet, Genotype = map$Genotype, Arsenic = map$Arsenic) #4-cholesten-3-one
c = ggplot(pt, aes(ASV, chm)) + geom_smooth(method = "lm",se = FALSE ) + 
  geom_point(aes(color = Arsenic, shape = Genotype, size = Diet)) + theme_classic() +
  ggtitle("4-cholesten-3-one ~ Bacteroidetes")+
  scale_color_social() + scale_size_manual(values = c(1.5,3))
c #https://pubmed.ncbi.nlm.nih.gov/31815740/
cor.test(Bacteroidetes,as.numeric(metab_tab_d["38125",]), method = "spearman")







#######Imidazole propinate
Clostridium = rownames(bac.otu.tax)[bac.otu.tax$Genus =="Clostridium_sensu_stricto_1" ]
Clostridium = rowSums(bac.otu.rel[,colnames(bac.otu.rel)%in%Clostridium])



pt2 = data.frame(chm = as.numeric(met.dist_mic["40730",]), 
                 ASV = Clostridium,
                 chem = "Imidazole propionate",Diet = map$Diet, Genotype = map$Genotype, Arsenic = map$Arsenic)  #succinate

d = ggplot(pt2, aes(ASV,chm )) + geom_smooth(method = "lm",se = FALSE ) + 
  geom_point(aes(color = Arsenic, shape = Genotype, size = Diet)) + theme_classic() +ggtitle("Imidazole propinate ~ Clostridium")+
  scale_color_social() + scale_size_manual(values = c(1.5,3))
d
cor.test(Clostridium,as.numeric(met.dist_mic["40730",]), method = "spearman")














###succinate
Clostridium= rownames(bac.otu.tax)[bac.otu.tax$Genus =="Clostridium_sensu_stricto_1" ]
Clostridium = rowSums(bac.otu.rel[,colnames(bac.otu.rel)%in%Clostridium])




pt3 = data.frame(chm = as.numeric(met.dist_mic["1437",]), 
                 ASV = Clostridium,
                 chem = "succinate",Diet = map$Diet, Genotype = map$Genotype, Arsenic = map$Arsenic)  #succinate
e = ggplot(pt3, aes(ASV,chm )) + geom_smooth(method = "lm",se = FALSE,formula = y ~ log(x+0.0001) ) + 
  geom_point(aes(color = Arsenic, shape = Genotype, size = Diet)) + theme_classic() +ggtitle("Succinate ~ Clostridium") +
  scale_color_social() + scale_size_manual(values = c(1.5,3))
e
cor.test(Clostridium,as.numeric(met.dist_mic["1437",]), method = "spearman")


###succinate
Prevotellaceae= rownames(bac.otu.tax)[bac.otu.tax$Family =="Prevotellaceae" ]
Prevotellaceae = rowSums(bac.otu.rel[,colnames(bac.otu.rel)%in%Prevotellaceae])




pt3 = data.frame(chm = as.numeric(met.dist_mic["1437",]), 
                 ASV = Prevotellaceae,
                 chem = "succinate",Diet = map$Diet, Genotype = map$Genotype, Arsenic = map$Arsenic)  #succinate
f = ggplot(pt3, aes(ASV,chm )) + geom_smooth(method = "lm",se = FALSE, formula = y ~ log(x+0.0001)) + 
  geom_point(aes(color = Arsenic, shape = Genotype, size = Diet)) + theme_classic() +ggtitle("Succinate ~ Prevotellaceae") +
  scale_color_social() + scale_size_manual(values = c(1.5,3))
f
cor.test(Prevotellaceae,as.numeric(met.dist_mic["1437",]), method = "spearman")













####tauroursodeoxycholate
Bacteroidaceae = rownames(bac.otu.tax)[bac.otu.tax$Family =="Bacteroidaceae" ]
Bacteroidaceae = rowSums(bac.otu.rel[,colnames(bac.otu.rel)%in%Bacteroidaceae])

pt5 = data.frame(chm = as.numeric(met.dist_mic["39378",]), 
                 ASV = Bacteroidaceae,
                 chem = "tauroursodeoxycholate",Diet = map$Diet, Genotype = map$Genotype, Arsenic = map$Arsenic) #trimethylamine N-oxide

g = ggplot(pt5, aes(ASV,chm )) + #geom_smooth(method = "lm",se = FALSE ) + 
  geom_point(aes(color = Arsenic, shape = Genotype, size = Diet)) + theme_classic() +
  ggtitle("tauroursodeoxycholate ~ Bacteroidaceae")+scale_color_social() + scale_size_manual(values = c(1.5,3)) ##add removed points
g
cor.test(Bacteroidaceae,as.numeric(met.dist_mic["39378",]), method = "spearman")






Akkermansia = rownames(bac.otu.tax)[bac.otu.tax$Genus =="Akkermansia" ]
Akkermansia = rowSums(bac.otu.rel[,colnames(bac.otu.rel)%in%Akkermansia])

pt5 = data.frame(chm = as.numeric(met.dist_mic["39378",]), 
                 ASV = Akkermansia,
                 chem = "tauroursodeoxycholate",Diet = map$Diet, Genotype = map$Genotype, Arsenic = map$Arsenic) #trimethylamine N-oxide

h = ggplot(pt5, aes(ASV,chm )) + #geom_smooth(method = "lm",se = FALSE ) + 
  geom_point(aes(color = Arsenic, shape = Genotype, size = Diet)) + theme_classic() +
  ggtitle("tauroursodeoxycholate ~ Akkermansia")+scale_color_social() + scale_size_manual(values = c(1.5,3)) ##add removed points
h
cor.test(Akkermansia,as.numeric(met.dist_mic["39378",]), method = "spearman")





Prevotellaceae = rownames(bac.otu.tax)[bac.otu.tax$Family =="Prevotellaceae" ]
Prevotellaceae = rowSums(bac.otu.rel[,colnames(bac.otu.rel)%in%Prevotellaceae])

pt5 = data.frame(chm = as.numeric(met.dist_mic["39378",]), 
                 ASV = Prevotellaceae,
                 chem = "tauroursodeoxycholate",Diet = map$Diet, Genotype = map$Genotype, Arsenic = map$Arsenic) #trimethylamine N-oxide

i = ggplot(pt5, aes(ASV,chm )) + geom_smooth(method = "lm",se = FALSE, formula = y ~ log(x+0.0001) ) + 
  geom_point(aes(color = Arsenic, shape = Genotype, size = Diet)) + theme_classic() +
  ggtitle("tauroursodeoxycholate ~ Prevotellaceae")+scale_color_social() + scale_size_manual(values = c(1.5,3)) ##add removed points
i
cor.test(Prevotellaceae,as.numeric(met.dist_mic["39378",]), method = "spearman")




##Now try with #trimethylamine N-oxide"
Proteobacteria = rownames(bac.otu.tax)[bac.otu.tax$Phylum =="Proteobacteria" ]
Proteobacteria = rowSums(bac.otu.rel[,colnames(bac.otu.rel)%in%Proteobacteria])


pt7 = data.frame(chm = as.numeric(met.dist_mic["40406",]), 
                 ASV =Proteobacteria,
                 chem = "trimethylamine N-oxide",Diet = map$Diet, Genotype = map$Genotype, Arsenic = map$Arsenic ) #trimethylamine N-oxide

l = ggplot(pt7, aes(ASV, chm))+ #+ geom_smooth(method = "lm",se = FALSE ) + 
  geom_point(aes(color = Arsenic, shape = Genotype, size = Diet)) + theme_classic() +ggtitle("trimethylamine N-oxide ~ ASV43-Clostridiales")+
  scale_color_social() + scale_size_manual(values = c(1.5,3))

l
cor.test(Proteobacteria,as.numeric(met.dist_mic["40406",]), method = "spearman")


Firmicutes = rownames(bac.otu.tax)[bac.otu.tax$Phylum =="Firmicutes" ]
Firmicutes = rowSums(bac.otu.rel[,colnames(bac.otu.rel)%in%Firmicutes])


pt7 = data.frame(chm = as.numeric(met.dist_mic["40406",]), 
                 ASV =Firmicutes,
                 chem = "trimethylamine N-oxide",Diet = map$Diet, Genotype = map$Genotype, Arsenic = map$Arsenic ) #trimethylamine N-oxide

m = ggplot(pt7, aes(ASV, chm)) + geom_smooth(method = "lm",se = FALSE ) + 
  geom_point(aes(color = Arsenic, shape = Genotype, size = Diet)) + theme_classic() +ggtitle("trimethylamine N-oxide ~ Firmicutes")+
  scale_color_social() + scale_size_manual(values = c(1.5,3))

m
cor.test(Firmicutes,as.numeric(met.dist_mic["39378",]), method = "spearman")








Bacteroidetes = rownames(bac.otu.tax)[bac.otu.tax$Phylum =="Bacteroidetes" ]
Bacteroidetes = rowSums(bac.otu.rel[,colnames(bac.otu.rel)%in%Bacteroidetes])


pt7 = data.frame(chm = as.numeric(met.dist_mic["40406",]), 
                 ASV =Bacteroidetes,
                 chem = "trimethylamine N-oxide",Diet = map$Diet, Genotype = map$Genotype, Arsenic = map$Arsenic ) #trimethylamine N-oxide

n = ggplot(pt7, aes(ASV, chm)) + geom_smooth(method = "lm",se = FALSE ) + 
  geom_point(aes(color = Arsenic, shape = Genotype, size = Diet)) + theme_classic() +ggtitle("trimethylamine N-oxide ~ Bacteroidetes")+
  scale_color_social() + scale_size_manual(values = c(1.5,3))

n
cor.test(Bacteroidetes,as.numeric(met.dist_mic["39378",]), method = "spearman")



a = ggarrange(a,c,d,e,f,i,m,n, nrow = 2, ncol = 4, common.legend = TRUE)
ggsave("/Users/gabri/OneDrive - University of Arizona/Arsenic Project/First_part_amp_sec/plots/fig5_correlations.pdf", width = 7.4, height = 5)







1605


#####Ursodeoxycholate (bad correaltions, let's look at sparman)
Bacteroidaceae = rownames(bac.otu.tax)[bac.otu.tax$Family =="Bacteroidaceae" ]
Bacteroidaceae = rowSums(bac.otu.rel[,colnames(bac.otu.rel)%in%Bacteroidaceae])

pt5 = data.frame(chm = as.numeric(met.dist_mic["1605",]), 
                 ASV = Bacteroidaceae,
                 chem = "tauroursodeoxycholate",Diet = map$Diet, Genotype = map$Genotype, Arsenic = map$Arsenic) #trimethylamine N-oxide

g = ggplot(pt5, aes(ASV,chm )) + #geom_smooth(method = "lm",se = FALSE ) + 
  geom_point(aes(color = Arsenic, shape = Genotype, size = Diet)) + theme_classic() +
  ggtitle("tauroursodeoxycholate ~ Bacteroidaceae")+scale_color_social() + scale_size_manual(values = c(1.5,3)) ##add removed points
g
cor.test(Bacteroidaceae,as.numeric(met.dist_mic["1605",]), method = "spearman")






Akkermansia = rownames(bac.otu.tax)[bac.otu.tax$Genus =="Akkermansia" ]
Akkermansia = rowSums(bac.otu.rel[,colnames(bac.otu.rel)%in%Akkermansia])

pt5 = data.frame(chm = as.numeric(met.dist_mic["1605",]), 
                 ASV = Akkermansia,
                 chem = "tauroursodeoxycholate",Diet = map$Diet, Genotype = map$Genotype, Arsenic = map$Arsenic) #trimethylamine N-oxide

h = ggplot(pt5, aes(ASV,chm )) + #geom_smooth(method = "lm",se = FALSE ) + 
  geom_point(aes(color = Arsenic, shape = Genotype, size = Diet)) + theme_classic() +
  ggtitle("tauroursodeoxycholate ~ Akkermansia")+scale_color_social() + scale_size_manual(values = c(1.5,3)) ##add removed points
h
cor.test(Akkermansia,as.numeric(met.dist_mic["1605",]), method = "spearman")





Prevotellaceae = rownames(bac.otu.tax)[bac.otu.tax$Family =="Prevotellaceae" ]
Prevotellaceae = rowSums(bac.otu.rel[,colnames(bac.otu.rel)%in%Prevotellaceae])

pt5 = data.frame(chm = as.numeric(met.dist_mic["1605",]), 
                 ASV = Prevotellaceae,
                 chem = "tauroursodeoxycholate",Diet = map$Diet, Genotype = map$Genotype, Arsenic = map$Arsenic) #trimethylamine N-oxide

i = ggplot(pt5, aes(ASV,chm )) + geom_smooth(method = "lm",se = FALSE, formula = y ~ log(x+0.0001) ) + 
  geom_point(aes(color = Arsenic, shape = Genotype, size = Diet)) + theme_classic() +
  ggtitle("tauroursodeoxycholate ~ Prevotellaceae")+scale_color_social() + scale_size_manual(values = c(1.5,3)) ##add removed points
i
cor.test(Prevotellaceae,as.numeric(met.dist_mic["1605",]), method = "spearman")












