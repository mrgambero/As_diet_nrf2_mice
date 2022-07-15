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




models = data.frame()
variation_parts = data.frame()
mapped = map[,c(5,6,7)]
mapped1 = mapped

plot(as.numeric(metab_tab_d[1,]) ~ mapped$Diet)
metab_tab_d_s = metab_tab_d
#metab_tab_d_s = decostand(metab_tab_d, method = "standardize", MARGIN = 1)
plot(as.numeric(metab_tab_d_s[1,]) ~ mapped$Diet)

mapped1$chemical = as.numeric(metab_tab_d_s[1,])
model = lm(log(chemical,2) ~ Diet*Genotype*Arsenic, data = mapped1)

#model$coefficients



for (i in 1:nrow(metab_tab_d_s))
{
  mapped1$chemical = as.numeric(metab_tab_d_s[i,])
  model = lm(chemical ~ Diet*Genotype*Arsenic, data = mapped1)
  a = data.frame(anova(model))
  a$coefficients = model$coefficients[c(2,3,4,5,6,7,8,1)]
  a$std_err = summary(model)$coefficients[c(2,3,4,5,6,7,8,1),2]
  a$chemical = as.character(metab_metadata$BIOCHEMICAL[i])
  a$treatment = rownames(a)
  models = rbind(models,a)
  
}



models = models[!models$treatment == "Residuals",]
models$q_value = p.adjust(models$Pr..F., method = "BH")
unique(models$treatment)
models$treatment =  factor(models$treatment , levels = c("Diet", "Genotype", "Arsenic",
                                                         "Diet:Genotype","Diet:Arsenic",
                                                         "Genotype:Arsenic", "Diet:Genotype:Arsenic"))


models = models[models$q_value < 0.05,]


###Diet
models_diet = models[models$treatment == "Diet",]
models_diet = models_diet[order(models_diet$coefficients),]
write.table(models_diet, "models_diet.txt")
models_diet = models_diet[c(1:5, (nrow(models_diet)-4):nrow(models_diet)),]
models_diet$dec_inc[1:5] = "decreasing"
models_diet$dec_inc[6:10] = "increasing"
models_diet = models_diet[,c(8,11)]
models_diet$COMP_ID = metab_metadata$COMP_ID[match(models_diet$chemical, metab_metadata$BIOCHEMICAL )]
metab_tab_diet = metab_tab_d[rownames(metab_tab_d) %in% models_diet$COMP_ID,]
metab_tab_diet["ID"] <- rownames(metab_tab_diet)
melted = data.table::melt(metab_tab_diet,id.vars	= "ID")
melted$inc_dec =  models_diet$dec_inc[match(melted$ID, models_diet$COMP_ID )]
melted$ID =  models_diet$chemical[match(melted$ID, models_diet$COMP_ID )]
melted$variable = mapped1$Diet[match(melted$variable, rownames(mapped1) )]
p1 <- ggplot(melted, aes(x=ID, y=value, fill=variable)) +
  facet_wrap(~inc_dec, scale="free") + theme_classic() +
  theme(axis.text.x = element_text(angle = 45 , hjust = 1),legend.position="top")+
  #geom_point(position=position_jitterdodge(), alpha = 0.15)+
  scale_fill_manual(values = c("#1AFF1A", "#4B0092")) + 
    geom_boxplot(outlier.shape = NA,alpha=0.3) +
  xlab(NULL) +ggtitle("Diet")
p1

####Genotype
models_geno = models[models$treatment == "Genotype",]
models_geno = models_geno[order(models_geno$coefficients),]
write.table(models_geno, "models_geno.txt")

models_geno = models_geno[c(1:5, (nrow(models_geno)-4):nrow(models_geno)),]
models_geno$dec_inc[1:5] = "decreasing"
models_geno$dec_inc[6:10] = "increasing"
models_geno = models_geno[,c(8,11)]
models_geno$COMP_ID = metab_metadata$COMP_ID[match(models_geno$chemical, metab_metadata$BIOCHEMICAL )]
metab_tab_geno = metab_tab_d[rownames(metab_tab_d) %in% models_geno$COMP_ID,]
metab_tab_geno["ID"] <- rownames(metab_tab_geno)
melted = data.table::melt(metab_tab_geno,id.vars	= "ID")
melted$inc_dec =  models_geno$dec_inc[match(melted$ID, models_geno$COMP_ID )]
melted$ID =  models_geno$chemical[match(melted$ID, models_geno$COMP_ID )]
melted$variable = mapped1$Genotype[match(melted$variable, rownames(mapped1) )]
p2 <- ggplot(melted, aes(x=ID, y=value, fill=variable)) +
  facet_wrap(~inc_dec, scale="free") + theme_classic() +
  theme(axis.text.x = element_text(angle = 45 , hjust = 1),legend.position="top")+
  #geom_point(position=position_jitterdodge(), alpha = 0.15)+
  scale_fill_manual(values = c("#00AFBB", "#E7B800")) + 
  geom_boxplot(outlier.shape = NA,alpha=0.3) +
  xlab(NULL) +ggtitle("Genotype")
p2



####Arsenic
models_ars = models[models$treatment == "Arsenic",]
models_ars = models_ars[order(models_ars$coefficients),]
write.table(models_ars, "models_ars.txt")

models_ars = models_ars[c(1:5, (nrow(models_ars)-4):nrow(models_ars)),]
models_ars$dec_inc[1:5] = "decreasing"
models_ars$dec_inc[6:10] = "increasing"
models_ars = models_ars[,c(8,11)]
models_ars$COMP_ID = metab_metadata$COMP_ID[match(models_ars$chemical, metab_metadata$BIOCHEMICAL )]
metab_tab_ars = metab_tab_d[rownames(metab_tab_d) %in% models_ars$COMP_ID,]
metab_tab_ars["ID"] <- rownames(metab_tab_ars)
melted = data.table::melt(metab_tab_ars,id.vars	= "ID")
melted$inc_dec =  models_ars$dec_inc[match(melted$ID, models_ars$COMP_ID )]
melted$ID =  models_ars$chemical[match(melted$ID, models_ars$COMP_ID )]
melted$variable = mapped1$Arsenic[match(melted$variable, rownames(mapped1) )]
p3 <- ggplot(melted, aes(x=ID, y=value, fill=variable)) +
  facet_wrap(~inc_dec, scale="free") + theme_classic() +
  theme(axis.text.x = element_text(angle = 45 , hjust = 1),legend.position="top")+
  #geom_point(position=position_jitterdodge(), alpha = 0.15)+
  geom_boxplot(outlier.shape = NA,alpha=0.3) +
  xlab(NULL) +ggtitle("Arsenic")+ scale_fill_social() 
p3










####Interactions####
models_DG = models[models$treatment == "Diet:Genotype",]
models_DG = models_DG[order(models_DG$coefficients),]
write.table(models_DG, "models_DG.txt")


models_DG = models_DG[c(1:5, (nrow(models_DG)-4):nrow(models_DG)),]
models_DG$dec_inc[1:5] = "decreasing"
models_DG$dec_inc[6:10] = "increasing"
models_DG = models_DG[,c(8,11)]
models_DG$COMP_ID = metab_metadata$COMP_ID[match(models_DG$chemical, metab_metadata$BIOCHEMICAL )]

df = data.frame(Imidazole_propionate = as.numeric(metab_tab_d[rownames(metab_tab_d) == "15716",]),
                Diet = map$Diet, Arsenic = map$Arsenic, Genotype = map$Genotype)


p5 = ggline(df, x = "Genotype", y = "Imidazole_propionate", 
               add = c("mean_se", "jitter"),
               color = "Diet",shape="Diet", 
               palette =c("#1AFF1A", "#4B0092"), 
               ylab = "Normalized abundance", title = "Imidazole_propionate", xlab = "Genotype")
p5




###Interaction 2
models_DA =models[models$treatment == "Diet:Arsenic",]
models_DA = models_DA[order(models_DA$coefficients),]
write.table(models_DA, "models_DA.txt")

models_DA = models_DA[c(1:5, (nrow(models_DA)-4):nrow(models_DA)),]
models_DA$dec_inc[1:5] = "decreasing"
models_DA$dec_inc[6:10] = "increasing"
models_DA = models_DA[,c(8,11)]
models_DA$COMP_ID = metab_metadata$COMP_ID[match(models_DA$chemical, metab_metadata$BIOCHEMICAL )]

df = data.frame(sev_ketodeoxycholate = as.numeric(metab_tab_d[rownames(metab_tab_d) == models_DA$COMP_ID[10],]),
                Diet = map$Diet, Arsenic = map$Arsenic, Genotype = map$Genotype)



p6 = ggline(df, x = "Diet", y = "sev_ketodeoxycholate", 
               add = c("mean_se", "jitter"),
               color = "Arsenic",shape="Arsenic", 
               palette =c("#0077b4", "#ce201f"), 
               ylab = "Normalized abundance", title = "7-ketodeoxycholate", xlab = "Diet")
p6



###Interaction 3
models_GA = models[models$treatment == "Genotype:Arsenic",]
models_GA = models_GA[order(models_GA$coefficients),]
write.table(models_GA, "models_GA.txt")

models_GA = models_GA[c(1:5, (nrow(models_GA)-4):nrow(models_GA)),]
models_GA$dec_inc[1:5] = "decreasing"
models_GA$dec_inc[6:10] = "increasing"
models_GA = models_GA[,c(8,11)]
models_GA$COMP_ID = metab_metadata$COMP_ID[match(models_GA$chemical, metab_metadata$BIOCHEMICAL )]

df = data.frame(maltose = as.numeric(metab_tab_d[rownames(metab_tab_d) == models_GA$COMP_ID[1],]),
                Diet = map$Diet, Arsenic = map$Arsenic, Genotype = map$Genotype)



p7 = ggline(df, x = "Arsenic", y = "maltose", 
             add = c("mean_se", "jitter"),
             color = "Genotype",shape="Genotype", 
             palette =c("#00AFBB", "#E7B800"), 
             ylab = "Normalized abundance", title = "Maltose", xlab = "Arsenic")
p7


###Interaction 4
models_DAG = models[models$treatment == "Diet:Genotype:Arsenic",]
models_DAG = models_DAG[order(models_DAG$coefficients),]
write.table(models_DAG, "models_DAG.txt")


models_DAG = models_DAG[c(1:5, (nrow(models_DAG)-4):nrow(models_DAG)),]
models_DAG$dec_inc[1:5] = "decreasing"
models_DAG$dec_inc[6:10] = "increasing"
models_DAG = models_DAG[,c(8,11)]
models_DAG$COMP_ID = metab_metadata$COMP_ID[match(models_DAG$chemical, metab_metadata$BIOCHEMICAL )]

df = data.frame(TMAVA = as.numeric(metab_tab_d[rownames(metab_tab_d) == models_DAG$COMP_ID[5],]),
                Diet = map$Diet, Arsenic = map$Arsenic, Genotype = map$Genotype)


p8 = ggline(df, x = "Arsenic", y = "TMAVA", 
              add = c("mean_se", "jitter"),
              color = "Genotype",shape="Genotype", 
              palette =c("#40B0A6", "#E1BE6A"), 
              ylab = "Relative Abundance", title = "TMAVA", xlab = "As intake", facet.by = "Diet")
p8



ggarrange(p1,p2,p3, common.legend = TRUE, nrow = 1)
#ggsave("/Users/gabri/OneDrive - University of Arizona/Arsenic Project/First_part_amp_sec/plots/figure_rel_ab_met/rel_ab_11.pdf", height = 4)
ggarrange(p5,p6,p7,p8, common.legend = TRUE, nrow = 1)
#ggsave("/Users/gabri/OneDrive - University of Arizona/Arsenic Project/First_part_amp_sec/plots/figure_rel_ab_met/rel_ab_22.pdf", height = 3, width = 8)














#############Subpathway
#############Subpathway
#############Subpathway



metab_tab_d_s = decostand(metab_tab_d, method = "standardize", MARGIN = 1)
###aggregate by other column 
metab_tab_d_s$sub_pathway = metab_metadata$SUB_PATHWAY
grouped = aggregate(.~sub_pathway, metab_tab_d_s, FUN = "mean") %>% tibble::column_to_rownames(var="sub_pathway") 


  
  
models = data.frame()
plot(as.numeric(grouped[1,]) ~ mapped$Diet)
metab_tab_d_s = grouped
mapped1$chemical = as.numeric(metab_tab_d_s[1,])

#model$coefficients



for (i in 1:nrow(metab_tab_d_s))
{
  mapped1$chemical = as.numeric(metab_tab_d_s[i,])
  model = lm(chemical ~ Diet*Genotype*Arsenic, data = mapped1)
  a = data.frame(anova(model))
  a$coefficients = model$coefficients[c(2,3,4,5,6,7,8,1)]
  a$std_err = summary(model)$coefficients[c(2,3,4,5,6,7,8,1),2]
  a$chemical = as.character(rownames(metab_tab_d_s)[i])
  a$treatment = rownames(a)
  models = rbind(models,a)
  
}



models = models[!models$treatment == "Residuals",]
models$q_value = p.adjust(models$Pr..F., method = "BH")
unique(models$treatment)
models$treatment =  factor(models$treatment , levels = c("Diet", "Genotype", "Arsenic",
                                                         "Diet:Genotype","Diet:Arsenic",
                                                         "Genotype:Arsenic", "Diet:Genotype:Arsenic"))


models = models[models$q_value < 0.05,]


###Diet
models_diet = models[models$treatment == "Diet",]
models_diet = models_diet[order(models_diet$coefficients),]
#write.table(models_diet, "models_diet.txt")
models_diet = models_diet[c(1:5, (nrow(models_diet)-4):nrow(models_diet)),]
models_diet$dec_inc[1:5] = "decreasing"
models_diet$dec_inc[6:10] = "increasing"
models_diet = models_diet[,c(8,11)]

metab_tab_diet = metab_tab_d[rownames(metab_tab_d) %in% models_diet$COMP_ID,]
metab_tab_diet["ID"] <- rownames(metab_tab_diet)


melted = data.table::melt(metab_tab_diet,id.vars	= "ID")
  melted$inc_dec =  models_diet$dec_inc[match(melted$ID, models_diet$COMP_ID )]
melted$ID =  models_diet$chemical[match(melted$ID, models_diet$COMP_ID )]
melted$variable = mapped1$Diet[match(melted$variable, rownames(mapped1) )]
p1 <- ggplot(melted, aes(x=ID, y=value, fill=variable)) +
  facet_wrap(~inc_dec, scale="free") + theme_classic() +
  theme(axis.text.x = element_text(angle = 45 , hjust = 1),legend.position="top")+
  #geom_point(position=position_jitterdodge(), alpha = 0.15)+
  scale_fill_manual(values = c("#1AFF1A", "#4B0092")) + 
  geom_boxplot(outlier.shape = NA,alpha=0.3) +
  xlab(NULL) +ggtitle("Diet")
p1

####Genotype
models_geno = models[models$treatment == "Genotype",]
models_geno = models_geno[order(models_geno$coefficients),]
write.table(models_geno, "models_geno.txt")

models_geno = models_geno[c(1:5, (nrow(models_geno)-4):nrow(models_geno)),]
models_geno$dec_inc[1:5] = "decreasing"
models_geno$dec_inc[6:10] = "increasing"
models_geno = models_geno[,c(8,11)]
models_geno$COMP_ID = metab_metadata$COMP_ID[match(models_geno$chemical, metab_metadata$BIOCHEMICAL )]
metab_tab_geno = metab_tab_d[rownames(metab_tab_d) %in% models_geno$COMP_ID,]
metab_tab_geno["ID"] <- rownames(metab_tab_geno)
melted = data.table::melt(metab_tab_geno,id.vars	= "ID")
melted$inc_dec =  models_geno$dec_inc[match(melted$ID, models_geno$COMP_ID )]
melted$ID =  models_geno$chemical[match(melted$ID, models_geno$COMP_ID )]
melted$variable = mapped1$Genotype[match(melted$variable, rownames(mapped1) )]
p2 <- ggplot(melted, aes(x=ID, y=value, fill=variable)) +
  facet_wrap(~inc_dec, scale="free") + theme_classic() +
  theme(axis.text.x = element_text(angle = 45 , hjust = 1),legend.position="top")+
  #geom_point(position=position_jitterdodge(), alpha = 0.15)+
  scale_fill_manual(values = c("#00AFBB", "#E7B800")) + 
  geom_boxplot(outlier.shape = NA,alpha=0.3) +
  xlab(NULL) +ggtitle("Genotype")
p2



####Arsenic
models_ars = models[models$treatment == "Arsenic",]
models_ars = models_ars[order(models_ars$coefficients),]
write.table(models_ars, "models_ars.txt")

models_ars = models_ars[c(1:5, (nrow(models_ars)-4):nrow(models_ars)),]
models_ars$dec_inc[1:5] = "decreasing"
models_ars$dec_inc[6:10] = "increasing"
models_ars = models_ars[,c(8,11)]
models_ars$COMP_ID = metab_metadata$COMP_ID[match(models_ars$chemical, metab_metadata$BIOCHEMICAL )]
metab_tab_ars = metab_tab_d[rownames(metab_tab_d) %in% models_ars$COMP_ID,]
metab_tab_ars["ID"] <- rownames(metab_tab_ars)
melted = data.table::melt(metab_tab_ars,id.vars	= "ID")
melted$inc_dec =  models_ars$dec_inc[match(melted$ID, models_ars$COMP_ID )]
melted$ID =  models_ars$chemical[match(melted$ID, models_ars$COMP_ID )]
melted$variable = mapped1$Arsenic[match(melted$variable, rownames(mapped1) )]
p3 <- ggplot(melted, aes(x=ID, y=value, fill=variable)) +
  facet_wrap(~inc_dec, scale="free") + theme_classic() +
  theme(axis.text.x = element_text(angle = 45 , hjust = 1),legend.position="top")+
  #geom_point(position=position_jitterdodge(), alpha = 0.15)+
  geom_boxplot(outlier.shape = NA,alpha=0.3) +
  xlab(NULL) +ggtitle("Arsenic")+ scale_fill_social() 
p3










####Interactions####
models_DG = models[models$treatment == "Diet:Genotype",]
models_DG = models_DG[order(models_DG$coefficients),]
write.table(models_DG, "models_DG.txt")


models_DG = models_DG[c(1:5, (nrow(models_DG)-4):nrow(models_DG)),]
models_DG$dec_inc[1:5] = "decreasing"
models_DG$dec_inc[6:10] = "increasing"
models_DG = models_DG[,c(8,11)]
models_DG$COMP_ID = metab_metadata$COMP_ID[match(models_DG$chemical, metab_metadata$BIOCHEMICAL )]

df = data.frame(Imidazole_propionate = as.numeric(metab_tab_d[rownames(metab_tab_d) == "15716",]),
                Diet = map$Diet, Arsenic = map$Arsenic, Genotype = map$Genotype)


p5 = ggline(df, x = "Genotype", y = "Imidazole_propionate", 
            add = c("mean_se", "jitter"),
            color = "Diet",shape="Diet", 
            palette =c("#1AFF1A", "#4B0092"), 
            ylab = "Normalized abundance", title = "Imidazole_propionate", xlab = "Genotype")
p5




###Interaction 2
models_DA =models[models$treatment == "Diet:Arsenic",]
models_DA = models_DA[order(models_DA$coefficients),]
write.table(models_DA, "models_DA.txt")

models_DA = models_DA[c(1:5, (nrow(models_DA)-4):nrow(models_DA)),]
models_DA$dec_inc[1:5] = "decreasing"
models_DA$dec_inc[6:10] = "increasing"
models_DA = models_DA[,c(8,11)]
models_DA$COMP_ID = metab_metadata$COMP_ID[match(models_DA$chemical, metab_metadata$BIOCHEMICAL )]

df = data.frame(sev_ketodeoxycholate = as.numeric(metab_tab_d[rownames(metab_tab_d) == models_DA$COMP_ID[10],]),
                Diet = map$Diet, Arsenic = map$Arsenic, Genotype = map$Genotype)



p6 = ggline(df, x = "Diet", y = "sev_ketodeoxycholate", 
            add = c("mean_se", "jitter"),
            color = "Arsenic",shape="Arsenic", 
            palette =c("#0077b4", "#ce201f"), 
            ylab = "Normalized abundance", title = "7-ketodeoxycholate", xlab = "Diet")
p6



###Interaction 3
models_GA = models[models$treatment == "Genotype:Arsenic",]
models_GA = models_GA[order(models_GA$coefficients),]
write.table(models_GA, "models_GA.txt")

models_GA = models_GA[c(1:5, (nrow(models_GA)-4):nrow(models_GA)),]
models_GA$dec_inc[1:5] = "decreasing"
models_GA$dec_inc[6:10] = "increasing"
models_GA = models_GA[,c(8,11)]
models_GA$COMP_ID = metab_metadata$COMP_ID[match(models_GA$chemical, metab_metadata$BIOCHEMICAL )]

df = data.frame(maltose = as.numeric(metab_tab_d[rownames(metab_tab_d) == models_GA$COMP_ID[1],]),
                Diet = map$Diet, Arsenic = map$Arsenic, Genotype = map$Genotype)



p7 = ggline(df, x = "Arsenic", y = "maltose", 
            add = c("mean_se", "jitter"),
            color = "Genotype",shape="Genotype", 
            palette =c("#00AFBB", "#E7B800"), 
            ylab = "Normalized abundance", title = "Maltose", xlab = "Arsenic")
p7


###Interaction 4
models_DAG = models[models$treatment == "Diet:Genotype:Arsenic",]
models_DAG = models_DAG[order(models_DAG$coefficients),]
write.table(models_DAG, "models_DAG.txt")


models_DAG = models_DAG[c(1:5, (nrow(models_DAG)-4):nrow(models_DAG)),]
models_DAG$dec_inc[1:5] = "decreasing"
models_DAG$dec_inc[6:10] = "increasing"
models_DAG = models_DAG[,c(8,11)]
models_DAG$COMP_ID = metab_metadata$COMP_ID[match(models_DAG$chemical, metab_metadata$BIOCHEMICAL )]

df = data.frame(TMAVA = as.numeric(metab_tab_d[rownames(metab_tab_d) == models_DAG$COMP_ID[5],]),
                Diet = map$Diet, Arsenic = map$Arsenic, Genotype = map$Genotype)


p8 = ggline(df, x = "Arsenic", y = "TMAVA", 
            add = c("mean_se", "jitter"),
            color = "Genotype",shape="Genotype", 
            palette =c("#40B0A6", "#E1BE6A"), 
            ylab = "Relative Abundance", title = "TMAVA", xlab = "As intake", facet.by = "Diet")
p8



  
  



