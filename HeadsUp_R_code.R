
#### prospective association analysis ####
library(tidyverse)
phenotype_data <- read.csv("C:/Users/Stephanie/OneDrive/S code/phenotype_data_pro.csv")

#load metabololomics data #
metabolomics_long0 <- read.csv("C:/Users/Stephanie/OneDrive/S code/Primary Metabolomics.Long.Form.csv")
metabolomics_long_y0 <- metabolomics_long0 %>%
  separate(treatment, c("group", "visit"), sep=" - ") %>%
  filter(visit == "Initial")
metabolomics_long_y1 <- metabolomics_long0 %>%
  separate(treatment, c("group", "visit"), sep=" - ") %>%
  filter(visit == "1 Year Later")

stopifnot(metabolomics_long_y0$label == metabolomics_long_y1$label)
stopifnot(metabolomics_long_y0$label == phenotype_data$Sample_Name)

#standardize the meta data
metabolomics_long_y0_data <- metabolomics_long_y0 %>%
  select(-c(`file id` : `BinBase.name`)) %>%
  log2() 

metabolomics_long_y1_data <- metabolomics_long_y1 %>%
  select(-c(`file id` : `BinBase.name`)) %>%
  log2()

metabolomics_long_y0_data <- as.data.frame(metabolomics_long_y0_data)
metabolomics_long_y1_data <- as.data.frame(metabolomics_long_y1_data)

pvalue_diffCor <- numeric(ncol(metabolomics_long_y0_data))
names(pvalue_diffCor) <- colnames(metabolomics_long_y0_data)
pvalue_diff <- pvalue_diffCor

pvalue_diffCor_BAND_IMI <- pvalue_diffCor
pvalue_diffCor_RYGB_IMI <- pvalue_diffCor
pvalue_diffCor_RYGB_BAND <- pvalue_diffCor


es_BAND_IMI <- pvalue_diff
es_RYGB_IMI <- pvalue_diff
es_RYGB_BAND <- pvalue_diff

se_BAND_IMI <- pvalue_diff
se_RYGB_IMI <- pvalue_diff
se_RYGB_BAND <- pvalue_diff

black.bold.text <- element_text(face = "bold", color = "black", size=20)

covs <- c("gluc", "hba1c")
types <- c("BAND_IMI", "BYPASS_IMI", "BYPASS_BAND")

contents <- c("p1_all", "p1_known", "p005_all", "p005_known", "q005_all", "q005_known")
summary_diff <- matrix(0,nrow=length(contents),ncol=length(covs))
rownames(summary_diff) <- contents
colnames(summary_diff) <- covs

summary_diff_list <- replicate(length(types), summary_diff, simplify=F)
names(summary_diff_list) <- types
phenotype_data$INSULIN <- ifelse(phenotype_data$INSULIN=="Yes", 2,
                                 ifelse(phenotype_data$INSULIN=="No", 1,NA))

phenotype_data$DMMED <- ifelse(phenotype_data$DMMED=="Yes", 2,
                               ifelse(phenotype_data$DMMED=="No", 1,NA))
table(phenotype_data$DM_drug)

## select variables 
#acov <- "gluc"
acov <- "hba1c"
for(acov in covs){
  cat("working on", acov,"\n")
  acov_y0 <- paste0(acov, "_v0")
  acov_y1 <- paste0(acov, "_v1")
  this_phenotype_data <- phenotype_data %>%
    mutate(value_base = phenotype_data[["hba1c_v0"]], diff=phenotype_data[["hba1c_v0"]] - phenotype_data[["hba1c_v1"]],
           BMI_v0 =  phenotype_data[["BMI_v0"]], 
           INSULIN =phenotype_data[["INSULIN"]], 
           DMMED=phenotype_data[["DMMED"]], 
           weightLoss = phenotype_data[["weight_v1"]] - phenotype_data[["weight_v0"]], 
           Reference2 = factor(Reference, levels=c("IMI", "BAND", "RYGB")), 
           Reference3 = factor(Reference, levels=c("RYGB", "BAND", "IMI"))) %>%
    select("Sample_Name", "Reference", "Reference2", "Reference3", "age", "sex", "value_base", 
           "weightLoss", "diff","BMI_v0","INSULIN","DMMED") 
}  
#table(this_phenotype_data$Reference)


## linear regrssion 
for(i in 1:ncol(metabolomics_long_y0_data)){
  aname <- colnames(metabolomics_long_y0_data)[i]
  ameta <- metabolomics_long_y0_data[[i]]
  that_phenotype_data <- this_phenotype_data %>%
    mutate(metabolite=ameta)
  
  that_phenotype_data$age <- scale(that_phenotype_data$age)
  that_phenotype_data$BMI_v0 <-scale(that_phenotype_data$BMI_v0,center = TRUE, scale = TRUE) 
  that_phenotype_data$weightLoss <- scale(that_phenotype_data$weightLoss,center = TRUE, scale = TRUE)
  that_phenotype_data$value_base <- scale(that_phenotype_data$value_base,center = TRUE, scale = TRUE)
  alm <- lm(diff ~ metabolite * Reference + age + sex + value_base + BMI_v0+weightLoss+metabolite+Reference, data=that_phenotype_data)
  blm <- lm(diff ~ metabolite * Reference2 + age + sex + value_base + BMI_v0+weightLoss+metabolite+Reference2, data=that_phenotype_data)
  
  asumm <- summary(alm)
  bsumm <- summary(blm)
  
  es_BAND_IMI[aname] <- asumm$coefficients["metabolite:ReferenceIMI", "Estimate"]
  es_RYGB_BAND[aname] <- bsumm$coefficients["metabolite:Reference2RYGB", "Estimate"]
  es_RYGB_IMI[aname] <- asumm$coefficients["metabolite:ReferenceRYGB", "Estimate"]
  
  se_BAND_IMI[aname] <- asumm$coefficients["metabolite:ReferenceIMI", "Std. Error"]
  se_RYGB_IMI[aname] <- bsumm$coefficients["metabolite:Reference2RYGB", "Std. Error"]
  se_RYGB_BAND[aname] <- asumm$coefficients["metabolite:ReferenceRYGB", "Std. Error"]
  
  pvalue_diffCor_BAND_IMI[aname] <- asumm$coefficients["metabolite:ReferenceIMI", "Pr(>|t|)"]
  pvalue_diffCor_RYGB_IMI[aname] <- bsumm$coefficients["metabolite:Reference2RYGB", "Pr(>|t|)"]
  pvalue_diffCor_RYGB_BAND[aname] <- asumm$coefficients["metabolite:ReferenceRYGB", "Pr(>|t|)"]
  print(i)	
}


qvalue_diffCor_BAND_IMI <- p.adjust(pvalue_diffCor_BAND_IMI, "BH")
qvalue_diffCor_RYGB_IMI <- p.adjust(pvalue_diffCor_RYGB_IMI, "BH")
qvalue_diffCor_RYGB_BAND <- p.adjust(pvalue_diffCor_RYGB_BAND, "BH")

table(qvalue_diffCor_RYGB_IMI < 0.05)
table(qvalue_diffCor_BAND_IMI < 0.05)
table(qvalue_diffCor_RYGB_BAND < 0.05)


setwd(WD_result)
result_diffCor <- tibble(
  metabolite=colnames(metabolomics_long_y0_data), 
  es_BAND_IMI = es_BAND_IMI, 
  es_RYGB_IMI = es_RYGB_IMI, 
  es_RYGB_BAND = es_RYGB_BAND, 
  se_BAND_IMI = se_BAND_IMI, 
  se_RYGB_IMI_ = se_RYGB_IMI, 
  se_RYGB_BAND = se_RYGB_BAND, 
  
  pvalue_BAND_IMI=pvalue_diffCor_BAND_IMI,
  qvalue_BAND_IMI=qvalue_diffCor_BAND_IMI, 
  pvalue_RYGB_IMI=pvalue_diffCor_RYGB_IMI,
  qvalue_RYGB_IMI=qvalue_diffCor_RYGB_IMI,
  pvalue_RYGB_BAND=pvalue_diffCor_RYGB_BAND,
  qvalue_RYGB_BAND=qvalue_diffCor_RYGB_BAND
) %>%
  arrange(pvalue_BAND_IMI)
write.csv(result_diffCor,"C:/Users/Stephanie/OneDrive/S code/gluc_prospective_new.csv")

#### Repeated measurement analysis ####
phenotype_data <- read.csv("C:/Users/Stephanie/OneDrive/S code/phenotype_data_pro.csv")

metabolomics_long0 <- read.csv("C:/Users/Stephanie/OneDrive/S code/Primary Metabolomics.Long.Form.csv")
metabolomics_long_y0 <- metabolomics_long0 %>%
  separate(treatment, c("group", "visit"), sep=" - ") %>%
  filter(visit == "Initial")
metabolomics_long_y1 <- metabolomics_long0 %>%
  separate(treatment, c("group", "visit"), sep=" - ") %>%
  filter(visit == "1 Year Later")

stopifnot(metabolomics_long_y0$label == metabolomics_long_y1$label)
stopifnot(metabolomics_long_y0$label == phenotype_data$Sample_Name)

#standardize the lipid data
metabolomics_long_y0_data <- metabolomics_long_y0 %>%
  select(-c(`file id` : `BinBase.name`)) %>%
  log2() 

metabolomics_long_y0_data <- as.data.frame(metabolomics_long_y0_data)
metabolomics_long_y1_data <- metabolomics_long_y1 %>%
  select(-c(`file id` : `BinBase.name`)) %>%
  log2() 


metabolomics_long_y0_data <- as.data.frame(metabolomics_long_y0_data)
metabolomics_long_y1_data <- as.data.frame(metabolomics_long_y1_data)  

if(F){
  substring(colnames(metabolomics_long_y0_data),1,1) <- toupper(substring(colnames(metabolomics_long_y0_data),1,1))
  substring(colnames(metabolomics_long_y1_data),1,1) <- toupper(substring(colnames(metabolomics_long_y1_data),1,1))	
}


types <- c("BAND_IMI", "RYGB_IMI", "RYGB_BAND")
types2 <- c("IMI", "BAND", "RYGB")

contents <- c("p1_all", "p1_known", "p005_all", "p005_known", "q005_all", "q005_known")
summary_diffChangeByGroup <- matrix(0,nrow=length(contents),ncol=length(types))
rownames(summary_diffChangeByGroup) <- contents
summary_visit <- summary_diffChangeByGroup
colnames(summary_diffChangeByGroup) <- types
colnames(summary_visit) <- types2


# covs <- c("BMI", "weight", "gluc", "hba1c")
## select variables 
this_phenotype_data <- phenotype_data %>%
  mutate(Reference2 = factor(Reference, levels=c("IMI", "BAND", "RYGB"))) %>%
  mutate(Reference3 = factor(Reference, levels=c("RYGB", "BAND", "IMI"))) %>%
  select("Sample_Name", "Reference", "Reference2", "Reference3", "age", "sex", 
         "BMI_v0","BMI_v1", "weight_v0", "weight_v1","gluc_v0", "gluc_v1","hba1c_v0", "hba1c_v1")

pvalue_diff <- numeric(ncol(metabolomics_long_y0_data))
names(pvalue_diff) <- colnames(metabolomics_long_y0_data)
pvalue_diff_BAND_IMI <- pvalue_diff
pvalue_diff_RYGB_IMI <- pvalue_diff
pvalue_diff_RYGB_BAND <- pvalue_diff

es_diff_BAND_IMI <- pvalue_diff
es_diff_RYGB_IMI <- pvalue_diff
es_diff_RYGB_BAND <- pvalue_diff

se_diff_BAND_IMI<- pvalue_diff
se_diff_RYGB_IMI <- pvalue_diff
se_diff_RYGB_BAND <- pvalue_diff

black.bold.text <- element_text(face = "bold", color = "black", size=20)

covs <- c("gluc", "hba1c")
types <- c("BAND_IMI", "RYGB_IMI", "RYGB_BAND")

## linear mixed models
for(i in 1:ncol(metabolomics_long_y0_data)){
  aname <- colnames(metabolomics_long_y0_data)[i]
  stopifnot(aname == colnames(metabolomics_long_y1_data)[i])
  
  ameta_y0 <- metabolomics_long_y0_data[[i]]   
  ameta_y1 <- metabolomics_long_y1_data[[i]]   
  
  that_phenotype_data <- this_phenotype_data %>% 
    mutate(base_meta = ameta_y0, meta_y0 = ameta_y0, meta_y1=ameta_y1)
  
  that_weight_data_long <- that_phenotype_data %>%
    gather(visit0,weight, c(weight_v0, weight_v1)) %>%
    separate(col = visit0, into = c("none", "visit"), sep = "_v")%>%
    mutate(visit=as.numeric(visit))%>% select(-none)
  
  that_BMI_data_long <- that_phenotype_data %>%
    gather(visit0,BMI, c(BMI_v0, BMI_v1)) %>%
    separate(col = visit0, into = c("none", "visit"), sep = "_v")%>%
    mutate(visit=as.numeric(visit))%>% select(-none)
  
  #that_glucose_data_long <- that_phenotype_data %>%
  #gather(visit0,gluc,c(gluc_v0, gluc_v1)) %>%
  #separate(col = visit0, into = c("none", "visit"), sep = "_v")%>%
  #mutate(visit=as.numeric(visit))%>% select(-none)
  
  that_hba1c_data_long <- that_phenotype_data %>%
    gather(visit0,hba1c,c(hba1c_v0, hba1c_v1)) %>%
    separate(col = visit0, into = c("none", "visit"), sep = "_v")%>%
    mutate(visit=as.numeric(visit))%>% select(-none)
  
  that_phenotype_data_WB <- full_join(that_weight_data_long, that_BMI_data_long)
  that_phenotype_data_short <- full_join(that_phenotype_data_WB,that_hba1c_data_long)
  
  that_meta_data_long <- that_phenotype_data %>%
    gather(visit0,metabolite,c(meta_y0, meta_y1)) %>%
    separate(col = visit0, into = c("none", "visit"), sep = "_y") %>%
    mutate(visit=as.numeric(visit)) %>% select(-none)
  
  that_phenotype_data_long <- full_join(that_meta_data_long, that_phenotype_data_short)
  
  #almer <- lmer(gluc ~ Reference2 *metabolite + age + sex + weight + BMI +visit+ (1|Sample_Name), data = that_phenotype_data_long)
  #blmer <- lmer(gluc ~ Reference  *metabolite + age + sex + weight + BMI +visit+(1|Sample_Name), data = that_phenotype_data_long)
  
  almer <- lmer(hba1c ~ Reference2 *metabolite + age + sex + weight + BMI +visit+ (1|Sample_Name), data = that_phenotype_data_long)
  blmer <- lmer(hba1c ~ Reference  *metabolite + age + sex + weight + BMI +visit+(1|Sample_Name), data = that_phenotype_data_long)
  
  
  <- summary(almer)
  bsumm <- summary(blmer)
  
  
  es_diff_BAND_IMI[aname] <- asumm$coefficients["Reference2BAND:metabolite", "Estimate"]
  es_diff_RYGB_IMI[aname] <- asumm$coefficients["Reference2RYGB:metabolite", "Estimate"]
  es_diff_RYGB_BAND[aname] <- bsumm$coefficients["ReferenceRYGB:metabolite", "Estimate"]
  
  se_diff_BAND_IMI[aname] <- asumm$coefficients["Reference2BAND:metabolite", "Std. Error"]
  se_diff_RYGB_IMI[aname] <- asumm$coefficients["Reference2RYGB:metabolite", "Std. Error"]
  se_diff_RYGB_BAND[aname] <- bsumm$coefficients["ReferenceRYGB:metabolite", "Std. Error"]
  
  pvalue_diff_BAND_IMI[aname] <- asumm$coefficients["Reference2BAND:metabolite", "Pr(>|t|)"]
  pvalue_diff_RYGB_IMI[aname] <- asumm$coefficients["Reference2RYGB:metabolite", "Pr(>|t|)"]
  pvalue_diff_RYGB_BAND[aname] <- bsumm$coefficients["ReferenceRYGB:metabolite", "Pr(>|t|)"]
  print(i)
}


qvalue_diff_BAND_IMI <- p.adjust(pvalue_diff_BAND_IMI, "BH")
qvalue_diff_RYGB_IMI <- p.adjust(pvalue_diff_RYGB_IMI, "BH")
qvalue_diff_RYGB_BAND <- p.adjust(pvalue_diff_RYGB_BAND, "BH")
table(qvalue_diff_BAND_IMI < 0.05)
table(qvalue_diff_RYGB_IMI < 0.05)
table(qvalue_diff_RYGB_BAND < 0.05)


result_diff <- tibble(metabolite=colnames(metabolomics_long_y0_data), 
                      es_diff_BAND_IMI = es_diff_BAND_IMI, 
                      se_diff_BAND_IMI = se_diff_BAND_IMI, 
                      es_diff_RYGB_IMI = es_diff_RYGB_IMI, 
                      se_diff_RYGB_IMI = se_diff_RYGB_IMI, 
                      es_diff_RYGB_BAND= es_diff_RYGB_BAND, 
                      se_diff_RYGB_BAND = se_diff_RYGB_BAND, 
                      pvalue_BAND_IMI=pvalue_diff_BAND_IMI,
                      pvalue_RYGB_IMI=pvalue_diff_RYGB_IMI,
                      pvalue_RYGB_BAND=pvalue_diff_RYGB_BAND,
                      qvalue_RYGB_IMI=qvalue_diff_RYGB_IMI,
                      qvalue_BAND_IMI=qvalue_diff_BAND_IMI, 
                      qvalue_RYGB_BAND=qvalue_diff_RYGB_BAND					  
) %>%
  arrange(pvalue_diff_RYGB_IMI)
afile <- paste0("result_all.csv")
result_diff2 <- round(result_diff[,c(2:7)],2)
result_diff[,c(2:7)] <- result_diff2[,c(1:6)]

write.csv(result_diff,"C:/Users/Stephanie/OneDrive/S code/hba1c_longitudinal_new.csv")


#### network analysis 
library(WGCNA)
#require(flashClust)
#library(magrittr)


#bsize=5000
#sft <- pickSoftThreshold(bypass,blockSize=bsize,dataIsExpr = TRUE,corFnc = bicor,networkType = "unsigned")
#net = blockwiseModules(bypass, maxBlockSize=25000,power = 4, corFnc = bicor,maxPOutliers = 0.05,
#networkType = "unsigned", minModuleSize = 30, # 20
#reassignThreshold = 0, mergeCutHeight = 0.1,deepSplit=4,
#numericLabels = TRUE, pamRespectsDendro = FALSE,saveTOMs = F,verbose = 3)
#table(net$colors)

TOM=TOMsimilarityFromExpr(bypass, corType='bicor',maxPOutliers=0.05, networkType='unsigned', power = 5)
diss=1-TOM
Tree= hclust(as.dist(diss),method="average")
diag(diss) = NA
#cutreeDynamic(dendro = Tree, distM = as.matrix(1-TOM),deepSplit = 4, pamRespectsDendro = FALSE,minClusterSize = 50)
tree = cutreeHybrid(dendro = Tree, pamStage=FALSE,minClusterSize = 30,# cutHeight = 0.9, ## the maximum joining heights for dendrogram
                    deepSplit = 4, distM = as.matrix(1-TOM))

mergedColors = labels2colors(tree$labels)
# Calculate eigengenes
MEList = moduleEigengenes(bypass, colors = mergedColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs)
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average")
# Call an automatic merging function
merge = mergeCloseModules(bypass, mergedColors, cutHeight = 0.15, verbose = 3)
# The merged module colors
mergedColors = merge$colors
table(mergedColors)
modules=unique(mergedColors)
mergedMEs = merge$newMEs
plotEigengeneNetworks(mergedMEs, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2))
gene_membership <- data.frame(mergedColors )

names(gene_membership) <- "Metabolite group"


sizeGrWindow(10,5)
plotDendroAndColors(Tree,
                    colors=gene_membership,
                    dendroLabels = FALSE, hang = 0.02,
                    addGuide = TRUE, guideHang = 0.01,
                    main="Cluster dendrogram of metabolites")
# save module info
write.table(mergedMEs,'network_eigenvalues_bypass.csv',sep=',',row.names = F, col.names=T) 
dt0=cbind(colnames(bypass),mergedColors)
colnames(dt0)=c('metabolite','module')
dt0=dt0[order(dt0[,'module'],decreasing = F),]
write.table(dt0,'log2_network_bypass.csv',sep=',',row.names = F, col.names=T)


######################################IMI network
TOM=TOMsimilarityFromExpr(imi, corType='bicor',maxPOutliers=0.05, networkType='unsigned', power = 5)
diss=1-TOM
Tree= hclust(as.dist(diss),method="average")
diag(diss) = NA
#cutreeDynamic(dendro = Tree, distM = as.matrix(1-TOM),deepSplit = 4, pamRespectsDendro = FALSE,minClusterSize = 50)
tree = cutreeHybrid(dendro = Tree, pamStage=FALSE,minClusterSize = 30,# cutHeight = 0.9, ## the maximum joining heights for dendrogram
                    deepSplit = 4, distM = as.matrix(1-TOM))

mergedColors = labels2colors(tree$labels)
# Calculate eigengenes
MEList = moduleEigengenes(imi, colors = mergedColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs)
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average")
# Call an automatic merging function
merge = mergeCloseModules(imi, mergedColors, cutHeight = 0.15, verbose = 3)
# The merged module colors
mergedColors = merge$colors
table(mergedColors)
modules=unique(mergedColors)
mergedMEs = merge$newMEs
plotEigengeneNetworks(mergedMEs, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2))
gene_membership <- data.frame(mergedColors )

names(gene_membership) <- "Metabolite group"


sizeGrWindow(10,5)
plotDendroAndColors(Tree,
                    colors=gene_membership,
                    dendroLabels = FALSE, hang = 0.02,
                    addGuide = TRUE, guideHang = 0.01,
                    main="Cluster dendrogram of metabolites")
# save module info
write.table(mergedMEs,'log2_network_eigenvalues_imi.csv',sep=',',row.names = F, col.names=T) 
dt0=cbind(colnames(imi),mergedColors)
colnames(dt0)=c('metabolite','module')
dt0=dt0[order(dt0[,'module'],decreasing = F),]
write.table(dt0,'log2_network_imi.csv',sep=',',row.names = F, col.names=T)


######################################BAND network
TOM=TOMsimilarityFromExpr(band, corType='bicor',maxPOutliers=0.05, networkType='unsigned', power = 5)
diss=1-TOM
Tree= hclust(as.dist(diss),method="average")
diag(diss) = NA
#cutreeDynamic(dendro = Tree, distM = as.matrix(1-TOM),deepSplit = 4, pamRespectsDendro = FALSE,minClusterSize = 50)
tree = cutreeHybrid(dendro = Tree, pamStage=FALSE,minClusterSize = 30,# cutHeight = 0.9, ## the maximum joining heights for dendrogram
                    deepSplit = 4, distM = as.matrix(1-TOM))

mergedColors = labels2colors(tree$labels)
# Calculate eigengenes
MEList = moduleEigengenes(band, colors = mergedColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs)
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average")
# Call an automatic merging function
merge = mergeCloseModules(band, mergedColors, cutHeight = 0.15, verbose = 3)
# The merged module colors
mergedColors = merge$colors
table(mergedColors)
modules=unique(mergedColors)
mergedMEs = merge$newMEs
gene_membership <- data.frame(mergedColors )

names(gene_membership) <- "Metabolite group"

sizeGrWindow(10,5)
plotDendroAndColors(Tree,
                    colors=gene_membership,
                    dendroLabels = FALSE, hang = 0.02,
                    addGuide = TRUE, guideHang = 0.01,
                    main="Cluster dendrogram of metabolites")

# save module info
write.table(mergedMEs,'log2_network_eigenvalues_band.csv',sep=',',row.names = F, col.names=T) 
dt0=cbind(colnames(band),mergedColors)
colnames(dt0)=c('metabolite','module')
dt0=dt0[order(dt0[,'module'],decreasing = F),]
write.table(dt0,'log2_network_band.csv',sep=',',row.names = F, col.names=T)


#Calculate connectivity RYGB VS IMI
network_par<-function(dt1,dt2,module){
  lipids_sub1=dt1[,mergedColors==module]
  lipids_sub2=dt2[,mergedColors==module]
  power=5 #both power 5
  adj1=adjacency(lipids_sub1, power = power, type = "unsigned",corFnc = 'bicor') # "signed hybrid"
  adj2=adjacency(lipids_sub2, power = power, type = "unsigned",corFnc = 'bicor') # "signed hybrid"
  diag(adj1)=0
  diag(adj2)=0
  ### Fundamental Network Concepts
  Size=dim(adj1)[1]
  Connectivity1=apply(adj1, 2, sum) # Within Module Connectivities
  Connectivity2=apply(adj2, 2, sum) 
  Connectivity=sum(Connectivity1)-sum(Connectivity2) #cor(Connectivity1,Connectivity2)
  Density1=sum(Connectivity1)/(Size*(Size-1))
  Density2=sum(Connectivity2)/(Size*(Size-1))
  Density=Density1-Density2
  vect_1=c(Size,Connectivity,Density)
  return(vect_1)
}

data_1=rbind(bypass,imi)
iter=1000
ot=NULL
for (i in 1:length(modules)){
  vect=NULL
  for (j in 1:iter){
    set.seed(j)
    case_id=sample(seq(1,100,1),size=50)
    #part1=sample(depression[depression[,'depression_progress']==1,'ID'],size=100)
    #control_id=c(part1,sample(depression[!depression[,'ID']%in%part1,'ID'],size=(nrow(lipids_control)-100)))
    dt1=data_1[case_id,]
    dt2=data_1[!seq(1,100,1)%in%case_id,]
    net_par=network_par(dt1,dt2,modules[i])
    vect=rbind(vect,c(modules[i],net_par))
  }
  #rownames(vect)=rep(modules[i],iter)
  ot=rbind(ot,vect)
  print(i)
}
colnames(ot)=c("Module","Size","Connectivity","Density")
write.table(ot,'log2_network comparison_bypass_imi.csv',sep=',',row.names = F, col.names=T)



#####
data <- read.csv("adj_residual.csv")
col_name=scan('adj_residual.csv', what = "",sep = ",",nlines = 1, quiet = TRUE, skip = 0, strip.white = TRUE)
colnames(data)=col_name
table(data[,'group'])

band=data[data[,'group']=='Band',2:ncol(data)]
bypass=data[data[,'group']=='Bypass',2:ncol(data)]
imi=data[data[,'group']=='IMI',2:ncol(data)]

library(WGCNA)

#######imi
TOM=TOMsimilarityFromExpr(imi, corType='bicor',maxPOutliers=0.05, networkType='unsigned', power = 5)
diss=1-TOM
Tree= hclust(as.dist(diss),method="average")
diag(diss) = NA
#cutreeDynamic(dendro = Tree, distM = as.matrix(1-TOM),deepSplit = 4, pamRespectsDendro = FALSE,minClusterSize = 50)
tree = cutreeHybrid(dendro = Tree, pamStage=FALSE,minClusterSize = 30,# cutHeight = 0.9, ## the maximum joining heights for dendrogram
                    deepSplit = 4, distM = as.matrix(1-TOM))

mergedColors = labels2colors(tree$labels)
# Calculate eigengenes
MEList = moduleEigengenes(imi, colors = mergedColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs)
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average")
# Call an automatic merging function
merge = mergeCloseModules(imi, mergedColors, cutHeight = 0.15, verbose = 3)
# The merged module colors
mergedColors = merge$colors
table(mergedColors)
modules=unique(mergedColors)
chooseTopHubInEachModule(imi,mergedColors,power=5)

# hub lipid
lipids_select1=imi[,1:153]
mergedColors1=mergedColors[1:153]
chooseTopHubInEachModule(lipids_select1,mergedColors1,power=5) # hub gene in each module


# hub lipids in each lipid module
restGenes= (mergedColors == "turquoise") 
adj = adjacency(bypass[,restGenes], power = 5,type = "unsigned",corFnc = 'bicor')
Connectivity=apply(adj, 2, sum)
sort(Connectivity,decreasing=T)
# delete unknown
restGenes= (mergedColors1 == "turquoise")
adj = adjacency(lipids_select1[,restGenes], power = 5,type = "unsigned",corFnc = 'bicor')
Connectivity=apply(adj, 2, sum)
sort(Connectivity,decreasing=T)


###bypass
TOM=TOMsimilarityFromExpr(bypass, corType='bicor',maxPOutliers=0.05, networkType='unsigned', power = 5)
diss=1-TOM
Tree= hclust(as.dist(diss),method="average")
diag(diss) = NA
#cutreeDynamic(dendro = Tree, distM = as.matrix(1-TOM),deepSplit = 4, pamRespectsDendro = FALSE,minClusterSize = 50)
tree = cutreeHybrid(dendro = Tree, pamStage=FALSE,minClusterSize = 30,# cutHeight = 0.9, ## the maximum joining heights for dendrogram
                    deepSplit = 4, distM = as.matrix(1-TOM))

mergedColors = labels2colors(tree$labels)
# Calculate eigengenes
MEList = moduleEigengenes(bypass, colors = mergedColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs)
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average")
# Call an automatic merging function
merge = mergeCloseModules(bypass, mergedColors, cutHeight = 0.15, verbose = 3)
# The merged module colors
mergedColors = merge$colors
table(mergedColors)
modules=unique(mergedColors)
chooseTopHubInEachModule(bypass,mergedColors,power=5)

# hub lipid
lipids_select1=bypass[,1:153]
mergedColors1=mergedColors[1:153]
chooseTopHubInEachModule(lipids_select1,mergedColors1,power=5) # hub gene in each module


network_par<-function(dt1,dt2){
  vect=NULL
  for (i in 1:length(modules)){
    lipids_sub1=dt1[,mergedColors==modules[i]]
    lipids_sub2=dt2[,mergedColors==modules[i]]
    power=5
    adj1=adjacency(lipids_sub1, power = power, type = "unsigned",corFnc = 'bicor') # "signed hybrid"
    adj2=adjacency(lipids_sub2, power = power, type = "unsigned",corFnc = 'bicor') # "signed hybrid"
    diag(adj1)=0
    diag(adj2)=0
    ### Fundamental Network Concepts
    Size=dim(adj1)[1]
    Connectivity1=apply(adj1, 2, sum) # Within Module Connectivities
    Connectivity2=apply(adj2, 2, sum) 
    Connectivity=sum(Connectivity1)-sum(Connectivity2) #cor(Connectivity1,Connectivity2)
    Density1=sum(Connectivity1)/(Size*(Size-1))
    Density2=sum(Connectivity2)/(Size*(Size-1))
    Density=Density1-Density2
    vect_1=c(Size,Connectivity,Density)
    vect=rbind(vect,vect_1)
  }
  return(vect)
}

net_par=network_par(bypass,imi)#net_par=network_par(band,imi)
rownames(net_par)=modules
sum(net_par[,1]) # 364


#### analyze permutation results
net_0=net_par
net_0=cbind(modules,net_0)
network=read.csv('log2_network comparison_bypass_imi.csv',sep=',',header=T) 
#network=read.csv('log2_network comparison_bypass_band.csv',sep=',',header=T) 
# case female vs. case male: network comparison_female_male.csv  network comparison_lipids_female_male.csv
# control female vs. case male: network comparison_female_male_control.csv  network comparison_lipids_female_male_control.csv
# female case vs. control: network comparison_female.csv network comparison_lipids_female.csv
# male case vs. control: network comparison_male.csv network comparison_lipids_male.csv
table(network[,1])
colnames(net_0)=colnames(network)
head(network)
network[,'Density']=abs(as.numeric(network[,'Density'])) 
network[,'Connectivity']=abs(as.numeric(network[,'Connectivity'])) #as.numeric(network[,'Connectivity']) 
net_0[,'Density']=abs(as.numeric(net_0[,'Density'])) 
net_0[,'Connectivity']=abs(as.numeric(net_0[,'Connectivity'])) #as.numeric(net_0[,'Connectivity']) 

# H0: case and control networks have similar structures: null connectivity correlation close to 0
FDR=matrix(0,nrow(net_0),(ncol(net_0)-2)) # each row is for each modules
for (i in 1:nrow(net_0)){
  FDR[i,1]=sum(as.numeric(network[network[,'Module']==net_0[i,1],'Connectivity'])>as.numeric(net_0[i,'Connectivity']))/1000
  #FDR[i,1]=ifelse(as.numeric(net_0[i,'Connectivity'])>0,sum(as.numeric(network[network[,'Module']==net_0[i,1],'Connectivity'])>as.numeric(net_0[i,'Connectivity']))/1000,
  #sum(as.numeric(network[network[,'Module']==net_0[i,1],'Connectivity'])<as.numeric(net_0[i,'Connectivity']))/1000)  
  FDR[i,2]=sum(as.numeric(network[network[,'Module']==net_0[i,1],'Density'])>as.numeric(net_0[i,'Density']))/1000
}
rownames(FDR)=modules
colnames(FDR)=colnames(net_0)[3:4]

new=cbind(modules,net_par,rep(0,nrow(net_par)))
colnames(new)=c(colnames(network),'cat')
new[,'cat']=''
new[FDR[,'Connectivity']<0.05&as.numeric(new[,'Connectivity'])>0,'cat']='GOC'   # P<0.05 or P<0.1
new[FDR[,'Connectivity']<0.05&as.numeric(new[,'Connectivity'])<0,'cat']='LOC'
new[FDR[,'Connectivity']>=0.05,'cat']='Conserved'

library(ggplot2)
require("ggrepel")
new=data.frame(new)
sizeGrWindow(6,6)
new[,'Size']=as.numeric(new[,'Size'])
new[,'Connectivity']=as.numeric(new[,'Connectivity'])
ggplot(new, aes(Size, Connectivity, color = cat)) + labs(x = "Module size",y='MDC')+geom_point(size = 3)  +
  scale_color_manual(values = c("#388ECC","#F68B33",'#00BFC4'))+theme_minimal()+
  theme(legend.text = element_text(size=12),axis.title.x = element_text( size=12),axis.title.y = element_text(size=12),axis.text.y = element_text(size=10),
        axis.text.x = element_text(size=10),legend.position = "right",panel.background = element_rect(fill = "white"),
        panel.grid.minor = element_blank(),panel.grid.major = element_blank())+ guides(color=guide_legend(title=""))+
  geom_text_repel(aes(label = Module),min.segment.length = 0,size=3,box.padding = 0.4,position = position_jitter(seed = 10),max.overlaps = Inf) # ,max.overlaps = 8


################ plot showing the overlap & difference between case modules & control modules
network_case=read.table('log2_network_bypass.csv',sep=',',header=T)
data.frame(table(network_case[,2]))
#network_control=read.table('log2_network_imi.csv',sep=',',header=T)
network_control=read.table('log2_network_band.csv',sep=',',header=T)
data.frame(table(network_control[,2]))

network_control=network_control[match(network_case[,1],network_control[,1]),]
sum(network_case[,1]==network_control[,1])
###### create dataset for the plot 1: module
dt=data.frame(case=network_case[,'module'],control=network_control[,'module'])
library(plyr)
library(dplyr)
counts <- ddply(dt, .(dt$case, dt$control), nrow)
names(counts) <- c("case", "control", "Freq")
count=counts[order(counts[,'case']),]
count[,'case']=paste0(count[,'case'],' (RYGB)')   # paste0(count[,'case'],'1')
count[,'control']=paste0(count[,'control'],' (BAND)') # paste0(count[,'control'],'0')
#count[,'control']=paste0(count[,'control'],' (IMI)') # paste0(count[,'control'],'0')
count[,'case']=paste(toupper(substr(count[,'case'], 1, 1)), substr(count[,'case'], 2, nchar(count[,'case'])), sep="")
count[,'control']=paste(toupper(substr(count[,'control'], 1, 1)), substr(count[,'control'], 2, nchar(count[,'control'])), sep="")
library(circlize)
label=unique(c(count[,1],count[,2]))
grid.col=sub("\\ .*", "", label)
library(stringr)
grid.col=str_to_sentence(grid.col)
names(grid.col)=label
#grid.col = c(Purple1='Purple',Red1='Red',Turquoise1='Turquoise',Blue1='Blue',Magenta1='Magenta',Black1='Black',
#Yellow1='Yellow',Brown1='Brown',Green1='Green',Pink1='Pink',Grey1='Grey',Purple0='Purple',Red0='Red',Turquoise0='Turquoise',Blue0='Blue',Magenta0='Magenta',Black0='Black',
#Yellow0='Yellow',Brown0='Brown',Green0='Green',Pink0='Pink',Grey0='Grey')
sizeGrWindow(6,5)
chordDiagram(count, annotationTrack = "grid", preAllocateTracks = 1,grid.col = grid.col, transparency = 0.35)
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, # cex=1.4,
              facing = "clockwise", niceFacing = TRUE, adj = c(-0.1, 0.5))}, bg.border = NA)
circos.clear()
dev.off()
# figure dimension: 927*850

############################ hypergeometric test for all pair-wise module comparison between case network & control network
network_case=read.table('log2_network_bypass.csv',sep=',',header=T)
data.frame(table(network_case[,2]))
network_control=read.table('log2_network_imi.csv',sep=',',header=T)
network_control=read.table('log2_network_band.csv',sep=',',header=T)
data.frame(table(network_control[,2]))

module_case=unique(network_case[,2])
n_case=length(module_case) # No. of case modules
module_control=unique(network_control[,2])
n_control=length(module_control) # No. of case modules

P=matrix(NA,n_case,n_control)
n_overlap=matrix(0,n_case,n_control)
for (i in 1:n_case){
  for (j in 1:n_control){
    n=sum(network_case[network_case[,2]==module_case[i],1]%in%network_control[network_control[,2]==module_control[j],1])
    n1=length(network_case[network_case[,2]==module_case[i],1])
    n2=length(network_control[network_control[,2]==module_control[j],1])
    n_overlap[i,j]=n
    if (n>0){
      p=fisher.test(matrix(c(n,n1-n,n2-n,nrow(network_case)-n1-n2+n),nr = 2))$p.value
      P[i,j]=p
    }
  }
}
library(qvalue)
q_1=qvalue(as.numeric(P[!is.na(P)]),lambda=0)$qvalues
q=as.numeric(P)
q[!is.na(P)]=q_1
q[is.na(P)]=1
q
####### adjust q value cut-off threshold to make sure each module only has at most 2 boxes with color "green"
textMatrix=ifelse(q<1e-20,'***','')
colorMatrix=ifelse(q<1e-20,'lightgreen','white')
textMatrix[textMatrix!='***']=ifelse(q[textMatrix!='***']<1e-10,'**','')
colorMatrix[colorMatrix!='lightgreen']=ifelse(q[colorMatrix!='lightgreen']<1e-10,'lightskyblue','white')
textMatrix[!textMatrix%in%c('***','**')]=ifelse(q[!textMatrix%in%c('***','**')]<0.00001,'*','')
colorMatrix[!colorMatrix%in%c('lightgreen','lightskyblue')]=ifelse(q[!colorMatrix%in%c('lightgreen','lightskyblue')]<0.00001,'orange','white')
textMatrix[!textMatrix%in%c('***','**','*')]=ifelse(q[!textMatrix%in%c('***','**','*')]<0.05,'+','')
colorMatrix[!colorMatrix%in%c('lightgreen','lightskyblue','orange')]=ifelse(q[!colorMatrix%in%c('lightgreen','lightskyblue','orange')]<0.05,'pink','white')
colorMatrix[!colorMatrix%in%c('lightgreen','lightskyblue','orange','pink')]=ifelse(n_overlap[!colorMatrix%in%c('lightgreen','lightskyblue','orange','pink')]>0,'grey','white')
q=matrix(q,n_case,n_control)
textMatrix=matrix(textMatrix,n_case,n_control)
FDR=-log10(q)
library(WGCNA)
# Extend margins to fit all labels
par(mar = c(6, 6, 4, 6)) # dimension: 700*600

labeledHeatmap(Matrix = FDR, yLabels = paste0('ME',module_case),ySymbols = module_case, xLabels = paste0('ME',module_control), 
               xSymbols = module_control, colorLabels = FALSE, colorMatrix=colorMatrix,#colors = blueWhiteRed(50), 
               textMatrix = textMatrix, setStdMargins = FALSE, cex.text = 1.4, zlim = c(0,90), 
               main = paste(""),yLabelsPosition = "left",xLabelsPosition = "bottom",
               xLabelsAdj = 0.7,x.adj.lab.y=0.9,xColorOffset =0.012, plotLegend = F,
               yColorOffset = 0.01,cex.lab.y = 1,xColorWidth =  strheight("M")-0.02,
               yColorWidth = strheight("M")-0.026)

