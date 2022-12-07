############################################################
############ Written by Jennifer Currenti ##################
############################################################
#                                                          #
#                          Libraries                       #
#                                                          #
############################################################

library(Seurat)
library(ggplot2)
library(tidyverse)
library(ggsankey)

############################################################
#                                                          #
#                   Genomic data - Fig 1                   #
#                                                          #
############################################################

#Sankey plot
dat.sankey <- data.frame(Sample2 = c("other_genomic","other_genomic",
                                     "targeted_genomic_interval","targeted_genomic_interval",
                                     "percentage_total_transcriptome"),
                         Sample2_depleted = c("other_genomic","percentage_total_transcriptome",
                                              "targeted_genomic_interval","percentage_total_transcriptome",
                                              "percentage_total_transcriptome"),
                         freq = c(7,16,
                                  5,44,
                                  28))

dat.sankey2 <- dat.sankey %>% 
  uncount(freq) %>% 
  make_long(Sample2,Sample2_depleted) %>%
  mutate(node = fct_relevel(node, "other_genomic", "targeted_genomic_interval", "percentage_total_transcriptome"), 
         next_node = fct_relevel(next_node, "percentage_total_transcriptome", "other_genomic", "targeted_genomic_interval"))

ggplot(dat.sankey2, aes(x = x,
                        next_x = next_x,
                        node = node,
                        next_node = next_node,
                        fill = factor(node))) +
  geom_sankey(flow.alpha = .6) +
  theme_alluvial(base_size = 8) +
  scale_fill_manual(values = c("other_genomic" = "#440154FF",
                               "targeted_genomic_interval" = "#39568CFF",
                               "percentage_total_transcriptome" = "#1F968BFF"))


#Library complexity
dat.lib <- data.frame(downsampling = c(1,2,3,4,5,6,
                                       7),
                      Sample2 = c(500,1120,3234,5963,8870,
                                  12215,17405),
                      Sample2_depleted = c(6435,9159,12074,13730,14991,
                                           16377,17652)) 


ggplot(dat.lib, aes(x=downsampling)) + 
  geom_line(aes(y = Sample2), color="steelblue", linetype="twodash") + 
  geom_line(aes(y = Sample2_depleted), color="steelblue", linetype="twodash") + theme_classic()



#Decrease in ribomito
dat.mito <- data.frame(Sample = c("Sample1","Sample1","Sample2","Sample2"),
                       Depl = c("Non-depleted","Depleted","Non-depleted","Depleted"),
                       Fold = c((0.6/0.26)+1,1,(3.34/1.24)+1,1))
dat.mito$Depl <- factor(dat.mito$Depl, levels = c("Non-depleted","Depleted"))

ggplot(dat.mito, mapping = aes(x=Sample, y=Fold, fill=Depl)) +
  geom_bar(stat="identity", position = "dodge") +
  theme_classic()


############################################################
#                                                          #
#                           Data                           #
#                                                          #
############################################################

syd4.1 <- readRDS("data/sample_Sydney4_saw_v2_bin50_seurat.rds")
hk.syd4.1 <- readRDS("data/sample_Sydney4_HK_saw_v2_bin50_seurat.rds") 


############################################################
#                                                          #
#                        Data QC                           #
#                                                          #
############################################################

#remove genes in less than 5 bins
keep <- rownames(syd4.1@assays$Spatial@counts[rowSums(syd4.1@assays$Spatial@counts)>5,])
syd4.1 <- subset(syd4.1, features = keep)

#remove bins with less than 5 genes
keep <- colnames(syd4.1@assays$Spatial@counts[,colSums(syd4.1@assays$Spatial@counts)>5])
syd4.1 <- subset(syd4.1, cells = keep)


syd4.1[["percent.mt"]] <- PercentageFeatureSet(syd4.1, pattern = "^MT-")
syd4.1[["percent.rpl"]] <- PercentageFeatureSet(syd4.1, pattern = "^RPL")
syd4.1[["percent.rps"]] <- PercentageFeatureSet(syd4.1, pattern = "^RPS")
syd4.1[["percent.ribo"]] <- syd4.1[["percent.rpl"]] + syd4.1[["percent.rps"]]
syd4.1[["ribo.mito"]] <- syd4.1[["percent.mt"]] + syd4.1[["percent.ribo"]]


syd4.1@images$slice1@coordinates$row <- trunc(syd4.1@images$slice1@coordinates$row/50) + 1
syd4.1@images$slice1@coordinates$col <- trunc(syd4.1@images$slice1@coordinates$col/50) + 1
syd4.1@images$slice1@coordinates$imagerow <- trunc(syd4.1@images$slice1@coordinates$imagerow/50) + 1
syd4.1@images$slice1@coordinates$imagecol <- trunc(syd4.1@images$slice1@coordinates$imagecol/50) + 1
syd4.1@images$slice1@image <- syd4.1@images$slice1@image[1:max(syd4.1@images$slice1@coordinates$imagerow),
                                                         1:max(syd4.1@images$slice1@coordinates$imagecol)]

syd4.1 <- SCTransform(syd4.1, assay = "Spatial",  verbose = TRUE, return.only.var.genes = FALSE, variable.features.n = 2000)

var <- c('gmean', 'variance', 'residual_variance')
syd4.1[["SCT"]]@meta.features <- SCTResults(syd4.1[["SCT"]], slot = "feature.attributes")[, var]


#remove genes in less than 5 bins
keep <- rownames(hk.syd4.1@assays$Spatial@counts[rowSums(hk.syd4.1@assays$Spatial@counts)>5,])
hk.syd4.1 <- subset(hk.syd4.1, features = keep)

#remove bins with less than 5 genes
keep <- colnames(hk.syd4.1@assays$Spatial@counts[,colSums(hk.syd4.1@assays$Spatial@counts)>5])
hk.syd4.1 <- subset(hk.syd4.1, cells = keep)


hk.syd4.1[["percent.mt"]] <- PercentageFeatureSet(hk.syd4.1, pattern = "^MT-")
hk.syd4.1[["percent.rpl"]] <- PercentageFeatureSet(hk.syd4.1, pattern = "^RPL")
hk.syd4.1[["percent.rps"]] <- PercentageFeatureSet(hk.syd4.1, pattern = "^RPS")
hk.syd4.1[["percent.ribo"]] <- hk.syd4.1[["percent.rpl"]] + hk.syd4.1[["percent.rps"]]
hk.syd4.1[["ribo.mito"]] <- hk.syd4.1[["percent.mt"]] + hk.syd4.1[["percent.ribo"]]


hk.syd4.1@images$slice1@coordinates$row <- trunc(hk.syd4.1@images$slice1@coordinates$row/50) + 1
hk.syd4.1@images$slice1@coordinates$col <- trunc(hk.syd4.1@images$slice1@coordinates$col/50) + 1
hk.syd4.1@images$slice1@coordinates$imagerow <- trunc(hk.syd4.1@images$slice1@coordinates$imagerow/50) + 1
hk.syd4.1@images$slice1@coordinates$imagecol <- trunc(hk.syd4.1@images$slice1@coordinates$imagecol/50) + 1
hk.syd4.1@images$slice1@image <- hk.syd4.1@images$slice1@image[1:max(hk.syd4.1@images$slice1@coordinates$imagerow),
                                                               1:max(hk.syd4.1@images$slice1@coordinates$imagecol)]

hk.syd4.1 <- SCTransform(hk.syd4.1, assay = "Spatial",  verbose = TRUE, return.only.var.genes = FALSE, variable.features.n = 2000)

var <- c('gmean', 'variance', 'residual_variance')
hk.syd4.1[["SCT"]]@meta.features <- SCTResults(hk.syd4.1[["SCT"]], slot = "feature.attributes")[, var]



############################################################
#                                                          #
#                        HVGs - Fig S1                     #
#                                                          #
############################################################

#histogram for cutoff
data <- data.frame('residual_variance' = c(syd4.1[["SCT"]]@meta.features$sct.residual_variance,
                                           hk.syd4.1[["SCT"]]@meta.features$sct.residual_variance,
                                           rep(NA, (length(syd4.1[["SCT"]]@meta.features$sct.residual_variance)-length(hk.syd4.1[["SCT"]]@meta.features$sct.residual_variance)))),
                   'state' = c(rep('non-depleted', length(syd4.1[["SCT"]]@meta.features$sct.residual_variance)),
                               rep('depleted', length(syd4.1[["SCT"]]@meta.features$sct.residual_variance))))


ggplot(data, aes(residual_variance, fill = state)) + 
  geom_histogram(bins = 100, alpha=.5, position="identity") + 
  theme_classic() + 
  scale_fill_manual(values=c("#330066", "#000033")) + 
  ylab("Count") + xlab("Residual variance") + xlim(c(0,3)) + 
  geom_vline(xintercept = 1.4, linetype = 'dotted')



#cut-off by residual variance of 1.4
nrow(syd4.1[["SCT"]]@meta.features[syd4.1[["SCT"]]@meta.features$sct.residual_variance>1.4,])
nrow(hk.syd4.1[["SCT"]]@meta.features[hk.syd4.1[["SCT"]]@meta.features$sct.residual_variance>1.4,]) 



table(rownames(syd4.1[["SCT"]]@meta.features[syd4.1[["SCT"]]@meta.features$sct.residual_variance>1.4,]) %in% rownames(hk.syd4.1[["SCT"]]@meta.features[hk.syd4.1[["SCT"]]@meta.features$sct.residual_variance>1.4,]))
nondepl <- rownames(syd4.1[["SCT"]]@meta.features[syd4.1[["SCT"]]@meta.features$sct.residual_variance>1.4,])
depl <- rownames(hk.syd4.1[["SCT"]]@meta.features[hk.syd4.1[["SCT"]]@meta.features$sct.residual_variance>1.4,])


overlap <- depl[depl %in% nondepl]
non.depleted <- nondepl[!nondepl %in% depl]
depleted <- depl[!depl %in% nondepl] 





#non-depleted HVGs
dat <- data.frame('genes' = rownames(syd4.1@assays$SCT@data[rownames(syd4.1@assays$SCT@data) %in% non.depleted,]),
                  'exp' = rowSums(syd4.1@assays$SCT@data[rownames(syd4.1@assays$SCT@data) %in% non.depleted,]))
dat.hk <- data.frame('genes' = rownames(hk.syd4.1@assays$SCT@data[rownames(hk.syd4.1@assays$SCT@data) %in% non.depleted,]),
                     'exp' = rowSums(hk.syd4.1@assays$SCT@data[rownames(hk.syd4.1@assays$SCT@data) %in% non.depleted,]))
dat.all <- merge(x = dat, y = dat.hk, by = 'genes')
colnames(dat.all) <- c('Genes','Non-depleted','Depleted')


highlight_mt <- dat.all %>% filter(grepl("^MT-",dat.all$Genes))
highlight_rpl <- dat.all %>% filter(grepl("^RPL",dat.all$Genes))
highlight_rps <- dat.all %>% filter(grepl("^RPS",dat.all$Genes))
highlight <- rbind(highlight_mt,highlight_rpl,highlight_rps)


ggplot(dat.all, aes(dat.all$`Non-depleted`,dat.all$Depleted)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, size = 0.5) +
  scale_x_sqrt() + scale_y_sqrt() +
  ggrepel::geom_text_repel(aes(label = ifelse(dat.all$`Non-depleted`>1000,as.character(dat.all$Genes),'')),
                           min.segment.length = 0) + theme_classic() +
  geom_point(data=highlight, 
             aes(highlight$`Non-depleted`,highlight$Depleted), 
             color='red')


#depleted HVGs
dat <- data.frame('genes' = rownames(syd4.1@assays$SCT@data[rownames(syd4.1@assays$SCT@data) %in% depleted,]),
                  'exp' = rowSums(syd4.1@assays$SCT@data[rownames(syd4.1@assays$SCT@data) %in% depleted,]))
dat.hk <- data.frame('genes' = rownames(hk.syd4.1@assays$SCT@data[rownames(hk.syd4.1@assays$SCT@data) %in% depleted,]),
                     'exp' = rowSums(hk.syd4.1@assays$SCT@data[rownames(hk.syd4.1@assays$SCT@data) %in% depleted,]))
dat.all <- merge(x = dat, y = dat.hk, by = 'genes')
colnames(dat.all) <- c('Genes','Non-depleted','Depleted')

highlight_mt <- dat.all %>% filter(grepl("^MT-",dat.all$Genes))
highlight_rpl <- dat.all %>% filter(grepl("^RPL",dat.all$Genes))
highlight_rps <- dat.all %>% filter(grepl("^RPS",dat.all$Genes))
highlight <- rbind(highlight_mt,highlight_rpl,highlight_rps)


ggplot(dat.all, aes(dat.all$`Non-depleted`,dat.all$Depleted)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, size = 0.5) +
  scale_x_sqrt() + scale_y_sqrt() +
  ggrepel::geom_text_repel(aes(label = ifelse(dat.all$`Non-depleted`>1000,as.character(dat.all$Genes),'')),
                           min.segment.length = 0) + theme_classic() +
  geom_point(data=highlight, 
             aes(highlight$`Non-depleted`,highlight$Depleted), 
             color='red')


#shared HVGs
dat <- data.frame('genes' = rownames(syd4.1@assays$SCT@data[rownames(syd4.1@assays$SCT@data) %in% overlap,]),
                  'exp' = rowSums(syd4.1@assays$SCT@data[rownames(syd4.1@assays$SCT@data) %in% overlap,]))
dat.hk <- data.frame('genes' = rownames(hk.syd4.1@assays$SCT@data[rownames(hk.syd4.1@assays$SCT@data) %in% overlap,]),
                     'exp' = rowSums(hk.syd4.1@assays$SCT@data[rownames(hk.syd4.1@assays$SCT@data) %in% overlap,]))
dat.all <- merge(x = dat, y = dat.hk, by = 'genes')
colnames(dat.all) <- c('Genes','Non-depleted','Depleted')


highlight_mt <- dat.all %>% filter(grepl("^MT-",dat.all$Genes))
highlight_rpl <- dat.all %>% filter(grepl("^RPL",dat.all$Genes))
highlight_rps <- dat.all %>% filter(grepl("^RPS",dat.all$Genes))
highlight <- rbind(highlight_mt,highlight_rpl,highlight_rps)


ggplot(dat.all, aes(dat.all$`Non-depleted`,dat.all$Depleted)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, size = 0.5) +
  scale_x_sqrt() + scale_y_sqrt() +
  ggrepel::geom_text_repel(aes(label = ifelse(dat.all$`Non-depleted`>1000,as.character(dat.all$Genes),'')),
                           min.segment.length = 0) + theme_classic() +
  geom_point(data=highlight, 
             aes(highlight$`Non-depleted`,highlight$Depleted), 
             color='red')



############################################################
#                                                          #
#                   Annotation - Fig 2                     #
#                                                          #
############################################################

lymphoid.genes <- c("CD3E", 'CD3D', 'CD8A', 'CD8B', 'CD4', 'IL7R', 'CXCL13', 'TRAV', 'TRBV','GNLY', 'NCAM1', 'MZB1', 'CD79A', 'CD19')
myeloid.genes <- c("LYZ", 'CD163', 'FOLR2', 'TREM2', 'KIT', 'SIGLEC','S100A8', 'S100A9', 'C1QA', 'CLEC10A')
epithelial.genes <- c("ALB",'KRT8', 'KRT18', 'KRT19', 'KRT5', 'KRT14')
fibro.genes <- c('ACTA2', 'TAGLN', 'THY1', 'MYH11')
endo.genes <- c('PECAM1', 'VWF','PLVAP', 'ACKR1')



#subset to have the same bins
bins.overlap <- gsub(".*:","",colnames(syd4.1@assays$Spatial@counts))[gsub(".*:","",colnames(syd4.1@assays$Spatial@counts)) %in% gsub(".*:","",colnames(hk.syd4.1@assays$Spatial@counts))]


syd4.1.bins <- subset(syd4.1, cells = colnames(syd4.1@assays$Spatial@counts)[gsub(".*:","",colnames(syd4.1@assays$Spatial@counts)) %in% bins.overlap])
hk.syd4.1.bins <- subset(hk.syd4.1, cells = colnames(hk.syd4.1@assays$Spatial@counts)[gsub(".*:","",colnames(hk.syd4.1@assays$Spatial@counts)) %in% bins.overlap])





#Lymphoid
syd4.1.bins[['lymphoid']] <- colSums(syd4.1.bins@assays$SCT@counts[rownames(syd4.1.bins@assays$SCT@counts) %in% lymphoid.genes,])
hk.syd4.1.bins[['lymphoid']] <- colSums(hk.syd4.1.bins@assays$SCT@counts[rownames(hk.syd4.1.bins@assays$SCT@counts) %in% lymphoid.genes,])

syd4.1.bins <- SetIdent(syd4.1.bins,value='lymphoid')
SpatialFeaturePlot(syd4.1.bins, features = 'lymphoid', pt.size.factor = 2, max.cutoff = max(syd4.1.bins[['lymphoid']]), stroke = NA) + 
  ggplot2::theme(legend.position = "right")

hk.syd4.1.bins <- SetIdent(hk.syd4.1.bins,value='lymphoid')
SpatialFeaturePlot(hk.syd4.1.bins, features = 'lymphoid', pt.size.factor = 2, max.cutoff = max(syd4.1.bins[['lymphoid']]), stroke = NA) + 
  ggplot2::theme(legend.position = "right")


x <- as.data.frame(table(syd4.1.bins[['lymphoid']]))
x$Var1 <- as.integer(x$Var1)-1
y <- as.data.frame(table(hk.syd4.1.bins[['lymphoid']]))
y$Var1 <- as.integer(y$Var1)-1
x2 <- c()
for(i in 2:(nrow(x))){
  dat <- rep(x[i,1],x[i,2])
  x2 <- c(x2,dat)
}
y2 <- c()
for(i in 2:(nrow(y))){
  dat <- rep(y[i,1],y[i,2])
  y2 <- c(y2,dat)
}
tab <- data.frame('non-depleted' = c(x2, rep(NA, length(y2)-length(x2))),
                  'depleted' = y2)
tab2 <- tidyr::gather(tab, key = "Status",value = "count","non.depleted":"depleted")
ggplot(tab2, aes(tab2$count, fill = Status)) + geom_histogram(binwidth=1, alpha=.5, position="identity") + 
  scale_y_sqrt() + theme_classic() + scale_fill_manual(values=c("#330066", "#000033")) + ylab("Number of bins") + xlab("Sum of nCount") + 
  scale_x_continuous(breaks=seq(1, max(tab2$count[which(!tab2$count %in% c(NA))]), 1))




#Myeloid
syd4.1.bins[['myeloid']] <- colSums(syd4.1.bins@assays$SCT@counts[rownames(syd4.1.bins@assays$SCT@counts) %in% myeloid.genes,])
hk.syd4.1.bins[['myeloid']] <- colSums(hk.syd4.1.bins@assays$SCT@counts[rownames(hk.syd4.1.bins@assays$SCT@counts) %in% myeloid.genes,])

syd4.1.bins <- SetIdent(syd4.1.bins,value='myeloid')
SpatialFeaturePlot(syd4.1.bins, features = 'myeloid', pt.size.factor = 2, max.cutoff = max(syd4.1.bins[['myeloid']]), stroke = NA) + 
  ggplot2::theme(legend.position = "right")

hk.syd4.1.bins <- SetIdent(hk.syd4.1.bins,value='myeloid')
SpatialFeaturePlot(hk.syd4.1.bins, features = 'myeloid', pt.size.factor = 2, max.cutoff = max(syd4.1.bins[['myeloid']]), stroke = NA) + 
  ggplot2::theme(legend.position = "right")


x <- as.data.frame(table(syd4.1.bins[['myeloid']]))
x$Var1 <- as.integer(x$Var1)-1
y <- as.data.frame(table(hk.syd4.1.bins[['myeloid']]))
y$Var1 <- as.integer(y$Var1)-1
x2 <- c()
for(i in 2:(nrow(x))){
  dat <- rep(x[i,1],x[i,2])
  x2 <- c(x2,dat)
}
y2 <- c()
for(i in 2:(nrow(y))){
  dat <- rep(y[i,1],y[i,2])
  y2 <- c(y2,dat)
}
tab <- data.frame('non-depleted' = c(x2, rep(NA, length(y2)-length(x2))),
                  'depleted' = y2)
tab2 <- tidyr::gather(tab, key = "Status",value = "count","non.depleted":"depleted")
ggplot(tab2, aes(tab2$count, fill = Status)) + 
  geom_histogram(binwidth=1, alpha=.5, position="identity") + scale_y_sqrt() + theme_classic() + 
  scale_fill_manual(values=c("#330066", "#000033")) + ylab("Number of bins") + xlab("Sum of nCount") + 
  scale_x_continuous(breaks=seq(1, max(tab2$count[which(!tab2$count %in% c(NA))]), 1))




#Epithelial
syd4.1.bins[['epithelial']] <- colSums(syd4.1.bins@assays$SCT@counts[rownames(syd4.1.bins@assays$SCT@counts) %in% epithelial.genes,])
hk.syd4.1.bins[['epithelial']] <- colSums(hk.syd4.1.bins@assays$SCT@counts[rownames(hk.syd4.1.bins@assays$SCT@counts) %in% epithelial.genes,])

x <- as.data.frame(table(syd4.1.bins[['epithelial']]))
x$Var1 <- as.integer(x$Var1)-1
y <- as.data.frame(table(hk.syd4.1.bins[['epithelial']]))
y$Var1 <- as.integer(y$Var1)-1
x2 <- c()
for(i in 2:(nrow(x))){
  dat <- rep(x[i,1],x[i,2])
  x2 <- c(x2,dat)
}
y2 <- c()
for(i in 2:(nrow(y))){
  dat <- rep(y[i,1],y[i,2])
  y2 <- c(y2,dat)
}
tab <- data.frame('non-depleted' = c(x2, rep(NA, length(y2)-length(x2))),
                  'depleted' = y2)
tab2 <- tidyr::gather(tab, key = "Status",value = "count","non.depleted":"depleted")
ggplot(tab2, aes(tab2$count, fill = Status)) + geom_histogram(binwidth=1, alpha=.5, position="identity") + scale_y_sqrt() + theme_classic() + 
  scale_fill_manual(values=c("#330066", "#000033")) + ylab("Number of bins") + xlab("Sum of nCount") + 
  scale_x_continuous(breaks=seq(1, max(tab2$count[which(!tab2$count %in% c(NA))]), 10))




#Fibroblast
syd4.1.bins[['fibroblast']] <- colSums(syd4.1.bins@assays$SCT@counts[rownames(syd4.1.bins@assays$SCT@counts) %in% fibro.genes,])
hk.syd4.1.bins[['fibroblast']] <- colSums(hk.syd4.1.bins@assays$SCT@counts[rownames(hk.syd4.1.bins@assays$SCT@counts) %in% fibro.genes,])

syd4.1.bins <- SetIdent(syd4.1.bins,value='fibroblast')
SpatialFeaturePlot(syd4.1.bins, features = 'fibroblast', pt.size.factor = 2, max.cutoff = max(syd4.1.bins[['fibroblast']]), stroke = NA) + 
  ggplot2::theme(legend.position = "right")

hk.syd4.1.bins <- SetIdent(hk.syd4.1.bins,value='fibroblast')
SpatialFeaturePlot(hk.syd4.1.bins, features = 'fibroblast', pt.size.factor = 2, max.cutoff = max(syd4.1.bins[['fibroblast']]), stroke = NA) + 
  ggplot2::theme(legend.position = "right")


x <- as.data.frame(table(syd4.1.bins[['fibroblast']]))
x$Var1 <- as.integer(x$Var1)-1
y <- as.data.frame(table(hk.syd4.1.bins[['fibroblast']]))
y$Var1 <- as.integer(y$Var1)-1
x2 <- c()
for(i in 2:(nrow(x))){
  dat <- rep(x[i,1],x[i,2])
  x2 <- c(x2,dat)
}
y2 <- c()
for(i in 2:(nrow(y))){
  dat <- rep(y[i,1],y[i,2])
  y2 <- c(y2,dat)
}
tab <- data.frame('non-depleted' = c(x2, rep(NA, length(y2)-length(x2))),
                  'depleted' = y2)
tab2 <- tidyr::gather(tab, key = "Status",value = "count","non.depleted":"depleted")
ggplot(tab2, aes(tab2$count, fill = Status)) + geom_histogram(binwidth=1, alpha=.5, position="identity") + scale_y_sqrt() + theme_classic() + 
  scale_fill_manual(values=c("#330066", "#000033")) + ylab("Number of bins") + xlab("Sum of nCount") + 
  scale_x_continuous(breaks=seq(1, max(tab2$count[which(!tab2$count %in% c(NA))]), 1))




#Endothelial
syd4.1.bins[['endothelial']] <- colSums(syd4.1.bins@assays$SCT@counts[rownames(syd4.1.bins@assays$SCT@counts) %in% endo.genes,])
hk.syd4.1.bins[['endothelial']] <- colSums(hk.syd4.1.bins@assays$SCT@counts[rownames(hk.syd4.1.bins@assays$SCT@counts) %in% endo.genes,])

syd4.1.bins <- SetIdent(syd4.1.bins,value='endothelial')
SpatialFeaturePlot(syd4.1.bins, features = 'endothelial', pt.size.factor = 2, max.cutoff = max(syd4.1.bins[['endothelial']]), stroke = NA) + 
  ggplot2::theme(legend.position = "right")

hk.syd4.1.bins <- SetIdent(hk.syd4.1.bins,value='endothelial')
SpatialFeaturePlot(hk.syd4.1.bins, features = 'endothelial', pt.size.factor = 2, max.cutoff = max(syd4.1.bins[['endothelial']]), stroke = NA) + 
  ggplot2::theme(legend.position = "right")


x <- as.data.frame(table(syd4.1.bins[['endothelial']]))
x$Var1 <- as.integer(x$Var1)-1
y <- as.data.frame(table(hk.syd4.1.bins[['endothelial']]))
y$Var1 <- as.integer(y$Var1)-1
x2 <- c()
for(i in 2:(nrow(x))){
  dat <- rep(x[i,1],x[i,2])
  x2 <- c(x2,dat)
}
y2 <- c()
for(i in 2:(nrow(y))){
  dat <- rep(y[i,1],y[i,2])
  y2 <- c(y2,dat)
}
tab <- data.frame('non-depleted' = c(x2, rep(NA, length(y2)-length(x2))),
                  'depleted' = y2)
tab2 <- tidyr::gather(tab, key = "Status",value = "count","non.depleted":"depleted")
ggplot(tab2, aes(tab2$count, fill = Status)) + geom_histogram(binwidth=1, alpha=.5, position="identity") + scale_y_sqrt() + theme_classic() + 
  scale_fill_manual(values=c("#330066", "#000033")) + ylab("Number of bins") + xlab("Sum of nCount") + 
  scale_x_continuous(breaks=seq(1, max(tab2$count[which(!tab2$count %in% c(NA))]), 1))







############################################################
#                                                          #
#             Spatial gene expression - Fig 3              #
#                                                          #
############################################################

SpatialFeaturePlot(syd4.1.bins, slot = 'counts', features = "AFP", pt.size.factor = 1.5, stroke = NA) + 
  ggplot2::theme(legend.position = "right")
SpatialFeaturePlot(hk.syd4.1.bins, slot = 'counts', features = "AFP", pt.size.factor = 1.5, stroke = NA, 
                   max.cutoff = max(syd4.1.bins@assays$SCT@counts[rownames(syd4.1.bins@assays$SCT@counts) == 'AFP',])) + 
  ggplot2::theme(legend.position = "right")


afp <- data.frame("non-depleted" = syd4.1.bins@assays$SCT@counts[rownames(syd4.1.bins@assays$SCT@counts) == 'AFP',],
                  "depleted" = hk.syd4.1.bins@assays$SCT@counts[rownames(hk.syd4.1.bins@assays$SCT@counts) == 'AFP',])
afp.2 <- tidyr::gather(afp, key = "Status",value = "count","non.depleted":"depleted")
afp.2[afp.2$count == 0,]$count <- NA

ggplot(afp.2, aes(count, fill = Status)) + 
  geom_histogram(binwidth=1, alpha=.5, position="identity") + scale_y_sqrt() + 
  theme_classic() + scale_fill_manual(values=c("#330066", "#000033")) + 
  ylab("Number of bins") + xlab("nCount")




SpatialFeaturePlot(syd4.1.bins, slot = 'counts', features = "HP", pt.size.factor = 1.5, stroke = NA, 
                   max.cutoff = max(hk.syd4.1.bins@assays$SCT@counts[rownames(hk.syd4.1.bins@assays$SCT@counts) == 'HP',])) + 
  ggplot2::theme(legend.position = "right")
SpatialFeaturePlot(hk.syd4.1.bins, slot = 'counts', features = "HP", pt.size.factor = 1.5, stroke = NA) + 
  ggplot2::theme(legend.position = "right")


hp <- data.frame("non-depleted" = syd4.1.bins@assays$SCT@counts[rownames(syd4.1.bins@assays$SCT@counts) == 'HP',],
                 "depleted" = hk.syd4.1.bins@assays$SCT@counts[rownames(hk.syd4.1.bins@assays$SCT@counts) == 'HP',])
hp.2 <- tidyr::gather(hp, key = "Status",value = "count","non.depleted":"depleted")
hp.2[hp.2$count == 0,]$count <- NA

ggplot(hp.2, aes(count, fill = Status)) + 
  geom_histogram(binwidth=1, alpha=.5, position="identity") + scale_y_sqrt() + 
  theme_classic() + scale_fill_manual(values=c("#330066", "#000033")) + 
  ylab("Number of bins") + xlab("nCount")




SpatialFeaturePlot(syd4.1.bins, slot = 'counts', features = "APOC1", pt.size.factor = 1.5, stroke = NA) + 
  ggplot2::theme(legend.position = "right")
SpatialFeaturePlot(hk.syd4.1.bins, slot = 'counts', features = "APOC1", pt.size.factor = 1.5, stroke = NA, 
                   max.cutoff = max(syd4.1.bins@assays$SCT@counts[rownames(syd4.1.bins@assays$SCT@counts) == 'APOC1',])) + 
  ggplot2::theme(legend.position = "right")


apoc1 <- data.frame("non-depleted" = syd4.1.bins@assays$SCT@counts[rownames(syd4.1.bins@assays$SCT@counts) == 'APOC1',],
                    "depleted" = hk.syd4.1.bins@assays$SCT@counts[rownames(hk.syd4.1.bins@assays$SCT@counts) == 'APOC1',])
apoc1.2 <- tidyr::gather(apoc1, key = "Status",value = "count","non.depleted":"depleted")
apoc1.2[apoc1.2$count == 0,]$count <- NA

ggplot(apoc1.2, aes(count, fill = Status)) + 
  geom_histogram(binwidth=1, alpha=.5, position="identity") + scale_y_sqrt() + 
  theme_classic() + scale_fill_manual(values=c("#330066", "#000033")) + 
  ylab("Number of bins") + xlab("nCount")



