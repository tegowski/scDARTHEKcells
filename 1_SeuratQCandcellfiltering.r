library(Seurat)
library(tidyverse)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)

for (file in c("WT1", "WT2", "WT3", "WT4", "WT5", "mut1", "mut2", "mut3", "mut4", "mut5")){
        seurat_data <- Read10X(data.dir = paste0("$WORKDIR/Data/", file))
        seurat_obj <- CreateSeuratObject(counts = seurat_data,  
											project = file)
        assign(file, seurat_obj)
}

SeuratSMARTcombined <- merge(x = WT1, y = c(WT2, WT3, WT4, WT5, mut1, mut2, mut3))
SeuratSMARTcombined$log10GenesPerRead <- log10(SeuratSMARTcombined$nFeature_RNA) / log10(SeuratSMARTcombined$nCount_RNA)
SeuratSMARTcombined$mitoRatio <- PercentageFeatureSet(object = SeuratSMARTcombined, pattern = "^MT-")
SeuratSMARTcombined$mitoRatio <- SeuratSMARTcombined@meta.data$mitoRatio / 100
metadata <- SeuratSMARTcombined@meta.data
metadata$cells <- rownames(metadata)
metadata <- metadata %>%
        dplyr::rename(seq_folder = orig.ident,
                      nRead = nCount_RNA,
                      nGene = nFeature_RNA)
metadata$sample <- NA
metadata$sample[which(str_detect(metadata$cells, "^WT"))] <- "WT"					  
metadata$sample[which(str_detect(metadata$cells, "^mut"))] <- "MUT"
SeuratSMARTcombined@meta.data <- metadata
save(SeuratSMARTcombined, file="Data/293SMARTDART_seurat.RData")

metadata %>% 
  	ggplot(aes(x=sample, fill=sample)) + 
  	geom_bar() +
  	theme_classic() +
  	theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  	theme(plot.title = element_text(hjust=0.5, face="bold")) +
  	ggtitle("NCells")

metadata %>% 
  	ggplot(aes(color=sample, x=nUMI, fill= sample)) + 
  	geom_density(alpha = 0.2) + 
  	scale_x_log10() + 
  	theme_classic() +
  	ylab("Cell density") +
  	geom_vline(xintercept = 500)
	
metadata %>% 
  	ggplot(aes(color=sample, x=nGene, fill= sample)) + 
  	geom_density(alpha = 0.2) + 
  	theme_classic() +
  	scale_x_log10() + 
  	geom_vline(xintercept = 300)

metadata %>% 
  	ggplot(aes(color=sample, x=mitoRatio, fill=sample)) + 
  	geom_density(alpha = 0.2) + 
  	scale_x_log10() + 
  	theme_classic() +
  	geom_vline(xintercept = 0.2)

metadata %>%
  	ggplot(aes(x=log10GenesPerUMI, color = sample, fill=sample)) +
  	geom_density(alpha = 0.2) +
  	theme_classic() +
  	geom_vline(xintercept = 0.8)
	
filtered_Seurat <- subset(x = SeuratSMARTcombined, 
                         subset= (nUMI >= 1000000) & 
                           (nGene >= 9000) & 
                           (log10GenesPerUMI > 0.58) & 
                           (mitoRatio < 0.10))

counts <- GetAssayData(object = filtered_Seurat, slot = "counts")
nonzero <- counts > 0
keep_genes <- Matrix::rowSums(nonzero) >= 10
filtered_counts <- counts[keep_genes, ]
filtered_Seurat <- CreateSeuratObject(filtered_counts, meta.data = filtered_Seurat@meta.data)
metadata_clean <- filtered_Seurat@meta.data

metadata_clean %>% 
  	ggplot(aes(x=sample, fill=sample)) + 
  	geom_bar() +
  	theme_classic() +
  	theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  	theme(plot.title = element_text(hjust=0.5, face="bold")) +
  	ggtitle("NCells")

metadata_clean %>% 
  	ggplot(aes(color=sample, x=nUMI, fill= sample)) + 
  	geom_density(alpha = 0.2) + 
  	scale_x_log10() + 
  	theme_classic() +
  	ylab("Cell density") +
  	geom_vline(xintercept = 500)

metadata_clean %>% 
  	ggplot(aes(color=sample, x=nGene, fill= sample)) + 
  	geom_density(alpha = 0.2) + 
  	theme_classic() +
  	scale_x_log10() + 
  	geom_vline(xintercept = 300)

metadata_clean %>% 
  	ggplot(aes(x=sample, y=log10(nGene), fill=sample)) + 
  	geom_boxplot() + 
  	theme_classic() +
  	theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  	theme(plot.title = element_text(hjust=0.5, face="bold")) +
  	ggtitle("NCells vs NGenes")

metadata_clean %>% 
  	ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
  	geom_point() + 
	scale_colour_gradient(low = "gray90", high = "black") +
  	stat_smooth(method=lm) +
  	scale_x_log10() + 
  	scale_y_log10() + 
  	theme_classic() +
  	geom_vline(xintercept = 500) +
  	geom_hline(yintercept = 250) +
  	facet_wrap(~sample)

metadata_clean %>% 
  	ggplot(aes(x=mitoRatio, color=sample, fill=sample)) + 
  	geom_density(alpha = 0.2) + 
  	scale_x_log10() + 
  	theme_classic() +
  	geom_vline(xintercept = 0.2)

metadata_clean %>%
  	ggplot(aes(x=log10GenesPerUMI, color = sample, fill=sample)) +
  	geom_density(alpha = 0.2) +
  	theme_classic() +
  	geom_vline(xintercept = 0.8)

save(filtered_Seurat, file="Data/293SMARTDART_filtered.RData")

length(filtered_Seurat@meta.data$nGene)