install.packages("BiocManager")

BiocManager::install("edgeR")

BiocManager::install("DESeq2")

install.packages("pheatmap")

install.packages("ggplot2")

?pheatmap # help page for pheatmap

library(edgeR)
library(DESeq2)
library(pheatmap)
library(ggplot2)


#setting up data
setwd('C:/Users/Asus/Documents/FOR-271 RNAseq Project')
Coinfection.targets<-read.delim("C:/Users/Asus/Documents/FOR-271 RNAseq Project/data/fileDesc.txt")
rownames(Coinfection.targets)<-c("Ha1","Ha2","Ha3","Ctr1","Ctr2","Ctr3")
Coinfection.orig <- readDGE(Coinfection.targets, header=F)
head(Coinfection.orig)
Coinfection.rawCount <- Coinfection.orig$count

#examining data
dim(Coinfection.rawCount)
head(Coinfection.rawCount)

#setting up metadata
sampletype <- factor(c(rep("Ha",3), rep("Ctr", 3)))
meta <- data.frame(sampletype, row.names = colnames(Coinfection.orig$count))
colnames(Coinfection.orig$count)
rownames(meta)
all(colnames(Coinfection.orig$count) %in% rownames(meta)) & all(rownames(meta) %in% colnames(Coinfection.orig$count)) #final check
  
dds <- DESeqDataSetFromMatrix(Coinfection.orig, colData = meta, design = ~ sampletype)
head(counts(dds))
dds <- estimateSizeFactors(dds)
sizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)
dir.create("results") #added this in so less "manual handling" needed
write.csv(normalized_counts, file="./results/coinfection_normalized_counts_DESeq2.csv")

rld <- rlog(dds, blind=TRUE)
plotPCA(rld, intgroup="sampletype")
pdf("./results/PlotPCA_dds.pdf")
plotPCA(rld, intgroup="sampletype")
dev.off()


rld_mat <- assay(rld)
rld_cor <- cor(rld_mat) 
head(rld_cor)
head(meta)

pheatmap(rld_cor, annotation = meta)

pdf("./results/PlotHeatmap_dds.pdf")
heat.colors <- RColorBrewer::brewer.pal(6, "Blues")
pheatmap(rld_cor, annotation = meta, color = heat.colors, border_color=NA, fontsize = 10, 
         fontsize_row = 10, height=20)

#getwd()
#setwd("C:/Users/username/Desktop/RNA-seq_DEA")
options(digits=3)
infection.targets<-read.delim("./data/fileDesc.txt")
infection.targets
rownames(infection.targets)<-c("Ha1","Ha2","Ha3","Ctr1","Ctr2","Ctr3")
infection.targets
infection <- readDGE(infection.targets, header=F)

dim(infection)
head(infection)
infection.rawCount <- infection$count
head(infection.rawCount)

ggplot(infection.rawCount) +
  geom_histogram(aes(x = Ha1), stat = "bin", bins = 200) +
  xlab("Raw expression counts") +
  ylab("Number of genes")

png("./results/count distribution.png", res=300, height=1800, width=1800)
ggplot(infection.rawCount) +
  geom_histogram(aes(x = Ha1), stat = "bin", bins = 200) +
  xlab("Raw expression counts") +
  ylab("Number of genes")
dev.off()
write.csv(infection.rawCount, file="./results/infection.rawCounts.csv")
infection.normCPM <- cpm(calcNormFactors(infection)) #counts per mil
dim(infection.normCPM)

head(infection.normCPM)
write.csv(infection.normCPM, file="./results/infection.normCPM.csv")
infection.filtered <- rowSums(cpm(infection)>1) >=3
table(infection.filtered)
infection$samples$lib.size
Infection <- infection[infection.filtered,]
dim(Infection)

