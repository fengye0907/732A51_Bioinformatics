library(GEOquery)
x = getGEOSuppFiles("GSE20986")
# 解压
untar("GSE20986/GSE20986_RAW.tar", exdir = "data")

cels = list.files("data/", pattern = "[gz]")
sapply(paste("data", cels, sep = "/"), gunzip)

phenodata = matrix(rep(list.files("data"), 2), ncol =2)
class(phenodata)
phenodata <- as.data.frame(phenodata)
colnames(phenodata) <- c("Name", "FileName")
phenodata$Targets <- c("iris", 
                       "retina", 
                       "retina", 
                       "iris", 
                       "retina", 
                       "iris", 
                       "choroid", 
                       "choroid", 
                       "choroid", 
                       "huvec", 
                       "huvec", 
                       "huvec")
write.table(phenodata, "data/phenodata.txt", quote = F, sep = "\t", row.names = F)

# source("https://bioconductor.org/biocLite.R")
# biocLite("simpleaffy")
library(simpleaffy)

celfiles <- read.affy(covdesc = "phenodata.txt", path = "data")
boxplot(celfiles)

library(RColorBrewer)
cols = brewer.pal(8, "Set1")
eset <- exprs(celfiles)
samples <- celfiles$Targets
colnames(eset)
# colnames(eset) <- samples

boxplot(celfiles, col = cols, las = 2)
distance <- dist(t(eset), method = "maximum")
clusters <- hclust(distance)
plot(clusters)

require(simpleaffy)
# devtools::install.github("bmbolstad/affyPLM")
require(affyPLM)
celfiles.gcrma = gcrma(celfiles)
par(mfrow=c(1,2))
boxplot(celfiles.gcrma, col = cols, las = 2, main = "Post-Normalization");
boxplot(celfiles, col = cols, las = 2, main = "Pre-Normalization")

dev.off()
# Cluster Dendrogram based on post-normalization
eset2 <- exprs(celfiles.gcrma)
# colnames(eset2) <- samples
distance2 <- dist(t(eset2), method = "maximum")
clusters2 <- hclust(distance2)
plot(clusters2)



#1.2##############################################################
library("gplots")
library("RColorBrewer")
eset <- exprs(celfiles)
samples <- celfiles$Targets
colnames(eset)
colnames(eset) <- samples
distance <- dist(t(eset), method = "maximum")
clusters <- hclust(distance)
plot(clusters)

eset2 <- exprs(celfiles.gcrma)
colnames(eset2) <- samples
distance2 <- dist(t(eset2), method = "maximum")
clusters2 <- hclust(distance2)
plot(clusters2)

colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
hc <- hclust(sampleDists)
heatmap.2( as.matrix(distance), Rowv=as.dendrogram(clusters),
           symm=TRUE, trace="none", col=colors,
           margins=c(2,10), labCol=FALSE )

heatmap.2( as.matrix(distance2), Rowv=as.dendrogram(clusters2),
           symm=TRUE, trace="none", col=colors,
           margins=c(2,10), labCol=FALSE )

#1.3##############################################################
library(limma)
phenodata

samples <- as.factor(samples)
design <- model.matrix(~0+samples)
colnames(design)
samples <- as.factor(samples)
design <- model.matrix(~0+samples)
colnames(design)
colnames(design) <- c("choroid", "huvec", "iris", "retina")
design

contrast.matrix = makeContrasts(
  huvec_choroid = huvec - choroid, 
  huvec_retina = huvec - retina, 
  huvec_iris = huvec - iris, 
  levels = design)

fit = lmFit(celfiles.gcrma, design)
huvec_fit <- contrasts.fit(fit, contrast.matrix)
huvec_ebay <- eBayes(huvec_fit)
# source("https://bioconductor.org/biocLite.R")
# biocLite("hgu133plus2.db")
library(hgu133plus2.db)
# biocLite("annotate")
library(annotate)
probenames.list <- rownames(topTable(huvec_ebay, number = 100000))
getsymbols <- getSYMBOL(probenames.list, "hgu133plus2")


results <- topTable(huvec_ebay, number = 100000, coef = "huvec_choroid")
results <- cbind(results, getsymbols)
summary(results)
results$threshold <- "1"
a <- subset(results, adj.P.Val < 0.05 & logFC > 5)
results[rownames(a), "threshold"] <- "2"
b <- subset(results, adj.P.Val < 0.05 & logFC < -5)
results[rownames(b), "threshold"] <- "3"
table(results$threshold)

library(ggplot2)
volcano <- ggplot(data = results, 
                  aes(x = logFC, y = -1*log10(adj.P.Val), 
                      colour = threshold, 
                      label = getsymbols))

volcano <- volcano + 
  geom_point() + 
  scale_color_manual(values = c("black", "red", "green"), 
                     labels = c("Not Significant", "Upregulated", "Downregulated"), 
                     name = "Key/Legend")

volcano + 
  geom_text(data = subset(results, logFC > 5 & -1*log10(adj.P.Val) > 5), 
            aes(x = logFC, y = -1*log10(adj.P.Val), colour = threshold, label = getsymbols)  )+
  ggtitle("Volcano plot for huvec-choroid pair")+theme_bw()


# huvec_retina pair
results <- topTable(huvec_ebay, number = 100000, coef = "huvec_retina")
results <- cbind(results, getsymbols)
results$threshold <- "1"
a <- subset(results, adj.P.Val < 0.05 & logFC > 5)
results[rownames(a), "threshold"] <- "2"
b <- subset(results, adj.P.Val < 0.05 & logFC < -5)
results[rownames(b), "threshold"] <- "3"
volcano <- ggplot(data = results, 
                  aes(x = logFC, y = -1*log10(adj.P.Val), 
                      colour = threshold, 
                      label = getsymbols))

volcano <- volcano + 
  geom_point() + 
  scale_color_manual(values = c("black", "red", "green"), 
                     labels = c("Not Significant", "Upregulated", "Downregulated"), 
                     name = "Key/Legend")

volcano + 
  geom_text(data = subset(results, logFC > 5 & -1*log10(adj.P.Val) > 5),
            aes(x = logFC, y = -1*log10(adj.P.Val), colour = threshold, label = getsymbols))+
  ggtitle("Volcano plot for huvec-retina pair")+theme_bw()


# huvec_iris pair
results <- topTable(huvec_ebay, number = 100000, coef = "huvec_iris")
results <- cbind(results, getsymbols)
results$threshold <- "1"
a <- subset(results, adj.P.Val < 0.05 & logFC > 5)
results[rownames(a), "threshold"] <- "2"
b <- subset(results, adj.P.Val < 0.05 & logFC < -5)
results[rownames(b), "threshold"] <- "3"
volcano <- ggplot(data = results, 
                  aes(x = logFC, y = -1*log10(adj.P.Val), 
                      colour = threshold, 
                      label = getsymbols))

volcano <- volcano + 
  geom_point() + 
  scale_color_manual(values = c("black", "red", "green"), 
                     labels = c("Not Significant", "Upregulated", "Downregulated"), 
                     name = "Key/Legend")

volcano + 
  geom_text(data = subset(results, logFC > 5 & -1*log10(adj.P.Val) > 5), 
            aes(x = logFC, y = -1*log10(adj.P.Val), 
                colour = threshold, label = getsymbols))+
  ggtitle("Volcano plot for huvec-iris pair")+theme_bw()


# huvec-retina MA 
esetd <- eset[, c("huvec","retina","iris","choroid")]
plot(esetd2[,c(1,4)])
plot(esetd2[,c(1,2)])
esetd2 <- eset2[, c("huvec","retina","iris","choroid")]


y3 <- (eset[, c("retina", "huvec")])
x3 <- (eset2[, c("retina", "huvec")])

library(affy)
split.screen(c(1,2))
ma.plot( rowMeans(log2(y3)), log2(y3[, 1])-log2(y3[, 2]), cex=1 )
title("Pre Norm (huvec Vs retina)")
screen(2) 
ma.plot( rowMeans(log2(x3)), log2(x3[, 1])-log2(x3[, 2]), cex=1 )
title("Post Norm (huvec Vs retina)")
