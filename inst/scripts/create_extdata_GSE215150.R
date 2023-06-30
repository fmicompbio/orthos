# Download and unzip gene expression tables
# GSE215150_4T1.gene.expression.txt.gz
# GSE215150_LM2.gene.expression.txt.gz
# from Supplementary Files @ https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE215150

LM2 <- read.delim("GSE215150_LM2.gene.expression.txt",comment.char = "#")
rownames(LM2) <- LM2[,1]
LM2 <- LM2[,-c(1:6)]
colnames(LM2) <- gsub(".*STAR.","",colnames(LM2) )
colnames(LM2) <- gsub("Aligned.sortedByCoord.out.bam","",colnames(LM2) )
LM2 <- cbind(rowSums(LM2[,1:3]), rowSums(LM2[,4:6]), rowSums(LM2[,7:9])  )
colnames(LM2) <- c("Ctrl","MKL1","caMKL1")
saveRDS(LM2, "inst/extdata/GSE215150_MKL1_Human.rds")

T1 <- read.delim("GSE215150_4T1.gene.expression.txt",comment.char = "#")
rownames(T1) <- T1[,1]
T1 <- T1[,-c(1:6)]
colnames(T1) <- gsub(".*STAR.","",colnames(T1) )
colnames(T1) <- gsub("Aligned.sortedByCoord.out.bam","",colnames(T1) )
T1 <- cbind(rowSums(T1[,1:3]), rowSums(T1[,4:6]), rowSums(T1[,7:9])  )
colnames(T1) <- c("Ctrl","MKL1","caMKL1")
saveRDS(T1, "inst/extdata/GSE215150_MKL1_Mouse.rds")
