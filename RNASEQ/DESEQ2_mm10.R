library('DESeq2')
library('data.table')

directory<-getwd()
files = grep("*mm10.genes.results", list.files(directory), value = T)

countData = data.frame(fread(files[1]))[c(1,5)]

# Loop and read the 5th column remaining files
for(i in 2:length(files)) {
        countData = cbind(countData, data.frame(fread(files[i]))[5])
}

colnames(countData) = c("GeneID", gsub(".genes.results", "", files))
rownames(countData) = countData$GeneID

countData = countData[,c(2:ncol(countData))]
countData = round(countData)  

sampleData<-read.table("sample_info.txt", header=T)

#Remove G841-05, G841-12 (outliers)
countData = countData[,-c(5, 12)]
sampleData = sampleData[-c(5, 12),]
rownames(sampleData)<-sampleData$ID

dds = DESeqDataSetFromMatrix(countData = countData, 
colData = sampleData, design = ~ GROUP)

dds = DESeq(dds)
res<-results(dds)
write.table(res,file="DEG_DRD_vs_CTRL_DESEQ2_mm10.txt",quote=FALSE,sep="\t")
resOrdered<-res[order(res$padj),]
res_sig = subset(resOrdered, padj<0.1)
write.table(res_sig,file="DEG_DRD_vs_CTRL_DESEQ2_SIG_mm10.txt",quote=FALSE,sep="\t")
