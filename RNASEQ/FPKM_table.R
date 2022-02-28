library('data.table')

directory<-getwd()
files = grep("*genes.results", list.files(directory), value = T)

countData = data.frame(fread(files[1]))[c(1,7)]

# Loop and read the 7th column (FPKM) remaining files
for(i in 2:length(files)) {
        countData = cbind(countData, data.frame(fread(files[i]))[7])
}

colnames(countData) = c("GeneID", gsub(".genes.results", "", files))
rownames(countData) = countData$GeneID

countData = countData[,c(2:ncol(countData))]

write.table(countData,file="FPKM_GeneID.txt",quote=FALSE,sep="\t")
