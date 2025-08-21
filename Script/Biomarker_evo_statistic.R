BioKA<-read.csv("Biomarker\\BioKA\\cancer_or_benign_tumor.csv")
BioKA_gentic<-BioKA[BioKA$biomarkerType %in%c("gene","protein"),]
BioKA_gentic_ID<-data.frame(Name=unique(BioKA_gentic$biomarkerName))
BioKA_gentic_ID$ID<-""
BioKA_gentic_ID$wuyu <- sub("\\[.*", "", BioKA_gentic_ID$Name )
BioKA_gentic_ID$wuyu <- sub("\\(Canis lupus dingo\\)", "", BioKA_gentic_ID$wuyu )
BioKA_gentic_ID$wuyu <- sub("\\(Canis lupus familiaris\\]", "", BioKA_gentic_ID$wuyu )
BioKA_gentic_ID$wuyu <- sub("\\s+$", "", BioKA_gentic_ID$wuyu)#删除末尾的空格
BioKA_gentic_ID$wuyu <- sub(".*\\(([^]]+)\\).*", "\\1", BioKA_gentic_ID$wuyu)#用括号里的代替，例如：stromal cell-derived factor-1 (SDF-1)	



library(org.Hs.eg.db)#org.Hs.eg.db 是用于geneID转换的包。物种为人类。Bioconductor上还有提供其他物种的。
eg2Symbol=toTable(org.Hs.egSYMBOL)##将包中gene_symbol转换成数据框
eg2name=toTable(org.Hs.egGENENAME)##将包中GENENAME转换成数据框
anno=merge(eg2Symbol,eg2name,by='gene_id')#根据gene_id合并两个数据框


# 预处理函数，用于移除空格和"-"，并转换为小写
preprocess <- function(x) {
  x <- gsub(" ", "", x) # 移除空格
  x <- gsub("-", "", x) # 移除"-"
  x <- gsub("[[:punct:]]", "", x)#移除标点
  tolower(x) # 转换为小写
}
# 对 anno 数据框的 symbol 和 gene_name 进行预处理
anno$symbol_processed <- sapply(anno$symbol, preprocess)
anno$gene_name_processed <- sapply(anno$gene_name, preprocess)

# 对 BioKA_gentic_ID 数据框的 wuyu 进行预处理
BioKA_gentic_ID$wuyu_processed <- sapply(BioKA_gentic_ID$wuyu, preprocess)
library(dplyr)
# 对每个处理后的 wuyu 值进行处理
BioKA_gentic_ID$ID <- sapply(BioKA_gentic_ID$wuyu_processed, function(wuyu_value) {
  # 首先尝试在处理后的 symbol 中查找
  symbol_match <- anno$symbol[anno$symbol_processed == wuyu_value]
  
  if (length(symbol_match) > 0) {
    # 如果在 symbol 中找到了匹配项
    return(symbol_match)
  } else {
    # 如果在 symbol 中没有找到，尝试在处理后的 gene_name 中查找
    gene_name_match <- anno$symbol[anno$gene_name_processed == wuyu_value]
    if (length(gene_name_match) > 0) {
      # 如果在 gene_name 中找到了匹配项
      return(gene_name_match)
    } else {
      # 如果都没有找到
      return("")
    }
  }
})
write.csv(BioKA_gentic_ID,"Biomarker\\BioKA\\cancer_biomarkerID.csv",fileEncoding='UTF-8')
write.table(BioKA_gentic_ID,"Biomarker\\BioKA\\cancer_biomarkerID.txt",fileEncoding='UTF-8',sep = "\t",quote = F,row.names = F)

BioKA_gentic_ID <- subset(BioKA_gentic_ID, !duplicated(BioKA_gentic_ID$wuyu))
a<-read.table("Biomarker\\BioKA\\protein_gene.txt",sep=" ",check.names = F)
a<-data.frame(a[,-1])
write.table(a,"Biomarker\\BioKA\\protein_gene.tsv",quote = F,row.names = F)

BioKA_cancer<-read.csv("Biomarker\\BioKA\\cancer_biomarkerID.csv",header = F)
BioKA_cancer<-data.frame(biomarker=unique(BioKA_cancer[,1]))
library(stringr)
BioKA_cancer$biomarker <- str_to_upper(BioKA_cancer$biomarker)
BioKA_cancer$biomarker <- gsub("-", "", BioKA_cancer$biomarker)
BioKA_cancer<-data.frame(biomarker=unique(BioKA_cancer[,1]))
evo<-read.table("进化信息\\main_HUMAN-gene-name.txt")
evo<-evo[,-1]
evo_filter<-evo[evo$V3%in%BioKA_cancer$biomarker,]


#------------------------------------------------------MarkDB--------------------------------------------------------------------------
MarkDB<-read.table("Biomarker\\MarkerDB\\biomarker_nore.tsv")
evo<-read.table("Evolutionary\\main_HUMAN-gene-name.txt")
evo<-evo[,-1]
evo_filter<-evo[evo$V3%in%MarkDB$V1,]

# 使用 table 函数计算每个 evo_stage 的频数
evo_stage_table <- as.data.frame(table(evo_filter$V2))

# 计算每个 evo_stage 的比例
evo_stage_proportions <- as.data.frame(prop.table(table(evo_filter$V2)))



# 打印结果
print(evo_stage_proportions)

#cancer
MarkDB<-read.csv("Biomarker\\MarkerDB\\selected_columns.csv",header = F)
disease<-data.frame(unique(MarkDB$V3))
#write.csv(disease,"Biomarker\\MarkerDB\\disese.csv")
cancer<-read.csv("Biomarker\\MarkerDB\\cancer.csv",header = F)
cancer_biomarker<-MarkDB[MarkDB$V3%in%cancer$V1,]
cancer_biomarker<-data.frame(unique(cancer_biomarker$V1))
evo_filter<-evo[evo$V3%in%cancer_biomarker$unique.cancer_biomarker.V1.,]
# 使用 table 函数计算每个 evo_stage 的频数
evo_stage_table <- as.data.frame(table(evo_filter$V2))

# 计算每个 evo_stage 的比例
evo_stage_proportions <- as.data.frame(prop.table(table(evo_filter$V2)))

eukaryota<- matrix(c(25,37,5217,12636), nrow = 2)
fisher.test(eukaryota)

Eumetazoa<- matrix(c(18,44,4550,13303), nrow = 2)
fisher.test(Eumetazoa)
Opisthokonta<- matrix(c(6,56,1024,16829), nrow = 2)
fisher.test(Opisthokonta)
