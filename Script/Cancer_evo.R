#------------------------------------------------------Biomarker--------------------------------------------------------------------------
library(dplyr)
library(stringr)
evo<-read.table("Evolutionary\\main_HUMAN-gene-name.txt",header = T)
library(readxl)

background_counts <- evo$evo_stage %>%
  data.frame(Period = .) %>%
  group_by(Period) %>%
  summarise(Total = n())
#BioKA
BioKA_cancer<-read.csv("Biomarker\\BioKA\\BioKa_cancer_biomarkerID.csv",header = F)
BioKA_cancer<-data.frame(biomarker=unique(BioKA_cancer[,1]))
BioKA_biomarker<-evo[evo$gene_name%in%BioKA_cancer$biomarker,]
biomarker_counts <- BioKA_biomarker$evo_stage %>%
  data.frame(Period = .) %>%
  group_by(Period) %>%
  summarise(Count = n())
biomarker_counts$prob<-biomarker_counts$Count / nrow(BioKA_biomarker)
enrichment_data <- left_join(biomarker_counts, background_counts, by = "Period")
results <- enrichment_data %>%
  rowwise() %>%
  mutate(
    Not_Biomarker = Total - Count,
    N_N = sum(biomarker_counts$Count)-Count,
    Not_Background = sum(background_counts$Total) - Not_Biomarker - Count -N_N,
    Fisher_p_value = fisher.test(matrix(c(Count, N_N, Not_Biomarker, Not_Background), nrow = 2))$p.value,
    Fisher_odd=fisher.test(matrix(c(Count, N_N, Not_Biomarker, Not_Background), nrow = 2))$estimate,
  )
results=data.frame(results$Period,results$Fisher_p_value,results$Fisher_odd)

#MarkDB
MarkDB<-read.csv("Biomarker\\MarkerDB\\selected_columns.csv",header = F)
MarkDB_cancer<-read.csv("Biomarker\\MarkerDB\\cancer.csv",header = F)
MarkDB_cancer_biomarker<-MarkDB[MarkDB$V3%in%MarkDB_cancer$V1,]
MarkDB_cancer_biomarker<-evo[evo$gene_name%in%MarkDB_cancer_biomarker$V1,]

biomarker_counts <- MarkDB_cancer_biomarker$evo_stage %>%
  data.frame(Period = .) %>%
  group_by(Period) %>%
  summarise(Count = n())
biomarker_counts$prob<-biomarker_counts$Count / nrow(MarkDB_cancer_biomarker)
enrichment_data <- left_join(biomarker_counts, background_counts, by = "Period")
results <- enrichment_data %>%
  rowwise() %>%
  mutate(
    Not_Biomarker = Total - Count,
    N_N = sum(biomarker_counts$Count)-Count,
    Not_Background = sum(background_counts$Total) - Not_Biomarker - Count -N_N,
    Fisher_p_value = fisher.test(matrix(c(Count, N_N, Not_Biomarker, Not_Background), nrow = 2))$p.value,
    Fisher_odd=fisher.test(matrix(c(Count, N_N, Not_Biomarker, Not_Background), nrow = 2))$estimate,
  )
results=data.frame(results$Period,results$Fisher_p_value,results$Fisher_odd)

#TTD
library(data.table)
TTD_biomarker<-fread("Biomarker\\TTD/Biomarker_disease.txt",sep="\t")
#根据ICD-11代码筛选 2开的是肿瘤
TTD_cancer_biomarker<- TTD_biomarker[grepl("ICD-11: 2", TTD_biomarker$ICD11), ]
TTD_biomarker_uni<-fread("Biomarker\\TTD/TTD-P1-07-Biomarker_disease-web-uniprot.txt",header=F)
TTD_biomarker_uni<-TTD_biomarker_uni[TTD_biomarker_uni$V1%in%TTD_cancer_biomarker$BiomarkerID,]
TTD_biomarker_gene<-fread("Biomarker\\TTD/Biomarker_gene.tsv",header=T)
TTD_biomarker_gene<-TTD_biomarker_gene[TTD_biomarker_gene$From%in%TTD_biomarker_uni$V2,]
TTD_cancer_biomarker<-evo[evo$gene_name%in%TTD_biomarker_gene$To,]

biomarker_counts <- TTD_cancer_biomarker$evo_stage %>%
  data.frame(Period = .) %>%
  group_by(Period) %>%
  summarise(Count = n())
biomarker_counts$prob<-biomarker_counts$Count / nrow(TTD_cancer_biomarker)
enrichment_data <- left_join(biomarker_counts, background_counts, by = "Period")
results <- enrichment_data %>%
  rowwise() %>%
  mutate(
    Not_Biomarker = Total - Count,
    N_N = sum(biomarker_counts$Count)-Count,
    Not_Background = sum(background_counts$Total) - Not_Biomarker - Count -N_N,
    Fisher_p_value = fisher.test(matrix(c(Count, N_N, Not_Biomarker, Not_Background), nrow = 2))$p.value,
    Fisher_odd=fisher.test(matrix(c(Count, N_N, Not_Biomarker, Not_Background), nrow = 2))$estimate,
  )
results=data.frame(results$Period,results$Fisher_p_value,results$Fisher_odd)



#FDA
FDA_cancer<-read.table("Biomarker/FDA/FDA-approved-biomarkers-tumour.txt",sep="\t",row.names = 1)
FDA_biomarker_cancer<-evo[evo$gene_name%in%FDA_cancer$V2,]
biomarker_counts <- FDA_biomarker_cancer$evo_stage %>%
  data.frame(Period = .) %>%
  group_by(Period) %>%
  summarise(Count = n())
biomarker_counts$prob<-biomarker_counts$Count / nrow(FDA_biomarker_cancer)
enrichment_data <- left_join(biomarker_counts, background_counts, by = "Period")
results <- enrichment_data %>%
  rowwise() %>%
  mutate(
    Not_Biomarker = Total - Count,
    N_N = sum(biomarker_counts$Count)-Count,
    Not_Background = sum(background_counts$Total) - Not_Biomarker - Count -N_N,
    Fisher_p_value = fisher.test(matrix(c(Count, N_N, Not_Biomarker, Not_Background), nrow = 2))$p.value,
    Fisher_odd=fisher.test(matrix(c(Count, N_N, Not_Biomarker, Not_Background), nrow = 2))$estimate,
  )
results=data.frame(results$Period,results$Fisher_p_value,results$Fisher_odd)



#union
Union_cancer_biomarker <- rbind(BioKA_biomarker,MarkDB_cancer_biomarker,TTD_cancer_biomarker,FDA_biomarker_cancer)
Union_cancer_biomarker<-distinct(Union_cancer_biomarker)

biomarker_counts <- Union_cancer_biomarker$evo_stage %>%
  data.frame(Period = .) %>%
  group_by(Period) %>%
  summarise(Count = n())
biomarker_counts$prob<-biomarker_counts$Count / nrow(Union_cancer_biomarker)
enrichment_data <- left_join(biomarker_counts, background_counts, by = "Period")
results <- enrichment_data %>%
  rowwise() %>%
  mutate(
    Not_Biomarker = Total - Count,
    N_N = sum(biomarker_counts$Count)-Count,
    Not_Background = sum(background_counts$Total) - Not_Biomarker - Count -N_N,
    Fisher_p_value = fisher.test(matrix(c(Count, N_N, Not_Biomarker, Not_Background), nrow = 2))$p.value,
    Fisher_odd=fisher.test(matrix(c(Count, N_N, Not_Biomarker, Not_Background), nrow = 2))$estimate,
  )
results=data.frame(results$Period,results$Fisher_p_value,results$Fisher_odd)


#------------------------------------------------------target--------------------------------------------------------------------------
evo<-read.table("Evolutionary\\main_HUMAN-gene-name.txt",header = T)
library(readxl)
evo_stage_proportions <- data.frame(matrix(ncol = 0, nrow = 8))
rownames(evo_stage_proportions)<-as.data.frame(table(evo$evo_stage))[,1]
evo_stage_proportions$all_human <- as.data.frame(table(evo$evo_stage))[,2]
evo_stage_proportions$all_human_prob<-as.data.frame(prop.table(table(evo$evo_stage)))[,2]
background_counts <- evo$evo_stage %>%
  data.frame(Period = .) %>%
  group_by(Period) %>%
  summarise(Total = n())
#TTD
TTD_disease <- read.delim("/Target/TTD/P1-06-Target_disease.txt", header = FALSE, sep = "\t")
TTD_cancer<- TTD_disease[grepl("ICD-11: 2", TTD_disease$V5), ]
TTD_cancer<- TTD_cancer[grepl("Approv", TTD_cancer$V3), ]
TTD_target_uni<-read.delim("/Target/TTD/P2-01-TTD_uniprot_all.txt", header = FALSE, sep = "\t")
library(tidyverse)
TTD_target_uni <- TTD_target_uni %>%
  mutate(Group = cumsum(V1 == "TARGETID")) %>%
  filter(V1 %in% c("TARGETID", "UNIPROID", "TARGNAME", "TARGTYPE")) %>%
  pivot_wider(names_from = V1, values_from = V2) %>%
  select(-Group)
TTD_target_uni<-TTD_target_uni[TTD_target_uni$TARGETID%in%TTD_cancer$V1,]
TTD_target_gene<- strsplit(TTD_target_uni$TARGNAME, "(", fixed= T)
TTD_target_gene <- sapply(TTD_target_gene, function(x) x[2])
TTD_target_gene <- gsub("\\)", "", TTD_target_gene)
TTD_target_gene <- trimws(TTD_target_gene)
TTD_target_gene <- evo[evo$gene_name%in%TTD_target_gene,]
write.csv(TTD_target_uni,"/Target/TTD/TTD_cancer_app.csv",row.names = F)

target_counts <- TTD_target_gene$evo_stage %>%
  data.frame(Period = .) %>%
  group_by(Period) %>%
  summarise(Count = n())
target_counts$prob<-target_counts$Count / nrow(TTD_target_gene)


enrichment_data <- left_join(target_counts, background_counts, by = "Period")
results <- enrichment_data %>%
  rowwise() %>%
  mutate(
    Not_Biomarker = Total - Count,
    N_N = sum(target_counts$Count)-Count,
    Not_Background = sum(background_counts$Total) - Not_Biomarker - Count -N_N,
    Fisher_p_value = fisher.test(matrix(c(Count, N_N, Not_Biomarker, Not_Background), nrow = 2))$p.value,
    Fisher_odd=fisher.test(matrix(c(Count, N_N, Not_Biomarker, Not_Background), nrow = 2))$estimate,
  )
results=data.frame(results$Period,results$Fisher_p_value,results$Fisher_odd)

#DrugBank
library(data.table)
DrugBank<-fread("Target/drugbank/drug_target.txt",sep="\t",header = F)
colnames(DrugBank)<-c("drug_id","drug_name","drug_stage","drug_type","drug_description","drug_indication","target_id","uniprot_id","target_name","target_gene_name")
library(dplyr)
keywords <- c("Cancer", "Tumour", "Carcinoma", "Sarcoma", 
              "Malignancy", "Neoplasm", "Oncology", 
              "Metastasis", "Lump", "Leukemia", 
              "Lymphoma", "Melanoma")


pattern <- paste(keywords, collapse = "|")


DrugBank_cancer <- DrugBank %>%
  filter(grepl(pattern, drug_description,ignore.case = TRUE) | grepl(pattern, drug_indication,ignore.case = TRUE))
length(DrugBank_cancer$target_gene_name)
table(DrugBank_cancer$drug_stage)
DrugBank_cancer_app<-DrugBank_cancer %>%
  filter(grepl('approv', drug_stage,ignore.case = TRUE))
DrugBank_cancer<-DrugBank_cancer_app
#write.csv(DrugBank_cancer_app,"/Target/drugbank/DrugBank_cancer_app.csv",row.names = F)

DrugBank_target<-evo[evo$gene_name%in%DrugBank_cancer$target_gene_name,]
target_counts <- DrugBank_target$evo_stage %>%
  data.frame(Period = .) %>%
  group_by(Period) %>%
  summarise(Count = n())
target_counts$prob<-target_counts$Count / nrow(DrugBank_target)

enrichment_data <- left_join(target_counts, background_counts, by = "Period")
results <- enrichment_data %>%
  rowwise() %>%
  mutate(
    Not_Biomarker = Total - Count,
    N_N = sum(target_counts$Count)-Count,
    Not_Background = sum(background_counts$Total) - Not_Biomarker - Count -N_N,
    Fisher_p_value = fisher.test(matrix(c(Count, N_N, Not_Biomarker, Not_Background), nrow = 2))$p.value,
    Fisher_odd=fisher.test(matrix(c(Count, N_N, Not_Biomarker, Not_Background), nrow = 2))$estimate,
  )
results=data.frame(results$Period,results$Fisher_p_value,results$Fisher_odd)

#dgidb
dgidb<-fread("/Target/1_DGIDB_2022-Feb/interactions.tsv")
dgidb_cancer<-dgidb[grepl("cancer", dgidb$interaction_source_db_name, ignore.case = T),]
dgidb_cancer<-dgidb_cancer[dgidb_cancer$approved==TRUE,]
#write.csv(dgidb_cancer,"/drug/Target/1_DGIDB_2022-Feb/Dgidb_cancer_app.csv",row.names = F)
dgidb_cancer_target<-evo[evo$gene_name%in%dgidb_cancer$gene_name,]
target_counts <- dgidb_cancer_target$evo_stage %>%
  data.frame(Period = .) %>%
  group_by(Period) %>%
  summarise(Count = n())
target_counts$prob<-target_counts$Count / nrow(dgidb_cancer_target)

enrichment_data <- left_join(target_counts, background_counts, by = "Period")
results <- enrichment_data %>%
  rowwise() %>%
  mutate(
    Not_Biomarker = Total - Count,
    N_N = sum(target_counts$Count)-Count,
    Not_Background = sum(background_counts$Total) - Not_Biomarker - Count -N_N,
    Fisher_p_value = fisher.test(matrix(c(Count, N_N, Not_Biomarker, Not_Background), nrow = 2))$p.value,
    Fisher_odd=fisher.test(matrix(c(Count, N_N, Not_Biomarker, Not_Background), nrow = 2))$estimate,
  )
results=data.frame(results$Period,results$Fisher_p_value,results$Fisher_odd)

#union
Union_cancer_target <- rbind(DrugBank_target,TTD_target_gene,dgidb_cancer_target)
Union_cancer_target<-distinct(Union_cancer_target)
target_counts <- Union_cancer_target$evo_stage %>%
  data.frame(Period = .) %>%
  group_by(Period) %>%
  summarise(Count = n())
target_counts$prob<-target_counts$Count / nrow(Union_cancer_target)
enrichment_data <- left_join(target_counts, background_counts, by = "Period")
results <- enrichment_data %>%
  rowwise() %>%
  mutate(
    Not_Biomarker = Total - Count,
    N_N = sum(target_counts$Count)-Count,
    Not_Background = sum(background_counts$Total) - Not_Biomarker - Count -N_N,
    Fisher_p_value = fisher.test(matrix(c(Count, N_N, Not_Biomarker, Not_Background), nrow = 2))$p.value,
    Fisher_odd=fisher.test(matrix(c(Count, N_N, Not_Biomarker, Not_Background), nrow = 2))$estimate,
  )
results=data.frame(results$Period,results$Fisher_p_value,results$Fisher_odd)


result <- data.frame(target = character())
for (target in gene$V1) {
  matches <- DrugBank_cancer_app[grepl(target, DrugBank_cancer_app$target_gene_name, ignore.case = TRUE), ]
  
  if (nrow(matches) > 0) {
    result <- rbind(result, matches)
  }
}




# -------------------------------Target--------------------------
DrugBank<-fread("Target/drugbank/drug_target.txt",sep="\t",header = F)
colnames(DrugBank)<-c("drug_id","drug_name","drug_stage","drug_type","drug_description","drug_indication","target_id","uniprot_id","target_name","target_gene_name")
library(dplyr)
keywords <- c("Cancer", "Tumour", "Carcinoma", "Sarcoma", 
              "Malignancy", "Neoplasm", "Oncology", 
              "Metastasis", "Lump", "Leukemia", 
              "Lymphoma", "Melanoma")

pattern <- paste(keywords, collapse = "|")
DrugBank_cancer <- DrugBank %>%
  filter(grepl(pattern, drug_description,ignore.case = TRUE) | grepl(pattern, drug_indication,ignore.case = TRUE))
length(DrugBank_cancer$target_gene_name)
table(DrugBank_cancer$drug_stage)
DrugBank_cancer_app<-DrugBank_cancer %>%
  filter(grepl('approv', drug_stage,ignore.case = TRUE))
DrugBank_cancer<-DrugBank_cancer_app
#write.csv(DrugBank_cancer_app,"F:/postgraduate/database/drug/Target/drugbank/DrugBank_cancer_app.csv",row.names = F)

DrugBank_target<-evo[evo$gene_name%in%DrugBank_cancer$target_gene_name,]
DrugBank<-DrugBank[DrugBank$target_gene_name%in%DrugBank_target$gene_name,]

target_summary <- DrugBank %>%
  group_by(target_gene_name) %>%
  summarise(drug_count = n_distinct(drug_name)) %>%
  arrange(desc(drug_count))
a<-merge(target_summary,DrugBank_target,by.x="target_gene_name",by.y = "gene_name")
target_counts$prob<-target_counts$Count / nrow(DrugBank_target)