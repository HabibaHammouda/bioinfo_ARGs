#Bioinformatics Final Project - ARGs

#set working directory
setwd("C:/Users/user/Desktop/Bioinformatics_FinalProject_team9")

#load libraries (packages)
library(Biostrings)
library(seqinr)
library(readr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(ShortRead)
library(stringr)
library(pheatmap)

#convert fastq to fasta
fq <- readFastq("SRR797943.fastq")
writeXStringSet(sread(fq), "SRR797943.fasta")

#check the fasta file
reads <- readDNAStringSet("SRR797943.fasta")
reads

#building the database
system('"C:/Program Files/NCBI/blast-2.17.0+/bin/makeblastdb.exe" -in nucleotide_fasta_protein_homolog_model.fasta -dbtype nucl -out card_db')
file.exists("card_db.nsq")

#run BLAST to detect ARGs
system('"C:/Program Files/NCBI/blast-2.17.0+/bin/blastn.exe" -query SRR797943.fasta -db card_db -out results_SRR797943.txt -outfmt 6')
file.exists("results_SRR797943.txt")

file.info("results_SRR797943.txt")  # check file size
readLines("results_SRR797943.txt", n = 5)  # see first few lines

#reading blast results in R (tabulating in a data frame)
blast_results <- read.table("results_SRR797943.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(blast_results) <- c("query_id","subject_id","perc_identity","align_length",
                             "mismatches","gap_opens","q_start","q_end",
                             "s_start","s_end","evalue","bit_score")
head(blast_results)

#preprocessing
#filtering
filtered <- blast_results %>%
  filter(perc_identity >= 80,
         evalue < 1e-5,
         align_length / nchar(query_id) >= 0.7)  # approximate coverage
head(filtered)

#annotating and data manipulation
metadata <- read.delim("aro_categories_index.tsv", stringsAsFactors=FALSE) # Load CARD metadata

# extract the DNA accession (second field between pipes)
filtered <- filtered %>%
  mutate(DNA_accession = str_split(subject_id, "\\|", simplify = TRUE)[,2])
head(filtered$DNA_accession)

#merging
annotated <- merge(filtered, metadata, by.x = "DNA_accession", by.y = "DNA.Accession", all.x = TRUE)
head(annotated)

#Summaries
# Count the number of hits per drug class
drug_summary <- annotated %>%
  group_by(Drug.Class) %>%
  summarize(count = n()) %>%
  arrange(desc(count))
drug_summary

# Count the number of hits per gene family
gene_summary <- annotated %>%
  group_by(AMR.Gene.Family) %>%
  summarize(count = n()) %>%
  arrange(desc(count))
gene_summary

# Count the number of hits per % identity
identity_summary <- annotated %>%
  group_by(AMR.Gene.Family) %>%
  summarize(mean_identity = mean(perc_identity),
            min_identity = min(perc_identity),
            max_identity = max(perc_identity)) %>%
  arrange(desc(mean_identity))
identity_summary

#Visualization (bar chart and heat map)
# plotting the bar chart 
drug_summary <- drug_summary %>%
  mutate(Drug.Class_wrapped = str_wrap(Drug.Class, width = 20))

#to shorten the drug class name
#drug_summary <- drug_summary %>%
#mutate(Drug.Class_short = str_trunc(Drug.Class, width = 20, side = "right"))

ggplot(drug_summary, aes(x = reorder(Drug.Class_wrapped, count), y = count)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +  # horizontal bars
  theme_minimal() +
  labs(title = "ARG abundance by Drug Class",
       x = "Drug Class", y = "Number of Hits") +
  theme(axis.text.y = element_text(size = 10))

#######################################

#preparing info for heatmap
# Aggregate: take the max percent identity per query and gene family
identity_agg <- annotated %>%
  group_by(query_id, AMR.Gene.Family) %>%
  summarize(max_identity = max(perc_identity), .groups = "drop")
# Pivot to wide format
identity_matrix <- identity_agg %>%
  pivot_wider(names_from = AMR.Gene.Family,
              values_from = max_identity,
              values_fill = 0)
# Convert to matrix for heatmap
identity_mat <- as.matrix(identity_matrix[,-1])  # remove query_id column
rownames(identity_mat) <- identity_matrix$query_id

# Shorten Drug.Class or Gene Family names to 15 characters
colnames(identity_mat) <- str_trunc(colnames(identity_mat), width = 15, side = "right")
rownames(identity_mat) <- str_trunc(rownames(identity_mat), width = 15, side = "right")

top_families <- annotated %>%
  group_by(AMR.Gene.Family) %>%
  summarize(count = n()) %>%
  arrange(desc(count)) %>%
  slice_head(n = 20) %>%
  pull(AMR.Gene.Family)

identity_top <- identity_mat[, colnames(identity_mat) %in% top_families]

identity_mat[is.na(identity_mat)] <- 0
identity_mat <- apply(identity_mat, 2, as.numeric)
rownames(identity_mat) <- identity_matrix$query_id  # restore row names

#plotting the heatmap
pheatmap(identity_mat,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         fontsize_row = 6,
         fontsize_col = 6,
         main = "Percent Identity of ARG Hits")