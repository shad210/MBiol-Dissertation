#install packages
install.packages("ggplot2")
install.packages("ggpubr")
install.packages("scales")
install.packages("ggrepel")
install.packages("dplyr")
install.packages("BiocManager")
BiocManager::install("DESeq2", force = TRUE)
BiocManager::install("biomaRt")
install.packages("reshape2")
install.packages("car")

#loading packages
library(ggplot2)
library(ggpubr)
library(scales)
library(ggrepel)
library(dplyr)
library(BiocManager)
library(DESeq2)
library(biomaRt)
options(timeout = 60)
library(reshape2)
library(car)


#reading counts .rds files
bulk_CD8_counts <- readRDS("bulk_CD8_counts.rds")
healthydonor_meta <- readRDS("healthydonor_meta.rds")
healthydonor_mixcr <- readRDS("healthydonor_mixcr.rds")
patient_meta <- readRDS("patient_meta.rds")
patient_mixcr <- readRDS("patient_mixcr.rds")

#add age data to patient and healthy donor MixCR

#extract age
patient_age <- patient_meta %>%
  distinct(patient, age)

colnames(healthydonor_meta)[colnames(healthydonor_meta) == "patient"] <- "donor"
healthydonor_age <- healthydonor_meta %>%
  distinct(donor, age)

#ensure age is a character
patient_age <- patient_age %>%
  mutate(patient = as.character(patient))

healthydonor_age <- healthydonor_age %>%
  mutate(donor = as.character(donor))

#append age data onto MixCR dataframes
patient_mixcr <- left_join(patient_mixcr, patient_age, by = "patient")

healthydonor_mixcr <- left_join(healthydonor_mixcr, healthydonor_age, by = "donor")

#determine which cells are MAITs
patient_mixcr <- patient_mixcr %>%
  mutate(MAITness = ifelse(bestVHit == "TRAV1-2" & bestJHit %in% c("TRAJ12", "TRAJ20", "TRAJ33"), "MAIT", "nonMAIT"))

healthydonor_mixcr <- healthydonor_mixcr %>%
  mutate(MAITness = ifelse(bestVHit == "TRAV1-2" & bestJHit %in% c("TRAJ12", "TRAJ20", "TRAJ33"), "MAIT", "nonMAIT"))

#subset MAITs into a new dataframe
patient_MAITs <- patient_mixcr %>%
  filter(MAITness == "MAIT")

healthydonor_MAITs <- healthydonor_mixcr %>%
  filter(MAITness == "MAIT")

#subset pre-treatment MAITs into a new dataframe for patients
patient_MAITs_C1 <- patient_MAITs %>%
  filter(cycle == "C1")

#group the patient MAIT data pretreatment by total clone count per patient, and plot compared to age
##first get MAIT totals
patient_MAITs_C1_totals <- patient_MAITs_C1 %>%
  group_by(patient, age) %>%
  summarise(total_clone_count = sum(cloneCount), .groups = "drop")

##fit a linear model to calculate R-squared
model_patient_MAITs_C1 <- lm(total_clone_count ~ age, data = patient_MAITs_C1_totals)
r_squared_patient_MAITs_C1 <- summary(model_patient_MAITs_C1)$r.squared

##scatterplot with trendline and R-squared
ggplot(patient_MAITs_C1_totals, aes(x = age, y = total_clone_count)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  annotate("text", x = 80, y = 1750,
           label = paste("R² =", round(r_squared_patient_MAITs_C1, 2)), size = 5, color = "darkred") +
  labs(title = "MAIT counts vs Age", x = "Age", y = "Total Clone Count")

#get MAIT proportions for each patient
patient_mixcr_C1 <- patient_mixcr %>%
  filter(cycle == "C1")

patient_C1_proportions <- patient_mixcr_C1 %>%
  group_by(patient, age) %>%
  summarise(
    total_clone_count = sum(cloneCount),
    mait_clone_count = sum(cloneCount[MAITness == "MAIT"]),
    .groups = "drop"
  )

patient_C1_proportions <- patient_C1_proportions %>%
  mutate(MAIT_proportion = mait_clone_count / total_clone_count)

#fit a linear model to this too
model_patient_C1_proportions <- lm(MAIT_proportion ~ age, data = patient_C1_proportions)
cor.test(log(patient_C1_proportions$MAIT_proportion), patient_C1_proportions$age, method = "spearman")
r_squared_patient_C1_proportions <- summary(model_patient_C1_proportions)$r.squared

#############

#plot MAIT proportion vs age
ggplot(patient_C1_proportions, aes(x = age, y = MAIT_proportion)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  scale_y_log10() +
  annotate("text", x = 80, y = 0.15,
           label = paste("R² =", round(r_squared_patient_C1_proportions, 2)), size = 5, color = "darkred") +
  labs(title = "MAIT proportion vs Age", x = "Age", y = "MAIT proportion")


#let's try doing number of MAIT rows compared to total rows per patient - read depths could be causing the lack of significance we are witnessing
patient_mixcr_C1_rows <- patient_mixcr_C1 %>%
  group_by(patient, age) %>%
  summarise(
    total_rows = n(),
    mait_rows = sum(MAITness == "MAIT"),
    .groups = "drop"
  )

patient_mixcr_C1_rows <- patient_mixcr_C1_rows %>%
  mutate(MAIT_proportion = mait_rows / total_rows)

patient_mixcr_C1_rows_proportions <- patient_mixcr_C1_rows %>%
  select(patient, age, MAIT_proportion)

#now plot this with linear regression
model_patient_C1_rows_proportions <- lm(MAIT_proportion ~ age, data = patient_mixcr_C1_rows_proportions)
r_squared_patient_C1_rows_proportions <- summary(model_patient_C1_rows_proportions)$r.squared

ggplot(patient_mixcr_C1_rows_proportions, aes(x = age, y = MAIT_proportion)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  annotate("text", x = 80, y = 0.04,
           label = paste("R² =", round(r_squared_patient_C1_rows_proportions, 2)), size = 5, color = "darkred") +
  labs(title = "MAIT row proportion vs Age", x = "Age", y = "MAIT row proportion")


#######DESeq
ENSG <- bulk_CD8_counts$ENSG
bulk_CD8_counts <- bulk_CD8_counts[, c(ncol(bulk_CD8_counts), 1:(ncol(bulk_CD8_counts) - 1))]
patient_meta$sample <- paste0(patient_meta$sample, "_")
colnames(bulk_CD8_counts) <- gsub("CD8_None", "", colnames(bulk_CD8_counts))
numeric_colnames <- as.numeric(gsub("_.*", "", colnames(bulk_CD8_counts)[-1]))
remaining_cols <- bulk_CD8_counts[-1][order(numeric_colnames)]
bulk_CD8_counts <- cbind(ENSG, remaining_cols)

HD_cols <- grepl("HD", colnames(bulk_CD8_counts))
HD_counts <- bulk_CD8_counts[, HD_cols]
HD_counts <- cbind(ENSG = bulk_CD8_counts$ENSG, HD_counts)
patient_counts <- bulk_CD8_counts[, !HD_cols]

C1_samples <- grep("_C1_", colnames(bulk_CD8_counts), value = TRUE)
C1_CD8_data <- bulk_CD8_counts[, C1_samples]
C1_CD8_data <- C1_CD8_data[, order(as.numeric(colnames(C1_CD8_data)))]
C1_CD8_data <- C1_CD8_data[, 1:(ncol(C1_CD8_data) - 4)]
C1_CD8_data <- cbind(ENSG = bulk_CD8_counts$ENSG, C1_CD8_data)

##make the metadata (this is for C1)
C1_meta <- patient_meta %>% filter(grepl("_C1_", sample))
C1_meta_filtered <- C1_meta[C1_meta$sample %in% colnames(C1_CD8_data), ]

patient_meta_filtered <- patient_meta[patient_meta$sample %in% colnames(patient_counts), ]

samples_in_patient_meta <- patient_meta_filtered$sample

patient_counts_filtered <- patient_counts[, colnames(patient_counts) %in% samples_in_patient_meta]
patient_counts_filtered <- cbind(ENSG, patient_counts_filtered)
row.names(patient_counts_filtered) <- patient_counts_filtered[, 1]
patient_counts_filtered <- patient_counts_filtered[, -1]
#convert ENSG to gene names
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
gene_mapping <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                      filters = "ensembl_gene_id",
                      values = ENSG,
                      mart = ensembl)
head(gene_mapping)

patient_C1_proportions <- patient_C1_proportions[!grepl("\\*", patient_C1_proportions$patient), ] #cleaning
patient_C1_proportions <- patient_C1_proportions[order(as.numeric(gsub("_.*", "", patient_C1_proportions$patient))), ] #ordering numerically

patient_C1_proportions_filtered <- patient_C1_proportions[patient_C1_proportions$patient %in% 
                                                            C1_meta_filtered$patient, ] #more cleaning to ensure same number of rows and columns for DESeq

C1_meta_filtered_combined <- cbind(C1_meta_filtered, 
                                   patient_C1_proportions_filtered$total_clone_count,
                                   patient_C1_proportions_filtered$mait_clone_count,
                                   patient_C1_proportions_filtered$MAIT_proportion)

C1_meta_filtered_combined <- C1_meta_filtered_combined %>%
  mutate(age_group = case_when(
    age < 59 ~ "young",
    age >= 59 & age < 72 ~ "mid",
    age >= 72 ~ "old"
  ))

C1_meta_filtered_combined <- C1_meta_filtered_combined %>%
  rename(`patient_C1_proportions_filtered$total_clone_count` = "total_clone_count", 
         `patient_C1_proportions_filtered$mait_clone_count` = "mait_clone_count", 
         `patient_C1_proportions_filtered$MAIT_proportion` = "MAIT_proportion")
C1_meta_filtered_combined$nonMAIT_clone_count <- C1_meta_filtered_combined$total_clone_count - C1_meta_filtered_combined$mait_clone_count
C1_meta_filtered_combined$nonMAIT_proportion <- 1 - C1_meta_filtered_combined$MAIT_proportion
C1_meta_filtered_combined$log_MAIT <- log(C1_meta_filtered_combined$MAIT_proportion)
C1_meta_filtered_combined$log_nonMAIT <- log(C1_meta_filtered_combined$nonMAIT_proportion)


matching_columns_C1 <- colnames(C1_CD8_data) %in% C1_meta_filtered_combined$sample

C1_CD8_filtered <- C1_CD8_data[, c(TRUE, matching_columns_C1[-1])]
rownames(C1_CD8_filtered) <- C1_CD8_filtered$ENSG
C1_CD8_filtered <- C1_CD8_filtered[, -which(names(C1_CD8_filtered) == "ENSG")]

dds_C1_age <- DESeqDataSetFromMatrix(countData = C1_CD8_filtered,
                                     colData = C1_meta_filtered_combined,
                                     design = ~log(nonMAIT_proportion) + log(MAIT_proportion)*scale(age))

dds_C1_age <- DESeq(dds_C1_age)

res_dds_C1_age <- results(dds_C1_age)

res_C1_age_df <- as.data.frame(res_dds_C1_age)

res_C1_age_df <- res_C1_age_df[!is.na(res_C1_age_df$padj) & !is.na(res_C1_age_df$log2FoldChange), ]

res_C1_age_df$significance <- "Not Sig"
res_C1_age_df$significance[res_C1_age_df$padj < 0.05] <- "Sig"

res_C1_age_df$ensembl_gene_id <- rownames(res_C1_age_df)

res_C1_age_df$hgnc_symbol <- gene_mapping$hgnc_symbol[match(res_C1_age_df$ensembl_gene_id, gene_mapping$ensembl_gene_id)]

View(res_C1_age_df)

##same thing but with age as age group factor instead of scaled numerical age
dds_C1_age_group <- DESeqDataSetFromMatrix(countData = C1_CD8_filtered,
                                           colData = C1_meta_filtered_combined,
                                           design = ~log(MAIT_proportion)*as.factor(age_group))

dds_C1_age_group <- DESeq(dds_C1_age_group)

res_dds_C1_age_group <- results(dds_C1_age_group)

res_C1_age_group_df <- as.data.frame(res_dds_C1_age_group)

res_C1_age_group_df <- res_C1_age_group_df[!is.na(res_C1_age_group_df$padj) & !is.na(res_C1_age_group_df$log2FoldChange), ]

res_C1_age_group_df$significance <- "Not Sig"
res_C1_age_group_df$significance[res_C1_age_group_df$padj < 0.05 & abs(res_C1_age_group_df$log2FoldChange) > 0.58] <- "Sig"

View(res_C1_age_group_df)

###Healthy donors
healthydonor_proportions <- healthydonor_mixcr %>%
  group_by(donor, age) %>%
  summarise(
    total_clone_count = sum(cloneCount),
    mait_clone_count = sum(cloneCount[MAITness == "MAIT"]),
    .groups = "drop"
  )

healthydonor_proportions <- healthydonor_proportions %>%
  mutate(MAIT_proportion = mait_clone_count / total_clone_count)

rownames(HD_counts) <- HD_counts[ ,1]
HD_counts <- HD_counts[, -1]

healthydonor_proportions <- healthydonor_proportions %>%
  mutate(sample = paste0(donor, "_HD_"))

matching_columns_HD <- colnames(HD_counts) %in% healthydonor_proportions$sample
healthydonor_proportions_filtered <- healthydonor_proportions[healthydonor_proportions$sample %in% colnames(HD_counts), ]

healthydonor_proportions_filtered <- healthydonor_proportions_filtered %>%
  mutate(age_group = case_when(
    age < 59 ~ "young",
    age >= 59 & age < 72 ~ "mid",
    age >= 72 ~ "old"
  ))

#DESeq for healthy
dds_HD_age <- DESeqDataSetFromMatrix(countData = HD_counts,
                                     colData = healthydonor_proportions_filtered,
                                     design = ~log(MAIT_proportion)*scale(age))

dds_HD_age <- DESeq(dds_HD_age)

res_dds_HD_age <- results(dds_HD_age)

res_HD_age_df <- as.data.frame(res_dds_HD_age)

res_HD_age_df <- res_HD_age_df[!is.na(res_HD_age_df$padj) & !is.na(res_HD_age_df$log2FoldChange), ]

res_HD_age_df$significance <- "Not Sig"
res_HD_age_df$significance[res_HD_age_df$padj < 0.05] <- "Sig"

res_HD_age_df$ensembl_gene_id <- rownames(res_HD_age_df)

res_HD_age_df$hgnc_symbol <- gene_mapping$hgnc_symbol[match(res_HD_age_df$ensembl_gene_id, gene_mapping$ensembl_gene_id)]

View(res_HD_age_df)

#DESeq for healthy with age groups
dds_HD_age_group <- DESeqDataSetFromMatrix(countData = HD_counts,
                                     colData = healthydonor_proportions_filtered,
                                     design = ~log(MAIT_proportion)*as.factor(age_group))

dds_HD_age_group <- DESeq(dds_HD_age_group)

res_dds_HD_age_group <- results(dds_HD_age_group)

res_HD_age_group_df <- as.data.frame(res_dds_HD_age_group)

res_HD_age_group_df <- res_HD_age_group_df[!is.na(res_HD_age_group_df$padj) & !is.na(res_HD_age_group_df$log2FoldChange), ]

res_HD_age_group_df$significance <- "Not Sig"
res_HD_age_group_df$significance[res_HD_age_group_df$padj < 0.05 & abs(res_HD_age_group_df$log2FoldChange) > 0.58] <- "Sig"



View(res_HD_age_group_df) #no significance?


#####################DESeq with healthy and C1 combined
colnames(healthydonor_proportions_filtered)[colnames(healthydonor_proportions_filtered) == "donor"] <- "patient" #rename donor to patient so that the column names will match
healthydonor_proportions_filtered$patient <- as.integer(healthydonor_proportions_filtered$patient)
healthydonor_meta$sample <- paste0(healthydonor_meta$sample, "_")
healthydonor_proportions_filtered <- merge(healthydonor_proportions_filtered, healthydonor_meta, by = "sample")
healthydonor_proportions_filtered$age <- healthydonor_proportions_filtered$age.x
healthydonor_proportions_filtered$age.x <- NULL
healthydonor_proportions_filtered$age.y <- NULL
C1_meta_filtered_combined$sex <- ifelse(C1_meta_filtered_combined$sex == 1, "M", "F")


combined_filtered <- bind_rows(C1_meta_filtered_combined, healthydonor_proportions_filtered)
combined_filtered$status <- ifelse(combined_filtered$patient > 525, "HD", "Patient")
combined_filtered$nonMAIT_proportion <- 1 - combined_filtered$MAIT_proportion

matching_columns_combined <- colnames(bulk_CD8_counts) %in% combined_filtered$sample

CD8_combined_filtered <- bulk_CD8_counts[, c(TRUE, matching_columns_combined[-1])]
rownames(CD8_combined_filtered) <- CD8_combined_filtered[,1]  # Set the first column as row names
CD8_combined_filtered <- CD8_combined_filtered[,-1]

dds_combined <- DESeqDataSetFromMatrix(countData = CD8_combined_filtered,
                                     colData = combined_filtered,
                                     design = ~as.factor(status) + log(nonMAIT_proportion) + 
                                       log(MAIT_proportion)*as.factor(age_group))

dds_combined <- DESeq(dds_combined)

res_dds_combined <- results(dds_combined)

res_dds_combined_df <- as.data.frame(res_dds_combined)

res_dds_combined_df <- res_dds_combined_df[!is.na(res_dds_combined_df$padj) & !is.na(res_dds_combined_df$log2FoldChange), ]

res_dds_combined_df$significance <- "Not Sig"
res_dds_combined_df$significance[res_dds_combined_df$padj < 0.05] <- "Sig"

res_dds_combined_df$ensembl_gene_id <- rownames(res_dds_combined_df)

res_dds_combined_df$hgnc_symbol <- gene_mapping$hgnc_symbol[match(res_dds_combined_df$ensembl_gene_id, gene_mapping$ensembl_gene_id)]

View(res_dds_combined_df)

####just MAITs
dds_MAIT <- DESeqDataSetFromMatrix(countData = CD8_combined_filtered,
                                       colData = combined_filtered,
                                       design = ~log(MAIT_proportion))

dds_MAIT <- DESeq(dds_MAIT)

res_dds_MAIT <- results(dds_MAIT)

res_dds_MAIT_df <- as.data.frame(res_dds_MAIT)

res_dds_MAIT_df <- res_dds_MAIT_df[!is.na(res_dds_MAIT_df$padj) & !is.na(res_dds_MAIT_df$log2FoldChange), ]

res_dds_MAIT_df$significance <- "Not Sig"
res_dds_MAIT_df$significance[res_dds_MAIT_df$padj < 0.05] <- "Sig"

res_dds_MAIT_df$ensembl_gene_id <- rownames(res_dds_MAIT_df)

res_dds_MAIT_df$hgnc_symbol <- gene_mapping$hgnc_symbol[match(res_dds_MAIT_df$ensembl_gene_id, gene_mapping$ensembl_gene_id)]

View(res_dds_MAIT_df)

################Age plot
model_proportions <- lm(MAIT_proportion ~ age, data = combined_filtered)
spearman_test <- cor.test(log(combined_filtered$MAIT_proportion), combined_filtered$age, method = "spearman")

rho_spearman <- round(spearman_test$estimate, 3)
p_value_spearman <- signif(spearman_test$p.value, 3)

#ggplot of MAIT vs age
ggplot(combined_filtered, aes(x = age, y = MAIT_proportion)) + 
  geom_point(alpha = 0.7, color = "blue") +
  geom_smooth(method = "lm", se = TRUE, color = "black", fill = "grey70", linetype = "dashed") +
  scale_y_log10(labels = label_number()) +
  labs(title = "MAIT proportion vs Age",
       x = "Age",
       y = "MAIT proportion") +
  annotate("text", x = max(combined_filtered$age) * 0.8, y = max(combined_filtered$MAIT_proportion) * 0.5,
           label = paste0("Rho = ", rho_spearman, "\n p = ", p_value_spearman),
           color = "black", size = 5, hjust = 0) +
  theme_minimal()

###comparison
combined_filtered <- combined_filtered %>%
  mutate(sample_type = ifelse(grepl("_C1_", sample), "C1", "HD"))

summary_combined_filtered <- combined_filtered %>%
  group_by(sample_type, age_group) %>%
  summarise(
    mean_MAIT = mean(MAIT_proportion, na.rm = TRUE),
    se_MAIT = sd(MAIT_proportion, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

summary_combined_filtered <- summary_combined_filtered %>%
  mutate(group = paste(sample_type, age_group, sep = " - "))

ggplot() +
  geom_bar(data = summary_combined_filtered,
           aes(x = age_group, y = mean_MAIT, fill = sample_type),
           stat = "identity",
           position = position_dodge(width = 0.7),
           width = 0.6,
           alpha = 0.3) +
  geom_errorbar(data = summary_combined_filtered,
                aes(x = age_group, ymin = mean_MAIT - se_MAIT, ymax = mean_MAIT + se_MAIT, group = sample_type),
                position = position_dodge(width = 0.7),
                width = 0.2) +
  geom_jitter(data = combined_filtered,
              aes(x = age_group, y = MAIT_proportion, color = sample_type),
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.7),
              size = 2, alpha = 0.7) +
  labs(
    title = "MAIT proportion Cycle 1 patients vs Healthy Donors",
    x = "Age group",
    y = "MAIT proportion"
  ) +
  scale_fill_brewer(palette = "Set3") +
  scale_color_brewer(palette = "Set3") +
  coord_cartesian(ylim = c(0, 0.1)) +
  theme_minimal()
###volcano plots

#MAIT
top_genes_MAIT <- res_dds_MAIT_df %>%
  arrange(padj) %>%
  head(12)

res_dds_MAIT_df <- res_dds_MAIT_df %>%
  mutate(diffexpressed = ifelse(log2FoldChange > 0, "Upregulated", "Downregulated"))

res_dds_MAIT_df$diffexpressed <- "No"
res_dds_MAIT_df$diffexpressed[res_dds_MAIT_df$log2FoldChange > 0.58 & res_dds_MAIT_df$padj < 0.05] <- "Upregulated"
res_dds_MAIT_df$diffexpressed[res_dds_MAIT_df$log2FoldChange < -0.58 & res_dds_MAIT_df$padj < 0.05] <- "Downregulated"

ggplot(data = res_dds_MAIT_df, aes(x = log2FoldChange, y = -log10(padj), color = diffexpressed)) +
  geom_point(alpha = 0.7) +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue")) +
  geom_text_repel(data = top_genes_MAIT, aes(label = hgnc_symbol, color = "black"), size = 4, max.overlaps = 10) +
  theme_minimal() +
  geom_vline(xintercept = c(-0.58, 0.58), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
  labs(title = "MAIT differential gene expression",
       x = "Log2FoldChange",
       y = "-Log10padj",
       color = "Gene regulation") +
  theme(legend.position = "top")

#Age-related
top_genes_age <- res_dds_combined_df %>%
  arrange(padj) %>%
  head(12)

res_dds_combined_df <- res_dds_combined_df %>%
  mutate(diffexpressed = ifelse(log2FoldChange > 0, "Upregulated", "Downregulated"))

res_dds_combined_df$diffexpressed <- "No"
res_dds_combined_df$diffexpressed[res_dds_combined_df$log2FoldChange > 0.01 & res_dds_combined_df$padj < 0.05] <- "Upregulated"
res_dds_combined_df$diffexpressed[res_dds_combined_df$log2FoldChange < -0.01 & res_dds_combined_df$padj < 0.05] <- "Downregulated"

ggplot(data = res_dds_combined_df, aes(x = log2FoldChange, y = -log10(padj), color = diffexpressed)) +
  geom_point(alpha = 0.7) +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue")) +
  geom_text_repel(data = top_genes_age, aes(label = hgnc_symbol, color = "black"), size = 4, max.overlaps = 20) +
  theme_minimal() +
  geom_vline(xintercept = c(-0.01, 0.01), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
  labs(title = "MAIT age-related differential gene expression",
       x = "Log2FoldChange",
       y = "-Log10padj",
       color = "Gene regulation") +
  theme(legend.position = "top")

#########expression graphs
bulk_subset <- bulk_CD8_counts[, colnames(bulk_CD8_counts) == "ENSG" | grepl("_C1_|_HD_", colnames(bulk_CD8_counts)), drop = FALSE]
bulk_subset <- bulk_subset[, -((ncol(bulk_subset)-3):ncol(bulk_subset))]
rownames(bulk_subset) <- bulk_subset$ENSG
bulk_subset <- bulk_subset[, !(colnames(bulk_subset) %in% c("ENSG"))]
cols_to_keep <- colnames(bulk_subset) %in% combined_filtered$sample
bulk_subset_filtered <- bulk_subset[, cols_to_keep]

baselineNormaliseTransform <- function(count_data, sample_metadata){
  allddsHTSeq<-DESeqDataSetFromMatrix(countData = count_data,
                                      colData = sample_metadata,
                                      design=~1)
  allddsHTSeq<-estimateSizeFactors(allddsHTSeq)
  allddsHTSeq<-estimateDispersions(allddsHTSeq)
  data  <- counts(allddsHTSeq, normalized = TRUE)
  vsd<-getVarianceStabilizedData(allddsHTSeq)
  toKeep  <- rowMeans(data) > 10
  data <- vsd[toKeep,]
}

normalized_bulk_subset_filtered <- baselineNormaliseTransform(count_data = bulk_subset_filtered, sample_metadata = combined_filtered)
normalized_bulk_df <- as.data.frame(normalized_bulk_subset_filtered)

normalized_bulk_df$ENSG <- rownames(normalized_bulk_df)
normalized_bulk_df$hgnc_symbol <- gene_mapping$hgnc_symbol[match(normalized_bulk_df$ENSG, gene_mapping$ensembl_gene_id)]

normalized_bulk_df_transposed <- t(normalized_bulk_df)

expression_df <- cbind(combined_filtered, normalized_bulk_df_transposed)

####ggplots for this
#RORC
ggplot(expression_df, aes(x = MAIT_proportion, y = log2(ENSG00000143365))) +
         geom_point(alpha = 0.7) +
         scale_x_log10(labels = label_number()) +
         scale_y_log10(labels = label_number()) +
         geom_smooth(method = "lm", se = TRUE, linetype = "dashed", color = "black") +
         stat_cor(method = "spearman", label.x = log10(0.03), label.y = 0.05) +
         theme_minimal() +
         labs(title = "RORC",
              x = "MAIT proportion",
              y = "Log2Expression") +
         theme(legend.position = "none")

#SLC4A10
ggplot(expression_df, aes(x = MAIT_proportion, y = log2(ENSG00000144290))) +
  geom_point(alpha = 0.7) +
  scale_x_log10(labels = label_number()) +
  scale_y_log10(labels = label_number()) +
  geom_smooth(method = "lm", se = TRUE, linetype = "dashed", color = "black") +
  stat_cor(method = "spearman", label.x = log10(0.03), label.y = 0.05) +
  theme_minimal() +
  labs(title = "SLC4A10",
       x = "MAIT proportion",
       y = "Log2Expression") +
  theme(legend.position = "none")

#IL23R
ggplot(expression_df, aes(x = MAIT_proportion, y = log2(ENSG00000162594))) +
  geom_point(alpha = 0.7) +
  scale_x_log10(labels = label_number()) +
  scale_y_log10(labels = label_number()) +
  geom_smooth(method = "lm", se = TRUE, linetype = "dashed", color = "black") +
  stat_cor(method = "spearman", label.x = log10(0.03), label.y = 0.05) +
  theme_minimal() +
  labs(title = "IL23R",
       x = "MAIT proportion",
       y = "Log2Expression") +
  theme(legend.position = "none")

#ME1
ggplot(expression_df, aes(x = MAIT_proportion, y = log2(ENSG00000065833))) +
  geom_point(alpha = 0.7) +
  scale_x_log10(labels = label_number()) +
  scale_y_log10(labels = label_number()) +
  geom_smooth(method = "lm", se = TRUE, linetype = "dashed", color = "black") +
  stat_cor(method = "spearman", label.x = log10(0.03), label.y = 0.05) +
  theme_minimal() +
  labs(title = "ME1",
       x = "MAIT proportion",
       y = "Log2Expression") +
  theme(legend.position = "none")

#PRSS35
ggplot(expression_df, aes(x = MAIT_proportion, y = log2(ENSG00000146250))) +
  geom_point(alpha = 0.7) +
  scale_x_log10(labels = label_number()) +
  geom_smooth(method = "lm", se = TRUE, linetype = "dashed", color = "black") +
  stat_cor(method = "spearman", label.x = log10(0.03), label.y = 5) +
  theme_minimal() +
  labs(title = "PRSS35",
       x = "MAIT proportion",
       y = "Log2Expression") +
  theme(legend.position = "none")

#LTK
ggplot(expression_df, aes(x = MAIT_proportion, y = log2(ENSG00000062524))) +
  geom_point(alpha = 0.7) +
  scale_x_log10(labels = label_number()) +
  geom_smooth(method = "lm", se = TRUE, linetype = "dashed", color = "black") +
  stat_cor(method = "spearman", label.x = log10(0.03), label.y = 5) +
  theme_minimal() +
  labs(title = "LTK",
       x = "MAIT proportion",
       y = "Log2Expression") +
  theme(legend.position = "none")

#MATN2
ggplot(expression_df, aes(x = MAIT_proportion, y = log2(ENSG00000132561))) +
  geom_point(alpha = 0.7) +
  scale_x_log10(labels = label_number()) +
  geom_smooth(method = "lm", se = TRUE, linetype = "dashed", color = "black") +
  stat_cor(method = "spearman", label.x = log10(0.03), label.y = 5) +
  theme_minimal() +
  labs(title = "MATN2",
       x = "MAIT proportion",
       y = "Log2Expression") +
  theme(legend.position = "none")

#COL5A1
ggplot(expression_df, aes(x = MAIT_proportion, y = log2(ENSG00000130635))) +
  geom_point(alpha = 0.7) +
  scale_x_log10(labels = label_number()) +
  geom_smooth(method = "lm", se = TRUE, linetype = "dashed", color = "black") +
  stat_cor(method = "spearman", label.x = log10(0.03), label.y = 5) +
  theme_minimal() +
  labs(title = "COL5A1",
       x = "MAIT proportion",
       y = "Log2Expression") +
  theme(legend.position = "none")

#CCR6
ggplot(expression_df, aes(x = MAIT_proportion, y = log2(ENSG00000112486))) +
  geom_point(alpha = 0.7) +
  scale_x_log10(labels = label_number()) +
  scale_y_log10(labels = label_number()) +
  geom_smooth(method = "lm", se = TRUE, linetype = "dashed", color = "black") +
  stat_cor(method = "spearman", label.x = log10(0.03), label.y = 0.05) +
  theme_minimal() +
  labs(title = "CCR6",
       x = "MAIT proportion",
       y = "Log2Expression") +
  theme(legend.position = "none")

#ADAM12
ggplot(expression_df, aes(x = MAIT_proportion, y = log2(ENSG00000148848))) +
  geom_point(alpha = 0.7) +
  scale_x_log10(labels = label_number()) +
  scale_y_log10(labels = label_number()) +
  geom_smooth(method = "lm", se = TRUE, linetype = "dashed", color = "black") +
  stat_cor(method = "spearman", label.x = log10(0.03), label.y = 0.05) +
  theme_minimal() +
  labs(title = "ADAM12",
       x = "MAIT proportion",
       y = "Log2Expression") +
  theme(legend.position = "none")

#NTN4
ggplot(expression_df, aes(x = MAIT_proportion, y = log2(ENSG00000074527))) +
  geom_point(alpha = 0.7) +
  scale_x_log10(labels = label_number()) +
  scale_y_log10(labels = label_number()) +
  geom_smooth(method = "lm", se = TRUE, linetype = "dashed", color = "black") +
  stat_cor(method = "spearman", label.x = log10(0.03), label.y = 0.05) +
  theme_minimal() +
  labs(title = "NTN4",
       x = "MAIT proportion",
       y = "Log2Expression") +
  theme(legend.position = "none")

#GPR15
ggplot(expression_df, aes(x = MAIT_proportion, y = log2(ENSG00000154165))) +
  geom_point(alpha = 0.7) +
  scale_x_log10(labels = label_number()) +
  scale_y_log10(labels = label_number()) +
  geom_smooth(method = "lm", se = TRUE, linetype = "dashed", color = "black") +
  stat_cor(method = "spearman", label.x = log10(0.03), label.y = 0.05) +
  theme_minimal() +
  labs(title = "GPR15",
       x = "MAIT proportion",
       y = "Log2Expression") +
  theme(legend.position = "none")

#USP13
ggplot(expression_df, aes(x = MAIT_proportion, y = log2(ENSG00000058056))) +
  geom_point(alpha = 0.7) +
  scale_x_log10(labels = label_number()) +
  scale_y_log10(labels = label_number()) +
  geom_smooth(method = "lm", se = TRUE, linetype = "dashed", color = "black") +
  stat_cor(method = "spearman", label.x = log10(0.03), label.y = 0.1) +
  theme_minimal() +
  labs(title = "USP13",
       x = "MAIT proportion",
       y = "Log2Expression") +
  theme(legend.position = "none")