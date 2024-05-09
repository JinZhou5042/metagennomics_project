library(tibble)
library(stringr)
library(phyloseq)
library(Biostrings)
library(dplyr)
library(vegan)
library(ggplot2)

output_data_path <- "/scratch365/zhuang8/R16S_data/workspace/output_data"
output_plot_path <- "/scratch365/zhuang8/R16S_data/workspace/plot_data"
seqtab.nochim <- readRDS(file.path(output_data_path, "seqtab.nochim.Rds"))
taxa <- readRDS(file.path(output_data_path, "taxa.Rds"))

# check basic info
# dim(seqtab.nochim)
# rownames(seqtab.nochim)     # row names of seqtab.nochim
# colnames(seqtab.nochim)     # column names of seqtab.nochim
# seqtab.nochim[, 1]          # first column of seqtab.nochim
# seqtab.nochim[1, ]          # first row of seqtab.nochim


# Extract the row names from a data table `seqtab.nochim` which likely contains sequence data
samples.out <- rownames(seqtab.nochim)
metadata <- tibble(Sample_names = samples.out) %>%
  mutate(Group = ifelse(str_detect(Sample_names, "UG"), "Group_UG", "Group_G")) %>%
  column_to_rownames(var = "Sample_names")
#       > metadata
#                             Group
#       24-E5-G-B1-surf   Group_G
#       24-E5-G-B3-50cm   Group_G
#       24-E5-G-B3-surf   Group_G
#       24-E5-G-UB-50cm   Group_G
#       24-E5-G-UB-surf   Group_G
#       24-E5-UG-B1-surf Group_UG
#       24-E5-UG-B3-50cm Group_UG
#       24-E5-UG-B3-surf Group_UG
#       24-E5-UG-UB-50cm Group_UG
#       24-E5-UG-UB-surf Group_UG

ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(metadata), 
               tax_table(taxa))
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq_along(taxa_names(ps)))
#        > ps
#        phyloseq-class experiment-level object
#        otu_table()   OTU Table:         [ 1543 taxa and 10 samples ]
#        sample_data() Sample Data:       [ 10 samples by 1 sample variables ]
#        tax_table()   Taxonomy Table:    [ 1543 taxa by 7 taxonomic ranks ]
#        refseq()      DNAStringSet:      [ 1543 reference sequences ]

# Show available ranks in the dataset
rank_names(ps)
# Create table, number of features for each phyla
table(tax_table(ps)[, "Phylum"], exclude = NULL)
# Filter out NA and unwanted phyla
ps <- subset_taxa(ps, !is.na(Phylum) & !Phylum %in% c(""))

# Compute prevalence of each feature, store as data.frame
prevdf <- apply(X = otu_table(ps),
               MARGIN = ifelse(taxa_are_rows(ps), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevdf <- data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(ps),
                    tax_table(ps))
# Compute the total and average prevalences of each feature
plyr::ddply(prevdf, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})
# Define phyla to filter (low counts)
filterPhyla <- c("")
# Filter entries with unidentified Phylum.
ps <- subset_taxa(ps, !Phylum %in% filterPhyla)
# Subset to the remaining phyla
prevdf1 <- subset(prevdf, Phylum %in% get_taxa_unique(ps, "Phylum"))

# Plot
plot <- ggplot(prevdf1, aes(TotalAbundance, Prevalence / nsamples(ps), color = Phylum)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) +
  scale_x_log10() +
  xlab("Total Abundance") + 
  ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) +
  theme(legend.position = "none")

# Save plot to a file
ggsave(file.path(output_plot_path, "prevalence_plot.png"), plot, width = 10, height = 6, units = "in", dpi = 300)

# Define prevalence threshold as 5% of total samples
prevalenceThreshold <- 0.05 * nsamples(ps)

# Execute prevalence filter, using `prune_taxa()` function
keepTaxa <- rownames(prevdf1)[(prevdf1$Prevalence >= prevalenceThreshold)]
ps1 <- prune_taxa(keepTaxa, ps)
saveRDS(ps1, file.path(output_data_path, "ps1.Rds"))
# Merge everything to the phylum level
ps1_phylum <- tax_glom(ps1, "Phylum", NArm = TRUE)

# Transform Taxa counts to relative abundance
ps1_phylum_relabun <- transform_sample_counts(ps1_phylum, function(OTU) OTU/sum(OTU) * 100)
# Extract the data from the phyloseq object
taxa_abundance_table_phylum <- psmelt(ps1_phylum_relabun)

# box plot
BoxPlot_phylum <- taxa_abundance_table_phylum %>% 
  ggplot(aes(x = Phylum, y = Abundance, fill = Phylum)) +
  geom_boxplot() +
  labs(x = "",
       y = "Relative Abundance",
       title = "Phylum Relative Abundance") +
  facet_grid(~ Group, scales = "free") +
  theme(
    axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 10),
    strip.text = element_text(size = 12)
  )
ggsave(file.path(output_plot_path, "boxplot_phylum.png"), BoxPlot_phylum, width = 10, height = 6, units = "in", dpi = 300)


# Phylum Relative Abundance
ps1_phylum <- tax_glom(ps1, "Phylum", NArm = TRUE)  #!!!!!
# Get top 10 genera
top10_genera <- names(sort(taxa_sums(ps1_phylum), decreasing=TRUE))[1:10]
# Transform Taxa counts to relative abundance
ps1_phylum_relabun <- transform_sample_counts(ps1_phylum, function(OTU) OTU/sum(OTU) * 100)
# Extract the top 10 taxa and Regular Diet Samples
ps1_phylum_top10 <- prune_taxa(top10_genera, ps1_phylum_relabun)
# Convert into dataframe
taxa_abundance_table_phylum <- psmelt(ps1_phylum_top10)

stacked_barplot_phylum <- taxa_abundance_table_phylum %>% 
  ggplot(aes(x = Sample, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity") +
  labs(x = "",
       y = "Relative Abundance",
       title = "Phylum Relative Abundance") +
  facet_grid(~ Group, scales = "free") +
  theme(
    axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 10),
    strip.text = element_text(size = 12)
  )
ggsave(file.path(output_plot_path, "stacked_barplot_phylum.png"), stacked_barplot_phylum, width = 10, height = 6, units = "in", dpi = 300)


# vegan for rarifaction curve 
seqtab.nochim <- readRDS(file.path(output_data_path, "seqtab.nochim.Rds"))
otutab <- t(seqtab.nochim)
otu_matrix <- as.matrix(apply(otutab, 2, as.integer))
rarecurve(otu_matrix, step=50, col="blue", cex=0.6, xlab="Sampling Effort", ylab="OTU Richness")

png(file.path(output_plot_path, "rarefaction_curve.png"), width=700, height=500)
rarecurve(otu_matrix, step=50, col="blue", cex=0.6, xlab="Sampling Size", ylab="OTU Richness")
dev.off()

# vegan for MDS analysis
seqtab.nochim <- readRDS(file.path(output_data_path, "seqtab.nochim.Rds"))
otutab <- seqtab.nochim
dist_matrix <- vegdist(otutab, method="jaccard")
mds_result <- cmdscale(dist_matrix, eig=TRUE, k=3)

set.seed(123)
mds_df <- as.data.frame(mds_result$points)
clusters <- kmeans(mds_df[, 1:2], centers = 3)
mds_df$Cluster <- as.factor(clusters$cluster)
ggplot(mds_df, aes(x = V1, y = V2, color = Cluster)) +
  geom_point(alpha = 0.6, size = 2) +
  ggtitle("MDS Plot with K-means Clustering") +
  xlab("MDS1") + ylab("MDS2") +
  scale_color_brewer(palette = "Set1")

ggsave(file.path(output_plot_path, "MDS_plot.png"), width=10, height=8)

total_variation = sum(mds_result$eig)
explained_variation = sum(mds_result$eig[1:6]) / total_variation * 100
print(paste("Explained Variation:", explained_variation, "%"))

# cumulative variation
eigenvalues <- mds_result$eig
total_variation <- sum(eigenvalues)
cumulative_variation <- cumsum(eigenvalues) / total_variation

cumulative_data <- data.frame(Dimension = 1:length(cumulative_variation), CumulativeVariation = cumulative_variation)

p <- ggplot(cumulative_data, aes(x = Dimension, y = CumulativeVariation)) +
  geom_line() + geom_point() +
  scale_y_continuous(labels = scales::percent_format()) +
  ggtitle("Cumulative Variation Explained by Each Dimension") +
  xlab("Dimension") +
  ylab("Cumulative Variation (%)") +
  theme(
    axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 10),
    strip.text = element_text(size = 12)
  )

ggsave(file.path(output_plot_path, "cumulative_variation.png"), plot = p, width = 10, height = 6, units = "in", dpi = 300)
