###LOAD LIBRARIES
library(phyloseq)
library(readxl)
library(ape)
library(microbiome)
library(ggplot2)
library(dplyr)
library(stringr)
library(tidyr)
library(scales)
library(reshape2)
library(ggpubr)
library(vegan)
library(dendextend)
library(tidyverse)
library(forcats)
library(pheatmap)
library(RColorBrewer)
library(readr)
library(tibble)
# Set working directory
#setwd("C:/Users/.....")

######################### PHYLOSEQ OBJECT CREATION #################

# 1. Import OTU table
# Load OTU table
otu <- read.table("feature-table.tsv",
                  header = TRUE,
                  sep = "\t",
                  row.names = 1,
                  skip = 1,               # Skip the "# Constructed from biom file" line
                  comment.char = "",     # Prevent "#" from removing header
                  check.names = FALSE)   # Keep sample names as-is (e.g. "sample1")

# Optional: Remove taxonomy column if it exists
if ("taxonomy" %in% colnames(otu)) {
  otu <- otu[, -which(colnames(otu) == "taxonomy")]
}

# Convert to matrix
otu_mat <- as.matrix(otu)

# Create phyloseq OTU table object
OTU <- otu_table(otu_mat, taxa_are_rows = TRUE)


tax <- read.table("exported-taxonomy/taxonomy.tsv", 
                  header = TRUE, 
                  sep = "\t", 
                  row.names = 1, 
                  comment.char = "", 
                  quote = "", 
                  stringsAsFactors = FALSE)

# Split taxonomy into levels
tax_split <- strsplit(tax$Taxon, ";\\s*")
max_ranks <- max(sapply(tax_split, length))
tax_matrix <- t(sapply(tax_split, function(x) c(x, rep(NA, max_ranks - length(x)))))

# Assign biological rank names
colnames(tax_matrix) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")[1:max_ranks]

# ❗ Assign original ASV hashes as row names
rownames(tax_matrix) <- rownames(tax)

# Create taxonomy table
TAX <- tax_table(tax_matrix)


# 3. Import tree
tree <- read_tree("exported-tree/tree.nwk")


# Step 1: Read metadata
metadata <- read_excel("metadata_milk.xlsx")

# Step 2: Convert to data frame
metadata <- as.data.frame(metadata)

# Step 3: Set row names to "SAMPLE ID"
rownames(metadata) <- metadata$`SAMPLE ID`

# Step 4: Remove the "SAMPLE ID" column now that it's in the rownames
metadata$`SAMPLE ID` <- NULL

# Step 5: Convert to phyloseq sample_data object
META <- sample_data(metadata)


# 7. Create phyloseq object
ps <- phyloseq(OTU, TAX, META, tree)

# Save the phyloseq object
saveRDS(ps, file = "phyloseq_object.rds")

ps <- readRDS("phyloseq_object.rds")

############ CORE MICROBIAL VISUALIZATION ##################
#####################PHYLUM#################################

# Normalize to relative abundance
physeq_rel <- transform_sample_counts(ps, function(x) x / sum(x))

# Agglomerate to Phylum level
physeq_phylum <- tax_glom(physeq_rel, taxrank = "Phylum")

# Convert to long format for ggplot
df <- psmelt(physeq_phylum)

# Convert to percentage
df$Abundance <- df$Abundance * 100

# Optionally: round for clean labels
df$label <- ifelse(df$Abundance > 1, paste0(round(df$Abundance, 1), "%"), "")

# Plot
phy <- ggplot(df, aes(x = Sample, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 3, color = "black") +
  facet_wrap(~ Product, scales = "free_x") +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.text = element_text(size = 13, face = "bold"),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12)
  ) +
  labs(
    title = "Phylum-Level Composition by Product (with Percentages)",
    x = "Samples (Grouped by Product)",
    y = "Relative Abundance (%)"
  )

ggsave("Phylum_Level_Composition.png", plot = phy, width = 10, height = 8, dpi = 800, bg = "white")

####################GENUS#######################
####################GENUS#######################

# Transform to relative abundance
ps_rel <- transform(ps, "compositional")

# Agglomerate at Genus level
ps_genus <- tax_glom(ps_rel, taxrank = "Genus")

# Melt to long format
df <- psmelt(ps_genus)

# Add Product info from sample_data
df$Product <- sample_data(ps_genus)$Product[match(df$Sample, rownames(sample_data(ps_genus)))]

# Identify top 10 genera by overall abundance
top10_genera <- df %>%
  group_by(Genus) %>%
  summarise(TotalAbundance = sum(Abundance)) %>%
  arrange(desc(TotalAbundance)) %>%
  slice_head(n = 10) %>%
  pull(Genus)

# Identify top 5 for labeling
top5_genera <- df %>%
  group_by(Genus) %>%
  summarise(TotalAbundance = sum(Abundance)) %>%
  arrange(desc(TotalAbundance)) %>%
  slice_head(n = 5) %>%
  pull(Genus)

# Filter and group
df_top10 <- df %>%
  mutate(Genus = as.character(Genus),
         Genus = ifelse(Genus %in% top10_genera, Genus, "Other")) %>%
  group_by(Product, Genus) %>%
  summarise(Abundance = mean(Abundance), .groups = "drop") %>%
  mutate(Percentage = Abundance * 100,
         Label = ifelse(Genus %in% top5_genera & Percentage > 1,
                        paste0(round(Percentage, 1), "%"), ""))

# Optional: wrap long product names
df_top10$Product <- str_wrap(df_top10$Product, width = 10)

# Plot with percentage labels for top 5 genera
gn <- ggplot(df_top10, aes(x = Product, y = Percentage, fill = Genus)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(label = Label), position = position_stack(vjust = 0.5),
            size = 3, color = "black") +
  theme_minimal() +
  labs(title = "Top 10 Genera by Product (with % for Top 5)",
       x = "Product", y = "Mean Relative Abundance (%)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "bold")) +
  scale_fill_brewer(palette = "Set3")

ggsave("Genus_Level_Composition.png", plot = gn, width = 10, height = 8, dpi = 800, bg = "white")



#################SPECIES PLOT WITHOUT USING PHYLOSEQ OBJECT#############
# Step 1: Load Excel file
#setwd("C:/Users/.....")
data <- read_excel("taxa_freq_table7.xlsx")

# Step 2: Convert to long format
df_long <- data %>%
  pivot_longer(cols = -Species, names_to = "Product", values_to = "Abundance") %>%
  filter(!is.na(Abundance))

# Step 3: Calculate relative abundance (percentage) per product
df_long <- df_long %>%
  group_by(Product) %>%
    mutate(TotalProduct = sum(Abundance),
         Percentage = (Abundance / TotalProduct) * 100) %>%
  ungroup()

# Step 4: Identify top 10 species per product
top10_per_product <- df_long %>%
  group_by(Product) %>%
  slice_max(order_by = Abundance, n = 10, with_ties = FALSE) %>%
  mutate(IsTop10 = TRUE)

# Step 5: Identify top 3 species per product for labeling
top3_labels <- df_long %>%
  group_by(Product) %>%
  slice_max(order_by = Abundance, n = 3, with_ties = FALSE) %>%
  select(Species, Product)

# Step 6: Keep only top 10 per product (plus “Other”) and label top 3
df_labeled <- df_long %>%
  left_join(top10_per_product %>% select(Species, Product, IsTop10),
            by = c("Species", "Product")) %>%
  mutate(Species = ifelse(is.na(IsTop10), "Other", Species)) %>%
  select(-IsTop10) %>%
  group_by(Product, Species) %>%
  summarise(Percentage = sum(Percentage), .groups = "drop") %>%
  mutate(Label = ifelse(paste(Species, Product) %in%
                          paste(top3_labels$Species, top3_labels$Product),
                        paste0(round(Percentage, 1), "%"), "")) %>%
  filter(Species %in% top10_per_product$Species | Species == "Other")

# Step 7: Custom color palette (ensure enough colors for unique species)
custom_palette <- c(
  "#fb8072", "#80b1d3", "#fdb462", "#b3de69", "#fccde5",
  "#fdbf6f", "#cab2d6", "#ff7f00", "#6a3d9a", "#1b9e77",
  "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02",
  "#a6761d", "#666666", "#a6cee3", "#fb9a99", "#b15928",
  "#ffffb3", "#bebada"
)

# Step 8: Plot pie charts with facets and percentage labels
p <- ggplot(df_labeled, aes(x = "", y = Percentage, fill = Species)) +
  geom_bar(stat = "identity", width = 1) +
  geom_text(aes(label = Label), position = position_stack(vjust = 0.5), size = 5) +
  coord_polar(theta = "y") +
  facet_wrap(~ Product, scales = "free") +
  scale_fill_manual(values = custom_palette) +
  theme_minimal() +
  labs(x = NULL, y = "Relative Abundance (%)") +
  theme(
    strip.text = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

# Step 9: Save high-resolution image
ggsave("top_species_pie_facets.png", plot = p, width = 12, height = 8, dpi = 800, bg = "white")






################### ALPHA DIVERSITY PHYLUM LEVEL#############################


ps <- readRDS("phyloseq_object.rds")

# ----------------------------------------------------
# 1. Agglomerate to Phylum level
# ----------------------------------------------------
ps_phylum <- tax_glom(ps, taxrank = "Phylum")

# ----------------------------------------------------
# 2. Estimate alpha diversity metrics
# ----------------------------------------------------
alpha_div <- estimate_richness(ps_phylum, measures = c("Observed", "Shannon", "Chao1"))
alpha_div$SampleID <- rownames(alpha_div)

# ----------------------------------------------------
# 3. Merge with Product info
# ----------------------------------------------------
meta <- data.frame(sample_data(ps_phylum))
meta$SampleID <- rownames(meta)
alpha_merged <- left_join(alpha_div, meta, by = "SampleID")

# ----------------------------------------------------
# 4. Save alpha diversity values to CSV
# ----------------------------------------------------
write.csv(alpha_merged, "phylum_level_alpha_diversity_values.csv", row.names = FALSE)
write.csv(alpha_merged, "phylum_level_alpha_diversity_values.csv", row.names = FALSE)
# ----------------------------------------------------
# 5. Perform Kruskal-Wallis tests
# ----------------------------------------------------
kw_observed <- kruskal.test(Observed ~ Product, data = alpha_merged)
kw_shannon  <- kruskal.test(Shannon ~ Product, data = alpha_merged)
kw_chao1    <- kruskal.test(Chao1 ~ Product, data = alpha_merged)

# Save test results
kw_results <- data.frame(
  Metric = c("Observed", "Shannon", "Chao1"),
  p_value = c(kw_observed$p.value, kw_shannon$p.value, kw_chao1$p.value)
)

write.csv(kw_results, "kruskal_test_results_phylum_level.csv", row.names = FALSE)

# ----------------------------------------------------
# 6. Plot with p-values
# ----------------------------------------------------
alpha_long <- melt(alpha_merged, id.vars = c("SampleID", "Product"),
                   measure.vars = c("Observed", "Shannon", "Chao1"))

colors <- RColorBrewer::brewer.pal(n = 5, name = "Set2")
names(colors) <- unique(alpha_long$Product)

ggplot(alpha_long, aes(x = Product, y = value, fill = Product)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(aes(color = Product), width = 0.2, alpha = 0.8, size = 4) +  # Larger dots
  facet_wrap(~variable, scales = "free_y") +
  stat_compare_means(method = "kruskal.test", label = "p.format", label.y.npc = "top") +
  theme_minimal() +
  labs(title = "Alpha Diversity by Product (Phylum Level) with Kruskal-Wallis Tests",
       x = "Product", y = "Diversity Index") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors)



################ALPHA DIVERSITY GENUS LEVEL#################################
############################################################################



# ----------------------------------------------------
# 1. Agglomerate to Genus level
# ----------------------------------------------------
ps_genus <- tax_glom(ps, taxrank = "Genus")

# ----------------------------------------------------
# 2. Estimate alpha diversity metrics
# ----------------------------------------------------
alpha_div <- estimate_richness(ps_genus, measures = c("Observed", "Shannon", "Chao1"))
alpha_div$SampleID <- rownames(alpha_div)

# ----------------------------------------------------
# 3. Merge with Product info
# ----------------------------------------------------
meta <- data.frame(sample_data(ps_genus))
meta$SampleID <- rownames(meta)
alpha_merged <- left_join(alpha_div, meta, by = "SampleID")

# ----------------------------------------------------
# 4. Save alpha diversity values to CSV
# ----------------------------------------------------
write.csv(alpha_merged, "genus_level_alpha_diversity_values.csv", row.names = FALSE)

# ----------------------------------------------------
# 5. Perform Kruskal-Wallis tests
# ----------------------------------------------------
kw_observed <- kruskal.test(Observed ~ Product, data = alpha_merged)
kw_shannon  <- kruskal.test(Shannon ~ Product, data = alpha_merged)
kw_chao1    <- kruskal.test(Chao1 ~ Product, data = alpha_merged)

# Save test results
kw_results <- data.frame(
  Metric = c("Observed", "Shannon", "Chao1"),
  p_value = c(kw_observed$p.value, kw_shannon$p.value, kw_chao1$p.value)
)

write.csv(kw_results, "kruskal_test_results_genus_level.csv", row.names = FALSE)

# ----------------------------------------------------
# 6. Plot with p-values
# ----------------------------------------------------
alpha_long <- melt(alpha_merged, id.vars = c("SampleID", "Product"),
                   measure.vars = c("Observed", "Shannon", "Chao1"))

colors <- RColorBrewer::brewer.pal(n = 5, name = "Set2")  # or any palette you prefer
names(colors) <- unique(alpha_long$Product)

ggplot(alpha_long, aes(x = Product, y = value, fill = Product)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(aes(color = Product), width = 0.2, alpha = 0.8, size = 4) +  # Increased dot size
  facet_wrap(~variable, scales = "free_y") +
  stat_compare_means(method = "kruskal.test", label = "p.format", label.y.npc = "top") +
  theme_minimal() +
  labs(title = "Alpha Diversity by Product (Genus Level) with Kruskal-Wallis Tests",
       x = "Product", y = "Diversity Index") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors)







#################DENDROGRAMS NO COLOUR##########################
#######scritp to extract unifract distnace######################


# Clean object
ps_clean <- prune_samples(sample_sums(ps) > 0, ps)
ps_clean <- prune_taxa(taxa_sums(ps_clean) > 0, ps_clean)

# Compute Weighted UniFrac distance
unifrac_dist <- distance(ps_clean, method = "unifrac", weighted = TRUE)

# ✅ Save UniFrac distance matrix as CSV
unifrac_matrix <- as.matrix(unifrac_dist)
write.csv(unifrac_matrix, file = "weighted_unifrac_distance_matrix.csv")
cat("✅ UniFrac distance matrix saved as 'weighted_unifrac_distance_matrix.csv'\n")

# Perform hierarchical clustering
hc <- hclust(as.dist(unifrac_dist), method = "average")

# Extract metadata
meta <- data.frame(sample_data(ps_clean))
product_labels <- as.character(meta$Product)
names(product_labels) <- rownames(meta)

# Convert hc to dendrogram
dend <- as.dendrogram(hc)

# Match product labels to dendrogram tip order
labels_in_dend <- labels(dend)
tip_labels <- product_labels[labels_in_dend]

# Create consistent color map
unique_products <- unique(tip_labels)
product_colors <- setNames(RColorBrewer::brewer.pal(length(unique_products), "Set2"), unique_products)

# Apply product labels and colors to dendrogram
dend <- dend %>%
  set("labels", tip_labels) %>%
  set("labels_colors", product_colors[tip_labels])  # exact color mapping

# Save plot
png("products_dendrogram_unifrac_named.png", width = 1000, height = 800)
par(cex = 1.8)

plot(dend,
     main = "Dendrogram (Weighted UniFrac) Labeled by Product",
     ylab = "Weighted UniFrac Distance")

legend("topright",
       legend = names(product_colors),
       col = product_colors,
       pch = 15, cex = 1.2)

dev.off()

cat("✅ Plot saved as 'products_dendrogram_unifrac_named.png'\n")


############################################################################
################################bray-curtis analysis########################

# Clean phyloseq object: remove empty samples and taxa
ps_clean <- prune_samples(sample_sums(ps) > 0, ps)
ps_clean <- prune_taxa(taxa_sums(ps_clean) > 0, ps_clean)

# Compute Bray-Curtis distance matrix
dist_all <- distance(ps_clean, method = "bray")

bray_matrix <- as.matrix(dist_all)
write.csv(bray_matrix, file = "bray_curtis_distance_matrix.csv")
cat("✅ Bray distance matrix saved as 'bray_distance_matrix.csv'\n")


# Perform hierarchical clustering
hc <- hclust(as.dist(dist_all), method = "average")

# Extract sample metadata
sample_info <- data.frame(sample_data(ps_clean))
product_labels <- as.character(sample_info$Product)
names(product_labels) <- rownames(sample_info)

# Get labels in dendrogram order
labels_in_dend <- labels(as.dendrogram(hc))
tip_labels <- product_labels[labels_in_dend]

# Create consistent color mapping
unique_products <- unique(tip_labels)
product_colors <- setNames(RColorBrewer::brewer.pal(length(unique_products), "Set2"), unique_products)

# Apply labels and colors
dend <- as.dendrogram(hc)
dend <- dend %>%
  set("labels", tip_labels) %>%
  set("labels_colors", product_colors[tip_labels])

# Save dendrogram to PNG
png("products_dendrogram_bray_named.png", width = 1000, height = 800)
par(cex = 1.8)  # Increase font size of tip labels

plot(dend,
     main = "Dendrogram of Samples (Bray-Curtis) Labeled by Product",
     ylab = "Bray-Curtis Dissimilarity")

legend("topright",
       legend = names(product_colors),
       col = product_colors,
       pch = 15, cex = 1.2)  # Increase legend font size

dev.off()

cat("✅ Dendrogram saved as 'products_dendrogram_bray_named.png'\n")




##################FUNCTIONAL ANALYSIS###############################
#################set working directory##############################
# Read PICRUSt2 output files
ec_unstra <- read_excel("C:/Users/DAVID/Desktop/phyloseq/trial/EC_metagenome_out/pred_metagenome_unstrat.xlsx")
ko_unstra <- read_excel("C:/Users/DAVID/Desktop/phyloseq/trial/KO_metagenome_out/pred_metagenome_unstrat.xlsx")
pathways_unstra <- read_excel("C:/Users/DAVID/Desktop/phyloseq/trial/pathways_out/path_abun_unstrat.xlsx")

# Convert from wide to long format for ggplot2
# Explicitly rename first column to "function"
colnames(ec_unstra)[1] <- "function"
colnames(ko_unstra)[1] <- "function"
colnames(pathways_unstra)[1] <- "function"

ec_unstra_long <- ec_unstra %>%
  pivot_longer(cols = -`function`, names_to = "SampleID", values_to = "Abundance")

ko_unstra_long <- ko_unstra %>%
  pivot_longer(cols = -`function`, names_to = "SampleID", values_to = "Abundance")

pathways_unstra_long <- pathways_unstra %>%
  pivot_longer(cols = -`function`, names_to = "SampleID", values_to = "Abundance")


metadata <- data.frame(sample_data(ps)) %>%
  tibble::rownames_to_column("SampleID")
ko_unstra_joined <- left_join(ko_unstra_long, metadata, by = "SampleID")
ec_unstra_joined <- left_join(ec_unstra_long, metadata, by = "SampleID")
pathways_unstra_joined <- left_join(pathways_unstra_long, metadata, by = "SampleID")



################################################################################
###########PATHWAY COMPARATIVE ANALYSIS#########################################

# Get top 10 pathways
top_pathways_unstra <- pathways_unstra_long %>%
  group_by(`function`) %>%
  summarise(Total = sum(Abundance)) %>%
  top_n(10, Total) %>%
  pull(`function`)

# Define custom colors
custom_colors <- c(
  "Ghee" = "#E69F00",
  "Nono" = "#56B4E9",
  "Nunu" = "#009E73",
  "Wara" = "#F0E442",
  "Kwerionik" = "#D55E00"
)

# Prepare data
plot_data <- pathways_unstra_long %>%
  filter(`function` %in% top_pathways_unstra) %>%
  left_join(metadata, by = "SampleID") %>%
  group_by(Product, `function`) %>%
  summarise(mean_abundance = mean(Abundance), .groups = "drop") %>%
  mutate(Product = factor(Product, levels = names(custom_colors)))

# Create the plot
p <- ggplot(plot_data, aes(x = mean_abundance, y = fct_reorder(`function`, mean_abundance), fill = Product)) +
  geom_col(position = position_dodge(width = 0.9), width = 0.7) +
  scale_fill_manual(values = custom_colors) +
  theme_minimal(base_size = 15) +  # Slightly larger base size
  labs(
    title = "Top 10 Predicted Pathways by Product",
    x = "Mean Abundance",
    y = "Pathway",
    fill = "Product"
  ) +
  theme(
    axis.text.y = element_text(face = "bold", size = 12),
    axis.text.x = element_text(size = 12),
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 12),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

# Save as high-resolution PNG
ggsave("top10_pathways_barplot.png", plot = p, width = 10, height = 6, dpi = 400)

readr::write_csv(plot_data, paste0("pathway_plot_data.csv"))

# Optional: also display in R session
print(p)


#############heatmap_ko Comparative Analysis ####################
#############heatmap-ko Comparative Analysis ###################


# Define custom product colors
custom_colors <- c(
  "Ghee" = "#E69F00",
  "Nono" = "#56B4E9",
  "Nunu" = "#009E73",
  "Wara" = "#F0E442",
  "Kwerionik" = "#D55E00"
)

# Ensure Product is a factor with desired order
ko_unstra_joined <- ko_unstra_joined %>%
  mutate(Product = factor(Product, levels = names(custom_colors)))

# Get top 10 most abundant KOs
top_kos <- ko_unstra_joined %>%
  group_by(`function`) %>%
  summarise(Total = sum(Abundance), .groups = "drop") %>%
  top_n(10, Total) %>%
  pull(`function`)

# Filter for top 10 KOs
top_ko_data <- ko_unstra_joined %>%
  filter(`function` %in% top_kos)

# Pivot to KO × Sample matrix (mean across replicates)
heatmap_data <- top_ko_data %>%
  group_by(SampleID, `function`) %>%
  summarise(mean_abundance = mean(Abundance), .groups = "drop") %>%
  pivot_wider(names_from = SampleID, values_from = mean_abundance, values_fill = 0) %>%
  column_to_rownames(var = "function")

# Prepare annotation for columns (samples)
sample_metadata <- ko_unstra_joined %>%
  select(SampleID, Product) %>%
  distinct() %>%
  column_to_rownames("SampleID")

# Plot heatmap with clustering
pheatmap(
  mat = as.matrix(heatmap_data),
  annotation_col = sample_metadata,
  annotation_colors = list(Product = custom_colors),
  color = colorRampPalette(brewer.pal(9, "YlOrRd"))(100),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  fontsize = 11,
  
  main = "Heatmap of Top 10 Predicted KOs Across Samples"
)


# Step 1: Calculate total KO abundance per sample to compute relative abundance
relative_ko_data <- top_ko_data %>%
  group_by(SampleID) %>%
  mutate(TotalSampleAbundance = sum(Abundance)) %>%
  ungroup() %>%
  mutate(RelativeAbundance = Abundance / TotalSampleAbundance * 100)  # in percent

# Step 2: Calculate mean abundance and mean relative abundance per Product and KO
ko_summary <- relative_ko_data %>%
  group_by(Product, `function`) %>%
  summarise(
    mean_abundance = mean(Abundance),
    mean_relative_abundance = mean(RelativeAbundance),
    .groups = "drop"
  )

# Step 3: Save to CSV
write.csv(ko_summary, "mean_abundance_relative_abundance_top_kos.csv", row.names = FALSE)

# Optional: View part of the summary
head(ko_summary)



#############################################################################
############ EC EC HEATMAP Comparative Analysis##############################

# Define custom product colors
custom_colors <- c(
  "Ghee" = "#E69F00",
  "Nono" = "#56B4E9",
  "Nunu" = "#009E73",
  "Wara" = "#F0E442",
  "Kwerionik" = "#D55E00"
)

# Step 1: Join EC predictions with metadata
ec_unstra_joined <- left_join(ec_unstra_long, metadata, by = "SampleID") %>%
  mutate(Product = factor(Product, levels = names(custom_colors)))

# Step 2: Identify top 10 most abundant ECs
top_ecs <- ec_unstra_joined %>%
  group_by(`function`) %>%
  summarise(Total = sum(Abundance), .groups = "drop") %>%
  top_n(10, Total) %>%
  pull(`function`)

# Step 3: Filter for top 10 ECs
top_ec_data <- ec_unstra_joined %>%
  filter(`function` %in% top_ecs)

# Step 4: Prepare wide matrix (EC x SampleID)
heatmap_data <- top_ec_data %>%
  group_by(`function`, SampleID) %>%
  summarise(mean_abundance = mean(Abundance), .groups = "drop") %>%
  pivot_wider(names_from = SampleID, values_from = mean_abundance, values_fill = 0) %>%
  column_to_rownames("function")

# Step 5: Prepare annotation data for samples
sample_metadata <- ec_unstra_joined %>%
  select(SampleID, Product) %>%
  distinct() %>%
  column_to_rownames("SampleID")

# Step 6: Save top ECs summary (optional)
top_ec_table <- ec_unstra_joined %>%
  filter(`function` %in% top_ecs) %>%
  group_by(Product, `function`) %>%
  summarise(mean_abundance = mean(Abundance), .groups = "drop")

write_csv(top_ec_table, "top10_predicted_ECs_by_product.csv")

# Step 7: Plot heatmap with clustering and dendrogram
pheatmap(
  mat = as.matrix(heatmap_data),
  annotation_col = sample_metadata,
  annotation_colors = list(Product = custom_colors),
  scale = "row",  # normalize across rows (ECs)
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "ward.D2",
  color = colorRampPalette(c("white", "red"))(100),
  fontsize_row = 10,
  fontsize_col = 10,
  main = "Heatmap of Top 10 Predicted EC Numbers Across Samples"
)


################################################################################
######################### SAMPLE SPECIFIC FUNCTIONAL ANALYSIS###################
############### SAMPLE NAME = NONO ############################

# Load taxonomy
taxonomy <- read_tsv("C:/Users/DAVID/Desktop/phyloseq/trial/exported-taxonomy/taxonomy.tsv")

# Split taxonomic levels
taxonomy_clean <- taxonomy %>%
  separate(Taxon, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
           sep = ";\\s*", fill = "right", remove = FALSE) %>%
  rename(ASV = `Feature ID`)

ko_contrib <- read_tsv("C:/Users/DAVID/Desktop/phyloseq/trial/KO_metagenome_cont_out/pred_metagenome_contrib.tsv", 
                       show_col_types = FALSE)

ec_contrib <- read_tsv("C:/Users/DAVID/Desktop/phyloseq/trial/EC_metagenome_cont_out/pred_metagenome_contrib.tsv", 
                       show_col_types = FALSE)
pathway_contrib <- read_tsv("C:/Users/DAVID/Desktop/phyloseq/trial/pathways_cont_out/path_abun_contrib.tsv", 
                            show_col_types = FALSE)
# Join KO contributions to taxonomic classification
ko_taxa_contrib <- left_join(ko_contrib, taxonomy_clean, by = c("taxon" = "ASV"))
ec_taxa_contrib <- left_join(ec_contrib, taxonomy_clean, by = c("taxon" = "ASV"))
pathway_taxa_contrib <- left_join(pathway_contrib, taxonomy_clean, by = c("taxon" = "ASV"))

# === STEP 0: SELECT SAMPLE ===
target_sample <- "sample1"  # change to your actual sample name
target_data <- ko_taxa_contrib %>% filter(sample == target_sample)

# === STEP 1: Top 10 KOs in selected sample ===
top_ko <- target_data %>%
  group_by(`function`) %>%
  summarise(total = sum(taxon_function_abun), .groups = "drop") %>%
  slice_max(total, n = 10)

# === STEP 2: Filter to top KOs ===
filtered_data <- target_data %>%
  filter(`function` %in% top_ko$`function`, !is.na(Genus))

# === STEP 3: Top 10 Genera in that sample and those KOs ===
top_genus <- filtered_data %>%
  group_by(Genus) %>%
  summarise(total = sum(taxon_function_abun), .groups = "drop") %>%
  slice_max(total, n = 10)

# === STEP 4: Filter again ===
filtered_data <- filtered_data %>%
  filter(Genus %in% top_genus$Genus)

# === STEP 5: Build function × genus matrix ===
heatmap_df <- filtered_data %>%
  group_by(`function`, Genus) %>%
  summarise(total = sum(taxon_function_abun), .groups = "drop") %>%
  pivot_wider(names_from = Genus, values_from = total, values_fill = 0)

heatmap_matrix <- heatmap_df %>%
  column_to_rownames("function") %>%
  as.matrix()

# === STEP 6: Plot heatmap ===
pheatmap(heatmap_matrix,
         main = paste("Top KO Functions in Nono"),
         color = colorRampPalette(brewer.pal(9, "YlGnBu"))(100),
         scale = "row",
         fontsize_row = 10,
         fontsize_col = 10,
         border_color = NA,
         cluster_rows = TRUE,
         cluster_cols = TRUE)

# === OPTIONAL: Save CSV (counts + % per function) ===
percent_mat <- sweep(heatmap_matrix, 1, rowSums(heatmap_matrix), FUN = "/") * 100
count_df <- as.data.frame(heatmap_matrix)
percent_df <- as.data.frame(percent_mat)
colnames(count_df) <- paste0(colnames(count_df), "_count")
colnames(percent_df) <- paste0(colnames(percent_df), "_percent")
count_df$KO <- rownames(count_df)
percent_df$KO <- rownames(percent_df)
combined <- left_join(count_df, percent_df, by = "KO")
write.csv(combined, paste0("KO_Function_vs_Genus_sample1", target_sample, ".csv"), row.names = FALSE)



############### SAMPLE NAME = WARA ############################

# === STEP 0: SELECT SAMPLE 2 ===
target_sample_2 <- "sample2"  # replace with your actual sample ID
target_data_2 <- ko_taxa_contrib %>% filter(sample == target_sample_2)

# === STEP 1: Top 10 KOs in Sample 2 ===
top_ko_2 <- target_data_2 %>%
  group_by(`function`) %>%
  summarise(total = sum(taxon_function_abun), .groups = "drop") %>%
  slice_max(total, n = 10)

# === STEP 2: Filter to top KOs ===
filtered_data_2 <- target_data_2 %>%
  filter(`function` %in% top_ko_2$`function`, !is.na(Genus))

# === STEP 3: Top 10 Genera contributing to top KOs ===
top_genus_2 <- filtered_data_2 %>%
  group_by(Genus) %>%
  summarise(total = sum(taxon_function_abun), .groups = "drop") %>%
  slice_max(total, n = 10)

# === STEP 4: Filter again for top genera ===
filtered_data_2 <- filtered_data_2 %>%
  filter(Genus %in% top_genus_2$Genus)

# === STEP 5: Create KO × Genus matrix ===
heatmap_df_2 <- filtered_data_2 %>%
  group_by(`function`, Genus) %>%
  summarise(total = sum(taxon_function_abun), .groups = "drop") %>%
  pivot_wider(names_from = Genus, values_from = total, values_fill = 0)

heatmap_matrix_2 <- heatmap_df_2 %>%
  column_to_rownames("function") %>%
  as.matrix()

# === STEP 6: Plot Heatmap for Sample 2 ===
pheatmap(heatmap_matrix_2,
         main = paste("Top KO Functions in Wara"),
         color = colorRampPalette(c("black", "yellow", "red"))(100),
         scale = "row",
         fontsize_row = 10,
         fontsize_col = 10,
         border_color = NA,
         cluster_rows = TRUE,
         cluster_cols = TRUE)

# === OPTIONAL: Save CSV for Sample 2 (Counts + % per Function) ===
percent_mat_2 <- sweep(heatmap_matrix_2, 1, rowSums(heatmap_matrix_2), FUN = "/") * 100
count_df_2 <- as.data.frame(heatmap_matrix_2)
percent_df_2 <- as.data.frame(percent_mat_2)
colnames(count_df_2) <- paste0(colnames(count_df_2), "_count")
colnames(percent_df_2) <- paste0(colnames(percent_df_2), "_percent")
count_df_2$KO <- rownames(count_df_2)
percent_df_2$KO <- rownames(percent_df_2)
combined_2 <- left_join(count_df_2, percent_df_2, by = "KO")
write.csv(combined_2, paste0("KO_Function_vs_Genus_", target_sample_2, ".csv"), row.names = FALSE)



############### SAMPLE NAME = KWERIONIK ############################

# === STEP 0: SELECT SAMPLE 3 ===
target_sample_3 <- "sample3"  # replace with actual sample ID
target_data_3 <- ko_taxa_contrib %>% filter(sample == target_sample_3)

# === STEP 1: Top 10 KOs in Sample 3 ===
top_ko_3 <- target_data_3 %>%
  group_by(`function`) %>%
  summarise(total = sum(taxon_function_abun), .groups = "drop") %>%
  slice_max(total, n = 10)

# === STEP 2: Filter to top KOs ===
filtered_data_3 <- target_data_3 %>%
  filter(`function` %in% top_ko_3$`function`, !is.na(Genus))

# === STEP 3: Top 10 Genera contributing to top KOs ===
top_genus_3 <- filtered_data_3 %>%
  group_by(Genus) %>%
  summarise(total = sum(taxon_function_abun), .groups = "drop") %>%
  slice_max(total, n = 10)

# === STEP 4: Filter again for top genera ===
filtered_data_3 <- filtered_data_3 %>%
  filter(Genus %in% top_genus_3$Genus)

# === STEP 5: Create KO × Genus matrix ===
heatmap_df_3 <- filtered_data_3 %>%
  group_by(`function`, Genus) %>%
  summarise(total = sum(taxon_function_abun), .groups = "drop") %>%
  pivot_wider(names_from = Genus, values_from = total, values_fill = 0)

heatmap_matrix_3 <- heatmap_df_3 %>%
  column_to_rownames("function") %>%
  as.matrix()

# === STEP 6: Plot Heatmap for Sample 3 ===
pheatmap(heatmap_matrix_3,
         main = paste("Top KO Functions in Kwerionik"),
         color = colorRampPalette(c("pink", "white", "red"))(100),
         scale = "row",
         fontsize_row = 10,
         fontsize_col = 10,
         border_color = NA,
         cluster_rows = TRUE,
         cluster_cols = TRUE)

# === STEP 7: Save Combined CSV (Counts + Percentages) ===
percent_mat_3 <- sweep(heatmap_matrix_3, 1, rowSums(heatmap_matrix_3), FUN = "/") * 100
count_df_3 <- as.data.frame(heatmap_matrix_3)
percent_df_3 <- as.data.frame(percent_mat_3)
colnames(count_df_3) <- paste0(colnames(count_df_3), "_count")
colnames(percent_df_3) <- paste0(colnames(percent_df_3), "_percent")
count_df_3$KO <- rownames(count_df_3)
percent_df_3$KO <- rownames(percent_df_3)
combined_3 <- left_join(count_df_3, percent_df_3, by = "KO")
write.csv(combined_3, paste0("KO_Function_vs_Genus_", target_sample_3, ".csv"), row.names = FALSE)





############### SAMPLE NAME = GHEE ############################


# === STEP 0: SELECT SAMPLE 4 ===
target_sample_4 <- "sample4"  # replace with actual sample ID
target_data_4 <- ko_taxa_contrib %>% filter(sample == target_sample_4)

# === STEP 1: Top 10 KOs in Sample 4 ===
top_ko_4 <- target_data_4 %>%
  group_by(`function`) %>%
  summarise(total = sum(taxon_function_abun), .groups = "drop") %>%
  slice_max(total, n = 10)

# === STEP 2: Filter to top KOs ===
filtered_data_4 <- target_data_4 %>%
  filter(`function` %in% top_ko_4$`function`, !is.na(Genus))

# === STEP 3: Top 10 Genera contributing to those KOs ===
top_genus_4 <- filtered_data_4 %>%
  group_by(Genus) %>%
  summarise(total = sum(taxon_function_abun), .groups = "drop") %>%
  slice_max(total, n = 10)

# === STEP 4: Filter again to only top genera ===
filtered_data_4 <- filtered_data_4 %>%
  filter(Genus %in% top_genus_4$Genus)

# === STEP 5: Create KO × Genus matrix ===
heatmap_df_4 <- filtered_data_4 %>%
  group_by(`function`, Genus) %>%
  summarise(total = sum(taxon_function_abun), .groups = "drop") %>%
  pivot_wider(names_from = Genus, values_from = total, values_fill = 0)

heatmap_matrix_4 <- heatmap_df_4 %>%
  column_to_rownames("function") %>%
  as.matrix()

# === STEP 6: Plot Heatmap ===
pheatmap(heatmap_matrix_4,
         main = paste("Top KO Functions in Ghee"),
         color = colorRampPalette(c("black", "white", "red"))(100),
         scale = "row",
         fontsize_row = 10,
         fontsize_col = 10,
         border_color = NA,
         cluster_rows = TRUE,
         cluster_cols = TRUE)

# === STEP 7: Save Combined CSV (Counts + Percentages) ===
percent_mat_4 <- sweep(heatmap_matrix_4, 1, rowSums(heatmap_matrix_4), FUN = "/") * 100
count_df_4 <- as.data.frame(heatmap_matrix_4)
percent_df_4 <- as.data.frame(percent_mat_4)
colnames(count_df_4) <- paste0(colnames(count_df_4), "_count")
colnames(percent_df_4) <- paste0(colnames(percent_df_4), "_percent")
count_df_4$KO <- rownames(count_df_4)
percent_df_4$KO <- rownames(percent_df_4)
combined_4 <- left_join(count_df_4, percent_df_4, by = "KO")
write.csv(combined_4, paste0("KO_Function_vs_Genus_", target_sample_4, ".csv"), row.names = FALSE)



############### SAMPLE NAME = NUNU ############################

# === STEP 0: SELECT SAMPLE 5 ===
target_sample_5 <- "sample5"  # replace with your actual sample ID
target_data_5 <- ko_taxa_contrib %>% filter(sample == target_sample_5)

# === STEP 1: Top 10 KOs in Sample 5 ===
top_ko_5 <- target_data_5 %>%
  group_by(`function`) %>%
  summarise(total = sum(taxon_function_abun), .groups = "drop") %>%
  slice_max(total, n = 10)

# === STEP 2: Filter to top KOs ===
filtered_data_5 <- target_data_5 %>%
  filter(`function` %in% top_ko_5$`function`, !is.na(Genus))

# === STEP 3: Top 10 Genera contributing to those KOs ===
top_genus_5 <- filtered_data_5 %>%
  group_by(Genus) %>%
  summarise(total = sum(taxon_function_abun), .groups = "drop") %>%
  slice_max(total, n = 10)

# === STEP 4: Filter again to only top genera ===
filtered_data_5 <- filtered_data_5 %>%
  filter(Genus %in% top_genus_5$Genus)

# === STEP 5: Create KO × Genus matrix ===
heatmap_df_5 <- filtered_data_5 %>%
  group_by(`function`, Genus) %>%
  summarise(total = sum(taxon_function_abun), .groups = "drop") %>%
  pivot_wider(names_from = Genus, values_from = total, values_fill = 0)

heatmap_matrix_5 <- heatmap_df_5 %>%
  column_to_rownames("function") %>%
  as.matrix()

# === STEP 6: Plot Heatmap ===
pheatmap(heatmap_matrix_5,
         main = paste("Top KO Functions in Nunu"),
         color = colorRampPalette(c("blue", "white", "red"))(100),
         scale = "row",
         fontsize_row = 10,
         fontsize_col = 10,
         border_color = NA,
         cluster_rows = TRUE,
         cluster_cols = TRUE)

# === STEP 7: Save Combined CSV (Counts + Percentages) ===
percent_mat_5 <- sweep(heatmap_matrix_5, 1, rowSums(heatmap_matrix_5), FUN = "/") * 100
count_df_5 <- as.data.frame(heatmap_matrix_5)
percent_df_5 <- as.data.frame(percent_mat_5)
colnames(count_df_5) <- paste0(colnames(count_df_5), "_count")
colnames(percent_df_5) <- paste0(colnames(percent_df_5), "_percent")
count_df_5$KO <- rownames(count_df_5)
percent_df_5$KO <- rownames(percent_df_5)
combined_5 <- left_join(count_df_5, percent_df_5, by = "KO")
write.csv(combined_5, paste0("KO_Function_vs_Genus_", target_sample_5, ".csv"), row.names = FALSE)


############################################################################
################## EC EC EC SAMPEL SPECIFIC ################################
############### SAMPLE NAME = NONO ############################

# === STEP 0: SELECT SAMPLE 1 ===
target_sample_ec1 <- "sample1"  # Replace with your actual sample ID
target_data_ec1 <- ec_taxa_contrib %>% filter(sample == target_sample_ec1)

# === STEP 1: Top 10 EC functions in Sample 1 ===
top_ec_1 <- target_data_ec1 %>%
  group_by(`function`) %>%
  summarise(total = sum(taxon_function_abun), .groups = "drop") %>%
  slice_max(total, n = 10)

# === STEP 2: Filter to top ECs ===
filtered_data_ec1 <- target_data_ec1 %>%
  filter(`function` %in% top_ec_1$`function`, !is.na(Genus))

# === STEP 3: Top 10 Genera contributing to top ECs ===
top_genus_ec1 <- filtered_data_ec1 %>%
  group_by(Genus) %>%
  summarise(total = sum(taxon_function_abun), .groups = "drop") %>%
  slice_max(total, n = 10)

# === STEP 4: Filter again for top genera ===
filtered_data_ec1 <- filtered_data_ec1 %>%
  filter(Genus %in% top_genus_ec1$Genus)

# === STEP 5: Create EC × Genus matrix ===
heatmap_df_ec1 <- filtered_data_ec1 %>%
  group_by(`function`, Genus) %>%
  summarise(total = sum(taxon_function_abun), .groups = "drop") %>%
  pivot_wider(names_from = Genus, values_from = total, values_fill = 0)

heatmap_matrix_ec1 <- heatmap_df_ec1 %>%
  column_to_rownames("function") %>%
  as.matrix()

# === STEP 6: Plot EC Heatmap ===
pheatmap(heatmap_matrix_ec1,
         main = paste("Top EC Functions in Nono"),
         color = colorRampPalette(brewer.pal(9, "YlOrRd"))(100),
         scale = "row",
         fontsize_row = 10,
         fontsize_col = 10,
         border_color = NA,
         cluster_rows = TRUE,
         cluster_cols = TRUE)

# === STEP 7: Save Combined CSV (Counts + Percentages) ===
percent_mat_ec1 <- sweep(heatmap_matrix_ec1, 1, rowSums(heatmap_matrix_ec1), FUN = "/") * 100
count_df_ec1 <- as.data.frame(heatmap_matrix_ec1)
percent_df_ec1 <- as.data.frame(percent_mat_ec1)
colnames(count_df_ec1) <- paste0(colnames(count_df_ec1), "_count")
colnames(percent_df_ec1) <- paste0(colnames(percent_df_ec1), "_percent")
count_df_ec1$EC <- rownames(count_df_ec1)
percent_df_ec1$EC <- rownames(percent_df_ec1)
combined_ec1 <- left_join(count_df_ec1, percent_df_ec1, by = "EC")
write.csv(combined_ec1, paste0("EC_Function_vs_Genus_", target_sample_ec1, ".csv"), row.names = FALSE)



############### SAMPLE NAME = WARA ############################
# === STEP 0: SELECT SAMPLE 2 ===
target_sample_ec2 <- "sample2"  # Replace with your actual sample name
target_data_ec2 <- ec_taxa_contrib %>% filter(sample == target_sample_ec2)

# === STEP 1: Top 10 EC functions in Sample 2 ===
top_ec_2 <- target_data_ec2 %>%
  group_by(`function`) %>%
  summarise(total = sum(taxon_function_abun), .groups = "drop") %>%
  slice_max(total, n = 10)

# === STEP 2: Filter to top ECs ===
filtered_data_ec2 <- target_data_ec2 %>%
  filter(`function` %in% top_ec_2$`function`, !is.na(Genus))

# === STEP 3: Top 10 Genera contributing to top ECs ===
top_genus_ec2 <- filtered_data_ec2 %>%
  group_by(Genus) %>%
  summarise(total = sum(taxon_function_abun), .groups = "drop") %>%
  slice_max(total, n = 10)

# === STEP 4: Filter again for top genera ===
filtered_data_ec2 <- filtered_data_ec2 %>%
  filter(Genus %in% top_genus_ec2$Genus)

# === STEP 5: Create EC × Genus matrix ===
heatmap_df_ec2 <- filtered_data_ec2 %>%
  group_by(`function`, Genus) %>%
  summarise(total = sum(taxon_function_abun), .groups = "drop") %>%
  pivot_wider(names_from = Genus, values_from = total, values_fill = 0)

heatmap_matrix_ec2 <- heatmap_df_ec2 %>%
  column_to_rownames("function") %>%
  as.matrix()

# === STEP 6: Plot EC Heatmap ===
pheatmap(heatmap_matrix_ec2,
         main = paste("Top EC Functions in Wara"),
         color = colorRampPalette(brewer.pal(9, "YlGnBu"))(100),
         scale = "row",
         fontsize_row = 10,
         fontsize_col = 10,
         border_color = NA,
         cluster_rows = TRUE,
         cluster_cols = TRUE)

# === STEP 7: Save Combined CSV (Counts + Percentages) ===
percent_mat_ec2 <- sweep(heatmap_matrix_ec2, 1, rowSums(heatmap_matrix_ec2), FUN = "/") * 100
count_df_ec2 <- as.data.frame(heatmap_matrix_ec2)
percent_df_ec2 <- as.data.frame(percent_mat_ec2)
colnames(count_df_ec2) <- paste0(colnames(count_df_ec2), "_count")
colnames(percent_df_ec2) <- paste0(colnames(percent_df_ec2), "_percent")
count_df_ec2$EC <- rownames(count_df_ec2)
percent_df_ec2$EC <- rownames(percent_df_ec2)
combined_ec2 <- left_join(count_df_ec2, percent_df_ec2, by = "EC")
write.csv(combined_ec2, paste0("EC_Function_vs_Genus_", target_sample_ec2, ".csv"), row.names = FALSE)



############### SAMPLE NAME = KWERIONIK ############################

# === STEP 0: SELECT SAMPLE 3 ===
target_sample_ec3 <- "sample3"  # Replace with your actual sample ID
target_data_ec3 <- ec_taxa_contrib %>% filter(sample == target_sample_ec3)

# === STEP 1: Top 10 EC functions in Sample 3 ===
top_ec_3 <- target_data_ec3 %>%
  group_by(`function`) %>%
  summarise(total = sum(taxon_function_abun), .groups = "drop") %>%
  slice_max(total, n = 10)

# === STEP 2: Filter to top ECs ===
filtered_data_ec3 <- target_data_ec3 %>%
  filter(`function` %in% top_ec_3$`function`, !is.na(Genus))

# === STEP 3: Top 10 Genera contributing to top ECs ===
top_genus_ec3 <- filtered_data_ec3 %>%
  group_by(Genus) %>%
  summarise(total = sum(taxon_function_abun), .groups = "drop") %>%
  slice_max(total, n = 10)

# === STEP 4: Filter again for top genera ===
filtered_data_ec3 <- filtered_data_ec3 %>%
  filter(Genus %in% top_genus_ec3$Genus)

# === STEP 5: Create EC × Genus matrix ===
heatmap_df_ec3 <- filtered_data_ec3 %>%
  group_by(`function`, Genus) %>%
  summarise(total = sum(taxon_function_abun), .groups = "drop") %>%
  pivot_wider(names_from = Genus, values_from = total, values_fill = 0)

heatmap_matrix_ec3 <- heatmap_df_ec3 %>%
  column_to_rownames("function") %>%
  as.matrix()

# === STEP 6: Plot EC Heatmap ===
pheatmap(heatmap_matrix_ec3,
         main = paste("Top EC Functions in Kwerionik"),
         color = colorRampPalette(c("navy", "yellow", "red"))(100),
         scale = "row",
         fontsize_row = 10,
         fontsize_col = 10,
         border_color = NA,
         cluster_rows = TRUE,
         cluster_cols = TRUE)

# === STEP 7: Save Combined CSV (Counts + Percentages) ===
percent_mat_ec3 <- sweep(heatmap_matrix_ec3, 1, rowSums(heatmap_matrix_ec3), FUN = "/") * 100
count_df_ec3 <- as.data.frame(heatmap_matrix_ec3)
percent_df_ec3 <- as.data.frame(percent_mat_ec3)
colnames(count_df_ec3) <- paste0(colnames(count_df_ec3), "_count")
colnames(percent_df_ec3) <- paste0(colnames(percent_df_ec3), "_percent")
count_df_ec3$EC <- rownames(count_df_ec3)
percent_df_ec3$EC <- rownames(percent_df_ec3)
combined_ec3 <- left_join(count_df_ec3, percent_df_ec3, by = "EC")
write.csv(combined_ec3, paste0("EC_Function_vs_Genus_", target_sample_ec3, ".csv"), row.names = FALSE)




############### SAMPLE NAME = GHEE ############################
# === STEP 0: SELECT SAMPLE 4 ===
target_sample_ec4 <- "sample4"  # Replace with actual sample name
target_data_ec4 <- ec_taxa_contrib %>% filter(sample == target_sample_ec4)

# === STEP 1: Identify Top 10 EC functions ===
top_ec_4 <- target_data_ec4 %>%
  group_by(`function`) %>%
  summarise(total = sum(taxon_function_abun), .groups = "drop") %>%
  slice_max(total, n = 10)

# === STEP 2: Filter data to top ECs ===
filtered_data_ec4 <- target_data_ec4 %>%
  filter(`function` %in% top_ec_4$`function`, !is.na(Genus))

# === STEP 3: Identify Top 10 Genera ===
top_genus_ec4 <- filtered_data_ec4 %>%
  group_by(Genus) %>%
  summarise(total = sum(taxon_function_abun), .groups = "drop") %>%
  slice_max(total, n = 10)

# === STEP 4: Filter data to top genera ===
filtered_data_ec4 <- filtered_data_ec4 %>%
  filter(Genus %in% top_genus_ec4$Genus)

# === STEP 5: Create EC × Genus matrix ===
heatmap_df_ec4 <- filtered_data_ec4 %>%
  group_by(`function`, Genus) %>%
  summarise(total = sum(taxon_function_abun), .groups = "drop") %>%
  pivot_wider(names_from = Genus, values_from = total, values_fill = 0)

heatmap_matrix_ec4 <- heatmap_df_ec4 %>%
  column_to_rownames("function") %>%
  as.matrix()

# === STEP 6: Plot Heatmap ===
pheatmap(heatmap_matrix_ec4,
         main = paste("Top EC Functions in Ghee"),
         color = colorRampPalette(c("pink", "white", "red"))(100),
         scale = "row",
         fontsize_row = 10,
         fontsize_col = 10,
         border_color = NA,
         cluster_rows = TRUE,
         cluster_cols = TRUE)

# === STEP 7: Save Combined CSV (Counts + Percentages) ===
percent_mat_ec4 <- sweep(heatmap_matrix_ec4, 1, rowSums(heatmap_matrix_ec4), FUN = "/") * 100
count_df_ec4 <- as.data.frame(heatmap_matrix_ec4)
percent_df_ec4 <- as.data.frame(percent_mat_ec4)
colnames(count_df_ec4) <- paste0(colnames(count_df_ec4), "_count")
colnames(percent_df_ec4) <- paste0(colnames(percent_df_ec4), "_percent")
count_df_ec4$EC <- rownames(count_df_ec4)
percent_df_ec4$EC <- rownames(percent_df_ec4)
combined_ec4 <- left_join(count_df_ec4, percent_df_ec4, by = "EC")
write.csv(combined_ec4, paste0("EC_Function_vs_Genus_", target_sample_ec4, ".csv"), row.names = FALSE)




############### SAMPLE NAME = NUNU ############################
# === STEP 0: SELECT SAMPLE 5 ===
target_sample_ec5 <- "sample5"  # Replace with your actual sample name
target_data_ec5 <- ec_taxa_contrib %>% filter(sample == target_sample_ec5)

# === STEP 1: Identify Top 10 EC functions ===
top_ec_5 <- target_data_ec5 %>%
  group_by(`function`) %>%
  summarise(total = sum(taxon_function_abun), .groups = "drop") %>%
  slice_max(total, n = 10)

# === STEP 2: Filter data to top ECs ===
filtered_data_ec5 <- target_data_ec5 %>%
  filter(`function` %in% top_ec_5$`function`, !is.na(Genus))

# === STEP 3: Identify Top 10 Genera ===
top_genus_ec5 <- filtered_data_ec5 %>%
  group_by(Genus) %>%
  summarise(total = sum(taxon_function_abun), .groups = "drop") %>%
  slice_max(total, n = 10)

# === STEP 4: Filter data to top genera ===
filtered_data_ec5 <- filtered_data_ec5 %>%
  filter(Genus %in% top_genus_ec5$Genus)

# === STEP 5: Create EC × Genus matrix ===
heatmap_df_ec5 <- filtered_data_ec5 %>%
  group_by(`function`, Genus) %>%
  summarise(total = sum(taxon_function_abun), .groups = "drop") %>%
  pivot_wider(names_from = Genus, values_from = total, values_fill = 0)

heatmap_matrix_ec5 <- heatmap_df_ec5 %>%
  column_to_rownames("function") %>%
  as.matrix()

# === STEP 6: Plot EC Heatmap ===
pheatmap(heatmap_matrix_ec5,
         main = paste("Top EC Functions in Nunu"),
         color = colorRampPalette(c("orange", "white", "red"))(100),
         scale = "row",
         fontsize_row = 10,
         fontsize_col = 10,
         border_color = NA,
         cluster_rows = TRUE,
         cluster_cols = TRUE)

# === STEP 7: Save Combined CSV (Counts + Percentages) ===
percent_mat_ec5 <- sweep(heatmap_matrix_ec5, 1, rowSums(heatmap_matrix_ec5), FUN = "/") * 100
count_df_ec5 <- as.data.frame(heatmap_matrix_ec5)
percent_df_ec5 <- as.data.frame(percent_mat_ec5)
colnames(count_df_ec5) <- paste0(colnames(count_df_ec5), "_count")
colnames(percent_df_ec5) <- paste0(colnames(percent_df_ec5), "_percent")
count_df_ec5$EC <- rownames(count_df_ec5)
percent_df_ec5$EC <- rownames(percent_df_ec5)
combined_ec5 <- left_join(count_df_ec5, percent_df_ec5, by = "EC")
write.csv(combined_ec5, paste0("EC_Function_vs_Genus_", target_sample_ec5, ".csv"), row.names = FALSE)


################################################################################
############ PATHWAYS PATHWAYS SAMPLE SPECIFIC ################################
############### SAMPLE NAME = NONO ############################

# === STEP 0: SELECT SAMPLE 1 ===
target_sample_pw1 <- "sample1"
target_data_pw1 <- pathway_taxa_contrib %>% filter(sample == target_sample_pw1)

# === STEP 1: Identify Top 10 Pathways ===
top_pathway_pw1 <- target_data_pw1 %>%
  group_by(`function`) %>%
  summarise(total = sum(taxon_function_abun), .groups = "drop") %>%
  slice_max(total, n = 10)

# === STEP 2: Filter data to top pathways ===
filtered_data_pw1 <- target_data_pw1 %>%
  filter(`function` %in% top_pathway_pw1$`function`, !is.na(Genus))

# === STEP 3: Identify Top 10 Genera ===
top_genus_pw1 <- filtered_data_pw1 %>%
  group_by(Genus) %>%
  summarise(total = sum(taxon_function_abun), .groups = "drop") %>%
  slice_max(total, n = 10)

# === STEP 4: Filter to top genera ===
filtered_data_pw1 <- filtered_data_pw1 %>%
  filter(Genus %in% top_genus_pw1$Genus)

# === STEP 5: Build Pathway × Genus matrix ===
heatmap_df_pw1 <- filtered_data_pw1 %>%
  group_by(`function`, Genus) %>%
  summarise(total = sum(taxon_function_abun), .groups = "drop") %>%
  pivot_wider(names_from = Genus, values_from = total, values_fill = 0)

heatmap_matrix_pw1 <- heatmap_df_pw1 %>%
  column_to_rownames("function") %>%
  as.matrix()

# === STEP 6: Plot Heatmap ===
pheatmap(heatmap_matrix_pw1,
         main = paste("Top Pathways in Nono"),
         color = colorRampPalette(brewer.pal(9, "BuPu"))(100),
         scale = "row",
         fontsize_row = 10,
         fontsize_col = 10,
         border_color = NA,
         cluster_rows = TRUE,
         cluster_cols = TRUE)

# === STEP 7: Save Combined CSV (Counts + Percentages) ===
percent_mat_pw1 <- sweep(heatmap_matrix_pw1, 1, rowSums(heatmap_matrix_pw1), FUN = "/") * 100
count_df_pw1 <- as.data.frame(heatmap_matrix_pw1)
percent_df_pw1 <- as.data.frame(percent_mat_pw1)
colnames(count_df_pw1) <- paste0(colnames(count_df_pw1), "_count")
colnames(percent_df_pw1) <- paste0(colnames(percent_df_pw1), "_percent")
count_df_pw1$Pathway <- rownames(count_df_pw1)
percent_df_pw1$Pathway <- rownames(percent_df_pw1)
combined_pw1 <- left_join(count_df_pw1, percent_df_pw1, by = "Pathway")
write.csv(combined_pw1, paste0("Pathway_Function_vs_Genus_", target_sample_pw1, ".csv"), row.names = FALSE)



############### SAMPLE NAME = NONO ############################

# === STEP 0: SELECT SAMPLE 2 ===
target_sample_pw2 <- "sample2"
target_data_pw2 <- pathway_taxa_contrib %>% filter(sample == target_sample_pw2)

# === STEP 1: Identify Top 10 Pathways ===
top_pathway_pw2 <- target_data_pw2 %>%
  group_by(`function`) %>%
  summarise(total = sum(taxon_function_abun), .groups = "drop") %>%
  slice_max(total, n = 10)

# === STEP 2: Filter data to top pathways ===
filtered_data_pw2 <- target_data_pw2 %>%
  filter(`function` %in% top_pathway_pw2$`function`, !is.na(Genus))

# === STEP 3: Identify Top 10 Genera ===
top_genus_pw2 <- filtered_data_pw2 %>%
  group_by(Genus) %>%
  summarise(total = sum(taxon_function_abun), .groups = "drop") %>%
  slice_max(total, n = 10)

# === STEP 4: Filter to top genera ===
filtered_data_pw2 <- filtered_data_pw2 %>%
  filter(Genus %in% top_genus_pw2$Genus)

# === STEP 5: Build Pathway × Genus matrix ===
heatmap_df_pw2 <- filtered_data_pw2 %>%
  group_by(`function`, Genus) %>%
  summarise(total = sum(taxon_function_abun), .groups = "drop") %>%
  pivot_wider(names_from = Genus, values_from = total, values_fill = 0)

heatmap_matrix_pw2 <- heatmap_df_pw2 %>%
  column_to_rownames("function") %>%
  as.matrix()

# === STEP 6: Plot Heatmap ===
pheatmap(heatmap_matrix_pw2,
         main = paste("Top Pathways in Wara"),
         color = colorRampPalette(brewer.pal(9, "Reds"))(100),
         scale = "row",
         fontsize_row = 10,
         fontsize_col = 10,
         border_color = NA,
         cluster_rows = TRUE,
         cluster_cols = TRUE)

# === STEP 7: Save Combined CSV (Counts + Percentages) ===
percent_mat_pw2 <- sweep(heatmap_matrix_pw2, 1, rowSums(heatmap_matrix_pw2), FUN = "/") * 100
count_df_pw2 <- as.data.frame(heatmap_matrix_pw2)
percent_df_pw2 <- as.data.frame(percent_mat_pw2)
colnames(count_df_pw2) <- paste0(colnames(count_df_pw2), "_count")
colnames(percent_df_pw2) <- paste0(colnames(percent_df_pw2), "_percent")
count_df_pw2$Pathway <- rownames(count_df_pw2)
percent_df_pw2$Pathway <- rownames(percent_df_pw2)
combined_pw2 <- left_join(count_df_pw2, percent_df_pw2, by = "Pathway")
write.csv(combined_pw2, paste0("Pathway_Function_vs_Genus_", target_sample_pw2, ".csv"), row.names = FALSE)



############### SAMPLE NAME = KWERIONIK ############################

# === STEP 0: SELECT SAMPLE 3 ===
target_sample_pw3 <- "sample3"
target_data_pw3 <- pathway_taxa_contrib %>% filter(sample == target_sample_pw3)

# === STEP 1: Identify Top 10 Pathways ===
top_pathway_pw3 <- target_data_pw3 %>%
  group_by(`function`) %>%
  summarise(total = sum(taxon_function_abun), .groups = "drop") %>%
  slice_max(total, n = 10)

# === STEP 2: Filter data to top pathways ===
filtered_data_pw3 <- target_data_pw3 %>%
  filter(`function` %in% top_pathway_pw3$`function`, !is.na(Genus))

# === STEP 3: Identify Top 10 Genera ===
top_genus_pw3 <- filtered_data_pw3 %>%
  group_by(Genus) %>%
  summarise(total = sum(taxon_function_abun), .groups = "drop") %>%
  slice_max(total, n = 10)

# === STEP 4: Filter to top genera ===
filtered_data_pw3 <- filtered_data_pw3 %>%
  filter(Genus %in% top_genus_pw3$Genus)

# === STEP 5: Build Pathway × Genus matrix ===
heatmap_df_pw3 <- filtered_data_pw3 %>%
  group_by(`function`, Genus) %>%
  summarise(total = sum(taxon_function_abun), .groups = "drop") %>%
  pivot_wider(names_from = Genus, values_from = total, values_fill = 0)

heatmap_matrix_pw3 <- heatmap_df_pw3 %>%
  column_to_rownames("function") %>%
  as.matrix()

# === STEP 6: Plot Heatmap ===
pheatmap(heatmap_matrix_pw3,
         main = paste("Top Pathways in Kwerionik"),
         color = colorRampPalette(brewer.pal(9, "Purples"))(100),
         scale = "row",
         fontsize_row = 10,
         fontsize_col = 10,
         border_color = NA,
         cluster_rows = TRUE,
         cluster_cols = TRUE)

# === STEP 7: Save Combined CSV (Counts + Percentages) ===
percent_mat_pw3 <- sweep(heatmap_matrix_pw3, 1, rowSums(heatmap_matrix_pw3), FUN = "/") * 100
count_df_pw3 <- as.data.frame(heatmap_matrix_pw3)
percent_df_pw3 <- as.data.frame(percent_mat_pw3)
colnames(count_df_pw3) <- paste0(colnames(count_df_pw3), "_count")
colnames(percent_df_pw3) <- paste0(colnames(percent_df_pw3), "_percent")
count_df_pw3$Pathway <- rownames(count_df_pw3)
percent_df_pw3$Pathway <- rownames(percent_df_pw3)
combined_pw3 <- left_join(count_df_pw3, percent_df_pw3, by = "Pathway")
write.csv(combined_pw3, paste0("Pathway_Function_vs_Genus_", target_sample_pw3, ".csv"), row.names = FALSE)




############### SAMPLE NAME = GHEE ############################
# === STEP 0: SELECT SAMPLE 4 ===
target_sample_pw4 <- "sample4"
target_data_pw4 <- pathway_taxa_contrib %>% filter(sample == target_sample_pw4)

# === STEP 1: Identify Top 10 Pathways ===
top_pathway_pw4 <- target_data_pw4 %>%
  group_by(`function`) %>%
  summarise(total = sum(taxon_function_abun), .groups = "drop") %>%
  slice_max(total, n = 10)

# === STEP 2: Filter to top pathways ===
filtered_data_pw4 <- target_data_pw4 %>%
  filter(`function` %in% top_pathway_pw4$`function`, !is.na(Genus))

# === STEP 3: Identify Top 10 Genera ===
top_genus_pw4 <- filtered_data_pw4 %>%
  group_by(Genus) %>%
  summarise(total = sum(taxon_function_abun), .groups = "drop") %>%
  slice_max(total, n = 10)

# === STEP 4: Filter to top genera ===
filtered_data_pw4 <- filtered_data_pw4 %>%
  filter(Genus %in% top_genus_pw4$Genus)

# === STEP 5: Build Pathway × Genus matrix ===
heatmap_df_pw4 <- filtered_data_pw4 %>%
  group_by(`function`, Genus) %>%
  summarise(total = sum(taxon_function_abun), .groups = "drop") %>%
  pivot_wider(names_from = Genus, values_from = total, values_fill = 0)

heatmap_matrix_pw4 <- heatmap_df_pw4 %>%
  column_to_rownames("function") %>%
  as.matrix()

# === STEP 6: Plot Heatmap ===
pheatmap(heatmap_matrix_pw4,
         main = paste("Top Pathways in Ghee"),
         color = colorRampPalette(brewer.pal(9, "Oranges"))(100),
         scale = "row",
         fontsize_row = 10,
         fontsize_col = 10,
         border_color = NA,
         cluster_rows = TRUE,
         cluster_cols = TRUE)

# === STEP 7: Save Combined CSV (Counts + Percentages) ===
percent_mat_pw4 <- sweep(heatmap_matrix_pw4, 1, rowSums(heatmap_matrix_pw4), FUN = "/") * 100
count_df_pw4 <- as.data.frame(heatmap_matrix_pw4)
percent_df_pw4 <- as.data.frame(percent_mat_pw4)
colnames(count_df_pw4) <- paste0(colnames(count_df_pw4), "_count")
colnames(percent_df_pw4) <- paste0(colnames(percent_df_pw4), "_percent")
count_df_pw4$Pathway <- rownames(count_df_pw4)
percent_df_pw4$Pathway <- rownames(percent_df_pw4)
combined_pw4 <- left_join(count_df_pw4, percent_df_pw4, by = "Pathway")
write.csv(combined_pw4, paste0("Pathway_Function_vs_Genus_", target_sample_pw4, ".csv"), row.names = FALSE)


############### SAMPLE NAME = NUNU ############################
# === STEP 0: SELECT SAMPLE 5 ===

target_sample_pw5 <- "sample5"
target_data_pw5 <- pathway_taxa_contrib %>% filter(sample == target_sample_pw5)

# === STEP 1: Identify Top 10 Pathways ===
top_pathway_pw5 <- target_data_pw5 %>%
  group_by(`function`) %>%
  summarise(total = sum(taxon_function_abun), .groups = "drop") %>%
  slice_max(total, n = 10)

# === STEP 2: Filter to top pathways ===
filtered_data_pw5 <- target_data_pw5 %>%
  filter(`function` %in% top_pathway_pw5$`function`, !is.na(Genus))

# === STEP 3: Identify Top 10 Genera ===
top_genus_pw5 <- filtered_data_pw5 %>%
  group_by(Genus) %>%
  summarise(total = sum(taxon_function_abun), .groups = "drop") %>%
  slice_max(total, n = 10)

# === STEP 4: Filter to top genera ===
filtered_data_pw5 <- filtered_data_pw5 %>%
  filter(Genus %in% top_genus_pw5$Genus)

# === STEP 5: Build Pathway × Genus matrix ===
heatmap_df_pw5 <- filtered_data_pw5 %>%
  group_by(`function`, Genus) %>%
  summarise(total = sum(taxon_function_abun), .groups = "drop") %>%
  pivot_wider(names_from = Genus, values_from = total, values_fill = 0)

heatmap_matrix_pw5 <- heatmap_df_pw5 %>%
  column_to_rownames("function") %>%
  as.matrix()

# === STEP 6: Plot Heatmap ===
pheatmap(heatmap_matrix_pw5,
         main = paste("Top Pathways in Nunu"),
         color = colorRampPalette(brewer.pal(9, "RdPu"))(100),
         scale = "row",
         fontsize_row = 10,
         fontsize_col = 10,
         border_color = NA,
         cluster_rows = TRUE,
         cluster_cols = TRUE)

# === STEP 7: Save Combined CSV (Counts + Percentages) ===
percent_mat_pw5 <- sweep(heatmap_matrix_pw5, 1, rowSums(heatmap_matrix_pw5), FUN = "/") * 100
count_df_pw5 <- as.data.frame(heatmap_matrix_pw5)
percent_df_pw5 <- as.data.frame(percent_mat_pw5)
colnames(count_df_pw5) <- paste0(colnames(count_df_pw5), "_count")
colnames(percent_df_pw5) <- paste0(colnames(percent_df_pw5), "_percent")
count_df_pw5$Pathway <- rownames(count_df_pw5)
percent_df_pw5$Pathway <- rownames(percent_df_pw5)
combined_pw5 <- left_join(count_df_pw5, percent_df_pw5, by = "Pathway")
write.csv(combined_pw5, paste0("Pathway_Function_vs_Genus_", target_sample_pw5, ".csv"), row.names = FALSE)


########################################################################
#######################################################################
################ GRID GRID ko ko ####################################

# Load required libraries
library(pheatmap)
library(gridExtra)
library(grid)

# === Grab each heatmap as grob ===

# Sample 1 – Nono
heatmap_plot_1 <- grid.grabExpr(pheatmap(heatmap_matrix,
                                         main = "Top KO Functions in Nono",
                                         color = colorRampPalette(brewer.pal(9, "YlGnBu"))(100),
                                         scale = "row", fontsize_row = 10, fontsize_col = 10,
                                         border_color = NA, cluster_rows = TRUE, cluster_cols = TRUE))

# Sample 2 – Wara
heatmap_plot_2 <- grid.grabExpr(pheatmap(heatmap_matrix_2,
                                         main = "Top KO Functions in Wara",
                                         color = colorRampPalette(c("black", "yellow", "red"))(100),
                                         scale = "row", fontsize_row = 10, fontsize_col = 10,
                                         border_color = NA, cluster_rows = TRUE, cluster_cols = TRUE))

# Sample 3 – Kwerionik
heatmap_plot_3 <- grid.grabExpr(pheatmap(heatmap_matrix_3,
                                         main = "Top KO Functions in Kwerionik",
                                         color = colorRampPalette(c("pink", "white", "red"))(100),
                                         scale = "row", fontsize_row = 10, fontsize_col = 10,
                                         border_color = NA, cluster_rows = TRUE, cluster_cols = TRUE))

# Sample 4 – Ghee
heatmap_plot_4 <- grid.grabExpr(pheatmap(heatmap_matrix_4,
                                         main = "Top KO Functions in Ghee",
                                         color = colorRampPalette(c("black", "white", "red"))(100),
                                         scale = "row", fontsize_row = 10, fontsize_col = 10,
                                         border_color = NA, cluster_rows = TRUE, cluster_cols = TRUE))

# Sample 5 – Nunu
heatmap_plot_5 <- grid.grabExpr(pheatmap(heatmap_matrix_5,
                                         main = "Top KO Functions in Nunu",
                                         color = colorRampPalette(c("blue", "white", "red"))(100),
                                         scale = "row", fontsize_row = 10, fontsize_col = 10,
                                         border_color = NA, cluster_rows = TRUE, cluster_cols = TRUE))

# Create an empty grob for slot 6 (optional)
empty_plot <- grid.rect(gp = gpar(col = NA))  # Transparent placeholder


# Combine plots 1, 2, and 5 into a single frame
frame_1 <- grid.arrange(
  heatmap_plot_1,  # Nono
  heatmap_plot_2,  # Wara
  heatmap_plot_5,  # Nunu
  ncol = 3
)


# Combine plots 3 and 4 into another frame
frame_2 <- grid.arrange(
  heatmap_plot_3,  # Kwerionik
  heatmap_plot_4,  # Ghee
  ncol = 2
)



#######################################################################
################################## grid grid ec ec #####################


heatmap_ec_plot_1 <- grid.grabExpr(
  pheatmap(heatmap_matrix_ec1,
           main = "Top EC Functions in Nono",
           color = colorRampPalette(brewer.pal(9, "YlOrRd"))(100),
           scale = "row",
           fontsize_row = 10,
           fontsize_col = 10,
           angle_col = 45,
           border_color = NA,
           cluster_rows = TRUE,
           cluster_cols = TRUE)
)


heatmap_ec_plot_2 <- grid.grabExpr(
  pheatmap(heatmap_matrix_ec2,
           main = "Top EC Functions in Wara",
           color = colorRampPalette(brewer.pal(9, "YlGnBu"))(100),
           scale = "row",
           fontsize_row = 10,
           fontsize_col = 10,
           angle_col = 45,
           border_color = NA,
           cluster_rows = TRUE,
           cluster_cols = TRUE)
)


heatmap_ec_plot_3 <- grid.grabExpr(
  pheatmap(heatmap_matrix_ec3,
           main = "Top EC Functions in Kwerionik",
           color = colorRampPalette(c("navy", "white", "red"))(100),
           scale = "row",
           fontsize_row = 10,
           fontsize_col = 10,
           angle_col = 45,
           border_color = NA,
           cluster_rows = TRUE,
           cluster_cols = TRUE)
)


heatmap_ec_plot_4 <- grid.grabExpr(
  pheatmap(heatmap_matrix_ec4,
           main = "Top EC Functions in Ghee",
           color = colorRampPalette(c("black", "yellow", "red"))(100),
           scale = "row",
           fontsize_row = 10,
           fontsize_col = 10,
           angle_col = 45,
           border_color = NA,
           cluster_rows = TRUE,
           cluster_cols = TRUE)
)

heatmap_ec_plot_5 <- grid.grabExpr(
  pheatmap(heatmap_matrix_ec5,
           main = "Top EC Functions in Nunu",
           color = colorRampPalette(c("orange", "white", "red"))(100),
           scale = "row",
           fontsize_row = 10,
           fontsize_col = 10,
           angle_col = 45,
           border_color = NA,
           cluster_rows = TRUE,
           cluster_cols = TRUE)
)

# Frame 1
grid.arrange(heatmap_ec_plot_1, heatmap_ec_plot_2, heatmap_ec_plot_5, ncol = 3)

# Frame 2
grid.arrange(heatmap_ec_plot_3, heatmap_ec_plot_4, ncol = 2)



#############################################################
########################### grid pathway 



# === Create all 5 pathway heatmaps and capture gtables ===

# Sample 1 - Nono
p1 <- pheatmap(heatmap_matrix_pw1,
               main = "Top Pathways in Nono",
               color = colorRampPalette(brewer.pal(9, "BuPu"))(100),
               scale = "row",
               fontsize_row = 10,
               fontsize_col = 10,
               border_color = NA,
               cluster_rows = TRUE,
               cluster_cols = TRUE)
g1 <- p1$gtable

# Sample 2 - Wara
p2 <- pheatmap(heatmap_matrix_pw2,
               main = "Top Pathways in Wara",
               color = colorRampPalette(brewer.pal(9, "Reds"))(100),
               scale = "row",
               fontsize_row = 10,
               fontsize_col = 10,
               border_color = NA,
               cluster_rows = TRUE,
               cluster_cols = TRUE)
g2 <- p2$gtable

# Sample 3 - Kwerionik
p3 <- pheatmap(heatmap_matrix_pw3,
               main = "Top Pathways in Kwerionik",
               color = colorRampPalette(brewer.pal(9, "Purples"))(100),
               scale = "row",
               fontsize_row = 10,
               fontsize_col = 10,
               border_color = NA,
               cluster_rows = TRUE,
               cluster_cols = TRUE)
g3 <- p3$gtable

# Sample 4 - Ghee
p4 <- pheatmap(heatmap_matrix_pw4,
               main = "Top Pathways in Ghee",
               color = colorRampPalette(brewer.pal(9, "Oranges"))(100),
               scale = "row",
               fontsize_row = 10,
               fontsize_col = 10,
               border_color = NA,
               cluster_rows = TRUE,
               cluster_cols = TRUE)
g4 <- p4$gtable

# Sample 5 - Nunu
p5 <- pheatmap(heatmap_matrix_pw5,
               main = "Top Pathways in Nunu",
               color = colorRampPalette(brewer.pal(9, "RdPu"))(100),
               scale = "row",
               fontsize_row = 10,
               fontsize_col = 10,
               border_color = NA,
               cluster_rows = TRUE,
               cluster_cols = TRUE)
g5 <- p5$gtable

# === Arrange into frames ===

# Frame 1: Nono, Wara, Nunu
frame1 <- grid.arrange(g1, g2, g5, ncol = 3)

# Frame 2: Kwerionik, Ghee
frame2 <- grid.arrange(g3, g4, ncol = 2)

# === Optional: Save to file ===
# ggsave("frame1_pathways.png", plot = frame1, width = 18, height = 6, dpi = 300)
# ggsave("frame2_pathways.png", plot = frame2, width = 12, height = 6, dpi = 300)




























