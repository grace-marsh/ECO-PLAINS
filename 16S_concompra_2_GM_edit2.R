


library(dplyr)
library(tidyverse)
library(phyloseq)
library(ggplot2)
library(biomformat)
library(microViz)
library(vegan)
library(microbiome)
library(reshape)
library(wesanderson)
library(pheatmap)
library(RColorBrewer)
library(svglite)
library(ape)
library(ggtree) # phylogeny
library(AICcmodavg) # anova analysis
library(multcompView) 
library(ggpubr)
library(compositions)
library(pairwiseAdonis) # permanova
library(ggcorrplot)
library(factoextra)
library(ANCOMBC) # da analysis
library(DT)
library(WGCNA) # taxa clusters
library(htmlwidgets)
library(iCAMP) # assembly processes
library(rstatix)
library(picante) # faith pd



setwd("C:/Users/gmarsh/OneDrive - epfl.ch/grace_shared/greenland_2023_16S_analysis/Grace_analysis/filter_1000")



# Import and Format Data -------------------------------------------------------------


# import metadata and format
metadata <- read.csv("gl_2023_metadata_edit.csv", sep = ",")
metadata$Sample <- gsub("T", "_T", metadata$Sample)
metadata$Sample <- gsub("L", "_L", metadata$Sample)
row.names(metadata) <- metadata$Sample
metadata[metadata$Depth == "Top" ,]$Soil_temp <- metadata[metadata$Depth == "Top" ,]$Surface_temp # format soil temp as one variable
metadata$Surface_temp = NULL

# import otu table
otu_table <- read.csv("../otu_table_concompra.csv", header = T, check.names = F, sep = ",", row.names = 1) # change *OTU ID to OTU_ID header in file
sum(otu_table) # 1109441

# format col names otu
colnames(otu_table) <- gsub(".CONCOMPRA", "", colnames(otu_table))
barcodes <- read.csv("../sample_barcodes.tsv", header = T, sep = "\t")
otu_table <- otu_table[, barcodes$Barcode] # select only sample barcodes
for (i in 1:length(colnames(otu_table))){
  colnames(otu_table)[i] <- barcodes[barcodes$Barcode == colnames(otu_table)[i], ]$X.SampleID
}



# Filtering ---------------------------------------------------------------


# filter samples
otumat_filt <- otu_table[(otu_table$extraction_control == 0) & (otu_table$pcr_water_control == 0),] # remove contaminants from control
otumat_filt$extraction_control = NULL # remove controls
otumat_filt$pcr_water_control = NULL
otumat_filt$'3C_L' = NULL # remove FC outlier
otumat_filt <- otumat_filt[, colSums(otumat_filt) > 1000] # remove samples with few reads (49 samples)

# format taxonomy
taxa <- read.csv("../clustered_consensus_silva_concompra.tsv", header = T, sep = "\t", row.names = 1)
colnames(taxa)[1] <- "Taxonomy"

taxonomy <- taxa %>% dplyr::select(Taxonomy) %>% 
  separate(Taxonomy, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), ";")  %>% 
  as.matrix()
taxonomy <- gsub("d__","",
                 gsub("p__","",
                      gsub("o__","",
                           gsub("c__","",
                                gsub("g__","",
                                     gsub("s__","",
                                          gsub("f__","", taxonomy)))))))  %>% 
  as.matrix()


# create phyloseq object
otu <- otu_table(otumat_filt, taxa_are_rows = TRUE)
sample_data <- sample_data(metadata)
tax_table <- tax_table(taxonomy)
otu_all <- merge_phyloseq(otu, sample_data, tax_table) 
# 2841 otus

# filter chloroplast, mitochondria, eukaryote
otu_clean <- subset_taxa(otu_all, (Order != "Chloroplast") | is.na(Order))
otu_clean <- subset_taxa(otu_clean, (Family != "Mitochondria") | is.na(Family)) 
otu_clean <- subset_taxa(otu_clean, (Kingdom == "Bacteria") | is.na(Kingdom)) 
# 2767 otus

# Remove singletons
otu_filt <- filter_taxa(otu_clean, function(x) sum(x) > 1, TRUE) # remove otus with less than 1 read
colSums(otu_table(otu_filt))/colSums(otu_table(otu_clean)) # check proportion after filtering
otu_filt # 2402 otus

# microViz check
otu_check <- phyloseq_validate(otu_filt, remove_undetected = TRUE) # removes 0 counts

# rename unassigned phyla
tax_table(otu_filt)[tax_table(otu_filt)[,"Phylum"]== "","Phylum"] <- "Unassigned" # rename NAs
tax_table(otu_filt)[is.na(tax_table(otu_filt)[,"Phylum"]),"Phylum"] <- "Unassigned" # rename NAs

# filter again samples with more than 1000 reads after this filtering
otu_table(otu_filt) <- otu_table(otu_filt)[,colSums(otu_table(otu_filt)) > 1000]

# stats
sum(otu_table(otu_filt)) # 848846 reads after fitlering
length(unique(tax_table(otu_filt)[,"Genus"])) # 288 unique genera
otu_filt # 49 samples

# Rarefaction curve after filtering
otus_filtered <- t(otu_table(otu_filt))
class(otus_filtered) <- "matrix"
rarecurve(otus_filtered, step=1000, ylab="OTU", label=T, cex.axis=1.5, cex.lab=1.5)





# Phylogeny ---------------------------------------------------------------


# for phylogeny tree construction
# input fasta sequences and subset filtered otus
#  library("Biostrings")
#  otu_seqs <- readDNAStringSet("../clustered_consensus.fasta")
#  length(names(otu_seqs)) # check length
# #
#  otu_seqs_filtered <- otu_seqs[names(otu_seqs) %in% row.names(otu_table(otu_filt))] # gets seqs for filtered otus only
#  length(names(otu_seqs_filtered)) # check length
# 
#  writeXStringSet(otu_seqs_filtered, "otu_concompra_filtered1000.fasta", append = FALSE, compress = FALSE, format = "fasta") # write out fasta file

# phylogeny analysis
tree <- read.tree("rooted_tree.nwk")
otu_all <- merge_phyloseq(otu_filt, tree) # add tree to phyloseq object
ntaxa(otu_all) # check taxa match


# plot circular tree with site labels
df1 <- data.frame(row.names = c("0","1","1_VEG", "2","3","4","5"), No = c(0:6)) # create numerical values for each site to plot x value

gg_color_hue <- function(n) { # default colours from colourwheel (for phyla)
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

ggtree(otu_all, branch.length='none', layout = "circular", aes(colour = Phylum)) +
  geom_point(aes(x=x + df1[Site,], colour = Site), na.rm=TRUE) + scale_colour_manual(values = c(wes_palette("Zissou1", 7, type = "continuous"), gg_color_hue(100)))
# ggsave("phyla_tree.png", plot = last_plot(), height = 10, width = 10)




# Abundance ---------------------------------------------------------------


# Read abundance plot
filt_data <- psmelt(otu_filt) # create dataframe from phyloseq object

# order phyla
filt_phyla <- aggregate_taxa(otu_filt, level = "Phylum") # phylum level data
filt_data$Phylum <- factor(filt_data$Phylum, levels = top_taxa(filt_phyla)) # order by phylum abundance

# top 10 phyla only
other <- subset(filt_data, !(filt_data$Phylum %in% top_taxa(filt_phyla, 10))) # 
other$Phylum <- "Other" # create other category
top_data <- rbind(other, filt_data[filt_data$Phylum %in% top_taxa(filt_phyla, 10), ]) # combine other and top 10 phyla
top_data$Phylum <- factor(top_data$Phylum, levels = c(top_taxa(filt_phyla, 10), "Other")) # order

ggplot(top_data, aes(x = Sample, y = Abundance, fill = Phylum)) + geom_bar(aes(fill = Phylum), stat="identity", position="stack") + 
  scale_fill_manual(values = wes_palette("Zissou1", 11, type = "continuous")) + theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
#ggsave("top_phyla_abundance.png", plot = last_plot(), height = 7, width = 10)


# relative abundance plot
top_data_formatted <- top_data %>% group_by(Sample, Site, Depth, Phylum) %>% summarise(relative_abundance = sum(Abundance))
ggplot(top_data_formatted, aes(x = Sample, y = relative_abundance, fill = Phylum)) + geom_bar(aes(fill = Phylum), stat="identity", position="fill", colour = "black") + 
  scale_fill_manual(values = wes_palette("Zissou1", 11, type = "continuous")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + coord_flip() + facet_grid(rows = vars(Site, Depth), scales = "free", space = "free") 


# Average for each depth and site
top_data_formatted_2 <- top_data_formatted %>% group_by(Site, Depth, Phylum) %>% summarise(relative_abundance_2 = median(relative_abundance))

top_data_formatted_2[top_data_formatted_2$Site == "0", ]$Site <- "Glacier" # rename sites for plotting
top_data_formatted_2[top_data_formatted_2$Site == "1", ]$Site <- "Site 1"
top_data_formatted_2[top_data_formatted_2$Site == "1_VEG", ]$Site <- "Site 1 Veg"
top_data_formatted_2[top_data_formatted_2$Site == "2", ]$Site <- "Site 2"
top_data_formatted_2[top_data_formatted_2$Site == "3", ]$Site <- "Site 3"
top_data_formatted_2[top_data_formatted_2$Site == "4", ]$Site <- "Site 4"
top_data_formatted_2[top_data_formatted_2$Site == "5", ]$Site <- "Site 5"

# new phyla names
levels(top_data_formatted_2$Phylum)[match("Proteobacteria",levels(top_data_formatted_2$Phylum))] <- "Pseudomonadota"
levels(top_data_formatted_2$Phylum)[match("Acidobacteria",levels(top_data_formatted_2$Phylum))] <- "Acidobacteriota"
levels(top_data_formatted_2$Phylum)[match("Cyanobacteria",levels(top_data_formatted_2$Phylum))] <- "Cyanobacteriota"
levels(top_data_formatted_2$Phylum)[match("Firmicutes",levels(top_data_formatted_2$Phylum))] <- "Bacillota"
levels(top_data_formatted_2$Phylum)[match("Chloroflexi",levels(top_data_formatted_2$Phylum))] <- "Chloroflexota"
levels(top_data_formatted_2$Phylum)[match("Actinobacteriota",levels(top_data_formatted_2$Phylum))] <- "Actinomycetota"

# plot relative abundance average for sites
rel_abund <- ggplot(top_data_formatted_2, aes(x = Depth, y = relative_abundance_2, fill = Phylum)) + geom_bar(aes(fill = Phylum), stat="identity", position="fill", colour = "black") + 
  scale_fill_manual(values = c(wes_palette("Zissou1", 10, type = "continuous"), "darkgrey")) + theme_bw() + ylab("Relative Abundance") + xlab("") +
  coord_flip() + facet_grid(rows = vars(Site), scales = "free", space = "free") 








# Flow cytometry transformation
# formula = otumat_filt * ( fc_values / colSums(otumat_filt) ) solution from: https://github.com/joey711/phyloseq/issues/497
otu_all_norm <- transform_sample_counts(otu_all, function(x){x/sum(x)}) # normalise count per sample (relative count)
otu_all_fc  <- otu_all_norm

for(n in 1:nsamples(otu_all_norm)){ # for each sample
  otu_table(otu_all_fc)[,n] <- otu_table(otu_all_norm)[,n] * sample_data(otu_all_norm)$FC_count [n] # transform count by total sample FC cells
}

# check the sample total equals the cell count
colSums(otu_table(otu_all_fc)) 
sample_data[,c("FC_count", "Sample")]

# plot fc transformed counts
fc <- sample_data[,c("FC_count", "Sample")]
ggplot(fc[colnames(otu_table(otu_all)),], aes(x = Sample, y = FC_count)) + geom_bar(aes(fill = Sample), stat="identity", position="dodge") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# remove outlier sample 3CL

# plot fc transfromed with ordered phyla
filt_data <- psmelt(otu_all_fc) # create dataframe from phyloseq object
filt_phyla <- aggregate_taxa(otu_all_fc, level = "Phylum") # phylum level data
filt_data$Phylum <- factor(filt_data$Phylum, levels = top_taxa(filt_phyla)) # order by phylum abundance

other <- subset(filt_data, !(filt_data$Phylum %in% top_taxa(filt_phyla, 10))) # get top 10 phyla
other$Phylum <- "other" # create other category
top_data <- rbind(other, filt_data[filt_data$Phylum %in% top_taxa(filt_phyla, 10), ]) # combine other and top 10 phyla
top_data$Phylum <- factor(top_data$Phylum, levels = c(top_taxa(filt_phyla, 10), "other")) # order

ggplot(top_data, aes(x = Sample, y = Abundance, fill = Phylum)) + geom_bar(aes(fill = Phylum), stat="identity", position="stack") + 
  scale_fill_manual(values = wes_palette("Zissou1", 11, type = "continuous")) + theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# ggsave("top_phyla_abundance_FC_trans_log.png", plot = last_plot(), height = 7, width = 10, dpi = 900)






# Alpha Diversity ---------------------------------------------------------



# plot shannon diversity
plot_richness(otu_all, measures=c("Shannon"), x = "Site", color = "Site") + geom_boxplot() + theme_bw() + facet_wrap(~Depth, nrow = 2) + scale_color_manual(values = wes_palette("Zissou1", 7, type = "continuous"))  # shannon index only
#ggsave(filename = "shannon_alphadiversity_depth.png", plot = last_plot(), height = 5, width = 5)


# format for anova
shannon <- estimate_richness(otu_filt, measures = "Shannon") # get shannon values
shannon$Sample <- row.names(shannon)
shannon <- merge(metadata, shannon, by = "Sample")


# kruskal wallis + dunn
k_1 <- kruskal.test(Shannon ~ Site, data = shannon)
kruskal.test(Shannon ~ Depth, data = shannon)
dunn <- dunn_test(Shannon ~ Site, data = shannon)

pvals <- dunn$p.adj
names <- paste(dunn$group1, dunn$group2, sep = "-")
names(pvals) <- names
clust_groups <- as.data.frame.list(multcompLetters(pvals))
clust_groups$Site <- row.names(clust_groups)

# format groups for plot
shannon_sum <- shannon %>% group_by(Site) %>% summarise(mean = mean(Shannon), quant = quantile(Shannon, probs = 0.75)) %>% arrange(desc(mean))
shannon_sum <- merge(shannon_sum, clust_groups, by = "Site")

alpha_diversity <- ggplot(shannon) + geom_boxplot(aes(x = Site, y = Shannon, fill = Site)) + 
  geom_text(data = shannon_sum, aes(x = Site, y = quant, label = clust_groups$Letters),  position = position_dodge(width = 0.8), size = 5, vjust=-2) + 
  scale_fill_manual(values = wes_palette("Zissou1", 7, type = "continuous")) + 
  scale_colour_manual(values = c("black", "black", "black", "black", "black", "black", "black"  )) + theme_bw() + 
  scale_x_discrete(labels = c("Glacier", "Site 1", "Site 1 \nVegetated", "Site 2", "Site 3", "Site 4", "Site 5")) +
  theme(legend.position="none", axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + xlab("")




# Faith's phylogenetic diversity (PD)
faith <- pd(t(as.matrix(otu_table(otu_all))), phy_tree(otu_all), include.root = T)

# format
faith$Sample <- row.names(faith)
faith <- merge(metadata, faith, by = "Sample")

# pd is sum of branch legnths
ggplot(faith) + geom_boxplot(aes(x = Site, y = PD, fill = Site)) + 
  scale_fill_manual(values = wes_palette("Zissou1", 7, type = "continuous")) + theme_bw() + ylab("Faith's Phylogenetic Diversity") +
  scale_x_discrete(labels = c("Glacier", "Site 1", "Site 1 \nVegetated", "Site 2", "Site 3", "Site 4", "Site 5")) +
  theme(legend.position="none", axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 


# faith kruskal wallis + dunn
k_2 <- kruskal.test(PD ~ Site, data = faith)
dunn_2 <- dunn_test(PD ~ Site, data = faith)

pvals2 <- dunn_2$p.adj
names2 <- paste(dunn_2$group1, dunn_2$group2, sep = "-")
names(pvals2) <- names2
clust_groups2 <- as.data.frame.list(multcompLetters(pvals2))
clust_groups2$Site <- row.names(clust_groups2)

# format groups for plot
faith_sum <- faith %>% group_by(Site) %>% summarise(mean = mean(PD), quant = quantile(PD, probs = 0.75)) %>% arrange(desc(mean))
faith_sum$groups <- merge(faith_sum, clust_groups2, by = "Site")

ggplot(faith) + geom_boxplot(aes(x = Site, y = PD, fill = Site)) + 
  geom_text(data = faith_sum, aes(x = Site, y = quant, label = clust_groups2$Letters),  position = position_dodge(width = 0.8), size = 5, vjust=-1.5) + 
  scale_fill_manual(values = wes_palette("Zissou1", 7, type = "continuous")) + theme_bw() + ylab("Faith's Phylogenetic Diversity") +
  scale_x_discrete(labels = c("Glacier", "Site 1", "Site 1 \nVegetated", "Site 2", "Site 3", "Site 4", "Site 5")) +
  theme(legend.position="none", axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + xlab("")
# ggsave(filename = "faith_alphadiversity_sig.png", plot = last_plot(), height = 5, width = 5, dpi = 900)





# UniFrac Distance --------------------------------------------------------


# Weighted unifrac appoach
# uses phylogenetic dist information too
unifrac_dist <- UniFrac(otu_all_fc, weighted = TRUE) %>% as.matrix() # for FC normalised
image(unifrac_dist)
sample_meta <- metadata[metadata$Sample %in% colnames(as.matrix(unifrac_dist)), c(2:3)] # format sample metadata
sample_colours <- list(Depth = c("Lower" = "grey", "Top" = "white" ),
                       Site = c("0" = "#3A9AB2", 
                                "1" = "#85B7B9", 
                                "1_VEG" = "#ADC397", 
                                "2" = "#DCCB4E", 
                                "3" = "#E5A208", 
                                "4" = "#ED6E04", 
                                "5" = "#F11B00")) # colours
pheatmap(as.matrix(unifrac_dist), cluster_rows = T, cluster_cols = T,
         color = colorRampPalette(brewer.pal(3,"GnBu"))(100),
         clustering_distance_cols = 'euclidean',
         clustering_distance_rows = 'euclidean',
         clustering_method = 'ward.D',
         annotation_col = sample_meta,
         annotation_colors = sample_colours,
)


# plot within vs between site differences
dist_df <- melt(as.matrix(unifrac_dist))
dist_df = dist_df[dist_df[1] != dist_df[2],] # remove dif between same subsite
dist_df$site1 <- substr(dist_df[[1]], 1, 1)
dist_df$site2 <- substr(dist_df[[2]], 1, 1)
dist_df$comparison <- "between"
dist_df[dist_df$site1 == dist_df$site2, ]$comparison <- "within" # dif within sites

# plot differences
ggplot(dist_df, aes(x = comparison, y = value, colour = comparison)) + geom_boxplot() + theme_bw() + stat_compare_means(label.y=0.1) +
  stat_compare_means(label="p.signif",ref.group = ".all.",
                     hide.ns = TRUE,comparisons = list(c("between", "within")))
# ggsave("unifrac_dist.png", plot = last_plot(), height = 5, width = 6)





# Distance decay ----------------------------------------------------------



## distances vs geographic distance - distance decay
# unifrac vs geographic distances
unifrac_dist <- UniFrac(otu_all_fc, weighted = TRUE) %>% as.matrix() # for FC normalised
dist_df <- melt(as.matrix(unifrac_dist))
dist_df = dist_df[dist_df[1] != dist_df[2],] # remove dif between same subsite
colnames(dist_df) <- c("sample1", "sample2", "unifrac")

# similar approach, combine glacier distance with unifrac distances
s_data <- sample_data(otu_all)
s_data <- as.matrix(s_data[,c(1,19)]) # get glacier dist info

g_dist <- merge(dist_df, s_data, by.x = "sample1", by.y = "Sample") # get glacier distance for sample1
colnames(g_dist)[4] <- "sample1_glacier_dist" # rename variable for sample1
g_dist <- merge(g_dist, s_data, by.x = "sample2", by.y = "Sample") # get glacier distance for sample1
colnames(g_dist)[5] <- "sample2_glacier_dist" # rename variable for sample1

g_dist$dist_df <- abs(as.numeric(g_dist$sample1_glacier_dist) - as.numeric(g_dist$sample2_glacier_dist)) # the distance between samples

fit1 <- lm(unifrac ~ dist_df, data = g_dist)
summary(fit1)
ggplot(g_dist, aes(x = dist_df, y = unifrac)) + geom_point() +
  stat_smooth(method = "lm", col = "red") + xlab("Geographic distance between sites (km)")  + 
  labs(title = paste("Adj R2 = ",signif(summary(fit1)$adj.r.squared, 5),
                     ", p =",signif(summary(fit1)$coef[2,4], 5))) + theme_bw() + ylab("Weighted Unifrac Distance")
#ggsave("unifrac_dist_geographic_dist.png", plot = last_plot(), height = 4, width = 4)



# Bray curtis vs geographic
# distance matrix
set.seed(4102023)
bray_dist <- vegdist(t(otu_table(otu_all_fc)), method ="bray") # bray-curtis dissimilarity
bray_dist_df <- melt(as.matrix(bray_dist))
bray_dist_df = bray_dist_df[bray_dist_df[1] != bray_dist_df[2],] # remove dif between same subsite
colnames(bray_dist_df) <- c("sample1", "sample2", "bray_curtis")
# image(as.matrix(bray_dist)) 

g_dist_b <- merge(bray_dist_df, s_data, by.x = "sample1", by.y = "Sample") # get glacier distance for sample1
colnames(g_dist_b)[4] <- "sample1_glacier_dist" # rename variable for sample1
g_dist_b <- merge(g_dist_b, s_data, by.x = "sample2", by.y = "Sample") # get glacier distance for sample1
colnames(g_dist_b)[5] <- "sample2_glacier_dist" # rename variable for sample1

g_dist_b$dist_df <- abs(as.numeric(g_dist_b$sample1_glacier_dist) - as.numeric(g_dist_b$sample2_glacier_dist)) # the distance between samples

fit2 <- lm(bray_curtis ~ dist_df, data = g_dist_b)
summary(fit2)
ggplot(g_dist_b, aes(x = dist_df, y = bray_curtis)) + geom_point() +
  stat_smooth(method = "lm", col = "red") + xlab("Geographic distance between sites (km)")  + 
  labs(title = paste("Adj R2 = ",signif(summary(fit2)$adj.r.squared, 5),
                     ", p =",signif(summary(fit2)$coef[2,4], 5))) + theme_bw() + ylab("Bray-curtis")
#ggsave("braycurtis_dist_geographic.png", plot = last_plot(), height = 4, width = 4)



# unifrac and bray-curtis vs geographic dist
# reformat to merge
g_dist$distance_m <- "unifrac"
g_dist_b$distance_m <- "bray_curtis"
colnames(g_dist)[3] <- "distance"
colnames(g_dist_b)[3] <- "distance"

all_dist <- rbind(g_dist, g_dist_b) # merge distance metrics
stat_values <- data.frame(distance_m = c("unifrac", "bray_curtis"),
                          p_val = c(signif(summary(fit1)$coef[2,4], 3), signif(summary(fit2)$coef[2,4], 3)),
                          r2 = c(round(summary(fit1)$adj.r.squared, 2), round(summary(fit2)$adj.r.squared, 2)),
                          slope = c(round(fit1$coefficients[2], 3), round(fit2$coefficients[2], 3)) ,
                          intercept = c(round(fit1$coefficients[1],2), round(fit2$coefficients[1], 2)))

#all_dist[all_dist$distance_m == "unifrac", ]$distance_m <- "weighted_unifrac" # rename
all_dist[all_dist$distance_m == "unifrac", ]$distance_m <- "Weighted UniFrac" # rename
all_dist[all_dist$distance_m == "bray_curtis", ]$distance_m <- "Bray-Curtis" # rename

dist_decay <- ggplot(all_dist, aes(x = dist_df, y = distance, colour = distance_m, group = distance_m)) + geom_point(alpha = 0.3) +
  stat_smooth(method = "lm") + xlab("Geographic distance between sites (km)")  + 
  annotate("text", y = c(0.1, 0.16), x = 3, label = paste0("slope==", stat_values$slope, "~~r^2==", stat_values$r2), parse = TRUE,
           colour = c("salmon", "#3A9AB2")) +
  scale_colour_manual(values = c("#3A9AB2", "salmon")) + 
  theme_bw() + ylab("Dissimilarity") + labs(colour = NULL)
#ggsave("distance_decay.png", plot = last_plot(), height = 4, width = 5.5)



# check for sig difference
summary(lm(distance ~ dist_df + distance_m + dist_df:distance_m, all_dist))
# the estimate is the slope, use the interaction slope
# this is saying that unifrac has a slope 0.004 lower than bray, and is significant







# mantel test for comparison
unifrac_dist_mat <- UniFrac(otu_all_fc, weighted = TRUE) %>% as.matrix() # for FC normalised
bray_dist_mat <- vegdist(t(otu_table(otu_all_fc)), method ="bray") %>% as.matrix()

# geographic distance mat
geo_wide <- g_dist %>% select(sample1, sample2, dist_df) %>%
  pivot_wider(names_from = sample2, values_from = dist_df) %>% as.matrix
rownames(geo_wide) <- geo_wide[,1]
rownames(geo_wide) <- geo_wide[,1]
geo_wide <- geo_wide[,-1]
geo_wide[is.na(geo_wide)] <- 0 # assign 0 comparisons
# order the same as distance mat
geo_wide <- geo_wide[rownames(unifrac_dist_mat), colnames(unifrac_dist_mat)]


# mantel test
mantel_unifrac <- mantel(unifrac_dist_mat, geo_wide, method = "spearman")
mantel_bray <- mantel(bray_dist_mat, geo_wide, method = "spearman")

colnames(unifrac_dist_mat)
colnames(geo_wide)



ggarrange(rel_abund, ggarrange(alpha_diversity, dist_decay, labels = c("B", "C"), align = "h"), nrow = 2, labels = "A", legend = "right", heights = c(1.5,1))
# ggsave("fig3.png", plot = last_plot(), height = 10, width = 10, dpi = 900)





## Top only
unifrac_dist_t <- UniFrac(subset_samples(otu_all_fc, Depth == "Top"), weighted = TRUE) %>% as.matrix()
dist_df_t <- melt(as.matrix(unifrac_dist_t))
dist_df_t = dist_df_t[dist_df_t[1] != dist_df_t[2],] # remove dif between same subsite
colnames(dist_df_t) <- c("sample1", "sample2", "unifrac")

g_dist_t <- merge(dist_df_t, s_data, by.x = "sample1", by.y = "Sample") # get glacier distance for sample1
colnames(g_dist_t)[4] <- "sample1_glacier_dist" # rename variable for sample1
g_dist_t <- merge(g_dist_t, s_data, by.x = "sample2", by.y = "Sample") # get glacier distance for sample1
colnames(g_dist_t)[5] <- "sample2_glacier_dist" # rename variable for sample1

g_dist_t$dist_df <- abs(as.numeric(g_dist_t$sample1_glacier_dist) - as.numeric(g_dist_t$sample2_glacier_dist)) # the distance between samples

fit1 <- lm(unifrac ~ dist_df, data = g_dist_t)
summary(fit1)


# mantel test
# geographic distance mat
geo_wide <- g_dist_t %>% select(sample1, sample2, dist_df) %>%
  pivot_wider(names_from = sample2, values_from = dist_df) %>% as.matrix
rownames(geo_wide) <- geo_wide[,1]
rownames(geo_wide) <- geo_wide[,1]
geo_wide <- geo_wide[,-1]
geo_wide[is.na(geo_wide)] <- 0 # assign 0 comparisons
# order the same as distance mat
geo_wide <- geo_wide[rownames(unifrac_dist_t), colnames(unifrac_dist_t)]

mantel_top <- mantel(unifrac_dist_t, geo_wide, method = "spearman")





## Lower only
unifrac_dist_l <- UniFrac(subset_samples(otu_all_fc, Depth == "Lower"), weighted = TRUE) %>% as.matrix()
dist_df_l <- melt(as.matrix(unifrac_dist_l))
dist_df_l = dist_df_l[dist_df_l[1] != dist_df_l[2],] # remove dif between same subsite
colnames(dist_df_l) <- c("sample1", "sample2", "unifrac")

g_dist_l <- merge(dist_df_l, s_data, by.x = "sample1", by.y = "Sample") # get glacier distance for sample1
colnames(g_dist_l)[4] <- "sample1_glacier_dist" # rename variable for sample1
g_dist_l <- merge(g_dist_l, s_data, by.x = "sample2", by.y = "Sample") # get glacier distance for sample1
colnames(g_dist_l)[5] <- "sample2_glacier_dist" # rename variable for sample1

g_dist_l$dist_df <- abs(as.numeric(g_dist_l$sample1_glacier_dist) - as.numeric(g_dist_l$sample2_glacier_dist)) # the distance between samples

fit2 <- lm(unifrac ~ dist_df, data = g_dist_l)
summary(fit2)

# mantel test
# geographic distance mat
geo_wide <- g_dist_l %>% select(sample1, sample2, dist_df) %>%
  pivot_wider(names_from = sample2, values_from = dist_df) %>% as.matrix
rownames(geo_wide) <- geo_wide[,1]
rownames(geo_wide) <- geo_wide[,1]
geo_wide <- geo_wide[,-1]
geo_wide[is.na(geo_wide)] <- 0 # assign 0 comparisons
# order the same as distance mat
geo_wide <- geo_wide[rownames(unifrac_dist_l), colnames(unifrac_dist_l)]

mantel_low <- mantel(unifrac_dist_l, geo_wide, method = "spearman")




# reformat to merge depths
g_dist_l$depth <- "Lower"
g_dist_t$depth <- "Top"

all_dist_d <- rbind(g_dist_t, g_dist_l) # merge distance metrics for depths

stat_values_d <- data.frame(depth = c("Top", "Lower"),
                            p_val = c(signif(summary(fit1)$coef[2,4], 3), signif(summary(fit2)$coef[2,4], 3)),
                            r2 = c(round(summary(fit1)$adj.r.squared, 2), round(summary(fit2)$adj.r.squared, 2)),
                            slope = c(round(fit1$coefficients[2], 3), round(fit2$coefficients[2], 3)) ,
                            intercept = c(round(fit1$coefficients[1],2), round(fit2$coefficients[1], 2)))

all_dist_d$depth <- factor(all_dist_d$depth, levels = c("Top", "Lower"))
dist_decay_depth <- ggplot(all_dist_d, aes(x = dist_df, y = unifrac, colour = depth, group = depth)) + geom_point() +
  stat_smooth(method = "lm") + xlab("Geographic distance between sites (km)")  + 
  annotate("text", y = c(0.1, 0.14), x = 3, label = paste0("slope==", stat_values_d$slope, "~~r^2==", stat_values_d$r2), parse = T,
           colour = c("darkgrey", "black")) +
  scale_colour_manual(values = c("darkgrey", "black")) + labs(colour = NULL) +
  theme_bw() + ylab("UniFrac Dissimilarity")
#ggsave("distance_decay_depth.png", plot = last_plot(), height = 4, width = 5, dpi = 900)


summary(lm(unifrac ~ dist_df + depth + dist_df:depth, all_dist_d))
# the estimate is the slope, this is saying that lower has a slope 0.001 higher than top, not significant








# NMDS --------------------------------------------------------------------


# NMDS
nmds <- metaMDS(unifrac_dist) # new NMDS object
nmds$stress # 0.184
plot(nmds)
metadata$Sample <- row.names(metadata)
scores_nmds <- scores(nmds) %>% as_tibble(rownames = "Sample") %>% inner_join(., metadata, by = "Sample") #combined metadata


# sig env with taxa envfit

# add env vars
otu_df <- as.data.frame(otu_table(otu_all_fc))
otu_df <- as.data.frame(otu_table(otu_all_fc))
env_vars <- metadata
env_vars <- subset(env_vars, env_vars$Sample %in% colnames(otu_df)) # subset env variables for only samples analysed
env_vars <- subset(env_vars, env_vars$Sample %in% scores_nmds$Sample) # subset env variables for only samples analysed
env_vars$Sample <- factor(env_vars$Sample, levels = colnames(otu_df)) #reorder
env_vars <- env_vars[order(env_vars$Sample), ] #reorder
env_trans <- decostand(env_vars[c(4:19)], "standardize") # transform values by standardising
env_vars <- cbind(env_vars[c(1:3)], env_trans)

# phyla abundance wide
filt_data <- psmelt(otu_all_fc) 
phyla_data <- filt_data[,c(2,3,5,6, 26)] # select only phylum data
phyla_sums <- phyla_data %>% group_by(Sample, Phylum) %>% summarise(Abundance_sum = sum(Abundance)) # total phyla abundance in sample
phyla <- pivot_wider(phyla_sums, names_from = Phylum, values_from = Abundance_sum) # create wide dataframe

phyla_trans <- cbind(phyla[,1], as.data.frame(clr(phyla[,-1]))) # CLR transformation, equivalent to log transformation for FC counts

metadata_p <- merge(env_vars, phyla_trans, by = "Sample") # merge with metadata 

#reorder
metadata_p$Sample <- factor(metadata_p$Sample, levels = colnames(otu_df))
metadata_p <- metadata_p[order(metadata_p$Sample), ]

env_ord <- envfit(nmds, metadata_p, permutations = 999, na.rm = T)
env_ord # see sig vectors
plot(env_ord)
en_coord_cont <- as.data.frame(scores(env_ord, "vectors")) * ordiArrowMul(env_ord)
en_coord_cat <- as.data.frame(scores(env_ord, "factors")) * ordiArrowMul(env_ord)

# sig env vectors
pvals <- as.data.frame(env_ord$vectors[4])
pvals$padj <- p.adjust(pvals$pvals, method = "fdr")
sig_env <- rownames(pvals[pvals$padj < 0.05, ])
en_coord_cont_sig <- en_coord_cont %>% filter(rownames(en_coord_cont) %in% c(sig_env)) # select sig env vectors and main env variables

scores_nmds$Depth <- factor(scores_nmds$Depth, levels = c("Top", "Lower")) # order depth

# rename taxa
rownames(en_coord_cont_sig)[rownames(en_coord_cont_sig) == "Proteobacteria"] <- "Pseudomonadota"
rownames(en_coord_cont_sig)[rownames(en_coord_cont_sig) == "Actinobacteriota"] <- "Actinomycetota"
rownames(en_coord_cont_sig)[rownames(en_coord_cont_sig) == "Chloroflexi"] <- "Chloroflexota"
rownames(en_coord_cont_sig)[rownames(en_coord_cont_sig) == "Cyanobacteria"] <- "Cyanobacteriota"


# rename vars
row.names(en_coord_cont_sig)[row.names(en_coord_cont_sig) == "Water_content"] <- "'Water Content'"
row.names(en_coord_cont_sig)[row.names(en_coord_cont_sig) == "Glacier_dist"] <- "'Glacier Dist'"
row.names(en_coord_cont_sig)[row.names(en_coord_cont_sig) == "CFU_count"] <- "'CFU Count'"

row.names(en_coord_cont_sig)[row.names(en_coord_cont_sig) == "Na"] <- "Na^'+'"
row.names(en_coord_cont_sig)[row.names(en_coord_cont_sig) == "NH4"] <- "NH[4]^'+'"
row.names(en_coord_cont_sig)[row.names(en_coord_cont_sig) == "K"] <- "K^'+'"
row.names(en_coord_cont_sig)[row.names(en_coord_cont_sig) == "Cl"] <- "Cl^'-'"
row.names(en_coord_cat)[row.names(en_coord_cat) == "DepthLower"] <- "Lower"
row.names(en_coord_cat)[row.names(en_coord_cat) == "DepthTop"] <- "Top"

nmds_all <- ggplot(data = scores_nmds, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(data = scores_nmds, aes(colour = Site, shape = Depth), size = 4)  + 
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), data = en_coord_cont_sig, colour = "grey") +
  ggrepel::geom_text_repel(data = en_coord_cont_sig, aes(x=NMDS1, y=NMDS2, label = row.names(en_coord_cont_sig)), direction = "both", parse = T, segment.size = 0.1, max.overlaps = 10) +
  scale_shape_manual(values=c(19,1)) + theme_bw() + scale_color_manual(labels = c("Glacier", "Site 1", "Site 1 Vegetated", "Site 2", "Site 3", "Site 4", "Site 5"), values = wes_palette("Zissou1", 7, type = "continuous"))
#ggsave("nmds_unifrac_env_sig_phlya_FC.png", plot = last_plot(), height = 5, width = 6)





# PERMANOVA ---------------------------------------------------------------



# PERMANOVA + pairwise
# standardised ordered env vars
adonis2(unifrac_dist ~ Site, metadata_p, permutations = 999) # r2 0.269, p < 0.001
adonis2(unifrac_dist ~ Depth, metadata_p, permutations = 999) # r2 0.113, p< 0.001
adonis2(unifrac_dist ~ Site * Depth, metadata_p, permutations = 999) # r2 0.453, p < 0.001, F = 2.2377

# install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
pw.adonis <- pairwise.adonis2(unifrac_dist ~ Site * Depth, metadata_p, permutations = 999) # does not require method as unifrac dist already calc
summary(pw.adonis)

# p value for each comparison
pairwise_results <- data.frame()

# Loop through each comparison in pw.adonis and extract the p-values
for (comparison in names(pw.adonis)[-1]) { # iterate through each comparison (avoiding the parent call)
  print(comparison)
  p_value <- pw.adonis[[comparison]]$`Pr(>F)`[1]  # get p-value
  comp1 <- gsub("_vs_.*", "", comparison) # get site comparison names
  comp2 <- gsub(".*_vs_", "", comparison)
  pairwise <- data.frame(Site1 = min(comp1, comp2), Site2 = max(comp1,comp2), P_Value = p_value) # associate site comparison (ordered)
  pairwise_results <- rbind(pairwise_results, pairwise) # store the result
}

# p adjust
pairwise_results$P_adj <- p.adjust(pairwise_results$P_Value, method = "fdr") 

# Create the heatmap with labels and color scale that splits at P_value = 0.05
ggplot(pairwise_results, aes(x = as.factor(Site1), y = as.factor(Site2))) +
  geom_tile(aes(fill = ifelse(P_adj <= 0.05, "p-adj<0.05", "NS")), color = "black") +  # Apply conditional fill
  scale_fill_manual(values = c("p-adj<0.05" = "pink", "NS" = "white")) +  # Manually set fill colors
  geom_text(aes(label = round(P_adj, 3)), color = "black", size = 3) +  # Add P_value labels
  theme_minimal() +
  labs(title = "Pairwise PERMANOVA",
    x = "", y = "", fill = "P-adj Value") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_discrete(labels = c("Glacier", "Site 1", "Site 1 \nVegetated", "Site 2", "Site 3", "Site 4")) +
  scale_y_discrete(labels = c("Site 1", "Site 1 \nVegetated", "Site 2", "Site 3", "Site 4", "Site 5")) 
#ggsave(filename = "pairwise_permanova_pvalues.png", plot = last_plot(), height = 4, width = 5, dpi = 900)






# Env Correlation ---------------------------------------------------------



# Spearman correlations
# non parametric (not normal data) and use raw
env_vars <- metadata[-c(3,5,6,20,21)] # remove depth, air temp, light and soil fluxes, remove ions

# hist(env_trans$FC_count)

# spearman correlation
meta_cor <- env_vars %>% filter(!is.na(Soil_temp) & !is.na(Grain_size)) %>% 
  dplyr::select(where(is.numeric)) %>% cor(method = "spearman") 

meta_cor_p <- env_vars %>% filter(!is.na(Soil_temp) & !is.na(Grain_size)) %>% 
  dplyr::select(where(is.numeric)) %>% cor_pmat(method = "spearman") # pvalue from correlation

# extract upper traingle (not duplicates) and apply p adjust to total tests (not column by column)
meta_cor_p_adjust <- meta_cor_p
meta_cor_p_upper <- meta_cor_p[upper.tri(meta_cor_p)] # upper triangle
meta_cor_p_adj_upper <- p.adjust(meta_cor_p_upper, method = "fdr") # padj of upper triangle
meta_cor_p_adjust[upper.tri(meta_cor_p_adjust)] <- meta_cor_p_adj_upper # reassign back into dataframe
meta_cor_p_adjust[lower.tri(meta_cor_p_adjust)] <- t(meta_cor_p_adjust)[lower.tri(meta_cor_p_adjust)] # mirror in lower triangle
diag(meta_cor_p_adjust) <- 0 # set diagonal


# plot
meta_cor_plot <- ggcorrplot(meta_cor,
           hc.order = F,
           type = "upper",
           outline.color = "white",
           p.mat = meta_cor_p_adjust,
           insig = "blank",
           col = c("#3A9AB2", "white", "#F11B00"),
           legend.title = "Spearman Correlation")

# add pvalues 
labs.function = function(x){
  case_when(x > 0.055 ~ "",
            x <= 0.052 & x >= 0.01 ~ "*",
            x < 0.01 & x >= 0.001 ~ "**",
            x < 0.001 ~ "***")}

# Get asterisk matrix based on p-values
meta_cor_p_labs <- as.data.frame(meta_cor_p_adjust) %>% mutate_all(labs.function) # determine asterisks from values
meta_cor_p_labs$Var1 <- as.factor(rownames(meta_cor_p_labs)) # add col before long transformation
meta_cor_p_labs <- melt(meta_cor_p_labs, id.vars = "Var1", variable.name = "Var2", value.name = "value")
colnames(meta_cor_p_labs)[2] <- "Var2"
meta_cor_p_labs$in.df <- ifelse(is.na(match(paste0(meta_cor_p_labs$Var1, meta_cor_p_labs$Var2),
                                            paste0(meta_cor_plot[["data"]]$Var1, meta_cor_plot[["data"]]$Var2))),
                                "No", "Yes")
meta_cor_p_labs <- dplyr::select(filter(meta_cor_p_labs, in.df == "Yes"), -in.df)

meta_cor_plot +
  geom_text(aes(x = meta_cor_p_labs$Var1,
                y = meta_cor_p_labs$Var2),
            label = meta_cor_p_labs$value,
            #nudge_y = 0.25,
            size = 4) + 
  scale_y_discrete(labels = c("Glacier_dist" = "Glacier Dist", "Soil_temp" = "Soil Temp", "FC_count" = "FC Count", "CFU_count" = "CFU Count", "Grain_size" = "Grain Size", "Water_content" = "Water Content",
                                                  "NH4" = expression(NH[4]^"+"), "Na" = expression(Na^'+'), "Mg" = expression(Mg^'2+'), "K" = expression(K^'+'), 
                                                  "Ca" = expression(Ca^'2+') , "Cl" = expression(Cl^'-'))) +
  scale_x_discrete(labels = c("Glacier_dist" = "Glacier Dist", "Soil_temp" = "Soil Temp", "FC_count" = "FC Count", "CFU_count" = "CFU Count", "Grain_size" = "Grain Size", "Water_content" = "Water Content",
                                                                                                                                         "NH4" = expression(NH[4]^"+"), "Na" = expression(Na^'+'), "Mg" = expression(Mg^'2+'), "K" = expression(K^'+'), 
                                                                                                                                         "Ca" = expression(Ca^'2+') , "Cl" = expression(Cl^'-')))
#ggsave("env_spearman_cor_sig.png", plot = last_plot(), height = 5, width = 7, dpi = 900)






# Env PCA -----------------------------------------------------------------



# ENV PCA


# format
env_trans <- metadata
env_trans[env_trans$Depth == "Top", ]$Depth <- 1 # make depth numeric
env_trans[env_trans$Depth == "Lower", ]$Depth <- 0
env_trans$Depth <- as.numeric(env_trans$Depth)
env_trans <- env_trans[-c(5,6, 20, 21)] # remove light and air temp, gas flux

#env_trans <- env_trans[env_trans$Sample %in% row.names(sample_data(otu_all)), ] # select only analysed samples

hist(env_trans$FC_count)
hist(env_trans$CFU_count)
env_trans$FC_count <- log(env_trans$FC_count) # log transform
env_trans$CFU_count <- log(env_trans$CFU_count + 1) # log transform

env_trans <- decostand(env_trans[c(4:17)], "standardize") # standardise

# check if data is linear (using standardised env data)
lm_fc <- lm(FC_count ~ Soil_temp + Grain_size + Water_content + pH + EC , env_trans)
summary(lm_fc) # linear regression check

par(mfrow=c(2,2))
plot(lm_fc) # check least square error
par(mfrow=c(1,1))

# not normal or linear data use PCA (NMDS for complex datasets)
fc_pca <- princomp(na.omit(env_trans))
plot(fc_pca )


# contribution of each variable to the 1 and 2 dimensions
fviz_cos2(fc_pca, choice = "var", axes = 1:2)
fviz_pca_var(fc_pca, col.var = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE)


# create dataframe of pca values
fc_pca_df <- as.data.frame(fc_pca$loadings[,1:2])
fc_pca_vars <- get_pca_var(fc_pca) # access cos2 variable
# is how well a variable is represented on a principle component, tutorial: https://biosakshat.github.io/pca.html
fc_pca_df$cos2 <- rowSums(fc_pca_vars$cos2[,1:2]) # order of co2 info is same as dataframe

# edit variable names for plot
fc_loadings <- as.data.frame(fc_pca$loadings[,1:2])
rownames(fc_loadings)[rownames(fc_loadings) == "Soil_temp"] <- "'Soil Temp'"
rownames(fc_loadings)[rownames(fc_loadings) == "Grain_size"] <- "'Grain Size'"
rownames(fc_loadings)[rownames(fc_loadings) == "Water_content"] <- "'Water Content'"
rownames(fc_loadings)[rownames(fc_loadings) == "Glacier_dist"] <- "'Glacier Dist'"
rownames(fc_loadings)[rownames(fc_loadings) == "FC_count"] <- "'FC Count'"
rownames(fc_loadings)[rownames(fc_loadings) == "CFU_count"] <- "'CFU Count'"

rownames(fc_loadings)[rownames(fc_loadings) == "Na"] <- "Na^'+'"
rownames(fc_loadings)[rownames(fc_loadings) == "NH4"] <- "NH[4]^'+'"
rownames(fc_loadings)[rownames(fc_loadings) == "Mg"] <- "Mg^'2+'"
rownames(fc_loadings)[rownames(fc_loadings) == "K"] <- "K^'+'"
rownames(fc_loadings)[rownames(fc_loadings) == "Ca"] <- "Ca^'2+'"
rownames(fc_loadings)[rownames(fc_loadings) == "Cl"] <- "Cl^'-'"



# in this case cos2 shows how well both components represent a variable
env_pca <- ggplot(fc_pca_df, aes(x=Comp.1, y=Comp.2)) + geom_point(aes(colour = cos2)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey") + geom_vline(xintercept = 0, linetype = "dashed", colour = "grey") +
  geom_segment(aes(x = 0, y = 0, xend = Comp.1, yend = Comp.2, colour = cos2)) +
  ggrepel::geom_text_repel(aes(x=Comp.1, y=Comp.2, label = row.names(fc_loadings)), parse = T) +
  theme_bw() + xlab("PCA1 (34.7%)") + ylab("PCA2 (14.7%)") + scale_colour_gradient2(low = "#00AFBB", mid = "#E7B800", high = "#FC4E07", midpoint = 0.5)
#ggsave("ENV PCA.png", plot = last_plot(), height = 4, width = 5)


# add point measurements
# create dataframe of pca values
fc_pca_scores <- as.data.frame(fc_pca$scores )
fc_pca_scores$Sample <- rownames(fc_pca_scores)
sample_data <- metadata[,c(1:3)]
fc_pca_scores_meta <- inner_join(fc_pca_scores, sample_data, by = "Sample") # add pca scores to data
fc_pca_scores_meta$Depth <- factor(fc_pca_scores_meta$Depth, levels = c("Top", "Lower")) # order depths

env_pca_point <- ggplot(data = fc_pca_scores_meta, aes(x=Comp.1, y=Comp.2)) +
  geom_point(data = fc_pca_scores_meta, aes(colour = Site, shape = Depth)) +
  geom_segment(aes(x = 0, y = 0, xend = Comp.1*10, yend = Comp.2*10), data = fc_pca_df, colour = "grey") +
  ggrepel::geom_text_repel(data = fc_pca_df, aes(x=Comp.1*10, y=Comp.2*10, label = row.names(as.data.frame(fc_pca$loadings[,1:2]))), max.overlaps = 20, segment.size = 0.1) +
  scale_shape_manual(values=c(16, 1)) + theme_bw() + scale_color_manual(labels = c("Glacier", "Site 1", "Site 1 Vegetated", "Site 2", "Site 3", "Site 4", "Site 5"), values = wes_palette("Zissou1", 7, type = "continuous")) +
  xlab("PCA1 (34.7%)") + ylab("PCA2 (14.7%)")
#ggsave("ENV PCA2.png", plot = last_plot(), height = 4, width = 5)


env_pca <- ggplot(fc_pca_df, aes(x=Comp.1*10, y=Comp.2*10)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey") + geom_vline(xintercept = 0, linetype = "dashed", colour = "grey") +
  geom_segment(aes(x = 0, y = 0, xend = Comp.1 *10, yend = Comp.2*10, colour = cos2)) +
  ggrepel::geom_text_repel(aes(x=Comp.1*10.5, y=Comp.2*10.5, label = row.names(fc_loadings)), parse = T) +
  theme_bw() + xlab("PCA1 (34.7%)") + ylab("PCA2 (14.7%)") + 
  scale_colour_gradient2(low = "#00AFBB", mid = "#E7B800", high = "#FC4E07", midpoint = 0.5) + 
  geom_point(data = fc_pca_scores_meta, aes(x=Comp.1, y=Comp.2, shape = Depth)) + scale_shape_manual(values=c(16, 1))






# Env Raw -----------------------------------------------------------------



# Environmental Data Plots

# plot ch4 data
ch4 <- ggplot(metadata[metadata$Site %in% c(1, "1_VEG",2,3,4),], aes(x = Site, y = flux_CH4)) + geom_boxplot() + geom_point() + theme_bw() +
  geom_hline(yintercept = 0, linetype = "dashed") + labs(colour = "Site") + xlab("") + ylab(expression(CH[4]~nmol~m^{-2}~s^{-1})) + 
  theme(legend.position="none", axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  scale_x_discrete(labels = c("Site 1", "Site 1 \nVegetated", "Site 2", "Site 3", "Site 4")) 
#ggsave("ch4.png", plot = last_plot(), height = 5, width = 5)
min(na.omit(metadata$flux_CH4)) # -0.8409409
max(na.omit(metadata$flux_CH4)) # 0.1501967

cor.test(metadata$flux_CH4, metadata$Glacier_dist, method = "spearman")
cor.test(metadata$flux_CO2, metadata$Glacier_dist, method = "spearman")


co2 <- ggplot(metadata[metadata$Site %in% c(1, "1_VEG",2,3,4),], aes(x = Site, y = flux_CO2)) + geom_boxplot() + geom_point() + theme_bw() +
  geom_hline(yintercept = 0, linetype = "dashed") + labs(colour = "Site") + xlab("") + ylab(expression(CO[2]~umol~m^{-2}~s^{-1})) + 
  theme(legend.position="none", axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  scale_x_discrete(labels = c("Site 1", "Site 1 \nVegetated", "Site 2", "Site 3", "Site 4")) 
#ggsave("co2.png", plot = last_plot(), height = 5, width = 5)
min(na.omit(metadata$flux_CO2)) # -0.02200417
max(na.omit(metadata$flux_CO2)) # 0.8290738




# Metadata (all samples)
metadata_long <- gather(metadata, key = env_variable, value = value, -c(1:3, 5, 19, 20, 21)) # only use samples for otu analysis and keep glacier distance

# order variables
metadata_long$env_variable <- factor(metadata_long$env_variable, levels = c("Soil_temp", "Air_temp", "Water_content", "pH", "EC", "Grain_size", "FC_count", "CFU_count",
                                                                            "Ca", "K", "Mg", "Na", "NH4", "Cl"))
# rename
env_labels <- c(`Soil_temp` = "Soil Temp (°C)", `Air_temp` = "Air Temp (°C)", `Water_content` = "Water content (%)", `pH` = "pH",
                          `EC` = "EC (dS/m)", `Grain_size` = "Grain size (um)", `FC_count` = "FC count (cells/g)", `CFU_count` = "CFU count (CFU/g)",
                `Ca` = "Ca2+", `K` = "K+", `Mg` = "Mg2+", `Na` = "Na+", `NH4` = "NH4+", `Cl` = "Cl-")

options(scipen = -3)
metadata_long$Depth <- factor(metadata_long$Depth, levels = c("Top", "Lower"))

metadata_long[metadata_long$env_variable != "Air_temp",] %>% 
  ggplot(aes(x = Site, y= value, fill = Site, alpha = Depth, alpha = 0.5)) + geom_boxplot(outliers = F)  + 
  facet_wrap(~env_variable, nrow = 5, scales = "free_y", labeller = as_labeller(env_labels)) + 
  theme_bw() + ylab("") + xlab("") + scale_fill_discrete(labels = c("Glacier", "Site 1", "Site 1 Vegetated", "Site 2", "Site 3", "Site 4", "Site 5") ) +
  scale_x_discrete(labels = c("Glacier", "Site 1", "Site 1 \nVegetated", "Site 2", "Site 3", "Site 4", "Site 5"))  + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
#ggsave("env_vars.png", plot = last_plot(), height = 10, width = 10, dpi = 900)

options(scipen = 100) # reset back





# plot cfu vs fc
max(metadata$FC_count)
max(na.omit(metadata$CFU_count))

min(metadata$FC_count)
order(na.omit(metadata$CFU_count))


# linear regression
# format
env_trans <- metadata
env_trans$FC_count <- log(env_trans$FC_count) # log transform
env_trans$CFU_count <- log(env_trans$CFU_count + 1) # log transform

fit1 <- lm( CFU_count ~ FC_count , data = env_trans) # log transformed
summary(fit1)


options(scipen = -3) 

# log transformed
ggplot(env_trans, aes(x = FC_count, y = CFU_count, colour = Site))  +
  geom_point() + geom_abline(linetype = "dashed") +
  theme_bw() + xlab(Flow~Cytometry~Counts~(log(cells~g^{-1}))) + ylab(expression(Culture-based~Counts~(log(CFU~g^{-1})))) + 
  scale_colour_discrete(labels = c("Glacier", "Site 1", "Site 1 Vegetated", "Site 2", "Site 3", "Site 4", "Site 5")) +
  xlim(0,15)
#ggsave("cell_abundance.png", plot = last_plot(), height = 5, width = 7, dpi = 900)

# t test for non normal data
wilcox.test(env_trans$FC_count, 
            env_trans$CFU_count, alternative = "greater", paired = T)



# test cell abundance between depths
wilcox.test(env_trans[env_trans$Depth == "Top", ]$FC_count, 
            env_trans[env_trans$Depth == "Lower", ]$FC_count, paired = F)

wilcox.test(env_trans[env_trans$Depth == "Top", ]$CFU_count, 
            env_trans[env_trans$Depth == "Lower", ]$CFU_count, paired = F)





# Paired depth ------------------------------------------------------------




# DEPTH COMPARISON 

# env vars
shannon_trans <- decostand(shannon[-c(1:3)], "standardize") # standardise
shannon_trans <- cbind(shannon[c(1:3)], shannon_trans)
shannon_trans <- shannon_trans[-c(5, 6, 19, 20, 21)] # remove light, air temp and glacier dist, fluxes
shannon_trans$Subsite <- gsub("_[TL]",  "", shannon_trans$Sample)

# split layers
top_data <- shannon_trans[shannon_trans$Depth == "Top", ]
low_data <- shannon_trans[shannon_trans$Depth == "Lower", ]

# remove samples without two depths
top_data <- top_data[top_data$Subsite %in% low_data$Subsite, ]
low_data <- low_data[low_data$Subsite %in% top_data$Subsite, ]
low_data$Subsite == top_data$Subsite # check true order

top_data_num <- top_data %>% dplyr::select(where(is.numeric)) # subset only numeric vars
low_data_num <- low_data %>% dplyr::select(where(is.numeric))

# get difference dataframe
dif_data <- top_data_num - low_data_num
dif_data <- cbind(dif_data, top_data[c(2,18)]) # add site and subsite

# format long
dif_data_long <- gather(dif_data, key = "Env_Variable", value = "Top_Lower_Diff", -c(15,16))

# test difference
wilcox_p_all <- data.frame()
for(env_var in colnames(top_data_num)){
  wilcox_p <- wilcox.test(top_data_num[[env_var]], low_data_num[[env_var]], paired = T)$p.value # non parametric paired test
  wilcox_p_df <- data.frame(Env_Variable = env_var, p_val = wilcox_p)
  wilcox_p_all <- rbind(wilcox_p_all, wilcox_p_df)
}

# adjust pvalues
wilcox_p_all$p_adj <- p.adjust(wilcox_p_all$p_val, method = "fdr")
dif_data_long <- merge(dif_data_long, wilcox_p_all, by = "Env_Variable")

# plot
# order by difference values
dif_data_mean <- na.omit(dif_data_long) %>% group_by(Env_Variable) %>% summarise(mean_dif = mean(Top_Lower_Diff), mean_pval = mean(p_adj))
dif_data_long$Env_Variable <- factor(dif_data_long$Env_Variable, levels = dif_data_mean[order(dif_data_mean$mean_dif, decreasing = F),]$Env_Variable)
dif_data_long$Env_Variable <- factor(dif_data_long$Env_Variable, levels = c("Ca", "Mg", "K", "NH4", "Cl", "Na", "EC", "Water_content","Grain_size", "CFU_count", "FC_count", "Shannon",  "pH", "Soil_temp"))


ggplot(dif_data_long, aes(x = Top_Lower_Diff, y = Env_Variable, fill = p_adj)) + geom_violin(scale = "width", draw_quantiles = 0.5) + theme_bw() + 
   annotate("text", y = wilcox_p_all$Env_Variable, x = 3.6, label = paste("p.adj = ", signif(wilcox_p_all$p_adj, digits = 2), sep = ""), size = 2.5) + 
   geom_vline(xintercept =  0, linetype = "dashed") + scale_fill_gradient(high = "white", low = "#F11B00") + ylab("") + theme(legend.position = "none") + xlim(-2.5, 4) + xlab("Top - Lower (LFC)")
# ggsave("top_lower_diff.png", plot = last_plot(), height = 4, width = 5)


ggplot(data = dif_data_long) + geom_violin(aes(x = Top_Lower_Diff, y = Env_Variable), scale = "width", draw_quantiles = 0.5) + 
  geom_violin(data = dif_data_long[dif_data_long$Env_Variable == "pH" | dif_data_long$Env_Variable == "Soil_temp",], aes(x = Top_Lower_Diff, y = Env_Variable), fill = "red", colour = "red", alpha = 0.2, scale = "width", draw_quantiles = 0.5) + 
  theme_bw() + 
  annotate("text", y = wilcox_p_all$Env_Variable, x = 3.6, label = paste("p.adj = ", signif(wilcox_p_all$p_adj, digits = 2), sep = ""), size = 2.5) + 
  ylab("") + theme(legend.position = "none") + xlab("Top - Lower (LFC)") + 
  scale_y_discrete(labels = c("Soil_temp" = "Soil Temp", "FC_count" = "FC Count", "CFU_count" = "CFU Count", "Grain_size" = "Grain Size", "Water_content" = "Water Content",
                              "NH4" = expression(NH[4]^"+"), "Na" = expression(Na^'+'), "Mg" = expression(Mg^'2+'), "K" = expression(K^'+'), 
                              "Ca" = expression(Ca^'2+') , "Cl" = expression(Cl^'-'))) + geom_vline(xintercept =  0, linetype = "dashed") + xlim(-2.5, 4) 



dif_data_long$colour <- "none"
dif_data_long[dif_data_long$p_adj < 0.05, ]$colour <- "red"

depth_env <- ggplot(data = dif_data_long, aes(x = Top_Lower_Diff, y = Env_Variable, colour = colour, fill = colour)) + 
  geom_violin(scale = "width", draw_quantiles = 0.5, alpha = 0.2) +
  #geom_point(aes(x = Top_Lower_Diff, y = Env_Variable), alpha = 0.2) +
  scale_color_manual(values = c("black", "red")) +
  scale_fill_manual(values = c(NA, "red")) +
  theme_bw() + 
  ylab("") + theme(legend.position = "none") + xlab("Difference (z-score)") + ggtitle("Low                                                            Top") +
  scale_y_discrete(labels = c("Soil_temp" = "Soil Temp", "FC_count" = "FC Count", "CFU_count" = "CFU Count", "Grain_size" = "Grain Size", "Water_content" = "Water Content",
                              "NH4" = expression(NH[4]^"+"), "Na" = expression(Na^'+'), "Mg" = expression(Mg^'2+'), "K" = expression(K^'+'), 
                              "Ca" = expression(Ca^'2+') , "Cl" = expression(Cl^'-'))) + geom_vline(xintercept =  0, linetype = "dashed") 



ggarrange(ggarrange(env_pca, depth_env, nrow = 2, labels = c("A", "B"), align = "hv"), ggarrange(ch4, co2, nrow = 2, labels = c("C", "D"), align = "hv"), ncol = 2, widths = c(6,4))

ggarrange(ggarrange(env_pca, depth_env, nrow = 2, labels = c("A", "B"), align = "v"), ggarrange(ch4, co2, nrow = 2, labels = c("C", "D"), align = "hv"), ncol = 2, widths = c(6,4))
# ggsave("fig2.png", plot = last_plot(), height = 9, width = 10, dpi = 900)
# ggsave("fig2.svg", plot = last_plot(), height = 9, width = 10, dpi = 900, device = svglite::svglite)




# Depth DA ----------------------------------------------------------------




## DA analysis - ANCOMB ANALYSIS

# using phyloseq object
otu_all
sample_data(otu_all)

# using the Holm-Bonferroni method
set.seed(123)


# run DA
# output_interact = ancombc2(data = otu_all, tax_level = "Genus",
#                            fix_formula = "Depth + Site", rand_formula = NULL,
#                            p_adj_method = "holm", pseudo_sens = TRUE,
#                            prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05,
#                            group = "Depth", struc_zero = TRUE, neg_lb = F,
#                            alpha = 0.05, n_cl = 2, verbose = TRUE
#                            )
# 
# write.csv(output_interact$res, "depth_ancomb_FC.csv", row.names = F)

# load results
output_interact <- read.csv("depth_ancomb_FC.csv", header = T)

# select fams with sig difference
sig <- output_interact[output_interact$diff_DepthTop == TRUE,]
sig$taxon <-  gsub("Genus:", "", output_interact[output_interact$diff_DepthTop == TRUE,]$taxon) # format family names

sig$taxon <-  gsub("Family:", "", sig$taxon) # format family names
# removes uncultured groups

taxa <- as.data.frame(tax_table(otu_all), row.names = F) # load full taxonomies as df

taxa_sig <- taxa[taxa$Genus %in% sig$taxon, ] # get full taxa for sig fams
taxa_sig <- taxa_sig[-c(7)] # remove species col
taxa_sig <- unique(taxa_sig) # remove duplicates
taxa_sig <- taxa_sig[taxa_sig$Genus != "uncultured", ] # remove uncultured
taxa_sig <- merge(taxa_sig, sig, by.x = "Genus", by.y = "taxon") # merge full taxa with DA results

# order by lfc
taxa_sig$Genus <- factor(taxa_sig$Genus, levels = taxa_sig[order(taxa_sig$lfc_DepthTop, decreasing = F), ]$Genus)

# rename phyla
taxa_sig[taxa_sig$Phylum == "Chloroflexi", ]$Phylum <- "Chloroflexota"
taxa_sig[taxa_sig$Phylum == "Cyanobacteria", ]$Phylum <- "Cyanobacteriota"
taxa_sig[taxa_sig$Phylum == "Firmicutes", ]$Phylum <- "Bacillota"
taxa_sig[taxa_sig$Phylum == "Proteobacteria", ]$Phylum <- "Pseudomonadota"
taxa_sig[taxa_sig$Phylum == "Actinobacteriota", ]$Phylum <- "Actinomycetota"


da_plot <- ggplot(taxa_sig, aes(x = lfc_DepthTop, y = Genus, fill = Phylum)) + geom_bar(stat = "identity") + theme_bw() + 
  facet_grid( vars(Phylum), scales = "free_y", space = "free_y") + ylab("Differentially Abundant Genera") +
  theme(strip.background = element_blank(), strip.text.y = element_blank()) + xlab("LFC") + ggtitle("Low                                                                    Top")
# ggsave("DA_ancomb_sig_genus_FC.png", plot = last_plot(), height = 5, width = 5)

# # plot individual taxa to check
# da_data <- filt_data[filt_data$Genus %in% taxa_sig$Genus,]
# da_data_mean <- da_data %>% group_by(Site, Depth, OTU, Family, Genus) %>% summarise(mean_abundance = mean(Abundance))
# ggplot(da_data_mean[da_data_mean$Genus == "Pir4_lineage", ], aes(x = Site, y = mean_abundance, fill = Genus, colour = Genus)) + facet_wrap(~Depth) +
#   geom_bar(stat = "identity", position="dodge") + theme_bw() + scale_y_sqrt() + ylab("OTU Abundance (FC)")

ggarrange(nmds_all, da_plot, labels = c("A", "B"), ncol = 1, heights = c(2,3))
# ggsave("fig4.png", plot = last_plot(), height = 13, width = 8, dpi = 900)





# WGCNA clusters ----------------------------------------------------------




## WGCNA ANALYSIS
# https://www.polarmicrobes.org/weighted-gene-correlation-network-analysis-wgcna-applied-to-microbial-communities/

# clr transformation
otu_norm <- otu_all %>% microbiome::transform(transform = "clr") # relative data with clr transformation

# only core taxa
otu_core_taxa <- filter_taxa(otu_all, function(x) sum(x >= 1) > (0.05*length(x)), TRUE) # selects otus in at least 5% of the sample
# 2083 OTUS
otu_core <- otu_norm
otu_table(otu_core) <- otu_table(otu_norm)[row.names(otu_table(otu_core_taxa))] # select these otus from the transformed data above
ntaxa(otu_core) #  2083





# percentage abundance
sum(otu_table(otu_all))
otu_core_rel <- otu_table(otu_all)[row.names(otu_table(otu_core_taxa))]

core_data <- filt_data[filt_data$OTU %in% row.names(otu_table(otu_core_taxa)),]

# relative abundances
otu_rel <- transform_sample_counts(otu_all, function(x) x / sum(x))
core_rel <- otu_table(otu_rel)[row.names(otu_table(otu_core_taxa))]

# sum otu per sample
otu_sum <- rowSums(otu_table(otu_rel))

sum(otu_sum)
median(otu_sum)
sd(otu_sum)

core_sum <- rowSums(core_rel)

sum(core_sum) / sum(otu_sum) * 100
median(core_sum)
sd(core_sum)





# analysis
otu_core_t <- as.data.frame(t(as.data.frame(otu_table(otu_core)))) # transpose

names(otu_core_t) # OTUS
rownames(otu_core_t) # samples

gsg <- goodSamplesGenes(otu_core_t, verbose = 3)
gsg$allOK # TRUE indicates that all otu had sufficient abundance (not less 0)

# cluster samples based on OTU
sample_tree <- hclust(dist(otu_core_t), method = "average")

sizeGrWindow(12,9)
par(cex = 0.6);
par(mar = c(0,4,2,0))

# plot clustered samples
plot(sample_tree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)

# add env data
sample_data(otu_core)[,c(4,6:19)]
traitColour <- numbers2colors(sample_data(otu_core)[,c(4,6:19)], signed = FALSE);

# plot env values for each sample
plotDendroAndColors(sample_tree, traitColour,
                    groupLabels = names(sample_data(otu_core)[,c(4,6:19)]),
                    main = "Sample dendrogram and trait heatmap")

# allow multithreading
options(stringsAsFactors = FALSE);
enableWGCNAThreads()

# select thresholds
powers <- c(c(1:10), seq(from = 11, to=30, by=1))

# network topology analysis
sft <- pickSoftThreshold(otu_core_t, powerVector = powers, verbose = 5, networkType = "signed")

# plot results
sizeGrWindow(9, 5)
par(mfrow = c(1,2))
cex1 = 0.9

plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");

abline(h=0.8,col="red")

plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")


# calculate adjacency
softPower = 6 # based on above plot
adjacency = adjacency(otu_core_t, power = softPower, type = "signed")

# adjacency matrix into a Topological Overlap Matrix (TOM) and calculate dissimilarity
TOM = TOMsimilarity(adjacency, TOMType = "signed")
dissTOM = 1-TOM

TaxaTree = hclust(as.dist(dissTOM), method = "average")

# clustering of all subsetted OTUs based on the TOM dissimilarity index
sizeGrWindow(12,9)
plot(TaxaTree, xlab="", sub="", main = "Taxa clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04);

# module size
minModuleSize = 20

# Module identification using dynamic tree cut
dynamicMods = cutreeDynamic(dendro = TaxaTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)

# colours
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors) # this shows how many modules

# plot modules
sizeGrWindow(8,6)
plotDendroAndColors(TaxaTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Taxa dendrogram and module colors")

# quantify co-expression of the entire modules 
MEList = moduleEigengenes(otu_core_t, colors = dynamicColors)
MEs = MEList$eigengenes

MEDiss = 1-cor(MEs) # dissimilarity 

# cluster modules
METree = hclust(as.dist(MEDiss), method = "average") 

# plot module clusters
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")

# check to merge modules
MEDissThres = 0.30 # cutoff 0.3

abline(h=MEDissThres, col = "red") # plot cutoff

# merge
merge = mergeCloseModules(otu_core_t, dynamicColors, cutHeight = MEDissThres, verbose = 3)
mergedColors = merge$colors
mergedMEs = merge$newMEs # eigengenes of new merged modules

# playing
mergedMEDiss = 1-cor(mergedMEs) # dissimilarity 
# cluster modules
mergedMETree = hclust(as.dist(mergedMEDiss), method = "average") 
# plot module clusters
sizeGrWindow(7, 6)
plot(mergedMETree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")

# replot
sizeGrWindow(12, 9)

plotDendroAndColors(TaxaTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

# rename
moduleColors = mergedColors

# create labels
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder)-1
MEs = mergedMEs


# save results
#write.csv(MEs, "MEs.csv", row.names = T)
#write.csv(moduleColors, "moduleColors.csv", row.names = F)




# WGCNA heatmaps ----------------------------------------------------------






# reload
# MEs2 <- read.csv("MEs.csv", header = T, row.names = 1) # load saved Module eigenvectors
# moduleColors <- read.csv("moduleColors.csv", header = T)
# moduleColors <- moduleColors$x
# 
# otu_norm <- otu_all %>% microbiome::transform(transform = "clr") 
# otu_core_taxa <- filter_taxa(otu_all, function(x) sum(x >= 1) > (0.05*length(x)), TRUE) 
# otu_core <- otu_norm
# otu_table(otu_core) <- otu_table(otu_norm)[row.names(otu_table(otu_core_taxa))] 
# otu_core_t <- as.data.frame(t(as.data.frame(otu_table(otu_core)))) 
# 
# # load merged me tree





# associate modules with env variables
nTaxa = ncol(otu_core_t)
nSamples = nrow(otu_core_t)

# first edit sample data for figure
sample_data_edit <- sample_data(otu_core)
sample_data_edit[sample_data_edit$Depth == "Top", ]$Depth <- 1 # make depth binary
sample_data_edit[sample_data_edit$Depth == "Lower", ]$Depth <- 0
sample_data_edit$Depth <- as.numeric(sample_data_edit$Depth) # make site 1 veg numeric
sample_data_edit_f <- sample_data_edit[,-c(1,2,5,6,12:17,20,21)] # reassign sample metadata, removing ions

MEs0 = moduleEigengenes(otu_core_t, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, sample_data_edit_f, use = "p", method = "spearman")
#moduleTraitCor = cor(MEs, sample_data_edit_f, use = "p") # pearson
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)


# plot modules with env data
sizeGrWindow(10,6)
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));

labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(sample_data_edit_f),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

# replot cuter
moduleTraitCor_df <- as.data.frame(moduleTraitCor)
moduleTraitCor_df$module <- row.names(moduleTraitCor_df)
moduleTraitCor_long <- gather(moduleTraitCor_df, key = Variable, value = value, -10)

textMatrix_df <- as.data.frame(textMatrix) # add p values
textMatrix_df <- gather(textMatrix_df, key = var, value = sig)
textMatrix_df$pvalue <- gsub(".*\\(", "", textMatrix_df$sig)
textMatrix_df$pvalue <- as.numeric(gsub("\\)", "", textMatrix_df$pvalue))
moduleTraitCor_long <- cbind(moduleTraitCor_long, textMatrix_df)


# add sig
labs.function = function(x){
  case_when(x > 0.055 ~ "",
            x <= 0.055 & x >= 0.01 ~ "*",
            x < 0.01 & x >= 0.001 ~ "**",
            x < 0.001 ~ "***")}

# adjust pvalue for multiple testing
moduleTraitCor_long$pvalue_adj <- p.adjust(moduleTraitCor_long$pvalue, method = "BH")
moduleTraitCor_long$symbol <- labs.function(moduleTraitCor_long$pvalue_adj)

#write.csv(moduleTraitCor_long, "moduleTraitCor_long.csv", row.names = F)
#moduleTraitCor_long <- read.csv("moduleTraitCor_long.csv", header = T) # load saved results for wgcna heatmap

# plot
ggplot(moduleTraitCor_long, aes(x = Variable, y = gsub("ME", "", module), fill = value)) + geom_tile() + scale_fill_gradient2(
  low = "blue",
  high = "red",
  mid = "white",
  midpoint = 0,
  limit = c(-1,1)) + theme_bw() +
  theme(axis.text.x = element_text(angle=90)) +
  labs(y = "", fill="corr") +
  geom_text(aes(x = moduleTraitCor_long$Variable, 
                y = gsub("ME", "", moduleTraitCor_long$module)), 
            label = moduleTraitCor_long$symbol, 
            #nudge_y = 0.25, 
            size = 5)





# add tree to heatmap
# ggtree(mergedMETree) + geom_treescale() +
#   geom_tiplab(align = TRUE, linesize = 0, size = 3)
module_tree <- ggtree(mergedMETree, branch.length = "none") 

# get module order in tree
module_order <- na.omit(arrange(module_tree$data, by = module_tree$data["y"])["label"])
# rearrange tile plot to match
moduleTraitCor_long$module <- gsub("ME", "", moduleTraitCor_long$module) # remove ME string
module_order$label <- gsub("ME", "", module_order$label)
moduleTraitCor_long$module <- factor(moduleTraitCor_long$module, levels = module_order$label)

module_heatmap <- ggplot(moduleTraitCor_long, aes(x = Variable, y = module, fill = value)) + geom_tile() + scale_fill_gradient2(
  low = "#3A9AB2", 
  high = "#F11B00",
  mid = "white",
  midpoint = 0,
  limit = c(-1,1)) + theme_bw() +
  theme(axis.text.x = element_text(angle=90)) +
  labs(y = "", fill="Spearman \nCorrelation") +
  geom_text(aes(x = moduleTraitCor_long$Variable, 
                y = moduleTraitCor_long$module), 
            label = moduleTraitCor_long$symbol, 
            #nudge_y = 0.25, 
            size = 5) + xlab("Total") + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  scale_y_discrete(labels = c("Module 11", "Module 10", "Module 9", "Module 8", "Module 7", "Module 6", 
                              "Module 5", "Module 4", "Module 3", "Module 2", "Module 1")) + 
  scale_x_discrete(labels = c("CFU Count","Depth", "EC", "FC Count", "Glacier Dist", "Grain Size", "pH", "Soil Temp", "Water Content"))





# arrange plots together
ggarrange(module_tree, module_heatmap, widths = c(1,10), labels = "A", align = "h") # add tree plot
# ggsave("wgcna_heatmap_spearman.png", plot = last_plot(), height = 5, width = 7)






# depth analysis for flux of wgcna
# top
mergedMEs_T <- mergedMEs[grep("_T", rownames(mergedMEs)),]

sample_data_edit_f2 <- sample_data_edit[,c(20,21)] # reassign sample metadata, removing ions
sample_data_edit_f2 <- sample_data_edit_f2[grep("_T", rownames(sample_data_edit_f2))]

rownames(sample_data_edit_f2) == rownames(mergedMEs_T)

moduleTraitCor2 = cor(mergedMEs_T, sample_data_edit_f2, use = "p", method = "spearman")
moduleTraitPvalue2 = corPvalueStudent(moduleTraitCor2, length(mergedMEs_T$MEturquoise))

textMatrix2 = paste(signif(moduleTraitCor2, 2), "\n(",
                    signif(moduleTraitPvalue2, 1), ")", sep = "");

labeledHeatmap(Matrix = moduleTraitCor2,
               xLabels = names(sample_data_edit_f2),
               yLabels = names(mergedMEs_T),
               ySymbols = names(mergedMEs_T),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix2,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

# replot cuter
moduleTraitCor2_df <- as.data.frame(moduleTraitCor2)
moduleTraitCor2_df$module <- row.names(moduleTraitCor2_df)
moduleTraitCor2_long <- gather(moduleTraitCor2_df, key = Variable, value = value, -3)

textMatrix2_df <- as.data.frame(textMatrix2) # add p values
textMatrix2_df <- gather(textMatrix2_df, key = var, value = sig)
textMatrix2_df$pvalue <- gsub(".*\\(", "", textMatrix2_df$sig)
textMatrix2_df$pvalue <- as.numeric(gsub("\\)", "", textMatrix2_df$pvalue))
moduleTraitCor2_long <- cbind(moduleTraitCor2_long, textMatrix2_df)

# multiple testing p adjustment
# adjust pvalue for multiple testing
moduleTraitCor2_long$pvalue_adj <- p.adjust(moduleTraitCor2_long$pvalue, method = "BH")
moduleTraitCor2_long$symbol <- labs.function(moduleTraitCor2_long$pvalue_adj)

# reorder to match phylogeny order
moduleTraitCor2_long$module <- gsub("ME", "", moduleTraitCor2_long$module) # remove ME string
moduleTraitCor2_long$module <- factor(moduleTraitCor2_long$module, levels = module_order$label)

# plot
ggplot(moduleTraitCor2_long, aes(x = Variable, y = module, fill = value)) + geom_tile() + scale_fill_gradient2(
  low = "#3A9AB2", 
  high = "#F11B00",
  mid = "white",
  midpoint = 0,
  limit = c(-1,1)) + theme_bw() +
  theme(axis.text.x = element_text(angle=90)) +
  labs(y = "", fill="Spearman Correlation") +
  geom_text(aes(x = moduleTraitCor2_long$Variable, 
                y = moduleTraitCor2_long$module), 
            label = moduleTraitCor2_long$symbol, 
            #nudge_y = 0.25, 
            size = 5) + xlab("Top") + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  scale_y_discrete(labels = c("Module 11", "Module 10", "Module 9", "Module 8", "Module 7", "Module 6", 
                              "Module 5", "Module 4", "Module 3", "Module 2", "Module 1"))

wgcna_top <- ggplot(moduleTraitCor2_long, aes(x = Variable, y = module, fill = value)) + geom_tile() + scale_fill_gradient2(
  low = "#3A9AB2", 
  high = "#F11B00",
  mid = "white",
  midpoint = 0,
  limit = c(-1,1)) + theme_bw() +
  theme(axis.text.x = element_text(angle=90)) +
  labs(y = "", fill="Spearman \nCorrelation") +
  geom_text(aes(x = moduleTraitCor2_long$Variable, 
                y = moduleTraitCor2_long$module), 
            label = moduleTraitCor2_long$symbol, 
            #nudge_y = 0.25, 
            size = 5) + xlab("Top") + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), 
                                            axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  scale_x_discrete(labels = c(expression(CH[4]~Flux), expression(CO[2]~Flux)))
# ggsave("wgcna_flux_heatmap_top.png", plot = last_plot(), height = 5, width = 3)






# depth analysis for flux of wgcna
# Lower
mergedMEs_L <- mergedMEs[grep("_L", rownames(mergedMEs)),]

sample_data_edit_f3 <- sample_data_edit[,c(20,21)] # reassign sample metadata, removing ions
sample_data_edit_f3 <- sample_data_edit_f3[grep("_L", rownames(sample_data_edit_f3))]

rownames(sample_data_edit_f3) == rownames(mergedMEs_L)

moduleTraitCor3 = cor(mergedMEs_L, sample_data_edit_f3, use = "p")
moduleTraitPvalue3 = corPvalueStudent(moduleTraitCor3, length(mergedMEs_L$MEturquoise))

textMatrix3 = paste(signif(moduleTraitCor3, 2), "\n(",
                    signif(moduleTraitPvalue3, 1), ")", sep = "");

labeledHeatmap(Matrix = moduleTraitCor3,
               xLabels = names(sample_data_edit_f3),
               yLabels = names(mergedMEs_L),
               ySymbols = names(mergedMEs_L),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix3,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

# replot cuter
moduleTraitCor3_df <- as.data.frame(moduleTraitCor3)
moduleTraitCor3_df$module <- row.names(moduleTraitCor3_df)
moduleTraitCor3_long <- gather(moduleTraitCor3_df, key = Variable, value = value, -3)

textMatrix3_df <- as.data.frame(textMatrix3) # add p values
textMatrix3_df <- gather(textMatrix3_df, key = var, value = sig)
textMatrix3_df$pvalue <- gsub(".*\\(", "", textMatrix3_df$sig)
textMatrix3_df$pvalue <- as.numeric(gsub("\\)", "", textMatrix3_df$pvalue))
moduleTraitCor3_long <- cbind(moduleTraitCor3_long, textMatrix3_df)

# multiple testing p adjustment
# adjust pvalue for multiple testing
moduleTraitCor3_long$pvalue_adj <- p.adjust(moduleTraitCor3_long$pvalue, method = "BH")
moduleTraitCor3_long$symbol <- labs.function(moduleTraitCor3_long$pvalue_adj)

# reorder to match phylogeny order
moduleTraitCor3_long$module <- gsub("ME", "", moduleTraitCor3_long$module) # remove ME string
moduleTraitCor3_long$module <- factor(moduleTraitCor3_long$module, levels = module_order$label)

# plot
ggplot(moduleTraitCor3_long, aes(x = Variable, y = module, fill = value)) + geom_tile() + scale_fill_gradient2(
  low = "#3A9AB2", 
  high = "#F11B00",
  mid = "white",
  midpoint = 0,
  limit = c(-1,1)) + theme_bw() +
  theme(axis.text.x = element_text(angle=90)) +
  labs(y = "", fill="Spearman Correlation") +
  geom_text(aes(x = moduleTraitCor3_long$Variable, 
                y = moduleTraitCor3_long$module), 
            label = moduleTraitCor3_long$symbol, 
            #nudge_y = 0.25, 
            size = 5) + xlab("Lower") + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  scale_y_discrete(labels = c("Module 11", "Module 10", "Module 9", "Module 8", "Module 7", "Module 6", 
                              "Module 5", "Module 4", "Module 3", "Module 2", "Module 1"))
# ggsave("wgcna_flux_heatmap_low.png", plot = last_plot(), height = 5, width = 3)

wgcna_low <- ggplot(moduleTraitCor3_long, aes(x = Variable, y = module, fill = value)) + geom_tile() + scale_fill_gradient2(
  low = "#3A9AB2", 
  high = "#F11B00",
  mid = "white",
  midpoint = 0,
  limit = c(-1,1)) + theme_bw() +
  theme(axis.text.x = element_text(angle=90)) +
  labs(y = "", fill="Spearman \nCorrelation") +
  geom_text(aes(x = moduleTraitCor3_long$Variable, 
                y = moduleTraitCor3_long$module), 
            label = moduleTraitCor3_long$symbol, 
            #nudge_y = 0.25, 
            size = 5) + xlab("Lower") + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), 
                                              axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  scale_x_discrete(labels = c(expression(CH[4]~Flux), expression(CO[2]~Flux)))

wgcna_heatmap_all <- ggarrange(module_tree, module_heatmap, wgcna_top, wgcna_low, widths = c(2,11,3,3), nrow = 1, labels = c("A"), common.legend = T, legend = "right", align = "h") # add tree plot that exactly fits nodes
# ggsave("fig4A.png", plot = last_plot(), height = 5, width = 8, dpi = 900)
# ggsave("fig4A.svg", plot = last_plot(), height = 5, width = 10, dpi = 900, device = svglite::svglite)






# WGCNA taxa enrichment ---------------------------------------------------




## Enrichment analysis of WGCNA modules

# using modules above
module_taxa <- as.data.frame(tax_table(otu_core)[row.names(tax_table(otu_core)) %in% names(otu_core_t), ])
module_taxa$Taxon <- row.names(module_taxa)
module_taxa$moduleColor <- moduleColors # add module group (this is order based!)
# write.csv(module_taxa, "WGCNA_module_taxa.csv", row.names = F)

# create column with full name so to not confound unique genera (e.g. uncultured)
module_taxa$name <- paste(module_taxa$Kingdom, module_taxa$Phylum, module_taxa$Class, module_taxa$Order, module_taxa$Family, sep = "_")

otu_core_counts <- as.data.frame(otu_table(otu_all_fc)[row.names(tax_table(otu_all_fc)) %in% names(otu_core_t), ])
otu_core_counts$sum <- rowSums(otu_core_counts) # sum fc count across all samples (maybe ave would be better here?)
otu_core_counts$Taxon <- rownames(otu_core_counts) # create barcode column
module_taxa <- merge(module_taxa, otu_core_counts[, c(50, 51)]) # associate modules + taxa with barcodes

ModCountsM <- module_taxa %>% count(name, Family, name="m") # family total counts
ModCountsN <- ModCountsM %>% mutate(n = nTaxa - m) %>% select (-m) # difference from total count


# get enrichment of taxa in each module
# get abundance of each otu in module and sum
ModFCounts <- data.frame()
for(module in unique(module_taxa$moduleColor)){
  
  DFModuleAnnot = (module_taxa) %>% filter(moduleColor==module) # subset module
  k = nrow(DFModuleAnnot) # total number of OTUs in Module

  n0 = length(unique(DFModuleAnnot$name)) # Count of uniq family in module (aka no. tests performed)
  
  ModCountsSum <- DFModuleAnnot %>% group_by(name, Phylum, Family) %>% summarise(family_mod_sum = sum(sum)) # get family sum of otu abundance in module
  
  ModCountsQ <- DFModuleAnnot %>% count(name, Family, name = "q") # each family count in module
  ModCounts <- left_join(ModCountsQ, ModCountsM,  by = "name")
  ModCounts <- right_join(ModCountsN, ModCounts, by = "name")
  ModCounts <- mutate(ModCounts,k=k)
  ModCounts <- mutate(ModCounts, p = phyper(q, m, n, k, lower.tail = FALSE)) # density distribution (sampling without replacement)
  Mod1Counts <- mutate(Module = module, ModCounts, pAdj = p.adjust(p, n0, method = "BH"))
  Mod1Counts <- merge(Mod1Counts, ModCountsSum[,-3], by = "name" ) # add module otu fam sum abundance to output, has to remove fam col which is duplicated
  
  ModFCounts <- rbind(ModFCounts, Mod1Counts)
}


# filter for significance
ModFCounts_filt <- ModFCounts %>% filter(q>0) %>% filter(pAdj<0.01)

# filter by percentage abundance of module
module_abundance <- ModFCounts %>% group_by(Module) %>% summarise(module_sum = sum(family_mod_sum)) # sum of all families sum in each module (before filtering)
ModFCounts_filt <- merge(ModFCounts_filt, module_abundance, by = "Module") # add to enrichment data
ModFCounts_filt <- ModFCounts_filt[ModFCounts_filt$family_mod_sum > (0.05 * ModFCounts_filt$module_sum), ] # families with abundance greater than 5% of the module total abundance

# order module
ModFCounts_filt$Module <- factor(ModFCounts_filt$Module , levels = gsub("ME", "", rev(module_order$label)))
ModFCounts_filt$name <- factor(ModFCounts_filt$name, levels = unique(ModFCounts_filt[order(ModFCounts_filt$Phylum), ]$name))

# rename taxa
ModFCounts_filt[ModFCounts_filt$Phylum == "Actinobacteriota", ]$Phylum <- "Actinomycetota"
ModFCounts_filt[ModFCounts_filt$Phylum == "Cyanobacteria", ]$Phylum <- "Cyanobacteriota"
ModFCounts_filt[ModFCounts_filt$Phylum == "Proteobacteria", ]$Phylum <- "Pseudomonadota"
#ModFCounts_filt[ModFCounts_filt$Phylum == "Actinobacteriota", ]$Phylum <- "Actinomycetota"

ModFCounts_filt[ModFCounts_filt$Family == "Hydrogenophilaceae", ]$Family <- "Thiobacillacaeae" # Thiobacillus genus from Hydrogenophilaceae to Thiobacillacaeae
ModFCounts_filt[ModFCounts_filt$Family == "uncultured", ]$Family <- "Gaiellales_uncultured" # change uncultured name

# rename modules as numbers
mod_col_num <- data.frame(mod = module_order, mod_num = c(11:1))
colnames(mod_col_num)[1] <- "Module"
ModFCounts_filt <- merge(ModFCounts_filt, mod_col_num, by = "Module")

# plot enriched fam abundance
ggplot(ModFCounts_filt, aes(y = Family, x =  family_mod_sum )) + geom_bar(aes(fill = Phylum), alpha = 0.5, stat="identity", position="stack", colour = "black") + 
  theme_bw() + ylab("Family") + xlab("Total OTU Read Count (FC)") + facet_grid( vars(mod_num), scales = "free_y", space = "free_y") + scale_x_sqrt()
# ggsave("module_phyla_enriched3.png", plot = last_plot(), height = 8, width = 6)

# plot enriched fam percentage
wgcna_taxa <- ggplot(ModFCounts_filt, aes(y = Family, x =  family_mod_sum/module_sum * 100)) + geom_bar(aes(fill = Phylum), alpha = 0.5, stat="identity", position="stack", colour = "black") + 
  theme_bw() + ylab("Enriched Families in Modules") + xlab("Proportion of Module Abundance (%)") + facet_grid( vars(mod_num), scales = "free_y", space = "free_y") + 
  theme(legend.position = "bottom", legend.justification = c(0, 0)) + guides(fill = guide_legend(nrow = 4, title.position = "top"))
# ggsave("module_phyla_enriched2.png", plot = last_plot(), height = 8, width = 6, dpi = 900)
# write.csv(ModFCounts_filt, "WGCNA_enriched_fam.csv", row.names = F)






# WGCNA cluster env ------------------------------------------------


# # investigate taxa module with one env factor
# env_factor <- sample_data_edit_f[,9] # glacier distance
# 
# # get correlation
# modNames = substring(names(MEs), 3)
# TaxaModuleMembership = as.data.frame(cor(otu_core_t, MEs, use = "p")); # pearson correlation
# MMPvalue = as.data.frame(corPvalueStudent(as.matrix(TaxaModuleMembership), nSamples)); # student asymptotic p-value for correlation
# names(TaxaModuleMembership) = paste("MM", modNames, sep="");
# names(MMPvalue) = paste("p.MM", modNames, sep="");
# 
# TaxaTraitSignificance = as.data.frame(cor(otu_core_t, env_factor, use = "p"));
# GSPvalue = as.data.frame(corPvalueStudent(as.matrix(TaxaTraitSignificance), nSamples));
# names(TaxaTraitSignificance) = paste("GS.env_var");
# names(GSPvalue) = paste("p.GS.env_var");
# 
# # get all module taxa data
# TaxaInfo0 = as.data.frame(tax_table(otu_core)[row.names(tax_table(otu_core)) %in% names(otu_core_t), ])
# TaxaInfo0$Taxon <- row.names(TaxaInfo0)
# TaxaInfo0$moduleColor <- moduleColors # add module group (this is order based!)
# 
# # add env correlation data (order is correct)
# TaxaInfo0$TaxaTraitSignificance <- TaxaTraitSignificance$GS.env_var
# TaxaInfo0$GSPvalue <- GSPvalue$p.GS.env_var
# 
# # add module membership after reformatting
# TaxaModuleMembership_clean <- TaxaModuleMembership
# TaxaModuleMembership_clean$Taxon <- row.names(TaxaModuleMembership)
# TaxaModuleMembership_clean <- gather(TaxaModuleMembership_clean, key = moduleColor, value = module_membership, -12)
# TaxaModuleMembership_clean$moduleColor <- gsub("MM", "", TaxaModuleMembership_clean$moduleColor)
# 
# TaxaInfo <- merge(TaxaInfo0, TaxaModuleMembership_clean, by = c("Taxon", "moduleColor")) # merges only membership for matching otu module
# 
# # order based on module and then env var
# TaxaOrder = order(TaxaInfo$moduleColor, -abs(TaxaInfo$TaxaTraitSignificance));
# TaxaInfo = TaxaInfo[TaxaOrder, ]
# 
# # note: 
# # TaxaInfo$TaxaTraitSignificance = correlation to env variable
# # TaxaInfo$GSPvalue = pval for the correlation
# # TaxaInfo$MM.colour = module membership (each otu belongs to one module but has values for all)
# # TaxaInfo$p.MM.colour = significance of module membership
# 
# # identify top modules for env correlation
# cor(MEs, env_factor, use = "p")[order(-abs(cor(MEs, env_factor, use = "p"))),] # check most correlated modules
# 
# ggplot(TaxaInfo, aes(x = module_membership, y = TaxaTraitSignificance, colour = moduleColor)) + geom_point() +
#   ylab("Correlation") + xlab("Module Membership")  + theme_bw() + facet_wrap(~moduleColor, ncol = 3)





# module abundance plot with env sig

# for all modules
top_mod_env_all <- data.frame()
for(mod in c(unique(TaxaInfo$moduleColor))){
  top_mod_env <- TaxaInfo[TaxaInfo$moduleColor == mod, ]
  top_mod_env <- top_mod_env[order(top_mod_env$module_membership, decreasing = T),][c(1:10),] # most significantly correlated
  top_mod_env_all <- rbind(top_mod_env_all, top_mod_env)
}

# get module average abundance pattern from top members
env_taxa_fc <- otu_all_fc
otu_table(env_taxa_fc) <- otu_table(env_taxa_fc)[row.names(otu_table(env_taxa_fc)) %in% top_mod_env_all$Taxon,] # select these otus from the transformed data above
ntaxa(env_taxa_fc) 
env_taxa_fc <- psmelt(env_taxa_fc) # df
env_taxa_fc <- merge(env_taxa_fc, top_mod_env_all, by.x = "OTU", by.y = "Taxon") # combine module data back with abundance
env_taxa_fc_ave <- env_taxa_fc %>% group_by(Glacier_dist, moduleColor, Genus.x, OTU) %>% summarise(Mean_abundance = mean(Abundance))


# resassign modules as numbers using mod_col_num
mod_col_num <- data.frame(mod = module_order, mod_num = c(11:1))
colnames(mod_col_num)[1] <- "moduleColor" # rename to match this dataframe
env_taxa_fc_ave <- merge(env_taxa_fc_ave, mod_col_num, by = "moduleColor") # add to dataframe
env_taxa_fc_ave$mod_num <- paste("Module ", env_taxa_fc_ave$mod_num, sep = "") #add  module prefix
env_taxa_fc_ave$mod_num <- factor(env_taxa_fc_ave$mod_num, levels = c("Module 1", "Module 2", "Module 3", "Module 4", 
                                                                      "Module 5", "Module 6", "Module 7", "Module 8", 
                                                                      "Module 9", "Module 10", "Module 11"))



# create categories for sig correlations
env_taxa_fc_ave$Sig <- "Other"
env_taxa_fc_ave[env_taxa_fc_ave$mod_num %in% c("Module 1", "Module 3", "Module 4", "Module 5",
                                               "Module 6", "Module 9", "Module 10"), ]$Sig <- "Glacier Distance"
env_taxa_fc_ave[env_taxa_fc_ave$mod_num %in% c("Module 4", "Module 10"), ]$Sig <- "Glacier Distance & CO2"
env_taxa_fc_ave$Sig <- factor(env_taxa_fc_ave$Sig, levels = c("Other", "Glacier Distance", "Glacier Distance & CO2")) # order


cluster_glacier <- ggplot(env_taxa_fc_ave, aes(x = Glacier_dist, y = Mean_abundance, group = OTU, colour = Sig, fill = Sig)) + 
  geom_area(alpha = 0.1, position = "dodge") + theme_bw() + 
  scale_fill_manual(labels = c("Other", "Glacier Distance", expression(Glacier~Distance~and~CO[2])), values = c("darkgrey","#85B7B9", "#3A9AB2" )) + 
  scale_colour_manual(labels = c("Other", "Glacier Distance", expression(Glacier~Distance~and~CO[2])), values = c("grey","#85B7B9", "#3A9AB2" )) +
  facet_wrap(~mod_num, ncol = 2, scales = "free_y") + ylab("Mean OTU Absolute Abundance") + xlab("Glacier Distance (km)") + 
  labs(fill = "Signifcant Correlation", colour = "Signifcant Correlation") + theme(legend.position = "bottom", legend.justification = c(0, 0)) +
  guides(fill = guide_legend(nrow = 4, title.position = "top"))



enriched_WGCNA <- ggarrange(wgcna_taxa, cluster_glacier, labels = c("B", "C"), ncol = 2, legend = "bottom", widths= c(6.5,4), common.legend = F, align = "h")

ggarrange(wgcna_heatmap_all, enriched_WGCNA, nrow = 2, heights = c(5, 10))
# ggsave("fig6.png", plot = last_plot(), height = 16, width = 10, dpi = 900)











# Other: WGCNA cluster gas ------------------------------------------------



# # METHANE
# # set var
# env_factor <- sample_data_edit[,20]
# 
# # get correlation
# modNames = substring(names(MEs), 3)
# TaxaModuleMembership = as.data.frame(cor(otu_core_t, MEs, use = "p"));
# MMPvalue = as.data.frame(corPvalueStudent(as.matrix(TaxaModuleMembership), nSamples));
# names(TaxaModuleMembership) = paste("MM", modNames, sep="");
# names(MMPvalue) = paste("p.MM", modNames, sep="");
# 
# TaxaTraitSignificance = as.data.frame(cor(otu_core_t, env_factor, use = "p"));
# GSPvalue = as.data.frame(corPvalueStudent(as.matrix(TaxaTraitSignificance), nSamples));
# names(TaxaTraitSignificance) = paste("GS.env_var");
# names(GSPvalue) = paste("p.GS.env_var");
# 
# # get all module taxa data
# TaxaInfo0 = as.data.frame(tax_table(otu_core)[row.names(tax_table(otu_core)) %in% names(otu_core_t), ])
# TaxaInfo0$Taxon <- row.names(TaxaInfo0)
# TaxaInfo0$moduleColor <- moduleColors # add module group (this is order based!)
# 
# # add env correlation data (order is correct)
# TaxaInfo0$TaxaTraitSignificance <- TaxaTraitSignificance$GS.env_var
# TaxaInfo0$GSPvalue <- GSPvalue$p.GS.env_var
# 
# # add module membership after reformatting
# TaxaModuleMembership_clean <- TaxaModuleMembership
# TaxaModuleMembership_clean$Taxon <- row.names(TaxaModuleMembership)
# TaxaModuleMembership_clean <- gather(TaxaModuleMembership_clean, key = moduleColor, value = module_membership, -12)
# TaxaModuleMembership_clean$moduleColor <- gsub("MM", "", TaxaModuleMembership_clean$moduleColor)
# 
# TaxaInfo <- merge(TaxaInfo0, TaxaModuleMembership_clean, by = c("Taxon", "moduleColor")) # merges only membership for matching otu module
# 
# # order based on module and then env var
# TaxaOrder = order(TaxaInfo$moduleColor, -abs(TaxaInfo$TaxaTraitSignificance));
# TaxaInfo = TaxaInfo[TaxaOrder, ]
# 
# # identify top modules for env correlation
# cor(MEs, env_factor, use = "p")[order(-abs(cor(MEs, env_factor, use = "p"))),] # check most correlated modules
# 
# ggplot(TaxaInfo, aes(x = module_membership, y = TaxaTraitSignificance, colour = moduleColor)) + geom_point() +
#   ylab("Correlation") + xlab("Module Membership")  + theme_bw() + facet_wrap(~moduleColor)
# 
# ggplot(TaxaInfo[TaxaInfo$moduleColor %in% c("black", "magenta", "green", "blue"), ], aes(x = module_membership, y = TaxaTraitSignificance, colour = moduleColor)) + geom_point() +
#   ylab("Correlation") + xlab("Module Membership")  + theme_bw() + facet_wrap(~moduleColor) + geom_smooth(method = "lm") + 
#   scale_colour_manual(values = c("darkgrey", "lightblue", "palegreen3", "deeppink")) 
# # ggsave("module_ch4_membership.png", plot = last_plot(), height = 5, width = 5)
# 
# otu_abund <- as.data.frame(otu_table(otu_all_fc)[row.names(otu_table(otu_all_fc)) %in% names(otu_core_t), ]) # select core otus with FC counts
# otu_abund$OTU <- rownames(otu_abund)
# otu_abund_l <- gather(otu_abund, key = Sample, value = Abundance, -50) # format otu counts long
# 
# otu_abund_module <- merge(otu_abund_l, TaxaInfo[c(1,2)], by.x = "OTU", by.y = "Taxon") # merge counts with module info
# otu_abund_module <- merge(otu_abund_module, data.frame(sample_data_edit[,c(1,3,20,21)]), by = "Sample") # merge this with env vars
# otu_abund_module[otu_abund_module$Depth == 1, ]$Depth <- "Top" # reformat depth back
# otu_abund_module[otu_abund_module$Depth == 0, ]$Depth <- "Lower"
# 
# # plot abundance of each otu against sample methane
# # select only module of interest and in top layer
# ggplot(na.omit(otu_abund_module[otu_abund_module$moduleColor %in% c("black", "magenta", "green", "blue") & otu_abund_module$Depth == "Top", ]), aes(x = flux_CH4, y = Abundance, colour = moduleColor)) + geom_point() +
#   ylab("Abundance (CLR)") + xlab("CH4 (nmol/m2/s)")  + theme_bw() + facet_wrap(~moduleColor, scales = "free") + geom_smooth(method = "lm", position = "identity") + 
#   scale_colour_manual(values = c("darkgrey", "lightblue", "palegreen3", "deeppink")) 
# # ggsave("module_ch4_top_otu.png", plot = last_plot(), height = 5, width = 5)
# 
# 
# # sum otus per sample
# otu_abund_module_sum <- na.omit(otu_abund_module) %>% group_by(Sample, Depth, moduleColor, flux_CH4, flux_CO2) %>% summarise(Total_abundance = sum(Abundance))
# 
# # plot module abundance in each sample against sample methane (top only)
# ggplot(otu_abund_module_sum[otu_abund_module_sum$moduleColor %in% c("black", "magenta", "green", "blue") & otu_abund_module_sum$Depth == "Top", ], aes(x = flux_CH4, y = Total_abundance, colour = moduleColor)) + geom_point() +
#   ylab("Abundance (FC)") + xlab("CH4 (nmol/m2/s)")  + theme_bw() + facet_wrap(~moduleColor, nrow = 2, scales = "free") + geom_smooth(method = "lm", position = "identity") + 
#   scale_colour_manual(values = c("darkgrey", "lightblue", "palegreen3", "deeppink")) 
# # ggsave("module_ch4_top_abund_FC.png", plot = last_plot(), height = 5, width = 5)
# 
# 
# 
# 
# 
# 
# # CO2 module
# # set var
# env_factor <- sample_data_edit[,21]
# 
# # get correlation
# modNames = substring(names(MEs), 3)
# TaxaModuleMembership = as.data.frame(cor(otu_core_t, MEs, use = "p"));
# MMPvalue = as.data.frame(corPvalueStudent(as.matrix(TaxaModuleMembership), nSamples));
# names(TaxaModuleMembership) = paste("MM", modNames, sep="");
# names(MMPvalue) = paste("p.MM", modNames, sep="");
# 
# TaxaTraitSignificance = as.data.frame(cor(otu_core_t, env_factor, use = "p"));
# GSPvalue = as.data.frame(corPvalueStudent(as.matrix(TaxaTraitSignificance), nSamples));
# names(TaxaTraitSignificance) = paste("GS.env_var");
# names(GSPvalue) = paste("p.GS.env_var");
# 
# # get all module taxa data
# TaxaInfo0 = as.data.frame(tax_table(otu_core)[row.names(tax_table(otu_core)) %in% names(otu_core_t), ])
# TaxaInfo0$Taxon <- row.names(TaxaInfo0)
# TaxaInfo0$moduleColor <- moduleColors # add module group (this is order based!)
# 
# # add env correlation data (order is correct)
# TaxaInfo0$TaxaTraitSignificance <- TaxaTraitSignificance$GS.env_var
# TaxaInfo0$GSPvalue <- GSPvalue$p.GS.env_var
# 
# # add module membership after reformatting
# TaxaModuleMembership_clean <- TaxaModuleMembership
# TaxaModuleMembership_clean$Taxon <- row.names(TaxaModuleMembership)
# TaxaModuleMembership_clean <- gather(TaxaModuleMembership_clean, key = moduleColor, value = module_membership, -12)
# TaxaModuleMembership_clean$moduleColor <- gsub("MM", "", TaxaModuleMembership_clean$moduleColor)
# 
# TaxaInfo <- merge(TaxaInfo0, TaxaModuleMembership_clean, by = c("Taxon", "moduleColor")) # merges only membership for matching otu module
# 
# # order based on module and then env var
# TaxaOrder = order(TaxaInfo$moduleColor, -abs(TaxaInfo$TaxaTraitSignificance));
# TaxaInfo = TaxaInfo[TaxaOrder, ]
# 
# # identify top modules for env correlation
# cor(MEs, env_factor, use = "p")[order(-abs(cor(MEs, env_factor, use = "p"))),] # check most correlated modules
# 
# ggplot(TaxaInfo, aes(x = module_membership, y = TaxaTraitSignificance, colour = moduleColor)) + geom_point() +
#   ylab("Correlation") + xlab("Module Membership")  + theme_bw() + facet_wrap(~moduleColor)
# 
# ggplot(TaxaInfo[TaxaInfo$moduleColor %in% c("green", "yellow"), ], aes(x = module_membership, y = TaxaTraitSignificance, colour = moduleColor)) + geom_point() +
#   ylab("Correlation") + xlab("Module Membership")  + theme_bw() + facet_wrap(~moduleColor) + geom_smooth(method = "lm") + 
#   scale_colour_manual(values = c("palegreen3", "orange")) 
# # ggsave("module_co2_membership.png", plot = last_plot(), height = 3, width = 5)
# 
# 
# # using dataframe from ch4 section above
# ggplot(na.omit(otu_abund_module[otu_abund_module$moduleColor %in% c("green", "yellow") & otu_abund_module$Depth == "Top", ]), aes(x = flux_CO2, y = Abundance, colour = moduleColor)) + geom_point() +
#   ylab("Abundance (CLR)") + xlab("CO2 (umol/m2/s)")  + theme_bw() + facet_wrap(~moduleColor, scales = "free") + geom_smooth(method = "lm", position = "identity") + scale_colour_manual(values = c("palegreen3", "orange")) 
# # ggsave("module_co2_top_otu.png", plot = last_plot(), height = 3, width = 5)
# 
# 
# ggplot(otu_abund_module_sum[otu_abund_module_sum$moduleColor%in% c("green", "yellow") & otu_abund_module_sum$Depth == "Top", ], aes(x = flux_CO2, y = Total_abundance, colour = moduleColor)) + geom_point() +
#   ylab("Abundance (FC)") + xlab("CO2 (umol/m2/s)")  + theme_bw() + facet_wrap(~moduleColor, scales = "free") + geom_smooth(method = "lm", position = "identity") + scale_colour_manual(values = c("palegreen3", "orange")) 
# # ggsave("module_co2_top_abund_FC.png", plot = last_plot(), height = 3, width = 5)
# 
# 
# 
# 
# 
# 
# # SOIL DEPTH
# env_factor <- sample_data_edit_f[,1] # depth
# 
# # combine depth module comparison
# top_mod_env_all <- data.frame()
# for(mod in c("yellow", "brown", "greenyellow")){
#   top_mod_env <- TaxaInfo[TaxaInfo$moduleColor == mod, ]
#   top_mod_env <- top_mod_env[order(top_mod_env$module_membership, decreasing = T),][c(1:10),] # most significantly correlated
#   top_mod_env_all <- rbind(top_mod_env_all, top_mod_env)
# }
# # save table of top members for modules sig for glacier dist
# # write.csv(top_mod_env_all, "top_module_membership_depth.csv", row.names = F)
# 
# # get module average abundance pattern from top members
# depth_taxa_fc <- otu_all_fc
# otu_table(depth_taxa_fc) <- otu_table(depth_taxa_fc)[row.names(otu_table(depth_taxa_fc)) %in% top_mod_env_all$Taxon,] # select these otus from the transformed data above
# ntaxa(depth_taxa_fc)
# depth_taxa_fc <- psmelt(depth_taxa_fc) # df
# depth_taxa_fc <- merge(depth_taxa_fc, top_mod_env_all, by.x = "OTU", by.y = "Taxon") # combine module data back with abundance
# depth_taxa_fc_ave <- depth_taxa_fc %>% group_by(Glacier_dist, Depth, moduleColor, Genus.x, OTU) %>% summarise(Mean_abundance = mean(Abundance))
# 
# ggplot(depth_taxa_fc_ave, aes(x = Depth, y = Mean_abundance, colour = moduleColor, fill = moduleColor, group = Genus.x)) +
#   geom_area(alpha = 0.1, position = "dodge") +
#   theme_bw() + scale_fill_manual(values = c("palevioletred2", "lightgreen", "orange")) +
#   scale_colour_manual(values = c("palevioletred2", "lightgreen", "orange")) +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + scale_y_sqrt() +
#   facet_wrap(~moduleColor , ncol = 1, scales = "free_y") + ylab("Mean OTU FC Abundance") + xlab("Depth")
# 
# depth_taxa_fc_ave[depth_taxa_fc_ave$Depth == "Lower", ]$Mean_abundance <- depth_taxa_fc_ave[depth_taxa_fc_ave$Depth == "Lower", ]$Mean_abundance * (-1) # make negative
# 
# ggplot(depth_taxa_fc_ave, aes(x = Glacier_dist, y = Mean_abundance, colour = moduleColor, fill = moduleColor, group = OTU)) +
#   geom_area(alpha = 0.1, position = "dodge") +
#   theme_bw() + scale_fill_manual(values = c("palevioletred2", "lightgreen", "orange")) +
#   scale_colour_manual(values = c("palevioletred2", "lightgreen", "orange")) +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
#   facet_wrap(~moduleColor , nrow = 2, scales = "free_y") + ylab("Mean OTU FC Abundance") + xlab("Depth")
# # ggsave("module_depth_sig.png", plot = last_plot(), height = 5, width = 4)






# Methylotrophs -----------------------------------------------------------



# Methanotrophs + methylotrophs
methanotroph_genus <- c("Methyloversatilis",
                        "Methylosinus",
                        "Methylocystis",
                        "Methylocella",
                        "Methylocapsa",
                        "Methyloferula",
                        "Methylovirgula",
                        "Methyloterricola",
                        "Methylospira",
                        "Methylobacter",
                        "Methylogaea",
                        "Methylomarinum",
                        "Methylomicrobium",
                        "Methylocucumis",
                        "Methylomonas",
                        "Methylosarcina",
                        "Methylotuvimicrobium",
                        "Methylosoma",
                        "Methylomagnum",
                        "Methylosphaera",
                        "Methylovulum",
                        "LS7-MC",
                        "Methylocaldum",
                        "Methylohalobius",
                        "Methylothermus",
                        "Methylomarinovum",
                        "Methyloglobulus",
                        "Methyloparacoccus",
                        "Methyloprofundus",
                        "Methylococcus",
                        "Methylotetracoccus",
                        "Methylacidiphilum",
                        "Methylacidimicrobium",
                        "Candidatus Methylomirabilis",
                        "Crenothrix",
                        "Clonothrix",
                        
                        "Methylotenera",
                        "Methylobacterium-Methylorubrum",
                        "OM43_clade",
                        "Methyloceanibacter",
                        "Methylorosula",
                        "Hyphomicrobium",
                        "RCP2-54")

# check for more methylotrophs
methanotroph_data <- filt_data[filt_data$Genus %in% methanotroph_genus,]
methanotroph_data$Depth <- factor(methanotroph_data$Depth, levels = c("Top", "Lower"))


# sum otus for each genus
methanotroph_data_sum <- methanotroph_data %>% group_by(Sample, Site, Depth, flux_CH4, Glacier_dist, Genus) %>% summarise(Total_abundance = sum(Abundance))
# ggplot(methanotroph_data_sum, aes(x = Site, y = Total_abundance, colour = Genus)) + geom_boxplot()  + facet_wrap(~Depth + Genus, scales = "free", nrow = 2) + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + scale_y_sqrt()





# methanotroph lm
ch4_lm_top <- lm(Total_abundance ~ flux_CH4, methanotroph_data_sum[methanotroph_data_sum$Genus == "Methylocapsa" & methanotroph_data_sum$Depth == "Top",])
coef(summary(ch4_lm_top))

ch4_lm_lower <- lm(Total_abundance ~ flux_CH4, methanotroph_data_sum[methanotroph_data_sum$Genus == "Methylocapsa" & methanotroph_data_sum$Depth == "Lower",])
coef(summary(ch4_lm_lower))

# dataframe with lm values
stat_text <- data.frame(label = c(paste0("r^2==", signif(summary(ch4_lm_top)$adj.r.squared, 3), "~~p==", signif(summary(ch4_lm_top)$coef[2,4], 3)),
                                  paste0("r^2==", signif(summary(ch4_lm_lower)$adj.r.squared, 3), "~~p==", signif(summary(ch4_lm_lower)$coef[2,4], 3))),
                        Depth = c("Top", "Lower"))
stat_text$Depth <- factor(stat_text$Depth, levels = c("Top", "Lower")) # order depth

# plot methanotroph lm
ggplot(na.omit(methanotroph_data_sum[methanotroph_data_sum$Genus == "Methylocapsa",]), aes(x = flux_CH4, y = Total_abundance))  + geom_smooth(method = "lm", colour = "grey", alpha = 0.2) +
  geom_point() + theme_bw() + ylab("Methylocapsa Abundance (FC)") + xlab(expression(CH[4]~nmol~m^{-2}~s^{-1})) + geom_text(data = stat_text, aes(x = -0.6, y = -500, label = label), parse = T) +
  facet_wrap(~Depth, nrow = 2) + geom_vline(xintercept = 0, linetype = "dashed")
#ggsave("methylocapsa_lm.png", plot = last_plot(), height = 5, width = 5, dpi = 900)








# depth wilcox test
methyl_depth <- data.frame()
for(genera in unique(methanotroph_data_sum$Genus)){
  test_p <- wilcox.test(methanotroph_data_sum[methanotroph_data_sum$Genus == genera & methanotroph_data_sum$Depth == "Top", ]$Total_abundance, 
                   methanotroph_data_sum[methanotroph_data_sum$Genus == genera & methanotroph_data_sum$Depth == "Lower", ]$Total_abundance, 
                   paired = F, na.action = na.omit)$p.value # non parametric nonpaired test
  test_p_df <- data.frame(Genus = genera, p_val = test_p)
  methyl_depth <- rbind(methyl_depth, test_p_df)
}

# methanotroph_data_depth <- methanotroph_data_sum %>% group_by(Depth, Genus) %>% summarise(Mean_abundance = mean(Total_abundance))
methanotroph_data_depth <- methanotroph_data_sum %>% group_by(Depth, Genus) %>% summarise(Mean_abundance = median(Total_abundance))

# add significance symbols
labs.function = function(x){
  case_when(x > 0.055 ~ "ns",
            x <= 0.055 & x >= 0.01 ~ "*",
            x < 0.01 & x >= 0.001 ~ "**",
            x < 0.001 ~ "***")}

# depth_meth_p <- data.frame(Genus = methyl_depth$Genus, group1 = c("Lower", "Lower", "Lower", "Lower", "Lower", "Lower", "Lower"), group2 = c("Lower", "Lower", "Lower", "Lower", "Lower", "Lower", "Lower"), 
#                            p = as.numeric(methyl_depth$p_val),
#                            y.position = c(1500,1350,1200,900,810,700,400))
depth_meth_p <- data.frame(Genus = methyl_depth$Genus, group1 = c("Top", "Top", "Top", "Top", "Top", "Top", "Top"), group2 = c("Top", "Top", "Top", "Top", "Top", "Top", "Top"), 
                           p = as.numeric(methyl_depth$p_val),
                           y.position = c(1100,1350,1100,900,500,700,75))

depth_meth_p$p <- labs.function(depth_meth_p$p)

methylotroph_depth <- ggplot(methanotroph_data_depth, aes(x = Depth, y = Mean_abundance, colour = Genus)) + 
  geom_bar(aes(fill=Genus), alpha = 0.2, stat="identity", position="stack") + theme_bw() + ylab("Mean OTU Abundance (FC)") + xlab("") + 
  scale_x_discrete(labels=c("1" = "Site 1 \n Unvegetated", "1_VEG" = "Site 1 \n Vegetated")) + labs(colour = "Genera", fill = "Genera") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  stat_pvalue_manual(depth_meth_p, remove.bracket = T, hide.ns = T)

# plot only sig boxplots
methanotroph_data_depth <- methanotroph_data_sum

methylotroph_depth2 <- ggplot(methanotroph_data_depth[methanotroph_data_depth$Genus %in% c("Hyphomicrobium", "Methylotenera", "RCP2-54"), ], aes(x = Genus, y = Total_abundance, fill = Depth)) + 
  geom_boxplot(outliers = F, position = position_dodge(), alpha = 0.5) + theme_bw()+ ylab("OTU Absolute Abundance") + xlab("") + 
  labs(colour = "Depth") + scale_fill_manual(values = wes_palette("Zissou1", 4, type = "continuous")[c(3,4)]) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.title = element_blank()) + geom_signif(y_position = c(2000, 4000, 1000), xmin = c(0.8,1.8, 2.8), xmax = c(1.2,2.2,3.2), 
                                                                                  annotation = c("*", "**", "*"), tip_length = 0, colour = "black") + ylim(0,4250)



# vegetation wilcox test
methyl_veg <- data.frame()
for(genera in unique(methanotroph_data_sum$Genus)){
  print(genera)
  test_p <- wilcox.test(methanotroph_data_sum[methanotroph_data_sum$Genus == genera & methanotroph_data_sum$Site == "1", ]$Total_abundance, 
                   methanotroph_data_sum[methanotroph_data_sum$Genus == genera & methanotroph_data_sum$Site == "1_VEG", ]$Total_abundance, 
                   paired = F, na.action = na.omit)$p.value # non parametric nonpaired test
  test_p_df <- data.frame(Genus = genera, p_val = test_p)
  methyl_veg <- rbind(methyl_veg, test_p_df)
}

# methanotroph_data_site <- methanotroph_data_sum %>% group_by(Site, Glacier_dist, Genus) %>% summarise(Mean_abundance = mean(Total_abundance))
methanotroph_data_site <- methanotroph_data_sum %>% group_by(Site, Glacier_dist, Genus) %>% summarise(Mean_abundance = median(Total_abundance))

# veg_meth_p <- data.frame(Genus = methyl_veg$Genus, group1 = c("1_VEG", "1_VEG", "1_VEG", "1_VEG", "1_VEG", "1_VEG", "1_VEG"), group2 = c("1_VEG", "1_VEG", "1_VEG", "1_VEG", "1_VEG", "1_VEG", "1_VEG"), 
#                            p = as.numeric(methyl_veg$p_val),
#                            y.position = c(5600,5300,5100,5000,3000,200,100))

veg_meth_p <- data.frame(Genus = methyl_veg$Genus, group1 = c("1", "1", "1", "1", "1", "1", "1"), group2 = c("1", "1", "1", "1", "1", "1", "1"), 
                         p = as.numeric(methyl_veg$p_val),
                         y.position = c(1200,5300,5100,5000,3000,200,100))
veg_meth_p$p <- labs.function(veg_meth_p$p)

methylotroph_veg <- ggplot(methanotroph_data_site[methanotroph_data_site$Site %in% c("1", "1_VEG"), ], aes(x = Site, y = Mean_abundance, colour = Genus)) + 
  geom_bar(aes(fill=Genus), alpha = 0.2, stat="identity", position="stack") + theme_bw()+ ylab("") + xlab("") + 
  scale_x_discrete(labels=c("1" = "Site 1 \n Unvegetated", "1_VEG" = "Site 1 \n Vegetated")) + labs(colour = "Genera", fill = "Genera") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))  +
  stat_pvalue_manual(veg_meth_p, remove.bracket = T, hide.ns = T)

# only sig boxplot
methylotroph_veg2 <- ggplot(methanotroph_data_sum[methanotroph_data_sum$Site %in% c("1", "1_VEG") & methanotroph_data_sum$Genus == "Hyphomicrobium", ], aes(x = Genus, y = Total_abundance, fill = Site)) + 
  geom_boxplot(outliers = F, alpha = 0.5) + theme_bw()+ ylab("OTU Absolute Abundance") + xlab("") + 
  labs(colour = "Site") +  scale_fill_manual(labels = c("Unvegetated", "Vegetated"), values = wes_palette("Zissou1", 4, type = "continuous")[c(1,2)]) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.title = element_blank()) + geom_signif(y_position = c(2500), xmin = c(0.8), xmax = c(1.2), 
                                                                                  annotation = c("**"), tip_length = 0, colour = "black") + ylim(0, 2750)




# across all sites
methylotroph_dist <- ggplot(methanotroph_data_site[methanotroph_data_site$Site %in% c(0:5), ], aes(x = Site, y = Mean_abundance, colour = Genus)) + 
  geom_bar(aes(fill=Genus), alpha = 0.2, stat="identity", position="stack") + theme_bw() + ylab("") + xlab("") + 
  scale_x_discrete(labels=c("0" = "Glacier", "1" = "Site 1", "2" = "Site 2", "3" = "Site 3", "4" = "Site 4", "5" = "Site 5")) + labs(colour = "Genera", fill = "Genera") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

methylotroph_dist2 <- ggplot(methanotroph_data_sum[methanotroph_data_sum$Site %in% c(0:5) & methanotroph_data_sum$Genus != "OM43_clade" , ], aes(x = Site, y = Total_abundance, colour = Genus)) + 
  geom_boxplot(outliers = F) + theme_bw() + ylab("OTU Absolute Abundance") + xlab("") + 
  scale_x_discrete(labels=c("0" = "Glacier", "1" = "Site 1", "2" = "Site 2", "3" = "Site 3", "4" = "Site 4", "5" = "Site 5")) + labs(colour = "Genera", fill = "Genera") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position = "none") + facet_wrap(~Genus, scales = "free_y")
#ggsave("methylotroph_box.png", plot = last_plot(), height = 6, width = 10, dpi = 900)


ggarrange(methylotroph_depth, methylotroph_veg, methylotroph_dist,  labels = c("A", "B", "C"), 
          common.legend = T, nrow = 1, legend = "right", widths = c(1,1,2), align = "h")
#ggsave("fig6.png", plot = last_plot(), height = 6, width = 10, dpi = 900)
#ggsave("fig6.svg", plot = last_plot(), height = 6, width = 10, dpi = 900, device = svglite::svglite)

ggarrange(ggarrange(methylotroph_depth2, methylotroph_veg2, labels = c("A", "B"), widths = c(6,4.5)), methylotroph_dist2,  labels = c("", "C"),
          nrow = 2,  align = "h")
#ggsave("methylotroph_boxplots_all.png", plot = last_plot(), height = 8, width = 8, dpi = 900)





# anova between sites
site_meth <- methanotroph_data_sum[methanotroph_data_sum$Site %in% c(0:5), ]

# test different anovas
aov_1 <- aov(Total_abundance ~ Site, data = site_meth)

# compare anova models
model.set <- list(aov_1, aov_2, aov_3)
model.names <- c("aov_1", "aov_2", "aov_3")
aictab(model.set, modnames = model.names) # model with lowest AIC score is best fit, model 1

summary(aov_1)

# check stats
par(mfrow=c(2,2))
plot(aov_1)
par(mfrow=c(1,1))

# check normality
hist(aov_1$residuals) 

# post hoc test to identify what these differences are
tukey.two.way <- TukeyHSD(aov_1)
tukey.two.way

# label sig groups
# https://stackoverflow.com/questions/70483632/how-to-show-tukey-groups-using-tidyverse-and-rstatix tutorial
clust_groups <- multcompLetters4(aov_1, tukey.two.way) # input is anova and tukey
clust_groups <- as.data.frame.list(clust_groups$Site) # create a dataframe for plotting, using clust_groups$variable





# ice nucleation ----------------------------------------------------------

## ICE NUCLEATION 
IN_species <- c("Pseudomonas syringae",
                "Pseudomonas fluorescens",
                "Pseudomonas trivalis",
                "Pseudomonas costantinii",
                "Pseudomonas orientalis",
                "Pseudomonas mandelii",
                "Pseudomonas frederiksbergensis",
                "Pseudomonas lini",
                "Pseudomonas abietaniphila",
                "Pseudomonas_E_647464 abietaniphila",
                "Pseudomonas graminis",
                "Pseudomonas_E_647464 graminis",
                "Pseudomonas putida",
                "Pseudomonas aeruginosa",
                "Pseudomonas viridiflava",
                "Xanthomonas campestris",
                "Erwinia Herbicola",
                "Lysinibacillus parviboronicapiens",
                "Lysinibacillus fusiformis",
                "Lysinibacillus sphaericus")

IN_genus <- c("Pseudomonas",
              "Pseudoxanthomonas",
              "Erwinia", "Pantoea",
              "Lysinibacillus", "Xanthomonas")

IN_data <- filt_data[filt_data$Genus %in% IN_genus,]
ggplot(IN_data, aes(x = Sample, y = Abundance, fill = Genus)) + geom_bar(aes(fill=Genus), stat="identity", position="identity") + facet_grid(~Site, scales = "free_x") + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

# summarise across samples
IN_data_total <- IN_data %>% group_by(Site, Depth, Sample, Glacier_dist, Family, Genus) %>% summarise(total_abundance = sum(Abundance))


# plot site and depth
IN_data_mean <- IN_data_total %>% group_by(Depth, Site, Glacier_dist, Family, Genus) %>% summarise(mean_abundance = mean(total_abundance), sd = sd(total_abundance))
IN_data_mean$Depth <- factor(IN_data_mean$Depth, levels = c("Top", "Lower")) # order

ggplot(IN_data_mean, aes(x = Site, y = mean_abundance, fill = Genus, colour = Genus)) + 
  geom_bar(aes(fill=Genus), alpha = 0.2, stat="identity", position="stack") + theme_bw() + 
  ylab("Mean Cell Abundance") + xlab("")  +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + facet_wrap(~Depth, nrow = 2) +
  scale_x_discrete(labels=c("0" = "Glacier", "1" = "Site 1",  "1_VEG" = "Site 1 \n Vegetated", "2" = "Site 2", "3" = "Site 3", "4" = "Site 4", "5" = "Site 5")) 
#ggsave("IN_genera_bar.png", plot = last_plot(), height = 7, width = 7, dpi = 900)


wilcox.test(IN_data_total[IN_data_total$Genus == "Pseudomonas" & IN_data_total$Depth == "Top", ]$total_abundance, 
            IN_data_total[IN_data_total$Genus == "Pseudomonas" & IN_data_total$Depth == "Lower", ]$total_abundance, 
            paired = F, na.action = na.omit)$p.value

wilcox.test(IN_data_total[IN_data_total$Genus == "Pseudoxanthomonas" & IN_data_total$Depth == "Top", ]$total_abundance, 
            IN_data_total[IN_data_total$Genus == "Pseudoxanthomonas" & IN_data_total$Depth == "Lower", ]$total_abundance, 
            paired = F, na.action = na.omit)$p.value

# not sig difference between depths



# Icamp Assembly ----------------------------------------------------------



# icamp
# identify community assembly processes
comm <- as.data.frame(t(as.data.frame(otu_table(otu_all)))) # transposed otu_table

## identify community assembly processes importance between samples
# separates taxa into bins 
# processes governing each bin is tested using beta Net Relatedness Index (betaNRI) and taxonomic beta diversities (Raup-Crick)
# fraction of pairwise comparisons with betaNRI < -1.96 = homogeneous selection %
# betaNRI > +1.96 = heterogeneous selection %
# for remaining abs(NRI) <= 1.96, RC analysis is used
# fraction of pairwise comparisons with RC < -0.95 = homogenizing dispersal %
# RC > 0.95 as dispersal limitation
# for remaining abs(NRI) <= 1.96 and abs(RC) <= 0.95 = drift, diversification, weak selection and/or weak dispersal %

# % of individual processes across all bins are weighted by relative abundance of each bin, and summarized to estimate relative importance of individual processes at community level
# null model significance can also be inferred by direct test based on null model distribution, which should be a preferred choice when the null model simulated values do not follow normal distribution

# icamp.out <- icamp.big(comm = comm,
#                     tree = tree, # using rooted tree
#                     pd.wd = paste0(tempdir(),"/pdbig"), # working directory for program
#                     rand = 1000, # recommended
#                     nworker = 2, # thread number
#                     bin.size.limit = 12) # taxa bin size of 12, 24, or 48 recommended
# "1.719275 hours"

# load saved results
load("iCAMP.iCAMP.Confidence.detail.rda")
results <- res$CbMPDiCBraya # process importance between samples

# test normality of null values
# reqiures icamp.big in which detail.null must be TRUE, to save all null values
# norm_test <- null.norm(icamp.output = icamp.out, rand.list = NULL, p.norm.cut = 0.05, detail.out = T)
# results <- icamp.out$CbMPDiCBraya # process importance between samples

results$sample <- row.names(res$CbMPDiCBraya)
ggplot(results, aes(x= sample1, y = sample2, fill = Heterogeneous.Selection)) + geom_tile() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 


# importance of each process in each bin and turnover
#icamp.bins.out <- icamp.bins(icamp.out$detail, treat = as.matrix(s_data), clas = NULL, silent = FALSE,
#           boot = FALSE, rand.time = 1000, between.group = F) # quantify importance


# bootstrapping to estimate variation of relative importance of each process within each group and compare between groups
# between sites
s_data <- sample_data(otu_all)
s_data <- s_data[,c(2)] # dataframe of sites
s_data <- s_data[s_data$Site != "1_VEG", ]  # remove vegetated site

# set between.group = TRUE for between site turnovers
icamp.boot.out <- icamp.boot(res$CbMPDiCBraya, as.matrix(s_data), rand.time = 1000, compare = TRUE,
                             silent = FALSE, between.group = T, ST.estimation = FALSE)

site_processes <- icamp.boot.out$summary[,c(1:10)] # select relative importance of assembly processes

site_processes_filt <- site_processes[site_processes$Group %in% c("0","1","2","3","4","5"), ]

ggplot(site_processes_filt, aes(x = Group, y = as.numeric(Observed), fill = Process)) + geom_bar(stat="identity", position="stack") + 
  scale_fill_manual(values = wes_palette("Zissou1", 5, type = "continuous")) + theme_bw() + ylab("Mean Relative Importance") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

process_details <- data.frame()


for(group in c("0", "1", "2", "3", "4", "5")){
  print(group)
  process_details_group <- as.data.frame(icamp.boot.out$boot.detail[[group]]) # get raw data from group
  process_details_group$Group <- group
  process_details <- rbind(process_details, process_details_group)
}

# perform ANOVA
library(AICcmodavg)
library(multcompView) # https://stackoverflow.com/questions/70483632/how-to-show-tukey-groups-using-tidyverse-and-rstatix tutorial

# for one process example
# aov_1 <- aov(process_details[["Heterogeneous.Selection"]] ~ Group, data = process_details)
# summary(aov_1)
# 
# par(mfrow=c(2,2))
# plot(aov_1)
# par(mfrow=c(1,1))
# 
# hist(aov_1$residuals) # check normality
# 
# # post hoc test to identify what these differences are
# tukey.two.way <- TukeyHSD(aov_1)
# tukey.two.way
# 
# # label sig groups
# clust_groups <- multcompLetters4(aov_1, tukey.two.way)
# clust_groups <- as.data.frame.list(clust_groups$Group) # create a dataframe for plotting

# # perform anova for each assembly process
# aov_process_groups <- data.frame()
# for(process in colnames(process_details[-6])){ # excluding group col
#   print(process)
#   aov_1 <- aov(process_details[[process]] ~ Group, data = process_details)
#   hist(aov_1$residuals) # check normality
#   print(summary(aov_1))
#   
#   tukey.two.way <- TukeyHSD(aov_1)
#   print(tukey.two.way)
#   
#   clust_groups <- multcompLetters4(aov_1, tukey.two.way)
#   clust_groups <- as.data.frame.list(clust_groups$Group) # create a dataframe for plotting
#   aov_groups <- clust_groups["Letters"]
#   aov_groups$Process <- process
#   aov_groups$Group <- row.names(aov_groups)
#   aov_process_groups <- rbind(aov_process_groups, aov_groups)
# }
# # check printout here for anova F and P values and tukey P values and CI
# 
# # combine with bar plot
# site_process_filt_aov <- merge(site_processes_filt, aov_process_groups, by = c("Group", "Process") )
# 
# # order
# site_processes_filt$Process <- factor(site_processes_filt$Process, levels = c("Dispersal.Limitation",  "Drift.and.Others", "Heterogeneous.Selection", "Homogeneous.Selection", "Homogenizing.Dispersal"))
# site_processes_filt <- site_processes_filt[order(site_processes_filt$Process), ]
# 
# site_process_filt_aov$Process <- factor(site_process_filt_aov$Process, levels = c("Dispersal.Limitation",  "Drift.and.Others", "Heterogeneous.Selection", "Homogeneous.Selection", "Homogenizing.Dispersal"))
# site_process_filt_aov <- site_process_filt_aov[order(site_process_filt_aov$Process), ]
# 
# assembly_site <- ggplot(site_processes_filt, aes(x = Group, y = as.numeric(Observed), fill = Process)) + geom_col(position="stack") + 
#   scale_fill_manual(values = wes_palette("Zissou1", 5, type = "continuous")) + theme_bw() + ylab("") + xlab("") +
#   geom_text(aes(x = site_process_filt_aov$Group, label = site_process_filt_aov$Letters), position = position_stack(vjust = .5)) + 
#   scale_x_discrete(labels = c("Glacier", "Site 1", "Site 2", "Site 3", "Site 4", "Site 5")) + 
#   theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1), strip.text.x = element_blank())
# #ggsave("assembly_processes_site_aov.png", plot = last_plot(), height = 6, width = 5)
# 
# # separate bar plot for each process
# ggplot(site_process_filt_aov, aes(x = Group, y = as.numeric(Observed), fill = Process)) + geom_col(position="stack", colour = "black") + 
#   scale_fill_manual(values = wes_palette("Zissou1", 5, type = "continuous")) + theme_bw() + xlab("") + ylab("Mean Relative Importance") + 
#   theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1), strip.text.x = element_blank()) + 
#   geom_errorbar(aes(ymin = as.numeric(Observed), ymax = as.numeric(Observed) + as.numeric(Stdev)), width = 0.5) + 
#   facet_wrap(~Process, scales = "free", nrow = 1) +
#   geom_text(aes(x = site_process_filt_aov$Group, label = site_process_filt_aov$Letters), position = position_stack(vjust = .5))
# #ggsave("assembly_processes_site_aov2.png", plot = last_plot(), height = 4, width = 12)



# compute statistical differences between sites
site_processes_p <- icamp.boot.out$compare[icamp.boot.out$compare$Group1 %in% c("0","1","2","3","4","5") & icamp.boot.out$compare$Group2 %in% c("0","1","2","3","4","5"), ]

# select columns for effect size and pvalue
site_processes_p_cohen <- select(site_processes_p, c("Group1", "Group2", "Heterogeneous.Selection_Cohen.d", "Homogeneous.Selection_Cohen.d",
                                                     "Dispersal.Limitation_Cohen.d", "Homogenizing.Dispersal_Cohen.d", "Drift.and.Others_Cohen.d"))

site_processes_p_val <- select(site_processes_p, c("Group1", "Group2", "Heterogeneous.Selection_P.value", "Homogeneous.Selection_P.value",
                                                   "Dispersal.Limitation_P.value", "Homogenizing.Dispersal_P.value", "Drift.and.Others_P.value"))

# make long format
site_processes_p_cohen <- gather(site_processes_p_cohen, key = "Process", "Cohen_d", -c(1,2))
site_processes_p_cohen$Process <- gsub("_Cohen.d", "", site_processes_p_cohen$Process)

site_processes_p_val <- gather(site_processes_p_val, key = "Process", "p_val", -c(1,2))
site_processes_p_val$Process <- gsub("_P.value", "", site_processes_p_val$Process)

# merge
site_processes_p_all <- merge(site_processes_p_cohen, site_processes_p_val, by = c("Group1", "Group2", "Process"))

# padjustment for multiple testing
site_processes_p_all$p_adj <- p.adjust(site_processes_p_all$p_val, method = "fdr")

# use multcompletters to annotate sig groups
# create a matrix of site comparisons and associated pvalues for each process from bootstrapping
site_process_sig <- data.frame()
for(process in unique(site_processes_p_all$Process)){
  print(process)
  site_process_data <- site_processes_p_all[site_processes_p_all$Process == process, ] # dataframe of one process
  groups <- c("0", "1", "2", "3", "4", "5")
  pmat <- matrix(1, length(groups), length(groups), dimnames=list(groups, groups)) # create matrix
  cmat <- matrix(1, length(groups), length(groups), dimnames=list(groups, groups)) # create matrix
  
  for(i in 1:nrow(site_process_data)) { # fill matrix with pvals
    pmat[site_process_data$Group1[i], site_process_data$Group2[i]] <- site_process_data$p_adj[i]
    pmat[site_process_data$Group2[i], site_process_data$Group1[i]] <- site_process_data$p_adj[i]
  }
  
  for(i in 1:nrow(site_process_data)) { # fill matrix with effect vals
    cmat[site_process_data$Group1[i], site_process_data$Group2[i]] <- site_process_data$Cohen_d[i]
    cmat[site_process_data$Group2[i], site_process_data$Group1[i]] <- site_process_data$Cohen_d[i]
  }
  
  print(pmat)
  print(cmat)
  sig_sites <- pmat < 0.05 & abs(as.numeric(cmat)) >= 0.5  # find sig for pvalue and cohen effect
  print(sig_sites)

  # assign letter groups
  if(process == "Drift.and.Others"){
    letters <- multcompLetters(sig_sites)$Letters
  }else{
    letters <- multcompLetters(sig_sites, reversed = T)$Letters # letters are not arranged by descending so have to manually check this
  }

  # associate letter with group and process and save
  process_site_result <- data.frame(
    Group = names(letters),
    Letter = letters,
    Process = process,
    row.names = NULL
  )

  site_process_sig <- rbind(site_process_sig, process_site_result)
}

# combine with bar plot
site_process_filt_group <- merge(site_processes_filt, site_process_sig, by = c("Group", "Process") )

# order
site_processes_filt$Process <- factor(site_processes_filt$Process, levels = c("Dispersal.Limitation",  "Drift.and.Others", "Heterogeneous.Selection", "Homogeneous.Selection", "Homogenizing.Dispersal"))
site_processes_filt <- site_processes_filt[order(site_processes_filt$Process), ]

site_process_filt_group$Process <- factor(site_process_filt_group$Process, levels = c("Dispersal.Limitation",  "Drift.and.Others", "Heterogeneous.Selection", "Homogeneous.Selection", "Homogenizing.Dispersal"))
site_process_filt_group <- site_process_filt_group[order(site_process_filt_group$Process), ]

assembly_site <- ggplot(site_processes_filt, aes(x = Group, y = as.numeric(Observed), fill = Process)) + geom_col(position="stack") +
  scale_fill_manual(values = wes_palette("Zissou1", 5, type = "continuous"), labels = c("Disperal Limitation", "Drift and Others", "Heterogeneous Selection", "Homogeneous Selection", "Homogenizing Dispersal")) + theme_bw() + ylab("") + xlab("") +
  geom_text(aes(x = site_process_filt_group$Group, label = site_process_filt_group$Letter), position = position_stack(vjust = 0.5)) +
  scale_x_discrete(labels = c("Glacier", "Site 1", "Site 2", "Site 3", "Site 4", "Site 5")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), strip.text.x = element_blank())
#ggsave("assembly_processes_site_aov.png", plot = last_plot(), height = 6, width = 5)

# # separate bar plot for each process
ggplot(site_processes_filt, aes(x = Group, y = as.numeric(Observed), fill = Process)) + geom_col(position="stack", colour = "black") +
  scale_fill_manual(values = wes_palette("Zissou1", 5, type = "continuous")) + theme_bw() + xlab("") + ylab("Mean Relative Importance") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1), strip.text.x = element_blank()) +
  geom_errorbar(aes(ymin = as.numeric(Observed), ymax = as.numeric(Observed) + as.numeric(Stdev)), width = 0.5) +
  facet_wrap(~Process, scales = "free", nrow = 1) +
  geom_text(aes(x = site_process_filt_group$Group, label = site_process_filt_group$Letter), position = position_stack(vjust = .5))
#ggsave("assembly_processes_site_aov2.png", plot = last_plot(), height = 4, width = 12)









# total outwash plain
s_data <- sample_data(otu_all)
s_data$Site <- "Total" # assign same group to all samples
s_data <- s_data[,c(2)] # dataframe of site
icamp.boot.out_t <- icamp.boot(res$CbMPDiCBraya, as.matrix(s_data), rand.time = 1000, compare = FALSE,
                               silent = FALSE, between.group = FALSE, ST.estimation = FALSE)

total_processes <- icamp.boot.out_t$summary[,c(1:10)]

# plot total outwashplain assembly processes
assembly_total <- ggplot(total_processes, aes(x = Group, y = as.numeric(Observed), fill = Process)) + geom_bar(stat="identity", position="stack") + 
  scale_fill_manual(values = wes_palette("Zissou1", 5, type = "continuous"), labels = c("Disperal Limitation", "Drift and Others", "Heterogeneous Selection", "Homogeneous Selection", "Homogenizing Dispersal")) + theme_bw() + ylab("Mean Relative Importance") + xlab("") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + scale_x_discrete(labels = c("Total"))
# ggsave("assembly_processes_total.png", plot = last_plot(), height = 5, width = 4)





# Depth
s_data <- sample_data(otu_all)
# s_data$Site_depth <- paste(s_data$Site, s_data$Depth, sep = "_")
s_data <- s_data[, c(3)]
icamp.boot.out_d <- icamp.boot(res$CbMPDiCBraya, as.matrix(s_data), rand.time = 1000, compare = T,
                               silent = FALSE, between.group = F, ST.estimation = FALSE)

depth_processes <- icamp.boot.out_d$summary[,c(1:10)]
depth_processes$Group <- factor(depth_processes$Group, levels = c("Top", "Lower"))

ggplot(depth_processes, aes(x = Group, y = as.numeric(Observed), fill = Process)) + geom_bar(stat="identity", position="stack") + 
  scale_fill_manual(values = wes_palette("Zissou1", 5, type = "continuous")) + theme_bw() + ylab("Mean Relative Importance") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
# ggsave("assembly_processes_depth0.png", plot = last_plot(), height = 6, width = 4)

# per process
ggplot(depth_processes, aes(x = Group, y = as.numeric(Mean), fill = Process)) + geom_bar(stat="identity", position="stack", colour = "black") + 
  scale_fill_manual(values = wes_palette("Zissou1", 5, type = "continuous")) + theme_bw() + ylab("Mean Relative Importance") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  geom_errorbar(aes(ymin = as.numeric(Mean), ymax = as.numeric(Mean) + as.numeric(Stdev)), width = 0.5) + 
  facet_wrap(~Process, scales = "free", nrow = 1) 
# ggsave("assembly_processes_depth.png", plot = last_plot(), height = 4, width = 9)

# statistical difference between depths
depth_processes$Mean <- as.numeric(depth_processes$Mean)

# format pvalue from bootstrapping
depth_process_sig <-icamp.boot.out_d$compare # original bootstrap comparison
depth_process_sig <- as.data.frame(t(depth_process_sig))
depth_process_sig$parameter <- row.names(depth_process_sig)
depth_process_sig$stat <- gsub(".*_", "", depth_process_sig$parameter) # get stats
depth_process_sig$process <- gsub("_.*", "", depth_process_sig$parameter) # get process # check cohen d here
depth_process_sig <- depth_process_sig[depth_process_sig$stat == "P.value",] # select only pval stat
colnames(depth_process_sig)[1] <- "pval" # rename pval col
depth_process_sig$group1 <- paste(depth_process_sig$process, "Lower", sep = "_")
depth_process_sig$group2 <- paste(depth_process_sig$process, "Top", sep = "_")

depth_process_p <- data.frame(group1 = depth_process_sig$group1, group2 = depth_process_sig$group2,
                              p = as.numeric(depth_process_sig$pval),
                              y.position = c(0.2, 0.05, 0.65, 0.1, 0.3))

labs.function = function(x){
  case_when(x > 0.055 ~ "ns",
            x <= 0.055 & x >= 0.01 ~ "*",
            x < 0.01 & x >= 0.001 ~ "**",
            x < 0.001 ~ "***")}

depth_process_p$p <- labs.function(depth_process_p$p)
depth_processes$Group2 <- paste(depth_processes$Process, depth_processes$Group, sep = "_")

# plot differences
ggplot(depth_processes, aes(x = Group2, y = as.numeric(Mean))) + geom_bar(stat="identity", position="stack", colour = "black", aes(fill = Process)) + 
  scale_fill_manual(values = wes_palette("Zissou1", 5, type = "continuous")) + theme_bw() + ylab("Mean Relative Importance") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  geom_errorbar(aes(ymin = as.numeric(Mean), ymax = as.numeric(Mean) + as.numeric(Stdev)), width = 0.5) +
  stat_pvalue_manual(depth_process_p, tip.length = 0.01, hide.ns = F)
# ggsave("assembly_processes_depth.png", plot = last_plot(), height = 6, width = 6)

depth_process_p <- data.frame(Process = depth_process_sig$process, group1 = c("Lower", "Lower", "Lower", "Lower", "Lower"), group2 = c("Lower", "Lower", "Lower", "Lower", "Lower"), 
                              p = as.numeric(depth_process_sig$pval),
                              y.position = c(0.1, 0.055, 0.65, 0.02, 0.3))
depth_process_p$p <- labs.function(depth_process_p$p)
assembly_depth <- ggplot(depth_processes, aes(x = Group, y = as.numeric(Observed), fill = Process)) + geom_bar(stat="identity", position="stack") + 
  scale_fill_manual(values = wes_palette("Zissou1", 5, type = "continuous"), labels = c("Disperal Limitation", "Drift and Others", "Heterogeneous Selection", "Homogeneous Selection", "Homogenizing Dispersal")) + theme_bw() + ylab("") + xlab("") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  stat_pvalue_manual(depth_process_p, remove.bracket = T, hide.ns = F)
# ggsave("assembly_processes_depth0.png", plot = last_plot(), height = 6, width = 4)

ggplot(depth_processes[depth_processes$Process == "Heterogeneous.Selection"  | depth_processes$Process == "Dispersal.Limitation" , ], aes(x = Group2, y = as.numeric(Mean))) + geom_bar(stat="identity", position="stack", colour = "black", aes(fill = Process)) + 
  scale_fill_manual(values = wes_palette("Zissou1", 5, type = "continuous")[c(1,3)]) + theme_bw() + ylab("Mean Relative Importance") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) + xlab("") +
  geom_errorbar(aes(ymin = as.numeric(Mean), ymax = as.numeric(Mean) + as.numeric(Stdev)), width = 0.5) +
  stat_pvalue_manual(depth_process_p[c(1,3),], tip.length = 0.01, hide.ns = F) + 
  scale_x_discrete(labels=c("Dispersal.Limitation_Lower" = "Lower", "Dispersal.Limitation_Top" = "Top", "Heterogeneous.Selection_Lower" = "Lower", "Heterogeneous.Selection_Top" = "Top"))
# ggsave("assembly_processes_depth.png", plot = last_plot(), height = 5, width = 5)




## vegetation
s_data <- sample_data(otu_all)
# s_data$Site_depth <- paste(s_data$Site, s_data$Depth, sep = "_")
s_data <- s_data[, c(2)]
s_data <- s_data[s_data$Site == "1" | s_data$Site == "1_VEG", ]
icamp.boot.out_v <- icamp.boot(res$CbMPDiCBraya, as.matrix(s_data), rand.time = 1000, compare = T,
                               silent = FALSE, between.group = F, ST.estimation = FALSE)

veg_processes <- icamp.boot.out_v$summary[,c(1:10)]

ggplot(veg_processes, aes(x = Group, y = as.numeric(Observed), fill = Process)) + geom_bar(stat="identity", position="stack") + 
  scale_fill_manual(values = wes_palette("Zissou1", 5, type = "continuous")) + theme_bw() + ylab("Mean Relative Importance") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

# per process
ggplot(veg_processes, aes(x = Group, y = as.numeric(Mean), fill = Process)) + geom_bar(stat="identity", position="stack", colour = "black") + 
  scale_fill_manual(values = wes_palette("Zissou1", 5, type = "continuous")) + theme_bw() + ylab("Mean Relative Importance") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  geom_errorbar(aes(ymin = as.numeric(Mean), ymax = as.numeric(Mean) + as.numeric(Stdev)), width = 0.5) + 
  facet_wrap(~Process, scales = "free", nrow = 1) 
# ggsave("assembly_processes_veg.png", plot = last_plot(), height = 4, width = 9)

# check significance
veg_processes$Group2 <- paste(veg_processes$Process, veg_processes$Group, sep = "_")

ggplot(veg_processes, aes(x = Group2, y = as.numeric(Mean), fill = Process)) + geom_bar(stat="identity", position="stack", colour = "black") + 
  scale_fill_manual(values = wes_palette("Zissou1", 5, type = "continuous")) + theme_bw() + ylab("Mean Relative Importance") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  geom_errorbar(aes(ymin = as.numeric(Mean), ymax = as.numeric(Mean) + as.numeric(Stdev)), width = 0.5) 

# statistical difference between depths
veg_processes$Mean <- as.numeric(veg_processes$Mean)

# format pvalue from bootstrapping
veg_process_sig <-icamp.boot.out_v$compare # original bootstrap comparison
veg_process_sig <- as.data.frame(t(veg_process_sig))
veg_process_sig$parameter <- row.names(veg_process_sig)
veg_process_sig$stat <- gsub(".*_", "", veg_process_sig$parameter) # get stats
veg_process_sig$process <- gsub("_.*", "", veg_process_sig$parameter) # get process
veg_process_sig <- veg_process_sig[veg_process_sig$stat == "P.value",] # select only pval stat
colnames(veg_process_sig)[1] <- "pval" # rename pval col
veg_process_sig$group1 <- paste(veg_process_sig$process, "1", sep = "_")
veg_process_sig$group2 <- paste(veg_process_sig$process, "1_VEG", sep = "_")

veg_process_p <- data.frame(group1 = veg_process_sig$group1, group2 = veg_process_sig$group2, 
                            p = as.numeric(veg_process_sig$pval),
                            y.position = c(0.2, 0.05, 0.65, 0.1, 0.3))

veg_process_p$p <- labs.function(veg_process_p$p)
ggplot(veg_processes, aes(x = Group2, y = as.numeric(Mean))) + geom_bar(stat="identity", position="stack", colour = "black", aes(fill = Process)) + 
  scale_fill_manual(values = wes_palette("Zissou1", 5, type = "continuous")) + theme_bw() + ylab("Mean Relative Importance") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + xlab("Vegetation") +
  geom_errorbar(aes(ymin = as.numeric(Mean), ymax = as.numeric(Mean) + as.numeric(Stdev)), width = 0.5)  +
  stat_pvalue_manual(veg_process_p[3,], tip.length = 0.02, hide.ns = F)

veg_process_p <- data.frame(Process = veg_process_sig$process, group1 = c("1_VEG", "1_VEG", "1_VEG", "1_VEG", "1_VEG"), group2 = c("1_VEG", "1_VEG", "1_VEG", "1_VEG", "1_VEG"), 
                            p = as.numeric(veg_process_sig$pval),
                            y.position = c(0.12, 0.03, 0.65, 0.08, 0.3))
veg_process_p$p <- labs.function(veg_process_p$p)
assembly_veg <- ggplot(veg_processes, aes(x = Group, y = as.numeric(Observed), fill = Process)) + geom_bar(stat="identity", position="stack") + 
  scale_fill_manual(values = wes_palette("Zissou1", 5, type = "continuous"), labels = c("Disperal Limitation", "Drift and Others", "Heterogeneous Selection", "Homogeneous Selection", "Homogenizing Dispersal")) + theme_bw() + ylab("") + xlab("") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  stat_pvalue_manual(veg_process_p, remove.bracket = T, hide.ns = F) +  scale_x_discrete(labels=c("1" = "Site 1 \n Unvegetated", "1_VEG" = "Site 1 \n Vegetated"))
# ggsave("assembly_processes_veg0.png", plot = last_plot(), height = 6, width = 4)


icamp_comps <- ggarrange(assembly_total, assembly_depth, assembly_veg, assembly_site, labels = c("A", "B", "C", "D"), nrow = 1, widths = c(1,1.5,1.5,3), common.legend =  T, legend = "right", align = "h")
# ggsave("fig5A.png", plot = last_plot(), height = 4, width = 10, dpi = 900)





# Icamp and Distance ------------------------------------------------------



# Assembly processes vs unifrac distance
# dist_df from unifrac dist above
unifrac_dist <- UniFrac(otu_all_fc, weighted = TRUE) %>% as.matrix() # for FC normalised
dist_df <- melt(as.matrix(unifrac_dist))
dist_df = dist_df[dist_df[1] != dist_df[2],] # remove dif between same subsite
colnames(dist_df) <- c("sample1", "sample2", "unifrac")
process_dist <- merge(results, dist_df, by = c("sample1", "sample2")) # add unifrac dist to assembly processes

# check how assembly process correlates with unifrac distance between samples
fit_ddisp <- lm(process_dist$Dispersal.Limitation ~ unifrac, data = process_dist)
summary(fit_ddisp)
# plot example
ggplot(process_dist, aes(x = unifrac, y = Dispersal.Limitation)) + geom_point() +
  stat_smooth(method = "lm", col = "red") +
  labs(title = paste("Adj R2 = ",signif(summary(fit_ddisp)$adj.r.squared, 5),
                     ", p =",signif(summary(fit_ddisp)$coef[2,4], 5)))

fit_dhet <- lm(process_dist$Heterogeneous.Selection ~ unifrac, data = process_dist) # fit lm
summary(fit_dhet)

fit_dhom <- lm(process_dist$Homogeneous.Selection ~ unifrac, data = process_dist) # fit lm
summary(fit_dhom)

fit_ddrift <- lm(process_dist$Drift.and.Others ~ unifrac, data = process_dist) # fit lm
summary(fit_ddrift)

fit_dhomd <- lm(Homogenizing.Dispersal ~ unifrac, data = process_dist)
summary(fit_dhomd)

# reformat long
# sometimes a col named samples needs to be removed first (col 8)
process_dist_all <- gather(process_dist[-c(8)], key = assembly_process, value = rel_importance, -c(1,2, 8)) # except sample names and unifrac

# extract lm stats
process_values <- data.frame(assembly_process = c("Dispersal.Limitation", "Drift.and.Others", "Heterogeneous.Selection", "Homogeneous.Selection",  "Homogenizing.Dispersal"),
                                  p_val = c(signif(summary(fit_ddisp)$coef[2,4], 3),
                                            signif(summary(fit_ddrift)$coef[2,4], 3),
                                            signif(summary(fit_dhet)$coef[2,4], 3), 
                                            signif(summary(fit_dhom)$coef[2,4], 3),
                                            signif(summary(fit_dhomd)$coef[2,4], 3)
                                            ),
                                  r2 = c(signif(summary(fit_ddisp)$adj.r.squared, 2), 
                                         signif(summary(fit_ddrift)$adj.r.squared, 2),
                                         signif(summary(fit_dhet)$adj.r.squared, 2), 
                                         signif(summary(fit_dhom)$adj.r.squared, 2),
                                         signif(summary(fit_dhomd)$adj.r.squared, 2) 
                                         ),
                             slope = c(signif(coef(fit_ddisp)[2], 2),
                                       signif(coef(fit_ddrift)[2], 2),
                                       signif(coef(fit_dhet)[2], 2),
                                       signif(coef(fit_dhom)[2], 2),
                                       signif(coef(fit_dhomd)[2], 2)
                                       )
                             )

# plot
assembly_unifrac <- ggplot(process_dist_all, aes(x = unifrac, y = rel_importance, colour = assembly_process, group = assembly_process)) + geom_point(alpha = 0.3) +
  stat_smooth(method = "lm") + xlab("Weighted UniFrac")  + 
  annotate("text", y = c(0.8, 0.77, 0.74, 0.71, 0.68), x = 0.15, label = paste0("slope==", process_values$slope, "~~r^2==", process_values$r2), parse = T,
           colour = c("#3A9AB2" ,"#9BBDAC" ,"#DCCB4E" ,"#E98905" , "#F11B00"), size = 3, hjust = "left") +
  scale_colour_manual(values = c("#3A9AB2" ,"#9BBDAC" ,"#DCCB4E" ,"#E98905" , "#F11B00"), labels = c("Disperal Limitation", "Drift and Others", "Heterogeneous Selection", "Homogeneous Selection", "Homogenizing Dispersal")) + 
  theme_bw() + ylab("Relative Importance") + labs(colour = "Process")
#ggsave("assembly processes unifrac distance.png", plot = last_plot(), height = 5, width = 6)


## Assembly process distance decay

# calculate how geographic distance correlates with assembly processes
# results <- icamp.out$CbMPDiCBraya from icamp
s_data <- sample_data(otu_all)
s_data <- as.matrix(s_data[,c(1,19)]) # get glacier dist info

process_dist_g <- merge(results, s_data, by.x = "sample1", by.y = "Sample") # get glacier distance for sample1
process_dist_g <- process_dist_g[-8] # remove weird sample column (ISSUE)
colnames(process_dist_g)[8] <- "sample1_glacier_dist" # rename variable for sample1
process_dist_g <- merge(process_dist_g, s_data, by.x = "sample2", by.y = "Sample") # get glacier distance for sample1
colnames(process_dist_g)[9] <- "sample2_glacier_dist" # rename variable for sample1

# the distance between samples
process_dist_g$sample_distg <- abs(as.numeric(process_dist_g$sample1_glacier_dist) - as.numeric(process_dist_g$sample2_glacier_dist)) 

# linear regressions
fit_disp <- lm(Dispersal.Limitation ~ sample_distg, data = process_dist_g)
summary(fit_disp)
# plotting example
ggplot(process_dist_g, aes(x = sample_distg, y = Dispersal.Limitation)) + geom_point() +
  stat_smooth(method = "lm", col = "red") + xlab("Geographic distance between sites (km)")   + 
  labs(title = paste("Adj R2 = ",signif(summary(fit_disp)$adj.r.squared, 5),
                     ", p =",signif(summary(fit_disp)$coef[2,4], 5)))

fit_het <- lm(process_dist_g$Heterogeneous.Selection ~ sample_distg, data = process_dist_g)
summary(fit_het)

fit_hom <- lm(process_dist_g$Homogeneous.Selection ~ sample_distg, data = process_dist_g)
summary(fit_hom)

fit_drift <- lm(process_dist_g$Drift.and.Others ~ sample_distg, data = process_dist_g)
summary(fit_drift)

fit_homd <- lm(process_dist_g$Homogenizing.Dispersal ~ sample_distg, data = process_dist_g)
summary(fit_homd)


# reformat long
process_dist_g_all <- gather(process_dist_g, key = assembly_process, value = rel_importance, -c(1,2, 8:10))

# extract lm stats (order matters here)
process_stat_values <- data.frame(assembly_process = c("Dispersal.Limitation", "Drift.and.Others", "Heterogeneous.Selection", "Homogeneous.Selection",  "Homogenizing.Dispersal"),
                             p_val = c(signif(summary(fit_disp)$coef[2,4], 3),
                                       signif(summary(fit_drift)$coef[2,4], 3),
                                       signif(summary(fit_het)$coef[2,4], 3), 
                                       signif(summary(fit_hom)$coef[2,4], 3),
                                       signif(summary(fit_homd)$coef[2,4], 3)),
                             r2 = c(signif(summary(fit_disp)$adj.r.squared, 2), 
                                    signif(summary(fit_drift)$adj.r.squared, 2),
                                    signif(summary(fit_het)$adj.r.squared, 2), 
                                    signif(summary(fit_hom)$adj.r.squared, 2),
                                    signif(summary(fit_homd)$adj.r.squared, 2)),
                             slope = c(signif(coef(fit_disp)[2], 2),
                                       signif(coef(fit_drift)[2], 2),
                                       signif(coef(fit_het)[2], 2),
                                       signif(coef(fit_hom)[2], 2),
                                       signif(coef(fit_homd)[2], 2)
                             ))

# plot
assembly_distance <- ggplot(process_dist_g_all, aes(x = sample_distg, y = rel_importance, colour = assembly_process, group = assembly_process)) + geom_point(alpha = 0.3) +
  stat_smooth(method = "lm") + xlab("Geographic distance between sites (km)")  + 
  annotate("text", y = c(0.8, 0.77, 0.74, 0.71, 0.68), x = 0, label = paste0("slope==", process_stat_values$slope, "~~r^2==", process_stat_values$r2), parse = T,
           colour = c("#3A9AB2" ,"#9BBDAC" ,"#DCCB4E" ,"#E98905" , "#F11B00"), size = 3, hjust = "left") +
  scale_colour_manual(values = c("#3A9AB2" ,"#9BBDAC" ,"#DCCB4E" ,"#E98905" , "#F11B00"), labels = c("Disperal Limitation", "Drift and Others", "Heterogeneous Selection", "Homogeneous Selection", "Homogenizing Dispersal")) + 
  theme_bw() + ylab("Relative Importance") + labs(colour = "Process")  
#ggsave("assembly processes geo distance.png", plot = last_plot(), height = 5, width = 6)

icamp_dist <- ggarrange(assembly_unifrac, assembly_distance, common.legend = T, legend = "right", labels = c("E", "F"))
#ggsave("fig5B.png", plot = last_plot(), height = 4, width = 10, dpi = 900)


ggarrange(icamp_comps, icamp_dist, nrow = 2)
# ggsave("fig5.png", plot = last_plot(), height = 8, width = 10, dpi = 900)



# compare assembly processes between unifrac vs geographic dist
# format
process_dist_all$distance_m <- "unifrac"
process_dist_g_all$distance_m <- "geographic"
process_dist_g_all <- process_dist_g_all[-c(3,4)]
colnames(process_dist_g_all) <- colnames(process_dist_all) # check colnames and rename them to match

process_dist_all_merged <- rbind(process_dist_all, process_dist_g_all)
colnames(process_dist_all_merged)[3] <- "dist"

summary(lm(rel_importance ~ dist + distance_m + dist:distance_m, process_dist_all_merged[process_dist_all_merged$assembly_process == "Heterogeneous.Selection", ]))
# this mean that unifrac distance has a slop 0.164 higher than geographic distance for Heterogeneous.Selection importance, this is a sig difference
summary(lm(rel_importance ~ dist + distance_m + dist:distance_m, process_dist_all_merged[process_dist_all_merged$assembly_process == "Homogeneous.Selection", ]))
summary(lm(rel_importance ~ dist + distance_m + dist:distance_m, process_dist_all_merged[process_dist_all_merged$assembly_process == "Dispersal.Limitation", ]))
summary(lm(rel_importance ~ dist + distance_m + dist:distance_m, process_dist_all_merged[process_dist_all_merged$assembly_process == "Homogenizing.Dispersal", ]))
summary(lm(rel_importance ~ dist + distance_m + dist:distance_m, process_dist_all_merged[process_dist_all_merged$assembly_process == "Drift.and.Others", ]))

unique(process_dist_all_merged$assembly_process)






