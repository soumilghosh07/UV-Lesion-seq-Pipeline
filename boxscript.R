library(GenomicRanges)
library(rtracklayer)
library(ggplot2)
bedgraph_file <- "17dec_wi38"
core15_bed <- "E017_15_coreMarks_dense.bed.bgz"
out <- "./"

# Import core15 bed file
df_core15 <- import.bed(core15_bed)
mcols(df_core15)$score <- NULL
mcols(df_core15)$itemRgb <- NULL
mcols(df_core15)$thick <- NULL
df_core15 <- sort(df_core15)

# Import the bedgraph file
dat <- import(bedgraph_file, format = "bedGraph")

# Calculate the coverage RLE
dat_rle <- coverage(dat, weight = "score")

# Ensure seqlevels are identical
seqlevels(df_core15) <- names(dat_rle)

# Perform binnedAverage
wi381 <- binnedAverage(bins = df_core15, numvar = dat_rle, varname = "wi381_mean_FC", na.rm = TRUE)

# Sort seqlevels and the wi381 object
wi381 <- sortSeqlevels(wi381)
wi381 <- sort(wi381)

# Save the result
save(wi381, file = paste0(out, "wi38_17dec_binned_to_chrom15.RData"))

wi38df <- wi381
chromatin_states <- unique(mcols(wi381)$name)
wi38df <- data.frame(states = chromatin_states, min = NA, q25 = NA, median = NA, q75 = NA, max = NA)
for (state in chromatin_states) {
  state_data <- mcols(wi381)$wi381_mean_FC[mcols(wi381)$name == state]
  wi38df[wi38df$states == state, ] <- c(
    state,
    min(state_data, na.rm = TRUE),
    quantile(state_data, 0.25, na.rm = TRUE),
    median(state_data, na.rm = TRUE),
    quantile(state_data, 0.75, na.rm = TRUE),
    max(state_data, na.rm = TRUE)
  )}
whole_genome_data <- mcols(wi381)$wi381_mean_FC
whole_genome_stats <- data.frame(
  states = "Whole genome",
  min = min(whole_genome_data, na.rm = TRUE),
  q25 = quantile(whole_genome_data, 0.25, na.rm = TRUE),
  median = median(whole_genome_data, na.rm = TRUE),
  q75 = quantile(whole_genome_data, 0.75, na.rm = TRUE),
  max = max(whole_genome_data, na.rm = TRUE)
)
wi38df <- rbind(wi38df, whole_genome_stats)


print(wi38df)
wi38df$min <- as.numeric(wi38df$min)
wi38df$q25 <- as.numeric(wi38df$q25)
wi38df$median <- as.numeric(wi38df$median)
wi38df$q75 <- as.numeric(wi38df$q75)
wi38df$max <- as.numeric(wi38df$max)

BarColors <- c("Red", "Indian Red", "Dark Salmon", "Dark Khaki", "grey", "Green", "Dark Green",
               "Green Yellow", "Yellow", "Medium Aquamarine", "Pale Turquoise", "Orange Red", 
               "Lime Green", "Gainsboro", "Black", "White")


wi38df$states <- factor(wi38df$states, levels = c(chromatin_states, "Whole genome"), 
                        labels = c(chromatin_states, "Whole genome"))

plottitle <- "Ranked wi381 Signal"
xlab <- "Chromatin States "
ylab <- "wi381 (ranked signal)"

wi38_Plot <- ggplot(wi38df, aes(x = reorder(states, median))) + 
  geom_boxplot(size = 0.25, aes(ymin = min, lower = q25, middle = median, upper = q75, max = max, fill = states), stat = "identity") +
  scale_fill_manual(values = BarColors) +  # Apply colors
  ggtitle(label = plottitle) + 
  xlab(xlab) + 
  ylab(ylab) + 
  guides(fill = "none") +
  theme_bw() + 
  theme(
    plot.title = element_text(hjust = 0.5),
    line = element_line(size = 0.25),
    panel.border = element_rect(fill = NA, size = 0.25),
    axis.title.x = element_text(size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.title.y = element_text(size = 10),
    axis.text.y = element_text(size = 8)
  ) +
  coord_cartesian(ylim = c(-2, 4))  # Adjust y-axis for closer view of signal ranges


ggsave(filename = "Ranked_wi381_17nov_signal_in_chrom_states.pdf", path = out, plot = wi38_Plot, device = "pdf", width = 5, height = 3, units = "in")
