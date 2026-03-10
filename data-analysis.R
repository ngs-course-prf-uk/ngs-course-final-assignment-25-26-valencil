library(tidyverse)
setwd("~/projects/ngs-course-final-assignment-25-26-valencil")

d <- read_tsv(
  "cols-all.tsv",
  col_names = c("DP", "TYPE")
)

# Histogram
ggplot(d, aes(x = DP, fill = TYPE)) +
  geom_histogram(alpha = 0.5, position = "identity", bins = 50) +
  scale_x_log10() +
  ggtitle("Distribution of read depth (DP) qualities INDELS vs. SNPs") +
  scale_fill_manual(values = c("SNP" = "steelblue", "INDEL" = "orange"))

ggsave(filename = "~/projects/ngs-course-final-assignment-25-26-valencil/results/histogram.jpeg", plot = last_plot(), width = 6, height = 4, units = "in", dpi = 300)

# Boxplot
d %>%
  filter(TYPE %in% c("SNP","INDEL"), DP > 0) %>%
  ggplot(aes(TYPE, DP, fill = TYPE)) +
  geom_boxplot(outlier.size = 0.01) +
  ggtitle("Distribution of read depth (DP) qualities INDELS vs. SNPs") +
  scale_fill_manual(values = c("SNP" = "steelblue", "INDEL" = "orange")) +
  scale_y_log10()

ggsave(filename = "~/projects/ngs-course-final-assignment-25-26-valencil/results/boxplot.jpeg", plot = last_plot(), width = 6, height = 4, units = "in", dpi = 300)

# General statistics
d %>%
  filter(TYPE %in% c("SNP", "INDEL")) %>%
  group_by(TYPE) %>%
  summarise(
    n = n(),
    mean_DP = mean(DP, na.rm = TRUE),
    median_DP = median(DP, na.rm = TRUE),
    sd_DP = sd(DP, na.rm = TRUE),
    IQR_DP = IQR(DP, na.rm = TRUE)
  )

# Wilcox test
wilcox.test(DP ~ TYPE, data = d %>% filter(TYPE %in% c("SNP","INDEL")))