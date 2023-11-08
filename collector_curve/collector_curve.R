# Collector Curves - Figure 6

library(readr)
library(dplyr)
library(ggplot2)
library(stringr)
library(scales)
library(vegan)
library(tidyverse)

# Upload table with taxonomic hiearchy that match the OTU ID's
taxa <- read_delim("Sabree_L7/KC_L7_Filter_OTU_Name.txt", 
                   delim = "\t", escape_double = FALSE, 
                   trim_ws = TRUE)

taxa2 <- taxa[, 1:2]
names(taxa2) <- c("taxa", "otu")

split <- as.data.frame(str_split_fixed(taxa2$taxa, ";", n = 7))
names(split) <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")

taxa3 <- cbind.data.frame(split, taxa2[, 2])

taxa_fam <- taxa3[c(5, 8)]

# OTU contributing to each
otu_mc <- read_delim("Sabree_L7/PICRUST2_all/picrust2_out_pipeline_2022/pathways_out/path_abun_strat.tsv", 
                     delim = "\t", escape_double = FALSE, 
                     trim_ws = TRUE)
# convert to 0s and 1s
otu_01 <- as.data.frame(apply(otu_mc[, -c(1:2)], 2, function(x) ifelse(x == 0, 0, 1)))
otu_01 <- cbind.data.frame(otu_mc[, 1:2], otu_01)

# split table up by diet
balanced <- otu_01[, 1:57]
protein <- otu_01[, c(1:2, 58:99)]
cellulose <- otu_01[, c(1:2, 100:125)]
# removing cellulose weeks above 4
celluloseB <- cellulose[, -c(13:14, 19:22)]

# BALANCED DIET

# split up each diet into week 1 and week 4
# last number is the week #

balanced_week1 <- balanced[, c(1:2, 3, 11, 19, 27, 35, 43)]
balanced_week1_count <- cbind.data.frame(
  balanced_week1[, 1:2],
  as.data.frame(apply(balanced_week1[, c(3:8)], 1, sum))
)
names(balanced_week1_count)[3] <- "count_wk1"

balanced_week4 <- balanced[, c(1:2, 6, 14, 22, 30, 38, 46)]
balanced_week4_count <- cbind.data.frame(
  balanced_week4[, 1:2],
  as.data.frame(apply(balanced_week4[, c(3:8)], 1, sum))
)
names(balanced_week4_count)[3] <- "count_wk4"

identical(balanced_week1_count$sequence, balanced_week4_count$sequence) #TRUE

balanced_wk1.4 <- cbind.data.frame(balanced_week1_count, balanced_week4_count$count_wk4)
names(balanced_wk1.4)[4] <- "count_wk4"

# remove rows that are 0
balanced_week1_count_rm0 <- balanced_week1_count[apply(balanced_week1_count, 1, function(x) all(x!=0)), ]
balanced_week4_count_rm0 <- balanced_week4_count[apply(balanced_week4_count, 1, function(x) all(x!=0)), ]

# use table "balanced_week1_count_rm0"

# remove counts
balanced_week1_otumc <- balanced_week1_count_rm0[, c(2, 1)]

# group pathways by otu

b1otu <- balanced_week1_otumc %>%
  arrange(sequence) %>%
  as.data.frame()

b1otu_PA <- b1otu %>%
  mutate(presence = 1) %>%
  pivot_wider(names_from = sequence, 
              values_from = presence, 
              values_fill = 0) %>%
  as.data.frame()

rownames(b1otu_PA) <- b1otu_PA$pathway

k <- b1otu_PA[, -1]
k2 <- t(k)

# species accumulation calculation
b1_accum_curve <- specaccum(k2, method = "random")

plot(b1_accum_curve, col = "blue",
     xlab = "Number of pathways",
     ylab = "Number of otus")

curve_mc <- b1_accum_curve[["richness"]]
curve_otu <- b1_accum_curve[["sites"]]

x <- cbind.data.frame(curve_otu, curve_mc)

ggplot(x, aes(x = curve_otu, y = curve_mc, color = "blue")) + geom_line()

k2_taxa <- as.data.frame(rownames(k2))
names(k2_taxa) <- "otu"

k2_taxa2 <- merge(k2_taxa, taxa_fam, by = "otu")

t2 <- k2_taxa2 %>%
  sample_n(n()) %>%
  mutate(observation = row_number()) %>%
  arrange(family, observation) %>%
  group_by(family) %>%
  mutate(distinct = row_number() == 1) %>%
  ungroup() %>%
  arrange(observation) %>%
  mutate(s = cumsum(distinct))

x2 <- cbind.data.frame(x, t2$s)

zero <- c(0, 0, 0, 0)

x3 <- rbind.data.frame(zero, x2)

names(x3)[3] <- "taxa"

# I divided taxa axis by the number of families that the OTU's match to
# Pathways are divided by 337, which is the total number of pathways

# Balanced diet Week 1 - Collector Curve | Y axis is Percent
ggplot(x3, aes(x = curve_otu)) + geom_line(aes(y = ((taxa/103) * 100), color = "family"), color = "red") +
  geom_line(aes(y = (curve_mc/337 * 100), color = "pathway"), color = "blue") + 
  labs(title = "Balanced Week 1: barcharts are cumulative percent family profile \n line is cumulative percent unique families", 
       x = "OTU Count", 
       y = "Percent") + theme_classic()

# Balanced diet Week 1 - Collector Curve | Y axis is Count
ggplot(x3, aes(x = curve_otu)) + geom_line(aes(y = (taxa), color = "family"), color = "red") +
  geom_line(aes(y = (curve_mc), color = "pathway"), color = "blue") + 
  labs(title = "Balanced Week 1: barcharts are cumulative percent family profile \n line is cumulative percent unique families", 
       x = "OTU Count", 
       y = "Number") + theme_classic()


# remove counts
balanced_week4_otumc <- balanced_week4_count_rm0[, c(2, 1)]

# group pathways by otu

b4otu <- balanced_week4_otumc %>%
  arrange(sequence) %>%
  as.data.frame()

b4otu_PA <- b4otu %>%
  mutate(presence = 1) %>%
  pivot_wider(names_from = sequence, 
              values_from = presence, 
              values_fill = 0) %>%
  as.data.frame()

rownames(b4otu_PA) <- b4otu_PA$pathway

k_b4 <- b4otu_PA[, -1]
k2_b4 <- t(k_b4)


# species accumulation calculation
b4_accum_curve <- specaccum(k2_b4, method = "random")

plot(b4_accum_curve, col = "blue",
     xlab = "Number of pathways",
     ylab = "Number of otus")

b4_curve_mc <- b4_accum_curve[["richness"]]
b4_curve_otu <- b4_accum_curve[["sites"]]

x_b4 <- cbind.data.frame(b4_curve_otu, b4_curve_mc)

ggplot(x_b4, aes(x = b4_curve_otu, y = b4_curve_mc, color = "blue")) + geom_line()

b4_k2_taxa <- as.data.frame(rownames(k2_b4))
names(b4_k2_taxa) <- "otu"

b4_k2_taxa2 <- merge(b4_k2_taxa, taxa_fam, by = "otu")

b4_t2 <- b4_k2_taxa2 %>%
  sample_n(n()) %>%
  mutate(observation = row_number()) %>%
  arrange(family, observation) %>%
  group_by(family) %>%
  mutate(distinct = row_number() == 1) %>%
  ungroup() %>%
  arrange(observation) %>%
  mutate(s = cumsum(distinct))

b4_x2 <- cbind.data.frame(x_b4, b4_t2$s)

zero <- c(0, 0, 0, 0)

b4_x3 <- rbind.data.frame(zero, b4_x2)

b4_x3$taxa <- as.numeric(b4_x3$`b4_t2$s`)


# Balanced diet Week 4 - Collector Curve | Y axis is Percent
ggplot(b4_x3, aes(x = b4_curve_otu)) + geom_line(aes(y = ((taxa/103) * 100), color = "family"), color = "red") +
  geom_line(aes(y = (b4_curve_mc/337 * 100), color = "pathway"), color = "blue") + 
  labs(title = "Balanced Week 4", 
       x = "OTU Count", 
       y = "Percent") + theme_classic()

# Balanced diet Week 4 - Collector Curve | Y axis is Count
ggplot(b4_x3, aes(x = b4_curve_otu)) + geom_line(aes(y = (taxa), color = "family"), color = "red") +
  geom_line(aes(y = (b4_curve_mc), color = "pathway"), color = "blue") + 
  labs(title = "Balanced Week 4", 
       x = "OTU Count", 
       y = "Count") + theme_classic()


# PROTEIN DIET
# split up each diet into week 1 and week 4

protein_week1 <- protein[, c(1:2, 3, 21, 29, 37)]
protein_week1_count <- cbind.data.frame(
  protein_week1[, 1:2],
  as.data.frame(apply(protein_week1[, c(3:6)], 1, sum))
)
names(protein_week1_count)[3] <- "count_wk1"

protein_week4 <- protein[, c(1:2, 6, 24, 32, 40)] # trying without individuals 4 and 5 so that samples match up between week 1 and 4
protein_week4_count <- cbind.data.frame(
  protein_week4[, 1:2],
  as.data.frame(apply(protein_week4[, c(3:6)], 1, sum))
)
names(protein_week4_count)[3] <- "count_wk4"

identical(protein_week1_count$sequence, protein_week4_count$sequence) #TRUE

protein_wk1.4 <- cbind.data.frame(protein_week1_count, protein_week4_count$count_wk4)
names(protein_wk1.4)[4] <- "count_wk4"

# remove rows that are 0
protein_week1_count_rm0 <- protein_week1_count[apply(protein_week1_count, 1, function(x) all(x!=0)), ]
protein_week4_count_rm0 <- protein_week4_count[apply(protein_week4_count, 1, function(x) all(x!=0)), ]

# remove counts
protein_week1_otumc <- protein_week1_count_rm0[, c(2, 1)]

# group pathways by otu

p1otu <- protein_week1_otumc %>%
  arrange(sequence) %>%
  as.data.frame()

p1otu_PA <- p1otu %>%
  mutate(presence = 1) %>%
  pivot_wider(names_from = sequence, 
              values_from = presence, 
              values_fill = 0) %>%
  as.data.frame()

rownames(p1otu_PA) <- p1otu_PA$pathway

k_p1 <- p1otu_PA[, -1]
k2_p1 <- t(k_p1)


# species accumulation calculation
p1_accum_curve <- specaccum(k2_p1, method = "random")

plot(p1_accum_curve, col = "blue",
     xlab = "Number of pathways",
     ylab = "Number of otus")

p1_curve_mc <- p1_accum_curve[["richness"]]
p1_curve_otu <- p1_accum_curve[["sites"]]

x_p1 <- cbind.data.frame(p1_curve_otu, p1_curve_mc)

ggplot(x_p1, aes(x = p1_curve_otu, y = p1_curve_mc, color = "blue")) + geom_line()

p1_k2_taxa <- as.data.frame(rownames(k2_p1))
names(p1_k2_taxa) <- "otu"

p1_k2_taxa2 <- merge(p1_k2_taxa, taxa_fam, by = "otu")

p1_t2 <- p1_k2_taxa2 %>%
  sample_n(n()) %>%
  mutate(observation = row_number()) %>%
  arrange(family, observation) %>%
  group_by(family) %>%
  mutate(distinct = row_number() == 1) %>%
  ungroup() %>%
  arrange(observation) %>%
  mutate(s = cumsum(distinct))

p1_x2 <- cbind.data.frame(x_p1, p1_t2$s)

zero <- c(0, 0, 0, 0)

p1_x3 <- rbind.data.frame(zero, p1_x2)

p1_x3$taxa <- as.numeric(p1_x3$`p1_t2$s`)

# Protein diet Week 1 - Collector Curve | Y axis is Percent
ggplot(p1_x3, aes(x = p1_curve_otu)) + geom_line(aes(y = ((taxa/103) * 100), color = "family"), color = "red") +
  geom_line(aes(y = (p1_curve_mc/337 * 100), color = "pathway"), color = "blue") + 
  labs(title = "Protein Week 1", 
       x = "OTU Count", 
       y = "Percent") + theme_classic()

# Protein diet Week 1 - Collector Curve | Y axis is Percent
ggplot(p1_x3, aes(x = p1_curve_otu)) + geom_line(aes(y = (taxa), color = "family"), color = "red") +
  geom_line(aes(y = (p1_curve_mc), color = "pathway"), color = "blue") + 
  labs(title = "Protein Week 1", 
       x = "OTU Count", 
       y = "Count") + theme_classic()


# remove counts
protein_week4_otumc <- protein_week4_count_rm0[, c(2, 1)]

# group pathways by otu

p4otu <- protein_week4_otumc %>%
  arrange(sequence) %>%
  as.data.frame()

p4otu_PA <- p4otu %>%
  mutate(presence = 1) %>%
  pivot_wider(names_from = sequence, 
              values_from = presence, 
              values_fill = 0) %>%
  as.data.frame()

rownames(p4otu_PA) <- p4otu_PA$pathway

k_p4 <- p4otu_PA[, -1]
k2_p4 <- t(k_p4)


# species accumulation calculation
p4_accum_curve <- specaccum(k2_p4, method = "random")

plot(p4_accum_curve, col = "blue",
     xlab = "Number of pathways",
     ylab = "Number of otus")

p4_curve_mc <- p4_accum_curve[["richness"]]
p4_curve_otu <- p4_accum_curve[["sites"]]

x_p4 <- cbind.data.frame(p4_curve_otu, p4_curve_mc)

ggplot(x_p4, aes(x = p4_curve_otu, y = p4_curve_mc, color = "blue")) + geom_line()

p4_k2_taxa <- as.data.frame(rownames(k2_p4))
names(p4_k2_taxa) <- "otu"

p4_k2_taxa2 <- merge(p4_k2_taxa, taxa_fam, by = "otu")

p4_t2 <- p4_k2_taxa2 %>%
  sample_n(n()) %>%
  mutate(observation = row_number()) %>%
  arrange(family, observation) %>%
  group_by(family) %>%
  mutate(distinct = row_number() == 1) %>%
  ungroup() %>%
  arrange(observation) %>%
  mutate(s = cumsum(distinct))

p4_x2 <- cbind.data.frame(x_p4, p4_t2$s)

zero <- c(0, 0, 0, 0)

p4_x3 <- rbind.data.frame(zero, p4_x2)

p4_x3$taxa <- as.numeric(p4_x3$`p4_t2$s`)

# Protein diet Week 4 - Collector Curve | Y axis is Percent
ggplot(p4_x3, aes(x = p4_curve_otu)) + geom_line(aes(y = ((taxa/102) * 100), color = "family"), color = "red") +
  geom_line(aes(y = (p4_curve_mc/337 * 100), color = "pathway"), color = "blue") + 
  labs(title = "Protein Week 4", 
       x = "OTU Count (418)", 
       y = "Percent") + theme_classic()

# Protein diet Week 4 - Collector Curve | Y axis is Count
ggplot(p4_x3, aes(x = p4_curve_otu)) + geom_line(aes(y = (taxa), color = "family"), color = "red") +
  geom_line(aes(y = (p4_curve_mc), color = "pathway"), color = "blue") + 
  labs(title = "Protein Week 4", 
       x = "OTU Count (418)", 
       y = "Count") + theme_classic()



# Cellulose Diet
# split up each diet into week 1 and week 4

cellulose_week1 <- celluloseB[, c(1:2, 6, 10, 14)]
cellulose_week1_count <- cbind.data.frame(
  cellulose_week1[, 1:2],
  as.data.frame(apply(cellulose_week1[, c(3:5)], 1, sum))
)
names(cellulose_week1_count)[3] <- "count_wk1"

cellulose_week4 <- celluloseB[, c(1:2, 9, 12, 16)]
cellulose_week4_count <- cbind.data.frame(
  cellulose_week4[, 1:2],
  as.data.frame(apply(cellulose_week4[, c(3:5)], 1, sum))
)
names(cellulose_week4_count)[3] <- "count_wk4"

identical(cellulose_week1_count$sequence, cellulose_week4_count$sequence) #TRUE

cellulose_wk1.4 <- cbind.data.frame(cellulose_week1_count, cellulose_week4_count$count_wk4)
names(cellulose_wk1.4)[4] <- "count_wk4"

# remove rows that are 0
cellulose_week1_count_rm0 <- cellulose_week1_count[apply(cellulose_week1_count, 1, function(x) all(x!=0)), ]
cellulose_week4_count_rm0 <- cellulose_week4_count[apply(cellulose_week4_count, 1, function(x) all(x!=0)), ]


# remove counts
cellulose_week1_otumc <- cellulose_week1_count_rm0[, c(2, 1)]

# group pathways by otu

c1otu <- cellulose_week1_otumc %>%
  arrange(sequence) %>%
  as.data.frame()

c1otu_PA <- c1otu %>%
  mutate(presence = 1) %>%
  pivot_wider(names_from = sequence, 
              values_from = presence, 
              values_fill = 0) %>%
  as.data.frame()

rownames(c1otu_PA) <- c1otu_PA$pathway

k_c1 <- c1otu_PA[, -1]
k2_c1 <- t(k_c1)


# species accumulation calculation
c1_accum_curve <- specaccum(k2_c1, method = "random")

plot(c1_accum_curve, col = "blue",
     xlab = "Number of pathways",
     ylab = "Number of otus")

c1_curve_mc <- c1_accum_curve[["richness"]]
c1_curve_otu <- c1_accum_curve[["sites"]]

x_c1 <- cbind.data.frame(c1_curve_otu, c1_curve_mc)

ggplot(x_c1, aes(x = c1_curve_otu, y = c1_curve_mc, color = "blue")) + geom_line()

c1_k2_taxa <- as.data.frame(rownames(k2_c1))
names(c1_k2_taxa) <- "otu"

c1_k2_taxa2 <- merge(c1_k2_taxa, taxa_fam, by = "otu")

c1_t2 <- c1_k2_taxa2 %>%
  sample_n(n()) %>%
  mutate(observation = row_number()) %>%
  arrange(family, observation) %>%
  group_by(family) %>%
  mutate(distinct = row_number() == 1) %>%
  ungroup() %>%
  arrange(observation) %>%
  mutate(s = cumsum(distinct))

c1_x2 <- cbind.data.frame(x_c1, c1_t2$s)

zero <- c(0, 0, 0, 0)

c1_x3 <- rbind.data.frame(zero, c1_x2)

c1_x3$taxa <- as.numeric(c1_x3$`c1_t2$s`)

# Cellulose diet Week 1 - Collector Curve | Y axis is Percent
ggplot(c1_x3, aes(x = c1_curve_otu)) + geom_line(aes(y = ((taxa/102) * 100), color = "family"), color = "red") +
  geom_line(aes(y = (c1_curve_mc/337 * 100), color = "pathway"), color = "blue") + 
  labs(title = "Cellulose Week 1", 
       x = "OTU Count 460", 
       y = "Percent") + theme_classic()

# Cellulose diet Week 1 - Collector Curve | Y axis is Count
ggplot(c1_x3, aes(x = c1_curve_otu)) + geom_line(aes(y = (taxa), color = "family"), color = "red") +
  geom_line(aes(y = (c1_curve_mc), color = "pathway"), color = "blue") + 
  labs(title = "Cellulose Week 1", 
       x = "OTU Count 460", 
       y = "Count") + theme_classic()

# remove counts
cellulose_week4_otumc <- cellulose_week4_count_rm0[, c(2, 1)]

# group pathways by otu

c4otu <- cellulose_week4_otumc %>%
  arrange(sequence) %>%
  as.data.frame()

c4otu_PA <- c4otu %>%
  mutate(presence = 1) %>%
  pivot_wider(names_from = sequence, 
              values_from = presence, 
              values_fill = 0) %>%
  as.data.frame()

rownames(c4otu_PA) <- c4otu_PA$pathway

k_c4 <- c4otu_PA[, -1]
k2_c4 <- t(k_c4)


# species accumulation calculation
c4_accum_curve <- specaccum(k2_c4, method = "random")

plot(c4_accum_curve, col = "blue",
     xlab = "Number of pathways",
     ylab = "Number of otus")

c4_curve_mc <- c4_accum_curve[["richness"]]
c4_curve_otu <- c4_accum_curve[["sites"]]

x_c4 <- cbind.data.frame(c4_curve_otu, c4_curve_mc)

ggplot(x_c4, aes(x = c4_curve_otu, y = c4_curve_mc, color = "blue")) + geom_line()

c4_k2_taxa <- as.data.frame(rownames(k2_c4))
names(c4_k2_taxa) <- "otu"

c4_k2_taxa2 <- merge(c4_k2_taxa, taxa_fam, by = "otu")

c4_t2 <- c4_k2_taxa2 %>%
  sample_n(n()) %>%
  mutate(observation = row_number()) %>%
  arrange(family, observation) %>%
  group_by(family) %>%
  mutate(distinct = row_number() == 1) %>%
  ungroup() %>%
  arrange(observation) %>%
  mutate(s = cumsum(distinct))

c4_x2 <- cbind.data.frame(x_c4, c4_t2$s)

zero <- c(0, 0, 0, 0)

c4_x3 <- rbind.data.frame(zero, c4_x2)

c4_x3$taxa <- as.numeric(c4_x3$`c4_t2$s`)

# Cellulose diet Week 4 - Collector Curve | Y axis is Percent
ggplot(c4_x3, aes(x = c4_curve_otu)) + geom_line(aes(y = ((taxa/91) * 100), color = "family"), color = "red") +
  geom_line(aes(y = (c4_curve_mc/337 * 100), color = "pathway"), color = "blue") + 
  labs(title = "Cellulose Week 4", 
       x = "OTU Count (396)", 
       y = "Percent") + theme_classic()

# Cellulose diet Week 4 - Collector Curve | Y axis is Percent
ggplot(c4_x3, aes(x = c4_curve_otu)) + geom_line(aes(y = (taxa), color = "family"), color = "red") +
  geom_line(aes(y = (c4_curve_mc), color = "pathway"), color = "blue") + 
  labs(title = "Cellulose Week 4", 
       x = "OTU Count (396)", 
       y = "Percent") + theme_classic()


