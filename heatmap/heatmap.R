# Heatmap of all bacterial families for frass paper
library(vegan)
library(ComplexHeatmap)
library(circlize)
library(stringr)
library(ggplot2)
library(scales)
library(readxl)

clean_data <- read_excel("~/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Frass_Paper/Sabree_L5/KC_L5_Filter_Data_20210324.xlsx")
clean_data <- clean_data[, -1]
# remove archaea 
clean_data <- clean_data[-c(1:3), ]

redata <- t(clean_data)
colnames(redata) <- redata[1, ]
redata <- redata[-1, ]

order <- as.data.frame(redata[order(row.names(redata)), ])

split_order <- str_split_fixed(rownames(order), '_', 3)

colnames(split_order) <- c("Diet", "Ind", "Week")

split_order <- cbind.data.frame(split_order, order)
order_split_order <- split_order[order(split_order$Week), ]

numeric <- as.data.frame(apply(order_split_order[, -c(1:3)], 2, as.numeric))
rownames(numeric) <- rownames(order_split_order)

order_split_order2 <- cbind.data.frame(order_split_order[, 1:3], numeric)

# separate by diet
DF.2 <- order_split_order2[grep("DF", order_split_order2$Diet), ]
HN.2 <- order_split_order2[grep("HN", order_split_order2$Diet), ]
LN.2 <- order_split_order2[grep("LN", order_split_order2$Diet), ]

# separate by week
### DF week
DF_Week1 <- as.data.frame(DF.2[which(DF.2$Week == '1'), 4:147])
DF_Week2 <- as.data.frame(DF.2[which(DF.2$Week == '2'), 4:147])
DF_Week3 <- as.data.frame(DF.2[which(DF.2$Week == '3'), 4:147])
DF_Week4 <- as.data.frame(DF.2[which(DF.2$Week == '4'), 4:147])
DF_Week5 <- as.data.frame(DF.2[which(DF.2$Week == '5'), 4:147])
DF_Week6 <- as.data.frame(DF.2[which(DF.2$Week == '6'), 4:147])
DF_Week7 <- as.data.frame(DF.2[which(DF.2$Week == '7'), 4:147])
DF_Week8 <- as.data.frame(DF.2[which(DF.2$Week == '8'), 4:147])

### HN week
HN_Week1 <- as.data.frame(HN.2[which(HN.2$Week == '1'), 4:147])
HN_Week2 <- as.data.frame(HN.2[which(HN.2$Week == '2'), 4:147])
HN_Week3 <- as.data.frame(HN.2[which(HN.2$Week == '3'), 4:147])
HN_Week4 <- as.data.frame(HN.2[which(HN.2$Week == '4'), 4:147])
HN_Week5 <- as.data.frame(HN.2[which(HN.2$Week == '5'), 4:147])
HN_Week6 <- as.data.frame(HN.2[which(HN.2$Week == '6'), 4:147])
HN_Week7 <- as.data.frame(HN.2[which(HN.2$Week == '7'), 4:147])
HN_Week8 <- as.data.frame(HN.2[which(HN.2$Week == '8'), 4:147])

### LN week
LN_Week1 <- as.data.frame(LN.2[which(LN.2$Week == '1'), 4:147])
LN_Week2 <- as.data.frame(LN.2[which(LN.2$Week == '2'), 4:147])
LN_Week3 <- as.data.frame(LN.2[which(LN.2$Week == '3'), 4:147])
LN_Week4 <- as.data.frame(LN.2[which(LN.2$Week == '4'), 4:147])


# Convert to Relative Abundance with Vegan Package
apply(DF_Week1, 1, sum)
DF_1_rel <- decostand(DF_Week1, method = "total")
apply(DF_1_rel, 1, sum)

# 2
apply(DF_Week2, 1, sum)
DF_2_rel <- decostand(DF_Week2, method = "total")
apply(DF_2_rel, 1, sum)
# 3
apply(DF_Week3, 1, sum)
DF_3_rel <- decostand(DF_Week3, method = "total")
apply(DF_3_rel, 1, sum)
# 4
apply(DF_Week4, 1, sum)
DF_4_rel <- decostand(DF_Week4, method = "total")
apply(DF_4_rel, 1, sum)
# 5
apply(DF_Week5, 1, sum)
DF_5_rel <- decostand(DF_Week5, method = "total")
apply(DF_5_rel, 1, sum)
# 6
apply(DF_Week6, 1, sum)
DF_6_rel <- decostand(DF_Week6, method = "total")
apply(DF_6_rel, 1, sum)
# 7
apply(DF_Week7, 1, sum)
DF_7_rel <- decostand(DF_Week7, method = "total")
apply(DF_7_rel, 1, sum)

# 8
apply(DF_Week8, 1, sum)
DF_8_rel <- decostand(DF_Week8, method = "total")
apply(DF_8_rel, 1, sum)

# HN:convert to Relative Abundance with vegan package
# 1
apply(HN_Week1, 1, sum)
HN_1_rel <- decostand(HN_Week1, method = "total")
apply(HN_1_rel, 1, sum)
# 2
apply(HN_Week2, 1, sum)
HN_2_rel <- decostand(HN_Week2, method = "total")
apply(HN_2_rel, 1, sum)
# 3
apply(HN_Week3, 1, sum)
HN_3_rel <- decostand(HN_Week3, method = "total")
apply(HN_3_rel, 1, sum)
# 4
apply(HN_Week4, 1, sum)
HN_4_rel <- decostand(HN_Week4, method = "total")
apply(HN_4_rel, 1, sum)
# 5
apply(HN_Week5, 1, sum)
HN_5_rel <- decostand(HN_Week5, method = "total")
apply(HN_5_rel, 1, sum)
# 6
apply(HN_Week6, 1, sum)
HN_6_rel <- decostand(HN_Week6, method = "total")
apply(HN_6_rel, 1, sum)

# 7
apply(HN_Week7, 1, sum)
HN_7_rel <- decostand(HN_Week7, method = "total")
apply(HN_7_rel, 1, sum)

# 8
apply(HN_Week8, 1, sum)
HN_8_rel <- decostand(HN_Week8, method = "total")
apply(HN_8_rel, 1, sum)


# LN:convert to Relative Abundance with vegan package
# 1
apply(LN_Week1, 1, sum)
LN_1_rel <- decostand(LN_Week1, method = "total")
apply(LN_1_rel, 1, sum)
# 2
apply(LN_Week2, 1, sum)
LN_2_rel <- decostand(LN_Week2, method = "total")
apply(LN_2_rel, 1, sum)
# 3
apply(LN_Week3, 1, sum)
LN_3_rel <- decostand(LN_Week3, method = "total")
apply(LN_3_rel, 1, sum)
# 4
apply(LN_Week4, 1, sum)
LN_4_rel <- decostand(LN_Week4, method = "total")
apply(LN_4_rel, 1, sum)


# DF: Means of relative abundance
# 1
DF_1_rel_means <- colMeans(DF_1_rel)
# 2
DF_2_rel_means <- colMeans(DF_2_rel)
# 3
DF_3_rel_means <- colMeans(DF_3_rel)
# 4
DF_4_rel_means <- colMeans(DF_4_rel)
# 5
DF_5_rel_means <- colMeans(DF_5_rel)
# 6
DF_6_rel_means <- colMeans(DF_6_rel)
# 7
DF_7_rel_means <- colMeans(DF_7_rel)
# 8
DF_8_rel_means <- colMeans(DF_8_rel)

DF_1_rel_means -> DF_Week1
DF_2_rel_means -> DF_Week2
DF_3_rel_means -> DF_Week3
DF_4_rel_means -> DF_Week4
DF_5_rel_means -> DF_Week5
DF_6_rel_means -> DF_Week6
DF_7_rel_means -> DF_Week7
DF_8_rel_means -> DF_Week8

# HN: Means of relative abundance
# 1
HN_1_rel_means <- colMeans(HN_1_rel)
# 2
HN_2_rel_means <- colMeans(HN_2_rel)
# 3
HN_3_rel_means <- colMeans(HN_3_rel)
# 4
HN_4_rel_means <- colMeans(HN_4_rel)
# 5
HN_5_rel_means <- colMeans(HN_5_rel)
# 6
HN_6_rel_means <- colMeans(HN_6_rel)
# 7
HN_7_rel_means <- colMeans(HN_7_rel)
# 8 
HN_8_rel_means <- colMeans(HN_8_rel)

HN_1_rel_means -> HN_Week1
HN_2_rel_means -> HN_Week2
HN_3_rel_means -> HN_Week3
HN_4_rel_means -> HN_Week4
HN_5_rel_means -> HN_Week5
HN_6_rel_means -> HN_Week6
HN_7_rel_means -> HN_Week7
HN_8_rel_means -> HN_Week8

# LN: Means of relative abundance
# 1
LN_1_rel_means <- colMeans(LN_1_rel)
# 2
LN_2_rel_means <- colMeans(LN_2_rel)
# 3
LN_3_rel_means <- colMeans(LN_3_rel)
#4
LN_4_rel_means <- colMeans(LN_4_rel)
LN_1_rel_means -> LN_Week1
LN_2_rel_means -> LN_Week2
LN_3_rel_means -> LN_Week3
LN_4_rel_means -> LN_Week4

# combine all DF, HN, LN rel means into a table 
DF_rbind <- rbind(DF_Week1, DF_Week2, DF_Week3, DF_Week4, DF_Week5, DF_Week6, 
                  DF_Week7, DF_Week8)
HN_rbind <- rbind(HN_Week1, HN_Week2, HN_Week3, HN_Week4, HN_Week5, HN_Week6, 
                  HN_Week7, HN_Week8)
LN_rbind <- rbind(LN_Week1, LN_Week2, LN_Week3, LN_Week4)

DF_names <- as.data.frame(colnames(DF_rbind))
DF_names_split <- str_split_fixed(DF_names$`colnames(DF_rbind)`, ';', 5)
DF_rownames <- paste(DF_names_split[,3], DF_names_split[, 4], DF_names_split[,5], sep = " ")

D = t(DF_rbind)
D_bind <- cbind.data.frame(DF_rownames, D)
D_bind <- D_bind[, -1]
D_bind[D_bind < 0.001] <- NA

HN_names <- as.data.frame(colnames(HN_rbind))
HN_names_split <- str_split_fixed(HN_names$`colnames(HN_rbind)`, ';', 5)

H = t(HN_rbind)
H_bind <- cbind.data.frame(HN_names_split[,5], H)
H_bind <- H_bind[,-1]
H_bind[H_bind < 0.001] <- NA

LN_names <- as.data.frame(colnames(LN_rbind))
LN_names_split <- str_split_fixed(LN_names$`colnames(LN_rbind)`, ';', 5)

L = t(LN_rbind)
L_bind <- cbind.data.frame(LN_names_split[,5], L)
L_bind <- L_bind[,-1]
L_bind[L_bind < 0.001] <- NA

D_bind <- as.matrix(D_bind)
H_bind <- as.matrix(H_bind)
L_bind <- as.matrix(L_bind)

col_fun = colorRamp2(c(0.0025, 0.15, 0.3), c("gold","red", "black"))

DF_H = Heatmap(D_bind, name = "DF", col = col_fun, na_col = "white", cluster_rows = FALSE, 
               cluster_columns = FALSE)
HN_H = Heatmap(H_bind, name = "H;.N", col = col_fun, na_col = "white", cluster_rows = FALSE, 
               cluster_columns = FALSE)
LN_H = Heatmap(L_bind, name = "LN", col = col_fun, na_col = "white", cluster_rows = FALSE, 
               cluster_columns = FALSE)
Ht_List = DF_H + HN_H + LN_H

lgd = Legend(title = "Relative Abundance", col_fun = col_fun, 
             at = c(0.0025, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3), labels = c("0.0025", "0.05", "0.1", "0.15", "0.2", "0.25", "0.3"), 
             direction = "horizontal")

draw(Ht_List, annotation_legend_list = lgd, annotation_legend_side = "top", 
     show_heatmap_legend = FALSE)

# double check on why there is 4:147 families here but in OTU pathway count there was 105 ...
# this was because families in OTU pathway count combined all families that were "other" or the same at the family level
