# ALDEx2 on Balanced v Protein weeks 3 and 4 
# Balanced v Cellulose weeks 3 and 4
# Comparison between balance and treatment diet used weeks three and four due to 
# (a) low sample size in the cellulose enriched diet during weeks five to eight and 
# (b) provide time for bacterial taxa to respond to the shift in diet type from balanced to week four of the nutrient-poor diet

library(ALDEx2)

MC = read_tsv("/Users/kaylacross/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Frass_Paper/Sabree_L7/PICRUST2_all/picrust2_out_pipeline2/pathways_out/path_abun_unstrat.tsv")
MC = as.data.frame(MC)

MC2 = MC

rownames(MC2) = MC2$pathway
MC2 = MC2[, -1]

MC2_t = t(MC2)
MC2_t = MC2_t[order(row.names(MC2_t)), ]

# seperate out diet, ind, and week
library(stringr)
diet <- str_split_fixed(rownames(MC2_t), "_", 3)
colnames(diet) <- c("diet", "ind", "week")
diet_week <- paste(diet[,1], diet[,3], sep = "_")
diet_week <- as.data.frame(diet_week)
diet2 <- cbind.data.frame(diet_week, diet)

pathway <- cbind.data.frame(diet_week, MC2_t)

DF_1 <- pathway[which(pathway$diet_week == "DF_1"), ]
DF_2 <- pathway[which(pathway$diet_week == "DF_2"), ]
DF_3 <- pathway[which(pathway$diet_week == "DF_3"), ]
DF_4 <- pathway[which(pathway$diet_week == "DF_4"), ]

HN_1 <- pathway[which(pathway$diet_week == "HN_1"), ]
HN_2 <- pathway[which(pathway$diet_week == "HN_2"), ]
HN_3 <- pathway[which(pathway$diet_week == "HN_3"), ]
HN_4 <- pathway[which(pathway$diet_week == "HN_4"), ]

DF_HN_3_4 <- rbind.data.frame(DF_3, DF_4, 
                              HN_3, HN_4)

data <- as.data.frame(t(DF_HN_3_4[, -1]))

conds.hn <- c(rep("DF", 13), rep("HN", 12))


df.hn.3.4 = aldex(round(data), conds.hn, effects = T, denom = "iqlr")
Wi_Effect = df.hn.3.4[which(df.hn.3.4$wi.eBH < 0.05 & abs(df.hn.3.4$effect) >= 0.8), ] 
# note: The number of significant features will be different each time aldex2 is ran
# this is due to the sampling procedure. See section 5.1 of the aldex2 vingette for
# more information


LN_1 <- pathway[which(pathway$diet_week == "LN_1"), ]
LN_2 <- pathway[which(pathway$diet_week == "LN_2"), ]
LN_3 <- pathway[which(pathway$diet_week == "LN_3"), ]
LN_4 <- pathway[which(pathway$diet_week == "LN_4"), ]

DF_LN_3_4 <- rbind.data.frame(DF_3, DF_4, 
                              LN_3, LN_4)

data3 <- as.data.frame(t(DF_LN_3_4[, -1]))

conds.ln <- c(rep("DF", 13), rep("LN", 10))


df.ln.3.4 = aldex(round(data3), conds.ln, effects = T, denom = "iqlr")
Wi_Effect.ln = df.ln.3.4[which(df.ln.3.4$wi.eBH < 0.05 & abs(df.ln.3.4$effect) >= 0.8), ] 
# note: The number of significant features will be different each time aldex2 is ran
# this is due to the sampling procedure. See section 5.1 of the aldex2 vingette for
# more information
