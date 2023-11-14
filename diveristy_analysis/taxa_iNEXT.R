library(knitr)
library(iNEXT)
library(ggplot2)
library(cowplot)
library(xtable)

location="/users/PAS0688/cross582/check"
file.name=c("KC_L5_Filter_Data_20210324.txt")

data=read.table(file.path(location,file.name),header=TRUE,sep="\t")

# Remove first column (first column are phylotypes names)
redata=data[,-1]
redata=redata[,-1]

order <- as.data.frame(t(redata))

order <- order[order(rownames(order)), ]

x <- order[-c(109, 108, 114:117), ]

y <- as.data.frame(t(x))

# Hill Numbers: q - 0, 1, 2
out0_new <- iNEXT(y,q=0,datatype="abundance", knots = 10000, nboot = 200)
save(out0_new,file=file.path(location,paste0("diversity-order0-",format(Sys.Date(),"%Y%m%d"),"family-knots10000.Rdata")))


diversity_min_size = estimateD(y,datatype="abundance", base = "size", level = min(out0_new$DataInfo$n)*2)
save(diversity_min_size,file=file.path(location,paste0("diversity-pointestimate-minsize",format(Sys.Date(),"%Y%m%d"),"family-knots10000.Rdata")))

diversity_max_size = estimateD(y,datatype="abundance", base = "size", level = max(out0_new$DataInfo$n))
save(diversity_max_size,file=file.path(location,paste0("diversity-pointestimate-maxsize",format(Sys.Date(),"%Y%m%d"),"family-knots10000.Rdata")))

diversity_min_coverage = estimateD(y,datatype="abundance", base = "coverage")
save(diversity_min_coverage,file=file.path(location,paste0("diversity-pointestimate-mincoverage",format(Sys.Date(),"%Y%m%d"),"family-knots10000.Rdata")))
