### A Harkey
### 6-10-2020

# This script is intended to take data from datasets which have been clustered and create a data frame that contains unique gene IDs and the clusters to which they belong

library(tidyr)

##### Brady et al. dataset #####
# Root cell type clustering
# These are the clusters shown as a supplemental figure in the Harkey et al. (2020) TF DEACoN manuscript

dat <- read.csv("G:/My Drive/WFU/01 Research/06 TF DEACoN/Code/input/clusterdatasets/bradyclusters.csv", row.names = 1)
dat$locus <- as.character(dat$locus)
dat <- dat[, -which(colnames(dat) %in% c("cluster", "order"))]
colnames(dat)[which(colnames(dat) == "newclust")] <- "cluster"

# (max(nchar(dat$locus)) + 1) / 10
# The locus with the most gene IDs has 15 of them

dat.sing <- dat[which(nchar(dat$locus) == 9),]
dat.multi <- dat[-which(nchar(dat$locus) == 9),]

dat.multi <- separate(dat.multi, "locus",  paste0("locus", 1:15))

dat.split <- NA

for (i in 1:15)
{
  
}
