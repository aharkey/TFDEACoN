### A Harkey
### 6-10-2020

# This script is intended to take data from datasets which have been clustered and create a data frame that contains unique gene IDs and the clusters to which they belong

##### Brady et al. dataset #####
# Root cell type clustering
# These are the clusters shown as a supplemental figure in the Harkey et al. (2020) TF DEACoN manuscript

dat <- read.csv("G:/My Drive/WFU/01 Research/06 TF DEACoN/Code/input/clusterdatasets/bradyclusters.csv", row.names = 1)
