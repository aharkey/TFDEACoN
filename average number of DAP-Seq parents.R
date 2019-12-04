setwd("C:/Users/aharkey/Desktop/DAP-Seq reverse search")

alltarg <- read.csv("input/alltarg.txt",
                    header = TRUE,
                    row.names = 1)

uniqtarg <- data.frame(as.character(unique(alltarg$target)))
colnames(uniqtarg) <- "target"

uniqtarg$TFcount <- NA

for (i in 4001:5000)
{
  uniqtarg$TFcount[i] <- length(which(alltarg$target == uniqtarg$target[i]))
}

length(unique(alltarg$TF))

# Each gene has an average of
# 2848928/32605=
# 87.38 TFs bound to it
# out of 387 TFs
# So, 22.58% of all DAP-Seq TFs bind to any given gene

# Each TF binds to an average of
# 2848928/387=
# 7361.57 targets
# out of 32605 genes
# So, 22.58 % of all targets are bound by any given TF




