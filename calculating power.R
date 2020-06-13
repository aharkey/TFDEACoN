### A Harkey
### 6-12-2020
### Calculating the sample size needed for TF DEACoN to detect enrichment

library(statmod)

# The TFs in the DAP-Seq dataset have 199-19,976 targets
# I will calculate the power for both extremes, and the average

# The known ratios to test against
tot <- 27655
maxratio <- 19976/tot
minratio <- 199/tot
averatio <- 7362/tot

# The effect size to look for
logfccut <- 0.5
effectsize <- 2^logfccut

maxeffect <- maxratio * effectsize
# maxeffect is > 1. Set to 0.8
maxeffect <- 0.8
mineffect <- minratio * effectsize
aveeffect <- averatio * effectsize

# The power cutoff
power <- 0.80

# The p-value cutoff
pval <- 0.05

# Max calculation
maxpower <- data.frame(n = c(seq(100, 1000, 100), seq(200, 300, 10)), power = as.numeric(NA))
for (i in 1:nrow(maxpower))
{
  maxpower$power[i] <- power.fisher.test(maxeffect, maxratio, maxpower$n[i], 27655, alpha = pval)
}

# You need an n of 240 to detect significance for the TF with the most targets

# Min calculation
minpower <- data.frame(n = seq(1000, 10000, 1000), power = as.numeric(NA))
for (i in 1:nrow(minpower))
{
  minpower$power[i] <- power.fisher.test(mineffect, minratio, minpower$n[i], 27655, alpha = pval, nsim = 500)
}

# You need an n of about 10,000 to detect significance for the TF with the least targets

# Ave calculation
avepower <- data.frame(n = seq(130, 160, 2), power = as.numeric(NA))
for (i in 1:nrow(avepower))
{
  avepower$power[i] <- power.fisher.test(aveeffect, averatio, avepower$n[i], 27655, alpha = pval, nsim = 500)
}

# You need an n of about 144 to detect significance for the TF with the least targets