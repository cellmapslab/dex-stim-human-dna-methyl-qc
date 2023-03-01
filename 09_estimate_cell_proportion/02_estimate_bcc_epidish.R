# 09 Cell Composition

# 1. Set up env and load data

BiocManager::install("EpiDISH")

library(EpiDISH)
library(data.table)
library(dplyr)
library(ggplot2)
library(RColorBrewer)

beta.mtrx.fn <- "~/bio/code/mpip/dex-stim-human-array/data/methylation/dex_methyl_beta_combat_mtrx.rds"
beta.mtrx    <- readRDS(beta.mtrx.fn)

pheno.fn     <- "~/bio/code/mpip/dex-stim-human-array/data/pheno/pheno_full_for_kimono.csv"
pheno        <- fread(pheno.fn)
veh.dnam.ids <- pheno[Dex == 0][!is.na(DNAm_ID)][,DNAm_ID]

beta.mtrx.veh <- beta.mtrx[, colnames(beta.mtrx) %in% veh.dnam.ids]

# 2. Estimate cell type proportion

bcc.epidish.rpc.mdl <- epidish(beta.m = beta.mtrx, ref.m = centDHSbloodDMC.m, method = "RPC")
bcc.epidish.rpc     <- bcc.epidish.rpc.mdl$estF %>% data.frame() %>% setDT(keep.rownames = T)

bcc.epidish.cbs.mdl <- epidish(beta.m = beta.mtrx, ref.m = centDHSbloodDMC.m, method = "CBS")
bcc.epidish.cbs     <- bcc.epidish.cbs.mdl$estF %>% data.frame() %>% setDT(keep.rownames = T)

bcc.epidish.cp.mdl <- epidish(beta.m = beta.mtrx, ref.m = centDHSbloodDMC.m, method = "CP")
bcc.epidish.cp     <- bcc.epidish.cp.mdl$estF %>% data.frame() %>% setDT(keep.rownames = T)

head(bcc.epidish.rpc)
head(bcc.epidish.cbs)
head(bcc.epidish.cp)

bcc.epidish <- data.frame(epidish_rpc = bcc.epidish.rpc, epidish_cbs = bcc.epidish.cbs[, -1], epidish_cp = bcc.epidish.cp[, -1])
colnames(bcc.epidish)[1] <- "DNAm_ID"

write.csv2(bcc.epidish, 
           file = "~/bio/code/mpip/dex-stim-human-array/output/data/methylation/dex_stim_array_human_epidish_bcc.csv", 
           row.names = F, quote = F)

# 3. Estimate cell type proportion: new ref data, 02.05.2022

# the new ref blood data set
load("~/bio/code/mpip/dex-stim-human-array/data/public_data/FlowSorted.BloodExtended.EPIC.compTable.rda")

bcc.epidish.rpc.mdl <- epidish(
  beta.m = beta.mtrx.veh,
  ref.m = FlowSorted.BloodExtended.EPIC.compTable,
  method = "RPC"
)

bcc.epidish.rpc              <- bcc.epidish.rpc.mdl$estF %>% data.frame() %>% setDT(keep.rownames = T)
colnames(bcc.epidish.rpc)[1] <- "DNAm_ID"

# store the results
write.csv2(bcc.epidish.rpc, 
           file = "~/bio/code/mpip/dex-stim-human-array/output/data/methylation/dex_stim_array_human_epidish_salas_bcc_rpc.csv", 
           row.names = F, quote = F)

# 
# summary(bcc.epidish)
# cor.mtrx <- cor(bcc.epidish[, -1])
# 
# bcc.epidish.by.meth <- rbind(bcc.epidish.rpc[,method := "RPC"],
#                              bcc.epidish.cbs[,method := "CBS"],
#                              bcc.epidish.cp[,method := "CP"])
# 
# dat.m <- melt(bcc.epidish.by.meth, id.vars = c('rn', 'method'))
# 
# ggplot(dat.m, aes(x = variable, y = value, fill =  method)) +
#   geom_boxplot() 
# 
