args                <- commandArgs(T)
input.parameters.fn <- as.character(args[1])

wd <- "/home/ahryhorzhevska/mpip/code/dex-stim-human-dna-methyl-qc/"
setwd(wd)

input.parameters    <- read.csv("input_parameters.csv", sep = ";", header = F )
input.parameters    <- as.data.frame(input.parameters)

for (row in 1:nrow(input.parameters))
  assign(input.parameters[row, 1], input.parameters[row, 2])

# Load utils

source("packages.R")
source("03_normalization/01_help_files/BMIQ_1.6_Teschendorff.R")

# Load RGSet
src.data.dir.pre <- "/binder/mgp/datasets/2020_DexStim_Array_Human/methylation/10_final_qc_data/"
rg.set.fn        <- paste0(src.data.dir.pre , "dex_methyl_rgset_final.rds")
rg.set           <- readRDS(rg.set.fn)

# Extract beta matrix and annotate RGSet
# annot     <- getAnnotation(rg.set)
beta.mtrx <- getBeta(rg.set)
beta.mtrx <- na.omit(beta.mtrx)

mset      <- preprocessRaw(rg.set)
manifest  <- getManifest(mset)
design    <- ifelse(rownames(beta.mtrx) %in% manifest@data$TypeI$Name, "1", "2")

# BMIQ Normalization

bmiq.beta.mtrx      <- apply(beta.mtrx[, 1:length(colnames(beta.mtrx))], 2, 
				             function(a) BMIQ(a, design, plots = FALSE)$nbeta)

rownames(bmiq.beta.mtrx) <- rownames(beta.mtrx)
colnames(bmiq.beta.mtrx) <- colnames(beta.mtrx)

# probeType           <- annot[rownames(beta.mtrx), c("Name","Type")]
# probeType           <- na.omit(probeType)
# probeType           <- as.data.frame(probeType)  
# probeType$probeType <- ifelse(probeType$Type %in% "I", 1, 2)

# bmiq.beta.mtrx      <- apply(beta.mtrx[, 1:length(colnames(beta.mtrx))], 2, 
# 				             function(a) BMIQ(a,probeType$probeType,plots = FALSE)$nbeta)

length(which(is.nan(bmiq.beta.mtrx))) # should be 0
saveRDS(bmiq.beta.mtrx, 
        file = paste0(src.data.dir.pre, "dex_methyl_bmiq_beta_mtrx.rds"))
