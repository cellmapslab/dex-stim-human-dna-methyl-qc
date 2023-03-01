# 06 Surrogate Variable Analysis (SVA)

# Reference: https://bioconductor.org/packages/release/bioc/vignettes/sva/inst/doc/sva.pdf

setwd("~/bio/code/mpip/dex-stim-human-dna-methyl-qc/")

library(sva)
library(lumi) # convert beta values to M-values 

# input.parameters.fn <- "input_parameters.csv"
# input.parameters    <- read.csv(input.parameters.fn, sep = ";", header = F )
# input.parameters    <- as.data.frame(input.parameters)
# 
# for (row in 1:nrow(input.parameters))
#   assign(as.character(input.parameters[row, 1]), input.parameters[row, 2])

# 1. Setting up the data

# m.mtrx.fn <- "~/bio/datasets/methylation/10_final_qc_data/dex_methyl_mval_mtrx_final.rds"
# pheno.fn <- "~/bio/datasets/methylation/10_final_qc_data/dex_methyl_phenotype.rds"

m.mtrx.fn <- "/binder/mgp/datasets/2020_DexStim_Array_Human/methylation/10_final_qc_data/dex_methyl_mval_mtrx_final.rds"
pheno.fn  <- "/binder/mgp/datasets/2020_DexStim_Array_Human/methylation/10_final_qc_data/dex_methyl_phenotype.rds"

m.mtrx    <- readRDS(m.mtrx.fn)
pheno     <- readRDS(pheno.fn)
pheno     <- pheno[pheno$Sample_Name %in% colnames(m.mtrx),]

#2. Convert beta values to M-values

# m.mtrx <- beta2m(beta.mtrx)

# 3. Create the full and null model matrices

pheno$person       <- factor(pheno$person)
pheno$Sample_Group <- factor(pheno$Sample_Group)
pheno$status       <- factor(pheno$status)

# Model 1

# mod0 <- model.matrix(~ 1, data = pheno)
# mod  <- model.matrix(~ Sample_Group, data = pheno)

# Model 2

mod0 <- model.matrix(~ 1 + person, data = pheno)
mod  <- model.matrix(~ person + Sample_Group, data = pheno)

# Model 3
# mod0 <- model.matrix(~ 1 + person, data = pheno)
# mod  <- model.matrix(~ person + Sample_Group:status, data = pheno)

# Model 4

# mod0 <- model.matrix(~ 1 + person, data = pheno)
# mod  <- model.matrix(~ person + Sample_Group:status, data = pheno)
# mod  <- mod[,!grepl("status0", colnames(mod))]

# 4. Applying the 'sva' to estimate surrogate variables
# n.sv  <- num.sv(m.mtrx, mod, method = "leek")
svobj <- sva(m.mtrx, mod, mod0) #, n.sv = n.sv)

# Save the results
methyl.svs.df <- as.data.frame(cbind(pheno$Sample_Name, svobj$sv))
colnames(methyl.svs.df) <- c("DNAm_ID", paste0("DNAm_SV" , 1:ncol(svobj$sv)))

write.csv2(methyl.svs.df, 
           file = paste0("/binder/mgp/datasets/2020_DexStim_Array_Human/methylation/15_DNAm_sva/", "methyl_svs.csv"),
           row.names = F, quote = F)

