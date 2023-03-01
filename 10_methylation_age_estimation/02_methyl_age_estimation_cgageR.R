# 10 Methylation age estimation

# 1. Set up env and load data

library(minfi)
library(cgageR)
library(SummarizedExperiment)
library(ggplot2)
library(tidyverse)

input.parameters.fn <- "input_parameters.csv"

input.parameters    <- read.csv(input.parameters.fn, sep = ";", header = F )
input.parameters    <- as.data.frame(input.parameters)


for (row in 1:nrow(input.parameters))
  assign(input.parameters[row, 1], input.parameters[row, 2])

src.final.data.dir <- "~/bio/datasets/methylation/10_final_qc_data/"

beta.mtrx.fn <- paste0(src.final.data.dir, "dex_methyl_qn_beta_mtrx.rds") # "dex_methyl_beta_combat_mtrx.rds")
beta.mtrx    <- readRDS(beta.mtrx.fn)

pheno.fn     <- paste0(src.final.data.dir, "dex_methyl_phenotype.rds")
pheno        <- readRDS(pheno.fn)

pheno$Status  <- as.factor(pheno$Status)
pheno$Sex     <- as.factor(pheno$Sex)

levels(pheno$Sex) <- c("Male", "Female")
levels(pheno$Status) <- c("Control", "MDD")

# 2. Estimate methylation age

meth.age.obj  <- getAgeR(beta.mtrx, epitoc=FALSE, horvath=T, hannum=FALSE, drift=FALSE, driftcg=NULL, chrage=NULL, 
                         showStatusHannum=FALSE, keepres=FALSE,keepcpgs.hannum=TRUE,
                         keepcpgs.epitoc=TRUE, keepcpgs.horvath=TRUE)

meth.age <- meth.age.obj$HorvathClock.output$Horvath.Est
meth.age[["SampleID"]] <- rownames(meth.age)

# 3. Merge estimated age with real and write out results

real.age.tbl   <- data.frame(SampleID = pheno$DNAm_ID, realAge = pheno$Age, mAge_Hovath = pheno$mAge_Hovath, Individual = pheno$Sample_ID, 
                             Dex = pheno$Dex, Sex = pheno$Sex, Status = pheno$Status)
merged.age.tbl <- merge(meth.age, real.age.tbl)
merged.age.tbl[["DeltaDNAmAge"]] <- merged.age.tbl$Horvath.Est - merged.age.tbl$mAge_Hovath
merged.age.tbl <- merged.age.tbl[merged.age.tbl$Dex == 0, ] %>% setDT()

pheno.full.new <- left_join(pheno.full, merged.age.tbl[, .(Individual, Horvath.Est)], by = c("Sample_ID" = "Individual"))
pheno.full.new[["mAge_Hovath"]] <- pheno.full.new[["Horvath.Est"]]
pheno.full.new <- pheno.full.new[, -155]
# write.table(merged.age.tbl, file = paste0(report.dir, "DEX-stim-array-human-meth-age.csv"), col.names = T, row.names = T, quote = F, sep = ";")

# 4. Relate DNAm age to chronological age

MedianAbsDevHelp <- function(x, y) median(abs(x-y), na.rm=TRUE)
MedianAbsDev     <- function(dnam.age.method, real.age) signif(MedianAbsDevHelp(dnam.age.method, real.age), 2)


pdf(file = paste0(report.dir, "DNAm_Age_and_Chronological_Age_Relation_cgageR.pdf"))

err <- MedianAbsDev (merged.age.tbl$Horvath.Est, merged.age.tbl$realAge)
ggplot(merged.age.tbl, aes(x = mAge_Hovath, y = realAge,  color = pheno$Sample_Group)) + 
  geom_point() +
  geom_smooth() +
  labs(title = paste("DNAm age (Hovath) vs Chronological age, err = ", err),
       x = "DNAm Age Hovath", y = "Chronological Age")
dev.off()

# 5. DNAm age Accelaration

delta.age.df <- merged.age.tbl[, .(Horvath.Est, mAge_Hovath, realAge, Status, Sex)] 
delta.age.df[["delta_Age"]] <- delta.age.df$Horvath.Est - delta.age.df$realAge 

delta.age.melted.df <- delta.age.df[, .(delta_Age, Status, Sex)]
delta.age.df$Sex    <- "All"
delta.age.melted.df <- rbind(delta.age.melted.df, delta.age.df[, .(delta_Age, Status, Sex)])

ggplot(delta.age.melted.df, aes(y = delta_Age, x = Status, fill = Sex)) +
  geom_violin(trim = F, width = 0.8) +
  geom_boxplot(width = 0.1, color = "black", fill = "white") +
  scale_fill_viridis(discrete = TRUE, alpha = 0.5) +
  facet_wrap(~ Sex, ncol = 3) +
  labs(title = " ", y = "DNAm Age Hovath - Chronological Age", x = " ") + 
  #  theme_ipsum() +
  theme( panel.background = element_blank(),
         plot.title = element_text(size = 10),
         axis.title = element_text(size = 10),
         axis.text.x = element_text(angle = 0, hjust = 1), 
         legend.position = "none")
