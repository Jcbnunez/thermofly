#### Mine Outliers

library(tidyverse)
library(data.table)
library(magrittr)

dat <- fread("/gpfs2/scratch/jcnunez/thermofly/basisetae/analyses/D.basisetae_annotated_SNPs.txt")

dat %>%
  filter(nAlleles == 2) %>%
  filter(Annotation_Impact %in%
  c("HIGH","MODERATE")) -> hits

write.table(hits,
            file = "top.hits.dbas.txt", 
            append = FALSE, quote = TRUE, 
            sep = "\t",
            eol = "\n", na = "NA", 
            dec = ".", row.names = FALSE,
            col.names = TRUE, 
            qmethod = c("escape", "double"),
            fileEncoding = "")

wr