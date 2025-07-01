
args = commandArgs(trailingOnly=TRUE)
i=as.numeric(args[1])

####

library(tidyverse)
library(magrittr)
library(data.table)
library(SeqArray)
library(gdsfmt)
library(SNPRelate)
library(adegenet)
library(reshape2)
library(FactoMineR)
require(gtools)
require(foreach)
library(lme4)

#####
samps <- fread("/netfiles/thermofly/METADATA/Thermofly_metadata.vNov11.2024.tsv")

#####
genofile <- seqOpen("/gpfs2/scratch/jcnunez/thermofly/basisetae/mapping/D.basisetae.GATK.pipe.gds", allow.duplicate=T)
seqResetFilter(genofile)

#### Load the filtering SNP object -- which JCBN created
snps.dt <- data.table(chr=seqGetData(genofile, "chromosome"),
                      pos=seqGetData(genofile, "position"),
                      variant.id=seqGetData(genofile, "variant.id"),
                      nAlleles=seqNumAllele(genofile),
                      missing=seqMissing(genofile))
####

## prepare windows
### define windows
win.bp <- 1.5e5
step.bp <- win.bp+1

chrs <- snps.dt$chr %>% unique
wins <- foreach(chr.i=chrs,
                .combine="rbind", 
                .errorhandling="remove")%dopar%{
                  
                  tmp <- snps.dt %>%
                    filter(chr == chr.i)
                  
                  data.table(chr=chr.i,
                             start=seq(from=min(tmp$pos), to=max(tmp$pos)-win.bp, by=step.bp),
                             end=seq(from=min(tmp$pos), to=max(tmp$pos)-win.bp, by=step.bp) + win.bp)
                }

wins[,i:=1:dim(wins)[1]]

dim(wins)
######

####
#4,532,545 snps
snps.dt %>%
  filter(nAlleles == 2, missing < 0.05) ->
  snps.dt.biallelic_lowmiss

seqSetFilter(genofile, 
             variant.id=snps.dt.biallelic_lowmiss$variant.id)

snps.dt.biallelic_lowmiss[,global_af:=seqGetData(genofile, "annotation/info/AF")$data]

snps.dt.biallelic_lowmiss %>%
  filter(global_af >= 0) ->
  snps.dt.biallelic_lowmiss.common

seqSetFilter(genofile, 
             variant.id=snps.dt.biallelic_lowmiss.common$variant.id)

dp <- seqGetData(genofile, "annotation/format/DP")
ad <- seqGetData(genofile, "annotation/format/AD")

col.min <- apply(dp, MARGIN = 2, FUN = min)
snps.dt.biallelic_lowmiss.common %>%
  mutate(min_cov = col.min) %>%
  filter(min_cov >= 4) ->
  snps.dt.biallelic_lowmiss.common.4X

#20,649 SNPs
seqResetFilter(genofile)
seqSetFilter(genofile, 
             variant.id=snps.dt.biallelic_lowmiss.common.4X$variant.id)


snpgdsGetGeno(genofile, 
              snp.id=snps.dt.biallelic_lowmiss.common.4X$variant.id, 
              verbose=TRUE, with.id=TRUE) ->
  GENO.matrix

MATRIX <- GENO.matrix$genotype
colnames(MATRIX)  <- GENO.matrix$snp.id
rownames(MATRIX)  <- GENO.matrix$sample.id

load("/gpfs2/scratch/jcnunez/thermofly/basisetae/analyses/Dbas.pca.object.Rdata")

####

scope = wins[i]
snps.dt %>%
  filter(chr == scope$chr) %>%
  filter(pos >= scope$start & pos >= scope$end ) ->
  selected.snps

######
sub.matrix = MATRIX[,which(colnames(MATRIX) %in% selected.snps$variant.id)]
REPS=100

results=        
  foreach(k=1:dim(sub.matrix)[1],
          .combine = "rbind",
        .errorhandling = "remove")%do%{
          
          data.frame(
            sampleId = row.names(sub.matrix),
            af=sub.matrix[,k],
            pc1=pca.object$ind$coord[,1],
            selected.snps[k]
          ) %>% left_join(samps) -> tmp
          
          chr.t=unique(tmp$chr)
          pos.t=unique(tmp$pos)
          variant.t=unique(tmp$variant.id)
          
          message(paste(k,chr.t,
                        pos.t,variant.t,
                        sep = "_"))
          
    foreach(j=0:REPS,
                  .combine = "rbind",
                  .errorhandling = "remove")%do%{
                    
          if(j==0){
            
            t0.real <- glm(af/2 ~ pc1, data = tmp, 
                           family= binomial)
            t1.real <- glm(af/2 ~ pc1 + as.factor(locale), data = tmp, 
                           family= binomial)
            
            p_lrt=anova(t1.real, t0.real, test="Chisq")[2,5]
            hab.AIC = extractAIC(t1.real)[2]
            null.AIC = extractAIC(t0.real)[2]
          } #if J ==0
          if(j>0){
                      
           t0.real <- glm(af/2 ~ pc1, data = tmp, 
                          family= binomial)
           t1.real <- glm(af/2 ~ pc1 + sample(as.factor(locale)), data = tmp, 
                          family= binomial)
           
           p_lrt=anova(t1.real, t0.real, test="Chisq")[2,5]
           hab.AIC = extractAIC(t1.real)[2]
           null.AIC = extractAIC(t0.real)[2]
          } # if j>0
                    obs <-
                      data.table(
                        chr=chr.t,
                        pos=pos.t,
                        perm=j,
                        p_lrt=p_lrt,
                        b_hab=summary(t1.real)$coef[2,1], 
                        se_hab=summary(t1.real)$coef[2,2],
                        n0=sum(tmp$af==0), 
                        n1=sum(tmp$af==1),
                        af=mean(tmp$af),
                        hab.AIC = hab.AIC,
                        null.AIC = null.AIC)
                    return(obs)
          }#j
          }#k

results %>%
  group_by(chr, pos,
           perm==0) %>%
  summarize(
  m.p_lrt.95 = quantile(p_lrt, 0.95, na.rm = T),
  m.p_lrt.05 = quantile(p_lrt, 0.05, na.rm = T),
  m.b_hab.95 = quantile(b_hab, 0.95, na.rm = T),
  m.b_hab.05 = quantile(b_hab, 0.05, na.rm = T),
  m.se_hab.95 = quantile(se_hab, 0.95, na.rm = T),
  m.se_hab.05 = quantile(se_hab, 0.05, na.rm = T),
  na = sum(is.na(p_lrt)),
  n0 = mean(n0, na.rm = T),
  n1 = mean(n1, na.rm = T),
  af = mean(af, na.rm = T)
  ) ->
  summaries

summaries %>%
  filter(`perm == 0` == TRUE) %>%
  dplyr::select(-`perm == 0`, 
         p_ltr=m.p_lrt.95,
         b_hab=m.b_hab.95,
         se_hab=m.se_hab.95,
         -m.p_lrt.05, -m.b_hab.05, -m.se_hab.05
         ) -> real

summaries %>%
  filter(`perm == 0` == FALSE) %>%
  .[,c("chr","pos","m.p_lrt.95","m.p_lrt.05",
       "m.b_hab.95","m.b_hab.05","m.se_hab.95",
       "m.se_hab.05")] -> permutations_df

left_join(real, permutations_df) ->
  summaries.final

file_name = paste(scope$chr,
                  scope$start,
                  scope$end,
                  sep = "_"
                  )
root = "/gpfs2/scratch/jcnunez/thermofly/basisetae/analyses/GLM/"

save(summaries.final,
     file = paste(root,file_name, "_dat.Rdata",sep = ""))

