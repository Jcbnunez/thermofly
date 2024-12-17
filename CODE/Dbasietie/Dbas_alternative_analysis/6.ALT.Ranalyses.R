#### Hawaiian Drosophila ablysis

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

#library(BEDASSLE)

#####
samps <- fread("/netfiles/thermofly/METADATA/Thermofly_metadata.tsv")

#####
genofile <- seqOpen("/gpfs2/scratch/jcnunez/thermofly/basisetae/mapping/D.basisetae.GATK.pipe.gds", allow.duplicate=T)
seqResetFilter(genofile)

#### Load the filtering SNP object -- which JCBN created
snps.dt <- data.table(chr=seqGetData(genofile, "chromosome"),
                      pos=seqGetData(genofile, "position"),
                      variant.id=seqGetData(genofile, "variant.id"),
                      nAlleles=seqNumAllele(genofile),
                      missing=seqMissing(genofile))
#4,532,545 snps
snps.dt %>%
  filter(nAlleles == 2, missing < 0.05) ->
  snps.dt.biallelic_lowmiss

seqSetFilter(genofile, 
             variant.id=snps.dt.biallelic_lowmiss$variant.id)

snps.dt.biallelic_lowmiss[,global_af:=seqGetData(genofile, "annotation/info/AF")$data]

snps.dt.biallelic_lowmiss %>%
  filter(global_af > 0.05) ->
  snps.dt.biallelic_lowmiss.common

seqSetFilter(genofile, 
             variant.id=snps.dt.biallelic_lowmiss.common$variant.id)

dp <- seqGetData(genofile, "annotation/format/DP")
ad <- seqGetData(genofile, "annotation/format/AD")

col.min <- apply(dp, MARGIN = 2, FUN = min)
snps.dt.biallelic_lowmiss.common %>%
  mutate(min_cov = col.min) %>%
  filter(min_cov >= 6) ->
  snps.dt.biallelic_lowmiss.common.6X

#20,649 SNPs
seqResetFilter(genofile)
seqSetFilter(genofile, 
             variant.id=snps.dt.biallelic_lowmiss.common.6X$variant.id)


snpgdsGetGeno(genofile, 
              snp.id=snps.dt.biallelic_lowmiss.common.6X$variant.id, 
              verbose=TRUE, with.id=TRUE) ->
  GENO.matrix

MATRIX <- GENO.matrix$genotype
colnames(MATRIX)  <- GENO.matrix$snp.id
rownames(MATRIX)  <- GENO.matrix$sample.id

samps.DBAS <- samps %>% filter(sampleId %in% GENO.matrix$sample.id)

v <- snpgdsFst(genofile, 
               snp.id=snps.dt.biallelic_lowmiss.common.6X$variant.id,
               sample.id=samps.DBAS$sampleId, 
               population=as.factor(samps.DBAS$locale),
               autosome.only=FALSE,
               method="W&C84")

v$MeanFst
#0.001883866

### RELATEDNESS
ibd.TOMS <- snpgdsIBDMoM(genofile, 
                    sample.id=filter(samps.DBAS, locale == "US-HI-Tom")$sampleId, 
                    snp.id=snps.dt.biallelic_lowmiss.common.6X$variant.id,
                    maf=0.05, missing.rate=0.05, num.thread=2,  autosome.only=FALSE)
ibd.coeff.TOMS <- snpgdsIBDSelection(ibd.TOMS)

ibd.Ola <- snpgdsIBDMoM(genofile, 
                         sample.id=filter(samps.DBAS, locale == "US-HI-Ola")$sampleId, 
                         snp.id=snps.dt.biallelic_lowmiss.common.6X$variant.id,
                         maf=0.05, missing.rate=0.05, num.thread=2,  autosome.only=FALSE)
ibd.coeff.OLA <- snpgdsIBDSelection(ibd.Ola)

rbind(mutate(ibd.coeff.TOMS, pop = "Toms" ), 
      mutate(ibd.coeff.OLA, pop = "Olaa")) ->
  kin.analysis

kin.analysis %>%
  group_by(pop) %>%
  summarise(m.kinship = mean(kinship))

#pop     m.kim
#<chr>   <dbl>
#  1 Olaa  0.00342
#  2 Toms  0.00783

#####
MATRIX %>% 
  PCA(graph = F, ncp = 5) ->
  pca.object
save(pca.object, file = "Dbas.pca.object.Rdata")
####
load("Dbas.pca.object.Rdata")

pca.object$ind$coord %>%
  as.data.frame() %>% 
  mutate(sampleId = rownames(.)) %>%
  left_join(samps) ->
  pca.meta.dim

#### PLOT PCA
pca.meta.dim %>%
  ggplot(aes(
    x=Dim.1,
    y=Dim.2,
    color = locale
  )) +
  geom_point() + 
  scale_color_manual(values = c("darkblue","firebrick2") ) +
  theme_bw() ->
  PCA12.dim

ggsave(PCA12.dim, file = "PCA12.dim.pdf", w = 5, h =4)

#### DO DISTANCE ANALYSIS
CalculateEuclideanDistance <- function(vect1, vect2) sqrt(sum((vect1 - vect2)^2))

pca.meta.dim %>%
  filter(locale == "US-HI-Ola") %>%
  .$sampleId -> OLA_SAMPS
permutations(n = length(OLA_SAMPS), r = 2, v = OLA_SAMPS) %>%
  as.data.frame()->
  OLA_PERMS

pca.meta.dim %>%
  filter(locale == "US-HI-Tom") %>%
  .$sampleId -> TOM_SAMPS
permutations(n = length(TOM_SAMPS), r = 2, v = TOM_SAMPS) %>%
  as.data.frame()->
  TOM_PERMS

OLA.D =
foreach(i=1:dim(OLA_PERMS)[1], 
        .combine = "rbind")%do%{

          i1 = pca.meta.dim %>% filter(sampleId == OLA_PERMS$V1[i])
          i2 = pca.meta.dim %>% filter(sampleId == OLA_PERMS$V2[i])
          
          x1=i1$Dim.1
          x2=i2$Dim.1
          
          y1=i1$Dim.2
          y2=i2$Dim.2
          
          euD = CalculateEuclideanDistance(c(x1,x2),c(y1,y2)) 
          
          data.frame(comp = "ola", euD)
          
}

TOM.D =
  foreach(i=1:dim(TOM_PERMS)[1], 
          .combine = "rbind")%do%{
            
            i1 = pca.meta.dim %>% filter(sampleId == TOM_PERMS$V1[i])
            i2 = pca.meta.dim %>% filter(sampleId == TOM_PERMS$V2[i])
            
            x1=i1$Dim.1
            x2=i2$Dim.1
            
            y1=i1$Dim.2
            y2=i2$Dim.2
            
            euD = CalculateEuclideanDistance(c(x1,x2),c(y1,y2)) 
            
            data.frame(comp = "Tom", euD)
            
          }

mean(OLA.D$euD)
sd(OLA.D$euD)

mean(TOM.D$euD)
sd(TOM.D$euD)

#### DAPC ANALYSIS
N_r <- length(GENO.matrix$sample.id)
dapc_r <- dapc(MATRIX, grp = pca.meta.dim$locale , n.da=100, n.pca=N_r/3)
optmin <- optim.a.score(dapc_r)
dapc_opt <- dapc(MATRIX, grp = pca.meta.dim$locale , n.da=100, n.pca=optmin$best )

pdf("comp.plot.pdf",width = 9, height = 4)
compoplot(dapc_opt)
dev.off()

dapc_opt$posterior %>%
  as.data.frame() %>% 
  mutate(sampleId = rownames(.)) %>%
  left_join(samps) ->
  dapc.meta.POST

dapc.meta.POST %>% 
  filter(locale == "US-HI-Ola") %>%
  filter(`US-HI-Ola` < 0.99) 

dapc.meta.POST %>% 
  filter(locale == "US-HI-Tom") %>%
  filter(`US-HI-Ola` < 0.99) 


## scatter
pdf("scatter.pdf")
scatter(dapc_opt)
dev.off()

dapc_opt$ind.coord %>%
  as.data.frame() %>% 
  mutate(sampleId = rownames(.)) %>%
  left_join(samps) ->
  dapc.meta.dim

dapc.meta.dim %>%
  ggplot(aes(
    x=LD1,
    fill = locale
  )) +
  geom_density() +
  theme_bw() ->
  dapc.dim

ggsave(dapc.dim, file = "dapc.dim.pdf", w = 5, h =4)
#### Compare PCA / DAPC

left_join(
pca.meta.dim[,c("sampleId", "Dim.2")],
dapc.meta.dim[,c("sampleId", "LD1")] ) %>%
  cor.test(~Dim.2+LD1, data = .)

left_join(
  pca.meta.dim[,c("sampleId", "Dim.1")],
  dapc.meta.dim[,c("sampleId", "LD1")] ) %>%
  cor.test(~Dim.1+LD1, data = .)


dapc_opt$var.contr %>%
  as.data.frame() %>%
  mutate(variant.id =  as.integer(rownames(.))) %>%
  left_join(snps.dt.biallelic_lowmiss.common.6X) ->
  var.contrib.snp6x
  
var.contrib.snp6x %>%
  arrange(-LD1) %>% head(100)


