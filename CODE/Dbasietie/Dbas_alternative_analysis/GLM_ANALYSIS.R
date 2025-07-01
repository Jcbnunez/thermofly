#### Hawaiian Drosophila ablysis - pt 2

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

####
###
files <- system(paste("ls /gpfs2/scratch/jcnunez/thermofly/basisetae/analyses/GLM/*"), intern = T)

glm.out = 
  foreach(i=files,
          .combine = "rbind")%do%{
            
            tmp <- get(load(i))    
            
          }

glm.out %>%
  group_by(chr,pos) %>%
  arrange(chr,pos) %>%
  ungroup() %>%
  mutate(pos.id = 1:dim(glm.out)[1]) ->
  glm.out

glm.out %>%
  filter(p_ltr < 0.01) %>%
  filter(p_ltr < m.p_lrt.05) %>%
  dim

glm.out %>%
  filter(p_ltr < 0.01) %>%
  filter(p_ltr < m.p_lrt.05) %>%
  .$chr %>% unique() -> chrs.focal


glm.out %>%
  arrange(p_ltr) %>%
  mutate(p_rank = 1:dim(.)) %>%
  dplyr::select(chr, pos, p_rank) ->
  p_rank_df

glm.out %>%
  arrange(m.p_lrt.05) %>%
  mutate(p_random_r = 1:dim(.)) %>%
  dplyr::select(chr, pos, p_random_r) ->
  p_norm_df

dim(glm.out)[1] -> L

left_join(glm.out, p_rank_df) %>%
  left_join(p_norm_df) %>%
  mutate(p_rnpv = p_rank/L,
         p_rand_rnpv = p_random_r/L
         ) -> glm.out
glm.out %<>% mutate(SNP_id = paste(chr, pos, sep = "_"))


#####
samps <- fread("/netfiles/thermofly/METADATA/Thermofly_metadata.vNov11.2024.tsv")

#####
genofile <- seqOpen("/netfiles/thermofly/ANALYSES/Dbas_paper/D.basisetae.GATK.pipe.gds", allow.duplicate=T)
seqResetFilter(genofile)


samps.DBAS <- samps %>% filter(sampleId %in% seqGetData(genofile, "sample.id"))
#### Load the filtering SNP object -- which JCBN created
snps.dt <- data.table(chr=seqGetData(genofile, "chromosome"),
                      pos=seqGetData(genofile, "position"),
                      variant.id=seqGetData(genofile, "variant.id"),
                      nAlleles=seqNumAllele(genofile),
                      missing=seqMissing(genofile))
snps.dt %<>% mutate(SNP_id = paste(chr, pos, sep = "_"))
###

glm.out$SNP_id -> glm_snps
snps.dt %>% filter(SNP_id %in% glm_snps) %>% .$variant.id ->
  anchors

snpgdsGetGeno(genofile, 
              snp.id=anchors,
              verbose=TRUE, with.id=TRUE) ->
  GENO.matrix.all
GENO.matrix.all$genotype -> geno_matrix
colnames(geno_matrix) = paste(seqGetData(genofile, "chromosome"), seqGetData(genofile, "position"), sep = "_")  
rownames(geno_matrix) = seqGetData(genofile, "sample.id")  
geno_matrix = data.frame(geno_matrix)

data.frame(
  SNP_id = colnames(geno_matrix),
  NAS = colSums(is.na(geno_matrix))
) -> missing_data

left_join(glm.out, missing_data) ->
  glm.out.miss


## # FST
fst.snp = 
foreach(i=1:dim(glm.out)[1], .combine = "rbind",
        .errorhandling = "remove")%do%{
  
  message(paste(i, dim(glm.out)[1], sep = "|"))

  snps.dt %>%
    filter(chr == glm.out[i,]$chr) %>%
    filter(pos == glm.out[i,]$pos) %>%
    .$variant.id -> anchor
  
  v <- snpgdsFst(genofile, 
                 snp.id=anchor,
                 sample.id=samps.DBAS$sampleId, 
                 population=as.factor(samps.DBAS$locale),
                 autosome.only=FALSE,
                 method="W&C84") 
  
  data.frame(
    chr = glm.out[i,]$chr,
    pos = glm.out[i,]$pos,
    fst = v$Fst
  )
  
}

save(fst.snp, file = "fst.snp.Rdata")
load("fst.snp.Rdata")

glm.out.miss %>%
  filter(NAS <= 1) %>%
  full_join(fst.snp) ->
  glm.out.fst

glm.out.fst$fst[which(glm.out.fst$fst < 0)] = 0

glm.out.fst %<>%
  mutate(p_ltr.cor = p.adjust(p_ltr, "fdr"))

glm.out.fst %>%
  ggplot(aes(
    x=-log10(p_ltr),
    y=fst,
    color = p_ltr < m.p_lrt.05
  )) + 
  geom_point(size = 2.3) +
  geom_hline(yintercept = quantile(glm.out.fst$fst, 
                                   0.90, na.rm = T),
             linetype = "dashed") + 
  geom_vline(xintercept = -log10(0.01),
             linetype = "dashed") + 
  theme_bw() + ylim(0, 0.6) ->
  fst.lrt

ggsave(fst.lrt, file = "fst.lrt.pdf")
####
glm.out.fst %>%
  filter(p_ltr < m.p_lrt.05) %>%
  filter(p_ltr < 0.01) %>%
  filter(fst > quantile(glm.out.fst$fst, 
                        0.90, na.rm = T)) -> 
  FSTGLM_SNPS.outlier

write.table(FSTGLM_SNPS.outlier, 
            file = "FSTGLM_SNPS.outlier.Feb7.V2.txt", 
            append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,
            col.names = TRUE, 
            qmethod = c("escape", "double"),
            fileEncoding = "")
####
FSTGLM_SNPS.outlier %<>% mutate(SNP_id = paste(chr, pos, sep = "_"))

annot_file<-"/gpfs2/scratch/jcnunez/thermofly/basisetae/analyses/D.basisetae_annotated_SNPs.txt"

annot <- fread(annot_file)

names(annot)[3] = "SNP_id"
FSTGLM_SNPS.outlier %>%
  left_join(annot) ->
  FSTGLM_SNPS.outlier.annot

FSTGLM_SNPS.outlier.annot %>%
  group_by(SNP_id) %>%
  slice_head ->
  FSTGLM_SNPS.outlier.annot.uniq
  

write.table(FSTGLM_SNPS.outlier.annot.uniq, 
            file = "FSTGLM_SNPS.outlier.annot.uniq.TABLE.txt", 
            append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,
            col.names = TRUE, 
            qmethod = c("escape", "double"),
            fileEncoding = "")


FSTGLM_SNPS.outlier.annot.uniq$Annotation  %>% table
FSTGLM_SNPS.outlier.annot.uniq$chr  %>% table
FSTGLM_SNPS.outlier.annot.uniq$Gene_Name  %>% table


FSTGLM_SNPS.outlier.annot.uniq %>% filter(Gene_Name == "g193") %>% .$chr -> chr.i
FSTGLM_SNPS.outlier.annot.uniq %>% filter(Gene_Name == "g193") %>% .$pos -> pos.i

#####

snps.dt %>%
  filter(SNP_id %in% FSTGLM_SNPS.outlier$SNP_id) %>%
  .$variant.id ->
  target

snpgdsGetGeno(genofile, 
              snp.id=target, 
              verbose=TRUE, with.id=TRUE) ->
  GENO.matrix.targ

snps.dt %>%
  filter(SNP_id %in% FSTGLM_SNPS.outlier$SNP_id)

FSTGLM_SNPS.outlier.annot.uniq %>%
  filter()

snps.dt %>%
  filter(variant.id %in% GENO.matrix.targ$snp.id)


GENO.matrix.targ$genotype %>%
  as.data.frame() %>%
  mutate(samp = GENO.matrix.targ$sample.id) ->
  forplot

colnames(forplot) = c(GENO.matrix.targ$snp.id, "samp")

forplot %>%
  melt(id = "samp", variable.name = "variant.id") %>%
  mutate(variant.id = as.numeric(as.character(variant.id))) %>%
  left_join(snps.dt) %>%
  ggplot(aes(
    x=SNP_id,
    y=samp,
    fill = as.character(value)
  )) + coord_flip() + geom_point(shape = 22, size = 5) ->
  genotypes_tmp

ggsave(genotypes_tmp, file ="genotypes_tmp.pdf",
       w= 9, h =  5)


###
glm.out.fst %>%
  filter(chr == "JAWNLB010000114.1") %>%
  filter(!is.na(p_ltr < m.p_lrt.05)) %>%
  ggplot(aes(
    x=pos,
    y=-log10(p_ltr),
    color = p_ltr < m.p_lrt.05,
  )) +
  geom_point(size = 2.3) +
  theme_bw() +
  facet_grid(~chr, scale = "free_x")-> plot.114.manh

ggsave(plot.114.manh, file = "plot.114.manh.pdf",
       w = 4, h = 3)



glm.out.fst %>%
  filter(chr == "JAWNLB010000114.1") %>%
  filter(!is.na(p_ltr < m.p_lrt.05)) %>%
  ggplot(aes(
    x=pos,
    y=fst,
    color = p_ltr < m.p_lrt.05,
  )) +
  geom_point(size = 2.3) +
  theme_bw() +
  facet_grid(~chr, scale = "free_x")-> plot.114.mfst

ggsave(plot.114.mfst, file = "plot.114.mfst.pdf",
       w = 4, h = 3)


glm.out.fst %>%
  filter(chr == "JAWNLB010000114.1") %>%
  filter(!is.na(p_ltr < m.p_lrt.05)) %>%
  arrange(-fst) %>% as.data.frame() %>% head


