##### Visualize mapping efficiency

library(tidyverse)
library(magrittr)
library(data.table)
library(reshape2)
library(FactoMineR)


samps = fread("/netfiles/thermofly/METADATA/Thermofly_metadata.vNov11.2024.tsv")
##
in_dat <- fread("bam_extract.txt")
names(in_dat) = c("chr","L","mappedReads","mean_cov","sd_cov","sampleId")

aff = c(
"JAZHFB010000004.1",
"JAZHFB010000003.1",
"CM074799.1",
"CM074800.1",
"CM074801.1",
"CM074802.1",
"JAZHFB010000007.1",
"CM074803.1"
)

in_dat %<>%
mutate(sp_chrs =
case_when(chr %in%  aff ~ "Affinis",)
)

grep("JAWNKY", in_dat$chr) -> alg_chrs

in_dat$sp_chrs[alg_chrs] = "Algonquin"
in_dat$sp_chrs[is.na(in_dat$sp_chrs)] = "Athabasca"


in_dat %>%
filter(L > 500000) %>%
group_by(sp_chrs, sampleId) %>%
summarize(Mean = mean(mean_cov)) ->
summarized_covs

summarized_covs %>%
dcast(sampleId~sp_chrs, value.var = "Mean")

summarized_covs %>%
group_by(sampleId) %>%
slice_max(Mean) ->
heuristic_species_id


## cast wide for PCA
in_dat %>%
filter(L > 500000) %>%
dcast(sampleId~chr+sp_chrs, value.var = "mean_cov") ->
in_dat_cast

in_dat_cast[-1] ->
in_dat_cast.redux

rownames(in_dat_cast.redux) = in_dat_cast[,1]

in_dat_cast.redux %>%
PCA(graph = F) ->
PCA.calc.obj

PCA.calc.obj$ind$coord %>%
as.data.frame %>%
mutate(sampleId = rownames(.)) %>%
left_join(heuristic_species_id) %>%
mutate(sampleId = paste("D_aff.wild.",sampleId, sep = "")) %>%
left_join(samps)
->
PCA.coords

PCA.coords %>%
ggplot(aes(
x=Dim.1,
y=Dim.2,
color=altitude,
shape=sp_chrs,
size=Mean
)) +
geom_point() ->
f1stpass.PCA.plot

ggsave(f1stpass.PCA.plot, file = "f1stpass.PCA.plot.pdf")







