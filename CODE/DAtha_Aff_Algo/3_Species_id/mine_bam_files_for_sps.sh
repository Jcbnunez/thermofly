### Loop for mining coverage data from BAMQC files
### Started Dec 6, 2024
### Jcbn

files=$( ls /netfiles/thermofly/ANALYSES/Affinis_Altitude/bamQCs/bams_qualimap/Qualimap_D_aff.wild.*/genome_results.txt)
for i in $files
do
echo $i
cat $i | sed -n '129,3862p;3862q' | sed -e 's/^\t//g' > $i.out.txt

awk 'BEGIN{OFS="\t"} {print $0, (FNR>1 ? FILENAME : FILENAME)}' $i.out.txt >> bam_extract.txt

done




