######################### make annotation ######################### 
chrs <- 1:22
for chr in $(seq 1 22)
do
    for ct in $( cat ${celltype_pt} )
    do
        python ./ldsc/make_annot.py \
               --bed-file ${bed_dir}/${ct} \
               --bimfile 1000G.EUR.QC.${chr}.bim \
               --annot-file 01.Anno.split/${ct}.chr${chr}.annot.gz
    done
done


######################### estimate LD score ######################### 
for ct in $( cat ${celltype_pt} )
do
        for chr in {1..22}
        do
                ct_chr_annot_pt="01.Anno.split/${ct}.chr${chr}.annot.gz"
                echo "Processing Celltype: ${ct}; Chrom: chr${chr}; Annot File: ${ct_chr_annot_pt}"
                output_pt="02.ld_score/${ct}.chr${chr}"
                python2 ldsc/ldsc.py \
                        --l2 \
                        --bfile 1000G.EUR.QC.${chr} \
                        --print-snps listHM3.txt \
                        --ld-wind-cm 1 \
                        --thin-annot \
                        --annot ${ct_chr_annot_pt}\
                        --out ${output_pt}
        done
done


######################### Partition heritability ######################### 
for ct in $( cat ${celltype_pt} )
do
        ct_qtl_pt="02.ld_score/${ct}.chr"
        for disease in "${diseaseV[@]}"
        do
                output_pt="output/${disease}.${ct}"
                python2  /home/wangwenhui/apps/ldsc/ldsc.py \
                        --h2 ${disease}.sumstats.gz \
                        --ref-ld-chr 1000G_Phase3_baselineLD_v2.2_ldscores/baselineLD.,${ct_qtl_pt} \
                        --w-ld-chr 1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
                        --frqfile-chr 1000G_Phase3_frq/1000G.EUR.QC. \
                        --overlap-annot \
                        --print-coefficients \
                        --out ${output_pt}
        done
done
