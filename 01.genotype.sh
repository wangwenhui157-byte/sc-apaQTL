# fillter_chrom
plink --file genotype.sex --exclude snp2del.txt --recode --out genotype.sex.qc1


# filter_sex
plink --file genotype.sex.qc1 --remove iid2del.txt --recode --out genotype.sex.qc2


# snp_missingness
mkdir snp_missingness
plink --bfile genotype.sex.qc2 --make-bed --geno 0.03 --out snp_missingness/snp_missingness --noweb


# indiv_missingness
mkdir indiv_missingness
plink --bfile snp_missingness/snp_missingness --make-bed --mind 0.03 --out indiv_missingness/indiv_missingness --noweb


# check_sex
mkdir check_sex
plink --bfile indiv_missingness/indiv_missingness --check-sex --out check_sex/check_sex --noweb


# maf
mkdir maf
plink --bfile indiv_missingness/indiv_missingness --freq --out maf/freq --noweb
plink --bfile indiv_missingness/indiv_missingness --maf 0.01 --allow-extra-chr --make-bed  --out maf/maf
printf "s/X/23/\ns/Y/24/\ns/MT/26/\n" > maf/chr_coding
paste <(sed -f maf/chr_coding <(cut -f1 maf/maf.bim)) <(cut -f2- maf/maf.bim) > maf/maf_tmp.bim


# Hardy-Weinberg
mkdir hwe
plink --bfile maf/maf --hardy --hwe 0.001 --make-bed --out hwe2/hwe2 --noweb


# het
mkdir het
plink --bfile hwe2/hwe2 --het --out het/het --noweb
python filter_het.py --input het/het.het --output het/het.inds


# het filter
mkdir het_filter
plink --bfile hwe2/hwe2 --remove het/het.inds --make-bed --out het_filter/het_filter --noweb


# pca
mkdir pca 
plink --bfile het_filter/het_filter --recode --out pca/pca --noweb
printf "genotypename:    pca.ped\nsnpname:         pca.map\nindivname:       pca.ped\nevecoutname:     pca.evec\nevaloutname: pca.eval\naltnormstyle:    NO\nnumoutevec:      5\nfamilynames:     NO\nnumoutlieriter:  0\n" > pca/par.smartpca
cd pca
awk '{if($6==-9){$6=1; print}}' pca.ped > pca_new.ped
mv pca.ped pca_old.ped
mv pca_new.ped pca.ped
smartpca -p par.smartpca | tee pca.pca


# ibd
mkdir ibd
plink --bfile het_filter/het_filter --genome --out ibd/ibd --noweb


# grm
mkdir grm
awk '$1 == "2324"' het_filter/het_filter.bim | cut -f2 > grm/y_par.snps
plink --bfile het_filter/het_filter --exclude grm/y_par.snps --make-bed --out grm/grm --noweb
gcta64 --bfile grm/grm --make-grm --out grm/grm


# grm_singleton
mkdir grm_singleton
gcta64 --grm grm/grm --grm-singleton 0.05 --out grm_singleton/grm_singleton


# grm_filter
mkdir grm_filter
gcta64 --grm grm/grm --grm-cutoff 0.125 --make-grm --out grm_filter/grm_filter
gcta64 --grm grm_filter/grm_filter --pca 6 --out grm_filter/grm_filter


# grm_subset
mkdir grm_subset
plink --bfile het_filter/het_filter --keep grm_filter/grm_filter.grm.id --make-bed --out  grm_subset/grm_subset


# forward_strand
mkdir forward_strand
snpflip -b grm_subset/grm_subset.bim -f ../refGenome/Homo_sapiens.GRCh37.dna.primary_assembly.fa -o forward_strand/forward_strand_snpflip
plink --bfile grm_subset/grm_subset --exclude forward_strand/forward_strand_snpflip.ambiguous --make-bed --out forward_strand/forward_strand_noambiguous --noweb
plink --bfile forward_strand/forward_strand_noambiguous --flip forward_strand/forward_strand_snpflip.reverse --make-bed --out forward_strand/forward_strand --noweb
mkdir forward_strand/dup
mv forward_strand/forward_strand* forward_strand/dup/
plink --bfile forward_strand/dup/forward_strand --recode vcf --out forward_strand/dup/dedup
python oneK1Kqc3.py --input forward_strand/dup/dedup.vcf --GenomeVersion hg19
plink --bfile forward_strand/dup/forward_strand --exclude forward_strand/dup/dedup.dupId --make-bed --out forward_strand/forward_strand


# id_map
mkdir idMap
bash idMap.sh > idMap.log 2>&1 &
python -u /data1/weiyihu/apa/sc/scripts/idMap.py \
    --path . \
    --bim forward_strand/forward_strand.bim \
    --out idMap/forward_strand.idMap
plink --bfile forward_strand/forward_strand --update-name idMap/forward_strand.idMap --make-bed --out idMap/forward_strand


# makeGrmLoading
mkdir makeGrmLoading
cd makeGrmLoading
mkdir 1000GIntersect
gcta64 \
    --bfile 1000Genome/bfile/intersect/merge \
    --maf 0.01 \
    --autosome \
    --make-grm \
    --out 1000GIntersect/REF
gcta64 \
    --grm 1000GIntersect/REF \
    --pca 20 \
    --out 1000GIntersect/REF_pca20
gcta64 \
    --bfile 1000Genome/bfile/intersect/merge \
    --pc-loading 1000GIntersect/REF_pca20 \
    --out 1000GIntersect/REF_snp_loading

# Michigan_Imputation_Server
mkdir upload
for chr in {1..22}; do
    plink --bfile ../GSE196829/forward_strand/forward_strand --chr ${chr} --recode vcf --out upload/chr${chr}.vcf
    bcftools convert split/chr${chr}.vcf -Oz -o upload/chr${chr}.vcf.gz
done

python michiganImputationServer.py \
    --path upload \
    --mode imputation \
    --refpanel 1000g-phase-3-v5 \
    --phasing eagle \
    --population eur \
    --build hg19 \
    --r2filter 0 > michiganImputationServer.log 2>&1 &

# qc
mkdir qc
for chr in {1..22}; do
    echo ${chr}
    file=download/chr${chr}.info.gz
    zcat ${file} | awk -F "\t" '{if($7>=0.8){print $1}}' > qc/chr${chr}.id
done

for chr in {1..22}; do
    plink --vcf download/chr${chr}.dose.vcf.gz \
        --threads 2 \
        --extract qc/chr${chr}.id \
        --const-fid \
        --keep-allele-order \
        --make-bed \
        --out qc/chr${chr}
done


# merge
for i in {1..22}; do
    echo chr${i}.bed chr${i}.bim chr${i}.fam >> qc/merge.list
done
cd qc
plink2 --pmerge-list merge.list --make-bed --out merge
cd ..


# maf filter
cd qc
mkdir nogeno
plink2 --bfile merge --maf 0.05 --make-bed --out nogeno/qc2
