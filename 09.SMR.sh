######################### run ######################### 
for disease in "${diseaseV[@]}"
do
    for cellType in "${cellTypeV[@]}"
    do
        echo "disease: ${disease}; CellType: ${cellType}"
        smr \
                --bfile "fromPlink/extractId" \
                --gwas-summary "${disease}.ma" \
                --beqtl-summary "${cellType}" \
                --out "smr_result/${disease}/${cellType}" \
                --thread-num 1 \
                --diff-freq-prop 0.99
    done
done
