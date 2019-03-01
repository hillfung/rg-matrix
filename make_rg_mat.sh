#!/bin/bash

module load python/2.7.9

if [[ $( find ./*_rg.log -maxdepth 0 2> /dev/null | wc -l ) != 0 ]];then
	echo -e "WARNING: directory contains log files. These were removed"
	rm *_rg.log
fi
if [[ $( find names.list -maxdepth 0 2> /dev/null | wc -l ) != 0 ]];then
	echo -e 'WARNING: found "names.list". This was removed'
	rm names.list
fi
if [[ $( find sumstats.list -maxdepth 0 2> /dev/null | wc -l ) != 0 ]];then
	echo -e 'WARNING: found "sumstats.list". This was removed'
	rm sumstats.list
fi

ls *.sumstats.gz | sed 's/.sumstats.gz//g' > names.list

for i in $(seq 1 $(wc -l < names.list));do
    if [[ ${i} == $(wc -l < names.list) ]];then
        echo -en $(tail -n1 names.list)".sumstats.gz" > sumstats.list
    else
        for j in $(seq ${i} $(wc -l < names.list));do
            if [[ ${j} == $(wc -l < names.list) ]];then
                echo -en $(tail -n1 names.list)".sumstats.gz" >> sumstats.list
            else
                echo -en $(head -n${j} names.list | tail -n1)".sumstats.gz," >> sumstats.list
            fi
        done
    fi
    python ~/GWAS/resources/software/ldsc/ldsc.py \
    --rg $(head -n${i} names.list | tail -n1).sumstats.gz,$(cat sumstats.list) \
    --ref-ld-chr ~/GWAS/resources/LDScores/eur_w_ld_chr/ \
    --w-ld-chr ~/GWAS/resources/LDScores/eur_w_ld_chr/ \
    --out ${i}_rg
    rm sumstats.list
done

for i in $(seq 1 $(wc -l < names.list));do
    START=$(grep -En '^p1 +p2' ${i}_rg.log | cut -f1 -d':')
    tail -n+${START} ${i}_rg.log | grep '.sumstats.gz' > dum
    awk 'BEGIN{OFS="\t"};{print $2,$3 > "dum.rg"};{print $2,$4 > "dum.se"};{print $2,$11 > "dum.cti"}' < dum
    if [[ ${i} == 1 ]];then
        mv dum.rg rg.mat.1
        mv dum.se se.mat.1
        mv dum.cti cti.mat.1
    else
        join -a1 -t $'\t' -1 1 -2 1 rg.mat.1 -o auto dum.rg > foo && mv foo rg.mat.1
        join -a1 -t $'\t' -1 1 -2 1 se.mat.1 -o auto dum.se > foo && mv foo se.mat.1
        join -a1 -t $'\t' -1 1 -2 1 cti.mat.1 -o auto dum.cti > foo && mv foo cti.mat.1
        rm dum.*
    fi
    rm dum
    mv ${i}_rg.log $(head -n${i} names.list | tail -n1)_rg.log
done

for i in $( echo 'rg.mat se.mat cti.mat' );do
    sed -i 's![ \t]*$!!g' ${i}.1
    cut -f1 ${i}.1 > files
    cut -f2- ${i}.1 > matrix
    for j in $(seq 1 $(wc -l < ${i}.1));do
    cut -f${j} matrix | tail -n+$((j + 1)) | rs -c$'\t' -C$'\t' -T > dum
    head -n${j} matrix | tail -n1 | paste - dum -d $'\t' >> ${i}.2
    rm dum
    done
    paste files ${i}.2 -d $'\t' > ${i}
    rm ${i}.? files matrix
done
rm names.list