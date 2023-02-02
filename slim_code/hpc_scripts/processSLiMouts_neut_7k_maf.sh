dir=$1
mig=$2

module load bcftools

cd ${dir}/K7000_m${mig}_neutral/

it=1
while [ $it -le 20 ]
do
	ls ${dir}/K7000_m${mig}_neutral/i${it}_*.vcf.gz > ${dir}/K7000_m${mig}_neutral/i${it}_mergelist.txt
	bcftools merge -l ${dir}/K7000_m${mig}_neutral/i${it}_mergelist.txt -Oz -o ${dir}/K7000_m${mig}_neutral/i${it}_merge.vcf.gz
	bcftools index ${dir}/K7000_m${mig}_neutral/i${it}_merge.vcf.gz
	zcat ${dir}/K7000_m${mig}_neutral/i${it}_merge.vcf.gz | sed 's#./.#0|0#g' | bgzip -c > ${dir}/K7000_m${mig}_neutral/i${it}_merge_mtoref.vcf.gz
	bcftools index ${dir}/K7000_m${mig}_neutral/i${it}_merge_mtoref.vcf.gz
	vcftools --gzvcf ${dir}/K7000_m${mig}_neutral/i${it}_merge_mtoref.vcf.gz --maf 0.05 --max-maf 0.95 --min-alleles 2 --max-alleles 2 --recode --out ${dir}/K7000_m${mig}_neutral/i${it}_merge_mtoref_maf
	vcftools --vcf ${dir}/K7000_m${mig}_neutral/i${it}_merge_mtoref_maf.recode.vcf --weir-fst-pop ${dir}/popfiles/Can13.txt --weir-fst-pop ${dir}/popfiles/Can40.txt --out ${dir}/K7000_m${mig}_neutral/i${it}_Can_tempfst_maf
	vcftools --vcf ${dir}/K7000_m${mig}_neutral/i${it}_merge_mtoref_maf.recode.vcf --weir-fst-pop ${dir}/popfiles/Lof07.txt --weir-fst-pop ${dir}/popfiles/LofContemp.txt --out ${dir}/K7000_m${mig}_neutral/i${it}_Lof_tempfst_maf
	vcftools --vcf ${dir}/K7000_m${mig}_neutral/i${it}_merge_mtoref_maf.recode.vcf --weir-fst-pop ${dir}/popfiles/Can13.txt --weir-fst-pop ${dir}/popfiles/LofContemp.txt --out ${dir}/K7000_m${mig}_neutral/i${it}_CanLof_contempfst_maf
	vcftools --vcf ${dir}/K7000_m${mig}_neutral/i${it}_merge_mtoref_maf.recode.vcf --weir-fst-pop ${dir}/popfiles/Can40.txt --weir-fst-pop ${dir}/popfiles/Lof07.txt --out ${dir}/K7000_m${mig}_neutral/i${it}_CanLof_historicfst_maf
	vcftools --vcf ${dir}/K7000_m${mig}_neutral/i${it}_merge_mtoref_maf.recode.vcf --weir-fst-pop ${dir}/popfiles/Lof11.txt --weir-fst-pop ${dir}/popfiles/Lof14.txt --out ${dir}/K7000_m${mig}_neutral/i${it}_Lof_sampfst_maf
	echo "m1e-3" > ${dir}/K7000_m${mig}_neutral/i${it}_maf_fsts
	echo "qtl1e-1" >> ${dir}/K7000_m${mig}_neutral/i${it}_maf_fsts
	echo "i1" >> ${dir}/K7000_m${mig}_neutral/i${it}_maf_fsts
	grep "weighted" ${dir}/K7000_m${mig}_neutral/i${it}_Can_tempfst_maf.log | cut -d' ' -f 7 >> ${dir}/K7000_m${mig}_neutral/i${it}_maf_fsts
	grep "weighted" ${dir}/K7000_m${mig}_neutral/i${it}_Lof_tempfst_maf.log | cut -d' ' -f 7 >> ${dir}/K7000_m${mig}_neutral/i${it}_maf_fsts
	grep "weighted" ${dir}/K7000_m${mig}_neutral/i${it}_CanLof_contempfst_maf.log | cut -d' ' -f 7 >> ${dir}/K7000_m${mig}_neutral/i${it}_maf_fsts
	grep "weighted" ${dir}/K7000_m${mig}_neutral/i${it}_CanLof_historicfst_maf.log | cut -d' ' -f 7 >> ${dir}/K7000_m${mig}_neutral/i${it}_maf_fsts
	grep "weighted" ${dir}/K7000_m${mig}_neutral/i${it}_Lof_sampfst_maf.log | cut -d' ' -f 7 >> ${dir}/K7000_m${mig}_neutral/i${it}_maf_fsts
	${dir}/plink2 --allow-extra-chr --vcf ${dir}/K7000_m${mig}_neutral/i${it}_merge_mtoref_maf.recode.vcf --freq --pheno ${dir}/popfiles/allpops.txt --loop-cats population --out ${dir}/K7000_m${mig}_neutral/i${it}_maf
	it=$(( $it+1 ))
done

paste ${dir}/popfiles/fstcats ${dir}/K7000_m${mig}_neutral/i*_maf_fsts > ${dir}/fstouts/K7000_m${mig}_neutral_maf_fsts
