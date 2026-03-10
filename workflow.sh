cd ~/projects/ngs-course-final-assignment-25-26-valencil

</data-shared/vcf_examples/luscinia_vars_flags.vcf.gz zcat |
  grep -v '^#' \
> no-headers.vcf


IN=no-headers.vcf
<$IN grep -E -o 'DP=([^;]+)' | sed 's/DP=//' > col-dp.tsv
<$IN awk '{if($0 ~ /INDEL/) print "INDEL"; else print "SNP"}' > col-type.tsv


wc -l *.tsv
paste col-dp.tsv col-type.tsv > cols-all.tsv
