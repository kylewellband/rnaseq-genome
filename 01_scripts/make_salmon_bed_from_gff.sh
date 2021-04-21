awk '($3 == gene){print $1 "\t" $4 "\t" $5 "\t" $9 "\t.\t" $7}' 02_reference/genes.gff #| perl -pe 's/(ID=)(.*)(;Name.*\w)(.*)/\2\4/' > 02_reference/salmon_genes.bed
