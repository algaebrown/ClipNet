# split into file
awk '{ if ($3 == "exon"){ print } }' $1 > $1_exon.gff3
awk '{split($9, fields, ";"); if (fields[5] == "gene_type=protein_coding"){print}}' $1_exon.gff3 > $1_exon_cds.gff3
awk '{split($9, fields, ";"); if (fields[5] == "gene_type=lncRNA"){print}}' $1_exon.gff3 > $1_exon_lncRNA.gff3
awk '{split($9, fields, ";"); if (fields[5] == "gene_type=miRNA"){print}}' $1_exon.gff3 > $1_exon_miRNA.gff3
awk '{split($9, fields, ";"); if (fields[5] == "gene_type=snoRNA"){print}}' $1_exon.gff3 > $1_exon_snoRNA.gff3

awk '{ if ($3 == "transcript"){ print } }' $1> $1_transcript.gff3
awk '{split($9, fields, ";"); if (fields[5] == "gene_type=protein_coding"){print}}' $1_transcript.gff3 > $1_intron_cds.gff3


cat $1 | grep three_prime_UTR > $1_3utr.gff3
cat $1 | grep five_prime_UTR > $1_5utr.gff3

# extract transcript id for exon
awk '{split ($9, fields, ":"); {print fields[2]}}' $1_exon_cds.gff3 | sort | uniq > $1_exon_cds.esnt
awk '{split ($9, fields, ":"); {print fields[2]}}' $1_exon_lncRNA.gff3 | sort | uniq > $1_exon_lncRNA.esnt

awk '{split($9, fields, ";"){print substr(fields[1],9)}}' $1_3utr.gff3 | sort | uniq > $1_3utr.esnt
awk '{split($9, fields, ";"){print substr(fields[1],9)}}' $1_5utr.gff3 | sort | uniq > $1_5utr.esnt

awk '{split ($9, fields, ";"); {print substr(fields[1],4)}}' $1_intron_cds.gff3 | sort | uniq > $1_intron_cds.esnt