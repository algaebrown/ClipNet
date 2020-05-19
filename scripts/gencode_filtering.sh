### This whole script is to filter for canonical transcripts. output is gencode_combine_correct

# filter for canonical transcript (knownCanonical) in the whole gencode
awk 'ARGIND == 1 {save[$5] = 1} ARGIND == 2 {split($9,fields, ":"); if(fields[2] in save) {print}}' knownCanonical.txt gencode.v33.annotation.gff3 > gencode_filtered.gff3

# find exon
awk '{ if ($3 == "exon"){ print } }' gencode_filtered.gff3 > gencode_exon.gff3

# remove intermediate
rm gencode_filtered.gff3

# but transcript has id in different place
awk 'ARGIND == 1 {save[$5] = 1} ARGIND == 2 {split($9,fields, ";"); a = fields[1];b = gensub(/ID=/, "", 1, a);if(b in save) {print}}' knownCanonical.txt gencode.v33.annotation.gff3 > gencode_filter_for_transcript.gff3

awk '{ if ($3 == "transcript"){ print } }' gencode_filter_for_transcript.gff3> gencode_transcript.gff3

# remove intermediate
rm gencode_filter_for_transcript.gff3

# transcript - exon = intron (but need to make intron id for yourself
bedtools subtract -a gencode_transcript.gff3 -b gencode_exon.gff3 > gencode_intron.gff3

# find three prime and five prime utr
cat gencode.v33.annotation.gff3 | grep prime_UTR > gencode_utr.gff3


# and filter for canonical transcript only
awk 'ARGIND == 1 {save[$5] = 1} ARGIND == 2 {split($9,fields, ";"); a = fields[1];b = substr(a,9) ;if(b in save) {print}}' knownCanonical.txt gencode_utr.gff3 > gencode_utr_filtered.gff3

# remove intermediate
rm gencode_utr.gff3

# combine all features into gencode_combine_correct
cat gencode_exon.gff3 > gencode_combine.gff3
cat gencode_intron.gff3 >> gencode_combine.gff3
cat gencode_utr_filtered.gff3 >> gencode_combine.gff3

# sort file
sort -k1,1 -k4,4n ~/projects/gencode_combine.gff3 > gencode_combine_sorted.gff3

# remove intermediate
rm gencode_combine.gff3

# but these exons/intron include lncRNA and other non-coding RNAs
# use intron_exon_coordinate.sh to filter for protein-coding ones;
# use intron_annotation to get intron ids

# filter for protein_coding only
awk '{split($9, fields, ";"); if (fields[5] == "gene_type=protein_coding"){print}}' gencode_intron.gff3 > gencode_cds.gff3

awk '{split($9, fields, ";"); if (fields[5] == "gene_type=protein_coding"){print}}' gencode_exon.gff3 >> gencode_cds.gff3

# sort 
cat gencode_cds.gff3 | sort -k1,1 -k4,4n > gencode_cds_sorted.gff3

rm gencode_cds.gff3
