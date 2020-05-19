exon=$(awk '{ if ($3 == "exon"){ print } }' $1 | wc -l)
intron=$(awk '{ if ($3 == "transcript"){ print } }' $1 | wc -l)
fiveutr=$(awk '{ if ($3 == "five_prime_UTR"){ print } }' $1 | wc -l)
threeutr=$(awk '{ if ($3 == "three_prime_UTR"){ print } }' $1 | wc -l)

echo $exon, $intron, $fiveutr, $threeutr