#$1 bam file, $2 portein name,
#getting read and coverge
BASEDIR='/home/hsher/seqdata/eclip_raw/'
OUTDIR='/home/hsher/projects/'
bedtools coverage -b $BASEDIR$1 -a /home/hsher/projects/gencode_combine_sorted.gff3 > $OUTDIR$2.coverage
# count by transcript type
cat $OUTDIR$2.coverage | cut -d $'\t' -f 3,10 | sort | uniq -c > $OUTDIR$2.readcount
# filter for key feature
awk '{ if ($3=="transcript" && $10 > 500 && $13 > 0.005){print}}' $OUTDIR$2.coverage >$OUTDIR$2.key_transcript
# merge repeat transcript
awk '{ if ($3=="exon" && $10 > 500 && $13 > 0.005){print}}' $OUTDIR$2.coverage >$OUTDIR$2.key_exon