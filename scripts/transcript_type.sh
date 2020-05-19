echo transcript > $1_transcript_type
cut -d $'\t' -f 9 $1_transcript.gff3 | cut -d ";" -f 5 | sort | uniq -c >> $1_transcript_type
echo exon >> $1_transcript_type
cut -d $'\t' -f 9 $1_exon.gff3 | cut -d ";" -f 5 | sort | uniq -c >> $1_transcript_type