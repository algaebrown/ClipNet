bedtools window -a $1 -b $2 -w 30 -c | sort | uniq -c