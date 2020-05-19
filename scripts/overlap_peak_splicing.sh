awk 'ARGIND == 1 {save[$5] = 1} ARGIND == 2 {split($9,fields, ";"); b = substr(fields[3],9); if(b in save) {print}}' $1 $2 > $3
