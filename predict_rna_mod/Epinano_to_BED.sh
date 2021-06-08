input=$1
awk -F',' '{if (!/^#/) print $1"\t"$2"\t"$2"\t"$4"\t"$0}' $input | perl -ne 'chomp; @a = split /\s+/,$_; print join("\t",@a),"\n";' > $input.bed