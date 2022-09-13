#!/bin/bash

#It appends the name of that genome to the header.
for i in *.fasta
do
sed -i -e "s/>/\>${i%.fasta}-/g" $i
done

#It makes BLAST database.
for i in *.fasta
do
makeblastdb -in $i -input_type fasta -dbtype nucl -title ${i%.fasta} -out ${i%.fasta}
done

#It generates an alias for all genome in the folder.
blastdb_aliastool -dblist "dana dere dgri dmoj dper dpse dsec dsim dvir dwil dyak" -dbtype nucl -out drosophilidae -title "drosophilidae"

