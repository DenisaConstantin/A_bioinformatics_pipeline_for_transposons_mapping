#!/bin/bash

echo "Folder"
read folder

singleline="singleline_"
echo "Genomic flanking length"
read LengthFL1

echo "Transposon length"
read LengthFL2

file_base_name="flankings_"
extension=".fa"

cd $folder

#It accesses each folder one at a time.
ls -d */ |
while read dirname 
do
 cd $dirname

#IT BUILDS COORD.TXT FILE.
#The variable "coordinate" takes the value of the first coordinate in the header of the respective remnant.
sed -i "1s/\../ ../" *.fa
coordinate=$(awk -F ":| " '{print$4}' *.fa)
if [[ $coordinate == "complement("* ]]
then
sed -i "s/(/\( /" *.fa
coordinate=$(awk -F " " '{print$4}' *.fa)
fi

#The variable "loc" is the chromosome in which the remnant is inserted.
loc=$(awk -F "=|:|;" '{print$4}' *.fa)
#I works on the BLAST results file, it searches for the highest value for column 12 (bitscore), if column 3 (pident) = 100 and column 12 has the maximum value > send the data to coord1.txt
awk 'BEGIN {max = 0} {if ($12>max) max=$12} END {echo max} {if ($3 == 100 && $12 == max) print $2,$9,$10,$13}' *.txt > coord1.txt
cat coord1.txt | while read chr start end strand
do
if [ $chr == $loc ] && [ $coordinate == $start ] || [ $coordinate == $end ]
then
echo $start $end $strand > $start.coord.txt
fi
done
cat *.coord.txt > coord.txt
rm -r *.coord.txt coord1.txt
#THE FILE COORD.TXT WAS BUILT.

#JUNCTION SEQUENCE EXTRACTION.
cat coord.txt | while read coordinate1 coordinate2 strand
do

if [ $strand = plus ]
then
#It calculates the values for START and STOP.
START=$( expr $coordinate2 - $LengthFL2 + 1 )
if [ $( expr $coordinate2 - $LengthFL2 + 1 ) -lt $coordinate1 ] 
then
START=$coordinate1
fi
STOP=$( expr $coordinate2 + $LengthFL1 )

#The command that performs the extraction.
seqkit subseq -r $START:$STOP ../../Genome/dmel.fasta > $file_base_name$coordinate1.$strand$extension
#It transforms the file into a single line file so you can count the nucleotides.
awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' $file_base_name$coordinate1.$strand$extension > $singleline$file_base_name$coordinate1.$strand$extension
Length=$( echo $singleline$file_base_name$coordinate1.$strand$extension | sed -e '1d' $singleline$file_base_name$coordinate1.$strand$extension | wc -c )
STOPh=$( expr $START + $Length - 1 )

rm $singleline$file_base_name$coordinate1.$strand$extension

#It edits the header of the file containing the junction sequence by adding the coordinate of interest, the strand, the extraction interval, the total length.
sed -i -e "s/>/\>Coordinate: $coordinate2. Strand: $strand. Range: $START:$STOP Length: /g" $file_base_name$coordinate1.$strand$extension
fi

#If the reference strand of the transposon is identified in the non-reference strand of the genome.
if [ $strand = minus ]
then
  if [  $coordinate1 -gt $coordinate2 ]
  then
  START=1 
  if ! [ $coordinate2 -le $LengthFL1 ]
  then

  START=$( expr $coordinate2 - $LengthFL1 )
  fi

  STOP=$( expr $coordinate2 + $LengthFL2 - 1 )
  if [ $( expr $coordinate2 + $LengthFL2 - 1 ) -gt $coordinate1 ]
  then
  STOP=$coordinate1
  fi

seqkit subseq -r $START:$STOP ../../Genome/dmel.fasta > $file_base_name$coordinate1.$strand$extension
awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' $file_base_name$coordinate1.$strand$extension > $singleline$file_base_name$coordinate1.$strand$extension
Length=$( echo $singleline$file_base_name$coordinate1.$strand$extension | sed -e '1d' $singleline$file_base_name$coordinate1.$strand$extension | wc -c )
STOPh=$( expr $START + $Length - 1 )

rm $singleline$file_base_name$coordinate1.$strand$extension

sed -i -e "s/>/\>Coordinate: $coordinate2. Strand: $strand. Range: $START:$STOP Length: /g" $file_base_name$coordinate1.$strand$extension

  else
  START=1
  if ! [ $coordinate1 -le $LengthFL1 ]
  then
  START=$( expr $coordinate1 - $LengthFL1 )
  fi

  STOP=$( expr $coordinate1 + $LengthFL2 - 1 )
  if [ $( expr $coordinate1 + $LengthFL2 - 1 ) -gt $coordinate2 ]
  then
  STOP=$coordinate2
  fi

  seqkit subseq -r $START:$STOP ../../Genome/dmel.fasta > $file_base_name$coordinate2.$strand$extension
  awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' $file_base_name$coordinate2.$strand$extension > $singleline$file_base_name$coordinate2.$strand$extension
  Length=$( echo $singleline$file_base_name$coordinate2.$strand$extension | sed -e '1d' $singleline$file_base_name$coordinate2.$strand$extension | wc -c )
  STOPh=$( expr $START + $Length - 1 )

  rm $singleline$file_base_name$coordinate2.$strand$extension

  sed -i -e "s/>/\>Coordinate: $coordinate1. Strand: $strand. Range: $START:$STOP Length: /g" $file_base_name$coordinate2.$strand$extension
  fi

fi

done

#It concatenates all files that contain "minus" in the name.
cat *.minus.fa > in.reverse.fa
#Reverse complement.
seqkit seq -r -p -t DNA in.reverse.fa > out.reverse.fa
#It concatenates all files with sequences extracted to the flankings.fa file.
cat *.plus.fa out.reverse.fa > flankings.fa

#It deletes intermediate files.
rm -r *.plus.fa *.minus.fa *.reverse.fa


#It removes the jonction sequences extracted from all the chromosomes except the one in which the transposon is identified.
awk -F "=|:" '/^>/ {s=$8 ".flankings.fa"}; {print > s}' flankings.fa
for i in *.flankings.fa
do
if ! [[ $i == $loc.flankings.fa ]]
then
rm -r $i
fi
done
 
#Because there were several sequences in the previous file, the total length was miscalculated.
#It calculates the length of the extracted sequence correctly.
for i in *.flankings.fa
do
awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' $i > $i.singleline 
Length1=$( echo $i.singleline | sed -e '1d' $i.singleline | wc -c )
#It appends the correct length to the header.
sed -i -e "s/Length: /\ Length: $Length1 /g " $i
done 
rm -r *.singleline flankings.fa

#Implemented step for sorting results using BLAST.
for i in *.flankings.fa
do
blastn -query $i -word_size 13 -db ../../../pipeline_extragere/Genome/drosophilidae -out ${i%.fa}.drosophilidae.blast.txt -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sstrand"
done

average=$( expr $LengthFL2 + $LengthFL1 / 2 )

awk -v average2="$average" '{ if ( $2 ~/^ *dana/ && $12 >= average2 ) print $0;}' *.flankings.drosophilidae.blast.txt > dana.results.txt
awk -v average2="$average" '{ if ( $2 ~/^ *dere/ && $12 >= average2 ) print $0;}' *.flankings.drosophilidae.blast.txt > dere.results.txt
awk -v average2="$average" '{ if ( $2 ~/^ *dgri/ && $12 >= average2 ) print $0;}' *.flankings.drosophilidae.blast.txt > dgri.results.txt
awk -v average2="$average" '{ if ( $2 ~/^ *dmoj/ && $12 >= average2 ) print $0;}' *.flankings.drosophilidae.blast.txt > dmoj.results.txt
awk -v average2="$average" '{ if ( $2 ~/^ *dper/ && $12 >= average2 ) print $0;}' *.flankings.drosophilidae.blast.txt > dper.results.txt
awk -v average2="$average" '{ if ( $2 ~/^ *dpse/ && $12 >= average2 ) print $0;}' *.flankings.drosophilidae.blast.txt > dpse.results.txt
awk -v average2="$average" '{ if ( $2 ~/^ *dsec/ && $12 >= average2 ) print $0;}' *.flankings.drosophilidae.blast.txt > dsec.results.txt
awk -v average2="$average" '{ if ( $2 ~/^ *dsim/ && $12 >= average2 ) print $0;}' *.flankings.drosophilidae.blast.txt > dsim.results.txt
awk -v average2="$average" '{ if ( $2 ~/^ *dvir/ && $12 >= average2 ) print $0;}' *.flankings.drosophilidae.blast.txt > dvir.results.txt
awk -v average2="$average" '{ if ( $2 ~/^ *dwil/ && $12 >= average2 ) print $0;}' *.flankings.drosophilidae.blast.txt > dwil.results.txt
awk -v average2="$average" '{ if ( $2 ~/^ *dyak/ && $12 >= average2 ) print $0;}' *.flankings.drosophilidae.blast.txt > dyak.results.txt

mkdir BLAST
mv *.drosophilidae.blast.txt *.results.txt ./BLAST/
#It deletes the empty files.
find . /tmp -size 0 -delete 

#It appends the remnant name to the header.
for i in *.flankings.fa 
do
awk -v dirname2="$dirname" '/^>/{sub(">", ">" dirname2 )}1' $i > ${i%.flankings.fa}.flankings.fasta
done
rm -r *.flankings.fa

 cd -
done

#It merges all the files which contain junction sequences.
shopt -s globstar | cat **/*.flankings.fasta > FLANKINGS.fa | shopt -u globstar
 
#It performs BLAST with XML results to be loaded into Kablammo.
blastn -query FLANKINGS.fa -word_size 13 -db ../../pipeline_extragere/Genome/drosophilidae -outfmt 5 -out KABLAMMO.xml
