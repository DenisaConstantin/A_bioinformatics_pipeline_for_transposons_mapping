#!/bin/bash

singleline2=".singleline"
echo "Folder"
read folder

echo "Remnant"
read remnant

#It splits the MULTIFASTA file into FASTA files, using the name of the remnant and its location (the chromosome in which it is identified).
awk -v folder2="$folder" -F "=|:|;" '/^>/ {s=folder2$7"-"$4 ".fa"}; {print > s}' dmel-all-transposon-r6.39.fasta

#It changes the folder to the user-specified working folder.
cd $folder

#It browses the files and it deletes the files whose name does not start with the user-specified remnant name (for hobo-H{}).
for i in *.fa
 do
 if ! [[ $i == "$remnant"* ]]
 then
 rm -r $i
 fi
done

#It deletes the files that have less than 30 nucleotides.
for i in *.fa

do

#It transforms the files into single-line ones so you can count nucleotides.
awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' $i > $i$singleline2

#It counts the nucleotides.
Length=$( sed -e '1d' $i$singleline2 | wc -c )
if [ $Length -lt 30 ]
then
rm -r $i
fi
rm -r $i$singleline2
done

cd ..

#It aligns each remnant with the genome of Drosophila melanogaster.
for i in $folder*.fa

do

#It performs the alignment and returns the result in tabular form.
blastn -query $i -word_size 13 -db ../pipeline_extragere/Genome/dmel -out ${i%.fasta}.blast.txt -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sstrand"
done 

#It creates a folder for each .txt it finds.
for i in $folder*.fa.blast.txt

do   
  
mkdir ${i%.fa.blast.txt}
 
done

#It moves each .txt file to the homonymous folder.
for i in $folder*.txt

do 

directory=${i%.fa.blast.txt}   
mv "$i" "$directory"

done 

cd $folder

#It moves .fa to the folder of the same name.
for i in *.fa

do

directory=${i%.fa}
mv "$i" "$directory"

done
