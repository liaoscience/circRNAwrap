
#conditions chose,  acfs 2, find_circ 2, circRNA_finder 2, mapsplice 2, CIRI2 2, CIRCexplorer 2, KNIFE 2. 

#for sample in SRR1049826 SRR1049827 SRR1049828 SRR1049829 SRR1049830 SRR1049831 SRR1049832 SRR1049833;do
for sample in circRNA linear mix;do
cd /d/work_dir/OneDrive/circRNA_tools_compare/result/raw_result
cd ./$sample""_test
rm -r circRNA_validate
mkdir circRNA_validate

########################################################################### 1 several tools together
#acfs 2
awk -F '[|\t]' '{if(/^chr/ && $7>=2 ) print $1"\t"$2"\t"$3"\t"$4"_"$8"\t"$7"\t"$9}' ./$sample""_acfs/circle_candidates_CBR.bed12 > ./circRNA_validate/$sample""_CBR.txt
awk -F '[|\t]' '{if(/^chr/ && $7>=2 ) print $1"\t"$2"\t"$3"\t"$4"_"$8"\t"$7"\t"$9}' ./$sample""_acfs/circle_candidates_MEA.bed12 > ./circRNA_validate/$sample""_MEA.txt
cat ./circRNA_validate/$sample""_MEA.txt ./circRNA_validate/$sample""_CBR.txt > ./circRNA_validate/$sample""_acfs.txt
rm ./circRNA_validate/$sample""_CBR.txt ./circRNA_validate/$sample""_MEA.txt
sed -i '1i\chr\tstart\tend\tsample\tacfs\tstrand' ./circRNA_validate/$sample""_acfs.txt

#find_circ 2
awk '{if($7>=2 && !/norm/) print $1"\t"$2"\t"$3"\t"$1"_"$2"_"$3"\t"$5"\t"$6}' ./$sample""_fc/$sample.sites.bed > ./circRNA_validate/$sample.find_circ.txt #unique mapped 2
sed -i '1d' ./circRNA_validate/$sample.find_circ.txt
sed -i '1i\chr\tstart\tend\tsample\tfind_circ\tstrand' ./circRNA_validate/$sample.find_circ.txt

#circRNA_finder STAR 2
awk '{if($5>=2) print $1"\t"$2"\t"$3"\t"$1"_"$2"_"$3"\t"$5"\t"$6}' ./$sample.circRNA_finder.bed > ./circRNA_validate/$sample.circRNA_finder.txt #unique mapped 3
sed -i '1i\chr\tstart\tend\tsample\tcircRNA_finder\tstrand' ./circRNA_validate/$sample.circRNA_finder.txt

#mapsplice 2
awk -F '[~\t]' '{ if($1==$2 && $3>$4 && $3<=$4+1000000 && $6>=2 && ($7~/\+\+/ || $7~/\-\-/)) print $1"\t"$4-1"\t"$3"\t"$1"_"$4-1"_"$3"\t"$6"\t"$7; if($1==$2 && $3<$4 && $3>=$4-1000000 && $6>=2&& ($7~/\+\+/ || $7~/\-\-/)) print $1"\t"$3-1"\t"$4"\t"$1"_"$3-1"_"$4"\t"$6"\t"$7}' ./fusions_candidates.txt | sed 's/++/+/g' | sed 's/--/-/g' > ./circRNA_validate/$sample.mapsplice.txt
sed -i '1i\chr\tstart\tend\tsample\tmapsplice\tstrand' ./circRNA_validate/$sample.mapsplice.txt

#KNIFE 2
awk '{if($5>=0.9 && /rev/) print $1"\t"$3"\t"$2"\t"$4"\t"$5"\t"$6"\t"$8}' ./combinedReports/${sample}_1__circJuncProbs.txt | awk -F '[:|\t]' '{print $1"\t"$3"\t"$5"\t"$1"_"$3"_"$5"\t"$12"\t"$7}' | sed '1d' | awk '{if($2>$3) print $1"\t"$3"\t"$2"\t"$4"\t"$5"\t"$6; if($2<$3) print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' > ./circRNA_validate/$sample.KNIFE.txt

awk '{if($4>0 && /rev/) print $1"\t"$3"\t"$2"\t"$4"\t"$5"\t"$6"\t"$8}' ./combinedReports/naive${sample}_1_report.txt | awk -F '[:|\t]' '{print $1"\t"$3"\t"$5"\t"$1"_"$3"_"$5"\t"$8"\t"$7}' | sed '1d' | awk '{if($2>$3) print $1"\t"$3"\t"$2"\t"$4"\t"$5"\t"$6; if($2<$3) print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' >> ./circRNA_validate/$sample.KNIFE.txt
sed -i '1i\chr\tstart\tend\tsample\tKNIFE\tstrand' ./circRNA_validate/$sample.KNIFE.txt

#CIRI2
awk '{print $2"\t"$3-1"\t"$4"\t"$2"_"$3-1"_"$4"\t"$5"\t"$11}' ./${sample}_CIRI/${sample}.CIRI.txt | sed '1d' > ./circRNA_validate/$sample.CIRI2.txt
awk '{print $2"\t"$3-1"\t"$4"\t"$2"_"$3-1"_"$4"\t"$5"\t"$11}' ./${sample}.circRNA | sed '1d' > ./circRNA_validate/$sample.CIRI2.txt
sed -i '1i\chr\tstart\tend\tsample\tCIRI\tstrand' ./circRNA_validate/$sample.CIRI2.txt

#CIRCexplorer
awk '{if($13>=2) print $1"\t"$2"\t"$3"\t"$1"_"$2"_"$3"\t"$13"\t"$6}' ./${sample}_CIRCexplorer/${sample}_circ.txt > ./circRNA_validate/$sample.CIRCexplorer.txt
awk '{if($13>=2) print $1"\t"$2"\t"$3"\t"$1"_"$2"_"$3"\t"$13"\t"$6}' ./CIRCexplorer/${sample}_circ.txt > ./circRNA_validate/$sample.CIRCexplorer.txt
sed -i '1i\chr\tstart\tend\tsample\tCIRCexplorer\tstrand' ./circRNA_validate/$sample.CIRCexplorer.txt

#CIRCexplorer2 
awk '{if($13>=2) print $1"\t"$2"\t"$3"\t"$1"_"$2"_"$3"\t"$13"\t"$6}' ./${sample}.CIRCexplorer2.txt > ./circRNA_validate/$sample.CIRCexplorer2.txt
sed -i '1i\chr\tstart\tend\tsample\tCIRCexplorer\tstrand' ./circRNA_validate/$sample.CIRCexplorer2.txt

#DCC 
awk '{if($13>=2) print $1"\t"$2"\t"$3"\t"$1"_"$2"_"$3"\t"$13"\t"$6}' ./${sample}.CIRCexplorer2.txt > ./circRNA_validate/$sample.CIRCexplorer2.txt
sed -i '1i\chr\tstart\tend\tsample\tCIRCexplorer\tstrand' ./circRNA_validate/$sample.CIRCexplorer2.txt

cat ./circRNA_validate/*.txt | awk '!a[$1"\t"$2"\t"$3]++' > ./circRNA_validate/$sample.7_tools.txt
awk '{if(!/start/) print }' ./circRNA_validate/$sample.7_tools.txt > ./circRNA_validate/$sample.7_tools.bed
awk '{i++; print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"i }' ./circRNA_validate/$sample.7_tools.bed > ./circRNA_validate/rest.7_tools.bed

done
