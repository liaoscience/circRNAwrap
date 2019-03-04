################# content, detected circRNAs, 
################# use the backsplice site to realign, abundance, fraction, unique, transcripts
################# 1 merge all tools output
################# 2 build index
################# 3 realign the unmapped reads and get the abundance reads, meanwhile, use sailfish for abundance estimate
################# 4 remove duplication reads and get the abundance
################# 5 circRNA backsplice site cis splice reads abundance
################# 6 circRNA backsplice site best linear transcripts and circRNA fasta sequence
################# 7 circRNA transcripts prediction and support paired-end reads
################# 8 linear RNA abundance analysis

#add the samples
sample=$1
threads=$2
dir=/home/lilin/workdir/data/circRNA/data
cd $dir/${sample}

#software
put_together_ws="/home/lilin/workdir/data/circRNA/put_together_ws.pl"
bedtools=/home/lilin/miniconda2/bin/bedtools
hisat2=/home/lilin/miniconda2/bin/hisat2
hisat2-build=/home/lilin/miniconda2/bin/hisat2-build
samtools="/home/lilin/workdir/git/samtools-1.4/samtools"
mosdepth=/home/lilin/miniconda2/bin/mosdepth

#reference
reference=/home/lilin/workdir/reference/gatk4/hg19_ref/annotation
genome=/home/lilin/workdir/reference/gatk4/hg19_ref/ucsc.hg19.fasta
GTF=/home/lilin/workdir/reference/gatk4/hg19_ref/hg19_genes.gtf

# 1 merge all tools output

#conditions chose,  acfs 2, find_circ 2, circRNA_finder 2, mapsplice 2, CIRI2 2, CIRCexplorer 2, KNIFE 2, DCC, 2

mkdir circRNA_validate

########################################################################### 1 several tools together

#KNIFE 2
awk '{if($5>=0.9 && /rev/) print $1"\t"$3"\t"$2"\t"$4"\t"$5"\t"$6"\t"$8}' ./${sample}_KNIFE/combinedReports/${sample}_1__circJuncProbs.txt | awk -F '[:|\t]' '{print $1"\t"$3"\t"$5"\t"$1"_"$3"_"$5"\t"$12"\t"$7}' | sed '1d' | awk '{if($2>$3) print $1"\t"$3"\t"$2"\t"$4"\t"$5"\t"$6; if($2<$3) print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' > ./circRNA_validate/${sample}.KNIFE.txt
#awk '{if($4>0 && /rev/) print $1"\t"$3"\t"$2"\t"$4"\t"$5"\t"$6"\t"$8}' ./combinedReports/naive${sample}_1_report.txt | awk -F '[:|\t]' '{print $1"\t"$3"\t"$5"\t"$1"_"$3"_"$5"\t"$8"\t"$7}' | sed '1d' | awk '{if($2>$3) print $1"\t"$3"\t"$2"\t"$4"\t"$5"\t"$6; if($2<$3) print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' >> ./circRNA_validate/${sample}.KNIFE.txt
sed -i '1i\chr\tstart\tend\tsample\tKNIFE\tstrand' ./circRNA_validate/${sample}.KNIFE.txt

#find_circ 2
awk '{if($7>=2 && !/norm/) print $1"\t"$2"\t"$3"\t"$1"_"$2"_"$3"\t"$5"\t"$6}' ./${sample}""_find_circ/${sample}.sites.bed > ./circRNA_validate/${sample}.find_circ.txt #unique mapped 2
sed -i '1d' ./circRNA_validate/${sample}.find_circ.txt
sed -i '1i\chr\tstart\tend\tsample\tfind_circ\tstrand' ./circRNA_validate/${sample}.find_circ.txt

#acfs 2
awk -F '[|\t]' '{if(/^chr/ && $7>=2 ) print $1"\t"$2"\t"$3"\t"$4"_"$8"\t"$7"\t"$9}' ./${sample}""_acfs/circle_candidates_CBR.bed12 > ./circRNA_validate/${sample}""_CBR.txt
awk -F '[|\t]' '{if(/^chr/ && $7>=2 ) print $1"\t"$2"\t"$3"\t"$4"_"$8"\t"$7"\t"$9}' ./${sample}""_acfs/circle_candidates_MEA.bed12 > ./circRNA_validate/${sample}""_MEA.txt
cat ./circRNA_validate/${sample}""_MEA.txt ./circRNA_validate/${sample}""_CBR.txt > ./circRNA_validate/${sample}""_acfs.txt
rm ./circRNA_validate/${sample}""_CBR.txt ./circRNA_validate/${sample}""_MEA.txt
sed -i '1i\chr\tstart\tend\tsample\tacfs\tstrand' ./circRNA_validate/${sample}""_acfs.txt

#CIRI 2
awk '{if($5>2) print $2"\t"$3-1"\t"$4"\t"$2"_"$3-1"_"$4"\t"$5"\t"$11}' ./${sample}_CIRI/${sample}.CIRI.txt | sed '1d' > ./circRNA_validate/${sample}.CIRI2.txt
sed -i '1i\chr\tstart\tend\tsample\tCIRI\tstrand' ./circRNA_validate/${sample}.CIRI2.txt

#CIRCexplorer 2
awk '{if($13>=2) print $1"\t"$2"\t"$3"\t"$1"_"$2"_"$3"\t"$13"\t"$6}' ./${sample}_CIRCexplorer2/${sample}_circ.txt > ./circRNA_validate/${sample}.CIRCexplorer2.txt
sed -i '1i\chr\tstart\tend\tsample\tCIRCexplorer\tstrand' ./circRNA_validate/${sample}.CIRCexplorer2.txt

#circRNA_finder STAR 2
awk '{if($5>=2) print $1"\t"$2"\t"$3"\t"$1"_"$2"_"$3"\t"$5"\t"$6}' ./${sample}_circRNA_finder/${sample}.circRNA_finder.bed > ./circRNA_validate/${sample}.circRNA_finder.txt #unique mapped 3
sed -i '1i\chr\tstart\tend\tsample\tcircRNA_finder\tstrand' ./circRNA_validate/${sample}.circRNA_finder.txt

#mapsplice 2
awk -F '[~\t]' '{ if($1==$2 && $3>$4 && $3<=$4+100000 && $6>=2 && ($7~/\+\+/ || $7~/\-\-/)) print $1"\t"$4-1"\t"$3"\t"$1"_"$4-1"_"$3"\t"$6"\t"$7; if($1==$2 && $3<$4 && $3>=$4-100000 && $6>=2&& ($7~/\+\+/ || $7~/\-\-/)) print $1"\t"$3-1"\t"$4"\t"$1"_"$3-1"_"$4"\t"$6"\t"$7}' ./${sample}_mapsplice/fusions_candidates.txt | sed 's/++/+/g' | sed 's/--/-/g' > ./circRNA_validate/${sample}.mapsplice.txt
sed -i '1i\chr\tstart\tend\tsample\tmapsplice\tstrand' ./circRNA_validate/${sample}.mapsplice.txt

#DCC 2
awk '{if($5>=2) print $1"\t"$2"\t"$3"\t"$1"_"$2"_"$3"\t"$13"\t"$6}' ./${sample}_DCC/CircCoordinates > ./circRNA_validate/${sample}.DCC.txt
sed -i '1i\chr\tstart\tend\tsample\tCIRCexplorer\tstrand' ./circRNA_validate/${sample}.DCC.txt


cat ./circRNA_validate/*.txt | awk '!a[$1"\t"$2"\t"$3]++' > ./circRNA_validate/${sample}.tools.txt
awk '{if(!/start/) print }' ./circRNA_validate/${sample}.tools.txt > ./circRNA_validate/${sample}1.tools.bed
awk '{i++; print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"i }' ./circRNA_validate/${sample}1.tools.bed > ./circRNA_validate/${sample}.tools.bed



########################################################################### 2 build index

cd ./circRNA_validate
# short circRNA
awk '{if($3-$2<=125) print }' ${sample}.circ.bed > short125.bed
#reference 
fastaFromBed -fi $reference -bed short125.bed -fo circRNA_short.fa -s
#repeat quadrupling
paste circRNA_short.fa circRNA_short.fa circRNA_short.fa circRNA_short.fa > short.fa

# long circRNA
awk '{if($3-$2>125) print }' ${sample}.circ.bed > long125.bed
#reference
fastaFromBed -fi $reference -bed long125.bed -fo circRNA_long.fa -s
#repeat twice
paste circRNA_long.fa circRNA_long.fa > long.fa
cat long.fa short.fa > total.fa
rm circRNA_long.fa circRNA_short.fa short.fa long.fa
mkdir index
mv total.fa ./index/
cd ./index/
sed 's/\t>.*//g' total.fa | sed 's/\t//g' > total.repeat1.fa
fold total.repeat1.fa > total.repeat.fa
samtools faidx total.repeat.fa
rm total.repeat1.fa
rm total.fa
${hisat2-build} -p $threads total.repeat.fa total.repeat

########################################################################### 3 RAISE abundance
cd $dir/${sample}
awk 'NR%4==1' ${sample}.unmap.fq > unmapped.title
$bedtools bamtobed -split -bed12 -i ${sample}.hisat2.bam > accepted_hits.bed
awk -F '[@_\t/]' 'NR==FNR{a[$2]=$0;next}NR>FNR{if($4 in a)print $0}' unmapped.title accepted_hits.bed > unmapped.mapped.bed
rm accepted_hits.bed
rm unmapped.title
$hisat2 -p $threads --mm --dta -x ./circRNA_validate/index/total.repeat -U ${sample}.unmap.fq -S ./circRNA_validate/${sample}.unmap.sam 2> ${sample}.circmap.rate 
$samtools sort -@ $threads -o ./circRNA_validate/${sample}.unmap.bam ./circRNA_validate/${sample}.unmap.sam
$samtools rmdup -S ./circRNA_validate/${sample}.unmap.bam --reference ./circRNA_validate/index/total.repeat ./circRNA_validate/rmdup.bam
$bedtools bamtobed -split -bed12 -i ./circRNA_validate/rmdup.bam > ./circRNA_validate/rmdup.bed
#sed -i 's/__/\//g' ./circRNA_validate/rmdup.bed
$bedtools bamtobed -split -bed12 -i ./circRNA_validate/${sample}.unmap.bam > ./circRNA_validate/accepted_hits.bed
sed -i 's/__/\//g' ./circRNA_validate/accepted_hits.bed
$bedtools intersect -split -wa -wb -abam ./circRNA_validate/${sample}.unmap.bam -b ./circRNA_validate/index/ref2bp.bed -bed | awk '!a[$0]++' > circ_splice.reads
sed -i 's/__/\//g' circ_splice.reads   #get whole circRNA splice reads in the splice site
awk -F '[_\t/]' 'NR==FNR{a[$4]=a[$4]$0" ";next}NR>FNR{if($4 in a)print $0"\t"a[$4]}' unmapped.mapped.bed circ_splice.reads > tmp.pair.circ_splice.reads
awk '{split($1,a,"[-:()/+]"); if(a[1]==$19 && a[2]-20 <= $20 && a[3]+20 >= $21 && $4!=$22 ) print }' tmp.pair.circ_splice.reads > tmp.pair_proper.circ_splice.reads
awk 'ARGIND==1{a[$0]}ARGIND>1&&!($0 in a ){print $0}' tmp.pair_proper.circ_splice.reads tmp.pair.circ_splice.reads > tmp.pair_unproper.circ_splice.reads
#count
#both paired reads in circRNA sequence and they are not align to genome
awk -F '[_\t]' 'NR==FNR{a[$1"\t"$4]=$0;next}NR>FNR{if($1"\t"$4 in a)print a[$1"\t"$4]"\t"$0}' circ_splice.reads ./circRNA_validate/accepted_hits.bed | awk -F '[_\t]' '{if($4==$24 && $6!~$26) print }' > tmp.both_circ_splice.reads
awk '{a[$4]++;b[$4]=$0}END{for(i in a){print i,a[i] | "sort -k 1"}}' tmp.both_circ_splice.reads > tmp.both_circ_splice.reads1
awk 'NR==FNR{a[$1]=$2;next}NR>FNR{if($4 in a)print a[$4]"\t"$0}' tmp.both_circ_splice.reads1 circ_splice.reads > tmp.both_circ_splice.reads2
awk '{a[$2]+=1/$1;}END{for(i in a){print i,a[i] | "sort -k 1"}}' tmp.both_circ_splice.reads2 > ${sample}.both.txt
sed -i '1i\circrna both' ${sample}.both.txt
#count total 
awk '{a[$4]++;b[$4]=$0}END{for(i in a){print i,a[i] | "sort -k 1"}}' circ_splice.reads > tmp.circ_validate.reads
awk 'NR==FNR{a[$1]=$2;next}NR>FNR{if($4 in a)print a[$4]"\t"$0}' tmp.circ_validate.reads circ_splice.reads > tmp.circ_splice.reads2
awk '{a[$2]+=1/$1;}END{for(i in a){print i,a[i] | "sort -k 1"}}' tmp.circ_splice.reads2 > ${sample}.total.txt
sed -i '1i\circrna '${sample}'' ${sample}.total.txt
#count proper paired
awk '{a[$4]++;b[$4]=$0}END{for(i in a){print i,a[i] | "sort -k 1"}}' tmp.pair_proper.circ_splice.reads > tmp.circ_validate.reads1
awk 'NR==FNR{a[$1]=$2;next}NR>FNR{if($4 in a)print a[$4]"\t"$0}' tmp.circ_validate.reads1 tmp.pair_proper.circ_splice.reads > tmp.circ_splice.reads21
awk '{a[$2]+=1/$1;}END{for(i in a){print i,a[i] | "sort -k 1"}}' tmp.circ_splice.reads21 > ${sample}.pp.txt
sed -i '1i\circrna pp' ${sample}.pp.txt
# count unproper paired
awk '{a[$4]++;b[$4]=$0}END{for(i in a){print i,a[i] | "sort -k 1"}}' tmp.pair_unproper.circ_splice.reads > tmp.circ_unvalidate.reads1
awk 'NR==FNR{a[$1]=$2;next}NR>FNR{if($4 in a)print a[$4]"\t"$0}' tmp.circ_unvalidate.reads1 tmp.pair_unproper.circ_splice.reads > tmp.circ_unsplice.reads21
awk '{a[$2]+=1/$1;}END{for(i in a){print i,a[i] | "sort -k 1"}}' tmp.circ_unsplice.reads21 > ${sample}.up.txt
sed -i '1i\circrna up' ${sample}.up.txt
mkdir tmp
cp ${sample}.total.txt ./tmp
mv ${sample}.both.txt ./tmp
mv ${sample}.up.txt ./tmp
mv ${sample}.pp.txt ./tmp
perl $put_together_ws ./tmp | sort -r -k2 > ${sample}.circ.txt
rm tmp.*
rm -r ./tmp

################# 4 remove PCR duplication reads and get the abundance
cd $dir/${sample}
$bedtools intersect -split -wa -wb -abam ./circRNA_validate/rmdup.bam -b ./circRNA_validate/index/ref2bp.bed -bed | awk '!a[$0]++' > unique.circ_splice.reads     
sed -i 's/__/\//g' unique.circ_splice.reads   #get whole circRNA splice reads in the splice site
awk -F '[_\t/]' 'NR==FNR{a[$4]=a[$4]$0" ";next}NR>FNR{if($4 in a)print $0"\t"a[$4]}' unmapped.mapped.bed unique.circ_splice.reads > tmp.pair.unique.circ_splice.reads
awk '{split($1,a,"[-:()/+]"); if(a[1]==$19 && a[2]-20 <= $20 && a[3]+20 >= $21 && $4!=$22 ) print }' tmp.pair.unique.circ_splice.reads > tmp.pair_proper.unique.circ_splice.reads
awk 'ARGIND==1{a[$0]}ARGIND>1&&!($0 in a ){print $0}' tmp.pair_proper.unique.circ_splice.reads tmp.pair.unique.circ_splice.reads > tmp.pair_unproper.unique.circ_splice.reads
#count
#both paired reads in circRNA sequence and they are not align to genome
awk -F '[_\t]' 'NR==FNR{a[$1"\t"$4]=$0;next}NR>FNR{if($1"\t"$4 in a)print a[$1"\t"$4]"\t"$0}' unique.circ_splice.reads ./circRNA_validate/rmdup.bed | awk -F '[_\t]' '{if($4==$24 && $6!~$26) print }' > tmp.both_unique.circ_splice.reads
awk '{a[$4]++;b[$4]=$0}END{for(i in a){print i,a[i] | "sort -k 1"}}' tmp.both_unique.circ_splice.reads > tmp.both_unique.circ_splice.reads1
awk 'NR==FNR{a[$1]=$2;next}NR>FNR{if($4 in a)print a[$4]"\t"$0}' tmp.both_unique.circ_splice.reads1 unique.circ_splice.reads > tmp.both_unique.circ_splice.reads2
awk '{a[$2]+=1/$1;}END{for(i in a){print i,a[i] | "sort -k 1"}}' tmp.both_unique.circ_splice.reads2 > ${sample}.both.txt
sed -i '1i\circrna both' ${sample}.both.txt
#count total 
awk '{a[$4]++;b[$4]=$0}END{for(i in a){print i,a[i] | "sort -k 1"}}' unique.circ_splice.reads > tmp.circ_validate.reads
awk 'NR==FNR{a[$1]=$2;next}NR>FNR{if($4 in a)print a[$4]"\t"$0}' tmp.circ_validate.reads unique.circ_splice.reads > tmp.unique.circ_splice.reads2
awk '{a[$2]+=1/$1;}END{for(i in a){print i,a[i] | "sort -k 1"}}' tmp.unique.circ_splice.reads2 > ${sample}.unique.total.txt
sed -i '1i\circrna '${sample}'' ${sample}.unique.total.txt
#count proper paired
awk '{a[$4]++;b[$4]=$0}END{for(i in a){print i,a[i] | "sort -k 1"}}' tmp.pair_proper.unique.circ_splice.reads > tmp.circ_validate.reads1
awk 'NR==FNR{a[$1]=$2;next}NR>FNR{if($4 in a)print a[$4]"\t"$0}' tmp.circ_validate.reads1 tmp.pair_proper.unique.circ_splice.reads > tmp.unique.circ_splice.reads21
awk '{a[$2]+=1/$1;}END{for(i in a){print i,a[i] | "sort -k 1"}}' tmp.unique.circ_splice.reads21 > ${sample}.pp.txt
sed -i '1i\circrna pp' ${sample}.pp.txt
# count unproper paired
awk '{a[$4]++;b[$4]=$0}END{for(i in a){print i,a[i] | "sort -k 1"}}' tmp.pair_unproper.unique.circ_splice.reads > tmp.circ_unvalidate.reads1
awk 'NR==FNR{a[$1]=$2;next}NR>FNR{if($4 in a)print a[$4]"\t"$0}' tmp.circ_unvalidate.reads1 tmp.pair_unproper.unique.circ_splice.reads > tmp.circ_unsplice.reads21
awk '{a[$2]+=1/$1;}END{for(i in a){print i,a[i] | "sort -k 1"}}' tmp.circ_unsplice.reads21 > ${sample}.up.txt
sed -i '1i\circrna up' ${sample}.up.txt
mkdir tmp
cp ${sample}.unique.total.txt ./tmp
mv ${sample}.both.txt ./tmp
mv ${sample}.up.txt ./tmp
mv ${sample}.pp.txt ./tmp
perl $put_together_ws ./tmp | sort -r -k2 > ${sample}.unique.circ.txt
rm tmp.*
rm -r ./tmp


################# 5 circRNA backsplice site cis splice reads abundance
cd $dir/${sample}
awk -v sample="${sample}" 'BEGIN{print "circrna\t"sample"\tboth\tpp\tup"}NR==1{for(i=0;i++<NF;)a[$i]=i;next}{print $a["circrna"]"\t"$a[sample]"\t"$a["both"]"\t"$a["pp"]"\t"$a["up"]}' ${sample}.circ.txt > ${sample}.tmp.circ.txt
cat ${sample}.tmp.circ.txt | sed 's/(-)/ nega/g' | sed 's/(+)/ posi/g' | awk -F '[ \t:-]' '{print $1,$2,$3,$0}' | sed 's/nega/-/g' | sed 's/posi/+/g' | sed 's/ /\t/g' | awk '{print $1"\t"$2"\t"$3"\t"$4"("$5")\t"$6"\t"$5}' | sed '1d' | sort -k 1,1 -k2,2n > ${sample}.circ1.bed
sortBed -i ${sample}.circ1.bed > ${sample}.circ.bed

awk '{print $1,$2-1,$2+1,$4,"up",$6"\n"$1,$3-1,$3+1,$4,"down",$6 }' ${sample}.circ.bed | sed 's/ /\t/g' > ${sample}.tmp.2.circ.bed
$mosdepth -t $threads --by ${sample}.tmp.2.circ.bed ${sample}.tmp.mos ${sample}.hisat2.rmdup.bam
gunzip ${sample}.tmp.mos.regions.bed.gz
awk '{a[$4]=a[$4](a[$4]?"\t":"")$5;}END{for (j in a) print j"\t"a[j]}' ${sample}.tmp.mos.regions.bed | sed '1i\circrna\tleft\tright' > ${sample}.tmp.mos
awk 'NR==FNR{a[$1]=$2"\t"$3;next}NR>FNR{if($1 in a)print $0"\t"a[$1]}' ${sample}.tmp.mos ${sample}.tmp.circ.txt > ${sample}.circ.fraction

################# 6 circRNA backsplice site annotation
cd $dir/${sample}

#annotation
# get circRNA inner exon sequence
closestBed -s -a ${sample}.circ.bed -b $anno/gencode.v19.annotation.refflat.gtf -d | awk '{ if($NF==0) print }' > tmp.1  # circ and exon same as acfs gtf file, inner exon of circRNAs
awk '{print $7"\t"$10"\t"$11"\t"$15"\t"$12"\t"$13}' tmp.1 | awk '{$2=$2; print }' | sed 's/ /\t/g' > tmp.2 #circ related exon
fastaFromBed -name $4 -s -tab -fi $genome -bed tmp.2 -fo tmp.3 #circ related exon fasta
awk -F '[\t\(]' 'NR==FNR{a[$1]=$0;next}NR>FNR{if($16 in a)print $0"\t"a[$16]}' tmp.3 tmp.1 | sort -n -k17 > tmp.4 #match the exon to related circ acording to the transcript and exon name
sed 's/___/\t/g' tmp.4 | awk -F "[@\t]" '{b[$1,$2,$3,$4,$5,$6,$14,$15]=b[$1,$2,$3,$4,$5,$6,$14,$15]$NF""; c[$1,$2,$3,$4,$5,$6,$14,$15]=c[$1,$2,$3,$4,$5,$6,$14,$15]$(NF-2)"___";}END{for(i in b){split(i,m,SUBSEP); len=length(b[i]); print ">"m[1]"\t"m[2]"\t"m[3]"\t"m[4]"\t"m[5]"\t"m[6]"\t"m[7]"\t"m[8]"\t"len"\t"c[i]"\t"b[i]}}' | sort -nr -k9 | awk '{if(($3-$2)>=$9) print }' | awk '!a[$4]++' > ${sample}.tmp.tab2 #circRNA related exon length  make sure each circRNA are output according this condition

sed 's/___/\t/g' tmp.4 | awk -F "[@\t]" '{b[$1,$2,$3,$4,$5,$6,$14,$15]=b[$1,$2,$3,$4,$5,$6,$14,$15]$NF""; c[$1,$2,$3,$4,$5,$6,$14,$15]=c[$1,$2,$3,$4,$5,$6,$14,$15]$(NF-2)"___";}END{for(i in b){split(i,m,SUBSEP); len=length(b[i]); print ">"m[1]"\t"m[2]"\t"m[3]"\t"m[4]"\t"m[5]"\t"m[6]"\t"m[7]"\t"m[8]"\t"len"\t"c[i]"\t"b[i]}}' | sort -nr -k9 | awk '!a[$4]++' > ${sample}.tmp.tab1   #keep the longest

awk 'ARGIND==1{a[$4]}ARGIND>1&&!($4 in a ){print $0}' ${sample}.tmp.tab2 ${sample}.tmp.tab1 | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"($3-$2)"\t"$10"\t"$11}' > ${sample}.tmp.tab3   #filter some annotation length longer than circRNA spans

cat ${sample}.tmp.tab3 ${sample}.tmp.tab2 > ${sample}.circ_candidates.tab

sed 's/___/\t/g' tmp.4 | awk -F "[@\t]" '{b[$1,$2,$3,$4,$5,$6,$14,$15]=b[$1,$2,$3,$4,$5,$6,$14,$15]$NF""; c[$1,$2,$3,$4,$5,$6,$14,$15]=c[$1,$2,$3,$4,$5,$6,$14,$15]$(NF-2)"___";}END{for(i in b){split(i,m,SUBSEP); len=length(b[i]); print ">"m[1]"\t"m[2]"\t"m[3]"\t"m[4]"\t"m[5]"\t"m[6]"\t"m[7]"\t"m[8]"\t"len"\t"c[i]"\t"b[i]}}' | sort -nr -k9 | awk '{if(($3-$2)>$9) print }' > ${sample}.circ_candidates_multi.tab

awk 'BEGIN{n=10}{for(i=1;i<n;i++)printf $i"\t";print $i}' ${sample}.circ_candidates.tab | sed 's/>//g' > tmp.anno #circRNA related annotation

#add the utr information
bedtools sort -i tmp.2 > tmp1.2
mv tmp1.2 tmp.2
closestBed -s -a tmp.2 -b $anno/utr.bed -d | awk '{ if($NF==0) print }' | sed 's/NM_/NM/g' > tmp.utr_in_exon # circRNA related exon which located in utr region
sed 's/___/\t/g' tmp.utr_in_exon | awk '{if($4==$17) print }' | awk '!a[$4"\t"$5]++' > tmp.utr_in_exon1
awk 'NR==FNR{a[$4]=$0;next}NR>FNR{if($7 in a)print $0"\t"a[$7]}' tmp.utr_in_exon1 tmp.anno | awk '{if($2 < $20 && $3 < $20 ) next; if($2 > $21 && $3 > $21 ) next; print }' > tmp.utr.anno #utr contain circRNAs
awk '{if($2 > $20 && $3 < $21) print }' tmp.utr.anno > ${sample}.tmp.circ_in_utr # only in utr region
awk 'ARGIND==1{a[$4]}ARGIND>1&&!($4 in a ){print $0}' ${sample}.tmp.circ_in_utr tmp.utr.anno > ${sample}.tmp.utr_cds.anno

#class 1, utr, 2, only exon, 3, not in exon (intron, intergenic, antisense)
awk 'ARGIND==1{a[$4]}ARGIND>1&&!($4 in a ){print $0}' tmp.1 ${sample}.circ.bed > tmp.not_in_exon.bed 
bedtools sort -i tmp.not_in_exon.bed > tmp.not_in_exon.bed1
mv tmp.not_in_exon.bed1 tmp.not_in_exon.bed
closestBed -s -a tmp.not_in_exon.bed -b $anno/gencode.v19.annotation.transcript.bed -d | awk '{ if($NF==0) print }' | awk '!a[$4]++' > ${sample}.tmp.circ_in_intron
awk 'ARGIND==1{a[$4]}ARGIND>1&&!($4 in a ){print $0}' ${sample}.tmp.circ_in_intron tmp.not_in_exon.bed > tmp.circ_not_intron
bedtools sort -i tmp.circ_not_intron > tmp.circ_not_intron1
mv tmp.circ_not_intron1 tmp.circ_not_intron
closestBed -S -a tmp.circ_not_intron -b $anno/gencode.v19.annotation.transcript.bed -d | awk '{ if($NF==0) print }' > ${sample}.circ_in_antisense 
awk 'ARGIND==1{a[$4]}ARGIND>1&&!($4 in a ){print $0}' ${sample}.circ_in_antisense tmp.circ_not_intron > ${sample}.tmp.circ_in_intergenic
awk 'ARGIND==1{a[$4]}ARGIND>1&&!($4 in a ){print $0}' tmp.utr.anno tmp.anno > ${sample}.tmp.only_exon.bed

# genetype only focus mRNA or noncoding RNA, NM or NR

#combine together 
#chr start end name abudance strand unique length splice_length annotation_gene annotation_tranccript type(intergenic antisense intron exon 3utr 5utr)
awk '{ leng=($3-$2); print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"leng"\tNA\tNA\tNA\tintergenic"}' ${sample}.tmp.circ_in_intergenic > ${sample}.tmp.circ_in_intergenic1
awk '{ leng=($3-$2); print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"leng"\tNA\t"$10"\t"$11"\tantisense"}' ${sample}.circ_in_antisense | awk '!a[$4]++' > ${sample}.circ_in_antisense1
awk '{ leng=($3-$2); print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"leng"\tNA\t"$10"\t"$11"\tintron"}' ${sample}.tmp.circ_in_intron > ${sample}.tmp.circ_in_intron1
awk '{ m=split($10, a, "___" ); if(m>2) type="exon intron "; if(m==2) type="exon "; leng=($3-$2); print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"leng"\t"$9"\t"$7"\t"$8"\t"type""$10}' ${sample}.tmp.only_exon.bed > ${sample}.tmp.only_exon.bed1
awk '{ m=split($10, a, "___" ); if(m>2) type="exon intron "; if(m==2) type="exon"; leng=($3-$2); print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"leng"\t"$9"\t"$7"\t"$8"\t"type" "$26" "$10}' ${sample}.tmp.circ_in_utr > ${sample}.tmp.circ_in_utr1
awk '{ m=split($10, a, "___" ); if(m>2) type="exon cds intron"; if(m==2) type="exon cds"; leng=($3-$2); print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"leng"\t"$9"\t"$7"\t"$8"\t"type" "$26" "$10}' ${sample}.tmp.utr_cds.anno > ${sample}.tmp.utr_cds.anno1
cat ${sample}.tmp.circ_in_intergenic1 ${sample}.circ_in_antisense1 ${sample}.tmp.circ_in_intron1 ${sample}.tmp.only_exon.bed1 ${sample}.tmp.utr_cds.anno1 ${sample}.tmp.circ_in_utr1 > ${sample}.annotation.txt

rm *tmp*

################# 7 circRNA transcripts prediction and support paired-end reads
cd $dir/${sample}
#according to the paired end reads, to build the circRNA inner splicing site
awk -F '[_\t/]' 'NR==FNR{a[$4]=a[$4]$0" ";next}NR>FNR{if($4 in a)print $0"\t"a[$4]}' unmapped.mapped.bed circ_splice.reads > tmp.pair.circ_splice.reads
awk '{split($1,a,"[-:()/+]"); if(a[1]==$19 && a[2]-2 <= $20 && a[3]+2 >= $21 && $4!=$22 ) print }' tmp.pair.circ_splice.reads > pair_proper_tmp.circ_splice # get the circRNA paired reads
awk 'BEGIN{n=30}{for(i=1;i<30;i++)printf $i"\t";print $i}' pair_proper_tmp.circ_splice | awk '{$23=$22; $22=$1; print }' | sed 's/ /\t/g' | awk 'BEGIN{n=11}{for(i=(NF-n);i<NF;i++)printf $i"\t";print $i}' > pair_proper_tmp.circ_splice.tmp.1 # convert it to bed12
awk '{if($10==2) {n11=split($11,t11,",");n12=split($12,t12,","); print $1"\t"($7+t11[1])"\t"($8-t11[2])"\t"$4"\t"$5"\t"$6} }' pair_proper_tmp.circ_splice.tmp.1 > pair_proper_tmp.circ_splice.tmp.2 # transform the split information to splice site, 2 split
awk '{if($10==3) {n11=split($11,t11,",");n12=split($12,t12,","); for (i = 1; i < 3; i++) {print $1"\t"($7+t11[i]+t12[i])"\t"($7+t12[i+1])"\t"$4"\t"$5"\t"$6}} }' pair_proper_tmp.circ_splice.tmp.1 > pair_proper_tmp.circ_splice.tmp.3 # transform the split information to splice site, 3 split
#get the circRNA inner splice site
cat pair_proper_tmp.circ_splice.tmp.2 pair_proper_tmp.circ_splice.tmp.3 | awk '!a[$1"\t"$2"\t"$3"\t"$4"\t"$5]++' > pair_proper_tmp.circ_splice.inner  # inner has splice site other there no detected splice site. maybe owning to low abundance
awk '{a[$1"\t"$2"\t"$3"\t"$4]++}END{for(i in a){print i"\t"a[i] | "sort -k 1"}}' pair_proper_tmp.circ_splice.inner | awk '{m=substr($4,length($4)-1,1); print $0"\t"m }' > pair_proper.circ_splice.intro  #inner of circRNA splicing site format same as intron
#test whether the left and right are overlap for each circRNAs positive is origin from right, negative is origin from left
awk '{a[$13]=1; if($24~/\-/ && $26>left[$13] ) left[$13]=$26; if($24~/\+/ && (!length(right[$13]) || $25<right[$13])) right[$13]=$25;}END{for(i in a){print i,right[i],left[i]}}' pair_proper_tmp.circ_splice | awk '{if($2<$3) print }' > pair_proper.circ_splice_spanned
#bulid a transcripts acording to the span circRNAs inner splicing site, detection of splicing site. maybe not all, need high abundance
awk -F '[ \t]' 'NR==FNR{a[$1]=$0;next}NR>FNR{if($4 in a)print $0}' pair_proper.circ_splice_spanned pair_proper.circ_splice.intro | sort -n -k3 | sort -n -k2 | awk '{a[$4]=a[$4](a[$4]?",":"")$2; b[$4]=b[$4](b[$4]?",":"")$3; c[$4]=c[$4](c[$4]?",":"")$5; d[$4]+=1; e[$4]+=$5; }END{for (j in a) print j"\t"d[j]"\t"e[j]"\t"c[j]"\t"a[j]"\t"b[j]}' > ${sample}.circ.transcripts
#compare with intron, extract the intron bed format
#awk '{ n9 = split($9, t9, ",");n10 = split($10, t10, ","); for (i = 0; ++i < n9-1;) { print $2"\t"t10[i]"\t"t9[i + 1]"\t"i "I@" $1"\t"$8"\t"$3 }}' gencode.v24.all.refFlat.txt > gencode.v24.all.refFlat.intron.bed
#awk '{ if($8==1) print $2"\t"$4"\t"$5"\t"$1"\t"$8"\t"$3}' gencode.v24.all.refFlat.txt > gencode.v24.all.refFlat.single.bed
$bedtools intersect -wa -wb -b $reference/gencode.v19.annotation.intron.bed -a pair_proper.circ_splice.intro > pair_proper_tmp.circ_splice.in_intron  #overlap with intron
awk '{if($2==$8 && $3==$9) print }' pair_proper_tmp.circ_splice.in_intron > pair_proper_tmp.circ_splice.in_intron_match
awk 'ARGIND==1{a[$1"\t"$2"\t"$3"\t"$4]}ARGIND>1&&!($1"\t"$2"\t"$3"\t"$4 in a ){print $0}' pair_proper_tmp.circ_splice.in_intron_match pair_proper_tmp.circ_splice.in_intron > pair_proper_tmp.circ_splice.in_intron_mis
awk 'ARGIND==1{a[$1"\t"$2"\t"$3"\t"$4]}ARGIND>1&&!($1"\t"$2"\t"$3"\t"$4 in a ){print $0}' pair_proper_tmp.circ_splice.in_intron pair_proper.circ_splice.intro > pair_proper_tmp.circ_splice.in_not_intron # not overlap with intron
awk '!a[$4]++' pair_proper_tmp.circ_splice.in_intron_mis > pair_proper.circ_splice.in_intron_mis
#not in intron region
awk 'ARGIND==1{a[$4]}ARGIND>1&&!($4 in a ){print $0}' pair_proper_tmp.circ_splice.in_intron pair_proper.circ_splice.intro > pair_proper_tmp.circ_splice.in_out
$bedtools intersect -wa -wb -b $reference/gencode.v19.annotation.single.bed -a pair_proper_tmp.circ_splice.in_out > pair_proper.circ_splice.in_single
awk 'ARGIND==1{a[$4]}ARGIND>1&&!($4 in a ){print $0}' pair_proper.circ_splice.in_single pair_proper_tmp.circ_splice.in_out > pair_proper.circ_splice.intergenic

#extract some discordant alignments
#circRNA bed format
sed 's/(-)/(n)/g' ${sample}.circ.fraction | sed 's/(+)/(p)/g' | awk '{ print $1,$2"_"$4 }' | awk -F '[ :\(\)-]' '{ print $1,$2,$3,$0,$4}' | sed 's/p/+/g' | sed 's/n/-/g' | sed 's/ /\t/g' | sed '1d' > ${sample}.tmp.validate.circ.bed
$bedtools intersect -f 0.99 -split -abam ${sample}.hisat2.bam -b ${sample}.tmp.validate.circ.bed > ${sample}.tmp.validate.circ.bam   #extract the reads in the circRNA region
$samtools view ${sample}.tmp.validate.circ.bam | awk '{if(($2==161 && $9<0) || ($2==81 && $9>0) || ($2==97 && $9>0) || ($2==145 && $9<0) ) print }' > ${sample}.tmp.validate.circ.support #get the support paired circRNAs
awk 'a[$1]++' ${sample}.tmp.validate.circ.support > ${sample}.tmp.validate.circ.support.title #filter some one end reads. kept both paired end reads
awk -F '[\t/]' 'NR==FNR{a[$1]=1;next}NR>FNR{if($1 in a)print $0}' ${sample}.tmp.validate.circ.support.title ${sample}.tmp.validate.circ.support | sort -k 1 > ${sample}.tmp.validate.circ.support1 #get the proper paired reads
$samtools view -H ${sample}.hisat2.bam > tmp.head.sam
cat tmp.head.sam ${sample}.tmp.validate.circ.support1 | $samtools view -bS > ${sample}.tmp.validate.circ.support.bam  # return to bam for further analysis bed split
$bedtools bamtobed -split -bed12 -i ${sample}.tmp.validate.circ.support.bam > ${sample}.tmp.validate.circ.support.bed
$bedtools intersect -wa -wb -a ${sample}.tmp.validate.circ.support.bed -b ${sample}.tmp.validate.circ.bed | awk '{if($2>$14 && $3<$15) print }' > ${sample}.tmp.validate.circ.match
# match the circRNA with related match paired reads 
awk 'ARGIND==1{a[$1]}ARGIND>1&&!($1 in a ){print $0}' ${sample}.tmp.validate.circ.match ${sample}.tmp.validate.circ.support.bed > ${sample}.tmp.validate.circ.support.test # test it whether it is match with circRNAs it should be no
awk '{if($10==2) {n11=split($11,t11,",");n12=split($12,t12,","); print $1"\t"($7+t11[1])"\t"($8-t11[2])"\t"$16"\t"$4"\t"$6} }' ${sample}.tmp.validate.circ.match > ${sample}.tmp.validate.circ.match.2  # transform the reads to splice site
awk '{a[$1"\t"$2"\t"$3"\t"$4]++}END{for(i in a){print i"\t"a[i] | "sort -k 1"}}' ${sample}.tmp.validate.circ.match.2 | awk '{m=substr($4,length($4)-1,1); print $0"\t"m }' > ${sample}.validate.circ.match #circRNA related splicing site and support paired reads number
#comapre the support with circRNA other end.
awk 'ARGIND==1{a[$1"\t"$2"\t"$3"\t"$4]=$0}ARGIND>1&&($1"\t"$2"\t"$3"\t"$4 in a ){print $0"\t"a[$1"\t"$2"\t"$3"\t"$4]}' ${sample}.validate.circ.match pair_proper.circ_splice.intro > pair_proper_tmp.circ_splice.paired_support #different
#compare with intron
$bedtools intersect -wa -wb -b $reference/gencode.v19.annotation.intron.bed -a ${sample}.validate.circ.match > ${sample}.validate.circ.match_intron
awk '{if($2==$8 && $3==$9) print }' ${sample}.validate.circ.match_intron > ${sample}.validate.circ.match_intron_match
awk 'ARGIND==1{a[$1"\t"$2"\t"$3"\t"$4]}ARGIND>1&&!($1"\t"$2"\t"$3"\t"$4 in a ){print $0}' ${sample}.validate.circ.match_intron_match ${sample}.validate.circ.match_intron > ${sample}.validate.circ.match_intron_mis
awk '!a[$4]++' ${sample}.validate.circ.match_intron_mis | wc -l
awk '{if($2==$8 || $3==$9) print }' pair_proper_tmp.circ_splice.in_intron_mis | awk '!a[$4]++' | wc -l
#not in intron region
awk 'ARGIND==1{a[$4]}ARGIND>1&&!($4 in a ){print $0}' ${sample}.validate.circ.match_intron ${sample}.validate.circ.match > ${sample}.validate.circ.match.in_out
$bedtools intersect -wa -wb -b $reference/gencode.v19.annotation.single.bed -a ${sample}.validate.circ.match.in_out > ${sample}.validate.circ.match.in_single
awk 'ARGIND==1{a[$4]}ARGIND>1&&!($4 in a ){print $0}' ${sample}.validate.circ.match.in_single ${sample}.validate.circ.match.in_out > ${sample}.validate.circ.match.in_intergenic
rm *tmp*
