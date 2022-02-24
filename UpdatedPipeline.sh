#aligning genomes in hisat 
idx=mm10
spliceFile=gencode.vM24.splice.txt
processors=3

#this is the one that I used without rna-strandness
for i in *PE1.fastq.gz
do
name=$(echo $i | awk -F"_PE." '{print $1}')
echo $name
hisat2 -p ${processors} -x ${idx} --known-splicesite-infile ${spliceFile} -1 ${i} -2 ${name}_PE2.fastq.gz | samtools view -bS - | samtools sort -n - -o $name.sorted.bam
done

#remove the junk that doesnt align to the chromosomes
#mice only have 19 chromosomes and the 2 sex chromosomes vs humans have 23 chromosomes and the two sex chromosomes 
for i in *sorted.bam
do
    name=$(echo $i | awk -F".sorted." '{print $1}')
    echo $name
    samtools sort $i -o $i
    samtools index $i
    samtools view -bh $i chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chrX chrY chrM > $name.noJunk.bam
done

#remove PCR duplicates
for i in *.noJunk.bam
do
    name=$(echo $i | awk -F".noJunk." '{print $1}')
    echo $name
    samtools sort -n $i -o $name.noJunk.bam
    samtools fixmate -m $name.noJunk.bam $name.fixmate.bam
    samtools sort $name.fixmate.bam -o $name.fixmate.bam
    samtools markdup -rs $name.fixmate.bam $name.noDups.bam
    rm $name.fixmate.bam
done

#sort reads into genes for unstranded samples 
6

#the following three chunks are some troubleshooting code for when there was a corrupted bam 
for i in *sorted2.bam
do
    name=$(echo $i | awk -F".sorted2." '{print $1}')
    echo $name
    samtools sort $i -o $i
    samtools index $i
    samtools view -bh $i chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chrX chrY chrM > $name.noJunk2.bam
done

for i in *.noJunk2.bam
do
    name=$(echo $i | awk -F".noJunk2." '{print $1}')
    echo $name
    samtools sort -n $i -o $name.noJunk.bam
    samtools fixmate -m $name.noJunk.bam $name.fixmate.bam
    samtools sort $name.fixmate.bam -o $name.fixmate.bam
    samtools markdup -rs $name.fixmate.bam $name.noDups2.bam
    rm $name.fixmate.bam
done

for i in *.noDups3.bam
do
    name=$(echo $i | awk -F".noDups3." '{print $1}')
    echo $name
    samtools view -bf 1 $i | htseq-count -r pos -f bam --stranded=no - gencode.vM24.annotation.gtf > $name.gene.counts.txt
done

#this code is to get some quality control measures

echo -e "name\tall.counts\twithout.junk\tunique.counts" >> qc_metrics.txt

for bam in *sorted.bam
do
    name=$(echo $bam | awk -F".sorted.bam" '{print $1}')
    echo $name
    ALL_COUNTS=`samtools view -c $name.sorted.bam`
    WITHOUT_JUNK=`samtools view -c $name.noJunk.bam`
    UNIQUE_COUNTS=`samtools view -c $name.noDups.bam`
    echo -e "${name}\t$ALL_COUNTS\t$WITHOUT_JUNK\t$UNIQUE_COUNTS" >> qc_metrics.txt
done
