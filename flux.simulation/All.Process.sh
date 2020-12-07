###Step 0
###Please make sure the following programms are added in $PATH.
###kallisto salmon rsem-prepare-reference rsem-calculate-expression STAR cufflinks

test -d hg19 || mkdir hg19
test -d simulation.run1 || mkdir simulation.run1
test -d simulation.run2 || mkdir simulation.run2
test -d simulation.run3 || mkdir simulation.run3
test -d simulation.run4 || mkdir simulation.run4

cd ./hg19

###Step 1
wget -c ftp://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz  &&  gunzip Homo_sapiens.GRCh37.75.gtf.gz
wget -c ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh37.75.cdna.all.fa.gz  &&  gunzip Homo_sapiens.GRCh37.75.cdna.all.fa.gz
wget -c ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.75.dna_rm.primary_assembly.fa.gz  &&  gunzip Homo_sapiens.GRCh37.75.dna_rm.primary_assembly.fa.gz

mkdir chr_fa
faSplit byName Homo_sapiens.GRCh37.75.dna_rm.primary_assembly.fa chr_fa/


###Step 2
awk '$2=="protein_coding" && ($3=="exon" || $3=="start_codon" || $3=="stop_codon")' Homo_sapiens.GRCh37.75.gtf | awk '(length($1) <= 2 && $1!="MT" ) || $1~/^GL/' > Homo_sapiens.GRCh37.75.sorted.gtf
awk 'BEGIN{FS=OFS="\t"}{gsub(/\"/,"\t");print $10,$12}' Homo_sapiens.GRCh37.75.sorted.gtf | uniq > Homo_sapiens.GRCh37.75.sorted.names
awk 'FILENAME==ARGV[1]{array[">"$2]=1}FILENAME==ARGV[2]{if($1~/^>/){if($1 in array){a=1}else{a=0}};if(a==1){print $1}}' Homo_sapiens.GRCh37.75.sorted.names Homo_sapiens.GRCh37.75.cdna.all.fa > Homo_sapiens.GRCh37.75.sorted.transcript.fa


###Step 3
kallisto index  Homo_sapiens.GRCh37.75.sorted.transcript.fa -i Homo_sapiens.GRCh37.75.sorted.transcript.kallisto.index

salmon index --transcripts Homo_sapiens.GRCh37.75.sorted.transcript.fa --index Homo_sapiens.GRCh37.75.sorted.transcript.salmon.index

mkdir rsem.star
rsem-prepare-reference --gtf Homo_sapiens.GRCh37.75.sorted.gtf --star -p 8 Homo_sapiens.GRCh37.75.dna_rm.primary_assembly.fa rsem.star/rsem.star

mkdir star-genome
STAR --runThreadN 12 --runMode genomeGenerate --genomeDir star-genome --genomeFastaFiles Homo_sapiens.GRCh37.75.dna_rm.primary_assembly.fa  --sjdbGTFfile Homo_sapiens.GRCh37.75.sorted.gtf  --sjdbOverhang 99


###Step 4
export FLUX_MEM="12G"

for i in `seq 1 3`;do
cd ../simulation.run${i}

###4.1 Simulation
flux-simulator -p hg19.ver75.hydrolysis.par

awk 'NR%8 >=1 && NR%8 <=4' hg19.ver75.hydrolysis.fastq > hg19.ver75.hydrolysis.1.fastq
awk 'NR%8 ==0 || NR%8 >=5' hg19.ver75.hydrolysis.fastq > hg19.ver75.hydrolysis.2.fastq

###4.2 quantification
kallisto quant -i ../hg19/Homo_sapiens.GRCh37.75.sorted.transcript.kallisto.index -t 4 hg19.ver75.hydrolysis.1.fastq hg19.ver75.hydrolysis.2.fastq  -o kallisto.output

salmon quant -i ../hg19/Homo_sapiens.GRCh37.75.sorted.transcript.salmon.index -l IU -1 hg19.ver75.hydrolysis.1.fastq -2 hg19.ver75.hydrolysis.2.fastq -o salmon.output -p 4

STAR --runThreadN 12 --runMode alignReads  --genomeDir ../hg19/star-genome  --readFilesIn hg19.ver75.hydrolysis.1.fastq hg19.ver75.hydrolysis.2.fastq --sjdbOverhang 99 --outFileNamePrefix star.output --outSAMtype BAM SortedByCoordinate --outFilterMultimapNmax 100 --quantMode TranscriptomeSAM  --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical

cufflinks -p 8 -u -G ../hg19/Homo_sapiens.GRCh37.75.sorted.gtf -o cufflinks star.outputAligned.sortedByCoord.out.bam

rsem-calculate-expression --paired-end -p 8 --star hg19.ver75.hydrolysis.1.fastq hg19.ver75.hydrolysis.2.fastq ../hg19/rsem.star/rsem.star rsem.star

###4.3 comparison
awk 'BEGIN{FS=OFS="\t"}FILENAME==ARGV[1]{array0[$2]=$5 * 1000000}FILENAME==ARGV[2]{array1[$1]=$4}FILENAME==ARGV[3]{array2[$1]=$5}FILENAME==ARGV[4]{array3[$1]=$6}FILENAME==ARGV[5]{array4[$1]=$9}FILENAME==ARGV[6]{print $0,array0[$2],array1[$2],array2[$2],array3[$2],array4[$2]}' hg19.ver75.hydrolysis.pro  salmon.output/quant.sf kallisto.output/abundance.tsv rsem.star.isoforms.results cufflinks/isoforms.fpkm_tracking ../hg19/Homo_sapiens.GRCh37.75.sorted.names  > simulation.run${i}.compare.result
done


for i in `seq 4 4`;do
cd ../simulation.run${i}

###4.1 Simulation
~/program/flux-simulator-1.2.1/bin/flux-simulator -p hg19.ver75.hydrolysis.par

###4.2 quantification
kallisto quant --single -l 100 -s 2 -i ../hg19/Homo_sapiens.GRCh37.75.sorted.transcript.kallisto.index -t 4 --single hg19.ver75.hydrolysis.fastq -o kallisto.output

salmon quant -i ../hg19/Homo_sapiens.GRCh37.75.sorted.transcript.salmon.index -l IU -r hg19.ver75.hydrolysis.fastq  -o salmon.output -p 4

STAR --runThreadN 12 --runMode alignReads  --genomeDir ../hg19/star-genome  --readFilesIn hg19.ver75.hydrolysis.fastq --sjdbOverhang 99 --outFileNamePrefix star.output --outSAMtype BAM SortedByCoordinate --outFilterMultimapNmax 100 --quantMode TranscriptomeSAM  --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical

cufflinks -p 8 -u -G ../hg19/Homo_sapiens.GRCh37.75.sorted.gtf -o cufflinks star.outputAligned.sortedByCoord.out.bam

rsem-calculate-expression  -p 8  --star  hg19.ver75.hydrolysis.fastq  ../hg19/rsem.star/rsem.star rsem.star

###4.3 comparison
awk 'BEGIN{FS=OFS="\t"}FILENAME==ARGV[1]{array0[$2]=$5 * 1000000}FILENAME==ARGV[2]{array1[$1]=$4}FILENAME==ARGV[3]{array2[$1]=$5}FILENAME==ARGV[4]{array3[$1]=$6}FILENAME==ARGV[5]{array4[$1]=$9}FILENAME==ARGV[6]{print $0,array0[$2],array1[$2],array2[$2],array3[$2],array4[$2]}' hg19.ver75.hydrolysis.pro  salmon.output/quant.sf kallisto.output/abundance.tsv rsem.star.isoforms.results cufflinks/isoforms.fpkm_tracking ../hg19/Homo_sapiens.GRCh37.75.sorted.names  > simulation.run${i}.compare.result
done


###Step 5 Summarize all simulation results
cd ..
R --vanilla < Plot.Rscript

