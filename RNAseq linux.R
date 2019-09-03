# unzip fastq.gz
gunzip *.gz

# under public/ ,there are resource/, work/ dirs.
# build genome/, genes/, bdgp6_tran/ dirs under resource/ dir.
# bdgp6_tran/ contains .ht2 files, which are index files, download from ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/data/bdgp6_tran.tar.gz
# genome/ contains .fa files, which are reference files, download from ftp://ftp.ensemblgenomes.org/pub/metazoa/release-44/fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.BDGP6.22.dna.toplevel.fa.gz
# genes/ contains .gtf file, which is annotation, download from ftp://ftp.ensemblgenomes.org/pub/metazoa/release-44/gtf/drosophila_melanogaster/Drosophila_melanogaster.BDGP6.22.44.gtf.gz

# There are data/, results/, fastqc/ dirs under the work/ dirs.
# data/ contains clean.fastq files.
# results/ contains prepDE.py and gtflist.
# fastqc/ contains fastqc results and flagstat results.

# fastqc, in data/ dir
fastqc *.fq -t 8 -o /public/work/fastqc/ &

# hisat2, in data/ dir
for i in *R1_paired.fq ; 
do
i=${i%R1_paired.fq*}; 
nohup hisat2 -p 4--dta -x /public/resource/bdgp6_tran/genome_tran -1 ${i}R1_paired.fq -2 ${i}R2_paired.fq -S  /public/work/results/${i}align.sam 2>/public/work/fastqc/${i}align.log  & 
done

# samtools, in results/ dir
for i in *.sam;
do 
i=${i%.align.sam*}; nohup samtools  view  -bS ${i}.align.sam > ${i}_unsorted.bam & 
done

# sam to sotrted.bam (in results/ dir) >1.3 version
for i in *.sam;
do 
i=${i%.align.sam*}; nohup samtools sort -@ 8 -o ${i}_sorted.bam ${i}.align.sam &
done

# Flagstat in results/ dir
for i in *.bam
do
i=${i%.bam*}; nohup samtools flagstat ${i}.bam > /public/work/fastqc/${i}.flagstat &
done

# bam Index, in results/ dir
for i in *.bam
do
i=${i%.bam*}; nohup samtools index ${i}.bam &
done

# igvtools (after building index), in results/ dir
for i in *.bam
do
i=${i%.bam*}; nohup igvtools count -z 5 -w 25  ${i}.bam ${i}.bam.tdf  /public/resource/genome/Drosophila_melanogaster.BDGP6.22.dna.toplevel.fa &
done

# stringtie, in results/ dir
# build gtf/ and tab/ dirs under results/ dir
for i in *.bam; 
do 
i=${i%.bam*}; nohup stringtie -p 8 -G /public/resource/gene/Drosophila_melanogaster.BDGP6.22.44.gtf -o /public/work/results/gtf/${i}.gtf -A /public/work/results/tab/${i}.tab -B -e -l ${i} ${i}.bam &
done

# under gtf/, build dir DESeq2/
# prepare gtflist.txt which contain sample gtf information and its path.
python /public/resource/prepDE.py -i gtflist.txt -g DESeq2/gene_count.csv -t DESeq2/transcript.csv