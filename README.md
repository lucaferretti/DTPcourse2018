# Start of the session:

first open a terminal, then connect to the Pirbright cluster via SSH

> ssh dtcstudents@mallorn.pirbright.ac.uk

and use the password that was given during the course. Then ask one of the demonstrators for the numeric code to enter.

----------------------

The first time that you connect, a few preparatory steps to work remotely:

choose a (unique!) username, containing only letters and numbers and underscores, then type

> tmux new -s MYUSERNAME

where MYUSERNAME is the username you chose. This will open a virtual session on the cluster that will avoid losing your work if you lose your connection!

Then create your own directory

> mkdir MYUSERNAME

copy all data there:

> cp Data/* MYUSERNAME/

then move there to work

> cd MYUSERNAME

Please do not work outside this directory!

Finally, reserve a computational unit for yourself by running

> srun -c 1 --tasks-per-node 1 --pty bash -i

> source /data/local/bash_profile_init

This will open you a shell on a computational node and load all defaults for this shell.

From now on, please avoid typing "exit" or anything similar. In case you want to get out, just close the window of the terminal.

Some tricks:

- to create another window: "CTRL-b c"

- to switch between windows: "CTRL-b 0" or "CTRL-b 1" or "CTRL-b n" or "CTRL-b p"

----------------------

After the first time, you can simply open a terminal again, reconnect via SSH as described above, then type

> tmux attach -t MYUSERNAME

and you'll get to the same point you were when you left.

----------------------
----------------------

# Alignment of short reads from viral deep sequencing:

Viral reads in the Fastq format can be found in the files "ViralReads_R1.fastq" and "ViralReads_R2.fastq".

These are two files because the sequencing was paired-end: the reads in the second file are paired with reads in the first.

A reference genome for a closely related FMDV strain can be found in Fasta format in "KR108954_1.fasta".

----------------------

The first step is read trimming, to remove adapter and low-quality bases at the ends of the reads:

> trim_galore --paired --fastqc -o . ViralReads_R1.fastq ViralReads_R2.fastq

and to make things easy, renaming the output files containing the trimmed reads:

> mv ViralReads_R1_val_1.fq ViralReads_trimmed_R1.fastq

> mv ViralReads_R2_val_2.fq ViralReads_trimmed_R2.fastq

----------------------

In preparation for the alignment with BWA, the reference should be indexed:

> bwa index KR108954_1.fasta

then the alignment can start:

> bwa mem KR108954_1.fasta ViralReads_R1.fastq ViralReads_R2.fastq > ViralReads.sam

The SAM output file now contains aligned reads, that should be converted to BAM format

> samtools view -b -q 20 -F 12 ViralReads.sam > ViralReads.bam

then sorted by position along the genome

> samtools sort ViralReads.bam -T tempname -O 'bam' > ViralReads.sorted.bam

and the final BAM file should be indexed for later purposes:

> samtools index ViralReads.sorted.bam

The content and position of the aligned read against the genome can be seen graphically from

> samtools tview ViralReads.sorted.bam KR108954_1.fasta

where each "." corresponds to a base that is identical to the reference. Press q or CTRL-c to exit.

Note how noisy the data are!

----------------------

# SNP calling from viral deep sequencing:

To call variants from these data, first extract the information base-by-base (pileup format) with samtools:

> samtools mpileup -f KR108954_1.fasta ViralReads.sorted.bam > ViralReads.pileup

then run Varscan2 on the pileup file to generate a tabular file with all the variants and related information

> java -jar /data/local/libexec/varscan/VarScan.jar pileup2snp ViralReads.pileup --min-coverage 100 --min-reads2 50 --min-avg-qual20 --min-var-freq 0.001 --p-value 0.05 > listVariants.tabular

You can also try to change the above parameters (see manual at http://dkoboldt.github.io/varscan/using-varscan.html#v2.3_pileup2snp) according to reasonable criteria, andd have have a look at the difference in the variants in the output file.

Then examine the resulting output file containing only variants that passed the filter, by doing

> less listVariants.tabular

The final question of this tutorial is:

these sequences belong to a single viral 'population', or viral swarm/quasispecies. In this context, we expect all sequences to form a 'cloud' of genotypes differing by a few random mutations. However, things are not so easy. Have a look at the frequencies of the variants. How are these frequencies distributed? Can you spot something that looks odd?

----------------------
----------------------

# Alignment and SNP calling from DNA sequencing of mosquitoes

List of files/datasets (Anopheles gambiae and stephensi):

- Anopheles-stephensi-Indian_SCAFFOLDS_AsteI2.fa : reference genome of A.stephensi

- Anopheles-gambiae-PEST_CHROMOSOMES_AgamP4.fa : reference genome of A.gambiae

- Anopheles-gambiae-PEST_BASEFEATURES_AgamP4.10.gtf : annotation of A.gambiae genome

- Anopheles-gambiae-PEST_BASEFEATURES_AgamP4.10.genes.bed : annotation of genes in BED format, generated by the command

> cat Anopheles-gambiae-PEST_BASEFEATURES_AgamP4.10.gtf | sed '1,5d' | grep -v transcript_id | cut -f 1,4,5,9 | cut -f 1 -d ";" | sed 's/gene_id //' | sed 's/"//g' > Anopheles-gambiae-PEST_BASEFEATURES_AgamP4.10.genes.bed

- SRR830370_exome_R1.fq : reads from exome of A.stephensi (first of paired-end pair)

- SRR830370_exome_R2.fq : reads from exome of A.stephensi (second of paired-end pair)

First thing to do is to reduce the data to a small subset, e.g. of a million lines, since it is more practical for a first test (every million lines, i.e. 250000 reads, will take about 4 minutes to align)

> head -n 1000000 SRR830370_exome_R1.fq > SRR830370_smallsubset_R1.fq

> head -n 1000000 SRR830370_exome_R2.fq > SRR830370_smallsubset_R2.fq

----------------------

Alignment of DNA reads from A.stephensis to different reference genomes:

Index reference genomes for various programs (takes about 10 minutes)

> samtools faidx Anopheles-gambiae-PEST_CHROMOSOMES_AgamP4.fa

> samtools faidx Anopheles-stephensi-Indian_SCAFFOLDS_AsteI2.fa

> bwa index Anopheles-gambiae-PEST_CHROMOSOMES_AgamP4.fa

> bwa index Anopheles-stephensi-Indian_SCAFFOLDS_AsteI2.fa

Then align A.stephensi reads both to A.gambiae and A.stephensi reference genomes with BWA

> bwa mem Anopheles-gambiae-PEST_CHROMOSOMES_AgamP4.fa SRR830370_smallsubset_R1.fq SRR830370_smallsubset_R2.fq | gzip -3 > SRR830370_smallsubset.bwa.Agambiae.sam.gz

> bwa mem Anopheles-stephensi-Indian_SCAFFOLDS_AsteI2.fa SRR830370_smallsubset_R1.fq SRR830370_smallsubset_R2.fq | gzip -3 > SRR830370_smallsubset.bwa.Astephensi.sam.gz

then transform into BAM

> zcat SRR830370_smallsubset.bwa.Agambiae.sam.gz | samtools view -bS - > SRR830370_smallsubset.bwa.Agambiae.bam

> zcat SRR830370_smallsubset.bwa.Astephensi.sam.gz | samtools view -bS - > SRR830370_smallsubset.bwa.Astephensi.bam

then sort

> samtools sort  -O 'bam' -T anytemp -o SRR830370_smallsubset.bwa.Agambiae.sorted.bam SRR830370_smallsubset.bwa.Agambiae.bam

> samtools sort  -O 'bam' -T anytemp -o SRR830370_smallsubset.bwa.Astephensi.sorted.bam SRR830370_smallsubset.bwa.Astephensi.bam

then index

> samtools index SRR830370_smallsubset.bwa.Agambiae.sorted.bam

> samtools index SRR830370_smallsubset.bwa.Astephensi.sorted.bam

and finally check the mapping statistics:

> samtools flagstat SRR830370_smallsubset.bwa.Agambiae.sorted.bam

> samtools flagstat SRR830370_smallsubset.bwa.Astephensi.sorted.bam

How do they compare? How many reads would you lose by mapping to the "wrong" reference genome? And how this depends on choosing to consider all mapped reads or properly mapped pairs only?

You can use

> samtools mpileup SRR830370_smallsubset.bwa.Astephensi.sorted.bam | less

or

> samtools tview SRR830370_smallsubset.bwa.Astephensi.sorted.bam Anopheles-stephensi-Indian_SCAFFOLDS_AsteI2.fa

to explore the data by hand. You can move around using arrows on the keyboard; press q to exit.

----------------------

SNP calling from the previous alignments using SAMtools:

> samtools mpileup -uf Anopheles-gambiae-PEST_CHROMOSOMES_AgamP4.fa SRR830370_smallsubset.bwa.Agambiae.sorted.bam | bcftools call -mO v > SRR830370_smallsubset.bwa.Agambiae.sorted.samtools.vcf

then filter for low quality polymorphisms and select either polymorphic variants

> bcftools filter -s LowQual -e '%QUAL<20 || DP<10 || AC>1 || AC<1' SRR830370_smallsubset.bwa.Agambiae.sorted.samtools.vcf | grep -v LowQual > SRR830370_smallsubset.GEM.Agambiae.sorted.samtools.filtered.variants.vcf

or select fixed differences with the reference

> bcftools filter -s LowQual -e '%QUAL<20 || DP<10 || AC<2' SRR830370_smallsubset.bwa.Agambiae.sorted.samtools.vcf | grep -v LowQual > SRR830370_smallsubset.GEM.Agambiae.sorted.samtools.filtered.diffs.vcf

To find the number of fixed differences per gene (including introns), you can use

> bedtools coverage -counts -a Anopheles-gambiae-PEST_BASEFEATURES_AgamP4.10.genes.bed -b SRR830370_smallsubset.GEM.Agambiae.sorted.samtools.filtered.diffs.vcf > genes.SRR830370_smallsubset.GEM.Agambiae.sorted.samtools.filtered.diffs.vcf.counts.txt

and similarly for the number of variants.

----------------------
----------------------

# Alignment of short reads and sex differences in gene expression from RNA-seq data of Anopheles gambiae

List of files/datasets (Anopheles gambiae):

- Anopheles-gambiae-PEST_CHROMOSOMES_AgamP4.fa : reference genome of A.gambiae

- Anopheles-gambiae-PEST_BASEFEATURES_AgamP4.10.gtf : annotation of A.gambiae genome

- SRR1805322_sub_R1.fastq : reads from transcriptome of adult male of A.gambiae (first of paired-end pair)

- SRR1805322_sub_R2.fastq : reads from transcriptome of adult male of A.gambiae (second of paired-end pair)

- SRR1805328_sub_R1.fastq : reads from transcriptome of adult female of A.gambiae (first of paired-end pair)

- SRR1805328_sub_R2.fastq : reads from transcriptome of adult female of A.gambiae (second of paired-end pair)

----------------------

Index reference transcriptome for GEM (takes about 15 minutes)

> gemtools index -i Anopheles-gambiae-PEST_CHROMOSOMES_AgamP4.fa

> gemtools t-index -a Anopheles-gambiae-PEST_BASEFEATURES_AgamP4.10.gtf -i Anopheles-gambiae-PEST_CHROMOSOMES_AgamP4.gem

Align with GEM (takes about 5 minutes per dataset)

> gemtools rna-pipeline -i Anopheles-gambiae-PEST_CHROMOSOMES_AgamP4.gem -a Anopheles-gambiae-PEST_BASEFEATURES_AgamP4.10.gtf -f SRR1805322_sub_R1.fastq SRR1805322_sub_R2.fastq -q 33 -t 2 --no-bam

then convert the resulting file to SAM

> zcat SRR1805322_sub_R1.map.gz | gem-2-sam -q offset-33 -o SRR1805322_sub_R1.sam -I Anopheles-gambiae-PEST_CHROMOSOMES_AgamP4.gem -l

and follow the steps above to generate a sorted, indexed bam file. Then, you could call SNPs from these data as well, if you would like to try.

Repeat the above steps for SRR1805328 as well.

----------------------

The files *.gtf.counts.txt contain the read counts for each gene therein.

You can use Excel/Openoffice or R to compare and analyse these data.

In case you would have problems with the ordering, you can use the "sort" command in Ubuntu before opening the files in Excel/Openoffice:

> sort -k 1 FILE1.gtf.counts.txt > FILE1.gtf.counts.sorted.txt

replacing FILE1 with the correct name.

You could also join the files for SRR1805322 and SRR1805328 in a single table before analysing them, using e.g. the join command:

> join -o 1.1,2.1,1.2,2.2 -a 1 -a 2 FILE1.gtf.counts.sorted.txt FILE2.gtf.counts.sorted.txt > jointcount.sorted.csv

replacing again FILE1 and FILE2 with the correct names.

This way you can compare the difference in expression in male vs female with the absolute expression of the gene.
Which genes are more expressed? Which genes have the largest sex differences in expression levels, compared with their average expression? Which genes are only expressed in males or females?

Finally, you can run a list of the most interesting genes that you found through easy-to-use webservers for functional annotation and gene enrichment analysis such as http://www.pantherdb.org or http://www.flymine.org in order to understand which biological functions are overrepresented among these genes.

----------------------

You could also join one of the previous files with the file genes.SRR830370_smallsubset.GEM.Agambiae.sorted.samtools.filtered.diffs.vcf.counts.txt as well.

If you do this, you could analyse in Excel/Openoffice or R the resulting table to answer an interesting question:

is there a correlation between the absolute expression level of a gene and its evolutionary properties?

More specifically, do highly expressed genes tend to diverge more between different species, or less?

----------------------

If you are interested in visualizing the data as expression profiles on the genome:

Bam files can be converted to traces in bedGraph format via commands like

> bedtools genomecov -trackline -trackopts autoscale=on -trackopts graphType=bar -bg -ibam SRR1805322_sub_1.sorted.bam > SRR1805322.bg
but the size of the file should be less than 20 Megabytes, e.g.

> head -n 500000 SRR1805322.bg > SRR1805322_partial.bg

then they can be visualised at https://www.vectorbase.org/Anopheles_gambiae/Location/View?r=2L%3A2686000-2746000 by clicking on Custom Tracks on the left and loading your bedGraph file.





