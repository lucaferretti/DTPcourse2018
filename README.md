# Start of the session:

first open a terminal, then connect to the Pirbright cluster via SSH

> ssh DTCstudents@mallorn.pirbright.ac.uk

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

This will open you a shell on a computational node.

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

where each "." corresponds to a base that is identical to the reference.

Note how noisy the data are!

----------------------
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

these sequences belong to a single viral population. However, things are not so easy. Have a look at the frequencies of the variants. How are these frequencies distributed? Can you spot something that looks odd?
