#Set input and parameters
round=2
threads=20 
read1=reads1_R1.fastq.gz
read2=reads1_R1.fastq.gz
input=input.genome.fa
for ((i=1; i<=${round};i++)); do
#step 1:
	#index the genome file and do alignment
	bwa index ${input};
	bwa mem -t ${threads} ${input} ${read1} ${read2}|samtools view --threads 3 -F 0x4 -b - |samtools sort - -m 2g --threads 5 -o sgs.sort.bam;
	#index bam and genome files
	samtools index -@ ${threads} sgs.sort.bam;
	samtools faidx ${input};
	#polish genome file
	python NextPolish/lib/nextPolish1.py -g ${input} -t 1 -p ${threads} -s sgs.sort.bam > genome.polishtemp.fa;
	input=genome.polishtemp.fa;
#step2:
	#index genome file and do alignment
	bwa index ${input};
	bwa mem -t ${threads} ${input} ${read1} ${read2}|samtools view --threads 3 -F 0x4 -b - |samtools sort - -m 2g --threads 5 -o sgs.sort.bam;
	#index bam and genome files
	samtools index -@ ${threads} sgs.sort.bam;
	samtools faidx ${input};
	#polish genome file
	python NextPolish/lib/nextPolish1.py -g ${input} -t 2 -p ${threads} -s sgs.sort.bam > genome.nextpolish.fa;
	input=genome.nextpolish.fa;
done;
#Finally polished genome file: genome.nextpolish.fa
