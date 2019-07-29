## Performance comparison between NextPolish, Pilon and Racon using simulation data  

* **REQUIREMENT**
    * [ART v2.5.8](https://www.niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm)
    * [PBSIM v1.0.4](https://github.com/pfaucon/PBSIM-PacBio-Simulator)
    * [CANU v1.8](https://github.com/marbl/canu)
    * [Pilon v1.23](https://github.com/broadinstitute/pilon)
    * [Racon v1.3.3](https://github.com/isovic/racon)
    * [NextPolish v1.0.3](https://github.com/Nextomics/NextPolish)
    * [Quast v5.0.2](https://github.com/ablab/quast)
    
* **Download reference**   
```
curl -SL ftp://ftp.ensembl.org/pub/release-96/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.1.fa.gz | gunzip - > chr01.fa
```

* **Simulate PacBio data**  
```
pbsim --data-type CLR --model_qc /PBSIM-PacBio-Simulator/data/model_qc_clr --depth 50 --length-mean 10000 --accuracy-mean 0.85 --prefix pacbio chr01.fa
```

* **Simulate Illumina data**  
```
art_illumina -ss HS25 -i chr01.fa -p -l 150 -f 50 -m 300 -s 10 -o NGS_
```

* **Assemble reference**  
```
canu -pacbio-raw pacbio_0001.fastq -p asm -d canu-pb useGrid=True genomeSize=250m gridEngineMemoryOption="-l vf=MEMORY"
```

* **Run Pilon**  
    + work.sh
    ```
    genome=asm.contigs.fasta  
    reads1=NGS_1.fq    
    reads2=NGS_1.fq    
    input=${genome}    
    for i in {1..4};do   
        NextPolish/bin/bwa index ${input};  
        NextPolish/bin/bwa mem -t 25 ${input} ${reads1} ${reads2} |NextPolish/bin/samtools view  -b - |NextPolish/bin/samtools fixmate -m --threads 5 - - |NextPolish/bin/samtools sort -m 5g --threads 5 - -o ${input}.sort.bam;   
        NextPolish/bin/samtools index ${input}.sort.bam;  
        time -p java -Xmx50G -jar /home/huj/software/pilon-1.23.jar --genome ${input} --frags ${input}.sort.bam  --output ${genome}.pilon.v${i} --threads 5 --fix bases;  
        input=${genome}.pilon.v${i}.fasta;  
    done
    ```
    + Run  
    `nohup sh work.sh > pilon.log &`
    + CPU time used for polishing  
    `egrep 'user|sys' pilon.log|awk '{x+=$2}END{print x}'`

* **Run Racon**   
    + work.sh
    ```
    awk '{if (NR%4==1){print $0"1"}else{print $0}}' NGS_1.fq > NGS_1.rn.fq;  
    awk '{if (NR%4==1){print $0"1"}else{print $0}}' NGS_2.fq > NGS_2.rn.fq;  
    cat NGS_1.rn.fq NGS_2.rn.fq > NGS.rn.fq;  
    genome=asm.contigs.fasta     
    reads1=NGS_1.rn.fq  
    reads2=NGS_2.rn.fq  
    input=${genome}  
    for i in {1..4};do  
        NextPolish/bin/minimap2 -ax sr ${input} ${reads1} ${reads2} > input.sam
        time -p racon NGS.rn.fq input.sam ${input} --include-unpolished --threads 5 > ${genome}.racon.v${i}.fasta;  
        input=${genome}.racon.v${i}.fasta;  
    done
    ```
    + Run  
    `nohup sh work.sh > racon.log &`
    + CPU time used for polishing  
    `egrep 'user|sys' racon.log|awk '{x+=$2}END{print x}'`

* **Run NextPolish**  
    + run.cfg  
    ```
    [General]
    job_type = local
    job_prefix = nextPolish
    task = 1212
    rewrite = yes
    rerun = 3
    parallel_jobs = 1
    multithread_jobs = 5
    genome = asm.contigs.fasta
    genome_size = auto
    workdir = ./01_rundir
    polish_options = -p {multithread_jobs}

    [sgs_option]
    sgs_fofn = sgs.fofn
    sgs_options = -max_depth 100 -bwa
    ```
    + Run  
    ```
    ls NGS_1.fq NGS_2.fq > sgs.fofn;  
    nextPolish run.cfg;
    ```
    + CPU time used for polishing  
    ```
    egrep 'user|sys' 01_rundir/*/0*.polish.ref.sh.work/polish_genome*/nextPolish.sh.e|awk '{print $2}'|sed 's/m/\t/' |sed 's/s//' |awk '{x+=$1*60+$2}END{print x}'
    ```

* **Run Quast**  
    + Input    
        * Pilon x 1  
        `asm.contigs.pilonv1.fasta`
        * Pilon x 2  
        `asm.contigs.pilonv2.fasta`
        * Pilon x 3  
        `asm.contigs.pilonv3.fasta`
        * Pilon x 4  
        `asm.contigs.pilonv4.fasta`
        * Racon x 1  
        `asm.contigs.raconv1.fasta`
        * Racon x 2  
        `asm.contigs.raconv2.fasta`
        * Racon x 3  
        `asm.contigs.raconv3.fasta`
        * Racon x 4  
        `asm.contigs.raconv4.fasta`
        * NextPolish x 1  
        ```
        cat 01_rundir/01.kmer_count/*.polish.ref.sh.work/polish_genome*/genome.nextpolish.part*.fasta > asm.contigs.nextpolishv1.fasta
        ```
        * NextPolish x 2  
        ```
        cat 01_rundir/03.kmer_count/*mar.polish.ref.sh.work/polish_genome*/genome.nextpolish.part*.fasta > asm.contigs.nextpolishv2.fasta
        ```
    + Run    
    ```
    quast/quast-5.0.2/quast.py -e --min-contig 1000000 --min-alignment 50000 --extensive-mis-size 7000 -r chr01.fa asm.contigs.fasta asm.contigs.nextpolishv1.fasta asm.contigs.nextpolishv2.fasta asm.contigs.pilonv1.fasta asm.contigs.pilonv2.fasta asm.contigs.pilonv3.fasta asm.contigs.pilonv4.fasta asm.contigs.raconv1.fasta asm.contigs.raconv2.fasta asm.contigs.raconv3.fasta asm.contigs.raconv4.fasta
    ```
    + **Click [here](TEST1.pdf) for result**
    

