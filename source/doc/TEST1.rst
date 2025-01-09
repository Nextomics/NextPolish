.. _simulate_short_reads:

.. title:: Simulated short reads

Performance comparison between NextPolish, Pilon and Racon using simulated short reads
--------------------------------------------------------------------------------------

**REQUIREMENT**

   -  `ART v2.5.8 <https://www.niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm>`__
   -  `PBSIM v1.0.4 <https://github.com/pfaucon/PBSIM-PacBio-Simulator>`__
   -  `CANU v1.8 <https://github.com/marbl/canu>`__
   -  `Pilon v1.23 <https://github.com/broadinstitute/pilon>`__
   -  `Racon v1.3.3 <https://github.com/isovic/racon>`__
   -  `NextPolish v1.0.3 <https://github.com/Nextomics/NextPolish>`__
   -  `Quast v5.0.2 <https://github.com/ablab/quast>`__

1. **Download reference**
   
  .. code-block:: shell

    curl -SL ftp://ftp.ensembl.org/pub/release-96/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.1.fa.gz | gunzip - > chr01.fa

2. **Simulate PacBio data**

  .. code-block:: shell

    pbsim --data-type CLR --model_qc /PBSIM-PacBio-Simulator/data/model_qc_clr --depth 50 --length-mean 10000 --accuracy-mean 0.85 --prefix pacbio chr01.fa

3. **Simulate Illumina data**

  .. code-block:: shell

    art_illumina -ss HS25 -i chr01.fa -p -l 150 -f 50 -m 300 -s 10 -o NGS_

4. **Assemble reference**

  .. code-block:: shell

    canu -pacbio-raw pacbio_0001.fastq -p asm -d canu-pb useGrid=True genomeSize=250m gridEngineMemoryOption="-l vf=MEMORY"

5. **Run Pilon**

  - work.sh

  .. code-block:: shell

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

  - Run
    
  .. code-block:: shell

    nohup sh work.sh > pilon.log &

  - CPU time used for polishing
  
  .. code-block:: shell

    egrep 'user|sys' pilon.log|awk '{x+=$2}END{print x}'

6. **Run Racon**

  - work.sh

  .. code-block:: shell

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

  - Run
    
  .. code-block:: shell

    nohup sh work.sh > racon.log &

  - CPU time used for polishing

  .. code-block:: shell

    egrep 'user|sys' racon.log|awk '{x+=$2}END{print x}'

7. **Run NextPolish**

   -  run.cfg

  .. code-block:: shell

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

  - Run  
  
  .. code-block:: shell

   ls NGS_1.fq NGS_2.fq > sgs.fofn
   nextPolish run.cfg

  - CPU time used for polishing  

  .. code-block:: shell
    
    egrep 'user|sys' 01_rundir/*/0*.polish.ref.sh.work/polish_genome*/nextPolish.sh.e|awk '{print $2}'|sed 's/m/\t/' |sed 's/s//' |awk '{x+=$1*60+$2}END{print x}'

8. **Run Quast**

  -  Input

    -  Pilon x 1: ``asm.contigs.pilonv1.fasta``
    -  Pilon x 2: ``asm.contigs.pilonv2.fasta``
    -  Pilon x 3: ``asm.contigs.pilonv3.fasta``
    -  Pilon x 4: ``asm.contigs.pilonv4.fasta``
    -  Racon x 1: ``asm.contigs.raconv1.fasta``
    -  Racon x 2: ``asm.contigs.raconv2.fasta``
    -  Racon x 3: ``asm.contigs.raconv3.fasta``
    -  Racon x 4: ``asm.contigs.raconv4.fasta``
    -  NextPolish x 1::
        
        cat 01_rundir/01.kmer_count/*.polish.ref.sh.work/polish_genome*/genome.nextpolish.part*.fasta > asm.contigs.nextpolishv1.fasta

    -  NextPolish x 2::

        cat 01_rundir/03.kmer_count/*mar.polish.ref.sh.work/polish_genome*/genome.nextpolish.part*.fasta > asm.contigs.nextpolishv2.fasta

  -  Run

  .. code-block:: shell

    quast/quast-5.0.2/quast.py -e --min-contig 1000000 --min-alignment 50000 --extensive-mis-size 7000 -r chr01.fa asm.contigs.fasta asm.contigs.nextpolishv1.fasta asm.contigs.nextpolishv2.fasta asm.contigs.pilonv1.fasta asm.contigs.pilonv2.fasta asm.contigs.pilonv3.fasta asm.contigs.pilonv4.fasta asm.contigs.raconv1.fasta asm.contigs.raconv2.fasta asm.contigs.raconv3.fasta asm.contigs.raconv4.fasta

  .. object:: Quast result

  +------------------------+-----------+------------------------+------------------------+-----------------+---------------------+-------------------+------------------+--------------------+-------------------+-------------------+-------------------+
  |                        |asm.contigs|asm.contigs.nextpolishv1|asm.contigs.nextpolishv2|asm.contigs.pilonv1|asm.contigs.pilonv2|asm.contigs.pilonv3|asm.contigs.pilonv4|asm.contigs.raconv1|asm.contigs.raconv2|asm.contigs.raconv3|asm.contigs.raconv4|
  +========================+===========+========================+========================+===================+===================+===================+==================+====================+===================+===================+===================+
  |Total length (>= 0 bp)  |224780032  |224716364               |215224152               |215223160          |215223131          |215223143          |215223109         |215217457           |215212057          |215209603          |215208478          |
  +------------------------+-----------+------------------------+------------------------+-------------------+-------------------+-------------------+------------------+--------------------+-------------------+-------------------+-------------------+
  |Reference length        |248956422  |248956422               |248956422               |248956422          |248956422          |248956422          |248956422         |248956422           |248956422          |248956422          |248956422          |
  +------------------------+-----------+------------------------+------------------------+-------------------+-------------------+-------------------+------------------+--------------------+-------------------+-------------------+-------------------+
  |Unaligned length        |56553      |61272                   |61269                   |62646              |61699              |61703              |61275             |163683              |177973             |176917             |193791             |
  +------------------------+-----------+------------------------+------------------------+-------------------+-------------------+-------------------+------------------+--------------------+-------------------+-------------------+-------------------+
  |Largest alignment       |38684842   |38657142                |38657130                |38657017           |38656999           |38657033           |38657014          |38554506            |38537001           |38535515           |38523009           |
  +------------------------+-----------+------------------------+------------------------+-------------------+-------------------+-------------------+------------------+--------------------+-------------------+-------------------+-------------------+
  |# mismatches per 100 kbp|17.82      |2.38                    |2.26                    |2.92               |2.39               |2.31               |2.31              |3.08                |2.91               |2.87               |2.64               |
  +------------------------+-----------+------------------------+------------------------+-------------------+-------------------+-------------------+------------------+--------------------+-------------------+-------------------+-------------------+
  |# indels per 100 kbp    |121.45     |0.81                    |0.71                    |1.60               |1.36               |1.28               |1.25              |1.97                |1.16               |1.08               |1.00               |
  +------------------------+-----------+------------------------+------------------------+-------------------+-------------------+-------------------+------------------+--------------------+-------------------+-------------------+-------------------+
  |# mismatches            |38286      |5107                    |4863                    |6275               |5134               |4974               |4957              |6625                |6249               |6177               |5684               |
  +------------------------+-----------+------------------------+------------------------+-------------------+-------------------+-------------------+------------------+--------------------+-------------------+-------------------+-------------------+
  |# indels                |261011     |1736                    |1527                    |3447               |2917               |2754               |2684              |4242                |2494               |2312               |2148               |
  +------------------------+-----------+------------------------+------------------------+-------------------+-------------------+-------------------+------------------+--------------------+-------------------+-------------------+-------------------+


  .. note:: The complete result of Quast can be seen from :download:`here <./TEST1.pdf>`.


