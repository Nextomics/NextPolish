.. _tutorial:

Tutorial
~~~~~~~~

.. contents::
    :local:

Polishing using short reads only
--------------------------------

#. Prepare sgs_fofn

   .. code-block:: shell

      ls reads1_R1.fq reads1_R2.fq reads2_R1.fq reads2_R2.fq > sgs.fofn

#. Create run.cfg

   .. code-block:: shell

    [General]
    job_type = local
    job_prefix = nextPolish
    task = best
    rewrite = yes
    rerun = 3
    parallel_jobs = 6
    multithread_jobs = 5
    genome = ./raw.genome.fasta #genome file
    genome_size = auto
    workdir = ./01_rundir
    polish_options = -p {multithread_jobs}

    [sgs_option]
    sgs_fofn = ./sgs.fofn
    sgs_options = -max_depth 100 -bwa

#. Run

   .. code-block:: shell

      nextPolish run.cfg

#. Finally polished genome

   - Sequence: ``/path_to_work_directory/genome.nextpolish.fasta``
   - Statistics: ``/path_to_work_directory/genome.nextpolish.fasta.stat``

.. tip:: User defined alignment pipeline, which will be faster than the default pipeline when runing on a local system. The accuracy of the polished genome is the same as the default.
    
    .. code-block:: shell

       #Set input and parameters
       round=2
       threads=20
       read1=reads_R1.fastq.gz
       read2=reads_R2.fastq.gz
       input=input.genome.fa
       for ((i=1; i<=${round};i++)); do
       #step 1:
          #index the genome file and do alignment
          bwa index ${input};
          bwa mem -t ${threads} ${input} ${read1} ${read2}|samtools view --threads 3 -F 0x4 -b -|samtools fixmate -m --threads 3  - -|samtools sort -m 2g --threads 5 -|samtools markdup --threads 5 -r - sgs.sort.bam
          #index bam and genome files
          samtools index -@ ${threads} sgs.sort.bam;
          samtools faidx ${input};
          #polish genome file
          python NextPolish/lib/nextpolish1.py -g ${input} -t 1 -p ${threads} -s sgs.sort.bam > genome.polishtemp.fa;
          input=genome.polishtemp.fa;
       #step2:
          #index genome file and do alignment
          bwa index ${input};
          bwa mem -t ${threads} ${input} ${read1} ${read2}|samtools view --threads 3 -F 0x4 -b -|samtools fixmate -m --threads 3  - -|samtools sort -m 2g --threads 5 -|samtools markdup --threads 5 -r - sgs.sort.bam
          #index bam and genome files
          samtools index -@ ${threads} sgs.sort.bam;
          samtools faidx ${input};
          #polish genome file
          python NextPolish/lib/nextpolish1.py -g ${input} -t 2 -p ${threads} -s sgs.sort.bam > genome.nextpolish.fa;
          input=genome.nextpolish.fa;
       done;
       #Finally polished genome file: genome.nextpolish.fa

.. _long_read_polish:

Polishing using long reads only
-------------------------------------

#. Prepare lgs_fofn

   .. code-block:: shell

      ls reads1.fq reads2.fa.gz > lgs.fofn

#. Create run.cfg

   .. code-block:: shell

    [General]
    job_type = local
    job_prefix = nextPolish
    task = best
    rewrite = yes
    rerun = 3
    parallel_jobs = 6
    multithread_jobs = 5
    genome = ./raw.genome.fasta #genome file
    genome_size = auto
    workdir = ./01_rundir
    polish_options = -p {multithread_jobs}

    [lgs_option]
    lgs_fofn = ./lgs.fofn
    lgs_options = -min_read_len 1k -max_depth 100
    lgs_minimap2_options = -x map-ont

#. Run

   .. code-block:: shell

      nextPolish run.cfg

#. Finally polished genome

   - Sequence: ``/path_to_work_directory/genome.nextpolish.fasta``
   - Statistics: ``/path_to_work_directory/genome.nextpolish.fasta.stat``

.. tip:: User defined alignment pipeline, which will be faster than the default pipeline when runing on a local system. The accuracy of the polished genome is the same as the default.
    
    .. code-block:: shell

        #Set input and parameters
        round=2
        threads=20
        read=read.fasta.gz
        read_type=ont #{clr,hifi,ont}, clr=PacBio continuous long read, hifi=PacBio highly accurate long reads, ont=NanoPore 1D reads
        declare -A mapping_option=(["clr"]="map-pb" ["hifi"]="asm20" ["ont"]="map-ont")
        input=input.genome.fa

        for ((i=1; i<=${round};i++)); do
            minimap2 -ax ${mapping_option[$read_type]} -t ${threads} ${input} ${read}|samtools sort - -m 2g --threads ${threads} -o lgs.sort.bam;
            samtools index lgs.sort.bam;
            ls `pwd`/lgs.sort.bam > lgs.sort.bam.fofn;
            python NextPolish/lib/nextpolish2.py -g ${input} -l lgs.sort.bam.fofn -r ${read_type} -p ${threads} -sp -o genome.nextpolish.fa;
            if ((i!=${round}));then
                mv genome.nextpolish.fa genome.nextpolishtmp.fa;
                input=genome.nextpolishtmp.fa;
            fi;
        done;
        # Finally polished genome file: genome.nextpolish.fa

Polishing using short reads and long reads
------------------------------------------------

#. Prepare sgs_fofn

   .. code-block:: shell

      ls reads1_R1.fq reads1_R2.fq reads2_R1.fq reads2_R2.fq > sgs.fofn

#. Prepare lgs_fofn

   .. code-block:: shell

      ls reads1.fq reads2.fa.gz > lgs.fofn

#. Create run.cfg

   .. code-block:: shell

    [General]
    job_type = local
    job_prefix = nextPolish
    task = best
    rewrite = yes
    rerun = 3
    parallel_jobs = 6
    multithread_jobs = 5
    genome = ./raw.genome.fasta
    genome_size = auto
    workdir = ./01_rundir
    polish_options = -p {multithread_jobs}

    [sgs_option]
    sgs_fofn = ./sgs.fofn
    sgs_options = -max_depth 100 -bwa

    [lgs_option]
    lgs_fofn = ./lgs.fofn
    lgs_options = -min_read_len 1k -max_depth 100
    lgs_minimap2_options = -x map-ont

#. Run

   .. code-block:: shell

      nextPolish run.cfg

#. Finally polished genome

   - Sequence: ``/path_to_work_directory/genome.nextpolish.fasta``
   - Statistics: ``/path_to_work_directory/genome.nextpolish.fasta.stat``

Polishing using short reads and hifi reads
------------------------------------------------

#. Prepare sgs_fofn

   .. code-block:: shell

      ls reads1_R1.fq reads1_R2.fq reads2_R1.fq reads2_R2.fq > sgs.fofn

#. Prepare hifi_fofn

   .. code-block:: shell

      ls reads1.fq reads2.fa.gz > hifi.fofn

#. Create run.cfg

   .. code-block:: shell

    [General]
    job_type = local
    job_prefix = nextPolish
    task = best
    rewrite = yes
    rerun = 3
    parallel_jobs = 6
    multithread_jobs = 5
    genome = ./raw.genome.fasta
    genome_size = auto
    workdir = ./01_rundir
    polish_options = -p {multithread_jobs}

    [sgs_option]
    sgs_fofn = ./sgs.fofn
    sgs_options = -max_depth 100 -bwa

    [hifi_option]
    hifi_fofn = ./hifi.fofn
    hifi_options = -min_read_len 1k -max_depth 100
    hifi_minimap2_options = -x map-pb

#. Run

   .. code-block:: shell

      nextPolish run.cfg

#. Finally polished genome

   - Sequence: ``/path_to_work_directory/genome.nextpolish.fasta``
   - Statistics: ``/path_to_work_directory/genome.nextpolish.fasta.stat``
