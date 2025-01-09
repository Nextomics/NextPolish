.. _actual_short_reads:

.. title:: Actual biological data

Performance comparison between NextPolish and Pilon using actual biological data
--------------------------------------------------------------------------------

**REQUIREMENT**

  -  `Miniasm v0.2 <https://github.com/lh3/miniasm>`__
  -  `Falcon v1.8.7 <https://github.com/PacificBiosciences/FALCON>`__
  -  `Pilon v1.23 <https://github.com/broadinstitute/pilon>`__
  -  `Racon v1.3.3 <https://github.com/isovic/racon>`__
  -  `NextPolish v1.0.3 <https://github.com/Nextomics/NextPolish>`__
  -  `Seqkit v0.10.1 <https://github.com/shenwei356/seqkit>`__
  -  `Gmap v2017-01-14 <http://research-pub.gene.com/gmap/>`__
  -  `Freebayes v1.2.0-10 <https://github.com/ekg/freebayes>`__

1. **Download data**

  -  Sequencing data

    -  `Arabidopsis thaliana <https://www.nature.com/articles/s41467-018-03016-2>`__
    -  `Homo sapiens <https://www.nature.com/articles/ncomms12065>`__

  -  Genes

    -  `Arabidopsis thaliana <https://www.arabidopsis.org/download/index-auto.jsp?dir=%2Fdownload_files%2FGenes%2FTAIR10_genome_release%2FTAIR10_gff3>`__
    -  `Homo sapiens <https://asia.ensembl.org/Homo_sapiens/Info/Index>`__

2. **Assembly**

  - *Arabidopsis thaliana*

    - PacBio data

    .. code-block:: shell

     minimap2 -x ava-pb pb.reads.fq pb.reads.fq | gzip -1 > overlaps.paf.gz
     miniasm -f pb.reads.fq overlaps.paf.gz > miniasm.gfa
     awk '{if($1=="S"){print ">"$2;print $3}}' miniasm.gfa > miniasm.fasta

    - NanoPore data

    .. code-block:: shell

     minimap2 -x ava-ont ont.reads.fq ont.reads.fq | gzip -1 > overlaps.paf.gz
     miniasm -f ont.reads.fq overlaps.paf.gz > miniasm.gfa
     awk '{if($1=="S"){print ">"$2;print $3}}' miniasm.gfa > miniasm.fasta

  - *Homo sapiens*

    - fc_run.cfg
      
    .. code-block:: shell

        job_type = sge
        input_fofn = input.fofn
        input_type = raw

        length_cutoff = 11000
        length_cutoff_pr = 12000

        stop_all_jobs_on_failure = False
        target = assembly

        job_queue = all.q
        sge_option_da = -pe smp 4 -q %(job_queue)s
        sge_option_la = -pe smp 4 -q %(job_queue)s
        sge_option_pda = -pe smp 4 -q %(job_queue)s
        sge_option_pla = -pe smp 4 -q %(job_queue)s
        sge_option_fc = -pe smp 10 -q %(job_queue)s
        sge_option_cns = -pe smp 4 -q %(job_queue)s

        pa_concurrent_jobs = 499
        ovlp_concurrent_jobs = 499
        cns_concurrent_jobs = 499

        pa_HPCdaligner_option =  -v -B256 -t12 -w8 -e0.75 -k18 -h260 -l2000 -s1000 -T4
        ovlp_HPCdaligner_option = -v -B128 -t12 -k20 -h360 -e.96 -l1800 -s1000 -T4

        pa_DBsplit_option = -x1000 -s200 -a
        ovlp_DBsplit_option = -x1000 -s200

        falcon_sense_option = --output_multi --min_idt 0.75 --min_cov 4  --max_n_read 200 --n_core 4
        overlap_filtering_setting = --max_diff 70 --max_cov 100  --min_cov 2 --bestn 10 --n_core 10

    - Run
      
    .. code-block:: shell
      
      fc_run.py fc_run.cfg

3. **Run Racon**

  -  *Arabidopsis thaliana*

    - PacBio data

    .. code-block:: shell

     genome=miniasm.fasta  
     reads=pb.reads.fq        
     input=${genome}  
     for i in {1..4};do
         minimap2 -x map-pb ${input} ${reads} > align.paf;
         racon -t 10 ${reads} align.paf ${input} > ${genome}.racon.v${i}.fasta;
         input=${genome}.racon.v${i}.fasta;
     done;

    - NanoPore data

    .. code-block:: shell

     genome=miniasm.fasta  
     reads=ont.reads.fq        
     input=${genome}  
     for i in {1..4};do
         minimap2 -x map-ont ${input} ${reads} > align.paf;
         racon -t 10 ${reads} align.paf ${input} > ${genome}.racon.v${i}.fasta;
         input=${genome}.racon.v${i}.fasta;
     done;

4. **Run Pilon**

  -  *Arabidopsis thaliana*

    -  work.sh
       
    .. code-block:: shell

      genome=miniasm.racon.v4.fasta
      reads1=NGS_1.fq    
      reads2=NGS_1.fq    
      input=${genome}    
      for i in {1..4};do   
         NextPolish/bin/bwa index ${input};  
         NextPolish/bin/bwa mem -t 25 ${input} ${reads1} ${reads2} |NextPolish/bin/samtools view  -b - |NextPolish/bin/samtools fixmate -m --threads 5 - - |NextPolish/bin/samtools sort -m 5g --threads 5 - -o ${input}.sort.bam;   
         NextPolish/bin/samtools index ${input}.sort.bam;  
         time -p java -Xmx50G -jar /home/huj/software/pilon-1.23.jar --genome ${input} --frags ${input}.sort.bam --output ${genome}.pilon.v${i} --threads 5 --fix bases;  
         input=${genome}.pilon.v${i}.fasta;  
      done

  -  *Homo sapiens*

    -  work.sh

    .. code-block:: shell

        genome=p_ctg.fa 
        reads1=NGS_1.fq    
        reads2=NGS_1.fq    
        input=${genome}    
        for i in {1..4};do   
           NextPolish/bin/bwa index ${input};  
           NextPolish/bin/bwa mem -t 25 ${input} ${reads1} ${reads2} |NextPolish/bin/samtools view  -b - |NextPolish/bin/samtools fixmate -m --threads 5 - - |NextPolish/bin/samtools sort -m 5g --threads 5 - -o ${input}.sort.bam;   
           NextPolish/bin/samtools index ${input}.sort.bam;
           seqkit split2 -p 20 ${input};
           ls ${input}.split|while read line;do time -p java -Xmx120G -jar /home/huj/software/pilon-1.23.jar --genome ${line} --frags ${input}.sort.bam --output ${line}.pilon --threads 5 --fix bases;done;
           cat ${input}.split/*.pilon.fasta > ${genome}.pilon.v${i}.fasta;
           input=${genome}.pilon.v${i}.fasta;  
        done

  -  Run
     
    .. code-block:: shell

      nohup sh work.sh > pilon.log &

  -  CPU time used for polishing
    
    .. code-block:: shell

      egrep 'user|sys' pilon.log|awk '{x+=$2}END{print x}'

5. **Run NextPolish**

  -  run.cfg
      
  .. code-block:: shell
    
    [General]
    job_type = local
    job_prefix = nextPolish
    task = 1212
    rewrite = yes
    rerun = 3
    parallel_jobs = 5
    multithread_jobs = 5
    genome = p_ctg.fa #miniasm.racon.v4.fasta
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

    egrep 'user|sys' 01_rundir/ */0*.polish.ref.sh.work/polish_genome */nextPolish.sh.e|awk '{print $2}'|sed 's/m//' |sed 's/s//' |awk '{x+=$1* 60+$2}END{print x}'

6. **Run Gmap**
  
  .. code-block:: shell

    genome=miniasm.racon.v4.pilon.v4.fasta # p_ctg.pilon.v4.fasta
    gmap_build -d ./${genome}.gmap ${genome}     
    gmap -D ./ -d ${genome}.gmap Homo_sapiens.GRCh38.cds.all.filter.fa -F -n 1 -i 0 -t 10 -A > ${genome}.gmap.blast

7. **Run Freebayes**
  
  .. code-block:: shell

    genome=miniasm.racon.v4.pilon.v4.fasta # p_ctg.pilon.v4.fasta     
    reads1=NGS_1.fq         
    reads2=NGS_1.fq      
    NextPolish/bin/bwa index ${genome};       
    NextPolish/bin/bwa mem -t 10 ${genome} ${reads1} ${reads2}|NextPolish/bin/samtools view -b - |NextPolish/bin/samtools sort -m 5g --threads 5 - -o ${genome}.bwa.sort.bam;     
    NextPolish/bin/samtools index -@ 10 ${genome}.bwa.sort.bam     
    freebayes -p 2 -b ${genome}.bwa.sort.bam -v ${genome}.sort.bam.vcf -f ${genome}

8. **Count mapped reads**
   
  .. code-block:: python
  
    #!/usr/bin/env python

    import sys
    import pysam

    bam_file = sys.argv[1]
    mapped = full_length_mapped = 0
    for i in pysam.AlignmentFile(bam_file, "r"):
        if i.is_unmapped or i.is_supplementary or i.is_secondary:
            continue
        qseq = i.query_sequence.upper()
        rseq = i.get_reference_sequence().upper()
        mapped += 1
        if qseq == rseq:
            full_length_mapped += 1

    print 'mapped: %d full_length_mapped: %d' % (mapped, full_length_mapped)

9. **Count SNP/Indel**
  
  .. code-block:: shell

    #!/bin/bash

    vcf=$1
    homosnp=$(grep -v '#' ${vcf}|grep snp|grep "1/1"|wc -l)
    echo homosnp: $homosnp

    homoindel=$(grep -v '#' ${vcf}|egrep 'ins|del'|grep "1/1"|wc -l)
    echo homoindel: $homoindel

    hetererrors=$(grep -v '#' ${vcf}|cut -f 10 |sed 's/:/\t/g' |awk '$4==0'|grep -v 1/1 |wc -l)
    echo hetererrors: $hetererrors

10. **Count mapped genes** 
  
  .. code-block:: python
    
    #!/usr/bin/env python

    import sys

    gmap_result_file = sys.argv[1]
    total_gene_count = int(sys.argv[2])
    maps = unmaps = truncate_maps = 0

    names = []
    name = cov = aa = qlen  = ''
    with open(gmap_result_file) as IN:
        for line in IN:
            line = line.strip()
            if not line:
                continue
            lines = line.strip().split()
            if line.startswith('>'):
                if qlen:
                    if int(aa) < int(qlen) * 0.95:
                        truncate_maps += 1
                    qlen = ''
                elif name in names:
                    names.remove(name)

                name = line[1:]
                if name in names:
                    print >>sys.stderr, 'deplicate name: ' + name
                    sys.exit(1)
                else:
                    names.append(name)
            elif line.startswith('Coverage'):
                            qlen = str(int(lines[-2])/3)
            elif line.startswith('Translation'):
                aa = lines[-2][1:]

    if qlen:
        if int(aa) < int(qlen) * 0.95:
            truncate_maps += 1
    elif name in names:
        names.remove(name)

    maps = len(names)
    unmaps = total_gene_count - maps

    print "\t".join(['#','unmap','truncate_map'])
    print "\t".join(map(str, ('#',unmaps,truncate_maps)))

11. **Result can be seen from** `NextPolish paper <https://doi.org/10.1093/bioinformatics/btz891>`__.
