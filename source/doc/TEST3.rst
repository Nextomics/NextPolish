.. _simulate_long_noisy_reads:

.. title:: Simulated long noisy reads

Performance comparison between NextPolish and Racon using simulated long noisy reads
------------------------------------------------------------------------------------

**REQUIREMENT**

  -  `PBSIM v1.0.4 <https://github.com/pfaucon/PBSIM-PacBio-Simulator>`__
  -  `NanoSim v2.6.0 <https://github.com/bcgsc/NanoSim>`__
  -  `minimap2 v2.15-r915-dirty <https://github.com/lh3/minimap2>`__
  -  `miniasm v0.3-r179 <https://github.com/lh3/miniasm>`__
  -  `gfatools v0.4-r179-dirty <https://github.com/lh3/gfatools>`__
  -  `samtools v1.9 <https://github.com/samtools/samtools>`__
  -  `Racon v1.3.3 <https://github.com/isovic/racon>`__
  -  `NextPolish v1.2.2 <https://github.com/Nextomics/NextPolish>`__
  -  `Quast v5.0.2 <https://github.com/ablab/quast>`__

1. **Download reference**

  .. code-block:: shell

    curl -SL ftp://ftp.ensembl.org/pub/release-96/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.1.fa.gz | gunzip - > chr01.fa

2. **Simulate PacBio data**

  .. code-block:: shell

    pbsim --data-type CLR --model_qc /PBSIM-PacBio-Simulator/data/model_qc_clr --depth 50 --length-mean 10000 --accuracy-mean 0.85 --prefix pacbio chr01.fa

3. **Simulate NanoPore data**

  .. code-block:: shell

    python NanoSim/src/simulator.py genome -rg chr01.fa -c NanoSim/pre-trained_models/human_NA12878_DNA_FAB49712_guppy/training -n 1631727 -b guppy
    cat simulated_aligned_reads.fasta simulated_unaligned_reads.fasta > ont.sumulated.reads.fa

4. **Assemble reference**
  
  - PacBio data

  .. code-block:: shell

    minimap2 -t 30 -x ava-pb pb.sumulated.reads.fa pb.sumulated.reads.fa > pb.asm.paf 
    miniasm -f pb.sumulated.reads.fa pb.asm.paf > pb.asm.gfa 
    gfatools gfa2fa pb.asm.gfa > pb.asm.fa

  - NanoPore data
  
  .. code-block:: shell

    minimap2 -t 30 -x ava-ont ont.sumulated.reads.fa ont.sumulated.reads.fa > ont.asm.paf 
    miniasm -f ont.sumulated.reads.fa ont.asm.paf > ont.asm.gfa 
    gfatools gfa2fa ont.asm.gfa > ont.asm.fa


5. **Run Racon**   

  - PacBio data

  .. code-block:: shell

    minimap2 -x map-pb -t 20 pb.asm.fa pb.sumulated.reads.fa > pb.map.paf
    racon -t 20 pb.sumulated.reads.fa pb.map.paf pb.asm.fa > pb.asm.racon1.fa

  - NanoPore data

  .. code-block:: shell

    minimap2 -x map-ont -t 20 ont.asm.fa ont.sumulated.reads.fa > ont.map.paf 
    racon -t 20 ont.sumulated.reads.fa ont.map.paf ont.asm.fa > ont.asm.racon1.fa

6. **Run NextPolish**  
  
  - PacBio data

  .. code-block:: shell

    minimap2 -ax map-pb -t 20 pb.asm.fa pb.sumulated.reads.fa|samtools sort - -m 2g --threads 20 -o pb.map.bam 
    samtools index pb.map.bam 
    ls `pwd`/pb.map.bam > pb.map.bam.fofn 
    python NextPolish/lib/nextpolish2.py -g pb.asm.fa -l pb.map.bam.fofn -r clr -p 20 -sp -o pb.asm.nextpolish1.fa

  - NanoPore data
  
  .. code-block:: shell

    minimap2 -ax map-ont -t 20 ont.asm.fa ont.sumulated.reads.fa|samtools sort - -m 2g --threads 20 -o ont.map.bam 
    samtools index ont.map.bam 
    ls `pwd`/ont.map.bam > ont.map.bam.fofn 
    python NextPolish/lib/nextpolish2.py -g ont.asm.fa -l ont.map.bam.fofn -r ont -p 20 -sp -o ont.asm.nextpolish1.fa

  .. note:: Here we use a custom alignment pipeline and then use NextPolish to polish the genome. The genome accuracy after polishing is the same as using NextPolish pipeline to do alignment, see :ref:`Tutorial <long_read_polish>`.

7. **Run Quast**

  - Input
  
    - PacBio data

      -  ``pb.asm.fa``
      -  ``pb.asm.nextpolish1.fa``
      -  ``pb.asm.racon1.fa``

    - NanoPore data

      -  ``ont.asm.fa``
      -  ``ont.asm.nextpolish1.fa``
      -  ``ont.asm.racon1.fa``

  - Run

  .. code-block:: shell

    quast.py --eukaryote --large --threads 25 --min-identity 85 -r chr01.fa pb.asm.fa pb.asm.nextpolish1.fa pb.asm.racon1.fa ont.asm.fa  ont.asm.nextpolish1.fa ont.asm.racon1.fa

  .. object:: Quast result

  +------------------------+----------+------------------+-------------+----------+-------------------+--------------+
  |                        |pb.asm    |pb.asm.nextpolish1|pb.asm.racon1|ont.asm   |ont.asm.nextpolish1|ont.asm.racon1|
  +========================+==========+==================+=============+==========+===================+==============+
  |Total length (>= 0 bp)  |238893883 |229392481         |231583305    |221739507 |231851442          |231932961     |
  +------------------------+----------+------------------+-------------+----------+-------------------+--------------+
  |Reference length        |248956422 |248956422         |248956422    |248956422 |248956422          |248956422     |          
  +------------------------+----------+------------------+-------------+----------+-------------------+--------------+
  |Unaligned length        |1002739   |307941            |70526        |6235359   |6163688            |6431927       |            
  +------------------------+----------+------------------+-------------+----------+-------------------+--------------+
  |Largest alignment       |26588612  |25515573          |25771470     |30803348  |32268337           |32271759      |        
  +------------------------+----------+------------------+-------------+----------+-------------------+--------------+
  |# mismatches per 100 kbp|5425.25   |165.25            |115.42       |4973.49   |30.79              |34.63         |
  +------------------------+----------+------------------+-------------+----------+-------------------+--------------+
  |# indels per 100 kbp    |7127.93   |631.97            |1233.12      |4126.88   |43.39              |83.87         |
  +------------------------+----------+------------------+-------------+----------+-------------------+--------------+
  |# mismatches            |12141134  |370583            |258809       | 11129037 |68890              |77504         |
  +------------------------+----------+------------------+-------------+----------+-------------------+--------------+
  |# indels                |15951531  |1417256           |2765093      |9234603   |97088              |187713        |
  +------------------------+----------+------------------+-------------+----------+-------------------+--------------+

  .. note:: The complete result of Quast can be seen from :download:`here <./TEST3.pdf>`.
