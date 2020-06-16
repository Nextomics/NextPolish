## Performance comparison between NextPolish and Racon using simulation data  

* **REQUIREMENT**
    * [PBSIM v1.0.4](https://github.com/pfaucon/PBSIM-PacBio-Simulator)
    * [NanoSim v2.6.0](https://github.com/bcgsc/NanoSim)
    * [minimap2 v2.15-r915-dirty](https://github.com/lh3/minimap2)
    * [miniasm v0.3-r179](https://github.com/lh3/miniasm)
    * [gfatools v0.4-r179-dirty](https://github.com/lh3/gfatools)
    * [samtools v1.9](https://github.com/samtools/samtools)
    * [Racon v1.3.3](https://github.com/isovic/racon)
    * [NextPolish v1.2.2](https://github.com/Nextomics/NextPolish)
    * [Quast v5.0.2](https://github.com/ablab/quast)
    
* **Download reference**   
```
curl -SL ftp://ftp.ensembl.org/pub/release-96/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.1.fa.gz | gunzip - > chr01.fa
```

* **Simulate PacBio data**  
```
pbsim --data-type CLR --model_qc /PBSIM-PacBio-Simulator/data/model_qc_clr --depth 50 --length-mean 10000 --accuracy-mean 0.85 --prefix pacbio chr01.fa
```

* **Simulate NanoPore data**  
```
python NanoSim/src/simulator.py genome -rg chr01.fa -c NanoSim/pre-trained_models/human_NA12878_DNA_FAB49712_guppy/training -n 1631727 -b guppy
cat simulated_aligned_reads.fasta simulated_unaligned_reads.fasta > ont.sumulated.reads.fa
```

* **Assemble reference**  
```
minimap2 -t 30 -x ava-pb pb.sumulated.reads.fa pb.sumulated.reads.fa > pb.asm.paf
miniasm -f pb.sumulated.reads.fa pb.asm.paf > pb.asm.gfa
gfatools gfa2fa pb.asm.gfa > pb.asm.fa

minimap2 -t 30 -x ava-ont ont.sumulated.reads.fa ont.sumulated.reads.fa > ont.asm.paf
miniasm -f ont.sumulated.reads.fa ont.asm.paf > ont.asm.gfa
gfatools gfa2fa ont.asm.gfa > ont.asm.fa
```

* **Run Racon**   
```
minimap2 -x map-pb -t 20 pb.asm.fa pb.sumulated.reads.fa > pb.map.paf
racon -t 20 pb.sumulated.reads.fa pb.map.paf pb.asm.fa > pb.asm.racon1.fa

minimap2 -x map-ont -t 20 ont.asm.fa ont.sumulated.reads.fa > ont.map.paf
racon -t 20 ont.sumulated.reads.fa ont.map.paf ont.asm.fa > ont.asm.racon1.fa
```

* **Run NextPolish**  
```
minimap2 -ax map-pb -t 20 pb.asm.fa pb.sumulated.reads.fa|samtools sort - -m 2g --threads 20 -o pb.map.bam
samtools index pb.map.bam
ls `pwd`/pb.map.bam > pb.map.bam.fofn
python NextPolish/lib/nextpolish2.py -g pb.asm.fa -l pb.map.bam.fofn -r clr -p 20 -a -s -o pb.asm.nextpolish1.fa

minimap2 -ax map-ont -t 20 ont.asm.fa ont.sumulated.reads.fa|samtools sort - -m 2g --threads 20 -o ont.map.bam
samtools index ont.map.bam
ls `pwd`/ont.map.bam > ont.map.bam.fofn
python NextPolish/lib/nextpolish2.py -g ont.asm.fa -l ont.map.bam.fofn -r ont -p 20 -s -a -o ont.asm.nextpolish1.fa

```

* **Run Quast**  
    + Input    
        * pb.asm.fa
        * pb.asm.nextpolish1.fa
        * pb.asm.racon1.fa
        * ont.asm.fa
        * ont.asm.nextpolish1.fa
        * ont.asm.racon1.fa
    + Run    
    ```
    quast.py --eukaryote --large --threads 25 --min-identity 85 -r chr01.fa pb.asm.fa pb.asm.nextpolish1.fa pb.asm.racon1.fa ont.asm.fa  ont.asm.nextpolish1.fa ont.asm.racon1.fa
    ```
    + **Click [here](TEST3.pdf) for result**
    

