[![Downloads](https://img.shields.io/github/downloads/Nextomics/NextPolish/total.svg)](https://github.com/Nextomics/NextPolish/releases/download/v1.0.5/NextPolish.tgz)
[![Release](https://img.shields.io/github/release/Nextomics/NextPolish.svg)](https://github.com/Nextomics/NextPolish/releases)
[![Issues](https://img.shields.io/github/issues/Nextomics/NextPolish.svg)](https://github.com/Nextomics/NextPolish/issues)

# NextPolish
NextPolish is used to fix base errors (SNV/Indel) in the genome generated by noisy long reads, it can be used with short read data only or long read data only or a combination of both. It contains two core modules, and use a stepwise fashion to correct the error bases in reference genome. To correct the raw third-generation sequencing (TGS) long reads with approximately 15-10% sequencing errors, please use [NextDenovo](https://github.com/Nextomics/NextDenovo).

* **DOWNLOAD**  
click [here](https://github.com/Nextomics/NextPolish/releases/download/v1.0.5/NextPolish.tgz) or use the following command:  
`wget https://github.com/Nextomics/NextPolish/releases/download/v1.0.5/NextPolish.tgz`  

* **REQUIREMENT**
	* [Python 2.7](https://www.python.org/download/releases/2.7/)
	* [Psutil](https://psutil.readthedocs.io/en/latest/)
	* [Drmaa](https://github.com/pygridtools/drmaa-python) (Only required by running under non-local system)

* **INSTALL**  
`tar -vxzf NextPolish.tgz && cd NextPolish && make`

* **UNINSTALL**  
`cd NextPolish && make clean`

* **TEST**  
	`nextPolish test_data/run.cfg`

* **QUICK RUN**  
	1. Prepare sgs_fofn  
	`ls reads1_R1.fq reads1_R2.fq reads2_R1.fq reads2_R2.fq > sgs.fofn`
	2. Create run.cfg  
	`genome=input.genome.fa`  
	`echo -e "task = 1212\ngenome = $genome\nsgs_fofn = sgs.fofn" > run.cfg`
	3. Run  
	`nextPolish run.cfg`
	4. Finally polished genome  
	`cat 03.kmer_count/*polish.ref.sh.work/polish_genome*/genome.nextpolish.part*.fasta > input.genome.nextpolish.v2.fa`  

>***Optional:*** You can also use your own alignment pipeline and then use NextPolish to polish the genome, which will faster than the default NextPolish pipeline when runing on a local system, see [here](./doc/bwa.sh) for an example (using bwa to do alignment).

>***Note:*** If the raw genome generated without a consensus step, such as [miniasm](https://github.com/lh3/miniasm), please run the following command or [racon](https://github.com/isovic/racon) 2-3 rounds using long reads before running [NextPolish](https://github.com/Nextomics/NextPolish) to avoid incorrect mapping of shortgun reads due to the high error rate in the genome assembly.

```bash
    threads=20  
    genome=input.genome.fa
    lgsreads=input.lgs.reads.fq.gz
    bin/minimap2 -a -t ${threads} -x map-ont/map-pb ${genome} ${lgsreads}|bin/samtools view -F 0x4 -b - |bin/samtools sort - -m 2g -@ ${threads} -o genome.lgs.bam;  
    bin/samtools index -@ ${threads} genome.lgs.bam;
    bin/samtools faidx ${genome};
    python lib/nextPolish.py -g ${genome} -t 5 --bam_lgs genome.lgs.bam -p ${threads} > genome.lgspolish.fa;
```
* **USAGE**    
Please see [doc/OPTION.md](doc/OPTION.md) for options introduction.

* **PERFORMANCE COMPARISION**   
	+ [Simulation data](./doc/TEST1.md)
	+ [Actual biological data](./doc/TEST2.md)

* **HELP**   
Please raise an issue at the [issue page](https://github.com/Nextomics/NextPolish/issues/new).

* **CONTACT INFORMATION**    
For additional help, please send an email to huj_at_grandomics_dot_com.

* **COPYRIGHT**    
NextPolish is freely available for academic use and other non-commercial use. 

* **PLEASE STAR AND THANKS**    

* **FAQ**  
	1. What is the difference between [NextPolish](https://github.com/Nextomics/NextPolish) and [Pilon](https://github.com/broadinstitute/pilon)?  
	Currently, NextPolish is focuses on genome correction using shotgun reads, which is also one of the most important steps (typically the last step) to accomplish a genome assembly, while Pilon can be used to make other improvements. For genome correction, NextPolish consumes considerable less time and has a higher correction accuracy for genomes with same sizes and such an advantage becomes more and more significant when the genome size of targeted assemblies increased compared to Pilon. See PERFORMANCE COMPARISION section for more details.
	2. Which job scheduling systems are supported by NextPolish?  
	NextPolish use [DRMAA](https://en.wikipedia.org/wiki/DRMAA) to submit, control, and monitor jobs, so in theory, support all DRMAA-compliant system, such as LOCAL, SGE, PBS, SLURM.
	3. How to continue running unfinished tasks?  
	No need to make any changes, simply run the same command again.
	4. Is it necessary to run steps 3 and 4?  
	In most cases, you can only run steps 1 and 2, steps 3 and 4 are experimental, and we do not currently recommend running on a actual project.
	5. How many iterations to run NextPolish cyclically to get the best result?   
	Our test shown that run NextPolish with 2 iterations, and most of the bases with effectively covered by SGS data can be corrected. Please set task = best to get the best result. Set task = best means NextPolish will cyclically run steps 1 and 2 with 2 iterations. Of course, you can require NextPolish to run with more iterations to get a better result, such as set task=12121212, which means NextPolish will cyclically run steps 1 and 2 with 4 iterations.
	6. What is the difference between bwa or minimap2 to do SGS data mapping?  
	Our test shown Minimap2 is about 3 times faster than bwa, but the accuracy of polished genomes using minimap2 or bwa is tricky, depending on the error rate of genomes and SGS data, see [here](https://lh3.github.io/2018/04/02/minimap2-and-the-future-of-bwa) for more details.
	7. How to specify the queue name/cpu/memory/bash to submit jobs?  
	Please use cluster_options, NextPolish will replace {vf}, {cpu}, {bash} with specific values needed for each jobs.
	8. RuntimeError: Could not find drmaa library.  Please specify its full path using the environment variable DRMAA_LIBRARY_PATH.   
	Please setup the environment variable: DRMAA_LIBRARY_PATH, see [here](https://github.com/pygridtools/drmaa-python) for more details.
	9. ERROR: drmaa.errors.DeniedByDrmException: code 17: error: no suitable queues.    
	This is usually caused by a wrong setting of cluster_options, please check cluster_options first. If you use SGE, you also can add '-w n' to cluster_options, it will switch off validation for invalid resource requests. Please add a similar option for other job scheduling systems. 
	10. OSError: /path/lib64/libc.so.6: version `GLIBC_2.14' not found (required by /path/NextPolish/lib/calgs.so).  
	Please download [this version](https://github.com/Nextomics/NextPolish/releases/download/v1.0.3/NextPolish-CentOS6.9.tgz) and try again.
