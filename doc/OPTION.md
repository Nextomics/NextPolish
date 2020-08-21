# NextPolish

NextPolish requires at least one assembly file (option: genome) and one read files list (option: sgs_fofn or lgs_fofn) as input, it works with gzip'd FASTA and FASTQ formats and uses a config file to pass options, see [here](./run.cfg) for an example.

* **INPUT**    
    - genome file  
    `genome=/path/to/need_to_be_polished_assembly_file`
    - reads files list (one file one line, paired-end files should be interleaved)  
    `ls reads1_R1.fq reads1_R2.fq reads2_R1.fq.gz reads2_R2.fq.gz ... > sgs.fofn`

* **OUTPUT**    
    - genome.nextpolish.fasta with fasta format, the fasta header includes primary seqID, length.
    - genome.nextpolish.fasta.stat, some basic statistical information of the polished genome.

* **OPTION** 

<pre>
    [General]                # global options
    job_type = sge           # [local, sge, pbs...]. (default: sge)
    job_prefix = nextPolish  # prefix tag for jobs. (default: nextPolish)
    task = best              # task need to run [all, default, best, 1, 2, 5, 12, 1212...], 1, 2 are different algorithm modules for short reads, while 5 is the algorithm module for long reads, all=[5]1234, default=[5]12, best=[55]1212. (default: best)
    rewrite = no             # overwrite existed directory [yes, no]. (default: no)
    rerun = 3                # re-run unfinished jobs untill finished or reached ${rerun} loops, 0=no. (default: 3)
    parallel_jobs = 6        # number of tasks used to run in parallel. (default: 6)
    multithread_jobs = 5     # number of threads used to in a task. (default: 5)
    cluster_options = auto   # a template to define the resource requirements for each job, which will pass to DRMAA as the nativeSpecification field.
    genome = genome.fa       # genome file need to be polished. (<b>required</b>)
    genome_size = auto       # genome size, auto = calculate genome size using the input ${genome} file. (default: auto)
    workdir = 01_rundir      # work directory. (default: ./)
<!--    round_count = 1          # number of iterations to run NextPolish cyclically. (default: 1)
    round_mode = 2           # preset mode of iterations to run NextPolish cyclically, 1 = 1234[1234], 2 = 12[12]34, 3 = 123[123]4. (default: 2) -->
    polish_options = -p {multithread_jobs}
                             # options used in polished step, see below.

    [sgs_option]             # options of short reads. (<b>optional</b>)
    sgs_fofn = ./sgs.fofn    # input short reads file, one file one line, paired-end files should be interleaved.
    sgs_options = -max_depth 100 -bwa
                             # -use_duplicate_reads, use duplicate pair-end reads in the analysis. (default: False)
                             # -unpaired, unpaired input files. (default: False)
                             # -max_depth, use up to ${max_depth} fold reads data to polish. (default: 100)
                             # -bwa, use bwa to do mapping. (default: -bwa) 
                             # -minimap2, use minimap2 to do mapping, which is much faster than bwa. 

    [lgs_option]             # options of long reads. (<b>optional</b>)
    lgs_fofn = ./lgs.fofn    # input long reads file, one file one line.             
    lgs_options = -min_read_len 1k -max_depth 100
                             # -min_read_len, filter reads with length shorter than ${min_read_len}. (default: 1k)
                             # -max_read_len, filter reads with length longer than $ {max_read_len}, ultra-long reads usually contain lots of errors, and the mapping step requires significantly more memory and time, 0=disable (default: 0)
                             # -max_depth, use up to ${max_depth} fold reads data to polish, 0=disable. (default: 100)
    lgs_minimap2_options = -x map-pb -t {multithread_jobs}
                             # minimap2 options, used to set PacBio/Nanopore read overlap. (<b>required</b>)
    
    [polish_options]         # options used in polished step.
    -debug                   # output details of polished bases to stderr. (default: False)

    algorithm arguments:
    -count_read_ins_sgs      # read ${count_read_ins_sgs} reads to estimate the insert size of paired-end reads. (default: 10000)
    -min_map_quality         # skip the mapped read with mapping quality < ${min_map_quality}. (default: 0)
    -max_ins_len_sgs         # skip the paired-end read with insert size > ${max_ins_len_sgs}. (default: 10000)
    -max_ins_fold_sgs        # skip the paired-end read with insert size > ${max_ins_fold_sgs} * estimated_average_insert_size. (default: 5)
    -max_clip_ratio_sgs      # skip the mapped read with clipped length > ${max_clip_ratio_sgs} * full_length, used for bam_sgs. (default: 0.15)
    -max_clip_ratio_lgs      # skip the mapped read with clipped length > ${max_clip_ratio_lgs} * full_length, used for bam_lgs. (default: 0.4)
    -trim_len_edge           # trimed length at the two edges of an alignment. (default: 2)
    -ext_len_edge            # extened length at the two edges of a low quality region. (default: 2)

    score_chain:
    -indel_balance_factor_sgs 
                             # a factor to control the ratio between indels, larger factor will produced more deletions, and vice versa. (default: 0.5)
    -min_count_ratio_skip    # skip a site if the fraction of the most genotype > ${min_count_ratio_skip}. (default: 0.8)

    kmer_count:
    -min_len_ldr             # minimum length requirement of a low depth region, which will be further processed using bam_lgs. (default: 3)
    -max_len_kmer            # maximum length requirement of a polished kmer, longer kmers will be splited. (default: 50)
    -min_len_inter_kmer      # minimum interval length between two adjacent kmers, shorter interval length will be merged. (default: 5)
    -max_count_kmer          # read up to this count of observed kmers for a polished kmer. (default: 50)

    snp_phase:
    -ploidy                  # set the ploidy of the sample of this genome. (default: 2)
    -max_variant_count_lgs   # exclude long reads with more than ${max_variant_count_lgs} variable sites, it is approximately equivalent to total error bases in the long read. (default: 150k)
    -indel_balance_factor_lgs 
                             # a factor to control the ratio between indels, larger factor will produced more deletions, and vice versa. (default: 0.33)
    -min_depth_snp           # recall snps using bam_lgs if the total depth of this site in bam_sgs < ${min_depth_snp}. (default: 3)
    -min_count_snp           # recall snps using bam_lgs if the count of this snp in bam_sgs < ${min_count_snp}. (default: 5)
    -min_count_snp_link      # find a snp linkage using bam_lgs if the count of this linkage in bam_sgs < ${min_count_snp_link}. (default: 5)
    -max_indel_factor_lgs    # recall indels with bam_sgs if the count of the second most genotype > ${max_indel_factor_lgs} * the count of the most genotype when the most genotype is different with ref in bam_lgs. (default: 0.21)
    -max_snp_factor_lgs      # recall snps with bam_lgs if the count of the second most genotype > ${max_snp_factor_lgs} * the count of the most genotype when the most genotype is different with ref. (default: 0.53)
    -min_snp_factor_sgs      # skip a snp if the count of the second most genotype < ${min_snp_factor_sgs} * the count of the most genotype. (default: 0.34)
</pre>
