.. _parameterreference:

NextPolish Parameter Reference
==============================

NextPolish requires at least one assembly file (option: ``genome``) and one read file list (option: ``sgs_fofn`` or ``lgs_fofn`` or ``hifi_fofn``) as input, it works with gzip'd FASTA and FASTQ formats and uses a ``config file`` to pass options.

Input
-----

- ``genome file``
  
  .. code-block:: shell
    
    genome=/path/to/need_to_be_polished_assembly_file

- ``read file list`` (one file one line, paired-end files should be interleaved)
  
  .. code-block:: shell
    
    ls reads1_R1.fq reads1_R2.fq reads2_R1.fq.gz reads2_R2.fq.gz ... > sgs.fofn

- ``config file``
  
  A config file is a text file that contains a set of parameters (key=value pairs) to set runtime parameters for NextPolish. The following is a typical config file, which is also located in ``doc/run.cfg``.

  .. code-block:: shell

    [General]
    job_type = local
    job_prefix = nextPolish
    task = best
    rewrite = yes
    deltmp = yes
    rerun = 3
    parallel_jobs = 6
    multithread_jobs = 5
    genome = ./raw.genome.fasta
    genome_size = auto
    workdir = ./01_rundir
    polish_options = -p {multithread_jobs}

    [sgs_option] #optional
    sgs_fofn = ./sgs.fofn
    sgs_options = -max_depth 100 -bwa

    [lgs_option] #optional
    lgs_fofn = ./lgs.fofn
    lgs_options = -min_read_len 1k -max_depth 100
    lgs_minimap2_options = -x map-ont

    [hifi_option] #optional
    hifi_fofn = ./hifi.fofn
    hifi_options = -min_read_len 1k -max_depth 100
    hifi_minimap2_options = -x asm20

Output
------

- ``genome.nextpolish.fasta`` 
  
  Polished genome with fasta format, the fasta header includes primary seqID, length. A lowercase letter indicates a low quality base after polishing, this usually caused by heterozygosity.
- ``genome.nextpolish.fasta.stat``

  Some basic statistical information of the polished genome. 

.. _options:

Options
-------

Global options
##############

  .. option:: job_type = sge           

    local, sge, pbs... (default: sge)
  .. option:: job_prefix = nextPolish  

    prefix tag for jobs. (default: nextPolish)
  .. option:: task = best              

    task need to run [all, default, best, 1, 2, 5, 12, 1212...], 1, 2 are different algorithm modules for short reads, while 5 is the algorithm module for long reads, all=[5]1234, default=[5]12, best=[55]1212. (default: best)
  .. option:: rewrite = no             

    overwrite existed directory [yes, no]. (default: no)
  .. option:: deltmp = yes      

    delete intermediate results. (default: yes)
  .. option:: rerun = 3                

    re-run unfinished jobs untill finished or reached ${rerun} loops, 0=no. (default: 3)
  .. option:: parallel_jobs = 6        

    number of tasks used to run in parallel. (default: 6)
  .. option:: multithread_jobs = 5     

    number of threads used to in a task. (default: 5)
  .. option:: submit = auto   

    command to submit a job, auto = automatically set by `Paralleltask <https://github.com/moold/ParallelTask>`__.
  .. option:: kill = auto   

    command to kill a job, auto = automatically set by `Paralleltask <https://github.com/moold/ParallelTask>`__.
  .. option:: check_alive = auto   

    command to check a job status, auto = automatically set by `Paralleltask <https://github.com/moold/ParallelTask>`__.
  .. option:: job_id_regex = auto   

    the job-id-regex to parse the job id from the out of ``submit``, auto = automatically set by `Paralleltask <https://github.com/moold/ParallelTask>`__.
  .. option:: use_drmaa = no   

    use drmaa to submit and control jobs.
  .. option:: genome = genome.fa       

    genome file need to be polished. (**required**)
  .. option:: genome_size = auto       

    genome size, auto = calculate genome size using the input ${genome} file. (default: auto)
  .. option:: workdir = 01_rundir      

    work directory. (default: ./)
  .. option:: polish_options = -p {multithread_jobs}

    ::
    
      -p, number of processes used for polishing.
      -u, output uppercase sequences. (default: False)
      -debug, output details of polished bases to stderr, only useful in short read polishing. (default: False)

Options for short reads
#######################

  .. option:: sgs_fofn = ./sgs.fofn    

    input short read files list, one file one line, paired-end files should be interleaved.
  .. option:: sgs_options = -max_depth 100 -bwa

    ::

      -N, don't discard a read/pair if the read contains N base.
      -use_duplicate_reads, use duplicate pair-end reads in the analysis. (default: False)
      -unpaired, unpaired input files. (default: False)
      -max_depth, use up to ${max_depth} fold reads data to polish. (default: 100)
      -bwa, use bwa to do mapping. (default: -bwa) 
      -minimap2, use minimap2 to do mapping, which is much faster than bwa. 

Options for long reads
#######################

  .. option:: lgs_fofn = ./lgs.fofn    

    input long read files list, one file one line.             
  .. option:: lgs_options = -min_read_len 1k -max_depth 100

    ::

      -min_read_len, filter reads with length shorter than ${min_read_len}. (default: 1k)
      -max_read_len, filter reads with length longer than $ {max_read_len}, ultra-long reads usually contain lots of errors, and the mapping step requires significantly more memory and time, 0=disable (default: 0)
      -max_depth, use up to ${max_depth} fold reads data to polish, 0=disable. (default: 100)  
  .. option:: lgs_minimap2_options = -x map-pb -t {multithread_jobs}
      
    minimap2 options, used to set PacBio/Nanopore reads mapping. (**required**)

Options for hifi reads
#######################

  .. option:: hifi_fofn = ./hifi.fofn    

    input hifi read files list, one file one line.             
  .. option:: hifi_options = -min_read_len 1k -max_depth 100

    ::

      -min_read_len, filter reads with length shorter than ${min_read_len}. (default: 1k)
      -max_read_len, filter reads with length longer than $ {max_read_len}, ultra-long reads usually contain lots of errors, and the mapping step requires significantly more memory and time, 0=disable (default: 0)
      -max_depth, use up to ${max_depth} fold reads data to polish, 0=disable. (default: 100)  
  .. option:: hifi_minimap2_options = -x map-pb -t {multithread_jobs}
      
    minimap2 options, used to set hifi reads mapping. (**required**)