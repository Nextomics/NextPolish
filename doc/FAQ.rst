.. _faq:

Frequently Asked Questions
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. contents::
  :local:

What is the difference between `NextPolish <https://github.com/Nextomics/NextPolish>`__ and `Pilon <https://github.com/broadinstitute/pilon>`__?
----------------------------------------------------------------------------------------------------------------------------------------------------------

Currently, NextPolish is focuses on genome correction using shotgun reads, which is also one of the most important steps (typically the last step) to accomplish a genome assembly, while Pilon can be used to make other improvements. For genome correction, NextPolish consumes considerable less time and has a higher correction accuracy for genomes with same sizes and such an advantage becomes more and more significant when the genome size of targeted assemblies increased compared to Pilon. See BENCHMARK section for more details.

Which job scheduling systems are supported by NextPolish?
------------------------------------------------------------

NextPolish use `Paralleltask <https://github.com/moold/ParallelTask>`__ to submit, control, and monitor jobs, so in theory, support all Paralleltask-compliant system, such as LOCAL, SGE, PBS, SLURM.

How to continue running unfinished tasks?
--------------------------------------------

No need to make any changes, simply run the same command again.

How to set the `task` parameter?
-------------------------------------

The ``task`` parameter is used to set the polishing algorithm logic, 1, 2, 3, 4 are different algorithm modules for short reads, while 5 is the algorithm module for long reads. BTW, steps 3 and 4 are experimental, and we do not currently recommend running on a actual project. Set ``task=551212`` means NextPolish will cyclically run steps 5, 1 and 2 with 2 iterations.

How many iterations to run NextPolish cyclically to get the best result?
---------------------------------------------------------------------------

Our test shown that run NextPolish with 2 iterations, and most of the bases with effectively covered by SGS data can be corrected. Please set ``task=best`` to get the best result. ``task = best`` means NextPolish will cyclically run steps [5], 1 and 2 with 2 iterations. Of course, you can require NextPolish to run with more iterations to get a better result, such as set ``task=555512121212``, which means NextPolish will cyclically run steps 5, 1 and 2 with 4 iterations.

Why does the contig N50 of polished genome become shorter or why does the polished genome contains some extra ``N``?
--------------------------------------------------------------------------------------------------------------------------

In some cases, if the short reads contain ``N``, some error bases will be fixed by ``N`` (the global score of a kmer with ``N`` is the largest and be selected), and remove ``N`` in short reads will avoid this.

What is the difference between bwa or minimap2 to do SGS data mapping?
--------------------------------------------------------------------------

Our test shown Minimap2 is about 3 times faster than bwa, but the accuracy of polished genomes using minimap2 or bwa is tricky, depending on the error rate of genomes and SGS data, see `here <https://lh3.github.io/2018/04/02/minimap2-and-the-future-of-bwa>`__ for more details.

How to specify the queue/cpu/memory/bash to submit jobs?   
------------------------------------------------------------
See `here <https://github.com/moold/ParallelTask#configuration>`__ to edit the `Paralleltask <https://github.com/moold/ParallelTask>`__ configure template file ``cluster.cfg``, or use the ``submit`` parameter.