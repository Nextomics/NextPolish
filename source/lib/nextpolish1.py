#!/usr/bin/env python
from __future__ import print_function

import os
import sys
import argparse
from multiprocessing import Pool
from kit import *
from ctypes import *

log = plog()

class HelpFormatter(
		argparse.RawDescriptionHelpFormatter,
		argparse.ArgumentDefaultsHelpFormatter):
	pass

class faidx_t (Structure):
	pass

class hts_idx_t (Structure):
	pass

class bam_hdr_t (Structure):
	pass

class Configure (Structure):
	_fields_ = [
		("trim_len_edge", c_uint8),
		("ext_len_edge", c_uint8),
		("min_map_quality", c_uint8),
		("indel_balance_factor_sgs", c_double),
		("min_count_ratio_skip", c_double),

		("min_len_ldr", c_uint8),
		("min_len_inter_kmer", c_uint8),
		("max_len_kmer", c_uint8),
		("max_count_kmer", c_uint8),

		("min_depth_snp", c_uint8),
		("min_count_snp", c_uint8),
		("min_count_snp_link", c_int8),
		("ploidy", c_double),
		("indel_balance_factor_lgs", c_double),
		("max_indel_factor_lgs", c_double),
		("max_snp_factor_lgs", c_double),
		("min_snp_factor_sgs", c_double),

		("region_count", c_int32),
		("count_read_ins_sgs", c_uint32),
		("max_ins_len_sgs", c_uint32),
		("max_ins_fold_sgs", c_int32),
		("max_variant_count_lgs", c_int32),

		("max_clip_ratio_sgs", c_double),
		("max_clip_ratio_lgs", c_double),
		
		("trace_polish_open", c_int32),

		("read_tlen", c_int32),
		("read_len", c_int32),
		("fastafn", c_char_p),
		("bamfn", c_char_p),
		("thirdbamfn", c_char_p)
	]

class PolishPoint (Structure):
	_fields_ = [
		("pos", c_int32),
		("index", c_int16),
		("curbase", c_char),
		("base", c_char)
	]

class PolishResult (Structure):
	_fields_ = [
		("contig", c_char_p),
		("data", POINTER(PolishPoint)),
		("length", c_int32),
		("datalength", c_int32)
	]

global P, CFG, FUN
P = CDLL(os.path.dirname(os.path.realpath(__file__)) + '/' + "nextpolish1.so")
P.config_init.argtypes = [c_char_p, c_char_p, c_char_p]
P.config_init.restype = POINTER(Configure)
P.config_destory.argtypes = [POINTER(Configure)]

P.polishresult_destory.argtypes = [POINTER(PolishResult)]

P.score_chain.argtypes = [c_char_p, POINTER(Configure)]
P.score_chain.restype = POINTER(PolishResult)
P.kmer_count.argtypes = [c_char_p, POINTER(Configure)]
P.kmer_count.restype = POINTER(PolishResult)
P.snp_phase.argtypes = [c_char_p, POINTER(Configure)]
P.snp_phase.restype = POINTER(PolishResult)
P.snp_valid.argtypes = [c_char_p, POINTER(Configure)]
P.snp_valid.restype = POINTER(PolishResult)
P.lgspolish.argtypes = [c_char_p, POINTER(Configure)]
P.lgspolish.restype = POINTER(PolishResult)

def update_cfg(args):
	global CFG
	CFG[0].trim_len_edge = args.trim_len_edge
	CFG[0].ext_len_edge = args.ext_len_edge
	CFG[0].min_map_quality = args.min_map_quality
	CFG[0].indel_balance_factor_sgs = args.indel_balance_factor_sgs
	CFG[0].min_count_ratio_skip = args.min_count_ratio_skip

	CFG[0].min_len_ldr = args.min_len_ldr
	CFG[0].min_len_inter_kmer = args.min_len_inter_kmer
	CFG[0].max_len_kmer = args.max_len_kmer
	CFG[0].max_count_kmer = args.max_count_kmer

	CFG[0].min_depth_snp = args.min_depth_snp
	CFG[0].min_count_snp = args.min_count_snp
	CFG[0].min_count_snp_link = args.min_count_snp_link
	CFG[0].ploidy = args.ploidy
	CFG[0].indel_balance_factor_lgs = args.indel_balance_factor_lgs
	CFG[0].max_indel_factor_lgs = args.max_indel_factor_lgs
	CFG[0].max_snp_factor_lgs = args.max_snp_factor_lgs
	CFG[0].min_snp_factor_sgs = args.min_snp_factor_sgs

	CFG[0].region_count = 10000
	CFG[0].count_read_ins_sgs = args.count_read_ins_sgs
	CFG[0].max_ins_len_sgs = args.max_ins_len_sgs
	CFG[0].max_ins_fold_sgs = args.max_ins_fold_sgs
	CFG[0].max_variant_count_lgs = args.max_variant_count_lgs

	CFG[0].max_clip_ratio_sgs = args.max_clip_ratio_sgs
	CFG[0].max_clip_ratio_lgs = args.max_clip_ratio_lgs

	CFG[0].trace_polish_open = 1 if args.debug else 0

def correct(seq_name, genome, sgs_bam, lgs_bam):
	seq = ''
	lable = 0
	with open(genome) as IN:
		for line in IN:
			if line.startswith('>'):
				if line[1:].strip().split()[0] == seq_name:
					lable = 1
				elif lable:
					return seq
			elif lable:
				seq += line.strip()

def read_unpolished_seqs(infile, index, polished_seqs):
	unpolished_seqs = set()
	if args.block_index != 'all':
		with open(infile) as IN:
			for line in IN:
				lines = line.strip().split()
				if lines and lines[0].split('_np')[0] not in polished_seqs and lines[1] == index:
					unpolished_seqs.add(str2byte(lines[0]))
	else:
		with open(infile) as IN:
			for line in IN:
				if line.startswith('>'):
					unpolished_seqs.add(str2byte(line.strip().split()[0][1:]))
	return unpolished_seqs

def read_polished_seqs(infile, polished_seqs):
	last_seq = ''
	cur_seq_offset = last_seq_position = 0
	with open(infile) as IN:
		for line in IN:
			if line.startswith('>'):
				last_seq_position += cur_seq_offset
				cur_seq_offset = len(line)
				last_seq = seq_name = line.split()[0].split('_np')[0][1:]
				log.warning('Skip polished seq: ' + seq_name)
				polished_seqs.add(seq_name)
			else:
				cur_seq_offset += len(line)

	if last_seq:
		polished_seqs.remove(last_seq)
	return last_seq_position

def worker(seed_name):
	c_seq = FUN(seed_name, CFG)

	seq = string_at(c_seq.contents.contig)
	seq_len = c_seq.contents.length
	pdata = c_seq.contents.data
	pbases = [(pdata[p].pos, pdata[p].index, byte2str(pdata[p].curbase), byte2str(pdata[p].base)) for p in range(c_seq.contents.datalength)]
	P.polishresult_destory(c_seq)
	return (byte2str(seed_name), byte2str(seq), seq_len, pbases)

def start():
	log.info(
		'Start a polished worker in %d from parent %d' %
		(os.getpid(), os.getppid()))

def main(args):
	log.info('Polished step options:')
	log.info(args)

	OUT = sys.stdout
	polished_seqs = set()
	last_seq_position = 0
	if args.out != 'stdout':
		if os.path.exists(args.out):
			last_seq_position = read_polished_seqs(args.out, polished_seqs)
			OUT = open(args.out, 'r+')
			OUT.seek(last_seq_position, os.SEEK_SET)
			# OUT.seek(-1 * last_seq_position, 2)
		else:
			OUT = open(args.out, 'w')
	
	blockfile = args.block
	if args.block_index == 'all' or not args.block:
		args.block_index = 'all'
		blockfile = args.genome
	unpolished_seqs = read_unpolished_seqs(blockfile, args.block_index, polished_seqs)
	
	global CFG, FUN
	CFG = P.config_init(str2byte(args.genome), str2byte(args.bam_sgs), str2byte(args.bam_lgs))
	FUN = {1 : P.score_chain, 2 : P.kmer_count, 3 : P.snp_phase, 4 : P.snp_valid, 5: P.lgspolish}[args.task]
	update_cfg(args)

	pool = Pool(args.process, initializer=start)
	for seq_pname, seq, seq_len, pbases in pool.imap_unordered(worker, unpolished_seqs, chunksize=1):
		# if seq_len >= args.min_len_seq:
		if args.uppercase:
			seq = seq.upper()
		seq_name = seq_pname + (str(args.task) if seq_pname.split('_')[-1].startswith('np') else ('_np' + str(args.task)))
		print('>%s %d\n%s' % (seq_name, seq_len, seq), file=OUT)
		for pbase in pbases:
			print(seq_pname + ' %d %d %c %c' % pbase, file=sys.stderr)
		# else:
			# log.warning('discard short seq %s with length %d', seq_name, seq_len)
	pool.close()
	pool.join()
	if args.out != 'stdout':
		OUT.close()
	P.config_destory(CFG)

if __name__ == '__main__':
	parser = argparse.ArgumentParser(
		formatter_class=HelpFormatter,
		description='''
%(prog)s:
	Polish the genome using multi-processor.

exmples:
	%(prog)s -f input.idxs -g genome.fa -s sgs.sort.bam -l lgs.sort.bam -p 10 -t 3
'''
	)

	io = parser.add_argument_group('Input/Output arguments')
	io.add_argument('-g', '--genome', metavar = 'FILE', required=True, type=str,
		help='genome file, the reference of bam alignments, require index file [samtools faidx].')
	io.add_argument('-s', '--bam_sgs', metavar = 'FILE', type=str,
		help='sorted bam file of short reads, require index file [samtools index].')
	io.add_argument('-l', '--bam_lgs', metavar = 'FILE', type=str,
		help='sorted bam file of long reads, require index file [samtools index].')
	io.add_argument('-b', '--block', metavar = 'FILE', type=str,
		help='genome block file, each line includes [seq_id, index].')
	io.add_argument('-i', '--block_index', type=str, default = 'all',
		help='index of seqs need to be polished in genome block file.')
	io.add_argument('-u', '--uppercase', action='store_true', default=False,
		help='output uppercase sequences.')
	io.add_argument('-debug', action='store_true', default=False,
		help='output details of polished bases to stderr.')
	io.add_argument('-o', '--out', metavar = 'FILE', default='stdout',
		help='output file, polished seqs in output file will be skipped.')
	
	algorithm = parser.add_argument_group('Algorithm arguments')
	algorithm.add_argument('-t', '--task', metavar = 'N', type=int, required=True, choices=[1,2,3,4,5],
		help='task need to run, [1=score_chain, 2=kmer_count, 3=snp_phase, 4=snp_valid]')
	algorithm.add_argument('-p', '--process', metavar = 'N', type=int, default=10,
		help='number of processes used for polishing.')
	# algorithm.add_argument('-min_len_seq', metavar = 'N', type=str, default = '10k',
	# 	help='minimum length requirement of a seq, shorter seqs will be skipped.')
	algorithm.add_argument('-count_read_ins_sgs', metavar = 'N', type=int, default=10000,
		help='read N reads to estimate the insert size of paired-end reads.')
	algorithm.add_argument('-min_map_quality', metavar = 'N', type=int, default=0,
		help='skip the mapped read with mapping quality < N.')
	algorithm.add_argument('-max_ins_len_sgs', metavar = 'N', type=int, default=10000,
		help='skip the paired-end read with insert size > N.')
	algorithm.add_argument('-max_ins_fold_sgs', metavar = 'N', type=int, default=5,
		help='skip the paired-end read with insert size > N * estimated average insert size.')
	algorithm.add_argument('-max_clip_ratio_sgs', metavar = 'F', type=float, default=0.15,
		help='skip the mapped read with clipped length > F * full length, used for bam_sgs.')
	algorithm.add_argument('-max_clip_ratio_lgs', metavar = 'F', type=float, default=0.4,
		help='skip the mapped read with clipped length > F * full length, used for bam_lgs.')
	algorithm.add_argument('-trim_len_edge', metavar = 'N', type=int, default=2,
		help='trimed length at the two edges of a alignment.')
	algorithm.add_argument('-ext_len_edge', metavar = 'N', type=int, default=2,
		help='extened length at the two edges of a low quality region.')

	score_chain = parser.add_argument_group('score_chain')
	score_chain.add_argument('-indel_balance_factor_sgs', metavar = 'F', type=float, default=0.5,
		help='a factor to control the ratio between indels, larger factor will produced more deletions, and vice versa.')
	score_chain.add_argument('-min_count_ratio_skip', metavar = 'F', type=float, default=0.8,
		help='skip a site if the fraction of the most genotype > F.')

	kmer_count = parser.add_argument_group('kmer_count')
	kmer_count.add_argument('-min_len_ldr', metavar = 'N', type=int, default=3,
		help='minimum length requirement of a low depth region, which will be further processed using bam_lgs.')
	kmer_count.add_argument('-max_len_kmer', metavar = 'N', type=int, default=50,
		help='maximum length requirement of a polished kmer, longer kmers will be splited.')
	kmer_count.add_argument('-min_len_inter_kmer', metavar = 'N', type=int, default=5,
		help='minimum interval length between two adjacent kmers, shorter interval length will be merged.')
	kmer_count.add_argument('-max_count_kmer', metavar = 'N', type=int, default=50,
		help='read up to this count of observed kmers for a polished kmer.')

	snp_phase = parser.add_argument_group('snp_phase')
	snp_phase.add_argument('-ploidy', metavar = 'N', type=int, default=2,
		help='set the ploidy of the sample of this genome.')
	snp_phase.add_argument('-max_variant_count_lgs', metavar = 'N', type=str, default = '150k',
		help='exclude long reads with more than N variable sites, it is approximately equivalent \
			to total error bases in the long read.')
	snp_phase.add_argument('-indel_balance_factor_lgs', metavar = 'F', type=float, default=0.33,
		help='a factor to control the ratio between indels, larger factor will produced more deletions, and vice versa.')
	snp_phase.add_argument('-min_depth_snp', metavar = 'N', type=int, default=3,
		help='recall snps using bam_lgs if the total depth of this site in bam_sgs < N.')
	snp_phase.add_argument('-min_count_snp', metavar = 'N', type=int, default=5,
		help='recall snps using bam_lgs if the count of this snp in bam_sgs < N.')
	snp_phase.add_argument('-min_count_snp_link', metavar = 'N', type=int, default=5,
		help='find a snp linkage using bam_lgs if the count of this linkage in bam_sgs < N.')
	snp_phase.add_argument('-max_indel_factor_lgs', metavar = 'F', type=float, default=0.21,
		help='recall indels with bam_sgs if the count of the second most genotype > F * \
			the count of the most genotype when the most genotype is different with ref in bam_lgs.')
	snp_phase.add_argument('-max_snp_factor_lgs', metavar = 'F', type=float, default=0.53,
		help='recall snps with bam_lgs if the count of the second most genotype > F * \
			the count of the most genotype when the most genotype is different with ref.')
	snp_phase.add_argument('-min_snp_factor_sgs', metavar = 'F', type=float, default=0.34,
		help='skip a snp if the count of the second most genotype < F * \
			the count of the most genotype.')
	# snp_valid = parser.add_argument_group('snp_valid arguments')
	args, unknown = parser.parse_known_args()

	# args.min_len_seq = parse_num_unit(args.min_len_seq)
	args.max_variant_count_lgs = parse_num_unit(args.max_variant_count_lgs)
	if args.task == 5:
		log.error("\033[35m Please use %s/nextpolish2.py to polish the genome with long reads. \033[35m" % os.path.dirname(os.path.realpath(__file__)))
		sys.exit(1)
	main(args)
