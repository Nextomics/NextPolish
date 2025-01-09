"""Configuration file parser."""
import re, os, sys
from kit import *

__all__ = ["ConfigParser"]

log = plog()

def _bool(v):
	return 1 if str(v).lower() not in ['no', '0', 'false'] and bool(v) else 0

class ConfigParser:
	def __init__(self, cfgfile):
		self.cfg = self._defaultcfg()
		self.cfgdir = os.path.dirname(os.path.abspath(cfgfile))
		self.read(cfgfile)
		self._check()
		self._settask()

	def _defaultcfg(self):
		cfg = {}
		cfg['job_type'] = 'local'
		cfg['job_prefix'] = 'nextPolish'
		cfg['task'] = 'best'
		cfg['rewrite'] = '0'
		cfg['cleantmp'] = '0'
		cfg['rerun'] = '3'
		cfg['parallel_jobs'] = '6'
		cfg['multithread_jobs'] = '5'
		cfg['use_drmaa'] = False
		cfg['deltmp'] = False
		cfg['submit'] = None
		cfg['kill'] = None
		cfg['check_alive'] = None
		cfg['job_id_regex'] = None
		# cfg['cluster_options'] = ''

		cfg['genome_size'] = 'auto'
		cfg['workdir'] = os.getcwd()
		# cfg['round_count'] = '1'
		# cfg['round_mode'] = '2'
		# cfg['min_len_seq'] = '10k'
		cfg['polish_options'] = ''

		cfg['sgs_options'] = '-max_depth 100'
		cfg['sgs_use_duplicate_reads'] = 0
		cfg['sgs_unpaired'] = '0'
		cfg['sgs_max_depth'] = '100'
		cfg['sgs_block_size'] = '500M'
		cfg['sgs_rm_nread'] = 1
		# cfg['sgs_align_options'] = 'minimap2 --split-prefix tmp -a -x sr'
		cfg['sgs_align_options'] = 'bwa mem'

		cfg['lgs_options'] = '-min_read_len 1k -max_depth 100'
		cfg['lgs_block_size'] = '500M'
		cfg['lgs_min_read_len'] = '1k'
		cfg['lgs_max_read_len'] = '0'
		cfg['lgs_max_depth'] = '100'
		cfg['lgs_minimap2_options'] = '-x map-ont'
		cfg['lgs_read_type'] = ''

		cfg['hifi_options'] = '-min_read_len 1k -max_depth 100'
		cfg['hifi_block_size'] = '500M'
		cfg['hifi_min_read_len'] = '1k'
		cfg['hifi_max_read_len'] = '0'
		cfg['hifi_max_depth'] = '100'
		cfg['hifi_minimap2_options'] = '-x map-pb'

		return cfg

	def read(self, cfgfile):
		with open(cfgfile) as IN:
			for line in IN:
				line = line.strip()
				if not line or line[0].startswith('#'):
					continue
				group = re.search(r'([^;\s]+)\s*[=:]\s*([^;#\n]+)(\s*|#.*)$', line) # a option = value1 value2 # annotation
				if group and group.groups()[1].strip():
					self.cfg[group.groups()[0]] = group.groups()[1].strip()

	def _settask(self):
		import re
		task = re.sub(r'[\s,;]+','', self.cfg['task'])
		if task == 'all':
			task = '561234'
		elif task == 'default':
			task = '5612'
		elif task == 'best':
			task = '55661212'

		task = list(map(int, list(task)))
		if 'sgs_fofn' not in self.cfg:
			for i in [1, 2, 3, 4]:
				while i in task:
					task.remove(i)
					log.warning('Delete task: %d due to missing sgs_fofn.' % i)

		if 'lgs_fofn' not in self.cfg:
			for i in [3, 5]:
				while i in task:
					task.remove(i)
					log.warning('Delete task: %d due to missing lgs_fofn.' % i)

		if 'hifi_fofn' not in self.cfg:
			for i in [6]:
				while i in task:
					task.remove(i)
					log.warning('Delete task: %d due to missing hifi_fofn.' % i)

		for i in range(len(task)):
			if task[i] == 2 and task[i-1] != 1:
				log.error('Error, task 2 must follow task 1.')
				exit(1)
			elif task[i] == 3 and task[i-1] != 2:
				log.error('Error, task 3 must follow task 2.')
				exit(1)
			if task[i] == 4 and task[i-1] != 3:
				log.error('Error, task 4 must follow task 3.')
				exit(1)

		self.cfg['task'] = task
		log.info('scheduled tasks:\n' + str(self.cfg['task']))

	def _check(self):
		
		self.cfg['workdir'] = self.cfg['workdir'] if self.cfg['workdir'].startswith('/') else self.cfgdir + '/' + self.cfg['workdir']
		self.cfg['lgs_polish'] = self.cfg['workdir'] + '/%02d.lgs_polish'
		self.cfg['hifi_polish'] = self.cfg['workdir'] + '/%02d.hifi_polish'
		self.cfg['score_chain'] = self.cfg['workdir'] + '/%02d.score_chain'
		self.cfg['kmer_count'] = self.cfg['workdir'] + '/%02d.kmer_count'
		self.cfg['snp_phase'] = self.cfg['workdir'] + '/%02d.snp_phase'
		self.cfg['snp_valid'] = self.cfg['workdir'] + '/%02d.snp_valid'
		self.cfg['cleantmp'] = _bool(self.cfg['cleantmp'])
		self.cfg['rewrite'] = _bool(self.cfg['rewrite'])
		self.cfg['use_drmaa'] = _bool(self.cfg['use_drmaa'])
		
		# if self.cfg['job_type'] not in ['sge', 'local']:
		# 	log.error('Error, job_type only accept: sge|local')
		# 	sys.exit(1)

		if self.cfg['rewrite']:
			log.warning('Re-write workdir')
		if _bool(self.cfg['deltmp']):
			self.cfg['deltmp'] = True
		if _bool(self.cfg['rerun']):
			self.cfg['rerun'] = min(int(self.cfg['rerun']), 10)
		else:
			self.cfg['rerun'] = False

		if 'genome' not in self.cfg:
			log.error('Error, can not find genome option')
			sys.exit(1)
		else:
			self.cfg['genome'] = self.cfg['genome'] if self.cfg['genome'].startswith('/') else self.cfgdir + '/' + self.cfg['genome']
		if not os.path.exists(self.cfg['genome']):
			log.error('Error, can not find genome %s.' % self.cfg['genome'])
			sys.exit(1)
		if self.cfg['genome'].endswith('gz'):
			log.error('Error, samtools can not index gzip-ed files, please decompress %s.' % self.cfg['genome'])
			sys.exit(1)
		else:
			self.cfg['genome'] = self.cfg['genome'] if self.cfg['genome'].startswith('/') else self.cfgdir + '/' + self.cfg['genome']
			self.cfg['genome_size'] = calgs(self.cfg['genome']) if self.cfg['genome_size'] == 'auto' else parse_num_unit(self.cfg['genome_size'])
		
		# if 'min_len_seq' in self.cfg['polish_options']:
		# 	self.cfg['min_len_seq'] = parse_options_value(self.cfg['polish_options'], '-min_len_seq')
		# self.cfg['min_len_seq'] = parse_num_unit(self.cfg['min_len_seq'])

		if 'use_duplicate_reads' in self.cfg['sgs_options']:
			self.cfg['sgs_use_duplicate_reads'] = 1
		if 'unpaired' in self.cfg['sgs_options']:
			self.cfg['sgs_unpaired'] = '1'
		if '-N' in self.cfg['sgs_options']:
			self.cfg['sgs_rm_nread'] = 0

		if '-max_depth' in self.cfg['sgs_options']:
			self.cfg['sgs_max_depth'] = parse_options_value(self.cfg['sgs_options'], '-max_depth')
		if '-block_size' in self.cfg['sgs_options']:
			self.cfg['sgs_block_size'] = parse_options_value(self.cfg['sgs_options'], '-block_size')
		else:
			self.cfg['sgs_block_size'] = min(parse_num_unit(self.cfg['sgs_block_size']), self.cfg['genome_size'] * int(self.cfg['sgs_max_depth'])/int(self.cfg['parallel_jobs']))
		
		if 'sgs_fofn' not in self.cfg and 'lgs_fofn' not in self.cfg and 'hifi_fofn' not in self.cfg:
			log.error('Error, can not find sgs_fofn, lgs_fofn and hifi_fofn options')
			sys.exit(1)
		
		self.cfg['align_threads'] = self.cfg['multithread_jobs']

		if 'sgs_fofn' in self.cfg:
			self.cfg['sgs_fofn'] = self.cfg['sgs_fofn'] if self.cfg['sgs_fofn'].startswith('/') else self.cfgdir + '/' + self.cfg['sgs_fofn']
			if not os.path.exists(self.cfg['sgs_fofn']):
				log.error('Error, can not find sgs_fofn %s.' % self.cfg['sgs_fofn'])
				sys.exit(1)
			if 'minimap2' in self.cfg['sgs_options']:
				self.cfg['sgs_align_options'] = 'minimap2 --split-prefix tmp -a -x sr'
			else:
				self.cfg['sgs_align_options'] = 'bwa mem'
				if not _bool(self.cfg['sgs_unpaired']):
					self.cfg['sgs_align_options'] += ' -p '
			self.cfg['sgs_align_options'] += ' -t ' + self.cfg['multithread_jobs']

		if 'lgs_fofn' in self.cfg:
			self.cfg['lgs_fofn'] = self.cfg['lgs_fofn'] if self.cfg['lgs_fofn'].startswith('/') else self.cfgdir + '/' + self.cfg['lgs_fofn']
			if not os.path.exists(self.cfg['lgs_fofn']):
				log.error('Error, can not find lgs_fofn %s.' % self.cfg['lgs_fofn'])
				sys.exit(1)

			if 'min_read_len' in self.cfg['lgs_options']:
				self.cfg['lgs_min_read_len'] = parse_options_value(self.cfg['lgs_options'], '-min_read_len')
			if 'max_read_len' in self.cfg['lgs_options']:
				self.cfg['lgs_max_read_len'] = parse_options_value(self.cfg['lgs_options'], '-max_read_len')
			if 'max_depth' in self.cfg['lgs_options']:
				self.cfg['lgs_max_depth'] = parse_options_value(self.cfg['lgs_options'], '-max_depth')				
			if '-t' not in self.cfg['lgs_minimap2_options']:
				self.cfg['lgs_minimap2_options'] += ' -t ' + self.cfg['multithread_jobs']
			elif 'multithread_jobs' in self.cfg['lgs_minimap2_options']:
				self.cfg['lgs_minimap2_options'] = self.cfg['lgs_minimap2_options'].format(multithread_jobs = self.cfg['multithread_jobs'])
			if '-x' not in self.cfg['lgs_minimap2_options']:
				log.error('Error, failed find \'-x\' option in lgs_minimap2_options')
				sys.exit(1)
			elif 'map-ont' in self.cfg['lgs_minimap2_options']:
				self.cfg['lgs_read_type'] = 'ont'
			elif 'map-pb' in self.cfg['lgs_minimap2_options']:
				self.cfg['lgs_read_type'] = 'clr'
			else:
				log.error('Error, can not detect read_type for lgs_fofn.')
				sys.exit(1)

			if '-block_size' in self.cfg['lgs_options']:
				self.cfg['lgs_block_size'] = parse_options_value(self.cfg['lgs_options'], '-block_size')
			else:
				self.cfg['lgs_block_size'] = min(parse_num_unit(self.cfg['lgs_block_size']), self.cfg['genome_size'] * int(self.cfg['lgs_max_depth'])/int(self.cfg['parallel_jobs']))
			self.cfg['align_threads'] = max(int(self.cfg['align_threads']), int(parse_options_value(self.cfg['lgs_minimap2_options'], '-t')))
		
		if 'hifi_fofn' in self.cfg:
			self.cfg['hifi_fofn'] = self.cfg['hifi_fofn'] if self.cfg['hifi_fofn'].startswith('/') else self.cfgdir + '/' + self.cfg['hifi_fofn']
			if not os.path.exists(self.cfg['hifi_fofn']):
				log.error('Error, can not find hifi_fofn %s.' % self.cfg['hifi_fofn'])
				sys.exit(1)

			if 'min_read_len' in self.cfg['hifi_options']:
				self.cfg['hifi_min_read_len'] = parse_options_value(self.cfg['hifi_options'], '-min_read_len')
			if 'max_read_len' in self.cfg['hifi_options']:
				self.cfg['hifi_max_read_len'] = parse_options_value(self.cfg['hifi_options'], '-max_read_len')
			if 'max_depth' in self.cfg['hifi_options']:
				self.cfg['hifi_max_depth'] = parse_options_value(self.cfg['hifi_options'], '-max_depth')				
			if '-t' not in self.cfg['hifi_minimap2_options']:
				self.cfg['hifi_minimap2_options'] += ' -t ' + self.cfg['multithread_jobs']
			elif 'multithread_jobs' in self.cfg['hifi_minimap2_options']:
				self.cfg['hifi_minimap2_options'] = self.cfg['hifi_minimap2_options'].format(multithread_jobs = self.cfg['multithread_jobs'])
			if '-x' not in self.cfg['hifi_minimap2_options']:
				log.error('Error, failed find \'-x\' option in hifi_minimap2_options')
				sys.exit(1)

			if '-block_size' in self.cfg['hifi_options']:
				self.cfg['hifi_block_size'] = parse_options_value(self.cfg['hifi_options'], '-block_size')
			else:
				self.cfg['hifi_block_size'] = min(parse_num_unit(self.cfg['hifi_block_size']), self.cfg['genome_size'] * int(self.cfg['lgs_max_depth'])/int(self.cfg['parallel_jobs']))
			self.cfg['align_threads'] = max(int(self.cfg['align_threads']), int(parse_options_value(self.cfg['hifi_minimap2_options'], '-t')))

		if '-p' not in self.cfg['polish_options'] and '--process' not in self.cfg['polish_options']:
			self.cfg['polish_options'] += ' -p ' + self.cfg['multithread_jobs']
		elif 'multithread_jobs' in self.cfg['polish_options']:
			self.cfg['polish_options'] = self.cfg['polish_options'].format(multithread_jobs = self.cfg['multithread_jobs'])

		for opt in ['submit', 'kill', 'check_alive', 'job_id_regex']:
			if self.cfg[opt] and self.cfg[opt].lower() == 'auto':
				self.cfg[opt] = None

		del self.cfg['sgs_options']
		del self.cfg['lgs_options']
		del self.cfg['hifi_options']
