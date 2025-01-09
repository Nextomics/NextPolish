#!/usr/bin/env python
from __future__ import print_function

import os, re
import time
import signal
import logging

import sys
if sys.version_info[0] == 2:
	PYTHON_VERSION = 2
elif sys.version_info[0] == 3:
	PYTHON_VERSION = 3
else:
	raise Exception("Unknown Python version")


__all__ = ['plog', 'pmkdir', 'write2file', 'parse_num_unit', 'parse_options_value', 
	'calgs', 'readfile', 'getver', 'pypath', 'which', 'str2byte', 'byte2str', 'cal_n50_info',
	'remove_option', 'remove_options']

class TimedOutExc(Exception):
	pass

class deco(object):
	@staticmethod
	def _deadline(timeout):
		def _decorate(f):
			def handler(signum, frame):
				raise TimedOutExc
			def new_f(*args):
				signal.signal(signal.SIGALRM, handler)
				signal.alarm(timeout)
				try:
					return f(*args)
				except TimedOutExc:
					return 0
				finally:
					signal.alarm(0)
			return new_f
		return _decorate

class ExitOnCritical(logging.StreamHandler):
	def emit(self, record):
		if isinstance(record.msg, list):
			record.msg = ' '.join(record.msg)
		elif isinstance(record.msg, dict):
			record.msg ="\n" +  "\n".join(("%-30s%s" % (str(k).strip() + ':', str(v).strip()) for k, v in \
				sorted(record.msg.items(), key = lambda x: len(str(x[0]) + str(x[1])))))
		elif hasattr(record.msg, '__dict__'):
			record.msg ="\n" + "\n".join(("%-30s%s" % (str(k).strip() + ':', str(v).strip()) for k, v in \
				sorted(vars(record.msg).items(), key = lambda x: len(str(x[0]) + str(x[1])))))

		if record.levelno >= logging.ERROR:
			msg = record.msg
			record.msg = '\033[35m%s\033[0m' % msg
			super(ExitOnCritical, self).emit(record)
			record.msg = msg
		else:
			super(ExitOnCritical, self).emit(record)
		if record.levelno >= logging.CRITICAL:
			need_emit = False
			for handler in logging.getLogger().handlers:
				if need_emit:
					handler.emit(record)
				elif handler == self:
					need_emit = True
			raise Exception(record.msg)

def plog(path=None):
	formatter = logging.Formatter('[%(process)d %(levelname)s] %(asctime)s %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
	log = logging.getLogger()

	if not log.handlers:
		log.setLevel(logging.INFO)
		# console_handler = logging.StreamHandler()
		console_handler = ExitOnCritical()
		console_handler.setFormatter(formatter)
		log.addHandler(console_handler)

	if path:
		has_path_logger = False
		for logger in log.handlers:
			if isinstance(logger, logging.FileHandler) and logger.baseFilename == path:
				has_path_logger = True
				break

		if not has_path_logger:
			fileHandler = logging.FileHandler(path, mode='w')
			fileHandler.setFormatter(formatter)
			log.addHandler(fileHandler)
	return log

def pmkdir(path):
	if not os.path.exists(path):
		os.makedirs(path)
		return True
	else:
		return False

def ptime(init_time = 0):
	if init_time:
		return time.strftime("%H:%M:%S",time.gmtime((time.time() - init_time)))
	else:
		return time.strftime("%Y-%m-%d %H:%M:%S",time.gmtime((time.time())))

def write2file(content, path, mode = 'w'):
	with open(path, mode) as OUT:
		print(content, file=OUT)

def readfile(path):
	with open(path) as IN:
		return IN.read()

def parse_options_value(content, option, last = True, first = False, index = False):
	'''-p1 [4] -p2 [smp [5]] -p3 [ad] -d4 [f] '''
	contents = content.strip().split()
	contents_len = len(contents)
	for i in range(contents_len):
		if contents[i] == option:
			i += 1
			break
	j = i
	while (j <= contents_len):
		if j == contents_len or contents[j].startswith('-') or contents[j].startswith('>') or contents[j].startswith('1>'):
			if first:
				return i if index else contents[i]
			elif last:
				return j - 1 if index else contents[j - 1]
			else:
				return (i, j) if index else ' '.join(contents[i:j])
		j += 1

def remove_option(content, option):
	'''-p1 [4] -p2 [smp [5]] -p3 [ad] =>  -p1 [4] -p3 [ad]'''
	contents = content.strip().split()
	contents_len = len(contents)
	for i in range(contents_len):
		if contents[i] == option:
			break
	else:
		return content
	v = i + 1
	while v < contents_len and not contents[v].startswith('-'):
		v += 1
	return " ".join(contents[:i] + contents[v:])

def remove_options(content, options):
	for option in options:
		content = remove_option (content, options)
	return content

def parse_num_unit(content, base=1000):
	'''2Gb 2kb 2.3 kb 3.5g'''
	def parse_unit(unit, base=1000):
		if unit[0] in ['k', 'K']:
			return base
		elif unit[0] in ['m', 'M']:
			return base * base
		elif unit[0] in ['g', 'G']:
			return base * base * base

	if str(content)[-1].isdigit():
		return int(content)
	value, unit = 1, 1
	contents = str(content).strip().split()
	if len(contents) != 1:
		value = float(contents[0])
		unit = parse_unit(contents[1], base)
	else:
		if contents[0][-2].isdigit():
			value = float(contents[0][:-1])
			unit = parse_unit(contents[0][-1], base)
		else:
			value = float(contents[0][:-2])
			unit = parse_unit(contents[0][-2:], base)
	return int(value * unit + .499)

def calgs(infile):
	import ctypes
	calgs = ctypes.CDLL(os.path.dirname(os.path.realpath(__file__)) + '/' + "calgs.so")
	calgs.calgs.argtypes = [ctypes.c_char_p]
	calgs.calgs.restype = ctypes.c_uint64
	return calgs.calgs(infile.encode('utf-8'))

def latestver(url):
	if PYTHON_VERSION == 3:
		from urllib.request import urlopen
	else:
		from urllib2 import urlopen

	try:
		g = re.search(r'download/(.*?)/', urlopen(url, timeout=1).read().decode("utf-8"))
		if g:
			return g.group(1)
		else:
			return 'Unknown'
	except Exception as e:
		#raise e
		return 'Unknown'

def getver(path):
	ver = 'Unknown'
	readme = path + '/README.md'
	if os.path.exists(readme):
		with open(readme, 'r') as IN:
			g = re.search(r'download/(.*?)/', IN.read())
			if g:
				ver = g.group(1)
	latest = latestver('https://api.github.com/repos/Nextomics/NextPolish/releases/latest')
	if latest != 'Unknown' and ver != latest:
		print(('\033[35mPlease update to the latest version: %s, current version: %s \033[0m') % (latest, ver))
	return ver

def pypath():
	return sys.executable

def cal_n50_info(stat, outfile = None):
	stat.sort(key = int, reverse = True)
	gs = sum(stat)
	k = 1
	l = j = 0
	out = "%-5s %20s %20s\n" % ("Type", "Length (bp)", "Count (#)")
	for i in stat:
		l += i
		j += 1
		while l >= gs * 0.1 * k and k < 10:
			out += "N%d0 %20d%20d\n" % (k, i, j)
			k += 1
	out += '\n'
	out += "%-5s %18d%20s\n" % ("Min.", stat[-1], '-')
	out += "%-5s %18d%20s\n" % ("Max.", stat[0], '-')
	out += "%-5s %18d%20s\n" % ("Ave.", gs/j, '-')
	out += "%-5s %18d%20d\n" % ("Total", gs, j)
	if outfile:
		write2file(out, outfile)
	return out

def which(program):
	def is_exe(fpath):
		return os.path.isfile(fpath) and os.access(fpath, os.X_OK)
	fpath, fname = os.path.split(program)
	if fpath:
		if is_exe(program):
			return program
	else:
		for path in os.environ["PATH"].split(os.pathsep):
			exe_file = os.path.join(path, program)
			if is_exe(exe_file):
				return exe_file
	return None

def str2byte(string, ignore = True):
	if ignore:
		try:
			string = string.encode("UTF-8")
		except Exception:
			pass
		finally:
			return string
	else:
		return string.encode("UTF-8")

def byte2str(byte, ignore = True):
	if ignore:
		try:
			byte = byte.decode("UTF-8")
		except Exception:
			pass
		finally:
			return byte
	else:
		return byte.decode("UTF-8")
