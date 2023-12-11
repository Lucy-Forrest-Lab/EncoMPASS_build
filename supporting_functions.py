# Name: supporting_functions.py
# Language: python3
# Description: 
# Author: Edoardo Sarti
# Date: 2018/10/28

import re
import os
import gc
import io
import ast
import sys
import copy
import glob
import time
import json
import math
import gzip
import shutil
import codecs
import pickle
import random
import datetime
import pathlib
import requests
import string
import collections.abc
import subprocess
import numpy as np
import pandas as pd
import scipy.stats
import scipy.spatial.distance
import multiprocessing
import urllib.request
from Bio.PDB import PDBParser
from Bio.PDB import Superimposer
from Bio import Align
from Bio.Align import substitution_matrices
from functools import reduce
from pympler import muppy, summary
from selenium import webdriver
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.firefox.options import Options
from selenium.webdriver.common.by import By

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib import cm

# python3 -m pip install --index-url https://test.pypi.org/simple/ --no-deps locusts
import locusts.swarm
import uniprot_requests

import config  # Uses datetime

# SECTION 0. Handling error and warning messages ------------------------------------------------------
def print_log(message, thr_log_status='ERROR'):	#, log_filename=''):
	"""
	Function to print log messages.
	"message" is a 3-tuple (STATUS, FUNC_NAME, TEXT). 
	STATUS must be one of the following:
	  DEBUG for debugging purposes
	  NOTICE for tricky things happening under the hood
	  WARNING if something could create an error
	  ERROR if a process fails
	  CRITICAL if there is a contradiction and the program has to stop
	FUNC_NAME is the name of the function calling this function
	TEXT is a string of text, that can contain any kind of printable
	structure and any number of lines.
	"thr_log_status" determines which stati are also printed in stdout.
	If thr_log_status="ERROR", stati ERROR and CRITICAL are printed.
	"""

	if config.log_filename is None:
		if multiprocessing.current_process().name == 'MainProcess':
			config.log_filename = 'DEBUG_' + str(os.getpid()) + '_log.txt'
		else:
			config.log_filename = 'DEBUG_' + str(os.getppid()) + '_log.txt'
	
	log_stati = ['DEBUG', 'NOTICE', 'WARNING', 'ERROR', 'CRITICAL']
	
	# Checks over print_log itself and checks fixed arguments
	if thr_log_status not in log_stati:
		print_log(('CRITICAL', print_log.__name__, 'thr_log_status \'{0}\' is not recognized. List of available log_stati: {1}'.format(thr_log_status, log_stati)))
	if len(message) != 3:
		print_log(('CRITICAL', print_log.__name__, 'message argument should be a 3-tuple: (log_status, function_name, message_text). Here, message was {0}'.format(message)))
	
	# Defines variables
	log_status = message[0]
	function_logged = message[1]
	message_texts = message[2].split('\n')
	ts = time.time()
	file_time_stamp = datetime.datetime.fromtimestamp(ts).strftime('%Y%m%d_%H%M%S')
	time_stamp = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
	if len(function_logged) > 30:
		formatted_function_logged = function_logged[:27] + '...'
	else:
		formatted_function_logged = function_logged
	
	# Checks arguments 
	if log_status not in log_stati:
		print_log(('CRITICAL', print_log.__name__, 'message log_status argument \'{0}\' is not recognized. List of available log_stati: {1}'.format(log_status, log_stati)))
		exit(1)

	# Formats the output
	thr_ilog = log_stati.index(thr_log_status)
	msg_ilog = log_stati.index(log_status)
	formatted_info = '{0:20}\t{1:10}\t{2:30}\t'.format(time_stamp, log_status, formatted_function_logged)
	init = True
	if config.log_filename is None:
		config.log_filename = './log__' + function_logged + '__' + file_time_stamp + '.txt'
	log_file = open(config.log_filename, 'a')
	for message_text in message_texts:
		nc = 0
		while nc < len(message_text):
			jump = 0
			if len(message_text[nc:]) > 60:
				last_space_i = message_text[nc:nc+60].rfind(' ')
				if last_space_i == -1:
					new_nc = nc+60
				else:
					jump = 1
					new_nc = nc+last_space_i
			else:
				new_nc = nc + len(message_text[nc:])
			formatted_text = formatted_info + message_text[nc:new_nc].lstrip() + '\n'
			if init:
				formatted_info = ' '*20 + '\t' + ' '*10 + '\t' + ' '*30 + '\t'
				init = False
			
			log_file.write(formatted_text)
			if msg_ilog >= thr_ilog:
				print(formatted_text[:-1])      # It already has a '\n'
			nc = new_nc + jump

	if log_status == 'CRITICAL':
		exit(1)


# SECTION 1. Processing python data structures from strings -----------------------------------------------
def string_is_float(s, thr_log_status='ERROR'):
	this_name = string_is_float.__name__

	if not type(s) == str:
		raise ValueError
		stp_err = ('CRITICAL', this_name, "input type {0}, while str expected".format(type(s)))
		print_log(stp_err)
	try:
		float(s)
		return True
	except ValueError:
		return False

def string_is_int(s, thr_log_status='ERROR'):
	this_name = string_is_int.__name__

	if not type(s) == str:
		stp_err = ('CRITICAL', this_name, "input type {0}, while str expected".format(type(s)))
		print_log(stp_err)
	try:
		int(s)
		return True
	except ValueError:
		return False


def string_isnot_dict(element):
	inquotes = ""
	for c in element:
		if c == '"' or c == "'":
			if c == inquotes:
				inquotes = ""
			elif not inquotes:
				inquotes = c
		elif c == ":" and not inquotes:
			return False
	return True


def string_is_properstring(element):
	if len(element) < 2 or (not ((element[0] == "'" and element[-1] == "'" ) or (element[0] == '"' and element[-1] == '"'))):
		return False
	for nc, c in enumerate(element[1:-1]):
		if c == element[0] and (nc == 0 or element[nc-1] != "\\" ):
			return False
	return True


def string_decode_element(element, is_immutable=False, permissive=False, thr_log_status='ERROR'):
	"""WARNING: Cannot decode following types:
	decimal, complex, range, bytes, bytearrary,
	and any mutable or immutable user-defined class
	or non-standard library class (e.g. np.array())
	"""
	this_name = string_decode_element.__name__

	nop_err = ()
	element = element.strip()
	if string_is_properstring(element):
		if len(element) == 2:
			return ""
		else:
			return element[1:-1]
	elif element == "None":
		return None
	elif element == "False":
		return False
	elif element == "True":
		return True
	elif string_is_int(element):
		return int(element)
	elif string_is_float(element):
		return float(element)
	elif (len(element) > 1 and element[0] == "(" and element[-1] == ")") or (len(element) > 6 and element[:6] == "tuple(" and element[-1] == ")"):
		return string_decode_list(element, is_tuple=True)
	elif len(element) > 10 and element[:10] == "frozenset(" and element[-1] == ")":
		return string_decode_list(element, is_fset=True)
	elif not is_immutable:
		if (len(element) > 1 and element[0] == "{" and element[-1] == "}" and string_isnot_dict(element)) or (len(element) > 4 and element[:4] == "set(" and element[-1] == ")"):
			return string_decode_list(element, is_set=True)
		elif (len(element) > 1 and element[0] == "[" and element[-1] == "]") or (len(element) > 5 and element[:5] == "list(" and element[-1] == ")"):
			return string_decode_list(element)
		elif (len(element) > 1 and element[0] == "{" and element[-1] == "}") or (len(element) > 5 and element[:5] == "dict(" and element[-1] == ")"):
			return string_decode_dict(element)
		elif permissive:
			return element
		else:
			nop_err = ('CRITICAL', this_name, "could not process element {0}".format(element))
	elif permissive:
		return element
	else:
		nop_err = ('CRITICAL', this_name, "could not process immutable element {0}".format(element))

	if nop_err:
		print_log(nop_err)


def string_decode_list(s, is_tuple=False, is_set=False, is_fset=False, thr_log_status='ERROR'):
	this_name = string_decode_list.__name__

	if not type(s) == str:
		stp_err = ('CRITICAL', this_name, "input type {0}, while str expected".format(type(s)))
		print_log(stp_err)

	s = s.strip()
	is_cast = False
	if is_set:
		if (len(s) > 1 and s[0] == "{" and s[-1] == "}"):
			s = s[1:-1]
		elif (len(s) > 4 and s[:4] == "set(" and s[-1] == ")"):
			s = s[4:-1]
			is_cast = True
		else:
			fst_err = ('CRITICAL', this_name, "string does not appear to be a set: {0}".format(s))
	elif is_fset:
		if len(s) > 10 and s[:10] == "frozenset(" and s[-1] == ")":
			s = s[10:-1]
			is_cast = True
		else:
			fst_err = ('CRITICAL', this_name, "string does not appear to be a frozenset: {0}".format(s))
	elif is_tuple:
		if len(s) > 1 and s[0] == "(" and s[-1] == ")":
			s = s[1:-1]
		elif len(s) > 6 and s[:6] == "tuple(" and s[-1] == ")":
			s = s[6:-1]
			is_cast = True
		else:
			fst_err = ('CRITICAL', this_name, "string does not appear to be a tuple: {0}".format(s))
	else:
		if len(s) > 1 and s[0] == "[" and s[-1] == "]":
			s = s[1:-1]
		elif len(s) > 5 and s[:5] == "list(" and s[-1] == ")":
			s = s[5:-1]
			is_cast = True
		else:
			fst_err = ('CRITICAL', this_name, "string does not appear to be a list: {0}".format(s))

	if not s:
		if is_tuple:
			return ()
		elif is_set:
			return set()
		elif is_fset:
			return frozenset()
		else:
			return []

	d_list = []
	element = None
	inquotes = ""
	inbracks = ""
	for nc, c in enumerate(s):
		if c == "," and not inquotes and not inbracks:
			if type(element) == str:
				d_list.append(string_decode_element(element))
				element = None
		else:
			if c == '"' or c == "'":
				if c == inquotes:
					inquotes = ""
				elif not inquotes:
					inquotes = c
			elif not inquotes:
				if c == '{' or c == "[" or c == "(":
					inbracks += c
				elif c == '}':
					if inbracks and inbracks[-1] == '{':
						inbracks = inbracks[:-1]
				elif c == ']':
					if inbracks and inbracks[-1] == '[':
						inbracks = inbracks[:-1]
				elif c == ')':
					if inbracks and inbracks[-1] == '(':
						inbracks = inbracks[:-1]

			if type(element) == type(None):
				element = c
			elif type(element) == str:
				element += c
			
	if type(element) == str:
		d_list.append(string_decode_element(element))

	if is_cast:
		d_list = d_list[0]
	
	if is_tuple:
		return tuple(d_list)
	elif is_set:
		return set(d_list)
	elif is_fset:
		
		return frozenset(d_list)
	else:
		return d_list


def string_decode_dict(s, thr_log_status='ERROR'):
	this_name = string_decode_dict.__name__

	if not type(s) == str:
		stp_err = ('CRITICAL', this_name, "input type {0}, while str expected".format(type(s)))
		print_log(stp_err)
	s = s.strip()
	if not (len(s) > 1 and s[0] == "{" and s[-1] == "}"):
		fst_err = ('CRITICAL', this_name, "string does not appear to be a dict: {0}".format(s))
		print_log(fst_err)
	s = s[1:-1]
	if not s:
		return {}
	
	d_dict = {}
	key = None
	substr = None
	inquotes = ""
	inbracks = ""
	for c in s:
		if c == "," and not inquotes and not inbracks:
			if type(key) == type(None):
				bfe_err = ('CRITICAL', this_name, "no key found for dictionary element: {0}".format(key))
				print_log(bfe_err)

			d_dict[key] = string_decode_element(substr)
			key = None
			substr = None
		elif c == ":" and not inquotes and not inbracks:
			key = string_decode_element(substr, is_immutable=True)
			substr = None
		else:
			if c == '"' or c == "'":
				if c == inquotes:
					inquotes = ""
				elif not inquotes:
					inquotes = c
			elif not inquotes:
				if c == '{' or c == "[" or c == "(":
					inbracks += c
				elif c == '}':
					if inbracks and inbracks[-1] == '{':
						inbracks = inbracks[:-1]
				elif c == ']':
					if inbracks and inbracks[-1] == '[':
						inbracks = inbracks[:-1]
				elif c == ')':
					if inbracks and inbracks[-1] == '(':
						inbracks = inbracks[:-1]

			if type(substr) == type(None):
				substr = c
			elif type(substr) == str:
				substr += c
		
	if type(key) == type(None):
		bfe_err = ('CRITICAL', this_name, "no key found for dictionary element: {0}".format(key))
		print_log(bfe_err)

	d_dict[key] = string_decode_element(substr)

	return d_dict
			

def stringconversion_test():
	ans = {}

	s = "{'this' : 'is', 'a' : {'nested' : {'set' : '!'}}, 'with' : ['also', 1, 'list']}"
	ans[s] = {'this' : 'is', 'a' : {'nested' : {'set' : '!'}}, 'with' : ['also', 1, 'list']}

	s = "this should not be parsed"
	ans[s] = False

	s = "'False'"
	ans[s] = False

	s = "set([False, 'True', (1, 3.4, -1000, frozenset([1, '1.2', 'lol']))])"
	ans[s] = set([False, 'True', (1, 3.4, -1000, frozenset([1, '1.2', 'lol']))])

	s = "[True, True, True, 'False, True', False]"
	ans[s] = [True, True, True, 'False, True', False]

	s = "{'a : 'b'}"
	ans[s] = False

	s = "[[[1]]]"
	ans[s] = [[[1]]]

	s = "[{2,3]]"
	ans[s] = False

	for k in ans:
		try:
			e = string_decode_element(k)
		except:
			e = False

def is_number(s):
	try:	
		float(str(s))
		return True
	except ValueError:
		return False


class FixedValDict(collections.abc.MutableMapping):
	"""
	Dictionary collection of Fixed Dictionaries (that must be all equal)
	In a Fixed Value Dictionary, keys are arbitrary, whereas values follow a Fixed Dictionary template
	1. Initialization
		If FixedValDict is initialized with a FixedDict, it will interpret it as its template,
		and it will return an empty structure. If it is initialized with a dictionary of FixedDict
		it will interpret each FixedDict as one of its values. All FixedDict must be compatible
		(i.e. all entries must be subsets of the most complete one), and the structure of the largest
		FixedDict is taken as template.
	2. Adding an entry
		To add an entry to an existing FixedValDict, the new(k) method must be called. k must contain
		the key to add. Once the method has been called, the new FixedDict can be filled
	"""
	def __init__(self, dictionary):
		self._dictionary = {}
		if isinstance(dictionary, FixedDict):	# It's the template. Must be FD
			self.elemdictionary = FixedDict(dictionary.get_fdict_template())	# It might be a full FD: you have to nullify
		elif type(dictionary) == dict:		# It's a full FVD. All valuse must be FD
			self.elemdictionary = FixedDict({})
			for x in dictionary:
				if not isinstance(dictionary[x], FixedDict):
					raise TypeError
				if self.elemdictionary.show_dictionary(quiet=True, nullify=True, deeptransparency=True, annotate=True) != dictionary[x].show_dictionary(quiet=True, nullify=True, deeptransparency=True, annotate=True):
					is_contained = True
					for y in self.elemdictionary:
						if y not in dictionary[x]:
							is_contained = False
					if is_contained:
						self.elemdictionary = FixedDict(dictionary[x].show_dictionary(quiet=True, nullify=True))
					else:
						raise TypeError
				self._dictionary[x] = dictionary[x]
			for x in self._dictionary:
				for y in self.elemdictionary:
					if y not in self._dictionary[x]:
						self._dictionary[x][y] = self.elemdictionary[y]
		else:
			raise TypeError
	def __setitem__(self, key, item):
		a = item.show_dictionary(quiet=True, nullify=True, deeptransparency=True, annotate=True)
		b = self.elemdictionary.show_dictionary(quiet=True, nullify=True, deeptransparency=True, annotate=True)
		if not isinstance(item, FixedDict) or a != b:
			print(">>>>>>>>>>>>>>>>> HERE2 <<<<<<<<<<<<<<<<<")
			raise TypeError
		self._dictionary[key] = item
	def __getitem__(self, key):
		return self._dictionary[key]
	def __delitem__(self, key):
		self._dictionary.__delitem__(key)
	def __contains__(self, key):
		return key in self._dictionary
	def __len__(self):
		return len(self._dictionary)
	def __iter__(self):
		return iter(self._dictionary)
	def __eq__(self, other):
		if isinstance(other, FixedValDict):
			for k in self._dictionary:
				if not (k in other._dictionary and self._dictionary[k] == other._dictionary[k]):
					return False
			for k in other._dictionary:
				if not (k in self._dictionary and self._dictionary[k] == other._dictionary[k]):
					return False
			return True
		else:
			return False
	def __steq__(self, other, selfdict={}):
		# Equality of structures, but not of contents
		return self.show_elemdictionary(quiet=True, nullify=True, deeptransparency=True, annotate=True) == other.show_elemdictionary(quiet=True, nullify=True, deeptransparency=True, annotate=True)
	def show_elemdictionary(self, t='', selfcontainer=None, quiet=False, transparency=False, deeptransparency=False, nullify=False, annotate=False):
		self.elemdictionary.show_dictionary(t=t, selfcontainer=selfcontainer, quiet=quiet, transparency=transparency, deeptransparency=deeptransparency, nullify=nullify, annotate=annotate)
	def show_dictionary(self, t='', selfcontainer=None, quiet=False, transparency=False, deeptransparency=False, nullify=False, annotate=False):
		sc = {}
		if len(self._dictionary) == 0 or nullify:
			if not quiet:
				print(t+'template_entry')
			if deeptransparency:
				self.elemdictionary.show_dictionary(t=t+'\t', quiet=quiet, nullify=True, transparency=transparency, deeptransparency=deeptransparency, annotate=annotate)
			elif transparency:
				self.elemdictionary.show_dictionary(t=t+'\t', quiet=quiet, nullify=True, transparency=False, deeptransparency=deeptransparency, annotate=annotate)
			else:
				sc = FixedValDict(FixedDict(self.elemdictionary.show_dictionary(t=t+'\t', quiet=quiet, nullify=True, transparency=False, deeptransparency=False, annotate=annotate)))
		else:
			for k in self._dictionary:
				if not quiet:
					print(t+str(k), "FixedDict")
				if annotate:
					y = "FixedDict::" + k
				else:
					y = k
				sc[y] = self._dictionary[k].show_dictionary(t=t+'\t', quiet=quiet, transparency=transparency, deeptransparency=deeptransparency, nullify=False, annotate=annotate)
		return sc
	def get_fdict_template(self, t='', quiet=True):
		return self.elemdictionary.show_dictionary(t=t, quiet=quiet, nullify=True)
	def new(self, key):
		self._dictionary[key] = FixedDict(self.elemdictionary)
	

class FixedDict(collections.abc.MutableMapping):
	def __init__(self, dictionary):
		if isinstance(dictionary, FixedDict):
			self._dictionary = dictionary.show_dictionary(quiet=True)
		elif type(dictionary) == dict:
			self._dictionary = dictionary
		else:
			print(">>>>>>>>>>>>>>>>> HERE3 <<<<<<<<<<<<<<<<<")
			raise TypeError
	def __setitem__(self, key, item):
		if key not in self._dictionary:
			raise KeyError("The key {} is not defined.".format(key))
		self._dictionary[key] = item
	def __getitem__(self, key):
		return self._dictionary[key]
	def __delitem__(self, key):
		raise NotImplementedError
	def __contains__(self, key):
		return key in self._dictionary
	def __len__(self):
		return len(self._dictionary)
	def __iter__(self):
		return iter(self._dictionary)
	def __eq__(self, other):
		if isinstance(other, FixedDict):
			for x in self._dictionary:
				if not (x in other._dictionary and self._dictionary[x] == other._dictionary[x]):
					return False
			for x in other._dictionary:
				if not (x in self._dictionary and self._dictionary[x] == other._dictionary[x]):
					return False
			return True
		else:
			return False
	def __steq__(self, other, selfdict={}):
		# Equality of structures, but not of contents
		od = other.show_dictionary(quiet=True, nullify=True, deeptransparency=True, annotate=True)
		sd = self.show_dictionary(quiet=True, nullify=True, deeptransparency=True, annotate=True)
		return od == sd
	def show_dictionary(self, t='', selfcontainer=None, quiet=False, transparency=False, deeptransparency=False, nullify=False, annotate=False):
		contlist = [dict, list, set, frozenset]
		emptycontlist = [{}, [], set(), frozenset()]
		if selfcontainer == None:
			con = self._dictionary
		else:
			con = selfcontainer
		if type(con) == dict:
			sc = {}
			if not con and not quiet:
				print(t+"{}")
			for x in con:
				if type(con[x]) in contlist:
					if not quiet:
						print(t+str(x), str(type(con[x]).__name__))
					if not nullify:
						fd = FixedDict({})
						sc[x] = fd.show_dictionary(t=t+'\t', selfcontainer=con[x], quiet=quiet, transparency=False, deeptransparency=deeptransparency, nullify=False, annotate=annotate)
					else:
						sc[x] = emptycontlist[contlist.index(type(con[x]))]
				elif isinstance(con[x], FixedList):
					if not quiet:
						print(t+str(x), "FixedList")
					if annotate:
						y = "FixedList::" + x
					else:
						y = x
					sc[y] = con[x].show_list(t=t+'\t', quiet=quiet, deeptransparency=deeptransparency, nullify=nullify, annotate=annotate)
				elif isinstance(con[x], FixedDict):
					if not quiet:
						print(t+str(x), "FixedDict")
					if annotate:
						y = "FixedDict::" + x
					else:
						y = x
					sc[y] = con[x].show_dictionary(t=t+'\t', quiet=quiet, transparency=False, deeptransparency=deeptransparency, nullify=nullify, annotate=annotate)
				elif isinstance(con[x], FixedValDict):
					if not quiet:
						print(t+str(x), "FixedValDict")
					if annotate:
						y = "FixedValDict::" + x
					else:
						y = x
					sc[y] = con[x].show_dictionary(t=t+'\t', quiet=quiet, transparency=False, deeptransparency=deeptransparency, nullify=nullify, annotate=annotate)
				else:
					if nullify:
						sc[x] = ''
					else:
						sc[x] = con[x]
					if not quiet:
						s = con[x]
						a = ''
						if (s == '') or (not is_number(s)):
							a = '\''
						print(t+str(x), a+str(s)+a)
			if (selfcontainer == None) and (not ((not t) or transparency or deeptransparency)):
				return FixedDict(sc)
			else:
				return sc
		elif type(con) == list:
			sc = []
			if nullify or not con:
				if not quiet:
					print('{0}{1}'.format(t,sc))
				return sc
			for x in con:
				if type(x) in contlist:
					if not quiet:
						print('{0}{1}'.format(t,str(type(x).__name__)))
					sc.append(self.show_dictionary(t=t+'\t', selfcontainer=x, quiet=quiet, transparency=False, deeptransparency=deeptransparency, nullify=nullify, annotate=annotate))
				elif isinstance(x, FixedList):
					if not quiet:
						print('{0}{1}'.format(t,"FixedList"))
					sc.append(x.show_list(t=t+'\t', quiet=quiet, deeptransparency=deeptransparency, nullify=nullify, annotate=annotate))
				elif isinstance(x, FixedDict):
					if not quiet:
						print('{0}{1}'.format(t,"FixedDict"))
					sc.append(x.show_dictionary(t=t+'\t', quiet=quiet, transparency=False, deeptransparency=deeptransparency, nullify=nullify, annotate=annotate))
				else:
					sc.append(x)
					if not quiet:
						a = ''
						if (x == '')  or not is_number(x):
							a = '\''
						print('{0}{1}'.format(t,a+str(x)+a))
			return sc
		elif type(con) == set:
			sc = set()
			if nullify or not con:
				if not quiet:
					print('{0}{1}'.format(t,sc))
				return sc
			for x in con:
				if type(x) == frozenset:
					if not quiet:
						print('{0}{1}'.format(t,str(type(x).__name__)))
					sc.add(self.show_dictionary(t=t+'\t', selfcontainer=x, quiet=quiet, transparency=False, deeptransparency=deeptransparency, nullify=nullify, annotate=annotate))
				elif isinstance(x, FixedDict) or isinstance(x, FixedList) or type(x) in contlist:
					print(">>>>>>>>>>>>>>>>> HERE4 <<<<<<<<<<<<<<<<<")
					raise TypeError
				else:
					sc.add(x)
					if not quiet:
						a = ''
						if not is_number(x):
							a = '\''
						print('{0}{1}'.format(t,a+str(x)+a))
			return sc
		elif type(con) == frozenset:
			sc = set()
			if nullify or not con:
				if not quiet:
					print('{0}{1}'.format(t,frozenset(sc)))
				return frozenset(sc)
			for x in con:
				if type(x) == frozenset:
					if not quiet:
						print('{0}{1}'.format(t,str(type(x).__name__)))
					sc.add(self.show_dictionary(t=t+'\t', selfcontainer=x, quiet=quiet, transparency=False, deeptransparency=deeptransparency, nullify=nullify, annotate=annotate))
				elif isinstance(x, FixedDict) or isinstance(x, FixedList) or type(x) in contlist:
					print(">>>>>>>>>>>>>>>>> HERE5 <<<<<<<<<<<<<<<<<")
					raise TypeError
				else:
					sc.add(x)
					if not quiet:
						a = ''
						if not is_number(x):
							a = '\''
						print('{0}{1}'.format(t,a+str(x)+a))
			return frozenset(sc)
	def get_fdict_template(self, t='', quiet=True):
		return self.show_dictionary(t=t, quiet=quiet, nullify=True)

def FDtemplate_test(d_templ, is_fd=False):
    print("\n### FDtemplate_test log ###")
    if not is_fd:
        print("\n> Template dictionary")
        print(d_templ)
        print("\n> FD construct")
        fd_templ = FixedDict(d_templ)
    else:
        fd_templ = d_templ
    print("\n> FD.show_dictionary(transparency=True)")
    print(fd_templ.show_dictionary(transparency=True))
    print("\n> FD.show_dictionary(deeptransparency=True)")
    print(fd_templ.show_dictionary(deeptransparency=True))
    print("\n> FD.show_dictionary(deeptransparency=True, annotate=True)")
    print(fd_templ.show_dictionary(deeptransparency=True, annotate=True))
    print("\n> FD.get_fdict_template()")
    print(fd_templ.get_fdict_template())
    print("\n### FDtemplate_test log end ###")
    return fd_templ


def FD_test():
    print("1. Create a template FD")
    print("    1.1 From an empty dictionary")
    d_template = {"key_1" : "", "key_2" : [], "key_3" : {}, "key_4" : set(), "key_5" : frozenset()}
    FDtemplate_test(d_template)
    print("\n    1.2 From an empty nested dictionary")
    d_template = {"key_1" : "", "key_2" : [], "key_3" : FixedDict({"key_31" : "", "key_32" : [], "key_33" : {}, "key_34" : set(), "key_35" : frozenset()}), "key_4" : set(), "key_5" : frozenset()}
    FDtemplate_test(d_template)
    print("\n    1.3 From a full dictionary")
    d_template = {"key_1" : "oh, mama!", "key_2" : [1,2,3], "key_3" : {1:2, 2:3, 3:1}, "key_4" : set({1,2,3}), "key_5" : frozenset({1,2,3})}
    FDtemplate_test(d_template)
    print("\n    1.2 From a full nested dictionary")
    d_template = {"key_1" : "oh, mama!", "key_2" : [1,2,3], "key_3" : FixedDict({"key_31" : "oh, papa!", "key_32" : [4,5,6], "key_33" : {4:5, 5:6, 6:4}, "key_34" : set({4,5,6}), "key_35" : frozenset({4,5,6})}), "key_4" : set({4,5,6}), "key_5" : frozenset({4,5,6})}
    FDtemplate_test(d_template)

    print("2. Read/Write a template FD")
    print("   Note: requires json package")
    print("   Note: set() and frozenset() not implemented")
    print("   WARNING: assumes ALL dictionary keys are strings!")
    print("    2.1 Write a full nested TD in FD_test_1.txt")
    d_template = {"key_1" : "oh, mama!", "key_2" : [1,2,3], "key_3" : FixedDict({"key_31" : "oh, papa!", "key_32" : [4,5,6], "key_33" : {'4':5, '5':6, '6':4}})}
    fd_templ = FDtemplate_test(d_template)
    text = json.dumps(fd_templ.show_dictionary(quiet=True, deeptransparency=True, annotate=True))
    with open("FD_test_1.txt", 'w') as wf:
        print("\n> Transcript")
        print(text)
        wf.write(text)
    fd_templ_transp = json.loads(open("FD_test_1.txt", 'r').read(), object_hook=object_decoder_transp)
    print("\n> Loaded transparent")
    print(fd_templ_transp)
    print("\n> Loaded")
    fd_templ_afterread = FixedDict(json.loads(open("FD_test_1.txt", 'r').read(), object_hook=object_decoder))
    FDtemplate_test(fd_templ_afterread, is_fd=True)
    print("\n> Check: is it the same as before the write/read cycle?")
    print(fd_templ == fd_templ_afterread)

    print("3. Fit a normal data structure to a template FD")
    print("   Note: requires json package")
    print("   Note: set() and frozenset() not implemented")
    print("   WARNING: assumes ALL dictionary keys are strings!")
    print("    3.1 fit an incomplete data structure to an empty template FD")
    d_template = {"key_1" : "", "key_2" : [], "key_3" : FixedDict({"key_31" : "", "key_32" : [], "key_33" : {}})}
    fd_templ = FDtemplate_test(d_template)
    print("\t> Data structure")
    data_structure = {"key_1" : "hello, there!", "key_3" : {"key_33" : {1:2, 2:3, 3:4}}}
    fitted_fd = fit_into_template(fd_templ, data_structure)
    FDtemplate_test(fitted_fd, is_fd=True)


def FVDtemplate_test(fd_templ, is_fvd=False):
    print("\n### FVDtemplate_test log ###")
    if not is_fvd:
        print("\n> Template FD.show_dictionary(deeptransparency=True, annotate=True)")
        print(fd_templ.show_dictionary(deeptransparency=True, annotate=True))
        print("\n> FVD construct")
        fvd_templ = FixedValDict(fd_templ)
    else:
        fvd_templ = FixedValDict(fd_templ)
    print("\n> FVD.show_dictionary(transparency=True)")
    print(fvd_templ.show_dictionary(transparency=True))
    print("\n> FVD.show_dictionary(deeptransparency=True)")
    print(fvd_templ.show_dictionary(deeptransparency=True))
    print("\n> FVD.show_dictionary(deeptransparency=True, annotate=True)")
    print(fvd_templ.show_dictionary(deeptransparency=True, annotate=True))
    print("\n> FVD.get_fdict_template()")
    print(fvd_templ.get_fdict_template())
    print("\n### FVDtemplate_test log end ###")
    return fvd_templ


def FVD_test():
    print("1. Create a template FVD")
    print("    1.1 From an empty dictionary")
    fd_template = FixedDict({"key_1" : "", "key_2" : [], "key_3" : {}, "key_4" : set(), "key_5" : frozenset()})
    FVDtemplate_test(fd_template)
    print("\n    1.2 From a full multiple FD")
    fd_template = {"entry_1" : FixedDict({"key_1" : "something", "key_2" : [1,2,3], "key_3" : {1:2, 2:3, 3:1}, "key_4" : set({1,2,3}), "key_5" : frozenset({1,2,3})}), "entry_2" : FixedDict({"key_1" : "something else", "key_2" : [4,5,6], "key_3" : {4:5, 5:6, 6:4}, "key_4" : set({4,5,6}), "key_5" : frozenset({4,5,6})})}
    FVDtemplate_test(fd_template, is_fvd=True)
    print("\n    1.3 From an empty nested dictionary")
    fd_template = FixedDict({"key_1" : "", "key_2" : [], "key_3" : FixedValDict(FixedDict({"key_31" : "", "key_32" : [], "key_33" : {}, "key_34" : set(), "key_35" : frozenset(), "key_36" : FixedValDict(FixedDict({"key_361" : ""}))})), "key_4" : set(), "key_5" : frozenset()})
    print("uno")
    uno = FixedValDict(FixedDict({"key_361" : ""}))
    print("duominu")
    duominu = FixedValDict(FixedDict({"key_31" : "", "key_32" : [], "key_33" : {}, "key_34" : set(), "key_35" : frozenset(), "key_36" : FixedDict({"key_361" : ""})}))
    duominu.new("fvd_key_1")
    duominu["fvd_key_1"]["key_36"]["key_361"] = "ciao"
    print("duo")
    duo = FixedDict({"key_31" : "", "key_32" : [], "key_33" : {}, "key_34" : set(), "key_35" : frozenset(), "key_36" : uno})
    duo["key_36"].new("fvd_key_11")
    duo["key_36"]["fvd_key_11"]["key_361"] = "ciao"
    print("SHOW DUO", duo.show_dictionary(nullify=True))
    truo = FixedValDict(duo)
    fd_template = FixedDict({"key_1" : "", "key_2" : [], "key_3" : truo, "key_4" : set(), "key_5" : frozenset()})
    fd_template.show_dictionary(deeptransparency=True)
    fd_template["key_3"].new("fvd_key_1")
    fd_template["key_3"]["fvd_key_1"]["key_36"].new("fvd_key_11")
    FVDtemplate_test(fd_template)

    print("\n    1.3 From a full dictionary")
    fd_template = FixedDict({"key_1" : "oh, mama!", "key_2" : [1,2,3], "key_3" : {1:2, 2:3, 3:1}, "key_4" : set({1,2,3}), "key_5" : frozenset({1,2,3})})
    FVDtemplate_test(fd_template)
    print("\n    1.4 From a full nested dictionary")
    fd_template = FixedDict({"key_1" : "oh, mama!", "key_2" : [1,2,3], "key_3" : FixedDict({"key_31" : "oh, papa!", "key_32" : [4,5,6], "key_33" : {4:5, 5:6, 6:4}, "key_34" : set({4,5,6}), "key_35" : frozenset({4,5,6})}), "key_4" : set({4,5,6}), "key_5" : frozenset({4,5,6})})
    FVDtemplate_test(fd_template)
    print("\n    1.3 From a full FVD")


# A list of fixed dictionaries
class FixedList(collections.abc.MutableSequence):
	def __init__(self, fdict):
		if isinstance(fdict, FixedDict):
			self._fixed_dict = FixedDict(fdict.show_dictionary(quiet=True, nullify=True))
			self._list = list()
		elif isinstance(fdict, FixedList):
			self._fixed_dict = FixedDict(fdict.show_dictionary(quiet=True, nullify=True))
			self._list = list()
			for x in fdict:
				self._list.append(FixedDict(x.show_dictionary(quiet=True)))
		elif type(fdict) == list:
			if len(fdict) == 0:
				print("ERROR: cannot initialize FixedList with an empty list")
				exit(1)
			elif isinstance(fdict[0], FixedDict):
				self._fixed_dict = FixedDict(fdict[0])
				self._list = list()
				for x in fdict:
					if isinstance(x, FixedDict) and x.__steq__(self._fixed_dict):
						self._list.append(FixedDict(x.show_dictionary(quiet=True)))
			else:
				print("ERROR: initalization list must contain FixedDict elements")
				exit(1)
		else:
			print(">>>>>>>>>>>>>>>>> HERE6 <<<<<<<<<<<<<<<<<")
			raise TypeError
	def __len__(self):
		return len(self._list)
	def __delitem__(self, index):
		self._list.__delitem__(index)
	def insert(self, index, value):
		if isinstance(value, FixedDict) and value.__steq__(self._fixed_dict):
			self._list.insert(index, FixedDict(value.show_dictionary(quiet=True)))
		else:
			print(">>>>>>>>>>>>>>>>> HERE7 <<<<<<<<<<<<<<<<<")
			raise TypeError
	def __setitem__(self, index, value):
		if isinstance(value, FixedDict) and value.__steq__(self._fixed_dict):
			self._list.__setitem__(index, FixedDict(value.show_dictionary(quiet=True)))
		else:
			print(">>>>>>>>>>>>>>>>> HERE8 <<<<<<<<<<<<<<<<<")
			raise TypeError
	def __getitem__(self, index):
		return self._list.__getitem__(index)
	def append(self, value):
#		print(value, value.show_dictionary(), isinstance(value, FixedDict), value.__steq__(self._fixed_dict))
		if isinstance(value, FixedDict) and value.__steq__(self._fixed_dict):
			self._list.append(FixedDict(value.show_dictionary(quiet=True)))
		else:
			print(">>>>>>>>>>>>>>>>> HERE9 <<<<<<<<<<<<<<<<<")
			raise TypeError
	def get_fdict_template(self, t='', quiet=True):
		return self._fixed_dict.show_dictionary(t=t, quiet=quiet, nullify=True)
	def show_list(self, t='', quiet=False, deeptransparency=False, nullify=False, annotate=False):
		if not t:
			transparency = True
		else:
			transparency = False
		if (not t) or deeptransparency:
			l = []
		else:
			l = FixedList(self.get_fdict_template(t=t))
		if len(self._list) == 0 or nullify:
			if not quiet:
				print(t+'template_entry')
				self._fixed_dict.show_dictionary(t=t+'\t', quiet=quiet, transparency=transparency, deeptransparency=deeptransparency, nullify=True, annotate=annotate)
		else:
			for i,x in enumerate(self._list):
				if not quiet:
					print(t+'entry_'+str(i))
				li = x.show_dictionary(t=t+'\t', quiet=quiet, transparency=transparency, deeptransparency=deeptransparency, nullify=nullify, annotate=annotate)
				l.append(li)
		return l
	show_dictionary = show_list


### NEW REFACTORING
def passport_entry(code_location, pdbi, message, short_message=""):
    tmp_str_data = define_str_data_entry(config.str_data_template_fn)['PASSPORT'].get_fdict_template()
    tmp_str_data['location'] = code_location
    tmp_str_data['date'] = datetime.datetime.now().strftime("%Y/%m/%d %H:%M:%S")
    tmp_str_data['short_message'] = short_message
    tmp_str_data['message'] = message
    print_log((
        'NOTICE', 
        code_location, 
        "[passport {0}] ".format(pdbi) + message))
    return FixedDict(tmp_str_data)
### END NEW REFACTORING


def object_decoder(obj, deeptransparency=False):
	new_obj = {}
	for x in obj:
		if 'FixedDict::' in x:
			new_key = ''.join(x.split('::')[1:])
			if deeptransparency:
				new_obj[new_key] = obj[x]
			else:
#				print(x, obj[x])
				new_obj[new_key] = FixedDict(obj[x])
		elif 'FixedList::' in x:
			new_key = ''.join(x.split('::')[1:])
			if deeptransparency:
				new_obj[new_key] = obj[x]
			else:
				new_obj[new_key] = FixedList(FixedDict(obj[x][0]))
		elif 'FixedValDict::' in x:
			new_key = ''.join(x.split('::')[1:])
			if deeptransparency:
				new_obj[new_key] = obj[x]
			else:
				new_obj[new_key] = FixedValDict(FixedDict(obj[x]))
		else:
			new_obj[x] = obj[x]
	return new_obj


def object_decoder_transp(obj):
	return object_decoder(obj, deeptransparency=True)


def define_str_data_entry(str_data_entry_fn):
    with open(str_data_entry_fn) as str_data_entry_f:
        str_data_entry = FixedDict(json.load(str_data_entry_f, 
            object_hook=object_decoder))
    return str_data_entry

def write_checkpoint(str_data, logfilename):
	transp_str_data = {}
	for pdbi in str_data:
		transp_str_data[pdbi] = str_data[pdbi].show_dictionary(quiet=True, deeptransparency=True, annotate=True)
	transp_str_data_str = json.dumps(transp_str_data)
	pickle.dump(transp_str_data_str, open(logfilename, 'wb'))


def fit_into_template(template_str_data_entry, data_structure, recur=0):
	if isinstance(template_str_data_entry, FixedDict):
		if type(data_structure) != dict:
			raise TypeError
		if data_structure == {}:
			return FixedDict(template_str_data_entry.show_dictionary(quiet=True, nullify=True))
		else:
			nulltemp = FixedDict(template_str_data_entry.show_dictionary(quiet=True, nullify=True))
			found_in_data_structure = set()
			for x in template_str_data_entry:
				if x in data_structure:
					nulltemp[x] = fit_into_template(template_str_data_entry[x], data_structure[x], recur=recur+1)
					found_in_data_structure.add(x)
				else:
					nulltemp[x] = template_str_data_entry[x]
			if {x for x in data_structure} != found_in_data_structure:
				print("ERROR: ARGUMENTS IN DATA STRUCTURE NOT IN TEMPLATE:", {x for x in data_structure} - found_in_data_structure)
				raise TypeError
			return nulltemp
	if isinstance(template_str_data_entry, FixedValDict):
		if type(data_structure) != dict:
			raise TypeError
		if data_structure == {}:
			return FixedValDict(FixedDict(template_str_data_entry.get_fdict_template()))
		else:
			nulltemp = FixedValDict(FixedDict(template_str_data_entry.get_fdict_template()))
			for x in data_structure:
				nulltemp.new(x)
				nulltemp[x] = fit_into_template(FixedDict(template_str_data_entry.get_fdict_template()), data_structure[x], recur=recur+1)
			return nulltemp
	elif isinstance(template_str_data_entry, FixedList):
		if type(data_structure) != list:
			raise TypeError
		if data_structure == []:
			return FixedList(FixedDict(template_str_data_entry.get_fdict_template()))
		else:
			nulltemp = FixedList(FixedDict(template_str_data_entry.get_fdict_template()))
			for x in data_structure:
				nulltemp.append(fit_into_template(FixedDict(template_str_data_entry.get_fdict_template()), x, recur=recur+1))
			return nulltemp
	else:
		if type(data_structure) != type(template_str_data_entry):
			if type(template_str_data_entry) == str:
				return str(data_structure)
			elif type(template_str_data_entry) == int:
				return int(data_structure)
			elif type(template_str_data_entry) == float:
				return float(data_structure)
			elif type(template_str_data_entry) == bool:
				return bool(data_structure)
			else:
				print("TYPE ERROR: ", type(data_structure), type(template_str_data_entry))
				raise TypeError
		else:
			return data_structure


def write_locusts_parfile(options, locusts_filename, job_code, only_gather=False):
    text =  (
        "### Generic\n"
        "{0:30}\t{1}\n\n"
    ).format('run_on_hpc', options['RUN'][('hpc', 'run_on_hpc')])

    if options['RUN'][('hpc', 'run_on_hpc')]:
        text += "### HPC\n"
       
        for k in ['host_name', 'partition', 'requested_nodes', 'cpus_per_node', 'local_shared_dir', 'data_transfer_protocol', 'nodewise_scratch_folder', 'nodewise_scratch_memory', 'email_address', 'extra_outer_statements', 'min_stack_per_core', 'walltime']:
            if options['ALL'][k]:
                text += "{0:30}\t{1}\n".format(k, options['ALL'][k])
        text += "{0:30}\t{1}\n".format('hpc_exec_dir', options['ALL']['hpc_exec_path'] + '/' + job_code)  # EncoMPASS hpc_exec_path is the parent of each locusts hpc_exec_dir
        text += "min_stack_per_core\t1\n"
        text += "\n### Singularity\n"
        for k in ['singularity', 'singularity_container', 'singularity_modload']:
            text += "{0:30}\t{1}\n".format(k, options['ALL']['hpc_' + k])
        text += "exclusive\tTrue\n"
    else:
        text += (
            "### Local multithreading\n"
            "{0:30}\t{1}\n\n"
        ).format('number_of_processors', options['ALL']['number_of_processors'])
        text += "\n### Singularity\n"
        for k in ['singularity', 'singularity_container', 'singularity_modload']:
            text += "{0:30}\t{1}\n".format(k, options['ALL'][k])

    if only_gather:
        text += "only_gather\tTrue\n"

    print("LOCUSTS FILE")
    print(text)

    with open(locusts_filename, 'w') as outf:
        outf.write(text)

    return

def read_checkpoint(logfilename, template_filename):
	template_str_data_entry = define_str_data_entry(template_filename)
	transp_str_data_str = pickle.load(open(logfilename, 'rb'))
	tmp_data = json.loads(transp_str_data_str, object_hook=object_decoder_transp)
	str_data = {}
	for pdbi in tmp_data:
		str_data[pdbi] = fit_into_template(template_str_data_entry, tmp_data[pdbi])
	return str_data


def listest():
	l = FixedList(FixedDict({'fl' : FixedList(FixedDict({1 : '', 2 : ''})), 2 : ''}))
	ltry = [{}, {'fl' : [], 2 : 'ciao'}, {'fl' : [{1 : 'bao'}], 2 : 'ciao'}]
	fitl = fit_into_template(l, ltry)
	tl = fitl.show_list(deeptransparency=True, annotate=True)
	print(tl)

	ll = FixedList(FixedDict({1 : '', 2 : ''}))
	l.append(FixedDict({'fl' : ll, 2 : 4}))
	tl = l.show_list(deeptransparency=True, annotate=True)
	print(tl)
	l2 = FixedList(FixedDict({'flfl' : l, 22 : ''}))
	l2.show_list()

def fixedvaltest():
	fv = FixedValDict(FixedDict({'questo' : '', 'quello' : FixedValDict(FixedDict({'uno' : '', 'due' : {}, 'tre' : []}))}))
	print("TEMPLATE")
	print(fv.show_dictionary(deeptransparency=True, annotate=True))

	fvdata = FixedValDict({'uff' : FixedDict({'questo' : ''}), 'uno' : FixedDict({'questo' : 'heheheh'})})
	print("DATA")
	print(fvdata.show_dictionary(deeptransparency=True, annotate=True))
	pickle.dump(json.dumps(fvdata.show_dictionary(quiet=True, deeptransparency=True, annotate=True)), open("provauno.pkl", 'wb'))
	print("WRITTEN ON FILE:", pickle.load(open("provauno.pkl", 'rb')))
	tmpfv = json.loads(pickle.load(open("provauno.pkl", 'rb')), object_hook=object_decoder_transp)
	print("AFTER LOAD")
	print(tmpfv)
	print("FITTING PROCESS")
	newfv = fit_into_template(fv, tmpfv)
	print(newfv)
	print(newfv.show_dictionary(transparency=True, annotate=True))
	print(newfv.show_dictionary(deeptransparency=True, annotate=True))
	print("DONE")
	
	fvdata = {'A' : {'quello' : {1: {"uno" : "bo"}}}, 'B' : {'questo' : '', 'quello' : {'1' : {'uno' : 'u', 'due' : {1:2}, 'tre' : [1,2,3]}}}}
	fvdatayep = fit_into_template(fv, fvdata)
	print(fvdatayep.show_dictionary(quiet=True, deeptransparency=True, annotate=True))

def fixedvaltest2():
	fv = FixedDict(json.load(open("test_data_entry.json"), object_hook=object_decoder))
	print(fv.show_dictionary(transparency=True, annotate=True))
	print("########################################")
	print(fv.show_dictionary(deeptransparency=True, annotate=True))


def fixedtest():
	f = FixedDict({'1': {}, '2':FixedDict({'a' : '', 'b' : []})})
	fdata = {'1' : {'1a' : 'aa', '1b' : 'bb'}, '2':{'a':'', 'b':[1,2,3,4]}}
	print(f.show_dictionary())
	print(f.show_dictionary(transparency=True, quiet=False, annotate=True))


def read_transf_matrices_json(json_matrices):
    # WARNING: any change on output format MUST be reflected in the mmCIF header reader!
    #  They must agree on output format

    res = []
    keys = [] # List with all key names. It could NOT have been a dict key -> list idx 
    for ijm, json_mx in enumerate(json_matrices):
        tr = {}
        w = []
        for i in range(3):
            v = []
            for j in range(3):
                v.append(json_mx['matrix_{0}_{1}'.format(i+1, j+1)])
            w.append(v)
        tr['matrix'] = w
        w = []
        for i in range(3):
            w.append(json_mx['vector_{0}'.format(i+1)])
        tr['vector'] = w
        res.append(tr)
        keys.append(json_mx['id'])
    return res, keys

def test_roe():
    text = "(1-2)"
    mtx_ids = ['1', '49', '2', '3', '1.5']
    print(read_oper_expression(text, mtx_ids))

    text = "1,2,3"
    mtx_ids = ['1', '3', '2', '5']
    print(read_oper_expression(text, mtx_ids))


def read_oper_expression(text, mtx_ids):
    # text contains the expression
    # mtx_ids is a dict list of matrix ids (strings) in the order of the data structure 
    #  (either mmCIF - containing one absolute list - or json - containing redundant lists,
    #  one per operator, but with file-unique IDs)

    # Remove quotes
    text = text.replace("'", "").replace('"', "")

    if "(" in text:
        # Split by parentheses 
        prodind = text.replace("(", "").split(")")
    else:
        # Jus consider all
        prodind = [text]

    factind = []
    for fact in prodind:
        if not fact:
            continue
        preindexes = [x for x in fact.split(",")]
        indexes = []
        for ind in preindexes:
            if "-" in ind:
                a, b = ind.split("-")
                # collect all IDs from and to two certain IDs (IMPLIES good order in mtx_ids)
                indexes += [[mtx_ids[x]] for x in range(mtx_ids.index(a), mtx_ids.index(b)+1)]
            else:
                indexes += [[ind]]
        factind.append(indexes)
    if len(factind) == 1:
        return factind[0] # List of lists: transformations of the innermost list must be combined. Ex: [[1,2], [3]] == apply combination 1*2, apply 3
    else:
        for ifact in range(len(factind)-1):
            newline = []
            for i in factind[ifact]:
                for j in factind[ifact+1]:
                    newline.append(i + j)
            factind[ifact+1] = newline
        return newline   # List of lists: transformations of the innermost list must be combined. Ex: [[1,2], [3]] == apply combination 1*2, apply 3
            

def parse_mmCIF(mmcif_filename, thr_log_status='ERROR'):
    # WARNING: any change on output format of transformations MUST be reflected in read_transf_matrices_json!
    #  They must agree on output format

    if not os.path.exists(mmcif_filename):
        return {}, {}, {}

    all_loops, loop_keys = parse_mmCIF_format(mmcif_filename)
    header, hreport = parse_mmCIF_header(all_loops, loop_keys)
    pdb_dict, creport = parse_mmCIF_coords(all_loops, loop_keys)

    return header, pdb_dict, {**hreport, **creport}


def parse_mmCIF_header_wrap(mmcif_filename, thr_log_status='ERROR'):

    if not os.path.exists(mmcif_filename):
        return {}, {}
    all_loops, loop_keys = parse_mmCIF_format(mmcif_filename)
    header, hreport = parse_mmCIF_header(all_loops, loop_keys)

    return header, hreport


def parse_mmCIF_coords_wrap(mmcif_filename, thr_log_status='ERROR', scan=False, options={}):

    if not os.path.exists(mmcif_filename):
        return {}, {}
    all_loops, loop_keys = parse_mmCIF_format(mmcif_filename)
    pdb_dict, creport = parse_mmCIF_coords(all_loops, loop_keys, scan=scan, options=options)

    return pdb_dict, creport


def mmCIF_header_stats(stats_location, headers):

    viral_proteins = 0
    depo_dates = []
    resolutions = []
    experiments = []
    biological_assemblies = []
    for header in headers:
        depo_dates.append([int(x) for x in header['deposition_date'].split('-')])
        if type(header['resolution']) is float:
            resolutions.append(header['resolution'])
        else:
            resolutions.append(-1)
        experiments.append(header['experiment'])
        if 'virus' in header['title'].lower() or 'viral' in header['title'].lower():
            viral_proteins += 1
        
        biological_assemblies.append(header['biological_assemblies'])
    
    # Deposition dates
    if len(depo_dates) > 2:
        depo_years = [x[0] for x in depo_dates]
        min_year, max_year = min(depo_years), max(depo_years)
        nbins = 20
        binsize = (max_year - min_year) / nbins
    
        plt.clf()
        fig, ax = plt.subplots(figsize=(30, 8))
        bins = np.arange(min_year, max_year + 1.5*binsize, binsize)
        ax.hist(depo_years, bins=bins)
        ax.set_xlabel('Deposition date (year)')
        ax.set_ylabel('#coordinate files from PDB')
        figname = stats_location + 'PDB_depo_dates.png'
        plt.savefig(figname)
        plt.savefig(figname[:-4]+'.pdf')
    
    # Resolution stacked experiment
    if len(set(resolutions)) > 2:
        exp_types = set(experiments)
        res_stacked = {et : [] for et in exp_types}
        for i in range(len(resolutions)):
            res_stacked[experiments[i]].append(resolutions[i])

        min_res, max_res = int(min(resolutions)), int(max(resolutions))
        binsize = 1

        legend = sorted(list(exp_types))
 
        plt.clf()
        fig, ax = plt.subplots(figsize=(30, 8))
        bins = np.arange(min_res, max_res + 1.5*binsize, binsize)
        ax.hist([res_stacked[et] for et in legend], bins=bins, stacked=True, label=legend, edgecolor='black', linewidth=0.5)
        plt.legend()
        ax.set_xlabel('Resolution (angstrom)')
        ax.set_ylabel('#coordinate files from PDB')
        figname = stats_location + 'resolution.png'
        plt.savefig(figname)
        plt.savefig(figname[:-4]+'.pdf')

    # Biological assemblies
    n_ba = [len(x) for x in biological_assemblies]
    if len(set(n_ba)) > 2:
        n_transf_in_ba1 = [len(x[0]) for x in biological_assemblies]
        min_transf, max_transf = min(n_transf_in_ba1), max(n_transf_in_ba1)

        plt.clf()
        fig, ax = plt.subplots(figsize=(30, 8))
        bins = np.arange(min_transf, max_transf + 1.5, 1)
        ax.hist(n_transf_in_ba1, bins=bins)
        ax.set_xlabel('#matrix transformations')
        ax.set_ylabel('#coordinate files from PDB')
        figname = stats_location + 'transformations_in_ba1.png'
        plt.savefig(figname)
        plt.savefig(figname[:-4]+'.pdf')

    return


def parse_mmCIF_header(all_loops, loop_keys):
    report = {}

    header = FixedDict({
        'name' : '',
        'deposition_date' : '',
        'identifier' : '',
        'title' : '',
        'translrots' : [],
        'ktranslrots' : [],
        'biological_assemblies' : [],
        'experiment' : '',
        'classification' : '',
        'reference' : '',
        'sheet_range' : '',
        'polymers' : '',
        'resolution' : '',
        'space_group' : '',
        'helix_range' : '',
        'chemicals' : '',
        'version' : '',
        'authors' : '',
        'dbref' : [],
        'ss_class' : '' #all-alpha, all-beta, mixed
    })

    # Name
    # Experimental method
    if '_exptl.method' in loop_keys:
        il, ic = loop_keys['_exptl.method']
        header['experiment'] = all_loops[il][1][ic]
    # Resolution
    if '_reflns.d_resolution_high' in loop_keys:
        il, ic = loop_keys['_reflns.d_resolution_high']
        if all_loops[il][1][ic]:
            header['resolution'] = float(all_loops[il][1][ic])
    elif '_em_3d_reconstruction.resolution' in loop_keys:
        il, ic = loop_keys['_em_3d_reconstruction.resolution']
        if all_loops[il][1][ic]:
            header['resolution'] = float(all_loops[il][1][ic])
    # Deposition date
    if '_pdbx_database_status.recvd_initial_deposition_date' in loop_keys:
        il, ic = loop_keys['_pdbx_database_status.recvd_initial_deposition_date']
        header['deposition_date'] = all_loops[il][1][ic]
    # Title
    if '_citation.title' in loop_keys:
        il, ic = loop_keys['_citation.title']
        header['title'] = all_loops[il][1][ic]

    # Secondary structure type
    presence_alpha = False
    presence_beta = False
    if '_struct_conf.conf_type_id' in loop_keys:
        il, ic = loop_keys['_struct_conf.conf_type_id']
        for line in all_loops[il]:
            if 'HELX' in line[ic]:
                presence_alpha = True
            elif 'STRN' in line[ic]:
                presence_beta = True
    if '_struct_sheet_range.sheet_id' in loop_keys:
        presence_beta = True
    ss_class = 'mixed'
    if presence_alpha and (not presence_beta):
        ss_class = 'all-alpha'
    elif presence_beta and (not presence_alpha):
        ss_class = 'all-beta'
    header['ss_class'] = ss_class


    # Transformations
    # There might be multiple "oper_expression", each giving a different complex
    # We always take the first one (that should be the most complete one)
    if '_pdbx_struct_assembly_gen.oper_expression' in loop_keys:
        transf = {'chains', 'operators'}
        oeil, oeic = loop_keys['_pdbx_struct_assembly_gen.oper_expression']
        _, ailic = loop_keys['_pdbx_struct_assembly_gen.asym_id_list']
        _, aidic = loop_keys['_pdbx_struct_assembly_gen.assembly_id']

        mil, mic = loop_keys['_pdbx_struct_oper_list.matrix[1][1]']
        _, midic = loop_keys['_pdbx_struct_oper_list.id']

        # Translrot
        # IMPLIES that the key order in loop in  mmCIF file is: 
        # matrix[1][1], matrix[1][2], matrix[1][3], vector[1], matrix[2][1], ..., vector[3]
        matrices = []
        matrix = []
        row = []
        vector = []
        kmatrices = []
        for im in range(1, len(all_loops[mil])):
            for i in range(mic, mic+12):
                if i%4<3:
                    row.append(all_loops[mil][im][i])
                else:
                    vector.append(all_loops[mil][im][i])
                    matrix.append(row)
                    row = []
            matrices.append({'matrix' : matrix, 'vector' : vector})
            kmatrices.append(all_loops[mil][im][midic])
            matrix = []
            vector = []
        header['translrots'] = matrices
        header['ktranslrots'] = kmatrices
        
        # Operators
        assembly = []
        assid = ''
        for iop in range(1, len(all_loops[oeil])):
            transf = {}
            transf['operators'] = read_oper_expression(all_loops[oeil][iop][oeic], kmatrices)
            transf['chains'] = all_loops[oeil][iop][ailic].split(',')
            if all_loops[oeil][iop][aidic] != assid:
                if assembly:
                    header['biological_assemblies'].append(assembly)
                    assembly = []
                assid = all_loops[oeil][iop][aidic]
            assembly.append(transf)
        if assembly:
            header['biological_assemblies'].append(assembly)

    # Database reference
    if '_struct_ref.db_name' in loop_keys:
        dbil, dbic = loop_keys['_struct_ref.db_name']
        _, eic = loop_keys['_struct_ref.entity_id']
        _, accic = loop_keys['_struct_ref.pdbx_db_accession']

        for dbi in range(1, len(all_loops[dbil])):
            dbr = {'entity' : '', 'accession' : '', 'database' : ''}
            dbr['entity'] = all_loops[dbil][dbi][eic]
            dbr['accession'] = all_loops[dbil][dbi][accic]
            dbr['database'] = all_loops[dbil][dbi][dbic]
            header['dbref'].append(dbr)

    return header.show_dictionary(quiet=True), report


def parse_mmCIF_coords(all_loops, loop_keys, scan=False, model_to_consider=1, options={}):
    record = {
        "unk_res" : set(),
        "unrecognized_res" : set(),
        "mse_res" : set(),
        "het_res" : set(),
        "double_record" : set(),
        "double_record_auth" : set(),
        "altlocs" : {},
        "pdb_format" : None,
        "lines" : None,
        "chains" : [],
        "new_chains" : [],
        "merged_chains" : []   # NOT IMPLEMENTED
    }

    polymer_entities = set()
    if '_entity.id' in loop_keys:
        eidil, eidic = loop_keys['_entity.id']
        _, etic = loop_keys['_entity.type'] 
        no_entity_info = False
        for il in range(1,len(all_loops[eidil])):
            if all_loops[eidil][il][etic] == 'polymer':
                polymer_entities.add(all_loops[eidil][il][eidic])
    else:
        no_entity_info = True

    atil, groupic = loop_keys['_atom_site.group_PDB']  # ATOM/HETATM
    _, atidic = loop_keys['_atom_site.id'] # ATOM ID
    _, symic = loop_keys['_atom_site.id'] # ELEMENT
    _, nameic = loop_keys['_atom_site.label_atom_id'] # ATOM NAME
    _, altic = loop_keys['_atom_site.label_alt_id'] # ALTLOC
    _, resnic = loop_keys['_atom_site.label_comp_id'] # 3-LETTER RESIDUE NAME
    _, chic = loop_keys['_atom_site.label_asym_id'] # CHAIN (mmCIF CONVENTION)
    _, entic = loop_keys['_atom_site.label_entity_id'] # ENTITY ID (FOR CONVERSION TABLE)
    _, residic = loop_keys['_atom_site.label_seq_id'] # RESIDUE ID
    _, inscic = loop_keys['_atom_site.pdbx_PDB_ins_code'] # INSERTION CODE (ALPHABETICAL COMPLEMENT TO RESID)
    _, xic = loop_keys['_atom_site.Cartn_x'] # X
    _, yic = loop_keys['_atom_site.Cartn_y'] # Y
    _, zic = loop_keys['_atom_site.Cartn_z'] # Z
    _, occic = loop_keys['_atom_site.occupancy'] # OCCUPANCY
    _, authatic = loop_keys['_atom_site.auth_atom_id'] # AUTHOR ATOM NAME
    _, authidic = loop_keys['_atom_site.auth_seq_id'] # AUTHOR RESID (FOR HOLE CHECKS)
    _, authresnic = loop_keys['_atom_site.auth_comp_id'] # AUTHOR 3-LETTER RESIDUE NAME
    _, authchic = loop_keys['_atom_site.auth_asym_id'] # AUTHOR CHAIN (FOR CONVERSION TABLE)
    _, modif = loop_keys['_atom_site.pdbx_PDB_model_num'] # MODEL NUMBER (FOR KEEPING TH FIRST)
    _, bfic = loop_keys['_atom_site.B_iso_or_equiv'] # B FACTOR

    bb_atoms = ['N', 'CA', 'C', 'O', 'OXT', 'CB'] # This is only used for UNKs and UNKs also keep CBs
    groups = ["ATOM", "HETATM"]

    golines = []
    identif = set()
    newidentif = set()
    chains = []
    new_chains = []
    altlocs = {}
    occupancies = {}
    normocc = {}
    for il in range(1,len(all_loops[atil])):
        line = all_loops[atil][il]

        if int(line[modif]) != model_to_consider:
            continue

        entity = line[entic] 
        if (not no_entity_info) and entity not in polymer_entities:
            continue

        group = line[groupic]

        # Residue identification
        resname = line[authresnic] # Author's convention
        resid = int(line[authidic]) # Author's convention
        inscode = line[inscic]
        altloc = line[altic]
        res = ((altloc, resid, inscode), resname)

        # Keep only ATOM and HETATM
        if group not in groups:
            #print("NOT ATOM NOR HETATM")
            continue

        atname = line[authatic] # Author's convention used
        ch = line[authchic] # Author's convention
        newch = line[chic]

        # Compile list of chains and list of new chains
        #  (they are kept in the same order)
        if ch not in chains:
            chains.append(ch)
            new_chains.append(newch)
            altlocs[ch] = {}

        # QUI METTI CODICE PER CONSECUTIVE AUTHOR CHAINS E CONTROLLA QUANDO SUCCEDE

        # Eliminate double lines
        if (atname, res, ch) in identif:
            record["double_record"].add((atname, res, ch))
            continue
        else:
            identif.add((atname, res, ch))

        newatname = line[nameic]
        newresid = int(line[residic])
        newresname = line[resnic]
        newres = ((altloc, newresid, inscode), newresname)
        newch = line[chic]
        # Eliminate double lines wrt new notation
        if (newatname, newres, newch) in newidentif:
            record["double_record_new"].add((newatname, newres, newch))
            continue
        else:
            newidentif.add((newatname, newres, newch))
       
        # Classify AltLocs
        res = ((resid, inscode), resname)
        occupancy = float(line[occic])
        if altloc:
            if (res, ch) not in occupancies:
                occupancies[(res, ch)] = 0
                normocc[(res, ch)] = 0
            occupancies[(res, ch)] += occupancy
            normocc[(res, ch)] += 1
            if altloc not in altlocs[ch]:
                altlocs[ch][altloc] = set()
            altlocs[ch][altloc].add((res, ch))
        
        golines.append(il)


    al = {}
    for ch in altlocs:
        al[ch] = {x for x in altlocs[ch]}
    record["altlocs"] = al

    # Decide what AltLoc to keep: highest occupancy, then most abundant
    best_altlocs = {}
    for ch in altlocs:
        best_altloc, best_occ, best_norm = '', 0, 0
        for altloc in altlocs[ch]:
            tot_occupancy, tot_norm = 0, 0
            for resch in altlocs[ch][altloc]:
                if normocc[resch] > 0:
                    tot_occupancy += occupancies[resch] / normocc[resch]
                tot_norm += 1
            if tot_norm > 0:
                tot_occupancy /= tot_norm
            if (tot_occupancy > best_occ) or (tot_occupancy == best_occ and tot_norm > best_norm):
                best_altloc = altloc
                best_occ = tot_occupancy
                best_norm = tot_norm
        best_altlocs[ch] = best_altloc

    #ils = [x for x in range(1,len(all_loops[atil])) if x not in dellines]
    golines2 = []
    for il in golines:
        line = all_loops[atil][il]
        altloc = line[altic]
        ch = line[authchic] # Author's convention
        if altloc != '' and (ch not in best_altlocs or altloc != best_altlocs[ch]):
            continue
        golines2.append(il)

    # Can write PDB?
    # line number <100000 & mono-character chain
    if len(golines2) < 100000 and max([len(x) for x in chains]) == 1:
        gopdb = True
    else:
        gopdb = False

    record["pdb_format"] = gopdb
    record["lines"] = len(golines2)
    record["chains"] = chains
    record["new_chains"] = new_chains

    old_ch = ''
    coords = {}
    resids = {}
    pdbseq = {}
    b_factor = {}
    for il in golines2:
        line = all_loops[atil][il]
        #print("LINE", line)

        group = line[groupic]
        # Chain information
        ch = line[authchic] # Author's convention
        # Residue identification
        resname = line[authresnic] # Author's convention
        resid = int(line[authidic]) # Author's convention
        inscode = line[inscic]
        res = ((resid, inscode), ch)
        # Atom information
        atname = line[authatic] # Author's convention
        bfact = float(line[bfic]) if string_is_float(line[bfic]) else -1

        x, y, z = float(line[xic]), float(line[yic]), float(line[zic])

        if group == "ATOM":
            # UNK are kept
            if resname == "UNK":
                record["unk_res"].add((res, ch))
            # ??? -> UNK (because we are in polymers!)
            elif from3to1(resname) == "X":
                record["unrecognized_res"].add((res, ch))
                resname = "UNK"
            # Sidechain atoms of UNK are removed (except CB)
            if atname not in bb_atoms:
                continue
        elif group == "HETATM":
            resname = line[authresnic] # Author's convention
            # MSE -> MET
            if resname == "MSE":
                record["mse_res"].add((res, ch))
                if atname == "SE":
                    atname = "S"
                resname = "MET"
            # Other HETATMs are deleted (can be glycans, or stuff like that)
            else:
                record["het_res"].add((res, ch))
                continue

        # Compile data structure
        res = ((resid, inscode), resname)
        if ch not in coords:
            coords[ch] = {}
            resids[ch] = []
            pdbseq[ch] = []
            b_factor[ch] = []
            if old_ch and b_factor[ch] and b_norm > 0:
                b_factor[ch][-1] /= b_norm
            old_ch = ch
        if res not in coords[ch]:
            coords[ch][res] = {}
            if b_factor[ch] and b_norm > 0:
                b_factor[ch][-1] /= b_norm
            b_norm = 0
            b_factor[ch].append(0)
        resids[ch].append((resid, inscode))
        pdbseq[ch].append(from3to1(resname))
        coords[ch][res][atname] = (x, y, z)
        if bfact != -1:
            b_factor[ch][-1] += bfact
            b_norm += 1

    if old_ch and ch in b_factor and b_factor[ch] and b_norm > 0:
        b_factor[ch][-1] /= b_norm

    # Until OPM does not convert to mmCIF, there will never be DUM residues in mmCIF structures
    pdb_dict = {
        'COORDS' : coords,  # (Per chain) nested dictionary <nested tuple>res -> <str>atname -> <3-tuple>coords
        'DUMMIES' : {},  # No dummies in mmCIF until OPM converts to mmCIF
        'PDBSEQ' : { ch : pdbseq[ch] for ch in pdbseq }, # (Per chain) List of 1-character residues
        'AVG_BFACTOR' : b_factor,    # (Per chain) List of B factor averages over residues
        'RESIDS' : resids            # (Per chain) List of pairs (<int>resid, <str>inscode)
    }

    if scan:
        if options:
            report = scan_struct(options, coords, {})
            for x in report:
                if x not in record:
                    record[x] = report[x]
        else:
            print("scan option can only be used with options data structure!")
            exit(1)

    return pdb_dict, record
         

def parse_mmCIF_format(mmcif_filename):    
    """Main "mindless" parser of mmCIF format
    Returns a list of tables with all the processed keys/values
    Each table has a 1-line header with keys
    """

    if not os.path.exists(mmcif_filename):
        print(f"ERROR: file {mmcif_filename} does not exist")
        return [], {}

    all_loops = []
    with open(mmcif_filename) as f:
        record_type = ""
        this_loop = []
        for line in f:
            if not line.strip() or line.startswith("data_"):
                continue
            if line.startswith("#"):
                if record_type == "def":
                    this_loop = list(zip(*this_loop))
                if this_loop:
                    all_loops.append(this_loop)
                this_loop = []
                record_type = ""
                continue
            elif record_type != "longline" and line.startswith("_") and record_type != "loop":
                record_type = "def"
                n_vals = 2
                vals = []
                longfield = ""
            elif line.startswith("loop_"):
                record_type = "loop"
                keys = []
                vals = []
                n_vals = 0
                longfield = ""
                continue

            # Define keys
            if line.startswith("_"):
                if record_type == "def":
                    keys = [line.split()[0]]
                elif record_type == "loop":
                    keys.append(line.strip())
                    n_vals += 1
            # Record vals
            if record_type in {"def", "loop"}:
                if line.startswith("_") and record_type == "loop":
                    continue
                # If we are in a long field (starts with ;, then several lines follow, ends with ; and newline)
                if line.startswith(";"):
                    if line.strip() == ";":
                        vals.append(longfield)
                        longfield = ""
                    else:
                        longfield = line[1:].strip()
                elif longfield:
                    longfield += line.strip()
                # Otherwise count the fields on this line
                else:
                    vals_on_line = []
                    wsfield = False
                    for w in line.split():
                        if w[0] == "'":
                            if w[-1] != "'":
                                wsfield = True
                                field = w[1:]
                                continue
                            else:
                                field = w[1:-1] 
                        elif w[-1] == "'" and wsfield:
                            field += " " + w[:-1]
                            wsfield = False
                        elif wsfield:
                             field += " " + w
                        elif w in ['.', '?']:  # mmCIF placeholders for empty fields
                             field = ''
                        else:
                             field = w
                        if not wsfield:
                            vals_on_line.append(field)
                    vals += vals_on_line
                if len(vals) == n_vals:
                    if record_type == "def":
                        this_loop.append(vals)
                    elif record_type == "loop":
                        if not this_loop:
                            this_loop.append(keys)
                        this_loop.append(vals)
                    vals = []

    key_d = {}
    for iloop, loop in enumerate(all_loops):
        for ik, k in enumerate(loop[0]):
            key_d[k] = (iloop, ik)

    return all_loops, key_d


# DEPRECATED - WE NO LONGER PARSE PDB HEADERS (ONLY MMCIF HEADERS)
def parse_PDB_header(pdb_filename, thr_log_status='ERROR'):

    def years(yy):
        if int(yy) < 50:
            return int('20'+yy)
        else:
            return int('19'+yy)

    this_name = parse_PDB_header.__name__

    header = FixedDict({
        'name' : '',
        'deposition_date' : '',
        'identifier' : '',
        'title' : '',
        'translrots' : [],
        'ktranslrots' : [],
        'biological_assemblies' : [],
        'experiment' : '',
        'classification' : '',
        'reference' : '',
        'sheet_range' : '',
        'polymers' : '',
        'resolution' : '',
        'space_group' : '',
        'helix_range' : '',
        'chemicals' : '',
        'version' : '',
        'authors' : '',
        'dbref' : [],
        'ss_class' : '' #all-alpha, all-beta, mixed
    })

    """
    header = FixedDict({
        'name' : '', 
        'deposition_date' : '', 
        'identifier' : '', 
        'title' : '', 
        'biomoltrans' : '', 
        'experiment' : '', 
        'classification' : '', 
        'reference' : '', 
        'sheet_range' : '', 
        'polymers' : '', 
        'resolution' : '', 
        'space_group' : '', 
        'helix_range' : '', 
        'chemicals' : '', 
        'version' : '', 
        'authors' : '', 
        'dbref' : '',
        'ss_class' : ''
    })
    """

    # Check existence of file
    if not os.path.exists(pdb_filename):
        if pdb_filename[-4:] == ".pdb" and os.path.exists(pdb_filename+'.gz'):
            os.system("cd {0}; gunzip {1}".format(os.path.dirname(pdb_filename), os.path.basename(pdb_filename+'.gz')))
        else:
            nof_err = ('CRITICAL', this_name, "File {0} not found".format(pdb_filename))
            print_log(nof_err)

    months = ['JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN', 'JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC']
    biomoltrans = {}
    biomoltrans_i = []
    title = ""
    experiment = ""
    dbref = {}
    saveline = ''
    with open(pdb_filename) as pdb_file:
        for nl, line in enumerate(pdb_file):
            if line.startswith("ATOM"):
                break
            # HEADER
            if line[:6] == 'HEADER':
                header['name'] = line[10:50].strip()
                dd, mm, yy = line[50:59].split("-")
                header['deposition_date'] = [int(dd), months.index(mm), years(yy)]
                header['identifier'] = line[62:66].lower()
            # TITLE
            elif line[:5] == 'TITLE':
                title += line[10:].strip() + " "
            # EXPDTA
            elif line[:6] == 'EXPDTA':
                experiment = line[10:].strip()
            # DBREF
            elif line[:5] == 'DBREF':
                db_name = line[26:29].strip()
                # Only UniProt references are parsed
                if db_name != "UNP":
                    continue
                # DBREF can have two line or just one
                # If two lines:
                if line[5] == '1':
                    saveline = line
                    continue
                elif line[5] == '2':
                    ch_pdb = saveline[12]
                    db_access = saveline[42:].strip()
                    db_id = line[14:].strip()
                    pdb_s = int(saveline[14:18])
                    pdb_e = int(saveline[20:24])
                    db_s = int(line[33:].split()[0])
                    db_e = int(line[33:].split()[1])
                # If just one line:
                else:
                    ch_pdb = line[12]
                    db_access = line[42:54].strip()
                    db_id = line[33:41].strip()
                    pdb_s = int(line[14:18])
                    pdb_e = int(line[20:24])
                    db_s = int(line[55:60])
                    db_e = int(line[62:67])
                dbref[ch_pdb] = ((db_id, db_access), (pdb_s, pdb_e), (db_s, db_e))
            # REMARK
            elif line[:6] == 'REMARK' and len(line.split()) > 2:
                # REMARK 2 (resolution)
                if int(line[6:10]) == 2:
                    if 'RESOLUTION' in line and experiment and ("RAY" in experiment or "ELECTRON" in experiment):
                        try:
                            header['resolution'] = float(line.split()[3])
                        except:
                            pass
                # REMARK 350 (biomatrices)
                # TO DO Thes matrices do not follow the new format!
                elif int(line[6:10]) == 350:
                    if 'BIOMOLECULE' in line.split()[2]:
                        if biomoltrans_i:
                            biomoltrans[biomolecule_id] = biomoltrans_i  # ERROR! Implies only one transformation per biomolecule!
                        biomolecule_id = int(line.split()[3])
                        biomoltrans_i = [] # Matrix (rotation and translation)
                    elif 'CHAINS:' in line:
                        if "APPLY THE FOLLOWING TO CHAINS:" in line:
                            biomolecule_chains = [x.strip() for x in line[line.index(':')+1:].split(',') if x.strip()]
                        elif 'AND CHAINS' in line:
                            biomolecule_chains += [x.strip() for x in line[line.index(':')+1:].split(',') if x.strip()]				
                    elif 'BIOMT' in line.split()[2]:
                        if 'BIOMT1' == line.split()[2]:
                            biomoltrans_i.append(biomolecule_chains)
                        biomoltrans_i.append(line[23:].strip())
        if biomoltrans_i:
            biomoltrans[biomolecule_id] = biomoltrans_i
        header['title'] = title
        header['biomoltrans'] = biomoltrans  # INCORRECT DATA
        header['experiment'] = experiment
        header['dbref'] = dbref
        hreport = {}

    return header.show_dictionary(quiet=True), hreport


def write_ENC(str_data_pdbi, c_dict_pdbi, out_filename, ignoredpc_filename, pad_fn, cache_filename='', thr_log_status='ERROR', no_pass=False, only_pass=False, chain_selection=None):
	this_name = write_ENC.__name__

	order = pdb_atom_dictionary(pad_fn)

	ignorepc = []
	with open(ignoredpc_filename) as ig_file:
		for line in ig_file:
			if line.strip():
				ignorepc.append(line.strip())
	old_title = None
	with open(out_filename, 'w') as enc_file:
		enc_file.write("HEADER    {0:40}{1}   {2}\n".format(str_data_pdbi['FROM_PDB']['name'], str_data_pdbi['FROM_PDB']['deposition_date'], str_data_pdbi['FROM_PDB']['pdb_id'].upper()))
		if str_data_pdbi['ENCOMPASS']['name']:
			title = str_data_pdbi['ENCOMPASS']['name']
		else:
			title = str_data_pdbi['FROM_PDB']['title']
		while title.strip():
			if len(title)>70 and title[69]!=" " and title[70]!=" ":
				words = " ".join([x for x in title[:70].split()[:-1]])
			else:
				words=title[:70]
			title = title[len(words):]
			if old_title is not None and title==old_title: # Prevents infinite loops
				break
			old_title = title
			enc_file.write("TITLE     {0:70}\n".format(words))
		enc_file.write("EXPDTA    {0}\n".format(str_data_pdbi['FROM_PDB']['experimental_method'].upper()))
		linelength = 60
		enc_file.write("REMARK    0 EncoMPASS version 2.0\n")
		ie = 1
		for entry in str_data_pdbi['PASSPORT']:
			if no_pass:
				continue
			ignore = False
			for ipc in ignorepc:
				if '_'+ipc in entry['location']:
					ignore = True
			if ignore:
				continue
			messages = entry['message']
			first_line = True
			ind = ""
			for message in messages.split('\n'):
				words = message.split()
				if not words[0]:
					continue
				mpart = ""
				a_capo = False
				endline = False
				proceed = True
				w = words[0]
				while words or w or mpart:
					if len(mpart) + len(w) + 1 > linelength or (mpart and (a_capo or endline or (w and w[0] in ['[', '{']))):
						if not ('{' in mpart or '[' in mpart) and ('}' in mpart or ']' in mpart):
							ind = ind[1:]
						if first_line:
							enc_file.write("REMARK    0 MSG{0:4d}: {1}\n".format(ie, ind+mpart))
							ie += 1
						else:
							enc_file.write("REMARK    0          {0}\n".format(ind+mpart))
						first_line = False
						a_capo = False
						if ('{' in mpart or '[' in mpart) and not ('}' in mpart or ']' in mpart):
							ind += " "
						mpart = ""
					ipar1a, ipar1b = -1, -1
					if ']' in w and '}' in w:
						ipar1a = w.index(']')
						ipar1b = w.index('}')
					if w and w[0] in ['[', '{'] and w.strip() != mpart.strip():
						mpart = w[0]
						w = w[1:]
						a_capo= False
						proceed = False
					elif ']' in w and (ipar1a+ipar1b == -2 or ipar1a < ipar1b):
						ipar = w.index(']')
						proceed = False
						ipar2a, ipar2b = -1, -1
						if ']' in w[ipar+1:] and '}' in w[ipar+1:]:
							ipar2a = w[ipar+1:].index(']')
							ipar2b = w[ipar+1:].index('}')
						if ']' in w[ipar+1:] and (ipar2a+ipar2b == -2 or ipar2a < ipar2b):
							ipar2 = w[ipar+1:].index(']')
							mpart += w[:ipar+1+ipar2]
							w = w[ipar+1+ipar2:]
						elif '}' in w[ipar+1:] and (ipar2a+ipar2b == -2 or ipar2a > ipar2b):
							ipar2 = w[ipar+1:].index('}')
							mpart += w[:ipar+1+ipar2]
							w = w[ipar+1+ipar2:]
						else:
							proceed = True
						a_capo = True
					elif '}' in w and (ipar1a+ipar1b == -2 or ipar1a > ipar1b):
						ipar = w.index('}')
						proceed = False
						ipar2a, ipar2b = -1, -1
						if ']' in w[ipar+1:] and '}' in w[ipar+1:]:
							ipar2a = w[ipar+1:].index(']')
							ipar2b = w[ipar+1:].index('}')
						if ']' in w[ipar+1:] and (ipar2a+ipar2b == -2 or ipar2a < ipar2b):
							ipar2 = w[ipar+1:].index(']')
							mpart += w[:ipar+1+ipar2]
							w = w[ipar+1+ipar2:]
						elif '}' in w[ipar+1:] and (ipar2a+ipar2b == -2 or ipar2a > ipar2b):
							ipar2 = w[ipar+1:].index('}')
							mpart += w[:ipar+1+ipar2]
							w = w[ipar+1+ipar2:]
						else:
							proceed = True
						a_capo = True
					if proceed:
						mpart += w + " "
						words = words[1:]
						if words:
							w = words[0]
						elif endline:
							break
						else:
							endline = True
					else:
						proceed = True
	
		if str_data_pdbi['FROM_PDB']['resolution'] and type(str_data_pdbi['FROM_PDB']['resolution']) == str and string_is_float(str_data_pdbi['FROM_PDB']['resolution']):
			enc_file.write("REMARK    2 RESOLUTION.  {0:6.2f} ANGSTROMS.\n".format(float(str_data_pdbi['FROM_PDB']['resolution'])))
		if str_data_pdbi['ENCOMPASS']['structure']['UniProt_stoichiometry']:
			enc_file.write("REMARK  999 STOICHIOMETRY QUANTITY                 CORRESP_CHAINS                      UNP_CODES\n")
			noup_d = {}
			for k in str_data_pdbi['ENCOMPASS']['structure']['UniProt_stoichiometry']:
				if str_data_pdbi['ENCOMPASS']['structure']['UniProt_stoichiometry'][k]:
					chs = str_data_pdbi['ENCOMPASS']['structure']['UniProt_stoichiometry'][k]
					enc_file.write("REMARK  999 STOICHIOMETRY {0:8d} {1:>30} {2:>30}\n".format(len(chs), ",".join(chs), ",".join(k.split(","))))
				else:
					for c in chs:
						ic = str_data_pdbi['ENCOMPASS']['structure']['kchains'].index(c) 
						s = str_data_pdbi['ENCOMPASS']['structure'][ic]['sequence']
						if s not in noup_d:
							noup_d[s] = set()
						noup_d[s].add(c)
			for seq in noup_d:
				chs = noup_d[seq]
				enc_file.write("REMARK  999 STOICHIOMETRY {0:8d} {1:>30} {2:>30}\n".format(len(chs), ",".join(chs), "N/A"))

					
		for ich, ch in enumerate(str_data_pdbi['ENCOMPASS']['structure']['kchains']):
			unp_acc_list = str_data_pdbi['ENCOMPASS']['structure']['chains'][ich]['UniProt_acc']
			for acc, (bn, en), (strb, stre) in unp_acc_list:
				enc_file.write("DBREF  {0} {1} {2:4d}{3} {4:4d}{5} UNP  {6:>8} {7:>11}    {8:5d}  {9:5d}\n".format(str_data_pdbi['FROM_PDB']['pdb_id'].upper(), ch, strb, ' ', stre, ' ', acc, acc, bn+1, en+1))

		if only_pass:
			return
		atid = 0
		for chain in sorted(list(c_dict_pdbi['COORDS'].keys())):
			if chain_selection and chain not in chain_selection:
				continue
			for resid_name in c_dict_pdbi['COORDS'][chain]:
				resid, resname = resid_name
				for atname in order[resname]:
					if atname in c_dict_pdbi['COORDS'][chain][resid_name]:
						atid += 1
						cx, cy, cz = c_dict_pdbi['COORDS'][chain][resid_name][atname]
						if len(atname) < 4:
							atname_real = ' ' + atname
						else:
							atname_real = atname
						enc_file.write('ATOM  {0:5d} {1:<4} {2:3} {3:1}{4:4d}{5:1}   {6: 8.3f}{7: 8.3f}{8: 8.3f}{9:24}\n'.format(atid, atname_real, resname, chain, resid[0], resid[1], float(cx), float(cy), float(cz), ' '*24))
			atid += 1
			if len(atname) < 4:
				atname_real = ' ' + atname
			else:
				atname_real = atname
			enc_file.write("TER   {0:5d} {1:<4} {2:3} {3:1}{4:4d}{5:1}\n".format(atid, "", resname, chain, resid[0], resid[1]))
		for atname, coords in c_dict_pdbi['DUMMIES']:
			cx, cy, cz = coords
			atid += 1
			if len(atname.strip()) < 4:
				atname_real = ' ' + atname
			else:
				atname_real = atname
			enc_file.write('HETATM{0:5d} {1:<4} {2:3} {3:1}{4:4}    {5: 8.3f}{6: 8.3f}{7: 8.3f}{8:24}\n'.format(atid, atname, 'DUM', '0', int(resid[0])+1, float(cx), float(cy), float(cz), ' '*24))
		enc_file.write("END\n")

	if cache_filename:
		write_mmCIF(str_data_pdbi, coord_dict, cache_filename, pad_fn)
	

def clash_spheres():
    """Central atom and radius enclosing at all times
    a freely-moving residue of a certain type
    """

    spheres = { 
        'ALA' : ('CA', 4),
        'ARG' : ('CA', 10), #('CG', 7)
        'ASN' : ('CA', 7), #('CB', 5.5)
        'ASP' : ('CA', 7), #('CB', 5.5)
        'CYS' : ('CA', 5.5), #('CB', 4)
        'GLN' : ('CA', 7), #('CB', 5.5)
        'GLU' : ('CA', 7), #('CB', 5.5)
        'GLY' : ('CA', 2.5),
        'HIS' : ('CA', 7), #('CB', 5.5)
        'ILE' : ('CA', 5.5),
        'LEU' : ('CA', 5.5),
        'LYS' : ('CA', 8.5), #('CG', 5.5)
        'MET' : ('CA', 7), #('CB', 5.5)
        'MSE' : ('CA', 7), #('CB', 5.5)
        'PHE' : ('CA', 8.5), #('CB', 7)
        'PRO' : ('CA', 4),
        'SER' : ('CA', 4),
        'THR' : ('CA', 2),
        'TRP' : ('CA', 10), #('CG', 7)
        'TYR' : ('CA', 10), #('CB', 8.5)
        'VAL' : ('CA', 4),
        'UNK' : ('CA', 4)
    }

    return spheres


def pdb_atom_dictionary(pad_fn):
    this_name = pdb_atom_dictionary.__name__

    if not os.path.exists(pad_fn):
        print_log((
            'CRITICAL',
             this_name,
             f"Could not open atom dict file {pad_fn}"
        ))

    order = {}
    with open(pad_fn) as in_file:
        for line in in_file:
            resname, atname = line.split()
            if resname not in order:
                order[resname] = []
            order[resname].append(atname)
    for resname in order:
        order[resname].append('OXT')
    return order


def write_mmCIF(str_data, coord_dict, output_filename, pad_fn):
	if not str_data or not coord_dict or not coord_dict['COORDS']:
		print(not str_data, not coord_dict, not coord_dict['COORDS'], output_filename)
		print(str_data['status'])
		print("ERROR777")
		exit(1)

	order = pdb_atom_dictionary(pad_fn)

	pdbi = str_data['FROM_PDB']['pdb_id']
	with open(output_filename, 'w') as out_file:
		out_file.write('data_{0}\n'.format(pdbi.upper()))

		out_file.write('#\n_entry.id {0}\n#\n'.format(pdbi.upper()))

		out_file.write('loop_\n')
		out_file.write('_atom_site.group_PDB\n')
		out_file.write('_atom_site.id\n')
		out_file.write('_atom_site.label_atom_id\n')
		out_file.write('_atom_site.label_comp_id\n')
		out_file.write('_atom_site.label_asym_id\n')
		out_file.write('_atom_site.label_seq_id\n')
		out_file.write('_atom_site.pdbx_PDB_ins_code\n')
		out_file.write('_atom_site.Cartn_x\n')
		out_file.write('_atom_site.Cartn_y\n')
		out_file.write('_atom_site.Cartn_z\n')
		atid = 0
		for chain in sorted(list(coord_dict['COORDS'].keys())):
			for resid_name in coord_dict['COORDS'][chain]:
				resid, resname = resid_name
				for atname in order[resname]:
					if atname in coord_dict['COORDS'][chain][resid_name]:
						atid += 1
						cx, cy, cz = coord_dict['COORDS'][chain][resid_name][atname]
						inscode = '?' if not resid[1] else resid[1]
						out_file.write('{0:6} {1:<8d} {2:<4} {3:3} {4:1} {5:<4d}  {6:1} {7: 10.3f} {8: 10.3f} {9: 10.3f}\n'.format('ATOM', atid, atname, resname, chain, resid[0], inscode, float(cx), float(cy), float(cz)))
		if 'DUMMIES' in coord_dict and coord_dict['DUMMIES']:
			for atname, coords in coord_dict['DUMMIES']:
				cx, cy, cz = coords
				atid += 1
				out_file.write('{0:6} {1:<8d} {2:<4} {3:3} {4:1} {5:<4}  ? {6: 10.3f} {7: 10.3f} {8: 10.3f}\n'.format('HETATM', atid, atname, 'DUM', '.', '.', float(cx), float(cy), float(cz)))
		out_file.write('#\n')


def retrieve_coords_from_CIF(cif_filename, four_outs=False):
	loop_header = False
	atom_lines = False
	c_dict = {'COORDS' : {}, 'DUMMIES' : []}
	#print(cif_filename)
	with open(cif_filename) as cif_file:
		for line in cif_file:
			if line.startswith('loop_'):
				loop_header = True
				argcount = []
				continue
			if loop_header:
				if line.startswith('_'):
					if line.startswith('_atom_site.label_atom_id'):
						argcount.append('atname')
					elif line.startswith('_atom_site.label_comp_id'):
						argcount.append('resname')
					elif line.startswith('_atom_site.label_asym_id'):
						argcount.append('chain')
					elif line.startswith('_atom_site.label_seq_id'):
						argcount.append('r')
					elif line.startswith('_atom_site.pdbx_PDB_ins_code'):
						argcount.append('ichar')
					elif line.startswith('_atom_site.Cartn_x'):
						argcount.append('cx')
					elif line.startswith('_atom_site.Cartn_y'):
						argcount.append('cy')
					elif line.startswith('_atom_site.Cartn_z'):
						argcount.append('cz')
					else:
						argcount.append(None)
				elif line.startswith('ATOM'):
					loop_header = False
					atom_lines = True
			if atom_lines and (line.startswith('ATOM') or line.startswith('HETATM')):
				fields = line.split()
				coordinates = [float(fields[argcount.index('cx')]), float(fields[argcount.index('cy')]), float(fields[argcount.index('cz')])]
				resname = fields[argcount.index('resname')]
				atname = fields[argcount.index('atname')]
				if line.startswith('HETATM'):
					if resname == "DUM":
						c_dict['DUMMIES'].append((atname, coordinates))
					continue
				chain = fields[argcount.index('chain')]
				r = int(fields[argcount.index('r')])
				ichar = '' if fields[argcount.index('ichar')].strip() == '?' else fields[argcount.index('ichar')]
				resid = (r, ichar)
				if chain not in c_dict['COORDS']:
					c_dict['COORDS'][chain] = {}
				if (resid, resname) not in c_dict['COORDS'][chain]:
					c_dict['COORDS'][chain][(resid, resname)] = {}
				c_dict['COORDS'][chain][(resid, resname)][atname] = coordinates
#				if resname == "UNK": 
#					print(resname)
	return c_dict


def read_coord_from_PDB_text(text):
    coords = {}	# Coordinates dictionary: (chain)(resid, resname)(atname) -> coordinates
    dummies = []	# Dummy atoms list: (atname, coordinates)
    for nl, line in enumerate(text.split('\n')):
        if not line:
            continue
        if not (line.startswith("ATOM") or line.startswith("HETATM")):
            continue
        if not line[21]:
            ch_name = '_'
        else:
            ch_name = line[21]
        resname = line[17:20]
        atname = line[12:16].strip()
        coordinates = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        if line.startswith("HETATM"):
            if resname == "DUM":
                dummies.append((atname, coordinates))
            continue

        r = int(line[22:26].strip())
        ichar = line[27].strip()
        if ch_name not in coords:
            coords[ch_name] = {}
        if ((r, ichar), resname) not in coords[ch_name]:
            coords[ch_name][((r, ichar), resname)] = {}
        if atname not in coords[ch_name][((r, ichar), resname)]:
            coords[ch_name][((r, ichar), resname)][atname] = coordinates

    c_dict = {'COORDS' : coords, 'DUMMIES' : dummies}
    return c_dict


def scan_struct_from_PDB_text(options, text='', coord_filename='', with_coords=False):

    if text:
        c_dict = read_coord_from_PDB_text(text)
    elif coord_filename:
        coord_file = open(coord_filename)
        text = coord_file.read()
        c_dict = read_coord_from_PDB_text(text)

    report = scan_struct(options, coords=c_dict['COORDS'], dummies=c_dict['DUMMIES'])
    if with_coords:
        return c_dict, report
    else:
        return report


def scan_struct(options, coords={}, dummies={}):
    """Scans main data structure for structural coordinates
    and dummy atoms, and returns a report with useful statistics
    on possible critical features (holes, clashes, etc)
    """

    def examine_res(res_pair):
        """Finds clashes between pairs of residues"""

        # It is an internal function because it does not
        #  use the typical res-coord format
        r1, r2, cd1, cd2, clash_thr = res_pair
        atnames_r1 = [x for x in cd1]
        atnames_r1.sort()
        atcoords_r1 = np.array([cd1[x] for x in atnames_r1])
        atnames_r2 = [x for x in cd2]
        atnames_r2.sort()
        atcoords_r2 = np.array([cd2[x] for x in atnames_r2])
        at_dX = scipy.spatial.distance.cdist(atcoords_r1, atcoords_r2, 'euclidean')
        clash_arr = np.argwhere(at_dX < clash_thr)
        if clash_arr.size > 0:
            clash_d = at_dX[at_dX < clash_thr]
            return (ch_name1, r1, ch_name2, r2, np.amin(clash_d))
        else:
            return ''

    def CA_data(coords):
        """Returns two lists with the CA coords and labels of all residues"""
        cX = {}
        cXlegend = {}
        for ch in coords:
            if ch not in cX:
                cX[ch] = []
                cXlegend[ch] = []
            for resid_name in coords[ch]:
                cXlegend[ch].append((ch, resid_name))
                if 'CA' in coords[ch][resid_name]:
                    cX[ch].append(coords[ch][resid_name]['CA'])
                else:
                    cX[ch].append((-9999, -9999, -9999))
        return cX, cXlegend


    # C_alpha coords
    cX, cXlegend = CA_data(coords)

    # 1. Check for: UNK, structure holes, disordered residues
    #  A disordered residue is a residue whose index does not follow incremental
    #  numbering (i.e. presence of insertion codes, or index gaps not reflected
    #  in structural holes).
    #  This analysis is agnostic wrt residue order in the coordinate file
    hole_thr = int(options['PARAMETERS'][('', 'hole_thr')])
    holes = {}
    holes_frac = {}
    disordered = {}
    unk_ids = {}
    for ch_name in coords:
        unk_ids[ch_name] = []  # List of UNK resids (not the same list as in parsers! This is generally after-ENC filter)
        holes[ch_name] = []  # List of 2-tuples with beginning and end resids of holes
        holes_frac[ch_name] = 0
        disordered[ch_name] = []  # List of disordered resids
        total_hole = 0
        resid_list = []
        for x, x2 in coords[ch_name]: # x is (<int>resid, <str>inscode), x2 is <str>3-resname
            resid_list.append(x)
            if x2 == "UNK":
                unk_ids[ch_name].append(x)

        # Check for disorder
        ref_res_disorder = -9999
        is_disordered = False
        resid_list = sorted(resid_list) # Does not check for really disordered residues (i.e. resid 120 coming before 31)
        for i, resid in [(ix, x) for ix, x in enumerate(resid_list)][1:]: # resid is (<int>resid, <str>inscode)
            resid_d = resid[0] - resid_list[i-1][0]
            # If there is an index gap
            if resid_d > 1:
                if list(cX[ch_name][i-1]) == [-9999, -9999, -9999]:
                    is_disordered = True
                    ref_res_disorder = resid_list[i-1][0]
                elif list(cX[ch_name][i]) == [-9999, -9999, -9999]:
                    is_disordered = True
                    ref_res_disorder = resid_list[i][0]
                else:
                    CA_dist = np.linalg.norm(np.array(cX[ch_name][i]) - np.array(cX[ch_name][i-1]))
                    # It is a hole if CA_dist > 4 A (this parameter is fixed)
                    if CA_dist > 4.0:
                        # Only index gaps of more than hole_thr get recorded
                        if resid_d > hole_thr:
                            holes[ch_name].append((resid_list[i-1], resid))
                        # ...but all gaps contribute to the total gap estimation
                        total_hole += resid_d
                    # If the index gap is not reflected in a structure hole, then residues are disordered
                    else:
                        is_disordered = True
                        ref_res_disorder = resid_list[i-1][0]
            else:
                is_disordered = True
                ref_res_disorder = resid_list[i-1][0]
            if is_disordered:
                if resid[0] - ref_res_disorder != 1:  # It can still be 0 (different inscode) 
                    disordered[ch_name].append(resid)
                else:
                    is_disordered = False

        # Fraction of holes counted as gaps over resids in the chain
        holes_frac[ch_name] = total_hole / (resid_list[-1][0] - resid_list[0][0] + 1)


    # 2. Checks for clashes
    clash_thr = float(options['PARAMETERS'][('', 'clash_distance')])
    clashes = {}
    chains = sorted([x for x in coords])
    ch_contacts = {}
    bound_chains = set()		
    spheres = clash_spheres()  # (<str>atname, <float>radius)  where atname is the label of the atom that's the center of the sphere
    for ich, ch_name1 in enumerate(chains):
        resids_ch1 = sorted([x for x in coords[ch_name1] if (ch_name1, x) in cXlegend[ch_name1]])
        spheres_ch1 = [spheres[x[1]][1] for x in resids_ch1]
        for ch_name2 in chains[ich:]:
            # Use sphere intersection method to filter out all trivial non-clashes
            dX = scipy.spatial.distance.cdist(cX[ch_name1], cX[ch_name2], 'euclidean')
            resids_ch2 = sorted([x for x in coords[ch_name2] if (ch_name2, x) in cXlegend[ch_name2]])
            spheres_ch2 = [spheres[x[1]][1] for x in resids_ch2]
            spheres_dX = np.array([np.array(spheres_ch2)+x for x in spheres_ch1]) # Matrix of same dimension of dX containing all thresholds
            res_pairs = np.argwhere(dX <= spheres_dX) # Threshold check
            new_res_pairs = []
            for x1, x2 in res_pairs:
               # Only record if clash is between two chains or two non-consecutive residues
               if ch_name1 != ch_name2:
                   new_res_pairs.append((resids_ch1[x1], resids_ch2[x2], coords[ch_name1][resids_ch1[x1]], coords[ch_name2][resids_ch2[x2]], clash_thr))
               else:
                   if x2 >= x1 + 2:
                       new_res_pairs.append((resids_ch1[x1], resids_ch2[x2], coords[ch_name1][resids_ch1[x1]], coords[ch_name2][resids_ch2[x2]], clash_thr))
            res_pairs = sorted(new_res_pairs, key=lambda x: (x[0], x[1]))

            # Refine clash report with the all-atom method examine_res (much slower)
            tmpclashes = map(examine_res, res_pairs)
            clashes[(ch_name1, ch_name2)] = clashes[(ch_name2, ch_name1)] = [x for x in tmpclashes if x != '']

            contact_found = False
            if not clashes[(ch_name1, ch_name2)]:
                # Use sphere intersection method again, this time with an offset of 3 A for catching contacts
                res_pairs_contact = np.argwhere(dX <= spheres_dX + 3)
                res_pairs_contact = [(x[0], x[1], coords[ch_name1][x[0]], coords[ch_name2][x[1]], 4.5) for x in res_pairs if ch_name1 != ch_name2 or x[1][0][0] >= x[0][0][0]+2]
                # We only need to find one contact...
                for x in res_pairs_contact:
                    if examine_res(x) != '':
                        contact_found = True
                        break
            else:
                # A clash is a contact
                contact_found = True

            ch_contacts[(ch_name1, ch_name2)] = contact_found
            ch_contacts[(ch_name2, ch_name1)] = contact_found
			
    report = {
        'clashes' : clashes,
        'clash_frac' : {x : 2 * len(clashes[x]) / (len(coords[x[0]]) + len(coords[x[1]])) for x in clashes},
        'holes' : holes,
        'holes_frac' : holes_frac,
        'disordered' : disordered,
        'chains_in_contact' : ch_contacts,
        'UNK' : unk_ids
    }

    return report


def write_biopython_input(original_filename, biopython_filename, text=[]):
		if (not original_filename) and (not text):
			# INCOMPLETE print_log
			exit(1)
		if (not original_filename) and text:
			orig_text = text
		else:
			orig_text = []
			with open(original_filename) as origf:
				for line in origf:
					orig_text.append(line)

		with open(biopython_filename, 'w') as biopython_file:
			for nl, line in enumerate(orig_text):
				# For this try, remove all HETATM: even good chains could lose a few residues, but we are only interested in
				# counting chains, so it is ok.
				if line.startswith("HETATM"):
					continue
				elif line.startswith("ATOM") and line[17:20] == "UNK":
					biopython_file.write(line[:17] + 'XXX' + line[20:]+'\n')
				else:
					biopython_file.write(line+'\n')


def write_representation_chart(str_data, out_filename):
    with open(out_filename, 'w') as outf:
        for pdbi in sorted([x for x in str_data]):
            if str_data[pdbi]['FROM_OPM']['primary_representation']:
                outf.write("{0}\t".format(pdbi))
                sec_pdbis = [x['pdbid'] for x in str_data[pdbi]['FROM_OPM']\
                    ['secondary_representations']]
                for sec_pdbi in sorted(sec_pdbis):
                    outf.write("{0} ".format(sec_pdbi))
                outf.write("\n")


def make_DSSP(dssp_path, coord_filename, dssp_out_filename, pdbi, c_dict={}, cache_coord_filename='', opm_check=True, leave_original_dssp=True, container=False):
	# Run DSSP
	fnull = open("/dev/null", "w")
	dssp_clist = dssp_path.split()
	p = subprocess.Popen(dssp_clist + ['-i', coord_filename, '-o', dssp_out_filename], stdout=fnull, stderr=fnull)
	p.wait()

	# If no output file is created, return empty
	if not os.path.exists(dssp_out_filename):
		return [], []

	if not c_dict:
		if cache_coord_filename:
			c_dict = retrieve_coords_from_CIF(cache_coord_filename)
		else:
			c_dict, _ = parse_PDB_coords(coord_filename)

	# Parse DSSP file and create a handy dictionary
	# Create a DSSP_segment list with elements such as ('A', 10, 34), where 'A' = helix and 'B' = sheet. DSSP: B, E = beta, H, G, I = alpha.
	dssp_to_ab = {'H':'A', 'G':'A', 'I':'A', 'B':'B', 'E':'B'}
	residues = []
	ss_set = set()
	ss_dict = {}
	Tres = {}
	start_next = False
	with open(dssp_out_filename, 'r') as dssp_file:
		for line in dssp_file:
			if not line:
				continue
			fields = line.split()
			if start_next and fields[1] != '!':
				if line[16].strip() in list(dssp_to_ab.keys()):
					ss_type = dssp_to_ab[line[16]].strip()
					ch_resid = (line[11], (int(line[5:10].strip()), line[10].strip()))
					ss_set.add(ch_resid)
					if leave_original_dssp:
						ss_dict[ch_resid] = line[16].strip()
					else:
						ss_dict[ch_resid] = ss_type
				elif line[16].strip() == "T":
					if line[11] not in Tres:
						Tres[line[11]] = set()
					resid = (int(line[5:10].strip()), line[10].strip())
					Tres[line[11]].add(resid)
	
			if fields[0].strip() == '#':
				start_next = True

	dssp_segments = {}
	dssp_seg = []
	newTres = {}
	for ch in c_dict['COORDS']:
		dssp_segments[ch] = []
		for resid, resname in [x for x in c_dict['COORDS'][ch]]:
			ch_resid = (ch, resid)
			if ch in Tres and resid in Tres[ch]:
				if ch not in newTres:
					newTres[ch] = []
				newTres[ch].append((resid, resname))
			if ch_resid in ss_set:
				if not dssp_seg:
					dssp_seg.append(ss_dict[ch_resid])
					dssp_seg.append((resid, resname))
				elif dssp_seg[0] == ss_dict[ch_resid]:
					dssp_seg.append((resid, resname))
				elif dssp_seg[0] != ss_dict[ch_resid]:
					dssp_segments[ch].append(dssp_seg)
					dssp_seg = []
					dssp_seg.append(ss_dict[ch_resid])
					dssp_seg.append((resid, resname))
			else:
				if dssp_seg:
					dssp_segments[ch].append(dssp_seg)
					dssp_seg = []
	# Residues marked with T are disordered but with a signal on the same H-bond pattern as helices. They can be important
	return dssp_segments, newTres


def parse_PDB_coords(struct_filename, output_filename='', scan=False, biopython_compatible=True, model_to_consider=1, with_UNK=True, options={}, thr_log_status="ERROR"):
    # Checks for altlocs
    # for each altloc residue, if altloc contains the CA, choose the altloc with the best placed CA.
    # If the CAs are equally placed or the CA is not contained in the altloc, look at the occupancies.
    # If the occupancies are equal, choose the first one.

    this_name = parse_PDB_coords.__name__

    record = {
        "unk_res" : set(),
        "unrecognized_res" : set(),
        "mse_res" : set(),
        "het_res" : set(),
        "hetunk_res" : set(),
        "double_record" : set(),
        "altlocs" : {},
        "lines" : None,
        "chains" : [],
        "merged_chains" : set(),
        "pdb_format" : True  # We are reading a PDB file, thus this entry is PDB-compatible
    }

    if not struct_filename:
        if scan:
            return {}, {}
        else:
            return {}

    struct_file = open(struct_filename,'r')
    text = struct_file.read().split('\n')
    struct_file.close()

    bb_atoms = ['N', 'CA', 'C', 'O', 'OXT', 'CB']
    get_rid_list = []

    # Cycle zero: considers only one model
    # It is vital for the good functioning of the first cycle (otherwise there will be double countings due to different models)
    monomodel_text = []
    this_is_it = False
    there_is_no_model = True
    for nl, line in enumerate(text):
        if not line:
            continue

        if line.startswith("MODEL"):
            fields = line.split()
            there_is_no_model = False
            if int(fields[1]) == model_to_consider:
               this_is_it = True
        elif line.startswith("ENDMDL"):
            this_is_it = False

        if line.startswith("ATOM") or line.startswith("HETATM") or line.startswith("TER ") or line.startswith("END ") or line.startswith("ANISOU"):
            if this_is_it or there_is_no_model:
                monomodel_text.append(line)
        else:
            monomodel_text.append(line)
    text = monomodel_text


    # First cycle:
    #  Gets rid of double atoms
    labels = set()
    res_at_with_altlocs = {}
    last_atom = 0
    modified = False
    chains = []
    old_ch = " "
    old_iresid = -1000
    continuity_check = set()
    continuity_hook = {}
    for nl, line in enumerate(text):
        if not line:
            continue
        # Gets rid of double atoms both in ATOM and HETATM
        if line.startswith("ATOM") or line.startswith("HETATM"):
            label = line[12:27]
            ch = line[21]
            s_resid = line[22:26].strip()

            # Set resid to -999 for DU<
            if "DUM" in line:
                resid = -999
            # Eliminate residues with non-numeric resID (without inscode)
            elif not string_is_int(s_resid):
                resname = line[17:20]
                if resname == "DUM" and s_resid == "****" and old_iresid > 0: # Out-of-range resid 
                    resid = old_iresid + 1
                else:
                    get_rid_list.append(nl)
                    continue
            else:
                resid = int(s_resid)

            if old_ch != ch:
                if ch.strip():
                    # If the resID does not reset from one chain to the next, activate continuity check
                    if abs(old_iresid - resid) < 2:
                        continuity_check.add((old_ch, ch))
                        continuity_hook[ch] = nl
                old_ch = ch
            if ch.strip() and old_iresid != resid:
                old_iresid = resid

            # Eliminates double lines
            if label in labels:
                get_rid_list.append(nl)
                modified = True
                continue
            last_atom = nl

        # Eliminate all HETATM which are not: DUM, MSE, HOH, or UNK
        if line.startswith("HETATM"):
            ch = line[21]
            s_resid = line[22:26].strip()
            # Eliminate residues with non-numeric resID (without inscode)
            if not string_is_int(s_resid):
                resname = line[17:20]
                if resname == "DUM" and s_resid == "****" and old_iresid > 0: # Out-of-range resid 
                    resid = old_iresid + 1
                else:
                    get_rid_list.append(nl)
                    continue
            else:
                resid = int(line[22:26].strip())
            if ch.strip() and old_iresid != resid:
                old_iresid = resid

            inscode = line[26].strip()
            res = (resid, inscode)        
            record["het_res"].add((res, ch))
            if (line[17:20]!='MSE' and line[17:20]!='DUM' and line[17:20]!='HOH' and (line[17:20]!='UNK' or not with_UNK)):
                if line[17:20] =='UNK': # UNK residues classified as HETATM
                    record["unk_res"].add((res, ch))
                    record["hetunk_res"].add((res, ch))
                get_rid_list.append(nl)
                modified = True
                continue
            elif line[17:20]=='MSE':
                record["mse_res"].add((res, ch))

        # Eliminate all TER atoms whose residue is X (if UNK option is not there)
        if line.startswith("TER"):
            # Gets rid of strange residue names. If with_UNK = True, only of 'X' names that are not 'UNK'
            # No residuewise statistics are taken here because TER is only just an atom of an ATOM residue
            if from3to1(line[17:20].strip()) == 'X' and (not with_UNK or (with_UNK and line[17:20].strip() != 'UNK')):
                get_rid_list.append(nl)
                modified = True
                continue
            last_atom = nl

        # Eliminate ANISOU, MASTER, CONECT lines
        if line.startswith("ANISOU") or line.startswith("MASTER") or line.startswith("CONECT"):
            get_rid_list.append(nl)
            continue

        if line.startswith("ATOM") or line.startswith("HETATM"):
            ch = line[21]
            resid = -999 if "DUM" in line else int(line[22:26].strip())
            inscode = line[26].strip()
            res = (resid, inscode)
            atname = line[12:16].strip()
            # Gets rid of strange residue names (like UNK), but not of mislabeled DUM (will be transformed into HETATM)
            if from3to1(line[17:20].strip()) == 'X' and line[17:20]!='DUM':
                if line[17:20].strip() == 'UNK':
                    record['unk_res'].add((res, ch))
                    if (not with_UNK) or line[12:16].strip() not in bb_atoms:
                        get_rid_list.append(nl)
                        modified = True
                        continue
                else:
                    # ALL unrecognized atoms are eliminated. This is a difference with mmCIF parser (because here we can't check the entity)
                    if line.startswith("ATOM"):
                        record['unrecognized_res'].add((resid, ch))
                    get_rid_list.append(nl)
                    modified = True
                    continue
            # Records atoms having AltLoc (so they can be recognized even if one AltLoc version has empty AltLoc)
            if line[16:17]!=' ':
                if (res, ch) not in res_at_with_altlocs:
                    res_at_with_altlocs[(res, ch)] = [] #[atname]
                res_at_with_altlocs[(res, ch)].append(atname)
            # Eliminates intertwined chains
            if ch in chains and ch != chains[-1]:
                get_rid_list.append(nl)
            elif ch not in chains:
                chains.append(ch) # Only chains with at least an ATOM record are recorded


    # Continuity check: if chain has changed name but resID has not reset
    #  then check the distance between the C_alphas of the two extrema.
    #  If it is < 5A, the two chains are joined
    if continuity_check:
        confirmed_cont = []
#       continuity_check = merge(continuity_check)
        for s in continuity_check:
            ch1, ch2 = s
            CAold_coords, CA_coords = None, None
            for line in text[continuity_hook[ch2]:continuity_hook[ch2]-21:-1]:
                if (line.startswith("ATOM") or line.startswith("HETATM")) and line[12:16].strip() == "CA":
                    CAold_coords = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
                    break
            for line in text[continuity_hook[ch2]:continuity_hook[ch2]+21]:
                if (line.startswith("ATOM") or line.startswith("HETATM")) and line[12:16].strip() == "CA":
                    CA_coords = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
                    break
            if type(CAold_coords) != type(None) and type(CA_coords) != type(None) and np.linalg.norm(CAold_coords - CA_coords) < 5.0:
                confirmed_cont.append((ch1, ch2))
        continuity_check = merge(confirmed_cont, sorted_lists=True) # If there are 3 consecutive chains, the second AND third will be replaced by the first
        record['merged_chains'] = continuity_check
		
    # Cycle 2
    #  Records occupancies and atom numbers of the various AltLocs, so that the AltLocs that are not chosen can be
    #  readily added to the get_rid_list
    occupancies = {}
    norm_occupancies = {}
    altloc_atoms = {}
    atnums = {}	# All atom numbers corresponding to a given chain+res
    res_with_altlocs = {}
    delete_these_residues = []
    DELappend = delete_these_residues.append
    res_ch_old = ((-9999, 'XXX'), 'XXX')
    ca_exists = False
    chain_list = []
    def VARAPP(arr):
        return arr.append
    def VARADD(arr):
        return arr.add
    for nl, line in enumerate(text):
        if not line:
            continue
        if line.startswith("ATOM"):
            if nl in get_rid_list:
                continue
            if line[17:20]=='DUM':
                continue
            ch = line[21]
            if ch not in occupancies:
                occupancies[ch] = {}
                norm_occupancies[ch] = {}
                altloc_atoms[ch] = {}
            resid = int(line[22:26].strip())
            inscode = line[26].strip()
            res = (resid, inscode)
            atname = line[12:16].strip()
            res_ch = (res, ch)
            try:
                VARAPP(atnums[res_ch])(nl)
            except KeyError:
                atnums[res_ch] = []
                VARAPP(atnums[res_ch])(nl)
            if res_ch in res_at_with_altlocs and atname in res_at_with_altlocs[res_ch]:
                altloc = line[16].strip()
                try:
                    norm_occupancies[ch][altloc] += 1
                except KeyError:
                    occupancies[ch][altloc] = 0
                    norm_occupancies[ch][altloc] = 0
                    altloc_atoms[ch][altloc] = []
                    norm_occupancies[ch][altloc] += 1
                occupancy = line[54:60].strip()
                if occupancy:
                    occupancies[ch][altloc] += float(occupancy)
                norm_occupancies[ch][altloc] += 1
                VARAPP(altloc_atoms[ch][altloc])(nl)
                try:
                    VARADD(res_with_altlocs[res_ch])(altloc)
                except KeyError:
                    res_with_altlocs[res_ch] = set()
                    VARADD(res_with_altlocs[res_ch])(altloc)
            if res_ch != res_ch_old:
                if res_ch_old[1] != 'XXX':
                    if not ca_exists:
                        DELappend(res_ch_old)
                        modified = True
                    ca_exists = False
                bb_list = []
                res_ch_old = res_ch
            if atname == 'CA':
                ca_exists = True
                if ch not in chain_list:
                    chain_list.append(ch)
    if not ca_exists and res_ch_old[1] != 'XXX':
        DELappend(res_ch_old)
        modified = True

    # Chooses the best AltLoc basing on average occupancy throughout the structure
#    print("Chooses best Altloc and gets rid of the others")
    best_altloc = {}
    for ch in occupancies:
        max_occupancy = -1
        for c in occupancies[ch]:
            if norm_occupancies[ch][c] > 0:
                occupancies[ch][c] = occupancies[ch][c] / norm_occupancies[ch][c]
                if max_occupancy < occupancies[ch][c]:
                    max_occupancy = occupancies[ch][c]
                    best_altloc[ch] = c
        record['altlocs'][ch] = {c for c in occupancies[ch]}

    # Gets rid of other AltLocs
    for ch in occupancies:
        for c in occupancies[ch]:
            # AltLocs that are not chosen are added to the get_rid_list
            if c != best_altloc[ch]:
                get_rid_list += altloc_atoms[ch][c]
            # If there are residues with AltLocs that do not contain the chosen AltLoc, they are eliminated as well
            else:
                for cr in res_with_altlocs:
                    if c not in res_with_altlocs[cr]:
                        get_rid_list += atnums[cr]
    for cr in delete_these_residues:
        get_rid_list += atnums[cr]

    # Write new structure (on a string)
    new_text_l = []
    NTappend = new_text_l.append
    LINELEM="{:76s}{:>2s}  \n"
    LINfmt = LINELEM.format
    HETfmt = 'HETATM{0} {1} {2}'.format
    UNKfmt = '{0} UNK{1} {2}'.format
    FREEfmt = '{0} {1} {2}'.format
    ELfmt = '{0:>2}'.format
    ATMET1fmt = 'ATOM  {0} SD  MET{1} {2}'.format
    ATMET2fmt = 'ATOM  {0} MET{1} {2}'.format
    FREELINfmt = '{0}\n'.format

    there_is_END = False
    atoms_begin = False
    continuity_choice = {}
    for i,s in enumerate(continuity_check):
        for c in s:
            continuity_choice[c] = sorted(list(s))[0]
		
    for nl, line in enumerate(text):
        if not line:
            continue
        if nl in get_rid_list:
            continue
        if (line.startswith("ATOM") or line.startswith("TER ")):
            # Enforce continuity between two adjacent chains...
            if line[21] in continuity_choice:
                line = line[:21] + continuity_choice[line[21]] + line[22:]
            # Write lines with their format
            if line.startswith("TER "):
                NTappend(FREELINfmt(line))
            else:
                atoms_begin = True
                atname = line[12:16].strip()
                elem = atname[0]
                if line[17:20]=='DUM':
                    NTappend(LINfmt(HETfmt(line[6:16],line[17:26],line[27:76]),ELfmt(elem)))
                elif from3to1(line[17:20]) == 'X':	# All unrecognized amino acids are transformed into UNK (and only the backbone and CB is kept)	
                    if with_UNK:
                        NTappend(LINfmt(UNKfmt(line[:16],line[20:26],line[27:76]),ELfmt(elem)))
                    else:
                        continue
                else:
                    NTappend(LINfmt(FREEfmt(line[:16],line[17:26],line[27:76]),ELfmt(elem)))
        elif line.startswith("HETATM") and line[17:20]=='MSE':
            atoms_begin = True
            atname = line[12:16].strip()
            if atname == 'SE':
                NTappend(LINfmt(ATMET1fmt(line[6:12],line[20:26],line[27:76]),ELfmt('S')))
            else:
                elem = atname[0]
                NTappend(LINfmt(ATMET2fmt(line[6:16],line[20:26],line[27:76]),ELfmt(elem)))
        elif line.startswith("HETATM") and line[17:20]=='UNK' and with_UNK:
            atoms_begin = True
            atname = line[12:16].strip()
            elem = atname[0]
            NTappend(LINfmt("ATOM  "+UNKfmt(line[6:16],line[20:26],line[27:76]),ELfmt(elem)))
        else:
            NTappend(FREELINfmt(line))
        if line.startswith("END"):
            there_is_END = True
    if not there_is_END:
        NTappend('END\n')
    new_text = "".join(new_text_l)
    # NOTE: from here on, 'new_text' contains the structure

    # Properties of the new PDB file
    new_chains_with_coords = set()
    res_ids = {}
    b_factor = {}
    for nl, line in enumerate(new_text.split('\n')):
        if not line:
            continue
        if (line.startswith("ATOM")):
            if not line[21]:
                ch_name = '_'
            else:
                ch_name = line[21]
            if ch_name not in new_chains_with_coords:
                res_ids[ch_name] = []
                b_factor[ch_name] = []
                bf, bn = 0, 0
                new_chains_with_coords.add(ch_name)

            resid = int(line[22:26])
            inscode = line[26].strip()
            res = (resid, inscode)
            if ch_name in res_ids and (not res_ids[ch_name] or  res_ids[ch_name][-1][0] != resid):
                resname = line[17:20]
                res_ids[ch_name].append((res, resname))
                if bn:
                    b_factor[ch_name].append(bf/bn)
                bf, bn = 0, 0

            if line[60:66].strip():
                bf += float(line[60:66])
                bn += 1

    record['lines'] = len(new_text.split('\n'))
    record['chains'] = sorted(list(new_chains_with_coords))

    pdb_dict = {
        'COORDS' : {},  # (Per chain) nested dictionary <nested tuple>res -> <str>atname -> <3-tuple>coords
        'DUMMIES' : {}, 
        'PDBSEQ' : { ch : [from3to1(x[1]) for x in res_ids[ch]] for ch in new_chains_with_coords }, # (Per chain) List of 1-character residues
        'AVG_BFACTOR' : { ch : b_factor[ch] for ch in new_chains_with_coords },    # (Per chain) List of B factor averages over residues
        'RESIDS' : { ch : [x[0] for x in res_ids[ch]] for ch in new_chains_with_coords }            # (Per chain) List of pairs (<int>resid, <str>inscode)
    }

    # Scan results are integrated in record data structure
    if scan:
        if not [line.startswith("ATOM") for line in new_text.split()]:
            noat_msg = ('FATAL', this_name, "No ATOM lines after filtering!")
            print_log(noat_msg)
        if not options:
            print("scan option can only be used with options data structure!")
            exit(1)
        c_dict, report = scan_struct_from_PDB_text(options, text=new_text, with_coords=True)
        for r in report:
            record[r] = report[r]
    else:
        c_dict = read_coord_from_PDB_text(new_text)

    pdb_dict['COORDS'] = c_dict['COORDS']
    pdb_dict['DUMMIES'] = c_dict['DUMMIES']

    return pdb_dict, record


def simple_download(url, local_filename):
	try:
		response = urllib.request.urlopen(url)
	except:
		return False

	local_file = open(local_filename, 'wb')
	shutil.copyfileobj(response, local_file)
	local_file.close()
	return True


def reiterated_simple_download(url, local_filename, thr_log_status="ERROR"):
	this_name = reiterated_simple_download.__name__

	max_r = 0.5#3
	error_found = False
	downloaded = False
	reiterations = 0
	while (error_found or not downloaded):
		if reiterations >= max_r:
			return False
		if reiterations > 0:
			time.sleep(reiterations**2)
			try_msg = ('NOTICE', this_name, "Connect to {0} \t Try number: {1}".format(url, reiterations))
			print_log(try_msg)
		downloaded = simple_download(url, local_filename)
		if downloaded:
			with codecs.open(local_filename, "r",encoding='utf-8', errors='ignore') as local_file:
				for line in local_file:
					if "Error: Could not connect to database. Please try again later." in line:
						error_found = True
		reiterations += 0.25
	if error_found:
		return False
	return True


def gzip_download(url, local_filename):
	try:
		response = urllib.request.urlopen(url)
	except:
		print("gzip_download not done!", local_filename)
		return False

	with gzip.GzipFile(fileobj=response) as uncompressed, open(local_filename, 'wb') as local_file:
		shutil.copyfileobj(uncompressed, local_file)
		local_file.close()
	return True


def reiterated_gzip_download(url, local_filename):
	this_name = 'reiterated_gzip_download'
	max_r = 0.5#3
	error_found = False
	downloaded = False
	reiterations = 0
	while (error_found or not downloaded):
		if reiterations >= max_r:
			return False
		if reiterations > 0:
			time.sleep(reiterations**2)
		downloaded = gzip_download(url, local_filename)
		if downloaded:
			with codecs.open(local_filename, "r",encoding='utf-8', errors='ignore') as local_file:
				for line in local_file:
					if "Error: Could not connect to database. Please try again later." in line:
						error_found = True
		reiterations += 0.25
	if error_found:
		return False
	return True


def from3to1(resname, include_MSEs=False):
	f3t1 = {'ALA' : 'A',
	        'ARG' : 'R',
	        'ASN' : 'N',
	        'ASP' : 'D',
	        'CYS' : 'C',
	        'GLN' : 'Q',
	        'GLU' : 'E',
	        'GLY' : 'G',
	        'HIS' : 'H',
	        'ILE' : 'I',
	        'LEU' : 'L',
	        'LYS' : 'K',
	        'MET' : 'M',
	        'PHE' : 'F',
	        'PRO' : 'P',
	        'SER' : 'S',
	        'THR' : 'T',
	        'TRP' : 'W',
	        'TYR' : 'Y',
	        'VAL' : 'V'}

	if resname in list(f3t1.keys()):
		return f3t1[resname]
	elif include_MSEs and resname == 'MSE':
		return 'M'
	else:
		return 'X'

def create_whole_png_text(pdb_filename, fig_filename):
	text = "".join(("delete all\n",
	                "reset\n",
	                "\n",
	                "load {0}\n".format(pdb_filename),
	                "viewport 350,350\n",
	                "hide\n",
	                "bg_color white\n",
	                "\n",
	                "util.cbc\n",
	                "color blue, (resn DUM and name O)\n",
	                "color blue, (resn DUM and name N)\n",
	                "set antialias=1\n",
	                "set line_smooth=1.00000\n",
	                "set depth_cue=1\n",
	                "set specular=1.00000\n",
	                "set surface_quality=1\n",
	                "set stick_quality=25\n",
	                "set sphere_quality=2\n",
	                "set cartoon_sampling=14\n",
	                "set ribbon_sampling=10\n",
	                "set ray_trace_fog=1\n",
	                "\n",
	                "rotate [1,0,0], angle=-90\n",
	                "show cartoon\n",
	                "set sphere_transparency=0.75\n",
	                "show spheres, (resname DUM)\n",
	                "\n",
	                "zoom\n",
	                "png {0}\n".format(fig_filename),
		        "\n"
		        "quit\n"))
	return text


def test_parse_PDB_coords(str_list=[]):
	if not str_list:
		str_list = ["4v8k", "5zf0"]

	test_folder = "./test_parser/"
	if not os.path.exists(test_folder):
		os.mkdir(test_folder)

	for pdbi in str_list:
		struct_filename = test_folder + pdbi + '.pdb'


def merge(lsts, sorted_lists=False):
    """From a list of lists, merges lists with common elements
    """

    sets = [set(lst) for lst in lsts if lst]
    merged = 1
    while merged:
        merged = 0
        results = []
        while sets:
            common, rest = sets[0], sets[1:]
            sets = []
            for x in rest:
                if x.isdisjoint(common):
                    sets.append(x)
                else:
                    merged = 1
                    common |= x
            results.append(frozenset([i for i in common]))
        sets = results
    if sorted_lists:
        main_list = []
        for s in sets:
            main_list.append(sorted(list(s)))
        return sorted(main_list)
    return sets


def find_gcd(l):
    x = reduce(math.gcd, l)
    return x


def reduce_chains(str_data):
    """Applies redundant_chains selection to each pdbi in str_data, removing those chains.
    Also removes each secpdbi in secondary_representations with same chains and same
    redundant_chains.
    """
    rm_set = set()
    rmsec_set = set()
    lock_set = set()
    all_set = set()
    for pdbi in str_data:
        all_set |= {pdbi + "_" + k for k in str_data[pdbi]['ENCOMPASS']['structure']['kchains']}
    isecreps = {}
    for pdbi in str_data:
        if pdbi in rmsec_set:
            continue
        for k in str_data[pdbi]['ENCOMPASS']['structure']['redundant_chains']:
            rm_set.add(pdbi + "_" + k)
        for secpdbi in [x['pdbid'] for x in str_data[pdbi]['FROM_OPM']['secondary_representations']]:
            if secpdbi in lock_set:
                continue
            test = (str_data[pdbi]['ENCOMPASS']['structure']['kchains'] == str_data[secpdbi]['ENCOMPASS']['structure']['kchains'])
            if test:
                str_data[pdbi]['ENCOMPASS']['structure']['redundant_chains']
                if str_data[pdbi]['ENCOMPASS']['structure']['redundant_chains'] == str_data[secpdbi]['ENCOMPASS']['structure']['redundant_chains']:
                    lock_set.add(pdbi)
                    rmsec_set.add(secpdbi)
                    rm_set |= {secpdbi + "_" + k for k in str_data[secpdbi]['ENCOMPASS']['structure']['kchains']}
                else:
                    pass

    return sorted(list(all_set - rm_set))


def colorblind_palette_OLD():
    d = {
        'cbf_dark_teal' : [1,   110, 130, 1],
        'cbf_purple'    : [125,  45, 145, 1],
        'cbf_softblue'  : [47,   94, 171, 1],
        'cbf_skyblue'   : [68,  152, 211, 1],
        'cbf_violet'    : [205, 133, 185, 1],
        'cbf_cyan'      : [70,  195, 208, 1],
        'cbf_redrusset' : [170,  29, 63,  1],
        'cbf_peach'     : [244, 119, 82,  1],
        'cbf_green'     : [25,  179, 90,  1],
        'cbf_lime'      : [237, 232, 59,  1],
        'cbf_peagreen'  : [171, 211, 122, 1],
        'cbf_wheat'     : [249, 229, 190, 1]
    }
    l = [
        'cbf_peach',
        'cbf_dark_teal',
        'cbf_lime',
        'cbf_purple',
        'cbf_softblue',
        'cbf_skyblue',
        'cbf_violet',
        'cbf_cyan',
        'cbf_redrusset',
        'cbf_green',
        'cbf_peagreen',
        'cbf_wheat'
    ]
    for i in d:
        d[i][0:3] = [x/256 for x in d[i][0:3]]
    return l, d


def colorblind_palette(n=0):
    import cbf_colors
    if n<9:
        l = [
            "cbf_orange", 
            "cbf_sky_blue",
            "cbf_bluish_green",
            "cbf_yellow",
            "cbf_blue",
            "cbf_vermillion",
            "cbf_reddish_purple"
        ]
    else:
        l = [
            'cbf_peach',
            'cbf_skyblue',
            'cbf_green',
            'cbf_wheat',
            'cbf_softblue',
            'cbf_redrusset',
            'cbf_purple',
            'cbf_lime',
            'cbf_dark_teal',
            'cbf_cyan',
            'cbf_violet',
            'cbf_peagreen'
        ]
    return l


def init_pairwise_alignment():
    aligner = Align.PairwiseAligner()
    aligner.mode = 'global'
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -2.0
    aligner.extend_gap_score = -1.0
    aligner.target_end_gap_score = 0.0
    aligner.query_end_gap_score = 0.0

    return aligner

def gap_function(x, y):  # x is gap position in seq, y is gap length
	if y == 0:  # No gap
		return 0
	elif y == 1:  # Gap open penalty
		return -2
	return - (2 + y/4.0 + math.log(y)/2.0)


def write_OPM_representation_chart(ids, secrep_d, OPMreprchart):
    """Writes an ORDERED list of all primary and secondary entries used
    Secondary entries always follow their corresponding primary entry
    """
    with open(OPMreprchart, 'w') as f:
        f.write("#OPM {0}\n".format(datetime.date.today()))
        for i, pdbi, is_temp in ids:
            t = "temporary" if is_temp else ""
            f.write(f"{i}\t{pdbi}\t  primary\t{t}\n")
            if pdbi in secrep_d:
                for sec_pdbi in secrep_d[pdbi]:
                    f.write(f"?\t{sec_pdbi}\tsecondary\t\n")

def write_OPM_representation_chart_2(str_data, OPMreprchart):
    """Writes an ORDERED list of all primary and secondary entries used
    Secondary entries always follow their corresponding primary entry
    """
    seclist = []
    with open(OPMreprchart, 'w') as f:
        f.write("#OPM {0}\n".format(datetime.date.today()))
        for i, pdbi in enumerate(str_data):
            if str_data[pdbi]['FROM_OPM']['primary_representation']:
                f.write(f"{i}\t{pdbi}\t  primary\t\n")
                for sec_pdbi in str_data[pdbi]['FROM_OPM']['secondary_representations']:
                    if sec_pdbi not in seclist:
                        f.write("?\t{0}\tsecondary\t\n".format(sec_pdbi['pdbid']))
                        seclist.append(sec_pdbi)


def read_OPM_representation_chart(OPMreprchart):
    """Reads the representation chart, ASSUMING it's ORDERED and thus
    that secondary entries refer to the preceding primary entry
    """
    ids = []
    secrep_d = {}
    with open(OPMreprchart) as f:
        last_primary = ""
        for line in f:
            if line.startswith("#"):
                continue
            i, pdbi, order, t = [x.strip() for x in line.split("\t")]
            if order == "primary":
                last_primary = pdbi
                if t == "temporary":
                    is_temp = True
                elif not t:
                    is_temp = False
                ids.append((int(i), pdbi, is_temp)) # (<int>ID, <str>PDBID, <bool>is_temp)
                secrep_d[pdbi] = []
            elif order == "secondary":
                secrep_d[last_primary].append(pdbi)
    # [(<int>ID, <str>PDBID, <bool>is_temp), ...], <str>PDBID->[<str>PDBID, ...]
    return ids, secrep_d



def write_OPM_representation_chart(ids, secrep_d, OPMreprchart):
    """Writes an ORDERED list of all primary and secondary entries used
    Secondary entries always follow their corresponding primary entry
    """
    with open(OPMreprchart, 'w') as f:
        f.write("#OPM {0}\n".format(datetime.date.today()))
        for i, pdbi, is_temp in ids:
            t = "temporary" if is_temp else ""
            f.write(f"{i}\t{pdbi}\t  primary\t{t}\n")
            if pdbi in secrep_d:
                for sec_pdbi in secrep_d[pdbi]:
                    f.write(f"?\t{sec_pdbi}\tsecondary\t\n")

def write_OPM_representation_chart_2(str_data, OPMreprchart):
    """Writes an ORDERED list of all primary and secondary entries used
    Secondary entries always follow their corresponding primary entry
    """
    seclist = []
    with open(OPMreprchart, 'w') as f:
        f.write("#OPM {0}\n".format(datetime.date.today()))
        for i, pdbi in enumerate(str_data):
            if str_data[pdbi]['FROM_OPM']['primary_representation']:
                f.write(f"{i}\t{pdbi}\t  primary\t\n")
                for sec_pdbi in str_data[pdbi]['FROM_OPM']['secondary_representations']:
                    if sec_pdbi not in seclist:
                        f.write("?\t{0}\tsecondary\t\n".format(sec_pdbi['pdbid']))
                        seclist.append(sec_pdbi)


def read_OPM_representation_chart(OPMreprchart):
    """Reads the representation chart, ASSUMING it's ORDERED and thus
    that secondary entries refer to the preceding primary entry
    """
    ids = []
    secrep_d = {}
    with open(OPMreprchart) as f:
        last_primary = ""
        for line in f:
            if line.startswith("#"):
                continue
            i, pdbi, order, t = [x.strip() for x in line.split("\t")]
            if order == "primary":
                last_primary = pdbi
                if t == "temporary":
                    is_temp = True
                elif not t:
                    is_temp = False
                ids.append((int(i), pdbi, is_temp)) # (<int>ID, <str>PDBID, <bool>is_temp)
                secrep_d[pdbi] = []
            elif order == "secondary":
                secrep_d[last_primary].append(pdbi)
    # [(<int>ID, <str>PDBID, <bool>is_temp), ...], <str>PDBID->[<str>PDBID, ...]
    return ids, secrep_d


### TEST MODULE ---------------------------------------------------------

def test_merge():
    l = {}
    l[1] = [[0, 1], [1, 2], [2, 3], [1, 4], [5, 6]]
    l[2] = [[0], [1, 2, 3, 4, 5], [0, 1]]
    l[3] = [[0, 2, 4], [1, 3, 5], []]
    for k in l:
        print(k, l[k])
        print(merge(l[k]))


def test_gcd():
    print(find_gcd([1,3,5,7]))
    print(find_gcd([2,4,6,8]))
    print(find_gcd([18, 24, 48]))
    print(find_gcd([2,6,8,2,6,6,10,6,6,14,12]))


def test_read_oper_expression():
    mtx_ind = [str(x) for x in range(1,29)]
    print(read_oper_expression("(1-2)", mtx_ind))
    print(read_oper_expression("(1-20)", mtx_ind))
    print(read_oper_expression("1,3,5-15", mtx_ind))
    print(read_oper_expression("(1-20)(21-25)(26,28)", mtx_ind))


def test_mmCIF_format_parser():
    pdbi_list = ['3dh4']
    for pdbi in pdbi_list:
        print('PDB_CODE', pdbi)
        cif_filename = "{0}.cif".format(pdbi.upper())
        if not os.path.exists(cif_filename):
            subprocess.run(["wget", "https://files.rcsb.org/view/{0}".format(cif_filename)])
        all_loops, loop_keys = parse_mmCIF_format(cif_filename)
        for i, loop in enumerate(all_loops):
            print("LOOP", i)
            for line in loop:
                print(",\t".join(line))
        if pdbi_list.index(pdbi) != len(pdbi_list)-1:
            print('\n#################################################################\n\n')


def test_header_parsers():
    pdbi_list = ['3dh4', '1bcc']
    for pdbi in pdbi_list:
        print('PDB_CODE', pdbi)
        cif_filename = "{0}.cif".format(pdbi.upper())
        pdb_filename = "{0}.pdb".format(pdbi.upper())
        if not os.path.exists(cif_filename):
            subprocess.run(["wget", "https://files.rcsb.org/view/{0}".format(cif_filename)])
        if not os.path.exists(pdb_filename):
            subprocess.run(["wget", "https://files.rcsb.org/view/{0}".format(pdb_filename)])
        print('\n########### mmCIF\n')
        all_loops, loop_keys = parse_mmCIF_format(cif_filename)
        print(parse_mmCIF_header(all_loops, loop_keys))
        print('\n########### PDB\n')
        res, rep = parse_PDB_header(pdb_filename)
        print(res.show_dictionary(quiet=True), rep)
        if pdbi_list.index(pdbi) != len(pdbi_list)-1:
            print('\n#################################################################\n\n')


def test_coord_parsers():
    pdbi_list = ['3dh4', '1bcc', '6bsb']
    for pdbi in pdbi_list:
        print('PDB_CODE', pdbi)
        cif_filename = "{0}.cif".format(pdbi.upper())
        pdb_filename = "{0}.pdb".format(pdbi.upper())
        if not os.path.exists(cif_filename):
            subprocess.run(["wget", "https://files.rcsb.org/view/{0}".format(cif_filename)])
        if not os.path.exists(pdb_filename):
            subprocess.run(["wget", "https://files.rcsb.org/view/{0}".format(pdb_filename)])
        print('\n########### mmCIF\n')
        all_loops, loop_keys = parse_mmCIF_format(cif_filename)
        res, rep = parse_mmCIF_coords(all_loops, loop_keys)
        for k in sorted(rep):
            print(k)
            if type(rep[k]) is set:
                print(sorted(list(rep[k])))
            else:
                print(rep[k])
        print('\n########### PDB\n')
        res, rep = parse_PDB_coords(pdb_filename)
        for k in sorted(rep):
            print(k)
            if type(rep[k]) is set:
                print(sorted(list(rep[k])))
            else:
                print(rep[k])
        
        if pdbi_list.index(pdbi) != len(pdbi_list)-1:
            print('\n#################################################################\n\n')


if __name__ == "__main__":
    config.log_filename = 'DEBUG.txt'
    test_roe()
