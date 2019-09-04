#!/usr/bin/env python3

import numpy as np

def read(fname):
	
	"""
	This module reads all of the information in an NWChem .movecs file,
	in unformatted FORTRAN binary format, and returns it as a python
	object with properties corresponding to each entry in the .movecs
	file.
	The only differences are that the lengths of strings in the header are 
	not stored, since these can be obtained from the strings themselves,
	that the object is labeled with the source file name, and that the 
	Python array indices of the HOMO orbital(s) are included as an extra
	entry, for convenience.
	The module only depends on Numerical Python (NumPy).
	Data types are respected throughout, with 64-bit floats, and 32- / 64-bit
	integers as in the .movecs file. Compiled-in FORTRAN integer types are
	automatically detected.
	
	Input:  the NWChem .movecs file name
	Output: the single object containing all information
	
	Erik R. McNellis 2017
	"""
	
	# Class of object containing all information in an NWChem .movecs file
	class movecs_info:
		def __init__(self, fname):
			self.sourcefile = fname
			self.geomsum    = ''
			self.basissum   = ''
			self.bqsum      = ''
			self.scftype20  = ''
			self.date       = ''
			self.scftype    = ''
			self.title      = ''
			self.basisname  = ''
			self.Etot       = 0.0
			self.Enrep      = 0.0
	
		# Function initializes arrays when dimensions are known
		def init_arrays(self,s,nbf,nmos):
			self.spins      = s
			self.nbf        = nbf
			self.nmos       = nmos
			# Copy 'nmos' to homo_index for correct size and type
			self.homo_index	= nmos.copy()
			self.occa       = np.empty([nbf], dtype=np.float64)
			self.evalsa     = np.empty([nbf], dtype=np.float64)
			self.psia       = np.empty([nbf,nmos[0]], dtype=np.float64)
			if (s == 2):
				self.occb   = np.empty([nbf], dtype=np.float64)
				self.evalsb = np.empty([nbf], dtype=np.float64)
				self.psib   = np.empty([nbf,nmos[1]], dtype=np.float64)
			
	# Read string of fixed length
	def read_str(f):
		l = np.fromfile(f, dtype=np.int32, count=1)[0]
		s = np.fromfile(f, dtype=np.dtype('S' + str(l)), count=1)[0]
		f.seek(4,1)
		return s.decode()
	
	# Skips a record without reading
	def skip_rec(f):
		c = np.fromfile(f, dtype=np.int32, count=1)[0]
		f.seek(c + 4,1)
		
	# Reads array of integers from file, with variable integer type
	def read_ints(f,dt):
		r = np.fromfile(f, dtype=np.int32, count=1)[0]
		r = int(r / np.dtype(dt).itemsize)
		ints = np.fromfile(f, dtype=dt, count=r)
		f.seek(4,1)
		return ints
	
	# Reads vector of 64-bit floats from file
	def read_float_vec(f,n):
		f.seek(4,1)
		v = np.fromfile(f, dtype=np.float64, count=n)
		f.seek(4,1)
		return v
	
	# Fills a matrix with values, avoids excess memory by using
	# only the memory already allocated to the matrix
	def fill_float_matrix(f,a,m,n):
		for i in range(n): a[:,i] = read_float_vec(f,m)
	
	# Initialize object
	mov = movecs_info(fname)
	# Open file in read-only binary mode
	f = open(fname, 'rb')
	# Read strings in header, skipping past the integer lengths
	# of the redundant scftype and title strings
	header        = read_str(f)
	mov.geomsum   = header[0:32]
	mov.basissum  = header[32:64]
	mov.bqsum     = header[64:96]
	mov.scftype20 = header[96:116]
	mov.date      = header[116:142]
	mov.scftype   = read_str(f); skip_rec(f)
	mov.title     = read_str(f); skip_rec(f)
	mov.basisname = read_str(f)
	
	# Now we're at the number of spins entry, which is the first with
	# variable integer type (4 or 8 bytes) - determine that type and rewind
	l = np.fromfile(f, dtype=np.int32, count=1)[0]
	itype = np.int32
	if ( l == 8 ): itype = np.int64
	f.seek(-4,1)
	
	# Read number of spins, basis functions and MOs with correct integer type
	s    = read_ints(f,itype)[0]
	nbf  = read_ints(f,itype)[0]
	nmos = read_ints(f,itype)
	# Initialize arrays depending on this information
	mov.init_arrays(s,nbf,nmos)
	
	# Now read orbital occupations, -eigenvalues and MOs for the alpha spin,
	# determine the array index of the HOMO orbital
	mov.occa          = read_float_vec(f,nbf)
	mov.homo_index[0] = np.where(mov.occa > np.float64(0.0))[0][-1]
	mov.evalsa        = read_float_vec(f,nbf)
	fill_float_matrix(f,mov.psia,nbf,nmos[0])
	
	# If there's a second spin channel, repeat for beta spin
	if (s == 2):
		mov.occb          = read_float_vec(f,nbf)
		mov.homo_index[1] = np.where(mov.occb > np.float64(0.0))[0][-1]
		mov.evalsb        = read_float_vec(f,nbf)
		fill_float_matrix(f,mov.psib,nbf,nmos[1])
	
	# Finally, read total- and nuclear repulsion energy, close file, exit
	E = read_float_vec(f,2)
	mov.Etot  = E[0]
	mov.Enrep = E[1]
	f.close()
	
	return mov

