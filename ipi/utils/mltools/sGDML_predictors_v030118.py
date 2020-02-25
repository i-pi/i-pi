#!/usr/bin/python

# GDL Machine Learning Model
# Author: Stefan Chmiela (stefan@chmiela.com)

import sys
import scipy.spatial.distance
import scipy.optimize
from scipy.linalg import norm
import numpy as np
import multiprocessing as mp
import ctypes

import time

VERSION = 20118
ENABLE_FAST_MATCHING = True

def shared_array(arr,typ):
	shpe = arr.shape
	arr_base = mp.Array(typ, arr.flatten(), lock=False)
	arr = np.ctypeslib.as_array(arr_base)
	arr = arr.reshape(shpe)
	return arr

def shared_value(val,typ):
	val_base = mp.Value(typ, val, lock=False)
	val = val_base.value
	return val

def find_optimal_permutation_worker(params):

	global trainPts_D,train_eig_abs
	r,z,v_abs,pdist,start,stop = params

	best_score = float("inf")
	best_perm = range(0,r.shape[0])
	for i in range(start,stop):

		# Cost for matching each atom pair.
		UU = train_eig_abs[i,:,:].dot(v_abs.T)
		cost = np.amax(UU) - UU

		# Find assignment with minimal cost.
		_,col_idx = scipy.optimize.linear_sum_assignment(cost)

		# Generate input descriptor.
		pdist_perm = pdist[:,col_idx]
		pdist_perm = pdist_perm[col_idx,:]
		#r_desc = self._r_to_desc(r,pdist_perm) # r does not need to be permuted because it not used in the funciton
		r_desc  = 1 / pdist_perm[np.tril_indices(r.shape[0],-1)] # HACK

		cdist = scipy.spatial.distance.cdist(r_desc.reshape(-1,1).T,trainPts_D,'euclidean')
		score = np.sum(cdist)

		if score < best_score and not np.sum(abs(z[col_idx] - z)) != 0: # matched_wrong_species
			best_score = score
			best_perm = col_idx

	return (best_score,best_perm)

def predict_worker(params):

	global trainPts_D,r_d_desc_alpha,train_perms_lin,sig
	r_desc,r_d_desc,start,stop = params

	n_perms = train_perms_lin.shape[0] / trainPts_D.shape[1]

	mat52_base_fact = 5/(3*sig**3)
	diag_scale_fact = 5/sig
	sqrt5 = np.sqrt(5)

	# Predict forces and energy.
	E = 0.0
	F = np.zeros((r_d_desc.shape[1],))

	for j in range(start,stop):

		# Create permutated variants of 'trainPts_D_j' and 'r_d_desc_alpha_j'.
		rj_desc_perms = np.reshape(np.tile(trainPts_D[j,:], n_perms)[train_perms_lin], (n_perms,-1), order='F') # 6,36
		rj_d_desc_alpha_perms = np.reshape(np.tile(r_d_desc_alpha[j,:], n_perms)[train_perms_lin], (n_perms,-1), order='F') # 6,36

		diff_ab_perms = r_desc - rj_desc_perms # 6,36
		norm_ab_perms = sqrt5 * np.linalg.norm(diff_ab_perms, axis=1) # 6,

		mat52_base = np.exp(-norm_ab_perms / sig) * mat52_base_fact # 6,
		a_x2 = np.einsum('ij,ij->i', diff_ab_perms, rj_d_desc_alpha_perms) # 6,

		F += np.einsum('ji,j->i', diff_ab_perms.dot(r_d_desc), a_x2 * mat52_base) * diag_scale_fact # 27,

		mat52_base *= (norm_ab_perms + sig) # 6,

		F += np.einsum('ij,i', rj_d_desc_alpha_perms.dot(r_d_desc), -mat52_base)
		E += np.sum(a_x2 * mat52_base)

	return np.append(F, E)

def predict_worker_batch(params):

	global trainPts_D,r_d_desc_alpha,train_perms_lin,sig,b_size
	r_desc,r_d_desc,wkr_start,wkr_stop = params

	n_perms = train_perms_lin.shape[0] / trainPts_D.shape[1]

	mat52_base_fact = 5/(3*sig**3)
	diag_scale_fact = 5/sig
	sqrt5 = np.sqrt(5)

	# Predict forces and energy.
	E = 0.0
	F = np.zeros((r_d_desc.shape[1],))

	#worker_len = wkr_stop - wkr_start
	#b_size = 125

	b_start = wkr_start
	for b_stop in range(wkr_start+b_size,wkr_stop,b_size) + [wkr_stop]:

		b_len = b_stop - b_start

		rj_desc_perms = np.reshape(np.tile(trainPts_D[b_start:b_stop,:], n_perms)[:,train_perms_lin], (b_len*n_perms,-1), order='F') # b_len,6,36
		rj_d_desc_alpha_perms = np.reshape(np.tile(r_d_desc_alpha[b_start:b_stop,:], n_perms)[:,train_perms_lin], (b_len*n_perms,-1), order='F') # b_len,6,36

		#diff_ab_perms = r_desc - rj_desc_perms # b_len*6,36
		np.subtract(r_desc, rj_desc_perms, out=rj_desc_perms); diff_ab_perms = rj_desc_perms # b_len*6,36
		norm_ab_perms = sqrt5 * np.linalg.norm(diff_ab_perms, axis=1) # b_len*6,

		mat52_base = np.exp(-norm_ab_perms / sig) * mat52_base_fact # b_len*6,
		a_x2 = np.einsum('ij,ij->i', diff_ab_perms, rj_d_desc_alpha_perms) # b_len*6,

		F += np.einsum('ji,j->i', diff_ab_perms.dot(r_d_desc), a_x2 * mat52_base) * diag_scale_fact # 27,

		mat52_base *= (norm_ab_perms + sig) # b_len*6,

		F += np.einsum('ij,i', rj_d_desc_alpha_perms.dot(r_d_desc), -mat52_base)
		E += np.sum(a_x2 * mat52_base)

		b_start = b_stop

	return np.append(F, E)

def predict_worker_cached(params):

	global trainPts_D,r_d_desc_alpha,R_desc_perms,R_d_desc_alpha_perms,sig,b_size
	r_desc,r_d_desc,wkr_start,wkr_stop = params

	n_perms = train_perms_lin.shape[0] / trainPts_D.shape[1]

	mat52_base_fact = 5/(3*sig**3)
	diag_scale_fact = 5/sig
	sqrt5 = np.sqrt(5)

	# Predict forces and energy.
	E = 0.0
	F = np.zeros((r_d_desc.shape[1],))

	wkr_start *= n_perms
	wkr_stop *= n_perms

	b_start = wkr_start
	for b_stop in range(wkr_start+b_size*n_perms,wkr_stop,b_size*n_perms) + [wkr_stop]:

		rj_desc_perms = R_desc_perms[b_start:b_stop,:]
		rj_d_desc_alpha_perms = R_d_desc_alpha_perms[b_start:b_stop,:]

		diff_ab_perms = r_desc - rj_desc_perms # b_len*6,36
		#np.subtract(r_desc, rj_desc_perms, out=rj_desc_perms); diff_ab_perms = rj_desc_perms # b_len*6,36
		norm_ab_perms = sqrt5 * np.linalg.norm(diff_ab_perms, axis=1) # b_len*6,

		mat52_base = np.exp(-norm_ab_perms / sig) * mat52_base_fact # b_len*6,
		a_x2 = np.einsum('ij,ij->i', diff_ab_perms, rj_d_desc_alpha_perms) # b_len*6,

		F += np.einsum('ji,j->i', diff_ab_perms.dot(r_d_desc), a_x2 * mat52_base) * diag_scale_fact # 27,

		mat52_base *= (norm_ab_perms + sig) # b_len*6,

		F += np.einsum('ij,i', rj_d_desc_alpha_perms.dot(r_d_desc), -mat52_base)
		E += np.sum(a_x2 * mat52_base)

		b_start = b_stop

	return np.append(F, E)

def predict_worker_new(params):

	global trainPts_D,r_d_desc_alpha,train_perms_lin,sig
	R_desc,R_d_desc,start,stop = params

	n_r = R_desc.shape[1]
	n_perms = train_perms_lin.shape[0] / trainPts_D.shape[1]

	mat52_base_fact = 5/(3*sig**3)
	diag_scale_fact = 5/sig
	sqrt5 = np.sqrt(5)

	# Predict forces and energy.
	E = np.zeros((n_r,))
	F = np.zeros((R_d_desc.shape[1], n_r))

	for j in range(start,stop):

		# Create permutated variants of 'trainPts_D_j' and 'r_d_desc_alpha_j'.
		rj_desc_perms = np.reshape(np.tile(trainPts_D[j,:], n_perms)[train_perms_lin], (n_perms,-1), order='F') # 6,36
		rj_d_desc_alpha_perms = np.reshape(np.tile(r_d_desc_alpha[j,:], n_perms)[train_perms_lin], (n_perms,-1), order='F') # 6,36
		diff_ab_perms = R_desc - np.tile(rj_desc_perms[:,:,None],n_r) # 6,36,2
		norm_ab_perms = sqrt5 * np.linalg.norm(diff_ab_perms, axis=1) # 6,2
		mat52_base = np.exp(-norm_ab_perms / sig) * mat52_base_fact # 6,2
		a_x2 = np.einsum('ijk,ij->ik', diff_ab_perms, rj_d_desc_alpha_perms) # 6,2
		F += np.einsum('jik,jk->ik', np.einsum('ijl,jkl->ikl', diff_ab_perms, R_d_desc), a_x2 * mat52_base) * diag_scale_fact # 27,2
		mat52_base *= (norm_ab_perms + sig) # 6,2
		F += np.einsum('ijk,ik->jk', np.einsum('ij,jkl->ikl', rj_d_desc_alpha_perms, R_d_desc), -mat52_base)
		E += np.einsum('ij,ij->j', a_x2, mat52_base)


	return (E,F)


class GDL_Predictor:

	def __init__(self,model,batch_size=125,num_workers=-1):

		global trainPts_D,r_d_desc_alpha,R_desc_perms,R_d_desc_alpha_perms,train_perms_lin,train_eig_abs,sig,b_size

		# Batch size (number of training samples summed up in prediction process) that a worker processes at once.
		b_size = shared_value(batch_size, ctypes.c_int)

		trainPts_D = shared_array(model['trainPts_D'].T,ctypes.c_double)
		r_d_desc_alpha = shared_array(np.squeeze(model['r_d_desc_alpha']),ctypes.c_double)

		self.trainPts_D = trainPts_D
		self.r_d_desc_alpha = r_d_desc_alpha

		sig = shared_value(model['sig'],ctypes.c_double)

		self.c = model['c']
		self.n_train = trainPts_D.shape[0]

		# Precompute batch permutations.
		self.has_matching_support = 'has_matching_support' in model
		if self.has_matching_support:
			train_eig_abs = shared_array(abs(model['train_eig']),ctypes.c_double)
			train_perms = model['train_perms']
		else:
			train_perms = np.arange(0,trainPts_D.shape[1]).reshape(1,-1) # HACK: non-matching models lack 'train_perms'
		perm_offsets = np.arange(train_perms.shape[0])[:,None] * train_perms.shape[1]
		train_perms_lin = shared_array((train_perms + perm_offsets).flatten('F'),ctypes.c_int)

		# Precompute everything permutated training descriptors and its first derivatives multiplied with the coefficients (only needed for cached variant).
		n_perms = perm_offsets.shape[0]
		R_desc_perms = shared_array(np.reshape(np.tile(trainPts_D, n_perms)[:,train_perms_lin], (self.n_train*n_perms,-1), order='F'),ctypes.c_double)
		R_d_desc_alpha_perms = shared_array(np.reshape(np.tile(r_d_desc_alpha, n_perms)[:,train_perms_lin], (self.n_train*n_perms,-1), order='F'),ctypes.c_double)

		# Precompute indices for nonzero entries in desriptor derivatices.
		n_atoms = model['z'].shape[0]
		self.d_desc_nnz_idx = np.zeros((n_atoms,n_atoms-1), dtype=np.int)
		for a in range(0,n_atoms): # for each partial deriavative
			rows,cols = np.tril_indices(n_atoms,-1)
			self.d_desc_nnz_idx[a,:] = np.concatenate([np.where( rows == a)[0], np.where( cols == a)[0]])

		# Set up multiprocessing
		if num_workers == -1:
			num_workers = mp.cpu_count() / 2

		self.n_procs = num_workers
		self.pool = mp.Pool(processes=self.n_procs)

		# Data ranges for processes
		self.wkr_starts = range(0,self.n_train,int(np.ceil(float(self.n_train)/self.n_procs)))
		self.wkr_stops = self.wkr_starts[1:] + [self.n_train]

	def __del__(self):
		self.pool.close()


	## Public ##

	def predict_new(self,R,z):

		global trainPts_D,train_perms_lin
		n_perms = train_perms_lin.shape[0] / trainPts_D.shape[1]

		# Retain backwards compatibility with old input format.
		is_old_input_format = R.ndim == 2
		if is_old_input_format:
			R = R[None,:,:]

		n_pred = R.shape[0]
		n_atoms = R.shape[1]
		dim_i = n_atoms * 3
		dim_d = (n_atoms**2 - n_atoms) / 2

		R_desc = np.empty((dim_d,n_pred))
		R_d_desc = np.empty((dim_d,dim_i,n_pred))
		for i,r in enumerate(R):
			pdist = scipy.spatial.distance.pdist(r,'euclidean')
			pdist = scipy.spatial.distance.squareform(pdist)

			# Generate input descriptor and its gradient.
			R_desc[:,i] = self._r_to_desc(r,pdist)
			R_d_desc[:,:,i] = self._r_to_d_desc(r,pdist)

		params = zip(R_desc[None,:,:].repeat(self.n_procs,0),\
					 R_d_desc[None,:,:,:].repeat(self.n_procs,0),\
					 self.wkr_starts,self.wkr_stops)

		results = self.pool.map(predict_worker_new, params)
		E = sum(E for E,_ in results) + self.c
		F = sum(F for _,F in results).reshape(-1,3,n_pred).transpose(2,0,1)

		# Retain backwards compatibility with old input format.
		if is_old_input_format:
			E = np.squeeze(E)
			F = np.squeeze(F)

		return (E,F)

	def predict(self,r,z):

		global trainPts_D,train_perms_lin
		n_perms = train_perms_lin.shape[0] / trainPts_D.shape[1]

		r = np.squeeze(r)

		pdist = scipy.spatial.distance.pdist(r,'euclidean')
		pdist = scipy.spatial.distance.squareform(pdist)

		# Generate input descriptor and its gradient.
		r_desc = self._r_to_desc(r,pdist)
		r_d_desc = self._r_to_d_desc_new(r,pdist)

		params = zip(r_desc[None,:].repeat(self.n_procs,0),\
					 r_d_desc[None,:,:].repeat(self.n_procs,0),\
					 self.wkr_starts,self.wkr_stops)

		res = sum(self.pool.map(predict_worker_cached, params))
		E = res[-1] + self.c
		F = res[:-1].reshape(-1,3)

		return (E,F)
		

	# Console output

	def print_info(self,model):

		theory_level_str = ''
		if 'theory_level' in model:
			theory_level_str = ' (' + unicode(model['theory_level']) + ')'

		print
		print "GDML Force Field Estimator %s" % ('(RF)' if model['is_rf_model'] else '')
		print '-------------------------------------------------------'
		print "Version %d, Author: Stefan Chmiela\n" % VERSION
		print "Dataset: '%s'%s" % (unicode(model['dataset']),theory_level_str)
		print "Symmetry matching: %s\n" % ('supported' if 'has_matching_support' in model else 'not supported (Check atom ordering.)')
		print 'Cross-validated prediction errors (MAE,RMSE):'
		print "|  Total Energies: %.3f/%.3f kcal/mol" % (model['e_mae'],model['e_rmse'])
		print "|  Atomic Forces:  %.3f/%.3f kcal/mol/Ang" % (model['f_mae'],model['f_rmse'])
		print '-------------------------------------------------------\n'

	def print_r(self,r,z_str):
		n_atoms = r.shape[0]

		print 'Input geometry [Ang]:'
		print '|  Atom       x            y            z'
		for a in range(0,n_atoms):
			print "|%3d: %s" % (a + 1,z_str[a]),
			for e in range(0,3):
				print "%12.4f" % r[a,e],
			print
		print

	def print_E(self,E):
		print "Total energy [kcal/mol]: %21.4f" % E

	def print_F(self,F):
		n_atoms = F.shape[0]

		print 'Atomic forces [kcal/mol/Ang]:'
		print '|  Atom       x            y            z'
		for a in range(0,n_atoms):
			print "|%3d:  " % (a + 1),
			for e in range(0,3):
				print "%12.4f" % F[a,e],
			print
		print


	## Protected ##

	# Input preprocessing

	_z_str_to_z_dict = {'O':8,'N':7,'C':6,'B':5,'H':1}
	def _z_str_to_z(self,z_str):
		return np.array([self._z_str_to_z_dict[x] for x in z_str])

	def _r_to_desc(self,r,pdist):
		n_atoms = r.shape[0]
		return 1 / pdist[np.tril_indices(n_atoms,-1)]

	def _r_to_d_desc(self,r,pdist):

		n_atoms = r.shape[0]
		d_dim = (n_atoms**2 - n_atoms)/2

		np.seterr(divide='ignore', invalid='ignore') # ignore division by zero below
		grad = np.zeros((d_dim,3*n_atoms))
		for a in range(0,n_atoms):
			d_dist = np.zeros((n_atoms,n_atoms))
			for e in range(0,3):
				d_dist[a,:] = (r[:,e] - r[a,e]) / pdist[a,:]**3
				d_dist[:,a] = d_dist[a,:]
				grad[:,3*a+e] = d_dist[np.tril_indices(n_atoms,-1)]

		return grad

	def _r_to_d_desc_new(self,r,pdist):

		n_atoms = r.shape[0]
		d_dim = (n_atoms**2 - n_atoms)/2

		np.seterr(divide='ignore', invalid='ignore') # ignore division by zero below
		grad = np.zeros((d_dim,3*n_atoms))
		for a in range(0,n_atoms):

			d_dist = (r - r[a,:]) / (pdist[a,:]**3)[:,None]

			idx = self.d_desc_nnz_idx[a,:]
			grad[idx,(3*a):(3*a+3)] = np.delete(d_dist, a, axis=0)

		return grad

	def find_optimal_permutation(self,r,z,pdist):

		w,v = np.linalg.eig(pdist)
		v_abs = abs(v[:,w.argsort()[::-1]])

		params = zip(r[None,:,:].repeat(self.n_procs,0),\
					 z[None,:].repeat(self.n_procs,0),\
					 v_abs[None,:].repeat(self.n_procs,0),\
					 pdist[None,:].repeat(self.n_procs,0),\
					 self.wkr_starts,self.proc_stops)

		results = self.pool.map(find_optimal_permutation_worker,params)

		scores = [result[0] for result in results]
		best_score_idx = scores.index(max(scores))
		best_perm = results[best_score_idx][1]

		return best_perm
