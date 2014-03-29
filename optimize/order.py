'''
solves the order ILP
'''
'''
Created on Mar 21, 2011

solves matching using CPLEX

@author: jlindsay
'''

# system imports.
import sys
import logging
import numpy as np
import cplex
import math
import networkx as nx

from cplex.exceptions import CplexSolverError

class OrderIlp(object):
	'''
	solves bi-partite matching using ILP
	'''
	
	def __init__(self, log_file, err_file):
		'''
		constructor
		'''	
		
		# save file ptrs.
		self._log_file = log_file
		self._err_file = err_file
		
		# clear logs.
		tmp = open(self._log_file, "w")
		tmp.close()
		tmp = open(self._err_file, "w")
		tmp.close()
		
		# set loaded var.
		self._loaded = False
		self._solved = False
		
	def load(self, matching_type, DG, card_val=False):
		''' loads variables from flow graph '''
		
		# sanity check.
		if self._loaded == True:
			logging.error("ILP already loaded.")
			sys.exit(1)

		# save reference ot graph.
		self._graph = DG
		
		# initiate cplex object.
		self._cpx = cplex.Cplex()
		
		# set log files.
		self._cpx.set_log_stream(self._log_file)
		self._cpx.set_results_stream(self._log_file)
		self._cpx.set_warning_stream(self._err_file)
		self._cpx.set_error_stream(self._err_file)
		
		# prepare lookup structures.
		self._var_defined = set()
		
		# add Xij variables.
		self._add_pair_vars()
		
		# constrain paths.
		self._constrain_paths()
		
		# build the objective.
		if matching_type == "weight":
			self._obj_weight()
		elif matching_type == "card":
			self._obj_card()
		elif matching_type == "mcmw":
			self._constrain_card(card_val)
			self._obj_weight()
		
		# set loaded.
		self._loaded = True			

	def solve(self, file_path=False):
		''' runs the ilp on loaded info '''

		# sanity check.
		if self._loaded == False:
			logging.error("ILP not loaded.")
			sys.exit(1)			
			
		# sanity check.
		if self._solved == True:
			logging.error("shouldn't solve ILP twice.")
			sys.exit(1)			
			
		# write ILP to file.
		if file_path != False:
			self._cpx.write(file_path, filetype="lp")
			
		# call the solve code.
		try:
			
			# call the solve method.	
			self._cpx.solve()
			
			# populate solution.
			self._DG = self._populate_sol()
			
		except CplexSolverError, e:
			
			# if no solution found return empty sol and -1.
			self._sol = None
			
			logging.error("exception raised during solve: " + str(e))
			sys.exit(1)
	
		# set solved to true.
		self._solved = True
		
		# return solution.
		#return self._sol, self._cpx.solution.get_objective_value()
		return self._DG

	def clear(self):
		''' resets ILP completely '''
		
		# sanity.
		if self._cpx == None:
			logging.error("ILP already deleted")
			sys.exit(1)
			
		# sanity.
		if self._solved == False:
			logging.error("ILP not solved")
			sys.exit(1)
		
		# remove cplex and other vars.
		del self._cpx
		del self._var_defined
		self._cpx = None	
	
		# clear loaded.
		self._loaded = False
		self._solved = False

	def _obj_weight(self):
		''' sets objective '''
		
		# loop over bundles.
		for e in self._graph.edges():
			
			# simplify.
			idxa = e[0]
			idxb = e[1]
			
			# build vars.
			Xij = "X#%s#%s" % (str(idxa), str(idxb))
			
			# get weight.
			Wij = self._graph[idxa][idxb]['bcnts']
            
			# add to objective.
			self._cpx.objective.set_linear(Xij, Wij)
				
		# set objective type.
		self._cpx.objective.set_sense(self._cpx.objective.sense.maximize)

	def _obj_card(self):
		''' sets objective '''
		
		# loop over bundles.
		for e in self._graph.edges():
			
			# simplify.
			idxa = e[0]
			idxb = e[1]
			
			# build vars.
			Xij = "X#%s#%s" % (str(idxa), str(idxb))
			
			# set simple weight.
			self._cpx.objective.set_linear(Xij, 1)
			
				
		# set objective type.
		self._cpx.objective.set_sense(self._cpx.objective.sense.maximize)
		
	def _constrain_paths(self):
		''' ensures each variable has in/out degree at most 1'''
		
		# loop over each node.
		for p in self._graph.nodes():
			
			# constrain in degree.
			inds = list()
			for q in self._graph.predecessors(p):	
				inds.append("X#%s#%s" % (str(q), str(p)))
			vals = [1] * len(inds)


			# build constraint.
			c1 = cplex.SparsePair( ind = inds, val = vals )
			
			# constrain out degree.
			inds = list()
			for q in self._graph.successors(p):	
				inds.append("X#%s#%s" % (str(p), str(q)))
			vals = [1] * len(inds)
			
			# build constraint.
			c2 = cplex.SparsePair( ind = inds, val = vals )
						
			# add them.
			self._cpx.linear_constraints.add( \
				lin_expr = [c1, c2],\
				senses = ["L", "L"],\
				rhs = [1, 1],\
				names = ['pair', 'pair']\
			)

	def _constrain_card(self, card_val):
		''' ensures paths have atleast a certain cardinality'''
		
		# loop over bundles.
		inds = list()
		for e in self._graph.edges():
			
			# simplify.
			idxa = e[0]
			idxb = e[1]
			
			# build vars.
			inds.append("X#%s#%s" % (str(idxa), str(idxb)))

		# build constraint.
		c = cplex.SparsePair( ind = inds, val = [1] * len(inds) )
									
		# add it.
		self._cpx.linear_constraints.add( \
			lin_expr = [c],\
			senses = ["E"],\
			rhs = [card_val],\
			names = ['card_val']\
		)

		
	def _add_pair_vars(self):
		''' adds pair variables '''
		
		# loop over bundles.
		for e in self._graph.edges():
			
			# simplify.
			idxa = e[0]
			idxb = e[1]
			
			# build vars.
			Xij = "X#%s#%s" % (str(idxa), str(idxb))
			
			# add the variables.
			self._cpx.variables.add( lb = [0], ub = [1], types = ["B"], names = [Xij] )	
			
			# populate lookup.
			self._var_defined.add(Xij)
	
	def _populate_sol(self):
		''' populates solution object after running '''
		
		# loop over bundles.
		elist = list()
		for e in self._graph.edges():
			
			# simplify.
			idxa = e[0]
			idxb = e[1]
			
			# build vars.
			Xij = "X#%s#%s" % (str(idxa), str(idxb))
			
			# get result.
			val = int(self._cpx.solution.get_values(Xij))
			
			# add to list if chosen.
			if val == 1:
				elist.append((idxa, idxb))
				
		# create a new graph.
		DG = nx.DiGraph()
		for n in self._graph.nodes():
			DG.add_node(n, self._graph.node[n])
			
		for p, q in elist:
			DG.add_edge(p, q, self._graph[p][q])
				
		# return the new graph.
		return DG
		#self._sol = elist
