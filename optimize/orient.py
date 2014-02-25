'''
solves the orientation ILP
'''

# system imports.
import sys
import logging
import numpy as np
import math

import cplex
from cplex.exceptions import CplexSolverError


class OrientIlp(object):
    '''
    implements SPQR ILP using CPLEX
    '''

    def __init__(self, log_file, err_file, prg_file, sol_file, weight_mode):
        '''
        constructor
        '''

        # save the weighting mode.
        self.weight_mode = weight_mode
        

        # save file ptrs.
        self._log_file = log_file
        self._err_file = err_file
        self._prg_file = prg_file
        self._sol_file = sol_file

        # clear logs.
        tmp = open(self._log_file, "w")
        tmp.close()
        tmp = open(self._err_file, "w")
        tmp.close()
        tmp = open(self._prg_file, "w")
        tmp.close()
        tmp = open(self._sol_file, "w")
        tmp.close()

        # set loaded var.
        self._loaded = False
        self._solved = False

    def solve(self):
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
        self._cpx.write(self._prg_file, filetype="lp")

        # call the solve code.
        try:

            # call the solve method.
            self._cpx.set_problem_type(self._cpx.problem_type.MILP)
            self._cpx.solve()

            # populate solution.
            sol = self._populate_sol()

        except CplexSolverError, e:

            logging.error("exception raised during solve: " + str(e))
            sys.exit(1)

        # set solved to true.
        self._solved = True

        if sol['obj'] == None:
            logging.error("got none for soluti9on???")
            sys.exit()

        # return the dictionary.
        return sol

    def load(self, G):
        ''' loads the node, bundle and triangle lists'''

        # sanity check.
        if self._loaded == True:
            logging.error("ILP already loaded.")
            sys.exit(1)

        # save pointer to graph.
        self._G = G

        # initiate cplex object.
        self._cpx = cplex.Cplex()

        # set log files.
        self._cpx.set_log_stream(self._log_file)
        self._cpx.set_results_stream(self._log_file)
        self._cpx.set_warning_stream(self._err_file)
        self._cpx.set_error_stream(self._err_file)

        # limit resources.
        self._cpx.parameters.threads.set(10)
        self._cpx.parameters.mip.polishing.time.set(1800)

        # set loaded.
        self._loaded = True

        # add variables.
        self._add_variables()

        # add constraints.
        self._orientation_pairing()
        self._state_coupling()
        self._forbid_2()
        self._forbid_3()

        # add objective.
        self._add_objective()

    def fix(self, fix):
        ''' called after load to fix orientation variables '''
        
        # sanity.
        if self._loaded == False:
            logging.error("problem with fixing logic")
            sys.exit()
            
        # sanity.
        if fix == None: return
            
        # fix values.
        keys = fix[0]
        vals = fix[1]
        for n, v in zip(keys, vals):
            Si = "S#%s" % str(n)
            c = cplex.SparsePair(ind=[Si], val=[1])
            self._cpx.linear_constraints.add(lin_expr = [c], senses = ['E'], rhs = [v], names = ['fix'])
            
    def objmod1(self, objlist):
        ''' called after load to modify objective '''
        
        # sanity.
        if self._loaded == False:
            logging.error("problem with fixing logic")
            sys.exit()
            
        # sanity.
        if objlist == []: return
            
        # fix values.
        for n, v0, v1 in objlist:
            
            # create dummy var constrain to be 1.
            Di = "D#%s" % str(n)
            self._cpx.variables.add(lb=[1],ub=[1],types = ["B"],names = [Di])
            
            # add to objective.
            Si = "S#%s" % str(n)
            self._cpx.objective.set_linear(Di, v0)
            self._cpx.objective.set_linear(Si, v1-v0)

    def objmod2(self, objlist):
        ''' called after load to modify objective '''
        
        # sanity.
        if self._loaded == False:
            logging.error("problem with fixing logic")
            sys.exit()
            
        # sanity.
        if objlist == []: return
            
        # fix values.
        for n0, n1, v00, v01, v10, v11 in objlist:

            # create dummyvar
            Di = "D#%s#%s" % (str(n0), str(n1))
            self._cpx.variables.add(lb=[1],ub=[1],types = ["B"],names = [Di])
                        
            # add to objective.
            if (n0, n1) in self._Oij:
                Oij = self._Oij[(n0, n1)]
            elif (n1, n0) in self._Oij:
                Oij = self._Oij[(n1, n0)]
            else:
                # create new dummy var with no weights
                Oij = "O#%s#%s" % (str(n0), str(n1))
                self._Oij[(n0,n1)] = Oij
                self._cpx.variables.add(lb=[0],ub=[1],types = ["B"],names = [Oij])
                
            self._cpx.objective.set_linear(Di, v00)
            self._cpx.objective.set_linear(Oij, v01)
            self._cpx.objective.set_linear(Oij, -1*v00)

    def _populate_sol(self):
        ''' populates solution object after running '''

        # create dictionary.
        sol = {'orien':dict(), 'state':dict(), 'obj':None}

        # loop over nodes.
        for n in self._G.nodes():
            sol['orien'][n] = self._cpx.solution.get_values(self._Si[n]) >= 0.5

        # loop over edges.
        for p, q in self._G.edges():

            # get result.
            valA = self._cpx.solution.get_values(self._Aij[(p,q)]) >= 0.5
            valB = self._cpx.solution.get_values(self._Bij[(p,q)]) >= 0.5
            valC = self._cpx.solution.get_values(self._Cij[(p,q)]) >= 0.5
            valD = self._cpx.solution.get_values(self._Dij[(p,q)]) >= 0.5

            # sanity check.
            if sum([valA, valB, valC, valD]) > 1:
                logging.error("more states chosen")
                sys.exit(1)

            # save state.
            if valA == 1:
                sol['state'][tuple(sorted([p,q]))] = 0
            elif valB == 1:
                sol['state'][tuple(sorted([p,q]))] = 1
            elif valC == 1:
                sol['state'][tuple(sorted([p,q]))] = 2
            elif valD == 1:
                sol['state'][tuple(sorted([p,q]))] = 3
            else:
                logging.error("no state chosen")
                sys.exit(1)

        # save objective.
        sol['obj'] = self._cpx.solution.get_objective_value()
        return sol

    def _forbid_3(self):
        ''' no three cycles '''

        # build neighbor sets.
        nsets = dict()
        for n in self._G.nodes():
            nsets[n] = set(self._G.neighbors(n))

        # loop over each edge.
        triangles = set()
        for p, q in self._G.edges():

            # look for intersection in their neighbors more than 1.
            isec = nsets[p].intersection(nsets[q])

            # found a triangle.
            if len(isec) > 1:

                # add triangles to set.
                for z in isec:
                    triangles.add(tuple(sorted([p, q, z])))

        # get existing variables.
        vrs =  set(self._cpx.variables.get_names())

        # make contraints.
        for i, j, k in triangles:

            # make vars.
            Aij = "A#%s#%s" % (str(i), str(j))
            Ajk = "A#%s#%s" % (str(j), str(k))
            Aki = "A#%s#%s" % (str(k), str(i))

            Bij = "B#%s#%s" % (str(i), str(j))
            Bjk = "B#%s#%s" % (str(j), str(k))
            Bki = "B#%s#%s" % (str(k), str(i))

            Cij = "C#%s#%s" % (str(i), str(j))
            Cjk = "C#%s#%s" % (str(j), str(k))
            Cki = "C#%s#%s" % (str(k), str(i))

            Dij = "D#%s#%s" % (str(i), str(j))
            Djk = "D#%s#%s" % (str(j), str(k))
            Dki = "D#%s#%s" % (str(k), str(i))

            # make list.
            vlist = list()
            vlist.append((Aij, Ajk, Aki))#[(Aij, Ajk, Aki, Bij, Bjk, Bki, Cij, Cjk, Cki, Dij, Djk, Dki]

            # check each constraint and add.
            for i in range(len(vlist)):

                # skip if variable not defined.
                x, y, z = vlist[i]
                if x not in vrs: continue
                if y not in vrs: continue
                if z not in vrs: continue

                # create constraint.
                c = cplex.SparsePair(ind=[x,y,z], val=[ 1, 1, 1 ])

                # add to cplex.
                self._cpx.linear_constraints.add(\
                    lin_expr = [c],\
                    senses = ['L'],\
                    rhs = [2],\
                    names = ['4%i' % i]\
                )

    def _forbid_2(self):
        ''' no two cycles '''

        # loop over each pair.
        for p, q in self._G.edges():

            # make variables.
            Oij = self._Oij[(p,q)]
            Aij = self._Aij[(p,q)]
            Bij = self._Bij[(p,q)]
            Cij = self._Cij[(p,q)]
            Dij = self._Dij[(p,q)]

            # simplify.
            wtA = self._G[p][q]['bcnts'][0]
            wtB = self._G[p][q]['bcnts'][1]
            wtC = self._G[p][q]['bcnts'][2]
            wtD = self._G[p][q]['bcnts'][3]

            # make constraints.
            lin_expr = list()
            senses = list()
            rhs = list()
            names = list()

            # check both.
            if wtA != 0 and wtD != 0:
                lin_expr.append(cplex.SparsePair(ind=[Aij,Dij,Oij], val=[ 1, 1, 1 ]))
                senses.append("L")
                rhs.append(1)
                names.append("3a")
            if wtB != 0 and wtC != 0:
                lin_expr.append(cplex.SparsePair(ind=[Bij,Cij,Oij], val=[ 1, 1,-1 ]))
                senses.append("L")
                rhs.append(0)
                names.append("3b")

            # add to cplex.
            self._cpx.linear_constraints.add(\
                lin_expr = lin_expr,\
                senses = senses,\
                rhs = rhs,\
                names = names\
            )

    def _state_coupling(self):
        ''' state coupling '''

        # loop over each pair.
        for p, q in self._G.edges():

            # get variables.
            Si = self._Si[p]
            Sj = self._Si[q]
            Aij = self._Aij[(p,q)]
            Bij = self._Bij[(p,q)]
            Cij = self._Cij[(p,q)]
            Dij = self._Dij[(p,q)]

            # make constraints.
            c2a = cplex.SparsePair(ind=[ Aij, Si, Sj ], val=[ 2, 1, 1 ])
            c2b = cplex.SparsePair(ind=[ Bij, Si, Sj ], val=[ 2, 1,-1 ])
            c2c = cplex.SparsePair(ind=[ Cij, Si, Sj ], val=[ 2,-1, 1 ])
            c2d = cplex.SparsePair(ind=[ Dij, Si, Sj ], val=[ 2,-1,-1 ])

            self._cpx.linear_constraints.add(lin_expr = [c2a], senses = ['L'], rhs = [2], names = ['2a'])
            self._cpx.linear_constraints.add(lin_expr = [c2b], senses = ['L'], rhs = [1], names = ['2b'])
            self._cpx.linear_constraints.add(lin_expr = [c2c], senses = ['L'], rhs = [1], names = ['2c'])
            self._cpx.linear_constraints.add(lin_expr = [c2d], senses = ['L'], rhs = [0], names = ['2d'])

            # add sum constraint.
            csum = cplex.SparsePair(ind=[ Aij, Bij, Cij, Dij ], val=[ 1, 1, 1, 1 ])
            self._cpx.linear_constraints.add(lin_expr = [csum], senses = ['E'], rhs = [1], names = ['2e'])

    def _orientation_pairing(self):
        ''' orientation pairing '''

        # loop over each pair.
        for p, q in self._G.edges():

            # make variables.
            Si = "S#%s" % str(p)
            Sj = "S#%s" % str(q)
            Oij = "O#%s#%s" % (str(p), str(q))

            # make constraints.
            c1 = cplex.SparsePair(ind=[Oij,Si,Sj], val=[ 1,-1,-1])
            c2 = cplex.SparsePair(ind=[Oij,Si,Sj], val=[ 1, 1, 1])
            c3 = cplex.SparsePair(ind=[Oij,Si,Sj], val=[ 1, 1,-1])
            c4 = cplex.SparsePair(ind=[Oij,Si,Sj], val=[ 1,-1, 1])
            names = ["1a", "1b", "1c", "1d"]

            # add to cplex.
            self._cpx.linear_constraints.add(\
                lin_expr = [c1, c2, c3, c4],\
                senses = ["L", "L", "G", "G"],\
                rhs = [0, 2, 0, 0],\
                names = names\
            )

    def _add_objective(self):
        ''' maximize weighted edges '''

        # add to objective.
        for p, q in self._G.edges():

            # convert counts to weights.
            counts = self._G[p][q]['bcnts']
            weights = [x * self._G[p][q]['cov'] for x in counts]
            uniqs = self._G[p][q]['u']
            combed = [weights[i]*uniqs[i] for i in range(4)]
            
            # get constants.
            if self.weight_mode == 0:
                wtA, wtB, wtC, wtD = counts
            elif self.weight_mode == 1:
                wtA, wtB, wtC, wtD = weights
            elif self.weight_mode == 2:
                wtA, wtB, wtC, wtD = uniqs
            elif self.weight_mode == 3:
                wtA, wtB, wtC, wtD = combed
            else:
                logging.error("unknown weight")
                sys.exit(1)

            # get variables.
            Aij = self._Aij[(p,q)]
            Bij = self._Bij[(p,q)]
            Cij = self._Cij[(p,q)]
            Dij = self._Dij[(p,q)]

            self._cpx.objective.set_linear(Aij, wtA)
            self._cpx.objective.set_linear(Bij, wtB)
            self._cpx.objective.set_linear(Cij, wtC)
            self._cpx.objective.set_linear(Dij, wtD)


        # set objective type.
        self._cpx.objective.set_sense(self._cpx.objective.sense.maximize)

    def _add_variables(self):
        ''' adds all variables '''

        # create orientation vars.
        si_list = list()
        self._Si = dict()
        for n in self._G.nodes():
            Si = "S#%s" % str(n)
            si_list.append(Si)
            self._Si[n] = Si

        # add orientation vars.
        self._cpx.variables.add(\
            lb = [0] * len(si_list),\
            ub = [1] * len(si_list),\
            types = ["B"] * len(si_list),\
            names = si_list\
        )

        # create consistency vars.
        oij_list = list()
        self._Oij = dict()
        for p, q in self._G.edges():
            Oij = "O#%s#%s" % (str(p), str(q))
            oij_list.append(Oij)
            self._Oij[(p,q)] = Oij

        # add consistency vars.
        self._cpx.variables.add(\
            lb = [0] * len(oij_list),\
            ub = [1] * len(oij_list),\
            types = ["B"] * len(oij_list),\
            names = oij_list\
        )

        # create state vars.
        states_list = list()
        self._Aij = dict()
        self._Bij = dict()
        self._Cij = dict()
        self._Dij = dict()
        for p,q in self._G.edges():
            Aij = "A#%s#%s" % (str(p), str(q))
            Bij = "B#%s#%s" % (str(p), str(q))
            Cij = "C#%s#%s" % (str(p), str(q))
            Dij = "D#%s#%s" % (str(p), str(q))
            states_list.append(Aij)
            states_list.append(Bij)
            states_list.append(Cij)
            states_list.append(Dij)
            self._Aij[(p,q)] = Aij
            self._Bij[(p,q)] = Bij
            self._Cij[(p,q)] = Cij
            self._Dij[(p,q)] = Dij

        # add consistency vars.
        self._cpx.variables.add(\
            lb = [0] * len(states_list),\
            ub = [1] * len(states_list),\
            types = ["B"] * len(states_list),\
            names = states_list\
        )



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
        self._cpx = None

        # zero out dummy vars.
        self._dummy_idx = 0

        # clear loaded.
        self._loaded = False
        self._solved = False

