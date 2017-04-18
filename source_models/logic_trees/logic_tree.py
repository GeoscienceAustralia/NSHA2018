"""Object for storing logic tree parameters
"""

import os, sys
import numpy as np

class LogicTreeBranch(object):
    """Class for storing logic tree parameters and weights
    for an individual branch
    """
    
    def __init__(self, branch_value, branch_weight):
        """Instantiate class
        """
#        self.branch_class = branch_class
        self.value = branch_value
        self.weight = branch_weight

class LogicTreeClass(object):
    """Class for storing multiple logic tree branches 
    related to similar parameters
    """
    def __init__(self,class_name):
        self.class_name = class_name
        self.branches = {}

    def add_branch(self, branch_value, weight):
        self.branches[branch_value] = LogicTreeBranch(branch_value, weight)

        
class LogicTreeSet(object):
    """Class for storing multiple logic tree branches 
    related to similar parameters
    """
    def __init__(self,set_name):#, set_classes, set_values, set_weights):
        self.set_name = set_name
        self.set_classes = {}
    
    def add_class(self, class_name):
        self.set_classes[class_name] = LogicTreeClass(class_name)
        
class LogicTree(object):
    """Class for storing logic tree parameters and weights
    for all branches and accessing these
    """
    def __init__(self, filename=None):
        self.sets = {}
        if filename is not None:
            self.lt_from_file(filename)
        #else:
        #    return

    def add_set(self, set_name):
        self.sets[set_name] = LogicTreeSet(set_name)

    def lt_from_file(self,filename):
        """Build logic trees from parameters in a file
        """
        
        f_in = open(filename,'r')
        header = f_in.readline()
        for line in f_in.readlines():
            row = line.rstrip().split(',')
            branch_set = row[1]
            branch_class = row[2]
            branch_value = row[3]
            branch_weight = row[4]
            if branch_set in self.sets:
                if branch_class in self.sets[branch_set].set_classes:
                    pass
                else:
                   self.sets[branch_set].add_class(branch_class) 
            else:
                self.add_set(branch_set)
                self.sets[branch_set].add_class(branch_class)
            # Add individual branch values and weightings
#            print self.sets[branch_set]
            self.sets[branch_set].set_classes[branch_class].add_branch(branch_value, branch_weight)
        print self.sets
            
if __name__ == "__main__":
    lt = LogicTree('../../shared/seismic_source_model_weights_rounded_p0.4.csv')
    print lt.sets['Mmax'].set_classes['Proterozoic'].branches
    for branch, lt_b in lt.sets['Mmax'].set_classes['Proterozoic'].branches.iteritems():
        print lt_b.value, lt_b.weight
