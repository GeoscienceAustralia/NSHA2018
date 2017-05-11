"""Object for storing and accessing logic tree parameters
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
    def __init__(self,set_name):
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
            branch_weight = float(row[4])
            if branch_set in self.sets:
                if branch_class in self.sets[branch_set].set_classes:
                    pass
                else:
                   self.sets[branch_set].add_class(branch_class) 
            else:
                self.add_set(branch_set)
                self.sets[branch_set].add_class(branch_class)
            self.sets[branch_set].set_classes[branch_class].add_branch(branch_value, branch_weight)

    def get_weights(self, set_name, class_name, branch_value = None):
        """Returns all branches and weights for a given set
        and class"""
        branch_values = []
        branch_weights = []
        for branch, lt_b in self.sets[set_name].set_classes[class_name].branches.iteritems():
            if branch_value is not None:
                if lt_b.value == branch_value:
                     branch_values.append(lt_b.value)
                     branch_weights.append(lt_b.weight)
            else:
                branch_values.append(lt_b.value)
                branch_weights.append(lt_b.weight)
        return branch_values, branch_weights

    def list_sets(self):
        """Print all set names
        """
        for set_name, obj in self.sets.iteritems():
            print set_name
            
    def list_classes(self, set_name=None):
        """Print class names of given set, or
        all if set_name is not specified"""
        if set_name is not None:
            print 'Set name:', set_name
            for class_name, obj in self.sets[set_name].set_classes.iteritems():
                print class_name
        else:
            for set_name, obj in self.sets.iteritems():
                 print 'Set name:', set_name
                 for class_name, obj in self.sets[set_name].set_classes.iteritems():
                     print class_name
            
if __name__ == "__main__":
    lt = LogicTree('../../shared/seismic_source_model_weights_rounded_p0.4.csv')
    lt.list_sets()
    lt.list_classes('FSM_MFD')
    branch_values, branch_weights = lt.get_weights('Mmax', 'Proterozoic')
    print branch_values
    print branch_weights
    branch_values, branch_weights = lt.get_weights('Mmax', 'Proterozoic', '7.3')
    print branch_values
    print branch_weights
    branch_values, branch_weights = lt.get_weights('FSM_MFD', 'Cratonic')
    print branch_values
    print branch_weights
