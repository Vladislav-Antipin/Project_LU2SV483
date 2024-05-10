import sys
sys.path.append('Project_library')

import embl_analysis as EMBL
import msa_analysis as MSA
from copy import deepcopy

def closest_nodes(DistMat):
    '''
    List[List[float]] -> int, int
    Returns row and column of a distance matrix corresponding to
    a  nodes with the least distance
    '''
    # min_dist : float ; minimal distance
    min_dist = DistMat[1][0]
    # i_of_min, j_of_min : int ; a row and column corresponding to
    # a pair of nodes with the least distance
    i_of_min = 1
    j_of_min = 0
    # i, j : int ; row and column of a distance matrix
    for i in range(1,len(DistMat)):
        for j in range(len(DistMat[i])-1):
            if DistMat[i][j] < min_dist :
                min_dist = DistMat[i][j]
                i_of_min = i
                j_of_min = j
    return i_of_min, j_of_min

def distance_to_parent(child_i, child_j,node_names, DistMat,tree):
    '''
    int * int * List[str] * List[List[float]] * Dict{ str : Tuple(int, Tuple(str, float), Tuple(str, float) ) }
    -> List[float]
    Assumption: the tree is organized as a dictionary {parent : (weight, (child_i_name, dist_i), (child_j_name, dist_j))},
    child_i and child_j are indices in distance matrix corresponding for child nodes of interest
    Returns a list containing distances from child_i and child_j to their parent
    '''
    # dist_i_j : List[int] ; list containing distances from child_i and child_j to their parent
    dist_i_j = []
    # parent_to_leaf : float ; distance from parent to any leaf
    parent_to_leaf = DistMat[child_i][child_j]/2
    # child_index : int ; child index that could be either child_i or child_j
    for child_index in child_i,child_j:
        # child_to_leaf : float ; a distance from a child to leaf
        child_to_leaf = 0
        # current_child_name : str
        current_child_name = node_names[child_index]
        # while current child is not a leaf <=>
        while current_child_name in tree.keys():
             # current_dist : float ; distance from a current child to it's descendant
             current_child_name, current_dist = tree[current_child_name][1]
             child_to_leaf += current_dist
        # parent_to_child : float ; distance from child_i or child_j to their parent
        parent_to_child = parent_to_leaf - child_to_leaf
        dist_i_j.append(parent_to_child)

    return dist_i_j

def add_node(tree,child_i,dist_i, child_j,dist_j,node_names, leaf_names):
    '''
    Dict{ str : Tuple(int, Tuple(str, float), Tuple(str, float) ) } * int * float * int * float * List[str] * List[str]
    -> str, int, int
    Assumption: the tree is organized as a dictionary {parent : (weight, (child_i_name, dist_i), (child_j_name, dist_j))},
    child_i and child_j are indices in distance matrix corresponding for child nodes of interest,
    dist_i and dist_j are distances from a corresponding child to their parent
    node_names is a header of the current distance matrix,
    leaf_names is a header of the original distance matrix
    Returns a generated parent name of the format node#<int> and weights (number of descendant leaves) of its children;
    Modifies a tree dictionary by adding a parent node
    '''
    # n : int ; the index of the node
    n = 1
    # parent_name : str ; a generated name of parent
    parent_name = 'node#1'
    while parent_name in tree.keys() or parent_name in leaf_names:
        n += 1
        parent_name = 'node#'+ str(n)
    # weights_i_j : List[int] ; weights (nb of leaves) of children with indices child_i and child_j
    weights_i_j = []
    # child_index : int ; child index that could be either child_i or child_j
    for child_index in child_i, child_j:
        if node_names[child_index] in tree :
            weights_i_j.append(tree[node_names[child_index]][0])
        else:
            weights_i_j.append(1)
    # weight_i, weight_j : int ; weights (nb of leaves) of children with indices child_i and child_j
    weight_i, weight_j = weights_i_j

    tree[parent_name] = (weight_i+weight_j, (node_names[child_i], dist_i), (node_names[child_j], dist_j))

    return parent_name, weight_i, weight_j


def update_nodes(child_i,child_j,parent_name,node_names):
    '''
    int * int * str * List[str] -> None
    Assumption: child_i and child_j are indices in distance matrix corresponding for child nodes of interest,
    parent_name is different from other names already present among nodes
    Modifies a header of a distance matrix (node names) by deleting the names of child_i and child_j
    and adding a parent name at the end
    '''

    # child_i_name, child_j_name : str ; names of children corresponding to indicies child_i and child_j
    child_i_name, child_j_name = node_names[child_i], node_names[child_j]
    # Hardcoding the removal of elements with 'pop' method : -------------------

    # old_node_names : List[str] ; header of distance matrix before modification
    #old_node_names = node_names[:]
    # popped : int ; counter of popped elements
    #popped = 0
    #for i in range(len(old_node_names)):
    #    if old_node_names[i] == child_i_name or old_node_names[i]==child_j_name:
    #        node_names.pop(i-popped)
    #        popped += 1

    # Another way, with 'remove' method: ---------------------------------------

    node_names.remove(child_i_name)
    node_names.remove(child_j_name)
    #---------------------------------------------------------------------------
    node_names.append(parent_name)


def update_matrix(DistMat,child_i, child_j, weight_i,weight_j):
    '''
    List[List[float]] * int * int * Dict{ str : Tuple(int, Tuple(str, float), Tuple(str, float) ) } * int * int * List[str]
    -> List[List[float]]
    Assumption: the DistMat is a distance matrix,
    the tree is organized as a dictionary {parent : (weight, (child_i_name, dist_i), (child_j_name, dist_j))},
    child_i and child_j are indices in distance matrix corresponding for child nodes of interest,
    weight_i and weight_j are weights (nb of leaves) of corresponding child nodes
    node_names is a header of the current distance matrix
    Returns the new distance matrix generated by including all distances from the old distance matrix
    ,except for those corresponding to children, and by calculating and adding the row of distances to
    a new parent node as the last row
    '''
    # NewDistMat : List[List[float]] ; new distance matrix
    NewDistMat = []
    # ParentDistRow : List[float] ;  a row in a new distance matrix corresponding to a new parent node
    ParentDistRow = []
    # i : int ; a row of distance matrix
    for i in range(len(DistMat)):
        if i != child_i and i != child_j:
            NewDistMat.append([])
            # j : int ; a column of distance matrix
            for j in range(len(DistMat[i])):
                if j != child_i and j != child_j:
                    NewDistMat[-1].append(DistMat[i][j])
            # dist_to_new_node : float ; a distance from a ith node to a new parent node
            dist_to_new_node = (MSA.triangular_to_rectangular(DistMat)[child_i][i]*weight_i + MSA.triangular_to_rectangular(DistMat)[child_j][i]*weight_j)/(weight_i+weight_j)
            ParentDistRow.append(dist_to_new_node)
    ParentDistRow.append(0.0)
    NewDistMat.append(ParentDistRow)
    return NewDistMat

def get_tree_upgma(OriginalDistMat, leaf_names):
    '''
    List[List[float]] * List[str] -> Dict{ str : Tuple(int, Tuple(str, float), Tuple(str, float) ) }, str
    Assumptions: OriginalDistMat is a distance matrix with a header leaf_names
    Returns a tree as a dictionary {parent : (weight, (child_i_name, dist_i), (child_j_name, dist_j))}
    obtained by  Unweighted Pair Group Method with Arithmetic mean (UPMGA) algorithm;
    and the name of the root node
    '''
    # CurrentDistMat : List[List[float]] ; current distance matrix
    CurrentDistMat = deepcopy(OriginalDistMat)
    # node_names : list[str] ; a header of a current distance matrix with node names
    node_names = leaf_names[:]
    # tree : Dict{ str : Tuple(int, Tuple(str, float), Tuple(str, float) ) }
    tree = {}
    while len(CurrentDistMat) > 1:

        # child_i, child_j : int ; indicies of closest nodes (currently considered as child nodes)
        child_i, child_j = closest_nodes(CurrentDistMat)

        #dist_i, dist_j : float ; distance from each of child nodes to their parent
        [dist_i, dist_j] = distance_to_parent(child_i, child_j,node_names, CurrentDistMat, tree)

        # parent_name : str ; generated name for the new parent node
        # weight_i, weight_j : int ; weights of child nodes (nb of leaves)
        parent_name, weight_i, weight_j = add_node(tree,child_i,dist_i, child_j,dist_j,node_names,leaf_names)

        update_nodes(child_i,child_j,parent_name,node_names)

        CurrentDistMat = update_matrix(CurrentDistMat,child_i, child_j, weight_i,weight_j)

    # root_name : str ; name of the root node (the last parent)
    root_name = parent_name

    return tree, root_name

def get_ordered_names_from_tree(tree, root_name):
    '''
    Dict{ str : Tuple(int, Tuple(str, float), Tuple(str, float) ) } * str -> List[str]
    Assumption: the tree is organized as a dictionary {parent : (weight, (child_i_name, dist_i), (child_j_name, dist_j))},
    root_name is the name of the root of the tree
    Returns an ordered list of leaf names (the order of appearance in a parenthesized form)
    '''
    # stack : List[str] ; a stack needed to get ordered list of leaves
    stack = [root_name]
    # leaves_in_order : List[str] ; list of ordered leaf names
    leaves_in_order = []
    while len(stack) > 0:
        # node : str ; a name of the popped node
        node = stack.pop()     # pops and returns the last element
        if node in tree.keys():
            stack.append(tree[node][2][0])
            stack.append(tree[node][1][0])
        else:
            leaves_in_order.append(node)
    return leaves_in_order

def convert_tree_from_dict_to_newik(tree, root_name):
    '''
    Dict{ str : Tuple(int, Tuple(str, float), Tuple(str, float) ) } * str -> str
    Assumption: the tree is organized as a dictionary {parent : (weight, (child_i_name, dist_i), (child_j_name, dist_j))},
    root_name is the name of the root of the tree
    Returns a tree written in a newik format
    '''
    # stack : List[str] ; a stack needed to get ordered list of leaves and parenthesis
    stack = [root_name]
    # newik : str ; a tree written in a newik format
    newik = ''
    while len(stack) > 0:
        # node : str ; a name of the popped node
        node = stack.pop()     # pops and returns the last element
        if node in tree.keys():
            stack.append(')')
            stack.append(':'+str(tree[node][2][1]))
            stack.append(tree[node][2][0])
            stack.append(',')
            stack.append(':'+str(tree[node][1][1]))
            stack.append(tree[node][1][0])
            stack.append('(')
        else:
            newik+=node
    return newik


if __name__ == '__main__':

     # Test_seqs : Dict[str:str] ; test sequences {id:sequence}
    Test_seqs = EMBL.read_fasta('Files_for_Project/SeqTest.fasta')
    # Score_dict : Dict{(str,str):float} ; {(aminoacid1, aminoacid2):score}
    Score_dict = MSA.read_score_matrix('Files_for_Project/blosum62.mat')
    
    # MSA_seqs : List[str] ; aligned sequences {id:sequence} 
    MSA_seqs = MSA.msa.star_align(Test_seqs, Score_dict)[1]
    

    # DisMat : List[List[float]] ; dissimilarity matrix
    # seq_names : List[str] ; sequences' names (header of dissimilarity matrix)
    DisMat , seq_names = MSA.compute_dissimilarity_matrix(MSA_seqs)

    # tree : Dict{ str : Tuple(int, Tuple(str, float), Tuple(str, float) ) } ; tree as a 
    #        dictionary {parent : (weight, (child_i_name, dist_i), (child_j_name, dist_j))}
    # root_name : str ; name of the root node
    tree, root_name = get_tree_upgma(DisMat, seq_names)

    print('Aligned test sequences:')
    print(*['>'+name+'\n'+ MSA_seqs[name] for name in MSA_seqs.keys()], sep='\n')
    print('Tree :',tree, sep='\n')
    print('Root :',root_name)
    print('Ordered sequences:', get_ordered_names_from_tree(tree,root_name))
    print('Newik format', convert_tree_from_dict_to_newik(tree, root_name))