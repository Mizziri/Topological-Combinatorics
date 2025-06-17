from itertools import permutations
from more_itertools import set_partitions, powerset_of_sets
from sympy.combinatorics import Permutation, PermutationGroup
from time import time
from typing import Iterable


# https://math.stackexchange.com/questions/5037772/list-of-all-non-homeomorphic-t-0-topologies-on-a-finite-set
# A kaba (per stackexchange user Kaba [kaba.hilvi.org]) is a hexidecimal representation of a topology on a finite set X
#   wherein a set U is open iff the ith bit in the decimal expansion of the kaba is 1
#       where i is the index of U in the canonical ordering of the powerset of X.

all_kabas_zero_through_five: list[list[int]] = [
    # n = 0
    [0x00000000,],
    # n = 1
    [0x00000003,],
    # n = 2
    [0x00000009, 0x0000000b, 0x0000000f,],
    # n = 3
    [0x00000081, 0x00000083, 0x00000089, 0x0000008b, 0x0000008f, 0x00000099, 0x000000ab, 0x000000af, 0x000000ff,],
    # n = 4
    [
        0x00008001, 0x00008003, 0x00008009, 0x0000800b, 0x0000800f, 0x00008081, 0x00008083, 0x00008089, 0x0000808b, 
        0x0000808f, 0x00008099, 0x000080ab, 0x000080af, 0x000080ff, 0x00008181, 0x00008241, 0x00008283, 0x000082c3, 
        0x00008383, 0x00008889, 0x0000888b, 0x0000888f, 0x00008899, 0x000088ab, 0x000088af, 0x000088bb, 0x000088ff, 
        0x00008acf, 0x00009999, 0x0000aaab, 0x0000aaaf, 0x0000aaff, 0x0000ffff,
    ],
    # n = 5
    [
        0x80000001, 0x80000003, 0x80000009, 0x8000000b, 0x8000000f, 0x80000081, 0x80000083, 0x80000089, 0x8000008b, 
        0x8000008f, 0x80000099, 0x800000ab, 0x800000af, 0x800000ff, 0x80008001, 0x80008003, 0x80008009, 0x8000800b, 
        0x8000800f, 0x80008081, 0x80008083, 0x80008089, 0x8000808b, 0x8000808f, 0x80008099, 0x800080ab, 0x800080af, 
        0x800080ff, 0x80008181, 0x80008241, 0x80008283, 0x800082c3, 0x80008383, 0x80008889, 0x8000888b, 0x8000888f, 
        0x80008899, 0x800088ab, 0x800088af, 0x800088bb, 0x800088ff, 0x80008acf, 0x80009999, 0x8000aaab, 0x8000aaaf, 
        0x8000aaff, 0x8000ffff, 0x80018001, 0x80024001, 0x80028003, 0x8002c003, 0x80038003, 0x80082003, 0x80083003, 
        0x80088009, 0x8008800b, 0x8008800f, 0x80089009, 0x8008a00b, 0x8008a00f, 0x8008b00b, 0x8008f00f, 0x80098009, 
        0x800a800b, 0x800a800f, 0x800ac00f, 0x800b800b, 0x800f800f, 0x80808081, 0x80808083, 0x80808089, 0x8080808b,
        0x8080808f, 0x80808099, 0x808080ab, 0x808080af, 0x808080ff, 0x80808181, 0x80808283, 0x808082c3, 0x80808383,
        0x80808889, 0x8080888b, 0x8080888f, 0x80808899, 0x808088ab, 0x808088af, 0x808088bb, 0x808088ff, 0x80808989, 
        0x80808a8b, 0x80808a8f, 0x80808acf, 0x80808b8b, 0x80808f8f, 0x80809999, 0x8080aaab, 0x8080aaaf, 0x8080aaff, 
        0x8080abab, 0x8080afaf, 0x8080ffff, 0x8082c0c3, 0x8088a0ab, 0x8088a0af, 0x8088a0ff, 0x8088b0bb, 0x8088f0ff, 
        0x81818181, 0x82418241, 0x82828283, 0x828282c3, 0x82828383, 0x8282c3c3, 0x83838383, 0x88888889, 0x8888888b, 
        0x8888888f, 0x88888899, 0x888888ab, 0x888888af, 0x888888bb, 0x888888ff, 0x88888acf, 0x88889999, 0x8888aaab, 
        0x8888aaaf, 0x8888aabb, 0x8888aaff, 0x8888bbbb, 0x8888ffff, 0x888acccf, 0x88aaccff, 0x99999999, 0xaaaaaaab, 
        0xaaaaaaaf, 0xaaaaaaff, 0xaaaaffff, 0xffffffff,
    ]
]


def bit_idxs(x: int) -> list[int]:
    """Return a list of indices n where the nth bit of x is set to 1 (0-indexed from the right)."""
    indices = []
    n = 0
    while x:
        if x & 1:
            indices.append(n)
        x >>= 1
        n += 1
    return indices


def int_to_bits(x: int, n:int) -> str:
    return format(x,f'0{n}b')

def bits_to_fset(bits: str) -> frozenset[int]:
    # out = set()
    # for bit, idx in zip(bits, range(len(bits))):
    #     if bit == '1':
    #         out.add(idx)
    # print(bits)
    # fset = frozenset({idx for bit, idx in zip(bits, range(len(bits)+1)) if bit == '1'})
    # print(fset)
    return frozenset({idx for bit, idx in zip(bits, range(len(bits)+1)) if bit == '1'})


def unkaba(all_kabas_on_n: list[int], n: int) -> frozenset[frozenset[int]]:
    return frozenset({bits_to_fset(int_to_bits(idx, n)) for idx in bit_idxs(all_kabas_on_n)})


def is_connected(subset: frozenset[int], topology: tuple[frozenset[int]]) -> bool:
    #check for {0,1,2}
    if len(subset) > 1:
        subspaceTopology = {subset & u for u in topology}
        for partition in set_partitions(subset, k=2):
            partition_of_sets = {frozenset(sub) for sub in partition}
            if partition_of_sets.issubset(subspaceTopology):
                return False
    return True

# Generate equivalent structures under relabeling, GPT nonsense fuckery
# Should return a static datatype, i.e. tuple or frozenset
#   Tuple cast is currently not expected to work since elements of tuple need to be static
# def canonical_form(synaptology, space: Iterable[int]) -> tuple[frozenset[int]]:
#     all_perms = PermutationGroup([Permutation(p) for p in permutations(list(space))])
#     # what fucking type is this please god help me
#     return tuple(min((((P(i) for i in subset) for subset in synaptology) for P in all_perms))) 

def canonical_form(synaptology, space: Iterable[int]) -> frozenset[frozenset[int]]:
    perms = permutations(space)
    forms = []

    for p in perms:
        mapping = dict(zip(space, p))
        transformed = frozenset(frozenset(mapping[i] for i in subset) for subset in synaptology)
        forms.append(transformed)

    return frozenset(min(forms))

def unique_synaptologies(synaptologies, label_permutations):
    seen, unique = set(), set()
    # Loop over the set of topologies or synaptologies
    for syn in synaptologies:
        if syn not in seen:
            unique.add(syn)
            equiv_structures = set()
            # Look at all possible relabellings of current structure
            for labels in label_permutations:
                new_structure = set()
                # Loop over each neighborhood/connected subspace
                for subset in syn: # Potiential issue here with unorderedness
                    new_structure.add(frozenset({labels[element] for element in subset}))
                equiv_structures.add(frozenset(new_structure))
            seen |= frozenset(equiv_structures)
    return frozenset(unique)


def connected_sets(topology: frozenset[frozenset[int]], powerset: list[set[int]]) -> frozenset[frozenset[int]]:
    return frozenset({subset for subset in powerset if is_connected(subset, topology)})


# def unique_synaptologies(n: int) -> set[tuple[frozenset[int]]]:
    
#     # tuple of ints
#     space: Iterable[int] = range(n)

#     label_permutations = list(permutations(space))

#     # This needs to be ordered to un-kaba properly, so list or tuple
#     #   Elements need to be static, so frozensets or tuples (frozensets should do)
#     powerset: frozenset[frozenset[int]] = {frozenset(s) for s in powerset_of_sets(space)}
 
#     topologies: set[frozenset[frozenset[int]]] = {unkaba(kaba) for kaba in all_kabas_zero_through_five[n]}
 
#     # print(f"Topologies:")
#     # for top in topologies:
#     #     print(top)                                
#     #print(len(topologies))
#     synaptologies: set[frozenset[frozenset[int]]] = {connected_sets(top, powerset) for top in topologies}
#     # print(f"Synaptologies:")
#     # for top in synaptologies:
#     #     print(top)                                
#     # print(len(synaptologies))

#     # for top in topologies:
#     #     if connectedSets(top,powerset) == frozenset({frozenset({2}), frozenset({0, 1, 2}), frozenset({1}), frozenset(), frozenset({0})}):
#     #         print(top)
    
#     #print(f"Synaptologies: {synaptologies}")   

#     # Uniquify by making into a set
#     #   Therefore, elements need to be static, i.e. tuples or frozensets
#     #unique_synaptologies: set[frozenset[frozenset[int]]] = {canonical_form(syn, space) for syn in synaptologies}
#     unique_synaptologies = unique_synaptologies(synaptologies, label_permutations)
#     #print(f"Unique synaptologies: {unique_synaptologies}")  
#     print(f"Unique synaptologies:")
#     for top in unique_synaptologies:
#         print(top)                                
#     print(len(unique_synaptologies))

#     return unique_synaptologies


def main() -> int:
    for n in [5]:
        start = time()
        # tuple of ints
        space: Iterable[int] = range(n)

        label_permutations = list(permutations(space))

        # This needs to be ordered to un-kaba properly, so list or tuple
        #   Elements need to be static, so frozensets or tuples (frozensets should do)
        powerset: frozenset[frozenset[int]] = {frozenset(s) for s in powerset_of_sets(space)}
    
        topologies: set[frozenset[frozenset[int]]] = {unkaba(kaba, n) for kaba in all_kabas_zero_through_five[n]}
    
        # print(f"Topologies:")
        # for top in topologies:
        #     print(top)                                
        #print(len(topologies))
        synaptologies: set[frozenset[frozenset[int]]] = {connected_sets(top, powerset) for top in topologies}
        # print(f"Synaptologies:")
        # for top in synaptologies:
        #     print(top)                                

        # for top in topologies:
        #     if connectedSets(top,powerset) == frozenset({frozenset({2}), frozenset({0, 1, 2}), frozenset({1}), frozenset(), frozenset({0})}):
        #         print(top)
        
        #print(f"Synaptologies: {synaptologies}")   

        # Uniquify by making into a set
        #   Therefore, elements need to be static, i.e. tuples or frozensets
        #unique_synaptologies: set[frozenset[frozenset[int]]] = {canonical_form(syn, space) for syn in synaptologies}
        unique = unique_synaptologies(synaptologies, label_permutations)
        #print(f"Unique synaptologies: {unique_synaptologies}")  
        print(f"Unique synaptologies:")
        for u in unique:
            pairs = []
            syn = []
            for s in u:
                if len(s) >1:
                    syn.append([el for el in s])
                if len(s) == 2:
                    pairs.append([el for el in s])
            if len(pairs) == 5:
                print("Synaptology pairs:")
                print(pairs)
                print("Synaptology:")
                print(syn)
                print("Topology:")
                for top in topologies:
                    connected = connected_sets(top, powerset)
                    syn_list = []
                    for s in connected:
                        if len(s) >1:
                            syn_list.append([el for el in s])
                    if syn_list == syn:
                        print([list(t) for t in top])        

                # IDEA: DETERMINE WHICH (N) PAIR BASES CORRESPOND TO THE SAME SYNAPTOLOGY
                # BY COMPARING TO NUMBER OF SIMPLE PLANAR GRAPHS WITH N EDGES
        print(len(unique))
        print(f"Took {time() - start} seconds.")
        
    return 0

main()



