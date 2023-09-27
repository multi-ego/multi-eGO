import numpy as np

def create_matrix_mask(set1, set2, types, symmetrize=False, inner_op=lambda x, y: x == y, outer_op=lambda x, y: x * y):
    # symmetrize type selection
    if symmetrize: types = list(set(types + [(t[1], t[0]) for t in types]))

    mask = np.full((set1.shape[0], set2.shape[0]), False)
    for (type1, type2) in types:
        mask |= outer_op(inner_op(set1, type1), inner_op(set2, type2)[:,np.newaxis])

    return mask

def create_array_mask(set1, set2, types, symmetrize=False, inner_op=lambda x, y: x == y, outer_op=lambda x, y: x * y):
    mask = create_matrix_mask(set1, set2, types, symmetrize, inner_op, outer_op)

    return mask.flatten()

def create_linearized_mask(set1, set2, types, symmetrize=False, inner_op=lambda x, y: x == y, outer_op=lambda x, y: x * y):
    if symmetrize: types = list(set(types + [(t[1], t[0]) for t in types]))

    mask = np.full(set1.shape[0], False)
    for (type1, type2) in types:
        mask |= outer_op(inner_op(set1, type1), inner_op(set2, type2))

    return mask

def map_c12_mask(types, mask, standard_c12_dict, special_c12_dict):
    translator = lambda types, c12s_dict: np.vectorize(my_dict.__getitem__)(a)
    standard_c12 = np.where(np.logical_not(mask), translator(types, standard_c12_dict), 0.0)
    special_c12 = np.where(mask, translator(types, special_c12_dict), 0.0)

    all_c12 = standard_c12 + special_c12

    return all_c12
