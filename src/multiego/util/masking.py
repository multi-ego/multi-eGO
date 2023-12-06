import numpy as np

def create_matrix_mask(set1, set2, types, symmetrize=False, inner_op=lambda x, y: x == y, outer_op=lambda x, y: x * y):
    """
    Creates a boolean matrix mask based on comparison operations between two sets.

    Args:
    - set1 (numpy.ndarray): First set of values.
    - set2 (numpy.ndarray): Second set of values.
    - types (list): List of tuples representing types for comparison.
    - symmetrize (bool): Flag to determine whether to symmetrize type selection.
    - inner_op (function): Operation for element-wise comparison within each set.
    - outer_op (function): Operation for the combination of comparisons between sets.

    Returns:
    - numpy.ndarray: A boolean matrix mask reflecting the comparison operations.

    Note:
    - `inner_op` should be a function that compares an element of `set1` with a type.
    - `outer_op` should be a function that combines comparisons between `set1` and `set2`.
    """
    # symmetrize type selection
    if symmetrize:
        types = list(set(types + [(t[1], t[0]) for t in types]))

    mask = np.full((set1.shape[0], set2.shape[0]), False)
    for (type1, type2) in types:
        mask |= outer_op(inner_op(set1, type1), inner_op(set2, type2)[:, np.newaxis]).T

    return mask


def create_array_mask(set1, set2, types, symmetrize=False, inner_op=lambda x, y: x == y, outer_op=lambda x, y: x * y):
    """
    Creates a linearized boolean array mask from a matrix mask based on comparison operations between two sets.

    Args:
    - set1 (numpy.ndarray): First set of values.
    - set2 (numpy.ndarray): Second set of values.
    - types (list): List of tuples representing types for comparison.
    - symmetrize (bool): Flag to determine whether to symmetrize type selection.
    - inner_op (function): Operation for element-wise comparison within each set.
    - outer_op (function): Operation for the combination of comparisons between sets.

    Returns:
    - numpy.ndarray: A linearized boolean array mask reflecting the comparison operations.

    Note:
    - This function is derived from create_matrix_mask() to provide a flattened mask array.
    """
    mask = create_matrix_mask(set1, set2, types, symmetrize, inner_op, outer_op)

    return mask.flatten()


def create_linearized_mask(set1, set2, types, symmetrize=False, inner_op=lambda x, y: x == y, outer_op=lambda x, y: x * y):
    """
    Creates a linearized boolean mask based on comparison operations between two sets.

    Args:
    - set1 (numpy.ndarray): First set of values.
    - set2 (numpy.ndarray): Second set of values.
    - types (list): List of tuples representing types for comparison.
    - symmetrize (bool): Flag to determine whether to symmetrize type selection.
    - inner_op (function): Operation for element-wise comparison within each set.
    - outer_op (function): Operation for the combination of comparisons between sets.

    Returns:
    - numpy.ndarray: A linearized boolean mask reflecting the comparison operations.

    Note:
    - This function does not provide a flattened mask array but operates on a 1D array.
    """
    if symmetrize:
        types = list(set(types + [(t[1], t[0]) for t in types]))

    mask = np.full(set1.shape[0], False)
    for (type1, type2) in types:
        mask |= outer_op(inner_op(set1, type1), inner_op(set2, type2))

    return mask


def map_c12_mask(types, mask, standard_c12_dict, special_c12_dict):
    """
    Maps the c12 values based on a mask to derive a combined c12 array.

    Args:
    - types (numpy.ndarray): Array of types used for the mask.
    - mask (numpy.ndarray): Boolean mask array.
    - standard_c12_dict (dict): Dictionary mapping types to standard c12 values.
    - special_c12_dict (dict): Dictionary mapping types to special c12 values.

    Returns:
    - numpy.ndarray: A combined array of c12 values based on the mask and dictionaries.
    """
    translator = lambda types, c12s_dict: np.vectorize(c12s_dict.__getitem__)(types)
    standard_c12 = np.where(np.logical_not(mask), translator(types, standard_c12_dict), 0.0)
    special_c12 = np.where(mask, translator(types, special_c12_dict), 0.0)

    all_c12 = standard_c12 + special_c12

    return all_c12
