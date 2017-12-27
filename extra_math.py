# -*- coding: utf-8 -*-
"""
EXTERNAL SOURCE
"""
from copy import deepcopy


# From here:
def mat_vec_mult(mat_a, vec_b):
    """
    Multiplying matrix with a vector, giving the result as a vector

    Source:
    https://stackoverflow.com/questions/10508021/matrix-multiplication-in-python
    """
    vec_c = [0.] * len(mat_a)
    for i, row in enumerate(mat_a):
        for j, elem in enumerate(vec_b):
            vec_c[i] += row[j] * elem
    return vec_c


def invert(mat_x):
    """
    Invert a matrix X according to gauss-jordan elimination
    In gauss-jordan elimination, we perform basic row operations to turn a matrix into
    row-echelon form.  If we concatenate an identity matrix to our input
    matrix during this process, we will turn the identity matrix into our inverse.
    X - input list of lists where each list is a matrix row
    output - inverse of X

    Source:
    http://www.vikparuchuri.com/blog/inverting-your-very-own-matrix/
    """
    # copy X to avoid altering input
    mat_x = deepcopy(mat_x)

    # Get dimensions of X
    rows = len(mat_x)
    cols = len(mat_x[0])

    # Get the identity matrix and append it to the right of mat_x
    # This is done because our row operations will make the identity into the inverse
    identity = make_identity(rows, cols)
    for i in range(0, rows):
        mat_x[i] += identity[i]

    i = 0
    for j in range(0, cols):
        # print("On col {0} and row {1}".format(j, i))
        # Check to see if there are any nonzero values below the current row in the current column
        zero_sum, first_non_zero = check_for_all_zeros(mat_x, i, j)
        # If everything is zero, increment the columns
        if zero_sum == 0:
            if j == cols:
                return mat_x
            raise Exception("Matrix is singular")
        # If mat_x[i][j] is 0, and there is a nonzero value below it, swap the two rows
        if first_non_zero != i:
            mat_x = swap_row(mat_x, i, first_non_zero)
        # Divide mat_x[i] by mat_x[i][j] to make mat_x[i][j] equal 1
        mat_x[i] = [m / mat_x[i][j] for m in mat_x[i]]

        # Rescale all other rows to make their values 0 below mat_x[i][j]
        for k in range(0, rows):
            if k != i:
                scaled_row = [mat_x[k][j] * m for m in mat_x[i]]
                mat_x[k] = [mat_x[k][m] - scaled_row[m] for m in range(0, len(scaled_row))]
        # If either of these is true, we have iterated through the matrix, and are done
        if i == rows or j == cols:
            break
        i += 1

    # Get just the right hand matrix, which is now our inverse
    for i in range(0, rows):
        mat_x[i] = mat_x[i][cols:len(mat_x[i])]
    return mat_x


def check_for_all_zeros(mat_x, i, j):
    """
    Check matrix mat_x to see if only zeros exist at or below row i in column j
    mat_x - a list of lists
    i - row index
    j - column index
    returns -
        zero_sum - the count of non zero entries
        first_non_zero - index of the first non value
    """
    non_zeros = []
    first_non_zero = -1
    for k in range(i, len(mat_x)):
        non_zero = mat_x[k][j] != 0
        non_zeros.append(non_zero)
        if first_non_zero == -1 and non_zero:
            first_non_zero = k
    zero_sum = sum(non_zeros)
    return zero_sum, first_non_zero


def swap_row(mat_x, i, j):
    """
    Swap row i and row j in a list of lists
    mat_x - list of lists
    i - row index
    j - row index
    returns- modified matrix
    """
    mat_x[j], mat_x[i] = mat_x[i], mat_x[j]
    return mat_x


def swap_col(mat_x, i, j):
    """
    Swap column i and column j in a list of lists
    mat_x - list of lists
    i - column index
    j - column index
    returns- modified matrix
    """
    for item in mat_x:
        item[i], item[j] = item[j], item[i]
    return mat_x


def make_identity(row_num, col_num):
    """
    Make an identity matrix with dimensions rxc
    row_num - number of rows
    col_num - number of columns
    returns - list of lists corresponding to  the identity matrix
    """
    identity = []
    for i in range(0, row_num):
        row = []
        for j in range(0, col_num):
            elem = 0
            if i == j:
                elem = 1
            row.append(elem)
        identity.append(row)
    return identity
