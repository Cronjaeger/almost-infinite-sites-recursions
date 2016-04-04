import numpy as np

def csv2psi(fileName):
    """
    Takes the path of a .csv-file as input. The last column is presumed to
    encode the "n"-vector enumerating the occurrance of different alleles; the first collumn is presumed to
    enumerate the multiplicity of each column. The other columns are presuemd to encode haplotypes.

    If the below were the contents of mycsv:
    "
    5, 1, 1, 1, 1, 1, 1, 1, 1, 2, 0
    0, 1, 1, 0, 0, 2, 0, 0, 3, 3, 22
    0, 0, 2, 1, 0, 3, 0, 0, 2, 1, 18
    1, 0, 1, 0, 0, 2, 0, 0, 0, 0, 12
    0, 1, 1, 0, 2, 2, 0, 0, 3, 3, 6
    0, 1, 1, 0, 0, 2, 0, 0, 0, 3, 6
    1, 2, 1, 0, 0, 0, 2, 0, 0, 0, 4
    1, 1, 1, 3, 0, 3, 0, 0, 3, 3, 3
    1, 1, 1, 0, 0, 2, 0, 0, 3, 3, 3
    0, 1, 1, 2, 0, 2, 0, 0, 3, 3, 3
    1, 0, 1, 1, 2, 2, 0, 0, 0, 0, 3
    0, 1, 1, 0, 2, 2, 0, 1, 3, 3, 3
    0, 3, 1, 0, 0, 2, 0, 0, 3, 3, 3
    3, 1, 1, 0, 0, 2, 0, 2, 3, 3, 3
    0, 3, 1, 1, 0, 3, 0, 0, 0, 1, 2
    1, 0, 1, 0, 0, 2, 1, 0, 0, 0, 2
    0, 0, 1, 1, 2, 3, 0, 0, 2, 1, 2
    0, 0, 1, 1, 0, 3, 0, 0, 0, 1, 1
    0, 0, 2, 1, 0, 3, 0, 0, 0, 1, 1
    1, 0, 1, 1, 2, 0, 0, 0, 0, 0, 1
    1, 1, 1, 3, 0, 3, 0, 3, 3, 3, 1
    3, 2, 2, 1, 0, 3, 0, 0, 0, 1, 1
    "
    Then psiFromCSV("mycsv") should return S,nr,nc whereby:

    >>> S
    array([[0, 1, 1, 0, 0, 2, 0, 0, 3, 3],
       [0, 0, 2, 1, 0, 3, 0, 0, 2, 1],
       [1, 0, 1, 0, 0, 2, 0, 0, 0, 0],
       [0, 1, 1, 0, 2, 2, 0, 0, 3, 3],
       [0, 1, 1, 0, 0, 2, 0, 0, 0, 3],
       [1, 2, 1, 0, 0, 0, 2, 0, 0, 0],
       [1, 1, 1, 3, 0, 3, 0, 0, 3, 3],
       [1, 1, 1, 0, 0, 2, 0, 0, 3, 3],
       [0, 1, 1, 2, 0, 2, 0, 0, 3, 3],
       [1, 0, 1, 1, 2, 2, 0, 0, 0, 0],
       [0, 1, 1, 0, 2, 2, 0, 1, 3, 3],
       [0, 3, 1, 0, 0, 2, 0, 0, 3, 3],
       [3, 1, 1, 0, 0, 2, 0, 2, 3, 3],
       [0, 3, 1, 1, 0, 3, 0, 0, 0, 1],
       [1, 0, 1, 0, 0, 2, 1, 0, 0, 0],
       [0, 0, 1, 1, 2, 3, 0, 0, 2, 1],
       [0, 0, 1, 1, 0, 3, 0, 0, 0, 1],
       [0, 0, 2, 1, 0, 3, 0, 0, 0, 1],
       [1, 0, 1, 1, 2, 0, 0, 0, 0, 0],
       [1, 1, 1, 3, 0, 3, 0, 3, 3, 3],
       [3, 2, 2, 1, 0, 3, 0, 0, 0, 1]])
    >>> nr
    array([22, 18, 12,  6,  6,  4,  3,  3,  3,  3,  3,  3,  3,  2,  2,  2,  1,
        1,  1,  1,  1])
    >>> nc
    array([5, 1, 1, 1, 1, 1, 1, 1, 1, 2])
    """
    raw = np.genfromtxt(fileName, delimiter=',')
    # print raw

    nr = np.array(raw[1:,-1],  dtype = int)
    nc = np.array(raw[0,:-1],  dtype = int)
    S  = np.array(raw[1:,:-1], dtype = int)

    S,nr = removeDuplicateRows(S,nr)
    St, nc = removeDuplicateRows(np.transpose(S),nc)
    S = np.transpose(St)

    return S,nr,nc


def removeDuplicateRows(S,n):

    rows,columns = S.shape
    row_counts = {}

    rows = S.shape[0]

    for row_index in range(rows):
        row = S[row_index]
        rowAsTuple = tuple(row)
        if rowAsTuple in row_counts:
            row_counts[rowAsTuple] += n[row_index]
        else:
            row_counts[rowAsTuple] = n[row_index]

    rowList = row_counts.keys()
    rowCounts = row_counts.values()

    S_new = np.array(np.r_["0,2",rowList],dtype = int)
    n_vec = np.array(rowCounts,dtype = int)

    return S_new,n_vec
