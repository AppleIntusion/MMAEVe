import numpy as np

x = np.array([[0,  1,  2],
              [3,  4,  5],
              [6,  7,  8],
              [9,  10, 11],
              [12, 13, 14]])
y = np.array([[0,  1,  2],
              [3,  4,  5],
              [6,  7,  8],
              [9,  10, 11],
              [12, 13, 14]])

def range_mean(lower, upper, data):
    return(np.mean(x[lower:upper], axis = 0))

range_mean = np.frompyfunc(range_mean, 3, 1)

lower = [0, 1, 2]
upper = [2, 3, 4]

test = range_mean(lower, upper, x)

#
#test = multiply(x, y)


