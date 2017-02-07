import scipy
import scipy.linalg as linalg
import numpy
import matplotlib.pyplot as pyplot


fig = pyplot.figure(0)
fig.clear()
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])

X = numpy.random.randn(50,2) * 10.0

X = numpy.array(X.tolist() + (numpy.random.randn(50,2)*10.0 + numpy.array([500.0, -500])).tolist() )
"""
equation is Ax + B = C

A = [ [ 1.3, 0.0002]
      [ -2.3, 1.2  ] ]


B = [ 20.3, -245.8 ]

"""

A = numpy.array([[1.3, 0.0002], [-2.3, 1.2]])
B = numpy.array([20.3, -245.8])
B = numpy.array([-10.0, 10.0])


C = A.dot(X.T).T + B

meanX = numpy.mean(X, axis=0)
meanC = numpy.mean(C, axis=0)
M = linalg.pinv(X-meanX).dot(C-meanC).T

blah = scipy.matrix([[1.0, x[0], x[1]] for x in X])

answer = linalg.pinv(blah.T.dot(blah)).dot(blah.T).dot(C)

A = answer.T.dot(blah.T).T

for x, c, a in zip(X, C, A[:,:]):
    ax.plot([x[0], c[0]], [x[1], c[1]], color = 'g')
    ax.plot([x[0], M.dot(x-meanX)[0]+meanC[0]], [x[1], M.dot(x-meanX)[1]+meanC[1]], color = 'b')
    ax.plot([x[0], a.flat[0]], [x[1], a.flat[1]], color = 'r')

fig.show()
