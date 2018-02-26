import numpyClone as np

print(np.arange(1,2,0.1))
print(np.linspace(1,2,9))
print(np.pi)
print(np.zeros(5))
print(np.zeros((5,3)))
print(np.ones((3,5)))
print(np.concatenate((range(10), range(2,12))))

x = np.arange(9,0,-1)
print(x)
print(x[2:7:2])
print(x[[3, 3, 1, 8]])
x[2:7] = 1
x[2:7] = np.arange(5)
A = x.reshape((3,3))
print(A)
B = A.copy()
B[1,1]=0
print(B)
print(A[1,1])
print(A,B)
#A.A #has to be type other than A to become A
A.T
print(A)
AA=A/2
print(A)
print(AA)
AA=A[A>2]
print(AA)
print(A)
m = np.arange(12)
m = m.reshape((3,4))
v = np.arange(10,31,10)
print(m)
print(v)
print(m+2)
print(v+2)

M = np.Mat([[1,3],[2,3]])
print(M)
print(M.__str__())
print(M.__repr__())
M= np.Array(range(9))
print(M)
M=np.reshape(M,(3,3))
print(M) # tutaj wykona typowe mnożenie dla typu 'array', czyli element przez element.
M=np.Mat(np.reshape(M,(3,3)))
print(M)
M*M
print(M*M) # tutaj wykona typowe mnożenie dla typu 'matrix', czyli wymnoży jedną macierz przez drugą
M=np.Mat(np.reshape(M,(3,3)))
print(M,"M")
M*M
print(M*M,"M*M")
v = np.Mat(range(1,4))
print(v)
v = v.T
print(A,"A")
print(v,"v")
print(A*v)
print(v*A)
print(M)
R=np.Mat([[1.1,3.9],[2.4,3.5]])
for i in R: print(i)
