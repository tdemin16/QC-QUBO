import random as r

n = 3

A = [ [ 0 for i in range(n) ] for j in range(n) ] 
B = [ [ 0 for i in range(n) ] for j in range(n) ]

for i in range(n):
    for j in range(i, n):
        if(i != j):
            k = r.randint(1, 100)
            A[i][j] = k
            A[j][i] = k

            k = r.randint(1, 100)
            B[i][j] = k
            B[j][i] = k
        pass
    pass

f = open("test_7.txt", "w")
f.write(str(n) + '\n'+'\n')

for i in range(n):
    for j in range(n):
        f.write(str(A[i][j]))
        if(j < n-1):
            f.write(' ')
            pass
        pass
    
    f.write('\n')

f.write('\n')

for i in range(n):
    for j in range(n):
        f.write(str(B[i][j]))
        if(j < n-1):
            f.write(' ')
            pass
        pass
    
    if i < n-1:
        f.write('\n')

f.close()

print(A)
print(B)