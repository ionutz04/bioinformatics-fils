'''
make a software implementation in which the native design of the global sequance aligment is used 2 sequances.
S1= "ACGTGCATGA"
S2 = "ACGT"


a[i,j] = max of:
    a[i-1, j-1] + S[xi,xj]
    a[i-1, j] + g
    a[i,j-1] + g

S :
    +2, xi = xj
    -1, xi != xj


        A    C   G   T   G   C   A   T   G   A
    0  -2  -4   -6   -8  -10 -12 -14 -16 -18 -20
A   -2  2  -1   -1  

C   -4 -1    2   -1

G   -6 -1   -1   2

T   -8 -1   -1   -1
'''
def S(a, b):
    if a == b:
        return 2
    else: return -1
def main():
    import numpy as np

    S1= "ACGTGCATGA"
    S2 = "ACGT" 
    n = len(S1)+1
    m = len(S2)+1
    a = np.zeros([n,m])
    g=-2
    
    for i in range(n):
        a[i,0] = i * -2

    for j in range(m):
        a[0,j] = j * -2
    # print(a)


    for i in range(1, n):
        for j in range(1, m):
            a[i,j] = max(a[i-1,j-1] + S(S1[i-1], S2[j-1]), a[i-1,j] + g, a[i,j-1] + g)

    print(a.T)



if __name__ == "__main__":
    main()
