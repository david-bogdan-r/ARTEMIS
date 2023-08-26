def globalAlign(
    rseq:'str', qseq:'str', scoremat
) -> 'tuple[str, str]':

    n = len(rseq)
    m = len(qseq)

    F = [
        [0]*(m+1) 
        for i in range(n+1)
    ]

    T = [
        [0, *[2]*m], 
        *[
            [1, *[0]*m] 
            for i in range(n)
        ]
    ]

    for i in range(1, n+1):

        fi   = F[i]
        ti   = T[i]
        fi_1 = F[i-1]
        si_1 = scoremat[i-1]

        for j  in range(1, m+1):

            mi = 0
            mf =  fi_1[j-1] + si_1[j-1]

            if fi_1[j] > mf:
                mf = fi_1[j]
                mi = 1

            if fi[j-1] > mf:
                mf = fi[j-1]
                mi = 2

            fi[j] = mf
            ti[j] = mi

    rali = ''
    qali = ''

    i = n
    j = m

    while i or j:

        t = T[i][j]

        if i and j and t == 0:
            i -= 1
            j -= 1
            rali += rseq[i]
            qali += qseq[j]

        elif i and t == 1:
            i -= 1
            rali += rseq[i]
            qali += '-'

        else:
            j -= 1
            rali += '-'
            qali += qseq[j]

    rali = rali[::-1]
    qali = qali[::-1]

    return rali, qali

# ---

def seqScore(seq1, seq2, i, j, weight=0.5):

    X = [s1[i] for s1 in seq1]
    Y = [s2[j] for s2 in seq2]

    res = 0

    for x in X:
        for y in Y:
            if x!='-' and x==y:
                res += weight

    return res/(len(X)*len(Y))

def gapGlobalAlign(v, w, sigma, eps, scoremat):

    lenv, lenw = len(scoremat),len(scoremat[0])

    s = {}
    backtrack = {}
    s[(0,0,0)] = 0
    s[(1,0,0)] = 0
    s[(2,0,0)] = 0

    for i in range(1,lenv+1):

        s[(0,i,0)] = 0 #-sigma-eps*(max(i-2,0))
        s[(1,i,0)] = 0 #-sigma-eps*(max(i-2,0))
        s[(2,i,0)] = 0 #-sigma-eps*(max(i-2,0))
        backtrack[(0,i,0)] = 'v'
        backtrack[(1,i,0)] = 'l'
        backtrack[(2,i,0)] = 'l'

    for j in range(1,lenw+1):

        s[(0,0,j)] = 0 #-sigma-eps*(max(j-2,0))
        s[(1,0,j)] = 0 #-sigma-eps*(max(j-2,0))
        s[(2,0,j)] = 0 #-sigma-eps*(max(j-2,0))
        backtrack[(0,0,j)] = 'u'
        backtrack[(1,0,j)] = 'u'
        backtrack[(2,0,j)] = 'h'

    for i in range(1,lenv+1):
        for j in range(1,lenw+1):

            if i==lenv or j==lenw:
                cursigma, cureps = 0, 0 
            else:
                cursigma, cureps = sigma, eps

            score = scoremat[i-1][j-1]# + SeqScore(v,w,i-1,j-1)

            if s[(0,i-1,j)]-cureps >= s[(1,i-1,j)]-cursigma:
                s[(0,i,j)] = s[(0,i-1,j)]-cureps
                backtrack[(0,i,j)] = 'v'
            else:
                s[(0,i,j)] = s[(1,i-1,j)]-cursigma
                backtrack[(0,i,j)] = 'vu'

            if s[(2,i,j-1)]-cureps >= s[(1,i,j-1)]-cursigma:
                s[(2,i,j)] = s[(2,i,j-1)]-cureps
                backtrack[(2,i,j)] = 'h'
            else:
                s[(2,i,j)] = s[(1,i,j-1)]-cursigma
                backtrack[(2,i,j)] = 'hl'

            if s[(1,i-1,j-1)]+score >= s[(0,i,j)] and s[(1,i-1,j-1)]+score >= s[(2,i,j)]:
                s[(1,i,j)] = s[(1,i-1,j-1)]+score
                backtrack[(1,i,j)] = 'd'
            elif s[(0,i,j)] >= s[(1,i-1,j-1)]+score and s[(0,i,j)] >= s[(2,i,j)]:
                s[(1,i,j)] = s[(0,i,j)]
                backtrack[(1,i,j)] = 'l'
            elif s[(2,i,j)] >= s[(1,i-1,j-1)]+score and s[(2,i,j)] >= s[(0,i,j)]:
                s[(1,i,j)] = s[(2,i,j)]
                backtrack[(1,i,j)] = 'u'

    return backtrack, s[(1,lenv,lenw)]


def outputGapAlign1(backtrack,v,i,j):

    lenv = len(v)

    res = ['' for _ in range(lenv)]

    r = 1

    while i or j:

        if   (r,i,j) in backtrack and backtrack[(r,i,j)] in ('v','vu'):
            for k in range(lenv):
                res[k] = v[k][i-1] + res[k]
            if backtrack[(r,i,j)] == 'vu': r += 1
            i -= 1
        elif (r,i,j) in backtrack and backtrack[(r,i,j)] in ('h','hl'):
            for k in range(lenv):
                res[k] = '-' + res[k]
            if backtrack[(r,i,j)] == 'hl': r -= 1
            j -= 1
        elif (r,i,j) in backtrack and backtrack[(r,i,j)] == 'l':
            r -= 1
        elif (r,i,j) in backtrack and backtrack[(r,i,j)] == 'u':
            r += 1
        else:
            for k in range(lenv):
                res[k] = v[k][i-1] + res[k]
            i,j = i-1, j-1

    return res

def outputGapAlign2(backtrack,w,i,j):

    lenw = len(w)

    res = ['' for _ in range(lenw)]

    r = 1

    while i or j:

        if   (r,i,j) in backtrack and backtrack[(r,i,j)] in ('v','vu'):
            for k in range(lenw):
                res[k] = '-' + res[k]
            if backtrack[(r,i,j)] == 'vu': r += 1
            i -= 1
        elif (r,i,j) in backtrack and backtrack[(r,i,j)] in ('h','hl'):
            for k in range(lenw):
                res[k] = w[k][j-1] + res[k]
            if backtrack[(r,i,j)] == 'hl': r -= 1
            j -= 1
        elif (r,i,j) in backtrack and backtrack[(r,i,j)] == 'l':
            r -= 1
        elif (r,i,j) in backtrack and backtrack[(r,i,j)] == 'u':
            r += 1
        else:
            for k in range(lenw):
                res[k] = w[k][j-1] + res[k]
            i,j = i-1, j-1

    return res


if __name__ == "__main__":

    sigma   = 0
    epsilon = 0

    seq1 = ("AACGGU",) # tuple of sequences of length lenv
    seq2 = ("ACGU",)   # tuple of sequences of length lenw

    mat = [[2,0,0,0],
           [0,0,0,0],
           [0,5,0,0],
           [0,0,5,0],
           [0,0,0,7],
           [0,0,0,0]]
    
    backtrack, score = gapGlobalAlign(seq1,seq2,sigma,epsilon,mat)

    aseq1 = outputGapAlign1(backtrack,seq1,len(seq1[0]),len(seq2[0]))
    aseq2 = outputGapAlign2(backtrack,seq2,len(seq1[0]),len(seq2[0]))

    for x in aseq1:
        print(x)
    for x in aseq2:
        print(x)
