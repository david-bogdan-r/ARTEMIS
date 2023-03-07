def SeqScore(seq1,seq2,i,j,weight=0.5):

    X = [s1[i] for s1 in seq1]
    Y = [s2[j] for s2 in seq2]

    res = 0

    for x in X:
        for y in Y:
            if x!='-' and x==y:
                res += weight

    return res/(len(X)*len(Y))
                

def GapGlobalAlign(v,w,sigma,eps,scoremat):

    lenv, lenw = len(scoremat),len(scoremat[0])

    s = {}
    backtrack = {}
    s[(0,0,0)] = 0
    s[(1,0,0)] = 0
    s[(2,0,0)] = 0

    for i in range(1,lenv+1):
        
        s[(0,i,0)] = -sigma-eps*(max(i-2,0))
        s[(1,i,0)] = -sigma-eps*(max(i-2,0))
        s[(2,i,0)] = -sigma-eps*(max(i-2,0))
        backtrack[(0,i,0)] = 'v'
        backtrack[(1,i,0)] = 'l'
        backtrack[(2,i,0)] = 'l'
        
    for j in range(1,lenw+1):
        
        s[(0,0,j)] = -sigma-eps*(max(j-2,0))
        s[(1,0,j)] = -sigma-eps*(max(j-2,0))
        s[(2,0,j)] = -sigma-eps*(max(j-2,0))
        backtrack[(0,0,j)] = 'u'
        backtrack[(1,0,j)] = 'u'
        backtrack[(2,0,j)] = 'h'


    for i in range(1,lenv+1):
        for j in range(1,lenw+1):

            score = scoremat[i-1][j-1] + SeqScore(v,w,i-1,j-1)

            if s[(0,i-1,j)]-eps >= s[(1,i-1,j)]-sigma:
                s[(0,i,j)] = s[(0,i-1,j)]-eps
                backtrack[(0,i,j)] = 'v'
            else:
                s[(0,i,j)] = s[(1,i-1,j)]-sigma
                backtrack[(0,i,j)] = 'vu'

            if s[(2,i,j-1)]-eps >= s[(1,i,j-1)]-sigma:
                s[(2,i,j)] = s[(2,i,j-1)]-eps
                backtrack[(2,i,j)] = 'h'
            else:
                s[(2,i,j)] = s[(1,i,j-1)]-sigma
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

def OutputGapAlign1(backtrack,v,i,j):

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

def OutputGapAlign2(backtrack,w,i,j):

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
    
    backtrack, score = GapGlobalAlign(seq1,seq2,sigma,epsilon,mat)

    aseq1 = OutputGapAlign1(backtrack,seq1,len(seq1[0]),len(seq2[0]))
    aseq2 = OutputGapAlign2(backtrack,seq2,len(seq1[0]),len(seq2[0]))

    for x in aseq1:
        print(x)
    for x in aseq2:
        print(x)