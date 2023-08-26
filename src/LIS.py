from bisect import bisect_left

def LIS(seq):

    path = [-1] * len(seq)
    ind  = [0]
    sub  = [seq[0]]

    for i, x in enumerate(seq[1:], 1):

        if sub[-1] < x:
            path[i] = ind[-1]
            sub.append(x)
            ind.append(i)

        else:
            ins = bisect_left(sub, x)
            path[i] = ind[ins - 1] if ins else -1
            sub[ins] = x
            ind[ins] = i

    ans = []
    i = ind[-1]

    while i >= 0:

        # ans.append(seq[i])
        ans.append(i)
        i = path[i]

    return ans[::-1]
