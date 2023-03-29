
def one_to_one_matches(xset, yset, col, eps):
    from itertools import count
    X = sorted(list(zip(xset[col], count())), key=lambda x:x[0])
    Y = sorted(list(zip(yset[col], count())), key=lambda x:x[0])
    ix, iy = 0, 0
    matchx, matchy = [], []
    while ix<len(X) and iy<len(Y):
        if abs(X[ix][0]-Y[iy][0]) < eps:
            matchx.append(X[ix][1])
            matchy.append(Y[iy][1])
            ix += 1
            iy += 1
        elif X[ix][0] < Y[iy][0]:
            ix += 1
        else:
            iy += 1
    return xset.iloc[matchx], yset.iloc[matchy] 
