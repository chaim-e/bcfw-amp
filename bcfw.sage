#  Copyright (c) 2022 Chaim Even-Zohar, Tsviqa Lakrec, Ran J Tessler
#
#  This file is the Sage script "bcfw.sage".
#
#  "bcfw.sage" is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  "bcfw.sage" is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with "bcfw.sage".  If not, see <https://www.gnu.org/licenses/>.

import itertools

MARKERS = '?123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ'

SIMPLIFY = True

# chord diagram for a BCFW cell and matrix construction
class CD:

    # constructor
    def __init__(self, n, chords):

        # verify chord diagrams rules
        assert all([isinstance(x,tuple) for x in chords])
        assert all([len(x) == 4 for x in chords])
        assert all([x[0]+1==x[1] and x[2]+1==x[3] for x in chords])
        assert all([not x[0]<y[0]<x[2]<y[2] for x in chords for y in chords])
        assert len(set([x[0] for x in chords])) == len(chords)
        assert min(map(min,chords+[[1]])) >= 1
        assert max(map(max,chords+[[n-1]])) <= n-1

        # initialize fields
        self.n = n
        self.k = len(chords)
        self.chords = chords = sorted(chords)
        self.tails = dict(zip([x[0] for x in chords], range(len(chords))))

    # textual representation and hash keys
    def __repr__(self):
        return 'CD(%d:%s)' % (self.n,','.join([('%s%s-%s%s' % 
            tuple([MARKERS[x] for x in c])) for c in self.chords])) 
    def __hash__(self):
        return hash(self.__repr__())
    def __eq__(self, other):
        return self.__repr__() == other.__repr__()

    # compute domino matrix
    def construct_matrix(self, params, ring = QQ):

        # setup variables 
        n = self.n
        k = self.k
        tails = self.tails
        chords = self.chords
        assert len(params) == k
        assert all([len(x)==4 for x in params])
        
        # init
        M = {n : vector([0]*k, ring)}
        rot = -1

        # main loop
        for m in range(n-1,0,-1):
            
            # fill
            if m not in M:
                M[m] = vector([0]*k, ring)
            
            # tail
            if m in tails and m-1 not in tails:
                for i in range(0,n):
                    if m+i not in tails:
                        break
                    M[m+i] += M[m+i+1]*params[tails[m+i]][1]   # y

            # head
            for c in chords:
                if c[2] == m:
                    for i in M:
                        M[i][tails[c[0]]:] *= -1
                        if i > c[1]: M[i] *= -1
                    M[c[1]] = vector([int(j==tails[c[0]]) for j in range(k)],ring)   # inc
                    rot *= -1
                    M[c[2]] += M[c[1]] * params[tails[c[0]]][2]   # x
                    M[c[3]] += M[c[2]] * params[tails[c[0]]][3]   # x
                    before = [x for x in M if x<c[1]]
                    before = max(before) if before else n
                    M[before] += M[c[1]] * params[tails[c[0]]][0]*rot**(before==n)   # y

        # domino form is ready
        return Matrix([M[i] for i in range(1,n+1)]).T

# iterator over chord diagrams given n,k 
def CDnk(n,k):
    for tails in itertools.combinations(range(1,n-3),k):
        for heads in itertools.combinations_with_replacement(
                [x for x in range(3,n-1) if x-1 not in tails],k):
            try:
                chords = []
                stack = []
                for i in range(1,n-1):
                    for j in range(heads.count(i)):
                        chords.append(stack.pop() + (i,i+1))
                    if i in tails:
                        stack.append((i,i+1))
                yield CD(n,sorted(chords))
            except IndexError:
                pass

# iterator over chord diagrams for all k given n
def CDn(n):
    for k in range(0,n-3):
        for cd in CDnk(n,k):
            yield cd

# a functionary: a formula in twistors            
class functionary:
    def __init__(self, twistor = None, terms = None, quadratic = None,
                 simplify = SIMPLIFY):
        if terms:
            self.twistor = False
            self.terms = terms
        elif twistor:
            self.terms = False
            assert isinstance(twistor, tuple)
            assert len(twistor) == 4
            assert twistor[0] <= twistor[1] <= twistor[2] <= twistor[3]
            self.twistor = twistor
        elif quadratic:
            self.twistor = False
            a,b,c,d,e,f,g = quadratic
            self.terms = [(+1,[functionary((a,c,d,g)),
                               functionary((b,e,f,g))]),
                          (-1,[functionary((a,e,f,g)),
                               functionary((b,c,d,g))])]
        if simplify:
            self.simplify()

    def simplify(self):
        if self.terms:
            self.terms = [(reduce(lambda x,y:x*y,[f.terms[0][0] for f in factors 
                                                  if f.terms and len(f.terms) == 1], coef),
                           reduce(lambda x,y:x+y, [f.terms[0][1] for f in factors 
                                                   if f.terms and len(f.terms) == 1], []) + 
                           [f for f in factors if f.twistor or len(f.terms) > 1]) 
                          for coef,factors in self.terms 
                          if coef != 0 
                          and not any([factor.twistor and len(set(factor.twistor))<4 
                                       for factor in factors])
                          and not any([factor.terms and len(factor.terms)==0 
                                       for factor in factors])
                         ]
        
    def __repr__(self):
        if self.twistor:
            a,b,c,d = self.twistor
            return '<%s%s%s%s>' % (
                MARKERS[a],MARKERS[b],MARKERS[c],MARKERS[d])
        return '(%s)' % ''.join([''.join(['+' if coef==+1 else
                                          '-' if coef==-1 else
                                          '+(%s)' % coef] + 
                                         [repr(factor) for factor in factors])
                           for coef,factors in self.terms])
        
    def pre(self, i):
        if self.twistor:
            return functionary(tuple([j if j<i else j+1 
                                      for j in self.twistor]))
        return functionary(terms = [(coef,[factor.pre(i) 
                                           for factor in factors]) 
                                    for coef,factors in self.terms])
    
    def upper(self, i, n):
        if self.twistor:
            a,b,c,d = self.twistor
            if d not in [n-1,n] or not a<b<c<d:
                return functionary((a,b,c,d))
            elif d == n and c < n-1:
                return functionary(terms =
                    [(+1,[functionary((a,b,c,n)),functionary((i,i+1,n-2,n-1))]),
                     (-1,[functionary((a,b,c,n-1)),functionary((i,i+1,n-2,n))]),
                     (+1,[functionary((a,b,c,n-2)),functionary((i,i+1,n-1,n))])])
            elif d == n-1 and c < n-1:
                return functionary(terms =
                    [(+1,[functionary((a,b,c,n-1)),functionary((i,i+1,n-2,n))]),
                     (-1,[functionary((a,b,c,n-2)),functionary((i,i+1,n-1,n))])])
            elif d == n and c == n-1:
                return functionary(terms = 
                  [(+1,[functionary((i,i+1,n-2,n-1)),functionary(terms =
                    [(+1,[functionary((a,b,n-1,n)),functionary((i,i+1,n-2,n))]),
                     (-1,[functionary((a,b,n-2,n)),functionary((i,i+1,n-1,n))])])
                       ])])                                    
        return functionary(terms = [(coef,[factor.upper(i,n) 
                                           for factor in factors]) 
                                    for coef,factors in self.terms])
    
    def lower(self, i, n):
        if self.twistor:
            a,b,c,d = self.twistor
            if i+1 not in [c,d] or not a<b<c<d:
                return functionary((a,b,c,d))
            elif d == i+1:
                return functionary(terms =
                    [(+1,[functionary((a,b,c,i+1)),functionary((i,n-2,n-1,n))]),
                     (-1,[functionary((a,b,c,i)),functionary((i+1,n-2,n-1,n))])])
            elif c == i+1 and d == n:
                return functionary(terms =
                    [(+1,[functionary((a,b,i+1,n)),functionary((i,n-2,n-1,n))]),
                     (-1,[functionary((a,b,i,n)),functionary((i+1,n-2,n-1,n))])])
        return functionary(terms = [(coef,[factor.lower(i,n) 
                                           for factor in factors]) 
                                    for coef,factors in self.terms])
    
# generate a functionary separating between two chord diagrams
def separate(cd1, cd2):

    # initialize variables
    c1 = cd1.chords
    c2 = cd2.chords
    assert c1 != c2
    assert cd1.n == cd2.n
    n = cd1.n
    markers1 = {i for x in c1 for i in x}
    absent1 = set(range(1,n)) - markers1
    markers2 = {i for x in c2 for i in x}
    absent2 = set(range(1,n)) - markers2

    # case (A)
    if not c1:
        j = max(markers2)-1
        c = min([c for c in c2 if j+1 in c])
        return functionary(twistor = (c[0],c[1],c[2],n))
    if not c2:
        j = max(markers1)-1
        c = min([c for c in c1 if j+1 in c])
        return functionary(twistor = (c[0],c[1],c[2],n))        
    
    # case (B)
    absent = set(absent1) & set(absent2)
    if absent:
        i = min(absent)
        cd1r = CD(n-1,[tuple([j if j<i else j-1 for j in c]) for c in c1])
        cd2r = CD(n-1,[tuple([j if j<i else j-1 for j in c]) for c in c2])
        sep = separate(cd1r, cd2r)
        return sep.pre(i)
        
    # case (C)
    if n-1 in markers2 and n-1 not in markers1:
        c = min([c for c in c2 if n-1 in c])
        return functionary(twistor = (c[0],c[1],c[2],n))
    if n-1 in markers1 and n-1 not in markers2:
        c = min([c for c in c1 if n-1 in c])
        return functionary(twistor = (c[0],c[1],c[2],n))
    
    # not (A)/(B)/(C)
    assert n-1 in markers1 and n-1 in markers2
    i = min([c[0] for c in c1 if n-1 in c])
    j = min([c[0] for c in c2 if n-1 in c])
    
    # case (D)
    if i != j:
        #print('case (D)',i,j)
        i,j = min(i,j),max(i,j)
        return functionary(quadratic=(i,i+1,j,j+1,n-2,n-1,n))
    
    # case (E)
    c1r = [c for c in c1 if c[0] > i] 
    c1l = [c for c in c1 if c[0] < i] 
    c2r = [c for c in c2 if c[0] > i] 
    c2l = [c for c in c2 if c[0] < i] 
    if c1r != c2r:
        sep = separate(CD(n,c1r), CD(n,c2r))
        return sep.upper(i,n)        
    
    # case (F)
    if c1l != c2l:
        sep = separate(CD(n,c1l), CD(n,c2l))
        return sep.lower(i,n)
    
    # not (A)/(B)/(C)/(D)/(E)/(F)
    assert False

# list all separator functionaries for a given n            
def list_separators(n):
    cds = list(CDn(n))
    for cd1,cd2 in itertools.combinations(cds, 2):
        sep = separate(cd1, cd2)
        print('%s %s %s' % (str(cd1).ljust((n-4)*6),
                            str(cd2).ljust((n-3)*6),sep))
