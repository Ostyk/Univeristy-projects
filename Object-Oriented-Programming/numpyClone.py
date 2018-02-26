# -*- coding: utf-8 -*-
#Michal Ostyk-Narbtt
class Array(object):
    def __init__(self,lista):
        if type(lista) == range:
            self.lista = list(lista)
            self.dim=len(lista)
        elif type(lista)==float:
            self.lista = lista
        else:
            self.lista = lista
            self.dim=len(lista)
    def __len__(self):
         return len(self.lista)
    def __str__(self):
        return self.__class__.__name__+"("+str(self.lista)+')'
    def __repr__(self):
        return Array.__str__(self)
    depth = lambda L: isinstance(L, list) and max(map(depth, L))+1
    def __add__(self,other):
        temp = []
        if type(other) is int:
            if depth(self.lista)==1:
                for j in range(self.dim):
                    temp.append(self.lista[j] + other)
                return Array(temp)
            else:
                for i in range(self.dim):
                    for j in range(len(self.lista[i])):
                        self.lista[i][j]=self.lista[i][j]+other
                    return Array(self.lista)
        else:
            if depth(other)==1:
                for j in range(len(self.lista)):
                    temp.append(self.lista[j] + other.lista[j])
                return Array(temp)
            else:
                for i in range(self.dim):
                    for j in range(len(self.lista[i])):
                        self.lista[i][j]=self.lista[i][j]+other[i][j]
                return Array(self.lista)
    def __mul__(self,other):
        if type(self) and type(other)!=Array:
            out=[]
            for i in range(len(other)): out.append(other[i])
            other=out
            return Mat([[sum(ele_a*ele_b for ele_a, ele_b in zip(row_a, col_b))
                     for col_b in other] for row_a in self.lista])
        else:
            if type(other) is  int or type(other) is float:
                l=[]
                for i in range(self.dim):
                    l.append(self.lista[i]*other)
                return Array(l)
            else:
                l=[]
                for i in range(self.dim):
                    l.append(self.lista[i]*other[i])
                return Array(l)
    def __truediv__(self,other):
        l=[]
        if type(other) is  int or type(other) is float:
            if depth(self.lista)!=1:
                l=ones((self.dim,len(self.lista[0])))
                for i in range(self.dim):
                    for j in range(len(self.lista[i])):
                        l[i][j]=self.lista[i][j]/other
                return l
            else:
                for i in range(self.dim):
                    l.append(self.lista[i]/other)
                return Array(l)
        else:
            for i in range(self.dim):
                l.append(self.lista[i]//other[i])
            return Array(l)
    def __setitem__(self,key,value):
        if isinstance(key,tuple):
            self.lista[key[0]][key[1]]=value
        elif isinstance(key,list):
            self.lista[key[0]:key[1]]=value
        elif isinstance(key,slice):
            el=list(range(len(self.lista)))[key]
            if type(value) is int:
                for i in range(len(el)):
                    self.lista[el[i]]=value
            else:
                for i in range(len(el)):
                    self.lista[el[i]]=value[i]
        else:
            self.lista[key]=value
    def __getitem__(self,key):
        if isinstance(key,slice):
            return Array([self[ii] for ii in range(*key.indices(len(self)))])
        elif isinstance(key,int):
            if key < 0 :
                key += len(self)
            if key < 0 or key >= len(self) :
                raise ValueError
            return self.lista[key]
        elif isinstance(key,tuple):
            a,b=int(key[0]),int(key[1])
            return [self[i][a:b+1] for i in range(a,b+1)][0][0]
        elif isinstance(key,list):
            return Array([self.lista[key[i]] for i in range(len(key))])
        else:
            return key[:]
    def reshape(self,t):
        a,b=t[0],t[1]
        return Array([self.lista[a*i : a*(i+1)] for i in range(b)])
    @property
    def T(self):
        x=self.lista
        if depth(x)==1:
            return Array(self.lista)
        else:
            y=[list(tup) for tup in zip(*x)]
            return self.__class__((Array(y).lista))
    def copy(self):
        return Array(self.lista)
    def __gt__(self,other):
        flat=[j for i in self.lista for j in i]
        new=sorted([i for i in flat if i>other])
        return Array(new)

import math as m
pi = m.pi
depth = lambda L: isinstance(L, list) and max(map(depth, L))+1
def zeros(n):
    if type(n) is tuple:
        if len(n)==2:
            return Array([[0]*n[0] for i in range(0,n[1])])
        elif len(n)==1:
            return Array([[0]*n for i in range(0,n)])
        else:
            raise TypeError
    elif type(n) is int:
        return Array([0 for i in range(n)])
def ones(n):
    if type(n) is tuple:
        if len(n)==2:
            return Array([[1]*n[0] for i in range(0,n[1])])
        elif len(n)==1:
            return Array([[1]*n for i in range(0,n)])
        else:
            raise TypeError
    elif type(n) is int:
        return Array([1 for i in range(n)])
def concatenate(a):
    return Array(sorted(list(a[0])+list(a[1])))
def arange(start, stop=None, step=1):
    if stop!=None:
        n = int(round((stop - start)/float(step)))
        if n > 1:
            return Array([float("{:3g}".format(start+step*i)) for i in range(n)])
        else:
            return Array([])
    else:
        return Array([i for i in range(0,start)])
def reshape(lista,t):
    a,b=t[0],t[1]
    out=[]
    for i in range(len(lista)): out.append(lista[i])
    if depth(out)!=1:
        return out
    else:
        x=[out[a*i : a*(i+1)] for i in range(b)]
        return Array(x)
def linspace(start, stop, num=50, endpoint=True):
    if num == 1:
        return Array([start])
    else:
        if endpoint==False:
            step = (stop - start) / (num - 1)
        else:
            step = (stop - start) / num
        return Array([float("{:3g}".format(start+step*i)) for i in range(num+1)])

class Mat(Array):
    def __init__(self,lista):
        if type(lista)==Array:
            x=[]
            for i in range(len(lista)):
                x.append(lista[i])
            self.lista=x
            super(Mat,self).__init__(x)
        else:
            self.lista=list(lista)
            super(Mat,self).__init__(lista)
    def __getitem__(self,key):
        return self.lista[key]
    def __iter__(self):
        q=[]
        for i in range(len(self.lista)):
            for j in range(len(self.lista[i])):
                q.append((j,i))
        return iter(q)
    def __next__(self):
        return next(self.lista)
    def __mul__(self,other):
        out=[]
        for i in range(len(other)): out.append(other[i])
        other=out
        return Mat([[sum(ele_a*ele_b for ele_a, ele_b in zip(row_a, col_b))
                 for col_b in other] for row_a in self.lista])
    @property
    def T(self):
        x=self.lista
        if depth(x)==1:
            return Mat([(self.lista)])
        else:
            y=[list(tup) for tup in zip(*x)]
            return self.__class__((Mat(y).lista))
    @property
    def A(self):
        return Array(self.lista)
