#!/usr/bin/env python
# coding: utf-8

# In[1]:


import math
from interval import interval, imath, fpu


# In[2]:


class Poly:
    def __init__(self, coef:list, order:int):
        self.order=order
        if len(coef)<(order+1):
            for i in range(order+1-len(coef)):
                coef.append(0) 
        self.coef=coef

    def __repr__(self):
        return 'Poly(%s, %d)'%(self.coef, self.order)
   
    def fixshape(x, y):
        order=max(x.order,y.order)
        x1=x.coef[:]
        y1=y.coef[:]
        if x.order<order:
            for i in range(order-x.order):
                x1.append(0)
        else:
            for i in range(order-y.order):
                y1.append(0)
        return Poly(x1,order), Poly(y1,order), order
        
    def fixshape1(x, y):
        order=x.order+y.order
        x1=x.coef[:]
        y1=y.coef[:]
        for i in range(order-x.order):
            x1.append(0)
        for j in range(order-y.order):
            y1.append(0)
        return Poly(x1,order), Poly(y1,order), order   
    
    def firstnonzero(self):
        order=self.order
        nonzero=order+1
        for i in range(order+1):
            if self.coef[i]!=0:
                nonzero=i
                break
        return nonzero
    
    def __add__(self, other):
        coef=[]
        if isinstance(other,(int,float)):
            other=Poly([other],0)
        self1, other1, order=self.fixshape(other)
        for i in range(order+1):
            coef.append(self1.coef[i]+other1.coef[i])
        return Poly(coef,order)
    
    def __radd__(self,other):
        return self+other
    
    def __sub__(self, other):
        coef=[]
        if isinstance(other,(int,float)):
            other=Poly([other],0)
        self1, other1, order=self.fixshape(other)
        for i in range(order+1):
            coef.append(self1.coef[i]-other1.coef[i])
        return Poly(coef,order)
    
    def __rsub__(self,other):
        if isinstance(other,(int,float)):
            other=Poly([other],0)
        return other-self
    
    def __mul__(self, other):
        coef=[]
        if isinstance(other,(int,float)):
            if other==0:
                return 0
            else:
                other=Poly([other],0)
        self1, other1, order=Poly.fixshape1(self, other)
        coef.append(self1.coef[0]*other1.coef[0])
        for i in range(1,order+1):
            s=0
            for j in range(i+1):
                s+=self1.coef[j]*other1.coef[i-j]
            coef.append(s)
        return Poly(coef,order)
    
    def __rmul__(self,other):
        return self*other
    
    def __truediv__(self, other):
        if isinstance(other,(int,float)):
            if other==0:
                raise Exception("Cannot dividied by zero!")
            else:
                coef=[]
                for i in range(self.order+1):
                    coef.append(self.coef[i]/other)
                return Poly(coef, self.order)
        else:
            raise Exception("Can only do scalar division!")

    def zero(n):
        coef=[]
        coef.append(0)
        for i in range(n):
            coef.append(0)
        return Poly(coef,n)
    
    def one(n):
        coef=[]
        coef.append(1)
        for i in range(n):
            coef.append(0)
        return Poly(coef,n)
    
    def __eq__(self,other):
        self1, other1, order=Poly.fixshape(self, other)
        return self1.coef==other1.coef
    
    def square(self):
        order=self.order*2
        for a in range(self.order):
            self.coef.append(0)
        coef=[]
        coef.append((self.coef[0])**2)
        for i in range(1,order+1):
            coef.append(self.squarecoef(i,self.coef))
        return Poly(coef,order)        
    
    def squarecoef(self,i,coef):
        s=0
        iend=int((i-2+i%2)/2)
        for j in range(iend+1):
            s+=coef[j]*coef[i-j]
        s=2*s
        if i%2==0:  
            s+=coef[int(i/2)]**2
        return s
    
    def __pow__(self,n=int):
        order=self.order*n
        for a in range(self.order):
            self.coef.append(0)
        if n<0:
            raise Exception("n must ≥ 0")
        elif n==0:
            return Poly.one(order)
        elif n%2==0:
            if n==2:
                return self.square()
            else:
                p=int(n/2)
                return Poly(self.__pow__(p).square().coef,order)
        elif n==1:
            return self
        else:
            p=int((n-1)/2)
            return Poly((self*self.__pow__(p).square()).coef,order)
    
    def exp(self):
        order=self.order
        coef=[]
        coef.append(math.exp(self.coef[0]))
        for i in range(1,order+1):
            s=0
            for j in range(i):
                s+=(i-j)*self.coef[i-j]*coef[j]
            coef.append(s/i)
        return Poly(coef,order)
    
    
    def log(self):
        order=self.order
        coef=[]
        coef.append(math.log(self.coef[0]))
        for i in range(1,order+1):
            s=0
            for j in range(1,i):
                s+=(i-j)*self.coef[j]*coef[i-j]
            coef.append((self.coef[i]-s/i)/self.coef[0])
        return Poly(coef,order)
    
    def sincos(self):
        order=self.order
        scoef=[];ccoef=[]
        scoef.append(math.sin(self.coef[0]))
        ccoef.append(math.cos(self.coef[0]))
        for i in range(1,order+1):
            s=0;c=0
            for j in range(1,i+1):
                s+=j*self.coef[j]*ccoef[i-j]
                c-=j*self.coef[j]*scoef[i-j]
            scoef.append(s/i)
            ccoef.append(c/i)
        return Poly(scoef,order), Poly(ccoef,order)
    
    def sin(self):
        return self.sincos()[0]
    
    def cos(self):
        return self.sincos()[1]
    
    def tan(self):
        order=self.order
        coef1=[]
        coef2=[]
        t=math.tan(self.coef[0])
        coef1.append(t)
        coef2.append(t**2)
        for i in range(1,order+1):
            coef1.append(self.tancoef(i,self.coef,coef2))
            coef2.append(self.squarecoef(i,coef1))
        return Poly(coef1,order)
    
    def tancoef(self,i,coef,coef2):
        s=0
        for j in range(i):
            s+=(i-j)*coef[i-j]*coef2[j]
        s=self.coef[i]+s/i
        return s
    
    def diffPoly(self):
        order=self.order
        coef=[]
        coef.append(self.coef[1])
        for i in range(2,order+1):
            coef.append(i*self.coef[i])
        coef.append(0)
        return Poly(coef,order)
        
    def intePoly(self,x):
        order=self.order
        coef=[]
        coef.append(x)
        for i in range(1,order+1):
            coef.append(self.coef[i-1]/i)
        return Poly(coef,order)
    
    #Horner's rule
    def __call__(self,dx):
        order=self.order
        s=self.coef[order]
        for i in range(order-1,-1,-1):
            s=s*dx+self.coef[i]
        return s
    
    # f(g(x)),self=f,other=g
    def compose(self,other):
        return self(other)

    def truncation(self,n):
        return Poly(self.coef[:n+1],n)    
    
    def bound_naive(self,x):
        B=self.coef[0]
        n=self.order
        for i in range(1,n+1):
            B+=self.coef[i]*(x**i)
        return B
    
    def bound(self,x,c='s'):
        if c=='naive' or c=='n':
            return self.bound_naive(x)
        elif c=='Horner' or c=='H':
            return self(x)
        elif c=='sub' or c=='s':
            return self.subdivision(x)
        elif c=='better' or c=='b':
            return self(x)&self.bound_naive(x)
        else:
            raise Exception("The options mush be naive(n), Horner(H), sub(s) or better(b)")
    
    #subdivision k times
    def binarychop(x):
        I=[]
        for i in x:
            sup=i[0].sup
            inf=i[0].inf
            l=(sup-inf)/2
            I.append(interval([inf,inf+l]))
            I.append(interval([inf+l,sup]))
        return I
    
    def subdivision(self,x):
        a=self(x) 
        I=Poly.binarychop([x])
        b=self(I[0])|self(I[1])    
        k=(a[0].sup-a[0].inf)-(b[0].sup-b[0].inf)
        while k>=0.001:
            a=b
            I=Poly.binarychop(I)
            b=self(I[0])
            for i in I:
                b=b|self(i)
            k=(a[0].sup-a[0].inf)-(b[0].sup-b[0].inf)
        return b  






