#!/usr/bin/env python
# coding: utf-8

# In[1]:


import math
from interval import interval, imath
from Polynomials import *


# In[2]:


class Taylor:
    def __init__(self, poly:Poly, remainder:interval, domain:interval):
        self.P=poly
        self.I=interval(remainder)
        self.D=interval(domain)

    @property
    def order(self):
        return self.P.order
    
    def __repr__(self):
        return 'Taylor({0}, {1}, {2})'.format(self.P, self.I, self.D)
    
    def __add__(self, other):
        if isinstance(other,(int,float)):
            other=Taylor(Poly([other],0),0,self.D)
        D=self.D&other.D
        return Taylor(self.P+other.P,self.I+other.I,D)
    
    def __radd__(self,other):
        return self+other
    
    def __sub__(self, other):
        if isinstance(other,(int,float)):
            other=Taylor(Poly([other],0),0,self.D)
        D=self.D&other.D
        return Taylor(self.P-other.P,self.I-other.I,D)
    
    def __rsub__(self,other):
        if isinstance(other,(int,float)):
            other=Taylor(Poly([other],0),0,self.D)
        return other-self
    
    # full multiplication, n is the result order
    def fullmul(self,other,n='full',c='s'):
        if isinstance(other,(int,float)):
            if other==0:
                return 0
            else:
                return Taylor(other*self.P,other*self.I,self.D)
        D=self.D&other.D
        P=self.P*other.P
        I1=self.P.bound(D,c)*other.I
        I2=self.I*other.P.bound(D,c)
        I=I1+I2+self.I*other.I
        if n=='full':
            return Taylor(P,I,D)
        elif isinstance(n,int) and n>0 and n<=(self.order+other.order):
            return Taylor(P,I,D).truncation(n,c)
        else:
            raise Exception("n must be 'full' or integer between 0 ~ self.order+other.order")
    
    # do not use subdivision
    def __mul__(self,other):
        if isinstance(other,(int,float)):
            if other==0:
                return 0
            else:
                return Taylor(other*self.P,other*self.I,self.D)
        return self.fullmul(other,n=max(self.order,other.order),c='s')
    
    
    def __rmul__(self,other):
        return self*other
    
    def bound(self,c='s'):
        return self.P.bound(self.D,c)+self.I
    
    def __call__(self,c='s'):
        return self.bound(c)
    
    def truncation(self,n,c='s'):
        P=self.P.truncation(n)
        I=(self.P-P).bound(self.D,c)+self.I
        return Taylor(P,I,self.D)
    
    def __truediv__(self, other):
        if isinstance(other,(int,float)):
            if other==0:
                raise Exception("Cannot dividied by zero!")
            else:
                return Taylor(self.P/other,self.I/other,self.D)
        else:
            raise Exception("Can only do scalar division!")
    
    def one(self):
        return Taylor(Poly.one(self.order),self.I,self.D)
    
    def zero(self):
        return Taylor(Poly.zero(self.order),self.I,self.D)
    
    def __eq__(self,other):
        if self.P==other.P and self.I==other.I and self.D==other.D:
            return True
    
    def square(self):
        return self*self
    
    def __pow__(self,n):
        if n==0:
            return self.one()
        elif n%2==0:
            if n==2:
                return self.square()
            else:
                p=int(n/2)
                return self.__pow__(p).square()
        elif n%2==1:
            if n==1:
                return self
            else:
                p=int((n-1)/2)
                return self*self.__pow__(p).square()
        else:
            raise Exception("n must be a positive integer")
    
    # 0~n+1 coefficients of exp sin cos
    def coeff(f,n,x,flag=0):
        if flag==1: # n+1 coefficients
            h=imath
        else: h=math
        if f=='exp':
            return h.exp(x)/math.factorial(n)
        if f=='sin':
            i=n%4
        if f=='cos':
            i=(n+1)%4
        if i==0:
            return h.sin(x)/math.factorial(n)
        elif i==1:
            return h.cos(x)/math.factorial(n)
        elif i==2:
            return -h.sin(x)/math.factorial(n)
        elif i==3:
            return -h.cos(x)/math.factorial(n)
    
    # elementary functions with lagrange remainder
    def element(f,n,D):
        coef=[]
        for i in range(n+1):
            coef.append(Taylor.coeff(f,i,0))
        coef1=Taylor.coeff(f,n+1,D,1)
        I=coef1*(D**(n+1))
        return Taylor(Poly(coef,n),I,D)
            
    #composition of a function f with a Taylor model T, f(T(x))
    def compose(self,f,n=10):
        if isinstance(f,Poly):
            return f(self)
        elif f=='exp' or f=='sin' or f=='cos':   
            Bf=self.bound()
            Mg=Taylor.element(f,n,Bf)
            M1=self.zero()
            M1.P.coef[1:]=self.P.coef[1:]
            M=Mg.P(M1)
            return Taylor(M.P,M.I+Mg.I,self.D)
    





