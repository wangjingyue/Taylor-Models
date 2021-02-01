#!/usr/bin/env python
# coding: utf-8

# In[1]:


import math
from interval import interval, imath
from Polynomials import *
import time


# In[2]:


class Taylor:
    def __init__(self, poly:Poly, remainder:interval, x0:interval, domain:interval):
        poly.coef=[interval(x) for x in poly.coef]
        self.P=poly
        self.x0=interval(x0)
        self.I=interval(remainder)
        self.D=interval(domain)

    @property
    def order(self):
        return self.P.order
    
    def __repr__(self):
        return 'Taylor({0}, {1}, {2}, {3})'.format(self.P, self.I, self.x0, self.D)
    
    def normalize(self):
        order=self.order
        r=(self.D[0].sup-self.D[0].inf)/2
        c=(self.D-self.x0).midpoint
        P=self.P(Poly([0,1],1)*r+c)
        P=Poly(P.coef,order)
        D=interval([-1,1])
        return Taylor(P,self.I,0,D)
    
    def meanvalue(self):
        c=self.D.midpoint
        diff=self.P.diffPoly()
        return self.P(c)+diff(self.D)*(self.D-c)+self.I
    
    def reexpand(self,x0):
        order=self.order
        P=self.P(Poly([0,1],1)+x0)
        P=Poly(P.coef,order)
        D=self.D-x0
        return Taylor(P,self.I,self.x0,D)
    
    def affine(d:interval):
        x0=interval((d[0].sup+d[0].inf)/2)
        x1=interval((d[0].sup-d[0].inf)/2)
        coef=[x0,x1]
        return Taylor(Poly(coef,1),0,0,interval[-1,1])
    
    def set_order(self,n):
        return Taylor(Poly(self.P.coef,n),self.I,self.x0,self.D)
    
    def set_x0(self,x0):
        return Taylor(self.P,self.I,x0,self.D)
    
    def set_domain(self,D):
        return Taylor(self.P,self.I,self.x0,D)
    
    def __neg__(self):
        return 0-self
    
    def __add__(self, other):
        if isinstance(other,(int,float,interval)):
            other=Taylor(Poly([other],0),0,0,self.D)
        D=self.D&other.D
        return Taylor(self.P+other.P,self.I+other.I,self.x0,D)
    
    def __radd__(self,other):
        return self+other
    
    def __sub__(self, other):
        if isinstance(other,(int,float,interval)):
            other=Taylor(Poly([other],0),0,0,self.D)
        D=self.D&other.D
        return Taylor(self.P-other.P,self.I-other.I,self.x0,D)
    
    def __rsub__(self,other):
        if isinstance(other,(int,float,interval)):
            other=Taylor(Poly([other],0),0,0,self.D)
        return other-self
    
    # full multiplication, n is the result order
    def fullmul(self,other,n='full',c='r'):
        if isinstance(other,(int,float,interval)):
            if other==0:
                return 0
            else:
                return Taylor(other*self.P,other*self.I,self.x0,self.D)
        D=self.D&other.D
        P=self.P*other.P
        I1=self.P.bound(D-self.x0,c)*other.I
        I2=self.I*other.P.bound(D-self.x0,c)
        I=I1+I2+self.I*other.I
        if n=='full':
            return Taylor(P,I,self.x0,D)
        elif isinstance(n,int) and n>0 and n<=(self.order+other.order):
            return Taylor(P,I,self.x0,D).truncation(n)
        else:
            raise Exception("n must be 'full' or integer between 0 ~ self.order+other.order")
    
    # do not use subdivision
    def __mul__(self,other):
        if isinstance(other,(int,float,interval)):
            if other==0:
                return 0
            else:
                return Taylor(other*self.P,other*self.I,self.x0,self.D)
        return self.fullmul(other,n=max(self.order,other.order),c='r')
    
    
    def __rmul__(self,other):
        return self*other
    
    def bound(self,c='r'):
        if c=='LDB':
            return self.LDB()
        elif c=='QFB':
            return self.QFB()
        elif c=='m' or c=='meanvalue':
            return self.meanvalue()
        elif c=='BNB':
            return self.BNB()
        else:
            return self.P.bound(self.D-self.x0,c)+self.I
    
    def __call__(self,c='H'):
        return self.bound(c)
    
    def truncation(self,n):
        P=self.P.truncation(n)
        Z=self.zero()
        Z.P.coef[n+1:]=self.P.coef[n+1:]
        I=Z.P.bound(self.D-self.x0,'r')+self.I
        return Taylor(P,I,self.x0,self.D)
    
    def __truediv__(self, other):
        if isinstance(other,(int,float,interval)):
            if other==0:
                raise Exception("Cannot dividied by zero!")
            else:
                return Taylor(self.P/other,self.I/other,self.x0,self.D)
        elif isinstance(other,Taylor):
            return other.inv()*self
    
    def one(self):
        return Taylor(Poly.one(self.order),self.I,self.x0,self.D)
    
    def zero(self):
        return Taylor(Poly.zero(self.order),self.I,self.x0,self.D)
    
    def __eq__(self,other):
        if self.P==other.P and self.x0==other.x0 and self.I==other.I and self.D==other.D:
            return True
    
    def doublefactorial(n):
         if n <= 0:
            return 1
         else:
            return n*Taylor.doublefactorial(n-2)
    
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
    def coeff(f,n,x,flag=1):
        if flag==1: # n+1 coefficients
            h=imath
        elif flag==0:
            h=math
        if f=='exp':
            return h.exp(x)/math.factorial(n)
        if f=='log':
            if n==0:
                return h.log(x)
            else:
                return (-1)**(n-1)/((x**n)*n)
        if f=='sqrt':
            if n==0:
                return h.sqrt(x)
            else:
                return h.sqrt(x)*((-1)**(n-1))*Taylor.doublefactorial(2*n-3)/(math.factorial(n)*(2**(n))*x**(n))
        if f=='inv':
            return (-1)**n/x**(n+1)
        if f=='sin' or f=='cos':
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
        if f=='sinh' or f=='cosh':
            if f=='sinh':
                i=n%2
            if f=='cosh':
                i=(n+1)%2
            if i==0:
                return h.sinh(x)/math.factorial(n)
            elif i==1:
                return h.cosh(x)/math.factorial(n)
        if f=='atan':
            if x==0 or x==interval[0]:
                if n%2==0:
                    return 0
                else: 
                    return (-1)**((n-1)/2)/n
            else:
                if n==0:
                    return h.atan(x)
                elif n%2==0:
                    return (-1)**(n-1)/((1+x**2)**(int(n/2))*n)*h.sin(n*(h.atan(1/x)))
                else:
                    return (-1)**(n-1)/(((1+x**2)**(int(n/2))*h.sqrt(1+x**2))*n)*h.sin(n*(h.atan(1/x)))
            
        
    
    # elementary functions with lagrange remainder   
    def element(f,n,x0,D):
        coef=[]
        for i in range(n+1):
            coef.append(Taylor.coeff(f,i,x0,1))
        P=Poly(coef,n)
        coef1=Taylor.coeff(f,n+1,D,1)
        f=eval('imath.'+f)
        if coef1[0].inf>=0 or coef1[0].sup<=0:
            a=interval(D[0].inf)
            b=interval(D[0].sup)
            I1=f(a)-P.bound_best(a-x0,['n','H'])
            I2=f(b)-P.bound_best(b-x0,['n','H'])
            I0=f(interval(x0))-P.bound(interval(0))
            I=(interval.hull((I1,I2,I0)))&(coef1*((D-x0)**(n+1)))
        else:
            I=coef1*((D-x0)**(n+1))
        return Taylor(P,I,x0,D)
    
    def element_inv(n,x0,D):
        coef=[]
        for i in range(n+1):
            coef.append(Taylor.coeff('inv',i,x0,1))
        P=Poly(coef,n)
        coef1=Taylor.coeff('inv',n+1,D,1)
        if coef1[0].inf>=0 or coef1[0].sup<=0:
            a=interval(D[0].inf)
            b=interval(D[0].sup)
            I1=interval(1)/a-P.bound(a-x0)
            I2=interval(1)/b-P.bound(b-x0)
            I0=interval(1)/interval(x0)-P.bound(interval(0))
            I=(interval.hull((I1,I2,I0)))&(coef1*((D-x0)**(n+1)))
        else:
            I=coef1*((D-x0)**(n+1))
        return Taylor(P,I,x0,D)
    
    #composition of a function f with a Taylor model T, f(T(x))
    def compose(self,f,n=10):
        if isinstance(f,Poly):
            return f(self)
        else:  
            x0=self.P.coef[0]
            Bf=self.bound_best(['n','H'])
            if f=='inv':
                Mg=Taylor.element_inv(n,x0,Bf)
            else:
                Mg=Taylor.element(f,n,x0,Bf)
            M1=self.zero()
            M1.P.coef[1:]=self.P.coef[1:]
            M=Mg.P(M1)
            return Taylor(M.P,M.I+Mg.I,self.x0,self.D)
    
    def exp(self):
        return self.compose('exp',self.order)
    
    def log(self):
        return self.compose('log',self.order)
    
    def sqrt(self):
        return self.compose('sqrt',self.order)
    
    def inv(self):
        return self.compose('inv',self.order)
    
    def sin(self):
        return self.compose('sin',self.order)
    
    def cos(self):
        return self.compose('cos',self.order)
    
    def tan(self):
        return self.sin()*(self.cos().inv())
    
    def sinh(self):
        return self.compose('sinh',self.order)
    
    def cosh(self):
        return self.compose('cosh',self.order)
    
    def tanh(self):
        return self.sinh()*(self.cosh().inv())
    
    def atan(self):
        return self.compose('atan',self.order)
    
    def asin(self):
        return (x/((1-self.square()).sqrt())).atan()
    
    def acos(self):
        return imath.pi/2-self.asin()
    
    def LDB_MIN(self,e=1e-3,max_iter=5):
        d=1
        Dn=self.D
        Dm=self.x0
        Pm=self.P
        bound=interval(0)
        n_iter=0
        while d>e and n_iter<max_iter:
            x0=(Dn-Dm).midpoint
            c=Dn.midpoint
            Pm=Pm(Poly([0,1],1)+x0)
            linear=Poly(Pm.coef[0:2],Pm.order)
            non_linear=Pm-linear
            L1=linear.coef[1].midpoint[0].inf
            I1=linear.bound(Dn-c,'n')
            Ih=non_linear.bound_best(Dn-c,['n','H'])
            bound=I1[0].inf+Ih
            d=bound[0].sup-bound[0].inf
            n_iter+=1
            if L1==0:
                break
            elif L1>0:
                new_hi=min(Dn[0].inf+(d/abs(L1)),Dn[0].sup)
                Dm=Dn
                Dn=interval[Dn[0].inf,new_hi]
            else:
                new_lo=max(Dn[0].sup-(d/abs(L1)),Dn[0].inf)
                Dm=Dn
                Dn=interval[new_lo,Dn[0].sup]
        hi=self.P.bound_best(self.D-self.x0,['n','H'])[0].sup
        return interval[bound[0].inf,hi]+self.I
    
    def LDB(self,e=1e-3,max_iter=5):
        a=self.LDB_MIN()
        b=-(0-self).LDB_MIN()
        return a&b
    
    def QFB_MIN(self):
        P=self.P
        bound_tm=self.bound_best(['n','H'])
        if P.coef[2]<=interval(0):
            return bound_tm
        else:
            x0=-P.coef[1]/(2*P.coef[2])
            x=Poly([0,1],self.order)
            Qx0=(x-x0)*P.coef[2]*(x-x0)
            bound_qfb=(P-Qx0).bound_best(self.D-self.x0,['n','H'])
            hi=P.bound_best(self.D-self.x0,['n','H'])[0].sup
            bound_qfb=interval[bound_qfb[0].inf,hi]+self.I
            bound=bound_qfb&bound_tm
            return bound
    
    def QFB(self):
        a=self.QFB_MIN()
        b=-(0-self).QFB_MIN()
        return a&b
    
    def subdomain(self,n):
        D=[]
        sup=self.D[0].sup
        inf=self.D[0].inf
        l=(sup-inf)/n
        for i in range(n):
            D.append(interval[inf+i*l,inf+(i+1)*l])
        return D
    
    def subbound(self,n,b,normalize=True):
        D=self.subdomain(n)
        print(D)
        truncation=self.order
        if normalize==True:
            bound=[self.set_domain(i).normalize().bound(b) for i in D] 
        else:
            bound=[self.set_domain(i).bound(b) for i in D]
        print(bound)
        resnew=bound[0]
        for i in bound:
            resnew=resnew|i
        return resnew
    
    def bound_best(self,methods):
        bound=self.bound(methods[0])
        for m in methods:
            bound=bound&self.bound(m)
        return bound
    
    def BNB(self,h=1e-5):
        e=h
        res0=self.bound_best(['n','H','LDB','QFB'])
        dom=Poly.binarychop([self.D])
        bound=[self.set_domain(i).bound_best(['n','H','LDB','QFB']) for i in dom]
        res=bound[0]|bound[1]
        loc=abs(res0[0].inf-res[0].inf)
        hic=abs(res0[0].sup-res[0].sup)
        while loc>e or hic>e:
            hi=[]
            lo=[]
            for i in bound:
                hi.append(i[0].sup)
                lo.append(i[0].inf)
            a=min(lo)
            minD=[i for i,x in enumerate(lo) if x==a]
            hii=[]
            for i in range(len(minD)):
                hii.append(hi[minD[i]])
            minD=[minD[hii.index(max(hii))]]
            b=max(hi)
            maxD=[i for i,x in enumerate(hi) if x==b]
            loo=[]
            for i in range(len(maxD)):
                loo.append(lo[maxD[i]])
            maxD=[maxD[loo.index(min(loo))]]
            minmax=list(set(minD+maxD))
            minmax.sort()
            j=0
            for i in range(len(minmax)):
                chop=Poly.binarychop([dom[minmax[i]+j]])
                dom[minmax[i]+j:minmax[i]+1+j]=Poly.binarychop([dom[minmax[i]+j]])
                bound[minmax[i]+j:minmax[i]+1+j]=[self.set_domain(i).bound_best(['n','H','LDB','QFB']) for i in chop]
                j=j+1
            resnew=bound[0]
            for i in bound:
                resnew=resnew|i
            loc=abs(resnew[0].inf-res[0].inf)
            hic=abs(res[0].sup-resnew[0].sup)
            res=resnew
            e=h*(res[0].sup-res[0].inf)
        return res






