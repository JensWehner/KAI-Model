#!/usr/bin/env python
import numpy as np
import sys

#a0=n a1=m
def Txx(n,m):
        a=np.arctan((m+0.5)/(n-0.5))-np.arctan((m-0.5)/(n-0.5))   
        return a


def Txy(n,m):
        a=0.5*np.log(((m-0.5)**2+(n+0.5)**2)/((m-0.5)**2+(n-0.5)**2))
        return a

f=1.0-2.0/np.pi
g=1.0+2.0/np.pi-16.0/(np.pi)**2


Nrange=[100,300,500,700]


for N in Nrange:
       

        n0=np.arange(-N,N,1.0)
        m0=np.arange(-N,N,1.0)
        nn,mm=np.meshgrid(n0,m0)
        n=nn.flatten()
        m=mm.flatten()
        h=Txx(n,m)*(Txx(n,m)-0.5*(Txx(n+1.0,m)+Txx(n-1.0,m)))
        i=Txy(n,m)*(Txy(n,m)-0.5*(Txy(n,m+1.0)+Txy(n,m-1.0)))
        deltaEx=np.sum(f*h+g*i)
        deltaEy=np.sum(g*h+f*i)
        print "deltaEx",deltaEx,"deltaEy",deltaEy
                        
