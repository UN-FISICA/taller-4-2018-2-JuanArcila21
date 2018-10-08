# -*- coding: utf-8 -*-
import math as m
from scipy import optimize
class Derivada:
    def __init__(self,f,metodo="adelante",dx=0.001):
        self.f=f
        self.metodo=metodo
        self.dx=dx
    def calc(self,x):
        if self.metodo=="adelante":
            return (self.f(x+self.dx)-self.f(x))/self.dx
        elif self.metodo=="central":
            return ((self.f(x+(self.dx)/2)-self.f(x-(self.dx)/2))/self.dx)
        elif self.metodo=="extrapolada":
            f1=(self.f(x+(self.dx)/2)-self.f(x-(self.dx)/2))/self.dx
            f2=(self.f(x+(self.dx)/4)-self.f(x-(self.dx)/4))/(self.dx/2)
            return (4*f2-f1)/3
        elif self.metodo=="segunda":
            return ((self.f(x+self.dx)+self.f(x-self.dx)-2*self.f(x))/(self.dx)**2)
        else:
            return "elija un metodo valido"
        
class Zeros:
    def __init__(self,f,metodo,error=1e-4,max_iter=100):
        self.f=f
        self.metodo=metodo
        self.error=error
        self.max=max_iter
    def zero(self,vi):
        if self.metodo=="newton":
            fprima=Derivada(self.f,"extrapolada")
            x=float(vi)
            i=0
            while i<=self.max:
                x=x-(self.f(x)/fprima.calc(x))
                if abs(self.f(x))<self.error:
                    i=self.max+1
                else:
                    i+=1
            i=0
            return x
        elif self.metodo=="bisectriz":
            a=min(vi[0],vi[1])
            b=max(vi[0],vi[1])
            if (self.f(a)>=0 and self.f(b)>=0) or (self.f(a)<=0 and self.f(b)<=0):
                print ("Se necesita otra tupla para este metodo")
                return
            else:
                i=0
                while i<=self.max:
                    if self.f((a+b)/2)<0:
                        a=(a+b)/2
                    elif self.f((a+b)/2)>0:
                        b=(a+b)/2
                    elif abs(self.f((a+b)/2))<self.error:
                        i=self.max+1 
                    i+=1
                i=0
                return (a+b)/2
        elif self.metodo=="interpolacion":
            i=0
            x=vi[0]
            x2=vi[1]
            y=self.f(x)
            y2=self.f(x2)
            if abs(self.f(x))<self.error:
                i=self.max+1 
            while i<=self.max:
                x=(((x2-x)/(y-y2))*y)+x
                y=self.f(x)
                if abs(self.f(x))<self.error:
                    i=self.max+1 
                i+=1
            return x
        elif self.metodo=="newton-sp":
            return optimize.newton(self.f,vi,maxiter=self.max,tol=self.error)
        elif self.metodo=="fsolve-sp":
            return optimize.fsolve(self.f,vi,xtol=self.error)
        elif self.metodo=="brentq-sp":
            return optimize.brentq(self.f,vi[0],vi[1],maxiter=self.max,xtol=self.error)

       
if __name__=="__main__":
    print ("La derivada del seno evaluada en pi/2 es {} por el metodo adelante,{} por el metodo central, {} por el metodo extrapolado. Con {} el valor de su segunda derivada en pi/2 ".format(Derivada(m.sin,"adelante").calc(m.pi/2),Derivada(m.sin,"central").calc(m.pi/2),Derivada(m.sin,"extrapolada").calc(m.pi/2),Derivada(m.sin,"segunda").calc(m.pi/2)))
    print ("Estas derivadas se calcularon con un dx de 0.001")
    print ("El cero de la funcion seno al rededor del intervalo [-pi/2,pi/2] para metodos de dos datos inciales o al rededor de pi/4 para metodos de un solo dato inicial son: {} por el metodo de newton, {} por el metodo de la bisectriz, {} por el metodo de interpolacion lineal, {} por el metodo de newton en scipy, {} por el metodo fsolve y {} por el metodo brentq".format(Zeros(m.sin,"newton").zero(m.pi/4),Zeros(m.sin,"bisectriz").zero((-m.pi/2,m.pi/2)),Zeros(m.sin,"interpolacion").zero((-m.pi/2,m.pi/2)),Zeros(m.sin,"newton-sp").zero(m.pi/4),Zeros(m.sin,"fsolve-sp").zero(m.pi/4),Zeros(m.sin,"brentq-sp").zero((-m.pi/2,m.pi/2))))
    print ("Los anteriores calculos fueron hechos con un máximo de 100 iteraciones, un error mínimo de 1e-4")