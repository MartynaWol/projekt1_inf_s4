# -*- coding: utf-8 -*-
"""
Created on Thu Apr 20 14:15:21 2023

@author: Martyna
"""

from math import * 
import numpy as np


class XYZ_BLH:
    def __init__(self,X,Y,Z):
        self.X = X
        self.Y = Y
        self.Z = Z
    def zamianaxyz_blh(self):
        x = self.X
        y = self.Y
        z = self.Z
        
        def Np(f,a,e2): #promien krzywizny N
            N = a / np.sqrt(1-e2 * np.sin(f)**2)
            return(N)

        def Mp (f, a, e2):
            M = a * (1-e2)/np.sqrt((1-e2*(np.sin(f))**2)**3)
            return(M)

        
        p = np.sqrt(X**2 + Y**2)
        f = np.arctan(Z / (p * (1-e2)))
        while True:
            N = Np(f,a,e2)
            h = p/np.cos(f) - N
            fs = f
            f = np.arctan(Z / (p * (1 - e2 * N / (N + h))))
            if abs(fs - f) < (0.000001/206265):
                break 
        l = np.arctan2(Y,X)
        xyz_blh = f"{f}, {l}, {h}"
        np.savetxt("xyz_blh.txt", xyz_blh, delimiter=",", fmt="%s")
        print("f", f)
        print("l", l)
        print("h", h)

class  BLH_XYZ:
    def __init(self, f,l,h):
        self.f = f
        self.l = l
        self.h = h
    def zamianablh(self):
        f = self.f
        l = self.l
        h = self.h
        
        def Np(f,a,e2): #promien krzywizny N
            N = a / np.sqrt(1-e2 * np.sin(f)**2)
            return(N)
        
        while True:
            N = Np(f,a,e2)
            X = (N + h) * np.cos(f) * np.cos(l)
            Xp = X
            Y = (N + h) * np.cos(f) * np.sin(l)
            Z =  (N * (1 - e2) + h) * np.sin(f)
            if abs(Xp - X) < (0.000001/206265):
                break
        blh_xyz = f"{X}, {Y}, {Z}"
        np.savetxt("blh_xyz.txt", blh_xyz, delimiter=",", fmt="%s")
        print("X", X)
        print("Y", Y)
        print("Z", Z)
        
        
        

class BL:
    
    def __init__(self, B, L, A, E2):
        self.B = B
        self.L = L
        self.A = A
        self.E2 = E2
        
    def zamiana_na_2000(self):
        fia = self.B
        lambdaa = self.L
        a = self.A
        e2 = self.E2
        if (lambdaa>(13.5*np.pi/180) and (lambdaa<16.5*np.pi/180)) or (lambdaa==(13.5*np.pi/180)):
            l0 = 5
        if (lambdaa>(16.5*np.pi/180) and lambdaa<(19.5*np.pi/180)) or (lambdaa==(16.5*np.pi/180)):
            l0 = 6
        if (lambdaa>(19.5*np.pi/180) and lambdaa<(22.5*np.pi/180)) or (lambdaa==(19.5*np.pi/180)):
            l0 = 7
        if (lambdaa>(22.5*np.pi/180) and lambdaa<(25.5*np.pi/180)) or (lambdaa==(22.5*np.pi/180)):
            l0 = 8
        L0=0.3141592653589793
        
        def deg2dms(dd): #dziesietne na sto min sec
            deg=np.trunc(dd) #stopnie
            mnt=np.trunc((dd-deg)*60) #minuty
            sec=((dd-deg) *60-mnt)*60 #sekundy
            #abs jest po to aby nie było - przed minutami i sekundami
            return(deg, mnt, sec)
        
        b2=a**2*(1-e2)
        er2=(((a**2)-(b2))/(b2))
        t=np.tan(fia)
        ni2=(er2)*((np.cos(fia))**2)
        l=lambdaa-L0    
        A0=1-(e2/4)-((3*(e2)**2)/64)-((5*(e2)**3)/256)
        A2=(3/8)*((e2+((e2**2)/4)+((15*(e2)**3))/128))
        A4=(15/256)*(e2**2+((3*((e2)**3))/4))
        A6=(35*((e2)**3))/3072
        
        sigma=a*(A0*fia-A2*np.sin(2*fia)+A4*np.sin(4*fia)-A6*np.sin(6*fia))
        N=a/(np.sqrt(1-e2*(np.sin(fia))**2))
        Xgk=sigma+((l**2)/2)*N*np.sin(fia)*np.cos(fia)*(1+((l**2)/12)*((np.cos(fia))**2)*(5-(t**2)+9*(ni2)+4*((ni2)**2))+((l**4)/360)*((np.cos(fia))**4)*(61-58*(t**2)+(t**4)+270*(ni2)-330*(ni2)*(t**2)))
        Ygk=l*N*np.cos(fia)*(1+((l**2)/6)*((np.cos(fia))**2)*(1-(t**2)+(ni2))+((l**4)/120)*((np.cos(fia))**4)*(5-18*(t**2)+(t**4)+14*(ni2)-58*(ni2)*(t**2)))
        
        phi=Xgk/(a*A0)
        while True:
            sig=a*(A0*phi-A2*np.sin(2*phi)+A4*np.sin(4*phi)-A6*np.sin(6*phi))
            phi_poprzednie=phi
            phi=phi_poprzednie+(Xgk-sig)/(a*A0)
            if abs(phi-phi_poprzednie)<0.000001/206265:
                break
        N1=a/np.sqrt(1-e2*np.sin(phi)**2)
        M1=(a*(1-e2))/(np.sqrt(1-(e2*(np.sin(phi))**2))**3)
        t1=np.tan(phi)
        ni=er2*np.cos(phi)*np.cos(phi)
        fiap=phi-((Ygk**2*t1)/(2*M1*N1)*(1-(Ygk**2/(12*N1**2))*(5+3*(t1**2)+(ni**2)-9*(ni**2)*(t1**2)-4*(ni**4))+(Ygk**4)/(360*N1**4)*(61+90*(t1**2)+45*(t1**4))))
        lambdaap=L0+(Ygk)/(N1*(np.cos(phi)))*(1-(((Ygk)**2)/(6*(N1**2)))*(1+2*(t1**2)+(ni))+((Ygk**4)/(120*(N1**4)))*(5+28*(t1**2)+24*(t1**4)+6*(ni)+8*(ni)*(t1**2)))
        fiapowrot=deg2dms(np.rad2deg(fiap))
        lambdaapowrot=deg2dms(np.rad2deg(lambdaap))      
        
        Xgka2000=0.999923*Xgk
        Ygka2000=0.999923*Ygk+l0*1000000+500000
        
        print('Xgka2000=', Xgka2000, 'Ygka2000=', Ygka2000)
        
        
    def zamiana_na_1992(self):
        fia = self.B
        lambdaa = self.L
        a = 6378137
        e2 = 0.00669438002290
        L0 = 19* np.pi / 180
        
        def deg2dms(dd): #dziesietne na sto min sec
            deg=np.trunc(dd) #stopnie
            mnt=np.trunc((dd-deg)*60) #minuty
            sec=((dd-deg) *60-mnt)*60 #sekundy
            #abs jest po to aby nie było - przed minutami i sekundami
            return(deg, mnt, sec)
        
        b2=a**2*(1-e2)
        er2=(((a**2)-(b2))/(b2))
        t=np.tan(fia)
        ni2=(er2)*((np.cos(fia))**2)
        l=lambdaa-L0    
        A0=1-(e2/4)-((3*(e2)**2)/64)-((5*(e2)**3)/256)
        A2=(3/8)*((e2+((e2**2)/4)+((15*(e2)**3))/128))
        A4=(15/256)*(e2**2+((3*((e2)**3))/4))
        A6=(35*((e2)**3))/3072
        
        sigma=a*(A0*fia-A2*np.sin(2*fia)+A4*np.sin(4*fia)-A6*np.sin(6*fia))
        N=a/(np.sqrt(1-e2*(np.sin(fia))**2))
        Xgk=sigma+((l**2)/2)*N*np.sin(fia)*np.cos(fia)*(1+((l**2)/12)*((np.cos(fia))**2)*(5-(t**2)+9*(ni2)+4*((ni2)**2))+((l**4)/360)*((np.cos(fia))**4)*(61-58*(t**2)+(t**4)+270*(ni2)-330*(ni2)*(t**2)))
        Ygk=l*N*np.cos(fia)*(1+((l**2)/6)*((np.cos(fia))**2)*(1-(t**2)+(ni2))+((l**4)/120)*((np.cos(fia))**4)*(5-18*(t**2)+(t**4)+14*(ni2)-58*(ni2)*(t**2)))
        
        phi=Xgk/(a*A0)
        while True:
            sig=a*(A0*phi-A2*np.sin(2*phi)+A4*np.sin(4*phi)-A6*np.sin(6*phi))
            phi_poprzednie=phi
            phi=phi_poprzednie+(Xgk-sig)/(a*A0)
            if abs(phi-phi_poprzednie)<0.000001/206265:
                break
        N1=a/np.sqrt(1-e2*np.sin(phi)**2)
        M1=(a*(1-e2))/(np.sqrt(1-(e2*(np.sin(phi))**2))**3)
        t1=np.tan(phi)
        ni=er2*np.cos(phi)*np.cos(phi)
        fiap=phi-((Ygk**2*t1)/(2*M1*N1)*(1-(Ygk**2/(12*N1**2))*(5+3*(t1**2)+(ni**2)-9*(ni**2)*(t1**2)-4*(ni**4))+(Ygk**4)/(360*N1**4)*(61+90*(t1**2)+45*(t1**4))))
        lambdaap=L0+(Ygk)/(N1*(np.cos(phi)))*(1-(((Ygk)**2)/(6*(N1**2)))*(1+2*(t1**2)+(ni))+((Ygk**4)/(120*(N1**4)))*(5+28*(t1**2)+24*(t1**4)+6*(ni)+8*(ni)*(t1**2)))
        fiapowrot=deg2dms(np.rad2deg(fiap))
        lambdaapowrot=deg2dms(np.rad2deg(lambdaap))
        
        Xgka1992 = Xgk * 0.9993 - 5300000
        Ygka1992 = Ygk * 0.9993 + 500000
        bl = f"{Xgka1992}, {Ygka1992}"
        np.savetxt("bl.txt", bl, delimiter=",", fmt="%s")
        print('Xgka1992', Xgka1992, 'Ygka1992', Ygka1992)
       
        
f = (52 + 50/60 + 00/3600)* np.pi / 180
l = (18 + 40/60 + 00/3600)* np.pi / 180
a = 6378137.000
e2 = 0.00669438002290
p1 = BL(f,l,a,e2)
BL.zamiana_na_2000(p1)
BL.zamiana_na_1992(p1)

class XYZtoNEU:
    
    def __init__(self, Xa, Ya, Za, Xb, Yb, Zb):
        self.Xa = Xa
        self.Ya = Ya
        self.Za = Za
        self.Xb = Xb
        self.Yb = Yb
        self.Zb = Zb
        
    def zamiana_na_NEU(self):
        xa = self.Xa
        ya = self.Ya
        za = self.Za
        xb = self.Xb
        yb = self.Yb
        zb = self.Zb
        a = 6378137.000
        e2 = 0.00669438002290
        
        def Np(f,a,e2):
            N = a / np.sqrt(1 - e2*np.sin(f)**2)
            return(N)
        
        #zamiana xyz na fl
        def Hirvonen(X,Y,Z,a,e2):
            p = np.sqrt(X**2 + Y**2)
            f = np.arctan(Z / (p * (1 - e2)))
            while True:
                N = Np(f,a,e2)
                h = (p/np.cos(f)) - N
                fs = f
                f = np.arctan(Z / (p* (1 - e2 * N / (N + h))))
                if abs(fs - f) < (0.000001/206265):
                    break
            l = np.arctan2(Y,X)
            return(f,l,h)
        fa, la, ha = Hirvonen(xa, ya, za, a, e2)
        fb, lb, hb = Hirvonen(xb, yb, zb, a, e2)
        
        dXYZ = [xb-xa, yb-ya, zb-za]
        
        R = np.array([[-np.sin(fa) * np.cos(la), -np.sin(la), np.cos(fa) * np.cos(la)],
                     [ -np.sin(fa) * np.sin(la),  np.cos(la), np.cos(fa) * np.sin(la)],
                     [np.cos(fa), 0, np.sin(fa)]])
        dneu = R.T @ dXYZ
        neu = f"{dneu}"
        np.savetxt("dneu.txt", neu, delimiter=",", fmt="%s")
        print('dneu = ', dneu)
        
        
X_b = 3658578.459
Y_b = 1235988.621
Z_b = 5059678.594
X_c =  3730357.4886301295
Y_c =  1215619.309205316
Z_c =  5011791.223440444

z1 = XYZtoNEU(X_c, Y_c, Z_c, X_b, Y_b, Z_b)
XYZtoNEU.zamiana_na_NEU(z1)


'''
Aby przetransformowac wspł. zapisz je jako:
    XYZ->NEU:
        X_a, Y_a, Z_a, X_b, Y_b, Z_b
    BL->U2000:
        f2, l2, a, e2
    BL->U1992:
        f92, l92, a, e2
'''


