# -*- coding: utf-8 -*-
"""
Created on Thu Apr 20 14:15:21 2023

@author: Martyna
"""

from math import * 
import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('dane', type=argparse.FileType('r'))
args = parser.parse_args()
lines = args.dane.readlines()


class XYZ_BLH:
    def __init__(self, x, y, z):
        self.X = x
        self.Y = y
        self.Z = z
    
    def zamianaxyz_blh(self):
        X = self.X
        Y = self.Y
        Z = self.Z
        a = 6378137
        e2 = 0.00669438002290
    
        def Np(f,a,e2):
            N = a / np.sqrt(1 - e2*np.sin(f)**2)
            return(N)
        
        def Mp(f,a,e2):
            M = (a*(1-e2))/ (np.sqrt(1 - e2*np.sin(f)**2))**3
            return(M)
        
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
            
        with open("wyniki-XYZ_BLH.txt", "w") as file:
            file.write(f"f = {f}, l = {l}, h = {h}")
        


class  BLH_XYZ:
    def __init__(self, fi,lam,ha):
        self.f = fi
        self.l = lam
        self.h = ha
        
    def zamianablh(self):
        f = self.f
        l = self.l
        h = self.h
        a = 6378137
        e2 = 0.00669438002290
       
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
            
        with open("wyniki-BLH_XYZ.txt", "w") as file:
            file.write(f"X = {X}, Y = {Y}, Z = {Z}")
       

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
        
        with open("wyniki-BL_U2000.txt", "w") as file:
            file.write(f"Xgka2000 = {Xgka2000}, Ygka2000 = {Ygka2000}")
        
        
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
       
        with open("wyniki-BL_U1992.txt", "w") as file:
            file.write(f"Xgka1992 = {Xgka1992}, Ygka1992 = {Ygka1992}")
       

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
        
        with open("wyniki-XYZ_NEU.txt", "w") as file:
            file.write(f"dneu = {dneu}")
        

 
try:
    if args.dane.name =="dane-XYZ_BLH.txt":
        for line in lines:
            data = line.strip().split()
            X = float(data[0])
            Y = float(data[1])
            Z = float(data[2])
        obl_XYZ_BLH = XYZ_BLH(X, Y, Z)
        XYZ_BLH.zamianaxyz_blh(obl_XYZ_BLH)
        
    if args.dane.name =="dane-XYZ_NEU.txt":
        for line in lines:
            data = line.strip().split()
            Xa = float(data[0])
            Ya = float(data[1])
            Za = float(data[2])
            Xb = float(data[3])
            Yb = float(data[4])
            Zb = float(data[5])
            
        blh = XYZtoNEU(Xa, Ya, Za, Xb, Yb, Zb)
        XYZtoNEU.zamiana_na_NEU(blh)
        
    if args.dane.name =="dane-BLH_XYZ.txt":
        for line in lines:
            data = line.strip().split()
            f = float(data[0])
            l = float(data[1])
            h = float(data[2])
            
        obl = BLH_XYZ(f, l, h)
        BLH_XYZ.zamianablh(obl)
        
    if args.dane.name =="dane-BL_U2000.txt":
        for line in lines:
            args = line.strip().split() 
            b = float(args[0])
            l = float(args[1])
            a = float(args[2])
            e2 = float(args[3])
            
        bl = BL(b,l,a,e2)
        BL.zamiana_na_2000(bl)

    if args.dane.name =="dane-BL_U1992.txt":
        for line in lines:
            args = line.strip().split() 
            b = float(args[0])
            l = float(args[1])
            a = float(args[2])
            e2 = float(args[3])
        
        bl2 = BL(b,l,a,e2)
        BL.zamiana_na_1992(bl2)
        
except AttributeError:
    print('Podałes zły atrubut')
except ValueError:
    print('Źle zapisałes plik z danymi')
except RuntimeWarning:
    print('Sprawdź czy wpisałe odpowiednie dane')
except IndexError:
    print('Sprawdz czy wprowadziłe odpowiednią ilosć danych')
#except:
 #   print('Cos poszło nie tak')
finally:
    print('Do zobaczenia!')