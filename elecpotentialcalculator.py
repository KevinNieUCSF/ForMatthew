import scipy
import numpy as np
import sys
import math

def main():
    CaIon,Carbon = arrayinit()
    CaIondist, Carbondist, Kdist = distcalc(CaIon, Carbon)
    CaPot, TotPot = potential(CaIondist, Carbondist, Kdist)
    print("CaPot: " + str(CaPot) + "kJ/mol")
    print("TotPot: " + str(TotPot) + "kJ/mol")
    return

def distcalc(CaIon, Carbon):
    a = -1.65
    b = -4.96
    K = ([0,0,a],[0,0,b])
    CaIondist = []
    Carbondist = []
    for i in CaIon:
        for center in K:
            dist = math.sqrt((i[0]-center[0])**2+(i[1]-center[1])**2+(i[2]-center[2])**2)
            CaIondist.append(dist)
    for i in Carbon:
        for center in K:
            dist = math.sqrt((i[0]-center[0])**2+(i[1]-center[1])**2+(i[2]-center[2])**2)
            Carbondist.append(dist)
    Kdist = math.sqrt((a - b)**2)
    return CaIondist, Carbondist, Kdist

def arrayinit():
    IonVert = (0,-3.3,-6.6,3.3,6.6)
    CaIon = []
    Carbon = []
    for i in IonVert:
        CaIon.append([-1.65,-1.65,i])
        CaIon.append([-1.65,1.65,i])
        CaIon.append([1.65,-1.65,i])
        CaIon.append([1.65,1.65,i])
    for i in IonVert:
        Carbon.append([3.31,0,i])
        Carbon.append([-3.31,0,i])
        Carbon.append([0,3.31,i])
        Carbon.append([0,-3.31,i])
    return CaIon, Carbon

def potential(CaIondist, Carbondist, Kdist):
    # U = k * q1q2 / (distance)
    # k = 8.976E9 in units of Nm**2/C**2
    k = (8.976E9)/4 #because of dielectric of 4
    A = 1e-10
    Avo = 6.022E23 / (1000*0.239)
    C = 1.602E-19
    CaPot = 0
    CarbPot = 0
    KPot = k * 1 * 1 * C**2 / (Kdist * A)
    TotPot = 0
    for dist in CaIondist:
        CaPot += k * 1 * -0.51 * C**2 / (dist * A)
    for dist in Carbondist:
        CarbPot += k * 1 * 0.51 * C**2 / (dist * A)
    print("Distance Between Ions: " + str(Kdist) + "A")
    CaPot = CaPot * Avo
    CaPot += KPot * Avo
    print(KPot * Avo)
    TotPot = CaPot + CarbPot * Avo
    return CaPot, TotPot


main()
