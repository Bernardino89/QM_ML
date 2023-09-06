#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri October 2  2020

@author: Bernardino Tirri
@contact:bernardino.tirri@chimieparistech.psl.eu
"""
from gthingsCTM import gjob
from gthingsCTM import gcomBuilder
# import decimal
import numpy as np
import random

# 0 put the exponent in exact appendix position
# 1 subimit a calculation for the dimer
# 2 send a calculation for the monomer
# 3 exstract EJ, EX, ED for the dimer;
#   ED is the Dimer Energy;
#   extracts EM for the  monomer;
#   EM is the Monomer Energy
# 4 Varandas constraint
# 6 search alpha which minimise Varandas


comstructure = gcomBuilder.gcomBuilder()


def varandas0(initial_guess):
    '''
    Whith this function is possible to subimit our jobs (dimer and monomer)
    We have to minimize this function, varing the exponent in the basis set
    beacause : ED{ALPHA}, EM{ALPHA}, EJ{ALPHA}, EX{ALPHA}.
    '''

    comstructure.readout("starter.com")

    for num in range(1):
        comstructure.link0 = '%CHK=dimer_' + str(num) +\
            '.chk \n%MEM=100GB \n%Nprocs=32'
        comstructure.title = str(num) +\
            "_dimer DH-SVPD --> Hydrogen Def2-SVPD --> Oxygen"

        # if ALPHA in comstructure.appendix:
        comstructure.appendix = comstructure.appendix.replace(
            "ALPHA", '{:1.11f}'.format(abs(initial_guess[0])))
        # print ("appendix :", comstructure.appendix)
        comstructure.appendix = comstructure.appendix.replace(
            "BETA", '{:1.11f}'.format(abs(initial_guess[1])))

        # print ("appendix :", comstructure.appendix)
        job = gjob.gjob()
        job.prepare(comname="dimer_" + str(num) + ".com",
                    logname="dimer_" + str(num) + ".log",
                    linejob=comstructure.combyline()
                    )
        job.launch('g16')
        # print("=" * 80)
        # debuginfo = comstructure.appendix.split('\n')[-10:]
        # print("ALPHA <", debuginfo[-8], ">")
        # print("BETA  <", debuginfo[-4], ">")
        # print("=" * 80)
        Infodimer = job.extract(["E(PBEQIDH)", "EJ",
                                 "Ex", "Error termination"])
        # print("Infodimer\n",Infodimer)

        for row in Infodimer:
            # if some problem happend in the execution
            # of the g16 job under examination
            if ("Error termination" in row):
                print("Error termination found in " + job.logfile)
                quit()

            # grep Coulomb energy EJ for the dimer
            if ("EJ=" in row):
                energy = row.split("=")
                EJ = float(energy[-3].split()[0])
                # print("*** EJ :", EJ)
                ej.append(EJ)
                #print("*** EJ :", ej)

            # grep Exchange energy EX  for the dimer
            if ("ENTVJ=" in row):
                energy = row.split("=")
                EX = float(energy[-4].split()[0])
                # print("*** EX :", EX)
                ex.append(EX)

            # grep total energy for the dimer
            if ("E(PBEQIDH)" in row):
                split_line = row.strip().split()
                ED = float(split_line[5].replace("D", "E"))
                ed.append(ED)

            if ("E(PBEQIDH)" in row):
                split_line = row.strip().split()
                mp2_dimer = float(split_line[2].replace("D", "E"))
                MP2_D.append(mp2_dimer)
                #print("*** MP2 :", MP2_D)



        job.status = 'void'

        comstructure.link0 = '%CHK=monomer_A' \
            '.chk \n%MEM=100GB \n%Nprocs=32'
        comstructure.title = str(num) +\
            "_monomer_A  DH-SVPD --> Hydrogen Def2-SVPD --> Oxygen"
        tml = comstructure.combyline_conf()
        tml1 = comstructure.combyline_conf()
        #print ("Temporanly list:, tml")
        comstructure.cleanatoms()  # possibilita di generalizzarlo
        comstructure.addatoms(tml[0:4])
        job.prepare(comname="monomer_A.com",
                    logname="monomer_A.log",
                    linejob=comstructure.combyline()
                    )
        job.launch('g16')
        Infomonomer_A = job.extract(["E(PBEQIDH)"])
        # print("Infomonomer\n",Infomonomer)

        for row in Infomonomer_A:

            # if some problem happend in the execution
            # of the g16 job under examination
            if (row.find("Error termination") != -1):
                print("Error termination found in monomer_A.log")
                quit()

            # grep total energy for the dimer
            if (row.find("E(PBEQIDH)") != -1):
                energy = row.split("=")
                split_line = row.strip().split()
                EM_A = float(split_line[5].replace("D", "E"))
                # print("*** EM_A :", EM_A)
                ema.append(EM_A)

            if ("E(PBEQIDH)" in row):
                split_line = row.strip().split()
                mp2_mono_A = float(split_line[2].replace("D", "E"))
                MP2_A.append(mp2_mono_A)


        job.status = 'void'
        # print ("appendix :", comstructure.appendix)


        comstructure.link0 = '%CHK=monomer_B' \
            '.chk \n%MEM=100GB \n%Nprocs=32'
        comstructure.title = str(num) +\
            "_monomer_B  DH-SVPD --> Hydrogen Def2-SVPD --> Oxygen"
        #tml = comstructure.combyline_conf()
        comstructure.cleanatoms()  # possibilita di generalizzarlo
        comstructure.addatoms(tml1[4:8])
        job.prepare(comname="monomer_B.com",
                    logname="monomer_B.log",
                    linejob=comstructure.combyline()
                    )
        job.launch('g16')
        Infomonomer_B = job.extract(["E(PBEQIDH)"])
        # print("Infomonomer\n",Infomonomer)

        for row in Infomonomer_B:

            # if some problem happend in the execution
            # of the g16 job under examination
            if (row.find("Error termination") != -1):
                print("Error termination found in monomer_B.log")
                quit()

            # grep total energy for the dimer
            if (row.find("E(PBEQIDH)") != -1):
                energy = row.split("=")
                split_line = row.strip().split()
                EM_B = float(split_line[5].replace("D", "E"))
                # print("*** EM_B :", EM_B)
                emb.append(EM_B)

            if ("E(PBEQIDH)" in row):
                split_line = row.strip().split()
                mp2_mono_B = float(split_line[2].replace("D", "E"))
                MP2_B.append(mp2_mono_B)


        job.status = 'void'
        # print ("appendix :", comstructure.appendix)


    return varandas([ED, EM_A, EJ, EX, EM_B])


def varandas(x):
    ''' this function performs the operation described in Eq.1 in
    "Small Basis Set Allowing the Recovery of Dispersion Interactions
    with Double-Hybrid Functionals"  J. Chem. Theory Comput. 2019, 15, 29442953
           (ED - EM)  -  (EJ + EX)
    F =  (-------------------------)^2
           (ED - EM)  +  (EJ + EX)
    '''

    ED, EM_A, EJ, EX, EM_B = x[0], x[1], x[2], x[3], x[4]
    delta = ED - (EM_A + EM_B)
    excor = EJ + EX
    return ((delta - excor) / (delta + excor))**2


def gradient(varandas0, initial_guess, h):
    '''
    Build the gradient using the finite differences (gradient_Varandas) :
        1)
        2)
        3)
    '''
    n = initial_guess.size
    Varandas0_backward = (initial_guess * np.ones((n, n))) -\
        (np.identity(n) * h)
    # print("initial_guess * np.ones((n,n))) :",
    #       '{:1.11f} {:1.11f} '.format(initial_guess[0],
    #       initial_guess[1]),
    #       initial_guess * np.ones((n,n)))
    Varandas0_forward = (initial_guess * np.ones((n, n))) +\
        (np.identity(n) * h)
    gradient_Varandas = np.zeros(n, dtype=np.double)

    # print("#" * 80)
    # print("GRADIENT VARANDAS")
    for i in range(n):

        C1 = varandas0(Varandas0_backward[i])
        # print ("C1 :", C1)
        C2 = varandas0(Varandas0_forward[i])
        # print ("C2 :", C2)
        gradient_Varandas[i] = (C2 - C1) / (2 * h)
        # print ("i , C1, C2, grad", i , C1, C2, gradient_Varandas[i] )
    # print("#" * 80)
    return gradient_Varandas
    # print ("GRADIENT VARANDAS",  gradient_Varandas)


def convergencecriteria(varandas0, initial_guess, exponent):
    '''
    This convergence criteria allows us to describe the errors by
    a single number using vector norm.
    Because a norm of a vector or matrix is a numerical measure of its size.

    convcrit = || [alpha_old, beta_old] - [alpha, beta]||/ || [alpha,  beta] ||

            = (((alpha-alpha_old)**2 + (beta-beta_old)**2)**0.5) /
                ((alpha + alpha_old)**2)**0.5
    '''

    convcrit = (varandas0(exponent) - varandas0(initial_guess)) /\
        varandas0(exponent)

    return np.sum(convcrit)


def smart_convergencecriteria(vec):
    '''
    This convergence criteria allows us to describe the errors by
    a single number using vector norm.
    Because a norm of a vector or matrix is a numerical measure of its size.

    convcrit = || [alpha_old, beta_old] - [alpha, beta]||/ || [alpha,  beta] ||

            = (((alpha-alpha_old)**2 + (beta-beta_old)**2)**0.5) /
                ((alpha + alpha_old)**2)**0.5
    '''
    convcrit = np.square(vec)
    return np.sqrt(convcrit.sum())

#print('{0:.13s}{1:>15}{2:>17}{3:>19}{4:>21}{5:>23}   {6:25}'.format(
#    'Iteration', 'ALPHA', 'BETA', 'Convergence', 'Varandas', 'Gradient(ALPHA)',
#    'Gradient(BETA)'))

def downward(varandas0, gradient_Varandas, initial_guess, epsilon,
             MaxIter, tollerance, h):
    '''
    descendent equation :

        exponent^(Iter) = initial_guess^(Iter - 1) -
                                epsilon *  gradient_Varandas^(Iter - 1)
        the error shuld be the difference between the predicted and
        the true value
    '''
    megaloop = True
    for I in range(0, MaxIter):

        valueV0 = varandas0(initial_guess)
        GV = gradient_Varandas(varandas0, initial_guess, h)
         
        #----
        loop = True
        #info = []
        subcounter = 0
        while (loop):

            # print("+" * 80)
            # print("[[",subcounter,"]]")
            # print("initial  alpha:", initial_guess[0], "    beta:", initial_guess[1])
            exponent = initial_guess - epsilon * GV
            # print("exponent alpha:", exponent[0], "    beta:", exponent[1], "    GV[alpha]:", GV[0], "    GV[beta]:", GV[1])
            valueV1 = varandas0(exponent)
            #print ("convcrit: ", convcrit)
            #print("varandas old :", valueV0, "varandas new :", valueV1)
            # print("+" * 80)
            
            if (valueV1 > valueV0):
                loop = False                
            convergence  = smart_convergencecriteria(GV)
            convergenceV = smart_convergencecriteria(valueV0)
            convergenceV1= smart_convergencecriteria(valueV1)
            diff_var = convergenceV- convergenceV1 
            #print ("diff_var", diff_var)
            #print ("convergenceV", convergenceV)
            #print ("convergenceV1", convergenceV1)
            #print ("valueV0 - valueV", convergenceV- convergenceV1)

            #if (convergence <= tollerance) and  ( diff_var <= 1.e-9):
            #if (convergence <= tollerance) or  ( diff_var == 1.e-11):
            #    megaloop = False
            #    break 
            #if (convergence <= tollerance) or  ( diff_var == 1.e-11):
                #megaloop = False
                #loop = True
            #    megaloop = False
            #    break
            # print("(loop, megaloop) = (",loop,",",megaloop,
            #       ")   convergence ==> ", convergence)

            initial_guess = np.copy(abs(exponent))
            if (exponent[0] < 0.01e-1 or exponent[0] > 0.999) or ( exponent[1] < 0.01e-1 or exponent[1] > 0.999):
                print ("the exponent used doesn't satisfies the basis set: change the random values!")
                print ("now the random values are:")
                #exit()
                a=random.uniform(0.110238297 ,0.854761703)
                b=random.uniform(0.110238297 ,0.854761703)
                print ("New ALPHA = ",a)
                print ("New BETA  = ",b)
                initial_guess = np.array([a,b], dtype = np.double)
            valueV0 = valueV1
            # if convergence <= tollerance:
            #     megaloop = False
            # if varandas0(exponent) <= 1.00028089:
            #    megaloop = False
            subcounter += 1
            VAR_FUN.append(varandas0(exponent))
            ALPHA_EXP.append(exponent[0])
            BETA_EXP.append(exponent[1])
            ED_UPDATE.insert(len(ed),ed[-1])
            EMA_UPDATE.insert(len(ema),ema[-1])
            EMB_UPDATE.insert(len(emb),emb[-1])
            EJ_UPDATE.insert(len(ej),ej[-1])
            EX_UPDATE.insert(len(ex),ex[-1])
            MP2D_UPDATE.insert(len(MP2_D),MP2_D[-1])
            MP2A_UPDATE.insert(len(MP2_A),MP2_A[-1])
            MP2B_UPDATE.insert(len(MP2_B),MP2_B[-1])
            print('{:4d}    {:>1.11f}    {:>1.11f}    {:>1.11f}    {:>1.16f}    {:>1.16f}    {:>1.16f}'.format(
                    subcounter, exponent[0], exponent[1], convergence,
                    varandas0(exponent), GV[0], GV[1]))

            if (subcounter == 20):
                break

            if (convergence <= tollerance) or  ( diff_var == 1.e-11):
                #megaloop = False
                #loop = True
                megaloop = False
                break
            
        if not megaloop:
            break
        # print("*** initial_guess :", initial_guess)
        # print('{:4d} {:>20f} {:>20f} {:>20f} {:>20f} {:>20f} {:>20f}'.format(
        #        I, exponent[0], exponent[1] ,convergence, varandas0(exponent),
        #        GV[0], GV[1]))
        #----
        #print(' {:4d}        {:>1.11f}        {:>1.11f}        {:>1.11f}        {:>1.11f}        {:>1.20f}        {:>1.20f}'.format(
        #            I, exponent[0], exponent[1], convergence,
        #            valueV0, GV[0], GV[1]))
    # else:
    #     exponent = None

    return abs(exponent)




if __name__ == "__main__":

    ED = []
    EM_A = []
    EM_B = []
    VAR_FUN = []
    ALPHA_EXP = []
    BETA_EXP = []	
    ed = []	
    ema = []	
    emb = []	
    ej = []	
    ex = []	
    ED_UPDATE = []	
    EMA_UPDATE = []	
    EMB_UPDATE = []	
    EJ_UPDATE = []	
    EX_UPDATE = []	
    MP2_D  = []
    MP2_A  = []
    MP2_B  = []

    MP2D_UPDATE = []	
    MP2A_UPDATE = []	
    MP2B_UPDATE = []	

    #----
    print("Legend:")
    print(" -------")
    print("\n")
    print("Legend:")
    print(" -------")
    print(" 1st col.: COUNTER")
    print(" 2nd col.: ALPHA")
    print(" 3rd col.: BETA")
    print(" 4th col.: CONVERGENCE")
    print(" 5th col.: VARANDAS")
    print(" 6th col.: GRADIENT ALPHA")
    print(" 7th col.: GRADIENT BETA")
    print(" -----------------------------")

    initial_guess = np.array(['0.21954348034','0.16697708112'], dtype = np.double ) 


    exponent = downward(varandas0, gradient, abs(initial_guess), epsilon=20,
                            MaxIter=50, tollerance= 0.000005, h=0.05)                      
 
    print("-"* 125,"\n")
    for i, (x,j,p) in enumerate(zip(VAR_FUN, ALPHA_EXP, BETA_EXP)):
        
        if x == min(VAR_FUN):	
            print ("Minimum VARANDAS value found in the whole set, ALPHA, and BETA:", x, j, p,"\n")
			
    #----
    print("Legend:")
    print(" -------")
    print("# Z-Axis:  VARANDAS") 
    print("# X-Axis:  ALPHA") 
    print("# Y-Axis:  BETA") 
    print("# Z-Axis2: DIMER ENERGY (Ha)") 
    print("# X-Axis2: MONOMER_A ENERGY (Ha)") 
    print("# Y-Axis2: MONOMER_B ENERGY (Ha)") 
    print("# Y-Axis3: DIMER COULOMB TERM (Ha)") 
    print("# Y-Axis3: DIMER EXCHANGE TERM (Ha)") 
    print(" ------------------------------- ")

    for i, (x, j, p, n, m, k, f, y) in enumerate(zip(VAR_FUN, ALPHA_EXP, BETA_EXP, ED_UPDATE, EMA_UPDATE, EMB_UPDATE, EJ_UPDATE, EX_UPDATE)):
        
        print ('{:>1.16f}   {:>1.11f}   {:>1.11f}   {:>1.16f}   {:>1.11f}   {:>1.11f}   {:>1.7f}   {:>1.7f}'.format(x, j, p, n, m, k, f, y))
	
    print("-"* 140,"\n")



    #----
    print("Legend:")
    print(" -------")
    print("# Z-Axis:  MP2 DIMER (Ha)") 
    print("# Y-Axis:  MP2 MONOMER A (Ha) ") 
    print("# Y-Axis:  MP2 MONOMER B (Ha) ") 
    print("# Z-Axis2: DIMER ENERGY (Ha)") 
    print("# X-Axis2: MONOMER_A ENERGY (Ha)") 
    print("# Y-Axis2: MONOMER_B ENERGY (Ha)") 
    print("# Y-Axis3: DIMER COULOMB TERM (Ha)") 
    print("# Y-Axis3: DIMER EXCHANGE TERM (Ha)") 
    print(" ------------------------------- ")

    for i, (x, j, p, n, m, k, f, y) in enumerate(zip(MP2D_UPDATE, MP2A_UPDATE, MP2B_UPDATE, ED_UPDATE, EMA_UPDATE, EMB_UPDATE, EJ_UPDATE, EX_UPDATE)):
        
        print ('{:>1.16f}   {:>1.11f}   {:>1.11f}   {:>1.16f}   {:>1.11f}   {:>1.11f}   {:>1.7f}   {:>1.7f}'.format(x, j, p, n, m, k, f, y))
	
    print("-"* 140,"\n")
		
