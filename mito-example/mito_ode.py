import numpy as np
import pickle
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter,ScalarFormatter
import matplotlib as mpl
from matplotlib.lines import Line2D

'''
Adapted from code written by Guadalupe Garcia
'''

from scipy.integrate import odeint

def mito_ode_calc(timeSpan, geoParam, ant_rate):
    vmito = geoParam[0]*1e-15
    vims = geoParam[1]*1e-15
    vcube = geoParam[2]*1e-15

    def f(y,t,vmito,vims,vcube, ant_rate):
        Na = 6.02214e23
        # ant_rate = 1e10 # rel ANT rate

        no_ant = 16471 #number of ants
        no1_atp = 267 #number of atphase
        kp = ant_rate * 9.2
        kcp = ant_rate * 0.35
        kt = ant_rate * 0.58
        kd = ant_rate * 0.48

        # ATP synthase rates
        k12 = 100.32#100 #s-1
        k21 = 40.0 #s-1
        k23 = 5.0 #uM-1s-1
        k32 = 5e3 #s-1
        k65 = 3969             #s-1
        k56 = 1002#1e3              #s-1,
        k16 = 147956#146516           #s-1
        k61 = 33658.7#33989            #s-1
        k25 = 5.85e-30         #s-1
        k52 = 1e-20            #s-1
        k54 = 1e2              #s-1
        k45 = 1e2              #s-1
        k43 = 2                #uM-1s-1
        k34 = 1e2

        # VDAC
        kvdac  = 1.0
        n_porin = 6340

        # binding and unbinding to ANTs
        kDi_on = 4.0 #          #KDo
        kDi_off = 100 #        #KDo
        kTm_on = 6.4#         #KTi
        kTm_off = 4e4#         #KTi
        kTi_on = 0.4           #KTo
        kTi_off = 200.0        #KTo
        kDm_on = 4.0           #KDi
        kDm_off = 4e4          #KDi
        fa1 = 0.5
        fa2 = 1.0

        D_matrix, T_matrix, T_cyto, T_IMS, L, DL, LD, LT, TL, TLD, DLD, DLDp, TLT, TLTp, E_IMS, E_Mat, E_Mat_H3Star, E_Mat_H3S, E_Mat_H3, D_IMS = y
        ac = D_IMS#nr_do_i
        dydt = [(-kDm_on*(1e6*D_matrix/(Na*vmito))*TL + kDm_off*TLD - kDm_on*(1e6*D_matrix/(Na*vmito))*L + 
                 kDm_off*LD - 2.0*kDm_on*fa1*(1e6*D_matrix/(Na*vmito))*DL + kDm_off*DLD + kDm_off*DLDp + 
                 k34*E_Mat_H3S - k43*(1e6*D_matrix/(Na*vmito))*E_Mat_H3Star),

                (-kTm_on*(1e6*T_matrix/(Na*vmito))*DL + kTm_off*(no_ant-L-DL-LD-LT-TL-TLD-DLD-TLT-DLDp-TLTp) - 
                 kTm_on*(1e6*T_matrix/(Na*vmito))*L + kTm_off*LT - 2.0*kTm_on*fa1*(1e6*T_matrix/(Na*vmito))*TL + 
                 kTm_off*TLT + kTm_off*TLTp -k23*E_Mat_H3*(1e6*T_matrix/(Na*vmito)) + k32*E_Mat_H3S),

                (-kvdac*(1e6*T_cyto/(Na*vcube))*n_porin + kvdac*(1e6*T_IMS/(Na*vims))*n_porin),

                (-kTi_on*(1e6*T_IMS/(Na*vims))*L + kTi_off*TL - kTi_on*(1e6*T_IMS/(Na*vims))*LD + kTi_off*TLD - 
                 2.0*kTi_on*fa1*(1e6*T_IMS/(Na*vims))*LT + kTi_off*TLT + kTi_off*TLTp - 
                 kvdac*(1e6*T_IMS/(Na*vims))*n_porin + kvdac*(1e6*T_cyto/(Na*vcube))*n_porin),

                (-kDi_on*(1e6*ac/(Na*vims))*L + kDi_off*DL- kTm_on*(1e6*T_matrix/(Na*vmito))*L + kTm_off*LT - 
                 kTi_on*(1e6*T_IMS/(Na*vims))*L + kTi_off*TL - kDm_on*(1e6*D_matrix/(Na*vmito))*L + kDm_off*LD),

                (kDi_on*(1e6*ac/(Na*vims))*L - kDi_off*DL - kTm_on*(1e6*T_matrix/(Na*vmito))*DL + 
                 kTm_off*(no_ant-L-DL-LD-LT-TL-TLD-DLD-TLT-DLDp-TLTp) - kDm_on*fa1*(1e6*D_matrix/(Na*vmito))*DL + 
                 kDm_off*DLD - kDm_on*fa1*(1e6*D_matrix/(Na*vmito))*DL + kDm_off*DLDp),

                (kDm_on*(1e6*D_matrix/(Na*vmito))*L - kDm_off*LD  - kTi_on*(1e6*T_IMS/(Na*vims))*LD + kTi_off*TLD - 
                 kDi_on*fa1*(1e6*ac/(Na*vims))*LD + kDi_off*DLD - kDi_on*fa1*(1e6*ac/(Na*vims))*LD + kDi_off*DLDp),

                (kTm_on*(1e6*T_matrix/(Na*vmito))*L - kTm_off*LT - kDi_on*(1e6*ac/(Na*vims))*LT + 
                 kDi_off*(no_ant-L-DL-LD-LT-TL-TLD-DLD-TLT-DLDp-TLTp) - kTi_on*fa1*(1e6*T_IMS/(Na*vims))*LT + 
                 kTi_off*TLT - kTi_on*fa1*(1e6*T_IMS/(Na*vims))*LT + kTi_off*TLTp),

                (kTi_on*(1e6*T_IMS/(Na*vims))*L - kTi_off*TL- kDm_on*(1e6*D_matrix/(Na*vmito))*TL + kDm_off*TLD - 
                 kTm_on*fa1*(1e6*T_matrix/(Na*vmito))*TL + kTm_off*TLT- kTm_on*fa1*(1e6*T_matrix/(Na*vmito))*TL + kTm_off*TLTp),

                (kDm_on*(1e6*D_matrix/(Na*vmito))*TL - kDm_off*TLD - kcp*TLD + kp*(no_ant-L-DL-LD-LT-TL-TLD-DLD-TLT-DLDp-TLTp) + 
                 kTi_on*(1e6*T_IMS/(Na*vims))*LD - kTi_off*TLD),

                (kDi_on*fa2*(1e6*ac/(Na*vims))*LD - kDi_off*DLD + kDm_on*fa2*(1e6*D_matrix/(Na*vmito))*DL - 
                 kDm_off*DLD),# -kd*DLD + kd*DLDp),

                0*(kDi_on*fa2*(1e6*ac/(Na*vims))*LD - kDi_off*DLDp + kDm_on*fa2*(1e6*D_matrix/(Na*vmito))*DL - 
                 kDm_off*DLDp - kd*DLDp + kd*DLD),

                (kTm_on*fa2*(1e6*T_matrix/(Na*vmito))*TL - kTm_off*TLT + kTi_on*fa2*(1e6*T_IMS/(Na*vims))*LT - 
                 kTi_off*TLT),# + kt*TLTp - kt*TLT),

                0*(kTm_on*fa2*(1e6*T_matrix/(Na*vmito))*TL - kTm_off*TLTp + kTi_on*fa2*(1e6*T_IMS/(Na*vims))*LT - 
                 kTi_off*TLTp -kt*TLTp +kt*TLT),

                (-k65*E_IMS + k56*(no1_atp-E_IMS-E_Mat-E_Mat_H3Star-E_Mat_H3S-E_Mat_H3)+ k16*E_Mat - k61*E_IMS),

                (-k16*E_Mat + k61*E_IMS - k12*E_Mat + k21*E_Mat_H3),

                (-k45*E_Mat_H3Star + k54*(no1_atp-E_IMS-E_Mat-E_Mat_H3Star-E_Mat_H3S-E_Mat_H3) + k34*E_Mat_H3S - k43*E_Mat_H3Star*(1e6*D_matrix/(Na*vmito))),

                (k43*(1e6*D_matrix/(Na*vmito))*E_Mat_H3Star - k34*E_Mat_H3S + k23*E_Mat_H3*(1e6*(T_matrix)/(Na*vmito)) - k32*E_Mat_H3S),

                (-k23*E_Mat_H3*(1e6*(T_matrix)/(Na*vmito)) + k32*E_Mat_H3S - k25*E_Mat_H3 + k52*(no1_atp-E_IMS-E_Mat-E_Mat_H3Star-E_Mat_H3S-E_Mat_H3) + 
                 k12*E_Mat - k21*E_Mat_H3),

                0*(-kDi_on*L*(1e6*D_IMS/(Na*vims)) + kDi_off*DL - kDi_on*LT*(1e6*D_IMS/(Na*vims)) + 
                kDi_off*(no_ant-L-DL-LD-LT-TL-TLD-DLD-TLT-DLDp-TLTp)
                - 2*fa1*kDi_on*LD*(1e6*D_IMS/(Na*vims)) + kDi_off*DLD + kDi_off*DLDp)
        ]
        return dydt

    cdm_i = 2 #10mM
    ctm_i = 13 #mM
    cdo_i = 0.1 #10mM
    cto_i = 6.5 #mM
    nr_dm_m = 0.45*0.8*6.02*cdm_i*vmito*1e20 #number of molecules Dm
    nr_tm_m = 0.05*6.02*ctm_i*vmito*1e20
    nr_do_i = 0.45*6.02*cdo_i*vims*1e20
    nr_to_i = 0.05*6.02*cto_i*vims*1e20
    nr_to_c = 0.05*6.02*cto_i*vcube*1e20
    no_ant = 16471 #number of ants
    no1_atp = 267 #number of atphase
    # initial conditions
    # order: D_matrix, T_matrix, T_cyto, T_IMS, L, DL, LD, LT, TL, TLD, DLD, DLDp, TLT, TLTp, E_IMS, E_Mat, E_Mat_H3Star, E_Mat_H3S, E_Mat_H3, D_IMS
    ant_frac = 1.0
    rem = (1-ant_frac)/10
    yinit  = [nr_dm_m, nr_tm_m, nr_to_c, nr_to_i, ant_frac*no_ant, 
              rem*no_ant, rem*no_ant, rem*no_ant, rem*no_ant, rem*no_ant, rem*no_ant, 
              rem*no_ant, rem*no_ant, rem*no_ant, no1_atp, 0, 0, 0, 0, nr_do_i]

    t = np.linspace(timeSpan[0], timeSpan[1], 1000)
    sol = odeint(f, yinit, t, args=(vmito,vims,vcube,ant_rate))

    return t, sol