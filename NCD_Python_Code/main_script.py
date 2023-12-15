# Python code to model non-cyanobacterial diazotrophs
# associated with sinking marine particles

# File: main_script.py

import numpy as np
from parameters import *
from scipy.io import savemat

# Time range
t_start = 0
t_end = 16
t_div = 2e-06 
t_length = int((t_end-t_start)/t_div)+1
t_range = np.linspace(t_start, t_end, t_length)
tspan = (t_start, t_end)

# psi range
psi_start = 0
psi_end = 1
length_psi_range = 50 
psi_range = np.linspace(psi_start,psi_end,length_psi_range) # cm
psi_div = (psi_end-psi_start)/(length_psi_range-1)

# Temperature 
temp_current = 17

# Division inside a single particle 
xstart = 1e-06 # cm
x_length = 50 

# Calculating reference carbon content and radius of a particle
xend = init_rad # cm % chosing a large partile to calculate reference carbon content
x = np.linspace(xstart,xend,x_length) # cm
dx = (xend-xstart)/(x_length-1) # cm

C_labile_conc = (f_G_C*C_init + f_A_C*P_init)*np.ones(x_length) # mug C/L
C_tot_labile_ref = 4*pi*np.sum(np.power(x, 2)*C_labile_conc)*dx*10**(-3) # *0.2*0.005; % mug C/particle;

C_nonlabile_conc = (f_G_C*((1-f_PS)/f_PS)*C_init + f_A_C*((1-f_PP)/f_PP)*P_init)*np.ones(x_length) # mug C/L
C_tot_nonlabile_ref = 4*pi*np.sum(np.power(x, 2)*C_nonlabile_conc)*dx*10**(-3) # *0.2*0.005; % mug C/particle;

C_ref = C_tot_labile_ref+C_tot_nonlabile_ref # mug C/particle
rp_ref = xend # fixed reference size

x_init = x
dx_init = dx

# Creating empty array for saving variables
y = np.zeros((7*x_length))
J_N2_per_day_1 = np.zeros(t_length+1)
J_N2_per_day_1[0] = 0

# Initial concentrations
y = np.concatenate((C_init*np.ones(x_length), P_init*np.ones(x_length), G_init*np.ones(x_length), A_init*np.ones(x_length), B_init*np.ones(x_length), O2_init*np.ones(x_length), NO3_init*np.ones(x_length)), axis=None)

# Calculating viscocity at the current temperature by interpolation method
vis_all = interpolate.interp1d(vis_temp_data, vis_data, fill_value = "extrapolate")
vis = vis_all(temp_current) 

# modification of parameters according to the Q10 rule
# hydrolysis
h_C_Q10 = h_C*Q10_h**((temp_current-ref_temp_hC)/10)
h_P_Q10 = h_P*Q10_h**((temp_current-ref_temp_hP)/10)
A_C_Q10 = A_C*Q10_h**((temp_current-ref_temp_AC)/10)
A_P_Q10 = A_P*Q10_h**((temp_current-ref_temp_AP)/10)

# uptake
J_max_G_Q10 = J_max_G*Q10_u**((temp_current-ref_temp_JmaxG)/10)
J_max_A_Q10 = J_max_A*Q10_u**((temp_current-ref_temp_JmaxA)/10)
alpha_G_Q10 = alpha_G*Q10_u**((temp_current-ref_temp_alphaG)/10)
alpha_A_Q10 = alpha_A*Q10_u**((temp_current-ref_temp_alphaA)/10)
   
J_max_NO3_Q10 = J_max_NO3*Q10_u**((temp_current-ref_temp_JmaxNO3)/10)
alpha_NO3_Q10 = alpha_NO3*Q10_u**((temp_current-ref_temp_alphaNO3)/10)
    
J_max_SO4_Q10 = J_max_SO4*Q10_u**((temp_current-ref_temp_JmaxSO4)/10)
    
M_N2_Q10 = M_N2*Q10_u**((temp_current-ref_temp_MN2)/10)
    
# respiration
R_B_Q10 = R_B*Q10_r**((temp_current-ref_temp_RB)/10)
R_G_Q10 = R_G*Q10_r**((temp_current-ref_temp_RG)/10)
R_A_Q10 = R_A*Q10_r**((temp_current-ref_temp_RA)/10)
R_N2_Q10 = R_N2*Q10_r**((temp_current-ref_temp_RN2)/10)
R_NO3_Q10 = R_NO3*Q10_r**((temp_current-ref_temp_RNO3)/10)
R_SO4_Q10 = R_SO4*Q10_r**((temp_current-ref_temp_RSO4)/10)
    
# diffusion
D_O2_cell_Q10 = D_O2_cell*ref_vis_DO2*(K+temp_current)/(vis*(K+ref_temp_DO2))
D_G_Q10 = D_G*ref_vis_DM*(K+temp_current)/(vis*(K+ref_temp_DM))
D_A_Q10 = D_A*ref_vis_DM*(K+temp_current)/(vis*(K+ref_temp_DM))
D_OP_Q10 = D_OP*ref_vis_DO2*(K+temp_current)/(vis*(K+ref_temp_DO2))
D_NP_Q10 = D_NP*ref_vis_DNO3*(K+temp_current)/(vis*(K+ref_temp_DNO3))

# memory allocation for future saving
J_N2_per_day_1 = np.zeros((t_length+1))
xend_1 = np.zeros((t_length+1))
    
J_N2_1 = np.zeros((t_length+1,x_length))
R_bar_1 = np.zeros((t_length+1,x_length))
mu_1 = np.zeros((t_length+1,x_length))
J_C_1 = np.zeros((t_length+1,x_length))
J_N_1 = np.zeros((t_length+1,x_length))
EA_1 = np.zeros((t_length+1,x_length))
R_B_1 = np.zeros((t_length+1,x_length))
R_E_1 = np.zeros((t_length+1,x_length))
R_N2_1 = np.zeros((t_length+1,x_length))
R_O2_1 = np.zeros((t_length+1,x_length))
J_C_L_1 = np.zeros((t_length+1,x_length))
J_P_L_1 = np.zeros((t_length+1,x_length))
J_G_1 = np.zeros((t_length+1,x_length))
J_A_1 = np.zeros((t_length+1,x_length))

DGDx_1 = np.zeros((t_length+1,x_length))
DADx_1 = np.zeros((t_length+1,x_length))
DODx_1 = np.zeros((t_length+1,x_length))
DNDx_1 = np.zeros((t_length+1,x_length))
V_O2_pot_1 = np.zeros((t_length+1,x_length))
V_O2_1 = np.zeros((t_length+1,x_length))
V_NO3_1 = np.zeros((t_length+1,x_length))
V_SO4_1 = np.zeros((t_length+1,x_length))

y_1 = np.zeros((t_length+1,7*x_length))
y_1[0,:] = np.concatenate((C_init*np.ones(x_length), P_init*np.ones(x_length), G_init*np.ones(x_length), A_init*np.ones(x_length), B_init*np.ones(x_length), O2_init*np.ones(x_length), NO3_init*np.ones(x_length)), axis=None)

# creating new time array considering time after 1000 time steps
dist = 1000
t_range_save = len(t_range[0:-1:dist])+1

t_range_1_1 = np.zeros(t_range_save)
xend_1_1 = np.zeros(t_range_save)

y_1_1 = np.zeros((t_range_save,x_length))

J_N2_1_1 = np.zeros((t_range_save,x_length))
R_bar_1_1 = np.zeros((t_range_save,x_length))
mu_1_1 = np.zeros((t_range_save,x_length))
J_C_1_1 = np.zeros((t_range_save,x_length))
J_N_1_1 = np.zeros((t_range_save,x_length))
EA_1_1 = np.zeros((t_range_save,x_length))
R_B_1_1 = np.zeros((t_range_save,x_length))
R_E_1_1 = np.zeros((t_range_save,x_length))
R_N2_1_1 = np.zeros((t_range_save,x_length))
R_O2_1_1 = np.zeros((t_range_save,x_length))
J_C_L_1_1 = np.zeros((t_range_save,x_length))
J_P_L_1_1 = np.zeros((t_range_save,x_length))
J_G_1_1 = np.zeros((t_range_save,x_length))
J_A_1_1 = np.zeros((t_range_save,x_length))

DGDx_1_1 = np.zeros((t_range_save,x_length))
DADx_1_1 = np.zeros((t_range_save,x_length))
DODx_1_1 = np.zeros((t_range_save,x_length))
DNDx_1_1 = np.zeros((t_range_save,x_length))
V_O2_pot_1_1 = np.zeros((t_range_save,x_length))
V_O2_1_1 = np.zeros((t_range_save,x_length))
V_NO3_1_1 = np.zeros((t_range_save,x_length))
V_SO4_1_1 = np.zeros((t_range_save,x_length))


#------------------------------------------------------------------
#------------------------- Main function --------------------------
#------------------------------------------------------------------

def nfix(y,psi_range,length_psi_range,x_length,xstart,C_tot_nonlabile_ref,C_ref,rp_ref, 
         h_C_Q10,h_P_Q10,A_C_Q10,A_P_Q10,J_max_G_Q10,J_max_A_Q10,alpha_G_Q10,alpha_A_Q10,J_max_NO3_Q10,
         alpha_NO3_Q10,J_max_SO4_Q10,M_N2_Q10,R_B_Q10,R_G_Q10,R_A_Q10,R_N2_Q10,R_NO3_Q10,R_SO4_Q10,
         D_O2_cell_Q10,D_G_Q10,D_A_Q10,D_OP_Q10,D_NP_Q10,f_G_C,f_A_C,f_A_N,alpha,x_B,Q_CN_B,
         rho_CO,rho_CNO3,rho_CSO4,m_B,G_bdry,A_bdry,O2_bdry,NO3_bdry,x,dx): 
    
    # -------- memory allocation ---------
    
    DyDt =  np.empty(7*x_length) #np.zeros(7*x_length)

    DCDt = np.empty(x_length)
    DPDt = np.empty(x_length)
    DGDt = np.empty(x_length)
    DADt = np.empty(x_length)
    DBDt = np.empty(x_length)
    DODt = np.empty(x_length)
    DNDt = np.empty(x_length)
    
    CF = np.empty(x_length) 
    PF = np.empty(x_length) 
    GF = np.empty(x_length) 
    AF = np.empty(x_length) 
    BF = np.empty(x_length) 
    OF = np.empty(x_length) 
    NF = np.empty(x_length) 
    
    DGDx = np.empty(x_length)
    DADx = np.empty(x_length)
    DODx = np.empty(x_length)
    DNDx = np.empty(x_length)

    psi = np.empty(x_length)
    
    C = y[0*x_length:1*x_length]
    P = y[1*x_length:2*x_length]
    G = y[2*x_length:3*x_length]
    A = y[3*x_length:4*x_length]
    B = y[4*x_length:5*x_length]
    O2 = y[5*x_length:6*x_length]
    NO3 = y[6*x_length:7*x_length]
    
    # -------- division of particle ---------
    
    C_labile = f_G_C*C+ f_A_C*P # mug C/L
    C_tot = C_tot_nonlabile_ref + 4*np.pi*np.sum(np.power(x,2)*C_labile)*dx*10**(-3) # *0.2*0.005; % mug C/particle;
    xend = np.power(C_tot/C_ref,1/alpha)*rp_ref
    x = np.linspace(xstart,xend,x_length) # cm % division of particle
    dx = (xend-xstart)/(x_length-1) # cm
    
    # -------- Cellular processes --------- 
    # ------ Hydrolysis of polymer ------
    
    J_C_L = (h_C_Q10*A_C_Q10*C)/(h_C_Q10+A_C_Q10*C) # mug G/cell/d
    J_P_L = (h_P_Q10*A_P_Q10*P)/(h_P_Q10+A_P_Q10*P) # mug A/cell/d
    
    # ----- Uptake of monomer ------
    
    J_G = (J_max_G_Q10*alpha_G_Q10*G)/(J_max_G_Q10+alpha_G_Q10*G) # mug G/cell/d
    J_A = (J_max_A_Q10*alpha_A_Q10*A)/(J_max_A_Q10+alpha_A_Q10*A) # mug A/cell/d
    
    # ----- Amount of C and N ------
    
    J_DOC = f_G_C*J_G+f_A_C*J_A # mug C/cell/d 
    J_DON =     0*J_G+f_A_N*J_A # mug N/cell/d 
    
    #------ Uptake of NO3 ------
    
    J_NO3 = (J_max_NO3_Q10*alpha_NO3_Q10*NO3)/(J_max_NO3_Q10+alpha_NO3_Q10*NO3) # mug N/cell/d

    # ------- Optimal psi ---------
    
    for j in range(0, x_length, 1):
        psi[j] = opt_psi_func(O2[j],length_psi_range,
              J_G[j],J_A[j],J_NO3[j],J_DOC[j],J_DON[j],x_B,Q_CN_B,rho_CO,rho_CNO3,
              rho_CSO4,J_max_SO4_Q10,M_N2_Q10,R_B_Q10,R_G_Q10,R_A_Q10,R_N2_Q10,R_NO3_Q10,R_SO4_Q10,D_O2_cell_Q10)

    # ------- N fixation rate ---------
    
    J_N2 = psi*M_N2_Q10 # mug N/cell/d
    J_N2_per_day = 0.01*4*np.pi*np.sum(np.power(x, 2)*B*J_N2)*dx*10**(-3) # mug N/particle/day
    
    # ------- Total incoming C and N into a cell ---------
    
    J_C = J_DOC # mug C/cell/d
    J_N = J_DON+J_N2 # mug N/cell/d
    
    # ------- Respiration rate ---------
    
    R1 = R_B_Q10*x_B + R_G_Q10*J_G + R_A_Q10*J_A + R_N2_Q10*Q_CN_B*J_N2 # mug C/cell/d  
    R_B = R_B_Q10*x_B
    R_E = 0.6*Q10_r**((temp_current-ref_temp_RB)/10)*x_B
    R_N2 = R_N2_Q10*Q_CN_B*J_N2
    
    # ------- Incoming O2 flux ---------
    
    O2in = np.maximum(0,O2 - R1/(rho_CO*D_O2_cell_Q10)) # mumol O2/L
    
    # ------- Cost of O2 removal ---------
    
    R_O2 = rho_CO*D_O2_cell_Q10*O2in # mug C/cell/d
    R_O2[psi==0] = 0 # effective only when N fixation is happening
            
    # ------- Sum of costs of basal respiration and O2 removal ---------
    
    R = R1 + R_O2 # mug C/cell/d
    
    # ------- Incoming flux of O2, NO3, SO4 ---------
    
    F_O2_max = rho_CO*D_O2_cell_Q10*O2 #  mug C/d/cell % maximum rate at which O2 can diffuse into the cell (no O2 in cell, all are used)
    F_O2 = np.maximum(0,rho_CO*D_O2_cell_Q10*(O2-O2in)) # mug C/d/cell % actual O2 diffusion rate into the cell (Cost of O2 removal)
    F_NO3_max = rho_CNO3*J_NO3 #  mug C/d/cell % maximum rate at which NO3 can be taken by the cell
    F_NO3 = np.minimum(F_NO3_max,np.maximum(0,(R-F_O2_max))) # mug C/d/cell % actual rate at which NO3 can be taken by the cell
    F_SO4_max = rho_CSO4*J_max_SO4_Q10 # mug C/cell/d % maximum SO4 uptake rate by cell
    F_SO4 = np.minimum(F_SO4_max,np.maximum(0,(R-(F_O2_max+F_NO3_max)))) # mug C/cell/d % actual amount of sulfate uptake
    
    # ------- Uptake cost of NO3 and SO4 ---------
    
    R_N_uptk = R_NO3_Q10*F_NO3 # mug C/d/cell % Cost of nitrate uptake
    R_S_uptk = R_SO4_Q10*F_SO4 # mug C/cell/d % Cost of sulfate uptake
    
    # ------- Total respiration rate ---------
    
    R_bar = R + R_N_uptk + R_S_uptk # mug C/cell/d
    
    #------- Available electron acceptor for growth ---------
    
    EA = F_O2+F_NO3+F_SO4 # mug C/cell/d
    
    # ------- Alailable fluxes for growth and growth rate ---------
    
    J_tot = np.array([J_C-R_bar, Q_CN_B*J_N, EA]) # mug C/cell/d
    mu = J_tot.min(axis=0)/x_B # 1/d % division rate
    
    #------- Used O2 and NO3 to calculate remaining concentrations inside particle ---------
    
    V_O2_pot = np.maximum(0,F_O2_max)/rho_CO #  mumol O2/cell/d
    V_O2 = np.maximum(0,F_O2)/rho_CO #  mumol O2/d/cell
    V_NO3 = np.maximum(0,F_NO3)/rho_CNO3 # mug NO3/d/cell
    V_SO4 = np.maximum(0,F_SO4)/rho_CSO4 # mug SO4/d/cell
    
    # ----------- Main equations ------------
    # ________ Changes due to reaction part ________
    
    CF = -J_C_L*B # mug C/d/L
    PF = -J_P_L*B # mug C/d/L
    GF = J_C_L*B-J_G*B # mug C/d/L
    AF = J_P_L*B-J_A*B # mug C/d/L
    BF = mu*B-m_B*B # mug C/d/L
    OF = -V_O2*B # mumol O2/d/L
    NF = -V_NO3*B # mug NO3/d/L
    
    # ________ Changes due to diffusion part ________
    # ---- Left boundary (center of particle) ----
    
    DGDx[0] = D_G_Q10*(1/dx**2)*(2*G[1] - 2*G[0]) #  du/dr = 0;
    DADx[0] = D_A_Q10*(1/dx**2)*(2*A[1] - 2*A[0])  
    DODx[0] = D_OP_Q10*(1/dx**2)*(2*O2[1] - 2*O2[0])  
    DNDx[0] = D_NP_Q10*(1/dx**2)*(2*NO3[1] - 2*NO3[0])  
    
    DCDt[0] = CF[0]
    DPDt[0] = PF[0]
    DGDt[0] = GF[0] + DGDx[0]
    DADt[0] = AF[0] + DADx[0]
    DBDt[0] = BF[0]
    DODt[0] = OF[0] + DODx[0]
    DNDt[0] = NF[0] + DNDx[0]
    
    # ---- Interior ----
    
    DGDx[1:-1] = (D_G_Q10/dx**2)*(G[:-2] - 2*G[1:-1] + G[2:]) + D_G_Q10*(2/x[1:-1])*(1/(2*dx))*(G[2:] - G[:-2])
    DADx[1:-1] = (D_A_Q10/dx**2)*(A[:-2] - 2*A[1:-1] + A[2:]) + D_A_Q10*(2/x[1:-1])*(1/(2*dx))*(A[2:] - A[:-2])
    DODx[1:-1] = (D_OP_Q10/dx**2)*(O2[:-2] - 2*O2[1:-1] + O2[2:]) + D_OP_Q10*(2/x[1:-1])*(1/(2*dx))*(O2[2:] - O2[:-2])
    DNDx[1:-1] = (D_NP_Q10/dx**2)*(NO3[:-2] - 2*NO3[1:-1] + NO3[2:]) + D_NP_Q10*(2/x[1:-1])*(1/(2*dx))*(NO3[2:] - NO3[:-2])

    
    DCDt[1:-1] = CF[1:-1]
    DPDt[1:-1] = PF[1:-1]
    DGDt[1:-1] = GF[1:-1] + DGDx[1:-1]
    DADt[1:-1] = AF[1:-1] + DADx[1:-1]
    DBDt[1:-1] = BF[1:-1]
    DODt[1:-1] = OF[1:-1] + DODx[1:-1]
    DNDt[1:-1] = NF[1:-1] + DNDx[1:-1]
    
    # ---- Right boundary (surface of particle) ----
    
    DGDx[-1] = D_G_Q10*(1/dx**2)*(G[-2] - 2*G[-1] + G_bdry) + D_G_Q10*(2/x[-1])*(1/(2*dx))*(G_bdry - G[-2])
    DADx[-1] = D_A_Q10*(1/dx**2)*(A[-2] - 2*A[-1] + A_bdry) + D_A_Q10*(2/x[-1])*(1/(2*dx))*(A_bdry - A[-2])
    DODx[-1] = D_OP_Q10*(1/dx**2)*(O2[-2] - 2*O2[-1] + O2_bdry) + D_OP_Q10*(2/x[-1])*(1/(2*dx))*(O2_bdry - O2[-2])
    DNDx[-1] = D_NP_Q10*(1/dx**2)*(NO3[-2] - 2*NO3[-1] + NO3_bdry) + D_NP_Q10*(2/x[-1])*(1/(2*dx))*(NO3_bdry - NO3[-2])
    
    DCDt[-1] = CF[-1]
    DPDt[-1] = PF[-1]
    DGDt[-1] = GF[-1] + DGDx[-1]
    DADt[-1] = AF[-1] + DADx[-1]
    DBDt[-1] = BF[-1]
    DODt[-1] = OF[-1] + DODx[-1]
    DNDt[-1] = NF[-1] + DNDx[-1]
    
    DyDt = np.concatenate((DCDt, DPDt, DGDt, DADt, DBDt, DODt, DNDt), axis=None)

    return (y,DyDt,x,dx,J_N2_per_day,J_N2,R_bar,mu,xend,J_C,J_N,EA,R_B,R_E,R_N2,R_O2,DGDx,DADx,DODx,DNDx,V_O2_pot,V_O2,V_NO3,V_SO4,J_C_L,J_P_L,J_G,J_A)





#------------------------------------------------------------------
#------------------- Function for optimal psi ---------------------
#------------------------------------------------------------------


def opt_psi_func(OO2,length_psi_range,
    JJ_G,JJ_A,JJ_NO3,JJ_DOC,JJ_DON,x_B,Q_CN_B,rho_CO,rho_CNO3,
    rho_CSO4,J_max_SO4_Q10,M_N2_Q10,R_B_Q10,R_G_Q10,R_A_Q10,R_N2_Q10,
    R_NO3_Q10,R_SO4_Q10,D_O2_cell_Q10):
    # -------------------------------------------
    psi_range = np.linspace(0,1,length_psi_range)
    JJ_N2 = psi_range*M_N2_Q10
    # -------------------------------------------
    JJ_C = JJ_DOC
    JJ_N = JJ_DON+JJ_N2
    # -------------------------------------------
    RR1 = R_B_Q10*x_B + R_G_Q10*JJ_G + R_A_Q10*JJ_A + R_N2_Q10*Q_CN_B*JJ_N2
    # -------------------------------------------
    OO2in = np.maximum(0,OO2-RR1/(rho_CO*D_O2_cell_Q10))
    RR_O2 = rho_CO*D_O2_cell_Q10*OO2in # mug C/cell/d
    RR_O2[psi_range==0] = 0
    RR = RR1 + RR_O2
    # -------------------------------------------
    FF_O2_max = rho_CO*D_O2_cell_Q10*OO2 # mugC/cell/d % maximum O2 diffuse rate into the cell
    FF_O2 = np.maximum(0,rho_CO*D_O2_cell_Q10*(OO2-OO2in)) # actual O2 diffuse rate into the cell
    FF_NO3_max = rho_CNO3*JJ_NO3 # maximum NO3 uptake rate by cell
    FF_NO3 = np.minimum(FF_NO3_max,np.maximum(0,(RR-FF_O2_max))) # actual amount of nitrate uptake
    FF_SO4_max = rho_CSO4*J_max_SO4_Q10 # maximum SO4 uptake rate by cell
    FF_SO4 = np.minimum(FF_SO4_max,np.maximum(0,(RR-(FF_O2_max+FF_NO3_max)))) # actual amount of sulfate uptake
    # -------------------------------------------
    RR_N_uptk = R_NO3_Q10*FF_NO3 # Cost of nitrate uptake
    RR_S_uptk = R_SO4_Q10*FF_SO4 # Cost of sulfate uptake
    RR_bar = RR + RR_N_uptk + RR_S_uptk
    # -------------------------------------------
    EEA = FF_O2+FF_NO3+FF_SO4
    # -------------------------------------------
    JJ_C_1 = JJ_C-RR_bar
    JJ_N_1 = Q_CN_B*JJ_N
    mu1_arr = np.array([JJ_C_1, JJ_N_1,EEA])
    mu1 = mu1_arr.min(axis=0)/x_B # just need min value, so no need to divide by mass
    # -------------------------------------------
    mu1[(JJ_N_1>JJ_C_1) | (JJ_N_1>EEA)] = -10**10;
    # -------------------------------------------
    idx = np.argmax(mu1)
    psi2 = psi_range[idx]
    return psi2





  
# changing time --------------------------

for t in range(0, t_length):
    y,DyDt,x,dx,J_N2_per_day,J_N2,R_bar,mu,xend,J_C,J_N,EA,R_B,R_E,R_N2,R_O2,DGDx,DADx,DODx,DNDx,V_O2_pot,V_O2,V_NO3,V_SO4,J_C_L,J_P_L,J_G,J_A = nfix(y,psi_range,length_psi_range,x_length,xstart,C_tot_nonlabile_ref,C_ref,rp_ref, 
                                    h_C_Q10,h_P_Q10,A_C_Q10,A_P_Q10,J_max_G_Q10,J_max_A_Q10,alpha_G_Q10,alpha_A_Q10,J_max_NO3_Q10,
                                    alpha_NO3_Q10,J_max_SO4_Q10,M_N2_Q10,R_B_Q10,R_G_Q10,R_A_Q10,R_N2_Q10,R_NO3_Q10,R_SO4_Q10,
                                    D_O2_cell_Q10,D_G_Q10,D_A_Q10,D_OP_Q10,D_NP_Q10,f_G_C,f_A_C,f_A_N,alpha,x_B,Q_CN_B,rho_CO,
                                    rho_CNO3,rho_CSO4,m_B,G_bdry,A_bdry,O2_bdry,NO3_bdry,x,dx)
    y = y + t_div*DyDt
    y[y < 0] = 0    

    # Saving variables and rates at time t+1     
    J_N2_per_day_1[t+1] = J_N2_per_day # mug N/particle/day
    
    xend_1[t+1] = xend
    
    J_N2_1[t+1,:] = J_N2
    R_bar_1[t+1,:] = R_bar
    mu_1[t+1,:] = mu
    J_C_1[t+1,:] = J_C
    J_N_1[t+1,:] = J_N
    EA_1[t+1,:] = EA
    R_B_1[t+1,:] = R_B
    R_E_1[t+1,:] = R_E
    R_N2_1[t+1,:] = R_N2
    R_O2_1[t+1,:] = R_O2
    J_C_L_1[t+1,:] = J_C_L
    J_P_L_1[t+1,:] = J_P_L
    J_G_1[t+1,:] = J_G
    J_A_1[t+1,:] = J_A
    
    DGDx_1[t+1,:] = DGDx
    DADx_1[t+1,:] = DADx
    DODx_1[t+1,:] = DODx
    DNDx_1[t+1,:] = DNDx
    V_O2_pot_1[t+1,:] = V_O2_pot
    V_O2_1[t+1,:] = V_O2
    V_NO3_1[t+1,:] = V_NO3
    V_SO4_1[t+1,:] = V_SO4
    
    y_1[t+1,:] = y
        
    
    
    
            
# Calculating N2 fixation per particle

J_N2_tot = t_div*np.sum(J_N2_per_day_1) # mug N/particle

print(J_N2_tot)

# To reduce the file size while saving, 
# we consider only values at each 1000 time points

t_range_1_1 = t_range[0::dist]
xend_1_1 = xend_1[0::dist]

J_N2_1_1 = J_N2_1[0::dist,:]
R_bar_1_1 = R_bar_1[0::dist,:]
mu_1_1 = mu_1[0::dist,:]
J_C_1_1 = J_C_1[0::dist,:]
J_N_1_1 = J_N_1[0::dist,:]
EA_1_1 = EA_1[0::dist,:]
R_B_1_1 = R_B_1[0::dist,:]
R_E_1_1 = R_E_1[0::dist,:]
R_N2_1_1 = R_N2_1[0::dist,:]
R_O2_1_1 = R_O2_1[0::dist,:]
J_C_L_1_1 = J_C_L_1[0::dist,:]
J_P_L_1_1 = J_P_L_1[0::dist,:]
J_G_1_1 = J_G_1[0::dist,:]
J_A_1_1 = J_A_1[0::dist,:]

DGDx_1_1 = DGDx_1[0::dist,:]
DADx_1_1 = DADx_1[0::dist,:]
DODx_1_1 = DODx_1[0::dist,:]
DNDx_1_1 = DNDx_1[0::dist,:]
V_O2_pot_1_1 = V_O2_pot_1[0::dist,:]
V_O2_1_1 = V_O2_1[0::dist,:]
V_NO3_1_1 = V_NO3_1[0::dist,:]
V_SO4_1_1 = V_SO4_1[0::dist,:]

y_1_1 = y_1[0::dist,:]

yC = y_1_1[:,0*x_length:1*x_length]
yP = y_1_1[:,1*x_length:2*x_length]
yG = y_1_1[:,2*x_length:3*x_length]
yA = y_1_1[:,3*x_length:4*x_length]
yB = y_1_1[:,4*x_length:5*x_length]
yO2 = y_1_1[:,5*x_length:6*x_length]
yNO3 = y_1_1[:,6*x_length:7*x_length]



# # ------- save in .mat format -------

time_dat = {"time_dat": t_range_1_1}
savemat('time_dat_at_temp_'+str(Temp)+'.mat', time_dat)

rad_dat = {"rad_dat": xend_1_1}
savemat('rad_dat_at_temp_'+str(Temp)+'.mat', rad_dat)

N_fix_dat = {"N_fix_dat": J_N2_1_1}
savemat('N_fix_dat_at_temp_'+str(Temp)+'.mat', N_fix_dat)

resp_dat = {"resp_dat": R_bar_1_1}
savemat('resp_dat_at_temp_'+str(Temp)+'.mat', resp_dat)

mu_dat = {"mu_dat": mu_1_1}
savemat('mu_dat_at_temp_'+str(Temp)+'.mat', mu_dat)

J_C_dat = {"J_C_dat": J_C_1_1}
savemat('J_C_dat_at_temp_'+str(Temp)+'.mat', J_C_dat)

J_N_dat = {"J_N_dat": J_N_1_1}
savemat('J_N_dat_at_temp_'+str(Temp)+'.mat', J_N_dat)

EA_dat = {"EA_dat": EA_1_1}
savemat('EA_dat_at_temp_'+str(Temp)+'.mat', EA_dat)

Cost_basal_dat = {"Cost_basal_dat": R_B_1_1}
savemat('Cost_basal_dat_at_temp_'+str(Temp)+'.mat', Cost_basal_dat) # mug C/cell/d

Cost_enzyme_dat = {"Cost_enzyme_dat": R_E_1_1}
savemat('Cost_enzyme_dat_at_temp_'+str(Temp)+'.mat', Cost_enzyme_dat)

Cost_direct_dat = {"Cost_direct_dat": R_N2_1_1}
savemat('Cost_direct_dat_at_temp_'+str(Temp)+'.mat', Cost_direct_dat)

Cost_O2_removal_dat = {"Cost_O2_removal_dat": R_O2_1_1}
savemat('Cost_O2_removal_dat_at_temp_'+str(Temp)+'.mat', Cost_O2_removal_dat)

J_C_L_dat = {"J_C_L_dat": J_C_L_1_1}
savemat('J_C_L_dat_at_temp_'+str(Temp)+'.mat', J_C_L_dat)

J_P_L_dat = {"J_P_L_dat": J_P_L_1_1}
savemat('J_P_L_dat_at_temp_'+str(Temp)+'.mat', J_P_L_dat)

J_G_dat = {"J_G_dat": J_G_1_1}
savemat('J_G_dat_at_temp_'+str(Temp)+'.mat', J_G_dat)

J_A_dat = {"J_A_dat": J_A_1_1}
savemat('J_A_dat_at_temp_'+str(Temp)+'.mat', J_A_dat)



diff_G_dat = {"diff_G_dat": DGDx_1_1}
savemat('diff_G_dat_at_temp_'+str(Temp)+'.mat', diff_G_dat) # mug C/L/d

diff_A_dat = {"diff_A_dat": DADx_1_1}
savemat('diff_A_dat_at_temp_'+str(Temp)+'.mat', diff_A_dat) # mug C/L/d

diff_O_dat = {"diff_O_dat": DODx_1_1}
savemat('diff_O_dat_at_temp_'+str(Temp)+'.mat', diff_O_dat) # mumol O2/L/d

diff_N_dat = {"diff_N_dat": DNDx_1_1}
savemat('diff_N_dat_at_temp_'+str(Temp)+'.mat', diff_N_dat) # mumol NO3/L/d

diff_cell_O_pot_dat = {"diff_cell_O_pot_dat": V_O2_pot_1_1}
savemat('diff_cell_O_pot_dat_at_temp_'+str(Temp)+'.mat', diff_cell_O_pot_dat) # mumol O2/d/cell

diff_cell_O_dat = {"diff_cell_O_dat": V_O2_1_1}
savemat('diff_cell_O_dat_at_temp_'+str(Temp)+'.mat', diff_cell_O_dat) # mumol O2/d/cell

diff_cell_N_dat = {"diff_cell_N_dat": V_NO3_1_1}
savemat('diff_cell_N_dat_at_temp_'+str(Temp)+'.mat', diff_cell_N_dat) # mumol NO3/d/cell

diff_cell_S_dat = {"diff_cell_S_dat": V_SO4_1_1}
savemat('diff_cell_S_dat_at_temp_'+str(Temp)+'.mat', diff_cell_S_dat) # mumol SO4/d/cell



var_C_dat = {"var_C_dat": yC}
savemat('var_C_dat_at_temp_'+str(Temp)+'.mat', var_C_dat)

var_P_dat = {"var_P_dat": yP}
savemat('var_P_dat_at_temp_'+str(Temp)+'.mat', var_P_dat)

var_G_dat = {"var_G_dat": yG}
savemat('var_G_dat_at_temp_'+str(Temp)+'.mat', var_G_dat)

var_A_dat = {"var_A_dat": yA}
savemat('var_A_dat_at_temp_'+str(Temp)+'.mat', var_A_dat)

var_B_dat = {"var_B_dat": yB}
savemat('var_B_dat_at_temp_'+str(Temp)+'.mat', var_B_dat)

var_O2_dat = {"var_O2_dat": yO2}
savemat('var_O2_dat_at_temp_'+str(Temp)+'.mat', var_O2_dat)

var_NO3_dat = {"var_NO3_dat": yNO3}
savemat('var_NO3_dat_at_temp_'+str(Temp)+'.mat', var_NO3_dat)

