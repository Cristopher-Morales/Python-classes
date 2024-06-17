# -*- coding: utf-8 -*-
"""
Created on Tue Jun  9 15:51:14 2020

@author: 31626
"""

# -*- coding: utf-8 -*-
"""
Created on Sun Jun  7 00:15:55 2020

@author: 31626
"""

import numpy as np
import matplotlib.pyplot as plt
import math
from math import cos as cos
from scipy.optimize import curve_fit as cf

pi = math.pi
def fit(x,alpha,C):
    return C*x**alpha

#    N_x = [60,120,240];
#N_x = [60,90,120,150,180]; # 240
N_x = [90]
delta_x = [];
for i_x1 in range(np.size(N_x)):
    delta_x.append( 4*pi/(8*N_x[i_x1]));
delta_x = np.array(delta_x);
#delta_t = [pi/240,pi/360,pi/480,pi/600,pi/720];

delta_t = [pi/1080];



e_n1_LW1 = [];
e_n1_Up1 = [];
e_n1_vL1 = [];
e_n1_LW2 = [];
e_n1_Up2 = [];
e_n1_vL2 = [];
e_n1_LW5 = [];
e_n1_Up5 = [];
e_n1_vL5 = [];
e_n2_LW1 = [];
e_n2_Up1 = [];
e_n2_vL1 = [];
e_n2_LW2 = [];
e_n2_Up2 = [];
e_n2_vL2 = [];
e_n2_LW5 = [];
e_n2_Up5 = [];
e_n2_vL5 = [];  
e_ninf_LW1 = [];
e_ninf_Up1 = [];
e_ninf_vL1 = [];
e_ninf_LW2 = [];
e_ninf_Up2 = [];
e_ninf_vL2 = [];
e_ninf_LW5 = [];
e_ninf_Up5 = [];
e_ninf_vL5 = [];       

class LW_sys:
    def __init__(self, N, dt):  # N is the number of cell grading
        self.N = 8*N # Number of Cells 8 16 ...
        self.T = 10*pi # Total Time
        self.dt = dt # Time Step
        self.Period = round(self.T/self.dt)+1
        self.xmin = 0 # Left Boundary of x  DELETE
        self.xmax = 4*pi # Right Boundary of x DELETE
        self.dx = (self.xmax - self.xmin)/self.N # Spatial Interval
        self.CFL = 2*self.dt/self.dx
        self.A = np.array([[0 , -2],[-1 ,  0]]);
        self.R = np.array([[2 , -2],[1,1]])
        self.R_1 = np.linalg.inv(self.R);
        self.abseig = np.array([[2,0],[0,2]]);
        self.abseig_m = np.array([[-2,0],[0,0]]);
        self.abseig_p = np.array([[0,0],[0,2]]);
        print("CFL NUMBER IS % 1.3f" % self.CFL)
        self.x = np.arange(self.xmin-self.dx, self.xmax+(2*self.dx), self.dx)
        self.iniCond()
        

    
    def iniCond(self):
        
        self.P_0 = np.zeros(self.x.size);


        for i in range(1+round(pi/(2*self.dx)),1+round((3*pi/2 + self.dx)/self.dx), 1):
            self.P_0[i] = 1;
        for i in range(1+round(5*pi/(2*self.dx)),1+round((7*pi/2+self.dx)/self.dx),1):
            self.P_0[i] = (1+cos(2*(i-1)*self.dx))/2;
        self.U_0 = np.ones(self.x.size);
        self.P_n_LW = self.P_0.copy()
        self.P_n_Up = self.P_0.copy()
        self.P_n_vL = self.P_0.copy()
        self.U_n_LW = self.U_0.copy()
        self.U_n_Up = self.U_0.copy()
        self.U_n_vL = self.U_0.copy()
        self.P_n1_LW = self.P_0.copy()
        self.P_n1_Up = self.P_0.copy()
        self.P_n1_vL = self.P_0.copy()
        self.U_n1_LW = self.U_0.copy()
        self.U_n1_Up = self.U_0.copy()
        self.U_n1_vL = self.U_0.copy()
        


    def solution(self):
        tc = 0
#        
        def phi(a,b):

            
                with np.errstate(divide='ignore',invalid='ignore'):
                    theta = np.divide(a,b);
                    if np.isinf(theta):
                        c = 2;
                    elif np.isnan(theta):
                        c = 1;
                    elif theta>10^300:
                        c = 2;
                    else:   
                        c = (theta+abs(theta))/(1+abs(theta)); # -inf inf NaN
        
                return c

            
        for i in range(self.Period):
            plt.clf()
            
            # Lax-Wendroff
#            [self.w1_n_LW , self.w2_n_LW] = np.dot(self.R_1,[self.P_n_LW , self.U_n_LW]);
            for j in range(self.N+2):
                
#                self.P_n1_LW[j] = self.P_n_LW[j] + self.CFL*2*(self.U_n_LW[j+1]-self.U_n_LW[j-1])+((self.CFL**2)/2)*4*(self.P_n_LW[j+1]-2*self.P_n_LW[j]+self.P_n_LW[j-1])
#
#                self.U_n1_LW[j] = self.U_n_LW[j] + 0.5*self.CFL*(self.P_n_LW[j+1]-self.P_n_LW[j-1])+((self.CFL**2)/2)*4*(self.U_n_LW[j+1]-2*self.U_n_LW[j]+self.U_n_LW[j-1])
                Q_n1_LW = np.zeros([2,1]);
                Q_n_LW = np.zeros([2,1]);
                Q_n_LW[0,0] = self.P_n_LW[j];
                Q_n_LW[1,0] = self.U_n_LW[j];
                Q_n_LWm = np.zeros([2,1]);
                Q_n_LWm[0,0] = self.P_n_LW[j-1];
                Q_n_LWm[1,0] = self.U_n_LW[j-1];
                Q_n_LWp = np.zeros([2,1]);
                Q_n_LWp[0,0] = self.P_n_LW[j+1];
                Q_n_LWp[1,0] = self.U_n_LW[j+1];
                
#                self.P_n1_Up[j] = self.P_n_Up[j] + 4*(self.dt/self.dx)* (self.U_n_Up[j]-self.U_n_Up[j-1])
#
#                self.U_n1_Up[j] = self.U_n_Up[j] + (self.dt/self.dx)*(self.P_n_Up[j]-self.P_n_Up[j-1])
                Q_n1_LW = Q_n_LW - (self.dt/(2*self.dx))*self.A*(Q_n_LWp-Q_n_LWm)+0.5*((self.dt/(self.dx))**2)*self.A*self.A*(Q_n_LWm-2*Q_n_LW+Q_n_LWp)
                
                self.P_n1_LW[j]= Q_n1_LW[0,0];
                self.U_n1_LW[j]= Q_n1_LW[1,0];
            
##            
                
                
            self.P_n_LW = self.P_n1_LW.copy()
            self.U_n_LW = self.U_n1_LW.copy()
 
          
#            # Upwind
#            
            for j in range(self.N+2):
                
                Q_n1_Up = np.zeros([2,1]);
                Q_n_Up = np.zeros([2,1]);
                Q_n_Up[0,0] = self.P_n_Up[j];
                Q_n_Up[1,0] = self.U_n_Up[j];
                Q_n_Upm = np.zeros([2,1]);
                Q_n_Upm[0,0] = self.P_n_Up[j-1];
                Q_n_Upm[1,0] = self.U_n_Up[j-1];
                Q_n_Upp = np.zeros([2,1]);
                Q_n_Upp[0,0] = self.P_n_Up[j+1];
                Q_n_Upp[1,0] = self.U_n_Up[j+1];
                
#                self.P_n1_Up[j] = self.P_n_Up[j] + 4*(self.dt/self.dx)* (self.U_n_Up[j]-self.U_n_Up[j-1])
#
#                self.U_n1_Up[j] = self.U_n_Up[j] + (self.dt/self.dx)*(self.P_n_Up[j]-self.P_n_Up[j-1])
                Q_n1_Up = Q_n_Up - (self.dt/self.dx)*((self.R*self.abseig*self.R_1)*Q_n_Up+(self.R*self.abseig_m*self.R_1)*Q_n_Upp-(self.R*self.abseig_p*self.R_1)*Q_n_Upm)
                
                self.P_n1_Up[j]= Q_n1_Up[0,0];
                self.U_n1_Up[j]= Q_n1_Up[1,0];
            
##            
#            
            self.P_n_Up = self.P_n1_Up.copy()
            self.U_n_Up = self.U_n1_Up.copy()
#            
#            self.w1_n_Up = self.w1_n1_Up.copy()
#            self.w2_n_Up = self.w2_n1_Up.copy()
#            
#            [self.P_n_Up, self.U_n_Up] = np.dot(self.R,[self.w1_n_Up, self.w2_n_Up]);
#            
#             van Leer
            

            
            
            for j in range(self.N+1):
                
                
                a1_m = (1/4)*((self.P_n_vL[j]-self.P_n_vL[j-1])+2*(self.U_n_vL[j]-self.U_n_vL[j-1]))
                a1_p = (1/4)*((self.P_n_vL[j+1]-self.P_n_vL[j])+2*(self.U_n_vL[j+1]-self.U_n_vL[j]))
                a1_pp = (1/4)*((self.P_n_vL[j+2]-self.P_n_vL[j+1])+2*(self.U_n_vL[j+2]-self.U_n_vL[j+1]))# +0.5
                a2_m = (1/4)*(-1*(self.P_n_vL[j]-self.P_n_vL[j-1])+2*(self.U_n_vL[j]-self.U_n_vL[j-1]))
                a2_p = (1/4)*(-1*(self.P_n_vL[j+1]-self.P_n_vL[j])+2*(self.U_n_vL[j+1]-self.U_n_vL[j]))
                a2_mm = (1/4)*(-1*(self.P_n_vL[j-1]-self.P_n_vL[j-2])+2*(self.U_n_vL[j-1]-self.U_n_vL[j-2]))
#                th1_m = a1_p/a1_m
#                th2_m = a2_p/a2_m
                
                Q_n1_vL = np.zeros([2,1]);
                Q_n_vL = np.zeros([2,1]);
                Q_n_vLp = np.zeros([2,1]);
                Q_n_vLm = np.zeros([2,1]);
                Q_n_vL[0,0] = self.P_n_vL[j];
                Q_n_vL[1,0] = self.U_n_vL[j];
                Q_n_vLp[0,0] = self.P_n_vL[j+1];
                Q_n_vLp[1,0] = self.U_n_vL[j+1];
                Q_n_vLm[0,0] = self.P_n_vL[j-1];
                Q_n_vLm[1,0] = self.U_n_vL[j-1];
                r_1 = np.array([[2],[1]]);
                r_2 = np.array([[-2],[1]]);
#                if j == 0:
#                    self.w1_n_vL[j-1] = self.w1_n_vL[self.N+1]
                
                Q_n1_vL = Q_n_vL-(self.dt/self.dx)*((self.R*self.abseig*self.R_1)*Q_n_vL+self.R*self.abseig_m*self.R_1*Q_n_vLp-self.R*self.abseig_p*self.R_1*Q_n_vLm + 0.5*(self.R*self.abseig*self.R_1)*(np.eye(2)-(self.dt/self.dx)*(self.R*self.abseig*self.R_1))*(a1_p*phi(a1_pp,a1_p)-a1_m*phi(a1_p,a1_m))*r_1+(a2_p*phi(a2_m,a2_p)-a2_m*phi(a2_mm,a2_m))*r_2)
                
                self.P_n1_vL[j]= Q_n1_vL[0,0];
                self.U_n1_vL[j]= Q_n1_vL[1,0];
##                print(F_im)
#                
##            [self.P_n1_vL, self.U_n1_vL] = np.dot(self.R,[self.w1_n1_vL, self.w2_n1_vL]);
#            
            self.P_n_vL = self.P_n1_vL.copy()
            self.U_n_vL = self.U_n1_vL.copy()
##            
#            self.w1_n_vL = self.w1_n1_vL.copy()
#            self.w2_n_vL = self.w1_n1_vL.copy()
#            
#            [self.P_n_vL, self.U_n_vL] = np.dot(self.R,[self.w1_n_vL, self.w2_n_vL]);
#            print(self.w2_n_vL)
            # Periodic boundary conditions
#            self.w1_n_LW[0] = self.w1_n_LW[self.N+1]
#            self.w1_n_LW[self.N+2] = self.w1_n_LW[1]
#            self.w2_n_LW[0] = self.w2_n_LW[self.N+1]
#            self.w2_n_LW[self.N+2] = self.w2_n_LW[1]
#            self.w1_n_Up[0] = self.w1_n_Up[self.N+1]
#            self.w1_n_Up[self.N+2] = self.w1_n_Up[1]
#            self.w2_n_Up[0] = self.w2_n_Up[self.N+1]
#            self.w2_n_Up[self.N+2] = self.w2_n_Up[1]
#            self.w1_n_vL[0] = self.w1_n_vL[self.N+1]
#            self.w1_n_vL[self.N+2] = self.w1_n_vL[1]
#            self.w2_n_vL[0] = self.w2_n_vL[self.N+1]
#            self.w2_n_vL[self.N+2] = self.w2_n_vL[1]
#            self.w1_n_LW[1] = self.w1_n_LW[0]
#            self.w2_n_LW[self.N+2] = self.w2_n_LW[1]
#            self.w1_n_Up[0] = self.w1_n_Up[self.N+1]
#            self.w2_n_Up[self.N+2] = self.w2_n_Up[1]
#            self.w1_n_vL[0] = self.w1_n_vL[self.N+1]
#            self.w2_n_vL[self.N+2] = self.w2_n_vL[1]
            
            
            

#            self.w1_n_LW[0] = self.w1_n_LW[self.N+2]
#            self.w2_n_LW[self.N+2] = self.w2_n_LW[0]
#            self.w1_n_Up[0] = self.w1_n_Up[self.N+2]
#            self.w2_n_Up[self.N+2] = self.w2_n_Up[0]
#            self.w1_n_vL[0] = self.w1_n_vL[self.N+2]
#            self.w2_n_vL[self.N+2] = self.w2_n_vL[0]
#            
            
            
            self.P_n_LW[0] = self.P_n_LW[self.N+1]
            self.P_n_LW[self.N+2] = self.P_n_LW[1]
            self.U_n_LW[0] = self.U_n_LW[self.N+1]
            self.U_n_LW[self.N+2] = self.U_n_LW[1]
            
            self.P_n_Up[0] = self.P_n_Up[self.N+1]
            self.P_n_Up[self.N+2] = self.P_n_Up[1]
            self.U_n_Up[0] = self.U_n_Up[self.N+1]
            self.U_n_Up[self.N+2] = self.U_n_Up[1]     
            
            self.P_n_vL[0] = self.P_n_vL[self.N+1]
            self.P_n_vL[self.N+2] = self.P_n_vL[1]
            self.U_n_vL[0] = self.U_n_vL[self.N+1]
            self.U_n_vL[self.N+2] = self.U_n_vL[1] 
#            self.P_n_LW[-1] = self.P_n_LW[self.N+1]
#            self.U_n_LW[-1] = self.U_n_LW[self.N+1]
#            self.P_n_Up[0] = self.P_n_Up[self.N+1]
#            self.U_n_Up[self.N+2] = self.U_n_Up[1]
#            self.P_n_vL[0] = self.P_n_vL[self.N+1]
#            self.U_n_vL[self.N+2] = self.P_n_vL[1]
#            
            
            # Exact Solution
            
            
#            for n in range(5,1,-1):
#                while self.dt*i-n*4*pi>0:
#                    nt = n;
#                    
#            tp = self.dt*i-nt*4*pi
#            q_ex = np.zeros(self.x.size);
#            for i in range(1+round(pi/(2*self.dx)+tc/self.dx),1+round((3*pi/2 + self.dx)/self.dx+tc/self.dx), 1):
#                q_ex[i] = 1;
#            for i in range(1+round(5*pi/(2*self.dx)+tc/self.dx),1+round((7*pi/2+self.dx)/self.dx+tc/self.dx),1):
#                q_ex[i] = (1+cos(2*(i-1)*self.dx-2*tc))/2;
#            
#            q_0 = np.zeros(self.x.size);
#            for m in range(1+round(pi/(2*self.dx)),1+round((3*pi/2 + self.dx)/self.dx), 1):
#                q_0[m] = 1;
#            for m in range(1+round(5*pi/(2*self.dx)),1+round((7*pi/2+self.dx)//self.dx),1):
#                q_0[m] = (1+cos(2*(m-1)*self.dx))/2;
         #######################################   
#            for n in [1,2,5]:
#                k = 1
#                
##                for j in range(self.N+2):
#    
#                
#                if abs(i*self.dt - (n)*2*pi)<0.001:
#                    [self.P_n_LW, self.U_n_LW] = np.dot(self.R,[self.w1_n_LW, self.w2_n_LW]);
#                    E_n_LW = self.q_n_LW - self.q_0;
#                    E_n_Up = self.q_n_Up - self.q_0;
#                    E_n_vL = self.q_n_vL - self.q_0;
                    
                    # Calculate Error
#                    if n == 1:
#                        
#                        e_n1_LW1.append(np.linalg.norm(delta_x[i_x]*E_n_LW,1));
#                        e_n1_Up1.append(np.linalg.norm(delta_x[i_x]*E_n_Up,1));
#                        e_n1_vL1.append(np.linalg.norm(delta_x[i_x]*E_n_vL,1));
#                        e_n2_LW1.append(np.linalg.norm(delta_x[i_x]*E_n_LW,2));
#                        e_n2_Up1.append(np.linalg.norm(delta_x[i_x]*E_n_Up,2));
#                        e_n2_vL1.append(np.linalg.norm(delta_x[i_x]*E_n_vL,2));
#                        e_ninf_LW1.append(np.linalg.norm(delta_x[i_x]*E_n_LW,np.inf));
#                        e_ninf_Up1.append(np.linalg.norm(delta_x[i_x]*E_n_Up,np.inf));
#                        e_ninf_vL1.append(np.linalg.norm(delta_x[i_x]*E_n_vL,np.inf));
#                    if n == 2:
#
#                        e_n1_LW2.append(np.linalg.norm(delta_x[i_x]*E_n_LW,1));
#                        e_n1_Up2.append(np.linalg.norm(delta_x[i_x]*E_n_Up,1));
#                        e_n1_vL2.append(np.linalg.norm(delta_x[i_x]*E_n_vL,1));
#                        e_n2_LW2.append(np.linalg.norm(delta_x[i_x]*E_n_LW,2));
#                        e_n2_Up2.append(np.linalg.norm(delta_x[i_x]*E_n_Up,2));
#                        e_n2_vL2.append(np.linalg.norm(delta_x[i_x]*E_n_vL,2));
#                        e_ninf_LW2.append(np.linalg.norm(delta_x[i_x]*E_n_LW,np.inf));
#                        e_ninf_Up2.append(np.linalg.norm(delta_x[i_x]*E_n_Up,np.inf));
#                        e_ninf_vL2.append(np.linalg.norm(delta_x[i_x]*E_n_vL,np.inf));
#                    if n == 5:
#                        
#                        e_n1_LW5.append(np.linalg.norm(delta_x[i_x]*E_n_LW,1));
#                        e_n1_Up5.append(np.linalg.norm(delta_x[i_x]*E_n_Up,1));
#                        e_n1_vL5.append(np.linalg.norm(delta_x[i_x]*E_n_vL,1));
#                        e_n2_LW5.append(np.linalg.norm(delta_x[i_x]*E_n_LW,2));
#                        e_n2_Up5.append(np.linalg.norm(delta_x[i_x]*E_n_Up,2));
#                        e_n2_vL5.append(np.linalg.norm(delta_x[i_x]*E_n_vL,2));
#                        e_ninf_LW5.append(np.linalg.norm(delta_x[i_x]*E_n_LW,np.inf));
#                        e_ninf_Up5.append(np.linalg.norm(delta_x[i_x]*E_n_Up,np.inf));
#                        e_ninf_vL5.append(np.linalg.norm(delta_x[i_x]*E_n_vL,np.inf));
#            print(self.P_n_LW)        
            plt.figure(figsize=(8,6))
            plt.plot(self.x, self.P_0, 'r-', label="Exact solution")
            plt.plot(self.x, self.P_n_LW, 'b', label="Lax-Wendroff")
            plt.plot(self.x, self.P_n_Up, 'g', label="Upwind")
            plt.plot(self.x, self.P_n_vL, 'k', label='van Leer')
            plt.axis((self.xmin-0.12, self.xmax+0.12, -0.2, 1.4))
            plt.grid(True)
            plt.xlabel("x")
            plt.ylabel("q")
            plt.legend(loc=1, fontsize=12) #
#            plt.title("Time = %d Period Pressure Plot" % n)
#            k += 1
            plt.pause(0.01)
            plt.show()
            
            plt.figure(figsize=(8,6))
            plt.plot(self.x, self.U_0, 'r-', label="Exact solution")
            plt.plot(self.x, self.U_n_LW, 'b', label="Lax-Wendroff")
            plt.plot(self.x, self.U_n_Up, 'g', label="Upwind")
            plt.plot(self.x, self.U_n_vL, 'k', label='van Leer')
            plt.axis((self.xmin-0.12, self.xmax+0.12, -0.2, 1.4))
            plt.grid(True)
            plt.xlabel("x")
            plt.ylabel("q")
            plt.legend(loc=1, fontsize=12) #
#            plt.title("Time = %d Period Velocity Plot" % n)
#            k += 1
            plt.pause(0.01)
            plt.show()
        tc += self.dt
    

    
##    N_x = [60,120,240];
#N_x = [150,300,600]; # 240
#
#delta_x = [];
#for i_x in range(np.size(N_x)):
#    delta_x.append( 4*pi/(8*N_x[i_x]));
#delta_x = np.array(delta_x);
##delta_t = [pi/400,pi/800,pi/1200];
#for i_x in range(np.size(N_x)):
#    sim = LW_adv(N_x[i_x], pi/1500)
#    sim.solution()

for i_x in range(np.size(N_x)):
    sim = LW_sys(N_x[i_x], delta_t[i_x])
    sim.solution()    
    
#    
#e_n1_LW1 = np.array(e_n1_LW1)
#e_n1_Up1 = np.array(e_n1_Up1);
#e_n1_vL1 = np.array(e_n1_vL1);
#e_n1_LW2 = np.array(e_n1_LW2);
#e_n1_Up2 = np.array(e_n1_Up2);
#e_n1_vL2 = np.array(e_n1_vL2);
#e_n1_LW5 = np.array(e_n1_LW5);
#e_n1_Up5 = np.array(e_n1_Up5);
#e_n1_vL5 = np.array(e_n1_vL5);
#e_n2_LW1 = np.array(e_n2_LW1);
#e_n2_Up1 = np.array(e_n2_Up1);
#e_n2_vL1 = np.array(e_n2_vL1);
#e_n2_LW2 = np.array(e_n2_LW2);
#e_n2_Up2 = np.array(e_n2_Up2);
#e_n2_vL2 = np.array(e_n2_vL2);
#e_n2_LW5 = np.array(e_n2_LW5);
#e_n2_Up5 = np.array(e_n2_Up2);
#e_n2_vL5 = np.array(e_n2_vL2);  
#e_ninf_LW1 = np.array(e_ninf_LW1);
#e_ninf_Up1 = np.array(e_ninf_Up1);
#e_ninf_vL1 = np.array(e_ninf_vL1);
#e_ninf_LW2 = np.array(e_ninf_LW2);
#e_ninf_Up2 = np.array(e_ninf_Up2);
#e_ninf_vL2 = np.array(e_ninf_vL2);
#e_ninf_LW5 = np.array(e_ninf_LW5);
#e_ninf_Up5 = np.array(e_ninf_Up5);
#e_ninf_vL5 = np.array(e_ninf_vL5);
#
## 1-norm
#
#para_n1_LW1,conv_n1_LW1 = cf(fit,delta_x,e_n1_LW1);
#para_n1_Up1,conv_n1_Up1 = cf(fit,delta_x,e_n1_Up1);
#para_n1_vL1,conv_n1_vL1 = cf(fit,delta_x,e_n1_vL1);
#C_n1_LW1,Alpha_n1_LW1 = para_n1_LW1[0],para_n1_LW1[1];
#C_n1_Up1,Alpha_n1_Up1 = para_n1_Up1[0],para_n1_Up1[1];
#C_n1_vL1,Alpha_n1_vL1 = para_n1_vL1[0],para_n1_vL1[1];
#plt.figure(figsize=(8,6))
#plt.plot(np.log10(delta_x), np.log10(fit(delta_x, *para_n1_LW1)), 'b-',\
#         np.log10(delta_x), np.log10(e_n1_LW1), 'rx' ,label="Lax-Wendroff，slope=%1.3f" %C_n1_LW1)
#plt.plot(np.log10(delta_x), np.log10(fit(delta_x, *para_n1_Up1)), 'g-',\
#         np.log10(delta_x), np.log10(e_n1_Up1), 'rx' ,label="Upwind，slope=%1.3f" %C_n1_Up1)
#plt.plot(np.log10(delta_x), np.log10(fit(delta_x, *para_n1_vL1)), 'm-',\
#         np.log10(delta_x), np.log10(e_n1_vL1), 'rx' ,label="van Leer，slope=%1.3f" %C_n1_vL1)
#plt.grid(True)
#plt.xlabel("Grid Size")
#plt.ylabel("Error")
#plt.legend(loc=1, fontsize=12) #
#plt.title("1-norm at Second Period.")
#plt.show()
#
## 2-norm
#
#para_n2_LW1,conv_n2_LW1 = cf(fit,delta_x,e_n2_LW1);
#para_n2_Up1,conv_n2_Up1 = cf(fit,delta_x,e_n2_Up1);
#para_n2_vL1,conv_n2_vL1 = cf(fit,delta_x,e_n2_vL1);
#C_n2_LW1,Alpha_n2_LW1 = para_n2_LW1[0],para_n2_LW1[1];
#C_n2_Up1,Alpha_n2_Up1 = para_n2_Up1[0],para_n2_Up1[1];
#C_n2_vL1,Alpha_n2_vL1 = para_n2_vL1[0],para_n2_vL1[1];
#plt.figure(figsize=(8,6))
#plt.plot(np.log10(delta_x), np.log10(fit(delta_x, *para_n2_LW1)), 'b-',\
#         np.log10(delta_x), np.log10(e_n2_LW1), 'rx' ,label="Lax-Wendroff，slope=%1.3f" %C_n2_LW1)
#plt.plot(np.log10(delta_x), np.log10(fit(delta_x, *para_n2_Up1)), 'g-',\
#         np.log10(delta_x), np.log10(e_n2_Up1), 'rx' ,label="Upwind，slope=%1.3f" %C_n2_Up1)
#plt.plot(np.log10(delta_x), np.log10(fit(delta_x, *para_n2_vL1)), 'm-',\
#         np.log10(delta_x), np.log10(e_n2_vL1), 'rx' ,label="van Leer，slope=%1.3f" %C_n2_vL1)
#
#plt.grid(True)
#plt.xlabel("Grid Size")
#plt.ylabel("Error")
#plt.legend(loc=1, fontsize=12) #
#plt.title("2-norm at Second Period.")
#plt.show()
#
#
## inf-norm
#
#para_ninf_LW1,conv_ninf_LW1 = cf(fit,delta_x,e_ninf_LW1);
#para_ninf_Up1,conv_ninf_Up1 = cf(fit,delta_x,e_ninf_Up1);
#para_ninf_vL1,conv_ninf_vL1 = cf(fit,delta_x,e_ninf_vL1);
#C_ninf_LW1,Alpha_ninf_LW1 = para_ninf_LW1[0],para_ninf_LW1[1];
#C_ninf_Up1,Alpha_ninf_Up1 = para_ninf_Up1[0],para_ninf_Up1[1];
#C_ninf_vL1,Alpha_ninf_vL1 = para_ninf_vL1[0],para_ninf_vL1[1];
#plt.figure(figsize=(8,6))
#plt.plot(np.log10(delta_x), np.log10(fit(delta_x, *para_ninf_LW1)), 'b-',\
#         np.log10(delta_x), np.log10(e_ninf_LW1), 'rx' ,label="Lax-Wendroff，slope=%1.3f" %C_ninf_LW1)
#plt.plot(np.log10(delta_x), np.log10(fit(delta_x, *para_ninf_Up1)), 'g-',\
#         np.log10(delta_x), np.log10(e_ninf_Up1), 'rx' ,label="Upwind，slope=%1.3f" %C_ninf_Up1)
#plt.plot(np.log10(delta_x), np.log10(fit(delta_x, *para_ninf_vL1)), 'm-',\
#         np.log10(delta_x), np.log10(e_ninf_vL1), 'rx' ,label="van Leer，slope=%1.3f" %C_ninf_vL1)
#plt.grid(True)
#plt.xlabel("Grid Size")
#plt.ylabel("Error")
#plt.legend(loc=1, fontsize=12) #
#plt.title("inf-norm at Second Period.")
#plt.show()

#fig, ax = plt.subplots(1)
#ax.plot(np.log10(delta_x),np.log10(e_n1_LW1))
#plt.show()
