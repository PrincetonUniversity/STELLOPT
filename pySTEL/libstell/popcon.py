"""
This library provides a python class for creating stellarator POPCON plots
"""

import numpy as np
import sys

# Constants
EC = 1.602176634E-19 # Electron charge [C]

# popcon Class
class POPCON:
    
    def __init__(self, B, a, R, iota, plasma_classes, popcon_title = 'POPCON'):
        # navg is an array given in units of m^-3
        # Tavg_array is an array in units of eV
        
        from scipy import integrate
        
        self.plasma_list = plasma_classes      
        
        self.popcon_title = popcon_title
        
        # Differential volume. This works well for large aspect-ratio; not so good otherwise
        self.dVdrho = lambda rho: 4*np.pi*np.pi*R*rho*a*a
        self.Volume, _ = integrate.quad(lambda x: self.dVdrho(x),0.0,1.0)
        
        #n_avg and T_avg have the same size as plasma_classes
        # note that they are the averages of the electron profile
        self.n_avg = self.get_averaged_density()
        self.T_avg = self.get_averaged_temperature()
        
        # set tauiss04
        self.tauiss04 = lambda P: 0.134*a**2.28*R**0.64*(P/1e6)**-0.61*(self.n_avg/1e19)**0.54*B**0.84*iota**0.41
        
        # set Sudo limit (this will be used only for plotting)
        self.sudo_max = lambda P: 1.25*0.25E20*np.sqrt(P*B/(a*a*R*1e6))
        
        # sets self.RHS, which is a 2D array w/ same shape as n_avg and T_avg
        # RHS has:
        # - Bremstrahlung radiation losses
        # - alpha heating
        self.set_RHS()
        
        # finds Pext for each point of self.RHS
        self.find_Pext()      

        # Plot popcon
        self.plot_popcon()
    
    def get_averaged_density(self):
        
        from scipy.integrate import trapezoid
        
        rho = np.linspace(0,1,100)
        
        n_avg = np.zeros_like(self.plasma_list,dtype=float)
        
        for idx,plasma in np.ndenumerate(self.plasma_list):
            
            int_n_dV = trapezoid( self.dVdrho(rho)*plasma.get_density('electrons',rho), rho )
                
            n_avg[idx] = int_n_dV / self.Volume
            
        return n_avg
    
    def get_averaged_temperature(self):
        
        from scipy.integrate import trapezoid
        
        rho = np.linspace(0,1,100)
        
        T_avg = np.zeros_like(self.plasma_list,dtype=float)
        
        for idx,plasma in np.ndenumerate(self.plasma_list):
                
            int_T_dV = trapezoid( self.dVdrho(rho)*plasma.get_temperature('electrons',rho), rho )
                
            T_avg[idx] = int_T_dV / self.Volume
            
        return T_avg      
        
    
    def set_RHS(self):
        
        # Bremstrahlung
        RHS_1 = self.get_Bremsstrahlung()
        
        # Alpha power (notice the minus sign)
        RHS_2 = -self.get_alpha_power()
        
        self.RHS = RHS_1 + RHS_2
        
        
    def get_Bremsstrahlung(self):
        
        from scipy.integrate import trapezoid
            
        rho = np.linspace(0,1,100)
        
        # c = 3e8
        # h = 6.626e-34
        # eps0 = 8.8541878188E-12
        # me = 9.1e-31
        # cte = (np.sqrt(2)/(3*np.pi**2.5)) * EC**6 / (eps0**3 * c**3 * h *me**1.5) * (1e20)**2 * np.sqrt(EC*1e3)
        # print(f'cte={cte}')
        
        # Bremsstrahlung power should have the same dimension as plasma_list
        PB = np.zeros_like(self.plasma_list,dtype=float)
        
        for idx,plasma in np.ndenumerate(self.plasma_list):
            
            ne = plasma.get_density('electrons',rho)
            Te = plasma.get_temperature('electrons',rho)
            
            # Compute Zeff
            Zeff = 0.0
            for ion_s in plasma.ion_species:
                Zeff += plasma.get_density(ion_s,rho) * plasma.Zcharge[ion_s]**2
            Zeff = Zeff / ne
                
            # formula according to Eq. (3.43), page 56, in Freidberg
            CB = 5.35e3
            n20 = ne / 1e20
            Tk = Te / 1e3
            
            SB = CB * Zeff * n20**2 * np.sqrt(Tk) # W/m^3
            
            #integrate in volume
            PB[idx] = trapezoid(self.dVdrho(rho) * SB, rho) # W
            
            self.PB = PB
            
        return PB
    
    def set_Ethermal_plasma(self):
        # sets plasma thermal energy [J] for each plasma_class in self.plasma_list
        
        from scipy.integrate import trapezoid
        
        rho = np.linspace(0,1,100)
        
        Ethermal = np.zeros_like(self.plasma_list,dtype=float)
        
        for idx,plasma in np.ndenumerate(self.plasma_list):
            
            dEdV = 0.0
            for species in plasma.list_of_species:
                
                dEdV += plasma.get_density(species,rho) * plasma.get_temperature(species,rho) * EC # Joule / m^3
                
            Ethermal[idx] = trapezoid(self.dVdrho(rho) * dEdV, rho) # Joule
            
        self.Eplasma_thermal = Ethermal
    
    def get_alpha_power(self):
        # P_alpha = E_alpha * integral(dV * nD * nT *sigmav )
        
        from fusion import FUSION
        from scipy.integrate import trapezoid
        
        # Fusion Class
        fusion = FUSION()   
        
        rho = np.linspace(0,1,100)
        P_alpha = np.zeros_like(self.plasma_list,dtype=float)
        
        for idx,plasma in np.ndenumerate(self.plasma_list):
        
            nD = plasma.get_density('deuterium', rho)
            nT = plasma.get_density('tritium', rho)
            
            Ti = 0.5* ( plasma.get_temperature('deuterium', rho) + plasma.get_temperature('tritium', rho) )
            
            sigmav = [fusion.sigmaBH(ti,'DT') for ti in Ti]
            
            S_alpha = nD * nT * sigmav *  fusion.E_DT_He # W/m^3 
            
            P_alpha[idx] = trapezoid(self.dVdrho(rho) * S_alpha,rho) # W
        
        self.P_alpha = P_alpha #this will be needed to computed tauISS04, hence why I keep it
            
        return P_alpha
    
    def find_Pext(self):
        # for each plasma_class (i.e., each point in the POPCON plot) solve the eq:
        # Pext - Eplasma/tauISS04(P=Pext+Palpha) = RHS
        
        from scipy.optimize import fsolve
        import matplotlib.pyplot as plt
        
        self.set_Ethermal_plasma()
        
        P_ext = np.zeros_like(self.plasma_list,dtype=float)
        
        for idx,plasma in np.ndenumerate(self.plasma_list):
            
            Eth = self.Eplasma_thermal[idx]
            Palpha = self.P_alpha[idx]
            
            rhs = self.RHS[idx]
            
            f_zero = lambda x: x - Eth/self.tauiss04(x+Palpha)[idx] - rhs
            
            # tt = np.linspace(0,185e6,1000)
            # plt.plot(tt,f_zero(tt),'.-')
            # plt.title(f'n_avg={self.n_avg[0]},   T_avg={self.T_avg[0]}')
            # plt.grid()
            # plt.show()
            
            Pext = fsolve(f_zero,10E6)
            
            P_ext[idx] = Pext[0]
            
        self.P_ext = P_ext    
        
            
    def plot_popcon(self):
        
        import matplotlib.pyplot as plt
        
        plt.rc('font', size=18)
        fig, ax = plt.subplots(figsize=(11,8))
        
        # sets negative values of P_ext to 0 and converts to MW
        P_MW = self.P_ext.clip(min=0) / 1e6
        
        # convert density to n20 and temperature to Tk
        n20 = self.n_avg / 1e20
        Tk  = self.T_avg / 1e3
        
        # fusion power = P_alpha + P_neutron = P_alpha + (E_neutron/E_alpha)*P_alpha = 5*P_alpha
        P_fusion_GW = 5*self.P_alpha / 1e9
        
        cntrf = ax.pcolor(Tk.transpose(),n20.transpose(),P_MW.transpose(),cmap='hot_r')#,levels=50)
        cntr = ax.contour(Tk.transpose(),n20.transpose(),P_MW.transpose(),levels=[0,10,20,30,50,70],linestyles='dashed')
        ax.contour(Tk.transpose(),n20.transpose(),P_fusion_GW.transpose(),levels=[3.0],linestyles='solid',colors='red')
        ax.clabel(cntr, inline=True, fontsize=17)
        
        fig.colorbar(cntrf,label='P_ext [MW]')
        ax.set_xlabel(r'$\left<T_e\right>$ [keV]')
        ax.set_ylabel(r'$\left<n_e\right>$ (x10$^{20}$ m$^{-3}$)')
        ax.set_title(f'{self.popcon_title}')
        ax.text(6.75, 1.05, r'$P_{\text{fusion}}=3$GW', color='red', fontsize=16)
        
        # overlay Sudo limit
        n_max_20 = self.sudo_max(P_MW*1e6 + self.P_alpha) / 1e20
        #get peak values for all points in the plot
        n0_20 = [[plasma.get_density('electrons', 0.0)/1e20 for plasma in row] for row in self.plasma_list]
        ax.contour(Tk.transpose(),n20.transpose(),(n_max_20-n0_20).transpose(),levels=[0.0],linestyles='solid',colors='green')
        ax.text(6.5, 1.7, r'$n_0/n_{\text{Sudo}}=1.25$', color='green', fontsize=16)
        
        # plot star
        #ax.scatter(1.02, 0.075, s=320, marker='*', color='blue', zorder=3)
        
        # fig, ax = plt.subplots(figsize=(11,8))
        # cntrf=ax.pcolor(Tk.transpose(),n20.transpose(),(n_max_20-n0_20).transpose(),cmap='seismic')#,levels=50)
        # fig.colorbar(cntrf)
        
        plt.show()
        
        
        
        
        
        
            
            

        
        
        
        
            
        
        
        
        
        
        




# Main routine
if __name__=="__main__":
	import sys
	sys.exit(0)