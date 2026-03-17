import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

# SPK Parameter
rho_0 = 1.225  # kg/m^3
c = 343  # m/s
S_D = 2.65e-4  # m^2
Bl = 1.927  # Tm
R_MS = 0.317  # Ns/m
R_E = 8.91  # ohm
C_MS = 0.22e-3  # m/N
M_MS = 0.405e-3  # kg
d = 0.1  # m
a = np.sqrt(S_D / np.pi)

# Further Parameter
V_AS = rho_0 * c**2 * (S_D**2) * C_MS
M_AS=M_MS/(S_D**2)
C_AS=(S_D**2)*C_MS
w_s=1 / np.sqrt(M_AS * C_AS)
f_s = w_s / (2 * np.pi)
Q_MS=(1/R_MS)*np.sqrt(M_MS / C_MS)
Q_ES=(R_E/Bl**2)*np.sqrt(M_MS / C_MS)
Q_TS=Q_MS*Q_ES/(Q_MS+Q_ES)

# Frequency & Ztransform Define
frequencies = np.logspace(np.log10(100), np.log10(10000), 1000)
omega = 2 * np.pi * frequencies
s = 1j * omega
p_ref = 2e-5  # 參考聲壓

# G(s) Transfunction
s_omega_ratio = s / w_s
G_s = (s_omega_ratio**2) / (s_omega_ratio**2 + (1/Q_TS)*s_omega_ratio + 1)

# Input Voltage
base_voltage_rated = 2.83  # Vrms

#  SPL & x_D 
def compute_SPL_and_displacement(E_g):
    
    p = (rho_0 / (4 * np.pi*d)) * (Bl * E_g) / (S_D * R_E * M_AS) * G_s
    SPL = 20 * np.log10(np.abs(p) / p_ref)

    x_D = E_g * np.sqrt(V_AS / (rho_0 * c**2 * S_D**2 * R_E * w_s * Q_ES)) * \
         (1 / ((s / w_s)**2 + (1 / Q_TS) * (s / w_s) + 1))
    x_D_mm = 1000 * x_D
    return SPL, x_D_mm

# Calac SPL 與 xD
SPL_rated, xD_rated = compute_SPL_and_displacement(2.83)


# SPL Plot
plt.figure(figsize=(10, 6))
plt.plot(frequencies, SPL_rated, label='Sim: Rated Voltage ', linestyle='--')

plt.xscale('log')
plt.xlim(80, 15000)
plt.ylim(40, 110)
plt.title('SPL: Simulated vs Measured')
plt.xlabel('Frequency (Hz)')
plt.ylabel('SPL (dB)')
plt.grid(True, which='both', linestyle='--', linewidth=0.5)
plt.xticks([100, 1000, 10000], ['100', '1000', '10000'])
plt.legend()
plt.tight_layout()
plt.show()

# --- Excursion Plot ---
plt.figure(figsize=(10, 6))
plt.plot(frequencies, np.abs(xD_rated), label="Simulated x_D (Rated)", color='b')
plt.axhline(0.4, color='orange', linestyle=':', linewidth=1.2, label='Xmax = 0.4 mm')
plt.xscale('log')
plt.xlim(80, 10000)
plt.ylim(0, 0.45)
plt.xlabel("Frequency (Hz)")
plt.ylabel("Displacement (mm)")
plt.title("Diaphragm Displacement: Simulated vs Measured")
plt.grid(True, which='both', linestyle='--', linewidth=0.5)
plt.legend()
plt.tight_layout()
plt.show()

