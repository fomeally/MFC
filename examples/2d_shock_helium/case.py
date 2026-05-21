#!/usr/bin/env python3
import math, json


p0 = 101325.0
ps = 159080.98
vs = -144.0
rho_s = 1.66
T0 = 293.0
Rbar = 8314.462618


#Air props
gamma_air = 1.4
W_air = 28.97
rho_air = p0*W_air / (Rbar * T0)
rho_air = 1.2062
Gamma_air = 1.0 / (gamma_air - 1.0)
cv_air = Rbar * Gamma_air / W_air
cv_air = 717.5
cp_air = gamma_air*cv_air
T0_air = 0.0
h0_air = 0.0
mu_air = 1.0 / (1.9e-5)
kappa_air = 1.0 / (0.7e-5)


#Helium props
gamma_he = 1.6451
W_he = 4.0026
rho_he = p0*W_he / (Rbar * T0)
rho_he = 0.2204
Gamma_he = 1.0 / (gamma_he - 1.0)
cv_he = Rbar * Gamma_he / W_he
cv_he = 2430.35
cp_he = gamma_he*cv_he
T0_he = 0.0
h0_he = 0.0
mu_he = 1.0 / (1.9e-5)
kappa_he = 1.0 / (1.0e-12)

# Binary diffusion coefficients
D11 = 0.0e0
D12 = 73.35e-6
D21 = D12
D22 = 0.0e0

Lx = 22.25e-2
Ly = 8.9e-2
xs = 16.8e-2

Nx = 2000
Ny = 356


# Configuring case dictionary
print(
    json.dumps(
        {
            # Logistics
            "run_time_info": "T",
            # Computational Domain Parameters
            "x_domain%beg": 0.0,
            "x_domain%end": Lx,
            "y_domain%beg": 0.0,
            "y_domain%end": Ly,
            "stretch_x": "F",
            "cyl_coord": "F",
            "m": int(Nx),
            "n": int(Ny),
            "p": 0,
            "dt": 5.0e-8,
            "t_step_start": 0,
            "t_step_stop": 1,
            "t_step_save": 1,
            # "t_step_stop": 13720,
            # "t_step_save": 20,
            # "t_step_print": 137,
            # "t_step_stop": 10000000,
            # "t_step_save": 50000,
            # Simulation Algorithm Parameters
            "num_patches": 3,
            "model_eqns": 2,
            "alt_soundspeed": "F",
            "diffusion": "F",
            "num_fluids": 2,
            "mpp_lim": "T",
            "mixture_err": "T",
            "time_stepper": 3,
            "weno_order": 5,
            "weno_eps": 1.0e-16,
            "weno_Re_flux": "F",
            "weno_Dif_flux": "F",
            "Dif_fv" : "T",
	        "dif_order": 2,
            "small_num_dif": 1.0e-8,
            "weno_avg": "F",
            "mapped_weno": "T",
            "null_weights": "F",
            "mp_weno": "F",
            "riemann_solver": 2,
            "wave_speeds": 1,
            "avg_state": 2,
            "bc_x%beg": -6,
            "bc_x%end": -6,
            "bc_y%beg": -6,
            "bc_y%end": -6,
       	    "viscous": "T",
            # Formatted Database Files Structure Parameters
            "format": 1,
            "precision": 2,
            "prim_vars_wrt": "T",
            "parallel_io": "T",
            "fd_order": 2,
            'schlieren_wrt'                :'T',
            "schlieren_alpha(1)": 0.5,
            "schlieren_alpha(2)": 0.5,
            "probe_wrt": "F",
            
            # Patch 1 Background Air
            "patch_icpp(1)%geometry": 3,
            "patch_icpp(1)%x_centroid": 0.50*Lx,
            "patch_icpp(1)%y_centroid": 0.50*Ly,
            "patch_icpp(1)%length_x": Lx,
            "patch_icpp(1)%length_y": Ly,
            "patch_icpp(1)%vel(1)": 0.0,
            "patch_icpp(1)%vel(2)": 0.0,
            "patch_icpp(1)%pres": p0,
            "patch_icpp(1)%alpha_rho(1)": rho_air,
            "patch_icpp(1)%alpha_rho(2)": 0.0,
            "patch_icpp(1)%alpha(1)": 1.0,
            "patch_icpp(1)%alpha(2)": 0.0,
            
            # Patch 2 Shocked Air
            "patch_icpp(2)%geometry": 3,
            "patch_icpp(2)%x_centroid": (Lx + xs) / 2.0,
            "patch_icpp(2)%y_centroid": 0.5*Ly,
            "patch_icpp(2)%length_x": Lx - xs,
            "patch_icpp(2)%length_y": Ly,
            "patch_icpp(2)%alter_patch(1)": "T",
            "patch_icpp(2)%vel(1)": vs,
            "patch_icpp(2)%vel(2)": 0.0,
            "patch_icpp(2)%pres": ps,
            "patch_icpp(2)%alpha_rho(1)": rho_s,
            "patch_icpp(2)%alpha_rho(2)": 0.0,
            "patch_icpp(2)%alpha(1)": 1.0,
            "patch_icpp(2)%alpha(2)": 0.0,

            # Patch 3 Circular Helium Droplet
            "patch_icpp(3)%geometry": 2,
            "patch_icpp(3)%x_centroid": 13.8e-2,
            "patch_icpp(3)%y_centroid": 0.5*Ly,
            "patch_icpp(3)%radius": 2.5e-2,
            "patch_icpp(3)%alter_patch(1)": "T",
            "patch_icpp(3)%vel(1)": 0.0,
            "patch_icpp(3)%vel(2)": 0.0,
            "patch_icpp(3)%pres": p0,
            "patch_icpp(3)%alpha_rho(1)": 0.0,
            "patch_icpp(3)%alpha_rho(2)": rho_he,
            "patch_icpp(3)%alpha(1)": 0.0,
            "patch_icpp(3)%alpha(2)": 1.0,

            # Fluids Physical Parameters
            # air
            "fluid_pp(1)%gamma": Gamma_air,
            "fluid_pp(1)%pi_inf": 0.0,
	        "fluid_pp(1)%Re(1)" : mu_air,
            "fluid_pp(1)%Re(2)" : kappa_air,

            # Helium
            "fluid_pp(2)%gamma": Gamma_he,
            "fluid_pp(2)%pi_inf": 0.0,
	        "fluid_pp(2)%Re(1)" : mu_he,
            "fluid_pp(2)%Re(2)" : kappa_he,

        }
    )
)
