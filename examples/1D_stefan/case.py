#!/usr/bin/env python3
import math, json


p0 = 101325.0
T0 = 298.0

#water props
# density - kg/m3
rho0w = 1000.0
# gama
gamw = 4.4
# pi infty - Pa
piw = 6.0e08
cvw = (p0 + piw) / (gamw - 1.0) / rho0w / T0
cpw = gamw * cvw
# qv
qvwl = -1167000
# qv'
qvpwl = 0.0e0
lambda_wl = 0.60

Rbar = 8314.462618
small_fluid = 1.0e-3


# water vapor properties
## Vapor water

Lv = 2.5e6
# Lv = 0.0

W_wv = 18.015
gamma_wv = 1.33
rho_wv = p0*W_wv / (Rbar * T0)
Gamma_wv = 1.0 / (gamma_wv - 1.0)
cv_wv = Rbar * Gamma_wv / W_wv
cp_wv = gamma_wv*cv_wv
qvwv = qvwl + Lv + (cpw - cp_wv)*T0
# qv'
# qvpwv = -21884.3  # J/(kg K)
qvpwv = -24349.7
# qvpwv = -32752.0
lambda_wv = 0.018

# air props
gamma_air = 1.4
W_air = 28.97
rho_air = p0*W_air / (Rbar * T0)
Gamma_air = 1.0 / (gamma_air - 1.0)
cv_air = Rbar * Gamma_air / W_air
cp_air = gamma_air*cv_air
lambda_air = 0.026

# # c02 props
# gamma_co2 = 1.28
# W_co2 = 44.01
# rho_co2 = p0*W_co2 / (Rbar * T0)
# Gamma_co2 = 1.0 / (gamma_co2 - 1.0)
# cv_co2 = Rbar * Gamma_co2 / W_co2
# cp_co2 = gamma_co2*cv_co2



# # N2 props
# gamma_n2 = 1.4
# #rho_n2 = 1.2506
# W_n2 = 28.02
# rho_n2 = p0*W_n2 / (Rbar * T0)
# Gamma_n2 = 1.0 / (gamma_n2 - 1.0)
# cv_n2 = Rbar * Gamma_n2 / W_n2
# cp_n2 = gamma_n2*cv_n2

# # H2 props
# gamma_h2 = 1.405
# W_h2 = 2.02
# rho_h2 = p0*W_h2 / (Rbar * T0)
# Gamma_h2 = 1.0 / (gamma_h2 - 1.0)
# cv_h2 = Rbar * Gamma_h2 / W_h2
# cp_h2 = gamma_h2*cv_h2

# Binary diffusion coefficients
D11 = 0.0e0
D12 = 2.5e-5
D13 = 0.674e-2
D21 = D12
D22 = 0.0e0
D23 = 1.0e-1
D31 = D13
D32 = D23
D33 = 0.0e0

Lx = 10.0e-6

Nx = 199


# Configuring case dictionary
print(
    json.dumps(
        {
            # Logistics
            "run_time_info": "T",
            # Computational Domain Parameters
            "x_domain%beg": 0.0,
            "x_domain%end": Lx,
            "stretch_x": "F",
            "cyl_coord": "F",
            "m": int(Nx),
            "n": 0,
            "p": 0,
            "dt": 0.5e-11,
            "t_step_start": 0,
            "t_step_stop": 400000,
            "t_step_save": 8000,
            # "t_step_save":  800,
            # "t_step_stop": 326000,
            # "t_step_save": 3260,
            "t_step_print": 4000,
            # Simulation Algorithm Parameters
            "num_patches": 3,
            "model_eqns": 2,
            "alt_soundspeed": "T",
            "diffusion": "T",
            "conduction": "F",
            "num_fluids": 3,
            "mpp_lim": "T",
            "mixture_err": "T",
            "relax": "T",
            "relax_model": 6,
            "palpha_eps": 1.0e-6,
            "ptgalpha_eps": 1.0e-2,
            "time_stepper": 3,
            "weno_order": 5,
            "weno_eps": 1.0e-16,
            "weno_Re_flux": "F",
            "small_num_dif": 1.0e-8,
            "weno_avg": "F",
            "mapped_weno": "T",
            "null_weights": "F",
            "mp_weno": "F",
            "riemann_solver": 2,
            "wave_speeds": 1,
            "avg_state": 2,
            "bc_x%beg": -2,
            "bc_x%end": -6,
            # Formatted Database Files Structure Parameters
            "format": 1,
            "precision": 2,
            "prim_vars_wrt": "T",
            "parallel_io": "F",
            "fd_order": 2,
            'schlieren_wrt'                :'F',
            "probe_wrt": "F",
            
            # Patch 1 Background Air
            "patch_icpp(1)%geometry": 1,
            "patch_icpp(1)%x_centroid": 0.50*Lx,
            "patch_icpp(1)%length_x": Lx,
            "patch_icpp(1)%vel(1)": 0.0,
            "patch_icpp(1)%pres": p0,
            "patch_icpp(1)%alpha_rho(1)": 0.0,
            "patch_icpp(1)%alpha_rho(2)": 0.0,
            "patch_icpp(1)%alpha_rho(3)": rho_air,
            "patch_icpp(1)%alpha(1)": 0.0,
            "patch_icpp(1)%alpha(2)": 0.0,
            "patch_icpp(1)%alpha(3)": 1.0,

            # Patch 2 Water 
            "patch_icpp(2)%geometry": 1,
            "patch_icpp(2)%x_centroid": 0.10*Lx,
            "patch_icpp(2)%length_x": 0.2*Lx,
            "patch_icpp(2)%alter_patch(1)": "T",
            "patch_icpp(2)%vel(1)": 0.0,
            "patch_icpp(2)%pres": p0,
            "patch_icpp(2)%alpha_rho(1)": rho0w,
            "patch_icpp(2)%alpha_rho(2)": 0.0,
            "patch_icpp(2)%alpha_rho(3)": 0.0,
            "patch_icpp(2)%alpha(1)": 1.0,
            "patch_icpp(2)%alpha(2)": 0.0,
            "patch_icpp(2)%alpha(3)": 0.0,

            # Patch 3 Smeared Region
            "patch_icpp(3)%geometry": 23,
            "patch_icpp(3)%x_centroid": 0.20*Lx,
            "patch_icpp(3)%length_x": 0.4*Lx,
            "patch_icpp(3)%alter_patch(1)": "T",
            "patch_icpp(3)%alter_patch(2)": "T",
            "patch_icpp(3)%vel(1)": 0.0,
            "patch_icpp(3)%pres": p0,
            "patch_icpp(3)%alpha_rho(1)": 0.0,
            "patch_icpp(3)%alpha_rho(2)": small_fluid * rho_wv,
            "patch_icpp(3)%alpha_rho(3)": (1.0 - small_fluid) * rho_air,
            "patch_icpp(3)%alpha(1)": 0.0,
            "patch_icpp(3)%alpha(2)": small_fluid,
            "patch_icpp(3)%alpha(3)": 1.0 - small_fluid,

            # Fluids Physical Parameters
            # Liquid Water
            "fluid_pp(1)%gamma": 1.0e0 / (gamw - 1.0e0),
            "fluid_pp(1)%pi_inf": gamw*piw / (gamw - 1.0e0),
	        "fluid_pp(1)%cv": cvw,
            "fluid_pp(1)%qv": qvwl,
            "fluid_pp(1)%qvp": qvpwl,
            # "fluid_pp(1)%k": lambda_wl,
            "fluid_pp(1)%gas_mixture" : "F",

            # Water Vapor
            "fluid_pp(2)%gamma": Gamma_wv,
            "fluid_pp(2)%pi_inf": 0.0,
	        "fluid_pp(2)%W": W_wv,
            "fluid_pp(2)%cv": cv_wv,
	        "fluid_pp(2)%cp": cp_wv,
            "fluid_pp(2)%qv": qvwv,
            "fluid_pp(2)%qvp": qvpwv,
            "fluid_pp(2)%D(1)": D22,
            "fluid_pp(2)%D(2)": D12,
            # "fluid_pp(2)%k": lambda_wv,
            "fluid_pp(2)%gas_mixture" : "T",

            # Air
            "fluid_pp(3)%gamma": Gamma_air,
            "fluid_pp(3)%pi_inf": 0.0,
            "fluid_pp(3)%W": W_air,
            "fluid_pp(3)%cp": cp_air,
            "fluid_pp(3)%cv": cv_air,
            "fluid_pp(3)%D(1)": D21,
            "fluid_pp(3)%D(2)": D22,
            # "fluid_pp(3)%k": lambda_air,
            "fluid_pp(3)%gas_mixture" : "T",

        }
    )
)
