#!/usr/bin/env python3
import math, json


p0 = 101325.0




Rbar = 8314.3
small_fluid = 1.0e-13
# c02 props
gamma_co2 = 1.28
rho_co2 = 1.977
W_co2 = 44.01
Gamma_co2 = 1.0 / (gamma_co2 - 1.0)
cv_co2 = Rbar * Gamma_co2 / W_co2



# N2 props
gamma_n2 = 1.4
rho_n2 = 1.2506
W_n2 = 28.02
Gamma_n2 = 1.0 / (gamma_n2 - 1.0)
cv_n2 = Rbar * Gamma_n2 / W_n2

D_ab = 1.0e-1

Nx = 199

# Configuring case dictionary
print(
    json.dumps(
        {
            # Logistics
            "run_time_info": "T",
            # Computational Domain Parameters
            "x_domain%beg": 0.0,
            "x_domain%end": 1.0,
            "stretch_x": "F",
            "cyl_coord": "F",
            "m": Nx,
            "n": 0,
            "p": 0,
            "dt": 1.0e-5,
            "t_step_start": 0,
            "t_step_stop": 200000,
            "t_step_save": 250,
            # Simulation Algorithm Parameters
            "num_patches": 2,
            "model_eqns": 2,
            "alt_soundspeed": "T",
            "diffusion": "T",
            "num_fluids": 2,
            "mpp_lim": "F",
            "mixture_err": "F",
            "time_stepper": 3,
            "weno_order": 5,
            "weno_eps": 1.0e-16,
            "weno_Re_flux": "F",
            "weno_Dif_flux": "F",
            "weno_avg": "F",
            "mapped_weno": "T",
            "null_weights": "F",
            "mp_weno": "F",
            "riemann_solver": 2,
            "wave_speeds": 1,
            "avg_state": 2,
            "bc_x%beg": -2,
            "bc_x%end": -2,
            # Formatted Database Files Structure Parameters
            "format": 1,
            "precision": 2,
            "prim_vars_wrt": "T",
            "parallel_io": "F",
            "fd_order": 2,
            'schlieren_wrt'                :'F',
            "probe_wrt": "F",
            
            # Patch 1 C02
            "patch_icpp(1)%geometry": 1,
            "patch_icpp(1)%x_centroid": 0.5,
            "patch_icpp(1)%length_x": 1.0,
            "patch_icpp(1)%vel(1)": 0.0,
            "patch_icpp(1)%pres": p0,
            "patch_icpp(1)%alpha_rho(1)": rho_co2*(1.0 - small_fluid),
            "patch_icpp(1)%alpha_rho(2)": rho_n2*small_fluid,
            "patch_icpp(1)%alpha(1)": 1.0 - small_fluid,
            "patch_icpp(1)%alpha(2)": small_fluid,
            
            # Patch 2 N2
            "patch_icpp(2)%geometry": 1,
            "patch_icpp(2)%x_centroid": 0.75,
            "patch_icpp(2)%length_x": 0.5,
            "patch_icpp(2)%alter_patch(1)": "T",
            "patch_icpp(2)%vel(1)": 0.0,
            "patch_icpp(2)%pres": p0,
            "patch_icpp(2)%alpha_rho(1)": rho_co2*small_fluid,
            "patch_icpp(2)%alpha_rho(2)": rho_n2*(1.0 - small_fluid),
            "patch_icpp(2)%alpha(1)": small_fluid,
            "patch_icpp(2)%alpha(2)": 1.0 - small_fluid,
            
            # Fluids Physical Parameters
            # CO2
            "fluid_pp(1)%gamma": Gamma_co2,
            "fluid_pp(1)%pi_inf": 0.0,
            "fluid_pp(1)%cv": cv_co2,
            "fluid_pp(1)%D": D_ab,
            
            # N2
            "fluid_pp(2)%gamma": Gamma_n2,
            "fluid_pp(2)%pi_inf": 0.0,
            "fluid_pp(2)%cv": cv_n2,
            "fluid_pp(2)%D": D_ab,
            
        }
    )
)
