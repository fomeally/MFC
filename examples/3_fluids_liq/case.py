#!/usr/bin/env python3
import math, json


p0 = 101325.0
T0 = 298.0

#water props
# density - kg/m3
rho0w = 1000
# gama
gamw = 6.12
# pi infty - Pa
piw = 3.43e08

Rbar = 8314.462618
small_fluid = 1.0e-13
# c02 props
gamma_co2 = 1.28
W_co2 = 44.01
rho_co2 = p0*W_co2 / (Rbar * T0)
Gamma_co2 = 1.0 / (gamma_co2 - 1.0)
cv_co2 = Rbar * Gamma_co2 / W_co2
cp_co2 = gamma_co2*cv_co2
T0_co2 = 0.0
h0_co2 = 0.0



# N2 props
gamma_n2 = 1.4
#rho_n2 = 1.2506
W_n2 = 28.02
rho_n2 = p0*W_n2 / (Rbar * T0)
Gamma_n2 = 1.0 / (gamma_n2 - 1.0)
cv_n2 = Rbar * Gamma_n2 / W_n2
cp_n2 = gamma_n2*cv_n2
T0_n2 = 0.0
h0_n2 = 0.0

# H2 props
gamma_h2 = 1.405
W_h2 = 2.02
rho_h2 = p0*W_h2 / (Rbar * T0)
Gamma_h2 = 1.0 / (gamma_h2 - 1.0)
cv_h2 = Rbar * Gamma_h2 / W_h2
cp_h2 = gamma_h2*cv_h2
T0_h2 = 0.0
h0_h2 = 0.0

# Binary diffusion coefficients
D11 = 0.0e0
D12 = 1.0e-1
D13 = 0.674e-2
D21 = D12
D22 = 0.0e0
D23 = 0.74e-2
D31 = D13
D32 = D23
D33 = 0.0e0

Lx = 1.0

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
            "m": Nx,
            "n": 0,
            "p": 0,
            "dt": 1.0e-5,
            "t_step_start": 0,
            # "t_step_stop": 10000,
            # "t_step_save": 20,
            "t_step_stop": 500000,
            "t_step_save": 5000,
            "t_step_print": 5000,
            # "t_step_stop": 10000000,
            # "t_step_save": 50000,
            # Simulation Algorithm Parameters
            "num_patches": 2,
            "model_eqns": 2,
            "alt_soundspeed": "T",
            "diffusion": "T",
            "num_fluids": 3,
            "mpp_lim": "F",
            "mixture_err": "F",
            "time_stepper": 3,
            "weno_order": 5,
            "weno_eps": 1.0e-16,
            "weno_Re_flux": "F",
            "weno_Dif_flux": "F",
            "Dif_fv" : "T",
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
            "schlieren_wrt": "F",
            "probe_wrt": "F",
            
            # Patch 1 N2
            "patch_icpp(1)%geometry": 1,
            "patch_icpp(1)%x_centroid": 0.50*Lx,
            "patch_icpp(1)%length_x": Lx,
            "patch_icpp(1)%vel(1)": 0.0,
            "patch_icpp(1)%pres": p0,
            "patch_icpp(1)%alpha_rho(1)": rho_n2,
            "patch_icpp(1)%alpha_rho(2)": 0.0,
            "patch_icpp(1)%alpha_rho(3)": 0.0,
            "patch_icpp(1)%alpha(1)": 1.0,
            "patch_icpp(1)%alpha(2)": 0.0,
            "patch_icpp(1)%alpha(3)": 0.0,
            
            # Patch 2 CO2
            "patch_icpp(2)%geometry": 1,
            "patch_icpp(2)%x_centroid": 0.75*Lx,
            "patch_icpp(2)%length_x": 0.50*Lx,
            "patch_icpp(2)%alter_patch(1)": "T",
            "patch_icpp(2)%vel(1)": 0.0,
            "patch_icpp(2)%pres": p0,
            "patch_icpp(2)%alpha_rho(1)": 0.0,
            "patch_icpp(2)%alpha_rho(2)": rho_co2,
            "patch_icpp(2)%alpha_rho(3)": 0.0,
            "patch_icpp(2)%alpha(1)": 0.0,
            "patch_icpp(2)%alpha(2)": 1.0,
            "patch_icpp(2)%alpha(3)": 0.0,

            # Fluids Physical Parameters
            # N2
            "fluid_pp(1)%gamma": Gamma_n2,
            "fluid_pp(1)%pi_inf": 0.0,
	        "fluid_pp(1)%W": W_n2,
	        "fluid_pp(1)%cp": cp_n2,
            "fluid_pp(1)%h0": h0_n2,
            "fluid_pp(1)%T0": T0_n2,
            "fluid_pp(1)%D(1)": D11,
            "fluid_pp(1)%D(2)": D12,
            "fluid_pp(1)%gas_mixture" : "T",

            # CO2
            "fluid_pp(2)%gamma": Gamma_co2,
            "fluid_pp(2)%pi_inf": 0.0,
	        "fluid_pp(2)%W": W_co2,
	        "fluid_pp(2)%cp": cp_co2,
            "fluid_pp(2)%h0": h0_co2,
            "fluid_pp(2)%T0": T0_co2,
            "fluid_pp(2)%D(1)": D21,
            "fluid_pp(2)%D(2)": D22,
            "fluid_pp(2)%gas_mixture" : "T",
    
            # Water
            "fluid_pp(3)%gamma": gamw,
            "fluid_pp(3)%pi_inf": piw,
            "fluid_pp(3)%gas_mixture" : "F",
        }
    )
)
