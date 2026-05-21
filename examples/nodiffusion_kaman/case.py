#!/usr/bin/env python3
import math, json


p0 = 0.956e5
M = 1.11
M2 = M*M
Rbar = 8314.462618
#Interface locations
shock_interface = 74.0e-2
material_interface = 72.0e-2


#Air props
gamma_air = 1.276
W_air = 34.76
rho_air = 1.351
T0 = p0*W_air / (rho_air * Rbar)
Gamma_air = 1.0 / (gamma_air - 1.0)
cv_air = Rbar * Gamma_air / W_air
cp_air = gamma_air*cv_air
T0_air = 0.0
h0_air = 0.0
mu_air = 1.0 / (2.0e-5)     #mixture air and acetone 75/2 by volume
kappa_air = 1.0 / (0.7e-5)
k_air = 2.2e-2

#shock air props
rho_shock = rho_air*(gamma_air + 1.0)*M2 / ((gamma_air - 1.0)*M2 + 2.0)
ps =  p0*(1.0 + 2.0*gamma_air/(gamma_air + 1.0)*(M2 - 1.0))
c_air = math.sqrt(gamma_air * p0 / rho_air)
vs = 2.0*c_air*(M2 - 1.0) / ((gamma_air + 1.0)*M)

# SF6 props
gamma_sf6 = 1.093
W_sf6 = 146.06
rho_sf6 = p0*W_sf6 / (Rbar * T0)
Gamma_sf6 = 1.0 / (gamma_sf6 - 1.0)
cv_sf6 = Rbar * Gamma_sf6 / W_sf6
cp_sf6 = gamma_sf6*cv_sf6
T0_sf6 = 0.0
h0_sf6 = 0.0
mu_sf6 = 1.0 / (1.5e-5)
kappa_sf6 = 1.0 / (1.0e-5)
k_sf6 = 0.013


# Binary diffusion coefficients
D11 = 0.0e0
D12 = 1.0e-5
D21 = D12
D22 = 0.0e0

Lx = 8.9e-2
Ly = 75.0e-2


Nx = 384
Ny = 3236


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
            "dt": 2.5e-7,
            "t_step_start": 0,
            "t_step_stop": 1,
            "t_step_save": 1,
            # "t_step_stop": 36000,
            # "t_step_save": 200,
            # "t_step_print": 360,
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
            #"Dif_fv" : "F",
	    #    "dif_order": 2,
            "small_num_dif": 1.0e-8,
            "weno_avg": "F",
            "mapped_weno": "T",
            "null_weights": "F",
            "mp_weno": "F",
            "riemann_solver": 2,
            "wave_speeds": 1,
            "avg_state": 2,
            "bc_x%beg": -17,
            "bc_x%end": -17, #symmetry without flipping normal velocity
            "bc_y%beg": -7, #subsonic inflow
            "bc_y%end": -2,  #symmetry with flipping normal velocity
       	    "viscous": "T",
            "conduction": "T",
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
            
            # Patch 1 Background air
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
            
            # Patch 2 Shocked air 
            "patch_icpp(2)%geometry": 3,
            "patch_icpp(2)%x_centroid": 0.5*Lx,
            "patch_icpp(2)%y_centroid": (Ly + shock_interface) / 2.0,
            "patch_icpp(2)%length_x": Lx,
            "patch_icpp(2)%length_y": Ly - shock_interface,
            "patch_icpp(2)%alter_patch(1)": "T",
            "patch_icpp(2)%vel(1)": 0.0,
            "patch_icpp(2)%vel(2)": -vs,
            "patch_icpp(2)%pres": ps,
            "patch_icpp(2)%alpha_rho(1)": rho_shock,
            "patch_icpp(2)%alpha_rho(2)": 0.0,
            "patch_icpp(2)%alpha(1)": 1.0,
            "patch_icpp(2)%alpha(2)": 0.0,

            # Patch 3 Sinusoidal SF6 interface
            "patch_icpp(3)%geometry": 26,
            "patch_icpp(3)%smoothen": "T",
            "patch_icpp(3)%smooth_patch_id": 1,
            "patch_icpp(3)%smooth_coeff": 0.5e-2,
            "patch_icpp(3)%x_centroid": 0.5*Lx,
            "patch_icpp(3)%y_centroid": 0.5*material_interface,
            "patch_icpp(3)%length_x": Lx,
            "patch_icpp(3)%length_y": material_interface,
            "patch_icpp(3)%alter_patch(1)": "T",
            "patch_icpp(3)%vel(1)": 0.0,
            "patch_icpp(3)%vel(2)": 0.0,
            "patch_icpp(3)%pres": p0,
            "patch_icpp(3)%alpha_rho(1)": 0.0,
            "patch_icpp(3)%alpha_rho(2)": rho_sf6,
            "patch_icpp(3)%alpha(1)": 0.0,
            "patch_icpp(3)%alpha(2)": 1.0,

            # Fluids Physical Parameters
            # air
            "fluid_pp(1)%gamma": Gamma_air,
            "fluid_pp(1)%pi_inf": 0.0,
	        "fluid_pp(1)%Re(1)" : mu_air,
            "fluid_pp(1)%Re(2)" : kappa_air,
	        "fluid_pp(1)%W": W_air,
	        "fluid_pp(1)%cp": cp_air,
            "fluid_pp(1)%cv": cv_air,
            #"fluid_pp(1)%h0": h0_air,
            #"fluid_pp(1)%T0": T0_air,
            #"fluid_pp(1)%D(1)": D11,
            #"fluid_pp(1)%D(2)": D12,
            #"fluid_pp(1)%gas_mixture" : "T",
	        "fluid_pp(1)%k": k_air,

            # SF6
            "fluid_pp(2)%gamma": Gamma_sf6,
            "fluid_pp(2)%pi_inf": 0.0,
	        "fluid_pp(2)%Re(1)" : mu_sf6,
            "fluid_pp(2)%Re(2)" : kappa_sf6,
	        "fluid_pp(2)%W": W_sf6,
	        "fluid_pp(2)%cp": cp_sf6,
            "fluid_pp(2)%cv": cv_sf6,
            #"fluid_pp(2)%h0": h0_sf6,
            #"fluid_pp(2)%T0": T0_sf6,
            #"fluid_pp(2)%D(1)": D21,
            #"fluid_pp(2)%D(2)": D22,
            #"fluid_pp(2)%gas_mixture" : "T",
            "fluid_pp(2)%k": k_sf6,

        }
    )
)
