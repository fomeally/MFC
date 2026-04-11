!>
!! @file m_conduction.f90
!! @brief Contains module m_conduction

#:include 'macros.fpp'

!> @brief This module is used to compute flux terms for heat conduction
module m_conduction

    ! Dependencies =============================================================

    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters

    use m_finite_differences   !< Finite difference module

    use m_mpi_proxy            !< Message passing interface (MPI) module proxy

    ! use m_weno                 !< WENO module

    use m_helper              !< Helper functions

    ! ==========================================================================

    implicit none

        private; public :: s_initialize_conduction_module, &
s_fill_face_weights_1d_2nd_cond, &
s_compute_conduction_rhs, &
s_finalize_conduction_module

    type(int_bounds_info) :: isc1, isc2, isc3
    !$acc declare create(isc1, isc2, isc3)

    real(wp), allocatable, dimension(:) :: ks
    !$acc declare create(ks)

    real(wp), allocatable, dimension(:, :, :) :: w_interp2_cond, w_grad2_cond
    !$acc declare create(w_interp2_cond, w_grad2_cond)

contains

    subroutine s_initialize_conduction_module

        integer :: i !< generic loop iterator

        @:ALLOCATE(ks(1:num_fluids))
        !$acc loop seq
        do i = 1, num_fluids
            ks(i) = fluid_pp(i)%k
        end do
        !$acc update device(ks)





        @:ALLOCATE(w_interp2_cond(0:1, -1:max(m, n, p), 1:num_dims))
        @:ALLOCATE(w_grad2_cond(0:1, -1:max(m, n, p), 1:num_dims))

        call s_fill_face_weights_1d_2nd_cond(-1, m, 1)
        if (n > 0) then
            call s_fill_face_weights_1d_2nd_cond(-1, n, 2)

            if (p > 0) then
                call s_fill_face_weights_1d_2nd_cond(-1, p, 3)
            end if
        end if
        !$acc update device(w_interp2_cond, w_grad2_cond)

    end subroutine s_initialize_conduction_module

    subroutine s_fill_face_weights_1d_2nd_cond(ifbeg, ifend, idir)

        integer,  intent(in) :: ifbeg, ifend, idir

        integer :: i
        real(wp) :: xf, dx_loc

        select case (idir)
        case (1)
            do i = ifbeg, ifend

                xf     = x_cb(i)
                dx_loc = x_cc(i+1) - x_cc(i)

                w_interp2_cond( 0, i, idir) = (x_cc(i+1) - xf) / dx_loc
                w_interp2_cond( 1, i, idir) = (xf - x_cc(i))   / dx_loc

                w_grad2_cond( 0, i, idir) = -1._wp / dx_loc
                w_grad2_cond( 1, i, idir) =  1._wp / dx_loc

            end do
        
        case (2)
            do i = ifbeg, ifend

                xf     = y_cb(i)
                dx_loc = y_cc(i+1) - y_cc(i)

                w_interp2_cond( 0, i, idir) = (y_cc(i+1) - xf) / dx_loc
                w_interp2_cond( 1, i, idir) = (xf - y_cc(i))   / dx_loc

                w_grad2_cond( 0, i, idir) = -1._wp / dx_loc
                w_grad2_cond( 1, i, idir) =  1._wp / dx_loc

            end do

        case (3)
            do i = ifbeg, ifend

                xf     = z_cb(i)
                dx_loc = z_cc(i+1) - z_cc(i)

                w_interp2_cond( 0, i, idir) = (z_cc(i+1) - xf) / dx_loc
                w_interp2_cond( 1, i, idir) = (xf - z_cc(i))   / dx_loc

                w_grad2_cond( 0, i, idir) = -1._wp / dx_loc
                w_grad2_cond( 1, i, idir) =  1._wp / dx_loc

            end do

        end select

    end subroutine s_fill_face_weights_1d_2nd_cond

    subroutine s_compute_conduction_rhs(idir, q_src_n, q_prim_vf, irx, iry, irz)

        integer, intent(in) :: idir
        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf
        type(scalar_field), dimension(sys_size), intent(inout) :: q_src_n
        type(int_bounds_info), intent(in) :: irx, iry, irz

        integer :: i, k, l, q, r !< Loop variables
        real(wp) :: R_univ
        real(wp) :: grid_spacing
        real(wp) :: rho_L, rho_R, rhog_L, rhog_R
        real(wp) :: rho_liq_L, rho_liq_R
        real(wp) :: alpha_m_L, alpha_m_R
        real(wp) :: P_L, P_R
        real(wp) :: W_L, W_R
        real(wp) :: ds_L, ds_R
        real(wp) :: kg_f
        real(wp) :: Tg_L, Tg_R, dTg_ds_f
        integer, dimension(3) :: offsets

        real(wp) :: alpharho_L(Dif_size), alpharho_R(Dif_size)
        real(wp) :: Y_L(Dif_size), Y_R(Dif_size)
        real(wp) :: S_L, S_R
        real(wp) :: alpha_i_f(num_fluids)
        real(wp) :: T_L(num_fluids), T_R(num_fluids)
        real(wp) :: rho_i_L(num_fluids), rho_i_R(num_fluids)
        real(wp) :: dT_ds_f(num_fluids)

        R_univ = 8314.462618_wp
        isc1 = irx; isc2 = iry; isc3 = irz

        ! Set offsets based on direction using array indexing
        offsets = 0
        offsets(idir) = 1
        
        ! Finite Volume with q_src_n Approach
        ! #########################################################################
        ! #########################################################################
        do q = isc3%beg, isc3%end
            do l = isc2%beg, isc2%end
                do k = isc1%beg, isc1%end

                    q_src_n(E_idx)%sf(k, l, q) = 0._wp

                    ! Calculate grid spacing using direction-based indexing
                    select case (idir)
                    case (1)
                        r = k
                        grid_spacing = x_cc(k+1) - x_cc(k)
                        ds_R = x_cc(k+1) - x_cb(k)
                        ds_L = x_cb(k) - x_cc(k)
                    case (2)
                        r = l
                        grid_spacing = y_cc(l+1) - y_cc(l)
                        ds_R = y_cc(l+1) - y_cb(l)
                        ds_L = y_cb(l) - y_cc(l)
                    case (3)
                        r = q
                        grid_spacing = z_cc(q+1) - z_cc(q)
                        ds_R = z_cc(q+1) - z_cb(q)
                        ds_L = z_cb(q) - z_cc(q)
                    end select

                    print *, "ds_L, ds_R, grid_spacing: ", ds_L, ds_R, grid_spacing


                    P_L = q_prim_vf(E_idx)%sf(k, l, q)
                    P_R = q_prim_vf(E_idx)%sf(k + offsets(1), l + offsets(2), q + offsets(3))

                    if ( abs(P_R - P_L) / max(P_L, P_R, 1.0e-8_wp) > 0.05_wp ) cycle !crude shock sensor (dont calculate conduction across shocks)

                    do i = 1, num_fluids
                        alpha_i_f(i) = 0.0_wp
                        dT_ds_f(i)   = 0.0_wp
                        T_L(i) = 0.0_wp
                        T_R(i) = 0.0_wp
                    end do

                    ! calculate temperatures depending on if there is a mixture gas or not
                    if (diffusion) then
                        alpha_m_L = q_prim_vf(advg_idx)%sf(k, l, q)
                        alpha_m_R = q_prim_vf(advg_idx)%sf(k + offsets(1), l + offsets(2), q + offsets(3))
                        if (alpha_m_L > small_num_dif .and. alpha_m_R > small_num_dif) then !(conduction at face for mixture gas)
                        
                            rho_L = 0.0_wp
                            rho_R = 0.0_wp
                            do i = 1, Dif_size
                                alpharho_L(i) = q_prim_vf(Dif_idx(i))%sf(k, l, q)
                                alpharho_R(i) = q_prim_vf(Dif_idx(i))%sf(k + offsets(1), l + offsets(2), q + offsets(3))
                            end do

                            do i = 1, Dif_size
                                rho_L = rho_L + alpharho_L(i)
                                rho_R = rho_R + alpharho_R(i)
                            end do

                            do i = 1, Dif_size
                                Y_L(i) = alpharho_L(i) / rho_L
                                Y_R(i) = alpharho_R(i) / rho_R
                            end do

                            rhog_L = rho_L / alpha_m_L
                            rhog_R = rho_R / alpha_m_R

                            W_L = 0._wp
                            W_R = 0._wp
                            do i = 1, Dif_size  
                                W_L = W_L + Y_L(i)/fluid_pp(Dif_idx(i))%W
                                W_R = W_R + Y_R(i)/fluid_pp(Dif_idx(i))%W        
                            end do

                            W_L = 1._wp / W_L
                            W_R = 1._wp / W_R

                            Tg_L = P_L * W_L / (rhog_L * R_univ)
                            Tg_R = P_R * W_R / (rhog_R * R_univ)
                            dTg_ds_f = (Tg_R - Tg_L) / grid_spacing

                            S_L = 0.0_wp
                            S_R = 0.0_wp   
                            do i = 1, Dif_size
                                S_L = S_L + (alpha_m_L * Y_L(i) * W_L / fluid_pp(Dif_idx(i))%W) * ks(i)
                                S_R = S_R + (alpha_m_R * Y_R(i) * W_R / fluid_pp(Dif_idx(i))%W) * ks(i)
                            end do
                            kg_f = (ds_L + ds_R) * S_L * S_R / (ds_L*S_R + ds_R*S_L)
                            q_src_n(E_idx)%sf(k, l, q) = q_src_n(E_idx)%sf(k, l, q) - kg_f * dTg_ds_f 

                        end if

                        if (num_fluids > Dif_size) then !there is liquid, use SG EOS
                            if (q_prim_vf(E_idx + liq_idx)%sf(k, l, q) > small_num_dif .and. q_prim_vf(E_idx + liq_idx)%sf(k + offsets(1), l + offsets(2), q + offsets(3)) > small_num_dif) then !(conduction at face for liquid)

                                alpha_i_f(liq_idx) = (ds_L + ds_R) * q_prim_vf(E_idx + liq_idx)%sf(k, l, q) * q_prim_vf(E_idx + liq_idx)%sf(k + offsets(1), l + offsets(2), q + offsets(3)) / &
                                    ( ds_L*q_prim_vf(E_idx + liq_idx)%sf(k + offsets(1), l + offsets(2), q + offsets(3)) + ds_R*q_prim_vf(E_idx + liq_idx)%sf(k, l, q) )
                                rho_liq_L = q_prim_vf(liq_idx)%sf(k, l, q) / q_prim_vf(E_idx + liq_idx)%sf(k, l, q)
                                rho_liq_R = q_prim_vf(liq_idx)%sf(k + offsets(1), l + offsets(2), q + offsets(3)) / q_prim_vf(E_idx + liq_idx)%sf(k + offsets(1), l + offsets(2), q + offsets(3))
                                T_L(liq_idx) = gammas(liq_idx) / (cvs(liq_idx) * rho_liq_L) * (P_L + pi_infs(liq_idx) / (gammas(liq_idx) + 1.0_wp))
                                T_R(liq_idx) = gammas(liq_idx) / (cvs(liq_idx) * rho_liq_R) * (P_R + pi_infs(liq_idx) / (gammas(liq_idx) + 1.0_wp))
                                dT_ds_f(liq_idx) = (T_R(liq_idx) - T_L(liq_idx)) / grid_spacing

                                q_src_n(E_idx)%sf(k, l, q) = q_src_n(E_idx)%sf(k, l, q) - ks(liq_idx)*alpha_i_f(liq_idx)*dT_ds_f(liq_idx)
                            end if
                        end if

                    else ! no diffusion, if present on both sides, compute temp of each component with SG EOS

                        do i = 1, num_fluids
                            if (q_prim_vf(E_idx + i)%sf(k, l, q) > small_num_dif .and. q_prim_vf(E_idx + i)%sf(k + offsets(1), l + offsets(2), q + offsets(3)) > small_num_dif) then
                                
                                rho_i_R(i) = q_prim_vf(i)%sf(k + offsets(1), l + offsets(2), q + offsets(3)) / q_prim_vf(i + E_idx)%sf(k + offsets(1), l + offsets(2), q + offsets(3))
                                rho_i_L(i) = q_prim_vf(i)%sf(k, l, q) / q_prim_vf(i + E_idx)%sf(k, l, q)
                                T_R(i) = gammas(i) / (cvs(i) * rho_i_R(i)) * (P_R + pi_infs(i) / (gammas(i) + 1.0_wp))
                                T_L(i) = gammas(i) / (cvs(i) * rho_i_L(i)) * (P_L + pi_infs(i) / (gammas(i) + 1.0_wp))
                                dT_ds_f(i) = (T_R(i) - T_L(i)) / grid_spacing
                                alpha_i_f(i) = (ds_L + ds_R) * q_prim_vf(E_idx + i)%sf(k, l, q) * q_prim_vf(E_idx + i)%sf(k + offsets(1), l + offsets(2), q + offsets(3)) / &
                                    ( ds_L*q_prim_vf(E_idx + i)%sf(k + offsets(1), l + offsets(2), q + offsets(3)) + ds_R*q_prim_vf(E_idx + i)%sf(k, l, q) )

                                q_src_n(E_idx)%sf(k, l, q) = q_src_n(E_idx)%sf(k, l, q) - ks(i)*alpha_i_f(i)*dT_ds_f(i)
                            end if

                        end do

                    end if
                        
                end do
            end do
        end do

        ! #########################################################################
        ! #########################################################################
            

        
    end subroutine s_compute_conduction_rhs

    subroutine s_finalize_conduction_module

        @:DEALLOCATE(ks)
        @:DEALLOCATE(w_interp2_cond, w_grad2_cond)

    end subroutine s_finalize_conduction_module

end module m_conduction