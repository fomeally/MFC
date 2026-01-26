!>
!! @file m_diffusion.f90
!! @brief Contains module m_diffusion

#:include 'macros.fpp'

!> @brief This module is used to compute flux terms for binary diffusion
module m_diffusion

    ! Dependencies =============================================================

    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters

    use m_finite_differences   !< Finite difference module

    use m_mpi_proxy            !< Message passing interface (MPI) module proxy

    ! use m_weno                 !< WENO module

    use m_helper              !< Helper functions

    ! ==========================================================================

    implicit none

    private; public :: s_initialize_diffusion_module, &
s_compute_sum_alpha_g, &
s_compute_diffusion_rhs, &
s_correct_volume_fractions, &
s_correct_riemann_volume_fractions, &
s_calculate_multicomponent_diffusion_flux, &
s_finalize_diffusion_module

    real(wp), allocatable, dimension(:, :) :: fd_coeff_x_d
    real(wp), allocatable, dimension(:, :) :: fd_coeff_y_d
    real(wp), allocatable, dimension(:, :) :: fd_coeff_z_d
    !$acc declare create(fd_coeff_x_d,fd_coeff_y_d,fd_coeff_z_d)

    type(int_bounds_info) :: isd1, isd2, isd3
    !$acc declare create(isd1, isd2, isd3)

    real(wp), allocatable, dimension(:, :) :: Ds
    !$acc declare create(Ds)

    real(wp), allocatable, dimension(:) :: Ws
    !$acc declare create(Ws)

    real(wp), allocatable, dimension(:) :: cps
    !$acc declare create(cp)

    real(wp), allocatable, dimension(:) :: T0s
    !$acc declare create(T0s)

    real(wp), allocatable, dimension(:) :: h0s
    !$acc declare create(h0s)

    real(wp), allocatable, dimension(:, :, :, :) :: dj_dx, dj_dy, dj_dz, djh_dx, djh_dy, djh_dz, dY_dx, dY_dy, dY_dz, alpha_K_dif, alpharho_K_dif, Y_dif, h_dif
    !$acc declare create(dj_dx, dj_dy, dj_dz, djh_dx, djh_dy, djh_dz, dY_dx, dY_dy, dY_dz, alpha_K_dif, alpharho_K_dif, Y_dif, h_dif)

    real(wp), allocatable, dimension(:, :, :) :: rho_dif, alpha_dif, dvel_dx, dvel_dy, dvel_dz, denom, rhogcg2, rho1c12, kdivu, W_dif, T_dif
    !$acc declare create(rho_dif, alpha_dif, dvel_dx, dvel_dy, dvel_dz, denom, rhogcg2, rho1c12, kdivu, W_dif, T_dif)

contains

    subroutine s_initialize_diffusion_module

        integer :: i, j !< generic loop iterators
        integer :: m_end, n_end, p_end, m_end_Y, n_end_Y, p_end_Y
        type(int_bounds_info) :: offset_s(1:3)

        offset_s(1)%beg = fd_number; offset_s(1)%end = fd_number
        if (n > 0) then
            offset_s(2)%beg = fd_number 
            offset_s(2)%end = fd_number
        else
            offset_s(2)%beg = 0
            offset_s(2)%end = 0
        end if
        if (p > 0) then
            offset_s(3)%beg = fd_number
            offset_s(3)%end = fd_number
        else
            offset_s(3)%beg = 0
            offset_s(3)%end = 0
        end if
        m_end = m + fd_number; n_end = n + fd_number; p_end = p + fd_number
        m_end_Y = m + 2*fd_number; n_end_Y = n + 2*fd_number; p_end_Y = p + 2*fd_number;

        @:ALLOCATE(Ds(1:Dif_size, 1:Dif_size))
        !$acc loop seq
        do i = 1, Dif_size
            do j = 1, Dif_size
                Ds(i, j) = fluid_pp(Dif_idx(i))%D(j)
            end do
        end do
        !$acc update device(Ds)

        @:ALLOCATE(Ws(1:Dif_size))
        @:ALLOCATE(cps(1:Dif_size))
        @:ALLOCATE(T0s(1:Dif_size))
        @:ALLOCATE(h0s(1:Dif_size))
        !$acc loop seq
        do i = 1, Dif_size
            Ws(i) = fluid_pp(Dif_idx(i))%W
            cps(i) = fluid_pp(Dif_idx(i))%cp
            T0s(i) = fluid_pp(Dif_idx(i))%T0
            h0s(i) = fluid_pp(Dif_idx(i))%h0
        end do
        !$acc update device(Ws, cps, T0s, h0s)
        ! Allocate arrays

        @:ALLOCATE(Y_dif(-2*fd_number:m_end_Y, -2*fd_number:n_end_Y, -2*fd_number:p_end_Y, 1:Dif_size))
        @:ALLOCATE(h_dif(-2*fd_number:m_end_Y, -2*fd_number:n_end_Y, -2*fd_number:p_end_Y, 1:Dif_size))
        @:ALLOCATE(dj_dx(-fd_number:m_end, 0:n, 0:p, 1:Dif_size))
        @:ALLOCATE(djh_dx(-fd_number:m_end, 0:n, 0:p, 1:Dif_size))
        @:ALLOCATE(dY_dx(-2*fd_number:m_end_Y, 0:n, 0:p, 1:Dif_size))
        @:ALLOCATE(dvel_dx(0:m, 0:n, 0:p))
        if (n > 0) then
            @:ALLOCATE(dj_dy(0:m, -fd_number:n_end, 0:p, 1:Dif_size))
            @:ALLOCATE(djh_dy(0:m, -fd_number:n_end, 0:p, 1:Dif_size))
            @:ALLOCATE(dY_dy(0:m, -fd_number:n_end, 0:p, 1:Dif_size))
            @:ALLOCATE(dvel_dy(0:m, 0:n, 0:p))
            if (p > 0) then
                @:ALLOCATE(dj_dz(0:m, 0:n, -fd_number:p_end, 1:Dif_size))
                @:ALLOCATE(djh_dz(0:m, 0:n, -fd_number:p_end, 1:Dif_size))
                @:ALLOCATE(dY_dz(0:m, 0:n, -fd_number:p_end, 1:Dif_size))
                @:ALLOCATE(dvel_dz(0:m, 0:n, 0:p))
            end if
        end if

        @:ALLOCATE(alpha_K_dif(-2*fd_number:m_end_Y, -2*fd_number:n_end_Y, -2*fd_number:p_end_Y, 1:Dif_size))
        @:ALLOCATE(alpharho_K_dif(-2*fd_number:m_end_Y, -2*fd_number:n_end_Y, -2*fd_number:p_end_Y, 1:Dif_size))
        @:ALLOCATE(rho_dif(-2*fd_number:m_end_Y, -2*fd_number:n_end_Y, -2*fd_number:p_end_Y))
        @:ALLOCATE(alpha_dif(-2*fd_number:m_end_Y, -2*fd_number:n_end_Y, -2*fd_number:p_end_Y))
        @:ALLOCATE(T_dif(-2*fd_number:m_end_Y, -2*fd_number:n_end_Y, -2*fd_number:p_end_Y))
        @:ALLOCATE(W_dif(-2*fd_number:m_end_Y, -2*fd_number:n_end_Y, -2*fd_number:p_end_Y))
        @:ALLOCATE(denom(0:m, 0:n, 0:p))
        @:ALLOCATE(rhogcg2(0:m, 0:n, 0:p))
        @:ALLOCATE(rho1c12(0:m, 0:n, 0:p))
        @:ALLOCATE(kdivu(0:m, 0:n, 0:p))

        @:ALLOCATE(fd_coeff_x_d(-fd_number:fd_number,-fd_number:m_end))
        if (n > 0) then
            @:ALLOCATE(fd_coeff_y_d(-fd_number:fd_number, -fd_number:n_end))
        end if
        if (p > 0) then
            @:ALLOCATE(fd_coeff_z_d(-fd_number:fd_number, -fd_number:p_end))
        end if


        ! Computing centered finite difference coefficients
        call s_compute_finite_difference_coefficients(m, x_cc, fd_coeff_x_d, buff_size, &
                                                      fd_number, fd_order, offset_s(1))
        !$acc update device(fd_coeff_x_d)
        if (n > 0) then
            call s_compute_finite_difference_coefficients(n, y_cc, fd_coeff_y_d, buff_size, &
                                                          fd_number, fd_order, offset_s(2))
            !$acc update device(fd_coeff_y_d)
        end if
        if (p > 0) then
            call s_compute_finite_difference_coefficients(p, z_cc, fd_coeff_z_d, buff_size, &
                                                          fd_number, fd_order, offset_s(3))
            !$acc update device(fd_coeff_z_d)
        end if

    end subroutine s_initialize_diffusion_module

    subroutine s_compute_sum_alpha_g(q_cons_vf, bounds)

        ! From the volume fractions of each mixture gas component, compute the
        ! total gas volume fraction field.

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf
        type(int_bounds_info), dimension(1:3), intent(in) :: bounds

        integer :: x, y, z, i
        real(wp) :: sum_alpha_g
 

        sum_alpha_g = 0.0_wp
        do z = bounds(3)%beg, bounds(3)%end
            do y = bounds(2)%beg, bounds(2)%end
                do x = bounds(1)%beg, bounds(1)%end
                    !$acc loop seq
                    sum_alpha_g = 0.0_wp
                    do i = 1, Dif_size
                        sum_alpha_g = sum_alpha_g + q_cons_vf(advxb + Dif_idx(i) - 1)%sf(x, y, z)
                    end do
                    q_cons_vf(advg_idx)%sf(x, y, z) = sum_alpha_g
                end do
            end do
        end do

    end subroutine s_compute_sum_alpha_g

    subroutine s_compute_diffusion_rhs(idir, j_src_n, rhs_vf, q_prim_vf, irx, iry, irz)

        integer, intent(in) :: idir
        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf
        type(scalar_field), dimension(sys_size), intent(inout) :: j_src_n, rhs_vf
        type(int_bounds_info), intent(in) :: irx, iry, irz

        integer :: i, k, l, q, r !< Loop variables
        real(wp) :: W1, W2, W3, D12, D13, D23
        real(wp) :: R_univ
        real(wp) :: grid_spacing
        real(wp) :: rho_L, rho_R, rho_f, rhog_f
        real(wp) :: alpha_m_L, alpha_m_R, alpha_m_f
        real(wp) :: g_f
        real(wp) :: P_L, P_R, P_f
        real(wp) :: T_f, W_f
        real(wp) :: sum_jflux
        integer, dimension(3) :: offsets

        real(wp) :: alpha_L(Dif_size), alpha_R(Dif_size), alpha_f(Dif_size)
        real(wp) :: alpharho_L(Dif_size), alpharho_R(Dif_size), alpharho_f(Dif_size)
        real(wp) :: Y_f(Dif_size), Y_L(Dif_size), Y_R(Dif_size)
        real(wp) :: dY_ds_f(Dif_size)
        real(wp) :: h_f(Dif_size)
        real(wp) :: j_flux(Dif_size)



        ! real(wp), allocatable :: alpha_L(:), alpha_R(:), alpha_f(:)
        ! real(wp), allocatable :: alpharho_L(:), alpharho_R(:), alpharho_f(:)
        ! real(wp), allocatable :: Y_f(:), Y_L(:), Y_R(:)
        ! real(wp), allocatable :: dY_ds_f(:)
        ! real(wp), allocatable :: h_f(:)
        ! real(wp), allocatable :: j_flux(:)





        ! allocate(alpha_L(Dif_size), alpha_R(Dif_size), alpha_f(Dif_size))
        ! allocate(alpharho_L(Dif_size), alpharho_R(Dif_size), alpharho_f(Dif_size))
        ! allocate(Y_f(Dif_size), Y_L(Dif_size), Y_R(Dif_size))
        ! allocate(dY_ds_f(Dif_size))
        ! allocate(h_f(Dif_size))
        ! allocate(j_flux(Dif_size))

        R_univ = 8314.462618_wp

        isd1 = irx; isd2 = iry; isd3 = irz

        ! Set offsets based on direction using array indexing
        offsets = 0
        offsets(idir) = 1
        
        if (Dif_fv) then
            ! Finite Volume with j_src_n Approach
            ! #########################################################################
            ! #########################################################################
            do q = isd3%beg, isd3%end
                do l = isd2%beg, isd2%end
                    do k = isd1%beg, isd1%end

                        do i = 1, Dif_size
                            j_src_n(Dif_idx(i))%sf(k, l, q) = 0._wp
                        end do
                        j_src_n(E_idx)%sf(k, l, q) = 0._wp

                        ! Calculate grid spacing using direction-based indexing
                        select case (idir)
                        case (1)
                            grid_spacing = x_cc(k + 1) - x_cc(k)
                        case (2)
                            grid_spacing = y_cc(l + 1) - y_cc(l)
                        case (3)
                            grid_spacing = z_cc(q + 1) - z_cc(q)
                        end select

                        do i = 1, Dif_size
                            alpha_L(i) = q_prim_vf(advxb + Dif_idx(i) - 1)%sf(k, l, q)
                            alpha_R(i) = q_prim_vf(advxb + Dif_idx(i) - 1)%sf(k + offsets(1), l + offsets(2), q + offsets(3))
                            alpharho_L(i) = q_prim_vf(Dif_idx(i))%sf(k, l, q)
                            alpharho_R(i) = q_prim_vf(Dif_idx(i))%sf(k + offsets(1), l + offsets(2), q + offsets(3))
                        end do

                        alpha_m_L = q_prim_vf(advg_idx)%sf(k, l, q)
                        alpha_m_R = q_prim_vf(advg_idx)%sf(k + offsets(1), l + offsets(2), q + offsets(3))
                        alpha_m_f = 0.5_wp * (alpha_m_L + alpha_m_R)

                        do i = 1, Dif_size
                            alpha_f(i) = 0.5_wp * (alpha_L(i) + alpha_R(i))
                            alpharho_f(i) = 0.5_wp * (alpharho_L(i) + alpharho_R(i))
                        end do

                        rho_L = 0._wp
                        rho_R = 0._wp
                        rho_f = 0._wp

                        do i = 1, Dif_size
                            rho_L = rho_L + alpharho_L(i)
                            rho_R = rho_R + alpharho_R(i)
                            rho_f = rho_f + alpharho_f(i)
                        end do
                    
                        if (alpha_m_L < small_num_dif .or. alpha_m_R < small_num_dif) cycle

                        g_f = 2._wp*alpha_m_L*alpha_m_R / (alpha_m_L + alpha_m_R)
                        g_f = min(alpha_m_R, alpha_m_L)

                        ! g_f = 1.0_wp

                        ! Total gas density at face
                        rhog_f = rho_f / alpha_m_f

                        P_L = q_prim_vf(E_idx)%sf(k, l, q)
                        P_R = q_prim_vf(E_idx)%sf(k + offsets(1), l + offsets(2), q + offsets(3))
                        P_f = 0.5_wp * (P_L + P_R)

                        
                        do i = 1, Dif_size
                            Y_L(i) = alpharho_L(i) / rho_L
                        end do
                        

                        do i = 1, Dif_size
                            Y_R(i) = alpharho_R(i) / rho_R
                        end do
             
                        do i = 1, Dif_size
                            Y_f(i) = alpharho_f(i) / rho_f
                        end do
                        

                        do i = 1, Dif_size
                            dY_ds_f(i) = (Y_R(i) - Y_L(i)) / grid_spacing
                        end do

                        W_f = 0._wp
                        do i = 1, Dif_size  
                            W_f = W_f + Y_f(i)/Ws(i)              
                        end do

                        W_f = 1._wp / W_f

                        T_f = P_f * W_f / (rhog_f * R_univ)

                        do i = 1, Dif_size
                            h_f(i) = h0s(i) + cps(i)*(T_f - T0s(i))
                        end do

                        ! Compute diffusion fluxes
                        if (Dif_size == 2) then
                            j_flux(1) = -rhog_f*Ds(1,2)*dY_ds_f(1)
                            j_flux(2) = -j_flux(1)
                        else if (Dif_size == 3) then
                            j_flux(1) = -rhog_f / (Y_f(1)*Ds(2,3) + Y_f(2)*Ds(3,1) + Y_f(3)*Ds(1,2)) * &
                                            ( Ds(1,2)*Ds(1,3)*dY_ds_f(1)*(1._wp - Y_f(1)) - Y_f(1)*Ds(2,3)*(Ds(1,2)*dY_ds_f(2) + Ds(1,3)*dY_ds_f(3)) )

                            j_flux(2) = -rhog_f / (Y_f(1)*Ds(2,3) + Y_f(2)*Ds(3,1) + Y_f(3)*Ds(1,2)) * &
                                            ( Ds(2,1)*Ds(2,3)*dY_ds_f(2)*(1._wp - Y_f(2)) - Y_f(2)*Ds(3,1)*(Ds(2,1)*dY_ds_f(1) + Ds(2,3)*dY_ds_f(3)) )

                            j_flux(3) = -sum(j_flux(1:2))
                        else
                            call s_calculate_multicomponent_diffusion_flux(Dif_size, rhog_f, Y_f, dY_ds_f, j_flux)
                        end if

                        ! Enforce mass conservation of diffusion fluxes
                        ! sum_jflux = 0.0_wp
                        ! if (alpha_m_f > small_num_dif) then
                        !     do i = 1, Dif_size
                        !         sum_jflux = sum_jflux + j_flux(i)
                        !     end do

                        !     do i = 1, Dif_size
                        !         j_flux(i) = j_flux(i) - Y_f(i)*sum_jflux
                        !     end do
                        ! end if

                        do i = 1, Dif_size
                            j_src_n(Dif_idx(i))%sf(k, l, q) = j_src_n(Dif_idx(i))%sf(k, l, q) + g_f*j_flux(i)
                            j_src_n(E_idx)%sf(k, l, q) = j_src_n(E_idx)%sf(k, l, q) + g_f*h_f(i)*j_flux(i)
                        end do
                    end do
                end do
            end do
        ! #########################################################################
        ! #########################################################################
        else
            ! Pure Finite Difference Approach 
            ! #########################################################################
            ! #########################################################################
            if (idir == 1) then
                !$acc parallel loop collapse(4) gang vector default(present)
                do q = 0, p
                    do l = 0, n
                        do k = -2*fd_number, m + 2*fd_number
                            do i = 1, Dif_size
                                alpha_K_dif(k, l, q, i) = q_prim_vf(advxb + Dif_idx(i) - 1)%sf(k, l, q)
                                alpharho_K_dif(k, l, q, i) = q_prim_vf(Dif_idx(i))%sf(k, l, q)
                            end do
                            rho_dif(k, l, q) = 0._wp
                            alpha_dif(k, l, q) = 0._wp
                        end do
                    end do
                end do
                !$acc end parallel loop

                !$acc parallel loop collapse(4) gang vector default(present)
                do q = 0, p
                    do l = 0, n
                        do k = -2*fd_number, m + 2*fd_number
                            do i = 1, Dif_size
                                rho_dif(k, l, q) = rho_dif(k, l, q) + alpharho_K_dif(k, l, q, i)
                                alpha_dif(k, l, q) = alpha_dif(k, l, q) + alpha_K_dif(k, l, q, i)
                            end do
                        end do
                    end do
                end do
                !$acc end parallel loop

                !$acc parallel loop collapse(3) gang vector default(present)
                do q = 0, p
                    do l = 0, n
                        do k = -2*fd_number, m + 2*fd_number
                            W_dif(k, l, q) = 0._wp
                            do i = 1, Dif_size
                                W_dif(k, l, q) = W_dif(k, l, q) + alpha_K_dif(k, l, q, i)*Ws(i)
                            end do
                        end do
                    end do
                end do
                !$acc end parallel loop

                !$acc parallel loop collapse(4) gang vector default(present)
                do q = 0, p
                    do l = 0, n
                        do k = -2*fd_number, m + 2*fd_number
                            ! gas cell
                            if (alpha_dif(k, l, q) > small_num_dif) then
                                T_dif(k, l, q) = q_prim_vf(E_idx)%sf(k, l, q) * W_dif(k, l, q) /( rho_dif(k, l, q)*R_univ )
                            else
                                T_dif(k, l, q) = 0._wp
                            end if
                        end do
                    end do
                end do
                !$acc end parallel loop

                !$acc parallel loop collapse(4) gang vector default(present)
                do q = 0, p
                    do l = 0, n
                        do k = -2*fd_number, m + 2*fd_number
                            if (alpha_dif(k, l, q) > small_num_dif) then
                                do i = 1, Dif_size
                                    
                                    Y_dif(k, l, q, i) = alpharho_K_dif(k, l, q, i) / rho_dif(k, l, q)
                                    !h_dif(k, l, q, i) = (q_prim_vf(E_idx)%sf(k, l, q) * (gammas(Dif_idx(i)) + 1._wp)) * alpha_K_dif(k, l, q, i) / alpharho_K_dif(k, l, q, i)
                                    h_dif(k, l, q, i) = h0s(i) + cps(i)*(T_dif(k, l, q) - T0s(i))
                                
                                end do
                            end if
                        end do
                    end do
                end do
                !$acc end parallel loop

                !$acc parallel loop collapse(4) gang vector default(present)
                do q = 0, p
                    do l = 0, n
                        do k = -fd_number, m + fd_number
                            do i = 1, Dif_size
                                dj_dx(k, l, q, i) = 0._wp
                                djh_dx(k, l, q, i) = 0._wp
                                dY_dx(k, l, q, i) = 0._wp
                            end do
                        end do
                    end do
                end do
                !$acc end parallel loop

                !set ghost cell for buffer region equal to k point (to enforce del dot j = 0)
                !$acc parallel loop collapse(5) gang vector default(present)
                do q = 0, p
                    do l = 0, n
                        do k = -fd_number, m + fd_number
                            if (alpha_dif(k, l, q) > small_num_dif) then
                                do i = 1, Dif_size
                                    do r = -fd_number, fd_number
                                        if (alpha_dif(k + r, l, q) > small_num_dif) then
                                            dY_dx(k, l, q, i) = dY_dx(k, l, q, i) &
                                                + Y_dif(k + r, l, q, i)*fd_coeff_x_d(r, k)
                                        else
                                            dY_dx(k, l, q, i) = dY_dx(k, l, q, i) &
                                                + Y_dif(k, l, q, i)*fd_coeff_x_d(r, k)
                                        end if
                                    end do
                                end do
                            end if
                        end do
                    end do
                end do
                !$acc end parallel loop

                if (Dif_size == 2) then
                    !$acc parallel loop collapse(5) gang vector default(present)
                    do q = 0, p
                        do l = 0, n
                            do k = 0, m
                                if (alpha_dif(k, l, q) > small_num_dif) then
                                    do i = 1, Dif_size
                                        do r = -fd_number, fd_number
                                            if (alpha_dif(k + r, l, q) > small_num_dif) then
                                                dj_dx(k, l, q, i) = dj_dx(k, l, q, i) &
                                                    + dY_dx(k + r, l, q, i)*rho_dif(k + r, l, q)*Ds(1,2)*fd_coeff_x_d(r, k)
                                                djh_dx(k, l, q, i) = djh_dx(k, l, q, i) &
                                                    + h_dif(k + r, l, q, i)*dY_dx(k + r, l, q, i)*rho_dif(k + r, l, q)*Ds(1,2)*fd_coeff_x_d(r, k)

                                            else
                                                dj_dx(k, l, q, i) = dj_dx(k, l, q, i) &
                                                    + dY_dx(k, l, q, i)*rho_dif(k, l, q)*Ds(1,2)*fd_coeff_x_d(r, k)
                                                djh_dx(k, l, q, i) = djh_dx(k, l, q, i) &
                                                    + h_dif(k, l, q, i)*dY_dx(k, l, q, i)*rho_dif(k, l, q)*Ds(1,2)*fd_coeff_x_d(r, k)
                                            end if
                                        end do
                                    end do
                                end if
                            end do
                        end do
                    end do
                    !$acc end parallel loop
                else if (Dif_size == 3) then
                    !$acc parallel loop collapse(5) gang vector default(present)
                    do q = 0, p
                        do l = 0, n
                            do k = 0, m
                                if (alpha_dif(k, l, q) > small_num_dif) then
                                    do r = -fd_number, fd_number
                                        do i = 1, Dif_size
                                            if (alpha_dif(k + r, l, q) > small_num_dif) then
                                                select case (i)
                                                    case (1)
                                                        dj_dx(k, l, q, i) = dj_dx(k, l, q, i) &
                                                            + rho_dif(k + r, l, q)*(Ds(1,2)*Ds(1,3)*(1._wp - Y_dif(k + r, l, q, 1))*dY_dx(k + r, l, q, 1) &
                                                            - Y_dif(k + r, l, q, 1)*Ds(2,3)*(Ds(1,2)*dY_dx(k + r, l, q, 2) + Ds(1,3)*dY_dx(k + r, l, q, 3))) &
                                                            / ( Y_dif(k + r, l, q, 1)*Ds(2,3) + Y_dif(k + r, l, q, 2)*Ds(1,3) + Y_dif(k + r, l, q, 3)*Ds(1,2) )*fd_coeff_x_d(r, k)

                                                        djh_dx(k, l, q, i) = djh_dx(k, l, q, i) &
                                                            + h_dif(k + r, l, q, i)*rho_dif(k + r, l, q)*(Ds(1,2)*Ds(1,3)*(1._wp - Y_dif(k + r, l, q, 1))*dY_dx(k + r, l, q, 1) &
                                                            - Y_dif(k + r, l, q, 1)*Ds(2,3)*(Ds(1,2)*dY_dx(k + r, l, q, 2) + Ds(1,3)*dY_dx(k + r, l, q, 3))) &
                                                            / ( Y_dif(k + r, l, q, 1)*Ds(2,3) + Y_dif(k + r, l, q, 2)*Ds(1,3) + Y_dif(k + r, l, q, 3)*Ds(1,2) )*fd_coeff_x_d(r, k)

                                                    case (2)
                                                        dj_dx(k, l, q, i) = dj_dx(k, l, q, i) &
                                                            + rho_dif(k + r, l, q)*(Ds(1,2)*Ds(2,3)*(1._wp - Y_dif(k + r, l, q, 2))*dY_dx(k + r, l, q, 2) &
                                                            - Y_dif(k + r, l, q, 2)*Ds(1,3)*(Ds(1,2)*dY_dx(k + r, l, q, 1) + Ds(2,3)*dY_dx(k + r, l, q, 3))) &
                                                            / ( Y_dif(k + r, l, q, 1)*Ds(2,3) + Y_dif(k + r, l, q, 2)*Ds(1,3) + Y_dif(k + r, l, q, 3)*Ds(1,2) )*fd_coeff_x_d(r, k)

                                                        djh_dx(k, l, q, i) = djh_dx(k, l, q, i) &
                                                            + h_dif(k + r, l, q, i)*rho_dif(k + r, l, q)*(Ds(1,2)*Ds(1,3)*(1._wp - Y_dif(k + r, l, q, 1))*dY_dx(k + r, l, q, 1) &
                                                            - Y_dif(k + r, l, q, 1)*Ds(2,3)*(Ds(1,2)*dY_dx(k + r, l, q, 2) + Ds(1,3)*dY_dx(k + r, l, q, 3))) &
                                                            / ( Y_dif(k + r, l, q, 1)*Ds(2,3) + Y_dif(k + r, l, q, 2)*Ds(1,3) + Y_dif(k + r, l, q, 3)*Ds(1,2) )*fd_coeff_x_d(r, k)

                                                    case (3)
                                                        dj_dx(k, l, q, i) = dj_dx(k, l, q, i) &
                                                            + rho_dif(k + r, l, q)*(Ds(2,3)*Ds(1,3)*(1._wp - Y_dif(k + r, l, q, 3))*dY_dx(k + r, l, q, 3) &
                                                            - Y_dif(k + r, l, q, 3)*Ds(1,2)*(Ds(1,3)*dY_dx(k + r, l, q, 1) + Ds(2,3)*dY_dx(k + r, l, q, 2))) &
                                                            / ( Y_dif(k + r, l, q, 1)*Ds(2,3) + Y_dif(k + r, l, q, 2)*Ds(1,3) + Y_dif(k + r, l, q, 3)*Ds(1,2) )*fd_coeff_x_d(r, k)

                                                        djh_dx(k, l, q, i) = djh_dx(k, l, q, i) &
                                                            + h_dif(k + r, l, q, i)*rho_dif(k + r, l, q)*(Ds(2,3)*Ds(1,3)*(1._wp - Y_dif(k + r, l, q, 3))*dY_dx(k + r, l, q, 3) &
                                                            - Y_dif(k + r, l, q, 3)*Ds(1,2)*(Ds(1,3)*dY_dx(k + r, l, q, 1) + Ds(2,3)*dY_dx(k + r, l, q, 2))) &
                                                            / ( Y_dif(k + r, l, q, 1)*Ds(2,3) + Y_dif(k + r, l, q, 2)*Ds(1,3) + Y_dif(k + r, l, q, 3)*Ds(1,2) )*fd_coeff_x_d(r, k)
                                                end select
                                        
                                            else
                                                select case (i)
                                                    case (1)
                                                        dj_dx(k, l, q, i) = dj_dx(k, l, q, i) &
                                                            + rho_dif(k, l, q)*(Ds(1,2)*Ds(1,3)*(1._wp - Y_dif(k, l, q, 1))*dY_dx(k, l, q, 1) &
                                                            - Y_dif(k, l, q, 1)*Ds(2,3)*(Ds(1,2)*dY_dx(k, l, q, 2) + Ds(1,3)*dY_dx(k, l, q, 3))) &
                                                            / ( Y_dif(k, l, q, 1)*Ds(2,3) + Y_dif(k, l, q, 2)*Ds(1,3) + Y_dif(k, l, q, 3)*Ds(1,2) )*fd_coeff_x_d(r, k)

                                                        djh_dx(k, l, q, i) = djh_dx(k, l, q, i) &
                                                            + h_dif(k, l, q, i)*rho_dif(k, l, q)*(Ds(2,3)*Ds(1,3)*(1._wp - Y_dif(k, l, q, 3))*dY_dx(k, l, q, 3) &
                                                            - Y_dif(k, l, q, 3)*Ds(1,2)*(Ds(1,3)*dY_dx(k, l, q, 1) + Ds(2,3)*dY_dx(k, l, q, 2))) &
                                                            / ( Y_dif(k, l, q, 1)*Ds(2,3) + Y_dif(k, l, q, 2)*Ds(1,3) + Y_dif(k, l, q, 3)*Ds(1,2) )*fd_coeff_x_d(r, k)

                                                    case (2)
                                                        dj_dx(k, l, q, i) = dj_dx(k, l, q, i) &
                                                            + rho_dif(k, l, q)*(Ds(1,2)*Ds(2,3)*(1._wp - Y_dif(k, l, q, 2))*dY_dx(k, l, q, 2) &
                                                            - Y_dif(k, l, q, 2)*Ds(1,3)*(Ds(1,2)*dY_dx(k, l, q, 1) + Ds(2,3)*dY_dx(k, l, q, 3))) &
                                                            / ( Y_dif(k, l, q, 1)*Ds(2,3) + Y_dif(k, l, q, 2)*Ds(1,3) + Y_dif(k, l, q, 3)*Ds(1,2) )*fd_coeff_x_d(r, k)

                                                        djh_dx(k, l, q, i) = djh_dx(k, l, q, i) &
                                                            + h_dif(k, l, q, i)*rho_dif(k, l, q)*(Ds(2,3)*Ds(1,3)*(1._wp - Y_dif(k, l, q, 3))*dY_dx(k, l, q, 3) &
                                                            - Y_dif(k, l, q, 3)*Ds(1,2)*(Ds(1,3)*dY_dx(k, l, q, 1) + Ds(2,3)*dY_dx(k, l, q, 2))) &
                                                            / ( Y_dif(k, l, q, 1)*Ds(2,3) + Y_dif(k, l, q, 2)*Ds(1,3) + Y_dif(k, l, q, 3)*Ds(1,2) )*fd_coeff_x_d(r, k)

                                                    case (3)
                                                        dj_dx(k, l, q, i) = dj_dx(k, l, q, i) &
                                                            + rho_dif(k, l, q)*(Ds(2,3)*Ds(1,3)*(1._wp - Y_dif(k, l, q, 3))*dY_dx(k, l, q, 3) &
                                                            - Y_dif(k, l, q, 3)*Ds(1,2)*(Ds(1,3)*dY_dx(k, l, q, 1) + Ds(2,3)*dY_dx(k, l, q, 2))) &
                                                            / ( Y_dif(k, l, q, 1)*Ds(2,3) + Y_dif(k, l, q, 2)*Ds(1,3) + Y_dif(k, l, q, 3)*Ds(1,2) )*fd_coeff_x_d(r, k)

                                                        djh_dx(k, l, q, i) = djh_dx(k, l, q, i) &
                                                            + h_dif(k, l, q, i)*rho_dif(k, l, q)*(Ds(2,3)*Ds(1,3)*(1._wp - Y_dif(k, l, q, 3))*dY_dx(k, l, q, 3) &
                                                            - Y_dif(k, l, q, 3)*Ds(1,2)*(Ds(1,3)*dY_dx(k, l, q, 1) + Ds(2,3)*dY_dx(k, l, q, 2))) &
                                                            / ( Y_dif(k, l, q, 1)*Ds(2,3) + Y_dif(k, l, q, 2)*Ds(1,3) + Y_dif(k, l, q, 3)*Ds(1,2) )*fd_coeff_x_d(r, k)
                                                end select

                                            end if
                                        end do
                                        
                                    end do
                                end if 
                            end do
                        end do
                    end do
                    !$acc end parallel loop
                end if

                !Valid for any number of species
                ! species continuity
                !$acc parallel loop collapse(4) gang vector default(present)
                do q = 0, p
                    do l = 0, n
                        do k = 0, m
                            if (alpha_dif(k, l, q) > small_num_dif) then
                                do i = 1, Dif_size
                                    rhs_vf(Dif_idx(i))%sf(k, l, q) = rhs_vf(Dif_idx(i))%sf(k, l, q) &
                                        + dj_dx(k, l, q, i)
                                end do
                            end if
                        end do
                    end do
                end do
                !$acc end parallel loop

                !Valid for any number of species
                !energy
                !$acc parallel loop collapse(4) gang vector default(present)
                do q = 0, p
                    do l = 0, n
                        do k = 0, m
                            if (alpha_dif(k, l, q) > small_num_dif) then
                                do i = 1, Dif_size
                                    rhs_vf(E_idx)%sf(k, l, q) = rhs_vf(E_idx)%sf(k, l, q) &
                                        + djh_dx(k, l, q, i)
                                end do
                            end if
                        end do
                    end do
                end do
                !$acc end parallel loop
            end if
        end if
        ! #########################################################################
        ! #########################################################################
        
    end subroutine s_compute_diffusion_rhs

    subroutine s_correct_volume_fractions(q_cons_vf, q_prim_vf)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf
        type(scalar_field), dimension(sys_size), intent(inout) :: q_prim_vf

        integer :: x, y, z, i
        real(wp) :: rho, W
        real(wp) :: alpharho(Dif_size), Y_s(Dif_size)

        do z = 0, p
            do y = 0, n
                do x = 0, m 
                    if (q_cons_vf(advg_idx)%sf(x, y, z) < small_num_dif) then
                        do i = 1, Dif_size
                            q_cons_vf(advxb + Dif_idx(i) - 1)%sf(x, y, z) = 0._wp
                        end do
                    else
                        rho = 0._wp
                        W = 0._wp

                        do i = 1, Dif_size
                            alpharho(i) = q_cons_vf(Dif_idx(i))%sf(x, y, z)
                        end do
 
                        do i = 1, Dif_size
                            rho = rho + alpharho(i)
                        end do
                        
                        do i = 1, Dif_size
                            Y_s(i) = alpharho(i) / rho
                        end do

                        do i = 1, Dif_size
                            W = W + Y_s(i)/Ws(i)
                        end do

                        W = 1._wp / W

                        do i = 1, Dif_size
                            q_cons_vf(advxb + Dif_idx(i) - 1)%sf(x, y, z) = q_cons_vf(advg_idx)%sf(x, y, z) * Y_s(i) * W / Ws(i)
                        end do

                        do i = 1, Dif_size
                            q_prim_vf(advxb + Dif_idx(i) - 1)%sf(x, y, z) = q_cons_vf(advxb + Dif_idx(i) - 1)%sf(x, y, z)
                        end do

                        q_prim_vf(advg_idx)%sf(x, y, z) = q_cons_vf(advg_idx)%sf(x, y, z)
                    end if
                end do
            end do
        end do

    end subroutine s_correct_volume_fractions

    subroutine s_correct_riemann_volume_fractions(q_rs_vf, bounds)

        real(wp), dimension(startx:, starty:, startz:, 1:), intent(inout) :: q_rs_vf
        type(int_bounds_info), dimension(1:3), intent(in) :: bounds


        integer :: x, y, z, i
        real(wp) :: rho, W
        real(wp) :: alpharho(Dif_size), Y_s(Dif_size)

        do z = bounds(3)%beg, bounds(3)%end
            do y = bounds(2)%beg, bounds(2)%end
                do x = bounds(1)%beg, bounds(1)%end
                    if (q_rs_vf(x, y, z, advg_idx) < small_num_dif) then
                        do i = 1, Dif_size
                            q_rs_vf(x, y, z, advxb + Dif_idx(i) - 1) = 0._wp
                        end do
                    else
                        rho = 0._wp
                        W = 0._wp

                        do i = 1, Dif_size
                            alpharho(i) = q_rs_vf(x, y, z, Dif_idx(i))
                        end do
 
                        do i = 1, Dif_size
                            rho = rho + alpharho(i)
                        end do
                        
                        do i = 1, Dif_size
                            Y_s(i) = alpharho(i) / rho
                        end do

                        do i = 1, Dif_size
                            W = W + Y_s(i)/Ws(i)
                        end do

                        W = 1._wp / W

                        do i = 1, Dif_size
                            q_rs_vf(x, y, z, advxb + Dif_idx(i) - 1) = q_rs_vf(x, y, z, advg_idx) * Y_s(i) * W / Ws(i)
                        end do

                    end if
                end do
            end do
        end do

    end subroutine s_correct_riemann_volume_fractions

    subroutine s_calculate_multicomponent_diffusion_flux(Dif_size, rho, Y, dY_ds, j_flux)

        integer, intent(in) :: Dif_size
        real(wp), intent(in) :: rho
        real(wp), dimension(Dif_size), intent(in) :: Y
        real(wp), dimension(Dif_size), intent(in) :: dY_ds
        real(wp), dimension(Dif_size), intent(out) :: j_flux

        integer :: i, j, k, Nm1, info
        integer :: ipiv(Dif_size-1)
        real(wp) :: A(Dif_size - 1, Dif_size - 1)
        real(wp) :: b(Dif_size - 1)
        real(wp) :: sum_diag, detA

        Nm1 = Dif_size - 1

        do i = 1, Nm1
            b(i) = rho * dY_ds(i)
        end do

        do i = 1, Nm1

            ! Diagonal A(i,i):
            !  - sum_{k!=i, k=1..N-1} Y_k / D(i,k)
            !  - (Y_i + Y_N) / D(i,N)
            sum_diag = 0.0_wp  

            do k = 1, Nm1
                if (k == i) cycle
                sum_diag = sum_diag + Y(k) / Ds(i, k)
            end do

            sum_diag = sum_diag + (Y(i) + Y(Dif_size)) / Ds(i, Dif_size)

            A(i, i) = -sum_diag

            ! Off-diagonals A(i,j) = Y_i * (1/D(i,j) - 1/D(i,N))
            do j = 1, Nm1
                if (j == i) cycle
                A(i, j) = Y(i) * ( 1.0_wp/Ds(i, j) - 1.0_wp/Ds(i, Dif_size) )
            end do

        end do

        ! if (Dif_size == 3) then

        ! else if (Dif_size == 4) then
        
        ! else

        ! Solve the linear system A * j = b
        if (Nm1 == 2) then
            ! For 2x2 system, use explicit formula
            detA = A(1,1)*A(2,2) - A(1,2)*A(2,1)
            j_flux(1) = ( A(2,2)*b(1) - A(1,2)*b(2) ) / detA
            j_flux(2) = ( -A(2,1)*b(1) + A(1,1)*b(2) ) / detA
        
        elseif (Nm1 == 3) then
            ! For 3x3 system, use explicit formula
            detA = A(1,1)*(A(2,2)*A(3,3) - A(2,3)*A(3,2)) - A(1,2)*(A(2,1)*A(3,3) - A(2,3)*A(3,1)) + A(1,3)*(A(2,1)*A(3,2) - A(2,2)*A(3,1))
            j_flux(1) = ( (A(2,2)*A(3,3) - A(2,3)*A(3,2))*b(1) - (A(1,2)*A(3,3) - A(1,3)*A(3,2))*b(2) + (A(1,2)*A(2,3) - A(1,3)*A(2,2))*b(3) ) / detA
            j_flux(2) = ( -(A(2,1)*A(3,3) - A(2,3)*A(3,1))*b(1) + (A(1,1)*A(3,3) - A(1,3)*A(3,1))*b(2) - (A(1,1)*A(2,3) - A(1,3)*A(2,1))*b(3) ) / detA
            j_flux(3) = ( (A(2,1)*A(3,2) - A(2,2)*A(3,1))*b(1) - (A(1,1)*A(3,2) - A(1,2)*A(3,1))*b(2) + (A(1,1)*A(2,2) - A(1,2)*A(2,1))*b(3) ) / detA
        else
            ! For larger systems, use LAPACK DGESV
            call dgesv(Nm1, 1, A, Nm1, ipiv, b, Nm1, info)
            do i = 1, Nm1
                j_flux(i) = b(i)
            end do
        end if 

        

        j_flux(Dif_size) = -sum(j_flux(1:Nm1))

    end subroutine s_calculate_multicomponent_diffusion_flux

    subroutine s_finalize_diffusion_module

        @:DEALLOCATE(Ds)
        @:DEALLOCATE(Ws)
        @:DEALLOCATE(cps)
        @:DEALLOCATE(T0s)
        @:DEALLOCATE(h0s)
        @:DEALLOCATE(alpha_K_dif, alpharho_K_dif, Y_dif, h_dif, T_dif)
        @:DEALLOCATE(rho_dif, alpha_dif, denom, rhogcg2, rho1c12, kdivu, W_dif)
        @:DEALLOCATE(fd_coeff_x_d)
        @:DEALLOCATE(dj_dx)
        @:DEALLOCATE(dY_dx)
        @:DEALLOCATE(djh_dx)
        @:DEALLOCATE(dvel_dx)
        if (n > 0) then
            @:DEALLOCATE(fd_coeff_y_d)
            @:DEALLOCATE(dj_dy)
            @:DEALLOCATE(djh_dy)
            @:DEALLOCATE(dY_dy)
            @:DEALLOCATE(dvel_dy)
            if (p > 0) then
                @:DEALLOCATE(fd_coeff_z_d)
                @:DEALLOCATE(dj_dz)
                @:DEALLOCATE(djh_dz)
                @:DEALLOCATE(dY_dz)
                @:DEALLOCATE(dvel_dz)
            end if
        end if

    end subroutine s_finalize_diffusion_module


end module m_diffusion