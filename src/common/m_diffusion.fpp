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

    use m_weno                 !< WENO module

    use m_helper              !< Helper functions

    ! ==========================================================================

    implicit none

    private; public :: s_compute_sum_alpha_g, &
s_initialize_diffusion_module, &
s_compute_diffusion_rhs, &
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

    subroutine s_compute_sum_alpha_g(q_cons_vf, bounds)

        ! From the volume fractions of each mixture gas component, compute the
        ! total gas volume fraction field.

        type(scalar_field), dimension(sys_size), intent(in) :: q_cons_vf
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

    subroutine s_compute_diffusion_rhs(idir, j_src_n, rhs_vf, q_prim_vf, irx, iry, irz)

        integer, intent(in) :: idir
        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf
        type(scalar_field), dimension(sys_size), intent(inout) :: j_src_n, rhs_vf
        type(int_bounds_info), intent(in) :: irx, iry, irz

        integer :: i, k, l, q, r !< Loop variables
        real(wp) :: W1, W2, W3, D12, D13, D23
        real(wp) :: R_univ, small_number
        real(wp) :: grid_spacing
        real(wp) :: rho_L, rho_R, rho_f
        real(wp) :: alpha_m_L, alpha_m_R, alpha_m_f
        real(wp) :: P_L, P_R
        real(wp) :: T_L, T_R, T_f, W_L, W_R, W_f
        real(wp) :: sum_jflux
        real(wp), allocatable :: alpha_L(:), alpha_R(:), alpha_f(:)
        real(wp), allocatable :: alpharho_L(:), alpharho_R(:), alpharho_f(:)
        real(wp), allocatable :: Y_L(:), Y_R(:), Y_f(:)
        real(wp), allocatable :: dY_ds_f(:)
        real(wp), allocatable :: h_f(:)
        real(wp), allocatable :: K_L(:), K_R(:), K_eff(:)
        real(wp), allocatable :: j_flux(:)
        real(wp), allocatable :: alpha_flux(:)
        real(wp), allocatable :: alpha_nonconserv(:)
        real(wp) :: d, s, c, sigma_max, sigma, gamma
        integer, dimension(3) :: offsets

        allocate(alpha_L(Dif_size), alpha_R(Dif_size), alpha_f(Dif_size))
        allocate(alpharho_L(Dif_size), alpharho_R(Dif_size), alpharho_f(Dif_size))
        allocate(Y_L(Dif_size), Y_R(Dif_size), Y_f(Dif_size))
        allocate(dY_ds_f(Dif_size))
        allocate(h_f(Dif_size))
        allocate(K_L(Dif_size), K_R(Dif_size), K_eff(Dif_size))
        allocate(j_flux(Dif_size))
        allocate(alpha_flux(Dif_size))
        allocate(alpha_nonconserv(Dif_size))

        R_univ = 8314.3_wp

        small_number = 1.0e-8_wp

        isd1 = irx; isd2 = iry; isd3 = irz

        ! Set offsets based on direction using array indexing
        offsets = 0
        offsets(idir) = 1

        do q = isd3%beg, isd3%end
            do l = isd2%beg, isd2%end
                do k = isd1%beg, isd1%end

                    do i = 1, Dif_size
                        j_src_n(Dif_idx(i))%sf(k, l, q) = 0._wp
                        j_src_n(advxb + Dif_idx(i) - 1)%sf(k, l, q) = 0._wp
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

                    do i = 1, Dif_size
                        alpha_f(i) = 0.5_wp * (alpha_L(i) + alpha_R(i))
                        alpharho_f(i) = 0.5_wp * (alpharho_L(i) + alpharho_R(i))

                        !alpha_f(i) = 2*alpha_L(i)*alpha_R(i) / (alpha_L(i) + alpha_R(i))
                        !alpharho_f(i) = 2*alpharho_L(i)*alpharho_R(i) / (alpharho_L(i) + alpharho_R(i))
                    end do

                    rho_L = 0._wp
                    rho_R = 0._wp
                    rho_f = 0._wp
                    alpha_m_L = 0._wp
                    alpha_m_R = 0._wp
                    alpha_m_f = 0._wp

                    do i = 1, Dif_size
                        rho_L = rho_L + alpharho_L(i)
                        rho_R = rho_R + alpharho_R(i)
                        !rho_f = rho_f + alpharho_f(i)
                        alpha_m_L = alpha_m_L + alpha_L(i)
                        alpha_m_R = alpha_m_R + alpha_R(i)
                        alpha_m_f = alpha_m_f + alpha_f(i)
                    end do
                    rho_f = 0.5_wp * (rho_L + rho_R)

                    P_L = q_prim_vf(E_idx)%sf(k, l, q)
                    P_R = q_prim_vf(E_idx)%sf(k + offsets(1), l + offsets(2), q + offsets(3))

                    if (alpha_m_L > small_number) then
                        do i = 1, Dif_size
                            Y_L(i) = alpharho_L(i) / rho_L
                        end do
                    else
                        do i = 1, Dif_size
                            Y_L(i) = 0._wp
                        end do
                    end if

                    if (alpha_m_R > small_number) then
                        do i = 1, Dif_size
                            Y_R(i) = alpharho_R(i) / rho_R
                        end do
                    else
                        do i = 1, Dif_size
                            Y_R(i) = 0._wp
                        end do
                    end if

                    if (alpha_m_f > small_number) then
                        do i = 1, Dif_size
                            Y_f(i) = alpharho_f(i) / rho_f
                        end do
                    else
                        do i = 1, Dif_size
                            Y_f(i) = 0._wp
                        end do
                    end if

                    do i = 1, Dif_size
                        dY_ds_f(i) = (Y_R(i) - Y_L(i)) / grid_spacing
                    end do

                    W_L = 0._wp
                    W_R = 0._wp
                    W_f = 0._wp
                    do i = 1, Dif_size
                        W_L = W_L + Y_L(i)/Ws(i)
                        W_R = W_R + Y_R(i)/Ws(i)      
                        W_f = W_f + Y_f(i)/Ws(i)
              
                    end do

                    W_L = 1._wp / W_L
                    W_R = 1._wp / W_R
                    W_f = 1._wp / W_f

                    if (alpha_m_L > small_number) then
                        T_L = P_L * W_L / (rho_L * R_univ)
                    else
                        T_L = 0._wp
                    end if

                    if (alpha_m_R > small_number) then
                        T_R = P_R * W_R / (rho_R * R_univ)
                    else
                        T_R = 0._wp
                    end if

                    T_f = 0.5_wp * (T_L + T_R)
                    !T_f = 298._wp
                    !T_f = 2*T_L*T_R / (T_L + T_R)

                    if (alpha_m_f > small_number) then
                        do i = 1, Dif_size
                            h_f(i) = h0s(i) + cps(i)*(T_f - T0s(i))
                        end do
                    else
                        do i = 1, Dif_size
                            h_f(i) = 0._wp
                        end do
                    end if 
                    K_L(1) = rho_L*Ds(1,2)
                    K_L(2) = rho_L*Ds(2,1)
                    K_R(1) = rho_R*Ds(1,2)
                    K_R(2) = rho_R*Ds(2,1)
                    K_eff(1) = 2.0_wp / (1.0_wp/K_L(1) + 1.0_wp/K_R(1))
                    K_eff(2) = 2.0_wp / (1.0_wp/K_L(2) + 1.0_wp/K_R(2))
                    !j_flux(1) = -K_eff(1)*dY_ds_f(1)
                    !j_flux(2) = -j_flux(1)
                    j_flux(1) = -rho_f*Ds(1,2)*dY_ds_f(1)
                    j_flux(2) = -j_flux(1)

                    ! Enforce mass conservation of diffusion fluxes
                    !sum_jflux = 0.0_wp
                    !if (alpha_m_f > small_number) then
                        !do i = 1, Dif_size
                            !sum_jflux = sum_jflux + j_flux(i)
                        !end do

                        !do i = 1, Dif_size
                            !j_flux(i) = j_flux(i) - Y_f(i)*sum_jflux
                        !end do
                    !end if

                    alpha_flux(1) = j_flux(1)*W_f**2._wp / (rho_f*Ws(1)*Ws(2)*alpha_m_f)
                    alpha_flux(2) = j_flux(2)*W_f**2._wp / (rho_f*Ws(1)*Ws(2)*alpha_m_f)

                    alpha_nonconserv(1) = j_flux(1)*(W_R**2._wp/(rho_R*Ws(1)*Ws(2)*alpha_m_R) - W_L**2._wp/(rho_L*Ws(1)*Ws(2)*alpha_m_L)) / grid_spacing
                    alpha_nonconserv(2) = j_flux(2)*(W_R**2._wp/(rho_R*Ws(1)*Ws(2)*alpha_m_R) - W_L**2._wp/(rho_L*Ws(1)*Ws(2)*alpha_m_L)) / grid_spacing

                    d = min(abs(x_cc(k)), abs(x_cc(m) - x_cc(k)))
                    s = max(0.0_wp, 1._wp - d / x_cc(m))
                    gamma = 0.0_wp
                    do i = 1, Dif_size
                        gamma = gamma + alpha_L(i)*(gammas(Dif_idx(i)) + 1._wp) / gammas(Dif_idx(i))
                    end do
                    c = sqrt(gamma*R_univ*T_L/W_L)
                    sigma_max = 3._wp*c/x_cc(m)
                    sigma = sigma_max*s**2._wp*(3._wp - 2._wp*s)

                    do i = 1, Dif_size
                        j_src_n(Dif_idx(i))%sf(k, l, q) = j_src_n(Dif_idx(i))%sf(k, l, q) + j_flux(i)
                        j_src_n(E_idx)%sf(k, l, q) = j_src_n(E_idx)%sf(k, l, q) + h_f(i)*j_flux(i)
                        j_src_n(advxb + Dif_idx(i) - 1)%sf(k, l, q) = j_src_n(advxb + Dif_idx(i) - 1)%sf(k, l, q) + alpha_flux(i)
                        if ((k == 0 .or. offsets(1)*k < m) .and. (l == 0 .or. offsets(2)*l < n) .and. (q == 0 .or. offsets(3)*q < p)) then
                            rhs_vf(advxb + Dif_idx(i) - 1)%sf(k + offsets(1), l + offsets(2), q + offsets(3)) = &
                                rhs_vf(advxb + Dif_idx(i) - 1)%sf(k + offsets(1), l + offsets(2), q + offsets(3)) + 0.5_wp*alpha_nonconserv(i)
                        end if

                        if (k > -1 .and. l > -1 .and. q > -1) then
                            rhs_vf(advxb + Dif_idx(i) - 1)%sf(k, l, q) = &
                                rhs_vf(advxb + Dif_idx(i) - 1)%sf(k, l, q) + 0.5_wp*alpha_nonconserv(i)
                        end if
                    end do
                    if (k > -1 .and. l > -1 .and. q > -1) then
                        rhs_vf(momxb)%sf(k, l, q) = rhs_vf(momxb)%sf(k, l, q) - sigma*rho_L*q_prim_vf(momxb)%sf(k, l, q)
                    end if
                end do
            end do
        end do
        
        ! compute kdivu term
        !if (num_fluids > Dif_size) then

            !$acc parallel loop collapse(3) gang vector default(present)
            !do q = 0, p
                !do l = 0, n
                    !do k = 0, m
                        !dvel_dx(k, l, q) = 0._wp
                        !do r = -fd_number, fd_number
                            !dvel_dx(k, l, q) = dvel_dx(k, l, q) + q_prim_vf(momxb + idir - 1)%sf(k + r, l, q)*fd_coeff_x_d(r, k)
                        !end do
                    !end do
                !end do
            !end do
            !$acc end parallel loop

            !$acc parallel loop collapse(3) gang vector default(present)
            !do q = 0, p
                !do l = 0, n
                    !do k = 0, m
                        !denom(k, l, q) = 0._wp
                        !if (alpha_dif(k, l, q) > small_number) then
                            !do i = 1, Dif_size
                                !denom(k, l, q) = denom(k, l, q) + alpha_K_dif(k, l, q, i) / ( q_prim_vf(E_idx)%sf(k, l, q)*(1._wp + gammas(Dif_idx(i))) / gammas(Dif_idx(i)) )
                            !end do
                        
                            !rho1c12(k, l, q) = ( (gammas(liq_idx) + 1._wp)*q_prim_vf(E_idx)%sf(k, l, q) + pi_infs(liq_idx) )/gammas(liq_idx)
                            !rhogcg2(k, l, q) = 1._wp / denom(k, l, q)
                            !kdivu(k, l, q) = alpha_dif(k, l, q)*(1._wp - alpha_dif(k, l, q))*(rhogcg2(k, l, q) - rho1c12(k, l, q)) / ( (1._wp - alpha_dif(k, l, q))*rhogcg2(k, l, q) + alpha_dif(k, l, q)*rho1c12(k, l, q) )
                        !end if
                    !end do
                !end do
            !end do
            !$acc end parallel loop

            !$acc parallel loop collapse(3) gang vector default(present)
            !do q = 0, p
                !do l = 0, n
                    !do k = 0, m
                        !if (alpha_dif(k, l, q) > small_number) then
                            !do i = 1, Dif_size
                                !rhs_vf(advxb + Dif_idx(i) - 1)%sf(k, l, q) = rhs_vf(advxb + Dif_idx(i) - 1)%sf(k, l, q) - alpha_K_dif(k, l, q, i) / alpha_dif(k, l, q) &
                                    !* kdivu(k, l, q)*dvel_dx(k, l, q)
                            !end do
                            !rhs_vf(advxb + liq_idx - 1)%sf(k, l, q) = rhs_vf(advxb + liq_idx - 1)%sf(k, l, q) + kdivu(k, l, q)*dvel_dx(k, l, q)
                        !end if
                    !end do
                !end do
            !end do
            !$acc end parallel loop
        !end if
        
        
        

    end subroutine s_compute_diffusion_rhs

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