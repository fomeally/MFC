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

    private; public :: s_initialize_diffusion_module, &
s_compute_diffusion_rhs, &
s_finalize_diffusion_module

    real(wp), allocatable, dimension(:, :) :: fd_coeff_x_d
    real(wp), allocatable, dimension(:, :) :: fd_coeff_y_d
    real(wp), allocatable, dimension(:, :) :: fd_coeff_z_d
    !$acc declare create(fd_coeff_x_d,fd_coeff_y_d,fd_coeff_z_d)

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

        integer :: i !< generic loop iterators
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

    subroutine s_compute_diffusion_rhs(idir, j_prim_vf, q_prim_vf, rhs_vf)

        integer, intent(in) :: idir
        !type(scalar_field), dimension(sys_size), intent(in) :: dj_prim_dx_qp, dj_prim_dy_qp, dj_prim_dz_qp
        type(scalar_field), dimension(sys_size), intent(inout) :: j_prim_vf, q_prim_vf, rhs_vf

        integer :: i, k, l, q, r !< Loop variables
        real(wp), dimension(2) :: dif_flg
        real(wp) :: W1, W2, W3, D12, D13, D23
        real(wp) :: R_univ, small_number

        dif_flg(1) = 1._wp; dif_flg(2) = -1._wp
        R_univ = 8314.3_wp

        small_number = 1.0e-8_wp


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

            !print *, "261"
            do q = 0, p
                do l = 0, n
                    do k = -2*fd_number, m + 2*fd_number
                        ! gas cell
                        if (alpha_dif(k, l, q) > small_number) then
                            T_dif(k, l, q) = q_prim_vf(E_idx)%sf(k, l, q) * W_dif(k, l, q) /( rho_dif(k, l, q)*R_univ )
                        else
                            T_dif(k, l, q) = 0._wp
                        end if
                    end do
                end do
            end do
        
            !print *, "297"
            !$acc parallel loop collapse(4) gang vector default(present)
            do q = 0, p
                do l = 0, n
                    do k = -2*fd_number, m + 2*fd_number
                        if (alpha_dif(k, l, q) > small_number) then
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
            !print *, "311"
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
            !print *, "325"
            !set ghost cell for buffer region equal to k point (to enforce del dot j = 0)
            !$acc parallel loop collapse(5) gang vector default(present)
            do q = 0, p
                do l = 0, n
                    do k = -fd_number, m + fd_number
                        if (alpha_dif(k, l, q) > small_number) then
                            do i = 1, Dif_size
                                do r = -fd_number, fd_number
                                    if (alpha_dif(k + r, l, q) > small_number) then
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
            !print *, "346"
            if (Dif_size == 2) then
                !$acc parallel loop collapse(5) gang vector default(present)
                do q = 0, p
                    do l = 0, n
                        do k = 0, m
                            if (alpha_dif(k, l, q) > small_number) then
                                do i = 1, Dif_size
                                    do r = -fd_number, fd_number
                                        if (alpha_dif(k + r, l, q) > small_number) then
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
                            if (alpha_dif(k, l, q) > small_number) then
                                do r = -fd_number, fd_number
                                    do i = 1, Dif_size
                                        if (alpha_dif(k + r, l, q) > small_number) then
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
            !print *, "460"
            !Valid for any number of species
            ! species continuity
            !$acc parallel loop collapse(4) gang vector default(present)
            do q = 0, p
                do l = 0, n
                    do k = 0, m
                        if (alpha_dif(k, l, q) > small_number) then
                            do i = 1, Dif_size
                                rhs_vf(Dif_idx(i))%sf(k, l, q) = rhs_vf(Dif_idx(i))%sf(k, l, q) &
                                    + dj_dx(k, l, q, i)
                            end do
                        end if
                    end do
                end do
            end do
            !$acc end parallel loop
            !print *, "475"
            !Valid for any number of species
            !energy
            !$acc parallel loop collapse(4) gang vector default(present)
            do q = 0, p
                do l = 0, n
                    do k = 0, m
                        if (alpha_dif(k, l, q) > small_number) then
                            do i = 1, Dif_size
                                rhs_vf(E_idx)%sf(k, l, q) = rhs_vf(E_idx)%sf(k, l, q) &
                                    + djh_dx(k, l, q, i)
                            end do
                        end if
                    end do
                end do
            end do
            !$acc end parallel loop
            !print *, "490"
            !volume fraction
            if (Dif_size == 2) then
                !$acc parallel loop collapse(3) gang vector default(present)
                do q = 0, p
                    do l = 0, n
                        do k = 0, m
                            if (alpha_dif(k, l, q) > small_number) then
                                do i = 1, Dif_size
                                    rhs_vf(advxb + Dif_idx(i) - 1)%sf(k, l, q) = rhs_vf(advxb + Dif_idx(i) - 1)%sf(k, l, q) &
                                        + dj_dx(k, l, q, i)*W_dif(k, l, q)*W_dif(k, l, q) / ( rho_dif(k, l, q)*Ws(1)*Ws(2)*alpha_dif(k, l, q))
                                end do
                            end if
                        end do
                    end do
                end do
                !$acc end parallel loop
            elseif (Dif_size == 3) then
                !$acc parallel loop collapse(4) gang vector default(present)
                do q = 0, p
                    do l = 0, n
                        do k = 0, m
                            if (alpha_dif(k, l, q) > small_number) then
                                do i = 1, Dif_size
                                    select case (i)
                                        case (1)
                                            rhs_vf(advxb + Dif_idx(i) - 1)%sf(k, l, q) = rhs_vf(advxb + Dif_idx(i) - 1)%sf(k, l, q) &
                                                + W_dif(k, l, q) / ( Ws(1)*Ws(2)*Ws(3)*rho_dif(k, l, q)*alpha_dif(k, l, q) ) &
                                            * ( dj_dx(k, l, q, 1)*(alpha_dif(k, l, q) - alpha_K_dif(k, l, q, 1))*Ws(2)*Ws(3) - alpha_K_dif(k, l, q, 1)*Ws(1)*( Ws(2)*dj_dx(k, l, q, 3) + Ws(3)*dj_dx(k, l, q, 2) ) )

                                        case (2)
                                            rhs_vf(advxb + Dif_idx(i) - 1)%sf(k, l, q) = rhs_vf(advxb + Dif_idx(i) - 1)%sf(k, l, q) &
                                                + W_dif(k, l, q) / ( Ws(1)*Ws(2)*Ws(3)*rho_dif(k, l, q)*alpha_dif(k, l, q) ) &
                                                * ( dj_dx(k, l, q, 2)*(alpha_dif(k, l, q) - alpha_K_dif(k, l, q, 2))*Ws(1)*Ws(3) - alpha_K_dif(k, l, q, 2)*Ws(2)*( Ws(1)*dj_dx(k, l, q, 3) + Ws(3)*dj_dx(k, l, q, 1) ) )

                                        case (3)
                                            rhs_vf(advxb + Dif_idx(i) - 1)%sf(k, l, q) = rhs_vf(advxb + Dif_idx(i) - 1)%sf(k, l, q) &
                                                + W_dif(k, l, q) / ( Ws(1)*Ws(2)*Ws(3)*rho_dif(k, l, q)*alpha_dif(k, l, q) ) &
                                                * ( dj_dx(k, l, q, 3)*(alpha_dif(k, l, q) - alpha_K_dif(k, l, q, 3))*Ws(1)*Ws(2) - alpha_K_dif(k, l, q, 3)*Ws(3)*( Ws(1)*dj_dx(k, l, q, 2) + Ws(2)*dj_dx(k, l, q, 1) ) )
                                    end select
                                end do
                            end if
                        end do
                    end do
                end do
                !$acc end parallel loop
            end if
            !print *, "533"
            !print *, liq_idx
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
        end if
        
        

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