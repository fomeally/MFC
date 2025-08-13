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
    $:GPU_DECLARE(create='[fd_coeff_x_d,fd_coeff_y_d,fd_coeff_z_d]')

    real(wp), allocatable, dimension(:, :) :: Ds
    $:GPU_DECLARE(create='[Ds]')

    real(wp), allocatable, dimension(:) :: Ws
    $:GPU_DECLARE(create='[Ws]')

    real(wp), allocatable, dimension(:, :, :, :) :: dj_dx, dj_dy, dj_dz, djh_dx, djh_dy, djh_dz, dY_dx, dY_dy, dY_dz, alpha_K_dif, alpharho_K_dif, dF_KdP, F_K_dif
    $:GPU_DECLARE(create='[dj_dx, dj_dy, dj_dz, djh_dx, djh_dy, djh_dz, dY_dx, dY_dy, dY_dz, alpha_K_dif, alpharho_K_dif, dF_KdP, F_K_dif]')

    real(wp), allocatable, dimension(:, :, :) :: Gamma_dif, Pi_inf_dif, F_dif, dYda, dYdP, dPdt, rho_dif, dvel_dx, dvel_dy, dvel_dz
    $:GPU_DECLARE(create='[Gamma_dif, Pi_inf_dif, F_dif, dYda, dYdP, dPdt, rho_dif, dvel_dx, dvel_dy, dvel_dz]')

contains

    subroutine s_initialize_diffusion_module

        integer :: i !< generic loop iterators
        integer :: m_end, n_end, p_end
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

        @:ALLOCATE(Ds(1:num_fluids, 1:num_fluids))

        $:GPU_LOOP(parallelism='[seq]')
        do i = 1, num_fluids
            do j = 1, num_fluids
                Ds(i, j) = fluid_pp(i)%D(j)
            end do
        end do

        $:GPU_UPDATE(device='[Ds]')


        @:ALLOCATE(Ws(1:num_fluids))

        $:GPU_LOOP(parallelism='[seq]')
        do i = 1, num_fluids
            Ws(i) = fluid_pp(i)%W
        end do
        $:GPU_UPDATE(device='[Ws]')

        @:ALLOCATE(dj_dx(-fd_number:m_end, 0:n, 0:p, 1:num_fluids))
        @:ALLOCATE(djh_dx(-fd_number:m_end, 0:n, 0:p, 1:num_fluids))
        @:ALLOCATE(dY_dx(-fd_number:m_end, 0:n, 0:p, 1:num_fluids))
        @:ALLOCATE(dvel_dx(0:m, 0:n, 0:p))
        if (n > 0) then
            @:ALLOCATE(dj_dy(0:m, -fd_number:n_end, 0:p, 1:num_fluids))
            @:ALLOCATE(djh_dy(0:m, -fd_number:n_end, 0:p, 1:num_fluids))
            @:ALLOCATE(dY_dy(0:m, -fd_number:n_end, 0:p, 1:num_fluids))
            @:ALLOCATE(dvel_dy(0:m, 0:n, 0:p))
            if (p > 0) then
                @:ALLOCATE(dj_dz(0:m, 0:n, -fd_number:p_end, 1:num_fluids))
                @:ALLOCATE(djh_dz(0:m, 0:n, -fd_number:p_end, 1:num_fluids))
                @:ALLOCATE(dY_dz(0:m, 0:n, -fd_number:p_end, 1:num_fluids))
                @:ALLOCATE(dvel_dz(0:m, 0:n, 0:p))
            end if
        end if

        @:ALLOCATE(alpha_K_dif(-fd_number:m_end, -fd_number:n_end, -fd_number:p_end, 1:num_fluids))
        @:ALLOCATE(alpharho_K_dif(-fd_number:m_end, -fd_number:n_end, -fd_number:p_end, 1:num_fluids))
        @:ALLOCATE(dF_KdP(-fd_number:m_end, -fd_number:n_end, -fd_number:p_end, 1:num_fluids))
        @:ALLOCATE(F_K_dif(-fd_number:m_end, -fd_number:n_end, -fd_number:p_end, 1:num_fluids))
        @:ALLOCATE(Pi_inf_dif(-fd_number:m_end, -fd_number:n_end, -fd_number:p_end))
        @:ALLOCATE(Gamma_dif(-fd_number:m_end, -fd_number:n_end, -fd_number:p_end))
        @:ALLOCATE(F_dif(-fd_number:m_end, -fd_number:n_end, -fd_number:p_end))
        @:ALLOCATE(rho_dif(-fd_number:m_end, -fd_number:n_end, -fd_number:p_end))
        @:ALLOCATE(dYda(0:m, 0:n, 0:p))
        @:ALLOCATE(dYdP(0:m, 0:n, 0:p))
        @:ALLOCATE(dPdt(0:m, 0:n, 0:p))

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
        $:GPU_UPDATE(device='[fd_coeff_x_d]')
        if (n > 0) then
            call s_compute_finite_difference_coefficients(n, y_cc, fd_coeff_y_d, buff_size, &
                                                          fd_number, fd_order, offset_s(2))
            $:GPU_UPDATE(device='[fd_coeff_y_d]')
        end if
        if (p > 0) then
            call s_compute_finite_difference_coefficients(p, z_cc, fd_coeff_z_d, buff_size, &
                                                          fd_number, fd_order, offset_s(3))
            $:GPU_UPDATE(device='[fd_coeff_z_d]')
        end if

    end subroutine s_initialize_diffusion_module


    subroutine s_compute_diffusion_rhs(idir, j_prim_vf, q_prim_vf, rhs_vf)

        integer, intent(in) :: idir
        type(scalar_field), dimension(sys_size), intent(inout) :: j_prim_vf, q_prim_vf, rhs_vf

        integer :: i, k, l, q, r !< Loop variables
        real(wp), dimension(2) :: dif_flg
        real(wp) :: W1, W2, W3, D12, D13, D23, W_dif, W_fac, inv_denom
        dif_flg(1) = 1._wp; dif_flg(2) = -1._wp

        if (cyl_coord) then
            if (idir == 1) then
                $:GPU_PARALLEL_LOOP(collapse=3)
                do q = 0, p
                    do l = 0, n
                        do k = -fd_number, m + fd_number
                            do i = 1, num_fluids
                                alpha_K_dif(k, l, q, i) = q_prim_vf(E_idx + i)%sf(k, l, q)
                                alpharho_K_dif(k, l, q, i) = q_prim_vf(i)%sf(k, l, q)
                                F_K_dif(k, l, q, i) = ( q_prim_vf(E_idx)%sf(k, l, q)*gammas(i) + &
                                                        (pi_infs(i)*gammas(i) / ( 1._wp + gammas(i) )) ) / cvs(i)
                                dF_KdP(k, l, q, i) = gammas(i) / cvs(i)
                            end do
                        end do
                    end do
                end do

                $:GPU_PARALLEL_LOOP(collapse=3)
                do q = 0, p
                    do l = 0, n
                        do k = -fd_number, m + fd_number
                            do i = 1, num_fluids
                                Gamma_dif(k, l, q) = 0._wp
                                Pi_inf_dif(k, l, q) = 0._wp
                                F_dif(k, l, q) = 0._wp
                                rho_dif(k, l, q) = 0._wp
                            end do
                        end do
                    end do
                end do

                $:GPU_PARALLEL_LOOP(collapse=3)
                do q = 0, p
                    do l = 0, n
                        do k = -fd_number, m + fd_number
                            do i = 1, num_fluids
                                Gamma_dif(k, l, q) = Gamma_dif(k, l, q) + alpha_K_dif(k, l, q, i)*gammas(i)
                                Pi_inf_dif(k, l, q) = Pi_inf_dif(k, l, q) + alpha_K_dif(k, l, q, i)*pi_infs(i)
                                F_dif(k, l, q) = F_dif(k, l, q) + alpha_K_dif(k, l, q, i)*F_K_dif(k, l, q, i)
                                rho_dif(k, l, q) = rho_dif(k, l, q) + alpharho_K_dif(k, l, q, i)
                            end do
                        end do
                    end do
                end do

                ! For now, only computes for BINARY diffusion. Computed dPdt without viscous term and velocity divergence
                $:GPU_PARALLEL_LOOP(collapse=3)
                do q = 0, p
                    do l = 0, n
                        do k = 0, m
                            dYda(k, l, q) = F_K_dif(k, l, q, 1)*F_K_dif(k, l, q, 2) / (F_dif(k, l, q) ** 2._wp)
                            dYdP(k, l, q) = ( alpha_K_dif(k, l, q, 1)*alpha_k_dif(k, l, q, 2) ) * (F_K_dif(k, l, q, 2)*dF_KdP(k, l, q, 1) - &
                                                F_K_dif(k, l, q, 1)*dF_KdP(k, l, q, 2)) / (F_dif(k, l, q) ** 2._wp)
                            dPdt(k, l, q) = -( q_prim_vf(E_idx)%sf(k, l, q)*(Gamma_dif(k, l, q) + 1._wp) + Pi_inf_dif(k, l, q) ) / Gamma_dif(k, l, q)                   
                        end do
                    end do
                end do

                $:GPU_PARALLEL_LOOP(collapse=3)
                do q = 0, p
                    do l = 0, n
                        do k = -fd_number, m + fd_number
                            do i = 1, num_fluids
                                dj_dx(k, l, q, i) = 0._wp
                                djh_dx(k, l, q, i) = 0._wp
                                dY_dx(k, l, q, i) = 0._wp
                            end do
                        end do
                    end do
                end do

                $:GPU_PARALLEL_LOOP(collapse=3)
                do q = 0, p
                    do l = 0, n
                        do k = 0, m
                            dvel_dx(k, l, q) = 0._wp               
                        end do
                    end do
                end do


                $:GPU_PARALLEL_LOOP(collapse=3)
                do q = 0, p
                    do l = 0, n
                        do k = 0, m
                            $:GPU_LOOP(parallelism='[seq]')
                            do r = -fd_number, fd_number
                                dvel_dx(k, l, q) = dvel_dx(k, l, q) + q_prim_vf(momxb + idir - 1)%sf(k + r, l, q)*fd_coeff_x_d(r, k)               
                            end do               
                        end do
                    end do
                end do

                $:GPU_PARALLEL_LOOP(collapse=3)
                do q = 0, p
                    do l = 0, n
                        do k = -fd_number, m + fd_number
                            do i = 1, num_fluids
                                $:GPU_LOOP(parallelism='[seq]')
                                do r = -fd_number, fd_number
                                    dY_dx(k, l, q, i) = dY_dx(k, l, q, i) &
                                        + j_prim_vf(i)%sf(k + r, l, q)*fd_coeff_x_d(r, k)
                                end do
                            end do
                        end do
                    end do
                end do

                $:GPU_PARALLEL_LOOP(collapse=3)
                do q = 0, p
                    do l = 0, n
                        do k = 0, m
                            do i = 1, num_fluids
                                $:GPU_LOOP(parallelism='[seq]')
                                do r = -fd_number, fd_number
                                    dj_dx(k, l, q, i) = dj_dx(k, l, q, i) &
                                        + dY_dx(k + r, l, q, i)*rho_dif(k + r, l, q)*Ds(1, 2)*fd_coeff_x_d(r, k)
                                    djh_dx(k, l, q, i) = djh_dx(k, l, q, i) &
                                        + dY_dx(k + r, l, q, i)*rho_dif(k + r, l, q)*j_prim_vf(advxb + i - 1)%sf(k + r, l, q)*Ds(1, 2)*fd_coeff_x_d(r, k)
                                end do
                            end do
                        end do
                    end do
                end do

                !Valid for any number of species
                ! species continuity
                $:GPU_PARALLEL_LOOP(collapse=3)
                do q = 0, p
                    do l = 0, n
                        do k = 0, m
                            do i = 1, num_fluids                           
                                rhs_vf(i)%sf(k, l, q) = rhs_vf(i)%sf(k, l, q) &
                                    + dj_dx(k, l, q, i)
                            end do
                        end do
                    end do
                end do

                !Only valid for binary diffusion
                !volume fraction
                $:GPU_PARALLEL_LOOP(collapse=3)
                do q = 0, p
                    do l = 0, n
                        do k = 0, m
                            do i = 1, num_fluids
                                rhs_vf(advxb + i - 1)%sf(k, l, q) = rhs_vf(advxb + i - 1)%sf(k, l, q) &
                                    + ( dj_dx(k, l, q, i) / rho_dif(k, l, q) - dif_flg(i)*dYdP(k, l, q)*dPdt(k, l, q)*dvel_dx(k, l, q) )/ dYda(k, l, q)
                            end do
                        end do
                    end do
                end do

                !Valid for any number of species
                !energy
                $:GPU_PARALLEL_LOOP(collapse=3)
                do q = 0, p
                    do l = 0, n
                        do k = 0, m
                            do i = 1, num_fluids
                                rhs_vf(E_idx)%sf(k, l, q) = rhs_vf(E_idx)%sf(k, l, q) &
                                    + djh_dx(k, l, q, i)
                            end do
                        end do
                    end do
                end do
            elseif (idir == 2) then
                $:GPU_PARALLEL_LOOP(collapse=3)
                do q = 0, p
                    do l = -fd_number, n + fd_number
                        do k = 0, m
                            do i = 1, num_fluids
                                alpha_K_dif(k, l, q, i) = q_prim_vf(E_idx + i)%sf(k, l, q)
                                alpharho_K_dif(k, l, q, i) = q_prim_vf(i)%sf(k, l, q)
                                F_K_dif(k, l, q, i) = ( q_prim_vf(E_idx)%sf(k, l, q)*gammas(i) + &
                                                        (pi_infs(i)*gammas(i) / ( 1._wp + gammas(i) )) ) / cvs(i)
                                dF_KdP(k, l, q, i) = gammas(i) / cvs(i)
                            end do
                        end do
                    end do
                end do

                $:GPU_PARALLEL_LOOP(collapse=3)
                do q = 0, p
                    do l = -fd_number, n + fd_number
                        do k = 0, m
                            do i = 1, num_fluids
                                Gamma_dif(k, l, q) = 0._wp
                                Pi_inf_dif(k, l, q) = 0._wp
                                F_dif(k, l, q) = 0._wp
                                rho_dif(k, l, q) = 0._wp
                            end do
                        end do
                    end do
                end do

                $:GPU_PARALLEL_LOOP(collapse=3)
                do q = 0, p
                    do l = -fd_number, n + fd_number
                        do k = 0, m
                            do i = 1, num_fluids
                                Gamma_dif(k, l, q) = Gamma_dif(k, l, q) + alpha_K_dif(k, l, q, i)*gammas(i)
                                Pi_inf_dif(k, l, q) = Pi_inf_dif(k, l, q) + alpha_K_dif(k, l, q, i)*pi_infs(i)
                                F_dif(k, l, q) = F_dif(k, l, q) + alpha_K_dif(k, l, q, i)*F_K_dif(k, l, q, i)
                                rho_dif(k, l, q) = rho_dif(k, l, q) + alpharho_K_dif(k, l, q, i)
                            end do
                        end do
                    end do
                end do

                ! For now, only computes for BINARY diffusion. Computed dPdt without viscous term and velocity divergence
                $:GPU_PARALLEL_LOOP(collapse=3)
                do q = 0, p
                    do l = 0, n
                        do k = 0, m
                            dYda(k, l, q) = F_K_dif(k, l, q, 1)*F_K_dif(k, l, q, 2) / (F_dif(k, l, q) ** 2._wp)
                            dYdP(k, l, q) = ( alpha_K_dif(k, l, q, 1)*alpha_k_dif(k, l, q, 2) ) * (F_K_dif(k, l, q, 2)*dF_KdP(k, l, q, 1) - &
                                                F_K_dif(k, l, q, 1)*dF_KdP(k, l, q, 2)) / (F_dif(k, l, q) ** 2._wp)
                            dPdt(k, l, q) = -( q_prim_vf(E_idx)%sf(k, l, q)*(Gamma_dif(k, l, q) + 1._wp) + Pi_inf_dif(k, l, q) ) / Gamma_dif(k, l, q)                   
                        end do
                    end do
                end do

                $:GPU_PARALLEL_LOOP(collapse=3)
                do q = 0, p
                    do l = -fd_number, n + fd_number
                        do k = 0, m
                            do i = 1, num_fluids
                                dj_dy(k, l, q, i) = 0._wp
                                djh_dy(k, l, q, i) = 0._wp
                                dY_dy(k, l, q, i) = 0._wp
                            end do
                        end do
                    end do
                end do

                $:GPU_PARALLEL_LOOP(collapse=3)
                do q = 0, p
                    do l = 0, n
                        do k = 0, m
                            dvel_dy(k, l, q) = 0._wp               
                        end do
                    end do
                end do


                $:GPU_PARALLEL_LOOP(collapse=3)
                do q = 0, p
                    do l = 0, n
                        do k = 0, m
                            $:GPU_LOOP(parallelism='[seq]')
                            do r = -fd_number, fd_number
                                dvel_dy(k, l, q) = dvel_dy(k, l, q) + q_prim_vf(momxb + idir - 1)%sf(k, l + r, q)*fd_coeff_y_d(r, l) + &
                                    q_prim_vf(momyb + idir - 1)%sf(k, l, q) / y_cc(l)         
                            end do               
                        end do
                    end do
                end do
                !$acc end parallel loop

                $:GPU_PARALLEL_LOOP(collapse=3)
                do q = 0, p
                    do l = -fd_number, n + fd_number
                        do k = 0, m
                            do i = 1, num_fluids
                                $:GPU_LOOP(parallelism='[seq]')
                                do r = -fd_number, fd_number
                                    dY_dy(k, l, q, i) = dY_dy(k, l, q, i) &
                                        + j_prim_vf(i)%sf(k, l + r, q)*fd_coeff_y_d(r, l)
                                end do
                            end do
                        end do
                    end do
                end do

                $:GPU_PARALLEL_LOOP(collapse=3)
                do q = 0, p
                    do l = 0, n
                        do k = 0, m
                            do i = 1, num_fluids
                                $:GPU_LOOP(parallelism='[seq]')
                                do r = -fd_number, fd_number
                                    dj_dy(k, l, q, i) = dj_dy(k, l, q, i) &
                                        + dY_dy(k, l + r, q, i)*rho_dif(k, l + r, q)*Ds(1, 2)*fd_coeff_y_d(r, l)
                                    djh_dy(k, l, q, i) = djh_dy(k, l, q, i) &
                                        + dY_dy(k, l + r, q, i)*rho_dif(k, l + r, q)*j_prim_vf(advxb + i - 1)%sf(k, l + r, q)*Ds(1, 2)*fd_coeff_y_d(r, l)
                                end do
                            end do
                        end do
                    end do
                end do

                !Valid for any number of species
                ! species continuity
                $:GPU_PARALLEL_LOOP(collapse=3)
                do q = 0, p
                    do l = 0, n
                        do k = 0, m
                            do i = 1, num_fluids                           
                                rhs_vf(i)%sf(k, l, q) = rhs_vf(i)%sf(k, l, q) &
                                    + dj_dy(k, l, q, i) + rho_dif(k, l, q)*Ds(1, 2)*dY_dy(k, l, q, i) / y_cc(l)
                            end do
                        end do
                    end do
                end do

                !Only valid for binary diffusion
                !volume fraction
                $:GPU_PARALLEL_LOOP(collapse=3)
                do q = 0, p
                    do l = 0, n
                        do k = 0, m
                            do i = 1, num_fluids
                                rhs_vf(advxb + i - 1)%sf(k, l, q) = rhs_vf(advxb + i - 1)%sf(k, l, q) &
                                    + ( (dj_dy(k, l, q, i) + rho_dif(k, l, q)*Ds(1, 2)*dY_dy(k, l, q, i) / y_cc(l)) / rho_dif(k, l, q) - dif_flg(i)*dYdP(k, l, q)*dPdt(k, l, q)*dvel_dy(k, l, q) )/ dYda(k, l, q)
                            end do
                        end do
                    end do
                end do

                !Valid for any number of species
                !energy
                $:GPU_PARALLEL_LOOP(collapse=3)
                do q = 0, p
                    do l = 0, n
                        do k = 0, m
                            do i = 1, num_fluids
                                rhs_vf(E_idx)%sf(k, l, q) = rhs_vf(E_idx)%sf(k, l, q) &
                                    + djh_dy(k, l, q, i) + rho_dif(k, l, q)*Ds(1, 2)*dY_dy(k, l, q, i) / y_cc(l)
                            end do
                        end do
                    end do
                end do    

            elseif (idir == 3) then
                $:GPU_PARALLEL_LOOP(collapse=3)
                do q = -fd_number, p + fd_number
                    do l = 0, n
                        do k = 0, m
                            do i = 1, num_fluids
                                alpha_K_dif(k, l, q, i) = q_prim_vf(E_idx + i)%sf(k, l, q)
                                alpharho_K_dif(k, l, q, i) = q_prim_vf(i)%sf(k, l, q)
                                F_K_dif(k, l, q, i) = ( q_prim_vf(E_idx)%sf(k, l, q)*gammas(i) + &
                                                        (pi_infs(i)*gammas(i) / ( 1._wp + gammas(i) )) ) / cvs(i)
                                dF_KdP(k, l, q, i) = gammas(i) / cvs(i)
                            end do
                        end do
                    end do
                end do

                $:GPU_PARALLEL_LOOP(collapse=3)
                do q = -fd_number, p + fd_number
                    do l = 0, n
                        do k = 0, m
                            do i = 1, num_fluids
                                Gamma_dif(k, l, q) = 0._wp
                                Pi_inf_dif(k, l, q) = 0._wp
                                F_dif(k, l, q) = 0._wp
                                rho_dif(k, l, q) = 0._wp
                            end do
                        end do
                    end do
                end do

                $:GPU_PARALLEL_LOOP(collapse=3)
                do q = -fd_number, p + fd_number
                    do l = 0, n
                        do k = 0, m
                            do i = 1, num_fluids
                                Gamma_dif(k, l, q) = Gamma_dif(k, l, q) + alpha_K_dif(k, l, q, i)*gammas(i)
                                Pi_inf_dif(k, l, q) = Pi_inf_dif(k, l, q) + alpha_K_dif(k, l, q, i)*pi_infs(i)
                                F_dif(k, l, q) = F_dif(k, l, q) + alpha_K_dif(k, l, q, i)*F_K_dif(k, l, q, i)
                                rho_dif(k, l, q) = rho_dif(k, l, q) + alpharho_K_dif(k, l, q, i)
                            end do
                        end do
                    end do
                end do

                ! For now, only computes for BINARY diffusion. Computed dPdt without viscous term and velocity divergence
                $:GPU_PARALLEL_LOOP(collapse=3)
                do q = 0, p
                    do l = 0, n
                        do k = 0, m
                            dYda(k, l, q) = F_K_dif(k, l, q, 1)*F_K_dif(k, l, q, 2) / (F_dif(k, l, q) ** 2._wp)
                            dYdP(k, l, q) = ( alpha_K_dif(k, l, q, 1)*alpha_k_dif(k, l, q, 2) ) * (F_K_dif(k, l, q, 2)*dF_KdP(k, l, q, 1) - &
                                                F_K_dif(k, l, q, 1)*dF_KdP(k, l, q, 2)) / (F_dif(k, l, q) ** 2._wp)
                            dPdt(k, l, q) = -( q_prim_vf(E_idx)%sf(k, l, q)*(Gamma_dif(k, l, q) + 1._wp) + Pi_inf_dif(k, l, q) ) / Gamma_dif(k, l, q)                   
                        end do
                    end do
                end do

                $:GPU_PARALLEL_LOOP(collapse=3)
                do q = -fd_number, p + fd_number
                    do l = 0, n
                        do k = 0, m
                            do i = 1, num_fluids
                                dj_dz(k, l, q, i) = 0._wp
                                djh_dz(k, l, q, i) = 0._wp
                                dY_dz(k, l, q, i) = 0._wp
                            end do
                        end do
                    end do
                end do

                $:GPU_PARALLEL_LOOP(collapse=3)
                do q = 0, p
                    do l = 0, n
                        do k = 0, m
                            dvel_dz(k, l, q) = 0._wp               
                        end do
                    end do
                end do

                $:GPU_PARALLEL_LOOP(collapse=3)
                do q = 0, p
                    do l = 0, n
                        do k = 0, m
                            do r = -fd_number, fd_number
                                dvel_dz(k, l, q) = dvel_dz(k, l, q) + q_prim_vf(momxb + idir - 1)%sf(k, l, q + r)*fd_coeff_z_d(r, q) / y_cc(l)               
                            end do               
                        end do
                    end do
                end do

                $:GPU_PARALLEL_LOOP(collapse=3)
                do q = -fd_number, p + fd_number
                    do l = 0, n
                        do k = 0, m
                            do i = 1, num_fluids
                                $:GPU_LOOP(parallelism='[seq]')
                                do r = -fd_number, fd_number
                                    dY_dz(k, l, q, i) = dY_dz(k, l, q, i) &
                                        + j_prim_vf(i)%sf(k, l, q + r)*fd_coeff_z_d(r, q) / y_cc(l)
                                end do
                            end do
                        end do
                    end do
                end do

                $:GPU_PARALLEL_LOOP(collapse=3)
                do q = 0, p
                    do l = 0, n
                        do k = 0, m
                            do i = 1, num_fluids
                                $:GPU_LOOP(parallelism='[seq]')
                                do r = -fd_number, fd_number
                                    dj_dz(k, l, q, i) = dj_dz(k, l, q, i) &
                                        + dY_dz(k, l, q + r, i)*rho_dif(k, l, q + r)*Ds(1, 2)*fd_coeff_z_d(r, p)
                                    djh_dz(k, l, q, i) = djh_dz(k, l, q, i) &
                                        + dY_dz(k, l, q + r, i)*rho_dif(k, l, q + r)*j_prim_vf(advxb + i - 1)%sf(k, l, q + r)*Ds(1, 2)*fd_coeff_z_d(r, p)
                                end do
                            end do
                        end do
                    end do
                end do

                !Valid for any number of species
                ! species continuity
                $:GPU_PARALLEL_LOOP(collapse=3)
                do q = 0, p
                    do l = 0, n
                        do k = 0, m
                            do i = 1, num_fluids
                                rhs_vf(i)%sf(k, l, q) = rhs_vf(i)%sf(k, l, q) &
                                    + dj_dz(k, l, q, i) / y_cc(l)
                            end do
                        end do
                    end do
                end do

                !Only valid for binary diffusion
                !volume fraction
                $:GPU_PARALLEL_LOOP(collapse=3)
                do q = 0, p
                    do l = 0, n
                        do k = 0, m
                            do i = 1, num_fluids
                                rhs_vf(advxb + i - 1)%sf(k, l, q) = rhs_vf(advxb + i - 1)%sf(k, l, q) &
                                    + ( ( dj_dz(k, l, q, i) / y_cc(l) )/ rho_dif(k, l, q) - dif_flg(i)*dYdP(k, l, q)*dPdt(k, l, q)*dvel_dz(k, l, q) )/ dYda(k, l, q)
                            end do
                        end do
                    end do
                end do

                !Valid for any number of species
                !energy
                $:GPU_PARALLEL_LOOP(collapse=3)
                do q = 0, p
                    do l = 0, n
                        do k = 0, m
                            do i = 1, num_fluids
                                rhs_vf(E_idx)%sf(k, l, q) = rhs_vf(E_idx)%sf(k, l, q) &
                                    + djh_dz(k, l, q, i) / y_cc(l)
                            end do
                        end do
                    end do
                end do


            end if

        else !cartesian coordinates
            if (num_fluids == 2) then
                if (idir == 1) then
                    $:GPU_PARALLEL_LOOP(collapse=3)
                    do q = 0, p
                        do l = 0, n
                            do k = -fd_number, m + fd_number
                                do i = 1, num_fluids
                                    alpha_K_dif(k, l, q, i) = q_prim_vf(E_idx + i)%sf(k, l, q)
                                    alpharho_K_dif(k, l, q, i) = q_prim_vf(i)%sf(k, l, q)
                                    F_K_dif(k, l, q, i) = ( q_prim_vf(E_idx)%sf(k, l, q)*gammas(i) + &
                                                            (pi_infs(i)*gammas(i) / ( 1._wp + gammas(i) )) ) / cvs(i)
                                    dF_KdP(k, l, q, i) = gammas(i) / cvs(i)
                                end do
                            end do
                        end do
                    end do

                    $:GPU_PARALLEL_LOOP(collapse=4)
                    do q = 0, p
                        do l = 0, n
                            do k = -fd_number, m + fd_number
                                do i = 1, num_fluids
                                    Gamma_dif(k, l, q) = 0._wp
                                    Pi_inf_dif(k, l, q) = 0._wp
                                    F_dif(k, l, q) = 0._wp
                                    rho_dif(k, l, q) = 0._wp
                                end do
                            end do
                        end do
                    end do

                    $:GPU_PARALLEL_LOOP(collapse=4)
                    do q = 0, p
                        do l = 0, n
                            do k = -fd_number, m + fd_number
                                do i = 1, num_fluids
                                    Gamma_dif(k, l, q) = Gamma_dif(k, l, q) + alpha_K_dif(k, l, q, i)*gammas(i)
                                    Pi_inf_dif(k, l, q) = Pi_inf_dif(k, l, q) + alpha_K_dif(k, l, q, i)*pi_infs(i)
                                    F_dif(k, l, q) = F_dif(k, l, q) + alpha_K_dif(k, l, q, i)*F_K_dif(k, l, q, i)
                                    rho_dif(k, l, q) = rho_dif(k, l, q) + alpharho_K_dif(k, l, q, i)
                                end do
                            end do
                        end do
                    end do

                    ! For now, only computes for BINARY diffusion. Computed dPdt without viscous term and velocity divergence
                    $:GPU_PARALLEL_LOOP(collapse=3)
                    do q = 0, p
                        do l = 0, n
                            do k = 0, m
                                dYda(k, l, q) = F_K_dif(k, l, q, 1)*F_K_dif(k, l, q, 2) / (F_dif(k, l, q) ** 2._wp)
                                dYdP(k, l, q) = ( alpha_K_dif(k, l, q, 1)*alpha_k_dif(k, l, q, 2) ) * (F_K_dif(k, l, q, 2)*dF_KdP(k, l, q, 1) - &
                                                    F_K_dif(k, l, q, 1)*dF_KdP(k, l, q, 2)) / (F_dif(k, l, q) ** 2._wp)
                                dPdt(k, l, q) = -( q_prim_vf(E_idx)%sf(k, l, q)*(Gamma_dif(k, l, q) + 1._wp) + Pi_inf_dif(k, l, q) ) / Gamma_dif(k, l, q)                   
                            end do
                        end do
                    end do

                    $:GPU_PARALLEL_LOOP(collapse=3)
                    do q = 0, p
                        do l = 0, n
                            do k = -fd_number, m + fd_number
                                do i = 1, num_fluids
                                    dj_dx(k, l, q, i) = 0._wp
                                    djh_dx(k, l, q, i) = 0._wp
                                    dY_dx(k, l, q, i) = 0._wp
                                end do
                            end do
                        end do
                    end do

                    $:GPU_PARALLEL_LOOP(collapse=3)
                    do q = 0, p
                        do l = 0, n
                            do k = 0, m
                                dvel_dx(k, l, q) = 0._wp               
                            end do
                        end do
                    end do


                    $:GPU_PARALLEL_LOOP(collapse=3)
                    do q = 0, p
                        do l = 0, n
                            do k = 0, m
                                $:GPU_LOOP(parallelism='[seq]')
                                do r = -fd_number, fd_number
                                    dvel_dx(k, l, q) = dvel_dx(k, l, q) + q_prim_vf(momxb + idir - 1)%sf(k + r, l, q)*fd_coeff_x_d(r, k)               
                                end do               
                            end do
                        end do
                    end do

                    $:GPU_PARALLEL_LOOP(collapse=3)
                    do q = 0, p
                        do l = 0, n
                            do k = -fd_number, m + fd_number
                                do i = 1, num_fluids
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do r = -fd_number, fd_number
                                        dY_dx(k, l, q, i) = dY_dx(k, l, q, i) &
                                            + j_prim_vf(i)%sf(k + r, l, q)*fd_coeff_x_d(r, k)
                                    end do
                                end do
                            end do
                        end do
                    end do

                    $:GPU_PARALLEL_LOOP(collapse=3)
                    do q = 0, p
                        do l = 0, n
                            do k = 0, m
                                do i = 1, num_fluids
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do r = -fd_number, fd_number
                                        dj_dx(k, l, q, i) = dj_dx(k, l, q, i) &
                                            + dY_dx(k + r, l, q, i)*rho_dif(k + r, l, q)*Ds(1, 2)*fd_coeff_x_d(r, k)
                                        djh_dx(k, l, q, i) = djh_dx(k, l, q, i) &
                                            + dY_dx(k + r, l, q, i)*rho_dif(k + r, l, q)*j_prim_vf(advxb + i - 1)%sf(k + r, l, q)*Ds(1, 2)*fd_coeff_x_d(r, k)
                                    end do
                                end do
                            end do
                        end do
                    end do

                    !Valid for any number of species
                    ! species continuity
                    $:GPU_PARALLEL_LOOP(collapse=3)
                    do q = 0, p
                        do l = 0, n
                            do k = 0, m
                                do i = 1, num_fluids
                                    rhs_vf(i)%sf(k, l, q) = rhs_vf(i)%sf(k, l, q) &
                                        + dj_dx(k, l, q, i)
                                end do
                            end do
                        end do
                    end do

                    !Only valid for binary diffusion
                    !volume fraction
                    $:GPU_PARALLEL_LOOP(collapse=3)
                    do q = 0, p
                        do l = 0, n
                            do k = 0, m
                                do i = 1, num_fluids

                                    rhs_vf(advxb + i - 1)%sf(k, l, q) = rhs_vf(advxb + i - 1)%sf(k, l, q) &
                                        + ( dj_dx(k, l, q, i) / rho_dif(k, l, q) - dif_flg(i)*dYdP(k, l, q)*dPdt(k, l, q)*dvel_dx(k, l, q) )/ dYda(k, l, q)
                                end do
                            end do
                        end do
                    end do

                    !Valid for any number of species
                    !energy
                    $:GPU_PARALLEL_LOOP(collapse=3)
                    do q = 0, p
                        do l = 0, n
                            do k = 0, m
                                do i = 1, num_fluids
                                    rhs_vf(E_idx)%sf(k, l, q) = rhs_vf(E_idx)%sf(k, l, q) &
                                        + djh_dx(k, l, q, i)
                                end do
                            end do
                        end do
                    end do
                elseif (idir == 2) then
                    $:GPU_PARALLEL_LOOP(collapse=3)
                    do q = 0, p
                        do l = -fd_number, n + fd_number
                            do k = 0, m
                                do i = 1, num_fluids
                                    alpha_K_dif(k, l, q, i) = q_prim_vf(E_idx + i)%sf(k, l, q)
                                    alpharho_K_dif(k, l, q, i) = q_prim_vf(i)%sf(k, l, q)
                                    F_K_dif(k, l, q, i) = ( q_prim_vf(E_idx)%sf(k, l, q)*gammas(i) + &
                                                            (pi_infs(i)*gammas(i) / ( 1._wp + gammas(i) )) ) / cvs(i)
                                    dF_KdP(k, l, q, i) = gammas(i) / cvs(i)
                                end do
                            end do
                        end do
                    end do

                    $:GPU_PARALLEL_LOOP(collapse=3)
                    do q = 0, p
                        do l = -fd_number, n + fd_number
                            do k = 0, m
                                do i = 1, num_fluids
                                    Gamma_dif(k, l, q) = 0._wp
                                    Pi_inf_dif(k, l, q) = 0._wp
                                    F_dif(k, l, q) = 0._wp
                                    rho_dif(k, l, q) = 0._wp
                                end do
                            end do
                        end do
                    end do
                    
                    $:GPU_PARALLEL_LOOP(collapse=3)
                    do q = 0, p
                        do l = -fd_number, n + fd_number
                            do k = 0, m
                                do i = 1, num_fluids
                                    Gamma_dif(k, l, q) = Gamma_dif(k, l, q) + alpha_K_dif(k, l, q, i)*gammas(i)
                                    Pi_inf_dif(k, l, q) = Pi_inf_dif(k, l, q) + alpha_K_dif(k, l, q, i)*pi_infs(i)
                                    F_dif(k, l, q) = F_dif(k, l, q) + alpha_K_dif(k, l, q, i)*F_K_dif(k, l, q, i)
                                    rho_dif(k, l, q) = rho_dif(k, l, q) + alpharho_K_dif(k, l, q, i)
                                end do
                            end do
                        end do
                    end do

                    ! For now, only computes for BINARY diffusion. Computed dPdt without viscous term and velocity divergence
                    $:GPU_PARALLEL_LOOP(collapse=3)
                    do q = 0, p
                        do l = 0, n
                            do k = 0, m
                                dYda(k, l, q) = F_K_dif(k, l, q, 1)*F_K_dif(k, l, q, 2) / (F_dif(k, l, q) ** 2._wp)
                                dYdP(k, l, q) = ( alpha_K_dif(k, l, q, 1)*alpha_k_dif(k, l, q, 2) ) * (F_K_dif(k, l, q, 2)*dF_KdP(k, l, q, 1) - &
                                                    F_K_dif(k, l, q, 1)*dF_KdP(k, l, q, 2)) / (F_dif(k, l, q) ** 2._wp)
                                dPdt(k, l, q) = -( q_prim_vf(E_idx)%sf(k, l, q)*(Gamma_dif(k, l, q) + 1._wp) + Pi_inf_dif(k, l, q) ) / Gamma_dif(k, l, q)                   
                            end do
                        end do
                    end do

                    $:GPU_PARALLEL_LOOP(collapse=3)
                    do q = 0, p
                        do l = -fd_number, n + fd_number
                            do k = 0, m
                                do i = 1, num_fluids
                                    dj_dy(k, l, q, i) = 0._wp
                                    djh_dy(k, l, q, i) = 0._wp
                                    dY_dy(k, l, q, i) = 0._wp
                                end do
                            end do
                        end do
                    end do

                    $:GPU_PARALLEL_LOOP(collapse=3)
                    do q = 0, p
                        do l = 0, n
                            do k = 0, m
                                dvel_dy(k, l, q) = 0._wp               
                            end do
                        end do
                    end do

                    $:GPU_PARALLEL_LOOP(collapse=3)
                    do q = 0, p
                        do l = 0, n
                            do k = 0, m
                                $:GPU_LOOP(parallelism='[seq]')
                                do r = -fd_number, fd_number
                                    dvel_dy(k, l, q) = dvel_dy(k, l, q) + q_prim_vf(momxb + idir - 1)%sf(k, l + r, q)*fd_coeff_y_d(r, l)               
                                end do               
                            end do
                        end do
                    end do

                    $:GPU_PARALLEL_LOOP(collapse=3)
                    do q = 0, p
                        do l = -fd_number, n + fd_number
                            do k = 0, m
                                do i = 1, num_fluids
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do r = -fd_number, fd_number
                                        dY_dy(k, l, q, i) = dY_dy(k, l, q, i) &
                                            + j_prim_vf(i)%sf(k, l + r, q)*fd_coeff_y_d(r, l)
                                    end do
                                end do
                            end do
                        end do
                    end do

                    $:GPU_PARALLEL_LOOP(collapse=3)
                    do q = 0, p
                        do l = 0, n
                            do k = 0, m
                                do i = 1, num_fluids
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do r = -fd_number, fd_number
                                        dj_dy(k, l, q, i) = dj_dy(k, l, q, i) &
                                            + dY_dy(k, l + r, q, i)*rho_dif(k, l + r, q)*Ds(1, 2)*fd_coeff_y_d(r, l)
                                        djh_dy(k, l, q, i) = djh_dy(k, l, q, i) &
                                            + dY_dy(k, l + r, q, i)*rho_dif(k, l + r, q)*j_prim_vf(advxb + i - 1)%sf(k, l + r, q)*Ds(1, 2)*fd_coeff_y_d(r, l)
                                    end do
                                end do
                            end do
                        end do
                    end do

                    !Valid for any number of species
                    ! species continuity
                    $:GPU_PARALLEL_LOOP(collapse=3)
                    do q = 0, p
                        do l = 0, n
                            do k = 0, m
                                do i = 1, num_fluids
                                    rhs_vf(i)%sf(k, l, q) = rhs_vf(i)%sf(k, l, q) &
                                        + dj_dy(k, l, q, i)
                                end do
                            end do
                        end do
                    end do

                    !Only valid for binary diffusion
                    !volume fraction
                    $:GPU_PARALLEL_LOOP(collapse=3)
                    do q = 0, p
                        do l = 0, n
                            do k = 0, m
                                do i = 1, num_fluids

                                    rhs_vf(advxb + i - 1)%sf(k, l, q) = rhs_vf(advxb + i - 1)%sf(k, l, q) &
                                        + ( dj_dy(k, l, q, i) / rho_dif(k, l, q) - dif_flg(i)*dYdP(k, l, q)*dPdt(k, l, q)*dvel_dy(k, l, q) )/ dYda(k, l, q)
                                end do
                            end do
                        end do
                    end do

                    !Valid for any number of species
                    !energy
                    $:GPU_PARALLEL_LOOP(collapse=3)
                    do q = 0, p
                        do l = 0, n
                            do k = 0, m
                                do i = 1, num_fluids
                                    rhs_vf(E_idx)%sf(k, l, q) = rhs_vf(E_idx)%sf(k, l, q) &
                                        + djh_dy(k, l, q, i)
                                end do
                            end do
                        end do
                    end do        
                elseif (idir == 3) then
                    $:GPU_PARALLEL_LOOP(collapse=3)
                    do q = -fd_number, p + fd_number
                        do l = 0, n
                            do k = 0, m
                                do i = 1, num_fluids
                                    alpha_K_dif(k, l, q, i) = q_prim_vf(E_idx + i)%sf(k, l, q)
                                    alpharho_K_dif(k, l, q, i) = q_prim_vf(i)%sf(k, l, q)
                                    F_K_dif(k, l, q, i) = ( q_prim_vf(E_idx)%sf(k, l, q)*gammas(i) + &
                                                            (pi_infs(i)*gammas(i) / ( 1._wp + gammas(i) )) ) / cvs(i)
                                    dF_KdP(k, l, q, i) = gammas(i) / cvs(i)
                                end do
                            end do
                        end do
                    end do

                    $:GPU_PARALLEL_LOOP(collapse=3)
                    do q = -fd_number, p + fd_number
                        do l = 0, n
                            do k = 0, m
                                do i = 1, num_fluids
                                    Gamma_dif(k, l, q) = 0._wp
                                    Pi_inf_dif(k, l, q) = 0._wp
                                    F_dif(k, l, q) = 0._wp
                                    rho_dif(k, l, q) = 0._wp
                                end do
                            end do
                        end do
                    end do

                    $:GPU_PARALLEL_LOOP(collapse=3)
                    do q = -fd_number, p + fd_number
                        do l = 0, n
                            do k = 0, m
                                do i = 1, num_fluids
                                    Gamma_dif(k, l, q) = Gamma_dif(k, l, q) + alpha_K_dif(k, l, q, i)*gammas(i)
                                    Pi_inf_dif(k, l, q) = Pi_inf_dif(k, l, q) + alpha_K_dif(k, l, q, i)*pi_infs(i)
                                    F_dif(k, l, q) = F_dif(k, l, q) + alpha_K_dif(k, l, q, i)*F_K_dif(k, l, q, i)
                                    rho_dif(k, l, q) = rho_dif(k, l, q) + alpharho_K_dif(k, l, q, i)
                                end do
                            end do
                        end do
                    end do

                    ! For now, only computes for BINARY diffusion. Computed dPdt without viscous term and velocity divergence
                    $:GPU_PARALLEL_LOOP(collapse=3)
                    do q = 0, p
                        do l = 0, n
                            do k = 0, m
                                dYda(k, l, q) = F_K_dif(k, l, q, 1)*F_K_dif(k, l, q, 2) / (F_dif(k, l, q) ** 2._wp)
                                dYdP(k, l, q) = ( alpha_K_dif(k, l, q, 1)*alpha_k_dif(k, l, q, 2) ) * (F_K_dif(k, l, q, 2)*dF_KdP(k, l, q, 1) - &
                                                    F_K_dif(k, l, q, 1)*dF_KdP(k, l, q, 2)) / (F_dif(k, l, q) ** 2._wp)
                                dPdt(k, l, q) = -( q_prim_vf(E_idx)%sf(k, l, q)*(Gamma_dif(k, l, q) + 1._wp) + Pi_inf_dif(k, l, q) ) / Gamma_dif(k, l, q)                   
                            end do
                        end do
                    end do

                    $:GPU_PARALLEL_LOOP(collapse=3)
                    do q = -fd_number, p + fd_number
                        do l = 0, n
                            do k = 0, m
                                do i = 1, num_fluids
                                    dj_dz(k, l, q, i) = 0._wp
                                    djh_dz(k, l, q, i) = 0._wp
                                    dY_dz(k, l, q, i) = 0._wp
                                end do
                            end do
                        end do
                    end do

                    $:GPU_PARALLEL_LOOP(collapse=3)
                    do q = 0, p
                        do l = 0, n
                            do k = 0, m
                                dvel_dz(k, l, q) = 0._wp               
                            end do
                        end do
                    end do


                    $:GPU_PARALLEL_LOOP(collapse=3)
                    do q = 0, p
                        do l = 0, n
                            do k = 0, m
                                $:GPU_LOOP(parallelism='[seq]')
                                do r = -fd_number, fd_number
                                    dvel_dz(k, l, q) = dvel_dz(k, l, q) + q_prim_vf(momxb + idir - 1)%sf(k, l, q + r)*fd_coeff_z_d(r, q)               
                                end do               
                            end do
                        end do
                    end do

                    $:GPU_PARALLEL_LOOP(collapse=3)
                    do q = -fd_number, p + fd_number
                        do l = 0, n
                            do k = 0, m
                                do i = 1, num_fluids
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do r = -fd_number, fd_number
                                        dY_dz(k, l, q, i) = dY_dz(k, l, q, i) &
                                            + j_prim_vf(i)%sf(k, l, q + r)*fd_coeff_z_d(r, q)
                                    end do
                                end do
                            end do
                        end do
                    end do

                    $:GPU_PARALLEL_LOOP(collapse=3)
                    do q = 0, p
                        do l = 0, n
                            do k = 0, m
                                do i = 1, num_fluids
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do r = -fd_number, fd_number
                                        dj_dz(k, l, q, i) = dj_dz(k, l, q, i) &
                                            + dY_dz(k, l, q + r, i)*rho_dif(k, l, q + r)*Ds(1, 2)*fd_coeff_z_d(r, p)
                                        djh_dz(k, l, q, i) = djh_dz(k, l, q, i) &
                                            + dY_dz(k, l, q + r, i)*rho_dif(k, l, q + r)*j_prim_vf(advxb + i - 1)%sf(k, l, q + r)*Ds(1, 2)*fd_coeff_z_d(r, p)
                                    end do
                                end do
                            end do
                        end do
                    end do

                    !Valid for any number of species
                    ! species continuity
                    $:GPU_PARALLEL_LOOP(collapse=3)
                    do q = 0, p
                        do l = 0, n
                            do k = 0, m
                                do i = 1, num_fluids
                                    rhs_vf(i)%sf(k, l, q) = rhs_vf(i)%sf(k, l, q) &
                                        + dj_dz(k, l, q, i)
                                end do
                            end do
                        end do
                    end do

                    !Only valid for binary diffusion
                    !volume fraction
                    $:GPU_PARALLEL_LOOP(collapse=3)
                    do q = 0, p
                        do l = 0, n
                            do k = 0, m
                                do i = 1, num_fluids
                                    rhs_vf(advxb + i - 1)%sf(k, l, q) = rhs_vf(advxb + i - 1)%sf(k, l, q) &
                                        + ( dj_dz(k, l, q, i) / rho_dif(k, l, q) - dif_flg(i)*dYdP(k, l, q)*dPdt(k, l, q)*dvel_dz(k, l, q) )/ dYda(k, l, q)
                                end do
                            end do
                        end do
                    end do

                    !Valid for any number of species
                    !energy
                    $:GPU_PARALLEL_LOOP(collapse=3)
                    do q = 0, p
                        do l = 0, n
                            do k = 0, m
                                do i = 1, num_fluids
                                    rhs_vf(E_idx)%sf(k, l, q) = rhs_vf(E_idx)%sf(k, l, q) &
                                        + djh_dz(k, l, q, i)
                                end do
                            end do
                        end do
                    end do
                end if
            elseif (num_fluids == 3) then
                ! For 3 fluids, we need to compute the diffusion terms for each pair of species
                ! 1 is N2, 2 is CO2, 3 is H2
                W1 = Ws(1)
                W2 = Ws(2)
                W3 = Ws(3)
                D12 = Ds(1,2)
                D13 = Ds(1,3)
                D23 = Ds(2,3)
                if (idir == 1) then
                    $:GPU_PARALLEL_LOOP(collapse=3)
                    do q = 0, p
                        do l = 0, n
                            do k = -fd_number, m + fd_number
                                do i = 1, num_fluids
                                    alpha_K_dif(k, l, q, i) = q_prim_vf(E_idx + i)%sf(k, l, q)
                                    alpharho_K_dif(k, l, q, i) = q_prim_vf(i)%sf(k, l, q)
                                end do
                            end do
                        end do
                    end do

                    $:GPU_PARALLEL_LOOP(collapse=3)
                    do q = 0, p
                        do l = 0, n
                            do k = -fd_number, m + fd_number
                                do i = 1, num_fluids
                                    rho_dif(k, l, q) = 0._wp
                                end do
                            end do
                        end do
                    end do


                    $:GPU_PARALLEL_LOOP(collapse=3)
                    do q = 0, p
                        do l = 0, n
                            do k = -fd_number, m + fd_number
                                do i = 1, num_fluids
                                    rho_dif(k, l, q) = rho_dif(k, l, q) + alpharho_K_dif(k, l, q, i)
                                end do
                            end do
                        end do
                    end do

                    $:GPU_PARALLEL_LOOP(collapse=3)
                    do q = 0, p
                        do l = 0, n
                            do k = -fd_number, m + fd_number
                                do i = 1, num_fluids
                                    dj_dx(k, l, q, i) = 0._wp
                                    djh_dx(k, l, q, i) = 0._wp
                                    dY_dx(k, l, q, i) = 0._wp
                                end do
                            end do
                        end do
                    end do

                    $:GPU_PARALLEL_LOOP(collapse=3)
                    do q = 0, p
                        do l = 0, n
                            do k = -fd_number, m + fd_number
                                do i = 1, num_fluids
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do r = -fd_number, fd_number
                                        dY_dx(k, l, q, i) = dY_dx(k, l, q, i) &
                                            + j_prim_vf(i)%sf(k + r, l, q)*fd_coeff_x_d(r, k)
                                    end do
                                end do
                            end do
                        end do
                    end do

                    $:GPU_PARALLEL_LOOP(collapse=3)
                    do q = 0, p
                        do l = 0, n
                            do k = 0, m
                                do r = -fd_number, fd_number
                                    inv_denom = 1._wp / ( j_prim_vf(1)%sf(k + r, l, q)*D23 + j_prim_vf(2)%sf(k + r, l, q)*D13 + j_prim_vf(3)%sf(k + r, l, q)*D12 )
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, num_fluids
                                        select case (i)
                                            case (1)
                                                dj_dx(k, l, q, i) = dj_dx(k, l, q, i) &
                                                    + rho_dif(k + r, l, q)*(D12*D13*(1._wp - j_prim_vf(1)%sf(k + r, l, q))*dY_dx(k + r, l, q, 1) &
                                                    - j_prim_vf(1)%sf(k + r, l, q)*D23*(D12*dY_dx(k + r, l, q, 2) + D13*dY_dx(k + r, l, q, 3))) &
                                                    * inv_denom*fd_coeff_x_d(r, k)

                                                djh_dx(k, l, q, i) = djh_dx(k, l, q, i) &
                                                    + j_prim_vf(advxb)%sf(k + r, l, q)*rho_dif(k + r, l, q)*(D12*D13*(1._wp - j_prim_vf(1)%sf(k + r, l, q))*dY_dx(k + r, l, q, 1) &
                                                    - j_prim_vf(1)%sf(k + r, l, q)*D23*(D12*dY_dx(k + r, l, q, 2) + D13*dY_dx(k + r, l, q, 3))) &
                                                    * inv_denom*fd_coeff_x_d(r, k)

                                            case (2)
                                                dj_dx(k, l, q, i) = dj_dx(k, l, q, i) &
                                                    + rho_dif(k + r, l, q)*(D12*D23*(1._wp - j_prim_vf(2)%sf(k + r, l, q))*dY_dx(k + r, l, q, 2) &
                                                    - j_prim_vf(2)%sf(k + r, l, q)*D13*(D12*dY_dx(k + r, l, q, 1) + D23*dY_dx(k + r, l, q, 3))) &
                                                    * inv_denom*fd_coeff_x_d(r, k)

                                                djh_dx(k, l, q, i) = djh_dx(k, l, q, i) &
                                                    + j_prim_vf(advxb + 1)%sf(k + r, l, q)*rho_dif(k + r, l, q)*(D12*D23*(1._wp - j_prim_vf(2)%sf(k + r, l, q))*dY_dx(k + r, l, q, 2) &
                                                    - j_prim_vf(2)%sf(k + r, l, q)*D13*(D12*dY_dx(k + r, l, q, 1) + D23*dY_dx(k + r, l, q, 3))) &
                                                    * inv_denom*fd_coeff_x_d(r, k)

                                            case (3)
                                                dj_dx(k, l, q, i) = dj_dx(k, l, q, i) &
                                                    + rho_dif(k + r, l, q)*(D23*D13*(1._wp - j_prim_vf(3)%sf(k + r, l, q))*dY_dx(k + r, l, q, 3) &
                                                    - j_prim_vf(3)%sf(k + r, l, q)*D12*(D13*dY_dx(k + r, l, q, 1) + D23*dY_dx(k + r, l, q, 2))) &
                                                    * inv_denom*fd_coeff_x_d(r, k)

                                                djh_dx(k, l, q, i) = djh_dx(k, l, q, i) &
                                                    + j_prim_vf(advxb + 2)%sf(k + r, l, q)*rho_dif(k + r, l, q)*(D23*D13*(1._wp - j_prim_vf(3)%sf(k + r, l, q))*dY_dx(k + r, l, q, 3) &
                                                    - j_prim_vf(3)%sf(k + r, l, q)*D12*(D13*dY_dx(k + r, l, q, 1) + D23*dY_dx(k + r, l, q, 2))) &
                                                    * inv_denom*fd_coeff_x_d(r, k)
                                        end select
                                    end do
                                    
                                end do 
                            end do
                        end do
                    end do

                    !Valid for any number of species
                    ! species continuity
                    $:GPU_PARALLEL_LOOP(collapse=3)
                    do q = 0, p
                        do l = 0, n
                            do k = 0, m
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_fluids
                                    rhs_vf(i)%sf(k, l, q) = rhs_vf(i)%sf(k, l, q) &
                                        + dj_dx(k, l, q, i)
                                end do
                            end do
                        end do
                    end do
                    

                    !volume fraction
                    !ideal gas only (3-component diffusion)
                    $:GPU_PARALLEL_LOOP(collapse=3)
                    do q = 0, p
                        do l = 0, n
                            do k = 0, m                        
                                W_dif = alpha_K_dif(k, l, q, 1)*W1 + alpha_K_dif(k, l, q, 2)*W2 + alpha_K_dif(k, l, q, 3)*W3
                                W_fac = W_dif / (W1*W2*W3*rho_dif(k, l, q))
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_fluids
                                    select case (i)
                                        case (1)
                                            rhs_vf(advxb + i - 1)%sf(k, l, q) = rhs_vf(advxb + i - 1)%sf(k, l, q) &
                                                + W_fac & 
                                                * ( dj_dx(k, l, q, 1)*(1._wp - alpha_K_dif(k, l, q, 1))*W2*W3 - alpha_K_dif(k, l, q, 1)*W1*( W2*dj_dx(k, l, q, 3) + W3*dj_dx(k, l, q, 2) ) )

                                        case (2)
                                            rhs_vf(advxb + i - 1)%sf(k, l, q) = rhs_vf(advxb + i - 1)%sf(k, l, q) &
                                                + W_fac &
                                                * ( dj_dx(k, l, q, 2)*(1._wp - alpha_K_dif(k, l, q, 2))*W1*W3 - alpha_K_dif(k, l, q, 2)*W2*( W1*dj_dx(k, l, q, 3) + W3*dj_dx(k, l, q, 1) ) )

                                        case (3)
                                            rhs_vf(advxb + i - 1)%sf(k, l, q) = rhs_vf(advxb + i - 1)%sf(k, l, q) &
                                                + W_fac &
                                                * ( dj_dx(k, l, q, 3)*(1._wp - alpha_K_dif(k, l, q, 3))*W1*W2 - alpha_K_dif(k, l, q, 3)*W3*( W1*dj_dx(k, l, q, 2) + W2*dj_dx(k, l, q, 1) ) )
                                    end select
                                end do
                            end do
                        end do
                    end do
                    

                    !Valid for any number of species
                    !energy
                    $:GPU_PARALLEL_LOOP(collapse=3)
                    do q = 0, p
                        do l = 0, n
                            do k = 0, m
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_fluids
                                    rhs_vf(E_idx)%sf(k, l, q) = rhs_vf(E_idx)%sf(k, l, q) &
                                        + djh_dx(k, l, q, i)
                                end do
                            end do
                        end do
                    end do
                    
                end if
            end if
        end if


    end subroutine s_compute_diffusion_rhs

    subroutine s_finalize_diffusion_module

        @:DEALLOCATE(Ds)
        @:DEALLOCATE(Ws)
        @:DEALLOCATE(alpha_K_dif, alpharho_K_dif, dF_KdP, F_K_dif, F_dif)
        @:DEALLOCATE(Gamma_dif, Pi_inf_dif, dYda, dYdP, dPdt, rho_dif)
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