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

    private; public s_get_diffusion, & 
s_compute_fd_gradient_diffusion, &
s_apply_scalar_divergence_theorem_diffusion, &
s_reconstruct_cell_boundary_values_diff, &
s_reconstruct_cell_boundary_values_diff_deriv, & 
s_initialize_diffusion_module, &
s_compute_diffusion_rhs, &
s_finalize_diffusion_module
    
    type(int_bounds_info) :: iv
    type(int_bounds_info) :: is1_diffusion, is2_diffusion, is3_diffusion
    !$acc declare create(is1_diffusion, is2_diffusion, is3_diffusion, iv)

    

    real(wp), allocatable, dimension(:, :) :: fd_coeff_x_d
    real(wp), allocatable, dimension(:, :) :: fd_coeff_y_d
    real(wp), allocatable, dimension(:, :) :: fd_coeff_z_d
    !$acc declare create(fd_coeff_x_d,fd_coeff_y_d,fd_coeff_z_d)

    real(wp), allocatable, dimension(:, :) :: Ds
    !$acc declare create(Ds)

    real(wp), allocatable, dimension(:) :: Ws
    !$acc declare create(Ws)

    real(wp), allocatable, dimension(:, :, :, :) :: dj_dx, dj_dy, dj_dz, djh_dx, djh_dy, djh_dz, dY_dx, dY_dy, dY_dz, alpha_K_dif, alpharho_K_dif, dF_KdP, F_K_dif
    !$acc declare create(dj_dx, dj_dy, dj_dz, djh_dx, djh_dy, djh_dz, dY_dx, dY_dy, dY_dz, alpha_K_dif, alpharho_K_dif, dF_KdP, F_K_dif)

    real(wp), allocatable, dimension(:, :, :) :: Gamma_dif, Pi_inf_dif, F_dif, dYda, dYdP, dPdt, rho_dif, dvel_dx, dvel_dy, dvel_dz
    !$acc declare create(Gamma_dif, Pi_inf_dif, F_dif, dYda, dYdP, dPdt, rho_dif, dvel_dx, dvel_dy, dvel_dz)

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

        !$acc loop seq
        do i = 1, num_fluids
            do j = 1, num_fluids
                Ds(i, j) = fluid_pp(i)%D(j)
            end do
        end do

        !$acc update device(Ds)
        !$acc enter data copyin(is1_diffusion, is2_diffusion, is3_diffusion, iv)


        @:ALLOCATE(Ws(1:num_fluids))

        !$acc loop seq
        do i = 1, num_fluids
            Ws(i) = fluid_pp(i)%W
        end do
        !$acc update device(Ws)

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
    

    !>  Computes diffusion terms
    !!  @param q_cons_vf Cell-averaged conservative variables
    !!  @param q_prim_vf Cell-averaged primitive variables
    !!  @param rhs_vf Cell-averaged RHS variables
    subroutine s_get_diffusion(qL_prim_rsx_vf, qL_prim_rsy_vf, qL_prim_rsz_vf, &
                             dqL_prim_dx_n, dqL_prim_dy_n, dqL_prim_dz_n, &
                             qL_prim, &
                             qR_prim_rsx_vf, qR_prim_rsy_vf, qR_prim_rsz_vf, &
                             dqR_prim_dx_n, dqR_prim_dy_n, dqR_prim_dz_n, &
                             qR_prim, &
                             q_prim_qp, &
                             dq_prim_dx_qp, dq_prim_dy_qp, dq_prim_dz_qp, &
                             ix, iy, iz)

        real(kind(0d0)), dimension(startx:, starty:, startz:, 1:), &
            intent(inout) :: qL_prim_rsx_vf, qR_prim_rsx_vf, &
                             qL_prim_rsy_vf, qR_prim_rsy_vf, &
                             qL_prim_rsz_vf, qR_prim_rsz_vf

        type(vector_field), dimension(num_dims), intent(inout) :: qL_prim, qR_prim

        type(vector_field), intent(in) :: q_prim_qp

        type(vector_field), dimension(1:num_dims), &
            intent(inout) :: dqL_prim_dx_n, dqR_prim_dx_n, &
                             dqL_prim_dy_n, dqR_prim_dy_n, &
                             dqL_prim_dz_n, dqR_prim_dz_n

        type(vector_field), dimension(1), intent(inout) :: dq_prim_dx_qp, dq_prim_dy_qp, dq_prim_dz_qp
        type(int_bounds_info), intent(in) :: ix, iy, iz

        integer :: i, j, k, l

        do i = 1, num_dims

            iv%beg = cont_idx%beg; iv%end = cont_idx%end

            !$acc update device(iv)

            call s_reconstruct_cell_boundary_values_diff( &
                q_prim_qp%vf(iv%beg:iv%end), &
                qL_prim_rsx_vf, qL_prim_rsy_vf, qL_prim_rsz_vf, &
                qR_prim_rsx_vf, qR_prim_rsy_vf, qR_prim_rsz_vf, &
                i, qL_prim(i)%vf(iv%beg:iv%end), qR_prim(i)%vf(iv%beg:iv%end), &
                ix, iy, iz)
        end do

        if (weno_Dif_flux) then
            ! Compute velocity gradient at cell centers using scalar
            ! divergence theorem
            do i = 1, num_dims
                if (i == 1) then
                    call s_apply_scalar_divergence_theorem_diffusion( &
                        qL_prim(i)%vf(iv%beg:iv%end), &
                        qR_prim(i)%vf(iv%beg:iv%end), &
                        dq_prim_dx_qp(1)%vf(iv%beg:iv%end), i, &
                        ix, iy, iz, iv, dx, m, buff_size)
                elseif (i == 2) then
                    call s_apply_scalar_divergence_theorem_diffusion( &
                        qL_prim(i)%vf(iv%beg:iv%end), &
                        qR_prim(i)%vf(iv%beg:iv%end), &
                        dq_prim_dy_qp(1)%vf(iv%beg:iv%end), i, &
                        ix, iy, iz, iv, dy, n, buff_size)
                else
                    call s_apply_scalar_divergence_theorem_diffusion( &
                        qL_prim(i)%vf(iv%beg:iv%end), &
                        qR_prim(i)%vf(iv%beg:iv%end), &
                        dq_prim_dz_qp(1)%vf(iv%beg:iv%end), i, &
                        ix, iy, iz, iv, dz, p, buff_size)
                end if
            end do

        else ! Compute velocity gradient at cell centers using finite differences

            iv%beg = cont_idx%beg; iv%end = cont_idx%end
            !$acc update device(iv)

            is1_diffusion = ix; is2_diffusion = iy; is3_diffusion = iz

            !$acc update device(is1_diffusion, is2_diffusion, is3_diffusion)

            !$acc parallel loop collapse(3) gang vector default(present)
            do l = is3_diffusion%beg, is3_diffusion%end
                do k = is2_diffusion%beg, is2_diffusion%end
                    do j = is1_diffusion%beg + 1, is1_diffusion%end
                        !$acc loop seq
                        do i = iv%beg, iv%end
                            dqL_prim_dx_n(1)%vf(i)%sf(j, k, l) = &
                                (q_prim_qp%vf(i)%sf(j, k, l) - &
                                 q_prim_qp%vf(i)%sf(j - 1, k, l))/ &
                                (x_cc(j) - x_cc(j - 1))
                        end do
                    end do
                end do
            end do

            !$acc parallel loop collapse(3) gang vector default(present)
            do l = is3_diffusion%beg, is3_diffusion%end
                do k = is2_diffusion%beg, is2_diffusion%end
                    do j = is1_diffusion%beg, is1_diffusion%end - 1
                        !$acc loop seq
                        do i = iv%beg, iv%end
                            dqR_prim_dx_n(1)%vf(i)%sf(j, k, l) = &
                                (q_prim_qp%vf(i)%sf(j + 1, k, l) - &
                                 q_prim_qp%vf(i)%sf(j, k, l))/ &
                                (x_cc(j + 1) - x_cc(j))
                        end do
                    end do
                end do
            end do

            if (n > 0) then

                !$acc parallel loop collapse(3) gang vector default(present)
                do l = is3_diffusion%beg, is3_diffusion%end
                    do j = is2_diffusion%beg + 1, is2_diffusion%end
                        do k = is1_diffusion%beg, is1_diffusion%end
                            !$acc loop seq
                            do i = iv%beg, iv%end
                                dqL_prim_dy_n(2)%vf(i)%sf(k, j, l) = &
                                    (q_prim_qp%vf(i)%sf(k, j, l) - &
                                     q_prim_qp%vf(i)%sf(k, j - 1, l))/ &
                                    (y_cc(j) - y_cc(j - 1))
                            end do
                        end do
                    end do
                end do

                !$acc parallel loop collapse(3) gang vector default(present)
                do l = is3_diffusion%beg, is3_diffusion%end
                    do j = is2_diffusion%beg, is2_diffusion%end - 1
                        do k = is1_diffusion%beg, is1_diffusion%end
                            !$acc loop seq
                            do i = iv%beg, iv%end
                                dqR_prim_dy_n(2)%vf(i)%sf(k, j, l) = &
                                    (q_prim_qp%vf(i)%sf(k, j + 1, l) - &
                                     q_prim_qp%vf(i)%sf(k, j, l))/ &
                                    (y_cc(j + 1) - y_cc(j))
                            end do
                        end do
                    end do
                end do

                !$acc parallel loop collapse(3) gang vector default(present)
                do l = is3_diffusion%beg, is3_diffusion%end
                    do j = is2_diffusion%beg + 1, is2_diffusion%end
                        do k = is1_diffusion%beg + 1, is1_diffusion%end - 1
                            !$acc loop seq
                            do i = iv%beg, iv%end
                                dqL_prim_dx_n(2)%vf(i)%sf(k, j, l) = &
                                    (dqL_prim_dx_n(1)%vf(i)%sf(k, j, l) + &
                                     dqR_prim_dx_n(1)%vf(i)%sf(k, j, l) + &
                                     dqL_prim_dx_n(1)%vf(i)%sf(k, j - 1, l) + &
                                     dqR_prim_dx_n(1)%vf(i)%sf(k, j - 1, l))

                                dqL_prim_dx_n(2)%vf(i)%sf(k, j, l) = 25d-2* &
                                                                     dqL_prim_dx_n(2)%vf(i)%sf(k, j, l)
                            end do
                        end do
                    end do
                end do

                !$acc parallel loop collapse(3) gang vector default(present)
                do l = is3_diffusion%beg, is3_diffusion%end
                    do j = is2_diffusion%beg, is2_diffusion%end - 1
                        do k = is1_diffusion%beg + 1, is1_diffusion%end - 1
                            !$acc loop seq
                            do i = iv%beg, iv%end
                                dqR_prim_dx_n(2)%vf(i)%sf(k, j, l) = &
                                    (dqL_prim_dx_n(1)%vf(i)%sf(k, j + 1, l) + &
                                     dqR_prim_dx_n(1)%vf(i)%sf(k, j + 1, l) + &
                                     dqL_prim_dx_n(1)%vf(i)%sf(k, j, l) + &
                                     dqR_prim_dx_n(1)%vf(i)%sf(k, j, l))

                                dqR_prim_dx_n(2)%vf(i)%sf(k, j, l) = 25d-2* &
                                                                     dqR_prim_dx_n(2)%vf(i)%sf(k, j, l)

                            end do
                        end do
                    end do
                end do

                !$acc parallel loop collapse(3) gang vector default(present)
                do l = is3_diffusion%beg, is3_diffusion%end
                    do k = is2_diffusion%beg + 1, is2_diffusion%end - 1
                        do j = is1_diffusion%beg + 1, is1_diffusion%end
                            !$acc loop seq
                            do i = iv%beg, iv%end
                                dqL_prim_dy_n(1)%vf(i)%sf(j, k, l) = &
                                    (dqL_prim_dy_n(2)%vf(i)%sf(j, k, l) + &
                                     dqR_prim_dy_n(2)%vf(i)%sf(j, k, l) + &
                                     dqL_prim_dy_n(2)%vf(i)%sf(j - 1, k, l) + &
                                     dqR_prim_dy_n(2)%vf(i)%sf(j - 1, k, l))

                                dqL_prim_dy_n(1)%vf(i)%sf(j, k, l) = 25d-2* &
                                                                     dqL_prim_dy_n(1)%vf(i)%sf(j, k, l)

                            end do
                        end do
                    end do
                end do

                !$acc parallel loop collapse(3) gang vector default(present)
                do l = is3_diffusion%beg, is3_diffusion%end
                    do k = is2_diffusion%beg + 1, is2_diffusion%end - 1
                        do j = is1_diffusion%beg, is1_diffusion%end - 1
                            !$acc loop seq
                            do i = iv%beg, iv%end
                                dqR_prim_dy_n(1)%vf(i)%sf(j, k, l) = &
                                    (dqL_prim_dy_n(2)%vf(i)%sf(j + 1, k, l) + &
                                     dqR_prim_dy_n(2)%vf(i)%sf(j + 1, k, l) + &
                                     dqL_prim_dy_n(2)%vf(i)%sf(j, k, l) + &
                                     dqR_prim_dy_n(2)%vf(i)%sf(j, k, l))

                                dqR_prim_dy_n(1)%vf(i)%sf(j, k, l) = 25d-2* &
                                                                     dqR_prim_dy_n(1)%vf(i)%sf(j, k, l)

                            end do
                        end do
                    end do
                end do

                if (p > 0) then

                    !$acc parallel loop collapse(3) gang vector default(present)
                    do j = is3_diffusion%beg + 1, is3_diffusion%end
                        do l = is2_diffusion%beg, is2_diffusion%end
                            do k = is1_diffusion%beg, is1_diffusion%end
                                !$acc loop seq
                                do i = iv%beg, iv%end

                                    dqL_prim_dz_n(3)%vf(i)%sf(k, l, j) = &
                                        (q_prim_qp%vf(i)%sf(k, l, j) - &
                                         q_prim_qp%vf(i)%sf(k, l, j - 1))/ &
                                        (z_cc(j) - z_cc(j - 1))
                                end do
                            end do
                        end do
                    end do

                    !$acc parallel loop collapse(3) gang vector default(present)
                    do j = is3_diffusion%beg, is3_diffusion%end - 1
                        do l = is2_diffusion%beg, is2_diffusion%end
                            do k = is1_diffusion%beg, is1_diffusion%end
                                !$acc loop seq
                                do i = iv%beg, iv%end

                                    dqR_prim_dz_n(3)%vf(i)%sf(k, l, j) = &
                                        (q_prim_qp%vf(i)%sf(k, l, j + 1) - &
                                         q_prim_qp%vf(i)%sf(k, l, j))/ &
                                        (z_cc(j + 1) - z_cc(j))
                                end do
                            end do
                        end do
                    end do

                    !$acc parallel loop collapse(3) gang vector default(present)
                    do l = is3_diffusion%beg + 1, is3_diffusion%end - 1
                        do k = is2_diffusion%beg, is2_diffusion%end
                            do j = is1_diffusion%beg + 1, is1_diffusion%end
                                !$acc loop seq
                                do i = iv%beg, iv%end

                                    dqL_prim_dz_n(1)%vf(i)%sf(j, k, l) = &
                                        (dqL_prim_dz_n(3)%vf(i)%sf(j, k, l) + &
                                         dqR_prim_dz_n(3)%vf(i)%sf(j, k, l) + &
                                         dqL_prim_dz_n(3)%vf(i)%sf(j - 1, k, l) + &
                                         dqR_prim_dz_n(3)%vf(i)%sf(j - 1, k, l))

                                    dqL_prim_dz_n(1)%vf(i)%sf(j, k, l) = 25d-2* &
                                                                         dqL_prim_dz_n(1)%vf(i)%sf(j, k, l)

                                end do
                            end do
                        end do
                    end do

                    !$acc parallel loop collapse(3) gang vector default(present)
                    do l = is3_diffusion%beg + 1, is3_diffusion%end - 1
                        do k = is2_diffusion%beg, is2_diffusion%end
                            do j = is1_diffusion%beg, is1_diffusion%end - 1
                                !$acc loop seq
                                do i = iv%beg, iv%end

                                    dqR_prim_dz_n(1)%vf(i)%sf(j, k, l) = &
                                        (dqL_prim_dz_n(3)%vf(i)%sf(j + 1, k, l) + &
                                         dqR_prim_dz_n(3)%vf(i)%sf(j + 1, k, l) + &
                                         dqL_prim_dz_n(3)%vf(i)%sf(j, k, l) + &
                                         dqR_prim_dz_n(3)%vf(i)%sf(j, k, l))

                                    dqR_prim_dz_n(1)%vf(i)%sf(j, k, l) = 25d-2* &
                                                                         dqR_prim_dz_n(1)%vf(i)%sf(j, k, l)

                                end do
                            end do
                        end do
                    end do

                    !$acc parallel loop collapse(3) gang vector default(present)
                    do l = is3_diffusion%beg + 1, is3_diffusion%end - 1
                        do j = is2_diffusion%beg + 1, is2_diffusion%end
                            do k = is1_diffusion%beg, is1_diffusion%end
                                !$acc loop seq
                                do i = iv%beg, iv%end

                                    dqL_prim_dz_n(2)%vf(i)%sf(k, j, l) = &
                                        (dqL_prim_dz_n(3)%vf(i)%sf(k, j, l) + &
                                         dqR_prim_dz_n(3)%vf(i)%sf(k, j, l) + &
                                         dqL_prim_dz_n(3)%vf(i)%sf(k, j - 1, l) + &
                                         dqR_prim_dz_n(3)%vf(i)%sf(k, j - 1, l))

                                    dqL_prim_dz_n(2)%vf(i)%sf(k, j, l) = 25d-2* &
                                                                         dqL_prim_dz_n(2)%vf(i)%sf(k, j, l)

                                end do
                            end do
                        end do
                    end do

                    !$acc parallel loop collapse(3) gang vector default(present)
                    do l = is3_diffusion%beg + 1, is3_diffusion%end - 1
                        do j = is2_diffusion%beg, is2_diffusion%end - 1
                            do k = is1_diffusion%beg, is1_diffusion%end
                                !$acc loop seq
                                do i = iv%beg, iv%end

                                    dqR_prim_dz_n(2)%vf(i)%sf(k, j, l) = &
                                        (dqL_prim_dz_n(3)%vf(i)%sf(k, j + 1, l) + &
                                         dqR_prim_dz_n(3)%vf(i)%sf(k, j + 1, l) + &
                                         dqL_prim_dz_n(3)%vf(i)%sf(k, j, l) + &
                                         dqR_prim_dz_n(3)%vf(i)%sf(k, j, l))

                                    dqR_prim_dz_n(2)%vf(i)%sf(k, j, l) = 25d-2* &
                                                                         dqR_prim_dz_n(2)%vf(i)%sf(k, j, l)

                                end do
                            end do
                        end do
                    end do

                    !$acc parallel loop collapse(3) gang vector default(present)
                    do j = is3_diffusion%beg + 1, is3_diffusion%end
                        do l = is2_diffusion%beg + 1, is2_diffusion%end - 1
                            do k = is1_diffusion%beg, is1_diffusion%end
                                !$acc loop seq
                                do i = iv%beg, iv%end

                                    dqL_prim_dy_n(3)%vf(i)%sf(k, l, j) = &
                                        (dqL_prim_dy_n(2)%vf(i)%sf(k, l, j) + &
                                         dqR_prim_dy_n(2)%vf(i)%sf(k, l, j) + &
                                         dqL_prim_dy_n(2)%vf(i)%sf(k, l, j - 1) + &
                                         dqR_prim_dy_n(2)%vf(i)%sf(k, l, j - 1))

                                    dqL_prim_dy_n(3)%vf(i)%sf(k, l, j) = 25d-2* &
                                                                         dqL_prim_dy_n(3)%vf(i)%sf(k, l, j)

                                end do
                            end do
                        end do
                    end do

                    !$acc parallel loop collapse(3) gang vector default(present)
                    do j = is3_diffusion%beg, is3_diffusion%end - 1
                        do l = is2_diffusion%beg + 1, is2_diffusion%end - 1
                            do k = is1_diffusion%beg, is1_diffusion%end
                                !$acc loop seq
                                do i = iv%beg, iv%end

                                    dqR_prim_dy_n(3)%vf(i)%sf(k, l, j) = &
                                        (dqL_prim_dy_n(2)%vf(i)%sf(k, l, j + 1) + &
                                         dqR_prim_dy_n(2)%vf(i)%sf(k, l, j + 1) + &
                                         dqL_prim_dy_n(2)%vf(i)%sf(k, l, j) + &
                                         dqR_prim_dy_n(2)%vf(i)%sf(k, l, j))

                                    dqR_prim_dy_n(3)%vf(i)%sf(k, l, j) = 25d-2* &
                                                                         dqR_prim_dy_n(3)%vf(i)%sf(k, l, j)

                                end do
                            end do
                        end do
                    end do
                    !$acc parallel loop collapse(3) gang vector default(present)
                    do j = is3_diffusion%beg + 1, is3_diffusion%end
                        do l = is2_diffusion%beg, is2_diffusion%end
                            do k = is1_diffusion%beg + 1, is1_diffusion%end - 1
                                !$acc loop seq
                                do i = iv%beg, iv%end

                                    dqL_prim_dx_n(3)%vf(i)%sf(k, l, j) = &
                                        (dqL_prim_dx_n(1)%vf(i)%sf(k, l, j) + &
                                         dqR_prim_dx_n(1)%vf(i)%sf(k, l, j) + &
                                         dqL_prim_dx_n(1)%vf(i)%sf(k, l, j - 1) + &
                                         dqR_prim_dx_n(1)%vf(i)%sf(k, l, j - 1))

                                    dqL_prim_dx_n(3)%vf(i)%sf(k, l, j) = 25d-2* &
                                                                         dqL_prim_dx_n(3)%vf(i)%sf(k, l, j)

                                end do
                            end do
                        end do
                    end do
                    !$acc parallel loop collapse(3) gang vector default(present)
                    do j = is3_diffusion%beg, is3_diffusion%end - 1
                        do l = is2_diffusion%beg, is2_diffusion%end
                            do k = is1_diffusion%beg + 1, is1_diffusion%end - 1
                                !$acc loop seq
                                do i = iv%beg, iv%end
                                    dqR_prim_dx_n(3)%vf(i)%sf(k, l, j) = &
                                        (dqL_prim_dx_n(1)%vf(i)%sf(k, l, j + 1) + &
                                         dqR_prim_dx_n(1)%vf(i)%sf(k, l, j + 1) + &
                                         dqL_prim_dx_n(1)%vf(i)%sf(k, l, j) + &
                                         dqR_prim_dx_n(1)%vf(i)%sf(k, l, j))

                                    dqR_prim_dx_n(3)%vf(i)%sf(k, l, j) = 25d-2* &
                                                                         dqR_prim_dx_n(3)%vf(i)%sf(k, l, j)

                                end do
                            end do
                        end do
                    end do

                    do i = iv%beg, iv%end
                        call s_compute_fd_gradient_diffusion(q_prim_qp%vf(i), &
                                                   dq_prim_dx_qp(1)%vf(i), &
                                                   dq_prim_dy_qp(1)%vf(i), &
                                                   dq_prim_dz_qp(1)%vf(i))
                    end do

                else

                    do i = iv%beg, iv%end
                        call s_compute_fd_gradient_diffusion(q_prim_qp%vf(i), &
                                                   dq_prim_dx_qp(1)%vf(i), &
                                                   dq_prim_dy_qp(1)%vf(i), &
                                                   dq_prim_dy_qp(1)%vf(i))
                    end do

                end if

            else

                do i = iv%beg, iv%end
                    call s_compute_fd_gradient_diffusion(q_prim_qp%vf(i), &
                                               dq_prim_dx_qp(1)%vf(i), &
                                               dq_prim_dx_qp(1)%vf(i), &
                                               dq_prim_dx_qp(1)%vf(i))
                end do

            end if

        end if

    end subroutine s_get_diffusion

    !>  Computes the scalar gradient fields via finite differences
        !!  @param var Variable to compute derivative of
        !!  @param grad_x First coordinate direction component of the derivative
        !!  @param grad_y Second coordinate direction component of the derivative
        !!  @param grad_z Third coordinate direction component of the derivative
        !!  @param norm Norm of the gradient vector
    subroutine s_compute_fd_gradient_diffusion(var, grad_x, grad_y, grad_z)

        type(scalar_field), intent(in) :: var
        type(scalar_field), intent(inout) :: grad_x
        type(scalar_field), intent(inout) :: grad_y
        type(scalar_field), intent(inout) :: grad_z
        type(int_bounds_info) :: ix, iy, iz

        integer :: j, k, l !< Generic loop iterators

        ix%beg = 1 - buff_size; ix%end = m + buff_size - 1
        if (n > 0) then
            iy%beg = 1 - buff_size; iy%end = n + buff_size - 1
        else
            iy%beg = 0; iy%end = 0
        end if

        if (p > 0) then
            iz%beg = 1 - buff_size; iz%end = p + buff_size - 1
        else
            iz%beg = 0; iz%end = 0
        end if

        is1_diffusion = ix; is2_diffusion = iy; is3_diffusion = iz

        !$acc update device(is1_diffusion, is2_diffusion, is3_diffusion)

        !$acc parallel loop collapse(3) gang vector default(present)
        do l = is3_diffusion%beg, is3_diffusion%end
            do k = is2_diffusion%beg, is2_diffusion%end
                do j = is1_diffusion%beg, is1_diffusion%end
                    grad_x%sf(j, k, l) = &
                        (var%sf(j + 1, k, l) - var%sf(j - 1, k, l))/ &
                        (x_cc(j + 1) - x_cc(j - 1))
                end do
            end do
        end do

        if (n > 0) then
            !$acc parallel loop collapse(3) gang vector
            do l = is3_diffusion%beg, is3_diffusion%end
                do k = is2_diffusion%beg, is2_diffusion%end
                    do j = is1_diffusion%beg, is1_diffusion%end
                        grad_y%sf(j, k, l) = &
                            (var%sf(j, k + 1, l) - var%sf(j, k - 1, l))/ &
                            (y_cc(k + 1) - y_cc(k - 1))
                    end do
                end do
            end do
        end if

        if (p > 0) then
            !$acc parallel loop collapse(3) gang vector
            do l = is3_diffusion%beg, is3_diffusion%end
                do k = is2_diffusion%beg, is2_diffusion%end
                    do j = is1_diffusion%beg, is1_diffusion%end
                        grad_z%sf(j, k, l) = &
                            (var%sf(j, k, l + 1) - var%sf(j, k, l - 1))/ &
                            (z_cc(l + 1) - z_cc(l - 1))
                    end do
                end do
            end do
        end if

        !$acc parallel loop collapse(2) gang vector default(present)
        do l = idwbuff(3)%beg, idwbuff(3)%end
            do k = idwbuff(2)%beg, idwbuff(2)%end
                grad_x%sf(idwbuff(1)%beg, k, l) = &
                    (-3d0*var%sf(idwbuff(1)%beg, k, l) + 4d0*var%sf(idwbuff(1)%beg + 1, k, l) - var%sf(idwbuff(1)%beg + 2, k, l))/ &
                    (x_cc(idwbuff(1)%beg + 2) - x_cc(idwbuff(1)%beg))
                grad_x%sf(idwbuff(1)%end, k, l) = &
                    (+3d0*var%sf(idwbuff(1)%end, k, l) - 4d0*var%sf(idwbuff(1)%end - 1, k, l) + var%sf(idwbuff(1)%end - 2, k, l))/ &
                    (x_cc(idwbuff(1)%end) - x_cc(idwbuff(1)%end - 2))
            end do
        end do
        if (n > 0) then
            !$acc parallel loop collapse(2) gang vector default(present)
            do l = idwbuff(3)%beg, idwbuff(3)%end
                do j = idwbuff(1)%beg, idwbuff(1)%end
                    grad_y%sf(j, idwbuff(2)%beg, l) = &
                        (-3d0*var%sf(j, idwbuff(2)%beg, l) + 4d0*var%sf(j, idwbuff(2)%beg + 1, l) - var%sf(j, idwbuff(2)%beg + 2, l))/ &
                        (y_cc(idwbuff(2)%beg + 2) - y_cc(idwbuff(2)%beg))
                    grad_y%sf(j, idwbuff(2)%end, l) = &
                        (+3d0*var%sf(j, idwbuff(2)%end, l) - 4d0*var%sf(j, idwbuff(2)%end - 1, l) + var%sf(j, idwbuff(2)%end - 2, l))/ &
                        (y_cc(idwbuff(2)%end) - y_cc(idwbuff(2)%end - 2))
                end do
            end do
            if (p > 0) then
                !$acc parallel loop collapse(2) gang vector default(present)
                do k = idwbuff(2)%beg, idwbuff(2)%end
                    do j = idwbuff(1)%beg, idwbuff(1)%end
                        grad_z%sf(j, k, idwbuff(3)%beg) = &
                            (-3d0*var%sf(j, k, idwbuff(3)%beg) + 4d0*var%sf(j, k, idwbuff(3)%beg + 1) - var%sf(j, k, idwbuff(3)%beg + 2))/ &
                            (z_cc(idwbuff(3)%beg + 2) - z_cc(idwbuff(3)%beg))
                        grad_z%sf(j, k, idwbuff(3)%end) = &
                            (+3d0*var%sf(j, k, idwbuff(3)%end) - 4d0*var%sf(j, k, idwbuff(3)%end - 1) + var%sf(j, k, idwbuff(3)%end - 2))/ &
                            (z_cc(idwbuff(3)%end) - z_cc(idwbuff(3)%end - 2))
                    end do
                end do
            end if
        end if

        if (bc_x%beg <= -3) then
            !$acc parallel loop collapse(2) gang vector default(present)
            do l = idwbuff(3)%beg, idwbuff(3)%end
                do k = idwbuff(2)%beg, idwbuff(2)%end
                    grad_x%sf(0, k, l) = (-3d0*var%sf(0, k, l) + 4d0*var%sf(1, k, l) - var%sf(2, k, l))/ &
                                         (x_cc(2) - x_cc(0))
                end do
            end do
        end if
        if (bc_x%end <= -3) then
            !$acc parallel loop collapse(2) gang vector default(present)
            do l = idwbuff(3)%beg, idwbuff(3)%end
                do k = idwbuff(2)%beg, idwbuff(2)%end
                    grad_x%sf(m, k, l) = (3d0*var%sf(m, k, l) - 4d0*var%sf(m - 1, k, l) + var%sf(m - 2, k, l))/ &
                                         (x_cc(m) - x_cc(m - 2))
                end do
            end do
        end if
        if (n > 0) then
            if (bc_y%beg <= -3 .and. bc_y%beg /= -13) then
                !$acc parallel loop collapse(2) gang vector default(present)
                do l = idwbuff(3)%beg, idwbuff(3)%end
                    do j = idwbuff(1)%beg, idwbuff(1)%end
                        grad_y%sf(j, 0, l) = (-3d0*var%sf(j, 0, l) + 4d0*var%sf(j, 1, l) - var%sf(j, 2, l))/ &
                                             (y_cc(2) - y_cc(0))
                    end do
                end do
            end if
            if (bc_y%end <= -3) then
                !$acc parallel loop collapse(2) gang vector default(present)
                do l = idwbuff(3)%beg, idwbuff(3)%end
                    do j = idwbuff(1)%beg, idwbuff(1)%end
                        grad_y%sf(j, n, l) = (3d0*var%sf(j, n, l) - 4d0*var%sf(j, n - 1, l) + var%sf(j, n - 2, l))/ &
                                             (y_cc(n) - y_cc(n - 2))
                    end do
                end do
            end if
            if (p > 0) then
                if (bc_z%beg <= -3) then
                    !$acc parallel loop collapse(2) gang vector default(present)
                    do k = idwbuff(2)%beg, idwbuff(2)%end
                        do j = idwbuff(1)%beg, idwbuff(1)%end
                            grad_z%sf(j, k, 0) = &
                                (-3d0*var%sf(j, k, 0) + 4d0*var%sf(j, k, 1) - var%sf(j, k, 2))/ &
                                (z_cc(2) - z_cc(0))
                        end do
                    end do
                end if
                if (bc_z%end <= -3) then
                    !$acc parallel loop collapse(2) gang vector default(present)
                    do k = idwbuff(2)%beg, idwbuff(2)%end
                        do j = idwbuff(1)%beg, idwbuff(1)%end
                            grad_z%sf(j, k, p) = &
                                (3d0*var%sf(j, k, p) - 4d0*var%sf(j, k, p - 1) + var%sf(j, k, p - 2))/ &
                                (z_cc(p) - z_cc(p - 2))
                        end do
                    end do
                end if
            end if
        end if

    end subroutine s_compute_fd_gradient_diffusion

    subroutine s_apply_scalar_divergence_theorem_diffusion(vL_vf, vR_vf, &
                                                 dv_ds_vf, &
                                                 norm_dir, &
                                                 ix, iy, iz, iv_in, &
                                                 dL, dim, buff_size_in)

        type(scalar_field), &
            dimension(iv%beg:iv%end), &
            intent(in) :: vL_vf, vR_vf

        type(scalar_field), &
            dimension(iv%beg:iv%end), &
            intent(inout) :: dv_ds_vf

        integer, intent(in) :: norm_dir
        type(int_bounds_info), intent(in) :: ix, iy, iz, iv_in
        integer, intent(in) :: dim, buff_size_in

        ! arrays of cell widths
        real(kind(0d0)), dimension(-buff_size_in:dim + buff_size_in), intent(in) :: dL

        integer :: i, j, k, l !< Generic loop iterators

        is1_diffusion = ix
        is2_diffusion = iy
        is3_diffusion = iz
        iv = iv_in

        !$acc update device(is1_diffusion, is2_diffusion, is3_diffusion, iv)

        ! First-Order Spatial Derivatives in x-direction ===================
        if (norm_dir == 1) then

            ! A general application of the scalar divergence theorem that
            ! utilizes the left and right cell-boundary integral-averages,
            ! inside each cell, or an arithmetic mean of these two at the
            ! cell-boundaries, to calculate the cell-averaged first-order
            ! spatial derivatives inside the cell.

            !$acc parallel loop collapse(3) gang vector default(present)
            do l = is3_diffusion%beg, is3_diffusion%end
                do k = is2_diffusion%beg, is2_diffusion%end
                    do j = is1_diffusion%beg + 1, is1_diffusion%end - 1
                        !$acc loop seq
                        do i = iv%beg, iv%end
                            dv_ds_vf(i)%sf(j, k, l) = &
                                1d0/((1d0 + wa_flg)*dL(j)) &
                                *(wa_flg*vL_vf(i)%sf(j + 1, k, l) &
                                  + vR_vf(i)%sf(j, k, l) &
                                  - vL_vf(i)%sf(j, k, l) &
                                  - wa_flg*vR_vf(i)%sf(j - 1, k, l))
                        end do
                    end do
                end do
            end do

            ! END: First-Order Spatial Derivatives in x-direction ==============

            ! First-Order Spatial Derivatives in y-direction ===================
        elseif (norm_dir == 2) then

            ! A general application of the scalar divergence theorem that
            ! utilizes the left and right cell-boundary integral-averages,
            ! inside each cell, or an arithmetic mean of these two at the
            ! cell-boundaries, to calculate the cell-averaged first-order
            ! spatial derivatives inside the cell.

            !$acc parallel loop collapse(3) gang vector default(present)
            do l = is3_diffusion%beg, is3_diffusion%end
                do k = is2_diffusion%beg + 1, is2_diffusion%end - 1
                    do j = is1_diffusion%beg, is1_diffusion%end
                        !$acc loop seq
                        do i = iv%beg, iv%end
                            dv_ds_vf(i)%sf(j, k, l) = &
                                1d0/((1d0 + wa_flg)*dL(k)) &
                                *(wa_flg*vL_vf(i)%sf(j, k + 1, l) &
                                  + vR_vf(i)%sf(j, k, l) &
                                  - vL_vf(i)%sf(j, k, l) &
                                  - wa_flg*vR_vf(i)%sf(j, k - 1, l))
                        end do
                    end do
                end do
            end do

            ! END: First-Order Spatial Derivatives in y-direction ==============

            ! First-Order Spatial Derivatives in z-direction ===================
        else

            ! A general application of the scalar divergence theorem that
            ! utilizes the left and right cell-boundary integral-averages,
            ! inside each cell, or an arithmetic mean of these two at the
            ! cell-boundaries, to calculate the cell-averaged first-order
            ! spatial derivatives inside the cell.

            !$acc parallel loop collapse(3) gang vector default(present)
            do l = is3_diffusion%beg + 1, is3_diffusion%end - 1
                do k = is2_diffusion%beg, is2_diffusion%end
                    do j = is1_diffusion%beg, is1_diffusion%end
                        !$acc loop seq
                        do i = iv%beg, iv%end
                            dv_ds_vf(i)%sf(j, k, l) = &
                                1d0/((1d0 + wa_flg)*dL(l)) &
                                *(wa_flg*vL_vf(i)%sf(j, k, l + 1) &
                                  + vR_vf(i)%sf(j, k, l) &
                                  - vL_vf(i)%sf(j, k, l) &
                                  - wa_flg*vR_vf(i)%sf(j, k, l - 1))
                        end do
                    end do
                end do
            end do

        end if
        ! END: First-Order Spatial Derivatives in z-direction ==============

    end subroutine s_apply_scalar_divergence_theorem_diffusion

 
    subroutine s_reconstruct_cell_boundary_values_diff(v_vf, vL_x, vL_y, vL_z, vR_x, vR_y, vR_z, &
                                                       norm_dir, vL_prim_vf, vR_prim_vf, ix, iy, iz)

        type(scalar_field), dimension(iv%beg:iv%end), intent(in) :: v_vf
        type(scalar_field), dimension(iv%beg:iv%end), intent(inout) :: vL_prim_vf, vR_prim_vf

        real(kind(0d0)), dimension(startx:, starty:, startz:, 1:), intent(inout) :: vL_x, vL_y, vL_z, vR_x, vR_y, vR_z
        integer, intent(in) :: norm_dir
        type(int_bounds_info), intent(in) :: ix, iy, iz

        integer :: weno_dir !< Coordinate direction of the WENO reconstruction

        integer :: i, j, k, l

        ! Reconstruction in s1-direction ===================================

        if (norm_dir == 1) then
            is1_diffusion = ix; is2_diffusion = iy; is3_diffusion = iz
            weno_dir = 1; is1_diffusion%beg = is1_diffusion%beg + weno_polyn
            is1_diffusion%end = is1_diffusion%end - weno_polyn

        elseif (norm_dir == 2) then
            is1_diffusion = iy; is2_diffusion = ix; is3_diffusion = iz
            weno_dir = 2; is1_diffusion%beg = is1_diffusion%beg + weno_polyn
            is1_diffusion%end = is1_diffusion%end - weno_polyn

        else
            is1_diffusion = iz; is2_diffusion = iy; is3_diffusion = ix
            weno_dir = 3; is1_diffusion%beg = is1_diffusion%beg + weno_polyn
            is1_diffusion%end = is1_diffusion%end - weno_polyn

        end if

        !$acc update device(is1_diffusion, is2_diffusion, is3_diffusion, iv)

        if (n > 0) then
            if (p > 0) then
                call s_weno(v_vf(iv%beg:iv%end), &
                            vL_x(:, :, :, iv%beg:iv%end), vL_y(:, :, :, iv%beg:iv%end), vL_z(:, :, :, iv%beg:iv%end), vR_x(:, :, :, iv%beg:iv%end), vR_y(:, :, :, iv%beg:iv%end), vR_z(:, :, :, iv%beg:iv%end), &
                            norm_dir, weno_dir, &
                            is1_diffusion, is2_diffusion, is3_diffusion)
            else
                call s_weno(v_vf(iv%beg:iv%end), &
                            vL_x(:, :, :, iv%beg:iv%end), vL_y(:, :, :, iv%beg:iv%end), vL_z(:, :, :, :), vR_x(:, :, :, iv%beg:iv%end), vR_y(:, :, :, iv%beg:iv%end), vR_z(:, :, :, :), &
                            norm_dir, weno_dir, &
                            is1_diffusion, is2_diffusion, is3_diffusion)
            end if
        else
            call s_weno(v_vf(iv%beg:iv%end), &
                        vL_x(:, :, :, iv%beg:iv%end), vL_y(:, :, :, :), vL_z(:, :, :, :), vR_x(:, :, :, iv%beg:iv%end), vR_y(:, :, :, :), vR_z(:, :, :, :), &
                        norm_dir, weno_dir, &
                        is1_diffusion, is2_diffusion, is3_diffusion)
        end if

        if (diffusion) then
            if (weno_Dif_flux) then
                if (norm_dir == 2) then
                    !$acc parallel loop collapse(4) gang vector default(present)
                    do i = iv%beg, iv%end
                        do l = is3_diffusion%beg, is3_diffusion%end
                            do j = is1_diffusion%beg, is1_diffusion%end
                                do k = is2_diffusion%beg, is2_diffusion%end
                                    vL_prim_vf(i)%sf(k, j, l) = vL_y(j, k, l, i)
                                    vR_prim_vf(i)%sf(k, j, l) = vR_y(j, k, l, i)
                                end do
                            end do
                        end do
                    end do
                elseif (norm_dir == 3) then
                    !$acc parallel loop collapse(4) gang vector default(present)
                    do i = iv%beg, iv%end
                        do j = is1_diffusion%beg, is1_diffusion%end
                            do k = is2_diffusion%beg, is2_diffusion%end
                                do l = is3_diffusion%beg, is3_diffusion%end
                                    vL_prim_vf(i)%sf(l, k, j) = vL_z(j, k, l, i)
                                    vR_prim_vf(i)%sf(l, k, j) = vR_z(j, k, l, i)
                                end do
                            end do
                        end do
                    end do
                elseif (norm_dir == 1) then
                    !$acc parallel loop collapse(4) gang vector default(present)
                    do i = iv%beg, iv%end
                        do l = is3_diffusion%beg, is3_diffusion%end
                            do k = is2_diffusion%beg, is2_diffusion%end
                                do j = is1_diffusion%beg, is1_diffusion%end
                                    vL_prim_vf(i)%sf(j, k, l) = vL_x(j, k, l, i)
                                    vR_prim_vf(i)%sf(j, k, l) = vR_x(j, k, l, i)
                                end do
                            end do
                        end do
                    end do
                end if
            end if
        end if

        ! ==================================================================

    end subroutine s_reconstruct_cell_boundary_values_diff

    subroutine s_reconstruct_cell_boundary_values_diff_deriv(v_vf, vL_x, vL_y, vL_z, vR_x, vR_y, vR_z, &
                                                             norm_dir, vL_prim_vf, vR_prim_vf, ix, iy, iz)

        type(scalar_field), dimension(iv%beg:iv%end), intent(in) :: v_vf
        real(kind(0d0)), dimension(startx:, starty:, startz:, iv%beg:), intent(inout) :: vL_x, vL_y, vL_z, vR_x, vR_y, vR_z
        type(scalar_field), dimension(iv%beg:iv%end), intent(inout) :: vL_prim_vf, vR_prim_vf
        type(int_bounds_info), intent(in) :: ix, iy, iz

        integer, intent(IN) :: norm_dir

        integer :: weno_dir !< Coordinate direction of the WENO reconstruction

        integer :: i, j, k, l
        ! Reconstruction in s1-direction ===================================

        if (norm_dir == 1) then
            is1_diffusion = ix; is2_diffusion = iy; is3_diffusion = iz
            weno_dir = 1; is1_diffusion%beg = is1_diffusion%beg + weno_polyn
            is1_diffusion%end = is1_diffusion%end - weno_polyn

        elseif (norm_dir == 2) then
            is1_diffusion = iy; is2_diffusion = ix; is3_diffusion = iz
            weno_dir = 2; is1_diffusion%beg = is1_diffusion%beg + weno_polyn
            is1_diffusion%end = is1_diffusion%end - weno_polyn

        else
            is1_diffusion = iz; is2_diffusion = iy; is3_diffusion = ix
            weno_dir = 3; is1_diffusion%beg = is1_diffusion%beg + weno_polyn
            is1_diffusion%end = is1_diffusion%end - weno_polyn

        end if

        !$acc update device(is1_diffusion, is2_diffusion, is3_diffusion, iv)

        if (n > 0) then
            if (p > 0) then

                call s_weno(v_vf(iv%beg:iv%end), &
                            vL_x(:, :, :, iv%beg:iv%end), vL_y(:, :, :, iv%beg:iv%end), vL_z(:, :, :, iv%beg:iv%end), vR_x(:, :, :, iv%beg:iv%end), vR_y(:, :, :, iv%beg:iv%end), vR_z(:, :, :, iv%beg:iv%end), &
                            norm_dir, weno_dir, &
                            is1_diffusion, is2_diffusion, is3_diffusion)
            else
                call s_weno(v_vf(iv%beg:iv%end), &
                            vL_x(:, :, :, iv%beg:iv%end), vL_y(:, :, :, iv%beg:iv%end), vL_z(:, :, :, :), vR_x(:, :, :, iv%beg:iv%end), vR_y(:, :, :, iv%beg:iv%end), vR_z(:, :, :, :), &
                            norm_dir, weno_dir, &
                            is1_diffusion, is2_diffusion, is3_diffusion)
            end if
        else

            call s_weno(v_vf(iv%beg:iv%end), &
                        vL_x(:, :, :, iv%beg:iv%end), vL_y(:, :, :, :), vL_z(:, :, :, :), vR_x(:, :, :, iv%beg:iv%end), vR_y(:, :, :, :), vR_z(:, :, :, :), &
                        norm_dir, weno_dir, &
                        is1_diffusion, is2_diffusion, is3_diffusion)
        end if

        if (diffusion) then
            if (weno_Dif_flux) then
                if (norm_dir == 2) then
                    !$acc parallel loop collapse(4) gang vector default(present)
                    do i = iv%beg, iv%end
                        do l = is3_diffusion%beg, is3_diffusion%end
                            do j = is1_diffusion%beg, is1_diffusion%end
                                do k = is2_diffusion%beg, is2_diffusion%end
                                    vL_prim_vf(i)%sf(k, j, l) = vL_y(j, k, l, i)
                                    vR_prim_vf(i)%sf(k, j, l) = vR_y(j, k, l, i)
                                end do
                            end do
                        end do
                    end do
                elseif (norm_dir == 3) then
                    !$acc parallel loop collapse(4) gang vector default(present)
                    do i = iv%beg, iv%end
                        do j = is1_diffusion%beg, is1_diffusion%end
                            do k = is2_diffusion%beg, is2_diffusion%end
                                do l = is3_diffusion%beg, is3_diffusion%end
                                    vL_prim_vf(i)%sf(l, k, j) = vL_z(j, k, l, i)
                                    vR_prim_vf(i)%sf(l, k, j) = vR_z(j, k, l, i)
                                end do
                            end do
                        end do
                    end do
                elseif (norm_dir == 1) then
                    !$acc parallel loop collapse(4) gang vector default(present)
                    do i = iv%beg, iv%end
                        do l = is3_diffusion%beg, is3_diffusion%end
                            do k = is2_diffusion%beg, is2_diffusion%end
                                do j = is1_diffusion%beg, is1_diffusion%end
                                    vL_prim_vf(i)%sf(j, k, l) = vL_x(j, k, l, i)
                                    vR_prim_vf(i)%sf(j, k, l) = vR_x(j, k, l, i)
                                end do
                            end do
                        end do
                    end do
                end if
            end if
        end if
        ! ==================================================================

    end subroutine s_reconstruct_cell_boundary_values_diff_deriv

    subroutine s_compute_diffusion_rhs(idir, j_prim_vf, q_prim_vf, rhs_vf)

        integer, intent(in) :: idir
        !type(scalar_field), dimension(sys_size), intent(in) :: dj_prim_dx_qp, dj_prim_dy_qp, dj_prim_dz_qp
        type(scalar_field), dimension(sys_size), intent(inout) :: j_prim_vf, q_prim_vf, rhs_vf

        integer :: i, k, l, q, r !< Loop variables
        real(wp), dimension(2) :: dif_flg
        real(wp) :: W1, W2, W3, D12, D13, D23, W_dif, W_fac, inv_denom
        dif_flg(1) = 1._wp; dif_flg(2) = -1._wp

        if (cyl_coord) then
            if (idir == 1) then
                !$acc parallel loop collapse(4) gang vector default(present)
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
                !$acc end parallel loop

                !$acc parallel loop collapse(4) gang vector default(present)
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
                !$acc end parallel loop
                
                !$acc parallel loop collapse(4) gang vector default(present)
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
                !$acc end parallel loop

                ! For now, only computes for BINARY diffusion. Computed dPdt without viscous term and velocity divergence
                !$acc parallel loop collapse(3) gang vector default(present)
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
                !$acc end parallel loop

                !$acc parallel loop collapse(4) gang vector default(present)
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
                !$acc end parallel loop

                !$acc parallel loop collapse(3) gang vector default(present)
                do q = 0, p
                    do l = 0, n
                        do k = 0, m
                            dvel_dx(k, l, q) = 0._wp               
                        end do
                    end do
                end do
                !$acc end parallel loop


                !$acc parallel loop collapse(4) gang vector default(present)
                do q = 0, p
                    do l = 0, n
                        do k = 0, m
                            do r = -fd_number, fd_number
                                dvel_dx(k, l, q) = dvel_dx(k, l, q) + q_prim_vf(momxb + idir - 1)%sf(k + r, l, q)*fd_coeff_x_d(r, k)               
                            end do               
                        end do
                    end do
                end do
                !$acc end parallel loop

                !$acc parallel loop collapse(5) gang vector default(present)
                do q = 0, p
                    do l = 0, n
                        do k = -fd_number, m + fd_number
                            do i = 1, num_fluids
                                do r = -fd_number, fd_number
                                    dY_dx(k, l, q, i) = dY_dx(k, l, q, i) &
                                        + j_prim_vf(i)%sf(k + r, l, q)*fd_coeff_x_d(r, k)
                                end do
                            end do
                        end do
                    end do
                end do
                !$acc end parallel loop

                !$acc parallel loop collapse(5) gang vector default(present)
                do q = 0, p
                    do l = 0, n
                        do k = 0, m
                            do i = 1, num_fluids
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
                !$acc end parallel loop

                !Valid for any number of species
                ! species continuity
                !$acc parallel loop collapse(4) gang vector default(present)
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
                !$acc end parallel loop

                !Only valid for binary diffusion
                !volume fraction
                !$acc parallel loop collapse(4) gang vector default(present)
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
                !$acc end parallel loop

                !Valid for any number of species
                !energy
                !$acc parallel loop collapse(4) gang vector default(present)
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
                !$acc end parallel loop
            elseif (idir == 2) then
                !$acc parallel loop collapse(4) gang vector default(present)
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
                !$acc end parallel loop
                
                !$acc parallel loop collapse(4) gang vector default(present)
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
                !$acc end parallel loop
                
                !$acc parallel loop collapse(4) gang vector default(present)
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
                !$acc end parallel loop

                ! For now, only computes for BINARY diffusion. Computed dPdt without viscous term and velocity divergence
                !$acc parallel loop collapse(3) gang vector default(present)
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
                !$acc end parallel loop

                !$acc parallel loop collapse(4) gang vector default(present)
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
                !$acc end parallel loop

                !$acc parallel loop collapse(3) gang vector default(present)
                do q = 0, p
                    do l = 0, n
                        do k = 0, m
                            dvel_dy(k, l, q) = 0._wp               
                        end do
                    end do
                end do
                !$acc end parallel loop


                !$acc parallel loop collapse(4) gang vector default(present)
                do q = 0, p
                    do l = 0, n
                        do k = 0, m
                            do r = -fd_number, fd_number
                                dvel_dy(k, l, q) = dvel_dy(k, l, q) + q_prim_vf(momxb + idir - 1)%sf(k, l + r, q)*fd_coeff_y_d(r, l) + &
                                    q_prim_vf(momyb + idir - 1)%sf(k, l, q) / y_cc(l)         
                            end do               
                        end do
                    end do
                end do
                !$acc end parallel loop

                !$acc parallel loop collapse(5) gang vector default(present)
                do q = 0, p
                    do l = -fd_number, n + fd_number
                        do k = 0, m
                            do i = 1, num_fluids
                                do r = -fd_number, fd_number
                                    dY_dy(k, l, q, i) = dY_dy(k, l, q, i) &
                                        + j_prim_vf(i)%sf(k, l + r, q)*fd_coeff_y_d(r, l)
                                end do
                            end do
                        end do
                    end do
                end do
                !$acc end parallel loop

                !$acc parallel loop collapse(5) gang vector default(present)
                do q = 0, p
                    do l = 0, n
                        do k = 0, m
                            do i = 1, num_fluids
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
                !$acc end parallel loop

                !Valid for any number of species
                ! species continuity
                !$acc parallel loop collapse(4) gang vector default(present)
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
                !$acc end parallel loop

                !Only valid for binary diffusion
                !volume fraction
                !$acc parallel loop collapse(4) gang vector default(present)
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
                !$acc end parallel loop

                !Valid for any number of species
                !energy
                !$acc parallel loop collapse(4) gang vector default(present)
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
                !$acc end parallel loop     

            elseif (idir == 3) then
                !$acc parallel loop collapse(4) gang vector default(present)
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
                !$acc end parallel loop
                
                !$acc parallel loop collapse(4) gang vector default(present)
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
                !$acc end parallel loop
                
                !$acc parallel loop collapse(4) gang vector default(present)
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
                !$acc end parallel loop

                ! For now, only computes for BINARY diffusion. Computed dPdt without viscous term and velocity divergence
                !$acc parallel loop collapse(3) gang vector default(present)
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
                !$acc end parallel loop

                !$acc parallel loop collapse(4) gang vector default(present)
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
                !$acc end parallel loop

                !$acc parallel loop collapse(3) gang vector default(present)
                do q = 0, p
                    do l = 0, n
                        do k = 0, m
                            dvel_dz(k, l, q) = 0._wp               
                        end do
                    end do
                end do
                !$acc end parallel loop


                !$acc parallel loop collapse(4) gang vector default(present)
                do q = 0, p
                    do l = 0, n
                        do k = 0, m
                            do r = -fd_number, fd_number
                                dvel_dz(k, l, q) = dvel_dz(k, l, q) + q_prim_vf(momxb + idir - 1)%sf(k, l, q + r)*fd_coeff_z_d(r, q) / y_cc(l)               
                            end do               
                        end do
                    end do
                end do
                !$acc end parallel loop

                !$acc parallel loop collapse(5) gang vector default(present)
                do q = -fd_number, p + fd_number
                    do l = 0, n
                        do k = 0, m
                            do i = 1, num_fluids
                                do r = -fd_number, fd_number
                                    dY_dz(k, l, q, i) = dY_dz(k, l, q, i) &
                                        + j_prim_vf(i)%sf(k, l, q + r)*fd_coeff_z_d(r, q) / y_cc(l)
                                end do
                            end do
                        end do
                    end do
                end do
                !$acc end parallel loop
    
                !$acc parallel loop collapse(5) gang vector default(present)
                do q = 0, p
                    do l = 0, n
                        do k = 0, m
                            do i = 1, num_fluids
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
                !$acc end parallel loop

                !Valid for any number of species
                ! species continuity
                !$acc parallel loop collapse(4) gang vector default(present)
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
                !$acc end parallel loop

                !Only valid for binary diffusion
                !volume fraction
                !$acc parallel loop collapse(4) gang vector default(present)
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
                !$acc end parallel loop

                !Valid for any number of species
                !energy
                !$acc parallel loop collapse(4) gang vector default(present)
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
                !$acc end parallel loop

            end if

        else !cartesian coordinates
            if (num_fluids == 2) then
                if (idir == 1) then
                    !$acc parallel loop collapse(4) gang vector default(present)
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
                    !$acc end parallel loop
                    
                    !$acc parallel loop collapse(4) gang vector default(present)
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
                    !$acc end parallel loop
                    
                    !$acc parallel loop collapse(4) gang vector default(present)
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
                    !$acc end parallel loop

                    ! For now, only computes for BINARY diffusion. Computed dPdt without viscous term and velocity divergence
                    !$acc parallel loop collapse(3) gang vector default(present)
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
                    !$acc end parallel loop

                    !$acc parallel loop collapse(4) gang vector default(present)
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
                    !$acc end parallel loop

                    !$acc parallel loop collapse(3) gang vector default(present)
                    do q = 0, p
                        do l = 0, n
                            do k = 0, m
                                dvel_dx(k, l, q) = 0._wp               
                            end do
                        end do
                    end do
                    !$acc end parallel loop


                    !$acc parallel loop collapse(4) gang vector default(present)
                    do q = 0, p
                        do l = 0, n
                            do k = 0, m
                                do r = -fd_number, fd_number
                                    dvel_dx(k, l, q) = dvel_dx(k, l, q) + q_prim_vf(momxb + idir - 1)%sf(k + r, l, q)*fd_coeff_x_d(r, k)               
                                end do               
                            end do
                        end do
                    end do
                    !$acc end parallel loop

                    !$acc parallel loop collapse(5) gang vector default(present)
                    do q = 0, p
                        do l = 0, n
                            do k = -fd_number, m + fd_number
                                do i = 1, num_fluids
                                    do r = -fd_number, fd_number
                                        dY_dx(k, l, q, i) = dY_dx(k, l, q, i) &
                                            + j_prim_vf(i)%sf(k + r, l, q)*fd_coeff_x_d(r, k)
                                    end do
                                end do
                            end do
                        end do
                    end do
                    !$acc end parallel loop
        
                    !$acc parallel loop collapse(5) gang vector default(present)
                    do q = 0, p
                        do l = 0, n
                            do k = 0, m
                                do i = 1, num_fluids
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
                    !$acc end parallel loop

                    !Valid for any number of species
                    ! species continuity
                    !$acc parallel loop collapse(4) gang vector default(present)
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
                    !$acc end parallel loop

                    !Only valid for binary diffusion
                    !volume fraction
                    !$acc parallel loop collapse(4) gang vector default(present)
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
                    !$acc end parallel loop

                    !Valid for any number of species
                    !energy
                    !$acc parallel loop collapse(4) gang vector default(present)
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
                    !$acc end parallel loop
                elseif (idir == 2) then
                    !$acc parallel loop collapse(4) gang vector default(present)
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
                    !$acc end parallel loop
                    
                    !$acc parallel loop collapse(4) gang vector default(present)
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
                    !$acc end parallel loop
                    
                    !$acc parallel loop collapse(4) gang vector default(present)
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
                    !$acc end parallel loop

                    ! For now, only computes for BINARY diffusion. Computed dPdt without viscous term and velocity divergence
                    !$acc parallel loop collapse(3) gang vector default(present)
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
                    !$acc end parallel loop

                    !$acc parallel loop collapse(4) gang vector default(present)
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
                    !$acc end parallel loop

                    !$acc parallel loop collapse(3) gang vector default(present)
                    do q = 0, p
                        do l = 0, n
                            do k = 0, m
                                dvel_dy(k, l, q) = 0._wp               
                            end do
                        end do
                    end do
                    !$acc end parallel loop


                    !$acc parallel loop collapse(4) gang vector default(present)
                    do q = 0, p
                        do l = 0, n
                            do k = 0, m
                                do r = -fd_number, fd_number
                                    dvel_dy(k, l, q) = dvel_dy(k, l, q) + q_prim_vf(momxb + idir - 1)%sf(k, l + r, q)*fd_coeff_y_d(r, l)               
                                end do               
                            end do
                        end do
                    end do
                    !$acc end parallel loop

                    !$acc parallel loop collapse(5) gang vector default(present)
                    do q = 0, p
                        do l = -fd_number, n + fd_number
                            do k = 0, m
                                do i = 1, num_fluids
                                    do r = -fd_number, fd_number
                                        dY_dy(k, l, q, i) = dY_dy(k, l, q, i) &
                                            + j_prim_vf(i)%sf(k, l + r, q)*fd_coeff_y_d(r, l)
                                    end do
                                end do
                            end do
                        end do
                    end do
                    !$acc end parallel loop
        
                    !$acc parallel loop collapse(5) gang vector default(present)
                    do q = 0, p
                        do l = 0, n
                            do k = 0, m
                                do i = 1, num_fluids
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
                    !$acc end parallel loop

                    !Valid for any number of species
                    ! species continuity
                    !$acc parallel loop collapse(4) gang vector default(present)
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
                    !$acc end parallel loop

                    !Only valid for binary diffusion
                    !volume fraction
                    !$acc parallel loop collapse(4) gang vector default(present)
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
                    !$acc end parallel loop

                    !Valid for any number of species
                    !energy
                    !$acc parallel loop collapse(4) gang vector default(present)
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
                    !$acc end parallel loop           
                elseif (idir == 3) then
                    !$acc parallel loop collapse(4) gang vector default(present)
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
                    !$acc end parallel loop
                    
                    !$acc parallel loop collapse(4) gang vector default(present)
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
                    !$acc end parallel loop
                    
                    !$acc parallel loop collapse(4) gang vector default(present)
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
                    !$acc end parallel loop

                    ! For now, only computes for BINARY diffusion. Computed dPdt without viscous term and velocity divergence
                    !$acc parallel loop collapse(3) gang vector default(present)
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
                    !$acc end parallel loop

                    !$acc parallel loop collapse(4) gang vector default(present)
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
                    !$acc end parallel loop

                    !$acc parallel loop collapse(3) gang vector default(present)
                    do q = 0, p
                        do l = 0, n
                            do k = 0, m
                                dvel_dz(k, l, q) = 0._wp               
                            end do
                        end do
                    end do
                    !$acc end parallel loop


                    !$acc parallel loop collapse(4) gang vector default(present)
                    do q = 0, p
                        do l = 0, n
                            do k = 0, m
                                do r = -fd_number, fd_number
                                    dvel_dz(k, l, q) = dvel_dz(k, l, q) + q_prim_vf(momxb + idir - 1)%sf(k, l, q + r)*fd_coeff_z_d(r, q)               
                                end do               
                            end do
                        end do
                    end do
                    !$acc end parallel loop

                    !$acc parallel loop collapse(5) gang vector default(present)
                    do q = -fd_number, p + fd_number
                        do l = 0, n
                            do k = 0, m
                                do i = 1, num_fluids
                                    do r = -fd_number, fd_number
                                        dY_dz(k, l, q, i) = dY_dz(k, l, q, i) &
                                            + j_prim_vf(i)%sf(k, l, q + r)*fd_coeff_z_d(r, q)
                                    end do
                                end do
                            end do
                        end do
                    end do
                    !$acc end parallel loop
        
                    !$acc parallel loop collapse(5) gang vector default(present)
                    do q = 0, p
                        do l = 0, n
                            do k = 0, m
                                do i = 1, num_fluids
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
                    !$acc end parallel loop

                    !Valid for any number of species
                    ! species continuity
                    !$acc parallel loop collapse(4) gang vector default(present)
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
                    !$acc end parallel loop

                    !Only valid for binary diffusion
                    !volume fraction
                    !$acc parallel loop collapse(4) gang vector default(present)
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
                    !$acc end parallel loop

                    !Valid for any number of species
                    !energy
                    !$acc parallel loop collapse(4) gang vector default(present)
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
                    !$acc end parallel loop
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
                    !$acc parallel loop collapse(5) gang vector default(present)
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
                    !$acc end parallel loop
                    
                    !$acc parallel loop collapse(5) gang vector default(present)
                    do q = 0, p
                        do l = 0, n
                            do k = -fd_number, m + fd_number
                                do i = 1, num_fluids
                                    rho_dif(k, l, q) = 0._wp
                                end do
                            end do
                        end do
                    end do
                    !$acc end parallel loop
                    
                    !$acc parallel loop collapse(4) gang vector default(present)
                    do q = 0, p
                        do l = 0, n
                            do k = -fd_number, m + fd_number
                                do i = 1, num_fluids
                                    rho_dif(k, l, q) = rho_dif(k, l, q) + alpharho_K_dif(k, l, q, i)
                                end do
                            end do
                        end do
                    end do
                    !$acc end parallel loop

                    !$acc parallel loop collapse(4) gang vector default(present)
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
                    !$acc end parallel loop

                    !$acc parallel loop collapse(5) gang vector default(present)
                    do q = 0, p
                        do l = 0, n
                            do k = -fd_number, m + fd_number
                                do i = 1, num_fluids
                                    do r = -fd_number, fd_number
                                        dY_dx(k, l, q, i) = dY_dx(k, l, q, i) &
                                            + j_prim_vf(i)%sf(k + r, l, q)*fd_coeff_x_d(r, k)
                                    end do
                                end do
                            end do
                        end do
                    end do
                    !$acc end parallel loop
        
                    !$acc parallel loop collapse(4) gang vector default(present)
                    do q = 0, p
                        do l = 0, n
                            do k = 0, m
                                do r = -fd_number, fd_number
                                    inv_denom = 1._wp / ( j_prim_vf(1)%sf(k + r, l, q)*D23 + j_prim_vf(2)%sf(k + r, l, q)*D13 + j_prim_vf(3)%sf(k + r, l, q)*D12 )
                                    do i = 1, num_fluids
                                        select case (i)
                                            case (1)
                                                dj_dx(k, l, q, i) = dj_dx(k, l, q, i) &
                                                    + rho_dif(k + r, l, q)*(D12*D13*(1._wp - j_prim_vf(1)%sf(k + r, l, q))*dY_dx(k + r, l, q, 1) &
                                                    - j_prim_vf(1)%sf(k + r, l, q)*D23*(D12*dY_dx(k + r, l, q, 2) + D13*dY_dx(k + r, l, q, 3))) &
                                                    * inv_denom*fd_coeff_x_d(r, k)

                                                djh_dx(k, l, q, i) = djh_dx(k, l, q, 1) &
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
                    !$acc end parallel loop

                    !Valid for any number of species
                    ! species continuity
                    !$acc parallel loop collapse(4) gang vector default(present)
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
                    !$acc end parallel loop

                    !volume fraction
                    !ideal gas only (3-component diffusion)
                    !$acc parallel loop collapse(3) gang vector default(present)
                    do q = 0, p
                        do l = 0, n
                            do k = 0, m                        
                                W_dif = alpha_K_dif(k, l, q, 1)*W1 + alpha_K_dif(k, l, q, 2)*W2 + alpha_K_dif(k, l, q, 3)*W3
                                W_fac = W_dif / (W1*W2*W3*rho_dif(k, l, q))
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
                    !$acc end parallel loop

                    !Valid for any number of species
                    !energy
                    !$acc parallel loop collapse(4) gang vector default(present)
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
                    !$acc end parallel loop
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