!>
!! @file m_diffusion.f90
!! @brief Contains module m_diffusion

#:include 'macros.fpp'

!> @brief This module is used to compute flux terms for binary diffusion
module m_diffusion

    ! Dependencies =============================================================

    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters

    use m_mpi_proxy            !< Message passing interface (MPI) module proxy

    use m_weno                 !< WENO module

    use m_helper              !< Helper functions

    ! ==========================================================================

    private; public s_get_diffusion, & 
s_compute_fd_gradient_diffusion, &
s_apply_scalar_divergence_theorem_diffusion, &
s_reconstruct_cell_boundary_values_diff, & 
s_initialize_diffusion_module

    type(int_bounds_info) :: iv
    type(int_bounds_info) :: is1, is2, is3
 !$acc declare create(is1, is2, is3, iv)   
    
    type(int_bounds_info) :: iv
    type(int_bounds_info) :: is1_diffusion, is2_diffusion, is3_diffusion
    !$acc declare create(is1_diffusion, is2_diffusion, is3_diffusion, iv)

    real(kind(0d0)), allocatable, dimension(:) :: Ds_diffusion
    !$acc declare create(Ds_diffusion)

contains

    subroutine s_initialize_diffusion_module

        integer :: i !< generic loop iterators

        @:ALLOCATE(Ds_diffusion(1:Dif_size))

        
        do i = 1, Dif_size
            Ds_diffusion(i) = fluid_pp(Dif_idx(i))%D
        end do
        
        !$acc update device(Ds_diffusion, Dif_idx, Dif_size)
        !$acc enter data copyin(is1_diffusion, is2_diffusion, is3_diffusion, iv)

    end subroutine s_initialize_diffusion_module
    

    contains

    !>  Computes diffusion terms
    !!  @param q_cons_vf Cell-averaged conservative variables
    !!  @param q_prim_vf Cell-averaged primitive variables
    !!  @param rhs_vf Cell-averaged RHS variables
    subroutine s_get_diffusion(jL_prim_rsx_vf, jL_prim_rsy_vf, jL_prim_rsz_vf, &
                             djL_prim_dx_n, djL_prim_dy_n, djL_prim_dz_n, &
                             jL_prim, & 
                             jR_prim_rsx_vf, jR_prim_rsy_vf, jR_prim_rsz_vf, &
                             djR_prim_dx_n, djR_prim_dy_n, djR_prim_dz_n, &
                             jR_prim, &
                             j_vf_qp, &
                             dj_prim_dx_qp, dj_prim_dy_qp, dj_prim_dz_qp,  &
                             ix, iy, iz)

        real(kind(0d0)), dimension(startx:, starty:, startz:, 1:), &
             intent(INOUT) :: jL_prim_rsx_vf, jR_prim_rsx_vf, &
                              jL_prim_rsy_vf, jR_prim_rsy_vf, &
                              jL_prim_rsz_vf, jR_prim_rsz_vf

        type(vector_field), dimension(1:num_dims) :: jL_prim, jR_prim

        type(vector_field) :: j_vf_qp

        type(vector_field), dimension(1), &
            intent(INOUT) :: djL_prim_dx_n, djR_prim_dx_n, &
                             djL_prim_dy_n, djR_prim_dy_n, &
                             djL_prim_dz_n, djR_prim_dz_n

        type(vector_field) :: dj_prim_dx_qp, dj_prim_dy_qp, dj_prim_dz_qp

        type(int_bounds_info), intent(IN) :: ix, iy, iz

        integer :: i, j, k, l

        do i = 1, num_dims

            iv%beg = cont_idx%beg; iv%end = cont_idx%beg

            !$acc update device(iv)

            call s_reconstruct_cell_boundary_values_diff( &
                j_vf_qp%vf(iv%beg:iv%end), &
                jL_prim_rsx_vf, jL_prim_rsy_vf, jL_prim_rsz_vf, &
                jR_prim_rsx_vf, jR_prim_rsy_vf, jR_prim_rsz_vf, &
                i, jL_prim(i)%vf(iv%beg:iv%end), jR_prim(i)%vf(iv%beg:iv%end), &
                ix, iy, iz)
        end do

        if (weno_Re_flux) then
            ! Compute mass fraction gradient at cell centers using scalar
            ! divergence theorem
            do i = 1, num_dims
                if (i == 1) then
                    call s_apply_scalar_divergence_theorem_diffusion( &
                        jL_prim(i)%vf(iv%beg:iv%end), &
                        jR_prim(i)%vf(iv%beg:iv%end), &
                        dj_prim_dx_qp%vf(iv%beg:iv%end), i, &
                        ix, iy, iz, iv, dx, m, buff_size)
                elseif (i == 2) then
                    call s_apply_scalar_divergence_theorem_diffusion( &
                        jL_prim(i)%vf(iv%beg:iv%end), &
                        jR_prim(i)%vf(iv%beg:iv%end), &
                        dj_prim_dy_qp%vf(iv%beg:iv%end), i, &
                        ix, iy, iz, iv, dy, n, buff_size)
                else
                    call s_apply_scalar_divergence_theorem_diffusion( &
                        jL_prim(i)%vf(iv%beg:iv%end), &
                        jR_prim(i)%vf(iv%beg:iv%end), &
                        dj_prim_dz_qp%vf(iv%beg:iv%end), i, &
                        ix, iy, iz, iv, dz, p, buff_size)
                end if
            end do

        else ! Compute velocity gradient at cell centers using finite differences

            iv%beg = cont_idx%beg; iv%end = cont_idx%beg
            !$acc update device(iv)

    !$acc parallel loop collapse(3) gang vector default(present)
            do l = iz%beg, iz%end
                do k = iy%beg, iy%end
                    do j = ix%beg + 1, ix%end
    !$acc loop seq
                        do i = iv%beg, iv%end
                            djL_prim_dx_n(1)%vf(i)%sf(j, k, l) = &
                                (j_vf_qp%vf(i)%sf(j, k, l) - &
                                j_vf_qp%vf(i)%sf(j - 1, k, l))/ &
                                (x_cc(j) - x_cc(j - 1))
                        end do
                    end do
                end do
            end do

    !$acc parallel loop collapse(3) gang vector default(present)
            do l = iz%beg, iz%end
                do k = iy%beg, iy%end
                    do j = ix%beg, ix%end - 1
    !$acc loop seq
                        do i = iv%beg, iv%end
                            djR_prim_dx_n(1)%vf(i)%sf(j, k, l) = &
                                (j_vf_qp%vf(i)%sf(j + 1, k, l) - &
                                j_vf_qp%vf(i)%sf(j, k, l))/ &
                                (x_cc(j + 1) - x_cc(j))
                        end do
                    end do
                end do
            end do

            if (n > 0) then

    !$acc parallel loop collapse(3) gang vector default(present)
                do l = iz%beg, iz%end
                    do j = iy%beg + 1, iy%end
                        do k = ix%beg, ix%end
    !$acc loop seq
                            do i = iv%beg, iv%end
                                djL_prim_dy_n(1)%vf(i)%sf(k, j, l) = &
                                    (j_vf_qp%vf(i)%sf(k, j, l) - &
                                    j_vf_qp%vf(i)%sf(k, j - 1, l))/ &
                                    (y_cc(j) - y_cc(j - 1))
                            end do
                        end do
                    end do
                end do

    !$acc parallel loop collapse(3) gang vector default(present)
                do l = iz%beg, iz%end
                    do j = iy%beg, iy%end - 1
                        do k = ix%beg, ix%end
                            !$acc loop seq
                            do i = iv%beg, iv%end
                                djR_prim_dy_n(1)%vf(i)%sf(k, j, l) = &
                                    (j_vf_qp%vf(i)%sf(k, j + 1, l) - &
                                    j_vf_qp%vf(i)%sf(k, j, l))/ &
                                    (y_cc(j + 1) - y_cc(j))
                            end do
                        end do
                    end do
                end do

   
    
    

                if (p > 0) then

    !$acc parallel loop collapse(3) gang vector default(present)
                    do j = iz%beg + 1, iz%end
                        do l = iy%beg, iy%end
                            do k = ix%beg, ix%end
    !$acc loop seq
                                do i = iv%beg, iv%end

                                    djL_prim_dz_n(1)%vf(i)%sf(k, l, j) = &
                                        (j_vf_qp%vf(i)%sf(k, l, j) - &
                                        j_vf_qp%vf(i)%sf(k, l, j - 1))/ &
                                        (z_cc(j) - z_cc(j - 1))
                                end do
                            end do
                        end do
                    end do

    !$acc parallel loop collapse(3) gang vector default(present)
                    do j = iz%beg, iz%end - 1
                        do l = iy%beg, iy%end
                            do k = ix%beg, ix%end
    !$acc loop seq
                                do i = iv%beg, iv%end

                                    djR_prim_dz_n(1)%vf(i)%sf(k, l, j) = &
                                        (j_vf_qp%vf(i)%sf(k, l, j + 1) - &
                                        j_vf_qp%vf(i)%sf(k, l, j))/ &
                                        (z_cc(j + 1) - z_cc(j))
                                end do
                            end do
                        end do
                    end do

                    do i = iv%beg, iv%end
                        call s_compute_fd_gradient_diffusion(j_vf_qp%vf(i), &
                                                dj_prim_dx_qp%vf(i), &
                                                dj_prim_dy_qp%vf(i), &
                                                dj_prim_dz_qp%vf(i), &
                                                ix, iy, iz, buff_size)
                    end do

                else

                    do i = iv%beg, iv%end
                        call s_compute_fd_gradient_diffusion(j_vf_qp%vf(i), &
                                                dj_prim_dx_qp%vf(i), &
                                                dj_prim_dy_qp%vf(i), &
                                                dj_prim_dy_qp%vf(i), &
                                                ix, iy, iz, buff_size)
                    end do

                end if

            else
                do i = iv%beg, iv%end
                    call s_compute_fd_gradient_diffusion(j_vf_qp%vf(i), &
                                            dj_prim_dx_qp%vf(i), &
                                            dj_prim_dx_qp%vf(i), &
                                            dj_prim_dx_qp%vf(i), &
                                            ix, iy, iz, buff_size)
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
    subroutine s_compute_fd_gradient_diffusion(var, grad_x, grad_y, grad_z, &
                                     ix, iy, iz, buff_size_in)

        type(scalar_field), intent(IN) :: var
        type(scalar_field), intent(INOUT) :: grad_x
        type(scalar_field), intent(INOUT) :: grad_y
        type(scalar_field), intent(INOUT) :: grad_z

        integer, intent(IN) :: buff_size_in

        integer :: j, k, l !< Generic loop iterators

        type(int_bounds_info) :: ix, iy, iz

        ix%beg = -buff_size_in; ix%end = m + buff_size_in; 
        if (n > 0) then
            iy%beg = -buff_size_in; iy%end = n + buff_size_in
        else
            iy%beg = -1; iy%end = 1
        end if

        if (p > 0) then
            iz%beg = -buff_size_in; iz%end = p + buff_size_in
        else
            iz%beg = -1; iz%end = 1
        end if

        !$acc update device(ix, iy, iz)

    !$acc parallel loop collapse(3) gang vector default(present)
        do l = iz%beg + 1, iz%end - 1
            do k = iy%beg + 1, iy%end - 1
                do j = ix%beg + 1, ix%end - 1
                    grad_x%sf(j, k, l) = &
                        (var%sf(j + 1, k, l) - var%sf(j - 1, k, l))/ &
                        (x_cc(j + 1) - x_cc(j - 1))
                end do
            end do
        end do

        if (n > 0) then
    !$acc parallel loop collapse(3) gang vector
            do l = iz%beg + 1, iz%end - 1
                do k = iy%beg + 1, iy%end - 1
                    do j = ix%beg + 1, ix%end - 1
                        grad_y%sf(j, k, l) = &
                            (var%sf(j, k + 1, l) - var%sf(j, k - 1, l))/ &
                            (y_cc(k + 1) - y_cc(k - 1))
                    end do
                end do
            end do
        end if

        if (p > 0) then
    !$acc parallel loop collapse(3) gang vector
            do l = iz%beg + 1, iz%end - 1
                do k = iy%beg + 1, iy%end - 1
                    do j = ix%beg + 1, ix%end - 1
                        grad_z%sf(j, k, l) = &
                            (var%sf(j, k, l + 1) - var%sf(j, k, l - 1))/ &
                            (z_cc(l + 1) - z_cc(l - 1))
                    end do
                end do
            end do
        end if

        ix%beg = -buff_size_in; ix%end = m + buff_size_in; 
        if (n > 0) then
            iy%beg = -buff_size_in; iy%end = n + buff_size_in
        else
            iy%beg = 0; iy%end = 0
        end if

        if (p > 0) then
            iz%beg = -buff_size_in; iz%end = p + buff_size_in
        else
            iz%beg = 0; iz%end = 0
        end if

        !$acc update device(ix, iy, iz)

    !$acc parallel loop collapse(2) gang vector default(present)
        do l = iz%beg, iz%end
            do k = iy%beg, iy%end
                grad_x%sf(ix%beg, k, l) = &
                    (-3d0*var%sf(ix%beg, k, l) + 4d0*var%sf(ix%beg + 1, k, l) - var%sf(ix%beg + 2, k, l))/ &
                    (x_cc(ix%beg + 2) - x_cc(ix%beg))
                grad_x%sf(ix%end, k, l) = &
                    (3d0*var%sf(ix%end, k, l) - 4d0*var%sf(ix%end - 1, k, l) + var%sf(ix%end - 2, k, l))/ &
                    (x_cc(ix%end) - x_cc(ix%end - 2))
            end do
        end do
        if (n > 0) then
    !$acc parallel loop collapse(2) gang vector default(present)
            do l = iz%beg, iz%end
                do j = ix%beg, ix%end
                    grad_y%sf(j, iy%beg, l) = &
                        (-3d0*var%sf(j, iy%beg, l) + 4d0*var%sf(j, iy%beg + 1, l) - var%sf(j, iy%beg + 2, l))/ &
                        (y_cc(iy%beg + 2) - y_cc(iy%beg))
                    grad_y%sf(j, iy%end, l) = &
                        (3d0*var%sf(j, iy%end, l) - 4d0*var%sf(j, iy%end - 1, l) + var%sf(j, iy%end - 2, l))/ &
                        (y_cc(iy%end) - y_cc(iy%end - 2))
                end do
            end do
            if (p > 0) then
    !$acc parallel loop collapse(2) gang vector default(present)
                do k = iy%beg, iy%end
                    do j = ix%beg, ix%end
                        grad_z%sf(j, k, iz%beg) = &
                            (-3d0*var%sf(j, k, iz%beg) + 4d0*var%sf(j, k, iz%beg + 1) - var%sf(j, k, iz%beg + 2))/ &
                            (z_cc(iz%beg + 2) - z_cc(iz%beg))
                        grad_z%sf(j, k, iz%end) = &
                            (3d0*var%sf(j, k, iz%end) - 4d0*var%sf(j, k, iz%end - 1) + var%sf(j, k, iz%end - 2))/ &
                            (z_cc(iz%end) - z_cc(iz%end - 2))
                    end do
                end do
            end if
        end if

        if (bc_x%beg <= -3) then
    !$acc parallel loop collapse(2) gang vector default(present)
            do l = iz%beg, iz%end
                do k = iy%beg, iy%end
                    grad_x%sf(0, k, l) = (-3d0*var%sf(0, k, l) + 4d0*var%sf(1, k, l) - var%sf(2, k, l))/ &
                                        (x_cc(2) - x_cc(0))
                end do
            end do
        end if
        if (bc_x%end <= -3) then
    !$acc parallel loop collapse(2) gang vector default(present)
            do l = iz%beg, iz%end
                do k = iy%beg, iy%end
                    grad_x%sf(m, k, l) = (3d0*var%sf(m, k, l) - 4d0*var%sf(m - 1, k, l) + var%sf(m - 2, k, l))/ &
                                        (x_cc(m) - x_cc(m - 2))
                end do
            end do
        end if
        if (n > 0) then
            if (bc_y%beg <= -3 .and. bc_y%beg /= -14) then
    !$acc parallel loop collapse(2) gang vector default(present)
                do l = iz%beg, iz%end
                    do j = ix%beg, ix%end
                        grad_y%sf(j, 0, l) = (-3d0*var%sf(j, 0, l) + 4d0*var%sf(j, 1, l) - var%sf(j, 2, l))/ &
                                            (y_cc(2) - y_cc(0))
                    end do
                end do
            end if
            if (bc_y%end <= -3) then
    !$acc parallel loop collapse(2) gang vector default(present)
                do l = iz%beg, iz%end
                    do j = ix%beg, ix%end
                        grad_y%sf(j, n, l) = (3d0*var%sf(j, n, l) - 4d0*var%sf(j, n - 1, l) + var%sf(j, n - 2, l))/ &
                                            (y_cc(n) - y_cc(n - 2))
                    end do
                end do
            end if
            if (p > 0) then
                if (bc_z%beg <= -3) then
    !$acc parallel loop collapse(2) gang vector default(present)
                    do k = iy%beg, iy%end
                        do j = ix%beg, ix%end
                            grad_z%sf(j, k, 0) = &
                                (-3d0*var%sf(j, k, 0) + 4d0*var%sf(j, k, 1) - var%sf(j, k, 2))/ &
                                (z_cc(2) - z_cc(0))
                        end do
                    end do
                end if
                if (bc_z%end <= -3) then
    !$acc parallel loop collapse(2) gang vector default(present)
                    do k = iy%beg, iy%end
                        do j = ix%beg, ix%end
                            grad_z%sf(j, k, p) = &
                                (3d0*var%sf(j, k, p) - 4d0*var%sf(j, k, p - 1) + var%sf(j, k, p - 2))/ &
                                (z_cc(p) - z_cc(p - 2))
                        end do
                    end do
                end if
            end if
        end if

    end subroutine s_compute_fd_gradient_diffusion ! --------------------------------------

            !>  The purpose of this subroutine is to employ the inputted
        !!      left and right cell-boundary integral-averaged variables
        !!      to compute the relevant cell-average first-order spatial
        !!      derivatives in the x-, y- or z-direction by means of the
        !!      scalar divergence theorem.
        !!  @param vL_vf Left cell-boundary integral averages
        !!  @param vR_vf Right cell-boundary integral averages
        !!  @param dv_ds_vf Cell-average first-order spatial derivatives
        !!  @param norm_dir Splitting coordinate direction
    subroutine s_apply_scalar_divergence_theorem_diffusion(vL_vf, vR_vf, & ! --------
                                                 dv_ds_vf, &
                                                 norm_dir, &
                                                 ix, iy, iz, iv, &
                                                 dL, dim, buff_size_in)

        type(int_bounds_info) :: ix, iy, iz, iv
            
        integer :: buff_size_in, dim

        real(kind(0d0)), dimension(-buff_size_in:dim + buff_size_in) :: dL
        ! arrays of cell widths

        type(scalar_field), &
            dimension(iv%beg:iv%end), &
            intent(IN) :: vL_vf, vR_vf

        type(scalar_field), &
            dimension(iv%beg:iv%end), &
            intent(INOUT) :: dv_ds_vf

        integer, intent(IN) :: norm_dir

        integer :: i, j, k, l !< Generic loop iterators

        !$acc update device(ix, iy, iz, iv)

        ! First-Order Spatial Derivatives in x-direction ===================
        if (norm_dir == 1) then

            ! A general application of the scalar divergence theorem that
            ! utilizes the left and right cell-boundary integral-averages,
            ! inside each cell, or an arithmetic mean of these two at the
            ! cell-boundaries, to calculate the cell-averaged first-order
            ! spatial derivatives inside the cell.

!$acc parallel loop collapse(3) gang vector default(present)
            do l = iz%beg, iz%end
                do k = iy%beg, iy%end
                    do j = ix%beg + 1, ix%end - 1
!$acc loop seq
                        do i = iv%beg, iv%end
                            dv_ds_vf(i)%sf(j,k,l) = &
                                        1d0/((1d0+wa_flg)*dL(j)) &
                                       *( wa_flg*vL_vf(i)%sf(j+1,k,l) &
                                        +        vR_vf(i)%sf( j ,k,l) &
                                        -        vL_vf(i)%sf( j ,k,l) &
                                        - wa_flg*vR_vf(i)%sf(j-1,k,l) )
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

            do l = iz%beg, iz%end
                do k = iy%beg + 1, iy%end - 1
                    do j = ix%beg, ix%end
!$acc loop seq
                        do i = iv%beg, iv%end
                            dv_ds_vf(i)%sf(j,k,l) = &
                                             1d0/((1d0+wa_flg)*dL(k)) &
                                       *( wa_flg*vL_vf(i)%sf(j,k+1,l) &
                                        +        vR_vf(i)%sf(j, k ,l) &
                                        -        vL_vf(i)%sf(j, k ,l) &
                                        - wa_flg*vR_vf(i)%sf(j,k-1,l) )
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
            do l = iz%beg + 1, iz%end - 1
                do k = iy%beg, iy%end
                    do j = ix%beg, ix%end
!$acc loop seq
                        do i = iv%beg, iv%end
                            dv_ds_vf(i)%sf(j,k,l) = &
                                             1d0/((1d0+wa_flg)*dL(l)) &
                                       *( wa_flg*vL_vf(i)%sf(j,k,l+1) &
                                        +        vR_vf(i)%sf(j,k, l ) &
                                        -        vL_vf(i)%sf(j,k, l ) &
                                        - wa_flg*vR_vf(i)%sf(j,k,l-1) )
                        end do
                    end do
                end do
            end do

        end if
        ! END: First-Order Spatial Derivatives in z-direction ==============

    end subroutine s_apply_scalar_divergence_theorem_diffusion ! ---------------------


    subroutine s_reconstruct_cell_boundary_values_diff(v_vf, vL_x, vL_y, vL_z, vR_x, vR_y, vR_z, & ! -
                                                       norm_dir, vL_prim_vf, vR_prim_vf, ix, iy, iz)

        type(scalar_field), dimension(iv%beg:iv%end), intent(IN) :: v_vf
        type(scalar_field), dimension(iv%beg:iv%end), intent(INOUT) :: vL_prim_vf, vR_prim_vf

        real(kind(0d0)), dimension(startx:, starty:, startz:, 1:), intent(INOUT) :: vL_x, vL_y, vL_z, vR_x, vR_y, vR_z 

        integer, intent(IN) :: norm_dir

        integer :: weno_dir !< Coordinate direction of the WENO reconstruction

        integer :: i, j, k, l

        type(int_bounds_info) :: ix, iy, iz
        ! Reconstruction in s1-direction ===================================

        if (norm_dir == 1) then
            is1 = ix; is2 = iy; is3 = iz
            weno_dir = 1; is1%beg = is1%beg + weno_polyn
            is1%end = is1%end - weno_polyn

        elseif (norm_dir == 2) then
            is1 = iy; is2 = ix; is3 = iz
            weno_dir = 2; is1%beg = is1%beg + weno_polyn
            is1%end = is1%end - weno_polyn

        else
            is1 = iz; is2 = iy; is3 = ix
            weno_dir = 3; is1%beg = is1%beg + weno_polyn
            is1%end = is1%end - weno_polyn

        end if

        !$acc update device(is1, is2, is3, iv)

        if (n > 0) then
            if (p > 0) then

                call s_weno(v_vf(iv%beg:iv%end), &
                    vL_x(:, :, :, iv%beg:iv%end), vL_y(:, :, :, iv%beg:iv%end), vL_z(:, :, :, iv%beg:iv%end), vR_x(:, :, :, iv%beg:iv%end), vR_y(:, :, :, iv%beg:iv%end), vR_z(:, :, :, iv%beg:iv%end), &
                                norm_dir, weno_dir, &
                                is1, is2, is3)
            else
                call s_weno(v_vf(iv%beg:iv%end), &
                    vL_x(:, :, :, iv%beg:iv%end), vL_y(:, :, :, iv%beg:iv%end), vL_z(:, :, :, :), vR_x(:, :, :, iv%beg:iv%end), vR_y(:, :, :, iv%beg:iv%end), vR_z(:, :, :, :), &
                                norm_dir, weno_dir, &
                                is1, is2, is3)
            end if
        else

            call s_weno(v_vf(iv%beg:iv%end), &
                        vL_x(:, :, :, iv%beg:iv%end), vL_y(:, :, :, :), vL_z(:, :, :, :), vR_x(:, :, :, iv%beg:iv%end), vR_y(:, :, :, :), vR_z(:, :, :, :), &
                            norm_dir, weno_dir, &
                            is1, is2, is3)
        end if

        if (diffusion) then
            if (weno_Re_flux) then
                if (norm_dir == 2) then
!$acc parallel loop collapse(4) gang vector default(present)
                    do i = iv%beg, iv%end
                        do l = is3%beg, is3%end
                            do j = is1%beg, is1%end
                                do k = is2%beg, is2%end
                                    vL_prim_vf(i)%sf(k, j, l) = vL_y(j, k, l, i)
                                    vR_prim_vf(i)%sf(k, j, l) = vR_y(j, k, l, i)
                                end do
                            end do
                        end do
                    end do
                elseif (norm_dir == 3) then
!$acc parallel loop collapse(4) gang vector default(present)
                    do i = iv%beg, iv%end
                        do j = is1%beg, is1%end
                            do k = is2%beg, is2%end
                                do l = is3%beg, is3%end
                                    vL_prim_vf(i)%sf(l, k, j) = vL_z(j, k, l, i)
                                    vR_prim_vf(i)%sf(l, k, j) = vR_z(j, k, l, i)
                                end do
                            end do
                        end do
                    end do
                elseif (norm_dir == 1) then
!$acc parallel loop collapse(4) gang vector default(present)
                    do i = iv%beg, iv%end
                        do l = is3%beg, is3%end
                            do k = is2%beg, is2%end
                                do j = is1%beg, is1%end
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

    end subroutine s_reconstruct_cell_boundary_values_diff ! --------------------


end module m_diffusion