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
s_initialize_fv4_weights_dir, &
s_solve_4x4, &
s_fill_face_weights_1d_2nd, &
s_fill_face_weights_1d_4th, &
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

    real(wp), allocatable, dimension(:, :, :) :: w_interp4, w_grad4, w_interp2, w_grad2
    !$acc declare create(w_interp4, w_grad4, w_interp2, w_grad2)

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

        @:ALLOCATE(w_interp2(0:1, -1:max(m, n, p), 1:num_dims))
        @:ALLOCATE(w_grad2(0:1, -1:max(m, n, p), 1:num_dims))

        call s_fill_face_weights_1d_2nd(-1, m, 1)
        if (n > 0) then
            call s_fill_face_weights_1d_2nd(-1, n, 2)

            if (p > 0) then
                call s_fill_face_weights_1d_2nd(-1, p, 3)
            end if
        end if
        !$acc update device(w_interp2, w_grad2)

        if (dif_order == 4) then
            @:ALLOCATE(w_interp4(-1:2, -1:max(m, n, p), 1:num_dims))
            @:ALLOCATE(w_grad4(-1:2, -1:max(m, n, p), 1:num_dims))

            ! call s_fill_face_weights_1d_4th(-1, m, 1)
            ! if (n > 0) then
            !     call s_fill_face_weights_1d_4th(-1, n, 2)

            !     if (p > 0) then
            !         call s_fill_face_weights_1d_4th(-1, p, 3)
            !     end if
            ! end if

            call s_initialize_fv4_weights_dir(-1, m, 1)
            if (n > 0) then
                call s_initialize_fv4_weights_dir(-1, n, 2)

                if (p > 0) then
                    call s_initialize_fv4_weights_dir(-1, p, 3)
                end if
            end if
            !$acc update device(w_interp4, w_grad4)
        end if

    end subroutine s_initialize_diffusion_module

    subroutine s_initialize_fv4_weights_dir(ifbeg, ifend, idir)

        implicit none

        integer, intent(in) :: ifbeg, ifend, idir

        integer :: i, col, n, info
        real(wp) :: xf_stencil(0:4)
        real(wp) :: A(4,4), rhs(4)
        real(wp) :: xL, xR, dx_loc, xf

        do i = ifbeg, ifend

            !------------------------------------------------------------
            ! Build local face stencil from cb arrays
            !------------------------------------------------------------
            select case (idir)

            case (1)
                xf_stencil(0) = x_cb(i-1)
                xf_stencil(1) = x_cb(i  )
                xf_stencil(2) = x_cb(i+1)
                xf_stencil(3) = x_cb(i+2)
                xf_stencil(4) = x_cb(i+3)

            case (2)
                xf_stencil(0) = y_cb(i-1)
                xf_stencil(1) = y_cb(i  )
                xf_stencil(2) = y_cb(i+1)
                xf_stencil(3) = y_cb(i+2)
                xf_stencil(4) = y_cb(i+3)

            case (3)
                xf_stencil(0) = z_cb(i-1)
                xf_stencil(1) = z_cb(i  )
                xf_stencil(2) = z_cb(i+1)
                xf_stencil(3) = z_cb(i+2)
                xf_stencil(4) = z_cb(i+3)

            end select

            xf = xf_stencil(2)   ! target face

            !------------------------------------------------------------
            ! Build moment matrix A
            !------------------------------------------------------------
            do col = 1, 4
                xL = xf_stencil(col-1)
                xR = xf_stencil(col)
                dx_loc = xR - xL

                do n = 0, 3
                    A(n+1,col) = (xR**(n+1) - xL**(n+1)) / ((n+1._wp) * dx_loc)
                end do
            end do

            !------------------------------------------------------------
            ! Interpolation weights
            !------------------------------------------------------------
            do n = 0, 3
                rhs(n+1) = xf**n
            end do

            call s_solve_4x4(A, rhs, w_interp4(:, i, idir), info)
            if (info /= 0) call s_mpi_abort('Error computing FV interp weights')

            !------------------------------------------------------------
            ! Rebuild A (solver overwrites it)
            !------------------------------------------------------------
            do col = 1, 4
                xL = xf_stencil(col-1)
                xR = xf_stencil(col)
                dx_loc = xR - xL

                do n = 0, 3
                    A(n+1,col) = (xR**(n+1) - xL**(n+1)) / ((n+1._wp) * dx_loc)
                end do
            end do

            !------------------------------------------------------------
            ! Gradient weights
            !------------------------------------------------------------
            rhs(1) = 0._wp
            do n = 1, 3
                rhs(n+1) = real(n, wp) * xf**(n-1)
            end do

            call s_solve_4x4(A, rhs, w_grad4(:, i, idir), info)
            if (info /= 0) call s_mpi_abort('Error computing FV grad weights')

        end do

    end subroutine s_initialize_fv4_weights_dir

    subroutine s_solve_4x4(Ain, b, w, info)

        implicit none

        real(wp), intent(inout) :: Ain(4,4)
        real(wp), intent(inout) :: b(4)
        real(wp), intent(out)   :: w(-1:2)
        integer,  intent(out)   :: info

        integer :: i, j, k, piv
        real(wp) :: maxval, factor, tmp
        real(wp) :: rowtmp(4)

        info = 0

        !------------------------------------------------------------
        ! Forward elimination with partial pivoting
        !------------------------------------------------------------
        do k = 1, 4

            ! Find pivot row
            piv = k
            maxval = abs(Ain(k,k))

            do i = k+1, 4
                if (abs(Ain(i,k)) > maxval) then
                    maxval = abs(Ain(i,k))
                    piv = i
                end if
            end do

            ! Check for singular matrix
            if (maxval < 1.0e-14_wp) then
                info = 1
                return
            end if

            ! Swap rows if needed
            if (piv /= k) then
                rowtmp(:) = Ain(k,:)
                Ain(k,:)  = Ain(piv,:)
                Ain(piv,:) = rowtmp(:)

                tmp   = b(k)
                b(k)  = b(piv)
                b(piv)= tmp
            end if

            ! Eliminate below pivot
            do i = k+1, 4
                factor = Ain(i,k) / Ain(k,k)

                do j = k, 4
                    Ain(i,j) = Ain(i,j) - factor * Ain(k,j)
                end do

                b(i) = b(i) - factor * b(k)
            end do

        end do

        !------------------------------------------------------------
        ! Back substitution
        !------------------------------------------------------------
        do i = 4, 1, -1
            tmp = b(i)

            do j = i+1, 4
                tmp = tmp - Ain(i,j) * b(j)
            end do

            b(i) = tmp / Ain(i,i)
        end do

        !------------------------------------------------------------
        ! Map solution to stencil indexing (-1:2)
        !------------------------------------------------------------
        w(-1) = b(1)
        w(0) = b(2)
        w(1) = b(3)
        w(2) = b(4)

    end subroutine s_solve_4x4

    subroutine s_fill_face_weights_1d_2nd(ifbeg, ifend, idir)

        integer,  intent(in) :: ifbeg, ifend, idir

        integer :: i
        real(wp) :: xf, dx_loc

        select case (idir)
        case (1)
            do i = ifbeg, ifend

                xf     = x_cb(i)
                dx_loc = x_cc(i+1) - x_cc(i)

                w_interp2( 0, i, idir) = (x_cc(i+1) - xf) / dx_loc
                w_interp2( 1, i, idir) = (xf - x_cc(i))   / dx_loc

                w_grad2( 0, i, idir) = -1._wp / dx_loc
                w_grad2( 1, i, idir) =  1._wp / dx_loc

            end do
        
        case (2)
            do i = ifbeg, ifend

                xf     = y_cb(i)
                dx_loc = y_cc(i+1) - y_cc(i)

                w_interp2( 0, i, idir) = (y_cc(i+1) - xf) / dx_loc
                w_interp2( 1, i, idir) = (xf - y_cc(i))   / dx_loc

                w_grad2( 0, i, idir) = -1._wp / dx_loc
                w_grad2( 1, i, idir) =  1._wp / dx_loc

            end do

        case (3)
            do i = ifbeg, ifend

                xf     = z_cb(i)
                dx_loc = z_cc(i+1) - z_cc(i)

                w_interp2( 0, i, idir) = (z_cc(i+1) - xf) / dx_loc
                w_interp2( 1, i, idir) = (xf - z_cc(i))   / dx_loc

                w_grad2( 0, i, idir) = -1._wp / dx_loc
                w_grad2( 1, i, idir) =  1._wp / dx_loc

            end do

        end select

    end subroutine s_fill_face_weights_1d_2nd

    subroutine s_fill_face_weights_1d_4th(ifbeg, ifend, idir)

        integer,  intent(in) :: ifbeg, ifend, idir

        integer :: i, a, b, c
        integer :: offs(4)
        real(wp) :: x(4), xf
        real(wp) :: prod, denom

        offs = (/ -1, 0, 1, 2 /)

        select case (idir)

        case (1)

            do i = ifbeg, ifend

                xf = x_cb(i)
                x(1) = x_cc(i - 1); 
                x(2) = x_cc(i) 
                x(3) = x_cc(i + 1) 
                x(4) = x_cc(i + 2)

                ! ------------------------------------------------------------
                ! Interpolation weights:
                ! q_f = sum_{s=-1}^{2} w_interp(s,i,idir) * q(i+s)
                ! ------------------------------------------------------------
                do a = 1, 4
                    prod = 1._wp
                    do b = 1, 4
                        if (b /= a) then
                            prod = prod * (xf - x(b)) / (x(a) - x(b))
                        end if
                    end do
                    w_interp4(offs(a), i, idir) = prod
                end do

                ! w_interp4(-1, i, idir) = -1._wp / 12._wp
                ! w_interp4( 0, i, idir) =  7._wp / 12._wp
                ! w_interp4( 1, i, idir) =  7._wp / 12._wp
                ! w_interp4( 2, i, idir) = -1._wp / 12._wp

                ! ------------------------------------------------------------
                ! Gradient weights:
                ! dqdx_f = sum_{s=-1}^{2} w_grad(s,i,idir) * q(i+s)
                ! ------------------------------------------------------------
                do a = 1, 4
                    w_grad4(offs(a), i, idir) = 0._wp

                    do b = 1, 4
                        if (b /= a) then
                            denom = x(a) - x(b)
                            prod  = 1._wp / denom

                            do c = 1, 4
                                if (c /= a .and. c /= b) then
                                    prod = prod * (xf - x(c)) / (x(a) - x(c))
                                end if
                            end do

                            w_grad4(offs(a), i, idir) = w_grad4(offs(a), i, idir) + prod
                        end if
                    end do
                end do

                ! w_grad4(-1, i, idir) =  1._wp / (12._wp*(x(2) - x(1)))
                ! w_grad4( 0, i, idir) = -15._wp / (12._wp*(x(2) - x(1)))
                ! w_grad4( 1, i, idir) =  15._wp / (12._wp*(x(2) - x(1)))
                ! w_grad4( 2, i, idir) = -1._wp / (12._wp*(x(2) - x(1)))
            end do
        case (2)

            do i = ifbeg, ifend

                xf = y_cb(i)
                x(1) = y_cc(i - 1); 
                x(2) = y_cc(i) 
                x(3) = y_cc(i + 1) 
                x(4) = y_cc(i + 2)

                ! ------------------------------------------------------------
                ! Interpolation weights:
                ! q_f = sum_{s=-1}^{2} w_interp(s,i,idir) * q(i+s)
                ! ------------------------------------------------------------
                do a = 1, 4
                    prod = 1._wp
                    do b = 1, 4
                        if (b /= a) then
                            prod = prod * (xf - x(b)) / (x(a) - x(b))
                        end if
                    end do
                    w_interp4(offs(a), i, idir) = prod
                end do

                ! ------------------------------------------------------------
                ! Gradient weights:
                ! dqdx_f = sum_{s=-1}^{2} w_grad(s,i,idir) * q(i+s)
                ! ------------------------------------------------------------
                do a = 1, 4
                    w_grad4(offs(a), i, idir) = 0._wp

                    do b = 1, 4
                        if (b /= a) then
                            denom = x(a) - x(b)
                            prod  = 1._wp / denom

                            do c = 1, 4
                                if (c /= a .and. c /= b) then
                                    prod = prod * (xf - x(c)) / (x(a) - x(c))
                                end if
                            end do

                            w_grad4(offs(a), i, idir) = w_grad4(offs(a), i, idir) + prod
                        end if
                    end do
                end do
            end do
        
        case (3)

            do i = ifbeg, ifend

                xf = z_cb(i)
                x(1) = z_cc(i - 1); 
                x(2) = z_cc(i) 
                x(3) = z_cc(i + 1) 
                x(4) = z_cc(i + 2)

                ! ------------------------------------------------------------
                ! Interpolation weights:
                ! q_f = sum_{s=-1}^{2} w_interp(s,i,idir) * q(i+s)
                ! ------------------------------------------------------------
                do a = 1, 4
                    prod = 1._wp
                    do b = 1, 4
                        if (b /= a) then
                            prod = prod * (xf - x(b)) / (x(a) - x(b))
                        end if
                    end do
                    w_interp4(offs(a), i, idir) = prod
                end do

                ! ------------------------------------------------------------
                ! Gradient weights:
                ! dqdx_f = sum_{s=-1}^{2} w_grad(s,i,idir) * q(i+s)
                ! ------------------------------------------------------------
                do a = 1, 4
                    w_grad4(offs(a), i, idir) = 0._wp

                    do b = 1, 4
                        if (b /= a) then
                            denom = x(a) - x(b)
                            prod  = 1._wp / denom

                            do c = 1, 4
                                if (c /= a .and. c /= b) then
                                    prod = prod * (xf - x(c)) / (x(a) - x(c))
                                end if
                            end do

                            w_grad4(offs(a), i, idir) = w_grad4(offs(a), i, idir) + prod
                        end if
                    end do
                end do
            end do

        end select


    end subroutine s_fill_face_weights_1d_4th

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
        real(wp) :: S_1, S_2, S_3 !< Shock sensor values
        real(wp) :: grid_spacing
        real(wp) :: rho_L, rho_LL, rho_R, rho_RR, rho_f, rhog_f
        real(wp) :: alpha_m_L, alpha_m_LL, alpha_m_R, alpha_m_RR, alpha_m_f
        real(wp) :: n_gate
        real(wp) :: g_f
        real(wp) :: P_L, P_LL, P_R, P_RR, P_f
        real(wp) :: T_f, W_f
        real(wp) :: sum_jflux
        integer, dimension(3) :: offsets

        real(wp) :: alpharho_L(Dif_size), alpharho_LL(Dif_size), alpharho_R(Dif_size), alpharho_RR(Dif_size), alpharho_f(Dif_size)
        real(wp) :: Y_f(Dif_size), Y_L(Dif_size), Y_LL(Dif_size), Y_R(Dif_size), Y_RR(Dif_size)
        real(wp) :: dY_ds_f(Dif_size)
        real(wp) :: h_f(Dif_size)
        real(wp) :: j_flux(Dif_size)

        R_univ = 8314.462618_wp
        n_gate = 2.0_wp

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

                        ! ! Calculate grid spacing using direction-based indexing
                        ! select case (idir)
                        ! case (1)
                        !     grid_spacing = x_cc(k + 1) - x_cc(k)
                        ! case (2)
                        !     grid_spacing = y_cc(l + 1) - y_cc(l)
                        ! case (3)
                        !     grid_spacing = z_cc(q + 1) - z_cc(q)
                        ! end select

                        ! Calculate grid spacing using direction-based indexing
                        select case (idir)
                        case (1)
                            r = k
                        case (2)
                            r = l
                        case (3)
                            r = q
                        end select

                        alpha_m_L = q_prim_vf(advg_idx)%sf(k, l, q)
                        alpha_m_R = q_prim_vf(advg_idx)%sf(k + offsets(1), l + offsets(2), q + offsets(3))
                        alpha_m_LL = q_prim_vf(advg_idx)%sf(k - offsets(1), l - offsets(2), q - offsets(3))
                        alpha_m_RR = q_prim_vf(advg_idx)%sf(k + 2*offsets(1), l + 2*offsets(2), q + 2*offsets(3))

                        if (alpha_m_L < small_num_dif .or. alpha_m_R < small_num_dif) cycle
                        
                        if (alpha_m_LL > small_num_dif .and. alpha_m_RR > small_num_dif .and. dif_order == 4) then

                            do i = 1, Dif_size
                                alpharho_LL(i) = q_prim_vf(Dif_idx(i))%sf(k - offsets(1), l - offsets(2), q - offsets(3))
                                alpharho_L(i) = q_prim_vf(Dif_idx(i))%sf(k, l, q)
                                alpharho_R(i) = q_prim_vf(Dif_idx(i))%sf(k + offsets(1), l + offsets(2), q + offsets(3))
                                alpharho_RR(i) = q_prim_vf(Dif_idx(i))%sf(k + 2*offsets(1), l + 2*offsets(2), q + 2*offsets(3))
                                alpharho_f(i) = w_interp4(-1, r, idir)*alpharho_LL(i) + w_interp4(0, r, idir)*alpharho_L(i) + &
                                                w_interp4(1, r, idir)*alpharho_R(i) + w_interp4(2, r, idir)*alpharho_RR(i)
                            end do

                            alpha_m_f = w_interp4(-1, r, idir)*alpha_m_LL + w_interp4(0, r, idir)*alpha_m_L + &
                                        w_interp4(1, r, idir)*alpha_m_R + w_interp4(2, r, idir)*alpha_m_RR

                            P_LL = q_prim_vf(E_idx)%sf(k - offsets(1), l - offsets(2), q - offsets(3))
                            P_L = q_prim_vf(E_idx)%sf(k, l, q)
                            P_R = q_prim_vf(E_idx)%sf(k + offsets(1), l + offsets(2), q + offsets(3))
                            P_RR = q_prim_vf(E_idx)%sf(k + 2*offsets(1), l + 2*offsets(2), q + 2*offsets(3))

                            S_1 = abs(P_R - P_L) / max(P_L, P_R, small_num_dif)
                            S_2 = abs(P_RR - P_R) / max(P_R, P_RR, small_num_dif)
                            S_3 = abs(P_L - P_LL) / max(P_L, P_LL, small_num_dif)
                            if (max(S_1, S_2, S_3) > 0.05_wp) cycle !crude shock sensor (dont calculate diffusion across shocks)
                            P_f = w_interp4(-1, r, idir)*P_LL + w_interp4(0, r, idir)*P_L + &
                                  w_interp4(1, r, idir)*P_R + w_interp4(2, r, idir)*P_RR
                            rho_LL = 0.0_wp
                            rho_L = 0.0_wp
                            rho_R = 0.0_wp
                            rho_RR = 0.0_wp
                            rho_f = 0.0_wp

                            do i = 1, Dif_size
                                rho_LL = rho_LL + alpharho_LL(i)
                                rho_L = rho_L + alpharho_L(i)
                                rho_R = rho_R + alpharho_R(i)
                                rho_RR = rho_RR + alpharho_RR(i)
                                rho_f = rho_f + alpharho_f(i)
                            end do

                            do i = 1, Dif_size
                                Y_LL(i) = alpharho_LL(i) / rho_LL
                                Y_L(i) = alpharho_L(i) / rho_L
                                Y_R(i) = alpharho_R(i) / rho_R
                                Y_RR(i) = alpharho_RR(i) / rho_RR
                                Y_f(i) = alpharho_f(i) / rho_f
                            end do
                            
                            do i = 1, Dif_size
                                dY_ds_f(i) = w_grad4(-1, r, idir)*Y_LL(i) + w_grad4(0, r, idir)*Y_L(i) + &
                                     w_grad4(1, r, idir)*Y_R(i) + w_grad4(2, r, idir)*Y_RR(i)
                            end do

                        else !use 2nd order

                            do i = 1, Dif_size
                                alpharho_L(i) = q_prim_vf(Dif_idx(i))%sf(k, l, q)
                                alpharho_R(i) = q_prim_vf(Dif_idx(i))%sf(k + offsets(1), l + offsets(2), q + offsets(3))
                                alpharho_f(i) = w_interp2(0, r, idir)*alpharho_L(i) + w_interp2(1, r, idir)*alpharho_R(i)
                            end do

                            alpha_m_f = w_interp2(0, r, idir)*alpha_m_L + w_interp2(1, r, idir)*alpha_m_R
                            P_L = q_prim_vf(E_idx)%sf(k, l, q)
                            P_R = q_prim_vf(E_idx)%sf(k + offsets(1), l + offsets(2), q + offsets(3))

                            if (abs(P_R - P_L) / max(P_L, P_R, small_num_dif) > 0.05_wp ) cycle !crude shock sensor (dont calculate diffusion across shocks)
                            P_f = w_interp2(0, r, idir)*P_L + w_interp2(1, r, idir)*P_R

                            rho_L = 0.0_wp
                            rho_R = 0.0_wp
                            rho_f = 0.0_wp
                            do i = 1, Dif_size
                                rho_L = rho_L + alpharho_L(i)
                                rho_R = rho_R + alpharho_R(i)
                                rho_f = rho_f + alpharho_f(i)
                            end do

                            do i = 1, Dif_size
                                Y_L(i) = alpharho_L(i) / rho_L
                                Y_R(i) = alpharho_R(i) / rho_R
                                Y_f(i) = alpharho_f(i) / rho_f
                            end do
                            
                            do i = 1, Dif_size
                                dY_ds_f(i) = w_grad2(0, r, idir)*Y_L(i) + w_grad2(1, r, idir)*Y_R(i)
                            end do
                        end if

                        g_f = 2.0_wp * (alpha_m_L**n_gate) * (alpha_m_R**n_gate) / ( (alpha_m_L**n_gate) + (alpha_m_R**n_gate) )
                        ! g_f = 1.0_wp
                        ! Total gas density at face
                        rhog_f = rho_f / alpha_m_f
                        ! rhog_f = 101325.0_wp * 28.02_wp / (R_univ * 298.0_wp)

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

    subroutine s_correct_riemann_volume_fractions(q_rs_vf, bounds1, bounds2, bounds3)

        type(int_bounds_info), intent(in) :: bounds1, bounds2, bounds3
        real(wp), dimension(bounds1%beg:, bounds2%beg:, bounds3%beg:, 1:), intent(inout) :: q_rs_vf
        type(int_bounds_info) :: weno_bounds1, weno_bounds2, weno_bounds3

        integer :: x, y, z, i
        real(wp) :: rho, W
        real(wp) :: alpharho(Dif_size), Y_s(Dif_size)

        weno_bounds1%beg = bounds1%beg + weno_polyn
        weno_bounds1%end = bounds1%end - weno_polyn
        weno_bounds2 = bounds2
        weno_bounds3 = bounds3


        
        do z = weno_bounds3%beg, weno_bounds3%end
            do y = weno_bounds2%beg, weno_bounds2%end
                do x = weno_bounds1%beg, weno_bounds1%end
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
        @:DEALLOCATE(w_interp2, w_grad2)
        if (dif_order == 4) then
            @:DEALLOCATE(w_interp4, w_grad4)
        end if

    end subroutine s_finalize_diffusion_module


end module m_diffusion