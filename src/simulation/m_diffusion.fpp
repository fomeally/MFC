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

    use m_constants, only: pi

    use m_precision_select, only: mpi_p

#ifdef MFC_MPI
    use mpi                    !< Message passing interface (MPI) module
#endif

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
s_write_diffusion_tbc_restart, &
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
    !$acc declare create(cps)

    real(wp), allocatable, dimension(:, :, :, :) :: dj_dx, dj_dy, dj_dz, djh_dx, djh_dy, djh_dz, dY_dx, dY_dy, dY_dz, alpha_K_dif, alpharho_K_dif, Y_dif, h_dif
    !$acc declare create(dj_dx, dj_dy, dj_dz, djh_dx, djh_dy, djh_dz, dY_dx, dY_dy, dY_dz, alpha_K_dif, alpharho_K_dif, Y_dif, h_dif)

    real(wp), allocatable, dimension(:, :, :) :: rho_dif, alpha_dif, dvel_dx, dvel_dy, dvel_dz, denom, rhogcg2, rho1c12, kdivu, W_dif, T_dif
    !$acc declare create(rho_dif, alpha_dif, dvel_dx, dvel_dy, dvel_dz, denom, rhogcg2, rho1c12, kdivu, W_dif, T_dif)

    real(wp), allocatable, dimension(:, :, :) :: w_interp4, w_grad4, w_interp2, w_grad2
    !$acc declare create(w_interp4, w_grad4, w_interp2, w_grad2)

    !> Transparent (infinite-domain) diffusion BC state for one boundary face.
    !! The last two indices are the tangential cell indices of the face; U and
    !! f_s carry one ghost layer (MPI seams, corners) for the wall Laplacian.
    type :: t_dif_tbc
        logical :: active = .false.
        real(wp), allocatable, dimension(:, :, :) :: U    !< History per exponential, (1:n_soe, -1:n1+1, -1:n2+1)
        real(wp), allocatable, dimension(:, :, :) :: f_s  !< f = X1 - X1_inf at the current stage, (1, -1:n1+1, -1:n2+1)
        real(wp), allocatable, dimension(:, :) :: f_a    !< f at t_{N-1}
        real(wp), allocatable, dimension(:, :) :: f_b    !< f at t_N
        real(wp), allocatable, dimension(:, :) :: G_a    !< u_t.grad_t f - D Lap_t f at t_{N-1}
        real(wp), allocatable, dimension(:, :) :: G_b    !< u_t.grad_t f - D Lap_t f at t_N
        real(wp), allocatable, dimension(:, :) :: X_inf  !< Far-field mole fraction of gas species 1 (< 0: not set yet)
        real(wp), allocatable, dimension(:, :) :: V_out  !< Outward volume flux of the boundary diffusive flux (last call)
    end type t_dif_tbc

    type(t_dif_tbc), dimension(1:3, 1:2) :: dif_tbc     !< (direction, 1 = beg / 2 = end)
    integer :: n_soe                                    !< Number of exponentials approximating t^(-1/2)
    real(wp), allocatable, dimension(:) :: soe_w, soe_s !< Weights and rates: t^(-1/2) ~ sum w exp(-s t)
    real(wp) :: tbc_d1, tbc_d2   !< Lengths of [t_{N-1}, t_N] and [t_{N-2}, t_{N-1}]
    real(wp) :: tbc_dt_last      !< dt of the previous time step
    real(wp) :: tbc_tau          !< Time offset of the current RK stage from t_N
    integer :: tbc_t_step = -huge(1), tbc_stage = 0
    integer :: tbc_first_dir = 0 !< First direction with a transparent face (stage detection)
    integer, parameter :: dif_tbc_restart_tag = 20261011 !< Format tag of the TBC restart file

    !> Test switch: one-sided-derivative diffusive flux at extrapolation/CBC
    !! boundaries instead of the transparent BC (order 1 or 2)
    logical, parameter :: dif_bc_one_sided = .false.
    integer, parameter :: dif_bc_os_order = 2

    !> Fixed wall composition at extrapolation/CBC boundaries (classic Stefan
    !! tube): the mole fraction of gas species 1 is held at dif_bc_X1_wall on the
    !! boundary face, imposed through the diffusive face flux only (ghost cells
    !! and inviscid fluxes untouched). Takes precedence over the transparent BC;
    !! dif_bc_os_order sets the order of the face gradient (1 or 2).
    logical, parameter :: dif_bc_fixed = .true.
    real(wp), parameter :: dif_bc_X1_wall = 0._wp

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
        !$acc loop seq
        do i = 1, Dif_size
            Ws(i) = fluid_pp(Dif_idx(i))%W
            cps(i) = fluid_pp(Dif_idx(i))%cp
        end do
        !$acc update device(Ws, cps)

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

        call s_initialize_diffusion_tbc()
        call s_read_diffusion_tbc_restart()

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

    subroutine s_compute_diffusion_rhs(idir, j_src_n, rhs_vf, q_prim_vf, irx, iry, irz, t_step)

        integer, intent(in) :: idir
        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf
        type(scalar_field), dimension(sys_size), intent(inout) :: j_src_n, rhs_vf
        type(int_bounds_info), intent(in) :: irx, iry, irz
        integer, intent(in) :: t_step

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
                        h_f(i) = qvs(Dif_idx(i)) + cps(i)*T_f
                    end do

                    ! Compute diffusion fluxes
                    call s_compute_face_species_flux(rhog_f, Y_f, dY_ds_f, j_flux)

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

        ! Diffusive flux at extrapolation/CBC boundaries: either the one-sided
        ! derivative test (dif_bc_one_sided) or the transparent (infinite-domain) BC
        if (dif_bc_one_sided) then
            call s_apply_diffusion_one_sided(idir, j_src_n, q_prim_vf)
        elseif (dif_bc_fixed) then
            call s_apply_diffusion_fixed(idir, j_src_n, q_prim_vf)
        elseif (any(dif_tbc(idir, :)%active)) then
            call s_apply_diffusion_tbc(idir, t_step, j_src_n, q_prim_vf)
        end if

    end subroutine s_compute_diffusion_rhs

    !> Species diffusion fluxes at a face from the face gas density, mass
        !! fractions and mass-fraction gradients: Fick (2 species), closed-form
        !! Stefan-Maxwell (3 species) or the general Stefan-Maxwell solve.
    subroutine s_compute_face_species_flux(rhog_f, Y_f, dY_ds_f, j_flux)

        real(wp), intent(in) :: rhog_f
        real(wp), dimension(Dif_size), intent(in) :: Y_f, dY_ds_f
        real(wp), dimension(Dif_size), intent(out) :: j_flux

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

    end subroutine s_compute_face_species_flux

    !> Test BC: diffusive flux on the physical extrapolation/CBC boundary faces
        !! (bc -3, -5..-13) of direction idir from a one-sided derivative of the
        !! mass fractions taken towards the domain, instead of the ghost cells.
        !!   order 1: (Y1 - Y0)/(s1 - s0), i.e. the first interior face gradient
        !!   order 2: derivative at the wall face of the quadratic through the
        !!            first three cell centres (-2, 3, -1)/ds on a uniform grid
        !! Face density, composition and pressure are the boundary cell values and
        !! the gate is alpha_g**2 (two equal neighbours). MPI seams, periodic,
        !! symmetry, walls and the axis are not touched.
    subroutine s_apply_diffusion_one_sided(idir, j_src_n, q_prim_vf)

        integer, intent(in) :: idir
        type(scalar_field), dimension(sys_size), intent(inout) :: j_src_n
        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf

        integer :: iside, i, j, t1, t2, n1, n2, nn
        integer, dimension(2) :: td
        integer, dimension(0:2) :: kk, ll, qq !< The three cells used, starting at the boundary cell
        integer :: kf, lf, qf                 !< Boundary face
        real(wp) :: s_f, alpha_g, rho_c, rhog, P_c, W_c, T_c, g_c, R_univ
        real(wp), dimension(0:2) :: s, w
        real(wp), dimension(Dif_size) :: Y_c, dY_ds, h_c, j_flux
        real(wp), dimension(0:2, Dif_size) :: Ys

        R_univ = 8314.462618_wp
        call s_tbc_face_dims(idir, td, n1, n2)

        select case (idir)
        case (1)
            nn = m
        case (2)
            nn = n
        case default
            nn = p
        end select
        if (nn < 2) return

        do iside = 1, 2
            if (f_dif_tbc_bc(idir, iside) < -13 .or. f_dif_tbc_bc(idir, iside) > -3) cycle

            ! Wall-face and cell-centre coordinates along idir, inward from the wall
            do j = 0, 2
                select case (idir)
                case (1)
                    s(j) = x_cc(merge(j, m - j, iside == 1))
                case (2)
                    s(j) = y_cc(merge(j, n - j, iside == 1))
                case default
                    s(j) = z_cc(merge(j, p - j, iside == 1))
                end select
            end do
            select case (idir)
            case (1)
                s_f = x_cb(merge(-1, m, iside == 1))
            case (2)
                s_f = y_cb(merge(-1, n, iside == 1))
            case default
                s_f = z_cb(merge(-1, p, iside == 1))
            end select

            ! Derivative weights (d/ds, s increasing along +idir) at the wall face
            w = 0._wp
            if (dif_bc_os_order == 1) then
                w(0) = -1._wp/(s(1) - s(0))
                w(1) = 1._wp/(s(1) - s(0))
            else
                w(0) = ((s_f - s(1)) + (s_f - s(2)))/((s(0) - s(1))*(s(0) - s(2)))
                w(1) = ((s_f - s(0)) + (s_f - s(2)))/((s(1) - s(0))*(s(1) - s(2)))
                w(2) = ((s_f - s(0)) + (s_f - s(1)))/((s(2) - s(0))*(s(2) - s(1)))
            end if

            do t2 = 0, n2
                do t1 = 0, n1
                    ! Cells 0, 1, 2 inward from the wall and the wall face
                    do j = 0, 2
                        select case (idir)
                        case (1)
                            kk(j) = merge(j, m - j, iside == 1); ll(j) = t1; qq(j) = t2
                        case (2)
                            kk(j) = t1; ll(j) = merge(j, n - j, iside == 1); qq(j) = t2
                        case default
                            kk(j) = t1; ll(j) = t2; qq(j) = merge(j, p - j, iside == 1)
                        end select
                    end do
                    kf = kk(0); lf = ll(0); qf = qq(0)
                    select case (idir)
                    case (1)
                        kf = merge(-1, m, iside == 1)
                    case (2)
                        lf = merge(-1, n, iside == 1)
                    case default
                        qf = merge(-1, p, iside == 1)
                    end select

                    do i = 1, Dif_size
                        j_src_n(Dif_idx(i))%sf(kf, lf, qf) = 0._wp
                    end do
                    j_src_n(E_idx)%sf(kf, lf, qf) = 0._wp

                    ! All three cells must contain gas
                    if (any((/(q_prim_vf(advg_idx)%sf(kk(j), ll(j), qq(j)), j=0, 2)/) < small_num_dif)) cycle

                    do j = 0, 2
                        rho_c = 0._wp
                        do i = 1, Dif_size
                            rho_c = rho_c + q_prim_vf(Dif_idx(i))%sf(kk(j), ll(j), qq(j))
                        end do
                        do i = 1, Dif_size
                            Ys(j, i) = q_prim_vf(Dif_idx(i))%sf(kk(j), ll(j), qq(j))/rho_c
                        end do
                    end do

                    do i = 1, Dif_size
                        dY_ds(i) = w(0)*Ys(0, i) + w(1)*Ys(1, i) + w(2)*Ys(2, i)
                        Y_c(i) = Ys(0, i)
                    end do

                    ! Face state = boundary cell state
                    alpha_g = q_prim_vf(advg_idx)%sf(kk(0), ll(0), qq(0))
                    rho_c = 0._wp
                    do i = 1, Dif_size
                        rho_c = rho_c + q_prim_vf(Dif_idx(i))%sf(kk(0), ll(0), qq(0))
                    end do
                    rhog = rho_c/alpha_g
                    P_c = q_prim_vf(E_idx)%sf(kk(0), ll(0), qq(0))
                    W_c = 0._wp
                    do i = 1, Dif_size
                        W_c = W_c + Y_c(i)/Ws(i)
                    end do
                    W_c = 1._wp/W_c
                    T_c = P_c*W_c/(rhog*R_univ)
                    do i = 1, Dif_size
                        h_c(i) = qvs(Dif_idx(i)) + cps(i)*T_c
                    end do

                    call s_compute_face_species_flux(rhog, Y_c, dY_ds, j_flux)

                    g_c = alpha_g**2
                    do i = 1, Dif_size
                        j_src_n(Dif_idx(i))%sf(kf, lf, qf) = g_c*j_flux(i)
                        j_src_n(E_idx)%sf(kf, lf, qf) = j_src_n(E_idx)%sf(kf, lf, qf) + g_c*h_c(i)*j_flux(i)
                    end do
                end do
            end do
        end do

    end subroutine s_apply_diffusion_one_sided

    !> Fixed wall composition on the extrapolation/CBC boundary faces of direction
        !! idir (classic Stefan tube): X1 = dif_bc_X1_wall on the face. Written in the
        !! mole fraction, whose binary flux is exactly Fickian with constant C*D, and
        !! converted to the mass flux with grad Y1 = W1 W2/W**2 grad X1:
        !!   j_out = alpha_g**2 rho_g W1 W2/W_b**2 D12 dX1/ds,   s inward from the face
        !! dX1/ds from (X1_wall, X1 of the first one or two cells): order 1 or 2
        !! (dif_bc_os_order). Species 2 carries -j_out and energy (h1 - h2)*j_out, so
        !! the net mass flux is zero and the inviscid/CBC treatment is unaffected.
        !! Binary gas only.
    subroutine s_apply_diffusion_fixed(idir, j_src_n, q_prim_vf)

        integer, intent(in) :: idir
        type(scalar_field), dimension(sys_size), intent(inout) :: j_src_n
        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf

        integer :: iside, i, t1, t2, n1, n2, nn, bc
        integer, dimension(2) :: td
        integer, dimension(3) :: e
        integer :: kc, lc, qc, kf, lf, qf, k1, l1, q1
        logical :: use_second
        real(wp) :: R_univ, sgn, D12, s_f, s0, s1, h0, h1, dXds
        real(wp) :: alpha_g, rho_c, rhog, Y1, X0, X1n, W_b, T_b, J_out, rho_1, alpha_1
        real(wp), dimension(2) :: h_e

        if (Dif_size /= 2) return

        R_univ = 8314.462618_wp
        D12 = Ds(1, 2)
        call s_tbc_face_dims(idir, td, n1, n2)
        e = 0; e(idir) = 1

        select case (idir)
        case (1)
            nn = m
        case (2)
            nn = n
        case default
            nn = p
        end select

        do iside = 1, 2
            bc = f_dif_tbc_bc(idir, iside)
            if (bc < -13 .or. bc > -3) cycle

            ! j_src_n is positive along +s; the outward normal is -s at the beginning
            sgn = merge(-1._wp, 1._wp, iside == 1)

            do t2 = 0, n2
                do t1 = 0, n1
                    call s_tbc_cell(idir, iside, t1, t2, kc, lc, qc, kf, lf, qf)

                    alpha_g = q_prim_vf(advg_idx)%sf(kc, lc, qc)
                    rho_c = q_prim_vf(Dif_idx(1))%sf(kc, lc, qc) + q_prim_vf(Dif_idx(2))%sf(kc, lc, qc)

                    if (.not. (alpha_g > small_num_dif .and. rho_c > 0._wp)) then
                        j_src_n(Dif_idx(1))%sf(kf, lf, qf) = 0._wp
                        j_src_n(Dif_idx(2))%sf(kf, lf, qf) = 0._wp
                        j_src_n(E_idx)%sf(kf, lf, qf) = 0._wp
                        cycle
                    end if

                    rhog = rho_c/alpha_g
                    Y1 = q_prim_vf(Dif_idx(1))%sf(kc, lc, qc)/rho_c
                    W_b = 1._wp/(Y1/Ws(1) + (1._wp - Y1)/Ws(2))
                    X0 = Y1*W_b/Ws(1)

                    ! Second cell inward (for the order-2 face gradient)
                    k1 = kc + nint(-sgn)*e(1); l1 = lc + nint(-sgn)*e(2); q1 = qc + nint(-sgn)*e(3)
                    use_second = dif_bc_os_order >= 2 .and. nn >= 1
                    if (use_second) then
                        alpha_1 = q_prim_vf(advg_idx)%sf(k1, l1, q1)
                        rho_1 = q_prim_vf(Dif_idx(1))%sf(k1, l1, q1) + q_prim_vf(Dif_idx(2))%sf(k1, l1, q1)
                        use_second = alpha_1 > small_num_dif .and. rho_1 > 0._wp
                    end if

                    ! Distances from the face to the first two cell centres
                    select case (idir)
                    case (1)
                        s_f = x_cb(kf); s0 = x_cc(kc)
                        if (use_second) s1 = x_cc(k1)
                    case (2)
                        s_f = y_cb(lf); s0 = y_cc(lc)
                        if (use_second) s1 = y_cc(l1)
                    case default
                        s_f = z_cb(qf); s0 = z_cc(qc)
                        if (use_second) s1 = z_cc(q1)
                    end select
                    h0 = abs(s0 - s_f)

                    if (use_second) then
                        h1 = abs(s1 - s_f)
                        X1n = q_prim_vf(Dif_idx(1))%sf(k1, l1, q1)/rho_1
                        X1n = X1n/(X1n/Ws(1) + (1._wp - X1n)/Ws(2))/Ws(1)
                        ! Derivative at the face of the quadratic through (0, X_wall), (h0, X0), (h1, X1n)
                        dXds = -(1._wp/h0 + 1._wp/h1)*dif_bc_X1_wall &
                               + h1/(h0*(h1 - h0))*X0 &
                               - h0/(h1*(h1 - h0))*X1n
                    else
                        dXds = (X0 - dif_bc_X1_wall)/h0
                    end if

                    ! Outward mass flux of species 1, gated like the interior faces
                    J_out = alpha_g**2*rhog*Ws(1)*Ws(2)/W_b**2*D12*dXds

                    T_b = q_prim_vf(E_idx)%sf(kc, lc, qc)*W_b/(rhog*R_univ)
                    do i = 1, 2
                        h_e(i) = qvs(Dif_idx(i)) + cps(i)*T_b
                    end do

                    j_src_n(Dif_idx(1))%sf(kf, lf, qf) = sgn*J_out
                    j_src_n(Dif_idx(2))%sf(kf, lf, qf) = -sgn*J_out
                    j_src_n(E_idx)%sf(kf, lf, qf) = sgn*(h_e(1) - h_e(2))*J_out
                end do
            end do
        end do

    end subroutine s_apply_diffusion_fixed

    !> Sets up the transparent (infinite-domain) diffusion BC on every physical
        !! boundary that uses ghost-cell extrapolation (-3) or a CBC (-5 to -13).
        !! MPI interior edges, periodic, symmetry, walls and the axis are excluded
        !! because their bc values fall outside that range.
        !! The kernel t^(-1/2) of the half-order time derivative is approximated
        !! on [delta, T] by a sum of exponentials, from the trapezoidal rule on
        !!   t^(-1/2) = (1/sqrt(pi)) * int_{-inf}^{inf} exp(-t e^u) e^(u/2) du
    subroutine s_initialize_diffusion_tbc

        integer :: idir, iside, bc, i, n1, n2
        integer, dimension(2) :: td
        logical :: any_tbc
        real(wp) :: delta, t_span, u_min, u_max, h_u, tol

        ! Not used when another wall treatment is selected
        if (dif_bc_one_sided .or. dif_bc_fixed) return

        any_tbc = .false.
        do idir = 1, num_dims
            do iside = 1, 2
                bc = f_dif_tbc_bc(idir, iside)
                if (bc >= -13 .and. bc <= -3) any_tbc = .true.
            end do
        end do
        if (.not. any_tbc) return

        if (Dif_size /= 2) then
            if (proc_rank == 0) print *, 'WARNING: transparent diffusion BC requires exactly 2 gas species; ', &
                'extrapolation/CBC boundaries keep zero diffusive flux.'
            return
        end if

        ! Kernel range: smallest lag is a fraction of a time step, largest twice the
        ! absolute end time (the history starts at t = 0, also after a restart)
        tol = 1.e-9_wp
        h_u = 0.4_wp
        if (cfl_dt) then
            t_span = t_stop
            delta = 1.e-9_wp*t_stop
        else
            t_span = real(t_step_stop, wp)*dt
            delta = 0.1_wp*dt
        end if
        t_span = max(2._wp*t_span, 100._wp*delta)
        u_max = log(40._wp/delta)
        u_min = 2._wp*log(tol*sqrt(pi)/(2._wp*sqrt(t_span)))
        n_soe = ceiling((u_max - u_min)/h_u) + 1

        allocate (soe_w(1:n_soe), soe_s(1:n_soe))
        do i = 1, n_soe
            soe_s(i) = exp(u_min + real(i - 1, wp)*h_u)
            soe_w(i) = h_u*sqrt(soe_s(i))/sqrt(pi)
        end do

        do idir = 1, num_dims
            call s_tbc_face_dims(idir, td, n1, n2)

            do iside = 1, 2
                bc = f_dif_tbc_bc(idir, iside)
                if (bc < -13 .or. bc > -3) cycle

                dif_tbc(idir, iside)%active = .true.
                if (tbc_first_dir == 0) tbc_first_dir = idir

                ! Arrays that need tangential neighbours carry one ghost layer
                allocate (dif_tbc(idir, iside)%U(1:n_soe, -1:n1 + 1, -1:n2 + 1))
                allocate (dif_tbc(idir, iside)%f_s(1:1, -1:n1 + 1, -1:n2 + 1))
                allocate (dif_tbc(idir, iside)%f_a(0:n1, 0:n2))
                allocate (dif_tbc(idir, iside)%f_b(0:n1, 0:n2))
                allocate (dif_tbc(idir, iside)%G_a(0:n1, 0:n2))
                allocate (dif_tbc(idir, iside)%G_b(0:n1, 0:n2))
                allocate (dif_tbc(idir, iside)%X_inf(0:n1, 0:n2))
                allocate (dif_tbc(idir, iside)%V_out(0:n1, 0:n2))
                dif_tbc(idir, iside)%U = 0._wp
                dif_tbc(idir, iside)%f_s = 0._wp
                dif_tbc(idir, iside)%f_a = 0._wp
                dif_tbc(idir, iside)%f_b = 0._wp
                dif_tbc(idir, iside)%G_a = 0._wp
                dif_tbc(idir, iside)%G_b = 0._wp
                dif_tbc(idir, iside)%X_inf = -1._wp
                dif_tbc(idir, iside)%V_out = 0._wp
            end do
        end do

    end subroutine s_initialize_diffusion_tbc

    !> Restart file of the transparent diffusion BC state of this rank:
        !! p_all/p<rank>/<step>/dif_tbc.dat (serial I/O) or
        !! restart_data/dif_tbc_<step>_<rank>.dat (parallel I/O).
    function f_dif_tbc_restart_path(step) result(path)

        integer, intent(in) :: step
        character(len=path_len) :: path

        if (parallel_io) then
            write (path, '(A,I0,A,I0,A)') trim(case_dir)//'/restart_data/dif_tbc_', step, '_', proc_rank, '.dat'
        else
            write (path, '(A,I0,A,I0,A)') trim(case_dir)//'/p_all/p', proc_rank, '/', step, '/dif_tbc.dat'
        end if

    end function f_dif_tbc_restart_path

    !> Writes the transparent diffusion BC state (kernel, step lengths and, per
        !! active face, the histories U, f and G at t_{N-1} and t_N, X_inf and V_out)
        !! next to the flow restart data of save index step. After a completed time
        !! step U covers [0, t_{N-1}], so this is all the next step needs.
        !! Must be called after the flow data of the same step are written (the
        !! serial writer recreates the step directory).
    subroutine s_write_diffusion_tbc_restart(step)

        integer, intent(in) :: step

        integer :: idir, iside, iu, act
        logical :: any_active

        if (dif_bc_one_sided) return
        any_active = .false.
        do idir = 1, 3
            do iside = 1, 2
                if (dif_tbc(idir, iside)%active) any_active = .true.
            end do
        end do
        if (.not. any_active) return

        open (newunit=iu, file=trim(f_dif_tbc_restart_path(step)), form='unformatted', &
              access='stream', status='replace', action='write')

        write (iu) dif_tbc_restart_tag, n_soe
        write (iu) soe_s(1:n_soe), soe_w(1:n_soe)
        write (iu) tbc_d1, tbc_dt_last

        do idir = 1, 3
            do iside = 1, 2
                act = merge(1, 0, dif_tbc(idir, iside)%active)
                write (iu) act
                if (act == 0) cycle
                write (iu) size(dif_tbc(idir, iside)%f_a, 1), size(dif_tbc(idir, iside)%f_a, 2)
                write (iu) dif_tbc(idir, iside)%U, &
                    dif_tbc(idir, iside)%f_a, dif_tbc(idir, iside)%f_b, &
                    dif_tbc(idir, iside)%G_a, dif_tbc(idir, iside)%G_b, &
                    dif_tbc(idir, iside)%X_inf, dif_tbc(idir, iside)%V_out
            end do
        end do

        close (iu)

    end subroutine s_write_diffusion_tbc_restart

    !> Restores the transparent diffusion BC state on a restart. The saved kernel
        !! replaces the one built in s_initialize_diffusion_tbc so the histories stay
        !! consistent with their rates. Without a file the BC starts cold (warning).
        !! The MPI decomposition must match the run that wrote the file.
    subroutine s_read_diffusion_tbc_restart

        integer :: step, idir, iside, iu, act, tag, n_soe_f, n1f, n2f, ios
        integer :: n1, n2
        integer, dimension(2) :: td
        logical :: any_active, file_exist
        character(len=path_len) :: path

        if (dif_bc_one_sided) return
        any_active = .false.
        do idir = 1, 3
            do iside = 1, 2
                if (dif_tbc(idir, iside)%active) any_active = .true.
            end do
        end do
        if (.not. any_active) return

        step = merge(n_start, t_step_start, cfl_dt)
        if (step <= 0) return

        path = f_dif_tbc_restart_path(step)
        inquire (file=trim(path), exist=file_exist)
        if (.not. file_exist) then
            print '(A,I0,A)', 'WARNING: rank ', proc_rank, ': no transparent diffusion BC restart file '// &
                trim(path)//'; the BC history starts from zero.'
            return
        end if

        open (newunit=iu, file=trim(path), form='unformatted', access='stream', &
              status='old', action='read', iostat=ios)
        if (ios /= 0) call s_mpi_abort('Cannot open '//trim(path))

        read (iu) tag, n_soe_f
        if (tag /= dif_tbc_restart_tag) call s_mpi_abort('Unknown format of '//trim(path))

        ! Use the kernel the histories were built with
        if (n_soe_f /= n_soe) then
            deallocate (soe_s, soe_w)
            allocate (soe_s(1:n_soe_f), soe_w(1:n_soe_f))
            n_soe = n_soe_f
            do idir = 1, num_dims
                do iside = 1, 2
                    if (.not. dif_tbc(idir, iside)%active) cycle
                    call s_tbc_face_dims(idir, td, n1, n2)
                    deallocate (dif_tbc(idir, iside)%U)
                    allocate (dif_tbc(idir, iside)%U(1:n_soe, -1:n1 + 1, -1:n2 + 1))
                end do
            end do
        end if
        read (iu) soe_s(1:n_soe), soe_w(1:n_soe)
        read (iu) tbc_d1, tbc_dt_last

        do idir = 1, 3
            do iside = 1, 2
                read (iu) act
                if ((act == 1) .neqv. dif_tbc(idir, iside)%active) then
                    call s_mpi_abort('Transparent BC faces in '//trim(path)// &
                                     ' do not match this run (same MPI decomposition required)')
                end if
                if (act == 0) cycle
                read (iu) n1f, n2f
                if (n1f /= size(dif_tbc(idir, iside)%f_a, 1) .or. n2f /= size(dif_tbc(idir, iside)%f_a, 2)) then
                    call s_mpi_abort('Transparent BC face size in '//trim(path)// &
                                     ' does not match this run (same MPI decomposition required)')
                end if
                read (iu) dif_tbc(idir, iside)%U, &
                    dif_tbc(idir, iside)%f_a, dif_tbc(idir, iside)%f_b, &
                    dif_tbc(idir, iside)%G_a, dif_tbc(idir, iside)%G_b, &
                    dif_tbc(idir, iside)%X_inf, dif_tbc(idir, iside)%V_out
            end do
        end do

        close (iu)

        ! Next call is stage 1 of a new step that continues the restored history
        ! (and not the cold-start branch, which would reset the step lengths)
        tbc_t_step = -huge(1) + 1

        if (proc_rank == 0) print '(A,I0)', 'Restored transparent diffusion BC history from step ', step

    end subroutine s_read_diffusion_tbc_restart

    !> Boundary condition code of face (idir, iside), iside = 1 beginning / 2 end.
        !! On MPI interior edges (and periodic ones) this is the neighbour's rank (>= 0).
    function f_dif_tbc_bc(idir, iside) result(bc)

        integer, intent(in) :: idir, iside
        integer :: bc

        select case (idir)
        case (1)
            bc = merge(bc_x%beg, bc_x%end, iside == 1)
        case (2)
            bc = merge(bc_y%beg, bc_y%end, iside == 1)
        case default
            bc = merge(bc_z%beg, bc_z%end, iside == 1)
        end select

    end function f_dif_tbc_bc

    !> Tangential directions td of a face normal to idir and the local index
        !! ranges 0:n1, 0:n2 of the face cells along them.
    subroutine s_tbc_face_dims(idir, td, n1, n2)

        integer, intent(in) :: idir
        integer, dimension(2), intent(out) :: td
        integer, intent(out) :: n1, n2

        select case (idir)
        case (1)
            td = (/2, 3/); n1 = n; n2 = p
        case (2)
            td = (/1, 3/); n1 = m; n2 = p
        case default
            td = (/1, 2/); n1 = m; n2 = n
        end select

    end subroutine s_tbc_face_dims

    !> Cell centres c and widths w along direction d, including one ghost cell.
    subroutine s_tbc_coords(d, nn, c, w)

        integer, intent(in) :: d, nn
        real(wp), dimension(-1:nn + 1), intent(out) :: c, w

        select case (d)
        case (1)
            c = x_cc(-1:nn + 1); w = dx(-1:nn + 1)
        case (2)
            c = y_cc(-1:nn + 1); w = dy(-1:nn + 1)
        case default
            c = z_cc(-1:nn + 1); w = dz(-1:nn + 1)
        end select

    end subroutine s_tbc_coords

    !> Fills the tangential ghost layer of a face array A(1:nv, -1:n1+1, -1:n2+1).
        !! Along each tangential direction, per side:
        !!   neighbour rank (MPI seam or periodic) -> exchanged with that rank
        !!   periodic on a single rank (-1)        -> wrapped
        !!   another transparent boundary (corner) -> linear extrapolation
        !!   anything else (wall, symmetry, axis)  -> mirror (zero gradient)
        !! All ranks along a wall own the same face, so the exchanges pair up.
    subroutine s_tbc_fill_ghosts(A, nv, n1, n2, td, tag)

        integer, intent(in) :: nv, n1, n2, tag
        integer, dimension(2), intent(in) :: td
        real(wp), dimension(1:nv, -1:n1 + 1, -1:n2 + 1), intent(inout) :: A

        integer :: ax, nn, nt, side, bc, j, dest, src
        integer, dimension(2) :: bcs
        real(wp), allocatable, dimension(:) :: c, w
        real(wp), allocatable, dimension(:, :) :: sbuf, rbuf
#ifdef MFC_MPI
        integer :: ierr
#endif

        do ax = 1, 2
            if (td(ax) > num_dims) cycle

            if (ax == 1) then
                nn = n1; nt = n2
            else
                nn = n2; nt = n1
            end if

            allocate (c(-1:nn + 1), w(-1:nn + 1))
            call s_tbc_coords(td(ax), nn, c, w)
            bcs(1) = f_dif_tbc_bc(td(ax), 1)
            bcs(2) = f_dif_tbc_bc(td(ax), 2)

#ifdef MFC_MPI
            if (any(bcs >= 0)) then
                allocate (sbuf(1:nv, 0:nt), rbuf(1:nv, 0:nt))
                do side = 1, 2
                    ! side 1: send first row to beg, receive end ghost from end
                    ! side 2: send last row to end, receive beg ghost from beg
                    dest = merge(bcs(1), bcs(2), side == 1)
                    src = merge(bcs(2), bcs(1), side == 1)
                    if (dest < 0) dest = MPI_PROC_NULL
                    if (src < 0) src = MPI_PROC_NULL

                    if (ax == 1) then
                        sbuf = A(:, merge(0, nn, side == 1), 0:nt)
                    else
                        sbuf = A(:, 0:nt, merge(0, nn, side == 1))
                    end if

                    call MPI_SENDRECV(sbuf, nv*(nt + 1), mpi_p, dest, tag + 2*(ax - 1) + side, &
                                      rbuf, nv*(nt + 1), mpi_p, src, tag + 2*(ax - 1) + side, &
                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

                    if (src /= MPI_PROC_NULL) then
                        if (ax == 1) then
                            A(:, merge(nn + 1, -1, side == 1), 0:nt) = rbuf
                        else
                            A(:, 0:nt, merge(nn + 1, -1, side == 1)) = rbuf
                        end if
                    end if
                end do
                deallocate (sbuf, rbuf)
            end if
#endif

            do side = 1, 2
                bc = bcs(side)
                if (bc >= 0) cycle

                do j = 0, nt
                    if (bc == -1) then
                        call s_tbc_set_ghost(A, nv, n1, n2, ax, j, side, &
                                             f_tbc_get(A, nv, n1, n2, ax, j, merge(nn, 0, side == 1)))
                    elseif (bc >= -13 .and. bc <= -3 .and. nn >= 1) then
                        if (side == 1) then
                            call s_tbc_set_ghost(A, nv, n1, n2, ax, j, side, &
                                                 f_tbc_get(A, nv, n1, n2, ax, j, 0) + &
                                                 (f_tbc_get(A, nv, n1, n2, ax, j, 0) - f_tbc_get(A, nv, n1, n2, ax, j, 1))* &
                                                 (c(0) - c(-1))/(c(1) - c(0)))
                        else
                            call s_tbc_set_ghost(A, nv, n1, n2, ax, j, side, &
                                                 f_tbc_get(A, nv, n1, n2, ax, j, nn) + &
                                                 (f_tbc_get(A, nv, n1, n2, ax, j, nn) - f_tbc_get(A, nv, n1, n2, ax, j, nn - 1))* &
                                                 (c(nn + 1) - c(nn))/(c(nn) - c(nn - 1)))
                        end if
                    else
                        call s_tbc_set_ghost(A, nv, n1, n2, ax, j, side, &
                                             f_tbc_get(A, nv, n1, n2, ax, j, merge(0, nn, side == 1)))
                    end if
                end do
            end do

            deallocate (c, w)
        end do

    end subroutine s_tbc_fill_ghosts

    !> Values A(:, i, j) along tangential axis 1, or A(:, j, i) along axis 2.
    function f_tbc_get(A, nv, n1, n2, ax, j, i) result(v)

        integer, intent(in) :: nv, n1, n2, ax, j, i
        real(wp), dimension(1:nv, -1:n1 + 1, -1:n2 + 1), intent(in) :: A
        real(wp), dimension(1:nv) :: v

        if (ax == 1) then
            v = A(:, i, j)
        else
            v = A(:, j, i)
        end if

    end function f_tbc_get

    !> Sets the ghost value at side (1 = beginning, 2 = end) of tangential axis ax.
    subroutine s_tbc_set_ghost(A, nv, n1, n2, ax, j, side, v)

        integer, intent(in) :: nv, n1, n2, ax, j, side
        real(wp), dimension(1:nv, -1:n1 + 1, -1:n2 + 1), intent(inout) :: A
        real(wp), dimension(1:nv), intent(in) :: v

        if (ax == 1) then
            A(:, merge(-1, n1 + 1, side == 1), j) = v
        else
            A(:, j, merge(-1, n2 + 1, side == 1)) = v
        end if

    end subroutine s_tbc_set_ghost

    !> Tangential Laplacian of A (ghosts filled) at the face cells: lap(1:nv, 0:n1, 0:n2).
    subroutine s_tbc_lap(A, nv, n1, n2, td, lap)

        integer, intent(in) :: nv, n1, n2
        integer, dimension(2), intent(in) :: td
        real(wp), dimension(1:nv, -1:n1 + 1, -1:n2 + 1), intent(in) :: A
        real(wp), dimension(1:nv, 0:n1, 0:n2), intent(out) :: lap

        integer :: t1, t2
        real(wp), dimension(-1:n1 + 1) :: c1, w1
        real(wp), dimension(-1:n2 + 1) :: c2, w2

        lap = 0._wp

        if (td(1) <= num_dims) then
            call s_tbc_coords(td(1), n1, c1, w1)
            do t2 = 0, n2
                do t1 = 0, n1
                    lap(:, t1, t2) = lap(:, t1, t2) + &
                                     ((A(:, t1 + 1, t2) - A(:, t1, t2))/(c1(t1 + 1) - c1(t1)) - &
                                      (A(:, t1, t2) - A(:, t1 - 1, t2))/(c1(t1) - c1(t1 - 1)))/w1(t1)
                end do
            end do
        end if

        if (td(2) <= num_dims) then
            call s_tbc_coords(td(2), n2, c2, w2)
            do t2 = 0, n2
                do t1 = 0, n1
                    lap(:, t1, t2) = lap(:, t1, t2) + &
                                     ((A(:, t1, t2 + 1) - A(:, t1, t2))/(c2(t2 + 1) - c2(t2)) - &
                                      (A(:, t1, t2) - A(:, t1, t2 - 1))/(c2(t2) - c2(t2 - 1)))/w2(t2)
                end do
            end do
        end if

    end subroutine s_tbc_lap

    !> Tangential advection sum_ax ut(ax) * dA/dtau_ax at the face cells, with
        !! central (upwind = .false.) or first-order upwind differences. Ghosts filled.
    subroutine s_tbc_tangent_adv(A, nv, n1, n2, td, ut, upwind, adv)

        integer, intent(in) :: nv, n1, n2
        integer, dimension(2), intent(in) :: td
        real(wp), dimension(1:nv, -1:n1 + 1, -1:n2 + 1), intent(in) :: A
        real(wp), dimension(1:2, 0:n1, 0:n2), intent(in) :: ut
        logical, intent(in) :: upwind
        real(wp), dimension(1:nv, 0:n1, 0:n2), intent(out) :: adv

        integer :: t1, t2
        real(wp), dimension(-1:n1 + 1) :: c1, w1
        real(wp), dimension(-1:n2 + 1) :: c2, w2

        adv = 0._wp

        if (td(1) <= num_dims) then
            call s_tbc_coords(td(1), n1, c1, w1)
            do t2 = 0, n2
                do t1 = 0, n1
                    if (.not. upwind) then
                        adv(:, t1, t2) = adv(:, t1, t2) + ut(1, t1, t2)* &
                                         (A(:, t1 + 1, t2) - A(:, t1 - 1, t2))/(c1(t1 + 1) - c1(t1 - 1))
                    elseif (ut(1, t1, t2) >= 0._wp) then
                        adv(:, t1, t2) = adv(:, t1, t2) + ut(1, t1, t2)* &
                                         (A(:, t1, t2) - A(:, t1 - 1, t2))/(c1(t1) - c1(t1 - 1))
                    else
                        adv(:, t1, t2) = adv(:, t1, t2) + ut(1, t1, t2)* &
                                         (A(:, t1 + 1, t2) - A(:, t1, t2))/(c1(t1 + 1) - c1(t1))
                    end if
                end do
            end do
        end if

        if (td(2) <= num_dims) then
            call s_tbc_coords(td(2), n2, c2, w2)
            do t2 = 0, n2
                do t1 = 0, n1
                    if (.not. upwind) then
                        adv(:, t1, t2) = adv(:, t1, t2) + ut(2, t1, t2)* &
                                         (A(:, t1, t2 + 1) - A(:, t1, t2 - 1))/(c2(t2 + 1) - c2(t2 - 1))
                    elseif (ut(2, t1, t2) >= 0._wp) then
                        adv(:, t1, t2) = adv(:, t1, t2) + ut(2, t1, t2)* &
                                         (A(:, t1, t2) - A(:, t1, t2 - 1))/(c2(t2) - c2(t2 - 1))
                    else
                        adv(:, t1, t2) = adv(:, t1, t2) + ut(2, t1, t2)* &
                                         (A(:, t1, t2 + 1) - A(:, t1, t2))/(c2(t2 + 1) - c2(t2))
                    end if
                end do
            end do
        end if

    end subroutine s_tbc_tangent_adv

    !> Moments of the shifted kernel over the local interval [0, L]:
        !!   M0 = int_0^L exp(-c r) r^(-1/2) dr,  M1 = int_0^L exp(-c r) r^(1/2) dr
    subroutine s_tbc_moments(c, L, M0, M1)

        real(wp), intent(in) :: c, L
        real(wp), intent(out) :: M0, M1

        real(wp) :: cL

        cL = c*L
        if (cL < 1.e-3_wp) then
            M0 = 2._wp*sqrt(L)*(1._wp - cL/3._wp + cL*cL/10._wp)
            M1 = 2._wp/3._wp*L*sqrt(L)*(1._wp - 0.6_wp*cL + 3._wp/14._wp*cL*cL)
        else
            M0 = sqrt(pi/c)*erf(sqrt(cL))
            M1 = (0.5_wp*M0 - sqrt(L)*exp(-cL))/c
        end if

    end subroutine s_tbc_moments

    !> Overwrites the diffusive fluxes on the transparent boundary faces of
        !! direction idir. The exterior problem is posed in the mole fraction X1:
        !! for a binary ideal gas at constant p and T the molar concentration is
        !! constant and the molar-average velocity u* = u + V is divergence free,
        !! so dX1/dt + u*.grad X1 = D Lap X1 is exactly linear (V = W/rho_g
        !! sum_i j_i/W_i is the volume flux of the diffusive fluxes, taken from the
        !! previous call). Outside the domain u* is a uniform stream with the
        !! boundary cell's value (u*_n along the outward normal, u_t along the wall)
        !! and the composition is uniform, equal to the boundary cell's composition
        !! when gas is first seen there. Solving that exterior half-space exactly
        !! and using grad Y1 = W1 W2/W**2 grad X1 gives the outgoing mass flux
        !!   j_out = rho_g*W1*W2/W**2*[sqrt(D) (d/dt + u_t.grad_t - D Lap_t + c)^(1/2) - u*_n/2] f
        !! with f = X1_b - X1_inf and c = u*_n**2/(4 D). For u* = 0 this is the
        !! quiescent half-space condition; for strong outflow j_out -> 0 and for
        !! strong inflow j_out -> rho_g*|u_n|*f (Danckwerts).
        !! Writing the kernel of (d/dt + ... + c)^(-1/2) as a sum of exponentials
        !! gives one history Z_l per exponential, obeying
        !!   dZ_l/dt = -(s_l + c) Z_l - u_t.grad_t Z_l + D Lap_t Z_l + h,
        !!   h = df/dt + c f + G,  G = u_t.grad_t f - D Lap_t f
        !! and the operator equals (1/sqrt(pi)) sum_l w_l Z_l. The history covers
        !! [0, t_{N-1}]; the part from t_{N-1} to the current RK stage is
        !! integrated exactly for f and G linear in time.
        !! The exterior velocity is frozen per step at the boundary cell value.
        !! Species 2 carries -j_out and energy (h1 - h2)*j_out.
    subroutine s_apply_diffusion_tbc(idir, t_step, j_src_n, q_prim_vf)

        integer, intent(in) :: idir, t_step
        type(scalar_field), dimension(sys_size), intent(inout) :: j_src_n
        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf

        integer :: iside, i, ax, t1, t2, n1, n2
        integer, dimension(2) :: td
        integer :: kc, lc, qc, kf, lf, qf !< Boundary cell and boundary face indices
        logical :: gas_present
        real(wp) :: R_univ, sgn, L, x, e_x, phi_x, psi_x, D12, c, fp, h_a, h_b, h_s
        real(wp) :: hist, M0, M1, op_f, J_out
        real(wp) :: alpha_g, rho_c, rhog, Y1, X1, W_b, T_b
        real(wp), dimension(2) :: h_e
        real(wp), allocatable, dimension(:, :, :) :: G, lapf, advf, lapU, advU, ut
        real(wp), allocatable, dimension(:, :) :: un

        R_univ = 8314.462618_wp
        D12 = Ds(1, 2)

        ! RK stage bookkeeping, done once per stage: the first call of a new
        ! time step is stage 1, whose state is the completed time level t_N
        if (idir == tbc_first_dir) then
            if (t_step /= tbc_t_step) then
                if (tbc_t_step == -huge(1)) then
                    tbc_dt_last = dt; tbc_d1 = dt
                end if
                tbc_t_step = t_step
                tbc_stage = 1
                tbc_d2 = tbc_d1
                tbc_d1 = tbc_dt_last
                tbc_dt_last = dt
            else
                tbc_stage = tbc_stage + 1
            end if

            ! Stage times relative to t_N (SSP-RK2: 0, dt; SSP-RK3: 0, dt, dt/2)
            tbc_tau = 0._wp
            if (time_stepper == 2 .and. tbc_stage == 2) tbc_tau = dt
            if (time_stepper == 3 .and. tbc_stage == 2) tbc_tau = dt
            if (time_stepper == 3 .and. tbc_stage == 3) tbc_tau = 0.5_wp*dt
        end if

        call s_tbc_face_dims(idir, td, n1, n2)
        L = tbc_d1 + tbc_tau
        allocate (G(1:1, 0:n1, 0:n2), lapf(1:1, 0:n1, 0:n2), advf(1:1, 0:n1, 0:n2))
        allocate (ut(1:2, 0:n1, 0:n2), un(0:n1, 0:n2))

        do iside = 1, 2
            if (.not. dif_tbc(idir, iside)%active) cycle

            ! j_src_n is positive along +s; the outward normal is -s at the beginning
            sgn = merge(-1._wp, 1._wp, iside == 1)

            ! Deviation f from the far field at the current stage (held while no
            ! gas) and the exterior velocity, taken from the boundary cell
            do t2 = 0, n2
                do t1 = 0, n1
                    call s_tbc_cell(idir, iside, t1, t2, kc, lc, qc, kf, lf, qf)
                    alpha_g = q_prim_vf(advg_idx)%sf(kc, lc, qc)
                    rho_c = q_prim_vf(Dif_idx(1))%sf(kc, lc, qc) + q_prim_vf(Dif_idx(2))%sf(kc, lc, qc)

                    ! Molar-average (= volume-average) normal velocity u* = u + V
                    un(t1, t2) = sgn*q_prim_vf(momxb + idir - 1)%sf(kc, lc, qc) + dif_tbc(idir, iside)%V_out(t1, t2)
                    do ax = 1, 2
                        ut(ax, t1, t2) = 0._wp
                        if (td(ax) <= num_dims) ut(ax, t1, t2) = q_prim_vf(momxb + td(ax) - 1)%sf(kc, lc, qc)
                    end do

                    dif_tbc(idir, iside)%f_s(1, t1, t2) = dif_tbc(idir, iside)%f_b(t1, t2)
                    if (alpha_g > small_num_dif .and. rho_c > 0._wp) then
                        Y1 = q_prim_vf(Dif_idx(1))%sf(kc, lc, qc)/rho_c
                        X1 = (Y1/Ws(1))/(Y1/Ws(1) + (1._wp - Y1)/Ws(2))
                        if (dif_tbc(idir, iside)%X_inf(t1, t2) < 0._wp) then
                            dif_tbc(idir, iside)%X_inf(t1, t2) = X1
                            dif_tbc(idir, iside)%f_a(t1, t2) = 0._wp
                            dif_tbc(idir, iside)%f_b(t1, t2) = 0._wp
                            dif_tbc(idir, iside)%G_a(t1, t2) = 0._wp
                            dif_tbc(idir, iside)%G_b(t1, t2) = 0._wp
                            dif_tbc(idir, iside)%U(:, t1, t2) = 0._wp
                        end if
                        dif_tbc(idir, iside)%f_s(1, t1, t2) = X1 - dif_tbc(idir, iside)%X_inf(t1, t2)
                    end if
                end do
            end do

            ! G = u_t.grad_t f - D Lap_t f at the current stage
            call s_tbc_fill_ghosts(dif_tbc(idir, iside)%f_s, 1, n1, n2, td, 1000*idir + 100*iside)
            call s_tbc_lap(dif_tbc(idir, iside)%f_s, 1, n1, n2, td, lapf)
            call s_tbc_tangent_adv(dif_tbc(idir, iside)%f_s, 1, n1, n2, td, ut, .false., advf)
            G = advf - D12*lapf

            ! Stage 1: advance the history over [t_{N-2}, t_{N-1}], then shift
            if (tbc_stage == 1) then
                allocate (lapU(1:n_soe, 0:n1, 0:n2), advU(1:n_soe, 0:n1, 0:n2))
                call s_tbc_fill_ghosts(dif_tbc(idir, iside)%U, n_soe, n1, n2, td, 1000*idir + 100*iside + 10)
                call s_tbc_lap(dif_tbc(idir, iside)%U, n_soe, n1, n2, td, lapU)
                call s_tbc_tangent_adv(dif_tbc(idir, iside)%U, n_soe, n1, n2, td, ut, .true., advU)

                do t2 = 0, n2
                    do t1 = 0, n1
                        c = un(t1, t2)**2/(4._wp*D12)
                        fp = (dif_tbc(idir, iside)%f_b(t1, t2) - dif_tbc(idir, iside)%f_a(t1, t2))/tbc_d2
                        h_a = fp + c*dif_tbc(idir, iside)%f_a(t1, t2) + dif_tbc(idir, iside)%G_a(t1, t2)
                        h_b = fp + c*dif_tbc(idir, iside)%f_b(t1, t2) + dif_tbc(idir, iside)%G_b(t1, t2)

                        do i = 1, n_soe
                            x = (soe_s(i) + c)*tbc_d2
                            e_x = 0._wp
                            if (x < 700._wp) e_x = exp(-x)
                            if (x < 1.e-8_wp) then
                                phi_x = 1._wp - 0.5_wp*x
                            else
                                phi_x = (1._wp - e_x)/x
                            end if
                            if (x < 1.e-4_wp) then
                                psi_x = 0.5_wp - x/6._wp
                            else
                                psi_x = (x - 1._wp + e_x)/(x*x)
                            end if
                            dif_tbc(idir, iside)%U(i, t1, t2) = &
                                e_x*(dif_tbc(idir, iside)%U(i, t1, t2) + tbc_d2*(D12*lapU(i, t1, t2) - advU(i, t1, t2))) + &
                                tbc_d2*(h_a*phi_x + (h_b - h_a)*psi_x)
                        end do

                        dif_tbc(idir, iside)%f_a(t1, t2) = dif_tbc(idir, iside)%f_b(t1, t2)
                        dif_tbc(idir, iside)%f_b(t1, t2) = dif_tbc(idir, iside)%f_s(1, t1, t2)
                        dif_tbc(idir, iside)%G_a(t1, t2) = dif_tbc(idir, iside)%G_b(t1, t2)
                        dif_tbc(idir, iside)%G_b(t1, t2) = G(1, t1, t2)
                    end do
                end do

                deallocate (lapU, advU)
            end if

            ! Boundary fluxes
            do t2 = 0, n2
                do t1 = 0, n1
                    call s_tbc_cell(idir, iside, t1, t2, kc, lc, qc, kf, lf, qf)
                    alpha_g = q_prim_vf(advg_idx)%sf(kc, lc, qc)
                    rho_c = q_prim_vf(Dif_idx(1))%sf(kc, lc, qc) + q_prim_vf(Dif_idx(2))%sf(kc, lc, qc)
                    gas_present = alpha_g > small_num_dif .and. rho_c > 0._wp

                    if (.not. gas_present) then
                        j_src_n(Dif_idx(1))%sf(kf, lf, qf) = 0._wp
                        j_src_n(Dif_idx(2))%sf(kf, lf, qf) = 0._wp
                        j_src_n(E_idx)%sf(kf, lf, qf) = 0._wp
                        dif_tbc(idir, iside)%V_out(t1, t2) = 0._wp
                        cycle
                    end if

                    ! (d/dt + u_t.grad_t - D Lap_t + c)^(1/2) f: history plus exact local part
                    c = un(t1, t2)**2/(4._wp*D12)
                    hist = 0._wp
                    do i = 1, n_soe
                        x = (soe_s(i) + c)*L
                        if (x < 700._wp) hist = hist + soe_w(i)*exp(-x)*dif_tbc(idir, iside)%U(i, t1, t2)
                    end do

                    fp = (dif_tbc(idir, iside)%f_s(1, t1, t2) - dif_tbc(idir, iside)%f_a(t1, t2))/L
                    h_s = fp + c*dif_tbc(idir, iside)%f_s(1, t1, t2) + G(1, t1, t2)
                    h_a = fp + c*dif_tbc(idir, iside)%f_a(t1, t2) + dif_tbc(idir, iside)%G_a(t1, t2)
                    call s_tbc_moments(c, L, M0, M1)
                    op_f = (hist + h_s*M0 - (h_s - h_a)*M1/L)/sqrt(pi)

                    ! Outward mass flux of species 1 from the mole-fraction flux,
                    ! grad Y1 = W1 W2/W**2 grad X1, gated by the gas fraction like
                    ! the interior faces (g_f = alpha_g**2 for two equal neighbours)
                    rhog = rho_c/alpha_g
                    Y1 = q_prim_vf(Dif_idx(1))%sf(kc, lc, qc)/rho_c
                    W_b = 1._wp/(Y1/Ws(1) + (1._wp - Y1)/Ws(2))
                    J_out = alpha_g**2*rhog*Ws(1)*Ws(2)/W_b**2* &
                            (sqrt(D12)*op_f - 0.5_wp*un(t1, t2)*dif_tbc(idir, iside)%f_s(1, t1, t2))

                    ! Outward volume flux V = W/rho_g sum_i j_i/W_i, used for u* at the next call
                    dif_tbc(idir, iside)%V_out(t1, t2) = J_out*W_b/rhog*(1._wp/Ws(1) - 1._wp/Ws(2))
                    T_b = q_prim_vf(E_idx)%sf(kc, lc, qc)*W_b/(rhog*R_univ)
                    do i = 1, 2
                        h_e(i) = qvs(Dif_idx(i)) + cps(i)*T_b
                    end do

                    j_src_n(Dif_idx(1))%sf(kf, lf, qf) = sgn*J_out
                    j_src_n(Dif_idx(2))%sf(kf, lf, qf) = -sgn*J_out
                    j_src_n(E_idx)%sf(kf, lf, qf) = sgn*(h_e(1) - h_e(2))*J_out
                end do
            end do
        end do

        deallocate (G, lapf, advf, ut, un)

    end subroutine s_apply_diffusion_tbc

    !> Boundary cell (kc, lc, qc) and boundary face (kf, lf, qf) of face cell
        !! (t1, t2) on face (idir, iside). j_src_n(k) is the flux at face k+1/2.
    subroutine s_tbc_cell(idir, iside, t1, t2, kc, lc, qc, kf, lf, qf)

        integer, intent(in) :: idir, iside, t1, t2
        integer, intent(out) :: kc, lc, qc, kf, lf, qf

        select case (idir)
        case (1)
            kc = merge(0, m, iside == 1); lc = t1; qc = t2
            kf = merge(-1, m, iside == 1); lf = lc; qf = qc
        case (2)
            kc = t1; lc = merge(0, n, iside == 1); qc = t2
            kf = kc; lf = merge(-1, n, iside == 1); qf = qc
        case default
            kc = t1; lc = t2; qc = merge(0, p, iside == 1)
            kf = kc; lf = lc; qf = merge(-1, p, iside == 1)
        end select

    end subroutine s_tbc_cell

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

                        ! No gas mass: no composition to split alpha_g by
                        if (rho <= 0._wp) then
                            do i = 1, Dif_size
                                q_cons_vf(advxb + Dif_idx(i) - 1)%sf(x, y, z) = 0._wp
                            end do
                            cycle
                        end if

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

                        ! No gas mass: no composition to split alpha_g by
                        if (rho <= 0._wp) then
                            do i = 1, Dif_size
                                q_rs_vf(x, y, z, advxb + Dif_idx(i) - 1) = 0._wp
                            end do
                            cycle
                        end if

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

        integer :: idir, iside

        do idir = 1, 3
            do iside = 1, 2
                if (allocated(dif_tbc(idir, iside)%U)) deallocate (dif_tbc(idir, iside)%U)
                if (allocated(dif_tbc(idir, iside)%f_a)) deallocate (dif_tbc(idir, iside)%f_a)
                if (allocated(dif_tbc(idir, iside)%f_b)) deallocate (dif_tbc(idir, iside)%f_b)
                if (allocated(dif_tbc(idir, iside)%f_s)) deallocate (dif_tbc(idir, iside)%f_s)
                if (allocated(dif_tbc(idir, iside)%G_a)) deallocate (dif_tbc(idir, iside)%G_a)
                if (allocated(dif_tbc(idir, iside)%G_b)) deallocate (dif_tbc(idir, iside)%G_b)
                if (allocated(dif_tbc(idir, iside)%X_inf)) deallocate (dif_tbc(idir, iside)%X_inf)
                if (allocated(dif_tbc(idir, iside)%V_out)) deallocate (dif_tbc(idir, iside)%V_out)
            end do
        end do
        if (allocated(soe_w)) deallocate (soe_w)
        if (allocated(soe_s)) deallocate (soe_s)

        @:DEALLOCATE(Ds)
        @:DEALLOCATE(Ws)
        @:DEALLOCATE(cps)
        @:DEALLOCATE(w_interp2)
        @:DEALLOCATE(w_grad2)
        if (dif_order == 4) then
            @:DEALLOCATE(w_interp4)
            @:DEALLOCATE(w_grad4)
        end if

    end subroutine s_finalize_diffusion_module


end module m_diffusion