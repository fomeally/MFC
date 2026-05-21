#:def Hardcoded1DVariables()
    ! Place any declaration of intermediate variables here
    real(wp) :: xloc, rho0, p0, dp, Ldom
    real(wp) :: dxloc, sinc_arg, sinc_fac
    real(wp) :: theta, ppert
#:enddef

#:def Hardcoded1D()

    select case (patch_icpp(patch_id)%hcid)
    case (100)
        xloc = x_cc(i)
        Ldom = x_boundary%end - x_boundary%beg
        dxloc = Ldom/real(m + 1, wp)

        rho0 = 1.0_wp
        p0   = 1.0e5_wp
        dp   = 1.0e2_wp

        sinc_arg = pi*dxloc/Ldom
        sinc_fac = sin(sinc_arg)/sinc_arg

        theta = 2.0_wp*pi*(xloc - x_boundary%beg)/Ldom
        ppert = dp*sinc_fac*sin(theta)

        q_prim_vf(1)%sf(i,0,0) = rho0
        q_prim_vf(2)%sf(i,0,0) = 0.0_wp
        q_prim_vf(3)%sf(i,0,0) = p0 + ppert
        q_prim_vf(4)%sf(i,0,0) = 1.0_wp

    case default
        call s_int_to_str(patch_id, iStr)
        call s_mpi_abort("Invalid hcid specified for patch "//trim(iStr))
    end select

#:enddef
