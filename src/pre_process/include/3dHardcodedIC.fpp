#:def Hardcoded3DVariables()
    ! Place any declaration of intermediate variables here

    real(wp) :: rhoH, rhoL, pRef, pInt, h, lam, wl, amp, intH, alph

    real(wp) :: eps

    real(wp) :: r_d2, B_d2, phi_d2, Y_d2, T_d2, rhov_d2, rhoa_d2, alpv_d2, eta_d2, shift_d2

    eps = 1e-9_wp
#:enddef

#:def Hardcoded3D()

    select case (patch_icpp(patch_id)%hcid)
    case (300) ! Rayleigh-Taylor instability
        rhoH = 3._wp
        rhoL = 1._wp
        pRef = 1.e5_wp
        pInt = pRef
        h = 0.7_wp
        lam = 0.2_wp
        wl = 2._wp*pi/lam
        amp = 0.025_wp/wl

        intH = amp*(sin(2._wp*pi*x_cc(i)/lam - pi/2._wp) + sin(2._wp*pi*z_cc(k)/lam - pi/2._wp)) + h

        alph = 5e-1_wp*(1._wp + tanh((y_cc(j) - intH)/2.5e-3_wp))

        if (alph < eps) alph = eps
        if (alph > 1._wp - eps) alph = 1._wp - eps

        if (y_cc(j) > intH) then
            q_prim_vf(advxb)%sf(i, j, k) = alph
            q_prim_vf(advxe)%sf(i, j, k) = 1._wp - alph
            q_prim_vf(contxb)%sf(i, j, k) = alph*rhoH
            q_prim_vf(contxe)%sf(i, j, k) = (1._wp - alph)*rhoL
            q_prim_vf(E_idx)%sf(i, j, k) = pref + rhoH*9.81_wp*(1.2_wp - y_cc(j))
        else
            q_prim_vf(advxb)%sf(i, j, k) = alph
            q_prim_vf(advxe)%sf(i, j, k) = 1._wp - alph
            q_prim_vf(contxb)%sf(i, j, k) = alph*rhoH
            q_prim_vf(contxe)%sf(i, j, k) = (1._wp - alph)*rhoL
            pInt = pref + rhoH*9.81_wp*(1.2_wp - intH)
            q_prim_vf(E_idx)%sf(i, j, k) = pInt + rhoL*9.81_wp*(intH - y_cc(j))
        end if

    case (301) ! (3D lung geometry in X direction, |sin(*)+sin(*)|)
        h = 0.0_wp
        lam = 1.0_wp
        amp = patch_icpp(patch_id)%a(2)
        intH = amp*abs((sin(2*pi*y_cc(j)/lam - pi/2) + sin(2*pi*z_cc(k)/lam - pi/2)) + h)
        if (x_cc(i) > intH) then
            q_prim_vf(contxb)%sf(i, j, k) = patch_icpp(1)%alpha_rho(1)
            q_prim_vf(contxe)%sf(i, j, k) = patch_icpp(1)%alpha_rho(2)
            q_prim_vf(E_idx)%sf(i, j, k) = patch_icpp(1)%pres
            q_prim_vf(advxb)%sf(i, j, k) = patch_icpp(1)%alpha(1)
            q_prim_vf(advxe)%sf(i, j, k) = patch_icpp(1)%alpha(2)
        end if

    case (302, 303) ! Droplet of radius patch_icpp%radius at the patch centroid in the quasi-steady
        ! d^2-law gas (Palmore & Desjardins 2019, Eqs. 82-83). With phi = (1 + B)**(1 - R/r),
        ! 1 - Y = (1 - Y_s)*phi and T = T_s + (T_inf - T_s)*(phi - 1)/B, B = (Y_s - Y_inf)/(1 - Y_s);
        ! inside r <= R the gas takes the surface state (T_s, Y_s). Fluid 1 is the liquid, fluid 2
        ! the vapor and fluid 3 the air, the gases ideal at the patch pressure, rho_i = p*gamma_i/(cv_i*T).
        ! Liquid and gas are blended with the sphere patch's tanh, eta = 1 in the liquid; it is
        ! built here because a smoothed sphere patch would refill this patch's cells with its constants.
        ! a(2) = T_s, a(3) = Y_s, a(4) = T_inf, a(5) = Y_inf, a(6) = liquid density, a(7) = smooth_coeff
        ! hcid 303 is the same, but the tanh is moved inward so the liquid volume fraction is a(8) at r = R
        ! (the outer edge of the band, not its centre, sits at R); the gas profile stays anchored at R
        shift_d2 = 0._wp
        if (patch_icpp(patch_id)%hcid == 303) then
            shift_d2 = atanh(1._wp - 2._wp*patch_icpp(patch_id)%a(8))*min(dx, dy, dz)/patch_icpp(patch_id)%a(7)
        end if
        r_d2 = sqrt((x_cc(i) - x_centroid)**2 + (cart_y - y_centroid)**2 + (cart_z - z_centroid)**2)
        B_d2 = (patch_icpp(patch_id)%a(3) - patch_icpp(patch_id)%a(5))/(1._wp - patch_icpp(patch_id)%a(3))
        phi_d2 = (1._wp + B_d2)**(1._wp - min(patch_icpp(patch_id)%radius/max(r_d2, 1e-30_wp), 1._wp))
        Y_d2 = 1._wp - (1._wp - patch_icpp(patch_id)%a(3))*phi_d2
        T_d2 = patch_icpp(patch_id)%a(2) + (patch_icpp(patch_id)%a(4) - patch_icpp(patch_id)%a(2))*(phi_d2 - 1._wp)/B_d2

        rhov_d2 = q_prim_vf(E_idx)%sf(i, j, k)*fluid_pp(2)%gamma/(fluid_pp(2)%cv*T_d2)
        rhoa_d2 = q_prim_vf(E_idx)%sf(i, j, k)*fluid_pp(3)%gamma/(fluid_pp(3)%cv*T_d2)
        alpv_d2 = (Y_d2/rhov_d2)/(Y_d2/rhov_d2 + (1._wp - Y_d2)/rhoa_d2)

        eta_d2 = 0.5_wp - 0.5_wp*tanh(patch_icpp(patch_id)%a(7)/min(dx, dy, dz)*(r_d2 - patch_icpp(patch_id)%radius + shift_d2))

        q_prim_vf(contxb)%sf(i, j, k) = eta_d2*patch_icpp(patch_id)%a(6)
        q_prim_vf(contxb + 1)%sf(i, j, k) = (1._wp - eta_d2)*alpv_d2*rhov_d2
        q_prim_vf(contxb + 2)%sf(i, j, k) = (1._wp - eta_d2)*(1._wp - alpv_d2)*rhoa_d2
        q_prim_vf(advxb)%sf(i, j, k) = eta_d2
        q_prim_vf(advxb + 1)%sf(i, j, k) = (1._wp - eta_d2)*alpv_d2
        q_prim_vf(advxb + 2)%sf(i, j, k) = (1._wp - eta_d2)*(1._wp - alpv_d2)

        ! Put your variable assignments here
    case default
        call s_int_to_str(patch_id, iStr)
        call s_mpi_abort("Invalid hcid specified for patch "//trim(iStr))
    end select

#:enddef
