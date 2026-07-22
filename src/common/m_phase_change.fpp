!> energies (6-eqn to 4-eqn) equilibrium through an infinitely fast (algebraic)
!> procedure.

#:include 'macros.fpp'

module m_phase_change

#ifndef MFC_POST_PROCESS

    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters

    use m_mpi_proxy            !< Message passing interface (MPI) module proxy

    use m_variables_conversion !< State variables type conversion procedures

    use ieee_arithmetic

    implicit none

    private; 
    public :: s_initialize_phasechange_module, &
              s_relaxation_solver, &
              s_infinite_relaxation_k, &
              s_get_Ysat, &
              s_compute_residual, &
              s_solve_pT_alpha_from_masses, &
              s_pS_residual, &
              s_finalize_relaxation_solver_module

    !> @name Parameters for the first order transition phase change
    !> @{
    integer, parameter :: max_iter = 1e8_wp        !< max # of iterations
    real(wp), parameter :: pCr = 4.94e7_wp   !< Critical water pressure
    real(wp), parameter :: TCr = 385.05_wp + 273.15_wp  !< Critical water temperature
    real(wp), parameter :: mixM = 1.0e-8_wp !< threshold for 'mixture cell'. If Y < mixM, phase change does not happen
    integer, parameter :: lp = 1    !< index for the liquid phase of the reacting fluid
    integer, parameter :: vp = 2    !< index for the vapor phase of the reacting fluid
    !> @}

    !> @name Gibbs free energy phase change parameters
    !> @{
    real(wp) :: A, B, C, D
    !> @}

    !$acc declare create(max_iter,pCr,TCr,mixM,lp,vp,A,B,C,D)

contains

    !> This subroutine should dispatch to the correct relaxation solver based
        !!      some parameter. It replaces the procedure pointer, which CCE
        !!      is breaking on.
    subroutine s_relaxation_solver(q_cons_vf)
        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf
        ! This is empty because in current master the procedure pointer
        ! was never assigned
        @:ASSERT(.false., "s_relaxation_solver called but it currently does nothing")
    end subroutine s_relaxation_solver

    !>  The purpose of this subroutine is to initialize the phase change module
        !!      by setting the parameters needed for phase change and
        !!      selecting the phase change module that will be used
        !!      (pT- or pTg-equilibrium)
    subroutine s_initialize_phasechange_module
        ! variables used in the calculation of the saturation curves for fluids 1 and 2
        A = (gs_min(lp)*cvs(lp) - gs_min(vp)*cvs(vp) &
             + qvps(vp) - qvps(lp))/((gs_min(vp) - 1.0_wp)*cvs(vp))

        B = (qvs(lp) - qvs(vp))/((gs_min(vp) - 1.0_wp)*cvs(vp))

        C = (gs_min(vp)*cvs(vp) - gs_min(lp)*cvs(lp)) &
            /((gs_min(vp) - 1.0_wp)*cvs(vp))

        D = ((gs_min(lp) - 1.0_wp)*cvs(lp)) &
            /((gs_min(vp) - 1.0_wp)*cvs(vp))

    end subroutine s_initialize_phasechange_module

    !>  This subroutine is created to activate either the pT- (N fluids) or the
        !!      pTg-equilibrium (2 fluids for g-equilibrium)
        !!      model, also considering mass depletion, depending on the incoming
        !!      state conditions.
        !!  @param q_cons_vf Cell-average conservative variables
    subroutine s_infinite_relaxation_k(q_cons_vf)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf
        real(wp) :: pS, pSOV, pSSL !< equilibrium pressure for mixture, overheated vapor, and subcooled liquid
        real(wp) :: TS, TSOV, TSSL, TSatOV, TSatSL !< equilibrium temperature for mixture, overheated vapor, and subcooled liquid. Saturation Temperatures at overheated vapor and subcooled liquid
        real(wp) :: rhoe, dynE, rhos !< total internal energy, kinetic energy, and total entropy
        real(wp) :: rho, rM, m1, m2, MCT !< total density, total reacting mass, individual reacting masses
        real(wp) :: TvF !< total volume fraction
        real(wp) :: m_g, dm, W_g, Y_sat, alpha_g !< auxiliary variables for the mass depletion procedure
        real(wp), dimension(Dif_size) :: Y_g !< auxiliary variable for the mass depletion procedure
        logical :: valid_Ysat !< auxiliary variable for the mass depletion procedure

        !$acc declare create(pS, pSOV, pSSL, TS, TSOV, TSatOV, TSatSL, TSSL, rhoe, dynE, rhos, rho, rM, m1, m2, MCT, TvF)

        real(wp), dimension(num_fluids) :: p_infOV, p_infpT, p_infSL, sk, hk, gk, ek, rhok, m_i, alpha_i !< auxiliary variables for the pT- and pTg-equilibrium solvers
        real(wp), dimension(num_fluids) :: m_trial
        real(wp) :: R, R_plus, R_minus, m_prime, dm_plus, dm_minus, dRdm !< auxiliary variables for the mass depletion procedure
        !< Generic loop iterators
        integer :: i, j, k, l, iter, max_iter

        !$acc declare create(p_infOV, p_infpT, p_infSL, sk, hk, gk, ek, rhok)

        max_iter = 20
        ! starting equilibrium solver
        !$acc parallel loop collapse(3) gang vector default(present) private(p_infOV, p_infpT, p_infSL, sk, hk, gk, ek, rhok,pS, pSOV, pSSL, TS, TSOV, TSatOV, TSatSL, TSSL, rhoe, dynE, rhos, rho, rM, m1, m2, MCT, TvF)
        do j = 0, m
            do k = 0, n
                do l = 0, p

                    rho = 0.0_wp; TvF = 0.0_wp
                    !$acc loop seq
                    do i = 1, num_fluids

                        ! Mixture density
                        rho = rho + q_cons_vf(i + contxb - 1)%sf(j, k, l)

                        ! Total Volume Fraction
                        TvF = TvF + q_cons_vf(i + advxb - 1)%sf(j, k, l)

                    end do

                    ! calculating the total reacting mass for the phase change process. By hypothesis, this should not change
                    ! throughout the phase-change process.
                    rM = q_cons_vf(lp + contxb - 1)%sf(j, k, l) + q_cons_vf(vp + contxb - 1)%sf(j, k, l)

                    ! correcting negative (reacting) mass fraction values in case they happen
                    call s_correct_partial_densities(MCT, q_cons_vf, rM, j, k, l)

                    ! fixing m1 and m2 AFTER correcting the partial densities. Note that these values must be stored for the phase
                    ! change process that will happen a posteriori
                    m1 = q_cons_vf(lp + contxb - 1)%sf(j, k, l)

                    m2 = q_cons_vf(vp + contxb - 1)%sf(j, k, l)

                    ! kinetic energy as an auxiliary variable to the calculation of the total internal energy
                    dynE = 0.0_wp
                    !$acc loop seq
                    do i = momxb, momxe

                        dynE = dynE + 5.0e-1_wp*q_cons_vf(i)%sf(j, k, l)**2/rho

                    end do

                    ! calculating the total energy that MUST be preserved throughout the pT- and pTg-relaxation procedures
                    ! at each of the cells. The internal energy is calculated as the total energy minus the kinetic
                    ! energy to preserved its value at sharp interfaces
                    rhoe = q_cons_vf(E_idx)%sf(j, k, l) - dynE

                    ! Calling pT-equilibrium for either finishing phase-change module, or as an IC for the pTg-equilibrium
                    ! for this case, MFL cannot be either 0 or 1, so I chose it to be 2
                    call s_infinite_pt_relaxation_k(j, k, l, 2, pS, p_infpT, rM, q_cons_vf, rhoe, TS)

                    ! check if pTg-equilibrium is required
                    ! NOTE that NOTHING else needs to be updated OTHER than the individual partial densities
                    ! given the outputs from the pT- and pTg-equilibrium solvers are just p and one of the partial masses
                    ! (pTg- case)
                    if ((relax_model == 6) .and. ((q_cons_vf(lp + contxb - 1)%sf(j, k, l) > mixM*rM) &
                                                  .and. (q_cons_vf(vp + contxb - 1)%sf(j, k, l) > mixM*rM)) &
                        .and. (pS < pCr) .and. (TS < TCr)) then

                        ! Checking if phase change is needed, by checking whether the final solution is either subcoooled
                        ! liquid or overheated vapor.

                        ! overheated vapor case
                        ! depleting the mass of liquid
                        q_cons_vf(lp + contxb - 1)%sf(j, k, l) = mixM*rM

                        ! tranferring the total mass to vapor
                        q_cons_vf(vp + contxb - 1)%sf(j, k, l) = (1.0_wp - mixM)*rM

                        ! calling pT-equilibrium for overheated vapor, which is MFL = 0
                        call s_infinite_pt_relaxation_k(j, k, l, 0, pSOV, p_infOV, rM, q_cons_vf, rhoe, TSOV)

                        ! calculating Saturation temperature
                        call s_TSat(pSOV, TSatOV, TSOV)

                        ! subcooled liquid case
                        ! tranferring the total mass to liquid
                        q_cons_vf(lp + contxb - 1)%sf(j, k, l) = (1.0_wp - mixM)*rM

                        ! depleting the mass of vapor
                        q_cons_vf(vp + contxb - 1)%sf(j, k, l) = mixM*rM

                        ! calling pT-equilibrium for subcooled liquid, which is MFL = 1
                        call s_infinite_pt_relaxation_k(j, k, l, 1, pSSL, p_infSL, rM, q_cons_vf, rhoe, TSSL)

                        ! calculating Saturation temperature
                        call s_TSat(pSSL, TSatSL, TSSL)

                        ! checking the conditions for overheated vapor and subcooled liquide
                        if (TSOV > TSatOV) then

                            ! Assigning pressure
                            pS = pSOV

                            ! Assigning Temperature
                            TS = TSOV

                            ! correcting the liquid partial density
                            q_cons_vf(lp + contxb - 1)%sf(j, k, l) = mixM*rM

                            ! correcting the vapor partial density
                            q_cons_vf(vp + contxb - 1)%sf(j, k, l) = (1.0_wp - mixM)*rM

                        elseif (TSSL < TSatSL) then

                            ! Assigning pressure
                            pS = pSSL

                            ! Assigning Temperature
                            TS = TSSL

                            ! correcting the liquid partial density
                            q_cons_vf(lp + contxb - 1)%sf(j, k, l) = (1.0_wp - mixM)*rM

                            ! correcting the vapor partial density
                            q_cons_vf(vp + contxb - 1)%sf(j, k, l) = mixM*rM

                        else
                            ! returning partial pressures to what they were from the homogeneous solver
                            ! liquid
                            q_cons_vf(lp + contxb - 1)%sf(j, k, l) = m1

                            ! vapor
                            q_cons_vf(vp + contxb - 1)%sf(j, k, l) = m2

                            ! calling the pTg-equilibrium solver
                            call s_infinite_ptg_relaxation_k(j, k, l, pS, p_infpT, rhoe, q_cons_vf, TS)

                        end if

                    else if (relax_model == 7) then
                        if (q_cons_vf(lp + advxb - 1)%sf(j, k, l) > mixM &
                            .and. q_cons_vf(advg_idx)%sf(j, k, l) > mixM) then

                                

                                !$acc loop seq
                                do i = 1, num_fluids
                                    m_i(i) = q_cons_vf(i + contxb - 1)%sf(j, k, l)
                                    alpha_i(i) = q_cons_vf(i + advxb - 1)%sf(j, k, l)
                                end do

                                alpha_g = q_cons_vf(advg_idx)%sf(j, k, l)

                                m_g = 0.0_wp
                                !$acc loop seq
                                do i = 1, Dif_size
                                    m_g = m_g + m_i(Dif_idx(i))
                                end do

                                call s_get_Ysat(m_i, TS, pS, m_g, Y_sat, valid_Ysat)

                                if (.not. valid_Ysat) cycle

                                dm = ( Y_sat*m_g - q_cons_vf(vp + contxb - 1)%sf(j, k, l) )/ (1.0_wp - Y_sat)

                                ! with this initial dm guess, run newton solver to find the real dm
                                do i = 1, max_iter

                                    
                                    call s_compute_residual(m_i, alpha_i, dm, rhoe, R)

                                    if (abs(R) < 1.0e-8_wp) then

                                        ! print *, "converged", i, " iterations"
                                        exit
                                    
                                    end if
                                    

                                    m_prime = 1.0e-6_wp*max(abs(dm), m_i(lp), m_i(vp), 1.0e-20_wp)
                                    dm_plus  = min(max(dm + m_prime, -m_i(vp)), m_i(lp))
                                    dm_minus = min(max(dm - m_prime, -m_i(vp)), m_i(lp))

                                    call s_compute_residual(m_i, alpha_i, dm_plus, rhoe, R_plus)

                                    call s_compute_residual(m_i, alpha_i, dm_minus, rhoe, R_minus)

                                    dRdm = (R_plus - R_minus)/(dm_plus - dm_minus)

                                    dm = dm - R/dRdm

                                end do

                                ! we converged to a solution, now we can update the partial densities and enthalpies accordingly

                                dm = max(dm, -m_i(vp))
                                dm = min(dm, m_i(lp))

                                m_trial(:) = m_i(:)

                                m_trial(lp) = m_i(lp) - dm
                                m_trial(vp) = m_i(vp) + dm

                                q_cons_vf(lp + contxb - 1)%sf(j, k, l) = m_trial(lp)
                                q_cons_vf(vp + contxb - 1)%sf(j, k, l) = m_trial(vp)

                                ! compute the new temperature and pressure and vol fracs after the mass depletion procedure
                                call s_solve_pT_alpha_from_masses(m_trial, rhoe, pS, TS, alpha_i)

                                alpha_g = 1.0_wp - alpha_i(lp)

                                q_cons_vf(advg_idx)%sf(j, k, l) = alpha_g

                                do i = 1, num_fluids
                                    q_cons_vf(i + advxb - 1)%sf(j, k, l) = alpha_i(i)
                                end do                           
                                
                        end if
                    end if

                    ! Calculations AFTER equilibrium

                    if ( relax_model /= 7) then

                        ! entropy
                        sk(1:num_fluids) = cvs(1:num_fluids)*log((TS**gs_min(1:num_fluids)) &
                                                                /((pS + ps_inf(1:num_fluids))**(gs_min(1:num_fluids) - 1.0_wp))) + qvps(1:num_fluids)

                        ! enthalpy
                        hk(1:num_fluids) = gs_min(1:num_fluids)*cvs(1:num_fluids)*TS &
                                        + qvs(1:num_fluids)

                        ! Gibbs-free energy
                        gk(1:num_fluids) = hk(1:num_fluids) - TS*sk(1:num_fluids)

                        ! densities
                        rhok(1:num_fluids) = (pS + ps_inf(1:num_fluids)) &
                                            /((gs_min(1:num_fluids) - 1)*cvs(1:num_fluids)*TS)

                        ! internal energy
                        ek(1:num_fluids) = (pS + gs_min(1:num_fluids) &
                                            *ps_inf(1:num_fluids))/(pS + ps_inf(1:num_fluids)) &
                                        *cvs(1:num_fluids)*TS + qvs(1:num_fluids)

                        ! calculating volume fractions, internal energies, and total entropy
                        rhos = 0.0_wp
                        !$acc loop seq
                        do i = 1, num_fluids

                            ! volume fractions
                            q_cons_vf(i + advxb - 1)%sf(j, k, l) = q_cons_vf(i + contxb - 1)%sf(j, k, l)/rhok(i)

                            ! alpha*rho*e
                            q_cons_vf(i + intxb - 1)%sf(j, k, l) = q_cons_vf(i + contxb - 1)%sf(j, k, l)*ek(i)

                            ! Total entropy
                            rhos = rhos + q_cons_vf(i + contxb - 1)%sf(j, k, l)*sk(i)

                        end do

                    end if
                end do
            end do
        end do

    end subroutine s_infinite_relaxation_k

    !>  This auxiliary subroutine is created to activate the pT-equilibrium for N fluids
        !!  @param j generic loop iterator for x direction
        !!  @param k generic loop iterator for y direction
        !!  @param l generic loop iterator for z direction
        !!  @param MFL flag that tells whether the fluid is pure gas (0), pure liquid (1), or a mixture (2)
        !!  @param pS equilibrium pressure at the interface
        !!  @param p_infpT stiffness for the participating fluids under pT-equilibrium
        !!  @param rM sum of the reacting masses
        !!  @param q_cons_vf Cell-average conservative variables
        !!  @param rhoe mixture energy
        !!  @param TS equilibrium temperature at the interface
    subroutine s_infinite_pt_relaxation_k(j, k, l, MFL, pS, p_infpT, rM, q_cons_vf, rhoe, TS)

#ifdef _CRAYFTN
        !DIR$ INLINEALWAYS s_infinite_pt_relaxation_k
#else
        !$acc routine seq
#endif

        ! initializing variables
        integer, intent(in) :: j, k, l, MFL
        real(wp), intent(out) :: pS
        real(wp), dimension(num_fluids), intent(out) :: p_infpT
        real(wp), intent(in) :: rM
        type(scalar_field), dimension(sys_size), intent(in) :: q_cons_vf
        real(wp), intent(in) :: rhoe
        real(wp), intent(out) :: TS
        real(wp) :: gp, gpp, hp, pO, mCP, mQ, alphaGam, alphaPi_inf !< variables for the Newton Solver
        real(wp), dimension(Dif_size) :: Y_g !< mass fractions of the gas mixture
        real(wp) :: W_g, m_g, rho_g, R_univ


        integer :: i, ns !< generic loop iterators

        R_univ = 8314.462618_wp

        ! auxiliary variables for the pT-equilibrium solver
        mCP = 0.0_wp; mQ = 0.0_wp; p_infpT = ps_inf; alphaGam = 0.0_wp; alphaPi_inf = 0.0_wp
        ! Performing tests before initializing the pT-equilibrium
        !$acc loop seq
        do i = 1, num_fluids

            ! sum of the total alpha*rho*cp of the system
            mCP = mCP + q_cons_vf(i + contxb - 1)%sf(j, k, l)*cvs(i)*gs_min(i)

            ! sum of the total alpha*rho*q of the system
            mQ = mQ + q_cons_vf(i + contxb - 1)%sf(j, k, l)*qvs(i)

            alphaGam = alphaGam + q_cons_vf(i + advxb - 1)%sf(j, k, l)*gammas(i)

            alphaPi_inf = alphaPi_inf + q_cons_vf(i + advxb - 1)%sf(j, k, l)*pi_infs(i)

        end do

        ! Checking energy constraint
        if ((rhoe - mQ - minval(p_infpT)) < 0.0_wp) then

            if ((MFL == 0) .or. (MFL == 1)) then

                ! Assigning zero values for mass depletion cases
                ! pressure
                pS = 0.0_wp

                ! temperature
                TS = 0.0_wp

                return
            end if

        end if

        ! calculating initial estimate for pressure in the pT-relaxation procedure. I will also use this variable to
        ! iterate over the Newton's solver
        pO = 0.0_wp

        ! Maybe improve this condition afterwards. As long as the initial guess is in between -min(ps_inf)
        ! and infinity, a solution should be able to be found.
        pS = 1.0e4_wp

        ! Newton Solver for the pT-equilibrium
        ns = 0

        if (relax_model /= 7) then
            ! change this relative error metric. 1e4_wp is just arbitrary
            do while ((abs(pS - pO) > palpha_eps) .and. (abs((pS - pO)/pO) > palpha_eps/1e4_wp) .or. (ns == 0))

                ! increasing counter
                ns = ns + 1

                ! updating old pressure
                pO = pS

                ! updating functions used in the Newton's solver
                gpp = 0.0_wp; gp = 0.0_wp; hp = 0.0_wp
                !$acc loop seq
                do i = 1, num_fluids

                    gp = gp + (gs_min(i) - 1.0_wp)*q_cons_vf(i + contxb - 1)%sf(j, k, l)*cvs(i) &
                        *(rhoe + pS - mQ)/(mCP*(pS + p_infpT(i)))

                    gpp = gpp + (gs_min(i) - 1.0_wp)*q_cons_vf(i + contxb - 1)%sf(j, k, l)*cvs(i) &
                        *(p_infpT(i) - rhoe + mQ)/(mCP*(pS + p_infpT(i))**2)

                end do

                hp = 1.0_wp/(rhoe + pS - mQ) + 1.0_wp/(pS + minval(p_infpT))

                ! updating common pressure for the newton solver
                pS = pO + ((1.0_wp - gp)/gpp)/(1.0_wp - (1.0_wp - gp + abs(1.0_wp - gp)) &
                                            /(2.0_wp*gpp)*hp)
            end do
        else ! vapor saturation model
            pS = ( rhoe - alphaPi_inf - mQ )/ alphaGam
        end if

        if (relax_model /= 7) then
            ! common temperature
            TS = (rhoe + pS - mQ)/mCP
        ! elseif (relax_model == 7 .and. q_cons_vf(advg_idx)%sf(j, k, l) > mixM*1e5_wp) then

        !     m_g = 0.0_wp
        !     do i = 1, Dif_size
        !         m_g = m_g + q_cons_vf(Dif_idx(i) + contxb - 1)%sf(j, k, l)
        !     end do

        !     rho_g = m_g / q_cons_vf(advg_idx)%sf(j, k, l)

        !     Y_g = 0.0_wp
        !     do i = 1, Dif_size
        !         Y_g(i) = q_cons_vf(Dif_idx(i) + contxb - 1)%sf(j, k, l) / m_g
        !     end do

        !     W_g = 0.0_wp
        !     do i = 1, Dif_size
        !         W_g = W_g + Y_g(i) / fluid_pp(Dif_idx(i))%W
        !     end do

        !     W_g = 1.0_wp / W_g

        !     TS = pS * W_g / (rho_g * R_univ)
        !     if (TS > 300.0_wp) then
        !         print *, "TS = ", TS, " pS = ", pS, " rho_g = ", rho_g, " W_g = ", W_g
        !     end if
        else 
            TS = (rhoe + pS - mQ)/mCP
            ! TS = 298.0_wp
        end if

    end subroutine s_infinite_pt_relaxation_k

    !>  This auxiliary subroutine is created to activate the pTg-equilibrium for N fluids under pT
        !!      and 2 fluids under pTg-equilibrium. There is a final common p and T during relaxation
        !!  @param j generic loop iterator for x direction
        !!  @param k generic loop iterator for y direction
        !!  @param l generic loop iterator for z direction
        !!  @param pS equilibrium pressure at the interface
        !!  @param p_infpT stiffness for the participating fluids under pT-equilibrium
        !!  @param rhoe mixture energy
        !!  @param q_cons_vf Cell-average conservative variables
        !!  @param TS equilibrium temperature at the interface
    subroutine s_infinite_ptg_relaxation_k(j, k, l, pS, p_infpT, rhoe, q_cons_vf, TS)

#ifdef _CRAYFTN
        !DIR$ INLINEALWAYS s_infinite_ptg_relaxation_k
#else
        !$acc routine seq
#endif

        integer, intent(in) :: j, k, l
        real(wp), intent(inout) :: pS
        real(wp), dimension(num_fluids), intent(in) :: p_infpT
        real(wp), intent(in) :: rhoe
        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf
        real(wp), intent(inout) :: TS

        real(wp), dimension(num_fluids) :: p_infpTg !< stiffness for the participating fluids for pTg-equilibrium
        real(wp), dimension(2, 2) :: Jac, InvJac, TJac !< matrices for the Newton Solver
        real(wp), dimension(2) :: R2D, DeltamP !< residual and correction array
        real(wp) :: Om ! underrelaxation factor
        real(wp) :: mCP, mCPD, mCVGP, mCVGP2, mQ, mQD ! auxiliary variables for the pTg-solver

        !< Generic loop iterators
        integer :: i, ns
        ! pTg-equilibrium solution procedure
        ! Newton Solver parameters
        ! counter
        ns = 0

        ! Relaxation factor
        Om = 1.0e-3_wp

        p_infpTg = p_infpT

        if (((pS < 0.0_wp) .and. ((q_cons_vf(lp + contxb - 1)%sf(j, k, l) &
                                   + q_cons_vf(vp + contxb - 1)%sf(j, k, l)) > ((rhoe &
                                                                                 - gs_min(lp)*ps_inf(lp)/(gs_min(lp) - 1))/qvs(lp)))) .or. &
            ((pS >= 0.0_wp) .and. (pS < 1.0e-1_wp))) then

            ! improve this initial condition
            pS = 1.0e4_wp

        end if

        ! Loop until the solution for F(X) is satisfied
        ! Check whether I need to use both absolute and relative values
        ! for the residual, and how to do it adequately.
        ! Dummy guess to start the pTg-equilibrium problem.
        ! improve this initial condition
        R2D(1) = 0.0_wp; R2D(2) = 0.0_wp
        DeltamP(1) = 0.0_wp; DeltamP(2) = 0.0_wp
        do while (((sqrt(R2D(1)**2 + R2D(2)**2) > ptgalpha_eps) &
                   .and. ((sqrt(R2D(1)**2 + R2D(2)**2)/rhoe) > (ptgalpha_eps/1e6_wp))) &
                  .or. (ns == 0))

            ! Updating counter for the iterative procedure
            ns = ns + 1

            ! Auxiliary variables to help in the calculation of the residue
            mCP = 0.0_wp; mCPD = 0.0_wp; mCVGP = 0.0_wp; mCVGP2 = 0.0_wp; mQ = 0.0_wp; mQD = 0.0_wp
            ! Those must be updated through the iterations, as they either depend on
            ! the partial masses for all fluids, or on the equilibrium pressure
            !$acc loop seq
            do i = 1, num_fluids

                ! sum of the total alpha*rho*cp of the system
                mCP = mCP + q_cons_vf(i + contxb - 1)%sf(j, k, l) &
                      *cvs(i)*gs_min(i)

                ! sum of the total alpha*rho*q of the system
                mQ = mQ + q_cons_vf(i + contxb - 1)%sf(j, k, l)*qvs(i)

                ! These auxiliary variables now need to be updated, as the partial densities now
                ! vary at every iteration
                if ((i /= lp) .and. (i /= vp)) then

                    mCVGP = mCVGP + q_cons_vf(i + contxb - 1)%sf(j, k, l) &
                            *cvs(i)*(gs_min(i) - 1)/(pS + ps_inf(i))

                    mCVGP2 = mCVGP2 + q_cons_vf(i + contxb - 1)%sf(j, k, l) &
                             *cvs(i)*(gs_min(i) - 1)/((pS + ps_inf(i))**2)

                    mQD = mQD + q_cons_vf(i + contxb - 1)%sf(j, k, l)*qvs(i)

                    ! sum of the total alpha*rho*cp of the system
                    mCPD = mCPD + q_cons_vf(i + contxb - 1)%sf(j, k, l)*cvs(i) &
                           *gs_min(i)

                end if

            end do

            ! calculating the (2D) Jacobian Matrix used in the solution of the pTg-quilibrium model
            call s_compute_jacobian_matrix(InvJac, j, Jac, k, l, mCPD, mCVGP, mCVGP2, pS, q_cons_vf, TJac)

            ! calculating correction array for Newton's method
            DeltamP = -1.0_wp*matmul(InvJac, R2D)

            ! updating two reacting 'masses'. Recall that inert 'masses' do not change during the phase change
            ! liquid
            q_cons_vf(lp + contxb - 1)%sf(j, k, l) = q_cons_vf(lp + contxb - 1)%sf(j, k, l) + Om*DeltamP(1)

            ! gas
            q_cons_vf(vp + contxb - 1)%sf(j, k, l) = q_cons_vf(vp + contxb - 1)%sf(j, k, l) - Om*DeltamP(1)

            ! updating pressure
            pS = pS + Om*DeltamP(2)

            ! calculating residuals, which are (i) the difference between the Gibbs Free energy of the gas and the liquid
            ! and (ii) the energy before and after the phase-change process.
            call s_compute_pTg_residue(j, k, l, mCPD, mCVGP, mQD, q_cons_vf, pS, rhoe, R2D)

        end do

        ! common temperature
        TS = (rhoe + pS - mQ)/mCP
    end subroutine s_infinite_ptg_relaxation_k

    !>  This auxiliary subroutine corrects the partial densities of the REACTING fluids in case one of them is negative
        !!      but their sum is positive. Inert phases are not corrected at this moment
        !!  @param MCT partial density correction parameter
        !!  @param q_cons_vf Cell-average conservative variables
        !!  @param rM sum of the reacting masses
        !!  @param j generic loop iterator for x direction
        !!  @param k generic loop iterator for y direction
        !!  @param l generic loop iterator for z direction
    subroutine s_correct_partial_densities(MCT, q_cons_vf, rM, j, k, l)

#ifdef _CRAYFTN
        !DIR$ INLINEALWAYS s_correct_partial_densities
#else
        !$acc routine seq
#endif

        !> @name variables for the correction of the reacting partial densities
        !> @{
        real(wp), intent(out) :: MCT
        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf
        real(wp), intent(inout) :: rM
        integer, intent(in) :: j, k, l
        !> @}
        if (rM < 0.0_wp) then

            if ((q_cons_vf(lp + contxb - 1)%sf(j, k, l) >= -1.0_wp*mixM) .and. &
                (q_cons_vf(vp + contxb - 1)%sf(j, k, l) >= -1.0_wp*mixM)) then

                q_cons_vf(lp + contxb - 1)%sf(j, k, l) = 0.0_wp

                q_cons_vf(vp + contxb - 1)%sf(j, k, l) = 0.0_wp

                rM = q_cons_vf(lp + contxb - 1)%sf(j, k, l) + q_cons_vf(vp + contxb - 1)%sf(j, k, l)

            end if

        end if

        ! Defining the correction in terms of an absolute value might not be the best practice.
        ! Maybe a good way to do this is to partition the partial densities, giving a small percentage of the total reacting density
        MCT = 2*mixM

        ! correcting the partial densities of the reacting fluids. What to do for the nonreacting ones?
        if (q_cons_vf(lp + contxb - 1)%sf(j, k, l) < 0.0_wp) then

            q_cons_vf(lp + contxb - 1)%sf(j, k, l) = MCT*rM

            q_cons_vf(vp + contxb - 1)%sf(j, k, l) = (1.0_wp - MCT)*rM

        elseif (q_cons_vf(vp + contxb - 1)%sf(j, k, l) < 0.0_wp) then

            q_cons_vf(lp + contxb - 1)%sf(j, k, l) = (1.0_wp - MCT)*rM

            q_cons_vf(vp + contxb - 1)%sf(j, k, l) = MCT*rM

        end if
    end subroutine s_correct_partial_densities

    !>  This auxiliary subroutine calculates the 2 x 2 Jacobian and, its inverse and transpose
        !!      to be used in the pTg-equilibirium procedure
        !!  @param InvJac Inverse of the Jacobian Matrix
        !!  @param j generic loop iterator for x direction
        !!  @param Jac Jacobian Matrix
        !!  @param k generic loop iterator for y direction
        !!  @param l generic loop iterator for z direction
        !!  @param mCPD  sum of the total alpha*rho*cp
        !!  @param mCVGP auxiliary variable for the calculation of the matrices: alpha*rho*cv*(g-1)/press
        !!  @param mCVGP2 auxiliary variable for the calculation of the matrices: alpha*rho*cv*(g-1)/press^2
        !!  @param pS equilibrium pressure at the interface
        !!  @param q_cons_vf Cell-average conservative variables
        !!  @param TJac Transpose of the Jacobian Matrix
    subroutine s_compute_jacobian_matrix(InvJac, j, Jac, k, l, mCPD, mCVGP, mCVGP2, pS, q_cons_vf, TJac)

#ifdef _CRAYFTN
        !DIR$ INLINEALWAYS s_compute_jacobian_matrix
#else
        !$acc routine seq
#endif

        real(wp), dimension(2, 2), intent(out) :: InvJac
        integer, intent(in) :: j
        real(wp), dimension(2, 2), intent(out) :: Jac
        integer, intent(in) :: k, l
        real(wp), intent(in) :: mCPD, mCVGP, mCVGP2, pS
        type(scalar_field), dimension(sys_size), intent(in) :: q_cons_vf
        real(wp), dimension(2, 2), intent(out) :: TJac

        real(wp) :: ml, mT, TS, dFdT, dTdm, dTdp ! mass of the reacting fluid, total reacting mass, and auxiliary variables

        ! mass of the reacting liquid
        ml = q_cons_vf(lp + contxb - 1)%sf(j, k, l)

        ! mass of the two participating fluids
        mT = q_cons_vf(lp + contxb - 1)%sf(j, k, l) &
             + q_cons_vf(vp + contxb - 1)%sf(j, k, l)

        TS = 1/(mT*cvs(vp)*(gs_min(vp) - 1)/(pS + ps_inf(vp)) &
                + ml*(cvs(lp)*(gs_min(lp) - 1)/(pS + ps_inf(lp)) &
                      - cvs(vp)*(gs_min(vp) - 1)/(pS + ps_inf(vp))) &
                + mCVGP)

        dFdT = &
            -(cvs(lp)*gs_min(lp) - cvs(vp)*gs_min(vp))*log(TS) &
            - (qvps(lp) - qvps(vp)) &
            + cvs(lp)*(gs_min(lp) - 1)*log(pS + ps_inf(lp)) &
            - cvs(vp)*(gs_min(vp) - 1)*log(pS + ps_inf(vp))

        dTdm = -(cvs(lp)*(gs_min(lp) - 1)/(pS + ps_inf(lp)) &
                 - cvs(vp)*(gs_min(vp) - 1)/(pS + ps_inf(vp)))*TS**2

        dTdp = (mT*cvs(vp)*(gs_min(vp) - 1)/(pS + ps_inf(vp))**2 &
                + ml*(cvs(lp)*(gs_min(lp) - 1)/(pS + ps_inf(lp))**2 &
                      - cvs(vp)*(gs_min(vp) - 1)/(pS + ps_inf(vp))**2) &
                + mCVGP2)*TS**2

        ! F = (F1,F2) is the function whose roots we are looking for
        ! x = (m1, p) are the independent variables. m1 = mass of the first participant fluid, p = pressure
        ! F1 = 0 is the Gibbs free energy quality
        ! F2 = 0 is the enforcement of the thermodynamic (total - kinectic) energy
        ! dF1dm
        Jac(1, 1) = dFdT*dTdm

        ! dF1dp
        Jac(1, 2) = dFdT*dTdp + TS &
                    *(cvs(lp)*(gs_min(lp) - 1)/(pS + ps_inf(lp)) &
                      - cvs(vp)*(gs_min(vp) - 1)/(pS + ps_inf(vp)))

        ! dF2dm
        Jac(2, 1) = (qvs(vp) - qvs(lp) &
                     + (cvs(vp)*gs_min(vp) - cvs(lp)*gs_min(lp)) &
                     /(ml*(cvs(lp)*(gs_min(lp) - 1)/(pS + ps_inf(lp)) &
                           - cvs(vp)*(gs_min(vp) - 1)/(pS + ps_inf(vp))) &
                       + mT*cvs(vp)*(gs_min(vp) - 1)/(pS + ps_inf(vp)) + mCVGP) &
                     - (ml*(cvs(vp)*gs_min(vp) - cvs(lp)*gs_min(lp)) &
                        - mT*cvs(vp)*gs_min(vp) - mCPD) &
                     *(cvs(lp)*(gs_min(lp) - 1)/(pS + ps_inf(lp)) &
                       - cvs(vp)*(gs_min(vp) - 1)/(pS + ps_inf(vp))) &
                     /((ml*(cvs(lp)*(gs_min(lp) - 1)/(pS + ps_inf(lp)) &
                            - cvs(vp)*(gs_min(vp) - 1)/(pS + ps_inf(vp))) &
                        + mT*cvs(vp)*(gs_min(vp) - 1)/(pS + ps_inf(vp)) + mCVGP)**2))/1
        ! dF2dp
        Jac(2, 2) = (1 + (ml*(cvs(vp)*gs_min(vp) - cvs(lp)*gs_min(lp)) &
                          - mT*cvs(vp)*gs_min(vp) - mCPD) &
                     *(ml*(cvs(lp)*(gs_min(lp) - 1)/(pS + ps_inf(lp))**2 &
                           - cvs(vp)*(gs_min(vp) - 1)/(pS + ps_inf(vp))**2) &
                       + mT*cvs(vp)*(gs_min(vp) - 1)/(pS + ps_inf(vp))**2 + mCVGP2) &
                     /(ml*(cvs(lp)*(gs_min(lp) - 1)/(pS + ps_inf(lp)) &
                           - cvs(vp)*(gs_min(vp) - 1)/(pS + ps_inf(vp))) &
                       + mT*cvs(vp)*(gs_min(vp) - 1)/(pS + ps_inf(vp)) + mCVGP)**2)/1

        ! intermediate elements of J^{-1}
        InvJac(1, 1) = Jac(2, 2)
        InvJac(1, 2) = -1.0_wp*Jac(1, 2)
        InvJac(2, 1) = -1.0_wp*Jac(2, 1)
        InvJac(2, 2) = Jac(1, 1)

        ! elements of J^{T}
        TJac(1, 1) = Jac(1, 1)
        TJac(1, 2) = Jac(2, 1)
        TJac(2, 1) = Jac(1, 2)
        TJac(2, 2) = Jac(2, 2)

        ! dividing by det(J)
        InvJac = InvJac/(Jac(1, 1)*Jac(2, 2) - Jac(1, 2)*Jac(2, 1))

    end subroutine s_compute_jacobian_matrix

    !>  This auxiliary subroutine computes the residue of the pTg-equilibrium procedure
        !!  @param j generic loop iterator for x direction
        !!  @param k generic loop iterator for y direction
        !!  @param l generic loop iterator for z direction
        !!  @param mCPD  sum of the total alpha*rho*cp
        !!  @param mCVGP auxiliary variable for the calculation of the matrices: alpha*rho*cv*(g-1)/press
        !!  @param mQD sum of the total alpha*rho*qv
        !!  @param q_cons_vf Cell-average conservative variables
        !!  @param pS equilibrium pressure at the interface
        !!  @param rhoe mixture energy
        !!  @param R2D (2D) residue array
    subroutine s_compute_pTg_residue(j, k, l, mCPD, mCVGP, mQD, q_cons_vf, pS, rhoe, R2D)

#ifdef _CRAYFTN
        !DIR$ INLINEALWAYS s_compute_pTg_residue
#else
        !$acc routine seq
#endif

        integer, intent(in) :: j, k, l
        real(wp), intent(in) :: mCPD, mCVGP, mQD
        type(scalar_field), dimension(sys_size), intent(in) :: q_cons_vf
        real(wp), intent(in) :: pS, rhoe
        real(wp), dimension(2), intent(out) :: R2D

        real(wp) :: ml, mT, TS !< mass of the reacting liquid, total reacting mass, equilibrium temperature

        ! mass of the reacting liquid
        ml = q_cons_vf(lp + contxb - 1)%sf(j, k, l)

        ! mass of the two participating fluids
        mT = q_cons_vf(lp + contxb - 1)%sf(j, k, l) &
             + q_cons_vf(vp + contxb - 1)%sf(j, k, l)

        TS = 1/(mT*cvs(vp)*(gs_min(vp) - 1)/(pS + ps_inf(vp)) &
                + ml*(cvs(lp)*(gs_min(lp) - 1)/(pS + ps_inf(lp)) &
                      - cvs(vp)*(gs_min(vp) - 1)/(pS + ps_inf(vp))) &
                + mCVGP)

        ! Gibbs Free Energy Equality condition (DG)
        R2D(1) = TS*((cvs(lp)*gs_min(lp) - cvs(vp)*gs_min(vp)) &
                     *(1 - log(TS)) - (qvps(lp) - qvps(vp)) &
                     + cvs(lp)*(gs_min(lp) - 1)*log(pS + ps_inf(lp)) &
                     - cvs(vp)*(gs_min(vp) - 1)*log(pS + ps_inf(vp))) &
                 + qvs(lp) - qvs(vp)

        ! Constant Energy Process condition (DE)
        R2D(2) = (rhoe + pS &
                  + ml*(qvs(vp) - qvs(lp)) - mT*qvs(vp) - mQD &
                  + (ml*(gs_min(vp)*cvs(vp) - gs_min(lp)*cvs(lp)) &
                     - mT*gs_min(vp)*cvs(vp) - mCPD) &
                  /(ml*(cvs(lp)*(gs_min(lp) - 1)/(pS + ps_inf(lp)) &
                        - cvs(vp)*(gs_min(vp) - 1)/(pS + ps_inf(vp))) &
                    + mT*cvs(vp)*(gs_min(vp) - 1)/(pS + ps_inf(vp)) + mCVGP))/1

    end subroutine s_compute_pTg_residue

    !>  This auxiliary subroutine finds the Saturation temperature for a given
        !!      saturation pressure through a newton solver
        !!  @param pSat Saturation Pressure
        !!  @param TSat Saturation Temperature
        !!  @param TSIn equilibrium Temperature
    subroutine s_TSat(pSat, TSat, TSIn)

#ifdef _CRAYFTN
        !DIR$ INLINEALWAYS s_TSat
#else
        !$acc routine seq
#endif

        real(wp), intent(in) :: pSat
        real(wp), intent(out) :: TSat
        real(wp), intent(in) :: TSIn

        real(wp) :: dFdT, FT, Om !< auxiliary variables

        ! Generic loop iterators
        integer :: ns

        if ((pSat == 0.0_wp) .and. (TSIn == 0.0_wp)) then

            ! assigning Saturation temperature
            TSat = 0.0_wp

        else

            ! calculating initial estimate for temperature in the TSat procedure. I will also use this variable to
            ! iterate over the Newton's solver
            TSat = TSIn

            ! iteration counter
            ns = 0

            ! underrelaxation factor
            Om = 1.0e-3_wp
            do while ((abs(FT) > ptgalpha_eps) .or. (ns == 0))
                ! increasing counter
                ns = ns + 1

                ! calculating residual
                FT = TSat*((cvs(lp)*gs_min(lp) - cvs(vp)*gs_min(vp)) &
                           *(1 - log(TSat)) - (qvps(lp) - qvps(vp)) &
                           + cvs(lp)*(gs_min(lp) - 1)*log(pSat + ps_inf(lp)) &
                           - cvs(vp)*(gs_min(vp) - 1)*log(pSat + ps_inf(vp))) &
                     + qvs(lp) - qvs(vp)

                ! calculating the jacobian
                dFdT = &
                    -(cvs(lp)*gs_min(lp) - cvs(vp)*gs_min(vp))*log(TSat) &
                    - (qvps(lp) - qvps(vp)) &
                    + cvs(lp)*(gs_min(lp) - 1)*log(pSat + ps_inf(lp)) &
                    - cvs(vp)*(gs_min(vp) - 1)*log(pSat + ps_inf(vp))

                ! updating saturation temperature
                TSat = TSat - Om*FT/dFdT

            end do

        end if

    end subroutine s_TSat

    subroutine s_get_Ysat(m_i, TS, pS, m_g, Y_sat, valid_Ysat)

        real(wp), dimension(num_fluids), intent(in) :: m_i
        real(wp), intent(in) :: TS, pS
        real(wp), intent(in) :: m_g
        real(wp), intent(out) :: Y_sat
        logical, intent(out) :: valid_Ysat

        integer :: i

        ! real(wp) :: rho, dynE, rhoe
        real(wp) :: pSat, Plog10_pSat_mmHg
        real(wp) :: X_sat
        real(wp) :: Y_g(Dif_size)
        real(wp) :: W_nc

        ! mCP = 0.0_wp; mQ = 0.0_wp; alphaGam = 0.0_wp; alphaPi_inf = 0.0_wp
        ! ! Performing tests before initializing the pT-equilibrium
        ! !$acc loop seq
        ! do i = 1, num_fluids

        !     ! sum of the total alpha*rho*cp of the system
        !     mCP = mCP + q_cons_vf(i + contxb - 1)%sf(j, k, l)*cvs(i)*gs_min(i)

        !     ! sum of the total alpha*rho*q of the system
        !     mQ = mQ + q_cons_vf(i + contxb - 1)%sf(j, k, l)*qvs(i)

        !     alphaGam = alphaGam + q_cons_vf(i + advxb - 1)%sf(j, k, l)*gammas(i)

        !     alphaPi_inf = alphaPi_inf + q_cons_vf(i + advxb - 1)%sf(j, k, l)*pi_infs(i)

        ! end do

        ! ! rho = 0.0_wp; dynE = 0.0_wp
        ! ! !$acc loop seq
        ! ! do i = 1, num_fluids
        ! !     rho = rho + q_cons_vf(i + contxb - 1)%sf(j, k, l)
        ! ! end do

        ! ! !$acc loop seq
        ! ! do i = 1, num_fluids
        ! !     dynE = dynE + 5.0e-1_wp*q_cons_vf(i)%sf(j, k, l)**2/rho
        ! ! end do

        ! ! rhoe = q_cons_vf(E_idx)%sf(j, k, l) - dynE

        ! pS = ( rhoe - alphaPi_inf - mQ )/ alphaGam
        
        ! TS = (rhoe + pS - mQ)/mCP

        Plog10_pSat_mmHg = 8.07131_wp - 1730.63_wp/(233.426_wp + TS - 273.15_wp)
        pSat = (10.0_wp**Plog10_pSat_mmHg)*133.322_wp

        if (pSat >= pS) then
            valid_Ysat = .false.
            return
        else
            valid_Ysat = .true.
        end if

        X_sat = pSat / pS

        !$acc loop seq
        do i = 1, Dif_size
            Y_g(i) = m_i(Dif_idx(i)) / m_g
        end do

        W_nc = 0.0_wp
        !$acc loop seq
        do i = 2, Dif_size
            W_nc = W_nc + Y_g(i) / fluid_pp(Dif_idx(i))%W
        end do

        W_nc = (1.0_wp - Y_g(1)) / W_nc

        Y_sat = X_sat * fluid_pp(Dif_idx(1))%W / (X_sat * fluid_pp(Dif_idx(1))%W + (1.0_wp - X_sat) * W_nc)

    end subroutine s_get_Ysat

    subroutine s_compute_residual(m_i, alpha_i, dm, rhoe, R)

        real(wp), intent(in) :: m_i(num_fluids), alpha_i(num_fluids)
        real(wp), intent(in) :: dm, rhoe
        real(wp), intent(out) :: R

        real(wp) :: Y_g(Dif_size)
        real(wp) :: m_g
        real(wp) :: dm_eff
        real(wp) :: mCP, mQ, alphaGam, alphaPi_inf
        real(wp) :: pS, TS
        real(wp) :: pSat, Plog10_pSat_mmHg
        real(wp) :: X_sat
        real(wp) :: W_nc
        real(wp) :: W_g
        real(wp) :: Y_sat

        real(wp), dimension(num_fluids) :: m_trial, alpha_trial

        integer :: i

        dm_eff = max(dm, -m_i(vp))
        dm_eff = min(dm_eff, m_i(lp))

        m_trial(:) = m_i(:)
        m_trial(lp) = m_i(lp) - dm_eff
        m_trial(vp) = m_i(vp) + dm_eff

        alpha_trial(:) = alpha_i(:)

        ! 2. Given trial masses and fixed rhoe, solve p,T,alpha
        call s_solve_pT_alpha_from_masses(m_trial, rhoe, pS, TS, alpha_trial)

        m_g = 0.0_wp
        !$acc loop seq
        do i = 1, Dif_size
            m_g = m_g + m_trial(Dif_idx(i))
        end do
        
        !$acc loop seq
        do i = 1, Dif_size
            Y_g(i) = m_trial(Dif_idx(i)) / m_g
        end do

        mCP = 0.0_wp; mQ = 0.0_wp; alphaGam = 0.0_wp; alphaPi_inf = 0.0_wp
        !$acc loop seq
        do i = 1, num_fluids
            
            ! sum of the total alpha*rho*cp of the system
            mCP = mCP + m_trial(i)*cvs(i)*gs_min(i)

            ! sum of the total alpha*rho*q of the system
            mQ = mQ + m_trial(i)*qvs(i)

            alphaGam = alphaGam + alpha_trial(i)*gammas(i)

            alphaPi_inf = alphaPi_inf + alpha_trial(i)*pi_infs(i)

        end do

        !compute new Psat
        Plog10_pSat_mmHg = 8.07131_wp - 1730.63_wp/(233.426_wp + TS - 273.15_wp)
        pSat = (10.0_wp**Plog10_pSat_mmHg)*133.322_wp

        X_sat = pSat / pS

        W_nc = 0.0_wp
        !$acc loop seq
        do i = 2, Dif_size
            W_nc = W_nc + Y_g(i) / fluid_pp(Dif_idx(i))%W
        end do

        W_nc = (1.0_wp - Y_g(1)) / W_nc

        Y_sat = X_sat * fluid_pp(Dif_idx(1))%W / (X_sat * fluid_pp(Dif_idx(1))%W + (1.0_wp - X_sat) * W_nc)

        R = Y_g(1) - Y_sat

    end subroutine s_compute_residual

    subroutine s_solve_pT_alpha_from_masses(m_trial, rhoe, pS, TS, alpha_trial)

        real(wp), dimension(num_fluids), intent(in) :: m_trial
        real(wp), intent(in) :: rhoe
        real(wp), intent(out) :: pS, TS
        real(wp), dimension(num_fluids), intent(inout) :: alpha_trial

        real(wp) :: Fp, Fp_plus, Fp_minus, dFpdp

        real(wp) :: mCP, mQ, alphaGam, alphaPi_inf
        real(wp) :: rhok(num_fluids)
        real(wp) :: Y_g(Dif_size)
        real(wp) :: W_g, alpha_g
        real(wp) :: m_g
        real(wp) :: pS_trial, dp
        real(wp) :: tol_p
        integer :: iter, max_p_iter, i

        max_p_iter = 20
        tol_p = 1.0e-6_wp

        mCP = 0.0_wp
        mQ  = 0.0_wp
        alphaGam = 0.0_wp
        alphaPi_inf = 0.0_wp

        do i = 1, num_fluids
            mCP = mCP + m_trial(i)*cvs(i)*gs_min(i)
            mQ  = mQ  + m_trial(i)*qvs(i)
            alphaGam = alphaGam + alpha_trial(i)*gammas(i)
            alphaPi_inf = alphaPi_inf + alpha_trial(i)*pi_infs(i)
        end do

        !initial guess for pS
        pS_trial = ( rhoe - alphaPi_inf - mQ )/ alphaGam

        do iter = 1, max_p_iter

            call s_pS_residual(pS_trial, rhoe, mQ, mCP, m_trial, Fp)

            if (abs(Fp)/max(abs(pS_trial), 1.0e5_wp) < tol_p) exit

            dp = 1.0e-6_wp*max(abs(pS_trial), 1.0e5_wp)

            call s_pS_residual(pS_trial + dp, rhoe, mQ, mCP, m_trial, Fp_plus)

            call s_pS_residual(pS_trial - dp, rhoe, mQ, mCP, m_trial, Fp_minus)

            dFpdp = (Fp_plus - Fp_minus)/(2.0_wp*dp)

            pS_trial = pS_trial - Fp/dFpdp

        end do

        pS = pS_trial

        TS = (rhoe + pS_trial - mQ)/mCP

        rhok(lp) = (pS_trial + ps_inf(lp)) &
                /((gs_min(lp) - 1.0_wp)*cvs(lp)*TS)

        alpha_trial(lp) = m_trial(lp) / rhok(lp)

        alpha_g = 1 - alpha_trial(lp)

        m_g = 0.0_wp
        !$acc loop seq
        do i = 1, Dif_size
            m_g = m_g + m_trial(Dif_idx(i))
        end do

        !$acc loop seq
        do i = 1, Dif_size
            Y_g(i) = m_trial(Dif_idx(i)) / m_g
        end do

        W_g = 0.0_wp
        !$acc loop seq
        do i = 1, Dif_size
            W_g = W_g + Y_g(i) / fluid_pp(Dif_idx(i))%W
        end do

        W_g = 1.0_wp / W_g

        !$acc loop seq
        do i = 1, Dif_size
            alpha_trial(Dif_idx(i)) = Y_g(i)*alpha_g*W_g / fluid_pp(Dif_idx(i))%W
        end do

    end subroutine s_solve_pT_alpha_from_masses


    subroutine s_pS_residual(pS, rhoe, mQ, mCP, m_i, Fp)

        real(wp), intent(in) :: pS
        real(wp), intent(in) :: rhoe
        real(wp), intent(in) :: mQ, mCP
        real(wp), dimension(num_fluids), intent(in) :: m_i
        real(wp), intent(out) :: Fp

        real(wp) :: TS
        real(wp) :: W_g, m_g, alpha_g
        real(wp) :: Y_g(Dif_size)
        real(wp) :: rhok(num_fluids), alpha_i(num_fluids)
        real(wp) :: alphaGam, alphaPi_inf
        real(wp) :: pS_update

        integer :: i



        TS = (rhoe + pS - mQ)/mCP

        rhok(lp) = (pS + ps_inf(lp)) &
                    /((gs_min(lp) - 1.0_wp)*cvs(lp)*TS)

        alpha_i(lp) = m_i(lp) / rhok(lp)

        alpha_g = 1 - alpha_i(lp)

        m_g = 0.0_wp
        !$acc loop seq
        do i = 1, Dif_size
            m_g = m_g + m_i(Dif_idx(i))
        end do

        !$acc loop seq
        do i = 1, Dif_size
            Y_g(i) = m_i(Dif_idx(i)) / m_g
        end do

        W_g = 0.0_wp
        !$acc loop seq
        do i = 1, Dif_size
            W_g = W_g + Y_g(i) / fluid_pp(Dif_idx(i))%W
        end do

        W_g = 1.0_wp / W_g

        !$acc loop seq
        do i = 1, Dif_size
            alpha_i(Dif_idx(i)) = Y_g(i)*alpha_g*W_g / fluid_pp(Dif_idx(i))%W
        end do
        
        alphaGam = 0.0_wp
        alphaPi_inf = 0.0_wp
        !$acc loop seq
        do i = 1, num_fluids
            alphaGam = alphaGam + alpha_i(i)*gammas(i)
            alphaPi_inf = alphaPi_inf + alpha_i(i)*pi_infs(i)
        end do

        pS_update = ( rhoe - alphaPi_inf - mQ )/ alphaGam

        Fp = pS_update - pS

    end subroutine s_pS_residual

    !>  This subroutine finalizes the phase change module
    subroutine s_finalize_relaxation_solver_module
    end subroutine s_finalize_relaxation_solver_module

#endif

end module m_phase_change
