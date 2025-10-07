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

    use m_helper              !< Helper functions

    ! ==========================================================================

    implicit none

    private; public :: s_compute_sum_alpha_g

contains

    subroutine s_compute_sum_alpha_g(q_cons_vf, bounds)

            ! From the volume fractions of each mixture gas component, compute the
            ! total gas volume fraction field.

            type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf
            type(int_bounds_info), dimension(1:3), intent(in) :: bounds

            integer :: x, y, z, i
            real(wp) :: sum_alpha_g
    

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

end module m_diffusion