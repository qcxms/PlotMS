!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Generate random normal derivate by Box-Muller distribution

! G.E.P. Box and M.E. Muller
! Ann. Meth. Statist. 29(2): 610-611, 1958
! DOI: https://doi.org/10.1214/aoms/1177706645
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module qcxms_boxmuller
  use xtb_mctc_accuracy
  implicit none

  contains 

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Do variation of the starting collision-energies
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function vary_energies(E_in, E_Distr) result(E_Scale)

      real(wp) :: E_in,sigma,E_Array(2)
      real(wp) :: dum,dum2 !dummys
      real(wp) :: E_Scale ! End-result
      real(wp) :: E_Distr ! How strongly the energy is distributed
      real(wp) :: z0,z1

      real(wp), parameter :: pi = 3.1415926535897932384626433832795029_wp
      
      call random_seed()
      call random_number(dum)
      call random_number(dum2)

      sigma = E_in * dble(E_Distr)  !dble(0.25)

      z0 = sqrt(-2.0_wp * LOG(dum)) * cos(2.0_wp * pi * dum2)
      z1 = sqrt(-2.0_wp * LOG(dum)) * sin(2.0_wp * pi * dum2)

      E_Array(1) = z0 * sigma + E_in
      E_Array(2) = z1 * sigma + E_in

      if (dum > 0.5) then
        E_Scale = E_Array(1)
      else
        E_Scale = E_Array(2)
      endif

    end function vary_energies

end module qcxms_boxmuller
