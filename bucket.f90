module bucket
  use xtb_mctc_accuracy, only: wp
  implicit none

  contains

    subroutine check_bucket(xmass, list_masses, index_mass, counter)
    
      integer :: loop, index_mass, counter
    
      real(wp) :: xmass, list_masses(counter)
    
      logical :: there = .true.
    
      if ( xmass > 0 ) then 
        do loop = 1, counter
          if     ( list_masses(loop) == xmass ) then
            there = .true.
            exit
          elseif ( list_masses(loop) /= xmass ) then
            there = .false.
          endif
        enddo
        if ( .not. there ) then
          index_mass = index_mass + 1
          list_masses(index_mass) = xmass
        endif
      endif
   
    end subroutine check_bucket


    !subroutine count_bucket (counter, save_mass, list_masses, index_mass)

    !  integer :: loop, counter, index_mass

    !  real(wp) :: save_mass(counter), list_masses(index_mass)

    !  do outer = 1, counter

    !    do inner = 1, index_mass

    !      if (save_mass(outer) == list_masses(inner)) then
    !        intensity(inner) = intensity(inner) + 1
    !      endif

    !    enddo

    !  enddo

    !end subroutine count_bucket 

end module bucket
