module bucket
  use xtb_mctc_accuracy, only: wp
  implicit none

  contains

    subroutine check_bucket(counter, xmass, list_masses, index_mass, intensity, mzmin)
    
      integer :: loop, index_mass, counter
      integer :: intensity(counter), mzmin
    
      real(wp) :: xmass, list_masses(counter)
    
      logical :: there = .true.
    
      if ( xmass > mzmin ) then 
        !> loop over all entries in the list
        do loop = 1, counter
          !>> true if already in list, end
          if     ( list_masses(loop) == xmass ) then
            there = .true.
            exit
          !>> false if not in list, store
          elseif ( list_masses(loop) /= xmass ) then
            there = .false.
          endif
        enddo

        !> if it is not in the list, add it
        if ( .not. there ) then
          index_mass = index_mass + 1
          list_masses(index_mass) = xmass
        endif
      endif

      !> count the number of signals
      do loop = 1, index_mass

        if (list_masses(loop) == xmass) then
          intensity(loop) = intensity(loop) + 1
        endif

      enddo
   
    end subroutine check_bucket


end module bucket
