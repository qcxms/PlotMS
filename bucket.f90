module bucket
  use xtb_mctc_accuracy, only: wp
  implicit none

  contains

    subroutine check_bucket(index_mass, isotope_masses, exact_intensity, &
        list_masses,  intensity, count_mass)
    
      integer :: loop, loop2, index_mass
      integer :: count_mass, entry
    
      real(wp) :: intensity(index_mass)
      real(wp) :: xmass, list_masses(index_mass)
      real(wp) :: isotope_masses(index_mass)
      real(wp) :: exact_intensity(index_mass)
    
      logical :: there = .true.

      xmass = 0
      entry = count_mass
      count_mass = count_mass + index_mass

      write(*,*) 'INDEX MASS', index_mass

outer: do loop = 1, index_mass
        !do loop2 = 1, count_mass
          write(*,*) isotope_masses(loop), exact_intensity(loop)
          write(*,*) 

          write(*,*) loop
         ! write(*,*) loop2
          xmass = isotope_masses(loop)
          write(*,*) 'INDEX MASS', xmass

          !if ( xmass > mzmin ) then 
          !> loop over all entries in the list
         ! do loop = 1, counter
          !>> true if already in list, end
          do loop2 = 1, count_mass
            if     ( list_masses(loop2) == xmass ) then
              there = .true.
              write(*,*) 'TRUE'
            !>> false if not in list, store
            elseif ( list_masses(loop2) /= xmass ) then
              there = .false.
              write(*,*) 'FALSE'
              entry = loop2
              exit 
            endif
          enddo
            !> if it is not in the list, add it
            if ( .not. there ) then
            !  index_mass = index_mass + 1
              write(*,*) 'entry', entry
              list_masses(entry) = xmass
              write(*,*) list_masses
            endif
        !enddo outer

          !> if it is not in the list, add it
          !if ( .not. there ) then
          !!  index_mass = index_mass + 1
          !  write(*,*) 'LP2', loop
          !  list_masses(loop) = xmass
          !  write(*,*) list_masses
          !endif
        !endif
        !enddo

        !> count the number of signals
        do loop2 = 1, count_mass

          if (list_masses(loop2) == xmass) then
            intensity(loop2) = (intensity(loop2) + 1) * exact_intensity(loop)
            write(*,*) 'INTENSITY', isotope_masses(loop), intensity(loop)
          endif

        enddo
        if ( .not. there ) then
          count_mass = count_mass + 1
          there = .true.
        endif
        write(*,*) 
      enddo outer
   
    end subroutine check_bucket


end module bucket
