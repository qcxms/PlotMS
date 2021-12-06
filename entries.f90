module count_entries
  use xtb_mctc_accuracy, only: wp
  implicit none

  contains

    subroutine check_entries(index_mass, isotope_masses, exact_intensity, &
        list_masses,  intensity,  count_mass, chrg)
   
      integer :: loop, loop2, index_mass
      integer :: count_mass
      integer :: sum_index, i
      integer :: check
    
      !real(wp) :: intensity(index_mass)
      real(wp) :: xmass
      real(wp) :: isotope_masses(index_mass)
      real(wp) :: exact_intensity(index_mass)
      real(wp) :: chrg

      real(wp) :: list_masses(10000)
      real(wp) :: intensity(10000)

      real(wp) :: mass_diff
    
      !real(wp), allocatable :: list_masses(:)
      !real(wp), allocatable :: save_list(:)
      !real(wp), allocatable :: save_int(:)
      !real(wp), allocatable :: intensity(:)
      !real (wp) :: save_list(1000)
      !real (wp) :: save_list(1000)

      logical :: there = .true.

      xmass = 0
      loop = 0
      loop2 = 0
      check = 0
      !count_mass = count_mass + index_mass
      if(count_mass == 0)then
        sum_index =  index_mass
      else
        sum_index = count_mass + index_mass
      endif


      !> loop over the list from the isotope subroutine
outer: do loop2 = 1, index_mass 
        loop = 0

        !> loop over all entries of new (check) list 
inner:  do 
          loop = loop + 1
          mass_diff = abs(list_masses(loop) - isotope_masses(loop2))

          if ( list_masses(loop) == 0.0_wp ) then
            there = .false.
             !write(*,*) 'NULL'
            exit inner

          !elseif  ( list_masses(loop) == isotope_masses(loop2) ) then
          elseif  ( mass_diff < 1.0d0-10 .or. mass_diff == 0.0_wp) then
            there = .true.
            intensity(loop) = intensity(loop) +  1 * exact_intensity(loop2) &
              * abs(chrg)
            if (loop2 < index_mass ) cycle outer 
            if (loop2 == index_mass ) exit inner


          !>> false if not in list, store
          !elseif ( list_masses(loop) /= isotope_masses(loop2) ) then
          elseif  ( mass_diff > 1.0d0-10 ) then
            there = .false.
            if (loop == sum_index ) exit inner
          endif
        enddo inner


        if ( .not. there ) then
          count_mass = count_mass + 1
          list_masses(count_mass) = isotope_masses(loop2)
          intensity(count_mass) = intensity(count_mass) + 1 * exact_intensity(loop2) &
            * abs(chrg)
        endif

      enddo outer

   
    end subroutine check_entries
    
end module count_entries 
