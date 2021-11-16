module bucket
  use xtb_mctc_accuracy, only: wp
  implicit none

  contains

    subroutine check_bucket(index_mass, isotope_masses, exact_intensity, &
        list_masses,  intensity,  count_mass, chrg)
   
      integer :: loop, loop2, index_mass
      integer :: count_mass
      integer :: sum_index, i
      integer :: check
    
      !real(wp) :: intensity(index_mass)
      real(wp) :: xmass
      real(wp) :: isotope_masses(index_mass)
      real(wp) :: exact_intensity(index_mass)
      real(wp) :: maxi
      real(wp) :: chrg

      real(wp) :: list_masses(10000)
      real(wp) :: intensity(10000)
    
      !real(wp), allocatable :: list_masses(:)
      !real(wp), allocatable :: save_list(:)
      !real(wp), allocatable :: save_int(:)
      !real(wp), allocatable :: intensity(:)
      !real (wp) :: save_list(1000)
      !real (wp) :: save_list(1000)

      logical :: there = .true.
      logical :: last  = .false.

      xmass = 0
      loop = 0
      loop2 = 0
      check = 0
      !count_mass = count_mass + index_mass
      if(count_mass == 0)then
        sum_index =  index_mass
      !  allocate ( list_masses(index_mass))
      !  allocate ( intensity(index_mass))
      else
        sum_index = count_mass + index_mass
      !  if(allocated(save_list)) write(*,*) 'JOP'
      !  if(allocated(save_list))deallocate (save_list)
      !  if(allocated(save_int))deallocate (save_int)
      !  allocate ( save_list(count_mass))
      !  allocate ( save_int (count_mass))

      ! !write(*,*) 'COUNT MASS', count_mass
      ! do i = 1, count_mass 
      !  save_list(i) = list_masses(i)
      !  save_int(i)  = intensity(i)
      ! enddo

      !  deallocate (list_masses)
      !  deallocate (intensity)
      !  allocate ( list_masses(sum_index))
      !  allocate ( intensity(sum_index))

      !   !write(*,*) 'save', save_list
      !   !write(*,*) 'list', list_masses

      !  do i = 1, count_mass
      !    list_masses(i)  = save_list(i)
      !    intensity(i)    = save_int(i)
      !  enddo
      endif


      !if (.not. allocated(list_masses)) allocate ( list_masses(sum_index))

      !write(*,*) 'BUCKET'
      !write(*,*) 'index', index_mass
      !write(*,*) 'sum', sum_index
      !write(*,*) 'count', count_mass

outer: do loop2 = 1, index_mass 
        loop = 0
        !write(*,*) isotope_masses(loop2)

        !write(*,*) 'LOOP2', loop2
        do 
          loop = loop + 1

          if  ( list_masses(loop) == isotope_masses(loop2) ) then
          !write(*,*) 'LOOP', loop, 'TRUE'
            there = .true.
            !write(*,*) 'TRUE'
            !write(*,*) 'intensity', loop,intensity(loop)
            !intensity(loop) = intensity(loop) + 1
            intensity(loop) = intensity(loop) +  1 * exact_intensity(loop2) &
              * abs(chrg)
            if (loop2 < index_mass ) cycle outer 
            if (loop2 == index_mass ) exit 

          elseif ( list_masses(loop) == 0.0_wp ) then
            there = .false.
             !write(*,*) 'NULL'
            exit

          !>> false if not in list, store
          elseif ( list_masses(loop) /= isotope_masses(loop2) ) then
          !write(*,*) 'LOOP', loop, 'FALSE'
            there = .false.
            !write(*,*) 'FALSE'
            !write(*,*) loop, loop2
            !if (loop == count_mass ) exit
            if (loop == sum_index ) exit
          endif
        enddo


        if ( .not. there ) then
          count_mass = count_mass + 1
          list_masses(count_mass) = isotope_masses(loop2)
          intensity(count_mass) = intensity(count_mass) + 1 * exact_intensity(loop2) &
            * abs(chrg)
        endif

      enddo outer

      !if(allocated(intensity) .and. allocated(save_int)) then
      !  deallocate (save_int)
      !  allocate (save_int(count_mass))

      !  do i = 1, count_mass 
      !    save_int(i)=intensity(i)
      !   enddo

      !    deallocate (intensity)
      !    allocate ( intensity(count_mass))

      !    do i = 1, count_mass 
      !      intensity(i)    = save_int(i)
      !    enddo
      !endif

      !if(allocated(list_masses) .and. allocated(save_list)) then
      !  deallocate (save_list)
      !  allocate (save_list(count_mass))

      !  do i = 1, count_mass 
      !    save_list(i)=list_masses(i)
      !   enddo

      !    deallocate (list_masses)
      !    allocate ( list_masses(count_mass))

      !    do i = 1, count_mass 
      !      list_masses(i)    = save_list(i)
      !    enddo
      !  write(*,*) 
      !endif
      !deallocate (save_list)
      !deallocate (save_int)

      !    do i = 1, count_mass 
        !    write(*,*) i,list_masses(i),intensity(i)
            !write(*,*) i,intensity(i)
       !   enddo
   
    end subroutine check_bucket


end module bucket
