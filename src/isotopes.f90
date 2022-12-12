!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  ! compute the isotopic mass distribution for
  ! given chemical formula (ntot atoms with OZ iat_save())
  ! returns nsig masses in mass() with probability
  ! this is a quick and dirty MC algorithm and its run-time depends
  ! critically on the number of trials nrnd
  
  ! TK See: ATOMIC WEIGHTS OF THE ELEMENTS: REVIEW 2000 (IUPAC Technical Report)
  ! TK Pure Appl. Chem., Vol. 75, No. 6, pp. 683800, 2003.
 
module isotope_pattern
  use xtb_mctc_accuracy, only: wp
  implicit none

  contains

subroutine isotope(counter, mzmin, ntot, iat_save, maxatm, rnd, &
    no_isotopes, index_mass, exact_intensity, isotope_masses, z_chrg)

  integer :: ntot,iat_save(*),maxatm
  integer :: niso(200)
  integer :: n,i,j,iso,iti
  integer :: counter
  integer,parameter :: nrnd = 50000 
  integer :: loop, index_mass
!  integer :: tmp_intensity
  integer :: store_int(1000)
  integer :: mzmin
  integer :: z_chrg

  real(wp) :: rnd(nrnd,maxatm)
  real(wp) :: r,sum_prob
  real(wp) :: prob(200,10),massiso(200,10),p1,p2,x,xmass
  real(wp) :: list_masses(nrnd)
  real(wp),allocatable :: isotope_masses(:)

  real(wp) :: mipmax, iipmax ! mass and intensity of highest isotope peak
  integer :: indexipmax ! index of peak with highest intensity
  real(wp) :: current_mass

  real(wp), allocatable :: exact_intensity(:)
  !real(wp) :: exact_intensity(1000)

  logical  :: no_isotopes
  logical  :: there

  niso=0
  prob=0
  massiso=0

  
   !  1 H  (Hydrogen)
        niso(1)           = 2
        prob(1,1)         = 0.0115_wp
        prob(1,2)         = 99.9885_wp
        massiso(1,1)      = 2.014101_wp
        massiso(1,2)      = 1.007825_wp

   !  5 B  (Boron)
        niso(5)           = 2
        prob(5,1)         = 19.9_wp
        prob(5,2)         = 80.1_wp
        massiso(5,1)      = 10.0129370_wp
        massiso(5,2)      = 11.0093055_wp

   !  6 C  (Carbon)
        niso(6)           = 2
        prob(6,1)         = 1.07_wp
        prob(6,2)         = 98.93_wp
        massiso(6,1)      = 13.003354_wp
        massiso(6,2)      = 12.00_wp
  
   !  7 N  (Nitrogen)
        niso(7)           = 2
        prob(7,1)         = 0.368_wp
        prob(7,2)         = 99.632_wp
        massiso(7,1)      = 15.000108_wp
        massiso(7,2)      = 14.003074_wp
  
   !  8 O  (Oxygen)
        niso(8)           = 3
        prob(8,1)         = 0.038_wp
        prob(8,2)         = 0.205_wp
        prob(8,3)         = 99.757_wp
        massiso(8,1)      = 16.999131_wp
        massiso(8,2)      = 17.999159_wp
        massiso(8,3)      = 15.994914_wp
  
   !  9 F  (Fluorine)
        niso(9)           = 1
        prob(9,1)         = 100.0_wp
        massiso(9,1)      = 18.998403_wp
  
  !  13 Al (Aluminium)                                                                                 
          niso(13)        = 1
          prob(13,1)      = 100.0_wp
          massiso(13,1)   = 26.981538_wp
  
   ! 14 Si (Silicon)
        niso(14)          = 3
        prob(14,1)        = 92.223_wp
        prob(14,2)        = 4.685_wp
        prob(14,3)        = 3.092_wp
        massiso(14,1)     = 27.976926_wp
        massiso(14,2)     = 28.976494_wp
        massiso(14,3)     = 29.973770_wp
  
   ! 15 P  (Phosphorus)
        niso(15)          = 1
        prob(15,1)        = 100.0_wp
        massiso(15,1)     = 30.973761_wp
  
   ! 16 S  (Sulfur)
        niso(16)          = 4
        prob(16,1)        = 0.02_wp
        prob(16,2)        = 0.76_wp
        prob(16,3)        = 4.29_wp
        prob(16,4)        = 94.93_wp
        massiso(16,1)     = 35.967080_wp
        massiso(16,2)     = 32.971458_wp
        massiso(16,3)     = 33.967867_wp
        massiso(16,4)     = 31.972071_wp
  
   ! 17 Cl (Chlorine) 
        niso(17)          = 2
        prob(17,1)        = 75.76_wp
        prob(17,2)        = 24.24_wp
        massiso(17,1)     = 34.968852_wp
        massiso(17,2)     = 36.965902_wp

   ! 20 Ca (Calcium)
        niso(20)          = 5
        prob(20,1)        = 96.94_wp
        prob(20,2)        = 0.65_wp
        prob(20,3)        = 0.14_wp
        prob(20,4)        = 2.09_wp
        prob(20,5)        = 0.18_wp
        massiso(20,1)     = 39.962591_wp
        massiso(20,2)     = 41.958618_wp
        massiso(20,3)     = 42.958766_wp
        massiso(20,4)     = 43.955482_wp
        massiso(20,5)     = 45.953688_wp      
  
   ! 22 Ti (Titanium) 
        niso(22)          = 5           
        prob(22,1)        = 8.0_wp
        prob(22,2)        = 7.3_wp
        prob(22,3)        = 73.8_wp
        prob(22,4)        = 5.5_wp
        prob(22,5)        = 5.4_wp
        massiso(22,1)     = 45.952632_wp
        massiso(22,2)     = 46.951763_wp
        massiso(22,3)     = 47.947946_wp
        massiso(22,4)     = 48.947870_wp
        massiso(22,5)     = 49.944791_wp
                               
   ! 24 Cr (Chromium)                
        niso(24)          = 4           
        prob(24,1)        = 4.35_wp
        prob(24,2)        = 83.79_wp
        prob(24,3)        = 9.50_wp
        prob(24,4)        = 2.37_wp
        massiso(24,1)     = 49.946044_wp
        massiso(24,2)     = 51.940507_wp
        massiso(24,3)     = 52.940649_wp
        massiso(24,4)     = 53.938880_wp
                               
   ! 25 Mn (Manganese)
        niso(25)          = 1           
        prob(25,1)        = 100.0_wp
        massiso(25,1)     = 54.938045_wp
  
  !  26 Fe (Iron)
        niso(26)          = 4
        prob(26,1)        = 5.845_wp
        prob(26,2)        = 91.754_wp
        prob(26,3)        = 2.119_wp
        prob(26,4)        = 0.2_wp
        massiso(26,1)     = 53.939_wp
        massiso(26,2)     = 55.934_wp
        massiso(26,3)     = 56.935_wp
        massiso(26,4)     = 57.933_wp
  
   ! 27 Co (Cobalt)
        niso(27)          = 1         
        prob(27,1)        = 100.000_wp
        massiso(27,1)     = 58.933195_wp
                             
  !  28 Ni (Nickel)
        niso(28)          = 5         
        prob(28,1)        = 68.08_wp
        prob(28,2)        = 26.22_wp
        prob(28,3)        = 1.14_wp
        prob(28,4)        = 3.63_wp
        prob(28,5)        = 0.93_wp
        massiso(28,1)     = 57.935343_wp
        massiso(28,2)     = 59.930786_wp
        massiso(28,3)     = 60.931056_wp
        massiso(28,4)     = 61.928345_wp
        massiso(28,5)     = 63.927966_wp
                             
   ! 29 Cu (Copper)                
        niso(29)          = 2         
        prob(29,1)        = 69.15_wp
        prob(29,2)        = 30.85_wp
        massiso(29,1)     = 62.929597_wp
        massiso(29,2)     = 64.927789_wp
                             
   ! 30 Zn (Zinc)
        niso(30)          = 5         
        prob(30,1)        = 48.6_wp
        prob(30,2)        = 27.9_wp
        prob(30,3)        = 4.1_wp
        prob(30,4)        = 18.8_wp
        prob(30,5)        = 0.6_wp
        massiso(30,1)     = 63.929142_wp
        massiso(30,2)     = 65.926033_wp
        massiso(30,3)     = 66.927127_wp
        massiso(30,4)     = 67.924884_wp
        massiso(30,5)     = 69.925319_wp
  
   ! 32 Ge (Germanium)
        niso(32)          = 5          
        prob(32,1)        = 20.38_wp
        prob(32,2)        = 27.31_wp
        prob(32,3)        = 7.76_wp
        prob(32,4)        = 36.72_wp
        prob(32,5)        = 7.83_wp
        massiso(32,1)     = 69.924247_wp
        massiso(32,2)     = 71.922076_wp
        massiso(32,3)     = 72.923459_wp
        massiso(32,4)     = 73.921178_wp
        massiso(32,5)     = 75.921402_wp
                                 
   ! 33 As (Arsenic)
        niso(33)          = 1          
        prob(33,1)        = 100.0_wp
        massiso(33,1)     = 74.921596_wp
                                 
   ! 34 Se (Selenium)
        niso(34)          = 6          
        prob(34,1)        = 0.89_wp
        prob(34,2)        = 9.37_wp
        prob(34,3)        = 7.63_wp
        prob(34,4)        = 23.77_wp
        prob(34,5)        = 49.61_wp
        prob(34,6)        = 8.73_wp
        massiso(34,1)     = 73.922476_wp
        massiso(34,2)     = 75.919213_wp
        massiso(34,3)     = 76.919914_wp
        massiso(34,4)     = 77.917309_wp
        massiso(34,5)     = 79.916521_wp
        massiso(34,6)     = 81.916699_wp
   
   ! 35 Br (Bromine)
        niso(35)          = 2
        prob(35,1)        = 50.69_wp
        prob(35,2)        = 49.31_wp
        massiso(35,1)     = 78.91833_wp
        massiso(35,2)     = 80.91629_wp
  
   ! 40 Zr (Zirconium)
          niso(40)        = 5          
          prob(40,1)      = 51.45_wp
          prob(40,2)      = 11.22_wp
          prob(40,3)      = 17.15_wp
          prob(40,4)      = 17.38_wp
          prob(40,5)      = 2.80_wp
          massiso(40,1)   = 89.904704_wp
          massiso(40,2)   = 90.905646_wp
          massiso(40,3)   = 91.905040_wp
          massiso(40,4)   = 93.906315_wp
          massiso(40,5)   = 95.908273_wp
                                 
   ! 44 Ru (Ruthenium)
          niso(44)        = 7          
          prob(44,1)      = 5.554_wp
          prob(44,2)      = 1.873_wp
          prob(44,3)      = 12.761_wp
          prob(44,4)      = 12.607_wp
          prob(44,5)      = 17.062_wp
          prob(44,6)      = 31.551_wp
          prob(44,7)      = 18.623_wp
          massiso(44,1)   = 95.907598_wp
          massiso(44,2)   = 97.905287_wp
          massiso(44,3)   = 98.905939_wp
          massiso(44,4)   = 99.904219_wp
          massiso(44,5)   = 100.905582_wp
          massiso(44,6)   = 101.906323_wp
          massiso(44,7)   = 103.905433_wp
                                 
   ! 46 Pd (Palladium)
          niso(46)        = 6          
          prob(46,1)      = 1.02_wp
          prob(46,2)      = 11.15_wp
          prob(46,3)      = 22.34_wp
          prob(46,4)      = 27.33_wp
          prob(46,5)      = 26.47_wp
          prob(46,6)      = 11.73_wp
          massiso(46,1)   = 101.905609_wp
          massiso(46,2)   = 103.904036_wp
          massiso(46,3)   = 104.905085_wp
          massiso(46,4)   = 105.903486_wp
          massiso(46,5)   = 107.903892_wp
          massiso(46,6)   = 109.905153_wp
  
   ! 50 Sn (Tin) 
          niso(50)        = 10         
          prob(50,1)      = 0.97_wp
          prob(50,2)      = 0.66_wp
          prob(50,3)      = 0.34_wp
          prob(50,4)      = 14.54_wp
          prob(50,5)      = 7.68_wp
          prob(50,6)      = 24.22_wp
          prob(50,7)      = 8.59_wp
          prob(50,8)      = 32.58_wp
          prob(50,9)      = 4.63_wp
          prob(50,10)     = 5.79_wp
          massiso(50,1)   = 111.904848_wp
          massiso(50,2)   = 113.902779_wp
          massiso(50,3)   = 114.903342_wp
          massiso(50,4)   = 115.901741_wp
          massiso(50,5)   = 116.902952_wp
          massiso(50,6)   = 117.901603_wp
          massiso(50,7)   = 118.903308_wp
          massiso(50,8)   = 119.902194_wp
          massiso(50,9)   = 121.903439_wp
          massiso(50,10)  = 123.905274_wp
                               
   ! 51 Sb (Antimony)
          niso(51)        = 2          
          prob(51,1)      = 57.21_wp
          prob(51,2)      = 42.79_wp
          massiso(51,1)   = 120.903816_wp
          massiso(51,2)   = 122.904214_wp
                               
   ! 52 Te (Tellurium)
          niso(52)        = 8          
          prob(52,1)      = 0.09_wp
          prob(52,2)      = 2.55_wp
          prob(52,3)      = 0.89_wp
          prob(52,4)      = 4.74_wp
          prob(52,5)      = 7.07_wp
          prob(52,6)      = 18.84_wp
          prob(52,7)      = 31.75_wp
          prob(52,8)      = 34.09_wp
          massiso(52,1)   = 119.90402_wp
          massiso(52,2)   = 121.903044_wp
          massiso(52,3)   = 122.904270_wp
          massiso(52,4)   = 123.902817_wp
          massiso(52,5)   = 124.904431_wp
          massiso(52,6)   = 125.903312_wp
          massiso(52,7)   = 127.904463_wp
          massiso(52,8)   = 129.906224_wp
  
   ! 53  I (Iodine)
          niso(53)        = 1          
          prob(53,1)      = 100.00_wp
          massiso(53,1)   = 126.904473_wp

   ! 57 La (Lanthanum)
          niso(57)        = 2          
          prob(57,1)      = 0.09_wp
          prob(57,2)      = 99.91_wp
          massiso(57,1)   = 137.907112_wp
          massiso(57,2)   = 138.909477_wp

   ! 74 W (Tungsten)
          niso(74)        = 5          
          prob(74,1)      = 0.12_wp
          prob(74,2)      = 26.50_wp
          prob(74,3)      = 14.31_wp
          prob(74,4)      = 30.64_wp
          prob(74,5)      = 28.43_wp
          massiso(74,1)   = 179.946704_wp
          massiso(74,2)   = 181.948204_wp
          massiso(74,3)   = 182.950223_wp
          massiso(74,4)   = 183.950931_wp
          massiso(74,5)   = 185.954364_wp
                               
   ! 78 Pt (Platinum)
          niso(78)        = 5          
          prob(78,1)      = 0.782_wp
          prob(78,2)      = 32.97_wp
          prob(78,3)      = 33.83_wp
          prob(78,4)      = 25.24_wp
          prob(78,5)      = 7.16_wp
          massiso(78,1)   = 191.961038_wp
          massiso(78,2)   = 193.962680_wp
          massiso(78,3)   = 194.964791_wp
          massiso(78,4)   = 195.964951_wp
          massiso(78,5)   = 197.967893_wp
                               
   ! 80 Hg (Mercury)
          niso(80)        = 7          
          prob(80,1)      = 0.15_wp
          prob(80,2)      = 10.0_wp
          prob(80,3)      = 16.87_wp
          prob(80,4)      = 23.10_wp
          prob(80,5)      = 13.19_wp
          prob(80,6)      = 29.86_wp
          prob(80,7)      = 6.87_wp
          massiso(80,1)   = 195.965833_wp
          massiso(80,2)   = 197.966769_wp
          massiso(80,3)   = 198.968279_wp
          massiso(80,4)   = 199.968326_wp
          massiso(80,5)   = 200.970302_wp
          massiso(80,6)   = 201.970643_wp
          massiso(80,7)   = 203.973494_wp
                               
   ! 81 Tl (Thallium)
          niso(81)        = 2          
          prob(81,1)      = 29.52_wp
          prob(81,2)      = 70.48_wp
          massiso(81,1)   = 202.972344_wp
          massiso(81,2)   = 204.974427_wp
  
   ! 82 Pb (Lead)
          niso(82)        = 4          
          prob(82,1)      = 1.3_wp
          prob(82,2)      = 24.1_wp
          prob(82,3)      = 22.1_wp
          prob(82,4)      = 52.4_wp
          massiso(82,1)   = 203.973043_wp
          massiso(82,2)   = 205.974465_wp
          massiso(82,3)   = 206.975897_wp
          massiso(82,4)   = 207.976652_wp
                               
   ! 83 Bi (Bismuth)
          niso(83)        = 1          
          prob(83,1)      = 100.0_wp
          massiso(83,1)   = 208.980398_wp
           

  ! postprocessing starts here 
  
  prob = prob * 0.01_wp
  list_masses = 0.0_wp
  index_mass = 0
  loop = 0
  store_int = 0.0_wp

  !  Has to be done only once to check for correct isotope entries
  !> loop over all mases to element 83 (Bismuth)
  if (counter == 1) then
    do i = 1, 83
       sum_prob = 0
       do j = 1, niso(i)
          sum_prob = sum_prob + prob(i,j)
       enddo
       if ( niso(i) > 0 .and. abs(sum_prob - 1.0_wp) > 0.01_wp ) stop 'internal isotope error 1'
    enddo
  endif

!> this checks if the fragment has atms larger than 100 (no clue why this is checked)
  !write(*,*) 'ntot', ntot
  !do i = 1, ntot
  !  write(*,*) i, iat_save(i)
  !  if ( iat_save(i) > 100 ) then
  !    niso   (iat_save(i))   = 1
  !    prob   (iat_save(i),1) = 1.0_wp
  !    massiso(iat_save(i),1) = iat_save(i) - 100.0_wp
  !  endif
  !enddo
 
  !> if number of isotopes == 0 <- often errors if wrong .res file is read (is (hopefully) fixed)
  do i=1,ntot
    if(niso(iat_save(i)) == 0) stop 'internal isotope error 2'
  enddo


  
  ! Calculate isotope pattern                     

    ! loop over random number runs
    do n = 1, nrnd
      xmass = 0
      do i = 1, ntot
         iti = iat_save(i)
         r = rnd(n,i)
         p1 = 0.0_wp
         p2 = prob(iti,1)
         do iso = 1, niso(iti)
            if ( r >= p1 .and. r <= p2 ) then
               x = massiso(iti,iso)
               exit
            endif
            p1 = p2
            p2 = p2 + prob(iti,iso+1)
         enddo
         xmass = xmass + x
      enddo

      current_mass = xmass / float(abs(z_chrg))

      there = .true.
      !> only values larger than user input 
      if ( current_mass > mzmin ) then 

        !> loop over all entries in the list
inner:  do 
          loop = loop + 1

          !>> write if there is no entry
          if ( list_masses(loop) == 0.0_wp ) then
            there = .false.
            exit inner
          !>> true if already in list, end
          elseif  ( abs(list_masses(loop) - current_mass) < 1.0d-10 ) then
            there = .true.
            store_int(loop) = store_int(loop) + 1
            exit inner
          !>> false if not in list, store
          elseif     ( abs(list_masses(loop) - current_mass) > 1.0d-10 ) then
            there = .false.
            cycle inner
          endif
        enddo inner
        loop = 0

        !> if it is not in the list, add it
        if ( .not. there ) then
          index_mass = index_mass + 1
          list_masses(index_mass) = current_mass
          store_int(index_mass) = store_int(index_mass) + 1
          !write(*,*) list_masses(index_mass)
        endif
      endif


    enddo


  allocate(exact_intensity(index_mass))
  allocate(isotope_masses(index_mass))


  do loop = 1, index_mass
    isotope_masses(loop) = list_masses(loop) 
    !exact_intensity(loop)  = real(store_int(loop),wp) / sum(real(store_int,wp))
    exact_intensity(loop)  = float(store_int(loop)) / sum(float(store_int))
    !write(*,*) isotope_masses(loop), store_int(loop)
    !write(*,*) exact_intensity(loop)
  enddo

   


  ! if no isotopes we take here the peak of the isotope pattern with the highest intensity
  if ( no_isotopes ) then 
  iipmax = maxval(exact_intensity)
  indexipmax = maxloc(exact_intensity, dim = 1)
  mipmax = isotope_masses(indexipmax)
  deallocate(exact_intensity,isotope_masses)
  index_mass = 1
  allocate(exact_intensity(1),isotope_masses(1))
  exact_intensity(1)=iipmax
  isotope_masses(1)=mipmax
  endif


end subroutine isotope


end module isotope_pattern
