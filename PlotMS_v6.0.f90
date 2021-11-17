! ============================================================================
! Name           : plotms.f90
!! Authors       : S. Grimme
!! Co-authors    : J.Koopman, C.Bauer
!! Contributions : T. Kind (FiehnLab 2013)
!! Version       : 5.1  (Aug 06 2021)
!! Copyright     : S. Grimme
!! Description   : plot mass spectra from QCxMS
!! ============================================================================
!
!!====================================================================
!!     Original PlotMS Program for extraction of mass spectra
!!     from quantum mechanical simulations used within QCxMS
!!
!!                      *********************************************
!!                      *                S. Grimme                  *
!                      * Mulliken Center for Theoretical Chemistry *
!                      *             Universitaet Bonn             *
!                      *                  2008-22                  *
!                      *********************************************
!
!     Please cite as:
!     S.Grimme, Angew.Chem.Int.Ed. 52 (2013) 6306-6312.
!
!     Involves:
!      - Matching score -- C.Bauer
!      - JCAMP-DX Files -- T.Kind
!
!====================================================================

! 5.0: Changing name from QCEIMS to QCxMS 
! 5.1: Changing from pseudo-random distr. to BoxMuller distribution


program plotms
  use bucket, only: check_bucket
  use isotope_pattern, only: isotope
  use qcxms_boxmuller, only: vary_energies
  use qcxms_readcommon
  use xtb_mctc_accuracy, only: wp
  use xtb_mctc_symbols, only: toSymbol
  implicit none

 !treat i,j,k,l,m,n as integers and all other variables as real (?!?!)
  integer :: n,m,i,j,k,kk,kkk,kkkk,nn,ntot,nsig
  integer :: atm_types
  integer :: maxatm,imax,imin,nagrfile,spec,io
  integer :: iat  (10000)
  integer :: nat  (10000)
  integer :: iat_save (10000)
  integer :: mass (1000)! maximum mass = 10000
  integer :: isec,jcoll,jsec,ial,jal(0:10),kal(0:30),nf,irun,mzmin
  integer :: maxrun
  ! TK number of peaks in in-silico spectra, needed for JCAMP-DX export
  integer :: numspec
  integer :: io_spec, io_raw, io_mass, io_exp, io_jcamp
  integer :: counter
  integer :: z_chrg
  integer :: cnt
  integer :: index_mass, count_mass, sum_index
  integer :: sorted_index

  real(wp) :: xx(100),tmax,r,rms,norm,cthr,cthr2
  real(wp) :: chrg,chrg2,dum,checksum,score
  real(wp) :: iexp (2,10000)
  real(wp) :: mint (1000)
  real(wp) :: tmass(10000) ! tmass contains abundances (real) for each unit mass (i)
  real(wp) :: chrg_wdth
  real(wp) :: total_charge
  real(wp) :: lowest
  real(wp), allocatable :: sorted_masses(:), sorted_intensity(:)
  real(wp), allocatable :: exact_intensity(:)
  real(wp), allocatable :: isotope_masses(:)
  real(wp), allocatable :: rnd(:,:)
  real(wp), allocatable :: checksum2(:), save_mass(:)
  !real(wp), allocatable :: list_masses(:)
  !real(wp), allocatable :: intensity(:)
  real(wp) :: list_masses(10000)
  real(wp) :: intensity(10000)

  logical  :: sel,echo,exdat,mpop,small
  logical  :: ex,ex1,ex2,ex3,ex4
  logical  :: noIso, Plasma
  
  ! fname=<qcxms.res> or result file, xname contains the mass.agr plot file
  character(len=80)  :: arg(10)
  character(len=80)  :: xname
  character(len=:), allocatable  :: fname,fname1,fname2,fname3,fname4
  !character(len=:), allocatable  :: Frag_Name
  character(len=20), dimension(20)  :: Frag_Name

  mpop          = .false.
  echo          = .false.
  small         = .false.
  isec          = 0
  norm          = 1.0
  nagrfile      = 410
  mzmin         = 10
  spec          = 0
  noIso         = .false.
  Plasma        = .false.
 
  !> setup overall (initial) charge as used in QCxMS (has to be set manually)
  cthr          = 1.d-3
  chrg_wdth     = 0.0_wp
  total_charge  = 1.0_wp
  cthr2         = 0.01 ! only used for information, not in program

  z_chrg = 1


 ! edit this path name to some standard xmgrace plot file
 ! JK changed to universal location/file
 
 !xname='~/.mass_raw.agr'
 !xname='~/.mass_jay.agr'
  fname=''
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! start loop reading arguments
  do i=1,9
     arg(i)=' '
    call get_command_argument(i,arg(i))
  enddo
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  ! start loop processing aruments
  
  !  comand line parameters
  !  -a count charges from -1000 (cthr) to cthr2
  !  -v print spectra "mass % intensity  counts   Int. exptl" to stdout;
  !     with "Int. exptl" (experimental) taken from exp.dat but not all exp peaks are exported
  !     if no theoretical counterpart exists
  !  -f filename or  -f <name_of_res_file>
  !  -t couting ions with charge x to y (give the value, e.g. "-t 1 2" for charge 1 to 2)
  !  -w broadening the charges by an SD 
  !  -s Take only secondary and tertiary fragmentations (give the value, e.g. "-s 2" for secondary)
  !  -m set minimum value for m/z, so 100% value will be calc. for higher values (CID)
  !  -i calculate NO isotope pattern
  
  do i = 1, 9
     if(index(arg(i),'-a') /= 0)  cthr   = -1000.0_wp
     if(index(arg(i),'-v') /= 0)  echo   = .true.
     if(index(arg(i),'-f') /= 0)  fname  = arg(i+1)
     if(index(arg(i),'-i') /= 0)  noIso  = .true.
     if(index(arg(i),'-p') /= 0)  Plasma = .true.

     ! if = 0, unity intensities are regarded
     if(index(arg(i),'-t') /= 0) then
        call readl(arg(i+1),xx,nn)
        cthr=xx(1)
     endif

     ! set width distr. of boltz sampling
     if(index(arg(i),'-w') /= 0) then
        call readl(arg(i+1),xx,nn)
        chrg_wdth=xx(1)
     endif

     ! set number of cascading runs
     if(index(arg(i),'-s') /= 0) then
        call readl(arg(i+1),xx,nn)
        isec=int(xx(1))
     endif

     ! set minimum mass
     if(index(arg(i),'-m') /= 0)then
        call readl(arg(i+1),xx,nn)
        mzmin=int(xx(1))
     endif

     ! get initial charge of system
     if(index(arg(i),'-c') /= 0)then
        call readl(arg(i+1),xx,nn)
        total_charge = real(xx(1),wp)
     endif
  enddo

  
!  if(index(arg(1),'-k') /= 0) then
!     xname='~/.mass_jay_klein.agr'
!     write(*,*) ' Using small Version'
!     small = .True.
!  else
! !    xname='~/.mass_jay.agr'
!     xname='~/.mass_raw.agr'
!  endif
  
  xname = '~/.mass_raw.agr'
  
  ! fname contains the results from each calculation or the temporary result tmpqcxms.res
  ! xname contains the xmgrace plot file

  if ( Plasma ) then
    write(*,*) 'Print neutral PLASMA spectrum'
    fname = 'qcxms_neutrals.res'
    spec  = 1
  
  elseif ( fname == '' ) then
    fname1 = 'qcxms.res'
    fname2 = 'qcxms_cid.res'
    fname3 = 'tmpqcxms.res'
    fname4 = 'tmpqcxms_cid.res'

    ! check if either EI or CID res file exist
    inquire ( file = fname1, exist = ex1 )
    inquire ( file = fname2, exist = ex2 )
    inquire ( file = fname3, exist = ex3 )
    inquire ( file = fname4, exist = ex4 )

    if ( ex1 .and. ex2 ) stop 'cannot process both EI and CID .res files!'
    if ( ex3 .and. ex4 ) stop 'cannot process both EI and CID tmp.res files!'

    ! none exist
    if ( .not. ex1 .and. .not. ex2 .and. .not. ex3 .and. .not. ex4 ) stop '! res file does not exist !' 

    ! if EI exists
    if ( ex1 .or. ex3 ) then
      spec = 1
      if ( ex1 ) then
        fname = 'qcxms.res'
      else
        fname = 'tmpqcxms.res'
      endif
    endif

    ! if CID exists
    if ( ex2 .or. ex4 ) then
      spec = 2
      if ( ex2 ) then
        fname = 'qcxms_cid.res'
      else
        fname = 'tmpqcxms_cid.res'
      endif
    endif
  endif

  
  write(*,*)
  write(*,*) '*******************************'
  write(*,*) '* QCxMS output reader PLOTMS  *'
  write(*,*) '*     v. 6.0 Nov 2021         *'
  write(*,*) '*         S. Grimme           *'
  write(*,*) '*         J. Koopman          *'
  write(*,*) '*******************************'
  write(*,*) 'Contributor:  C.Bauer, T.Kind '
  write(*,*)
  write(*,*) 'xmgrace file body ',trim(xname)
  write(*,*) 'Reading ... ', trim(fname)
  write(*,*)
  
  ! ------------------------------------------------------------------------------------------------------!
  ! execute arguments
  
  write(*,*) ' '
  ! -w
  if(chrg_wdth > 0)then
     write(*,'('' broadening the charges by an SD, wdth :'',f4.1)')chrg_wdth*100.
  endif
  
  ! -s
  if(isec /= 0)then
     write(*,*) 'Taking only secondary, tertiary ... fragmentations ',isec
  endif
  
  ! -m
  if(mzmin > 10)then
     write(*,*) 'Only m/z values greater than',mzmin,'are considered'
  endif

  ! -c
  !if ( total_charge > 0 ) then ! always (as information)
  !   write(*,*) 'The total charge of the system is ', total_charge
  !endif

  ! -t (choose if unity intensities or normal)
  if(cthr >= 0)then
   !  write(*,'('' couting ions with charge from '',f4.1,'' to '',f4.1)')   cthr,cthr2 + total_charge
  else
     write(*,'('' counting all fragments with unity charge (frag. overview)'')')
  endif

  write(*,*) ' '
  
  ! ------------------------------------------------------------------------------------------------------!
  tmass = 0
  maxrun = 0
  i=1
  maxatm=0
  counter = 0
  ! ------------------------------------------------------------------------------------------------------!
 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! read .res file for the first time:
  ! charge(chrg2), trajectory(irun),number of collision event(icoll), numbers of (secondary-)fragments (jsec,nf), 
  ! no. of different atom types in the fragment (atm_types), atomic number (iat), 
  ! amount of this atom typ (nat) <== from kk = 1 to k
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  open ( file=fname, newunit=io_spec, status='old' )
    
  do

    counter = counter + 1
    !write(*,*) 'COUNTER', counter

    if (spec == 1) then    !EI
      read(io_spec,*,iostat=io) chrg2, z_chrg, irun, jsec, nf, atm_types, (iat(kk),nat(kk), kk = 1, atm_types)
      if(io<0) exit !EOF
      if(io>0) stop 'Fail in read' !fail
  
    elseif(spec == 2) then !CID
      read(io_spec,*,iostat=io) chrg2, z_chrg, irun, jcoll, jsec, nf, atm_types, (iat(kk), nat(kk), kk = 1, atm_types)
      if(io<0) exit !EOF
      if(io>0) stop 'Fail in read' !fail

    else 
      write(*,*) 'S T O P - Something wrong in plotms - no spec'
      stop
    endif


    if (isec > 0 .and. isec /= jsec) cycle !check tertiary,etc fragmentation (isec)
  
    if ( z_chrg > 0 ) then
      if (chrg2 > cthr) then  !default: chrg2 > 0.001
        ntot = sum(nat(1:atm_types))
        if ( ntot > maxatm ) maxatm = ntot  !get highest number of atoms in fragment
        if ( z_chrg > int(total_charge) ) total_charge = real( z_chrg, wp )
      endif
    elseif ( z_chrg < 0 ) then
      if (chrg2 < -1.0_wp*cthr) then  !default: chrg2 > 0.001
        ntot = sum(nat(1:atm_types))
        if ( ntot > maxatm ) maxatm = ntot  !get highest number of atoms in fragment
        if ( z_chrg < int(total_charge) ) total_charge = real( z_chrg, wp )
      endif

    endif

    i = i + 1 !count number of single fragments with charge > chrt
    if ( irun > maxrun ) maxrun = irun !save highest run number
  
  enddo
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  write(*,*) 'The maximum charge of the system is ', total_charge

  !> allocate the variables
  allocate (checksum2(maxrun), &
    &       save_mass(counter))

  n = i - 1  ! n = no. of fragments with amount maxatm of atoms
  
  write(*,*) n,' fragments with ',maxatm,' atoms max.'
  
  !close file
  close(io_spec)

  !seed=42                                                                                     
  !call random_seed(seed)                                                
  ! initialize the random number array (efficiency)
  allocate(rnd(50000,maxatm))  
  do i = 1, 50000
    do j = 1, maxatm
      call random_number(r)
      rnd(i,j) = r
    enddo
  enddo
  
  
  ! read it again
  write(*,*) 'Computing ...'
  
  ! contains qcxms.res as standard option
  open ( file=fname, newunit=io_spec, status='old')
  imin=100000
  i         = 1
  counter   = 0
  ial       = 0
  jal       = 0
  kal       = 0
  checksum  = 0.0_wp
  checksum2 = 0.0_wp
  save_mass = 0.0_wp
  index_mass= 0
  intensity = 0
  count_mass = 0
  list_masses = 0.0_wp
  sum_index = 0

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> start the loop, processing the qcxms.res file, line-by-line
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  do

    counter = counter + 1

    !> read the fragment entry from qcxms.res line by line
    if (spec == 1) then    !EI
       read(io_spec,*,iostat=io) chrg2, z_chrg, irun, jsec, nf, atm_types, (iat(kk),nat(kk),kk=1,atm_types)
       if(io<0) exit !EOF
       if(io>0) stop !fail
  
    elseif(spec == 2) then !CID
       read(io_spec,*,iostat=io) chrg2, z_chrg, irun, jcoll, jsec, nf, atm_types, (iat(kk),nat(kk),kk=1,atm_types)
       if(io<0) exit !EOF
       if(io>0) stop !fail

    endif


    !do kk = 1, atm_types
    ! 
    !  syb_frag(kk) = iat(kk)
    !  nmb_frag(kk) = nat(kk)
    !  !write(*,'(a2,i2)') toSymbol(iat(kk)), nat(kk)
    !  !write(frag_name,'(a2,i2)') toSymbol(iat(kk)), nat(kk)
    !  !write(*,*) frag_name(kk)
    !enddo

    ! write(*,*) frag_name
  
    sel = .false.
    chrg = chrg2

    !> sum up the charges
    checksum2(irun) = checksum2(irun) + chrg2

    !> manual setting (see above)
    if ( cthr < 0 ) chrg = total_charge

    !> count the moment where fragments were created (first/second MD or collision...)
    !if ( (z_chrg > 0 .and. chrg > cthr) .or. (z_chrg < 0 .and. chrg < -1* cthr) ) then   ! default: yes
    if ( abs(chrg) > cthr )  then   ! default: yes
      sel = .true.
      ial = ial + 1              !count total amount of fragmentations
      jal(jsec) = jal(jsec) + 1  !count amount of first, sec., tert., ... fragmentations

      !> if not fragmented (i.e. smaller/larger than total_charge) count it to kal(1)
      !if ( abs(chrg) < abs(total_charge) ) then
      !  kal(jcoll) = kal(jcoll) + 1 !count the collisions
      !else
      !  kal(1) = kal(jcoll) + 1     !no fragmentation
      !endif

      !> count the collision at which most framentation occured
      if ( spec == 2 ) then
        if ( chrg == total_charge ) then
          kal(1)     = kal(1) + 1 !count the collisions
        else
          kal(jcoll) = kal(jcoll) + 1     !no fragmentation
        endif
      endif

    endif
    !----
    if ( isec > 0 .and. isec /= jsec ) sel = .false.  !Only if asked, not default
    !----

    !> get the specifics of the current fragment entry
    if ( sel ) then ! default
      ! sum all atoms of all types (total atoms)
      ntot = sum(nat(1:atm_types))
      ! all types of atoms
      kkkk = 0
      do kk = 1, atm_types
        ! all atoms of this type
        do kkk = 1, nat(kk)
          kkkk = kkkk + 1
          iat_save(kkkk) = iat(kk) 
        enddo
      enddo

      checksum = checksum + chrg   !Check total charge
    
      !> compute isotope patterns via monte carlo and save intensites 
      !  of all possible combinations 
      call isotope (counter, mzmin, ntot, iat_save, maxatm, rnd, mass, mint, &
        nsig, noIso, index_mass, exact_intensity, isotope_masses, z_chrg)

      !>> distribute charge (if set as input)
      if(chrg_wdth > 1.0d-6) chrg = vary_energies(chrg,chrg_wdth) 

      !> sort the single fragment intensities over the entire list of frags
       if (index_mass > 0) then
      call check_bucket( index_mass, isotope_masses, exact_intensity, &
        list_masses, intensity, count_mass, chrg)
      endif
 
      deallocate(exact_intensity)
      deallocate(isotope_masses)


      !> calculate the total mass
      !do atm_types = 1, nsig
      !  tmass(mass(atm_types)) = tmass(mass(atm_types)) + mint(atm_types) * abs(chrg)
      !enddo




    endif

  enddo
  !> end the loop
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


     !write(*,*) frag_name

  !> sort the list with increasing masses
  allocate(sorted_masses(count_mass), sorted_intensity(count_mass))
  lowest = 0.0_wp
  cnt=1

  !>> find the lowest entry and save it in new lists, repeat with all values
  !   larger than the last value (lowest)
  do i = 1, count_mass
    sorted_index = minloc(list_masses,1, mask = list_masses > lowest)
    lowest = list_masses(sorted_index)
    sorted_masses(i) = list_masses(sorted_index)
    sorted_intensity(i) = intensity(sorted_index)
    !write(*,*) i, sorted_index, sorted_masses(i), sorted_intensity(i)
  enddo
  

  !> write out the sum of charges
  write(*,*) n,' (charged) fragments done.'
  write(*,*)
  write(*,*) 'checksum of charge :', checksum 
  write(*,*)
  
  k=0
  !> count absolute values larger than 1e-6
  do i = 1, maxrun 
    if ( abs(checksum2(i) ) >  1.0d-6) k = k + 1
    if ( abs(checksum2(i) ) > 1.0d-6 .and. abs(checksum2(i) - total_charge ) >  1.0d-3)then
      write(*,*) 'checksum error for trj', i,' chrg=',checksum2(i)
    endif
  enddo
  
  write(*,*) k,' successfull runs.'
  
  write(*,*)
  write(*,*) 'ion sources:'
  write(*,'(''% inital run'',f6.1)') 100.0_wp * float(jal(1)) / float(ial)
  
  do j=2,8
     write(*,'(''%        '',i1,''th'',f6.1)') j ,100.0_wp * float(jal(j))/float(ial)
  enddo 

  if ( spec == 2 ) then
    write(*,*)
    write(*,'(''% collisions '',f6.1)') 100.0_wp * float(kal(1)) / float(ial)
    
    do j=2,16
       write(*,'(''%        '',i2,''th'',f6.1)') j ,100.0_wp * float(kal(j))/float(ial)
    enddo 
  endif
  
  write(*,*)


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! read the provided experimental file (exp.dat) 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  imax=0
  imin=0
  
  inquire(file='exp.dat',exist=exdat)
  
  if(exdat) then
    write(*,*) 'Reading exp.dat ...'
    open(file='exp.dat', newunit=io_exp, status='old')
  
    !> read the experimental (NIST) file 
rd: do
      read(io_exp,'(a)',iostat=iocheck)line
      if (iocheck < 0) exit rd

  
      if(index(line,'##MAXY=') /= 0)then
        line(7:7)=' '
        call readl(line,xx,nn)
        norm=xx(nn)/100.0_wp
      endif
  
      if(index(line,'##PEAK TABLE') /= 0)then
        kk=0
        do
          read(io_exp,'(a)',iostat=iocheck)line
          if (iocheck < 0) exit rd
  
          ! TK JCAMP DX for MS data has "##END="
          if(index(line,'##END') /= 0) cycle  
  
          do k=1,80
            if(line(k:k) == ',')line(k:k)=' '
          enddo
  
          call readl(line,xx,nn)
  
          do k = 1, nn/2
             kk = kk + 1
             iexp(1,kk) = xx(2 * k - 1)
             iexp(2,kk) = xx(2 * k)
             if(iexp(1,kk) > imax) imax = iexp(1,kk)
             if(iexp(1,kk) < imin) imin = iexp(1,kk)
          enddo
        enddo
      endif
    enddo rd
  
    close(io_exp)
  endif
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
  !> get the minimum m/z value that is to be plotted
  j = 0
  !>> sum the total mass only over all values higher than provided (default: 10) 
  !do i = mzmin, 10000
  !  if ( tmass(i) /= 0 ) j = i
  !enddo
  !do i = 1, index_mass
  !  if ( list_masses(i) <= mzmin ) then
  !    final_masses(i) = 0 
  !  endif
  !enddo

  !final_masses = pack(list_masses, list_masses > mzmin)

  imin = mzmin

  !imax = max(j,imax)
  imax = nint(maxval(list_masses))
  !> find the highest entry
  tmax = maxval(intensity)
  
  !tmax = maxval(tmass(mzmin:10000))
  !tmax = maxval(tmass)
  !write(*,*) 'tmax', tmax

  ! echt nicht sicher ob tmass auch das tmass sein soll, aber schauen
  
  ! counts of structures in the 100% peak
  write(*,*) 'Theoretical counts in 100 % signal:', idint(tmax)
 
  !> verbose information
  if(echo)then
    write(*,*)'mass % intensity  counts   Int. exptl'
    do i = mzmin, 10000
       if (tmass(i) /= 0) then
          dum = 0
          do j = 1, kk
             if(int(iexp(1,j)) == i) dum = iexp(2,j)
          enddo

          write(*,'(i8,F8.2,i8,F8.2)') &
          i, 100.0_wp *tmass(i)/tmax, idint(tmass(i)), dum/norm
       endif
    enddo
  endif
  
  ! when the template exists
  
  inquire(file=xname, exist=ex)
  if(ex)then
     open(file=xname, newunit=io_raw)
   
  ! Original file called .mass_raw.agr
     if(       small ) open( file='small_mass.agr', newunit=io_mass)
     if( .not. small ) open( file='mass.agr',       newunit=io_mass)
  !   open(unit=9,file='intensities.txt')
  
  ! my xmgrace file mass.agr has nagrfile lines
     write(*,*)
     write(*,*) 'Writing mass.agr ...'
     do i = 1, nagrfile
        read (io_raw,'(a)') line
        if( .not. small )then 
           if(index(line,'world 10') /= 0)then
              write(line,'(''@    world '',i3,'', -105,'',i3,'', 100'')') imin, imax + 5
           endif
           if(index(line,'xaxis  tick major') /= 0)then
              line='@    xaxis  tick major 20'
                 if(imax > 200) line='@    xaxis  tick major 50' 
           endif
  
           if(index(line,'@    s1 symbol size') /= 0)then
              line='@    s1 symbol size 0.160000'
              if(imax > 200) line='@    s1 symbol size 0.100000'
           endif
  
           if(index(line,'@    s2 symbol size') /= 0)then
              line='@    s2 symbol size 0.160000'
              if(imax > 200) line='@    s2 symbol size 0.100000'
           endif
        !kleines Format
        else
           if(index(line,'world 10') /= 0)then
              write(line,'(''@    world '',i3,'', 0,'',i3,'', 100'')') imin, imax + 5
           endif  
  !         if(index(line,'xaxis  tick major') /= 0)then
  !            line='@    xaxis  tick major 10'
  !!               if(imax > 200) line='@    xaxis  tick major 50' 
  !         endif
           
           if(index(line,'@    s1 symbol size') /= 0)then
              line='@    s1 symbol size 0.320000'
  !            if(imax > 200) line='@    s1 symbol size 0.100000'
           endif
        endif
        write(io_mass,'(a)')line
     enddo
  
  ! write only masses > mzmin
     !do i = mzmin, 10000
     !   if(tmass(i) /= 0)then
     !      write(io_mass,*) i, 100.0_wp * tmass(i)/tmax
!    !       write(9,*) i,100.*tmass(i)/tmax   !ce50values
     !   endif
     !enddo
     !write(io_mass,*)'&'
!     do i = 1, index_mass
     !do i = 1, count_mass
     !   if(tmass(i) /= 0)then
     !      write(*,*) list_masses(i) , 100.0_wp * tmass(i) !/tmax
     !      write(io_mass,*) list_masses(i) , 100.0_wp * tmass(i)!/tmax
!    !       write(9,*) i,100.*tmass(i)/tmax   !ce50values
     !   endif
     !enddo
    do i = 1, count_mass 
!      write(*,*) i, sorted_masses(i), sorted_intensity(i) / tmax
      write(io_mass,*)  sorted_masses(i), 100.0_wp * (sorted_intensity(i) / tmax)

    enddo
     write(io_mass,*)'&'
  
!     close(9)
  
  ! TK establish again if exdat (experimental JCAMPDX) exists
  ! TK @target G0.S1 means upper theory spectrum in current implementation
  ! TK @target G0.S2 means lower experimental spectrum
  ! TK iexp(1,i) - contains the exp masses and  iexp(2,i) contains the abundance
  ! TK kk contains number of experimental spectra
  
     if(exdat)then
       write(io_mass,*)'@target G0.S2'
       write(io_mass,*)'@type bar'
       do i=1,kk
          write(io_mass,*) iexp(1,i),-iexp(2,i)/norm
       enddo
       write(io_mass,*)'&'
     endif
  endif
  
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
  ! TK here we export the spectrum to JCAMP-DX, subroutine would be better for
  ! TK allowing later code merges or inclusion of new and missing elements such as
  ! TK currently unit mass only
  
  write(*,*) 'Writing JCAMP file: result.jdx ...'
  write(*,*)
  ! TK open file as JCAMP-DX MS result file, open new or replace
  open( file='result.jdx', STATUS="REPLACE", newunit=io_jcamp)
  
  ! TK minimum JCAMP-DX header
  write(io_jcamp,"(A)")'##TITLE=Theoretical in-silico spectrum (QCxMS)'
  write(io_jcamp,"(A)")'##JCAMP-DX=4.24'
  write(io_jcamp,"(A)")'##DATA TYPE=MASS SPECTRUM'
  
  !write(io_jcamp,*)'##MOLFORM=C13 H22 O3 Si2'
  !write(io_jcamp,*)'##MW=282'
  write(io_jcamp,"(A)")'##XUNITS=M/Z'
  write(io_jcamp,"(A)")'##YUNITS=RELATIVE INTENSITY'
  !write(io_jcamp,*)'##XFACTOR=1'
  !write(io_jcamp,*)'##YFACTOR=1'
  !write(io_jcamp,*)'##FIRSTX=15'
  !write(io_jcamp,*)'##LASTX=281'
  !write(io_jcamp,*)'##FIRSTY=20'
  !write(io_jcamp,*)'##MAXX=281'
  !write(io_jcamp,*)'##MINX=15'
  !write(io_jcamp,*)'##MAXY=9999'
  !write(io_jcamp,*)'##MINY=10'
  
  ! TK calculate number of in-silico spectra
  ! TK new: numspec = idint(tmax) ** not related to number ofr spectral peaks
  numspec = 0
  do i=10,10000
    if(tmass(i) /= 0)then
      numspec = numspec + 1
    endif
  enddo
  
  
  write(io_jcamp,'(A, I0)')'##NPOINTS=' ,numspec
  
  ! TK ##XYDATA=(XY..XY) 1 designates one line per m/z abd pair
  ! TK separated by comma
  ! TK 1.500000e+001 ,5.000000e+000
  
  write(io_jcamp,"(A)")'##PEAK TABLE=(XY..XY) 1'
  !15,20 26,10 27,20 29,50
  
  
  do i=10,10000
       if(tmass(i) /= 0)then
         ! unit mass to exact mass
         ! write(io_jcamp,"(I4, I8)") i, int(100.*tmass(i)/tmax)
         write(io_jcamp,"(F12.6, I8)") real(i), int(100.0_wp * tmass(i)/tmax)
       endif
  enddo
  
  write(io_jcamp,"(A)")'##END='
  
  ! TK close the JCAMP-DX file
  close(io_jcamp)
  
 ! compute deviation exp-theor.
  
  if(exdat)then

    rms  = 0
    nn   = 0
    kkk  = 0
    kkkk = 0

    do i = 1, kk
      k = iexp(1,i)
      if ( iexp(2,i)/norm > 5.0_wp )then
        r = 100.0_wp * tmass(k)/tmax - iexp(2,i)/norm
        kkk = kkk + 1
        if ( 100.0_wp * tmass(k)/tmax > 2.5_wp) kkkk = kkkk + 1
        rms = rms + abs(r)
      endif
    enddo
    write(*,*)'MAD(exptl./theor.) = ', rms / kkk
    write(*,*)'# exptl. > 5 %     = ', kkk
    write(*,*)'% correctly found  = ', 100 * float(kkkk)/float(kkk)

  endif
  
  close(io_spec)
  close(io_raw)
  close(io_mass)
  close(io_jcamp)
  
  ! compute mass spectral match score
  if(exdat)then 
    score = 0.0_wp
    ! bring iexp to the right order of magnitude 
    do i = 1, 10000
      iexp(2,i) = iexp(2,i) / norm
 enddo

    call match(tmass,iexp,score,tmax)
    write(*,*)
    write(*,*)"!!!!!!!!!!!!!!!!!!!!!!! "
    write(*,*)"  Matching score:  "
    write(*,'(6F10.3)') score   
    write(*,*)"!!!!!!!!!!!!!!!!!!!!!!! "
    write(*,*)
    write(*,*)"Composite match score, see "
    write(*,*)"Stein, S. E.; Scott, D. R. J. Am. Soc. Mass Spectrom. 1994, 5, 859-866" 
    write(*,*)"For our implementation, see "
    write(*,*)"Bauer, C. A.; Grimme,S. J. Phys. Chem. A 2014, 118, 11479-11484"
    
    ! Overlap score tried...but failes
    !      call match2(tmass,iexp,score,1.0d0)
    !      write(*,'(''overlap score'')')
    !      write(*,'(6F10.3)') score
  
  endif 
   
end program plotms

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine match2(tmass,iexp,score,alpha)
  use xtb_mctc_accuracy, only: wp
  implicit none

  integer  :: i,j

  real(wp) :: tmass(10000)
!  real(wp) :: int_exp(10000)
  real(wp) :: iexp(2,10000)
  real(wp) :: score
  real(wp) :: alpha

  ! needed storage variables
  real(wp) :: norm_experiment
  real(wp) :: norm_calculation
  real(wp) :: overlap
  
  !      int_exp = 0.0d0
  !      do i=1,10000
  !          print*,i,iexp(1,i),iexp(2,i)
  !        if(i == iexp(1,i)) then
  !          int_exp(int(iexp(1,i))) = iexp(2,i)
  !        endif
  !      enddo
  
  norm_experiment = 0.0d0
  norm_calculation = 0.0d0
  overlap = 0.0d0
  do i = 1, 10000
    do j = 1, 10000
      norm_experiment = norm_experiment + (-1)**(i+j) * iexp(2,j) * iexp(2,i) * exp(-alpha*(i-j)**2)
      norm_calculation = norm_calculation + (-1)**(i+j) * tmass(j) * tmass(i) * exp(-alpha*(i-j)**2)
      overlap = overlap + (-1)**(i+j) * iexp(2,j) * tmass(i) * exp(-alpha * (i-j)**2)
    enddo
  enddo
  
  score = overlap / sqrt(norm_experiment * norm_calculation)
  
end subroutine match2
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine match(tmass,iexp,score,tmax)
  use xtb_mctc_accuracy, only: wp
  implicit none

  integer :: i,j
  integer :: pp !peak pair number pp
  integer :: numcalc ! number of calculated peaks

  real(wp) :: tmass(10000)
  real(wp) :: iexp(2,10000)
  real(wp) :: score
  real(wp) :: tmax
  real(wp) :: w(2,10000) !weighted vectors
  real(wp) :: sum1,sum2,sum3,sum4,dot,fr
  real(wp) :: norm(2)
  real(wp) :: p(2,10000) ! peak pair matrix 

  numcalc=0
  
  do i = 1, 10000
    if ( tmass(i) /= 0.0_wp ) numcalc = numcalc + 1
  enddo
  
  !weighted spectral vectors, m**3, int**0.6 scaling
  w = 0.0_wp
    do i = 1, 10000
      do j = 1, 10000
        if ( j == iexp(1, i) ) then
          w(2,j) = j**2 * ( 1000.0_wp * tmass(j) / tmax )**0.6_wp ! changed this to 1000 <- check it!
          w(1,j) = iexp(1,i)**2 * iexp(2,i)**0.6_wp
        endif
      enddo
    enddo
  
  norm = 0.0_wp
  do i = 1, 10000
     norm(1) = norm(1) + w(1,i)**2
     norm(2) = norm(2) + w(2,i)**2
  enddo
  
  norm = sqrt(norm)
  
  do i = 1, 10000
     w(1,i) = w(1,i) / norm(1)
     w(2,i) = w(2,i) / norm(2)
  enddo
  
  dot  =  0.0_wp
  sum1 =  0.0_wp
  sum2 =  0.0_wp
  sum3 =  0.0_wp

  do i = 1, 10000
     sum1 = sum1 + w(1,i) * w(2,i)
     sum2 = sum2 + (w(1,i))**2
     sum3 = sum3 + (w(2,i))**2
  enddo
  
  dot = sum1**2 / (sum2 * sum3)
  
  ! masses and intensity scalement for the second term
  ! m**0
  ! int**1
  
  w=0.0_wp

  do i = 1, 10000
     do j = 1, 10000
        if ( j == iexp(1,i) ) then
           w(2,j) = 100.0_wp * tmass(j) / tmax
           w(1,j) = iexp(2,i)
        endif
     enddo
  enddo
  
  !calculate the norm
  
  norm=0.0_wp
  do i=1, 10000
     norm(1) = norm(1) + w(1,i)**2
     norm(2) = norm(2) + w(2,i)**2
  enddo
  
  norm = sqrt(norm)
  
  do i = 1, 10000
     w(1,i) = w(1,i) / norm(1)
     w(2,i) = w(2,i) / norm(2)
  enddo
  
  pp   = 0
  sum4 = 0.0_wp
  fr   = 0.0_wp
  p    = 0.0_wp

  do i = 1, 10000
     if ( (w(1,i) * w(2,i) ) /= 0) then
        pp = pp + 1
        p(1,pp) = w(1,i)
        p(2,pp) = w(2,i)
  !         print*,p(1,i),p(2,i)
     endif
  enddo
  
  ! ardous loop 
  call calcfr(pp,p,sum4)
  if (pp == numcalc) sum4 = sum4 + 1.0_wp
  
  fr = sum4 / float(pp)
  score = (numcalc * dot + pp * fr) / (numcalc + pp)
  return
end subroutine match

 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! computes the second expression in the mass spec match score 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine calcfr(pp,pair,sum4)
  use xtb_mctc_accuracy, only: wp
  implicit none

  integer :: i
  integer :: pp

  real(wp) :: pair(2,pp)
  real(wp) :: sum4
  
  sum4 = 0.0_wp
  
  do i = 1, pp
     if     (abs((pair(1,i) * pair(2,i-1)) / (pair(1,i-1) * pair(2,i)))   < 1.0_wp)then
  
        sum4 = sum4 + (pair(1,i) * pair(2,i-1)) / (pair(1,i-1) * pair(2,i))
     
     elseif (abs((pair(1,i) * pair(2,i-1)) / (pair(1,i-1) * pair(2,i)))  > 1.0d0)then
  
        sum4 = sum4 + 1 / ((pair(1,i) * pair(2,i-1)) / (pair(1,i-1) * pair(2,i))) 
     
     elseif (abs((pair(1,i) * pair(2,i-1)) / (pair(1,i-1) * pair(2,i)))  == 1.0d0)then
        sum4 = sum4 + 1.0_wp
     endif
  enddo
  
end subroutine calcfr
