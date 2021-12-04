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
  use count_entries, only: check_entries
  use isotope_pattern, only: isotope
  use qcxms_boxmuller, only: vary_energies
  use qcxms_readcommon
  use xtb_mctc_accuracy, only: wp
  use xtb_mctc_symbols, only: toSymbol
  implicit none

  integer :: n,i,j,k,kk,kkk,kkkk,nn,ntot
  integer :: atm_types
  integer :: maxatm,imax,imin,nagrfile
  integer :: spec
  integer :: iat  (10000)
  integer :: nat  (10000)
  integer :: iat_save (10000)
  integer :: isec,jcoll,jsec,ial,jal(0:10),kal(0:30),nf,irun,mzmin
  integer :: maxrun
  ! TK number of peaks in in-silico spectra, needed for JCAMP-DX export
  integer :: io, io_spec, io_raw, io_mass, io_csv, io_exp, io_jcamp
  integer :: counter
  integer :: z_chrg
  integer :: index_mass, count_mass, sum_index
  integer :: sorted_index
  integer :: list_length
  integer :: exp_entries

  real(wp) :: xx(100),tmax,r,rms,norm,cthr,cthr2
  real(wp) :: chrg,chrg2,checksum,score
  !real(wp) :: dum
  real(wp) :: chrg_wdth
  real(wp) :: total_charge
  real(wp) :: lowest
  real(wp) :: min_intensity 
  real(wp), allocatable :: exp_mass (:)
  real(wp), allocatable :: exp_int  (:)
  real(wp), allocatable :: reduced_intensities(:), reduced_masses(:)
  real(wp), allocatable :: sorted_masses(:), sorted_intensities(:)
  real(wp), allocatable :: exact_intensity(:)
  real(wp), allocatable :: isotope_masses(:)
  real(wp), allocatable :: rnd(:,:)
  real(wp), allocatable :: checksum2(:)
  !real(wp), allocatable :: list_masses(:)
  !real(wp), allocatable :: intensity(:)
  real(wp) :: list_masses(10000)
  real(wp) :: intensity(10000)

  logical  :: sel,verbose,exdat,small
  logical  :: ex,ex1,ex2,ex3,ex4
  logical  :: noIso, Plasma
  logical  :: args = .false.
  
  ! fname=<qcxms.res> or result file, xname contains the mass.agr plot file
  character(len=80)  :: arg(10)
  character(len=80)  :: xname
  character(len=:), allocatable  :: fname,fname1,fname2,fname3,fname4

  verbose          = .false.
  small         = .false.
  isec          = 0
  norm          = 1.0
  nagrfile      = 410
  mzmin         = 10
  spec          = 0
  noIso         = .false.
  Plasma        = .false.
 
  cthr          = 1.d-3
  chrg_wdth     = 0.0_wp
  total_charge  = 1.0_wp
  cthr2         = 0.01 ! only used for information, not in program
  min_intensity = 0.0_wp

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
    if (arg(i) /= ' ') args = .true. 
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
  !  -m set minimum value for m/z, so 100% value will be calc. for higher values (x-axis)
  !  -i set the minimum rel. intensity at which the signals are counted          (y-axis)
  !  -p DO NOT calculate isotope pattern
  
  do i = 1, 9
     if(index(arg(i),'-a') /= 0)  cthr   = -1000.0_wp
     if(index(arg(i),'-v') /= 0)  verbose   = .true.
     if(index(arg(i),'-f') /= 0)  fname  = arg(i+1)
     if(index(arg(i),'-p') /= 0)  noIso  = .true.

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

     ! set the minimum intensity
     if(index(arg(i),'-i') /= 0)then
        call readl(arg(i+1),xx,nn)
        min_intensity = real(xx(1),wp)
     endif
  enddo

  
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
  write(*,'(6x,''* * * * * * * * * * * * * * * * * *'')')
  write(*,'(6x,''*           P l o t M S           *'')')
  write(*,'(6x,''* * * * * * * * * * * * * * * * * *'')')
  write(*,'(6x,''*    QCxMS spectra plotting tool  *'')')
  write(*,'(6x,''*         -  v. 6.0.1  -          *'')')
  write(*,'(6x,''*          04. Dec 2021           *'')')
  write(*,'(6x,''*           *   *   *             *'')')
  write(*,'(6x,''*           S. Grimme             *'')')
  write(*,'(6x,''*           J. Koopman            *'')')
  write(*,'(6x,''* * * * * * * * * * * * * * * * * *'')')
  write(*,'(6x,''   Contributor:  C.Bauer, T.Kind   '')')
  write(*,*)
  write(*,'(6x,''   -> Reading file: '',(a)         )')trim(fname)
  if (verbose) write(*,'(6x ''   -> xmgrace file body '',(a)         )')trim(xname)
  write(*,*)
  write(*,*)
  
  ! ------------------------------------------------------------------------------------------------------!
  ! execute arguments
 
  if (args) then
    write(*,'(50(''!''))')
    write(*,'('' The following settings are being used :  '')')
    write(*,*)

    ! -w
    if(chrg_wdth > 0)then
      write(*,'('' - Broadening of the charges by an SD, wdth :'',f4.1)')chrg_wdth*100.
    endif
    
    ! -s
    if(isec /= 0)then
      write(*,'('' - Taking only secondary, tertiary ... fragmentations '')') isec
    endif
    
    ! -m
    if(mzmin > 10)then
      write(*,'('' - Only m/z values greater than '',i4,'' are considered'')') mzmin
    endif

    ! -i
    if ( min_intensity > 0.0_wp ) then 
      write(*,'('' - Minimum signal intensity: '',f7.4,'' % (relative)'')') min_intensity
    endif

    ! -p 
    if ( noIso ) then 
      write(*,'('' - No isotope pattern calculated '')') 
    endif

    ! -t (choose if unity intensities or normal)
    if(cthr >= 0)then
      write(*,'('' - couting ions with charge from '',f4.1,'' to '',f4.1)')   cthr,cthr2 + total_charge
    else
      write(*,'('' - counting all fragments with unity charge (frag. overview)'')')
    endif
    write(*,*)
    write(*,'(50(''!''))')
  endif

  ! ------------------------------------------------------------------------------------------------------!
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
      read(io_spec,*,iostat=io) chrg2, irun, jsec, nf, atm_types, (iat(kk),nat(kk), kk = 1, atm_types)
      if(io<0) exit !EOF
      if(io>0) stop ' -- Fail in read -- ' !fail
  
    elseif(spec == 2) then !CID
      read(io_spec,*,iostat=io) chrg2, irun, jcoll, jsec, nf, atm_types, (iat(kk), nat(kk), kk = 1, atm_types)
      if(io<0) exit !EOF
      if(io>0) stop ' -- Fail in read -- ' !fail

    else 
      write(*,*) 'S T O P - Something wrong in plotms - no spec'
      stop
    endif


    if (isec > 0 .and. isec /= jsec) cycle !check tertiary,etc fragmentation (isec)
  
    if (chrg2 > cthr) then  !default: chrg2 > 0.001
      ntot = sum(nat(1:atm_types))
      if ( ntot > maxatm ) maxatm = ntot  !get highest number of atoms in fragment
    endif

    i = i + 1 !count number of single fragments with charge > chrt
    if ( irun > maxrun ) maxrun = irun !save highest run number
  
  enddo
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !write(*,*) 'The maximum charge of the system is ', total_charge

  !> allocate the variables
  allocate (checksum2(maxrun))

  n = i - 1  ! n = no. of fragments with amount maxatm of atoms
  
  write(*,*)
  write(*,'(''--------------------------------------- '')')
  write(*,'(i4, '' fragments with '', i3 ,'' atoms max.'')') n, maxatm
  write(*,'(''--------------------------------------- '')')
  write(*,*)
  
  !close file
  close(io_spec)

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
  write(*,*)
  
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
       read(io_spec,*,iostat=io) chrg2,  irun, jsec, nf, atm_types, (iat(kk),nat(kk),kk=1,atm_types)
       if(io<0) exit !EOF
       if(io>0) stop !fail
  
    elseif(spec == 2) then !CID
       read(io_spec,*,iostat=io) chrg2, irun, jcoll, jsec, nf, atm_types, (iat(kk),nat(kk),kk=1,atm_types)
       if(io<0) exit !EOF
       if(io>0) stop !fail

    endif


    sel = .false.
    chrg = chrg2

    !> sum up the charges
    checksum2(irun) = checksum2(irun) + chrg2

    !> manual setting (see above)
    if ( cthr < 0 ) chrg = total_charge

    !> count the moment where fragments were created (first/second MD or collision...)
    if ( abs(chrg) > cthr )  then   ! default: yes
      sel = .true.
      ial = ial + 1              !count total amount of fragmentations
      jal(jsec) = jal(jsec) + 1  !count amount of first, sec., tert., ... fragmentations


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
      call isotope (counter, mzmin, ntot, iat_save, maxatm, rnd, &
        noIso, index_mass, exact_intensity, isotope_masses)

      !>> distribute charge (if set as input)
      if(chrg_wdth > 1.0d-6) chrg = vary_energies(chrg,chrg_wdth) 

      !> sort the single fragment intensities over the entire list of frags
       if (index_mass > 0) then
        call check_entries( index_mass, isotope_masses, exact_intensity, &
          list_masses, intensity, count_mass, chrg)
      endif
 
      deallocate(exact_intensity)
      deallocate(isotope_masses)

    endif

  enddo
  !> end the loop
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  list_length=0

  !> find the highest intensity
  tmax = maxval(intensity)

  intensity(:) = intensity(:)/tmax * 1000.0_wp

  allocate(reduced_masses(count_mass), reduced_intensities(count_mass))

  reduced_masses      = 0.0_wp
  reduced_intensities = 0.0_wp

  !> sort out everything we don't like
  do i = 1, count_mass
    if ( (intensity(i)) >= min_intensity ) then
      list_length = list_length + 1
      reduced_intensities(list_length) = intensity(i)
      reduced_masses(list_length)      = list_masses(i)
    endif
  enddo

  !> sort the list with increasing masses
  allocate(sorted_masses(list_length), sorted_intensities(list_length))
  lowest = 0.0_wp

  !>> find the lowest entry and save it in new lists, repeat with all values
  !   larger than the last value (lowest)
  do i = 1, list_length
    sorted_index = minloc(reduced_masses,1, mask = reduced_masses > lowest)
    lowest = reduced_masses(sorted_index)
    sorted_masses(i)    = reduced_masses(sorted_index)
    sorted_intensities(i) = reduced_intensities(sorted_index)
    !write(*,*) i, sorted_index, sorted_masses(i), sorted_intensities(i)
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
  
  write(*,'(''* * * * * * * * * * * * * * * * * * * * '')')
  write(*,*) k,' successfull runs!! '
  write(*,'(''* * * * * * * * * * * * * * * * * * * * '')')
  
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
! read the provided experimental file (exp.dat) ! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  imax=0
  imin=0
  
  inquire(file='exp.dat',exist=exdat)
  
  if(exdat) then
    write(*,*) 'Reading exp.dat ...'
    open( file = 'exp.dat', newunit = io_exp, status = 'old')
  
    !> read the experimental (NIST) file 
rd: do
      read(io_exp,'(a)',iostat=iocheck)line
      if (iocheck < 0) exit rd

  
      !> read the maximum intensity of exp.
      if(index(line,'##MAXY=') /= 0)then
        line(7:7)=' '
        call readl(line,xx,nn)
        norm=xx(nn)/100.0_wp
      endif

      !> allocate the number of entries
      if(index(line,'##NPOINTS=') /= 0)then
        line(10:10)=' '
        call readl(line,xx,nn)
        exp_entries = nint(xx(nn))
        !write(*,*) 'EXP ENTRIES', exp_entries
        allocate (exp_mass(exp_entries), &
          &       exp_int(exp_entries))
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
             exp_mass(kk) = xx(2 * k - 1)
             exp_int (kk) = xx(2 * k)
             if(exp_mass(kk) > imax) imax = exp_mass(kk)
             if(exp_int (kk) < imin) imin = exp_int (kk)
          enddo
        enddo
      endif
    enddo rd
  
    close(io_exp)
  endif
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
  !> get the minimum m/z value that is to be plotted
  imin = mzmin

  !> get the maximum m/z value that is to be plotted
  imax = nint(maxval(sorted_masses))

  ! counts of structures in the 100% peak
  write(*,*) 'Theoretical counts in 100 % signal:', idint(tmax)
 
  !----------------------------------------------------------------------------- 
  !> verbose information -> not yet available anymore with acc masses
  !if(verbose)then
  !  write(*,*)'mass % intensity  counts   Int. exptl'
  !  do i = 1, list_length !count_mass
  !     if (sorted_masses(i) > mzmin) then
  !        dum = 0
  !        do j = 1, exp_entries
  !           !if(int(iexp(1,j)) == i) dum = iexp(2,j)
  !           if(int(iexp(1,j)) == i) dum = iexp(2,j)
  !        enddo

  !        write(*,'(i8,F8.2,i8,F8.2)') &
  !        !i, 100.0_wp *sorted_masses(i)/tmax, idint(tmass(i)), dum/norm
  !        !i, 100.0_wp *sorted_intensities(i)/tmax, sorted_masses(i), dum/norm
  !        i, sorted_intensities(i), sorted_masses(i), dum/norm
  !     endif
  !  enddo
  !endif
  !----------------------------------------------------------------------------- 


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! P  R  I  N  T  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  write(*,*)
  write(*,'(''  Writing ...'')')
  !> print out comma-separated-value (csv) file

  write(*,'('' - result.csv file'')')
  open( file = 'result.csv', newunit = io_csv,   status='replace')

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> print out JCAMP file

  write(*,'('' - result.jdx file'')')
  open( file = 'result.jdx', newunit = io_jcamp, status='replace')

  ! TK minimum JCAMP-DX header
  write(io_jcamp,"(A)")'##TITLE=Theoretical in-silico spectrum (QCxMS)'
  write(io_jcamp,"(A)")'##JCAMP-DX=4.24'
  write(io_jcamp,"(A)")'##DATA TYPE=MASS SPECTRUM'
  
  write(io_jcamp,"(A)")'##XUNITS=M/Z'
  write(io_jcamp,"(A)")'##YUNITS=RELATIVE INTENSITY'
  
  !> number of in-silico spectra
  write(io_jcamp,'(A, I0)')'##NPOINTS=' , list_length   
  write(io_jcamp,"(A)")'##PEAK TABLE=(XY..XY) 1'
  
  !       >>>> Write out the results into the files <<<<         !
  do i = 1, list_length 
    write(io_csv,  '(f12.6, 1x, a1, 3x, f26.18)')  &
      sorted_masses(i),',', sorted_intensities(i) 
    write(io_jcamp,'(f12.6,4x, f26.18)')           &
      sorted_masses(i),     sorted_intensities(i) 
  enddo

  write(io_jcamp,"(A)")'##END='


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> print out XMGRACE file

  !> when the template exists
  inquire(file=xname, exist=ex)
  if(ex)then
    open(file=xname, newunit=io_raw)
   
    ! Original file called .mass_raw.agr
    if(       small ) open( file='small_mass.agr', newunit=io_mass)
    if( .not. small ) open( file='mass.agr',       newunit=io_mass)
  
    write(*,'('' - mass.agr file'')')
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
        !if(index(line,'xaxis  tick major') /= 0)then
        !line='@    xaxis  tick major 10'
        !if(imax > 200) line='@    xaxis  tick major 50' 
        !endif
         
        if(index(line,'@    s1 symbol size') /= 0)then
          line='@    s1 symbol size 0.320000'
        !if(imax > 200) line='@    s1 symbol size 0.100000'
        endif
      endif
      write(io_mass,'(a)')line
    enddo
 
    !> Write out the results into the mass.agr file !
    do i = 1, list_length ! count_mass 
      write(io_mass,*) sorted_masses(i), sorted_intensities(i)/10.0_wp 
    enddo

    write(io_mass,*)'&'
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! TK establish again if exdat (experimental JCAMPDX) exists
    ! TK exp_mass(i) - contains the exp masses and  exp_int(i) contains the abundance
    ! TK  exp_entries contains number of experimental spectra
  
    if(exdat)then
      write(io_mass,*)'@target G0.S2'
      write(io_mass,*)'@type bar'
      do i=1,exp_entries
         write(io_mass,*) exp_mass(i),-exp_int(i)/norm
      enddo
      write(io_mass,*)'&'
    endif

  endif !ex

  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
 ! compute deviation exp-theor.
  if(exdat)then

    rms  = 0
    nn   = 0
    kkk  = 0
    kkkk = 0

ol: do i = 1, exp_entries
il:   do j = 1, list_length
        if ( exp_int(i)/norm > 5.0_wp )then !.and. sorted_intensities(j) > 1.0_wp )then
          if ( exp_mass(i) == nint(sorted_masses(j))) then
            r = sorted_intensities(j) - (exp_int(i) / norm)
            !write(*,*) i, exp_mass(i), exp_int(i)
            !write(*,*) j, sorted_masses(j), sorted_intensities(j) 
            !write(*,*) r
            kkk = kkk + 1
            if ( sorted_intensities(j) > 2.5_wp) kkkk = kkkk + 1
            rms = rms + abs(r)
          endif
        endif
      enddo il
    enddo ol


    write(*,*)'MAD(exptl./theor.) = ', rms / kkk
    write(*,*)'# exptl. > 5 %     = ', kkk
    write(*,*)'% correctly found  = ', 100 * float(kkkk)/float(kkk)

  endif
  
  close(io_spec)
  close(io_raw)
  close(io_mass)
  close(io_csv)
  close(io_jcamp)
  
  ! compute mass spectral match score
  if(exdat)then 
    score = 0.0_wp
    ! bring exp to the right order of magnitude 
    do i = 1, exp_entries
      exp_int(i) = exp_int(i) / norm
    enddo

    write(*,*) 
    write(*,*) ' !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! '
    write(*,*) ' - NIST mass comparison not supported in v6.0 - '
    write(*,*) '   -    no automatic score comparison       -'
    stop       ' !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! '


  !  call match(sorted_masses, sorted_intensities, list_length, &
  !              exp_entries, exp_mass, exp_int,score,tmax)
    write(*,*)
    write(*,*)"!!!!!!!!!!!!!!!!!!!!!!! "
    write(*,*)"  Matching score:  "
   ! write(*,'(6F10.3)') score   
    write(*,*)"!!!!!!!!!!!!!!!!!!!!!!! "
    write(*,*)
    write(*,*)"Composite match score, see "
    write(*,*)"Stein, S. E.; Scott, D. R. J. Am. Soc. Mass Spectrom. 1994, 5, 859-866" 
    write(*,*)"For our implementation, see "
    write(*,*)"Bauer, C. A.; Grimme,S. J. Phys. Chem. A 2014, 118, 11479-11484"
    
  
  endif 
   
end program plotms

 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine match(sorted_masses, sorted_intensities, list_length, &
                  exp_entries,  exp_mass, exp_int,score,tmax)
  use xtb_mctc_accuracy, only: wp
  implicit none

  integer :: i,j
  integer :: pp !peak pair number pp
  integer :: numcalc ! number of calculated peaks
  integer :: list_length, exp_entries

  !real(wp) :: tmass(10000)
  !real(wp) :: iexp(2,10000)
  real(wp) :: exp_mass(exp_entries)
  real(wp) :: exp_int(exp_entries)
  real(wp) :: score
  real(wp) :: tmax
  real(wp) :: w(2,10000) !weighted vectors
  real(wp) :: sum1,sum2,sum3,sum4,dot,fr
  real(wp) :: norm(2)
  real(wp) :: p(2,10000) ! peak pair matrix 
  real(wp) :: sorted_masses(list_length), sorted_intensities(list_length)

  numcalc=0
  
  !do i = 1, 10000
  !  if ( tmass(i) /= 0.0_wp ) numcalc = numcalc + 1
  !enddo

  numcalc = list_length
  
  !weighted spectral vectors, m**3, int**0.6 scaling
  w = 0.0_wp
    !do i = 1, 10000
    !  do j = 1, 10000
    !    if ( j == iexp(1, i) ) then
    !      w(2,j) = j**2 * ( 1000.0_wp * tmass(j) / tmax )**0.6_wp ! changed this to 1000 <- check it!
    !      w(1,j) = iexp(1,i)**2 * iexp(2,i)**0.6_wp
    !    endif
    !  enddo
    !enddo
  !do i = 1, 10000
ol: do i = 1, exp_entries
il:  do j = 1, list_length
       if ( sorted_masses(j) == exp_mass(i) ) then
        !w(2,j) = j**2 * ( 100.0_wp * tmass(j) / tmax )**0.6_wp ! changed this to 1000 <- check it!
        w(2,j) = j**2 * ( sorted_intensities(j) )**0.6_wp ! changed this to 1000 <- check it!
        w(1,j) = exp_mass(i)**2 * exp_int(i)**0.6_wp
      endif
    enddo il
  enddo ol
  
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
     !do j = 1, 10000
     do j = 1, list_length
        if ( j == exp_mass(i) ) then
           w(2,j) = sorted_intensities(j) 
           w(1,j) = exp_int(i)
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
