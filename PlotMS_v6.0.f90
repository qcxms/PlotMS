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
!                      *                  2008-21                  *
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
  use qcxms_boxmuller, only: vary_energies
  use xtb_mctc_accuracy, only:wp
  implicit none

 !treat i,j,k,l,m,n as integers and all other variables as real (?!?!)
  integer :: n,i,j,k,kk,kkk,kkkk,nn,ntot,nsig
  integer :: atm_types
  integer :: maxatm,imax,imin,nagrfile,spec,io
  integer :: iat  (10000)
  integer :: nat  (10000)
  integer :: iat_sum (10000)
  integer :: mass (1000)! maximum mass = 10000
  integer :: isec,jcoll,jsec,ial,jal(0:10),nf,irun,seed,mzmin
  integer :: total_charge
  ! TK number of peaks in in-silico spectra, needed for JCAMP-DX export
  integer :: numspec
  integer :: io_spec

  real(wp) ::  iexp (2,10000)
  real(wp) ::  mint (1000)
  real(wp) ::  tmass(10000) ! tmass contains abundances (real) for each unit mass (i)
  real(wp) ::  checksum2(50000)

  ! cthr and cthr contain ions (counts) with different charges
  ! tmax the maximum abundance factor over the whole spectrum
  real(wp) :: xx(100),tmax,r,rms,norm,chrg,cthr,cthr2
  real(wp) :: chrg1,chrg2,dum,checksum,chrg_wdh,score
  real(wp),allocatable :: rnd(:,:)

  logical  :: sel,echo,exdat,mpop,small
  logical  :: ex,ex1,ex2,ex3,ex4
  logical  :: noIso, Plasma
  
  ! fname=<qcxms.res> or result file, xname contains the mass.agr plot file
  character(len=80)  :: arg(10),line
  character(len=80)  :: xname
  character(len=:), allocatable  :: fname,fname1,fname2,fname3,fname4
  character(len=255) :: atmp

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
  chrg_wdh      = 0
  total_charge  = 1
  cthr2         = 0.01 ! only used for information, not in program


 ! edit this path name to some standard xmgrace plot file
 ! JK changed to universal location/file
 
 !xname='~/.mass_raw.agr'
 !xname='~/.mass_jay.agr'
  fname=''
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! start loop reading arguments
  do i=1,9
     arg(i)=' '
  !   call getarg(i,arg(i))
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
        chrg_wdh=xx(1)
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
  write(*,*) '* QCxMS output reader PLOTMS *'
  write(*,*) '*     v. 6.0 Sep 2021        *'
  write(*,*) '*         S. Grimme          *'
  write(*,*) '*         J. Koopman         *'
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
  if(chrg_wdh > 0)then
     write(*,'('' broadening the charges by an SD, wdth :'',f4.1)')chrg_wdh*100.
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
  if(total_charge > 0)then ! always (as information)
     write(*,*) 'The total charge of the system is ', total_charge
  endif

  ! -t (choose if unity intensities or normal)
  if(cthr >= 0)then
     write(*,'('' couting ions with charge from '',f4.1,'' to '',f4.1)')   cthr,cthr2 + total_charge
  else
     write(*,'('' counting all fragments with unity charge (frag. overview)'')')
  endif

  write(*,*) ' '
  
  ! ------------------------------------------------------------------------------------------------------!
  tmass = 0
  
  
  
  ! read file once
  i=1
  maxatm=0
  
  ! read .res file:
  ! charge(chrg2), trajectory(irun),number of collision event(icoll), numbers of (secondary-)fragments (jsec,nf), 
  ! no. of different atom types in the fragment (atm_types), atomic number (iat), amount of this atom typ (nat) <== from kk = 1 to k
  
  open ( file=fname, newunit=io_spec, status='old' )
  
  do
    if (spec == 1) then    !EI
       read(io_spec,*,iostat=io) chrg2, irun, jsec, nf, atm_types, (iat(kk),nat(kk), kk = 1, atm_types)
       if(io<0) exit !EOF
       if(io>0) stop 'Fail in read' !fail
  
    elseif(spec == 2) then !CID
       read(io_spec,*,iostat=io) chrg2, irun, jcoll, jsec, nf, atm_types, (iat(kk), nat(kk), kk = 1, atm_types)
       if(io<0) exit !EOF
       if(io>0) stop !fail

     else 
       write(*,*) 'S T O P - Something wrong in plotms - no spec'
       stop
     endif
  
     if (isec > 0 .and. isec /= jsec) cycle !check tertiary,etc fragmentation (isec)
  
     if (chrg2 > cthr) then  !default: chrg2 > 0.001
        ntot = sum(nat(1:atm_types))
        if(ntot > maxatm) maxatm = ntot  !get highest number of atoms in fragment
     endif
  
     i=i+1 !count number of single fragments with charge > chrt
  
  enddo

  n = i - 1  ! n = no. of fragments with amount maxatm of atoms
  
  write(*,*) n,' fragments with ',maxatm,' atoms max.'
  
  !close file
  close(io_spec)
  
  
  ! read it again
  write(*,*) 'Computing ...'
  
  ! contains qcxms.res as standard option
  open ( file=fname, newunit=io_spec, status='old')
  imin=100000
  i         = 1
  ial       = 0
  jal       = 0
  checksum  = 0
  checksum2 = 0


    !seed=42                                                                                     
    !call random_seed(seed)                                                
    !      I do not know why one should really randomize the results             
    !      so set the seed to 42. Now the pseudo-randomized                     
    !      computed spectrum (isotope postprocessing) is consistent every time. 
    !      Before, the score would fluctuate a lot. nrnd=50000 untouched.       
    !        Christoph Bauer Dec. 16, 2014. 

    ! initialize the random number array (efficiency)
    allocate(rnd(50000,maxatm))  
    do i = 1, 50000
      do j = 1, maxatm
        call random_number(r)
        rnd(i,j) = r
      enddo
    enddo


  do

    if (spec == 1) then    !EI
       read(io_spec,*,iostat=io) chrg2, irun, jsec, nf, atm_types, (iat(kk),nat(kk),kk=1,atm_types)
       if(io<0) exit !EOF
       if(io>0) stop !fail
  
    elseif(spec == 2) then !CID
       read(io_spec,*,iostat=io) chrg2, irun, jcoll, jsec, nf, atm_types, (iat(kk),nat(kk),kk=1,atm_types)
       if(io<0) exit !EOF
       if(io>0) stop !fail

    endif
  
    sel = .false.
    chrg = chrg2
    checksum2(irun) = checksum2(irun) + chrg2 !Check if charge sums up to 1

    if ( cthr < 0 ) chrg = total_charge

    if ( chrg > cthr ) then   ! default: yes
       sel = .true.
       ial = ial + 1              !count total amount of fragmentations
       jal(jsec) = jal(jsec) + 1  !count amount of first, sec., tert., ... fragmentations
    endif
    !----
    if ( isec > 0 .and. isec /= jsec ) sel = .false.  !Only if asked, not default
    !----

    if ( sel ) then ! default
      ntot = sum(nat(1:atm_types))
      ! all types of atoms
      kkkk = 0
      do kk = 1, atm_types
        ! all atoms of this type
        do kkk = 1, nat(kk)
          kkkk = kkkk + 1
          iat_sum(kkkk) = iat(kk)
        enddo
      enddo
      checksum = checksum + chrg   !Check total charge
    
      ! compute pattern, nsig signals at masses mass with int mint
      call isotope (ntot, iat_sum, maxatm, rnd, mass, mint, nsig, noIso)
      if(chrg_wdh > 1.0d-6) chrg = vary_energies(chrg,chrg_wdh) 

      do atm_types = 1, nsig
        tmass(mass(atm_types)) = tmass(mass(atm_types)) + mint(atm_types) * chrg
      enddo

      i = i + 1
    endif

  enddo 
  
  write(*,*) n,' (charged) fragments done.'
  write(*,*)
  write(*,*) 'checksum of charge :', checksum
  write(*,*)
  
  k=0
  do i = 1, 50000
     if (abs(checksum2(i) ) > 1.d-6) k = k + 1
  
     if ( abs(checksum2(i) ) > 1.d-6 .and. abs(checksum2(i) - total_charge ) > 1.d-3)then
        write(*,*) 'checksum error for trj', i,' chrg=',checksum2(i)
     endif
  enddo
  
  write(*,*) k,' successfull runs.'
  
  write(*,*)
  write(*,*) 'ion sources:'
  write(*,'(''% inital run'',f6.1)')100.*float(jal(1))/float(ial)
  
  do j=2,8
     write(*,'(''%        '',i1,''th'',f6.1)')j,100.*float(jal(j))/float(ial)
  enddo 
  
  write(*,*)
  
  ! read exp.
  imax=0
  imin=0
  
  inquire(file='exp.dat',exist=exdat)
  
  if(exdat) then
     write(*,*) 'Reading exp.dat ...'
     open(unit=4,file='exp.dat')
  
  ! TK (a) is edit descriptor for character strings
  20 read(4,'(a)',end=200)line
  
  if(index(line,'##MAXY=') /= 0)then
      line(7:7)=' '
      call readl(line,xx,nn)
      norm=xx(nn)/100.0d0
  endif
  
  if(index(line,'##PEAK TABLE') /= 0)then
      kk=0
  30  read(4,'(a)',end=200)line
  
     ! TK JCAMP DX for MS data has "##END="
     if(index(line,'##END') /= 0)goto 200
  
     do k=1,80
        if(line(k:k) == ',')line(k:k)=' '
     enddo
  
     call readl(line,xx,nn)
  
     do k=1,nn/2
        kk=kk+1
        iexp(1,kk)=xx(2*k-1)
        iexp(2,kk)=xx(2*k)
        if(iexp(1,kk) > imax) imax=iexp(1,kk)
        if(iexp(1,kk) < imin) imin=iexp(1,kk)
     enddo
     goto 30
  
     endif
     goto 20
  
  200 continue
     close(4)
  endif
  
  imin=max(imin,mzmin)
  j=0
  do i=mzmin,10000
     if(tmass(i) /= 0)j=i
  enddo
  imax=max(j,imax)
  
  tmax=maxval(tmass(mzmin:10000))
  
  ! counts of structures in the 100% peak
  write(*,*)'Theoretical counts in 100 % signal:',idint(tmax)
  
  if(echo)then
  write(*,*)'mass % intensity  counts   Int. exptl'
  do i=mzmin,10000
     if(tmass(i) /= 0)then
        dum=0
        do j=1,kk
           if(int(iexp(1,j)) == i)dum=iexp(2,j)
        enddo
        write(*,'(i8,F8.2,i8,F8.2)') &
        i,100.*tmass(i)/tmax,idint(tmass(i)),dum/norm
     endif
  enddo
  endif
  
  ! when the template exists
  
  inquire(file=xname,exist=ex)
  if(ex)then
     open(unit=2,file=xname)
   
  ! Original file called .mass_raw.agr
     if(small == .false.)open(unit=7,file='mass.agr')
     if(small == .true.)open(unit=7,file='small_mass.agr')
  !   open(unit=7,file='mass.agr')
  !   open(unit=9,file='intensities.txt')
  
  ! my xmgrace file mass.agr has nagrfile lines
     write(*,*)
     write(*,*) 'Writing mass.agr ...'
     do i=1,nagrfile
        read (2,'(a)')line
        if(small == .false.)then 
           if(index(line,'world 10') /= 0)then
              write(line,'(''@    world '',i3,'', -105,'',i3,'', 100'')')imin,imax+5
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
              write(line,'(''@    world '',i3,'', 0,'',i3,'', 100'')')imin,imax+5
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
        write(7,'(a)')line
     enddo
  
  ! write only masses > mzmin
     do i=mzmin,10000
        if(tmass(i) /= 0)then
           write(7,*) i,100.*tmass(i)/tmax
!           write(9,*) i,100.*tmass(i)/tmax   !ce50values
        endif
     enddo
     write(7,*)'&'
  
!     close(9)
  
  ! TK establish again if exdat (experimental JCAMPDX) exists
  ! TK @target G0.S1 means upper theory spectrum in current implementation
  ! TK @target G0.S2 means lower experimental spectrum
  ! TK iexp(1,i) - contains the exp masses and  iexp(2,i) contains the abundance
  ! TK kk contains number of experimental spectra
  
     if(exdat)then
       write(7,*)'@target G0.S2'
       write(7,*)'@type bar'
       do i=1,kk
          write(7,*) iexp(1,i),-iexp(2,i)/norm
       enddo
       write(7,*)'&'
     endif
  endif
  
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
  ! TK here we export the spectrum to JCAMP-DX, subroutine would be better for
  ! TK allowing later code merges or inclusion of new and missing elements such as
  ! TK currently unit mass only
  
  write(*,*) 'Writing JCAMP file: result.jdx ...'
  write(*,*)
  ! TK open file as JCAMP-DX MS result file, open new or replace
  open(unit=11, file='result.jdx', STATUS="REPLACE")
  
  ! TK minimum JCAMP-DX header
  write(11,"(A)")'##TITLE=Theoretical in-silico spectrum (QCxMS)'
  write(11,"(A)")'##JCAMP-DX=4.24'
  write(11,"(A)")'##DATA TYPE=MASS SPECTRUM'
  
  !write(11,*)'##MOLFORM=C13 H22 O3 Si2'
  !write(11,*)'##MW=282'
  write(11,"(A)")'##XUNITS=M/Z'
  write(11,"(A)")'##YUNITS=RELATIVE INTENSITY'
  !write(11,*)'##XFACTOR=1'
  !write(11,*)'##YFACTOR=1'
  !write(11,*)'##FIRSTX=15'
  !write(11,*)'##LASTX=281'
  !write(11,*)'##FIRSTY=20'
  !write(11,*)'##MAXX=281'
  !write(11,*)'##MINX=15'
  !write(11,*)'##MAXY=9999'
  !write(11,*)'##MINY=10'
  
  ! TK calculate number of in-silico spectra
  ! TK new: numspec = idint(tmax) ** not related to number ofr spectral peaks
  numspec = 0
  do i=10,10000
    if(tmass(i) /= 0)then
      numspec = numspec + 1
    endif
  enddo
  
  
  write(11,'(A, I0)')'##NPOINTS=' ,numspec
  
  ! TK ##XYDATA=(XY..XY) 1 designates one line per m/z abd pair
  ! TK separated by comma
  ! TK 1.500000e+001 ,5.000000e+000
  
  write(11,"(A)")'##PEAK TABLE=(XY..XY) 1'
  !15,20 26,10 27,20 29,50
  
  
  do i=10,10000
       if(tmass(i) /= 0)then
         ! unit mass to exact mass
         ! write(11,"(I4, I8)") i, int(100.*tmass(i)/tmax)
         write(11,"(F12.6, I8)") real(i), int(100.*tmass(i)/tmax)
       endif
  enddo
  
  write(11,"(A)")'##END='
  
  ! TK close the JCAMP-DX file
  close(11)
  
 ! compute deviation exp-theor.
  
  if(exdat)then

    rms=0
    nn=0
    kkk=0
    kkkk=0

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
  close(2)
  close(3)
  close(7)
  
  ! compute mass spectral match score                                                                 
  if(exdat)then                                                                                 
  score=0.                                                                                      
  ! bring iexp to the right order of magnitude                                                        
    do i=1,10000                                                                                  
      iexp(2,i)=iexp(2,i)/norm                                                                    
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

  integer  :: i,j,k

  real(wp) :: tmass(10000)
  real(wp) :: int_exp(10000)
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
  do i=1,10000
     do j=1,10000
         norm_experiment = norm_experiment + (-1)**(i+j)*iexp(2,j)*iexp(2,i)*exp(-alpha*(i-j)**2)
         norm_calculation = norm_calculation + (-1)**(i+j)*tmass(j)*tmass(i)*exp(-alpha*(i-j)**2)
         overlap = overlap + (-1)**(i+j)*iexp(2,j)*tmass(i)*exp(-alpha*(i-j)**2)
     enddo
  enddo
  
  score = overlap/sqrt(norm_experiment*norm_calculation)
  
  end subroutine match2
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine match(tmass,iexp,score,tmax)
  use xtb_mctc_accuracy, only: wp
  implicit none

  integer :: i,j,k
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
          w(2,j) = j**2 * ( 100.0_wp * tmass(j) / tmax )**0.6_wp
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
  
  dot  =0.0_wp
  sum1 =0.0_wp
  sum2 =0.0_wp
  sum3 =0.0_wp

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
