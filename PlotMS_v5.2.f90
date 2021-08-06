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
!      - CE50 values    -- J.Koopman
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
  integer :: nrnd,maxatm,imax,imin,nagrfile,spec,io
  integer :: iat  (10000)
  integer :: nat  (10000)
  integer :: idum (10000)
  integer :: mass (1000)! maximum mass = 10000
  integer :: isec,jcoll,jsec,ial,jal(0:10),nf,irun,seed,mzmin
  ! TK number of peaks in in-silico spectra, needed for JCAMP-DX export
  integer :: numspec

  real(wp) ::  iexp (2,10000)
  real(wp) ::  mint (1000)
  real(wp) ::  tmass(10000) ! tmass contains abundances (real) for each unit mass (i)
  real(wp) ::  ce_val,sum_ce,CE50
  real(wp) ::  checksum2(50000)
  real(wp), allocatable :: rnd(:,:)

  ! cthr and cthr contain ions (counts) with different charges
  ! tmax the maximum abundance factor over the whole spectrum
  real(wp) :: xx(100),tmax,r,rms,norm,chrg,cthr,cthr2
  real(wp) :: chrg1,chrg2,dum,checksum,chrg_wdh,score
  real(wp) :: snorm

  logical  :: sel,echo,exdat,mpop,small
  logical  :: ex,ex1,ex2,ex3,ex4
  logical  :: noIso
  
  ! fname=<qcxms.res> or result file, xname contains the mass.agr plot file
  character(len=80)  :: arg(10),line
  character(len=80)  :: xname
  character(len=80)  :: fname,fname1,fname2,fname3,fname4
  character(len=255) :: atmp
 
  iexp =0
  mpop =.false.
  echo =.false.
  small =.false.
  isec =0
  cthr =1.d-3
  chrg_wdh   =0
  norm =1.0
  cthr2=1.01
  nagrfile=410
  mzmin = 10
  spec = 0
  noIso = .false.
  
  ! edit this path name to some standard xmgrace plot file
  ! JK changed to universal location/file
  
  !xname='~/.mass_raw.agr'
  !xname='~/.mass_jay.agr'
  fname=''
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! start loop reading arguments
  do i=1,9
     arg(i)=' '
     call getarg(i,arg(i))
  enddo
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  ! start loop processing aruments
  
  !  comand line parameters
  !  -a count charges from -1000 (cthr) to cthr2
  !  -p print spectra "mass % intensity  counts   Int. exptl" to stdout;
  !     with "Int. exptl" (experimental) taken from exp.dat but not all exp peaks are exported
  !     if no theoretical counterpart exists
  !  -f filename or  -f <name_of_res_file>
  !  -t couting ions with charge x to y (give the value, e.g. "-t 1 2" for charge 1 to 2)
  !  -w broadening the charges by an SD 
  !  -s Take only secondary and tertiary fragmentations (give the value, e.g. "-s 2" for secondary)
  !  -m set minimum value for m/z, so 100% value will be calc. for higher values (CID)
  !  -i calculate NO isotope pattern
  
  do i=1,9
     if(index(arg(i),'-a') /= 0)  cthr=-1000.
     if(index(arg(i),'-p') /= 0)  echo=.true.
     if(index(arg(i),'-f') /= 0)  fname=arg(i+1)
     if(index(arg(i),'-i') /= 0)  noIso=.true.
     if(index(arg(i),'-t') /= 0) then
        call readl(arg(i+1),xx,nn)
        cthr=xx(1)
     endif
     if(index(arg(i),'-w') /= 0) then
        call readl(arg(i+1),xx,nn)
        chrg_wdh=xx(1)
     endif
     if(index(arg(i),'-s') /= 0) then
        call readl(arg(i+1),xx,nn)
        isec=int(xx(1))
     endif
     if(index(arg(i),'-m') /= 0)then
        call readl(arg(i+1),xx,nn)
        mzmin=int(xx(1))
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
  
  xname='~/.mass_raw.agr'
  
  ! fname contains the results from each calculation or the temporary result tmpqcxms.res
  ! xname contains the xmgrace plot file
  
  if ( fname == '' ) then
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
  write(*,*) '*     v. 5.1 Aug 2021        *'
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
  ! -t
  if(cthr >= 0)then
     write(*,'('' couting ions with charge from '',f4.1,'' to '',f4.1)')   cthr,cthr2
  else
     write(*,'('' counting all fragments with unity charge (frag. overview)'')')
  endif
  
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
  write(*,*) ' '
  
  ! ------------------------------------------------------------------------------------------------------!
  tmass = 0
  sum_ce = 0
  
  seed=42                                                                                     
  call random_seed(seed)                                                
  !      I do not know why one should really randomize the results             
  !      so set the seed to 42. Now the pseudo-randomized                     
  !      computed spectrum (isotope postprocessing) is consistent every time. 
  !      Before, the score would fluctuate a lot. nrnd=50000 untouched.       
  !        Christoph Bauer Dec. 16, 2014. 
  
  
  ! read file once
  i=1
  maxatm=0
  
  ! read .res file:
  ! charge(chrg2), trajectory(irun),number of collision event(icoll), numbers of (secondary-)fragments (jsec,nf), 
  ! types of atoms in fragment (k), atomic number (iat), amount of this atom typ (nat) <== from kk = 1 to k
  
  open(unit=1,file=fname)
  
  do
    if (spec == 1) then    !EI
       read(1,*,iostat=io)chrg2,irun,jsec,nf,k,(iat(kk),nat(kk),kk=1,k)
       if(io<0) exit !EOF
       if(io>0) stop 'Fail in read' !fail
  
    elseif(spec == 2) then !CID
       read(1,*,iostat=io)chrg2,irun,jcoll,jsec,nf,k,(iat(kk),nat(kk),kk=1,k)
       if(io<0) exit !EOF
       if(io>0) stop !fail
  
     else !not yet decided what to do with tmpqcxms
       write(*,*) 'S T O P - Something wrong in plotms - no spec'
       stop
     endif
  
     if(isec > 0.and.isec /= jsec) cycle !check tertiary,etc fragmentation (isec)
  
     if(chrg2 > cthr) then  !default: chrg2 > 0.001
        ntot=sum(nat(1:k))
        if(ntot > maxatm)maxatm=ntot  !get highest number of atoms in fragment
     endif
  
     i=i+1 !count number of single fragments with charge > chrt
  
  enddo

  n=i-1  ! n = no. of fragments with amount maxatm of atoms
  
  write(*,*) n,' fragments with ',maxatm,' atoms max.'
  
  !close file
  close(1)
  
  ! initialize the random number array (efficiency)
  nrnd=50000
  allocate(rnd(maxatm,nrnd))
  do i=1,nrnd
     do j=1,maxatm
        call random_number(r)
        rnd(j,i)=r
     enddo
  enddo
  
  ! read it again
  write(*,*) 'Computing ...'
  
  ! contains qcxms.res as standard option
  open(unit=1,file=fname)
  imin=100000
  i=1
  ial=0
  jal=0
  checksum =0
  checksum2=0
 
  do
    if (spec == 1) then    !EI
       read(1,*,iostat=io)chrg2,irun,jsec,nf,k,(iat(kk),nat(kk),kk=1,k)
       if(io<0) exit !EOF
       if(io>0) stop !fail
  
    elseif(spec == 2) then !CID
       read(1,*,iostat=io)chrg2,irun,jcoll,jsec,nf,k,(iat(kk),nat(kk),kk=1,k)
       if(io<0) exit !EOF
       if(io>0) stop !fail
    endif
  
  sel = .false.
  chrg = chrg2
  checksum2(irun) = checksum2(irun) + chrg2 !Check if charge sums up to 1

  if(cthr < 0)chrg=1.0                  

  if(chrg > cthr)then   ! default: yes
     sel=.true.
     ial=ial+1              !count total amount of fragmentations
     jal(jsec)=jal(jsec)+1  !count amount of first, sec., tert., ... fragmentations
  endif
  !----
  if(isec > 0.and.isec /= jsec) sel=.false.  !Only if asked, not default
  !----

  if(sel)then
    ntot=sum(nat(1:k))
  ! all types of atoms
     kkkk=0
     do kk=1,k
  ! all atoms of this type
        do kkk=1,nat(kk)
           kkkk=kkkk+1
           idum(kkkk)=iat(kk)
        enddo
      enddo
      checksum=checksum+chrg   !Check total charge
  
  ! compute pattern, nsig signals at masses mass with int mint
      call isotope(ntot,idum,maxatm,nrnd,rnd,mass,mint,nsig,noIso)
      if(chrg_wdh > 1.0d-6) chrg = vary_energies(chrg,chrg_wdh) 
         do k=1,nsig
           tmass(mass(k))=tmass(mass(k))+mint(k)*chrg
         enddo
         i=i+1
      endif
  enddo 
  
  write(*,*) n,' (charged) fragments done.'
  write(*,*)
  write(*,*) 'checksum of charge :',checksum
  write(*,*)
  
  k=0
  do i=1,50000
     if(abs(checksum2(i)) > 1.d-6)k=k+1
  
     if(abs(checksum2(i)) > 1.d-6.and.abs(checksum2(i)-1.) > 1.d-3)then
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
  do i=1,kk
     k=iexp(1,i)
     if(iexp(2,i)/norm > 5.0)then
        r=100.*tmass(k)/tmax-iexp(2,i)/norm
        kkk=kkk+1
        if(100.*tmass(k)/tmax > 2.5)kkkk=kkkk+1
        rms=rms+abs(r)
     endif
  enddo
  write(*,*)'MAD(exptl./theor.) = ',rms/kkk
  write(*,*)'# exptl. > 5 %     = ',kkk
  write(*,*)'% correctly found  = ',100*float(kkkk)/float(kkk)
  endif
  
  close(1)
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
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  ! compute the isotopic mass distribution for
  ! given chemical formula (nat atoms with OZ it())
  ! returns nsig masses in mass() with probability
  ! mint()
  ! this is a quick and dirty MC algorithm and its run-time depends
  ! critically on the number of trials nrnd
  
  ! TK See: ATOMIC WEIGHTS OF THE ELEMENTS: REVIEW 2000 (IUPAC Technical Report)
  ! TK Pure Appl. Chem., Vol. 75, No. 6, pp. 683–800, 2003.
  
  subroutine isotope(nat,it,maxatm,nrnd,rnd,mass,mint,nsig,no_isotopes)
  use xtb_mctc_accuracy, only: wp
  implicit none

  integer :: nat,it(*),nsig,nrnd,maxatm
  integer :: mass(*)
  integer :: niso(200)
  integer :: nmass(10000)
  integer :: n,i,j,iso,imass,isum,k,iti

  real(wp) ::  mint(*)
  real(wp) ::  rnd(maxatm,nrnd)
  real(wp) :: r,xm
  real(wp) ::  prob(200,10),massiso(200,10),p1,p2,x,xmass

  logical  :: no_isotopes

  
  niso=0
  prob=0
  massiso=0
  
   !  1 H  (Hydrogen)
        niso(1)           = 2
        prob(1,1)         = 0.0115_wp
        prob(1,2)         = 99.9885_wp
        massiso(1,1)      = 2.014101_wp
        massiso(1,2)      = 1.007825_wp

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
  
  ! mass currently loops to element 83 (Bismuth)
  do i = 1,83
     xm = 0
     do j = 1, niso(i)
        xm=xm+prob(i,j)
     enddo
     if ( niso(i) > 0 .and. abs(xm-1.) > 0.01) stop 'internal isotope error 1'
  enddo
  
  do i = 1, nat
     if ( it(i) > 100 ) then
        niso   (it(i))   = 1
        prob   (it(i),1) = 1.0_wp
        massiso(it(i),1) = it(i) - 100.0_wp
     endif
  enddo
 
  ! if number of isotopes == 0 <- often errors if wrong .res file is read (is (hopefully) fixed)
  do i=1,nat
     if(niso(it(i)) == 0) stop 'internal isotope error 2'
  enddo

  !Calculate NO isotope pattern
  if (no_isotopes == .true. ) then 
    nmass = 0
    !do n=1,nrnd ! I commented this out, because I see no effects on the spectrum  
    xmass=0
    do i = 1, nat
      iti = it(i)
      p1 = 0.0
      do iso = 1, niso(iti)
         if ( prob(iti,iso) > p1 ) then
           x = massiso(iti,iso)
         endif
         p1 = prob(iti,iso)
      enddo
      xmass = xmass + x
    enddo
    imass = nint(xmass)                  ! here it gets to integers - changed to nearest int nint
    nmass(imass) = nmass(imass)+1
    !enddo

  ! Calculate isotope pattern
  else                      
    nmass = 0
    do n = 1, nrnd
      xmass = 0
      do i = 1, nat
         iti = it(i)
         r = rnd(i,n)
         p1 = 0.0
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
      imass = nint(xmass)                  ! here it gets to integers - changed to nearest int nint
      nmass(imass) = nmass(imass) + 1
    enddo
  endif
  
  isum = sum(nmass)
  k = 0
  do i = 1, 2000
     if ( nmass(i) > 0 ) then
        k = k + 1
        mass(k) = i
        mint(k) = float(nmass(i)) / float(isum)
     endif
  enddo
  
  nsig = k

end subroutine isotope

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine match2(tmass,iexp,score,alpha)
  implicit none
  real*8 tmass(10000)
  real*8 int_exp(10000)
  real*8 iexp(2,10000)
  real*8 score
  real*8 alpha
  integer i,j,k
  !     needed storage variables
  real*8 norm_experiment
  real*8 norm_calculation
  real*8 overlap
  
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
  
  ! ardous loop external file calcfr.f90
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
  
  return 
  end subroutine calcfr
