! ============================================================================
! Name          : plotms.f90
! Author        : S. Grimme
! Co-authors    : J.Koopman, C.Bauer
! Contributions : T. Kind (FiehnLab 2013)
! Version       : 5.0  (Jun 11 2021)
! Copyright     :
! Description   : plot mass spectra from QCxMS
! ============================================================================

!====================================================================
!     Original PlotMS Program for extraction of mass spectra
!     from quantum mechanical simulations used within QCxMS
!
!                      *********************************************
!                      *                S. Grimme                  *
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

!!! 5.0: Changing name from QCEIMS to QCxMS !!! 


program plotms
!  use xtb_mctc_accuracy, only:wp
  !treat i,j,k,l,m,n as integers and all other cariables as real
  implicit none
  integer n,i,j,k,kk,kkk,kkkk,nn,ntot,nsig
  integer nrnd,maxatm,imax,imin,nagrfile,spec,io
  ! maximum mass = 10000
  real*8  iexp (2,10000)
  integer iat  (10000)
  integer nat  (10000)
  integer idum (10000)
  integer mass (1000)
  integer isec,jcoll,jsec,ial,jal(0:10),nf,irun,seed,mzmin
  real*8  mint (1000)
  ! tmass contains abundances (real) for each unit mass (i)
  real*8  tmass(10000)
  real*8  ce_val,sum_ce,CE50
  real*8  checksum2(50000)
  real*8, allocatable :: rnd(:,:)
  ! cthr and cthr contain ions (counts) with different charges
  ! tmax the maximum abundance factor over the whole spectrum
  real*8 xx(100),tmax,r,rms,norm,chrg,cthr,cthr2
  real*8 chrg1,chrg2,dum,checksum,cw,score
  real   snorm
  logical ex,sel,echo,exdat,mpop,small
  logical noIso
  
  ! fname=<qcxms.res> or result file, xname contains the mass.agr plot file
  character*80 arg(10),line,fname,xname
  character*255 atmp
  
  ! TK number of peaks in in-silico spectra, needed for JCAMP-DX export
  integer numspec
  
  iexp =0
  mpop =.false.
  echo =.false.
  small =.false.
  isec =0
  cthr =1.d-3
  cw   =0
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
     if(index(arg(i),'-a').ne.0)  cthr=-1000.
     if(index(arg(i),'-p').ne.0)  echo=.true.
     if(index(arg(i),'-f').ne.0)  fname=arg(i+1)
     if(index(arg(i),'-i').ne.0)  noIso=.true.
     if(index(arg(i),'-t').ne.0) then
        call readl(arg(i+1),xx,nn)
        cthr=xx(1)
     endif
     if(index(arg(i),'-w').ne.0) then
        call readl(arg(i+1),xx,nn)
        cw=xx(1)
     endif
     if(index(arg(i),'-s').ne.0) then
        call readl(arg(i+1),xx,nn)
        isec=int(xx(1))
     endif
     if(index(arg(i),'-m').ne.0)then
        call readl(arg(i+1),xx,nn)
        mzmin=int(xx(1))
     endif
  enddo

  
!  if(index(arg(1),'-k').ne.0) then
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
  
  if(fname.eq.'')then
    fname='qcxms.res'
    inquire(file=fname,exist=ex)
    if(ex) spec = 1
    if(.not.ex)then
      fname='qcxms_cid.res'
      spec = 2
    endif
  
    inquire(file=fname,exist=ex)
    if(.not.ex) fname='tmpqcxms.res' ! not solved yet
  endif
  
  inquire(file=fname,exist=ex)
  if(.not.ex) stop 'res file does not exist'
  
  write(*,*)
  write(*,*) '*******************************'
  write(*,*) '* QCxMS output reader PLOTMS *'
  write(*,*) '*      V 5.0 Jun 2021        *'
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
  if(cthr.ge.0)then
     write(*,'('' couting ions with charge from '',f4.1,'' to '',f4.1)')   cthr,cthr2
  else
     write(*,'('' counting all fragments with unity charge (frag. overview)'')')
  endif
  
  write(*,*) ' '
  ! -w
  if(cw.gt.0)then
     write(*,'('' broadening the charges by an SD, wdth :'',f4.1)')cw*100.
  endif
  
  ! -s
  if(isec.ne.0)then
     write(*,*) 'Taking only secondary, tertiary ... fragmentations ',isec
  endif
  
  ! -m
  if(mzmin.gt.10)then
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
    if (spec.eq.1) then    !EI
       read(1,*,iostat=io)chrg2,irun,jsec,nf,k,(iat(kk),nat(kk),kk=1,k)
       if(io<0) exit !EOF
       if(io>0) stop !fail
  
    elseif(spec.eq.2) then !CID
       read(1,*,iostat=io)chrg2,irun,jcoll,jsec,nf,k,(iat(kk),nat(kk),kk=1,k)
       if(io<0) exit !EOF
       if(io>0) stop !fail
  
     else !not yet decided what to do with tmpqcxms
       write(*,*) 'S T O P - Something wrong in plotms - no spec'
       stop
     endif
  
     if(isec.gt.0.and.isec.ne.jsec) cycle !check tertiary,etc fragmentation (isec)
  
     if(chrg2.gt.cthr) then  !default: chrg2 > 0.001
        ntot=sum(nat(1:k))
        if(ntot.gt.maxatm)maxatm=ntot  !get highest number of atoms in fragment
     endif
  
     i=i+1 !count number of single fragments with charge > chrt
  
  enddo

  n=i-1
  ! n = no. of fragments with amount maxatm of atoms
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
    if (spec.eq.1) then    !EI
       read(1,*,iostat=io)chrg2,irun,jsec,nf,k,(iat(kk),nat(kk),kk=1,k)
       if(io<0) exit !EOF
       if(io>0) stop !fail
  
    elseif(spec.eq.2) then !CID
       read(1,*,iostat=io)chrg2,irun,jcoll,jsec,nf,k,(iat(kk),nat(kk),kk=1,k)
       if(io<0) exit !EOF
       if(io>0) stop !fail
    endif
  
  sel=.false.
  chrg=chrg2
  checksum2(irun)=checksum2(irun)+chrg2 !Check if charge sums up to 1

  if(cthr.lt.0)chrg=1.0                  

  if(chrg.gt.cthr)then   ! default: yes
     sel=.true.
     ial=ial+1              !count total amount of fragmentations
     jal(jsec)=jal(jsec)+1  !count amount of first, sec., tert., ... fragmentations
  endif
  !----
  if(isec.gt.0.and.isec.ne.jsec) sel=.false.  !Only if asked, not default
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
      if(cw.gt.1.d-6)chrg=chrg+cw*chrg*snorm()
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
     if(abs(checksum2(i)).gt.1.d-6)k=k+1
  
     if(abs(checksum2(i)).gt.1.d-6.and.abs(checksum2(i)-1.).gt.1.d-3)then
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
  
  if(index(line,'##MAXY=').ne.0)then
      line(7:7)=' '
      call readl(line,xx,nn)
      norm=xx(nn)/100.0d0
  endif
  
  if(index(line,'##PEAK TABLE').ne.0)then
      kk=0
  30  read(4,'(a)',end=200)line
  
     ! TK JCAMP DX for MS data has "##END="
     if(index(line,'##END').ne.0)goto 200
  
     do k=1,80
        if(line(k:k).eq.',')line(k:k)=' '
     enddo
  
     call readl(line,xx,nn)
  
     do k=1,nn/2
        kk=kk+1
        iexp(1,kk)=xx(2*k-1)
        iexp(2,kk)=xx(2*k)
        if(iexp(1,kk).gt.imax) imax=iexp(1,kk)
        if(iexp(1,kk).lt.imin) imin=iexp(1,kk)
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
     if(tmass(i).ne.0)j=i
  enddo
  imax=max(j,imax)
  
  tmax=maxval(tmass(mzmin:10000))
  
  ! counts of structures in the 100% peak
  write(*,*)'Theoretical counts in 100 % signal:',idint(tmax)
  
  if(echo)then
  write(*,*)'mass % intensity  counts   Int. exptl'
  do i=mzmin,10000
     if(tmass(i).ne.0)then
        dum=0
        do j=1,kk
           if(int(iexp(1,j)).eq.i)dum=iexp(2,j)
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
     if(small.eq..false.)open(unit=7,file='mass.agr')
     if(small.eq..true.)open(unit=7,file='small_mass.agr')
  !   open(unit=7,file='mass.agr')
  !   open(unit=9,file='intensities.txt')
  
  ! my xmgrace file mass.agr has nagrfile lines
     write(*,*)
     write(*,*) 'Writing mass.agr ...'
     do i=1,nagrfile
        read (2,'(a)')line
        if(small.eq..false.)then 
           if(index(line,'world 10').ne.0)then
              write(line,'(''@    world '',i3,'', -105,'',i3,'', 100'')')imin,imax+5
           endif
           if(index(line,'xaxis  tick major').ne.0)then
              line='@    xaxis  tick major 20'
                 if(imax.gt.200) line='@    xaxis  tick major 50' 
           endif
  
           if(index(line,'@    s1 symbol size').ne.0)then
              line='@    s1 symbol size 0.160000'
              if(imax.gt.200) line='@    s1 symbol size 0.100000'
           endif
  
           if(index(line,'@    s2 symbol size').ne.0)then
              line='@    s2 symbol size 0.160000'
              if(imax.gt.200) line='@    s2 symbol size 0.100000'
           endif
        !kleines Format
        else
           if(index(line,'world 10').ne.0)then
              write(line,'(''@    world '',i3,'', 0,'',i3,'', 100'')')imin,imax+5
           endif  
  !         if(index(line,'xaxis  tick major').ne.0)then
  !            line='@    xaxis  tick major 10'
  !!               if(imax.gt.200) line='@    xaxis  tick major 50' 
  !         endif
           
           if(index(line,'@    s1 symbol size').ne.0)then
              line='@    s1 symbol size 0.320000'
  !            if(imax.gt.200) line='@    s1 symbol size 0.100000'
           endif
        endif
        write(7,'(a)')line
     enddo
  
  ! write only masses > mzmin
     do i=mzmin,10000
        if(tmass(i).ne.0)then
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
      if(tmass(i).ne.0)then
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
       if(tmass(i).ne.0)then
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
     if(iexp(2,i)/norm.gt.5.0)then
        r=100.*tmass(k)/tmax-iexp(2,i)/norm
        kkk=kkk+1
        if(100.*tmass(k)/tmax.gt.2.5)kkkk=kkkk+1
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
   
        end
  
  ! test code
  !     real*8 mint(100)
  !     integer mass(100),it(15),nrnd
  !     real*4 rnd(15,50000)
  
  !     nrnd=50000
  !     do i=1,nrnd
  !        do j=1,15
  !           call random_number(r)
  !           rnd(j,i)=r
  !       enddo
  !     enddo
  
  !     do i=1,6
  !        it(i)=6
  !     enddo
  !     do i=7,12
  !        it(i)=1
  !     enddo
  !     it(13)=16
  !     it(14)=7
  !     it(15)=17
  !     call isotope(15,it,nrnd,rnd,mass,mint,nsig)
  !     do i=1,nsig
  !        write(*,*) mass(i),mint(i)*100.
  !     enddo
  !     end
  
  ! compute the isotopic mass distribution for
  ! given chemical formula (nat atoms with OZ it())
  ! returns nsig masses in mass() with probability
  ! mint()
  ! this is a quick and dirty MC algorithm and its run-time depends
  ! critically on the number of trials nrnd
  ! TK increased mass acuracy values, currently unit masses                          <-- Unit masses are correct, otherwise it doesn't count right, because it makes integers in isotop subroutine. Change maybe?!
  
  ! TK See: ATOMIC WEIGHTS OF THE ELEMENTS: REVIEW 2000 (IUPAC Technical Report)
  ! TK Pure Appl. Chem., Vol. 75, No. 6, pp. 683–800, 2003.
  
  subroutine isotope(nat,it,maxatm,nrnd,rnd,mass,mint,nsig,no_isotopes)
  implicit none

  integer nat,it(*),nsig,nrnd,maxatm
  integer mass(*)
  integer niso(200)
  integer nmass(10000)
  integer n,i,j,iso,imass,isum,k,iti

  real*8  mint(*)
  real*8  rnd(maxatm,nrnd)
  real*8 r,xm
  real*8  prob(200,10),massiso(200,10),p1,p2,x,xmass

  logical no_isotopes

  
  niso=0
  prob=0
  massiso=0
  
   !  1 H  (Hydrogen)
        niso(1)=2
        prob(1,1)=0.0115
        prob(1,2)=99.9885
        massiso(1,1)=2.014101 
        massiso(1,2)=1.007825
   !  6 C  (Carbon)
        niso(6)=2
        prob(6,1)=1.07
        prob(6,2)=98.93
        massiso(6,1)=13.003354
        massiso(6,2)=12.00
  
   !  7 N  (Nitrogen)
        niso(7)=2
        prob(7,1)=0.368
        prob(7,2)=99.632
        massiso(7,1)=15.000108
        massiso(7,2)=14.003074
  
   !  8 O  (Oxygen)
        niso(8)=3
        prob(8,1)=0.038
        prob(8,2)=0.205
        prob(8,3)=99.757
        massiso(8,1)=16.999131
        massiso(8,2)=17.999159
        massiso(8,3)=15.994914
  
   !  9 F  (Fluorine)
        niso(9)=1
        prob(9,1)=100.
        massiso(9,1)=18.998403
  
  !  13 Al (Aluminium)                                                                                
          niso(13)=1                                                                                  
          prob(13,1)=100.0                                                                            
          massiso(13,1)=26.981538
  
   ! 14 Si (Silicon)
        niso(14)=3
        prob(14,1)=92.223
        prob(14,2)=4.685
        prob(14,3)=3.092
        massiso(14,1)=27.976926
        massiso(14,2)=28.976494
        massiso(14,3)=29.973770
  
   ! 15 P  (Phosphorus)
        niso(15)=1
        prob(15,1)=100.
        massiso(15,1)=30.973761
  
   ! 16 S  (Sulfur)
        niso(16)=4
        prob(16,1)=0.02
        prob(16,2)=0.76
        prob(16,3)=4.29
        prob(16,4)=94.93
        massiso(16,1)=35.967080
        massiso(16,2)=32.971458
        massiso(16,3)=33.967867
        massiso(16,4)=31.972071
  
   ! 17 Cl (Chlorine) 
        niso(17)=2
        prob(17,1)=75.76
        prob(17,2)=24.24
        massiso(17,1)=34.968852 !35.0 !
        massiso(17,2)=36.965902 !37.0 !
  
   ! 22 Ti (Titanium) 
        niso(22)=5           
        prob(22,1)=8.0       
        prob(22,2)=7.3       
        prob(22,3)=73.8      
        prob(22,4)=5.5       
        prob(22,5)=5.4       
        massiso(22,1)=45.952632
        massiso(22,2)=46.951763
        massiso(22,3)=47.947946
        massiso(22,4)=48.947870
        massiso(22,5)=49.944791
                            
   ! 24 Cr (Chromium)             
        niso(24)=4           
        prob(24,1)=4.35      
        prob(24,2)=83.79     
        prob(24,3)=9.50      
        prob(24,4)=2.37      
        massiso(24,1)=49.946044
        massiso(24,2)=51.940507
        massiso(24,3)=52.940649
        massiso(24,4)=53.938880
                            
   ! 25 Mn (Manganese)
        niso(25)=1           
        prob(25,1)=100.0     
        massiso(25,1)=54.938045
  
  !  26 Fe (Iron)
        niso(26)=4
        prob(26,1)=5.845
        prob(26,2)=91.754
        prob(26,3)=2.119
        prob(26,4)=0.2
        massiso(26,1)=53.939
        massiso(26,2)=55.934
        massiso(26,3)=56.935
        massiso(26,4)=57.933
  
   ! 27 Co (Cobalt)
        niso(27)=1         
        prob(27,1)=100.000 
        massiso(27,1)=58.933195 
                          
  !  28 Ni (Nickel)
        niso(28)=5         
        prob(28,1)=68.08   
        prob(28,2)=26.22   
        prob(28,3)=1.14    
        prob(28,4)=3.63    
        prob(28,5)=0.93    
        massiso(28,1)=57.935343
        massiso(28,2)=59.930786
        massiso(28,3)=60.931056
        massiso(28,4)=61.928345
        massiso(28,5)=63.927966
                          
   ! 29 Cu (Copper)             
        niso(29)=2         
        prob(29,1)=69.15   
        prob(29,2)=30.85   
        massiso(29,1)=62.929597
        massiso(29,2)=64.927789
                          
   ! 30 Zn (Zinc)
        niso(30)=5         
        prob(30,1)=48.6    
        prob(30,2)=27.9    
        prob(30,3)=4.1     
        prob(30,4)=18.8    
        prob(30,5)=0.6     
        massiso(30,1)=63.929142
        massiso(30,2)=65.926033
        massiso(30,3)=66.927127
        massiso(30,4)=67.924884
        massiso(30,5)=69.925319
  
   ! 32 Ge (Germanium)
        niso(32)=5          
        prob(32,1)=20.38    
        prob(32,2)=27.31    
        prob(32,3)=7.76     
        prob(32,4)=36.72    
        prob(32,5)=7.83     
        massiso(32,1)=69.924247 
        massiso(32,2)=71.922076  
        massiso(32,3)=72.923459  
        massiso(32,4)=73.921178  
        massiso(32,5)=75.921402  
                              
   ! 33 As (Arsenic)
        niso(33)=1          
        prob(33,1)=100.0    
        massiso(33,1)=74.921596
                              
   ! 34 Se (Selenium)
        niso(34)=6          
        prob(34,1)=0.89     
        prob(34,2)=9.37     
        prob(34,3)=7.63     
        prob(34,4)=23.77    
        prob(34,5)=49.61    
        prob(34,6)=8.73     
        massiso(34,1)=73.922476
        massiso(34,2)=75.919213
        massiso(34,3)=76.919914
        massiso(34,4)=77.917309
        massiso(34,5)=79.916521
        massiso(34,6)=81.916699
   
   ! 35 Br (Bromine)
        niso(35)=2
        prob(35,1)=50.69
        prob(35,2)=49.31
        massiso(35,1)=78.91833
        massiso(35,2)=80.91629
  
   ! 40 Zr (Zirconium)
          niso(40)=5          
          prob(40,1)=51.45    
          prob(40,2)=11.22    
          prob(40,3)=17.15    
          prob(40,4)=17.38    
          prob(40,5)=2.80     
          massiso(40,1)=89.904704
          massiso(40,2)=90.905646
          massiso(40,3)=91.905040
          massiso(40,4)=93.906315
          massiso(40,5)=95.908273
                              
   ! 44 Ru (Ruthenium)
          niso(44)=7          
          prob(44,1)=5.554    
          prob(44,2)=1.873    
          prob(44,3)=12.761   
          prob(44,4)=12.607   
          prob(44,5)=17.062   
          prob(44,6)=31.551   
          prob(44,7)=18.623   
          massiso(44,1)=95.907598
          massiso(44,2)=97.905287
          massiso(44,3)=98.905939
          massiso(44,4)=99.904219
          massiso(44,5)=100.905582
          massiso(44,6)=101.906323
          massiso(44,7)=103.905433
                              
   ! 46 Pd (Palladium)
          niso(46)=6          
          prob(46,1)=1.02     
          prob(46,2)=11.15    
          prob(46,3)=22.34    
          prob(46,4)=27.33    
          prob(46,5)=26.47    
          prob(46,6)=11.73    
          massiso(46,1)=101.905609
          massiso(46,2)=103.904036
          massiso(46,3)=104.905085
          massiso(46,4)=105.903486
          massiso(46,5)=107.903892
          massiso(46,6)=109.905153
  
   ! 50 Sn (Tin) 
          niso(50)=10         
          prob(50,1)=0.97     
          prob(50,2)=0.66     
          prob(50,3)=0.34     
          prob(50,4)=14.54    
          prob(50,5)=7.68     
          prob(50,6)=24.22    
          prob(50,7)=8.59     
          prob(50,8)=32.58    
          prob(50,9)=4.63     
          prob(50,10)=5.79    
          massiso(50,1)=111.904848
          massiso(50,2)=113.902779
          massiso(50,3)=114.903342
          massiso(50,4)=115.901741
          massiso(50,5)=116.902952
          massiso(50,6)=117.901603
          massiso(50,7)=118.903308
          massiso(50,8)=119.902194
          massiso(50,9)=121.903439
          massiso(50,10)=123.905274
                              
   ! 51 Sb (Antimony)
          niso(51)=2          
          prob(51,1)=57.21    
          prob(51,2)=42.79    
          massiso(51,1)=120.903816
          massiso(51,2)=122.904214
                              
   ! 52 Te (Tellurium)
          niso(52)=8          
          prob(52,1)=0.09     
          prob(52,2)=2.55     
          prob(52,3)=0.89     
          prob(52,4)=4.74     
          prob(52,5)=7.07     
          prob(52,6)=18.84    
          prob(52,7)=31.75    
          prob(52,8)=34.09    
          massiso(52,1)=119.90402
          massiso(52,2)=121.903044
          massiso(52,3)=122.904270
          massiso(52,4)=123.902817
          massiso(52,5)=124.904431
          massiso(52,6)=125.903312
          massiso(52,7)=127.904463
          massiso(52,8)=129.906224
  
   ! 53  I (Iodine)
          niso(53)=1          
          prob(53,1)=100.00   
          massiso(53,1)=126.904473
                              
   ! 78 Pt (Platinum)
          niso(78)=5          
          prob(78,1)=0.782    
          prob(78,2)=32.97    
          prob(78,3)=33.83    
          prob(78,4)=25.24    
          prob(78,5)=7.16     
          massiso(78,1)=191.961038
          massiso(78,2)=193.962680
          massiso(78,3)=194.964791
          massiso(78,4)=195.964951
          massiso(78,5)=197.967893
                              
   ! 80 Hg (Mercury)
          niso(80)=7          
          prob(80,1)=0.15     
          prob(80,2)=10.0     
          prob(80,3)=16.87    
          prob(80,4)=23.10    
          prob(80,5)=13.19    
          prob(80,6)=29.86    
          prob(80,7)=6.87     
          massiso(80,1)=195.965833
          massiso(80,2)=197.966769
          massiso(80,3)=198.968279
          massiso(80,4)=199.968326
          massiso(80,5)=200.970302
          massiso(80,6)=201.970643
          massiso(80,7)=203.973494
                              
   ! 81 Tl (Thallium)
          niso(81)=2          
          prob(81,1)=29.52    
          prob(81,2)=70.48    
          massiso(81,1)=202.972344
          massiso(81,2)=204.974427
  
   ! 82 Pb (Lead)
          niso(82)=4          
          prob(82,1)=1.3      
          prob(82,2)=24.1     
          prob(82,3)=22.1     
          prob(82,4)=52.4     
          massiso(82,1)=203.973043
          massiso(82,2)=205.974465
          massiso(82,3)=206.975897
          massiso(82,4)=207.976652
                              
   ! 83 Bi (Bismuth)
          niso(83)=1          
          prob(83,1)=100.     
          massiso(83,1)=208.980398
           
                                                                                                                                                                                                    
  ! postprocessing starts here 
  
  prob = prob * 0.01
  
  ! JK mass currently loops to element 83 (Bismuth)
  do i=1,83
     xm=0
     do j=1,niso(i)
        xm=xm+prob(i,j)
     enddo
     if(niso(i).gt.0.and.abs(xm-1.).gt.0.01) &
     stop 'internal isotope error 1'
  enddo
  
  do i=1,nat
     if(it(i).gt.100)then
        niso   (it(i))  =1
        prob   (it(i),1)=1.
        massiso(it(i),1)=it(i)-100.
     endif
  enddo
  
  do i=1,nat
     if(niso(it(i)).eq.0) stop 'internal isotope error 2'
  enddo

  !Calculate NO isotope pattern
  if (no_isotopes.eq..true.)then 
    nmass=0
    !do n=1,nrnd
    xmass=0
    do i=1,nat
      iti=it(i)
      p1=0.0
      do iso=1,niso(iti)
         if (prob(iti,iso).gt.p1)then
           x=massiso(iti,iso)
         endif
         p1=prob(iti,iso)
      enddo
      xmass = xmass + x
    enddo
    imass=nint(xmass)                  ! here it gets to integers - changed to nearest int nint
    nmass(imass)=nmass(imass)+1
    !enddo

  ! Calculate isotope pattern
  else                      
    nmass=0
    do n=1,nrnd
      xmass=0
      do i=1,nat
         iti=it(i)
         r=rnd(i,n)
         p1=0.0
         p2=prob(iti,1)
         do iso=1,niso(iti)
            if(r.ge.p1.and.r.le.p2)then
               x=massiso(iti,iso)
               exit
            endif
            p1=p2
            p2=p2+prob(iti,iso+1)
         enddo
         xmass=xmass+x
      enddo
      imass=nint(xmass)                  ! here it gets to integers - changed to nearest int nint
      nmass(imass)=nmass(imass)+1
    enddo
  endif
  
  isum=sum(nmass)
  k=0
  do i=1,2000
     if(nmass(i).gt.0)then
        k=k+1
        mass(k)=i
        mint(k)=float(nmass(i))/float(isum)
     endif
  enddo
  
  nsig=k
  end
  
  
  ! TK sub function as explained below
        REAL FUNCTION snorm()
  
  !C**********************************************************************C
  !C                                                                      C
  !C                                                                      C
  !C     (STANDARD-)  N O R M A L  DISTRIBUTION                           C
  !C                                                                      C
  !C                                                                      C
  !C**********************************************************************C
  !C**********************************************************************C
  !C                                                                      C
  !C     FOR DETAILS SEE:                                                 C
  !C                                                                      C
  !C               AHRENS, J.H. AND DIETER, U.                            C
  !C               EXTENSIONS OF FORSYTHE'S METHOD FOR RANDOM             C
  !C               SAMPLING FROM THE NORMAL DISTRIBUTION.                 C
  !C               MATH. COMPUT., 27,124 (OCT. 1973), 927 - 937.          C
  !C                                                                      C
  !C     ALL STATEMENT NUMBERS CORRESPOND TO THE STEPS OF ALGORITHM 'FL'  C
  !C     (M=5) IN THE ABOVE PAPER     (SLIGHTLY MODIFIED IMPLEMENTATION)  C
  !C                                                                      C
  !C     Modified by Barry W. Brown, Feb 3, 1988 to use RANF instead of   C
  !C     SUNIF.  The argument IR thus goes away.                          C
  !C                                                                      C
  !C**********************************************************************C
  !C
        DIMENSION a(32),d(31),t(31),h(31)
  !C
  !C     THE DEFINITIONS OF THE CONSTANTS A(K), D(K), T(K) AND
  !C     H(K) ARE ACCORDING TO THE ABOVEMENTIONED ARTICLE
  !C
  
        DATA a/0.0,.3917609E-1,.7841241E-1,.1177699,.1573107,.1970991, &
             .2372021,.2776904,.3186394,.3601299,.4022501,.4450965, &
             .4887764,.5334097,.5791322,.6260990,.6744898,.7245144, &
             .7764218,.8305109,.8871466,.9467818,1.009990,1.077516, &
             1.150349,1.229859,1.318011,1.417797,1.534121,1.675940, &
             1.862732,2.153875/
        DATA d/5*0.0,.2636843,.2425085,.2255674,.2116342,.1999243, &
             .1899108,.1812252,.1736014,.1668419,.1607967,.1553497, &
             .1504094,.1459026,.1417700,.1379632,.1344418,.1311722, &
             .1281260,.1252791,.1226109,.1201036,.1177417,.1155119, &
             .1134023,.1114027,.1095039/
        DATA t/.7673828E-3,.2306870E-2,.3860618E-2,.5438454E-2, &
             .7050699E-2,.8708396E-2,.1042357E-1,.1220953E-1,.1408125E-1, &
             .1605579E-1,.1815290E-1,.2039573E-1,.2281177E-1,.2543407E-1, &
             .2830296E-1,.3146822E-1,.3499233E-1,.3895483E-1,.4345878E-1, &
             .4864035E-1,.5468334E-1,.6184222E-1,.7047983E-1,.8113195E-1, &
             .9462444E-1,.1123001,.1364980,.1716886,.2276241,.3304980, &
             .5847031/
        DATA h/.3920617E-1,.3932705E-1,.3950999E-1,.3975703E-1, &
             .4007093E-1,.4045533E-1,.4091481E-1,.4145507E-1,.4208311E-1, &
             .4280748E-1,.4363863E-1,.4458932E-1,.4567523E-1,.4691571E-1, &
             .4833487E-1,.4996298E-1,.5183859E-1,.5401138E-1,.5654656E-1, &
             .5953130E-1,.6308489E-1,.6737503E-1,.7264544E-1,.7926471E-1, &
             .8781922E-1,.9930398E-1,.1155599,.1404344,.1836142,.2790016, &
             .7010474/
  
  !C
     10 call random_number(u)
        s = 0.0
        IF (u.GT.0.5) s = 1.0
        u = u + u - s
     20 u = 32.0*u
        i = int(u)
        IF (i.EQ.32) i = 31
        IF (i.EQ.0) GO TO 100
  !C
  !C                                START CENTER
  !C
     30 ustar = u - float(i)
        aa = a(i)
     40 IF (ustar.LE.t(i)) GO TO 60
        w = (ustar-t(i))*h(i)
  !C
  !C                                EXIT   (BOTH CASES)
  !C
     50 y = aa + w
        snorm = y
        IF (s.EQ.1.0) snorm = -y
        RETURN
  !C
  !C                                CENTER CONTINUED
  !C
     60 call random_number(u)
        w = u* (a(i+1)-aa)
        tt = (0.5*w+aa)*w
        GO TO 80
  
     70 tt = u
        call random_number(ustar)
     80 IF (ustar.GT.tt) GO TO 50
     90 call random_number(u)
        IF (ustar.GE.u) GO TO 70
        call random_number(ustar)
        GO TO 40
  !C
  !C                                START TAIL
  !C
    100 i = 6
        aa = a(32)
        GO TO 120
  
    110 aa = aa + d(i)
        i = i + 1
    120 u = u + u
        IF (u.LT.1.0) GO TO 110
    130 u = u - 1.0
    140 w = u*d(i)
        tt = (0.5*w+aa)*w
        GO TO 160
  
    150 tt = u
    160 call random_number(ustar)
        IF (ustar.GT.tt) GO TO 50
    170 call random_number(u)
        IF (ustar.GE.u) GO TO 150
        call random_number(u)
        GO TO 140
  
        END
  
  !     *****************************************************************
  ! TK  This routine handles
  ! TK  Called by: Main program
  ! TK  USES:      READAA
  
        SUBROUTINE READL(A1,X,N)
        IMPLICIT REAL*8 (A-H,O-Z)
        CHARACTER*(*) A1
        DIMENSION X(*)
        I=0
        IS=1
    10  I=I+1
        X(I)=READAA(A1,IS,IB,IE)
        IF(IB.GT.0 .AND. IE.GT.0) THEN
                                  IS=IE
                                  GOTO 10
        ENDIF
        N=I-1
        RETURN
        END
  
  !     *****************************************************************
  ! TK  This function handles
  ! TK  Called by: Main program
  ! TK  USES:      no other sub
  
        FUNCTION READAA(A,ISTART,IEND,IEND2)
        IMPLICIT REAL*8 (A-H,O-Z)
        REAL*8 READAA
        CHARACTER*(*) A
        NINE=ICHAR('9')
        IZERO=ICHAR('0')
        MINUS=ICHAR('-')
        IDOT=ICHAR('.')
        ND=ICHAR('D')
        NE=ICHAR('E')
        IBL=ICHAR(' ')
        IEND=0
        IEND2=0
        IDIG=0
        C1=0
        C2=0
        ONE=1.D0
        X = 1.D0
        NL=LEN(A)
        DO 10 J=ISTART,NL-1
           N=ICHAR(A(J:J))
           M=ICHAR(A(J+1:J+1))
           IF(N.LE.NINE.AND.N.GE.IZERO .OR.N.EQ.IDOT)GOTO 20
           IF(N.EQ.MINUS.AND.(M.LE.NINE.AND.M.GE.IZERO .OR. M.EQ.IDOT)) GOTO 20
     10 CONTINUE
        READAA=0.D0
        RETURN
     20 CONTINUE
        IEND=J
        DO 30 I=J,NL
           N=ICHAR(A(I:I))
           IF(N.LE.NINE.AND.N.GE.IZERO) THEN
              IDIG=IDIG+1
              IF (IDIG.GT.10) GOTO 60
              C1=C1*10+N-IZERO
           ELSEIF(N.EQ.MINUS.AND.I.EQ.J) THEN
              ONE=-1.D0
           ELSEIF(N.EQ.IDOT) THEN
              GOTO 40
           ELSE
              GOTO 60
           ENDIF
     30 CONTINUE
     40 CONTINUE
        IDIG=0
        DO 50 II=I+1,NL
           N=ICHAR(A(II:II))
           IF(N.LE.NINE.AND.N.GE.IZERO) THEN
              IDIG=IDIG+1
              IF (IDIG.GT.10) GOTO 60
              C2=C2*10+N-IZERO
              X = X /10
           ELSEIF(N.EQ.MINUS.AND.II.EQ.I) THEN
              X=-X
           ELSE
              GOTO 60
           ENDIF
     50 CONTINUE
  
  !C
  !C PUT THE PIECES TOGETHER
  !C
  
     60 CONTINUE
        READAA= ONE * ( C1 + C2 * X)
        DO 55 J=IEND,NL
           N=ICHAR(A(J:J))
           IEND2=J
           IF(N.EQ.IBL)RETURN
     55 IF(N.EQ.ND .OR. N.EQ.NE)GOTO 57
        RETURN
  
     57 C1=0.0D0
        ONE=1.0D0
        DO 31 I=J+1,NL
           N=ICHAR(A(I:I))
           IEND2=I
           IF(N.EQ.IBL)GOTO 70
           IF(N.LE.NINE.AND.N.GE.IZERO) C1=C1*10.0D0+N-IZERO
           IF(N.EQ.MINUS)ONE=-1.0D0
     31 CONTINUE
     61 CONTINUE
     70 READAA=READAA*10**(ONE*C1)
        RETURN
        END
  
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
  !        if(i.eq.iexp(1,i)) then
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
  
  subroutine match(tmass,iexp,score,tmax)
  implicit none
  real*8 tmass(10000)
  real*8 iexp(2,10000)
  real*8 score
  real*8 tmax
  integer i,j,k
  !     weighted vectors
  real*8 w(2,10000)
  !     needed storage variables
  real*8 sum1,sum2,sum3,sum4,dot,fr
  real*8 norm(2)
  !     peak pair number pp
  integer pp
  !     peak pair matrix 
  real*8 p(2,10000)
  !     number of calculated peaks
  integer numcalc
  
  numcalc=0
  do i=1,10000
    if(tmass(i).ne.0.)numcalc=numcalc+1
  enddo
  
  !     weighted spectral vectors, m**3, int**0.5 scaling ! Not sure why 0.5
  !     weighted spectral vectors, m**3, int**0.6 scaling
  w=0.
    do i=1,10000
    do j=1,10000
       if(j.eq.iexp(1,i)) then
       w(2,j)=j**2*(100.*tmass(j)/tmax)**0.6
       w(1,j)=iexp(1,i)**2*iexp(2,i)**0.6
       endif
    enddo
    enddo
  
  norm=0.
  do i=1,10000
     norm(1)=norm(1)+w(1,i)**2
     norm(2)=norm(2)+w(2,i)**2
  enddo
  
  norm=sqrt(norm)
  
  do i=1,10000
     w(1,i)=w(1,i)/norm(1)
     w(2,i)=w(2,i)/norm(2)
  enddo
  
  dot=0.
  sum1=0.
  sum2=0.
  sum3=0.
  do i=1,10000
     sum1=sum1+w(1,i)*w(2,i)
     sum2=sum2+(w(1,i))**2
     sum3=sum3+(w(2,i))**2
  enddo
  
  dot=sum1**2/(sum2*sum3)
  
  ! masses and intensity scalement for the second term
  ! m**0
  ! int**1
  
  w=0.
  do i=1,10000
     do j=1,10000
        if(j.eq.iexp(1,i)) then
           w(2,j)=100.*tmass(j)/tmax
           w(1,j)=iexp(2,i)
        endif
     enddo
  enddo
  
  !calculate the norm
  
  norm=0.
  do i=1,10000
     norm(1)=norm(1)+w(1,i)**2
     norm(2)=norm(2)+w(2,i)**2
  enddo
  
  norm=sqrt(norm)
  
  do i=1,10000
     w(1,i)=w(1,i)/norm(1)
     w(2,i)=w(2,i)/norm(2)
  enddo
  
  pp=0
  sum4=0.
  fr=0.
  p=0.
  do i=1,10000
     if((w(1,i)*w(2,i)).ne.0) then
        pp=pp+1
        p(1,pp)=w(1,i)
        p(2,pp)=w(2,i)
  !         print*,p(1,i),p(2,i)
     endif
  enddo
  
  !     ardous loop external file calcfr.f90
  call calcfr(pp,p,sum4)
  if (pp.eq.numcalc) sum4=sum4+1.0d0
  
  fr=sum4/dble(pp)
  score=(numcalc*dot+pp*fr)/(numcalc+pp)
  return
  end 
    
  
  subroutine calcfr(pp,pair,sum4)
  !     computes the second expression in the mass spec match score 
  implicit none
  integer pp
  real*8 pair(2,pp)
  real*8 sum4
  integer i
  
  sum4=0.d0
  
  do i=1,pp
     if     (abs((pair(1,i)*pair(2,i-1)) / (pair(1,i-1)*pair(2,i)))  .lt.1.0d0)then
  
        sum4=sum4+(pair(1,i)*pair(2,i-1))/(pair(1,i-1)*pair(2,i))
     
     elseif (abs((pair(1,i)*pair(2,i-1)) / (pair(1,i-1)*pair(2,i))) .gt.1.0d0)then
  
       sum4=sum4+1/((pair(1,i)*pair(2,i-1))/(pair(1,i-1)*pair(2,i))) 
     
     elseif (abs((pair(1,i)*pair(2,i-1)) / (pair(1,i-1)*pair(2,i))) .eq.1.0d0)then
       sum4=sum4+1.0d0
     endif
  enddo
  
  return 
  end
  
               
