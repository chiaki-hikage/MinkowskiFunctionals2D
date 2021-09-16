module mf_tools

  double precision, parameter :: pi=3.141592653589793

contains

!***************************************************

  subroutine calcmf(f,m,npix,mf,numin,dnu,nbin)

!***************************************************

    implicit none

    integer :: n,umpix,ib,nb
    integer, intent(in) :: npix,nbin
    real, pointer :: m(:)
    double precision, pointer :: f(:,:)
    double precision, intent(in) :: numin,dnu
    double precision, intent(out) :: mf(nbin,0:2)
    double precision :: u0,u1,u2,u11,u12,u22,g2,nu,b,numax

    numax=numin+dnu*nbin
    mf=0.d0
    umpix=0
    
    do n=1,npix
       if (m(n)==0.) cycle
       umpix=umpix+1
       
       u0=f(n,1)
       u1=f(n,2)
       u2=f(n,3)
       u11=f(n,4)
       u12=f(n,5)
       u22=f(n,6)

       g2=u1*u1+u2*u2

       if (u0<numin) then
          cycle
       else if (u0>numax) then
          do nb=1,nbin
             mf(nb,0)=mf(nb,0)+1
          enddo
          cycle
       else
          b=(u0-numin)/dnu
          ib=int(b)+1
          do nb=1,ib-1
             mf(nb,0)=mf(nb,0)+1.d0
          enddo
          if (b-ib>-0.5) mf(ib,0)=mf(ib,0)+1.d0
          mf(ib,1)=mf(ib,1)+sqrt(g2)/4.d0/dnu
          mf(ib,2)=mf(ib,2)+(2*u1*u2*u12-u1*u1*u22-u2*u2*u11)/(2*pi*g2)/dnu
       endif
    enddo

    mf=mf/umpix
    
  end subroutine calcmf

!***************************************************

  subroutine calcmoment(f,m,npix,sig,sk)

!***************************************************

    implicit none

    integer :: n,umpix
    integer, intent(in) :: npix
    real, pointer :: m(:)
    double precision, pointer :: f(:,:)
    double precision, intent(out) :: sk(0:2)
    double precision, intent(inout) :: sig(0:1)
    double precision :: u0,u1,u2,u11,u12,u22,g2,lap,tau

    sk=0.d0
    umpix=0
    tau=0.d0
    
    do n=1,npix
       if (m(n)==0.) cycle
       umpix=umpix+1
       
       u0=f(n,1)
       u1=f(n,2)
       u2=f(n,3)
       u11=f(n,4)
       u12=f(n,5)
       u22=f(n,6)

       lap=u11+u22
       g2=u1*u1+u2*u2
       tau=tau+g2/2.

       sk(0)=sk(0)+u0*u0*u0
       sk(1)=sk(1)+u0*u0*lap
       sk(2)=sk(2)+2*g2*lap

    enddo
    tau=sqrt(tau/umpix)
    sig(1)=tau*sqrt(2.)*sig(0)
    
    sk=sk/umpix
    sk(1)=sk(1)/(2*tau*tau)
    sk(2)=sk(2)/(4*tau**4)

    write(*,*) 'sigma',real(sig)
    write(*,*) 'skewness:',real(sk)

  end subroutine calcmoment

!***************************************************

  subroutine out_mf(outf,mf,sig,sk,numin,dnu,nbin)

!       1st line: sigma, skewness values
!        After 2nd line:
!          1st col: threshold value
!          2-4th cols: MFs for the input data with mask
!          5-7th cols: Gaussian formulae using sigma for the input data
!          8-10th cols: 1st order perturbative component from skewness
!
!***************************************************

    implicit none

    integer :: n
    integer, intent(in) :: nbin
    double precision, intent(in) :: mf(nbin,0:2),sig(0:1),sk(0:2),numin,dnu
    double precision :: mfg(0:2),dmfp1(0:2),nu,tau
    character(100) :: outf
    tau=sig(1)/sig(0)/sqrt(2.)

    write(*,*) 'Output file of Minkowski Functionals: ',outf(1:len_trim(outf))
    open(2,file=outf,status='unknown')
    write(2,'(5(1pe15.8,1x))') real(sig),real(sk)
    do n=1,nbin
       nu=numin+(n-0.5)*dnu
       call mf_gauss(2,nu,tau,mfg)
       call dmf_pb1(2,nu,tau,sk,sig(0),dmfp1)
       write(2,'(10(1pe15.8,1x))') nu,mf(n,:),mfg,dmfp1
    enddo
    close(2)

  end subroutine out_mf

!***************************************************
  
  subroutine convert_covderiv(f,npix)

!    inclination angle of each pixel in HEALPix format
!    pixel numbering scheme is "Ring" 

!***************************************************

    implicit none

    integer :: n,nside,i,j
    integer, intent(in) :: npix
    double precision, pointer :: f(:,:)
    double precision :: z, theta(npix)

    nside=nint(sqrt(npix/12.))
    n=0
    do i=1,nside-1
       z=1.-i**2/3./float(nside)**2
       do j=1,4*i
          n=n+1
          theta(n)=acos(z)
       enddo
    enddo

    do i=nside,3*nside
       z=4./3.-2.*i/3./float(nside)
       do j=1,4*nside
          n=n+1
          theta(n)=acos(z)
       enddo
    enddo

    do i=3*nside+1,4*nside-1
       z=-1.+(4*nside-i)**2/3./float(nside)**2
       do j=1,4*(4*nside-i)
          n=n+1
          theta(n)=acos(z)
       enddo
    enddo

! convert covariant derivative 

    do n=1,npix
       f(n,5)=f(n,5)-f(n,3)/tan(theta(n))
       f(n,6)=f(n,6)+f(n,2)/tan(theta(n))
    enddo

  end subroutine convert_covderiv

!*******************************************

  subroutine norm(f,m,sig0,npix)

! normalize data by its standard deviation

!*******************************************

    implicit none

    integer  n,i
    integer, intent(in) :: npix
    real, pointer :: m(:)
    double precision, pointer :: f(:,:)
    double precision ::  mean,sig0

    mean=0.d0
    sig0=0.d0
    i=0
    do n=1,npix
       if (m(n)>0.) then
          i=i+1
          mean=mean+f(n,1)
          sig0=sig0+f(n,1)*f(n,1)
       endif
    enddo
    mean=mean/i
    sig0=sqrt((sig0-i*mean*mean)/(i-1))
    do n=1,npix
       f(n,1)=f(n,1)-mean
    enddo
    f=f/sig0

  end subroutine norm

!***************************************************

  subroutine readmask(inf,m,npix,fsky)

!  read mask file in binary
!     mask(i)=0:  i-th pixel is masked
!     mask(i)=1:  i-th pixel is unmasked and then used for MFs calculation
!  pixel numbering scheme is RING scheme
!
!***************************************************
    
    implicit none
    
    integer :: status,unit,readwrite,blocksize,hdutype,ntable,n
    integer :: nfound,irow,i,range(2),irowmax,tfields,nb,nf,npix2
    integer*8 :: felem,nelems
    integer, parameter :: nbite=4
    integer, intent(in) :: npix
    real, pointer :: m(:)
    real ::  nulle,nobs
    double precision, intent(out) :: fsky
    character(100) :: inf,comment
    logical :: anynull

    allocate(m(npix))
    
    if (inf(1:6)=='nomask') then 
       m=1.
       return
    endif

    status=0
    call ftgiou(unit,status)
    readwrite=0
    call ftopen(unit,inf,readwrite,blocksize,status)
    ntable=2
    call ftmahd(unit,ntable,hdutype,status)
    call ftgkyj(unit,'TFIELDS',tfields,comment,status)
    call ftgknj(unit,'NAXIS',1,2,range,nfound,status)

    nelems=range(1)/nbite/tfields
    irowmax=range(2)
    npix2=nelems*irowmax
    if (npix/=npix2) then
       write(*,*) 'pixel number does not agree'
       stop
    endif

    nulle=0.d0
    i=0
    do irow=1,irowmax
       call ftgcve(unit,1,irow,1,nelems,nulle,m(irow),anynull,status)
       if (m(irow)>0.) i=i+1
    enddo

    fsky=float(i)/float(irowmax)
    write(*,*) 'fraction of unmasked area:', fsky
    
    call ftclos(unit, status)
    call ftfiou(unit, status)

    if (status>0) call printerror(status)

    return

  end subroutine readmask

!*******************************************************

  double precision function mf_amp(d,k,tau)

!*******************************************************

!   amplitude of MFs for the Gaussian distribution in d-dimension
!   tau=(sig1/sig0)/sqrt(d)

    implicit none
    
    integer, intent(in) :: d,k
    double precision, intent(in) :: tau

    mf_amp=tau**k/(2.*pi)**((k+1)/2.)*om(d)/(om(d-k)*om(k))

    return
    
  end function mf_amp
  
!*******************************************************

  subroutine mf_gauss(d,nu,tau,mf)

!*******************************************************

!   MFs for the Gaussian distribution in d-dimension
!   tau=(sig1/sig0)/sqrt(d)

    implicit none
    
    integer, intent(in) :: d
    double precision, intent(in) :: nu,tau
    double precision, intent(out) :: mf(0:d)
    integer :: k

    do k=0,d
       mf(k)=mf_amp(d,k,tau)*exp(-nu**2/2.)*herm(k-1,nu)
    enddo

    return
    
  end subroutine mf_gauss

!*******************************************************

  subroutine dmf_pb1(d,nu,tau,sk,sig0,dmf)

!*******************************************************

!   MFs for the weakly non-Gaussian distribution in d-dimension

    implicit none
    
    integer, intent(in) :: d
    double precision, intent(in) :: nu,tau,sk(0:d),sig0
    double precision, intent(out) :: dmf(0:d)
    double precision :: dmf_sk
    integer :: k

    do k=0,d
       dmf_sk=sk(0)/6.*herm(k+2,nu)-k*sk(1)/4.*herm(k,nu)-k*(k-1)*sk(2)/4*herm(k-2,nu)
       dmf(k)=mf_amp(d,k,tau)*exp(-nu**2/2.)*dmf_sk
    enddo

    return
    
  end subroutine dmf_pb1

!*****************************************
  
  double precision function om(m)

!  om = pi**(k/2)/Gamma(k/2+1)

!*****************************************
    
    integer, intent(in) :: m
    
    if (m<0) then
       om=0
    else if (m==0) then
       om=1
    else if (m==1) then
       om=2
    else if (m==2) then
       om=pi
    else if (m==3) then
       om=4*pi/3
    else if (m==4) then
       om=pi*pi/2.
    else if (m==5) then
       om=8*pi*pi/15.
    else if (m==6) then
       om=pi**3/6.
    else if (m==7) then
       om=16*pi**3/105.
    else if (m==8) then
       om=pi**4/24.
    else if (m==9) then
       om=32*pi**4/945.
    else if (m==10) then
       om=pi**5/120.
    else
       write(*,*) 'omega with m>10 is not described' 
       stop
    endif
    
    return
    
  end function om

!*****************************************

  double precision function herm(m,x)

!   m-th order Hermite polynomials
!
!*****************************************

    implicit none
    
    integer, intent(in) ::  m
    double precision, intent(in) :: x

    if (m<-1) then
       herm=0.
    else if (m==-1) then
       herm=sqrt(pi/2.)*exp(x**2/2.)*erfcc(dble(x/sqrt(2.)))
    elseif (m==0) then
       herm=1.
    elseif (m==1) then
       herm=x
    elseif (m==2) then
       herm=x*x-1.
    elseif (m==3) then
       herm=x**3-3*x
    elseif (m==4) then
       herm=x**4-6*x*x+3
    elseif (m==5) then
       herm=x**5-10*x**3+15*x
    elseif (m==6) then
       herm=x**6-15*x**4+45*x*x-15
    elseif (m==7) then
       herm=x**7-21*x**5+105*x**3-105*x
    elseif (m==8) then
       herm=x**8-28*x**6+210*x**4-420*x*x+105
    elseif (m==9) then
       herm=x**9-36*x**7+378*x**5-1260*x**3+945*x
    elseif (m==10) then
       herm=x**10-45*x**8+630*x**6-3150*x**4+4725*x*x-945
    else
       write(*,*) '>10-th order hermite polynomials are not described'
       stop
    endif

    return

  end function herm

!*****************************************

  double precision FUNCTION erfcc(x)

! calculates the complimentary error function
! Taken from Numerical Recipes. 

!*****************************************

    DOUBLE PRECISION, intent(in) :: x
    DOUBLE PRECISION :: t,z
    z=abs(x)
    t=1./(1.+0.5*z)
    erfcc=t*exp(-z*z-1.26551223+t*(1.00002368+t*(.37409196+t* &
         (.09678418+t*(-.18628806+t*(.27886807+t*(-1.13520398+t* &
         (1.48851587+t*(-.82215223+t*.17087277)))))))))
    if (x<0.) erfcc=2.-erfcc
    return

  END FUNCTION erfcc

!***************************************************

  subroutine read_synfast(inf,f,npix)

! read Temperature + 1st & 2nd derivatives
! output of synfast routine of HEALPix
! pixel numbering scheme is RING scheme
!
!  1st col  T
!  2nd col  dT/dtheta
!  3rd col  dT/dphi/sin(theta)
!  4th col  d^2 T/dtheta^2, 
!  5th col  d^2T/dtheta/dphi/sin(theta)
!  6th col  d^2T/dphi^2/sin^2(theta)
!   theta: elevation angle, phi: azimuth angle

!***************************************************
    
    implicit none
    
    integer :: status,unit,readwrite,blocksize,hdutype,ntable,n
    integer :: nfound,irow,i,range(2),irowmax,tfields,nb,nf
    integer*8 :: felem,nelems
    integer, parameter :: nbite=4
    integer, intent(out) :: npix
    double precision, pointer :: f(:,:)
    double precision ::  nulld
    character(100) :: inf,comment
    logical :: anynull
    
    status=0
    call ftgiou(unit,status)
  
    readwrite=0
    call ftopen(unit,inf,readwrite,blocksize,status)

    ntable=2
    call ftmahd(unit,ntable,hdutype,status)

    call ftgkyj(unit,'TFIELDS',tfields,comment,status)
    call ftgknj(unit,'NAXIS',1,2,range,nfound,status)

    nelems=range(1)/nbite/tfields
    irowmax=range(2)
    npix=nelems*irowmax

    felem=1
    nulld=0.d0

    allocate(f(npix,tfields))
  
    do irow=1,irowmax
       nb=1+nelems*(irow-1)
       nf=nelems*irow
       do i=1,tfields
          call ftgcvd(unit,i,irow,felem,nelems,nulld,f(nb:nf,i),anynull,status)
       enddo
    enddo
    
    call ftclos(unit, status)
    call ftfiou(unit, status)

    if (status>0) call printerror(status)

    return

  end subroutine read_synfast

!********************************************************

  subroutine correct_binning(mf,sig,numin,dnu,nbin)

! correction of finite binning size (ref. Lim and Simon arXiv:1103.4300)

!********************************************************

    implicit none

    integer, intent(in) :: nbin
    double precision, intent(inout) :: mf(nbin,0:2)
    double precision, intent(in) :: numin,dnu,sig(0:1)
    double precision :: dmf(0:2),nu,tau
    integer :: n,k
    tau=sig(1)/sig(0)/sqrt(2.)

    do n=1,nbin
       nu=numin+(n-0.5)*dnu
       do k=1,2
          dmf(k)=dnu*dnu/24*herm(k+1,nu)*mf_amp(2,k,tau)*exp(-nu*nu/2)
          mf(n,k)=mf(n,k)-dmf(k)
       enddo
    enddo

  end subroutine correct_binning

! *************************************************************************
 
  subroutine printerror(status)

!  This subroutine prints out the descriptive text corresponding to the
!  error status value and prints out the contents of the internal
!  error message stack generated by FITSIO whenever an error occurs.

! *************************************************************************

    integer :: status
    character :: errtext*30,errmessage*80

!  Check if status is OK (no error); if so, simply return
    if (status<=0)return
      
!  The FTGERR subroutine returns a descriptive 30-character text string that
!  corresponds to the integer error status number.  A complete list of all
!  the error numbers can be found in the back of the FITSIO User's Guide.
    call ftgerr(status,errtext)
    write(*,*) 'FITSIO Error Status =',status,': ',errtext

!  FITSIO usually generates an internal stack of error messages whenever
!  an error occurs.  These messages provide much more information on the
!  cause of the problem than can be provided by the single integer error
!  status value.  The FTGMSG subroutine retrieves the oldest message from
!  the stack and shifts any remaining messages on the stack down one
!  position.  FTGMSG is called repeatedly until a blank message is
!  returned, which indicates that the stack is empty.  Each error message
!  may be up to 80 characters in length.  Another subroutine, called
!  FTCMSG, is available to simply clear the whole error message stack in
!  cases where one is not interested in the contents.
    call ftgmsg(errmessage)
    do while (errmessage .ne. ' ')
       print *,errmessage
       call ftgmsg(errmessage)
    end do

  end subroutine printerror

end module mf_tools
