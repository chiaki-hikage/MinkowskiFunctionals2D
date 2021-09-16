program mfcalc

!  input
!     inft: data file including temperature and its 1st and 2nd derivatives 
!     infm: mask file
!     nbin: binning number
!     numin: minimum threshold nu 
!     numax: maximum threshold nu 
!  output
!     outf: output file of Minkowski Functionals

  use mf_tools
  implicit none

  integer :: nbin,npix
  real, pointer :: m(:)
  double precision, pointer :: f(:,:),mf(:,:)
  double precision :: sig(0:1),sk(0:2),numin,numax,dnu,fsky
  character(100) :: inft,infm,outf
  character(10) :: cnbin,cnumin,cnumax

  call getarg(1,inft)
  call getarg(2,infm)
  call getarg(3,outf)
  call getarg(4,cnbin)
  read(cnbin,*) nbin
  call getarg(5,cnumin)
  read(cnumin,*) numin
  call getarg(6,cnumax)
  read(cnumax,*) numax
  allocate(mf(nbin,0:2))
  dnu=(numax-numin)/float(nbin)

  call read_synfast(inft,f,npix)
  call convert_covderiv(f,npix)
  call readmask(infm,m,npix,fsky)
  call norm(f,m,sig(0),npix)
  call calcmoment(f,m,npix,sig,sk)
  call calcmf(f,m,npix,mf,numin,dnu,nbin)
  deallocate(f,m)
  call correct_binning(mf,sig,numin,dnu,nbin)
  call out_mf(outf,mf,sig,sk,numin,dnu,nbin)
  deallocate(mf)

end program mfcalc
