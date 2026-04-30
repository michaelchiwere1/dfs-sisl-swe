module monit_spectrum

#ifdef _OPENMP
  use omp_lib, only : omp_get_max_threads, omp_get_thread_num
#endif
  use prm_phconst, only : PI
  use legendre, only : legendre__ini, legendre__g2w
  implicit none

  private
  public :: monit_spectrum__ini, monit_spectrum__output

  integer,save :: JCN_DFS
  integer,save :: NUMPE_OMP
  integer,save :: NMAX_SH
  integer,save :: NNUM_SH
  integer,save :: NMAX
  integer,save :: MMAX
  integer,save :: NNUM
  integer,save :: MNUM
  integer,save :: MNWAV
  integer,save :: MNWAV_UV
  integer,save :: IMAX
  integer,save :: JMAX_SH
  integer,save,allocatable :: MNSTART(:)
  integer,save,allocatable :: MNSTART_UV(:)
  real(8),save,allocatable :: CCC(:)
  integer,save,allocatable :: NTOTAL(:)
  integer,save,allocatable :: MWAVE(:)

  real(8),save,allocatable :: GW(:)
  real(8),save,allocatable :: ALAT(:)
  real(8),save,allocatable :: SINLAT(:)
  real(8),save,allocatable :: COSLAT(:)
  real(8),save,allocatable :: COSLAT_INV(:)

  real(8),save,allocatable :: cosnphi(:,:)
  real(8),save,allocatable :: sinnphi(:,:)

  integer,parameter :: IUNIT_CTL  = 71
  integer,parameter :: IUNIT_DATA = 72

contains

  subroutine monit_spectrum__ini( filename, intkt_monit, ntmax, jcn_dfs_in, nmax_in, mmax_in, earth )
    character(len=*),intent(in) :: filename   !! Output file name
    integer,intent(in) :: intkt_monit
    integer,intent(in) :: ntmax
    integer,intent(in) :: jcn_dfs_in
    integer,intent(in) :: nmax_in
    integer,intent(in) :: mmax_in
    real(8),intent(in) :: earth

    character(len=30) :: fdata
    character(len=30) :: fctl
    character(len=60) :: vars(2)
    real(8) :: an,colat
    integer :: n,m,mn,j

#ifdef _OPENMP
    NUMPE_OMP = omp_get_max_threads()
#else
    NUMPE_OMP = 1
#endif

    JCN_DFS = jcn_dfs_in

    NMAX = nmax_in
    MMAX = mmax_in

    NNUM = NMAX + 1
    MNUM = MMAX + 1

    NNUM_SH = max( NNUM, MNUM )
!xx    NNUM_SH = NNUM*2
    NMAX_SH = NNUM_SH - 1
    IMAX    = MNUM*2
    JMAX_SH = NNUM_SH + mod(NNUM_SH,2)  !! JMAX should be even.
    allocate( MNSTART(0:MNUM) )
    allocate( MNSTART_UV(0:MNUM) )
    allocate( GW(JMAX_SH) )
    allocate( ALAT(JMAX_SH) )
    allocate( SINLAT(JMAX_SH) )
    allocate( COSLAT(JMAX_SH) )
    allocate( COSLAT_INV(JMAX_SH) )

    write(6,*) "monit_spectrum__ini: call legendre__ini for monitoring power energy spectrum."
    write(6,*) "  NNUM, NNUM_SH, JMAX_SH=",NNUM, NNUM_SH, JMAX_SH
    call legendre__ini                          &
     &( IMAX, JMAX_SH, NNUM_SH, MNUM,              &!IN
     &  MNWAV, MNWAV_UV, MNSTART, MNSTART_UV,   &!OUT
     &  gw, alat, sinlat, coslat, coslat_inv )   !OUT

    allocate( NTOTAL(MNWAV) )
    allocate( MWAVE(MNWAV) )
!$OMP PARALLEL default(SHARED), private(m,mn)
 !$OMP DO schedule(DYNAMIC)
    do m=0,MMAX
      do mn=MNSTART(m),MNSTART(m+1)-1
        NTOTAL(mn) = m + mn - MNSTART(m)
        MWAVE(mn) = m
      end do
    end do
 !$OMP END DO
!$OMP END PARALLEL

    allocate( CCC(NMAX_SH) )
 !$OMP PARALLEL default(SHARED), private(n,an)
   !$OMP DO schedule(STATIC)
    do n=1,NMAX_SH
       an=n
       CCC(n) = 0.5d0*earth**2/( an*(an+1) )
    end do
   !$OMP END DO
 !$OMP END PARALLEL

    allocate( cosnphi(NNUM,JMAX_SH) )
    allocate( sinnphi(NNUM,JMAX_SH) )
 !$OMP PARALLEL default(SHARED), private(j,n,colat)
   !$OMP DO schedule(STATIC)
    do j=1,JMAX_SH
       do n=1,NNUM
          colat = PI/2 - ALAT(j)
          cosnphi(n,j) = cos( (n-1)*colat )
          sinnphi(n,j) = sin( n*colat )
       end do
    end do
   !$OMP END DO
 !$OMP END PARALLEL

    fctl  = trim(filename)//'.ctl'
    fdata = trim(filename)//'.dr'

    vars(1) = 'KE_VOR    0  0   Kinetic Energy (Vorticity)'
    vars(2) = 'KE_DIV    0  0   Kinetic Energy (Divergence)'

    open(IUNIT_CTL,file=fctl,form='formatted',access='sequential')
    write(IUNIT_CTL,'(a)') 'DSET ^'//fdata
    write(IUNIT_CTL,'(a)') 'TITLE '//filename
    write(IUNIT_CTL,'(a)') 'UNDEF -9.99E33'
    write(IUNIT_CTL,'(a,i5,a)') 'XDEF ', NMAX_SH, '  LINEAR 1.0 1.0'
    write(IUNIT_CTL,'(a)') 'YDEF 1  LEVELS  0.0'
    write(IUNIT_CTL,'(a)') 'ZDEF 1  LEVELS  1000.0'
    write(IUNIT_CTL,'(a,i6,a,i4,a)') 'TDEF ',ntmax,'  LINEAR  00Z01JAN2000 ',intkt_monit,'HR'
    write(IUNIT_CTL,'(a,i4)') 'VARS  ',2
    do n=1,2
       write(IUNIT_CTL,'(a)') trim(vars(n))
    end do
    write(IUNIT_CTL,'(a)') 'ENDVARS'
    close(IUNIT_CTL)
    !
    open( IUNIT_DATA, file=fdata, form='unformatted', access='direct', recl=4*NMAX_SH )

  end subroutine monit_spectrum__ini


  subroutine monit_spectrum__output( it, qrot, qdiv )
    integer,intent(in) :: it
    complex(8),intent(inout) :: qrot(NNUM,0:MMAX)
    complex(8),intent(inout) :: qdiv(NNUM,0:MMAX)

    real(8) :: qrot_sh(2,MNWAV)
    real(8) :: qdiv_sh(2,MNWAV)

    real(4) :: ke_rot(NMAX_SH)
    real(4) :: ke_div(NMAX_SH)
    real(8) :: work_rot(0:NMAX_SH,0:NUMPE_OMP-1)
    real(8) :: work_div(0:NMAX_SH,0:NUMPE_OMP-1)
    real(8) :: data_sh(2,0:MMAX,JMAX_SH)
    integer :: np,mn,n,m,irec
    
!    integer,parameter :: nn=20
!    real(8) :: work(nn)
!    integer :: iwork(nn)
!    real(8) :: ww
!    integer :: m


!    qrot(:,:) = 0.0d0
!    qrot(1,0:3) = 1.0d0
!    qdiv(:,:) = 0.0d0

    call w2g_gauss  &
     &( qrot,       &!IN
     &  data_sh   )  !OUT
    call legendre__g2w &
     &( data_sh,       &!IN
     &  qrot_sh  )      !OUT

!    m=3
!    write(6,*) "m=",m
!    write(6,*) "qrot(:,m)=",qrot(:,m)
!    write(6,*) "data_sh(1,m,:)=",data_sh(1,m,:)
!    write(6,*) "qrot_sh(1,MNSTART(m):MNSTART(m+1)-1)=",qrot_sh(1,MNSTART(m):MNSTART(m+1)-1)
!    stop 333

    call w2g_gauss  &
     &( qdiv,       &!IN
     &  data_sh   )  !OUT
    call legendre__g2w &
     &( data_sh,       &!IN
     &  qdiv_sh  )      !OUT

 !$OMP PARALLEL default(SHARED), private(np,mn,m,n)

#ifdef _OPENMP
    np=omp_get_thread_num()
#else
    np=0
#endif

    work_rot(:,np) = 0.0d0
    work_div(:,np) = 0.0d0

   !$OMP DO schedule(STATIC)
    do mn=1,MNWAV
       n = NTOTAL(mn)
       m = MWAVE(mn)
       if ( m == 0 ) then
          work_rot(n,np) = work_rot(n,np) + qrot_sh(1,mn)**2
          work_div(n,np) = work_div(n,np) + qdiv_sh(1,mn)**2
       else
          work_rot(n,np) = work_rot(n,np) + 2.0d0*( qrot_sh(1,mn)**2 + qrot_sh(2,mn)**2 )
          work_div(n,np) = work_div(n,np) + 2.0d0*( qdiv_sh(1,mn)**2 + qdiv_sh(2,mn)**2 )
       end if
    end do
   !$OMP END DO

   !$OMP DO schedule(STATIC)
    do n=1,NMAX_SH
       ke_rot(n) = CCC(n)*sum( work_rot(n,0:NUMPE_OMP-1) )
       ke_div(n) = CCC(n)*sum( work_div(n,0:NUMPE_OMP-1) )
    end do
   !$OMP END DO
 !$OMP END PARALLEL
    !
    irec = 2*(it-1) + 1
    write(IUNIT_DATA,rec=irec  ) ke_rot
    write(IUNIT_DATA,rec=irec+1) ke_div


  end subroutine monit_spectrum__output


  subroutine w2g_gauss &
   &( qdata,       &!IN
   &  data_sh   )   !OUT

    complex(8),intent(in) :: qdata(NNUM,0:MMAX)
    real(8),intent(out) :: data_sh(2,0:MMAX,JMAX_SH)

    integer :: n,m,j

 !$OMP PARALLEL default(SHARED), private(n,m,j)

   !$OMP DO schedule(STATIC)
    do j=1,JMAX_SH
       data_sh(:,:,j) = 0.0d0
    end do
   !$OMP END DO

   !$OMP DO schedule(STATIC)
    do m=0,MMAX
       if ( JCN_DFS >= 1 ) then
          if ( m == 0 ) then
             do j=1,JMAX_SH
                do n=1,NNUM
                   data_sh(1,m,j) = data_sh(1,m,j) + real( qdata(n,m), kind=8 )*cosnphi(n,j)
                   data_sh(2,m,j) = data_sh(2,m,j) + dimag( qdata(n,m) )*cosnphi(n,j)
                end do
             end do
          else if ( m == 1 ) then
             do j=1,JMAX_SH
                do n=1,NNUM
                   data_sh(1,m,j) = data_sh(1,m,j) + real( qdata(n,m), kind=8 )*sinnphi(1,j)*cosnphi(n,j)
                   data_sh(2,m,j) = data_sh(2,m,j) + dimag( qdata(n,m) )*sinnphi(1,j)*cosnphi(n,j)
                end do
             end do
          else if ( mod(m,2) == 0 ) then  !! m=2,4,6,8,...
             do j=1,JMAX_SH
                do n=1,NNUM
                   data_sh(1,m,j) = data_sh(1,m,j) + real( qdata(n,m), kind=8 )*sinnphi(1,j)*sinnphi(n,j)
                   data_sh(2,m,j) = data_sh(2,m,j) + dimag( qdata(n,m) )*sinnphi(1,j)*sinnphi(n,j)
                end do
             end do
          else                            !! m=3,5,7,9,...
             if ( jcn_dfs == 1 ) then
                do j=1,JMAX_SH
                   do n=1,NNUM
                      data_sh(1,m,j) = data_sh(1,m,j) + real( qdata(n,m), kind=8 )*sinnphi(1,j)*cosnphi(n,j)
                      data_sh(2,m,j) = data_sh(2,m,j) + dimag( qdata(n,m) )*sinnphi(1,j)*cosnphi(n,j)
                   end do
                end do
             else
                do j=1,JMAX_SH
                   do n=1,NNUM
                      data_sh(1,m,j) = data_sh(1,m,j) + real( qdata(n,m), kind=8 )*sinnphi(1,j)**2*sinnphi(n,j)
                      data_sh(2,m,j) = data_sh(2,m,j) + dimag( qdata(n,m) )*sinnphi(1,j)**2*sinnphi(n,j)
                   end do
                end do
             end if
          end if
       else
          if ( m == 0 ) then
             do j=1,JMAX_SH
                do n=1,NNUM
                   data_sh(1,m,j) = data_sh(1,m,j) + real( qdata(n,m), kind=8 )*cosnphi(n,j)
                   data_sh(2,m,j) = data_sh(2,m,j) + dimag( qdata(n,m) )*cosnphi(n,j)
                end do
             end do
          else if ( mod(m,2) == 1 ) then  !! m=1,3,5,7,...
             do j=1,JMAX_SH
                do n=1,NNUM
                   data_sh(1,m,j) = data_sh(1,m,j) + real( qdata(n,m), kind=8 )*sinnphi(n,j)
                   data_sh(2,m,j) = data_sh(2,m,j) + dimag( qdata(n,m) )*sinnphi(n,j)
                end do
             end do
          else                            !! m=2,4,6,8,...
             do j=1,JMAX_SH
                do n=1,NNUM
                   data_sh(1,m,j) = data_sh(1,m,j) + real( qdata(n,m), kind=8 )*sinnphi(1,j)*sinnphi(n,j)
                   data_sh(2,m,j) = data_sh(2,m,j) + dimag( qdata(n,m) )*sinnphi(1,j)*sinnphi(n,j)
                end do
             end do
          end if
       end if
    end do
   !$OMP END DO
 !$OMP END PARALLEL

  end subroutine w2g_gauss

end module monit_spectrum
