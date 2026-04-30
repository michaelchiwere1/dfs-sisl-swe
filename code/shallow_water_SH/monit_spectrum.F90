module monit_spectrum

#ifdef _OPENMP
  use omp_lib, only : omp_get_max_threads, omp_get_thread_num
#endif

  implicit none

  private
  public :: monit_spectrum__ini, monit_spectrum__output

  integer,save :: NUMPE_OMP
  integer,save :: MNWAV
  integer,save :: NMAX
  real(8),save,allocatable :: CCC(:)

  integer,parameter :: IUNIT_CTL  = 71
  integer,parameter :: IUNIT_DATA = 72

contains

  subroutine monit_spectrum__ini( filename, intkt_monit, ntmax, mnwav_in, nmax_in, earth )
    character(len=*),intent(in) :: filename   !! Output file 
    integer,intent(in) :: intkt_monit
    integer,intent(in) :: ntmax
    integer,intent(in) :: mnwav_in
    integer,intent(in) :: nmax_in
    real(8),intent(in) :: earth

    character(len=30) :: fdata
    character(len=30) :: fctl
    character(len=60) :: vars(2)
    real(8) :: an
    integer :: j,n

#ifdef _OPENMP
    NUMPE_OMP = omp_get_max_threads()
#else
    NUMPE_OMP = 1
#endif

    MNWAV = mnwav_in
    NMAX = nmax_in

    allocate( CCC(NMAX) )
! !$OMP PARALLEL default(SHARED), private(n,an)
!   !$OMP DO schedule(STATIC)
    do n=1,NMAX
       an=n
       CCC(n) = 0.5d0*earth**2/( an*(an+1) )
    end do
!   !$OMP END DO
! !$OMP END PARALLEL

    fctl  = trim(filename)//'.ctl'
    fdata = trim(filename)//'.dr'

    vars(1) = 'KE_VOR    0  0   Kinetic Energy (Vorticity)'
    vars(2) = 'KE_DIV    0  0   Kinetic Energy (Divergence)'

    open(IUNIT_CTL,file=fctl,form='formatted',access='sequential')
    write(IUNIT_CTL,'(a)') 'DSET ^'//fdata
    write(IUNIT_CTL,'(a)') 'TITLE '//filename
    write(IUNIT_CTL,'(a)') 'UNDEF -9.99E33'
    write(IUNIT_CTL,'(a,i5,a)') 'XDEF ', NMAX, '  LINEAR 1.0 1.0'
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
    open( IUNIT_DATA, file=fdata, form='unformatted', access='direct', recl=4*NMAX )

  end subroutine monit_spectrum__ini


  subroutine monit_spectrum__output( it, ntotal, mwave, mnstart, qrot, qdiv )
    integer,intent(in) :: it
    integer,intent(in) :: ntotal(MNWAV)
    integer,intent(in) :: mwave(MNWAV)
    integer,intent(in) :: mnstart(0:NMAX+1)
    real(8),intent(in) :: qrot(2,MNWAV)
    real(8),intent(in) :: qdiv(2,MNWAV)

    real(4) :: ke_rot(NMAX)
    real(4) :: ke_div(NMAX)
    real(8) :: work_rot(0:NMAX,0:NUMPE_OMP-1)
    real(8) :: work_div(0:NMAX,0:NUMPE_OMP-1)
    integer,parameter :: nn=20
    real(8) :: work(nn)
    integer :: iwork(nn)
    real(8) :: ww
    integer :: np,mn,n,m,irec

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
       n = ntotal(mn)
       m = mwave(mn)
       if ( m == 0 ) then
          work_rot(n,np) = work_rot(n,np) + qrot(1,mn)**2
          work_div(n,np) = work_div(n,np) + qdiv(1,mn)**2
       else
          work_rot(n,np) = work_rot(n,np) + 2.0d0*( qrot(1,mn)**2 + qrot(2,mn)**2 )
          work_div(n,np) = work_div(n,np) + 2.0d0*( qdiv(1,mn)**2 + qdiv(2,mn)**2 )
       end if
    end do
   !$OMP END DO

   !$OMP DO schedule(STATIC)
    do n=1,NMAX
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

end module monit_spectrum
