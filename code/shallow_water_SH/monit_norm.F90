module monit_norm

  implicit none

  private
  public :: monit_norm__ini, monit_norm__output

  integer,save :: IMAX
  integer,save :: JMAX
!  real(8),save :: SUM_WEIGHT
  real(8),save :: NORM_L1_INIT
  real(8),save :: NORM_L2_INIT
  real(8),save :: NORM_LMAX_INIT
  real(8),save,allocatable :: WEIGHT(:)
  real(8),save,allocatable :: INITIAL_DATA(:,:)

  integer,parameter :: IUNIT_CTL  = 61
  integer,parameter :: IUNIT_DATA = 62

contains

  subroutine monit_norm__ini( filename, intkt_monit, ntmax, imax_in, jmax_in, weight_in, initial_data_in )
    character(len=*),intent(in) :: filename   !! Output file name
    integer,intent(in) :: intkt_monit
    integer,intent(in) :: ntmax
    integer,intent(in) :: imax_in
    integer,intent(in) :: jmax_in
    real(8),intent(in) :: weight_in(jmax_in)
    real(8),intent(in) :: initial_data_in(imax_in,jmax_in)

    real(8) :: norm_l1_j(jmax_in), norm_l2_j(jmax_in), norm_lmax_j(jmax_in)
    character(len=30) :: fdata
    character(len=30) :: fctl
    character(len=60) :: vars(3)
    real(8) :: ww
    integer :: i,j,n

    IMAX = imax_in
    JMAX = jmax_in
    allocate( WEIGHT(JMAX) )
    allocate( INITIAL_DATA(IMAX,JMAX) )

   !$OMP PARALLEL default(SHARED), private(j)
   !$OMP DO schedule(STATIC)
    do j=1,JMAX
       WEIGHT(j) = weight_in(j)
       INITIAL_DATA(:,j) = initial_data_in(:,j)
    end do
   !$OMP END DO
   !$OMP END PARALLEL

!    SUM_WEIGHT = sum( WEIGHT(:) )*IMAX

    norm_l1_j(:)   = 0.0d0
    norm_l2_j(:)   = 0.0d0
    norm_lmax_j(:) = -9.99d33
    
   !$OMP PARALLEL default(SHARED), private(i,j,ww)
   !$OMP DO schedule(STATIC)
    do j=1,JMAX
       do i=1,IMAX
          ww = abs( INITIAL_DATA(i,j) )
          norm_l1_j(j)   = norm_l1_j(j) + ww
          norm_l2_j(j)   = norm_l2_j(j) + ww*ww
          norm_lmax_j(j) = max( norm_lmax_j(j), ww )
       end do
       norm_l1_j(j) = norm_l1_j(j)*WEIGHT(j)
       norm_l2_j(j) = norm_l2_j(j)*WEIGHT(j)
    end do
   !$OMP END DO
   !$OMP END PARALLEL
     
    NORM_L1_INIT   = sum( norm_l1_j(:) )
    NORM_L2_INIT   = sqrt( sum( norm_l2_j(:) ) )
    NORM_LMAX_INIT = maxval( norm_lmax_j(:) )
    
    write(6,*) "NORM_L1_INIT  =",NORM_L1_INIT
    write(6,*) "NORM_L2_INIT  =",NORM_L2_INIT
    write(6,*) "NORM_LMAX_INIT=",NORM_LMAX_INIT

    fctl  = trim(filename)//'.ctl'
    fdata = trim(filename)//'.dr'

    vars(1) = 'norm_L1    0  0   L1 norm'
    vars(2) = 'norm_L2    0  0   L2 norm'
    vars(3) = 'norm_Lmax  0  0   Lmax norm'

    open(IUNIT_CTL,file=fctl,form='formatted',access='sequential')
    write(IUNIT_CTL,'(a)') 'DSET ^'//fdata
    write(IUNIT_CTL,'(a)') 'TITLE '//filename
    write(IUNIT_CTL,'(a)') 'UNDEF -9.99E33'
    write(IUNIT_CTL,'(a)') 'XDEF 1  LEVELS  0.0 '
    write(IUNIT_CTL,'(a)') 'YDEF 1  LEVELS  0.0'
    write(IUNIT_CTL,'(a)') 'ZDEF 1  LEVELS  1000.0'
    write(IUNIT_CTL,'(a,i6,a,i4,a)') 'TDEF ',ntmax,'  LINEAR  00Z01JAN2000 ',intkt_monit,'HR'
    write(IUNIT_CTL,'(a,i4)') 'VARS  ',3
    do n=1,3
       write(IUNIT_CTL,'(a)') trim(vars(n))
    end do
    write(IUNIT_CTL,'(a)') 'ENDVARS'
    close(IUNIT_CTL)
    !
    open( IUNIT_DATA, file=fdata, form='unformatted', access='direct', recl=4 )

  end subroutine monit_norm__ini


  subroutine monit_norm__output( it, data )
    integer,intent(in) :: it
    real(8),intent(in) :: data(IMAX,JMAX)

    real(8) :: norm_l1_j(JMAX), norm_l2_j(JMAX), norm_lmax_j(JMAX)
    real(4) :: norm_l1, norm_l2, norm_lmax
    real(8) :: ww
    integer :: i,j,irec

    norm_l1_j(:)   = 0.0d0
    norm_l2_j(:)   = 0.0d0
    norm_lmax_j(:) = -9.99d33

   !$OMP PARALLEL default(SHARED), private(i,j,ww)
   !$OMP DO schedule(STATIC)
    do j=1,JMAX
       do i=1,IMAX
          ww = abs( data(i,j) - INITIAL_DATA(i,j) )
          norm_l1_j(j)   = norm_l1_j(j) + ww
          norm_l2_j(j)   = norm_l2_j(j) + ww*ww
          norm_lmax_j(j) = max( norm_lmax_j(j), ww )
       end do
       norm_l1_j(j) = norm_l1_j(j)*weight(j)
       norm_l2_j(j) = norm_l2_j(j)*weight(j)
    end do
   !$OMP END DO
   !$OMP END PARALLEL

    norm_l1   = sum( norm_l1_j(:) )/NORM_L1_INIT
    norm_l2   = sqrt( sum( norm_l2_j(:) ) )/NORM_L2_INIT
    norm_lmax = maxval( norm_lmax_j(:) )/NORM_LMAX_INIT
    !
    irec = 3*(it-1) + 1
    write(IUNIT_DATA,rec=irec  ) norm_l1
    write(IUNIT_DATA,rec=irec+1) norm_l2
    write(IUNIT_DATA,rec=irec+2) norm_lmax

  end subroutine monit_norm__output

end module monit_norm
