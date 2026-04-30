module monit

  implicit none

  private
  public :: monit__ini, monit__output
  
  integer,parameter :: NVAR=5
  integer,parameter :: IUNIT_DATA = 35
  integer,save :: IMAX
  integer,save :: JMAX

contains

  subroutine monit__ini( filename, intkt_monit, ntmax, imax_in, jmax_in, ylat ) !IN
    !
    integer,intent(in) :: intkt_monit     !! Interval of monitor output (hour)
    integer,intent(in) :: ntmax
    character(len=*),intent(in) :: filename   !! Output file name
    integer,intent(in) :: imax_in            !! imax
    integer,intent(in) :: jmax_in            !! jmax
    real(8),intent(in) :: ylat(jmax_in)      !! ylat(jmax)
    !
    character(len=30) :: fdata
    character(len=30) :: fctl
    character(len=60) :: vars(NVAR)
    integer :: iunit, n
    !
    IMAX = imax_in
    JMAX = jmax_in
    !
    fdata = trim(filename)//'.dr'
    fctl  = trim(filename)//'.ctl'
    !
    vars(1) = 'h      0  0   height              (m)'
    vars(2) = 'u      0  0   zonal wind          (m/s)'
    vars(3) = 'v      0  0   meridional wind     (m/s)'
    vars(4) = 'vor    0  0   vorticity           (1/s)'
    vars(5) = 'div    0  0   divergence          (1/s)'
    !
    iunit=31
    open(iunit,file=fctl,form='formatted',access='sequential')
!    write(iunit,'(a)') 'OPTIONS big_endian'
    write(iunit,'(a)') 'DSET ^'//fdata
    write(iunit,'(a)') 'TITLE '//filename
    write(iunit,'(a)') 'UNDEF -9.99E33'
    write(iunit,'(a,i6,a,f14.10)') 'XDEF ',IMAX,'  LINEAR  0.0 ',360.0d0/IMAX
    write(iunit,'(a,i6,a,f14.10,a,f14.10)') 'YDEF ',JMAX,'  LINEAR  ',ylat(JMAX),' ',ylat(JMAX-1)-ylat(JMAX)
    write(iunit,'(a,i4,a)') 'ZDEF 1     LEVELS  1000.0'
    write(iunit,'(a,i6,a,i4,a)') 'TDEF ',ntmax,'  LINEAR  00Z01JAN2000 ',intkt_monit,'HR'
    write(iunit,'(a,i4)') 'VARS  ',NVAR
    do n=1,NVAR
       write(iunit,'(a)') trim(vars(n))
    end do
    write(iunit,'(a)') 'ENDVARS'
    close(iunit)
    !
    open( IUNIT_DATA, file=fdata, form='unformatted', access='direct', &
     &    recl=4*imax*jmax )
    !
  end subroutine monit__ini


  subroutine monit__output( ivar, it, data ) !IN
    !
    !! Write 4 byte direct file
    !
    integer,intent(in) :: ivar
    integer,intent(in) :: it
    real(8),intent(in) :: data(IMAX,JMAX)   !! Data for output
    !
    real(4) :: data_r4(IMAX,JMAX)           !! 4 byte Data
    integer :: irec,j
    !
 !$OMP PARALLEL default(SHARED), private(j)
   !$OMP DO schedule(STATIC)
    do j=1,JMAX
       data_r4(1:IMAX,j) = data(1:IMAX,JMAX+1-j)
    end do
   !$OMP END DO
 !$OMP END PARALLEL
    !
    irec = nvar*(it-1) + ivar
    write(IUNIT_DATA,rec=irec) data_r4
    !
  end subroutine monit__output

end module monit
