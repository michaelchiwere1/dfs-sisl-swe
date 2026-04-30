module grads

  implicit none

  private
  public :: grads__outxy

contains

  subroutine grads__outxy( filename, n, imax, jmax, ylat, data ) !IN
    !
    !! Output xy-decomposed data to grads file
    !
    character(*),intent(in) :: filename   !! Output file name
    integer,intent(in) :: n               !! Record number
    integer,intent(in) :: imax            !! imax
    integer,intent(in) :: jmax            !! jmax
    real(8),intent(in) :: ylat(jmax)      !! ylat(jmax)
    real(8),intent(in) :: data(IMAX,JMAX) !! Local data for output
    !
    real(4) :: data_r4(IMAX,JMAX)     !! real(4) of data
    character(len=30) :: fdata
    character(len=30) :: fctl
    integer :: j
    !
    !
!$OMP PARALLEL default(SHARED), private(j)
!$OMP DO schedule(STATIC)
    do j=1,JMAX
       data_r4(:,j) = data(:,JMAX+1-j)
    end do
!$OMP END DO
!$OMP END PARALLEL
    !
    fdata = trim(filename)//'.dr'
    fctl  = trim(filename)//'.ctl.tmp'
    !
    !! Write 4-byte direct file
    call write_data( fdata, IMAX, JMAX, n, data_r4 ) !IN
      !
    if ( n == 1 ) then
      !! Write grads control file
      call write_ctlfile( fctl, fdata, filename, imax, jmax, n, &!IN
       &                  dlon=360.0d0/IMAX, ylat=YLAT(jmax:1:-1)     )  !IN(Opt.)
    end if
      !
!      write(6,*) 'grads__outxy: output( ',trim(fctl),', ',trim(fdata),' )'
    !
  end subroutine grads__outxy


  subroutine write_data( fdata, imax, jmax, n, data_r4 ) !IN
    !
    !! Write grads direct file
    !
    character(*),intent(in) :: fdata         !! Name of data file
    integer,intent(in) :: imax,jmax          !! Array size of data(:,:,:)
    integer,intent(in) :: n                  !! record number
    real(4),intent(in) :: data_r4(imax,jmax) !! Data for output
    !
    integer :: iunit
    !
    iunit=30
    open( iunit, file=fdata, form='unformatted', access='direct', &
     &    recl=4*imax*jmax )
    write(iunit,rec=n) data_r4
    close(iunit)
    !
  end subroutine write_data


  subroutine write_ctlfile( fctl, fdata, filename, imax, jmax, n, &!IN
   &                        dlon, ylat                          )  !IN(Opt.)
    !
    !! Write grads control file
    !
    character(len=*),intent(in) :: fctl      !! Name of grads control file
    character(len=*),intent(in) :: fdata     !! Name of grads data file
    character(len=*),intent(in) :: filename  !! File name
    integer,intent(in) :: imax, jmax       !! Array size of data(:,:,:)
    integer,intent(in) :: n                !! record number
    real(8),intent(in),optional :: dlon
    real(8),intent(in),optional :: ylat(jmax)
    !
    integer :: iunit
    !
    iunit=31
    open(iunit,file=fctl,form='formatted',access='sequential')
!    write(iunit,'(a)') 'OPTIONS big_endian'
    write(iunit,'(a)') 'DSET ^'//fdata
    write(iunit,'(a)') 'TITLE '//filename
    write(iunit,'(a)') 'UNDEF -9.99E33'
    if ( present(dlon) ) then
      write(iunit,'(a,i4,a,f8.4)') 'XDEF ',imax,'  LINEAR  0.0 ',dlon
    else
      write(iunit,'(a,i4,a)') 'XDEF ',imax,'  LINEAR   1  1'
    end if
    if ( present(ylat) ) then
      write(iunit,'(a,i4,a)') 'YDEF ',jmax,'  LEVELS'
      write(iunit,'( (5x,5(f8.3,1x)) )') ylat
    else
      write(iunit,'(a,i4,a)') 'YDEF ',jmax,'  LINEAR   1  1'
    end if
    write(iunit,'(a,i4,a)') 'ZDEF 1      LINEAR   1  1'
    write(iunit,'(a,i4,a)') 'TDEF 1000  LINEAR  JAN2000 1HR'
    write(iunit,'(a)') 'VARS    1'
    write(iunit,'(a,i4,a)') 'data  0  0  '//filename
    write(iunit,'(a)') 'ENDVARS'
    close(iunit)
    !
  end subroutine write_ctlfile

end module grads

