module e_time

  !! Measure elapsed wall-clock time
  
  implicit none

  private
  public :: e_time__start, e_time__end, e_time__output

  integer,parameter :: NMAX    = 100  !! Maximum e_time number
  integer,parameter :: NOUTMAX = 100  !! Maximum e_time number (Output)
  integer,parameter :: N_OUT(2) = (/ 1, NOUTMAX /)      !! Elapsed time from n=1 to n=NOUTMAX is outputted

  integer,parameter :: NNESTMAX = 50 !! Maximum nesting depth of e_time__start/e_time__end blocks

  logical,save,allocatable :: FIRST_LABEL(:)
  character(len=20),save,allocatable :: LABEL(:) !! e_time label for each e_time number and thread number
  integer,save,allocatable :: NCALL(:)           !! Calling number for each ..
  logical,save,allocatable :: ALREADY_START(:)   !! e_time__start is already called or not for each ..

  integer,save,allocatable :: N_NEST(:)     !! Nested e_time numbers for each thread number
  integer,save :: IPTR_N_NEST  !! pointer for n_nest(:)

  integer(8),save,allocatable :: ICOUNT_SUM(:)   !! Summation of timer count for each ..

  integer(8),save :: ICOUNT_BEFORE
  integer(8),save :: ICRATE
  integer(8),save :: ICMAX

  logical,save :: L_FIRST = .true.
  logical,save :: L_OUTPUT(NMAX)

contains


  subroutine e_time__start(n,label_in) !IN
    integer,intent(in) :: n                  !! e_time number
    character(len=*),intent(in) :: label_in  !! e_time label
    integer :: n_before
    integer(8) :: icount,icount_diff
    
!    write(0,*) "e_time__start: label=",label_in

    if ( n < 1 .or. n > NMAX ) then
      write(6,*) "Error: e_time__start (e_time.F90): n is out of range."
      write(6,*) "       n, NMAX =",n,NMAX
      stop 999
    end if

    if ( L_FIRST ) then
      L_FIRST = .false.
      allocate( FIRST_LABEL  (NMAX) )
      FIRST_LABEL(:) = .true.
      allocate( LABEL        (NMAX) )
      LABEL(:) = ' '
      allocate( NCALL        (NMAX) )
      allocate( ALREADY_START(NMAX) )
      allocate( N_NEST       (NNESTMAX) )

      allocate( ICOUNT_SUM   (NMAX) )
      call system_clock(icount,ICRATE,ICMAX)
      ICOUNT_SUM(:) = 0
      ICOUNT_BEFORE = icount

      NCALL(:) = 0
      ALREADY_START(:) = .false.
      IPTR_N_NEST = 0

      L_OUTPUT(:) = .false.
      if ( N_OUT(1) /= -999 ) L_OUTPUT(max(N_OUT(1),1):min(N_OUT(2),NMAX)) = .true.
    end if

    if ( .not. L_OUTPUT(n) ) then
      return
    end if

    if ( ALREADY_START(n) ) then
      write(6,*)               "Error: e_time__start (e_time.F90): "
      write(6,'(A,I3,A,I3,A)') "       e_time__start(",n,",'"//trim(label_in)//"') is called twice " &
       &                               //"before calling e_time__end(",n,",'"//trim(label_in)//"')."
      stop 999
    end if
    ALREADY_START(n) = .true.

    if ( FIRST_LABEL(n) ) then
      LABEL(n) = label_in
      FIRST_LABEL(n) = .false.
    end if

    call system_clock(icount)
    if ( IPTR_N_NEST > 0 ) then
      n_before = N_NEST(IPTR_N_NEST)
      icount_diff = icount - ICOUNT_BEFORE
      if ( icount_diff < 0 ) then
        icount_diff = icount_diff + ICMAX
      end if
      ICOUNT_SUM(n_before) = ICOUNT_SUM(n_before) + icount_diff
    end if
    ICOUNT_BEFORE = icount

    NCALL(n) = NCALL(n) + 1
    IPTR_N_NEST = IPTR_N_NEST + 1
    if ( IPTR_N_NEST > NNESTMAX ) then
      write(6,*) "Error: e_time_start (e_time.F90): IPTR_N_NEST > NNESTMAX"
      write(6,*) "       IPTR_N_NEST = ",IPTR_N_NEST
      write(6,*) "       NNESTMAX = ",NNESTMAX
      stop 999
    end if
    N_NEST(IPTR_N_NEST) = n

    return
  end subroutine e_time__start


  subroutine e_time__end(n,label_in) !IN
    integer,intent(in) :: n                  !! e_time number
    character(len=*),intent(in) :: label_in  !! e_time label
    integer :: n_before
    integer(8) :: icount,icount_diff
    
!    write(0,*) "e_time__end: label=",label_in

    if ( n < 1 .or. n > NMAX ) then
      write(6,*) "Error: e_time_end (e_time.F90): n is out of range."
      write(6,*) "       n, NMAX =",n,NMAX
      stop 999
    end if

    if ( .not. L_OUTPUT(n) ) then
      return
    end if

    if ( .not. ALREADY_START(n) .or. IPTR_N_NEST<=0 ) then
      write(6,'(A)')           "Error: e_time__end (e_time.F90):"
      write(6,'(A,I3,A,I3,A)') "       e_time__start(",n,",'"//trim(label_in)//"') is not called " &
       &                       //"before calling e_time__end(",n,",'"//trim(label_in)//"')."
      stop 999
    end if

    n_before = N_NEST(IPTR_N_NEST)
    if ( n_before /= n ) then
      write(6,'(A)')           "Error: e_time__end (e_time.F90):"
      write(6,'(A,I3,A,I3,A)') "       e_time__end(",n,",'"//trim(label_in)//"') is called "                 &
       &                       //"before calling e_time__end(",n_before,",'"//trim(LABEL(n_before))//"')."
      stop 999
    end if
    
    ALREADY_START(n) = .false.

    call system_clock(icount)
    icount_diff = icount - ICOUNT_BEFORE
    if ( icount_diff < 0 ) then
      icount_diff = icount_diff + ICMAX
    end if
    ICOUNT_SUM(n) = ICOUNT_SUM(n) + icount_diff
    ICOUNT_BEFORE = icount

    IPTR_N_NEST = IPTR_N_NEST - 1
    
  end subroutine e_time__end


  subroutine e_time__output
    real(8) :: time(NMAX)
    real(8) :: total
    integer :: n_before, i

    if ( IPTR_N_NEST > 0 ) then
      n_before = N_NEST(IPTR_N_NEST)
      write(6,'(A,I3,A)') "Error: e_time__output (e_time.F90): e_time__output is called "      &
       &                  //"before calling e_time__end(",n_before,",'"//trim(LABEL(n_before))//"')."
      stop 999
    end if

    do i=1,NMAX
      time(i) = 1.0d0*ICOUNT_SUM(i)/ICRATE
    end do

    total = 0.0d0
    do i=1,NMAX
      total = total + time(i)
    end do

    write(6,'(A)') ' '
    write(6,'(A)') '  ======== ELAPSED TIME ========'
    write(6,'(7X, A20, A10, A12, A10)')                                  &
     & '  SUBROUTINE        ', '    CPU(S)', '     RATE(%)',  '    NUMBER'
    do i=1,NMAX
      if ( time(i) > 0.0d0 ) then
        write(6,'(I5, 2X, A20, F10.2, F12.2, I10)')               &
         & i, LABEL(i), time(i), time(i)/total*100.0d0, NCALL(i)
      end if
    end do
    write(6,'(3X, A, 1X, F11.2)')       &
     & '--- TOTAL TIME (S) ---', total

    return
  end subroutine e_time__output

end module e_time

