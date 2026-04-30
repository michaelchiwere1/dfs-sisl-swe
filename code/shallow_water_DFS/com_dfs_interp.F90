module com_dfs_interp

  implicit none

  integer,save :: jcn_dfs
  integer,save :: jcn_grid
  
  integer,save :: NNUM
  integer,save :: NMAX
  integer,save :: MMAX
  integer,save :: IMAX
  integer,save :: JMAX
  integer,save :: NNUMHF1
  integer,save :: NNUMHF2

  integer,save :: N_TRUNC_M0_DFS
  integer,save :: N_TRUNC_M1_DFS
  integer,save :: N_TRUNC_M2_DFS
  integer,save :: N_TRUNC_M3_DFS

  integer,save :: N_L2TRUNC1_M0_DFS
  integer,save :: N_L2TRUNC2_M0_DFS

  integer,save :: N_L2TRUNC1_M1_DFS
  integer,save :: N_L2TRUNC2_M1_DFS

  integer,save :: N_L2TRUNC1_M2_DFS
  integer,save :: N_L2TRUNC2_M2_DFS

  integer,save :: N_L2TRUNC1_M3_DFS
  integer,save :: N_L2TRUNC2_M3_DFS
  
  integer,save,allocatable :: NNUM1_M_DFS(:)
  integer,save,allocatable :: NNUM2_M_DFS(:)
  integer,save,allocatable :: NNUM_M_DFS(:)

  private
  public :: com_dfs__ini, jcn_dfs, jcn_grid, MMAX, NNUM, IMAX, JMAX, NNUMHF1, NNUMHF2
  public :: N_TRUNC_M0_DFS, N_TRUNC_M1_DFS, N_TRUNC_M2_DFS, N_TRUNC_M3_DFS
  public :: N_L2TRUNC1_M0_DFS, N_L2TRUNC2_M0_DFS, N_L2TRUNC1_M1_DFS, N_L2TRUNC2_M1_DFS
  public :: N_L2TRUNC1_M2_DFS, N_L2TRUNC2_M2_DFS, N_L2TRUNC1_M3_DFS, N_L2TRUNC2_M3_DFS
  public :: NNUM1_M_DFS, NNUM2_M_DFS, NNUM_M_DFS

contains

  subroutine com_dfs__ini_interp( jcn_dfs_in, jcn_grid_in, imax_in, jmax_in, mmax_in, nmax_in )
    integer,intent(IN) :: jcn_dfs_in
    integer,intent(IN) :: jcn_grid_in
    integer,intent(IN) :: imax_in
    integer,intent(IN) :: jmax_in
    integer,intent(IN) :: mmax_in
    integer,intent(IN) :: nmax_in
    !
    integer :: m
    !
    jcn_dfs  = jcn_dfs_in
    jcn_grid = jcn_grid_in
    !
    MMAX = mmax_in
    NMAX = nmax_in
    NNUM = NMAX+1
    IMAX = imax_in
    JMAX = jmax_in
    NNUMHF1 = (NNUM+1)/2  
    NNUMHF2 = NNUM/2      !! NNUMHF1+NNUMHF2=NNUM
    !
    if ( jcn_grid == -1 .and. NMAX > JMAX ) then
       write(6,*) "Error: NMAX should be NMAX<=JMAX when jcn_grid = -1"
       write(6,*) "       NMAX, JMAX =", NMAX, JMAX
       stop 999
    end if
    if ( jcn_grid /= -1 .and. NMAX >= JMAX ) then
       write(6,*) "Error: NMAX should be NMAX<JMAX when jcn_grid = ",jcn_grid
       write(6,*) "       NMAX, JMAX =", NMAX, JMAX
       stop 999
    end if
    !
    !
    if ( jcn_dfs == 0 ) then  !! Old DFS
      N_TRUNC_M0_DFS = NNUM
      N_TRUNC_M1_DFS = NNUM
      N_TRUNC_M2_DFS = NNUM
      N_TRUNC_M3_DFS = NNUM
    else if ( jcn_dfs == 1 ) then
      if ( jcn_grid == -1 ) then
        N_TRUNC_M0_DFS = min(NNUM,JMAX)
        N_TRUNC_M1_DFS = min(NNUM-1,JMAX-1)
        N_TRUNC_M2_DFS = NNUM-2
        N_TRUNC_M3_DFS = NNUM-1
      else
        N_TRUNC_M0_DFS = NNUM
        N_TRUNC_M1_DFS = NNUM-1
        N_TRUNC_M2_DFS = NNUM-2
        N_TRUNC_M3_DFS = NNUM-1
      end if
    else  !! jcn_dfs == 2, New DFS
      if ( jcn_grid == -1 ) then
        N_TRUNC_M0_DFS = min(NNUM,JMAX)
        N_TRUNC_M1_DFS = min(NNUM-1,JMAX-1)
        N_TRUNC_M2_DFS = NNUM-2
        N_TRUNC_M3_DFS = NNUM-3
      else
        N_TRUNC_M0_DFS = NNUM
        N_TRUNC_M1_DFS = NNUM-1
        N_TRUNC_M2_DFS = NNUM-2
        N_TRUNC_M3_DFS = NNUM-3
      end if
    end if
    !
    N_L2TRUNC1_M0_DFS = (N_TRUNC_M0_DFS+1)/2
    N_L2TRUNC2_M0_DFS = N_TRUNC_M0_DFS/2
    !
    N_L2TRUNC1_M1_DFS = (N_TRUNC_M1_DFS+1)/2
    N_L2TRUNC2_M1_DFS = N_TRUNC_M1_DFS/2
    !
    N_L2TRUNC1_M2_DFS = (N_TRUNC_M2_DFS+1)/2
    N_L2TRUNC2_M2_DFS = N_TRUNC_M2_DFS/2
    !
    N_L2TRUNC1_M3_DFS = (N_TRUNC_M3_DFS+1)/2
    N_L2TRUNC2_M3_DFS = N_TRUNC_M3_DFS/2
    !
    allocate( NNUM1_M_DFS(0:MMAX) )
    allocate( NNUM2_M_DFS(0:MMAX) )
    allocate( NNUM_M_DFS(0:MMAX) )
    !
!$OMP PARALLEL default(SHARED), private(m)
!$OMP DO schedule(STATIC)
    do m=0,MMAX
      if ( m == 0 ) then
        NNUM1_M_DFS(m) = N_L2TRUNC1_M0_DFS
        NNUM2_M_DFS(m) = N_L2TRUNC2_M0_DFS
        NNUM_M_DFS(m)  = N_TRUNC_M0_DFS
      else if ( m == 1 ) then
        NNUM1_M_DFS(m) = N_L2TRUNC1_M1_DFS
        NNUM2_M_DFS(m) = N_L2TRUNC2_M1_DFS
        NNUM_M_DFS(m)  = N_TRUNC_M1_DFS
      else if ( mod(m,2) == 0 ) then
        NNUM1_M_DFS(m) = N_L2TRUNC1_M2_DFS
        NNUM2_M_DFS(m) = N_L2TRUNC2_M2_DFS
        NNUM_M_DFS(m)  = N_TRUNC_M2_DFS
      else
        NNUM1_M_DFS(m) = N_L2TRUNC1_M3_DFS
        NNUM2_M_DFS(m) = N_L2TRUNC2_M3_DFS
        NNUM_M_DFS(m)  = N_TRUNC_M3_DFS
      end if
    end do
!$OMP END DO
!$OMP END PARALLEL


!    write(0,*) 'N_TRUNC_M0_DFS=',N_TRUNC_M0_DFS
!    write(0,*) 'N_TRUNC_M1_DFS=',N_TRUNC_M1_DFS
!    write(0,*) 'N_TRUNC_M2_DFS=',N_TRUNC_M2_DFS
!    write(0,*) 'N_TRUNC_M3_DFS=',N_TRUNC_M3_DFS
!    write(0,*)
    !
!    write(0,*) 'N_L2TRUNC1_M0_DFS=',N_L2TRUNC1_M0_DFS
!    write(0,*) 'N_L2TRUNC2_M0_DFS=',N_L2TRUNC2_M0_DFS

!    write(0,*) 'N_L2TRUNC1_M1_DFS=',N_L2TRUNC1_M1_DFS
!    write(0,*) 'N_L2TRUNC2_M1_DFS=',N_L2TRUNC2_M1_DFS

!    write(0,*) 'N_L2TRUNC1_M2_DFS=',N_L2TRUNC1_M2_DFS
!    write(0,*) 'N_L2TRUNC2_M2_DFS=',N_L2TRUNC2_M2_DFS

!    write(0,*) 'N_L2TRUNC1_M3_DFS=',N_L2TRUNC1_M3_DFS
!    write(0,*) 'N_L2TRUNC2_M3_DFS=',N_L2TRUNC2_M3_DFS
!    write(0,*)

    !
  end subroutine com_dfs__ini_interp

end module com_dfs_interp

