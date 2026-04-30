module legendre

  use e_time, only : e_time__start, e_time__end

!---------------------------------------------------------------
  implicit none
!---------------------------------------------------------------
  private
  
  public :: legendre__ini, legendre__g2w, legendre__g2w_uv
  public :: legendre__w2g, legendre__w2g_uv, legendre__w2g_dy
!---------------------------------------------------------------
  
  integer(8) :: IMAX
  integer(8) :: JMAX
  integer(8) :: NNUM
  integer(8) :: MNUM
  integer(8) :: NTRNC
  integer(8) :: MTRNC
  real(8),save,allocatable :: PZ(:,:)
  real(8),save,allocatable :: RM(:,:)
  real(8),save,allocatable :: PM(:,:,:)
  integer(8),save,allocatable :: JC(:,:)
  real(8),save,allocatable :: CM(:,:)
  
  integer(8),save :: MNWAV
  integer(8),save :: MNWAV_UV
  integer(8),save,allocatable :: MNSTART(:)
  integer(8),save,allocatable :: MNSTART_UV(:)
  
  real(8),save,allocatable :: emn0(:)
  real(8),save,allocatable :: emn1(:)
!
  logical,save :: initialized = .false.
!
!******************************************************************
contains
!******************************************************************

subroutine legendre__ini                                   &
 &( imax_in, jmax_in, nnum_in, mnum_in,                    &!IN
 &  mnwav_out, mnwav_uv_out, mnstart_out, mnstart_uv_out,  &!OUT
 &  weight, alat, sinlat, coslat, coslat_inv )                  !OUT
!
  integer,intent(in) :: imax_in
  integer,intent(in) :: jmax_in
  integer,intent(in) :: nnum_in
  integer,intent(in) :: mnum_in
  integer,intent(out) :: mnwav_out
  integer,intent(out) :: mnwav_uv_out
  integer,intent(out) :: mnstart_out(0:mnum_in)
  integer,intent(out) :: mnstart_uv_out(0:mnum_in)
  real(8),intent(out) :: weight(jmax_in)
  real(8),intent(out) :: alat(jmax_in)
  real(8),intent(out) :: sinlat(jmax_in)
  real(8),intent(out) :: coslat(jmax_in)
  real(8),intent(out) :: coslat_inv(jmax_in)
!
  integer(8) :: m, num, j, mn_uv
  real(8) :: am, an
!
  if ( initialized ) then
     write( 6,* ) 'legendre_ini : already initialized.'
     stop 999
  end if
!
  initialized = .true.
  
  IMAX = imax_in
  JMAX = jmax_in
  NNUM = nnum_in
  MNUM = mnum_in
  MNWAV = MNUM*( NNUM + (NNUM-MNUM+1) )/2
  MNWAV_UV = MNWAV + MNUM

  mnwav_out    = MNWAV
  mnwav_uv_out = MNWAV_UV
!
!  allocate( pnm (mnwav,jmax/2) )
!  allocate( dpnm(mnwav,jmax/2) )
!  allocate( pnm1(mnum,jmax/2) )
!  allocate( weight  (jmax) )
!
!  pnm (:,:) = pnm_in (:,:)
!  dpnm(:,:) = dpnm_in(:,:)
!  pnm1(:,:) = pnm1_in(:,:)
!  weight(:)     = weight_in(:)
!

  MTRNC = MNUM-1
  NTRNC = NNUM-1
  write(6,*) "NTRNC,MTRNC,NNUM,MNUM=",NTRNC,MTRNC,NNUM,MNUM
  
  allocate( PZ(JMAX/2,5) )
  allocate( RM(NNUM/2*3+NNUM+1,0:MTRNC) )
  allocate( PM(JMAX/2,2,0:MTRNC) )
  allocate( JC(NNUM/8+1,0:MTRNC) )
  allocate( MNSTART(0:MTRNC+1) )
  allocate( MNSTART_UV(0:MTRNC+1) )
  allocate( CM(2*NTRNC+1,0:MTRNC) )
  
  call LXINIG( JMAX, PZ, 1_8 )
  
!  write(6,*) "PZ(:,1)=",PZ(:,1)
!  write(6,*) "PZ(:,2)=",PZ(:,2)
!  write(6,*) "PZ(:,3)=",PZ(:,3)
  
!$OMP PARALLEL default(SHARED), private(m)
 !$OMP DO schedule(DYNAMIC)
  do m=0,MTRNC
     call LXINIR( NNUM, m, RM(:,m) )
     call LXINIW( NNUM, JMAX, m, PZ, PM(:,:,m), RM(:,m), JC(:,m) )
     call LXINIC( NTRNC, m, CM(:,m) )
  end do
 !$OMP END DO
!$OMP END PARALLEL

  MNSTART(0) = 1
  MNSTART_UV(0) = 1
  do m=0,MTRNC
     num = NTRNC - m + 1            !! Number of m ~ NTRNC
     MNSTART(m+1) = MNSTART(m) + num
     MNSTART_UV(m+1) = MNSTART_UV(m) + num + 1
  end do
  
  do m=0,MTRNC+1
     mnstart_out(m) = MNSTART(m)
     mnstart_uv_out(m) = MNSTART_UV(m)
  end do
  
!!$OMP PARALLEL default(SHARED), private(j)
! !$OMP DO schedule(STATIC)
!  do j=1,JMAX/2
!     weight(JMAX/2+j)           = PZ(j,2)*2.0d0
!     weight(JMAX/2+1-j)         = weight(JMAX/2+j)
!     sinlat(JMAX/2+j)       = PZ(j,1)
!     sinlat(JMAX/2+1-j)     = -sinlat(JMAX/2+j)
!     coslat(JMAX/2+j)       = PZ(j,3)
!     coslat(JMAX/2+1-j)     = coslat(JMAX/2+j)
!     coslat_inv(JMAX/2+j)   = PZ(j,4)
!     coslat_inv(JMAX/2+1-j) = coslat_inv(JMAX/2+j)
!     alat(JMAX/2+j)         = atan( sinlat(JMAX/2+j), coslat(JMAX/2+j) )
!     alat(JMAX/2+1-j)       = -alat(JMAX/2+j)
!  end do
! !$OMP END DO
!!$OMP END PARALLEL
!  
!  weight(:)         = weight(JMAX:1:-1)
!  sinlat(:)     = sinlat(JMAX:1:-1)
!  coslat(:)     = coslat(JMAX:1:-1)
!  coslat_inv(:) = coslat_inv(JMAX:1:-1)
!  alat(:)       = alat(JMAX:1:-1)
  
  
  
!$OMP PARALLEL default(SHARED), private(j)
 !$OMP DO schedule(STATIC)
  do j=1,JMAX/2
!xx     weight(JMAX/2+1-j)         = PZ(j,2)*2.0d0
     weight(JMAX/2+1-j)         = PZ(j,2)
     weight(JMAX/2+j)           = weight(JMAX/2+1-j)
     sinlat(JMAX/2+1-j)     = PZ(j,1)
     sinlat(JMAX/2+j)       = -sinlat(JMAX/2+1-j)
     coslat(JMAX/2+1-j)     = PZ(j,3)
     coslat(JMAX/2+j)       = coslat(JMAX/2+1-j)
     coslat_inv(JMAX/2+1-j) = PZ(j,4)
     coslat_inv(JMAX/2+j)   = coslat_inv(JMAX/2+1-j)
     alat(JMAX/2+1-j)       = atan( sinlat(JMAX/2+1-j), coslat(JMAX/2+1-j) )
     alat(JMAX/2+j)         = -alat(JMAX/2+1-j)
  end do
 !$OMP END DO
!$OMP END PARALLEL
  
  
!  write(6,*) "weight(:)        =",weight(:)
!  write(6,*) "alat(:)      =",alat(:)
!  write(6,*) "sinlat(:)    =",sinlat(:)
!  write(6,*) "coslat(:)    =",coslat(:)
!  write(6,*) "coslat_inv(:)=",coslat_inv(:)

  
  allocate( emn0(mnwav+mnum) )
  allocate( emn1(mnwav+mnum) )
  
!$OMP PARALLEL default(SHARED), private(m,mn_uv,am,an)
 !$OMP DO schedule(DYNAMIC)
  do m=0,MTRNC
    do mn_uv=MNSTART_UV(m),MNSTART_UV(m+1)-1
      am=m
      an=m+mn_uv-MNSTART_UV(m)
      emn0(mn_uv) = (an+2)*sqrt(((an+1)**2-am**2)/(4*(an+1)**2-1))
      emn1(mn_uv) = (an-1)*sqrt((an**2-am**2)/(4*an**2-1))
    end do
  end do
 !$OMP END DO
!$OMP END PARALLEL

end subroutine legendre__ini


!******************************************************************


subroutine legendre__g2w &
 &( gdata,             &!IN
 &  qdata   )           !INOUT
!------------------------------------------------------
  real(8),intent(in) :: gdata(2,0:imax/2-1,jmax)
  real(8),intent(inout) :: qdata(2,mnwav)
  
  integer(8) :: m
!
!-------------------------------------------------------
!
  call e_time__start(21,"legendre__g2w")

!
  if ( .not. initialized ) then
    write( 6,* ) 'legendre_run_g2w : not initialized yet.'
    stop 999
  end if
  
!$OMP PARALLEL default(SHARED), private(m)
 !$OMP DO schedule(DYNAMIC)
  do m=0,MTRNC
     if ( m == 0 ) then
        call LXTGZS(NNUM,NTRNC,JMAX,qdata(1,MNSTART(m):MNSTART(m+1)-1), &
         &          gdata(1,m,JMAX:1:-1),PZ,RM(:,m),0_8 )
        qdata(2,MNSTART(m):MNSTART(m+1)-1) = 0.0d0
     else
        call LXTGWS(NNUM,NTRNC,JMAX,m,qdata(:,MNSTART(m):MNSTART(m+1)-1), &
         &          gdata(:,m,JMAX:1:-1),PZ,PM(:,:,m),RM(:,m),JC(:,m),0_8 )
     end if
  end do
 !$OMP END DO
!$OMP END PARALLEL
!
  call e_time__end(21,"legendre__g2w")
  
end subroutine legendre__g2w


!******************************************************************


subroutine legendre__g2w_uv &
 &( gdata,     &!IN
 &  qdata    )  !OUT
!------------------------------------------------------
  real(8),intent(in) :: gdata(2,0:imax/2-1,jmax)
  real(8),intent(inout) :: qdata(2,mnwav_uv)
  
  integer(8) :: m
!
!-------------------------------------------------------
!
  call e_time__start(22,"legendre__g2w_uv")
!
  if ( .not. initialized ) then
    write( 6,* ) 'legendre_run_g2w : not initialized yet.'
    stop 999
  end if
  
!$OMP PARALLEL default(SHARED), private(m)
 !$OMP DO schedule(DYNAMIC)
  do m=0,MTRNC
     if ( m == 0 ) then
        call LXTGZS(NNUM,NTRNC+1,JMAX,qdata(1,MNSTART_UV(m):MNSTART_UV(m+1)-1), &
         &          gdata(1,m,JMAX:1:-1),PZ,RM(:,m),0_8 )
        qdata(2,MNSTART_UV(m):MNSTART_UV(m+1)-1) = 0.0d0
     else
        call LXTGWS(NNUM,NTRNC+1,JMAX,m,qdata(:,MNSTART_UV(m):MNSTART_UV(m+1)-1), &
         &          gdata(:,m,JMAX:1:-1),PZ,PM(:,:,m),RM(:,m),JC(:,m),0_8 )
     end if
  end do
 !$OMP END DO
!$OMP END PARALLEL

  call e_time__end(22,"legendre__g2w_uv")
!
end subroutine legendre__g2w_uv


!******************************************************************


subroutine legendre__w2g &
 &( qdata,             &!IN
 &  gdata        )      !OUT

  real(8),intent(in) :: qdata(2,mnwav)
  real(8),intent(out) :: gdata(2,0:imax/2-1,jmax)
  
  integer(8) :: m,j
!
!-------------------------------------------------------
!
  call e_time__start(23,"legendre__w2g")
!
  if ( .not. initialized ) then
    write( 6,* ) 'legendre_run_w2g : not initialized yet.'
    stop 999
  end if
  
!$OMP PARALLEL default(SHARED), private(m,j)
 !$OMP DO schedule(DYNAMIC)
  do m=0,MTRNC
     if ( m == 0 ) then
        call LXTSZG(NNUM,NTRNC,JMAX,qdata(1,MNSTART(m):MNSTART(m+1)-1), &
         &          gdata(1,m,JMAX:1:-1),PZ,RM(:,m),0_8 )
        gdata(2,m,:) = 0.0d0
     else
        call LXTSWG(NNUM,NTRNC,JMAX,m,qdata(:,MNSTART(m):MNSTART(m+1)-1), &
        &          gdata(:,m,JMAX:1:-1),PZ,PM(:,:,m),RM(:,m),JC(:,m),0_8 )
     end if
  end do
 !$OMP END DO

 !$OMP DO schedule(STATIC)
  do j=1,JMAX
    gdata(:,MTRNC+1:imax/2-1,j) = 0.0d0
  end do
 !$OMP END DO
!$OMP END PARALLEL
!
  call e_time__end(23,"legendre__w2g")
!
end subroutine legendre__w2g


!******************************************************************


subroutine legendre__w2g_dy &
 &( qdata,             &!IN
 &  gdata   )           !OUT
 
  real(8),intent(in)  :: qdata(2,mnwav)
  real(8),intent(out) :: gdata(2,0:imax/2-1,jmax)

  real(8) :: qwork(2,0:NTRNC+1)
!
  integer :: j
  integer(8) :: m, mn, mn_uv, mn_uv_s
!

  call e_time__start(24,"legendre__w2g_dy")
!
  if ( .not. initialized ) then
    write( 6,* ) 'legendre_run_w2g_dy : not initialized yet.'
    stop 999
  end if
  
  
  
!  qdata=0.0d0 !!!!!!!!!!!!!!!!11111
!  qdata(:,1:10)=1.0d0
!  qdata(:,100:mnwav)=1.0d0
   
!
!  jj1  = jblk2/2 *(jl-1)+1
!  call w2gdd_hy_run                                   &
!   &( mnum, nend1, jend1, mnwav, kqs, kqe, kqmax,  &!IN
!   &  imax , jblk2/2, jblk2,                        &!IN
!   &  dpnm(1,jj1),                                  &!IN
!   &  qdata,                                        &!IN
!   &  gdata  )                                       !OUT
!  gdata( :,:, jmax/2+1:jmax, : ) = gdata( :,:, jmax:jmax/2+1:-1, : )
  
  
  
!#ifdef AAAAAAA
!  write(6,*) "1: gdata(1,0,1:10,1)=",gdata(1,0,1:10,1)
!  write(6,*) "1: gdata(1,1,1:10,1)=",gdata(1,1,1:10,1)
!  write(6,*) "1: gdata(1,:,1,1)=",gdata(1,:,1,1)
!  write(6,*) "1: gdata(2,:,1,1)=",gdata(2,:,1,1)
!  write(6,*) "1: gdata(2,:,jmax/2,1)=",gdata(2,:,jmax/2,1)

!$OMP PARALLEL default(SHARED), private(m,mn,mn_uv,mn_uv_s,qwork,j)
 !$OMP DO schedule(DYNAMIC)
  do m=0,MTRNC
     mn_uv_s = MNSTART_UV(m)
     mn_uv = MNSTART_UV(m)
     mn    = MNSTART(m)
!     if ( mn+1 <= MNSTART(m+1)-1 ) then
     if ( m /= MTRNC) then
        qwork(1,mn_uv-mn_uv_s) = emn0(mn_uv)*qdata(1,mn+1)
        qwork(2,mn_uv-mn_uv_s) = emn0(mn_uv)*qdata(2,mn+1)
     else
        qwork(1,mn_uv-mn_uv_s) = 0.0d0
        qwork(2,mn_uv-mn_uv_s) = 0.0d0
     end if
     do mn_uv=MNSTART_UV(m)+1,MNSTART_UV(m+1)-3
        mn = mn_uv - MNSTART_UV(m) + MNSTART(m)
        qwork(1,mn_uv-mn_uv_s) = emn0(mn_uv)*qdata(1,mn+1) - emn1(mn_uv)*qdata(1,mn-1)
        qwork(2,mn_uv-mn_uv_s) = emn0(mn_uv)*qdata(2,mn+1) - emn1(mn_uv)*qdata(2,mn-1)
     end do
     mn_uv = MNSTART_UV(m+1) - 2
!     if ( mn_uv >= MNSTART_UV(m)+1 ) then
     if ( m /= MTRNC) then
        mn = mn_uv - MNSTART_UV(m) + MNSTART(m)
        qwork(1,mn_uv-mn_uv_s) = - emn1(mn_uv)*qdata(1,mn-1)
        qwork(2,mn_uv-mn_uv_s) = - emn1(mn_uv)*qdata(2,mn-1)
     end if
     mn_uv = MNSTART_UV(m+1) - 1
!     if ( mn_uv >= MNSTART_UV(m)+1 ) then
        mn = mn_uv - MNSTART_UV(m) + MNSTART(m)
        qwork(1,mn_uv-mn_uv_s) = - emn1(mn_uv)*qdata(1,mn-1)
        qwork(2,mn_uv-mn_uv_s) = - emn1(mn_uv)*qdata(2,mn-1)
!     end if


     if ( m == 0 ) then
        call LXTSZG( NNUM, NTRNC+1, JMAX, qwork(1,:),            &
         &           gdata(1,m,JMAX:1:-1), PZ, RM(:,m), 0_8 )
        gdata(2,m,JMAX:1:-1) = 0.0d0
     else
        call LXTSWG( NNUM, NTRNC+1, JMAX, m, qwork,                               &
         &           gdata(:,m,JMAX:1:-1), PZ, PM(:,:,m), RM(:,m), JC(:,m), 0_8 )
     end if
  end do
 !$OMP END DO

 !$OMP DO schedule(STATIC)
  do j=1,JMAX
    gdata(:,MTRNC+1:imax/2-1,j) = 0.0d0
  end do
 !$OMP END DO
!$OMP END PARALLEL
  
!  gdata=-gdata        !!!!!!
  
!  write(6,*) "2: gdata(1,0,1:10,1)=",gdata(1,0,1:10,1)
!  write(6,*) "2: gdata(1,1,1:10,1)=",gdata(1,1,1:10,1)
!  write(6,*) "2: gdata(1,:,1,1)=",gdata(1,:,1,1)
!  write(6,*) "2: gdata(2,:,1,1)=",gdata(2,:,1,1)
!  write(6,*) "2: gdata(2,:,jmax/2,1)=",gdata(2,:,jmax/2,1)
!  stop 333
  
  
!#ifdef AAAAAAAaa  
!  write(6,*) "1: gdata(1,0,1:10,1)=",gdata(1,0,1:10,1)
!  write(6,*) "1: gdata(1,1,1:10,1)=",gdata(1,1,1:10,1)

!  do m=0,MTRNC
!     if ( m == 0 ) then
!        call LXCYZS( NTRNC, qdata(1,MNSTART(m):MNSTART(m+1)-1), qwork1, CM(:,m) )
!        call LXTSZG( NNUM, NTRNC+1, JMAX, qwork1,            &
!         &           gdata(1,m,JMAX:1:-1,1), PZ, RM(:,m), 0_8 )
!        gdata(2,m,JMAX:1:-1,1) = 0.0d0
!     else
!        call LXCYWS( NTRNC, m, qdata(:,MNSTART(m):MNSTART(m+1)-1), qwork2, CM(:,m) )
!        call LXTSWG( NNUM, NTRNC+1, JMAX, m, qwork2,                               &
!         &           gdata(:,m,JMAX:1:-1,1), PZ, PM(:,:,m), RM(:,m), JC(:,m), 0_8 )
!     end if
!  end do
!  gdata(:,MTRNC+1:imax/2-1,:,:) = 0.0d0

!  write(6,*) "2: gdata(1,0,1:10,1)=",gdata(1,0,1:10,1)
!  write(6,*) "2: gdata(1,1,1:10,1)=",gdata(1,1,1:10,1)
!  stop 333
!#endif
  
  call e_time__end(24,"legendre__w2g_dy")
!
end subroutine legendre__w2g_dy


!******************************************************************


subroutine legendre__w2g_uv &
 &( qdata,     &!IN
 &  gdata  )    !OUT

  real(8),intent(in) :: qdata(2,mnwav_uv)
  real(8),intent(out) :: gdata(2,0:imax/2-1,jmax)
!
!  integer :: jj1
!  integer :: nend1, jend1, kqmax, kqs, kqe, jblk2, jl, kqu
  
  integer(8) :: m,j
!
!-------------------------------------------------------
!
  call e_time__start(25,"legendre__w2g_uv")
  
!  nend1 = mnum
!  jend1 = mnum
!  kqmax = kmax*2
!  kqs = 1
!  kqe = kmax*2
!  jblk2 = jmax
!  jl = 1
!  kqu = 0
!
  if ( .not. initialized ) then
    write( 6,* ) 'legendre_run_w2g : not initialized yet.'
    stop 999
  end if
!


!  qdata(:,:) = 1.0d0 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11
!  qdata0=1.0d0
!  !qucor=0.0d0
!  qucor=1.0d0


!
!  jj1  = jblk2/2 *(jl-1)+1
!    call w2g_hy_run                                   &
!     &( mnum, nend1, jend1, mnwav, kqs, kqe, kqmax,  &!IN
!     &  imax , jblk2/2, jblk2,                        &!IN
!     &  pnm(1,jj1),                                   &!IN
!     &  qdata0,                                        &!IN
!     &  gdata,                                        &!OUT
!     &  kmax,                                         &!IN
!     &  pnm1(1,jj1), kqu, qucor )                      !IN,optional
!  gdata( :,:, jmax/2+1:jmax, : ) = gdata( :,:, jmax:jmax/2+1:-1, : )
!
  
!  write(6,*) "1: gdata(1,0,1:10,1)=",gdata(1,0,1:10,1)
!  write(6,*) "1: gdata(1,1,1:10,1)=",gdata(1,1,1:10,1)
!  write(6,*) "1: gdata(1,2,1:10,1)=",gdata(1,2,1:10,1)

!  write(6,*) "1: gdata(1,0,JMAX-10:JMAX,1)=",gdata(1,0,JMAX-10:JMAX,1)
!  write(6,*) "1: gdata(1,1,JMAX-10:JMAX,1)=",gdata(1,1,JMAX-10:JMAX,1)
!  write(6,*) "1: gdata(1,2,JMAX-10:JMAX,1)=",gdata(1,2,JMAX-10:JMAX,1

!  write(6,*) "1: gdata(1,0,:,1)=",gdata(1,0,:,1)
  
  
!  gdata(:,:,:,:) = 0.0d0

!$OMP PARALLEL default(SHARED), private(m,j)
 !$OMP DO schedule(DYNAMIC)
  do m=0,MTRNC
     if ( m == 0 ) then
        call LXTSZG(NNUM,NTRNC+1,JMAX,qdata(1,MNSTART_UV(m):MNSTART_UV(m+1)-1), &
         &          gdata(1,m,JMAX:1:-1),PZ,RM(:,m),0_8 )
        gdata(2,m,JMAX:1:-1) = 0.0d0
     else
        call LXTSWG(NNUM,NTRNC+1,JMAX,m,qdata(:,MNSTART_UV(m):MNSTART_UV(m+1)-1), &
        &          gdata(:,m,JMAX:1:-1),PZ,PM(:,:,m),RM(:,m),JC(:,m),0_8 )
     end if
  end do
 !$OMP END DO

 !$OMP DO schedule(STATIC)
  do j=1,JMAX
    gdata(:,MTRNC+1:imax/2-1,j) = 0.0d0
  end do
 !$OMP END DO
!$OMP END PARALLEL

  
!  write(6,*) "2: gdata(1,0,1:10,1)=",gdata(1,0,1:10,1)
!  write(6,*) "2: gdata(1,1,1:10,1)=",gdata(1,1,1:10,1)
!  write(6,*) "2: gdata(1,2,1:10,1)=",gdata(1,2,1:10,1)
  
!  write(6,*) "2: gdata(1,0,JMAX-10:JMAX,1)=",gdata(1,0,JMAX-10:JMAX,1)
!  write(6,*) "2: gdata(1,1,JMAX-10:JMAX,1)=",gdata(1,1,JMAX-10:JMAX,1)
!  write(6,*) "2: gdata(1,2,JMAX-10:JMAX,1)=",gdata(1,2,JMAX-10:JMAX,1

!  write(6,*) "2: gdata(1,0,:,1)=",gdata(1,0,:,1)
 
!  stop 333



!  do m=0,MTRNC
!     write(6,*) "m=",m
!     write(6,*) "qdata(1,MNSTART_UV(m):MNSTART_UV(m+1)-1)        =",qdata(1,MNSTART_UV(m):MNSTART_UV(m+1)-1)
!     write(6,*) "qdata0(1,MNSTART(m):MNSTART(m+1)-1),qucor(1,m+1)=",qdata0(1,MNSTART(m):MNSTART(m+1)-1),qucor(1,m+1)
!  end do
  
!  stop 333

  call e_time__end(25,"legendre__w2g_uv")
  
end subroutine legendre__w2g_uv


!******************************************************************


end module legendre
