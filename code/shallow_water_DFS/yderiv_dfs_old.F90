module yderiv_dfs_old

  use com_dfs, only : NNUM, MMAX,                               &
   &                  N_TRUNC_M0_DFS, N_TRUNC_M1_DFS, N_TRUNC_M2_DFS,                           &
   &                  N_L2TRUNC1_M1_DFS, N_L2TRUNC2_M1_DFS, N_L2TRUNC1_M2_DFS, N_L2TRUNC2_M2_DFS

  implicit none

  private
  public :: yderiv_dfs_old__run1, yderiv_dfs_old__run2

contains

  subroutine yderiv_dfs_old__run1(ne,wdata,wdatay) !IN,IN,OUT
    !
    integer(kind=4),intent(IN) :: ne                              !! Number of elements
    complex(8),intent(IN) :: wdata(NNUM,0:MMAX)  !! Wave data
    complex(8),intent(OUT) :: wdatay(NNUM,0:MMAX)!! cos(lat)*@wdata/@lat
    !
    complex(8) :: work(NNUM)
    complex(8) :: wyp1
    integer :: m,n,nn
    !
    !
!$OMP PARALLEL default(SHARED), private(m,n,nn,wyp1,work)
!$OMP DO schedule(STATIC)
      do m=0,MMAX
        if ( m == 0 ) then               !! mg is zero
           do n=1,NNUM-1               !! Wave number+1
             nn=n+1
             work(n) = -n*wdata(nn,m)  !! @(psi or chi)/@phi, cos(n*phi)=>sin(n*phi)
           end do
           work(NNUM)=0.0d0
           !
           !! sin(n*phi) => sin(phi)*cos(n*phi)
           n=NNUM
           wdatay(n,m) = work(n)*2.0
           n=NNUM-1
           wdatay(n,m) = work(n)*2.0
           do n=NNUM-2,2,-1  !! Wave number+1
             wdatay(n,m) = work(n)*2.0d0 + wdatay(n+2,m)
           end do
           n=1                  !! Wave number+1
           wdatay(n,m) = work(n) + wdatay(n+2,m)*0.5d0
        else if ( mod(m,2) == 0 ) then   !! m is even ( not zero )
          n=1   !! Wave number
          wdatay(n,m) = -n*wdata(n+1,m)*0.5d0
          do n=2,NNUM-1   !! Wave number
            wdatay(n,m) = n*( - wdata(n+1,m) + wdata(n-1,m) )*0.5d0
          end do
          n=NNUM   !! Wave number
          wdatay(n,m) = n*wdata(n-1,m)*0.5d0
          !
        else                              !! m is odd
          n=1       !! Wave number
          wdatay(n,m) = -(n+1)*wdata(n+1,m)*0.5d0
          do n=2,NNUM-1   !! Wave number
            wdatay(n,m) = ( - (n+1)*wdata(n+1,m) + (n-1)*wdata(n-1,m) )*0.5d0
          end do
          n=NNUM   !! Wave number
          wdatay(n,m) = (n-1)*wdata(n-1,m)*0.5d0
          !
        end if
        wdatay(:,m) = -wdatay(:,m)
      end do
!$OMP END DO
!$OMP END PARALLEL
    !
  end subroutine yderiv_dfs_old__run1


  subroutine yderiv_dfs_old__run2(ne,wdata,wdatay) !IN,IN,OUT
    !
    integer(kind=4),intent(IN) :: ne                              !! Number of elements
    complex(8),intent(IN) :: wdata(NNUM,0:MMAX)  !! Wave data
    complex(8),intent(OUT) :: wdatay(NNUM,0:MMAX)!! cos(lat)*@wdata/@lat
    !
    complex(8) :: wyp1, ww
    integer :: m,n,nn
    !
    !
!$OMP PARALLEL default(SHARED), private(m,n,nn,wyp1,ww)
!$OMP DO schedule(STATIC)
      do m=0,MMAX
        if ( m == 0 ) then               !! m is zero
          n=1
          nn = n-1       !! Wave number
          wdatay(n,m) = - (nn-1)*wdata(n+1,m)*0.5d0
          n=2
          nn = n-1       !! Wave number
          wdatay(n,m) = 2.0d0*wdata(1,m)
          do n=3,NNUM-1
            nn = n-1     !! Wave number
            wdatay(n,m) = ( - (nn-1)*wdata(n+1,m) + (nn+1)*wdata(n-1,m) )*0.5d0  
          end do
          n=NNUM
          nn = n-1       !! Wave number
          wdatay(n,m) = (nn+1)*wdata(n-1,m)*0.5d0
          n=NNUM+1
          nn = n-1       !! Wave number
          wyp1 = (nn+1)*wdata(n-1,m)*0.5d0
          !
        else if ( mod(m,2) == 0 ) then   !! m is even ( not zero )
          n=1
          nn = n     !! Wave number
          wdatay(n,m) = -nn*wdata(n+1,m)*0.5d0
          do n=2,NNUM-1
            nn = n   !! Wave number
            wdatay(n,m) = nn*( - wdata(n+1,m) + wdata(n-1,m) )*0.5d0
          end do
          n=NNUM
          nn = n     !! Wave number
          wdatay(n,m) = nn*wdata(n-1,m)*0.5d0
          !
        else                              !! m is odd
          n=1
          nn = n     !! Wave number
          wdatay(n,m) = -(nn+1)*wdata(n+1,m)*0.5d0
          do n=2,NNUM-1
            nn = n   !! Wave number
            wdatay(n,m) = ( - (nn+1)*wdata(n+1,m) + (nn-1)*wdata(n-1,m) )*0.5d0
          end do
          n=NNUM
          nn = n     !! Wave number
          wdatay(n,m) = (nn-1)*wdata(n-1,m)*0.5d0
          !
        end if    
        !
        wdatay(:,m) = -wdatay(:,m)
      end do
!$OMP END DO
!$OMP END PARALLEL
    !
  end subroutine yderiv_dfs_old__run2

end module yderiv_dfs_old
