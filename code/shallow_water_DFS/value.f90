 !$OMP PARALLEL default(SHARED), private(m,j,nt,am,an)
        !$OMP DO schedule(STATIC)
        do m=0,MNMAX
            am=m
            if ( m >= 3 .and. mod(m,2) == 1 ) then
              j=1
              an=j-1
              myqu(j,m) = ER_INV*0.5d0*dcmplx(0.0d0,1.0d0)*am*myqchi(j,m) & !! n=0
               &        - ER_INV*0.25d0*an*myqpsi(j+1,m)
              myqv(j,m) = ER_INV*0.5d0*dcmplx(0.0d0,1.0d0)*am*myqpsi(j,m)  & !! n=0
               &        + ER_INV*0.25d0*an*myqchi(j+1,m)
              j=2
              an=j-1
              myqu(j,m) = ER_INV*0.5d0*dcmplx(0.0d0,1.0d0)*am*myqchi(j,m)    & !! n=1
               &        + ER_INV*0.25d0*an*( 3.0d0*myqpsi(j-1,m) - myqpsi(j+1,m) )
              myqv(j,m) = ER_INV*0.5d0*dcmplx(0.0d0,1.0d0)*am*myqpsi(j,m)    & !! n=1
               &        + ER_INV*0.25d0*an*( -3.0d0*myqchi(j-1,m) + myqchi(j+1,m) )
              j=3
              an=j-1
              myqu(j,m) = ER_INV*0.5d0*dcmplx(0.0d0,1.0d0)*am*( -myqchi(j-2,m) + myqchi(j,m) ) &
               &        + ER_INV*0.25d0*an*( 2.0d0*myqpsi(j-1,m) - myqpsi(j+1,m) )
              myqv(j,m) = ER_INV*0.5d0*dcmplx(0.0d0,1.0d0)*am*( -myqpsi(j-2,m) + myqpsi(j,m) ) &
               &        + ER_INV*0.25d0*an*( -2.0d0*myqchi(j-1,m) + myqchi(j+1,m) )
              do j=4,MNUM-1
                an=j-1
                myqu(j,m) = ER_INV*0.5d0*dcmplx(0.0d0,1.0d0)*am*( -myqchi(j-2,m) + myqchi(j,m) ) &
                 &        + ER_INV*0.25d0*an*( -myqpsi(j-3,m) + 2.0d0*myqpsi(j-1,m) - myqpsi(j+1,m) )
                myqv(j,m) = ER_INV*0.5d0*dcmplx(0.0d0,1.0d0)*am*( -myqpsi(j-2,m) + myqpsi(j,m) ) &
                 &        + ER_INV*0.25d0*an*( myqchi(j-3,m) -2.0d0*myqchi(j-1,m) + myqchi(j+1,m) )
              end do
              j=MNUM
              an=j-1
              myqu(j,m) = ER_INV*0.5d0*dcmplx(0.0d0,1.0d0)*am*( -myqchi(j-2,m) + myqchi(j,m) ) &
               &        + ER_INV*0.25d0*an*( -myqpsi(j-3,m) + 2.0d0*myqpsi(j-1,m) )
              myqv(j,m) = ER_INV*0.5d0*dcmplx(0.0d0,1.0d0)*am*( -myqpsi(j-2,m) + myqpsi(j,m) ) &
               &        + ER_INV*0.25d0*an*( myqchi(j-3,m) -2.0d0*myqchi(j-1,m) )
            else if ( m == 0 ) then
              do j=1,MNUM-1
                an=j
                myqu(j,m) = -ER_INV*an*myqpsi(j+1,m) 
                myqv(j,m) =  ER_INV*an*myqchi(j+1,m) 
              end do
              j=MNUM
              myqu(j,m) = 0.0d0
              myqv(j,m) = 0.0d0
            else if ( mod(m,2) == 1 ) then
              j=1
              an=j-1
              myqu(j,m) = ER_INV*dcmplx(0.0d0,1.0d0)*am*myqchi(j,m)  !! n=0
              myqv(j,m) = ER_INV*dcmplx(0.0d0,1.0d0)*am*myqpsi(j,m)  !! n=0
              j=2
              an=j-1
              myqu(j,m) = ER_INV*dcmplx(0.0d0,1.0d0)*am*myqchi(j,m)           & !! n=1
               &        + ER_INV*0.5d0*an*( 2.0d0*myqpsi(j-1,m) - myqpsi(j+1,m) )
              myqv(j,m) = ER_INV*dcmplx(0.0d0,1.0d0)*am*myqpsi(j,m)           & !! n=1
               &        + ER_INV*0.5d0*an*( -2.0d0*myqchi(j-1,m) + myqchi(j+1,m) )
              do j=3,MNUM-1
                an=j-1
                myqu(j,m) = ER_INV*dcmplx(0.0d0,1.0d0)*am*myqchi(j,m)  &
                 &        + ER_INV*0.5d0*an*( myqpsi(j-1,m) - myqpsi(j+1,m) )
                myqv(j,m) = ER_INV*dcmplx(0.0d0,1.0d0)*am*myqpsi(j,m) &
                 &        + ER_INV*0.5d0*an*( -myqchi(j-1,m) + myqchi(j+1,m) )
              end do
              j=MNUM
              an=j-1
              myqu(j,m) = ER_INV*dcmplx(0.0d0,1.0d0)*am*myqchi(j,m)      &
               &        + ER_INV*0.5d0*an*myqpsi(j-1,m) 
              myqv(j,m) = ER_INV*dcmplx(0.0d0,1.0d0)*am*myqpsi(j,m)      &
               &        - ER_INV*0.5d0*an*myqchi(j-1,m)
            else  !! mod(m,2) == 0
              j=1
              an=j
              myqu(j,m) = ER_INV*dcmplx(0.0d0,1.0d0)*am*myqchi(j,m)  &
               &        - ER_INV*0.5d0*an*myqpsi(j+1,m)
              myqv(j,m) = ER_INV*dcmplx(0.0d0,1.0d0)*am*myqpsi(j,m)  &
               &        + ER_INV*0.5d0*an*myqchi(j+1,m)
              do j=2,MNUM-1
                an=j
                myqu(j,m) = ER_INV*dcmplx(0.0d0,1.0d0)*am*myqchi(j,m) &
                 &        + ER_INV*0.5d0*an*( myqpsi(j-1,m) - myqpsi(j+1,m) )
                myqv(j,m) = ER_INV*dcmplx(0.0d0,1.0d0)*am*myqpsi(j,m) &
                 &        + ER_INV*0.5d0*an*( -myqchi(j-1,m) + myqchi(j+1,m) )
              end do
              j=MNUM
              an=j 
              myqu(j,m) = ER_INV*dcmplx(0.0d0,1.0d0)*am*myqchi(j,m)              &
               &        + ER_INV*0.5d0*an*myqpsi(j-1,m)
              myqv(j,m) = ER_INV*dcmplx(0.0d0,1.0d0)*am*myqpsi(j,m)             &
               &        - ER_INV*0.5d0*an*myqchi(j-1,m)
            end if
          end do
  !$OMP END DO
  !$OMP END PARALLEL
        !
