module LU

implicit none

contains

subroutine linear(AA,bb,xx)
        implicit none
        integer, parameter :: dp = selected_real_kind(15, 307)
        integer :: ii,jj,kk,nn
        real(dp), dimension(:,:),intent(in) :: AA
        real(dp), dimension(:),intent(in) :: bb                !size(A,1)
        real(dp), dimension(:),intent(inout) :: xx             !size(A,1) and inout instead of out
        real(dp), dimension(size(AA,1),size(AA,2)) :: LL
        real(dp), dimension(size(AA,1),size(AA,2)) :: UU
        real(dp), dimension(size(AA,1)) :: dd
        real(dp) :: sum

        nn=size(AA,1)
        dd=0
        LL=0
        UU=0

        !---LU Doolitle---

        do ii=1,nn
            do kk=ii,nn
                sum=0
                    do jj=1,ii
                        sum=sum+(LL(ii,jj)*UU(jj,kk)) !U filling
                    enddo
                UU(ii,kk)=AA(ii,kk)-sum
            enddo

            do kk=ii,nn
                if (ii==kk) then
                    LL(ii,kk)=1
                else
                    sum=0
                        do jj=1,ii
                            sum=sum+(LL(kk,jj)*UU(jj,ii)) !L filling
                        enddo
                    LL(kk,ii)=(AA(kk,ii)-sum)/UU(ii,ii)
                endif
            enddo
        enddo
        
        !---Solving LUx=b

        do ii=1,nn
            sum=0
                do jj=1,ii
                    sum=sum+(LL(ii,jj)*dd(jj))
                enddo
            dd(ii)=(bb(ii)-sum)/LL(ii,ii)
        enddo

        xx(nn)=dd(nn)/UU(nn,nn)

        do ii=nn-1,1,-1
            sum=0
                do jj=ii+1,nn
                    sum=sum+(UU(ii,jj)*xx(jj))            
                enddo
            xx(ii)=(dd(ii)-sum)/UU(ii,ii)
        enddo

end subroutine linear

end module LU