module mod
    implicit none

 contains

    subroutine linear(A,b,x)
        implicit none
        integer, parameter :: dp = selected_real_kind(15, 307)
        integer :: i,j,k,n
        real(dp), dimension(:,:),intent(in) :: A
        real(dp), dimension(:),intent(in) :: b                !size(A,1)
        real(dp), dimension(:),intent(inout) :: x             !size(A,1) and inout instead of out
        real(dp), dimension(size(A,1),size(A,2)) :: L
        real(dp), dimension(size(A,1),size(A,2)) :: U
        real(dp), dimension(size(A,1)) :: d
        real(dp) :: sum

        n=size(A,1)
        d=0
        x=0
        L=0
        U=0

        

        !---LU Doolitle---

        do i=1,n
            do k=i,n
                sum=0
                    do j=1,i
                        sum=sum+(L(i,j)*U(j,k)) !U filling
                    enddo
                U(i,k)=A(i,k)-sum
            enddo

            do k=i,n
                if (i==k) then
                    L(i,k)=1
                else
                    sum=0
                        do j=1,i
                            sum=sum+(L(k,j)*U(j,i)) !L filling
                        enddo
                    L(k,i)=(A(k,i)-sum)/U(i,i)
                endif
            enddo
        enddo
        
        !---Solving LUx=b

        do i=1,n
            sum=0
                do j=1,i
                    sum=sum+(L(i,j)*d(j))
                enddo
            d(i)=(b(i)-sum)/L(i,i)
        enddo

        x(n)=d(n)/U(n,n)

        do i=n-1,1,-1
            sum=0
                do j=i+1,n
                    sum=sum+(U(i,j)*x(j))            
                enddo
            x(i)=(d(i)-sum)/U(i,i)
        enddo

    end subroutine linear

end module mod