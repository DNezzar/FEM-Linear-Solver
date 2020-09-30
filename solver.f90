module mod

    implicit none

 contains

    subroutine solver(conn,node,free,A,E,f,d,strain,stress)

    integer, parameter :: dp = selected_real_kind(15, 307)    ! 15 significant digits of precision and an exponent range of at least 307
    integer(dp) :: n,m,ndof,i                                 ! Number of nodes, Number of elements, Number of dof, index
    integer(dp) :: n1,n2                                      ! 2 Nodes of an element
    integer(dp), dimension(1:4) :: scat                       ! Vector location where the matrix K is to scatter
    integer(dp), dimension(:), intent(in) :: free             ! Vector of free dof
    real(dp), dimension(:,:), intent(in) :: node              ! Node: Matrix of node coordinates (x,y)
    integer(dp), dimension(:,:), intent(in) :: conn              ! Conn: Matrix of element connectivity (1st node, 2nd node)
    real(dp), dimension(:,:), allocatable :: K, Kfree         ! Global Stiffness Matrix, Sub Stiffness Matrix of free dof 
    real(dp), dimension(size(node,1)*2), intent(out) :: d     ! Displacment vector, Sub Displacment Vector of free dof
    real(dp), dimension(size(conn,1)), intent(out) :: stress
    real(dp), dimension(size(conn,1)), intent(out) :: strain
    real(dp), dimension(:), allocatable :: dfree
    real(dp), dimension(:), intent(inout) :: f                ! Force vector, Sub Force Vector of free dof
    real(dp), dimension(:), allocatable :: ffree              ! Force vector, Sub Force Vector of free dof
    real(dp), dimension(1:4) :: v,B                           ! Vector of cosinus
    real(dp), dimension(1:4,1:4) :: T,ke                      ! Transformation Matrix, Element Stiffness Matrix
    real(dp) :: L,x1,x2,y1,y2,c,s,c2,s2,cs                    ! A: Area, E: Young's Modulus, P: Force Magnitude, Length of an element,cosinus,sinus (and square), strain and stress of each element
    real(dp),intent(in) :: A,E

    
    n=size(node,1)
    m=size(conn,1)
    ndof=2*n
    allocate(K(1:ndof,1:ndof))                                !Allocation of size global Stiffness Matrix

    !---initialization
    K=0.0
    d=0.0

    allocate(Kfree(1:size(free),1:size(free)))                !Allocation of free dof reduced stiffness matrix
    allocate(dfree(1:size(free)))
    allocate(ffree(1:size(free)))
    Kfree=0.0
    dfree=0.0
    ffree=0.0

!----------------------END INPUT-----------------------
!---------Filling Global Stiffness Matrix
    do i=1,m                                         
        n1=conn(i,1) !ID of the first node
        n2=conn(i,2) !ID of the second node

        x1=node(n1,1)
        y1=node(n1,2)

        x2=node(n2,1)
        y2=node(n2,2)
        
        L=sqrt(((x2-x1)**2)+((y2-y1)**2))

        c=(x2-x1)/L
        s=(y2-y1)/L
        c2=c*c
        s2=s*s
        cs=c*s

        T(1,:)=(/c2,cs,-c2,-cs/)
        T(2,:)=(/cs,s2,-cs,-s2/)
        T(3,:)=(/-c2,-cs,c2,cs/)
        T(4,:)=(/-cs,-s2,cs,s2/)

        ke=(A*E/L)*T
        
        scat=(/2*n1-1,2*n1,2*n2-1,2*n2/)

        K(scat,scat)=K(scat,scat)+ke    !add ke into K
    enddo 
    Kfree=K(free,free)
    ffree=f(free)
    dfree=d(free)

!---------Solving Linear Systeme

    call linear(Kfree,ffree,dfree)
    
!---------Compute Nodale Displacment

   d(free)=dfree
   
!---------Compute Reaction Forces
    
    f=matmul(K,d)

!---------Compute Strain and Stress

   do i=1,m
        n1=conn(i,1) !ID of the first node
        n2=conn(i,2) !ID of the second node

        x1=node(n1,1)
        y1=node(n1,2)

        x2=node(n2,1)
        y2=node(n2,2)
        
        L=sqrt(((x2-x1)**2)+((y2-y1)**2))

        c=(x2-x1)/L
        s=(y2-y1)/L
        v=(/-c,-s, c, s/)
        B=(1/L)*v
        scat=(/2*n1-1,2*n1,2*n2-1,2*n2/)

        strain(i)=dot_product(B,d(scat))
        stress(i)=E*strain(i)
   enddo




!-----------------PRINT-----------

    print*, "N o d a l  D i s p l a c m e n t s"
    print*, "NID        X-DISP        Y-DISP"

    do i=1,n
        print "(i3,es17.3,es14.3)",i,d(2*i-1),d(2*i)
    enddo

    print*, "E x t e r n a l  F o r c e s"
    print*, "NID        X-FORCE       Y-FORCE"

    do i=1,n
        print "(i3,es17.3,es14.3)",i,f(2*i-1),f(2*i)
    enddo

    print*, "E l e m e n t   S t r e s s"
    print*, "EID         STRAIN       STRESS"

    do i=1,m
        print "(i3,es17.3,es14.3)",i,strain(i),stress(i)
    enddo

!---------------END PRINT


end subroutine solver

subroutine linear(AA,bb,xx)
        
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

end module mod