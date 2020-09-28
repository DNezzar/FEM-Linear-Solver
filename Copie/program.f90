program test

    use mod 

    implicit none

!------------------------INPUT-------------------------
    integer, parameter :: dp = selected_real_kind(15, 307)    ! 15 significant digits of precision and an exponent range of at least 307
    integer(dp) :: n,m,ndof,i                                 ! Number of nodes, Number of elements, Number of dof, index
    integer(dp) :: n1,n2                                      ! 2 Nodes of an element
    integer(dp), dimension(1:4) :: scat                       ! Vector location where the matrix K is to scatter
    integer(dp), dimension(:), allocatable :: free            ! Vector of free dof
    real(dp), dimension(:,:), allocatable :: Node             ! Node: Matrix of node coordinates (x,y)
    real(dp), dimension(:,:), allocatable :: Conn             ! Conn: Matrix of element connectivity (1st node, 2nd node)
    real(dp), dimension(:,:), allocatable :: K, Kfree         ! Global Stiffness Matrix, Sub Stiffness Matrix of free dof 
    real(dp), dimension(:), allocatable :: d,dfree            ! Displacment vector, Sub Displacment Vector of free dof
    real(dp), dimension(:), allocatable :: f,ffree            ! Force vector, Sub Force Vector of free dof
    real(dp), dimension(1:4) :: v,B                           ! Vector of cosinus
    real(dp), dimension(1:4,1:4) :: T,ke                      ! Transformation Matrix, Element Stiffness Matrix
    real(dp) :: A,E,P,L,x1,x2,y1,y2,c,s,c2,s2,cs,strain,stress! A: Area, E: Young's Modulus, P: Force Magnitude, Length of an element,cosinus,sinus (and square), strain and stress of each element

    print *, "Combien de noeuds ?"
    read *, n
    print *, "Combien d'elements ?"
    read *, m

    allocate(Node(1:n,1:2))                                   !Allocation of size of matrix of node coordinates in 2D
    allocate(Conn(1:m,1:2))                                   !Allocation of size of matrix of connectivity in 2D

    ndof=2*n
    allocate(K(1:ndof,1:ndof))                                !Allocation of size global Stiffness Matrix
    allocate(d(1:ndof))                                       !Allocation of size global displacment vector
    allocate(f(1:ndof))                                       !Allocation of size global force vector
     

    !---initialization
    K=0.0
    d=0.0
    f=0.0
    !-----Properties
    A=1.0
    E=10000000.0 

    !------BC------
    P=1000.0
    f(6)=P
    allocate(free(1:4))                                       !Input of the free dof
    free=(/3,4,5,6/)                                          !Taille max=ndof

    allocate(Kfree(1:size(free),1:size(free)))                !Allocation of free dof reduced stiffness matrix
    allocate(dfree(1:size(free)))
    allocate(ffree(1:size(free)))
    Kfree=0.0
    dfree=0.0
    ffree=0.0

    do i=1,n
        print "(a35,i1)", "entrez coordonnees (x,y) du point ",i
        read *, Node(i,1)
        read *, Node(i,2)
    enddo
    Node=Node*12
    do i=1,m
        print "(a34,i1)", "entrez les 2 noeuds de l'element ",i
        read*, Conn(i,1)
        read*, Conn(i,2)
    enddo
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
    d(free)=dfree
!---------Compute Nodale Displacment

    print*, "N o d a l  D i s p l a c m e n t s"
    print*, "NID        X-DISP        Y-DISP"

    do i=1,n
        print "(i3,es17.3,es14.3)",i,d(2*i-1),d(2*i)
    enddo
!---------Compute Reaction Forces
    f=matmul(K,d)

    print*, "E x t e r n a l  F o r c e s"
    print*, "NID        X-FORCE       Y-FORCE"

    do i=1,n
        print "(i3,es17.3,es14.3)",i,f(2*i-1),f(2*i)
    enddo
!---------Compute Strain and Stress
    
    print*, "E l e m e n t   S t r e s s"
    print*, "EID         STRAIN       STRESS"

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

        strain=dot_product(B,d(scat))
        stress=E*strain

        print "(i3,es17.3,es14.3)",i,strain,stress
   enddo

end program test