program test
    use mod

    implicit none
    integer, parameter :: dp = selected_real_kind(15, 307)    ! 15 significant digits of precision and an exponent range of at least 307
    integer(dp) :: i
    real(dp), dimension(1:4,1:2) :: coord
    real(dp), dimension(1:6,1:2) :: con
    real(dp), dimension(1:8) :: force
    real(dp), dimension(:), allocatable :: disp
    integer(dp), dimension(1:4) :: libre
    real(dp), dimension(:), allocatable :: stra, stre
    real(dp) :: S,Young 

    coord(1,:)=(/0,0/)
    coord(2,:)=(/20,0/)
    coord(3,:)=(/20,15/)
    coord(4,:)=(/0,15/)

    coord=coord*12

    con(1,:)=(/1,2/)
    con(2,:)=(/1,3/)
    con(3,:)=(/2,4/)
    con(4,:)=(/2,3/)
    con(5,:)=(/4,3/)
    con(6,:)=(/1,4/)

    libre=(/3,4,5,6/)
    force=(/0,0,0,0,0,1000,0,0/)
    S=1.0
    Young=10000000

    call solver(con,coord,libre,S,Young,force,disp,stra,stre)

    print*, "N o d a l  D i s p l a c m e n t s"
    print*, "NID        X-DISP        Y-DISP"

    do i=1,size(coord,1)
        print "(i3,es17.3,es14.3)",i,disp(2*i-1),disp(2*i)
    enddo

    print*, "E x t e r n a l  F o r c e s"
    print*, "NID        X-FORCE       Y-FORCE"

    do i=1,size(coord,1)
        print "(i3,es17.3,es14.3)",i,force(2*i-1),force(2*i)
    enddo

    print*, "E l e m e n t   S t r e s s"
    print*, "EID         STRAIN       STRESS"

    do i=1,size(con,1)
        print "(i3,es17.3,es14.3)",i,strain(i),stress(i)
    enddo
    
    


end program test