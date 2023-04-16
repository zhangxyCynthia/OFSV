MODULE ElementType 
    IMPLICIT REAL*8(A-H,O-Z)
   
    TYPE D_p 
        REAL(kind=8) ::node(2) !node
    END TYPE D_p 
    
    
    TYPE E_p  
         REAL(kind=8) :: initial(9)
         REAL(kind=8) :: du0(9),L1(9) ! initial coefficients of every time step
         REAL(kind=8) :: du1(9),L2(9) ! coefficients updated by first intermediate stage of rk3
         REAL(kind=8) :: du2(9),L3(9) ! coefficients updated by second intermediate stage of rk3
         REAL(kind=8) :: du3(9),L4(9) ! coefficients updated by second intermediate stage of rk3
         REAL(kind=8) :: fflux_x(9),fflux_y(9) ! element boundary integral
         REAL(kind=8) :: iflux(9) ! element integral
         REAL(kind=8) :: x(2) ! element center
         REAL(kind=8) :: distx(2) ! element size
         REAL(kind=8) :: subnode(4,4,2) ! subnode
         REAL(kind=8) :: xg(3,3,2) ! subelement center
         REAL(kind=8) :: distxg(3,3,2) ! subelement size
         REAL(kind=8) :: correction(9)  ! sharp correction
         REAL(kind=8) :: dtu_x(9),dtu_y(9)  ! sharp correction
         
         
    END TYPE E_p 



    TYPE F_p 
         REAL(kind=8) :: hflux(2) ! numerical flux of boundary
         REAL(kind=8) :: hgflux(3,3), hlgflux(4,3,3)  ! numerical flux of boundary at gauss points
         REAL(kind=8) :: partial_jump(6,2)
         
    END TYPE F_p
    

    END MODULE ElementType
!-------------------------------------------------------------------	 

!-------------------------------------------------------------------
    MODULE GlobalData     
    USE ElementType  
    IMPLICIT REAL*8(A-H,O-Z) 
    integer :: IT,NP1,NP2,stage
    real(kind=8) :: PI,CK,GAM        
    real(kind=8) :: tau,time,tstop,utime,sstop               
    real(kind=8) :: ep,dx,dy,stiff(9,9)

    

    type(E_p),allocatable,dimension(:,:)::pElem 
    type(F_p),allocatable,dimension(:,:)::pFace_x,pFace_y
    type(D_p),allocatable,dimension(:,:)::peNode
    END MODULE GlobalData
!-------------------------------------------------------------------

