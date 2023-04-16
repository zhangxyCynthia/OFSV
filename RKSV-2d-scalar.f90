! 'ElementType' manages data structures of every element
!-----------------------------------------------------------------
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



!------------------------------------------------------------------- 
      SUBROUTINE GETTIMESTEP 
      USE ElementType      
      USE GlobalData    
      IMPLICIT REAL*8(A-H,O-Z)  
      

 
      
      TAU=1.25D-3*DBLE(32.d0/(NP1-2))
      

      UBD=0.0  
      END SUBROUTINE
!-------------------------------------------------------------------
 
 
!periodic boundary condition
!-------------------------------------------------------------------
      SUBROUTINE BOUNDARY
      USE ElementType      
      USE GlobalData    
      IMPLICIT NONE 
      INTEGER :: i,j
      
      SELECT CASE(stage)
         CASE(1)
         DO J=2,NP2-1
         pElem(1,J)%du0(:)=pElem(NP1-1,J)%du0(:)
         pElem(NP1,J)%du0(:)=pElem(2,J)%du0(:)
         END DO
      
         DO I=2,NP1-1
         pElem(I,1)%du0(:)=pElem(I,NP2-1)%du0(:)
         pElem(I,NP2)%du0(:)=pElem(I,2)%du0(:)
         END DO
         
         CASE(2)
         DO J=2,NP2-1
         pElem(1,J)%du1(:)=pElem(NP1-1,J)%du1(:)
         pElem(NP1,J)%du1(:)=pElem(2,J)%du1(:)
         END DO
      
         DO I=2,NP1-1
         pElem(I,1)%du1(:)=pElem(I,NP2-1)%du1(:)
         pElem(I,NP2)%du1(:)=pElem(I,2)%du1(:)
         END DO
         
         CASE(3)
         DO J=2,NP2-1
         pElem(1,J)%du2(:)=pElem(NP1-1,J)%du2(:)
         pElem(NP1,J)%du2(:)=pElem(2,J)%du2(:)
         END DO
         
         DO I=2,NP1-1   
         pElem(I,1)%du2(:)=pElem(I,NP2-1)%du2(:)
         pElem(I,NP2)%du2(:)=pElem(I,2)%du2(:)
         END DO
         
         CASE(4)
         DO J=2,NP2-1
         pElem(1,J)%du3(:)=pElem(NP1-1,J)%du3(:)
         pElem(NP1,J)%du3(:)=pElem(2,J)%du3(:)
         END DO
         
         DO I=2,NP1-1   
         pElem(I,1)%du3(:)=pElem(I,NP2-1)%du3(:)
         pElem(I,NP2)%du3(:)=pElem(I,2)%du3(:)
         END DO
      
         END SELECT
          
      !  SELECT CASE(stage)
      !     CASE(1)
      !     DO J=2,NP2-1
      !     pElem(1,J)%du0(:,:)=pElem(2,J)%du0(:,:)
      !     pElem(NP1,J)%du0(:,:)=pElem(NP1-1,J)%du0(:,:)
      !     END DO
      !
      !     DO I=2,NP1-1
      !     pElem(I,1)%du0(:,:)=pElem(I,2)%du0(:,:)
      !     pElem(I,NP2)%du0(:,:)=pElem(I,NP2-1)%du0(:,:)
      !     END DO
      !    
      !     CASE(2)
      !     DO J=2,NP2-1
      !     pElem(1,J)%du1(:,:)=pElem(2,J)%du1(:,:)
      !     pElem(NP1,J)%du1(:,:)=pElem(NP1-1,J)%du1(:,:)
      !     END DO
      !
      !     DO I=2,NP1-1
      !     pElem(I,1)%du1(:,:)=pElem(I,2)%du1(:,:)
      !     pElem(I,NP2)%du1(:,:)=pElem(I,NP2-1)%du1(:,:)
      !     END DO
      !    
      !     CASE(3)
      !     DO J=2,NP2-1
      !     pElem(1,J)%du2(:,:)=pElem(2,J)%du2(:,:)
      !     pElem(NP1,J)%du2(:,:)=pElem(NP1-1,J)%du2(:,:)
      !     END DO
      !    
      !     DO I=2,NP1-1   
      !     pElem(I,1)%du2(:,:)=pElem(I,2)%du2(:,:)
      !     pElem(I,NP2)%du2(:,:)=pElem(I,NP2-1)%du2(:,:)
      !     END DO
      !    
      !END SELECT
      
    END SUBROUTINE
!-------------------------------------------------------------------

!-------------------------------------------------------------------
    SUBROUTINE GAUSS(xc,dist,pointg)
    USE ElementType      
    USE GlobalData   
    IMPLICIT REAL*8(A-H,O-Z) 
      real(kind=8) :: x,y,w
      real(kind=8) :: pointg(2,2),xc(2),dist(2)
      integer ::i,j
      
      w=DSQRT(3.d0)/3.d0
    
    pointg(1,1)=xc(1)-w*dist(1)/2.d0
    pointg(1,2)=xc(1)+w*dist(1)/2.d0
    
    
    pointg(2,1)=xc(2)-w*dist(2)/2.d0
    pointg(2,2)=xc(2)+w*dist(2)/2.d0
    
    END SUBROUTINE
!-------------------------------------------------------------------
!-------------------------------------------------------------------
    SUBROUTINE GAUSS3(xc,dist,pointg,weight)
    USE ElementType      
    USE GlobalData   
    IMPLICIT REAL*8(A-H,O-Z) 
      real(kind=8) :: x,y,g
      real(kind=8) :: pointg(2,3),xc(2),dist(2),weight(3)
      integer ::i,j
      
      g=DSQRT(3.d0/5.d0)
    
    pointg(1,1)=xc(1)-g*dist(1)/2.d0
    pointg(1,2)=xc(1)
    pointg(1,3)=xc(1)+g*dist(1)/2.d0
    
    
    pointg(2,1)=xc(2)-g*dist(2)/2.d0
    pointg(2,2)=xc(2)
    pointg(2,3)=xc(2)+g*dist(2)/2.d0
    
    weight(1)=5.d0/9.d0
    weight(2)=8.d0/9.d0
    weight(3)=5.d0/9.d0
    
    END SUBROUTINE
!-------------------------------------------------------------------
    
!-------------------------------------------------------------------
    SUBROUTINE BASISELEM(xc,xg,elem) 
    USE ElementType      
    USE GlobalData   
    IMPLICIT NONE 
    real(kind=8) :: xc(2),xg(2)
    real(kind=8) :: elem(9) 

    elem(1)=1.d0
    
    elem(2)=(xg(1)-xc(1))/(DX/2.d0)
    elem(3)=(xg(2)-xc(2))/(DY/2.d0)
    elem(4)=elem(2)*elem(3)
    
    elem(5)=3.d0/2.d0*elem(2)**2-1.d0/2.d0
    elem(6)=3.d0/2.d0*elem(3)**2-1.d0/2.d0
    
    elem(7)=elem(5)*elem(3)
    elem(8)=elem(2)*elem(6)
    elem(9)=elem(5)*elem(6)
    
    END SUBROUTINE
!-------------------------------------------------------------------
    
 SUBROUTINE BASISELEM_de(xc,xg,elemd) 
    USE ElementType      
    USE GlobalData   
    IMPLICIT NONE
    real(kind=8) :: xc(2),xg(2)
    real(kind=8) :: elemd(2,9)
  
    !\partial x
    elemd(1,1)=0.d0
    
    elemd(1,2)=2.d0/DX
    elemd(1,3)=0.d0
    elemd(1,4)=elemd(1,2)*(xg(2)-xc(2))/(DY/2.d0)
    
    elemd(1,5)=3.d0*(xg(1)-xc(1))/(DX/2.d0)*(2.d0/DX)
    elemd(1,6)=0.d0
    
    elemd(1,7)=elemd(1,5) * (xg(2)-xc(2))/(DY/2.d0)
    elemd(1,8)=elemd(1,2) * (3.d0/2.d0*((xg(2)-xc(2))/(DY/2.d0))**2-1.d0/2.d0)
    elemd(1,9)=elemd(1,5) * (3.d0/2.d0*((xg(2)-xc(2))/(DY/2.d0))**2-1.d0/2.d0)
    
    !\partial y
    elemd(2,1)=0.d0
    
    elemd(2,2)=0.d0
    elemd(2,3)=2.d0/DY
    elemd(2,4)=elemd(2,3)*((xg(1)-xc(1))/(DX/2.d0))
    
    elemd(2,5)=0.d0
    elemd(2,6)=3.d0*2.d0/DY*(xg(2)-xc(2))/(DY/2.d0)
    
    elemd(2,7)=elemd(2,3)*(3.d0/2.d0*((xg(1)-xc(1))/(DX/2.d0))**2-1.d0/2.d0)
    elemd(2,8)=elemd(2,6)*((xg(1)-xc(1))/(DX/2.d0))
    elemd(2,9)=elemd(2,6)*(3.d0/2.d0*((xg(1)-xc(1))/(DX/2.d0))**2-1.d0/2.d0)

    END SUBROUTINE

!-------------------------------------------------------------------
SUBROUTINE BASISELEM_de2(xc,xg,elemd2) 
    USE ElementType      
    USE GlobalData   
    IMPLICIT NONE
    real(kind=8) :: xc(2),xg(2)
    real(kind=8) :: elemd2(3,9)
  
    
     !\partial x^2
    elemd2(1,1)=0.d0
    
    elemd2(1,2)=0.d0
    elemd2(1,3)=0.d0
    elemd2(1,4)=0.d0
    
    elemd2(1,5)=3.d0*2.d0/DX*(2.d0/DX)
    elemd2(1,6)=0.d0
    
    elemd2(1,7)=elemd2(1,5) * ((xg(2)-xc(2))/(DY/2.d0))
    elemd2(1,8)=0.d0
    elemd2(1,9)=elemd2(1,5) * (3.d0/2.d0*((xg(2)-xc(2))/(DY/2.d0))**2-1.d0/2.d0)
    
    !\partial xy
    elemd2(2,1)=0.d0
    
    elemd2(2,2)=0.d0
    elemd2(2,3)=0.d0
    elemd2(2,4)=2.d0/DX*(2.d0/DY)
    
    elemd2(2,5)=0.d0
    elemd2(2,6)=0.d0
    
    elemd2(2,7)=3.d0*(2.d0/DX)*(xg(1)-xc(1))/(DX/2.d0)*(2.d0/DY)
    elemd2(2,8)=2.d0/DX * (3.d0*2.d0/DY*(xg(2)-xc(2))/(DY/2.d0))
    elemd2(2,9)=3.d0*2.d0/DX*(xg(1)-xc(1))/(DX/2.d0) * (3.d0*2.d0/DY*(xg(2)-xc(2))/(DY/2.d0))

    !\partial y^2
    elemd2(3,1)=0.d0
    
    elemd2(3,2)=0.d0
    elemd2(3,3)=0.d0
    elemd2(3,4)=0.d0
    
    elemd2(3,5)=0.d0
    elemd2(3,6)=3.d0*2.d0/DY*(2.d0/DY)
    
    elemd2(3,7)=0.d0
    elemd2(3,8)=(xg(1)-xc(1))/(DX/2.d0) * elemd2(3,6)
    elemd2(3,9)=(3.d0/2.d0*((xg(1)-xc(1))/(DX/2.d0))**2-1.d0/2.d0) * elemd2(3,6)


    END SUBROUTINE
!-------------------------------------------------------------------   
!-------------------------------------------------------------------


!-------------------------------------------------------------------   
      SUBROUTINE INITIAL 
      USE ElementType      
      USE GlobalData   
      IMPLICIT REAL*8(A-H,O-Z) 
      real(kind=8) :: con
      real(kind=8) :: elem(9),wei(3)
      real(kind=8) :: pointg(2,2),pointgi(2,3)
      real(kind=8) :: Mass(9),xc(2),xg(2)
      real(kind=8) :: xx,yy
      real(kind=8) :: vel0,vey0,p0,rho0
      integer ::i,j,n,ii,jj,a
      
      ! for 32+2
      ! eg: 32 is the quantity of the elements of computational domain
      !      2 is the quantity of boundary
      NP1=8+2 ! the quantity of elements of x direction
	  NP2=8+2 ! the quantity of elements of y direction
 
      PI=4.d0*atan(1.d0)
      CK=3.d0 
      GAM=1.4d0
      
      
      
	  TSTOP=1.9999999999d0 ! the whole time of computation
      !TSTOP=0.19999999999d0 ! the whole time of computation for sod
      !TSTOP=0.13999999999d0 ! the whole time of computation for lax
      !TSTOP=1.79999999999d0 ! the whole time of computation for shu
      !TSTOP=4.99999999999d0 ! the whole time of computation for Titarev 
      UTIME=TSTOP/5.d0 ! decide the quantity of pictures
	  SSTOP=0.d0 
      TIME=0.d0 
      IT=0

      ! the computation domain is [0,2]x[0,2]
      DX = 2.d0/DBLE(NP1-2) 
      DY = 2.d0/DBLE(NP2-2)
      !DX = 1.d0/DBLE(NP1-2) 
      !DY = 1.d0/DBLE(NP2-2)
      !DX = 10.d0/DBLE(NP1-2) 
      !DY = 10.d0/DBLE(NP2-2)
 
      ep=0.001

      ! allocate memory
      ALLOCATE(pElem(NP1,NP2))    
      ALLOCATE(pFace_x(NP1+1,NP2))
      ALLOCATE(pFace_y(NP1,NP2+1))   	
      ALLOCATE(peNode(NP1+1,NP2+1))
      
      
      
      ! stiff matrix
 
         
      stiff(1,:)=(/  0.25d0,                      0.25d0, 0.25d0, 0.25d0,                       0.25d0, 0.25d0, 0.25d0,                      0.25d0, 0.25d0 /)
      
      stiff(2,:)=(/ -0.75d0,                        0.d0, 0.75d0,-0.75d0,                         0.d0, 0.75d0,-0.75d0,                        0.d0, 0.75d0 /)
      stiff(3,:)=(/ -0.75d0,                     -0.75d0, -0.75d0, 0.d0, 0.d0, 0.d0,  0.75d0,                       0.75d0, 0.75d0 /)
      stiff(4,:)=(/  2.25d0,                        0.d0, -2.25d0, 0.d0, 0.d0, 0.d0, -2.25d0,                         0.d0, 2.25d0 /)
      
      stiff(5,:)=(/  0.75d0,3.d0/4.d0*(1.d0-DSQRT(3.d0)), 0.75d0, 0.75d0, 3.d0/4.d0*(1.d0-DSQRT(3.d0)), 0.75d0, 0.75d0,3.d0/4.d0*(1.d0-DSQRT(3.d0)), 0.75d0 /)
      stiff(6,:)=(/  0.75d0,                       0.75d0, 0.75d0,  3.d0/4.d0*(1.d0-DSQRT(3.d0)),3.d0/4.d0*(1.d0-DSQRT(3.d0)), 3.d0/4.d0*(1.d0-DSQRT(3.d0)),  0.75d0,                        0.75d0, 0.75d0 /)
      
      stiff(7,:)=(/ -2.25d0,9.d0/4.d0*(DSQRT(3.d0)-1.d0), -2.25d0, 0.d0, 0.d0, 0.d0,  2.25d0,-9.d0/4.d0*(DSQRT(3.d0)-1.d0), 2.25d0 /)
      stiff(8,:)=(/ -2.25d0,                         0.d0, 2.25d0,  9.d0/4.d0*(DSQRT(3.d0)-1.d0),                        0.d0,-9.d0/4.d0*(DSQRT(3.d0)-1.d0), -2.25d0,                          0.d0, 2.25d0 /)
      stiff(9,:)=(/  2.25d0,-9.d0/4.d0*(DSQRT(3.d0)-1.d0), 2.25d0, -9.d0/4.d0*(DSQRT(3.d0)-1.d0),9.d0/2.d0*(2.d0-DSQRT(3.d0)),-9.d0/4.d0*(DSQRT(3.d0)-1.d0),  2.25d0, -9.d0/4.d0*(DSQRT(3.d0)-1.d0), 2.25d0 /)
 
      
      !write(*,*) stiff
      !! the mass matrix of basis
      Mass(1)=DX*DY
      
      Mass(2)=1.d0/3.d0*DX*DY
      Mass(3)=1.d0/3.d0*DX*DY
      Mass(4)=1.d0/9.d0*DX*DY
      
      Mass(5)=1.d0/5.d0*DX*DY
      Mass(6)=1.d0/5.d0*DX*DY
      
      Mass(7)=1.d0/15.d0*DX*DY
      Mass(8)=1.d0/15.d0*DX*DY
      Mass(9)=1.d0/25.d0*DX*DY



      DO IY=1,NP2+1
      DO IX=1,NP1+1
          xx=DBLE(IX-2)*DX
          yy=DBLE(IY-2)*DY    
          peNode(IX,IY)%node(1) = xx
          peNode(IX,IY)%node(2) = yy
      ENDDO
      ENDDO
      
      !DO IY=1,NP2+1
      !DO IX=1,NP1+1
      !    xx=DBLE(IX-2)*DX
      !    yy=DBLE(IY-2)*DY    
      !    peNode(IX,IY)%node(1) = xx-5.d0
      !    peNode(IX,IY)%node(2) = yy-5.d0
      !ENDDO
      !ENDDO
      

      DO IY=1,NP2
      DO IX=1,NP1
          pElem(IX,IY)%x(:)=(peNode(IX,IY)%node(:) + peNode(IX+1,IY)%node(:) + peNode(IX,IY+1)%node(:) + peNode(IX+1,IY+1)%node(:))/4.d0

          pElem(IX,IY)%distx(1)=peNode(IX+1,IY)%node(1)- peNode(IX,IY)%node(1)
          pElem(IX,IY)%distx(2)=peNode(IX,IY+1)%node(2)- peNode(IX,IY)%node(2)
          pElem(IX,IY)%correction(:)=0.d0
          
      ENDDO
      ENDDO
      
      DO IY=1,NP2
      DO IX=1,NP1
      
            CALL GAUSS(pElem(IX,IY)%x(:),pElem(IX,IY)%distx(:),pointg)
            
             pElem(IX,IY)%subnode(1,:,1)=peNode(IX,IY)%node(1)
             pElem(IX,IY)%subnode(2,:,1)=pointg(1,1)
             pElem(IX,IY)%subnode(3,:,1)=pointg(1,2)
             pElem(IX,IY)%subnode(4,:,1)=peNode(IX+1,IY)%node(1)
             
             pElem(IX,IY)%subnode(:,1,2)=peNode(IX,IY)%node(2)
             pElem(IX,IY)%subnode(:,2,2)=pointg(2,1)
             pElem(IX,IY)%subnode(:,3,2)=pointg(2,2)
             pElem(IX,IY)%subnode(:,4,2)=peNode(IX,IY+1)%node(2)
             
             DO ii=1,3
             DO jj=1,3
                pElem(IX,IY)%xg(ii,jj,:)=(pElem(IX,IY)%subnode(ii,jj,:)+pElem(IX,IY)%subnode(ii+1,jj,:)+&
                                          pElem(IX,IY)%subnode(ii,jj+1,:)+pElem(IX,IY)%subnode(ii+1,jj+1,:))/4.d0
                
                pElem(IX,IY)%distxg(ii,jj,1)=pElem(IX,IY)%subnode(ii+1,jj,1)- pElem(IX,IY)%subnode(ii,jj,1)
                pElem(IX,IY)%distxg(ii,jj,2)=pElem(IX,IY)%subnode(ii,jj+1,2)- pElem(IX,IY)%subnode(ii,jj,2)
             ENDDO
             ENDDO
             
      ENDDO
      ENDDO

      !vel0=1.d0
      !vey0=0.d0 ! initial velocity of y direction
      !p0=1.d0 ! initial pressure
      

      vel0=1.d0
      vey0=1.d0 ! initial velocity of y direction
      p0=1.d0 ! initial pressure
      
      ! compute the initial coefficients
      DO IY=1,NP2
      DO IX=1,NP1  
       
      ! compute the Gauss points of element
          
          xc(:)= pElem(IX,IY)%x(:)
      
      CALL GAUSS3(xc,pElem(IX,IY)%distx(:),pointgi,wei)
      pElem(IX,IY)%du0(:)=0.d0
      
      
      ! compute initial conservative variables on Gauss points
      DO i=1,3
      DO j=1,3
      xg(1)=pointgi(1,i)
      xg(2)=pointgi(2,j)
      
      CALL BASISELEM(xc,xg,elem)    


      con= 1.d0+0.2d0*DSIN(PI*(xg(1)+xg(2))) ! density
      
      ! compute six coefficients(degrees of freedom) of conservative variables
      DO n=1,9
             pElem(IX,IY)%du0(n)=pElem(IX,IY)%du0(n)+con*wei(i)*wei(j)*elem(n)*DX*DY/4.d0/Mass(n)
      END DO   

      ENDDO
      ENDDO

      pElem(IX,IY)%initial(:)=pElem(IX,IY)%du0(:)

      ENDDO
      ENDDO 
     
    END SUBROUTINE      
!-------------------------------------------------------------------
!---------------------------------------Flux-----------------------------------
!subroutine flux(prim,flx)
!USE ElementType      
!USE GlobalData   
!IMPLICIT REAL*8(A-H,O-Z) 
!
!real*8,dimension(1:4)::prim,flx
!
!flx(1)=prim(1)*prim(2)
!flx(2)=prim(1)*prim(2)**2+prim(4)
!flx(3)=prim(1)*prim(2)*prim(3)
!flx(4)=(prim(1)*(prim(2)**2+prim(3)**2)/2.d0+prim(4)/(GAM-1.d0)+prim(4))*prim(2)
!
!return
!end subroutine
!    
  
!
!subroutine HLLC(statel,stater,flx)
!implicit none
!		
!real*8,dimension(1:4)::stateL,stateR,flx
!
!real*8,dimension(1:4)::pl,pr    !primitive
!
!integer::m
!real*8::al,ar,pvars,pstar,tmp1,tmp2,tmp3,qk, sl,sr,star,gamma
!real*8::flxtmp(1:4),qstar(1:4)
!
!gamma=1.4d0
!
!call exchange(pl,statel,"c2p")   !conservative to primitive
!call exchange(pr,stater,"c2p")
!
!al=sqrt(gamma*pl(4)/pl(1)); ar=sqrt(gamma*pr(4)/pr(1))
!tmp1=0.5d0*(al+ar)
!tmp2=0.5d0*(pl(1)+pr(1))
!
!pvars=0.5d0*(pl(4)+pr(4))-0.5d0*(pr(2)-pl(2))*tmp1*tmp2
!pstar=max(0.d0,pvars)
!
!if(pstar.le.pl(4))then
!  qk=1.d0
!else
!  tmp1=(gamma+1.d0)/(2.d0*gamma); tmp2=(pstar/pl(4)-1.d0)
!  qk=sqrt(1.d0+tmp1*tmp2)
!endif
!sl=pl(2)-al*qk
!
!if(pstar.le.pr(4))then
!  qk=1.d0
!else
!  tmp1=(gamma+1.d0)/(2.d0*gamma); tmp2=(pstar/pr(4)-1.d0)
!  qk=sqrt(1.d0+tmp1*tmp2)
!endif
!sr=pr(2)+ar*qk
!
!tmp1=pr(4)-pl(4)+pl(1)*pl(2)*(sl-pl(2))-pr(1)*pr(2)*(sr-pr(2))
!tmp2=pl(1)*(sl-pl(2))-pr(1)*(sr-pr(2))
!star=tmp1/tmp2
!
!if(sl.ge.0.d0)then
!    call flux(pl,flx)
!elseif(sr.le.0.d0)then
!    call flux(pr,flx)
!elseif((star.ge.0.d0) .and. (sl.le.0.d0))then
!    call flux(pL,flxtmp)
!    call ustarforHLLC(pL(1),pL(2),pL(3),pL(4),sL,star,gamma, qstar)
!    
!    do m=1,4
!    flx(m)=flxtmp(m)+sl*(qstar(m)-statel(m))
!    enddo
!elseif((star.le.0.d0) .and. (sr.ge.0.d0))then
!    call flux(pr,flxtmp)
!    call ustarforHLLC(pr(1),pr(2),pr(3),pr(4),sr,star,gamma, qstar)
!    do m=1,4
!    flx(m)=flxtmp(m)+sr*(qstar(m)-stater(m))
!    enddo
!endif
!
!return	
!end subroutine
!
!
!
!subroutine ustarforHLLC(d1,u1,v1,p1,s1,star1,gamma, ustar)
!implicit none
!
!real*8::d1,u1,v1,p1,s1,star1,gamma, ustar(1:4)
!real*8::tmp1,tmp2,tmp3
!
!tmp1=d1*(s1-u1)/(s1-star1)
!
!tmp2=0.5d0*(u1**2+v1**2)+p1/ ((gamma-1.d0)*d1) 
!tmp3=star1+p1/(d1*(s1-u1))
!
!ustar(1)=tmp1
!ustar(2)=tmp1*star1
!ustar(3)=tmp1*v1
!ustar(4)=tmp1*(tmp2+(star1-u1)*tmp3)
!
!return
!    end subroutine
!!-------------------------------------------------------------------
!    subroutine HLLC1(statel,stater,flx)
!implicit none
!		
!real(kind=8) ::statel,stater,flx
!
!   flx=statel
!   
!return	
!end subroutine
!-------------------------------------------------------------------
!-------------------------------------------------------------------
!-------------------------------------------------------------------
    
    SUBROUTINE GETHFLUX(i,j)
    USE ElementType
    use GlobalData    
    IMPLICIT REAL*8 (A-H,O-Z)
    real(kind=8) :: pointg(2,3),xc(2),xg(2),wei(3)
    real(kind=8) :: corr(2,2) ! the rotation matrix (global --> local)
    real(kind=8) :: statel, stater ! conservative variables(global) on Gauss points of leftside and rightside of element boundary
    real(kind=8) :: statelg, staterg ! local conservative variables
    real(kind=8) :: fluxg ! local flux computed by local conservative variables
    real(kind=8) :: face_l(9),face_r(9) ! the value of basis on Gauss points of leftside and rightside of element boundary
    real(kind=8) :: uu(3,9) ! coefficients of element and neighbor
    integer :: ii,jj,gg
    
    SELECT CASE(stage)
        
    CASE(1)
        uu(1,:)=pElem(i-1,j)%du0(:)
        uu(2,:)=pElem(i,j)%du0(:)
        uu(3,:)=pElem(i,j-1)%du0(:)
        
    CASE(2)
        uu(1,:)=pElem(i-1,j)%du1(:)
        uu(2,:)=pElem(i,j)%du1(:)
        uu(3,:)=pElem(i,j-1)%du1(:)
        
    CASE(3)
        uu(1,:)=pElem(i-1,j)%du2(:)
        uu(2,:)=pElem(i,j)%du2(:)
        uu(3,:)=pElem(i,j-1)%du2(:)
            
        
    !CASE(4)
    !    uu(1,:,:)=pElem(i-1,j)%du3(:,:)
    !    uu(2,:,:)=pElem(i,j)%du3(:,:)
    !    uu(3,:,:)=pElem(i,j-1)%du3(:,:)

    END SELECT
    
      ! x directon
          
              DO jj=1,3 ! for y -direction 3 subelements                  
                  
                   CALL GAUSS3(pElem(i,j)%xg(1,jj,:),pElem(i,j)%distxg(1,jj,:),pointg,wei)
                   
              DO gg=1,3 ! for 3 gauss points
               xc(1)=pElem(i-1,j)%x(1)
               xc(2)=pElem(i,j)%x(2) ! center of (i-1,j)
               xg(1)=peNode(i,j)%node(1)    ! for x_{i-1/2}
               xg(2)=pointg(2,gg)    ! for each gauss points
               
               CALL BASISELEM(xc,xg,face_l)
               
               xc(1)=pElem(i,j)%x(1) !center of (i,j)
               CALL BASISELEM(xc,xg,face_r)
               
           statel=0.d0
           stater=0.d0
           DO n=1,9
               statel=statel+uu(1,n)*face_l(n)
               stater=stater+uu(2,n)*face_r(n)  
           END DO
           
           ! HLLC flux
           !CALL HLLC1(statel,stater,pFace_x(i,j)%hgflux(jj,gg))   !flux for (x_{i-1/2},g_{gg})
           pFace_x(i,j)%hgflux(jj,gg) = statel
           END DO
              END DO
                
  
      ! y directon     

      
     DO ii=1,3 ! for x-direction 3 subelements
                  
                   CALL GAUSS3(pElem(i,j)%xg(ii,1,:),pElem(i,j)%distxg(ii,1,:),pointg,wei)
                   
           DO gg=1,3 ! for 3 gauss points
               xc(1)=pElem(i,j)%x(1)
               xc(2)=pElem(i,j-1)%x(2) !center of (i,j-1)
               xg(1)=pointg(1,gg)       ! for each gauss points
               xg(2)=peNode(i,j)%node(2)       ! for y_{j-1/2}
               
               CALL BASISELEM(xc,xg,face_l)
               
               xc(2)=pElem(i,j)%x(2)  !center of (i,j)
               CALL BASISELEM(xc,xg,face_r)
               
               statel=0.d0
               stater=0.d0
           DO n=1,9
               statel=statel+uu(3,n)*face_l(n)
               stater=stater+uu(2,n)*face_r(n)
           END DO
           
           !CALL HLLC1(statel,stater,fluxg)          
           pFace_y(i,j)%hgflux(ii,gg) = statel
           
           END DO
     END DO
       
    END SUBROUTINE
!-------------------------------------------------------------------
         
!-------------------------------------------------------------------
    SUBROUTINE GETHLGFLUX(i,j)
    USE ElementType
    use GlobalData    
    IMPLICIT REAL*8 (A-H,O-Z)
    real(kind=8) :: pointg(2,3),xc(2),xg(2),wei(3)
    real(kind=8) :: state ! conservative variables(global) on Gauss points 
    real(kind=8) :: prim ! local  prim variables
    real(kind=8) :: fluxg ! local flux computed by local conservative variables
    real(kind=8) :: face(9) ! the value of basis on Gauss points of leftside and rightside of element boundary
    real(kind=8) :: uu(9) ! coefficients of element and neighbor
    integer :: ii,jj,kk,gg
    
    
    
    
    SELECT CASE(stage)
    CASE(1)
        uu(:)=pElem(i,j)%du0(:)
    CASE(2)
        uu(:)=pElem(i,j)%du1(:)
    CASE(3)
        uu(:)=pElem(i,j)%du2(:)
    CASE(4)
        uu(:)=pElem(i,j)%du3(:)
    END SELECT
    
    ! x-direction
               xc(1)=pElem(i,j)%x(1)
               xc(2)=pElem(i,j)%x(2)
    DO jj=1,3 ! y-direction
                  
      CALL GAUSS3(pElem(i,j)%xg(2,jj,:),pElem(i,j)%distxg(2,jj,:),pointg,wei)
  
    DO gg=1,3 ! y-direction three gauss points in each subelement
       DO kk=2,3 ! x-direction last two subunits
               xg(1)=pElem(i,j)%subnode(kk,jj,1)
               xg(2)=pointg(2,gg)
               
               CALL BASISELEM(xc,xg,face)
           state=0.d0
           DO n=1,9
               state=state+uu(n)*face(n)
           END DO
           
            
            fluxg=state
 
            pFace_x(i,j)%hlgflux(kk,jj,gg)= fluxg
       END DO
   
    END DO
    END DO
    
    DO jj=1,3
        DO gg=1,3
            pFace_x(i,j)%hlgflux(1,jj,gg)=pFace_x(i,j)%hgflux(jj,gg) ! x_{i-1/2} y-> 3*3 gauss points
            pFace_x(i,j)%hlgflux(4,jj,gg)=pFace_x(i+1,j)%hgflux(jj,gg) ! x_{i+1/2}  y-> 3*3 gauss points
        END DO
    END DO
      
     ! y-direction
    
    DO ii=1,3 ! x-direction
                  
      CALL GAUSS3(pElem(i,j)%xg(ii,2,:),pElem(i,j)%distxg(ii,2,:),pointg,wei)
      
    DO kk=2,3 ! x-direction last two subunits
    DO gg=1,3 ! y-direction two gauss points in each subelement
               xg(1)=pointg(1,gg)
               xg(2)=pElem(i,j)%subnode(ii,kk,2)
               
               CALL BASISELEM(xc,xg,face)
               
           state=0.d0
           DO n=1,9
               state=state+uu(n)*face(n)
           END DO
           
          
            fluxg=state
 

            pFace_y(i,j)%hlgflux(kk,ii,gg)= fluxg
    END DO
    END DO
    END DO
    
    DO ii=1,3
        DO gg=1,3
            pFace_y(i,j)%hlgflux(1,ii,gg)=pFace_y(i,j)%hgflux(ii,gg) ! y_{j-1/2} x-> 9 gauss points
            pFace_y(i,j)%hlgflux(4,ii,gg)=pFace_y(i,j+1)%hgflux(ii,gg) ! y_{j+1/2}  x-> 9 gauss points
        END DO
    END DO
       
   END SUBROUTINE
!-------------------------------------------------------------------
  
    
!-------------------------------------------------------------------
    SUBROUTINE VERTICESJUMP(i,j)
    USE ElementType
    use GlobalData    
    IMPLICIT REAL*8 (A-H,O-Z)
    real(kind=8) :: pnode(2),pointg(2,2),xc(2),xg(2)
    real(kind=8) :: statel(6), stater(6),charist_statel(6),charist_stater(6) ! conservative variables(global) on Gauss points of leftside and rightside of element boundary
    real(kind=8) :: statel_1(6), stater_1(6),charist_statel_1(6),charist_stater_1(6) ! conservative variables(global) on Gauss points of leftside and rightside of element boundary
    real(kind=8) :: face_l(6,9),face_r(6,9) ! the value of basis on Gauss points of leftside and rightside of element boundary
    real(kind=8) :: partial_jump_x(6,2),partial_jump_y(6,2)
    real(kind=8) :: uu(3,9) ! coefficients of element and neighbor
    integer :: ii,m,n

    
    SELECT CASE(stage)
        
    CASE(1) 
        uu(1,:)=pElem(i-1,j)%du0(:)
        uu(2,:)=pElem(i,j)%du0(:)
        uu(3,:)=pElem(i,j-1)%du0(:)
        
    CASE(2)
        uu(1,:)=pElem(i-1,j)%du1(:)
        uu(2,:)=pElem(i,j)%du1(:)
        uu(3,:)=pElem(i,j-1)%du1(:)
        
    CASE(3)
        uu(1,:)=pElem(i-1,j)%du2(:)
        uu(2,:)=pElem(i,j)%du2(:)
        uu(3,:)=pElem(i,j-1)%du2(:)
        
    CASE(4)
        uu(1,:)=pElem(i-1,j)%du3(:)
        uu(2,:)=pElem(i,j)%du3(:)
        uu(3,:)=pElem(i,j-1)%du3(:)
            
    END SELECT
    
           
! compute for vertices jump

       !x directon     
        pointg(2,1)=peNode(i,j)%node(2)  !y_{j-1/2} 
        pointg(2,2)=peNode(i,j+1)%node(2)  !y_{j+1/2} 
       
           DO ii=1,2
           
               xc(1)=pElem(i-1,j)%x(1) !\tau_{i-1,j} x-center
               xc(2)=pElem(i,j)%x(2) !\tau_{i-1,j} \tau_{i,j} y-center
               xg(1)=peNode(i,j)%node(1)  !x_{i-1/2} boundary
               xg(2)=pointg(2,ii)         ! gauss point term to y not \eta
               

               ! partial^\alpha    
            CALL BASISELEM(xc,xg,face_l(1,:))  ! \alpha=0
            CALL BASISELEM_de(xc,xg,face_l(2:3,:))   ! \alpha=1
            CALL BASISELEM_de2(xc,xg,face_l(4:6,:))   ! \alpha=2
                 
               
               xc(1)=pElem(i,j)%x(1)             !\tau_{i,j} x-center
            CALL BASISELEM(xc,xg,face_r(1,:))  ! \alpha=0
            CALL BASISELEM_de(xc,xg,face_r(2:3,:))   ! \alpha=1
            CALL BASISELEM_de2(xc,xg,face_r(4:6,:))   ! \alpha=2
               
           statel(:)=0.d0
           stater(:)=0.d0
           DO m=1,6
           DO n=1,9
               statel(m)=statel(m)+uu(1,n)*face_l(m,n)     ! \alpha=0: m=1,\alpha=1:m=2,3,\alpha=2:m=4,5,6
               stater(m)=stater(m)+uu(2,n)*face_r(m,n)
           END DO
           END DO
           
            DO n = 1,6 
                !partial_jump_x(m,n,ii) = (stater(m,n) - statel(m,n))**2
                partial_jump_x(n,ii) = (stater(n) - statel(n))**2
            END DO
           
        END DO
            pFace_x(i,j)%partial_jump(:,:) = partial_jump_x(:,:) !local to global

      ! y directon  

        pointg(1,1)=peNode(i,j)%node(1)  !x_{i-1/2} 
        pointg(1,2)=peNode(i+1,j)%node(1)  !x_{i+1/2} 

           DO ii=1,2
               
               xc(1)=pElem(i,j)%x(1)       !\tau_{i,j} \tau_{i,j-1} x-center
               xc(2)=pElem(i,j-1)%x(2)       !\tau_{i,j-1} y-center
               xg(1)=pointg(1,ii)               !x_{i-1/2} or !x_{i+1/2} 
               xg(2)=peNode(i,j)%node(2)              !y_{j-1/2}

               ! partial^\alpha     \alpha=0

            CALL BASISELEM(xc,xg,face_l(1,:))  ! \alpha=0
            CALL BASISELEM_de(xc,xg,face_l(2:3,:))   ! \alpha=1
            CALL BASISELEM_de2(xc,xg,face_l(4:6,:))   ! \alpha=2
                       
               xc(2)=pElem(i,j)%x(2)       !\tau_{i,j}      y-center
            CALL BASISELEM(xc,xg,face_r(1,:))  ! \alpha=0
            CALL BASISELEM_de(xc,xg,face_r(2:3,:))   ! \alpha=1
            CALL BASISELEM_de2(xc,xg,face_r(4:6,:))   ! \alpha=2
               
   
          statel(:)=0.d0
          stater(:)=0.d0
    
          DO m=1,6
          DO n=1,9
               statel(m)=statel(m)+uu(3,n)*face_l(m,n)     ! \alpha=0
               stater(m)=stater(m)+uu(2,n)*face_r(m,n)
          END DO
          END DO
   
               
            DO n = 1,6
                partial_jump_y(n,ii) = (stater(n) - statel(n))**2
            END DO
                            
           END DO
           
            pFace_y(i,j)%partial_jump(:,:) = partial_jump_y(:,:) !local to global

     
    END SUBROUTINE
!-------------------------------------------------------------------
!-------------------------------------------------------------------
    SUBROUTINE ELEMFACE(i,j)
    USE ElementType
    use GlobalData    
    IMPLICIT REAL*8 (A-H,O-Z)
    real(kind=8) :: wei(3)
    integer :: ii,jj,gg
    
    wei(1)=5.d0/9.d0
    wei(2)=8.d0/9.d0
    wei(3)=5.d0/9.d0
    
    pElem(i,j)%fflux_x(:)=0.d0
    pElem(i,j)%fflux_y(:)=0.d0
    
    DO ii=1,3 !x-direction
        DO jj=1,3 !y-direction
            DO gg=1,3 !gauss points
               
                pElem(i,j)%fflux_x((jj-1)*3+ii)=pElem(i,j)%fflux_x((jj-1)*3+ii)+ (-pFace_x(i,j)%hlgflux(ii+1,jj,gg)+&
                                                                    pFace_x(i,j)%hlgflux(ii,jj,gg)) * wei(gg) *pElem(i,j)%distxg(ii,jj,2)/2.d0
                pElem(i,j)%fflux_y((jj-1)*3+ii)=pElem(i,j)%fflux_y((jj-1)*3+ii)+ (-pFace_y(i,j)%hlgflux(jj+1,ii,gg)+&
                                                                    pFace_y(i,j)%hlgflux(jj,ii,gg)) * wei(gg) * pElem(i,j)%distxg(ii,jj,1)/2.d0 
                
            END DO
        END DO
    END DO
   
    
    END SUBROUTINE
!-------------------------------------------------------------------
!-------------------------------------------------------------------
!-----------------------------------------------------------------
    SUBROUTINE SHARPCORRECTION(i,j)
    USE ElementType
    use GlobalData    
    IMPLICIT REAL*8 (A-H,O-Z)
    real(kind=8) :: uu(9)
    real(kind=8) :: xigema(3),jump_for_vertices,sharp_corection(9)
    
    integer :: ii,jj,ll

    SELECT CASE(stage)
    CASE(1)
        uu(:)=pElem(i,j)%du0(:)
    CASE(2)
        uu(:)=pElem(i,j)%du1(:)
    CASE(3)
        uu(:)=pElem(i,j)%du2(:)
    CASE(4)
        uu(:)=pElem(i,j)%du3(:)
    END SELECT
        
        ! calculate for sharp correction
        
        xigema(:) = 0.d0
        ! \alpha^l_k
                !l=0
            jump_for_vertices= pFace_x(i,j)%partial_jump(1,1) + pFace_x(i,j)%partial_jump(1,2)&
                                   + pFace_x(i+1,j)%partial_jump(1,1) + pFace_x(i+1,j)%partial_jump(1,2)&
                                   + pFace_y(i,j)%partial_jump(1,1) + pFace_y(i,j)%partial_jump(1,2)&
                                   + pFace_y(i,j+1)%partial_jump(1,1) + pFace_y(i,j+1)%partial_jump(1,2)

                xigema(1) = xigema(1) + 2.d0*1.d0/(3.d0*1.d0)* DSQRT(DX**2+DY**2)**(0-1) * DSQRT(1.d0/4.d0*jump_for_vertices)
                !xigema(1) = xigema(1) + 2.d0*1.d0/(3.d0*1.d0)* DSQRT(DX**2)**(0-1) * DSQRT(1.d0/4.d0*jump_for_vertices)

        
                !l=1
            DO ll = 2,3
        
            jump_for_vertices= pFace_x(i,j)%partial_jump(ll,1) + pFace_x(i,j)%partial_jump(ll,2)&
                                   + pFace_x(i+1,j)%partial_jump(ll,1) + pFace_x(i+1,j)%partial_jump(ll,2)&
                                   + pFace_y(i,j)%partial_jump(ll,1) + pFace_y(i,j)%partial_jump(ll,2)&
                                   + pFace_y(i,j+1)%partial_jump(ll,1) + pFace_y(i,j+1)%partial_jump(ll,2)

                xigema(2) = xigema(2) + 2.d0*3.d0/(3.d0*1.d0)* DSQRT(DX**2+DY**2)**(1-1) * DSQRT(1.d0/4.d0*jump_for_vertices)
                !xigema(2) = xigema(2) + 2.d0*3.d0/(3.d0*1.d0)* DSQRT(DX**2)**(1-1) * DSQRT(1.d0/4.d0*jump_for_vertices)

            
            END do
                !l=2
            DO ll = 4,6
            jump_for_vertices= pFace_x(i,j)%partial_jump(ll,1) + pFace_x(i,j)%partial_jump(ll,2)&
                                   + pFace_x(i+1,j)%partial_jump(ll,1) + pFace_x(i+1,j)%partial_jump(ll,2)&
                                   + pFace_y(i,j)%partial_jump(ll,1) + pFace_y(i,j)%partial_jump(ll,2)&
                                   + pFace_y(i,j+1)%partial_jump(ll,1) + pFace_y(i,j+1)%partial_jump(ll,2)

                xigema(3) = xigema(3) + 2.d0*5.d0/(3.d0*2.d0)* DSQRT(DX**2+DY**2) **(2-1) * DSQRT(1.d0/4.d0*jump_for_vertices)
                !xigema(3) = xigema(3) + 2.d0*5.d0/(3.d0*2.d0)* DSQRT(DX**2) **(2-1) * DSQRT(1.d0/4.d0*jump_for_vertices)

        
            END do
  
       
       xigema=1.0d0*xigema
        ! correction part
        
        sharp_corection(:) = 0.d0

                !l=0
        DO n=2,9
            sharp_corection(n)=sharp_corection(n) + xigema(1) * uu(n)
        END DO
        
                !l=1
        DO n=2,9
            sharp_corection(n)=sharp_corection(n) + xigema(2) * uu(n)
        END DO
        
                !l=2
        DO n=5,9
            sharp_corection(n)=sharp_corection(n) + xigema(3) * uu(n)
        END DO

        pElem(i,j)%correction(:) = sharp_corection(:)
        !write(*,*) i,j, pElem(i,j)%correction(1,3)
        
    END SUBROUTINE
!-------------------------------------------------------------------

!-------------------------------------------------------------------
    SUBROUTINE UPDATE(i,j)
    USE ElementType
    use GlobalData    
    IMPLICIT REAL*8 (A-H,O-Z)
    real(kind=8) :: dtu(9)
    
      pElem(i,j)%dtu_x(:)=0.d0
      pElem(i,j)%dtu_y(:)=0.d0
      
       dtu(:)=0.d0

      DO n=1,9
      DO k=1,9   
          dtu(n)= dtu(n)+stiff(n,k)*(pElem(i,j)%fflux_x(k)+pElem(i,j)%fflux_y(k))*4.d0/(DX*DY)
          !pElem(i,j)%dtu_x(m,n)= pElem(i,j)%dtu_x(m,n)+stiff(n,k)*(pElem(i,j)%fflux_x(m,k))*4.d0/(DX*DY)
          !pElem(i,j)%dtu_y(m,n)= pElem(i,j)%dtu_y(m,n)+stiff(n,k)*(pElem(i,j)%fflux_y(m,k))*4.d0/(DX*DY)
      END DO
      END DO

      
      dtu(:)= dtu(:) - pElem(i,j)%correction(:)

    SELECT CASE(stage)
    CASE(1)
        pElem(i,j)%du1(:)=pElem(i,j)%du0(:)+TAU*dtu(:)
    CASE(2)
        pElem(i,j)%du2(:)=3.d0/4.d0*pElem(i,j)%du0(:)+1.d0/4.d0*pElem(i,j)%du1(:)+1.d0/4.d0*TAU*dtu(:)
    CASE(3)
        pElem(i,j)%du0(:)=1.d0/3.d0*pElem(i,j)%du0(:)+2.d0/3.d0*pElem(i,j)%du2(:)+2.d0/3.d0*TAU*dtu(:)
    END SELECT

    !SELECT CASE(stage)
    !CASE(1)
    !    pElem(i,j)%L1(:,:)=dtu(:,:)
    !    !write(*,*) i,j,pElem(i,j)%L1(1,1)
    !    pElem(i,j)%du1(:,:)=pElem(i,j)%du0(:,:)+1.d0/2.d0*TAU*pElem(i,j)%L1(:,:)
    !CASE(2)
    !    pElem(i,j)%L2(:,:)=dtu(:,:)
    !    pElem(i,j)%du2(:,:)=pElem(i,j)%du0(:,:)+1.d0/2.d0*TAU*pElem(i,j)%L2(:,:)
    !CASE(3)
    !    pElem(i,j)%L3(:,:)=dtu(:,:)
    !    pElem(i,j)%du3(:,:)=pElem(i,j)%du0(:,:)+TAU*pElem(i,j)%L3(:,:)
    !CASE(4)
    !    pElem(i,j)%L4(:,:)=dtu(:,:)
    !    pElem(i,j)%du0(:,:)=pElem(i,j)%du0(:,:)+1.d0/6.d0*TAU*(pElem(i,j)%L1(:,:)+2.d0*pElem(i,j)%L2(:,:)+2.d0*pElem(i,j)%L3(:,:)+pElem(i,j)%L4(:,:))
    !END SELECT 
      
    END SUBROUTINE
!-------------------------------------------------------------------

      PROGRAM STRUCTURE_RKDG
      USE ElementType      
      USE GlobalData 
      use omp_lib
      IMPLICIT REAL*8 (a-h,o-z)   
      integer :: clock_count1,clock_count2,clock_max,clock_rate

      ! start counting time
      call system_clock ( clock_count1, clock_rate, clock_max  )
      
      ! set mesh & inital coefficients
	  CALL INITIAL   
      
      ! start computing
      DO WHILE(TIME.LE.TSTOP)

      ! get time step
      CALL GETTIMESTEP  
      
      ! output: numbers of computation, current time, current time step
      WRITE(*,*),IT,TIME,TAU 
     

	  TIME=TIME+TAU
      IT=IT+1
      
      ! output picture
      ! the number of picture depends on the setting of subroutine INITIAL
      ! eg: UTIME = TSTOP/5.d0 output 0 1 2 3 4 5
      IF(TIME.GE.SSTOP) THEN    
      CALL OUTPUT(SSTOP,UTIME)
      SSTOP=SSTOP+UTIME
      ENDIF  
     
          ! update the coefficients of this time step
          ! stage(global variable) manages the state of rk3
          DO stage=1,3
              
          ! set boundary
          CALL BOUNDARY
          
          ! 'omp parallel' is the sign of the beginning of parallel computation
          ! 'omp end parallel' is the sign of the end of parallel computation
          ! 'omp do' & 'omp end do nowait' is sign of do loop
          
          
          ! compute the numerical flux of boundaries
          ! boundary: for x direction, compute the leftside
          !           for y direction, compute the downside
          !$omp parallel
          !$omp do         
          DO i=2,NP1
          DO j=2,NP2
          CALL GETHFLUX(i,j)             
          END DO
          END DO
          !$omp end do nowait
          !$omp end parallel
           !write(*,*)  pFace_y(3,3)%hgflux(1,2,1)
          !$omp parallel
          !$omp do         
          DO i=2,NP1-1
          DO j=2,NP2-1
          CALL GETHLGFLUX(i,j)           
          END DO
          END DO
          !$omp end do nowait
          !$omp end parallel
           !write(*,*)  pFace_x(3,3)%hlgflux(2,1,1,:)
          ! compute the jump^2 at vertices in each element
          !$omp parallel
          !$omp do         
          DO i=2,NP1
          DO j=2,NP2
          CALL VERTICESJUMP(i,j)             
          END DO
          END DO
          !$omp end do nowait
          !$omp end parallel
          
          ! compute the boundary integral
          !$omp parallel
          !$omp do
          DO i=2,NP1-1    
          DO j=2,NP2-1
          CALL ELEMFACE(i,j)
          END DO
          END DO
          !$omp end do nowait
          !$omp end parallel
      
      
          ! compute the correction term (for 2D, surface)
          !$omp parallel
          !$omp do 
          DO i=2,NP1-1    
          DO j=2,NP2-1
          CALL SHARPCORRECTION(i,j)
          END DO
          END DO
          !$omp end do nowait
          !$omp end parallel 
          
          ! update coefficients (subroutine rk3)
          !$omp parallel
          !$omp do 
          DO i=2,NP1-1    
          DO j=2,NP2-1
          CALL UPDATE(i,j)
          END DO
          END DO
          !$omp end do nowait
          !$omp end parallel 
          END DO

      END DO 
      
                   !write(*,*) pElem(25,2)%du0(1,:)

      
    ! output errors eg: L^2 errors
    CALL OUTPUT_error
      
    ! stop counting time
    call system_clock ( clock_count2, clock_rate, clock_max  )

    
    write(*,*) " The program's calculation time is", DBLE((clock_count2-clock_count1)/clock_rate), "seconds"

    pause

      
    END PROGRAM   




!----------------------------------------------------------------------------  
!----------------------------------------------------------------------------  
  SUBROUTINE OUTPUT(STIM,UTIM) 
    USE ElementType      
    USE GlobalData  
    IMPLICIT REAL*8(A-H,O-Z)   
    INTEGER :: i,NUM,IDX
    REAL(kind=8) :: STIM, UTIM
    REAL(kind=8) :: x,y
    REAL(kind=8) :: pointg(2,3),xc(2),xg(2),elem(9),wei(3)
    REAL(kind=8) :: u_ave,rho,u1
    REAL(kind=8) :: u_g,rhog
    
    integer :: ii,jj
    CHARACTER CHR*26,CHR1*26,CHR2*26,CHR6*26,CHR7*26

    NUM=IDINT(STIM/UTIM)
    WRITE(CHR,'(I8)') NUM
    DO WHILE (INDEX(CHR,' ').EQ.1)
    CHR=CHR(2:)
    ENDDO
    IDX=INDEX(CHR,' ')-1       
    CHR2='RKSV-2d-sin-16-'//CHR(:IDX)//'.plt'   
    
    OPEN(8,FILE=CHR2, STATUS='UNKNOWN') 
 
    WRITE(8,*) 'title="contour"'
    WRITE(8,*) 'variables="x","y","den","u","v","press"' 
    WRITE(8,*) 'ZONE',' I=',(NP1-2)*1,' J=',(NP2-2)*1,' f=point'
 
    
    DO IY=2,NP2-1 
    DO IX=2,NP1-1 
    
    u_ave=0.d0
    u_g=0.d0
    
    xc(:) = pElem(IX,IY)%x(:)
      CALL GAUSS3(xc,pElem(IX,IY)%distx(:),pointg,wei)
    
    
    DO ii=1,3
    DO jj=1,3
        xg(1)=pointg(1,ii)
        xg(2)=pointg(2,jj)
        
        CALL BASISELEM(xc,xg,elem)
        
    
        
        rho=0.d0
        DO n=1,9
            rho=rho+pElem(IX,IY)%du0(n)*elem(n)
        END DO
        
        
        u_ave=u_ave+rho*wei(ii)*wei(jj)*DX*DY/4.d0
        
    
        
    END DO
    END DO
        u_ave=u_ave/DX/DY
    

    WRITE(8,11) xc(1),xc(2),pElem(IX,IY)%du0(1)
    
        
    ENDDO      
    ENDDO
 
11  FORMAT(8(F14.8, 2X))  
    
    END SUBROUTINE 
        
!----------------------------------------------------------------------------  
!!----------------------------------------------------------------------------  
!  !output full type
!    SUBROUTINE OUTPUT(STIM,UTIM) 
!    USE ElementType      
!    USE GlobalData  
!    IMPLICIT REAL*8(A-H,O-Z)   
!    INTEGER :: i,NUM,IDX
!    REAL(kind=8) :: STIM, UTIM
!    REAL(kind=8) :: x,y,pointy
!    REAL(kind=8) :: xnode(2),xc(2),xg(2),dist(2),elem(9),wei(3),pointx(21)
!    REAL(kind=8) :: u_ave(4),rho(4),u1(4)
!    REAL(kind=8) :: u_g(4),rhog(4)
!    
!    integer :: ii,jj
!    CHARACTER CHR*26,CHR1*26,CHR2*26,CHR6*26,CHR7*26
!
!    NUM=IDINT(STIM/UTIM)
!    WRITE(CHR,'(I8)') NUM
!    DO WHILE (INDEX(CHR,' ').EQ.1)
!    CHR=CHR(2:)
!    ENDDO
!    IDX=INDEX(CHR,' ')-1       
!    CHR2='RKSV-full-lax-2-128-'//CHR(:IDX)//'.plt'   
!    
!    OPEN(8,FILE=CHR2, STATUS='UNKNOWN') 
! 
!    WRITE(8,*) 'title="contour"'
!    WRITE(8,*) 'variables="x","y","den","u","v","press"' 
!    WRITE(8,*) 'ZONE',' I=',(NP1-2)*21,' J=',(NP2-2)*1,' f=point'
! 
!    
!    DO IY=2,NP2-1 
!    DO IX=2,NP1-1 
!    
!    u_ave=0.d0
!    u_g=0.d0
!    
!    xnode(:) = peNode(IX,IY)%node(:)
!    xc(:) = pElem(IX,IY)%x(:)
!    dist(:) = pElem(IX,IY)%distx(:)
!    
!    DO ii=1,21
!        pointx(ii) = xnode(1) + DBLE(ii)/21.d0*dist(1)
!    END DO
!    
!    pointy = xc(2)
!    
!    
!    
!      !CALL GAUSS3(xc,pElem(IX,IY)%distx(:),pointg,wei)
!    
!    
!    DO ii=1,21
!        xg(1)=pointx(ii)
!        xg(2)=pointy
!        
!        CALL BASISELEM(xc,xg,elem)
!        
!        rho=0.d0
!        
!        DO n=1,9
!            rho=rho+pElem(IX,IY)%du0(:,n)*elem(n)
!        END DO
!        
!        CALL exchange(u1,rho,"c2p")
!                
!    WRITE(8,11) xg(1),xg(2),u1(1),u1(2),u1(3),u1(4)
!        
!    END DO
!    
!    END DO
!    END DO      
! 
!11  FORMAT(8(F14.8, 2X))  
!    
!    END SUBROUTINE 
!        
!!----------------------------------------------------------------------------  
    
!----------------------------------------------------------------------------  
 
SUBROUTINE OUTPUT_error
USE ElementType      
USE GlobalData 
IMPLICIT REAL*8 (a-h,o-z)   
REAL(KIND=8) :: error1,error2,error3,error4 ! L^1 error, L^2 error
REAL(KIND=8) :: temp1,temp2,temp3,temp4
REAL(KIND=8) :: exact, den, exact_n, den_n ! exact solution, numerical solution
REAL(KIND=8) :: elem(9),elem_n(9)
REAL(KIND=8) :: pointg(2,3),xc(2),xg(2),wei(3),node(2)
integer :: ii,jj

error1=0.d0
error2=0.d0
error3=0.d0
error4=0.d0


DO i=2,NP1-1  
DO j=2,NP2-1    
   
   
   
    temp1=0.d0
    temp2=0.d0
    temp3=0.d0
    
    xc(:)=pElem(i,j)%x(:)
      CALL GAUSS3(xc,pElem(i,j)%distx(:),pointg,wei)
    
    DO ii=1,3
    DO jj=1,3
        xg(1)=pointg(1,ii)
        xg(2)=pointg(2,jj)
        
        CALL BASISELEM(xc,xg,elem)
        
        exact=0.d0
        den=0.d0
        DO n=1,9
            !exact=exact+pElem(i,j)%initial(1,n)*elem(n)
            den=den+pElem(i,j)%du0(n)*elem(n)
        END DO
    exact =  1.d0+0.2d0*DSIN(PI*(xg(1)+xg(2)))
    !exact =  1.d0+0.2d0*DSIN(PI*(xg(1)))
        
        temp1=temp1+ABS(den-exact)*wei(ii)*wei(jj)
        temp2=temp2+ABS(den-exact)**2*wei(ii)*wei(jj)
        temp3=temp3+(exact-den)*wei(ii)*wei(jj)
        
    END DO
    END DO
    temp1=temp1*DX*DY/4.d0
    temp2=temp2*DX*DY/4.d0
    temp3=temp3/4.d0
    temp3=temp3**2
    
    error1=error1+temp1 ! L^1 error
    error2=error2+temp2 ! L^2 error
    error3=error3+temp3 ! cell average
   
END DO
END DO

error3=error3/DBLE((NP1-2)*(NP2-2))


DO i=2,NP1-1  
DO j=2,NP2-1  
    
    xc(:)=pElem(i,j)%x(:)
    
    node(:)= peNode(i+1,j+1)%node(:)
    
    CALL BASISELEM(xc,node,elem_n)
    
     exact_n=0.d0
     den_n=0.d0

        DO n=1,9
            !exact_n=exact_n+pElem(i,j)%initial(1,n)*elem_n(n)
            den_n=den_n+pElem(i,j)%du0(n)*elem_n(n)
        END DO
        
    exact_n =  1.d0+0.2d0*DSIN(PI*(node(1)+node(2)))
    !exact_n =  1.d0+0.2d0*DSIN(PI*(node(1)))
    
    temp4=(exact_n-den_n)**2
     
    error4 = error4+temp4
        
END DO
END DO

error4=error4/DBLE((NP1-2)*(NP2-2))

PRINT *,'mesh:',NP1-2,NP2-2
PRINT *,'c-a error:', DSQRT(error3)
PRINT *,'L^2 error:', DSQRT(error2)
PRINT *,'node error:', DSQRT(error4)


      
END SUBROUTINE

!----------------------------------------------------------------------------  
!----------------------------------------------------------------------------  