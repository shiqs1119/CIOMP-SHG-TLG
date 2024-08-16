
Program workfile
    Implicit none
! some parameters    
    complex(8) , parameter :: i = (0.0, 1.0) , t0 = (0.0,0.0)
    real(8) , parameter :: PI = acos(-1.0) , a0 = 2.46D-010
    real(8) , parameter :: t1 = 3.12D0 , t2 = 0.29D0 , t3 = -0.0103D0 !expressed in eV
    real(8) , parameter :: h_bar = 1.054572663D-034 , r0 = sqrt(3.0)/3.0D0
! variables of the k value    
    integer , parameter :: q = 1000
    real(8) :: b1(2), b2(2) , a1(2) , a2(2) 
    real(8) , allocatable :: k(:,:,:)
    integer :: x , y
! the energy band and length
    real(8) :: E(6), lx
    real(8), parameter :: l = (4*sqrt(2.0d0)+2*sqrt(5.0d0))/6.0d0
   
    allocate(k(2,q,q))
    b1 = [ 2*PI/(a0*sqrt(3.0)) , 2*PI/a0 ]
    b2 = [ 2*PI/(a0*sqrt(3.0)) , -2*PI/a0 ]
    a1 = [ a0*sqrt(3.0)/2 , a0/2 ]
    a2 = [ a0*sqrt(3.0)/2 , -a0/2 ]
    Forall(x=1:size(k,2), y=1:size(k,3))
            k(:,x,y) = (b1*x+b2*y)/q
    End Forall

! about the file
open(1, file = 'energyband.xls', status = 'new')
! the energy band
do x = 1,999
   if(x<=500) then
               call eigen(x, x, E)
               lx = (real(x)*sqrt(2.0)/1000)/l
       elseif (x<=667) then
               call eigen((1000-x), x, E)
               lx = (real(x)*sqrt(2.0)/1000)/l
       elseif (x<=1000) then
               call eigen((1000-x),(2000-2*x), E)
               lx =(real(x)*sqrt(5.0d0)/1000-2.0d0*sqrt(5.0)/3.0d0+2.0d0*sqrt(2.0d0)/3.0)/l
   endif
   write(1,*)x,lx,E(1),E(2),E(3),E(4),E(5),E(6)
 enddo
 
 contains
! about the hamiltonian matrix
subroutine eigen(kx,ky,fW)
    implicit none
! the input and output 
    integer, parameter :: N = 6
    integer,intent(in)::kx,ky
    real(8), intent(out) :: fw(N)
! about the k value
    complex(8) :: m1 , m2
    complex(8) :: f1 , f2
! variables of the Eigenvalue and Eigenvector of the Hamiltonian matrix
    complex(8) :: H(N,N) 
 
    m1 = dot_product(a1 , k(:,kx,ky))
    m2 = dot_product(a2 , k(:,kx,ky))
    f1 = 1+exp(i*m1)+exp(i*m2)
    f2 = 1+exp(-i*m1)+exp(-i*m2)
    H = reshape([complex(8) :: t0, t1*f1, t0, t0, t0, t0,&
            t1*f2, t0, t0, t2, t0, t3,&
            t0, t0, t0, t1*f2, t0 ,t0,&
            t0, t2, t1*f1,t0, t0, t2, &
            t0, t0, t0, t0, t0, t1*f1,&
            t0, t3, t0, t2, t1*f2, t0],[6 , 6])
    call myzheev(N,fW,H)
endsubroutine

! eigenvalue of the zheev
subroutine myzheev(N,W,H)
    implicit none
    integer :: N
    real(8),intent(out) :: W(N)
    complex(8),intent(inout) :: H(N,N)
    integer:: INFO=0
    real(8) :: RWORK(3*N-2)
    complex(8) :: WORK(3*N)
    call ZHEEV('V', 'U', N, H, N, W, WORK, 3*N, RWORK, INFO)
 endsubroutine

    

end program workfile
