Program workfile
    Implicit none
! some parameters    
    complex(8) , parameter :: i = (0.0, 1.0) , t0 = (0.0,0.0)
    real(8) , parameter :: PI = acos(-1.0) , a0 = 2.4595D-010, e = 1.602176634D-19, eV = 1.6021766208D-19
    real(8) , parameter :: t1 = 3.12D0*eV , t2 = 0.29D0*eV , t3 = -0.0103D0*eV !expressed in eV
    real(8) , parameter :: h_bar = 1.054572663D-034 , r0 = sqrt(3.0)/3.0D0
    real(8) , parameter :: kB = 1.380649D-023, T = 273.15D0
! variables of the k value    
    integer , parameter :: q = 1000
    real(8) :: b1(2), b2(2) , a1(2) , a2(2) 
    real(8) , allocatable :: k(:,:,:)
    integer :: x , y
! variables of the count number    
    integer:: i1, i2, i3, i4, i5, i6=0
! varibales of polarization and conductivity and step-size
    complex(8) :: w2 ! which is unknown
    real(8) :: dk
! the final polarization and conductivity
    complex(8) :: S1_xx, S1_xy, S1_yx, S1_yy

! the k value
    allocate(k(2,q,q))
    b1 = [ 2*PI/(a0*sqrt(3.0)) , 2*PI/a0 ]
    b2 = [ 2*PI/(a0*sqrt(3.0)) , -2*PI/a0 ]
    a1 = [ a0*sqrt(3.0)/2 , a0/2 ]
    a2 = [ a0*sqrt(3.0)/2 , -a0/2 ]
    Forall(x=1:size(k,2), y=1:size(k,3))
            k(:,x,y) = (b1*x+b2*y)/q
    End Forall
! the step-size
    dk = 8d0*(pi**2)/((a0**2)*sqrt(3.0)*(10**6))

! output data into a file 
!open(1, file = 'sigma1.txt', status = 'new') 
! assign value to the w1 and w2
!do i5 = 0,100
    !w2 = complex((0.02d0)*eV*i5/h_bar ,0.03d0*eV/h_bar)
    w2 = complex((0.52d0)*eV/h_bar ,0.03d0*eV/h_bar)
! integral of k
    call integralk(w2,S1_xx,S1_xy,S1_yx,S1_yy)
    write(*,'(6es32.16,1x)')real(w2),real(S1_xx),real(S1_xy),real(S1_yx),real(S1_yy)
!end do
    !close(1)
    
    contains
! integration of k
subroutine integralk(w2,fS1_xx,fS1_xy,fS1_yx,fS1_yy)
    implicit none
    complex(8),intent(in) :: w2
    complex(8),intent(out) :: fS1_xx,fS1_xy,fS1_yx,fS1_yy
! definition of variables
    complex(8) :: m1 , m2
    complex(8) :: f1 , f2
! variables of the Eigenvalue and Eigenvector of the Hamiltonian matrix
    integer , parameter :: N=6
    real(8) fW(N)
    complex(8) H(N,N)
! variables of the velocity operator matrix
    complex(8) :: r_kx(6,6), r_ky(6,6), H_kx(6,6), H_ky(6,6),&
                  v_kx(6,6), v_ky(6,6), H_1(6,6),H_2(6,6), tH(6,6)
    complex(8) :: f1_kx, f1_ky, f2_kx, f2_ky,&
                  f1_kxx, f1_kxy, f1_kyy, f2_kxx, f2_kxy, f2_kyy
! variables of the assistant matrix and examine matrix
    complex(8) :: M_kxx(6,6), M_kxy(6,6), M_kyy(6,6), M_kyx(6,6),&
                  v_kxx(6,6), v_kxy(6,6), v_kyy(6,6), v_kyx(6,6)
! variables of eigenstate and differential
    real(8) :: f(6), f_nm(6,6), Omega(6,6)
    complex(8) :: r_x(6,6), r_y(6,6), v_nmx(6,6), v_nmy(6,6),&
                  M_xx(6,6) , M_xy(6,6), M_yx(6,6), M_yy(6,6)
! the final polarization and conductivity
    complex(8) ::  P1_x(6,6), P1_y(6,6),P1_xx(6,6),P1_xy(6,6),P1_yx(6,6),P1_yy(6,6),&
                P1_kxx(6,6), P1_kxy(6,6), P1_kyx(6,6), P1_kyy(6,6)
                   
   
! before the big loop
do i1 = 1,6
    do i2 = 1,6
        P1_kxx(i1,i2) = t0
        P1_kxy(i1,i2) = t0
        P1_kyx(i1,i2) = t0
        P1_kyy(i1,i2) = t0
    enddo
enddo
do i3 = 1,1000
    do i4 = 1,1000
! give a test number
    !i3 = 1
        !i4 = 2
    m1 = dot_product(a1 , k(:,i3,i4))
    m2 = dot_product(a2 , k(:,i3,i4))
    f1 = 1+exp(i*m1)+exp(i*m2)
    f2 = 1+exp(-i*m1)+exp(-i*m2)
    H = reshape([complex(8) :: t0, t1*f1, t0, t0, t0, t0,&
            t1*f2, t0, t0, t2, t0, t3,&
            t0, t0, t0, t1*f2, t0 ,t0,&
            t0, t2, t1*f1,t0, t0, t2, &
            t0, t0, t0, t0, t0, t1*f1,&
            t0, t3, t0, t2, t1*f2, t0],[6 , 6])
    H_1 = reshape([complex(8) :: &
            t0, t1*f1, t0, t0, t0, t0,&
            t1*f2, t0, t0, t2, t0, t3,&
            t0, t0, t0, t1*f2, t0 ,t0,&
            t0, t2, t1*f1,t0, t0, t2, &
            t0, t0, t0, t0, t0, t1*f1,&
            t0, t3, t0, t2, t1*f2, t0],[6 , 6])

! eigenvalue and eigenvector of the Hamiltonian matrix
call myzheev(N,fW,H)

! nonlinear polarization and conductivity
    H_2 = conjg(transpose(H))
    tH = transpose(H)     

! differential value of f1 and f2 on kx or ky
! differential of the Hamiltonian matrix on kx or ky
    f1_kx = (sqrt(3.0)/2)*a0*i*(exp(a0*i*(sqrt(3.0)*k(1,i3,i4)/2+k(2,i3,i4)/2))+&
            exp(a0*i*(sqrt(3.0)*k(1,i3,i4)/2-k(2,i3,i4)/2)))
    f2_kx = (sqrt(3.0)/2)*a0*i*(-exp(a0*i*(sqrt(3.0)*k(1,i3,i4)/2+k(2,i3,i4)/2))-&
            exp(-a0*i*(sqrt(3.0)*k(1,i3,i4)/2-k(2,i3,i4)/2)))
    f1_ky = (1.0/2.0)*a0*i*(exp(a0*i*(sqrt(3.0)*k(1,i3,i4)/2+k(2,i3,i4)/2))-&
            exp(a0*i*(sqrt(3.0)*k(1,i3,i4)/2-k(2,i3,i4)/2)))
    f2_ky = (1.0/2.0)*a0*i*(-exp(-a0*i*(sqrt(3.0)*k(1,i3,i4)/2+k(2,i3,i4)/2))+&
            exp(-a0*i*(sqrt(3.0)*k(1,i3,i4)/2-k(2,i3,i4)/2)))
! partial differential of the Hamiltonian matrix on kx and ky    
    f1_kxx = -(3.0/4.0)*(a0**2.0)*(exp(a0*i*(sqrt(3.0)*k(1,i3,i4)/2+k(2,i3,i4)/2))+&
            exp(a0*i*(sqrt(3.0)*k(1,i3,i4)/2-k(2,i3,i4)/2)))
    f1_kxy = (sqrt(3.0)/4)*(a0**2.0)*(-exp(a0*i*(sqrt(3.0)*k(1,i3,i4)/2+k(2,i3,i4)/2))+&
            exp(a0*i*(sqrt(3.0)*k(1,i3,i4)/2-k(2,i3,i4)/2)))
    f1_kyy = -(1.0/4.0)*(a0**2.0)*(exp(a0*i*(sqrt(3.0)*k(1,i3,i4)/2+k(2,i3,i4)/2))+&
            exp(a0*i*(sqrt(3.0)*k(1,i3,i4)/2-k(2,i3,i4)/2)))
    f2_kxx = -(3.0/4.0)*(a0**2.0)*(exp(-a0*i*(sqrt(3.0)*k(1,i3,i4)/2+k(2,i3,i4)/2))+&
            exp(-a0*i*(sqrt(3.0)*k(1,i3,i4)/2-k(2,i3,i4)/2)))
    f2_kxy = (sqrt(3.0)/4)*(a0**2.0)*(-exp(-a0*i*(sqrt(3.0)*k(1,i3,i4)/2+k(2,i3,i4)/2))+&
            exp(-a0*i*(sqrt(3.0)*k(1,i3,i4)/2-k(2,i3,i4)/2)))
    f2_kyy = -(1.0/4.0)*(a0**2.0)*(exp(-a0*i*(sqrt(3.0)*k(1,i3,i4)/2+k(2,i3,i4)/2))+&
            exp(-a0*i*(sqrt(3.0)*k(1,i3,i4)/2-k(2,i3,i4)/2)))
! differential matrix of Hamiltonian on ky and kx
H_kx = reshape([complex(8) :: t0, t1*f1_kx, t0, t0, t0, t0,&
        t1*f2_kx, t0, t0, t0, t0, t0,&
        t0, t0, t0, t1*f2_kx, t0, t0,&
        t0, t0, t1*f1_kx, t0, t0 ,t0,&
        t0, t0, t0, t0, t0, t1*f1_kx,&
        t0, t0, t0, t0, t1*f2_kx, t0],[6 , 6])
H_ky = reshape([complex(8) :: t0, t1*f1_ky, t0, t0, t0, t0,&
        t1*f2_ky, t0, t0, t0, t0, t0,&
        t0, t0, t0, t1*f2_ky, t0, t0,&
        t0, t0, t1*f1_ky, t0, t0 ,t0,&
        t0, t0, t0, t0, t0, t1*f1_ky,&
        t0, t0, t0, t0, t1*f2_ky, t0],[6 , 6])
! the position matrix in k space
r_kx = reshape([real(8) :: -a0*r0, t0, t0, t0, t0, t0,&
        t0, t0, t0, t0, t0, t0,&
        t0, t0, a0*r0, t0, t0 ,t0,&
        t0, t0, t0, t0, t0 ,t0,&
        t0, t0, t0, t0, -a0*r0, t0,&
        t0, t0, t0, t0, t0 ,t0],[6,6])
r_ky = reshape([real(8) :: t0, t0, t0, t0, t0, t0,&
        t0, t0, t0, t0, t0, t0,&
        t0, t0, t0, t0, t0 ,t0,&
        t0, t0, t0, t0, t0 ,t0,&
        t0, t0, t0, t0, t0 ,t0,&
        t0, t0, t0, t0, t0 ,t0],[6,6])
! velocity operator matrix
v_kx = (-i/h_bar)*(matmul(r_kx,H_1)-matmul(H_1,r_kx))+(1.0/h_bar)*H_kx  
v_ky = (-i/h_bar)*(matmul(r_ky,H_1)-matmul(H_1,r_ky))+(1.0/h_bar)*H_ky
! differential and partial differential matrix of Hamiltonian on ky and kx
v_kxx = reshape([complex(8) :: t0,(-i*sqrt(3.0)*t1/(3*h_bar))*f1_kx+(t1/h_bar)*f1_kxx, t0, t0, t0, t0,&
        (i*sqrt(3.0)*t1/(3*h_bar))*f2_kx+(t1/h_bar)*f2_kxx, t0, t0, t0, t0, t0,&
        t0, t0, t0, (i*sqrt(3.0)*t1/(3*h_bar))*f2_kx+(t1/h_bar)*f2_kxx, t0, t0,&
        t0, t0, (-i*sqrt(3.0)*t1/(3*h_bar))*f1_kx+(t1/h_bar)*f1_kxx, t0, t0, t0,&
        t0, t0, t0, t0, t0, (-i*sqrt(3.0)*t1/(3*h_bar))*f1_kx+(t1/h_bar)*f1_kxx,&
        t0, t0, t0, t0, (i*sqrt(3.0)*t1/(3*h_bar))*f2_kx+(t1/h_bar)*f2_kxx,t0],[6,6])
v_kxy = reshape([complex(8) :: t0,(-i*sqrt(3.0)*t1/(3*h_bar))*f1_ky+(t1/h_bar)*f1_kxy, t0, t0, t0, t0,&
        (i*sqrt(3.0)*t1/(3*h_bar))*f2_ky+(t1/h_bar)*f2_kxy, t0, t0, t0, t0, t0,&
        t0, t0, t0,(i*sqrt(3.0)*t1/(3*h_bar))*f2_ky+(t1/h_bar)*f2_kxy, t0, t0,&
        t0, t0, (-i*sqrt(3.0)*t1/(3*h_bar))*f1_ky+(t1/h_bar)*f1_kxy, t0, t0, t0,&
        t0, t0, t0, t0, t0,(-i*sqrt(3.0)*t1/(3*h_bar))*f1_ky+(t1/h_bar)*f1_kxy,&
        t0, t0, t0, t0, (i*sqrt(3.0)*t1/(3*h_bar))*f2_ky+(t1/h_bar)*f2_kxy, t0],[6,6]) 
v_kyx = reshape([complex(8) :: t0, (t1/h_bar)*f1_kxy, t0, t0, t0, t0,&
        (t1/h_bar)*f2_kxy, t0, t0, t0, t0, t0,&
        t0, t0, t0, (t1/h_bar)*f2_kxy, t0, t0,&
        t0, t0, (t1/h_bar)*f1_kxy, t0, t0, t0,&
        t0, t0, t0, t0, t0, (t1/h_bar)*f1_kxy,&
        t0, t0, t0, t0, (t1/h_bar)*f2_kxy,t0],[6,6])
v_kyy = reshape([complex(8) :: t0, (t1/h_bar)*f1_kyy, t0, t0, t0, t0,&
        (t1/h_bar)*f2_kyy, t0, t0, t0, t0, t0,&
        t0, t0, t0, (t1/h_bar)*f2_kyy, t0, t0,&
        t0, t0, (1.0/h_bar)*f1_kyy, t0, t0, t0,&
        t0, t0, t0, t0, t0, (t1/h_bar)*f1_kyy,&
        t0, t0, t0, t0, (t1/h_bar)*f2_kyy,t0],[6,6])
! the assistant matrix M
M_kxx = (-i/h_bar)*(matmul(r_kx,v_kx)-matmul(v_kx,r_kx))+(1.0/h_bar)*v_kxx
M_kyx = (-i/h_bar)*(matmul(r_kx,v_ky)-matmul(v_ky,r_kx))+(1.0/h_bar)*v_kyx
M_kxy = (-i/h_bar)*(matmul(r_ky,v_kx)-matmul(v_kx,r_ky))+(1.0/h_bar)*v_kxy
M_kyy = (-i/h_bar)*(matmul(r_ky,v_ky)-matmul(v_ky,r_ky))+(1.0/h_bar)*v_kyy

! the eigenstate matrix
    v_nmx = matmul(matmul(H_2,v_kx),H)
    v_nmy = matmul(matmul(H_2,v_ky),H)
    M_xx = matmul(matmul(H_2,M_kxx),H)
    M_xy = matmul(matmul(H_2,M_kxy),H)
    M_yx = matmul(matmul(H_2,M_kyx),H)
    M_yy = matmul(matmul(H_2,M_kyy),H)

! nonlinear polarization and conductivity
do i1 = 1, 6
    f(i1) = 1.0/(exp(fW(i1)/(kB*T))+1.0)
enddo
do i1 = 1,6
    do i2 = 1,6
    f_nm(i2,i1) = f(i1)-f(i2)
    ! the omega should be complex so that it can be calculate with w2, w2 is complex
    omega(i2,i1) = (fW(i2)-fW(i1))/h_bar
    if(i1==i2) then
            P1_x(i2,i1) = (i/(h_bar*w2))*(-1.0d0/(kB*T))*h_bar*v_nmx(i2,i1)*f(i1)*(1.0d0-f(i1))
            P1_y(i2,i1) = (i/(h_bar*w2))*(-1.0d0/(kB*T))*h_bar*v_nmy(i2,i1)*f(i1)*(1.0d0-f(i1))
    else
           r_x(i2,i1) = -i*v_nmx(i2,i1)/omega(i2,i1)
           r_y(i2,i1) = -i*v_nmy(i2,i1)/omega(i2,i1)
           P1_x(i2,i1) = r_x(i2,i1)*f_nm(i1,i2)/(h_bar*w2-h_bar*omega(i2,i1))
           P1_y(i2,i1) = r_y(i2,i1)*f_nm(i1,i2)/(h_bar*w2-h_bar*omega(i2,i1))
    endif
    P1_xx(i2,i1) = v_nmx(i2,i1)*P1_x(i2,i1)
    P1_xy(i2,i1) = v_nmx(i2,i1)*P1_y(i2,i1)
    P1_yx(i2,i1) = v_nmy(i2,i1)*P1_x(i2,i1)
    P1_yy(i2,i1) = v_nmy(i2,i1)*P1_y(i2,i1)
    enddo 
enddo
! summarize the n1 and n2
    P1_kxx = P1_kxx + (P1_xx/((2.0d0*pi)**2.0d0))*dk
    P1_kxy = P1_kxy + (P1_xy/((2.0d0*pi)**2.0d0))*dk
    P1_kyx = P1_kyx + (P1_yx/((2.0d0*pi)**2.0d0))*dk
    P1_kyy = P1_kyx + (P1_yy/((2.0d0*pi)**2.0d0))*dk
! first order conductivity
   enddo
enddo
! have the test
!do i1 = 1,6
    !print*,v_kx(i1,:)
    !enddo
fS1_xx = -2d0*(e**2.0d0)*sum(P1_kxx(1:6,1:6))
fS1_xy = -2d0*(e**2.0d0)*sum(P1_kxy(1:6,1:6))
fS1_yx = -2d0*(e**2.0d0)*sum(P1_kyx(1:6,1:6))
fS1_yy = -2d0*(e**2.0d0)*sum(P1_kyy(1:6,1:6))

endsubroutine

! eigenvalue of the zheev
subroutine myzheev(N,W,H)
    implicit none
    integer :: INFO=0
    real(8) :: RWORK(3*N-2)
    complex(8) :: WORK(3*N)
    integer,intent(in) :: N
    real(8),intent(out) :: W(N)
    complex(8),intent(inout) :: H(N,N)
    call ZHEEV('V', 'U', N, H, N, W, WORK, 3*N, RWORK, INFO)
 endsubroutine

End Program workfile
