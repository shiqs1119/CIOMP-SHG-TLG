module ABA
    implicit none
    
    complex(8) , parameter :: CI = (0.0, 1.0)
    real(8) , parameter :: PI = acos(-1.0) , a0 = 2.4595D-010, e = 1.602176634D-19, eV = 1.602176634D-19
    real(8) , parameter :: h_bar = 1.054572663D-034, r0 = sqrt(3.0)/3.0D0
    real(8) , parameter :: kB = 1.380649D-023, T = 273.15D0, sigma0 = (e**2.0d0)/(4.0d0*h_bar)
    real(8) , parameter :: t1 = 3.12D0*eV , t2 = 0.39D0*eV , t3 = 0.052D0*eV    
    
    real(8), parameter :: b1(2) = [2*PI/(a0*sqrt(3.0)) , 2*PI/a0 ], b2(2) = [ 2*PI/(a0*sqrt(3.0)) , -2*PI/a0 ]
    real(8), parameter :: a1(2) = [ a0*sqrt(3.0)/2 , a0/2 ], a2(2) = [ a0*sqrt(3.0)/2 , -a0/2]

    real(8), parameter :: kprimd(2, 2) = reshape([2*PI/(a0*sqrt(3.0)) , 2*PI/a0, 2*PI/(a0*sqrt(3.0)) , -2*PI/a0 ], [2,2])
    real(8), parameter :: rprimd(2,2 ) = reshape([a0*sqrt(3.0)/2 , a0/2 ,  a0*sqrt(3.0)/2 , -a0/2], [2,2])
    
    real(8) :: chemmu = 0.d0 * eV
    
contains

    function get_H(k) result(H)
        implicit none
        real(8), intent(in)  :: k(2)
        complex(8) :: H(6, 6)
        
        complex*16 :: f1, f2
        
        H = 0.d0

        f1 = 1+exp(-2.d0*Pi*CI*k(1))+exp(-2.d0*Pi*CI*k(2))

        H(1, 2) = t1*f1;        H(2, 1) = dconjg(H(1, 2))
        H(2, 3) = t2;            H(3, 2) = dconjg(H(2, 3))
        H(2, 6) = t3;            H(6, 2) = dconjg(H(2, 6))
        H(3, 4) = t1*f1;        H(4, 3) = dconjg(H(3, 4))
        H(3, 6) = t2;            H(6, 3) = dconjg(H(3, 6))
        H(5, 6) = t1*f1;        H(6, 5) = dconjg(H(5, 6))
         
    end function get_H
    
  
  function get_r(k) result(r)
      implicit none
      real(8) , intent(in)  :: k(2)
      complex(8) :: r(2, 6, 6) 
      ! r is a three dimension array
      r = 0.d0
      r(1, 1, 1) = -a0*r0
      r(1, 4, 4) = a0*r0
      r(1, 5, 5) = -a0*r0
        
  end function get_r
    
    
  function get_v(k) result(v)
        implicit none
        real(8), intent(in)  :: k(2)
        complex(8) :: v(2, 6, 6) 
        complex(8) :: f1(2)
        
        complex(8) :: H(6, 6), r(2, 6, 6)
        integer :: i
        
        v = 0.d0
        
        f1(:) = -CI*(rprimd(:, 1)*exp(-2.d0*Pi*CI*k(1))+rprimd(:, 2)*exp(-2.d0*Pi*CI*k(2)))

        v(:, 1, 2) = t1*f1;        v(:, 2, 1) = dconjg(v(:, 1, 2))
        v(:, 3, 4) = t1*f1;        v(:, 4, 3) = dconjg(v(:, 3, 4))
        v(:, 5, 6) = t1*f1;        v(:, 6, 5) = dconjg(v(:, 5, 6))
        
        H = get_H(k)
        
        r = get_r(k)
      
        do i = 1, 2, 1
          v(i, :, :) = (v(i, :, :) - CI*(Matmul(r(i, :, :), H) - MatMul(H, r(i, :, :))))/h_bar
        end do
        
  end function get_v
    
    
  subroutine get_sigma_atk(v, en, w, numw, gamma, sigmak)
      implicit none 
      complex(8), intent(in) :: v(2, 6, 6)
      real(8), intent(in)  :: en(6)
      integer, intent(in)  :: numw
      real(8), intent(in)  :: w(numw), gamma
      complex(8), intent(out) :: sigmak(2, 2, numw)
        
        
      real(8) :: f(6)
      integer :: i1, i2, iw
      complex(8) :: tempf(2,2), temp(2,2), omega, ctmp, df(6)


      do i1 = 1, 6
        f(i1) = 1.0/(exp((en(i1)-chemmu)/(kB*T))+1.0)
        df(i1) = -1.0d0/(kB*T)*f(i1)*(1.0d0-f(i1))
      enddo

      sigmak = 0.d0
      
        
      do i1 = 1,6
        do i2 = 1,6
          
          temp(1,:) = v(1,i1,i2)*v(:,i2,i1)
          temp(2,:) = v(2,i1,i2)*v(:,i2,i1)
          
          ! we multuply the w by the unit eV, that is the frequency
          
          if(i1==i2) then
            do iw  = 1, numw
              sigmak(:, :, iw) = sigmak(:, :, iw) + temp*(df(i1) / dcmplx(w(iw), gamma))
            end do
            
          else
            ctmp = -(f(i1)-f(i2))/(en(i2)-en(i1))
            
            do iw = 1, numw, 1
              sigmak(:, :, iw) = sigmak(:, :, iw) + temp * (ctmp / dcmplx(w(iw) - (en(i2)-en(i1)), gamma))
            end do
          endif
          
        enddo
      enddo
      
      sigmak = sigmak * (CI*h_bar)
  end subroutine get_sigma_atk
  
end module ABA



Program workfile
    use ABA, only : eV
    Implicit none
 
    integer , parameter :: numk = 5000 , numw = 400

    integer :: x , y

    integer:: i1, i2, i5, i6=0

    real(8) :: w(numw) ! which is unknown

    real(8) ::  gamma
  
    complex(8) :: Sigma(2,2,numw)

    gamma = 0.05
    

    do i1 = 1,numw
      w(i1) = i1*8.d0/numw
    enddo
    
    call integralk(numk, numw, w * eV, gamma * eV, Sigma)
    
    open(1, file = 'sigma2.txt' , status = 'new')
    do i1 = 1, numw
        write(1,'(100es32.16,1x)')w(i1), real(Sigma(1, 1, i1)), real(Sigma(1, 2, i1)), real(Sigma(2, 1, i1)), real(sigma(2, 2, i1))
    end do    
    close(1)

end Program workfile


! integration of k

subroutine integralk(numk, numw, w, gamma, sigma)

    use ABA 

    implicit none
    integer, intent(in) :: numw, numk
    real(8),intent(in) :: w(numw)
    real(8), intent(in) :: gamma
    complex(8),intent(out) ::sigma(2, 2, numw)


    integer :: N = 6, i1, i2, i3
    real(8) :: k(2), en(6), ds 
    complex(8) ::  sigmak(2, 2, numw)
    complex(8) :: H(6,6),  v(2, 6, 6)

    ds = 8d0*(PI**2)/((a0**2)*sqrt(3.0)*(numk**2.d0))

    sigma = 0.d0

    do i2 = 1,numk
      ! write(*, *)i2
      do i3 = 1,numk
        
        k(1) = real(i2)/numk; k(2) = real(i3)/numk

        H = get_H(k)
        v = get_v(k)
        
        call myzheev(N, en, H)

        do i1 = 1, 2, 1
            v(i1, :, :) = matmul(matmul(conjg(transpose(H)), v(i1, :, :)), H)
        enddo

        call get_sigma_atk(v, en, w, numw, gamma, sigmak)
!        write(10, '(2(I5,1x),100(e12.6,1x))')i2, i3, v(:, 4, 4), en(4)
        sigma = sigma + sigmak
      
      enddo
!      write(10, *)
    enddo
  
    sigma = -2*ds*(e**2.d0)*sigma/((2*PI)**2.d0)

endsubroutine integralk


subroutine myzheev(N,W,H)
    implicit none
    integer,intent(in) :: N
    real(8),intent(out) :: W(N)
    complex(8),intent(inout) :: H(N,N)
    integer :: INFO=0
    real(8) :: RWORK(3*N-2)
    complex(8) :: WORK(3*N)

    call ZHEEV('V', 'U', N, H, N, W, WORK, 3*N, RWORK, INFO)

endsubroutine myzheev
