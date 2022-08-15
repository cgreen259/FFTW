 program fftw2
   use, intrinsic :: ISO_C_BINDING

   implicit none
   include 'fftw3.f03'

   integer(C_INT) :: n, sign, flags
   type(C_PTR) :: fr
   complex(C_DOUBLE_COMPLEX), dimension(:), allocatable :: pout
   real (C_DOUBLE), dimension(:), allocatable :: pin
   real :: x
   integer :: i, j
   real, dimension (100,2) :: funct 

   n=16

   allocate(pin(1:n))
   allocate(pout(1:n))
   
   write(6,*)'This program calculates the fourier transform of'
   write(6,*)'given data.'
   write(6,*)'It does this by using the fastest fourier transform in the west'
   write(6,*)'library and the real data DFT from this library.'
   write(6,*)'It computes this up to sin(7x) and cos(7x).'
   write(6,*)'a_j and b_j are written to le fftw.dat'
   write(6,*)'the fourier expansion is written to function.dat'
   write(6,*)'These can be plotted with gnuplot using:'
   write(6,*)'plot fftw.dat with points'
   write(6,*)'replot function.dat with lines'
   
   !By setting flags to this it keeps the input array
   flags = FFTW_ESTIMATE

   !creating the input array in C
   pin = (/ &
        (-0.2_C_DOUBLE), &
        (-0.1_C_DOUBLE), &
        (0.3_C_DOUBLE), &
        (0.2_C_DOUBLE), &
        (0.4_C_DOUBLE), &
        (0.5_C_DOUBLE), &
        (0.0_C_DOUBLE), &
        (-0.4_C_DOUBLE), &
        (-0.4_C_DOUBLE), &
        (-0.2_C_DOUBLE), &
        (0.1_C_DOUBLE), &
        (0.2_C_DOUBLE), &
        (0.2_C_DOUBLE), &
        (0.1_C_DOUBLE), &
        (0.1_C_DOUBLE), &
        (-0.1_C_DOUBLE) &
        /)
  
   !setting up the fftw plan for real to complex
   fr=fftw_plan_dft_r2c_1d(n,pin,pout,flags)
   !calling the fftw plan
   call fftw_execute_dft_r2c(fr,pin,pout)

   pout = conjg(pout) /(n/2)
   write(6,'(2(a10,1x))') 'a_j','b_j'
   write(6,'(2(f10.6,1x))') pout(1:8)
   
   !writing the ouput form the FFTW to a .dat file
   !this out puts a_j and -b_j
   !hence the complex conjugate is taken above to get b_j
   open(unit=10, file='fftw.dat')
     write(10,'(2(f10.6,1x))') pout(1:8)
   close(10)
   
   !destrying the FFTW plan, as this is good practise
   call fftw_destroy_plan(fr)
   
   !setting the fourier part of the array to be initally equal to 0.5*a0
   funct(1:100,2)=0.5*pout(1)
   
   !calculating the fourier series to over a range of x=0.1 to x=10
   do j=0,99
	 x=real(j)/10.0
	 x=(x/100)*2*acos(0.0)
	 funct(i,1)=x
	 !do loop to calculate the fourier series for a given x value
	 do i=1,7
	   funct(j,2)=funct(j,2)+pout(2*i+1)*cos(real(i)*x)
	   funct(j,2)=funct(j,2)+pout(2*i+2)*sin(real(i)*x)
	   !f(x)=0.5*ao+sum(aj*cos(j*x)+bj*sin(j*x))  sum from j=1 to j=7
	 end do 
   end do 
	
   
    !writing the fourier series to a .dat file
   open(unit=11, file='function.dat')
     write(11,'(2(f10.6,1x))') funct
   close (11)

end program fftw2

 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 