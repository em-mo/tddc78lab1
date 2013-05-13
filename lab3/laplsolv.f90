program laplsolv
!-----------------------------------------------------------------------
! Serial program for solving the heat conduction problem 
! on a square using the Jacobi method. 
! Written by Fredrik Berntsson (frber@math.liu.se) March 2003
! Modified by Berkant Savas (besav@math.liu.se) April 2006
!-----------------------------------------------------------------------
  double precision omp_get_wtime
  integer omp_get_num_threads
  integer omp_get_thread_num
  integer, parameter                  :: n=1000, maxiter=1000
  double precision,parameter          :: tol=1.0E-3
  double precision,dimension(0:n+1,0:n+1) :: T
  double precision,dimension(n)       :: tmp1
  double precision,dimension(0:n+1)   :: tmp2
  double precision                    :: error,x
  real                                :: t1,t0
  integer                             :: i,j,k, iters
  character(len=20)                   :: str
  
  ! Set boundary conditions and initial values for the unknowns
  T=0.0D0
  T(0:n+1 , 0)     = 1.0D0
  T(0:n+1 , n+1)   = 1.0D0
  T(n+1   , 0:n+1) = 2.0D0
  

  ! Solve the linear system of equations using the Jacobi method
  t0 = omp_get_wtime()

  !$omp parallel private(i, j, k, tmp1, error, tmp2) shared(iter)

  numthreads = omp_get_num_threads()
  threadId = omp_get_thread_num()

  numcolumns = n/numthreads
  restCols = mod(numcolumns, numthreads)

  if (threadId < restcols) then
     numcolumns = numcolumns + 1
  end if

  startColumnIndex = numColumns * thread_id + min(restCols, threadId)
  lastColumnIndex = startColumnIndex + numColumns

  do k=1,maxiter
     tmp1=T(1:n,0)
     error=0.0D0
     
     !$omp do     
     do j=1,n
        
                
        T(1:n,j)=(T(0:n-1,j)+T(2:n+1,j)+T(1:n,j+1)+tmp1)/4.0D0
     
        !$omp critical
        error=max(error,maxval(abs(tmp2(1:n)-T(1:n,j))))
        !$omp end critical

        tmp1=tmp2(1:n)

     end do
     !$omp end do
     if (error<tol) then
        exit
     end if
     
  end do
  iter = k
  !$omp end parallel

  t1 = omp_get_wtime()

  write(unit=*,fmt=*) 'Time:',t1-t0,'Number of Iterations:',iter
  write(unit=*,fmt=*) 'Temperature of element T(1,1)  =',T(1,1)

  ! Uncomment the next part if you want to write the whole solution
  ! to a file. Useful for plotting. 
  
  !open(unit=7,action='write',file='result.dat',status='unknown')
  !write(unit=str,fmt='(a,i6,a)') '(',N,'F10.6)'
  !do i=0,n+1
  !   write (unit=7,fmt=str) T(i,0:n+1)  
  !end do
  !close(unit=7)
  
end program laplsolv
