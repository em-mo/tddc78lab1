program laplsolv
!-----------------------------------------------------------------------
! Serial program for solving the heat conduction problem 
! on a square using the Jacobi method. 
! Written by Fredrik Berntsson (frber@math.liu.se) March 2003
! Modified by Berkant Savas (besav@math.liu.se) April 2006
!-----------------------------------------------------------------------
  double precision omp_get_wtime
  integer omp_get_num_threads,omp_get_thread_num,numthreads,threadId, omp_get_max_threads
  integer, parameter                  :: n=100, maxiter=1000
  double precision,parameter          :: tol=1.0E-3
  double precision,dimension(0:n+1,0:n+1) :: T
  double precision,dimension(n)       :: tmp1, tmp2, startTmp, endTmp
  double precision                    :: error
  double precision                    :: t1,t0
  integer                             :: j,k,iters,startCol, endCol, numColumns
  character(len=20)                   :: str
  
  ! Set boundary conditions and initial values for the unknowns
  write(*,*) "Size", n, " threads ", omp_get_max_threads()

  T=0.0D0
  T(0:n+1 , 0)     = 1.0D0
  T(0:n+1 , n+1)   = 1.0D0
  T(n+1   , 0:n+1) = 2.0D0
  
  iter = 0;

  ! Solve the linear system of equations using the Jacobi method
  t0 = omp_get_wtime()

  !$omp parallel private(j, k, tmp1, tmp2, startCol, endCol, numthreads, threadId, numColumns, startTmp, endTmp) shared(error)

  numthreads = omp_get_num_threads()
  threadId = omp_get_thread_num()

  ! Calculate which columns that are needed to copy for each thread
  ! The allocation of columns follows the same allocation as omp do
  if (mod(n, numthreads) == 0) then
    numColumns = n / numthreads
    startCol = numColumns * threadId
    endCol = startCol + numColumns + 1
  else
    numColumns = ceiling(dble(n) / numthreads)
    if ((threadId+1) * numColumns <= n) then
       startCol = numColumns * threadId
       endCol = numColumns * (threadId + 1) + 1
    else if (threadId * numColumns < n .and. (threadId + 1) * numColumns > n) then
       startCol = numColumns * threadId
       endCol = n + 1
    else
       startCol = 0
       endCol = 0
    end if
  end if

  do k=1,maxiter
    tmp1=T(1:n,startCol)
    endTmp = T(1:n, endCol)
    
    error=0.0D0
 
    !$omp barrier
    !$omp do reduction(max:error)
    do j=1,n
      tmp2=T(1:n,j)

      ! Border case
      if (j == endCol - 1) then
        T(1:n,j)=(T(0:n-1,j)+T(2:n+1,j)+endTmp+tmp1)/4.0D0
      else        
        T(1:n,j)=(T(0:n-1,j)+T(2:n+1,j)+T(1:n,j+1)+tmp1)/4.0D0
      end if
      
      error=max(error,maxval(abs(tmp2-T(1:n,j))))
      
      numthreads = j
      tmp1=tmp2
    end do
    !$omp end do

    if (error<tol) then
      exit
    end if
    !$omp barrier
  end do

  ! To get the displaying of iterations correct
  !$omp critical
  if (iter < k) then
    iter = k
  end if
  !$omp end critical
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
