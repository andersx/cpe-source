
c     geometry optimization using bfgs (bfgs.f)

      subroutine bfgsgeom(n3,telec,coord,scftol,atomic,energy,eel,gr,
     & niter,evector,qmat,lgrad,lexc,nexc,nrem,geseatom,icycle,nbeweg, 
     & izp,xm,maxrun,fmax,name,namev,period,ntype,boxsiz,atnames,mode,
     & outfile,qzero,converge,newtries,ql,qlup,qldown,lmax,docm3,qcm3)
 
      implicit none    
      include 'maxima.inc'
      logical atomic,evector,lgrad,lexc,period,converge,docm3
      integer n3,niter,nexc,nrem,icycle,nbeweg,izp(*),maxrun,ntype,mode
      integer newtries,lmax(MAXTYP),lcount
      real*8  telec,coord(3,*),scftol,energy,eel,gr(3,*),qmat(*)
      real*8  qcm3(NNDIM),xm(MAXTYP)
      real*8  fmax,geseatom,boxsiz(3,3),qzero(MAXTYP,4)
      real*8  ql(*),qlup(*),qldown(*)
      character*65 name,namev,outfile
      character*70 atnames
 
      real*8  Edis,dxsmx,dxtmp,fsmx,ftmp,ranmar
      logical dispers,zerorun

      integer          nmax, mmax, lenwa
      parameter       (nmax  = 3*NNDIM, mmax = 5)
      parameter       (lenwa = 2*mmax*nmax +  4*nmax
     +                      + 11*mmax*mmax + 8*mmax)

c     nmax  is the dimension of the largest problem to be solved.
c     mmax  is the maximum number of limited memory corrections.
c     lenwa is the corresponding real workspace required.
 
c     Declare the variables needed by the code.
c     A description of all these variables is given at the end of the file
 
      character*60     task, csave
      logical          lsave(4)
      integer          n, m, iprint,
     +                 nbd(nmax), iwa(3*nmax), isave(44)
      double precision f, factr, pgtol, 
     +                 x(nmax), l(nmax), u(nmax), g(nmax), dsave(29), 
     +                 wa(lenwa),xold(nmax)

c     Declare a few additional variables for this sample problem.

      integer          i,j

      common /disper/ Edis, dispers

c     set number of bfgs-restarts to zero
      newtries=0
      zerorun=.true.

c     We wish to have output at every iteration.

      iprint = -1

c     We specify the tolerances in the stopping criteria.
c     we will have our own criteria - if (task(1:5) .eq. 'NEW_X')
      factr  = 0.0d0
      pgtol  = 0.0d0

c     We specify the dimension n of the sample problem and the number
c     m of limited memory corrections stored.  (n and m should not
c     exceed the limits nmax and mmax respectively.)
 
      n      = n3
      m      =  5
 
c     We now provide nbd which defines the bounds on the variables:
c     l   specifies the lower bounds,
c     u   specifies the upper bounds. 
 
c     We set no bounds (geometry constraints are handled within eglcao)

      do  i = 1, n
         nbd(i) = 0
      enddo

c     We now define the starting point.

      do i = 1,n3/3
        do j = 1,3
          x((i-1)*3+j) = coord(j,i)
        enddo
      enddo
      do i = 1,n3/3
        do j = 1,3
          xold((i-1)*3+j) = coord(j,i)
        enddo
      enddo

c     open files
      if (icycle.eq.1) then
        open(95,file='ENERGY.TMP',form='formatted',status='unknown')  
      else
        open(95,file='ENERGY.TMP',form='formatted',status='unknown'
     &       ,access='append')  
      endif

c     We start the iteration by initializing task.
 
      task = 'START'

c     ------- The beginning of the loop ----------
  111 continue
      
c     This is the call to the L-BFGS-B code.
 
      call setulb(n,m,x,l,u,nbd,f,g,factr,pgtol,wa,iwa,task,iprint,
     +            csave,lsave,isave,dsave)
 
      if (task(1:2) .eq. 'FG') then
 
c        The minimization routine has returned to request the
c        function f and gradient g values at the current x.

c        Compute function value f and gradient g for the sample problem.
  
         call eglcao(n3,telec,x,scftol,atomic,energy,eel,gr,niter,
     &               evector,qmat,ql,qlup,qldown,lgrad,lexc,nexc,nrem, 
     &               docm3,qcm3,xm)
         f = energy
         do i = 1,n3/3
           do j = 1,3
             g((i-1)*3+j)=gr(j,i)
           enddo
         enddo

c        output for first calculation at iteration 0
         if ((icycle.eq.1).AND.(zerorun)) then
           call outforces(nbeweg,izp,gr,xm)
           !calculate maximal force fsmx
           fsmx = 0.0d0
           do i = 1,nbeweg 
             ftmp = 0.0d0
             do j = 1,3 
               ftmp = ftmp + gr(j,i)**2
             enddo
             ftmp = dsqrt(ftmp)
             fsmx = max(fsmx,ftmp)
           enddo  
           ! write output
           write (*,13) icycle,0,niter,energy,eel,
     &        (energy-geseatom)*627.5095d0,Edis,dxsmx,fsmx
           write (95,13) icycle,0,niter,energy,eel,
     &        (energy-geseatom)*627.5095d0,Edis,dxsmx,fsmx
           write (name,'(a,i3.3)') outfile(:index(outfile,' ')-1),icycle
           call outcoord(name, namev, n3/3, period, ntype, izp, 
     &                   x, x, boxsiz, atnames,mode)
           if (fsmx.lt.fmax) then
             converge=.true.
           endif
           ! if only single-point calculation or converged, end bfgs
           if ((maxrun.eq.0).or.(converge)) then 
             goto 112
           endif
           zerorun=.false.
         endif

c        Go back to the minimization routine.
         goto 111

      elseif (task(1:5) .eq. 'NEW_X') then
 
c        The minimization routine has returned with a new iterate

c        calculate maximal displacement dxsmx
         dxsmx = 0.0d0
         do i = 1,nbeweg 
           dxtmp = 0.0d0
           do j = -2,0 
             dxtmp = dxtmp + (x(3*i+j)-xold(3*i+j))**2
             xold(3*i+j) = x(3*i+j)
           enddo
           dxtmp = dsqrt(dxtmp)
           dxsmx = max(dxsmx,dxtmp)
         enddo

c        calculate maximal force fsmx
         fsmx = 0.0d0
         do i = 1,nbeweg 
           ftmp = 0.0d0
           do j = 1,3 
             ftmp = ftmp + gr(j,i)**2
           enddo
           ftmp = dsqrt(ftmp)
           fsmx = max(fsmx,ftmp)
         enddo  

c        write some output
         call outforces(nbeweg,izp,gr,xm)

         write (*,13) icycle,isave(30),niter,energy,eel,
     &      (energy-geseatom)*627.5095d0,Edis,dxsmx,fsmx
         write (95,13) icycle,isave(30),niter,energy,eel,
     &      (energy-geseatom)*627.5095d0,Edis,dxsmx,fsmx

         write (name,'(a,i3.3)') outfile(:index(outfile,' ')-1), icycle
         call outcoord(name, namev, n3/3, period, ntype, izp, 
     &       x, xold, boxsiz, atnames,mode)
 
 
c        goto next iteration step if forcetolerance and maxrun is not reached
         if ((isave(30).lt.maxrun) .AND. (fsmx.ge.fmax)) then
             goto 111
         else
             if (fsmx.lt.fmax) converge=.true.
             do i=1,n3/3
               do j=1,3
                 coord(j,i)=x((i-1)*3+j)
               enddo
             enddo
             write(*,*) "bfgs relaxation: ",isave(34)," eglcao calls"
         endif
         
      else

c        Try again when task is neither FG nor NEW_X
         newtries=newtries+1
         if (newtries.eq.1) then
           lcount=0
           do i=1,n3/3
             qmat(i)=qzero(izp(i),4)
             do j=1,lmax(izp(i))
               lcount=lcount+1
               ql(lcount)     = qzero(izp(i),j)
               qlup(lcount)   = qzero(izp(i),j)/2.0d0
               qldown(lcount) = qzero(izp(i),j)/2.0d0
             enddo
           enddo
           write(*,*) "bfgs relaxation: ",isave(34)," eglcao calls"
           write(*,*) "bfgs error:",task
           write(*,*) "restart with qmat(i)=qzero(izp(i),4),i=1,nat"
           goto 111
         endif 
         ! if newtries==3 siehe dylcao.r (ierror)
         if (newtries.lt.3) then
           lcount=0
           do i=1,n3/3
             qmat(i)=qzero(izp(i),4)
             do j=1,lmax(izp(i))
               lcount=lcount+1
               ql(lcount)     = qzero(izp(i),j)
               qlup(lcount)   = qzero(izp(i),j)/2.0d0
               qldown(lcount) = qzero(izp(i),j)/2.0d0
             enddo
           enddo
           do i=1,n3
             x(1)=x(1)+0.01d0*(ranmar()-0.5d0)
           enddo
           write(*,*) "bfgs relaxation: ",isave(34)," eglcao calls"
           write(*,*) "bfgs error:",task
           write(*,*) "restart with displaced coordinates and",
     &                " qmat(i)=qzero(izp(i),4),i=1,nat"
           goto 111
         endif

c        We terminate execution when task is neither FG nor NEW_X.
c        We print the information contained in the string task
c        if the default output is not used and the execution is
c        not stopped intentionally by the user. 

         if (iprint .le. -1 .and. task(1:4) .ne. 'STOP') then
           write(*,*) "bfgs relaxation: ",isave(34)," eglcao calls"
           write(*,*) "BFGS geometry optimization ended abnormally"
           write(*,*) "Error message: ",task
         endif

      endif

c     ---------- The end of the loop -------------

112   continue
c     close open files
      close(95)

 13   format(i5,1x,i5,' / ',i2,3(1x,f14.6),1x,f14.6,2(1x,f14.6))

 
      end

c======================= The end of driver1 ============================

c      call setulb(n,m,x,l,u,nbd,f,g,factr,pgtol,wa,iwa,task,iprint,
c     +            csave,lsave,isave,dsave)
c     --------------------------------------------------------------
c             DESCRIPTION OF THE VARIABLES IN L-BFGS-B
c     --------------------------------------------------------------
c
c     n is an INTEGER variable that must be set by the user to the
c       number of variables.  It is not altered by the routine.
c
c     m is an INTEGER variable that must be set by the user to the
c       number of corrections used in the limited memory matrix.
c       It is not altered by the routine.  Values of m < 3  are
c       not recommended, and large values of m can result in excessive
c       computing time. The range  3 <= m <= 20 is recommended. 
c
c     x is a DOUBLE PRECISION array of length n.  On initial entry
c       it must be set by the user to the values of the initial
c       estimate of the solution vector.  Upon successful exit, it
c       contains the values of the variables at the best point
c       found (usually an approximate solution).
c
c     l is a DOUBLE PRECISION array of length n that must be set by
c       the user to the values of the lower bounds on the variables. If
c       the i-th variable has no lower bound, l(i) need not be defined.
c
c     u is a DOUBLE PRECISION array of length n that must be set by
c       the user to the values of the upper bounds on the variables. If
c       the i-th variable has no upper bound, u(i) need not be defined.
c
c     nbd is an INTEGER array of dimension n that must be set by the
c       user to the type of bounds imposed on the variables:
c       nbd(i)=0 if x(i) is unbounded,
c              1 if x(i) has only a lower bound,
c              2 if x(i) has both lower and upper bounds, 
c              3 if x(i) has only an upper bound.
c
c     f is a DOUBLE PRECISION variable.  If the routine setulb returns
c       with task(1:2)= 'FG', then f must be set by the user to
c       contain the value of the function at the point x.
c
c     g is a DOUBLE PRECISION array of length n.  If the routine setulb
c       returns with taskb(1:2)= 'FG', then g must be set by the user to
c       contain the components of the gradient at the point x.
c
c     factr is a DOUBLE PRECISION variable that must be set by the user.
c       It is a tolerance in the termination test for the algorithm.
c       The iteration will stop when
c
c        (f^k - f^{k+1})/max{|f^k|,|f^{k+1}|,1} <= factr*epsmch
c
c       where epsmch is the machine precision which is automatically
c       generated by the code. Typical values for factr on a computer
c       with 15 digits of accuracy in double precision are:
c       factr=1.d+12 for low accuracy;
c             1.d+7  for moderate accuracy; 
c             1.d+1  for extremely high accuracy.
c       The user can suppress this termination test by setting factr=0.
c
c     pgtol is a double precision variable.
c       On entry pgtol >= 0 is specified by the user.  The iteration
c         will stop when
c
c                 max{|proj g_i | i = 1, ..., n} <= pgtol
c
c         where pg_i is the ith component of the projected gradient.
c       The user can suppress this termination test by setting pgtol=0.
c
c     wa is a DOUBLE PRECISION  array of length 
c       (2mmax + 4)nmax + 11mmax^2 + 8mmax used as workspace.
c       This array must not be altered by the user.
c
c     iwa is an INTEGER  array of length 3nmax used as
c       workspace. This array must not be altered by the user.
c
c     task is a CHARACTER string of length 60.
c       On first entry, it must be set to 'START'.
c       On a return with task(1:2)='FG', the user must evaluate the
c         function f and gradient g at the returned value of x.
c       On a return with task(1:5)='NEW_X', an iteration of the
c         algorithm has concluded, and f and g contain f(x) and g(x)
c         respectively.  The user can decide whether to continue or stop
c         the iteration. 
c       When
c         task(1:4)='CONV', the termination test in L-BFGS-B has been 
c           satisfied;
c         task(1:4)='ABNO', the routine has terminated abnormally
c           without being able to satisfy the termination conditions,
c           x contains the best approximation found,
c           f and g contain f(x) and g(x) respectively;
c         task(1:5)='ERROR', the routine has detected an error in the
c           input parameters;
c       On exit with task = 'CONV', 'ABNO' or 'ERROR', the variable task
c         contains additional information that the user can print.
c       This array should not be altered unless the user wants to
c          stop the run for some reason.  See driver2 or driver3
c          for a detailed explanation on how to stop the run 
c          by assigning task(1:4)='STOP' in the driver.
c
c     iprint is an INTEGER variable that must be set by the user.
c       It controls the frequency and type of output generated:
c        iprint<0    no output is generated;
c        iprint=0    print only one line at the last iteration;
c        0<iprint<99 print also f and |proj g| every iprint iterations;
c        iprint=99   print details of every iteration except n-vectors;
c        iprint=100  print also the changes of active set and final x;
c        iprint>100  print details of every iteration including x and g;
c       When iprint > 0, the file iterate.dat will be created to
c                        summarize the iteration.
c
c     csave  is a CHARACTER working array of length 60.
c
c     lsave is a LOGICAL working array of dimension 4.
c       On exit with task = 'NEW_X', the following information is
c         available:
c       lsave(1) = .true.  the initial x did not satisfy the bounds;
c       lsave(2) = .true.  the problem contains bounds;
c       lsave(3) = .true.  each variable has upper and lower bounds.
c
c     isave is an INTEGER working array of dimension 44.
c       On exit with task = 'NEW_X', it contains information that
c       the user may want to access:
c         isave(30) = the current iteration number;
c         isave(34) = the total number of function and gradient
c                         evaluations;
c         isave(36) = the number of function value or gradient
c                                  evaluations in the current iteration;
c         isave(38) = the number of free variables in the current
c                         iteration;
c         isave(39) = the number of active constraints at the current
c                         iteration;
c
c       See the subroutine setulb.f for a description of other 
c       information contained in isave.
c
c     dsave is a DOUBLE PRECISION working array of dimension 29.
c       On exit with task = 'NEW_X', it contains information that
c       the user may want to access:
c         dsave(2) = the value of f at the previous iteration;
c         dsave(5) = the machine precision epsmch generated by the code;
c         dsave(13) = the infinity norm of the projected gradient;
c
c       See the subroutine setulb.f for a description of other 
c       information contained in dsave.
c
c     --------------------------------------------------------------
c           END OF THE DESCRIPTION OF THE VARIABLES IN L-BFGS-B
c     --------------------------------------------------------------
c
c     << An example of subroutine 'timer' for AIX Version 3.2 >>
c
c     subroutine timer(ttime)
c     double precision ttime
c     integer itemp, integer mclock
c
c     itemp = mclock()
c     ttime = dble(itemp)*1.0d-2
c     return
c     end
c-----------------------------------------------------------------------

