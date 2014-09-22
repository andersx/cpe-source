c     SUBROUTINE EGLCAO
c     =================
c
c     Copyright 1997 by Peter Blaudeck, Dirk Porezag, Michael Haugk,
c                       Joachim Elsner
c
c *********************************************************************
c
c     PROGRAM CHARACTERISTICS
c     -----------------------
c
c eglcao calculates energy and gradient for dylcao as shown by Seifert. 
c The determination of the occupation numbers has been changed to be
c also valid for metallic systems.
c
c PARAMETERS:
c n3      i  3 * number of atoms
c telec   r  electronic temperature (for Fermi distribution)
c x       r  coordinates (n3)
c scftol  r  convergence criterion for scf (if negative, no scf)
c atomic  l  should there be an atomic energy calculation
c e       r  total energy
c eel     r  electronic energy
c gr      r  gradient (n3)
c miter   i  number of scf-iterations performed
c evector l  prints out the eigenvectors
c lgrad   l  should gradients be calculated?
c
c *********************************************************************
c
       subroutine eglcao(n3,telec,x,scftol,atomic,e,eel,gr,
     &    miter,evector,qmat,ql,qlup,qldown,lgrad,lexc,nexc,nrem,
     &    docm3,qcm3,xm,excited,docpe)
c      =================
c
c maxima definitions
c
       implicit REAL*8 (A-H,O-Z)
       include 'maxima.inc'
c
       logical doscf, atomic, evector,constr,lgrad,lexc,excited
       integer izp(NNDIM), ind(NNDIM+1),lmax(MAXTYP)
       integer indj,indk,indj1,indk1,nstart,nexc,nrem
       integer conat(NNDIM+1),cai
       real*8 x(3,*), gr(3,*), hgrad(3,NNDIM),sgrad(3,NNDIM), scftol
       real*8 eel, ecoul, e, efermi,convec(3,NNDIM)
       real*8 au(LDIM,LDIM), bu(LDIM,LDIM), nel
       real*8 a(MDIM,MDIM), b(MDIM,MDIM), spro
c      real*8 h(34*MDIM), ev(MDIM), occ(MDIM)
c       modification for sgi-diagonalization
       real*8 ev(MDIM), occ(MDIM)
       real*8 qmulli(MDIM), qmold(NNDIM), qmat(NNDIM)
       real*8 qdiff(NNDIM),dipol(3),dipabs
       real*8 hamil(MDIM,MDIM),overl(MDIM,MDIM),shift(NNDIM,3)
       real*8 shift3(NNDIM,3),shift3A(NNDIM,3),shift3B(NNDIM)
       real*8 afudge(NNDIM),afuder(NNDIM),efudge
c       real*8 kspr
c      ME constraints like in CHARMM
c       integer lcc,jcc(nndim),atcc1(nndim),atcc2(nndim)
c       real*8 rccval,kccval,rccinc
c       common /cconstr/ lcc,jcc,atcc1,atcc2,rccval,kccval,rccinc
cccccccccccccccc  additional potential
c        real*8 ecorr
c        integer nco,nc(NNDIM,2)
cccccccccccccc CHARMM+ external charges
c      real*8 xE(3,4*NNDIM),ZE(4*NNDIM),shiftE(NNDIM)
       real*8 xE(3,EXTDIM),ZE(EXTDIM),shiftE(NNDIM)
       integer nE
       character*2 EXT
       character*1 sccmode
       logical fudge,dummyl
       real*8 a0fudge,bfudge,cfudge
       logical xhgamma,izpxh(MAXTYP)
       real*8 kl1,gammamat(NNDIM,NNDIM,3,3),gammader(NNDIM,NNDIM,3,3)
c      MG:120412 for collinear spin-formalism (JPCA2007,111,5622) 
       real*8  espin(MAXTYP)
       real*8  nelup,neldown,wspin(MAXTYP,3,3)
       logical spin
       integer mi,mj,mu,nu,ier,ier2,lldim,indl(NNDIM+1)
       real*8  evdown(MDIM),adown(MDIM,MDIM),spinshift(NNDIM,3)
       real*8  occdown(MDIM),qmulliup(MDIM),qmullidown(MDIM),phelp
       real*8  qmatdown(NNDIM),qmatup(NNDIM),ql(3*NNDIM),qlup(3*NNDIM)
       real*8  qldown(3*NNDIM),qlold(3*NNDIM),qlupold(3*NNDIM)
       real*8  qldownold(3*NNDIM),efermiup,efermidown
       real*8  dipolup(3),dipoldown(3)
c      MG:120504 l-shell dependent hubbard
       logical ldep
ccccccccccccccc
       external skhpar,skspar
c     DISPERSION
       logical dispers
       real*8 Edis
       common /disper/ Edis, dispers
c     DISPERSION END
cccccccccccccc CHARMM
       common /extchr/ xE, ZE, nE, nEiter ,icycle, EXT
cccccccccc
cccc additional potential
       common /correct/  r0yz , r0x ,kspr , nco, nc
cccccccccccc
       common /muell/ a,b
c      only for not taking space from stack (FORTRAN problem)
       common /concgr/ convec, constr, conat
c      constrains
       common /lmax/ lmax
c              s : 1, p : 2, u.s.w. for all atom types
       common /izp/ izp, nel, nbeweg
c              atom types, number of electrons, number of movable atoms
       common /spltab/ coeff,xr,efkt,cutoff,numint
c              spline data for repulsiv force

       include 'commonbox.inc'

       common /machine/ racc
c              machine accuracy
       common /hxmod/   kl1,xhgamma,izpxh
       common /mcharge/ qzero(MAXTYP,4), uhubb(MAXTYP,3), ldep
       common /uhq/ uhder(MAXTYP,3)
c
       common /scchelp/ sccmode

       common /fudgev/ fudge,dummyl,a0fudge,bfudge,cfudge
c      MG:120412 for collinear spin-formalism (JPCA2007,111,5622) 
       common /spin/   espin,nelup,neldown,wspin,spin

c      K.W. 2013-06-20 additional variables for CM3 calculation
       logical docm3
       real*8 qcm3(NNDIM),cm3dipol(3),cm3dipabs,order(NNDIM),qcm3tot
       integer cbnd,number1(NNDIM),number2(NNDIM)
       real*8 bord(nndim,nndim),xm(MAXTYP)
c      K.W. 2014-04-14 CPE calculation
       logical docpe
       real*8  E_cpe,C(3*NNDIM),shiftCPE(NNDIM),cpedipol(3),cpedipabs
       real*8  grad_cpe(3,NNDIM), CPEPar(NNDIM,4)
c
c setup
c
       ier=0
       ier2=0

c no gradients in mode 6
       if (atomic) lgrad= .false.
       if (scftol .lt. 0.0d0) then
         doscf = .false. 
       else
         doscf = .true.
       endif
       nn = n3 / 3
       
c calculation of indices for matrices H and S 
c
       ind(1) = 0
       do j = 1,nn 
         izpj = izp(j)
         ind(j+1) = ind(j)+lmax(izpj)**2
       end do
c
c actual dimension of matrix
       ndim = ind(nn+1)      
       if (ndim .gt. MDIM) then
         write(*,*) ' eglcao: ndim > ', MDIM
         stop
       endif

c maximum number of l-index
       lldim = 0
       do i=1,nn
         lldim=lldim+lmax(izp(i))
       enddo
c index for l-dependent charges
       indl(1)=0
       do i=1,nn
         indl(i+1) = indl(i)+lmax(izp(i))
       enddo
       
c setup of charge-independent part of H and S
c
       do j = 1,nn 
         do k = 1,j 
         
           call slkmatrices(j,k,x,au,bu)
           
           do n = 1,ind(k+1)-ind(k) 
             do m = 1,ind(j+1)-ind(j)
               hamil(ind(j)+m,ind(k)+n) = au(m,n)
               hamil(ind(k)+n,ind(j)+m) = au(m,n)       
               overl(ind(j)+m,ind(k)+n) = bu(m,n)
               overl(ind(k)+n,ind(j)+m) = bu(m,n)
             end do
           end do
         end do
       end do
c for external charges : extra shift shiftE 
      if ((EXT.eq.'CH') .or. (EXT.eq.'EF')) then
       call externalshift(nn,x,izp,shiftE)
      endif
c setup for SCF cycle
c
       if (doscf ) then
        if (atomic) then
         open(70,file='scfsum',form='formatted',status='unknown') 
        endif
c         close(70, status='delete')
         almix = 0.2d0
         maxiter = 70
         eelold = 1.0d10
       else
          maxiter = 1
       endif

c initialize oldcharges
         if (doscf) then
           if(ldep) then
             do i=1,lldim
               qlold(i) = ql(i)
             enddo
           else
             do i = 1,nn 
               qmold(i) = qmat(i)
             end do
           endif
           if (spin) then
             do i = 1,lldim
               qlupold(i)   = qlup(i)
               qldownold(i) = qldown(i)
             enddo
           endif
         endif

c precalculate gammamatrix and derivative for all atom-pairings
       if (doscf) then
         call get_gammamat(nn,x,izp,uhubb,uhder,kl1,boxsiz,period,
     &                     izpxh,ldep,lmax,gammamat,gammader)
       endif
       
c  *************************
c  SCF cycle starts here
c  *************************
       do niter = 1,maxiter 
         miter = niter
c        if ((niter.ne.1).and.(niter.eq.maxiter)) then
c          write(*,*) "ERROR in eglcao: SCC not converged!",
c    &       " maxiter=",maxiter,". Exit program"
c          stop
c        endif

c save old charges --> MG: not necessary because broyden takes care of that
c         if (doscf) then
c           do i = 1,nn 
c             qmold(i) = qmat(i)
c           end do
c         endif
       
c charge-independent part of H and S
c
       do j = 1,nn 
         indj  = ind(j)
         indj1 = ind(j+1)
         do k = 1,j 
           indk  = ind(k)
           indk1 = ind(k+1)
           
           do n = 1,indk1-indk 
             do m = 1,indj1-indj 
               a(indj+m,indk+n) = hamil(indj+m,indk+n)
               b(indj+m,indk+n) = overl(indj+m,indk+n)
             end do
           end do
         end do
       end do
       
c add charge dependent terms (Hubbard, etc. )
c start by defining qdiff and shift 
c
      if (doscf) then
           call HAMILSHIFT(nn,qmat,qzero,izp,qdiff,gammamat,gammader,
     &               sccmode,shift,shift3,shift3A,shift3B,ldep,ql,
     &               spin,qlup,qldown,wspin,indl,lmax,spinshift)

ccc add external charges! CHARMM
       if ((EXT.eq.'CH').or. (EXT.eq.'EF'))then 
         if (ldep) then
           do i=1,nn
           do li=1,lmax(izp(i))
             shift(i,li) = shift(i,li) - shiftE(i)
           enddo
           enddo
         else
           do i = 1,nn
             shift(i,1) = shift(i,1) - shiftE(i)
             ! minus, weil Z positiv sind!, positive DeltaQ meinen aber negative Ladung
           end do
         endif
       endif 

c K.W. 2014-04-14 SCF part of CPE calculation, 
c c/f Kaminski et al., JPC A 116 (2012), 9131. 
      if (docpe) then
        if (spin) then 
          write (*,*) 'Sorry, CPE implemented for closed shell
     &        for closed shell systems only at the moment'
          stop
        else if (ldep) then
          write (*,*) 'Sorry, CPE does not work with l-dependent
     &        hubbards at the moment' 
          stop
        else 
          call cpe_scf(nn,qmat,qzero,izp,x,xm,E_cpe,shiftCPE,
     &         C,dipol,cpedipol,cpedipabs,uhder,
     &         uhubb,CPEPar)


c Update the Hamiltonian shift with CPE potential
          do i = 1,nn
            shift(i,1) = shift(i,1) - shiftCPE(i)
          end do
        end if
      end if

c
c update hamiltonian matrix (shift3 =0 if sccmode.ne.3)
c
      if (ldep) then
        do i = 1,nn 
        do li = 1,lmax(izp(i))
        do mi = 1,2*li-1
           mu = ind(i)+(li-1)**2+mi
           do j = 1,i
           do lj = 1,lmax(izp(j))
           do mj = 1,2*lj-1
              nu = ind(j)+(lj-1)**2+mj
              a(mu,nu)     = a(mu,nu)+0.5d0*overl(mu,nu)
     &        *(shift(i,li)+shift(j,lj)+1.0d0/3.0d0*(shift3(i,li)+
     &          shift3A(i,li)+shift3(j,lj)+shift3A(j,lj)+
     &          shift3B(i)+shift3B(j)))
           end do
           end do
           end do
        end do
        end do
        end do
      else
        do i = 1,nn 
         do li = 1,lmax(izp(i))**2 
          do j = 1,i
           do lj = 1,lmax(izp(j))**2 
              a(ind(i)+li,ind(j)+lj) = a(ind(i)+li,ind(j)+lj) 
     &        +0.5d0*overl(ind(i)+li,ind(j)+lj)
     &        *(shift(i,1)+shift(j,1)+1.0d0/3.0d0*(2.0d0*shift3(i,1)+
     &          shift3A(i,1)+2.0d0*shift3(j,1)+shift3A(j,1)))
           end do
          end do
         end do
        end do
      endif

C Add high order fudge term (qdiff is updated in HAMILSHIFT)
c afuder(i) = d afudge(i) / d q(i)
c not for ldep as excluded in dylcao.r
       if (fudge) then
         do i = 1, nn 
           afudge(i) = a0fudge*dexp(-bfudge*(qdiff(i)-cfudge)**2)
           afuder(i) = -2.0d0*bfudge*(qdiff(i)-cfudge)*afudge(i)
         enddo    
c        update Hamiltonian
         do i = 1,nn
           do li= 1,lmax(izp(i))**2
             do j= 1,i
               do lj= 1,lmax(izp(j))**2
                 a(ind(i)+li,ind(j)+lj)=a(ind(i)+li,ind(j)+lj)
     &           +1.0d0/12.0d0*overl(ind(i)+li,ind(j)+lj)*(
     &           qdiff(i)**2*(3.0d0*afudge(i)+qdiff(i)*afuder(i))+
     &           qdiff(j)**2*(3.0d0*afudge(j)+qdiff(j)*afuder(j)) )
               end do  
             end do  
           end do  
         end do   
       endif   
        
      endif
!     end doscf      

c MG:120412 
      if(spin) then
        ! add spin contribution for spinup
        ! add spin contribution for spindown (additional matrix)
        do i = 1,nn 
        do li = 1,lmax(izp(i))
        do mi = 1,2*li-1
           mu = ind(i)+(li-1)**2+mi
           do j = 1,i
           do lj = 1,lmax(izp(j))
           do mj = 1,2*lj-1
              nu = ind(j)+(lj-1)**2+mj
              adown(mu,nu) = a(mu,nu)-0.5d0*overl(mu,nu)
     &          *(spinshift(i,li)+spinshift(j,lj))
              a(mu,nu)     = a(mu,nu)+0.5d0*overl(mu,nu)
     &          *(spinshift(i,li)+spinshift(j,lj))
           end do
           end do
           end do
        end do
        end do
        end do
      endif
         
c solution of eigenvalue problem 
c
         call ewevge(MDIM,MDIM,ndim,a,b,ev,h,1,-1,ier)
c MG: spinpolarized case: solve also for down-spin electrons
         if(spin) then
             do i=1,ndim
               do j=1,i
                 b(j,i) = 0.0d0
                 b(i,j) = overl(i,j)
               enddo
             enddo
             call ewevge(MDIM,MDIM,ndim,adown,b,evdown,h,1,-1,ier2)
         endif
       
         if ((ier.ne.0).or.(ier2.ne.0)) then 
           print *,' ewevge: ier =',ier,'niter=',niter
           stop
         endif

c Machine accuracy 
c
         dacc = 4.0d0*racc
c
c Calculate occupation (occ) and fermi energy (efermi)
c using fermi-distribution
c
         if (spin) then
           call FERMI(nelup,telec,ndim,ev,dacc,occ,efermiup,1.0d0)
           call FERMI(neldown,telec,ndim,evdown,dacc,occdown,
     &                                             efermidown,1.0d0)
         else if (excited) then
c          K.W. 2014-05-02 
c          return occupation numbers for a fixed HOMO/LUMO occupation
c          cf. fixed_fermi.f90
           call fixed_fermi(nel,ndim,occ,efermi,2.0d0)
         else
           call FERMI(nel,telec,ndim,ev,dacc,occ,efermi,2.0d0)
         endif

c sum of occupied eigenvalues
c
         eel = 0.0d0
         do i = 1,ndim 
           if (occ(i) .lt. dacc) goto 75
           eel = eel + occ(i)*ev(i)
         end do
75       continue
         if (spin) then
           do i = 1,ndim
             if (occdown(i).lt.dacc) goto 76
             eel = eel + occdown(i)*evdown(i)
           enddo
76         continue
         endif
c
c Lowest unoccupied level
c
c         nstart = i
       
c
c determine mulliken charges , the charge of the whole system
c and the mulliken 
c
         if(doscf) then
           ! spinpolarized - be careful: dipol will not be calculated!
           if (spin) then
             call MULLIKEN(nn,x,izp,qmatup,qzero,qmulliup,qlup,
     &                     qtotup,ndim,dacc,occ,a,
     &                     overl,lmax,ind,dipolup,dipabs)
             call MULLIKEN(nn,x,izp,qmatdown,qzero,qmullidown,qldown,
     &                     qtotdown,ndim,dacc,occdown,adown,
     &                     overl,lmax,ind,dipoldown,dipabs)
             qtot=qtotup+qtotdown
             do i=1,nn
               qmat(i) = qmatup(i)+qmatdown(i)
             enddo
             do i=1,lldim
               ql(i) = qlup(i)+qldown(i)
             enddo
           ! spinunpolarized
           else
             call MULLIKEN(nn,x,izp,qmat,qzero,qmulli,ql,qtot,
     &              ndim,dacc,occ,a,overl,lmax,ind,dipol,dipabs)
           endif
           if (dabs(qtot-nel).gt.1.0d-06) then
             write(*,*) "Warning: actual number of electrons is ",qtot
           endif 
         endif

       
c output spektra, occupation numbers, mulliken charges 
c
c        if (atomic) then
c         call outspec(nn,nbeweg,ndim,ind,lmax,izp,ev,occ,efermi,
c     &                qmat,qmulli,qtot,dipol,dipabs,
c     &                spin,evdown,occdown,qmatup,qmatdown,qmulliup,
c     &                qmullidown,qtotup,qtotdown,efermiup,efermidown)
c        endif 
c         if (evector) then
c           call outeigenvectors(a,ev,occ,ind,nn)
c         endif
c        do i=1,nn
c          ind1 = ind(i)+1; ind2 = ind1+lmax(izp(i))**2-1
c          write (*,'(1x,A6,i4,11(1x,f12.6))') 
c     &        "up",i,qmatup(i),(qmulliup(j), j= ind1,ind2)
c          write (*,'(1x,A6,i4,11(1x,f12.6))') 
c     &        "down",i,qmatdown(i),(qmullidown(j), j= ind1,ind2)
c        enddo
        
       
c complete calculation of electronic energy:
c charge-dependent energy contribution
c warning: this will only lead to the right result if convergence
c has been reached
c note: shift3=0 if sccmode.ne.3
         if (doscf) then
           ecoul = 0.0d0
           ecoul3=0.0d0
           eext = 0.0d0
           if (ldep) then
             do i=1,nn
               do li=1,lmax(izp(i))
                 indli=indl(i)+li
                 ecoul =ecoul +shift(i,li)*(ql(indli)+qzero(izp(i),li))
                 ecoul3=ecoul3+shift3(i,li)*qzero(izp(i),li)
     &                        +(shift3A(i,li)+shift3B(i))*ql(indli)
               enddo
             enddo
           else
             do i = 1,nn
              ecoul = ecoul + shift(i,1)*(qmat(i)+qzero(izp(i),4))
              ecoul3= ecoul3+ shift3(i,1)*(qmat(i)+qzero(izp(i),4))
     &                      + shift3A(i,1)*qmat(i)
             end do
           endif
ccccccc CHARMM
           if((EXT.eq.'CH').or. (EXT.eq.'EF'))then
           do i = 1,nn
            eext = eext + shiftE(i)*(qzero(izp(i),4)-qmat(i))
           end do
           endif

c        highorder fudge term
           efudge = 0.0d0
           if (fudge) then
             do i=1,nn
             efudge = efudge + afudge(i)*qdiff(i)**3
     &         -qmat(i)*qdiff(i)**2*(3.0d0*afudge(i)+qdiff(i)*afuder(i))
             enddo
             efudge = efudge/6.0d0
           endif

c        spinpolarization term
c        etot=sum_i n_i epsilon_i - ..... - 1/2*sum_a sum_l sum_l' p_al p_al' W_all'  = ... -1/2*sum_a sum_l p_al spinshift(a,l)
           espinpol=0.0d0
           if (spin) then
              do i=1,nn
                do li=1,lmax(izp(i))
                  indli=indl(i)+li
                  espinpol = espinpol + spinshift(i,li)*
     &                      ( qlup(indli) - qldown(indli) )
                 enddo
               enddo
           endif

ccccccc 
           eel = eel-0.5d0*ecoul -ecoul3/3.0d0 + 0.5d0*eext + efudge 
     &           -0.5d0*espinpol

         endif
c remark: eel contains shiftE aready via ev,
c shift also contains -shiftE, i.e. ecoul also
c contains contributions from EXT
       
c print energy(deleted), check convergence, call broyden mixing
c
         if (doscf ) then
           if (dabs(eel-eelold) .lt. scftol) goto 70
           eelold = eel
c
c   broyden mixing
c
           call mixer(niter,almix,nn,qmold,qmat,spin,ndim,lldim,ldep,
     &          lmax,izp,indl,qlold,ql,qlupold,qldownold,qlup,qldown)
         end if
       end do
70     continue
       
c end of SCF

c  update qdiff and shifts
       if (doscf) then
           call HAMILSHIFT(nn,qmat,qzero,izp,qdiff,gammamat,gammader,
     &               sccmode,shift,shift3,shift3A,shift3B,ldep,ql,
     &               spin,qlup,qldown,wspin,indl,lmax,spinshift)
       endif
       
c K.W. 2014-01-31 calculation of CM3 charges
      if (docm3) then
        if (spin) then
          write (*,*) "Sorry, CM3 charge calculation implemented
     &                for closed shell systems only at the moment"
          stop 
        else
          call cm3charges(nn,izp,qcm3,qmat,ndim,x,qzero,
     &         cm3dipol,cm3dipabs,ind,occ,dacc,a,overl,cbnd,
     &         number1,number2,order,xm,bord,qcm3tot,lmax)
        endif
      endif

c K.W. 2014-05-07 calculation of CPE gradient
      if (docpe) then
        call cpe_grad(nn,izp,x,qmat,qzero,CPEPar,C,grad_cpe,uhubb,uhder)
      end if

c output spektra, occupation numbers, mulliken charges 
c
c        if (atomic) then
         call outspec(nn,nbeweg,ndim,ind,lmax,izp,ev,occ,efermi,
     &                qmat,qmulli,qtot,dipol,dipabs,
     &                spin,evdown,occdown,qmatup,qmatdown,qmulliup,
     &                qmullidown,qtotup,qtotdown,efermiup,efermidown,
     &                qzero,docm3,qcm3,qcm3tot,cm3dipol,cm3dipabs,cbnd,
     &                order,number1,number2,bord,docpe,E_cpe,C,
     &                cpedipabs,cpedipol,grad_cpe,shiftCPE)

c        endif 
         if (evector) then
           if (spin) then
             write(*,*) "Warning: print eigenvectors for spin 
     &                   currently not implemented."
           else
             call outeigenvectors(a,ev,occ,ind,nn)
           endif
         endif
        
ccc add EXTERNAL CHARGES CHARMM
          if ((EXT.eq.'CH').or. (EXT.eq.'EF'))then
            if (ldep) then
              do i=1,nn
              do li=1,lmax(izp(i))
                shift(i,li) = shift(i,li) - shiftE(i)
              enddo
              enddo
            else
              do i = 1,nn
                shift(i,1) = shift(i,1) - shiftE(i)
                ! minus, weil Z positiv sind!, positive DeltaQ meinen aber negative Ladung  
              enddo
            endif
          endif
c CPE shifts (? again)
          if (docpe) then
              do i = 1,nn
                shift(i,1) = shift(i,1) - shiftCPE(i)
              end do
          end if

       if (atomic) then       
c
c calculation of atomic energies
c 
         if (sccmode.eq.'3') then 
           write(*,*) '3rd order currently not implemented for "atomic"
     &           calculation'
           stop
         endif
         if (spin) then 
           write(*,*) 'spin-polarization-formalism not implemented for
     &           "atomic" calculation'
         endif
         if (ldep) then
           write(*,*) 'l-dependent Hubbard formalism not implemented for
     &            "atomic" calculation'
         endif
         call atomenerg(nn,ndim,ind,izp,shift,shift3,qdiff,
     &      x,hamil,a,occ,doscf,period,boxsiz,xinvbox,nlat)
       endif
       
c start with repulsive energy
       erep   = 0.0d0
       
c
c Calculate repulsive energy
c
       call repulsive(nn,izp,period,x,boxsiz,xinvbox,nlat,erep)
     
c total energy 
c
       e = eel + erep 

c excitation energies via lin. response
c
c       if(lexc) then       
c         call lrespo(nn,ind,ndim,a,ev,overl,x,boxsiz,
c     &                         period,nexc,nrem,e)
c       endif

       
c
c start with gradient
        
       do j = 1,nn 
         do i = 1,3 
           gr(i,j) = 0.0d0
           hgrad(i,j) = 0.0d0
         end do
       end do

c  put high order fudge term on shift for forces
      if (fudge) then
        do i=1,nn
         shift(i,1) = shift(i,1) + 1.0d0/6.0d0*
     &     qdiff(i)**2*(3.0d0*afudge(i)+qdiff(i)*afuder(i))
        enddo
      endif


c DISERPSION
       if (dispers) then
         call  dispersion_egr(nn,x,hgrad,boxsiz,nlat,period)
        do j=1,nn
         do i=1,3
           gr(i,j) = gr(i,j) + hgrad(i,j)
         end do
        end do
        e = e + Edis
       endif
c END DISPERSION       
c BEGIN EXTERNAL CHARGES CHARMM
      if((EXT.eq.'CH').or. (EXT.eq.'EF'))then
       call externalchgrad(nn,x,izp,hgrad,qmat)
        do j=1,nn
         do i=1,3
           gr(i,j) = gr(i,j) + hgrad(i,j)
         end do
        end do
      endif
c END EXTERNAL CHRARGES
       if(lgrad) then
c

c gradient 
c
       call usualgrd(nn,nbeweg,ndim,izp,lmax,ind,period,doscf,
     &               x,ev,a,b,occ,shift,shift3,shift3A,shift3B,dacc,
     &               boxsiz,xinvbox,hgrad,ldep,spin,1.0d0,spinshift)
       do j=1,nn 
        do i=1,3
         gr(i,j) = gr(i,j) + hgrad(i,j)
        end do
       end do

       if (spin) then
           call usualgrd(nn,nbeweg,ndim,izp,lmax,ind,period,doscf,
     &       x,evdown,adown,b,occdown,shift,shift3,shift3A,shift3B,dacc,
     &       boxsiz,xinvbox,hgrad,ldep,spin,-1.0d0,spinshift)
           do j=1,nn 
            do i=1,3
             gr(i,j) = gr(i,j) + hgrad(i,j)
            end do
           end do
       endif
                    
        
c here are the contributions due to gamma if in scf mode
       if (doscf) then
         call GAMMAGRAD(nn,nbeweg,x,izp,uhubb,uhder,boxsiz,period,
     &            izpxh,kl1,qdiff,sccmode,ldep,ql,qzero,lmax,indl,hgrad)
         do j=1,nn 
           do i=1,3 
             gr(i,j) = gr(i,j) + hgrad(i,j)
           end do
         end do
       endif

c Calculate gradient of the repulsive energy
c
       call repulsivegrd(nn,nbeweg,izp,period,x,boxsiz,
     &                   xinvbox,nlat,hgrad)

c write "dftb.e" which is used for program erepfit
        open (unit=99, file='dftb.e',status='replace')
        write (99,'(2E25.16)') eel,erep
        do i=1,nn
            write(99,'(4E25.16)') (-gr(j,i),j=1,3)
        end do
        close(99)

c add gradient of repulsive energy to total gradient
       do j=1,nn 
         do i=1,3 
           gr(i,j) = gr(i,j) + hgrad(i,j)
         end do
       end do

c K.W. 2014-05-02; CPE gradient calculation
      if (docpe) then
        do j = 1,nn
          do i = 1,3
            gr(i,j) = gr(i,j) + grad_cpe(i,j)
          end do
        end do
      endif

c
c constraints?
c
       if(constr) then
        do i=1,conat(1)
         cai=conat(i+1)
         spro=convec(1,i)*gr(1,cai) 
         spro=spro+convec(2,i)*gr(2,cai) 
         spro=spro+convec(3,i)*gr(3,cai) 
c         gr(1,cai)=gr(1,cai)-spro*convec(1,i)
c         gr(2,cai)=gr(2,cai)-spro*convec(2,i)
c         gr(3,cai)=gr(3,cai)-spro*convec(3,i)
c other constraint! fix x,y,z individually
         gr(1,cai)=gr(1,cai)*convec(1,i)
         gr(2,cai)=gr(2,cai)*convec(2,i)
         gr(3,cai)=gr(3,cai)*convec(3,i)
        end do
       endif

c End of gradient calculation , only done when not in mode 6
c
       endif

       end
