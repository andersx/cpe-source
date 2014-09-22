subroutine outvib(flag,nn,atnames,x,gr,xm,izp,e,qzero,qmat,docm3,qcm3)

include 'maxima.h'

logical      flag,docm3
integer      nn,izp(NNDIM)
character*70 atnames
real*8       x(3,NNDIM),gr(3*NNDIM),xm(MAXTYP),e
real*8       qzero(MAXTYP,4),qmat(NNDIM),qcm3(NNDIM)
integer      i,j
real*8       xhlp(3),dipol(3),cm3dipol(3)

# K.W. 2013-06-20 added CM3 support for spectral calculation
  # open files
  if (flag){
    open (1,file='POSITION.0',status='replace')
    open (2,file='FORCE.0',status='replace')   
    open (3,file='DIPOLE.0',status='replace')   
    open (4,file='MASSES-E0',status='replace')  
    if (docm3) {
      open (5,file='CM3DIPOLE.0',status='replace')
    }
  }
  else{
    open (1,file='POSITION',status='unknown',access='append')
    open (2,file='FORCE',status='unknown',access='append')
    open (3,file='DIPOLE',status='unknown',access='append')
    open (4,file='ENERGY',status='unknown')
    #open (4,file='ENERGY',status='unknown',access='append')
    if (docm3) {
      open (5,file='CM3DIPOLE',status='unknown',access='append')
    }
  }

  # write POSITION*
  write(1,'(I4,1X,A1)') nn,'C'
  write(1,'(A)') atnames
  do i=1,nn {
    do j=1,3 {
     xhlp(j)=x(j,i)*0.529177d0
    }
    write(1,'(I4,1X,I3,3(1X,F12.6))') i,izp(i),(xhlp(j),j=1,3)
  }
  close(1)

  # write FORCE*
  write(2,'(A)') '# INTERATOMIC FORCES AT ATOM N0. I'
  write(2,'(A)') '# ATOM N0.'
  write(2,'(A)') '# FORCE_X, FORCE_Y, FORCE_Z'
  write(2,'(A)') '#'
  write(2,'(A)') '#'
  do i = 1,nn {
    write (2,'(i5,3(1x,e19.10))') i,(-gr(3*i+j-3), j = 1,3)
  }
  write(2,'(A)') 'rms/d:   0.0000   0.0000'
  close (2)
 
  # calculation of dipol moment and CM3
  do i = 1,3 {
    dipol(i) = 0.0d0
    cm3dipol(i) = 0.0d0
    do j = 1,nn {
      dipol(i) = dipol(i)+(qzero(izp(j),4)-qmat(j))*x(i,j)
      if (docm3) {
        cm3dipol(i) = cm3dipol(i) + (qzero(izp(j),4)-qcm3(j)) * x(i,j)
      }
    }
    # conversion to debye
    dipol(i) = dipol(i)*2.541765d0
  }

  # write DIPOLE*
  write (3,'(3(1x,f20.12))') (dipol(i), i= 1,3)
  close(3)
  
  if (docm3) {
    write (5,'(3(1x,f20.12))') (cm3dipol(i), i=1,3)
    close(5)
  }

  # write MASSES-E0
  if (flag){
    do i = 1,nn {
      write (4,'(1x,f24.15)') xm(izp(i))
    }
    write(4,'(1x,f24.15)') e  
    close(4)
  }
  # write ENERGY
  else {
    write(4,'(1x,f24.15)') e
  }

END

