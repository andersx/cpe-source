
      function dist(x,a,c)
        implicit none
        include 'maxima.inc'
        real*8  dist,x(3,nndim),rx,ry,rz
        integer a,c

        rx   = x(1,c)-x(1,a)
        ry   = x(2,c)-x(2,a)
        rz   = x(3,c)-x(3,a)
        dist = dsqrt(rx*rx+ry*ry+rz*rz)

        return
      end


      function gamz(r,u)
        ! gamz   = \int\int \frac{\rho(r_a)}{|r_a-R_b|}\text{d}r_a\text{d}r_b
        ! \rho(x)= \frac{\tau^3}{8\pi}\text{exp}(-\tau|x-r_a|)
        ! r      = |r_a - R_b|
        ! u      = 16/5*\tau_a
        implicit none
        real*8 gamz,r,u,zero
        parameter(zero=1.0d-4)

        if (r .lt. zero) then
          gamz=1.6d0*u
          return
        else
          gamz = 1.0d0/r - (1.6d0*u+1.0d0/r)*dexp(-3.2d0*u*r)
          return
        endif
      end


      function gamnew(r,ua,ub)
        ! gamnew = \int\int \frac{\rho(r_a)\rho(r_b)}{|r_a-r_b|}\text{d}r_a\text{d}r_b
        ! \rho(x)= \frac{\tau^3}{8\pi}\text{exp}(-\tau|x-r_a|)
        ! r      = |r_a - r_b|
        ! ua     = 16/5*\tau_a
        implicit none
        real*8 gamnew,r,ua,ub,zero,a,b,gams
        parameter(zero=1.0d-4)
    
        a=3.2d0*ua
        b=3.2d0*ub 

        if(a+b .lt. zero) then
          gamnew = 0.0d0
          return
        elseif (r .lt. zero) then
          gamnew = 0.5d0*( a*b/(a+b)+a*a*b*b/((a+b)**3) )
          return
        else 
          gamnew = 1.0d0/r - gams(r,ua,ub)
          return
        endif
      end


      function gams(r,ua,ub)
        implicit none
        real*8 gams,r,ua,ub
        real*8 a,b,ar,a2,a4,a6,b2,b4,b6,fab,fba

        a=3.2d0*ua
        b=3.2d0*ub

        if (dabs(a-b) .lt. 1.0d-5) then
          ar   = (a+b)/2.0d0*r
          gams = dexp(-ar)*(48.0d0+33.0d0*ar+9.0d0*ar*ar+ar*ar*ar)
     &                                                /(r*48.0d0)
          return
        else
          a2   = a*a
          a4   = a2*a2
          a6   = a4*a2
          b2   = b*b
          b4   = b2*b2
          b6   = b4*b2
          fab  = a*b4/(2.0d0*(a2-b2)**2)-(b6-3.0d0*a2*b4)/((a2-b2)**3*r)
          fba  = b*a4/(2.0d0*(b2-a2)**2)-(a6-3.0d0*b2*a4)/((b2-a2)**3*r)
          gams = dexp(-a*r)*fab + dexp(-b*r)*fba
          return
        endif
      end


      function gamko(r,ua,ub)
        ! Klopman-Ohno gamma
        implicit none
        real*8 gamko,r,ua,ub
  
        gamko=(r*r+0.25d0*(1.0d0/ua+1.0d0/ub)**2)**(-0.5d0)

        return
      end
   
      function gamhx(r,ua,ub,zeta)
        ! gamma^h as defined in paper Gaus,Cui,Elstner,JCTC2011
        ! only for HX-interactions, for other interactions use gamnew
        implicit none
        real*8 gamhx,r,ua,ub,zeta,zero,a,b,gams
        parameter(zero=1.0d-4)
    
        a=3.2d0*ua
        b=3.2d0*ub 

        if(a+b .lt. zero) then
          gamhx = 0.0d0
          return
        elseif (r .lt. zero) then
          gamhx = 0.5d0*( a*b/(a+b)+a*a*b*b/((a+b)**3) )
          return
        else 
          gamhx = 1.0d0/r -gams(r,ua,ub)*dexp(-0.5d0*r*r*(ua+ub)**zeta)
          return
        endif
      end

    
      function dist1(x,a,c,k,i)
        implicit none
        include 'maxima.inc'
        real*8  dist1,x(3,nndim),rx,ry,rz,dist
        integer a,c,k,i 
        ! dist1 = dr_{ac}/dR_{ki}   ;   i is x,y,or z

        if ((a.eq.k).and.(c.eq.k)) then
          dist1 = 0.0d0
        elseif (a.eq.k) then
          dist1 = -(x(i,c)-x(i,a))/dist(x,a,c)
        elseif (c.eq.k) then
          dist1 =  (x(i,c)-x(i,a))/dist(x,a,c)
        else 
          dist1 = 0.0d0
        endif

        return
      end


      function gamz1(r,u)
        implicit none
        real*8 gamz1,r,u,zero,a
        parameter(zero=1.0d-4)

        if (r .lt. zero) then
          gamz1=0.0d0
          return
        else
          a=3.2d0*u
          gamz1 = -1.0d0/(r*r)+dexp(-a*r)*(1.0d0/(r*r)+a*a*0.5d0+a/r)
          return
        endif
      end


      function gamnew1(r,ua,ub)
        implicit none
        real*8 gamnew1,r,ua,ub,zero,a,b,gams1
        parameter(zero=1.0d-4)

        a=3.2d0*ua
        b=3.2d0*ub

        if((a+b.lt.zero).or.(r.lt.zero)) then
          gamnew1 = 0.0d0
          return
        else
          gamnew1 = -1.0d0/(r*r) - gams1(r,ua,ub)
          return
        endif
      end

 
      function gams1(r,ua,ub)
        implicit none
        real*8 gams1,r,ua,ub
        real*8 a,b,z,z2,zr,g,dgdr,a2,a4,b2,b4,fab,fba
        real*8 dfabdr,dfbadr

        a=3.2d0*ua
        b=3.2d0*ub

        if (dabs(a-b) .lt. 1.0d-5) then
          z     = (a+b)/2.0d0
          z2    = z*z
          zr    = z*r
          g     = (48.0d0+33.0d0*zr+9.0d0*zr*zr+zr*zr*zr)/(48.0d0*r)
          dgdr  = -1.0d0/(r*r)+3.0d0*z2/16.0d0+z2*zr/24.0d0
          gams1 = dexp(-zr)*(dgdr-z*g)
        else
          a2 = a*a
          a4 = a2*a2
          b2 = b*b
          b4 = b2*b2
          fab=a*b4/(2.0d0*(a2-b2)**2)-(b4*b2-3.0d0*a2*b4)/((a2-b2)**3*r)
          fba=b*a4/(2.0d0*(b2-a2)**2)-(a4*a2-3.0d0*b2*a4)/((b2-a2)**3*r)
          dfabdr = (b4*b2-3.0d0*a2*b4)/((a2-b2)**3*r*r)
          dfbadr = (a4*a2-3.0d0*b2*a4)/((b2-a2)**3*r*r)
          gams1  = dexp(-a*r)*(dfabdr-a*fab)+dexp(-b*r)*(dfbadr-b*fba)
        endif

        return
      end


      function gamko1(r,ua,ub)
        ! derivative of Klopman-Ohno gamma w.r.t. r
        implicit none
        real*8 gamko1,r,ua,ub

        gamko1=-r*(r*r+0.25d0*(1.0d0/ua+1.0d0/ub)**2)**(-1.5d0)

        return
      end

