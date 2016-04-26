c ======================================================================
c      Preditor-multicorretor para derivada de ordem 1.                |
c ====================================================================== 
      subroutine predict(nnode,ndf,id,u,v)
c **********************************************************************
c *                                                                    *
c *   PREDICT                                                          *
c *                                                                    *
c **********************************************************************
      implicit none
      include 'transiente.fi'
      integer nnode,ndf,id(ndf,*),i,j,k
      real*8 u(ndf,*),v(ndf,*),b
c ......................................................................      
      b = dt*(1.d0-alfa)
      do 110 i = 1, nnode
      do 100 j = 1, ndf
         k = id(j,i)
         if(k .gt. 0) then         
            u(j,i) = u(j,i) + v(j,i)*b
            v(j,i) = 0.d0
         endif
  100 continue
  110 continue
c ......................................................................    
      return
      end
c *********************************************************************
c
c *********************************************************************
      subroutine update(nnode,ndf,id,u,v,x)
c **********************************************************************
c *                                                                    *
c *   UPDATE                                                           *
c *                                                                    *
c **********************************************************************
      implicit none
      include 'transiente.fi'      
      integer nnode,ndf,id(ndf,*),i,j,k
      real*8 u(ndf,*),v(ndf,*),x(*),b
c ......................................................................            
      b = alfa*dt
      do 110 i = 1, nnode
         do 100 j = 1, ndf
            k = id(j,i)
            if(k .gt. 0) then
                u(j,i) =  u(j,i) + x(k)*b
                v(j,i) =  v(j,i) + x(k)
            endif
  100    continue
  110 continue
c ......................................................................  
      return
      end
c *********************************************************************
c
c *********************************************************************
c ======================================================================
c      Preditor-multicorretor para derivada de ordem 2.                |
c ====================================================================== 
      subroutine predict2(nnode,ndf,id,u,v,a,dt,alfa,beta)
c **********************************************************************
c *                                                                    *
c *   PREDICT2                                                         *
c *                                                                    *
c **********************************************************************
      integer nnode,ndf,k
      integer id(ndf,*)
      real*8  u(ndf,*),v(ndf,*),a(ndf,*)
      real*8  dt,alfa,beta,c,d
c ......................................................................
      c = dt*dt*.5d0*(1.d0-2.d0*beta)
      d = dt*(1.d0-alfa)
      do 110 i = 1, nnode
        do 100 j = 1, ndf
          k = id(j,i)
          if (k .gt. 0) then
            u(j,i) = u(j,i) + dt*v(j,i) + c*a(j,i)
            v(j,i) = v(j,i) + d*a(j,i)
            a(j,i) = 0.d0
          endif
  100   continue
  110 continue
      return
      end
c *********************************************************************
c
c *********************************************************************
      subroutine update2(nnode,ndf,id,u,v,a,x,dt,alfa,beta)
c **********************************************************************
c *                                                                    *
c *   UPDATE2                                                          *
c *                                                                    *
c **********************************************************************
      integer nnode,ndf,k
      integer id(ndf,*)
      real*8  u(ndf,*),v(ndf,*),a(ndf,*),x(*)
      real*8  dt,alfa,beta,c,d
c ......................................................................
      c = beta*dt*dt
      d = alfa*dt
      do 110 i = 1, nnode
        do 100 j = 1, ndf
          k = id(j,i)
          if (k .gt. 0) then
            a(j,i) = a(j,i) + x(k)
            u(j,i) = u(j,i) + c*x(k)
            v(j,i) = v(j,i) + d*x(k)
          endif
  100   continue
  110 continue
      return
      end
c *********************************************************************
c
c *********************************************************************
      subroutine predict_pm(nnode,nnodev,ndf,u,u0)
c **********************************************************************
c *                                                                    *
c *   PREDICT _PM                                                      *
c *                                                                    *
c **********************************************************************
      implicit none
      integer nnode,nnodev,ndf,i,j
      real*8 u(ndf,*),u0(ndf,*),b
c ......................................................................      
c
c ...
      do 110 i = 1, nnode
        do 100 j = 1, ndf - 1
          u(j,i) = u(j,i) - u0(j,i)
  100   continue
  110 continue
c .....................................................................
c
c ...
      do 120 i = 1, nnodev
        u(ndf,i) = u(ndf,i) - u0(ndf,i)  
  120 continue
c ......................................................................    
      return
      end
c *********************************************************************
c
c *********************************************************************
      subroutine update_pm(nnode,nnodev,ndf,id,u,x)
c **********************************************************************
c *                                                                    *
c *   UPDATEPM :  atualizacao do graus de liberdade do poro mecanico   *
c *                                                                    *
c **********************************************************************
      integer nnode,nnodev,ndf,k
      integer id(ndf,*)
      real*8  u(ndf,*),x(*)
c ......................................................................
c
c ... atualizacao dos incrementos dos deslocamentos
      do 110 i = 1, nnode
        do 100 j = 1, ndf-1
          k = id(j,i)
          if (k .gt. 0) then
            u(j,i) =  u(j,i) + x(k)
          endif
  100   continue
  110 continue
c ......................................................................
c
c ... atualizao da deltap
      do 120 i = 1, nnodev
        k = id(ndf,i)
        if (k .gt. 0) then
          u(ndf,i) = u(ndf,i) + x(k)
        endif
  120 continue
c ......................................................................
c
c ...
      return
      end
c *********************************************************************
c
c *********************************************************************
      subroutine update_mec(nnode,ndf,id,u,x)
c **********************************************************************
c *                                                                    *
c *   UPDATE_MEC : atualizacao dos graus de liberdade do mecanico      *
c *   estatico                                                         *
c *                                                                    *
c **********************************************************************
      integer nnode,nnodev,ndf,k
      integer id(ndf,*)
      real*8  u(ndf,*),x(*)
c ......................................................................
c
c ... atualizacao dos incrementos dos deslocamentos
      do 110 i = 1, nnode
        do 100 j = 1, ndf
          k = id(j,i)
          if (k .gt. 0) then
            u(j,i) = x(k)
          endif
  100   continue
  110 continue
c ......................................................................
c
c ...
      return
      end
c *********************************************************************
