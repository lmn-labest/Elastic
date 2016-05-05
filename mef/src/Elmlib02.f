c *********************************************************************
c * Biblioteca de elementos mecanicos                                 * 
c * ----------------------------------------------------------------- *
c * Elastico-estatico:                                                *
c * ----------------------------------------------------------------- *
c * ------------------ Elementos lineares --------------------------- *
c *                                                                   *
c * ELMT02_MEC- triangulo de 3 nos para o problema mecanico elastico  *
c * estatico (Estado plano de deformacao)                             *
c *                                                                   *
c * ELMT03_MEC- triangulo de 3 nos para o problema mecanico elastico  *
c * estatico (Estado plano de deformacao)                             *
c *                                                                   *
c * ELMT04_MEC- quadrilateros de 4 nos para o problema mecanico       * 
c * elastico estatico (Estado plano de deformacao)                    *
c *                                                                   *
c * ELMT05_MEC- quadrilateros de 4 nos para o problema mecanico       *
c * elastico estatico (Estado plano de tensao)                        *
c *                                                                   *
c * ELMT06_MEC- tetraedros de 4 nos para o problema mecanico elastico *        
c * estatico                                                          *
c *                                                                   *
c * ELMT07_MEC- hexaedros de  8 nos para o problema mecanico  elastico* 
c * estatico                                                          *
c * ------------------ Elementos quadraticos ------------------------ *
c * ELMT12_MEC- tetraedros de 10 nos para o problema mecanico elastico*
c * estatico                                                          *
c *                                                                   *
c * ELMT13_MEC- hexaedros de 20 nos para o problema mecanico  elastico* 
c * estatico                                                          *
c *                                                                   *
c * ----------------------------------------------------------------- *
c * Plastico-estatico:                                                *
c * ----------------------------------------------------------------- *
c *********************************************************************
      subroutine elmt02_mec(e,iq,x,u,p,s,txn,ndm,nst,nel,isw)
c **********************************************************************
c * Data de criacao    : 25/04/2016                                    *
c * Data de modificaco : 00/00/0000                                    * 
c * ------------------------------------------------------------------ *
c * ELMT02_MEC: Elemento triangulares de 3 nos para problemas          *  
c * mecanico estaico-elasticos (Estado plano de deformacao)            *
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ * 
c * e(10) - constantes fisicas                                         *
c *           e(1) = modulo de Young                                   *
c *           e(2) = coeficiente de Poisson                            *
c *           e(3) =                                                   *
c *           e(4) =                                                   *
c *           e(5) =                                                   *
c *           e(6) =                                                   *
c *           e(7) =                                                   *
c * x(ndm,nem) - coordenadas nodais locais                             *
c * u(nst)     - graus de liberade por elemento (u)                    *
c * p(nst)     - nao definido                                          *
c * s(nst,nst) - nao definido                                          *
c * txn(6,nen) - tensoes nodais                                        *      
c * ndm - dimensao                                                     *
c * nst - numero de graus de liberdade por elemento (2*3)              *
c * nel - numero do elemento                                           *
c * isw - codigo de instrucao                                          *
c *  1 =                                                               *
c *  2 = matriz K e forcas internas Ku                                 *
c *  3 = tensoes                                                       *
c *  4 = forcas de volume e superficies                                *
c *  5 =                                                               *
c *  6 =                                                               *
c *  7 =                                                               *
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ * 
c * e - constantes fisicas                                             *
c * s - matriz de elemento                                             *
c * p - isw = 2  residuo                                               *
c *     isw = 3  tensao                                                *
c *     isw = 4  cargas de superfice e volume                          *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ * 
c **********************************************************************
      implicit none
      include 'load.fi'
      integer ndm,nst,nel,isw
      integer i,j,l1,l2,k1,k2,k,tp
      integer nen
      integer iq(*)
c ...
      real*8 face_f(2),ddum
c ...
      real*8 u(*)
      real*8 p(*),s(nst,*)
c ... funcoes de interpolacao
      real*8 h(3),hx(3),hy(3)
c ...
      real*8 xj(2,2),xji(2,2)
      real*8 w,wt1,wt2,wt3,det
      real*8 epsi(4),txn(4,*)
c ... matriz jacobiana
      real*8 xj11,xj12,xj21,xj22
      real*8 xji11,xji12,xji21,xji22
c ...
      real*8 e(*),x(ndm,*),xf(2,2)
c ...
      real*8 a,b,c,ym,ps
      real*8 a1,a2,a3,d11,d12,d22,d33
c ... 
      integer tria_side_node3(2,3),no
c ...
      data tria_side_node3 / 1, 2
     .                     , 2, 3
     .                     , 3, 1/
c ... 
      data nen/3/
c ......................................................................
      goto (100,200,200,400,500) isw
c ======================================================================
c
c ...                                          
  100 continue
c ...
c .....................................................................
      return
c ======================================================================
c
c ... Matriz de rigidez:
c
c ......................................................................
  200 continue
c
c     Matriz Jacobiana:
c
      xj11 = x(1,1)-x(1,3)
      xj12 = x(2,1)-x(2,3)
      xj21 = x(1,2)-x(1,3)
      xj22 = x(2,2)-x(2,3)
      det  = xj11*xj22-xj12*xj21
      if (det .le. 0.d0) go to 1100
c
c	Inversa da matriz Jacobiana:
c
      xji11 =  xj22/det
      xji12 = -xj12/det
      xji21 = -xj21/det
      xji22 =  xj11/det
c
c     Derivadas das funcoes de interpolacao:
c
      hx(1) =  xji11
      hx(2) =  xji12
      hx(3) = -xji11-xji12
      hy(1) =  xji21
      hy(2) =  xji22
      hy(3) = -xji21-xji22
c
c     Matriz constitutiva:
c
      ym       = e(1)
      ps       = e(2)
c ...
      a1       = 1.d0 - ps
      a2       = 1.d0 - 2.d0*ps
      a3       = 1.d0 + ps
c ...
      a        = (ym*a1)/(a3*a2)
      b        = ps/a1
      c        = 0.5d0*(a2/a1) 
c ....
      d11      = a
      d12      = a*b
      d22      = a
      d33      = a*c  
c .....................................................................
c
c ...
      if(isw .eq. 3) go to 300
c .....................................................................
c
c ... Matriz de rigidez:
      do i = 1, nst
        do j = 1, nst
          s(j,i) = 0.d0
        enddo
      enddo
c .....................................................................
c
c ...
      wt1 = 0.5d0*det
      do 220 j = 1, 3
c        k1 = (j-1)*2+1
         k1 = 2*j-1
         k2 = k1 + 1    
         do 210 i = 1, 3
c          l1 = (i-1)*2+1
           l1 = 2*i-1
           l2 = l1 + 1    
c                                                                    
           s(l1,k1) = ( hx(i)*d11*hx(j) + hy(i)*d33*hy(j) ) * wt1
c                                                                     
           s(l1,k2) = ( hx(i)*d12*hy(j) + hy(i)*d33*hx(j) ) * wt1
c                                                                      
           s(l2,k1) = ( hy(i)*d12*hx(j) + hx(i)*d33*hy(j) ) * wt1
c                                                                       
           s(l2,k2) = ( hy(i)*d22*hy(j) + hx(i)*d33*hx(j) ) * wt1
c                                                                      
  210    continue
  220 continue
c .....................................................................
c
c ... Forcas internas:
      call lku_m(s,u,p,nst)
c .....................................................................
c
c .....................................................................      
      return  
c ======================================================================
c
c ... Tensoes nodais :
c
c ......................................................................
  300 continue
c ... tensao nodal | sxx syy szz sxy|
      call deform2d(hx,hy,u,epsi,3)
      call stress2d_m(d11,d12,d22,d33,ps,.true.,epsi,p)
c ...
      p(5)  = p(1)
      p(6)  = p(2)
      p(7)  = p(3)
      p(8)  = p(4)
c ...
      p(9)  = p(1)
      p(10) = p(2)
      p(11) = p(3)
      p(12) = p(4)
c .....................................................................
  310 continue
c .....................................................................
c
c ......................................................................
      return
c ======================================================================
c
c ... Cargas distribuidas no volume e no contorno:
c
c ......................................................................
 400  continue
c ... forca e fluxo distribuida no contorno
c     iq(1) = 1 | no 1 2 |
c             2 | no 2 3 |
c             3 | no 3 4 |     
c             4 | no 4 1 |
c
c ... verifica se ha alguma lado com carga
      tp = 0
      do 435 i = 1, 3
        tp    = tp + iq(i)
  435 continue
c .....................................................................
c
c ...
      if( tp .gt. 0 ) then
        do 440 j = 1, 3
c ... lado 
          if(iq(j) .gt. 0 ) then
c ...
            do 445 i = 1, 2 
              no      = tria_side_node3(i,j)
              xf(1,i) = x(1,no)
              xf(2,i) = x(2,no)
  445       continue
c .....................................................................
c
c ... jacobiano
            det = (xf(1,2)-xf(1,1))*(xf(1,2)-xf(1,1))
     .          + (xf(2,2)-xf(2,1))*(xf(2,2)-xf(2,1))
            if (det .le. 0.d0) goto 1000
            det = dsqrt(det)
c ......................................................................
c
c ... forcas
            call tload(iq(j),0.d0,0.d0,ddum,face_f)
c.......................................................................
            wt1 = 0.5d0*det
            do 455 i = 1, 2
              no    = tria_side_node3(i,j)
c              l1  = (no-1)*3+1
c              wt1   = h(i)*w
              l1    = 2*no-1
              l2    = l1 + 1
              p(l1) = p(l1) - face_f(1)*wt1
              p(l2) = p(l2) - face_f(2)*wt1
  455       continue
c .....................................................................
          endif
c .....................................................................
  440   continue
c .....................................................................
      endif
c .....................................................................
      return
c ======================================================================
c
c ... Tensoes iniciais:                  
c
c ......................................................................
  500 continue
      return
c ======================================================================
c
c ... Matriz geometrica:
c
c ......................................................................
  600 continue
c ......................................................................
 1000 continue
      print*, '*** Subrotina ELMT02_MEC: determinante nulo ou negativo 
     . do elemento ',nel,' na carga da aresta',j
      stop
 1100 continue
      print*, '*** Subrotina ELMT02_mec: determinante nulo ou negativo 
     .do elemento ',nel
      stop
      end
c **********************************************************************
      subroutine elmt03_mec(e,iq,x,u,p,s,txn,ndm,nst,nel,isw)
c **********************************************************************
c * Data de criacao    : 27/04/2016                                    *
c * Data de modificaco : 00/00/0000                                    * 
c * ------------------------------------------------------------------ *
c * ELMT03_MEC: Elemento triangulares de 3 nos para problemas          *  
c * mecanico estaico-elasticos (Estado plano de tensao)                *
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ * 
c * e(10) - constantes fisicas                                         *
c *           e(1) = modulo de Young                                   *
c *           e(2) = coeficiente de Poisson                            *
c *           e(3) = espessura                                         *
c *           e(4) =                                                   *
c *           e(5) =                                                   *
c *           e(6) =                                                   *
c *           e(7) =                                                   *
c * x(ndm,nem) - coordenadas nodais locais                             *
c * u(nst)     - graus de liberade por elemento (u)                    *
c * p(nst)     - nao definido                                          *
c * s(nst,nst) - nao definido                                          *
c * txn(4,nen) - tensoes nodais                                        *      
c * ndm - dimensao                                                     *
c * nst - numero de graus de liberdade por elemento (2*3)              *
c * nel - numero do elemento                                           *
c * isw - codigo de instrucao                                          *
c *  1 =                                                               *
c *  2 = matriz K e forcas internas Ku                                 *
c *  3 = tensoes                                                       *
c *  4 = forcas de volume e superficies                                *
c *  5 =                                                               *
c *  6 =                                                               *
c *  7 =                                                               *
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ * 
c * e - constantes fisicas                                             *
c * s - matriz de elemento                                             *
c * p - isw = 2  residuo                                               *
c *     isw = 3  tensao                                                *
c *     isw = 4  cargas de superfice e volume                          *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ * 
c **********************************************************************
      implicit none
      include 'load.fi'
      integer ndm,nst,nel,isw
      integer i,j,l1,l2,k1,k2,k,tp,tp1
      integer nen
      integer iq(*)
c ...
      real*8 face_f(2),ddum
c ...
      real*8 u(*)
      real*8 p(*),s(nst,*)
c ... funcoes de interpolacao
      real*8 h(3),hx(3),hy(3)
c ...
      real*8 xj(2,2),xji(2,2)
      real*8 wt1,det
      real*8 epsi(4),txn(4,*)
c ... matriz jacobiana
      real*8 xj11,xj12,xj21,xj22
      real*8 xji11,xji12,xji21,xji22
c ...
      real*8 e(*),x(ndm,*),xf(2,2)
c ...
      real*8 a,b,c,ym,ps,thic
      real*8 a1,a2,d11,d12,d22,d33
c ... 
      integer tria_side_node3(2,3),no
c ...
      data tria_side_node3 / 1, 2
     .                     , 2, 3
     .                     , 3, 1/
c ... 
      data nen/3/
c ......................................................................
      goto (100,200,200,400,500) isw
c ======================================================================
c
c ...                                          
  100 continue
c ...
c .....................................................................
      return
c ======================================================================
c
c ... Matriz de rigidez:
c
c ......................................................................
  200 continue
c
c     Matriz Jacobiana:
c
      xj11 = x(1,1)-x(1,3)
      xj12 = x(2,1)-x(2,3)
      xj21 = x(1,2)-x(1,3)
      xj22 = x(2,2)-x(2,3)
      det  = xj11*xj22-xj12*xj21
      if (det .le. 0.d0) go to 1100
c
c	Inversa da matriz Jacobiana:
c
      xji11 =  xj22/det
      xji12 = -xj12/det
      xji21 = -xj21/det
      xji22 =  xj11/det
c
c     Derivadas das funcoes de interpolacao:
c
      hx(1) =  xji11
      hx(2) =  xji12
      hx(3) = -xji11-xji12
      hy(1) =  xji21
      hy(2) =  xji22
      hy(3) = -xji21-xji22
c
c     Matriz constitutiva:
c
      ym       = e(1)
      ps       = e(2)
      thic     = e(3)
      if( thic .eq. 0 ) thic = 1.d0
c ...
      a1       = 1.d0 - ps
      a2       = 1.d0 - ps*ps
c ...
      a        = ym/a2
      b        = ps
      c        = 0.5d0*a1 
c ....
      d11      = a
      d12      = a*b
      d22      = a
      d33      = a*c  
c .....................................................................
c
c ...
      if(isw .eq. 3) go to 300
c .....................................................................
c
c ... Matriz de rigidez:
      do i = 1, nst
        do j = 1, nst
          s(j,i) = 0.d0
        enddo
      enddo
c .....................................................................
c
c ...
      wt1 = 0.5d0*det
      do 220 j = 1, 3
c        k1 = (j-1)*2+1
         k1 = 2*j-1
         k2 = k1 + 1    
         do 210 i = 1, 3
c          l1 = (i-1)*2+1
           l1 = 2*i-1
           l2 = l1 + 1    
c                                                                    
           s(l1,k1) = ( hx(i)*d11*hx(j) + hy(i)*d33*hy(j) ) * wt1
c                                                                     
           s(l1,k2) = ( hx(i)*d12*hy(j) + hy(i)*d33*hx(j) ) * wt1
c                                                                      
           s(l2,k1) = ( hy(i)*d12*hx(j) + hx(i)*d33*hy(j) ) * wt1
c                                                                       
           s(l2,k2) = ( hy(i)*d22*hy(j) + hx(i)*d33*hx(j) ) * wt1
c                                                                      
  210    continue
  220 continue
c .....................................................................
c
c ... Forcas internas:
      call lku_m(s,u,p,nst)
c .....................................................................
c
c .....................................................................      
      return  
c ======================================================================
c
c ... Tensoes nodais :
c
c ......................................................................
  300 continue
c ... tensao nodal | sxx syy 0 sxy|
      call deform2d(hx,hy,u,epsi,3)
      call stress2d_m(d11,d12,d22,d33,ps,.false.,epsi,p)
c ...
      p(5)  = p(1)
      p(6)  = p(2)
      p(7)  = p(3)
      p(8)  = p(4)
c ...
      p(9)  = p(1)
      p(10) = p(2)
      p(11) = p(3)
      p(12) = p(4)
c .....................................................................
  310 continue
c .....................................................................
c
c ......................................................................
      return
c ======================================================================
c
c ... Cargas distribuidas no volume e no contorno:
c
c ......................................................................
 400  continue
c ... forca e fluxo distribuida no contorno
c     iq(1) = 1 | no 1 2 |
c             2 | no 2 3 |
c             3 | no 3 4 |     
c             4 | no 4 1 |
c
c ... verifica se ha alguma lado com carga
      tp = 0
      do 435 i = 1, 3
        tp    = tp + iq(i)
  435 continue
c .....................................................................
c
c ...
      if( tp .gt. 0 ) then
        do 440 j = 1, 3
c ... lado 
          if(iq(j) .gt. 0 ) then
c ...
            do 445 i = 1, 2 
              no      = tria_side_node3(i,j)
              xf(1,i) = x(1,no)
              xf(2,i) = x(2,no)
  445       continue
c .....................................................................
c
c ... jacobiano
            det = (xf(1,2)-xf(1,1))*(xf(1,2)-xf(1,1))
     .          + (xf(2,2)-xf(2,1))*(xf(2,2)-xf(2,1))
            if (det .le. 0.d0) goto 1000
            det = dsqrt(det)
c ......................................................................
c
c ... forcas
            call tload(iq(j),0.d0,0.d0,ddum,face_f)
c.......................................................................
            wt1 = 0.5d0*det
            do 455 i = 1, 2
              no    = tria_side_node3(i,j)
c              l1  = (no-1)*3+1
c              wt1   = h(i)*w
              l1    = 2*no-1
              l2    = l1 + 1
              p(l1) = p(l1) - face_f(1)*wt1
              p(l2) = p(l2) - face_f(2)*wt1
  455       continue
c .....................................................................
          endif
c .....................................................................
  440   continue
c .....................................................................
      endif
c .....................................................................
      return
c ======================================================================
c
c ... Tensoes iniciais:                  
c
c ......................................................................
  500 continue
      return
c ======================================================================
c
c ... Matriz geometrica:
c
c ......................................................................
  600 continue
c ......................................................................
 1000 continue
      print*, '*** Subrotina ELMT02_MEC: determinante nulo ou negativo 
     . do elemento ',nel,' na carga da aresta',j
      stop
 1100 continue
      print*, '*** Subrotina ELMT02_mec: determinante nulo ou negativo 
     .do elemento ',nel
      stop
      end
c **********************************************************************
      subroutine elmt04_mec(e,iq,x,u,p,s,txn,ndm,nst,nel,isw)
c **********************************************************************
c * Data de criacao    : 25/04/2016                                    *
c * Data de modificaco : 00/00/0000                                    * 
c * ------------------------------------------------------------------ *
c * ELMT04_MEC: Elemento quadrilateros 4 nos para problemas            *  
c * mecanico estaico-elasticos (Estado plano de deformacao)            *
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ * 
c * e(10) - constantes fisicas                                         *
c *           e(1) = modulo de Young                                   *
c *           e(2) = coeficiente de Poisson                            *
c *           e(3) =                                                   *
c *           e(4) =                                                   *
c *           e(5) =                                                   *
c *           e(6) =                                                   *
c *           e(7) =                                                   *
c * x(ndm,nem) - coordenadas nodais locais                             *
c * u(nst)     - graus de liberade por elemento (u)                    *
c * p(nst)     - nao definido                                          *
c * s(nst,nst) - nao definido                                          *
c * txn(6,nen) - tensoes nodais                                        *      
c * ndm - dimensao                                                     *
c * nst - numero de graus de liberdade por elemento (2*4)              *
c * nel - numero do elemento                                           *
c * isw - codigo de instrucao                                          *
c *  1 =                                                               *
c *  2 = matriz K e forcas internas Ku                                 *
c *  3 = tensoes                                                       *
c *  4 = forcas de volume e superficies                                *
c *  5 =                                                               *
c *  6 =                                                               *
c *  7 =                                                               *
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ * 
c * e - constantes fisicas                                             *
c * s - matriz de elemento                                             *
c * p - isw = 2  residuo                                               *
c *     isw = 3  tensao                                                *
c *     isw = 4  cargas de superfice e volume                          *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ * 
c **********************************************************************
      implicit none
      include 'load.fi'
      common /gauss/ pg, wg
      integer ndm,nst,nel,isw
      integer i,j,l1,l2,k1,k2,k,tp,tp1
      integer nen,nint,lx,ly,lz
      integer iq(*)
c ...
      real*8 face_f(2),ddum
c ...
      real*8 u(*)
      real*8 p(*),s(nst,*)
c ... funcoes de interpolacao
      real*8 h(4),hx(4),hy(4),hz(4)
c ...
      real*8 xj(2,2),xji(2,2),r(3,3)
      real*8 ri,si,ti,w,wt1,wt2,wt3,det
      real*8 rn(4),sn(4)
      real*8 pi,epsi(4),txn(4,*)
c ... integracao numerica de tetraedros          
      real*8 pg(10,10),wg(10,10)
c ...
      real*8 e(*),x(ndm,*),xf(2,2)
c ...
      real*8 a,b,c,ym,ps
      real*8 a1,a2,a3,d11,d12,d22,d33
c ... 
      integer quad_side_node4(2,4),no
c ...
      data quad_side_node4 / 1, 2
     .                     , 2, 3
     .                     , 3, 4
     .                     , 4, 1/
c ... 
c
      data rn / 1.d0,-1.d0,-1.d0, 1.d0/ ! r1, r2, r3, r4
c
      data sn / 1.d0, 1.d0,-1.d0,-1.d0/ ! s1, s2, s3, s4           
c
      data nen/4/
c ......................................................................
      goto (100,200,300,400,500) isw
c ======================================================================
c
c ...                                          
  100 continue
c ...
c .....................................................................
      return
c ======================================================================
c
c ... Matriz de rigidez:
c
c ......................................................................
  200 continue
c ... matriz constitutiva:
      ym       = e(1)
      ps       = e(2)
c ...
      a1       = 1.d0 - ps
      a2       = 1.d0 - 2.d0*ps
      a3       = 1.d0 + ps
c ...
      a        = (ym*a1)/(a3*a2)
      b        = ps/a1
      c        = 0.5d0*(a2/a1) 
c ....
      d11      = a
      d12      = a*b
      d22      = a
      d33      = a*c  
c .....................................................................
c
c ... Matriz de rigidez:
      do i = 1, nst
        do j = 1, nst
          s(j,i) = 0.d0
        enddo
      enddo
c .....................................................................
c
c ... 
      nint = 2 
      do 210 ly = 1, nint
        si = pg(ly,nint)
        do 215 lx = 1, nint
          ri = pg(lx,nint)
c ...                                               
          call sfquad4_m(h,hx,hy,ri,si,.false.,.true.)
          call jacob2d_m(hx,hy,xj,xji,x,det,4,ndm,nel,.true.)
c .....................................................................
c
c ...
          w   = wg(lx,nint)*wg(ly,nint)*det
c .....................................................................
c
c ... Kuu ( Int((Bt)*C*(B)*dV) )
          wt1 = w
          do 220 i = 1, 4
c             l1 = (i-1)*2+1
            l1 = 2*i-1
            l2 = l1 + 1
            do 225 j = 1, 4
c               k1 = (j-1)*2+1
              k1 = 2*j-1
              k2 = k1 + 1
c
              s(l1,k1) = s(l1,k1) 
     .                 + (hx(i)*d11*hx(j)+hy(i)*d33*hy(j))*wt1
c                                                                        
              s(l1,k2) = s(l1,k2) 
     .                 + (hx(i)*d12*hy(j)+hy(i)*d33*hx(j))*wt1
c                                                                      
              s(l2,k1) = s(l2,k1) 
     .                 + (hy(i)*d12*hx(j)+hx(i)*d33*hy(j))*wt1
c                                                                         
              s(l2,k2) = s(l2,k2)
     .                 + (hy(i)*d22*hy(j)+hx(i)*d33*hx(j))*wt1
c .....................................................................                   
  225       continue
c .....................................................................              
  220     continue 
c .....................................................................
  215   continue
c .....................................................................
  210 continue
c .....................................................................
c
c ... Forcas internas:
      call lku_m(s,u,p,nst)
c .....................................................................
c
c .....................................................................      
      return  
c ======================================================================
c
c ... Tensoes nodais :
c
c ......................................................................
  300 continue
c ... matriz constitutiva:
      ym        = e(1)
      ps        = e(2)
c ...
      a1        = 1.d0 - ps
      a2        = 1.d0 - 2.d0*ps
      a3        = 1.d0 + ps
c ...
      a         = (ym*a1)/(a3*a2)
      b         = ps/a1
      c         = 0.5d0*(a2/a1)
c .....................................................................
c
c ....
      d11      = a
      d12      = a*b
      d22      = a
      d33      = a*c  
c .....................................................................
c
c ... tensao nodal | sxx syy szz sxy |
      do 310 i = 1, 4
c       tp = (i-1)*4 + 1
        tp  = 4*i - 3
c ... calculo do terminante
        call sfquad4_m(h,hx,hy,rn(i),sn(i),.false.,.true.)
        call jacob2d_m(hx,hy,xj,xji,x,det,4,ndm,nel,.true.)
c .....................................................................
        call deform2d(hx,hy,u,epsi,4)
        call stress2d_m(d11,d12,d22,d33,ps,.true.,epsi,p(tp))
c .....................................................................
  310 continue
c .....................................................................
c
c ......................................................................
      return
c ======================================================================
c
c ... Cargas distribuidas no volume e no contorno:
c
c ......................................................................
 400  continue
c ... forca e fluxo distribuida no contorno
c     iq(1) = 1 | no 1 2 |
c             2 | no 2 3 |
c             3 | no 3 4 |     
c             4 | no 4 1 |
c
c ... verifica se ha alguma lado com carga
      tp = 0
      do 435 i = 1, 4
        tp    = tp + iq(i)
  435 continue
c .....................................................................
c
c ...
      if( tp .gt. 0 ) then
        do 440 j = 1, 4
c ... lado 
          if(iq(j) .gt. 0 ) then
c ...
            do 445 i = 1, 2 
              no      = quad_side_node4(i,j)
              xf(1,i) = x(1,no)
              xf(2,i) = x(2,no)
  445       continue
c .....................................................................
c
c ... jacobiano
            det = (xf(1,2)-xf(1,1))*(xf(1,2)-xf(1,1))
     .          + (xf(2,2)-xf(2,1))*(xf(2,2)-xf(2,1))
            if (det .le. 0.d0) goto 1000
            det = dsqrt(det)
c ......................................................................
c
c ... forcas
            call tload(iq(j),0.d0,0.d0,ddum,face_f)
c.......................................................................
            wt1 = 0.5d0*det
            do 455 i = 1, 2
              no      = quad_side_node4(i,j)
c              l1  = (no-1)*3+1
c              wt1   = h(i)*w
              l1    = 2*no-1
              l2    = l1 + 1
              p(l1) = p(l1) - face_f(1)*wt1
              p(l2) = p(l2) - face_f(2)*wt1
  455       continue
c .....................................................................
          endif
c .....................................................................
  440   continue
c .....................................................................
      endif
c .....................................................................
      return
c ======================================================================
c
c ... Tensoes iniciais:                  
c
c ......................................................................
  500 continue
      return
c ======================================================================
c
c ... Matriz geometrica:
c
c ......................................................................
  600 continue
c ......................................................................
 1000 continue
      print*, '*** Subrotina ELMT05_MEC: determinante nulo ou negativo 
     . do elemento ',nel,' na carga da aresta',j
      stop
      end
c **********************************************************************
      subroutine elmt05_mec(e,iq,x,u,p,s,txn,ndm,nst,nel,isw)
c **********************************************************************
c * Data de criacao    : 26/04/2016                                    *
c * Data de modificaco : 00/00/0000                                    * 
c * ------------------------------------------------------------------ *
c * ELMT04_MEC: Elemento quadrilateros 4 nos para problemas            *  
c * mecanico estaico-elasticos (Estado plano de tensao)                *
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ * 
c * e(10) - constantes fisicas                                         *
c *           e(1) = modulo de Young                                   *
c *           e(2) = coeficiente de Poisson                            *
c *           e(3) = espessura                                         *
c *           e(4) =                                                   *
c *           e(5) =                                                   *
c *           e(6) =                                                   *
c *           e(7) =                                                   *
c * x(ndm,nem) - coordenadas nodais locais                             *
c * u(nst)     - graus de liberade por elemento (u)                    *
c * p(nst)     - nao definido                                          *
c * s(nst,nst) - nao definido                                          *
c * txn(4,nen) - tensoes nodais                                        *      
c * ndm - dimensao                                                     *
c * nst - numero de graus de liberdade por elemento (2*4)              *
c * nel - numero do elemento                                           *
c * isw - codigo de instrucao                                          *
c *  1 =                                                               *
c *  2 = matriz K e forcas internas Ku                                 *
c *  3 = tensoes                                                       *
c *  4 = forcas de volume e superficies                                *
c *  5 =                                                               *
c *  6 =                                                               *
c *  7 =                                                               *
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ * 
c * e - constantes fisicas                                             *
c * s - matriz de elemento                                             *
c * p - isw = 2  residuo                                               *
c *     isw = 3  tensao                                                *
c *     isw = 4  cargas de superfice e volume                          *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ * 
c **********************************************************************
      implicit none
      include 'load.fi'
      common /gauss/ pg, wg
      integer ndm,nst,nel,isw
      integer i,j,l1,l2,k1,k2,k,tp,tp1
      integer nen,nint,lx,ly,lz
      integer iq(*)
c ...
      real*8 face_f(2),ddum
c ...
      real*8 u(*)
      real*8 p(*),s(nst,*)
c ... funcoes de interpolacao
      real*8 h(4),hx(4),hy(4),hz(4)
c ...
      real*8 xj(2,2),xji(2,2),r(3,3)
      real*8 ri,si,ti,w,wt1,wt2,wt3,det
      real*8 rn(4),sn(4)
      real*8 pi,epsi(4),txn(4,*)
c ... integracao numerica de tetraedros          
      real*8 pg(10,10),wg(10,10)
c ...
      real*8 e(*),x(ndm,*),xf(2,2)
c ...
      real*8 a,b,c,ym,ps,thic
      real*8 a1,a2,a3,d11,d12,d22,d33
c ... 
      integer quad_side_node4(2,4),no
c ...
      data quad_side_node4 / 1, 2
     .                     , 2, 3
     .                     , 3, 4
     .                     , 4, 1/
c ... 
c
      data rn / 1.d0,-1.d0,-1.d0, 1.d0/ ! r1, r2, r3, r4
c
      data sn / 1.d0, 1.d0,-1.d0,-1.d0/ ! s1, s2, s3, s4           
c
      data nen/4/
c ......................................................................
      goto (100,200,300,400,500) isw
c ======================================================================
c
c ...                                          
  100 continue
c ...
c .....................................................................
      return
c ======================================================================
c
c ... Matriz de rigidez:
c
c ......................................................................
  200 continue
c ... matriz constitutiva:
      ym       = e(1)
      ps       = e(2)
      thic     = e(3)
      if( thic .eq. 0 ) thic = 1.d0
c ...
      a1       = 1.d0 - ps
      a2       = 1.d0 - ps*ps
c ...
      a        = ym/a2
      b        = ps
      c        = 0.5d0*a1 
c ....
      d11      = a
      d12      = a*b
      d22      = a
      d33      = a*c  
c .....................................................................
c
c ... Matriz de rigidez:
      do i = 1, nst
        do j = 1, nst
          s(j,i) = 0.d0
        enddo
      enddo
c .....................................................................
c
c ... 
      nint = 2 
      do 210 ly = 1, nint
        si = pg(ly,nint)
        do 215 lx = 1, nint
          ri = pg(lx,nint)
c ...                                               
          call sfquad4_m(h,hx,hy,ri,si,.false.,.true.)
          call jacob2d_m(hx,hy,xj,xji,x,det,4,ndm,nel,.true.)
c .....................................................................
c
c ...
          w   = wg(lx,nint)*wg(ly,nint)*det*thic
c .....................................................................
c
c ... Kuu ( Int((Bt)*C*(B)*dV) )
          wt1 = w
          do 220 i = 1, 4
c             l1 = (i-1)*2+1
            l1 = 2*i-1
            l2 = l1 + 1
            do 225 j = 1, 4
c               k1 = (j-1)*2+1
              k1 = 2*j-1
              k2 = k1 + 1
c
              s(l1,k1) = s(l1,k1) 
     .                 + (hx(i)*d11*hx(j)+hy(i)*d33*hy(j))*wt1
c                                                                         
              s(l1,k2) = s(l1,k2) 
     .                 + (hx(i)*d12*hy(j)+hy(i)*d33*hx(j))*wt1
c                                                                      
              s(l2,k1) = s(l2,k1) 
     .                 + (hy(i)*d12*hx(j)+hx(i)*d33*hy(j))*wt1
c                                                                         
              s(l2,k2) = s(l2,k2)
     .                 + (hy(i)*d22*hy(j)+hx(i)*d33*hx(j))*wt1
c .....................................................................                   
  225       continue
c .....................................................................              
  220     continue 
c .....................................................................
  215   continue
c .....................................................................
  210 continue
c .....................................................................
c
c ... Forcas internas:
      call lku_m(s,u,p,nst)
c .....................................................................
c
c .....................................................................      
      return  
c ======================================================================
c
c ... Tensoes nodais :
c
c ......................................................................
  300 continue
c ... matriz constitutiva:
      ym       = e(1)
      ps       = e(2)
      thic     = e(3)
      if( thic .eq. 0 ) thic = 1.d0
c ...
      a1       = 1.d0 - ps
      a2       = 1.d0 - ps*ps
c ...
      a        = ym/a2
      b        = ps
      c        = 0.5d0*a1 
c ....
      d11      = a
      d12      = a*b
      d22      = a
      d33      = a*c  
c .....................................................................
c
c ... tensao nodal | sxx syy szz sxy |
      do 310 i = 1, 4
c       tp = (i-1)*4 + 1
        tp  = 4*i - 3
c ... calculo do terminante
        call sfquad4_m(h,hx,hy,rn(i),sn(i),.false.,.true.)
        call jacob2d_m(hx,hy,xj,xji,x,det,4,ndm,nel,.true.)
c .....................................................................
        call deform2d(hx,hy,u,epsi,4)
        call stress2d_m(d11,d12,d22,d33,ps,.false.,epsi,p(tp))
c .....................................................................
  310 continue
c .....................................................................
c
c ......................................................................
      return
c ======================================================================
c
c ... Cargas distribuidas no volume e no contorno:
c
c ......................................................................
 400  continue
c ... forca e fluxo distribuida no contorno
c     iq(1) = 1 | no 1 2 |
c             2 | no 2 3 |
c             3 | no 3 4 |     
c             4 | no 4 1 |
c
c ... verifica se ha alguma lado com carga
      tp = 0
      do 435 i = 1, 4
        tp    = tp + iq(i)
  435 continue
c .....................................................................
c
c ...
      if( tp .gt. 0 ) then
        do 440 j = 1, 4
c ... lado 
          if(iq(j) .gt. 0 ) then
c ...
            do 445 i = 1, 2 
              no      = quad_side_node4(i,j)
              xf(1,i) = x(1,no)
              xf(2,i) = x(2,no)
  445       continue
c .....................................................................
c
c ... jacobiano
            det = (xf(1,2)-xf(1,1))*(xf(1,2)-xf(1,1))
     .          + (xf(2,2)-xf(2,1))*(xf(2,2)-xf(2,1))
            if (det .le. 0.d0) goto 1000
            det = dsqrt(det)
c ......................................................................
c
c ... forcas
            call tload(iq(j),0.d0,0.d0,ddum,face_f)
c.......................................................................
            wt1 = 0.5d0*det
            do 455 i = 1, 2
              no      = quad_side_node4(i,j)
c             l1  = (no-1)*2+1
              l1    = 2*no-1
              l2    = l1 + 1
              p(l1) = p(l1) - face_f(1)*wt1
              p(l2) = p(l2) - face_f(2)*wt1
  455       continue
c .....................................................................
          endif
c .....................................................................
  440   continue
c .....................................................................
      endif
c .....................................................................
      return
c ======================================================================
c
c ... Tensoes iniciais:                  
c
c ......................................................................
  500 continue
      return
c ======================================================================
c
c ... Matriz geometrica:
c
c ......................................................................
  600 continue
c ......................................................................
 1000 continue
      print*, '*** Subrotina ELMT05_MEC: determinante nulo ou negativo 
     . do elemento ',nel,' na carga da aresta',j
      stop
      end
c **********************************************************************
      subroutine elmt06_mec(e,iq,x,u,p,s,txn,ndm,nst,nel,isw)
c **********************************************************************
c * Data de criacao    : 24/04/2016                                    *
c * Data de modificaco : 00/00/0000                                    * 
c * ------------------------------------------------------------------ *      
c * ELMT06_PM: Elemento tetraedrico de  4 nos para problemas           *  
c * mecanico elasticos                                                 *
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ * 
c * e(10) - constantes fisicas                                         *
c *           e(1) = modulo de Young                                   *
c *           e(2) = coeficiente de Poisson                            *
c *           e(3) = massa especifica do meio poroso                   *
c * x(ndm,nem) - coordenadas nodais locais                             *
c * u(nst)     - graus de liberade por elemento (u)                    *
c * p(nst)     - nao definido                                          *
c * s(nst,nst) - nao definido                                          *
c * txn(6,nen) - tensoes nodais                                        *
c * ndm - dimensao                                                     *
c * nst - numero de graus de liberdade por elemento (3*4)              *
c * nel - numero do elemento                                           *
c * isw - codigo de instrucao                                          *
c *  1 =                                                               *
c *  2 = matriz K e forcas internas Ku                                 *
c *  3 = tensoes                                                       *
c *  4 = forcas de volume e superficies                                *
c *  5 = Tensoes iniciais                                              *
c *  6 =                                                               *
c *  7 =                                                               *
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ * 
c * e - constantes fisicas                                             *
c * s - matriz de elemento                                             *
c * p - isw = 2  residuo                                               *
c *     isw = 3  tensao                                                *
c *     isw = 4  cargas de superfice e volume                          *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ * 
c **********************************************************************
      implicit none
      include 'gravity.fi'
      include 'load.fi'
      real*8  div6,div144       
      parameter ( div6   = 0.166666666666667d+00)
      parameter ( div144 = 0.364444444444444d-02)
      integer ndm,nst,nel,isw
      integer i,j,l1,l2,l3,l,k1,k2,k3,k,tp,tp1
      integer nen,nint,lx,ly,lz
      integer iq(*)
c ...
      real*8 face_f(3),ddum
c ...
      real*8 u(*)
      real*8 p(*),s(nst,*)
c ... funcoes de interpolacao
      real*8 h(4),hx(4),hy(4),hz(4)
c ...
      real*8 xj(3,3),xji(3,3),r(3,3)
      real*8 ri,si,ti,w,wt1,det
      real*8 rn(10),sn(10),tn(10)
      real*8 pi,epsi(6),txi(6),txn(6,*),pm
c ...
      real*8 e(*),x(ndm,*),xf(3,3)
c ...
      real*8 a,b,c,ym,ps
      real*8 density,gl(3)
      real*8 a1,a2,a3,tmp
c ... 
      integer tetra_face_node4(3,4),no
c ...
      data tetra_face_node4 / 2, 3, 4
     .                      , 1, 4, 3
     .                      , 1, 2, 4
     .                      , 1, 3, 2/
c
      data nen/4/    
c ......................................................................
      goto (100,200,300,400,500) isw
c ======================================================================
c
c.... calculo do delta t critico               
c
c ......................................................................
  100 continue
c .....................................................................
      return
c ======================================================================
c
c ... Matriz de rigidez:
c
c ......................................................................
  200 continue
c ... matriz constitutiva:
      ym       = e(1)
      ps       = e(2)
c ...
      a1       = 1.d0 - ps
      a2       = 1.d0 - 2.d0*ps
      a3       = 1.d0 + ps
c ...
      a        = (ym*a1)/(a3*a2)
      b        = ps/a1
      c        = 0.5d0*(a2/a1) 
c .....................................................................
c
c ... Matriz de rigidez:
      do i = 1, nst
        do j = 1, nst
          s(j,i) = 0.d0
        enddo
      enddo
c .....................................................................
c
c ... Matriz Jacobiana:
      call jtetra4(x,hx,hy,hz,det,.true.,nel)
c .....................................................................
c
c ... Kuu ( Int((Bt)*C*(B)*dV) )
      wt1 = a*det*div6
      do 215 i = 1, 4
c       l = (i-1)*3+1
        l1 = 3*i - 2
        l2 = l1 + 1
        l3 = l2 + 1
        do 205 j = 1, 4
c         k1 = (j-1)*3+1
          k1 = 3*j-2
          k2 = k1 + 1
          k3 = k2 + 1
c
          s(l1,k1) = (hx(i)*hx(j) + c*hy(i)*hy(j) + c*hz(i)*hz(j))*wt1
c
          s(l1,k2) = (b*hx(i)*hy(j) + c*hy(i)*hx(j))*wt1
c
          s(l1,k3) = (b*hx(i)*hz(j) + c*hz(i)*hx(j))*wt1
c
          s(l2,k1) = (b*hy(i)*hx(j) + c*hx(i)*hy(j))*wt1
c
          s(l2,k2) = (hy(i)*hy(j) + c*hx(i)*hx(j) + c*hz(i)*hz(j))*wt1
c
          s(l2,k3) = (b*hy(i)*hz(j) + c*hz(i)*hy(j))*wt1
c
          s(l3,k1) = (b*hz(i)*hx(j) + c*hx(i)*hz(j))*wt1
c
          s(l3,k2) = (b*hz(i)*hy(j) + c*hy(i)*hz(j))*wt1
c
          s(l3,k3) = (hz(i)*hz(j) + c*hy(i)*hy(j) + c*hx(i)*hx(j))*wt1
 205    continue
 215  continue
c
c ... Forcas internas:
      call lku_m(s,u,p,nst)
c     call lku_sym(s,u,p,nst)    
c ......................................................................
c
c ......................................................................   
      return  
c ======================================================================
c
c ... Tensoes nodais:
c
c ......................................................................
  300 continue
c ... matriz constitutiva:
      ym        = e(1)
      ps        = e(2)
c ...
      a1        = 1.d0 - ps
      a2        = 1.d0 - 2.d0*ps
      a3        = 1.d0 + ps
c ...
      a         = (ym*a1)/(a3*a2)
      b         = ps/a1
      c         = 0.5d0*(a2/a1)
c .....................................................................
c
c ... Matriz Jacobiana:
      call jtetra4(x,hx,hy,hz,det,.true.,nel)
      call deform3d(hx,hy,hz,u,epsi,4)
      call stress3d(a,b,c,epsi,p)
c .....................................................................
c
c ... tensao nodal total
      p(7)  = p(1)
      p(8)  = p(2)
      p(9)  = p(3)
      p(10) = p(4)
      p(11) = p(5)
      p(12) = p(6)
c
      p(13) = p(1)
      p(14) = p(2)
      p(15) = p(3)
      p(16) = p(4)
      p(17) = p(5)
      p(18) = p(6)
c
      p(19) = p(1)
      p(20) = p(2)
      p(21) = p(3)
      p(22) = p(4)
      p(23) = p(5)
      p(24) = p(6)
c .....................................................................
  310 continue
c .....................................................................
c
c ......................................................................
      return
c ======================================================================
c
c ... Cargas distribuidas no volume e no contorno:
c
c ......................................................................
 400  continue
c ... peso proprio
      if( density .eq. 0.0d0) goto 430
      density   = e(3)*1.0d-06 
c .....................................................................
c
c ... 
      gl(1)     =  gravity(1)*density
      gl(2)     =  gravity(2)*density
      gl(3)     =  gravity(3)*density
c .....................................................................
c
c ... Matriz Jacobiana:
      call jtetra4(x,hx,hy,hz,det,.false.,nel)
c .....................................................................          
c
c ...
c     wt1 = det*(1.0d0/6.0d0)*(1.0d0/24.0d0) 
      wt1 = det*div144 
      do 420 i = 1, 4
c         l1  = (i-1)*3+1
          l1  = 3*i-2
          l2  = l1 + 1
          l3  = l2 + 1
c ... Fu = int(pm_d*(N~T)*ge*dv)  
          p(l1) = - wt1*gl(1)
          p(l2) = - wt1*gl(2)
          p(l3) = - wt1*gl(3)
  420   continue    
c .....................................................................
c
 410  continue   
c .....................................................................
c
 430  continue
c ... forca e fluxo distribuida no contorno
c     iq(1) = 1 | no 2 3 4  |
c             2 | no 1 4 3  |
c             3 | no 1 2 4  |     
c             4 | no 1 3 2  |
c ... verifica se ha alguma face com carga      
      tp = 0
      do 435 i = 1, 4
        tp  = tp + iq(i)
  435 continue
c .....................................................................
c
c ...      
      if( tp .gt. 0 ) then
        do 440 j = 1, 4
c ... face  
          if(iq(j) .gt. 0 ) then
c ...
            do 445 i = 1, 3
              no       = tetra_face_node4(i,j)
              xf(1,i) = x(1,no) 
              xf(2,i) = x(2,no) 
              xf(3,i) = x(3,no)
  445       continue
c ...    
            call rotmatrix(xf,r)
            call rotate(xf(1,1),r,xf(1,1),.false.)
            call rotate(xf(1,2),r,xf(1,2),.false.)
            call rotate(xf(1,3),r,xf(1,3),.false.)
c ... forcas
            call tload(iq(j),0.d0,0.0d0,ddum,face_f)
c ...
            call jtria3(xf,hx,hy,hz,det,.false.,nel,3)
c ...
            wt1 = det*div6
            do 450 i = 1, 3 
                no  = tetra_face_node4(i,j)
c               l1  = (no-1)*3+1
                l1    = 3*no-2
                l2    = l1 + 1
                l3    = l2 + 1
                p(l1) = p(l1) - face_f(1)*wt1
                p(l2) = p(l2) - face_f(2)*wt1
                p(l3) = p(l3) - face_f(3)*wt1
  450       continue
c .....................................................................
          endif
c .....................................................................
  440   continue
c .....................................................................
      endif
c .....................................................................
      return
c ======================================================================
c
c ... Tensoes iniciais:                  
c
c ......................................................................
  500 continue
c ... Matriz Jacobiana:
      call jtetra4(x,hx,hy,hz,det,.true.,nel)
c .....................................................................
c
c ...
      wt1 = det*div6
      do 520 i = 1, 4
c       l1  = (i-1)*3+1
        l1  = 3*i-2
        l2  = l1 + 1
        l3  = l2 + 1
c ... Fu = int(BeT*sigma*dv)
        p(l1) = wt1*(hx(i)*txn(1,i) + hy(i)*txn(4,i) + hz(i)*txn(6,i))
c         
        p(l2) = wt1*(hx(i)*txn(4,i) + hy(i)*txn(2,i) + hz(i)*txn(5,i))
c          
        p(l3) = wt1*(hx(i)*txn(6,i) + hy(i)*txn(5,i) + hz(i)*txn(3,i))
c ....................................................................
  520 continue  
c .....................................................................
      return
c ======================================================================
c
c ... Matriz geometrica:
c
c ......................................................................
  600 continue
      stop
c ....................................................................
      end
c *********************************************************************
c
c *********************************************************************
      subroutine elmt07_mec(e,iq,x,u,p,s,txn,ndm,nst,nel,isw)
c *********************************************************************
c * Data de criacao    : 24/04/2016                                   *
c * Data de modificaco : 00/00/0000                                    * 
c * ------------------------------------------------------------------ *       
c * ELMT07_MEC: Elemento hexaedricos de 8 nos para problemas           *  
c * mecanico estaico-elasticos                                         *
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ * 
c * e(10) - constantes fisicas                                         *
c *           e(1) = modulo de Young                                   *
c *           e(2) = coeficiente de Poisson                            *
c *           e(3) = massa especifica do meio                          *
c *           e(4) =                                                   *
c *           e(5) =                                                   *
c *           e(6) =                                                   *
c *           e(7) =                                                   *
c * x(ndm,nem) - coordenadas nodais locais                             *
c * u(nst)     - graus de liberade por elemento (u)                    *
c * p(nst)     - nao definido                                          *
c * s(nst,nst) - nao definido                                          *
c * txn(6,nen) - tensoes nodais                                        *      
c * ndm - dimensao                                                     *
c * nst - numero de graus de liberdade por elemento (3*8)              *
c * nel - numero do elemento                                           *
c * isw - codigo de instrucao                                          *
c *  1 =                                                               *
c *  2 = matriz K e forcas internas Ku                                 *
c *  3 = tensoes                                                       *
c *  4 = forcas de volume e superficies                                *
c *  5 = Tensoes iniciais                                              *
c *  6 =                                                               *
c *  7 =                                                               *
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ * 
c * e - constantes fisicas                                             *
c * s - matriz de elemento                                             *
c * p - isw = 2  residuo                                               *
c *     isw = 3  tensao       o                                        *
c *     isw = 4  cargas de superfice e volume                          *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ * 
c **********************************************************************
      implicit none
      include 'gravity.fi'
      include 'load.fi'
      common /gauss/ pg, wg
      integer ndm,nst,nel,isw
      integer i,j,l1,l2,l3,l,k1,k2,k3,k,tp,tp1
      integer nen,nint,lx,ly,lz
      integer iq(*)
c ...
      real*8 face_f(3),ddum
c ...
      real*8 u(*)
      real*8 p(*),s(nst,*)
c ... funcoes de interpolacao
      real*8 h(8),hx(8),hy(8),hz(8)
c ...
      real*8 xj(3,3),xji(3,3),xj2D(2,2),xji2D(2,2),r(3,3)
      real*8 ri,si,ti,w,wt1,wt2,wt3,det
      real*8 rn(8),sn(8),tn(8)
      real*8 pi,epsi(6),txi(6),txn(6,*)
c ... integracao numerica de tetraedros          
      real*8 pg(10,10),wg(10,10)
c ...
      real*8 e(*),x(ndm,*),xf(3,4)
c ...
      real*8 a,b,c,ym,ps,density
      real*8 gl(3)
      real*8 a1,a2,a3
c ... 
      integer hexa_face_node8(4,6),no
      real*8 hexa_vol,volum
c ...
      data hexa_face_node8 / 1, 2, 3, 4
     .                     , 5, 6, 7, 8
     .                     , 1, 5, 6, 2
     .                     , 4, 3, 7, 8
     .                     , 1, 4, 8, 5
     .                     , 2, 6, 7, 3/
c
      data rn / 1.d0,-1.d0,-1.d0, 1.d0, ! r1, r2, r3, r4
     .          1.d0,-1.d0,-1.d0, 1.d0/ ! r5, r6, r7, r8              
     
c
      data sn / 1.d0, 1.d0,-1.d0,-1.d0, ! s1, s2, s3, s4           
     .          1.d0, 1.d0,-1.d0,-1.d0/ ! s5, s6, s7, s8             
     
c
      data tn / 1.d0, 1.d0, 1.d0, 1.d0, ! t1, t2, t3, t4
     .         -1.d0,-1.d0,-1.d0,-1.d0/ ! t5, t6, t7, t8
      
c
      data nen/8/
c ......................................................................
      goto (100,200,300,400,500) isw
c ======================================================================
c
c ...                                          
c
c ......................................................................
  100 continue
c ...
c .....................................................................
      return
c ======================================================================
c
c ... Matriz de rigidez:
c
c ......................................................................
  200 continue
c ... matriz constitutiva:
      ym       = e(1)
      ps       = e(2)
c ...
      a1       = 1.d0 - ps
      a2       = 1.d0 - 2.d0*ps
      a3       = 1.d0 + ps
c ...
      a        = (ym*a1)/(a3*a2)
      b        = ps/a1
      c        = 0.5d0*(a2/a1) 
c .....................................................................
c
c ... Matriz de rigidez:
      do i = 1, nst
        do j = 1, nst
          s(j,i) = 0.d0
        enddo
      enddo
c .....................................................................
c
c ... 
      nint = 2 
      do 205 lz = 1, nint
        ti = pg(lz,nint)
        do 210 ly = 1, nint
          si = pg(ly,nint)
          do 215 lx = 1, nint
            ri = pg(lx,nint)
c ...                                               
            call sfhexa8_m(h,hx,hy,hz,ri,si,ti,.false.,.true.)
            call jacob3d_m(hx,hy,hz,xj,xji,x,det,8,nel,.true.)
c .....................................................................
c
c ...
            w   = wg(lx,nint)*wg(ly,nint)*wg(lz,nint)*det
c .....................................................................
c
c ... Kuu ( Int((Bt)*C*(B)*dV) )
            wt1 = w*a
            do 220 i = 1, 8
c             l1 = (i-1)*3+1
              l1 = 3*i-2
              l2 = l1 + 1
              l3 = l2 + 1
              do 225 j = 1, 8
c               k1 = (j-1)*3+1
                k1 = 3*j-2
                k2 = k1 + 1
                k3 = k2 + 1
                s(l1,k1) = s(l1,k1) +
     .          (hx(i)*hx(j) + c*(hy(i)*hy(j) + hz(i)*hz(j)))*wt1
c
                s(l1,k2) = s(l1,k2) +
     .          (b*hx(i)*hy(j) + c*hy(i)*hx(j))*wt1
c
                s(l1,k3) = s(l1,k3) + 
     .          (b*hx(i)*hz(j) + c*hz(i)*hx(j))*wt1
c
                s(l2,k1) = s(l2,k1) + 
     .          (b*hy(i)*hx(j) + c*hx(i)*hy(j))*wt1
c
                s(l2,k2) = s(l2,k2) + 
     .          (hy(i)*hy(j) + c*(hx(i)*hx(j) + hz(i)*hz(j)))*wt1
c
                s(l2,k3) = s(l2,k3) +
     .          (b*hy(i)*hz(j) + c*hz(i)*hy(j))*wt1
c
                s(l3,k1) = s(l3,k1) +
     .          (b*hz(i)*hx(j) + c*hx(i)*hz(j))*wt1
c
                s(l3,k2) = s(l3,k2) + 
     .          (b*hz(i)*hy(j) + c*hy(i)*hz(j))*wt1
c
                s(l3,k3) = s(l3,k3) +
     .          (hz(i)*hz(j) + c*(hy(i)*hy(j) + hx(i)*hx(j)))*wt1
c .....................................................................                   
  225         continue
c .....................................................................              
  220       continue 
c .....................................................................
  215     continue
c .....................................................................
  210   continue
c .....................................................................
  205 continue
c .....................................................................
c
c ... Forcas internas:
      call lku_m(s,u,p,nst)
c .....................................................................
c
c .....................................................................      
      return  
c ======================================================================
c
c ... Tensoes nodais :
c
c ......................................................................
  300 continue
c ... matriz constitutiva:
      ym        = e(1)
      ps        = e(2)
c ...
      a1        = 1.d0 - ps
      a2        = 1.d0 - 2.d0*ps
      a3        = 1.d0 + ps
c ...
      a         = (ym*a1)/(a3*a2)
      b         = ps/a1
      c         = 0.5d0*(a2/a1)
c .....................................................................
c
c ... tensao nodal total
      do 310 i = 1, 8
c       tp = (i-1)*6 + 1
        tp  = 6*i - 5
c ... calculo do terminante
        call sfhexa8_m(h,hx,hy,hz,rn(i),sn(i),tn(i),.false.,.true.)
        call jacob3d_m(hx,hy,hz,xj,xji,x,det,8,nel ,.true.)
c .....................................................................
        call deform3d(hx,hy,hz,u,epsi,8)
        call stress3d(a,b,c,epsi,p(tp))
c .....................................................................
  310 continue
c .....................................................................
c
c ......................................................................
      return
c ======================================================================
c
c ... Cargas distribuidas no volume e no contorno:
c
c ......................................................................
 400  continue
c ... peso proprio  
      if( density .eq. 0.0d0) goto 430
      density   = e(3)*1.0d-06 
c .....................................................................
c
c ... 
      gl(1)     =  gravity(1)*density
      gl(2)     =  gravity(2)*density
      gl(3)     =  gravity(3)*density
c .....................................................................
c
c ...                            
      nint = 2 
      do 405 lz = 1, nint
        ti = pg(lz,nint)
        do 410 ly = 1, nint
          si = pg(ly,nint)
          do 415 lx = 1, nint
            ri = pg(lx,nint)
c ...                                               
            call sfhexa8_m(h,hx,hy,hz,ri,si,ti,.true.,.true.)
            call jacob3d_m(hx,hy,hz,xj,xji,x,det,8,nel,.true.)
c .....................................................................
c
c ...
            w   = wg(lx,nint)*wg(ly,nint)*wg(lz,nint)*det
c .....................................................................
c
c ...
            do 420 i = 1, 8
c             l1  = (i-1)*3+1
              l1  = 3*i-2
              l2  = l1 + 1
              l3  = l2 + 1
c ... Fu = int(pm_d*(N~T)*ge*dv)
              wt1   = h(i)*w  
              p(l1) = p(l1) - wt1*gl(1)
              p(l2) = p(l2) - wt1*gl(2)
              p(l3) = p(l3) - wt1*gl(3)
  420       continue
c .....................................................................    
  415     continue
c .....................................................................    
  410   continue
c .....................................................................    
  405 continue
c .....................................................................
c
c .....................................................................
c
 430  continue 
c ... forca e fluxo distribuida no contorno
c     iq(1) = 1 | no 1 2 3 4 |
c             2 | no 5 6 7 8 |
c             3 | no 1 5 6 2 |     
c             4 | no 4 3 7 8 |
c             5 | no 1 4 8 5 |
c             6 | no 2 6 7 3 |
c
c ... verifica se ha alguma face com carga
      tp = 0
      do 435 i = 1, 6
        tp    = tp + iq(i)
  435 continue
c .....................................................................
c
c ...
      if( tp .gt. 0 ) then
        do 440 j = 1, 6
c ... face 
          if(iq(j) .gt. 0 ) then
c ...
            do 445 i = 1, 4
              no       = hexa_face_node8(i,j)
              xf(1,i) = x(1,no) 
              xf(2,i) = x(2,no) 
              xf(3,i) = x(3,no)
  445       continue
c ...    
            call rotmatrix(xf,r)
            call rotate(xf(1,1),r,xf(1,1),.false.)
            call rotate(xf(1,2),r,xf(1,2),.false.)
            call rotate(xf(1,3),r,xf(1,3),.false.)
            call rotate(xf(1,4),r,xf(1,4),.false.)
c ... forcas
            call tload(iq(j),0.d0,0.d0,ddum,face_f)
            nint = 2
            do 450 ly = 1, nint
              si = pg(ly,nint)
              do 455 lx = 1, nint
                ri = pg(lx,nint)
c ...                                               
                call sfquad4_m(h,hx,hy,ri,si,.true.,.true.)
                call jacob2d_m(hx,hy,xj2D,xji2D,xf,det,4,ndm
     .                        ,nel,.true.)
c .....................................................................
c
c ...                                               
                w   = wg(lx,nint)*wg(ly,nint)*det
c .....................................................................
c
c ...
                do 460 i = 1, 4
                  no  = hexa_face_node8(i,j)
c                   l1  = (no-1)*3+1
                  wt1   = h(i)*w
                  l1    = 3*no-2
                  l2    = l1 + 1
                  l3    = l2 + 1
                  p(l1) = p(l1) - face_f(1)*wt1
                  p(l2) = p(l2) - face_f(2)*wt1
                  p(l3) = p(l3) - face_f(3)*wt1  
  460           continue
c .....................................................................
  455         continue
c .....................................................................             
  450       continue
c .....................................................................
          endif
c .....................................................................
  440   continue
c .....................................................................
      endif
c .....................................................................
      return
c ======================================================================
c
c ... Tensoes iniciais:                  
c
c ......................................................................
  500 continue
      nint = 4 
      do 510 lz = 1, nint
        ti = pg(lz,nint)
        do 520 ly = 1, nint
          si = pg(ly,nint)
          do 530 lx = 1, nint
            ri = pg(lx,nint)
c ...                                               
            call sfhexa8_m(h,hx,hy,hz,ri,si,ti,.true.,.true.)
            call jacob3d_m(hx,hy,hz,xj,xji,x,det,8,nel,.true.)
c .....................................................................
c
c ...
            w   = wg(lx,nint)*wg(ly,nint)*wg(lz,nint)*det
            wt1 = w  
c .....................................................................
c
c ...
            txi(1:6)= 0.d0
            do 540 j = 1,   8
              txi(1) = txi(1) + h(j)*txn(1,j)
              txi(2) = txi(2) + h(j)*txn(2,j) 
              txi(3) = txi(3) + h(j)*txn(3,j) 
              txi(4) = txi(4) + h(j)*txn(4,j) 
              txi(5) = txi(5) + h(j)*txn(5,j)
              txi(6) = txi(6) + h(j)*txn(6,j) 
  540       continue 
c ..................................................................... 
c
c ...
            do 560 i = 1, 8
c             l1  = (i-1)*3+1
              l1  = 3*i-2
              l2  = l1 + 1
              l3  = l2 + 1
c ... Fu = int(BeT*sigma*dv)
              p(l1) = p(l1) 
     .        + wt1*(hx(i)*txi(1) + hy(i)*txi(4) + hz(i)*txi(6))
c
              p(l2) = p(l2) 
     .        + wt1*(hx(i)*txi(4) + hy(i)*txi(2) + hz(i)*txi(5))
c
               p(l3) = p(l3) 
     .        + wt1*(hx(i)*txi(6) + hy(i)*txi(5) + hz(i)*txi(3))
c ...................................................................... 
  560       continue  
c .....................................................................
  530     continue
c .....................................................................
  520   continue
c .....................................................................
  510 continue
c .....................................................................
c
c ....
      return
c ======================================================================
c
c ... Matriz geometrica:
c
c ......................................................................
  600 continue
      stop
c ......................................................................
      end
c **********************************************************************
      subroutine elmt12_mec(e,iq,x,u,p,s,txn,ndm,nst,nel,isw)
c **********************************************************************
c * Data de criacao    : 09/04/2016                                    *
c * Data de modificaco : 00/00/0000                                    * 
c * ------------------------------------------------------------------ *      
c * ELMT12_MEC : Elemento tetraedrico de 10 nos para problemas         *  
c * mecanico elasticos                                                 *
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ * 
c * e(10) - constantes fisicas                                         *
c *           e(1) = modulo de Young                                   *
c *           e(2) = coeficiente de Poisson                            *
c *           e(3) = massa especifica do meio poroso                   *
c * x(ndm,nem) - coordenadas nodais locais                             *
c * u(nst)     - graus de liberade por elemento (u)                    *
c * p(nst)     - nao definido                                          *
c * s(nst,nst) - nao definido                                          *
c * txn(6,nen) - tensoes nodais                                        *
c * ndm - dimensao                                                     *
c * nst - numero de graus de liberdade por elemento (3*10)             *
c * nel - numero do elemento                                           *
c * isw - codigo de instrucao                                          *
c *  1 =                                                               *
c *  2 = matriz K e forcas internas Ku                                 *
c *  3 = tensoes                                                       *
c *  4 = forcas de volume e superficies                                *
c *  5 = Tensoes iniciais                                              *
c *  6 =                                                               *
c *  7 =                                                               *
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ * 
c * e - constantes fisicas                                             *
c * s - matriz de elemento                                             *
c * p - isw = 2  residuo                                               *
c *     isw = 3  tensao                                                *
c *     isw = 4  cargas de superfice e volume                          *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ * 
c **********************************************************************
      implicit none
      include 'gravity.fi'
      include 'load.fi'
      common /pint4/ pri4,psi4,pti4,wf4,npint4
      common /pint / pri,psi,wf,npint
      real*8  div6         
      parameter ( div6 = 0.166666666666667d0)
      integer ndm,nst,nel,isw
      integer i,j,l1,l2,l3,l,k1,k2,k3,k,tp,tp1
      integer nen,nint,lx,ly,lz
      integer iq(*)
c ...
      real*8 face_f(3),ddum
c ...
      real*8 u(*)
      real*8 p(*),s(nst,*)
c ... funcoes de interpolacao
      real*8 hu(10),hux(10),huy(10),huz(10)
      real*8 h(4),hx(4),hy(4),hz(4)
c ...
      real*8 xj(3,3),xji(3,3),xj2D(2,2),xji2D(2,2),r(3,3)
      real*8 ri,si,ti,w,wt1,det
      real*8 rn(10),sn(10),tn(10)
      real*8 pi,epsi(6),txi(6),txn(6,*),pm
c ... integracao numerica de tetraedros     
      real*8 pri4(5,3),psi4(5,3),pti4(5,3),wf4(5,3)
      real*8  pri(12,5),psi(12,5),wf(12,5)
      integer igrau,igrau_vol,igrau_face,npint(5),npint4(3)
c ...
      real*8 e(*),x(ndm,*),xf(3,3)
c ...
      real*8 a,b,c,ym,ps
      real*8 density,gl(3)
      real*8 a1,a2,a3,tmp
c ... 
      integer tetra_face_node10(6,4),no
      real*8 tetra_vol,volum
c ...
      data tetra_face_node10 / 2, 3, 4, 8, 9,10
     .                       , 1, 4, 3, 7, 9, 6
     .                       , 1, 2, 4, 5,10, 7
     .                       , 1, 3, 2, 6, 8, 5/
c
      data rn / 1.0d0, 0.0d0, 0.0d0 ! r1, r2, r3
     .       ,  0.0d0, 0.5d0, 0.5d0 ! r4, r5, r6              
     .       ,  0.5d0, 0.0d0, 0.0d0 ! r7, r8, r9              
     .       ,  0.0d0/              ! r10              
c
      data sn / 0.0d0, 1.0d0, 0.0d0 ! s1, s2, s3
     .        , 0.0d0, 0.5d0, 0.0d0 ! s4, s5, s6              
     .        , 0.0d0, 0.5d0, 0.0d0 ! s7, s8, s9              
     .        , 0.5d0/              ! s10              
c
      data tn / 0.0d0, 0.0d0, 0.0d0 ! t1, t2, t3
     .        , 1.0d0, 0.0d0, 0.0d0 ! t4, t5, t6              
     .        , 0.5d0, 0.0d0, 0.0d0 ! t7, t8, t9              
     .        , 0.5d0/              ! t10              
c
      data nen/10/,igrau_vol/2/,igrau_face/2/    
c ......................................................................
      goto (100,200,300,400,500) isw
c ======================================================================
c
c.... calculo do delta t critico               
c
c ......................................................................
  100 continue
c .....................................................................
      return
c ======================================================================
c
c ... Matriz de rigidez:
c
c ......................................................................
  200 continue
c ... matriz constitutiva:
      ym       = e(1)
      ps       = e(2)
c ...
      a1       = 1.d0 - ps
      a2       = 1.d0 - 2.d0*ps
      a3       = 1.d0 + ps
c ...
      a        = (ym*a1)/(a3*a2)
      b        = ps/a1
      c        = 0.5d0*(a2/a1) 
c .....................................................................
c
c ... Matriz de rigidez:
      do i = 1, nst
        do j = 1, nst
          s(j,i) = 0.d0
        enddo
      enddo
c .....................................................................
c
c ... 
      igrau = igrau_vol 
      nint = npint4(igrau) 
      do 205 lx = 1, nint
c ...
        ti = pti4(lx,igrau)
        si = psi4(lx,igrau)
        ri = pri4(lx,igrau)
c ...                                               
        call sftetra4(h,hx,hy,hz,ri,si,ti,.false.,.true.)
        call jacob3d_m(hx,hy,hz,xj,xji,x,det,4,nel,.true.)
c .....................................................................
c
c ...
        w   = wf4(lx,igrau)*det
        w   = w*div6
c .....................................................................
c
c ...
        call sftetra10(hu,hux,huy,huz,ri,si,ti,.false.,.true.)
        call jacob3d_m(hux,huy,huz,xj,xji,x,det,10,nel,.false.)
c .....................................................................
c
c ... Kuu ( Int((Bt)*C*(B)*dV) )
        wt1 = w*a
        do 210 i = 1, 10
c         l1 = (i-1)*3+1
          l1 = 3*i-2
          l2 = l1 + 1
          l3 = l2 + 1
          do 215 j = 1, 10
c           k1 = (j-1)*3+1
            k1 = 3*j-2
            k2 = k1 + 1
            k3 = k2 + 1
            s(l1,k1) = s(l1,k1) +
     .      (hux(i)*hux(j) + c*(huy(i)*huy(j) + huz(i)*huz(j)))*wt1
c
            s(l1,k2) = s(l1,k2) +
     .      (b*hux(i)*huy(j) + c*huy(i)*hux(j))*wt1
c
            s(l1,k3) = s(l1,k3) + 
     .      (b*hux(i)*huz(j) + c*huz(i)*hux(j))*wt1
c
            s(l2,k1) = s(l2,k1) + 
     .      (b*huy(i)*hux(j) + c*hux(i)*huy(j))*wt1
c
            s(l2,k2) = s(l2,k2) + 
     .      (huy(i)*huy(j) + c*(hux(i)*hux(j) + huz(i)*huz(j)))*wt1
c
            s(l2,k3) = s(l2,k3) +
     .      (b*huy(i)*huz(j) + c*huz(i)*huy(j))*wt1
c
            s(l3,k1) = s(l3,k1) +
     .      (b*huz(i)*hux(j) + c*hux(i)*huz(j))*wt1
c
            s(l3,k2) = s(l3,k2) + 
     .      (b*huz(i)*huy(j) + c*huy(i)*huz(j))*wt1
c
            s(l3,k3) = s(l3,k3) +
     .      (huz(i)*huz(j) + c*(huy(i)*huy(j) + hux(i)*hux(j)))*wt1
c .....................................................................            
  215     continue
c .....................................................................     
  210   continue 
c .....................................................................
 205  continue 
c .....................................................................
c
c ... Forcas internas:
      call lku_m(s,u,p,nst)
c     call lku_sym(s,u,p,nst)    
c ......................................................................
c
c ......................................................................   
      return  
c ======================================================================
c
c ... Tensoes nodais:
c
c ......................................................................
  300 continue
c ... matriz constitutiva:
      ym        = e(1)
      ps        = e(2)
c ...
      a1        = 1.d0 - ps
      a2        = 1.d0 - 2.d0*ps
      a3        = 1.d0 + ps
c ...
      a         = (ym*a1)/(a3*a2)
      b         = ps/a1
      c         = 0.5d0*(a2/a1)
c .....................................................................
c
c ... tensao nodal total
      do 310 i = 1, 4
c       tp = (i-1)*6 + 1
        tp  = 6*i - 5
c ... calculo do determinante
        call sftetra4(h,hx,hy,hz,rn(i),sn(i),tn(i),.false.,.true.)
        call jacob3d_m(hx,hy,hz,xj,xji,x,det,4,nel ,.true.)
c ... calculo da derivadas das funcoes de interpolacao
        call sftetra10(hu,hux,huy,huz,rn(i),sn(i),tn(i),.false.,.true.)
        call jacob3d_m(hux,huy,huz,xj,xji,x,det,10,nel ,.false.)
c .....................................................................
        call deform3d(hux,huy,huz,u,epsi,10)
        call stress3d(a,b,c,epsi,p(tp))
c .....................................................................
  310 continue
c .....................................................................
c
c ......................................................................
      return
c ======================================================================
c
c ... Cargas distribuidas no volume e no contorno:
c
c ......................................................................
 400  continue
c ... peso proprio
      if( density .eq. 0.0d0) goto 430
      density   = e(3)*1.0d-06 
c .....................................................................
c
c ... 
      gl(1)     =  gravity(1)*density
      gl(2)     =  gravity(2)*density
      gl(3)     =  gravity(3)*density
c .....................................................................
c
c ...     
      igrau = igrau_vol 
      nint  = npint4(igrau) 
      do 410 lx = 1, nint
c ...
        ti = pti4(lx,igrau)
        si = psi4(lx,igrau)
        ri = pri4(lx,igrau)
c ...                                               
        call sftetra4(h,hx,hy,hz,ri,si,ti,.false.,.true.)
        call jacob3d_m(hx,hy,hz,xj,xji,x,det,4,nel,.true.)
c .....................................................................
c
c ...
        w   = wf4(lx,igrau)*det
        w   = w*div6
c .....................................................................
c
c ...
       call sftetra10(hu,hux,huy,huz,ri,si,ti,.true.,.false.)
c .....................................................................
c
c ...
        do 420 i = 1, 10
c         l1  = (i-1)*3+1
          l1  = 3*i-2
          l2  = l1 + 1
          l3  = l2 + 1
c ... Fu = int(pm_d*(N~T)*ge*dv)  
          wt1   = hu(i)*w  
          p(l1) = p(l1) - wt1*gl(1)
          p(l2) = p(l2) - wt1*gl(2)
          p(l3) = p(l3) - wt1*gl(3)
  420   continue    
c .....................................................................
c
 410  continue   
c .....................................................................
c
 430  continue  
c ... forca e fluxo distribuida no contorno
c     iq(1) = 1 | no 2 3 4  8  9 10 |
c             2 | no 1 4 3  7  9  6 |
c             3 | no 1 2 4  5 10  7 |     
c             4 | no 1 3 2  6  8  5 |
c ... verifica se ha alguma face com carga      
      tp = 0
      do 435 i = 1, 4
        tp  = tp + iq(i)
  435 continue
c .....................................................................
c
c ...      
      if( tp .gt. 0 ) then
        do 440 j = 1, 4
c ... face  
          if(iq(j) .gt. 0 ) then
c ...
            do 445 i = 1, 3
              no       = tetra_face_node10(i,j)
              xf(1,i) = x(1,no) 
              xf(2,i) = x(2,no) 
              xf(3,i) = x(3,no)
  445      continue
c ...    
            call rotmatrix(xf,r)
            call rotate(xf(1,1),r,xf(1,1),.false.)
            call rotate(xf(1,2),r,xf(1,2),.false.)
            call rotate(xf(1,3),r,xf(1,3),.false.)
c ... forcas
            call tload(iq(j),0.d0,0.d0,ddum,face_f)
            igrau = igrau_face 
            nint = npint(igrau)
            do 450 lx = 1, nint
              ri = pri(lx,igrau)
              si = psi(lx,igrau)
c ...                                 
              call sftria3(h,hx,hy,ri,si,.false.,.true.)
              call jacob2d_m(hx,hy,xj2D,xji2D,xf,det,3,ndm
     .                        ,nel,.true.)
c .....................................................................
c
c ...                                               
              w   = wf(lx,igrau)*det
              w   = 0.5d0*w
c .....................................................................
c
c ...                   
              call sftria6(hu,hux,huy,ri,si,.true.,.false.)
c .....................................................................
c
c ...
              do 455 i = 1, 6 
                no  = tetra_face_node10(i,j)
c               l1  = (no-1)*3+1
                wt1   = hu(i)*w
                l1    = 3*no-2
                l2    = l1 + 1
                l3    = l2 + 1
                p(l1) = p(l1) - face_f(1)*wt1
                p(l2) = p(l2) - face_f(2)*wt1
                p(l3) = p(l3) - face_f(3)*wt1
  455         continue
c .....................................................................
  450       continue
c .....................................................................
          endif
c .....................................................................
  440   continue
c .....................................................................
      endif
c .....................................................................
      return
c ======================================================================
c
c ... Tensoes iniciais:                  
c
c ......................................................................
  500 continue
c ...     
      igrau = igrau_vol 
      nint  = npint4(igrau) 
      do 510 lx = 1, nint
c ...
        ti = pti4(lx,igrau)
        si = psi4(lx,igrau)
        ri = pri4(lx,igrau)
c ...                                               
        call sftetra4(h,hx,hy,hz,ri,si,ti,.true.,.true.)
        call jacob3d_m(hx,hy,hz,xj,xji,x,det,4,nel,.true.)
c .....................................................................
c
c ...
        w   = wf4(lx,igrau)*det
        wt1 = w*div6
c .....................................................................
c
c ...
        call sftetra10(hu,hux,huy,huz,ri,si,ti,.false.,.true.)
        call jacob3d_m(hux,huy,huz,xj,xji,x,det,10,nel,.false.)
c .....................................................................
c  
c ...
        txi(1:6)= 0.d0
        do 515 j = 1,   4
          txi(1) = txi(1) + h(j)*txn(1,j)
          txi(2) = txi(2) + h(j)*txn(2,j) 
          txi(3) = txi(3) + h(j)*txn(3,j) 
          txi(4) = txi(4) + h(j)*txn(4,j) 
          txi(5) = txi(5) + h(j)*txn(5,j)
          txi(6) = txi(6) + h(j)*txn(6,j) 
  515   continue
c .....................................................................
c        
c ...
        do 520 i = 1, 10
c         l1  = (i-1)*3+1
          l1  = 3*i-2
          l2  = l1 + 1
          l3  = l2 + 1
c ... Fu = int(BeT*sigma*dv)
          p(l1) = p(l1) 
     .    + wt1*(hux(i)*txi(1) + huy(i)*txi(4) + huz(i)*txi(6))
c         
          p(l2) = p(l2) 
     .    + wt1*(hux(i)*txi(4) + huy(i)*txi(2) + huz(i)*txi(5))
c          
          p(l3) = p(l3) 
     .    + wt1*(hux(i)*txi(6) + huy(i)*txi(5) + huz(i)*txi(3))
c .....................................................................
  520   continue  
c .....................................................................
  510 continue 
c .....................................................................
      return
c ======================================================================
c
c ... Matriz geometrica:
c
c ......................................................................
  600 continue
      stop
c ....................................................................
      end
c *********************************************************************
      subroutine elmt13_mec(e,iq,x,u,p,s,txn,ndm,nst,nel,isw)
c **********************************************************************
c * Data de criacao    : 09/04/2016                                    *
c * Data de modificaco : 00/00/0000                                    * 
c * ------------------------------------------------------------------ *       
c * ELMT013_MEC: Elemento hexaedricos de 20 nos para problemas         *  
c * mecanico estaico-elasticos                                         *
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ * 
c * e(10) - constantes fisicas                                         *
c *           e(1) = modulo de Young                                   *
c *           e(2) = coeficiente de Poisson                            *
c *           e(3) = massa especifica do meio poroso                   *
c *           e(4) =                                                   *
c *           e(5) =                                                   *
c *           e(6) =                                                   *
c *           e(7) =                                                   *
c * x(ndm,nem) - coordenadas nodais locais                             *
c * u(nst)     - graus de liberade por elemento (u)                    *
c * p(nst)     - nao definido                                          *
c * s(nst,nst) - nao definido                                          *
c * txn(6,nen) - tensoes nodais                                        *      
c * ndm - dimensao                                                     *
c * nst - numero de graus de liberdade por elemento (3*20)             *
c * nel - numero do elemento                                           *
c * isw - codigo de instrucao                                          *
c *  1 =                                                               *
c *  2 = matriz K e forcas internas Ku                                 *
c *  3 = tensoes                                                       *
c *  4 = forcas de volume e superficies                                *
c *  5 = Tensoes iniciais                                              *
c *  6 =                                                               *
c *  7 =                                                               *
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ * 
c * e - constantes fisicas                                             *
c * s - matriz de elemento                                             *
c * p - isw = 2  residuo                                               *
c *     isw = 3  tensao e fluxo                                        *
c *     isw = 4  cargas de superfice e volume                          *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ * 
c **********************************************************************
      implicit none
      include 'gravity.fi'
      include 'load.fi'
      common /gauss/ pg, wg
      integer ndm,nst,nel,isw
      integer i,j,l1,l2,l3,l,k1,k2,k3,k,tp,tp1
      integer nen,nint,lx,ly,lz
      integer iq(*)
c ...
      real*8 face_f(3),ddum
c ...
      real*8 u(*)
      real*8 p(*),s(nst,*)
c ... funcoes de interpolacao
      real*8 hu(20),hux(20),huy(20),huz(20)
      real*8 h(8),hx(8),hy(8),hz(8)
c ...
      real*8 xj(3,3),xji(3,3),xj2D(2,2),xji2D(2,2),r(3,3)
      real*8 ri,si,ti,w,wt1,wt2,wt3,det
      real*8 rn(20),sn(20),tn(20)
      real*8 pi,epsi(6),txi(6),txn(6,*)
c ... integracao numerica de tetraedros          
      real*8 pg(10,10),wg(10,10)
c ...
      real*8 e(*),x(ndm,*),xf(3,4)
c ...
      real*8 a,b,c,ym,ps,density
      real*8 gl(3)
      real*8 a1,a2,a3
c ... 
      integer hexa_face_node20(8,6),no
      real*8 hexa_vol,volum
c ...
      data hexa_face_node20 / 1, 2, 3, 4, 9,10,11,12
     .                      , 5, 6, 7, 8,13,14,15,16
     .                      , 1, 5, 6, 2,17,13,18, 9
     .                      , 4, 3, 7, 8,11,19,15,20
     .                      , 1, 4, 8, 5,12,20,16,17
     .                      , 2, 6, 7, 3,18,14,19,10/
c
      data rn / 1.d0,-1.d0,-1.d0, 1.d0, ! r1, r2, r3, r4
     .          1.d0,-1.d0,-1.d0, 1.d0, ! r5, r6, r7, r8              
     .          0.d0,-1.d0, 0.d0, 1.d0, ! r9,r10,r11,r12              
     .          0.d0,-1.d0, 0.d0, 1.d0, !r13,r14,r15,r16               
     .          1.d0,-1.d0,-1.d0, 1.d0/ !r17,r18,r19,r20  
c
      data sn / 1.d0, 1.d0,-1.d0,-1.d0, ! s1, s2, s3, s4           
     .          1.d0, 1.d0,-1.d0,-1.d0, ! s5, s6, s7, s8             
     .          1.d0, 0.d0,-1.d0, 0.d0, ! s9,s10,s11,s12             
     .          1.d0, 0.d0,-1.d0, 0.d0, !s13,s14,s15,s16             
     .          1.d0, 1.d0,-1.d0,-1.d0/ !s17,s18,s19,s20  
c
      data tn / 1.d0, 1.d0, 1.d0, 1.d0, ! t1, t2, t3, t4
     .         -1.d0,-1.d0,-1.d0,-1.d0, ! t5, t6, t7, t8
     .          1.d0, 1.d0, 1.d0, 1.d0, ! t9,t10,t11,t12
     .         -1.d0,-1.d0,-1.d0,-1.d0, !t13,t14,t15,t16      
     .          0.d0, 0.d0, 0.d0, 0.d0/ !t17,t18,t19,t20       
c
      data nen/20/
c ......................................................................
      goto (100,200,300,400,500) isw
c ======================================================================
c
c ...                                          
c
c ......................................................................
  100 continue
c ...
c .....................................................................
      return
c ======================================================================
c
c ... Matriz de rigidez:
c
c ......................................................................
  200 continue
c ... matriz constitutiva:
      ym       = e(1)
      ps       = e(2)
c ...
      a1       = 1.d0 - ps
      a2       = 1.d0 - 2.d0*ps
      a3       = 1.d0 + ps
c ...
      a        = (ym*a1)/(a3*a2)
      b        = ps/a1
      c        = 0.5d0*(a2/a1) 
c .....................................................................
c
c ... Matriz de rigidez:
      do i = 1, nst
        do j = 1, nst
          s(j,i) = 0.d0
        enddo
      enddo
c .....................................................................
c
c ... 
      nint = 4 
      do 205 lz = 1, nint
        ti = pg(lz,nint)
        do 210 ly = 1, nint
          si = pg(ly,nint)
          do 215 lx = 1, nint
            ri = pg(lx,nint)
c ...                                               
            call sfhexa8_m(h,hx,hy,hz,ri,si,ti,.true.,.true.)
            call jacob3d_m(hx,hy,hz,xj,xji,x,det,8,nel,.true.)
c .....................................................................
c
c ...
            w   = wg(lx,nint)*wg(ly,nint)*wg(lz,nint)*det
c .....................................................................
c
c ...
            call sfhexa20_m(hu,hux,huy,huz,ri,si,ti,.false.,.true.)
            call jacob3d_m(hux,huy,huz,xj,xji,x,det,20,nel,.false.)
c .....................................................................
c
c ... Kuu ( Int((Bt)*C*(B)*dV) )
            wt1 = w*a
            do 220 i = 1, 20
c             l1 = (i-1)*3+1
              l1 = 3*i-2
              l2 = l1 + 1
              l3 = l2 + 1
              do 225 j = 1, 20
c               k1 = (j-1)*3+1
                k1 = 3*j-2
                k2 = k1 + 1
                k3 = k2 + 1
                s(l1,k1) = s(l1,k1) +
     .          (hux(i)*hux(j) + c*(huy(i)*huy(j) + huz(i)*huz(j)))*wt1
c
                s(l1,k2) = s(l1,k2) +
     .          (b*hux(i)*huy(j) + c*huy(i)*hux(j))*wt1
c
                s(l1,k3) = s(l1,k3) + 
     .          (b*hux(i)*huz(j) + c*huz(i)*hux(j))*wt1
c
                s(l2,k1) = s(l2,k1) + 
     .          (b*huy(i)*hux(j) + c*hux(i)*huy(j))*wt1
c
                s(l2,k2) = s(l2,k2) + 
     .          (huy(i)*huy(j) + c*(hux(i)*hux(j) + huz(i)*huz(j)))*wt1
c
                s(l2,k3) = s(l2,k3) +
     .          (b*huy(i)*huz(j) + c*huz(i)*huy(j))*wt1
c
                s(l3,k1) = s(l3,k1) +
     .          (b*huz(i)*hux(j) + c*hux(i)*huz(j))*wt1
c
                s(l3,k2) = s(l3,k2) + 
     .          (b*huz(i)*huy(j) + c*huy(i)*huz(j))*wt1
c
                s(l3,k3) = s(l3,k3) +
     .          (huz(i)*huz(j) + c*(huy(i)*huy(j) + hux(i)*hux(j)))*wt1
c .....................................................................                   
  225         continue
c .....................................................................              
  220       continue 
c .....................................................................
  215     continue
c .....................................................................
  210   continue
c .....................................................................
  205 continue
c .....................................................................
c
c ... Forcas internas:
      call lku_m(s,u,p,nst)
c     call lku_sym(s,u,p,nst)    
c .....................................................................
c
c .....................................................................      
      return  
c ======================================================================
c
c ... Tensoes nodais :
c
c ......................................................................
  300 continue
c ... matriz constitutiva:
      ym        = e(1)
      ps        = e(2)
c ...
      a1        = 1.d0 - ps
      a2        = 1.d0 - 2.d0*ps
      a3        = 1.d0 + ps
c ...
      a         = (ym*a1)/(a3*a2)
      b         = ps/a1
      c         = 0.5d0*(a2/a1)
c .....................................................................
c
c ... tensao nodal total
      do 310 i = 1, 8
c       tp = (i-1)*6 + 1
        tp  = 6*i - 5
c ... calculo do terminante
        call sfhexa8_m(h,hx,hy,hz,rn(i),sn(i),tn(i),.false.,.true.)
        call jacob3d_m(hx,hy,hz,xj,xji,x,det,8,nel ,.true.)
c ... calculo da derivadas das funcoes de interpolacao
        call sfhexa20_m(hu,hux,huy,huz,rn(i),sn(i),tn(i),.false.,.true.)
        call jacob3d_m(hux,huy,huz,xj,xji,x,det,20,nel ,.false.)
c .....................................................................
        call deform3d(hux,huy,huz,u,epsi,20)
        call stress3d(a,b,c,epsi,p(tp))
c .....................................................................
  310 continue
c .....................................................................
c
c ......................................................................
      return
c ======================================================================
c
c ... Cargas distribuidas no volume e no contorno:
c
c ......................................................................
 400  continue
c ... peso proprio
      if( density .eq. 0.0d0) goto 430
      density   = e(3)*1.0d-06 
c .....................................................................
c
c ... 
      gl(1)     =  gravity(1)*density
      gl(2)     =  gravity(2)*density
      gl(3)     =  gravity(3)*density
c .....................................................................
c
c ...                            
      nint = 4 
      do 405 lz = 1, nint
        ti = pg(lz,nint)
        do 410 ly = 1, nint
          si = pg(ly,nint)
          do 415 lx = 1, nint
            ri = pg(lx,nint)
c ...                                               
            call sfhexa8_m(h,hx,hy,hz,ri,si,ti,.false.,.true.)
            call jacob3d_m(hx,hy,hz,xj,xji,x,det,8,nel,.true.)
c .....................................................................
c
c ...
            w   = wg(lx,nint)*wg(ly,nint)*wg(lz,nint)*det
c .....................................................................
c
c ...
            call sfhexa20_m(hu,hux,huy,huz,ri,si,ti,.true.,.false.)
c .....................................................................
c
c ...
            do 420 i = 1, 20
c             l1  = (i-1)*3+1
              l1  = 3*i-2
              l2  = l1 + 1
              l3  = l2 + 1
c ... Fu = int(pm_d*(N~T)*ge*dv)
              wt1   = hu(i)*w  
              p(l1) = p(l1) - wt1*gl(1)
              p(l2) = p(l2) - wt1*gl(2)
              p(l3) = p(l3) - wt1*gl(3)
  420       continue
c .....................................................................    
  415     continue
c .....................................................................    
  410   continue
c .....................................................................    
  405 continue
c .....................................................................
c
c .....................................................................
c
  430 continue
c ... forca e fluxo distribuida no contorno
c     iq(1) = 1 | no 1 2 3 4  9 10 11 12 |
c             2 | no 5 6 7 8 13 14 15 16 |
c             3 | no 1 5 6 2 17 13 18  9 |     
c             4 | no 4 3 7 8 11 19 15 20 |
c             5 | no 1 4 8 5 12 20 16 17 |
c             6 | no 2 6 7 3 18 14 19 10 |
c
c ... verifica se ha alguma face com carga
      tp = 0
      do 435 i = 1, 6
        tp    = tp + iq(i)
  435 continue
c .....................................................................
c
c ...
      if( tp .gt. 0 ) then
        do 440 j = 1, 6
c ... face 
          if(iq(j) .gt. 0 ) then
c ...
            do 445 i = 1, 4
              no       = hexa_face_node20(i,j)
              xf(1,i) = x(1,no) 
              xf(2,i) = x(2,no) 
              xf(3,i) = x(3,no)
  445       continue
c ...    
            call rotmatrix(xf,r)
            call rotate(xf(1,1),r,xf(1,1),.false.)
            call rotate(xf(1,2),r,xf(1,2),.false.)
            call rotate(xf(1,3),r,xf(1,3),.false.)
            call rotate(xf(1,4),r,xf(1,4),.false.)
c ... forcas
            call tload(iq(j),0.d0,0.d0,ddum,face_f)
            nint = 3
            do 450 ly = 1, nint
              si = pg(ly,nint)
              do 455 lx = 1, nint
                ri = pg(lx,nint)
c ...                                               
                call sfquad4_m(h,hx,hy,ri,si,.false.,.true.)
                call jacob2d_m(hx,hy,xj2D,xji2D,xf,det,4,ndm
     .                        ,nel,.true.)
c .....................................................................
c
c ...                                               
                w   = wg(lx,nint)*wg(ly,nint)*det
c .....................................................................
c
c ...                                               
                call sfquad8_m(hu,hux,huy,ri,si,.true.,.false.)
c .....................................................................
c
c ...
                do 460 i = 1, 8
                  no  = hexa_face_node20(i,j)
c                   l1  = (no-1)*3+1
                  wt1 = hu(i)*w
                  l1    = 3*no-2
                  l2    = l1 + 1
                  l3    = l2 + 1
                  p(l1) = p(l1) - face_f(1)*wt1
                  p(l2) = p(l2) - face_f(2)*wt1
                  p(l3) = p(l3) - face_f(3)*wt1  
  460           continue
c .....................................................................
  455         continue
c .....................................................................             
  450       continue
c .....................................................................
          endif
c .....................................................................
  440   continue
c .....................................................................
      endif
c .....................................................................
      return
c ======================================================================
c
c ... Tensoes iniciais:                  
c
c ......................................................................
  500 continue
      nint = 4 
      do 510 lz = 1, nint
        ti = pg(lz,nint)
        do 520 ly = 1, nint
          si = pg(ly,nint)
          do 530 lx = 1, nint
            ri = pg(lx,nint)
c ...                                               
            call sfhexa8_m(h,hx,hy,hz,ri,si,ti,.true.,.true.)
            call jacob3d_m(hx,hy,hz,xj,xji,x,det,8,nel,.true.)
c .....................................................................
c
c ...
            w   = wg(lx,nint)*wg(ly,nint)*wg(lz,nint)*det
            wt1 = w  
c .....................................................................
c
c ...
            call sfhexa20_m(hu,hux,huy,huz,ri,si,ti,.false.,.true.)
            call jacob3d_m(hux,huy,huz,xj,xji,x,det,20,nel,.false.)
c .....................................................................
c
c ...
            txi(1:6)= 0.d0
            do 540 j = 1,   8
              txi(1) = txi(1) + h(j)*txn(1,j)
              txi(2) = txi(2) + h(j)*txn(2,j) 
              txi(3) = txi(3) + h(j)*txn(3,j) 
              txi(4) = txi(4) + h(j)*txn(4,j) 
              txi(5) = txi(5) + h(j)*txn(5,j)
              txi(6) = txi(6) + h(j)*txn(6,j) 
  540       continue 
c ..................................................................... 
c
c ...
            do 560 i = 1, 20
c             l1  = (i-1)*3+1
              l1  = 3*i-2
              l2  = l1 + 1
              l3  = l2 + 1
c ... Fu = int(BeT*sigma*dv)
              p(l1) = p(l1) 
     .        + wt1*(hux(i)*txi(1) + huy(i)*txi(4) + huz(i)*txi(6))
c
              p(l2) = p(l2) 
     .        + wt1*(hux(i)*txi(4) + huy(i)*txi(2) + huz(i)*txi(5))
c
               p(l3) = p(l3) 
     .        + wt1*(hux(i)*txi(6) + huy(i)*txi(5) + huz(i)*txi(3))
c ...................................................................... 
  560       continue  
c .....................................................................
  530     continue
c .....................................................................
  520   continue
c .....................................................................
  510 continue
c .....................................................................
c
c ....
      return
c ======================================================================
c
c ... Matriz geometrica:
c
c ......................................................................
  600 continue
      stop
c ....................................................................
      end
c *********************************************************************
 
