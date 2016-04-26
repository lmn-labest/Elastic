       subroutine elmlib_mec(e ,iq  ,x   ,u      ,p  ,s,txn
     .                     ,dt,ndm ,nst ,nel    ,iel,isw
     .                     ,ma,nlit,ilib)
c **********************************************************************
c *                                                                    *
c *   ELMLIB: biblioteca de elementos.                                 *
c *   ------                                                           *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ---------------------                                            *
c *                                                                    *
c *     e(10) - constantes fisicas                                     *
c *     iq(7) - cargas nos elementos                                   *
c *     x(ndm,nen)- coordenadas nodais locais                          *
c *     u(nst)     - solucao anterior                                  *
c *     p(nst)     - nao definido                                      *
c *     s(nst,nst)   - nao definido                                    *
c *     txn(6,nen) - tensoes nodais                                    *
c *     ndm - dimensao                                                 *
c *     nst - numero de graus de liberdade por elemento                *
c *     nel - numero do elemento                                       *
c *     iel - tipo   do elemento                                       *
c *     isw - codigo de instrucao                                      *
c *     ma   -  numero de material do elemento                         *
c *     nlit -  numero da iteracao nao linear                          *
c *     ilib -  codigo da biblioteca                                   *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   -------------------                                              *
c *                                                                    *
c *     e - constantes fisicas                                         *
c *     s - matriz de elemento                                         *
c *     p - isw = 2  residuo                                           *
c *         isw = 3  tensao e fluxo                                    *
c *         isw = 4  cargas de superfice, volume e integras do passo   *
c *         de tempo anterior                                          *
c **********************************************************************
      implicit none
      integer iq(*),iel,nel,ndm,nst,isw,ilib,ma,nlit
      real*8 e(*),x(*),u(*),p(*),s(nst,*),txn(6,*)
      real*8 dt
c ......................................................................
      goto (100,200,300,400,500,600,700,800,900,1000,1100,1200,1300) iel
   10 write(*,2000) iel,nel
      stop
c ......................................................................
  100 continue
      goto 10
      return
c ......................................................................            
  200 continue
      goto 10
      return
c ......................................................................
  300 continue
      goto 10
      return
c ......................................................................
  400 continue
      if (ilib .eq. 1) then  
c     Elemento quadrilatero de 4 nos (mec-elastico)
        call elmt04_mec(e,iq,x,u,p,s,txn,ndm,nst,nel,isw)
      endif   
      return
c ......................................................................
  500 continue
      goto 10
      return
c ......................................................................
  600 continue
      if (ilib .eq. 1) then  
c     Elemento tetraedrico de 4 nos (mec-elastico)
        call elmt06_mec(e,iq,x,u,p,s,txn,ndm,nst,nel,isw)
      endif
      return
c ......................................................................
  700 continue
      if (ilib .eq. 1) then  
c     Elemento hexaedrico de 8 nos (mec-elastico)
        call elmt07_mec(e,iq,x,u,p,s,txn,ndm,nst,nel,isw)
      endif
      return
c ......................................................................
  800 continue
      goto 10
      return
c ......................................................................            
  900 continue
      goto 10
      return 
c ......................................................................
 1000 continue
      goto 10
      return 
c ......................................................................
 1100 continue
      goto 10
      return       
c ......................................................................
 1200 continue
      if (ilib .eq. 1) then  
c     Elemento tetraedro de 10 nos (mec-elastico)
        call elmt12_mec(e,iq,x,u,p,s,txn,ndm,nst,nel,isw)
      endif
      return       
c ......................................................................
c ......................................................................
 1300 continue
      if (ilib .eq. 1) then  
c     Elemento hexaedrico de 20 nos (mec-elastico)
        call elmt13_mec(e,iq,x,u,p,s,txn,ndm,nst,nel,isw)
      endif
      return       
c ......................................................................         
 2000 format(1x,'SUBROUTINE ELMLIBPMEC:'
     .,/,5x,'tipo de elemento ',i2,' nao existente, elemento ',i6,' !')
       end
c **********************************************************************
c
c **********************************************************************
      subroutine uhx(hx,hy,vx,vy,hxu,nen)
c **********************************************************************
c *                                                                    *
c *   Termos de estabilzacao CG (characteristic-Galerkin)              *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *                                                                    *
c *   hx(nen) - derivadas de h                                         *
c *   hy(nen) - derivadas de h                                         *
c *   vx,vy - componentes da velocidade                                *
c *   kd    - coeficiente de difusao                                   *
c *   hc    - tamanho caracteristico do elemento                       *
c *                                                                    *
c *   Parametros de saida:                                             *
c *                                                                    *
c *   hpg(nen) - termos de estabilizacao: hpg = v.grad(h)              *
c *                                                                    *
c **********************************************************************
      implicit none
      integer nen,i
      real*8  hx(*),hy(*),vx,vy,hxu(*)
c ......................................................................
      do 100 i = 1, nen
         hxu(i) = vx*hx(i) + vy*hy(i)
  100 continue
c ......................................................................      
      return
      end           
      subroutine hpgquad4(h,hx,hy,vx,vy,kd,hc,hpg)
c **********************************************************************
c *                                                                    *
c *                                                                    *
c *   HPGQUAD4: funcoes de interpolacao de Petrov-Galerkin             *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *                                                                    *
c *   h(4)  - funcoes de interpolacao                                  *
c *   hx(4) - derivadas de h                                           *
c *   hy(4) - derivadas de h                                           *
c *   vx,vy - componeentes da velocodade                               *
c *   kd    - coeficiente de difusao                                   *
c *   hc    - tamanho caracteristico do elemento                       *
c *                                                                    *
c *   Parametros de saida:                                             *
c *                                                                    *
c *   hpg(4) - funcoes de Petrov-Galerkin                              *
c *                                                                    *
c **********************************************************************
      implicit none
      integer i
      real*8  h(*),hx(*),hy(*),vx,vy,kd,hc,hpg(*)
      real*8 zero,peclet,alpha,xk,vm
      parameter (zero = 1.d-14)
c ......................................................................
c
c ... Modulo da velocidade:
c
      vm = dsqrt(vx*vx + vy*vy)
c
c ... Parametro de upwind:
c
      if(kd .gt. zero) then
           peclet = (hc*vm)/(2.d0*kd)
           alpha  = min((1.d0/3.d0)*peclet,1.d0)
      else
           alpha = 1.d0
      end if
c
c ... difusao artificial:
c
      xk = 0.d0
      if(vm .gt. zero) xk = alpha*hc*0.5d0/vm
c
c ... funcoes de Petrov-Galerkin:
c
      do 100 i = 1, 4
         hpg(i) = h(i) + xk*(vx*hx(i) + vy*hy(i))
  100 continue
c ......................................................................      
      return
      end
      subroutine hpghexa8(h,hx,hy,hz,vx,vy,vz,kd,hc,hpg)
c **********************************************************************
c *                                                                    *
c *                                                                    *
c *   HPGHEXA8: funcoes de interpolacao de Petrov-Galerkin             *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *                                                                    *
c *   h(8)  - funcoes de interpolacao                                  *
c *   hx(8) - derivadas de h                                           *
c *   hy(8) - derivadas de h                                           *
c *   hz(8) - derivadas de h                                           *
c *   vx,vy,vz - componeentes da velocodade                            *
c *   kd    - coeficiente de difusao                                   *
c *   hc    - tamanho caracteristico do elemento                       *
c *                                                                    *
c *   Parametros de saida:                                             *
c *                                                                    *
c *   hpg(4) - funcoes de Petrov-Galerkin                              *
c *                                                                    *
c **********************************************************************
      implicit none
      integer i
      real*8  h(*),hx(*),hy(*),hz(*),vx,vy,vz,kd,hc,hpg(*)
      real*8 zero,peclet,alpha,xk,vm
      parameter (zero = 1.d-14)
c ......................................................................
c
c ... Modulo da velocidade:
c
      vm = dsqrt(vx*vx + vy*vy + vz*vz)
c
c ... Parametro de upwind:
c
      if(kd .gt. zero) then
           peclet = (hc*vm)/(2.d0*kd)
           alpha  = min((1.d0/3.d0)*peclet,1.d0)
      else
           alpha = 1.d0
      end if
c
c ... difusao artificial:
c
      xk = 0.d0
      if(vm .gt. zero) xk = alpha*hc*0.5d0/vm
c
c ... funcoes de Petrov-Galerkin:
c
      do 100 i = 1, 8
         hpg(i) = h(i) + xk*(vx*hx(i) + vy*hy(i) + vz*hz(i))
  100 continue
c ......................................................................      
      return
      end      
      subroutine lku(s,u,p,nst)
c **********************************************************************
c *                                                                    *
c *   LKU: forcas internas K.U                                         *
c *   ---                                                              *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ---------------------                                            *
c *                                                                    *
c *     s    - matriz de elemento                                      *
c *     u    - deslocamentos                                           *
c *     nst  - numero de graus de liberdade por elemento               *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   -------------------                                              *
c *                                                                    *
c *     p = s.u - forcas internas                                      *
c *                                                                    *
c **********************************************************************
      implicit none
      integer nst,i,j
      real*8  s(nst,*),u(*),p(*)
c ......................................................................
      do 200 i = 1, nst
         do 100 j = 1, nst
            p(i) = p(i) + s(i,j)*u(j)
  100    continue
  200 continue
c ......................................................................  
      return
      end
      subroutine lku_m(s,u,p,nst)
c **********************************************************************
c * Data de criacao    : 12/12/2016                                    *
c * Data de modificaco :                                               * 
c * ------------------------------------------------------------------ *
c * LKU_M : forcas internas K.U                                        *  
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ * 
c * s    - matriz de elemento                                          *
c * u    - deslocamentos                                               *
c * nst  - numero de graus de liberdade por elemento                   *
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ * 
c * p = s.u - forcas internas                                          *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ * 
c * Efetuando o produto varrendo as linhas                             *     
c **********************************************************************
      implicit none
      integer nst,i,j
      real*8  s(nst,*),u(*),p(*)
c ......................................................................
      p(1:nst) = 0.d0
      do 200 j = 1, nst
         do 100 i = 1, nst
            p(i) = p(i) + s(i,j)*u(j)
  100    continue
  200 continue
c ......................................................................  
      return
      end
      subroutine lku_sym(s,u,p,nst)
c **********************************************************************
c * Data de criacao    : 02/04/2016                                    *
c * Data de modificaco :                                               * 
c * ------------------------------------------------------------------ *
c * LKU_SYM  : forcas internas K.U para matriz k simetria              *  
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ * 
c * s    - matriz de elemento                                          *
c * u    - deslocamentos                                               *
c * nst  - numero de graus de liberdade por elemento                   *
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ * 
c * p = s.u - forcas internas                                          *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ * 
c * Matriz K com armazenamento da parte triangular inferior            *     
c **********************************************************************
      implicit none
      integer nst,i,j
      real*8  s(nst,*),u(*),p(*),uj,sij
c ......................................................................
      p(1:nst) = 0.d0
      do 200 j = 1, nst
         uj   = u(j) 
         p(j) = p(j) + s(j,j)*uj  
         do 100 i = j+1, nst 
            sij  = s(i,j) 
            p(i) = p(i) + sij*uj 
            p(j) = p(j) + sij*u(i)
  100    continue
  200 continue
c ......................................................................  
      return
      end
      subroutine addlku(s,u,p,nst)
c **********************************************************************
c *                                                                    *
c *   ADDLKU: forcas internas K.U                                      *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *                                                                    *
c *     s    - matriz de elemento                                      *
c *     u    - deslocamentos                                           *
c *     nst  - numero de graus de liberdade por elemento               *
c *                                                                    *
c *   Parametros de saida:                                             *
c *                                                                    *
c *     p = s.u - forcas internas                                      *
c *                                                                    *
c **********************************************************************
      implicit none
      integer nst,i,j
      real*8  s(nst,*),u(*),p(*)
c ......................................................................
      do 200 i = 1, nst
         do 100 j = 1, nst
            p(i) = p(i) + s(i,j)*u(j)
  100    continue
  200 continue
c ......................................................................  
      return
      end      
      subroutine addlkum(s,u,p,nst)
c **********************************************************************
c *                                                                    *
c *   ADDLKU: forcas internas K.U                                      *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *                                                                    *
c *     s    - matriz de elemento                                      *
c *     u    - deslocamentos                                           *
c *     nst  - numero de graus de liberdade por elemento               *
c *                                                                    *
c *   Parametros de saida:                                             *
c *                                                                    *
c *     p = s.u - forcas internas                                      *
c *                                                                    *
c **********************************************************************
      implicit none
      integer nst,i,j
      real*8  s(nst,*),u(*),p(*)
c ......................................................................
      do 200 j = 1, nst
         do 100 i = 1, nst
            p(i) = p(i) + s(i,j)*u(j)
  100    continue
  200 continue
c ......................................................................  
      return
      end      
      subroutine sftria3(h,hr,hs,r,s,afl,bfl)
c **********************************************************************
c *                                                                    *
c *                                                    22/12/01        *
c *   SFTRIA3:                                                         *
c *   -------                                                          *
c *                                                                    *
c *   Calcula as funcoes de interpolacao e suas derivadas              *
c *   no ponto (r,s) do triangulo de 3 nos.                            *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ---------------------                                            *
c *                                                                    *
c *     r,s - coordenadas naturais                                     *
c *     afl = true ( calcula funcoes )                                 *
c *     bfl = true ( calcula as derivadas )                            *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   -------------------                                              *
c *                                                                    *
c *     h(6)    - funcoes de interolocao no ponto (r,s)                *
c *     hr(6)   - derivadas de h em relacao a r                        *
c *     hs(6)   - derivadas de h em relacao a s                        *
c *                                                                    *
c **********************************************************************
      implicit none
      real*8  h(*),hr(*),hs(*),r,s,t
      logical afl,bfl
c ......................................................................
      t = 1.d0 - r - s
      if (afl) then
c
c ...... Funcoes de interpolacao lineares standard:
c
         h(1) = r
         h(2) = s
         h(3) = t
      endif
      if (bfl) then
c
c ...... Derivadas em relacao a r :
c
         hr(1) =  1.d0
         hr(2) =  0.d0
         hr(3) = -1.d0
c
c ...... Derivadas em relacao a s :
c
         hs(1) =  0.d0
         hs(2) =  1.d0
         hs(3) = -1.d0
      endif
c ......................................................................                    
      return
      end
      subroutine sftria6(h,hr,hs,r,s,afl,bfl)
c **********************************************************************
c *                                                                    *
c *                                                    22/12/01        *
c *   SFTRIA6:                                                         *
c *   -------                                                          *
c *                                                                    *
c *   Calcula as funcoes de interpolacao e suas derivadas              *
c *   no ponto (r,s) do triangulo de 6 nos.                            *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ---------------------                                            *
c *                                                                    *
c *     r,s - coordenadas naturais                                     *
c *     afl = true ( calcula funcoes )                                 *
c *     bfl = true ( calcula as derivadas )                            *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   -------------------                                              *
c *                                                                    *
c *     h(6)    - funcoes de interolocao no ponto (r,s)                *
c *     hr(6)   - derivadas de h em relacao a r                        *
c *     hs(6)   - derivadas de h em relacao a s                        *
c *                                                                    *
c **********************************************************************
      implicit none
      real*8  h(*),hr(*),hs(*),r,s,t
      logical afl,bfl
c ......................................................................
      t = 1.d0 - r - s
      if (afl) then
c
c ...... Funcoes de interpolacao quadraticas standard:
c
         h(1) = r * (2.d0*r-1.d0)
         h(2) = s * (2.d0*s-1.d0)
         h(3) = t * (2.d0*t-1.d0)
         h(4) = 4.d0 * r * s
         h(5) = 4.d0 * s * t
         h(6) = 4.d0 * t * r
      endif
      if (bfl) then
c
c ...... Derivadas em relacao a r :
c
         hr(1) =  4.d0* r - 1.d0
         hr(2) =  0.d0
         hr(3) =  1.d0- 4.d0* t
         hr(4) =  4.d0* s
         hr(5) = -4.d0* s
         hr(6) =  4.d0*(t-r)
c
c ...... Derivadas em relacao a s :
c
         hs(1) =  0.d0
         hs(2) =  4.d0* s - 1.d0
         hs(3) =  1.d0 - 4.d0* t
         hs(4) =  4.d0* r
         hs(5) =  4.d0*(t-s)
         hs(6) = -4.d0* r
      endif
c ......................................................................      
      return
      end
      subroutine sftetra4(h,hr,hs,ht,r,s,t,afl,bfl)
c **********************************************************************
c *                                                                    *
c *                                                    16/03/2016      *
c *   SFTETRA4:                                                        *
c *   ---------                                                        *
c *                                                                    *
c *   Calcula as funcoes de interpolacao e suas derivadas              *
c *   no ponto (r,s,t) do tetraedro de  4 nos.                         *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ---------------------                                            *
c *                                                                    *
c *   r,s,t - coordenadas naturais                                     *
c *   afl = true ( calcula funcoes )                                   *
c *   bfl = true ( calcula as derivadas )                              *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   -------------------                                              *
c *                                                                    *
c *   h( 4)    - funcoes de interolocao no ponto (r,s)                 *
c *   hr( 4)   - derivadas de h em relacao a r                         *
c *   hs( 4)   - derivadas de h em relacao a s                         *
c *   ht( 4)   - derivadas de h em relacao a t                         *
c *                                                                    *
c * no  (   r,   s,   t,  u)                                           *
c * no1 (   1,   0,   0,  0)                                           * 
c * no2 (   0,   1,   0,  0)                                           * 
c * no3 (   0,   0,   1,  0)                                           * 
c * no4 (   0,   0,   0,  1)                                           * 
c **********************************************************************
      implicit none
      real*8  h(*),hr(*),hs(*),ht(*),r,s,t,u
      logical afl,bfl
c ......................................................................
      u = 1.d0 - r - s - t
      if (afl) then
c
c ...... Funcoes de interpolacao lineares standard:
c
         h(1) = r
         h(2) = s
         h(3) = t
         h(4) = u         
      endif
      if (bfl) then
c
c ...... Derivadas em relacao a r :
c
         hr(1) =   1.d0
         hr(2) =   0.d0
         hr(3) =   0.d0
         hr(4) =  -1.d0
c
c ...... Derivadas em relacao a s :
c
         hs(1) =  0.d0
         hs(2) =  1.d0           
         hs(3) =  0.d0
         hs(4) = -1.d0               
c
c ...... Derivadas em relacao a t :
c
         ht(1) =  0.d0
         ht(2) =  0.d0
         ht(3) =  1.d0               
         ht(4) = -1.d0             
      endif
c ......................................................................      
      return
      end      
      subroutine sftetra10(h,hr,hs,ht,r,s,t,afl,bfl)
c **********************************************************************
c *                                                                    *
c *                                                    22/12/01        *
c *   SFTETRA10:                                                       *
c *   ---------                                                        *
c *                                                                    *
c *   Calcula as funcoes de interpolacao e suas derivadas              *
c *   no ponto (r,s,t) do tetraedro de 10 nos.                         *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ---------------------                                            *
c *                                                                    *
c *   r,s,t - coordenadas naturais                                     *
c *   afl = true ( calcula funcoes )                                   *
c *   bfl = true ( calcula as derivadas )                              *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   -------------------                                              *
c *                                                                    *
c *   h(10)    - funcoes de interolocao no ponto (r,s)                 *
c *   hr(10)   - derivadas de h em relacao a r                         *
c *   hs(10)   - derivadas de h em relacao a s                         *
c *   ht(10)   - derivadas de h em relacao a t                         *
c *                                                                    *
c * no  (   r,   s,   t,  u)                                           *
c * no1 (   1,   0,   0,  0)                                           * 
c * no2 (   0,   1,   0,  0)                                           * 
c * no3 (   0,   0,   1,  0)                                           * 
c * no4 (   0,   0,   0,  1)                                           * 
c * no5 ( 1/2, 1/2,   0,  0)                                           * 
c * no6 ( 1/2,   0, 1/2,  0)                                           * 
c * no7 ( 1/2,   0,   0,1/2)                                           * 
c * no8 (   0, 1/2,   0,  0)                                           * 
c * no9 (   0,   0, 1/2,1/2)                                           * 
c * no10(   0, 1/2,   0,1/2)                                           * 
c **********************************************************************
      implicit none
      real*8  h(*),hr(*),hs(*),ht(*),r,s,t,u
      logical afl,bfl
c ......................................................................
      u = 1.d0 - r - s - t
      if (afl) then
c
c ...... Funcoes de interpolacao quadraticas standard:
c
         h(1) = r * (2.d0*r-1.d0)
         h(2) = s * (2.d0*s-1.d0)
         h(3) = t * (2.d0*t-1.d0)
         h(4) = u * (2.d0*u-1.d0)         
         h(5) = 4.d0 * r * s
         h(6) = 4.d0 * r * t
         h(7) = 4.d0 * r * u
         h(8) = 4.d0 * s * t
         h(9) = 4.d0 * t * u
         h(10)= 4.d0 * s * u         
      endif
      if (bfl) then
c
c ...... Derivadas em relacao a r :
c
         hr(1) =  4.d0*r - 1.d0
         hr(2) =  0.d0
         hr(3) =  0.d0
         hr(4) =  1.d0 - 4.d0*u
         hr(5) =  4.d0*s
         hr(6) =  4.d0*t
         hr(7) =  4.d0*(u-r)
         hr(8) =  0.d0
         hr(9) = -4.d0*t 
         hr(10)= -4.d0*s         
c
c ...... Derivadas em relacao a s :
c
         hs(1) =  0.d0
         hs(2) =  4.d0*s - 1.d0
         hs(3) =  0.d0
         hs(4) =  1.d0 - 4.d0*u
         hs(5) =  4.d0*r
         hs(6) =  0.d0
         hs(7) = -4.d0*r
         hs(8) =  4.d0*t
         hs(9) = -4.d0*t 
         hs(10)=  4.d0*(u-s)
c
c ...... Derivadas em relacao a t :
c
         ht(1) =  0.d0
         ht(2) =  0.d0
         ht(3) =  4.d0*t - 1.d0
         ht(4) =  1.d0 - 4.d0*u
         ht(5) =  0.d0
         ht(6) =  4.d0*r
         ht(7) = -4.d0*r
         ht(8) =  4.d0*s
         ht(9) =  4.d0*(u-t)
         ht(10)= -4.d0*s
      endif
c ......................................................................      
      return
      end      
      subroutine sfquad4(h,hx,hy,r,s,afl,bfl)
c **********************************************************************
c *                                                                    *
c *                                                    27/01/89        *
c *   SFQUAD4:                                                         *
c *                                                                    *
c *   Calcula as funcoes de interpolacao e suas derivadas              *
c *   no ponto (r,s) do quadrilatero de 4 nos.                         *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *                                                                    *
c *   r,s - coordenadas naturais                                       *
c *   afl = true ( calcula funcoes )                                   *
c *   bfl = true ( calcula as derivadas )                              *
c *                                                                    *
c *   Parametros de saida:                                             *
c *                                                                    *
c *   h(4)    - funcoes de interolocao no ponto (r,s)                  *
c *   hx(4)   - derivadas de h em relacao a r                          *
c *   hy(4)   - derivadas de h em relacao a s                          *
c * no  ( r, s)                                                        *
c * no1 ( 1, 1)                                                        * 
c * no2 (-1, 1)                                                        * 
c * no3 (-1,-1)                                                        * 
c * no4 ( 1,-1)                                                        * 
c **********************************************************************
      implicit none
      logical afl,bfl
      real*8  h(*),hx(*),hy(*),r,s
c ......................................................................
      if ( afl ) then
c
c     funcoes de interpolacao
c
      h(1)  = ( 1.d0+r ) * ( 1.d0+s ) / 4.d0
      h(2)  = ( 1.d0-r ) * ( 1.d0+s ) / 4.d0
      h(3)  = ( 1.d0-r ) * ( 1.d0-s ) / 4.d0
      h(4)  = ( 1.d0+r ) * ( 1.d0-s ) / 4.d0
c
      endif
      if ( bfl ) then
c
c     derivadas em relacao a r :
c
      hx(1)  =   (1.d0+s) / 4.d0
      hx(2)  = - (1.d0+s) / 4.d0
      hx(3)  = - (1.d0-s) / 4.d0
      hx(4)  =   (1.d0-s) / 4.d0
c
c     derivadas em relacao a s :
c
      hy(1)  =  (1.d0+r) / 4.d0
      hy(2)  =  (1.d0-r) / 4.d0
      hy(3)  = -(1.d0-r) / 4.d0
      hy(4)  = -(1.d0+r) / 4.d0
c
      endif
c ......................................................................      
      return
      end
      subroutine sfquad4_m(h,hx,hy,r,s,afl,bfl)
c **********************************************************************
c *                                                                    *
c *                                                    05/11/15        *
c *   SFQUAD4:                                                         *
c *                                                                    *
c *   Calcula as funcoes de interpolacao e suas derivadas              *
c *   no ponto (r,s) do quadrilatero de 4 nos.                         *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *                                                                    *
c *   r,s - coordenadas naturais                                       *
c *   afl = true ( calcula funcoes )                                   *
c *   bfl = true ( calcula as derivadas )                              *
c *                                                                    *
c *   Parametros de saida:                                             *
c *                                                                    *
c *   h(4)    - funcoes de interolocao no ponto (r,s)                  *
c *   hx(4)   - derivadas de h em relacao a r                          *
c *   hy(4)   - derivadas de h em relacao a s                          *
c * OBS:                                                               *
c * no  ( r, s)                                                        *
c * no1 ( 1, 1)                                                        * 
c * no2 (-1, 1)                                                        * 
c * no3 (-1,-1)                                                        * 
c * no4 ( 1,-1)                                                        * 
c **********************************************************************
      implicit none
      logical afl,bfl
      real*8 one_plus_r,one_plus_s,one_minus_r,one_minus_s
      real*8  h(*),hx(*),hy(*),r,s
c ......................................................................
c
c ...
      one_plus_r  = 1.d0 + r
      one_plus_s  = 1.d0 + s
      one_minus_r = 1.d0 - r
      one_minus_s = 1.d0 - s
c ......................................................................
c
c ...   
      if ( afl ) then
c
c     funcoes de interpolacao
c
        h(1)  = one_plus_r *one_plus_s *0.25d0
        h(2)  = one_minus_r*one_plus_s *0.25d0
        h(3)  = one_minus_r*one_minus_s*0.25d0
        h(4)  = one_plus_r *one_minus_s*0.25d0
c
      endif
      if ( bfl ) then
c
c     derivadas em relacao a r :
c
        hx(1)  =   one_plus_s *0.25d0
        hx(2)  = - one_plus_s *0.25d0
        hx(3)  = - one_minus_s*0.25d0
        hx(4)  =   one_minus_s*0.25d0
c
c     derivadas em relacao a s :
c
        hy(1)  =  one_plus_r *0.25d0
        hy(2)  =  one_minus_r*0.25d0
        hy(3)  = -one_minus_r*0.25d0
        hy(4)  = -one_plus_r *0.25d0
c
      endif
c ......................................................................      
      return
      end
c **********************************************************************
      subroutine sfquad8(h,hx,hy,r,s,afl,bfl)
c **********************************************************************
c *                                                                    *
c *                                                    27/01/89        *
c *   SFQUAD8:                                                         *
c *   -------                                                          *
c *                                                                    *
c *   Calcula as funcoes de interpolacao e suas derivadas              *
c *   no ponto (r,s) do quadrilatero quadratico de 8 nos.              *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ---------------------                                            *
c *                                                                    *
c *     r,s - coordenadas naturais                                     *
c *     afl = true ( calcula funcoes )                                 *
c *     bfl = true ( calcula as derivadas )                            *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   -------------------                                              *
c *                                                                    *
c *     h(8)    - funcoes de interolocao no ponto (r,s)                *
c *     hx(8)   - derivadas de h em relacao a r                        *
c *     hy(8)   - derivadas de h em relacao a s                        *
c *                                                                    *
c **********************************************************************
      implicit none
      logical afl,bfl
      real*8  h(*),hx(*),hy(*),r,s
c ......................................................................
      if ( afl ) then
c
c     funcoes de interpolacao
c
      h(5)  = ( 1.d0-r*r ) * ( 1.d0+s ) / 2.
      h(6)  = ( 1.d0-s*s ) * ( 1.d0-r ) / 2.
      h(7)  = ( 1.d0-r*r ) * ( 1.d0-s ) / 2.
      h(8)  = ( 1.d0-s*s ) * ( 1.d0+r ) / 2.
      h(1)  = ( 1.d0+r ) * ( 1.d0+s ) / 4.d0 - h(5) / 2.d0 - h(8) / 2.d0
      h(2)  = ( 1.d0-r ) * ( 1.d0+s ) / 4.d0 - h(5) / 2.d0 - h(6) / 2.d0
      h(3)  = ( 1.d0-r ) * ( 1.d0-s ) / 4.d0 - h(6) / 2.d0 - h(7) / 2.d0
      h(4)  = ( 1.d0+r ) * ( 1.d0-s ) / 4.d0 - h(7) / 2.d0 - h(8) / 2.d0
c
      endif
      if ( bfl ) then
c
c     derivadas em relacao a r :
c
      hx(5)  = - r * (1.d0+s)
      hx(6)  = -.5d0 * (1.d0-s*s)
      hx(7)  = - r * (1.d0-s)
      hx(8)  =  .5d0 * (1.d0-s*s)
      hx(1)  =   (1.d0+s) / 4.d0 - hx(5) / 2.d0 - hx(8) / 2.d0
      hx(2)  = - (1.d0+s) / 4.d0 - hx(5) / 2.d0 - hx(6) / 2.d0
      hx(3)  = - (1.d0-s) / 4.d0 - hx(6) / 2.d0 - hx(7) / 2.d0
      hx(4)  =   (1.d0-s) / 4.d0 - hx(7) / 2.d0 - hx(8) / 2.d0
c
c     derivadas em relacao a s :
c
      hy(5)  =  (1.d0-r*r) / 2.d0
      hy(6)  = - s * (1.d0-r)
      hy(7)  = -(1.d0-r*r) / 2.d0
      hy(8)  = - s * (1.+r)
      hy(1)  =  (1.d0+r) / 4.d0 - hy(5) / 2.d0 - hy(8) / 2.d0
      hy(2)  =  (1.d0-r) / 4.d0 - hy(5) / 2.d0 - hy(6) / 2.d0
      hy(3)  = -(1.d0-r) / 4.d0 - hy(6) / 2.d0 - hy(7) / 2.d0
      hy(4)  = -(1.d0+r) / 4.d0 - hy(7) / 2.d0 - hy(8) / 2.d0
c
      endif
c ......................................................................      
      return
      end
c **********************************************************************
      subroutine sfquad8_m(h,hx,hy,r,s,afl,bfl)
c **********************************************************************
c *                                                                    *
c *                                                    30/11/15        *
c *   SFQUAD8:                                                         *
c *   -------                                                          *
c *                                                                    *
c *   Calcula as funcoes de interpolacao e suas derivadas              *
c *   no ponto (r,s) do quadrilatero quadratico de 8 nos.              *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ---------------------                                            *
c *                                                                    *
c *     r,s - coordenadas naturais                                     *
c *     afl = true ( calcula funcoes )                                 *
c *     bfl = true ( calcula as derivadas )                            *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   -------------------                                              *
c *                                                                    *
c *     h(8)    - funcoes de interolocao no ponto (r,s)                *
c *     hx(8)   - derivadas de h em relacao a r                        *
c *     hy(8)   - derivadas de h em relacao a s                        *
c * OBS:                                                               *
c * no  ( r, s)                                                        *
c * no1 ( 1, 1)                                                        * 
c * no2 (-1, 1)                                                        * 
c * no3 (-1,-1)                                                        * 
c * no4 ( 1,-1)                                                        * 
c * no5 ( 0, 1)                                                        * 
c * no6 (-1, 0)                                                        * 
c * no7 ( 0,-1)                                                        * 
c * no8 ( 1, 0)                                                        * 
c **********************************************************************
      implicit none
      logical afl,bfl
      real*8 one_plus_r,one_plus_s,one_minus_r,one_minus_s
      real*8 one_minus_rr,one_minus_ss
      real*8  h(*),hx(*),hy(*),r,s
      real*8 ha5,ha6,ha7,ha8
c ...
      one_plus_r  = 1.d0 + r
      one_plus_s  = 1.d0 + s
      one_minus_r = 1.d0 - r
      one_minus_s = 1.d0 - s
      one_minus_rr= 1.d0 - r*r 
      one_minus_ss= 1.d0 - s*s
c ......................................................................
c
c ...                                                                   
      if ( afl ) then
c
c     funcoes de interpolacao
c
        h(5)  = one_minus_rr * one_plus_s * 0.5d0
        h(6)  = one_minus_ss * one_minus_r* 0.5d0
        h(7)  = one_minus_rr * one_minus_s* 0.5d0
        h(8)  = one_minus_ss * one_plus_r * 0.5d0
c
        ha5   =  h(5) * 0.5d0
        ha6   =  h(6) * 0.5d0
        ha7   =  h(7) * 0.5d0
        ha8   =  h(8) * 0.5d0
c       
        h(1)  = one_plus_r *one_plus_s * 0.25d0 - ha5 - ha8
        h(2)  = one_minus_r*one_plus_s * 0.25d0 - ha5 - ha6
        h(3)  = one_minus_r*one_minus_s* 0.25d0 - ha6 - ha7
        h(4)  = one_plus_r *one_minus_s* 0.25d0 - ha7 - ha8
c
      endif
      if ( bfl ) then
c
c     derivadas em relacao a r :
c
        hx(5)  = - r * one_plus_s
        hx(6)  = -0.5d0 *  one_minus_ss
        hx(7)  = - r * one_minus_s
        hx(8)  =  0.5d0 *  one_minus_ss
c
        ha5   =  hx(5)*0.5d0
        ha6   =  hx(6)*0.5d0
        ha7   =  hx(7)*0.5d0
        ha8   =  hx(8)*0.5d0
c      
        hx(1)  =   one_plus_s *0.25d0 - ha5 - ha8
        hx(2)  = - one_plus_s *0.25d0 - ha5 - ha6
        hx(3)  = - one_minus_s*0.25d0 - ha6 - ha7
        hx(4)  =   one_minus_s*0.25d0 - ha7 - ha8
c
c     derivadas em relacao a s :
c
        hy(5)  =  one_minus_rr * 0.5d0
        hy(6)  = - s * one_minus_r
        hy(7)  = -one_minus_rr * 0.5d0
        hy(8)  = - s * one_plus_r
c
        ha5   =  hy(5) * 0.5d0
        ha6   =  hy(6) * 0.5d0
        ha7   =  hy(7) * 0.5d0
        ha8   =  hy(8) * 0.5d0
c      
        hy(1)  =  one_plus_r * 0.25d0 - ha5 - ha8
        hy(2)  =  one_minus_r* 0.25d0 - ha5 - ha6 
        hy(3)  = -one_minus_r* 0.25d0 - ha6 - ha7
        hy(4)  = -one_plus_r * 0.25d0 - ha7 - ha8 
c
      endif
c ......................................................................      
      return
      end
      subroutine sfhexa8(h,hx,hy,hz,r,s,t,afl,bfl)
c **********************************************************************
c *                                                                    *
c *                                                    27/01/89        *
c *   SFHEXA8:                                                         *
c *   -------                                                          *
c *                                                                    *
c *   Calcula as funcoes de interpolacao e suas derivadas              *
c *   no ponto (r,s,t), do hexaedro linear de 8 nos.                   *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ---------------------                                            *
c *                                                                    *
c *     r,s,t - coordenadas naturais                                   *
c *     afl = true ( calcula funcoes )                                 *
c *     bfl = true ( calcula as derivadas )                            *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   -------------------                                              *
c *                                                                    *
c *     h(8)    - funcoes de interolocao no ponto (r,s,t)              *
c *     hx(8)   - derivadas de h em relacao a r                        *
c *     hy(8)   - derivadas de h em relacao a s                        *
c *     hz(8)   - derivadas de h em relacao a t                        *
c *                                                                    *
c **********************************************************************
      implicit none
      logical afl,bfl
      real*8  h(*),hx(*),hy(*),hz(*),r,s,t
c ......................................................................
      if ( afl ) then
c
c     funcoes de interpolacao
c
      h(1)  = ( 1.d0+r ) * ( 1.d0+s ) * (1.d0+t) / 8.d0
      h(2)  = ( 1.d0-r ) * ( 1.d0+s ) * (1.d0+t) / 8.d0
      h(3)  = ( 1.d0-r ) * ( 1.d0-s ) * (1.d0+t) / 8.d0
      h(4)  = ( 1.d0+r ) * ( 1.d0-s ) * (1.d0+t) / 8.d0
      h(5)  = ( 1.d0+r ) * ( 1.d0+s ) * (1.d0-t) / 8.d0
      h(6)  = ( 1.d0-r ) * ( 1.d0+s ) * (1.d0-t) / 8.d0
      h(7)  = ( 1.d0-r ) * ( 1.d0-s ) * (1.d0-t) / 8.d0
      h(8)  = ( 1.d0+r ) * ( 1.d0-s ) * (1.d0-t) / 8.d0
c
      endif
      if ( bfl ) then
c
c     derivadas em relacao a r :
c
      hx(1)  =   (1.d0+s) * (1.d0+t) / 8.d0
      hx(2)  = - (1.d0+s) * (1.d0+t) / 8.d0
      hx(3)  = - (1.d0-s) * (1.d0+t) / 8.d0
      hx(4)  =   (1.d0-s) * (1.d0+t) / 8.d0
      hx(5)  =   (1.d0+s) * (1.d0-t) / 8.d0
      hx(6)  = - (1.d0+s) * (1.d0-t) / 8.d0
      hx(7)  = - (1.d0-s) * (1.d0-t) / 8.d0
      hx(8)  =   (1.d0-s) * (1.d0-t) / 8.d0
c
c     derivadas em relacao a s :
c
      hy(1)  =  (1.d0+r) * (1.d0+t) / 8.d0
      hy(2)  =  (1.d0-r) * (1.d0+t) / 8.d0
      hy(3)  = -(1.d0-r) * (1.d0+t) / 8.d0
      hy(4)  = -(1.d0+r) * (1.d0+t) / 8.d0
      hy(5)  =  (1.d0+r) * (1.d0-t) / 8.d0
      hy(6)  =  (1.d0-r) * (1.d0-t) / 8.d0
      hy(7)  = -(1.d0-r) * (1.d0-t) / 8.d0
      hy(8)  = -(1.d0+r) * (1.d0-t) / 8.d0
c
c     derivadas em relacao a t :
c
      hz(1)  =  ( 1.d0+r ) * ( 1.d0+s ) / 8.d0
      hz(2)  =  ( 1.d0-r ) * ( 1.d0+s ) / 8.d0
      hz(3)  =  ( 1.d0-r ) * ( 1.d0-s ) / 8.d0
      hz(4)  =  ( 1.d0+r ) * ( 1.d0-s ) / 8.d0
      hz(5)  = -( 1.d0+r ) * ( 1.d0+s ) / 8.d0
      hz(6)  = -( 1.d0-r ) * ( 1.d0+s ) / 8.d0
      hz(7)  = -( 1.d0-r ) * ( 1.d0-s ) / 8.d0
      hz(8)  = -( 1.d0+r ) * ( 1.d0-s ) / 8.d0
      endif
c ......................................................................      
      return
      end
      subroutine sfhexa8_m(h,hx,hy,hz,r,s,t,afl,bfl)
c **********************************************************************
c *                                                                    *
c *                                                    06/11/15        *
c *   SFHEXA8:                                                         *
c *   -------                                                          *
c *                                                                    *
c *   Calcula as funcoes de interpolacao e suas derivadas              *
c *   no ponto (r,s,t), do hexaedro linear de 8 nos.                   *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ---------------------                                            *
c *                                                                    *
c *     r,s,t - coordenadas naturais                                   *
c *     afl = true ( calcula funcoes )                                 *
c *     bfl = true ( calcula as derivadas )                            *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   -------------------                                              *
c *                                                                    *
c *     h(8)    - funcoes de interolocao no ponto (r,s,t)              *
c *     hx(8)   - derivadas de h em relacao a r                        *
c *     hy(8)   - derivadas de h em relacao a s                        *
c *     hz(8)   - derivadas de h em relacao a t                        *
c * OBS:                                                               *
c * no  ( r, s, t)                                                     *
c * no1 ( 1, 1, 1)                                                     * 
c * no2 (-1, 1, 1)                                                     * 
c * no3 (-1,-1, 1)                                                     * 
c * no4 ( 1,-1, 1)                                                     * 
c * no5 ( 1, 1,-1)                                                     * 
c * no6 (-1, 1,-1)                                                     * 
c * no7 (-1,-1,-1)                                                     * 
c * no8 ( 1,-1,-1)                                                     * 
c **********************************************************************
      implicit none
      logical afl,bfl
      real*8 onePlusR,onePlusS,onePlusT
      real*8 oneMinusR,oneMinusS,oneMinusT
      real*8  h(*),hx(*),hy(*),hz(*),r,s,t
c ......................................................................
c
c ...
      onePlusR  = 1.d0 + r
      onePlusS  = 1.d0 + s
      onePlusT  = 1.d0 + t
      oneMinusR = 1.d0 - r
      oneMinusS = 1.d0 - s
      oneMinusT = 1.d0 - t
c ......................................................................
c
c ... funcoes de interpolacao
c
      if ( afl ) then
        h(1)  = onePlusR *onePlusS *onePlusT*0.125d0
        h(2)  = oneMinusR*onePlusS *onePlusT*0.125d0
        h(3)  = oneMinusR*oneMinusS*onePlusT*0.125d0
        h(4)  = onePlusR *oneMinusS*onePlusT*0.125d0
c
        h(5)  = onePlusR *onePlusS *oneMinusT*0.125d0
        h(6)  = oneMinusR*onePlusS *oneMinusT*0.125d0
        h(7)  = oneMinusR*oneMinusS*oneMinusT*0.125d0
        h(8)  = onePlusR *oneMinusS*oneMinusT*0.125d0
      endif
c .....................................................................
c
c ... derivadas
      if ( bfl ) then
c
c     derivadas em relacao a r :
c
        hx(1)  =  onePlusS *onePlusT*0.125d0
        hx(2)  = -onePlusS *onePlusT*0.125d0
        hx(3)  = -oneMinusS*onePlusT*0.125d0
        hx(4)  =  oneMinusS*onePlusT*0.125d0
c
        hx(5)  =  onePlusS *oneMinusT*0.125d0
        hx(6)  = -onePlusS *oneMinusT*0.125d0
        hx(7)  = -oneMinusS*oneMinusT*0.125d0
        hx(8)  =  oneMinusS*oneMinusT*0.125d0

c
c     derivadas em relacao a s :
c
        hy(1)  =  onePlusR *onePlusT*0.125d0
        hy(2)  =  oneMinusR*onePlusT*0.125d0
        hy(3)  = -oneMinusR*onePlusT*0.125d0
        hy(4)  = -onePlusR *onePlusT*0.125d0
c
        hy(5)  =  onePlusR *oneMinusT*0.125d0
        hy(6)  =  oneMinusR*oneMinusT*0.125d0
        hy(7)  = -oneMinusR*oneMinusT*0.125d0
        hy(8)  = -onePlusR *oneMinusT*0.125d0
c
c     derivadas em relacao a t :
c
        hz(1)  =  onePlusR *onePlusS *0.125d0
        hz(2)  =  oneMinusR*onePlusS *0.125d0
        hz(3)  =  oneMinusR*oneMinusS*0.125d0
        hz(4)  =  onePlusR *oneMinusS*0.125d0
c
        hz(5)  = -onePlusR *onePlusS *0.125d0
        hz(6)  = -oneMinusR*onePlusS *0.125d0
        hz(7)  = -oneMinusR*oneMinusS*0.125d0
        hz(8)  = -onePlusR *oneMinusS*0.125d0
      endif
c ......................................................................      
      return
      end
      subroutine sfhexa20_m(h,hx,hy,hz,r,s,t,afl,bfl)
c **********************************************************************
c *                                                                    *
c *                                                    22/09/15        *
c *   SFHEXA20:                                                        *
c *   -------                                                          *
c *                                                                    *
c *   Calcula as funcoes de interpolacao e suas derivadas              *
c *   no ponto (r,s,t), do hexaedro quadratico de 20 nos.              *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ---------------------                                            *
c *                                                                    *
c *     r,s,t - coordenadas naturais                                   *
c *     afl = true ( calcula funcoes )                                 *
c *     bfl = true ( calcula as derivadas )                            *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   -------------------                                              *
c *                                                                    *
c *     h(2)    - funcoes de interolocao no ponto (r,s,t)              *
c *     hx(20)  - derivadas de h em relacao a r                        *
c *     hy(20)  - derivadas de h em relacao a s                        *
c *     hz(20)  - derivadas de h em relacao a t                        *
c *                                                                    *
c * OBS:                                                               *
c * no  ( r, s, t)                                                     *
c * no1 ( 1, 1, 1)                                                     * 
c * no2 (-1, 1, 1)                                                     * 
c * no3 (-1,-1, 1)                                                     * 
c * no4 ( 1,-1, 1)                                                     * 
c * no5 ( 1, 1,-1)                                                     * 
c * no6 (-1, 1,-1)                                                     * 
c * no7 (-1,-1,-1)                                                     * 
c * no8 ( 1,-1,-1)                                                     * 
c * no9 ( 0, 1, 1)                                                     * 
c * no10(-1, 0, 1)                                                     * 
c * no11( 0,-1, 1)                                                     * 
c * no12( 1, 0, 1)                                                     * 
c * no13( 0, 1,-1)                                                     * 
c * no14(-1, 0,-1)                                                     * 
c * no15( 0,-1,-1)                                                     * 
c * no16( 1, 0,-1)                                                     * 
c * no17( 1, 1, 0)                                                     * 
c * no18(-1, 1, 0)                                                     * 
c * no19(-1,-1, 0)                                                     * 
c * no20( 1,-1, 0)                                                     * 
c **********************************************************************
      implicit none
      logical afl,bfl
      real*8 onePlusR,onePlusS,onePlusT
      real*8 oneMinusR,oneMinusS,oneMinusT
      real*8 oneMinusRr,oneMinusSs,oneMinusTt
      real*8  h(*),hx(*),hy(*),hz(*),r,s,t
c ......................................................................
c
c ...
      onePlusR  = 1.d0 + r
      onePlusS  = 1.d0 + s
      onePlusT  = 1.d0 + t
      oneMinusR = 1.d0 - r
      oneMinusS = 1.d0 - s
      oneMinusT = 1.d0 - t
      oneMinusRr= 1.d0-r*r
      oneMinusSs= 1.d0-s*s
      oneMinusTt= 1.d0-t*t
c ......................................................................
c
c ... funcoes de interpolacao
c
      if ( afl ) then
c     Nos dos Vertices
        h(1)  = onePlusR *onePlusS *onePlusT *( r+s+t-2.d0)*0.125d0
        h(2)  = oneMinusR*onePlusS *onePlusT *(-r+s+t-2.d0)*0.125d0
        h(3)  = oneMinusR*oneMinusS*onePlusT *(-r-s+t-2.d0)*0.125d0
        h(4)  = onePlusR *oneMinusS*onePlusT *( r-s+t-2.d0)*0.125d0
c
        h(5)  = onePlusR *onePlusS *oneMinusT*( r+s-t-2.d0)*0.125d0
        h(6)  = oneMinusR*onePlusS *oneMinusT*(-r+s-t-2.d0)*0.125d0
        h(7)  = oneMinusR*oneMinusS*oneMinusT*(-r-s-t-2.d0)*0.125d0
        h(8)  = onePlusR *oneMinusS*oneMinusT*( r-s-t-2.d0)*0.125d0
c     Nos pontos médios das arestas      
        h(9)  = oneMinusRr*onePlusS  *onePlusT  *0.25d0
        h(10) = oneMinusR *oneMinusSs*onePlusT  *0.25d0
        h(11) = oneMinusRr*oneMinusS *onePlusT  *0.25d0
        h(12) =  onePlusR *oneMinusSs*onePlusT  *0.25d0
c
        h(13) = oneMinusRr*onePlusS  *oneMinusT *0.25d0
        h(14) = oneMinusR *oneMinusSs*oneMinusT *0.25d0
        h(15) = oneMinusRr*oneMinusS *oneMinusT *0.25d0
        h(16) =  onePlusR *oneMinusSs*oneMinusT *0.25d0   
c 
        h(17) =  onePlusR *onePlusS  *oneMinusTt*0.25d0
        h(18) = oneMinusR *onePlusS  *oneMinusTt*0.25d0
        h(19) = oneMinusR *oneMinusS *oneMinusTt*0.25d0
        h(20) =  onePlusR *oneMinusS *oneMinusTt*0.25d0
c
      endif
c .....................................................................
c
c ... derivadas
      if ( bfl ) then
c
c     derivadas em relacao a r :
c
c     Nos dos Vertices
        hx(1)  = onePlusS *onePlusT *(2.d0*r+s+t-1.0d0)*0.125d0     
        hx(2)  = onePlusS *onePlusT *(2.d0*r-s-t+1.0d0)*0.125d0
        hx(3)  = oneMinusS*onePlusT *(2.d0*r+s-t+1.0d0)*0.125d0
        hx(4)  = oneMinusS*onePlusT *(2.d0*r-s+t-1.0d0)*0.125d0
c
        hx(5)  = onePlusS *oneMinusT*(2.d0*r+s-t-1.0d0)*0.125d0
        hx(6)  = onePlusS *oneMinusT*(2.d0*r-s+t+1.0d0)*0.125d0
        hx(7)  = oneMinusS*oneMinusT*(2.d0*r+s+t+1.0d0)*0.125d0
        hx(8)  = oneMinusS*oneMinusT*(2.d0*r-s-t-1.0d0)*0.125d0
c     Nos pontos médios das arestas      
        hx(9)  = -r*onePlusS  *onePlusT  *0.50d0
        hx(10) =   -oneMinusSs*onePlusT  *0.25d0
        hx(11) = -r*oneMinusS *onePlusT  *0.50d0
        hx(12) =    oneMinusSs*onePlusT  *0.25d0
c
        hx(13) = -r*onePlusS  *oneMinusT *0.50d0
        hx(14) =   -oneMinusSs*oneMinusT *0.25d0
        hx(15) = -r*oneMinusS *oneMinusT *0.50d0
        hx(16) =    oneMinusSs*oneMinusT *0.25d0   
c                  
        hx(17) =    onePlusS  *oneMinusTt*0.25d0                     
        hx(18) =   -onePlusS  *oneMinusTt*0.25d0                     
        hx(19) =   -oneMinusS *oneMinusTt*0.25d0                     
        hx(20) =    oneMinusS *oneMinusTt*0.25d0                     

c
c     derivadas em relacao a s :
        hy(1)  = onePlusR  *onePlusT *( r+2.d0*s+t-1.0d0)*0.125d0     
        hy(2)  = oneMinusR *onePlusT *(-r+2.d0*s+t-1.0d0)*0.125d0     
        hy(3)  = oneMinusR *onePlusT *( r+2.d0*s-t+1.0d0)*0.125d0     
        hy(4)  = onePlusR  *onePlusT *(-r+2.d0*s-t+1.0d0)*0.125d0    
c 
        hy(5)  = onePlusR  *oneMinusT*( r+2.d0*s-t-1.0d0)*0.125d0     
        hy(6)  = oneMinusR *oneMinusT*(-r+2.d0*s-t-1.0d0)*0.125d0     
        hy(7)  = oneMinusR *oneMinusT*( r+2.d0*s+t+1.0d0)*0.125d0     
        hy(8)  = onePlusR  *oneMinusT*(-r+2.d0*s+t+1.0d0)*0.125d0     
c     Nos pontos médios das arestas      
        hy(9)  =   oneMinusRr*onePlusT  *0.25d0
        hy(10) =-s*oneMinusR *onePlusT  *0.50d0
        hy(11) =  -oneMinusRr*onePlusT  *0.25d0
        hy(12) =-s*onePlusR  *onePlusT  *0.50d0
c
        hy(13) =   oneMinusRr*oneMinusT *0.25d0
        hy(14) =-s*oneMinusR *oneMinusT *0.50d0
        hy(15) =  -oneMinusRr*oneMinusT *0.25d0
        hy(16) =-s*onePlusR  *oneMinusT *0.50d0
c
        hy(17) =   onePlusR  *oneMinusTt*0.25d0
        hy(18) =   oneMinusR *oneMinusTt*0.25d0
        hy(19) =  -oneMinusR *oneMinusTt*0.25d0
        hy(20) =  -onePlusR  *oneMinusTt*0.25d0
c
c
c     derivadas em relacao a t :
c
        hz(1)  = onePlusR *onePlusS  *( r+s+2.d0*t-1.0d0)*0.125d0     
        hz(2)  = oneMinusR*onePlusS  *(-r+s+2.d0*t-1.0d0)*0.125d0     
        hz(3)  = oneMinusR*oneMinusS *(-r-s+2.d0*t-1.0d0)*0.125d0     
        hz(4)  = onePlusR *oneMinusS *( r-s+2.d0*t-1.0d0)*0.125d0     
c
        hz(5)  = onePlusR *onePlusS  *(-r-s+2.d0*t+1.0d0)*0.125d0     
        hz(6)  = oneMinusR*onePlusS  *( r-s+2.d0*t+1.0d0)*0.125d0     
        hz(7)  = oneMinusR*oneMinusS *( r+s+2.d0*t+1.0d0)*0.125d0     
        hz(8)  = onePlusR *oneMinusS *(-r+s+2.d0*t+1.0d0)*0.125d0     
c     Nos pontos médios das arestas      
        hz(9)  = oneMinusRr *onePlusS  *0.25d0
        hz(10) = oneMinusR  *oneMinusSs*0.25d0
        hz(11) = oneMinusRr *oneMinusS *0.25d0
        hz(12) =  onePlusR  *oneMinusSs*0.25d0
c
        hz(13) = -oneMinusRr*onePlusS  *0.25d0
        hz(14) = -oneMinusR *oneMinusSs*0.25d0
        hz(15) = -oneMinusRr*oneMinusS *0.25d0
        hz(16) =  -onePlusR *oneMinusSs*0.25d0 
c   
        hz(17) =-t*onePlusR *onePlusS  *0.50d0
        hz(18) =-t*oneMinusR*onePlusS  *0.50d0
        hz(19) =-t*oneMinusR*oneMinusS *0.50d0
        hz(20) =-t*onePlusR *oneMinusS *0.50d0
      endif
c ......................................................................      
      return
      end
c **********************************************************************
      subroutine sfhexa20(h,hx,hy,hz,r,s,t,afl,bfl)
c **********************************************************************
c *                                                                    *
c *                                                    27/01/11        *
c *   SFHEXA20:                                                        *
c *   -------                                                          *
c *                                                                    *
c *   Calcula as funcoes de interpolacao e suas derivadas              *
c *   no ponto (r,s,t), do hexaedro quadratico de 20 nos.              *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ---------------------                                            *
c *                                                                    *
c *     r,s,t - coordenadas naturais                                   *
c *     afl = true ( calcula funcoes )                                 *
c *     bfl = true ( calcula as derivadas )                            *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   -------------------                                              *
c *                                                                    *
c *     h(8)    - funcoes de interolocao no ponto (r,s,t)              *
c *     hx(8)   - derivadas de h em relacao a r                        *
c *     hy(8)   - derivadas de h em relacao a s                        *
c *     hz(8)   - derivadas de h em relacao a t                        *
c *                                                                    *
c **********************************************************************
      implicit none
      logical afl,bfl
      real*8  h(*),hx(*),hy(*),hz(*),r,s,t
c ......................................................................
      if ( afl ) then
c
c     funcoes de interpolacao
c
c     Nos dos Vertices
      h(1)  = ( 1.d0+r ) * ( 1.d0+s ) * (1.d0+t) *(r+s+t-2)/ 8.d0
      h(2)  = ( 1.d0-r ) * ( 1.d0+s ) * (1.d0+t) *(-r+s+t-2)/ 8.d0
      h(3)  = ( 1.d0-r ) * ( 1.d0-s ) * (1.d0+t) *(-r-s+t-2)/ 8.d0
      h(4)  = ( 1.d0+r ) * ( 1.d0-s ) * (1.d0+t) *(r-s+t-2)/ 8.d0
c
      h(5)  = ( 1.d0+r ) * ( 1.d0+s ) * (1.d0-t) *(r+s-t-2)/ 8.d0
      h(6)  = ( 1.d0-r ) * ( 1.d0+s ) * (1.d0-t) *(-r+s-t-2)/ 8.d0
      h(7)  = ( 1.d0-r ) * ( 1.d0-s ) * (1.d0-t) *(-r-s-t-2)/ 8.d0
      h(8)  = ( 1.d0+r ) * ( 1.d0-s ) * (1.d0-t) *(r-s-t-2)/ 8.d0
c     Nos pontos médios das arestas      
      h(9)  =  ( 1.d0-r**2 ) * ( 1.d0+s ) * (1.d0+t) / 4.d0
      h(10)  = ( 1.d0-r ) * ( 1.d0-s**2 ) * (1.d0+t) / 4.d0
      h(11)  = ( 1.d0-r**2 ) * ( 1.d0-s ) * (1.d0+t) / 4.d0
      h(12)  = ( 1.d0+r ) * ( 1.d0-s**2 ) * (1.d0+t) / 4.d0
c
      h(13)  =  ( 1.d0-r**2 ) * ( 1.d0+s ) * (1.d0-t) / 4.d0
      h(14)  =  ( 1.d0-r ) * ( 1.d0-s**2 ) * (1.d0-t) / 4.d0
      h(15)  =  ( 1.d0-r**2 ) * ( 1.d0-s ) * (1.d0-t) / 4.d0
      h(16)  =  ( 1.d0+r ) * ( 1.d0-s**2 ) * (1.d0-t) / 4.d0    
c
      h(17)  = ( 1.d0+r ) * ( 1.d0+s ) * (1.d0-t**2) / 4.d0
      h(18)  = ( 1.d0-r ) * ( 1.d0+s ) * (1.d0-t**2) / 4.d0
      h(19)  = ( 1.d0-r ) * ( 1.d0-s ) * (1.d0-t**2) / 4.d0
      h(20)  = ( 1.d0+r ) * ( 1.d0-s ) * (1.d0-t**2) / 4.d0
c
      endif
      if ( bfl ) then
c
c     derivadas em relacao a r :
c
      hx(1) = ( 1.d0+s ) * (1.d0+t) *(r+s+t-2)/ 8.d0+
     . ( 1.d0+r ) * ( 1.d0+s ) * (1.d0+t) / 8.d0     
      hx(2) = -( 1.d0+s ) * (1.d0+t) *(-r+s+t-2)/ 8.d0-
     .( 1.d0-r ) * ( 1.d0+s ) * (1.d0+t) / 8.d0
      hx(3) = - ( 1.d0-s ) * (1.d0+t) *(-r-s+t-2)/ 8.d0-
     .( 1.d0-r ) * ( 1.d0-s ) * (1.d0+t) / 8.d0
      hx(4) =  ( 1.d0-s ) * (1.d0+t) *(r-s+t-2)/ 8.d0+
     .( 1.d0+r ) * ( 1.d0-s ) * (1.d0+t) / 8.d0
c
      hx(5) = ( 1.d0+s ) * (1.d0-t) *(r+s-t-2)/ 8.d0+
     . ( 1.d0+r ) * ( 1.d0+s ) * (1.d0-t) / 8.d0
      hx(6) = -( 1.d0+s ) * (1.d0-t) *(-r+s-t-2)/ 8.d0-
     .( 1.d0-r ) * ( 1.d0+s ) * (1.d0-t)/ 8.d0
      hx(7) = -( 1.d0-s ) * (1.d0-t) *(-r-s-t-2)/ 8.d0-
     .( 1.d0-r ) * ( 1.d0-s ) * (1.d0-t) / 8.d0
      hx(8) =  ( 1.d0-s ) * (1.d0-t) *(r-s-t-2)/ 8.d0+
     .( 1.d0+r ) * ( 1.d0-s ) * (1.d0-t)/ 8.d0
c
      hx(9)  =  -2*r * ( 1.d0+s ) * (1.d0+t) / 4.d0
      hx(10)  = -( 1.d0-s**2 ) * (1.d0+t) / 4.d0
      hx(11)  = -2*r * ( 1.d0-s ) * (1.d0+t) / 4.d0
      hx(12)  =  ( 1.d0-s**2 ) * (1.d0+t) / 4.d0
c
      hx(13)  = -2*r * ( 1.d0+s ) * (1.d0-t) / 4.d0
      hx(14)  = -( 1.d0-s**2 ) * (1.d0-t) / 4.d0
      hx(15)  = -2*r * ( 1.d0-s ) * (1.d0-t) / 4.d0
      hx(16)  =  ( 1.d0-s**2 ) * (1.d0-t) / 4.d0  
c
      hx(17)  =  ( 1.d0+s ) * (1.d0-t**2) / 4.d0
      hx(18)  = -( 1.d0+s ) * (1.d0-t**2) / 4.d0
      hx(19)  = -( 1.d0-s ) * (1.d0-t**2) / 4.d0
      hx(20)  =  ( 1.d0-s ) * (1.d0-t**2) / 4.d0
c
c     derivadas em relacao a s :
c
      hy(1)  = ( 1.d0+r ) *  (1.d0+t) *(r+s+t-2)/ 8.d0+
     . ( 1.d0+r ) * ( 1.d0+s ) * (1.d0+t) / 8.d0
      hy(2)  = ( 1.d0-r ) * (1.d0+t) *(-r+s+t-2)/ 8.d0+
     .( 1.d0-r ) * ( 1.d0+s ) * (1.d0+t) / 8.d0
      hy(3)  = -( 1.d0-r ) *  (1.d0+t) *(-r-s+t-2)/ 8.d0-
     .( 1.d0-r ) * ( 1.d0-s ) * (1.d0+t) / 8.d0
      hy(4)  = -( 1.d0+r ) *  (1.d0+t) *(r-s+t-2)/ 8.d0-
     .( 1.d0+r ) * ( 1.d0-s ) * (1.d0+t)/ 8.d0
c
      hy(5)  = ( 1.d0+r ) * (1.d0-t) *(r+s-t-2)/ 8.d0+
     .( 1.d0+r ) * ( 1.d0+s ) * (1.d0-t) / 8.d0
      hy(6)  = ( 1.d0-r ) * (1.d0-t) *(-r+s-t-2)/ 8.d0+
     .( 1.d0-r ) * ( 1.d0+s ) * (1.d0-t) / 8.d0
      hy(7)  = -( 1.d0-r ) * (1.d0-t) *(-r-s-t-2)/ 8.d0-
     .( 1.d0-r ) * ( 1.d0-s ) * (1.d0-t) / 8.d0
      hy(8)  = -( 1.d0+r ) *  (1.d0-t) *(r-s-t-2)/ 8.d0-
     .( 1.d0+r ) * ( 1.d0-s ) * (1.d0-t) / 8.d0
c
      hy(9)  =  ( 1.d0-r**2 ) * (1.d0+t) / 4.d0
      hy(10)  = ( 1.d0-r ) * ( -2*s ) * (1.d0+t) / 4.d0
      hy(11)  = -( 1.d0-r**2 )  * (1.d0+t) / 4.d0
      hy(12)  = ( 1.d0+r ) * (-2*s ) * (1.d0+t) / 4.d0
c
      hy(17)  = ( 1.d0+r )  * (1.d0-t**2) / 4.d0
      hy(18)  = ( 1.d0-r )  * (1.d0-t**2) / 4.d0
      hy(19)  = -( 1.d0-r )  * (1.d0-t**2) / 4.d0
      hy(20)  = -( 1.d0+r )  * (1.d0-t**2) / 4.d0
c
      hy(13)  =  ( 1.d0-r**2 ) * (1.d0-t) / 4.d0
      hy(14)  =  ( 1.d0-r ) * ( -2*s ) * (1.d0-t) / 4.d0
      hy(15)  = -( 1.d0-r**2 )  * (1.d0-t) / 4.d0
      hy(16)  =  ( 1.d0+r ) * ( -2*s ) * (1.d0-t) / 4.d0  
c
c     derivadas em relacao a t :
c
      hz(1)  = ( 1.d0+r ) * ( 1.d0+s ) * (r+s+t-2)/ 8.d0+
     .( 1.d0+r ) * ( 1.d0+s ) * (1.d0+t) / 8.d0
      hz(2)  = ( 1.d0-r ) * ( 1.d0+s ) *(-r+s+t-2)/ 8.d0+
     .( 1.d0-r ) * ( 1.d0+s ) * (1.d0+t)/ 8.d0
      hz(3)  = ( 1.d0-r ) * ( 1.d0-s ) *(-r-s+t-2)/ 8.d0+
     . ( 1.d0-r ) * ( 1.d0-s ) * (1.d0+t) / 8.d0
      hz(4)  = ( 1.d0+r ) * ( 1.d0-s ) * (r-s+t-2)/ 8.d0+
     .( 1.d0+r ) * ( 1.d0-s ) * (1.d0+t) / 8.d0
c
      hz(5)  = -( 1.d0+r ) * ( 1.d0+s ) * (r+s-t-2)/ 8.d0-
     .( 1.d0+r ) * ( 1.d0+s ) * (1.d0-t)/ 8.d0
      hz(6)  = -( 1.d0-r ) * ( 1.d0+s ) * (-r+s-t-2)/ 8.d0-
     .( 1.d0-r ) * ( 1.d0+s ) * (1.d0-t)/ 8.d0
      hz(7)  = -( 1.d0-r ) * ( 1.d0-s ) * (-r-s-t-2)/ 8.d0-
     .( 1.d0-r ) * ( 1.d0-s ) * (1.d0-t) / 8.d0
      hz(8)  = -( 1.d0+r ) * ( 1.d0-s )  *(r-s-t-2)/ 8.d0-
     .( 1.d0+r ) * ( 1.d0-s ) * (1.d0-t) / 8.d0
c
      hz(9)  =  ( 1.d0-r**2 ) * ( 1.d0+s )  / 4.d0
      hz(10)  = ( 1.d0-r ) * ( 1.d0-s**2 )  / 4.d0
      hz(11)  = ( 1.d0-r**2 ) * ( 1.d0-s )  / 4.d0
      hz(12)  = ( 1.d0+r ) * ( 1.d0-s**2 )  / 4.d0
c
      hz(13)  =  -( 1.d0-r**2 ) * ( 1.d0+s )  / 4.d0
      hz(14)  =  -( 1.d0-r ) * ( 1.d0-s**2 )  / 4.d0
      hz(15)  =  -( 1.d0-r**2 ) * ( 1.d0-s )  / 4.d0
      hz(16)  =  -( 1.d0+r ) * ( 1.d0-s**2 )  / 4.d0 
c
      hz(17)  = ( 1.d0+r ) * ( 1.d0+s ) * (-2*t) / 4.d0
      hz(18)  = ( 1.d0-r ) * ( 1.d0+s ) * (-2*t) / 4.d0
      hz(19)  = ( 1.d0-r ) * ( 1.d0-s ) * (-2*t) / 4.d0
      hz(20)  = ( 1.d0+r ) * ( 1.d0-s ) * (-2*t) / 4.d0
      endif
c ......................................................................      
      return
      end
c **********************************************************************
      subroutine sfwedge(h,hx,hy,hz,r,s,t,afl,bfl)
c **********************************************************************
c *                                                                    *
c *                                                    26/03/03        *
c *   SFWEDGE:                                                         *
c *   -------                                                          *
c *                                                                    *
c *   Calcula as funcoes de interpolacao e suas derivadas              *
c *   no ponto (r,s,t) para o elemento cunha                           *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ---------------------                                            *
c *                                                                    *
c *     r,s,t - coordenadas naturais                                   *
c *     afl = true ( calcula funcoes )                                 *
c *     bfl = true ( calcula as derivadas )                            *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   -------------------                                              *
c *                                                                    *
c *     h(6)    - funcoes de interolocao no ponto (r,s,t)              *
c *     hr(6)   - derivadas de h em relacao a r                        *
c *     hs(6)   - derivadas de h em relacao a s                        *
c *                                                                    *
c **********************************************************************
      implicit none
      real*8  h(*),hx(*),hy(*),hz(*),r,s,t
      logical afl,bfl
c ......................................................................
      if (afl) then
c
c ...... Funcoes de interpolacao:
c
         h(1) = r * (1.d0+t) * 0.5d0
         h(2) = s * (1.d0+t) * 0.5d0
         h(3) =(1.d0-r-s) * (1.d0+t) * 0.5d0
         h(4) = r * (1.d0-t) * 0.5d0
         h(5) = s * (1.d0-t) * 0.5d0
         h(6) =(1.d0-r-s) * (1.d0-t) * 0.5d0
      endif
      if (bfl) then
c
c ...... Derivadas em relacao a r:
c
         hx(1) = (1.d0+t) * 0.5d0
         hx(2) =  0.d0
         hx(3) =-(1.d0+t) * 0.5d0
         hx(4) = (1.d0-t) * 0.5d0
         hx(5) =  0.d0
         hx(6) = (t-1.d0) * 0.5d0
c
c ...... Derivadas em relacao a s:
c
         hy(1) =  0.d0
         hy(2) = (1.d0+t) * 0.5d0
         hy(3) =-(1.d0+t) * 0.5d0
         hy(4) =  0.d0
         hy(5) = (1.d0-t) * 0.5d0
         hy(6) = (t-1.d0) * 0.5d0
c
c ...... Derivadas em relacao a t:
c
         hz(1) =  r * 0.5d0
         hz(2) =  s * 0.5d0
         hz(3) = (1.d0-r-s) * 0.5d0
         hz(4) = -r * 0.5d0
         hz(5) = -s * 0.5d0
         hz(6) = (r+s-1.d0) * 0.5d0
      endif
c ......................................................................      
      return
      end
      subroutine jacob2d(hx,hy,xj,xji,x,det,nen,ndm,nel)
c **********************************************************************
c *                                                                    *
c *                                                    20/03/03        *
c *                                                                    *
c *   JACOB2D: determinante do Jacobiano no ponto (r,s)                *
c *   -------                                                          *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ---------------------                                            *
c *                                                                    *
c *     hx(nen)   - derivadas de h em relacao a r                      *
c *     hy(nen)   - derivadas de h em relacao a s                      *
c *     x(ndm,nen)- coordenadas nodais do elemento                     *
c *     nen       - numero de nos do elemento                          *
c *     nel       - numero do elemento                                 *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   -------------------                                              *
c *                                                                    *
c *     hx(nen)  - derivadas de h em relacao a x                       *
c *     hy(nen)  - derivadas de h em relacao a y                       *
c *     xj(2,2)  - matriz jacobiana no ponto (r,s,t)                   *
c *     xji(2,2) - inversa de xj                                       *
c *     det      - determinante da matriz jacobiana no ponto (r,s)     *
c *                                                                    *
c *                                                                    *
c **********************************************************************
      implicit none
      integer nen,ndm,nel,j,k
      real*8  ZERO
      real*8  hx(*),hy(*),xj(2,*),xji(2,*),x(ndm,*)
      real*8  det,hxk,hyk
      parameter (ZERO = 1.d-14)
c ......................................................................
c
c ... Matriz Jacobiana:
c
      do 200 j = 1 , 2
         xj(1,j) = 0.d0
         xj(2,j) = 0.d0
         do 100 k = 1 , nen
            xj(1,j) = xj(1,j) + hx(k) * x(j,k)
            xj(2,j) = xj(2,j) + hy(k) * x(j,k)
  100    continue
  200 continue
c
c ... Determinante da matriz Jacobiana:  
c
c ......................................................................              
      det = xj(1,1)*xj(2,2)-xj(2,1)*xj(1,2)
      if (det .le. ZERO) then        
         print*,'*** Subrotina ELMT: determinante <= 0 ',nel
         stop
      endif
c ......................................................................                    
c
c ... Inversa da matriz Jacobiana:  
c
      xji(1,1) =  xj(2,2) / det
      xji(1,2) = -xj(1,2) / det
      xji(2,1) = -xj(2,1) / det
      xji(2,2) =  xj(1,1) / det
c
c ... Derivadas das funcoes de interpolacao:
c
      do 300 k = 1, nen
         hxk = hx(k)
         hyk = hy(k)
         hx(k) = xji(1,1)*hxk + xji(1,2)*hyk
         hy(k) = xji(2,1)*hxk + xji(2,2)*hyk
  300 continue
c ......................................................................                
      return
      end
      subroutine jacob2d_m(hx,hy,xj,xji,x,det,nen,ndm,nel,afl)
c **********************************************************************
c *                                                                    *
c *                                                    20/11/15        *
c *                                                                    *
c *   JACOB2D: determinante do Jacobiano no ponto (r,s)                *
c *   -------                                                          *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ---------------------                                            *
c *                                                                    *
c *     hx(nen)   - derivadas de h em relacao a r                      *
c *     hy(nen)   - derivadas de h em relacao a s                      *
c *     x(ndm,nen)- coordenadas nodais do elemento                     *
c *     nen       - numero de nos do elemento                          *
c *     nel       - numero do elemento                                 *
c *     afl       - true - calcula a inversa da matriz Jacobiana       *
c *                 false- calcula somente as derivadas das funcoes de *
c *                 interpolacao                                       *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   -------------------                                              *
c *                                                                    *
c *     hx(nen)  - derivadas de h em relacao a x                       *
c *     hy(nen)  - derivadas de h em relacao a y                       *
c *     xj(2,2)  - matriz jacobiana no ponto (r,s,t)                   *
c *     xji(2,2) - inversa de xj                                       *
c *     det      - determinante da matriz jacobiana no ponto (r,s)     *
c *                                                                    *
c *                                                                    *
c **********************************************************************
      implicit none
      integer nen,ndm,nel,j,k
      real*8  ZERO
      real*8  hx(*),hy(*),xj(2,*),xji(2,*),x(ndm,*)
      real*8  det,hxk,hyk
      logical afl
      parameter (ZERO = 1.d-14)
c ......................................................................
c
c ...
      if (afl) then 
c ... Matriz Jacobiana:
c
        do 200 j = 1 , 2
          xj(1,j) = 0.d0
          xj(2,j) = 0.d0
          do 100 k = 1 , nen
            xj(1,j) = xj(1,j) + hx(k) * x(j,k)
            xj(2,j) = xj(2,j) + hy(k) * x(j,k)
  100     continue
  200   continue
c
c ... Determinante da matriz Jacobiana:  
c
c ......................................................................
        det = xj(1,1)*xj(2,2)-xj(2,1)*xj(1,2)
        if (det .le. ZERO) then        
          print*,'*** Subrotina ELMT: determinante <= 0 ',nel
          stop
        endif
c ......................................................................                    
c
c ... Inversa da matriz Jacobiana:  
c
        xji(1,1) =  xj(2,2) / det
        xji(1,2) = -xj(1,2) / det
        xji(2,1) = -xj(2,1) / det
        xji(2,2) =  xj(1,1) / det
      endif
c ......................................................................                    
c
c ... Derivadas das funcoes de interpolacao:
c
      do 300 k = 1, nen
         hxk = hx(k)
         hyk = hy(k)
         hx(k) = xji(1,1)*hxk + xji(1,2)*hyk
         hy(k) = xji(2,1)*hxk + xji(2,2)*hyk
  300 continue
c ......................................................................                
      return
      end
c **********************************************************************
      subroutine jacob3d(hx,hy,hz,xj,xji,x,det,nen,nel)
c **********************************************************************
c *                                                                    *
c *                                                    10/10/01        *
c *                                                                    *
c *   JACOB3D: determinante do Jacobiano no ponto (r,s,t)              *
c *   -------                                                          *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ---------------------                                            *
c *                                                                    *
c *     hx(nen)   - derivadas de h em relacao a r                      *
c *     hy(nen)   - derivadas de h em relacao a s                      *
c *     hz(nen)   - derivadas de h em relacao a s                      *
c *     x(ndm,nen)- coordenadas nodais do elemento                     *
c *     nen       - numero de nos do elemento                          *
c *     nel       - numero do elemento                                 *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   -------------------                                              *
c *                                                                    *
c *     hx(nen)  - derivadas de h em relacao a x                       *
c *     hy(nen)  - derivadas de h em relacao a y                       *
c *     hz(nen)  - derivadas de h em relacao a z                       *
c *     xj(3,3)  - matriz jacobiana no ponto (r,s,t)                   *
c *     xji(3,3) - inversa de xj                                       *
c *     det      - determinante da matriz jacobiana no ponto (r,s,t)   *
c *                                                                    *
c *                                                                    *
c **********************************************************************
      implicit none
      integer nen,nel,j,k
      real*8  hx(*),hy(*),hz(*),xj(3,*),xji(3,*),x(3,*)
      real*8  det,hxk,hyk,hzk
      real*8  ZERO
      parameter (ZERO = 1.d-14)
c ......................................................................
c
c ... Matria Jacobiana:
c
      do 200 j = 1 , 3
         xj(1,j) = 0.d0
         xj(2,j) = 0.d0
         xj(3,j) = 0.d0
         do 100 k = 1 , nen
            xj(1,j) = xj(1,j) + hx(k) * x(j,k)
            xj(2,j) = xj(2,j) + hy(k) * x(j,k)
            xj(3,j) = xj(3,j) + hz(k) * x(j,k)
  100    continue
  200 continue
c
c ... Determinante da matriz Jacobiana:  
c
      det = xj(1,1)*xj(2,2)*xj(3,3) + xj(1,2)*xj(2,3)*xj(3,1) + xj(1,3)*
     .      xj(2,1)*xj(3,2) - xj(3,1)*xj(2,2)*xj(1,3) - xj(2,1)*xj(1,2)*
     .      xj(3,3) - xj(1,1)*xj(3,2)*xj(2,3)
c ......................................................................                   
      if (det .le. ZERO) then
         print*,'*** Subrotina ELMT__: determinante <= 0 ',nel
         stop
      endif
c ......................................................................                    
c
c ... Inversa da matriz Jacobiana:  
c
      xji(1,1) =  ( xj(2,2) * xj(3,3) - xj(2,3) * xj(3,2) ) / det
      xji(2,1) = -( xj(2,1) * xj(3,3) - xj(2,3) * xj(3,1) ) / det
      xji(3,1) =  ( xj(2,1) * xj(3,2) - xj(2,2) * xj(3,1) ) / det
      xji(1,2) = -( xj(1,2) * xj(3,3) - xj(1,3) * xj(3,2) ) / det
      xji(2,2) =  ( xj(1,1) * xj(3,3) - xj(1,3) * xj(3,1) ) / det
      xji(3,2) = -( xj(1,1) * xj(3,2) - xj(1,2) * xj(3,1) ) / det
      xji(1,3) =  ( xj(1,2) * xj(2,3) - xj(1,3) * xj(2,2) ) / det
      xji(2,3) = -( xj(1,1) * xj(2,3) - xj(1,3) * xj(2,1) ) / det
      xji(3,3) =  ( xj(1,1) * xj(2,2) - xj(1,2) * xj(2,1) ) / det
c
c ... Derivadas das funcoes de interpolacao:
c
      do 300 k = 1, nen
         hxk = hx(k)
         hyk = hy(k)
         hzk = hz(k)
         hx(k) = xji(1,1)*hxk + xji(1,2)*hyk + xji(1,3)*hzk
         hy(k) = xji(2,1)*hxk + xji(2,2)*hyk + xji(2,3)*hzk
         hz(k) = xji(3,1)*hxk + xji(3,2)*hyk + xji(3,3)*hzk
  300 continue
c ......................................................................              
      return
      end
c **********************************************************************
      subroutine jacob3d_m(hx,hy,hz,xj,xji,x,det,nen,nel,afl)
c **********************************************************************
c *                                                                    *
c *                                                    10/10/15        *
c *                                                                    *
c *   JACOB3D: determinante do Jacobiano no ponto (r,s,t)              *
c *   -------                                                          *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ---------------------                                            *
c *                                                                    *
c *     hx(nen)   - derivadas de h em relacao a r                      *
c *     hy(nen)   - derivadas de h em relacao a s                      *
c *     hz(nen)   - derivadas de h em relacao a s                      *
c *     x(ndm,nen)- coordenadas nodais do elemento                     *
c *     nen       - numero de nos do elemento                          *
c *     nel       - numero do elemento                                 *
c *     afl       - true - calcula a inversa da matriz Jacobiana       *
c *                 false- calcula somente as derivadas das funcoes de *
c *                 interpolacao                                       *
c *   Parametros de saida:                                             *
c *   -------------------                                              *
c *                                                                    *
c *     hx(nen)  - derivadas de h em relacao a x                       *
c *     hy(nen)  - derivadas de h em relacao a y                       *
c *     hz(nen)  - derivadas de h em relacao a z                       *
c *     xj(3,3)  - matriz jacobiana no ponto (r,s,t)                   *
c *     xji(3,3) - inversa de xj                                       *
c *     det      - determinante da matriz jacobiana no ponto (r,s,t)   *
c *                                                                    *
c *                                                                    *
c **********************************************************************
      implicit none
      integer nen,nel,j,k
      real*8  hx(*),hy(*),hz(*),xj(3,*),xji(3,*),x(3,*)
      real*8  det,hxk,hyk,hzk
      logical afl
      real*8  ZERO
      parameter (ZERO = 1.d-14)
c ......................................................................
c
c ...
      if (afl) then 
c ... Matria Jacobiana:
c
        do 200 j = 1 , 3
          xj(1,j) = 0.d0
          xj(2,j) = 0.d0
          xj(3,j) = 0.d0
          do 100 k = 1 , nen
            xj(1,j) = xj(1,j) + hx(k) * x(j,k)
            xj(2,j) = xj(2,j) + hy(k) * x(j,k)
            xj(3,j) = xj(3,j) + hz(k) * x(j,k)
  100     continue
  200   continue
c
c ... Determinante da matriz Jacobiana:  
c
        det = xj(1,1)*xj(2,2)*xj(3,3) + xj(1,2)*xj(2,3)*xj(3,1)+xj(1,3)*
     .        xj(2,1)*xj(3,2) - xj(3,1)*xj(2,2)*xj(1,3)-xj(2,1)*xj(1,2)*
     .        xj(3,3) - xj(1,1)*xj(3,2)*xj(2,3)
c ......................................................................                   
        if (det .le. ZERO) then
          print*,'*** Subrotina ELMT__: determinante <= 0 ',nel
          stop
        endif
c ......................................................................                    
c     
c ... Inversa da matriz Jacobiana:  
c
        xji(1,1) =  ( xj(2,2) * xj(3,3) - xj(2,3) * xj(3,2) ) / det
        xji(2,1) = -( xj(2,1) * xj(3,3) - xj(2,3) * xj(3,1) ) / det
        xji(3,1) =  ( xj(2,1) * xj(3,2) - xj(2,2) * xj(3,1) ) / det
        xji(1,2) = -( xj(1,2) * xj(3,3) - xj(1,3) * xj(3,2) ) / det
        xji(2,2) =  ( xj(1,1) * xj(3,3) - xj(1,3) * xj(3,1) ) / det
        xji(3,2) = -( xj(1,1) * xj(3,2) - xj(1,2) * xj(3,1) ) / det
        xji(1,3) =  ( xj(1,2) * xj(2,3) - xj(1,3) * xj(2,2) ) / det
        xji(2,3) = -( xj(1,1) * xj(2,3) - xj(1,3) * xj(2,1) ) / det
        xji(3,3) =  ( xj(1,1) * xj(2,2) - xj(1,2) * xj(2,1) ) / det
c ......................................................................          
      endif
c ......................................................................          
c
c ... Derivadas das funcoes de interpolacao:
c
      do 300 k = 1, nen
         hxk = hx(k)
         hyk = hy(k)
         hzk = hz(k)
         hx(k) = xji(1,1)*hxk + xji(1,2)*hyk + xji(1,3)*hzk
         hy(k) = xji(2,1)*hxk + xji(2,2)*hyk + xji(2,3)*hzk
         hz(k) = xji(3,1)*hxk + xji(3,2)*hyk + xji(3,3)*hzk
  300 continue
c ......................................................................              
      return
      end
      subroutine rotmatrix(x,r)
c **********************************************************************
c *                                                                    *
c *                                                    26/03/03        *
c *   ROTACAO:                                                         *
c *   -------                                                          *
c *                                                                    *
c **********************************************************************
      implicit none
      real*8  x(3,*),r(3,3)
      real*8 v1x,v1y,v1z,v2x,v2y,v2z,v3x,v3y,v3z,v1,v2,v3
c ......................................................................
c
c ... v1 (x2-x1)     
      v1x = x(1,2)-x(1,1)
      v1y = x(2,2)-x(2,1)
      v1z = x(3,2)-x(3,1)
c ... v2 (x3-x1)
      v2x = x(1,3)-x(1,1)
      v2y = x(2,3)-x(2,1)
      v2z = x(3,3)-x(3,1)
c ... v3 = v1 x v2
      v3x = v1y*v2z - v1z*v2y
      v3y = v1z*v2x - v1x*v2z
      v3z = v1x*v2y - v1y*v2x
c ... v2 = v1 x v3
      v2x = v1z*v3y - v1y*v3z
      v2y = v1x*v3z - v1z*v3x
      v2z = v1y*v3x - v1x*v3y
c ... norm
      v1 = dsqrt(v1x*v1x+v1y*v1y+v1z*v1z)
      v2 = dsqrt(v2x*v2x+v2y*v2y+v2z*v2z)
      v3 = dsqrt(v3x*v3x+v3y*v3y+v3z*v3z)
c ...
c     print*,v1x/v1,v1y/v1,v1z/v1
c     print*,v2x/v2,v2y/v2,v2z/v2
c     print*,v3x/v3,v3y/v3,v3z/v3
c ...
      r(1,1) = v1x/v1
      r(1,2) = v1y/v1
      r(1,3) = v1z/v1
      r(2,1) = v2x/v2
      r(2,2) = v2y/v2
      r(2,3) = v2z/v2
      r(3,1) = v3x/v3
      r(3,2) = v3y/v3
      r(3,3) = v3z/v3
c ......................................................................      
      return
      end
      subroutine rotmatrix2(x1,x2,x3,r)
c **********************************************************************
c *                                                                    *
c *                                                    26/03/03        *
c *   MATRIZ DE ROTACAO: 3D                                            *
c *   ------------------                                               *
c *                                                                    *
c **********************************************************************
      implicit none
      real*8  x1(3),x2(3),x3(3),r(3,3)
      real*8 v1x,v1y,v1z,v2x,v2y,v2z,v3x,v3y,v3z,v1,v2,v3
c ......................................................................
c
c ... v1 (x2-x1)     
      v1x = x2(1)-x1(1)
      v1y = x2(2)-x1(2)
      v1z = x2(3)-x1(3)
c ... v2 (x3-x1)
      v2x = x3(1)-x1(1)
      v2y = x3(2)-x1(2)
      v2z = x3(3)-x1(3)
c ... v3 = v1 x v2
      v3x = v1y*v2z - v1z*v2y
      v3y = v1z*v2x - v1x*v2z
      v3z = v1x*v2y - v1y*v2x
c ... v2 = v1 x v3
      v2x = v1z*v3y - v1y*v3z
      v2y = v1x*v3z - v1z*v3x
      v2z = v1y*v3x - v1x*v3y
c ... norm
      v1 = dsqrt(v1x*v1x+v1y*v1y+v1z*v1z)
      v2 = dsqrt(v2x*v2x+v2y*v2y+v2z*v2z)
      v3 = dsqrt(v3x*v3x+v3y*v3y+v3z*v3z)
c ...
      r(1,1) = v1x/v1
      r(1,2) = v1y/v1
      r(1,3) = v1z/v1
      r(2,1) = v2x/v2
      r(2,2) = v2y/v2
      r(2,3) = v2z/v2
      r(3,1) = v3x/v3
      r(3,2) = v3y/v3
      r(3,3) = v3z/v3
c ......................................................................      
      return
      end
      subroutine irotmatrix(x,r)
c **********************************************************************
c *                                                                    *
c *                                                    21/11/2010      *
c *   ROTACAO:                                                         *
c *   calcula matriz de rotacao sistema local para 2D                  *
c *                                                                    *
c **********************************************************************
      implicit none
      real*8  x(2,*),r(2,2)
      real*8 v1x,v1y,v1z,v2x,v2y,v2z,v3x,v3y,v3z,v1,v2,v3
c ......................................................................
      v1x = x(1,2)-x(1,1)
      v1y = x(2,2)-x(2,1)
      v2x = -v1y
      v2y = v1x
      v1 = dsqrt(v1x*v1x+v1y*v1y)
      v2 = dsqrt(v2x*v2x+v2y*v2y)
      r(1,1) = v1x/v1
      r(1,2) = v1y/v1
      r(2,1) = v2x/v2
      r(2,2) = v2y/v2
c ......................................................................      
      return
      end
      subroutine rotate(x,r,xl,transp)
c **********************************************************************
c *                                                                    *
c *                                                    26/03/03        *
c *   ROTACAO:                                                         *
c *   -------                                                          *
c *                                                                    *
c **********************************************************************
      implicit none
      real*8  x(3),r(3,3),xl(3),x1,x2,x3
      logical transp
c ......................................................................
      if(transp) then
         x1 = r(1,1)*x(1)+r(2,1)*x(2)+r(3,1)*x(3)
         x2 = r(1,2)*x(1)+r(2,2)*x(2)+r(3,2)*x(3)
         x3 = r(1,3)*x(1)+r(2,3)*x(2)+r(3,3)*x(3)      
      else
         x1 = r(1,1)*x(1)+r(1,2)*x(2)+r(1,3)*x(3)
         x2 = r(2,1)*x(1)+r(2,2)*x(2)+r(2,3)*x(3)
         x3 = r(3,1)*x(1)+r(3,2)*x(2)+r(3,3)*x(3)
      endif
      xl(1) = x1
      xl(2) = x2
      xl(3) = x3
c ......................................................................      
      return
      end
      subroutine irotate(x,r,xl,transp)
c **********************************************************************
c *                                                                    *
c *                                                   21/11/2010       *
c *   ROTACAO:                                                         *
c *   Rotaciona as coordenadas                                         *
c *                                                                    *
c **********************************************************************
      implicit none
      real*8  x(2),r(2,2),xl(2),x1,x2
      logical transp
c ......................................................................
      if(transp) then
         x1 = r(1,1)*x(1)+r(2,1)*x(2)
         x2 = r(1,2)*x(1)+r(2,2)*x(2)
      else
         x1 = r(1,1)*x(1)+r(1,2)*x(2)
         x2 = r(2,1)*x(1)+r(2,2)*x(2)
      endif
      xl(1) = x1
      xl(2) = x2
c ......................................................................      
      return
      end
      subroutine irotate2(x,r,xl,inv)
c **********************************************************************
c *                                                                    *
c *                                                   21/11/2010       *
c *   ROTACAO:                                                         *
c *   Rotaciona tensor de segunda ordem                                *
c *                                                                    *
c **********************************************************************
      implicit none
      real*8  x(3),r(2,2),xl(3),rt11,rt12,rt21,rt22,x1,x2,x3
      logical inv
c ......................................................................
      if(inv) then
         rt11 = r(1,1)*x(1)+r(2,1)*x(3)
         rt12 = r(1,1)*x(3)+r(2,1)*x(2)
         rt21 = r(1,2)*x(1)+r(2,2)*x(3)
         rt22 = r(1,2)*x(3)+r(2,2)*x(2)
         x1 = rt11*r(1,1)+rt12*r(2,1)
         x2 = rt21*r(1,2)+rt22*r(2,2)
         x3 = rt11*r(1,2)+rt12*r(2,2)
      else
         rt11 = r(1,1)*x(1)+r(1,2)*x(3)
         rt12 = r(1,1)*x(3)+r(1,2)*x(2)
         rt21 = r(2,1)*x(1)+r(2,2)*x(3)
         rt22 = r(2,1)*x(3)+r(2,2)*x(2)
         x1 = rt11*r(1,1)+rt12*r(1,2)
         x2 = rt21*r(2,1)+rt22*r(2,2)
         x3 = rt11*r(2,1)+rt12*r(2,2)
      endif
      xl(1) = x1
      xl(2) = x2
      xl(3) = x3
c ......................................................................      
      return
      end
      subroutine irotatek(s,r) 
c **********************************************************************
c *                                                                    *
c *                                                   21/11/2010       *
c *   ROTACAO:                                                         *
c *   Rotaciona matriz de rigidez                                      *
c *                                                                    *
c **********************************************************************
      implicit none 
      real*8 s(8,8),r(2,2),sl(8,8) 
      integer i,j 
c ......................................................................
      do j = 1,4 
         do i =1,4 
            sl(2*j-1,2*i-1)=r(1,1)*s(2*j-1,2*i-1)*r(1,1)+ 
     .                  r(2,1)*s(2*j,2*i)*r(2,1) 
            sl(2*j-1,2*i)=r(1,1)*s(2*j-1,2*i-1)*r(1,2)+ 
     .                  r(2,1)*s(2*j,2*i)*r(2,2) 
            sl(2*j,2*i-1)=r(1,2)*s(2*j-1,2*i-1)*r(1,1)+ 
     .                  r(2,2)*s(2*j,2*i)*r(2,1) 
            sl(2*j,2*i)=r(1,2)*s(2*j-1,2*i-1)*r(1,2)+ 
     .                  r(2,2)*s(2*j,2*i)*r(2,2) 
      enddo 
      enddo 
      do i = 1,8 
         do j = 1,8 
            s(i,j) = sl(i,j) 
         enddo 
      enddo 
c ...................................................................... 
      return 
      end
      subroutine irotatek3d(s,r,nst) 
c **********************************************************************
c *                                                                    *
c *                                                   21/11/2010       *
c *   ROTACAO:                                                         *
c *   Rotaciona matriz de rigidez 3d                                   *
c *                                                                    *
c **********************************************************************
      implicit none 
      integer i,j,k,l,i1,i2,i3,j1,j2,j3,nst
      real*8 s(nst,nst),r(3,3),sl(nst,nst),aux(nst,nst),raux(nst,nst)
c ......................................................................
      l = nst/3
      do i = 1, nst
         do j = 1, nst
            sl(i,j) = 0.d0
            aux(i,j) = 0.d0
            raux(i,j) = 0.d0
         enddo
      enddo
c
c      do i = 1, l 
c      i1 = 3*(i-1)+1
c      i2 = 3*(i-1)+2
c      i3 = 3*(i-1)+3
c         do j =1, l 
c            j1 = 3*(j-1)+1
c            j2 = 3*(j-1)+2
c            j3 = 3*(j-1)+3
c            do k = 1, 18
c               
c               aux(i1,j1) = aux(i1,j1) + r(1,k)*s(k,j1)
c               aux(i2,j2) = aux(i2,j2) + r(2,k)*s(k,j2)
c               aux(i3,j3) = aux(i3,j3) + r(3,k)*s(k,j3) 
c            enddo  
c         enddo 
c      enddo 
c      do i = 1, l 
c      i1 = 3*(i-1)+1
c      i2 = 3*(i-1)+2
c      i3 = 3*(i-1)+3
c         do j =1, l 
c            j1 = 3*(j-1)+1
c            j2 = 3*(j-1)+2
c            j3 = 3*(j-1)+3
c            do k = 1, 3
c               sl(i1,j1) = sl(i1,j1) + aux(i1,k)*r(1,k) 
c               sl(i2,j2) = sl(i2,j2) + aux(i2,k)*r(2,k) 
c               sl(i3,j3) = sl(i3,j3) + aux(i3,k)*r(3,k) 
c            enddo  
c         enddo 
c      enddo 
      do i = 1, l
         i1 = 3*(i-1)+1
         i2 = 3*(i-1)+2
         i3 = 3*(i-1)+3
            raux(i1,i1) = r(1,1)
            raux(i1,i2) = r(1,2)
            raux(i1,i3) = r(1,3)
            raux(i2,i1) = r(2,1)
            raux(i2,i2) = r(2,2)
            raux(i2,i3) = r(2,3)
            raux(i3,i1) = r(3,1)
            raux(i3,i2) = r(3,2)
            raux(i3,i3) = r(3,3)            
      enddo
      do i = 1, nst
         do j = 1, nst
            do k  = 1, nst
               aux(i,j) = aux(i,j) + raux(k,i)*s(k,j)
            enddo
         enddo
      enddo
      do i = 1, nst
         do j = 1, nst
            do k  = 1, nst
               sl(i,j) = sl(i,j) + aux(i,k)*raux(k,j)
            enddo
         enddo
      enddo      
      do i = 1,nst 
         do j = 1,nst 
            s(i,j) = sl(i,j) 
         enddo 
      enddo 
c ...................................................................... 
      return 
      end
c **********************************************************************
c
c **********************************************************************
      subroutine deform2d(hx,hy,u,eps,nen)
c **********************************************************************
c *                                                                    *
c *                                                    25/03/03        *
c *   DEFORM2D: calcula deformacoes 2D                                 *
c *   --------                                                         *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ---------------------                                            *
c *                                                                    *
c *     hx(nen)   - derivadas de h em relacao a x                      *
c *     hy(nen)   - derivadas de h em relacao a y                      *
c *     u(nst)    - deslocamentos                                      *
c *     nen       - numero de nos do elemento                          *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   -------------------                                              *
c *                                                                    *
c *     eps(3) - deformacoes                                           *
c *                                                                    *
c *                                                                    *
c **********************************************************************
      implicit none
      integer nen,i,j1,j2
      real*8  hx(*),hy(*),u(*),eps(*)
      real*8  ux,uy,vx,vy
c ......................................................................
      ux = 0.d0
      uy = 0.d0
      vx = 0.d0
      vy = 0.d0
      do 100 i = 1, nen
c       j1 = (i-1)*2+1
        j1 = 2*i-1    
        j2 = j1 + 1
        ux = ux + hx(i)*u(j1)
        uy = uy + hy(i)*u(j1)
        vx = vx + hx(i)*u(j2)
        vy = vy + hy(i)*u(j2)
  100 continue
      eps(1) = ux
      eps(2) = vy
      eps(3) = uy + vx
      return
      end
      subroutine deform3d(hx,hy,hz,u,eps,nen)
c **********************************************************************
c *                                                                    *
c *                                                    25/03/03        *
c *   DEFORM3D: calcula deformacoes 3D                                 *
c *   --------                                                         *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ---------------------                                            *
c *                                                                    *
c *     hx(nen)   - derivadas de h em relacao a x                      *
c *     hy(nen)   - derivadas de h em relacao a y                      *
c *     hz(nen)   - derivadas de h em relacao a z                      *
c *     u(nst)    - deslocamentos                                      *
c *     nen       - numero de nos do elemento                          *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   -------------------                                              *
c *                                                                    *
c *     eps(6) - deformacoes                                           *
c *                                                                    *
c *                                                                    *
c **********************************************************************
      implicit none
      integer nen,i,j1,j2,j3
      real*8  hx(*),hy(*),hz(*),u(*),eps(*)
      real*8  ux,uy,uz,vx,vy,vz,wx,wy,wz
c ......................................................................
      ux = 0.d0
      uy = 0.d0
      uz = 0.d0
      vx = 0.d0
      vy = 0.d0
      vz = 0.d0
      wx = 0.d0
      wy = 0.d0
      wz = 0.d0
      do 100 i = 1, nen
c        j1 = (i-1)*3+1
         j1 = 3*i-2
         j2 = j1 + 1
         j3 = j2 + 1
         ux = ux + hx(i)*u(j1)
         uy = uy + hy(i)*u(j1)
         uz = uz + hz(i)*u(j1)
         vx = vx + hx(i)*u(j2)
         vy = vy + hy(i)*u(j2)
         vz = vz + hz(i)*u(j2)
         wx = wx + hx(i)*u(j3)
         wy = wy + hy(i)*u(j3)
         wz = wz + hz(i)*u(j3)
  100 continue
      eps(1) = ux
      eps(2) = vy
      eps(3) = wz
      eps(4) = uy + vx
      eps(5) = vz + wy
      eps(6) = uz + wx
      return
      end
      subroutine stress2d(d11,d12,d22,d33,eps,t)
c **********************************************************************
c *                                                                    *
c *                                                    25/03/03        *
c *   STRESS2D: calcula tensoes 2D                                     *
c *   --------                                                         *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ---------------------                                            *
c *                                                                    *
c *     d11,d12,d22,d33 - coeficientes da matriz constitutiva          *
c *     eps(3) - deformacoes                                           *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   -------------------                                              *
c *                                                                    *
c *     t(3) - tensoes                                                 *
c *                                                                    *
c *                                                                    *
c **********************************************************************
      implicit none
      real*8 d11,d12,d22,d33,eps(*),t(*)
c ......................................................................
      t(1) = d11*eps(1) + d12*eps(2)
      t(2) = d12*eps(1) + d22*eps(2)
      t(3) = d33*eps(3)
      return
      end
c **********************************************************************
      subroutine stress2d_m(d11,d12,d22,d33,ps,ept,eps,t)
c **********************************************************************
c * Data de criacao    : 25/04/2016                                    *
c * Data de modificaco : 00/00/0000                                    * 
c * ------------------------------------------------------------------ * 
c * STRESS2D: calcula tensoes 2D                                       *
c * ------------------------------------------------------------------ *      
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ *                                                                  *
c * d11,d12,d22,d33 - coeficientes da matriz constitutiva              *
c * ps     - coeficiente de Poisson                                    *
c * ept    - estado plano de deformacai                                *
c * eps(3) - deformacoes                                               *
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *                                                                    *
c * ------------------------------------------------------------------ *                                              *
c * t(1) - sigmaxx                                                     *
c * t(2) - sigmayy                                                     *
c * t(3) - sigmazz                                                     *
c * t(4) - sigmaxy                                                     *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ * 
c **********************************************************************
      implicit none
      real*8 d11,d12,d22,d33,eps(*),t(*),ps
      logical ept
c ......................................................................
      t(1) = d11*eps(1) + d12*eps(2)
      t(2) = d12*eps(1) + d22*eps(2)
      if(ept) then 
        t(3) = ps*(t(1)+t(2)) 
      else
        t(3) = 0.0d0
      endif 
      t(4) = d33*eps(3)
      return
      end
      subroutine stress3d(a,b,c,eps,t)
c **********************************************************************
c *                                                                    *
c *                                                    25/03/03        *
c *   STRESS3D: calcula tensoes 3D                                     *
c *   --------                                                         *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ---------------------                                            *
c *                                                                    *
c *     a,b,c  - coeficientes da matriz constitutiva                   *
c *     eps(6) - deformacoes                                           *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   -------------------                                              *
c *                                                                    *
c *     t(6) - tensoes                                                 *
c *                                                                    *
c *                                                                    *
c **********************************************************************
      implicit none
      real*8  a,b,c,eps(*),t(*)
c ......................................................................
      t(1) = (  eps(1) + b*eps(2) + b*eps(3))*a
      t(2) = (b*eps(1) +   eps(2) + b*eps(3))*a
      t(3) = (b*eps(1) + b*eps(2) +   eps(3))*a
      t(4) = eps(4)*c*a
      t(5) = eps(5)*c*a
      t(6) = eps(6)*c*a
      return
      end
c **********************************************************************
c
c **********************************************************************
      subroutine jtria3(x,hx,hy,hz,det,afl,nel,ndm)
c **********************************************************************
c * Data de criacao    : 24/04/2016                                    *
c * Data de modificaco : 00/00/0000                                    *
c * ------------------------------------------------------------------ *   
c * JTRIA3 : Jacobiano de triangulo de 3 nos.                          *
c * ------------------------------------------------------------------ *  
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ * 
c * x(ndm,nen)- coordenadas nodais do elemento                         *
c * nel       - numero do elemento                                     *
c * afl       - .true. calcula das derivadas                           *
c * ndm       - dimensao                                               *
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ * 
c * hx(3) - derivadas de h em relacao a x                              *
c * hy(3) - derivadas de h em relacao a y                              *
c * det   - determinante da matriz jacobiana                           *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ *
c **********************************************************************
      implicit none
      logical afl
	integer nel,ndm
      real*8  x(ndm,*),hx(*),hy(*),hz(*),xj(3,3),det,deti,ZERO
	parameter (ZERO = 1.d-14)
c ......................................................................

c ... Matriz Jacobiana:

      xj(1,1) = x(1,1)-x(1,3)
      xj(1,2) = x(2,1)-x(2,3)
      xj(2,1) = x(1,2)-x(1,3)
      xj(2,2) = x(2,2)-x(2,3)

c ... Determinante da matriz Jacobiana:  
      det  = xj(1,1)*xj(2,2)-xj(1,2)*xj(2,1)
c ......................................................................
      if (det .le. ZERO) then
        print*,'*** Subrotina ELMT__: determinante <= 0 ',nel
        stop
      endif
c ......................................................................
c
c ... Derivadas das funcoes de interpolacao:

      if (afl) then
        deti  = 1.0d0/det
c           
        hx(1) =  xj(2,2)*deti
        hx(2) = -xj(1,2)*deti
        hx(3) = (xj(1,2)-xj(2,2))*deti
c
        hy(1) = -xj(2,1)*deti
        hy(2) =  xj(1,1)*deti
        hy(3) = (xj(2,1)-xj(1,1))*deti
      endif
      return
      end
c **********************************************************************
      subroutine jtetra4(x,hx,hy,hz,det,afl,nel)
c **********************************************************************
c *                                                                    *
c *                                                    22/01/05        *
c *                                                                    *
c *   JTETRA4: Jacobiano de tetraedros de 4 nos.                       *
c *   -------                                                          *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ---------------------                                            *
c *                                                                    *
c *     x(ndm,nen)- coordenadas nodais do elemento                     *
c *     nel       - numero do elemento                                 *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   -------------------                                              *
c *                                                                    *
c *     hx(4) - derivadas de h em relacao a x                          *
c *     hy(4) - derivadas de h em relacao a y                          *
c *     hz(4) - derivadas de h em relacao a z                          *
c *     det   - determinante da matriz jacobiana                       *
c *                                                                    *
c *                                                                    *
c **********************************************************************
      implicit none
      logical afl
	integer nel
      real*8  x(3,*),hx(*),hy(*),hz(*),xj(3,3),det,deti,ZERO
	parameter (ZERO = 1.d-14)
c ......................................................................

c ... Matriz Jacobiana:

      xj(1,1) = x(1,1)-x(1,4)
      xj(1,2) = x(2,1)-x(2,4)
	xj(1,3) = x(3,1)-x(3,4)
      xj(2,1) = x(1,2)-x(1,4)
	xj(2,2) = x(2,2)-x(2,4)
      xj(2,3) = x(3,2)-x(3,4)
      xj(3,1) = x(1,3)-x(1,4)
	xj(3,2) = x(2,3)-x(2,4)
      xj(3,3) = x(3,3)-x(3,4)

c ... Determinante da matriz Jacobiana:  
      
      det  = xj(1,1)*xj(2,2)*xj(3,3) + xj(1,2)*xj(2,3)*xj(3,1) +
     .       xj(1,3)*xj(2,1)*xj(3,2) - xj(3,1)*xj(2,2)*xj(1,3) -
     .       xj(1,2)*xj(2,1)*xj(3,3) - xj(1,1)*xj(3,2)*xj(2,3)
c ......................................................................
      if (det .le. ZERO) then
        print*,'*** Subrotina ELMT__: determinante <= 0 ',nel
        stop
      endif
c ......................................................................

c ... Derivadas das funcoes de interpolacao:

      if (afl) then
        deti  = 1.0d0/det
c           
        hx(1) = (xj(2,2)*xj(3,3)-xj(2,3)*xj(3,2))*deti
        hx(2) = (xj(1,3)*xj(3,2)-xj(1,2)*xj(3,3))*deti
        hx(3) = (xj(1,2)*xj(2,3)-xj(1,3)*xj(2,2))*deti
c
        hy(1) = (xj(2,3)*xj(3,1)-xj(2,1)*xj(3,3))*deti
        hy(2) = (xj(1,1)*xj(3,3)-xj(1,3)*xj(3,1))*deti
        hy(3) = (xj(1,3)*xj(2,1)-xj(1,1)*xj(2,3))*deti
c
        hz(1) = (xj(2,1)*xj(3,2)-xj(2,2)*xj(3,1))*deti
        hz(2) = (xj(1,2)*xj(3,1)-xj(1,1)*xj(3,2))*deti
        hz(3) = (xj(1,1)*xj(2,2)-xj(1,2)*xj(2,1))*deti
c
        hx(4) = -hx(1)-hx(2)-hx(3)
        hy(4) = -hy(1)-hy(2)-hy(3)
        hz(4) = -hz(1)-hz(2)-hz(3)
      endif
      return
      end
      subroutine invjtetra4(x,xji,det,nel)
c **********************************************************************
c *                                                                    *
c *                                                    22/01/05        *
c *                                                                    *
c *   INVJTETRA4: inversa do Jacobiano de tetraedros de 4 nos.         *
c *   ----------                                                       *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ---------------------                                            *
c *                                                                    *
c *     x(ndm,nen)- coordenadas nodais do elemento                     *
c *     nel       - numero do elemento                                 *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   -------------------                                              *
c *                                                                    *
c *     det      - determinante da matriz jacobiana                    *
c *     xji(3,3) - inversa da matriz jacobiana                         *
c *                                                                    *
c **********************************************************************
      implicit none
	integer nel
      real*8  xj(3,3),xji(3,3),x(3,*),det
	real*8  ZERO
	parameter (ZERO = 1.d-14)
c ......................................................................
c
c ... Matriz Jacobiana:

      xj(1,1) = x(1,1)-x(1,4)
      xj(1,2) = x(2,1)-x(2,4)
	xj(1,3) = x(3,1)-x(3,4)
      xj(2,1) = x(1,2)-x(1,4)
	xj(2,2) = x(2,2)-x(2,4)
      xj(2,3) = x(3,2)-x(3,4)
      xj(3,1) = x(1,3)-x(1,4)
 	xj(3,2) = x(2,3)-x(2,4)
      xj(3,3) = x(3,3)-x(3,4)
c
c ... Determinante da matriz Jacobiana:  
      
      det  = xj(1,1)*xj(2,2)*xj(3,3) + xj(1,2)*xj(2,3)*xj(3,1) +
     .       xj(1,3)*xj(2,1)*xj(3,2) - xj(3,1)*xj(2,2)*xj(1,3) -
     .       xj(1,2)*xj(2,1)*xj(3,3) - xj(1,1)*xj(3,2)*xj(2,3)
c ......................................................................
      if (det .le. ZERO) then
        print*,'*** Subrotina ELMT__: determinante <= 0 ',nel
        stop
      endif
c ......................................................................
	
c ..  Inversa da matriz Jacobiana:

      xji(1,1) = (xj(2,2)*xj(3,3)-xj(2,3)*xj(3,2))/det
      xji(1,2) = (xj(1,3)*xj(3,2)-xj(1,2)*xj(3,3))/det
      xji(1,3) = (xj(1,2)*xj(2,3)-xj(1,3)*xj(2,2))/det
      xji(2,1) = (xj(2,3)*xj(3,1)-xj(2,1)*xj(3,3))/det
      xji(2,2) = (xj(1,1)*xj(3,3)-xj(1,3)*xj(3,1))/det
      xji(2,3) = (xj(1,3)*xj(2,1)-xj(1,1)*xj(2,3))/det
      xji(3,1) = (xj(2,1)*xj(3,2)-xj(2,2)*xj(3,1))/det
      xji(3,2) = (xj(1,2)*xj(3,1)-xj(1,1)*xj(3,2))/det
      xji(3,3) = (xj(1,1)*xj(2,2)-xj(1,2)*xj(2,1))/det
c ......................................................................
      return
      end   
c **********************************************************************
c
c **********************************************************************
c * Data de criacao    : 02/12/2015                                    *
c * Data de modificaco :                                               * 
c * ------------------------------------------------------------------ * 
c * DARCY_FLUX_3D : fluxo de darcy e 3D                                *
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ *
c * perm   - coeficiente de permebilidade ( perm/peso especifico)      *
c * gl     - aceleracao da gravidade                                   *
c * fluid_d- massa especifica do fluido                                *
c * hx     - derivada das funcoes de interpolacao                      *
c * hy     - derivada das funcoes de interpolacao                      *
c * hz     - derivada das funcoes de interpolacao                      *
c * u      - pressao nodal                                             *
c * nen    - numero de pontos por elemento                             *
c * p      - nao definido                                              *
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ *                                                                     *
c * p(3) - fluxo de darcy( k( grad(P) + ro_fluid*g)                    *
c * ------------------------------------------------------------------ *       
c *   OBS:                                                             *
c * ------------------------------------------------------------------ * 
c **********************************************************************
      subroutine darcy_flux(perm,gl,fluid_d,hx,hy,hz,u,nen,p) 
      implicit none 
      real*8 perm,fluid_d,gl(3)
      real*8 hx(*),hy(*),hz(*),p(3),u(*) 
      integer i,nen
c ... ro_fluid*g 
      p(1) = fluid_d*gl(1)
      p(2) = fluid_d*gl(2)
      p(3) = fluid_d*gl(3)
c ......................................................................
c
c ... ro_fluid*g -grad(P)
      do i = 1, nen
        p(1) = p(1) - hx(i)*u(i)
        p(2) = p(2) - hy(i)*u(i)
        p(3) = p(3) - hz(i)*u(i)
      enddo
c .......................................................................
c
c ... 
      p(1) = perm*p(1)
      p(2) = perm*p(2)
      p(3) = perm*p(3)
c .......................................................................
      return
      end
c **********************************************************************
c
c **********************************************************************
c * Data de criacao    : 12/12/2015                                    *
c * Data de modificaco : 02/04/2016                                    * 
c * ------------------------------------------------------------------ * 
c * HEXA_VOL : calcula o volume do hexaedro                            *
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ * 
c * x - coordendas dos vertices do hexaedro                            *
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ * 
c * retorna o volume do hexaedro                                       *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ * 
c *                         | v1x v1y v1z |                            *
c * vol = v1 * ( v2 x v3) = | v2x v2y v3z |                            *
c *                         | v1x v1y v1z |                            *       
c **********************************************************************
      real*8 function hexa_vol(x)
      implicit none
      real*8 x(3,*),d(3,3),det
c ... v1 - no1 - no5
      d(1,1) = x(1,1)-x(1,5)
      d(1,2) = x(2,1)-x(2,5)
	d(1,3) = x(3,1)-x(3,5)
c ... v2 - no6 - no5     
      d(2,1) = x(1,6)-x(1,5)
	d(2,2) = x(2,6)-x(2,5)
      d(2,3) = x(3,6)-x(3,5)
c ... v3 - no8 - no5        
      d(3,1) = x(1,8)-x(1,5)
	d(3,2) = x(2,8)-x(2,5)
      d(3,3) = x(3,8)-x(3,5)
c
c ... Determinante:  
c      
      det  = d(1,1)*d(2,2)*d(3,3) + d(1,2)*d(2,3)*d(3,1) +
     .       d(1,3)*d(2,1)*d(3,2) - d(3,1)*d(2,2)*d(1,3) -
     .       d(1,2)*d(2,1)*d(3,3) - d(1,1)*d(3,2)*d(2,3)
c ...
      hexa_vol = det
c ......................................................................      
      return
      end
c **********************************************************************
c
c **********************************************************************
c * Data de criacao    : 27/03/2016                                    *
c * Data de modificaco :                                               * 
c * ------------------------------------------------------------------ * 
c * TETRA_VOL : calcula o volume do tetreadro                          *
c* ------------------------------------------------------------------- * 
c * x - coordendas dos vertices do hexaedro                            *
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ * 
c * retorna o volume do tetraedro                                      *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ * 
c *             | x1 - x4   y1 - y4   z1 - z4 |                        *
c * vol = (1/6) | x2 - x4   y2 - y4   z2 - z4 |                        *
c *             | x3 - x4   z3 - y4   z3 - z4 |                        *     
c **********************************************************************
      real*8 function tetra_vol(x)
      implicit none
      real*8  div6         
      parameter ( div6 = 0.166666666666667d0)
      real*8 x(3,*),d(3,3),det
c ... 
      d(1,1) = x(1,1)-x(1,4)
      d(1,2) = x(2,1)-x(2,4)
	d(1,3) = x(3,1)-x(3,4)
      d(2,1) = x(1,2)-x(1,4)
	d(2,2) = x(2,2)-x(2,4)
      d(2,3) = x(3,2)-x(3,4)
      d(3,1) = x(1,3)-x(1,4)
	d(3,2) = x(2,3)-x(2,4)
      d(3,3) = x(3,3)-x(3,4)
c
c ... Determinante da matriz Jacobiana:  
c      
      det  = d(1,1)*d(2,2)*d(3,3) + d(1,2)*d(2,3)*d(3,1) +
     .       d(1,3)*d(2,1)*d(3,2) - d(3,1)*d(2,2)*d(1,3) -
     .       d(1,2)*d(2,1)*d(3,3) - d(1,1)*d(3,2)*d(2,3)
c ...
      tetra_vol = det*div6
c ......................................................................      
      return
      end
c **********************************************************************
      block data
c **********************************************************************
c *                                                                    *
c *                                                    25/01/89        *
c *   - Descricao:                                                     *
c *                                                                    *
c *     Inicializa os pontos de integracao.                            *
c *                                                                    *
c *     Gauss-Legendre:                                                *
c *     --------------                                                 *
c *     do lz = 1, nint                                                *
c *       ti = pg(lz,nint)                                             *
c *       do ly = 1, nint                                              *
c *          si = pg(ly,nint)                                          *
c *          do lx = 1, nint                                           *
c *             ri = pg(lx,nint)                                       *
c *             wt = wg(lx,nint)*wg(ly,nint)*wg(lz,nint)*det           *
c *          enddo                                                     *
c *       enddo                                                        *
c *     enddo                                                          *
c *                                                                    *
c *     Integracao em triangulos:                                      *
c *     -------------------------                                      *
c *     nint = npint(igrau)                                            *
c *     do lx = 1, nint                                                *
c *        ri = pri(lx,igrau)                                          *
c *        si = psi(lx,igrau)                                          *
c *        ti = 1-ri-si                                                *
c *        call dgh(node,xl,ri,si,ti,det,noel,hfl)                     *
c *        wt = wf(lx,igrau) * 0.5d0 * det * thic                      *
c *     end do                                                         *
c *                                                                    *
c *     Integracao em tetraedros:                                      *
c *     -------------------------                                      *
c *     nint = npint4(grau)                                            *
c *     do lx = 1, nint                                                *
c *        ri = pri4(lx,igrau)                                         *
c *        si = psi4(lx,igrau)                                         *
c *        ti = pti4(lx,igrau)                                         *
c *        ui = 1 - ri - si - ti                                       *
c *        call dgh(node,xl,ri,si,ti,ui,det,noel,hfl)                  *
c *        wt = wf4(lx,igrau)                                          *
c *     end do                                                         *
c *                                                                    *
c **********************************************************************
      implicit none
      common /gauss/ pg, wg
      common /pint / pri,psi,wf,npint
      common /pint4/ pri4,psi4,pti4,wf4,npint4      
      real*8  pg(10,10), wg(10,10)
      real*8  pri(12,5),psi(12,5),wf(12,5)
      real*8  pri4(5,3),psi4(5,3),pti4(5,3),wf4(5,3)
      integer npint(5),npint4(3)
c ======================================================================
c
c ... Pontos de integracao de Gauss-Legendre:
c
      data pg /      0.d0,    0.d0, 0.d0,  0.d0,    0.d0, 0.d0,  0.d0,
     .               0.d0,    0.d0, 0.d0,
     .-.577350269189626d0,   .577350269189626d0,    0.d0, 0.d0,  0.d0,
     .  0.d0, 0.d0,  0.d0,    0.d0, 0.d0,
     .-.774596669241483d0,    0.d0,                .774596669241483d0,
     .  0.d0, 0.d0,  0.d0,    0.d0, 0.d0,  0.d0,    0.d0,
     .-.861136311594053d0,  -.339981043584856d0,   .339981043584856d0,
     . .861136311594053d0,    0.d0, 0.d0,  0.d0,    0.d0, 0.d0,  0.d0,
     .-.906179845938664d0,  -.538469310105683d0,    0.d0,
     . .538469310105683d0,   .906179845938664d0,    0.d0, 0.d0,  0.d0,
     .  0.d0, 0.d0,
     .-.932469514203152d0,  -.661209386466265d0,  -.238619186083197d0,
     . .238619186083197d0,   .661209386466265d0,   .932469514203152d0,
     .  0.d0, 0.d0,  0.d0,    0.d0,
     .-.949107912342759d0,  -.741531185599394d0,  -.405845151377397d0,
     .  0.d0,                .405845151377397d0,   .741531185599394d0,
     . .949107912342759d0,    0.d0, 0.d0,  0.d0,
     .-.960289856497536d0,  -.796666477413627d0,  -.525532409916329d0,
     .-.183434642495650d0,   .183434642495650d0,   .525532409916329d0,
     . .796666477413627d0,   .960289856497536d0,    0.d0, 0.d0,
     .-.968160239507626d0,  -.836031107326636d0,  -.613371432700590d0,
     .-.324253423403809d0,    0.d0,                .324253423403809d0,
     . .613371432700590d0,   .836031107326636d0,   .968160239507626d0,
     .  0.d0,
     .-.973906528517172d0,  -.865063366688985d0,  -.679409568299024d0,
     .-.433395394129247d0,  -.148874338981631d0,   .148874338981631d0,
     . .433395394129247d0,   .679409568299024d0,   .865063366688985d0,
     . .973906528517172d0 /
c
c ... Fatores de ponderacao de Gauss-Legendre:
c
      data wg /      2.d0,    0.d0, 0.d0,  0.d0,    0.d0, 0.d0,  0.d0,
     .               0.d0,    0.d0, 0.d0,
     .  1.d0, 1.d0,  0.d0,    0.d0, 0.d0,  0.d0,    0.d0, 0.d0,  0.d0,
     .  0.d0,
     . .555555555555556d0,   .888888888888889d0,   .555555555555556d0,
     .  0.d0, 0.d0,  0.d0,    0.d0, 0.d0,  0.d0,    0.d0,
     . .347854845137454d0,   .652145154862546d0,   .652145154862546d0,
     . .347854845137454d0,    0.d0, 0.d0,  0.d0,    0.d0, 0.d0,  0.d0,
     . .236926885056189d0,   .478628670499366d0,   .568888888888889d0,
     . .478628670499366d0,   .236926885056189d0,    0.d0, 0.d0,  0.d0,
     .  0.d0, 0.d0,
     . .171324492379170d0,   .360761573048139d0,   .467913934572691d0,
     . .467913934572691d0,   .360761573048139d0,   .171324492379170d0,
     .  0.d0, 0.d0,  0.d0,    0.d0,
     . .129484966168870d0,   .279705391489277d0,   .381830050505119d0,
     . .417959183673469d0,   .381830050505119d0,   .279705391489277d0,
     . .129484966168870d0,    0.d0, 0.d0,  0.d0,
     . .101228532690376d0,   .222381034453374d0,   .313706645877887d0,
     . .362683783378362d0,   .362683783378362d0,   .313706645877887d0,
     . .222381034453374d0,   .101228536290376d0,    0.d0, 0.d0,
     . .081274388361574d0,   .180648160694857d0,   .260610696402935d0,
     . .312347077040003d0,   .330239355001260d0,   .312347077040003d0,
     . .260610696402935d0,   .180648160694857d0,   .081274388361574d0,
     .  0.d0,
     . .066671344308688d0,   .149451349150581d0,   .219086362515982d0,
     . .269266719309996d0,   .295524224714753d0,   .295524224714753d0,
     . .269266719309996d0,   .219086362515982d0,   .149451349150581d0,
     . .066671344308688d0 /
c ======================================================================     
c
c ... Pontos de integracao em triangulos:
c ... Numero de pontos de integracao (grau = 1,..., 5):
      data npint /1,3,4,7,12/      
c ... Coordenadas r:
      data pri /.333333333333333d0,.0,.0,.0,.0,.0,.0,.0,.0,.0,.0,.0,
     .                      .5d0,.0,.5d0,.0,.0,.0,.0,.0,.0,.0,.0,.0,
     .    .333333333333333d0,.6d0,.2d0,.2d0,.0,.0,.0,.0,.0,.0,.0,.0,
     ..333333333333333d0,.059715871789770d0,.470142064105115d0,
     ..470142064105115d0,.797426985353087d0,.101286507323456d0,
     ..101286507323456d0,.0,.0,.0,.0,.0,
     ..873821971016996d0,.063089014491502d0,.063089014491502d0,
     ..501426509658179d0,.249286745170910d0,.249286745170910d0,
     ..636502499121399d0,.636502499121399d0,.310352451033785d0,
     ..310352451033785d0,.053145049844816d0,.053145049844816d0/
c ... Coordenadas s:
      data psi /.333333333333333d0,.0,.0,.0,.0,.0,.0,.0,.0,.0,.0,.0,
     .                      .5d0,.5d0,.0,.0,.0,.0,.0,.0,.0,.0,.0,.0,
     .    .333333333333333d0,.2d0,.6d0,.2d0,.0,.0,.0,.0,.0,.0,.0,.0,
     ..333333333333333d0,.470142064105115d0,.059715871789770d0,
     ..470142064105115d0,.101286507323456d0,.797426985353087d0,
     ..101286507323456d0,.0,.0,.0,.0,.0,
     ..063089014491502d0,.873821971016996d0,.063089014491502d0,
     ..249286745170910d0,.501426509658179d0,.249286745170910d0,
     ..310352451033785d0,.053145049844816d0,.636502499121399d0,
     ..053145049844816d0,.636502499121399d0,.310352451033785d0/
c ... Pesos para triangulos:
      data wf /                 1.d0,.0,.0,.0,.0,.0,.0,.0,.0,.0,.0,.0,
     . .333333333333333d0,.333333333333333d0,.333333333333333d0,.0,.0,
     ..0,.0,.0,.0,.0,.0,.0,
     .-.562500000000000d0,.520833333333333d0,.520833333333333d0,
     . .520833333333333d0,.0,.0,.0,.0,.0,.0,.0,.0,
     . .225030000300000d0,.132394152788506d0,.132394152788506d0,
     . .132394152788506d0,.125939180544827d0,.125939180544827d0,
     . .125939180544827d0,.0,.0,.0,.0,.0,
     . .050844906370207d0,.050844906370207d0,.050844906370207d0,
     . .116786275726379d0,.116786275726379d0,.116786275726379d0,
     . .082851075618374d0,.082851075618374d0,.082851075618374d0,
     . .082851075618374d0,.082851075618374d0,.082851075618374d0/
c=======================================================================     
c
c ... Pontos de integracao em tetraedros:
c
c ... Numero de pontos de integracao (grau = 1, 2, 3):
      data npint4 /1,4,5/    
c ... Coordenadas r:
      data pri4 / 0.25d0,0.d0,0.d0,0.d0,0.d0,
     .  0.58541020d0,0.13819660d0,0.13819660d0,0.13819660d0,0.d0,
     .  0.25d0,.5d0,.166666666666667d0,.166666666666667d0,
     .  .166666666666667d0/
c ... Coordenadas s:
      data psi4 / 0.25d0,0.d0,0.d0,0.d0,0.d0,
     .  0.13819660d0,0.58541020d0,0.13819660d0,0.13819660d0,0.d0,
     .  0.25d0,.166666666666667d0,.5d0,.166666666666667d0,
     .  .166666666666667d0/
c ... Coordenadas t:
      data pti4 / 0.25d0,0.d0,0.d0,0.d0,0.d0,
     . 0.13819660d0,0.13819660d0,0.58541020d0,0.13819660d0,0.d0,
     . 0.25d0,.166666666666667d0,.166666666666667d0,.5d0,
     ..166666666666667d0/
c ... Pesos para tetraedros:
      data wf4 /1.d0  , 0.d0  , 0.d0  , 0.d0  ,0.d0,
     .          0.25d0, 0.25d0, 0.25d0, 0.25d0,0.d0,
     .         -0.8d0 , 0.45d0, 0.45d0, 0.45d0,0.45d0/
c ......................................................................      
      end
