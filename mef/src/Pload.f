c **********************************************************************
      subroutine pload_mec(id,f,u,b,nload,nnode,ndf)
c **********************************************************************
c * Data de criacao    : 09/04/2016                                    *
c * Data de modificaco : 00/00/0000                                    * 
c * ------------------------------------------------------------------ *   
c * PLOAD: Monta o vetor de forcas do problema mecanico estatico       *
c * ------------------------------------------------------------------ *   
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ *   
c * id(ndf,nnode)    - numeracao das equacoes                          *
c * f(ndf,nnode)     - cargas concentradas e incognitas prescritas     *
c * u(ndf,nnode)     - graus de liberdade                              *
c * b(neq)           - nao definido                                    *
c * nload(ndf,nnode) - numero da carga nodal                           *
c * nnode            - numero de nos acessados na particao             *
c * ndf              - numero de graus de liberdade por no             *
c * ------------------------------------------------------------------ *      
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ *      
c * b - vetor de forcas                                                *
c **********************************************************************
      implicit none
      include 'transiente.fi'
      integer id(ndf,*),nload(ndf,*),nnode,ndf,i,j,k,nc
      real*8  f(ndf,*),u(ndf,*),b(*),c,vc(3),a
c.......................................................................
c
c ... Cargas nodais e desloc. prescritos variaveis no tempo:
c
      do 110 i = 1, nnode
        do 100 j = 1, ndf
          nc = nload(j,i)
          if(nc .gt. 0) then
            call tload(nc,t,u(j,i),c,vc)
            f(j,i) = c
          endif
  100   continue
  110 continue
c.......................................................................  
c
c ... Cargas nodais e deslocamentos prescritos no tempo t:
c
      a = alfa*dt
      do 200 i = 1, nnode
        do 200 j = 1, ndf
         k = id(j,i)
         if (k .gt. 0) then
            b(k) = f(j,i)
         else
            u(j,i) =  f(j,i)
         endif
  200 continue
c.......................................................................  
      return
      end
c **********************************************************************
c                                                                       
c **********************************************************************
      subroutine pload(x,id,f,u,v,b,nload,nnode,ndf)
c **********************************************************************
c *                                                                    *
c *   PLOAD: Monta o vetor de forcas                                   *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *                                                                    *
c *    id(ndf,nnode)   - numeracao das equacoes                        *
c *    f(ndf,nnode)    - cargas concentradas e incognitas prescritas   *
c *    b(neq) - nao definido                                           *
c *    nload(ndf,nnode) - numero da carga nodal                        *
c *    nnode  - numero de nos acessados na particao                    *
c *    ndf    - numero de graus de liberdade por no                    *
c *                                                                    *
c *   Parametros de saida:                                             *
c *                                                                    *
c *    b - vetor de forcas                                             *
c *                                                                    *
c **********************************************************************
      implicit none
      include 'transiente.fi'
      integer id(ndf,*),nload(ndf,*),nnode,ndf,i,j,k,nc
      real*8  x(3,*),f(ndf,*),u(ndf,*),v(ndf,*),b(*),c,vc(3),a
c.......................................................................
c
c ... Cargas nodais e desloc. prescritos variaveis no tempo:
c
      do 110 i = 1, nnode
        do 100 j = 1, ndf
          nc = nload(j,i)
          if(nc .gt. 0) then
            call tload(nc,t,u(j,i),c,vc)
            f(j,i) = c
          endif
  100   continue
  110 continue
c.......................................................................  
c
c ... Cargas nodais e deslocamentos prescritos no tempo t:
c
      a = alfa*dt
      do 200 i = 1, nnode
        do 200 j = 1, ndf
         k = id(j,i)
         if (k .gt. 0) then
            b(k) = f(j,i)
         else
c            v(j,i) = (f(j,i)-u(j,i))/a - v(j,i)*((1.d0-alfa)/alfa)
c            v(j,i) = 0.d0            
            v(j,i) = (f(j,i)-u(j,i))/dt
            u(j,i) =  f(j,i)
         endif
  200 continue
c.......................................................................  
      return
      end
c **********************************************************************
c                                                                        
c **********************************************************************
      subroutine tload(nc,t,u,c,vc)
c **********************************************************************
c *                                                                    *
c *   TLOAD:                                                           *
c *   -----                                                            *
c *                                                                    *
c *   Funcoes de cargas nodais                                         *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ---------------------                                            *
c *                                                                    *
c *    load(2,numload)  - load(1,j) = tipo da carga j                  *
c *                     - load(2,j) = numero de parcelas da carga j    *
c *    fload(numpar,numsum,numload) - definicao das cargas             *
c *                 fload(i,j,k) = parametro i da parcela j da carga j *
c *    nterms  - numero maximo de parcelas                             *
c *    npar    - numero maximo de parametros                           *
c *    nc      - numero da carga                                       *
c *    t       - tempo                                                 *
c *    u       - solucao                                               *
c *    c       - nao definido                                          *
c *    vc      - nao definido                                          *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   -------------------                                              *
c *                                                                    *
c *    c - valor da carga nodal no tempo t                             *
c *    vc- valor da carga vetorial nodal no tempo t                    *
c *                                                                    *
c *                                                                    *
c **********************************************************************
      implicit none
      include 'load.fi'
      integer nc,itype,nparc,i
      real*8  t,u,c,vc(*),a,w,t0,ti,DOISPI,tm,uext,b,h
      parameter (DOISPI = 6.283185307d0)
c.......................................................................
c      
c ... Cargas nodais variaveis no tempo:
c
      c  = 0.d0
      if(nc .le. 0) return
      itype = load(1,nc)
      nparc = load(2,nc)
c ...    Carga constante no tempo:
      if(itype .eq. 1) then
        c = fload(1,1,nc)
c ......................................................................
c
c ...    Variacao linear:            
      elseif(itype .eq. 2) then
        c = fload(1,1,nc) + fload(2,1,nc)*t
c ......................................................................
c
c ...    troca de calor de Newton:            
      elseif(itype .eq. 3) then
        c = (fload(1,1,nc)-u)*fload(2,1,nc)
c ......................................................................
c
c ...    Variacao senoidal:            
      elseif(itype .eq. 4) then
        do 200 i = 1, nparc
c ...       Amplitude:
          a = fload(1,i,nc)
c ...       Frequencia:
          w = 1.d0/fload(2,i,nc)
c ...       Fase:
          t0 = fload(3,i,nc)
          ti = max(t-t0,0.d0)
          c  = c + a*dsin(DOISPI*w*ti) 
  200   continue
c ......................................................................
c
c ...    Variacao cosssenoidal:            
      elseif(itype .eq. 5) then
        do 300 i = 1, nparc
c ...       Amplitude:
          a = fload(1,i,nc)
c ...       Frequencia:
          w = 1.d0/fload(2,i,nc)
c ...       Fase:
          t0 = fload(3,i,nc)
          ti = max(t-t0,0.d0)
          c  = c + a*dcos(DOISPI*w*ti)
  300   continue
c ......................................................................
c
c ...             
      elseif(itype .eq. 6) then
        do 400 i = 1, nparc  
          a  = fload(1,i,nc)
c           periodo:
          w  = 1.d0/fload(2,i,nc)
c           frequencia:
c		    w  = fload(2,i,nc)
          t0 = fload(3,i,nc)
c	      ti = max(t-t0,0.d0)
c           ti = t + t0
          c  = c + a*dcos(DOISPI*w*t+t0)
  400   continue
c ......................................................................
c
c ... kdu/dx = emiss * const(Stef-Boltz) *(uext4-u4)+H(uext-u)
      elseif(itype .eq. 7) then
c         call interpol(fload(1,2,nc),fload(1,3,nc),t,nparc,uext)
        call interpol(fload(1,2,nc),fload(1,3,nc),t,nparc,c)
c         a = fload(1,1,nc)
c         b = fload(2,1,nc)
c         h = fload(3,1,nc)
c         c = a*b*(uext**4-u**4)+h*(uext-u)
c ......................................................................
c
c ...  kdu/dx = emiss * const(Stef-Boltz) *(uext4-u4)+H(uext-u)
      elseif(itype .eq. 8) then
c         call interpol(fload(1,2,nc),fload(1,3,nc),t,nparc,uext)
c           if( t .ge. 172800.d0)   then
        call interpol(fload(1,1,nc),fload(1,2,nc),t,nparc,c)
c           else  
c              c = 0.d0
c           endif
c         a = fload(1,1,nc)
c         b = fload(2,1,nc)
c         h = fload(3,1,nc)
c         c = a*b*(uext**4-u**4)+h*(uext-u)
c ......................................................................
c
c ... forca distribuida constante no contorno
      elseif (itype .eq. 40) then
c ... numero de parcelas: 
        vc(1:nparc)  = fload(1:nparc,1,nc)
c .....................................................................
c
c ...
      else
        print*,' Subrotina TLOAD : carga invalida',itype
        call stop_mef()
      endif
c .....................................................................
      return
      end
c **********************************************************************
c
c **********************************************************************
      subroutine global(id,f,x,u,nnode,ndf)
c **********************************************************************
c *                                                                    *
c *   Subroutine GLOBAL                                                *
c *   -----------------                                                *
c *                                                                    *
c *                                                                    *
c **********************************************************************
      implicit none
      integer nnode,ndf,id(ndf,*),i,j,k
      real*8  f(ndf,*),x(*),u(ndf,*)
c ......................................................................
      do 100 i = 1, nnode
      do 100 j = 1, ndf
         k = id(j,i)
         if (k .gt. 0) then
            u(j,i) = x(k)
         else
            u(j,i) = f(j,i)
         endif
  100 continue
      return
      end
c **********************************************************************
c
c **********************************************************************
      subroutine updgeo(id,f,v,u,x,nnode,ndm,ndf)
c **********************************************************************
c *                                                                    *
c *   Subroutine UPDGEO                                                *
c *   -----------------                                                *
c *                                                                    *
c *                                                                    *
c **********************************************************************
      implicit none
      integer nnode,ndm,ndf,id(ndf,*),i,j,k
      real*8  x(ndm,*),u(ndf,*),v(*),f(ndf,*)
c ......................................................................
      do 100 i = 1, nnode
      do 100 j = 1, ndm
         k = id(j,i)
         if (k .gt. 0) then
            x(j,i) = x(j,i) + v(k)
         else
            x(j,i) = x(j,i) + f(j,i) - u(j,i)
         endif
  100 continue
      return
      end
c **********************************************************************
c
c **********************************************************************
      subroutine update_res(nnode,nnodev,ndf,u,u0,dp,pres0)
c **********************************************************************
c * Data de criacao    : 28/03/2016                                    *
c * Data de modificaco :                                               * 
c * ------------------------------------------------------------------ * 
c * UPDATE_RES : atualiza os resultados                                *
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ * 
c * nnode         - numero de nos acessados na particao                *
c * nnodev        - numero de nos de vertices acessados na particao    *
c * ndf           - numero de graus de liberdade por no                *
c * u0(ndf,nnode) - pressoes e deslocamentos em t                      *
c * u(ndf,nnode) - delta de pressoes e deslocamentos em t+dt           * 
c * dp(nnodev)    - nao definido                                       *
c * pres0(nnodev) - pressao inicial                                    *      
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ * 
c * u(ndf,nnode)   - pressoes e deslocamentos atualizados para o       * 
c *                  tempo t+dt                                        *
c * u0(ndf,nnode)  - pressoes e deslocamentos atualizados para o       * 
c *                  tempo t+dt                                        *
c * dp(nnodev)     - delta de pressao total                            * 
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ * 
c **********************************************************************
      implicit none
      integer nnode,nnodev,ndf
      integer i,j
      real*8  u(ndf,*),u0(ndf,*),dp(*),pres0(*)
c ......................................................................
      do 110 i = 1, nnode 
        do 100 j = 1, ndf-1 
           u(j,i)  = u0(j,i) + u(j,i)
           u0(j,i) = u(j,i)
  100   continue
  110 continue
c .......................................................................
c
c ...
      do 120 i = 1, nnodev
         u(ndf,i)  = u0(ndf,i) + u(ndf,i)
         dp(i)     = u(ndf,i)  - pres0(i) 
         u0(ndf,i) = u(ndf,i)
  120 continue
c .......................................................................
      return
      end
c **********************************************************************
c
c **********************************************************************
      subroutine update_res_v2(nnode,nnodev,ndf,u,u0,dp)
c **********************************************************************
c * Data de criacao    : 28/03/2016                                    *
c * Data de modificaco :                                               * 
c * ------------------------------------------------------------------ *     
c * UPDATE_RES_V2 : atualiza os resultados                             *
c * ------------------------------------------------------------------ *  
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ *  
c * nnode         - numero de nos acessados na particao                *
c * nnodev        - numero de nos de vertices acessados na particao    *
c * ndf           - numero de graus de liberdade por no                *
c * u0(ndf,nnode) - pressoes e deslocamentos em t                      *
c * u(ndf,nnode) - delta de pressoes e deslocamentos em t+dt           * 
c * dp(nnodev)    - nao definido                                       *
c * ------------------------------------------------------------------ *        
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ *          
c * u(ndf,nnode)   - pressoes e deslocamentos atualizados para o       * 
c *                  tempo t+dt                                        *
c * u0(ndf,nnode)  - pressoes e deslocamentos atualizados para o       * 
c *                  tempo t+dt                                        *
c * dp(nnodev)     - delta de pressao total                            * 
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ * 
c **********************************************************************
      implicit none
      integer nnode,nnodev,ndf
      integer i,j
      real*8  u(ndf,*),u0(ndf,*),dp(*)
c ......................................................................
      do 110 i = 1, nnode 
        do 100 j = 1, ndf-1 
           u(j,i)  = u0(j,i) + u(j,i)
           u0(j,i) = u(j,i)
  100   continue
  110 continue
c .......................................................................
c
c ...
      do 120 i = 1, nnodev
         dp(i)     = dp(i) + u(ndf,i) 
         u(ndf,i)  = u0(ndf,i) + u(ndf,i)
         u0(ndf,i) = u(ndf,i)
  120 continue
c .......................................................................
      return
      end
c **********************************************************************