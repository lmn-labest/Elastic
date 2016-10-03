c **********************************************************************
      subroutine pform_mec(ix      ,iq         ,ie       ,e
     .                    ,x       ,id         ,ia       ,ja
     .                    ,au      ,al         ,ad       ,b
     .                    ,u       ,tx0
     .                    ,xl      ,ul         ,pl   
     .                    ,sl      ,ld         ,txnl
     .                    ,numel   ,nen        ,nenv     ,ndf 
     .                    ,ndm     ,nst        ,neq      ,nad 
     .                    ,lhs     ,rhs        ,unsym 
     .                    ,stge    ,isw        ,ilib     ,nlit
     .                    ,i_colorg,i_elcolor  ,numcolors,fstress0)
c **********************************************************************
c * Data de criacao    : 09/04/2016                                    *
c * Data de modificaco : 24/04/2016                                    * 
c * ------------------------------------------------------------------ * 
c * PFORM_MEC:                                                         *
c * ------------------------------------------------------------------ *      
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ *                                                                  *
c * ix(nen+1,numel) - conetividades nodais                             *
c * iq(7,numel)     - cargas nos elementos                             *
c * ie(numat)       - tipo de elemento                                 *
c * e(10,numat)     - constantes fisicas dos materiais                 *
c * x(ndm,nnode)    - coordenadas nodais                               *
c * id(ndf,nnode)   - numeracao global das equacoes                    *
c * ia(neq+1) - ia(i) informa a posicao no vetor a do primeiro         *
c *             coeficiente nao-nulo da equacao i, no formato CSR      *
c * ja(nad)   - ja(k) informa a coluna do coeficiente que ocupa        *
c *             a posicao k no vetor a, no formato CSR (stge = 1)      *
c *           - ponteiro da diagonal do skyline (stge = 4)             *
c * ad(neq)   - nao definido                                           *
c * au(nad)   - nao definido                                           *
c * al(nad)   - nao definido                                           *
c * b(neq)    - vetor de forcas nodais equivalentes                    *
c * u(ndf,nnode) - solucao (com valores prescritos)                    *
c * xl(ndm,nen)  - nao definido                                        *
c * ul(nst)      - nao definido                                        *
c * pl(nst)      - nao definido                                        *
c * sl(nst,nst)  - nao definido                                        *
c * ld(nst)      - nao definido                                        *
c * numel - numero de elementos                                        *
c * nen   - numero de nos por elemento                                 *
c * nenv  - numero de nos de vertice por elemento                      *
c * ndf   - numero max. de graus de liberdade por no                   *
c * ndm   - dimensao                                                   *
c * nst   - nen*ndf                                                    *
c * neq   - numero de equacoes                                         *
c * lhs   - flag para montagem de ad (letf hand side)                  *
c * rhs   - flag para correcao de b  (right hand side)                 *
c * unsym - flag para matrizes nao simetricas                          *
c * stge  - = 1, armazenamento CSR                                     *
c *         = 2, armazenamento por arestas                             *
c *         = 3, armazenamento EBE                                     *
c *         = 4, armazenamento SKYLINE                                 *
c *         = 0, nao monta matriz                                      *
c * isw   - codigo de instrucao para as rotinas de elemento            *
c * ilib  - determina a biblioteca do elemento                         *
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ * 
c *     b(neq) - vetor de forcas corrigido     (rhs = true)            *
c *    ad(neq) - coeficientes da diagonal      (lhs = true)            *
c *    au(nad) - parte triangular sup. de A (lhs = true)               *
c *    al(nad) - parte triangular inf. de A (lhs = true, unsym = true) *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ * 
c **********************************************************************
      implicit none
      include 'openmp.fi'
      include 'omp_lib.h'
      include 'transiente.fi'
      include 'termprop.fi'
      integer numel,nen,nenv,ndf,ndm,nst,nad
      integer stge,isw,numat,nlit
      integer neq
      integer ix(nen+1,*),iq(7,*),ie(*),id(ndf,*),ld(nst)
      integer ia(*),ja(*)
      integer iel,ma,nel,no,i,j,k,kk,nad1,ilib
      integer i_colorg(2,*),i_elcolor(*),numcolors
      integer ic,jc
      real*8  e(prop,*),x(ndm,*),ad(*),au(*),al(*),b(*)
      real*8  u(ndf,*),tx0(6,*)
      real*8  xl(ndm,nenv),ul(nst),el(prop)
      real*8  pl(nst),sl(nst,nst),txnl(6,nen)
      logical lhs,rhs,unsym,fstress0
c .....................................................................
c
c ... Zera a matriz de coeficientes:
c
      if(lhs) then
        call azero(au,nad)
        call azero(al,nad)
        call azero(ad,neq)      
      endif
c ..................................................................... 
c
c ... openmp
      if(omp_elmt)then
c ... loop nos cores
c$omp parallel num_threads(nth_elmt)
c$omp.default(none)
c$omp.private(nel,ic,jc,no,k,ma,iel,i,j,kk,xl,ld,ul,pl,el,sl,txnl)
c$omp.shared(numcolors,i_colorg,i_elcolor,nen,nenv,ndf,ndm)
c$omp.shared(id,u,ix,tx0,ie,e)
c$omp.shared(stge,unsym,rhs,lhs,ilib,nlit,isw,nst,dt)
c$omp.shared(neq,nad,fstress0)
        do ic = 1, numcolors
c$omp do
c ... Loop nos elementos:
c     ------------------
          do jc = i_colorg(1,ic), i_colorg(2,ic)
            nel = i_elcolor(jc)
c           print*,omp_get_thread_num(),jc,nel,ic
            kk = 0
c ... loop nos nos de deslocamentos
            do i = 1, nen
              no = ix(i,nel)
              do j = 1, ndf
                kk     = kk + 1
                ld(kk) = id(j,no)
                ul(kk) = u(j,no)
                pl(kk) = 0.d0
              enddo
c ......................................................................              
           enddo   
c ......................................................................               
c
c ... loop nos nos vertices
            do i = 1, nenv
              no       = ix(i,nel)
              do j = 1, ndm
                xl(j,i) = x(j,no)
              enddo   
            enddo    
c ......................................................................
c            
c ... tensao inicial
          if(fstress0) then
            do i = 1, nen  
              no        = ix(i,nel)  
              txnl(1,i) = tx0(1,no)
              txnl(2,i) = tx0(2,no) 
              txnl(3,i) = tx0(3,no) 
              txnl(4,i) = tx0(4,no) 
              txnl(5,i) = tx0(5,no) 
              txnl(6,i) = tx0(6,no)
            enddo  
          endif  
c ......................................................................                 
c
c ...... Arranjos de elemento:
           ma  = ix(nen+1,nel)
           iel = ie(ma)
           do i = 1, prop
             el(i) = e(i,ma)
           enddo   
c ......................................................................
c
c ...... Chama biblioteca de elementos:
           call elmlib_mec(el,iq(1,nel),xl ,ul,pl ,sl ,txnl
     .                    ,dt,ndm     ,nst ,nel     ,iel,isw
     .                    ,ma,nlit    ,ilib)
c ......................................................................
c
c ...... Monta arrranjos locais em arranjos globais:
           call assbly(sl      ,pl   ,ld
     .                ,ia      ,ja   ,au
     .                ,al      ,ad   ,b    ,nst
     .                ,neq     ,lhs  ,rhs  ,unsym,stge)           
          enddo   
c$omp end do
c .....................................................................
        enddo
c$omp end parallel
c .....................................................................
c
c ... sequencial
      else  
c ... Loop nos elementos:
c     ------------------
        do nel = 1, numel
          kk = 0
c ... loop em todos os nos
          do i = 1, nen
            no = ix(i,nel)
            do j = 1, ndf
              kk     = kk + 1
              ld(kk) = id(j,no)
              pl(kk) = 0.d0
              ul(kk) = u(j,no)
            enddo
c ......................................................................
         enddo
c ......................................................................      
c
c ... loop nos nos vertices
          do i = 1, nenv
            no       = ix(i,nel)
            do j = 1, ndm
              xl(j,i) = x(j,no)
            enddo   
          enddo
c .....................................................................
c            
c ... tensao inicial
          if(fstress0) then
            do i = 1, nen  
              no        = ix(i,nel)  
              txnl(1,i) = tx0(1,no)
              txnl(2,i) = tx0(2,no) 
              txnl(3,i) = tx0(3,no) 
              txnl(4,i) = tx0(4,no) 
              txnl(5,i) = tx0(5,no) 
              txnl(6,i) = tx0(6,no)
            enddo  
          endif  
c ......................................................................      
c
c ...... Arranjos de elemento:
          ma  = ix(nen+1,nel)
          iel = ie(ma)
          do i = 1, prop
            el(i) = e(i,ma)
          enddo   
c ......................................................................
c
c ...... Chama biblioteca de elementos:
          call elmlib_mec(el,iq(1,nel),xl ,ul  ,pl ,sl ,txnl
     .                   ,dt,ndm     ,nst ,nel ,iel,isw
     .                   ,ma,nlit    ,ilib)
c ......................................................................
c
c ...... Monta arrranjos locais em arranjos globais:
          call assbly(sl      ,pl   ,ld
     .                ,ia      ,ja   ,au
     .                ,al      ,ad   ,b    ,nst
     .                ,neq     ,lhs  ,rhs  ,unsym,stge)  
        enddo 
c ......................................................................
      endif
c ......................................................................
c
c ...      
      return
      end
c **********************************************************************
c
c **********************************************************************
      subroutine tform_mec(ix    ,x    ,e   ,ie
     .                    ,ic    ,xl   ,ul  
     .                    ,pl    ,u    ,tx0 ,t   
     .                    ,nnodev,numel,nen ,nenv
     .                    ,ndm   ,ndf  ,nst ,ntn  
     .                    ,isw   ,ilib)
c **********************************************************************
c * Data de criacao    : 09/04/2016                                    *
c * Data de modificaco : 00/00/0000                                    * 
c * ------------------------------------------------------------------ * 
c * TFORM: caluclo das tensoes e dos fluxo nodal nos vertices          *                                                   *
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ * 
c * ix(nen+1,numel)  - conetividades nodais                            *
c * x(ndm,nnode)     - coordenadas nodais                              *
c * e(10,numat)      - constantes fisicas dos materiais                *
c * ie(numat)        - tipo de elemento                                *
c * ic(nnode)        - nao definido                                    *
c * xl(ndm,nen)      - nao definido                                    *
c * ul(ndf,nen)      - nao definido                                    *
c * pl(ntn*nen)      - nao definido                                    *
c * u(ndf,nnode)     - solucao corrente                                *
c * tx0(ntn,nnodev)  - tensao inicial                                  *
c * t(ntn,nnodev)    - nao definido                                    *
c * nnodev           - numero de nos de vertices                       *
c * numel            - numero de elementos                             *
c * nen              - numero max. de nos por elemento                 *
c * nenv             - numero de nos de vertice por elemento           *
c * ndf              - numero max. de graus de liberdade por no        *
c * nst              - nst = nen*ndf                                   *
c * ndm              - dimensao                                        *
c * ndf              - numero max. de graus de liberdade por no        *
c * ntn              - numero max. de derivadas por no                 *
c * ovlp             - chave indicativa de overlapping                 *
c * nprcs            - numero de processos                             *
c * isw              - codigo de instrucao para as rotinas             *
c *                    de elemento                                     *
c * ilib             - determina a biblioteca do elemento              *
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ * 
c * t(ntn   ,nnodev) - tensoes medias nodais                           *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ * 
c * Calcula as tensoes apenas nos vertices                             *     
c **********************************************************************
      use Malloc
      implicit none
c     include 'mpif.h'
      include 'parallel.fi'
      include 'termprop.fi'
      integer nnodev,numel,nen,nenv,ndf,nst,ndm,ntn
c ......................................................................      
      integer ix(nen+1,*),ie(*),ic(nnodev)
      integer nel,ma,iel,i,j,k,k1,no,kk
      integer ilib,isw
      real*8  xl(ndm,nenv),ul(nst),pl(nenv*ntn)
      real*8  x(ndm,*),e(prop,*)
      real*8  u(ndf,*),el(prop),tx0(ntn,*)
      real*8  t(ntn,*)
c ...
      logical ldum
      integer idum
      real*8 ddum
c ......................................................................
c
c ... deslocamentos do vetor local p
c ......................................................................
c
c ...
      do 30 i = 1, nnodev
        ic(i) = 0
c ... tensao
        do 10 j = 1, ntn
           t(j,i) = 0.d0
   10   continue 
c ......................................................................
c
   30 continue 
c ......................................................................
c
c ... Loop nos elementos:
      do 900 nel = 1, numel
        kk = 0
c ...
        do 300 i = 1, nenv*ntn
          pl(i) = 0.d0
  300   continue
c .......................................................................
c
c ... loop em todos os nos
        do 400 i = 1, nen
          no = ix(i,nel)
c ... loop nos deslocamentos
          do 410 j = 1, ndf
            kk     = kk + 1
            ul(kk) = u(j,no)
  410     continue
c ......................................................................
  400   continue
c ......................................................................
c
c ... loop nos nos vertices
        do 510 i = 1, nenv
          no     = ix(i,nel)
          kk     = kk + 1
          do 500 j = 1, ndm
            xl(j,i) = x(j,no)
  500     continue
  510   continue
c ......................................................................
c
c ...... form element array
        ma  = ix(nen+1,nel)
        iel = ie(ma)      
        do 610 i = 1, prop
          el(i) = e(i,ma)
  610   continue
c ......................................................................
c
c ...... Chama biblioteca de elementos:
        call elmlib_mec(el,idum,xl,ul,pl,ddum,ddum,ddum,ndm,nst,nel
     .                 ,iel,isw,ma,idum,ilib)
c ...... media do vetor global
        do 800 i = 1, nenv
           no = ix(i,nel)
           if (no .le. 0) go to 800
           ic(no) = ic(no) + 1  
           do 700 j = 1, ntn
c ... tensao total
              k       = (i-1)*ntn + j
              t(j,no) = t(j,no)  + pl(k)
  700      continue
  800   continue
  900 continue
c .......................................................................
c
c ... Comunica vetor de contagem de elementos por no'
c     if (novlp) call allgatheri(ic,i_xfi)
c .......................................................................
      do 1000 i = 1, nnodev
c ... tensao
        do 1010 j = 1, ntn
c ... tensao total            
         t(j,i)  = t(j,i)/ic(i)  + tx0(j,i)
 1010   continue 
c ......................................................................
 1000 continue
c ......................................................................
c
c ......................................................................
      return
      end      
c **********************************************************************
c
c **********************************************************************
c     subroutine pform(ix,iq,ie,e,x,id,ia,ja,au,al,ad,b,u,u0,ut,v,a,du,
c    .         hi,tx,txp,eps,w,xl,ul,vl,zl,wl,pl,sl,tl,ld,numel,nen,ndf,
c    .         ndm,npi,nst,neq,nad,lhs,rhs,unsym,stge,isw,ilib,
c    .         i_colorg,i_elcolor,numcolors,nlit,flaghidr)
c **********************************************************************
c *                                                                    *
c *   PFORM:                                                           *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *                                                                    *
c *    ix(nen+1,numel) - conetividades nodais                          *
c *    iq(7,numel)     - cargas nos elementos                          *
c *    ie(numat)       - tipo de elemento                              *
c *    e(10,numat)     - constantes fisicas dos materiais              *
c *    x(ndm,nnode)    - coordenadas nodais                            *
c *    id(ndf,nnode)   - numeracao global das equacoes                 *
c *    ia(neq+1) - ia(i) informa a posicao no vetor a do primeiro      *
c *                coeficiente nao-nulo da equacao i, no formato CSR   *
c *    ja(nad)   - ja(k) informa a coluna do coeficiente que ocupa     *
c *                a posicao k no vetor a, no formato CSR (stge = 1)   *
c *              - ponteiro da diagonal do skyline (stge = 4)          *
c *    ad(neq)   - nao definido                                        *
c *    au(nad)   - nao definido                                        *
c *    al(nad)   - nao definido                                        *
c *    b(neq)    - vetor de forcas nodais equivalentes                 *
c *    u(ndf,nnode) - solucao (com valores prescritos)                 *
c *    v(ndf,nnode) - derivada no tempo da solucao                     *
c *    a(ndf,nnode) - derivada segunda no tempo da solucao             *
c *    w(ndm,nnode) - campo de velocidades de conveccao                *
c *    xl(ndm,nen)  - nao definido                                     *
c *    ul(nst)      - nao definido                                     *
c *    vl(nst)      - nao definido                                     *
c *    zl(nst)      - nao definido                                     *
c *    wl(ndm,nen)  - nao definido                                     *
c *    pl(nst)      - nao definido                                     *
c *    sl(nst,nst)  - nao definido                                     *
c *    ld(nst)      - nao definido                                     *
c *    numel - numero de elementos                                     *
c *    nen   - numero de nos por elemento                              *
c *    ndf   - numero max. de graus de liberdade por no                *
c *    ndm   - dimensao                                                *
c *    nst   - nen*ndf                                                 *
c *    neq   - numero de equacoes                                      *
c *    nad   - numero de posicoes no CSR          (storage = 1)        *
c *    adfl  - flag para montagem de ad                                *
c *    bdfl  - flag para armazenamento da diagonal em blocos           *
c *    bfl   - flag para correcao de b                                 *
c *    unsym - flag para matrizes nao simetricas                       *
c *    stge  - = 1, armazenamento CSR                                  *
c *            = 2, armazenamento por arestas                          *
c *            = 3, armazenamento EBE                                  *
c *            = 4, armazenamento SKYLINE                              *
c *            = 0, nao monta matriz                                   *
c *    isw   - codigo de instrucao para as rotinas de elemento         *
c *    ilib  - determina a biblioteca do elemento                      *
c *                                                                    *
c *   Parametros de saida:                                             *
c *                                                                    *
c *     b(neq) - vetor de forcas corrigido     (rhs = true)            *
c *    ad(neq) - coeficientes da diagonal      (lhs = true)            *
c *    au(nad) - parte triangular sup. de A (lhs = true)               *
c *    al(nad) - parte triangular inf. de A (lhs = true, unsym = true) *
c *                                                                    *
c **********************************************************************
c     implicit none
c      include 'openmp.fi'
c     include 'omp_lib.h'
c     include 'transiente.fi'
c     include 'termprop.fi'
c     integer numel,nen,ndf,ndm,nst,neq,nad,stge,isw,numat,nlit,npi
c     integer ix(nen+1,*),iq(7,*),ie(*),id(ndf,*),ld(ndf,nen)
c     integer ia(*),ja(*)
c     integer iel,ma,nel,no,i,j,k,nad1,ilib
c     integer i_colorg(2,*),i_elcolor(*),numcolors
c     integer ic,jc
c     real*8  e(prop,*),x(ndm,*),ad(neq),au(*),al(*),b(neq),du(*)
c     real*8  u(ndf,*),v(ndf,*),a(ndm,*),w(ndm,*),ut(1,*),hi(3,*)
c     real*8  xl(ndm,nen),ul(ndf,nst),vl(ndf,nen),zl(ndf,nen),el(prop)
c     real*8  wl(ndm,nen),pl(ndf,nen),sl(nst,nst),tl(nen),utl(nen)
c     real*8  tx(npi,*),eps(npi,*),u0(ndf*nen,*),txp(npi,*)
c     logical lhs,rhs,unsym,flaghidr
c ----------------------------------------------------------------------
c ... Verifica se a matriz eh retangular e calc. o no. de coeficientes:
c     nad1 = 0
c     if(stge .eq. 1) nad1 = ja(nad+1)
c     if(nad1 .gt. 0) nad1 = ia(2*neq+2)-1
c ......................................................................      
c
c.... Zera a matriz de coeficientes:
c
c     if(lhs) then
c        call azero(au,nad)
c        call azero(al,nad+nad1)
c        call azero(ad,neq)      
c     endif
c ----------------------------------------------------------------------
c
c.... Loop nos elementos:
c     ------------------
cc$omp parallel 
cc$omp.private(nel,jc,no,k,ma,iel,i,j,xl,wl,ld,ul,vl,pl,sl,zl,tl,utl,el)
c     do 1000 ic = 1, numcolors
cc$omp do        
c       do 900 jc = i_colorg(1,ic), i_colorg(2,ic)
c          nel = i_elcolor(jc)      
c
c...... Arranjos locais:
c
c       do 600 i = 1, nen
c         no = ix(i,nel)
c         if (no .gt. 0) go to 300
c            do 100 j = 1, ndm
c              xl(j,i) = 0.d0
c              wl(j,i) = 0.d0
c 100        continue
c            do 200 j = 1, ndf
c              ld(j,i) = 0
c              ul(j,i) = 0.d0
c              vl(j,i) = 0.d0
c              zl(j,i) = 0.d0               
c              pl(j,i) = 0.d0
c              tl(i)   = 0.d0
c 200        continue
c            go to 600
c 300     continue
c         do 400 j = 1, ndm
c            xl(j,i) = x(j,no)
c            wl(j,i) = w(j,no)             
c 400     continue
c         do 500 j = 1, ndf
c            k = id(j,no)
c            ld(j,i) = k
c            pl(j,i) = 0.d0
c            ul(j,i) = u(j,no)
c            vl(j,i) = v(j,no)
c            zl(j,i) = a(j,no)
c 500     continue
c         tl(i) = du(no)
c         utl(i)= ut(1,no)
c 600   continue
c ......................................................................
c
c ...... Arranjos de elemento:
c
c       ma  = ix(nen+1,nel)
c       iel = ie(ma)
c       do 610 i = 1, prop
c          el(i) = e(i,ma)
c 610   continue
c ...... Chama biblioteca de elementos:
c       call elmlib(el,iq(1,nel),xl,ul,vl,zl,wl,pl,sl,tl,utl,hi(1,nel),
c    .              u0(1,nel),tx(1,nel),txp(1,nel),eps(1,nel),ndm,nst,
c    .              nel,iel,isw,ma,nlit,ilib,flaghidr,1.d0)
c ----------------------------------------------------------------------
c ...... Monta arrranjos locais em arranjos globais:
c       call assbly(sl,pl,ld,ia,ja,au,al,ad,b,nst,neq,lhs,rhs,unsym,
c    .              stge)
c 900 continue
cc$omp end do
c1000 continue
cc$omp end parallel
c ......................................................................
c     return
c     end
c     subroutine tform(ix,x,e,ie,ic,xl,ul,pl,tl,u,du,hi,tx,eps,t,nnode,
c    .             numel,nen,ndm,npi,ndf,nst,ntn,i_xfi,ilib,flaghidr)
c **********************************************************************
c *                                                                    *
c *   Subroutine TFORM                                                 *
c *   ----------------                                                 *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ----------------------                                           *
c *                                                                    *
c *    ix(nen+1,numel) - conetividades nodais                          *
c *    x(ndm,nnode)    - coordenadas nodais                            *
c *    e(10,numat)     - constantes fisicas dos materiais              *
c *    ie(numat)       - tipo de elemento                              *
c *    ic(nnode)   - nao definido                                      *
c *    xl(ndm,nen) - nao definido                                      *
c *    ul(ndf,nen) - nao definido                                      *
c *    pl(ntn*nen) - nao definido                                      *
c *    u(ndf,nnode)- solucao corrente                                  *
c *    stres(nte,numel) - tensoes por elemento                         *
c *    t(ntn,nnode)- nao definido                                      *
c *    nnode - numero de nos                                           *
c *    numel - numero de elementos                                     *
c *    nen   - numero max. de nos por elemento                         *
c *    ndf   - numero max. de graus de liberdade por no                *
c *    nst   - nst = nen*ndf                                           *
c *    ndm   - dimensao                                                *
c *    ntn   - numero max. de derivadas por no                         *
c *    ovlp  - chave indicativa de overlapping                         *
c *    nprcs - numero de processos                                     *
c *    ilib - determina a biblioteca do elemento                       *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   --------------------                                             *
c *                                                                    *
c *    t(ntn,nnode)- tensoes medias nodais                             *
c *                                                                    *
c **********************************************************************
c     use Malloc
c     implicit none
c     include 'mpif.h'
c     include 'parallel.fi'
c     include 'termprop.fi'
c     integer nnode,numel,nen,ndf,nst,ndm,ntn,npi
c ... ponteiro      
c     integer*8 i_xfi
c     logical flaghidr
c ......................................................................      
c     integer ix(nen+1,*),ie(*),ic(*),dum
c     integer nel,ma,iel,i,j,k,no
c     integer ilib
c     real*8  x(ndm,*),e(prop,*),xl(ndm,*),ul(ndf,*),pl(*),tl(*)
c     real*8  utl(nen),u(ndf,*),t(ntn,*),du(*),hi(3,*),tx(npi,*)
c     real*8  eps(npi,*)
c ......................................................................
c     do 20 i = 1, nnode
c        ic(i) = 0
c        do 10 j = 1, ntn
c           t(j,i) = 0.d0
c  10    continue
c  20 continue
c ......................................................................
c .... loop on elements:
c     do 900 nel = 1, numel
c ...... set up local arrays
c       do 600 i = 1, nen
c         no = ix(i,nel)
c         if (no .gt. 0) go to 300
c         do 100 j = 1, ndm
c           xl(j,i) = 0.d0
c 100     continue
c         do 200 j = 1, ndf
c           ul(j,i)  = 0.d0
c 200     continue
c         go to 600
c 300     continue
c         do 400 j = 1, ndm
c           xl(j,i) = x(j,no)
c 400     continue
c         do 500 j = 1, ndf
c           ul(j,i) = u(j,no)
c 500     continue
c         tl(i) = du(no)
c 600   continue
c ...... form element array
c       ma  = ix(nen+1,nel)
c        ma  = ix(21,nel)
c       iel = ie(ma)      
c ...... Biblioteca de elementos:
c       if (iel .eq. 10 .or. iel .eq. 11 .or. iel .eq. 12)  goto 900
c       call elmlib(e(1,ma),dum,xl,ul,pl,pl,pl,pl,pl,tl,utl,hi(1,nel),
c    .    pl,tx(1,nel),tx(1,nel),eps(1,nel),ndm,nst,nel,iel,3,ma,1,ilib,
c    .    flaghidr)
c ...... add to global array
c 650   do 800 i = 1, nen
c          no = ix(i,nel)
c          if (no .le. 0) go to 800
c          ic(no) = ic(no) + 1
c          do 700 j = 1, ntn
c             k = (i-1)*ntn + j
c             t(j,no) = t(j,no) + pl(k)
c 700      continue
c 800   continue
c 900 continue
c .......................................................................
c ... Comunica vetor de contagem de elementos por no'
c     if (novlp) call allgatheri(ic,i_xfi)
c .......................................................................
c     do 1000 i = 1, nnode
c     do 1000 j = 1, ntn
c        t(j,i) = t(j,i)/ic(i)
c1000 continue
c ......................................................................
c     return
c     end      
c     subroutine t_efetive(ix,x,e,ie,ic,xl,ul,pl,tl,u,du,hi,tx,eps,t
c    .         ,nnode,numel,nen,ndm,npi,ndf,nst,ntn,i_xfi,ilib,flaghidr,
c    .                                         plastic)
c **********************************************************************
c *                                                                    *
c *   Subroutine T_EFETIVE                                             *
c *   ----------------                                                 *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ----------------------                                           *
c *                                                                    *
c *    ix(nen+1,numel) - conetividades nodais                          *
c *    x(ndm,nnode)    - coordenadas nodais                            *
c *    e(10,numat)     - constantes fisicas dos materiais              *
c *    ie(numat)       - tipo de elemento                              *
c *    ic(nnode)   - nao definido                                      *
c *    xl(ndm,nen) - nao definido                                      *
c *    ul(ndf,nen) - nao definido                                      *
c *    pl(ntn*nen) - nao definido                                      *
c *    u(ndf,nnode)- solucao corrente                                  *
c *    stres(nte,numel) - tensoes por elemento                         *
c *    t(ntn,nnode)- nao definido                                      *
c *    nnode - numero de nos                                           *
c *    numel - numero de elementos                                     *
c *    nen   - numero max. de nos por elemento                         *
c *    ndf   - numero max. de graus de liberdade por no                *
c *    nst   - nst = nen*ndf                                           *
c *    ndm   - dimensao                                                *
c *    ntn   - numero max. de derivadas por no                         *
c *    ovlp  - chave indicativa de overlapping                         *
c *    nprcs - numero de processos                                     *
c *    ilib - determina a biblioteca do elemento                       *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   --------------------                                             *
c *                                                                    *
c *    t(ntn,nnode)- tensoes medias nodais                             *
c *                                                                    *
c **********************************************************************
c     use Malloc
c     implicit none
c     include 'mpif.h'
c     include 'parallel.fi'
c     include 'termprop.fi'
c     integer nnode,numel,nen,ndf,nst,ndm,ntn,npi
c ... ponteiro      
c     integer*8 i_xfi
c     logical flaghidr
c ......................................................................      
c     integer ix(nen+1,*),ie(*),ic(*),dum
c     integer nel,ma,iel,i,j,k,no
c     integer ilib
c     real*8 plastic(2,*)
c     real*8  x(ndm,*),e(prop,*),xl(ndm,*),ul(ndf,*),pl(*),tl(*)
c     real*8  utl(nen),u(ndf,*),du(*),hi(3,*),tx(npi,*)
c     real*8 t(ntn,*)
c     real*8  eps(npi,*)
c ......................................................................
c      do 20 i = 1, nnode
c         ic(i) = 0
c         do 10 j = 1, ntn
c            t(j,i) = 0.d0
c   10    continue
c   20 continue
c ......................................................................
c .... loop on elements:
c     do 900 nel = 1, numel
c ...... set up local arrays
c       do 600 i = 1, nen
c         no = ix(i,nel)
c         if (no .gt. 0) go to 300
c         do 100 j = 1, ndm
c           xl(j,i) = 0.d0
c 100     continue
c         do 200 j = 1, ndf
c           ul(j,i)  = 0.d0
c 200     continue
c         go to 600
c 300     continue
c         do 400 j = 1, ndm
c           xl(j,i) = x(j,no)
c 400     continue
c         do 500 j = 1, ndf
c           ul(j,i) = u(j,no)
c 500     continue
c         tl(i) = du(no)
c 600   continue
c ...... form element array
c       ma  = ix(nen+1,nel)
c        ma  = ix(21,nel)
c       iel = ie(ma)      
c ...... Biblioteca de elementos:
c        if (iel .eq. 11 .or. iel .eq. 12)  then
c         call elmlib(e(1,ma),dum,xl,ul,pl,pl,pl,pl,pl,tl,utl,hi(1,nel),
c    .    pl,tx(1,nel),tx(1,nel),eps(1,nel),ndm,nst,nel,iel,5,ma,1,ilib,
c    .    flaghidr,plastic(1,nel))
c        else
c           plastic(1,nel) = 0
c           plastic(2,nel) = 0
c        endif
c 900 continue
c .......................................................................
c     return
c     end 
c ***********************************************************************
c


