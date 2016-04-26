c **********************************************************************
c *                                                                    *
c *   CSR.F                                               31/08/2005   *
c *                                                                    *
c *   Este arquivo contem subrotinas para montagem da estrutura de     *
c *   dados CSR:                                                       *
c *                                                                    *
c *   csrstruct                                                        *
c *   csria                                                            *
c *   csrja                                                            *
c *                                                                    *
c **********************************************************************
      subroutine csrstruct(id,ix,num,nnode,numel,nen,ndf,neq,i2,i3,nad,
     .                     nad1,lower,diag,upper,right,ija,ja)
c **********************************************************************
c *                                                                    *
c *   CSRSTRUCT: monta os arranjos ia e ja do formato CSR.             *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *                                                                    *
c *    id(ndf,nnode)   - numeracao global das equacoes                 *
c *    ix(nen+1,numel) - conetividades nodais                          *
c *    num(nnode)      - numeracao RCM                                 *
c *    nnode - numero de nos total da particao                         *
c *    numel - numero de elementos                                     *
c *    nen   - numero de nos por elemento                              *
c *    ndf   - numero max. de graus de liberdade por no                *
c *    neq   - numero de equacoes                                      *
c *    lower = .true.  -> inclui a parte triangular inferior no csr    *
c *    diag  = .true.  -> inclui a diagonal no csr                     *
c *    upper = .true.  -> inclui a parte triangular superior no csr    *
c *                                                                    *
c *   Parametros de saida:                                             *
c *                                                                    *
c *    i2    - ponteiro para o arranjo ia(neq+1)                       *
c *    i3    - ponteiro para o arranjo ja(nad)                         *
c *    nad   - numero de coeficientes nao nulos                        *
c *    nad1  - numero de coeficientes nao nulos da parte retangular    *                                              *
c *                                                                    *
c **********************************************************************
      use Malloc
      implicit none
      integer id(ndf,*),ix(nen+1,*),num(*)
      integer nnode,numel,nen,ndf,neq,nad,nad1
c ... ponteiros      
      integer*8 i0,i1,i2,i3
c .....................................................................      
      integer n
      logical lower,diag,upper,right
      character*8 ija,ja
c ......................................................................
c
c ... Grafo da malha:
c
c     ia(i0) => ip(nnode+1) - ip(i) indica a posicao em ips do primeiro
c                             no vizinho ao no i.
c     ia(i1) => ips(ipos) - contem as conectividades nodais de cada no
      call graph(ix,nnode,numel,nen,i0,i1)
c ......................................................................      
c
c ... Montagem do arranjo ia(neq+1):
c
c ... ia(i2)=>ia(neq+1) - ia(i) informa a posicao no vetor a do primeiro
c                               coeficiente nao-nulo da equacao i   
c
      n = neq + 1
      if (right) n = 2*n
      i2 = alloc_4(ija,1,n)
      call csria(id,num,ia(i0),ia(i1),ia(i2),nnode,ndf,neq,nad,nad1,
     .           lower,diag,upper,right)
c ......................................................................     
c
c ... Montagem do arranjo ja(nad):
c
c ... ia(i3)=>ja(nad) - ja(k) informa a coluna do coeficiente que ocupa
c                       a posicao k no vetor a  
c
      i3 = alloc_4(ja,1,nad+nad1+1)
      ia(i3+nad) = 0
      call csrja(id,num,ia(i0),ia(i1),ia(i3),nnode,ndf,neq,nad,lower,
     .           diag,upper,right)
      call sortgraph(ia(i2),ia(i3),neq)
      if (right) then
         call sortgraph(ia(i2+neq+1),ia(i3+nad),neq)      
      endif
c ......................................................................
      i1 = dealloc('iaux1   ')
      i0 = dealloc('iaux0   ')
      i2 = locate (ija)
      i3 = locate (ja)
      return
      end
      subroutine csria(id,num,ip,ips,ia,nnode,ndf,neq,nad,nad1,lower,
     .                 diag,upper,right)
c **********************************************************************
c *                                                                    *
c *   CSRIA: monta o arranjo ia do formato CSR                         *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *                                                                    *
c *    id(ndf,nnode)- numeracao global das equacoes                    *
c *    num(nnode)   - renumeracao nodal                                *
c *    ip(nnode+1)  - ip(i) indica a posicao em ips do primeiro no     *
c *                   conectado ao no i                                *
c *    ips(ipos)    - contem as conetividades nodais de cada no        *
c *    nnode - numero de nos                                           *
c *    ndf   - numero max. de graus de liberdade por no                *
c *    neq   - numero de equacoes                                      *
c *    lower = .true.  -> inclui a parte triangular inferior no csr    *
c *    diag  = .true.  -> inclui a diagonal no csr                     *
c *    upper = .true.  -> inclui a parte triangular superior no csr    *
c *                                                                    *
c *   Parametros de saida:                                             *
c *                                                                    *
c *    ia(neq+1) - ia(i) informa a posicao no vetor a do primeiro      *
c *                      coeficiente nao-nulo da equacao i             *
c *    nad   - numero de coeficientes nao nulos                        *
c *                                                                    *
c **********************************************************************
      implicit none
      integer nnode,ndf,neq,nad,nad1
      integer id(ndf,nnode),ip(nnode+1),ips(*),ia(*),num(nnode)
      integer i,j,k,ii,jj,kk,neqi,neqj,no
      logical lower,diag,upper,right
c ......................................................................
c
c ... Inicializa o arranjo ia:
c
      do 50 i = 1, neq
         ia(i) = 0
         if (right) ia(i+neq+1) = 0
   50 continue
c ----------------------------------------------
c
c ... Loop nos vertices:
c
      do 140 no = 1, nnode
         i = num(no)
c
c ...    Loop nas equacoes do vertice i:
c
         do 130 ii = 1, ndf
            neqi = id(ii,i)
            if (neqi .gt. 0 .and. neqi .le. neq) then
c ----------------------------------------------
c
c ...          Loop nas equacoes do vertice i:
c
               do 100 kk = 1, ndf
                  neqj = id(kk,i)
                  if (neqj .le. 0) goto 100
                  if (neqj .lt. neqi) then
                     if (lower) ia(neqi) = ia(neqi) + 1
                  elseif (neqj .eq. neqi) then
                     if (diag) ia(neqi)  = ia(neqi) + 1
                  elseif (neqj .gt. neqi) then
                     if (upper) ia(neqi) = ia(neqi) + 1
                  endif
  100          continue
c
c ...          Loop nos vertices conectados ao vertice i:
c
               do 120 k = ip(i), ip(i+1)-1
                  j = ips(k)
c
c ...             Loop nas equacoes do vertice j:
c
                  do 110 jj = 1, ndf
                     neqj = id(jj,j)
                     if (neqj .le. 0) goto 110
                     if (neqj .lt. neqi) then
                        if (lower) ia(neqi) = ia(neqi) + 1
                     elseif (neqj .gt. neq) then
                        if (right) ia(neq+1+neqi) = ia(neq+1+neqi) + 1
                     elseif (neqj .gt. neqi) then
                        if (upper) ia(neqi) = ia(neqi) + 1
                     endif
  110             continue
  120          continue
c ----------------------------------------------
            endif
  130    continue
  140 continue
c ----------------------------------------------
      do 200 i = neq, 1, -1
         ia(i+1) = ia(i)
         if (right) ia(i+neq+2) = ia(i+neq+1)         
  200 continue
      ia(1) = 1
      if (right) ia(neq+2) = 1
      do 210 i = 1, neq
         ia(i+1) = ia(i+1) + ia(i)
         if (right) ia(i+neq+2) = ia(i+neq+2) + ia(i+neq+1)         
  210 continue
      nad  = ia(neq+1)-1
      nad1 = 0
      if (right) nad1 = ia(2*neq+2)-1      
c ......................................................................      
      return
      end
      subroutine csrja(id,num,ip,ips,ja,nnode,ndf,neq,nad,lower,diag,
     .                 upper,right)
c **********************************************************************
c *                                                                    *
c *   CSRJA: monta o arranjo ja do formato CSR (matrizes simetricas    *
c *                                             e nao simetricas)      *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *                                                                    *
c *    id(ndf,nnode)- numeracao global das equacoes                    *
c *    num(nnode)   - renumeracao nodal                                *
c *    ip(nnode+1)  - ip(i) indica a posicao em ips do primeiro        *
c *    ips(ipos)    - contem as conetividades nodais de cada no        *
c *    nnode - numero de nos                                           *
c *    ndf   - numero max. de graus de liberdade por no                *
c *    neq   - numero de equacoes                                      *
c *    nad   - numero de coeficientes nao nulos                        *
c *    lower = .true.  -> inclui a parte triangular inferior no csr    *
c *    diag  = .true.  -> inclui a diagonal no csr                     *
c *    upper = .true.  -> inclui a parte triangular superior no csr    *
c *                                                                    *
c *   Parametros de saida:                                             *
c *                                                                    *
c *    ja(nad) - ja(k) informa a coluna do coeficiente que ocupa       *
c *              a posicao k no vetor a                                *
c *                                                                    *
c **********************************************************************
      implicit none
      integer nnode,ndf,neq,nad
      integer id(ndf,*),num(*),ip(*),ips(*),ja(*)
      integer i,j,k,ii,jj,kk,neqi,neqj,no,n,m
      logical lower,diag,upper,right
c ......................................................................
      m = 0
      n = 0
c
c ... Loop nos vertices:
c
      do 140 no = 1, nnode
         i = num(no)
c
c ...    Loop nas equacoes do vertice i:
c
         do 130 ii = 1, ndf
            neqi = id(ii,i)
            if (neqi .gt. 0 .and. neqi .le. neq) then
c ----------------------------------------------------------
c
c ...          Loop nas equacoes do vertice i:
c
               do 100 kk = 1, ndf
                  neqj  = id(kk,i)
                  if (neqj .le. 0) goto 100
                  if (neqj .lt. neqi) then
                     if (lower) then
                        n = n + 1
                        ja(n) = neqj
                     endif
                  elseif (neqj .eq. neqi) then
                     if (diag) then
                        n = n + 1
                        ja(n) = neqj
                     endif
                  elseif (neqj .gt. neqi) then
                     if (upper) then
                        n = n + 1
                        ja(n) = neqj
                     endif
                  endif                  
  100          continue
c
c ...          Loop nos vertices conectados ao vertice i:
c
               do 120 k = ip(i), ip(i+1)-1
                  j = ips(k)
c
c ...             Loop nas equacoes do vertice j:
c
                  do 110 jj = 1, ndf
                     neqj  = id(jj,j)
                     if (neqj .le. 0) goto 110
                     if (neqj .lt. neqi) then
                        if (lower) then
                           n = n + 1
                           ja(n) = neqj
                        endif
                     elseif (neqj .gt. neq) then
                        if (right) then
                           m = m + 1
                           ja(nad+m) = neqj
                        endif
                     elseif (neqj .gt. neqi) then
                        if (upper) then
                           n = n + 1
                           ja(n) = neqj
                        endif
                     endif
  110             continue
  120          continue
c ----------------------------------------------------------
            endif
  130    continue
  140 continue
c ......................................................................
      return
      end
c **********************************************************************