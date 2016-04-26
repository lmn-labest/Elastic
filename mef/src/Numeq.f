      subroutine numeqpmec1(id,num,idr,nnode,nnodev,ndf,neq,nequ,neqp)
c **********************************************************************
c *                                                                    *
c *   Subroutine NUMEQPMEC1                                            *
c *                                                                    *
c *   Numeracao das equacoes problema poro-mecanico                    *
c *                                                                    *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *                                                                    *
c *    ix(nen+1,numel) - conetividades nodais dos elementos            *
c *    id(ndf,nnode)   - condicoes nodais (0 = livre, 1 = restringido) *
c *    num(nnode)      - numeracao original dos nos                    *
c *                      num(i) e o numero original do no i            *
c *    idr(ndf,nnode)  - nao definido                                  *
c *    nnode - numero de nos                                           *
c *    nnodev- numero de nos dos vertives                              *
c *    ndf   - numero max. de graus de liberdade por no                *
c *    neq   - nao definido                                            *
c *    nequ  - nao definido                                            *
c *    neqp  - nao definido                                            *
c *                                                                    *
c *   Parametros de saida:                                             *
c *                                                                    *
c *    idr   - numeracao nodal das equacoes                            *
c *    neq   - numero de equacoes                                      *
c *    nequ  - numero de equacoes em kuu                               *
c *    neqp  - numero de equacoes em kpp                               *
c *                                                                    *
c *  OBS: numera primeiro os deslocamentos                             *
c **********************************************************************
      implicit none
      integer nnode,nnodev,ndf,neq,nequ,neqp
      integer id(ndf,*),idr(ndf,*)
      integer num(*),i,j,k,n
c ......................................................................
c
c.... Numeracao nodal das equacoes:
c
c ... numera os deslocamentos 
      neq = 0
      do 110 i = 1, nnode
         n = num(i)
         do 100 j = 1, ndf - 1
            k = id(j,n)
            if (k .eq. 0) then
               neq = neq + 1
               idr(j,n) = neq
            else
               idr(j,n) = 0
            endif
  100    continue
  110 continue
      nequ = neq
c ......................................................................
c
c ... numera a pressao
      do 120 i = 1, nnode
         n = num(i)
         if( n .le. nnodev) then 
           k = id(ndf,n)
           if (k .eq. 0) then
             neq = neq + 1
             idr(ndf,n) = neq
           else
             idr(ndf,n) = 0
           endif
        endif
  120 continue
c ......................................................................
c
      neqp = neq - nequ
      return
      end
c **********************************************************************
      subroutine numeqpmec2(id,num,idr,nnode,nnodev,ndf,neq)
c **********************************************************************
c *                                                                    *
c *   Subroutine NUMEQPMEC2                                            *
c *                                                                    *
c *   Numeracao das equacoes problema poro-mecanico                    *
c *                                                                    *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *                                                                    *
c *    ix(nen+1,numel) - conetividades nodais dos elementos            *
c *    id(ndf,nnode)   - condicoes nodais (0 = livre, 1 = restringido) *
c *    num(nnode)      - numeracao original dos nos                    *
c *                      num(i) e o numero original do no i            *
c *    idr(ndf,nnode)  - nao definido                                  *
c *    nnode - numero de nos                                           *
c *    nnodev- numero de nos dos vertives                              *
c *    numel - numero de elementos                                     *
c *    nen   - numero de nos por elemento                              *
c *    ndf   - numero max. de graus de liberdade por no                *
c *    nst   - numero max. de graus de liberdade por elemento          *
c *    neq   - nao definido                                            *
c *                                                                    *
c *   Parametros de saida:                                             *
c *                                                                    *
c *    idr   - numeracao nodal das equacoes                            *
c *    neq   - numero de equacoes                                      *
c *                                                                    *
c *  OBS: numera simultaniamente os deslocamentos e as pressoes        *
c **********************************************************************
      implicit none
      integer nnode,nnodev,ndf,neq
      integer id(ndf,*),idr(ndf,*)
      integer num(*),i,j,k,n,ndfn
c ......................................................................
c
c.... Numeracao nodal das equacoes:
c
      neq = 0
      do 110 i = 1, nnode
         n = num(i)
c ...
         if( n .gt. nnodev) then
           ndfn = ndf - 1 
         else
           ndfn = ndf 
         endif
c .....................................................................
c
c ...
         do 100 j = 1, ndfn
           k = id(j,n)
           if (k .eq. 0) then
             neq = neq + 1
             idr(j,n) = neq
           else
             idr(j,n) = 0
           endif
  100    continue
c .....................................................................
  110 continue
c     do i = 1, nnode
c       n = num(i)
c       print*,n,idr(1,n),idr(2,n),idr(3,n),idr(4,n)
c     enddo
c ......................................................................  
      return
      end
      subroutine numeq(id,num,idr,nnode,ndf,neq)
c **********************************************************************
c *                                                                    *
c *   Subroutine NUMEQ                                                 *
c *                                                                    *
c *   Numeracao das equacoes                                           *
c *                                                                    *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *                                                                    *
c *    ix(nen+1,numel) - conetividades nodais dos elementos            *
c *    id(ndf,nnode)   - condicoes nodais (0 = livre, 1 = restringido) *
c *    num(nnode)      - numeracao original dos nos                    *
c *                      num(i) e o numero original do no i            *
c *    idr(ndf,nnode)  - nao definido                                  *
c *    nnode - numero de nos                                           *
c *    numel - numero de elementos                                     *
c *    nen   - numero de nos por elemento                              *
c *    ndf   - numero max. de graus de liberdade por no                *
c *    nst   - numero max. de graus de liberdade por elemento          *
c *    neq   - nao definido                                            *
c *    elflag- TRUE -> numera as equacoes por elemento                 *
c *                                                                    *
c *   Parametros de saida:                                             *
c *                                                                    *
c *    idr   - numeracao nodal das equacoes                            *
c *    neq   - numero de equacoes                                      *
c *                                                                    *
c **********************************************************************
      implicit none
      integer nnode,ndf,neq
      integer id(ndf,nnode),idr(ndf,nnode)
      integer num(nnode),i,j,k,n
c ......................................................................
c
c.... Numeracao nodal das equacoes:
c
      neq = 0
      do 110 i = 1, nnode
         n = num(i)
         do 100 j = 1, ndf
            k = id(j,n)
            if (k .eq. 0) then
               neq = neq + 1
               idr(j,n) = neq
            else
               idr(j,n) = 0
            endif
  100    continue
  110 continue
c ......................................................................  
      return
      end
      subroutine neqel(ix,id,ld,numel,nen,ndf,nst)
c **********************************************************************
c *                                                                    *
c *   Subroutine NEQEL: numeracao das equacoes por elemento            *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *                                                                    *
c *    ix(nen1,numel) - conetividades nodais dos elementos/arestas     *
c *    id(ndf,nnode)  - numeracao nodal das equacoes                   *
c *    ld(nst,numel)  - nao definido                                   *
c *    numel - numero de elementos                                     *
c *    nen   - numero de nos por elemento                              *
c *    ndf   - numero max. de graus de liberdade por no                *
c *    nst   - numero max. de graus de liberdade por elemento          *
c *                                                                    *
c *   Parametros de saida:                                             *
c *                                                                    *
c *    ld    - numeracao das equacoes por elemento                     *
c *                                                                    *
c **********************************************************************
      implicit none
      integer numel,nen,ndf,nst
      integer ix(nen+1,numel),id(ndf,*),ld(nst,numel)
      integer i,j,k,l,n
c ......................................................................
c
c.... Numeracao das equacoes por elemento:
c
      do 220 i = 1, numel
         l = 0
         do 210 j = 1, nen
            n = ix(j,i)
            do 200 k = 1, ndf
               l = l + 1
               ld(l,i) = 0
               if(n .gt. 0) ld(l,i) = id(k,n)
  200       continue
  210    continue
  220 continue
c ......................................................................              
      return
      end      
      subroutine neqedg(edge,id,ld,nedge,ndf,nst)
c **********************************************************************
c *                                                                    *
c *   NEQEDG                                                           *
c *                                                                    *
c *   Numeracao das equacoes por arestas                               *
c *                                                                    *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *                                                                    *
c *    id(ndf,nnode)   - numeracao nodal das equacoes                  *
c *    edge(2,nedge)   - conetividades nodais das arestas              *
c *    ld(nst,nedge)   - nao definido                                  *
c *    nedge - numero de arestas                                       *
c *    ndf   - numero max. de graus de liberdade por no                *
c *    nst   - numero max. de graus de liberdade por aresta            *
c *                                                                    *
c *   Parametros de saida:                                             *
c *                                                                    *
c *    ld    - numeracao das equacoes por aresta                       *
c *                                                                    *
c **********************************************************************
      implicit none
      integer nedge,ndf,nst,i,j,k,l,n
      integer edge(2,*),id(ndf,*),ld(nst,*)
c ......................................................................
c
c.... Numeracao das equacoes por aresta:
c
      do 220 i = 1, nedge
         l = 0
         do 210 j = 1, 2
            n = edge(j,i)
            do 200 k = 1, ndf
               l = l + 1
               ld(l,i) = 0
               if(n .gt. 0) ld(l,i) = id(k,n)
  200            continue
  210    continue
  220 continue
c ......................................................................              
      return
      end
