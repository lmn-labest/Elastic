c*****************************Svn***************************************      
c*$Date: 2011-03-16 15:32:53 -0300 (Wed, 16 Mar 2011) $                 
c*$Rev: 914 $                                                           
c*$Author: henrique $                                                   
c***********************************************************************
c **********************************************************************
c *                                                                    *
c *   GRAPH.F                                             31/08/2005   *
c *                                                                    *
c *   Este arquivo contem subrotinas para montagem do grafo da malha:  *
c *                                                                    *
c *   graph                                                            *
c *   aponta                                                           *
c *   graph0                                                           *
c *   graph1                                                           *
c *   sortgrpah                                                        *
c *   bubblesort                                                       *
c *                                                                    *
c **********************************************************************
      subroutine graph(ix,nnode,numel,nen,i0,i1)
c **********************************************************************
c *                                                                    *
c *   GRAPH: monta o grafo da malha.                                   *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *                                                                    *
c *    ix(nen+1,numel) - conetividades nodais                          *
c *    nnode - numero de nos                                           *
c *    numel - numero de elementos                                     *
c *    nen   - numero de nos por elemento                              *
c *                                                                    *
c *   Parametros de saida:                                             *
c *                                                                    *
c *    i0    - ponteiro para o arranjo ia(i0: nnode+1), apontador do   *
c *            do arranjo que contem o grafo                           *
c *    i1    - ponteiro para o arranjo ia(i1: ipos), que contem o      *
c *            o grafo                                                 *
c *                                                                    *
c *    Obs: esta rotina deve ser utilizada com a estrutura de dados    *
c *         definida no arquivo 'common.h'                             *
c *                                                                    *
c **********************************************************************
      use Malloc
      implicit none
      integer nnode,numel,nen
c ... ponteiros      
      integer*8 i0,i1,i2
c .....................................................................      
      integer ipos,ix(nen+1,*)
c ......................................................................
c
c ... Grafo da malha:
c
c     ia(i0) => ip(nnode+1) - ip(i) indica a posicao em ips do primeiro
c                             no vizinho ao no i.
c     ia(i1) => ips(ipos) - contem as conetividades nodais de cada no
c     ia(i2) => arranjo auxiliar
c
      i0 = alloc_4('iaux0   ',1,nnode+1)
      call mzero(ia(i0),nnode+1)
      call aponta(ix,ia(i0),nnode,numel,nen,ipos)
      i1 = alloc_4('iaux1   ',1,ipos)
      i2 = alloc_4('iaux2   ',1,ipos)
      call mzero(ia(i1),ipos)
      call mzero(ia(i2),ipos)
      call graph0(ix,ia(i0),ia(i1),ia(i2),nnode,numel,nen)
      call graph1(ia(i0),ia(i1),ia(i2),nnode,ipos)      
      i2 = dealloc('iaux2   ')
c ......................................................................      
      return
      end
      subroutine aponta(ix,ip,nnode,numel,nen,ipos)
c **********************************************************************
c *                                                                    *
c *   APONTA: monta o arranjo ip(nnode+1), apontador do grafo inicial. *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *                                                                    *
c *    ix(nen+1,numel) - conetividades nodais                          *
c *    nnode - numero de nos                                           *
c *    numel - numero de elementos                                     *
c *    nen   - numero de nos por elemento                              *
c *                                                                    *
c *   Parametros de saida:                                             *
c *                                                                    *
c *    ip(nnode+1) - ip(i) indica a posicao no grafo do primeiro       *
c *                  no vizinho ao no i                                *
c *    ipos - numero de posicoes do arranjo que contem o grafo         *
c *           inicial, que contem nos repetidos.                       *
c *                                                                    *
c **********************************************************************
      implicit none
      integer nnode,numel,nen,ipos,i,j,no,ix(nen+1,*),ip(nnode+1)
c ......................................................................
      ip(1) = 1
      do 100 i = 1, numel
      do 100 j = 1, nen
         no = ix(j,i)
         if (no .eq. 0 .or. no .gt. nnode) goto 100
         ip(no+1) = ip(no+1) + nen - 1
  100 continue
      do 200 i = 2, nnode+1
         ip(i) = ip(i) + ip(i-1)
  200 continue
      ipos = ip(nnode+1) - 1
c ......................................................................            
      return
      end
      subroutine graph0(ix,ip,ips,noviz,nnode,numel,nen)
c **********************************************************************
c *                                                                    *
c *   GRAPH0: monta o grafo inicial,armazenado no arranjo noviz.       *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *                                                                    *
c *    ix(nen+1,numel) - conetividades nodais                          *
c *    ip(nnode+1) - apontador do grafo inicial                        *
c *    numel - numero de elementos                                     *
c *    nen   - numero de nos por elemento                              *
c *                                                                    *
c *   Parametros de saida:                                             *
c *                                                                    *
c *    ips(nnode+1) = ip(nnode+1) - arranjo auxiliar                   *
c *    noviz(ipos)  - grafo inicial, contendo nos repetidos            *
c *                   noviz() = adjncy() no reordenador RCM            *
c *                                                                    *
c **********************************************************************
      implicit none
      integer ip(*),ips(*),noviz(*),nnode,numel,nen,ix(nen+1,*)
      integer i,j,k,noj,nok,ipos
c ......................................................................
      do 400 i = 1, numel
      do 300 j = 1, nen
         noj = ix(j,i)
         if (noj .eq. 0 .or. noj .gt. nnode) goto 300
         do 200 k = 1, nen
            nok = ix(k,i)
            if (nok .eq. 0 .or. nok .gt. nnode) goto 100
            if (nok .eq. noj) goto 100
            ipos = ip(noj) + ips(noj)
            ips(noj) = ips(noj) + 1
            noviz(ipos) = nok
  100       continue
  200    continue
  300 continue
  400 continue
c ......................................................................        
      return
      end
      subroutine graph1(ip,ips,noviz,nnode,ipos)
c **********************************************************************
c *                                                                    *
c *   GRAPH1: ordena os elementos do arranjo noviz, eliminado os nos   *
c *           repetidos, e armazenando o resultado em ips              *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *                                                                    *
c *    ip(nnode+1) - apontador do grafo inicial                        *
c *    ips(ipos)   = ip()                                              *
c *    noviz(ipos) - grafo inicial                                     *
c *    nnode - numero de nos                                           *
c *                                                                    *
c *   Parametros de saida:                                             *
c *                                                                    *
c *    ip(nnode+1)  - apontador do grafo final                         *
c *    ips()        - grafo final, ips()=adjncy() no reordenador RCM   *
c *                                                                    *
c **********************************************************************      
      implicit none
      integer nnode,ipos,i,j,k,l,no,ip(nnode+1),ips(ipos),noviz(*)
c ......................................................................
      call sortgraph(ip,noviz,nnode)
c ......................................................................
c
c ... Substitui os nos repetidos por zero e atualiza o no. de vizinhos
c
      do 110 i = 1, nnode
         l = ip(i)
         k = ip(i+1)-2
         do 100 j = l, k
            no = noviz(j)           
            if ((no .ne. 0) .and. (noviz(j+1) .eq. no)) then
               noviz(j) = 0
               ips(i) = ips(i) - 1
            endif
  100    continue
  110 continue
c ......................................................................
c
c ... Atualiza o vetor apontador:                                       
c 
      do 200 i = 1, nnode
         ip(i+1) = ip(i) + ips(i)         
  200 continue
c ......................................................................  
c
c ... Elimina os zeros da vizinhanca:                                   
c  
      j = 1
      do 210 i = 1, ipos
         if (noviz(i) .ne. 0) then
            ips(j) = noviz(i)
            j = j+1
         endif
  210 continue
c ......................................................................
      return
      end      
      subroutine sortgraph(ia,ja,neq)
c **********************************************************************
c *                                                                    *
c *   SORTGRAPH: ordena o grafo armazenado em ja() em ordem crescente  *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *                                                                    *
c *    ia(neq+1) - ia(i) informa a posicao no vetor a do primeiro      *
c *                      coeficiente nao-nulo da equacao i             *
c *    ja(nad) - ja(k) informa a coluna do coeficiente que ocupa       *
c *              a posicao k no vetor a                                *
c *                                                                    *
c *   Parametros de saida:                                             *
c *                                                                    *
c *    ja(nad) - grafo ordenado                                        *
c *                                                                    *
c **********************************************************************
      implicit none
      integer ia(*),ja(*),neq,i
c ......................................................................               
      do 100 i = 1,neq
         call bubblesort(ja(ia(i)),ia(i+1)-ia(i))
  100 continue
c ......................................................................           
      return
      end
      subroutine bubblesort(ja,n)
c **********************************************************************
c *                                                                    *
c *   BUBBLESORT: ordena o arranjo ja(n), em ordem crescente.          *
c *                                                                    *
c **********************************************************************
      implicit none
      integer ja(*),n,i,j,itroca
c ......................................................................               
   50 continue
      itroca = 0
      do 100 i = 2, n
         if(ja(i) .lt. ja(i-1)) then
            j = ja(i-1)
            ja(i-1) = ja(i)
            ja(i) = j
            itroca = 1
         endif
  100 continue
      if (itroca .eq. 1) goto 50
c ......................................................................               
      return
      end
