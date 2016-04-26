c*****************************Svn***************************************      
c*$Date: 2011-03-16 15:32:53 -0300 (Wed, 16 Mar 2011) $                 
c*$Rev: 914 $                                                           
c*$Author: henrique $                                                   
c***********************************************************************      
      subroutine edgstruct(ix,nnode,numel,nen,nedge,i_edg,i_ipedg)
c **********************************************************************
c *                                                                    *
c *   EDGSTRUCT: monta a estrutura de dados para o formato ARESTAS.    *
c *   ---------                                                        *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ---------------------                                            *
c *                                                                    *
c *    ix(nen+1,numel) - conetividades nodais                          *
c *    nnode - numero de nos                                           *
c *    numel - numero de elementos                                     *
c *    nen   - numero de nos por elemento                              *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   -------------------                                              *
c *                                                                    *
c *    nedge - numero de arestas                                       *
c *    i_edg - ponteiro para o arranjo edge(2,nedge)                   *
c *            edge(1,k) = no1 da aresta k                             *
c *            edge(2,k) = no2 da aresta k (no1 < no2)                 *
c *    i_ipedg - ponteiro para o arranjo ip(nnode+1)                   *
c *              ip(i) = primeira aresta para a qual no1 = i           *
c *                                                                    *
c **********************************************************************
      use Malloc
      implicit none
      integer nen,ix(nen+1,*),nnode,numel,nedge
c ... ponteiros      
      integer*8 i_edg,i_ipedg
      integer*8 i0,i1
c ......................................................................
c
c ... Grafo da malha:
c
c     ia(i0) => ip(nnode+1) - ip(i) indica a posicao em ips do primeiro
c                             no vizinho ao no i.
c     ia(i1) => ips(ipos) - contem as conetividades nodais de cada no
      call graph(ix,nnode,numel,nen,i0,i1)
c ......................................................................      
c
c ... Numero de arestas:
c
      call numedges(ia(i0),ia(i1),nnode,nedge)
c ......................................................................
c
c ... Determinacao das arestas:
c
      i_edg   = alloc_4('edge    ',2,nedge)
      i_ipedg = alloc_4('ipedge  ',1,nnode+1)
      call edges(ia(i0),ia(i1),ia(i_ipedg),ia(i_edg),nnode)     
      i1      = dealloc('iaux1   ')
      i0      = dealloc('iaux0   ')
      i_edg   =  locate('edge    ')
      i_ipedg =  locate('ipedge  ')
c ......................................................................      
      return
      end
      subroutine numedges(ip,ips,nnode,nedge)
c **********************************************************************
c *                                                                    *
c *   NUMEDGES: calcula o numero de arestas do grafo                   *
c *   --------                                                         *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ---------------------                                            *
c *                                                                    *
c *    ip(nnode+1) - apontador do grafo                                *
c *    ips(ipos)   - grafo da malha                                    *
c *    nnode - numero de nos                                           *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   -------------------                                              *
c *                                                                    *
c *    nedge - numero de arestas do grafo                              *
c *                                                                    *
c **********************************************************************      
      implicit none
      integer ip(*),ips(*),nnode,nedge,i,j,k
c ......................................................................
      nedge = 0
      do 110 i = 1, nnode
         do 100 k = ip(i), ip(i+1)-1
            j = ips(k)
            if (j .gt. i) then
               nedge = nedge + 1
            endif
  100    continue
  110 continue
c ......................................................................  
      return
      end
      subroutine edges(ip,ips,ipedge,edge,nnode)
c **********************************************************************
c *                                                                    *
c *   EDGES: determina as conetividades das arestas                    *
c *   -----                                                            *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ---------------------                                            *
c *                                                                    *
c *    ip(nnode+1) - apontador do grafo                                *
c *    ips(ipos)   - grafo da malha                                    *
c *    nnode       - numero de nos                                     *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   -------------------                                              *
c *                                                                    *
c *    ipedge(nnode+1) - apontador de arestas                          *
c *          ipedge(i) = primeira aresta para a qual no1 = i           *
c *    edge(2,nedge)   - conetividades das arestas                     *
c *          edge(1,k) = no1 da aresta k                               *
c *          edge(2,k) = no2 da aresta k (no1 < no2)                   *
c *                                                                    *
c **********************************************************************      
      implicit none
      integer ip(*),ips(*),ipedge(*),edge(2,*),nnode,i,j,k,n
c ......................................................................
      n = 0
      do 110 i = 1, nnode
         ipedge(i) = n+1
         do 100 k = ip(i), ip(i+1)-1
            j = ips(k)
            if (j .gt. i) then
               n = n + 1
               edge(1,n) = i
               edge(2,n) = j
            endif
  100    continue
  110 continue
      ipedge(nnode+1) = ipedge(nnode)
c ......................................................................  
      return
      end                              
