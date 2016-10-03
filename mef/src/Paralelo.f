      subroutine noeq(nnof,ndf,neqf,fmap0,fmap,idr,noGL,rcvs,rcvs0,
     .                dspl0,nviz)
c **********************************************************************
c *                                                                    *
c *   NOEQ                                                             *
c *   ----                                                             *
c *   Converte mapa Interface -> Global de no's                        *
c *        em  mapa Interface -> Local de equacoes                     *
c *                                                                    *
c *   Parâmetros de entrada:                                           *
c *   ----------------------                                           *
c *   nnof  -  numero de no's em comunicacao                           *
c *   ndf   -  numero de graus de liberdade por no'                    *
c *   nnode -  numero de no's local                                    *
c *   neqf  -  indefinido                                              *
c *   fmap0(nnof)     - mapa Interface->Global de no's                 *
c *   fmap(ndf,nnof)  - indefinido                                     *
c *   idr(ndf,nnode)  - numeracao local de equacoes                    *
c *   noGL(nnoG)      - mapa Global->Local de no's                     *
c *                                                                    *
c *   Parâmetros de saida:                                             *
c *   -------------------                                              *
c *   neqf - numero de equacoes em comunicacao                         *
c *   fmap(neqf) - mapa Interface->Local de equacoes                   *
c *                                                                    *
c **********************************************************************
      implicit none
      integer nnof,ndf,i,j,k,m,n,neqf,rcvs(*),rcvs0(*),dspl0(*)
      integer fmap0(*),idr(ndf,*),fmap(*),noGL(*),viz,nviz
c ......................................................................
      call mzero(fmap,nnof*ndf)
      rcvs(1:nviz) = 0
      neqf = 0
c ... viz = loop nos vizinhos
      do viz = 1, nviz
c ...    i = loop nos nohs de fronteira do vizinho
         do i = dspl0(viz)+1, dspl0(viz) + rcvs0(viz)
c ...       j = numeração global do nó
            j = fmap0(i)    
c ...       n = numeração local do nó
            n = noGL(j)
c ...       k = loop nos graus de liberdade do nó
            do k = 1, ndf
c ...          m = equação local do g.l.
               m = idr(k,n)
c ...          se o g.l. tiver equação
               if (m .gt. 0) then
c ...             neqf = numeração de equações de interface
                  neqf = neqf+1
c ...             mapa fronteira -> local
                  fmap(neqf) = m
c ...             estruturas para Irecv
                  rcvs(viz) = rcvs(viz) + 1
               endif
            enddo
         enddo
      enddo
c ......................................................................
      return
      end
c ......................................................................
      subroutine communicate(x,neqf1,neqf2,i_fmap,i_xf,
     .           i_rcvs,i_dspl)
c **********************************************************************
c *                                                                    *
c *   communicate: comunicação de vetores reais                        *
c *   -----------                                                      *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ---------------------                                            *
c *    x(0:neq)    -  vetor local destualizado                         *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   -------------------                                              *
c *    x(0:neq)    -  vetor local corrigido                            *
c *                                                                    *
c **********************************************************************
      use Malloc
      implicit none
      include 'mpif.h'
      include 'parallel.fi'                        
      include 'time.fi'
      real*8  x(*)
      integer i,j,k,l,p,d,c,e,status(MPI_STATUS_SIZE,2*nprcs)
c ......................................................................
      if(mpi .eqv. .false.) return  
c
c ... Comunicacao de vetores
c
      time0 = MPI_Wtime()
c
c ... Preenche buffer de envio:
      call mapfront(x,ia(i_xf+2*neqf1),ia(i_fmap+neqf1),neqf2)
      ovhtime = ovhtime + MPI_Wtime()-time0
c ......................................................................
      time0 = MPI_Wtime()
c
c ... Envia dados do ciclo atual:
c
c ... loop nos processos k, vizinhos a my_id:
      do k = 0, nviz1-1
         j = ia(i_ranks+k)
c
c ...    recebe do k-esimo vizinho:
         call MPI_Irecv(ia(i_xf+2*ia(i_dspl+k)),ia(i_rcvs+k),mdp,j,
     .                  my_id,mcw,ia(i_rreqs+k),ierr)
      enddo
      do k = nviz1, nviz1+nviz2-1
         j = ia(i_ranks+k)
c
c ...    envia para o k-esimo vizinho:
         call MPI_Isend(ia(i_xf+2*neqf1+2*ia(i_dspl+k)),ia(i_rcvs+k),
     .                  mdp,j,j,mcw,ia(i_sreqs+k-nviz1),ierr)
      enddo
c
c ... Aguarda recebimento dos dados:
      call MPI_waitall(nviz1,ia(i_rreqs),status,ierr)
      sendtime = sendtime + MPI_Wtime() - time0
c ......................................................................
      time0 = MPI_Wtime()
c
c ... Recupera valores recebidos:
      call maploc(x,ia(i_xf),ia(i_fmap),neqf1,ovlp)
c
c ... Libera requests dos Isend:
      call MPI_waitall(nviz2,ia(i_sreqs),status,ierr)
      ovhtime = ovhtime + MPI_Wtime()-time0
c ......................................................................
      return
      end
      subroutine maploc(x,xf,fmap,neqf,ovlp)
c **********************************************************************
c *                                                                    *
c *   MAPLOC_V   mapeia valores da fronteira para o local da particao  *
c *   --------                                                         *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ---------------------                                            *
c *    xf(neqf)    -  vetor de buffer                                  *
c *    fmap(neqf)  -  mapa de equações da fronteira p/ local           *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   -------------------                                              *
c *    x(0:neq)    -  vetor local corrigido                            *
c *                                                                    *
c **********************************************************************      
      implicit none
      real*8  x(*),xf(*)
      integer fmap(*),neqf,i,j
      logical ovlp
c ......................................................................
      if (ovlp) then
         do i = 1, neqf
            x(fmap(i)) = xf(i)
         enddo
      else
         do i = 1, neqf
            j = fmap(i)
            x(j) = x(j) + xf(i)
         enddo
      endif
c ......................................................................
      return
      end
      subroutine mapfront(x,xf,fmap,neqf)
c **********************************************************************
c *                                                                    *
c *   MAPFRONT   mapeia valores do local para a fronteira da particao  *
c *   --------                                                         *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ---------------------                                            *
c *    x(neq)    -  vetor local                                        *
c *    fmap(neqf)  -  mapa de equações de local p/ fronteira           *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   -------------------                                              *
c *    xf(neqf)    -  vetor de buffer                                  *
c *                                                                    *
c **********************************************************************
      implicit none
      real*8  x(*),xf(*)
c      integer fmap(*),neqf,neqf2,i,j
      integer fmap(*),neqf,i
c ......................................................................
      do i = 1, neqf
         xf(i) = x(fmap(i))
      enddo
c ......................................................................
      return
      end
      subroutine allgatheri(ic,i_xf)
c **********************************************************************
c *                                                                    *
c *   allgatherv0: comunicação de vetor de inteiros                    *
c *   -----------                                                      *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ---------------------                                            *
c *    ic(neq)    -  vetor no local destualizado                       *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   -------------------                                              *
c *    ic(neq)    -  vetor local corrigido                             *
c *                                                                    *
c **********************************************************************
      use Malloc
      implicit none
      include 'mpif.h'
      include 'parallel.fi'                        
      integer  ic(*),j
c ----------------------------------------------------------------------
c
c ... Comunica os {ic(ps)}
c
      call MPI_allgatherv(ic(nno1+1),nno2+nno3,MPI_INTEGER,
     .                   ia(i_xf),ia(i_rcvs0i),ia(i_dspl0i),MPI_INTEGER,
     .                   MPI_COMM_WORLD,ierr)
c ......................................................................
c
c ... Recupera valores recebidos
c
      call maploc_ic(ic,ia(i_xf),ia(i_fmap0i),nno1+1,nnofi,nno2+nno3)
c ......................................................................
      return
      end
      subroutine maploc_ic(ic,xf,fmap,k,n,nl)
c **********************************************************************
c *                                                                    *
c *   MAPLOC_IC  mapeia valores da fronteira para o local da particao  *
c *   ---------                                                        *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ---------------------                                            *
c *    xf(neqf)    -  vetor de buffer                                  *
c *    fmap(neqf)  -  mapa de equações da fronteira p/ local           *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   -------------------                                              *
c *    ic(neq)    -  vetor local corrigido                             *
c *                                                                    *
c **********************************************************************      
      implicit none
      integer ic(*),xf(*)
      integer fmap(*),k,n,nl,i,j
c ......................................................................
      call mzero(ic(k),nl)
      do i = 1, n
         j = fmap(i)
         if (j .gt. 0) ic(j) = ic(j) + xf(i)
      enddo                     
c ......................................................................      
      return
      end
      subroutine global_v(nl,ncl,i_l,i_g,name)
c **********************************************************************
c *                                                                    *
c *   GLOBAL_V   arranjo nodal real global para montagem de saidas     *
c *   --------                                                         *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ----------------------                                           *
c *      nl        -  número de linhas                                 *
c *      ncl       -  número de colunas local                          *
c *      i_l       -  ponteiro do arranjo local                        *
c *      name      -  string com nome para alocacao do arranjo global  *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   -------------------                                              *
c *      i_g  -  ponteiro do arranjo global (apenas processo 0)        *
c *                                                                    *
c **********************************************************************
      use Malloc
      implicit none
      include 'mpif.h'
      include 'parallel.fi'
      integer nl,ncl
c ... ponteiro      
      integer*8 i_l,i_g,i1
c ......................................................................      
      character*8 name
c ......................................................................
      if (nprcs .le. 1) then
         i_g = i_l
         return
      endif
c ......................................................................
c
c ... Cria arranjo global
      i_g = alloc_8(   name   , nl, nnoG)
      i1  = alloc_8('temp    ', nl, nnoG)
c ... Verificar com o iuri esta modificacao:
      call  azero  (ia(i_g),nl*nnoG)
      call  azero  (ia(i1) ,nl*nnoG)      
c ----------------------------------------------------------------------      
c
c     Passa valores do arranjo local i_l para o arranjo global i1
      call map_v(ia(i_l),nl,ncl,ia(i1),ia(i_noLG))
c
c     Soma i1's de todos os processos em i_g, no processo raiz
      call MPI_reduce(ia(i1),ia(i_g),nl*nnoG,MPI_DOUBLE_PRECISION,
     .                MPI_SUM,0,MPI_COMM_WORLD,ierr)
      i1 = dealloc('temp    ')
c ......................................................................
      return
      end
      subroutine map_v(vl,nl,ncl,vg,mapa)
c **********************************************************************
c *                                                                    *
c *   MAP_V     mapeia arranjo real*8 local em global                  *
c *   -----                                                            *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ----------------------                                           *
c *      vl(nl,ncl) -  vetor local                                     *
c *      nl         -  número de linhas                                *
c *      ncl        -  número de colunas local                         *
c *      ncg        -  número de colunas global                        *
c *      mapa(ncl)  -  relacao local->global                           *
c *                                                                    *
c *   Parametro de saida:                                              *
c *   ------------------                                               *
c *      vg(nl,ncg) -  vetor global                                    *
c *                                                                    *
c **********************************************************************
      implicit none
      integer  nl,ncl,mapa(*)
      real*8   vl(nl,*),vg(nl,*)
      integer i,j,k
c ......................................................................
      do i = 1, ncl
         k = mapa(i)
         do j = 1, nl
            vg(j,k) = vl(j,i)
         enddo
      enddo
c ......................................................................
      return
      end
      subroutine global_ix(nl,ncl,i_l,i_g,name)
c **********************************************************************
c *                                                                    *
c *   GLOBAL_IX    arranjo inteiro global para montagem de saidas      *
c *   ---------                                                        *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ----------------------                                           *
c *      nl        -  número de linhas                                 *
c *      ncl       -  número de colunas local                          *
c *      i_l       -  ponteiro do arranjo local                        *
c *      name      -  string com nome para alocacao do arranjo global  *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   -------------------                                              *
c *      i_g  -  ponteiro do arranjo global (apenas processo 0)        *
c *                                                                    *
c **********************************************************************
      use Malloc
      implicit none
      include 'mpif.h'
      include 'parallel.fi'
      integer nl,ncl
c ... ponteiro      
      integer*8 i_l,i_g,i1
c ......................................................................      
      character*8 name
c ......................................................................
      if (nprcs .le. 1) then
         i_g = i_l
         return
      endif
c ......................................................................
c
c ... Cria arranjo global
      i_g = alloc_4(   name   , nl, nelG)
      i1  = alloc_4('temp    ', nl, nelG)
      call  mzero(ia(i_g)  ,nl*nelG)
      call  mzero(ia(i1)   ,nl*nelG)      
c
c     Passa valores do arranjo local i_l para o arranjo global i1
      call map_ix(ia(i_l),nl,ncl,ia(i1),ia(i_elLG),ia(i_noLG))
c
c     Soma i1's de todos os processos em i_g, no processo raiz
      call MPI_reduce(ia(i1),ia(i_g),nl*nelG,MPI_INTEGER,
     .                MPI_SUM,0,MPI_COMM_WORLD,ierr)
      i1 = dealloc('temp    ')
c ......................................................................
      return
      end
c **********************************************************************
c
c **********************************************************************
c *                                                                    *
c *   GLOBAL_V_ELM : Monta o a saida global da variaveis por elemento  *
c *   -------------                                                    *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ----------------------                                           *
c *      nl        -  número de linhas                                 *
c *      ncl       -  número de colunas local                          *
c *      i_l       -  ponteiro do arranjo local                        *
c *      name      -  string com nome para alocacao do arranjo global  *
c *      nelG      -  numero global de elementos                       *
c *      cod       -  1 - inteiro de 4 bytes                           *
c *                   2 - real de 8 bytes                              *   
c *                                                                    *
c *   Parametros de saida:                                             *
c *   -------------------                                              *
c *      i_g  -  ponteiro do arranjo global (apenas processo 0)        *
c *                                                                    *
c **********************************************************************
      subroutine global_v_elm(nl,ncl,i_l,i_g,name,cod)
      use Malloc
      implicit none
      include 'mpif.h'
      include 'parallel.fi'
      integer nl,ncl
      integer cod
c ... dummy var      
      integer idum
      real*8  ddum
c ... ponteiro      
      integer*8 i_l,i_g,i1
c ......................................................................      
      character*8 name
c ......................................................................
      if (nprcs .le. 1) then
         i_g = i_l
         return
      endif
c ......................................................................
c
c ... Cria arranjo global
      if(cod .eq. 1 ) then
        i_g = alloc_4(   name   , nl, nelG)
        i1  = alloc_4('g_v_elm ', nl, nelG)
        call  mzero(ia(i_g)  ,nl*nelG)
        call  mzero(ia(i1)   ,nl*nelG)      
c
c     Passa valores do arranjo local i_l para o arranjo global i1
        call map_v_elm(ia(i_l),ddum,nl,ncl,ia(i1),ddum,ia(i_elLG),cod)
c
        call MPI_reduce(ia(i1),ia(i_g),nl*nelG,MPI_INTEGER,
     .                  MPI_SUM,0,MPI_COMM_WORLD,ierr)
      elseif( cod .eq. 2) then
        i_g = alloc_8(   name   , nl, nelG)
        i1  = alloc_8('g_v_elm ', nl, nelG)
        call  azero(ia(i_g)  ,nl*nelG)
        call  azero(ia(i1)   ,nl*nelG)      
c     Passa valores do arranjo local i_l para o arranjo global i1
        call map_v_elm(idum,ia(i_l),nl,ncl,idum,ia(i1),ia(i_elLG),cod)
c
        call MPI_reduce(ia(i1),ia(i_g),nl*nelG,MPI_DOUBLE_PRECISION,
     .                  MPI_SUM,0,MPI_COMM_WORLD,ierr)
      endif
c ......................................................................
c
      i1 = dealloc('g_v_elm ')
c ......................................................................
      return
      end
c **********************************************************************
c
c **********************************************************************
c *                                                                    *
c *   MAP_V_ELM : mapeia arranjo v  local em global de elementos       *
c *   ---------                                                        *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ----------------------                                           *
c *      v_l_i(ndf,numel) -  grandeza v local inteira                  *
c *      v_l_d(ndf,numel) -  grandeza v local real                     *
c *      ndf              -  numero de graus de liberdade              *
c *      numel            -  número de elmentos                        *
c *      v_G_i(ndf,numel) -  nao definido                              *
c *      v_G_d(ndf,numel) -  nao definido                              *
c *      elLG(numel)      -  relacao local->global de elementos        *
c *      cod              -  1 - inteiro de 4 bytes                    *
c *                          2 - real de 8 bytes                       *   
c *                                                                    *
c *   Parametro de saida:                                              *
c *   ------------------                                               *
c *      v_G_i(ndf,numel) -  v global  inteitro                        *
c *                            ou                                       *
c *      v_G_d(ndf,numel) -  v global  real                            *
c **********************************************************************
      subroutine map_v_elm(v_l_i,v_l_d,ndf,numel,v_G_i,v_G_d,elLG,cod)
      implicit none
      integer  ndf,numel,elLG(*)
      real*8   v_l_d(ndf,*),v_G_d(ndf,*)
      integer  v_l_i(ndf,*),v_G_i(ndf,*)
      integer i,j,k
      integer cod
c ......................................................................
c
c ... variaveis inteiras
      if( cod .eq. 1) then
        do i = 1, numel
           k = elLG(i)
           do j = 1, ndf
             v_G_i(j,k) = v_l_i(j,i)
           enddo
        enddo
c ... variaveis reais de dupla precisao
      elseif(cod .eq. 2) then
        do i = 1, numel
           k = elLG(i)
           do j = 1, ndf
             v_G_d(j,k) = v_l_d(j,i)
           enddo
        enddo
      endif  
c ......................................................................
      return
      end
c **********************************************************************
c
c **********************************************************************
      subroutine map_ix(ix_l,nen1,numel,ix_G,elLG,noLG)
c **********************************************************************
c *                                                                    *
c *   MAP_IX     mapeia arranjo ix local em global                     *
c *   ------                                                           *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ----------------------                                           *
c *      ix_l(nen1,numel) -  ix local                                  *
c *      nen1             -  número de linhas                          *
c *      numel            -  número de colunas local                   *
c *      ix_G(nen1,numel) -  relacao local->global                     *
c *      elLG(numel)      -  relacao local->global de elementos        *
c *      noLG(nnode)      -  relacao local->global de no's             *
c *                                                                    *
c *   Parametro de saida:                                              *
c *   ------------------                                               *
c *      ix_G(nen1,ncg)   -  ix global                                 *
c *                                                                    *
c **********************************************************************
      implicit none
      integer  nen1,numel,elLG(*),ix_l(nen1,*),ix_G(nen1,*),noLG(*)
      integer i,j,k
c ......................................................................
      do i = 1, numel
         k = elLG(i)
         do j = 1, nen1-1
            ix_G(j,k) = noLG(ix_l(j,i))
         enddo
         ix_G(nen1,k) = ix_l(nen1,i)
      enddo
c ......................................................................
      return
      end
c ......................................................................
      subroutine prt(v,n)
      real*8 v(*)
      integer n,i
      write(*,'(9999f6.0)') (v(i),i=1,n)
      return
      end
c ......................................................................
      subroutine twrite(text,ts,nprcs,nout)
      implicit none
      character*6 text
      real*8 ts(*)
      integer i,nprcs,nout
      write(nout,99) text,(ts(i),i=1,nprcs)
      return
   99 format(a6,128f18.6)      
      end
c ......................................................................
      subroutine itwrite(text,ts,nprcs,nout)
      implicit none
      character*6 text
      integer ts(*)
      integer i,nprcs,nout
      write(nout,99) text,(ts(i),i=1,nprcs)
      return
   99 format(a6,128i16)
      end
c **********************************************************************
c
c **********************************************************************
c *                                                                    *
c *   init_front : inicia a estrutura do fronte                        *
c *   ------                                                           *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ----------------------                                           *
c *      noLG(nnode)      -  arranjo de nos local -> global            *
c *      noGL             -  arranjo de nos global -> local            *
c *      nno1             -  numero de nos V1                          *
c *      nno2             -  numero de nos V2                          *
c *      nno3             -  numero de nos V3                          *
c *      nnofi            -  nao definido                              *
c *      nno_pload        -  nao definino                              *
c *      nnoG             -  nao definido                              *
c *      nelG             -  nao definido                              *
c *      ovlp             -  nao definido                              *
c *      novlp            -  nao definido                              *
c *      nprcs            -  numero de processo do Mpi                 *
c *      nviz1            -  numeros de vizinhos recebimento           *
c *      nviz2            -  numeros de vizinhos envio                 *
c *      i_rreqs          -  nao definido                              *
c *      i_rreqs          -  nao definido                              *
c *      i_sreqs          -  nao definido                              *
c *      i_rcvs0i         -  nao definido                              *
c *      i_dspl0i         -  nao definido                              *
c *      i_fmap0i         -  nao definido                              *
c *                                                                    *
c *   Parametro de saida:                                              *
c *   ------------------                                               *
c *      nno_pload        -  nnode (sem Mpi)                           *
c *      nnoG             -  nnode (sem Mpi)                           *
c *      nelG             -  numel (sem Mpi)                           *
c *      ovlp             -  .false.(sem Mpi)                          *
c *      novlp            -  .false. (sem Mpi)                         *
c *      i_rrqes          - estrutura de comunicacao do mpi            *
c *      i_srqes          - estrutura de comunicacao do mpi            *
c *      i_rcvs0i         - mapa para os nos                           *
c *      i_dspl0i         - mapa para os nos                           *
c *      i_fmap0i         - mapa para nos                              *
c *      nnofi            - numero de nos com comunicacao              *
c *      mpi              - excucao em mpi                             *
c **********************************************************************
      subroutine init_front(i_noLG,i_noGL,nno1,nno2,nno3,nnofi,nno_pload
     .                     ,nnoG,nelG,nnode,numel,ovlp,novlp,nprcs,nviz1
     .                     ,nviz2,i_rreqs,i_sreqs,i_rcvs0i,i_dspl0i
     .                     ,i_fmap0i,mpi) 
c ===     
      use Malloc
      implicit none
      include 'mpif.h'
      integer nno1,nno2,nno3,nnofi
      integer nno_pload,nnoG,nelG
      logical ovlp,novlp
      integer nprcs,nviz1,nviz2
      integer nnode,numel
      logical mpi
c ... ponteiros      
      integer*8 i_noLG,i_noGL
      integer*8 i_rreqs,i_sreqs
      integer*8 i_fmap0i,i_dspl0i,i_rcvs0i
      integer*8 i_fmap2
c .....................................................................
      integer ierr
      integer i,k
c =====================================================================
c
c === sem Mpi
      if (nprcs .eq. 1) then
        nno_pload = nnode
        nnoG      = nnode
        nelG      = numel
        ovlp      = .false.
        novlp     = .false.        
        mpi       = .false.        
        return
      endif
c =====================================================================
c
c === com Mpi
      mpi = .true.
      nno_pload = nno1 + nno2
c        Estruturas de comunicacao para send-receive 
c ... {rreqs} = receive requests:
      i_rreqs = alloc_4('rreqs   ', 1, nviz1)
c
c ... {sreqs} = send requests:
      i_sreqs = alloc_4('sreqs   ', 1, nviz2) 
c       Adapta arranjos para allgatheri
c
c ... {rcvs0}
      i_rcvs0i = alloc_4('rcvs0i  ', 1, nprcs)
      call mzero(ia(i_rcvs0i),nprcs)
      call MPI_allgather(nno2+nno3,1,MPI_INTEGER,ia(i_rcvs0i),1,
     .                   MPI_INTEGER,MPI_COMM_WORLD,ierr)
c ... {dspl0}
      i_dspl0i = alloc_4('dspl0i  ', 1, nprcs)
      ia(i_dspl0i) = 0
      do i = 1, nprcs-1
         ia(i_dspl0i+i) = ia(i_dspl0i+i-1) + ia(i_rcvs0i+i-1)
      enddo
c ... {fmap0}
      nnofi = ia(i_dspl0i+nprcs-1)+ia(i_rcvs0i+nprcs-1)
      i_fmap0i = alloc_4('fmap0i  ', 1, nnofi)
      i_fmap2  = alloc_4('fmap2   ', 1, nno2+nno3)
      call mzero(ia(i_fmap0i),nnofi)
      call mzero(ia(i_fmap2),nno2+nno3)
      do i=1,nno2+nno3
         ia(i_fmap2-1+i) = ia(i_noLG-1+nno1+i)
      enddo
      call MPI_allgatherv(ia(i_fmap2),nno2+nno3,MPI_INTEGER,
     .                    ia(i_fmap0i),ia(i_rcvs0i),ia(i_dspl0i),
     .                    MPI_INTEGER,MPI_COMM_WORLD,ierr)
      i_fmap2 = dealloc('fmap2   ')
      do i=1, nnofi
         ia(i_fmap0i-1+i+k) = ia(i_noGL-1+ia(i_fmap0i-1+i))
c         ia(i_fmap0i-1+i+k) = noGL(ia(i_fmap0i-1+i))
      enddo
c =====================================================================
c
c === 
      return
      end
c =====================================================================
c **********************************************************************
c     
c **********************************************************************
      subroutine frontb(ndf,idr,neqi,neq1i,neq2i,neq3i
     .                 ,neq4i,neq1ai,neqf1i,neqf2i,neq32i,neq_doti
     .                 ,i_fmapi,i_rcvsi,i_dspli,i_xfi
     .                 ,namefmap,namercvs,namedspl,namexf,cod)
c **********************************************************************
c *                                                                    *
c *   FRONTB:                                                          *
c *   -----                                                            *
c *       Calcula neqs por tipo de no' e                               *
c *       Transforma o mapa Interface->Global de no's em               *
c *       mapa Interface -> Local de equacoes                          *
c *                                                                    *
c *   Parâmetros de entrada:                                           *
c *   ----------------------                                           *
c *      ndf            - numero de graus de liberdade por no'         *
c *      idr(ndf,nnode) - numeracao local de equacoes                  *
c *      neqi           -                                              *
c *      neq1i          -                                              *
c *      neq2i          -                                              *
c *      neq3i          -                                              *
c *      neq4i          -                                              *
c *      neq1ai         -                                              *
c *      neqf1i         -                                              *
c *      neqf2i         -                                              *
c *      noGL(nnoG)     - mapa Global->Local de no's                   *
c *      neq            - numero de equacoes de numeq                  *
c *      namefmap       - nome do arranjo fmap                         *
c *      namercvs       - nome do arranjo rcvs                         *
c *      namedspl       - nome do arranjo dspl                         *
c *      namexf         - nome do arranjo xf                           *
c *      cod            - imprime o mapa em um arquivo auxiliar        *
c *                       1 - mecanico                                 *
c *                       2 - termico                                  *
c *   Parâmetros de saida:                                             *
c *   -------------------                                              *
c *      neqi  - numero de equacoes do sistema                         *
c *      neq1i - numero de equacoes em V1                              *
c *      neq2i - numero de equacoes em V2                              *
c *      neq3i - numero de equacoes em V3                              *
c *      neq4i - numero de equacoes em V4                              *
c *      neq1ai- numero de equacoes em V1a                             *
c *      neqf1i- numero de equacoes recebidas na comunicacao           *
c *      neqf2i- numero de equacoes enviadas na comunicacao            *
c *      fmap(neqf) - correlacao neq_interface -> neq_local            *
c *      rcvs(nprcs) - arranjo usado no allgatherv                     *
c *      dspl(nprcs) - arranjo usado no allgatherv                     *
c *      xf  (*    ) - buffer de comunicao                             *
c *                                                                    *
c **********************************************************************
      use Malloc
      implicit none
      include 'mpif.h'
      include 'parallel.fi'
c ... ponteiros      
      integer*8 i_fmapi,i_rcvsi,i_dspli,i_xfi,i_fmap2
c .....................................................................
c
c ... eqs
      integer neq1i,neq2i,neq3i,neq4i,neq1ai,neq32i,neq_doti,neqi
      integer neqf1i,neqf2i
c .....................................................................      
      character*8 namefmap,namedspl,namercvs,namexf
      integer ndf,i,j,k,n,cod
      integer idr(ndf,*)
c ......................................................................
      neq1i  = 0
      neq2i  = 0
      neq3i  = 0
      neq4i  = 0
      neq32i = 0
      neq1ai = 0
      if (nprcs .eq. 1) then
        neq_doti  = neqi
        return
      endif
c ......................................................................
c
c ... Numero de equacoes por grupo de no's
c
c
c ... neq1
      do i = 1, nno1
        do k = 1, ndf
          n = idr(k,i)
          if (n .ne. 0) neq1i = neq1i + 1
        enddo
      enddo
c
c ... neq2
      do i = nno1+1, nno1+nno2
        do k = 1, ndf
          n = idr(k,i)
          if (n .ne. 0) neq2i = neq2i + 1
        enddo
      enddo
c
c ... neq3
      do i = nno1+nno2+1, nno1+nno2+nno3
        do k = 1, ndf
          n = idr(k,i)
          if (n .ne. 0) neq3i = neq3i + 1
        enddo
      enddo
c.............................
      if (ovlp) then
c
c ...   neq4
        do i = nno1+nno2+nno3+1,nno1+nno2+nno3+nno4
          do k = 1, ndf
            n = idr(k,i)
            if (n .ne. 0) neq4i = neq4i + 1
          enddo
        enddo
c
c ...   neq1a
        do i = nno1-nno1a+1, nno1
          do k = 1, ndf
            n = idr(k,i)
            if (n .ne. 0) neq1ai = neq1ai + 1
          enddo
        enddo
      endif
c ......................................................................
c
c ... Parametros especificos de cada metodo a serem usados na solucao:
c
      neq_doti   = neq1i + neq2i
      neq32i     = neq3i
      if (ovlp) then
         neqi    = neq_doti
      else
         neq3i  = 0
      endif
c
      mdp = MPI_DOUBLE_PRECISION
      mcw = MPI_COMM_WORLD
      min = MPI_INTEGER
c ......................................................................
c
c        Estruturas de comunicacao para send-receive 
c
c ... {rcvs} = lista dos tamanhos dos blocos de eqs {Vfi}
      i_rcvsi = alloc_4(namercvs  , 1, nviz1+nviz2)
c
c ... {dspl} = Ponteiro p/ início de cada bloco em {Vf}
      i_dspli = alloc_4(namedspl  , 1, nviz1+nviz2)
c
c ......................................................................
c
c ... {fmap2} = mapa de equacoes Interface -> Local temporario (recv)
c ... {rcvs}
c     
      i_fmap2 = alloc_4('fmap2   ', ndf, nnof1+nnof2)
      call noeq(nnof1,ndf,neqf1i,ia(i_fmap0),ia(i_fmap2),idr,ia(i_noGL)
     .         ,ia(i_rcvsi),ia(i_rcvs0),ia(i_dspl0),nviz1)
c ... {dspl}
      ia(i_dspli) = 0
      do i = 1, nviz1-1
         ia(i_dspli+i) = ia(i_dspli+i-1) + ia(i_rcvsi+i-1)
      enddo
c ......................................................................
c
c ... {fmap2} = mapa de equacoes Interface -> Local temporario (send)
c ... {rcvs}
c
      if(ovlp)then
         call noeq(nnof2,ndf,neqf2i,ia(i_fmap0+nnof1),ia(i_fmap2+neqf1i)
     .            ,idr,ia(i_noGL),ia(i_rcvsi+nviz1),ia(i_rcvs0+nviz1)
     .            ,ia(i_dspl0+nviz1),nviz2)
c ...... {dspl}
         ia(i_dspli+nviz1) = 0
         do i = nviz1+1, nviz1+nviz2-1
            ia(i_dspli+i) = ia(i_dspli+i-1) + ia(i_rcvsi+i-1)
         enddo
      else
         call icopy(i_fmap2,i_fmap2+neqf1i,i_fmap2+neqf1i)
         call icopy(i_rcvsi,i_rcvsi+nviz1 ,i_rcvsi+nviz1)
         call icopy(i_dspli,i_dspli+nviz1 ,i_dspli+nviz1)
         call icopy(i_ranks,i_ranks+nviz1 ,i_ranks+nviz1)
         neqf2i = neqf1i
      endif
c ......................................................................
c
c ... {fmap} = mapa de equacoes Interface -> Local  (arranjo definitivo)
c
      j = neqf1i+neqf2i
      i_fmapi = alloc_4(namefmap  , 1, j)
      call icopy(i_fmap2,i_fmap2+j,i_fmapi)
c ......................................................................
c
c                Libera memoria de arranjos desnecessarios
c
      i_fmap2 = dealloc('fmap2   ')
c
c ... Localiza ponteiros
      i_fmapi  = locate(namefmap)
c ......................................................................
c
c ... Buffer de comunicacao:
c
      k = max(nnofi/2+1,neqf1i+neqf2i)
      i_xfi = alloc_8(namexf, 1, k)
      call azero(ia(i_xfi),k)
c ......................................................................
c
c ... escreve os mapas em um arquivo
c ... mecanico
      if( cod .eq. 1) call printmap(ia(i_fmapi),ia(i_rcvsi),ia(i_dspli)
     .              ,neqi,neq1i,neq2i,neq3i,neq4i,neq1ai,neq32i,neq_doti
     .              ,neqf1i,neqf2i,nviz1,nviz2,'mecanico ',16,my_id)
c ... termico     
      if( cod .eq. 2) call printmap(ia(i_fmapi),ia(i_rcvsi),ia(i_dspli)
     .              ,neqi,neq1i,neq2i,neq3i,neq4i,neq1ai,neq32i,neq_doti
     .              ,neqf1i,neqf2i,nviz1,nviz2,'termico  ',16,my_id)
c ......................................................................     
      return
      end
c *********************************************************************
c                                                                       
c *********************************************************************
c * PRINTMAP : saida do mapa criado pelo front                        *
c * ----------------------------------------------------------------- *
c * Paramentros de entrada:                                           *
c * ----------------------------------------------------------------- *
c * fmap(*) - mapa das equacaoes                                      *
c * rcvs(*) - numero de equacoes por vizinho                          *
c * dspl(*) -                                                         *
c * neq     - numero total de equacoes                                *
c * neq1    - numero V1                                               *
c * neq2    - numero V2                                               *
c * neq3    - numero V3                                               *
c * neq4    - numero V4                                               *
c * neq1a   -                                                         *
c * neq32   -                                                         *
c * neq_dot - numero de equacoes no produto escalar                   *
c * neqf1   - numero total de equacoes no buffer de recebimento       *
c * neqf2   - numero total de equacoes no buffer de envio             *
c * ----------------------------------------------------------------- *
c *********************************************************************
      subroutine printmap(fmap,rcvs,dspl,neq,neq1,neq2,neq3,neq4,neq1a,
     .           neq32,neq_dot,neqf1,neqf2,nviz1,nviz2,prefixo,
     .           cod,my_id)
      implicit none
      character*80 aux,name
      character prefixo*9
      integer my_id
      integer neq,neq1,neq2,neq3,neq4,neq1a,neq32,neq_dot
      integer neqf1,neqf2
      integer nviz1,nviz2
      integer fmap(*),rcvs(*),dspl(*)
      integer cod
      integer i
      integer nout
      data nout/1000/
c =====================================================================
c
c ===
      aux = prefixo
      aux = name(aux,my_id,cod) 
      open(nout,file=aux)
c ...      
      write(nout,'(a)')"#mapa gerado pelo fronte" 
      write(nout,'(a,i6)')"my_id",my_id 
c ...      
      write(nout,'(a6)')"numero eqs" 
      write(nout,100)"neq    ",neq 
      write(nout,100)"neq1   ",neq1 
      write(nout,100)"neq2   ",neq2 
      write(nout,100)"neq3   ",neq3 
      write(nout,100)"neq4   ",neq4 
      write(nout,100)"neq1a  ",neq1a 
      write(nout,100)"neq32  ",neq32 
      write(nout,100)"neq_dot",neq_dot 
      write(nout,100)"neqf1  ",neqf1 
      write(nout,100)"neqf2  ",neqf2 
c ....................................................................       
      write(nout,'(a)')"fmap"
      do i = 1, neqf1+neqf2
        write(nout,200)i,fmap(i)
      enddo 
      write(nout,'(a)')"end fmap"
      write(nout,'(a)')"rcvs"
      do i = 1, nviz1+nviz2
        write(nout,200)i,rcvs(i)
      enddo
      write(nout,'(a)')"end rcvs"
      write(nout,'(a)')"dspl"
      do i = 1, nviz1+nviz2
        write(nout,200)i,dspl(i)
      enddo
      write(nout,'(a)')"end dspl"
      close(nout)
c ====================================================================
c
c ===
100   format(a7,i8)
200   format(i8,i8)
      return
      end
c =====================================================================
c *********************************************************************
c
c *********************************************************************
c * STOP_MEF: aborta a excucao do mefpar encerrando devidamente o MPI *
c *********************************************************************
      subroutine stop_mef()
      implicit none
      include 'mpif.h'
      integer ierr
      call Mpi_finalize(ierr)
      stop
      return
      end
c *********************************************************************
