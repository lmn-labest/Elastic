c *********************************************************************
c * PARTITION_MATRIX : divide o trabalho do matevec entre as threads  *
c * por linhas e inicializa a estruturas do buffer do matvec          *
c * ----------------------------------------------------------------- *
c * parametros de entrada :                                           *
c * ----------------------------------------------------------------- *
c *   ia(neq+1) - ia(i) informa a posicao no vetor au do primeiro     *
c *                     coeficiente nao-nulo da linha   i             *
c *   ja(nad)   - ja(k) informa a coluna do coeficiente que ocupa     *
c *               a posicao k no vetor au                             *
c *   neq       - numero de equacoes                                  *
c *   ovlp      - overllaping( mpi)  
c * ----------------------------------------------------------------- *
c * parametros de saida                                               *
c * ----------------------------------------------------------------- *
c *  thread_begin(1:nthreads) - primeira linha do sistema da thread i *
c *  thread_end(1:nthreads)   - ultima linha do sistema da thread i   *
c *  thread_heigth(1:nthread) - altura onde a thread escreve no vetor *
c *                             y                                     *
c * ----------------------------------------------------------------- *
c *********************************************************************
      subroutine partition_matrix(ia,ja,neq,ovlp)
      implicit none
      include 'openmp.fi'
      integer ia(*),ja(*),neq,nequ,nad,nnzr
      integer i,idum
      logical ovlp
c
      nnzr = 0 
      if(ovlp) nnzr = ia(neq+2+neq)-1
      call partition_csrc_bynonzeros(ia,ia(neq+2),nnzr,neq,ovlp)
      call compute_effective_work(ia,ja,neq)
c .....................................................................      
c      
c ... checa se a divisao da matriz ocorreu sem problemas
c
      do i = 1 ,nth_solv
        if (thread_height(i)    .le. 0 ) goto 1000
        if ( thread_begin(i)    .le. 0 ) goto 1000
        if ( thread_end(i)      .le. 0 ) goto 1000
      enddo
c .....................................................................      
c
c ...
      return
c ... controle de erro       
1000  continue 
      print*,'***error: divisao da matrix para o openmp falhou!'
      print*,'Diminua o numero de threads usado ou desabilite o openmp.'
      print*,'Log da partition_matrix:'
      do i = 1 ,nth_solv
        print*,'thread id',i,'height' , thread_height(i)
        print*,'thread id',i,'begin ' , thread_begin(i)
        print*,'thread id',i,'end   ' , thread_end(i)
      enddo
      call stop_mef()
c .....................................................................      
c
c ... controle de erro       
1100  continue 
      print*,'***error: divisao da matrix para o openmp falhou!'
      print*,'Diminua o numero de threads usado ou desabilite o openmp.'
      print*,'Log da partition_matrix:'
      do i = 1 ,nth_solv
        print*,'thread id',i,'height  ' , thread_height(i)
        print*,'thread id',i,'begin   ' , thread_begin(i)
        print*,'thread id',i,'end     ' , thread_end(i)
        print*,'thread id',i,'begin pu' , thread_pu_begin(i)
        print*,'thread id',i,'end   pu' , thread_pu_end(i)
      enddo
      call stop_mef()
c .....................................................................      
      end
c *********************************************************************
c
c *********************************************************************
c * PARTITION_CSRC_BYNOZEROS :dividi o trabalho do matvec entre as    *
c * threads considerando valores nao nulos e inicializa a estruturas  *
c * do buffer do matvec                                               *
c * ----------------------------------------------------------------- *
c * parametros de entrada :                                           *
c * ----------------------------------------------------------------- *
c *   ia(neq+1) - ia(i) informa a posicao no vetor au do primeiro     *
c *                     coeficiente nao-nulo da linha   i             *
c *   iar(neq+1)- iar(i) informa a posicao no vetor ar do primeiro    *
c *                     coeficiente nao-nulo da linha i da parte      *
c *                     retangular                                    *
c *   neq       - numero de equacoes                                  *
c *   nnzr      - numero de nad zeros na matriz retangular(mpi)       *
c *   ovlp      - overllaping( mpi)                                   *
c * ----------------------------------------------------------------- *
c * parametros de saida                                               *
c * ----------------------------------------------------------------- *
c *  thread_begin(1:nthreads) - primeira linha do sistema da thread i *
c *  thread_end(1:nthreads)   - ultima linha do sistema da thread i   *
c * ----------------------------------------------------------------- *
c *********************************************************************
      subroutine partition_csrc_bynonzeros(ia,iar,nnzr,neq,ovlp)
      implicit none
      include 'omp_lib.h'
      include 'openmp.fi'
      integer ia(*),iar(*),nnzr,neq
      logical ovlp
      integer mean_variables,line,thread_size(max_num_threads),tam,i
      integer*8 nad  
c
      nad = ia(neq+1)-1
      mean_variables = (2*nad + neq)/nth_solv + 1
      if (ovlp) mean_variables = mean_variables + nnzr/nth_solv
      line = 2
      thread_begin(1) = 1
      do i = 1, nth_solv - 1
        thread_size(i) = 0
        tam = 0
c
 100    tam = 2*(ia(line) - ia(line - 1)) + 1
        if (ovlp) tam = tam + iar(line) - iar(line - 1)
        if (thread_size(i) + tam .le. mean_variables) then
           thread_size(i) = thread_size(i) + tam
           thread_end(i) = line - 1
           thread_begin(i+1) = line
           line = line + 1
           goto 100
        endif
      enddo
      thread_size(nth_solv) = 2*(ia(neq+1) -
     .   ia(thread_begin(nth_solv))) +
     .   neq + 1 - thread_begin(nth_solv)
      if (ovlp) thread_size(nth_solv) = thread_size(nth_solv) +
     .   iar(neq+1) - iar(thread_begin(nth_solv))
      thread_end(nth_solv) = neq
      return
      end
c *********************************************************************
c
c *********************************************************************
c * PARTITION_CSRC_EVENLY : dividi o trabalho entre as threads        *
c * considerando valores nulos e nao nulos  inicializa a estruturas   *
c * do buffer do matvec                                               *
c * ----------------------------------------------------------------- *
c * parametros de entrada :                                           *
c * ----------------------------------------------------------------- *
c *   neq       - numero de equacoes                                  *
c * ----------------------------------------------------------------- *
c * parametros de saida                                               *
c * ----------------------------------------------------------------- *
c *  thread_begin(1:nthreads) - primeira linha do sistema da thread i *
c *  thread_end(1:nthreads)   - ultima linha do sistema da thread i   *
c * ----------------------------------------------------------------- *
c *********************************************************************
      subroutine partition_csrc_evenly(neq)
      implicit none
      include 'omp_lib.h'
      include 'openmp.fi'
      integer neq,split
      integer i,j
c
c$omp parallel private(split)
!$    thread_id = omp_get_thread_num()
      split = mod(neq, nth_solv)
      thread_begin(thread_id+1) = thread_id*(neq/nth_solv)
c
      if (thread_id .lt. split) then
         thread_begin(thread_id+1) = thread_begin(thread_id+1) +
     .      thread_id + 1
         thread_end(thread_id+1) = thread_begin(thread_id+1) +
     .      neq/nth_solv
      else
         thread_begin(thread_id+1) = thread_begin(thread_id+1) +
     .      split + 1
         thread_end(thread_id+1) = thread_begin(thread_id+1) +
     .      neq/nth_solv - 1
      endif
c$omp end parallel
      return
      end
c *********************************************************************
c
c *********************************************************************
c * COMPUTE_EFFECTIVE_WORK: Calculo do trabalho efeitivo por thread   *
c * ----------------------------------------------------------------- *
c * parametros de entrada :                                           *
c * ----------------------------------------------------------------- *
c *   ia(neq+1) - ia(i) informa a posicao no vetor au do primeiro     *
c *                     coeficiente nao-nulo da linha   i             *
c *   ja(nad  ) - ja(k) informa a coluna do coeficiente que ocupa     *
c *               a posicao k no vetor au                             *
c *   neq       - numero de equacoes                                  *
c * ----------------------------------------------------------------- *
c * parametros de saida                                               *
c * ----------------------------------------------------------------- *
c *  thread_heigth(1:nthread) - altura onde a thread escreve no vetor *
c *                             y                                     *
c * ----------------------------------------------------------------- *
c *********************************************************************
      subroutine compute_effective_work(ia,ja,neq)
      implicit none
      include 'omp_lib.h'
      include 'openmp.fi'
      integer ia(*),ja(*),neq,h,i
c
c$omp parallel private(h)
!$    thread_id = omp_get_thread_num()
      h = thread_begin(thread_id+1)
      do i = thread_begin(thread_id+1), thread_end(thread_id+1)
         h = min(h, ja( ia(i) ))
      enddo
      thread_height(thread_id+1) = h
c$omp end parallel
      return
      end
c **********************************************************************
c
c *********************************************************************
c * PARTITION_CSRC_BYNOZEROS_PM :dividi o trabalho do matvec entre as *
c * threads comsiderando valores nao nulos e inicializa a estruturas  *
c * do buffer do matvec para matrizes blocadas Kuu, kpp, kpu.         *
c * ----------------------------------------------------------------- *
c * parametros de entrada :                                           *
c * ----------------------------------------------------------------- *
c *   ia(neq+1) - parte kuu e kpp                                     *
c *   iapu(neq+1)- parte kpu                                          *
c *   neq       - numero de equacoes total                            *
c *   nequ      - numero de equacoes kuu                              *
c * ----------------------------------------------------------------- *
c * parametros de saida                                               *
c * ----------------------------------------------------------------- *
c *  thread_begin(1:nthreads) - primeira linha do sistema da thread i *
c *  thread_end(1:nthreads)   - ultima linha do sistema da thread i   *
c * ----------------------------------------------------------------- *
c *********************************************************************
      subroutine partition_csrc_bynonzeros_pm(ia,iapu,neq,nequ)
      implicit none
      include 'omp_lib.h'
      include 'openmp.fi'
      integer ia(*),iapu(*),nnzr,neq,nequ,neqp
      integer mean_variables,line,thread_size(max_num_threads),tam,i
      integer*8 nad,nadpu  
c
      nad   = ia(neq+1)-1
c ...
      mean_variables = (2*nad + neq)/nth_solv + 1
      line = 2
      thread_begin(1) = 1
      do i = 1, nth_solv - 1
        thread_size(i) = 0
        tam = 0
c
 100    tam = 2*(ia(line) - ia(line - 1) ) + 1
        if (thread_size(i) + tam .le. mean_variables) then
           thread_size(i) = thread_size(i) + tam
           thread_end(i) = line - 1
           thread_begin(i+1) = line
           line = line + 1
           goto 100
        endif
      enddo
      thread_size(nth_solv) = 2*(ia(neq+1) -
     .   ia(thread_begin(nth_solv))) +
     .   neq + 1 - thread_begin(nth_solv)
      thread_end(nth_solv) = neq
c .....................................................................
c
c ...
      neqp  = neq - nequ
      nadpu = iapu(neqp+1)-1
      mean_variables = (2*nadpu + neqp)/nth_solv + 1
      line = 2
      thread_pu_begin(1) = 1
      do i = 1, nth_solv - 1
        thread_size(i) = 0
        tam = 0
c
 110    tam = 2*(iapu(line) - iapu(line - 1) ) + 1
        if (thread_size(i) + tam .le. mean_variables) then
           thread_size(i) = thread_size(i) + tam
           thread_pu_end(i) = line - 1
           thread_pu_begin(i+1) = line
           line = line + 1
           goto 110
        endif
      enddo
      thread_size(nth_solv) = 2*(iapu(neqp+1) -
     .   ia(thread_pu_begin(nth_solv))) +
     .   neqp + 1 - thread_pu_begin(nth_solv)
      thread_pu_end(nth_solv) = neqp
c ...
      do i = 1, nth_solv 
        thread_pm_begin(i)= min(thread_pu_begin(i)+nequ,thread_begin(i))
        thread_pm_end(i)  = max(thread_pu_end(i)+nequ,thread_end(i))
      enddo
    
c .....................................................................
      return
      end
c *********************************************************************
c
c *********************************************************************
c * COMPUTE_EFFECTIVE_WORK: Calculo do trabalho efeitivo por thread   *
c * ----------------------------------------------------------------- *
c * parametros de entrada :                                           *
c * ----------------------------------------------------------------- *
c * ia(neq+1) - matrix kuu e kpp                                      *
c * ja(nad  ) - matrix kuu e kpp                                      *
c *               a posicao k no vetor au                             *
c * iapu(neq+1) - matrix kpu                                          *
c * japu(nadpu) - matriz kpu                                          *
c * nequ        - numero de equacoes                                  *
c * ----------------------------------------------------------------- *
c * parametros de saida                                               *
c * ----------------------------------------------------------------- *
c *  thread_heigth(1:nthread) - altura onde a thread escreve no vetor *
c *                             y                                     *
c * ----------------------------------------------------------------- *
c *********************************************************************
      subroutine compute_effective_work_pm(ia,ja,iapu,japu,nequ)
      implicit none
      include 'omp_lib.h'
      include 'openmp.fi'
      integer ia(*),ja(*),iapu(*),japu(*)
      integer h,i,nequ
c
c$omp parallel private(h)
!$    thread_id = omp_get_thread_num()
c ... kuu e kpp
      h = thread_begin(thread_id+1)
      do i = thread_begin(thread_id+1), thread_end(thread_id+1)
         h = min(h, ja( ia(i) ))
      enddo
      thread_height(thread_id+1) = h
c .....................................................................
c
c ... kpu      
      h = thread_pu_begin(thread_id+1) + nequ
      do i = thread_pu_begin(thread_id+1), thread_pu_end(thread_id+1)
         h = min(h, japu( iapu(i) ))
      enddo
      thread_height(thread_id+1) = min(h,thread_height(thread_id+1))
c$omp end parallel
      return
      end
c **********************************************************************
c
c **********************************************************************
c * GET_BUFFER_SIZE: quantidade de memoria usada no buffer do matvec   *
c * ------------------------------------------------------------------ *
c * Parametro de entrada :                                             *
c * cod - 1 Bytes, 2 KBytes, 3 MBytes, 4 GBytes                        *
c * ------------------------------------------------------------------ *
c * ------------------------------------------------------------------ *
c **********************************************************************
      real*8 function get_buffer_size(tp,neq)
      implicit none
      include 'openmp.fi'
      integer nbytes
      integer neq
      character*2 tp
c ... 8 byts vector
      nbytes = 8*nth_solv*neq
c ......................................................................
c
c ...      
      if(tp .eq.' B') then
        get_buffer_size = nbytes
      else if (tp .eq. 'KB') then
        get_buffer_size = nbytes / 1024
      else if (tp .eq. 'MB') then
        get_buffer_size = nbytes / (1024**2)
      else if (tp .eq. 'GB') then
        get_buffer_size = nbytes / (1024**3)
      endif
c ......................................................................
c
c ...
      return
      end
c ......................................................................
c **********************************************************************
      
