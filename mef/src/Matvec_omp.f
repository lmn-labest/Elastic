c **********************************************************************
c *                                                                    *
c *   MATVEC.F                                           14/09/2011    *
c *                                                                    *
c *   Este arquivo contem subrotinas para produto matriz-vetor e       *
c *   produto escalar de vetores                                       *
c *                                                                    *
c *   matvec_csrc_omp                                                  *
c *   matvec_csrsym_omp                                                *
c *   matvec_csrcr_omp                                                 *
c *   matvec_csrcrsym_omp                                              *
c *   dot_par_omp                                                      *
c *   dot_par_omp_loopwise                                             *
c *                                                                    *
c *   1 = loop interno desenrolado                                     *
c *   2 = loops desenrolados                                           *
c *                                                                    *
c **********************************************************************
c
c **********************************************************************
c *                                                                    *
c *   MATVEC_CSRC_OMP: produto matriz-vetor y = Ax  (A nao-simetrica), *
c *                    coef. de A no formato CSRC e grafo simetrico.   *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *                                                                    *
c *   neq   - numero de equacoes                                       *
c *   ia(neq+1) - ia(i) informa a posicao no vetor au do primeiro      *
c *                     coeficiente nao-nulo da linha   i              *
c *   ja(neq+1) - ja(k) informa a coluna do coeficiente que ocupa      *
c *               a posicao k no vetor au                              *
c *   ad(neq)- diagonal da matriz A                                    *
c *   al(nad)- parte triangular inferior de A, no formato CSR, ou      *
c *            parte triangular superior de A, no formato CSC          *
c *   au(*)  - nao utilizado                                           *
c *   x(neq) - vetor a ser multiplicado                                *
c *   y(neq) - nao definido                                            *
c *   neqf1i - numero de equacoes no buffer de recebimento (MPI)       *
c *   neqf2i - numero de equacoes no buffer de envio (MPI)             *
c *   i_fmapi- ponteiro para o mapa de comunicacao  (MPI)              *
c *   i_xfi  - ponteiro para o buffer de valores    (MPI)              *
c *   i_rcvsi- ponteiro extrutura da comunicacao    (MPI)              *
c *   i_dspli- ponteiro extrutura da comunicacao    (MPI)              *
c *   thread_y - buffer de equacoes para o vetor y (openmp)            *
c *                                                                    *
c *   Parametros de saida:                                             *
c *                                                                    *
c *   y(neq) - vetor contendo o resultado do produto y = Ax            *
c *                                                                    *
c **********************************************************************
      subroutine matvec_csrc_omp(neq,ia,ja,dum0,dum1,ad,al,au,dum2,
     .                           x,y,neqf1i,neqf2i,i_fmapi,i_xfi,
     .                           i_rcvsi,i_dspli,thread_y)
      implicit none
      include 'mpif.h'
      include 'omp_lib.h'
      include 'parallel.fi'
      include 'time.fi'
      include 'openmp.fi'
      integer neq,ia(*),ja(*),dum0,dum1,dum2,i,j,k,jak,inc
      real*8  ad(*),al(*),au(*),x(*),y(*),t,xi
      integer neqf1i,neqf2i
c ... ponteiros      
      integer*8 i_fmapi,i_xfi,i_rcvsi,i_dspli
      real*8  thread_y(*)
c ......................................................................
c$omp single
      time0 = MPI_Wtime()
c$omp end single
!$    thread_id = omp_get_thread_num() + 1
      do i = 1, nth_solv
         inc = (i-1)*neq
c$omp do
         do j = thread_height(i)+inc, thread_begin(i)+inc-1
            thread_y(j) = 0.d0
         enddo
c$omp end do
      enddo
      inc = (thread_id - 1)*neq
c$omp barrier
      do 110 i = thread_begin(thread_id), thread_end(thread_id)
         y(i) = 0.d0
         xi = x(i)
         t  = ad(i)*xi
         do 100 k = ia(i), ia(i+1)-1
            jak = ja(k)
            t   = t + al(k)*x(jak)
            jak = jak + inc
            thread_y(jak) = thread_y(jak) + au(k)*xi
  100    continue
         thread_y(i+inc) = t
  110 continue
c$omp barrier
c
c ... Acumula thread_y(i) em y(i)
c
      do i = 1, nth_solv
         inc = (i-1)*neq
c$omp do
         do j = thread_height(i), thread_end(i)
            y(j) = y(j) + thread_y(j+inc)
         end do
c$omp end do
      end do
c$omp single
      matvectime = matvectime + MPI_Wtime() - time0
c ......................................................................
c
c ... Comunicacao do vetor y no sistema non-overlapping:
c
      if (novlp) call communicate(y,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,
     .                            i_dspli)
c$omp end single
c ......................................................................
      return
      end
c *********************************************************************
c
c **********************************************************************
c *                                                                    *
c *   MATVEC_CSRCSYM_OMP: produto matriz-vetor y = Ax  (A simetrica),  *
c *                    coef. de A no formato CSRC e grafo simetrico.   *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *                                                                    *
c *   neq   - numero de equacoes                                       *
c *   ia(neq+1) - ia(i) informa a posicao no vetor au do primeiro      *
c *                     coeficiente nao-nulo da linha   i              *
c *   ja(neq+1) - ja(k) informa a coluna do coeficiente que ocupa      *
c *               a posicao k no vetor au                              *
c *   ad(neq)- diagonal da matriz A                                    *
c *   al(nad)- parte triangular inferior de A, no formato CSR, ou      *
c *            parte triangular superior de A, no formato CSC          *
c *   au(*)  - nao utilizado                                           *
c *   x(neq) - vetor a ser multiplicado                                *
c *   y(neq) - nao definido                                            *
c *   neqf1i - numero de equacoes no buffer de recebimento (MPI)       *
c *   neqf2i - numero de equacoes no buffer de envio (MPI)             *
c *   i_fmapi- ponteiro para o mapa de comunicacao  (MPI)              *
c *   i_xfi  - ponteiro para o buffer de valores    (MPI)              *
c *   i_rcvsi- ponteiro extrutura da comunicacao    (MPI)              *
c *   i_dspli- ponteiro extrutura da comunicacao    (MPI)              *
c *   thread_y - buffer de equacoes para o vetor y (openmp)            *
c *                                                                    *
c *   Parametros de saida:                                             *
c *                                                                    *
c *   y(neq) - vetor contendo o resultado do produto y = Ax            *
c *                                                                    *
c **********************************************************************
      subroutine matvec_csrcsym_omp(neq,ia,ja,dum0,dum1,ad,al,dum2,dum3,
     .                              x,y,neqf1i,neqf2i,i_fmapi,i_xfi,
     .                              i_rcvsi,i_dspli,thread_y)
      implicit none
      include 'mpif.h'
      include 'omp_lib.h'
      include 'parallel.fi'
      include 'time.fi'
      include 'openmp.fi'
      integer neq,ia(*),ja(*),dum0,dum1,dum2,dum3,i,j,k,jak,inc
      real*8  ad(*),al(*),x(*),y(*),s,t,xi
      integer neqf1i,neqf2i
c ... ponteiros      
      integer*8 i_fmapi,i_xfi,i_rcvsi,i_dspli
      real*8  thread_y(*)
c ......................................................................
c$omp single
      time0 = MPI_Wtime()
c$omp end single
!$    thread_id = omp_get_thread_num() + 1
      do i = 1, nth_solv
         inc = (i-1)*neq
c$omp do
         do j = thread_height(i)+inc, thread_begin(i)+inc-1
            thread_y(j) = 0.d0
         enddo
c$omp end do
      enddo
      inc = (thread_id - 1)*neq
c$omp barrier
      do 110 i = thread_begin(thread_id), thread_end(thread_id)
         y(i) = 0.d0
         xi = x(i)
         t  = ad(i)*xi
         do 100 k = ia(i), ia(i+1)-1
            jak = ja(k)
            s   = al(k)
            t   = t + s*x(jak)
            jak = jak + inc
            thread_y(jak) = thread_y(jak) + s*xi
  100    continue
         thread_y(i+inc) = t
  110 continue
c$omp barrier
c
c ... Accumulate thread_y(i) into y(i)
c
      do i = 1, nth_solv
         inc = (i-1)*neq
c$omp do
         do j = thread_height(i), thread_end(i)
            y(j) = y(j) + thread_y(j+inc)
         end do
c$omp end do
      end do
c$omp single
      matvectime = matvectime + MPI_Wtime() - time0
      if (novlp) call communicate(y,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,
     .                            i_dspli)
c$omp end single
c ......................................................................
      return
      end
c **********************************************************************
c 
c **********************************************************************
c *                                                                    *
c *   MATVEC_CSRCR_OMP: produto matriz-vetor y = Ax  (A nao-simetrica),*
c *                    coef. de A no formato CSRCR.                    *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *                                                                    *
c *   neq   - numero de equacoes                                       *
c *   ia(neq+1) - ia(i) informa a posicao no vetor au do primeiro      *
c *                     coeficiente nao-nulo da linha   i              *
c *   ja(neq+1) - ja(k) informa a coluna do coeficiente que ocupa      *
c *               a posicao k no vetor au                              *
c *   ad(neq)- diagonal da matriz A                                    *
c *   al(nad)- parte triangular inferior de A, no formato CSR, ou      *
c *            parte triangular superior de A, no formato CSC          *
c *   au(*)  - nao utilizado                                           *
c *   x(neq) - vetor a ser multiplicado                                *
c *   y(neq) - nao definido                                            *
c *   neqf1i - numero de equacoes no buffer de recebimento (MPI)       *
c *   neqf2i - numero de equacoes no buffer de envio (MPI)             *
c *   i_fmapi- ponteiro para o mapa de comunicacao  (MPI)              *
c *   i_xfi  - ponteiro para o buffer de valores    (MPI)              *
c *   i_rcvsi- ponteiro extrutura da comunicacao    (MPI)              *
c *   i_dspli- ponteiro extrutura da comunicacao    (MPI)              *
c *   thread_y - buffer de equacoes para o vetor y (openmp)            *
c *                                                                    *
c *   Parametros de saida:                                             *
c *                                                                    *
c *   y(neq) - vetor contendo o resultado do produto y = Ax            *
c *                                                                    *
c **********************************************************************
      subroutine matvec_csrcr_omp(neq,ia,ja,ia1,ja1,ad,al,au,ar,x,y,
     .                            neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,
     .                            i_dspli,thread_y)
      implicit none
      include 'mpif.h'
      include 'omp_lib.h'
      include 'parallel.fi'
      include 'time.fi'
      include 'openmp.fi'
      integer neq,ia(*),ja(*),ia1(*),ja1(*),i,j,k,jak,inc
      real*8  ad(*),al(*),au(*),ar(*),x(*),y(*),t,xi
      integer neqf1i,neqf2i
c ... ponteiros      
      integer*8 i_fmapi,i_xfi,i_rcvsi,i_dspli
      real*8  thread_y(*)
c$omp single
c ......................................................................
c
c ... Comunicacao do vetor x no sistema overlapping:
      call communicate(x,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli)
c ......................................................................
      time0 = MPI_Wtime()
c$omp end single
!$    thread_id = omp_get_thread_num() + 1
      do i = 1, nth_solv
         inc = (i-1)*neq
c$omp do
         do j = thread_height(i)+inc, thread_begin(i)+inc-1
            thread_y(j) = 0.d0
         enddo
c$omp end do
      enddo
      inc = (thread_id - 1)*neq
c$omp barrier
      do 200 i = thread_begin(thread_id), thread_end(thread_id)
         y(i) = 0.d0
         xi = x(i)
         t  = ad(i)*xi
         do 100 k = ia(i), ia(i+1)-1
            jak = ja(k)
            t   = t + al(k)*x(jak)
            jak = jak + inc
            thread_y(jak) = thread_y(jak) + au(k)*xi
  100    continue
c 
         do 110 k = ia1(i), ia1(i+1)-1
            jak = ja1(k)
            t   = t + ar(k)*x(jak)
  110    continue
         thread_y(i+inc) = t         
  200 continue
c$omp barrier
c
c ... Accumulate thread_y(i) into y(i)
c
      do i = 1, nth_solv
         inc = (i-1)*neq
c$omp do
         do j = thread_height(i), thread_end(i)
            y(j) = y(j) + thread_y(j+inc)
         end do
c$omp end do
      end do
c$omp single
c ......................................................................
      matvectime = matvectime + MPI_Wtime() - time0
c$omp end single
      return
      end
c **********************************************************************
c 
c **********************************************************************
c *                                                                    *
c *   MATVEC_CSRCRSYM_OMP: produto matriz-vetor y = Ax                 *
c *                        (A simetrica) e grafo simetrico             *
c *                        coef. de A no formato CSRCR.                *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *                                                                    *
c *   neq   - numero de equacoes                                       *
c *   ia(neq+1) - ia(i) informa a posicao no vetor au do primeiro      *
c *                     coeficiente nao-nulo da linha   i              *
c *   ja(neq+1) - ja(k) informa a coluna do coeficiente que ocupa      *
c *               a posicao k no vetor au                              *
c *   ad(neq)- diagonal da matriz A                                    *
c *   al(nad)- parte triangular inferior de A, no formato CSR, ou      *
c *            parte triangular superior de A, no formato CSC          *
c *   au(*)  - nao utilizado                                           *
c *   x(neq) - vetor a ser multiplicado                                *
c *   y(neq) - nao definido                                            *
c *   neqf1i - numero de equacoes no buffer de recebimento (MPI)       *
c *   neqf2i - numero de equacoes no buffer de envio (MPI)             *
c *   i_fmapi- ponteiro para o mapa de comunicacao  (MPI)              *
c *   i_xfi  - ponteiro para o buffer de valores    (MPI)              *
c *   i_rcvsi- ponteiro extrutura da comunicacao    (MPI)              *
c *   i_dspli- ponteiro extrutura da comunicacao    (MPI)              *
c *   thread_y - buffer de equacoes para o vetor y (openmp)            *
c *                                                                    *
c *   Parametros de saida:                                             *
c *                                                                    *
c *   y(neq) - vetor contendo o resultado do produto y = Ax            *
c *                                                                    *
c **********************************************************************
      subroutine matvec_csrcrsym_omp(neq,ia,ja,ia1,ja1,ad,al,dum0,ar,x,
     .           y,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli,thread_y)
      implicit none
      include 'mpif.h'
      include 'omp_lib.h'
      include 'parallel.fi'
      include 'time.fi'
      include 'openmp.fi'
      integer neq,ia(*),ja(*),ia1(*),ja1(*),i,j,k,jak,dum0,inc
      real*8  ad(*),al(*),ar(*),x(*),y(*),t,xi,s,tm0
      real*8 thread_y(*)
      integer neqf1i,neqf2i
c ... ponteiros      
      integer*8 i_fmapi,i_xfi,i_rcvsi,i_dspli
c$omp single
c ......................................................................
c
c ... Comunicacao do vetor x no sistema overlapping:
c
      call communicate(x,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli)
c ......................................................................
      time0 = MPI_Wtime()
c$omp end single
!$    thread_id = omp_get_thread_num() + 1
      do i = 1, nth_solv
         inc = (i-1)*neq
c$omp do
         do j = thread_height(i)+inc, thread_begin(i)+inc-1
            thread_y(j) = 0.d0
         enddo
c$omp end do
      enddo
      inc = (thread_id - 1)*neq
c$omp barrier
      do 200 i = thread_begin(thread_id), thread_end(thread_id)
         y(i) = 0.d0
         xi = x(i)
         t  = ad(i)*xi
         do 100 k = ia(i), ia(i+1)-1
            jak = ja(k)
            s = al(k)
            t   = t + s*x(jak)
            jak = jak + inc
            thread_y(jak) = thread_y(jak) + s*xi
  100    continue
c
         do 110 k = ia1(i), ia1(i+1)-1
            jak = ja1(k)
            t   = t + ar(k)*x(jak)
  110    continue
         thread_y(i+inc) = t
  200 continue
c$omp barrier
c
c ... Accumulate thread_y(i) into y(i)
c
      do i = 1, nth_solv
         inc = (i-1)*neq
c$omp do
         do j = thread_height(i), thread_end(i)
            y(j) = y(j) + thread_y(j+inc)
         end do
c$omp end do
      end do
c$omp single
      matvectime = matvectime + MPI_Wtime() - time0
c$omp end single
c ......................................................................
      return
      end
c *********************************************************************
c
c **********************************************************************
c *                                                                    *
c *   DOT_PAR_OMP: produto escalar omp(versao de diretiva orfã,  funcao*
c *   chamada de uma regiao paralela, portanto nao ha a necessidade de *
c *   abertura de uma nova regia paralela).                            *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *                                                                    *
c *   a        - vetor                                                 *
c *   b        - vetor                                                 *
c *   neq_doti - dimensao dos vetores                                  *
c *                                                                    *
c *   Parametros de saida:                                             *
c *                                                                    *
c *   dot    - valor do produto escalar                                *
c *                                                                    *
c **********************************************************************
      real*8 function dot_par_omp(a,b,neq_doti)
      implicit none
      include 'mpif.h'
      include 'omp_lib.h'
      include 'parallel.fi'
      include 'time.fi'
      include 'openmp.fi'
      integer i,neq_doti
      real*8  a(*),b(*),tmp
c ......................................................................
c$omp single
      time0 = MPI_Wtime()
      omp_dot = 0.d0
c$omp end single
c$omp do reduction(+:omp_dot)
      do 100 i = 1, neq_doti
         omp_dot = omp_dot + a(i)*b(i)
  100 continue
c$omp end do
c$omp single
      dottime = dottime + MPI_Wtime() - time0
c ......................................................................
      if (nprcs .gt. 1) then
         call MPI_ALLREDUCE(omp_dot,tmp,1,
     .                      MPI_DOUBLE_PRECISION,
     .                      MPI_SUM,MPI_COMM_WORLD,ierr)
         omp_dot = tmp
      endif
c$omp end single
      dot_par_omp = omp_dot
c$omp barrier 
c ......................................................................
      return
      end
c *********************************************************************
c
c **********************************************************************
c *                                                                    *
c *   DOT_PAR_OMP_LLOPWISE: produto escalar omp(funcao chamada de uma  *
c *   regiao sequancial, portanto a necessidade da abertura de uma     *
c *   regiao paralela).                                                *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *                                                                    *
c *   a        - vetor                                                 *
c *   b        - vetor                                                 *
c *   neq_doti - dimensao dos vetores                                  *
c *                                                                    *
c *   Parametros de saida:                                             *
c *                                                                    *
c *   dot    - valor do produto escalar                                *
c *                                                                    *
c **********************************************************************
      real*8 function dot_par_omp_loopwise(a,b,neq_doti)
      implicit none
      include 'mpif.h'            
      include 'parallel.fi'
      include 'time.fi'
      integer neq_doti 
      integer i
      real*8  a(*),b(*),tmp
c ......................................................................
      time0 = MPI_Wtime()
      tmp = 0.d0
c$omp parallel do reduction(+:tmp)
      do 100 i = 1, neq_doti
         tmp = tmp + a(i)*b(i)
  100 continue
c$omp end parallel do
      dottime = dottime + MPI_Wtime() - time0
c ......................................................................
      if (nprcs .gt. 1) then
         call MPI_ALLREDUCE(tmp,dot_par_omp_loopwise,1,
     .                 MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      else
         dot_par_omp_loopwise = tmp
      endif
c ......................................................................    
      return
      end
