c *********************************************************************
c * Metodos iterativos para solucao de sistemas lineares (OpenMP)     *
c * ----------------------------------------------------------------- *
c * simetricos:                                                       *
c * ----------------------------------------------------------------- *
c * PCG - gradiente conjugados com precondicionador diagonal          *
c * ----------------------------------------------------------------- *
c * nao-simetricos:                                                   *
c * ----------------------------------------------------------------- *                                                       *
c * pbicgstab - gradiente bi-conjugados estabilizados  com            * 
c * precondicionador diagonal                                         *
c *                                                                   *
c * gmres(m) - GMRES com precondicionador diagonal                    *
c *                                                                   *
c * ----------------------------------------------------------------- *
c ********************************************************************* 
      subroutine pcg_omp(neq   ,nad ,ia ,ja
     .                  ,ad    ,au  ,al ,m  ,b    ,x     
     .                  ,z     ,r   ,p  ,tol,maxit
     .                  ,matvec,dot
     .                  ,my_id ,neqf1i ,neqf2i ,neq_doti,i_fmapi
     .                  ,i_xfi ,i_rcvsi,i_dspli,thread_y
     .                  ,fprint,flog   ,fnew)
c **********************************************************************
c * Data de criacao    : 00/00/0000                                    *
c * Data de modificaco : 25/04/2016                                    * 
c * ------------------------------------------------------------------ *   
c * PCG_OMP : Solucao de sistemas de equacoes pelo metodo dos          *
c * gradientes conjugados com precondicionador diagonal para matrizes  *
c * simetricas                                                         *
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ * 
c * neq      - numero de equacoes                                      *
c * nequ     - numero de equacoes no bloco Kuu                         *
c * nad      - numero de termos nao nulos no bloco Kuu e Kpu  ou K     *
c * ia(*)    - ponteiro do formato CSR                                 *
c * ja(*)    - ponteiro das colunas no formato CSR                     *
c * ad(neq)  - diagonal da matriz A                                    *
c * au(*)    - parte triangular superior de A                          *
c * al(*)    - parte triangular inferior de A                          *
c * m(*)     - precondicionador diagonal                               *
c * b(neq)   - vetor de forcas                                         *
c * x(neq)   - chute inicial                                           *
c * z(neq)   - arranjo local de trabalho                               *
c * r(neq)   - arranjo local de trabalho                               *
c * p(neq)   - arranjo local de trabalho                               *
c * tol      - tolerancia de convergencia                              *
c * maxit    - numero maximo de iteracoes                              *
c * matvec   - nome da funcao externa para o produto matrix-vetor      *
c * dot      - nome da funcao externa para o produto escalar           *
c * my_id    - id do processo(MPI)                                     *
c * neqf1i   - numero de equacoes no buffer de recebimento (MPI)       *
c * neqf2i   - numero de equacoes no buffer de envio (MPI)             *
c * neq_doti - numero de equacoes no produto interno (MPI)             *
c * i_fmapi  - ponteiro para o mapa de comunicacao  (MPI)              *
c * i_xfi    - ponteiro para o buffer de valores    (MPI)              *
c * i_rcvsi  - ponteiro extrutura da comunicacao    (MPI)              *
c * i_dspli  - ponteiro extrutura da comunicacao    (MPI)              *
c * thread_y - buffer de equacoes para o vetor y (openmp)              *
c * fprint   - saida na tela                                           *
c * flog     - log do arquivo de saida                                 *
c * fnew     - .true.  -> x0 igual a zero                              *
c *            .false. -> x0 dado                                      *
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ *
c * x(neq) - vetor solucao                                             *
c * b(neq) - modificado                                                *
c * ad(*),al(*),au(*) - inalterados                                    *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ *
c ********************************************************************** 
      implicit none
      include 'mpif.h'
      include 'omp_lib.h'
      include 'openmp.fi'
c ... Mpi      
      integer ierr
c .....................................................................      
      integer neqf1i,neqf2i,neq_doti
c ... ponteiros      
      integer*8 i_fmapi,i_xfi
      integer*8 i_rcvsi,i_dspli
c ......................................................................      
      integer neq,maxit,i,j,jj,nad
      integer ia(*),ja(*),my_id
      real*8  ad(*),au(*),al(*),b(*),m(*),x(*)
      real*8  r(*),z(*),p(*)
      real*8  dot,tol,conv,xkx,norm,d,di,alpha,beta,tmp
      real*8  time0,time
      real*8  thread_y(*)
      logical flog,fprint,fnew
      external matvec,dot
c ======================================================================
      time0 = MPI_Wtime()
c ......................................................................
c
c ...
      do 5 i = 1, neq
        if(ad(i) .eq. 0.d0 ) then
          write(*,1000) i
          call stop_mef()
        endif 
   5  continue
c ......................................................................
c$omp parallel default(none) 
c$omp.private(i,j,jj,xkx,norm,d,di,conv,alpha,beta,tmp)
c$omp.shared(neq,nad,ia,ja,al,ad,au,b,x,m,z,r,p,tol,maxit,thread_y)
c$omp.shared(neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli,neq_doti)
c$omp.shared(flog,fprint,fnew,my_id,time,time0)
c$omp.num_threads(nth_solv)                                          
c ......................................................................
c
c ... Chute inicial:
c
      if(fnew) then 
c$omp do
        do 10 i = 1, neq
           x(i) = 0.d0
   10   continue
c$omp end do
      endif
c ......................................................................
c
c ... conv = tol * |(M-1)b|
c$omp do
      do 15 i = 1, neq
         z(i) = b(i) * m(i)
   15 continue
c$omp end do
      d    = dot(z,z,neq_doti)
      conv = tol*dsqrt(dabs(d))
c .......................................................................
c
c ... Ax0
      call matvec(neq,ia,ja,ia(neq+2),ja(nad+1),ad,al,au,al(nad+1) 
     .           ,x,z
     .           ,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli,thread_y)
c .......................................................................
c
c ...
c$omp do
      do 100 i = 1, neq
c ... r0 = b - Ax0
         r(i) = b(i) - z(i)
c ... z0 = (M-1)r0
         z(i) = r(i) * m(i)
c ... p0 = r0
         p(i) = z(i)
  100 continue
c$omp end do
c ... ( r(0),z(0) ) = ( r(0), (M-1)r0 )
      d    = dot(r,z,neq_doti)
c ......................................................................
      jj = 1
      do 230 j = 1, maxit
c ... z = Ap(j)
         call matvec(neq,ia,ja,ia(neq+2),ja(nad+1),ad,al,au,al(nad+1) 
     .              ,p,z
     .              ,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli
     .              ,thread_y)
c .....................................................................
c
c ... alpha = ( r(j),z(j) ) / ( Ap(j), p(j) ))
         alpha = d / dot(z,p,neq_doti)
c .....................................................................
c
c ...
c$omp do
         do 210 i = 1, neq
c ... x(j+1) = x(j) + alpha*p
            x(i) = x(i) + alpha * p(i)
c ... r(j+1) = r(j) - alpha*Ap
            r(i) = r(i) - alpha * z(i)
c ... z  = (M-1)r0
            z(i) = r(i) * m(i)
  210    continue
c$omp end do
c .....................................................................
c
c ... ( r(j+1),(M-1)r(j+1) ) = ( r(j+1),z )
         di   = dot(r,z,neq_doti) 
c ... beta = ( r(j+1),(M-1)r(j+1) ) / ( r(j),r(j) ) 
         beta = di / d
c .....................................................................
c
c ...
c$omp do
         do 220 i = 1, neq
c ... p(j+1) = (M-1)r(j+1) + beta*p(j) = z + beta*p(j)
            p(i) = z(i) + beta * p(i)
  220    continue
c$omp end do
c .....................................................................
c
c ...
         d =  di
         if (dsqrt(dabs(d)) .lt. conv) goto 300
c ......................................................................
c$omp master
         if( jj .eq.500) then
           jj = 0
           write(*,1300),j,dsqrt(dabs(d)),conv 
         endif  
         jj = jj + 1
c$omp end master
c ......................................................................
  230 continue
c ----------------------------------------------------------------------
c$omp single
      write(*,1200) maxit
      if(flog) write(10,1200) maxit
      call stop_mef()
c$omp end single
  300 continue
c
c ... norm:  x*Kx
      call matvec(neq,ia,ja,ia(neq+2),ja(nad+1),ad,al,au,al(nad+1) 
     .           ,x,z
     .           ,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli
     .           ,thread_y)
      xkx = dot(x,z,neq_doti)
c ......................................................................
c
c ... norm-2 = || x ||
      norm = dsqrt(dot(x,x,neq_doti))
c ......................................................................
c
c ... r = b - Ax (calculo do residuo explicito)
c$omp do
      do 310 i = 1, neq
        r(i) = b(i) - z(i)
  310 continue
c$omp end do
      tmp  = dot(r,r,neq_doti)
      tmp = dsqrt(tmp)
c$omp single
      if( tmp .gt. 3.16d0*conv ) then
         if(my_id .eq.0 )then
           write(*,1400) tmp,conv
         endif 
      endif
c$omp end single
c ......................................................................
c$omp single
c ......................................................................
      time = MPI_Wtime()
      time = time-time0
c ......................................................................
      if(my_id .eq.0 .and. fprint )then
        write(*,1100)tol,conv,neq,nad,j,xkx,norm,time
      endif
c ......................................................................
c
c ... Controle de flops
      if(flog) then
        if(my_id.eq.0) then
          write(10,'(a,a,i9,a,d20.10,a,d20.10,a,d20.10,a,f20.2)')
     .       'PCG_OMP: ',' it ',j, ' x * Kx ',xkx,' ||x|| ',norm
     .      ,' tol ',tol,' time ',time
        endif
      endif
c ......................................................................
c$omp end single
c$omp end parallel
c ......................................................................
      return
c ======================================================================
 1000 format (//,5x,'SUBROTINA PCG:',/,5x,'Coeficiente da diagonal nulo'
     . '- equacao ',i9)
 1100 format(' (PCG_OMP) solver:'/
     . 5x,'Solver tol           = ',d20.6/
     . 5x,'tol * ||b||m         = ',d20.6/
     . 5x,'Number of equations  = ',i20/
     . 5x,'nad                  = ',i20/
     . 5x,'Number of iterations = ',i20/
     . 5x,'x * Kx               = ',d20.10/
     . 5x,'|| x ||              = ',d20.10/
     . 5x,'CPU time (s)         = ',f20.2/)
 1200 format (' *** WARNING: No convergence reached after ',i9,
     .        ' iterations !',/)
 1300 format (' PCG_OMP:',5x,'It',i7,5x,2d20.10)
1400  format (' PCG:',1x,'Residuo exato > 3.16d0*conv '
     .       ,1x,d20.10,1x,d20.10)
      end
c **********************************************************************
      subroutine gmres_omp(neq    ,nequ ,nad,ia ,ja
     .                    ,ad     ,au   ,al ,m  ,b
     .                    ,x      ,k    ,g  ,h  ,y
     .                    ,c      ,s    ,e  ,tol,maxit
     .                    ,matvec ,dot
     .                    ,neqovlp
     .                    ,my_id  ,neqf1i,neqf2i,neq_doti,i_fmapi
     .                    ,i_xfi  ,i_rcvsi,i_dspli,thread_y,flog)
c **********************************************************************
c * Data de criacao    : 00/00/0000                                    *
c * Data de modificaco : 12/12/2015                                    * 
c * ------------------------------------------------------------------ *   
c * GMRES_OMP:Solucao iterativa de sistemas simetricos e nao-simetricos*
c * pelo metodo GMRES com precondicionador diagonal.                   *
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ *
c * neq      - numero de equacoes                                      *
c * nequ     - numero de equacoes no bloco Kuu                         *
c * nad      - numero de termos nao nulos no bloco Kuu e Kpu  ou K     *
c * ia(*)    - ponteiro do formato CSR                                 *
c * ja(*)    - ponteiro das colunas no formato CSR                     *
c * ad(neq)  - diagonal da matriz A                                    *
c * au(*)    - parte triangular superior de A                          *
c * al(*)    - parte triangular inferior de A                          *
c * m(*)     - precondicionador diagonal                               *
c * b(neq)   - vetor de forcas                                         *
c * x(neq)   - chute inicial                                           *
c * k        - base de Krylov                                          *
c * z(neq)   - arranjo local de trabalho                               *
c * r(neq)   - arranjo local de trabalho                               *
c * tol      - tolerancia de convergencia                              *
c * maxit    - numero maximo de iteracoes                              *
c * matvec   - nome da funcao externa para o produto matrix-vetor      *
c * dot      - nome da funcao externa para o produto escalar           *
c * my_id    -                                                         *
c * neqf1i   -                                                         *
c * neqf2i   -                                                         *
c * neq_doti -                                                         *
c * i_fmap   -                                                         *
c * i_xfi    -                                                         *
c * i_rvcs   -                                                         *
c * i_dspli  -                                                         *
c * thread_y - buffer de equacoes para o vetor y (openmp)              *
c * flog     - log do arquivo de saida                                 *
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ *
c * x(neq) - vetor solucao                                             *
c * b(neq),ad(*),al(*),au(*) - modificados                             *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ *
c * Arranjos locais de trabalho:                                       *
c *                                                                    *
c * g(neq+1,k+1)                                                       *
c * h(k+1,k)                                                           *
c * y(k)                                                               *
c * c(k)                                                               *
c * s(k)                                                               *
c * e(k+1)                                                             *
c **********************************************************************
      implicit none
      include 'mpif.h'
      include 'omp_lib.h'
      include 'openmp.fi'
      integer neqf1i,neqf2i,neq_doti
c ... ponteiros      
      integer*8 i_fmapi,i_xfi
      integer*8 i_rcvsi,i_dspli
c .....................................................................      
      integer neq,nequ,k,maxit,ia(*),ja(*)
      integer neqovlp,nit,i,j,jj,l,ni,ic,nad,nad1
      real*8  ad(neq),au(*),al(*),m(*),b(*),x(*)
      real*8  g(neqovlp,1:k+1),h(k+1,k),y(k),c(k),s(k),e(k+1),tol
      real*8  energy,econv,norm,dot,ddot,r,aux1,aux2,beta
      real*8  time0,time
      real*8  thread_y(*)
      logical flog
      external matvec,dot
      integer my_id
c ......................................................................
      time0 = MPI_Wtime()
c ......................................................................
c
c.... Chute inicial:
c
c$omp parallel do num_threads(nth_solv)
      do 10 i = 1, neq
        x(i) = 0.d0
c ...    pre-condicionador diagonal:                           
        g(i,1) = b(i) * m(i)
   10 continue
c$omp end parallel do
c ----------------------------------------------------------------------
c
c ... Limite de convergencia:
c
c$omp parallel private(aux1) num_threads(nth_solv) 
      aux1  = dot(g(1,1),g(1,1),neq_doti)
c$omp single
      econv = tol*dsqrt(aux1)
c$omp end single
c$omp end parallel
c ----------------------------------------------------------------------
c
c ... Ciclos GMRES:
c
      nit = 0
      jj  = 0
      do 1000 l = 1, maxit
c
c ...... Residuo g(1) = b - A x:
c
c$omp parallel num_threads(nth_solv) 
         call matvec(neq,ia,ja,ia(neq+2),ja(nad+1),ad,al,au,al(nad+1)
     .              ,x,g(1,1)
     .              ,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli
     .              ,thread_y)
c$omp end parallel
c
c ...... Residuo com precondicionador diagonal:
c
c$omp parallel do num_threads(nth_solv) 
         do 200 i = 1, neq
            g(i,1) = (b(i) - g(i,1))*m(i)
  200    continue
c$omp end parallel do
c
c ...... Norma do residuo:
c
c$omp parallel private(aux1) num_threads(nth_solv) 
         aux1 = dot(g(1,1),g(1,1),neq_doti)
c$omp single
         e(1) = dsqrt(aux1)
c$omp end single
c$omp end parallel
c
c ...... Normalizacao de g1:
c
c$omp parallel do num_threads(nth_solv) 
         do 210 i = 1, neq
            g(i,1) = g(i,1)/e(1)
  210    continue
c$omp end parallel do
c
c ...... Iteracoes GMRES:
c
         ni = 0
         do 400 i = 1, k
            nit = nit + 1
            ni  = ni  + 1
c
c ......... Produto g(i+1) = A.g(i):
c
c$omp parallel num_threads(nth_solv) 
            call matvec(neq,ia,ja,ia(neq+2),ja(nad+1),ad,al,au,al(nad+1)
     .                 ,g(1,i) ,g(1,i+1)
     .                 ,neqf1i ,neqf2i
     .                 ,i_fmapi,i_xfi,i_rcvsi,i_dspli,thread_y)
c$omp end parallel
c
c ......... Precondicionador diagonal:
c
c$omp parallel do num_threads(nth_solv) 
            do 300 j = 1, neq
               g(j,i+1) = g(j,i+1)*m(j)
  300       continue
c$omp end parallel do
c
c ......... Ortogonalizacao (Gram-Schmidt modificado):
c
            do 320 j = 1, i
c$omp parallel private(aux1) num_threads(nth_solv) 
               aux1 = dot(g(1,i+1),g(1,j),neq_doti)
c$omp single
               beta = aux1
c$omp end single
c$omp do 
               do 310 ic = 1, neq
                  g(ic,i+1) = g(ic,i+1) - beta * g(ic,j)
  310          continue
c$omp end do
c$omp end parallel
               h(j,i) = beta
  320       continue
c
c ......... Norma de g(i+1):
c
c$omp parallel private(aux1) num_threads(nth_solv) 
            aux1 = dot(g(1,i+1),g(1,i+1),neq_doti)
c$omp single
            norm = dsqrt(aux1)
c$omp end single
c$omp end parallel
c
            h(i+1,i) = norm
c
c ......... Normalizacao de g(i+1):
c
c$omp parallel do num_threads(nth_solv) 
            do 330 ic = 1, neq
               g(ic,i+1) = g(ic,i+1)/norm
  330       continue
c$omp end parallel do
c
            do 340 j = 1, i-1
               aux1 =  c(j) * h(j,i) + s(j) * h(j+1,i)
               aux2 = -s(j) * h(j,i) + c(j) * h(j+1,i)
               h(j,i)   = aux1
               h(j+1,i) = aux2
  340       continue
            r = dsqrt(h(i,i)*h(i,i) + h(i+1,i)*h(i+1,i))
            c(i) = h(i,i)/r
            s(i) = h(i+1,i)/r
            h(i,i)   = r
            h(i+1,i) = 0.d0
            e(i+1) = -s(i) * e(i)
            e(i)   =  c(i) * e(i)
            if (dabs(e(i+1)) .le. econv) goto 500
  400    continue
  500    continue
c
c ...... Resolve o sistema h y = e :
c
         y(ni) = e(ni) / h(ni,ni)
         do 520 i = ni-1, 1, -1
            y(i) = 0.d0
            do 510 j = i+1, ni
               y(i) = y(i) - h(i,j)*y(j)
  510       continue
            y(i) = ( y(i) + e(i)) / h(i,i)
  520    continue
c
c ...... Atualizacao de x:
c
c$omp parallel do num_threads(nth_solv) 
         do 610 i = 1, neq
            do 600 j = 1, ni
               x(i) = x(i) + y(j) * g(i,j)
  600       continue
  610    continue
c$omp end parallel do
c ......................................................................
c
c ...
         jj = jj + 1
         if( jj .eq. 10) then
           jj = 0
           write(*,2300),l,nit,dabs(e(ni+1)),econv
         endif
c ......................................................................
c
c ...... Verifica a convergencia:
c
         if (dabs(e(ni+1)) .le. econv) goto 1100
 1000 continue
c ----------------------------------------------------------------------
 1100 continue
c
c ... Norma de energia da solucao:
c
c$omp parallel  num_threads(nth_solv) 
      call matvec(neq,ia,ja,ia(neq+2),ja(nad+1),ad,al,au,al(nad+1)
     .           ,x     ,g(1,1)  
     .           ,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli
     .           ,thread_y)
      aux1 = dot(x,g(1,1),neq_doti)
c$omp single
      energy = aux1
c$omp end single
c$omp end parallel
c ......................................................................
      time = MPI_Wtime()
      time = time-time0
c ......................................................................
      if (dabs(e(ni+1)) .gt. econv) then
         if(my_id .eq. 0) then
            write(*,2100) maxit,k,nit
           if(flog) write(10,2100) maxit,k,nit
         endif 
         call stop_mef()
      endif
c .....................................................................
      if(my_id.eq.0)write(*,2000) tol,neq,l,nit,dabs(e(ni+1)),energy
     .                           ,time
c ......................................................................
      if(flog) then
      if(my_id.eq.0) write(10,'(a,a,i9,a,d20.10,a,d20.10,a,f20.2)')
     .         "GMRES_OMP: "," it ",nit, " energy norm ",energy
     .        ," tol ",tol," time ",time
      endif
c .....................................................................
      return
c ----------------------------------------------------------------------
 2000 format(' (GMRES_OMP) solver:'/
     . 5x,'Solver tol           = ',d20.6/
     . 5x,'Number of equations  = ',i20/
     . 5x,'Number of cycles     = ',i20/
     . 5x,'Number of iterations = ',i20/
     . 5x,'Norm                 = ',d20.10/
     . 5x,'Energy norm          = ',d20.10/
     . 5x,'CPU time (s)         = ',f20.2/)
 2100 format(' *** WARNING: no convergence reached for '
     .      ,i9,' cycles !',5x,i7,' nKylov',5x,' It ',i7/)
 2300 format (' GMRES_OMP:',5x,'cycles',i7,5x,'It',i7,5x,2d20.10)
      end
c **********************************************************************
      subroutine pbicgstab_omp(neq   ,nequ   ,nad    ,ia,ja
     .                        ,ad    ,au     ,al     ,m ,b
     .                        ,x     ,t      ,v      ,r ,p ,z
     .                        ,tol   ,maxit
     .                        ,matvec,dot
     .                        ,my_id ,neqf1i ,neqf2i,neq_doti,i_fmapi
     .                        ,i_xfi ,i_rcvsi,i_dspli,thread_y,flog)
c **********************************************************************
c * Data de criacao    : 00/00/0000                                    *
c * Data de modificaco : 12/12/2015                                    * 
c * ------------------------------------------------------------------ *   
c * PBICGSTAB_OMP : Solucao de sistemas de equacoes pelo metodo dos    * 
c * gradientes biconjugados com precondicionador diagonal para         *
c * matrizes nao-simetricas.                                           *
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ *                                                                   *
c * neq      - numero de equacoes                                      *
c * nequ     - numero de equacoes no bloco Kuu                         *
c * nad      - numero de termos nao nulos no bloco Kuu e Kpu  ou K     *
c * ia(*)    - ponteiro do formato CSR                                 *
c * ja(*)    - ponteiro das colunas no formato CSR                     *
c * ad(neq)  - diagonal da matriz A                                    *
c * au(*)    - parte triangular superior de A                          *
c * al(*)    - parte triangular inferior de A                          *
c * m(*)     - precondicionador diagonal                               *
c * b(neq)   - vetor de forcas                                         *
c * x(neq)   - chute inicial                                           *
c * t(neq)   - arranjo local de trabalho                               *
c * v(neq)   - arranjo local de trabalho                               *
c * r(neq)   - arranjo local de trabalho                               *
c * p(neq)   - arranjo local de trabalho                               *
c * z(neq)   - arranjo local de trabalho                               *
c * tol      - tolerancia de convergencia                              *
c * maxit    - numero maximo de iteracoes                              *
c * matvec   - nome da funcao externa para o produto matrix-vetor      *
c * dot      - nome da funcao externa para o produto escalar           *
c * my_id    -                                                         *
c * neqf1i   -                                                         *
c * neqf2i   -                                                         *
c * neq_doti -                                                         *
c * i_fmap   -                                                         *
c * i_xfi    -                                                         *
c * i_rvcs   -                                                         *
c * i_dspli  -                                                         *
c * flog     - log do arquivo de saida                                 *
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ *
c * x(neq) - vetor solucao                                             *
c * b(neq) - modificado                                                *
c * ad(*),al(*),au(*) - inalterados                                    *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ *   
c **********************************************************************
      implicit none
      include 'mpif.h'
      include 'omp_lib.h'
      include 'openmp.fi'
c ... Mpi      
      integer ierr
c .....................................................................      
      integer neqf1i,neqf2i,neq_doti
c ... ponteiros      
      integer*8 i_fmapi,i_xfi
      integer*8 i_rcvsi,i_dspli      
c .....................................................................      
      integer neq,nequ,maxit,nad,i,j,jj,k
      integer ia(*),ja(*),my_id
      real*8 ad(*),au(*),al(*),m(*),x(*),r(*),p(*),b(*),t(*),v(*),z(*)
      real*8  dot,ddot,tol,conv,energy,d,alpha,beta,rr0,w
      real*8  time0,time,dottmp1,dottmp2
      real*8 thread_y(*) 
      logical flog
      external matvec,dot
c ======================================================================
      time0 = MPI_Wtime()
c ......................................................................
c$omp parallel default(none) 
c$omp.private(i,j,jj,d,conv,beta,alpha,rr0,w,energy) 
c$omp.shared(neq,nequ,nad,ia,ja,al,ad,au)
c$omp.shared(m,x,r,p,t,v,z,tol,maxit,thread_y)
c$omp.shared(neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli,neq_doti,flog)
c$omp.shared(my_id,time,time0)
c$omp.num_threads(nth_solv)                                          
c
c ... Chute inicial:
c
c$omp do
      do 10 i = 1, neq
         x(i) = 0.d0
   10 continue
c$omp end do
c ----------------------------------------------------------------------
      call matvec(neq,ia,ja,ia(neq+2),ja(nad+1),ad,al,au,al(nad+1) 
     .           ,x,z
     .           ,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli,thread_y)
c$omp do
      do 100 i = 1, neq
         r(i) = b(i) - z(i)
         p(i) = r(i)
         b(i) = p(i)
         z(i) = p(i)*m(i) 
  100 continue
c$omp end do
      d = dot(r,z,neq_doti)
      conv = tol*dsqrt(dabs(d))
c ----------------------------------------------------------------------
      jj = 1
      do 230 j = 1, maxit
         call matvec(neq,ia,ja,ia(neq+2),ja(nad+1),ad,al,au,al(nad+1) 
     .              ,z,v
     .              ,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli
     .              ,thread_y)
         rr0   = dot(b,r,neq_doti)
         alpha = rr0 / dot(v,r,neq_doti)
c$omp do          
         do 210 i = 1, neq
            x(i) = x(i) + alpha * z(i)
            b(i) = b(i) - alpha * v(i)
            z(i) = b(i) * m(i)
  210    continue
c$omp end do
         call matvec(neq,ia,ja,ia(neq+2),ja(nad+1),ad,al,au,al(nad+1) 
     .              ,z,t
     .              ,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli
     .              ,thread_y)
         w = dot(t,b,neq_doti)/dot(t,t,neq_doti)
c$omp do         
         do 220 i = 1, neq
             x(i) = x(i) + w*z(i)
             b(i) = b(i) - w*t(i)
  220    continue
c$omp end do
         d = dot(b,z,neq_doti)
         if (dsqrt(dabs(d)) .lt. conv) goto 300
         beta = (dot(r,b,neq_doti) / rr0)*(alpha/w)
c$omp do         
         do 225 i = 1, neq
              p(i) = b(i) + beta*(p(i)-w*v(i))
              z(i) = p(i)*m(i)
  225    continue
c$omp end do  
c ......................................................................
c$omp master
         if( jj .eq.500) then
           jj = 0
           write(*,1300),j,dsqrt(dabs(d)),conv 
         endif  
         jj = jj + 1
c$omp end master
c ......................................................................
  230 continue
c ----------------------------------------------------------------------
c$omp single      
      write(*,1200) maxit
      if(flog) write(10,1200) maxit
      call stop_mef()
c$omp end single      
  300 continue
c
c ... Energy norm:
c
      call matvec(neq,ia,ja,ia(neq+2),ja(nad+1),ad,al,au,al(nad+1) 
     .           ,x,z
     .           ,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli
     .           ,thread_y)
      energy   = dot(x,z,neq_doti)
c ......................................................................
c$omp single
      time = MPI_Wtime()
      time = time-time0
c ----------------------------------------------------------------------
      if(my_id.eq.0)write(*,1100) tol,neq,j,energy,time
c ......................................................................
c     Controle de flops
      if(my_id.eq.0) write(10,'(a,a,i9,a,d20.10,a,d20.10,a,f20.2)')
     .               "PBICGSTAB_OMP: ","it",j, " energy norm ",energy,
     .               " tol ",tol," time ",time
c ......................................................................
c$omp end single
c$omp end parallel
      return
c ======================================================================
 1000 format (//,5x,'SUBROTINA PBICGSTAB:',/,5x,'Coeficiente da diagonal
     . nulo ou negativo - equacao ',i7)
 1100 format(' (PBICGSTAB_OMP) solver:'/
     . 5x,'Solver tol           = ',d20.6/
     . 5x,'Number of equations  = ',i20/
     . 5x,'Number of iterations = ',i20/
     . 5x,'Energy norm          = ',d20.10/
     . 5x,'CPU time (s)         = ',f20.2/)
 1200 format (' *** WARNING: No convergence reached after ',i9,
     .        ' iterations !',/)
 1300 format (' BICGSTAB:',5x,'It',i7,5x,2d20.10)
      end
c **********************************************************************
      subroutine pcg_omp_loopwise(neq,ia,ja,ad,au,al,m,b,x,z,r,tol,
     .              maxit,matvec,dot,my_id,neqf1i,neqf2i,neq_doti,
     .              i_fmapi,i_xfi,i_rcvsi,i_dspli,thread_y)
c **********************************************************************
c *                                                                    *
c *   Subroutine PCG                                                   *
c *                                                                    *
c *   Solucao de sistemas de equacoes pelo metodo dos gradientes       *
c *   conjugados com precondicionador diagonal para matrizes           *
c *   simetricas.                                                      *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *                                                                    *
c *   neq    - numero de equacoes                                      *
c *   ia(*)  - ponteiro do formato CSR                                 *
c *   ja(*)  - ponteiro das colunas no formato CSR                     *
c *   ad(neq)- diagonal da matriz A                                    *
c *   au(*)  - parte triangular superior de A                          *
c *   al(*)  - parte triangular inferior de A                          *
c *   m(*)   - precondicionador diagonal                               *
c *   b(neq) - vetor de forcas                                         *
c *   x(neq) - chute inicial                                           *
c *   z(neq) - arranjo local de trabalho                               *
c *   r(neq) - arranjo local de trabalho                               *
c *   tol    - tolerancia de convergencia                              *
c *   maxit  - numero maximo de iteracoes                              *
c *   matvec - nome da funcao externa para o produto matrix-vetor      *
c *   dot    - nome da funcao externa para o produto escalar           *
c *   energy - nao definido                                            *
c *   thread_y - buffer de equacoes para o vetor y (openmp)            *
c *                                                                    *
c *   Parametros de saida:                                             *
c *                                                                    *
c *   x(neq) - vetor solucao                                           *
c *   b(neq) - modificado                                              *
c *   ad(*),al(*),au(*) - inalterados                                  *
c *                                                                    *
c **********************************************************************
      implicit none
      include 'mpif.h'
      include 'omp_lib.h'
      include 'openmp.fi'
      integer neqf1i,neqf2i,neq_doti
c ... ponteiros      
      integer*8 i_fmapi,i_xfi
      integer*8 i_rcvsi,i_dspli
c .....................................................................
      integer neq,maxit,i,j,nad
      integer ia(*),ja(*),my_id
      real*8  ad(*),au(*),al(*),m(*),x(*),r(*),z(*),b(*)
      real*8  dot,tol,conv,energy,d,alpha,beta
      real*8  time0,time
      real*8  dottmp
      real*8 thread_y(*)
      external matvec,dot
c ======================================================================
      time0 = MPI_Wtime()
c ......................................................................
      nad = ia(neq+1)-1
      if(my_id.eq.0)print *, 'nad :',nad
c ......................................................................
c
c ... Chute inicial:
c
c$omp parallel do
      do 10 i = 1, neq
         x(i) = 0.d0
   10 continue
c$omp end parallel do
c ----------------------------------------------------------------------
c$omp parallel
      call matvec(neq,ia,ja,ia(neq+2),ja(nad+1),ad,al,au,al(nad+1),x,z,
     .            neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli,thread_y)
c$omp end parallel
c$omp parallel do
      do 100 i = 1, neq
         r(i) = b(i) - z(i)
         z(i) = r(i) / m(i)
         b(i) = z(i)
  100 continue
c$omp end parallel do
      d    = dot(r(1),z(1),neq_doti)
      conv = tol*dsqrt(dabs(d))
c ----------------------------------------------------------------------
      do 230 j = 1, maxit
c$omp parallel
         call matvec(neq,ia,ja,ia(neq+2),ja(nad+1),ad,al,au,al(nad+1),
     .               b,z,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli,
     .               thread_y)
c$omp end parallel
         dottmp = dot(b,z,neq_doti)
         alpha = d / dottmp
c$omp parallel do
         do 210 i = 1, neq
            x(i) = x(i) + alpha * b(i)
            r(i) = r(i) - alpha * z(i)
            z(i) = r(i) / m(i)
  210    continue
c$omp end parallel do
         dottmp = dot(r,z,neq_doti)
         beta = dottmp / d
c$omp parallel do
         do 220 i = 1, neq
            b(i) = z(i) + beta * b(i)
  220    continue
c$omp end parallel do
         d = beta * d
         if (dsqrt(dabs(d)) .lt. conv) goto 300
  230 continue
c ----------------------------------------------------------------------
      write(*,1200) maxit
      call stop_mef()
  300 continue
c
c ... Energy norm:
c
c$omp parallel
      call matvec(neq,ia,ja,ia(neq+2),ja(nad+1),ad,al,au,al(nad+1),x,z,
     .            neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli,thread_y)
c$omp end parallel
      energy = dot(x(1),z(1),neq_doti)
c ......................................................................
      time = MPI_Wtime()
      time = time-time0
c ----------------------------------------------------------------------
      if(my_id.eq.0)write(*,1100) tol,neq,j,energy,time
c ......................................................................
c     Controle de flops
      if(my_id.eq.0) write(10,'(a,a,i9,a,d20.10,a,d20.10)')
     .               "PCG_OMP_LOOPWISE: ",
     .               "it",j, " energy norm ",energy," tol ",tol
c ......................................................................
      return
c ======================================================================
 1000 format (//,5x,'SUBROTINA PCG:',/,5x,'Coeficiente da diagonal nulo
     .ou negativo - equacao ',i7)
 1100 format(' (PCG_OMP_LOOPWISE) solver:'/
     . 5x,'Solver tol           = ',d20.6/
     . 5x,'Number of equations  = ',i20/
     . 5x,'Number of iterations = ',i20/
     . 5x,'Energy norm          = ',d20.10/
     . 5x,'CPU time (s)         = ',f20.2/)
 1200 format (' *** WARNING: No convergence reached after ',i9,
     .        ' iterations !',/)
      end
      subroutine gmres_omp_loopwise(neq,ia,ja,ad,au,al,m,b,x,k,g,h,y,c,
     .             s,e,tol,maxit,matvec,dot,neqovlp,my_id,
     .             neqf1i,neqf2i,neq_doti,i_fmapi,i_xfi,i_rcvsi,i_dspli,
     .             thread_y) 
c **********************************************************************
c *                                                                    *
c *   GMRES: Solucao iterativa de sistemas simetricos e nao-simetricos *
c *          pelo metodo GMRES com precondicionador diagonal.          *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *                                                                    *
c *   neq    - numero de equacoes                                      *
c *   ia(*)  - ponteiro do formato CSR                                 *
c *   ja(*)  - ponteiro das colunas no formato CSR                     *
c *   ad(neq)- diagonal da matriz A                                    *
c *   au(*)  - parte triangular superior de A                          *
c *   al(*)  - parte triangular inferior de A                          *
c *   m(*)   - precondicionador diagonal                               *
c *   b(neq) - vetor de forcas                                         *
c *   x(neq) - chute inicial                                           *
c *   z(neq) - arranjo local de trabalho                               *
c *   r(neq) - arranjo local de trabalho                               *
c *   tol    - tolerancia de convergencia                              *
c *   maxit  - numero maximo de iteracoes                              *
c *   matvec - nome da funcao externa para o produto matrix-vetor      *
c *   dot    - nome da funcao externa para o produto escalar           *
c *   thread_y - buffer de equacoes para o vetor y (openmp)            *
c *                                                                    *
c *   Arranjos locais de trabalho:                                     *
c *                                                                    *
c *      g(neq+1,k+1)                                                  *
c *      h(k+1,k)                                                      *
c *      y(k)                                                          *
c *      c(k)                                                          *
c *      s(k)                                                          *
c *      e(k+1)                                                        *
c *                                                                    *
c *   Parametros de saida:                                             *
c *                                                                    *
c *   x(neq) - vetor solucao                                           *
c *   b(neq),ad(*),al(*),au(*) - modificados                           *
c *                                                                    *
c **********************************************************************
      implicit none
      include 'mpif.h'
      include 'omp_lib.h'

      integer neqf1i,neqf2i,neq_doti
c ... ponteiros      
      integer*8 i_fmapi,i_xfi
      integer*8 i_rcvsi,i_dspli
c .....................................................................      
      include 'openmp.fi'
      integer neq,k,maxit,ia(*),ja(*),neqovlp,nit,i,j,l,ni,ic,nad,nad1
      real*8  ad(neq),au(*),al(*),m(*),b(*),x(*)
      real*8  g(neqovlp,1:k+1),h(k+1,k),y(k),c(k),s(k),e(k+1),tol
      real*8  energy,econv,norm,dot,r,aux1,aux2,beta, gnorm
      real*8  time0,time
      real*8 thread_y(*)
      external matvec,dot
      integer nii(maxit),my_id
      integer ns
c ......................................................................
      time0 = MPI_Wtime()
c ......................................................................
c     call omp_set_num_threads(num_threads)
      nad = ia(neq+1)-1
      if(my_id.eq.0)print*,nad
c ----------------------------------------------------------------------
c
c.... Chute inicial:
c
c$omp parallel do
      do 10 i = 1, neq
         x(i) = 0.d0
   10 continue
c$omp end parallel do
c ----------------------------------------------------------------------
c
c ... Limite de convergencia:
c
      norm  = dsqrt(dot(b(1),b(1),neq_doti))
      econv = tol*norm
c ----------------------------------------------------------------------
c
c ... Ciclos GMRES:
c
      nit = 0
      do 1000 l = 1, maxit
c
c ...... Residuo g(1) = b - A x:
c
c$omp parallel
         call matvec(neq,ia,ja,ia(neq+2),ja(nad+1),ad,al,au,al(nad+1),
     .               x,g(1,1),neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,
     .               i_dspli,thread_y)
c$omp end parallel
c
c ...... Residuo com precondicionador diagonal:
c
c$omp parallel do
         do 200 i = 1, neq
            g(i,1) = (b(i) - g(i,1))/m(i)
  200    continue
c$omp end parallel do
c
c ...... Norma do residuo:
c
         e(1) = dsqrt(dot(g(1,1),g(1,1),neq_doti))
c
c ...... Normalizacao de g1:
c
c$omp parallel do
         do 210 i = 1, neq
            g(i,1) = g(i,1)/e(1)
  210    continue
c$omp end parallel do
c
c ...... Iteracoes GMRES:
c
         ni = 0
         do 400 i = 1, k
            nit =nit + 1
            ni  = ni  + 1
c
c ......... Produto g(i+1) = A.g(i):
c
c$omp parallel
            call matvec(neq,ia,ja,ia(neq+2),ja(nad+1),ad,al,au,
     .                  al(nad+1),g(1,i),g(1,i+1),neqf1i,neqf2i,i_fmapi,
     .                  i_xfi,i_rcvsi,i_dspli,thread_y)
c$omp end parallel
c
c ......... Precondicionador diagonal:
c
c$omp parallel do
            do 300 j = 1, neq
               g(j,i+1) = g(j,i+1)/m(j)
  300       continue
c$omp end parallel do
c
c ......... Ortogonalizacao (Gram-Schmidt modificado):
c
            do 320 j = 1, i
               beta = dot(g(1,i+1),g(1,j),neq_doti)
c$omp parallel do
               do 310 ic = 1, neq
                  g(ic,i+1) = g(ic,i+1) - beta * g(ic,j)
  310          continue
c$omp end parallel do
               h(j,i) = beta
  320       continue
c
c ......... Norma de g(i+1):
c
            norm = dsqrt(dot(g(1,i+1),g(1,i+1),neq_doti))
c
            h(i+1,i) = norm
c
c ......... Normalizacao de g(i+1):
c
c$omp parallel do
            do 330 ic = 1, neq
               g(ic,i+1) = g(ic,i+1)/norm
  330       continue
c$omp end parallel do
  400    continue
c
c ..... Givens QR factorization of the Hessenberg matrix
c
         r = dsqrt(h(1,1)*h(1,1) + h(2,1)*h(2,1))
         c(1)   = h(1,1)/r
         s(1)   = h(2,1)/r
         h(1,1) = r
         h(2,1) = 0.d0
         e(2) = -s(1) * e(1)
         e(1) =  c(1) * e(1)
         ns = 1
         if (dabs(e(2)) .le. econv) goto 500
         do 450 i = 2, k
c$omp parallel do
            do 340 j = i, k
               aux1 =  c(i-1) * h(i-1,j) + s(i-1) * h(i,j)
               aux2 = -s(i-1) * h(i-1,j) + c(i-1) * h(i,j)
               h(i-1,j)   = aux1
               h(i,j)     = aux2
  340       continue
c$omp end parallel do
c ......... New Givens matrix elements
            r = dsqrt(h(i,i)*h(i,i) + h(i+1,i)*h(i+1,i))
            c(i) = h(i,i)/r
            s(i) = h(i+1,i)/r
            h(i,i)   = r
            h(i+1,i) = 0.d0
            e(i+1) = -s(i) * e(i)
            e(i)   =  c(i) * e(i)
            ns = i
            if (dabs(e(i+1)) .le. econv) goto 500
  450    continue
  500    continue
c
c ...... Resolve o sistema h y = e :
c
         y(ns) = e(ns) / h(ns,ns)
         do 520 i = ns-1, 1, -1
            y(i) = 0.d0
            do 510 j = i+1, ns
               y(i) = y(i) - h(i,j)*y(j)
  510       continue
            y(i) = (y(i) + e(i)) / h(i,i)
  520    continue
c
c ...... Atualizacao de x:
c
c TODO: junior dense matvec with columns arrangement
c$omp parallel do
         do 610 i = 1, neq
            do 600 j = 1, ns
               x(i) = x(i) + y(j) * g(i,j)
  600       continue
  610    continue
c$omp end parallel do
c
c ...... Verifica a convergencia:
c
         nii(l)=ni
         if (dabs(e(ns+1)) .le. econv) goto 1100
 1000 continue
c ----------------------------------------------------------------------
 1100 continue
c
c ... Norma de energia da solucao:
c
      energy = dot(x,b,neq_doti)
c ......................................................................
      time = MPI_Wtime()
      time = time-time0
c ----------------------------------------------------------------------
      if(my_id.eq.0)write(*,2000) tol,neq,l,nit,dabs(e(ns+1)),energy
     .                           ,time
      if (dabs(e(ns+1)) .gt. econv) then
         write(*,2100) maxit
c         stop
      endif
c ......................................................................
c     Controle de flops
      if(my_id.eq.0)write(10,'(999(i4,1x))') l,nit,(nii(j),j=1,l)
c ......................................................................
      return
c ----------------------------------------------------------------------
 2000 format(' (GMRES_OMP_LOOPWISE) solver:'/
     . 5x,'Solver tol           = ',d20.6/
     . 5x,'Number of equations  = ',i20/
     . 5x,'Number of cycles     = ',i20/
     . 5x,'Number of iterations = ',i20/
     . 5x,'Norm                 = ',d20.6/
     . 5x,'Energy norm          = ',d20.10/
     . 5x,'CPU time (s)         = ',f20.2/)
 2100 format(' *** WARNING: no convergence reached for ',i9,' cycles !',
     . /)
      end
c **********************************************************************
      subroutine pbicgstab_omp_loopwise(neq,ia,ja,ad,au,al,m,b,x,t,v,r,
     .            p,z,tol,maxit,matvec,dot,my_id,neqf1i,neqf2i,
     .            neq_doti,i_fmapi,i_xfi,i_rcvsi,i_dspli,thread_y)
c **********************************************************************
c *                                                                    *
c *   Subroutine PBICGSTAB                                             *
c *                                                                    *
c *   Solucao de sistemas de equacoes pelo metodo dos gradientes       *
c *   biconjugados com precondicionador diagonal para matrizes         *
c *   simetricas.                                                      *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *                                                                    *
c *   neq    - numero de equacoes                                      *
c *   ia(*)  - ponteiro do formato CSR                                 *
c *   ja(*)  - ponteiro das colunas no formato CSR                     *
c *   ad(neq)- diagonal da matriz A                                    *
c *   au(*)  - parte triangular superior de A                          *
c *   al(*)  - parte triangular inferior de A                          *
c *   m(*)   - precondicionador diagonal                               *
c *   b(neq) - vetor de forcas                                         *
c *   x(neq) - chute inicial                                           *
c *   t(neq) - arranjo local de trabalho                               *
c *   v(neq) - arranjo local de trabalho
c *   r(neq) - arranjo local de trabalho                               *
c *   p(neq) - arranjo local de trabalho
c *   z(neq) - arranjo local de trabalho                               *
c *   tol    - tolerancia de convergencia                              *
c *   maxit  - numero maximo de iteracoes                              *
c *   matvec - nome da funcao externa para o produto matrix-vetor      *
c *   dot    - nome da funcao externa para o produto escalar           *
c *   energy - nao definido                                            *
c *   thread_y - buffer de equacoes para o vetor y (openmp)            *
c *                                                                    *
c *   Parametros de saida:                                             *
c *                                                                    *
c *   x(neq) - vetor solucao                                           *
c *   b(neq) - modificado                                              *
c *   ad(*),al(*),au(*) - inalterados                                  *
c *                                                                    *
c **********************************************************************
      implicit none
      include 'mpif.h'
c ... Mpi      
      integer ierr
c .....................................................................      
      integer neqf1i,neqf2i,neq_doti
c ... ponteiros      
      integer*8 i_fmapi,i_xfi
      integer*8 i_rcvsi,i_dspli      
c .....................................................................      
      integer neq,maxit,nad,i,j,k
      integer ia(*),ja(*),my_id
      real*8  ad(*),au(*),al(*),m(*),x(*),r(*),p(*),b(*),t(*),v(*),z(*)
      real*8  dot,ddot,tol,conv,energy,d,alpha,beta,rr0,w
      real*8  time0,time
      real*8  thread_y(*)
      external matvec,dot
c ======================================================================
      time0 = MPI_Wtime()
c ......................................................................
      nad = ia(neq+1)-1
      if(my_id.eq.0)print *, 'nad :',nad
c ......................................................................
c
c ... Chute inicial:
c
c$omp parallel do private(i)
      do 10 i = 1, neq
         x(i) = 0.d0
   10 continue
c$omp end parallel do
c ----------------------------------------------------------------------
c$omp parallel
      call matvec(neq,ia,ja,ia(neq+2),ja(nad+1),ad,al,au,al(nad+1),x,z,
     .            neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli,thread_y)
c$omp end parallel
c
c$omp parallel do private(i)
      do 100 i = 1, neq
         r(i) = b(i) - z(i)
         p(i) = r(i)
         b(i) = p(i)
         z(i) = p(i)/m(i) 
  100 continue
c$omp end parallel do
      d    = dot(r(1),z(1),neq_doti)
      conv = tol*dsqrt(dabs(d))
c ----------------------------------------------------------------------
      do 230 j = 1, maxit
c$omp parallel      
         call matvec(neq,ia,ja,ia(neq+2),ja(nad+1),ad,al,au,al(nad+1),
     .               z,v,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli,
     .               thread_y)
c$omp end parallel     
         rr0 = dot(b,r,neq_doti)
         alpha = rr0/dot(v,r,neq_doti)
c$omp parallel do private(i)         
         do 210 i = 1, neq
            x(i) = x(i) + alpha * z(i)
            b(i) = b(i) - alpha * v(i)
            z(i) = b(i) / m(i)
  210    continue
c$omp end parallel do   
c$omp parallel  
         call matvec(neq,ia,ja,ia(neq+2),ja(nad+1),ad,al,au,al(nad+1),
     .               z,t, neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli,
     .               thread_y)
c$omp end parallel     
         w = dot(t,b,neq_doti) / dot(t,t,neq_doti)
c$omp parallel do private(i)        
         do 220 i = 1, neq
            x(i) = x(i) + w*z(i)
            b(i) = b(i) - w*t(i)
  220    continue
c$omp end parallel do
         d = dot(b,z,neq_doti)
         if (dsqrt(dabs(d)) .lt. conv) goto 300   
         beta = (dot(r,b,neq_doti) / rr0)*(alpha/w)
c$omp parallel do private(i)         
         do 225 i = 1, neq
             p(i) = b(i) + beta*(p(i)-w*v(i))
             z(i) = p(i)/m(i)
  225    continue
c$omp end parallel do  

  230 continue
c ----------------------------------------------------------------------
      write(*,1200) maxit
      call stop_mef()
  300 continue
c
c ... Energy norm:
c
c$omp parallel
      call matvec(neq,ia,ja,ia(neq+2),ja(nad+1),ad,al,au,al(nad+1),x,z,
     .            neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli,thread_y)
c$omp end parallel
      energy   = dot(x(1),z(1),neq_doti)
c ......................................................................
      time = MPI_Wtime()
      time = time-time0
c ----------------------------------------------------------------------
      if(my_id.eq.0)write(*,1100) tol,neq,j,energy,time
c ......................................................................
c     Controle de flops
      if(my_id.eq.0) write(10,'(a,a,i9,a,d20.10,a,d20.10)')
     .               "PBICGSTAB_OMP: ","it",j, " energy norm ",energy,
     .               " tol ",tol
c ......................................................................
      return
c ======================================================================
 1000 format (//,5x,'SUBROTINA PBICGSTAB:',/,5x,'Coeficiente da diagonal
     . nulo ou negativo - equacao ',i7)
 1100 format(' (PBICGSTAB_OMP_LOOPWISE) solver:'/
     . 5x,'Solver tol           = ',d20.6/
     . 5x,'Number of equations  = ',i20/
     . 5x,'Number of iterations = ',i20/
     . 5x,'Energy norm          = ',d20.10/
     . 5x,'CPU time (s)         = ',f20.2/)
 1200 format (' *** WARNING: No convergence reached after ',i9,
     .        ' iterations !',/)
      end
c *********************************************************************

