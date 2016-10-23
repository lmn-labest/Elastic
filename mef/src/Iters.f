c *********************************************************************
c * Metodos iterativos para solucao de sistemas lineares              *
c * ----------------------------------------------------------------- *
c * simetricos:                                                       *
c * ----------------------------------------------------------------- *
c * CG   - gradiente conjugados                                       *
c *                                                                   *
c * PCG  - gradiente conjugados com precondicionador diagonal         *
c *                                                                   *
c * ICCG - gradiente conjugados com precondicionador de fatoracoes    *              
c * incompletas                                                       *
c *                                                                   *
c * ----------------------------------------------------------------- *
c * nao-simetricos:                                                   *
c * ----------------------------------------------------------------- *
c * bicgstab - gradiente bi-conjugados estabilizados                  *
c *                                                                   *
c * pbicgstab - gradiente bi-conjugados estabilizados  com            * 
c * precondicionador diagonal                                         *
c *                                                                   *
c * icbicgstab - gradiente bi-conjugados estabilizados fatoracoes     *              
c * incompletas                                                       *   
c *                                                                   *
c * gmres(m) - GMRES com precondicionador diagonal                    *
c *                                                                   *
c * ----------------------------------------------------------------- *
c *********************************************************************  
      subroutine cg(neq   ,nad   ,ia    ,ja
     .             ,ad    ,au    ,al    ,b       ,x
     .             ,z     ,r     ,p     ,tol     ,maxit
     .             ,matvec,dot
     .             ,my_id ,neqf1i ,neqf2i,neq_doti,i_fmapi
     .             ,i_xfi ,i_rcvsi,i_dspli
     .             ,fprint,flog   ,fnew)
c **********************************************************************
c * Data de criacao    : 00/00/0000                                    *
c * Data de modificaco : 24/04/2016                                    * 
c * ------------------------------------------------------------------ *   
c * Subroutine CG : Solucao de sistemas de equacoes pelo metodo dos    *    
c * gradientes conjugados                                              *
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ * 
c * neq      - numero de equacoes                                      *
c * nad      - numero de termos nao nulos                              *
c * ia(*)    - ponteiro do formato CSR                                 *
c * ja(*)    - ponteiro das colunas no formato CSR                     *
c * ad(neq)  - diagonal da matriz A                                    *
c * au(*)    - parte triangular superior de A                          *
c * al(*)    - parte triangular inferior de A                          *
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
c * fprint   - saida na tela                                           *
c * flog     - log do arquivo de saida                                 *
c * fnew     - .true.  x0 igual a zero                                 *
c *            .false. x0 dado                                         *
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ *
c * x(neq) - vetor solucao                                             *
c * ad(*),al(*),au(*) e b - inalterados                                *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ * 
c **********************************************************************
      implicit none
      include 'mpif.h'
      integer neqf1i,neqf2i,neq_doti
c ... ponteiros      
      integer*8 i_fmapi,i_xfi
      integer*8 i_rcvsi,i_dspli
c .....................................................................      
      integer neq,nequ,nad,maxit,i,j,jj
      integer ia(*),ja(*),my_id
      real*8  ad(*),au(*),al(*),x(*),b(*)
      real*8  r(*),z(*),p(*)
      real*8  dot,tol,conv,xkx,norm,d,di,alpha,beta,tmp
      real*8  time0,time
      real*8 dum1
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
c
c ... Chute inicial:
c
      if(fnew) then  
        do 10 i = 1, neq
          x(i)  = 0.d0
   10   continue
      endif  
c .......................................................................
c
c ... conv = tol * |b|
      d    = dot(b,b,neq_doti)
      conv = tol*dsqrt(dabs(d))
c .......................................................................
c  
c ... Ax0
      call matvec(neq,ia,ja,ia(neq+2),ja(nad+1),ad,al,au,al(nad+1)  
     .           ,x,z 
     .           ,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli,dum1)
c .......................................................................
c
c ...
      do 100 i = 1, neq
c ... r0 = b - Ax0
         r(i) = b(i) - z(i)
c ... p0 = r0
         p(i) = r(i)
  100 continue
c ... ( r(0),r(0) )
      d    = dot(r,r,neq_doti)
c ......................................................................
      jj = 1
      do 230 j = 1, maxit
c ... z = Ap(j)
         call matvec(neq,ia,ja,ia(neq+2),ja(nad+1),ad,al,au,al(nad+1) 
     .              ,p,z,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli 
     .              ,dum1)
c .....................................................................
c
c ... alpha = ( r(j),r(j) ) / ( Ap(j), p(j) ))
         alpha = d / dot(z,p,neq_doti)
c .....................................................................
c
c ...
         do 210 i = 1, neq
c ... x(j+1) = x(j) + alpha*p
            x(i) = x(i) + alpha * p(i)
c ... r(j+1) = r(j) - alpha*Ap
            r(i) = r(i) - alpha * z(i)
  210    continue
c .....................................................................
c
c ...    
         di   = dot(r,r,neq_doti) 
c ... beta = ( r(j+1),r(j+1) ) / ( r(j),r(j) )
         beta = di / d
c .....................................................................
c
c ...
         do 220 i = 1, neq
c ... p(j+1) = r(j+1) + beta*p(j)
            p(i) = r(i) + beta * p(i)
  220    continue
c .....................................................................
c
c ...
         d =  di
         if (dsqrt(dabs(d)) .lt. conv) goto 300
c ......................................................................
c
c ...
         if( jj .eq.500) then
           jj = 0
           write(*,1300),j,dsqrt(dabs(d)),conv 
         endif  
         jj = jj + 1
c ......................................................................
  230 continue
c ......................................................................
      write(*,1200) maxit
      if(flog) write(10,1200) maxit
      call stop_mef()
  300 continue
c
c ... produto:  x*Kx
c
      call matvec(neq,ia,ja,ia(neq+2),ja(nad+1),ad,al,au,al(nad+1) 
     .           ,x,z 
     .           ,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli,dum1)
       xkx = dot(x,z,neq_doti)
c ......................................................................
c
c ... norm-2 = || x ||
      norm = dsqrt(dot(x,x,neq_doti))
c ......................................................................
c
c ... r = b - Ax (calculo do residuo explicito)
      do 310 i = 1, neq
        r(i) = b(i) - z(i)
  310 continue
      tmp  = dot(r,r,neq_doti)
      tmp = dsqrt(tmp)
      if( tmp .gt. conv ) then
        write(*,1400) tmp
      endif
c ......................................................................
      time = MPI_Wtime()
      time = time-time0
c ......................................................................
      if(my_id .eq.0 .and. fprint )then
        write(*,1100)tol,conv,neq,nad,j,xkx,norm,time
      endif
c ......................................................................
c     Controle de flops
      if(flog) then
        if(my_id.eq.0) then
          write(10,'(a,a,i9,a,d20.10,a,d20.10,a,d20.10,a,f20.2)')
     .       "CG: "," it ",j, " x * Kx ",xkx," ||x|| ",norm
     .      ," tol ",tol," time ",time
        endif
      endif
c ......................................................................
      return
c ======================================================================
 1000 format (//,5x,'SUBROTINA CG:',/,5x,'Coeficiente da diagonal ' 
     . '- equacao ',i9)
 1100 format(' (CG) solver:'/
     . 5x,'Solver tol           = ',d20.6/
     . 5x,'tol * || b ||        = ',d20.6/
     . 5x,'Number of equations  = ',i20/
     . 5x,'nad                  = ',i20/
     . 5x,'Number of iterations = ',i20/
     . 5x,'x * Kx               = ',d20.10/
     . 5x,'|| x ||              = ',d20.10/
     . 5x,'CPU time (s)         = ',f20.2/)
 1200 format (' *** WARNING: No convergence reached after ',i9,
     .        ' iterations !',/)
 1300 format (' CG:',5x,'It',i7,5x,2d20.10)
 1400 format (' CG:',1x,'Residuo exato > conv ',1x,d20.10)
      end
c *********************************************************************  
c
c *********************************************************************  
      subroutine pcg(neq   ,nad    ,ia    ,ja
     .              ,ad    ,au     ,al    ,m        ,b      
     .              ,x     ,z      ,r     ,p     
     .              ,tol   ,maxit
     .              ,matvec,dot
     .              ,my_id ,neqf1i ,neqf2i,neq_doti,i_fmapi
     .              ,i_xfi ,i_rcvsi,i_dspli
     .              ,fprint,flog   ,fnew,mpi,nprcs)
c **********************************************************************
c * Data de criacao    : 00/00/0000                                    *
c * Data de modificaco : 24/04/2016                                    * 
c * ------------------------------------------------------------------ *   
c * Subroutine PCG : Solucao de sistemas de equacoes pelo metodo dos   *
c * gradientes conjugados com precondicionador diagonal para matrizes  *
c * simetricas.                                                        *
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ * 
c * neq      - numero de equacoes                                      *
c * nad      - numero de termos nao nulos                              *
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
c * i_dspli  - ponteiro extrutura da comunicacao    (MPI)              *                                                      *
c * fprint   - saida na tela                                           *
c * flog     - log do arquivo de saida                                 *
c * fnew     - .true.  -> x0 igual a zero                              *
c *            .false. -> x0 dado                                      *
c * mpi      - true|false                                              *
c * nprcs    - numero de processos mpi                                 *  
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
      integer neqf1i,neqf2i,neq_doti,nprcs
c ... ponteiros      
      integer*8 i_fmapi,i_xfi
      integer*8 i_rcvsi,i_dspli
c .....................................................................      
      integer neq,nequ,nad,maxit,i,j,jj
      integer ia(*),ja(*),my_id
      real*8  ad(*),au(*),al(*),b(*),m(*),x(*)
      real*8  r(*),z(*),p(*)
      real*8  dot,tol,conv,xkx,norm,d,di,alpha,beta,tmp,norm_b
      real*8  norm_r,norm_m_r
      real*8  time0,time
      real*8 dum1
      logical flog,fprint,fnew,mpi
c ...
      integer*8 flop_cg
      real*8  mflops,vmean
c .....................................................................
      external matvec,dot, flop_cg
c ======================================================================
      time0 = MPI_Wtime()
c ......................................................................
c
c ...
      do 5 i = 1, neq
        if(ad(i) .eq. 0.d0 ) then
          write(*,1000) i,ad(i)
          call stop_mef()
        endif 
   5  continue
c ......................................................................
c
c ... Chute inicial:
c
      if(fnew) then  
        do 10 i = 1, neq
          x(i)  = 0.d0
   10   continue
      endif  
c .......................................................................
c
c ... conv = tol * |(M-1)b|m = tol *(b,M-1b)
      do 15 i = 1, neq
         z(i) = b(i) * m(i)
   15 continue
      d      = dot(b,z,neq_doti)
      norm_b = dsqrt(dabs(d))  
      conv   = tol*dsqrt(dabs(d))
c .......................................................................
c  
c ... Ax0                                                            
      call matvec(neq,ia,ja,ia(neq+2),ja(nad+1),ad,al,au,al(nad+1) 
     .           ,x,z,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli
     .           ,dum1)
c .......................................................................
c
c ...
      do 100 i = 1, neq
c ... r0 = b - Ax0
         r(i) = b(i) - z(i)
c ... z0 = (M-1)r0
         z(i) = r(i) * m(i)
c ... p0 = r0
         p(i) = z(i)
  100 continue
c ... ( r(0),z(0) ) = ( r(0), (M-1)r0 )
      d    = dot(r,z,neq_doti)
c ......................................................................
      jj = 1
      do 230 j = 1, maxit
c ... z = Ap(j)
         call matvec(neq,ia,ja,ia(neq+2),ja(nad+1),ad,al,au,al(nad+1)
     .              ,p,z,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli
     .              ,dum1)
c .....................................................................
c
c ... alpha = ( r(j),z(j) ) / ( Ap(j), p(j) ))
         alpha = d / dot(z,p,neq_doti)
c .....................................................................
c
c ...
         do 210 i = 1, neq
c ... x(j+1) = x(j) + alpha*p
            x(i) = x(i) + alpha * p(i)
c ... r(j+1) = r(j) - alpha*Ap
            r(i) = r(i) - alpha * z(i)
c ... z  = (M-1)r0
            z(i) = r(i) * m(i)
  210    continue
c .....................................................................
c
c ... ( r(j+1),(M-1)r(j+1) ) = ( r(j+1),z )
         di   = dot(r,z,neq_doti) 
c ... beta = ( r(j+1),(M-1)r(j+1) ) / ( r(j),r(j) ) 
         beta = di / d
c .....................................................................
c
c ...
         do 220 i = 1, neq
c ... p(j+1) = (M-1)r(j+1) + beta*p(j) = z + beta*p(j)
            p(i) = z(i) + beta * p(i)
  220    continue
c .....................................................................
c
c ...
         d =  di
         if (dsqrt(dabs(d)) .lt. conv) goto 300
c ......................................................................
         if( jj .eq.1000) then
           jj = 0
           if(my_id .eq.0) write(*,1300),j,dsqrt(dabs(d)),conv 
         endif  
         jj = jj + 1
c ......................................................................
  230 continue
c ----------------------------------------------------------------------
      write(*,1200) maxit
      if(flog) write(10,1200) maxit
      call stop_mef()
  300 continue
c
c ... norm:  x*Kx
      call matvec(neq,ia,ja,ia(neq+2),ja(nad+1),ad,al,au,al(nad+1)
     .           ,x,z,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli,dum1)
      xkx = dot(x,z,neq_doti)
c ......................................................................
c
c ... norm-2 = || x ||
      norm = dsqrt(dot(x,x,neq_doti))
c ......................................................................
c
c ... r =M(-1)(b - Ax) (calculo do residuo explicito)
      do 310 i = 1, neq
        r(i) = b(i) - z(i)
        z(i) = r(i)*m(i)
  310 continue
      norm_m_r = dot(r,z,neq_doti)
      norm_m_r = dsqrt(dabs(norm_m_r))
      norm_r   = dot(r,r,neq_doti)
      norm_r   = dsqrt(norm_r)
      if( norm_m_r .gt. conv ) then
         if(my_id .eq.0 )then
           write(*,1400) norm_m_r,conv
         endif 
      endif
c ......................................................................
      time = MPI_Wtime()
      time = time-time0
c ......................................................................
c 
c ...    
      if(mpi) then
        call mpi_mean(vmean,time,nprcs) 
        time   = vmean        
      endif    
      mflops = (flop_cg(neq,nad,j,2,mpi)/1000000)/time  
c ......................................................................
c
c ...
      if(my_id .eq.0 .and. fprint )then
        if(mpi) then
          write(*,1110)tol,conv,j,xkx,norm,norm_r,norm_m_r,mflops,time
        else
          write(*,1100)tol,conv,neq,nad,j,xkx,norm,norm_r,norm_m_r
     .                ,mflops,time
        endif
      endif
c ......................................................................
c
c ... Controle de flops
      if(flog) then
        if(my_id.eq.0) then
          write(10,'(a,a,i9,a,d20.10,a,d20.10,a,d20.10,a,f20.2)')
     .       'PCG: ',' it ',j, ' x * Kx ',xkx,' ||x|| ',norm
     .      ,' tol ',tol,' time ',time
        endif
      endif
c ......................................................................
      return
c ======================================================================
 1000 format (//,5x,'SUBROTINA PCG:',/,5x,'Diagonal coefficient ' 
     . '- equation ',i9,d20.6)
 1100 format(' (PCG) solver:'/
     . 5x,'Solver tol           = ',d20.6/
     . 5x,'tol * ||b||m         = ',d20.6/
     . 5x,'Number of equations  = ',i20/
     . 5x,'nad                  = ',i20/
     . 5x,'Number of iterations = ',i20/
     . 5x,'x * Kx               = ',d20.10/
     . 5x,'|| x ||              = ',d20.10/
     . 5x,'|| b - Ax ||         = ',d20.10/
     . 5x,'|| b - Ax ||m        = ',d20.10/
     . 5x,'Mflops               = ',f20.2/
     . 5x,'CPU time (s)         = ',f20.2/)
 1110 format(' (PCG_MPI) solver:'/
     . 5x,'Solver tol           = ',d20.6/
     . 5x,'tol * ||b||m         = ',d20.6/
     . 5x,'Number of iterations = ',i20/
     . 5x,'x * Kx               = ',d20.10/
     . 5x,'|| x ||              = ',d20.10/
     . 5x,'|| b - Ax ||         = ',d20.10/
     . 5x,'|| b - Ax ||m        = ',d20.10/
     . 5x,'Mflops               = ',f20.2/
     . 5x,'CPU time (s)         = ',f20.2/)
 1200 format (' *** WARNING: No convergence reached after ',i9,
     .        ' iterations !',/)
 1300 format (' PCG:',5x,'It',i7,5x,2d20.10)
 1400 format (' PCG:',1x,'Explicit residual > tol * ||b||| :'
     .       ,1x,d20.10,1x,d20.10)
 1500 format ( 'PCG: ',5x,i7,5x,2es20.10)
      end
c *********************************************************************  
c
c *********************************************************************  
      subroutine gmres(neq,nad,ia,ja,ad,au,al,m,b,x,k,g,h,y,c,s,e,
     .              tol,maxit,matvec,dot,neqovlp,my_id,neqf1i,neqf2i,
     .              neq_doti,i_fmapi,i_xfi,i_rcvsi,i_dspli,flog)
c **********************************************************************
c * Data de criacao    : 00/00/0000                                    *
c * Data de modificaco : 24/04/2016                                    * 
c * ------------------------------------------------------------------ *   
c * GMRES: Solucao iterativa de sistemas simetricos e nao-simetricos   *
c *        pelo metodo GMRES com precondicionador diagonal.            *
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ *
c * neq      - numero de equacoes                                      *
c * nad      - numero de termos nao nulos                              *
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
c * my_id    - id do processo(MPI)                                     *
c * neqf1i   - numero de equacoes no buffer de recebimento (MPI)       *
c * neqf2i   - numero de equacoes no buffer de envio (MPI)             *
c * neq_doti - numero de equacoes no produto interno (MPI)             *
c * i_fmapi  - ponteiro para o mapa de comunicacao  (MPI)              *
c * i_xfi    - ponteiro para o buffer de valores    (MPI)              *
c * i_rcvsi  - ponteiro extrutura da comunicacao    (MPI)              *
c * i_dspli  - ponteiro extrutura da comunicacao    (MPI)              *                             
c * flog     - log do arquivo de saida                                 *
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ *
c * x(neq) - vetor solucao                                             *
c * b(*),ad(*),al(*),au(*) - inalterados                               *
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
      integer neqf1i,neqf2i,neq_doti
c ... ponteiros      
      integer*8 i_fmapi,i_xfi
      integer*8 i_rcvsi,i_dspli
c .....................................................................      
      integer neq,nad
      integer k,maxit,ia(*),ja(*),neqovlp,nit,i,j,jj,l,ni,ic,nad1
      real*8  ad(*),au(*),al(*),m(*),b(*),x(*)
      real*8  g(neqovlp,1:k+1),h(k+1,k),y(k),c(k),s(k),e(k+1),tol
      real*8  xkx,econv,norm,dot,r,aux1,aux2,beta
      real*8  time0,time
      real*8 dum1
      logical flog
      external matvec,dot
      integer my_id
c ......................................................................
      time0 = MPI_Wtime()
c ......................................................................
c
c.... Chute inicial:
c
      do 10 i = 1, neq
         x(i) = 0.d0
c ...    pre-condicionador diagonal:                  
         g(i,1) = b(i)*m(i)
   10 continue
c ----------------------------------------------------------------------
c
c ... Limite de convergencia:
c
      norm  = dsqrt(dot(g(1,1),g(1,1),neq_doti))
      econv = tol*norm
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
         call matvec(neq,ia,ja,ia(neq+2),ja(nad+1),ad,al,au,al(nad+1) 
     .              ,x,g(1,1)
     .              ,neqf1i,neqf2i,i_fmapi,i_xfi 
     .              ,i_rcvsi,i_dspli,dum1)
c
c ...... Residuo com precondicionador diagonal:
c
         do 200 i = 1, neq
            g(i,1) = (b(i) - g(i,1))*m(i)
  200    continue
c
c ...... Norma do residuo:
c
         e(1) = dsqrt(dot(g(1,1),g(1,1),neq_doti))
c
c ...... Normalizacao de g1:
c
         do 210 i = 1, neq
            g(i,1) = g(i,1)/e(1)
  210    continue
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
            call matvec(neq,ia,ja,ia(neq+2),ja(nad+1),ad,al,au
     .                 ,al(nad+1)
     .                 ,g(1,i),g(1,i+1)
     .                 ,neqf1i,neqf2i 
     .                 ,i_fmapi,i_xfi,i_rcvsi,i_dspli,dum1)
c
c ......... Precondicionador diagonal:
c
            do 300 j = 1, neq
                g(j,i+1) = g(j,i+1)*m(j)
  300       continue
c
c ......... Ortogonalizacao (Gram-Schmidt modificado):
c
            do 320 j = 1, i
               beta = dot(g(1,i+1),g(1,j),neq_doti)
               do 310 ic = 1, neq
                  g(ic,i+1) = g(ic,i+1) - beta * g(ic,j)
  310          continue
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
            do 330 ic = 1, neq
               g(ic,i+1) = g(ic,i+1)/norm
  330       continue
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
            y(i) = (y(i) + e(i)) / h(i,i)
  520    continue
c
c ...... Atualizacao de x:
c
         do 610 i = 1, neq
            do 600 j = 1, ni
               x(i) = x(i) + y(j) * g(i,j)
  600       continue
  610    continue
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
c         nii(l)=ni
         if (dabs(e(ni+1)) .le. econv) goto 1100
c ......................................................................
 1000 continue
c ......................................................................
 1100 continue
c
c ... Norma da solucao: x*Kx
c
      call matvec(neq,ia,ja,ia(neq+2),ja(nad+1),ad,al,au,
     .           al(nad+1),x     ,g(1,1)  ,neqf1i,neqf2i,
     .           i_fmapi,i_xfi,i_rcvsi,i_dspli,dum1)
      xkx = dot(x,g,neq_doti)
c ......................................................................
c
c ... norm-2 = || x ||
      norm = dsqrt(dot(x,x,neq_doti))
c ......................................................................
c
c ... r = b - Ax (calculo do residuo explicito)
      do 1200 i = 1, neq
        g(i,2) = b(i) - g(i,1)
 1200 continue
      aux1 = dot(g(1,2),g(1,2),neq_doti)
      aux1 = dsqrt(aux1)
      if( aux1 .gt. 3.16d0*econv ) then
         if(my_id .eq.0 )then
           write(*,2400) aux1,econv
         endif 
      endif
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
c ......................................................................
      if(my_id.eq.0)write(*,2000) tol,neq,l,nit,dabs(e(ni+1)),xkx,norm
     .                           ,time
c ......................................................................
c     Controle de flops
      if(flog) then
        if(my_id.eq.0) then
          write(10,'(a,a,i9,a,d20.10,a,d20.10,a,i9,a,d20.10,a,f20.2)')
     .         'GMRES: ',' it ',nit, ' x * Kx ',xkx,' ||x|| ',norm,
     .         ' nKylov ',k,' tol ',tol,' time ',time
        endif
      endif
c ......................................................................
      return
c ----------------------------------------------------------------------
 2000 format(' (GMRES) solver:'/
     . 5x,'Solver tol           = ',d20.6/
     . 5x,'Number of equations  = ',i20/
     . 5x,'Number of cycles     = ',i20/
     . 5x,'Number of iterations = ',i20/
     . 5x,'Norm                 = ',d20.10/
     . 5x,'x * Kx               = ',d20.10/
     . 5x,'|| x ||              = ',d20.10/
     . 5x,'CPU time (s)         = ',f20.2/)
 2100 format(' *** WARNING: no convergence reached for '
     .      ,i9,' cycles !',5x,i7,' nKylov',5x,' It ',i7/)
 2300 format (' GMRES:',5x,'cycles',i7,5x,'It',i7,5x,2d20.10)
 2400 format (' GMRES:',1x,'Residuo exato > 3.16d0*conv '
     .       ,1x,d20.10,1x,d20.10)
      end
c **********************************************************************
c
c **********************************************************************
      subroutine bicgstab(neq      ,nad   ,ia,ja 
     .                   ,ad      ,au    ,al,b  ,x   
     .                   ,t       ,v     ,r ,p  ,r0
     .                   ,tol     ,maxit  
     .                   ,matvec  ,dot    
     .                   ,my_id   ,neqf1i,neqf2i 
     .                   ,neq_doti,i_fmapi,i_xfi,i_rcvsi,i_dspli
     .                   ,fprint  ,flog   ,fnew)
c **********************************************************************
c * Data de criacao    : 15/04/2016                                    *
c * Data de modificaco : 24/04/2016                                    * 
c * ------------------------------------------------------------------ *   
c * BICGSTAB  : Solucao de sistemas de equacoes pelo metodo dos        * 
c * gradientes biconjugados para matrizes nao-simetricas.              *                                      
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ *  
c * neq      - numero de equacoes                                      *
c * nad      - numero de termos nao nulos                              *
c * ia(*)    - ponteiro do formato CSR                                 *
c * ja(*)    - ponteiro das colunas no formato CSR                     *
c * ad(neq)  - diagonal da matriz A                                    *
c * au(*)    - parte triangular superior de A                          *
c * al(*)    - parte triangular inferior de A                          *
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
c * my_id    - id do processo(MPI)                                     *
c * neqf1i   - numero de equacoes no buffer de recebimento (MPI)       *
c * neqf2i   - numero de equacoes no buffer de envio (MPI)             *
c * neq_doti - numero de equacoes no produto interno (MPI)             *
c * i_fmapi  - ponteiro para o mapa de comunicacao  (MPI)              *
c * i_xfi    - ponteiro para o buffer de valores    (MPI)              *
c * i_rcvsi  - ponteiro extrutura da comunicacao    (MPI)              *
c * i_dspli  - ponteiro extrutura da comunicacao    (MPI)              *                                                     *
c * fprint   - saida na tela                                           *
c * flog     - log do arquivo de saida                                 *
c * fnew     - .true.  x0 igual a zero                                 *
c *            .false. x0 dado                                         *
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ *
c * x(neq) - vetor solucao                                             *
c * b(neq) - modificado                                                *
c * ad(*),al(*),au(*) e b - inalterados                                *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ **
c **********************************************************************
      implicit none
      include 'mpif.h'
      integer neqf1i,neqf2i,neq_doti
c ... ponteiros      
      integer*8 i_fmapi,i_xfi
      integer*8 i_rcvsi,i_dspli
c .....................................................................   
      integer neq,nad
      integer maxit,i,j,jj,k
      integer ia(*),ja(*),my_id
      real*8  ad(*),au(*),al(*),x(*),b(*)
      real*8  r(*),p(*),t(*),v(*),r0(*)
      real*8  dot,tol,conv,xkx,norm,d,alpha,beta,rr0,w,tmp
      real*8  time0,time
      real*8  dum1 
      logical flog,fprint,fnew
      external matvec,dot
c ======================================================================
      time0 = MPI_Wtime()
c ......................................................................
c     if(my_id.eq.0) print *, 'nad :',nad
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
c
c ... Chute inicial:
c
      if(fnew) then  
        do 10 i = 1, neq
          x(i)  = 0.d0
   10   continue
      endif  
c .......................................................................
c
c ... conv = tol * |b|
      d    = dot(b,b,neq_doti)
      conv = tol*dsqrt(dabs(d))
c .......................................................................
c
c ...
      call matvec(neq,ia,ja,ia(neq+2),ja(nad+1),ad,al,au,al(nad+1) 
     .           ,x,p
     .           ,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli,dum1)
c .......................................................................
c
c ...
      do 100 i = 1, neq
c ... r0 = b - Ax0
         r0(i) = b(i) - p(i)
c ... r = r0
         r(i)  = r0(i)
c ... p = r0
         p(i)  = r0(i)
  100 continue
c .......................................................................
c
c ...
      jj = 1
      do 230 j = 1, maxit
c ... v = Ap(j)
         call matvec(neq,ia,ja,ia(neq+2),ja(nad+1),ad,al,au,al(nad+1)
     .              ,p,v
     .              ,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli,dum1)
c .......................................................................
c
c ... alpha = ( r(j),r0 ) / ( Ap(j), r0 ))
         rr0   = dot(r,r0,neq_doti)
         alpha = rr0/dot(v,r0,neq_doti)
c .......................................................................
c
c ...
         do 210 i = 1, neq
c ... x(j+1) = x(j) + alpha*p
            x(i) = x(i) + alpha * p(i)
c ... s(j)   = r(j) - alpha*Ap(j)
            r(i) = r(i) - alpha * v(i)
  210    continue
c ........................................................................
c
c ... (s,s)
         d = dot(r,r,neq_doti)
c ...
         if (dsqrt(dabs(d)) .lt. conv) goto 300
c ........................................................................
c
c ... t = Ar(j)
         call matvec(neq,ia,ja,ia(neq+2),ja(nad+1),ad,al,au,al(nad+1) 
     .              ,r,t 
     .              ,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli,dum1)
c ........................................................................
c
c ... w = ( Ar(j),r(j) ) / ( Ar(j), Ar(j) ))
         w = dot(t,r,neq_doti) / dot(t,t,neq_doti)
c ........................................................................
c
c ... 
         do 220 i = 1, neq
c ... x(j+1) = x(j) + w*r(j)
            x(i) = x(i) + w*r(i)
c ... r(j+1) = s(j) - w*As(j)
            r(i) = r(i) - w*t(i)
  220    continue
c .......................................................................
c
c ... (r,r)
         d = dot(r,r,neq_doti)
c ...
         if (dsqrt(dabs(d)) .lt. conv) goto 300
c .......................................................................
c
c ... beta = ( r(j+1),r0 ) / ( r(j), r0 )) * (alpha/w) 
         beta = (dot(r,b,neq_doti) / rr0)*(alpha/w)
c .......................................................................
c
c ...
         do 225 i = 1, neq
             p(i) = r(i) + beta*(p(i)-w*v(i))
  225    continue
c .......................................................................
c
c ...
         if( jj .eq. 500) then
           jj = 0
           if(my_id.eq.0) write(*,1300),j,dsqrt(dabs(d)),conv 
         endif  
         jj = jj + 1
c ......................................................................
  230 continue
c ......................................................................
c
c ...
      if(my_id.eq.0) then
        write(*,1200) maxit
        if(fLog) write(10,1200) maxit
      endif
      call stop_mef()
  300 continue
c ......................................................................
c
c ... produto:  x*Kx
c
      call matvec(neq,ia,ja,ia(neq+2),ja(nad+1),ad,al,au,al(nad+1) 
     .           ,x,p 
     .           ,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli,dum1)
      xkx =  dot(x,p,neq_doti)
c .......................................................................
c
c ... norm-2 = || x ||
      norm = dsqrt(dot(x,x,neq_doti))
c ......................................................................
c
c ... r = b - Ax (calculo do residuo explicito)
      do 310 i = 1, neq
        r(i) = b(i) - p(i)
  310 continue
      tmp  = dot(r,r,neq_doti)
      tmp = dsqrt(tmp)
      if( tmp .gt. conv ) then
        write(*,1400) tmp
      endif
c ......................................................................
c
c ...
      time = MPI_Wtime()
      time = time-time0
c .......................................................................
c
c ...
      if(my_id .eq.0 .and. fprint )then
        write(*,1100)tol,conv,neq,nad,j,xkx,norm,time
      endif
c .......................................................................
c
c ... Controle de flops
      if(flog) then
        if(my_id.eq.0) then
          write(10,'(a,a,i9,a,d20.10,a,d20.10,a,d20.10,a,f20.2)')
     .       "BICGSTAB: "," it ",j, " x * Kx ",xkx," ||x|| ",norm
     .      ," tol ",tol," time ",time
        endif
      endif
c ......................................................................
      return
c ======================================================================
 1000 format (/,1x,'SUBROTINA BICCSTAB:',/,1x
     .       ,'Coeficiente da diagonal nulo ',i9)
 1100 format(' (BICGSTAB) solver:'/
     . 5x,'Solver tol           = ',d20.6/
     . 5x,'tol * || b ||        = ',d20.6/
     . 5x,'Number of equations  = ',i20/
     . 5x,'nad                  = ',i20/
     . 5x,'Number of iterations = ',i20/
     . 5x,'x * Kx               = ',d20.10/
     . 5x,'|| x ||              = ',d20.10/
     . 5x,'CPU time (s)         = ',f20.2/)
 1200 format (' *** WARNING: No convergence reached after ',i9,
     .        ' iterations !',/)
 1300 format (' BICGSTAB:',5x,'It',i7,5x,2d20.10)
 1400 format (' BICGSTAB:',1x,'Residuo exato > conv ',1x,d20.10)
      end
c *********************************************************************
c
c **********************************************************************
      subroutine pbicgstab(neq     ,nad   ,ia ,ja 
     .                    ,ad      ,au    ,al ,m  ,b ,x   
     .                    ,t       ,v     ,r  ,p  ,z ,r0 
     .                    ,tol     ,maxit  
     .                    ,matvec  ,dot    
     .                    ,my_id   ,neqf1i,neqf2i 
     .                    ,neq_doti,i_fmapi,i_xfi,i_rcvsi,i_dspli
     .                    ,fprint  ,flog   ,fnew)
c **********************************************************************
c * Data de criacao    : 00/00/0000                                    *
c * Data de modificaco : 24/04/2016                                    * 
c * ------------------------------------------------------------------ *   
c * PBICGSTAB : Solucao de sistemas de equacoes pelo metodo dos        * 
c * gradientes biconjugados com precondicionador diagonal para         *
c * matrizes nao-simetricas.                                           *
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ *
c * neq      - numero de equacoes                                      *
c * nad      - numero de termos nao nulos                              *
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
c * r0(neq)  - arranjo local de trabalho                               *
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
c * i_dspli  - ponteiro extrutura da comunicacao    (MPI)              *                                                         *
c * fprint   - saida na tela                                           *
c * flog     - log do arquivo de saida                                 *
c * fnew     - .true.  x0 igual a zero                                 *
c *            .false. x0 dado                                         *
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ *
c * x(neq) - vetor solucao                                             *
c * ad(*),al(*),au(*) e b(*)  inalterados                              *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * A(M-1)y=b precondicionador direita                                 *
c * ------------------------------------------------------------------ *
c **********************************************************************
      implicit none
      include 'mpif.h'
      integer neqf1i,neqf2i,neq_doti
c ... ponteiros      
      integer*8 i_fmapi,i_xfi
      integer*8 i_rcvsi,i_dspli
c .....................................................................   
      integer neq,nad
      integer maxit,i,j,jj,k
      integer ia(*),ja(*),my_id
      real*8  ad(*),au(*),al(*),m(*),x(*),b(*)
      real*8  r(*),p(*),t(*),v(*),z(*),r0(*)
      real*8  dot,tol,conv,d,alpha,beta,rr0,w,xkx,norm,tmp
      real*8  time0,time
      real*8  dum1 
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
c
c ... Chute inicial:
c 
      do 10 i = 1, neq
         x(i) = 0.d0
   10 continue
c .......................................................................
c
c ... conv = tol * |b|
      d    = dot(b,b,neq_doti)
      conv = tol*dsqrt(dabs(d))
c .......................................................................
c
c ... Ax0
      call matvec(neq,ia,ja,ia(neq+2),ja(nad+1),ad,al,au,al(nad+1) 
     .           ,x,z 
     .           ,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli,dum1)
c .......................................................................
c
c ...
      do 100 i = 1, neq
c ... r0 = b - Ax0
         r0(i) = b(i) - z(i)
c ... r = r0
         p(i)  = r0(i)
c ... p = r0
         r(i)  = r0(i)
c ... z = M(-1)p
         z(i)  = p(i)*m(i) 
  100 continue
c .......................................................................
c
c ...
      jj = 1
      do 230 j = 1, maxit
c ... v = Az(j) = AM(-1)p(j)
         call matvec(neq,ia,ja,ia(neq+2),ja(nad+1),ad,al,au,al(nad+1) 
     .              ,z,v 
     .              ,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli,dum1)
c .......................................................................
c
c ... alpha = ( r(j),r0 ) / ( AM(-1)p(j), r0 ))
         rr0   = dot(r,r0,neq_doti)
         alpha = rr0/dot(v,r0,neq_doti)
c .......................................................................
c
c ...
         do 210 i = 1, neq
c ... x(j+1) = x(j) + alpha*M(-1)p
            x(i) = x(i) + alpha * z(i)
c ... s(j)   = r(j) - alpha*AM(-1)p(j)
            r(i) = r(i) - alpha * v(i)
c ... z = M(-1)s
            z(i) = r(i) * m(i)
  210    continue
c ........................................................................
c
c ... (s,s)
         d = dot(r,r,neq_doti)
c ...
         if (dsqrt(dabs(d)) .lt. conv) goto 300
c ........................................................................
c
c ... t = Az = AM(-1)s(j)
         call matvec(neq,ia,ja,ia(neq+2),ja(nad+1),ad,al,au,al(nad+1) 
     .              ,z,t 
     .              ,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli,dum1)
c ........................................................................
c
c ... w = ( AM(-1)s(j) ,s(j) ) / ( AM(-1)s(j), AM(-1)s(j) )
         w = dot(t,r,neq_doti) / dot(t,t,neq_doti)
c ........................................................................
c
c ... 
         do 220 i = 1, neq
c ... x(j+1) = x(j) + w*M(-1)s
            x(i) = x(i) + w*z(i)
c ... r(j+1) = s(j) - w*AM(-1)s
            r(i) = r(i) - w*t(i)
  220    continue
c ........................................................................
c
c ... (r,r)
         d = dot(r,r,neq_doti)
c ...
         if (dsqrt(dabs(d)) .lt. conv) goto 300
c ........................................................................
c
c ... beta = ( r(j+1),r0 ) / ( r(j), r0 )) * (alpha/w) 
         beta = (dot(r,r0,neq_doti) / rr0)*(alpha/w)
c .......................................................................
c
c ...
         do 225 i = 1, neq
c ... p(j+1) = r(i) + beta*(p(j)-w*v(i))
             p(i) = r(i) + beta*(p(i)-w*v(i))
c ... z = M(-1)p
             z(i) = p(i)*m(i)
  225    continue
c .......................................................................
c
c ...
         if( jj .eq. 500) then
           jj = 0
           if(my_id.eq.0) write(*,1300),j,dsqrt(dabs(d)),conv 
         endif  
         jj = jj + 1
c ......................................................................
  230 continue
c ......................................................................
c
c ...
      if(my_id.eq.0) then
        write(*,1200) maxit
        if(fLog) write(10,1200) maxit
      endif
      call stop_mef()
  300 continue
c ......................................................................
c
c ... produto:  x*Kx
c
      call matvec(neq,ia,ja,ia(neq+2),ja(nad+1),ad,al,au,al(nad+1),
     .            x,z,
     .            neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli,dum1)
      xkx =  dot(x,z,neq_doti)
c .......................................................................
c
c ... norm-2 = || x ||
      norm = dsqrt(dot(x,x,neq_doti))
c ......................................................................
c
c ... r = b - Ax (calculo do residuo explicito)
      do 310 i = 1, neq
        r(i) = b(i) - z(i)
  310 continue
      tmp  = dot(r,r,neq_doti)
      tmp = dsqrt(tmp)
      if( tmp .gt. 3.16d0*conv ) then
        write(*,1400) tmp,conv
      endif
c ......................................................................
c
c ...
      time = MPI_Wtime()
      time = time-time0
c .......................................................................
c
c ...
      if(my_id .eq.0 .and. fprint )then
        write(*,1100)tol,conv,neq,nad,j,xkx,norm,time
      endif
c .......................................................................
c
c ... Controle de flops
      if(flog) then
        if(my_id.eq.0) then
          write(10,'(a,a,i9,a,d20.10,a,d20.10,a,d20.10,a,f20.2)')
     .       "BICGSTAB: "," it ",j, " x * Kx ",xkx," ||x|| ",norm
     .      ," tol ",tol," time ",time
        endif
      endif
c ......................................................................
      return
c ======================================================================
 1000 format (/,1x,'SUBROTINA PBICCSTAB:',/,1x
     .       ,'Coeficiente da diagonal nulo ',i9)
 1100 format(' (PBICGSTAB) solver:'/
     . 5x,'Solver tol           = ',d20.6/
     . 5x,'tol * || b ||        = ',d20.6/
     . 5x,'Number of equations  = ',i20/
     . 5x,'nad                  = ',i20/
     . 5x,'Number of iterations = ',i20/
     . 5x,'x * Kx               = ',d20.10/
     . 5x,'|| x ||              = ',d20.10/
     . 5x,'CPU time (s)         = ',f20.2/)
 1200 format (' *** WARNING: No convergence reached after ',i9,
     .        ' iterations !',/)
 1300 format (' PBICGSTAB:',5x,'It',i7,5x,2d20.10)
 1400 format (' PBICGSTAB:',1x,'Residuo exato > 3.16d0*conv '
     .       ,1x,d20.10,1x,d20.10)
      end
c *********************************************************************
c
c **********************************************************************
      subroutine icbicgstab(neq     ,nad   ,ia ,ja 
     .                     ,ad      ,au    ,al ,m  ,b     ,x   
     .                     ,t       ,v     ,r  ,p  ,z     ,r0
     .                     ,tol     ,maxit  
     .                     ,matvec  ,dot    
     .                     ,my_id   ,neqf1i,neqf2i 
     .                     ,neq_doti,i_fmapi,i_xfi,i_rcvsi,i_dspli
     .                     ,fprint  ,flog   ,fnew)
c **********************************************************************
c * Data de criacao    : 15/04/2016                                    *
c * Data de modificaco : 24/04/2016                                    * 
c * ------------------------------------------------------------------ *   
c * PBICGSTAB : Solucao de sistemas de equacoes pelo metodo dos        * 
c * gradientes biconjugados com precondicionador diagonal para         *
c * matrizes nao-simetricas.                                           *
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ *
c * neq      - numero de equacoes                                      *
c * nad      - numero de termos nao nulos                              *
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
c * r0(neq)  - arranjo local de trabalho                               *
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
c * i_dspli  - ponteiro extrutura da comunicacao    (MPI)              *                                                         *
c * fprint   - saida na tela                                           *
c * flog     - log do arquivo de saida                                 *
c * fnew     - .true.  x0 igual a zero                                 *
c *            .false. x0 dado                                         *
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ *
c * x(neq) - vetor solucao                                             *
c * b(neq) - modificado                                                *
c * ad(*),al(*),au(*) - inalterados                                    *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * A(M-1)y=b precondicionador direita                                 *
c * ------------------------------------------------------------------ *
c **********************************************************************
      implicit none
      include 'mpif.h'
      integer neqf1i,neqf2i,neq_doti
c ... ponteiros      
      integer*8 i_fmapi,i_xfi
      integer*8 i_rcvsi,i_dspli
c .....................................................................   
      integer neq,nad
      integer maxit,i,j,jj,k
      integer ia(*),ja(*),my_id
      real*8  ad(*),au(*),al(*),m(*),x(*),b(*)
      real*8  r(*),p(*),t(*),v(*),z(*),r0(*)
      real*8  dot,tol,conv,xkx,norm,d,alpha,beta,rr0,w,tmp
      real*8  time0,time
      real*8  dum1 
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
c
c ... Chute inicial:
c 
      do 10 i = 1, neq
         x(i) = 0.d0
   10 continue
c .......................................................................
c
c ... conv = tol * |b|
      d    = dot(b,b,neq_doti)
      conv = tol*dsqrt(dabs(d))
c .......................................................................
c
c ... Ax0
      call matvec(neq,ia,ja,ia(neq+2),ja(nad+1),ad,al,au,al(nad+1),
     .            x,z,
     .            neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli,dum1)
c .......................................................................
c
c ...
      do 100 i = 1, neq
c ... r0 = b - Ax0
         r0(i) = b(i) - z(i)
c ... r = r0
         p(i) = r0(i)
c ... p = r0
         r(i) = r0(i)
  100 continue
c .......................................................................
c
c ... Mz=p  
      call ildlt_solv(neq,ia,ja,m,m(neq+1),p,z)
c .......................................................................      
c
c ...
      jj = 1
      do 230 j = 1, maxit
c ... v = Az(j) = AM(-1)p(j)
         call matvec(neq,ia,ja,ia(neq+2),ja(nad+1),ad,al,au,al(nad+1),
     .               z,v,
     .               neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli,dum1)
c .......................................................................
c
c ... alpha = ( r(j),r0 ) / ( AM(-1)p(j), r0 ))
         rr0   = dot(r,r0,neq_doti)
         alpha = rr0/dot(v,r0,neq_doti)
c .......................................................................
c
c ...
         do 210 i = 1, neq
c ... x(j+1) = x(j) + alpha*M(-1)p
            x(i) = x(i) + alpha * z(i)
c ... s(j)   = r(j) - alpha*AM(-1)p(j)
            r(i) = r(i) - alpha * v(i)
  210    continue
c ........................................................................
c
c
c ... (s,s)
         d = dot(r,r,neq_doti)
c ...
         if (dsqrt(dabs(d)) .lt. conv) goto 300
c ........................................................................
c
c ... Mz=r  
         call ildlt_solv(neq,ia,ja,m,m(neq+1),r,z)
c .......................................................................
c
c ... (r,r)
         d = dot(r,r,neq_doti)
c ...
         if (dsqrt(dabs(d)) .lt. conv) goto 300
c ........................................................................
c
c ...  t = Az = AM(-1)s(j)
         call matvec(neq,ia,ja,ia(neq+2),ja(nad+1),ad,al,au,al(nad+1),
     .               z,t,
     .               neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli,dum1)
c ........................................................................
c
c ... w = ( AM(-1)s(j) ,s(j) ) / ( AM(-1)s(j), AM(-1)s(j) )
         w = dot(t,r,neq_doti) / dot(t,t,neq_doti)
c ........................................................................
c
c ... 
         do 220 i = 1, neq
c ... x(j+1) = x(j) + w*M(-1)s
            x(i) = x(i) + w*z(i)
c ... r(j+1) = s(j) - w*AM(-1)s
            r(i) = r(i) - w*t(i)
  220    continue
c ........................................................................
c
c ... (r,r)
         d = dot(r,r,neq_doti)
c ...
         if (dsqrt(dabs(d)) .lt. conv) goto 300
c ........................................................................
c
c ... beta = ( r(j+1),r0 ) / ( r(j), r0 )) * (alpha/w) 
         beta = (dot(r,r0,neq_doti) / rr0)*(alpha/w)
c .......................................................................
c
c ...
         do 225 i = 1, neq
c ... p(j+1) = r(i) + beta*(p(j)-w*v(i))
             p(i) = r(i) + beta*(p(i)-w*v(i))
  225    continue
c .......................................................................
c
c ... Mz=p  
         call ildlt_solv(neq,ia,ja,m,m(neq+1),p,z)
c .......................................................................
c
c ...
         if( jj .eq. 500) then
           jj = 0
           if(my_id.eq.0) write(*,1300),j,dsqrt(dabs(d)),conv 
         endif  
         jj = jj + 1
c ......................................................................
  230 continue
c ......................................................................
      if(my_id.eq.0) then
        write(*,1200) maxit
        if(fLog) write(10,1200) maxit
      endif
      call stop_mef()
  300 continue
c
c ... produto:  x*Kx
c
      call matvec(neq,ia,ja,ia(neq+2),ja(nad+1),ad,al,au,al(nad+1),
     .            x,z,
     .            neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli,dum1)
       xkx = dot(x,z,neq_doti)
c .......................................................................
c
c ... norm-2 = || x ||
      norm = dsqrt(dot(x,x,neq_doti))
c ......................................................................
c
c ... r = b - Ax (calculo do residuo explicito)
      do 310 i = 1, neq
        r(i) = b(i) - z(i)
  310 continue
      tmp  = dot(r,r,neq_doti)
      tmp = dsqrt(tmp)
      if( tmp .gt. conv ) then
        write(*,1400) tmp
      endif
c ......................................................................
c
c ...
      time = MPI_Wtime()
      time = time-time0
c .......................................................................
c
c ...
      if(my_id .eq.0 .and. fprint )then
        write(*,1100)tol,conv,neq,nad,j,xkx,norm,time
      endif
c .......................................................................
c
c ... Controle de flops
      if(flog) then
        if(my_id.eq.0) then
          write(10,'(a,a,i9,a,d20.10,a,d20.10,a,d20.10,a,f20.2)')
     .       "ICBICGSTAB: "," it ",j, " x * Kx ",xkx," ||x|| ",norm
     .      ," tol ",tol," time ",time
        endif

      endif
c ......................................................................
      return
c ======================================================================
 1000 format (/,1x,'SUBROTINA ICBICCSTAB:',/,1x
     .       ,'Coeficiente da diagonal nulo ',i9)
 1100 format(' (ICBICGSTAB) solver:'/
     . 5x,'Solver tol           = ',d20.6/
     . 5x,'tol * || b ||        = ',d20.6/
     . 5x,'Number of equations  = ',i20/
     . 5x,'nad                  = ',i20/
     . 5x,'Number of iterations = ',i20/
     . 5x,'x * Kx               = ',d20.10/
     . 5x,'|| x ||              = ',d20.10/
     . 5x,'CPU time (s)         = ',f20.2/)
 1200 format (' *** WARNING: No convergence reached after ',i9,
     .        ' iterations !',/)
 1300 format (' ICBICCSTAB:',5x,'It',i7,5x,2d20.10)
 1400 format (' ICBICCSTAB:',1x,'Residuo exato > conv ',1x,d20.10)
      end
c *********************************************************************
c
c *********************************************************************
      subroutine iccg(neq   ,nad    ,ia      ,ja
     .               ,ad    ,au     ,al      ,m       ,b       
     .               ,x     ,z      ,r       ,p   
     .               ,tol   ,maxit
     .               ,matvec,dot    ,triasolv
     .               ,my_id ,neqf1i ,neqf2i  ,neq_doti,i_fmapi
     .               ,i_xfi ,i_rcvsi,i_dspli
     .               ,fprint,flog   ,fnew)
c **********************************************************************
c * Data de criacao    : 10/04/2016                                    *
c * Data de modificaco : 24/04/2016                                    * 
c * ------------------------------------------------------------------ *   
c * IcCG : Solucao de sistemas de equacoes pelo metodo dos gradientes   *
c * conjugados com precondicionador incompleto                         *
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ * 
c * neq      - numero de equacoes                                      *
c * nad      - numero de termos nao nulos                              *
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
c * triasolv - nome da funcao externa para o o solver triangular       *
c * my_id    - id do processo(MPI)                                     *
c * neqf1i   - numero de equacoes no buffer de recebimento (MPI)       *
c * neqf2i   - numero de equacoes no buffer de envio (MPI)             *
c * neq_doti - numero de equacoes no produto interno (MPI)             *
c * i_fmapi  - ponteiro para o mapa de comunicacao  (MPI)              *
c * i_xfi    - ponteiro para o buffer de valores    (MPI)              *
c * i_rcvsi  - ponteiro extrutura da comunicacao    (MPI)              *
c * i_dspli  - ponteiro extrutura da comunicacao    (MPI)              *                                                        *
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
c * Arranjos jat,iat e kat são utilizados na retrosubstituizao do      *
c * solver iLDLt                                                       * 
c **********************************************************************
      implicit none
      include 'mpif.h'
      integer neqf1i,neqf2i,neq_doti
c ... ponteiros      
      integer*8 i_fmapi,i_xfi
      integer*8 i_rcvsi,i_dspli
c .....................................................................      
      integer neq,nad,maxit,i,j,jj
      integer ia(*),ja(*),my_id
      real*8  ad(*),au(*),al(*),x(*),b(*),r(*),z(*),p(*)
      real*8  m(*)
      real*8  dot,ddot,tol,conv,xkx,norm,d,di,alpha,beta,tmp
      real*8  time0,time
      real*8 dum1
      logical flog,fnew,fprint
      external matvec,dot,triasolv
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
c
c ... Chute inicial:
c
      if(fnew) then  
        do 10 i = 1, neq
          x(i)  = 0.d0
   10   continue
      endif  
c ......................................................................
c
c ... conv = tol * |(M-1)b|
      call triasolv(neq,ia,ja,m,m(neq+1),b,z)
      d    = dot(z,z,neq_doti)
      conv = tol*dsqrt(dabs(d))
c .......................................................................
c  
c ... Ax0
      call matvec(neq,ia,ja,ia(neq+2),ja(nad+1),ad,al,au,al(nad+1)  
     .           ,x,z 
     .           ,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli,dum1)
c .......................................................................      
c
      do 100 i = 1, neq
c ... r0 = b - Ax0
         r(i) = b(i) - z(i)
  100 continue
c .......................................................................      
c
c ... Mz=r  
      call triasolv(neq,ia,ja,m,m(neq+1),r,z)
c .......................................................................      
c
c ...
      do 105 i = 1, neq
c ... p0 = r0
         p(i) = z(i)
  105 continue
c .......................................................................      
c
c ... ( r(0),z(0) ) = ( r(0), (M-1)r0 )
      d    = dot(r,z,neq_doti) 
c .......................................................................
      jj = 1
      do 230 j = 1, maxit
c ... z = Ap(j)
         call matvec(neq,ia,ja,ia(neq+2),ja(nad+1),ad,al,au,al(nad+1) 
     .              ,p,z,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli 
     .              ,dum1)
c .....................................................................
c
c ... alpha = ( r(j),z(j) ) / ( Ap(j), p(j) ))
         alpha = d / dot(z,p,neq_doti)
c .....................................................................
c
c ...
         do 210 i = 1, neq
c ... x(j+1) = x(j) + alpha*p
            x(i) = x(i) + alpha * p(i)
c ... r(j+1) = r(j) - alpha*Ap
            r(i) = r(i) - alpha * z(i)
  210    continue
c ......................................................................
c
c ... Mz=r  
         call triasolv(neq,ia,ja,m,m(neq+1),r,z)
c .......................................................................

c ... ( r(j+1),(M-1)r(j+1) ) = ( r(j+1),z )
         di   = dot(r,z,neq_doti) 
c ... beta = ( r(j+1),(M-1)r(j+1) ) / ( r(j),r(j) ) 
         beta = di / d
c .....................................................................
c
c ...         
         do 220 i = 1, neq
            p(i) = z(i) + beta * p(i)
  220    continue
c .....................................................................
c
c ...
         d = di           
         if (dsqrt(dabs(d)) .lt. conv) goto 300
c ......................................................................
         if( jj .eq.500) then
           jj = 0
           write(*,1300),j,dsqrt(dabs(d)),conv 
         endif  
         jj = jj + 1
c ......................................................................
  230 continue
c ----------------------------------------------------------------------
      write(*,1200) maxit
      if(flog) write(10,1200) maxit
      call stop_mef()
  300 continue
c
c ... Energy norm: x*Kx
c
      call matvec(neq,ia,ja,ia(neq+2),ja(nad+1),ad,al,au,al(nad+1) 
     .          ,x,z 
     .          ,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli,dum1)
      xkx    = dot(x,z,neq_doti)
c ......................................................................
c
c ... norm-2 = || x ||
      norm = dsqrt(dot(x,x,neq_doti))
c ......................................................................
c
c ... r = b - Ax (calculo do residuo explicito)
      do 310 i = 1, neq
        r(i) = b(i) - z(i)
  310 continue
      tmp  = dot(r,r,neq_doti)
      tmp = dsqrt(tmp)
      if( tmp .gt. 3.16d0*conv ) then
        write(*,1400) tmp,conv
      endif
c ......................................................................
      time = MPI_Wtime()
      time = time-time0
c ----------------------------------------------------------------------
      if(my_id .eq.0 .and. fprint )then
        write(*,1100)tol,conv,neq,nad,j,xkx,norm,time
      endif
c ......................................................................
c
c ... Controle de flops
      if(flog) then
        if(my_id.eq.0) then
          write(10,'(a,a,i9,a,d20.10,a,d20.10,a,d20.10,a,f20.2)')
     .       "ICCG: "," it ",j, " x * Kx ",xkx," ||x|| ",norm
     .      ," tol ",tol," time ",time
        endif
      endif
c ......................................................................
      return
c ======================================================================
 1000 format (//,5x,'SUBROTINA ICCG:',/,5x,'Coeficiente da diagonal nulo
     . - equacao ',i9)
 1100 format(' (ICCG) solver:'/
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
 1300 format ( 'ICCG: ',5x,'It',i7,5x,2d20.10)
 1400 format (' ICCG:',1x,'Residuo exato > 3.16*conv ',1x,d20.10
     .       ,1x,d20.10)
      end
c **********************************************************************
      real*8 function smachn()
c **********************************************************************
c *                                                                    *
c *   SMACHN: calcula a precisao da maquina para real*8                *
c *                                                                    *
c **********************************************************************
      implicit none
      smachn = 1.0d0
100   smachn = smachn*0.5d0
      if(smachn+1.d0 .ne. 1.d0) go to 100
c     smachn = 2.0d0*smachn
      return
      end
c ***********************************************************************
c     subroutine bicgstab(neq,ia,ja,ad,au,al,m,b,x,y,z,p,r,s,tol,maxit,
c    .                    matvec,dot,my_id,neqf1i,neqf2i,neq_doti,
c    .                    i_fmapi,i_xfi,i_rcvsi,i_dspli)
c **********************************************************************
c *                                                                    *
c *   Subroutine BICGSTAB                                              *
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
c *                                                                    *
c *   Parametros de saida:                                             *
c *                                                                    *
c *   x(neq) - vetor solucao                                           *
c *   b(neq) - modificado                                              *
c *   ad(*),al(*),au(*) - modificado                                   *
c *                                                                    *
c **********************************************************************
c     implicit none
c     include 'mpif.h'
c     integer neqf1i,neqf2i,neq_doti
c ... ponteiros      
c     integer*8 i_fmapi,i_xfi
c     integer*8 i_rcvsi,i_dspli
c .....................................................................
c     integer neq,maxit,nad,i,j,k
c     integer ia(*),ja(*),my_id
c     real*8  ad(*),au(*),al(*),m(*),x(*),r(*),p(*),b(*),y(*),z(*),s(*)
c     real*8  dot,ddot,tol,conv,energy,d,alpha,beta,rr0,w
c     real*8  time0,time
c     real*8 dum1
c     external matvec,dot
c ======================================================================
c     time0 = MPI_Wtime() 
c ......................................................................
c     nad = ia(neq+1)-1
c     if(my_id.eq.0)print *, 'nad :',nad
c ......................................................................
c
c ... Chute inicial:
c
c     do 10 i = 1, neq
c        x(i) = 0.d0
c ...... pre-condicionador diagonal:         
c        b(i)  = b(i)/m(i)
c        ad(i) = ad(i)/m(i)
c        do 5 k = ia(i), ia(i+1)-1
c           j = ja(k)
c           al(k) = al(k) / m(i)
c           au(k) = au(k) / m(j)
c  5     continue      
c 10  continue
c ----------------------------------------------------------------------
c     call matvec(neq,ia,ja,ia(neq+2),ja(nad+1),ad,al,au,al(nad+1),x,r,
c    .            neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli,dum1)    
c     do 100 i = 1, neq
c        r(i) = b(i) - r(i)
c        p(i) = r(i)
c        b(i) = r(i) 
c 100 continue
c     d    = dot(r(1),r(1),neq)
c     conv = tol*dsqrt(dabs(d))
c ----------------------------------------------------------------------
c     do 230 j = 1, maxit
c        call matvec(neq,ia,ja,ia(neq+2),ja(nad+1),ad,al,au,al(nad+1),
c    .               p,z,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli,
c    .               dum1)
c        rr0 = dot(r,b,neq)
c        alpha = rr0/dot(z,b,neq)
c         do 210 i = 1, neq
c           s(i) = r(i) - alpha * z(i)
c 210    continue
c        call matvec(neq,ia,ja,ia(neq+2),ja(nad+1),ad,al,au,al(nad+1),
c    .               s,y, neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli,
c    .               dum1)
c        w = dot(y,s,neq) / dot(y,y,neq)
c        do 220 i = 1, neq
c           x(i) = x(i) + alpha*p(i) + w * s(i)
c           r(i) = s(i) - w*y(i)
c 220    continue
c        beta = (dot(r,b,neq) / rr0)*(alpha/w)
c        do 225 i = 1, neq
c            p(i) = r(i) + beta*(p(i)-w*z(i))
c 225    continue
c        d = dot(r,r,neq)  
c        if (dsqrt(dabs(d)) .lt. conv) goto 300
c 230 continue
c ----------------------------------------------------------------------
c     write(*,1200) maxit
c     call stop_mef()
c 300 continue
c
c ... Energy norm:
c
c      call matvec(neq,ia,ja,ia(neq+2),ja(nad+1),ad,al,au,al(nad+1),x,z,
c     .            neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli)
c     b(1:neq) = b(1:neq)*m(1:neq)  
c     energy   = dot(x(1),b(1),neq)
c ......................................................................
c     time = MPI_Wtime()
c     time = time-time0
c ----------------------------------------------------------------------
c     if(my_id.eq.0)write(*,1100) neq,j,energy,time
c ......................................................................
c     Controle de flops
c     if(my_id.eq.0)write(10,'(a,a,i9,a,d20.10,a,d20.10,f20.2)')
c    .              "BICGSTAB: ", "it",j, " energy norm ",energy,
c    .              " tol ",tol,"time",time
c ......................................................................
c     return
c ======================================================================
c1000 format (//,5x,'SUBROTINA BICGSTAB:',/,5x,'Coeficiente da diagonal
c    . nulo ou negativo - equacao ',i7)
c1100 format(' (BICGSTAB) solver:'/
c    . 5x,'Number of equations  = ',i20/
c    . 5x,'Number of iterations = ',i20/
c    . 5x,'Energy norm          = ',d20.6/
c    . 5x,'CPU time (s)         = ',f20.2/)
c1200 format (' *** WARNING: No convergence reached after ',i4,
c    .        ' iterations !',/)
c     end
