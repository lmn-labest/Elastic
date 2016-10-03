c **********************************************************************
c *                                                                    *
c *   SOLV.F                                             31/08/2005    *
c *                                                                    *
c *   Metodos iterativos de solucao:                                   *
c *                                                                    *
c *   pcg                                                              *
c *   gmres                                                            *
c *   prediag                                                          *
c *   pbcgstab                                                         *
c *                                                                    *
c **********************************************************************
      subroutine solv(neq,nad,ip,ja,ad,au,al,m,b,x,tol,maxit,ngram,
     .               unsym,solver,neqf1i,neqf2i,neq3i,neq4i,neq_doti,
     .               i_fmapi,i_xfi,i_rcvsi,i_dspli)
      use Malloc
      implicit none
      include 'mpif.h'
      include 'parallel.fi'
      include 'openmp.fi'
      include 'time.fi'
      integer neqf1i,neqf2i
c ... ponteiros      
      integer*8 i_fmapi,i_xfi,i_rcvsi,i_dspli
      integer*8 i_z,i_r,i_g,i_h,i_y,i_c,i_s
c ......................................................................
      integer neq3i,neq4i,neq_doti
      integer ip(*),ja(*),neq,nequ,nad
      integer maxit,ngram,solver
      real*8  ad(*),au(*),al(*),m(*),x(*),b(*),tol,energy
      logical unsym
      integer neqovlp
      external matvec_csrsym1
      external dot,dot_par
      external matvec_csrc,matvec_csrcr,matvec_csrcsym,matvec_csrcrsym
      external matvec_csrc1,matvec_csrcr1
      external matvec_csrcsym1,matvec_csrcrsym1
c     OpenMP'ed subroutines
      external dot_par_omp,dot_par_omp_loopwise
      external matvec_csrc_omp,matvec_csrcsym_omp,
     .         matvec_csrcr_omp,matvec_csrcrsym_omp
c ......................................................................
c
c ... numero total de equacoes na particao overlapping:
c    (neqovlp = neq, no sequencial e no non-overlapping)
      neqovlp = neq+neq3i+neq4i
      if (omp_solv) then
         pmatrixtime = Mpi_Wtime() - pmatrixtime 
         i_threads_y = alloc_8('buffer_y',nth_solv,neq)
         call partition_matrix(ip,ja,neq,ovlp)
         pmatrixtime = Mpi_Wtime() - pmatrixtime
      endif
c ......................................................................
c
c ... Gradientes conjugados com precondicionador diagonal:
      if (solver .eq. 1) then
         if (unsym) then
            print*,'Solver 1 nao disponivel para matriz nao-simetrica !'
            call stop_mef()
         endif
         i_z = alloc_8('z       ',1,neq)
         i_r = alloc_8('r       ',1,neq)
         i_s = alloc_8('s       ',1,neq)
c ...    precondicionador diagonal:
         call aequalb(m,ad,neq) 
c ...    Comunicacao da diagonal para o caso non-overlapping:
         if (novlp) call communicate(m,neqf1i,neqf2i,i_fmapi,i_xfi,
     .                   i_rcvsi,i_dspli)
         call pre_diag(m,m,neq,.false.)
c ......................................................................
         if (ovlp) then
c ......... Overlapping:
            if (omp_solv) then
c ... com OpenMp
               call pcg_omp(neq,nad,ip,ja,ad,au,al,m,b,x,
     .                  ia(i_z),ia(i_r),ia(i_s),
     .                  tol,maxit,
c ... matvec comum:
     .                  matvec_csrcrsym_omp,dot_par_omp,
c
     .                  my_id,neqf1i,neqf2i,neq_doti,
     .                  i_fmapi,i_xfi,i_rcvsi,i_dspli,ia(i_threads_y),
     .                  .true.,.true.,.true.)
c .......................................................................
c
c ... sem OpenMp
            else
               call pcg(neq ,nad  ,ip     ,ja,
     .                 ad ,au   ,al     ,m,
     .                 b  ,x    ,ia(i_z),ia(i_r),ia(i_s),
     .                 tol,maxit,
c ... matvec comum:
c    .                 matvec_csrcrsym,dot_par
c ... matvec desenrolado:
     .                 matvec_csrcrsym1,dot_par,
c
     .                 my_id,neqf1i,neqf2i,neq_doti,i_fmapi,
     .                 i_xfi,i_rcvsi,i_dspli,
     .                 .true.,.true.,.true.,mpi)
c ......................................................................
            endif
c ......................................................................
c
c ... non-overlapping:
        else
c ... com OpenMp
            if (omp_solv) then
              call pcg_omp(neq,nad,ip,ja,ad,au,al,m,b,x,
     .                 ia(i_z),ia(i_r),ia(i_s),
     .                 tol,maxit,
c ... matvec comum:
     .                 matvec_csrcsym_omp,dot_par_omp,
c
     .                 my_id,neqf1i,neqf2i,neq_doti,
     .                 i_fmapi,i_xfi,i_rcvsi,i_dspli,ia(i_threads_y),
     .                 .true.,.true.,.true.)
            else
c .......................................................................
c
c ... sem OpenMp 
              call pcg(neq ,nad  ,ip   ,ja,
     .                ad ,au   ,al     ,m,
     .                b  ,x    ,ia(i_z),ia(i_r),ia(i_s),
     .                tol,maxit,
c ... matvec comum:
c    .                matvec_csrcrsym,dot_par
c ... matvec desenrolado:
     .                matvec_csrcsym1,dot_par,
c
     .                my_id,neqf1i,neqf2i,neq_doti,i_fmapi,
     .                i_xfi,i_rcvsi,i_dspli,
     .                .true.,.true.,.true.,mpi)
            endif
         endif
c ......................................................................
         i_s = dealloc('s       ') 
         i_r = dealloc('r       ')
         i_z = dealloc('z       ')
c ......................................................................
c
c ... Gmres com precondicionador diagonal:
      elseif(solver .eq. 2) then
         i_g = alloc_8('g       ',neqovlp,ngram+1)         
         i_h = alloc_8('h       ',ngram+1,ngram)
         i_y = alloc_8('y       ',1,ngram)
         i_c = alloc_8('c       ',1,ngram)
         i_s = alloc_8('s       ',1,ngram)
         i_r = alloc_8('r       ',1,ngram+1)
c ...... precondicionador diagonal:
         call aequalb(m,ad,neq)
c ...... Comunicacao da diagonal para o caso non-overlapping:
         if (novlp) call communicate(m,neqf1i,neqf2i,i_fmapi,i_xfi,
     .                   i_rcvsi,i_dspli)
         if(unsym) then
c ......................................................................
            if(ovlp) then
c ............ Matriz nao-simetrica, overlapping:
               if (omp_solv) then
                  call gmres_omp(neq,ip,ja,ad,au,al,m,b,x,ngram,ia(i_g),
     .                       ia(i_h),ia(i_y),ia(i_c),ia(i_s),ia(i_r),
     .                       tol,maxit,matvec_csrcr_omp,dot_par_omp,
     .                       neqovlp,my_id,neqf1i,neqf2i,neq_doti,
     .                       i_fmapi,i_xfi,i_rcvsi,i_dspli,
     .                       ia(i_threads_y))
               else
                  call gmres(neq,ip,ja,ad,au,al,m,b,x,ngram,ia(i_g),
     .                       ia(i_h),ia(i_y),ia(i_c),ia(i_s),ia(i_r),
c ... matvec comum:
c    .                        maxit,matvec_csrcr,dot_par,neqovlp)
c ... matvec desenrolado:
     .                       tol,maxit,matvec_csrcr1,dot_par,neqovlp,
     .                       my_id,neqf1i,neqf2i,neq_doti,i_fmapi,i_xfi,
     .                       i_rcvsi,i_dspli)
               endif
               call communicate(x,neqf1i,neqf2i,i_fmapi,i_xfi,
     .                          i_rcvsi,i_dspli)
c ......................................................................
            else
c ............ Matriz nao-simetrica, sequencial e non-overlapping:
               if (omp_solv) then
                  call gmres_omp(neq,ip,ja,ad,au,al,m,b,x,ngram,ia(i_g),
     .                       ia(i_h),ia(i_y),ia(i_c),ia(i_s),ia(i_r),
     .                       tol,maxit,matvec_csrc_omp,dot_par_omp,
     .                       neqovlp,my_id,neqf1i,neqf2i,neq_doti,
     .                       i_fmapi,i_xfi,i_rcvsi,i_dspli,
     .                       ia(i_threads_y))
               else
                  call gmres(neq,ip,ja,ad,au,al,m,b,x,ngram,ia(i_g),
     .                       ia(i_h),ia(i_y),ia(i_c),ia(i_s),ia(i_r),
c ... matvec comum:
c     .                       maxit,matvec_csrc,dot_par,neqovlp)
c ... matvec desenrolado:
     .                       tol,maxit,matvec_csrc1,dot_par,neqovlp,
     .                       my_id,neqf1i,neqf2i,neq_doti,i_fmapi,i_xfi,
     .                       i_rcvsi,i_dspli)
              endif
            endif
c ......................................................................
         else
c ......................................................................
            if (ovlp) then
c ............ Matriz simetrica, overlapping:
               if (omp_solv) then
                  call gmres_omp(neq,ip,ja,ad,au,al,m,b,x,ngram,ia(i_g),
     .                       ia(i_h),ia(i_y),ia(i_c),ia(i_s),ia(i_r),
     .                       tol,maxit,matvec_csrcrsym_omp,dot_par_omp,
     .                       neqovlp,my_id,neqf1i,neqf2i,neq_doti,
     .                       i_fmapi,i_xfi,i_rcvsi,i_dspli,
     .                       ia(i_threads_y))
               else
                  call gmres(neq,ip,ja,ad,au,al,m,b,x,ngram,ia(i_g),
     .                       ia(i_h),ia(i_y),ia(i_c),ia(i_s),ia(i_r),
c ... matvec comum:
c    .                        maxit,matvec_csrcrsym,dot_par,neqovlp)
c ... matvec desenrolado:
     .                       tol,maxit,matvec_csrcrsym1,dot_par,neqovlp,
     .                       my_id,neqf1i,neqf2i,neq_doti,i_fmapi,i_xfi,
     .                       i_rcvsi,i_dspli)
               endif
               call communicate(x,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,
     .                          i_dspli)
c ......................................................................
            else
c ............ Matriz simetrica, sequencial e non-overlapping:
               if (omp_solv) then
                  call gmres_omp(neq,ip,ja,ad,au,al,m,b,x,ngram,ia(i_g),
     .                       ia(i_h),ia(i_y),ia(i_c),ia(i_s),ia(i_r),
     .                       tol,maxit,matvec_csrcsym_omp,dot_par_omp,
     .                       neqovlp,my_id,neqf1i,neqf2i,neq_doti,
     .                       i_fmapi,i_xfi,i_rcvsi,i_dspli,
     .                       ia(i_threads_y))
              else
                 call gmres(neq,ip,ja,ad,au,al,m,b,x,ngram,ia(i_g),
     .                       ia(i_h),ia(i_y),ia(i_c),ia(i_s),ia(i_r),
c ... matvec comum:
c     .                       maxit,matvec_csrcsym,dot_par,neqovlp)
c ... matvec desenrolado:
     .                       tol,maxit,matvec_csrcsym1,dot_par,neqovlp,
     .                       my_id,neqf1i,neqf2i,neq_doti,i_fmapi,i_xfi,
     .                       i_rcvsi,i_dspli)
               endif
            endif
c ......................................................................
         endif
c ......................................................................
         i_r = dealloc('r       ')
         i_s = dealloc('s       ')
         i_c = dealloc('c       ')
         i_c = dealloc('y       ')
         i_h = dealloc('h       ')
         i_g = dealloc('g       ')
c ......................................................................
c
c ... Gauss:
      elseif(solver .eq. 3) then
         time0 = MPI_Wtime()
         call dtri(ad,au,al,ja,neq,unsym)
         call dsolv(ad,au,al,b,ja,neq,energy,.true.)
         time = MPI_Wtime()
         print*,'CPU time (s) ',time-time0
         x(1:neq) = b(1:neq)
c ......................................................................
c
c ... BICGSTAB com precondicionador diagonal:
      elseif(solver .eq. 4) then
c ...    precondicionador diagonal:
         call aequalb(m,ad,neq)      
c ...    Comunicacao da diagonal para o caso non-overlapping:
         if (novlp) call communicate(m,neqf1i,neqf2i,i_fmapi,i_xfi,
     .                               i_rcvsi,i_dspli)
c .....................................................................
c
c ...     
         i_c = alloc_8('tsolver ',1,neq)
         i_h = alloc_8('hsolver ',1,neq)
         i_r = alloc_8('rsolver ',1,neq)
         i_s = alloc_8('psolver ',1,neq)
         i_z = alloc_8('zsolver ',1,neq)
c .....................................................................
c
c ............ Matriz nao-simetrica 
         if(unsym) then
c ............ Matriz nao-simetrica, overlapping:
           if(ovlp) then
c ... Openmp            
             if(omp_solv) then 
c              call pbicgstab_omp_loopwise(neq,ip,ja,ad,au,al,m,b,x,
               call pbicgstab_omp(neq,ip,ja,ad,au,al,m,b,x,           
     .                            ia(i_c),ia(i_h),   
     .                            ia(i_r),ia(i_s),ia(i_z),tol,maxit,
c ... dot_par_omp_loopwise:     
c    .                            matvec_csrcr_omp,dot_par_omp_loopwise,
c ... dot_par_omp      
     .                            matvec_csrcr_omp,dot_par_omp,     
     .                            my_id,neqf1i,neqf2i,neq_doti,i_fmapi,
     .                            i_xfi,i_rcvsi,i_dspli,
     .                            ia(i_threads_y))
c ... sem openmp     
             else
               call pbicgstab(neq,ip,ja,ad,au,al,m,b,x,
     .                      ia(i_c),ia(i_h),ia(i_r),ia(i_s),ia(i_z),
     .                      tol,maxit,
c ... matvec comum:
c    .                      matvec_csrcr,dot_par,
c ... matvec desenrolado:
     .                      matvec_csrcr1,dot_par,
     .                      my_id,neqf1i,neqf2i,neq_doti,i_fmapi,
     .                      i_xfi,i_rcvsi,i_dspli)
             endif
             call communicate(x,neqf1i,neqf2i,i_fmapi,i_xfi,
     .                        i_rcvsi,i_dspli)
c .....................................................................
c
c ............ Matriz nao-simetrica, sequencial e non-overlapping:  
           else
c ... Openmp            
             if(omp_solv) then
c               call pbicgstab_omp_loopwise(neq,ip,ja,ad,au,al,m,b,x,
                call pbicgstab_omp(neq,ip,ja,ad,au,al,m,b,x,           
     .                            ia(i_c),ia(i_h),   
     .                            ia(i_r),ia(i_s),ia(i_z),tol,maxit,
c ... dot_par_omp_loopwise:     
c    .                            matvec_csrc_omp,dot_par_omp_loopwise,
c ... dot_par_omp      
     .                            matvec_csrc_omp,dot_par_omp,     
     .                            my_id,neqf1i,neqf2i,neq_doti,i_fmapi,
     .                            i_xfi,i_rcvsi,i_dspli,
     .                            ia(i_threads_y))
c .....................................................................
c     
c ... sem openmp     
c            else
               call pbicgstab(neq,nequ,nad,ip,ja,ad,au,al,m,b,x,
c              call bicgstab(neq,ip,ja,ad,au,al,m,b,x,
     .                      ia(i_c),ia(i_h),ia(i_r),ia(i_s),ia(i_z),
     .                      tol,maxit,
c ... matvec comum:
c    .                      matvec_csrc,dot_par,
c ... matvec desenrolado:
     .                      matvec_csrc1,dot_par,
     .                      my_id,neqf1i,neqf2i,neq_doti,i_fmapi,
     .                      i_xfi,i_rcvsi,i_dspli)
            endif
c .....................................................................
c
c ....................................................................     
           endif
c ....................................................................             
c
c ............ Matriz simetrica 
         else       
c ............ Matriz simetrica, overlapping:
           if(ovlp) then
c ... Openmp       
             if(omp_solv) then
c              call pbicgstab_omp_loopwise(neq,ip,ja,ad,au,al,m,b,x,
               call pbicgstab_omp(neq,ip,ja,ad,au,al,m,b,x,           
     .                         ia(i_c),ia(i_h),   
     .                         ia(i_r),ia(i_s),ia(i_z),tol,maxit,
c ... dot_par_omp_loopwise:     
c    .                         matvec_csrcrsym_omp,dot_par_omp_loopwise,
c ... dot_par_omp      
     .                         matvec_csrcrsym_omp,dot_par_omp,     
     .                         my_id,neqf1i,neqf2i,neq_doti,i_fmapi,
     .                         i_xfi,i_rcvsi,i_dspli,
     .                         ia(i_threads_y))
c .....................................................................
c     
c ... sem openmp     
             else
               call pbicgstab(neq,ip,ja,ad,au,al,m,b,x,
     .                      ia(i_c),ia(i_h),ia(i_r),ia(i_s),ia(i_z),
     .                      tol,maxit,
c ... matvec comum:
c    .                      matvec_csrcrsym,dot_par,
c ... matvec desenrolado:
     .                      matvec_csrcrsym1,dot_par,
     .                      my_id,neqf1i,neqf2i,neq_doti,i_fmapi,
     .                      i_xfi,i_rcvsi,i_dspli)
             endif
c .....................................................................
c
c ...
             call communicate(x,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,
     .                        i_dspli)
c .....................................................................
c
c ............ Matriz simetrica, sequencial e non-overlapping:
           else
c ... Openmp       
             if(omp_solv) then
c              call pbicgstab_omp_loopwise(neq,ip,ja,ad,au,al,m,b,x, 
               call pbicgstab_omp(neq,ip,ja,ad,au,al,m,b,x,           
     .                          ia(i_c),ia(i_h),   
     .                          ia(i_r),ia(i_s),ia(i_z),tol,maxit,
c ... dot_par_omp_loopwise:     
c    .                          matvec_csrcsym_omp,dot_par_omp_loopwise,
c ... dot_par_omp      
     .                          matvec_csrcsym_omp,dot_par_omp,     
     .                          my_id,neqf1i,neqf2i,neq_doti,i_fmapi,
     .                          i_xfi,i_rcvsi,i_dspli,
     .                          ia(i_threads_y))
c .....................................................................
c     
c ... sem openmp     
             else
               call pbicgstab(neq   ,nequ   ,nad    ,ip     ,ja     ,
     .                       ad     ,au     ,al     ,m      ,b      ,x, 
     .                       ia(i_c),ia(i_h),ia(i_r),ia(i_s),ia(i_z),
     .                       tol    ,maxit,
c ... matvec comum:
c    .                       matvec_csrcb  ,dot_par,
c    .                     matvec_csrcsym,dot_par,
c ... matvec desenrolado:
     .                       matvec_csrcsym1,dot_par,
     .                       my_id   ,neqf1i  ,neqf2i,
     .                       neq_doti,i_fmapi,i_xfi  ,i_rcvsi,i_dspli)
             endif
c .....................................................................
           endif  
c .....................................................................       
         endif
c .....................................................................
c
c ...             
         i_z = dealloc('zsolver ')     
         i_s = dealloc('psolver ')
         i_r = dealloc('rsolver ')
         i_h = dealloc('hsolver ')
         i_c = dealloc('tsolver ')
c ......................................................................         
      endif
c ......................................................................
      if (omp_solv) then
         pmatrixtime = Mpi_Wtime() - pmatrixtime 
         i_threads_y = dealloc('buffer_y')
         pmatrixtime = Mpi_Wtime() - pmatrixtime
      endif
      return
      end
c **********************************************************************
