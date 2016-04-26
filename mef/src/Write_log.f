c *********************************************************************
c * WRITE_LOG : Escrever o arquivo de log de excuacao                 *
c * ----------------------------------------------------------------- *
c * parametros de entrada :                                           *
c * ----------------------------------------------------------------- *
c * nnode     - numero de nos                                         *
c * numel     - numero de elementos                                   *
c * numel_nov - numero de elementos non-overllaping                   *  
c * numel_ov  - numero de elementos overllaping                       *
c * ndf       - graus de liberdade mecanico                           *
c * ndft      - graus de liberdade termico                            *
c * neq       - numero de equacoes total                              *
c * nequ      - numero de equacoes u                                  *
c * neqp      - numero de equacoes p                                  *
c * neq1      - numero de equacoes V1 do mecanico                     *
c * neq2      - numero de equacoes V2 do mecanico                     *
c * neq32     - numero de equacoes V3 do mecanico                     *
c * neq4      - numero de equacoes V4 do mecanico                     *
c * neq1a     - numero de equacoes V1a do mecanico                    *
c * neqf1     - buffer de equacoes do mecanico                        *
c * neqf2     - buffer de equacoes do mecanico                        *
c * nad       - numero de elementos nao nulos fora da diag principal  *
c * nadu      - numero de coeficientes nao nulos do bloco u           *
c * nadu      - numero de coeficientes nao nulos do bloco p           *
c * nadup     - numero de coeficientes nao nulos do bloco up          *
c * nad1      - numero de elementos nao nulos di csrcr(overlaping)    *
c * neqt      - numero de equacoes do termico                         *
c * omp_elmt  - flag do openmp na fase de elemento                    *
c * nth_elmt  - numero de threads usado na fase de elemento           *
c * omp_solv  - flag do openmp na fase do solver                      *
c * nth_solv  - numero de threads usado do solver                     *
c * num_colors- numero de cores usado para colorir a malha            *
c * prename   - prefixo do nome do arquivo de saida                   *
c * my_id     - rank do porcesso mpi                                  *
c * nprcs     - numero de porcessos do mpi                            *
c * nlog      - numero do arquivo de saida                            *
c * ----------------------------------------------------------------- *
c * parametros de saida                                               *
c * ----------------------------------------------------------------- *
c * ----------------------------------------------------------------- *
c *********************************************************************
      subroutine write_log_file(nnode   ,numel,numel_nov ,numel_ov,ndf 
     .                         ,neq     ,neq1 ,neq2
     .                         ,neq32   ,neq4 ,neq1a,neqf1  ,neqf2 
     .                         ,nad     ,nad1
     .                         ,omp_elmt,nth_elmt,omp_solv  ,nth_solv
     .                         ,fmec    ,num_colors,prename
     .                         ,my_id ,nprcs      ,nlog)
      use Malloc
      implicit none
      include 'time.fi'
      include 'mpif.h'
c ... malha
      integer nnode,numel,numel_nov,numel_ov
c ... informacoes do sistema      
      integer neq,nequ,neqp,neq1,neq2,neq32,neq4,neq1a,neqf1,neqf2
      integer nad,nad1
c ...
      integer ndf
      logical fmec
c ... mpi      
      integer mcw,mi,mdp,ierr
      integer my_id,nprcs
c ... openmp
      integer nth_elmt,nth_solv,num_colors
      logical omp_elmt,omp_solv
c ... variaveis de arquivos      
      character*80 fname,name,prename
      integer nlog
c ...
      real*8 use_work_vector,get_buffer_size
      integer buf
c ... variaveis de controle do mpi
      mcw = MPI_COMM_WORLD
      mi  = MPI_INTEGER
      mdp = MPI_DOUBLE_PRECISION
c .....................................................................
c
c ... abre o arquivo de logs
      if(my_id .eq.0) then
        fname = name(prename,nprcs,14)
        open(nlog, file= fname)
        write(nlog,'(a)')"# Arquivo de log do poro mecanico"
        write(nlog,'(a)')"Tempos (seg):"
        write(nlog,'(a,es10.2)') "Resolucao do tempo:", mpi_wtick() 
      endif
c .....................................................................
c
c ... verifica se estamos usando o mpi     
c ... aloca a variavel onde serao armazenados os tempos      
      i_ts  = alloc_8('times   ',1,nprcs)
c .....................................................................
c
c ... Tempo levado na reord
      call MPI_GATHER(reordtime,1,mdp,ia(i_ts),1,mdp,0,mcw,ierr)
      if (my_id.eq.0) call twrite('REORD ',ia(i_ts),nprcs,nlog)
c
c
c ... Tempo levado na numeq
      call MPI_GATHER(numeqtime,1,mdp,ia(i_ts),1,mdp,0,mcw,ierr)
      if (my_id.eq.0) call twrite('NUMEQ ',ia(i_ts),nprcs,nlog)
c
c ... Tempo levado na front
c
      call MPI_GATHER(frontime,1,mdp,ia(i_ts),1,mdp,0,mcw,ierr)
      if (my_id.eq.0) call twrite('FRONT ',ia(i_ts),nprcs,nlog)
c
c ... Tempo levado na struc
c
      call MPI_GATHER(dstime,1,mdp,ia(i_ts),1,mdp,0,mcw,ierr)
      if (my_id.eq.0) call twrite('STRUC ',ia(i_ts),nprcs,nlog)
c
c ...                        
c
      call MPI_GATHER(vectime,1,mdp,ia(i_ts),1,mdp,0,mcw,ierr)
      if (my_id.eq.0) call twrite('VECTR ',ia(i_ts),nprcs,nlog)
c
c ... Tempo levado na pform
c
      call MPI_GATHER(elmtime,1,mdp,ia(i_ts),1,mdp,0,mcw,ierr)
      if (my_id.eq.0) call twrite('ELMT  ',ia(i_ts),nprcs,nlog)
c
c ... Tempo levado na tform
c
      call MPI_GATHER(tformtime,1,mdp,ia(i_ts),1,mdp,0,mcw,ierr)
      if (my_id.eq.0) call twrite('TFORM ',ia(i_ts),nprcs,nlog)
c
c ... Tempo levado no precondicionador
c
      call MPI_GATHER(precondtime,1,mdp,ia(i_ts),1,mdp,0,mcw,ierr)
      if (my_id.eq.0) call twrite('PCOND ',ia(i_ts),nprcs,nlog)
c
c ... Tempo levado no solver triagonal do solver fatorado
c
      call MPI_GATHER(ifatsolvtime,1,mdp,ia(i_ts),1,mdp,0,mcw,ierr)
      if (my_id.eq.0) call twrite('IFSOLV',ia(i_ts),nprcs,nlog)
c
c ... Tempo levado no solver
c
      call MPI_GATHER(soltime,1,mdp,ia(i_ts),1,mdp,0,mcw,ierr)
      if (my_id.eq.0) call twrite('SOLVR ',ia(i_ts),nprcs,nlog)
c
c ... Tempo levado no matvec
c
      call MPI_GATHER(matvectime,1,mdp,ia(i_ts),1,mdp,0,mcw,ierr)
      if (my_id.eq.0) call twrite('MATVC ',ia(i_ts),nprcs,nlog)
c
c ... Tempo levado no produto interno (dot)
c
      call MPI_GATHER(dottime,1,mdp,ia(i_ts),1,mdp,0,mcw,ierr)
      if (my_id.eq.0) call twrite('DOT   ',ia(i_ts),nprcs,nlog)
c
c ... Tempo levado no envia e recebimento de dados (MPI)
c
      call MPI_GATHER(sendtime,1,mdp,ia(i_ts),1,mdp,0,mcw,ierr)
      if (my_id.eq.0) call twrite('SNDRC ',ia(i_ts),nprcs,nlog)
c
c ... Tempo levado na geracao/atualizacao do buffer de dados
c     rebidos e enviados (MPI)
c
      call MPI_GATHER(ovhtime,1,mdp,ia(i_ts),1,mdp,0,mcw,ierr)
      if (my_id.eq.0) call twrite('OVERH ',ia(i_ts),nprcs,nlog)
c
c ... Tempo levado no colormesh
c
      call MPI_GATHER(colortime,1,mdp,ia(i_ts),1,mdp,0,mcw,ierr)
      if (my_id.eq.0) call twrite('COLOR ',ia(i_ts),nprcs,nlog)
c
c ... Tempo levado no matix partition
c
      if(omp_solv) then
        call MPI_GATHER(pmatrixtime,1,mdp,ia(i_ts),1,mdp,0,mcw,ierr)
        if (my_id.eq.0) call twrite('PMATRI ',ia(i_ts),nprcs,nlog)
      endif   
c
c ... Tempo levado na escrita dos resultados
c
      call MPI_GATHER(writetime,1,mdp,ia(i_ts),1,mdp,0,mcw,ierr)
      if (my_id.eq.0) call twrite('WRES   ',ia(i_ts),nprcs,nlog)
c
c ... Tempo Total
      call MPI_GATHER(totaltime,1,mdp,ia(i_ts),1,mdp,0,mcw,ierr)
      if (my_id.eq.0) call twrite('TOTAL ',ia(i_ts),nprcs,nlog)
c
c ... 
c
      if(my_id.eq.0) write(nlog,'(a)')"Memoria:"
c
c ... Total de memoria usado no vetor de trabalho
      call MPI_GATHER(use_work_vector('MB'),1,mdp,ia(i_ts),1,mdp,0,mcw
     .                               ,ierr)
      if (my_id.eq.0) call twrite('ia MB    ',ia(i_ts),nprcs,nlog)
c
c ...
c
      if(my_id.eq.0) write(nlog,'(a)')"Malha e sistema linear:"
c ... sequecianl
      if(nprcs.eq.1)then
c ... mecancio
        if(fmec)then
          write(nlog,'(a)')"Mecanico:"
          call itwrite('neq   ',neq  ,nprcs,nlog)
          call itwrite('nad   ',nad  ,nprcs,nlog)
        endif
c .....................................................................
c
c .....................................................................
        call itwrite('nnode ',nnode,nprcs,nlog)
        call itwrite('numel ',numel,nprcs,nlog)
c .... MPI       
      else
c .....................................................................
c
c ...
        if( ndf .gt. 0) then
          if(my_id.eq.0) write(nlog,'(a)')"Mecanico:"
c        
c ... numero de equacao 
          call MPI_GATHER(neq,1,mi,ia(i_ts),1,mi,0,mcw,ierr)
          if (my_id.eq.0) call itwrite('neq   ',ia(i_ts),nprcs,nlog)
c        
c ... numero de equacao no no V1
c
          call MPI_GATHER(neq1,1,mi,ia(i_ts),1,mi,0,mcw,ierr)
          if (my_id.eq.0) call itwrite('neq1  ',ia(i_ts),nprcs,nlog)
c
c ... numero de equacao no no V2
c
           call MPI_GATHER(neq2,1,mi,ia(i_ts),1,mi,0,mcw,ierr)
           if (my_id.eq.0) call itwrite('neq2  ',ia(i_ts),nprcs,nlog)
c
c ... numero de equacao no no V3          
c
           call MPI_GATHER(neq32,1,mi,ia(i_ts),1,mi,0,mcw,ierr)
           if (my_id.eq.0) call itwrite('neq3  ',ia(i_ts),nprcs,nlog)
c
c ... numero de equacao no no V4          
c
           call MPI_GATHER(neq4,1,mi,ia(i_ts),1,mi,0,mcw,ierr)
           if (my_id.eq.0) call itwrite('neq4  ',ia(i_ts),nprcs,nlog)
c
c ... numero de equacao no no V4a         
c
           call MPI_GATHER(neq1a,1,mi,ia(i_ts),1,mi,0,mcw,ierr)
           if (my_id.eq.0) call itwrite('neq1a ',ia(i_ts),nprcs,nlog)
c
c ... numero de coeficientes nao nulos    
c
           call MPI_GATHER(nad,1,mi,ia(i_ts),1,mi,0,mcw,ierr)
           if (my_id.eq.0) call itwrite('nad   ',ia(i_ts),nprcs,nlog)
c
c ... numero de coeficieno nao nulos overlapping
c
           call MPI_GATHER(nad1,1,mi,ia(i_ts),1,mi,0,mcw,ierr)
           if (my_id.eq.0) call itwrite('nad1  ',ia(i_ts),nprcs,nlog)
c
c ... numero de equacao no buffer de recebimento
c
           call MPI_GATHER(neqf1,1,mi,ia(i_ts),1,mi,0,mcw,ierr)
           if (my_id.eq.0) call itwrite('neqf1 ',ia(i_ts),nprcs,nlog)
c
c ... numero de equacao no buffer de envio
c
           call MPI_GATHER(neqf2,1,mi,ia(i_ts),1,mi,0,mcw,ierr)
           if (my_id.eq.0) call itwrite('neqf2 ',ia(i_ts),nprcs,nlog)
c
        endif
c .....................................................................
c
c ...
c
c ... numero de nos      
c
        call MPI_GATHER(nnode,1,mi,ia(i_ts),1,mi,0,mcw,ierr)
        if (my_id.eq.0) call itwrite('nnode ',ia(i_ts),nprcs,nlog)
c
c ... numero de elementos
 
        call MPI_GATHER(numel,1,mi,ia(i_ts),1,mi,0,mcw,ierr)
        if (my_id.eq.0) call itwrite('nel   ',ia(i_ts),nprcs,nlog)
 
c ... numero de elementos no-overlaping
c
        call MPI_GATHER(numel_nov,1,mi,ia(i_ts),1,mi,0,mcw,ierr)
        if (my_id.eq.0) call itwrite('nelnov',ia(i_ts),nprcs,nlog)
c
c ... numero de elementos overlaping
c
        call MPI_GATHER(numel_ov,1,mi,ia(i_ts),1,mi,0,mcw,ierr)
        if (my_id.eq.0) call itwrite('nel_ov',ia(i_ts),nprcs,nlog)
      endif
c .....................................................................
c
c ... openmp
c
      if(omp_elmt .or. omp_solv) then
        if(my_id.eq.0) write(nlog,'(a)')"Openmp:"
c
c ... numero de cores usado         
c
        call MPI_GATHER(num_colors,1,mi,ia(i_ts),1,mi,0,mcw,ierr)
        if (my_id.eq.0) call itwrite('ncolor',ia(i_ts),nprcs,nlog)
c 
c ... numero de nthreads na fase do elemento
c 
        if(omp_elmt) then
          call MPI_GATHER(nth_elmt,1,mi,ia(i_ts),1,mi,0,mcw,ierr)
          if (my_id.eq.0) call itwrite('nth_elmt',ia(i_ts),nprcs,nlog)
        endif
c ... numero de nthreads na fase do elemento
c 
        if(omp_solv)then
          call MPI_GATHER(nth_solv,1,mi,ia(i_ts),1,mi,0,mcw,ierr)
          if (my_id.eq.0) call itwrite('nth_solv',ia(i_ts),nprcs,nlog)
        endif
c
c ... poro-mecanico      
        if(ndf .gt. 0 .and. omp_solv) then
          buf = neq
          call MPI_GATHER(get_buffer_size('MB',buf),1,mdp,ia(i_ts),1,mdp
     .                   ,0,mcw,ierr)
          if (my_id.eq.0) call twrite('bfm MB',ia(i_ts),nprcs,nlog)
        endif
      endif  
c .....................................................................
c
      close(nlog) 
      return
      end
c *********************************************************************      
