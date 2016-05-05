c *********************************************************************
c * INIT_OPENMP: incia o openmp                                       *
c * ----------------------------------------------------------------- *
c * parametros de entrada :                                           *
c * ----------------------------------------------------------------- *
c * my_id       - identificacao do processo mpi                       *
c * nth_elmt    - numero de threads nos elementos                     *
c * num_solv    - numero de threads no solver                         *
c * omp_elmt    - openmp na fase de elementos                         *
c * omp_solv    - openmp na fase do solvertos                         *
c * ----------------------------------------------------------------- *
c * parametros de saida                                               *
c * ----------------------------------------------------------------- *
c * ----------------------------------------------------------------- *
c *********************************************************************
      subroutine init_openmp(omp_elmt,omp_solv,nth_elmt,nth_solv,my_id)
      implicit none
      include 'omp_lib.h'
      integer my_id
      integer max_num_threads_available
      integer nth_elmt,nth_solv
      logical omp_elmt,omp_solv
c ... checa se eh solicitado o uso do openmp sem o codigo ter sido
c     compilado com flag do opnemp    
      max_num_threads_available = 0                      
!$    max_num_threads_available      = omp_get_max_threads()
      if (max_num_threads_available .eq. 0 .and. 
     .   (omp_elmt .or. omp_solv) )then
        print*," Erro: Openmp disable in compile time,", 
     .         " but enable in run-time. "
        stop
      endif  
c .....................................................................
c
c ... openmp in compile time and run time
      if(my_id.eq.0 .and. max_num_threads_available .ne. 0)then
        print*,"**************  ELASTIC_MEC    ***********************"
        print*,"**************     OPENMP      ***********************"
        print*,"number of threads available: ",max_num_threads_available
        if(omp_elmt) print*,"number of threads in elmt  : ",nth_elmt
        if(omp_solv) print*,"number of threads in solv  : ",nth_solv
        print*,"******************************************************"
c ...  openmp disable in complile time     
      elseif(my_id.eq.0)then
        print*,"**************  ESLATIC_MEC    ***********************"
        print*,"**************     OPENMP      ***********************"
        print*,"Openmp disable in compile time.                      "
        print*,"******************************************************"
      endif
c .....................................................................
c
c ...      
      return
      end
* **********************************************************************
