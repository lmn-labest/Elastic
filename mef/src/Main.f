      program Mef
c **********************************************************************
c *                                                                    *
c *   Metodo dos elementos finitos para problemas mecanicos elasticos  *
c *                                                                    *
c **********************************************************************
      use Malloc
      implicit none
      include 'mpif.h'
      include 'omp_lib.h'
      include 'string.fi'
      include 'transiente.fi'
      include 'parallel.fi'
      include 'gravity.fi'
      include 'elementos.fi'
      include 'time.fi'
      include 'openmp.fi'
c ......................................................................
c
c ... Variaveis da estrutura interna de macro-comandos:
c
      character*8 mc,macro(40),lmacro(50)
      character*30 string
      character*80 str
      integer  iloop,imacro,nmacro,nmc,loop,j
      logical flag_macro_mesh
c ......................................................................
c
c ... Variaveis para controle de arquivos:
c
      character*80 prename,fname,name,filein
      character*80 pnodename
      integer nin,nplot,nout,nout_face
      integer logsolv,fconf,logsolvd
      integer totfiles,openflag
      integer*8 i_no,i_nfile
      integer print_nnode
      integer num_pnode
c ... arquivo de impressao nos nos ( pu,stress,stressE,stressB,flux,...)  
      integer nfiles,ifiles
      parameter ( nfiles = 5)
      logical new_file(nfiles),flag_pnd,print_flag(10)
c .....................................................................
c
c ... Variaveis de controle de solucao:
c
      integer maxit,maxnlit,tmaxnlit,ngram,stge,solver,istep
      integer ilib,ntn,code
      real*8  tol,solvtol,resid,resid0
      logical reordf,unsym
c......................................................................
c
c ... Variaveis descritivas do problema:
c
      integer nnodev,nnode,numel,numat,nen,nenv,ndf,ndm,nst
      logical fmec
c .....................................................................
c
c ... Variaveis do sistema de equacoes:
      integer neq,nad
      character*8 sia,sja,sau,sal,sad
      logical cal_neq
c .....................................................................
c
c ... precondicionador
      integer precond  
c .....................................................................
c
c ... Variaveis da interface de linha de comando
      integer nargs
      character arg*80
c .....................................................................
c
c ... Variaveis locais:
c
      integer i,k
      real*8  dot_par
c .....................................................................
c 
c ...
      logical bvtk
c ......................................................................
c
c ... Variaveis de medicao de tempo:
c
      real*8 timei      
c ----------------------------------------------------------------------
c
c ... Variaveis Colormesh:
c
      integer numcolors 
c ... Ponteiros:
      integer*8 i_colorg,i_elcolor    
c ......................................................................
c
c ... Ponteiros:
c
c ... malha
      integer*8 i_ix,i_id,i_ie,i_nload,i_eload,i_e,i_x,i_inum
      integer*8 i_ic,i_fnno
c ... arranjos locais ao elemento
      integer*8 i_xl,i_ul,i_pl,i_sl,i_ld,i_txl,i_txnl
c ... forcas e graus de liberdade 
      integer*8 i_f
      integer*8 i_u,i_tx0
      integer*8 i_tx
c ... sistema de equacoes
      integer*8 i_ia,i_ja,i_ad,i_au,i_al,i_b,i_b0,i_x0,i_bst0
c ... precondicionador
      integer*8 i_m
c ... arranjos globais (MPI - escrita)
      integer*8 i_g,i_g1,i_g2,i_g3,i_g4
c ......................................................................
c
c ... tensoes iniciais
      logical fstress0,fcstress0
c ......................................................................
c
c ... Variaveis de controle do MPI:
c
      integer status(MPI_STATUS_SIZE)      
c ......................................................................
c
c ... Macro-comandos disponiveis:
c
      data nmc /40/
      data macro/'loop    ','hextotet','mesh    ','solver  ','dt      ',
     .'pgeo    ','setprint','calneq  ','gravity ','        ','gmres   ',
     .'        ','        ','        ','        ','        ','        ',
     .'solvm   ','pmecres ','        ','        ','        ','        ',
     .'        ','        ','maxnlit ','        ','nltol   ','        ',
     .'        ','        ','setpnode','        ','        ','pnu     ',
     .'pns     ','config  ','maxit   ','solvtol ','stop    '/
c ......................................................................
c
c ... Arquivos de entrada e saida:
c
      data nin /1/, nplot /3/, logsolv /10/, nout /15/,
     .     logsolvd /16/, nout_face /17/ 
      data fconf /5/
      data flag_pnd /.false./ 
c     arquivo de impressao de nos associados aos seguintes inteiros
c     nfile = 50,51,52,...,num_pnode
c ......................................................................      
c
c ... Inicializacao MPI:
c
      call MPI_INIT( ierr )
      call MPI_COMM_SIZE( MPI_COMM_WORLD, nprcs, ierr )
      call MPI_COMM_RANK( MPI_COMM_WORLD, my_id, ierr )
c ......................................................................
c
c ... Inicializacao de variaveis da estrutura interna de macro-comandos:
c
      iloop   = 0
      imacro  = 0
      nmacro  = 0
c ......................................................................
c
c ... Inicializacao de variaveis de controle de solucao:
c
c ......................................................................
      istep   =  0
      t       =  0.d0
      dt      =  1.d0
      alfa    =  1.d0
      beta    =  1.d0
c ... tipo do problema
c ... fporomec  = problema poromecanico                    
c ... fmec      = problema mecanico              
      fmec     = .false.
c ... tensoes iniciais
c ... fstress0  = tensoes iniciais como condicao inicial
c ... fcstress0 = tensoes iniciais utilizadas       
      fstress0  = .false.
      fcstress0 = .false.
c ... reordf  =  true -> reordenacao Cuthill-McKee
      reordf    = .true.
c ... maxit   =  numero max. de iteracoes do solver iterativo
c ... solvtol =  tolerancia para o solver iterativo
c ... maxnlit =  numero max. de iteracoes nao-lineares
c ... tol     =  tolerancia do algoritmo nao-linear
c ... ngram   =  base de Krylov (utilizada somente no gmres)
c ... precond =  1 - NONE, 2 - diag, 3 - iLDLt(0), 4 - iC(0)
      maxit   =  10000
      solvtol =  1.d-11
      maxnlit =  2 
      tol     =  1.d-04
      ngram   =  50
      precond =  2
c ... unsym   =  true -> matriz de coeficientes nao-simetrica      
c ... solver  =  1 (pcg), 2 (gmres), 3 (gauss / LDU), 4 (bicgstab)
c ... stge    =  1 (csr), 2 (edges), 3 (ebe), 4 (skyline)
      unsym   = .false.
      solver  =  1
      stge    =  1
      resid0  =  0.d0
c ... ilib    =  1 define a biblioteca padrão ( default = mec )
      ilib    =  1
c ... cal_neq = calcula somente o numero de equacoes
      cal_neq = .false.
c ... escrita dos nohs quadraticos
      print_flag(1) = .false.
c ... desloc
      print_flag(2) = .true.
c ... stress
      print_flag(3) = .false.
c ... campo gravitacional (Padrao)
      gravity(1)  =   0.0d0
      gravity(2)  =   0.0d0
      gravity(3)  = -9.81d0
      gravity_mod = dsqrt(gravity(1)*gravity(1)
     .                   +gravity(2)*gravity(2)
     .                   +gravity(3)*gravity(3))
c ...
      flag_macro_mesh = .false.
c ... 
      bvtk = .false.
c ... OpenMP
      omp_elmt = .false.
      omp_solv = .false.
      nth_elmt = 1
      nth_solv = 1
c ......................................................................
c
c ... Abertura de arquivos:    
      nargs = iargc()
   10 continue
      if (my_id .eq. 0) then
c ... intervace de linha de comando        
        if(nargs .gt. 0) then
          call getarg(1,arg)
          filein = arg
        else
          print*, 'Arquivo de dados:'
          read(*,'(a)') filein
        endif 
      endif      
      if (nprcs .eq. 1) then
         open(nin, file= filein, status= 'old', err=15, action= 'read')
         goto 20
   15    continue
         print*, 'Arquivo nao existente !'
         nargs = 0
         goto 10   
      else
c ...    Passa nome do arquivo de entrada do processo 0 para os demais:
         call MPI_BCAST(filein,20,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
c ...    Nome do arquivo de dados do processo my_id:
         fname = name(filein,my_id,13)
c ...    Abre o arquivo de dados de cada processo:
         totfiles = 0
         openflag = 0
         open(nin, file= fname, status= 'old', err=11, action= 'read')
         openflag = 1
   11    continue
c ...    Testa se todos abriram sem erro:
         call MPI_ALLREDUCE(openflag,totfiles,1,MPI_INTEGER,MPI_SUM,
     .                      MPI_COMM_WORLD,ierr)
         if (totfiles .ne. nprcs) then
            if (my_id .eq. 0) print*, '*** Erro na abertura de arquivos'
            call stop_mef()
         endif
      endif
   20 continue
c ......................................................................   
      if (my_id .eq. 0) then   
c ... intervace de linha de comando        
        if(nargs .eq. 2) then
          call getarg(2,arg)
          prename = arg
        else
          print*, 'Arquivo de saida: '
          read(*,'(a)') prename
        endif 
        fname = name(prename,nprcs,15)
        open(logsolv,file=fname)
        write(logsolv,'(a)') 'Solver control flop.'
      endif
c ......................................................................      
      call MPI_barrier(MPI_COMM_WORLD,ierr)
c ......................................................................
c
c ... Leitura dos macro-comandos:
c
c ......................................................................
   50 continue
      if (iloop .eq. 0) then
         call readmacro(nin,.true.)
         write(mc,'(8a)') (word(i),i=1,8)
      else
         if (imacro .eq. 0 .or. imacro .eq. nmacro) then
            imacro = 1
         else
            imacro = imacro + 1
         endif
         mc  = lmacro(imacro)
         iloop = iloop - 1
      endif
      do 60 j = 1, nmc
         if (mc .eq. macro(j)) goto 70
   60 continue
      goto 50
   70 continue
      goto (100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,
     .     1400,1500,1600,1700,1800,1900,2000,2100,2200,2300,2400,2500,
     .     2600,2700,2800,2900,3000,3100,3200,3300,3400,3500,3600,3700,
     .     3800,3900,5000) j
c ......................................................................
c
c ... Execucao dos macro-comandos:
c
c ......................................................................
c
c ... Macro-comando LOOP:
  100 continue
      call readmacro(nin,.false.)
      write(string,'(12a)') (word(i),i=1,12)
      read(string,*,err = 120,end = 120) loop
      nmacro = 0
      imacro = 0
      iloop  = 0
  110 continue
      call readmacro(nin,.true.)
      write(mc,'(8a)') (word(i),i=1,8)        
      if (mc .eq. 'next ') goto 50
      nmacro = nmacro + 1
      iloop = loop*nmacro
      lmacro(nmacro) = mc
      goto 110
  120 continue
      print*,'Erro na leitura da macro (LOOP) !'
      goto 5000            
c ......................................................................
c
c ... Macro-comando HEXTOTET:
c
c ......................................................................
  200 continue
      if(my_id.eq.0)print*, 'Macro HEXTOTET'  
      call hexa_to_tetra(ia(i_ix),numel,nen,prename,nplot)
      goto 50
c ----------------------------------------------------------------------
c
c ... Macro-comando MESH:
c
c ......................................................................
  300 continue
      if(my_id.eq.0)print*, 'Macro MESH'  
c ... Inicializacao da estrutura de dados de gerenciamento de memoria:
c
      call init_malloc(maxmem)
c ......................................................................
c
c ...
      call init_openmp(omp_elmt,omp_solv,nth_elmt,nth_solv,my_id)
c ......................................................................
c
c
      flag_macro_mesh = .true.
c
c.... Leitura de dados:
c
c      call rdat(nnode,numel,numat,nen,ndf,ndm,nst,i_ix,i_id,i_ie,
c     .          i_nload,i_eload,i_inum,i_e,i_x,i_f,i_u,i_v,i_a,nin)
c
      call rdat(nnode     ,nnodev  ,numel      ,numat  
     .         ,nen       ,nenv   
     .         ,ndf       ,ndm     ,nst        ,i_ix  
     .         ,i_ie      ,i_inum  ,i_e        ,i_x
     .         ,i_id      ,i_nload ,i_eload    ,i_f  
     .         ,i_u       ,i_tx0   ,i_fnno
     .         ,fstress0  ,fmec    ,print_flag(1) ,nin ) 
c
c    -----------------------------------------------------------------
c    |  ix  | id | ie | nload | eload | inum | e | x | f | u | tx0 |
c
c    | fnno |    
c    -----------------------------------------------------------------
c ......................................................................      
c
c ... calculo das forcas internas devidos as tensoes inicias tensoes
c     habilitado
      if(fstress0) fcstress0 = .true.
c ......................................................................   
c
c.... Controle de tempos:
c
      soltime      = 0.d0 
      elmtime      = 0.d0
      vectime      = 0.d0
      dstime       = 0.d0
      numeqtime    = 0.d0
      reordtime    = 0.d0
      frontime     = 0.d0
      totaltime    = 0.d0
      matvectime   = 0.d0
      dottime      = 0.d0
      allgtime     = 0.d0
      gvtime       = 0.d0
      allrtime     = 0.d0
      sendtime     = 0.d0
      ovhtime      = 0.d0
      colortime    = 0.d0
      pmatrixtime  = 0.d0
      writetime    = 0.d0
      tformtime    = 0.d0
      precondtime  = 0.d0
      ifatsolvtime = 0.d0
      prebdiagtime = 0.d0
      totaltime    = MPI_Wtime()
c ......................................................................      
c
c.... Otimizacao da largura de banda:
c
      timei = MPI_Wtime()
      if (nprcs .gt. 1) then
       call reord(ia(i_ix),ia(i_inum),nno1-nno1a,nnode,numel,nen,reordf)
      else
       call reord(ia(i_ix),ia(i_inum),nnode,nnode,numel,nen,reordf)
      endif      
      reordtime = MPI_Wtime()-timei
c ......................................................................
c
c.... Numeracao nodal das equacoes:
c
      timei = MPI_Wtime()
c ... mecanico
      if(fmec) then
        call numeq(ia(i_id),ia(i_inum),ia(i_id),nnode,ndf,neq)
        if(cal_neq) print*,'Neq : ',neq
      endif
c ......................................................................
      numeqtime = MPI_Wtime()-timei
c ......................................................................
c
c ... Front Mpi
      timei = MPI_Wtime()
      call init_front(i_noLG,i_noGL,nno1,nno2,nno3,nnofi
     .                ,nno_pload,nnovG,nnoG
     .                ,nelG,nnodev,nnode,numel,ovlp,novlp,nprcs,nviz1
     .                ,nviz2,i_rreqs,i_sreqs,i_rcvs0i,i_dspl0i
     .                ,i_fmap0i,mpi)
c
c.... Mapa de equacoes de fronteira:
c
      if (ndf  .gt. 0) call frontb(ndf,ia(i_id),neq,neq1
     .                ,neq2,neq3,neq4,neq1a,neqf1,neqf2,neq32
     .                ,neq_dot
     .                ,i_fmap,i_rcvs,i_dspl,i_xf
     .                ,'fmap     ','rcvs    ','dspl    ','xf      ',0)
      frontime = MPI_Wtime()-timei
c ----------------------------------------------------------------------
c                | noLG | noGL | elLG | fmap | rcvs | dspl | 
c ----------------------------------------------------------------------
c
c ----------------------------------------------------------------------             
c                     | fmap0i | rcvs0i | dspl0i |
c ----------------------------------------------------------------------                      
c ......................................................................
c
c.... Arranjos locais de elemento:
c
c ... mecanico
      if(fmec) then
        i_xl   = alloc_8('xl      ',ndm,nenv)
        i_ul   = alloc_8('ul      ',1  ,nst)
        i_pl   = alloc_8('pl      ',1  ,nst)
c ...      
        i_txnl = alloc_8('txnl    ',  6,nen)
c ...      
        i_sl   = alloc_8('sl      ',nst,nst)
        i_ld   = alloc_4('ld      ',  1,nst)
      endif
c .....................................................................
c
c     -----------------------------------------
c     | xl | ul | pl | sl | ld |
c     -----------------------------------------
c
c ... Memoria para a estrutura de dados do sistema de equacoes:
c
      timei = MPI_Wtime()
      if (ndf .gt. 0) then
         sia = 'ia' 
         sja = 'ja' 
         sau = 'au' 
         sal = 'al' 
         sad = 'ad' 
         call datastruct(ia(i_ix),ia(i_id),ia(i_inum),nnode
     .                  ,numel   ,nen     ,ndf       ,nst
     .                  ,neq     ,stge    ,unsym     ,nad  ,nad1
     .                  ,i_ia    ,i_ja    ,i_au      ,i_al ,i_ad 
     .                  ,sia     ,sja     ,sau       ,sal  ,sad         
     .                  ,ovlp     )
      endif
      dstime = MPI_Wtime()-timei
c .....................................................................
c
c ... colorir a malha (openmp)
c
      colortime = MPI_Wtime()
      call coloredmesh(ia(i_ix),nnode,nnodev,numel,nenv,nen,numcolors
     .               ,i_colorg,i_elcolor)     
      colortime = MPI_Wtime()-colortime
c ......................................................................
c
c ......................................................................
c
c ... Memoria para o vetor de forcas e solucao:
c
c ......................................................................
      if (ndf .gt. 0) then
         i_x0  = alloc_8('x0      ',    1,neq+neq3+neq4)
         i_bst0= alloc_8('bstress0',    1,neq+neq3+neq4)
         i_b0  = alloc_8('b0      ',    1,neq+neq3+neq4)
         i_b   = alloc_8('b       ',    1,neq+neq3+neq4)
         call azero(ia(i_x0)  ,neq+neq3+neq4)
         call azero(ia(i_bst0),neq+neq3+neq4) 
         call azero(ia(i_b0)  ,neq+neq3+neq4)      
         call azero(ia(i_b )  ,neq+neq3+neq4)
c ...
         i_m   = 1
c ...  Memoria para o precondicionador diagonal:
         if(precond .eq. 2 ) then 
           i_m   = alloc_8('m       ',    1,neq)
           call azero(ia(i_m),neq)
c ......................................................................
c
c ...  Memoria para o precondicionador iLDLt e iLLT (cholesky)
         else if( precond .eq. 3 .or.  precond .eq. 4) then
           i_m   = alloc_8('m       ',    1,neq+nad)
           call azero(ia(i_m),neq+nad)
c ..................................................................... 
         endif       
      endif
c ......................................................................
      goto 50
c ----------------------------------------------------------------------
c
c ... Macro-comando : Solver
c
c ......................................................................
  400 continue
      if(my_id.eq.0)print*, 'Macro SOLVER'
      if(flag_macro_mesh) then
        print*,'Macro so pode ser utilizada antes da macro mesh'
        goto 5000
      endif
      call read_solver_config(solver,solvtol,maxit,precond,nin,my_id)
      goto 50
c ----------------------------------------------------------------------
c
c ... Macro-comando: DT
c
c ......................................................................
  500 continue
      if(my_id.eq.0)print*, 'Macro DT   '
      call readmacro(nin,.false.)
      write(string,'(30a)') (word(i),i=1,30)
      read(string,*,err = 510,end = 510) dt
      goto 50
  510 continue
      print*,'Erro na leitura da macro (DT) !'
      goto 5000      
c ----------------------------------------------------------------------
c
c ... Macro-comando: PGEO
c
c ......................................................................
  600 continue
c ...     
      ntn   = 6
c ...
      print_nnode = nnovG   
      if(print_flag(1)) print_nnode = nnoG
c ......................................................................
c
c ... Geometria:
      if(mpi) then
        
        writetime = writetime + MPI_Wtime()-timei 
        call global_ix(nen+1,numel_nov,i_ix,i_g,'ixg     ')
        call global_v(ndm   ,nno_pload,i_x ,i_g1,'xg      ')
c .....................................................................
c
c ...
        if( my_id .eq. 0 ) then
          print*, 'Macro PGEO'
          call write_mesh_geo(ia(i_g),ia(i_g1),print_nnode,nelG
     .                       ,nen    ,ndm     ,prename    ,bvtk
     .                       ,.true. ,nplot)
        endif
c .....................................................................
c
c ...        
        i_g1 = dealloc('xg      ')
        i_g  = dealloc('ixg     ')
        call MPI_barrier(MPI_COMM_WORLD,ierr)
        writetime = writetime + MPI_Wtime()-timei
c ......................................................................
c
c ...
      else
        print*, 'Macro PGEO'
        writetime = writetime + MPI_Wtime()-timei 
        call write_mesh_geo_bc(ia(i_ix) ,ia(i_x)    ,ia(i_ie)
     .                      ,ia(i_id)   ,ia(i_f)    ,ia(i_u) 
     .                      ,ia(i_tx0)  ,ia(i_nload),ia(i_eload)
     .                      ,print_nnode,numel      ,ndf     ,ntn
     .                      ,nen        ,ndm        ,prename
     .                      ,bvtk       ,macros     ,.true.
     .                      ,nplot      ,nout_face)
        writetime = writetime + MPI_Wtime()-timei
      endif
c .....................................................................
      goto 50
c ----------------------------------------------------------------------
c
c ... Macro-comando: SETPRINT
c
c ......................................................................
  700 continue
      if(my_id.eq.0) print*, 'Macro SETPRINT'
      if(flag_macro_mesh) then
        print*,'Macro so pode ser utilizada antes da macro mesh'
        goto 5000
      endif
      call set_print_vtk(print_flag,my_id,nin)      
      goto 50
c ----------------------------------------------------------------------
c
c ... Macro-comando: 
c
c ......................................................................
  800 continue
      if(my_id.eq.0)print*, 'Macro  CALNEQ'
      cal_neq = .true.
      goto 50
c ......................................................................
c
c ... Macro-comando: GRAVITY
c
c ......................................................................
  900 continue
      if (my_id .eq. 0)   print*, 'Macro Gravity'
c ... gx
      call readmacro(nin,.false.)
      write(string,'(30a)') (word(i),i=1,30)
      read(string,*,err =910,end =910) gravity(1)    
c ... gy
      call readmacro(nin,.false.)
      write(string,'(30a)') (word(i),i=1,30)
      read(string,*,err =920,end =920) gravity(2)    
c ... gz
      call readmacro(nin,.false.)
      write(string,'(30a)') (word(i),i=1,30)
      read(string,*,err =930,end =930) gravity(3)    
      gravity_mod = dsqrt(gravity(1)*gravity(1)
     .                   +gravity(2)*gravity(2)
     .                   +gravity(3)*gravity(3))
      goto 50
  910 continue
      print*,'Erro na leitura da macro (GRAVITY) gx !'
      goto 5000
  920 continue
      print*,'Erro na leitura da macro (GRAVITY) gy !'
      goto 5000
  930 continue
      print*,'Erro na leitura da macro (GRAVITY) gz !'
      goto 5000
c ----------------------------------------------------------------------
c
c ... Macro-comando:       
c
c ......................................................................
 1000 continue
      if(my_id.eq.0)print*, 'Macro      '
      goto 50
c ......................................................................
c
c ... Macro-comando:  
c
c ......................................................................
 1100 continue
      if(my_id.eq.0)print*, 'Macro      '
      goto 50
c ......................................................................
c
c ... Macro-comando: 
c
c ......................................................................
 1200 continue
      if(my_id.eq.0) print*, 'Macro '
      goto 50
c ----------------------------------------------------------------------
c
c ... Macro-comando: 
c
c ......................................................................
 1300 continue
      if(my_id.eq.0) print*, 'Macro '
      goto 50     
c ----------------------------------------------------------------------
c
c ... Macro-comando: BICGSTAB
c
c ......................................................................
 1400 continue
      if(my_id.eq.0)print*, 'Macro '
      goto 50
c ----------------------------------------------------------------------
c
c ... Macro-comando: PCG
c
c ......................................................................
 1500 continue
      goto 50 
c ----------------------------------------------------------------------
c
c ... Macro-comando:      
c
c ......................................................................
 1600 continue
      if(my_id.eq.0)print*,'Macro '
      goto 50     
c ......................................................................
c
c ... Macro-comando:         
c
c ......................................................................
 1700 continue
      goto 50
c ----------------------------------------------------------------------
c
c ... Macro-comando: SOLVM
c
c ......................................................................
 1800 continue
      if (my_id .eq. 0 ) print*, 'Macro  SOLVM'
c ...
      ilib   = 1
      i      = 1
      istep  = istep + 1
      resid0 = 0.d0
c .....................................................................
c
c ... Cargas nodais e valores prescritos no tempo t+dt:
      timei = MPI_Wtime()
      call pload_mec(ia(i_id),ia(i_f),ia(i_u),ia(i_b0),ia(i_nload)
     .              ,nnode   ,ndf)
      vectime = vectime + MPI_Wtime()-timei
c .....................................................................
c
c ... forcas internas devidos as tensoes inicias tensoes
      if(fstress0 .and. fcstress0) then
          timei = MPI_Wtime()
          call pform_mec(ia(i_ix)    ,ia(i_eload)  ,ia(i_ie) ,ia(i_e) 
     .                 ,ia(i_x)     ,ia(i_id)     ,ia(i_ia) ,ia(i_ja)
     .                 ,ia(i_au)    ,ia(i_al)     ,ia(i_ad) ,ia(i_bst0) 
     .                 ,ia(i_u )    ,ia(i_tx0)
     .                 ,ia(i_xl)    ,ia(i_ul)     ,ia(i_pl)
     .                 ,ia(i_sl)    ,ia(i_ld)     ,ia(i_txnl) 
     .                 ,numel       ,nen          ,nenv     ,ndf 
     .                 ,ndm         ,nst          ,neq      ,nad        
     .                 ,.false.     ,.true.       ,unsym 
     .                 ,stge        ,5            ,ilib     ,i
     .                 ,ia(i_colorg),ia(i_elcolor),numcolors,.true.)
          elmtime = elmtime + MPI_Wtime()-timei
          fcstress0= .false.
      endif 
c .....................................................................      
c
c ... forcas de volume e superficie do tempo t+dt :  
      timei = MPI_Wtime()
      call pform_mec(ia(i_ix)    ,ia(i_eload)  ,ia(i_ie) ,ia(i_e) 
     .              ,ia(i_x)     ,ia(i_id)     ,ia(i_ia) ,ia(i_ja)
     .              ,ia(i_au)    ,ia(i_al)     ,ia(i_ad) ,ia(i_b0) 
     .              ,ia(i_u )    ,ia(i_tx0)
     .              ,ia(i_xl)    ,ia(i_ul)     ,ia(i_pl)
     .              ,ia(i_sl)    ,ia(i_ld)     ,ia(i_txnl) 
     .              ,numel       ,nen          ,nenv     ,ndf 
     .              ,ndm         ,nst          
     .              ,neq         ,nad           
     .              ,.false.     ,.true.       ,unsym 
     .              ,stge,4      ,ilib         ,i
     .              ,ia(i_colorg),ia(i_elcolor),numcolors,.false.)
      elmtime = elmtime + MPI_Wtime()-timei
c .....................................................................
c
c ... tensao inicial
      if(fstress0) then
        call vsum(ia(i_b0),ia(i_bst0),neq,ia(i_b0))
      endif  
c .....................................................................
c
c ---------------------------------------------------------------------
c loop nao linear:
c ---------------------------------------------------------------------
 1810 continue
c ... Loop multi-corretor:      
      if(my_id.eq.0) print*,'nonlinear iteration ',i
c ...
      timei = MPI_Wtime()
      call aequalb(ia(i_b),ia(i_b0),neq)
      vectime = vectime + MPI_Wtime()-timei
c .....................................................................
c
c ... Residuo: b = F - K.u(n+1,i)
      timei = MPI_Wtime()
      call pform_mec(ia(i_ix)    ,ia(i_eload)  ,ia(i_ie) ,ia(i_e)
     .              ,ia(i_x)     ,ia(i_id)     ,ia(i_ia) ,ia(i_ja)
     .              ,ia(i_au)    ,ia(i_al)     ,ia(i_ad) ,ia(i_b)
     .              ,ia(i_u)     ,ia(i_tx0)
     .              ,ia(i_xl)    ,ia(i_ul)     ,ia(i_pl)
     .              ,ia(i_sl)    ,ia(i_ld)     ,ia(i_txnl)
     .              ,numel       ,nen          ,nenv     ,ndf
     .              ,ndm         ,nst          
     .              ,neq         ,nad        
     .              ,.true.      ,.true.       ,unsym
     .              ,stge        ,2            ,ilib     ,i
     .              ,ia(i_colorg),ia(i_elcolor),numcolors,.false.)
      elmtime = elmtime + MPI_Wtime()-timei
c .....................................................................
c
c ... Comunicacao do residuo para o caso non-overlapping:
      if (novlp) call communicate(ia(i_b),neqf1,neqf2,i_fmap,i_xf,
     .                            i_rcvs,i_dspl)
c ......................................................................
c
c ......................................................................
      resid = dsqrt(dot_par(ia(i_b),ia(i_b),neq_dot))
      if(i .eq. 1) resid0 = max(resid0,resid)
      if(my_id .eq. 0) print*,'resid/resid0',resid/resid0,'resid',resid
      if ((resid/resid0) .lt. tol) goto 1820     
c ......................................................................            
c
c ... solver (Ku(n+1,i+1) = b; u(t+dt) )
      timei = MPI_Wtime()
      call solv(neq     ,nad     ,ia(i_ia),ia(i_ja)
     .         ,ia(i_ad),ia(i_au),ia(i_al)
     .         ,ia(i_m) ,ia(i_b) ,ia(i_x0)
     .         ,solvtol ,maxit   ,ngram
     .         ,unsym   ,solver  ,neqf1   ,neqf2
     .         ,neq3    ,neq4    ,neq_dot
     .         ,i_fmap  ,i_xf    ,i_rcvs ,i_dspl)
      soltime = soltime + MPI_Wtime()-timei
c .....................................................................
c
c ... atualizacao :      u(n+1,i+1) = x(n+1,i)
      timei = MPI_Wtime()
      call update_mec(nnode,ndf,ia(i_id),ia(i_u),ia(i_x0))
      vectime = vectime + MPI_Wtime()-timei
c .....................................................................
c
c ...
      if (i .ge. maxnlit) goto 1820
      i = i + 1
      goto 1810 
c .....................................................................
c
c ---------------------------------------------------------------------
c fim do loop nao linear:
c ---------------------------------------------------------------------
c
c ...
 1820 continue 
c .....................................................................
      goto 50
c ----------------------------------------------------------------------
c
c ... Macro-comando: PMECRES
c
c ......................................................................
 1900 continue
      if(my_id.eq.0)print*,'Macro PMECRES'
c ... print_flag (true| false)
c     2 - desloc
c     3 - stress
c
c ... numero do tensor de tensoes
c ... | sxx syy szz sxy  0 0 0|
      if( ndm .eq. 2) then
        ntn = 4
c ... | sxx syy szz  sxy syz sxz |
      else if(ndm .eq. 3) then
        ntn = 6
      endif
c ......................................................................
c
c ... 
      i_tx = 1
      i_g3 = 1
      if(mpi) then
c ... comunicao
        call global_v(ndf   ,nno_pload,i_u   ,i_g ,'dispG   ')
        call global_ix(nen+1,numel_nov,i_ix  ,i_g1,'ixG     ')
        call global_v(ndm   ,nno_pload,i_x   ,i_g2,'xG      ')
        if(print_flag(3)) then
          call global_v(ntn   ,nno_pload,i_tx0 ,i_g3,'tx0G    ')
        endif
c ......................................................................
c
c
c ...
        if(my_id.eq.0) then
c ...
          if(print_flag(3)) then
            i_tx  = alloc_8('tx      ',  ntn,print_nnode)
            i_ic  = alloc_4('ic      ',    1,print_nnode)
            call azero(ia(i_tx)    ,print_nnode*ntn)
            call mzero(ia(i_ic)    ,print_nnode)
c .....................................................................
c
c ...
            timei = MPI_Wtime()
            call tform_mec(ia(i_g1) ,ia(i_g2),ia(i_e)   ,ia(i_ie)
     .                   ,ia(i_ic)  ,ia(i_xl),ia(i_ul) 
     .                   ,ia(i_txnl),ia(i_g) ,ia(i_g3),ia(i_tx) 
     .                   ,nnovG     ,nnoG
     .                   ,nelG      ,nenv      ,nen
     .                   ,ndm       ,ndf       ,nst  ,ntn
     .                   ,3         ,print_flag(1),ilib)
            tformtime = tformtime + MPI_Wtime()-timei
c ......................................................................
          endif
c ......................................................................
c
c ...
          fname = name(prename,istep,2)
          call write_mesh_res_mec(ia(i_g1) ,ia(i_g2) ,ia(i_g),ia(i_tx)
     .                         ,print_nnode,nelG
     .                         ,nen        ,ndm      ,ndf   ,ntn
     .                         ,fname      ,.false.,.true.  ,print_flag
     .                         ,nplot)
          close(nplot)  
c ......................................................................
c
c ...
          if(print_flag(3)) then
            i_ic  = dealloc('ic      ')
            i_tx  = dealloc('tx      ')
          endif
c ......................................................................
        endif
c ......................................................................
c
c ...
        if(print_flag(3)) i_g3  = dealloc('tx0G    ')
        i_g2  = dealloc('xG      ')
        i_g1  = dealloc('ixG     ')
        i_g   = dealloc('dispG   ')
c ......................................................................
c
c ...
      else
c ...
        if(print_flag(3)) then
          i_tx  = alloc_8('tx      ',  ntn,print_nnode)
          i_ic  = alloc_4('ic      ',    1,print_nnode)
          call azero(ia(i_tx)    ,print_nnode*ntn)
          call azero(ia(i_ic)    ,print_nnode)
c .....................................................................
c
c ...
         
          timei = MPI_Wtime()
          call tform_mec(ia(i_ix) ,ia(i_x)  ,ia(i_e)  ,ia(i_ie)
     .                  ,ia(i_ic)  ,ia(i_xl) ,ia(i_ul) 
     .                  ,ia(i_txnl),ia(i_u)  ,ia(i_tx0),ia(i_tx) 
     .                  ,nnodev    ,nnode  
     .                  ,numel     ,nenv      ,nen
     .                  ,ndm       ,ndf       ,nst       ,ntn
     .                  ,3         ,print_flag(1),ilib)
          tformtime = tformtime + MPI_Wtime()-timei
c ......................................................................
        endif
c ......................................................................
c
c ...
        fname = name(prename,istep,2)
        call write_mesh_res_mec(ia(i_ix) ,ia(i_x)  ,ia(i_u),ia(i_tx)
     .                       ,print_nnode,numel
     .                       ,nen        ,ndm      ,ndf   ,ntn
     .                       ,fname      ,.false.,.true.  ,print_flag
     .                       ,nplot)
        close(nplot)  
c ......................................................................
c
c ...
        if(print_flag(3)) then
          i_ic  = dealloc('ic      ')
          i_tx  = dealloc('tx      ')
        endif
c ......................................................................
      endif
c ......................................................................
      goto 50     
c ----------------------------------------------------------------------
c
c ... Macro-comando:
c
c ......................................................................
 2000 continue
      print*, 'Macro '
      goto 50
c ----------------------------------------------------------------------
c
c ... Macro-comando:
c
c ......................................................................
 2100 continue
      print*, 'Macro     '
      goto 50
c ----------------------------------------------------------------------
c
c ... Macro-comando:
c
c ......................................................................
 2200 continue
      print*, 'Macro     '
      goto 50
c ----------------------------------------------------------------------
c
c ... Macro-comando:
c
c ......................................................................
 2300 continue
      print*, 'Macro     '
      goto 50
c ----------------------------------------------------------------------
c
c ... Macro-comando:
c
c ......................................................................
 2400 continue
      print*, 'Macro     '
      goto 50
c ----------------------------------------------------------------------
c
c ... Macro-comando:
c
c ......................................................................
 2500 continue
      print*, 'Macro     '
      goto 50
c ----------------------------------------------------------------------
c
c ... Macro-comando:
c
c ......................................................................
 2600 continue
      goto 50
c ----------------------------------------------------------------------
c
c ... Macro-comando:
c
c ......................................................................
 2700 continue
      goto 50
c ----------------------------------------------------------------------
c
c ... Macro-comando: NLTOL
c
c ......................................................................
 2800 continue
      if(my_id.eq.0)print*, 'Macro NLTOL'
      call readmacro(nin,.false.)
      write(string,'(30a)') (word(i),i=1,30)
      read(string,*,err =2810,end =2810) tol
      write(*,'(a,d10.2)')' Set noliner tol for ',tol
      goto 50
 2810 continue
      print*,'Erro na leitura da macro (NLTOL) !'
      goto 5000
c ----------------------------------------------------------------------
c
c ... Macro-comando: 
c ......................................................................
 2900 continue
      print*, 'Macro     '
      goto 50
c ----------------------------------------------------------------------
c
c ... Macro-comando:
c
c ......................................................................
 3000 continue
      print*, 'Macro     '
      goto 50
c ----------------------------------------------------------------------
c
c ... Macro-comando:
c
c ......................................................................
 3100 continue
      print*, 'Macro     '
      goto 50
c ----------------------------------------------------------------------
c
c ... Macro-comando: SETPNODE impressao de grandezas por no no tempo
c
c ......................................................................
 3200 continue
      if(my_id.eq.0) print*, 'Macro SETPNODE'
      if(my_id.eq.0) then
        call readmacro(nin,.false.)
        write(str,'(80a)') (word(i),i=1,80)
        read(str,*,err=2310,end = 2310) pnodename
        goto 2320
c ... problema no arquivo auxiliar        
 2310   continue
        print*,'Erro na leitura da macro (SETPNODE)'
        flag_pnd = .false.
        goto 2330
c ... leitura normal 
 2320   continue     
        call readpnode(pnodename,i_no,i_nfile,num_pnode,flag_pnd,nout)
        new_file(1:nfiles) = .true.
 2330   continue
      endif
      call MPI_BCAST(flag_pnd,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
c ... erro na letura do nome do arquivo auxiliar      
      if( flag_pnd .eqv. .false.) call stop_mef()
      goto 50
c ----------------------------------------------------------------------
c
c ... Macro-comando:                                                    
c ......................................................................
 3300 continue
      if(my_id.eq.0) print*, 'Macro '
      goto 50
c ----------------------------------------------------------------------
c
c ... Macro-comando:                                               
c ......................................................................
 3400 continue
      if(my_id.eq.0) print*, 'Macro '
      goto 50
c ----------------------------------------------------------------------
c
c ... Macro-comando: PNDISP impressao do deslocamento por no no tempo
c     (SETPNODE)                                                   
c ......................................................................
 3500 continue
      if(my_id.eq.0) print*, 'Macro PNU'      
      if(flag_pnd.eqv..false.) then
        if(my_id.eq.0)print*,'Nemhum no de impressao para PNU!'  
        call stop_mef()
      endif
c ... codigo para o arquivo u_node.txt      
      code   = 30
      ifiles = 1
c .....................................................................
      call global_v(ndf,nno_pload,i_u,i_g1,'dispG   ')
      string = 'DeslocAndPress'
      if( my_id .eq. 0) then
        do j = 1, num_pnode
          call printnode(ia(i_g1),ia(i_no+j-1),ndf            ,istep,dt
     .                  ,string  ,prename     ,ia(i_nfile+j-1)
     .                  ,code    ,new_file(ifiles))
        enddo
        new_file(ifiles) = .false.
      endif
      if(mpi) then
        i_g1 = dealloc('dispG   ')
      endif
      call MPI_barrier(MPI_COMM_WORLD,ierr)
      goto 50
c ----------------------------------------------------------------------
c
c ... Macro-comando:                                                    
c ......................................................................
 3600 continue
      goto 50
c ----------------------------------------------------------------------
c
c ... Macro-comando: CONFIG
c
c ......................................................................
 3700 continue
      if(my_id) print*, 'Macro CONFIG    '
      if(flag_macro_mesh) then
        print*,'Macro so pode ser utilizada antes da macro mesh'
        goto 5000
      endif
      call read_config(maxmem
     .                ,omp_elmt,omp_solv
     .                ,nth_elmt,nth_solv
     .                ,reordf  ,bvtk 
     .                ,nin)
      goto 50
c ----------------------------------------------------------------------
c
c ... Macro-comando: 
c
c ......................................................................
 3800 continue
      goto 50
c ----------------------------------------------------------------------
c
c ... Macro-comando: 
c
c ......................................................................
 3900 continue
      goto 50
c ----------------------------------------------------------------------
c
c ... Macro-comando STOP:
c
c ......................................................................
 5000 continue
      if(my_id.eq.0)print*, 'Macro STOP'
c ... fecha o arquivo de entrada de dados
      close(nin)
c ... fecha o arquivo gid , view3d  e log do solv
      close(logsolv)
      if(solver .eq. 5 ) close(logsolvd)
c .....................................................................
c
c ...
      totaltime = MPI_Wtime() - totaltime
c .....................................................................
c
c ... arquivo de tempo      
      call write_log_file(nnode    ,numel   ,numel_nov,numel_ov,ndf 
     .                   ,neq      ,neq1    ,neq2     ,neq32   ,neq4  
     .                   ,neq1a    ,neqf1   ,neqf2 
     .                   ,nad      ,nad1
     .                   ,omp_elmt ,nth_elmt,omp_solv ,nth_solv
     .                   ,fmec     ,numcolors,prename
     .                   ,my_id    ,nprcs   ,nout)
c .....................................................................
c
c ... meida do tempo mpio  
      if(mpi) then    
        call mpi_log_mean_time(nnovG,nnoG,nelG
     .                        ,omp_elmt ,nth_elmt
     .                        ,omp_solv ,nth_solv
     .                        ,fmec     ,numcolors,prename
     .                        ,my_id    ,nprcs   ,nout)
      endif
c .....................................................................
c
c ... desalocando a memoria do vetor de trabalho  
      call common_finalize()
c .....................................................................
c
c ...
      call stop_mef()
c .....................................................................
      end
c **********************************************************************      
c      
      subroutine printv (v,n)
      implicit none 
      real*8 v(*)
      integer i,n
      do i = 1, n
      print*,i,v(i)
      enddo
      return
      end

