      subroutine rdat(nnode   ,nnodev    ,numel  ,numat   
     .               ,nen     ,nenv
     .               ,ndf     ,ndm       ,nst    ,i_ix 
     .               ,i_ie    ,i_inum    ,i_e    ,i_x 
     .               ,i_id    ,i_nload   ,i_eload,i_f
     .               ,i_u     ,i_tx0  
     .               ,fstress0,fmec      ,nin)
c **********************************************************************
c * Data de criacao    : 10/01/2016                                    *
c * Data de modificaco : 24/04/2016                                    *
c * ------------------------------------------------------------------ *    
c * RDAT: leitura de dados do problema poromecanico.                   *
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ * 
c * nin     - arquivo de entrada                                       *
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ * 
c * nnode - numero total de nos                                        *
c * nnodev- numero de nos dos vertices                                 *
c * numel - numero de elementos                                        *
c * numat - numero de materiais                                        *
c * nen   - numero max. de nos por elemento                            *
c * nenv  - numero max. de nos geometicos por elemento                 *
c * ndf   - numero max. de graus de liberdade por no                   *
c * ndm   - dimensao (1, 2 ou 3)                                       *
c * nst   - numero de graus de liberdade por elemento                  *
c * i_ix    - ponteiro para conetividades                              *
c * i_id    - ponteiro para restricoes nodais (mecanico)               *
c * i_ie    - ponteiro para materiais                                  *
c * i_nload - ponteiro para o arranjo nload (mecanico)                 *
c * i_eload - ponteiro para o arranjo eload (mecanico)                 *
c * i_inum  - ponteiro para o arranjo inum                             *
c * i_e     - ponteiro para o arranjo e                                *
c * i_x     - ponteiro para o arranjo x                                *
c * i_f     - ponteiro para o arranjo f (mecanico)                     *
c * i_u     - ponteiro para o arranjo u (mecanico)                     *
c * i_tx0   - ponteiro para o arranjo tx(mecanico)                     *
c * fstress0- leitura de tensoes iniciais (true/false)                 *      
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ * 
c * ix    - conetividades nodais dos elementos                         *
c * id    - restricoes nodais                                          *
c * ie    - tipo de elemento                                           *
c * e     - constantes fisicas                                         *
c * x     - coordenadas nodais                                         *
c * f     - forcas e prescricoes nodais                                *
c * u     - condicoes de contorno e inicial                            * 
c * tx0   - tensoes inicias                                            *
c * nload(i,j) - numero identificador da carga na direcao i do no j    *
c * load(1,n)  - tipo da carga n                                       *
c * load(2,n)  - numero de termos da carga n                           *
c * fload(i,j,k) - coeficiente i do termo j da carga k                 *
c **********************************************************************
      use Malloc
      implicit none
c ......................................................................      
      include 'elementos.fi'
      include 'load.fi'
      include 'string.fi'
      include 'parallel.fi'
      include 'termprop.fi'
c ......................................................................      
      integer nnodev,nnode,numel,numat,nen,nenv,ndf,ndm,nst,ntn
      integer maxgrade
c ... ponteiros      
      integer*8 i_e,i_x,i_f,i_nload,i_eload,i_inum
      integer*8 i_u,i_tx0
      integer*8 i_ix,i_id,i_ie
      integer*8 i_nelcon,i_nodcon,i_nincid,i_incid
c ......................................................................      
      integer nin
      integer j,nmc,totnel,nmacros,itmp
      character*15 macro(30),rc
      character*80 fname
      integer naux
      integer nincl /7/
      logical fstress0,fmec
      logical f_read_el /.false./
      logical el_quad  /.false./
c ......................................................................
      data macro/'materials      ','bar2           ','tria3          ',
     .           'quad4          ','tetra4         ','hexa8          ',
     .           'hexa20         ','               ','               ',
     .           'coordinates    ','               ','constraindisp  ',
     .           'nodalforces    ','elmtloads      ','nodalloads     ',
     .           '               ','               ','               ',
     .           'loads          ','               ','               ',
     .           'initialdisp    ','               ','initialstress  ',
     .           'parallel       ','insert         ','return         ',
     .           '               ','               ','end            '/
      data nmc /30/      
c ......................................................................
c
c ... Leitura dos parametros da malha: nnode,numel,numat,nen,ndf,ndm
      print*,'load parameters ...'
      call parameters(nnodev,numel,numat,nen,ndf,ndm,nin)
      print*,'load.'
      nnode  = nnodev
      nst    = 0
c ......................................................................
c
c ... tipo do problema
      fmec     = .false.
c ... mecancico 
      if( ndf .eq. 2 .or. ndf .eq. 3 ) then
        fmec = .true.
      endif
c ......................................................................
c
c ... numero do tensor de tensoes
c ... | sxx syy szz sxy |
      if( ndm .eq. 2) then
        ntn = 4
c ... | sxx syy szz  sxy syz sxz |
      else if(ndm .eq. 3) then
        ntn = 6
      endif
c ......................................................................
c   
c     Alocacao de arranjos na memoria: mecancio
c     ---------------------------------------------------------------
c     | ix | ie | e | x | eload |
c     ---------------------------------------------------------------
c
      if(fmec) then
        i_ix    = alloc_4('ix      ',nen+1,numel)
        i_ie    = alloc_4('ie      ',    1,numat)
        i_e     = alloc_8('e       ', prop,numat)
        i_x     = alloc_8('x       ',  ndm,nnodev)
        i_eload = alloc_4('eload   ',    7,numel)
        call mzero(ia(i_ix),numel*(nen+1))
        call mzero(ia(i_ie),numat) 
        call azero(ia(i_e),numat*10)
        call azero(ia(i_x),nnodev*ndm)
        call mzero(ia(i_eload),numel*7)
      endif 
c ......................................................................
c
c ......................................................................
      totnel  = 0            
      nmacros = 0
c ......................................................................                  
  100 continue
      call readmacro(nin,.true.)
      write(rc,'(15a)') (word(j),j=1,15)
      do 200 j = 1, nmc
         if (rc .eq. macro(j)) go to 300
  200 continue
      go to 100
  300 continue
c ......................................................................      
c
c ... Controle de Macros (PRE-processador):
c
      nmacros = nmacros + 1
      write(macros(nmacros),'(15a)') rc
c ......................................................................      
      go to (400, 450, 500,    ! materials ,bar2         ,tria3
     .       550, 600, 650,    !quad4      ,tetra4       ,hexa8
     .       700, 750, 800,    !hexa20     ,             ,
     .       850, 900, 950,    !coordinates,             ,constraindisp
     .      1000,1050,1100,    !nodalforces,elmtloads    ,nodalloads
     .      1150,1200,1250,    !           ,             ,
     .      1300,1350,1400,    !loads      ,             ,
     .      1450,1500,1550,    !initialdisp,             ,initialstress
     .      1600,1650,1700,    !parallel   ,insert       ,return
     .      1900,1900,2000) j  !           ,             ,end
c ......................................................................
c
c ... Propriedades dos materiais:
c
  400 continue
      call mate(ia(i_ie),ia(i_e),numat,nin)
      go to 100
c ......................................................................      
c
c ... Conetividades bar2:    
c
  450 continue
      nbar2(1) = 0
      call elconn(ia(i_ix),nen+1,2,nbar2(1),numel,nin)
      nbar2(2) = totnel+1
      totnel   = totnel + nbar2(1)      
      go to 100
c ......................................................................      
c
c ... Conetividades tria3:          
c
  500 continue
      print*,'load tria3 ...'
      f_read_el = .true.
      ntria3(1) = 0
      nenv      = 3
      nst       = nenv*ndf
      call elconn(ia(i_ix),nen+1,3,ntria3(1),numel,nin)
      ntria3(2) = totnel+1
      totnel    = totnel + ntria3(1)
c ... transforma os elementos lineares em quadraticos (3 nos)
c     if( nen .eq. 3) then
c     endif
c ...                                                                   
c     Alocacao de arranjos na memoria:
c     ---------------------------------------------------------------
c     | id  nload | inum | e | x | f | u | tx0 |
c     ---------------------------------------------------------------
      if (ndf .gt. 0) then
        i_inum  = alloc_4('inum    ',    1,nnode)  
        i_id    = alloc_4('id      ',  ndf,nnode)
        i_nload = alloc_4('nload   ',  ndf,nnode)
        i_f     = alloc_8('f       ',  ndf,nnode)
        i_u     = alloc_8('u       ',  ndf,nnode)
        i_tx0   = alloc_8('tx0     ',  ntn,nnode)
        call mzero(ia(i_inum) ,nnode)  
        call mzero(ia(i_id)   ,nnode*ndf)
        call mzero(ia(i_nload),nnode*ndf)
        call azero(ia(i_f)    ,nnode*ndf)
        call azero(ia(i_u)    ,nnode*ndf)
        call azero(ia(i_tx0)  ,nnode*ntn)
      endif
c .....................................................................
      print*,'load.'
      go to 100
c ......................................................................      
c
c ... Conetividades quad4:          
c
  550 continue
      print*,'load quad4 ...'
      f_read_el = .true.
      nquad4(1) = 0
      nenv      = 4
      nst       = nenv*ndf
      call elconn(ia(i_ix),nen+1,4,nquad4(1),numel,nin)
      nquad4(2) = totnel+1
      totnel    = totnel + nquad4(1)
c ... transforma os elementos lineares em quadraticos (8 nos)
c     if( nen .eq. 8) then
c     endif
c ...                                                                   
c     Alocacao de arranjos na memoria:
c     ---------------------------------------------------------------
c     | id  nload | inum | e | x | f | u | tx0 |
c     ---------------------------------------------------------------
      if (ndf .gt. 0) then
        i_inum  = alloc_4('inum    ',    1,nnode)  
        i_id    = alloc_4('id      ',  ndf,nnode)
        i_nload = alloc_4('nload   ',  ndf,nnode)
        i_f     = alloc_8('f       ',  ndf,nnode)
        i_u     = alloc_8('u       ',  ndf,nnode)
        i_tx0   = alloc_8('tx0     ',  ntn,nnode)
        call mzero(ia(i_inum) ,nnode)  
        call mzero(ia(i_id)   ,nnode*ndf)
        call mzero(ia(i_nload),nnode*ndf)
        call azero(ia(i_f)    ,nnode*ndf)
        call azero(ia(i_u)    ,nnode*ndf)
        call azero(ia(i_tx0)  ,nnode*ntn)
      endif
c .....................................................................
      print*,'load.'
      go to 100 
c ......................................................................      
c
c ... Conetividades tetra4:          
c
  600 continue
      print*,'load tetra4 ...'
      f_read_el  = .true.
      ntetra4(1) = 0
      nenv       = 4
      nst       = nenv*ndf  
      call elconn(ia(i_ix),nen+1,nenv,ntetra4(1),numel,nin)
      ntetra4(2) = totnel+1
      totnel    = totnel + ntetra4(1)
c ... transforma os elementos lineares em quadraticos (10 nos)
      if( nen .eq. 10) then
        if(fmec) then
          itmp =  nen*ndf 
        endif
        nst   = max(nst,itmp)
        el_quad   = .true.
c .....................................................................
c
c ... 
        i_nincid = alloc_4('nincid  ',1,nnodev) 
        call nodegrade(ia(i_ix),nnodev,numel,nenv,nen,ia(i_nincid)
     .                ,maxgrade) 
        i_incid  = alloc_4('incid   ',maxgrade,nnode)
        call elmincid(ia(i_ix),ia(i_incid),ia(i_nincid),nnodev,numel
     .               ,nenv    ,nen        ,maxgrade)
c ... gera a conectividade dos elementos quadraticos
        call mk_elconn_quad_v1(ia(i_ix),ia(i_incid),ia(i_nincid)
     .                        ,numel     ,nnode      ,nnodev
     .                        ,nen       ,nenv       ,maxgrade)
c .....................................................................
c
c ...
        i_incid     = dealloc('incid   ')
        i_nincid    = dealloc('nincid  ')
c .....................................................................
      endif
c ...                                                                   
c     Alocacao de arranjos na memoria:
c     ---------------------------------------------------------------
c     | id  nload | inum | e | x | f | u | tx0 |
c     ---------------------------------------------------------------
      if (ndf .gt. 0) then
        i_inum  = alloc_4('inum    ',    1,nnode)  
        i_id    = alloc_4('id      ',  ndf,nnode)
        i_nload = alloc_4('nload   ',  ndf,nnode)
        i_f     = alloc_8('f       ',  ndf,nnode)
        i_u     = alloc_8('u       ',  ndf,nnode)
        i_tx0   = alloc_8('tx0     ',  ntn,nnode)
        call mzero(ia(i_inum) ,nnode)  
        call mzero(ia(i_id)   ,nnode*ndf)
        call mzero(ia(i_nload),nnode*ndf)
        call azero(ia(i_f)    ,nnode*ndf)
        call azero(ia(i_u)    ,nnode*ndf)
        call azero(ia(i_tx0)  ,nnode*ntn)
      endif
c .....................................................................
      print*,'load.'
      go to 100
c ......................................................................      
c
c ... Conetividades hexa8:          
c
  650 continue
      print*,'load hexa8 ...'
      f_read_el = .true.
      nhexa8(1) = 0
      nenv      = 8
      nst       = nenv*ndf  
      call elconn(ia(i_ix),nen+1,nenv,nhexa8(1),numel,nin)
      nhexa8(2) = totnel + 1
      totnel    = totnel + nhexa8(1)
c ... transforma os elementos lineares em quadraticos (20 nos)
      if(nen .eq. 20) then
        if(fmec) then
          itmp =  nen*ndf 
        endif
        nst       = max(nst,itmp) 
        el_quad   = .true.
c .....................................................................
c
c ...  
        i_nincid = alloc_4('nincid  ',1,nnodev) 
        call nodegrade(ia(i_ix),nnodev,numel,nenv,nen,ia(i_nincid)
     .                ,maxgrade) 
        i_incid  = alloc_4('incid   ',maxgrade,nnode)
        call elmincid(ia(i_ix),ia(i_incid),ia(i_nincid),nnodev,numel
     .               ,nenv    ,nen        ,maxgrade)
c ... gera a conectividade dos elementos quadraticos
        call mk_elconn_quad_v1(ia(i_ix),ia(i_incid),ia(i_nincid)
     .                        ,numel   ,nnode      ,nnodev
     .                        ,nen     ,nenv       ,maxgrade)
c .....................................................................
c
c ...
        i_incid     = dealloc('incid   ')
        i_nincid    = dealloc('nincid  ')
c .....................................................................
      endif
c ...                                                                   
c     Alocacao de arranjos na memoria:
c     ---------------------------------------------------------------
c     | id  nload | inum | e | x | f | u | tx0 |
c     ---------------------------------------------------------------
      if (ndf .gt. 0) then
        i_inum  = alloc_4('inum    ',    1,nnode)  
        i_id    = alloc_4('id      ',  ndf,nnode)
        i_nload = alloc_4('nload   ',  ndf,nnode)
        i_f     = alloc_8('f       ',  ndf,nnode)
        i_u     = alloc_8('u       ',  ndf,nnode)
        i_tx0   = alloc_8('tx0     ',  ntn,nnode)
        call mzero(ia(i_inum) ,nnode)  
        call mzero(ia(i_id)   ,nnode*ndf)
        call mzero(ia(i_nload),nnode*ndf)
        call azero(ia(i_f)    ,nnode*ndf)
        call azero(ia(i_u)    ,nnode*ndf)
        call azero(ia(i_tx0)  ,nnode*ntn)
      endif
c .....................................................................
      print*,'load.'
      go to 100
c ......................................................................      
c
c ... Conetividades hexa20:          
c
  700 continue
      nhexa8(1) = 0
      call elconn(ia(i_ix),nen+1,20,nhexa8(1),numel,nin)
      nhexa8(2) = totnel+1
      totnel    = totnel + nhexa8(1)
      go to 100
c ......................................................................      
c
c ...                                                          
c
  750 continue
      go to 100
c ......................................................................      
c
c ...                                                                         
c      
  800 continue
      go to 100                                      
c ......................................................................    
c
c ... Coordenadas:
c
  850 continue
      print*,'load coordinates ...'
      call coord(ia(i_x),nnodev,ndm,nin)
      print*,'load.'
      go to 100
c ......................................................................      
c
c ... 
c
  900 continue
      go to 100
c ......................................................................      
c
c ... constraindisp - restricoes nodais (deslocamentos)
c
  950 continue
      print*,'load constraindisp ...'
      if(f_read_el) then
        call bound(ia(i_id),nnodev,ndf,nin,1)
c ...
        if(el_quad) then
          call mk_bound_quad(ia(i_id),ia(i_ix),numel,ndf,ndf,nen)
        endif
c ......................................................................
      else
        print*,'MACRO: constraindisp !! elementos nao lidos'
      endif
      print*,'load.'
      go to 100
c ......................................................................      
c
c ... nodalforces - forcas nodais:
c
 1000 continue
      if(f_read_el) then
        print*,'load nodalforces ...'
        call forces(ia(i_f),nnodev,ndf,nin)
c ...
        if(el_quad) then
          call mk_forces_quad(ia(i_id),ia(i_f),ia(i_ix),numel,ndf,nen)
        endif
c ......................................................................
      else
        print*,'MACRO: nodalforces !! elementos nao lidos'
      endif
      print*,'load.'
      go to 100
c ......................................................................      
c
c ... elmtloads - cargas nos elementos                                                    
c
 1050 continue
      print*,'load elmtloads ...'
      if(f_read_el) then
        call bound(ia(i_eload),numel,7,nin,3) 
      else
        print*,'MACRO: elmtloads !! elementos nao lidos'
      endif
      print*,'load.'
      go to 100
c ......................................................................      
c
c ... nodalloads - nos com cargas variaveis no tempo:
c
 1100 continue
      if(f_read_el) then
        call bound(ia(i_nload),nnodev,ndf,nin,2) 
      else
        print*,'MACRO: nodalloads !! elementos nao lidos'
      endif
      goto 100
c ......................................................................      
c
c ...                                                               
c
 1150 continue
      goto 100
c ......................................................................
c
c ...                                                                           
c      
 1200 continue
      go to 100
c ......................................................................
c
c ...                                                                                      
c      
 1250 continue
      go to 100
c ......................................................................
c
c ... Definicao das cargas variaveis no tempo:
c
 1300 continue
      print*,'load loads ...'
      call rload(nin)
      print*,'load.'
      goto 100
c ......................................................................
c
c ...                                                                     
c
 1350 continue
      goto 100
c ......................................................................
c
c ...                                     
c
 1400 continue
      goto 100
c ......................................................................      
c
c ... initialdisp - deslocamentos iniciais:
c
 1450 continue
c     if(fReadEl) then
c       call init_poro_mec(ia(i_u0),nnodev,ndf,1,ndf-1,nin)  
c     else
c       print*,'MACRO: initialdisp !! elementos nao lidos'
c     endif
      go to 100
c ......................................................................      
c
c ...                           
c
 1500 continue
      go to 100                                                       
c ......................................................................
c
c ... initialstress - tensao inicial                                                        
c      
 1550 continue
      print*,'load initialstress ...'
      if(f_read_el) then
        fstress0 = .true.  
        call init_poro_mec(ia(i_tx0),nnodev,ntn,1,ntn,nin)
        if(el_quad) then
          call mk_initial_quad(ia(i_tx0),ia(i_ix),numel,ntn,nen)
        endif
      else
        print*,'MACRO: initialstress !! elementos nao lidos'
      endif
      print*,'load.'
      go to 100
c ......................................................................
c
c ... Paralelo:                                         
c      
 1600 continue
      call read_par(nin,nnode,numel) 
      go to 100 
c ......................................................................
c
c ... (insert) Desvia leitura para arquivo auxiliar:
c      
 1650 continue
      nmacros = nmacros - 1
      naux = nin
      call readmacro(nin,.false.)
      write(fname,'(80a)') (word(j),j=1,strl)
      open(nincl, file= fname,status= 'old',err=1651,action='read')
      nin = nincl
      go to 100
 1651 continue
      print*, trim(fname), ' arquivo nao existente !'
      stop

c ......................................................................
c
c ... (return) Retorna leitura para arquivo de dados basico:
c      
 1700 continue
      nmacros = nmacros - 1
      close(nincl)
      nin = naux
      go to 100
c ......................................................................
c
c ...                                       
c      
 1900 continue
      go to 100                  
c ......................................................................      
c
c ... End:
c
 2000 continue
c
c ... Inclui macro de paralelo (PRE-processador):
c      
      if (nprcs .gt. 1) then
         write(macros(nmacros),'(15a)') 'parallel'
         nmacros = nmacros + 1
         write(macros(nmacros),'(15a)') 'end'
      endif
c
c ... Inicializa as condicoes de contorno no vetor u:
c
      if(ndf .gt.0) call boundc(nnode,ndf,ia(i_id),ia(i_f),ia(i_u))
      return
c ......................................................................      
      end
      subroutine elconn(ix,nen1,nen,nel,numel,nin)
c **********************************************************************
c *                                                                    *
c *   Conetividades nodais.                                            *
c *                                                                    *
c **********************************************************************
      implicit none
      include 'string.fi'      
      integer ix(nen1,*),nen1,nen,nel,numel,nin,j,k,m
      character*12 string
c ......................................................................
       do 100 j = 1, numel
          read(nin,*) k,(ix(m,k),m=1,nen),ix(nen1,k)
  100  continue
       call readmacro(nin,.true.)
       write(string,'(12a)') (word(j),j=1,12)  
       if (string .ne. 'end') goto 200
       nel=numel
       return
c ......................................................................
c      nel = 0
c      call readmacro(nin,.true.)
c      write(string,'(12a)') (word(j),j=1,12)
c      do 100 while(string .ne. 'end')
c        read(string,*,err = 200,end = 200) k
c        if(k .lt. 1 .or. k .gt. numel) goto 200      
c        do 10 j = 1, nen
c           call readmacro(nin,.false.)
c           write(string,'(12a)') (word(m),m=1,12)
c           read(string,*,err = 200,end = 200) ix(j,k)
c  10    continue
c        call readmacro(nin,.false.)
c        write(string,'(12a)') (word(m),m=1,12)
c        read(string,*,err = 200,end = 200) ix(nen1,k)   
c        nel = nel + 1
c        call readmacro(nin,.true.)
c        write(string,'(12a)') (word(j),j=1,12)
c  100 continue
c      return
c ......................................................................
  200 continue
      print*,'*** Erro na leitura dos elementos !'
      stop                   
      end
      subroutine coord(x,nnode,ndm,nin)
c **********************************************************************
c *                                                                    *
c *   Coordenadas.                                                     *
c *                                                                    *
c **********************************************************************
      implicit none
      include 'string.fi'
      integer nnode,ndm,nin,j,k,m
      real*8 x(ndm,nnode)
      character*30 string
c ......................................................................
      do 100 j = 1, nnode
        read(nin,*) k,(x(m,k),m=1,ndm)
  100 continue
      call readmacro(nin,.true.)
      write(string,'(12a)') (word(j),j=1,12)
      if (string .ne. 'end') goto 200
      return
c ......................................................................
c      call readmacro(nin,.true.)
c      write(string,'(12a)') (word(j),j=1,12)
cc     do 100 while(string .ne. 'end')
c         read(string,*,err = 200,end = 200) k
cc        if(k .lt. 1 .or. k .gt. nnode) goto 200
c         do 10 j = 1, ndm
c            call readmacro(nin,.false.)
c            write(string,'(30a)') (word(m),m=1,30)
c            read(string,*,err = 200,end = 200) x(j,k)         
c  10     continue
c         call readmacro(nin,.true.)
c         write(string,'(12a)') (word(j),j=1,12)   
c 100  continue
c      return
c ......................................................................                        
  200 continue
      print*,'*** Erro na leitura das coordenadas !'
      stop             
      end
      subroutine mate(ie,e,numat,nin)
c **********************************************************************
c *                                                                    *
c *   Materiais.                                                       *
c *                                                                    *
c **********************************************************************
      implicit none
      include 'string.fi'
      include 'termprop.fi'
      integer ie(*),numat,nin,i,j,m,ma
      real*8  e(prop,*)
      character*30 string
c ......................................................................
      call readmacro(nin,.true.)
      write(string,'(12a)') (word(j),j=1,12)
c .....................................................         
      do 110 while(string .ne. 'end')
         read(string,*,err = 200,end = 200) ma
         if(ma .lt. 1 .or. ma .gt. numat) goto 200
c .....................................................            
         call readmacro(nin,.false.)
         write(string,'(12a)') (word(m),m=1,12)
         read(string,*,err = 200,end = 200) ie(ma)
c .....................................................     
         call readmacro(nin,.false.)
         write(string,'(30a)') (word(m),m=1,30)
         i = 0
c
c ... linha de comando original:
c         do 100 while(string .ne. ' ')
c
c ... formato necessário ao funcionamento em linux:
         do 100 while ( (string .ne.   ' '    ) .and.
     .                  (string .ne. CHAR(13) ) )
            i = i + 1
            if(i .gt. prop) goto 200
            read(string,*,err = 200,end = 200) e(i,ma)
            call readmacro(nin,.false.)
            write(string,'(30a)') (word(m),m=1,30)               
  100    continue
         call readmacro(nin,.true.)
         write(string,'(12a)') (word(j),j=1,12)
  110 continue
      return
c ......................................................................                        
  200 continue
      print*,'*** Erro na leitura dos materiais !'
      stop             
      end      
      subroutine forces(f,nnode,ndf,nin)
c **********************************************************************
c *                                                                    *
c *   Forcas nodais e valores prescritos.                              *
c *                                                                    *
c **********************************************************************
      implicit none
      include 'string.fi'      
      integer nnode,ndf,nin,i,j,k
      real*8 f(ndf,*)
      character*30 string            
c ......................................................................
      call readmacro(nin,.true.)
      write(string,'(12a)') (word(i),i=1,12)
      do 100 while(string .ne. 'end')
         read(string,*,err = 200,end = 200) k
         if(k .lt. 1 .or. k .gt. nnode) goto 200      
         do 10 j = 1, ndf
            call readmacro(nin,.false.)
            write(string,'(30a)') (word(i),i=1,30)
            read(string,*,err = 200,end = 200) f(j,k)
   10    continue
         call readmacro(nin,.true.)
         write(string,'(12a)') (word(i),i=1,12)
  100 continue
      return
c ......................................................................      
  200 continue
      print*,'*** Erro na leitura das forcas nodais !'
      stop             
      end
      subroutine init_poro_mec(f,nnode,ndf,n1,n2,nin)
c **********************************************************************
c *                                                                    *
c *   Valores iniciais para o problema poro mecanico                   *
c *                                                                    *
c **********************************************************************
      implicit none
      include 'string.fi'      
      integer n1,n2,nnode,ndf,nin,i,j,k
      real*8 f(ndf,*)
      character*30 string            
c ......................................................................
      call readmacro(nin,.true.)
      write(string,'(12a)') (word(i),i=1,12)
      do 100 while(string .ne. 'end')
         read(string,*,err = 200,end = 200) k
         if(k .lt. 1 .or. k .gt. nnode) goto 200      
         do 10 j = n1, n2
            call readmacro(nin,.false.)
            write(string,'(30a)') (word(i),i=1,30)
            read(string,*,err = 200,end = 200) f(j,k)
   10    continue
         call readmacro(nin,.true.)
         write(string,'(12a)') (word(i),i=1,12)
  100 continue
      return
c ......................................................................      
  200 continue
      print*,'*** Erro na leitura das forcas nodais !'
      stop             
      end
      subroutine bound(id,nnode,ndf,nin,code)
c **********************************************************************
c *                                                                    *
c *   Coordenadas.                                                     *
c *                                                                    *
c **********************************************************************
      implicit none
      include 'string.fi'
      integer id(ndf,*),nnode,ndf,nin,code,i,j,k
      character*30 string
c ......................................................................
      call readmacro(nin,.true.)
      write(string,'(12a)') (word(j),j=1,12)
c .....................................................         
      do 110 while(string .ne. 'end')
         read(string,*,err = 200,end = 200) j
         if(j .lt. 1 .or. j .gt. nnode) goto 200
c .....................................................            
         call readmacro(nin,.false.)
         write(string,'(30a)') (word(k),k=1,30)
         i = 0
c
c ... linha de comando original:
c         do 100 while(string .ne. ' ')
c
c ... formato necessário ao funcionamento em linux:
         do 100 while ( (string .ne.   ' '    ) .and.
     .                  (string .ne. CHAR(13) ) )
            i = i + 1
            if(i .gt. ndf) goto 200
            read(string,*,err = 200,end = 200) id(i,j)
            call readmacro(nin,.false.)
            write(string,'(30a)') (word(k),k=1,30)               
  100    continue
         call readmacro(nin,.true.)
         write(string,'(12a)') (word(j),j=1,12)
  110 continue
      return
c ......................................................................                        
  200 continue
      if    (code .eq. 1) then
         print*,'*** Erro na leitura das restricoes !'
      elseif(code .eq. 2) then
         print*,'*** Erro na leitura das cargas nodais !'
      elseif(code .eq. 3) then
         print*,'*** Erro na leitura das cargas nos elementos !'
      endif
      stop             
      end            
      subroutine bound0(id,nnode,ndf,nin)
c **********************************************************************
c *                                                                    *
c *   Restricoes nodais.                                               *
c *                                                                    *
c **********************************************************************
      implicit none
      include 'string.fi'
      integer id(ndf,*),nnode,ndf,nin,i,j,k
      character*12 string      
c ......................................................................
      call readmacro(nin,.true.)
      write(string,'(12a)') (word(i),i=1,12)
      do 100 while(string .ne. 'end')
         read(string,*,err = 200,end = 200) k
         if(k .lt. 1 .or. k .gt. nnode) goto 200      
         do 10 j = 1, ndf
            call readmacro(nin,.false.)
            write(string,'(12a)') (word(i),i=1,12)
            read(string,*,err = 200,end = 200) id(j,k)
   10    continue
         call readmacro(nin,.true.)
         write(string,'(12a)') (word(i),i=1,12)
  100 continue
      return
c ......................................................................      
  200 continue
      print*,'*** Erro na leitura das restricoes nodais !'
      stop             
      end
      subroutine rload(nin)
c **********************************************************************
c *                                                                    *
c *   Leitura das cargas variaveis no tempo                            *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *                                                                    *
c *    idl(6) - variavel auxuliar                                      *
c *    nin - arquivo de entrada                                        *
c *                                                                    *
c *   Parametros de saida:                                             *
c *                                                                    *
c *    load(1,n) - tipo da carga n (n <= numload)                      *
c *    load(2,n) - numero de termos da carga n                         *
c *    fload(i,j,k) - coeficiente do termo j da carga k                *
c *                                                                    *
c **********************************************************************
      implicit none
      include 'string.fi'
      include 'load.fi'
      integer  nin,i,j,k,nc
      character*30 string            
c ......................................................................
      call readmacro(nin,.true.)
      write(string,'(12a)') (word(i),i=1,12)
      do 100 while(string .ne. 'end')
c ...    numero da carga:
         read(string,*,err = 200,end = 200) i
         if(i .lt. 1 .or. i .gt. numload) goto 200 
         call readmacro(nin,.false.)
         write(string,'(12a)') (word(j),j=1,12)
c ...    tipo da carga:
         read(string,*,err = 200,end = 200) load(1,i)
         load(2,i) = 1
c .....................................................................
c   
c ...   
         if     (load(1,i) .eq. 1) then
c ...       valor da carga:            
            call readmacro(nin,.false.)
            write(string,'(30a)') (word(j),j=1,30)
            read(string,*,err = 200,end = 200) fload(1,1,i)
c .....................................................................
c      
c ...       u(t) = a.t + b 
         elseif (load(1,i) .eq. 2) then
c ...       valor de b:            
            call readmacro(nin,.false.)
            write(string,'(30a)') (word(j),j=1,30)
            read(string,*,err = 200,end = 200) fload(1,1,i)
c ...       valor de a:            
            call readmacro(nin,.false.)
            write(string,'(30a)') (word(j),j=1,30)
            read(string,*,err = 200,end = 200) fload(2,1,i)
c .....................................................................
c      
c ...       -kdu/dx = -H.(u-uext) 
         elseif (load(1,i) .eq. 3) then
c ...       valor de uext:            
            call readmacro(nin,.false.)
            write(string,'(30a)') (word(j),j=1,30)
            read(string,*,err = 200,end = 200) fload(1,1,i)
c ...       valor de H:            
            call readmacro(nin,.false.)
            write(string,'(30a)') (word(j),j=1,30)
            read(string,*,err = 200,end = 200) fload(2,1,i)
c .....................................................................
c      
         elseif (load(1,i) .eq. 6) then
c ...       numero de parcelas:            
            call readmacro(nin,.false.)
            write(string,'(30a)') (word(j),j=1,30)
            read(string,*,err = 200,end = 200) load(2,i)
            do k = 1, load(2,i)
               call readmacro(nin,.true.)
               write(string,'(30a)') (word(j),j=1,30)
               read(string,*,err = 200,end = 200) fload(1,k,i)
               call readmacro(nin,.false.)
               write(string,'(30a)') (word(j),j=1,30)
               read(string,*,err = 200,end = 200) fload(2,k,i)
               call readmacro(nin,.false.)
               write(string,'(30a)') (word(j),j=1,30)
               read(string,*,err = 200,end = 200) fload(3,k,i)
            enddo            
c .....................................................................
c      
c ...       kdu/dx = emiss * const(Stef-Boltz) *(u4-uext4)+H(uext-u)
         elseif (load(1,i) .eq. 7) then
c ...       Emissividade:            
            call readmacro(nin,.false.)
            write(string,'(30a)') (word(j),j=1,30)
            read(string,*,err = 200,end = 200) fload(1,1,i)
c ...       Constante de Stefan-Boltzmann:            
            call readmacro(nin,.false.)
            write(string,'(30a)') (word(j),j=1,30)
            read(string,*,err = 200,end = 200) fload(2,1,i)
c ...       constante de conveccao (h):            
            call readmacro(nin,.false.)
            write(string,'(30a)') (word(j),j=1,30)
            read(string,*,err = 200,end = 200) fload(3,1,i)
c ...       numero de parcelas (uext no tempo):            
            call readmacro(nin,.false.)
            write(string,'(30a)') (word(j),j=1,30)
            read(string,*,err = 200,end = 200) load(2,i)
            do k = 1, load(2,i)
c ...       valor de tempo: 
               call readmacro(nin,.true.)
               write(string,'(30a)') (word(j),j=1,30)
               read(string,*,err = 200,end = 200) fload(k,2,i)
c ...       valor de uext:
               call readmacro(nin,.false.)
               write(string,'(30a)') (word(j),j=1,30)
               read(string,*,err = 200,end = 200) fload(k,3,i)
            enddo   
c .....................................................................
c   
c ...   
         elseif (load(1,i) .eq. 8) then
c ...       numero de parcelas:            
            call readmacro(nin,.false.)
            write(string,'(30a)') (word(j),j=1,30)
            read(string,*,err = 200,end = 200) load(2,i)
            do k = 1, load(2,i)
               call readmacro(nin,.true.)
               write(string,'(30a)') (word(j),j=1,30)
               read(string,*,err = 200,end = 200) fload(k,1,i)
               call readmacro(nin,.false.)
               write(string,'(30a)') (word(j),j=1,30)
               read(string,*,err = 200,end = 200) fload(k,2,i)
            enddo     
c .....................................................................
c      
c ... forca distribuida constante no contorno
         elseif (load(1,i) .eq. 40) then
c ...       numero de parcelas:            
            call readmacro(nin,.false.)
            write(string,'(30a)') (word(j),j=1,30)
            read(string,*,err = 200,end = 200) load(2,i)
            do k = 1, load(2,i)
               call readmacro(nin,.false.)
               write(string,'(30a)') (word(j),j=1,30)
               read(string,*,err = 200,end = 200) fload(k,1,i)
            enddo           
         endif
         call readmacro(nin,.true.)
         write(string,'(12a)') (word(j),j=1,12)
  100 continue
      return
c ......................................................................      
  200 continue
      print*,'*** Erro na leitura das cargas !'
      stop
      end    
      subroutine rtermprop(numat,nin)
c **********************************************************************
c *                                                                    *
c *   Leitura das propriedades variaveis com a temperatura             *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *                                                                    *
c *    idl(6) - variavel auxuliar                                      *
c *    nin - arquivo de entrada                                        *
c *                                                                    *
c *   Parametros de saida:                                             *
c *                                                                    *
c *    nprop (i,j,k) - valor i da propriedade j do material tipo k     *
c *    eprop (j,k) - numero de termos da propriedade j do material tipo*
c *                 k                                                  *
c *                                                                    *
c **********************************************************************
      implicit none
      include 'string.fi'
      include 'termprop.fi'
      integer  nin,i,j,k,ma, numat,m
      character*30 string            
c ......................................................................
      call readmacro(nin,.true.)
      write(string,'(12a)') (word(i),i=1,12)
      do 100 while(string .ne. 'end')
c ...    numero do material:
         read(string,*,err = 200,end = 200) ma
         if(ma .lt. 1 .or. ma .gt. numat) goto 200 
         call readmacro(nin,.false.)
         write(string,'(12a)') (word(j),j=1,12)
c ...    numero da propriedade variavel no material ma:
         read(string,*,err = 200,end = 200) i
         call readmacro(nin,.false.)
         write(string,'(30a)') (word(j),j=1,30)
c ...    numero de pontos da poligonal:
         read(string,*,err = 200,end = 200) eprop(i,ma)
         k = eprop(i,ma)
         do m = 1, k
c ...       valor da temperatura no ponto j:            
            call readmacro(nin,.true.)
            write(string,'(30a)') (word(j),j=1,30)
            read(string,*,err = 200,end = 200) nprop(m,i,ma)
            call readmacro(nin,.false.)
            write(string,'(30a)') (word(j),j=1,30)
c ...       valor da propriedade no ponto j:  
            read(string,*,err = 200,end = 200) nprop(m+k,i,ma)          
         enddo
      call readmacro(nin,.true.)
      write(string,'(12a)') (word(j),j=1,12)
  100 continue
      return
c ......................................................................      
  200 continue
      print*,'*** Erro na leitura das propriedades variaveis !'
      stop
      end  
      subroutine boundc(nnode,ndf,id,f,u)
c **********************************************************************
c *                                                                    *
c *   Valores prescritos.                                              *
c *                                                                    *
c **********************************************************************
      implicit none
      integer nnode,ndf,id(ndf,*),i,j,k
      real*8 f(ndf,*),u(ndf,*)
c ......................................................................
      do 200 i = 1, nnode
         do 100 j = 1, ndf
            k = id(j,i)
            if(k .gt. 0) then
               u(j,i) = f(j,i)
            endif
  100    continue
  200 continue
c ......................................................................  
      return
      end      
      subroutine initc(nnode,ndf,id,u0,u)
c **********************************************************************
c *                                                                    *
c *   Valores prescritos.                                              *
c *                                                                    *
c **********************************************************************
      implicit none
      integer nnode,ndf,id(ndf,*),i,j,k
      real*8 u0(ndf,*),u(ndf,*)
c ......................................................................
      do 200 i = 1, nnode
         do 100 j = 1, ndf
            k = id(j,i)
            if(k .le. 0) then
               u(j,i) = u0(j,i)
            endif
  100    continue
  200 continue
c ......................................................................  
      return
      end   
      subroutine parameters(nnode,numel,numat,nen,ndf,ndm,nin)
c **********************************************************************
c *                                                                    *
c *   Parameters                                                       *
c *                                                                    *
c **********************************************************************
      implicit none
      include 'string.fi'
      character*12 string
      integer nnode,numel,numat,nen,ndf,ndft,ndm,nin,n,j,npi
      logical flag(6)
      character*6 macro(6)
      integer i,nmc 
      data macro/'nnode ','numel ','numat '
     .          ,'maxno ','ndf   ','dim   '/
      data nmc /6/
c ......................................................................  
      flag(1:nmc) = .false.          
      n   = 0
      call readmacro(nin,.true.)
      write(string,'(6a)') (word(j),j=1,6)                        
      do while (strl .ne. 0)
         if     (string .eq. 'nnode') then
            call readmacro(nin,.false.)
            write(string,'(12a)') (word(j),j=1,12)            
            read(string,*,err = 100,end = 100) nnode
            flag(1) = .true.
            n = n + 1
         elseif (string .eq. 'numel') then
            call readmacro(nin,.false.)
            write(string,'(12a)') (word(j),j=1,12)
            read(string,*,err = 100,end = 100) numel
            flag(2) = .true.
            n = n + 1
         elseif (string .eq. 'numat') then
            call readmacro(nin,.false.)
            write(string,'(12a)') (word(j),j=1,12)
            read(string,*,err = 100,end = 100) numat
            flag(3) = .true.
            n = n + 1
         elseif (string .eq. 'maxno') then
            call readmacro(nin,.false.)
            write(string,'(12a)') (word(j),j=1,12)           
            read(string,*,err = 100,end = 100) nen
            flag(4) = .true.
            n = n + 1
         elseif (string .eq. 'ndf') then
            call readmacro(nin,.false.)
            write(string,'(12a)') (word(j),j=1,12)            
            read(string,*,err = 100,end = 100) ndf
            flag(5) = .true.
            n = n + 1
         elseif (string .eq. 'dim') then
            call readmacro(nin,.false.)
            write(string,'(12a)') (word(j),j=1,12)         
            read(string,*,err = 100,end = 100) ndm
            flag(6) = .true.
            n = n + 1
         endif
         call readmacro(nin,.false.)
         write(string,'(6a)') (word(j),j=1,6)                 
      end do
      if(n .lt. nmc) goto 110
      return
c ......................................................................                        
  100 continue
      print*,'*** Erro na leitura das variaveis de controle !'
      stop       
  110 continue
        print*,'*** Erro na leitura das variaveis de controle !'
        do i = 1, nmc
          if(flag(i) .eqv. .false.) print *,"falta a macro ",macro(i)
        enddo
      stop       
c ......................................................................                  
      end
      subroutine readmacro(nin,newline)
c **********************************************************************
c *                                                                    *
c *   Subroutine READMACRO: le uma macro numa linha nova ou a partir   *
c *                         da coluna line_col de uma linha lida       *
c *                         anteriormente e armazenada em line(*)      *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *                                                                    *
c *   nin - numero do arquivo de entrada                               *
c *   newline = .true.  - nova linha deve ser lida                     *
c *           = .false. - macro deve ser lida a partir da coluna       *
c *                       line_col em line(*)                          *
c *                                                                    *
c **********************************************************************
      implicit none
      include 'string.fi'
      integer j,k,nin
      logical newline
c ......................................................................
      if(newline) then
         line_col = 1
         read(nin,'(200a1)',err = 200,end = 200) (line(j),j=1,200)
      endif
c ......................................................................      
      do j = 1, maxstrl
         word(j) = ' '
      enddo
      strl = 0
      if(line_col .gt. maxstrl) return
c ......................................................................      
      do while (line(line_col) .eq. ' ')
         line_col = line_col + 1
         if (line_col .gt. maxstrl) return         
      end do
c ......................................................................      
      do while ( (line(line_col) .ne. ' ') .and.
     .           (line(line_col) .ne. CHAR(13)) )
         strl = strl + 1
         line_col = line_col + 1
         if (line_col .gt. maxstrl) goto 100
      end do
c ......................................................................      
  100 continue      
      k = line_col-strl
      do j = 1, strl
         write(word(j),'(a)') line(k)
         k = k + 1
      end do      
      return
c ......................................................................      
  200 continue
      print*,'*** Erro na leitura do arquivo de dados !'
      stop             
c ......................................................................      
      end 
c ======================================================================      
c
c ======================================================================
c
c     Leitura da estrutura de dados do paralelo
c
c ======================================================================      
      subroutine read_par(nin,nnode,numel)
c **********************************************************************
c *                                                                    *
c *   READ_PAR                                                         *
c *   --------                                                         *
c *                                                                    *
c *   Le informacoes da paralelizacao geradas pelo pre-processador     *
c *   atraves do macro-comando read_par em rdat e cria os arranjos:    *
c *                                                                    *
c *      noLG(nnode) - mapa Local->Global de nos                       *
c *      noGL(nnoG)  - mapa Global-<Local de nos                       *
c *      elLG(numel) - mapa Local->Global de elementos                 *
c *      fmap(nnof)  - mapa Interface->Global de nos                   *
c *      idG(ndf,nnoG) - restricoes nodais globais                     *
c *                                                                    *
c *   Parâmetros de entrada:                                           *
c *   ----------------------                                           *
c *                                                                    *
c *                                                                    *
c *   Parâmetros de saida:                                             *
c *   -------------------                                              *
c *                                                                    *
c *                                                                    *
c **********************************************************************
      use Malloc
      implicit none
c      include 'mpif.h'
      include 'parallel.fi'
      include 'elementos.fi'
      integer nin,nnode,numel,i,j,k
      character*160 comando,string
c ......................................................................     
      read(nin,*) string,nnoG      
      read(nin,*) string,nelG      
c ......................................................................
c
c ... Metodo de sub-divisao de dominio
c
      ovlp  = .true.
      novlp = .false.
      read(nin,*)  comando
      if(comando(1:3) .eq. 'non') then
         ovlp  = .false.
         novlp = .true.
      endif
c ......................................................................
c
c ... Numeros de no's do processo, parte non-overlapping
c
      read(nin,'(a80)')  comando
      read(comando,*) string,nno1,
     .                string,nno2,
     .                string,nno3,
     .                string,nno4,
     .                string,nno1a
      if (ovlp)read(nin,*) string,numel_ov
c ......................................................................
c
c ... Numeros de elementos do processo, parte non-overlapping
c
      numel_nov = numel - numel_ov
c ......................................................................
c
c ... {noLG} = Mapa Local -> Global de no's
c
      read(nin,*) string
      i_noLG = alloc_4('noLG    ', 1, nnode)
      call Levet_4(ia(i_noLG),nnode,nin)
      i_noGL = alloc_4('noGL    ', 1, nnoG)
      call mzero(ia(i_noGL),nnoG)      
      do i = 1, nnode
         j = ia(i_noLG+i-1)
         ia(i_noGL+j-1) = i
      enddo        
c ......................................................................
c
c ... {noGL} = Mapa Global -> Local de no's
c
c      read(nin,*) string,nnoG
c      i_noGL = alloc_4('noGL    ', 1, nnoG)
c      call Lemtx_4(ia(i_noGL),1,nin)
c ......................................................................
c
c ..... Numeros de elementos pelo tipo
c
      read(nin,'(a)')  comando
      read(comando,*) string, nbar2(1)  , nbar2(2)  ,
     .                string, ntria3(1) , ntria3(2) ,
     .                string, nquad4(1) , nquad4(2) ,
     .                string, ntetra4(1), ntetra4(2),
     .                string, nhexa8(1) , nhexa8(2)
c ......................................................................
c
c ... {elLG} = Mapa Local -> Global de elementos
c
c      read(nin,*) string,nelG
      read(nin,*) string
      i_elLG = alloc_4('elLG    ', 1, numel)
      call Levet_4(ia(i_elLG),numel,nin)
c ......................................................................
c
c ... {rcvs0} = lista dos tamanhos dos blocos de nós {Vfi}
c
c      i_rcvs0 = alloc_4('rcvs0   ', 1, nprcs)
c      if (ovlp) then
c        neqfi=nno1a+nno2
c      else
c        neqfi=nno2+nno3
c      endif
c      call MPI_allgather(neqfi,1,MPI_INTEGER,ia(i_rcvs0),1,MPI_INTEGER,
c     .                   MPI_COMM_WORLD,ierr)
c ......................................................................
c
c ... {dspl0} = Ponteiro p/ início de cada bloco em {Vf}
c
c      i_dspl0 = alloc_4('dspl0   ', 1, nprcs)
c      ia(i_dspl0) = 0
c      do i = 2, nprcs
c        ia(i_dspl0+i-1) = ia(i_dspl0+i-2) + ia(i_rcvs0+i-2)
c      enddo
c ......................................................................
c
c ... {fmap0} = Mapa Interface -> Global de no's
c
c      read(nin,*) string,nnof
c      i_fmap0 = alloc_4('fmap0   ', 1, nnof)
c      call Levet_4(ia(i_fmap0),nnof,nin)
      read(nin,*) string,nnof1,nnof2
      i_fmap0 = alloc_4('fmap0   ', 1, nnof1+nnof2)
      call Levet_4(ia(i_fmap0),nnof1,nin)
      if(ovlp)call Levet_4(ia(i_fmap0+nnof1),nnof2,nin)
c ......................................................................
c
c ...   Restricoes Globais
c
c      if(ndf .gt. 0) then
c         read(nin,*) string
c         i_idG = alloc_4('idG     ', ndf, nnoG)
c         call Lemtx_4(ia(i_idG),ndf,nin)
c      endif
c      if (ndft .gt. 0) then
c         read(nin,*) string
c         i_idGt = alloc_4('idGt    ', ndft, nnoG)
c         call Lemtx_4(ia(i_idGt),ndft,nin)    
c      endif
c ......................................................................
c
c ... Estruturas de comunicacao sendr
c
c      read(nin,*) string !nviz
c      read(nin,*) nviz
c      i_ranks  = alloc_4('ranks   ', 1, nprcs)
c      read(nin,*) string !ranks
c      read(nin,*) (ia(i_ranks+j), j = 0,nviz-1)
      nviz1 = 0
      nviz2 = 0
      read(nin,*) string !nviz
      read(nin,*) nviz1
      nviz2 = nviz1
      if(ovlp)read(nin,*) nviz2
c ....
      i_ranks  = alloc_4('ranks   ', 1, nviz1+nviz2)
      read(nin,*) string !ranks
      read(nin,*) (ia(i_ranks+j), j = 0,nviz1-1)
      if(ovlp)read(nin,*) (ia(i_ranks+nviz1+j), j = 0,nviz2-1)
c ....
      i_rcvs0 = alloc_4('rcvs0   ', 1, nviz1+nviz2)
      read(nin,*) string !rcvs
      read(nin,*) (ia(i_rcvs0+j), j = 0,nviz1-1)
      if(ovlp)read(nin,*) (ia(i_rcvs0+nviz1+j), j = 0,nviz2-1)
c ....
      i_dspl0 = alloc_4('dspl0   ', 1, nviz1+nviz2)
      read(nin,*) string !dspl
      read(nin,*) (ia(i_dspl0+j), j = 0,nviz1-1)
      if(ovlp)read(nin,*) (ia(i_dspl0+nviz1+j), j = 0,nviz2-1)
c ... end parallel
      read(nin,*) string
c ......................................................................
      return
      end
c ......................................................................
      subroutine Levet_4(array,lin,nin)
c **********************************************************************
c *                                                                    *
c *   Levet_4                                                          *
c *   -------                                                          *
c *                                                                    *
c *   Le arranjo inteiro                                               *
c *                                                                    *
c *                                                                    *
c *   Parâmetros de entrada:                                           *
c *   ----------------------                                           *
c *                                                                    *
c *                                                                    *
c *   Parâmetros de saida:                                             *
c *   -------------------                                              *
c *                                                                    *
c *                                                                    *
c **********************************************************************
      implicit none
      integer array(*),lin,nin,i
c ......................................................................
      do i=1,lin
        read(nin,*)array(i)
      enddo
c ......................................................................
      return
      end
c ......................................................................
      subroutine Lemtx_4(array,lin,nin)
c **********************************************************************
c *                                                                    *
c *   Lemtx_4                                                          *
c *   -------                                                          *
c *                                                                    *
c *   Le arranjo inteiro                                               *
c *                                                                    *
c *                                                                    *
c *   Parâmetros de entrada:                                           *
c *   ----------------------                                           *
c *                                                                    *
c *                                                                    *
c *   Parâmetros de saida:                                             *
c *   -------------------                                              *
c *                                                                    *
c *                                                                    *
c **********************************************************************
      implicit none
      integer array(lin,*),lin,nin,i,k
      character*80 comando
c ......................................................................
      read(nin,'(a80)') comando
      do while(comando(1:3) .ne. 'end')
        read(comando,*) k,(array(i,k),i=1,lin)
        read(nin,'(a80)') comando
      enddo
c ......................................................................
      return
      end
c *********************************************************************
c
c **********************************************************************
c * Data de criacao    : 28/03/2016                                    *
c * Data de modificaco : 09/04/2016                                    * 
c * ------------------------------------------------------------------ *
c * MK_BOUND_QUAD: gera as restricoes nos deslocamente nos pontos      *  
c *              intermediarios                                        *
c * ------------------------------------------------------------------ *
c * Parâmetros de entrada:                                             *
c * ------------------------------------------------------------------ *
c * id(ndf,*)  - restricoes nos vertices                               *
c * numel      - numero de elementos                                   *
c * ndf1       - grau de liberdade liberdade de deslocamento           *
c * ndf2       - grau de liberdade liberdade total do problema         *
c * nen        - numero de nos elementos quadraticos                   *
c * ------------------------------------------------------------------ *
c * Parâmetros de saida:                                               *
c * ------------------------------------------------------------------ *
c * id(ndf1,*)  - restricoes atualizadas                               *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ *
c * mecanico     : ndf1 = x, y e z  |  ndf2 = x, y e z                 *  
c * poromecanico : ndf1 = x, y,e z  |  ndf2 = x, y, z e p              *
c **********************************************************************
      subroutine mk_bound_quad(id,el,numel,ndf1,ndf2,nen)
      implicit none
      integer maxEdge
      parameter (maxEdge = 12) 
c ...
      integer i,j,k
      integer el(nen+1,*)
      integer iEdge(3,maxEdge)
      integer numel,ndf1,ndf2,nedge,no1,no2,no3,nen
      integer id(ndf2,*)
c ...
      nedge = 0 
c ... tetraedros de 10 nos 
      if( nen .eq. 10 ) then
        nedge =  6
        call tetra10edgeNod(iEdge) 
c ... hexaedros de 20 nos 
      else if( nen .eq. 20 ) then
        nedge = 12
        call hexa20edgeNod(iEdge) 
      endif
c ...
      do i = 1, numel
        do j = 1, nedge
c ... no vertices
          no1     = el(iEdge(1,j),i)
          no2     = el(iEdge(2,j),i)
c ... no central
          no3     = el(iEdge(3,j),i)
          do k = 1, ndf1
            if( id(k,no1) .eq. 1 .and. id(k,no2) .eq. 1) then
              id(k,no3) = 1
            endif
          enddo
        enddo
      enddo
c .....................................................................
c
c ...
c     do i = 1, 27
c       print*,i,id(1:4,i)
c     enddo
c .....................................................................  
      return
      end
c *********************************************************************
c
c **********************************************************************
c * Data de criacao    : 28/03/2016                                    *
c * Data de modificaco :                                               * 
c * ------------------------------------------------------------------ *
c * MK_FORCES_QUAD: gera as valores das cargas  nos pontos             *
c * intermediarios                                                     *
c * ------------------------------------------------------------------ *
c * Parâmetros de entrada:                                             *
c * ------------------------------------------------------------------ *                                                                  *
c * id(ndf,*)  - restricoes atualizadas                                *
c * f (ndf,*)  - valor das cargas nos vertices                         *
c * el(nen+1,*)- conectividade nodal                                   *
c * numel      - numero de elementos                                   *
c * ndf        - grau de liberdade                                     *
c * nen        - numero de nos elementos quadraticos                   *
c * ------------------------------------------------------------------ *
c * Parâmetros de saida:                                               *
c * ------------------------------------------------------------------ *                                          *
c * id(ndf,*)  - restricoes atualizadas                                *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ *                                                                    *
c **********************************************************************
      subroutine mk_forces_quad(id,f,el,numel,ndf,nen)
      implicit none
      integer maxEdge
      parameter (maxEdge = 12) 
c ...
      integer i,j,k
      integer el(nen+1,*)
      integer iEdge(3,maxEdge)
      integer id(ndf,*)
      integer numel,ndf,nedge,no1,no2,no3,nen
      real*8  f(ndf,*)
c ...
      nedge = 0  
c ... tetraedros de 10 nos 
      if( nen .eq. 10 ) then
        nedge =  6
        call tetra10edgeNod(iEdge) 
c ... hexaedros de 20 nos 
      else if( nen .eq. 20 ) then
        nedge = 12
        call hexa20edgeNod(iEdge) 
      endif
c ...
      do i = 1, numel
        do j = 1, nedge
c ... no vertices
          no1     = el(iEdge(1,j),i)
          no2     = el(iEdge(2,j),i)
c ... no central
          no3     = el(iEdge(3,j),i)
          do k = 1, ndf-1
            if( id(k,no1) .eq. 1 .and. id(k,no2) .eq. 1) then
              f(k,no3) = 0.5d0*(f(k,no1) +  f(k,no2))
            endif 
          enddo
        enddo
      enddo
c .....................................................................
c
c ...
c     do i = 1, 70
c       print*,i,f(1:4,i)
c     enddo
c .....................................................................  
      return
      end
c *********************************************************************
c
c **********************************************************************
c * Data de criacao    : 28/03/2016                                    *
c * Data de modificaco :                                               * 
c * ------------------------------------------------------------------ *
c * MK_INITIAL_QUAD:gera os valores iniciais nos pontos intermediarios *
c * ------------------------------------------------------------------ *
c * Parâmetros de entrada:                                             *
c * ------------------------------------------------------------------ *
c * f(ndf,*)   - valor das cargas nos vertices                         *
c * el(nen+1,*)- conectividade nodal                                   *
c * numel      - numero de elementos                                   *
c * ndf        - grau de liberdade                                     *
c * nen        - numero de nos elementos quadraticos                   *
c * ------------------------------------------------------------------ *
c * Parâmetros de saida:                                               *
c * ------------------------------------------------------------------ *
c * f(ndf,*)   - valor das cargas atualizadas                          *
c * ------------------------------------------------------------------ * 
c *  OBS:                                                              *
c * ------------------------------------------------------------------ * 
c **********************************************************************
      subroutine mk_initial_quad(f,el,numel,ndf,nen)
      implicit none
      integer maxEdge
      parameter (maxEdge = 12) 
c ...
      integer i,j,k
      integer el(nen+1,*)
      integer iEdge(3,maxEdge)
      integer numel,ndf,nedge,no1,no2,no3,nen
      real*8  f(ndf,*)
c ...
      nedge = 0  
c ... tetraedros de 10 nos 
      if( nen .eq. 10 ) then
        nedge =  6
        call tetra10edgeNod(iEdge) 
c ... hexaedros de 20 nos 
      else if( nen .eq. 20 ) then
        nedge = 12
        call hexa20edgeNod(iEdge) 
      endif
c ...
      do i = 1, numel
        do j = 1, nedge
c ... no vertices
          no1     = el(iEdge(1,j),i)
          no2     = el(iEdge(2,j),i)
c ... no central
          no3     = el(iEdge(3,j),i)
          do k = 1, ndf
            f(k,no3) = 0.5d0*(f(k,no1) +  f(k,no2))
          enddo
        enddo
      enddo
c .....................................................................
c
c ...
c     do i = 1, 70
c       print*,i,f(1:4,i)
c     enddo
c .....................................................................  
      return
      end
c **********************************************************************
c
c **********************************************************************
c *                                                                    *
c *   GET_PRES : obtem as presoes                                      *
c *                                                                    *
c *   Parâmetros de entrada:                                           *
c *   ----------------------                                           *
c *                                                                    *
c *   u (ndf,*)  - deslocamento e pressoes                             *
c *   pres(*)    - indefinido                                          *
c *   nnodev     - nos de vertices                                     *
c *   ndf        - grau de liberdade                                   *
c *                                                                    *
c *   Parâmetros de saida:                                             *
c *   -------------------                                              *
c *   pres(*)    - pressoes                                            *
c *                                                                    *
c **********************************************************************
      subroutine get_pres(u,pres,nnodev,ndf)
      implicit none
      integer ndf,i,nnodev
      real*8 u(ndf,*),pres(*)
c ...
      do i = 1, nnodev
        pres(i) = u(ndf,i)
      enddo
c .....................................................................  
      return
      end
c *********************************************************************
c
c *********************************************************************
c * read_config : leitura das configuracoes basicas de excucao        *
c * ------------------------------------------------------------------*
c * Parametros de entrada :                                           *
c * ----------------------------------------------------------------- *
c * maxmem        - memoria do vetor de trabalho                      *
c * omp_elmt  - flag do openmp na fase de elemento                    *
c * nth_elmt  - numero de threads usado na fase de elemento           *
c * omp_solv  - flag do openmp na fase do solver                      *
c * nth_solv  - numero de threads usado do solver                     *
c * reord     - reordanaco do sistema de equacao                      *
c * bvtk      - saida binario para o vtk legacy                       *
c * nin       - arquivo de entrada                                    *
c * ----------------------------------------------------------------- *
c * Parametros de saida :                                             *
c * ----------------------------------------------------------------- *
c * ----------------------------------------------------------------- *
c *********************************************************************
      subroutine read_config(maxmem  
     .                      ,omp_elmt,omp_solver
     .                      ,nth_elmt,nth_solver
     .                      ,reord   ,bvtk    
     .                      ,nin)
      implicit none
      include 'string.fi'
      character*15 string,macro(7)
      integer*8 maxmem
      integer nth_elmt,nth_solver
      logical omp_elmt,omp_solver,r(7),reord,bvtk
      integer n,j,nmacro
      integer nin
      data nmacro /7/
      data macro/'memory         ','omp_elmt       ','omp_solver     ',
     .           'nth_elmt       ','nth_solver     ','reord          ',
     .           'bvtk           '/
c .....................................................................
      n      = 0
      r(1:nmacro) = .false.
      call readmacro(nin,.true.)
      write(string,'(15a)') (word(j),j=1,12)
      do while (strl .ne. 0)
         if     (string .eq. macro(1)) then
            call readmacro(nin,.false.)
            write(string,'(15a)') (word(j),j=1,15)            
            read(string,*,err = 100,end = 100) maxmem
c ... convertendo de Mbytes para para numero de inteiros e 4 bytes
            maxmem = (maxmem*1024*1024)/4
            n = n + 1
            r(1) = .true.
         elseif (string .eq. macro(2)) then
            call readmacro(nin,.false.)
            write(string,'(15a)') (word(j),j=1,15)
            read(string,*,err = 100,end = 100) omp_elmt
            n = n + 1
            r(2) = .true.
         elseif (string .eq. macro(3)) then
            call readmacro(nin,.false.)
            write(string,'(15a)') (word(j),j=1,15)           
            read(string,*,err = 100,end = 100) omp_solver
            n = n + 1
            r(3) = .true.
         elseif (string .eq. macro(4)) then
            call readmacro(nin,.false.)
            write(string,'(15a)') (word(j),j=1,15)            
            read(string,*,err = 100,end = 100) nth_elmt
            n = n + 1
            r(4) = .true.
         elseif (string .eq. macro(5)) then
            call readmacro(nin,.false.)
            write(string,'(15a)') (word(j),j=1,15)            
            read(string,*,err = 100,end = 100) nth_solver
            n = n + 1
            r(5) = .true.
         elseif (string .eq. macro(6)) then
            call readmacro(nin,.false.)
            write(string,'(15a)') (word(j),j=1,15)            
            read(string,*,err = 100,end = 100) reord          
            n = n + 1
            r(6) = .true.
         elseif (string .eq. macro(7)) then
            call readmacro(nin,.false.)
            write(string,'(15a)') (word(j),j=1,15)            
            read(string,*,err = 100,end = 100) bvtk          
            n = n + 1
            r(7) = .true.
         endif 
         call readmacro(nin,.false.)
         write(string,'(15a)') (word(j),j=1,15)                 
      end do
c ...
      call readmacro(nin,.true.)
      call readmacro(nin,.false.)
c ......................................................................
c
c ...
      if(n .lt. 7) goto 100
      return
c ......................................................................                        
  100 continue
      print*,'*** Erro na leitura das variaveis de controle !'
      do j = 1, nmacro
        if(.NOT.r(j)) print*,macro(j), ' faltando!!'
      enddo
      stop       
c ......................................................................                  
      end
c **********************************************************************

