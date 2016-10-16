c *********************************************************************
c * Data de criacao    : 12/12/2015                                   *
c * Data de modificaco : 03/04/2016                                   * 
c * ------------------------------------------------------------------*    
c * WRITE_MESH_GEO: escreve a malha no vtk                            *
c * ------------------------------------------------------------------*
c * Parametros de entrada :                                           *
c * ----------------------------------------------------------------- *
c * el(nen+1,*) - conectividade com material                          *
c *  x(ndm,*)   - coordenadas                                         *
c * nnode       - numero de nos de vertices                           *
c * numel       - numero de elementos                                 *
c * nen         - numero de nos por elementos                         *
c * ndm         - numero de dimensoes                                 *
c * filein      - prefix do arquivo de saida                          *
c * bvtk        - true BINARY vtk false ASCII vtk                     *
c * legacy      - true (formato padrão .vtk) false (formato xlm .vtu)*
c * nout        - arquivo de saida                                    *
c * ----------------------------------------------------------------- *    
c * Parametros de saida :                                             *
c * ----------------------------------------------------------------- *
c * ----------------------------------------------------------------- * 
c * OBS:                                                              *
c * ----------------------------------------------------------------- *         
c *********************************************************************
      subroutine write_mesh_geo(el     ,x  ,nnode ,numel
     .                         ,nen    ,ndm,filein,bvtk
     .                         ,legacy ,nout)
c ===
      use Malloc 
      implicit none
c ... variaveis da malha      
      integer nnode,numel,nen,ndm
      integer el(nen+1,numel)
      real*8  x(ndm,nnode)
      integer nel,nno
c ... locais     
      integer*8 i_p
      data i_p/1/
      character*15 aux1
      character*30 aux
c ... variaveis dums
      real*8 ddum
      real*4 fdum
      integer idum
c ... arquivo      
      integer nout
      character*80 fileout,name,filein
      logical bvtk,legacy
      integer cod,cod2,gdl
c =====================================================================
c
c ===
      if(legacy) then
        fileout = name(filein,0,0)
      else  
        fileout = name(filein,0,1)
      endif  
      if(bvtk)then
        open(unit=nout,file=fileout,access='stream'
     .      ,form='unformatted',convert='big_endian')
      else
        open(unit=nout,file=fileout)
      endif  
c =====================================================================
c
c === cabecalho
      if(legacy) then
        write(aux,'(30a)')"Geometria" 
        call head_vtk(aux,bvtk,ddum,idum,.false.,nout)
      else  
        call head_vtu(nnode,numel,bvtk,nout) 
      endif  
c =====================================================================
c
c === Coordenadas
      if(legacy) then
        call coor_vtk(x,nnode,ndm,bvtk,nout)
      else  
        call coor_vtu(x,nnode,ndm,bvtk,nout)
      endif
c =====================================================================
c
c === Elementos
      if(legacy) then
        call elm_vtk(el,numel,nen,bvtk,nout)
      else  
        call elm_vtu(el,numel,nen,bvtk,nout)
      endif  
c =====================================================================
c
c === cell 
      if(legacy) then
        call cell_data_vtk(numel,bvtk,nout)
      else  
        call cell_data_vtu(bvtk,nout)
      endif  
c ... materiais      
      i_p = alloc_4('p       ', 1,numel)
      do nel = 1, numel
        ia(i_p+nel-1) = el(nen+1,nel)
      enddo
      write(aux1,'(15a)')"mat" 
c ... cod = 1 variaveis interias
      cod = 1
      gdl = 1
      if(legacy) then
        call cell_prop_vtk(ia(i_p),fdum,ddum,numel,aux1,cod,gdl,bvtk
     .                    ,nout)
      else
        call cell_prop_vtu(ia(i_p),fdum,ddum,numel,aux1,cod,gdl,bvtk
     .                    ,nout)
      endif 
      i_p = dealloc('p       ')
c .....................................................................
c
c ... elIdG
      i_p = alloc_4('p       ', 1,numel)
      do nel = 1, numel
        ia(i_p+nel-1) = nel
      enddo
      write(aux1,'(15a)')"elIdG" 
c ... cod = 1 variaveis interias
      cod = 1
      gdl = 1
      if(legacy) then
        call cell_prop_vtk(ia(i_p),fdum,ddum,numel,aux1,cod,gdl,bvtk
     .                    ,nout)
      else
        call cell_prop_vtu(ia(i_p),fdum,ddum,numel,aux1,cod,gdl,bvtk
     .                    ,nout)
      endif 
      i_p = dealloc('p       ')
c .....................................................................
c =====================================================================
c
c === nos  
c ...       
      if(legacy) then
        call point_data_vtk(nnode,bvtk,nout)
      else  
        call point_data_vtu(bvtk,nout)
      endif  
      i_p = alloc_4('p       ', 1,nnode)
      do nno = 1, nnode
        ia(i_p+nno-1) = nno    
      enddo
      write(aux1,'(15a)')"noIdG"
c ... gdb graus de liberdade
c     cod  1 escalar
c     cod2 1 int(4bytes) 
      gdl =  1
      cod =  1
      cod2 = 1
      if(legacy) then
        call point_prop_vtk(ia(i_p),fdum,ddum,nnode,aux1,ndm,gdl,cod,cod2
     .                    ,bvtk,nout)
      else
        call point_prop_vtu(ia(i_p),fdum,ddum,nnode,aux1,ndm,gdl,cod,cod2
     .                    ,bvtk,nout)
      endif
      i_p = dealloc('p       ')

      if(legacy .eqv. .false.) then
        call point_data_finalize_vtu(bvtk,nout)
        call finalize_vtu(bvtk,nout)
      endif
c .....................................................................
c =====================================================================
      close(nout)
      return
      end
c =====================================================================
c *********************************************************************
c
c *********************************************************************
c * Data de criacao    : 03/04/2016                                   *
c * Data de modificaco : 25/04/2016                                   * 
c * ------------------------------------------------------------------* 
c * WRITE_MESH_GEO_BC: escreve a malha global com carregamento        *
c * no formato vtk.                                                   *
c * ----------------------------------------------------------------- *
c * Parametros de entrada :                                           *
c * ----------------------------------------------------------------- *
c * el(nen+1,*)  - conectividade com material                         *
c *  x(ndm,*)    - coordenadas                                        *
c * ie           - tipo do elemento                                   *
c * id(ndf,*)    - restricoes mecanico                                *
c * f(ndf,*)     - carregamento mecanicos                             *
c * u0(ndf,*)    - valores iniciais                                   *
c * nload(ndf,*) - cargas nos nos                                     *
c * eload(7,*)   - cargas nos elementos                               *
c * nnode        - numero de nos de vertices                          *
c * numel        - numero de elementos                                *
c * ndf          - graus de liberdade mecanicas                       *
c * ntn          - numero de termos no tensor de tensoes              *
c * nen          - numero de nos por elementos                        *
c * ndm          - numero de dimensoes                                *
c * filein       - prefix do arquivo de saida                         *
c * bvtk         - true BINARY vtk false ASCII vtk                    *
c * macros       - macros lidas pela rdata                            *
c * legacy       - true (formato padrão .vtk) false (formato xlm .vtu)*
c * nout         - arquivo de saida principal                         *
c * nelemtload   - arquivo de saida para face com cargas(elloads)     *
c * ----------------------------------------------------------------- *
c * Parametros de saida :                                             *
c * ----------------------------------------------------------------- *
c *********************************************************************
      subroutine write_mesh_geo_bc(el   ,x      ,ie
     .                         ,id   ,f         ,u0    
     .                         ,tx0  ,nload     ,eload
     .                         ,nnode,numel     ,ndf   ,ntn
     .                         ,nen  ,ndm       ,filein
     .                         ,bvtk ,macros    ,legacy
     .                         ,nout ,nelemtload)
c ===
      use Malloc 
      implicit none
      include 'elementos.fi'
c ... variaveis da malha      
      integer nnode,numel,nen,ndm
      integer el(nen+1,numel),ie(*)
      real*8  x(ndm,nnode)
      integer nel,nno
c ... variaveis do problema
c     poro mecanico      
      integer ndf,ntn,ntn1
      integer id(ndf,*),nload(ndf,*),eload(*)
      real*8  f(ndf,*),u0(ndf,*),tx0(6,*)
c ... locais     
      integer*8 i_p,i_b1,i_b2
      data i_p/1/ i_b1/1/ i_b2/1/
      character*15 aux1
      character*30 aux
c ... faces                                          
      integer line_face,tria_face,quad_face,nface
c ... macrp lida pela rdat      
      character*15 macros(*)
c ... variaveis dums
      real*8 ddum
      real*4 fdum
      integer idum
c ...
      character*8 malloc_name
c ... arquivo      
      integer nout,nelemtload
      character*80 fileout,fileelmtload,name,filein
      logical bvtk,legacy
      integer cod,cod2,gdl
c =====================================================================
c ... auxiliares
      integer nmacros,nmc,j,i
      character*15 macro1(18),macro2(18),rc
c ......................................................................
      data macro1/'              ','nodalforces    ','nodalloads     ',
     .            'constraindisp ','               ','               ',
     .            '              ','               ','               ',
     .            '              ','               ','               ',
     .            '              ','               ','               ',
     .            '              ','initialstress  ','end            '/
c      
      data macro2/'elmtloads      ','               ','               ',
     .            '               ','               ','               ',
     .            '               ','               ','               ',
     .            '               ','               ','               ',
     .            '               ','               ','               ',
     .            '               ','               ','end            '/
      data nmc /18/
c ......................................................................

c
c ===
      if(legacy) then
        fileout = name(filein,0,5)
      else
        fileout = name(filein,0,6)
      endif
      if(bvtk)then
        open(unit=nout,file=fileout,access='stream'
     .      ,form='unformatted',convert='big_endian')
      else
        open(unit=nout,file=fileout)
      endif  
c =====================================================================
c
c === cabecalho
      if(legacy) then
        write(aux,'(30a)')"Geometria lida" 
        call head_vtk(aux,bvtk,ddum,idum,.false.,nout)
      else
        call head_vtu(nnode,numel,bvtk,nout)
      endif
c =====================================================================
c
c === Coordenadas
      if(legacy) then
        call coor_vtk(x,nnode,ndm,bvtk,nout)
      else  
        call coor_vtu(x,nnode,ndm,bvtk,nout)
      endif  
c =====================================================================
c
c === Elementos
      if(legacy) then
        call elm_vtk(el,numel,nen,bvtk,nout)
      else  
        call elm_vtu(el,numel,nen,bvtk,nout)
      endif  
c =====================================================================
c
c === point data
      if(legacy) then
        call point_data_vtk(nnode,bvtk,nout)
      else  
        call point_data_vtu(bvtk,nout)
      endif  
c .....................................................................
c
c ... nos com numeracao global     
      write(aux1,'(15a)')"noIdG"
      gdl =  1   
      cod =  1
      cod2 = 1
      malloc_name ='p'
      i_p = alloc_4(malloc_name, 1,nnode)
      do nno = 1, nnode
        ia(i_p+nno-1) = nno
      enddo
      if(legacy) then
        call point_prop_vtk(ia(i_p),fdum,ddum,nnode,aux1,ndf,gdl
     .                     ,cod,cod2,bvtk,nout)
      else
        call point_prop_vtu(ia(i_p),fdum,ddum,nnode,aux1,ndf,gdl
     .                     ,cod,cod2,bvtk,nout)
      endif
      i_p = dealloc(malloc_name)
c .....................................................................
c
c ...
      nmacros = 0 
  100 continue
      nmacros = nmacros + 1 
      write(rc,'(15a)')macros(nmacros)
      do 200 j = 1, nmc
        if(rc .eq. macro1(j)) go to 300
  200 continue
      goto 100
c .....................................................................
  300 continue
c .....................................................................
      goto( 400, 450, 500   !             ,nodalforces  ,nodalloads
     .    , 550, 600, 650   !constraindisp,             ,        
     .    , 700, 750, 800   !             ,             ,
     .    , 850, 900, 950   !             ,             ,             
     .    ,1000,1050,1100   !             ,             ,          
     .    ,1150,1200,1250)j !             ,initialstress,end
c ......................................................................
c
c ...                  
  400 continue
      goto 100
c .....................................................................
c
c ... nodalforces  
  450 continue
c ... gdb graus de liberdade
c     cod  2 vetorial
c     cod2 3 real(8bytes) 
      write(aux1,'(15a)')"forces"
      gdl =  ndf
      cod =  1
      cod2 = 3
      if(legacy) then
        call point_prop_vtk(idum,fdum,f,nnode,aux1,ndf,gdl,cod,cod2
     .                    ,bvtk,nout)
      else
        call point_prop_vtu(idum,fdum,f,nnode,aux1,ndf,gdl,cod,cod2
     .                    ,bvtk,nout)
      endif
      goto 100
c .....................................................................      
c 
c ... nodalloads
  500 continue
c ... gdb graus de liberdade
c     cod  2 vetorial
c     cod2 1 interio(4bytes) 
      write(aux1,'(15a)')"nodalloads"
      gdl =  ndf
      cod =  2
      cod2 = 1
      if(legacy) then
        call point_prop_vtk(nload,fdum,ddum,nnode,aux1,ndf,gdl,cod,cod2
     .                    ,bvtk,nout)
      else
        call point_prop_vtu(nload,fdum,ddum,nnode,aux1,ndf,gdl,cod,cod2
     .                    ,bvtk,nout)
      endif
      goto 100
c ......................................................................
c
c ... constriandisp  
  550 continue
c ... gdb graus de liberdade
c     cod  1 scalar      
c     cod2 1 int(4bytes)
      malloc_name ='id_bc'
      i_p = alloc_4(malloc_name,ndf,nnode)
      call mzero(ia(i_p),nnode*ndf)
c ... identifica os nos prescritos   
      call mk_dot_bc(id,ia(i_p),nnode,ndf)
c .....................................................................      
      write(aux1,'(15a)')"constraindisp"
      gdl =  ndf 
      cod =  1
      cod2 = 1
      if(legacy) then
        call point_prop_vtk(ia(i_p),fdum,ddum,nnode,aux1,ndf,gdl,cod
     .                     ,cod2   ,bvtk,nout)
      else
        call point_prop_vtu(ia(i_p),fdum,ddum,nnode,aux1,ndf,gdl,cod
     .                     ,cod2   ,bvtk,nout)
      endif
c .....................................................................
c
c ... dealloc      
      i_p = dealloc(malloc_name)
c .....................................................................
      goto 100
c ......................................................................
c
c ... 
  600 continue
  650 continue
  700 continue
  750 continue
  800 continue
  850 continue
  900 continue
  950 continue
 1000 continue
 1050 continue
 1100 continue
 1150 continue
      goto 100
c .....................................................................
c
c ...
 1200 continue
c ... gerando o tensor completo 
      if( ntn .eq. 6 ) then
        malloc_name ='tensor'
        ntn1        = 9
        i_p         = alloc_8(malloc_name,ntn1 ,nnode)
        call make_full_tensor(tx0,ia(i_p),nnode,ntn, ntn1 )
      endif
c .....................................................................
c
c ... initialstress          
      write(aux1,'(15a)')"initialstress"
c ... gdb graus de liberdade
c     cod  3 tensor     
c     cod2 3 real(8bytes) 
      gdl  =  ntn1              
      cod  =  3
      cod2 =  3
      if(legacy)then
        call point_prop_vtk(idum,fdum,ia(i_p),nnode,aux1,ndm,gdl
     .                    ,cod ,cod2,bvtk   ,nout)
      else
        call point_prop_vtu(idum,fdum,ia(i_p),nnode,aux1,ndm,gdl
     .                    ,cod ,cod2,bvtk   ,nout)
      endif
c .....................................................................
c
c ... dealloc
      i_p = dealloc(malloc_name)
c .....................................................................
      goto 100
c .....................................................................
c
c ... end      
 1250 continue
      if(legacy .eqv. .false.) call point_data_finalize_vtu(bvtk,nout)
c .....................................................................
c =====================================================================
c
c === cell data
      if(legacy) then
        call cell_data_vtk(numel,bvtk,nout)
      else
        call cell_data_vtu(bvtk,nout)
      endif
c ... materiais 
      malloc_name ='p'
      i_p = alloc_4(malloc_name, 1,numel)
      do nel = 1, numel
        ia(i_p+nel-1) = el(nen+1,nel)
      enddo
      write(aux1,'(15a)')"mat" 
c ... cod = 1 variaveis interias
      cod = 1
      if(legacy)then
        call cell_prop_vtk(ia(i_p),fdum,ddum,numel,aux1,cod,1,bvtk
     .                    ,nout)
      else
        call cell_prop_vtu(ia(i_p),fdum,ddum,numel,aux1,cod,1,bvtk
     .                    ,nout)
      endif
      i_p = dealloc(malloc_name)
c ... elementos com numeracao global     
      write(aux1,'(15a)')"elIdG"
      gdl =  1   
      cod =  1
      cod2 = 1
      malloc_name ='p'
      i_p = alloc_4(malloc_name, 1,numel)
      do nel = 1, numel
        ia(i_p+nel-1) = nel
      enddo
      if(legacy) then
        call point_prop_vtk(ia(i_p),fdum,ddum,numel,aux1,ndf,gdl
     .                     ,cod,cod2,bvtk,nout)
      else
        call point_prop_vtu(ia(i_p),fdum,ddum,numel,aux1,ndf,gdl
     .                     ,cod,cod2,bvtk,nout)
      endif
      i_p = dealloc(malloc_name)
c .....................................................................
c =====================================================================
c ...
      nmacros = 0 
  101 continue
      nmacros = nmacros + 1 
      write(rc,'(15a)')macros(nmacros)
      do 201 j = 1, nmc
        if(rc .eq. macro2(j)) go to 301
  201 continue
      goto 101
c .....................................................................
  301 continue
      goto( 401, 451, 501     !elmtloads    ,             ,             
     .    , 551, 601, 651     !             ,             ,        
     .    , 701, 751, 801     !             ,             ,
     .    , 851, 901, 951     !             ,             ,             
     .    ,1001,1051,1101     !             ,             ,          
     .    ,1151,1201,1251)j   !             ,             ,   
c .....................................................................
c
c ... elmtthermloads
  401 continue
c ... gdb graus de liberdade
c     cod  1 interio(4bytes)
      write(aux1,'(15a)')"elmtloads"
      gdl =  7   
      cod =  1
      if(legacy)then
        call cell_prop_vtk(eload,fdum,ddum,numel,aux1,cod,gdl,bvtk
     .                    ,nout)
      else
        call cell_prop_vtu(eload,fdum,ddum,numel,aux1,cod,gdl,bvtk
     .                    ,nout)
      endif
c .....................................................................
c
c ... arquivo auxiliar com as faces
      if(legacy) then
        fileelmtload = name(filein,0,7)
      else
        fileelmtload = name(filein,0,8)
      endif
      if(bvtk)then
        open(unit=nelemtload,file=fileelmtload,access='stream'
     .      ,form='unformatted',convert='big_endian')
      else
        open(unit=nelemtload,file=fileelmtload )
      endif
c .....................................................................
c
c ... criando a estrutura para as faces  
      malloc_name = 'face'
      i_p  = alloc_4(malloc_name, max_face*max_no_face,numel)
      malloc_name = 'carga'
      i_b1 = alloc_4(malloc_name,       1,numel)
      malloc_name = 'faceType'
      i_b2 = alloc_4(malloc_name,       1,numel)
      call make_face(el        ,ie      ,eload   ,ia(i_p)  ,ia(i_b1) 
     .              ,ia(i_b2)  ,numel   ,nen     ,max_face ,max_no_face
     .              ,line_face,tria_face,quad_face)
c .....................................................................
c
c ... cabecalho
      nface = line_face + tria_face + quad_face 
      write(aux,'(30a)')'face com carregamento'
      if(legacy) then
        call head_vtk(aux,bvtk,ddum,idum,.false.,nelemtload)
      else
        call head_vtu(nnode,nface,bvtk,nelemtload)
      endif
c .....................................................................
c
c ... Coordenadas
      if(legacy) then
        call coor_vtk(x,nnode,ndm,bvtk,nelemtload)
      else  
        call coor_vtu(x,nnode,ndm,bvtk,nelemtload)
      endif  
c .....................................................................
c
c ... faces
      if(legacy) then
        call face_vtk(ia(i_p)    ,ia(i_b2)
     .               ,max_no_face,line_face,tria_face
     .               ,quad_face  ,bvtk     ,nelemtload)
      else
        call face_vtu(ia(i_p)    ,ia(i_b1)  ,ia(i_b2)
     .               ,max_no_face,line_face,tria_face
     .               ,quad_face  ,bvtk     ,nelemtload)
      endif
c .....................................................................
c
c ... cargas 
      if(legacy) then
        call cell_data_vtk(nface,bvtk,nelemtload)
      else
        call cell_data_vtu(bvtk,nelemtload)
      endif
c
      write(aux1,'(15a)')"cargas"
      gdl =  1   
      cod =  1
      cod2 = 1
      if(legacy) then
        call point_prop_vtk(ia(i_b1),fdum,ddum,nface,aux1,ndf,gdl,cod
     .                    ,cod2    ,bvtk,nelemtload)
      else
        call point_prop_vtu(ia(i_b1),fdum,ddum,nface,aux1,ndf,gdl,cod
     .                    ,cod2    ,bvtk,nelemtload)
      endif
c .....................................................................
c
c ...
      if(legacy .eqv. .false.) then
        call cell_data_finalize_vtu(bvtk,nelemtload)
        call finalize_vtu(bvtk,nelemtload)
      endif  
c .....................................................................
c
c ...  
      malloc_name = 'faceType'
      i_b2 = dealloc(malloc_name)
      malloc_name = 'carga'
      i_b1 = dealloc(malloc_name)
      malloc_name = 'face'
      i_p  = dealloc(malloc_name)
c .....................................................................
      close(nelemtload)
c .....................................................................
      goto 101
c
c ...
  451 continue
c ... gdb graus de liberdade
c     cod  1 interio(4bytes)
c     write(aux1,'(15a)')"elmtloads"
c     gdl =  7   
c     cod =  1
c     if(legacy)then
c       call cell_prop_vtk(eload,fdum,ddum,numel,aux1,cod,gdl,bvtk
c    .                    ,nout)
c     else
c       call cell_prop_vtu(eload,fdum,ddum,numel,aux1,cod,gdl,bvtk
c    .                    ,nout)
c     endif
      goto 101
c .....................................................................      
c
c ...
  501 continue
  551 continue
  601 continue
  651 continue
  701 continue
  751 continue
  801 continue
  851 continue
  901 continue
  951 continue
 1001 continue
 1051 continue
 1101 continue
 1151 continue
 1201 continue
      goto 101 
c .....................................................................
c
c ... end      
 1251 continue
      if(legacy .eqv. .false.) then
        call cell_data_finalize_vtu(bvtk,nout)
        call finalize_vtu(bvtk,nout)
      endif  
c ..................................................................... 
      close(nout)
      return
      end
c =====================================================================
c *********************************************************************
c
c *********************************************************************
c * Data de criacao    : 09/04/2016                                   *
c * Data de modificaco : 16/10/2016                                   * 
c * ------------------------------------------------------------------*    
c * WRITE_MESH_RES_MEC: escreve a malha com os resultdados do problema*
c * mecanico no formato vtk                                           *
c * ------------------------------------------------------------------*
c * Parametros de entrada :                                           *
c * ----------------------------------------------------------------- *
c * el(nen+1,*) - conectividade com material                          *
c *  x(ndm,*)   - coordenadas                                         *
c *  u(ndf,*)   - deslocamento                                        *
c * tx(6,*)     - tensao total nodal                                  *
c * nnode       - numero de nos de vertices                           *
c * numel       - numero de elementos                                 *
c * nen         - numero de nos por elementos                         *
c * ndm         - numero de dimensoes                                 *
c * prename     - prefix do arquivo de saida                          *
c * istep       - passo tempo                                         *
c * bvtk        - true BINARY vtk false ASCII vtk                     *
c * legacy      - true (formato padrão .vtk) false (formato xlm .vtu)*
c * fprint      -                                                     *
c *                 deslocamento    (2)                               *
c *                 stress Total    (3)                               *
c * nout        - arquivo de saida                                    *
c * ----------------------------------------------------------------- *    
c * Parametros de saida :                                             *
c * ----------------------------------------------------------------- *
c * ----------------------------------------------------------------- * 
c * OBS:                                                              *
c * ----------------------------------------------------------------- *
c ********************************************************************* 
       subroutine write_mesh_res_mec(el     ,x     ,u     ,tx    
     .                             ,nnode  ,numel  
     .                             ,nen    ,ndm   ,ndf   ,ntn  
     .                             ,prename,istep
     ,                             ,bvtk   ,legacy,fprint,nout)
c ===
      use Malloc 
      implicit none
c ... variaveis da malha      
      integer nnode,numel,nen,ndm,ndf,ntn,ntn1
      integer el(nen+1,numel)
      real*8  x(ndm,*),u(ndf,*)
      real*8 tx(ntn,*)
      integer nel,nno
c ... locais     
      integer*8 i_p,i_tensor
      data i_p/1/,i_tensor/1/
      character*15 aux1
      character*30 aux
c ... tempo
      integer istep
c ... variaveis dums
      real*8 ddum
      real*4 fdum
      integer idum 
c ... arquivo      
      integer nout
      character*80 fileout,name,prename
      logical bvtk,legacy,fprint(*)
      integer cod,cod2,gdl
c =====================================================================
c
c ===
      if(legacy) then
        fileout = name(prename,istep,2)
      else  
        fileout = name(prename,istep,3)
      endif
      if(bvtk)then
        open(unit=nout,file=fileout,access='stream'
     .      ,form='unformatted',convert='big_endian')
      else
        open(unit=nout,file=fileout)
      endif  
c =====================================================================
c
c === cabecalho
      if(legacy) then
        write(aux,'(30a)')'Malha Res' 
        call head_vtk(aux,bvtk,ddum,idum,.false.,nout) 
      else  
        call head_vtu(nnode,numel,bvtk,nout) 
      endif  
c =====================================================================
c
c === Coordenadas
      if(legacy) then
        call coor_vtk(x,nnode,ndm,bvtk,nout)
      else  
        call coor_vtu(x,nnode,ndm,bvtk,nout)
      endif
c =====================================================================
c
c === Elementos
      if(legacy) then
        call elm_vtk(el,numel,nen,bvtk,nout)
      else  
        call elm_vtu(el,numel,nen,bvtk,nout)
      endif  
c =====================================================================
c
c === cell 
      if(legacy) then
        call cell_data_vtk(numel,bvtk,nout)
      else  
        call cell_data_vtu(bvtk,nout)
      endif  
c ... materiais      
      i_p = alloc_4('p       ', 1,numel)
      do nel = 1, numel
        ia(i_p+nel-1) = el(nen+1,nel)
      enddo
      write(aux1,'(15a)')'mat' 
c ... cod = 1 variaveis interias
      cod = 1
      gdl = 1
      if(legacy) then
        call cell_prop_vtk(ia(i_p),fdum,ddum,numel,aux1,cod,gdl,bvtk
     .                    ,nout)
      else
        call cell_prop_vtu(ia(i_p),fdum,ddum,numel,aux1,cod,gdl,bvtk
     .                    ,nout)
      endif 
      i_p = dealloc('p       ')
c .....................................................................
c
c ... elIdG
      i_p = alloc_4('p       ', 1,numel)
      do nel = 1, numel
        ia(i_p+nel-1) = nel
      enddo
      write(aux1,'(15a)')'elIdG' 
c ... cod = 1 variaveis interias
      cod = 1
      gdl = 1
      if(legacy) then
        call cell_prop_vtk(ia(i_p),fdum,ddum,numel,aux1,cod,gdl,bvtk
     .                    ,nout)
      else
        call cell_prop_vtu(ia(i_p),fdum,ddum,numel,aux1,cod,gdl,bvtk
     .                    ,nout)
      endif 
      i_p = dealloc('p       ')
c .....................................................................
c
c ...
      if(legacy .eqv. .false.) then
        call cell_data_finalize_vtu(bvtk,nout)
      endif 
c =====================================================================
c
c === nos  
c ...       
      if(legacy) then
        call point_data_vtk(nnode,bvtk,nout)
      else  
        call point_data_vtu(bvtk,nout)
      endif  
      i_p = alloc_4('p       ', 1,nnode)
      do nno = 1, nnode
        ia(i_p+nno-1) = nno    
      enddo
      write(aux1,'(15a)')'noIdG'
c ... gdb graus de liberdade
c     cod  1 escalar
c     cod2 1 int(4bytes) 
      gdl =  1
      cod =  1
      cod2 = 1
      if(legacy) then
        call point_prop_vtk(ia(i_p),fdum,ddum,nnode,aux1,ndm,gdl,cod
     .                    ,cod2   ,bvtk,nout)
      else
        call point_prop_vtu(ia(i_p),fdum,ddum,nnode,aux1,ndm,gdl,cod
     .                    ,cod2   ,bvtk,nout)
      endif
      i_p = dealloc('p       ')
c .....................................................................
c
c ... desloc
      if(fprint(2)) then     
        write(aux1,'(15a)')'desloc'
c ... gdb graus de liberdade
c     cod  1 escalar
c     cod2 3 real(8bytes) 
        gdl =  ndf
        cod =  2
        cod2 = 3
        if(legacy) then
          call point_prop_vtk(idum,fdum,u,nnode,aux1,ndm,gdl,cod
     .                      ,cod2,bvtk,nout)
        else
          call point_prop_vtu(idum,fdum,u,nnode,aux1,ndm,gdl,cod
     .                      ,cod2,bvtk,nout)
        endif
      endif
c .....................................................................
c
c ...
      if(fprint(3))then
c ... gerando o tensor completo
        if ( ntn .eq. 4 ) then
          i_tensor = alloc_8('tensor  ',9 ,nnode)
          ntn1     = 9
          call make_full_tensor(tx,ia(i_tensor),nnode,4, ntn1 ) 
        else if( ntn .eq. 6 ) then
          i_tensor = alloc_8('tensor  ',9 ,nnode)
          ntn1     = 9
          call make_full_tensor(tx,ia(i_tensor),nnode,6, ntn1 )
        endif
c .....................................................................
c    
c ...
        write(aux1,'(15a)')'stress'           
c ... gdb graus de liberdade
c     cod  3 tensor     
c     cod2 3 real(8bytes) 
        gdl  =  ntn1            
        cod  =  3
        cod2 =  3
        if(legacy) then
          call point_prop_vtk(idum,fdum,ia(i_tensor),nnode,aux1,ndm,gdl
     .                      ,cod ,cod2,bvtk,nout)
        else
          call point_prop_vtu(idum,fdum,ia(i_tensor),nnode,aux1,ndm,gdl
     .                       ,cod ,cod2,bvtk,nout)
        endif
c .....................................................................
c
c ...
        i_tensor = dealloc('tensor  ')
      endif
c .....................................................................
c    
c ...      
      if(legacy .eqv. .false.) then
        call point_data_finalize_vtu(bvtk,nout)
        call finalize_vtu(bvtk,nout)
      endif
c .....................................................................
c =====================================================================
      close(nout)
      return
      end
c =====================================================================
c *********************************************************************
c 
c *********************************************************************
c * MAKE_FACE: gera a conectividades das faces                        *
c * ----------------------------------------------------------------- *
c * Parametros de entrada :                                           *
c * ----------------------------------------------------------------- *
c * el         - conectividade                                        *
c * ie         - tipo do elemento                                     *
c * eload      - cargas nos elementos                                 *
c * face       - indefinido                                           *
c * carga      - indefinido                                           *
c * tipoface   - indefinido                                           *
c * numel      - numero de elementos                                  *
c * maxno      - numero maximo de no por elemento                     *
c * maxface    - numero maximo de faces                               *
c * lineface   - indefinido                                           *
c * triaface   - indefinido                                           *
c * quadface   - indefinido                                           *
c * ----------------------------------------------------------------- *
c * Parametros de saida :                                             *
c * face       - face onde as carga sao aplicadas                     *
c * carga      - carga na face                                        *
c * tipoface   - tipo da face ( 1 linha; 2 tria; 3 quad)              *
c * lineface   - numero de linha com carga                            *
c * triaface   - numero de face triangulares com carga                *
c * quadface   - numero de faca quandrigulares com carga              *
c * ----------------------------------------------------------------- *
c *********************************************************************
      subroutine make_face(el         ,ie      ,eload   ,face    ,carga
     .                    ,tipo_face  ,numel   ,maxno   ,maxface
     .                    ,max_no_face,line_face,tria_face,quad_face)
      implicit none
      integer numel,maxno
      integer i,j,k,c,maxface,nface,max_no_face,ty
      integer line_face,tria_face,quad_face
      integer carga(*),tipo_face(*)
      integer eload(maxface+1,*),face(max_no_face,*)
      integer el(maxno+1,*),ie(*)
      integer tria_side(2,3),quad_side(2,4)
      integer tetra_face(3,4),hexa_face(4,6)
c 
      data tria_side / 1, 2
     .               , 2, 3
     .               , 3, 1/
c 
      data quad_side / 1, 2
     .               , 2, 3
     .               , 3, 4
     .               , 4, 1/
c
      data tetra_face/2,3,4  
     .               ,1,4,3
     .               ,1,2,4
     .               ,1,3,2/
c
      data hexa_face/1,2,3,4  
     .              ,5,6,7,8
     .              ,1,5,6,2
     .              ,4,3,7,8
     .              ,1,4,8,5 
     .              ,2,6,7,3/
c .....................................................................
      line_face = 0 
      tria_face = 0 
      quad_face = 0 
      nface     = 0
c ...
      do i = 1, numel
        ty = ie(el(maxno+1,i))
c ... triangulo
        if( ty .eq. 2 .or. ty .eq. 3 ) then
c ... verifica se ha carga nas faces do elemento
          do j = 1, 3
            c = eload(j,i)
            if( c .ne. 0) then
              line_face         = line_face + 1
              nface             = nface + 1
              do k = 1, 2
                face(k,nface) = el(tria_side(k,j),i)
              enddo
              carga(nface)      = c
              tipo_face(nface)  = 1
            endif
          enddo
c ....................................................................
c
c ... quadrilatero
        else if(ty .eq. 4 .or. ty .eq. 5 ) then
c ... verifica se ha carga nas faces do elemento
          do j = 1, 4
            c = eload(j,i)
            if( c .ne. 0) then
              line_face         = line_face + 1
              nface             = nface + 1
              do k = 1, 2
                face(k,nface) = el(quad_side(k,j),i)
              enddo
              carga(nface)      = c
              tipo_face(nface)  = 1
            endif
          enddo
c ....................................................................
c
c ... Tetraedro
        else if( ty .eq. 6  .or. ty .eq. 16) then
c ... verifica se ha carga nas faces do elemento
          do j = 1, 4
            c = eload(j,i)
            if( c .ne. 0) then
              tria_face         = tria_face + 1
              nface             = nface + 1
              do k = 1, 3
                face(k,nface) = el(tetra_face(k,j),i)
              enddo
              carga(nface)      = c
              tipo_face(nface)  = 2
            endif
          enddo
c ....................................................................
c 
c ...  hexaedros      
        else if(ty .eq. 7 .or. ty .eq. 17) then
c ... verifica se ha carga nas faces do elemento
          do j = 1, 6
            c = eload(j,i)
            if( c .ne. 0) then
              quad_face         = quad_face + 1
              nface             = nface + 1
              do k = 1, 4
                face(k,nface) = el(hexa_face(k,j),i)
              enddo
              carga(nface)       = c
              tipo_face(nface)   = 3
            endif
          enddo
c ....................................................................
        endif
c ..................................................................... 
      enddo
c ..................................................................... 
      return
      end
c *********************************************************************
c
c *********************************************************************
c * SPLIT_U_P : pega a solução apenas nos vertices em divide em dois*
c *  vetores                                                          *
c * ----------------------------------------------------------------- *
c * Parametros de entrada :                                           *
c * ----------------------------------------------------------------- *
c * p      - indefinido                                               *
c * uv     - indefinido                                               *
c * u      - campo de pressao e deslocamente(u(4,i)=p;u(1:3)=desloc)  *
c * nnodev - numero de nos de vertices                                *
c * ndf    - numero de gaus de liberade                               *
c * ----------------------------------------------------------------- *
c * Parametros de saida :                                             *
c * ----------------------------------------------------------------- *
c * p(nnodev) - pressao                                               *
c * uv(nnodev)- deslocamento nos vertices                             *
c *********************************************************************
      subroutine split_u_p(p,uv,u,nnodev,ndf)
      implicit none
      integer nnodev,i,j,ndf
      real*8 uv(ndf-1,*),u(ndf,*),p(*)
c ...      
      do i = 1, nnodev
        do j = 1, ndf  
        enddo
      enddo
      do i = 1, nnodev
        do j = 1, ndf - 1  
          uv(j,i) = u(j,i)  
        enddo
      enddo
c ..................................................................... 
c
c ...      
      do i = 1, nnodev
        p(i) = u(ndf,i)   
      enddo
c ..................................................................... 
c
c ...
      return
      end
c *********************************************************************      
c       
c *********************************************************************
c * make_full_tensor : transforma um tensor simetrico em um tensor    *
c * geral                                                             *
c * ----------------------------------------------------------------- *
c * Parametros de entrada :                                           *
c * ----------------------------------------------------------------- *
c * tx     - tensor simetrico                                         *
c * tensor - indefinido                                               *
c * nnode  - numero de nos de vertices                                *
c * ntn    - numero total de termnos do tensor simetrico              *
c * n      - numero total de termnos do tensor geral                  *
c * ----------------------------------------------------------------- *
c * Parametros de saida :                                             *
c * ----------------------------------------------------------------- *
c * tensor - tensor geral                                             *
c * ----------------------------------------------------------------- *
c * OBS:                                                              *
c * tx    (Sxx,Syy,Szz, Sxy, Syz, Sxz)                                *
c * tensor(Sxx,Sxy,Sxz, Syx, Syy, Syz, Szx, Syz, Szz)                 *
c *********************************************************************      
      subroutine make_full_tensor(tx,tensor,nnode,ntn,n)
      implicit none
      real*8 tx(ntn,*),tensor(n,*)
      integer nnode,ntn,n,i
c ...
      if(ntn .eq. 4) then
        do i = 1, nnode
c ... sgima xx
          tensor(1,i) = tx(1,i)
c ... sgima xy
          tensor(2,i) = tx(4,i)
c ... sgima xz
          tensor(3,i) = 0.d0   
c ... sgima yx
          tensor(4,i) = tx(4,i)
c ... sgima yy
          tensor(5,i) = tx(2,i)
c ... sgima yz
          tensor(6,i) = 0.0d0  
c ... sgima zx
          tensor(7,i) = 0.0d0
c ... sgima zy
          tensor(8,i) = 0.0d0  
c ... sgima zz
          tensor(9,i) = tx(3,i)  
        enddo
c ...
      else if(ntn .eq. 6) then
        do i = 1, nnode
c ... sgima xx
          tensor(1,i) = tx(1,i)
c ... sgima xy
          tensor(2,i) = tx(4,i)
c ... sgima xz
          tensor(3,i) = tx(6,i)
c ... sgima yx
          tensor(4,i) = tx(4,i)
c ... sgima yy
          tensor(5,i) = tx(2,i)
c ... sgima yz
          tensor(6,i) = tx(5,i)
c ... sgima zx
          tensor(7,i) = tx(6,i)
c ... sgima zy
          tensor(8,i) = tx(5,i)
c ... sgima zz
          tensor(9,i) = tx(3,i)
        enddo
      endif
c .....................................................................
c
c ...
      return
      end
c **********************************************************************
c
c **********************************************************************  
c * Data de criacao    : 03/04/2016                                    *
c * Data de modificaco :                                               * 
c * ------------------------------------------------------------------ *      
c * MK_DOT_BC: Identifica o nos com valores prescritos                 *  
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ * 
c * id(ndf,*)    - numeracao da equacoes por no                        *
c * ib_bc(ndf,*) - nao definido                                        *
c * nnode        - numero de nos                                       *
c * ndf          - graus de liberdade                                  *
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ * 
c * ib_bc(ndf,*) - restricoes nodais                                   *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ *
c * id(ndf,*)    - 0 para nos restritos                                *
c * id_bc(ndf,*) - 1 para nos restritos e 0 para nos livres            *
c **********************************************************************
      subroutine mk_dot_bc(id,id_bc,nnode,ndf)
      implicit none
      integer nnode,ndf,i,j
      integer id(ndf,*),id_bc(ndf,*)

c ... identifica nos prescritos 
      do i = 1, nnode
        do j = 1, ndf
          if( id(j,i) .eq. 0) then
            id_bc(j,i)  = 1   
          endif
        enddo
      enddo
c ..................................................................... 
c
c ...
      return
      end
c *********************************************************************
