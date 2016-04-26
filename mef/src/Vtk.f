c **********************************************************************
c Depende das informacao da malha do mefpar, utiliando as variaveis do
c elementos.fi
c***********************************************************************
c
c **********************************************************************
c * Data de criacao    : 00/00/0000                                    *
c * Data de modificaco : 00/00/0000                                    * 
c * ------------------------------------------------------------------ *  
c * COOR_VTK: Escreve coordenadas no formato legacy                    *
c * -------------------------------------------------------------------*
c * Parametros de entrada:                                             *
c * -------------------------------------------------------------------*
c * coor              - coordenadas                                    *
c * nnode             - numero de nos                                  *
c * ndm               - numeros de dimensoes                           *
c * bvtk              - true formato binary false ascii                *
c * nfile             - numero associado ao arquivo de saida           *
c * -------------------------------------------------------------------*
c * Prametros de saida:                                                *
c * -------------------------------------------------------------------*
c * -------------------------------------------------------------------*
c **********************************************************************
c **********************************************************************
      subroutine coor_vtk(coor,nnode,ndm,bvtk,nfile)
      implicit none
      integer nnode,ndm,nfile
      integer i
      integer ierr
      logical bvtk,legacy
      character buffer*1024,str1*15,lf*1
      character buffer1*80,buffer2*80,buffer3*89
      Real*8 coor(ndm,*)
      Real*8 dum
      lf =char(10)
c ======================================================================
c
c ... legacy binario      
      if(bvtk)then
        write(str1(1:15),'(i15)') nnode
        buffer ='POINTS '//str1//' double'//lf
        write(nfile) trim(buffer)
c ... 
        if(ndm.eq.1)then
          dum = 0.d0
          do i = 1,nnode
            write(nfile)coor(1,i),dum,dum
          enddo
c ...
        elseif(ndm.eq.2) then
          dum = 0.d0
          do i=1,nnode
            write(nfile) coor(1,i),coor(2,i),dum
          enddo
c ...
        elseif(ndm.eq.3) then
          do i=1,nnode
            write(nfile)coor(1,i),coor(2,i),coor(3,i)
          enddo
        endif    
c ......................................................................
c
c ... legacy ascii
      else
        write(nfile,'(a,i10,a)') 'POINTS ', nnode ,' double'
c ...    
        if(ndm.eq.1) then
          do i=1,nnode
            write(nfile,'(3e15.5,3e15.5,3e15.5)') coor(1,i),0.0,0.0
          enddo
c ......................................................................
c
c ...
        elseif(ndm.eq.2) then
          do i=1,nnode
            write(nfile,'(3e15.5,3e15.5,3e15.5)') coor(1,i),coor(2,i),0.0
          enddo
c ......................................................................
c
c ...
        elseif(ndm.eq.3) then
          do i=1,nnode
            write(nfile,'(3e15.5,3e15.5,3e15.5)') coor(1,i),coor(2,i)
     .                ,coor(3,i)
          enddo
c ......................................................................
        endif
      endif  
c ======================================================================
      return
      end
c ======================================================================
c **********************************************************************
c
c **********************************************************************
c * Data de criacao    : 00/00/0000                                    *
c * Data de modificaco : 00/00/0000                                    * 
c * ------------------------------------------------------------------ *  
c * COOR_VTU: Escreve coordenadas no formato XML                       *
c * -------------------------------------------------------------------*
c * Parametros de entrada:                                             *
c * -------------------------------------------------------------------*
c * coor              - coordenadas                                    *
c * nnode             - numero de nos                                  *
c * ndm               - numeros de dimensoes                           *
c * bvtk              - true formato binary false ascii                *
c * nfile             - numero associado ao arquivo de saida           *
c * -------------------------------------------------------------------*
c * Prametros de saida:                                                *
c * -------------------------------------------------------------------*
c * -------------------------------------------------------------------*
c **********************************************************************
      subroutine coor_vtu(coor,nnode,ndm,bvtk,nfile)
      integer nnode,ndm,nfile
      integer i
      integer ierr
      logical bvtk,legacy
      character buffer*1024,str1*15,lf*1
      character buffer1*80,buffer2*80,buffer3*89
      Real*8 coor(ndm,*)
      Real*8 dum
      lf =char(10)
c ======================================================================
c ... xml binario      
      if(bvtk)then
        print*,'não implementado.'
        stop
c ... xml ascii
      else
        write(nfile,'(a)') '<Points>'
        buffer1='<DataArray type="Float64" '
        buffer2=' NumberOfComponents ="3" ' 
        buffer3=' format="ascii">' 
        buffer = trim(buffer1) // trim(buffer2) // trim(buffer3)
        write(nfile,'(a)') trim(buffer)
c ...          
        if(ndm.eq.1) then
          do i=1,nnode
            write(nfile,'(3e15.5,3e15.5,3e15.5)') coor(1,i),0.0,0.0
          enddo
c ......................................................................
c
c ...
        elseif(ndm.eq.2) then
          do i=1,nnode
            write(nfile,'(3e15.5,3e15.5,3e15.5)') coor(1,i),coor(2,i),0.0
          enddo
c ......................................................................
c
c ...
        elseif(ndm.eq.3) then
          do i=1,nnode
            write(nfile,'(3e15.5,3e15.5,3e15.5)') coor(1,i),coor(2,i)
     .                 ,coor(3,i)
          enddo
        endif  
c ......................................................................
        write(nfile,'(a)') '</DataArray>'
        write(nfile,'(a)') '</Points>'
      endif  
c ======================================================================
      return
      end
c ======================================================================
c **********************************************************************
c
c **********************************************************************
      subroutine elm_vtk(nos,numel,nen,bvtk,nfile)
c **********************************************************************
c * Data de criacao    : 00/00/0000                                    *
c * Data de modificaco : 00/00/0000                                    * 
c * ------------------------------------------------------------------ *  
c * ELM_VTK: escreve elementos nos formato vtk com os seus respectivos *
c *           materias                                                 *
c * -------------------------------------------------------------------*
c * Parametros de entrada:                                             *
c * -------------------------------------------------------------------*
c * nos(nen+1,numel)  - conetividade nodai                             *
c * numel             - numero de elementos                            *
c * ndm               - numeros de dimensoes                           *
c * nen               - numero maximo de nos por elemento              *
c * nfile             - arquivo de saida                               *
c * nbar2             - numeros de barras                              *
c * ntria3            - numero de triangulos                           *
c * nquad4            - numero de quaddrilateros                       *
c * ntetra4           - numero de tetraedros                           *
c * nhexa8            - numero de hexaedros                            *
c * bvtk              - true formato binary false ascii                *
c * nfile             - numero associado ao arquivo de saida           *
c * -------------------------------------------------------------------*
c * Prametros de saida:                                                *
c * -------------------------------------------------------------------*
c * -------------------------------------------------------------------*
c * OBS:                                                               *
c * -------------------------------------------------------------------*
c * Os prametors nbar2,ntria3,nquad4,ntetra4 e nhexa8 foram herdados do*
c *   mefpar                                                           *
c *--------------------------------------------------------------------*
c **********************************************************************
      implicit none
      include 'elementos.fi'
      integer numel,nen,nfile
      integer nnoel
      integer numet,nb2,nt3,nq4,nt4,nh8,nt10,nh20
      integer i,j
      integer nos(nen+1,*)
      character buffer*1024,lf*1,str1*15,str2*15
      logical bvtk

      lf =char(10)
c ======================================================================
c
c ... Calculo do numero de tipo de cada elemento
      nb2  = nbar2(1)
      nt3  = ntria3(1) 
      nq4  = nquad4(1) 
      nt4  = ntetra4(1)
      nh8  = nhexa8(1)
      nt10 = ntetra10(1)
      nh20 = nhexa20(1)
      numet =3*nb2 + 4*nt3 + 5*nq4 + 5*nt4 + 9*nh8 + 11*nt10
     .      +21*nh20
c ... elemento em overllaping 
      nt3  = ntria3(3) 
      nq4  = nquad4(3) 
      nt4  = ntetra4(3)
      nh8  = nhexa8(3)
      nt10 = ntetra10(3)
      nh20 = nhexa20(3)
      numet = 4*nt3 + 5*nq4 + 5*nt4 + 9*nh8 + 11*nt10 + 21*nh20 + numet
c ... total de elemntos e tamnhanho dos dados totais
      if(bvtk)then
        write(str1(1:15),'(i15)')numel
        write(str2(1:15),'(i15)')numet
        buffer = lf//'CELLS '//str1//str2//lf
        write(nfile) trim(buffer)  
      else
        write(nfile,'(a,i10,i10)') 'CELLS ',numel,numet  
      endif 
c ... escrevendo a malha
c     
c ... nos dos elementos
c      
c ......................................................................
      if (nbar2(1) .gt. 0) then
        nnoel = 2 
        do i = nbar2(2), nbar2(2) + nbar2(1)-1
          if(bvtk)then
             write(nfile) nnoel,(nos(j,i)-1,j=1,2)
          else
             write(nfile,'(10i10)') nnoel,(nos(j,i)-1,j=1,2)
          endif  
        enddo
      endif
c ......................................................................
c
c ......................................................................
      if (ntria3(1) .gt. 0) then
        nnoel = 3 
        do i = ntria3(2), ntria3(2) + ntria3(1)-1
          if(bvtk)then
            write(nfile) nnoel,(nos(j,i)-1,j=1,3)
          else
            write(nfile,'(10i10)') nnoel,(nos(j,i)-1,j=1,3)
          endif  
        enddo
      endif
c ......................................................................
c
c ......................................................................
      if (nquad4(1) .gt. 0) then
        nnoel = 4 
        do i = nquad4(2), nquad4(2) + nquad4(1)-1
          if(bvtk)then
            write(nfile) nnoel,(nos(j,i)-1,j=1,4)
          else
            write(nfile,'(10i10)') nnoel,(nos(j,i)-1,j=1,4)
          endif  
        enddo
      endif
c ......................................................................
      if (ntetra4(1) .gt. 0) then
        nnoel = 4 
        do i = ntetra4(2), ntetra4(2) + ntetra4(1)-1
          if(bvtk)then
            write(nfile) nnoel,(nos(j,i)-1,j=1,4)
          else
            write(nfile,'(10i10)') nnoel,(nos(j,i)-1,j=1,4)
          endif  
        enddo              
      endif
c ......................................................................
      if (nhexa8(1) .gt. 0) then
        nnoel = 8 
        do i = nhexa8(2), nhexa8(2) + nhexa8(1)-1
          if(bvtk)then
            write(nfile) nnoel,(nos(j,i)-1,j=1,8)
          else
            write(nfile,'(10i10)') nnoel,(nos(j,i)-1,j=1,8)
          endif  
        enddo 
      endif
c ......................................................................
c
c ......................................................................
      if (ntetra10(1) .gt. 0) then
        nnoel = 10 
        do i = ntetra10(2), ntetra10(2) + ntetra10(1)-1
          if(bvtk)then
            write(nfile) nnoel,(nos(j,i)-1,j=1,10)
          else
            write(nfile,'(11i10)') nnoel,nos(1 ,i)-1,nos(2,i)-1
     .                                  ,nos(3 ,i)-1,nos(4,i)-1
     .                                  ,nos(5 ,i)-1,nos(8,i)-1
     .                                  ,nos(6 ,i)-1,nos(7,i)-1
     .                                  ,nos(10,i)-1,nos(9,i)-1
          endif  
        enddo 
      endif
c ......................................................................
c
c ......................................................................
      if (nhexa20(1) .gt. 0) then
        nnoel = 20 
        do i = nhexa20(2), nhexa20(2) + nhexa20(1)-1
          if(bvtk)then
            write(nfile) nnoel,(nos(j,i)-1,j=1,20)
          else
            write(nfile,'(21i10)') nnoel
     .                  ,nos( 5,i)-1,nos( 6,i)-1,nos( 7,i)-1,nos( 8,i)-1
     .                  ,nos( 1,i)-1,nos( 2,i)-1,nos( 3,i)-1,nos( 4,i)-1
     .                  ,nos(13,i)-1,nos(14,i)-1,nos(15,i)-1,nos(16,i)-1
     .                  ,nos( 9,i)-1,nos(10,i)-1,nos(11,i)-1,nos(12,i)-1
     .                  ,nos(17,i)-1,nos(18,i)-1,nos(19,i)-1,nos(20,i)-1
          endif  
        enddo 
      endif
c ......................................................................
c
c ... elementos em overllaping
       if (ntria3(3) .gt. 0) then
        nnoel = 3 
        do i = ntria3(4), ntria3(4) + ntria3(3)-1
          if(bvtk)then
            write(nfile) nnoel,(nos(j,i)-1,j=1,3)
          else
            write(nfile,'(10i10)') nnoel,(nos(j,i)-1,j=1,3)
          endif  
        enddo
      endif
c ......................................................................
c
c ......................................................................
      if (nquad4(3) .gt. 0) then
        nnoel = 4 
        do i = nquad4(4), nquad4(4) + nquad4(3)-1
          if(bvtk)then
            write(nfile) nnoel,(nos(j,i)-1,j=1,4)
          else
            write(nfile,'(10i10)') nnoel,(nos(j,i)-1,j=1,4)
          endif  
        enddo
      endif
c ......................................................................
      if (ntetra4(3) .gt. 0) then
        nnoel = 4 
        do i = ntetra4(4), ntetra4(4) + ntetra4(3)-1
          if(bvtk)then
            write(nfile) nnoel,(nos(j,i)-1,j=1,4)
          else
            write(nfile,'(10i10)') nnoel,(nos(j,i)-1,j=1,4)
          endif  
        enddo              
      endif
c ......................................................................
      if (nhexa8(3) .gt. 0) then
        nnoel = 8 
        do i = nhexa8(4), nhexa8(4) + nhexa8(3)-1
          if(bvtk)then
            write(nfile) nnoel,(nos(j,i)-1,j=1,8)
          else
            write(nfile,'(10i10)') nnoel,(nos(j,i)-1,j=1,8)
          endif  
        enddo 
      endif
c ......................................................................
c ======================================================================
c ... tipo dos elementos
c
      if(bvtk)then
        write(str1(1:15),'(i15)')numel
        buffer = lf//'CELL_TYPES '//str1//lf
        write(nfile) trim(buffer)
      else
        write(nfile,'(a,i10)') 'CELL_TYPES ',numel
      endif
c .......................................................................      
c
c ...     
      if (nbar2(1) .gt. 0) then
        nnoel = 3
        do i = nbar2(2), nbar2(2) + nbar2(1)-1
          if(bvtk)then
            write(nfile)nnoel
          else  
            write(nfile,'(i3)')nnoel
          endif  
        enddo
      endif
c .......................................................................
c
c ...
      if (ntria3(1) .gt. 0) then
        nnoel = 5
        do i = ntria3(2), ntria3(2) + ntria3(1)-1
          if(bvtk)then
            write(nfile)nnoel
          else  
            write(nfile,'(i3)')nnoel
          endif  
        enddo
      endif
c ......................................................................
c
c ...
      if (nquad4(1) .gt. 0) then
        nnoel = 9
        do i = nquad4(2), nquad4(2) + nquad4(1)-1
          if(bvtk)then
            write(nfile)nnoel
          else  
            write(nfile,'(i3)')nnoel
          endif  
        enddo
      endif
c ......................................................................
c
c ...
      if (ntetra4(1) .gt. 0) then
        nnoel = 10
        do i = ntetra4(2), ntetra4(2) + ntetra4(1)-1
          if(bvtk)then
            write(nfile)nnoel
          else  
            write(nfile,'(i3)') nnoel
          endif  
        enddo
      endif
c ......................................................................
c
c ...
      if (nhexa8(1) .gt. 0) then
        nnoel = 12
        do i = nhexa8(2), nhexa8(2) + nhexa8(1)-1
          if(bvtk)then
            write(nfile)nnoel
          else  
            write(nfile,'(i3)')nnoel
          endif  
        enddo
      endif
c ......................................................................
c
c ...
      if (ntetra10(1) .gt. 0) then
        nnoel = 24
        do i = ntetra10(2), ntetra10(2) + ntetra10(1)-1
          if(bvtk)then
            write(nfile)nnoel
          else  
            write(nfile,'(i3)')nnoel
          endif  
        enddo
      endif
c ......................................................................
c
c ...
      if (nhexa20(1) .gt. 0) then
        nnoel = 25
        do i = nhexa20(2), nhexa20(2) + nhexa20(1)-1
          if(bvtk)then
            write(nfile)nnoel
          else  
            write(nfile,'(i3)')nnoel
          endif  
        enddo
      endif
c ......................................................................
c
c ... elmentos em overllaping    
      if (ntria3(3) .gt. 0) then
        nnoel = 5
        do i = ntria3(4), ntria3(4) + ntria3(3)-1
          if(bvtk)then
            write(nfile)nnoel
          else  
            write(nfile,'(i3)')nnoel
          endif  
        enddo
      endif
c ......................................................................
c
c ...
      if (nquad4(3) .gt. 0) then
        nnoel = 9
        do i = nquad4(4), nquad4(4) + nquad4(3)-1
          if(bvtk)then
            write(nfile)nnoel
          else  
            write(nfile,'(i3)')nnoel
          endif  
        enddo
      endif
c ......................................................................
c
c ...
      if (ntetra4(3) .gt. 0) then
        nnoel = 10
        do i = ntetra4(4), ntetra4(4) + ntetra4(3)-1
          if(bvtk)then
            write(nfile)nnoel
          else  
            write(nfile,'(i3)') nnoel
          endif  
        enddo
      endif
c ......................................................................
c
c ...
      if (nhexa8(3) .gt. 0) then
        nnoel = 12
        do i = nhexa8(4), nhexa8(4) + nhexa8(3)-1
          if(bvtk)then
            write(nfile)nnoel
          else  
            write(nfile,'(i3)')nnoel
          endif  
        enddo
      endif
c ......................................................................
c
      return
      end
c ======================================================================
c **********************************************************************
c
c **********************************************************************
      subroutine elm_vtu(nos,numel,nen,bvtk,nfile)
c **********************************************************************
c * Data de criacao    : 00/00/0000                                    *
c * Data de modificaco : 00/00/0000                                    * 
c * ------------------------------------------------------------------ *  
c * ELM_VTU: escreve elementos nos formato vtu (xml)                   *
c * -------------------------------------------------------------------*
c * Parametros de entrada:                                             *
c * -------------------------------------------------------------------*
c * nos(nen+1,numel)  - conetividade nodai                             *
c * numel             - numero de elementos                            *
c * ndm               - numeros de dimensoes                           *
c * nen               - numero maximo de nos por elemento              *
c * nfile             - arquivo de saida                               *
c * nbar2             - numeros de barras                              *
c * ntria3            - numero de triangulos                           *
c * nquad4            - numero de quaddrilateros                       *
c * ntetra4           - numero de tetraedros                           *
c * nhexa8            - numero de hexaedros                            *
c * bvtk              - true formato binary false ascii                *
c * nfile             - numero associado ao arquivo de saida           *
c * -------------------------------------------------------------------*
c * Prametros de saida:                                                *
c * -------------------------------------------------------------------*
c * -------------------------------------------------------------------*
c * OBS:                                                               *
c * -------------------------------------------------------------------*
c * Os prametors nbar2,ntria3,nquad4,ntetra4 e nhexa8 foram herdados do*
c *   mefpar                                                           *
c *--------------------------------------------------------------------*
c **********************************************************************
      implicit none
      include 'elementos.fi'
      integer numel,nen,nfile
      integer nnoel
      integer numet,nb2,nt3,nq4,nt4,nh8
      integer i,j
      integer nos(nen+1,*)
      integer*8 offset
      character buffer*1024,lf*1,str1*15,str2*15,buffer1*80,buffer2*80
      logical bvtk

      lf =char(10)
c ======================================================================
c
c ... Calculo do numero de tipo de cada elemento 
      nb2 = nbar2(1)
      nt3 = ntria3(1) 
      nq4 = nquad4(1) 
      nt4 = ntetra4(1)
      nh8 = nhexa8(1)
      numet = 3*nb2 + 4*nt3 + 5*nq4 + 5*nt4 + 9*nh8 
c
c ... total de elemntos e tamnhanho dos dados totais
      if(bvtk)then
        print*,'nao implementado.'
        stop
      else
        write(nfile,'(a)') '<Cells> '
        buffer1 ='<DataArray type="Int32" Name="connectivity" '
        buffer2 =' format="ascii">'
        buffer = trim(buffer1) // trim(buffer2)
        write(nfile,'(a)') trim(buffer)  
      endif 
c ... escrevendo a malha
c     
c ... nos dos elementos
c
c ......................................................................
      
      if (nbar2(1) .gt. 0) then
         do i = nbar2(2), nbar2(2) + nbar2(1)-1
           if(bvtk)then
           else
              write(nfile,'(10i10)') (nos(j,i)-1,j=1,2)
           endif  
         enddo
       endif
c ......................................................................
c
c ......................................................................
      if (ntria3(1) .gt. 0) then
        do i = ntria3(2), ntria3(2) + ntria3(1)-1
          if(bvtk)then
          else
            write(nfile,'(10i10)')(nos(j,i)-1,j=1,3)
          endif  
        enddo
      endif
c ......................................................................
c
c ......................................................................
      if (nquad4(1) .gt. 0) then
        do i = nquad4(2), nquad4(2) + nquad4(1)-1
          if(bvtk)then
          else
            write(nfile,'(10i10)')(nos(j,i)-1,j=1,4)
          endif  
        enddo
      endif
c ......................................................................
      if (ntetra4(1) .gt. 0) then
        do i = ntetra4(2), ntetra4(2) + ntetra4(1)-1
          if(bvtk)then
          else
            write(nfile,'(10i10)')(nos(j,i)-1,j=1,4)
          endif  
        enddo              
      endif
c ......................................................................
      if (nhexa8(1) .gt. 0) then
        do i = nhexa8(2), nhexa8(2) + nhexa8(1)-1
          if(bvtk)then
          else
            write(nfile,'(10i10)') (nos(j,i)-1,j=1,8)
          endif  
        enddo 
      endif
c ......................................................................
      if(bvtk)then
      else
        write(nfile,'(a)') '</DataArray>'
      endif
c ... offset no arranjo das conectvidade(funcionando para malha
c                         de um elemento apenas)
      if(bvtk)then
      else
        buffer1 ='<DataArray type="Int64" Name="offsets" '
        buffer2 =' format="ascii">'
        buffer = trim(buffer1) // trim(buffer2)
        write(nfile,'(a)') trim(buffer)
        if (nbar2(1) .gt. 0) then
          offset = 2
          do i = 1, numel
            write(nfile,'(i20)') offset
              offset = 2 + offset
            enddo  
        endif  
        if (ntria3(1) .gt. 0) then
          offset = 3
          do i = 1, numel
            write(nfile,'(i20)') offset
            offset = 3 + offset
          enddo  
        endif  
        if (nquad4(1) .gt. 0) then
          offset = 4
          do i = 1, numel
            write(nfile,'(i20)') offset
            offset = 4 + offset
          enddo  
        endif  
        if (ntetra4(1) .gt. 0) then
          offset = 4
          do i = 1, numel
            write(nfile,'(i20)') offset
            offset = 4 + offset
          enddo  
        endif  
        if (nhexa8(1) .gt. 0) then
          offset = 8
          do i = 1, numel
            write(nfile,'(i20)') offset
            offset = 8 + offset
          enddo  
        endif  
        write(nfile,'(a)') '</DataArray>'
      endif 
c ======================================================================
c ... tipo dos elementos
c
      if(bvtk)then
      else
        buffer1 ='<DataArray type="UInt8" Name="types" '
        buffer2 =' format="ascii">'
        buffer = trim(buffer1) // trim(buffer2)
        write(nfile,'(a)') trim(buffer)
      endif  
c
c ...     
      if (nbar2(1) .gt. 0) then
        nnoel = 3
        do i = nbar2(2), nbar2(2) + nbar2(1)-1
          if(bvtk)then
            write(nfile)nnoel
          else  
            write(nfile,'(i3)')nnoel
          endif  
        enddo
      endif
c .......................................................................
c
c ...
      if (ntria3(1) .gt. 0) then
        nnoel = 5
        do i = ntria3(2), ntria3(2) + ntria3(1)-1
          if(bvtk)then
            write(nfile)nnoel
          else  
            write(nfile,'(i3)')nnoel
          endif  
        enddo
      endif
c ......................................................................
c
c ...
      if (nquad4(1) .gt. 0) then
        nnoel = 9
        do i = nquad4(2), nquad4(2) + nquad4(1)-1
          if(bvtk)then
            write(nfile)nnoel
          else  
            write(nfile,'(i3)')nnoel
          endif  
        enddo
      endif
c ......................................................................
c
c ...
      if (ntetra4(1) .gt. 0) then
        nnoel = 10
        do i = ntetra4(2), ntetra4(2) + ntetra4(1)-1
          if(bvtk)then
            write(nfile)nnoel
          else  
            write(nfile,'(i3)') nnoel
          endif  
        enddo
      endif
c ......................................................................
c
c ...
      if (nhexa8(1) .gt. 0) then
        nnoel = 12
        do i = nhexa8(2), nhexa8(2) + nhexa8(1)-1
          if(bvtk)then
            write(nfile)nnoel
          else  
            write(nfile,'(i3)')nnoel
          endif  
        enddo
      endif
c ......................................................................
c
      if(bvtk)then
      else
        write(nfile,'(a)') '</DataArray>'
        write(nfile,'(a)') '</Cells>'
      endif  
c ......................................................................
c
c ====
      return
      end
c ======================================================================
c **********************************************************************
c
c **********************************************************************
c * Data de criacao    : 00/00/0000                                    *
c * Data de modificaco : 00/00/0000                                    * 
c * ------------------------------------------------------------------ *  
c * PONT_PROP_VTK : Escreve propriedades nodais                        *
c * -------------------------------------------------------------------*
c * Parametros de entrada:                                             *
c * -------------------------------------------------------------------*
c * res(ndf,nnode)    - campo escalar,vetorial ou tensorial            *
c * nnode             - numero de nos                                  *
c * ndm               - dimesao                                        *
c * gdl               - graus de liberdade                             *
c * nfile             - arquivo de saida                               *
c * cname             - nome da variavel                               *
c * cod1              - codigo de instrucao                            *
c *                   1 -> campo escalar                               *
c *                   2 -> campo vetorial                              *
c *                   3 -> campo tensorial                             *
c * cod2              - 1 - interio de  4 bytes                        *  
c *                   - 2 - real de 4 bytes                            *  
c *                   - 3 - real de 8 bytes                            *  
c * bvtk              - true BINARY vtk, false ASCII                   *
c * nfile - numero associado ao arquivo de saida                       *
c * -------------------------------------------------------------------*
c * Parmetros de saida:                                                *
c * -------------------------------------------------------------------*
c **********************************************************************
      subroutine point_prop_vtk(iprop,fprop,dprop,nnode,cname,ndm,gdl
     .                        ,cod1 
     .                        ,cod2,bvtk,nfile)
      implicit none
      integer nnode,ndm,gdl,cod1,cod2,nfile
      integer  i,j
      integer iprop(gdl,*)
      Real*4 fprop(gdl,*)
      Real*8 dprop(gdl,*)
      character buffer*1024,lf*1,cname*15
      integer idum
      real fdum
      real*8 ddum
      logical bvtk
      lf = char(10)
c ======================================================================
c
c === BINARY
      if(bvtk)then
c ... campo escalar
        if(cod1.eq.1) then
c .. escalar int BINARY     
          if(cod2.eq.1)then
            buffer =' SCALARS '//cname//' int'//lf
            write(nfile)trim(buffer)
            buffer =' LOOKUP_TABLE '//'default '//lf
            write(nfile)trim(buffer)
            do i=1,nnode
              write(nfile)iprop(1,i)
            enddo
c ... escalar float BINARY         
          elseif(cod2.eq.2)then  
            buffer = lf//'SCALARS '//cname//' float'//lf
            write(nfile)trim(buffer)
            buffer = lf//'LOOKUP_TABLE '//'default '//lf
            write(nfile)trim(buffer)
            do i=1,nnode
              write(nfile)fprop(1,i)
            enddo
c ... escalardouble BINARY         
          elseif(cod2.eq.3)then  
            buffer = lf//'SCALARS '//cname//' double'//lf
            write(nfile)trim(buffer)
            buffer = lf//'LOOKUP_TABLE '//'default '//lf
            write(nfile)trim(buffer)
            do i=1,nnode
              write(nfile)dprop(1,i)
            enddo
          endif  
c ......................................................................
c
c ... campo vetorial BINARY
        elseif(cod1.eq.2) then
          if(cod2.eq.1)then
            buffer =' VECTORS '//cname//' int'//lf
            write(nfile)trim(buffer)
            idum = 0
            do i=1,nnode
              if(ndm .eq. 1) then
                write(nfile)iprop(1,i),idum,idum
              else if(ndm .eq. 2) then
                write(nfile)iprop(1,i),iprop(2,i),idum
              elseif(ndm .eq. 3) then
                write(nfile)iprop(1,i),iprop(2,i),iprop(3,i)
              endif
            enddo
c ... escalar float BINARY         
          elseif(cod2.eq.2)then  
            buffer =' VECTORS '//cname//' float'//lf
            write(nfile)trim(buffer)
            fdum = 0.d0
            do i=1,nnode
              if(ndm .eq. 2) then
                write(nfile) fprop(1,i),fprop(2,i),fdum
              elseif(ndm .eq. 3) then
                write(nfile)fprop(1,i),fprop(2,i),fprop(3,i)
              endif
            enddo
c ... escalar double BINARY         
          elseif(cod2.eq.3)then  
            buffer =' VECTORS '//cname//' double'//lf
            write(nfile)trim(buffer)
            ddum = 0.d0
            do i=1,nnode
              if(ndm .eq. 1) then
                write(nfile) dprop(1,i),ddum,ddum
              else if(ndm .eq. 2) then
                write(nfile) dprop(1,i),dprop(2,i),ddum
              elseif(ndm .eq. 3) then
                write(nfile)dprop(1,i),dprop(2,i),dprop(3,i)
              endif
            enddo
          endif  
c ......................................................................
c 
c ... campo tensorial BINARY
        elseif(cod1.eq.3) then
          print*,'\nnao implementado '
        endif
c ......................................................................
c
c ======================================================================
c
c === ASCII
      else
c ... campo escalar ASCII
        if(cod1.eq.1) then
c .. escalar int ASCII 
          if(cod2.eq.1)then
            write(nfile,'(a,1x,a15,1x,a8,1x,i3)')'SCALARS'
     .           ,cname,'int    ',gdl
            write(nfile,'(a)')'LOOKUP_TABLE default'
            do i = 1, nnode
              write(nfile,'(99i10)')(iprop(j,i),j=1,gdl) 
            enddo
c .. escalar float ASCII 
          elseif(cod2.eq.2)then  
           write(nfile,'(a,1x,a15,1x,a8,1x,i3)')'SCALARS'
     .           ,cname,'float  ',gdl
            write(nfile,'(a)')'LOOKUP_TABLE default'
            do i = 1, nnode
               write(nfile,'(99e15.5e3)')(fprop(j,i),j=1,gdl) 
            enddo
c .. escalar double ASCII 
          elseif(cod2.eq.3)then  
            write(nfile,'(a,1x,a15,1x,a8,1x,i3)')'SCALARS'
     .           ,cname,'double ',gdl
            write(nfile,'(a)')'LOOKUP_TABLE default'
            do i = 1, nnode
               write(nfile,'(99e15.5e3)')(dprop(j,i),j=1,gdl) 
            enddo
          endif  
c ......................................................................
c
c ... campo vetorial ASCII
        elseif(cod1.eq.2) then
c ... int ASCII
          if(cod2.eq.1)then
            idum = 0
            write(nfile,'(a,15a,a)')'VECTORS ',cname,' int'
            do i=1,nnode
              if(ndm .eq. 1) then
                write(nfile,'(i10,i10,i10)') iprop(1,i),idum,idum
              elseif(ndm .eq. 2) then
                write(nfile,'(i10,i10,i10)') iprop(1,i),iprop(2,i),idum
              elseif(ndm .eq. 3) then
                write(nfile,'(i10,i10,i10)')
     .          iprop(1,i),iprop(2,i),iprop(3,i)
              endif
            enddo
c .. float ASCII 
          elseif(cod2.eq.2)then  
            fdum = 0.d0
            write(nfile,'(a,15a,a)')'VECTORS ',cname,' float'
            do i=1,nnode
              if(ndm .eq. 2) then
                write(nfile,'(i10,i10,i10)') fprop(1,i),fprop(2,i),fdum
              elseif(ndm .eq. 3) then
                write(nfile,'(i10,i10,i10)')
     .          fprop(1,i),fprop(2,i),fprop(3,i)
              endif
            enddo
c .. double ASCII 
          elseif(cod2.eq.3)then
            ddum = 0.d0
            write(nfile,'(a,15a,a)')'VECTORS ',cname,' double'
            do i=1,nnode
              if(ndm .eq. 1) then
                write(nfile,'(3e15.5e3)') dprop(1,i),ddum,ddum
              else if(ndm .eq. 2) then
                write(nfile,'(3e15.5e3)') dprop(1,i),dprop(2,i),ddum
              elseif(ndm .eq. 3) then
                write(nfile,'(3e15.5e3)')
     .          dprop(1,i),dprop(2,i),dprop(3,i)
              endif
            enddo
          endif  
c ......................................................................
c 
c ... campo tensorial ASCII
        elseif(cod1.eq.3) then
c ... integer ASCII
          if(cod2.eq.1)then
            print*,'\nnao implementado '
c ... float ASCII 
          elseif(cod2.eq.2)then  
            print*,'\nnao implementado '
c ... double ASCII 
          elseif(cod2.eq.3)then
            write(nfile,'(a,15a,a)')'TENSORS ',cname,' double'
            do i=1,nnode
              write(nfile,'(9e15.5e3)')dprop(1,i),dprop(2,i),dprop(3,i)
     .                                ,dprop(4,i),dprop(5,i),dprop(6,i)
     .                                ,dprop(7,i),dprop(8,i),dprop(9,i)
            enddo
c ......................................................................
          endif
c ......................................................................
        endif
c ......................................................................
      endif  
c ======================================================================
c
      return
      end
c ======================================================================
c **********************************************************************
c
c **********************************************************************
c * Data de criacao    : 00/00/0000                                    *
c * Data de modificaco : 00/00/0000                                    * 
c * ------------------------------------------------------------------ *  
c * PONT_PROP_VTU : Escreve propriedades nodais (xlm)                  *
c * -------------------------------------------------------------------*
c * Parametros de entrada:                                             *
c * -------------------------------------------------------------------*
c * res(ndf,nnode)    - campo escalar,vetorial ou tensorial            *
c * nnode             - numero de nos                                  *
c * ndm               - dimesao                                        *
c * gdl               - graus de liberdade                             *
c * nfile             - arquivo de saida                               *
c * cname             - nome da variavel                               *
c * cod1              - codigo de instrucao                            *
c *                   1 -> campo escalar                               *
c *                   2 -> campo vetorial                              *
c *                   3 -> campo tensorial                             *
c * cod2              - 1 - interio de  4 bytes                        *  
c *                   - 2 - real de 4 bytes                            *  
c *                   - 3 - real de 8 bytes                            *  
c * bvtk              - true BINARY vtk, false ASCII                   *
c * nfile - numero associado ao arquivo de saida                       *
c * -------------------------------------------------------------------*
c * Parmetros de saida:                                                *
c * -------------------------------------------------------------------*
c **********************************************************************
      subroutine point_prop_vtu(iprop,fprop,dprop,nnode,cname,ndm,gdl
     .                        ,cod1,cod2,bvtk,nfile)
      implicit none
      integer nnode,ndm,gdl,cod1,cod2,nfile
      integer  i,j
      integer iprop(gdl,*)
      Real*4 fprop(gdl,*)
      Real*8 dprop(gdl,*)
      character buffer*1024,cname*15,buffer1*80,buffer2*80,str1*15
      character buffer3*80
      integer idum
      real fdum
      real*8 ddum
      logical bvtk
c ======================================================================
c
c === BINARY
      if(bvtk)then
c ... campo escalar
        if(cod1.eq.1) then
c .. escalar int BINARY     
          if(cod2.eq.1)then
c ... escalar float BINARY         
          elseif(cod2.eq.2)then  
c ... escalardouble BINARY         
          elseif(cod2.eq.3)then  
          endif  
c ......................................................................
c
c ... campo vetorial BINARY
        elseif(cod1.eq.2) then
          if(cod2.eq.1)then
c ... escalar float BINARY         
          elseif(cod2.eq.2)then  
c ... escalardouble BINARY         
          elseif(cod2.eq.3)then  
          endif  
c ......................................................................
c 
c ... campo tensorial BINARY
        elseif(cod1.eq.3) then
          print*,'\nnao implementado '
        endif
c ......................................................................
c
c ======================================================================
c
c === ASCII
      else
c ... campo escalar ASCII
        if(cod1.eq.1) then
c ... escalar int ASCII 
          if(cod2.eq.1)then
            buffer1= '<DataArray type="Int32"'
            buffer2= ' Name="' // trim( cname) // '" format="ascii">'
            buffer=trim(buffer1)//trim(buffer2)
            write(nfile,'(a)')trim(buffer)
            do i=1,nnode
              write(nfile,'(i10)') iprop(1,i)
            enddo
c ... escalar float ASCII 
          elseif(cod2.eq.2)then  
c ... escalar double ASCII 
          elseif(cod2.eq.3)then  
            buffer1= '<DataArray type="Float64" '
            buffer2= ' Name="' // trim( cname) // '" format="ascii">'
            buffer=trim(buffer1)//trim(buffer2)
            write(nfile,'(a)')trim(buffer)
            do i=1,nnode
              write(nfile,'(e15.5e3)') dprop(1,i)
            enddo
          endif  
c ......................................................................
c
c ... campo vetorial ASCII
        elseif(cod1.eq.2) then
c ... vetorial int ASCII 
          if(cod2.eq.1)then
            write(str1,'(i1)')gdl
            buffer1= '<DataArray type="Int32" '
            buffer2= ' NumberOfComponents="'//trim(str1)//'"'
            buffer3= ' Name="' // trim(cname) // '" format="ascii">'
            buffer=trim(buffer1)//trim(buffer2)//trim(buffer3)
            write(nfile,'(a)')trim(buffer)
            do i=1,nnode
              write(nfile,'(10i10)') (iprop(j,i),j=1,gdl)
            enddo
c ... vetorial float ASCII 
          elseif(cod2.eq.2)then  
c ... vetorial double ASCII 
          elseif(cod2.eq.3)then
            write(str1,'(i1)')gdl
            buffer1= '<DataArray type="Float64" '
            buffer2= ' NumberOfComponents="'//trim(str1)//'"'
            buffer3= ' Name="' // trim(cname) // '" format="ascii">'
            buffer=trim(buffer1)//trim(buffer2)//trim(buffer3)
            write(nfile,'(a)')trim(buffer)
            do i=1,nnode
              write(nfile,'(10e15.5e3)')(dprop(j,i),j=1,gdl)
            enddo
          endif  
c ......................................................................
c ......................................................................
c 
c ... campo tensorial ASCII
        elseif(cod1.eq.3) then
          print*,'nao implementado '
c ......................................................................
        endif
        write(nfile,'(a)')'</DataArray>'
      endif  
c ======================================================================
c
      return
      end
c ======================================================================
c **********************************************************************
c
c **********************************************************************
c * Data de criacao    : 00/00/0000                                    *
c * Data de modificaco : 00/00/0000                                    * 
c * ------------------------------------------------------------------ *  
c * CELL_PROP_VTK : escreve propriedades das celulas                   *
c * -------------------------------------------------------------------*
c * Parametros de entrada:                                             *
c * -------------------------------------------------------------------*
c * dprop - propriedas real 8 bytes                                    *
c * iprop - propriedas interia de 8 bytes                              *
c * fprop - propriedas real de 4 bytes                                 *
c * cod   - 1 - interio de  4 bytes                                    *  
c *       - 2 - real de 4 bytes                                        *  
c *       - 3 - real de 8 bytes                                        *  
c * gdl   - graus de liberdade                                         *  
c * bvtk  - escrita no formato vtk binario                             *  
c * nfile - numero associado ao arquivo de saida                       *
c * -------------------------------------------------------------------*
c * Parmetros de saida:                                                *
c * -------------------------------------------------------------------*
c **********************************************************************
      subroutine cell_prop_vtk(iprop,fprop,dprop,numel,cname,cod,gdl 
     .                        ,bvtk,nfile)
c ===
      implicit none
      integer gdl
      integer iprop(gdl,*)
      real*4 fprop(gdl,*)
      real*8 dprop(gdl,*)
      integer numel
      integer nfile
      logical bvtk
      integer cod
      integer i,j
      character*15 cname
      character buffer*1024,lf*1,str1*8
      lf = char(10)
c ====================================================================== 
c === DADOS          
c ... VTk BINARY
      if(bvtk)then
c ... int      
        if(cod .eq. 1)then
          write(str1(1:8),'(i8)')gdl
          buffer = 'SCALARS '//cname//' int'//str1//lf
          write(nfile)trim(buffer)
          buffer = 'LOOKUP_TABLE default'//lf
          write(nfile)trim(buffer)
          do i = 1, numel
             write(nfile)(iprop(j,i),j=1,gdl)
          enddo
c ... float          
        else if(cod .eq. 2)then  
          write(str1(1:8),'(i8)')gdl
          buffer = 'SCALARS '//cname//'float'//str1//lf
          write(nfile)trim(buffer)
          buffer = 'LOOKUP_TABLE default'//lf
          write(nfile)trim(buffer)
          do i = 1, numel
            write(nfile)(fprop(j,i),j=1,gdl)
          enddo
c ... double         
        else if(cod .eq. 3)then  
          write(str1(1:8),'(i8)')gdl
          buffer = 'SCALARS '//cname//'double'//str1//lf
          write(nfile)trim(buffer)
          buffer = 'LOOKUP_TABLE default'//lf
          write(nfile)trim(buffer)
          do i = 1, numel
            write(nfile)(dprop(j,i),j=1,gdl)
          enddo
        endif
c ... Vtk ASCII        
      else 
c ... int      
        if(cod .eq. 1)then
         write(nfile,*)'SCALARS ',cname,'int',gdl
c          write(nfile,'(a,2x,15a,3a,i)')'SCALARS',cname,'int ',gdl
          write(nfile,'(a)')'LOOKUP_TABLE default'
          do i = 1, numel
            write(nfile,'(7i10)')(iprop(j,i),j=1,gdl)
          enddo
c ... float          
        else if(cod .eq. 2)then  
          write(nfile,'(a,1x,15a,a,1x,i4)')'SCALARS',cname,'float',gdl
          write(nfile,'(a)')'LOOKUP_TABLE default'
          do i = 1, numel
            write(nfile,'(7f10.5)')(fprop(j,i),j=1,gdl)
          enddo
c ... double         
        else if(cod .eq. 3)then  
          write(nfile,'(a,1x,15a,a,1x,i4)')'SCALARS',cname,'double',gdl
          write(nfile,'(a)')'LOOKUP_TABLE default'
          do i = 1, numel
            write(nfile,'(7f10.5)')(dprop(j,i),j=1,gdl)
          enddo
        endif
      endif  
c ====================================================================== 
c
c ===
      return
      end
c ====================================================================== 
c **********************************************************************
c
c **********************************************************************
c * Data de criacao    : 00/00/0000                                    *
c * Data de modificaco : 00/00/0000                                    * 
c * ------------------------------------------------------------------ *  
c * CELL_PROP_VTU : escreve propriedades das celulas (xlm)             *
c * -------------------------------------------------------------------*
c * Parametros de entrada:                                             *
c * -------------------------------------------------------------------*
c * dprop - propriedas real 8 bytes                                    *
c * iprop - propriedas interia de 8 bytes                              *
c * fprop - propriedas real de 4 bytes                                 *
c * cod   - 1 - interio de  4 bytes                                    *  
c *       - 2 - real de 4 bytes                                        *  
c *       - 3 - real de 8 bytes                                        *  
c * gdl   - graus de liberdade                                         *  
c * bvtk  - escrita no formato vtk binario                             *  
c * nfile - numero associado ao arquivo de saida                       *
c * -------------------------------------------------------------------*
c * Parmetros de saida:                                                *
c * -------------------------------------------------------------------*
c **********************************************************************
      subroutine cell_prop_vtu(iprop,fprop,dprop,numel,cname,cod,gdl 
     .                        ,bvtk,nfile)
c ===
      implicit none
      integer gdl
      integer iprop(gdl,*)
      real*4 fprop(gdl,*)
      real*8 dprop(gdl,*)
      integer numel
      integer nfile
      logical bvtk
      integer cod
      integer i,j
      character*15 cname
      character buffer*1024,buffer1*80,buffer2*80,buffer3*80,str1*8
c ====================================================================== 
c === DADOS          
c ... VTk BINARY
      if(bvtk)then
c ... int      
        if(cod .eq. 1)then
c ... float          
        else if(cod .eq. 2)then  
c ... double         
        else if(cod .eq. 3)then  
        endif
c ... Vtk ASCII        
      else 
c ... int      
        if(cod .eq. 1)then
          write(str1,'(i1)')gdl
          buffer1= '<DataArray type="Int32"'
          buffer2= ' NumberOfComponents="'//trim(str1)//'"'
          buffer3= ' Name="' // trim(cname) // '" format="ascii">'
          buffer=trim(buffer1)//trim(buffer2)//trim(buffer3)
          write(nfile,'(a)')trim(buffer)
          do i = 1, numel
            write(nfile,'(15i9)')(iprop(j,i),j=1,gdl)
          enddo
c ... float          
        else if(cod .eq. 2)then  
c ... double         
        else if(cod .eq. 3)then  
        endif
        write(nfile,'(a)')'</DataArray>'
      endif  
c ====================================================================== 
c
c ===
      return
      end
c ====================================================================== 
c **********************************************************************
c
c **********************************************************************
c * Data de criacao    : 00/00/0000                                    *
c * Data de modificaco : 09/04/2016                                    * 
c * ------------------------------------------------------------------ *  
c * HEAD_VTK : escreve o vabecalho do arquivo vtk                      *
c * -------------------------------------------------------------------*
c * Parametros de entrada:                                             *
c * -------------------------------------------------------------------*
c * headfile - cabecalho do arquivo vtk                                *
c * bvtk     - true formato binary false ascii                         *
c * t        - tempo real                                              *
c * istep    - passo de tempo                                          *
c * nfile    - numero associado ao arquivo de saida                    *
c * -------------------------------------------------------------------*
c * Parmetros de saida:                                                *
c * -------------------------------------------------------------------*
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ * 
c **********************************************************************
      subroutine head_vtk(headfile,bvtk,t,istep,time,nfile)
c ===
      implicit none
      character*30 headfile
      character lf*1, buffer*1024
      character buffer1*80,buffer2*80
      integer nfile
      integer istep
      real*8 t
      logical bvtk,time
      lf =char(10)
c ====================================================================== 
c
c ===
c ... legacy binario
      if(bvtk)then
        buffer ='# vtk DataFile Version 3.0'//lf
        write(nfile) trim(buffer)
        write(nfile) headfile
        buffer =lf//'BINARY'//lf
        write(nfile) trim(buffer)
        buffer = 'DATASET UNSTRUCTURED_GRID'//lf
        write(nfile) trim(buffer)
        buffer = 'FIELD FieldData 3'//lf
        write(nfile) trim(buffer)
        buffer = 'TIME(s) 1 1 double'//lf
        write(nfile) trim(buffer)
        write(nfile) t           
        buffer = 'TIME(h) 1 1 double'//lf
        write(nfile) trim(buffer)
        write(nfile) t/3600.0d0           
        buffer = 'CYCLE 1 1 int'//lf
        write(nfile) trim(buffer)
        write(nfile) istep 
c ... legacy ascii  
      else
        write(nfile,'(a)') '# vtk DataFile Version 3.0'
        write(nfile,'(a)')headfile 
        write(nfile,'(a)') 'ASCII'
        write(nfile,'(a)') 'DATASET UNSTRUCTURED_GRID'
        if(time) then
          write(nfile,'(a)') 'FIELD FieldData 3'
          write(nfile,'(a)') 'TIME(s) 1 1 double'
          write(nfile,'(f16.6)') t
          write(nfile,'(a)') 'TIME(h) 1 1 double'
          write(nfile,'(f16.6)') t/3600.0d0
          write(nfile,'(a)') 'CYCLE 1 1 int'
          write(nfile,'(i9)') istep
        endif  
      endif
c ======================================================================
      return
      end
c ======================================================================
c **********************************************************************
c
c **********************************************************************
c * Data de criacao    : 00/00/0000                                    *
c * Data de modificaco : 00/00/0000                                    * 
c * ------------------------------------------------------------------ *  
c * HEAD_VTU : escreve o cabecalho do arquivo vtu (xml)                *
c * -------------------------------------------------------------------*
c * Parametros de entrada:                                             *
c * -------------------------------------------------------------------*
c * nnode    - numero de nos                                           *
c * numel    - numero de elementos                                     *
c * bvtk     - true formato binary false ascii                         *
c * nfile    - numero associado ao arquivo de saida                    *
c * -------------------------------------------------------------------*
c * Parmetros de saida:                                                *
c * -------------------------------------------------------------------*
c **********************************************************************
      subroutine head_vtu(nnode,numel,bvtk,nfile)
c ===
      implicit none
      integer nnode,numel
      character lf*1, buffer*1024
      character buffer1*80,buffer2*80,str1*15,str2*15
      integer nfile
      logical bvtk
      lf =char(10)
c ====================================================================== 
c
c ===
c ... xml binario
      if(bvtk)then
        print*,'não implementado.'
        stop
c ... xml ascii  
      else
        write(nfile,'(a)') '<?xml version="1.0"?>'
        buffer1='<VTKFile type="UnstructuredGrid" '
        buffer2=' version="0.1" byte_order="BigEndian">'
        buffer = trim(buffer1) // trim(buffer2)
        write(nfile,'(a)') trim(buffer)
        buffer='<UnstructuredGrid>'
        write(nfile,'(a)') trim(buffer)
        write(str1(1:15),'(i15)') nnode
        write(str2(1:15),'(i15)') numel
        buffer1='<Piece NumberOfPoints="'//trim(str1)//'" '
        buffer2=' NumberOfCells="'//trim(str2)//'">'
        buffer = trim(buffer1) // trim(buffer2)
        write(nfile,'(a)') trim(buffer)
      endif
c ======================================================================
      return
      end
c ======================================================================
c **********************************************************************
c
c **********************************************************************
c * Data de criacao    : 00/00/0000                                    *
c * Data de modificaco : 00/00/0000                                    * 
c * ------------------------------------------------------------------ *  
c * CELL_DATA_VTK: incia uma secao CELL_DATA                           *
c * -------------------------------------------------------------------*
c * Parametros de entrada:                                             *
c * -------------------------------------------------------------------*
c * numel    - numero de elementos                                     *
c * bvtk     - true formato binary false ascii                         *
c * nfile    - numero associado ao arquivo de saida                    *
c * -------------------------------------------------------------------*
c * Parmetros de saida:                                                *
c * -------------------------------------------------------------------*
c **********************************************************************
      subroutine cell_data_vtk(numel,bvtk,nfile)
      implicit none
      integer numel
      integer nfile
      logical bvtk
      character str1*8,buffer*1024,lf*1
      lf = char(10)
c ... VTK BINARY
      if(bvtk)then
        write(str1(1:8),'(i8)')numel
        buffer = lf//'CELL_DATA '//str1//lf
        write(nfile) trim(buffer)
c ... VTK ASCII 
      else
        write(nfile,'(a,1x,i9)')'CELL_DATA',numel
      endif
c........................................................................
c 
c ... 
      return
      end
c **********************************************************************
c
c **********************************************************************
c * Data de criacao    : 00/00/0000                                    *
c * Data de modificaco : 00/00/0000                                    * 
c * ------------------------------------------------------------------ *  
c * CELL_DATA_VTU : inicia uma secao CELL_DATA (xlm)                   *
c * -------------------------------------------------------------------*
c * Parametros de entrada:                                             *
c * -------------------------------------------------------------------*
c * bvtk     - true formato binary false ascii                         *
c * nfile    - numero associado ao arquivo de saida                    *
c * -------------------------------------------------------------------*
c * Parmetros de saida:                                                *
c * -------------------------------------------------------------------*
c **********************************************************************
      subroutine cell_data_vtu(bvtk,nfile)
      implicit none
      integer nfile
      logical bvtk
      character str1*8,buffer*1024,lf*1
      lf = char(10)
c ... xlm BINARY
      if(bvtk)then
c ... xlm ASCII 
      else
        write(nfile,'(a)')'<CellData Scalars="scalars">'
      endif
c........................................................................
c 
c ... 
      return
      end
c **********************************************************************
c       
c **********************************************************************
c * Data de criacao    : 00/00/0000                                    *
c * Data de modificaco : 00/00/0000                                    * 
c * ------------------------------------------------------------------ *  
c * CELL_DATA_FINALIZE_VTU : fecha uma secao CELL_DATA (xlm)           *
c * -------------------------------------------------------------------*
c * Parametros de entrada:                                             *
c * -------------------------------------------------------------------*
c * bvtk     - true formato binary false ascii                         *
c * nfile    - numero associado ao arquivo de saida                    *
c * -------------------------------------------------------------------*
c * Parmetros de saida:                                                *
c * -------------------------------------------------------------------*
c **********************************************************************
      subroutine cell_data_finalize_vtu(bvtk,nfile)
      implicit none
      integer nfile
      logical bvtk
c ... xlm BINARY
      if(bvtk)then
c ... xlm ASCII 
      else
        write(nfile,'(a)')'</CellData>'
      endif
c........................................................................
c 
c ... 
      return
      end
c **********************************************************************
c
c **********************************************************************
c * Data de criacao    : 00/00/0000                                    *
c * Data de modificaco : 00/00/0000                                    * 
c * ------------------------------------------------------------------ *  
c * POINT_DATA_FINALIZE_VTU : fecha uma secao POINT_DATA (xlm)         *
c * -------------------------------------------------------------------*
c * Parametros de entrada:                                             *
c * -------------------------------------------------------------------*
c * bvtk     - true formato binary false ascii                         *
c * nfile    - numero associado ao arquivo de saida                    *
c * -------------------------------------------------------------------*
c * Parmetros de saida:                                                *
c * -------------------------------------------------------------------*
c **********************************************************************
      subroutine point_data_finalize_vtu(bvtk,nfile)
      implicit none
      integer nfile
      logical bvtk
c ... xlm BINARY
      if(bvtk)then
c ... xlm ASCII 
      else
        write(nfile,'(a)')'</PointData>'
      endif
c........................................................................
c 
c ... 
      return
      end
c **********************************************************************
c
c **********************************************************************
c * Data de criacao    : 00/00/0000                                    *
c * Data de modificaco : 00/00/0000                                    * 
c * ------------------------------------------------------------------ *  
c * FINALIZE_VTU : finaliza arquivo vtu                                *
c * -------------------------------------------------------------------*
c * Parametros de entrada:                                             *
c * -------------------------------------------------------------------*
c * bvtk     - true formato binary false ascii                         *
c * nfile    - numero associado ao arquivo de saida                    *
c * -------------------------------------------------------------------*
c * Parmetros de saida:                                                *
c * -------------------------------------------------------------------*
c **********************************************************************
      subroutine finalize_vtu(bvtk,nfile)
      implicit none
      integer nfile
      logical bvtk
c ... xlm BINARY
      if(bvtk)then
c ... xlm ASCII 
      else
        write(nfile,'(a)') '</Piece>'
        write(nfile,'(a)') '</UnstructuredGrid>'
        write(nfile,'(a)') '</VTKFile>'
      endif
c........................................................................
c 
c ... 
      return
      end
c **********************************************************************
c       
c       
c **********************************************************************
c * Data de criacao    : 00/00/0000                                    *
c * Data de modificaco : 00/00/0000                                    * 
c * ------------------------------------------------------------------ *  
c * POINT_DATA: incia um secao POINT_DATA                              *
c * -------------------------------------------------------------------*
c * Parametros de entrada:                                             *
c * -------------------------------------------------------------------*
c * nnode    - numero de pontos                                        *
c * bvtk     - true formato binary false ascii                         *
c * nfile    - numero associado ao arquivo de saida                    *
c * ------------------------------------------------------------------ *
c * Parmetros de saida:                                                *
c * -------------------------------------------------------------------*
c **********************************************************************
      subroutine point_data_vtk(nnode,bvtk,nfile)
      implicit none
      integer nnode
      integer nfile
      logical bvtk
      character str1*8,buffer*1024,lf*1
      lf = char(10)
c ... VTK BINARY
      if(bvtk)then
        write(str1(1:8),'(i8)')nnode
        buffer = lf//' POINT_DATA '//str1//lf
        write(nfile) trim(buffer)
c ... VTK ASCII 
      else
        write(nfile,'(a,1x,i9)')'POINT_DATA',nnode
      endif
c........................................................................
c 
c ... 
      return
      end
c **********************************************************************
c
c **********************************************************************
c * Data de criacao    : 00/00/0000                                    *
c * Data de modificaco : 00/00/0000                                    * 
c * ------------------------------------------------------------------ *  
c * POINT_DATA_VTU: incia uma secao POINT_DATA  (xlm)                  *
c * -------------------------------------------------------------------*
c * Parametros de entrada:                                             *
c * -------------------------------------------------------------------*
c * bvtk     - true formato binary false ascii                         *
c * nfile    - numero associado ao arquivo de saida                    *
c * -------------------------------------------------------------------*
c * Parmetros de saida:                                                *
c * -------------------------------------------------------------------*
c **********************************************************************
      subroutine point_data_vtu(bvtk,nfile)
      implicit none
      integer nfile
      logical bvtk
c ... xlm BINARY
      if(bvtk)then
c ... xlm ASCII 
      else
        write(nfile,'(a)')'<PointData>'
      endif
c........................................................................
c 
c ... 
      return
      end
c **********************************************************************
c  
c **********************************************************************
c * Data de criacao    : 00/00/0000                                    *
c * Data de modificaco : 27/03/2016                                    * 
c * ------------------------------------------------------------------ *  
c * FACE_VTK:escreve faces no formato vtk com as suas respectivos      *
c * cargas                                                             *
c * -------------------------------------------------------------------*
c * Parametros de entrada:                                             *
c * -------------------------------------------------------------------*
c * face              - conetividade das faces                         *
c * numel             - numero de elementos                            *
c * ndm               - numeros de dimensoes                           *
c * nen               - numero maximo de nos por elemento              *
c * bvtk              - true formato binary false ascii                *
c * nfile             - numero associado ao arquivo de saida           *
c * -------------------------------------------------------------------*
c * Prametros de saida:                                                *
c * -------------------------------------------------------------------*
c * -------------------------------------------------------------------*
c * OBS:                                                               *
c * -------------------------------------------------------------------*
c * Os prametors nbar2,ntria3,nquad4,ntetra4 e nhexa8 foram herdados do*
c *   mefpar                                                           *
c *--------------------------------------------------------------------*
c **********************************************************************
      subroutine face_vtk(face       ,tipo_face,max_no_face
     .                   ,line_face  ,tria_face,quad_face
     .                   ,bvtk       ,nfile)
      implicit none
c ...
      integer numet,nface,nnoel
      integer line_face,tria_face,quad_face,max_no_face
      integer nfile
      integer face(max_no_face,*),tipo_face(*)
      character buffer*1024,lf*1,str1*15,str2*15
      logical bvtk
      integer i,j  
      integer nno_face(3),tipo_face_vtk(3)
      data nno_face /2,3,4/
      data tipo_face_vtk /3,5,9/
c ......................................................................
c
c ...
      lf =char(10)
c ......................................................................
c
c ... Calculo do numero de tipo de cada elemento 
      numet = 3*line_face + 4*tria_face + 5*quad_face 
      nface = line_face + tria_face+ quad_face
c
c ... total de elemntos e tamnhanho dos dados totais
      if(bvtk)then
        write(str1(1:15),'(i15)')nface
        write(str2(1:15),'(i15)')numet
        buffer = lf//'CELLS '//str1//str2//lf
        write(nfile) trim(buffer)  
      else
        write(nfile,'(a,1x,i10,1x,i10)') 'CELLS ',nface,numet  
      endif 
c ... escrevendo a malha
c     
c ... nos dos elementos
      do i = 1, nface
        nnoel = nno_face(tipo_face(i))
        if(bvtk)then
          write(nfile) nnoel,(face(j,i)-1,j=1,nnoel)
        else
          write(nfile,'(10i10)') nnoel,(face(j,i)-1,j=1,nnoel)
        endif  
      enddo
c ......................................................................
c
c ... tipo do face
      if(bvtk)then
        write(str1(1:15),'(i15)')nface
        buffer = lf//'CELL_TYPES '//str1//lf
        write(nfile) trim(buffer)
      else
        write(nfile,'(a,1x,i10)') 'CELL_TYPES ',nface
      endif
c ...  
      do i = 1, nface
        nnoel = tipo_face_vtk(tipo_face(i))
        if(bvtk)then
          write(nfile)nnoel
        else  
          write(nfile,'(i3)')nnoel
        endif
      enddo  
c ......................................................................
c
c ...
      return
      end
c ......................................................................
c **********************************************************************
c
c **********************************************************************
c * Data de criacao    : 00/00/0000                                    *
c * Data de modificaco : 00/00/0000                                    * 
c * ------------------------------------------------------------------ *  
c * FACE_VTU: escreve faces no formato vtu com as suas respectivos     *
c * cargas                                                             *
c * -------------------------------------------------------------------*
c * Parametros de entrada:                                             *
c * -------------------------------------------------------------------*
c * face              - conetividade das faces                         *
c * numel             - numero de elementos                            *
c * ndm               - numeros de dimensoes                           *
c * nen               - numero maximo de nos por elemento              *
c * bvtk              - true formato binary false ascii                *
c * nfile             - numero associado ao arquivo de saida           *
c * -------------------------------------------------------------------*
c * Prametros de saida:                                                *
c * -------------------------------------------------------------------*
c * -------------------------------------------------------------------*
c * OBS:                                                               *
c * -------------------------------------------------------------------*
c * Os prametors nbar2,ntria3,nquad4,ntetra4 e nhexa8 foram herdados do*
c *   mefpar                                                           *
c *--------------------------------------------------------------------*
c **********************************************************************
      subroutine face_vtu(face     ,carga   ,tipoface
     .                   ,maxnoface,lineface,triaface
     .                   ,quadface ,bvtk    ,nfile)
      implicit none
c ...
      integer numet,nface,nnoel,offset
      integer lineface,triaface,quadface,maxnoface
      integer nfile
      integer face(maxnoface,*),tipoface(*),carga(*)
      character buffer*1024,lf*1,str1*15,str2*15,buffer1*80,buffer2*80
      logical bvtk
      integer i,j  
      integer nnoface(3),tipofacevtk(3)
      data nnoface /2,3,4/
      data tipofacevtk /3,5,9/
c ......................................................................
c
c ...
      lf =char(10)
c ......................................................................
c
c ... Calculo do numero de tipo de cada elemento 
      numet = 3*lineface + 4*triaface + 5*quadface 
      nface = lineface + triaface+ quadface
c
c ... total de elemntos e tamnhanho dos dados totais
      if(bvtk)then
        print*,'nao implementado.'
        stop
      else
        write(nfile,'(a)') '<Cells> '
        buffer1 ='<DataArray type="Int32" Name="connectivity" '
        buffer2 =' format="ascii">'
        buffer = trim(buffer1) // trim(buffer2)
        write(nfile,'(a)') trim(buffer)  
      endif 
c ... escrevendo a malha
c     
c ... nos dos elementos
      do i = 1, nface
        nnoel = nnoface(tipoface(i))
        if(bvtk)then
          print*,'nao implementado.'
          stop
        else
          write(nfile,'(10i10)') (face(j,i)-1,j=1,nnoel)
        endif  
      enddo
c ......................................................................
c
c ...
      if(bvtk)then
      else
        write(nfile,'(a)') '</DataArray>'
      endif
c ......................................................................
c
c ... offset no arranjo das conectvidade(funcionando para malha
c                         de um elemento apenas)
      offset = 0
      if(bvtk)then
        print*,'nao implementado.'
        stop
      else
        buffer1 ='<DataArray type="Int64" Name="offsets" '
        buffer2 =' format="ascii">'
        buffer = trim(buffer1) // trim(buffer2)
        write(nfile,'(a)') trim(buffer)
        do i = 1, nface
          nnoel  = nnoface(tipoface(i))
          offset = nnoel + offset
          write(nfile,'(i20)') offset
        enddo  
        write(nfile,'(a)') '</DataArray>'
      endif 
c .....................................................................
c
c ... tipo do face
      if(bvtk)then
        print*,'nao implementado.'
        stop
      else
        buffer1 ='<DataArray type="UInt8" Name="types" '
        buffer2 =' format="ascii">'
        buffer = trim(buffer1) // trim(buffer2)
        write(nfile,'(a)') trim(buffer)
      endif  
c ...  
      do i = 1, nface
        nnoel = tipofacevtk(tipoface(i))
        if(bvtk)then
          print*,'nao implementado.'
          stop
        else  
          write(nfile,'(i3)')nnoel
        endif
      enddo
c ......................................................................
      if(bvtk)then
        print*,'nao implementado.'
        stop
      else
        write(nfile,'(a)') '</DataArray>'
        write(nfile,'(a)') '</Cells>'
      endif  
c ......................................................................
c
c ...
      return
      end
c ......................................................................
c **********************************************************************

