      subroutine pdata(ix,x,u,nnode,numel,ndm,nen,ndf,istep,nplot,code)
c **********************************************************************
c *                                                                    *
c *   PDATA: Imprime arquivos de saida                                 *
c *   -----                                                            *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ----------------------                                           *
c *                                                                    *
c *    ix(nen+1,numel) - conetividades nodais                          *
c *    x(ndm,nnode)    - coordenadas nodais                            *
c *    u(ndf,nnode)    - campo escalar ou vetorial                     *
c *    nnode  - numero de nos                                          *
c *    numel  - numero de elementos                                    *
c *    nbar2  - numero de barras                                       *
c *    ntria3 - numero de triangulos                                   *
c *    nquad4 - numero de quadrilateros                                *
c *    ntetra4 - numero de tetraedros                                  *
c *    npenta5 - numero de elementos cunha                             *
c *    nhexa8  - numero de hexaedros                                   *
c *    ndm - dimensao                                                  *
c *    nen - numero maximo de nos por elementos                        *
c *    ndf - dimensao de u                                             *
c *    nplot - arquivo de saida                                        *
c *    cont1 - controle de impressao vtk                               *
c *    code - codigo de instrucao:                                     *
c *       ---  Geometria ( valores < 0 )  ---                          *
c *         = -1 - geometria View3D                                    *
c *         = -2 - geometria Gid                                       *
c *         = -3 - geometria ParaView                                  *
c *       ---  Solucao ( valores > 0 )  ---                            *
c *         =  1 - Temperatura View3D                                  *
c *         =  2 - Temperatura Gid                                     *
c *         =  3 - Temperatura ParaView                                *
c *         =  4 - Fluxo de Calor View3D                               *
c *         =  5 - Fluxo de Calor Gid                                  *
c *         =  6 - Fluxo de Calor ParaView                             *
c *         =  7 - Deslocamentos View3D                                *
c *         =  8 - Deslocamentos Gid                                   *
c *         =  9 - Deslocamentos ParaView                              *
c *         =  10 - Tensoes View3D                                     *
c *         =  11 - Tensoes Gid                                        *
c *         =  12 - Tensoes ParaView                                   *
c *                                                                    *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   --------------------                                             *
c *                                                                    *
c *                                                                    *
c **********************************************************************
      implicit none
      include 'elementos.fi'
      integer ix(nen+1,*),nnode,numel,ndm,nen,ndf,nplot,code,istep
c      integer ix(21,*),nnode,numel,ndm,nen,ndf,nplot,code,istep
      integer numelt
      integer i,j,k
      real*8  x(ndm,*),u(ndf,*),scale
      logical cont1
c ======================================================================
c
c ...    Arquivo de saida no formato VIEW3D v. 2.01
c
      scale = 0.0e+00
      if(code .eq. 3 .or. code .eq. 6 .or. code .eq. 9 
     .   .or. code .eq. 12 .or. code .eq. 17 .or. code .eq. 18)  goto 40
c
c...      Geometria:
   41 continue
      if(code .lt. 0) then
c ======================================================================
c
c ... Geometria para o View3D 2.01:
c     
         if(code .eq. -1) then
c ======================================================================
c        Coordenadas:
c ======================================================================
            write(nplot,'(a,i10)') 'coor ', nnode
            if (abs(code) .eq. 1) then
               if(ndm .eq. 1) then
                  do 5 i = 1, nnode
                     write(nplot,'(i10,3e15.5)') i,x(1,i),0.0,0.0
    5             continue         
               elseif(ndm .eq. 2) then
                  do 10 i = 1, nnode
                     write(nplot,'(i10,3e15.5)') i,(x(j,i),j=1,ndm),0.0
   10             continue
               elseif (ndm .eq. 3) then
                  do 11 i = 1, nnode
                       write(nplot,'(i10,3e15.5)') i,(x(j,i),j=1,ndm)
   11                  continue
               endif
c ...    Deformada:            
            elseif(abs(code) .eq. 2) then
               if(ndm .eq. 2) then
                  do 12 i = 1, nnode
                     write(nplot,'(i10,3e15.5)') i,x(1,i)+scale*u(1,i),
     .                                          x(2,i)+scale*u(2,i),0.0
   12             continue
               elseif (ndm .eq. 3) then
                  do 13 i = 1, nnode
                     write(nplot,'(i10,3e15.5)') i,x(1,i)+scale*u(1,i),
     .                                             x(2,i)+scale*u(2,i),
     .                                             x(3,i)+scale*u(3,i)
   13             continue
               endif
            endif
c ======================================================================
c ...    Conetividades dos elementos:
c ======================================================================
            if (nbar2(1) .gt. 0) then
               write(nplot,'(a,i10)') 'bar2 ', nbar2(1)
               do 14 i = nbar2(2), nbar2(2) + nbar2(1)-1
                  write(nplot,'(9i10)') i,(ix(j,i),j=1,2)
   14          continue
            endif   
c ......................................................................
            if (ntria3(1) .gt. 0) then
                  write(nplot,'(a,i10)') 'tria3', ntria3(1)
               do 15 i = ntria3(2), ntria3(2) + ntria3(1)-1
                  write(nplot,'(9i10)') i,(ix(j,i),j=1,3)
   15          continue
            endif
c ......................................................................   
            if (nquad4(1) .gt. 0) then
                  write(nplot,'(a,i10)') 'quad4', nquad4(1)
               do 16 i = nquad4(2), nquad4(2) + nquad4(1)-1
                  write(nplot,'(9i10)') i,(ix(j,i),j=1,4)
   16          continue
            endif
c ......................................................................            
            if (ntetra4(1) .gt. 0) then
                  write(nplot,'(a,i10)') 'tetra4', ntetra4(1)
               do 17 i = ntetra4(2), ntetra4(2) + ntetra4(1)-1
                  write(nplot,'(9i10)') i,(ix(j,i),j=1,4)
   17          continue            
            endif
c ......................................................................            
            if (nhexa8(1) .gt. 0) then
                  write(nplot,'(a,i10)') 'hexa8 ', nhexa8(1)
               do 18 i = nhexa8(2), nhexa8(2)+nhexa8(1)-1
                  write(nplot,'(9i10)') i,(ix(j,i),j=1,8)
   18          continue
            endif
c ======================================================================
c
c ... Geometria para o Gid:
c ======================================================================
         elseif(code .eq. -2) then
c ======================================================================
c        Coordenadas:
c ======================================================================
             if (nbar2(1) .gt. 0) then
               write(nplot,'(a,i5,a)') 'MESH " "  dimension ',ndm,
     .                                 ' ElemType Linear Nnode 2'
               write(nplot,'(a)')'Coordinates'
               do i = 1, nnode
                  write(nplot,'(i10,3e15.5)') i,(x(j,i),j=1,ndm)
               enddo
               write(nplot,'(a)') 'End coordinates'
               write(nplot,'(a)') 'Elements '
               do 19 i = nbar2(2), nbar2(2) + nbar2(1)-1
                  write(nplot,'(9i10)') i,(ix(j,i),j=1,3)
   19          continue
               write(nplot,'(a)') 'End elements '
             endif
c ......................................................................
            if (ntria3(1) .gt. 0) then
               write(nplot,'(a,i5,a)') 'MESH " "  dimension',ndm,
     .                             ' ElemType Triangle Nnode 3'
               write(nplot,'(a)') 'Coordinates'
               do i = 1, nnode
                  write(nplot,'(i10,3e15.5)') i,(x(j,i),j=1,ndm)
               enddo
               write(nplot,'(a)') 'End coordinates'
               write(nplot,'(a)') 'Elements '
               do 20 i = ntria3(2), ntria3(2) + ntria3(1)-1
                  write(nplot,'(10i10)') i,(ix(j,i),j=1,4)
   20          continue
               write(nplot,'(a)') 'End elements '
            endif
c ......................................................................
            if (nquad4(1) .gt. 0) then
               write(nplot,'(a,i5,a)') 'MESH " " dimension',ndm,
     .                             ' ElemType Quadrilateral Nnode 4'
               write(nplot,'(a)') 'Coordinates'
               do i = 1, nnode
                  write(nplot,'(i10,3e15.5)') i,(x(j,i),j=1,ndm)
               enddo
               write(nplot,'(a)') 'End coordinates'
               write(nplot,'(a)') 'Elements '
               do 21 i = nquad4(2), nquad4(2) + nquad4(1)-1
                  write(nplot,'(10i10)') i,(ix(j,i),j=1,5)
   21          continue
               write(nplot,'(a)') 'End elements '
            endif
c ......................................................................
            if (ntetra4(1) .gt. 0) then
               write(nplot,'(a,i5,a)') 'MESH " " dimension',ndm,
     .                             ' ElemType Tetrahedra Nnode 4'
               write(nplot,'(a)') 'Coordinates'
               do i = 1, nnode
                  write(nplot,'(i10,3e15.5)') i,(x(j,i),j=1,ndm)
               enddo
               write(nplot,'(a)') 'End coordinates'
               write(nplot,'(a)') 'Elements '
               do 22 i = ntetra4(2), ntetra4(2) + ntetra4(1)-1
                  write(nplot,'(10i10)') i,(ix(j,i),j=1,5)
   22          continue
               write(nplot,'(a)') 'End elements '
            endif
c ......................................................................
            if (nhexa8(1) .gt. 0) then
               write(nplot,'(a,i5,a)') 'MESH " " dimension',ndm,
     .                             ' ElemType Hexahedra Nnode 8'
               write(nplot,'(a)') 'Coordinates'
               do i = 1, nnode
                  write(nplot,'(i10,3e15.5)') i,(x(j,i),j=1,ndm)
               enddo
               write(nplot,'(a)') 'End coordinates'
               write(nplot,'(a)') 'Elements '
               do 23 i = nhexa8(2), nhexa8(2)+nhexa8(1)-1
                  write(nplot,'(10i10)') i,(ix(j,i),j=1,9)
   23          continue
               write(nplot,'(a)') 'End elements '
            endif
c ======================================================================
c
c ... Geometria para o ParaView:
c ======================================================================
         elseif(code .eq. -3)  then
            go to 40
         endif
         return
c ======================================================================
c ... Temperatura para o View3D 2.01:
c ======================================================================
      elseif(code .eq. 1) then
         write(nplot,'(a,i5)') 'nosc ', 1  
         write(nplot,'(a)') 'Temperatura'
         do 28 i = 1, nnode
            write(nplot,'(i10,7e15.5e3)') i,u(ndf,i)
   28    continue 
c ======================================================================
c ... Temperatura para o Gid:
c ======================================================================
      elseif (code .eq. 2) then
         write(nplot,'(a,i8,a)') 'Result "Temperatura" " "',
     .         istep,' Scalar onNodes'
         write(nplot,'(a)') 'Values'
         do 29 i = 1, nnode
            write(nplot,'(i10,7e15.5e3)') i,(u(j,i),j=1,ndf)
   29    continue
         write(nplot,'(a)') 'End values'
      return
c ======================================================================
c ... Temperatura para o ParaView:
c ======================================================================
      elseif (code .eq. 3) then
            write(nplot,'(a,i10)') 'POINT_DATA', nnode
         write(nplot,'(a)') 'SCALARS Temperaturas double'
         write(nplot,'(a)') 'LOOKUP_TABLE default'
         do i = 1, nnode
            if (ndm .eq. 2)   then
               write(nplot,'(1e15.5e3)') (u(j,i),j=1,ndf)
            elseif (ndm .eq. 3)   then
               write(nplot,'(1e15.5e3)') (u(j,i),j=1,ndf)
            endif
         enddo
      return
c ======================================================================
c ... Fluxo de Calor para o View3D 2.01:
c ======================================================================
      elseif(code .eq. 4) then
         write(nplot,'(a,i5)') 'novec '
         if (ndf .eq. 1) then
            do 32 i = 1, nnode
               write(nplot,'(i10,3e15.5e3)') i,u(1,i),0.d0,0.d0
   32       continue
         elseif (ndf .eq. 2) then
            do 35 i = 1, nnode
               write(nplot,'(i10,3e20.8e3)') i,(u(j,i),j=1,ndm),0.d0
   35       continue   
         elseif (ndf .eq. 3) then
            do 70 i = 1, nnode
               write(nplot,'(i10,3e20.8e3)') i,(u(j,i),j=1,ndm)
   70       continue
         endif
         return
c ======================================================================
c ... Fluxo de Calor para o Gid:
c ======================================================================
      elseif (code .eq. 5) then     
         write(nplot,'(a,i8,a)') 'Result "Fluxo" " "',
     .         istep,' Vector onNodes'
         write(nplot,'(a)') 'Values'
         do i = 1, nnode
            write(nplot,'(i10,7e15.5e3)') i,(u(j,i),j=1,ndm)
         enddo
         write(nplot,'(a)') 'End values'
         return
c ======================================================================
c ... Fluxo de Calor para o ParaView:
c ======================================================================
      elseif (code .eq. 6) then
            write(nplot,'(a,i8)') 'POINT_DATA', nnode
         write(nplot,'(a)') 'VECTORS FluxoDeCalor double'
         do i = 1, nnode
            if (ndm .eq. 2)   then
               write(nplot,'(3e15.5e3)') (u(j,i),j=1,ndm), 0.0e+00
            elseif (ndm .eq. 3)   then
               write(nplot,'(3e15.5e3)') (u(j,i),j=1,ndm)
            endif
         enddo
         return
c ======================================================================
c ... Deslocamentos para o View3D 2.01:
c ======================================================================
      elseif(code .eq. 7) then
         write(nplot,'(a,i5)') 'nosc ', ndf
         write(nplot,'(a)') 'u'
         if (ndf .gt. 1) write(nplot,'(a)') 'v'
         if(ndf .eq. 3) write (nplot,'(a)') 'w'
         do 31 i = 1, nnode
            write(nplot,'(i10,6e15.5e3)') i,(u(j,i),j=1,ndf)
   31    continue
         return
c ======================================================================
c ... Deslocamentos para o Gid:
c ======================================================================
      elseif (code .eq. 8) then
         write(nplot,'(a,i8,a)') 'Result "Deslocamentos" " "',
     .         istep,' Vector onNodes'
         write(nplot,'(a)') 'Values'
         do i = 1, nnode
            write(nplot,'(i10,7e15.5e3)') i,(u(j,i),j=1,ndf)
         enddo
         write(nplot,'(a)') 'End values'
         return
c ======================================================================
c ... Deslocamentos para o ParaView:
c ======================================================================
      elseif (code .eq. 9) then
         write(nplot,'(a,i8)') 'POINT_DATA', nnode
c         write(nplot,'(a,a,a)') 'VECTORS ', ' Deslocamentos ', ' double'
         write(nplot,'(a,a,a)') 'SCALARS ',' up            ',' double 4'
         write(nplot,'(a)') 'LOOKUP_TABLE default'
         do i = 1, nnode
            if (ndm .eq. 2)   then
               write(nplot,'(7e15.5e3)') (u(j,i),j=1,ndf), 0.0e+00
            elseif (ndm .eq. 3)   then
               write(nplot,'(7e15.5e3)') (u(j,i),j=1,ndf)
            endif
         enddo
         write(nplot,'(a,a,a)') 'VECTORS ', ' Deslocamentos ', ' double'
         do i = 1, nnode
            if (ndm .eq. 2)   then
               write(nplot,'(7e15.5e3)') (u(j,i),j=1,ndf), 0.0e+00
            elseif (ndm .eq. 3)   then
               write(nplot,'(7e15.5e3)') (u(j,i),j=1,ndf-1)
            endif
         enddo
         return
c ======================================================================
c ... Tensoes para o View3D 2.01:
c ======================================================================
      elseif(code .eq. 10) then
         write(nplot,'(a,i5)') 'nosc ', ndf
         write(nplot,'(a)') 'tx'
         write(nplot,'(a)') 'ty'
         if (ndf .eq. 3) then
            write(nplot,'(a)') 'txy'
         elseif (ndf .eq. 4) then
            write(nplot,'(a)') 'tz'
            write(nplot,'(a)') 'txy'
            write(nplot,'(a)') 'Tensao efetiva'
         elseif (ndf .eq. 6) then
            write(nplot,'(a)') 'tz'
            write(nplot,'(a)') 'txy'
            write(nplot,'(a)') 'tyz'
            write(nplot,'(a)') 'txz'
            write(nplot,'(a)') 'Tensao efetiva'
         endif
         do 50 i = 1, nnode
            write(nplot,'(i10,7e15.5e3)') i,(u(j,i),j=1,ndf)
   50    continue
         return
c ======================================================================
c ... Tensoes para o Gid:
c ======================================================================
      elseif (code .eq. 11) then
         write(nplot,'(a,i8,a)') 'Result "Tensao" " "',
     .         istep,' Vector onNodes'
         write(nplot,'(a)') 'Values'
         do i = 1, nnode
            write(nplot,'(i10,7e15.5e3)') i,(u(j,i),j=1,ndf)
         enddo
         write(nplot,'(a)') 'End values'
         return
c ======================================================================
c ... Tensoes para o ParaView:
c ======================================================================
      elseif (code .eq. 12) then
            write(nplot,'(a,i8)') 'POINT_DATA', nnode
         write(nplot,'(a,a,a)') 'TENSORS ', ' Tensoes ', ' double'
         do i = 1, nnode
            if (ndm .eq. 2)   then
               write(nplot,'(7e15.5e3)') u(1,i), u(4,i), 0.0e+00 , 
     .             u(4,i), u(2,i), 0.0e+00 , 0.0e+00, 0.0e+00 ,u(3,i) 
            elseif (ndm .eq. 3)   then
               write(nplot,'(7e15.5e3)') u(1,i), u(4,i), u(6,i), 
     .             u(4,i), u(2,i), u(5,i), u(6,i), u(5,i) , u(3,i)
            endif
         enddo
         return
c ======================================================================
c ... Grau de hidratacao para o ParaView:
c ======================================================================
      elseif (code .eq. 17) then
         write(nplot,'(a)') 'SCALARS hidratacao double'
         write(nplot,'(a)') 'LOOKUP_TABLE default'
         do i = 1, numel
            write(nplot,'(1e15.5e3)') u(3,i)
         enddo
      return
c ======================================================================
c ... Imprime tensao efetiva e criterio o ParaView:
c ======================================================================
      elseif (code .eq. 18) then
c         write(nplot,'(a,i8)') 'CELL_DATA', numel
         write(nplot,'(a,a,a)') 'VECTORS ', ' Tensao_efetiva', ' double'
         do i = 1, numel
            write(nplot,'(3e15.5e3)') u(1,i),u(2,i),0.d0
         enddo
      return
c ... Criterios de Ruptura por elemento:
c ======================================================================
c      elseif (code .eq. 18) then
c         write(nplot,'(a)') 'SCALARS plasticidade double'
c         write(nplot,'(a)') 'LOOKUP_TABLE default'
c         do i = 1, numel
c            if ( u(2,i) .gt. 0.d0 ) then
c               write(nplot,'(i8)')   2
c            elseif (u(1,i) .gt. 0.d0 ) then
c                write(nplot,'(i8)')   1
c            else   
c                write (nplot,'(i8)')   0
c            endif
c         enddo
c      return
c ======================================================================
c ... Escalares (tensoes principais) para o View3D 2.01:
c ======================================================================
      elseif(code .eq. 40) then
         write(nplot,'(a,i5)') 'nosc ', 3
         write(nplot,'(a)') 't1'
         write(nplot,'(a)') 't2'
         write(nplot,'(a)') 't3'
         do 60 i = 1, nnode
            write(nplot,'(i10,6e15.5e3)') i,(u(j,i),j=1,3)
   60    continue
c ======================================================================
c ... Escalares (deformacoes) para o View3D 2.01:
c ======================================================================
      elseif(code .eq. 50) then
         write(nplot,'(a,i5)') 'nosc ', ndf
         write(nplot,'(a)') 'ex'
         write(nplot,'(a)') 'ey'
         write(nplot,'(a)') 'ez'
         write(nplot,'(a)') 'exy'
         if (ndf .eq. 6) then
            write(nplot,'(a)') 'eyz'
            write(nplot,'(a)') 'exz'
         endif
         do 90 i = 1, nnode
            write(nplot,'(i10,7e15.5e3)') i,(u(j,i),j=1,ndf)
   90    continue
c ======================================================================
c ...    Arquivo de saida no formato VIEW3D v. 1.02
c ======================================================================
      elseif(code .eq. 111) then
         write (nplot,'(a)') 'Formato VIEW3D v. 1.02'
         write (nplot,'(a,i6)') 'coor ',nnode
         do 100 i = 1, nnode
            write (nplot,'(i10,3e15.5)') i,(x(j,i), j = 1,ndm)
  100    continue
         write (nplot,'(a,i6)') 'elem ',numel
         do 110 i = 1, numel
            write (nplot,'(8i10)') i,7,(ix(j,i), j = 1, nen)
  110    continue
         write(nplot,'(a,i2)') 'nosc ',ndf
         do 120 i = 1, nnode
            write (nplot,'(i10,3e15.5e3)') i,(u(j,i),j=1,ndf)
  120    continue
         write (nplot,'(a)') 'end '
c ======================================================================
c ... Transforma uma malha de hexaedros em tetraedros
c ======================================================================
      elseif(code .eq. 121) then
        write(nplot,'(a,i10)') 'tetra4 ', numel*6
        k = 1
        do 200 i = 1, nhexa8(1)
         write(nplot,'(6i10)') k,ix(2,i),ix(5,i),ix(6,i),ix(7,i),ix(9,i)
         k = k + 1
         write(nplot,'(6i10)') k,ix(2,i),ix(5,i),ix(7,i),ix(3,i),ix(9,i)
         k = k + 1
         write(nplot,'(6i10)') k,ix(2,i),ix(5,i),ix(3,i),ix(1,i),ix(9,i)
         k = k + 1
         write(nplot,'(6i10)') k,ix(3,i),ix(5,i),ix(7,i),ix(8,i),ix(9,i)
         k = k + 1
         write(nplot,'(6i10)') k,ix(3,i),ix(1,i),ix(8,i),ix(4,i),ix(9,i)
         k = k + 1
         write(nplot,'(6i10)') k,ix(3,i),ix(1,i),ix(5,i),ix(8,i),ix(9,i)
         k = k + 1
  200  continue
c ======================================================================
c ... Transforma uma malha de hexaedros em elementos cunha
c ======================================================================
      elseif(code .eq. 122) then
       k = 1
       do 300 i = 1, nhexa8(1)
         write(nplot,'(8i10)') k,ix(1,i),ix(2,i),ix(3,i),
     .                           ix(5,i),ix(6,i),ix(7,i),ix(9,i)
         k = k + 1
         write(nplot,'(8i10)') k,ix(4,i),ix(1,i),ix(3,i),
     .                           ix(8,i),ix(5,i),ix(7,i),ix(9,i)
         k = k + 1
  300  continue
      endif
      return
c ======================================================================
c        Comando de Geometria do Paraview:
c ======================================================================
   40 continue    
c ======================================================================
c        Coordenadas:
c ======================================================================
      write(nplot,'(a)') '# vtk DataFile Version 3.0' 
      write(nplot,'(a)') 'Malhas de elementos finitos'          
      write(nplot,'(a)') 'ASCII' 
      write(nplot,'(a)') 'DATASET UNSTRUCTURED_GRID' 
      write(nplot,'(a,i10,a)') 'POINTS  ',nnode, ' double'
      do i = 1, nnode
         if (ndm .eq. 2) then
            write(nplot,'(3e15.5)') (x(j,i),j=1,ndm), 0.0e+00
         else if (ndm .eq. 3) then
            write(nplot,'(3e15.5)') (x(j,i),j=1,ndm)
c                write(nplot,'(7e15.5e3)') u(1,i)*scale+x(1,i),
c     .          u(2,i)*scale+x(2,i),u(3,i)*scale+x(3,i)
         endif
      enddo
c =================================================================            
c Escrevendo os nos dos elementos
c =================================================================  
       numelt = 3*nbar2(1)   + 4*ntria3(1) + 5*nquad4(1) + 
     .         5*ntetra4(1) + 9*nhexa8(1) + 7 *nprism6(1)
c     .         5*ntetra4(1) + 21*nhexa8(1)  
      write(nplot,'(a, i10, i10)') 'CELLS ',numel, numelt
c ................................................................
       if (nbar2(1) .gt. 0) then
          do i = nbar2(2), nbar2(2) + nbar2(1)-1
            write(nplot,'(10i10)') 2,(ix(j,i)-1,j=1,2)
         enddo
      endif
c .................................................................
      if (ntria3(1) .gt. 0) then
          do 24 i = ntria3(2), ntria3(2) + ntria3(1)-1
            write(nplot,'(10i10)') 3,(ix(j,i)-1,j=1,3)
   24    continue
      endif
c ................................................................
      if (nquad4(1) .gt. 0) then
          do 25 i = nquad4(2), nquad4(2) + nquad4(1)-1
            write(nplot,'(10i10)') 4,(ix(j,i)-1,j=1,4)
   25    continue
      endif
c ................................................................
      if (ntetra4(1) .gt. 0) then
          do 26 i = ntetra4(2), ntetra4(2) + ntetra4(1)-1
            write(nplot,'(10i10)') 4,(ix(j,i)-1,j=1,4)
   26    continue
      endif
c ................................................................
      if (nhexa8(1) .gt. 0) then
          do 27 i = nhexa8(2), nhexa8(2) + nhexa8(1)-1
            write(nplot,'(10i10)') 8,ix(5,i)-1,ix(6,i)-1,ix(7,i)-1,
     .    ix(8,i)-1,ix(1,i)-1,ix(2,i)-1,ix(3,i)-1,ix(4,i)-1
   27    continue
      endif
c ................................................................
      if (nprism6(1) .gt. 0) then
          do 228 i = nprism6(2), nprism6(2) + nprism6(1)-1
            write(nplot,'(10i10)') 6,ix(4,i)-1,ix(5,i)-1,ix(6,i)-1,
     .    ix(1,i)-1,ix(2,i)-1,ix(3,i)-1
  228    continue
      endif
c ................................................................
c      if (nhexa8(1) .gt. 0) then
c          do i = nhexa8(2), nhexa8(2) + nhexa8(1)-1
c            write(nplot,'(10i10)') 20,ix(5,i)-1,ix(6,i)-1,ix(7,i)-1,
c     .       ix(8,i)-1,ix(1,i)-1,ix(2,i)-1,ix(3,i)-1,ix(4,i)-1,
c     .       ix(17,i)-1,ix(18,i)-1,ix(19,i)-1,ix(20,i)-1,
c     .       ix(9,i)-1,ix(10,i)-1,ix(11,i)-1,ix(12,i)-1,
c     .       ix(13,i)-1,ix(14,i)-1,ix(15,i)-1,ix(16,i)-1
c          enddo
c      endif
c ................................................................
c =================================================================            
c Escrevendo os tipos dos elementos
c =================================================================
      write(nplot,'(a, i10)') 'CELL_TYPES ', numel
c ................................................................
      if (nbar2(1) .gt. 0) then
          do i = nbar2(2), nbar2(2) + nbar2(1)-1
             write(nplot,'(a)') '3'
          enddo
      endif
c .................................................................
      if (ntria3(1) .gt. 0) then
          do i = ntria3(2), ntria3(2) + ntria3(1)-1
            write(nplot,'(a)') '5'
          enddo
      endif
c ................................................................
      if (nquad4(1) .gt. 0) then
          do i = nquad4(2), nquad4(2) + nquad4(1)-1
            write(nplot,'(a)') '9'
          enddo
      endif
c ................................................................
      if (ntetra4(1) .gt. 0) then
          do i = ntetra4(2), ntetra4(2) + ntetra4(1)-1
             write(nplot,'(a)') '10'
          enddo
      endif
c ................................................................
      if (nhexa8(1) .gt. 0) then
          do i = nhexa8(2), nhexa8(2) + nhexa8(1)-1
             write(nplot,'(a)') '12'
          enddo
       endif
c ................................................................
      if (nprism6(1) .gt. 0) then
          do i = nprism6(2), nprism6(2) + nprism6(1)-1
            write(nplot,'(a)') '13'
          enddo
      endif
c ................................................................
c      if (nhexa8(1) .gt. 0) then
c          do i = nhexa8(2), nhexa8(2) + nhexa8(1)-1
c             write(nplot,'(a)') '25'
c          enddo
c      endif
c ................................................................            
c =================================================================            
c Escrevendo os materiais
c =================================================================            
      write(nplot,'(a, i10)') 'CELL_DATA ', numel
      write(nplot,'(a)') 'SCALARS mesh_mat int'
      write(nplot,'(a)') 'LOOKUP_TABLE default'
c ................................................................
      if (nbar2(1) .gt. 0) then
          do i = nbar2(2), nbar2(2) + nbar2(1)-1
             write(nplot,'(i9)') ix(nen+1,i)
          enddo
      endif
c .................................................................
      if (ntria3(1) .gt. 0) then
          do i = ntria3(2), ntria3(2) + ntria3(1)-1
c            write(nplot,'(i5)') ix(4,i)
            write(nplot,'(i9)') ix(nen+1,i)
          enddo
      endif
c ................................................................
      if (nquad4(1) .gt. 0) then
          do i = nquad4(2), nquad4(2) + nquad4(1)-1
            write(nplot,'(i9)') ix(nen+1,i)
          enddo
      endif
c ................................................................
      if (ntetra4(1) .gt. 0) then
          do i = ntetra4(2), ntetra4(2) + ntetra4(1)-1
             write(nplot,'(i9)') ix(nen+1,i)
          enddo
      endif
c ................................................................
      if (nhexa8(1) .gt. 0) then
          do i = nhexa8(2), nhexa8(2) + nhexa8(1)-1
             write(nplot,'(i9)') ix(nen+1,i)
          enddo
      endif
c ................................................................
      if (nprism6(1) .gt. 0) then
          do i = nprism6(2), nprism6(2) + nprism6(1)-1
             write(nplot,'(i9)') ix(nen+1,i)
          enddo
      endif
c ................................................................
c      if (nhexa8(1) .gt. 0) then
c         do i = nhexa8(2), nhexa8(2) + nhexa8(1)-1
c             write(nplot,'(i5)') ix(21,i)
c          enddo
c      endif
      if(code .gt. 0) goto 41
c ...................................................................... 
      return 
      end
c *********************************************************************
c
c *********************************************************************
c * READPNODE : le o arquivo auxiliar com os nos que terao alguma     *
c * de suas grandeza impressas no tempo                               *
c * ----------------------------------------------------------------- *
c * parametros de entrada :                                           *
c * ----------------------------------------------------------------- *
c * fname   - nome do arquivo de entrada                              *
c * i_no    -                                                         *
c * i_nfile -                                                         *
c * num_node-                                                         *
c * nout    - arquivo de entrada                                      *
c * ----------------------------------------------------------------- *
c * parametros de saida                                               *
c * ----------------------------------------------------------------- *
c * i_no    -ponteiros para os nos                                    *
c * i_nfile -ponteiros para os arquivos                               *
c * num_node-numero de nos                                            *
c * flag    -verifica sucesso na abertura do arquivo                  *
c * ----------------------------------------------------------------- *
c *********************************************************************
      subroutine readpnode(fname,i_no,i_nfile,num_node,flag,nout)
      use Malloc
      implicit none
      include 'string.fi'
      integer*8 i_no,i_nfile
      integer num_node,i,j,no
      integer nout
      logical flag
      character*80 fname
      character*30 string
c .....................................................................
c
c ... 
      open(nout,file=fname,status= 'old' , err=1000 , action='read')
      flag = .true.
c ... numero total de nos      
      call readmacro(nout,.true.)
      write(string,'(30a)') (word(i),i=1,30)
      read(string,*) num_node
      print*,num_node,' print set nodes'
c .....................................................................
c
c ... memoria para os nos
      i_no        = alloc_4('no      ',1,num_node)
c ... memoria para os arquivos      
      i_nfile     = alloc_4('nnodew  ',1,num_node)
c ... lendo nos
      do i = 1, num_node
         call readmacro(nout,.true.)
         write(string,'(30a)') (word(j),j=1,30)
         read(string,*) no
         ia(i_no+i-1) = no
         ia(i_nfile+i-1) = 50 + (i-1)
      enddo
c .....................................................................
      close(nout)
      return
 1000 continue
      print*, '*** Erro na abertura de arquivos: ',trim(fname)
      flag = .false.
      end
c *********************************************************************
c
c *********************************************************************
c *  PRINTNODE: imprime uma grandeza do no pelo tempo                 *
c *  ---------------------------------------------------------------- *
c *  Parametro de Entrada :                                           *
c *  ---------------------------------------------------------------- *
c *  u       - vetor com a grandeza a ser escrita                     *
c *  no      - numero do no                                           *
c *  ndf     - graus de liberdade da grandeza                         *
c *  istep   - passo de tempo                                         *
c *  dt      - intercalo de tempo                                     *
c *  nameres - nome da gradeza a ser escrita                          *
c *  prename - nome do arquivo de saida                               *
c *  nout    - numero do arquivo a ser escrito                        *
c *  open    - abertura de um novo arquivo .true.                     *
c *  ---------------------------------------------------------------- *
c *  Parametro de Saida :                                             *
c *********************************************************************
      subroutine printnode(u,no,ndf,istep,dt,nameres,prename,nout,code
     .                    ,open)
      implicit none
c ===        
      character*30 nameres
      real*8 u(ndf,*)
      real*8 dt
      integer no,istep,nout,ndf,i,code
      character*80 fileout,prename,name
      logical open
c =====================================================================        
c
c ===
c ... abre um novo arquivo
      if(open) then
        fileout = name(prename,no,code)
        open(unit=nout,file=fileout)
        write(nout,'(a,a,a,i9)')
     .  '# Variacao ',trim(nameres),' no tempo no ponto',no
      else
c ... reabre o arquivo e add uma nova linha      
        fileout = name(prename,no,code)
        open(unit=nout,file=fileout,status ='old',access='append') 
      endif
c =====================================================================        
c
c === 
      write(nout,'(i9,f30.6,9f30.10)')istep,istep*dt,(u(i,no), i=1, ndf)
c      write(nout,*)istep,',',istep*dt,',',u(no)
      close(nout)
      return
      end
c =====================================================================
c
c *********************************************************************
c * Data de criacao    : 27/03/2016                                   *
c * Data de modificaco :                                              *
c * ----------------------------------------------------------------- *      
c * HEXA_TO_TETRA : transforma uma malha de hexaedros em tetraedros   *
c * ----------------------------------------------------------------- *
c * Parametro de Entrada :                                            *
c * ----------------------------------------------------------------- *
c * ix      - conectividades                                          *
c * numel   - numero do elmentos                                      *
c * nen     - numero de nos por elemento   a                          *
c * prename - nome do arquivo de saida                                *
c * nplot   - numero da unidade de saida                              *
c *  ---------------------------------------------------------------- *
c *  Parametro de Saida :                                             *
c * ----------------------------------------------------------------- *
c * OBS:                                                              *
c * ----------------------------------------------------------------- *      
c *********************************************************************      
        subroutine hexa_to_tetra(ix,numel,nen,prename,nplot)
        implicit none
        integer ix(nen+1,*),numel,nen,nplot,k,i,n
        character*80 prename,fname,name
c ...
        fname = name(prename,0,4)
        open(nplot,file=fname)
        n = nen + 1
c ...
        write(nplot,'(a,i10)') 'tetra4 ', numel*6
        k = 1
        do i = 1, numel
         write(nplot,'(6i10)') k,ix(2,i),ix(5,i),ix(7,i),ix(6,i),ix(n,i)
         k = k + 1
         write(nplot,'(6i10)') k,ix(2,i),ix(5,i),ix(3,i),ix(7,i),ix(n,i)
         k = k + 1
         write(nplot,'(6i10)') k,ix(2,i),ix(5,i),ix(1,i),ix(3,i),ix(n,i)
         k = k + 1
         write(nplot,'(6i10)') k,ix(3,i),ix(5,i),ix(8,i),ix(7,i),ix(n,i)
         k = k + 1
         write(nplot,'(6i10)') k,ix(3,i),ix(1,i),ix(4,i),ix(8,i),ix(n,i)
         k = k + 1
         write(nplot,'(6i10)') k,ix(3,i),ix(1,i),ix(8,i),ix(5,i),ix(n,i)
         k = k + 1
       enddo  
       return
      end 
c ***********************************************************************
c      
c *********************************************************************
      subroutine p_sint(ix,x,u,nnode,numel,ndm,nen,ndf,istep,nplot,code)
c **********************************************************************
c *                                                                    *
c *   PDATA: Imprime arquivos de saida de deslocamentos sem interface  *
c *   -----                                                            *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ----------------------                                           *
c *                                                                    *
c *    ix(nen+1,numel) - conetividades nodais                          *
c *    x(ndm,nnode)    - coordenadas nodais                            *
c *    u(ndf,nnode)    - campo escalar ou vetorial                     *
c *    nnode  - numero de nos                                          *
c *    numel  - numero de elementos                                    *
c *    nbar2  - numero de barras                                       *
c *    ntria3 - numero de triangulos                                   *
c *    nquad4 - numero de quadrilateros                                *
c *    ntetra4 - numero de tetraedros                                  *
c *    npenta5 - numero de elementos cunha                             *
c *    nhexa8  - numero de hexaedros                                   *
c *    ndm - dimensao                                                  *
c *    nen - numero maximo de nos por elementos                        *
c *    ndf - dimensao de u                                             *
c *    nplot - arquivo de saida                                        *
c *    cont1 - controle de impressao vtk                               *
c *    code - codigo de instrucao:                                     *
c *       ---  Geometria ( valores < 0 )  ---                          *
c *         = -1 - geometria View3D                                    *
c *         = -2 - geometria Gid                                       *
c *         = -3 - geometria ParaView                                  *
c *       ---  Solucao ( valores > 0 )  ---                            *
c *         =  1 - Temperatura View3D                                  *
c *         =  2 - Temperatura Gid                                     *
c *         =  3 - Temperatura ParaView                                *
c *         =  4 - Fluxo de Calor View3D                               *
c *         =  5 - Fluxo de Calor Gid                                  *
c *         =  6 - Fluxo de Calor ParaView                             *
c *         =  7 - Deslocamentos View3D                                *
c *         =  8 - Deslocamentos Gid                                   *
c *         =  9 - Deslocamentos ParaView                              *
c *         =  10 - Tensoes View3D                                     *
c *         =  11 - Tensoes Gid                                        *
c *         =  12 - Tensoes ParaView                                   *
c *                                                                    *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   --------------------                                             *
c *                                                                    *
c *                                                                    *
c **********************************************************************
      implicit none
      include 'elementos.fi'
      integer ix(nen+1,*),nnode,numel,ndm,nen,ndf,nplot,code,istep
c      integer ix(21,*),nnode,numel,ndm,nen,ndf,nplot,code,istep
      integer numelt
      integer i,j,k
      real*8  x(ndm,*),u(ndf,*),scale
      logical cont1
c        Comando de Geometria do Paraview:
c ======================================================================
c ======================================================================
c        Coordenadas:
c ======================================================================
      write(nplot,'(a)') '# vtk DataFile Version 3.0' 
      write(nplot,'(a)') 'Malhas de elementos finitos'          
      write(nplot,'(a)') 'ASCII' 
      write(nplot,'(a)') 'DATASET UNSTRUCTURED_GRID' 
      write(nplot,'(a,i10,a)') 'POINTS  ',nnode, ' double'
      do i = 1, nnode
         if (ndm .eq. 2) then
            write(nplot,'(3e15.5)') (x(j,i),j=1,ndm), 0.0e+00
         else if (ndm .eq. 3) then
            write(nplot,'(3e15.5)') (x(j,i),j=1,ndm)
c                write(nplot,'(7e15.5e3)') u(1,i)*scale+x(1,i),
c     .          u(2,i)*scale+x(2,i),u(3,i)*scale+x(3,i)
         endif
      enddo
c =================================================================            
c Escrevendo os nos dos elementos
c =================================================================  
c       numelt = 3*nbar2(1)   + 4*ntria3(1) + 5*nquad4(1) + 
c     .         5*ntetra4(1) + 9*nhexa8(1) + 7 *nprism6(1)
c     .         5*ntetra4(1) + 21*nhexa8(1)  
        numelt = (nen+1)*numel
      write(nplot,'(a, i10, i10)') 'CELLS ',numel, numelt
c ................................................................
       do i = 1, numel
          if (nbar2(1) .gt. 0) then
             write(nplot,'(10i10)') 2,(ix(j,i)-1,j=1,2)
          else if (ntria3(1) .gt. 0) then
             write(nplot,'(10i10)') 3,(ix(j,i)-1,j=1,3)
          else if (nquad4(1) .gt. 0) then
             write(nplot,'(10i10)') 4,(ix(j,i)-1,j=1,4)
          else if (ntetra4(1) .gt. 0) then
             write(nplot,'(10i10)') 4,(ix(j,i)-1,j=1,4)
          else if (nhexa8(1) .gt. 0) then
             write(nplot,'(10i10)') 8,ix(5,i)-1,ix(6,i)-1,ix(7,i)-1,
     .    ix(8,i)-1,ix(1,i)-1,ix(2,i)-1,ix(3,i)-1,ix(4,i)-1
          endif
       enddo
c ................................................................
c =================================================================            
c Escrevendo os tipos dos elementos
c =================================================================
      write(nplot,'(a, i10)') 'CELL_TYPES ', numel
c ................................................................
      if (nbar2(1) .gt. 0) then
          do i = 1, numel
             write(nplot,'(a)') '3'
          enddo
      endif
c .................................................................
      if (ntria3(1) .gt. 0) then
          do i = 1, numel
            write(nplot,'(a)') '5'
          enddo
      endif
c ................................................................
      if (nquad4(1) .gt. 0) then
          do i = 1, numel
            write(nplot,'(a)') '9'
          enddo
      endif
c ................................................................
      if (ntetra4(1) .gt. 0) then
          do i = 1, numel
             write(nplot,'(a)') '10'
          enddo
      endif
c ................................................................
      if (nhexa8(1) .gt. 0) then
          do i = 1, numel
             write(nplot,'(a)') '12'
          enddo
       endif
c ................................................................
      if (nprism6(1) .gt. 0) then
          do i = 1, numel
            write(nplot,'(a)') '13'
          enddo
      endif
c ................................................................
c      if (nhexa8(1) .gt. 0) then
c          do i = nhexa8(2), nhexa8(2) + nhexa8(1)-1
c             write(nplot,'(a)') '25'
c          enddo
c      endif
c ................................................................            
c =================================================================            
c Escrevendo os materiais
c =================================================================            
      write(nplot,'(a, i10)') 'CELL_DATA ', numel
      write(nplot,'(a)') 'SCALARS mesh_mat int'
      write(nplot,'(a)') 'LOOKUP_TABLE default'
c ................................................................
      if (nbar2(1) .gt. 0) then
          do i = 1, numel
             write(nplot,'(i9)') ix(nen+1,i)
          enddo
      endif
c .................................................................
      if (ntria3(1) .gt. 0) then
          do i = 1, numel
c            write(nplot,'(i5)') ix(4,i)
            write(nplot,'(i9)') ix(nen+1,i)
          enddo
      endif
c ................................................................
      if (nquad4(1) .gt. 0) then
          do i = 1, numel
            write(nplot,'(i9)') ix(nen+1,i)
          enddo
      endif
c ................................................................
      if (ntetra4(1) .gt. 0) then
          do i = 1, numel
             write(nplot,'(i9)') ix(nen+1,i)
          enddo
      endif
c ................................................................
      if (nhexa8(1) .gt. 0) then
          do i = 1, numel
             write(nplot,'(i9)') ix(nen+1,i)
          enddo
      endif
c ................................................................
      if (nprism6(1) .gt. 0) then
          do i = 1, numel
             write(nplot,'(i9)') ix(nen+1,i)
          enddo
      endif
c ................................................................
      write(nplot,'(a,i8)') 'POINT_DATA', nnode
      write(nplot,'(a,a,a)') 'VECTORS ', ' Deslocamentos ', ' double'
      do i = 1, nnode
         if (ndm .eq. 2)   then
            write(nplot,'(7e15.5e3)') (u(j,i),j=1,ndm), 0.0e+00
         elseif (ndm .eq. 3)   then
            write(nplot,'(7e15.5e3)') (u(j,i),j=1,ndm)
         endif
      enddo
      if (code .eq. 20)   then
         write(nplot,'(a)') 'SCALARS Temperaturas double'
         write(nplot,'(a)') 'LOOKUP_TABLE default'
         do i = 1, nnode
            if (ndm .eq. 2)   then
               write(nplot,'(1e15.5e3)') (u(j,i),j=1,ndf)
            elseif (ndm .eq. 3)   then
               write(nplot,'(1e15.5e3)') (u(j,i),j=1,ndf)
            endif
         enddo
      endif
      return
      end
