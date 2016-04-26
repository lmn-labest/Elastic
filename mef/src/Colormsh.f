c **********************************************************************
c *                                                                    *
c * SUBROUTINE: coloredmesh                                            *
c *                                                                    *
c * PURPOSE:    colororing a mesh in order to avoid memory conflicts   *
c *             during element assemblage phase                        *
c *                                                                    *
c **********************************************************************
      subroutine coloredmesh(ix      ,nnode ,nnodev 
     .                      ,numel   ,nenv  ,nen ,numcolors
     .                      ,i_colorg,i_elcolor)
c ----------------------------------------------------------------------      
c   Variable declarations
c ----------------------------------------------------------------------
      use Malloc
      implicit none
      include 'openmp.fi'
c 
c ... Parameters:
c 
      integer nnodev,nnode,numel,nenv,nen,numcolors
      integer ix(nen+1,*)
c 
c ... Local variables:
c 
      integer i,j,k,maxgrade
c ... Pointers
      integer*8 i_colorg,i_elcolor
      integer*8 i_icolor,i_nincid,i_incid
c ----------------------------------------------------------------------
c  Implementation of nodegrade
c ----------------------------------------------------------------------
      if (omp_elmt) then
c ...  Multicore finite element assembling:
        i_icolor = alloc_4('icolor  ',1,numel) 
        i_nincid = alloc_4('nincid  ',1,nnodev) 
c ... Compute the maxgrade of the mesh and element incidences:
        call nodegrade(ix,nnodev,numel,nenv,nen,ia(i_nincid),maxgrade) 
        i_incid  = alloc_4('incid   ',maxgrade,nnode)
        call elmincid(ix,ia(i_incid),ia(i_nincid),nnodev,numel,nenv,nen,
     .                maxgrade)
c ... Mesh coloring:
        call colormsh(ix,ia(i_nincid),ia(i_incid),ia(i_icolor),nnodev,
     .                numel,nenv,nen,maxgrade,numcolors)
        print *, 'Number of colors used : ', numcolors
        i_incid   = dealloc('incid   ')
        i_nincid  = dealloc('nincid  ')
        i_elcolor = alloc_4('elcolor ',1,numel)
        i_colorg  = alloc_4('colorg  ',2,numcolors)
      else
c ... Sequential finite element assembling:
        numcolors = 1
        i_icolor  = 1
        i_elcolor = alloc_4('elcolor ',1,numel)
        i_colorg  = alloc_4('colorg  ',2,numcolors)
      endif
      call setcolor(ix,ia(i_icolor),ia(i_elcolor),ia(i_colorg),numel,
     .              nenv,nen,numcolors)
      if (omp_elmt) i_icolor  = dealloc('icolor  ') 
      i_elcolor = locate('elcolor ')
      i_colorg  = locate('colorg  ')
      return
      end
c **********************************************************************
c *                                                                    *
c * SUBROUTINE: nodegrade                                              *
c *                                                                    *
c * PURPOSE:    determine the maxgrade of a given mesh                 *
c *                                                                    *
c **********************************************************************
      subroutine nodegrade(ix,nnodev,numel,nenv,nen,nincid,maxgrade)
c ----------------------------------------------------------------------     
c   Variable declarations
c ----------------------------------------------------------------------
      implicit none
c 
c ... Parameters:
c 
      integer nnodev,numel,nenv,nen,maxgrade
      integer ix(nen+1,*),nincid(*)
c
c ... Local variables:
c
      integer i,j,k,grade,node
c ----------------------------------------------------------------------
c  Implementation of nodegrade
c ----------------------------------------------------------------------
      maxgrade = 0
      grade    = 0
      do i = 1, nnodev
        nincid(i) = 0
      enddo
c ... Count the grade of a node:
      do i = 1, numel
        do j = 1, nenv
          node = ix(j,i)
          if(node .gt. 0)then
            nincid(node) = nincid(node) + 1 
            grade = nincid(node)
c         print *, 'grade - node : ', grade, node 
c ... Checking if the computed grade is the highest grade:
            if (grade .gt. maxgrade) maxgrade = grade
          endif  
        enddo
      enddo
      return
      end
c **********************************************************************
c *                                                                    *
c * SUBROUTINE: elmincid                                               *
c *                                                                    *
c * PURPOSE   : determines all incident elements of a node             *
c *                                                                    *
c **********************************************************************
      subroutine elmincid(ix,incid,nincid,nnodev,numel,nenv,nen
     .                   ,maxgrade)
c ---------------------------------------------------------------------
c  Variables
c ---------------------------------------------------------------------
      implicit none
c
c ... Parameters:
c
      integer nnodev,numel,nenv,nen,maxgrade
      integer ix(nen+1,*),incid(maxgrade,*),nincid(*)
c
c ... Local variables:
c
      integer i,j,k,node,ipos
c ----------------------------------------------------------------------
c  Implementation of elmincid
c ----------------------------------------------------------------------
c ... Initializing all elements of nincid:
      do i = 1, nnodev
        nincid(i) = 0
      enddo
c ----------------------------------------------------------------------
c  Computation of nodes incidences
c -----------------------------------------------------------------------
      do i = 1, numel
        do j = 1, nenv
          node = ix(j,i)
          if(node .gt. 0) then
            nincid(node) = nincid(node) + 1  
            ipos = nincid(node)
c ... Array incid stores the element incidences for each node: 
            incid(ipos,node) = i
          endif  
        enddo
      enddo
      return
      end
c **********************************************************************
c *                                                                    *
c * SUBROUTINE: colormesh                                              *
c *                                                                    *
c * PURPOSE   :                                                        *
c *                                                                    *
c **********************************************************************
      subroutine colormsh(ix,nincid,incid,icolor,nnodev,numel,nenv,nen 
     .                   ,maxgrade,numcolors)
c ----------------------------------------------------------------------
c  Variables
c ----------------------------------------------------------------------
      use Malloc
      implicit none
c
c ... Parameters:
c
      integer nnodev,numel,nenv,nen,maxgrade,numcolors
      integer ix(nen+1,*),incid(maxgrade,*),nincid(*),icolor(*)
c
c ... Local variables
c
      integer i,j,k,l,m,n,node,maxcolors,ielk,ielfirst,ixj,ik
      parameter (maxcolors = 200)
      logical found, quit
c ----------------------------------------------------------------------
c  Implementation of colormesh
c ----------------------------------------------------------------------
c ... Initializing all elements of nincid:
      do i = 1, numel
        icolor(i) = 0
      enddo
c ----------------------------------------------------------------------
c  Coloring loop (icolor(i) = color of ith element)
c ----------------------------------------------------------------------
      quit = .false.
      do i = 1, maxcolors
        do node = 1, nnodev
          found = .false.
          l = nincid(node)
c ----------------------------------------------------------------------
          do m = 1, l
            ik = incid(m,node)
            if (icolor(ik) .eq. i) found = .true.
          enddo
c ----------------------------------------------------------------------
          if (found .eqv. .false. ) then
            ielfirst = 0
            do m = 1, l
              ik = incid(m,node)
              if (icolor(ik) .eq. 0) ielfirst = ik
            enddo
              if (ielfirst .ne. 0) then
                icolor(ielfirst) = i
c ----------------------------------------------------------------------
c ... must check !
                do j = 1, nenv
                  ixj = ix(j,ielfirst)
                  if (ixj .gt. 0) then
                    do k = 1, nincid(ixj)
                      ik = incid(k,ixj)
                      ielk = ik
                      if (icolor(ielk) .eq. 0) icolor(ielk) = -1
                    enddo
                  endif
                enddo
c ----------------------------------------------------------------------
              endif
            endif
          enddo
        quit = .true. 
        do n = 1, numel
          if (icolor(n) .eq. -1) then
            icolor(n) = 0
            quit = .false.
          endif
        enddo
        numcolors = i
        if (quit) goto 100
      enddo
      print *, 'Numero maximo de cores excedido ', maxcolors
      stop
 100  continue
      return
      end
c **********************************************************************
c *                                                                    *
c * SUBROUTINE: setcolor                                               *
c *                                                                    *
c * PURPOSE   :                                                        *
c *                                                                    *
c **********************************************************************
      subroutine setcolor(ix,icolor,ielcolor,icolorg,numel,nenv,nen,
     .                    numcolors)
c ----------------------------------------------------------------------      
c  Variable declarations
c ----------------------------------------------------------------------
      implicit none
c
c ... Parameters:
c
      integer numel,nenv,nen,numcolors
      integer ix(nen+1,*),icolor(*),ielcolor(*),icolorg(2,*)
c
c ... Local variables:
c
      integer i,j,ipos,icount
      logical first
c ----------------------------------------------------------------------
c ... Implementation of setcolor
c ----------------------------------------------------------------------
      ipos = 1
      first = .true.
c     do i = 1, numel
c       ix(nen+1,i) = icolor(i)
c     enddo
      if (numcolors .gt. 1) then
        do i = 1, numcolors
          icount = 0
          do j = 1, numel
            if (icolor(j) .eq. i) then
              ielcolor(ipos) = j
              if (first) icolorg(1,i) = ipos
              first = .false.
              ipos = ipos + 1 
              icount = icount + 1
            endif
          enddo
          icolorg(2,i) = icolorg(1,i) + icount - 1
          first = .true.
        enddo
      else
        do i = 1, numel
          ielcolor(i) = i
        enddo
        icolorg(1,1) = 1
        icolorg(2,1) = numel
      endif
c     do i = 1, numel
c       print *, ' ielcolor ', ielcolor(i)
c     enddo
c     do i = 1, numcolors
c       print *, ' icolorg ', icolorg(i)
c     enddo
      return
      end
