c **********************************************************************
c * Data de criacao    : 10/04/2016                                    *
c * Data de modificaco : 00/00/0000                                    * 
c * ------------------------------------------------------------------ *
c * PRE_DIAG: precondicionador diagonal                                *
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ *
c * m   - indefinido                                                   *
c * ad  - coeficientes da diagonal principal                           *
c * neq - numero de equacoes                                           *
c * pd  - .true.  modulo da diagonal                                   *
c *       .false. valor com sinal da diagonal                          *
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ * 
c * m   - precondicionador diagonal ( M-1)                             *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ * 
c **********************************************************************
      subroutine pre_diag(m,ad,neq,pd)
      implicit none
      real*8 m(*),ad(*)
      integer i,neq
      logical pd
c ... positva definida
      if(pd) then
        do i = 1, neq
          m(i) = dabs(1.d0/ad(i))
        enddo
c ....................................................................
c
c ...
      else
        do i = 1, neq
          m(i) = 1.d0/ad(i)
        enddo
      endif
c ....................................................................
c
c ...
      return
      end
c **********************************************************************
c
c **********************************************************************
c * Data de criacao    : 10/04/2016                                    *
c * Data de modificaco : 00/00/0000                                    * 
c * ------------------------------------------------------------------ *
c * PRE_DIAG_SCHUR : precondicionador diagonal com complemento schur   *
c * aproximado ( problema poro mecanico)
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ *
c * ia  - csrc                                                         *
c * ja  - csrc                                                         *
c * m   - indefinido                                                   *
c * ad  - coeficientes da diagonal principal                           *
c * al  - coeficientes da parte inferior da matriz                     *
c * neq - numero de equacoes                                           *
c * nequ- numero de equacoes no bloco Kuu                              *
c * pd  - .true.  modulo da diagonal                                   *
c *       .false. valor com sinal da diagonal                          *
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ * 
c * m   - precondicionador diagonal ( M-1)                             *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ * 
c *     | diag(Kuu)                0                         |         *
c * M = |                                                    |         *  
c *     |    0         diag(Kpp + Kpu*(diag(Kuu))(-1)*(Kpu)T |         *
c *                                                                    *                
c *     | Kuu              (Kpu)T   |                                  *
c * K = |                           |                                  *  
c *     | Kpu              -Kpp     |                                  *          
c **********************************************************************
      subroutine pre_diag_schur(ia,ja,m,ad,al,w,neq,nequ)
      implicit none
      real*8 m(*),ad(*),al(*),w(*),tmp
      integer i,j,neq,nequ,jak,ia(*),ja(*)
c ... positva definida
      do i = 1, nequ
        m(i) = 1.d0/ad(i)
      enddo
c ....................................................................
c
c ...
      do i = nequ + 1, neq
        w(1:nequ) = 0.d0
        do j = ia(i), ia(i+1) - 1
          jak = ja(j)
          if( jak .le. nequ) w(jak) = al(j)
        enddo
c ... C + BD(-1)BT
        tmp = dabs(ad(i))
        do j = 1, nequ
          tmp = tmp + w(j)*w(j)*m(j)
        enddo
        m(i) = 1.0d0/tmp
      enddo
c ....................................................................
c
c ...
      return
      end
c **********************************************************************
c
c **********************************************************************
c * Data de criacao    : 15/06/2016                                    *
c * Data de modificaco : 00/00/0000                                    * 
c * ------------------------------------------------------------------ *
c * BLOCK_PRECOND : precondicionador bloco diagonal simetrico          *
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ *
c * ad     - coeficientes da diagonal principal                        *
c * al     - coeficientes da parte inferior da matriz                  *
c * ia     - csrc                                                      *
c * ja     - csrc                                                      *
c * m      - indefinido                                                *
c * n      - numero de linhas                                          *
c * aux(100)- arranjo auxiliar de trabalho                             *
c * iparam - parametros do bloco diagonal                              *
c *        - iparam(1) - nao definido                                  *
c *        - iparam(2) - nao definido                                  *
c *        - iparam(3) - numero de termos nos bloco                    *
c *        - iparam(4) - tamanho do bloco                              *
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ * 
c * m      - precondicionador bloco diagonal (M-1)                     *
c * iparam - parametros do bloco diagonal                              *
c *        - iparam(1) - numero de sub matriz em blocos                *
c *        - iparam(2) -  numero de inversos da diagonal simples       *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ * 
c * para tamanho de blocos 2 e 3 inversao direta M(-1)                 *        
c * para tamanho de blocos maior ou igual a 4 fatoracao LLt            *                                       c *                                                                    *
c *                                                                    *
c * bloco (2x2) - iparam(3) = 2                                        *
c *    | a11 a21 |        | a22 a21 |                                  *
c *    |         | =1/det |         |                                  *
c *    | a21 a22 |        |-a21 a11 |                                  *
c * m(3,*) = 1/det | a11 a22 -a21 |                                    *
c *                                                                    *
c * bloco (3x3) - iparam(3) = 3                                        *
c *    | a11 a21 a31 |        | c11 c12 c31|                           *
c *    | a21 a22 a32 | =1/det | c21 c22 c32|                           *  
c *    | a31 a32 a33 |        | c31 c32 c33|                           *
c * m(6,*) = 1/det | c11 c22 c33 c21 c31 c32 |                         *
c *                                                                    *
c * bloco (4x4) fatorado(LLt)                                          *
c *          | a11 a21 a31 a41 |                                       *
c *          | a21 a22 a32 a42 |                                       * 
c *          | a31 d32 a33 a43 |                                       *
c *          | a41 a42 d43 a44 |                                       *
c * m(10,1) = | a11 a22 a33 a44 a21 a31 a32 a41 a42 a43 |              *   
c **********************************************************************
      subroutine block_precond(ad,al,ia,ja,m,n,aux,param)
      implicit none
      include 'precond.fi'
      integer nblock,n,ndiv,rest,tmp1,tmp2,i,j,ii,jj,param(*)
      integer ia(*),ja(*),kk,k,l,c
      real*8 ad(*),al(*),m(param(3),*),idet
      real*8 a11,a21,a22,a33,a31,a32
      real*8 aux(param(4),*),v(max_block),w(max_block)
      logical bspd
c .....................................................................
c
c ...
      bspd = .true.
c .....................................................................
c
c ...
      nblock = param(4) 
c .....................................................................
c 
c ... numero de blocos completos
      ndiv = n/nblock
c .....................................................................
c
c ... 
      rest = mod(n,nblock)
c .....................................................................
c
c ...
      param(1) = ndiv
      param(2) = rest
c .....................................................................    
c
c .....................................................................
      do i = 1, ndiv + rest
        do j = 1, param(3)
          m(j,i) = 0.d0
        enddo
      enddo
c .....................................................................
c
c ... blocos inverso(2x2)
c      | a11 a21 |        | a22 -a21|
c      |         | =1/det |         | 
c      | a21 a22 |        |-a21  a11| 
      if( nblock .eq. 2) then
        do i = 1, ndiv
          ii = 2*i - 1
c ... 
          a11 = ad(ii)
          a22 = ad(ii+1)
c ......................................................................
c
c ...
          a21 = 0.d0
          do j = ia(ii+1) , ia(ii+2) - 1
            jj = ja(j)
            if( jj .eq. ii) a21 = al(j)
          enddo
c ......................................................................
c
c ...        
          idet = 1.d0/(a11*a22-a21*a21)
c ...
          m(1,i) =   idet*a22
          m(2,i) =   idet*a11
          m(3,i) =  -idet*a21
        enddo
c .....................................................................    
c
c ... blocos inverso(3x3)
c      | a11 a21 a31 |        | c11 c12 c31|
c      | a21 a22 a32 | =1/det | c21 c22 c32| 
c      | a31 a32 a33 |        | c31 c32 c33| 
      else if (nblock .eq. 3) then
        do i = 1, ndiv
          ii = 3*i - 2
c ... 
          a11 = ad(ii)
          a22 = ad(ii+1)
          a33 = ad(ii+2)
c ......................................................................
c
c ...
          a21 = 0.d0
          do j = ia(ii+1) , ia(ii+2) - 1
            jj = ja(j)
            if( jj .eq. ii) a21 = al(j)
          enddo
c ......................................................................
c
c ...
          a31 = 0.d0
          a32 = 0.d0
          do j = ia(ii+2) , ia(ii+3) - 1
            jj = ja(j)
            if( jj .eq. ii  ) then
              a31 =al(j)
            else if( jj .eq. ii+1) then
              a32 =al(j)
            endif
          enddo
c .....................................................................
c
c ...
          call inverse_matrix_3x3_sym(m,a11,a22,a33,a21,a31,a32)
c .....................................................................
c
c ...
          m(1,i) = a11
          m(2,i) = a22
          m(3,i) = a33
          m(4,i) = a21
          m(5,i) = a31
          m(6,i) = a32
c .....................................................................
        enddo 
c .....................................................................
c
c ... blocos LLT(4x4,5x5,6x6,...)
      else if (nblock .ge. 4) then
c ...
        do i = 1, ndiv
c ...
          do j = 1, nblock
            do k = 1, nblock
              aux(k,j) = 0.d0
            enddo
          enddo
c ......................................................................
c 
c ...
          ii = (i-1)*nblock + 1
c ... extraindo o bloco da csrc
          do j = 1, nblock 
            aux(j,j) = ad(ii+j-1)
          enddo
c
          do l = 1, nblock - 1
            kk = ii + l
            do j = ia(kk) , ia(kk+1) - 1
              jj = ja(j) 
              do c = 0, l - 1
                 if( jj .eq. ii + c) then
                   aux(l+1,c+1) = al(j)
                endif
              enddo
            enddo
          enddo
c .....................................................................
c
c ... fat LLT
          if(bspd) then
            call choleskyc(nblock,aux)
c ... fat LDLt
          else
            call ldltc(nblock,aux,v,w)
          endif
c .....................................................................
c
c ... armazenando L
          do j = 1, nblock
            m(j,i) = aux(j,j) 
          enddo
c
          kk = nblock
          do j = 2, nblock
            do k = 1, j - 1
              kk = kk + 1
              m(kk,i) = aux(j,k)
            enddo
          enddo
c .....................................................................

c .....................................................................
        enddo 
c .....................................................................     
      endif
c .....................................................................    
c
c ... 
      tmp1 = ndiv*nblock
      tmp2 = ndiv
      do i = 1, rest 
        ii      = tmp1 + i
        jj      = tmp2 + i
        m(1,jj) = 1.d0/ad(ii)
      enddo
c .....................................................................    
c
c ...
c     print*,ndiv,rest
c     do i = 1, ndiv+rest
c       write(*,'(6es9.2)')m(1:iparam(3),i)
c     enddo
      return
      end
c *********************************************************************
c
c **********************************************************************
c * Data de criacao    : 18/06/2016                                    *
c * Data de modificaco : 00/00/0000                                    *
c * ------------------------------------------------------------------ *
c * op_block_precond : operador do precondicionador bloco jacobi       *
c * diagonal                                                           *
c * ------------------------------------------------------------------ *
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ *
c * m       - matriz bloco diagonal                                    *
c * x       - vetor x                                                  *
c * y       - nao definido                                             *
c * aux(100)- arranjo auxiliar de trabalho                             *
c * iparam - parametros do bloco diagonal                              *
c *        - iparam(1) - numero de sub matriz em blocos                *
c *        - iparam(2) - numero de inversos da diagonal simples        *
c *        - iparam(3) - numero de termos nos bloco                    *
c *        - iparam(4) - tamanho do bloco                              *
c * ------------------------------------------------------------------ *
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ *
c * y - M(-1) * x                                                      *
c * ------------------------------------------------------------------ *
c * OBS:                                                               *
c * ------------------------------------------------------------------ *
c * para tamanho de blocos 2 e 3 inversao direta (y=M(-1)x             *        
c * para tamanho de blocos maior ou igual a 4 solucao do sistema por   *
c * metodo direto ( (LLt)y=x )                                         *
c *                                                                    *                                                  
c * bloco (2x2) invertido                                              *
c *          | a11 a21              |                                  *
c *          | a21 a22              |                                  * 
c * inv(M) = |         a33 a34      |                                  *
c *          |         a43 a44      |                                  *
c *          |                  a55 |                                  *
c * m(3,1) = | a11 a22 a21 |                                           *
c * m(3,2) = | a33 a44 a34 |                                           *
c * m(3,3) = | a55 0.0 0.0 |                                           *
c *                                                                    *
c * bloco (3x3) invertido                                              *
c *          | a11 a21 a31          |                                  *
c *          | a21 a22 a32          |                                  * 
c * inv(M) = | a31 d32 a33          |                                  *
c *          |             a44      |                                  *
c *          |                  a55 |                                  *
c * m(6,1) = | a11 a22 a33 a21 a31 a32 |                               *
c * m(6,2) = | a44 0.0 0.0 0.0 0.0 0.0 |                               *
c * m(6,3) = | a55 0.0 0.0 0.0 0.0 0.0 |                               *
c *                                                                    *
c * bloco (4x4) fatorado(LLt)                                          *
c *          | a11 a21 a31 a41      |                                  *
c *          | a21 a22 a32 a42      |                                  * 
c * LLt    = | a31 d32 a33 a43      |                                  *
c *          | a41 a42 d43 a44      |                                  *
c *          |                  a55 |                                  *
c * m(10,1) = | a11 a22 a33 a44 a21 a31 a32 a41 a42 a43 |              *   
c * m(10,2) = | a44 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 |              *   
c **********************************************************************
      subroutine op_block_precond(m,x,y,aux,iparam)
      implicit none
      include 'mpif.h'
      include 'time.fi'
      integer i,ii,j,jj,k,kk,iparam(*),tmp1,tmp2,ndiv,nblock,rest
      real*8 m(iparam(3),*),x(*),y(*),a1,a2,a3,a4,a5,a6,x1,x2,x3
      real*8 aux(iparam(4),*)
      logical bspd
c ...
      bspd = .true.
c .....................................................................
c
c ...
      prebdiagtime = Mpi_Wtime() - prebdiagtime
c .....................................................................
c
c
c ...
      nblock = iparam(4)
      ndiv   = iparam(1) 
      rest   = iparam(2)
c .....................................................................
c 
c ... bloco diagonal
      if (nblock .eq. 2) then  
        do i = 1, ndiv
          ii      = i*2 - 1
c
          a1      = m(1,i)
          a2      = m(2,i)
          a3      = m(3,i)
c
          x1      = x(ii)
          x2      = x(ii+1)
c ... y = M(-1)x
          y(ii)   = a1*x1 + a3*x2
          y(ii+1) = a3*x1 + a2*x2
        enddo
c .....................................................................
c
c ...
      else if( nblock .eq. 3) then
        do i = 1, ndiv
          ii      = i*3 - 2
c ... a11
          a1      = m(1,i)
c ... a22
          a2      = m(2,i)
c ... a33
          a3      = m(3,i)
c ... a21
          a4      = m(4,i)
c ... a31
          a5      = m(5,i)
c ... a32
          a6      = m(6,i)
c ...
          x1      = x(ii)
          x2      = x(ii+1)
          x3      = x(ii+2)
c ... y = M(-1)x
          y(ii)   = a1*x1 + a4*x2 + a5*x3
          y(ii+1) = a4*x1 + a2*x2 + a6*x3
          y(ii+2) = a5*x1 + a6*x2 + a3*x3
        enddo  
c .....................................................................
c
c ...
      else if( nblock .ge. 4) then
        do i = 1, ndiv
          ii = (i-1)*nblock + 1
c ...    
          do j = 1, nblock
            aux(j,j) = m(j,i)
          enddo
c .....................................................................
c
c ...
          kk = nblock
          do j = 2, nblock
            do k = 1, j - 1
              kk = kk + 1
              aux(j,k) =  m(kk,i)
              aux(k,j) = aux(j,k)  
            enddo
          enddo
c .....................................................................
c
c ... LLTy = x
          if(bspd) then
            call solv_cholesky(aux,x(ii),y(ii),nblock)
c ... LDLTy = x
          else
            call solv_ldlt(aux,x(ii),y(ii),nblock)
          endif
c .....................................................................
        enddo 
      endif
c .....................................................................
c
c ... diagonal simples
      tmp1 = ndiv*nblock
      tmp2 = ndiv
      do i = 1, rest
        ii      = tmp1 + i
        jj      = tmp2 + i  
c
        a1      = m(1,jj)
c
        y(ii)   = a1*x(ii) 
      enddo
c .....................................................................
c
c ...
      prebdiagtime = Mpi_Wtime() - prebdiagtime
c .....................................................................
c ... 
c     do i = 1, 13176
c       print*,i,x(i),y(i)
c     enddo
      return  
      end
c **********************************************************************  
c
c **********************************************************************
       subroutine ildlt1(n,ia,ja,al,ad,ldlt,v,shift,fcheck)  
c **********************************************************************
c * Data de criacao    : 10/04/2016                                    *
c * Data de modificaco : 00/00/0000                                    * 
c * ------------------------------------------------------------------ *   
c * ILDLT1:fatora imcompleta LDLt com matriz simetrica A no formato    *
c * CSRC.                                                              *
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ *
c * n - numero de equacoes                                             *
c * ia(1:n+1) - ponteiro das linhas, para o formato CSR                *
c * ja(1:nad) - ponteiro das colunas no formato CSR                    *
c * al(1:nad) - parte triangular inferior de A                         *
c * ad(1:n)   - diagonal da matriz A                                   *
c * ldlt(1:n+nad) - nao definido                                       *
c * v(1:n)    - arranjo auxiliar (nao inicializado)                    *
c * shift     - parametro de deslocamento para manter a fotaracao      * 
c *             estavel ( shift >= 0)                                  *   
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ * 
c * ldlt(n+1:n+nad) - fator L                                          *
c * ldlt(1:n)       - inverso da diagonal fatorada                     *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ * 
c * LDLT: fatora a matriz simetrica A na forma LDL. Os coeficientes    *
c * Aij (j<i) sao armazenados em a(1:nad), no formato CSR,             *
c * e a diagonal eh armazenada em ad(1:n). Os fatores L e D sao        *
c * armazenados em ldlt(n+1:n+1+nad) e ldlt(1:n), no mesmo formato CSR *
c * e com a mesma estrututa ia, ja. Pela propria definicao do CSR,     *
c * esta fatoracao eh incompleta. No entanto, pode tornar-se           *
c * completa utilizando-se a tecnica de fill-in.                       *
c * Esta implementacao eh uma adaptacao para o formato CSR do          *
c * algoritmo 5.1.2 em Matriz Computations, p. 84, Golub e Van Loan.   *
c **********************************************************************
      implicit none      
      integer n,ia(*),ja(*),i,j,k,kk,i1,i2,kj,nad
      real*8 al(*),ad(*),ldlt(*),v(*),t
      real*8 dmax,mod,shift
      logical fcheck  
      real*8 tol
      parameter (tol=1.d-14)
c ......................................................................
      nad      = ia(n+1)-1
      do i = 1, nad
        ldlt(i+n) = al(i)
      enddo
c ... A = A + shift*diag(A)
      ldlt(1:n)   = ad(1:n) + shift*ad(1:n)
c ......................................................................
      v(1:n)   = 0.d0
c ......................................................................
c
c ...
      dmax = 0.d0
      do 50 i = 1, n
        mod  = dabs(ad(i))
        dmax = max(mod,dmax)
  50  continue
c ......................................................................
c
c ......................................................................      
      do 1000 j = 1, n   
        i1 = ia(j)
        i2 = ia(j+1)-1
        do 100 i = i1, i2
          k    = ja(i)
          v(k) = ldlt(i+n)*ldlt(k)
  100   continue
c ......................................................................                    
        t = 0.d0
        do 200 i = i1, i2
          t = t + ldlt(i+n)*v(ja(i))
  200   continue
        v(j) = ldlt(j) - t
c ... pivo negativo ou muito pequeno
        if(v(j) .lt. tol*dmax .and. fcheck) then
          print*,'(iLDLt) Valor do Pivo invalido !!!',v(j)
          print*,'Linha: ',j
          call stop_mef()
        endif
        ldlt(j) = v(j)
c .....................................................................                  
        do 310 k = j+1, n
          kj = 0
          t  = 0.d0
          do 300 i = ia(k), ia(k+1)-1
            kk = ja(i)
            if (kk .lt. j) t  = t + ldlt(i+n)*v(kk)
            if (kk .eq. j) kj = i + n
  300     continue
          if (kj .gt. 0 ) ldlt(kj) = (ldlt(kj) - t)/v(j)
  310   continue
c .....................................................................
        v(j) = 0.d0  
        do 400 i = i1, i2             
          k    = ja(i) 
          v(k) = 0.d0   
  400   continue
c .....................................................................  
 1000 continue
c .....................................................................
c
c ...
      ldlt(1:n)   = 1.d0/ldlt(1:n)
c .....................................................................
c     open(15,file='ildlt.csr')
c     do i = 1, 1500
c       write(15,'(i9,1d15.2)')i,ldlt(i)
c     enddo 
c     write(15,*)'l'
c     do i = 1, nad
c       write(15,'(i9,2d15.2)')i,al(i),ldlt(i+n)
c     enddo 
c      close(15)
      return
      end
c **********************************************************************
c
c **********************************************************************
       subroutine ildlt2(n,ia,ja,al,ad,ldlt,v,shift,fcheck)  
c **********************************************************************
c * Data de criacao    : 10/04/2016                                    *
c * Data de modificaco : 20/04/2016                                    * 
c * ------------------------------------------------------------------ *   
c * ILDLT: fatora imcompleta LDLt com matriz simetrica A no formato    *
c * CSRC.                                                              *
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ *
c * n - numero de equacoes                                             *
c * ia(1:n+1) - ponteiro das linhas, para o formato CSR                *
c * ja(1:nad) - ponteiro das colunas no formato CSR                    *
c * al(1:nad) - parte triangular inferior de A                         *
c * ad(1:n)   - diagonal da matriz A                                   *
c * ldlt(1:n+nad) - nao definido                                       *
c * v(1:n)    - arranjo auxiliar (nao inicializado)                    *
c * shift     - parametro de deslocamento para manter a fotaracao      * 
c *             estavel ( shift >= 0)                                  *   
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ * 
c * ldlt(n+1:n+nad) - fator L                                          *
c * ldlt(1:n)       - inverso da diagonal fatorada                     *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ * 
c * LDLT: fatora a matriz simetrica A na forma LDL. Os coeficientes    *
c * Aij (j<i) sao armazenados em a(1:nad), no formato CSR,             *
c * e a diagonal eh armazenada em ad(1:n). Os fatores L e D sao        *
c * armazenados em ldlt(n+1:n+1+nad) e ldlt(1:n), no mesmo formato CSR *
c * e com a mesma estrututa ia, ja. Pela propria definicao do CSR,     *
c * esta fatoracao eh incompleta. No entanto, pode tornar-se           *
c * completa utilizando-se a tecnica de fill-in.                       *
c * Esta implementacao eh uma adaptacao para o formato CSR do          *
c * algoritmo 5.1.2 em Matriz Computations, p. 84, Golub e Van Loan.   *
c **********************************************************************
      implicit none      
      integer n,ia(*),ja(*),i,j,k,kk,i1,i2,kj,nad
      real*8 al(*),ad(*),ldlt(*),v(*),t
      real*8 dmax,mod,shift
      logical fcheck  
      real*8 tol
      parameter (tol=1.d-14)
c ......................................................................
      nad      = ia(n+1)-1
      do i = 1, nad
        ldlt(i+n) = al(i)
      enddo
c ... A = A + shift*diag(A)
      ldlt(1:n)   = ad(1:n) + shift*ad(1:n)
c ......................................................................
      v(1:n)   = 0.d0
c ......................................................................
c
c ...
      dmax = 0.d0
      do 50 i = 1, n
        mod  = dabs(ad(i))
        dmax = max(mod,dmax)
  50  continue
c ......................................................................
c
c ......................................................................      
      do 1000 j = 1, n   
        i1 = ia(j)
        i2 = ia(j+1)-1
        do 100 i = i1, i2
          k    = ja(i)
          v(k) = ldlt(i+n)*ldlt(k)
  100   continue
c ......................................................................                    
        t = 0.d0
        do 200 i = i1, i2
          t = t + ldlt(i+n)*v(ja(i))
  200   continue
        v(j) = ldlt(j) - t
c ... pivo negativo ou muito pequeno
        if(v(j) .lt. tol*dmax .and. fcheck) then
          print*,'(iLDLt) Valor do Pivo invalido !!!',v(j)
          print*,'Linha: ',j
          call stop_mef()
        endif
        ldlt(j) = v(j)
c .....................................................................                  
        do 310 k = j+1, n
          kj = 0
          t  = 0.d0
c ... somatorio da da coluna 1 ate j - 1
          do 300 i = ia(k), ia(k+1)-1
            kk = ja(i)
            if (kk .lt. j) then
              t = t + ldlt(i+n)*v(kk)
            else if (kk .eq. j) then
              kj = i + n
              go to 305
            else if (kk .gt. j) then
              go to 306
            endif
  300     continue
c ......................................................................
c
c ...  
  305     continue
          if(kj .ne. 0) ldlt(kj) = (ldlt(kj) - t)/v(j)
c ......................................................................
c
c ...
  306     continue
  310   continue
c .....................................................................
        v(j) = 0.d0  
        do 400 i = i1, i2             
          k    = ja(i) 
          v(k) = 0.d0   
  400   continue
c .....................................................................  
 1000 continue
c .....................................................................
c
c ...
      ldlt(1:n)   = 1.d0/ldlt(1:n)
c .....................................................................
c     open(15,file='ildlt.csr')
c     do i = 1, 1500
c       write(15,'(i9,1d15.2)')i,ldlt(i)
c     enddo 
c     write(15,*)'l'
c     do i = 1, nad
c       write(15,'(i9,2d15.2)')i,al(i),ldlt(i+n)
c     enddo 
c      close(15)
      return
      end
c **********************************************************************
c 
c **********************************************************************
       subroutine ichfat(n,ia,ja,al,ad,ich,v,shift,fcheck)  
c **********************************************************************
c * Data de criacao    : 22/04/2016                                    *
c * Data de modificaco : 00/00/0000                                    * 
c * ------------------------------------------------------------------ *   
c * ICHFAT : fatora imcompleta LLt (choleskyc) com matriz simetrica A  *
c * no formato CSRC.                                                   *
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ *
c * n - numero de equacoes                                             *
c * ia(1:n+1)    - ponteiro das linhas, para o formato CSR             *
c * ja(1:nad)    - ponteiro das colunas no formato CSR                 *
c * al(1:nad)    - parte triangular inferior de A                      *
c * ad(1:n)      - diagonal da matriz A                                *
c * ich(1:n+nad) - nao definido                                        *
c * v(1:n)       - arranjo auxiliar (nao inicializado)                 *
c * shift        - parametro de deslocamento para manter a fotaracao   * 
c *                estavel ( shift >= 0)                               *   
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ * 
c * ich(n+1:n+nad) - fator L                                           *
c * ich(1:n)       - inverso da diagonal fatorada                      *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ * 
c **********************************************************************
      implicit none      
      integer n,ia(*),ja(*),i,j,k,kk,i1,i2,kj,nad
      real*8 al(*),ad(*),ich(*),v(*),t
      real*8 shift
      logical fcheck  
      real*8 tol
      parameter (tol=1.d-14)
c ......................................................................
      nad      = ia(n+1)-1
      do i = 1, nad
        ich(i+n) = al(i)
      enddo
c ... A = A + shift*diag(A)
      ich(1:n)   = ad(1:n) + shift*ad(1:n)
c ......................................................................
      v(1:n)   = 0.d0
c ......................................................................
c
c ......................................................................      
      do 1000 j = 1, n   
        i1 = ia(j)
        i2 = ia(j+1)-1
c ... somatorio ajp (j = linha;p = 1 ... j -1)
        t = 0.d0
        do 200 i = i1, i2
          k    = ja(i)
          v(k) = ich(i+n)
          t    = t + ich(i+n)*ich(i+n)
  200   continue
c .....................................................................
        v(j) = ich(j) - t
c ... pivo negativo ou muito pequeno
        if(v(j) .le. 0.d0 .and. fcheck) then
          print*,'(iLLT) Valor do Pivo invalido !!!',v(j)
          print*,'Linha: ',j  
          call stop_mef()
        endif
        ich(j) = dsqrt(v(j)) 
c .....................................................................                  
        do 310 k = j+1, n
          kj = 0
          t  = 0.d0
c ... somatorio da da coluna 1 ate j - 1
          do 300 i = ia(k), ia(k+1)-1
            kk = ja(i)
            if (kk .lt. j) then
              t = t + ich(i+n)*v(kk)
            else if (kk .eq. j) then
              kj = i + n
              go to 305
            else if (kk .gt. j) then
              go to 306
            endif
  300     continue
c ......................................................................
c
c ...  
  305     continue
          if(kj .ne. 0) ich(kj) = (ich(kj) - t)/ich(j)
c ......................................................................
c
c ...
  306     continue
  310   continue
c .....................................................................
        v(j) = 0.d0  
        do 400 i = i1, i2             
          k    = ja(i) 
          v(k) = 0.d0
  400   continue
c .....................................................................  
 1000 continue
c .....................................................................
c
c ...
       ich(1:n)   = 1.0d0/ich(1:n)
c .....................................................................
c     open(15,file='ich.csr')
c     write(15,*)'d'
c     do i = 1, n
c       write(15,'(d15.2)')ich(i)
c     enddo 
c     write(15,*)'l'
c     do i = 1, nad 
c      do i = 1, 1500
c       write(15,'(d15.2)')ich(i+n)
c     enddo 
c     close(15)
      return
      end
c **********************************************************************
c
c **********************************************************************
      subroutine ildlt_solv(n,ia,ja,ad,a,b,x)
c **********************************************************************
c * Data de criacao    : 10/04/2016                                    *
c * Data de modificaco : 18/04/2016                                    *
c * ------------------------------------------------------------------ *  
c * ildlt_solv: resolve o sistema com a fatoracao imcompleta LDLt      *
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ * 
c * n - numero de equacoes                                             *
c * ia(1:n+1)  - ponteiro das linhas, para o formato CSR               *
c * ja(1:nad)  - ponteiro das colunas no formato CSR                   *
c * ad(1:nad)  - fatores D da matriz A ( D: 1/D )                      *
c * a(1:nad)   - fatores L da matriz A                                 *
c * b(1:n)     - vetor independente                                    *
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ * 
c * x(1:n) - vetor solucao.                                            *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ * 
c * Todas as informação da matriz são da parte inferior                * 
c **********************************************************************
      implicit none
      include 'mpif.h'
      include 'time.fi'
      integer n,ia(*),ja(*)
      integer i,k,kk
      real*8  a(*),ad(*),b(*),x(*),t
c ......................................................................
c
c ...
      ifatsolvtime = Mpi_Wtime() - ifatsolvtime
      x(1:n) = b(1:n)
c ......................................................................
c
c ... Forward substitution: Ly = b 
       do 110 i = 2, n
        t = x(i)
        do 100 k = ia(i), ia(i+1) - 1
           t = t - a(k)*x(ja(k))
  100   continue
        x(i) = t 
  110 continue
c ......................................................................
c
c ... Dz = y
      do 115 i = 1, n
        x(i) = x(i) * ad(i)
  115 continue
c ......................................................................
c
c ... Backward substitution: Ltx = z
      do 210 i = n, 1, -1 
        t = x(i)
        do 200 k = ia(i), ia(i+1)-1
           kk    = ja(k)
           x(kk) = x(kk) - t*a(k)
  200   continue
  210 continue
c ......................................................................
      ifatsolvtime = Mpi_Wtime() - ifatsolvtime    
      return
      end
c **********************************************************************
c
c **********************************************************************
      subroutine illt_solv(n,ia,ja,ad,a,b,x)
c **********************************************************************
c * Data de criacao    : 23/04/2016                                    *
c * Data de modificaco : 00/00/0000                                    *
c * ------------------------------------------------------------------ *  
c * illt_solv: resolve o sistema com a fatoracao imcompleta LLt        *
c * (cholesky)                                                         *
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ * 
c * n - numero de equacoes                                             *
c * ia(1:n+1)  - ponteiro das linhas, para o formato CSR               *
c * ja(1:nad)  - ponteiro das colunas no formato CSR                   *
c * ad(1:nad)  - fatores D da matriz A ( D: 1/D )                      *
c * a(1:nad)   - fatores L da matriz A                                 *
c * b(1:n)     - vetor independente                                    *
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ * 
c * x(1:n) - vetor solucao.                                            *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ * 
c * Todas as informação da matriz são da parte inferior                * 
c **********************************************************************
      implicit none
      include 'mpif.h'
      include 'time.fi'
      integer n,ia(*),ja(*)
      integer i,k,kk
      real*8  a(*),ad(*),b(*),x(*),t
c ......................................................................
c
c ...
      ifatsolvtime = Mpi_Wtime() - ifatsolvtime
      x(1:n) = b(1:n)
c ......................................................................
c
c ... Forward substitution: Ly = b 
       x(1) = x(1)*ad(1)  
       do 110 i = 2, n
        t = x(i)
        do 100 k = ia(i), ia(i+1) - 1
           t = t - a(k)*x(ja(k))
  100   continue
        x(i) = t * ad(i)
  110 continue
c ......................................................................
c
c ... Backward substitution: Lt x = y
      x(n) = x(n)*ad(n)
      do 210 i = n, 2, -1 
        t = x(i)
        do 200 k = ia(i), ia(i+1)-1
           kk    = ja(k)
           x(kk) = x(kk) - t*a(k)
  200   continue
        x(i-1) = x(i-1)*ad(i-1)
  210 continue
c ......................................................................
      ifatsolvtime = Mpi_Wtime() - ifatsolvtime    
      return
      end
c **********************************************************************
c
c **********************************************************************
      subroutine ildlt_csrc_aux(neq,ia,ja,jat,iat,at)
c **********************************************************************
c * Data de criacao    : 10/04/2016                                    *
c * Data de modificaco : 00/00/0000                                    *
c * ------------------------------------------------------------------ *  
c * ILDLT_CSRC_AUX: resolve o sistema com a fatoracao imcompleta LDLt  *
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ * 
c * neq        - numero de equacoes                                    *
c * ia(1:n+1)  - ponteiro das linhas, para o formato CSR               *
c * ja(1:nad)  - ponteiro das colunas no formato CSR                   *
c * jat(1:n+1) - nao definido                                          *
c * iat(1:nad) - nao definido                                          *
c * kat(1:nad  - nao definido                                          * 
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ * 
c * jat(1:n+1) - ponteiro das colunas, para o formato CSC              *
c * iat(1:nad) - ponteiro das linhas no formato CSC                    *
c * kat(1:nad  - ponteiro que relaciona o CSC com a matriz a do CSR    * 
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ *
c * Todas as informação da matriz são da parte inferior                * 
c **********************************************************************
 
      implicit none
      integer ia(*),ja(*),jat(*),iat(*),at(*)
      integer neq,nl,nc,nnc,k
c ...
      nnc          = 0
      jat(1:neq+1) = 0
c ... loop na colunas
      do nc = 1, neq
c ... loop nas linhas abaixo da diagonal
        do nl = nc + 1, neq
c ...
          do k = ia(nl), ia(nl+1) - 1
           
            if( nc .eq. ja(k)) then
              jat(nc)  = jat(nc) + 1 
              nnc      = nnc     + 1
              iat(nnc) = nl
              at(nnc)  = k
              go to 10
            endif
          enddo
c ....................................................................
  10      continue
        enddo
c ....................................................................
      enddo
c .....................................................................
c
      do k = neq, 1, -1
        jat(k+1) =  jat(k)    
      enddo
c     
      jat(1) = 1
      do k = 1, neq 
        jat(k+1) =  jat(k+1) + jat(k)   
      enddo
c .....................................................................
c
c ...
      return
      end
c **********************************************************************
c
c **********************************************************************
      subroutine ldltc(n,a,v,w)
c **********************************************************************
c * Data de criacao    : 20/06/2016                                    *
c * Data de modificaco : 00/00/0000                                    * 
c * ------------------------------------------------------------------ *  
c * LDLTC : fatoracao LDLt completa com a matriz cheia                 *
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ * 
c * n      - numero de equacoes                                        *
c * a(n,n) - matriz cheia                                              *
c * v(n)   - vetor auxiliar de trabalho                                *
c * w(n)   - vetor auxiliar de trabalho                                *
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ * 
c * a(n,n) - matriz fatorada LDLt                                      *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ *
c * parte triangular superior inalterada e diagonal invertida          *  
c **********************************************************************      
      implicit none
      integer n,i,j,k
      real*8 a(n,*),v(*),w(*),dot
c ......................................................................
c
c ...
      do j = 1, n
         do i = 1, j-1
            v(i) = a(j,i)*a(i,i)
         enddo
         w(1:j-1) = a(j,1:j-1)
         v(j) = a(j,j) - dot(w,v,j-1)
         if (v(j) .eq. 0.d0) then
           print*, 'Subrotina LDLTC: Pivot nulo !'
           stop
         endif 
         a(j,j) = v(j)
         do k = j+1, n
            w(1:j-1) = a(k,1:j-1)
            a(k,j) = (a(k,j) - dot(w,v,j-1))/v(j) 
         enddo
      enddo
c ......................................................................
c
c ...
      do j = 1, n
        a(j,j) = 1.d0/a(j,j)
      enddo
c ......................................................................
c
c ......................................................................
c     open(15,file='ldtl.matrix')
c     do i = 1, n
c       write(15,'(100d15.1)')(a(i,j),j=1,i)
c     enddo 
c     do i = 1, n
c        write(15,'(100d15.2)')(a(i,j),j=1,i)
c        write(15,'(100d15.2)')a(i,i)
c     enddo 
c     close(15)
c ...
      return
      end
c **********************************************************************
c
c **********************************************************************
      subroutine ildlt_full_matrix(n,a,v,w)
c **********************************************************************
c * Data de criacao    : 10/04/2016                                    *
c * Data de modificaco : 00/00/0000                                    * 
c * ------------------------------------------------------------------ *  
c * ILDLTC_FULL_MATRIX : fatoracao LDLt incompleta com a matriz cheia  *
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ * 
c * n      - numero de equacoes                                        *
c * a(n,n) - matriz cheia                                              *
c * v(n)   - vetor auxiliar de trabalho                                *
c * w(n)   - vetor auxiliar de trabalho                                *
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ * 
c * a(n,n) - matriz fatorada LDLt incompleta                           *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ *
c * parte triangular superior inalterada                               *  
c **********************************************************************      
      implicit none
      integer n,i,j,k
      real*8 a(n,*),v(*),w(*),dot
c ......................................................................
      do j = 1, n   
        do i = 1, j-1
          v(i) = a(j,i)*a(i,i)
        enddo
        w(1:j-1) = a(j,1:j-1)
        v(j) = a(j,j) - dot(w,v,j-1)
c
        if (v(j) .eq. 0.d0) then
          print*, 'Subrotina LDLTC: Pivot nulo !'
          stop
        endif 
c
        a(j,j) = v(j)
c
        do k = j+1, n
          w(1:j-1) = a(k,1:j-1)
          if( a(k,j) .ne. 0.0d0) then
            a(k,j) = (a(k,j) - dot(w,v,j-1))/v(j)
          endif 
        enddo
      enddo
c ......................................................................
c     open(15,file='ildtl.matrix')
c     do i = 1, 1500
c       write(15,'(i9,d15.2)')i,a(i,i)
c     enddo 
c     write(15,*)'l'
c     k = 1
c     do i = 1, n
c       do j = 1, i - 1
c         if( a(i,j) .ne. 0.0d0 ) then 
c           write(15,'(i9,d15.2)')k,a(i,j)
c           k = k + 1
c         endif
c       enddo
c     enddo 
c     close(15)
c ...
      return
      end
c **********************************************************************
c
c **********************************************************************
      subroutine choleskyc(n,a)
c **********************************************************************
c * Data de criacao    : 10/04/2016                                    *
c * Data de modificaco : 00/00/0000                                    * 
c * ------------------------------------------------------------------ *  
c * CHOLESKYC : fatoracao CHOLESKY completa com matriz cheia           *
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ * 
c * n      - numero de equacoes                                        *
c * a(n,n) - matriz cheia                                              *
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ * 
c * a(n,n) - matriz fatorada LDLt                                      *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ *
c * parte triangular superior inalterada e diagonal invertida          *  
c **********************************************************************      
      implicit none
      integer n,i,j,k
      real*8 a(n,*),t
c     open(15,file='ch.matrix')
c     do i = 1, n
c       write(15,'(100d15.1)')(a(i,j),j=1,n)
c     enddo 
c     write(15,'(a6)')'d'
c     do i = 1, n
c       write(15,'(d15.2)')a(i,i)
c     enddo
c     write(15,'(a6)')'d'
c     do i = 1, n
c       do j = 1, i - 1
c         write(15,'(d15.2)')a(i,j)
c        enddo 
c     enddo      
c     close(15)
c ......................................................................
      do j = 1, n
c ...
        t = 0.0d0
        do i = 1, j-1
          t = t + a(j,i)*a(j,i)
        enddo
c .....................................................................
c
c ...
        t = a(j,j) - t
        if (t .lt. 0.d0) then
          print*, 'Subrotina CHOLESKYC: Pivot negativo !',j
          stop
        endif 
        a(j,j) = dsqrt(t)
c .....................................................................
c
c ...
        do k = j+1, n
          t = 0.0d0
          do i = 1, j-1
            t = t +  a(k,i)*a(j,i)
          enddo
c ....................................................................
c
c ...
          a(k,j) = (a(k,j) - t) / a(j,j)
        enddo
      enddo
c ......................................................................
c
c ...
      do j = 1, n
        a(j,j) = 1.d0/a(j,j)
      enddo
c ......................................................................
c     open(15,file='ch.matrix')
c     do i = 1, n
c       write(15,'(100d15.1)')(a(i,j),j=1,i)
c     enddo 
c     write(15,'(a6)')'d'
c     do i = 1, n
c       write(15,'(d15.2)')a(i,i)
c     enddo
c     write(15,'(a6)')'d'
c     do i = 1, n
c       do j = 1, i - 1
c         write(15,'(d15.2)')a(i,j)
c        enddo 
c     enddo      
c     close(15)
c ...
      return
      end
c **********************************************************************
c
c **********************************************************************
      subroutine icholeskyc_full_matrix(n,a)
c **********************************************************************
c * Data de criacao    : 10/04/2016                                    *
c * Data de modificaco : 00/00/0000                                    * 
c * ------------------------------------------------------------------ *  
c * CHOLESKYC : fatoracao CHOLESKY completa com matriz cheia           *
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ * 
c * n      - numero de equacoes                                        *
c * a(n,n) - matriz cheia                                              *
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ * 
c * a(n,n) - matriz fatorada LDLt                                      *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ *
c * parte triangular superio inalterada                                *  
c **********************************************************************      
      implicit none
      integer n,i,j,k
      real*8 a(n,*),t
c ......................................................................
      do j = 1, 500
c ...
         t = 0.0d0
         do i = 1, j-1
            t = t +  a(j,i)*a(j,i)
         enddo
c .....................................................................
c
c ...
         t = a(j,j) - t
         if (t .lt. 0.d0) then
           print*, 'Subrotina CHOLESKYC: Pivot negativo !'
           stop
         endif 
         a(j,j) = dsqrt(t)
c .....................................................................
c
c ...
         do k = j+1, n
           t = 0.0d0
           do i = 1, j-1
             t = t +  a(k,i)*a(j,i)
           enddo
c .....................................................................
c
c ...      
           if( a(k,j) .ne. 0.d0 ) then 
             a(k,j) = (a(k,j) - t) / a(j,j)
           endif
c ..................................................................... 
         enddo
c .....................................................................
      enddo
c ......................................................................
c     open(15,file='ich.matrix')
c     do i = 1, n
c       write(15,'(100d15.1)')(a(i,j),j=1,i)
c     enddo 
c     write(15,*)'d'
c     do i = 1, 500
c       write(15,'(d15.2)')a(i,i)
c     enddo
c     write(15,*)'l' 
c     do i = 1, 500
c       do j = 1, i - 1
c         if ( a(i,j) .ne. 0.d0 ) then
c           write(15,'(d15.2)')a(i,j)
c         endif
c       enddo 
c     enddo      
c     close(15)
c ...
      return
      end
c **********************************************************************
c
c **********************************************************************
c * Data de criacao    : 18/04/2016                                    *
c * Data de modificaco : 00/00/0000                                    * 
c * ------------------------------------------------------------------ *  
c * SET_PRECOND : escolhe o precondicionandor                          *
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ * 
c * macro   - precondicionador escolhido                               *
c * precond - nao definido                                             *
c * nin     - aqruivo de entrada                                       *
c * my_id   - id do processo do mpi                                    *
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ * 
c * precond - precondicionador escolhido                               * 
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ *
c ********************************************************************** 
      subroutine set_precond(macro,precond,nin,my_id)
      implicit none
      include 'string.fi'
      include 'precond.fi'
      character macro(maxstrl)
      character*6 macros(7),string
      integer precond,size,nin,my_id
      integer i,nmc 
      data macros/'none  ','diag  ','ildlt '
     .           ,'illt  ','diagm ','bdiag '
     .           ,'diags '/
      data nmc /7/
c ...
      write(string,'(6a)') (word(i),i=1,6)
c ... nenhum precondicionador
      if( string .eq. macros(1)) then
        precond = 1
        if(my_id.eq.0) then
          write(*,'(1x,a10,1x,a6)')'Precond:',macros(1)  
        endif
c .....................................................................
c
c ... precondicionador diagonal
      elseif( string .eq. macros(2)) then
        precond = 2
        if(my_id.eq.0) then
          write(*,'(1x,a10,1x,a6)')'Precond:',macros(2)
        endif
c .....................................................................
c
c ... precondicionador ILDLT
      elseif( string .eq. macros(3)) then
        precond = 3
        if(my_id.eq.0) then
          write(*,'(1x,a10,1x,a6)')'Precond:',macros(3)
        endif
c .....................................................................
c
c ... precondicionador ILLT ( Cholesky )
      elseif( string .eq. macros(4)) then
        precond = 4
        if(my_id.eq.0) then
          write(*,'(1x,a10,1x,a6)')'Precond:',macros(4)
        endif
c .....................................................................
c
c ... precondicionador ILLT ( Cholesky )
      elseif( string .eq. macros(5)) then
        precond = 5
        if(my_id.eq.0) then
          write(*,'(1x,a10,1x,a6)')'Precond:',macros(5)
        endif
c .....................................................................
c
c ... precondicionador Block diagnal 
      elseif( string .eq. macros(6)) then
        precond = 6
c ... numero de blocos
        call size_Block(size,nin)
c .....................................................................
c
c ... 
        iparam(3) = size*(size+1)/2
        iparam(4) = size
c .....................................................................
c
c ...  
        if(my_id.eq.0) then
          write(*,'(1x,a10,1x,a6,i3)')'Precond:',macros(6),iparam(4)
        endif
c .....................................................................
c
c ... precondicionador Block diagnal 
      elseif( string .eq. macros(7)) then
        precond = 7
        if(my_id.eq.0) then
          write(*,'(1x,a10,1x,a6)')'Precond:',macros(7)
        endif
c .....................................................................
c
c ...                         
      else
        print*,'Erro na leitura da macro precond !'
        stop
      endif 
c .....................................................................
c
c ...
      return
      end
c ********************************************************************** 
c
c c ********************************************************************
c * Data de criacao    : 27/06/2016                                    *
c * Data de modificaco : 00/00/0000                                    * 
c * ------------------------------------------------------------------ *  
c * CAL_PRECOND : calculo o precondicionador                           *
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ * 
c * ia     - csrc                                                      *
c * ja     - csrc                                                      *
c * m      - indefinido                                                *
c * ad     - coeficientes da diagonal principal                        *
c * al     - coeficientes da parte inferior da matriz                  *
c * w(neq) - vetor auxiliar                                            *
c * precond- codigo para o precondicionador                            *
c * neq    - numero de equacoes                                        *
c * nequ   - numero de equacoes no bloco Kuu                           *
c * my_id  - id do processo do mpi                                     *
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ * 
c * m      - precondicionador                                          *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ *
c ********************************************************************** 
      subroutine cal_precond(ia,ja,m,ad,al,w,precond,neq,nequ,my_id)
      implicit none
      include 'mpif.h'
      include 'precond.fi'
      include 'time.fi'
      real*8 m(*),ad(*),al(*),w(*)
      integer ia(*),ja(*),neq,nequ
      integer precond,my_id
      real*8  max_block_a(max_block*max_block)
c ... sem precondicionador:
      if(precond .eq. 1) then
c .....................................................................
c
c ... precondicionador diagonal:
      else if(precond .eq. 2) then 
        precondtime = Mpi_Wtime() - precondtime  
        call pre_diag(m,ad,neq,.false.)  
        precondtime = Mpi_Wtime() - precondtime 
c .....................................................................
c
c ... precondicionador LDLT incompleto
      else if(precond .eq. 3) then
c ...
        precondtime = Mpi_Wtime() - precondtime 
        call ildlt2(neq,ia,ja,al,ad,m,w,0.0d0,.false.)
        precondtime = Mpi_Wtime() - precondtime 
c .....................................................................
c
c ... precondicionador Cholesky LLT incompleto
      else if(precond .eq. 4) then
c ...
        precondtime = Mpi_Wtime() - precondtime 
        call ichfat(neq,ia,ja,al,ad,m,w,0.0d0,.true.)
        precondtime = Mpi_Wtime() - precondtime 
c .....................................................................
c
c ... precondicionador modulo da diagonal:
      else if(precond .eq. 5) then
c ...
        precondtime = Mpi_Wtime() - precondtime  
        call pre_diag(m,ad,neq,.true.)
        precondtime = Mpi_Wtime() - precondtime 
c .....................................................................
c
c ... precondicionador modulo da diagonal:
      else if(precond .eq. 6) then
c ...
        precondtime = Mpi_Wtime() - precondtime 
        call block_precond(ad,al,ia,ja,m,neq,max_block_a,iparam) 
        precondtime = Mpi_Wtime() - precondtime
c .....................................................................
c
c ... precondicionador diagonal com complemento de schur:
      else if(precond .eq. 7) then
c ...
        precondtime = Mpi_Wtime() - precondtime  
        call pre_diag_schur(ia,ja,m,ad,al,w,neq,nequ)
        precondtime = Mpi_Wtime() - precondtime 
c .....................................................................
c
c ...
      else
        if( my_id.eq.0 ) then
          print*,"Precond invalido: !!",precond
          call stop_mef()
        endif 
      endif
c ......................................................................
      return      
      end
c **********************************************************************
c
c **********************************************************************
c * Data de criacao    : 18/06/2016                                    *
c * Data de modificaco : 00/00/0000                                    * 
c * ------------------------------------------------------------------ *  
c * inverse_matrix_3x3_sym : inverte um matriz simeterica 3x3          *                          *
c * ------------------------------------------------------------------ * 
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ * 
c * m(6,*)  - indefinido                                               *
c * a11     - coeficiente da posicao a(1,1)                            *
c * a22     - coeficiente da posicao a(2,2)                            *
c * a33     - coeficiente da posicao a(3,3)                            *
c * a21     - coeficiente da posicao a(2,1)                            *
c * a31     - coeficiente da posicao a(3,1)                            *
c * a32     - coeficiente da posicao a(3,2)                            *
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ * 
c * precond - precondicionador escolhido                               * 
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ *
c * Inversao feita atraves de cofatores                                *  
c **********************************************************************
      subroutine inverse_matrix_3x3_sym(m,a11,a22,a33,a21,a31,a32)
      implicit none
      real*8 m(6,*),a11,a22,a33,a21,a31,a32
      real*8 c11,c22,c33,c21,c31,c32,det
c ... cofatores
      c11 = a22*a33 - a32*a32
c
      c22 = a11*a33 - a31*a31
c     
      c33 = a11*a22 - a21*a21
c
      c21 = a31*a32 - a21*a33
c
      c31 = a21*a32 - a31*a22
c
      c32 = a31*a21 - a11*a32
c .....................................................................
c
c ... determinante 
      det = a11*a22*a33 + a21*a32*a31 + a31*a21*a32
     .    -(a22*a31*a31 + a11*a32*a32 + a33*a21*a21)
      det = 1.d0/det 
c .....................................................................
c
c ...
      a11 = c11*det
      a22 = c22*det
      a33 = c33*det
c ...
      a21 = c21*det
      a31 = c31*det
      a32 = c32*det
c ......................................................................
      return
      end
c **********************************************************************
c
c **********************************************************************
      subroutine size_Block(size,nin)
      implicit none
      include 'string.fi'
      character*30 string
      integer size,i,nin
c ...
      call readmacro(nin,.false.)
      write(string,'(30a)') (word(i),i=1,30)
      read(string,*,err =10,end =10) size
c .....................................................................
c
c ...
      return
   10 continue
      print*,'Erro na leitura da tamanho do bloco diagonal !'
      stop
      end
c **********************************************************************
c
c **********************************************************************
c * Data de criacao    : 19/06/2016                                    *
c * Data de modificaco : 00/00/0000                                    *
c * ------------------------------------------------------------------ *
c * solv_cholesky : solver direto para fatoracao LLt                   *
c * ------------------------------------------------------------------ *
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ *
c * cllt(n,n) - fatora LLt                                             *
c * b         - vetor de forcas                                        *
c * x         - nao definido                                           *
c * n         - numero de equacoes                                     *
c * ------------------------------------------------------------------ *
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ *
c * x         - solucao                                                *
c * ------------------------------------------------------------------ *
c * OBS:                                                               *
c * ------------------------------------------------------------------ *
c * Diagonal LLt ja invertida informacoes da parte inferior e superior *
c * loop interno por coluna                                            *
c **********************************************************************
      subroutine solv_cholesky(cllt,b,x,n)
      implicit none
      real*8 cllt(n,*),b(*),x(*),t
      integer i,j,n
c ...
      x(1:n) = b(1:n)
c .....................................................................
c
c ... Gy = b
      x(1) = x(1)*cllt(1,1)
      do i =2, n
        t = x(i)
        do j = 1, i - 1
c         t = t - cllt(i,j)*x(j) 
          t = t - cllt(j,i)*x(j) 
        enddo
        x(i) = t*cllt(i,i)
      enddo
c .....................................................................
c
c ... Gtx = y
      x(n) = x(n)*cllt(n,n)
      do i =n-1, 1, -1
        t = x(i)
        do j = i+1, n
          t = t - cllt(j,i)*x(j) 
        enddo
        x(i) = t*cllt(i,i)
      enddo
c .....................................................................
c
c ... 
      return
      end
c *********************************************************************
c
c **********************************************************************
c * Data de criacao    : 20/06/2016                                    *
c * Data de modificaco : 00/00/0000                                    *
c * ------------------------------------------------------------------ *
c * solv_ldtl : solver direto para fatoracao LDLt                      *
c * ------------------------------------------------------------------ *
c * Parametros de entrada:                                             *
c * ------------------------------------------------------------------ *
c * ldlt(n,n) - fatora LDLt                                            *
c * b         - vetor de forcas                                        *
c * x         - nao definido                                           *
c * n         - numero de equacoes                                     *
c * ------------------------------------------------------------------ *
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ *
c * x         - solucao                                                *
c * ------------------------------------------------------------------ *
c * OBS:                                                               *
c * ------------------------------------------------------------------ *
c * Diagonal LDLt ja invertida                                         *
c * informacoes da parte inferior e superior                           *
c * loop interno por coluna                                            *
c **********************************************************************
      subroutine solv_ldlt(cllt,b,x,n)
      implicit none
      real*8 cllt(n,*),b(*),x(*),t
      integer i,j,n
c ...
      x(1:n) = b(1:n)
c .....................................................................
c
c ... Gy = b
      do i =2, n
        t = x(i)
        do j = 1, i - 1
c         t = t - cllt(i,j)*x(j) 
          t = t - cllt(j,i)*x(j) 
        enddo
        x(i) = t
      enddo
c .....................................................................
c
c ...
      do i = 1, n
        x(i) = x(i)*cllt(i,i)
      enddo
c .....................................................................
c
c ... Gtx = y
      do i =n-1, 1, -1
        t = x(i)
        do j = i+1, n
          t = t - cllt(j,i)*x(j) 
        enddo
        x(i) = t
      enddo
c .....................................................................
c
c ... 
      return
      end
c *********************************************************************
      
       
