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
c * ------------------------------------------------------------------ * 
c * Parametros de saida:                                               *
c * ------------------------------------------------------------------ * 
c * m   - precondicionador diagonal ( M-1)                             *
c * ------------------------------------------------------------------ * 
c * OBS:                                                               *
c * ------------------------------------------------------------------ * 
c **********************************************************************
      subroutine pre_diag(m,ad,neq)
      implicit none
      real*8 m(*),ad(*)
      integer i,neq
      do i = 1, neq
        m(i) = 1.d0/ad(i)
      enddo
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
          print*,'Valor do Pivo invalido !!!',v(j)
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
          print*,'Valor do Pivo invalido !!!',v(j)
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
          print*,'Valor do Pivo invalido !!!',v(j)
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
c * Data de criacao    : 10/04/2016                                    *
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
c * parte triangular superio inalterada                                *  
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
c * parte triangular superio inalterada                                *  
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
      close(15)
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
c * parte triangular superio inalterada                                *  
c **********************************************************************      
      implicit none
      integer n,i,j,k
      real*8 a(n,*),t
c ......................................................................
      do j = 1, n
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
           a(k,j) = (a(k,j) - t) / a(j,j)
         enddo
      enddo
c ......................................................................
      open(15,file='ch.matrix')
c     do i = 1, n
c       write(15,'(100d15.1)')(a(i,j),j=1,i)
c     enddo 
      write(15,'(a6)')'d'
      do i = 1, n
        write(15,'(d15.2)')a(i,i)
      enddo
      write(15,'(a6)')'d'
      do i = 1, n
        do j = 1, i - 1
          write(15,'(d15.2)')a(i,j)
         enddo 
      enddo      
      close(15)
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
      open(15,file='ich.matrix')
c     do i = 1, n
c       write(15,'(100d15.1)')(a(i,j),j=1,i)
c     enddo 
      write(15,*)'d'
      do i = 1, 500
        write(15,'(d15.2)')a(i,i)
      enddo
      write(15,*)'l'
      do i = 1, 500
        do j = 1, i - 1
          if ( a(i,j) .ne. 0.d0 ) then
            write(15,'(d15.2)')a(i,j)
          endif
        enddo 
      enddo      
      close(15)
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
      character macro(maxstrl)
      character*6 macros(6),string
      integer precond,nin,my_id
      integer i,nmc 
      data macros/'none  ','diag  ','ildlt '
     .           ,'illt  ','      ','      '/
      data nmc /6/
c ...
      write(string,'(6a)') (word(i),i=1,6)
c ... nenhum precondicionador
      if( string .eq. macros(1)) then
        precond = 1
        if(my_id.eq.0) then
          write(*,'(1x,a25,1x,a6)')'precond:',macros(1)  
        endif
        print*,'Precond nap disponivel !!!',macros(1) 
        stop
c .....................................................................
c
c ... precondicionador diagonal
      elseif( string .eq. macros(2)) then
        precond = 2
        if(my_id.eq.0) then
          write(*,'(1x,a25,1x,a6)')'precond:',macros(2)
        endif
c .....................................................................
c
c ... precondicionador ILDLT
       elseif( string .eq. macros(3)) then
        precond = 3
        if(my_id.eq.0) then
          write(*,'(1x,a25,1x,a6)')'precond:',macros(3)
        endif
        print*,'Precond nap disponivel !!!',macros(3) 
        stop
c .....................................................................
c
c ... precondicionador ILLT ( Cholesky )
       elseif( string .eq. macros(4)) then
        precond = 4
        if(my_id.eq.0) then
          write(*,'(1x,a25,1x,a6)')'precond:',macros(4)
        endif
        print*,'Precond nap disponivel !!!',macros(4) 
        stop
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