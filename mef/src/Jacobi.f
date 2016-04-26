      subroutine jacobi(a,x,w,ip,n,tol,maxciclos)
c **********************************************************************
c *                                                                    *
c *   JACOBI                                                           *
c *   ------                                                           *
c *                                                                    *
c *   Calculo de autovalores e autovetores de matrizes simetricas      *
c *   pelo metodo de JACOBI.                                           *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ---------------------                                            *
c *                                                                    *
c *   a(n,n) - matriz simetrica                                        *
c *   x(n,n) - nao definido                                            *
c *   w(n)   - nao definido                                            *
c *   ip(n)  - arranjo auxiliar (nao definido)                         *
c *   n      - dimensao da matriz                                      *
c *   tol    - tolerancia de convergencia                              *
c *   maxciclos - nmumero maximo de ciclos                             *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   --------------------                                             *
c *                                                                    *
c *   a - matriz espectral: a(i,i) = i-esimo autovalor                 *
c *   x - matriz de autovetores: x(i,j),(i=1,n) = autovetor            *
c *                              associado ao j-esimo autovalor        *
c *   w - w(i) = (a(i,i) = i-esimo autovalor                           *
c *                                                                    *
c **********************************************************************
      implicit none
      integer k
      integer ip(n),n,p,q,maxciclos
      real*8 a(n,n),x(n,n),w(n),tol,diag,off,conv
      real*8 t,c,s,d,r,r1,r2
c ......................................................................
c
c ... Inicializa a matriz x (matriz identidade):
c
      call initx(x,n)
c
c ... Ciclos de Jacobi:
c
      call diagnorm(a,n,diag)
      do 1000 k = 1, maxciclos
        do  110 p = 1, n
          do  100 q = p+1, n
c
c ... Calculo da rotacao para o coeficiente a(p,q):
c
            c = 1.d0
            s = 0.d0
            if(a(p,q) .ne. 0.d0) then
              r = (a(q,q)-a(p,p))/(2.d0*a(p,q))
              d = dsqrt(1.d0+r*r)
              r1 = -r + d
              r2 = -r - d
              t = r1
              if(dabs(r2).lt.dabs(r1)) then
                t = r2
              endif
              c = (1.d0/dsqrt(1.d0+t*t))
              s = (t*c)
            endif
c                   t
c ... Produto X A X:
c
            call xax(a,c,s,n,p,q)
c
c ... Produto (Xi)(Xi+1):
c                   
            call xx(x,c,s,n,p,q)
  100     continue
  110   continue
        call offnorm(a,n,off)
        conv = off/diag
        if(conv .lt. tol) goto 2000
 1000 continue
      print*,'JACOBI: Sem convergecia apos ',maxciclos,' ciclos !'
      stop
 2000 continue
c
c ...    Ordena os autovalores e autovetores:
c
      call ordena(a,x,w,ip,n)
      return
      end
      subroutine xax(a,c,s,n,p,q)
c **********************************************************************
c *                                                                    *
c *      XAX                                                           *
c *   ---                                                              *
c *                                                                    *
c *                                                                    *
c **********************************************************************
      implicit none
      integer i,j
      integer n,p,q
      real*8  a(n,n),c,s,a1,a2
c ......................................................................
c
c ... Produto AX:
c
      do 100 i = 1, n
        a1 = a(i,p)
        a2 = a(i,q)
        a(i,p) = a1*c - a2*s
        a(i,q) = a1*s + a2*c
  100 continue
c                   t
c ... Produto X(AX):
c
      do 200 j = 1, n
        a1 = a(p,j)
        a2 = a(q,j)
        a(p,j) = c*a1 - s*a2
        a(q,j) = s*a1 + c*a2
  200 continue
      return
      end
      subroutine xx(x,c,s,n,p,q)
c **********************************************************************
c *                                                                    *
c *      XX                                                            *
c *   --                                                               *
c *                                                                    *
c *                                                                    *
c **********************************************************************
      implicit none
      integer i
      integer n,p,q
      real*8  x(n,n),c,s,a1,a2
c ......................................................................
      do 100 i = 1, n
        a1 = x(i,p)
        a2 = x(i,q)
        x(i,p) =  a1*c - a2*s
        x(i,q) =  a1*s + a2*c
  100 continue
      return
      end
      subroutine initx(x,n)
c **********************************************************************
c *                                                                    *
c *      INITX                                                         *
c *   -----                                                            *
c *                                                                    *
c *                                                                    *
c **********************************************************************
      implicit none
      integer i,j
      integer n
      real*8  x(n,n)
c ......................................................................
      do 110 i = 1, n
         x(i,i) = 1.d0
         do 100 j = i+1, n
            x(i,j) = 0.d0
            x(j,i) = 0.d0
  100    continue
  110 continue
      return
      end
      subroutine ordena(a,x,w,ip,n)
c **********************************************************************
c *                                                                    *
c *      ORDENA                                                        *
c *   ------                                                           *
c *                                                                    *
c *                                                                    *
c **********************************************************************
      implicit none
      integer i,j,k
      integer ip(n),n,iaux
      real*8  a(n,n),x(n,n),w(n),aux
      logical itroca
c ......................................................................
      do 100 i = 1, n
         ip(i) = i
         w(i)  = a(i,i)
  100 continue
  200      continue
      itroca = .false.
      do 210 i = 2, n
         j = i-1
         if(w(j).gt.w(i)) then
            aux   = w(i)
            w(i) = w(j)
            w(j) = aux
            iaux  = ip(i)
            ip(i) = ip(j)
            ip(j) = iaux
            do 110 k = 1, n
               aux = x(k,i)
               x(k,i) = x(k,j)
               x(k,j) = aux
  110       continue
            itroca = .true.
         endif
  210 continue
      if(itroca) goto 200
      return
      end
      subroutine diagnorm(a,n,diag)
c **********************************************************************
c *                                                                    *
c *      DIAGNORM                                                      *
c *   --------                                                         *
c *                                                                    *
c *                                                                    *
c **********************************************************************
      implicit none
      integer i
      integer n
      real*8  a(n,n),diag
c ......................................................................
      diag = 0.d0
      do 100 i = 1, n
         diag = diag + a(i,i)*a(i,i)
  100 continue
      return
      end
      subroutine offnorm(a,n,off)
c **********************************************************************
c *                                                                    *
c *      DIAGNORM                                                      *
c *   --------                                                         *
c *                                                                    *
c *                                                                    *
c **********************************************************************
      implicit none
      integer i,j
      integer n
      real*8  a(n,n),off
c ......................................................................
      off = 0.d0
      do 110 i = 1, n
        do 100 j = i+1,n
          off = off + a(i,j)*a(i,j)
  100   continue
  110 continue
      off = off*2.d0
      return
      end
      subroutine xaxold(a,c,s,n,p,q)
c **********************************************************************
c *                                                                    *
c *      XAX                                                           *
c *   ---                                                              *
c *                                                                    *
c *                                                                    *
c **********************************************************************
      implicit none
      integer i,j
      integer n,p,q
      real*8  a(n,n),c,s,a1,a2
c ......................................................................
c
c ... Produto AX:
c
      do 100 i = 1, n
        a1 = a(i,p)*c - a(i,q)*s
        a2 = a(i,p)*s + a(i,q)*c
        a(i,p) = a1
        a(i,q) = a2
  100 continue
c                   t
c ... Produto X(AX):
c
      do 200 j = 1, n
        a1 = c*a(p,j) - s*a(q,j)
        a2 = s*a(p,j) + c*a(q,j)
        a(p,j) = a1
        a(q,j) = a2
  200 continue
      return
      end
      subroutine xxold(x,c,s,n,p,q)
c **********************************************************************
c *                                                                    *
c *      XX                                                            *
c *   --                                                               *
c *                                                                    *
c *                                                                    *
c **********************************************************************
      implicit none
      integer i
      integer n,p,q
      real*8  x(n,n),c,s,a1,a2
c ......................................................................
      do 100 i = 1, n
        a1 =  x(i,p)*c - x(i,q)*s
        a2 =  x(i,p)*s + x(i,q)*c
        x(i,p) = a1
        x(i,q) = a2
  100 continue
      return
      end

