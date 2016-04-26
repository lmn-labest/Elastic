       subroutine interpol(x,y,xi,n,c)
c **********************************************************************
c *                                                                    *
c *   Interpol:                                                         *
c *   -----                                                            *
c *                                                                    *
c *   Interpola valores com base em uma poligonal                      *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ---------------------                                            *
c *                                                                    *
c *    x - coluna dos valores do eixo x                                *
c *    y - coluna dos valores do eixo y                                *
c *    xi - valor de x para o qual se quer determinar o y              *
c *    n - numero de pontos da poligonal                               *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   -------------------                                              *
c *                                                                    *
c *    c - valor de y correspondente ao xi                             *
c *                                                                    *
c *                                                                    *
c **********************************************************************  
      implicit none
      integer n,i
      real*8 x(*),y(*),xi,c,x1,x2,y1,y2
      if  (xi .lt. x(1)) then
         c = y(1)
         goto 1000
      elseif (xi .gt. x(n)) then
         c = y(n)
         goto 1100      
      endif
      do i = 1, n-1
         x1 = x(i)
         x2 = x(i+1)
         y1 = y(i)
         y2 = y(i+1)
         if (xi .ge. x1 .and. xi .le. x2) then
            c = y1 + (y2-y1)*(xi-x1)/(x2-x1)
            goto 100
         endif
      enddo
  100 continue
      return
 1000 continue
c      print*, '*** valor abaixo do limite de interpolacao ',xi,x(1)
c      stop
      return
 1100 continue
c      print*, '*** valor acima do limite de interpolacao ',xi,x(n)
c      stop
      return
      end
      
   