      subroutine check_var_term(a,b,c,d11,d12,d22,d33,h,u,ym,pr,alpha,
     .                        nen,ndm,ma,flaghidr,epd)
c **********************************************************************
c *                                                                    *
c *   CHECK_VAR_HIDR:                                                  *
c *   -----                                                            *
c *                                                                    *
c *   Verifica se as propriedades termicas variam com a temperatura    *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ---------------------                                            *
c *                                                                    *
c *    ma - material                                                   *
c *    hi - grau de hidratacao                                         *
c *    e - numero da propriedade                                       *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   -------------------                                              *
c *                                                                    *
c *    c - valor da propriedade                                        *
c *                                                                    *
c *                                                                    *
c **********************************************************************    
      implicit none
      integer ma,nen,ndm,i
      real*8 ym,pr,alpha,h(nen),u(nen),tm,a,b,c,d11,d12,d22,d33,b1,b2
     .       ,a1,a2,a3
      logical flaghidr,epd
      if (flaghidr)  then
         return
      else
         tm = 0.d0
         do i = 1, nen
            tm = tm + h(i)*u(i)
         enddo
         call updprop(ym,tm,ma,1)
         call updprop(pr,tm,ma,2)
         call updprop(alpha,tm,ma,4)
         if( ndm .eq. 2) then
            if(epd)  then
               b1 = 1.d0-pr
               b2 = 1.d0-2.d0*pr
               a = ym*b1/((1.+pr)*b2)
               b = pr/b1
               c = b2/(2.d0*b1)
               d11 = a
               d12 = a*b
               d22 = a
               d33 = a*c  
             else
                a = ym/(1.d0-pr*pr)
                d11 = a
                d12 = a*pr
                d22 = a
                d33 = ym/(2.d0*(1.d0+pr))
             endif         
         elseif(ndm .eq. 3) then
            a1 = 1.d0+pr
           a2 = 1.d0-pr
           a3 = 1.d0-2.d0*pr
           a = ym*a2/(a1*a3)
           b = pr/a2
           c = a3/(2.d0*a2)
         endif
      endif
      return
      end      
      subroutine updprop(c,tm,ma,e)
c **********************************************************************
c *                                                                    *
c *   UPDPROP:                                                         *
c *   -----                                                            *
c *                                                                    *
c *   Interpola valores com base em uma poligonal                      *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ---------------------                                            *
c *                                                                    *
c *    ma - material                                                   *
c *    tm - temperatura no ponto                                       *
c *    e - numero da propriedade                                       *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   -------------------                                              *
c *                                                                    *
c *    c - valor da propriedade                                        *
c *                                                                    *
c *                                                                    *
c **********************************************************************    
      implicit none
      include 'termprop.fi'
      integer ma,e,n
      real*8 c,tm
      n = eprop(e,ma)
      if (n .ne. 0) then
          call interpol(nprop(1,e,ma),nprop(n+1,e,ma),tm,n,c)
      endif
      return
      end