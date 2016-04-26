       subroutine plasticst2d(a1,b1,c1,deps,eps,tx,txp,sy,hl,phi,pr,
     .          iyied,c14,g,epd)
c **********************************************************************
c *                                                                    *
c *                                                                    *
c *   PLASTICST2D: calcula tensoes 2D, considerando plasticidade       *
c *   -----------                                                      *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ---------------------                                            *
c *                                                                    *
c *     eps(4)  - deformacao plastica + elastica                       *
c *     a1,b1,c1   - coeficientes da matriz constitutiva               *
c *     sy      - tensao de escoamento                                 *
c *     deps(4) - delta deformacao                                     *
c *     tx(4)   - tensoes anteriores                                   *
c *     tx(5)   - tensao efetiva no passo anterior                     *
c *     iyied   - numero correspondente ao criterio a ser utilizado    *
c *     txp(4)   - tensao plastica acumulada ate o passo anterior      *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   -------------------                                              *
c *                                                                    *
c *     tx(4)   - tensoes atualizadas: tx = tx + dtx                   *
c *     txp(4)  - tensoes plasticas atualizadas                        *
c *     tx(5)   - tensao efetiva atualizada                            *
c *                                                                    *
c *                                                                    *
c **********************************************************************
      implicit none
      logical epd
      integer i,iyied,mstep,j,ist
      real*8  eps(*),deps(*),tx(*),txp(*),a1,b1,c1,epe,g,c14,dtxp(4)
      real*8  tx0(4),dtx(4),r,re,ratio,ddeps(4),phi,snphi,pr,se
      real*8  preys,sy,hl,dste(4),steff,varj2,theta,sint3,devs(4)
      real*8  se0,escur,rfac,astep,stot(4),reduc,abeta,epcon,agash,bgash
      real*8  dlamb,a(4),d(4),delep,curys,bring,str(4)
c ......................................................................
c   
c      call radian(phi)
      snphi =dsin(phi)
      epe = eps(5)
      se0 = tx(5)
c ... calculate elastic stress increment
      call stress2d0(a1,b1,c1,deps,dtx)
      preys=sy+epe*hl
      dste(1) = dtx(1)
      dste(2) = dtx(4)
      dste(3) = dtx(2)
      dste(4) = pr*(dste(1)+dste(3))
c ...    t0 = t:
      call aequalb(tx0,tx,4)
c ...    Tensoes lineares, t = t0 + dt
      call vsum(tx0,dste,4,str)
      if (iyied .eq. 0) goto 4
c ... evaluate effective stress
      call effst(se,iyied,str,devs,steff,theta,varj2,sint3,snphi)     
      if ( (se0-preys) .ge. 0.) goto 3
c ... check for yielding durint this iterations
      escur = se - preys
      if(escur .le. 0.) goto 4
c ... compute factor R
      rfac = escur/(se-se0)
      goto 5
    3 escur = se - se0
      if(escur .le. 0.d0) goto 4
      rfac = 1.
c ... evaluate number of stress reduction steps
    5 mstep = escur*8./sy+1.
      astep = mstep
      reduc = 1.-rfac
      do 6 j = 1,4
         stot(j) = tx0(j)+reduc*dste(j)
    6    dste(j) = rfac*dste(j)/astep
      epcon = 0.d0
c ... loop over each stress reduction step
      do 8 ist = 1, mstep
c ... calculate vectors A and D
         call effst(se,iyied,stot,devs,steff,theta,varj2,sint3,snphi)
         call flow(a,d,devs,abeta,iyied,steff,theta,varj2,sint3,hl,
     .          snphi,c14,g,epd)
c ... compute plastic multiplier
         agash = 0.
         do 9 j = 1,4
    9       agash = agash + a(j) * dste(j)
         dlamb = agash *abeta
         if ( dlamb .lt. 0.d0) dlamb = 0.d0
c ...compute elastoplastic stresses
         bgash = 0.d0
         do 10 j = 1,4
            bgash = bgash + a(j) * stot(j)
   10       stot(j) = stot(j) + dste(j)-dlamb*d(j)
c ... calculate equivalent plastic strain increment
         delep = dlamb*bgash/se
         epcon = epcon+delep
c ... calculae equivalent plastic strain
         eps(5) = eps(5) + delep
    8  continue
c ... compute effective stress
      call effst(se,iyied,stot,devs,steff,theta,varj2,sint3,snphi)
c ... compute equivalent yield stress
      curys = sy + eps(5)*hl
c ... scale down stresses to yield surface
      bring = 1.d0
      if ( se .gt. curys)    bring = curys/se
      do 11 j = 1,4
   11    tx(j) = bring*stot(j)
      tx(5) = bring * se
c     compute residual stresses
      do 12 j = 1, 4
         dtxp(j) = str(j) - tx(j)
   12    txp(j) = txp(j) + dtxp(j)
      goto 2
c ... compute stresses for elastic nodes or points
    4 continue
      do 13 i = 1, 4
   13 tx(i) = str(i)
      tx(5) = se
    2 continue
      return
      end
      subroutine effst(se,iyied,tx,devs,steff,theta,varj2,sint3,snphi)
c **********************************************************************
c *                                                                    *
c *                                                                    *
c *   Effst: calculates the deviatoric stresses, effective stress and  *
c *          invariants                                                *
c *   -----------                                                      *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ---------------------                                            *
c *                                                                    *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   -------------------                                              *
c *                                                                    *
c *                                                                    *
c *                                                                    *
c **********************************************************************
      implicit none
      integer iyied
      real*8 root3,sm,tx(*),devs(*),varj2,varj3,steff,sint3,se
      real*8 snphi,theta,arsin
      root3 = 1.73205080757d0
      sm = (tx(1)+tx(3)+tx(4))/3.d0
      devs(1) = tx(1)-sm
      devs(3) = tx(3)-sm
      devs(4) = tx(4)-sm
      devs(2) = tx(2)
      varj2 = devs(2)*devs(2)+0.5d0*(devs(1)*devs(1)+devs(3)*devs(3)+
     . devs(4)*devs(4))
c      varj3 = devs(4)*(devs(4)*devs(4)-varj2)
      varj3 = devs(4)*(-devs(2)*devs(2)+devs(1)*devs(3))
      steff = sqrt(varj2)
      sint3 = -3.d0*root3*varj3/(2.d0*varj2*steff)
      if ( abs(sint3) .gt. 1.d0)  sint3 = sign(1.d0,sint3)
      theta = asin(sint3)/3.d0
      goto (1,2,3,4,5) iyied
c ... Tresca
    1 se = 2.d0*cos(theta)*steff
      return
c ... Von Mises
    2 se = root3*steff
      return
c ... Mohr-Coulomb
    3 se = sm*snphi+steff*(dcos(theta)-dsin(theta)*snphi/root3)
      return
c ... Drucker-Prager
    4 se = 6.d0*sm*snphi/(root3*(3.d0-snphi))+steff
      return
    5 se = 3.d0*sm*snphi/(root3*dsqrt(3.d0+snphi*snphi))+steff
      return
      end
      subroutine flow(a,d,devs,abeta,iyied,steff,theta,varj2,sint3,hl,
     .                              snphi,c14,g,epd)
c **********************************************************************
c *                                                                    *
c *                                                                    *
c *   Effst: calculates vectors A and D                                *
c *   -----------                                                      *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ---------------------                                            *
c *                                                                    *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   -------------------                                              *
c *                                                                    *
c *                                                                    *
c *                                                                    *
c **********************************************************************
      implicit none
      logical epd
      integer iyied,j,ipl
      real*8 a1(4),a2(4),a3(4),a(4),d(4),root3,devs(4),steff,varj2
      real*8 cons1,cons2,cons3,theta,snphi,cons,costh,tanth,tant3
      real*8 plumi,c14,xml,sinth,g,abeta,sint3,hl,denom
      root3 = 1.73205080757d0
c ... vector A1
      a1(1) = 1.d0
      a1(2) = 0.d0
      a1(3) = 1.d0
      a1(4) = 1.d0
c ... vector A2
      a2(1) = devs(1)/(2.d0*steff)
      a2(2) = devs(2)/steff
      a2(3) = devs(3)/(2.d0*steff)
      a2(4) = devs(4)/(2.d0*steff)
c ... vector A3
c      a3(1)=(2.d0*devs(3)*devs(4)-devs(1)*devs(4)-devs(1)*devs(3)
c     .       -devs(2)*devs(2))/3.d0
c      a3(2) = 2.d0*devs(2)*devs(4)
c      a3(3)=(2.d0*devs(1)*devs(4)-devs(3)*devs(4)-devs(1)*devs(3)
c     .       -devs(2)*devs(2))/3.d0
c      a3(4)=(2.d0*devs(1)*devs(3)-devs(3)*devs(4)-devs(4)*devs(1)
c     .       +2.d0*devs(2)*devs(2))/3.d0     
      a3(1) = devs(3)*devs(4)+varj2/3.d0
      a3(2) = -2.d0*devs(2)*devs(4)
      a3(3) = devs(1)*devs(4)+varj2/3.d0
      a3(4) = devs(1)*devs(3)-devs(2)*devs(2)+varj2/3.d0
      goto (1,2,3,4,5),iyied
c ... Tresca
    1 cons1 = 0.
      if (abs(theta)*57.29577951308d0 .lt. 29.d0) goto 20
      cons2 = root3
      cons3 = 0.d0
      goto 40
   20 sinth = sin(theta)
      cons2 = 2.d0*(cos(theta)+sinth*tan(3.d0*theta))
      cons3 = root3*sinth/(varj2*cos(3.d0*theta))
      goto 40
c ... Von Mises
    2 cons1 = 0.d0
      cons2 = root3
      cons3 = 0.d0
      goto 40
c ... Mohr-Coulomb
    3 cons1 = snphi/3.d0
      if (abs(theta)*57.29577951308d0 .lt. 29.d0) goto 30
      plumi = 1.d0
      if ( theta .gt. 0.d0) plumi = -1.d0
      cons2 = 0.5d0*(root3+plumi*snphi/root3)
      cons = 0.d0
      goto 40
   30 costh = cos(theta)
      tanth = tan(theta)
      tant3 = tan(3.d0*theta)
      cons2 = costh*((1.d0+tanth*tant3)+snphi*(tant3-tanth)/root3)
      cons3 =(root3*sin(theta)+snphi*costh)/(2.d0*varj2*cos(3.d0*theta))
      goto 40
c ... Drucker - Prager
    4 cons1 = 2.d0*snphi/(root3*(3.d0-snphi))
      cons2 = 1.d0
      cons3 = 0.d0
      goto 40
    5 cons1 = snphi/(root3*sqrt(3.d0+snphi**2))
      cons2 = 1.d0
      cons3 = 0.d0
c ... vector A
   40 do 50 j = 1,4
   50 a(j) = cons1*a1(j)+cons2*a2(j)+cons3*a3(j)
c ... vector D
      xml = c14*(a(1)+a(3))
      if ( epd .eqv. .true.) xml = xml + c14*a(4)
      d(1) = 2.d0*g*(a(1)+xml)
      d(2) = g*a(2)
      d(3) = 2.d0*g*(a(3)+xml)
      d(4) = 2.d0*g*(a(4)+xml)
      if(epd .neqv. .true.) d(4) =0.d0
c ... compute part of plastic multiplier
      denom = hl
      do 80  j = 1,4
   80 denom = denom + a(j)*d(j)
      abeta = 1.d0/denom
      return
      end
      subroutine plasticst3d(a1,b1,c1,deps,eps,tx,txp,sy,hl,phi,pr,
     .          iyied,c14,g)
c **********************************************************************
c *                                                                    *
c *                                                                    *
c *   PLASTICST3D: calcula tensoes 3D, considerando plasticidade       *
c *   -----------                                                      *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ---------------------                                            *
c *                                                                    *
c *     eps(6)  - deformacao plastica + elastica                       *
c *     a1,b1,c1   - coeficientes da matriz constitutiva               *
c *     sy      - tensao de escoamento                                 *
c *     deps(6) - delta deformacao                                     *
c *     tx(6)   - tensoes anteriores                                   *
c *     tx(7)   - tensao efetiva no passo anterior                     *
c *     iyied   - numero correspondente ao criterio a ser utilizado    *
c *     txp(6)   - tensao plastica acumulada ate o passo anterior      *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   -------------------                                              *
c *                                                                    *
c *     tx(6)   - tensoes atualizadas: tx = tx + dtx                   *
c *     txp(6)  - tensoes plasticas atualizadas                        *
c *     tx(7)   - tensao efetiva atualizada                            *
c *                                                                    *
c *                                                                    *
c **********************************************************************
      implicit none
      logical epd
      integer i,ist,iyied,mstep,j
      real*8  deps(*),tx(*),eps(*),txp(*),dtxp(6),a1,b1,c1,yield,epe,g
      real*8  tx0(6),dtx(6),r,phi,snphi,pr,se
      real*8  preys,sy,hl,dste(6),steff,varj2,theta,sint3,c14,devs(6)
      real*8  se0,escur,rfac,astep,stot(6),reduc,abeta,epcon,agash,bgash
      real*8  dlamb,a(6),d(6),delep,curys,bring,str(6)
c ......................................................................
c   
c      call radian(phi)
      snphi =dsin(phi)
      epe = eps(7)
      se0 = tx(7)
c ... calculate elastic stress increment
      call stress3d(a1,b1,c1,deps,dtx)
      preys=sy+epe*hl
      dste(1) = dtx(1)
      dste(2) = dtx(2)
      dste(3) = dtx(3)
      dste(4) = dtx(4)
      dste(5) = dtx(5)
      dste(6) = dtx(6)
c ...    t0 = t:
      call aequalb(tx0,tx,6)
c ...    Tensoes lineares, t = t0 + dt
      call vsum(tx0,dste,6,str)
      if (iyied .eq. 0) goto 4
c ... evaluate effective stress
      call effst3d(se,iyied,str,devs,steff,theta,varj2,sint3,snphi,0.d0)    
c      iyied = 0 
c      if (iyied .eq. 0) goto 4
      if ( (se0-preys) .ge. 0.d0) goto 3
c ... check for yielding durint this iterations
      escur = se - preys
      if(escur .le. 0.d0) goto 4
c ... compute factor R
      rfac = escur/(se-se0)
      goto 5
    3 escur = se - se0
      if(escur .le. 0.d0) goto 4
      rfac = 1.d0
c ... evaluate number of stress reduction steps
    5 mstep = escur*8.d0/sy+1.d0
      astep = mstep
      reduc = 1.d0-rfac
      do 6 j = 1,6
         stot(j) = tx0(j)+reduc*dste(j)
    6    dste(j) = rfac*dste(j)/astep
      epcon = 0.d0
c ... loop over each stress reduction step
      do 8 ist = 1, mstep
c ... calculate vectors A and D
         call effst3d(se,iyied,stot,devs,steff,theta,varj2,sint3,snphi,
     .0.d0)
         call flow3d(a,d,devs,abeta,steff,theta,varj2,sint3,hl,
     .          snphi,c14,g,pr,0.d0,iyied)
c ... compute plastic multiplier
         agash = 0.d0
         do 9 j = 1,6
    9       agash = agash + a(j) * dste(j)
         dlamb = agash *abeta
         if ( dlamb .lt. 0.d0) dlamb = 0.d0
c ...compute elastoplastic stresses
         bgash = 0.d0
         do 10 j = 1,6
            bgash = bgash + a(j) * stot(j)
   10       stot(j) = stot(j) + dste(j)-dlamb*d(j)
c ... calculate equivalent plastic strain increment
         delep = dlamb*bgash/se
         epcon = epcon+delep
c ... calculae equivalent plastic strain
         eps(7) = eps(7) + delep
    8  continue
c ... compute effective stress
      call effst3d(se,iyied,stot,devs,steff,theta,varj2,sint3,snphi,0.0)
c ... compute equivalent yield stress
      curys = sy + eps(7)*hl
c ... scale down stresses to yield surface
      bring = 1.d0
      if ( se .gt. curys)    bring = curys/se
      do 11 j = 1,6
   11    tx(j) = bring*stot(j)
      tx(7) = bring * se
c     compute residual stresses
      do 12 j = 1, 6
         dtxp(j) = str(j) - tx(j)
   12    txp(j) = txp(j) + dtxp(j)
      goto 2
c ... compute stresses for elastic nodes or points
    4 continue
      do 13 i = 1, 6
   13    tx(i) = str(i)
      tx(7) = se
    2 continue
      return
      end
      subroutine plastic3d_new(a,b1,c1,deps,eps,tx,txp,sy1,sy2,k,hl1,
     .          hl2,phi,pr,iyied,c14,g)
c **********************************************************************
c *                                                                    *
c *                                                                    *
c *   PLASTICST3D: calcula tensoes 3D, considerando plasticidade       *
c *   -----------                                                      *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ---------------------                                            *
c *                                                                    *
c *     eps(6)  - deformacao plastica + elastica                       *
c *     a1,b1,c1   - coeficientes da matriz constitutiva               *
c *     sy      - tensao de escoamento                                 *
c *     deps(6) - delta deformacao                                     *
c *     tx(6)   - tensoes anteriores                                   *
c *     tx(7)   - tensao efetiva no passo anterior                     *
c *     iyied   - numero correspondente ao criterio a ser utilizado    *
c *     txp(6)   - tensao plastica acumulada ate o passo anterior      *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   -------------------                                              *
c *                                                                    *
c *     tx(6)   - tensoes atualizadas: tx = tx + dtx                   *
c *     txp(6)  - tensoes plasticas atualizadas                        *
c *     tx(7)   - tensao efetiva atualizada                            *
c *                                                                    *
c *                                                                    *
c **********************************************************************
      implicit none
      logical epd
      integer i,ist,iyied,mstep,j,maxit
      real*8  deps(*),tx(*),eps(*),txp(*),dtxp(6),a,b1,c1,yield,epe1,g
      real*8  tx0(6),dtx(6),r,phi,snphi,pr,se1,se2,se02,tol1,tol2
      real*8  preys,sy,hl,dste(6),steff,varj2,theta,sint3,c14,devs(6)
      real*8  se01,escur,rfac,astep,stot(6),reduc,abeta,epcon,agash
      real*8  dlamb,d(6),delep,curys,bring,str(6),k,hl1,hl2,bgash
      real*8  epe2,preys1,preys2,sy1,sy2,escur1,escur2,tol,det
      real*8  dlamb1, dlamb2, abeta1,abeta2,a1(6),d1(6),a2(6),d2(6)
      real*8  lamb, lamb1,lamb2
c ......................................................................
c   
c    
      maxit = 100
      tol = 1.0e-03
      call radian(phi)
      snphi =dsin(phi)
      epe1 = eps(7)
c      se01 = tx(7)
      preys1=sy1+epe1*hl1
c ... calculate elastic stress increment
      call stress3d(a,b1,c1,deps,dtx)
      dste(1) = dtx(1)
      dste(2) = dtx(2)
      dste(3) = dtx(3)
      dste(4) = dtx(4)
      dste(5) = dtx(5)
      dste(6) = dtx(6)
c ...    t0 = t:
      call aequalb(tx0,tx,6)
c ...    Tensoes lineares, t = t0 + dt
      call vsum(tx0,dste,6,stot)
      if ( iyied .eq. 7)   goto 5
      if (iyied .eq. 0) goto 4
      if (iyied .eq. 6)   then
         epe2 = eps(8)
         preys1 = sy2+epe2*hl2
         hl1 = hl2
      endif
c ... evaluate effective stress
      call effst3d(se1,iyied,stot,devs,steff,theta,varj2,sint3,snphi,k) 
c      if ( (se01-preys1) .ge. 0.d0) goto 3
c ... check for yielding durint this iterations
      escur1 = se1 - preys1
      if(escur1 .le. 1.0e-10) goto 4
c    3 escur1 = se1 - se01
c      if(escur .le. 0.d0) goto 4
      goto 6
    5 epe2 = eps(8)
c      se02 = tx(8)
      preys2=sy2+epe2*hl2
      call aequalb(str,stot,6)
      call effst3d(se1,5,stot,devs,steff,theta,varj2,sint3,snphi,k) 
      call effst3d(se2,6,stot,devs,steff,theta,varj2,sint3,snphi,k)  
      escur1 = se1 - preys1
      escur2 = se2 - preys2
      if (escur1 .gt. 0.d0 .and. escur2 .gt. 0.d0)  goto 15
      if (escur1 .gt. 0.d0)   iyied = 5
      if (escur2 .gt. 0.d0)   then
         iyied = 6
         preys1 = preys2
         hl1 = hl2
      endif
      if(escur1 .le. 0.d0 .and. escur2 .le. 0.d0)    goto 4
c      do 6 j = 1,6
c    6     stot(j) = tx0(j)+dste(j)
    6 continue
      lamb = 0.d0
c ... loop Newton-Raphson
      do 8 i = 1, maxit
c ... calculate vectors A and D
         call effst3d(se1,iyied,stot,devs,steff,theta,varj2,sint3,snphi
     .   ,k)
         escur1 = se1 - preys1
         tol1 = preys1*tol
         if( abs(escur1) .lt. tol1)  goto 11
         call flow3d(a1,d1,devs,abeta,steff,theta,varj2,sint3,hl1,
     .          snphi,c14,g,pr,k,iyied)
c ... compute plastic multiplier
         dlamb = escur1 *abeta
         lamb = lamb + dlamb
c ...compute elastoplastic stresses
         do 10 j = 1,6
   10       stot(j) = stot(j)-dlamb*d1(j)
    8  continue
        write(*,*) 'numero maximo de iteracoes excedido - plasticidade',
     .          maxit
        stop
   11 continue
c ... calculate equivalent plastic strain
      if (iyied .eq. 6)   then
         eps(8) = eps(8) + lamb
      else
         eps(7) = eps(7) + lamb
      endif
      goto 14
   15 continue
      lamb1 = 0.d0
      lamb2 = 0.d0
c ..  Plasticity multi-surface
c ... loop Newton-Raphson
      do 17 i = 1, maxit
c ... calculate vectors A and D
         call effst3d(se1,5,stot,devs,steff,theta,varj2,sint3,snphi
     .   ,k)
         escur1 = se1 - preys1
         call effst3d(se2,6,stot,devs,steff,theta,varj2,sint3,snphi
     .   ,k)
         escur2 = se2 - preys2
         tol1 = preys1*tol
         tol2 = preys2*tol
         if(dabs(escur1) .lt. tol1 .and. dabs(escur2) .lt. tol2) goto 19
         call flow3d(a1,d1,devs,abeta1,steff,theta,varj2,sint3,hl1,
     .          snphi,c14,g,pr,k,5)
         call flow3d(a2,d2,devs,abeta2,steff,theta,varj2,sint3,hl1,
     .          snphi,c14,g,pr,k,6)
c ... compute determinant
         det = (1.d0/abeta1)*(1.d0/abeta2)
         do j = 1, 6
            det = det - d1(j)*a2(j)*d2(j)*a1(j)
         enddo
c ... compute plastic multiplier
         dlamb1 = escur1*(1.d0/abeta2)
         dlamb2 = escur2*(1.d0/abeta1)
         do j = 1 , 6
            dlamb1 = dlamb1 - escur2*d1(j)*a2(j)
            dlamb2 = dlamb2 - escur1*d2(j)*a1(j)
         enddo
         dlamb1 = dlamb1/det
         dlamb2 = dlamb2/det
         lamb1 = lamb1 + dlamb1
         lamb2 = lamb2 + dlamb2
c ...compute elastoplastic stresses
         do 18 j = 1,6
   18       stot(j) = stot(j)-dlamb1*d1(j)-dlamb2*d2(j)
   17   continue
        write(*,*) 'numero maximo de iteracoes excedido - plasticidade 
     .       multisurface',   maxit
        stop
   19 continue  
      if (lamb1 .lt. 0.d0)   then
         iyied = 6
         preys1 = preys2
         hl1 = hl2
         call aequalb(stot,str,6)
         goto 6
      elseif (lamb2 .lt. 0.d0)   then
         iyied = 5
         call aequalb(stot,str,6)
         goto 6
      else
         eps(7) = eps(7) + lamb1
         eps(8) = eps(8) + lamb2
      endif
   14 continue
      do 16 j = 1,6
   16    tx(j) = stot(j)
      tx(7) = se1
      if (iyied .eq. 6)   tx(8) = se2
c     compute residual stresses
c      do 12 j = 1, 6
c         dtxp(j) = str(j) - tx(j)
c   12    txp(j) = txp(j) + dtxp(j)
      goto 2
c ... compute stresses for elastic nodes or points
    4 continue
      do 13 i = 1, 6
   13    tx(i) = stot(i)
      tx(7) = se1
      if (iyied .eq. 6)   tx(8) = se2
    2 continue
      return
      end
      subroutine effst3d(se,iyied,tx,devs,steff,theta,varj2,sint3,snphi,
     .k)
c **********************************************************************
c *                                                                    *
c *                                                                    *
c *   Effst: calculates the deviatoric stresses, effective stress and  *
c *          invariants                                                *
c *   -----------                                                      *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ---------------------                                            *
c *                                                                    *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   -------------------                                              *
c *                                                                    *
c *                                                                    *
c *                                                                    *
c **********************************************************************
      implicit none
      integer iyied
      real*8 root3,sm,tx(*),devs(*),varj2,varj3,steff,sint3,se
      real*8 snphi,theta,arsin,beta,alfbet,k
      root3 = 1.73205080757d0
      sm = (tx(1)+tx(2)+tx(3))/3.d0
      devs(1) = tx(1)-sm
      devs(2) = tx(2)-sm
      devs(3) = tx(3)-sm
      devs(4) = tx(4)
      devs(5) = tx(5)
      devs(6) = tx(6)
      varj2 = devs(4)*devs(4)+devs(5)*devs(5)+devs(6)*devs(6)+
     .0.5d0*(devs(1)*devs(1)+devs(2)*devs(2)+ devs(3)*devs(3))
c      varj3 = devs(4)*(devs(4)*devs(4)-varj2)
      varj3 = devs(1)*devs(2)*devs(3)+2.d0*devs(4)*devs(5)*devs(6)
     .        -devs(1)*devs(5)*devs(5)-devs(2)*devs(6)*devs(6)
     .        -devs(3)*devs(4)*devs(4)
      steff = sqrt(varj2)
      sint3 = -3.d0*root3*varj3/(2.d0*varj2*steff)
      if ( abs(sint3) .gt. 1.d0)  sint3 = sign(1.d0,sint3)
      theta = asin(sint3)/3.d0
      goto (1,2,3,4,5,6) iyied
c ... Tresca
    1 se = 2.d0*cos(theta)*steff
      return
c ... Von Mises
    2 se = root3*steff
      return
c ... Mohr-Coulomb
    3 se = sm*snphi+steff*(cos(theta)-sin(theta)*snphi/root3)
      return
c ... Drucker-Prager
    4 se = 6.*sm*snphi/(root3*(3.-snphi))+steff
      return
C ... Tension-Cut-Off
    5 se = 3.d0*sm
      return
c ... Drucker-Prager para CONCRETO
    6 beta = sqrt(3.d0)*(2.d0*k-1.d0)/k
      alfbet = (k-1)/k
      se = 3.d0*sm*AlfBet+Beta*steff
      return
      end
      subroutine flow3d(a,d,devs,abeta,steff,theta,varj2,sint3,hl,snphi,
     .                  c14,g,pr,k,iyied)
c **********************************************************************
c *                                                                    *
c *                                                                    *
c *   Effst: calculates vectors A and D                                *
c *   -----------                                                      *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ---------------------                                            *
c *                                                                    *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   -------------------                                              *
c *                                                                    *
c *                                                                    *
c *                                                                    *
c **********************************************************************
      implicit none
      logical epd
      integer iyied,j,ipl
      real*8 a1(6),a2(6),a3(6),a(6),d(6),devs(6),root3,steff,varj2
      real*8 cons1,cons2,cons3,theta,snphi,cons,costh,tanth,tant3
      real*8 plumi,c14,xml,sinth,g,abeta,sint3,hl,denom,pr,k,alfbet,beta
      root3 = 1.73205080757d0
c ... vector A1:
      a1(1) = 1.d0
      a1(2) = 1.d0
      a1(3) = 1.d0
      a1(4) = 0.d0
      a1(5) = 0.d0
      a1(6) = 0.d0
c ... vector A2:
      a2(1) = devs(1)/(2.d0*steff)
      a2(2) = devs(2)/(2.d0*steff)
      a2(3) = devs(3)/(2.d0*steff)
      a2(4) = devs(4)/steff
      a2(5) = devs(5)/steff
      a2(6) = devs(6)/steff
c ... vector A3
      a3(1) =(2.d0/3.d0)*(devs(2)*devs(3)-devs(5)*devs(5))+(1.d0/3.d0)*
     .       (devs(6)*devs(6)+devs(4)*devs(4)-devs(1)*(devs(3)+devs(2)))
      a3(2) =(2.d0/3.d0)*(devs(1)*devs(3)-devs(6)*devs(6))+(1.d0/3.d0)*
     .       (devs(5)*devs(5)+devs(4)*devs(4)-devs(2)*(devs(1)+devs(3)))
      a3(3) =(2.d0/3.d0)*(devs(1)*devs(2)-devs(4)*devs(4))+(1.d0/3.d0)*
     .       (devs(5)*devs(5)+devs(6)*devs(6)-devs(3)*(devs(1)+devs(2)))
      a3(4) = -2.d0*devs(3)*devs(4)+2*devs(5)*devs(6)
      a3(5) = -2.d0*devs(1)*devs(5)+2*devs(4)*devs(6)
      a3(6) = -2.d0*devs(2)*devs(6)+2*devs(4)*devs(5)
      goto (1,2,3,4,5,6),iyied
c ... Tresca
    1 cons1 = 0.d0
      if (dabs(theta)*57.29577951308d0 .lt. 29.d0) goto 20
      cons2 = root3
      cons3 = 0.d0
      goto 40
   20 sinth = dsin(theta)
      cons2 = 2.d0*(dcos(theta)+sinth*dtan(3.d0*theta))
      cons3 = root3*sinth/(varj2*dcos(3.d0*theta))
      goto 40
c ... Von Mises
    2 cons1 = 0.d0
      cons2 = root3
      cons3 = 0.d0
      goto 40
c ... Mohr-Coulomb
    3 cons1 = snphi/3.d0
      if (abs(theta)*57.29577951308d0 .lt. 29.d0) goto 30
      plumi = 1.d0
      if ( theta .gt. 0.d0) plumi = -1.d0
      cons2 = 0.5d0*(root3+plumi*snphi/root3)
      cons = 0.d0
      goto 40
   30 costh = cos(theta)
      tanth = tan(theta)
      tant3 = tan(3.d0*theta)
      cons2 = costh*((1.d0+tanth*tant3)+snphi*(tant3-tanth)/root3)
      cons3 =(root3*sin(theta)+snphi*costh)/(2.d0*varj2*cos(3.d0*theta))
      goto 40
c ... Drucker - Prager
    4 cons1 = 2.d0*snphi/(root3*(3.d0-snphi))
      cons2 = 1.d0
      cons3 = 0.d0
      goto 40
C ... Tension-Cut-Off
    5 cons1 = 1.d0
      cons2 = 0.d0
      cons3 = 0.d0
      goto 40
c ... Drucker - Prager para CONCRETO
    6 beta = sqrt(3.d0)*(2.d0*k-1.d0)/k
      alfbet = (k-1.d0)/k
      cons1 = AlfBet
      cons2 = Beta
      cons3 = 0.d0
c ... vector A
   40 do 50 j = 1,6
   50 a(j) = cons1*a1(j)+cons2*a2(j)+cons3*a3(j)
c ... vector D
      xml = c14*(a(1)+a(2)+a(3))
c      if ( epd .eqv. .true.) xml = xml + c14*a(4)
      d(1) = 2.d0*g*(a(1)+xml)
      d(2) = 2.d0*g*(a(2)+xml)
      d(3) = 2.d0*g*(a(3)+xml)
      d(4) = g*a(4)
      d(5) = g*a(5)
      d(6) = g*a(6)
c      if(epd .neqv. .true.) d(4) =0.
c ... compute part of plastic multiplier
      denom = hl
      do 80  j = 1,6
   80 denom = denom + a(j)*d(j)
      abeta = 1.d0/denom
      return
      end      
      subroutine plasticity(deps,a,b,c,mu,lamb,yield,phi,coes,nsteps,
     .                       nen,tx,dtx,flag)
c **********************************************************************
c *                                                                    *
c *                                                                    *
c *   PLASTICST3D: calcula tensoes 3D, considerando plasticidade       *
c *   -----------                                                      *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ---------------------                                            *
c *                                                                    *
c *     deps(6) - incrementos de deformacoes                           *
c *     a,b,c   - coeficientes da matriz constitutiva linear           *
c *     mu,lamb - constantes de Lame                                   *
c *     yield   - tensao de escoamento                                 *
c *     phi     - angulo de atrito (Drucker-Prager)                    *
c *     coes    - coesao (Drucker-Prager                               *
c *     nsteps  - numero de sub-incrementos                            *
c *     nen     - numero de nos do elemento                            *
c *     tx(6)   - tensoes do passo anterior                            *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   -------------------                                              *
c *                                                                    *
c *     tx(6)   - tensoes atualizadas: tx = tx + dtx                   *
c *                                                                    *
c *                                                                    *
c **********************************************************************
      implicit none
      integer nsteps,nen,flag,vonmises,druckerprager,test0,test,i
      real*8  deps(*),tx(*),a,b,c,mu,lamb,yield,phi,coes
      real*8  ddeps(6),tx0(6),dtx(6),r,re,ratio0,ratio,ratiodp,d,rdp
c ......................................................................
c
      flag = 0
c
c ...    dtx = D.deps (incrementos lineares de tensao):
c
      call stress3d(a,b,c,deps,dtx)
c
c ...    tx0 = tx + dtx (tensoes candidatas):
c
      call vsum(tx,dtx,6,tx0)
c
c ...    Verifica se F <= 0, para tx0:
c
      test = druckerprager(tx0,phi,coes,mu,lamb)
      if(test .eq. 0) then
        call aequalb(tx,tx0,6)         
      elseif(test .eq. 1) then
        call aequalb(tx,tx0,6)
        flag = 1
      elseif(test .eq. 2) then
        call returnmap(tx0,deps,a,b,c,mu,lamb,phi,coes,dtx)
        call vsum(tx,dtx,6,tx)
        flag = 1
      elseif(test .eq. 3) then
        dtx(1) = tx0(1)-tx(1)
        dtx(2) = tx0(2)-tx(2)
        dtx(3) = tx0(3)-tx(3)
        dtx(4) = tx0(4)-tx(4)
        dtx(5) = tx0(5)-tx(5)
        dtx(6) = tx0(6)-tx(6)
        call aequalb(tx,tx0,6)
        flag = 0
      endif
      return
      end
      subroutine dtang2d(t,deps,a,b,c,dt,epd)           
c **********************************************************************
c *                                                                    *
c *                                                                    *
c *   DTANG3D: Calcula os incrementos de tensao atraves da             *
c *            matriz D tangente                                       *
c *   -------                                                          *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ---------------------                                            *
c *                                                                    *
c *     t(6)    - tensoes                                              *
c *     deps(6) - incrementos de deformacao                            *
c *     a,b,c   - coeficientes da matriz constitutiva elastica         *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   -------------------                                              *
c *                                                                    *
c *     dt(6)  - incrementos de tensao                                 *
c *                                                                    *
c *                                                                    *
c **********************************************************************
      implicit none
      logical epd
      real*8 ZERO
      parameter (ZERO = 1.d-14)
      real*8 t(4),deps(4),a,b,c,dt(4)
      real*8 tm,sx,sy,sz,j2,d11,d12,d13,d14,d22,d23,d24,d33,d34,d44
      real*8 alpha,beta,chi,H,HD,H1,d1,d2,d3,hc
      real*8 DF1,DF2,DF3,DF4,DJ31,DJ32,DJ33
c ......................................................................
      tm = (t(1)+t(2)+t(3)) / 3.d0
      sx = t(1) - tm
      sy = t(2) - tm
      sz = t(3) - tm
      j2 = (sx*sx+sy*sy+sz*sz)/2.d0 + t(4)*t(4)
      if (j2 .le. ZERO) then
         dt(1) = 0.d0
         dt(2) = 0.d0
         dt(3) = 0.d0
         dt(4) = 0.d0
         return
      endif
c ......................................................................
c
c ... Von Mises:
c
      alpha  = 0.d0
c     Verificar qual das duas linhas esta correta:
      beta = dsqrt(3.d0) / (2.d0*dsqrt(j2))
c     beta = dsqrt(3.d0) / dsqrt(j2)
      chi  = 0.d0
c
c ... Hardening:
c
      HD    = 0.d0
c ......................................................................
c
c ... Derivada da superficie de escoamento em relacao as tensoes:
c
      DJ31 = (2*sy*sz - sx*sz - sx*sy)/3.d0
      DJ32 = (2*sx*sz - sx*sy - sy*sz)/3.d0
      DJ33 = (2*sx*sy - sx*sz - sy*sz)/3.d0
c ......................................................................
      DF1 = alpha/3.d0 + beta*sx + chi*DJ31
      DF2 = alpha/3.d0 + beta*sy + chi*DJ32
      DF3 = alpha/3.d0 + beta*sz + chi*DJ32
      DF4 = 2.d0*beta*t(4)
c ......................................................................
      H1 = DF1*DF1*a + DF1*DF2*a*b + DF1*DF3*a*b + DF1*DF2*a*b +
     .     DF2*DF2*a + DF2*DF3*a*b + DF3*DF1*a*b + DF3*DF2*a*b +
     .     DF3*DF3*a + DF4*DF4*a*c
      H =  a*a /(HD + H1)
c ......................................................................
      d1 = (DF1   + DF2*b + DF3*b)
      d2 = (DF1*b + DF2   + DF3*b)
      d3 = (DF1*b + DF2*b + DF3  )
      hc = c*H
c ......................................................................
      d11 =   a - d1*d1*H
      d12 = a*b - d1*d2*H
      d13 = a*b - d1*d3*H
      d14 = -d1*DF4*hc
c ......................................................................
      d22 = a   - d2*d2*H
      d23 = a*b - d2*d3*H
      d24 = -d2*DF4*hc
c ......................................................................
      d33 = a - d3*d3*H
      d34 = -d3*DF4*hc
c ......................................................................
      d44 = (a - DF4*DF4*hc)*c
c ......................................................................
c
c ... Incrementos de tensoes:
c
      if (epd) then
c ...    Estado plano de deformacao:
         dt(1) = d11*deps(1)+d12*deps(2)+d13*deps(3)+d14*deps(4)
         dt(2) = d12*deps(1)+d22*deps(2)+d23*deps(3)+d24*deps(4)
         dt(3) = d13*deps(1)+d23*deps(2)+d33*deps(3)+d34*deps(4)
         dt(4) = d14*deps(1)+d24*deps(2)+d34*deps(3)+d44*deps(4)
      else
c ...    Estado plano de tensao:    
         deps(3) = -(d13*deps(1)+d23*deps(2)+d34*deps(4))/d33
         dt(1) = d11*deps(1)+d12*deps(2)+d13*deps(3)+d14*deps(4)
         dt(2) = d12*deps(1)+d22*deps(2)+d23*deps(3)+d24*deps(4)
         dt(3) = 0.d0
         dt(4) = d14*deps(1)+d24*deps(2)+d34*deps(3)+d44*deps(4)
      endif
      return
      end
      subroutine matrizd(a,b,c,t,phi,d,flag)               
c **********************************************************************
c *                                                                    *
c *                                                                    *
c *   DTANG3D: Calcula os incrementos de tensao atraves da             *
c *            matriz D tangente                                       *
c *   -------                                                          *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ---------------------                                            *
c *                                                                    *
c *     t(6)    - tensoes                                              *
c *     deps(6) - incrementos de deformacao                            *
c *     a,b,c   - coeficientes da matriz constitutiva elastica         *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   -------------------                                              *
c *                                                                    *
c *     dt(6)  - incrementos de tensao                                 *
c *                                                                    *
c *                                                                    *
c **********************************************************************
      implicit none
      integer flag
      real*8 ZERO
      parameter (ZERO = 1.d-14)
      real*8 t(6),d(6,6),a,b,c
      real*8 tm,sx,sy,sz,j2
      real*8 alpha,beta,chi,H,HD,H1,d1,d2,d3,hc
      real*8 DF1,DF2,DF3,DF4,DF5,DF6,DJ31,DJ32,DJ33,phi,sinphi
c ......................................................................
      tm = (t(1)+t(2)+t(3)) / 3.d0
      sx = t(1) - tm
      sy = t(2) - tm
      sz = t(3) - tm
      j2 = (sx*sx+sy*sy+sz*sz)/2.d0 + t(4)*t(4) + t(5)*t(5) + t(6)*t(6)
      if (flag .eq. 0 .or. dsqrt(j2) .lt. 1.d-10) then
         d(1,1) = a
         d(1,2) = b*a
         d(1,3) = b*a
         d(1,4) = 0.d0
         d(1,5) = 0.d0
         d(1,6) = 0.d0
c ..........................
         d(2,1) = d(1,2)
         d(2,2) = a
         d(2,3) = b*a
         d(2,4) = 0.d0
         d(2,5) = 0.d0
         d(2,6) = 0.d0
c ..........................
         d(3,1) = d(1,3)
         d(3,2) = d(2,3)
         d(3,3) = a
         d(3,4) = 0.d0
         d(3,5) = 0.d0
         d(3,6) = 0.d0
c ..........................
         d(4,1) = d(1,4)
         d(4,2) = d(2,4)
         d(4,3) = d(3,4)
         d(4,4) = c*a
         d(4,5) = 0.d0
         d(4,6) = 0.d0
c ..........................
         d(5,1) = d(1,5)
         d(5,2) = d(2,5)
         d(5,3) = d(3,5)
         d(5,4) = d(4,5)
         d(5,5) = c*a
         d(5,6) = 0.d0
c ..........................
         d(6,1) = d(1,6)
         d(6,2) = d(2,6)
         d(6,3) = d(3,6)
         d(6,4) = d(4,6)
         d(6,5) = d(5,6)
         d(6,6) = c*a
c ..........................
         return
      endif
      tm = (t(1)+t(2)+t(3)) / 3.d0
      sx = t(1) - tm
      sy = t(2) - tm
      sz = t(3) - tm
      j2 = (sx*sx+sy*sy+sz*sz)/2.d0 + t(4)*t(4) + t(5)*t(5) + t(6)*t(6)
      if (j2 .eq. 0.d0) then
         print*, 'J2 = 0 !'
c         stop
      endif
c ......................................................................
c
c ... Von Mises:
c
c      alpha  = 0.d0
c     Verificar qual das duas linhas esta correta:
c      beta = dsqrt(3.d0) / (2.d0*dsqrt(j2))
c      beta = dsqrt(3.d0) / dsqrt(j2)
c ......................................................................
c
c ... Drucker-Prager:
c
      sinphi = dsin(phi)
      alpha = 3.d0*(sinphi/dsqrt(9.d0+3.d0*sinphi*sinphi))
      beta  = 1.d0 / (2.d0*dsqrt(j2))
c .......................................................................
      chi  = 0.d0
c
c ... Hardening:
c
      HD    = 0.d0
c ......................................................................
c
c ... Derivada da superficie de escoamento em relacao as tensoes:
c
      DJ31 = (2*sy*sz - sx*sz - sx*sy)/3.d0
      DJ32 = (2*sx*sz - sx*sy - sy*sz)/3.d0
      DJ33 = (2*sx*sy - sx*sz - sy*sz)/3.d0
c ......................................................................
      DF1 = alpha/3.d0 + beta*sx + chi*DJ31
      DF2 = alpha/3.d0 + beta*sy + chi*DJ32
      DF3 = alpha/3.d0 + beta*sz + chi*DJ32
      DF4 = 2.d0*beta*t(4)
      DF5 = 2.d0*beta*t(5)
      DF6 = 2.d0*beta*t(6)
c ......................................................................
      H1 = DF1*DF1*a + DF1*DF2*a*b + DF1*DF3*a*b + DF1*DF2*a*b +
     .     DF2*DF2*a + DF2*DF3*a*b + DF3*DF1*a*b + DF3*DF2*a*b +
     .     DF3*DF3*a + DF4*DF4*a*c + DF5*DF5*a*c + DF6*DF6*a*c
      H =  a*a /(HD + H1)
c ......................................................................
      d1 = (DF1   + DF2*b + DF3*b)
      d2 = (DF1*b + DF2   + DF3*b)
      d3 = (DF1*b + DF2*b + DF3  )
      hc = c*H
c ......................................................................
      d(1,1) =   a - d1*d1*H
      d(1,2) = a*b - d1*d2*H
      d(1,3) = a*b - d1*d3*H
      d(1,4) = -d1*DF4*hc
      d(1,5) = -d1*DF5*hc
      d(1,6) = -d1*DF6*hc
c ......................................................................
      d(2,1) = d(1,2)
      d(2,2) = a   - d2*d2*H
      d(2,3) = a*b - d2*d3*H
      d(2,4) = -d2*DF4*hc
      d(2,5) = -d2*DF5*hc
      d(2,6) = -d2*DF6*hc
c ......................................................................
      d(3,1) = d(1,3)
      d(3,2) = d(2,3)
      d(3,3) = a - d3*d3*H
      d(3,4) = -d3*DF4*hc
      d(3,5) = -d3*DF5*hc
      d(3,6) = -d3*DF6*hc
c ......................................................................
      d(4,1) = d(1,4)
      d(4,2) = d(2,4)
      d(4,3) = d(3,4)
      d(4,4) = (a - DF4*DF4*hc)*c
      d(4,5) = -DF4*DF5*c*hc
      d(4,6) = -DF4*DF6*c*hc
c ......................................................................
      d(5,1) = d(1,5)
      d(5,2) = d(2,5)
      d(5,3) = d(3,5)
      d(5,4) = d(4,5)
      d(5,5) = (a - DF5*DF5*hc)*c
      d(5,6) = -DF5*DF6*c*hc
c ......................................................................
      d(6,1) = d(1,6)
      d(6,2) = d(2,6)
      d(6,3) = d(3,6)
      d(6,4) = d(4,6)
      d(6,5) = d(5,6)
      d(6,6) = (a - DF6*DF6*hc)*c
      return
      end
      real*8 function ratio(t,dt,n,yield)               
c **********************************************************************
c *                                                                    *
c *                                                    18/04/03        *
c *   RATIO:                                                           *
c *   -----                                                            *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ---------------------                                            *
c *                                                                    *
c *     t(n)  - tensoes                                                *
c *     dt(n) - tensoes incrementais                                   *
c *     yield - tensao de escoamento                                   *
c *                                                                    *
c *   Valores de retorno:                                              *
c *   -------------------                                              *
c *                                                                    *
c *     ratio                                                          *
c *                                                                    *
c **********************************************************************
      implicit none
      integer n
      real*8  ZERO
      parameter (ZERO = 1.d-14)
      real*8 t(*),dt(*),yield
      real*8 sx,sy,sz,tm,j2,dj2
c ......................................................................
      ratio = 0.d0
      tm = (t(1)+t(2)+t(3)) / 3.d0
      sx = t(1) - tm
      sy = t(2) - tm
      sz = t(3) - tm
      j2 = (sx*sx+sy*sy+sz*sz)/2.d0 + t(4)*t(4)
c ......................................................................
      tm  = (dt(1)+dt(2)+dt(3)) / 3.d0
      sx  = dt(1) - tm
      sy  = dt(2) - tm
      sz  = dt(3) - tm
      dj2 = (sx*sx+sy*sy+sz*sz)/2.d0 + dt(4)*dt(4)
      if(n .eq. 6) then
         j2  =  j2 +  t(5)*t(5) + t(6)*t(6)
         dj2 = dj2 + dt(5)*dt(5)+dt(6)*dt(6)
      endif
      if (dj2 .le. ZERO) return
c ......................................................................
      ratio = (yield/dsqrt(3.d0)-dsqrt(j2)) / dsqrt(dj2)
c ......................................................................
      return
      end
      subroutine deform2d0(hx,hy,u,b,eps,nen)
c **********************************************************************
c *                                                                    *
c *                                                                    *
c *   DEFORM2D_EPD: calcula deformacoes de estado plano de tensao      *
c *   ------------                                                     *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ---------------------                                            *
c *                                                                    *
c *     hx(nen)   - derivadas de h em relacao a x                      *
c *     hy(nen)   - derivadas de h em relacao a y                      *
c *     u(nst)    - deslocamentos                                      *
c *     nen       - numero de nos do elemento                          *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   -------------------                                              *
c *                                                                    *
c *     eps(4) - deformacoes                                           *
c *                                                                    *
c *                                                                    *
c **********************************************************************
      implicit none
      integer nen,i,j1,j2,j3
      real*8  hx(*),hy(*),u(*),eps(*),b
      real*8  ux,uy,vx,vy
c ......................................................................
      ux = 0.d0
      uy = 0.d0
      vx = 0.d0
      vy = 0.d0
      do 100 i = 1, nen
         j1 = (i-1)*2+1
         j2 = j1 + 1
         j3 = j2 + 1
         ux = ux + hx(i)*u(j1)
         uy = uy + hy(i)*u(j1)
         vx = vx + hx(i)*u(j2)
         vy = vy + hy(i)*u(j2)
  100 continue
      eps(1) = ux
      eps(2) = vy
      eps(3) = -b*(eps(1)+eps(2))
      eps(4) = uy + vx
      return
      end
      subroutine stress2d0(a,b,c,eps,t)
c **********************************************************************
c *                                                                    *
c *                                                                    *
c *   STRESS2D_EPD: calcula tensoes                                    *
c *   ------------                                                     *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ---------------------                                            *
c *                                                                    *
c *     a,b,c  - coeficientes da matriz constitutiva                   *
c *     eps(4) - deformacoes                                           *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   -------------------                                              *
c *                                                                    *
c *     t(4) - tensoes                                                 *
c *                                                                    *
c *                                                                    *
c **********************************************************************
      implicit none
      real*8 a,b,c,eps(*),t(*)
c ......................................................................
      t(1) = (  eps(1) + b*eps(2) + b*eps(3))*a
      t(2) = (b*eps(1) +   eps(2) + b*eps(3))*a
      t(3) = (b*eps(1) + b*eps(2) +   eps(3))*a
      t(4) = eps(4)*c*a
      return
      end
      real*8 function effectivet(t,phi,coes)               
c **********************************************************************
c *                                                                    *
c *   INVAR2:                                                          *
c *   ------                                                           *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ---------------------                                            *
c *                                                                    *
c *     t(n)   - tensoes                                               *
c *                                                                    *
c *   Valores de retorno:                                              *
c *   -------------------                                              *
c *                                                                    *
c *     invar2 - segundo invariante do tensor desviatorio              *
c *                                                                    *
c **********************************************************************
      implicit none
      real*8  t(*),phi,coes,sx,sy,sz,tm,j2,sinphi,alpha,tb,K,raiz
c ......................................................................
      tm = (t(1)+t(2)+t(3)) / 3.d0
      sx = t(1) - tm
      sy = t(2) - tm
      sz = t(3) - tm
      j2 = (sx*sx+sy*sy+sz*sz)/2.d0 + t(4)*t(4)+ t(5)*t(5)+t(6)*t(6)
      tb = dsqrt(j2)
      sinphi = dsin(phi)
c      alpha = sinphi/dsqrt(9.d0+3.d0*sinphi*sinphi)
      raiz = dsqrt(3.d0)*(3.d0-sinphi)
      alpha  = 2.d0*sinphi/raiz
c      K = (3.d0*coes*dcos(phi))/dsqrt(9.d0+3.d0*sinphi*sinphi)
      K = (6.d0*coes*dcos(phi))/raiz
      effectivet =  (3.d0*alpha*tm + tb)/K
c ......................................................................
      return
      end
      real*8 function invar2(t,n)               
c **********************************************************************
c *                                                                    *
c *   INVAR2:                                                          *
c *   ------                                                           *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ---------------------                                            *
c *                                                                    *
c *     t(n)   - tensoes                                               *
c *                                                                    *
c *   Valores de retorno:                                              *
c *   -------------------                                              *
c *                                                                    *
c *     invar2 - segundo invariante do tensor desviatorio              *
c *                                                                    *
c **********************************************************************
      implicit none
      integer n
      real*8  t(*),sx,sy,sz,tm
c ......................................................................
      tm = (t(1)+t(2)+t(3)) / 3.d0
      sx = t(1) - tm
      sy = t(2) - tm
      sz = t(3) - tm
      invar2 = (sx*sx+sy*sy+sz*sz)/2.d0 + t(4)*t(4)
      if(n .eq. 6) then
         invar2 = invar2 + t(5)*t(5)+t(6)*t(6)
      endif
c ......................................................................
      return
      end
      integer function vonmises(t,n,yield)               
c **********************************************************************
c *                                                                    *
c *                                                    18/04/03        *
c *   VONMISES: Criterio de plastificacao de Von Mises                 *
c *   --------                                                         *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ---------------------                                            *
c *                                                                    *
c *     t(n) - tensoes                                                 *
c *     yield  - tensao de escoamento                                  *
c *                                                                    *
c *   Valores de retorno:                                              *
c *   -------------------                                              *
c *                                                                    *
c *     vnomises = 0 (fase elastica linear)                            *
c *     vonmises = 1 (fase plastica)                                   *
c *                                                                    *
c **********************************************************************
      implicit none
      integer n
      real*8 t(*),yield
      real*8 sx,sy,sz,tm,j2
c ......................................................................
      vonmises = 0
      tm = (t(1)+t(2)+t(3)) / 3.d0
      sx = t(1) - tm
      sy = t(2) - tm
      sz = t(3) - tm
      j2 = (sx*sx+sy*sy+sz*sz)/2.d0 + t(4)*t(4)
      if(n .eq. 6) then
         j2 = j2 + t(5)*t(5) + t(6)*t(6)
      endif
      tm = dsqrt(j2)*dsqrt(3.d0)
      if (tm .gt. yield) then
         vonmises = 1
         return
      endif
c ......................................................................
      return
      end
      integer function druckerprager(t,phi,c,mu,lamb)               
c **********************************************************************
c *                                                                    *
c *                                                    18/04/03        *
c *   VONMISES: Criterio de plastificacao de Von Mises                 *
c *   --------                                                         *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ---------------------                                            *
c *                                                                    *
c *     t(n) - tensoes                                                 *
c *     yield  - tensao de escoamento                                  *
c *                                                                    *
c *   Valores de retorno:                                              *
c *   -------------------                                              *
c *                                                                    *
c *     vnomises = 0 (fase elastica linear)                            *
c *     vonmises = 1 (fase plastica)                                   *
c *                                                                    *
c **********************************************************************
      implicit none
      real*8  ZERO
      real*8 t(*),c,phi,mu,lamb
      real*8 sx,sy,sz,tm,j2,tb,alpha,K,F,F1,sinphi,raiz
      parameter (ZERO = 1.d-12)
c ......................................................................
      druckerprager = 0
      tm = (t(1)+t(2)+t(3)) / 3.d0
      sx = t(1) - tm
      sy = t(2) - tm
      sz = t(3) - tm
      j2 = (sx*sx+sy*sy+sz*sz)/2.d0 + t(4)*t(4)+ t(5)*t(5) + t(6)*t(6)
      tb = dsqrt(j2)
      sinphi = dsin(phi)
      raiz   = dsqrt(3.d0)*(3.d0-sinphi)
      alpha  = 2.d0*sinphi/raiz
c      alpha = sinphi/dsqrt(9.d0+3.d0*sinphi*sinphi)
c      K = (3.d0*c*dcos(phi))/dsqrt(9.d0+3.d0*sinphi*sinphi)
      K = (6.d0*c*dcos(phi))/raiz
      F =  3.d0*alpha*tm + tb - K
      if (F .lt. 0.d0) return
      druckerprager = 1
      if (F .le. ZERO) return
      druckerprager = 2
      F1 = -3.d0*tm + alpha*tb*(1.d0+3.d0*lamb/(2.d0*mu))+K/alpha
      if (F1 .lt. 0.d0) then
         druckerprager = 3
         t(1) = K/(3.d0*alpha)
         t(2) = K/(3.d0*alpha)
         t(3) = K/(3.d0*alpha)
         t(4) = 0.d0
         t(5) = 0.d0
         t(6) = 0.d0
      endif
c ......................................................................
      return
      end
      subroutine returnmap(t,deps,a,b,c,mu,lamb,phi,coes,dt)         
c **********************************************************************
c *                                                                    *
c *                                                                    *
c *   DTANG3D: Calcula os incrementos de tensao atraves da             *
c *            matriz D tangente                                       *
c *   -------                                                          *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ---------------------                                            *
c *                                                                    *
c *     t(6)    - tensoes                                              *
c *     deps(6) - incrementos de deformacao                            *
c *     a,b,c   - coeficientes da matriz constitutiva elastica         *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   -------------------                                              *
c *                                                                    *
c *     dt(6)  - incrementos de tensao                                 *
c *                                                                    *
c *                                                                    *
c **********************************************************************
      implicit none
      integer test,druckerprag0,druckerprager
      real*8 t(6),deps(6),a,b,c,mu,lamb,phi,coes,dt(6),dp(6)
      real*8 tm,sx,sy,sz,j2
      real*8 alpha,sinphi,F,K,rj2,raiz
c ......................................................................
      sinphi = dsin(phi)
c      alpha  = sinphi/dsqrt(9.d0+3.d0*sinphi*sinphi)
      raiz = dsqrt(3.d0)*(3.d0-sinphi)
      alpha  = 2.d0*sinphi/raiz
      tm = (t(1)+t(2)+t(3)) / 3.d0
      sx = t(1) - tm
      sy = t(2) - tm
      sz = t(3) - tm
      j2 = (sx*sx+sy*sy+sz*sz)/2.d0 + t(4)*t(4) + t(5)*t(5) + t(6)*t(6)
      if(j2 .le. 0.d0) then
         print*,'J2 <= 0 !'
         stop
      endif
c      K = (3.d0*coes*dcos(phi))/dsqrt(9.d0+3.d0*sinphi*sinphi)
      K = (6.d0*coes*dcos(phi))/raiz
      F = (3.d0*alpha*tm + dsqrt(j2) -K) / 
     .    (mu + 3.d0*alpha*alpha*(3.d0*lamb+2.d0*mu))
      rj2 = 1.d0/(2.d0*dsqrt(j2))
      dp(1) = deps(1) - F*(alpha + (2.d0*t(1)-t(2)-t(3))*rj2/3.d0)
      dp(2) = deps(2) - F*(alpha + (2.d0*t(2)-t(1)-t(3))*rj2/3.d0)
      dp(3) = deps(3) - F*(alpha + (2.d0*t(3)-t(1)-t(2))*rj2/3.d0)
      dp(4) = deps(4) - F*(t(4)*rj2*2.d0)
      dp(5) = deps(5) - F*(t(5)*rj2*2.d0)
      dp(6) = deps(6) - F*(t(6)*rj2*2.d0)
      call stress3d(a,b,c,dp,dt)
      return
      end
      integer function druckerprager3d(t,phi,fc)               
c **********************************************************************
c *                                                                    *
c *                                                    18/04/03        *
c *   DRUCKERPRAGER3d: Criterio de plastificacao de Drucker-Prager     *
c *   --------                                                         *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ---------------------                                            *
c *                                                                    *
c *     mu,lamb - constantes de Lame                                   *
c *     phi     - angulo de atrito (Drucker-Prager)                    *
c *     fc     - tensao de compressao                                  *
c *     c       - coesao (Drucker-Prager                               * 
c *                                                                    *
c *   Valores de retorno:                                              *
c *   -------------------                                              *
c *                                                                    *
c *     druckerprager = 0 (fase elastica linear)                       *
c *     druckerprager = 1 (fase plastica)                              *
c *                                                                    *
c **********************************************************************
      implicit none
      real*8  ZERO
      real*8 t(*),c,phi
      real*8 sx,sy,sz,tm,j2,tb,alpha,K,F,F1,sinphi,raiz,fc
      parameter (ZERO = 1.d-12)
c ......................................................................
      druckerprager3d = 0
      call radian(phi)
      c = fc*(1-dsin(phi))/(2*dcos(phi))
      tm = (t(1)+t(2)+t(3)) / 3.d0
      sx = t(1) - tm
      sy = t(2) - tm
      sz = t(3) - tm
      j2 = (sx*sx+sy*sy+sz*sz)/2.d0 + t(4)*t(4)+ t(5)*t(5) + t(6)*t(6)
      tb = dsqrt(j2)
      sinphi = dsin(phi)
      alpha = sinphi/dsqrt(9.d0+3.d0*sinphi*sinphi)
      K = (3.d0*c*dcos(phi))/dsqrt(9.d0+3.d0*sinphi*sinphi)
      F =  3.d0*alpha*tm + tb - K
      if (F .ge. 0.d0) then
          druckerprager3d = 1
      endif
c ......................................................................
      return
      end
      integer function rankine(sigma,yield)
c **********************************************************************
c *                                                                    *
c *                                                    18/04/03        *
c *   RANKINE : Criterio de Rankine                                    *
c *   --------                                                         *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ---------------------                                            *
c *     yield - tensao de escoamento                                   *
c *     sigma - tensao principal sigma1                                *
c *                                                                    *
c *   Valores de retorno:                                              *
c *   -------------------                                              *
c *                                                                    *
c *     rankine = 0 (fase elastica linear)                             *
c *     rankine = 1 (fase plastica)                                    *
c *                                                                    *
c **********************************************************************
      implicit none
      real*8  sigma,yield
      if (sigma .lt. yield) rankine = 0
      if (sigma .ge. yield) rankine = 1
      return
      end
      subroutine principal_stress(t,ntn)
c **********************************************************************
c *                                                                    *
c *   PSTRESS3D: tensoes principais                                    *
c *   ---------                                                        *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ----------------------                                           *
c *                                                                    *
c *    t(ntn,nnode)  - tensoes nodais referenciadas a xyz              *
c *    nnode - numero de nos                                           *
c *    ntn   - numero de tensoes por no                                *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   --------------------                                             *
c *                                                                    *
c *    t(ntn,nnode)  - tensoes principais                              *
c *    v1(3,nnode) - direcao principal 1                               *
c *    v2(3,nnode) - direcao principal 2                               *
c *    v3(3,nnode) - direcao principal 3                               *
c *                                                                    *
c **********************************************************************
      implicit none
      integer nnode,ntn,i,ip(3)
      real*8  t(ntn)
      real*8  a(3,3),x(3,3),w(3)
c ......................................................................
      a(1,1) = t(1)
      a(1,2) = t(4)
      a(1,3) = 0.d0
      a(2,1) = a(1,2)
      a(2,2) = t(2)
      a(2,3) = 0.d0
      a(3,1) = a(1,3)
      a(3,2) = a(2,3)
      a(3,3) = t(3)
      if(ntn .eq. 6) then
        a(1,3) = t(6)
        a(2,3) = t(5)
        a(3,1) = a(1,3)
        a(3,2) = a(2,3)
      endif
      call jacobi(a,x,w,ip,3,1.d-12,100)
      t(1) = w(1)
      t(2) = w(2)
      t(3) = w(3)
   20 continue
      return
      end
      subroutine radian(phi)
c **********************************************************************
c *                                                                    *
c *   RADIAN : converte angulo para radianos                           *
c *   ---------                                                        *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ----------------------                                           *
c *                                                                    *
c *    phi - angulo em graus                                           *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   --------------------                                             *
c *                                                                    *
c *    phi - angulo em radianos                                        *
c *                                                                    *
c **********************************************************************      
      implicit none
      real*8 phi, pi, temp
      pi = 4*atan(1.d0)
      temp = phi*pi/180.d0
      phi = temp
      return
      end
