C     ----------------------------------------------------------------
      program example
C     ----------------------------------------------------------------
C--   12-03-07  Basic QCDNUM example job
C     ----------------------------------------------------------------
C--   Michiel Botje   h24@nikhef.nl
      implicit double precision (a-h,o-z)
      data as0/0.118/, r20/8315.25D0/, iord/1/, nfin/0/ !alphas, NNLO, VFNS
      external func                                !input parton dists
      dimension def(-6:6,12)                       !flavor composition
      data def  /
C--   tb  bb  cb  sb  ub  db   g   d   u   s   c   b   t
C--   -6  -5  -4  -3  -2  -1   0   1   2   3   4   5   6  
     + 0., 0., 0., 0., 0.,-1., 0., 1., 0., 0., 0., 0., 0.,   !dval
     + 0., 0., 0., 0.,-1., 0., 0., 0., 1., 0., 0., 0., 0.,   !uval
     + 0., 0., 0.,-1., 0., 0., 0., 0., 0., 1., 0., 0., 0.,   !sval
     + 0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0.,   !dbar
     + 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0.,   !ubar
     + 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0.,   !sbar
     + 0., 0.,-1., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0.,   !cval
     + 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,   !cbar
     + 52*0.    /
      data xmin/1.D-5/, nxin/540/, iosp/2/            !x grid, splord
      dimension qq(2),wt(2)                           !mu2 grid
      data qq/2.0D0,1.0D8/, wt/1.D0,1.D0/, nqin/240/     !mu2 grid
      data q2c/3.0D0/, q2b/25.0D0/, q2t/28900.0D0/, q0/2.0/  !thresh and mu20
      data x/1.D-3/, q/4.8D0/, qmz2/8315.25D0/         !output scales
      dimension pdf(-6:6)                             !pdfout


!      call InitPDFsetByName('cteq5l.LHgrid')


C--   lun = 6 stdout, -6 stdout w/out banner page
      lun    = 6
      lunout = abs(lun)      
      call qcinit(lun,' ')                      !initialize
      call setval('elim',-1.d0)
      call gxmake(xmin,1,1,nxin,nx,iosp)        !x-grid
      call gqmake(qq,wt,2,nqin,nq)              !mu2-grid
      call setord(iord)                         !LO, NLO, NNLO
      call fillwt(1,id1,id2,nw)                 !calculate weights
      call setord(iord)                         !LO, NLO, NNLO
      call setalf(as0,r20)                      !input alphas
      iqc  = iqfrmq(q2c)                        !charm threshold
      iqb  = iqfrmq(q2b)                        !bottom threshold
      iqt  = iqfrmq(q2t)                        !bottom threshold
      call setcbt(nfin,iqc,iqb,iqt)             !thesholds in the vfns
      iq0  = iqfrmq(q0)                         !starting scale
      call evolfg(1,func,def,iq0,eps)           !evolve all pdf's

!      xx = 0.5d0
!      do qmy=2.d0, 40.d0, 1.d0
!      call fpdfxq(1,xx,qmy,pdf,1)                  !interpolate all pdf's
!      write(*,*) 'pdfC',qmy,pdf(4),pdf(5)
!      enddo

!      write(*,*)'alpha', asfunc(2.69d0,nfout,ierr)
!      write(*,*)'nfout', nfout
!       write(*,*)'alphaMz', asfunc(8315.25D0,nfout,ierr)
!       write(*,*)'nfout', nfout
!       write(*,*)'alpha100', asfunc(100.0D0,nfout,ierr)
!       write(*,*)'nfout', nfout
!       write(*,*)'alpha20', asfunc(20.00d0,nfout,ierr)
!       write(*,*)'nfout', nfout
!       write(*,*)'alpha2', asfunc(2.00d0,nfout,ierr)
!       write(*,*)'nfout', nfout
      
!      call GETCBT( nfixO, q2cO, q2bO, q2tO )
!      write(*,*) nfixO, q2cO, q2bO, q2tO
!      write(*,*)'alpha', asfunc(q2cO,nfout,ierr)
!      write(*,*)'nfout', nfout

!      return

!      xx = 120.d0/14000.d0
       call printTable(2.d0)
       call printTable(1.d1)
       call printTable(1.d3)
       call printTable(1.d5)
       call printTable(1.d7)
       call printTable(1.d8)
!      xx = 1.d-5;
!      do qloop=2.0,20.0, 1.0
!        call fpdfxq(1,xx,qloop,pdf,1)                  !interpolate all pdf's
!        write(*,*)'gluon ',xx,qloop, pdf(0)
!      enddo

    

!      csea = 2.D0*pdf(4)                        !charm sea at x,Q2
!      asmz = asfunc(qmz2,nfout,ierr)            !alphas(mz2)
!      write(lunout,'('' x, q, CharmSea ='',3E13.5)'), x,q,csea
!      write(lunout,'('' as(mz2)        ='', E13.5)'), asmz
      end
      
C     ----------------------------------------------------------------
      subroutine printTable(q2my)
      implicit none
      double precision logX,xx,q2my,pdf,step
      integer nsteps
      dimension pdf(-6:6)                             !pdfout
      nsteps=2000

      step = 5.d0/dble(nsteps)
      write(*,*) nsteps+1, q2my
      do logX=-5,0.d0+step/2.d0, step
         xx = 10**logX
         call fpdfxq(1,xx,q2my,pdf,1)                  !interpolate all pdf's
         write(*,*) xx, pdf(-6),pdf(-5),pdf(-4),pdf(-3),
     + pdf(-2),pdf(-1),pdf(0),pdf(1),pdf(2),pdf(3),pdf(4),pdf(5),pdf(6)
      enddo

      endsubroutine

C     ======================================
      double precision function func(ipdf,x)
C     ======================================

      implicit double precision (a-h,o-z)

!      double precision xf(-6:6)
!      q0 = 1.d0

!      call evolvePDF(x,q0,xf)

!                     func = 0.D0
!      if(ipdf.eq. 0) func = xf(0)
!      if(ipdf.eq. 1) func = xf(1)-xf(-1)
!      if(ipdf.eq. 2) func = xf(2)-xf(-2)
!      if(ipdf.eq. 3) func = 0.D0
!      if(ipdf.eq. 4) func = xf(-1)
!      if(ipdf.eq. 5) func = xf(-2)
!      if(ipdf.eq. 6) func = xf(-3)
!      if(ipdf.eq. 7) func = 0.D0
!      if(ipdf.eq. 8) func = 0.D0
!      if(ipdf.eq. 9) func = 0.D0
!      if(ipdf.eq.10) func = 0.D0
!      if(ipdf.eq.11) func = 0.D0
!      if(ipdf.eq.12) func = 0.D0


                     func = 0.D0
      if(ipdf.eq. 0) func = xglu(x)
      if(ipdf.eq. 1) func = xdnv(x)
      if(ipdf.eq. 2) func = xupv(x)
      if(ipdf.eq. 3) func = 0.D0
      if(ipdf.eq. 4) func = xdbar(x)
      if(ipdf.eq. 5) func = xubar(x)
      if(ipdf.eq. 6) func = xsbar(x)
      if(ipdf.eq. 7) func = 0.D0
      if(ipdf.eq. 8) func = 0.D0
      if(ipdf.eq. 9) func = 0.D0
      if(ipdf.eq.10) func = 0.D0
      if(ipdf.eq.11) func = 0.D0
      if(ipdf.eq.12) func = 0.D0










      return
      end
 
C     =================================
      double precision function xupv(x)
C     =================================

      implicit double precision (a-h,o-z)

      data au /5.107200D0/

      xupv = au * x**0.8D0 * (1.D0-x)**3.D0

      return
      end

C     =================================
      double precision function xdnv(x)
C     =================================

      implicit double precision (a-h,o-z)

      data ad /3.064320D0/

      xdnv = ad * x**0.8D0 * (1.D0-x)**4.D0

      return
      end

C     =================================
      double precision function xglu(x)
C     =================================

      implicit double precision (a-h,o-z)

      data ag    /1.7D0 /

      common /msum/ glsum, uvsum, dvsum

      xglu = ag * x**(-0.1D0) * (1.D0-x)**5.D0

      return
      end

C     ==================================
      double precision function xdbar(x)
C     ==================================

      implicit double precision (a-h,o-z)

      data adbar /0.1939875D0/

      xdbar = adbar * x**(-0.1D0) * (1.D0-x)**6.D0

      return
      end

C     ==================================
      double precision function xubar(x)
C     ==================================

      implicit double precision (a-h,o-z)

      xubar = xdbar(x) * (1.D0-x)

      return
      end

C     ==================================
      double precision function xsbar(x)
C     ==================================

      implicit double precision (a-h,o-z)

      xsbar = 0.2D0 * (xdbar(x)+xubar(x))

      return
      end







