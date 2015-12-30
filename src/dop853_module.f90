
module dop853_module

use dop853_constants

implicit none

    abstract interface
        subroutine func(n,x,y,f)
            implicit none
            integer,intent(in) :: n
            double precision,intent(in) :: x
            double precision,dimension(n),intent(in) :: y
            double precision,dimension(n),intent(out) :: f
        end subroutine func
    end interface

contains

! ----------------------------------------------------------
!     numerical solution of a system of first order
!     ordinary differential equations  y'=f(x,y).
!     this is an explicit runge-kutta method of order 8(5,3)
!     due to dormand & prince (with stepsize control and
!     dense output)
!
!     authors: e. hairer and g. wanner
!              universite de geneve, dept. de mathematiques
!              ch-1211 geneve 24, switzerland
!              e-mail:  ernst.hairer@unige.ch
!                       gerhard.wanner@unige.ch
!
!     this code is described in:
!         e. hairer, s.p. norsett and g. wanner, solving ordinary
!         differential equations i. nonstiff problems. 2nd edition.
!         springer series in computational mathematics,
!         springer-verlag (1993)
!
!     version of october 11, 2009
!      (new option iout=3 for sparse dense output)
!
!     input parameters
!     ----------------
!     n           dimension of the system
!
!     fcn         name (external) of subroutine computing the
!                 value of f(x,y):
!                    subroutine fcn(n,x,y,f)
!                    double precision x,y(n),f(n)
!                    f(1)=...   etc.
!
!     x           initial x-value
!
!     y(n)        initial values for y
!
!     xend        final x-value (xend-x may be positive or negative)
!
!     rtol,atol   relative and absolute error tolerances. they
!                 can be both scalars or else both vectors of length n.
!                 atol should be strictly positive (possibly very small)
!
!     itol        switch for rtol and atol:
!                   itol=0: both rtol and atol are scalars.
!                     the code keeps, roughly, the local error of
!                     y(i) below rtol*abs(y(i))+atol
!                   itol=1: both rtol and atol are vectors.
!                     the code keeps the local error of y(i) below
!                     rtol(i)*abs(y(i))+atol(i).
!
!     solout      name (external) of subroutine providing the
!                 numerical solution during integration.
!                 if iout.ge.1, it is called during integration.
!                 supply a dummy subroutine if iout=0.
!                 it must have the form
!                    subroutine solout (nr,xold,x,y,n,con,icomp,nd,
!                                       irtrn,xout)
!                    dimension y(n),con(8*nd),icomp(nd)
!                    ....
!                 solout furnishes the solution "y" at the nr-th
!                    grid-point "x" (thereby the initial value is
!                    the first grid-point).
!                 "xold" is the preceeding grid-point.
!                 "irtrn" serves to interrupt the integration. if irtrn
!                    is set <0, dop853 will return to the calling program.
!                    if the numerical solution is altered in solout,
!                    set  irtrn = 2
!                 "xout" can be used for efficient intermediate output
!                    if one puts iout=3
!                    when nr=1 define the first output point xout in solout.
!                      the subroutine solout will be called only when
!                      xout is in the interval [xold,x]; during this call
!                      a new value for xout can be defined, etc.
!
!          -----  continuous output: -----
!                 during calls to "solout", a continuous solution
!                 for the interval [xold,x] is available through
!                 the function
!                        >>>   contd8(i,s,con,icomp,nd)   <<<
!                 which provides an approximation to the i-th
!                 component of the solution at the point s. the value
!                 s should lie in the interval [xold,x].
!
!     iout        switch for calling the subroutine solout:
!                    iout=0: subroutine is never called
!                    iout=1: subroutine is called after every successful step
!                    iout=2: dense output is performed after every successful step
!                            (in this case iwork(5) must be specified)
!                    iout=3: dense output is performed in steps defined by the user
!                            (see "xout" above)
!
!     work        array of working space of length "lwork".
!                 work(1),...,work(20) serve as parameters for the code.
!                 for standard use, set them to zero before calling.
!                 "lwork" must be at least  11*n+8*nrdens+21
!                 where  nrdens = iwork(5)
!
!     lwork       declared lenght of array "work".
!
!     iwork       integer working space of lenght "liwork".
!                 iwork(1),...,iwork(20) serve as parameters for the code.
!                 for standard use, set them to zero before calling.
!                 "liwork" must be at least nrdens+21 .
!
!     liwork      declared lenght of array "iwork".
!
!-----------------------------------------------------------------------
!
!     sophisticated setting of parameters
!     -----------------------------------
!              several parameters (work(1),...,iwork(1),...) allow
!              to adapt the code to the problem and to the needs of
!              the user. for zero input, the code chooses default values.
!
!    work(2)   the safety factor in step size prediction,
!              default 0.9d0.
!
!    work(3), work(4)   parameters for step size selection
!              the new step size is chosen subject to the restriction
!                 work(3) <= hnew/hold <= work(4)
!              default values: work(3)=0.333d0, work(4)=6.d0
!
!    work(5)   is the "beta" for stabilized step size control
!              (see section iv.2). positive values of beta ( <= 0.04 )
!              make the step size control more stable.
!              negative work(5) provoke beta=0.
!              default 0.0d0.
!
!    work(6)   maximal step size, default xend-x.
!
!    work(7)   initial step size, for work(7)=0.d0 an initial guess
!              is computed with help of the function hinit
!
!    iwork(1)  this is the maximal number of allowed steps.
!              the default value (for iwork(1)=0) is 100000.
!
!    iwork(2)  switch for the choice of the coefficients
!              if iwork(2).eq.1  method dop853 of dormand and prince
!              (section ii.6).
!              the default value (for iwork(2)=0) is iwork(2)=1.
!
!    iwork(3)  switch for printing error messages
!              if iwork(3).lt.0 no messages are being printed
!              if iwork(3).gt.0 messages are printed with
!              write (iwork(3),*) ...
!              default value (for iwork(3)=0) is iwork(3)=6
!
!    iwork(4)  test for stiffness is activated after step number
!              j*iwork(4) (j integer), provided iwork(4).gt.0.
!              for negative iwork(4) the stiffness test is
!              never activated; default value is iwork(4)=1000
!
!    iwork(5)  = nrdens = number of components, for which dense output
!              is required; default value is iwork(5)=0;
!              for   0 < nrdens < n   the components (for which dense
!              output is required) have to be specified in
!              iwork(21),...,iwork(nrdens+20);
!              for  nrdens=n  this is done by the code.
!
!----------------------------------------------------------------------
!
!     output parameters
!     -----------------
!     x           x-value for which the solution has been computed
!                 (after successful return x=xend).
!
!     y(n)        numerical solution at x
!
!     h           predicted step size of the last accepted step
!
!     idid        reports on successfulness upon return:
!                   idid= 1  computation successful,
!                   idid= 2  comput. successful (interrupted by solout)
!                   idid=-1  input is not consistent,
!                   idid=-2  larger nmax is needed,
!                   idid=-3  step size becomes too small.
!                   idid=-4  problem is probably stiff (interrupted).
!
!   iwork(17)  nfcn    number of function evaluations
!   iwork(18)  nstep   number of computed steps
!   iwork(19)  naccpt  number of accepted steps
!   iwork(20)  nrejct  number of rejected steps (due to error test),
!                      (step rejections in the first step are not counted)
!-----------------------------------------------------------------------

      subroutine dop853(n,fcn,x,y,xend,rtol,atol,itol,solout,iout,work, &
                      & lwork,iwork,liwork,idid)

      implicit none

      double precision atol , beta , fac1 , fac2 , h , hmax ,     &
                     & rtol , safe , work , x ,&
                     & xend , y
      integer i , icomp , idid , ieco , iek1 , iek10 , iek2 , iek3 ,    &
            & iek4 , iek5 , iek6 , iek7 , iek8 , iek9 , iey1 , iout ,   &
            & iprint , istore , itol
      integer iwork , liwork , lwork , meth , n , naccpt , nfcn , nmax ,&
            & nrdens , nrejct , nstep , nstiff
      dimension y(n) , atol(*) , rtol(*) , work(lwork) , iwork(liwork)
      logical arret
      procedure(func) :: fcn  !! subroutine computing the value of f(x,y)
      external solout

! *** *** *** *** *** *** ***
!        setting the parameters
! *** *** *** *** *** *** ***
      nfcn = 0
      nstep = 0
      naccpt = 0
      nrejct = 0
      arret = .false.
! -------- iprint for monitoring the printing
      if ( iwork(3)==0 ) then
         iprint = 6
      else
         iprint = iwork(3)
      endif
! -------- nmax , the maximal number of steps -----
      if ( iwork(1)==0 ) then
         nmax = 100000
      else
         nmax = iwork(1)
         if ( nmax<=0 ) then
            if ( iprint>0 ) write (iprint,*) ' wrong input iwork(1)=' , &
                                 & iwork(1)
            arret = .true.
         endif
      endif
! -------- meth   coefficients of the method
      if ( iwork(2)==0 ) then
         meth = 1
      else
         meth = iwork(2)
         if ( meth<=0 .or. meth>=4 ) then
            if ( iprint>0 ) write (iprint,*)                            &
                                 & ' curious input iwork(2)=' , iwork(2)
            arret = .true.
         endif
      endif
! -------- nstiff   parameter for stiffness detection
      nstiff = iwork(4)
      if ( nstiff==0 ) nstiff = 1000
      if ( nstiff<0 ) nstiff = nmax + 10
! -------- nrdens   number of dense output components
      nrdens = iwork(5)
      if ( nrdens<0 .or. nrdens>n ) then
         if ( iprint>0 ) write (iprint,*) ' curious input iwork(5)=' ,  &
                              & iwork(5)
         arret = .true.
      else
         if ( nrdens>0 .and. iout<2 ) then
            if ( iprint>0 ) write (iprint,*)                            &
                     &' warning: put iout=2 or iout=3 for dense output '
         endif
         if ( nrdens==n ) then
            do i = 1 , nrdens
               iwork(i+20) = i
            enddo
         endif
      endif
! -------  safety factor -------------
      if ( work(2)==0.d0 ) then
         safe = 0.9d0
      else
         safe = work(2)
         if ( safe>=1.d0 .or. safe<=1.d-4 ) then
            if ( iprint>0 ) write (iprint,*)                            &
                            &' curious input for safety factor work(2)='&
                           & , work(2)
            arret = .true.
         endif
      endif
! -------  fac1,fac2     parameters for step size selection
      if ( work(3)==0.d0 ) then
         fac1 = 0.333d0
      else
         fac1 = work(3)
      endif
      if ( work(4)==0.d0 ) then
         fac2 = 6.d0
      else
         fac2 = work(4)
      endif
! --------- beta for step control stabilization -----------
      if ( work(5)==0.d0 ) then
         beta = 0.0d0
      elseif ( work(5)<0.d0 ) then
         beta = 0.d0
      else
         beta = work(5)
         if ( beta>0.2d0 ) then
            if ( iprint>0 ) write (iprint,*)                            &
                                  & ' curious input for beta: work(5)=' &
                                 & , work(5)
            arret = .true.
         endif
      endif
! -------- maximal step size
      if ( work(6)==0.d0 ) then
         hmax = xend - x
      else
         hmax = work(6)
      endif
! -------- initial step size
      h = work(7)
! ------- prepare the entry-points for the arrays in work -----
      iek1 = 21
      iek2 = iek1 + n
      iek3 = iek2 + n
      iek4 = iek3 + n
      iek5 = iek4 + n
      iek6 = iek5 + n
      iek7 = iek6 + n
      iek8 = iek7 + n
      iek9 = iek8 + n
      iek10 = iek9 + n
      iey1 = iek10 + n
      ieco = iey1 + n
! ------ total storage requirement -----------
      istore = ieco + 8*nrdens - 1
      if ( istore>lwork ) then
         if ( iprint>0 ) write (iprint,*)                               &
                          &' insufficient storage for work, min. lwork='&
                         & , istore
         arret = .true.
      endif
      icomp = 21
      istore = icomp + nrdens - 1
      if ( istore>liwork ) then
         if ( iprint>0 ) write (iprint,*)                               &
                        &' insufficient storage for iwork, min. liwork='&
                       & , istore
         arret = .true.
      endif
! -------- when a fail has occured, we return with idid=-1
      if ( arret ) then
         idid = -1
         return
      endif
! -------- call to core integrator ------------
      call dp86co(n,fcn,x,y,xend,hmax,h,rtol,atol,itol,iprint,solout,   &
                & iout,idid,nmax,meth,nstiff,safe,beta,fac1,fac2,&
                & work(iek1),work(iek2),work(iek3),work(iek4),work(iek5)&
                & ,work(iek6),work(iek7),work(iek8),work(iek9),         &
                & work(iek10),work(iey1),work(ieco),iwork(icomp),nrdens,&
                & nfcn,nstep,naccpt,nrejct)
      work(7) = h
      iwork(17) = nfcn
      iwork(18) = nstep
      iwork(19) = naccpt
      iwork(20) = nrejct

      end subroutine dop853

! ----------------------------------------------------------
!     core integrator for dop853
!     parameters same as in dop853 with workspace added

      subroutine dp86co(n,fcn,x,y,xend,hmax,h,rtol,atol,itol,iprint,    &
                      & solout,iout,idid,nmax,meth,nstiff,safe,  &
                      & beta,fac1,fac2,k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,  &
                      & y1,cont,icomp,nrd,nfcn,nstep,naccpt,  &
                      & nrejct)

      implicit none

      double precision :: atoli,beta,bspl,cont,deno,err
      double precision err2 , erri , expo1 , fac , fac1 , fac11 , fac2 ,&
                     & facc1 , facc2 , facold , h , hlamb ,     &
                     & hmax , hnew , hout , posneg , rtoli , safe
      double precision sk , stden , stnum , x , xend , xold ,  &
                     & xout , xph , ydiff
      integer i , iasti , icomp , idid , iord , iout , iprint ,  &
            & irtrn , itol , j , meth , n , naccpt , nfcn , nmax ,      &
            & nonsti , nrd , nrejct , nstep
      integer nstiff

    double precision,parameter :: c2 = 0.526001519587677318785587544488d-01
    double precision,parameter :: c3 = 0.789002279381515978178381316732d-01
    double precision,parameter :: c4 = 0.118350341907227396726757197510d+00
    double precision,parameter :: c5 = 0.281649658092772603273242802490d+00
    double precision,parameter :: c6 = 0.333333333333333333333333333333d+00
    double precision,parameter :: c7 = 0.25d+00
    double precision,parameter :: c8 = 0.307692307692307692307692307692d+00
    double precision,parameter :: c9 = 0.651282051282051282051282051282d+00
    double precision,parameter :: c10 = 0.6d+00
    double precision,parameter :: c11 = 0.857142857142857142857142857142d+00
    double precision,parameter :: c14 = 0.1d+00
    double precision,parameter :: c15 = 0.2d+00
    double precision,parameter :: c16 = 0.777777777777777777777777777778d+00
    double precision,parameter :: b1 = 5.42937341165687622380535766363d-2
    double precision,parameter :: b6 = 4.45031289275240888144113950566d0
    double precision,parameter :: b7 = 1.89151789931450038304281599044d0
    double precision,parameter :: b8 = -5.8012039600105847814672114227d0
    double precision,parameter :: b9 = 3.1116436695781989440891606237d-1
    double precision,parameter :: b10 = -1.52160949662516078556178806805d-1
    double precision,parameter :: b11 = 2.01365400804030348374776537501d-1
    double precision,parameter :: b12 = 4.47106157277725905176885569043d-2
    double precision,parameter :: bhh1 = 0.244094488188976377952755905512d+00
    double precision,parameter :: bhh2 = 0.733846688281611857341361741547d+00
    double precision,parameter :: bhh3 = 0.220588235294117647058823529412d-01
    double precision,parameter :: er1 = 0.1312004499419488073250102996d-01
    double precision,parameter :: er6 = -0.1225156446376204440720569753d+01
    double precision,parameter :: er7 = -0.4957589496572501915214079952d+00
    double precision,parameter :: er8 = 0.1664377182454986536961530415d+01
    double precision,parameter :: er9 = -0.3503288487499736816886487290d+00
    double precision,parameter :: er10 = 0.3341791187130174790297318841d+00
    double precision,parameter :: er11 = 0.8192320648511571246570742613d-01
    double precision,parameter :: er12 = -0.2235530786388629525884427845d-01
    double precision,parameter :: a21 = 5.26001519587677318785587544488d-2
    double precision,parameter :: a31 = 1.97250569845378994544595329183d-2
    double precision,parameter :: a32 = 5.91751709536136983633785987549d-2
    double precision,parameter :: a41 = 2.95875854768068491816892993775d-2
    double precision,parameter :: a43 = 8.87627564304205475450678981324d-2
    double precision,parameter :: a51 = 2.41365134159266685502369798665d-1
    double precision,parameter :: a53 = -8.84549479328286085344864962717d-1
    double precision,parameter :: a54 = 9.24834003261792003115737966543d-1
    double precision,parameter :: a61 = 3.7037037037037037037037037037d-2
    double precision,parameter :: a64 = 1.70828608729473871279604482173d-1
    double precision,parameter :: a65 = 1.25467687566822425016691814123d-1
    double precision,parameter :: a71 = 3.7109375d-2
    double precision,parameter :: a74 = 1.70252211019544039314978060272d-1
    double precision,parameter :: a75 = 6.02165389804559606850219397283d-2
    double precision,parameter :: a76 = -1.7578125d-2
    double precision,parameter :: a81 = 3.70920001185047927108779319836d-2
    double precision,parameter :: a84 = 1.70383925712239993810214054705d-1
    double precision,parameter :: a85 = 1.07262030446373284651809199168d-1
    double precision,parameter :: a86 = -1.53194377486244017527936158236d-2
    double precision,parameter :: a87 = 8.27378916381402288758473766002d-3
    double precision,parameter :: a91 = 6.24110958716075717114429577812d-1
    double precision,parameter :: a94 = -3.36089262944694129406857109825d0
    double precision,parameter :: a95 = -8.68219346841726006818189891453d-1
    double precision,parameter :: a96 = 2.75920996994467083049415600797d1
    double precision,parameter :: a97 = 2.01540675504778934086186788979d1
    double precision,parameter :: a98 = -4.34898841810699588477366255144d1
    double precision,parameter :: a101 = 4.77662536438264365890433908527d-1
    double precision,parameter :: a104 = -2.48811461997166764192642586468d0
    double precision,parameter :: a105 = -5.90290826836842996371446475743d-1
    double precision,parameter :: a106 = 2.12300514481811942347288949897d1
    double precision,parameter :: a107 = 1.52792336328824235832596922938d1
    double precision,parameter :: a108 = -3.32882109689848629194453265587d1
    double precision,parameter :: a109 = -2.03312017085086261358222928593d-2
    double precision,parameter :: a111 = -9.3714243008598732571704021658d-1
    double precision,parameter :: a114 = 5.18637242884406370830023853209d0
    double precision,parameter :: a115 = 1.09143734899672957818500254654d0
    double precision,parameter :: a116 = -8.14978701074692612513997267357d0
    double precision,parameter :: a117 = -1.85200656599969598641566180701d1
    double precision,parameter :: a118 = 2.27394870993505042818970056734d1
    double precision,parameter :: a119 = 2.49360555267965238987089396762d0
    double precision,parameter :: a1110 = -3.0467644718982195003823669022d0
    double precision,parameter :: a121 = 2.27331014751653820792359768449d0
    double precision,parameter :: a124 = -1.05344954667372501984066689879d1
    double precision,parameter :: a125 = -2.00087205822486249909675718444d0
    double precision,parameter :: a126 = -1.79589318631187989172765950534d1
    double precision,parameter :: a127 = 2.79488845294199600508499808837d1
    double precision,parameter :: a128 = -2.85899827713502369474065508674d0
    double precision,parameter :: a129 = -8.87285693353062954433549289258d0
    double precision,parameter :: a1210 = 1.23605671757943030647266201528d1
    double precision,parameter :: a1211 = 6.43392746015763530355970484046d-1
    double precision,parameter :: a141 = 5.61675022830479523392909219681d-2
    double precision,parameter :: a147 = 2.53500210216624811088794765333d-1
    double precision,parameter :: a148 = -2.46239037470802489917441475441d-1
    double precision,parameter :: a149 = -1.24191423263816360469010140626d-1
    double precision,parameter :: a1410 = 1.5329179827876569731206322685d-1
    double precision,parameter :: a1411 = 8.20105229563468988491666602057d-3
    double precision,parameter :: a1412 = 7.56789766054569976138603589584d-3
    double precision,parameter :: a1413 = -8.298d-3
    double precision,parameter :: a151 = 3.18346481635021405060768473261d-2
    double precision,parameter :: a156 = 2.83009096723667755288322961402d-2
    double precision,parameter :: a157 = 5.35419883074385676223797384372d-2
    double precision,parameter :: a158 = -5.49237485713909884646569340306d-2
    double precision,parameter :: a1511 = -1.08347328697249322858509316994d-4
    double precision,parameter :: a1512 = 3.82571090835658412954920192323d-4
    double precision,parameter :: a1513 = -3.40465008687404560802977114492d-4
    double precision,parameter :: a1514 = 1.41312443674632500278074618366d-1
    double precision,parameter :: a161 = -4.28896301583791923408573538692d-1
    double precision,parameter :: a166 = -4.69762141536116384314449447206d0
    double precision,parameter :: a167 = 7.68342119606259904184240953878d0
    double precision,parameter :: a168 = 4.06898981839711007970213554331d0
    double precision,parameter :: a169 = 3.56727187455281109270669543021d-1
    double precision,parameter :: a1613 = -1.39902416515901462129418009734d-3
    double precision,parameter :: a1614 = 2.9475147891527723389556272149d0
    double precision,parameter :: a1615 = -9.15095847217987001081870187138d0
    double precision,parameter :: d41 = -0.84289382761090128651353491142d+01
    double precision,parameter :: d46 = 0.56671495351937776962531783590d+00
    double precision,parameter :: d47 = -0.30689499459498916912797304727d+01
    double precision,parameter :: d48 = 0.23846676565120698287728149680d+01
    double precision,parameter :: d49 = 0.21170345824450282767155149946d+01
    double precision,parameter :: d410 = -0.87139158377797299206789907490d+00
    double precision,parameter :: d411 = 0.22404374302607882758541771650d+01
    double precision,parameter :: d412 = 0.63157877876946881815570249290d+00
    double precision,parameter :: d413 = -0.88990336451333310820698117400d-01
    double precision,parameter :: d414 = 0.18148505520854727256656404962d+02
    double precision,parameter :: d415 = -0.91946323924783554000451984436d+01
    double precision,parameter :: d416 = -0.44360363875948939664310572000d+01
    double precision,parameter :: d51 = 0.10427508642579134603413151009d+02
    double precision,parameter :: d56 = 0.24228349177525818288430175319d+03
    double precision,parameter :: d57 = 0.16520045171727028198505394887d+03
    double precision,parameter :: d58 = -0.37454675472269020279518312152d+03
    double precision,parameter :: d59 = -0.22113666853125306036270938578d+02
    double precision,parameter :: d510 = 0.77334326684722638389603898808d+01
    double precision,parameter :: d511 = -0.30674084731089398182061213626d+02
    double precision,parameter :: d512 = -0.93321305264302278729567221706d+01
    double precision,parameter :: d513 = 0.15697238121770843886131091075d+02
    double precision,parameter :: d514 = -0.31139403219565177677282850411d+02
    double precision,parameter :: d515 = -0.93529243588444783865713862664d+01
    double precision,parameter :: d516 = 0.35816841486394083752465898540d+02
    double precision,parameter :: d61 = 0.19985053242002433820987653617d+02
    double precision,parameter :: d66 = -0.38703730874935176555105901742d+03
    double precision,parameter :: d67 = -0.18917813819516756882830838328d+03
    double precision,parameter :: d68 = 0.52780815920542364900561016686d+03
    double precision,parameter :: d69 = -0.11573902539959630126141871134d+02
    double precision,parameter :: d610 = 0.68812326946963000169666922661d+01
    double precision,parameter :: d611 = -0.10006050966910838403183860980d+01
    double precision,parameter :: d612 = 0.77771377980534432092869265740d+00
    double precision,parameter :: d613 = -0.27782057523535084065932004339d+01
    double precision,parameter :: d614 = -0.60196695231264120758267380846d+02
    double precision,parameter :: d615 = 0.84320405506677161018159903784d+02
    double precision,parameter :: d616 = 0.11992291136182789328035130030d+02
    double precision,parameter :: d71 = -0.25693933462703749003312586129d+02
    double precision,parameter :: d76 = -0.15418974869023643374053993627d+03
    double precision,parameter :: d77 = -0.23152937917604549567536039109d+03
    double precision,parameter :: d78 = 0.35763911791061412378285349910d+03
    double precision,parameter :: d79 = 0.93405324183624310003907691704d+02
    double precision,parameter :: d710 = -0.37458323136451633156875139351d+02
    double precision,parameter :: d711 = 0.10409964950896230045147246184d+03
    double precision,parameter :: d712 = 0.29840293426660503123344363579d+02
    double precision,parameter :: d713 = -0.43533456590011143754432175058d+02
    double precision,parameter :: d714 = 0.96324553959188282948394950600d+02
    double precision,parameter :: d715 = -0.39177261675615439165231486172d+02
    double precision,parameter :: d716 = -0.14972683625798562581422125276d+03

      double precision y(n) , y1(n) , k1(n) , k2(n) , k3(n) , k4(n) , k5(n) , k6(n)
      double precision k7(n) , k8(n) , k9(n) , k10(n) , atol(*) , rtol(*)
      dimension cont(8*nrd) , icomp(nrd)
      logical reject , last , event
      procedure(func) :: fcn
      common /condo8/ xold , hout

! *** *** *** *** *** *** ***
!  initialisations
! *** *** *** *** *** *** ***
      facold = 1.d-4
      expo1 = 1.d0/8.d0 - beta*0.2d0
      facc1 = 1.d0/fac1
      facc2 = 1.d0/fac2
      posneg = sign(1.d0,xend-x)
! --- initial preparations
      atoli = atol(1)
      rtoli = rtol(1)
      last = .false.
      hlamb = 0.d0
      iasti = 0
      call fcn(n,x,y,k1)
      hmax = abs(hmax)
      iord = 8
      if ( h==0.d0 ) h = hinit(n,fcn,x,y,xend,posneg,k1,k2,k3,iord,hmax,&
                       & atol,rtol,itol)
      nfcn = nfcn + 2
      reject = .false.
      xold = x
      if ( iout/=0 ) then
         irtrn = 1
         hout = 1.d0
         call solout(naccpt+1,xold,x,y,n,cont,icomp,nrd,irtrn,&
                   & xout)
         if ( irtrn<0 ) goto 200
      endif
! --- basic integration step
 100  if ( nstep>nmax ) then
         if ( iprint>0 ) write (iprint,99001) x
         if ( iprint>0 ) write (iprint,*) ' more than nmax =' , nmax , 'steps are needed'
         idid = -2
         return
      elseif ( 0.1d0*abs(h)<=abs(x)*uround ) then
         if ( iprint>0 ) write (iprint,99001) x
         if ( iprint>0 ) write (iprint,*) ' step size too small, h=' , h
         idid = -3
         return
      else
         if ( (x+1.01d0*h-xend)*posneg>0.d0 ) then
            h = xend - x
            last = .true.
         endif
         nstep = nstep + 1
! --- the twelve stages
         if ( irtrn>=2 ) call fcn(n,x,y,k1)
         y1 = y + h*a21*k1
         call fcn(n,x+c2*h,y1,k2)
         y1 = y + h*(a31*k1+a32*k2)
         call fcn(n,x+c3*h,y1,k3)
         y1 = y + h*(a41*k1+a43*k3)
         call fcn(n,x+c4*h,y1,k4)
         y1 = y + h*(a51*k1+a53*k3+a54*k4)
         call fcn(n,x+c5*h,y1,k5)
         y1 = y + h*(a61*k1+a64*k4+a65*k5)
         call fcn(n,x+c6*h,y1,k6)
         y1 = y + h*(a71*k1+a74*k4+a75*k5+a76*k6)
         call fcn(n,x+c7*h,y1,k7)
         y1 = y + h*(a81*k1+a84*k4+a85*k5+a86*k6+a87*k7)
         call fcn(n,x+c8*h,y1,k8)
         y1 = y + h*(a91*k1+a94*k4+a95*k5+a96*k6+a97*k7+a98*k8)
         call fcn(n,x+c9*h,y1,k9)
         y1 = y + h*(a101*k1+a104*k4+a105*k5+a106*k6+a107*k7+a108*k8+a109*k9)
         call fcn(n,x+c10*h,y1,k10)
         y1 = y + h*(a111*k1+a114*k4+a115*k5+a116*k6+a117*k7+a118*k8+a119*k9+a1110*k10)
         call fcn(n,x+c11*h,y1,k2)
         xph = x + h
         y1 = y + h*(a121*k1+a124*k4+a125*k5+a126*k6+a127*k7+a128*k8+a129*k9+a1210*k10+a1211*k2)
         call fcn(n,xph,y1,k3)
         nfcn = nfcn + 11
         k4 = b1*k1 + b6*k6 + b7*k7 + b8*k8 + b9*k9 + b10*k10 + b11*k2 + b12*k3
         k5 = y + h*k4
! --- error estimation
         err = 0.d0
         err2 = 0.d0
         if ( itol==0 ) then
            do i = 1 , n
               sk = atoli + rtoli*max(abs(y(i)),abs(k5(i)))
               erri = k4(i) - bhh1*k1(i) - bhh2*k9(i) - bhh3*k3(i)
               err2 = err2 + (erri/sk)**2
               erri = er1*k1(i) + er6*k6(i) + er7*k7(i) + er8*k8(i)     &
                    & + er9*k9(i) + er10*k10(i) + er11*k2(i)            &
                    & + er12*k3(i)
               err = err + (erri/sk)**2
            enddo
         else
            do i = 1 , n
               sk = atol(i) + rtol(i)*max(abs(y(i)),abs(k5(i)))
               erri = k4(i) - bhh1*k1(i) - bhh2*k9(i) - bhh3*k3(i)
               err2 = err2 + (erri/sk)**2
               erri = er1*k1(i) + er6*k6(i) + er7*k7(i) + er8*k8(i)     &
                    & + er9*k9(i) + er10*k10(i) + er11*k2(i)            &
                    & + er12*k3(i)
               err = err + (erri/sk)**2
            enddo
         endif
         deno = err + 0.01d0*err2
         if ( deno<=0.d0 ) deno = 1.d0
         err = abs(h)*err*sqrt(1.d0/(n*deno))
! --- computation of hnew
         fac11 = err**expo1
! --- lund-stabilization
         fac = fac11/facold**beta
! --- we require  fac1 <= hnew/h <= fac2
         fac = max(facc2,min(facc1,fac/safe))
         hnew = h/fac
         if ( err<=1.d0 ) then
! --- step is accepted
            facold = max(err,1.0d-4)
            naccpt = naccpt + 1
            call fcn(n,xph,k5,k4)
            nfcn = nfcn + 1
! ------- stiffness detection
            if ( mod(naccpt,nstiff)==0 .or. iasti>0 ) then
               stnum = 0.d0
               stden = 0.d0
               do i = 1 , n
                  stnum = stnum + (k4(i)-k3(i))**2
                  stden = stden + (k5(i)-y1(i))**2
               enddo
               if ( stden>0.d0 ) hlamb = abs(h)*sqrt(stnum/stden)
               if ( hlamb>6.1d0 ) then
                  nonsti = 0
                  iasti = iasti + 1
                  if ( iasti==15 ) then
                     if ( iprint>0 ) write (iprint,*) ' the problem seems to become stiff at x = ', x
                     if ( iprint<=0 ) then
                        ! --- fail exit
                        idid = -4
                        return
                     endif
                  endif
               else
                  nonsti = nonsti + 1
                  if ( nonsti==6 ) iasti = 0
               endif
            endif
! ------- final preparation for dense output
            event = (iout==3) .and. (xout<=xph)
            if ( iout==2 .or. event ) then
! ----    save the first function evaluations
               do j = 1 , nrd
                  i = icomp(j)
                  cont(j) = y(i)
                  ydiff = k5(i) - y(i)
                  cont(j+nrd) = ydiff
                  bspl = h*k1(i) - ydiff
                  cont(j+nrd*2) = bspl
                  cont(j+nrd*3) = ydiff - h*k4(i) - bspl
                  cont(j+nrd*4) = d41*k1(i) + d46*k6(i) + d47*k7(i) + d48*k8(i) + d49*k9(i) + d410*k10(i) + d411*k2(i) + d412*k3(i)
                  cont(j+nrd*5) = d51*k1(i) + d56*k6(i) + d57*k7(i) + d58*k8(i) + d59*k9(i) + d510*k10(i) + d511*k2(i) + d512*k3(i)
                  cont(j+nrd*6) = d61*k1(i) + d66*k6(i) + d67*k7(i) + d68*k8(i) + d69*k9(i) + d610*k10(i) + d611*k2(i) + d612*k3(i)
                  cont(j+nrd*7) = d71*k1(i) + d76*k6(i) + d77*k7(i) + d78*k8(i) + d79*k9(i) + d710*k10(i) + d711*k2(i) + d712*k3(i)
               enddo
! ---     the next three function evaluations
               y1 = y + h*(a141*k1+a147*k7+a148*k8+a149*k9+a1410*k10+a1411*k2+a1412*k3+a1413*k4)
               call fcn(n,x+c14*h,y1,k10)
               y1 = y + h*(a151*k1+a156*k6+a157*k7+a158*k8+a1511*k2+a1512*k3+a1513*k4+a1514*k10)
               call fcn(n,x+c15*h,y1,k2)
               y1 = y + h*(a161*k1+a166*k6+a167*k7+a168*k8+a169*k9+a1613*k4+a1614*k10+a1615*k2)
               call fcn(n,x+c16*h,y1,k3)
               nfcn = nfcn + 3
! ---     final preparation
               do j = 1 , nrd
                  i = icomp(j)
                  cont(j+nrd*4) = h*(cont(j+nrd*4)+d413*k4(i)+d414*k10(i)+d415*k2(i)+d416*k3(i))
                  cont(j+nrd*5) = h*(cont(j+nrd*5)+d513*k4(i)+d514*k10(i)+d515*k2(i)+d516*k3(i))
                  cont(j+nrd*6) = h*(cont(j+nrd*6)+d613*k4(i)+d614*k10(i)+d615*k2(i)+d616*k3(i))
                  cont(j+nrd*7) = h*(cont(j+nrd*7)+d713*k4(i)+d714*k10(i)+d715*k2(i)+d716*k3(i))
               enddo
               hout = h
            endif
            k1 = k4
            y = k5
            xold = x
            x = xph
            if ( iout==1 .or. iout==2 .or. event ) then
               call solout(naccpt+1,xold,x,y,n,cont,icomp,nrd,irtrn,xout)
               if ( irtrn<0 ) goto 200
            endif
! ------- normal exit
            if ( last ) then
               h = hnew
               idid = 1
               return
            endif
            if ( abs(hnew)>hmax ) hnew = posneg*hmax
            if ( reject ) hnew = posneg*min(abs(hnew),abs(h))
            reject = .false.
         else
! --- step is rejected
            hnew = h/min(facc1,fac11/safe)
            reject = .true.
            if ( naccpt>=1 ) nrejct = nrejct + 1
            last = .false.
         endif
         h = hnew
         goto 100
      endif
 200  if ( iprint>0 ) write (iprint,99001) x
      idid = 2
99001 format (' exit of dop853 at x=',e18.4)

      end subroutine dp86co

! ----------------------------------------------------------
! ----  computation of an initial step size guess
! ----------------------------------------------------------
    function hinit(n,fcn,x,y,xend,posneg,f0,f1,y1,iord,hmax,atol,rtol,itol)

      implicit none

      double precision atol , atoli , der12 , der2 , dnf , dny , f0 ,   &
                     & f1 , h , h1 , hinit , hmax , posneg ,     &
                     & rtol , rtoli , sk , x , xend , y
      double precision y1
      integer i , iord , itol , n
      dimension y(n) , y1(n) , f0(n) , f1(n) , atol(*) , rtol(*)
      procedure(func) :: fcn

! ---- compute a first guess for explicit euler as
! ----   h = 0.01 * norm (y0) / norm (f0)
! ---- the increment for explicit euler is small
! ---- compared to the solution
      dnf = 0.0d0
      dny = 0.0d0
      atoli = atol(1)
      rtoli = rtol(1)
      if ( itol==0 ) then
         do i = 1 , n
            sk = atoli + rtoli*abs(y(i))
            dnf = dnf + (f0(i)/sk)**2
            dny = dny + (y(i)/sk)**2
         enddo
      else
         do i = 1 , n
            sk = atol(i) + rtol(i)*abs(y(i))
            dnf = dnf + (f0(i)/sk)**2
            dny = dny + (y(i)/sk)**2
         enddo
      endif
      if ( dnf<=1.d-10 .or. dny<=1.d-10 ) then
         h = 1.0d-6
      else
         h = sqrt(dny/dnf)*0.01d0
      endif
      h = min(h,hmax)
      h = sign(h,posneg)
! ---- perform an explicit euler step
      do i = 1 , n
         y1(i) = y(i) + h*f0(i)
      enddo
      call fcn(n,x+h,y1,f1)
! ---- estimate the second derivative of the solution
      der2 = 0.0d0
      if ( itol==0 ) then
         do i = 1 , n
            sk = atoli + rtoli*abs(y(i))
            der2 = der2 + ((f1(i)-f0(i))/sk)**2
         enddo
      else
         do i = 1 , n
            sk = atol(i) + rtol(i)*abs(y(i))
            der2 = der2 + ((f1(i)-f0(i))/sk)**2
         enddo
      endif
      der2 = sqrt(der2)/h
! ---- step size is computed such that
! ----  h**iord * max ( norm (f0), norm (der2)) = 0.01
      der12 = max(abs(der2),sqrt(dnf))
      if ( der12<=1.d-15 ) then
         h1 = max(1.0d-6,abs(h)*1.0d-3)
      else
         h1 = (0.01d0/der12)**(1.d0/iord)
      endif
      h = min(100*abs(h),h1,hmax)
      hinit = sign(h,posneg)

      end function hinit

! ----------------------------------------------------------
!     this function can be used for coninuous output in connection
!     with the output-subroutine for dop853. it provides an
!     approximation to the ii-th component of the solution at x.
! ----------------------------------------------------------
    function contd8(ii,x,con,icomp,nd)

      implicit none

      double precision con , conpar , contd8 , h , s , s1 , x , xold
      integer i , icomp , ii , j , nd
      dimension con(8*nd) , icomp(nd)
      common /condo8/ xold , h

! ----- compute place of ii-th component
      i = 0
      do j = 1 , nd
         if ( icomp(j)==ii ) i = j
      enddo
      if ( i==0 ) then
         write (6,*) ' no dense output available for comp.' , ii
      else
          s = (x-xold)/h
          s1 = 1.d0 - s
          conpar = con(i+nd*4) + s*(con(i+nd*5)+s1*(con(i+nd*6)+s*con(i+nd*7)))
          contd8 = con(i) + s*(con(i+nd)+s1*(con(i+nd*2)+s*(con(i+nd*3)+s1*conpar)))
      end if
    end function contd8

end module dop853_module
