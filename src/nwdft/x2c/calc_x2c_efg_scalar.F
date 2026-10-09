      subroutine calc_x2c_efg_scalar (rtdb, needu2c, nexc,
     &   g_dscf, g_dens_at, g_u, g_u2c, g_c2u)

c     -----------------------------------------------------------------
c     purpose: calculate one-electron property EFG using the 
c     two-component exact-decoupling Hamiltonian

c     For the time being we MUST use an uncontracted basis !

c     Integrals are supposed to be initialized by the calling
c     routine.

c     Driver (c) 2012, 2022 J. Autschbach, SUNY, jochena@buffalo.edu

c     Approach: The transformation matrices g_u etc. come from
c     calc_x2c_1e if certain options are set.

c     instead of calculating the operator in the uncontracted basis,
c     then re-contracting, and then calculating the property, we keep
c     the operator in the uncontracted basis and simply transform the
c     density matrix to the same basis, assuming the result to be
c     equivalent.
c
c     -----------------------------------------------------------------

      implicit none

#include "mafdecls.fh"
#include "rtdb.fh"
#include "apiP.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "bas.fh"
#include "util.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "dftpara.fh"
#include "cdft.fh"
#include "geom.fh"
c#include "rel_consts.fh"
#include "cosmo.fh"
#include "x2c_scf.fh"
#include "cint1cache.fh"

c     subroutine arguments:

      integer rtdb          ! [in] runtime database GA handle
      integer nexc          ! [in] no. of XC terms
      integer g_dscf(2)     ! [in] density matrix from SCF
      integer g_dens_at(2)  ! [in] atom densities for model pot.
      integer g_u(2)        ! [in] decoupling matrices
      integer g_u2c         ! [in] basis transformation (optional)
      integer g_c2u         ! [in] basis transformation (optional)
      logical needu2c       ! [in] switch for u2c/c2u

      
c     local variables:

      integer me, type, i, j
      double precision sum
      character*2  symbol
      character*16 element, at_tag

c     GA handles used for temp storage in this routine:
      integer g_P, g_efgint, g_temp, g_tmp2, g_efg, g_w(6),
     &  g_densu, g_s, g_wtmp

c     dimensioning arrays
      integer alo(3), ahi(3), blo(3), bhi(3), clo(3), chi(3)

c     other MA handles
      integer l_xyzpt, k_xyzpt, l_zanpt, k_zanpt

c     bq variables (MV)
      logical dobq
      integer bq_ncent
      integer i_cbq
      integer i_qbq
      double precision elpotbq

c     other local variables:

      integer ubas, cbas
      integer nubas, ncbas

      double precision ga_trace_diag
      external ga_trace_diag
      double precision tol2e, accval
      logical master, analytic_pvp, status, oskel_sav
      integer dbg
      
      integer ndmem, ntemp
      integer iat, nat
      
      double precision xp, yp, zp, xn, yn, zn, zan
      double precision rr, rr5
      double precision efgel(6), efgnuc(6), etmp(6)
      double precision origin(3)

      double precision pifac, pi, deg, third
      parameter (third=1.0d0/3.0d0)

      character*14 pname, ctmp

c     =================================================================

c     ---------------
c     initializations
c     ---------------

c     options / settings

      dbg=1

      me=ga_nodeid()
      master = me.eq.0

      pname = "x2c_efg_scalar"

      if (dbg>0.and.master) write(luout,*) 'entering '//pname

      call ga_sync()

c     src/NWints/dk/int_1edk_ga.F uses similar transformations as used
c     here. There, oskel is deliberately set to .false., so I'm guessing
c     we should do the same here but I'm not sure yet

      oskel_sav = oskel
      oskel = .false.

c     constants:

      pi  = acos(-1.0d0)
      deg = 180.0d0/pi
      pifac = (4.0d0/3.0d0)*pi

      tol2e=10.d0**(-itol2e)

      call int_acc_get(accval)

      if (dbg>0.and.master) write(luout,*) 'accval = ',accval
c      accval = 1d-25
c     call int_acc_set(accval)

c     ---------------
c     initializations
c     ---------------

c     the ncbas, nubas gymnastics is usually needed because
c     the x2c operator is always calculated in an uncontracted
c     basis

      ncbas = 0
      nubas = 0

      if (dbg>0.and.master) write(luout,*) pname//': needu2c =',needu2c

      if (needu2c) then

        call ga_inquire_name(g_u2c,ctmp)
        if(dbg>0.and.master) write (luout,*) 'u2c ga# and name: ',
     &    g_u2c, trim(ctmp)
        call ga_inquire_name(g_c2u,ctmp)
        if(dbg>0.and.master) write (luout,*) 'c2u ga# and name: ',
     &    g_c2u, trim(ctmp)

        ubas = x2c_bas

        if (.not. bas_numbf(ubas,nubas))
     &    call errquit(pname//': nbas x2c error',nubas, BASIS_ERR)
        cbas = ao_bas_han

        if (.not. bas_numbf(cbas,ncbas))
     &    call errquit(pname//': nbas ao error',nubas, BASIS_ERR)

        if (master) write(luout,*)
     &    pname//': working with uncontracted auxiliary basis'

c       make sure we have the correct u2c and c2u arrays

        call ga_inquire_name(g_u2c, ctmp)
        if (ctmp(1:3) .ne. 'u2c') call errquit(
     &    pname//': g_u2c inconsistency', 1, UNKNOWN_ERR)

        call ga_inquire_name(g_c2u, ctmp)
        if (ctmp(1:3) .ne. 'c2u') call errquit(
     &    pname//': g_c2u inconsistency', 1, UNKNOWN_ERR)

        if (dbg>0.and.master) write(luout,*) 'ncbas,nubas =',ncbas,nubas

      else

        ubas = ao_bas_han

        if (.not. bas_numbf(ubas,nubas))
     &    call errquit(pname//': nbas ao error',nubas, BASIS_ERR)

        cbas = ubas
        ncbas = nubas

        if (master) write(luout,*)
     &    pname//': the basis seems to be fully uncontracted ...'

      end if

c     make sure the input GAs for the transformation have the correct
c     dimension

      call ga_inquire(g_u(1),type,ntemp,ntemp)

      if (nubas.ne.ntemp) then
        call errquit(pname//': inconsistent dim. for g_u',nubas,GA_ERR)
      end if

      if (master.and.dbg>0) write (luout,*) pname//': nubas =',nubas

c     now that we have everything set for the contracted/uncontracted
c     scheme, let's calculate the density matrix in the uncontracted
c     basis if we do in fact have a contracted main basis, do some
c     tests, and perform the remaining calculation in the uncontracted
c     basis. The main complication is that the operator in the small
c     component basis is calculated numerically, and that gets done in
c     the contracted basis. Therefore, that step requires a similar step
c     as in calc_x2c_1e_scalar, as done in the part that deals with the
c     nummerical pVp integrals and the model potential, to deal with
c     this contribution.

      if (needu2c) then

        if (dbg>0.and.master) then
          write(luout,*) 'c2u*u2c test begin'
          write(luout,*) g_u2c, g_c2u
        end if

c       first, test whether g_c2u and g_u2c are present and are
c       pseudo-inverses (they should be, but it won't hurt checking
c       that). We're checking both diagonal elements and the
c       off-diagonal. The diagonal should equal the number of
c       contracted basis functions

        if (.not. ga_create(MT_DBL,ncbas,ncbas,'temp',1,1,g_temp))
     &    call errquit(pname//': error creating temp',56,GA_ERR)

        call ga_dgemm('n','n',ncbas,ncbas,nubas,1.0d0,g_c2u,g_u2c,
     &    0.0d0,g_temp)

        sum = ga_trace_diag(g_temp)
        if (master) write (luout,'(1x,a,1x,e15.7)')
     &    pname//': c2u*u2c: ',sum

        if (abs(sum-real(ncbas)).gt.1.0d-10) then
          if (master) write (luout,*) pname//
     &   ': WARNING: left-inverse in u/c trafo seems inaccurate !'
        end if

        call ga_add_todiag(g_temp, -1.0d0)
        call ga_maxelt(g_temp, sum)

        if (master .and. abs(sum).gt.1d-6) then
          write(luout,'(/1x,a,1x,f18.8)')
     &      pname//': largest element in c2u*u2c - 1: ',
     &      sum
          write(luout,'(1x,a/1x,a/)')
     &      'Significant deviations from zero may indicate',
     &      'elevated numerical noise in the basis de/re-contraction'
        end if

        if (.not. ga_destroy(g_temp))
     &    call errquit(pname//': error destroy temp',56, GA_ERR)

        if (dbg>0.and.master) write(luout,*) 'c2u*u2c test end'

c       second, calculate the density matrix g_densu in the uncontracted
c       basis. See calc_x2c_1e for further comments

        if (.not. ga_create(MT_DBL,ncbas,nubas,'temp',1,1,g_temp))
     &    call errquit(pname//': error creating temp',57,GA_ERR)

        if (.not. ga_create(MT_DBL,nubas,nubas,'densu',1,1,g_densu))
     &    call errquit(pname//': error creating densu',58,GA_ERR)

        if (master.and.dbg>0) write(luout,*) pname//': dens c2u start'

        call ga_zero(g_densu)
        call ga_dgemm('n','t',ncbas,nubas,ncbas, 1.0d0,
     &    g_dscf(1),g_u2c,0.0d0,g_temp)
        if (master.and.dbg>0) write(luout,*) pname//': dens c2u IM'
        call ga_dgemm('n','n',nubas,nubas,ncbas ,1.0d0,
     &    g_u2c, g_temp, 0.0d0,g_densu)
c        call ga_symmetrize(g_densu) ! not sure we need this

        if (.not. ga_destroy(g_temp))
     &    call errquit(pname//': error destroy temp',57, GA_ERR)

        if (master.and.dbg>0) write(luout,*) pname//': dens c2u done'

      else

c       if we aren't using a contracted basis elsewhere, we can
c       simply make the g_densu GA handle the same as that for
c       g_dscf(1)

        g_densu = g_dscf(1) ! sort of like the old F77 EQUIVALENCE ...

      end if                    ! needu2c

c     irrespective of whether the main basis is contracted or not, we
c     can now calculate N = tr[P S] with the array g_densu

      if (.not. ga_create(MT_DBL,nubas,nubas,'g_s',1,1,g_s)) call
     &  errquit(pname//': ga_create failed g_s',11,GA_ERR)
      call ga_zero(g_s)

      if (ocache) call int_1e_uncache_ga()
      call int_1e_ga(ubas, ubas , g_s, 'overlap', oskel)

!      call ga_symmetrize(g_s) ! not sure we need this

      call ga_sync() ! also not sure we need this ...

      if (.not. ga_create(MT_DBL,nubas,nubas,'g_temp',1,1,g_temp)) call
     &  errquit(pname//': ga_create g_temp failed',59,GA_ERR)

      call ga_zero(g_temp)
      call ga_dgemm('n','n',nubas,nubas,nubas, 1.0d0,
     &  g_densu, g_s, 0.0d0, g_temp)

      sum = ga_trace_diag(g_temp)

      if (master) write (luout,'(1x,a,1x,e15.7)')
     &  pname//': tr(P * S) : ',sum

      if (.not. ga_destroy(g_temp)) call errquit (pname//
     &  ': temp ga corrupt?',59, GA_ERR)
      if (.not. ga_destroy(g_s)) call errquit (pname//
     &  ': overlap ga corrupt?',11, GA_ERR)


c      ------------- Calculate EFG -------------

c     initialize integrals that were not already initialized

c      call hnd_giao_init(ubas,1)

c     ----- define points for calculation -----
c           1. nuclei: for now, we take ALL nuclei
c
      status=geom_ncent(geom,nat)
c
      if (.not. ma_push_get(mt_dbl,3*nat,'xyz pnt',l_xyzpt,k_xyzpt))
     &    call errquit(pname//': ma failed',911,MA_ERR)
      if (.not. ma_push_get(mt_dbl,nat,'zan pnt',l_zanpt,k_zanpt))
     &    call errquit(pname//': ma failed',911,MA_ERR)
c
      do  iat=1,nat
        status=geom_cent_get(geom,iat,at_tag,dbl_mb(k_xyzpt+3*(iat-1)),
     &     dbl_mb(k_zanpt+iat-1))
       end do

      alo(1) = nubas
      alo(2) = -1
      alo(3) = -1
      ahi(1) = nubas
      ahi(2) = nubas
      ahi(3) = 6*nat

      if (.not.nga_create(MT_DBL,3,ahi,'g_efgint',alo,g_efgint)) call
     &   errquit(pname//': nga_create failed g_efgint',11,GA_ERR)
      call ga_zero(g_efgint)

c     calculate field gradient integrals (same as FC+SD in NMR but
c     without spin)

      call int_giao_1ega(ubas,ubas,g_efgint,'sd+fc',
     &   dbl_mb(k_xyzpt),nat, oskel)

c     comment: The electronic integrals contain the isotropic
c     `contact' piece. However, when constructing the traceless
c     EFG tensor, those contributions disappear. We might as well
c     take them out explicitly

c     calculate FC integrals separately
      ahi(3) = nat
      if (.not.nga_create(MT_DBL,3,ahi,'temp',alo,g_temp)) call
     &   errquit('add_x2c: nga_create failed g_temp',31,GA_ERR)
      call ga_zero(g_temp)

      call int_giao_1ega(ubas,ubas,g_temp,'fc',
     &   dbl_mb(k_xyzpt),nat,oskel)

c
c     Take out FC part from SD+FC (i.e. from xx, yy and zz components)
c
      alo(1) = 1
      ahi(1) = nubas
      alo(2) = 1
      ahi(2) = nubas
      blo(1) = 1
      bhi(1) = nubas
      blo(2) = 1
      bhi(2) = nubas

      do iat = 1, nat
        alo(3) = 1+(iat-1)*6
        ahi(3) = 1+(iat-1)*6
        blo(3) = iat
        bhi(3) = iat
        do i = 1, 3 ! runs over xx, yy, zz
          call nga_add_patch(1.0d0,g_efgint,alo,ahi,
     &       pifac,g_temp,blo,bhi,
     &       g_efgint,alo,ahi)
          alo(3) = alo(3)+1
          ahi(3) = ahi(3)+1
         enddo
      enddo

      if (.not. ga_destroy(g_temp)) call errquit
     &   ('x2c_efg: temp ga corrupt?',31, GA_ERR)

c     the field gradient integrals carry a factor of 3
c     that we don't want:

      call ga_scale(g_efgint, third)

      if (dbg>1 .and.master) write (luout,*) pname//' before iat loop'

c     the following arrays are used throughout the following loop

      if (needu2c) then
        if(.not.ga_create(mt_dbl,nubas,nubas,'wtmp',0,0,g_wtmp))
     &    call errquit(pname//': error creating wtmp',12, GA_ERR)
        
        if(.not.ga_create(mt_dbl,ncbas,nubas,'tmp2',0,0,g_tmp2))
     &    call errquit(pname//': error creating tmp2',148, GA_ERR)
      end if                    ! needu2c

c     ---------------------------------------
c     calculate EFGs in a loop over the atoms
c     ---------------------------------------

      do iat = 1,nat

        xp = dbl_mb(k_xyzpt  +3*(iat-1))
        yp = dbl_mb(k_xyzpt+1+3*(iat-1))
        zp = dbl_mb(k_xyzpt+2+3*(iat-1))

        if (dbg>1 .and.master) write (luout,*) pname//': iat=',iat

c       --------------------
c       nuclear contribution
c       --------------------

        do i = 1,6
          efgnuc(i) = 0d0
        end do

        do i = 1,nat
          xn  = dbl_mb(k_xyzpt  +3*(i-1)) - xp
          yn  = dbl_mb(k_xyzpt+1+3*(i-1)) - yp
          zn  = dbl_mb(k_xyzpt+2+3*(i-1)) - zp
          zan = dbl_mb(k_zanpt+i-1)
          rr = sqrt(xn*xn + yn*yn + zn*zn)
          if (rr.ge.1.0d-3) then 
c           (skip if it is the same nucleus or nuclei too close)
            rr5=rr*rr*rr*rr*rr
            efgnuc(1) = efgnuc(1) - zan*xn*xn/rr5
            efgnuc(2) = efgnuc(2) - zan*yn*yn/rr5
            efgnuc(3) = efgnuc(3) - zan*zn*zn/rr5
            efgnuc(4) = efgnuc(4) - zan*xn*yn/rr5
            efgnuc(5) = efgnuc(5) - zan*xn*zn/rr5
            efgnuc(6) = efgnuc(6) - zan*yn*zn/rr5
          end if 
        end do                  ! i = loop over all other nuclei

c       create traceless version of nuclear contribution to EFG:
        etmp(1) = 2.0d0*efgnuc(1) - efgnuc(2) - efgnuc(3)
        etmp(2) = 2.0d0*efgnuc(2) - efgnuc(1) - efgnuc(3)
        etmp(3) = 2.0d0*efgnuc(3) - efgnuc(1) - efgnuc(2)
        etmp(4) = 3.0d0*efgnuc(4)
        etmp(5) = 3.0d0*efgnuc(5)
        etmp(6) = 3.0d0*efgnuc(6)
        do i = 1,6
          efgnuc(i) = etmp(i)
        end do
        
c       -----------------------
c       electronic contribution
c       -----------------------
        
c       First, we do the case without picture change corrections
        
c       copy electronic EFG integrals into temp array:
c       we also rearrange them such that we get the traceless
c       versions of the integrals right away
        
        alo(1) = nubas
        alo(2) = -1
        alo(3) = -1
        ahi(1) = nubas
        ahi(2) = nubas
        ahi(3) = 6
       
       if (.not.nga_create(MT_DBL,3,ahi,'efg',alo,g_efg)) call
     &     errquit('add_x2c: nga_create failed g_efg',32,GA_ERR)
        call ga_zero(g_efg)
 
        ahi(3) = 1
        if (.not.nga_create(MT_DBL,2,ahi,'temp',alo,g_temp)) call
     &     errquit(pname//': nga_create failed g_temp',31,GA_ERR)
        call ga_zero(g_temp)
               
        alo(1) = 1
        ahi(1) = nubas
        alo(2) = 1
        ahi(2) = nubas

        blo(1) = 1
        bhi(1) = nubas
        blo(2) = 1
        bhi(2) = nubas

        clo(1) = 1
        chi(1) = nubas
        clo(2) = 1
        chi(2) = nubas
        clo(3) = 1
        chi(3) = 1

c       create traceless version of EFG matrices in g_efg
c       for the selected reference point

c       step 1: accumulate 3xx, 3yy, 3zz, 
c       and accumulate xx+yy+zz in g_temp
        do i = 1,3
          alo(3) = (iat-1)*6 + i 
          ahi(3) = alo(3)
          blo(3) = i
          bhi(3) = blo(3)
          call nga_add_patch(3.0d0,g_efgint,alo,ahi,
     &       1d0,g_efg,blo,bhi,
     &       g_efg,blo,bhi)
          call nga_add_patch(1.0d0,g_efgint,alo,ahi,
     &       1d0,g_temp,clo,chi,
     &       g_temp,clo,chi)          
        end do
c       step 2: subtract xx+yy+zz from diagonal 
        do i = 1,3
          blo(3) = i
          bhi(3) = blo(3)
          call nga_add_patch(1.0d0,g_efg,blo,bhi,
     &       -1d0,g_temp,clo,chi,
     &       g_efg,blo,bhi)
        end do

c       indices 1, 2, 3 of g_efg are now the components
c       3xx-rr, 3yy-rr, 3zz-rr

c       step 3: assign 3xy, 3xz, 3yz to the remaining elements 4, 5, 6
        do i = 4,6
          alo(3) = (iat-1)*6 + i 
          ahi(3) = alo(3)
          blo(3) = i
          bhi(3) = blo(3)
          call nga_add_patch(3.0d0,g_efgint,alo,ahi,
     &       1d0,g_efg,blo,bhi,
     &       g_efg,blo,bhi)
        end do

c       calculate electronic EFG components:

        alo(3) = 1 ! for density matrix
        ahi(3) = 1
               
        do i = 1,6              ! loop over EFG components
          
          if (dbg>1 .and.master) write (luout,*) pname//': i=',i

          efgel(i) = 0d0

c         if we use array g_efg, use this:
          blo(3) = i            ! select this section from array g_efg
          bhi(3) = blo(3)

c         if we use g_efgint, use the boundaries below:
c$$$          blo(3) = i + (iat-1)*6
c$$$      bhi(3) = blo(3)

c         tr[P*EFG]

          alo(1) = 1
          ahi(1) = nubas
          alo(2) = 1
          ahi(2) = nubas

          call ga_zero(g_temp)
          call nga_matmul_patch('n','n',1d0,0d0,
     &       g_densu,alo,ahi,
     &       g_efg, blo,bhi,
     &       g_temp,clo,chi)
          sum = ga_trace_diag(g_temp)
          
          efgel(i) = efgel(i) + sum
          
        end do                  ! i = terms in EFG

c        print results:
        
        if (master) then
          write (luout,*)
          write (luout,*) 'X2C EFG: Before picture-change corr.'
          write (luout,'(1x,a,i4,6e15.7)') 'x2cefg: inuc, efgnuc =',
     &       iat, (efgnuc(i), i=1,6)
          write (luout,'(1x,a,i4,6e15.7)') 'x2cefg: inuc, efgel  =',
     &       iat, (efgel(i), i = 1,6)
          write (luout,'(1x,a,i4,6e15.7)') 'x2cefg: inuc, efgtot =',
     &       iat, (efgel(i)+efgnuc(i), i = 1,6)
          write (luout,*)
        end if

c       ---------------------------------------------------
c       calculate picture-change corrected EFG operator
c       for the selected reference point (nuclear position)
c       ---------------------------------------------------

c       allocate arrays to hold the AO matrices of the EFG
c       operator in the small-component basis
c       (we'll calculate all 6 matrices in one go)

        alo(1) = ncbas
        alo(2) = -1
        alo(3) = -1
        ahi(1) = ncbas
        ahi(2) = ncbas
        ahi(3) = 1
        do j = 1,6
          if (.not.nga_create(MT_DBL,2,ahi,'g_w',alo,g_w(j)))
     &       call errquit(
     &       pname//': nga_create failed g_w',1+10*j,GA_ERR)  
          call ga_zero(g_w(j))
        end do
        
        origin(1) = xp
        origin(2) = yp
        origin(3) = zp
        if (dbg>0.and.master) write (luout,*) 'calling getw'
        call x2c_getw_efg_scalar(
     &     rtdb, g_dens_at, origin,
     &     g_w,
     &     nexc)
        if (dbg>0.and.master) write (luout,*) 'return from getw'

c       --------------------
c       allocate operator GA
c       --------------------

        if(.not.ga_create(mt_dbl, nubas, nubas,'g_P',0,0,
     &     g_P))
     &     call errquit(
     &     pname//': error creating g_P',666, GA_ERR)

        
c       ------------------------------------------------------------
c       in a loop over the EFG components, transform the 
c       operator to two-component form, and calculate the electronic
c       EFG component again
c       ------------------------------------------------------------
        
        do i = 1,6              ! loop over EFG components
          
          call ga_zero(g_P)

c         array bounds needed to extract EFG component from GA:

          alo(1) = 1
          ahi(1) = nubas
          alo(2) = 1
          ahi(2) = nubas
          alo(3) = i
          ahi(3) = i

c         g_temp always uses the same dimensions clo, chi
c         we use the same ones for g_u, g_P
          clo(1) = 1
          chi(1) = nubas
          clo(2) = 1
          chi(2) = nubas

c         transform EFG operator

c         UL' EFG UL, store in g_P

          call nga_matmul_patch('t','n',1d0,0d0,
     &      g_u(1),clo,chi,
     &      g_efg,alo,ahi,
     &      g_temp,clo,chi)
          
          call nga_matmul_patch('n','n',1d0,0d0,
     &      g_temp,clo,chi,
     &      g_u(1),clo,chi,
     &      g_P,clo,chi)
          
c         US' pEFGp US, add to g_P
c
c         for US' pEFGp US, we might need to transform g_w(i)
c         to the uncontracted basis first:

          if (needu2c) then

            call ga_zero(g_tmp2)
            call ga_zero (g_wtmp)
            call ga_dgemm('n','n',ncbas,nubas,ncbas, 1.0d0,
     &        g_w(i),g_c2u,0.0d0,g_tmp2)
            call ga_dgemm('t','n',nubas,nubas,ncbas ,1.0d0,
     &        g_c2u, g_tmp2, 1.0d0,g_wtmp)
          else
            g_wtmp = g_w(i)
          end if ! needu2c
        
          
          call nga_matmul_patch('t','n',1d0,0d0,
     &      g_u(2),clo,chi,
     &      g_wtmp,clo,chi,
     &      g_temp,clo,chi)
          
          call nga_matmul_patch('n','n',1d0,1d0,
     &      g_temp,clo,chi,
     &      g_u(2),clo,chi,
     &      g_P,clo,chi)

c          call ga_print(g_P)

c         -----------------------------------------
c         calculate electronic EFG again, this time 
c         using the transformed operator
c         -----------------------------------------

          efgel(i) = 0.0d0
          
          alo(1) = 1                
          ahi(1) = nubas
          alo(2) = 1                
          ahi(2) = nubas
          alo(3) = 1
          ahi(3) = 1
          blo(1) = 1               
          bhi(1) = nubas
          blo(2) = 1
          bhi(2) = nubas
          blo(3) = 1
          bhi(3) = 1
          
          
c         tr[P*EFG]
          call ga_zero(g_temp)
          call nga_matmul_patch('n','n',1d0,0d0,
     &       g_densu,alo,ahi,
     &       g_P,alo,ahi,
     &       g_temp,alo,ahi)
          sum = ga_trace_diag(g_temp)
          efgel(i) = efgel(i) + sum          
          
        end do                  ! i = loop over EFG components

c       print results:
        
        if (master) then
          write (luout,*)
          write (luout,*) 'X2C EFG: After picture-change corr.'
          write (luout,'(1x,a,i4,6e15.7)') 'x2cefg: inuc, efgnuc =',
     &       iat, (efgnuc(i), i=1,6)
          write (luout,'(1x,a,i4,6e15.7)') 'x2cefg: inuc, efgel  =',
     &       iat, (efgel(i), i = 1,6)
          write (luout,'(1x,a,i4,6e15.7)') 'x2cefg: inuc, efgtot =',
     &       iat, (efgel(i)+efgnuc(i), i = 1,6)
          write (luout,*)
        end if

c       --------------------
c       clean up some GA memory
c       --------------------

        if (.not. ga_destroy(g_efg)) call errquit
     &     (pname//': ga corrupt?',32, GA_ERR)

        if (.not. ga_destroy(g_temp)) call errquit
     &     (pname//': ga corrupt?',31, GA_ERR)

        if (.not. ga_destroy(g_P)) call errquit
     &     (pname//': ga P corrupt?',666, GA_ERR)

        do j = 1,6
          if (.not. ga_destroy(g_w(j))) call errquit
     &       (pname//': failed destroying g_w',1+10*j, GA_ERR)
        end do

c       --------------------------------------------------------
c       end loop over refernce points at which EFG is calculated
c       (typically the list of nuclei)
c       --------------------------------------------------------

      end do ! iat, loop over nuclei

      if (.not. ga_destroy(g_efgint)) call errquit
     &  (pname//': ga corrupt?',11, GA_ERR)

      if (needu2c) then
        if (.not. ga_destroy(g_wtmp)) call errquit
     &    (pname//': ga wtmp corrupt?',12, GA_ERR)
        if (.not. ga_destroy(g_tmp2)) call errquit
     &    (pname//': ga tmp2 corrupt?',12, GA_ERR)
        if (.not. ga_destroy(g_densu)) call errquit
     &    (pname//': ga densu corrupt?',58, GA_ERR)
      end if                    ! needu2c
        

c     --------
c     all done
c     --------

      call ga_sync()     

      if (dbg>0.and.master) write(luout,*) 'leaving x2c_efg'

c     =================================================================

      return
      end

