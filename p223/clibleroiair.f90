!
! https://gitlab.com/jrh/wiglaf
!  murchison-35-g785ffede25b-dirty
!
! external/p223/clibleroiair.f90
!
! Copyright (C) 2015 CSIRO
!
! This file is part of wiglaf
!
! wiglaf is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License Version 2
! as published by the Free Software Foundation.
!
! wiglaf is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with wiglaf. If not, see <http://www.gnu.org/licenses/>.
!
! To contact CSIRO about this software you can e-mail
! juerg.hauser@csiro.au
!
! Wiglaf relies on several of the programs that form the P223 suite 
! to solve the forward problem. This file belongs to the version of 
! LeroiAir distributed as part of wiglaf. The original code has been
! re-organised into modules with unique names that can be compiled 
! into a library with a C interface using iso c bindings. It provides
! functions that can be called  to compute the model response for 
! specific AEM systems and is also linked to from a standalone 
! executable that is functionally equivalent to the original LeroiAir.
!
! leroiair.f90 	   - the main program
! libleroiair.f90  - original subroutines reorganised into different modules
! clibleroiair.f90 - functions to compute the model response for specific systems
! libfranken.f90   - leroi routines that are used if the plunge is not zero
! 
! The orignal copy right notice and documentation for LeroiAir is in libleroiair.f90

  module leroiair
  
  contains

  subroutine get_rotation_matrix_rx(alpha,rx)
     use leroiair_subroutines
    IMPLICIT NONE

    real alpha
    real rx(3,3)
    rx=0

    rx(1,1)=1

    rx(2,2)=cos(alpha)
    rx(2,3)=sin(alpha)

    rx(3,2)=-sin(alpha)
    rx(3,3)=cos(alpha)

  end subroutine get_rotation_matrix_rx

  subroutine get_rotation_matrix_ry(beta,ry)
     use leroiair_subroutines
    IMPLICIT NONE

    real beta
    real ry(3,3)
    ry=0

    ry(1,1)=cos(beta)
    ry(1,3)=-sin(beta)
    ry(2,2)=1
    ry(3,1)=sin(beta)
    ry(3,3)=cos(beta)

  end subroutine get_rotation_matrix_ry

  subroutine get_rotation_matrix_rz(gamma,rz)
     use leroiair_subroutines
    IMPLICIT NONE

    real gamma
    real rz(3,3)
    rz=0
    rz(1,1)=cos(gamma)
    rz(1,2)=sin(gamma)
    rz(2,1)=-sin(gamma)
    rz(2,2)=cos(gamma)
    rz(3,3)=1

  end subroutine get_rotation_matrix_rz

  subroutine get_rotation_matrix(plazm,pldip,plunj,r)
     use leroiair_subroutines
    IMPLICIT NONE


    real plazm,pldip,plunj
    real,dimension(3,3) :: rplng,rdip,razm,r,t

    call get_rotation_matrix_rz (plunj,rplng);
    call get_rotation_matrix_rx (pldip,rdip);
    call get_rotation_matrix_rz (plazm,razm);

    t=matmul(rdip,rplng)
    r=matmul(razm,t)

    return


  end subroutine get_rotation_matrix

  SUBROUTINE GET_ORIGIN_SHIFT(XCNTR,YCNTR,NPLT,NSHFT,ESHFT)
     use leroiair_subroutines

    INTEGER NPLT
    REAL XCNTR(NPLT),YCNTR(NPLT)
    REAL ESHFT,NSHFT
    ! Set up body centred origin. Shift origin by (D_NORTH, D_EAST)

    IF (NPLT > 0) THEN
       ESHFT = YCNTR(1)
       NSHFT = XCNTR(1)
       IF (MAX (ABS(ESHFT),ABS(NSHFT)) < 2000) THEN
          NSHFT = 0.0
          ESHFT = 0.0
       END IF
    ELSE  !// JRH
       XCNTR=0.0
       YCNTR=0.0
    END IF


  END SUBROUTINE GET_ORIGIN_SHIFT


  subroutine get_p223_plate_reference_point(nplt,xp,yp,zp,plazm,pldip,plunj,plngth1,plngth2,pwdth1,xcntr,ycntr,pltop)
    IMPLICIT NONE

    integer nplt,i
    real, dimension(nplt) :: xp,yp,zp,plazm,pldip,plunj,pwdth1,plngth1,plngth2,xcntr,ycntr,pltop
    real:: dzm
    real,dimension(3) :: p,pp
    real,dimension(3,3):: r

    do i=1,nplt
       p(1)=(plngth1(i)+plngth2(i))/2.0-plngth1(i)
       p(2)=-pwdth1(i)
       p(3)=0

       dzm=plazm(i)+1.0/2.0*ACOS(-1.0)
       call get_rotation_matrix(dzm,pldip(i),plunj(i),r)
       pp=matmul(r,p)
       xcntr(i)=pp(2)+xp(i)
       ycntr(i)=pp(1)+yp(i)
       pltop(i)=zp(i)-pp(3)
    end do

  end subroutine get_p223_plate_reference_point

  subroutine get_mxb_and_mxab(nplt,plngth,pwdth,cellw,mxb,mxab)

    implicit none

    integer mxb,mxab,nplt,jp,mxabt
    real cellw
    real plngth(nplt),pwdth(nplt)
    real tmpl ,da(nplt),db(nplt)
    integer na(nplt),nb(nplt)
    mxab=4

    do jp=1,nplt
       tmpl = min (plngth(jp), pwdth(jp)) / 2.
       tmpl = min (tmpl, cellw) + .01
       na(jp) = ceiling (plngth(jp) / tmpl) + 1
       nb(jp) = ceiling (pwdth(jp) / tmpl) + 1
       na(jp) = max (2,na(jp))
       nb(jp) = max (2,nb(jp))
       da(jp) = plngth(jp) / real (na(jp),4)
       db(jp) = pwdth(jp) / real (nb(jp),4)
       mxabt = (na(jp)) * (nb(jp))
       mxab = max (mxab, mxabt)
    end do

    mxb = max (2, maxval (nb))

        MXB = CEILING (1.5 * MXB)
       MXAB = CEILING (1.5 * MXAB)

  end subroutine get_mxb_and_mxab



  subroutine formod_tempest_data_and_jacobian(nlyr,res,pbres,thk,  &
       nplt,peast_,pnorth_,ptop,pres,plngth1,plngth2,pwdth1,pwdth2,pthk,pdzm_,pdip,plng, &
       nchnl,topn,tcls,nstat,tx_,ty_,tz,tazi,tincl,rx_,ry_,rz,trdx,trdy,trdz, &
       xmodl,a,ncmp,cmp,cellw,ijac,leroiair_failure_count) &
       bind(c,name="leroiair_formod_tempest_data_and_jacobian")

    ! TEMPEST system forward modelling
    !              nlyr - number of layers
    !         res(nlyr) - resistivity values
    !       thk(nlyr-1) - layer thickness
    !             nchnl - number of time domain channels
    !       topn(nchnl) - time at which receiver channel I opens.
    !       tcls(nchnl) - time at which receiver channel I closes.
    !                tz - transmitter elevation
    !            txincl - transmitter inclination
    !              trdx - in-line receiver offset for each frequency;  behind = positive
    !              trdy - transverse receiver offset for each frequency; Port = positive.
    !              trdz - vertical receiver offset for each frequency; below = positive
    !    xmodl(nchnl*2) - model response vertical component and then inline
    ! a(nchnl*2nlyr*2-1) - Jacobian
    use iso_c_binding
    use leroiair_subroutines
    implicit none
    integer (kind=c_int) nlyr
    real (kind=c_float) res(nlyr*nstat),thk((nlyr-1)*nstat)
    integer (kind=c_int)  nplt
    real (kind=c_float) peast(nplt),pnorth(nplt),ptop(nplt),pres(nplt)
    real (kind=c_float) plngth1(nplt),plngth2(nplt),pwdth1(nplt),pwdth2(nplt),pthk(nplt)
    real (kind=c_float) pdzm(nplt),pdip(nplt),plng(nplt)
    real (kind=c_float) peast_(nplt),pnorth_(nplt),pdzm_(nplt)
    integer (kind=c_int) nchnl
    real (kind=c_float) topn(nchnl),tcls(nchnl)
    integer (kind=c_int)  nstat
    real (kind=c_float) tx(nstat),ty(nstat),tz(nstat),tazi(nstat),tincl(nstat),trdx(nstat),trdy(nstat),trdz(nstat)
    real (kind=c_float) tx_(nstat),ty_(nstat),rx_(nstat),ry_(nstat)
    real (kind=c_float) xmodl(nchnl*ncmp*nstat)

    REAL (kind=c_float) a(nstat*((nlyr*2-1)*(nchnl*ncmp)+(nchnl*ncmp)*12*nplt))
    integer (kind=c_int) leroiair_failure_count

    integer (kind=c_int) ijac(15)
    integer (kind=c_int) ncmp,cmp
    INTEGER NDATA
    INTEGER I
    INTEGER,PARAMETER:: TDFD=1,STEP=1,IDER=4,NSX=1,NTYRP=47,NTYPLS=45,NPULS=1,&
         GSTRP=0,ASTRP=0,KNRM=3
    INTEGER NFRQ
    REAL,ALLOCATABLE :: FREQ(:)
    REAL,PARAMETER :: SWX(1)=0.0,SWY(3)=(/1.,0.,0./),PULSE=2.0E-2
    REAL,PARAMETER :: NORM(3) =(/1E6,1E6,1E6 /)


    REAL RMU(NLYR),REPS(NLYR),CALF(NLYR),CTAU(NLYR),CFREQ(NLYR)
    REAL CALFP(NPLT),CTAUP(NPLT),CFREQP(NPLT)
    REAL (kind=c_float) RX(NSTAT),RY(NSTAT),RZ(NSTAT)
    REAL (kind=c_float) CELLW
    INTEGER MXAB,MXB,MCHNL
    LOGICAL SAME_TX(NSTAT)
    real (kind=c_float) pbres

    REAL,PARAMETER :: TRP(47)= &
         (/6.09752306E-06,   7.38731796E-06,   8.94993991E-06,   1.08430995E-05,&
         1.31367151E-05,   1.59154933E-05,   1.92820607E-05,   2.33607498E-05,&
         2.83021946E-05,   3.42888925E-05,   4.15419418E-05,   5.03292104E-05,&
         6.09752315E-05,   7.38731760E-05,   8.94994009E-05,   1.08430999E-04,&
         1.31367153E-04,   1.59154937E-04,   1.92820604E-04,   2.33607498E-04,&
         2.83021946E-04,   3.42888932E-04,   4.15419403E-04,   5.03292133E-04,&
         6.09752315E-04,   7.38731760E-04,   8.94994009E-04,   1.08431000E-03,&
         1.31367147E-03,   1.59154937E-03,   1.92820607E-03,   2.33607506E-03,&
         2.83021946E-03,   3.42888921E-03,   4.15419415E-03,   5.03292121E-03,&
         6.09752303E-03,   7.38731772E-03,   8.94994009E-03,   1.08431000E-02,&
         1.31367156E-02,   1.59154944E-02,   1.92820616E-02,   2.33607497E-02,&
         2.83021946E-02,   3.42888907E-02,   4.15419415E-02/)

    LOGICAL,PARAMETER :: TXA90=.FALSE.
    LOGICAL,PARAMETER :: JCBN=.TRUE.
    REAL NSHFT,ESHFT

    ! HSBOSS_TD and LEROI_3D use different geometries
    CALL GET_ORIGIN_SHIFT(pnorth_,peast_,NPLT,NSHFT,ESHFT)

    tx=ty_-nshft
    ty=tx_-eshft

    rx=ry_-nshft
    ry=rx_-eshft

    peast=pnorth_-nshft
    pnorth=peast_-eshft

    do i=1,nplt
       pdzm(i)=pdzm_(i)-1.0/2.0*ACOS(-1.0)
    end do
    NDATA=nchnl*ncmp*nstat


    !       MCHNL = total number of readings per station to be inverted(TD or FD)
    !             = NCHNL for time-domain when CMP = 11, 13, 4, 42, 43
    !             = 2* NCHNL for time-domain when CMP = 2
    !             = 3* NCHNL for time-domain when CMP = 3
    if (cmp==2) then
       mchnl=nchnl*2
    else if (cmp==3) then
       mchnl=nchnl*3
    else
       mchnl=nchnl
    end if


    RMU=1.
    REPS=1.
    CALF=1.
    CTAU=0.
    CFREQ=1.

    CALFP=1.0
    CTAUP=0.0
    CFREQP=1.0

    call  get_mxb_and_mxab(nplt,plngth1+plngth2,pwdth1+pwdth2,cellw,mxb,mxab)

    SAME_TX=.TRUE.
    SAME_TX(1)=.FALSE.

    CALL SET_FRQ2(NFRQ,FREQ,TOPN,SWX,NSX)




    call FORJAC2 (nlyr,res,pbres,thk, &
         NPLT,pres,plngth1,plngth2,pwdth1,pwdth2,pthk,peast,pnorth,ptop,pdzm,pdip,plng, &
         NDATA,XMODL,JCBN,A, &
         TDFD, CMP,STEP,IDER,NSX,SWX,SWY,NTYRP,TRP,NPULS,PULSE,NTYPLS, &
         NCMP,NCHNL,MCHNL,TOPN, TCLS,GSTRP,ASTRP,NFRQ,FREQ,KNRM,NORM,TXA90, &
         NSTAT,trdx,trdy,trdz,tx,ty,tz,tincl,tazi,SAME_TX,RX,RY,RZ, &
         RMU,REPS,CALF,CTAU,CFREQ,CALFP,CTAUP,CFREQP, &
         CELLW,MXAB,MXB,ijac,leroiair_failure_count)




    DEALLOCATE(FREQ)


  end subroutine formod_tempest_data_and_jacobian

  subroutine formod_tempest_data(nlyr,res,pbres,thk,  &
       nplt,peast_,pnorth_,ptop,pres,plngth1,plngth2,pwdth1,pwdth2,pthk,pdzm_,pdip,plng,  &
       nchnl,topn,tcls,nstat,tx_,ty_,tz,tazi,tincl,rx_,ry_,rz,trdx,trdy,trdz, &
       xmodl,ncmp,cmp,cellw,leroiair_failure_count)  bind(c,name="leroiair_formod_tempest_data")

    ! TEMPEST system forward modelling
    !              nlyr - number of layers
    !         res(nlyr) - resistivity values
    !       thk(nlyr-1) - layer thickness
    !             nchnl - number of time domain channels
    !       topn(nchnl) - time at which receiver channel I opens.
    !       tcls(nchnl) - time at which receiver channel I closes.
    !                tz - transmitter elevation
    !            txincl - transmitter inclination
    !              trdx - in-line receiver offset for each frequency;  behind = positive
    !              trdy - transverse receiver offset for each frequency; Port = positive.
    !              trdz - vertical receiver offset for each frequency; below = positive
    !    xmodl(nchnl*2) - model response vertical component and then inline
    ! a(nchnl*2nlyr*2-1) - Jacobian
    use iso_c_binding
        use leroiair_subroutines

    implicit none
    integer (kind=c_int) nlyr
    real (kind=c_float) res(nlyr*nstat),thk((nlyr-1)*nstat)
    integer (kind=c_int) nplt
    real (kind=c_float) peast(nplt),pnorth(nplt),ptop(nplt),pres(nplt)
    real (kind=c_float) plngth1(nplt),plngth2(nplt),pwdth1(nplt),pwdth2(nplt),pthk(nplt)
    real (kind=c_float) pdzm(nplt),pdip(nplt),plng(nplt)
    real (kind=c_float) peast_(nplt),pnorth_(nplt),pdzm_(nplt)
    integer (kind=c_int) nchnl
    real (kind=c_float) topn(nchnl),tcls(nchnl)
    integer (kind=c_int) nstat
    real (kind=c_float) tx(nstat),ty(nstat),tz(nstat),tazi(nstat),tincl(nstat),trdx(nstat),trdy(nstat),trdz(nstat)
    real (kind=c_float) tx_(nstat),ty_(nstat),rx_(nstat),ry_(nstat)
    real (kind=c_float) xmodl(nchnl*ncmp*nstat)
    real (kind=c_float) pbres

    REAL (kind=c_float) a(nstat*(nlyr*2-1)*(nchnl*ncmp)+(nchnl*ncmp)*12*nplt)
    integer (kind=c_int) leroiair_failure_count

    integer (kind=c_int) ijac(15)
    integer (kind=c_int) ncmp,cmp
    INTEGER NDATA
    INTEGER I
    INTEGER,PARAMETER:: TDFD=1,STEP=1,IDER=4,NSX=1,NTYRP=47,NTYPLS=45,NPULS=1,&
         GSTRP=0,ASTRP=0,KNRM=3
    INTEGER NFRQ
    REAL,ALLOCATABLE :: FREQ(:)
    REAL,PARAMETER :: SWX(1)=0.0,SWY(3)=(/1.,0.,0./),PULSE=2.0E-2
    REAL,PARAMETER :: NORM(3) =(/1E6,1E6,1E6 /)


    REAL RMU(NLYR),REPS(NLYR),CALF(NLYR),CTAU(NLYR),CFREQ(NLYR)
    REAL CALFP(NPLT),CTAUP(NPLT),CFREQP(NPLT)
    REAL (kind=c_float) RX(NSTAT),RY(NSTAT),RZ(NSTAT)
    REAL (kind=c_float) CELLW
    INTEGER MXAB,MXB,MCHNL
    LOGICAL SAME_TX(NSTAT)
    REAL NSHFT,ESHFT

    REAL,PARAMETER :: TRP(47)= &
         (/6.09752306E-06,   7.38731796E-06,   8.94993991E-06,   1.08430995E-05,&
         1.31367151E-05,   1.59154933E-05,   1.92820607E-05,   2.33607498E-05,&
         2.83021946E-05,   3.42888925E-05,   4.15419418E-05,   5.03292104E-05,&
         6.09752315E-05,   7.38731760E-05,   8.94994009E-05,   1.08430999E-04,&
         1.31367153E-04,   1.59154937E-04,   1.92820604E-04,   2.33607498E-04,&
         2.83021946E-04,   3.42888932E-04,   4.15419403E-04,   5.03292133E-04,&
         6.09752315E-04,   7.38731760E-04,   8.94994009E-04,   1.08431000E-03,&
         1.31367147E-03,   1.59154937E-03,   1.92820607E-03,   2.33607506E-03,&
         2.83021946E-03,   3.42888921E-03,   4.15419415E-03,   5.03292121E-03,&
         6.09752303E-03,   7.38731772E-03,   8.94994009E-03,   1.08431000E-02,&
         1.31367156E-02,   1.59154944E-02,   1.92820616E-02,   2.33607497E-02,&
         2.83021946E-02,   3.42888907E-02,   4.15419415E-02/)

    LOGICAL,PARAMETER :: TXA90=.FALSE.
    LOGICAL,PARAMETER :: JCBN=.FALSE.

    ! HSBOSS_TD and LEROI_3D use different geometries!!!!!!

    CALL GET_ORIGIN_SHIFT(pnorth_,peast_,NPLT,NSHFT,ESHFT)

	ijac=0

    tx=ty_-nshft
    ty=tx_-eshft

    rx=ry_-nshft
    ry=rx_-eshft

    peast=pnorth_-nshft
    pnorth=peast_-eshft



    do i=1,nplt
       pdzm(i)=pdzm_(i)-1.0/2.0*ACOS(-1.0)
    end do

    !NRXST is equal NSTAT for Td systems.

    NDATA=nchnl*ncmp*nstat
    if (cmp==2) then
       mchnl=nchnl*2
    else if (cmp==3) then
       mchnl=nchnl*3
    else
       mchnl=nchnl
    end if

    RMU=1.
    REPS=1.
    CALF=1.
    CTAU=0.
    CFREQ=1.

    CALFP=1.0
    CTAUP=0.0
    CFREQP=1.0


    call  get_mxb_and_mxab(nplt,plngth1+plngth2,pwdth1+pwdth2,cellw,mxb,mxab)



    !NPAR=nlyr*2-1+nplt*9
    SAME_TX=.TRUE.
    SAME_TX(1)=.FALSE.

    CALL SET_FRQ2(NFRQ,FREQ,TOPN,SWX,NSX)

    !print *,"calling leroiair forjac2"
    ! Compute plate conductances

    call FORJAC2 (nlyr,res,pbres,thk, &
         NPLT,pres,plngth1,plngth2,pwdth1,pwdth2,pthk,peast,pnorth,ptop,pdzm,pdip,plng, &
         NDATA,XMODL,JCBN,A, &
         TDFD, CMP,STEP,IDER,NSX,SWX,SWY,NTYRP,TRP,NPULS,PULSE,NTYPLS, &
         NCMP,NCHNL,MCHNL,TOPN, TCLS,GSTRP,ASTRP,NFRQ,FREQ,KNRM,NORM,TXA90, &
         NSTAT,trdx,trdy,trdz,tx,ty,tz,tincl,tazi,SAME_TX,RX,RY,RZ, &
         RMU,REPS,CALF,CTAU,CFREQ,CALFP,CTAUP,CFREQP, &
         CELLW,MXAB,MXB,ijac,leroiair_failure_count)


    DEALLOCATE(FREQ)

  end subroutine formod_tempest_data

  subroutine formod_vtem_max_data_and_jacobian(nlyr,res,pbres,thk,  &
       nplt,peast_,pnorth_,ptop,pres,plngth1,plngth2,pwdth1,pwdth2,pthk,pdzm_,pdip,plng,  &
       ntrn,txarea,nchnl,topn,tcls,nsx,swx,ampt,waveform_,nstat,tx_,ty_,tz,tazi,tincl,rx_,ry_,rz,trdx,trdy,trdz, &
       xmodl,a,ncmp,cmp,cellw,ijac,leroiair_failure_count) bind(c,name="leroiair_formod_vtem_max_data_and_jacobian")

    ! Geotem system forward modelling
    !              nlyr - number of layers
    !         res(nlyr) - resistivity values
    !       thk(nlyr-1) - layer thickness
    !             nchnl - number of time domain channels
    !       topn(nchnl) - time at which receiver channel I opens.
    !       tcls(nchnl) - time at which receiver channel I closes.
    !                tz - transmitter elevation
    !            txincl - transmitter inclination
    !              trdx - in-line receiver offset for each frequency;  behind = positive
    !              trdy - transverse receiver offset for each frequency; Port = positive.
    !              trdz - vertical receiver offset for each frequency; below = positive
    !    xmodl(nchnl*2) - model response vertical component and then inline
    ! a(nchnl*2nlyr*2-1) - Jacobian
    use iso_c_binding
        use leroiair_subroutines

    implicit none
    integer (kind=c_int) nlyr
    integer (kind=c_int) ampt
    real (kind=c_float) res(nlyr*nstat),thk((nlyr-1)*nstat)
    integer (kind=c_int) nplt
    real (kind=c_float) peast(nplt),pnorth(nplt),ptop(nplt),pres(nplt)
    real (kind=c_float) plngth1(nplt),plngth2(nplt),pwdth1(nplt),pwdth2(nplt),pthk(nplt)
    real (kind=c_float) pdzm(nplt),pdip(nplt),plng(nplt)
    real (kind=c_float) peast_(nplt),pnorth_(nplt),pdzm_(nplt)
    integer (kind=c_int) nchnl
    real (kind=c_float) topn(nchnl),tcls(nchnl)
    integer (kind=c_int) nsx
    real (kind=c_float) swx(nsx),waveform(nsx),waveform_(nsx)
    real (kind=c_float) swy(nsx,3)
    integer (kind=c_int) nstat
    real (kind=c_float) tx(nstat),ty(nstat),tz(nstat),tazi(nstat),tincl(nstat),rx(nstat,1),ry(nstat,1),rz(nstat,1)
    real (kind=c_float) trdx(nstat),trdy(nstat),trdz(nstat)
    real (kind=c_float) tx_(nstat),ty_(nstat),rx_(nstat),ry_(nstat)
    real (kind=c_float) xmodl(nchnl*ncmp*nstat)
    real (kind=c_float) pbres
    integer (kind=c_int)ntrn
    real (kind=c_float) txarea
    REAL (kind=c_float) a((nchnl*ncmp)*(nlyr*2-1+12*nplt))
    integer (kind=c_int) leroiair_failure_count

    INTEGER NTYRP
    integer (kind=c_int)  ncmp,cmp
    INTEGER NDATA
    INTEGER I,NTYPLS
    INTEGER,PARAMETER:: TDFD=1,NPULS=5,GSTRP=0,ASTRP=0,KNRM=3
     INTEGER ISW,IDER,STEP
    REAL :: PULSE
    REAL :: BFFAC
    INTEGER :: NFRQ
    REAL,ALLOCATABLE :: FREQ(:)
    REAL,PARAMETER :: NORM(3) =(/1.,1.,1. /)

    REAL RMU(NLYR),REPS(NLYR),CALF(NLYR),CTAU(NLYR),CFREQ(NLYR)
    REAL CALFP(NPLT),CTAUP(NPLT),CFREQP(NPLT)

    REAL (kind=c_float) CELLW
    INTEGER MXAB,MXB,MCHNL
    LOGICAL SAME_TX(NSTAT)
    REAL :: PRM_TD(3)


    REAL,ALLOCATABLE :: TRP(:)

    LOGICAL,PARAMETER :: TXA90=.FALSE.
    LOGICAL,PARAMETER :: JCBN=.TRUE.

    integer (kind=c_int)  ijac(15)
    REAL, PARAMETER :: PI=3.141592654
    REAL NSHFT,ESHFT

     if (ampt==0) then
       STEP=0
       ISW = 1
       IDER = 0
       WAVEFORM = waveform_ * NTRN * TXAREA*1000.0
       BFFAC= 1
    else if (ampt==1) then
       STEP=0
       ISW =130
       IDER =1
       BFFAC=1
       WAVEFORM=waveform_
    endif

    PULSE=SWX(NSX)

    CALL SET_TRP2(NPULS,PULSE,nchnl,TOPN,nsx,SWX,NTYPLS,NTYRP,TRP)
    CALL DCPRM_TD (trdx(1),trdy(1),trdz(1),tincl(1)/PI*180,txarea,PRM_TD)

    ! HSBOSS_TD and LEROI_3D use different geometries!!!!!!
    CALL GET_ORIGIN_SHIFT(pnorth_,peast_,NPLT,NSHFT,ESHFT)

    tx=ty_ -nshft
    ty=tx_-eshft

    rx(:,1)=ry_-nshft
    ry(:,1)=rx_-eshft

    peast=pnorth_-nshft
    pnorth=peast_-eshft

    do i=1,nplt
       pdzm(i)=pdzm_(i)-1.0/2.0*ACOS(-1.0)
    end do

    !NRXST is equal NSTAT for Td systems.

    NDATA=nchnl*ncmp*nstat
    if (cmp==2) then
       mchnl=nchnl*2
    else if (cmp==3) then
       mchnl=nchnl*3
    else
       mchnl=nchnl
    end if

    RMU=1.
    REPS=1.
    CALF=1.
    CTAU=0.
    CFREQ=1.

    CALFP=1.0
    CTAUP=0.0
    CFREQP=1.0

    call  get_mxb_and_mxab(nplt,plngth1+plngth2,pwdth1+pwdth2,cellw,mxb,mxab)

    !NPAR=nlyr*2-1+nplt*9
    SAME_TX=.TRUE.
    SAME_TX(1)=.FALSE.


    CALL SET_FRQ2(NFRQ,FREQ,TOPN,SWX,NSX)
    CALL SET_SOURCE (STEP,ISW,BFFAC,WAVEFORM,NSX,SWX,SWY,PRM_TD)






    call FORJAC2 (nlyr,res,pbres,thk, &
         NPLT,pres,plngth1,plngth2,pwdth1,pwdth2,pthk,peast,pnorth,ptop,pdzm,pdip,plng, &
         NDATA,XMODL,JCBN,A, &
         TDFD, CMP,STEP,IDER,NSX,SWX,SWY,NTYRP,TRP,NPULS,PULSE,NTYPLS, &
         NCMP,NCHNL,MCHNL,TOPN, TCLS,GSTRP,ASTRP,NFRQ,FREQ,KNRM,NORM,TXA90, &
         NSTAT,trdx,trdy,trdz,tx,ty,tz,tincl,tazi,SAME_TX,RX,RY,RZ, &
         RMU,REPS,CALF,CTAU,CFREQ,CALFP,CTAUP,CFREQP, &
         CELLW,MXAB,MXB,ijac,leroiair_failure_count)




    DEALLOCATE(FREQ)

  end subroutine formod_vtem_max_data_and_jacobian

  subroutine formod_vtem_max_data(nlyr,res,pbres,thk,  &
       nplt,peast_,pnorth_,ptop,pres,plngth1,plngth2,pwdth1,pwdth2,pthk,pdzm_,pdip,plng,  &
       ntrn,txarea,nchnl,topn,tcls,nsx,swx,ampt,waveform_,nstat,tx_,ty_,tz,tazi,tincl,rx_,ry_,rz,trdx,trdy,trdz, &
       xmodl,ncmp,cmp,cellw,leroiair_failure_count) bind(c,name="leroiair_formod_vtem_max_data")

    ! Geotem system forward modelling
    !              nlyr - number of layers
    !         res(nlyr) - resistivity values
    !       thk(nlyr-1) - layer thickness
    !             nchnl - number of time domain channels
    !       topn(nchnl) - time at which receiver channel I opens.
    !       tcls(nchnl) - time at which receiver channel I closes.
    !                tz - transmitter elevation
    !            txincl - transmitter inclination
    !              trdx - in-line receiver offset for each frequency;  behind = positive
    !              trdy - transverse receiver offset for each frequency; Port = positive.
    !              trdz - vertical receiver offset for each frequency; below = positive
    !    xmodl(nchnl*2) - model response vertical component and then inline
    ! a(nchnl*2nlyr*2-1) - Jacobian
    use iso_c_binding
        use leroiair_subroutines

    implicit none
    integer (kind=c_int) nlyr
    integer (kind=c_int) ampt
    real (kind=c_float) res(nlyr*nstat),thk((nlyr-1)*nstat)
    integer (kind=c_int) nplt
    real (kind=c_float) peast(nplt),pnorth(nplt),ptop(nplt),pres(nplt)
    real (kind=c_float) plngth1(nplt),plngth2(nplt),pwdth1(nplt),pwdth2(nplt),pthk(nplt)
    real (kind=c_float) pdzm(nplt),pdip(nplt),plng(nplt)
    real (kind=c_float) peast_(nplt),pnorth_(nplt),pdzm_(nplt)
    integer (kind=c_int) nchnl
    real (kind=c_float)  topn(nchnl),tcls(nchnl)
    integer (kind=c_int) nsx
    real (kind=c_float) swx(nsx),waveform(nsx),waveform_(nsx)
    real (kind=c_float) swy(nsx,3)
    integer (kind=c_int) nstat
    real (kind=c_float) tx(nstat),ty(nstat),tz(nstat),tazi(nstat),tincl(nstat)
    real (kind=c_float) rx(nstat),ry(nstat),rz(nstat),trdx(nstat),trdy(nstat),trdz(nstat)
    real (kind=c_float) tx_(nstat),ty_(nstat),rx_(nstat),ry_(nstat)
    real (kind=c_float) xmodl(nchnl*ncmp*nstat)
    real (kind=c_float) pbres

    integer (kind=c_int) ntrn
    real (kind=c_float) txarea

    REAL (kind=c_float) a((nchnl*ncmp)*(nlyr*2-1+12*nplt))
    integer (kind=c_int) leroiair_failure_count

    INTEGER NTYRP
    integer (kind=c_int) ncmp,cmp
    INTEGER NDATA
    INTEGER I,NTYPLS
    INTEGER,PARAMETER:: TDFD=1,NPULS=5,GSTRP=0,ASTRP=0,KNRM=3
        INTEGER ISW,IDER,STEP
    REAL :: PULSE
    REAL :: BFFAC
    INTEGER :: NFRQ
    REAL,ALLOCATABLE :: FREQ(:)
    REAL,PARAMETER :: NORM(3) =(/1.,1.,1. /)

    REAL RMU(NLYR),REPS(NLYR),CALF(NLYR),CTAU(NLYR),CFREQ(NLYR)
    REAL CALFP(NPLT),CTAUP(NPLT),CFREQP(NPLT)
    integer (kind=c_int)  ijac(15)
    REAL (kind=c_float) CELLW
    INTEGER MXAB,MXB,MCHNL
    LOGICAL SAME_TX(NSTAT)
    REAL :: PRM_TD(3)

    REAL,ALLOCATABLE :: TRP(:)

    LOGICAL,PARAMETER :: TXA90=.FALSE.
    LOGICAL,PARAMETER :: JCBN=.FALSE.


    REAL, PARAMETER :: PI=3.141592654
    REAL NSHFT,ESHFT

	ijac=0

     if (ampt==0) then
       STEP=0
       ISW = 1
       IDER = 0
       WAVEFORM = waveform_ * NTRN * TXAREA*1000.0
       BFFAC= 1
    else if (ampt==1) then
       STEP=0
       ISW =130
       IDER =1
       BFFAC=1
       WAVEFORM=waveform_
    endif

    PULSE=SWX(NSX)

    CALL SET_TRP2(NPULS,PULSE,nchnl,TOPN,nsx,SWX,NTYPLS,NTYRP,TRP)
    CALL DCPRM_TD (trdx(1),trdy(1),trdz(1),tincl(1)/PI*180,txarea,PRM_TD)

    ! HSBOSS_TD and LEROI_3D use different geometries!!!!!!
    CALL GET_ORIGIN_SHIFT(pnorth_,peast_,NPLT,NSHFT,ESHFT)

    tx=ty_-nshft
    ty=tx_-eshft

    rx=ry_-nshft
    ry=rx_-eshft

    peast=pnorth_-nshft
    pnorth=peast_-eshft

    do i=1,nplt
       pdzm(i)=pdzm_(i)-1.0/2.0*ACOS(-1.0)
    end do

    !NRXST is equal NSTAT for Td systems.

    NDATA=nchnl*ncmp*nstat
    if (cmp==2) then
       mchnl=nchnl*2
    else if (cmp==3) then
       mchnl=nchnl*3
    else
       mchnl=nchnl
    end if

    RMU=1.
    REPS=1.
    CALF=1.
    CTAU=0.
    CFREQ=1.

    CALFP=1.0
    CTAUP=0.0
    CFREQP=1.0

    call  get_mxb_and_mxab(nplt,plngth1+plngth2,pwdth1+pwdth2,cellw,mxb,mxab)


    SAME_TX=.TRUE.
    SAME_TX(1)=.FALSE.


    CALL SET_FRQ2(NFRQ,FREQ,TOPN,SWX,NSX)
    CALL SET_SOURCE (STEP,ISW,BFFAC,WAVEFORM,NSX,SWX,SWY,PRM_TD)


    call FORJAC2 (nlyr,res,pbres,thk, &
         NPLT,pres,plngth1,plngth2,pwdth1,pwdth2,pthk,peast,pnorth,ptop,pdzm,pdip,plng, &
         NDATA,XMODL,JCBN,A, &
         TDFD, CMP,STEP,IDER,NSX,SWX,SWY,NTYRP,TRP,NPULS,PULSE,NTYPLS, &
         NCMP,NCHNL,MCHNL,TOPN, TCLS,GSTRP,ASTRP,NFRQ,FREQ,KNRM,NORM,TXA90, &
         NSTAT,trdx,trdy,trdz,tx,ty,tz,tincl,tazi,SAME_TX,RX,RY,RZ, &
         RMU,REPS,CALF,CTAU,CFREQ,CALFP,CTAUP,CFREQP, &
         CELLW,MXAB,MXB,ijac,leroiair_failure_count)




    DEALLOCATE(FREQ)
  end subroutine formod_vtem_max_data

  subroutine formod_geotem_data_and_jacobian(nlyr,res,pbres,thk,  &
       nplt,peast_,pnorth_,ptop,pres,plngth1,plngth2,pwdth1,pwdth2,pthk,pdzm_,pdip,plng,  &
       nchnl,topn,tcls,nsx,swx,waveform,nstat,tx_,ty_,tz,tazi,tincl,rx_,ry_,rz,trdx,trdy,trdz, &
       xmodl,a,ncmp,cmp,cellw,ijac,leroiair_failure_count) bind(c,name="leroiair_formod_geotem_data_and_jacobian")

    ! Geotem system forward modelling
    !              nlyr - number of layers
    !         res(nlyr) - resistivity values
    !       thk(nlyr-1) - layer thickness
    !             nchnl - number of time domain channels
    !       topn(nchnl) - time at which receiver channel I opens.
    !       tcls(nchnl) - time at which receiver channel I closes.
    !                tz - transmitter elevation
    !            txincl - transmitter inclination
    !              trdx - in-line receiver offset for each frequency;  behind = positive
    !              trdy - transverse receiver offset for each frequency; Port = positive.
    !              trdz - vertical receiver offset for each frequency; below = positive
    !    xmodl(nchnl*2) - model response vertical component and then inline
    ! a(nchnl*2nlyr*2-1) - Jacobian
    use iso_c_binding
        use leroiair_subroutines
    
    implicit none
    integer (kind=c_int) nlyr
    real (kind=c_float) res(nlyr*nstat),thk((nlyr-1)*nstat)
    integer (kind=c_int) nplt
    real (kind=c_float) peast(nplt),pnorth(nplt),ptop(nplt),pres(nplt)
    real (kind=C_float) plngth1(nplt),plngth2(nplt),pwdth1(nplt),pwdth2(nplt),pthk(nplt)
    real (kind=c_float) pdzm(nplt),pdip(nplt),plng(nplt)
    real (kind=c_float) peast_(nplt),pnorth_(nplt),pdzm_(nplt)
    integer (kind=c_int) nchnl
    real (kind=c_float) topn(nchnl),tcls(nchnl)
    integer (kind=c_int) nsx
    real (kind=c_float) swx(nsx),waveform(nsx)
    real (kind=c_float) swy(nsx,3)
    integer (kind=c_int) nstat
    real (kind=c_float)  tx(nstat),ty(nstat),tz(nstat),tazi(nstat),tincl(nstat),rx(nstat,1),ry(nstat,1),rz(nstat,1)
    real (kind=c_float) trdx(nstat),trdy(nstat),trdz(nstat)
    real (kind=c_float) tx_(nstat),ty_(nstat),rx_(nstat),ry_(nstat)
    real (kind=c_float) xmodl(nchnl*ncmp*nstat)
    real (kind=c_float) pbres

    REAL (kind=c_float) a((nchnl*ncmp)*(nlyr*2-1+12*nplt))
    integer (kind=c_int) leroiair_failure_count

    INTEGER NTYRP
    integer (kind=c_int) ncmp,cmp
    INTEGER NDATA
    INTEGER I,NTYPLS
    INTEGER,PARAMETER:: TDFD=1,STEP=0,IDER=1,NPULS=5,&
         GSTRP=1,ASTRP=0,KNRM=3,ISW=30
    REAL :: PULSE
    REAL,PARAMETER :: BFFAC=1.00000000
    INTEGER :: NFRQ
    REAL,ALLOCATABLE :: FREQ(:)
    REAL,PARAMETER :: NORM(3) =(/1.,1.,1. /)

    REAL RMU(NLYR),REPS(NLYR),CALF(NLYR),CTAU(NLYR),CFREQ(NLYR)
    REAL CALFP(NPLT),CTAUP(NPLT),CFREQP(NPLT)

    REAL (kind=c_float) CELLW
    INTEGER MXAB,MXB,MCHNL
    LOGICAL SAME_TX(NSTAT)
    REAL :: PRM_TD(3)
    REAL,PARAMETER :: TXAREA=1.0


    REAL,ALLOCATABLE :: TRP(:)

    LOGICAL,PARAMETER :: TXA90=.FALSE.
    LOGICAL,PARAMETER :: JCBN=.TRUE.

    integer (kind=c_int) ijac(15)
    REAL, PARAMETER :: PI=3.141592654
    REAL NSHFT,ESHFT

    PULSE=SWX(NSX)

    CALL SET_TRP2(NPULS,PULSE,nchnl,TOPN,nsx,SWX,NTYPLS,NTYRP,TRP)
    CALL DCPRM_TD (trdx(1),trdy(1),trdz(1),tincl(1)/PI*180,txarea,PRM_TD)

    ! HSBOSS_TD and LEROI_3D use different geometries!!!!!!
    CALL GET_ORIGIN_SHIFT(pnorth_,peast_,NPLT,NSHFT,ESHFT)

    tx=ty_ -nshft
    ty=tx_-eshft

    rx(:,1)=ry_-nshft
    ry(:,1)=rx_-eshft

    peast=pnorth_-nshft
    pnorth=peast_-eshft

    do i=1,nplt
       pdzm(i)=pdzm_(i)-1.0/2.0*ACOS(-1.0)
    end do

    !NRXST is equal NSTAT for Td systems.

    NDATA=nchnl*ncmp*nstat
    if (cmp==2) then
       mchnl=nchnl*2
    else if (cmp==3) then
       mchnl=nchnl*3
    else
       mchnl=nchnl
    end if

    RMU=1.
    REPS=1.
    CALF=1.
    CTAU=0.
    CFREQ=1.

    CALFP=1.0
    CTAUP=0.0
    CFREQP=1.0

    call  get_mxb_and_mxab(nplt,plngth1+plngth2,pwdth1+pwdth2,cellw,mxb,mxab)

    !NPAR=nlyr*2-1+nplt*9
    SAME_TX=.TRUE.
    SAME_TX(1)=.FALSE.


    CALL SET_FRQ2(NFRQ,FREQ,TOPN,SWX,NSX)
    CALL SET_SOURCE (STEP,ISW,BFFAC,WAVEFORM,NSX,SWX,SWY,PRM_TD)






    call FORJAC2 (nlyr,res,pbres,thk, &
         NPLT,pres,plngth1,plngth2,pwdth1,pwdth2,pthk,peast,pnorth,ptop,pdzm,pdip,plng, &
         NDATA,XMODL,JCBN,A, &
         TDFD, CMP,STEP,IDER,NSX,SWX,SWY,NTYRP,TRP,NPULS,PULSE,NTYPLS, &
         NCMP,NCHNL,MCHNL,TOPN, TCLS,GSTRP,ASTRP,NFRQ,FREQ,KNRM,NORM,TXA90, &
         NSTAT,trdx,trdy,trdz,tx,ty,tz,tincl,tazi,SAME_TX,RX,RY,RZ, &
         RMU,REPS,CALF,CTAU,CFREQ,CALFP,CTAUP,CFREQP, &
         CELLW,MXAB,MXB,ijac,leroiair_failure_count)




    DEALLOCATE(FREQ)

  end subroutine formod_geotem_data_and_jacobian

  subroutine formod_geotem_data(nlyr,res,pbres,thk,  &
       nplt,peast_,pnorth_,ptop,pres,plngth1,plngth2,pwdth1,pwdth2,pthk,pdzm_,pdip,plng,  &
       nchnl,topn,tcls,nsx,swx,waveform,nstat,tx_,ty_,tz,tazi,tincl,rx_,ry_,rz,trdx,trdy,trdz, &
       xmodl,ncmp,cmp,cellw,leroiair_failure_count) bind(c,name="leroiair_formod_geotem_data")

    ! Geotem system forward modelling
    !              nlyr - number of layers
    !         res(nlyr) - resistivity values
    !       thk(nlyr-1) - layer thickness
    !             nchnl - number of time domain channels
    !       topn(nchnl) - time at which receiver channel I opens.
    !       tcls(nchnl) - time at which receiver channel I closes.
    !                tz - transmitter elevation
    !            txincl - transmitter inclination
    !              trdx - in-line receiver offset for each frequency;  behind = positive
    !              trdy - transverse receiver offset for each frequency; Port = positive.
    !              trdz - vertical receiver offset for each frequency; below = positive
    !    xmodl(nchnl*2) - model response vertical component and then inline
    ! a(nchnl*2nlyr*2-1) - Jacobian
    use iso_c_binding
        use leroiair_subroutines

    implicit none
    integer (kind=c_int) nlyr
    real (kind=c_float) res(nlyr*nstat),thk((nlyr-1)*nstat)
    integer (kind=c_int) nplt
    real (kind=c_float) peast(nplt),pnorth(nplt),ptop(nplt),pres(nplt)
    real (kind=c_float) plngth1(nplt),plngth2(nplt),pwdth1(nplt),pwdth2(nplt),pthk(nplt)
    real (kind=c_float) pdzm(nplt),pdip(nplt),plng(nplt)
    real (kind=c_float) peast_(nplt),pnorth_(nplt),pdzm_(nplt)
    integer (kind=c_int) nchnl
    real (kind=c_float)  topn(nchnl),tcls(nchnl)
    integer (kind=c_int) nsx
    real (kind=c_float) swx(nsx),waveform(nsx)
    real (kind=c_float) swy(nsx,3)
    integer (kind=c_int) nstat
    real (kind=c_float) tx(nstat),ty(nstat),tz(nstat),tazi(nstat),tincl(nstat)
    real (kind=c_float) rx(nstat),ry(nstat),rz(nstat),trdx(nstat),trdy(nstat),trdz(nstat)
    real (kind=c_float) tx_(nstat),ty_(nstat),rx_(nstat),ry_(nstat)
    real (kind=c_float) xmodl(nchnl*ncmp*nstat)
    real (kind=c_float) pbres

    REAL (kind=c_float) a((nchnl*ncmp)*(nlyr*2-1+12*nplt))
    integer (kind=c_int) leroiair_failure_count

    INTEGER NTYRP
    integer (kind=c_int) ncmp,cmp
    INTEGER NDATA
    INTEGER I,NTYPLS
    INTEGER,PARAMETER:: TDFD=1,STEP=0,IDER=1,NPULS=5,&
         GSTRP=1,ASTRP=0,KNRM=3,ISW=30
    REAL :: PULSE
    REAL,PARAMETER :: BFFAC=1.00000000
    INTEGER :: NFRQ
    REAL,ALLOCATABLE :: FREQ(:)
    REAL,PARAMETER :: NORM(3) =(/1.,1.,1. /)

    REAL RMU(NLYR),REPS(NLYR),CALF(NLYR),CTAU(NLYR),CFREQ(NLYR)
    REAL CALFP(NPLT),CTAUP(NPLT),CFREQP(NPLT)
    integer (kind=c_int) ijac(15)
    REAL (kind=c_float) CELLW
    INTEGER MXAB,MXB,MCHNL
    LOGICAL SAME_TX(NSTAT)
    REAL :: PRM_TD(3)
    REAL,PARAMETER :: TXAREA=1.0

    REAL,ALLOCATABLE :: TRP(:)

    LOGICAL,PARAMETER :: TXA90=.FALSE.
    LOGICAL,PARAMETER :: JCBN=.FALSE.


    REAL, PARAMETER :: PI=3.141592654
    REAL NSHFT,ESHFT

	ijac=0
    PULSE=SWX(NSX)

    CALL SET_TRP2(NPULS,PULSE,nchnl,TOPN,nsx,SWX,NTYPLS,NTYRP,TRP)
    CALL DCPRM_TD (trdx(1),trdy(1),trdz(1),tincl(1)/PI*180,txarea,PRM_TD)

    ! HSBOSS_TD and LEROI_3D use different geometries!!!!!!
    CALL GET_ORIGIN_SHIFT(pnorth_,peast_,NPLT,NSHFT,ESHFT)

    tx=ty_-nshft
    ty=tx_-eshft

    rx=ry_-nshft
    ry=rx_-eshft

    peast=pnorth_-nshft
    pnorth=peast_-eshft

    do i=1,nplt
       pdzm(i)=pdzm_(i)-1.0/2.0*ACOS(-1.0)
    end do

    !NRXST is equal NSTAT for Td systems.

    NDATA=nchnl*ncmp*nstat
    if (cmp==2) then
       mchnl=nchnl*2
    else if (cmp==3) then
       mchnl=nchnl*3
    else
       mchnl=nchnl
    end if

    RMU=1.
    REPS=1.
    CALF=1.
    CTAU=0.
    CFREQ=1.

    CALFP=1.0
    CTAUP=0.0
    CFREQP=1.0

    call  get_mxb_and_mxab(nplt,plngth1+plngth2,pwdth1+pwdth2,cellw,mxb,mxab)


    SAME_TX=.TRUE.
    SAME_TX(1)=.FALSE.


    CALL SET_FRQ2(NFRQ,FREQ,TOPN,SWX,NSX)
    CALL SET_SOURCE (STEP,ISW,BFFAC,WAVEFORM,NSX,SWX,SWY,PRM_TD)


    call FORJAC2 (nlyr,res,pbres,thk, &
         NPLT,pres,plngth1,plngth2,pwdth1,pwdth2,pthk,peast,pnorth,ptop,pdzm,pdip,plng, &
         NDATA,XMODL,JCBN,A, &
         TDFD, CMP,STEP,IDER,NSX,SWX,SWY,NTYRP,TRP,NPULS,PULSE,NTYPLS, &
         NCMP,NCHNL,MCHNL,TOPN, TCLS,GSTRP,ASTRP,NFRQ,FREQ,KNRM,NORM,TXA90, &
         NSTAT,trdx,trdy,trdz,tx,ty,tz,tincl,tazi,SAME_TX,RX,RY,RZ, &
         RMU,REPS,CALF,CTAU,CFREQ,CALFP,CTAUP,CFREQP, &
         CELLW,MXAB,MXB,ijac,leroiair_failure_count)

    DEALLOCATE(FREQ)
  end subroutine formod_geotem_data

  subroutine formod_spectrem_data_and_jacobian(nlyr,res,pbres,thk,  &
       nplt,peast_,pnorth_,ptop,pres,plngth1,plngth2,pwdth1,pwdth2,pthk,pdzm_,pdip,plng,  &
        ntrn_,txarea_,nchnl,topn,tcls,nsx,swx,waveform_,nstat,tx_,ty_,tz,tazi,tincl,rx_,ry_,rz,trdx,trdy,trdz, &
       xmodl,a,ncmp,cmp,cellw,ijac,leroiair_failure_count) bind(c,name="leroiair_formod_spectrem_data_and_jacobian")

    ! Geotem system forward modelling
    !              nlyr - number of layers
    !         res(nlyr) - resistivity values
    !       thk(nlyr-1) - layer thickness
    !             nchnl - number of time domain channels
    !       topn(nchnl) - time at which receiver channel I opens.
    !       tcls(nchnl) - time at which receiver channel I closes.
    !                tz - transmitter elevation
    !            txincl - transmitter inclination
    !              trdx - in-line receiver offset for each frequency;  behind = positive
    !              trdy - transverse receiver offset for each frequency; Port = positive.
    !              trdz - vertical receiver offset for each frequency; below = positive
    !    xmodl(nchnl*2) - model response vertical component and then inline
    ! a(nchnl*2nlyr*2-1) - Jacobian
    use iso_c_binding
        use leroiair_subroutines

    implicit none
    integer (kind=c_int)  nlyr
    real (kind=c_float) res(nlyr*nstat),thk((nlyr-1)*nstat)
    integer (kind=c_int) nplt
    real (kind=c_float) peast(nplt),pnorth(nplt),ptop(nplt),pres(nplt)
    real (kind=c_float) plngth1(nplt),plngth2(nplt),pwdth1(nplt),pwdth2(nplt),pthk(nplt)
    real (kind=c_float) pdzm(nplt),pdip(nplt),plng(nplt)
    integer (kind=c_int) ntrn,ntrn_
    real (kind=c_float) txarea,txarea_
    real (kind=c_float) peast_(nplt),pnorth_(nplt),pdzm_(nplt)
    integer (kind=c_int) nchnl
    real (kind=c_float) topn(nchnl),tcls(nchnl)
    integer (kind=c_int) nsx
    real (kind=c_float) swx(nsx),waveform_(nsx),waveform(nsx)
    real (kind=c_float) swy(nsx,3)
    integer (kind=c_int) nstat
    real (kind=c_float) tx(nstat),ty(nstat),tz(nstat),tazi(nstat),tincl(nstat),rx(nstat,1),ry(nstat,1),rz(nstat,1)
    real (kind=c_float) trdx(nstat),trdy(nstat),trdz(nstat)
    real (kind=c_float) tx_(nstat),ty_(nstat),rx_(nstat),ry_(nstat)
    real (kind=c_float) xmodl(nchnl*ncmp*nstat)
    real (kind=c_float) pbres

    REAL (kind=c_float) a((nchnl*ncmp)*(nlyr*2-1+12*nplt))
    integer (kind=c_int) leroiair_failure_count

    INTEGER NTYRP
    integer (kind=c_int) ncmp,cmp
    INTEGER NDATA
    INTEGER I,NTYPLS
    INTEGER,PARAMETER:: TDFD=1,STEP=1,IDER=0,NPULS=5,&
         GSTRP=0,ASTRP=0,KNRM=4,ISW=1

    REAL :: PULSE
    REAL,PARAMETER :: BFFAC=1.00000000
    INTEGER :: NFRQ
    REAL,ALLOCATABLE :: FREQ(:)
    REAL :: NORM(KNRM)

    REAL RMU(NLYR),REPS(NLYR),CALF(NLYR),CTAU(NLYR),CFREQ(NLYR)
    REAL CALFP(NPLT),CTAUP(NPLT),CFREQP(NPLT)

    REAL (kind=c_float) CELLW
    INTEGER MXAB,MXB,MCHNL
    LOGICAL SAME_TX(NSTAT)
    REAL :: PRM_TD(3)


    REAL,ALLOCATABLE :: TRP(:)

    LOGICAL,PARAMETER :: TXA90=.FALSE.
    LOGICAL,PARAMETER :: JCBN=.TRUE.

    integer (kind=c_int) ijac(15)
    REAL, PARAMETER :: PI=3.141592654
    REAL NSHFT,ESHFT
    INTEGER JC
    REAL XRN
      REAL, PARAMETER ::  TOL=1.E-3
    PULSE=SWX(NSX)

    WAVEFORM = waveform_ * NTRN_ * TXAREA_

    CALL SET_TRP2(NPULS,PULSE,nchnl,TOPN,nsx,SWX,NTYPLS,NTYRP,TRP)
    txarea=txarea_
    ntrn=ntrn_

    CALL DCPRM_TD (trdx(1),trdy(1),trdz(1),tincl(1)/PI*180,txarea,PRM_TD)

    ! HSBOSS_TD and LEROI_3D use different geometries!!!!!!
    CALL GET_ORIGIN_SHIFT(pnorth_,peast_,NPLT,NSHFT,ESHFT)

    tx=ty_ -nshft
    ty=tx_-eshft

    rx(:,1)=ry_-nshft
    ry(:,1)=rx_-eshft

    peast=pnorth_-nshft
    pnorth=peast_-eshft

    do i=1,nplt
       pdzm(i)=pdzm_(i)-1.0/2.0*ACOS(-1.0)
    end do

    !NRXST is equal NSTAT for Td systems.

    NDATA=nchnl*ncmp*nstat
    if (cmp==2) then
       mchnl=nchnl*2
    else if (cmp==3) then
       mchnl=nchnl*3
    else
       mchnl=nchnl
    end if

    RMU=1.
    REPS=1.
    CALF=1.
    CTAU=0.
    CFREQ=1.

    CALFP=1.0
    CTAUP=0.0
    CFREQP=1.0

    call  get_mxb_and_mxab(nplt,plngth1+plngth2,pwdth1+pwdth2,cellw,mxb,mxab)

    !NPAR=nlyr*2-1+nplt*9
    SAME_TX=.TRUE.
    SAME_TX(1)=.FALSE.


    CALL SET_FRQ2(NFRQ,FREQ,TOPN,SWX,NSX)

    CALL SET_SOURCE (STEP,ISW,BFFAC,WAVEFORM,NSX,SWX,SWY,PRM_TD)


    NORM(1) = ABS (PRM_TD(1))
    NORM(2) = ABS (PRM_TD(2))
    NORM(3) = ABS (PRM_TD(3))
    NORM(4) = SQRT (NORM(1)**2 + NORM(2)**2 + NORM(3)**2)

    DO JC = 1,3
        XRN = NORM(JC) / NORM (4)
        IF (XRN < TOL) NORM(JC) = NORM(4)
            NORM(JC)=1000000.00/NORM(JC)
    END DO

    call FORJAC2 (nlyr,res,pbres,thk, &
         NPLT,pres,plngth1,plngth2,pwdth1,pwdth2,pthk,peast,pnorth,ptop,pdzm,pdip,plng, &
         NDATA,XMODL,JCBN,A, &
         TDFD, CMP,STEP,IDER,NSX,SWX,SWY,NTYRP,TRP,NPULS,PULSE,NTYPLS, &
         NCMP,NCHNL,MCHNL,TOPN, TCLS,GSTRP,ASTRP,NFRQ,FREQ,KNRM,NORM,TXA90, &
         NSTAT,trdx,trdy,trdz,tx,ty,tz,tincl,tazi,SAME_TX,RX,RY,RZ, &
         RMU,REPS,CALF,CTAU,CFREQ,CALFP,CTAUP,CFREQP, &
         CELLW,MXAB,MXB,ijac,leroiair_failure_count)




    DEALLOCATE(FREQ)

  end subroutine formod_spectrem_data_and_jacobian

   subroutine formod_spectrem_data(nlyr,res,pbres,thk,  &
       nplt,peast_,pnorth_,ptop,pres,plngth1,plngth2,pwdth1,pwdth2,pthk,pdzm_,pdip,plng,  &
        ntrn_,txarea_,nchnl,topn,tcls,nsx,swx,waveform_,nstat,tx_,ty_,tz,tazi,tincl,rx_,ry_,rz,trdx,trdy,trdz, &
       xmodl,ncmp,cmp,cellw,leroiair_failure_count) bind(c,name="leroiair_formod_spectrem_data")

    ! Geotem system forward modelling
    !              nlyr - number of layers
    !         res(nlyr) - resistivity values
    !       thk(nlyr-1) - layer thickness
    !             nchnl - number of time domain channels
    !       topn(nchnl) - time at which receiver channel I opens.
    !       tcls(nchnl) - time at which receiver channel I closes.
    !                tz - transmitter elevation
    !            txincl - transmitter inclination
    !              trdx - in-line receiver offset for each frequency;  behind = positive
    !              trdy - transverse receiver offset for each frequency; Port = positive.
    !              trdz - vertical receiver offset for each frequency; below = positive
    !    xmodl(nchnl*2) - model response vertical component and then inline
    ! a(nchnl*2nlyr*2-1) - Jacobian
    use iso_c_binding
        use leroiair_subroutines

    implicit none
    integer (kind=c_int) nlyr
    real (kind=c_float) res(nlyr*nstat),thk((nlyr-1)*nstat)
    integer(kind=c_int) nplt
    real (kind=c_float) peast(nplt),pnorth(nplt),ptop(nplt),pres(nplt)
    real (kind=c_float) plngth1(nplt),plngth2(nplt),pwdth1(nplt),pwdth2(nplt),pthk(nplt)
    real (kind=c_float) pdzm(nplt),pdip(nplt),plng(nplt)
    integer (kind=c_int) ntrn,ntrn_
    real (kind=c_float)  txarea,txarea_
    real (kind=c_float) peast_(nplt),pnorth_(nplt),pdzm_(nplt)
    integer (kind=c_int) nchnl
    real (kind=c_float) topn(nchnl),tcls(nchnl)
    integer (kind=c_int) nsx
    real (kind=c_float) swx(nsx),waveform_(nsx),waveform(nsx)
    real (kind=c_float) swy(nsx,3)
    integer (kind=c_int) nstat
    real (kind=c_float) tx(nstat),ty(nstat),tz(nstat),tazi(nstat),tincl(nstat),rx(nstat,1),ry(nstat,1),rz(nstat,1)
    real (kind=c_float) trdx(nstat),trdy(nstat),trdz(nstat)
    real (kind=c_float) tx_(nstat),ty_(nstat),rx_(nstat),ry_(nstat)
    real (kind=c_float) xmodl(nchnl*ncmp*nstat)
    real (kind=c_float) pbres

    REAL (kind=c_float) a((nchnl*ncmp)*(nlyr*2-1+12*nplt))
    integer  (kind=c_int) leroiair_failure_count

    INTEGER NTYRP
    integer (kind=c_int) ncmp,cmp
    INTEGER NDATA
    INTEGER I,NTYPLS
    INTEGER,PARAMETER:: TDFD=1,STEP=1,IDER=0,NPULS=5,&
         GSTRP=0,ASTRP=0,KNRM=4,ISW=1

    REAL :: PULSE
    REAL,PARAMETER :: BFFAC=1.00000000
    INTEGER :: NFRQ
    REAL,ALLOCATABLE :: FREQ(:)
    REAL :: NORM(KNRM)

    REAL RMU(NLYR),REPS(NLYR),CALF(NLYR),CTAU(NLYR),CFREQ(NLYR)
    REAL CALFP(NPLT),CTAUP(NPLT),CFREQP(NPLT)

    REAL (kind=c_float) CELLW
    INTEGER MXAB,MXB,MCHNL
    LOGICAL SAME_TX(NSTAT)
    REAL :: PRM_TD(3)


    REAL,ALLOCATABLE :: TRP(:)

    LOGICAL,PARAMETER :: TXA90=.FALSE.
    LOGICAL,PARAMETER :: JCBN=.FALSE.

    integer (kind=c_int) ijac(15)
    REAL, PARAMETER :: PI=3.141592654
    REAL NSHFT,ESHFT
    INTEGER JC
    REAL XRN
      REAL, PARAMETER ::  TOL=1.E-3
      
    ijac=0  
      
    PULSE=SWX(NSX)

    WAVEFORM = waveform_ * NTRN_ * TXAREA_

    CALL SET_TRP2(NPULS,PULSE,nchnl,TOPN,nsx,SWX,NTYPLS,NTYRP,TRP)
    txarea=txarea_
    ntrn=ntrn_

    CALL DCPRM_TD (trdx(1),trdy(1),trdz(1),tincl(1)/PI*180,txarea,PRM_TD)

    ! HSBOSS_TD and LEROI_3D use different geometries!!!!!!
    CALL GET_ORIGIN_SHIFT(pnorth_,peast_,NPLT,NSHFT,ESHFT)

    tx=ty_ -nshft
    ty=tx_-eshft

    rx(:,1)=ry_-nshft
    ry(:,1)=rx_-eshft

    peast=pnorth_-nshft
    pnorth=peast_-eshft

    do i=1,nplt
       pdzm(i)=pdzm_(i)-1.0/2.0*ACOS(-1.0)
    end do

    !NRXST is equal NSTAT for Td systems.

    NDATA=nchnl*ncmp*nstat
    if (cmp==2) then
       mchnl=nchnl*2
    else if (cmp==3) then
       mchnl=nchnl*3
    else
       mchnl=nchnl
    end if

    RMU=1.
    REPS=1.
    CALF=1.
    CTAU=0.
    CFREQ=1.

    CALFP=1.0
    CTAUP=0.0
    CFREQP=1.0

    call  get_mxb_and_mxab(nplt,plngth1+plngth2,pwdth1+pwdth2,cellw,mxb,mxab)

    !NPAR=nlyr*2-1+nplt*9
    SAME_TX=.TRUE.
    SAME_TX(1)=.FALSE.


    CALL SET_FRQ2(NFRQ,FREQ,TOPN,SWX,NSX)

    CALL SET_SOURCE (STEP,ISW,BFFAC,WAVEFORM,NSX,SWX,SWY,PRM_TD)


    NORM(1) = ABS (PRM_TD(1))
    NORM(2) = ABS (PRM_TD(2))
    NORM(3) = ABS (PRM_TD(3))
    NORM(4) = SQRT (NORM(1)**2 + NORM(2)**2 + NORM(3)**2)

    DO JC = 1,3
        XRN = NORM(JC) / NORM (4)
        IF (XRN < TOL) NORM(JC) = NORM(4)
            NORM(JC)=1000000.00/NORM(JC)
    END DO

    call FORJAC2 (nlyr,res,pbres,thk, &
         NPLT,pres,plngth1,plngth2,pwdth1,pwdth2,pthk,peast,pnorth,ptop,pdzm,pdip,plng, &
         NDATA,XMODL,JCBN,A, &
         TDFD, CMP,STEP,IDER,NSX,SWX,SWY,NTYRP,TRP,NPULS,PULSE,NTYPLS, &
         NCMP,NCHNL,MCHNL,TOPN, TCLS,GSTRP,ASTRP,NFRQ,FREQ,KNRM,NORM,TXA90, &
         NSTAT,trdx,trdy,trdz,tx,ty,tz,tincl,tazi,SAME_TX,RX,RY,RZ, &
         RMU,REPS,CALF,CTAU,CFREQ,CALFP,CTAUP,CFREQP, &
         CELLW,MXAB,MXB,ijac,leroiair_failure_count)




    DEALLOCATE(FREQ)

  end subroutine formod_spectrem_data


  subroutine formod_xtem_data_and_jacobian(nlyr,res,pbres,thk,  &
       nplt,peast_,pnorth_,ptop,pres,plngth1,plngth2,pwdth1,pwdth2,pthk,pdzm_,pdip,plng,  &
       ntrn,txarea,nchnl,topn,tcls,nsx,swx,ampt,waveform_,nstat,tx_,ty_,tz,tazi,tincl,rx_,ry_,rz,trdx,trdy,trdz, &
       xmodl,a,ncmp,cmp,cellw,ijac,leroiair_failure_count) bind(c,name="leroiair_formod_xtem_data_and_jacobian")

    ! Geotem system forward modelling
    !              nlyr - number of layers
    !         res(nlyr) - resistivity values
    !       thk(nlyr-1) - layer thickness
    !             nchnl - number of time domain channels
    !       topn(nchnl) - time at which receiver channel I opens.
    !       tcls(nchnl) - time at which receiver channel I closes.
    !                tz - transmitter elevation
    !            txincl - transmitter inclination
    !              trdx - in-line receiver offset for each frequency;  behind = positive
    !              trdy - transverse receiver offset for each frequency; Port = positive.
    !              trdz - vertical receiver offset for each frequency; below = positive
    !    xmodl(nchnl*2) - model response vertical component and then inline
    ! a(nchnl*2nlyr*2-1) - Jacobian
    use iso_c_binding
        use leroiair_subroutines

    implicit none
    integer (kind=c_int) nlyr
    integer (kind=c_int) ampt
    real (kind=c_float) res(nlyr*nstat),thk((nlyr-1)*nstat)
    integer (kind=c_int) nplt
    real (kind=c_float) peast(nplt),pnorth(nplt),ptop(nplt),pres(nplt)
    real (kind=c_float) plngth1(nplt),plngth2(nplt),pwdth1(nplt),pwdth2(nplt),pthk(nplt)
    real (kind=c_float) pdzm(nplt),pdip(nplt),plng(nplt)
    real (kind=c_float) peast_(nplt),pnorth_(nplt),pdzm_(nplt)
    integer (kind=c_int) nchnl
    real (kind=c_float) topn(nchnl),tcls(nchnl)
    integer (kind=c_int) nsx
    real (kind=c_float) swx(nsx),waveform(nsx),waveform_(nsx)
    real (kind=c_float) swy(nsx,3)
    integer (kind=c_int) nstat
    real (kind=c_float) tx(nstat),ty(nstat),tz(nstat),tazi(nstat),tincl(nstat),rx(nstat,1),ry(nstat,1),rz(nstat,1)
    real (kind=c_float) trdx(nstat),trdy(nstat),trdz(nstat)
    real (kind=c_float) tx_(nstat),ty_(nstat),rx_(nstat),ry_(nstat)
    real (kind=c_float) xmodl(nchnl*ncmp*nstat)
    real (kind=c_float) pbres
    integer (kind=c_int) ntrn
    real (kind=c_float) txarea
    REAL (kind=c_float) a((nchnl*ncmp)*(nlyr*2-1+12*nplt))
    integer (kind=c_int) leroiair_failure_count

    INTEGER NTYRP
    integer (kind=c_int) ncmp,cmp
    INTEGER NDATA
    INTEGER I,NTYPLS
    INTEGER,PARAMETER:: TDFD=1,NPULS=5,GSTRP=0,ASTRP=0,KNRM=3
     INTEGER ISW,IDER,STEP
    REAL :: PULSE
    REAL :: BFFAC
    INTEGER :: NFRQ
    REAL,ALLOCATABLE :: FREQ(:)
    REAL,PARAMETER :: NORM(3) =(/1.,1.,1. /)

    REAL RMU(NLYR),REPS(NLYR),CALF(NLYR),CTAU(NLYR),CFREQ(NLYR)
    REAL CALFP(NPLT),CTAUP(NPLT),CFREQP(NPLT)

    REAL (kind=c_float) CELLW
    INTEGER MXAB,MXB,MCHNL
    LOGICAL SAME_TX(NSTAT)
    REAL :: PRM_TD(3)


    REAL,ALLOCATABLE :: TRP(:)

    LOGICAL,PARAMETER :: TXA90=.FALSE.
    LOGICAL,PARAMETER :: JCBN=.TRUE.

    integer (kind=c_int) ijac(15)
    REAL, PARAMETER :: PI=3.141592654
    REAL NSHFT,ESHFT


       STEP=0
       ISW = 1
       IDER = 0
       WAVEFORM = waveform_ * NTRN * TXAREA
       BFFAC= 1.00000000


    PULSE=SWX(NSX)

    CALL SET_TRP2(NPULS,PULSE,nchnl,TOPN,nsx,SWX,NTYPLS,NTYRP,TRP)
    CALL DCPRM_TD (trdx(1),trdy(1),trdz(1),tincl(1)/PI*180,txarea,PRM_TD)

    ! HSBOSS_TD and LEROI_3D use different geometries!!!!!!
    CALL GET_ORIGIN_SHIFT(pnorth_,peast_,NPLT,NSHFT,ESHFT)

    tx=ty_ -nshft
    ty=tx_-eshft

    rx(:,1)=ry_-nshft
    ry(:,1)=rx_-eshft

    peast=pnorth_-nshft
    pnorth=peast_-eshft

    do i=1,nplt
       pdzm(i)=pdzm_(i)-1.0/2.0*ACOS(-1.0)
    end do

    !NRXST is equal NSTAT for Td systems.

    NDATA=nchnl*ncmp*nstat
    if (cmp==2) then
       mchnl=nchnl*2
    else if (cmp==3) then
       mchnl=nchnl*3
    else
       mchnl=nchnl
    end if

    RMU=1.
    REPS=1.
    CALF=1.
    CTAU=0.
    CFREQ=1.

    CALFP=1.0
    CTAUP=0.0
    CFREQP=1.0

    call  get_mxb_and_mxab(nplt,plngth1+plngth2,pwdth1+pwdth2,cellw,mxb,mxab)

    !NPAR=nlyr*2-1+nplt*9
    SAME_TX=.TRUE.
    SAME_TX(1)=.FALSE.


    CALL SET_FRQ2(NFRQ,FREQ,TOPN,SWX,NSX)
    CALL SET_SOURCE (STEP,ISW,BFFAC,WAVEFORM,NSX,SWX,SWY,PRM_TD)

    call FORJAC2 (nlyr,res,pbres,thk, &
         NPLT,pres,plngth1,plngth2,pwdth1,pwdth2,pthk,peast,pnorth,ptop,pdzm,pdip,plng, &
         NDATA,XMODL,JCBN,A, &
         TDFD, CMP,STEP,IDER,NSX,SWX,SWY,NTYRP,TRP,NPULS,PULSE,NTYPLS, &
         NCMP,NCHNL,MCHNL,TOPN, TCLS,GSTRP,ASTRP,NFRQ,FREQ,KNRM,NORM,TXA90, &
         NSTAT,trdx,trdy,trdz,tx,ty,tz,tincl,tazi,SAME_TX,RX,RY,RZ, &
         RMU,REPS,CALF,CTAU,CFREQ,CALFP,CTAUP,CFREQP, &
         CELLW,MXAB,MXB,ijac,leroiair_failure_count)

    DEALLOCATE(FREQ)

  end subroutine formod_xtem_data_and_jacobian

  subroutine formod_xtem_data(nlyr,res,pbres,thk,  &
       nplt,peast_,pnorth_,ptop,pres,plngth1,plngth2,pwdth1,pwdth2,pthk,pdzm_,pdip,plng,  &
       ntrn,txarea,nchnl,topn,tcls,nsx,swx,ampt,waveform_,nstat,tx_,ty_,tz,tazi,tincl,rx_,ry_,rz,trdx,trdy,trdz, &
       xmodl,ncmp,cmp,cellw,leroiair_failure_count) bind(c,name="leroiair_formod_xtem_data")

    ! Geotem system forward modelling
    !              nlyr - number of layers
    !         res(nlyr) - resistivity values
    !       thk(nlyr-1) - layer thickness
    !             nchnl - number of time domain channels
    !       topn(nchnl) - time at which receiver channel I opens.
    !       tcls(nchnl) - time at which receiver channel I closes.
    !                tz - transmitter elevation
    !            txincl - transmitter inclination
    !              trdx - in-line receiver offset for each frequency;  behind = positive
    !              trdy - transverse receiver offset for each frequency; Port = positive.
    !              trdz - vertical receiver offset for each frequency; below = positive
    !    xmodl(nchnl*2) - model response vertical component and then inline
    ! a(nchnl*2nlyr*2-1) - Jacobian
    use iso_c_binding
        use leroiair_subroutines

    implicit none
    integer (kind=c_int) nlyr
    integer (kind=c_int) ampt
    real (kind=c_float) res(nlyr*nstat),thk((nlyr-1)*nstat)
    integer (kind=c_int) nplt
    real (kind=c_float) peast(nplt),pnorth(nplt),ptop(nplt),pres(nplt)
    real (kind=c_float) plngth1(nplt),plngth2(nplt),pwdth1(nplt),pwdth2(nplt),pthk(nplt)
    real (kind=c_float) pdzm(nplt),pdip(nplt),plng(nplt)
    real (kind=c_float) peast_(nplt),pnorth_(nplt),pdzm_(nplt)
    integer (kind=c_int) nchnl
    real (kind=c_float) topn(nchnl),tcls(nchnl)
    integer (kind=c_int) nsx
    real (kind=c_float) swx(nsx),waveform(nsx),waveform_(nsx)
    real (kind=c_float) swy(nsx,3)
    integer (kind=c_int) nstat
    real (kind=c_float) tx(nstat),ty(nstat),tz(nstat),tazi(nstat),tincl(nstat)
    real(kind=c_float) rx(nstat),ry(nstat),rz(nstat),trdx(nstat),trdy(nstat),trdz(nstat)
    real (kind=c_float) tx_(nstat),ty_(nstat),rx_(nstat),ry_(nstat)
    real (kind=c_float) xmodl(nchnl*ncmp*nstat)
    real (kind=c_float) pbres

    integer (kind=c_int) ntrn
    real (kind=c_float) txarea

    REAL (kind=c_float)  a((nchnl*ncmp)*(nlyr*2-1+12*nplt))
    integer (kind=c_int) leroiair_failure_count

    INTEGER NTYRP
    integer (kind=c_int) ncmp,cmp
    INTEGER NDATA
    INTEGER I,NTYPLS
    INTEGER,PARAMETER:: TDFD=1,NPULS=5,GSTRP=0,ASTRP=0,KNRM=3
        INTEGER ISW,IDER,STEP
    REAL :: PULSE
    REAL :: BFFAC
    INTEGER :: NFRQ
    REAL,ALLOCATABLE :: FREQ(:)
    REAL,PARAMETER :: NORM(3) =(/1.,1.,1. /)

    REAL RMU(NLYR),REPS(NLYR),CALF(NLYR),CTAU(NLYR),CFREQ(NLYR)
    REAL CALFP(NPLT),CTAUP(NPLT),CFREQP(NPLT)
    integer (kind=c_int) ijac(15)
    REAL (kind=c_float) CELLW
    INTEGER MXAB,MXB,MCHNL
    LOGICAL SAME_TX(NSTAT)
    REAL :: PRM_TD(3)

    REAL,ALLOCATABLE :: TRP(:)
    LOGICAL,PARAMETER :: TXA90=.FALSE.
    LOGICAL,PARAMETER :: JCBN=.FALSE.


    REAL, PARAMETER :: PI=3.141592654
    REAL NSHFT,ESHFT

	ijac=0


       STEP=0
       ISW = 1
       IDER = 0
       WAVEFORM = waveform_ * NTRN * TXAREA
       BFFAC= 1.00000000


    PULSE=SWX(NSX)

    CALL SET_TRP2(NPULS,PULSE,nchnl,TOPN,nsx,SWX,NTYPLS,NTYRP,TRP)
    CALL DCPRM_TD (trdx(1),trdy(1),trdz(1),tincl(1)/PI*180,txarea,PRM_TD)

    ! HSBOSS_TD and LEROI_3D use different geometries!!!!!!
    CALL GET_ORIGIN_SHIFT(pnorth_,peast_,NPLT,NSHFT,ESHFT)

    tx=ty_-nshft
    ty=tx_-eshft

    rx=ry_-nshft
    ry=rx_-eshft

    peast=pnorth_-nshft
    pnorth=peast_-eshft

    do i=1,nplt
       pdzm(i)=pdzm_(i)-1.0/2.0*ACOS(-1.0)
    end do

    !NRXST is equal NSTAT for Td systems.

    NDATA=nchnl*ncmp*nstat
    if (cmp==2) then
       mchnl=nchnl*2
    else if (cmp==3) then
       mchnl=nchnl*3
    else
       mchnl=nchnl
    end if

    RMU=1.
    REPS=1.
    CALF=1.
    CTAU=0.
    CFREQ=1.

    CALFP=1.0
    CTAUP=0.0
    CFREQP=1.0

    call  get_mxb_and_mxab(nplt,plngth1+plngth2,pwdth1+pwdth2,cellw,mxb,mxab)

    SAME_TX=.TRUE.
    SAME_TX(1)=.FALSE.

    CALL SET_FRQ2(NFRQ,FREQ,TOPN,SWX,NSX)
    CALL SET_SOURCE (STEP,ISW,BFFAC,WAVEFORM,NSX,SWX,SWY,PRM_TD)

    call FORJAC2 (nlyr,res,pbres,thk, &
         NPLT,pres,plngth1,plngth2,pwdth1,pwdth2,pthk,peast,pnorth,ptop,pdzm,pdip,plng, &
         NDATA,XMODL,JCBN,A, &
         TDFD, CMP,STEP,IDER,NSX,SWX,SWY,NTYRP,TRP,NPULS,PULSE,NTYPLS, &
         NCMP,NCHNL,MCHNL,TOPN, TCLS,GSTRP,ASTRP,NFRQ,FREQ,KNRM,NORM,TXA90, &
         NSTAT,trdx,trdy,trdz,tx,ty,tz,tincl,tazi,SAME_TX,RX,RY,RZ, &
         RMU,REPS,CALF,CTAU,CFREQ,CALFP,CTAUP,CFREQP, &
         CELLW,MXAB,MXB,ijac,leroiair_failure_count)

    DEALLOCATE(FREQ)
  end subroutine formod_xtem_data

  subroutine formod_skytem_data_and_jacobian(nlyr,res,pbres,thk,  &
       nplt,peast_,pnorth_,ptop,pres,plngth1,plngth2,pwdth1,pwdth2,pthk,pdzm_,pdip,plng,  &
       ntrn,txarea,nchnl,topn,tcls,nsx,swx,ampt,waveform_,nstat,tx_,ty_,tz,tazi,tincl,rx_,ry_,rz,trdx,trdy,trdz, &
       xmodl,a,ncmp,cmp,cellw,ijac,leroiair_failure_count) bind(c,name="leroiair_formod_skytem_data_and_jacobian")

    ! Geotem system forward modelling
    !              nlyr - number of layers
    !         res(nlyr) - resistivity values
    !       thk(nlyr-1) - layer thickness
    !             nchnl - number of time domain channels
    !       topn(nchnl) - time at which receiver channel I opens.
    !       tcls(nchnl) - time at which receiver channel I closes.
    !                tz - transmitter elevation
    !            txincl - transmitter inclination
    !              trdx - in-line receiver offset for each frequency;  behind = positive
    !              trdy - transverse receiver offset for each frequency; Port = positive.
    !              trdz - vertical receiver offset for each frequency; below = positive
    !    xmodl(nchnl*2) - model response vertical component and then inline
    ! a(nchnl*2nlyr*2-1) - Jacobian
    use iso_c_binding
        use leroiair_subroutines

    implicit none
    integer (kind=c_int) nlyr
    integer (kind=c_int) ampt
    real (kind=c_float) res(nlyr*nstat),thk((nlyr-1)*nstat)
    integer (kind=c_int)  nplt
    real (kind=c_float) peast(nplt),pnorth(nplt),ptop(nplt),pres(nplt)
    real (kind=c_float) plngth1(nplt),plngth2(nplt),pwdth1(nplt),pwdth2(nplt),pthk(nplt)
    real (kind=c_float) pdzm(nplt),pdip(nplt),plng(nplt)
    real (kind=c_float) peast_(nplt),pnorth_(nplt),pdzm_(nplt)
    integer (kind=c_int) nchnl
    real (kind=c_float) topn(nchnl),tcls(nchnl)
    integer (kind=c_int) nsx
    real (kind=c_float) swx(nsx),waveform(nsx),waveform_(nsx)
    real (kind=c_float) swy(nsx,3)
    integer (kind=c_int) nstat
    real (kind=c_float) tx(nstat),ty(nstat),tz(nstat),tazi(nstat),tincl(nstat),rx(nstat,1),ry(nstat,1),rz(nstat,1)
    real (kind=c_float) trdx(nstat),trdy(nstat),trdz(nstat)
    real (kind=c_float) tx_(nstat),ty_(nstat),rx_(nstat),ry_(nstat)
    real (kind=c_float) xmodl(nchnl*ncmp*nstat)
    real (kind=c_float) pbres
    integer (kind=c_int) ntrn
    real (kind=c_float) txarea
    REAL (kind=c_float) a((nchnl*ncmp)*(nlyr*2-1+12*nplt))
    integer (kind=c_int)leroiair_failure_count

    INTEGER NTYRP
    integer (kind=c_int) ncmp,cmp
    INTEGER NDATA
    INTEGER I,NTYPLS
    INTEGER,PARAMETER:: TDFD=1,NPULS=5,GSTRP=0,ASTRP=0,KNRM=3
     INTEGER ISW,IDER,STEP
    REAL :: PULSE
    REAL :: BFFAC
    INTEGER :: NFRQ
    REAL,ALLOCATABLE :: FREQ(:)
    REAL,PARAMETER :: NORM(3) =(/1.,1.,1. /)

    REAL RMU(NLYR),REPS(NLYR),CALF(NLYR),CTAU(NLYR),CFREQ(NLYR)
    REAL CALFP(NPLT),CTAUP(NPLT),CFREQP(NPLT)

    REAL (kind=c_float) CELLW
    INTEGER MXAB,MXB,MCHNL
    LOGICAL SAME_TX(NSTAT)
    REAL :: PRM_TD(3)


    REAL,ALLOCATABLE :: TRP(:)

    LOGICAL,PARAMETER :: TXA90=.FALSE.
    LOGICAL,PARAMETER :: JCBN=.TRUE.

    integer (kind=c_int) ijac(15)
    REAL, PARAMETER :: PI=3.141592654
    REAL NSHFT,ESHFT

   if (ampt==1) then
       STEP=0
       ISW = 1
       IDER = 0
       WAVEFORM = waveform_ * NTRN * TXAREA
       BFFAC= 1.00000000
    endif

    PULSE=SWX(NSX)

    CALL SET_TRP2(NPULS,PULSE,nchnl,TOPN,nsx,SWX,NTYPLS,NTYRP,TRP)
    CALL DCPRM_TD (trdx(1),trdy(1),trdz(1),tincl(1)/PI*180,txarea,PRM_TD)

    ! HSBOSS_TD and LEROI_3D use different geometries!!!!!!
    CALL GET_ORIGIN_SHIFT(pnorth_,peast_,NPLT,NSHFT,ESHFT)

    tx=ty_ -nshft
    ty=tx_-eshft

    rx(:,1)=ry_-nshft
    ry(:,1)=rx_-eshft

    peast=pnorth_-nshft
    pnorth=peast_-eshft

    do i=1,nplt
       pdzm(i)=pdzm_(i)-1.0/2.0*ACOS(-1.0)
    end do

    !NRXST is equal NSTAT for Td systems.

    NDATA=nchnl*ncmp*nstat
    if (cmp==2) then
       mchnl=nchnl*2
    else if (cmp==3) then
       mchnl=nchnl*3
    else
       mchnl=nchnl
    end if

    RMU=1.
    REPS=1.
    CALF=1.
    CTAU=0.
    CFREQ=1.

    CALFP=1.0
    CTAUP=0.0
    CFREQP=1.0

    call  get_mxb_and_mxab(nplt,plngth1+plngth2,pwdth1+pwdth2,cellw,mxb,mxab)

    !NPAR=nlyr*2-1+nplt*9
    SAME_TX=.TRUE.
    SAME_TX(1)=.FALSE.


    CALL SET_FRQ2(NFRQ,FREQ,TOPN,SWX,NSX)
    CALL SET_SOURCE (STEP,ISW,BFFAC,WAVEFORM,NSX,SWX,SWY,PRM_TD)

    call FORJAC2 (nlyr,res,pbres,thk, &
         NPLT,pres,plngth1,plngth2,pwdth1,pwdth2,pthk,peast,pnorth,ptop,pdzm,pdip,plng, &
         NDATA,XMODL,JCBN,A, &
         TDFD, CMP,STEP,IDER,NSX,SWX,SWY,NTYRP,TRP,NPULS,PULSE,NTYPLS, &
         NCMP,NCHNL,MCHNL,TOPN, TCLS,GSTRP,ASTRP,NFRQ,FREQ,KNRM,NORM,TXA90, &
         NSTAT,trdx,trdy,trdz,tx,ty,tz,tincl,tazi,SAME_TX,RX,RY,RZ, &
         RMU,REPS,CALF,CTAU,CFREQ,CALFP,CTAUP,CFREQP, &
         CELLW,MXAB,MXB,ijac,leroiair_failure_count)

    DEALLOCATE(FREQ)

  end subroutine formod_skytem_data_and_jacobian

  subroutine formod_skytem_data(nlyr,res,pbres,thk,  &
       nplt,peast_,pnorth_,ptop,pres,plngth1,plngth2,pwdth1,pwdth2,pthk,pdzm_,pdip,plng,  &
       ntrn,txarea,nchnl,topn,tcls,nsx,swx,ampt,waveform_,nstat,tx_,ty_,tz,tazi,tincl,rx_,ry_,rz,trdx,trdy,trdz, &
       xmodl,ncmp,cmp,cellw,leroiair_failure_count) bind(c,name="leroiair_formod_skytem_data")

    ! Geotem system forward modelling
    !              nlyr - number of layers
    !         res(nlyr) - resistivity values
    !       thk(nlyr-1) - layer thickness
    !             nchnl - number of time domain channels
    !       topn(nchnl) - time at which receiver channel I opens.
    !       tcls(nchnl) - time at which receiver channel I closes.
    !                tz - transmitter elevation
    !            txincl - transmitter inclination
    !              trdx - in-line receiver offset for each frequency;  behind = positive
    !              trdy - transverse receiver offset for each frequency; Port = positive.
    !              trdz - vertical receiver offset for each frequency; below = positive
    !    xmodl(nchnl*2) - model response vertical component and then inline
    ! a(nchnl*2nlyr*2-1) - Jacobian
    use iso_c_binding
        use leroiair_subroutines

    implicit none
    integer (kind=c_int) nlyr
    integer (kind=c_int)  ampt
    real (kind=c_float) res(nlyr*nstat),thk((nlyr-1)*nstat)
    integer (kind=c_int) nplt
    real (kind=c_float) peast(nplt),pnorth(nplt),ptop(nplt),pres(nplt)
    real(kind=c_float) plngth1(nplt),plngth2(nplt),pwdth1(nplt),pwdth2(nplt),pthk(nplt)
    real (kind=c_float) pdzm(nplt),pdip(nplt),plng(nplt)
    real (kind=c_float) peast_(nplt),pnorth_(nplt),pdzm_(nplt)
    integer (kind=c_int) nchnl
    real (kind=c_float) topn(nchnl),tcls(nchnl)
    integer (kind=c_int) nsx
    real (kind=c_float) swx(nsx),waveform(nsx),waveform_(nsx)
    real (kind=c_float) swy(nsx,3)
    integer (kind=c_int)  nstat
    real (kind=c_float) tx(nstat),ty(nstat),tz(nstat),tazi(nstat),tincl(nstat)
    real (kind=c_float) rx(nstat),ry(nstat),rz(nstat),trdx(nstat),trdy(nstat),trdz(nstat)
    real (kind=c_float) tx_(nstat),ty_(nstat),rx_(nstat),ry_(nstat)
    real (kind=c_float) xmodl(nchnl*ncmp*nstat)
    real (kind=c_float) pbres

    integer (kind=c_int) ntrn
    real (kind=c_float) txarea

    REAL (kind=c_float)  a((nchnl*ncmp)*(nlyr*2-1+12*nplt))
    integer (kind=c_int) leroiair_failure_count

    INTEGER NTYRP
    integer (kind=c_int) ncmp,cmp
    INTEGER NDATA
    INTEGER I,NTYPLS
    INTEGER,PARAMETER:: TDFD=1,NPULS=5,GSTRP=0,ASTRP=0,KNRM=3
        INTEGER ISW,IDER,STEP
    REAL :: PULSE
    REAL :: BFFAC
    INTEGER :: NFRQ
    REAL,ALLOCATABLE :: FREQ(:)
    REAL,PARAMETER :: NORM(3) =(/1.,1.,1. /)

    REAL RMU(NLYR),REPS(NLYR),CALF(NLYR),CTAU(NLYR),CFREQ(NLYR)
    REAL CALFP(NPLT),CTAUP(NPLT),CFREQP(NPLT)
    integer (kind=c_int) ijac(15)
    REAL (kind=c_float )CELLW
    INTEGER MXAB,MXB,MCHNL
    LOGICAL SAME_TX(NSTAT)
    REAL :: PRM_TD(3)

    REAL,ALLOCATABLE :: TRP(:)
    LOGICAL,PARAMETER :: TXA90=.FALSE.
    LOGICAL,PARAMETER :: JCBN=.FALSE.


    REAL, PARAMETER :: PI=3.141592654
    REAL NSHFT,ESHFT

	ijac=0

   if (ampt==1) then
       STEP=0
       ISW = 1
       IDER = 0
       WAVEFORM = waveform_ * NTRN * TXAREA
       BFFAC= 1.00000000
    endif

    PULSE=SWX(NSX)

    CALL SET_TRP2(NPULS,PULSE,nchnl,TOPN,nsx,SWX,NTYPLS,NTYRP,TRP)
    CALL DCPRM_TD (trdx(1),trdy(1),trdz(1),tincl(1)/PI*180,txarea,PRM_TD)

    ! HSBOSS_TD and LEROI_3D use different geometries!!!!!!
    CALL GET_ORIGIN_SHIFT(pnorth_,peast_,NPLT,NSHFT,ESHFT)

    tx=ty_-nshft
    ty=tx_-eshft

    rx=ry_-nshft
    ry=rx_-eshft

    peast=pnorth_-nshft
    pnorth=peast_-eshft

    do i=1,nplt
       pdzm(i)=pdzm_(i)-1.0/2.0*ACOS(-1.0)
    end do

    !NRXST is equal NSTAT for Td systems.

    NDATA=nchnl*ncmp*nstat
    if (cmp==2) then
       mchnl=nchnl*2
    else if (cmp==3) then
       mchnl=nchnl*3
    else
       mchnl=nchnl
    end if

    RMU=1.
    REPS=1.
    CALF=1.
    CTAU=0.
    CFREQ=1.

    CALFP=1.0
    CTAUP=0.0
    CFREQP=1.0

    call  get_mxb_and_mxab(nplt,plngth1+plngth2,pwdth1+pwdth2,cellw,mxb,mxab)

    SAME_TX=.TRUE.
    SAME_TX(1)=.FALSE.

    CALL SET_FRQ2(NFRQ,FREQ,TOPN,SWX,NSX)
    CALL SET_SOURCE (STEP,ISW,BFFAC,WAVEFORM,NSX,SWX,SWY,PRM_TD)

    call FORJAC2 (nlyr,res,pbres,thk, &
         NPLT,pres,plngth1,plngth2,pwdth1,pwdth2,pthk,peast,pnorth,ptop,pdzm,pdip,plng, &
         NDATA,XMODL,JCBN,A, &
         TDFD, CMP,STEP,IDER,NSX,SWX,SWY,NTYRP,TRP,NPULS,PULSE,NTYPLS, &
         NCMP,NCHNL,MCHNL,TOPN, TCLS,GSTRP,ASTRP,NFRQ,FREQ,KNRM,NORM,TXA90, &
         NSTAT,trdx,trdy,trdz,tx,ty,tz,tincl,tazi,SAME_TX,RX,RY,RZ, &
         RMU,REPS,CALF,CTAU,CFREQ,CALFP,CTAUP,CFREQP, &
         CELLW,MXAB,MXB,ijac,leroiair_failure_count)

    DEALLOCATE(FREQ)
  end subroutine formod_skytem_data


  subroutine formod_helitem_data_and_jacobian(nlyr,res,pbres,thk,  &
       nplt,peast_,pnorth_,ptop,pres,plngth1,plngth2,pwdth1,pwdth2,pthk,pdzm_,pdip,plng,  &
       ntrn,txarea,nchnl,topn,tcls,nsx,swx,ampt,waveform_,nstat,tx_,ty_,tz,tazi,tincl,rx_,ry_,rz,trdx,trdy,trdz, &
       xmodl,a,ncmp,cmp,cellw,ijac,leroiair_failure_count) bind(c,name="leroiair_formod_helitem_data_and_jacobian")

    ! Geotem system forward modelling
    !              nlyr - number of layers
    !         res(nlyr) - resistivity values
    !       thk(nlyr-1) - layer thickness
    !             nchnl - number of time domain channels
    !       topn(nchnl) - time at which receiver channel I opens.
    !       tcls(nchnl) - time at which receiver channel I closes.
    !                tz - transmitter elevation
    !            txincl - transmitter inclination
    !              trdx - in-line receiver offset for each frequency;  behind = positive
    !              trdy - transverse receiver offset for each frequency; Port = positive.
    !              trdz - vertical receiver offset for each frequency; below = positive
    !    xmodl(nchnl*2) - model response vertical component and then inline
    ! a(nchnl*2nlyr*2-1) - Jacobian
    use iso_c_binding
        use leroiair_subroutines
    !   use print_leroiair_vars

    implicit none
    integer (kind=c_int)  nlyr
    integer (kind=c_int) ampt
    real (kind=c_float) res(nlyr*nstat),thk((nlyr-1)*nstat)
    integer (kind=c_int) nplt
    real (kind=c_float) peast(nplt),pnorth(nplt),ptop(nplt),pres(nplt)
    real(kind=c_float) plngth1(nplt),plngth2(nplt),pwdth1(nplt),pwdth2(nplt),pthk(nplt)
    real (kind=c_float) pdzm(nplt),pdip(nplt),plng(nplt)
    real (kind=c_float) peast_(nplt),pnorth_(nplt),pdzm_(nplt)
    integer (kind=c_int)  nchnl
    real (kind=c_float) topn(nchnl),tcls(nchnl)
    integer (kind=c_int) nsx
    real (kind=c_float) swx(nsx),waveform(nsx),waveform_(nsx)
    real (kind=c_float) swy(nsx,3)
    integer (kind=c_int) nstat
    real (kind=c_float) tx(nstat),ty(nstat),tz(nstat),tazi(nstat),tincl(nstat),rx(nstat,1),ry(nstat,1),rz(nstat,1)
    real (kind=c_float) trdx(nstat),trdy(nstat),trdz(nstat)
    real (kind=c_float) tx_(nstat),ty_(nstat),rx_(nstat),ry_(nstat)
    real (kind=c_float) xmodl(nchnl*ncmp*nstat)
    real (kind=c_float) pbres
    integer (kind=c_int) ntrn
    real (kind=c_float) txarea
    REAL (kind=c_float) a((nchnl*ncmp)*(nlyr*2-1+12*nplt))
    integer (kind=c_int) leroiair_failure_count

    INTEGER NTYRP
    integer (kind=c_int) ncmp,cmp
    INTEGER NDATA
    INTEGER I,NTYPLS
    INTEGER,PARAMETER:: TDFD=1,NPULS=5,GSTRP=0,ASTRP=0,KNRM=3
     INTEGER ISW,IDER,STEP
    REAL :: PULSE
    REAL :: BFFAC
    INTEGER :: NFRQ
    REAL,ALLOCATABLE :: FREQ(:)
    REAL,PARAMETER :: NORM(3) =(/1.,1.,1. /)

    REAL RMU(NLYR),REPS(NLYR),CALF(NLYR),CTAU(NLYR),CFREQ(NLYR)
    REAL CALFP(NPLT),CTAUP(NPLT),CFREQP(NPLT)

    REAL (kind=c_float) CELLW
    INTEGER MXAB,MXB,MCHNL
    LOGICAL SAME_TX(NSTAT)
    REAL :: PRM_TD(3)


    REAL,ALLOCATABLE :: TRP(:)

    LOGICAL,PARAMETER :: TXA90=.FALSE.
    LOGICAL,PARAMETER :: JCBN=.TRUE.

    integer (kind=c_int)  ijac(15)
    REAL, PARAMETER :: PI=3.141592654
    REAL NSHFT,ESHFT


   if (ampt==0) then
       STEP=0
       ISW = 1
       IDER = 0

         WAVEFORM = waveform_ * NTRN * TXAREA
        BFFAC= 1.00000000
    endif

    PULSE=SWX(NSX)

    CALL SET_TRP2(NPULS,PULSE,nchnl,TOPN,nsx,SWX,NTYPLS,NTYRP,TRP)
    CALL DCPRM_TD (trdx(1),trdy(1),trdz(1),tincl(1)/PI*180,txarea,PRM_TD)

    !call PRINT_DCPRM_TD_VARS(trdx(1),trdy(1),trdz(1),tincl(1)/PI*180,txarea,PRM_TD)
    ! HSBOSS_TD and LEROI_3D use different geometries!!!!!!
    CALL GET_ORIGIN_SHIFT(pnorth_,peast_,NPLT,NSHFT,ESHFT)

    tx=ty_ -nshft
    ty=tx_-eshft

    rx(:,1)=ry_-nshft
    ry(:,1)=rx_-eshft

    peast=pnorth_-nshft
    pnorth=peast_-eshft

    do i=1,nplt
       pdzm(i)=pdzm_(i)-1.0/2.0*ACOS(-1.0)
    end do

    !NRXST is equal NSTAT for Td systems.

    NDATA=nchnl*ncmp*nstat
    if (cmp==2) then
       mchnl=nchnl*2
    else if (cmp==3) then
       mchnl=nchnl*3
    else
       mchnl=nchnl
    end if

    RMU=1.
    REPS=1.
    CALF=1.
    CTAU=0.
    CFREQ=1.

    CALFP=1.0
    CTAUP=0.0
    CFREQP=1.0

    call  get_mxb_and_mxab(nplt,plngth1+plngth2,pwdth1+pwdth2,cellw,mxb,mxab)

    !NPAR=nlyr*2-1+nplt*9
    SAME_TX=.TRUE.
    SAME_TX(1)=.FALSE.


    CALL SET_FRQ2(NFRQ,FREQ,TOPN,SWX,NSX)
    CALL SET_SOURCE (STEP,ISW,BFFAC,WAVEFORM,NSX,SWX,SWY,PRM_TD)
 ! CALL PRINT_SET_SOURCE_VARS (STEP,ISW,BFFAC,WAVEFORM,NSX,SWX,SWY,PRM_TD)

    call FORJAC2 (nlyr,res,pbres,thk, &
         NPLT,pres,plngth1,plngth2,pwdth1,pwdth2,pthk,peast,pnorth,ptop,pdzm,pdip,plng, &
         NDATA,XMODL,JCBN,A, &
         TDFD, CMP,STEP,IDER,NSX,SWX,SWY,NTYRP,TRP,NPULS,PULSE,NTYPLS, &
         NCMP,NCHNL,MCHNL,TOPN, TCLS,GSTRP,ASTRP,NFRQ,FREQ,KNRM,NORM,TXA90, &
         NSTAT,trdx,trdy,trdz,tx,ty,tz,tincl,tazi,SAME_TX,RX,RY,RZ, &
         RMU,REPS,CALF,CTAU,CFREQ,CALFP,CTAUP,CFREQP, &
         CELLW,MXAB,MXB,ijac,leroiair_failure_count)

    DEALLOCATE(FREQ)

  end subroutine formod_helitem_data_and_jacobian

  subroutine formod_helitem_data(nlyr,res,pbres,thk,  &
       nplt,peast_,pnorth_,ptop,pres,plngth1,plngth2,pwdth1,pwdth2,pthk,pdzm_,pdip,plng,  &
       ntrn,txarea,nchnl,topn,tcls,nsx,swx,ampt,waveform_,nstat,tx_,ty_,tz,tazi,tincl,rx_,ry_,rz,trdx,trdy,trdz, &
       xmodl,ncmp,cmp,cellw,leroiair_failure_count) bind(c,name="leroiair_formod_helitem_data")

    ! Geotem system forward modelling
    !              nlyr - number of layers
    !         res(nlyr) - resistivity values
    !       thk(nlyr-1) - layer thickness
    !             nchnl - number of time domain channels
    !       topn(nchnl) - time at which receiver channel I opens.
    !       tcls(nchnl) - time at which receiver channel I closes.
    !                tz - transmitter elevation
    !            txincl - transmitter inclination
    !              trdx - in-line receiver offset for each frequency;  behind = positive
    !              trdy - transverse receiver offset for each frequency; Port = positive.
    !              trdz - vertical receiver offset for each frequency; below = positive
    !    xmodl(nchnl*2) - model response vertical component and then inline
    ! a(nchnl*2nlyr*2-1) - Jacobian
    use iso_c_binding
        use leroiair_subroutines

    implicit none
    integer (kind=c_int) nlyr
    integer (kind=c_int) ampt
    real (kind=c_float) res(nlyr*nstat),thk((nlyr-1)*nstat)
    integer (kind=c_int) nplt
    real (kind=c_float) peast(nplt),pnorth(nplt),ptop(nplt),pres(nplt)
    real(kind=c_float) plngth1(nplt),plngth2(nplt),pwdth1(nplt),pwdth2(nplt),pthk(nplt)
    real (kind=c_float) pdzm(nplt),pdip(nplt),plng(nplt)
    real (kind=c_float) peast_(nplt),pnorth_(nplt),pdzm_(nplt)
    integer (kind=c_int)  nchnl
    real (kind=c_float) topn(nchnl),tcls(nchnl)
    integer (kind=c_int) nsx
    real (kind=c_float) swx(nsx),waveform(nsx),waveform_(nsx)
    real (kind=c_float) swy(nsx,3)
    integer (kind=c_int) nstat
    real (kind=c_float) tx(nstat),ty(nstat),tz(nstat),tazi(nstat),tincl(nstat)
    real (kind=c_float) rx(nstat),ry(nstat),rz(nstat),trdx(nstat),trdy(nstat),trdz(nstat)
    real (kind=c_float) tx_(nstat),ty_(nstat),rx_(nstat),ry_(nstat)
    real (kind=c_float) xmodl(nchnl*ncmp*nstat)
    real (kind=c_float) pbres

    integer (kind=c_int) ntrn
    real (kind=c_float) txarea

    REAL (kind=c_float) a((nchnl*ncmp)*(nlyr*2-1+12*nplt))
    integer (kind=c_int) leroiair_failure_count

    INTEGER NTYRP
    integer (kind=c_int)  ncmp,cmp
    INTEGER NDATA
    INTEGER I,NTYPLS
    INTEGER,PARAMETER:: TDFD=1,NPULS=5,GSTRP=0,ASTRP=0,KNRM=3
        INTEGER ISW,IDER,STEP
    REAL :: PULSE
    REAL :: BFFAC
    INTEGER :: NFRQ
    REAL,ALLOCATABLE :: FREQ(:)
    REAL,PARAMETER :: NORM(3) =(/1.,1.,1. /)

    REAL RMU(NLYR),REPS(NLYR),CALF(NLYR),CTAU(NLYR),CFREQ(NLYR)
    REAL CALFP(NPLT),CTAUP(NPLT),CFREQP(NPLT)
    integer (kind=c_int)  ijac(15)
    REAL (kind=c_float)  CELLW
    INTEGER MXAB,MXB,MCHNL
    LOGICAL SAME_TX(NSTAT)
    REAL :: PRM_TD(3)

    REAL,ALLOCATABLE :: TRP(:)
    LOGICAL,PARAMETER :: TXA90=.FALSE.
    LOGICAL,PARAMETER :: JCBN=.FALSE.


    REAL, PARAMETER :: PI=3.141592654
    REAL NSHFT,ESHFT

	ijac=0

   if (ampt==0) then
       STEP=0
       ISW = 1
       IDER = 0
       WAVEFORM = waveform_ * NTRN * TXAREA
       BFFAC= 1.00000000
    endif

    PULSE=SWX(NSX)

    CALL SET_TRP2(NPULS,PULSE,nchnl,TOPN,nsx,SWX,NTYPLS,NTYRP,TRP)
    CALL DCPRM_TD (trdx(1),trdy(1),trdz(1),tincl(1)/PI*180,txarea,PRM_TD)

    ! HSBOSS_TD and LEROI_3D use different geometries!!!!!!
    CALL GET_ORIGIN_SHIFT(pnorth_,peast_,NPLT,NSHFT,ESHFT)

    tx=ty_-nshft
    ty=tx_-eshft

    rx=ry_-nshft
    ry=rx_-eshft

    peast=pnorth_-nshft
    pnorth=peast_-eshft

    do i=1,nplt
       pdzm(i)=pdzm_(i)-1.0/2.0*ACOS(-1.0)
    end do

    !NRXST is equal NSTAT for Td systems.

    NDATA=nchnl*ncmp*nstat
    if (cmp==2) then
       mchnl=nchnl*2
    else if (cmp==3) then
       mchnl=nchnl*3
    else
       mchnl=nchnl
    end if

    RMU=1.
    REPS=1.
    CALF=1.
    CTAU=0.
    CFREQ=1.

    CALFP=1.0
    CTAUP=0.0
    CFREQP=1.0

    call  get_mxb_and_mxab(nplt,plngth1+plngth2,pwdth1+pwdth2,cellw,mxb,mxab)

    SAME_TX=.TRUE.
    SAME_TX(1)=.FALSE.

    CALL SET_FRQ2(NFRQ,FREQ,TOPN,SWX,NSX)
    CALL SET_SOURCE (STEP,ISW,BFFAC,WAVEFORM,NSX,SWX,SWY,PRM_TD)

    call FORJAC2 (nlyr,res,pbres,thk, &
         NPLT,pres,plngth1,plngth2,pwdth1,pwdth2,pthk,peast,pnorth,ptop,pdzm,pdip,plng, &
         NDATA,XMODL,JCBN,A, &
         TDFD, CMP,STEP,IDER,NSX,SWX,SWY,NTYRP,TRP,NPULS,PULSE,NTYPLS, &
         NCMP,NCHNL,MCHNL,TOPN, TCLS,GSTRP,ASTRP,NFRQ,FREQ,KNRM,NORM,TXA90, &
         NSTAT,trdx,trdy,trdz,tx,ty,tz,tincl,tazi,SAME_TX,RX,RY,RZ, &
         RMU,REPS,CALF,CTAU,CFREQ,CALFP,CTAUP,CFREQP, &
         CELLW,MXAB,MXB,ijac,leroiair_failure_count)

    DEALLOCATE(FREQ)
  end subroutine formod_helitem_data


  SUBROUTINE SET_FRQ2(NFRQ,FREQ,TOPN,SWX,NSX)
    !  ------------------

    !***  Called by: MAIN

    !  For time-domain options:
    !
    !  Numerical experiments where the frequency-domain responses of layered
    !  layered 1/2 space models were transformed to time-domain for a wide
    !  range of resistivities has indicated that 6 points per decade is
    !  adequate frequency discretisation, even at the high end.  Moreover,
    !  for perfect frequency-domain data, only 3 frequencies are required
    !  from 1 to 10 Hz to maintain an accuracy of better than .1 percent for
    !  the models studied.  This means that 28 frequencies are needed for the
    !  1Hz to 100kHz range.
    !
    !  The need to go to 1 MHz depends upon receiver channels and conductivity.
    !  A 10,000 ohm-m 1/2 space requires extended range if the earliest window
    !  opening is within .28 ms of signal turn-off.  A 1 ohm-m 1/2 space allows
    !  normal range if the first window open time is > .002 ms.
    !
    !  In a 3D program, this is pretty hard to control because of non-uniform
    !  resistivity distribution so the safe option of going to 1 MHz is chosen
    !  if the first window opens earlier than .28 ms after signal turn-off.  This
    !  requires computation for 34 frequencies.
    !
    !  This can be over-ridden by setting TDFD = 0 in which case the user needs to
    !  specify five integers: KLO, KHI, KMD, PPD1, PPD2
    !
    !    KLO - Lowest frequency expressed as an integer power of 10
    !          For example setting LOW = -1, 0 or 1 would mean that the lowest
    !          frequency would be either 0.1, 1.0, or 10.0 Hz repectively.
    !
    !    KHI - Highest frequency expressed as an integer power of 10
    !          For example setting HIGH = 5 or 6 would mean that the highest
    !          frequency would be either 0.1 or 1.0 MHz repectively.
    !
    !    KMD - Frequencies below MID are spaced at PPD1 points per decade
    !          Frequencies above MID are spaced at PPD2 points per decade
    !
    !    The allowed values for both PPD1 and PPD2 are 3, 6 or 12.
    !    If PPD1 < 3, it is changed to 3
    !    If 3 < PPD1 < 6, it is changed to 6
    !    If PPD1 > 6, it is set to 12
    !
    !  Thus, depending upon the first channel time, for the default case of TDFD = 1,
    !
    !  KLO = 0,  KHI = 5,  KMD = 10,  PPD1 = 3 and PPD2 = 6 => NFRQ = 28   or
    !  KLO = 0,  KHI = 6,  KMD = 10,  PPD1 = 3 and PPD2 = 6 => NFRQ = 34
    use leroiair_subroutines
    IMPLICIT NONE
    INTEGER J,NFRQ,PPD1,PPD2,NSX, KLO,KHI,KMD
    INTEGER, PARAMETER :: QL=SELECTED_REAL_KIND(12,80)
    REAL T0
    REAL,ALLOCATABLE :: FREQ(:)
    REAL TOPN(:),SWX(:)
    REAL RLO,RHI,RMD
    REAL(KIND=QL) QLO,QHI,QMD,QFRQ1,QFRQ2,FQQ
    REAL(KIND=QL), ALLOCATABLE :: FDUM(:)

    ALLOCATE (FDUM(1000))

    KLO = 0
    KHI = 5
    KMD = 1
    PPD1 = 3
    PPD2 = 6
    T0 = MINVAL (TOPN) - SWX(NSX)
    IF (T0 < .28E-3) KHI = 6         ! Extend range to 1 MHz


    RLO = 10.**KLO
    RHI = 10.**KHI
    RMD = 10.**KMD
    QLO = REAL (RLO,QL)
    QHI = REAL (RHI,QL)
    QMD = REAL (RMD,QL)
    QFRQ1 = EXP ( LOG (10.D0) / REAL (PPD1,QL) )
    QFRQ2 = EXP ( LOG (10.D0) / REAL (PPD2,QL) )




    FDUM(1) = QLO
    NFRQ = 1
    FQQ = REAL (FDUM(1),KIND=QL)

    DO J = 2,1000
       IF (FDUM(J-1) * QFRQ1 < QMD) THEN
          FQQ = FDUM(J-1) * QFRQ1
       ELSE
          FQQ = FDUM(J-1) * QFRQ2
       END IF
       IF (FQQ > 1.001 * QHI) EXIT
       NFRQ = J
       FDUM(J) = FQQ
    END DO

    ALLOCATE (FREQ(NFRQ))
    FREQ(1:NFRQ) = REAL (FDUM(1:NFRQ))
    DEALLOCATE (FDUM)



  END SUBROUTINE SET_FRQ2


  SUBROUTINE SET_TRP2(NPULS,PULSE,nchnl,TOPN,nsx,SWX,NTYPLS,NTYRP,TRP)
    !  ------------------

    !***  Called by: MAIN

    !  Sets up interpolation times for FD -> TD transform which use the
    !  exact 6 points per decade frequency-domain data plus 6 per decade
    !  interpolated values.  These are based on a 12 point per decade
    !  cosine filter derived from the Niels Christensen routine FILCOA
    !  with OMEGA = .3 PI and shift 0.

    !             OUTPUT
    !             ------

    !        TRP - array of time values for FD -> TD transformations
    !      NTYRP - number of values in TRP
    !     EXTENT - the latest time for which time-domain output is required.
    !      PULSE - time length of one signal pulse
    !     NTYPLS - number of TRP values in 1 PULSE


    use leroiair_subroutines

    IMPLICIT NONE
    INTEGER NPULS,NTYPLS
    REAL PULSE
    integer nsx,nchnl
    real swx(nsx)
    real topn(nchnl)
    integer NTYRP
    REAL,ALLOCATABLE :: TRP(:)
    INTEGER MXTYM,J1
    REAL T0,EXTENT
    REAL,ALLOCATABLE :: QQQ(:)
    INTEGER, PARAMETER :: QL=SELECTED_REAL_KIND(12,80)
    REAL(KIND=QL) TBASE,QTYM, TQ

    REAL, PARAMETER :: TWOPI=6.2831853
    REAL, PARAMETER :: T0_MIN=1.E-7

    MXTYM=200
    ALLOCATE (QQQ(MXTYM))
    QQQ = 0.

    QTYM = LOG (10.D0) /12.D0
    QTYM = EXP (QTYM)
    EXTENT = 2.0 * NPULS * PULSE


    T0 = MINVAL (TOPN) - SWX(NSX)
    T0 = MAX (T0, T0_MIN)
    TBASE = 1.D0 / DBLE (TWOPI)
    DO J1 = 1,MXTYM
       IF (TBASE < T0) EXIT
       TBASE = TBASE / QTYM
    END DO

    TQ = TBASE
    QQQ(1) = REAL (TQ)
    DO J1 = 2, MXTYM
       NTYRP = J1
       TQ = TQ * QTYM
       QQQ(J1) = REAL(TQ)
       IF (QQQ(J1) < PULSE) NTYPLS = J1+2
       IF( QQQ(J1) > EXTENT) EXIT
    END DO

    ALLOCATE (TRP(NTYRP))
    TRP(1:NTYRP) = QQQ(1:NTYRP)
    DEALLOCATE (QQQ)

  END SUBROUTINE SET_TRP2




  SUBROUTINE SET_CELLS_adv (KP,NPLT,NLYR,MXAB,MXB,CELLW,NSTAT,THK_NSTAT,PLNGTH,PLWDTH,NA,NB,DA,DB, &
       XCNTR,YCNTR,PLTOP,PLAZM,PLDIP,XCELL,YCELL,ZCELL,PCNR)
    !----------------------------_-----------------------------------------------------

    !*** Called by: MAIN, FORJAC

    !   KP = 0 converts the location and dimension of all plates into cells.
    !          all plates extending into OB are pushed down into host
    !
    !   1 <= KP <= NPLT  converts the location and dimension of plate KP into cells.
    !                    This would be used during the Jacobian process when plates
    !                    are varied individually.  It wouldn't go into the OB because
    !                    location can only increase downwards in this case.
    !
    !
    !             OTHER INPUT
    !             -----------
    !          KP = 0 => push plate down into host and leave OB as is.
    !             = 1 => Leave the plate where it is but decrease OB thickness &
    !          KP = (1, 2, or 3) => shift cells (north, east, down) by DELQ
    !        NPLT = number of thin plates
    !        NLYR = number of layers
    !         THK = layer thicknesses
    !        MXAB - maximum number of cells per plate
    !         MXB - maximum number of cell rows down dip per plate
    !       CELLW = cell discretisation dimension
    !      PLNGTH = plate strike lengths
    !      PLWDTH = plate widths along dip
    !       XCNTR = north coordinates of plate reference points
    !       YCNTR = east coordinates of plate reference points
    !       PLAZM = plate azimuths in radians (0 = north)
    !       PLDIP = plate dip angle in radians
    !
    !               INPUT/OUTPUT
    !               -----------
    !           THK = layer thicknesses
    !        NA, NB = Number of cells along strike and down dip respectively
    !        DA, DB = Length of cells along strike and down dip respectively
    !         PLTOP = depth (+) from ground surface to reference point on top edge of plate
    !   XCELL(k,JP) - X (north) coordinate of centre of cell k of plate JP.
    !   YCELL(k,JP) - Y (east) coordinate of centre of cell k of plate JP.
    !   ZCELL(i,JP) - depth of cell centre in row i of plate nrelative to surface.
    !  PCNR(I,J,JP) = Ith component (x,y,z) of Jth corner (1 to 4) of plate JP

    use leroiair_subroutines
    IMPLICIT NONE
    INTEGER NPLT,NLYR,MXAB,MXB,MXA,NA(NPLT),NB(NPLT),KP,JP,JA,JB,JAB,KPF,KPS,NSTAT
    REAL THICK,THK_NSTAT((NLYR-1)*NSTAT),CELLW,TMPL,SDIP,CDIP,CSTR,SSTR,PL2,X0,Y0,Z0,XCTMP,YCTMP,HSTR,PLWC, &
         XCELL(MXAB,NPLT),YCELL(MXAB,NPLT),ZCELL(MXB,NPLT),PCNR(3,4,NPLT)
    REAL, DIMENSION(NPLT) :: XCNTR,YCNTR,PLTOP,PLNGTH,PLWDTH,PLAZM,PLDIP,DA,DB


    IF (KP > NPLT) KP = 0
    IF (KP > 0) THEN
       KPS = KP
       KPF = KP
    ELSE
       KPS = 1
       KPF = NPLT
    END IF

    !  Plates must be at least 1 metre below surface.  If necessary, adjust
    !  PLTOP so that all plates are in the basement.

    THICK = SUM (THK_NSTAT)/FLOAT(NSTAT)

    PLTOP = MAX (PLTOP,1.0)
    PLTOP = MAX (PLTOP, THICK)

    DO JP = KPS, KPF
       CDIP = COS (PLDIP(JP))
       SDIP = SIN (PLDIP(JP))
       IF (ABS (CDIP) < .02) CDIP = 0.
       IF (ABS (SDIP) < .02) SDIP = 0.

       CSTR = COS (PLAZM(JP))
       SSTR = SIN (PLAZM(JP))
       IF (ABS (CSTR) < 1.0E-3) CSTR = 0.
       IF (ABS (SSTR) < 1.0E-3) SSTR = 0.

       TMPL = MIN (PLNGTH(JP), PLWDTH(JP)) / 2.
       TMPL = MIN (TMPL, CELLW) + .01
       NB(JP) = CEILING (PLWDTH(JP) / TMPL)
       NB(JP) = MAX (2,NB(JP))
       NB(JP) = MIN (MXB, NB(JP))

       MXA = MXAB / NB(JP)
       NA(JP) = CEILING (PLNGTH(JP) / TMPL)
       NA(JP) = MAX (2,NA(JP))
       NA(JP) = MIN (MXA, NA(JP))

       DA(JP) = PLNGTH(JP) / REAL (NA(JP),4)
       DB(JP) = PLWDTH(JP) / REAL (NB(JP),4)
       PL2 = (NA(JP) + 1) * DA(JP) / 2.
       X0 = - PL2
       Y0 = -CDIP * DB(JP) /2.
       Z0 = PLTOP(JP) - SDIP * DB(JP) /2.

       DO JB = 1, NB(JP)
          ZCELL(JB,JP) = Z0 + JB* DB(JP) *SDIP
          DO JA = 1, NA(JP)
             JAB = JA + ((JB-1) * NA(JP))
             XCTMP = X0 + JA *DA(JP)
             YCTMP = Y0 + JB * DB(JP) * CDIP
             XCELL(JAB,JP) = XCNTR(JP) + (XCTMP * CSTR) - (YCTMP * SSTR)
             YCELL(JAB,JP) = YCNTR(JP) + (YCTMP * CSTR) + (XCTMP * SSTR)
          END DO
       END DO
       HSTR = PLNGTH(JP) /2.
       PLWC = PLWDTH(JP) * CDIP

       PCNR(3,1:2,JP) = PLTOP(JP)                        ! Z for top corners
       PCNR(3,3:4,JP) = PLTOP(JP) + PLWDTH(JP) * SDIP    ! Z for bottom corners

       PCNR(1,1,JP) = XCNTR(JP) - HSTR * CSTR     ! X,Y for top corners
       PCNR(2,1,JP) = YCNTR(JP) - HSTR * SSTR
       PCNR(1,2,JP) = XCNTR(JP) + HSTR * CSTR
       PCNR(2,2,JP) = YCNTR(JP) + HSTR * SSTR

       PCNR(1,3,JP) = PCNR(1,1,JP) - PLWC * SSTR  ! X,Y for bottom corners
       PCNR(2,3,JP) = PCNR(2,1,JP) + PLWC * CSTR
       PCNR(1,4,JP) = PCNR(1,2,JP) - PLWC * SSTR
       PCNR(2,4,JP) = PCNR(2,2,JP) + PLWC * CSTR





    END DO
    ! print * ,'east','       north'
    !print *,pcnr(1,1,1),pcnr(2,1,1)
    !print *,pcnr(1,2,1),pcnr(2,2,1)
    !print *,pcnr(1,3,1),pcnr(2,3,1)
    !print *,pcnr(1,4,1),pcnr(2,4,1)

    !print *,'north'
    !print *,xcell
    !
    !print *,'east'
    !print *,ycell




  END SUBROUTINE SET_CELLS_adv



  SUBROUTINE FORJAC2 (NLYR,RES,PBRES,THK, &
       NPLT,PRES,PLNGTH1,PLNGTH2,PWDTH1,PWDTH2,PTHK,XP,YP,ZP,PLAZM,PLDIP,PLUNJ, &
       NDATA,XMODL,JCBN,A, &
       TDFD, CMP,STEP,IDER,NSX,SWX,SWY,NTYRP,TRP,NPULS,PULSE,NTYPLS, &
       NCMP,NCHNL,MCHNL,TOPN, TCLS,GSTRP,ASTRP,NFRQ,FREQ,KNRM,NORM,TXA90, &
       NSTAT,XRX,YRX,ZRX,SX,SY,SZ,TXCLN,FANGLE,SAME_TX,RX,RY,RZ, &
       RMU,REPS,CALF,CTAU,CFREQ,CALFP,CTAUP,CFREQP, &
       CELLW,MXAB,MXB,IJAC,leroiair_failure_count)
    !----------------------------------------------------------------------------------------------

    !  Sets up and calls for model computation.
    !  It also calculates the Jacobian and error vector if required
    !  New convention: Dec, 2003: VERR is now VD - VM
    !  Thus DELPAR is now added rather than subtracted during updates.

    !*** Called by: NLSQ2
    !***     Calls: SET_CELLS, SHIFT_CELLS, GET_FWD_MODL, CNVRT2_MPAR

    !             General Inversion Input Variables
    !             ---------------------------------
    !
    !       NDATA - dimension of vector to be inverted: = NCMP * NCHNL or 2*NFRQ.
    !       DNORM(NDATA) : L1 norm of data for each channel or frequency, averaged over the survey
    !        INRM = 1 : point symmetric norm used for fitting error and sensitivity matrix
    !        INRM = 2 : survey norm used for fitting error and sensitivity matrix
    !       XDATA - data to be inverted in user-specified units
    !       XMODL - model data in user-specified units
    !        XWTS - weights for XDATA  (0 or 1)
    !        NPAR - number of parameters to be inverted (nominally 2*NLYR - 1)
    !        XPAR - array of transformed model parameters
    !       CXPAR = 0 => parameter is completely free to vary as dictated by inversion step
    !             = 1 => parameter is fixed
    !             = 2 => parameter is constrained by elasticity.
    !             = 3 => parameter bounds are buffered.
    !
    !             General Inversion Output Variables
    !             ----------------------------------
    !
    !       SUMSQ(1) - sum squared scaled error for point norm
    !       SUMSQ(2) - sum squared scaled error for survey norm
    !       JCBN    - true if Jacobian required
    !       A       - a large array which carries the Jacobian out of this routine
    !       VERR(J) - scaled error in channel J
    !
    !             AEM System Input Variables
    !             --------------------------
    !
    !        TDFD = 1 for TD; 2 for FD
    !         CMP - (time-domain only) component to be inverted
    !               11: in-line; 13: vertical; 2: joint vertical & in-line
    !               3 component; 4 total field
    !
    !        KNRM - dimension of NORM: = 3 for time-domain,  = NFRQ for frequency-domain
    !        NORM - PPM conversion
    !        FREQ - array of NFRQ frequencies
    !        STEP = 1 iff step response is to be computed
    !        IDER = 1 if source waveform was dB/dt; = 0 if amps pr B
    !         NSX - number of points used to discretise transmitter signal
    !         SWX - abscissae (seconds) of current waveform
    !         SWY - dI/dt * Tx moment & nanotesla conversion at times SWX
    !       NTYRP - number of values in TRP for total signal length: 2 * NPULS *PULSE
    !         TRP - array of time values for FD -> TD transformations
    !      NTYPLS - number of TRP values in 1 PULSE
    !       NCHNL = number of time domain channels
    !       MCHNL = total number of readings per station to be inverted(TD or FD)
    !             = NCHNL for time-domain when CMP = 11, 13, 4, 42, 43
    !             = 2* NCHNL for time-domain when CMP = 2
    !             = 3* NCHNL for time-domain when CMP = 3
    !
    !             = 2 * NFRQ for frequency-domain
    !        TOPN - time at which receiver channel I opens.
    !        TCLS - time at which receiver channel I closes.
    !       GSTRP = 1 => apply Questem-Geotem stripping algorithm
    !       ASTRP = 1 => apply Aerotem stripping algorithm
    !        FREQ = array of NFRQ frequencies
    !       TXCLN - angle in radians that TX dipole makes with vertical (climb = +)
    !       TXA90 - true for vertical co-planar briadside array
    !         NRX - number of receivers = NFRQ in FD; = 1 in TD
    !       NRXST - dimension for receiver offset & transmitter tilt
    !             = NFRQ in FD;  = NSTAT in TD
    !         ZRX - vertical receiver offset for each frequency;   below = positive
    !         XRX - in-line receiver offset for each frequency;    behind = positive
    !         YRX - transverse receiver offset for each frequency; left = positive.
    !
    !             AEM Survey Input Variables
    !             --------------------------
    !
    !       NSTAT - number of stations in survey line.
    !      FANGLE - flight path angle in radians. (North = 0; East = PI/2)
    !     SAME_TX - used to reduce redundant computations
    !       SX, SY, SZ: north, east and altitude (re gnd level) of transmitter
    !       RX, RY, RZ: north, east and altitude (re gnd level) of receiver(s)
    !
    !             Model Description Input Variables
    !             ---------------------------------
    !
    !        NLYR - number of layers (1 or 2)
    !         THK - layer thicknesses
    !         RMU - mu(i) / mu(0)
    !        REPS - relative dielectric constant
    !        CALF, CTAU & CFREQ are the layered earth Cole-Cole parameters.
    !
    !      XYNORM - survey dependent normalisation for plate X-Y location.
    !        NPLT - Number of plates
    !       CELLW - maximum cell dimension
    !        MXAB - Number of cells in biggest plate
    !         MXB - Maximim number of cell rows
    !        CALFP, CTAUP & CFREQP are the plate Cole-Cole parameters.
    !
    !
    !          PHYSICAL PARAMETER 0RDERING
    !          ---------------------------
    !
    !     KP = 9* (JP-1) where JP = plate index
    !
    !     KP + 1 - SIG_T = conductance of plate JP
    !        + 2 - PLTOP = depth of reference point of plate JP
    !        + 3 - PLNGTH = strike length of plate JP
    !        + 4 - PLWDTH = plate width of plate JP
    !        + 5 - YCNTR = east coordinate of reference point of plate JP
    !        + 6 - XCNTR = north coordinate of reference point of plate JP
    !        + 7 - PLAZM = strike angle of plate JP
    !        + 8 - PLDIP = dip angle of plate JP
    !        + 9 - PLUNJ = plunge angle of plate JP
    !
    !      9*NPLT + 1 - host resistivity
    !      9*NPLT + 2 - overburden resistivity (only if NLYR = 2)
    !      9*NPLT + 3 - overburden thickness   (only if NLYR = 2)
    !
    !
    !          TRANSFORMED PARAMETERS FOR INVERSION
    !          ------------------------------------
    !
    !   SIG_T, PLTOP, PLWDTH, PLNGTH, RES & THK are represented logarithmically
    !
    !   XCNTR * YCNTR are normalised to PLTOP
    !
    !   PLAZM, PLDIP & PLUNJ are normalised to PI

          use leroiair_subroutines
      use franken_subroutines
    !   use print_leroiair_vars

    IMPLICIT NONE
    REAL, PARAMETER :: PI=3.141592654
    INTEGER NDATA,NPAR,NPLT,NLYR,XWTS(NDATA),NCHNL,MCHNL,TDFD,CMP,INRM,KNRM, &
         NFRQ,MXAB,MXB,NA(NPLT),NB(NPLT),NSTAT,NRX,NRXST,STEP,IDER,NSX,NPULS,NTYPLS,  &
         NTYRP,GSTRP,ASTRP,JS,JD,JP,JL,LP,JF,JT,KS1,KS2,KP,K0
    REAL NORM(KNRM),FREQ(NFRQ),CELLW,SWX(NSX),SWY(NSX,3), &
         PULSE,TOPN(NCHNL),TCLS(NCHNL),TRP(NTYRP),X2,X3,VM,VJ,DENOM,DELTA,    &
         DELXY,DELZ,DELPHI,PCNR(3,4,NPLT),XP0,PARFAC,DNORM(NDATA)

    INTEGER NCMP,IJAC(15)
    REAL a(nstat*((nlyr*2-1)*(nchnl*ncmp)+(nchnl*ncmp)*12*nplt))

    REAL pbres
    REAL, DIMENSION(NDATA) ::  XMODL,XMODL0
    REAL, DIMENSION(NSTAT) :: TXCLN,XRX,YRX,ZRX
    REAL, DIMENSION(NSTAT) :: SX,SY,SZ,FANGLE
    REAL, DIMENSION(NSTAT,1) :: RX,RY,RZ
    REAL, DIMENSION(NLYR) :: REPS,RMU,CALF,CTAU,CFREQ
    REAL :: RES(NLYR*NSTAT),THK((NLYR-1)*NSTAT)
    REAL, DIMENSION(NPLT) :: PRES,PTHK,CALFP,CTAUP,CFREQP,PLTOP,PLWDTH,XCNTR,YCNTR,PLNGTH, &
         PLAZM,PLDIP,PLUNJ,DA,DB,PLNGTH1,PLNGTH2,PWDTH1,PWDTH2,&
         XP,YP,ZP
    REAL, DIMENSION(MXAB,NPLT) :: XCELL,YCELL,ZCELL
    REAL, DIMENSION(NCHNL,NSTAT,3) :: BTD,BTD_SCAT
    COMPLEX XBFD
    COMPLEX, DIMENSION(NFRQ,NSTAT,3) :: BFD,BFD_SCAT
    LOGICAL JCBN,TXA90,SAME_TX(NSTAT)
    integer I
    REAL ZP0(NSTAT)
    ! Compute initial model & compute error
    INTEGER leroiair_failure_count
    REAL PPARFAC,PARFAC0
    leroiair_failure_count=0






    NRXST=NSTAT
    NRX=1
    INRM=1
    XWTS=1
    NPAR=(nlyr*2-1+12*nplt)

    KP = 0
    PLNGTH=PLNGTH1+PLNGTH2
    PLWDTH=PWDTH1+PWDTH2
    CALL GET_P223_PLATE_REFERENCE_POINT(NPLT, XP,YP,ZP,PLAZM,PLDIP,PLUNJ,plngth1,plngth2,PWDTH1,XCNTR,YCNTR,PLTOP)
    CALL SET_CELLS (KP,NPLT,NLYR,MXAB,MXB,CELLW,THK,PLNGTH,PLWDTH,NA,NB,DA,DB, &
         XCNTR,YCNTR,PLTOP,PLAZM,PLDIP,XCELL,YCELL,ZCELL,PCNR)





    LP = 0;
    CALL GET_FWD_MODL(LP,PARFAC)
    XMODL0 = XMODL

    !XYNORM = (SUM (SZ) + SUM(RZ) )/ REAL (2*NSTAT)
    !XYNORM = sum(sqrt(xrx*xrx+yrx*yrx))/REAL (NSTAT)
    !  Initialise and then compute the Jacobian as the derivative of log(volts) wrt
    !  log(parameter) for a three percent step.  Skip over held parameters

    IF (JCBN) THEN
       A = 0.
       PARFAC0=0.03
       DELXY=0.025
       DELZ=-0.05
       DELPHI =PI*PARFAC0*0.25
       DO JP = 1,NPLT
          KP = JP
          K0 = 12*(JP-1)
          LP = K0 + 1

          if (ijac(1)==1) then
             XP0 = YP(JP)
             YP(JP) = YP(JP) + DELXY
             PLNGTH=PLNGTH1+PLNGTH2
             PLWDTH=PWDTH1+PWDTH2
             CALL GET_P223_PLATE_REFERENCE_POINT(NPLT, XP,YP,ZP,PLAZM,PLDIP,PLUNJ,plngth1,plngth2,PWDTH1,XCNTR,YCNTR,PLTOP)
             CALL SET_CELLS (KP,NPLT,NLYR,MXAB,MXB,CELLW,THK,PLNGTH,PLWDTH,NA,NB,DA,DB, &
                  XCNTR,YCNTR,PLTOP,PLAZM,PLDIP,XCELL,YCELL,ZCELL,PCNR)
             PPARFAC=DELXY
             CALL GET_FWD_MODL(LP,PPARFAC)
             if (leroiair_failure_count>0) return
             YP(JP) = XP0
          end if

          LP = K0 + 2
          if (ijac(2)==1) then
             XP0 = XP(JP)
             XP(JP) = XP(JP) + DELXY
             PLNGTH=PLNGTH1+PLNGTH2
             PLWDTH=PWDTH1+PWDTH2
             CALL GET_P223_PLATE_REFERENCE_POINT(NPLT, XP,YP,ZP,PLAZM,PLDIP,PLUNJ,plngth1,plngth2,PWDTH1,XCNTR,YCNTR,PLTOP)
             CALL SET_CELLS (KP,NPLT,NLYR,MXAB,MXB,CELLW,THK,PLNGTH,PLWDTH,NA,NB,DA,DB, &
                  XCNTR,YCNTR,PLTOP,PLAZM,PLDIP,XCELL,YCELL,ZCELL,PCNR)
             PPARFAC=DELXY
             CALL GET_FWD_MODL(LP,PPARFAC)
             if (leroiair_failure_count>0) return

             XP(JP) = XP0
          end if

          LP = K0 + 3
          if (ijac(3)==1) then
             XP0 = ZP(JP)
             ZP(JP) = ZP(JP) + DELZ
             CALL GET_P223_PLATE_REFERENCE_POINT(NPLT, XP,YP,ZP,PLAZM,PLDIP,PLUNJ,plngth1,plngth2,PWDTH1,XCNTR,YCNTR,PLTOP)
             CALL SHIFT_CELLS (JP,NPLT,MXB,NB,DELZ,ZCELL)
             PPARFAC=-DELZ
             CALL GET_FWD_MODL(LP,PPARFAC)
             ZP(JP) = XP0
          end if
          LP = K0 + 4
          if (ijac(4)==1) then
             XP0 = PRES(JP)
             PARFAC=PARFAC0
             PRES(JP) = (1. + PARFAC) * PRES(JP)
             CALL GET_FWD_MODL(LP,PARFAC)
             if (leroiair_failure_count>0) return
             PRES(JP) = XP0
          end if

          LP = K0 + 5
          if (ijac(5)==1) then
             XP0 = PLNGTH1(JP)
             PARFAC=PARFAC0
             PLNGTH1(JP) = (1. + PARFAC) * PLNGTH1(JP)
             PLNGTH=PLNGTH1+PLNGTH2
             PLWDTH=PWDTH1+PWDTH2
             CALL GET_P223_PLATE_REFERENCE_POINT(NPLT, XP,YP,ZP,PLAZM,PLDIP,PLUNJ,plngth1,plngth2,PWDTH1,XCNTR,YCNTR,PLTOP)
             CALL SET_CELLS (KP,NPLT,NLYR,MXAB,MXB,CELLW,THK,PLNGTH,PLWDTH,NA,NB,DA,DB, &
                  XCNTR,YCNTR,PLTOP,PLAZM,PLDIP,XCELL,YCELL,ZCELL,PCNR)
             CALL GET_FWD_MODL(LP,PARFAC)
             if (leroiair_failure_count>0) return

             PLNGTH1(JP) = XP0
          end if

          LP = K0 + 6
          if (ijac(6)==1) then
             XP0 = PLNGTH2(JP)
             PARFAC=PARFAC0
             PLNGTH2(JP) = (1. + PARFAC) * PLNGTH2(JP)
             PLNGTH=PLNGTH1+PLNGTH2
             PLWDTH=PWDTH1+PWDTH2
             CALL GET_P223_PLATE_REFERENCE_POINT(NPLT, XP,YP,ZP,PLAZM,PLDIP,PLUNJ,plngth1,plngth2,PWDTH1,XCNTR,YCNTR,PLTOP)
             CALL SET_CELLS (KP,NPLT,NLYR,MXAB,MXB,CELLW,THK,PLNGTH,PLWDTH,NA,NB,DA,DB, &
                  XCNTR,YCNTR,PLTOP,PLAZM,PLDIP,XCELL,YCELL,ZCELL,PCNR)
             CALL GET_FWD_MODL(LP,PARFAC)
             if (leroiair_failure_count>0) return

             PLNGTH2(JP) = XP0
          end if


          LP = K0 + 7
          if (ijac(7)==1) then
             XP0 = PWDTH1(JP)
             PARFAC=PARFAC0
             PWDTH1(JP) = (1. + PARFAC) * PWDTH1(JP)
             PLNGTH=PLNGTH1+PLNGTH2
             PLWDTH=PWDTH1+PWDTH2
             CALL GET_P223_PLATE_REFERENCE_POINT( NPLT,XP,YP,ZP,PLAZM,PLDIP,PLUNJ,plngth1,plngth2,PWDTH1,XCNTR,YCNTR,PLTOP)
             CALL SET_CELLS (KP,NPLT,NLYR,MXAB,MXB,CELLW,THK,PLNGTH,PLWDTH,NA,NB,DA,DB, &
                  XCNTR,YCNTR,PLTOP,PLAZM,PLDIP,XCELL,YCELL,ZCELL,PCNR)
             CALL GET_FWD_MODL(LP,PARFAC)
             if (leroiair_failure_count>0) return

             PWDTH1(JP) = XP0
          end if

          LP = K0 + 8
          if (ijac(8)==1) then
             XP0 = PWDTH2(JP)
             PARFAC=PARFAC0
             PWDTH2(JP) = (1. + PARFAC) * PWDTH2(JP)
             PLNGTH=PLNGTH1+PLNGTH2
             PLWDTH=PWDTH1+PWDTH2
             CALL GET_P223_PLATE_REFERENCE_POINT(NPLT, XP,YP,ZP,PLAZM,PLDIP,PLUNJ,plngth1,plngth2,PWDTH1,XCNTR,YCNTR,PLTOP)
             CALL SET_CELLS (KP,NPLT,NLYR,MXAB,MXB,CELLW,THK,PLNGTH,PLWDTH,NA,NB,DA,DB, &
                  XCNTR,YCNTR,PLTOP,PLAZM,PLDIP,XCELL,YCELL,ZCELL,PCNR)
             CALL GET_FWD_MODL(LP,PARFAC)
             if (leroiair_failure_count>0) return

             PWDTH2(JP) = XP0
          end if

          LP = K0 + 9
          if (ijac(9)==1) then
             XP0 = PTHK(JP)
             PARFAC=PARFAC0
             PTHK(JP) = (1. + PARFAC) * PTHK(JP)
             CALL GET_FWD_MODL(LP,PARFAC)
             if (leroiair_failure_count>0) return

             PTHK(JP) = XP0
          end if


          LP = K0 + 10
          if (ijac(10)==1) then
             XP0 = PLAZM(JP)
             PLAZM(JP) = PLAZM(JP) + DELPHI
             PLNGTH=PLNGTH1+PLNGTH2
             PLWDTH=PWDTH1+PWDTH2
             CALL GET_P223_PLATE_REFERENCE_POINT(NPLT, XP,YP,ZP,PLAZM,PLDIP,PLUNJ,plngth1,plngth2,PWDTH1,XCNTR,YCNTR,PLTOP)
             CALL SET_CELLS (KP,NPLT,NLYR,MXAB,MXB,CELLW,THK,PLNGTH,PLWDTH,NA,NB,DA,DB, &
                  XCNTR,YCNTR,PLTOP,PLAZM,PLDIP,XCELL,YCELL,ZCELL,PCNR)
             PPARFAC=DELPHI
             CALL GET_FWD_MODL(LP,PPARFAC)
             if (leroiair_failure_count>0) return

             PLAZM(JP) = XP0
          end if

          LP = K0 + 11
          if (ijac(11)==1) then
             XP0 = PLDIP(JP)
             PLDIP(JP) = PLDIP(JP) + DELPHI
             PLNGTH=PLNGTH1+PLNGTH2
             PLWDTH=PWDTH1+PWDTH2
             CALL GET_P223_PLATE_REFERENCE_POINT(NPLT, XP,YP,ZP,PLAZM,PLDIP,PLUNJ,plngth1,plngth2,PWDTH1,XCNTR,YCNTR,PLTOP)
             CALL SET_CELLS (KP,NPLT,NLYR,MXAB,MXB,CELLW,THK,PLNGTH,PLWDTH,NA,NB,DA,DB, &
                  XCNTR,YCNTR,PLTOP,PLAZM,PLDIP,XCELL,YCELL,ZCELL,PCNR)
             PPARFAC=DELPHI
             CALL GET_FWD_MODL(LP,PPARFAC)
             if (leroiair_failure_count>0) return

             PLDIP(JP) = XP0
          end if

          LP = K0 + 12
          if (ijac(12)==1) then
             XP0 = PLUNJ(JP)
             PLUNJ(JP) = PLUNJ(JP) + DELPHI
             PLNGTH=PLNGTH1+PLNGTH2
             PLWDTH=PWDTH1+PWDTH2
             CALL GET_P223_PLATE_REFERENCE_POINT(NPLT, XP,YP,ZP,PLAZM,PLDIP,PLUNJ,plngth1,plngth2,PWDTH1,XCNTR,YCNTR,PLTOP)
             CALL SET_CELLS (KP,NPLT,NLYR,MXAB,MXB,CELLW,THK,PLNGTH,PLWDTH,NA,NB,DA,DB, &
                  XCNTR,YCNTR,PLTOP,PLAZM,PLDIP,XCELL,YCELL,ZCELL,PCNR)
             PPARFAC=DELPHI
             CALL GET_FWD_MODL(LP,PPARFAC)
             if (leroiair_failure_count>0) return
             PLUNJ(JP) = XP0
          end if

          ! reset things
          PARFAC=PARFAC0
          PLNGTH=PLNGTH1+PLNGTH2
          PLWDTH=PWDTH1+PWDTH2
          CALL GET_P223_PLATE_REFERENCE_POINT(NPLT, XP,YP,ZP,PLAZM,PLDIP,PLUNJ,plngth1,plngth2,PWDTH1,XCNTR,YCNTR,PLTOP)

          CALL SET_CELLS (KP,NPLT,NLYR,MXAB,MXB,CELLW,THK,PLNGTH,PLWDTH,NA,NB,DA,DB, &
               XCNTR,YCNTR,PLTOP,PLAZM,PLDIP,XCELL,YCELL,ZCELL,PCNR)
       END DO

       K0 = 12*NPLT
       if (ijac(14)==1) then

          DO JL = 1,NLYR-1
             LP = K0 + JL
             DO JS=1,NSTAT
                ZP0(JS)=RES(JL+NLYR*(JS-1))
                RES(JL+NLYR*(JS-1))=(1. + PARFAC) * RES(JL+NLYR*(JS-1))
             END DO
             CALL GET_FWD_MODL(LP,PARFAC)
             DO JS=1,NSTAT
                RES(JL+NLYR*(JS-1))=ZP0(JS)
             END DO

          END DO
       end if

       if (ijac(15)==1) then
          JL=NLYR
          LP = K0 + JL
          DO JS=1,NSTAT
             ZP0(JS)=RES(JL+NLYR*(JS-1))
             RES(JL+NLYR*(JS-1))=(1. + PARFAC) * RES(JL+NLYR*(JS-1))
          END DO
          CALL GET_FWD_MODL(LP,PARFAC)
          DO JS=1,NSTAT
             RES(JL+NLYR*(JS-1))=ZP0(JS)
          END DO

       end if

       if (ijac(13)==1) then

          DO JL = 1,NLYR-1
             LP = K0 + NLYR + JL
             DO JS=1,NSTAT
                ZP0(JS)=THK(JL+(NLYR-1)*(JS-1))
                THK(JL+(NLYR-1)*(JS-1))=(1. + PARFAC) * THK(JL+(NLYR-1)*(JS-1))
             END DO
             CALL GET_FWD_MODL(LP,PARFAC)
             DO JS=1,NSTAT
                THK(JL+(NLYR-1)*(JS-1))=ZP0(JS)
             END DO
          END DO
       end if
    END IF




    XMODL = XMODL0
  CONTAINS

    SUBROUTINE GET_FWD_MODL(LP,PARFAC)
      !  ----------------------------------

      !***  Called by: FORJAC
      !***      Calls: LEROI_3D, TEM_3D, HSBOSS_TD, HSBOSS_FD

      !  If LP = 0, performs model computation using existing parameters
      !  If LP > 0, performs model computation where parameter LP is multiplied by
      !             PARFAC and Jacobian column LP isconstructed.



      INTEGER LP
      REAL PARFAC,CSTX(NRXST),SNTX(NRXST)
      LOGICAL INVERT
      DATA INVERT /.TRUE./
      REAL SIG_T(NPLT)
      LOGICAL INTRUDE
      DATA INTRUDE /.FALSE./
      ! convert plate thickness and resistivity into conductance
      do i=1,NPLT
         sig_t(i)=pthk(i)/pres(i)
      end do

      !print *,'>>>>',JP,LP
      !print *,">>>>",nstat

	   !print *,ijac(12),ABS(PLUNJ)
      if (ijac(12)==0 .AND. ALL(ABS(PLUNJ)<1e-3)) then

            CALL LEROI_3D_adv (TDFD,NFRQ,FREQ,NLYR,THK,RES,PBRES,RMU,REPS,CALF,CTAU,CFREQ,NPLT,MXAB,  &
                 MXB,NB,NA,DA,DB,XCELL,YCELL,ZCELL,PLTOP,PLWDTH,PLNGTH,XCNTR,YCNTR, &
                 PLAZM,PLDIP,SIG_T,CALFP,CTAUP,CFREQP,NSTAT,FANGLE,SX,SY,SZ,TXCLN,  &
                 TXA90,SAME_TX,NRX,NRXST,RX,RY,RZ,BFD_SCAT,leroiair_failure_count)

      else

         CALL FRNK_LEROI_3D (NFRQ,FREQ,NSTAT,SX,SY,SZ,FANGLE, &
              RX,RY,RZ,NLYR,THK,RES,REPS,CALF+1,CTAU,CFREQ, &
              NPLT,MXAB,CELLW,PLNGTH,PLWDTH,XCNTR,YCNTR,PLTOP,PLAZM,PLDIP,PLUNJ,  &
              INTRUDE,SIG_T,CALFP+1,CTAUP,CFREQP,BFD_SCAT,leroiair_failure_count)
      end if


!CALL PRINT_LEROI_3D_VARS (TDFD,NFRQ,FREQ,NLYR,THK,RES,RMU,REPS,CALF,CTAU,CFREQ,NPLT,MXAB,  &
!       MXB,NB,NA,DA,DB,XCELL,YCELL,ZCELL,PLTOP,PLWDTH,PLNGTH,XCNTR,YCNTR, &
!       PLAZM,PLDIP,SIG_T,CALFP,CTAUP,CFREQP,NSTAT,FANGLE,SX,SY,SZ,TXCLN,  &
!       TXA90,SAME_TX,NRX,NRXST,RX,RY,RZ,BFD_SCAT,leroiair_failure_count)




      if (leroiair_failure_count>0) return


      IF (TDFD < 2) THEN

         !       BTD(JT,JS,JC) - total field for channel JT, source position JS,component JC
         !  BTD_SCAT(JT,JS,JC) - scattered field: units = nT or nT/s for STEP = 1 or 0 respectively.
         !          JC = 1 => in-line component; 2 => transverse component;  3 => vertical component


         CALL TDEM_3D (STEP,IDER,NSX,SWX,SWY,NPULS,PULSE,NTYPLS,NTYRP,TRP,NCHNL, &
              TOPN,TCLS,FREQ,NFRQ,NSTAT,BFD_SCAT,GSTRP,ASTRP,BTD_SCAT)




         CALL HSBOSS_TD_adv (STEP,IDER,NSX,SWX,SWY,NPULS,PULSE,NTYPLS,NTYRP,TRP,NCHNL,TOPN, &
              TCLS,TXCLN,NSTAT,SZ,ZRX,XRX,YRX,NLYR,RES,REPS,RMU,THK, &
              CALF,CTAU,CFREQ,GSTRP,ASTRP,BTD)

        ! call PRINT_HSBOSS_TD_VARS (STEP,IDER,NSX,SWX,SWY,NPULS,PULSE,NTYPLS,NTYRP,TRP,NCHNL,TOPN, &
        ! TCLS,TXCLN,NSTAT,SZ,ZRX,XRX,YRX,NLYR,RES,REPS,RMU,THK, &
        ! CALF,CTAU,CFREQ,GSTRP,ASTRP,BTD)

         BTD = BTD + BTD_SCAT

         IF (CMP == 13 .OR. CMP == 2 .OR. CMP == 3) THEN
            DO JS = 1,NSTAT
               KS1 = (JS - 1) * MCHNL + 1
               KS2 = KS1 + NCHNL -1
               XMODL(KS1:KS2) = NORM(3) * BTD(1:NCHNL,JS,3)
            END DO
            IF (CMP == 2 .OR. CMP == 3) THEN
               DO JS = 1,NSTAT
                  KS1 = (JS - 1) * MCHNL + NCHNL + 1
                  KS2 = KS1 + NCHNL -1
                  XMODL(KS1:KS2) = NORM(1) * BTD(1:NCHNL,JS,1)
               END DO
               IF (CMP == 3) THEN
                  DO JS = 1,NSTAT
                     KS1 = (JS - 1) * MCHNL + 2*NCHNL + 1
                     KS2 = KS1 + NCHNL -1
                     XMODL(KS1:KS2) = NORM(2) * BTD(1:NCHNL,JS,2)
                  END DO
               END IF
            END IF
         ELSE IF (CMP == 11) THEN
            DO JS = 1,NSTAT
               KS1 = (JS - 1) * MCHNL + 1
               KS2 = KS1 + NCHNL -1
               XMODL(KS1:KS2) = NORM(1) * BTD(1:NCHNL,JS,1)
            END DO
         ELSE IF (CMP == 4 .OR. CMP == 42) THEN
            DO JS = 1,NSTAT
               KS1 = (JS - 1) * MCHNL
               DO JT = 1,NCHNL
                  X2 = BTD(JT,JS,1)**2 + BTD(JT,JS,3)**2
                  X3 = X2 + BTD(JT,JS,2)**2
                  IF (CMP == 42) THEN
                     XMODL(JT + KS1) = SQRT (X2)
                  ELSE
                     XMODL(JT + KS1) = SQRT (X3)
                  END IF
               END DO
            END DO
         END IF

      ELSE
         CALL HSBOSS_FD (NFRQ,FREQ,TXCLN,TXA90,NSTAT,SAME_TX,SZ,ZRX,XRX, &
              YRX,NLYR,RES,REPS,RMU,THK,CALF,CTAU,CFREQ,BFD)

         BFD = BFD - BFD_SCAT  !  Redefine BFD as the total response,

         !      The "-" sign is a consequence of using the sine transform for a +iwt
         !      sign convention.  It is thus consistant with the convention used for the
         !      layered half space and in TDEM for time-domain scattered fields.

         !  Store maximally coupled components in XMODL

         CSTX(1:NFRQ) = COS (TXCLN(1:NFRQ))
         SNTX(1:NFRQ) = SIN (TXCLN(1:NFRQ))



         DO JS = 1,NSTAT
            KS1 = (JS - 1) * 2*NFRQ
            DO JF = 1,NFRQ
               IF (TXA90) THEN
                  XBFD = NORM(JF) * BFD(JF,JS,2)
               ELSE
                  XBFD = NORM(JF) * (BFD(JF,JS,1) * SNTX(JF) + BFD(JF,JS,3) * CSTX(JF))
               END IF
               JD = JF + KS1
               XMODL(JD)        =  REAL (XBFD)
               XMODL(JD + NFRQ) = AIMAG (XBFD)
            END DO
         END DO
      END IF




      IF (LP > 0) THEN
         !   print *,'>>>>'
         DO JD = 1, NDATA
            VM = XMODL0(JD)
            VJ = XMODL(JD)
            DELTA = XWTS(JD) * (VJ - VM)
            IF (INRM == 1) THEN
               DENOM = SQRT ((VM**2 + VJ**2)/ 2.)
               !IF (ABS(DELTA) > 1.E-8 * DENOM) A(JD,LP) = DELTA / (PARFAC * DENOM)   //JRH
               !print *,JD+(LP-1)*NDATA
               ! print *,DELTA, DELTA / (PARFAC * DENOM)
               IF (ABS(DELTA) > 1.E-8 * DENOM) A(JD+(LP-1)*NDATA) = DELTA / (PARFAC * DENOM)
            ELSE IF (INRM == 2) THEN
               IF (ABS(DELTA) > 1.E-8 * DNORM(JD)) A(JD+(LP-1)*NDATA)  = DELTA / (PARFAC * DNORM(JD))
            END IF


         END DO
      END IF



    END SUBROUTINE GET_FWD_MODL

  END SUBROUTINE FORJAC2



  SUBROUTINE COLRES_adv(FRQ,NLYR,NPLT,NSTAT,RES_NSTAT,PBRES,REPS,CALF,CTAU,CFREQ,SIG_T,CALFP,CTAUP, &
       CFREQP,SIGL,RMU,KSQ_LYR,SIGT,KSQ_SHT)
    !---------------------------------------------------------------------------

    !  Computes SIGL, the complex conductivities of layers and SIGT, the complex
    !  conductances of plates, at frequency FRQ using the layered earth Cole-Cole
    !  parameters, CALF, CTAU, CFREQ, the plate Cole-Cole parameters: CALFP,
    !  CTAUP, CFREQP,
    !  the NLYR real layer resistivities, RES, and the real conductance, SIG_T,
    !  SIGL and SIGT include displacement currents.

    !  KSQ_LYR = iwu * SIGL = the layered earth propagation constants.
    !  KSQ_SHT = iwu * SIGT = the propagation constants for plates.

    !***  Called by: LEROI_3D
    use leroiair_subroutines

    IMPLICIT NONE
    REAL, PARAMETER :: TWOPI=6.2831853, MU0=12.56637E-7, EPS0=8.854156E-12
    COMPLEX, PARAMETER :: ONE=(1.,0.), CI=(0.,1.)
    INTEGER J,NLYR,NPLT,I
    REAL, DIMENSION(NLYR*NSTAT) :: RES_NSTAT
    REAL :: PBRES
    integer nstat
    REAL, DIMENSION(NLYR) :: RES,RMU,REPS,CALF,CTAU,CFREQ
    REAL OMEGA,MU,EPS,FRQ,SIG_T(NPLT),CALFP(NPLT),CTAUP(NPLT),CFREQP(NPLT)
    COMPLEX SIGT(NPLT),KSQ_SHT(NPLT),SIGL(NLYR),KSQ_LYR(NLYR,NSTAT),P

    INTENT (IN) FRQ,NLYR,NPLT,RMU,RES_NSTAT,REPS,CALF,CTAU,CFREQ,SIG_T,CALFP,CTAUP,CFREQP
    INTENT (OUT) SIGT,KSQ_SHT,SIGL,KSQ_LYR

    OMEGA = TWOPI * FRQ
    DO I=1,NSTAT
       RES(1:NLYR)=RES_NSTAT(1+(I-1)*NLYR:(I)*NLYR)
       RES(NLYR)=PBRES
       SIGL = CMPLX ((1. / RES), 0.)
       SIGT = CMPLX (SIG_T, 0.)

       ! Compute complex conductivity using Cole-Cole parameters if appropriate

       DO J = 1,NLYR                ! Layers
          P = (CI * OMEGA * CTAU(J) )**CFREQ(J)
          SIGL(J) = SIGL(J) * (ONE + P) / (ONE + CALF(J)*P)
          MU = MU0 * RMU(J)
          EPS = EPS0 * REPS(J)
          SIGL(J) = SIGL(J) + CI * OMEGA * EPS  !  Add in displacement term
          KSQ_LYR(J,I) = CI * OMEGA * MU * SIGL(J)
       END DO

    END DO  ! nstat loop

    DO J = 1,NPLT                ! Plates
       P = (CI * OMEGA * CTAUP(J) )**CFREQP(J)
       SIGT(J) = SIGT(J) * (ONE + P) / (ONE + CALFP(J)*P)
    END DO
    KSQ_SHT = CI * OMEGA * MU0 * RMU(NLYR) * SIGT

  END SUBROUTINE COLRES_adv



  SUBROUTINE PRM_BOSS_adv(TDFD,NRXST,TXCLN,TXA90,JF,FRQ,NSTAT,SAME_TX,SX,SY,SZ,FANGLE, &
       NLYR,RMU,KSQ_LYR_STAT,THK,NPLT,PLAZM,PLDIP,MXB,MXAB,NB,NA,XCELL,  &
       YCELL,ZCELL,NRPRM,RHOTRP,E_PRYM)
    ! --------------------------------------------------------------------------------

    !***  Called by: LEROI_3D
    !***      Calls: PRMHNK, CUBSPL, CUBVAL

    !  For a unit airborne dipole source of dip TXCLN, and azimuth, FANGLE, HSBOSS
    !  computes the tangential components of the primary electric field, E_PRYM,
    !  in NPLT plates lying in a uniform half-space or in the basement of a host
    !  beneath NLYR horizontal layers.  RES contains the resistivities of
    !  the overburden and basement.  The target (plate) is discretised along
    !  strike into NAL cells of length DA and down dip into NBL cells of width DB.
    !
    !  For a vertical plate along the X axis, the cells are numbered along strike
    !  from the South to the North, first along the top row and then in the same
    !  direction along lower rows.  The last cell would be in the bottom
    !  North corner.  In plan view, the X and Y components of E_PRYM are positive
    !  to the North and East respectively.
    !
    !      TDFD = 1 for time domain;   = 2 for frequency domain
    !     NRXST = NSTAT for time domain;   = NFRQ for frequency domain
    !     TXCLN - TX inclination angle in the vertical plane along flight path
    !     TXA90 - true for vertical co-planar briadside array
    !        JF - frequency index
    !       FRQ - frequency
    !   SAME_TX - used to eliminate repeat computations
    !     NSTAT - number of flight line transmitter positions
    !     SX,SY - North, East coordinates for transmitter station
    !        SZ - source altitude
    !    FANGLE - flight path angle in radians. (north = 0; east = PI/2)
    !      NLYR - number of layers including basement.
    !   KSQ_LYR - iwu * SIG for all layers
    !       THK - layer thicknesses
    !      NPLT - number of plates
    !     PLAZM - strike angle in radians
    !     PLDIP - dip angle (PI/2 for vertical plate)
    !       MXB - maximum number of cells down dip
    !      MXAB - maximum number of cells in one plate
    !    NA, NB - number of cells along strike & down dip respectively for each plate
    ! XCELL(k,*) - north coordinate of centre of cell k
    ! YCELL(k,*) - east coordinate of centre of cell k
    !     NRPRM - number of horizontal interpolation points needed for primary field.
    !    RHOTRP - 15 points per decade interpolation array
    !   ZCELL(i,*) - depth of cell centre in row i relative to surface.
    !  E_PRYM(1,*) - complex primary electric fields parallel to strike.
    !  E_PRYM(2,*) - complex primary electric fields perpendicular to strike.
    !
    !    The first NAB components of E_PRYM are the fields along strike and
    !    the second group of NA*NB components are along dip.  Thus for a strike
    !    along the X axis, the first NA*NB components would be the X component
    !    of the primary electric field and the second NA*NB components
    !    would be the Y component multiplied by the cosine of the dip angle.
    ! -------------------------------------------------------------------------

    USE LEROIAIR_FILTER_COEFFICIENTS
        use leroiair_subroutines

    IMPLICIT NONE
    REAL, PARAMETER :: MU0=12.56637E-7, PI2= 1.570796
    INTEGER TDFD,MXB,MXAB,NPLT,JP,JS,JF,NLYR,NBL,NAL,NSTAT,JB,JA,JAB,NRPRM,NA(NPLT), &
         NB(NPLT),NRXST,NINTG,JQ
    REAL FRQ,TXCLF,TXCLN(NRXST),CSTX,SNTX,SX(NSTAT),SY(NSTAT),SZ(NSTAT),FANGLE(NSTAT), &
         CSF,SNF,ALT,THK(NLYR),PLAZM(NPLT),PLDIP(NPLT),XCELL(MXAB,NPLT),YCELL(MXAB,NPLT), &
         ZCELL(MXB,NPLT),RHO,YR,RHOTRP(NRPRM),SSTRL,CSTRL,YR2,YR3,YI,YI2,YI3,XB2,XBAR,XBAR0, &
         YBAR,YBAR0,RMU(NLYR),DELT
    REAL, DIMENSION (4,NRPRM) :: QR1,QR2,QR3,QI1,QI2,QI3
    COMPLEX EFAC,KSQ_LYR_STAT(NLYR,NSTAT),KSQ_LYR(NLYR),EX,EY,E_PRYM(2,MXAB,NSTAT,NPLT),QQ(3)
    COMPLEX, DIMENSION(:,:,:), ALLOCATABLE :: FXX,KER
    LOGICAL TXA90,SAME_TX(NSTAT)

    INTENT (OUT) E_PRYM

    !  Set up transmitter sine and cosine

    NINTG = 1
    DELT = 0.
    IF (TXA90) DELT = PI2  ! For VCPB array rotate transmitter azimuth by PI/ 2

    ALLOCATE (FXX(NRPRM,3,MXB), KER(JNLO-NRPRM:JNHI,2,MXB) )
    FXX=(0.,0.); KER=(0.,0.)

    EFAC = (0.,1.) * FRQ * MU0 / 2.   ! iwu / 4 PI u0 used because source is in air.
    E_PRYM = (0.,0.)

    PLATE_LOOP: DO JP = 1, NPLT

       NAL = NA(JP)  ! Set up local dimensions for one-time plate operations
       NBL = NB(JP)

       !  Set up array FXX which contains values of the F-potential as a function of
       !  of horizontal distance and depth.  Step through different altitude levels.

       STN_LOOP: DO JS = 1,NSTAT
          KSQ_LYR=KSQ_LYR_STAT(:,JS)
          CSF = COS (FANGLE(JS) + DELT)
          SNF = SIN (FANGLE(JS) + DELT)
          SNTX = 0.
          CSTX = 1.
          JQ = JF
          IF (TDFD < 2) JQ = JS
          TXCLF = TXCLN(JQ)
          IF (ABS (TXCLF) .GT. 1.E-4) THEN
             NINTG = 3
             SNTX = SIN (TXCLF)
             CSTX = COS (TXCLF)
          END IF

          ALT = SZ(JS)
          IF (.NOT. SAME_TX(JS)) &
               CALL PRMHNK (NRPRM,RHOTRP,NPLT,JP,MXB,NBL,ZCELL,ALT,NLYR,THK,RMU,KSQ_LYR, &
               NINTG,KER,FXX)

          !  Step through depths, splining FXX for each depth as a function of RHO, the
          !  horizontal distance.  Add up contributions for each dipole for each
          !  transmitter.

          DEPTH_STEP: DO JB = 1, NBL
             QR1(1,1:NRPRM) = REAL (FXX(1:NRPRM,1,JB))
             QI1(1,1:NRPRM) = AIMAG (FXX(1:NRPRM,1,JB))
             CALL CUBSPL (RHOTRP,QR1,NRPRM,0,0)
             CALL CUBSPL (RHOTRP,QI1,NRPRM,0,0)
             IF (NINTG == 3) THEN
                QR2(1,1:NRPRM) = REAL (FXX(1:NRPRM,2,JB))
                QI2(1,1:NRPRM) = AIMAG (FXX(1:NRPRM,2,JB))
                QR3(1,1:NRPRM) = REAL (FXX(1:NRPRM,3,JB))
                QI3(1,1:NRPRM) = AIMAG (FXX(1:NRPRM,3,JB))
                CALL CUBSPL (RHOTRP,QR2,NRPRM,0,0)
                CALL CUBSPL (RHOTRP,QI2,NRPRM,0,0)
                CALL CUBSPL (RHOTRP,QR3,NRPRM,0,0)
                CALL CUBSPL (RHOTRP,QI3,NRPRM,0,0)
             END IF

             STRYK_STEP: DO JA = 1, NAL
                JAB = JA + (JB-1)*NAL
                EX = (0.,0.)
                EY = (0.,0.)
                QQ = (0.,0.)
                RHO = SQRT ( (XCELL(JAB,JP) - SX(JS) )**2 + (YCELL(JAB,JP) - SY(JS) )**2)

                !  Compute components in system with x along flight path centred at
                !  transmitter position.

                IF (RHO > 0.) THEN
                   XBAR0 = (XCELL(JAB,JP) - SX(JS)) / RHO
                   YBAR0 = (YCELL(JAB,JP) - SY(JS)) / RHO
                   XBAR = XBAR0 * CSF + YBAR0 * SNF
                   YBAR = YBAR0 * CSF - XBAR0 * SNF
                   XB2 = XBAR**2
                   YR = CUBVAL (RHOTRP,QR1,NRPRM,RHO)
                   YI = CUBVAL (RHOTRP,QI1,NRPRM,RHO)
                   QQ(1) = CMPLX (YR, YI)
                   EX =  CSTX * QQ(1) * YBAR
                   EY = -CSTX * QQ(1) * XBAR
                   IF (NINTG == 3) THEN
                      YR2 = CUBVAL (RHOTRP,QR2,NRPRM,RHO)
                      YI2 = CUBVAL (RHOTRP,QI2,NRPRM,RHO)
                      YR3 = CUBVAL (RHOTRP,QR3,NRPRM,RHO)
                      YI3 = CUBVAL (RHOTRP,QI3,NRPRM,RHO)
                      QQ(2) = CMPLX (YR2, YI2)
                      QQ(3) = CMPLX (YR3, YI3)
                      EX = EX + SNTX * XBAR * YBAR * (2.*QQ(3) - QQ(2))
                      EY = EY + SNTX * ((1. - 2.*XB2) * QQ(3) + XB2 * QQ(2))
                   END IF
                   IF (TXA90) THEN
                      QQ(1) = EX     ! Kluge to allow horizontal in-line dipole Tx results
                      QQ(2) = EY     ! to be used for horizontal transverse dipole Tx
                      EX = -QQ(2)
                      EY =  QQ(1)
                   END IF

                   !  Resolve the field components parallel and perpendicular to strike
                   !  Apportion X component along the strike direction
                   !  Return for cases where energy is attenuated due to depth and reset the
                   !  value of NB for all higher frequencies.

                   CSTRL = COS (PLAZM(JP) - FANGLE(JS))
                   SSTRL = SIN (PLAZM(JP) - FANGLE(JS))
                   E_PRYM(1,JAB,JS,JP) = EFAC * (EX*CSTRL + EY*SSTRL)
                   E_PRYM(2,JAB,JS,JP) = EFAC * (EY*CSTRL - EX*SSTRL) * COS (PLDIP(JP))
                END IF
             END DO STRYK_STEP
          END DO DEPTH_STEP
       END DO STN_LOOP
    END DO PLATE_LOOP

    DEALLOCATE (FXX,KER)

  END SUBROUTINE PRM_BOSS_adv


  SUBROUTINE LEROI_3D_adv (TDFD,NFRQ,FREQ,NLYR,THK_NSTAT,RES_NSTAT,PBRES,RMU,REPS,CALF,CTAU,CFREQ,NPLT,MXAB,  &
       MXB,NB,NA,DA,DB,XCELL,YCELL,ZCELL,PLTOP,PLWDTH,PLNGTH,XCNTR,YCNTR, &
       PLAZM,PLDIP,SIG_T,CALFP,CTAUP,CFREQP,NSTAT,FANGLE,SX,SY,SZ,TXCLN,  &
       TXA90,SAME_TX,NRX,NRXST,RX,RY,RZ,BFD_SCAT,leroiair_failure_count)
    !-----------------------------------------------------------------------------------
    !
    !***  Called by: MAIN, GET_FWD_MODL
    !***      Calls: COLRES, PRM_BOSS, SCAT_MAG, SCAT_MTRX-BOSS, SET_MGT, SET_NCELL2, SET_RHO,

    ! Main routine for LeroiAir Computation.
    ! Note that from 25.11.03 BTD_SCAT changes from: north & east components to
    !                --------          aircraft in-line & transverse components
    !
    !             OUTPUT
    !             ------
    !
    !   BFD_SCAT(JF,JS,JC) - scattered frequency-domain magnetic field in nT for
    !                        frequency, JF, station, JS, component JC
    !
    !                from SCAT_MAG  JC = 1 => north;   = 2 => east;       = 3 => vertical
    !                from LEROI_3D  JC = 1 => in-line; = 2 => transverse; = 3 => vertical
    !
    !             INPUT
    !             -----
    !
    !      TDFD - 1=> time-domain;  2 => frequency-domain
    !      FREQ - array of NFRQ frequencies
    !      NLYR - number of layers (1 or 2)
    !       THK - layer thicknesses
    !       RES - array of layer resistivities
    !       RMU - mu(i) / mu(0)
    !      REPS - relative dielectric constant
    !      CALF, CTAU, CFREQ are the layered earth Cole-Cole parameters.
    !
    !      NPLT - number of thin plates
    !     PLAZM - strike angle - 0 => due north
    !     PLDIP - 0 => horizontal,  90 => vertical
    !     XCNTR, YCNTR, PLTOP: north, east & depth coordinates of plate reference point
    !                           midpoint of top edge (left edge for flat-lying plate
    !     PLNGTH, PLWDTH: strike length and length down dip
    !     SIG_T, CALFP, CTAUP, CFREQP: conductance, C-C parameters for each plate

    !        NB,NA - number of cells down dip and along strike for each plate.
    !        DB,DA - down-dip and along-strike dimensions for these cells
    !         MXAB - Number of cells in biggest plate
    !          MXB - maximum number of cells down dip
    !   XCELL(k,*) - north coordinate of centre of cell k
    !   YCELL(k,*) - east coordinate of centre of cell k
    !   ZCELL(i,*) - depth of cell centre in row i relative to surface.
    !
    !        NSTAT - number of transmitter positions.
    !       FANGLE - flight path angle in radians. (North = 0; East = PI/2)
    !   SX, SY, SZ - north, east and altitude (re gnd level) of transmitter
    !       TXCLN  - transmitter inclination (> 0 => nose up)
    !        TXA90 - true for vertical co-planar briadside array
    !      SAME_TX - used to reduce redundant computations

    !        NRXST = TXCLN dimension: = NSTAT for time domain; = NFRQ for frequency domain
    !          NRX = receiver offset dimension: = 1 for time domain; = NFRQ for frequency domain
    !  RX,RY,RZ(JS,JR) - north, east and vertical receiver positions for station JS, offset JR

    use leroiair_subroutines

    IMPLICIT NONE
    INTEGER, PARAMETER :: NCPTS=3     !  Number of points for MGT gaussian integration
    COMPLEX, PARAMETER :: ZERO=(0.,0.)
    INTEGER TDFD,NFRQ,NLYR,NPLT,MXAB,MXB,MXRHO,NRPRM,NREGT,NRMGT,MXCL2,NSTAT,NRX, &
         IFIN,JF,JS,JQ,NA(NPLT),NB(NPLT),NCELL2(0:NPLT),MXGS,NRXST
    REAL FRQ,FREQ(NFRQ),THK(NLYR),MGTMX,XCELL(MXAB,NPLT),YCELL(MXAB,NPLT),ZCELL(MXB,NPLT), &
         XRM(NSTAT,NRX,NPLT),YRM(NSTAT,NRX,NPLT),TXCLN(NRXST),CSF,SNF
    REAL, DIMENSION(NLYR) :: RES,RMU,REPS,CALF,CTAU,CFREQ
    REAL, DIMENSION(NPLT) :: SIG_T,CALFP,CTAUP,CFREQP,PLTOP,PLWDTH,XCNTR,YCNTR, &
         PLNGTH,PLAZM,PLDIP,DA,DB,WXZY
    REAL, DIMENSION(NSTAT) :: SX,SY,SZ,FANGLE
    REAL, DIMENSION(NSTAT,NRX) :: RX,RY,RZ
    REAL, DIMENSION(:), ALLOCATABLE :: B,RHOTRP
    REAL, DIMENSION(:,:), ALLOCATABLE :: XGS,YGS,ZGS
    LOGICAL SAME_TX(NSTAT),TXA90
    COMPLEX SIGL(NLYR),KSQ_LYR(NLYR),KSQ_LYR_STAT(NLYR,NSTAT),SIGT(NPLT),KSQ_SHT(NPLT),E_PRYM(2,MXAB,NSTAT,NPLT), &
         J_SCAT(2,MXAB,NSTAT,NPLT),XLINE,YTRNS,BFD_SCAT(NFRQ,NSTAT,3)
    LOGICAL DCMP_FAIL
    REAL PBRES
    INTENT (OUT) BFD_SCAT

    integer leroiair_failure_count
    real res_nstat(nlyr*nstat),thk_nstat((nlyr-1)*nstat)


    RES(1:NLYR)=RES_NSTAT(1:NLYR)
    THK(1:NLYR-1)=THK_NSTAT(1:NLYR-1)
    RES(NLYR)=PBRES


    ! Set up horizontal distance descriptors and arrays for
    ! primary field, & electric & magnetic Green's functions.

    ! Set up Gausian integration for magnetic field computation

    !**************************
    IF (NPLT < 1) RETURN
    !**************************

    MXGS = NCPTS * MAX (MAXVAL(NA), MAXVAL(NB))
    ALLOCATE (XGS(MXGS,NPLT),YGS(MXGS,NPLT),ZGS(MXGS,NPLT))

    CALL SET_MGT (NSTAT,NRX,NPLT,NCPTS,NA,NB,DA,DB,PLTOP,PLWDTH,PLNGTH,PLDIP, &
         PLAZM,XCNTR,YCNTR,RX,RY,XRM,YRM,WXZY,MXGS,XGS,YGS,ZGS,MGTMX)

    CALL SET_NCELL2 (NPLT,NA,NB,NCELL2,MXCL2)

    ALLOCATE (B(1000) )
    CALL SET_RHO (MXAB,NPLT,NSTAT,MGTMX,NA,NB,DA,PLNGTH,XCELL,YCELL,SX,SY, &
         NRMGT,NRPRM,NREGT,MXRHO,B)
    ALLOCATE (RHOTRP(MXRHO))
    RHOTRP(1:MXRHO) = B(1:MXRHO)
    DEALLOCATE (B)

    BFD_SCAT = ZERO

    !! IF (.NOT. INVERT) WRITE(*,1) NFRQ  // JRH
    FREQUENCY_LOOP: DO JF = 1,NFRQ

       IFIN = NINT (100. * (JF-1) / REAL (NFRQ) )
       FRQ = FREQ(JF)
       !!  IF (.NOT. INVERT) WRITE(*,2) JF,FRQ,IFIN // JRH

       ! Set up complex conductivities.  Then compute primary fields, Green's tensors,
       ! scattering matrices and scattering currents for each plate individually.

       !CALL COLRES (FRQ,NLYR,NPLT,RES,REPS,CALF,CTAU,CFREQ,SIG_T,CALFP,CTAUP,CFREQP, &
       !             SIGL,RMU,KSQ_LYR,SIGT,KSQ_SHT)

       CALL COLRES_adv (FRQ,NLYR,NPLT,NSTAT,RES_NSTAT,PBRES,REPS,CALF,CTAU,CFREQ,SIG_T,CALFP,CTAUP,CFREQP, &
            SIGL,RMU,KSQ_LYR_STAT,SIGT,KSQ_SHT)
       KSQ_LYR=KSQ_LYR_STAT(:,1)
       !  Compute the layered earth electric fields on the target, E_PRYM.

       E_PRYM = ZERO
       JQ =1
       IF (TDFD == 2) JQ = JF            ! Tx-Rx offset index

       !CALL PRM_BOSS (TDFD,NRXST,TXCLN,TXA90,JF,FRQ,NSTAT,SAME_TX,SX,SY,SZ,FANGLE, &
       !               NLYR,RMU,KSQ_LYR,THK,NPLT,PLAZM,PLDIP,MXB,MXAB,NB,NA,XCELL,  &
       !               YCELL,ZCELL,NRPRM,RHOTRP,E_PRYM)

       CALL PRM_BOSS_adv (TDFD,NRXST,TXCLN,TXA90,JF,FRQ,NSTAT,SAME_TX,SX,SY,SZ,FANGLE, &
            NLYR,RMU,KSQ_LYR_STAT,THK,NPLT,PLAZM,PLDIP,MXB,MXAB,NB,NA,XCELL,  &
            YCELL,ZCELL,NRPRM,RHOTRP,E_PRYM)
       KSQ_LYR=KSQ_LYR_STAT(:,1)


       !  Set up scattering matrix SCAT_MTRX as an LU decomposition.


       CALL SCAT_MTRX_BOSS (NPLT,NLYR,MXAB,MXB,NCELL2,MXCL2,NA,NB,DA,DB,KSQ_SHT,SIGT, &
            SIGL,RMU,KSQ_LYR,THK,PLTOP,PLDIP,PLAZM,XCELL,YCELL,ZCELL, &
            NREGT,RHOTRP,NSTAT,E_PRYM,J_SCAT,DCMP_FAIL)


       IF (DCMP_FAIL) THEN
          !WRITE(NW,3)
          !WRITE(*,3)
          leroiair_failure_count=leroiair_failure_count+1
       !   CALL PRINT_FAILED_MODEL (NLYR,RES,THK,CALF,CTAU,CFREQ,RMU,REPS,NPLT,PLNGTH,PLWDTH, &
       !        SIG_T,CALFP,CTAUP,CFREQP,XCNTR,YCNTR,PLTOP,PLAZM,PLDIP)
          RETURN
       END IF

       !  Compute BFD_SCAT, the scattered frequency-domain magnetic fields.

       CALL SCAT_MAG_adv (JQ,JF,NFRQ,NPLT,MXAB,NA,NB,PLDIP,PLAZM,NLYR,RMU,KSQ_LYR_STAT,THK_NSTAT, &
            NSTAT,SAME_TX,NRX,XRM,YRM,RZ,NRMGT,RHOTRP,NCPTS,MXGS,XGS,YGS,  &
            ZGS,WXZY,J_SCAT,BFD_SCAT)

    END DO FREQUENCY_LOOP

    ! BFD_SCAT(*,*,JC) = the north, east, vertical components for JC = 1, 2, 3 respectively.
    ! Rotate fields by FANGLE into aircraft system so that
    ! BFD_SCAT(*,*,JC) = the in-line, transverse, vertical components for JC = 1, 2, 3 respectively.

    DO JS = 1,NSTAT
       CSF = COS (FANGLE(JS))
       SNF = SIN (FANGLE(JS))
       DO JF = 1,NFRQ
          XLINE = BFD_SCAT(JF,JS,1) * CSF + BFD_SCAT(JF,JS,2) * SNF
          YTRNS = BFD_SCAT(JF,JS,2) * CSF - BFD_SCAT(JF,JS,1) * SNF
          BFD_SCAT(JF,JS,1) = XLINE
          BFD_SCAT(JF,JS,2) = YTRNS
       END DO
    END DO
    !!  print *,'>>>',BFD_SCAT(1,1,1),BFD_SCAT(2,2,2)

    DEALLOCATE (XGS,YGS,ZGS,RHOTRP)


    ! 1 FORMAT(/T3,'A maximum of',I3,' 3D frequency-domain responses', &
    ! /T3,'will be computed initially.'/)
    ! 2 FORMAT(T3,'frequency',I3,'  =',G12.4,I8,' percent done')
    !   FORMAT (//T3,'An evil spirit has entered SCAT_MTRX_LU_DCMP causing the matrix to be singular.', &
    !         /T3,'The model leading to this crash is Described in LeroiAir.out.' &
    !         /T3,'COMPUTATION HALTED.  SEEK HELP.  (art.raiche@csiro.au)')
    !!
    !    CALL PRINT_LEROI_3D_VARS (TDFD,NFRQ,FREQ,NLYR,THK,RES,RMU,REPS,CALF,CTAU,CFREQ,NPLT,MXAB,  &
    !           MXB,NB,NA,DA,DB,XCELL,YCELL,ZCELL,PLTOP,PLWDTH,PLNGTH,XCNTR,YCNTR, &
    !           PLAZM,PLDIP,SIG_T,CALFP,CTAUP,CFREQP,NSTAT,FANGLE,SX,SY,SZ,TXCLN,  &
    !           TXA90,SAME_TX,NRX,NRXST,RX,RY,RZ,BFD_SCAT,leroiair_failure_count)


  END SUBROUTINE LEROI_3D_adv








  SUBROUTINE HSBOSS_TD_adv (STEP,IDER,NSX,SWX,SWY,NPULS,PULSE,NTYPLS,NTYRP,TRP,NCHNL,TOPN, &
       TCLS,TXCLN,NSTAT,SZ,ZRX,XRX,YRX,NLYR,RES_NSTAT,REPS,RMU,THK_NSTAT, &
       CALF,CTAU,CFREQ,GSTRP,ASTRP,BTD)
    !------------------------------------------------------------------------------------

    !***  Called by: MAIN, GET_FWD_MODL
    !***      Calls: COSTRN, CUBSPL HSMD_FD, FOLD_AND_CONVOLVE

    !  Computes BTD, the time-domain layered earth response convolved with the
    !  excitation waveform and the receiver channels per unit receiver area.
    !  For impulse response, it computes dB/dt in nT / s which is the same as
    !  nanovolts per unit area.

    !  For step response, it computes B in nanoteslas.

    !&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    !  SIGN CONVENTION:
    !  ----------------
    !  The normal layered earth field coordinate system used in this
    !  subroutine has X (JC=1) positive along the flight path, Y (JC=2)
    !  positive to starboard, and Z (JC=3) positive down.
    !&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

    !                             INPUT
    !                             -----
    !       STEP = 1 iff step response is to be computed
    !       IDER = 1 if source waveform was dB/dt; = 0 if amps pr B
    !        NSX - number of points used to discretise transmitter signal
    !        SWX - abscissae (seconds) of current waveform
    !        SWY - dI/dt * Tx moment & nanotesla conversion at times SWX
    !      NPULS - number of bipolar pulses of length PULSE
    !      PULSE - length of half-cycle on pulse plus off-time
    !      NTYRP - number of values in TRP for total signal length: 2 * NPULS *PULSE
    !        TRP - array of time values for FD -> TD transformations
    !     NTYPLS - number of TRP values in 1 PULSE
    !      NCHNL - number of channels
    !       TOPN - time at which receiver channel I opens.
    !       TCLS - time at which receiver channel I closes.
    !      TXCLN - angle in radians that TX dipole makes with vertical (climb = +)
    !      NSTAT - number of stations in survey line.
    !    SAME_TX - used to avoid repeat computations
    !         SZ - array of transmitter altitudes
    !        ZRX - vertical offset of RX at each station from transmitter  (below = +)
    !        XRX - in-line horizontal offset of RX at each station J;      (behind = +)
    !        YRX - transverse offset of RX at each station J               (left = +)
    !       NLYR - number of layers
    !        RES - array of layer resistivities
    !       REPS - relative dielectric constant
    !        RMU - mu(i) / mu(0)
    !        THK - array of layer thicknesses
    !     CALF, CTAU, CFREQ are the layered earth Cole-Cole parameters.
    !
    !                             OUTPUT
    !                             ------
    !     BTD(JT,JS,1) - the in-line component of the layered earth response at
    !                    time JT, station JS.
    !     BTD(JT,JS,2) - the horizontal transverse component
    !     BTD(JT,JS,3) - the vertical component

    USE LEROIAIR_FREQUENCY_SELECT
    use leroiair_subroutines

    IMPLICIT NONE
    INTEGER, PARAMETER :: NFRQ=NF_6PDE, NRXF=1, QL=SELECTED_REAL_KIND(12,80)
    REAL, PARAMETER :: TWOPI=6.283185307
    INTEGER STEP,IDER,NSX,NPULS,NTYPLS,NTYRP,NCHNL,NSTAT,NLYR,TDFD,GSTRP,ASTRP, &
         JS,JF,JT,JC
    REAL SWX(NSX),SWY(NSX,3),PULSE,TRP(NTYRP),SZ(NSTAT),ALT,T,YPRM(4,NTYRP),   &
         YCUM(NCHNL),YFRQ(4,NFRQ),FREQ(NFRQ),WF(NFRQ),BTD(NCHNL,NSTAT,3)
    REAL, DIMENSION(NCHNL) :: TOPN,TCLS
    REAL, DIMENSION(NLYR) :: RES,REPS,THK,CTAU,CFREQ,CALF,RMU
    REAL, DIMENSION(NSTAT) :: TXCLN,ZRX,XRX,YRX
    REAL, DIMENSION(NLYR*NSTAT) :: RES_NSTAT
    REAL, DIMENSION((NLYR-1)*NSTAT) :: THK_NSTAT

    REAL(KIND=QL), DIMENSION(NRXF) :: XRXD,YRXD,ZRXD,TXCLND
    COMPLEX(KIND=QL) BFDD(NFRQ,3)
    COMPLEX BFD
    LOGICAL TXA90

    INTENT (IN) STEP,IDER,NSX,SWX,SWY,NPULS,PULSE,NTYPLS,NTYRP,NCHNL,TOPN,TCLS,SZ, &
         TXCLN,NSTAT,ZRX,XRX,YRX,NLYR,THK_NSTAT,RMU,CALF,CTAU,CFREQ
    INTENT (OUT) BTD
    INTENT (INOUT) TRP

    BTD = 0.
    TDFD = 1
    TXA90 = .FALSE.

    FREQ(1:NFRQ) = FRQ_6PDE(1:NFRQ)
    WF(1:NFRQ) = LOG (TWOPI * FREQ(1:NFRQ))

    DO JS = 1, NSTAT
       RES=0.0
       THK=0.0
       RES(1:NLYR)=RES_NSTAT(1+(JS-1)*NLYR:(JS)*NLYR)
       THK(1:NLYR-1)=THK_NSTAT(1+(JS-1)*(NLYR-1):(JS)*(NLYR-1))

       !print *,res
       !print *,thk
       !IF (SAME_TX(JS)) THEN
       !  DO JC = 1,3
       !    BTD(1:NCHNL,JS,JC) = BTD(1:NCHNL,JS-1,JC)
       !  END DO
       !ELSE
       TXCLND(1) = REAL (TXCLN(JS), KIND=QL)
       ALT = SZ(JS)
       ZRXD(1) = REAL (ZRX(JS), KIND=QL)
       XRXD(1) = REAL (XRX(JS), KIND=QL)
       YRXD(1) = REAL (YRX(JS), KIND=QL)

       CALL HSMD_FD (NFRQ,FREQ,ALT,NRXF,TXCLND,TXA90,ZRXD,XRXD,YRXD,NLYR, &
            THK,RES,REPS,RMU,CALF,CTAU,CFREQ,TDFD,BFDD)

       !    Compute BTD, the 'observed' layered earth response by folding the BLEXT,
       !    the extended response over NPULS bipolar cycles into 1 PULSE and then
       !    convolving this with the TX waveform.  It is during the convolution that we
       !    shift from teslas to nanoteslas or nT/s.

       YFRQ = 0.
       DO JC = 1,3
          DO JF = 1,NFRQ
             BFD = CMPLX (BFDD(JF,JC) )
             YFRQ(1,JF) = AIMAG (BFD) / (TWOPI * FREQ(JF) )
          END DO
          CALL CUBSPL (WF,YFRQ,NFRQ,0,0)

          YPRM = 0.
          DO JT = 1, NTYRP   !  Convert to step-function time-domain.
             T = TRP(JT)
             YPRM(1,JT) = COSTRN (WF,YFRQ,NFRQ,NFRQ,T)
          END DO
          CALL FOLD_AND_CONVOLVE (STEP,IDER,NSX,SWX,SWY,NPULS,PULSE,TRP,NTYPLS, &
               NTYRP,NCHNL,TOPN,TCLS,YPRM,GSTRP,ASTRP,YCUM)

          BTD(1:NCHNL,JS,JC) = YCUM(1:NCHNL)
       END DO
       !END IF
    END DO

!CALL PRINT_HSBOSS_TD_VARS (STEP,IDER,NSX,SWX,SWY,NPULS,PULSE,NTYPLS,NTYRP,TRP,NCHNL,TOPN, &
!       TCLS,TXCLN,NSTAT,SZ,ZRX,XRX,YRX,NLYR,RES,REPS,RMU,THK, &
!       CALF,CTAU,CFREQ,GSTRP,ASTRP,BTD)


  END SUBROUTINE HSBOSS_TD_adv




  SUBROUTINE SCAT_MAG_adv (JQ,JF,NFRQ,NPLT,MXAB,NA,NB,PLDIP,PLAZM,NLYR,RMU,KSQ_LYR_STAT,THK_NSTAT,  &
       NSTAT,SAME_TX,NRX,XRM,YRM,RZ,NRMGT,RHOTRP,NCPTS,MXGS,XGS,YGS, &
       ZGS,WXZY,J_SCAT,BFD_SCAT)
    !---------------------------------------------------------------------------------

    !***  Called by: LEROI_3D
    !***      Calls: MGT_BOSS
    !
    !  Computes BFD_SCAT(JF,JS,JC), scattered frequency-domain magnetic field in nT
    !  for frequency, JF, station JS, component JC (1,2,3 => (north, east, vertical)
    !  J_SCAT are the scattered currents.
    !
    !  NFRQ, NPLT, NLYR, NSTAT, NRX: number of frequencies, plates, layers, stations
    !                                and receiver offsets respectively
    !           JQ - Tx-Rx offset index
    !           JF - frequency index
    !         MXAB - Number of cells in biggest plate
    !        NB,NA - number of cells down dip and along strike for each plate.
    ! PLAZM, PLDIP - dip and strike angle for each plate
    !          THK - layer thicknesses
    ! RMU, KSQ_LYR - relative permeability and propagation constant for each layer
    !      SAME_TX - used to reduce redundant computations
    !     XRM, YRM - transformed  receiver positions wrt plate centre for plate JP
    !       RHOTRP - horizontal interpolation array for MGT
    !        NRMGT - number of points in RHOTRP
    !          XGS - location of integration points along strike
    !     YGS, ZGS - coordinates for integration down dip
    !         MXGS - maximum of Gaussian integration points
    !        NCPTS - number of integration points per dimension per cell.
    !         WXZY - total integration weight for each point for plate JP.

    !  J_SCAT (J1,JCL,JS,JP) is the current in the direction J1 in the centre of
    !                        cell JCL of plate JP at station JS.
    !          J1 = 1 => along strike;  = 2 => down dip
    use leroiair_subroutines

    IMPLICIT NONE
    COMPLEX, PARAMETER :: ZERO=(0.,0.)
    INTEGER JQ,JF,NFRQ,JP,NPLT,MXAB,MXGS,NSTAT,NRX,NLYR,NCPTS,NRMGT,JS,J1,JCL
    INTEGER, DIMENSION (NPLT) :: NA,NB
    REAL RHOTRP(NRMGT),THK(NLYR),RZ(NSTAT,NRX),RMU(NLYR),CSTR,SSTR,THK_NSTAT((NLYR-1)*NSTAT)
    REAL, DIMENSION (NPLT) :: PLDIP,PLAZM,WXZY
    REAL, DIMENSION (NSTAT,NRX,NPLT) :: XRM,YRM
    REAL, DIMENSION (MXGS,NPLT) :: XGS,YGS,ZGS
    COMPLEX HA(MXAB,3,NSTAT),HB(MXAB,3,NSTAT),KSQ_LYR(NLYR),TMP(3),TMPR(3), &
         J_SCAT(2,MXAB,NSTAT,NPLT),BFD_SCAT(NFRQ,NSTAT,3),KSQ_LYR_STAT(NLYR,NSTAT)
    LOGICAL SAME_TX(NSTAT)

    !INTENT (INOUT) BFD_SCAT

    !  Compute the magnetic Green's tensor integrals for each plate and combine
    !  them with the scattering currents and sum.


    PLATE_LOOP: DO JP = 1,NPLT

       CSTR = COS (PLAZM(JP))
       SSTR = SIN (PLAZM(JP))


       TX_LOOP: DO JS = 1,NSTAT

          KSQ_LYR=KSQ_LYR_STAT(:,JS)
          THK(1:NLYR-1)=THK_NSTAT(1+(JS-1)*(NLYR-1):(JS)*(NLYR-1))

          CALL MGT_BOSS (JQ,JP,NPLT,MXAB,NA,NB,PLDIP,NLYR,RMU,KSQ_LYR,THK,NSTAT,SAME_TX, &
               NRX,XRM,YRM,RZ,NRMGT,RHOTRP,NCPTS,MXGS,XGS,YGS,ZGS,WXZY,HA,HB)


          !  Multiply the magnetic Green's tensor elements times the strike and
          !  downdip currents.  J1 = 1,2,3 represents the along strike, horizontal
          !  cross strike and vertical components respectively.  This has to be rotated
          !  so that J1 = 1,2,3 represent the North, East and vertical component
          !  respectively. Strike angle is defined as positive, clockwise from North.

          TMP = ZERO; TMPR = ZERO
          DO J1 = 1,3
             DO JCL = 1, NA(JP) * NB(JP)
                TMP(J1) = TMP(J1) + J_SCAT(1,JCL,JS,JP) * HA(JCL,J1,JS) &
                     + J_SCAT(2,JCL,JS,JP) * HB(JCL,J1,JS)
             END DO
          END DO

          !  Rotate the components back into the user specified coordinate system.

          TMPR(1) = TMP(1) * CSTR - TMP(2) * SSTR
          TMPR(2) = TMP(1) * SSTR + TMP(2) * CSTR
          TMPR(3) = TMP(3)

          BFD_SCAT(JF,JS,1:3) = BFD_SCAT(JF,JS,1:3) + TMPR(1:3)

       END DO TX_LOOP
    END DO PLATE_LOOP

  END SUBROUTINE SCAT_MAG_adv




end module leroiair



