# Copyright (C) 2017-2024 CSIRO
#
# This file is part of PyP223
#
# PyP223 is free software: you can redistribute it and/or modify it under the terms of
# the GNU General Public License as published by the Free Software Foundation, either
# version 3 of the License, or (at your option) any later version.
#
# PyP223 is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
# without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
# PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with PyP223.
# If not, see <https://www.gnu.org/licenses/>.
#
# To contact CSIRO about this software you can e-mail
# juerg.hauser@csiro.au
#

import os
import glob
import platform
import numpy
import numpy.ctypeslib
import ctypes

class LeroiAir():
    def __init__(self):
        self.p223=ctypes.cdll.LoadLibrary(glob.glob(os.path.dirname(__file__)+"/p223*.so")[0])

    def formod_tempest_data(self,nlyr,nstat,res,pbres,thk,nplt,peast,pnorth,ptop,pres,
                                         plngth1,plngth2,pwdth1,pwdth2,pthk,cellw,
                                         pdzm,pdip,plng,
                                         ncmp,cmp,
                                         nchnl,topn,tcls,
                                         tx,ty,tz,tazi,tincl,
                                         rx,ry,rz,trdx,trdy,trdz,
                                         xmodl,
                                         leroiair_failure_count):
        
        # Call to LeroiAir to comute the forward response and the Jacobian for
        # a thin plate target and a TEMPEST system
        # nlyr:  Number of layers
        # nstat: Number of stations
        # res:   Layer resistivities
        # pbres: Halfspace resistivity
        # thk: Layer thicknesses
        # nplt: Number of plates
        # peast: Easting of the plate reference point
        # pnorth: Northing of the plate reference point
        # pres: Plate resistivity
        # plnght1: Plate length in the positive direction using a right-hand rule
        # plnght2: Plate length in the negate direction using a right-hand rule
        # pwdith1: Plate width in the positive direction using a right-hand rule
        # pwdith2: Plate width in the negate direction using a right-hand rule
        # pthk: Plate thickness
        # cellw: Cell size used to discreties the thin plate
        # pdzm: Plate dip azimuth
        # pdip: Plate dip
        # plng: Plate plunge
        # ncmp: Number of components
        # cmp: Component selection integer
        # nhcnl: Number of channels/gates
        # topn: Gate start/opening time
        # tcls: Gate closing time
        # nstat: Number of stations/fiducials
        # tx: Transmitter easting
        # ty: Transmitter northing
        # tz: Transmitter elevation
        # tazi: Transmitter azimuth
        # tincl: Tranmitter inclination
        # rx_: Receiver easting
        # ry_: Receiver northing
        # rz_: Receiver elevation
        # trdx_: Transmiter receiver separation where trd(x,y,z) from a right handed coordinate sytem in azimuth direction
        # trdy: Transmiter receiver separation where trd(x,y,z) from a right handed coordinate sytem in azimuth direction
        # trdz: Transmiter receiver separation where trd(x,y,z) from a right handed coordinate sytem in azimuth direction
        # xmodl: Model response
        # leroiair_failure_count: keeping track of instances of LeroiAir failing
      
        nlyr_=ctypes.c_int(nlyr)
        res_=numpy.asfortranarray(res,dtype=numpy.float32)
        pbres_=ctypes.c_float(pbres)
        thk_=numpy.asfortranarray(thk,dtype=numpy.float32)
        nplt_=ctypes.c_int(nplt)
        peast_=numpy.asfortranarray(peast,dtype=numpy.float32)
        pnorth_=numpy.asfortranarray(pnorth,dtype=numpy.float32)
        ptop_=numpy.asfortranarray(ptop,dtype=numpy.float32)
        pres_=numpy.asfortranarray(pres,dtype=numpy.float32)
        plngth1_=numpy.asfortranarray(plngth1,dtype=numpy.float32)
        plngth2_=numpy.asfortranarray(plngth2,dtype=numpy.float32)
        pwdth1_=numpy.asfortranarray(pwdth1,dtype=numpy.float32)
        pwdth2_=numpy.asfortranarray(pwdth2,dtype=numpy.float32)
        pthk_=numpy.asfortranarray(pthk,dtype=numpy.float32)
        cellw_=ctypes.c_float(cellw)
        pdzm_=numpy.asfortranarray(pdzm,dtype=numpy.float32)
        pdip_=numpy.asfortranarray(pdip,dtype=numpy.float32)
        plng_=numpy.asfortranarray(plng,dtype=numpy.float32)
        ncmp_=ctypes.c_int(ncmp)
        cmp_=ctypes.c_int(cmp)
        nchnl_=ctypes.c_int(nchnl)
        topn_=numpy.asfortranarray(topn,dtype=numpy.float32)
        tcls_=numpy.asfortranarray(tcls,dtype=numpy.float32)
        nstat_=ctypes.c_int(nstat)
        tx_=numpy.asfortranarray(tx,dtype=numpy.float32)
        ty_=numpy.asfortranarray(ty,dtype=numpy.float32)
        tz_=numpy.asfortranarray(tz,dtype=numpy.float32)
        tazi_=numpy.asfortranarray(tazi,dtype=numpy.float32)
        tincl_=numpy.asfortranarray(tincl,dtype=numpy.float32)
        rx_=numpy.asfortranarray(rx,dtype=numpy.float32)
        ry_=numpy.asfortranarray(ry,dtype=numpy.float32)
        rz_=numpy.asfortranarray(rz,dtype=numpy.float32)
        trdx_=numpy.asfortranarray(trdx,dtype=numpy.float32)
        trdy_=numpy.asfortranarray(trdy,dtype=numpy.float32)
        trdz_=numpy.asfortranarray(trdz,dtype=numpy.float32)
        xmodl_=numpy.asfortranarray(numpy.zeros(nchnl*ncmp*nstat),dtype=numpy.float32)
        leroiair_failure_count_=ctypes.c_int(leroiair_failure_count)
        
        self.p223.leroiair_formod_tempest_data(
                                       ctypes.byref(nlyr_),
                                       res_.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                                       ctypes.byref(pbres_),
                                       thk_.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                                       ctypes.byref(nplt_),
                                       peast_.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                                       pnorth_.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                                       ptop_.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                                       pres_.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                                       plngth1_.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                                       plngth2_.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                                       pwdth1_.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                                       pwdth2_.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                                       pthk_.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                                       pdzm_.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                                       pdip_.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                                       plng_.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                                       ctypes.byref(nchnl_),
                                       topn_.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                                       tcls_.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                                       ctypes.byref(nstat_),
                                       tx_.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                                       ty_.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                                       tz_.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                                       tazi_.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                                       tincl_.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                                       rx_.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                                       ry_.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                                       rz_.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                                       trdx_.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                                       trdy_.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                                       trdz_.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                                       xmodl_.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                                       ctypes.byref(ncmp_),
                                       ctypes.byref(cmp_),
                                       ctypes.byref(cellw_),
                                       ctypes.byref(leroiair_failure_count_)
         )
        rsp=numpy.zeros([nchnl,ncmp])
        rsp[:,0]=xmodl_[0:nchnl]
        rsp[:,1]=xmodl_[nchnl:nchnl*2]
        return rsp

    def formod_vtem_max_data(self,nlyr,nstat,res,pbres,thk,nplt,peast,pnorth,ptop,pres,
                                         plngth1,plngth2,pwdth1,pwdth2,pthk,cellw,
                                         pdzm,pdip,plng,
                                         ncmp,cmp,
                                         ntrn,txarea,
                                         nchnl,topn,tcls,
                                         nsx,swx,ampt,waveform,
                                         tx,ty,tz,tazi,tincl,
                                         rx,ry,rz,trdx,trdy,trdz,
                                         xmodl,
                                         leroiair_failure_count):
        
        # Call to LeroiAir to comute the forward response and the Jacobian for
        # a thin plate target and a Vtem system
        # nlyr:  number of layers
        # nstat: number of stations
        # res:   resistivity vector
        # pbres: halfspace resistivity
        # thk: layer thickness vector
        # nplt: Number of plates
        # peast:
                
        nlyr_=ctypes.c_int(nlyr)
        res_=numpy.asfortranarray(res,dtype=numpy.float32)
        pbres_=ctypes.c_float(pbres)
        thk_=numpy.asfortranarray(thk,dtype=numpy.float32)
        nplt_=ctypes.c_int(nplt)
        peast_=numpy.asfortranarray(peast,dtype=numpy.float32)
        pnorth_=numpy.asfortranarray(pnorth,dtype=numpy.float32)
        ptop_=numpy.asfortranarray(ptop,dtype=numpy.float32)
        pres_=numpy.asfortranarray(pres,dtype=numpy.float32)
        plngth1_=numpy.asfortranarray(plngth1,dtype=numpy.float32)
        plngth2_=numpy.asfortranarray(plngth2,dtype=numpy.float32)
        pwdth1_=numpy.asfortranarray(pwdth1,dtype=numpy.float32)
        pwdth2_=numpy.asfortranarray(pwdth2,dtype=numpy.float32)
        pthk_=numpy.asfortranarray(pthk,dtype=numpy.float32)
        cellw_=ctypes.c_float(cellw)
        pdzm_=numpy.asfortranarray(pdzm,dtype=numpy.float32)
        pdip_=numpy.asfortranarray(pdip,dtype=numpy.float32)
        plng_=numpy.asfortranarray(plng,dtype=numpy.float32)
        ncmp_=ctypes.c_int(ncmp)
        cmp_=ctypes.c_int(cmp)
        
        ntrn_=ctypes.c_int(ntrn)
        txarea_=ctypes.c_float(txarea)
        
        
        nchnl_=ctypes.c_int(nchnl)
        topn_=numpy.asfortranarray(topn,dtype=numpy.float32)
        tcls_=numpy.asfortranarray(tcls,dtype=numpy.float32)
        nsx_=ctypes.c_int(nsx)
        swx_=numpy.asfortranarray(swx,dtype=numpy.float32)
        ampt_=ctypes.c_int(ampt)    
        waveform_=numpy.asfortranarray(waveform,dtype=numpy.float32)
        nstat_=ctypes.c_int(nstat)
        tx_=numpy.asfortranarray(tx,dtype=numpy.float32)
        ty_=numpy.asfortranarray(ty,dtype=numpy.float32)
        tz_=numpy.asfortranarray(tz,dtype=numpy.float32)
        tazi_=numpy.asfortranarray(tazi,dtype=numpy.float32)
        tincl_=numpy.asfortranarray(tincl,dtype=numpy.float32)
        rx_=numpy.asfortranarray(rx,dtype=numpy.float32)
        ry_=numpy.asfortranarray(ry,dtype=numpy.float32)
        rz_=numpy.asfortranarray(rz,dtype=numpy.float32)
        trdx_=numpy.asfortranarray(trdx,dtype=numpy.float32)
        trdy_=numpy.asfortranarray(trdy,dtype=numpy.float32)
        trdz_=numpy.asfortranarray(trdz,dtype=numpy.float32)
        xmodl_=numpy.asfortranarray(numpy.zeros(nchnl*ncmp*nstat),dtype=numpy.float32)
        leroiair_failure_count_=ctypes.c_int(leroiair_failure_count)
        
        self.p223.leroiair_formod_vtem_max_data(
                                       ctypes.byref(nlyr_),
                                       res_.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                                       ctypes.byref(pbres_),
                                       thk_.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                                       ctypes.byref(nplt_),
                                       peast_.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                                       pnorth_.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                                       ptop_.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                                       pres_.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                                       plngth1_.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                                       plngth2_.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                                       pwdth1_.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                                       pwdth2_.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                                       pthk_.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                                       pdzm_.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                                       pdip_.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                                       plng_.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                                       
                                                                              ctypes.byref(ntrn_),
                                                                              ctypes.byref(txarea_),
                                       ctypes.byref(nchnl_),
                                       topn_.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                                       tcls_.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                                       ctypes.byref(nsx_),
                                       swx_.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                                       ctypes.byref(ampt_),
                                       waveform_.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                                       ctypes.byref(nstat_),
                                       tx_.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                                       ty_.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                                       tz_.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                                       tazi_.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                                       tincl_.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                                       rx_.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                                       ry_.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                                       rz_.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                                       trdx_.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                                       trdy_.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                                       trdz_.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                                       xmodl_.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                                       ctypes.byref(ncmp_),
                                       ctypes.byref(cmp_),
                                       ctypes.byref(cellw_),
                                       ctypes.byref(leroiair_failure_count_)
         )
        rsp=numpy.zeros([nchnl,ncmp])
        rsp[:,0]=xmodl_[0:nchnl]
        rsp[:,1]=xmodl_[nchnl:nchnl*2]
        return rsp
