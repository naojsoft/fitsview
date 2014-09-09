#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Code for calculating the position of the Hyper-Suprime Cam CCD
# positions, accounting for various distortions of the field.
#
# Satoshi Kawanomoto (kawanomoto.satoshi@nao.ac.jp)
#
# Edited by
# Eric Jeschke (eric@naoj.org)
#
#
import numpy as np
from astro.subaru import SUBARU_LATITUDE_DEG

### fixed parameters
#latitude_d = 19.8   # Mauna Kea latitude in degree
latitude_d = SUBARU_LATITUDE_DEG
MJD_0 = 51544.5     # MJD at J2000.0 (2000/01/1.5)


class SCAGCCDPositions(object):

    def __init__(self):
        self.ccdpos = self.ccdpositions()
        self.dither_ccdpos = self.dithering_ccdpositions() 
        self.vignette_ccdpos = self.vignette_ccdpositions()
        
    def ccdpositions(self):
        ccdpos = np.zeros((4,4,2))

        # There are four guiding CCDs.  These are the positions in pa=0
        # 1
        ccdpos[0,0,0] =  -1.105
        ccdpos[0,0,1] =-270.080
        ccdpos[0,1,0] =  -1.105
        ccdpos[0,1,1] =-270.080 +30.720
        ccdpos[0,2,0] =  -1.105 -63.360
        ccdpos[0,2,1] =-270.080 +30.720
        ccdpos[0,3,0] =  -1.105 -63.360
        ccdpos[0,3,1] =-270.080

        # 2
        ccdpos[1,0,0] =  65.920
        ccdpos[1,0,1] =-270.080
        ccdpos[1,1,0] =  65.920
        ccdpos[1,1,1] =-270.080 +30.720
        ccdpos[1,2,0] =  65.920 -63.360
        ccdpos[1,2,1] =-270.080 +30.720
        ccdpos[1,3,0] =  65.920 -63.360
        ccdpos[1,3,1] =-270.080

        # 3 (AG218)
        ccdpos[2,0,0] =  68.465
        ccdpos[2,0,1] = 239.360
        ccdpos[2,1,0] =  68.465
        ccdpos[2,1,1] = 239.360 +7.68
        ccdpos[2,2,0] =  68.465 -63.360
        ccdpos[2,2,1] = 239.360 +7.68
        ccdpos[2,3,0] =  68.465 -63.360
        ccdpos[2,3,1] = 239.360

        # 4 (AG217)
        ccdpos[3,0,0] =   1.440
        ccdpos[3,0,1] = 239.360
        ccdpos[3,1,0] =   1.440
        ccdpos[3,1,1] = 239.360 +30.720
        ccdpos[3,2,0] =   1.440 -63.360
        ccdpos[3,2,1] = 239.360 +30.720
        ccdpos[3,3,0] =   1.440 -63.360
        ccdpos[3,3,1] = 239.360

        return ccdpos

    def vignette_ccdpositions(self):
        ccdpos = np.zeros((4,4,2))

        no_use = 2.2  # 12x6 arcmin ccd. 100pixel that is equivalent to 0.4 arcimin should not be used for guiding      

        # There are four guiding CCDs.  
        # the positions of dithering area in pa0  
        # 1
        ccdpos[0,0,0] = -1.105 - no_use  
        ccdpos[0,0,1] = -270.080 + no_use
        ccdpos[0,1,0] = -1.105 - no_use
        ccdpos[0,1,1] = -270.080 +30.720 - no_use
        ccdpos[0,2,0] = -1.105 -63.360 + no_use
        ccdpos[0,2,1] = -270.080 +30.720 - no_use
        ccdpos[0,3,0] = -1.105 -63.360 + no_use
        ccdpos[0,3,1] = -270.080 + no_use

        # 2
        ccdpos[1,0,0] = 65.920 - no_use
        ccdpos[1,0,1] = -270.080 + no_use
        ccdpos[1,1,0] = 65.920 - no_use
        ccdpos[1,1,1] = -270.080 +30.720 - no_use
        ccdpos[1,2,0] =  65.920 -63.360 + no_use
        ccdpos[1,2,1] = -270.080 +30.720 - no_use
        ccdpos[1,3,0] = 65.920 -63.360 + no_use
        ccdpos[1,3,1] = -270.080 + no_use

        # 3
        ccdpos[2,0,0] = 68.465 - no_use
        ccdpos[2,0,1] = 239.360 + no_use
        ccdpos[2,1,0] = 68.465 - no_use
        ccdpos[2,1,1] = 239.360 +7.68 - no_use
        ccdpos[2,2,0] = 68.465 -63.360 + no_use
        ccdpos[2,2,1] = 239.360 +7.68 - no_use
        ccdpos[2,3,0] = 68.465 -63.360 + no_use
        ccdpos[2,3,1] = 239.360 + no_use

        # 4
        ccdpos[3,0,0] = 1.440 - no_use
        ccdpos[3,0,1] = 239.360 + no_use
        ccdpos[3,1,0] = 1.440 - no_use
        ccdpos[3,1,1] = 239.360 +30.720 - no_use
        ccdpos[3,2,0] = 1.440 -63.360 + no_use
        ccdpos[3,2,1] = 239.360 +30.720 - no_use
        ccdpos[3,3,0] = 1.440 -63.360 + no_use
        ccdpos[3,3,1] = 239.360 + no_use

        return ccdpos


    def dithering_ccdpositions(self):
        ccdpos = np.zeros((4,4,2))

        no_use = 2 * 5.3  # 12x6 arcmin ccd. extract 2 arcmin for dithering area     

        # There are four guiding CCDs.  
        # the positions of dithering area in pa0  
        # 1
        ccdpos[0,0,0] = -1.105 - no_use  
        ccdpos[0,0,1] = -270.080 + no_use
        ccdpos[0,1,0] = -1.105 - no_use
        ccdpos[0,1,1] = -270.080 +30.720 - no_use
        ccdpos[0,2,0] = -1.105 -63.360 + no_use
        ccdpos[0,2,1] = -270.080 +30.720 - no_use
        ccdpos[0,3,0] = -1.105 -63.360 + no_use
        ccdpos[0,3,1] = -270.080 + no_use

        # 2
        ccdpos[1,0,0] = 65.920 - no_use
        ccdpos[1,0,1] = -270.080 + no_use
        ccdpos[1,1,0] = 65.920 - no_use
        ccdpos[1,1,1] = -270.080 +30.720 - no_use
        ccdpos[1,2,0] =  65.920 -63.360 + no_use
        ccdpos[1,2,1] = -270.080 +30.720 - no_use
        ccdpos[1,3,0] = 65.920 -63.360 + no_use
        ccdpos[1,3,1] = -270.080 + no_use

        # 3
        ccdpos[2,0,0] = 68.465 - no_use
        ccdpos[2,0,1] = 239.360 + no_use
        ccdpos[2,1,0] = 68.465 - no_use
        ccdpos[2,1,1] = 239.360 +7.68 - no_use
        ccdpos[2,2,0] = 68.465 -63.360 + no_use
        ccdpos[2,2,1] = 239.360 +7.68 - no_use
        ccdpos[2,3,0] = 68.465 -63.360 + no_use
        ccdpos[2,3,1] = 239.360 + no_use

        # 4
        ccdpos[3,0,0] = 1.440 - no_use
        ccdpos[3,0,1] = 239.360 + no_use
        ccdpos[3,1,0] = 1.440 - no_use
        ccdpos[3,1,1] = 239.360 +30.720 - no_use
        ccdpos[3,2,0] = 1.440 -63.360 + no_use
        ccdpos[3,2,1] = 239.360 +30.720 - no_use
        ccdpos[3,3,0] = 1.440 -63.360 + no_use
        ccdpos[3,3,1] = 239.360 + no_use

        return ccdpos

    def precessionMatrix(self, MJD):
        TT = (MJD - MJD_0)/36525.0
        zeta_A  = np.deg2rad((2306.2181*TT + 0.30188*TT**2.0 + 0.017998*TT**3.0)/3600.0)
        z_A     = np.deg2rad((2306.2181*TT + 1.09468*TT**2.0 + 0.018203*TT**3.0)/3600.0)
        theta_A = np.deg2rad((2004.3109*TT - 0.42665*TT**2.0 - 0.041833*TT**3.0)/3600.0)
        pM = np.matrix([[+np.cos(zeta_A)*np.cos(theta_A)*np.cos(z_A)-np.sin(zeta_A)*np.sin(z_A), \
                         -np.sin(zeta_A)*np.cos(theta_A)*np.cos(z_A)-np.cos(zeta_A)*np.sin(z_A), \
                         -np.sin(theta_A)*np.cos(z_A)],
                        [+np.cos(zeta_A)*np.cos(theta_A)*np.sin(z_A)+np.sin(zeta_A)*np.cos(z_A), \
                         -np.sin(zeta_A)*np.cos(theta_A)*np.sin(z_A)+np.cos(zeta_A)*np.cos(z_A), \
                         -np.sin(theta_A)*np.sin(z_A)],
                        [+np.cos(zeta_A)*np.sin(theta_A)                                       , \
                         -np.sin(zeta_A)*np.sin(theta_A)                                       , \
                         +np.cos(theta_A)]])
        return pM

    def distortion(self, omega_d):
        r = 319.99905369*omega_d \
            +15.37158789*omega_d**3.0 \
            + 2.61186992*omega_d**5.0 \
            + 4.27786427*omega_d**7.0
        return r

    def distortionInv(self, r):
        omega_d =  0.00312492 *r \
                  -1.447e-09  *r**3.0 \
                  -1.44128e-15*r**5.0 \
                  -1.10711e-20*r**7.0
        return omega_d

    def xy2rt(self, x, y, inr_d):
        r  = np.sqrt(x*x+y*y)
        t0 = np.rad2deg(np.arctan2(y,x))
        t  = inr_d + t0 - 90.0
        return r, t

    def r2rho(self, r, zd_d):
        z = np.deg2rad(zd_d)
        # atmospheric refraction parameter
        R_0 = 1.7325e-04
        # mean magnification parameter
        mm  = 1.0-0.5*R_0*(1.0+1.0/(np.cos(z)**2))
        rho = self.distortionInv(r)/mm
        return rho

    def sp2vec(self, ra_d, dec_d):
        ra  = np.deg2rad(ra_d)
        dec = np.deg2rad(dec_d)
        v   = np.array(([np.cos(ra)*np.cos(dec)],
                        [np.sin(ra)*np.cos(dec)],
                        [np.sin(dec)]))
        return v

    def vec2sp(self, v):
        x = v[0,0]
        y = v[1,0]
        z = v[2,0]
        ra  = np.arctan2(y,x)
        dec = np.arcsin(z)
        ra_d  = np.rad2deg(ra)
        dec_d = np.rad2deg(dec)
        return ra_d, dec_d

    def radeclst2hrzdpa(self, lat_d, ra_d, dec_d, lst_d):
        ra   = np.deg2rad(ra_d)
        dec  = np.deg2rad(dec_d)
        lat  = np.deg2rad(lat_d)
        hr_d = lst_d - ra_d
        hr   = np.deg2rad(hr_d)
        zd   = np.arccos(np.sin(lat)*np.sin(dec)+np.cos(lat)*np.cos(dec)*np.cos(hr)) 
        pa   = np.arctan2(np.cos(lat)*np.sin(hr),np.sin(lat)*np.cos(dec)-np.cos(lat)*np.sin(dec)*np.cos(hr))
        zd_d = np.rad2deg(zd)
        pa_d = np.rad2deg(pa)
        return hr_d, zd_d, pa_d

    def rhoth2hrdec(self, rho_d, th_d, hr_d, dec_d):
        rho = np.deg2rad(rho_d)
        th  = np.deg2rad(th_d)
        hr  = np.deg2rad(hr_d)
        dec = np.deg2rad(dec_d)
        dec2 = np.arcsin(np.cos(rho)*np.sin(dec)+np.sin(rho)*np.cos(dec)*np.cos(th))
        dhr  = np.arcsin(np.sin(th)*np.sin(rho)/np.cos(dec2))
        dec2_d = np.rad2deg(dec2)
        dhr_d  = np.rad2deg(dhr)
        hr2_d = hr_d - dhr_d
        return hr2_d, dec2_d

    def xy2radec(self, x, y, ra_tel_d, dec_tel_d, theta_inr_d, LST_d, MJD):
        pm = self.precessionMatrix(MJD)
        v_tel = self.sp2vec(ra_tel_d, dec_tel_d)
        v_tel_m = np.dot(pm, v_tel)
        ra_tel_m, dec_tel_m = self.vec2sp(v_tel_m)
        hr_tel, zd_tel, pa_tel = self.radeclst2hrzdpa(latitude_d, ra_tel_m, dec_tel_m, LST_d)
        r, t = self.xy2rt(x,y,theta_inr_d)
        rho = self.r2rho(r, zd_tel)
        th = pa_tel - t
        hr_p, dec_p = self.rhoth2hrdec(rho, th, hr_tel, dec_tel_m)
        ra_p = LST_d - hr_p
        v_p = self.sp2vec(ra_p, dec_p)
        v_p_2000 = np.dot(pm.I,v_p)
        ra_p, dec_p = self.vec2sp(v_p_2000)
        return ra_p, dec_p 
        
    def get_ccdpos(self, ra_deg, dec_deg, insrot_deg, lst, mjd):
        """
        Get the coordinates in WCS of the four corners of the HSC
        guiding CCDs.
        
        Takes the telescope pointing ra/dec (in deg), the instrument
        rotation (in deg) and the local sidereal time and modified
        julian time.
        
        Returns a list of lists of tuples, where each tuple is the
        (ra, dec) coordinates for a corner of one of the guiding CCDs.
        """
        res = []
        for i in range(4):
            corners = []
            for j in range(4):
                x = self.ccdpos[i, j, 0]
                y = self.ccdpos[i, j, 1]
    
                ra, dec = self.xy2radec(x, y, ra_deg, dec_deg,
                                        insrot_deg, lst, mjd)
                corners.append((ra, dec))
            res.append(corners)
        return res

    def get_dither_ccdpos(self, ra_deg, dec_deg, insrot_deg, lst, mjd):
        """
        Get the coordinates in WCS of the four corners of the HSC
        guiding CCDs for dithering area.
        
        Takes the telescope pointing ra/dec (in deg), the instrument
        rotation (in deg) and the local sidereal time and modified
        julian time.
        
        Returns a list of lists of tuples, where each tuple is the
        (ra, dec) coordinates for a corner of one of the guiding CCDs.
        """
        res = []
        for i in range(4):
            corners = []
            for j in range(4):
                x = self.dither_ccdpos[i, j, 0]
                y = self.dither_ccdpos[i, j, 1]
    
                ra, dec = self.xy2radec(x, y, ra_deg, dec_deg,
                                        insrot_deg, lst, mjd)
                corners.append((ra, dec))
            res.append(corners)
        return res

    def get_vignette_ccdpos(self, ra_deg, dec_deg, insrot_deg, lst, mjd):
        """
        Get the coordinates in WCS of the four corners of the HSC
        guiding CCDs for non vignetting area.
        
        Takes the telescope pointing ra/dec (in deg), the instrument
        rotation (in deg) and the local sidereal time and modified
        julian time.
        
        Returns a list of lists of tuples, where each tuple is the
        (ra, dec) coordinates for a corner of one of the guiding CCDs.
        """
        res = []
        for i in range(4):
            corners = []
            for j in range(4):
                x = self.vignette_ccdpos[i, j, 0]
                y = self.vignette_ccdpos[i, j, 1]
    
                ra, dec = self.xy2radec(x, y, ra_deg, dec_deg,
                                        insrot_deg, lst, mjd)
                corners.append((ra, dec))
            res.append(corners)
        return res


        
#####
if __name__ == "__main__":
    import pprint
    
    hscobj = SCAGCCDPositions()

    mjd = 51545.0    # Modified JD
    lst = 0          # Local sidereal time
    inr = 20         # instrument rotator angle
    dec = 0          # telescope decl.
    ra  = 0          # telescope R.A.

    ## ccdpos = func.ccdpositions()

    ## for i in range(4):
    ##     for j in range(4):
    ##         x = ccdpos[i, j, 0]
    ##         y = ccdpos[i, j, 1]
    
    ##         print func.xy2radec(x, y, ra, dec, inr, lst, mjd)

    coords = hscobj.get_ccdpos(ra, dec, inr, lst, mjd)
    pprint.pprint(coords)
