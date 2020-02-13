# -*- coding: utf-8 -*-
"""
This module is designed to do coordinate transformations using VMEC data

    Original version in MATLAB written by: Christopher A. Clark Creation Date: 06/2012
    Ported to PYTHON by Gavin M. Weir, Creation Date: 03/2016

     File Name: CalcFluxCart_VMEC.py
         CalcFluxCart_VMEC
         FindVMEC_Coords
         GetLabCoordsFromVMEC

    TODO:
    finish porting
    ----- check _dsi.spline, _dsi.spline_deriv, _dsi.ppval
    ----- verify mod( A, B )  == A%B in python

@author: gawe
"""
# ======================================================================== #
# ======================================================================== #

from __future__ import absolute_import, with_statement, absolute_import, \
                       division, print_function, unicode_literals

import numpy as _np
import os as _os
#import scipy.interpolate as _dsi   # broken
#import scipy.integrate as _int

try:
    import netCDF4 as netcdf
except:
    print('Unable to import netCDf4')
# end try

import utils as _ut

try:
    from .VMEC.read_vmec import read_vmec
except:
    from VMEC.read_vmec import read_vmec
# end try

# ========================================================================= #
# ------------------------------------------------------------------------- #
# ========================================================================= #
# GLOBALS
DEBUG = False
VMEC_Data = None
VMEC_DerivedQuant = None
VMEC_DataSource = None

# ========================================================================= #
# FUNCTIONS THAT USE GLOBAL VARIABLES

def Cart2Flux(cLabPos, VMEC_FilePath, ForceReadVMEC=False, PosTol=1e-3, verbose=True):
    """
     Description:
       This function finds the coordinates of a point in VMEC flux coordinates that
       corresponds with a point given in cartesian lab coordinates.  The
       algorithm for this routine is described in detail in the paper:
       Attenberger, Houlberg, and Hirshman, J Comp Phys 72 (1987) 435.

       This version is based on one created by John Schmitt, but is optimized
       for python (vectorized) and somewhat simplified.  It works out to be
       quite fast: 100ms for the first iteration and 10ms for subsequent
       iterations on a 12 year old Core Duo desktop.

       It is worth noting that much of the complexity of this code deals with
       points near the axis.

     Inputs:
       cLabPos - A 2D array describing the points that you want to have
           converted into VMEC coordinates.  The first index is the point
           number, and the second identifies the coordinate.  For example:
          (x1, y1, z1 ; x2, y2, z2;...).  Units are meters.

       VMEC_FilePath - The location and name of the VMEC output file
           that you want to use for the transformation (eithet netcdf or txt).

       ForceReadVMEC (optional) - Set this to one if you want to force the
           routine to read the VMEC file from disk instead of trying to reuse
           data stored in memory.  Otherwise set to zero. Default = 0

       PosTol (optional)  - The routine will complete once it finds a point in
           VMEC coordinates that is within this distance of the point specified
           in cLabPos.  Units are meters. Default is 10^-4.

     Outputs:
       Rho - The approximate normalized poloidal flux (sqrt(S)) at the input
           location
    """
    # Default arguments: We don't typically want to force the reading of the
    # VMEC files, and the tolerance on the position is limited by the radial
    # splining to somewhere around 1e-3, so there is no point in attempting to
    # converge much further than that.

    # ===================================================================== #
    # Fetch VMEC data

    # We'll store the relevant VMEC data in a persistent variables to improve
    # speed for the case where this function is called many times for the same
    # VMEC data.  Unless the user asks for a new file to be used as input, we'll
    # continue to use the persistent data.  Otherwise, load it from a file.
    global VMEC_Data
    global VMEC_DerivedQuant
    global VMEC_DataSource

    # note this test should be for whether the var field is present in VMEC_DataSource
    if ForceReadVMEC or (VMEC_Data is None) or (VMEC_DerivedQuant is None) or (VMEC_DataSource is None) \
        or (VMEC_DataSource != VMEC_FilePath):
        if 1:
#        if type(VMEC_FilePath)!=type('') or (VMEC_Data is None):
            VMEC_Data =  extract_data(VMEC_FilePath, ForceReadVMEC=ForceReadVMEC, verbose=verbose)
#        if VMEC_DerivedQuant is None:
            VMEC_DerivedQuant = spline_data(VMEC_Data, ForceReadVMEC=True, verbose=verbose)
        # end if
    # end
    VMEC_DataSource = VMEC_FilePath

    # ===================================================================== #
    # Use Newton's method to find the values of Rho and Theta that match our
    # position in real space.

    #Convert the user's cartesian data into polar form indexed as r,phi,z
    pLabPos = _np.zeros_like(cLabPos)
    pLabPos[:, 1], pLabPos[:, 0], pLabPos[:, 2] = _ut.cart2pol(cLabPos[:, 0], cLabPos[:, 1], cLabPos[:, 2] )
    Phi = pLabPos[:, 1]

    nrho = len(Phi)
    Rho = _np.zeros((nrho,), _np.float64)
    Theta = _np.zeros_like(Rho)
    for ii in range(nrho):
        Rho[ii], Theta[ii] = FindVMEC_Coords(pLabPos[ii, :], VMEC_Data, VMEC_DerivedQuant, PosTol=PosTol)
    #end for
    return Rho, Theta
#end def CalcFluxCart_VMEC

# ========================================================================= #


def Flux2Cart(sPos, VMEC_FilePath, ForceReadVMEC=False, kMinRho=2.0e-9, jacout=False, verbose=True):
    """
     Description:
       This function finds the coordinates of a point in cartesian lab coordinates
       that corresponds with a point in VMEC flux coordinates.

     Inputs:
       sPos - A 2D array describing the points that you want to have
           converted from VMEC coordinates.  The first index is the point
           number, and the second identifies the coordinate.  For example:
          (rho1, th1, fi1 ; rho2, th2, fi2;...).
              rho ~ normalized poloidal flux (sqrt(S))
              th ~ pseudo-poloidal angle [rad]
              fi ~ pseudo-toroidal angle [rad]

       VMEC_FilePath - The location and name of the VMEC output file
           that you want to use for the transformation (eithet netcdf or txt).

       ForceReadVMEC (optional) - Set this to one if you want to force the
           routine to read the VMEC file from disk instead of trying to reuse
           data stored in memory.  Otherwise set to zero. Default = 0

     Outputs:
       cLabPos - The cartesian laboratory coordinates at the approximate input location
    """

    out = Flux2Polar(sPos, VMEC_FilePath, ForceReadVMEC=ForceReadVMEC, kMinRho=kMinRho, jacout=jacout, verbose=verbose)

    if jacout:
        RR, fi, ZZ = tuple(out[0][:,0], out[0][:,1], out[0][:,2])
    else:
        RR, fi, ZZ = tuple(out[:,0], out[:,1], out[:,2])
    # end if

    XX, YY, ZZ = _ut.pol2cart(RR, fi, ZZ)
    cLabPos = _np.hstack((XX,YY,ZZ))
    if jacout:
        # Convert to configuration space coordinates from polar coordinates
        return cLabPos, _np.nan # TODO!:
    else:
        return cLabPos
    # end if
# end def


def Flux2Polar(sPos, VMEC_FilePath, ForceReadVMEC=False, kMinRho=2.0e-9, jacout=False, verbose=True):
    """
     Description:
       This function finds the coordinates of a point in polar lab coordinates
       that corresponds with a point in VMEC flux coordinates.

     Inputs:
       sPos - A 2D array describing the points that you want to have
           converted from VMEC coordinates.  The first index is the point
           number, and the second identifies the coordinate.  For example:
          (rho1, th1, fi1 ; rho2, th2, fi2;...).
              rho ~ normalized poloidal flux (sqrt(S))
              th ~ pseudo-poloidal angle [rad]
              fi ~ pseudo-toroidal angle [rad]

       VMEC_FilePath - The location and name of the VMEC output file
           that you want to use for the transformation (eithet netcdf or txt).

       ForceReadVMEC (optional) - Set this to one if you want to force the
           routine to read the VMEC file from disk instead of trying to reuse
           data stored in memory.  Otherwise set to zero. Default = 0

     Outputs:
       pLabPos - The polar laboratory coordinates at the approximate input location
    """
    # ===================================================================== #
    # Fetch VMEC data

    # We'll store the relevant VMEC data in a persistent variables to improve
    # speed for the case where this function is called many times for the same
    # VMEC data.  Unless the user asks for a new file to be used as input, we'll
    # continue to use the persistent data.  Otherwise, load it from a file.
    global VMEC_Data
    global VMEC_DerivedQuant
    global VMEC_DataSource

    # note this test should be for whether the var field is present in VMEC_DataSource
    if ForceReadVMEC or (VMEC_Data is None) or (VMEC_DerivedQuant is None) or (VMEC_DataSource is None) \
        or (VMEC_DataSource != VMEC_FilePath):
        if 1:
            VMEC_Data =  extract_data(VMEC_FilePath, ForceReadVMEC=ForceReadVMEC, verbose=verbose)
            VMEC_DerivedQuant = spline_data(VMEC_Data, ForceReadVMEC=True, verbose=verbose)
        # end if
    # end
    VMEC_DataSource = VMEC_FilePath

    # ===================================================================== #

    sPos = _np.atleast_2d(sPos)
    nrho = len(sPos[:,0])

    nargout = 2
    RR = _np.zeros((nrho,), _np.float64)
    ZZ = _np.zeros_like(RR)

    if jacout:
        nargout=8
        dRdRho, dRdTheta, dZdRho, dZdTheta, dRdZi, dZdZi \
            = tuple([_np.zeros_like(RR) for _ in range(6)])
    # end if

    for ii in range(nrho):
        out = GetLabCoordsFromVMEC(sPos[ii,0], sPos[ii,1], sPos[ii,2],
                    VMEC_Data, VMEC_DerivedQuant, kMinRho=kMinRho, nargout=nargout)

        if jacout:
            RR[ii], ZZ[ii], dRdRho[ii], dRdTheta[ii], dZdRho[ii], dZdTheta[ii], \
                dRdZi[ii], dZdZi[ii] = out
        else:
            RR[ii], ZZ[ii] = out
        # end if
    #end for

    pLabPos = _np.hstack((RR, sPos[:,2], ZZ)) # phi is not the cylindrical angle!!!!
    if jacout:
        return pLabPos, _np.hstack((dRdRho, dRdTheta, dZdRho, dZdTheta, dRdZi, dZdZi))
    return pLabPos
#end def Flux2Polar

# ========================================================================= #
# ------------------------------------------------------------------------- #
# ========================================================================= #
# FUNCTIONS THAT USE ONLY LOCAL VARIABLES BELOW HERE


def FindVMEC_Coords(pLabPos, VMEC_Data, VMEC_DerivedQuant, PosTol=1e-3, RhoGuess=0.5, ThetaGuess=1.0, verbose=True,
                    kMaxIters=150, kMinJacDet=1.0e-16, kMaxItersPerThetaStep=50, kMaxInitialThetaStep=0.1):
                    # kMaxIters=150, kMinJacDet=1.0e-16, kMaxItersPerThetaStep=30, kMaxInitialThetaStep=0.5):
    """
     Description:
       This routine takes a single point in lab space (polar coordinates),
       and uses Newton's method to find approximate the coresponding point
       in VMEC coordinates.  It is intended to be driven by the routine
       CalcFluxCart_VMEC, which handles the reading of the VMEC file and the
       splining of the relevant data contained within.

     Inputs:
           pLabPos - The position of a single pointlab space that is to be
               converted from polar coordinates to VMEC space.
               Ordered [R, Phi, Z].  Units are meters.
           VMEC_Data - This structure contains some raw data straight out of
               the VMEC results file.  We are using the following fields
                   .xm - Poloidal mode numbers
                   .xn - Toroidal mode numbers
           VMEC_DerivedQuant -
                   .rmnc_norm_spline - The "pp" form spline fitting of the
                       normalized spectral coefficients for the R coordinate:
                       \tilde{R_{mn}} in Attenberger's notation
                   .zmns_norm_spline - The "pp" form spline fitting of the
                       normalized spectral coefficients for the Z coordinate:
                       \tilde{Z_{mn}} in Attenberger's notation
                   .rmnc_norm_spline_derivs - The "pp" form spline fitting of the
                       radial derivative of the normalized spectral
                       coefficients for the R coordinate:
                       \frac{\partial}{\partial Rho} \tilde{R{mn}} in
                       Attenberger's notation
                   .zmns_norm_spline_derivs - The "pp" form spline fitting of the
                       radial derivative of the normalized spectral
                       coefficients for the Z coordinate:
                       \frac{\partial}{\partial Rho} \tilde{Z{mn}} in
                       Attenberger's notation
           PosTol - This value, in meters, sets when the routine will stop
               looking for a better fit in Rho, Theta space.  If the
           RhoGuess - A starting guess as to what the Rho value might be.
               This isn't very important.
           ThetaGuess - A starting guess as to what Theta might be to improve
               convergence times.  Also not very important.
           verbose - Boolean for printing to screen

            Because we expect the first derivative of either lab space coordinate to
            change sign for a theta step of pi/2, we can't accurately take a step
            on that scale.  However, capping the maximum step size too low will cause
            the code to take longer.  We attempt to strike a balance by trying to
            detect the failure to converge and adjusting the maximum theta step size.
                kMaxIters - Maximum number of iterations
                kMinJacDet - Minimum value of the Jacobian determinant
                kMaxItersPerThetaStep - Maximum number of iterations at a single theta value
                kMaxInitialThetaStep - Maximum initial step size in theta

     Outputs:
           Rho - converged estimate of the normalized effective radius (sqrt tor. flux / tor. flux @ LCFS)
           Theta - converged estimate of the pseudo-poloidal angle in flux-coordinates
    """
    Rho_tp = _np.copy(RhoGuess)
    Theta_tp = _np.copy(ThetaGuess)
    # Phi_tp = _np.copy(pLabPos[1])
    IterNum = 0
    Converged = False
    MaxThetaStep = _np.copy(kMaxInitialThetaStep)

    Dist_BestMatch = 0.0
    DistR_BestMatch = 0.0
    DistZ_BestMatch = 0.0
    JacDet_BestMatch = 0.0
    dR_dRho_BestMatch = 0.0
    dZ_dRho_BestMatch = 0.0
    dR_dTheta_BestMatch = 0.0
    dZ_dTheta_BestMatch = 0.0
    Theta_BestMatch = 0.0
    Rho_BestMatch = 0.0
    while (IterNum<kMaxIters) and (Converged is False):
        IterNum += 1

        # If we are having trouble converging, try a smaller maximum step in
        # theta.  This seems to help a lot, but is rarely needed.
        if ( IterNum % kMaxItersPerThetaStep == 0):
            MaxThetaStep /= 2.0
        # end if

        # Convert the \Rho, and \Theta for the current test point in to R, Z
        # using the splined coefficients calculated earlier.  This routine also
        # calculates components of the Jacobian matrix at the test point.
        R_tp, Z_tp, dR_dRho_tp, dR_dTheta_tp, dZ_dRho_tp, dZ_dTheta_tp = \
            GetLabCoordsFromVMEC(Rho_tp, Theta_tp, _np.copy(pLabPos[1]), VMEC_Data, VMEC_DerivedQuant)

        #How far off was our guess
        DistR_tp = _np.copy(pLabPos[0]) - R_tp  # in major radius
        DistZ_tp = _np.copy(pLabPos[2]) - Z_tp  # vertically
        Dist_tp = _np.sqrt(DistR_tp*DistR_tp + DistZ_tp*DistZ_tp) # linear distance to correct point

        #Calculate the Jacobian determinant at the test point
        JacDet_tp = dR_dTheta_tp*dZ_dRho_tp-dR_dRho_tp*dZ_dTheta_tp

        if (_np.abs(JacDet_tp)<kMinJacDet) and verbose:
            print('Find_VMEC_Coords: Jacobian determinant approximately zero. R=%6.3e, Z=%6.3e'%(R_tp, Z_tp))
        # end if

        # As mentioned around eq. 5a in the Attenberger paper, it is
        # possible (but rare) for our new test point to be worse than the previous
        # one.  If so, we can try successively smaller steps by cutting the
        # step size in half until we get an improved answer.  We will improve
        # the approximation of the Jacobian determinant using the information
        # gained by evaluating it at the over-step (the choice of weighting
        # appears to be arbitrary). In the notation of the paper, the k-th
        # iteration is the "Best Match" in the program and the k+1-th is the
        # "Test Point".  We still need a maximum step size in the theta
        # direction here, as in the normal step.
        if (IterNum!=1) and (Dist_tp>Dist_BestMatch):
          StepReduction = 1.0

          # Keep halving until you beat your last guest
          while (Dist_tp>Dist_BestMatch) and (IterNum<kMaxIters):
            StepReduction *= 2.0
            ApproxJacDet = 0.75*JacDet_BestMatch + 0.25*JacDet_tp

            ThetaStep = ( (dZ_dRho_BestMatch*DistR_BestMatch-dR_dRho_BestMatch*DistZ_BestMatch)
                           / (StepReduction*ApproxJacDet) )
            if (_np.abs(ThetaStep)>MaxThetaStep):
                ThetaStep = MaxThetaStep*_np.sign(ThetaStep)
            # end if
            Theta_tp = Theta_BestMatch + ThetaStep

            Rho_tp = (Rho_BestMatch +
                        ( (dR_dTheta_BestMatch*DistZ_BestMatch-dZ_dTheta_BestMatch*DistR_BestMatch)
                          / (StepReduction*ApproxJacDet) ) )

            if (Rho_tp<0):
                Rho_tp = -Rho_tp
                Theta_tp = Theta_BestMatch + _np.pi - ThetaStep
            # end if

            # Recalculate the test position coordinates in flux-coordinates
            R_tp, Z_tp, dR_dRho_tp, dR_dTheta_tp, dZ_dRho_tp, dZ_dTheta_tp = \
                GetLabCoordsFromVMEC(Rho_tp, Theta_tp, pLabPos[1], VMEC_Data, VMEC_DerivedQuant)

            # How far off was our test point
            DistR_tp = pLabPos[0] - R_tp   # in major radius
            DistZ_tp = pLabPos[2] - Z_tp   # vertically
            Dist_tp = _np.sqrt(DistR_tp*DistR_tp + DistZ_tp*DistZ_tp)  # linear distance

            # Recalculate the Jacobian determinant at the test point
            JacDet_tp = dR_dTheta_tp*dZ_dRho_tp-dR_dRho_tp*dZ_dTheta_tp

            IterNum += 1

            if (_np.abs(JacDet_tp)<kMinJacDet) and verbose:
                print('Find_VMEC_Coords: Jacobian determinant approximately zero. R=%6.3e, Z= %6.3e'%(R_tp, Z_tp))
            # end if jacobian determinant less than tolerance
          # end while distance greater than best match and iteration number below maximum
        # end if distance greater than best match and iteration number greater than 1

        # If necessary, pick the next test point using the typical Newton's
        # method. Otherwise, return the approximate location in VMEC coordinates.
        if (Dist_tp<PosTol):
            Converged = True
            Rho = _np.copy(Rho_tp)
            Theta = Theta_tp%(2.0*_np.pi)
            if (Theta<0):
                Theta += 2.0*_np.pi
            #end
        else:
            # This is the typical case: We have just identified an improved
            # estimate of Rho and Theta and about to take another normal step.
            # We only typically need to know the coordinates of the best match
            # so far, but the half-step case needs the Jacobian elements and the
            # distance that the test point was from its target.
            Rho_BestMatch = _np.copy(Rho_tp)
            Theta_BestMatch = _np.copy(Theta_tp)
            Dist_BestMatch = _np.copy(Dist_tp)
            DistR_BestMatch = _np.copy(DistR_tp)
            DistZ_BestMatch = _np.copy(DistZ_tp)
            dR_dRho_BestMatch = _np.copy(dR_dRho_tp)
            dR_dTheta_BestMatch = _np.copy(dR_dTheta_tp)
            dZ_dRho_BestMatch = _np.copy(dZ_dRho_tp)
            dZ_dTheta_BestMatch = _np.copy(dZ_dTheta_tp)
            JacDet_BestMatch = _np.copy(JacDet_tp)

            # The normal step, as taken from eq. 3,4 in the paper.  I'm not sure
            # why exactly they formulated the problem this way, but it works out
            # to be identical to what you would expect from Newton's method.
            # (Note that you have to use the fact that \Rho and \Theta are
            # orthagonal to get the algebra to work out).  Without a maximum
            # step size in the theta direction, convergence is very slow or
            # nonexistant.
            Rho_tp = Rho_tp + (dR_dTheta_tp*DistZ_tp-dZ_dTheta_tp*DistR_tp)/JacDet_tp
            ThetaStep = (dZ_dRho_tp*DistR_tp-dR_dRho_tp*DistZ_tp)/JacDet_tp
            if (_np.abs(ThetaStep)>MaxThetaStep):
                ThetaStep = MaxThetaStep*_np.sign(ThetaStep)
            #end
            Theta_tp = Theta_tp + ThetaStep

            # If the solver tries to look at the "negative rho", we must correct
            # it and the theta step that goes along with it.
            if (Rho_tp < 0):
                Rho_tp = -Rho_tp
                Theta_tp = Theta_BestMatch + _np.pi - ThetaStep
            #end
        #end
    #end

    if (IterNum == kMaxIters) or ('Rho' not in locals()):
        print('Attempt to find flux coordinate failed to converge:\n Phi=%6.4e\n R=%6.4e\n Z=%6.4e'
               %(pLabPos[1],pLabPos[0],pLabPos[2]))

        print('Best guess:\n Phi=%6.4e\n R=%6.4e\n Z=%6.4e'%(pLabPos[1],R_tp,Z_tp))
        Rho = _np.nan
        Theta = _np.nan

    # end if
    return Rho, Theta
#    [ttemp,rtemp,ztemp] = cart2pol( pLabPos[0],pLabPos[1],pLabPos[2] )
#    VMEC_DerivedQuant.tor = np.concatenate( (VMEC_DerivedQuant.tor, ttemp), axis=0 )
#    VMEC_DerivedQuant.rr  = np.concatenate( (VMEC_DerivedQuant.rr,  rtemp), axis=0 )
#    VMEC_DerivedQuant.zz  = np.concatenate( (VMEC_DerivedQuant.zz,  ztemp), axis=0 )
#    VMEC_DerivedQuant.rho = np.concatenate( (VMEC_DerivedQuant.rho, Rho_BestMatch), axis=0 )
#    VMEC_DerivedQuant.th  = np.concatenate( (VMEC_DerivedQuant.th,  Theta_BestMatch), axis=0 )
#    VMEC_DerivedQuant.drdth = np.concatenate( (VMEC_DerivedQuant.drdth, dR_dTheta_BestMatch), axis=0 )
#    VMEC_DerivedQuant.dzdth = np.concatenate( (VMEC_DerivedQuant.dzdth, dZ_dTheta_BestMatch), axis=0 )
#    VMEC_DerivedQuant.drdrho = np.concatenate( (VMEC_DerivedQuant.drdth, dR_dRho_BestMatch), axis=0 )
#    VMEC_DerivedQuant.dzdrho = np.concatenate( (VMEC_DerivedQuant.dzdth, dZ_dRho_BestMatch), axis=0 )
#    VMEC_DerivedQuant.jacdet = np.concatenate( (VMEC_DerivedQuant.jacdet, JacDet_BestMatch), axis=0 )


# ========================================================================= #


def GetLabCoordsFromVMEC(Rho, Theta, Phi, VMEC_Data, VMEC_DerivedQuant, kMinRho=2.0e-9, nargout=6):
    """
     File Name: GetLabCoordsFromVMEC.m

       Original version in MATLAB written by: Christopher A. Clark Creation Date: 06/2012
       Ported to PYTHON by Gavin M. Weir, Creation Date: 03/2016

     Description:
       This routine takes a single point in VMEC space, and finds the
       corresponding point in polar lab coordinates along with elements of the
       local Jacobian matrix of the transformation between the two systems.

     Inputs:
       Rho - VMEC normalized flux coordinate (sqrt(S))
       Theta - VMEC poloidal angle           (theta)
       Phi - VMEC (and lab) toroidal angle  (zi)
       VMEC_Data - A structure of parsed VMEC data, which must include:
                .xm - The poloidal mode numbers associated with the spectral
                      components of the lab space coordinates.
                .xn - The toroidal mode numbers associated with the spectral
                      components of the lab space coordinates.
       VMEC_DerivedQuant -
               .rmnc_norm_spline - The "pp" form spline fitting of the
                   normalized spectral coefficients for the R coordinate:
                   \tilde{R_{mn}} in Attenberger's notation
               .zmns_norm_spline - The "pp" form spline fitting of the
                   normalized spectral coefficients for the Z coordinate:
                   \tilde{Z_{mn}} in Attenberger's notation
               .rmnc_norm_spline_derivs - The "pp" form spline fitting of the
                   radial derivative of the normalized spectral
                   coefficients for the R coordinate:
                   \frac{\partial}{\partial Rho} \tilde{R{mn}} in
                   Attenberger's notation
               .zmns_norm_spline_derivs - The "pp" form spline fitting of the
                   radial derivative of the normalized spectral
                   coefficients for the Z coordinate:
                   \frac{\partial}{\partial Rho} \tilde{Z{mn}} in
                   Attenberger's notation
     kMinRho - We must set a minimum rho for the purposes of normalizing the spectral components.

     Outputs:
       r - The major radial coordinate of the input rho, theta, phi
       z - The vertical coordinate of the input rho, theta, phi

       dR_dRho - The partial derivative of the major radial coordinate with
           respect to VMEC Rho. (One term in the Jacobian Matrix)
       dZ_dRho - The partial derivative of the vertical coordinate with
           respect to VMEC Rho. (One term in the Jacobian Matrix)
       dR_dTheta - The partial derivative of the major radial coordinate with
           respect to VMEC poloidal angle. (One term of the Jacobian Matrix)
       dZ_dTheta - The partial derivative of the vertical coordinate with
           respect to VMEC poloidal angle. (One term of the Jacobian Matrix)
       if requested:
           dR_dZi - The partial derivative of the major radial coordinate with
               respect to VMEC toroidal angle.
           dZ_dZi - The partial derivative of the vertical coordinate with
               respect to VMEC toroidal angle.
    """
    # Get the normalized spectral coefficients and the radial derivative of the
    # normalized spectral coefficients from the splines evaluated at the
    # specific radial location of interest.  If the point of interest is
    # outside of the last closed flux surface, we will hold all of the spectral
    # terms except for m=1,n=0 at their rho=1 value.  That 1,0 term is extended
    # proportionally to rho, but that is handled in the de-normalization
    # section, below.
    #  .... note this is kind of like asymptoting things outside the LCFS to an ellipse
    if (Rho>1):
        RhoSpline = 1 # evaluate at 1 if rho is larger than 1
    else:
        RhoSpline = _np.copy(Rho)
    # end if rho>1

    rmnc_norm = _np.zeros((len(VMEC_DerivedQuant.rmnc_norm_spline),len(_np.atleast_1d(RhoSpline))), dtype=_np.float64)
    zmns_norm = _np.zeros_like(rmnc_norm)

    RhoDeriv_rmnc_norm = _np.zeros_like(rmnc_norm)
    RhoDeriv_zmns_norm = _np.zeros_like(rmnc_norm)

    for ii in range(len(VMEC_DerivedQuant.rmnc_norm_spline)):
        rmnc_norm[ii,:] = VMEC_DerivedQuant.rmnc_norm_spline[ii](RhoSpline)
        zmns_norm[ii,:] = VMEC_DerivedQuant.zmns_norm_spline[ii](RhoSpline)

        RhoDeriv_rmnc_norm[ii,:] = VMEC_DerivedQuant.rmnc_norm_spline[ii](RhoSpline, nu=1)
        RhoDeriv_zmns_norm[ii,:] = VMEC_DerivedQuant.zmns_norm_spline[ii](RhoSpline, nu=1)
#        RhoDeriv_rmnc_norm[ii,:] = VMEC_DerivedQuant.rmnc_norm_spline_derivs[ii](RhoSpline)
#        RhoDeriv_zmns_norm[ii,:] = VMEC_DerivedQuant.zmns_norm_spline_derivs[ii](RhoSpline)
    # end for

    # Find the cos(m\theta - n \phi) and sin(m\theta - n\phi) terms only once,
    # to save processor time.
    ModeAngle = VMEC_Data.xm*Theta-VMEC_Data.xn*Phi
    CosModeAngle = _np.cos(ModeAngle)
    SinModeAngle = _np.sin(ModeAngle)

    # We'll need to undo the normalization of the spectral magnitudes that was
    # performed prior to the splining. As described in the Attenberger paper,
    # this must be approximated by some finite, but small value as \Rho goes to
    # zero. To extend the VMEC coordinates outside of the LCFS, we grow the m=1,
    # n=0 modes linearly with \Rho, and hold the other terms constant at their
    # value evaluated at \Rho=1
    if (Rho>1):   # if rho is outside the LCFS
        SpectCoeffNorm = _np.ones((_np.size(VMEC_Data.xm,axis=0),1), _np.float64)
        SpectCoeffNorm[(VMEC_Data.xm == 1)*(VMEC_Data.xn == 0)] = _np.copy(Rho)

        dSpectCoeffNorm_dRho = _np.zeros((_np.size(VMEC_Data.xm,axis=0),1), _np.float64)
        dSpectCoeffNorm_dRho[(VMEC_Data.xm == 1)*(VMEC_Data.xn == 0)] = 1.0

    elif (Rho<kMinRho):   # if rho is near the axis
        SpectCoeffNorm = kMinRho**_np.abs(VMEC_Data.xm)
        dSpectCoeffNorm_dRho = _np.abs(VMEC_Data.xm) * kMinRho**(_np.abs(VMEC_Data.xm-1))
    else: # everywhere else in the plasma
        SpectCoeffNorm = Rho**_np.abs(VMEC_Data.xm)
        dSpectCoeffNorm_dRho = _np.abs(VMEC_Data.xm)*(Rho**(_np.abs(VMEC_Data.xm-1)))
    # end if rho>1 or rho small

    # Sum up all of the contributions of each mode to get both the positions in
    # major radius and vertical direction (r,z).
    r = _np.sum(rmnc_norm*CosModeAngle*SpectCoeffNorm)
    z = _np.sum(zmns_norm*SinModeAngle*SpectCoeffNorm)

    if nargout == 2:
        return r, z

    # The partial derivative in the theta direction is easily extracted from the
    # above equation, since only the cos (sin) term contains a theta dependence
    dR_dTheta = -1.0*_np.sum(VMEC_Data.xm*rmnc_norm*SinModeAngle*SpectCoeffNorm)
    dZ_dTheta = _np.sum(VMEC_Data.xm*zmns_norm*CosModeAngle*SpectCoeffNorm)

    # Properly renormalizing the _radial_ derivatives of r and z is more
    # complicated than is immediately obvious because we only have the radial
    # derivatives of the _normalized_ spectral coefficients.  Application of the
    # chain rule is necessary.  In the notation of Attenberger:
    # \frac{\partial}{\partial \Rho} R_mn = m \Rho^{m-1} \tilde{R_mn} + \Rho^m
    # \frac{partial}{\partial \Rho} \tilde{R_mn}
    dR_dRho = _np.sum( (RhoDeriv_rmnc_norm*SpectCoeffNorm + rmnc_norm*dSpectCoeffNorm_dRho)*CosModeAngle )
    dZ_dRho = _np.sum( (RhoDeriv_zmns_norm*SpectCoeffNorm + zmns_norm*dSpectCoeffNorm_dRho)*SinModeAngle )

    if nargout != 6:
        # if requested, also calculate the derivatives wrt to toroidal angle
        dR_dZi = _np.sum( VMEC_Data.xn*rmnc_norm*SinModeAngle*SpectCoeffNorm )
        dZ_dZi = -1.0*_np.sum( VMEC_Data.xn*zmns_norm*CosModeAngle*SpectCoeffNorm )
        return r, z, dR_dRho, dR_dTheta, dZ_dRho, dZ_dTheta, dR_dZi, dZ_dZi
    # end if
    return r, z, dR_dRho, dR_dTheta, dZ_dRho, dZ_dTheta


# ========================================================================= #
# ========================================================================= #


def extract_data(VMEC_FilePath, ForceReadVMEC=False, verbose=True):
    global VMEC_Data
    if VMEC_Data is None or ForceReadVMEC:
        VMEC_Data = _ut.Struct() # Instantiate an empty class of type structure

        VMEC_Data = read_vmec(VMEC_FilePath)

        if VMEC_FilePath.lower().find('.nc')>-1:
            VMEC_Data.ns = VMEC_Data.ns
            VMEC_Data.xm = VMEC_Data.xm   # (mn_mode) polidal mode numbers
            VMEC_Data.xn = VMEC_Data.xn   # (mn_mode) toridal mode numbers
            VMEC_Data.zmns = VMEC_Data.zmns.T # (radius mn_mode ) sinmn component of cylindrical Z, full mesh in [m]
            VMEC_Data.rmnc = VMEC_Data.rmnc.T # (radius mn_mode ) cosmn component of cylindrical R, full mesh in [m]

            # defines xn with the opposite sign as is usual in VMEC
            VMEC_Data.xn = -VMEC_Data.xn

            # VMEC_Data.xm = VMEC_Data.xm.T
            # VMEC_Data.xn = -VMEC_Data.xn.T
            # VMEC_Data.ns = VMEC_Data.ns.T
            # VMEC_Data.zmns = VMEC_Data.zmns.T
            # VMEC_Data.rmnc = VMEC_Data.rmnc.T
            # # VMEC_Data.xn = -VMEC_Data.xn
        #end
        VMEC_Data.loaded = True
        if verbose:
            print('Done Reading VMEC File')
        # end if
    # end if
    return VMEC_Data

# ========================================================================= #


def spline_data(VMEC_Data, ForceReadVMEC=True, verbose=True):
    global VMEC_DerivedQuant
    if VMEC_DerivedQuant is None or ForceReadVMEC:
        if verbose:
            print('Normalizing.')
        # Set up radial grids.  VMEC produces data that is equally spaced in 'S'
        # (normalized total flux) space.  We often want it in normalized minor
        # radius (rho), though.
        VMEC_DerivedQuant = _ut.Struct()
#        VMEC_DerivedQuant.GridS = ((1:VMEC_Data.ns)-1)/(VMEC_Data.ns-1)
        VMEC_DerivedQuant.GridS = (_np.asarray(range(1,VMEC_Data.ns+1))-1)/(VMEC_Data.ns-1.0)
        VMEC_DerivedQuant.GridRho = _np.sqrt(VMEC_DerivedQuant.GridS)

        # Normalizing the spectral coefficients to rho^m helps the splining get
        # the radial dependence correct near the axis, as described in the paper.
        # See equations 7a and 7b and the paragraphs directly above those
        # equations for more details. Notationally, rmnc_norm corresponds with
        # \tilde{R}_mn, for example. The values right at the axis are linearly
        # extrapolated: \tidle{R_(mn)}(\rho = 0) = \tilde{R_(mn)}(\rho = \Delta \rho) -
        # \Delta \rho \frac{\tilde{R_(mn)(2\Delta \rho) - \tilde{R_(mn)(\Delta \rho}{\Delta \rho}}
        TempExp = _np.abs(VMEC_Data.xm)*_np.ones((1,len(VMEC_DerivedQuant.GridRho)), float)
        TempBase = _np.ones((len(VMEC_Data.xm), 1), float)*_np.atleast_2d(VMEC_DerivedQuant.GridRho)
        NormalizationFactor = TempBase**TempExp

#        TempExp = _np.dot( _np.abs(VMEC_Data.xm), _np.ones((1,len(VMEC_DerivedQuant.GridRho)), float) )
#        TempBase = _np.dot( _np.ones((len(VMEC_Data.xm), 1), float), _np.atleast_2d(VMEC_DerivedQuant.GridRho))
#        NormalizationFactor = TempBase**TempExp

        rmnc_norm = _np.zeros_like(VMEC_Data.rmnc)
#        rmnc_norm[:,:] = _np.divide(VMEC_Data.rmnc[:, :], NormalizationFactor[:, :], where=NormalizationFactor[:, :]!=0)
        rmnc_norm[:,1:] = _np.divide(VMEC_Data.rmnc[:, 1:], NormalizationFactor[:, 1:], where=NormalizationFactor[:, 1:]!=0)
        rmnc_norm[:,0] = 2.0*rmnc_norm[:,1] - rmnc_norm[:,2]

        zmns_norm = _np.zeros_like(VMEC_Data.zmns)
#        zmns_norm[:,:] = _np.divide(VMEC_Data.zmns[:, :], NormalizationFactor[:, :], where=NormalizationFactor[:, :]!=0)
        zmns_norm[:,1:] = _np.divide(VMEC_Data.zmns[:, 1:], NormalizationFactor[:, 1:], where=NormalizationFactor[:, 1:]!=0)
        zmns_norm[:,0] = 2.0*zmns_norm[:,1] - zmns_norm[:,2]

        if verbose:
            print('Done normalizing.')
            print('Generating splines for spectral cooefficients.')
        # end if

        # We will capture the radial dependence of the spectral normalized
        # components of r and z using splines fit to 'double-sided' profiles.
        # This approximates the radial partial derivatives of the normalized
        # spectral components, too (but not directly the radial derivatives of r,z).
        Spline_Rhos = _ut.cylsym_odd(VMEC_DerivedQuant.GridRho, exclude_axis=True)
        rmnc_vals_to_spline = _ut.cylsym_even(rmnc_norm, axis=1, exclude_axis=True)
        zmns_vals_to_spline = _ut.cylsym_even(zmns_norm, axis=1, exclude_axis=True)
        # Spline_Rhos = np.concatenate((_np.flipud( -VMEC_DerivedQuant.GridRho[1:]), VMEC_DerivedQuant.GridRho), axis=1)
        # rmnc_vals_to_spline = np.concatenate( ( np.fliplr( rmnc_norm[ :, 1:-1 ] ) , rmnc_norm ), axis=1)
        # zmns_vals_to_spline = np.concatenate( ( np.fliplr( zmns_norm[ :, 1:-1 ] ) , zmns_norm ), axis=1)

        VMEC_DerivedQuant.rmnc_norm_spline = []
#        VMEC_DerivedQuant.rmnc_norm_spline_derivs = []
        VMEC_DerivedQuant.zmns_norm_spline = []
#        VMEC_DerivedQuant.zmns_norm_spline_derivs = []

        if DEBUG:
            nharm = _np.size(rmnc_vals_to_spline, axis=1)
            maxplots = 5
            iplot = nharm // (maxplots+1)
        # end if

        # Do not specify zero slope at convex hull
        splargs = {'end1':int(0), 'end2':int(0), 'slope1':int(0), 'slope2':int(0)}

        # Specify zero slope at convex hull
#        splargs = {'end1':int(1), 'end2':int(1), 'slope1':int(0), 'slope2':int(0)}

        for ii in range(rmnc_vals_to_spline.shape[0]):
            SplObj = _ut.Spline(Spline_Rhos, rmnc_vals_to_spline[ii,:], **splargs).spline() # homemade
#            SplObj = _dsi.UnivariateSpline(Spline_Rhos, rmnc_vals_to_spline[ii,:])  # dsi

            if DEBUG and (ii % (iplot//2) == 0):
                SplObj.show()
            # end def

            VMEC_DerivedQuant.rmnc_norm_spline.append( SplObj )
#            VMEC_DerivedQuant.rmnc_norm_spline_derivs.append( SplObj.derivative(1) ) # First derivative

            # ====== #

            SplObj = _ut.Spline(Spline_Rhos,zmns_vals_to_spline[ii,:], **splargs).spline() #  homemade
#            SplObj = _dsi.UnivariateSpline(Spline_Rhos,zmns_vals_to_spline[ii,:])

            if DEBUG and (ii % (iplot//2) == 0):
                SplObj.show()
            # end def

            VMEC_DerivedQuant.zmns_norm_spline.append( SplObj )
#            VMEC_DerivedQuant.zmns_norm_spline_derivs.append( SplObj.derivative(1) )

        # end for

        if verbose:
            print('Done splining the spectral coefficients')
        # end if
    # end if
    return VMEC_DerivedQuant

# ========================================================================= #
# ------------------------------------------------------------------------- #
# ========================================================================= #


class PYVMEC(_ut.Struct):
    loaded = False
    verbose = True
    def __init__(self, fname=None, **kwargs):
        kwargs['fname'] = fname
        self.__dict__.update(kwargs)

        if hasattr(self, 'input_extension') and fname is None:
            tstname = _os.path.join(_os.path.abspath(_os.path.curdir), self.input_extension)
            if _os.path.exists(tstname):
                self.fname = tstname
            # end if
        # end if

        # Check if the data was already loaded to this file through a dictionary
        if _os.path.exists(self.fname):
            self.read(self.fname, False, self.verbose)
        if hasattr(self, 'rmnc'):
            self.loaded = True
#            self.roa = _np.linspace(1e-5, 1, 200)
#            self.getB00()
#            self.setB00(Bfactor*self.B00)
#            self.getB00()
#            self.getamin()
#            self.getgradrho()
#            self.getVol()
#            self.getiota()
        # end if
    # end def __init__

    def read(self, VMEC_FilePath, ForceReadVMEC=True, verbose=False):
#        self.__dict__.update(extract_data(VMEC_FilePath, ForceReadVMEC, verbose))
        self.VMEC_Data = extract_data(VMEC_FilePath, ForceReadVMEC, verbose)
        self.VMEC_DerivedQuant = spline_data(self.VMEC_Data, True, verbose)
        self.VMEC_FilePath = VMEC_FilePath
        self.loaded = True
    # end def

    def cart2flux(self, cLabPos, verbose=True, **kwargs):
        # ===================================================================== #
        # Use Newton's method to find the values of Rho and Theta that match our
        # position in real space.
        PosTol = kwargs.setdefault('PosTol', 1e-3)

        #Convert the user's cartesian data into polar form indexed as r,phi,z
        pLabPos = _np.zeros_like(cLabPos)
        pLabPos[:, 1], pLabPos[:, 0], pLabPos[:, 2] = _ut.cart2pol(cLabPos[:, 0], cLabPos[:, 1], cLabPos[:, 2] )
        Phi = pLabPos[:, 1]

        nrho = len(Phi)
        Rho = _np.zeros((nrho,), _np.float64)
        Theta = _np.zeros_like(Rho)
        for ii in range(nrho):
            Rho[ii], Theta[ii] = FindVMEC_Coords(pLabPos[ii, :], self.VMEC_Data, self.VMEC_DerivedQuant, PosTol=PosTol)
        #end for
        return Rho, Theta
#        return Cart2Flux_VMEC(cLabPos, self.VMEC_FilePath, verbose=verbose, **kwargs)

    def flux2cart(self, sLabPos, **kwargs):
        out = self.flux2cart(sLabPos, self.VMEC_FilePath, **kwargs)

        RR, fi, ZZ = tuple(out[:,0], out[:,1], out[:,2])
        XX, YY, ZZ = _ut.pol2cart(RR, fi, ZZ)
        cLabPos = _np.hstack((XX,YY,ZZ))
        return cLabPos
#        return Flux2Cart(sLabPos, self.VMEC_FilePath, verbose=verbose, **kwargs)

    def flux2polar(self, sLabPos, **kwargs):
        kwargs.setdefault('kMinRho', 2.0e-9)
        kwargs.setdefault('verbose', True)
        return Flux2Polar(sLabPos, self.VMEC_FilePath, **kwargs)

    def __call__(self, cLabPos):
        return self.cart2flux(cLabPos, self.VMEC_FilePath)
    # end def

#    @property
#    def B00(self):
#        """
#         Get the B<SUB>00</SUB> value of the magnetic field on axis.
#         @return the B<SUB>00</SUB> value of the magnetic field on magnetic axis.
#        """
#        if not hasattr(self,"_B00"):
#            self._B00 = _np.copy(self.VMEC_Data.B00)
#        # endif
#        return self._B00
#    @B00.setter
#    def B00(self, value=None):
#        # Here we need to scale all the magnetic fields and associate the changes into the currents as well...
#        self.mconf.MCsetB00(self.MC, value)
#        self._B00 = value
#    @B00.deleter
#    def B00(self):
#        del self._B00

    # ============================================= #

    @staticmethod
    def _RootG(RR, dR_dRho, dR_dTheta, dZ_dRho, dZ_dTheta):
        """
            Return the jacobian for flux surface average integration (flux surface area element)
            ... this is the sqrt(g)
        """
        return RR*(dR_dTheta*dZ_dRho - dR_dRho*dZ_dTheta)

    @staticmethod
    def _dVdrho(M, Rho, VMEC_Data, VMEC_DerivedQuant, kMinRho=2.0e-9):
        """
        Return the incremental volume between \Rho and \Rho + d\Rho
        """
        def _sqrtG(Theta, Phi):
            RR, _, dR_dRho, dR_dTheta, dZ_dRho, dZ_dTheta = \
                GetLabCoordsFromVMEC(Rho, Theta, Phi, VMEC_Data, VMEC_DerivedQuant, kMinRho=kMinRho, nargout=6)
            return PYVMEC._RootG(RR, dR_dRho, dR_dTheta, dZ_dRho, dZ_dTheta)

        # ======== Returns the incremental volume between \Rho and \Rho + d\Rho ==== #
        # simpson's integration over toroidal and poloidal angle (one field period)
        # func = _sqrtG(Theta, Phi)
        dVdrho = M*_ut.aquad.dblquad_asr(_sqrtG, 0.0, 2.0*_np.pi, 0.0, 2.0*_np.pi/M, tol=1e-3)
#        dVdrho = M*_ut.aquad.dblquad_scipy(_sqrtG, 0.0, 2.0*_np.pi, 0.0, 2.0*_np.pi/M, tol=1e-3)
#        dVdrho = M*_ut.aquad.dblquadrature(_sqrtG, 0.0, 2.0*_np.pi, 0.0, 2.0*_np.pi/M, tol=1e-3)

#        # due to speed reasons, let's try to mesh it and calculate
#        ph = _np.linspace(0.0, 2.0*_np.pi, num=10)/M
##        func = lambda x: _sqrtG(x, ph)
##        dVdrho = M*_ut.aquad.quad_asr(func, 0.0, 2.0*_np.pi, tol=1e-3)
#        th = _np.linspace(0.0, 2.0*_np.pi, num=15)
##        func = lambda x: _sqrtG(th, x)
##        dVdrho = M*_ut.aquad.quad_asr(func, 0.0, 2.0*_np.pi/M, tol=1e-3)
#        res = 0.0
#        dphi = ph[1]-ph[0]
#        dth = th[1]-th[0]
#        for _ph in ph:
#            for _th in th:
#                res += _sqrtG(_th, _ph)
#            # end for
#        # end for
#        dVdrho = res*dphi*dth
        return dVdrho

#    @staticmethod
#    def GetdVdrho(M, Rho, VMEC_Data, VMEC_DerivedQuant, kMinRho=2.0e-9):
    def GetdVdrho(self, Rho, **kwargs):
        """
        Return the incremental volume between \Rho and \Rho + d\Rho for a vector of
        input flux-surfaces
        """
        kMinRho = kwargs.setdefault('kMinRho',None)
        if kMinRho is None and hasattr(self, 'kMinRho'):
            kMinRho = self.kMinRho
        else:
            kMinRho = 2.0e-9
        # end if
        self.kMinRho = kMinRho

        Rho = _np.atleast_1d(Rho).copy()
        dVdrho = _np.zeros_like(Rho)
        nrho = len(dVdrho)
        for ii in range(nrho):
            if Rho[ii] == 0:
                Rho[ii] = 1e-5
            dVdrho[ii] = PYVMEC._dVdrho(self.VMEC_Data.nfp, Rho[ii], self.VMEC_Data, self.VMEC_DerivedQuant, kMinRho=self.kMinRho)
        # end for
        self.Vol = _np.cumsum(dVdrho)
        self.dVdrho = dVdrho
        return self.dVdrho, self.Vol

    # ========================================================================= #

    @staticmethod
    def LocalGradRhoSquared(Rho, Theta, Phi, VMEC_Data, VMEC_DerivedQuant, kMinRho=2.0e-9):
        RR, ZZ, dR_dRho, dR_dTheta, dZ_dRho, dZ_dTheta, dR_dZi, dZ_dZi = \
            GetLabCoordsFromVMEC(Rho, Theta, Phi, VMEC_Data, VMEC_DerivedQuant, kMinRho=kMinRho, nargout=8)

        _gr2, sqg = PYVMEC._GradRhoSquared(RR, dR_dRho, dR_dTheta, dZ_dRho, dZ_dTheta, dR_dZi, dZ_dZi)
        return _gr2

    @staticmethod
    def LocalGradRho(Rho, Theta, Phi, VMEC_Data, VMEC_DerivedQuant, kMinRho=2.0e-9):
        RR, ZZ, dR_dRho, dR_dTheta, dZ_dRho, dZ_dTheta, dR_dZi, dZ_dZi = \
            GetLabCoordsFromVMEC(Rho, Theta, Phi, VMEC_Data, VMEC_DerivedQuant, kMinRho=kMinRho, nargout=8)

        _gr2, sqg = PYVMEC._GradRhoSquared(RR, dR_dRho, dR_dTheta, dZ_dRho, dZ_dTheta, dR_dZi, dZ_dZi)
        return _np.sqrt(PYVMEC.LocalGradRhoSquared(Rho, Theta, Phi, VMEC_Data, VMEC_DerivedQuant, kMinRho=kMinRho))


    @staticmethod
    def _GradRhoSquared(RR, dR_dRho, dR_dTheta, dZ_dRho, dZ_dTheta, dR_dZi, dZ_dZi):
        """
            Calculate the gradient of the normalized effective radius at a single-point in space
            ... gradrho(\Rho, \Theta, \Phi)

            ( \nabla rho )^2 = 1/g * [ g_(th,th) * g_(zi,zi) - g_(th,zi)^2 ]
                  = [(R_th^2+Z_th^2)*(R^2+R_zi^2+Z_zi^2) - (R_th*R_zi+Z_th*Z_zi)]
                    / [R^2 * (R_th*Z_rho - R_rho*Z_th)^2]

                sqrt(g) = R * (R_th*Z_rho - R_rho*Z_th)
                dV/drho = int( sqrt(g) dth dzi )
        """
        numerator = (dR_dTheta**2.0+dZ_dTheta**2.0)*(dR_dZi**2.0+RR**2.0+dZ_dZi**2.0)
        numerator -= (dR_dTheta*dR_dZi+dZ_dTheta*dZ_dZi)**2.0

        sqrtG = PYVMEC._RootG(RR, dR_dRho, dR_dTheta, dZ_dRho, dZ_dTheta)
        return numerator/(sqrtG**2.0), sqrtG

    @staticmethod
    def _GradRho(RR, dR_dRho, dR_dTheta, dZ_dRho, dZ_dTheta, dR_dZi, dZ_dZi):
        _gr2, sqg = PYVMEC._GradRhoSquared(RR, dR_dRho, dR_dTheta, dZ_dRho, dZ_dTheta, dR_dZi, dZ_dZi)
        return _np.sqrt(_gr2)*sqg

    @staticmethod
    def _FSAGradRho(M, Rho, VMEC_Data, VMEC_DerivedQuant, kMinRho=2.0e-9):
        """
            Calculate the gradient of the normalized effective radius at a single-point in space
            ... gradrho(\Rho, \Theta, \Phi)

            M - number of field periods
            Rho - sqrt(s), normalized effective radius
        """
        def _Gradrho(Theta, Phi):
            RR, ZZ, dR_dRho, dR_dTheta, dZ_dRho, dZ_dTheta, dR_dZi, dZ_dZi = \
                GetLabCoordsFromVMEC(Rho, Theta, Phi, VMEC_Data, VMEC_DerivedQuant, kMinRho=kMinRho, nargout=8)

            _gr2, sqg = PYVMEC._GradRhoSquared(RR, dR_dRho, dR_dTheta, dZ_dRho, dZ_dTheta, dR_dZi, dZ_dZi)
            return _np.sqrt(_gr2)*sqg

        # simpson's integration over toroidal and poloidal angle (one field period)
        dVdrho = PYVMEC._dVdrho(M, Rho, VMEC_Data, VMEC_DerivedQuant, kMinRho=kMinRho)

        # adaptive integration:
        gradrho = M*_ut.aquad.dblquad_asr(_Gradrho, 0.0, 2.0*_np.pi, 0.0, 2.0*_np.pi/M, tol=1e-3)
#        gradrho = M*_ut.aquad.dblquad_scipy(_Gradrho, 0.0, 2.0*_np.pi, 0.0, 2.0*_np.pi/M, tol=1e-3)
#        gradrho = M*_ut.aquad.dblquadrature(_Gradrho, 0.0, 2.0*_np.pi, 0.0, 2.0*_np.pi/M, tol=1e-3)

#        # due to speed reasons, let's try to mesh it and calculate
#        ph = _np.linspace(0.0, 2.0*_np.pi, num=10)/M
##        func = lambda x: _Gradrho(x, ph)
##        gradrho = M*_ut.aquad.quad_asr(func, 0.0, 2.0*_np.pi, tol=1e-3)
#        th = _np.linspace(0.0, 2.0*_np.pi, num=15)
##        func = lambda x: _Gradrho(th, x)
##        gradrho = M*_ut.aquad.quad_asr(func, 0.0, 2.0*_np.pi/M, tol=1e-3)
#        res = 0.0
#        dphi = ph[1]-ph[0]
#        dth = th[1]-th[0]
#        for _ph in ph:
#            for _th in th:
#                res += _Gradrho(_th, _ph)
#            # end for
#        # end for
#        gradrho = res*dphi*dth
        return gradrho/dVdrho, dVdrho

#    @staticmethod
#    def FSAGradRho(M, Rho, VMEC_Data, VMEC_DerivedQuant, kMinRho=2.0e-9):
    def FSAGradRho(self, Rho, **kwargs):
        """
            Calculate the flux-surface averaged gradient of the normalized
            effective radius on a single surface
            ... gradrho(\Rho)
        """
        kMinRho = kwargs.setdefault('kMinRho',None)
        if kMinRho is None and hasattr(self, 'kMinRho'):
            kMinRho = self.kMinRho
        else:
            kMinRho = 2.0e-9
        # end if
        self.kMinRho = kMinRho

        Rho = _np.atleast_1d(Rho).copy()
        gradrho = _np.zeros_like(Rho)
        dVoldrho = _np.zeros_like(Rho)
        nrho = len(gradrho)
        for ii in range(nrho):
            if Rho[ii] == 0:
                Rho[ii] = 1e-5
            gradrho[ii], dVoldrho[ii] = PYVMEC._FSAGradRho(self.VMEC_Data.nfp, Rho[ii], self.VMEC_Data, self.VMEC_DerivedQuant, kMinRho=self.kMinRho)
        # end for
        self.Vol = _np.cumsum(dVoldrho)
        self.gradrho = gradrho
        self.dVdrho = dVoldrho
        return self.gradrho, self.dVdrho, self.Vol


    # ========================================================================= #

    @staticmethod
    def _FSAGradRho2(M, Rho, VMEC_Data, VMEC_DerivedQuant, kMinRho=2.0e-9):
#    def _FSAGradRho2(self, Rho, kMinRho=2.0e-9):
        """
            Calculate the gradient of the normalized effective radius at a single-point in space
            ... gradrho(\Rho, \Theta, \Phi)
        """
#        M = self.VMEC_Data.nfp
        def _Gradrho2(Theta, Phi):
            RR, ZZ, dR_dRho, dR_dTheta, dZ_dRho, dZ_dTheta, dR_dZi, dZ_dZi = \
                GetLabCoordsFromVMEC(Rho, Theta, Phi, VMEC_Data, VMEC_DerivedQuant, kMinRho=kMinRho, nargout=8)

            _gr2, sqg = PYVMEC._GradRhoSquared(RR, dR_dRho, dR_dTheta, dZ_dRho, dZ_dTheta, dR_dZi, dZ_dZi)
            return _gr2*sqg

        # simpson's integration over toroidal and poloidal angle (one field period)
#        dVdrho = self._dVdrho(Rho, kMinRho=kMinRho)
        dVdrho = PYVMEC._dVdrho(M, Rho, VMEC_Data, VMEC_DerivedQuant, kMinRho=kMinRho)
        gradrho2 = M*_ut.aquad.dblquad_asr(_Gradrho2, 0.0, 2.0*_np.pi, 0.0, 2.0*_np.pi/M, tol=1e-3)
#        gradrho2 = M*_ut.aquad.dblquad_scipy(_Gradrho2, 0.0, 2.0*_np.pi, 0.0, 2.0*_np.pi/M, tol=1e-3)
#        gradrho2 = M*_ut.aquad.dblquadrature(_Gradrho2, 0.0, 2.0*_np.pi, 0.0, 2.0*_np.pi/M, tol=1e-3)

#        # due to speed reasons, let's try to mesh it and calculate
#        ph = _np.linspace(0.0, 2.0*_np.pi, num=10)/M
##        func = lambda x: _Gradrho2(x, ph)
##        gradrho2 = M*_ut.aquad.quad_asr(func, 0.0, 2.0*_np.pi, tol=1e-3)
#        th = _np.linspace(0.0, 2.0*_np.pi, num=15)
##        func = lambda x: _Gradrho2(th, x)
##        gradrho2 = M*_ut.aquad.quad_asr(func, 0.0, 2.0*_np.pi/M, tol=1e-3)
#        res = 0.0
#        dphi = ph[1]-ph[0]
#        dth = th[1]-th[0]
#        for _ph in ph:
#            for _th in th:
#                res += _Gradrho2(_th, _ph)
#            # end for
#        # end for
#        gradrho2 = res*dphi*dth
        return gradrho2/dVdrho, dVdrho

#    @staticmethod
#    def FSAGradRho2(M, Rho, VMEC_Data, VMEC_DerivedQuant, kMinRho=2.0e-9):
    def FSAGradRho2(self, Rho, **kwargs):
        """
            Calculate the flux-surface averaged gradient of the normalized
            effective radius on a single surface
            ... gradrho(\Rho)
        """
        kMinRho = kwargs.setdefault('kMinRho',None)
        if kMinRho is None and hasattr(self, 'kMinRho'):
            kMinRho = self.kMinRho
        else:
            kMinRho = 2.0e-9
        # end if
        self.kMinRho = kMinRho

        Rho = _np.atleast_1d(Rho).copy()
        gradrho2 = _np.zeros_like(Rho)
        dVoldrho = _np.zeros_like(Rho)
        nrho = len(gradrho2)
        for ii in range(nrho):
            if Rho[ii] == 0:
                Rho[ii] = 1e-5
            gradrho2[ii], dVoldrho[ii] = self._FSAGradRho2(self.VMEC_Data.nfp, Rho[ii], self.VMEC_Data, self.VMEC_DerivedQuant, kMinRho=kMinRho)
        # end for
        self.Vol = _np.cumsum(dVoldrho)
        self.gradrho2 = gradrho2
        self.dVdrho = dVoldrho
        return self.gradrho2, self.dVdrho, self.Vol
    # end def

 # end class VMEC_Data


# ========================================================================= #


if __name__=="__main__":
    import time
    start = time.time()

#    import os as _os
    rootdir = _os.path.join('d:/', 'Workshop', 'TRAVIS', 'MagnConfigs', 'W7X')
    if not _os.path.exists(rootdir):
        rootdir = _os.path.join('G:/', 'Workshop', 'TRAVIS_tree', 'MagnConfigs', 'W7X')
    # kVMECfile = _os.path.join(rootdir, 'wout_w7x.1000_1000_1000_1000_+0390_+0000.05.0144.txt')
#    kVMECfile = _os.path.join(rootdir, 'wout_w7x.1000_1000_1000_1000_+0000_+0000.01.00.txt')
    kVMECfile = _os.path.join(rootdir, 'wout_w7x.1000_1000_1000_1000_+0390_+0000.05.0144.nc')

#    import os as _os
#    kVMECfile = _os.path.join('G://', 'Workshop','TRAVIS_tree','MagnConfigs','W7X')
#    kVMECfile = _os.path.join(kVMECfile, 'wout_w7x.1000_1000_1000_1000_+0000_+0000.10.06.txt')

#    VMEC_Data = extract_data(kVMECfile, ForceReadVMEC=False, verbose=True)
#    VMEC_Data = extract_data(kVMECfile, ForceReadVMEC=True, verbose=True)

#    VMEC_DerivedQuant = spline_data(VMEC_Data, ForceReadVMEC=False, verbose=True)
#    VMEC_DerivedQuant = spline_data(VMEC_Data, ForceReadVMEC=True, verbose=True)

    # within the LCFS
    cLabPos = _np.atleast_2d([5.85, 0, 0]);
    sPos = _np.atleast_2d([_np.sqrt(0.11407577361150142), 180.0*_np.pi/180.0, 0.0])
    rho, th = Cart2Flux(cLabPos, kVMECfile, ForceReadVMEC=False, PosTol=1e-3, verbose=True)

    fiRZ = _ut.cart2pol(cLabPos[0,0],cLabPos[0,1],cLabPos[0,2])
    pLabOut = Flux2Polar(sPos, kVMECfile, verbose=True)

    # # outside the LCFS
    # cLabPos2 = _np.atleast_2d([5.5, 0, 0]);
    # rho2, th2 = Cart2Flux(cLabPos2, kVMECfile, ForceReadVMEC=False, PosTol=1e-3, verbose=True)

    vmc = PYVMEC(kVMECfile)
    # _rho, _th = vmc(cLabPos)
    # _rho2, _th2 = vmc(cLabPos2)

    VMEC_Data = vmc.VMEC_Data
    VMEC_DerivedQuant = vmc.VMEC_DerivedQuant
    VMEC_DataSource = vmc.fname

    print('Starting single-point calculation of gradrho')
    # _gradrho, _ = vmc._FSAGradRho(5, sPos[0,0], VMEC_Data, VMEC_DerivedQuant, kMinRho=2.0e-9)

    _localgradrho = vmc.LocalGradRho(sPos[0,0], sPos[0,1], sPos[0,2], VMEC_Data, VMEC_DerivedQuant, kMinRho=2.0e-9)

    theta = _np.linspace(0, 2*_np.pi, num=30, endpoint=True)
    phi = _np.linspace(0.0, 2*_np.pi/5.0, num=30, endpoint=True)
    lgradrho = _np.zeros((len(theta)*len(phi),), float)
    cylR  = _np.zeros((len(theta)*len(phi),), float)
    cylfi = _np.zeros((len(theta)*len(phi),), float)
    cylZ  = _np.zeros((len(theta)*len(phi),), float)
    kk = -1
    for ph in phi:
        for th in theta:
            kk += 1
            # pLabPos = Flux2Polar(_np.atleast_2d([sPos[0,0], th, ph]), kVMECfile, verbose=True)
            # lgradrho[kk] = vmc.LocalGradRho(sPos[0,0], th, ph, VMEC_Data, VMEC_DerivedQuant, kMinRho=2.0e-9)

            pLabPos = Flux2Polar(_np.atleast_2d([0.95, th, ph]), kVMECfile, verbose=True)
            lgradrho[kk] = vmc.LocalGradRho(0.95, th, ph, VMEC_Data, VMEC_DerivedQuant, kMinRho=2.0e-9)
            # print((kk, ph, th, lgradrho[kk]))

            cylR[kk], cylfi[kk], cylZ[kk] = tuple(pLabPos.copy())
            print((kk, ph, th, cylR[kk], cylfi[kk], cylZ[kk], lgradrho[kk]))
        # end for
    # end for

    # ======= #

    import matplotlib.pyplot as _plt
    from mpl_toolkits.mplot3d import Axes3D   # analysis:ignore

    XX, YY, ZZ, = _ut.pol2cart(cylR, cylfi, cylZ)
    _lgr = _np.copy(lgradrho)


    norm=_plt.Normalize(vmin=_lgr.min(), vmax=_lgr.max(), clip=False)
    colors = _plt.cm.viridis(norm(_lgr))

    # _plt.figure()
    # _ax = _plt.axes(projection='3d')
    # _ax.scatter(XX, YY, ZZ, c=_lgr, cmap='viridis')
    # _ax.set_xlabel(r'x [m]')
    # _ax.set_ylabel(r'y [m]')
    # _ax.set_zlabel(r'z [m]')


    XX = XX.reshape((len(theta), len(phi)))
    YY = YY.reshape((len(theta), len(phi)))
    ZZ = ZZ.reshape((len(theta), len(phi)))
    lgr = _lgr.reshape((len(theta), len(phi)))
    cm = colors.reshape((len(theta), len(phi),4))

    _plt.figure()
    _ax = _plt.axes(projection='3d')
    surf = _ax.plot_surface(XX, YY, ZZ, antialiased=True, rstride=1, cstride=1, facecolors=cm)
    # _ax.plot_surface(XX, YY, ZZ, facecolors=colors)
    # _ax.plot_trisurf(XX, YY, ZZ, facecolors=colors, edgecolor='none')
    _ax.set_xlabel(r'x [m]')
    _ax.set_ylabel(r'y [m]')
    _ax.set_zlabel(r'z [m]')
    # _ax.set_clabel(r'$\nabla\rho$')
    _ax.set_title(r'$|\nabla\rho|$ at $\rho$=0.95')

    m = _plt.cm.ScalarMappable(cmap=surf.cmap, norm=surf.norm)
    m.set_array(lgr)
    _plt.colorbar(m, shrink=0.5, aspect=5)

    # lim = (max(abs(max(_np.max(XX), _np.max(YY), _np.max(ZZ))), abs(min(_np.min(XX), _np.min(YY), _np.min(ZZ)))))
    # _ax.set_xlim(-lim, lim)
    # _ax.set_ylim(-lim, lim)
    # _ax.set_zlim(-lim, lim)

#     print('Starting single-point calculation of gradrho2')
#     _gradrho2, _dVoldrho = vmc._FSAGradRho2(5, sPos[0,0], VMEC_Data, VMEC_DerivedQuant, kMinRho=2.0e-9)

#     rho = _np.linspace(1e-3, 1.0, 21, endpoint=True)
#     print('Starting multi-point calculation of gradrho')
#     gradrho, dVoldrho, Vol = vmc.FSAGradRho(rho)
#     print('Starting multi-point calculation of gradrho2')
#     gradrho2, dVoldrho, Vol = vmc.FSAGradRho2(rho)

#     from pymconf import pyvmec as pymconf
#     pym = pymconf(kVMECfile)
#     pym.getgradrho(sPos[0,0])
#     pym.getVol(sPos[0,0])
#     print((_dVoldrho, pym.dVdrho.flatten(), 100*((_np.abs(_dVoldrho)-pym.dVdrho)/pym.dVdrho).flatten()))
#     print((_gradrho,  pym.gradrho.flatten(), 100*((_np.abs(_gradrho) -pym.gradrho)/pym.gradrho).flatten()))
#     print((_gradrho2, pym.gradrho2.flatten(), 100*((_np.abs(_gradrho2)-pym.gradrho2)/pym.gradrho2).flatten()))
#     pym.getgradrho(rho)
#     pym.getVol(rho)
#     print((Vol, pym.Vol.flatten(), 100*((_np.abs(Vol)-pym.Vol)/pym.Vol).flatten()))
#     print((dVoldrho, pym.dVdrho.flatten(), 100*((_np.abs(dVoldrho)-pym.dVdrho)/pym.dVdrho).flatten()))
#     print((gradrho,  pym.gradrho.flatten(), 100*((_np.abs(gradrho) -pym.gradrho)/pym.gradrho).flatten()))
#     print((gradrho2, pym.gradrho2.flatten(), 100*((_np.abs(gradrho2)-pym.gradrho2)/pym.gradrho2).flatten()))

# end if

# ========================================================================= #
# ========================================================================= #

