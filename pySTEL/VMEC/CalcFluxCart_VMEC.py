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
import scipy.interpolate as _dsi
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

# ========================================================================== #
# ========================================================================== #

VMEC_Data = None
VMEC_DerivedQuant = None
VMEC_DataSource = None

def __extract_data(VMEC_FilePath, ForceReadVMEC=False, verbose=True):
    global VMEC_Data
    if VMEC_Data is None or ForceReadVMEC:
        VMEC_Data = read_vmec(VMEC_FilePath)

        VMEC_Data = _ut.Struct() # Instantiate an empty class of type structure

        # If this is a netCDF VMEC file, read the output ourseleves, otherwise,
        # let Sam Lazerson's code handle it
        #  nFileName = len( VMEC_FilePath )
        # if VMEC_FilePath[nFileName-2:nFileName].lower()=='.nc':
        if VMEC_FilePath.lower().find('.nc')>-1:
            if verbose:   print('Checking file')  # end if
            ncID = netcdf.open(VMEC_FilePath, 'nc_nowrite')

            if verbose:   print('Loading data from VMEC wout file')  # end if
            VarID = netcdf.inqVarID(ncID, 'ns')
            VMEC_Data.ns = double(netcdf.getVar(ncID, VarID))

            VarID = netcdf.inqVarID(ncID, 'mnmax')
            VMEC_Data.mnmax = double(netcdf.getVar(ncID, VarID))

            # (mn_mode) toridal mode numbers
            VarID = netcdf.inqVarID(ncID, 'xn')
            VMEC_Data.xn = double(netcdf.getVar(ncID, VarID))

            # (mn_mode) polidal mode numbers
            VarID = netcdf.inqVarID(ncID, 'xm')
            VMEC_Data.xm = double(netcdf.getVar(ncID, VarID))

            # (radius mn_mode ) cosmn component of cylindrical R, full mesh in [m]
            VarID = netcdf.inqVarID(ncID, 'rmnc')
            VMEC_Data.rmnc = (netcdf.getVar(ncID, VarID)).T

            # (radius mn_mode ) sinmn component of cylindrical Z, full mesh in [m]
            VarID = netcdf.inqVarID(ncID, 'zmns')
            VMEC_Data.zmns = (netcdf.getVar(ncID, VarID)).T
        else:
            # Use Sam Lazerson's code to read the VMEC file, but note that he
            # defines xn with the opposite sign as is usual in VMEC
            VMEC_Data = read_vmec(VMEC_FilePath)

#            VMEC_Data.xm = VMEC_Data.xm.T
#            VMEC_Data.xn = -VMEC_Data.xn.T
#            VMEC_Data.ns = VMEC_Data.ns.T
#            VMEC_Data.zmns = VMEC_Data.zmns.T
#            VMEC_Data.rmnc = VMEC_Data.rmnc.T
            VMEC_Data.xn = -VMEC_Data.xn
        #end
        VMEC_Data.loaded = True
        if verbose:
            print('Done Reading VMEC File')
        # end if
    # end if
    return VMEC_Data

def __spline_data(VMEC_Data, ForceReadVMEC=False, verbose=True):
    global VMEC_DerivedQuant
    if VMEC_DerivedQuant is None or ForceReadVMEC:
        if verbose:
            print('Normalizing.')
        # Set up radial grids.  VMEC produces data that is equally spaced in 'S'
        # (normalized total flux) space.  We often want it in normalized minor
        # radius (rho), though.
        VMEC_DerivedQuant = _ut.Struct()
#        VMEC_DerivedQuant.GridS = ((1:VMEC_Data.ns:1)-1)/(VMEC_Data.ns-1)
        VMEC_DerivedQuant.GridS = (_np.asarray(range(1, VMEC_Data.ns))-1.0)/(VMEC_Data.ns-1.0)
        VMEC_DerivedQuant.GridRho = _np.sqrt(VMEC_DerivedQuant.GridS)

        # Normalizing the spectral coefficients to rho^m helps the splining get
        # the radial dependence correct near the axis, as described in the paper.
        # See equations 7a and 7b and the paragraphs directly above those
        # equations for more details. Notationally, rmnc_norm corresponds with
        # \tilde{R}_mn, for example. The values right at the axis are linearly
        # extrapolated: \tidle{R_(mn)}(\rho = 0) = \tilde{R_(mn)}(\rho = \Delta \rho) -
        # \Delta \rho \frac{\tilde{R_(mn)(2\Delta \rho) - \tilde{R_(mn)(\Delta \rho}{\Delta \rho}}
        TempExp = _np.dot( _np.abs(VMEC_Data.xm), _np.ones((1,len(VMEC_DerivedQuant.GridRho)), float) )
        TempBase = _np.dot( _np.ones((len(VMEC_Data.xm), 1), float), _np.atleast_2d(VMEC_DerivedQuant.GridRho))
        NormalizationFactor = TempBase**TempExp

        rmnc_norm = _np.zeros_like(VMEC_Data.rmnc)
        rmnc_norm[:,1:] = VMEC_Data.rmnc[1:, :].T/NormalizationFactor[:, 1:]
        rmnc_norm[:,0] = 2.0*rmnc_norm[:,1] - rmnc_norm[:,2]

        zmns_norm = _np.zeros_like(VMEC_Data.zmnc)
        zmns_norm[:,1:] = VMEC_Data.zmns[1:, :].T/NormalizationFactor[:, 1:]
        zmns_norm[:,0] = 2.0*zmns_norm[:,1] - zmns_norm[:,2]

        if verbose:
            print('Done normalizing.')
            print('Generating splines for spectral cooefficients.')
        # end if

        # We will capture the radial dependence of the spectral normalized
        # components of r and z using splines fit to 'double-sided' profiles.
        # This approximates the radial partial derivatives of the normalized
        # spectral components, too (but not directly the radial derivatives of r,z).
        Spline_Rhos = _ut.cylsym_odd(VMEC_DerivedQuant.GridRho)
        rmnc_vals_to_spline = _ut.cylsym_even(rmnc_norm, axis=1)
        zmns_vals_to_spline = _ut.cylsym_even(zmns_norm, axis=1)
        # Spline_Rhos = np.concatenate((_np.flipud( -VMEC_DerivedQuant.GridRho[1:]), VMEC_DerivedQuant.GridRho), axis=1)
        # rmnc_vals_to_spline = np.concatenate( ( np.fliplr( rmnc_norm[ :, 1:-1 ] ) , rmnc_norm ), axis=1)
        # zmns_vals_to_spline = np.concatenate( ( np.fliplr( zmns_norm[ :, 1:-1 ] ) , zmns_norm ), axis=1)

        VMEC_DerivedQuant.rmnc_norm_spline = _dsi.spline(Spline_Rhos, rmnc_vals_to_spline)
        VMEC_DerivedQuant.rmnc_norm_spline_derivs = _dsi.spline_deriv(VMEC_DerivedQuant.rmnc_norm_spline, order=1) #First derivative

        VMEC_DerivedQuant.zmns_norm_spline = _dsi.spline(Spline_Rhos,zmns_vals_to_spline)
        VMEC_DerivedQuant.zmns_norm_spline_derivs = _dsi.spline_deriv(VMEC_DerivedQuant.zmns_norm_spline, order=1)

        if verbose:
            print('Done splining the spectral coefficients')
        # end if
    # end if
    return VMEC_DerivedQuant

def CalcFluxCart_VMEC(cLabPos, VMEC_FilePath, ForceReadVMEC=False, PosTol=1e-4, verbose=True):
    """
     Description:
       This function finds the coordinates of a point in VMEC flux coordinates that
       corresponds with a point given in cartesian lab coordinates.  The
       algorithm for this routine is described in detail in the paper:
       Attenberger, Houlberg, and Hirshman, J Comp Phys 72 (1987) 435.

       This version is based on one created by John Schmitt, but is optimized
       for Matlab (vectorized) and somewhat simplified.  It works out to be
       quite fast: 100ms for the first iteration and 10ms for subsequent
       iterations on my 4 year old Core Duo desktop.

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
            VMEC_Data =  __extract_data(VMEC_FilePath, ForceReadVMEC=ForceReadVMEC, verbose=verbose)
#        if VMEC_DerivedQuant is None:
            VMEC_DerivedQuant = __spline_data(VMEC_Data, ForceReadVMEC=ForceReadVMEC, verbose=verbose)
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


def FindVMEC_Coords(pLabPos, VMEC_Data, VMEC_DerivedQuant, PosTol=1e-4, RhoGuess=0.5, ThetaGuess=1.0, verbose=True,
                    kMaxIters=150, kMinJacDet=1.0e-16, kMaxItersPerThetaStep=30, kMaxInitialThetaStep=0.5):
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
    MaxThetaStep = kMaxInitialThetaStep

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
        if (IterNum>1) and (Dist_tp>Dist_BestMatch):
          StepReduction = 1.0

          # Keep halving until you beat your last guest
          while (Dist_tp>Dist_BestMatch) and (IterNum<kMaxIters):
            StepReduction = 2.0*StepReduction
            ApproxJacDet = 0.75*JacDet_BestMatch + 0.25*JacDet_tp

            ThetaStep = ((dZ_dRho_BestMatch*DistR_BestMatch-dR_dRho_BestMatch*DistZ_BestMatch)
                            /StepReduction/ApproxJacDet)
            if (_np.abs(ThetaStep)>MaxThetaStep):
                ThetaStep = MaxThetaStep*_ut.sign(ThetaStep)
            # end if
            Theta_tp = Theta_BestMatch + ThetaStep

            Rho_tp = (Rho_BestMatch +
                        (dR_dTheta_BestMatch*DistZ_BestMatch-dZ_dTheta_BestMatch*DistR_BestMatch)
                          /StepReduction/ApproxJacDet)

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
        # end if distance grreater than best match and iteration number greater than 1

        # If necessary, pick the next test point using the typical Newton's
        # method. Otherwise, return the approximate location in VMEC coordinates.
        if (Dist_tp<PosTol):
            Converged = True
            Rho = Rho_tp
            Theta = Theta_tp%(2.0*_np.pi)
            if (Theta<0):
                Theta = 2.0*_np.pi + Theta
            #end
        else:
            # This is the typical case: We have just identified an improved
            # estimate of Rho and Theta and about to take another normal step.
            # We only typically need to know the coordinates of the best match
            # so far, but the half-step case needs the Jacobian elements and the
            # distance that the test point was from its target.
            Rho_BestMatch = Rho_tp
            Theta_BestMatch = Theta_tp
            Dist_BestMatch = Dist_tp
            DistR_BestMatch = DistR_tp
            DistZ_BestMatch = DistZ_tp
            dR_dRho_BestMatch = dR_dRho_tp
            dR_dTheta_BestMatch = dR_dTheta_tp
            dZ_dRho_BestMatch = dZ_dRho_tp
            dZ_dTheta_BestMatch = dZ_dTheta_tp
            JacDet_BestMatch = JacDet_tp

            # The normal step, as taken from eq. 3,4 in the paper.  I'm not sure
            # why exactly they formulated the problem this way, but it works out
            # to be identical to what you would expect from Newton's method.
            # (Note that you have to use the fact that \Rho and \Theta are
            # orthagonal to get the algebra to work out).  Without a maximum
            # step size in the theta direction, convergence is very slow or
            # nonexistant.
            Rho_tp += (dR_dTheta_tp*DistZ_tp-dZ_dTheta_tp*DistR_tp)/JacDet_tp
            ThetaStep = (dZ_dRho_tp*DistR_tp-dR_dRho_tp*DistZ_tp)/JacDet_tp
            if (_np.abs(ThetaStep)>MaxThetaStep):
                ThetaStep = MaxThetaStep*_ut.sign(ThetaStep)
            #end
            Theta_tp += ThetaStep

            # If the solver tries to look at the "negative rho", we must correct
            # it and the theta step that goes along with it.
            if (Rho_tp < 0):
                Rho_tp = -Rho_tp
                Theta_tp = Theta_BestMatch + _np.pi - ThetaStep
            #end
        #end
    #end

    if (IterNum == kMaxIters):
        print('Attempt to find flux coordinate failed to converge:\n Phi=%6.4e\n R=%6.4e\n Z=%6.4e'
               %(pLabPos[1],pLabPos[0],pLabPos[2]))
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
    rmnc_norm = _dsi.ppval(VMEC_DerivedQuant.rmnc_norm_spline, RhoSpline)
    zmns_norm = _dsi.ppval(VMEC_DerivedQuant.zmns_norm_spline, RhoSpline)

    RhoDeriv_rmnc_norm = _dsi.ppval(VMEC_DerivedQuant.rmnc_norm_spline_derivs, RhoSpline)
    RhoDeriv_zmns_norm = _dsi.ppval(VMEC_DerivedQuant.zmns_norm_spline_derivs, RhoSpline)

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
        SpectCoeffNorm = _np.ones((len(VMEC_Data.xm),1), _np.float64)
        dSpectCoeffNorm_dRho = _np.zeros((len(VMEC_Data.xm),1), _np.float64)

        tst = (VMEC_Data.xm==1) and (VMEC_Data.xn==0)
        if (tst).any():
            SpectCoeffNorm[tst] = _np.copy(Rho)
            dSpectCoeffNorm_dRho[tst] = 1
        # end if tst
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


def _RootG(RR, dR_dRho, dR_dTheta, dZ_dRho, dZ_dTheta):
    """
        Return the jacobian for flux surface average integration (flux surface area element)
        ... this is the sqrt(g)
    """
    return RR*(dR_dTheta*dZ_dRho - dR_dRho*dZ_dTheta)


def _dVdrho(M, Rho, VMEC_Data, VMEC_DerivedQuant, kMinRho=2.0e-9):
    """
    Return the incremental volume between \Rho and \Rho + d\Rho
    """
    def _sqrtG(Theta, Phi):
        RR, _, dR_dRho, dR_dTheta, dZ_dRho, dZ_dTheta = \
            GetLabCoordsFromVMEC(Rho, Theta, Phi, VMEC_Data, VMEC_DerivedQuant,
                                 kMinRho=kMinRho, nargout=6)
        return _RootG(RR, dR_dRho, dR_dTheta, dZ_dRho, dZ_dTheta)

    def zifunc(Phi):
        def thfunc(x, Phi):
            return _sqrtG(x, Phi)
        # simpson's integration over poloidal angle
        # ... note that this excludes the point theta=0 AND theta=2pi!!!
        return _ut.openpoints(thfunc, 0.0, 2.0*_np.pi, TOL=1e-6, verbose=True)

    # simpson's integration over toroidal angle (one field period)
    # ... note that this excludes the point Phi=0 AND Phi=2pi!!!
    # Returns the incremental volume between \Rho and \Rho + d\Rho
    dVdrho = M*_ut.openpoints(zifunc, 0.0, 2.0*_np.pi/M, TOL=1e-6, verbose=True)
    return dVdrho


def GetdVdrho(M, Rho, VMEC_Data, VMEC_DerivedQuant, kMinRho=2.0e-9):
    """
    Return the incremental volume between \Rho and \Rho + d\Rho for a vector of
    input flux-surfaces
    """
    dVdrho = _np.zeros_like(Rho)
    nrho = len(dVdrho)
    for ii in range(nrho):
        dVdrho[ii] = _dVdrho(M, Rho[ii], VMEC_Data, VMEC_DerivedQuant, kMinRho=kMinRho)
    # end for
    Vol = _np.cumsum(dVdrho)
    return dVdrho, Vol


# ========================================================================= #


def _GradRhoSquared(RR, dR_dRho, dR_dTheta, dZ_dRho, dZ_dTheta, dR_dZi, dZ_dZi):
    """
        Calculate the gradient of the normalized effective radius at a single-point in space
        ... gradrho(\Rho, \Theta, \Phi)
    """
    numerator = (dR_dTheta*dR_dTheta+dZ_dTheta*dZ_dTheta)*(dR_dZi*dR_dZi+RR*RR+dZ_dZi*dZ_dZi)
    numerator -= (dR_dTheta*dR_dZi+dZ_dTheta*dZ_dZi)**2.0

    sqrtG = _RootG(RR, dR_dRho, dR_dTheta, dZ_dRho, dZ_dTheta)
    return numerator/(sqrtG*sqrtG)

def _FSAGradRho(M, Rho, VMEC_Data, VMEC_DerivedQuant, kMinRho=2.0e-9):
    """
        Calculate the gradient of the normalized effective radius at a single-point in space
        ... gradrho(\Rho, \Theta, \Phi)
    """
    def _Gradrho(Theta, Phi):
        RR, ZZ, dR_dRho, dR_dTheta, dZ_dRho, dZ_dTheta, dR_dZi, dZ_dZi = \
            GetLabCoordsFromVMEC(Rho, Theta, Phi, VMEC_Data, VMEC_DerivedQuant,
                                 kMinRho=kMinRho, nargout=8)
        return (_RootG(RR, dR_dRho, dR_dTheta, dZ_dRho, dZ_dTheta),
                _np.sqrt(_GradRhoSquared(RR, dR_dRho, dR_dTheta, dZ_dRho, dZ_dTheta, dR_dZi, dZ_dZi)))

    def zifunc(Phi):
        def thfunc(x, Phi):
            _sqrtG, _gr = _Gradrho(x, Phi)
            return _sqrtG*_gr
        # simpson's integration over poloidal angle
        # ... note that this excludes the point theta=0 AND theta=2pi!!!
        return _ut.openpoints(thfunc, 0.0, 2.0*_np.pi, TOL=1e-6, verbose=True)

    # simpson's integration over toroidal angle (one field period)
    # ... note that this excludes the point Phi=0 AND Phi=2pi!!!
    gradrho = M*_ut.openpoints(zifunc, 0.0, 2.0*_np.pi/M, TOL=1e-6, verbose=True)
    gradrho /= _dVdrho(M, Rho, VMEC_Data, VMEC_DerivedQuant, kMinRho=kMinRho)
    return gradrho

def FSAGradRho(M, Rho, VMEC_Data, VMEC_DerivedQuant, kMinRho=2.0e-9):
    """
        Calculate the flux-surface averaged gradient of the normalized
        effective radius on a single surface
        ... gradrho(\Rho)
    """
    gradrho = _np.zeros_like(Rho)
    dVoldrho = _np.zeros_like(Rho)
    nrho = len(gradrho)
    for ii in range(nrho):
        gradrho[ii], dVoldrho[ii] = _FSAGradRho(M, Rho[ii], VMEC_Data, VMEC_DerivedQuant, kMinRho=kMinRho)
    # end for
    Vol = _np.cumsum(dVoldrho)
    return gradrho, dVoldrho, Vol


# ========================================================================= #


def _FSAGradRho2(M, Rho, VMEC_Data, VMEC_DerivedQuant, kMinRho=2.0e-9):
    """
        Calculate the gradient of the normalized effective radius at a single-point in space
        ... gradrho(\Rho, \Theta, \Phi)
    """
    def _Gradrho2(Theta, Phi):
        RR, ZZ, dR_dRho, dR_dTheta, dZ_dRho, dZ_dTheta, dR_dZi, dZ_dZi = \
            GetLabCoordsFromVMEC(Rho, Theta, Phi, VMEC_Data, VMEC_DerivedQuant,
                                 kMinRho=kMinRho, nargout=8)
        return (_RootG(RR, dR_dRho, dR_dTheta, dZ_dRho, dZ_dTheta),
                _GradRhoSquared(RR, dR_dRho, dR_dTheta, dZ_dRho, dZ_dTheta, dR_dZi, dZ_dZi))

    def zifunc(Phi):
        def thfunc(x, Phi):
            _sqrtG, _gr2 = _Gradrho2(x, Phi)
            return _sqrtG*_gr2
        # simpson's integration over poloidal angle
        # ... note that this excludes the point theta=0 AND theta=2pi!!!
        return _ut.openpoints(thfunc, 0.0, 2.0*_np.pi, TOL=1e-6, verbose=True)

    # simpson's integration over toroidal angle (one field period)
    # ... note that this excludes the point Phi=0 AND Phi=2pi!!!
    gradrho2 = M*_ut.openpoints(zifunc, 0.0, 2.0*_np.pi/M, TOL=1e-6, verbose=True)
    dVoldrho = _dVdrho(M, Rho, VMEC_Data, VMEC_DerivedQuant, kMinRho=kMinRho)
    gradrho2 /= dVoldrho
    return gradrho2, dVoldrho


def FSAGradRho2(M, Rho, VMEC_Data, VMEC_DerivedQuant, kMinRho=2.0e-9):
    """
        Calculate the flux-surface averaged gradient of the normalized
        effective radius on a single surface
        ... gradrho(\Rho)
    """
    # Theta, dTheta = _np.linspace(0.0, 2.0*_np.pi, num=50, endpoint=True, retstep=True)

    gradrho2 = _np.zeros_like(Rho)
    dVoldrho = _np.zeros_like(Rho)
    nrho = len(gradrho2)
    for ii in range(nrho):
        gradrho2[ii], dVoldrho[ii] = _FSAGradRho2(M, Rho[ii], VMEC_Data, VMEC_DerivedQuant, kMinRho=kMinRho)
    # end for
    Vol = _np.cumsum(dVoldrho)
    return gradrho2, dVoldrho, Vol


# ========================================================================= #
# ========================================================================= #

#class VMEC_Data(Struct):
#    def __init__(self):
#        loaded = False
#        self.ns = float()
#        self.mnmax = float()
#        self.xn = np.array([], float)
#        self.xm = np.array([], float)
#        self.rmnc = np.array([], float)
#        self.zmnc = np.array([], float)
#    # end def __init__
# # end class VMEC_Data


# ========================================================================= #


if __name__=="__main__":
    import os as _os
    kVMECfile = _os.path.join('G://', 'Workshop','TRAVIS_tree','MagnConfigs','W7X')
    kVMECfile = _os.path.join(kVMECfile, 'wout_w7x.1000_1000_1000_1000_+0000_+0000.10.06.txt')

    VMEC_Data = __extract_data(kVMECfile, ForceReadVMEC=False, verbose=True)
    VMEC_Data = __extract_data(kVMECfile, ForceReadVMEC=True, verbose=True)

    VMEC_DerivedQuant = __spline_data(VMEC_Data, ForceReadVMEC=False, verbose=True)
    VMEC_DerivedQuant = __spline_data(VMEC_Data, ForceReadVMEC=True, verbose=True)

    tstpos = _np.atleast_2d([5.5, 0, 0])
    rho, th = CalcFluxCart_VMEC(cLabPos, VMEC_FilePath, ForceReadVMEC=False, PosTol=1e-4, verbose=True)
# end if

# ========================================================================= #
# ========================================================================= #

