!-----------------------------------------------------------------------
!     Subroutine:    stellopt_write_header
!     Authors:       S. Lazerson (samuel.lazerson@gauss-fusion.com)
!     Date:          05/31/2024
!     Description:   This subroutine writest the STELLOPT header
!                    information.
!-----------------------------------------------------------------------
      SUBROUTINE stellopt_write_header
      USE stellopt_runtime, ONLY: bigno
      USE stellopt_vars, ONLY: equil_type
      USE stellopt_targets, ONLY: txport_proxy, sigma_orbit, &
         sigma_bootstrap, sigma_balloon, sigma_kink, sigma_ece, &
         sigma_coil_bnorm, sigma_regcoil_chi2_b, sigma_dkes, &
         sigma_dkes_Erdiff, sigma_dkes_alpha, sigma_fluxloop, &
         sigma_bprobe, sigma_segrog, sigma_neo, sigma_txport, &
         sigma_regcoil_current_density
      USE mpi_params
      USE diagno_runtime, ONLY: DIAGNO_VERSION
      USE beams3d_runtime, ONLY: BEAMS3D_VERSION
      USE vmec_params, ONLY: version_vmec=> version_
      USE vmec0, ONLY: version_bootsj => version_
!DEC$ IF DEFINED (GENE)
      USE par_other, ONLY: svn_gene => git_master, release_gene => git_branch
!DEC$ ENDIF
      IMPLICIT NONE
      ! Print code messages
      CALL tolower(equil_type)
      IF ((myid == master) .and. (TRIM(equil_type(1:4)) == 'vmec') ) THEN
         WRITE(6,*)        " Equilibrium calculation provided by: "
         WRITE(6,"(2X,A)") "================================================================================="
         WRITE(6,"(2X,A)") "=========   Parallel Variational Moments Equilibrium Code (v "//TRIM(version_vmec)//")       ========="
         WRITE(6,"(2X,A)") "=========                (S. Hirshman, J. Whitson)                      ========="
         WRITE(6,"(2X,A)") "=========         http://vmecwiki.pppl.wikispaces.net/VMEC              ========="
         WRITE(6,"(2X,A)") "================================================================================="
         WRITE(6,*)        "    "
      END IF
!DEC$ IF DEFINED (BEAMS3D_OPT)
      IF (myid == master .and. ANY(sigma_orbit < bigno) ) THEN
         WRITE(6,*)               " Energetic Particle calculation provided by: "
         WRITE(6,"(2X,A)")        "================================================================================="
         WRITE(6,"(2X,A,F5.2,A)") "=========                      BEAMS3D (v",BEAMS3D_VERSION,")                         ========="
         WRITE(6,"(2X,A)")        "=========                  (M. McMillan, S. Lazerson)                   ========="
         WRITE(6,"(2X,A)")        "=========                       lazerson@pppl.gov                       ========="
         WRITE(6,"(2X,A)")        "=========          http://vmecwiki.pppl.wikispaces.net/BEAMS3D          ========="
         WRITE(6,"(2X,A)")        "================================================================================="
         WRITE(6,*)        "    "
      END IF
!DEC$ ELSE
      IF (ANY(sigma_orbit < bigno)) THEN
         sigma_orbit = bigno
         IF (myid == master) THEN
            WRITE(6,*) '!!!!!!!!!!!!!!!!!!!! WARNING !!!!!!!!!!!!!!!!!!!!!!!!!'
            WRITE(6,*) '  STELLOPT has not been linked to the BEAMS3D code.   '
            WRITE(6,*) '  Optimization of particle orbits not possible.'
            WRITE(6,*) '  Disabling energetic particle targets.'
         END IF
      END IF
!DEC$ ENDIF
      IF (myid == master .and. ANY(sigma_bootstrap < bigno)) THEN
         WRITE(6,*)        " Bootstrap calculation provided by: "
         WRITE(6,"(2X,A)") "================================================================================="
         WRITE(6,"(2X,A,A,A)") "=========                    BOOTSJ (v",version_bootsj,")                             ========="
         WRITE(6,"(2X,A)") "=========            (J. Tolliver, K. Shaing, P. Moroz)                 ========="
         WRITE(6,"(2X,A)") "=========        http://vmecwiki.pppl.wikispaces.net/BOOTSJ             ========="
         WRITE(6,"(2X,A)") "================================================================================="
         WRITE(6,*)        "    "
      END IF
      IF (myid == master .and. ANY(sigma_balloon < bigno)) THEN
         WRITE(6,*)        " Ballooning stability calculation provided by: "
         WRITE(6,"(2X,A)") "================================================================================="
         WRITE(6,"(2X,A,F5.2,A)") "=========                    COBRAVMEC (v",4.10,")                         ========="
         WRITE(6,"(2X,A)") "=========                 (R. Sanchez, S. Hirshman)                     ========="
         WRITE(6,"(2X,A)") "=========                   raul.sanchez@uc3m.es                        ========="
         WRITE(6,"(2X,A)") "================================================================================="
         WRITE(6,*)        "    "
      END IF
!DEC$ IF DEFINED (TERPSICHORE)
      IF (myid == master .and. ANY(sigma_kink < bigno)) THEN
         WRITE(6,*)        " Kink stability calculation provided by: "
         WRITE(6,"(2X,A)") "================================================================================="
         WRITE(6,"(2X,A,F5.2,A)") "=========                    TERPSICHORE (v2016)                        ========="
         WRITE(6,"(2X,A)") "=========    (D. V. ANDERSON, W. A. COOPER, R. GRUBER AND U. SCHWENN)   ========="
         WRITE(6,"(2X,A)") "=========                   wilfred.cooper@epfl.ch                      ========="
         WRITE(6,"(2X,A)") "================================================================================="
         WRITE(6,*)        "    "         
      END IF
!DEC$ ELSE
      IF (ANY(sigma_kink < bigno)) THEN
         sigma_kink(:) = bigno
         IF (myid == master) THEN
            WRITE(6,*) '!!!!!!!!!!!!!!!!!!!! WARNING !!!!!!!!!!!!!!!!!!!!!!!!!'
            WRITE(6,*) '  STELLOPT has not been linked to the TERPSICHORE code.   '
            WRITE(6,*) '  Optimization of kink stability not possible.'
            WRITE(6,*) '  Disabling kink stability targets.'
         END IF
      END IF
!DEC$ ENDIF
!DEC$ IF DEFINED (TRAVIS)
      IF (myid == master .and. ANY(sigma_ece < bigno)) THEN
         WRITE(6,*)        " ECE Radiation calculation provided by: "
         WRITE(6,"(2X,A)")        "================================================================================="
         CALL printversion_sopt_f77
         !WRITE(6,"(2X,A,F5.2,A)") "=========                            TRAVIS                             ========="
         WRITE(6,"(2X,A)")        "=========                    (N. Marushchenko)                          ========="
         WRITE(6,"(2X,A)")        "=========              nikolai.marushchenko@ipp.mpg.de                  ========="
         WRITE(6,"(2X,A)")        "================================================================================="
         WRITE(6,*)               "    "
      END IF
!DEC$ ELSE
      IF (ANY(sigma_ece < bigno)) THEN
         sigma_ece = bigno
         IF (myid == master) THEN
            WRITE(6,*) '!!!!!!!!!!!!!!!!!!!! WARNING !!!!!!!!!!!!!!!!!!!!!!!!!'
            WRITE(6,*) '  STELLOPT has not been linked to the TRAVIS code.   '
            WRITE(6,*) '  Optimization of ECE Radiation not possible.'
            WRITE(6,*) '  Disabling ECE Radiation targets.'
         END IF
      END IF
!DEC$ ENDIF
!DEC$ IF DEFINED (COILOPTPP)
      IF (myid == master .and. (sigma_coil_bnorm < bigno)) THEN
         WRITE(6,*)        " Stellarator Coil Optimization provided by: "
         WRITE(6,"(2X,A)") "================================================================================="
         WRITE(6,"(2X,A)") "=========                            COILOPT++                          ========="
         WRITE(6,"(2X,A)") "=========                    (J. Breslau, S. Lazerson)                  ========="
         WRITE(6,"(2X,A)") "=========                        jbreslau@pppl.gov                      ========="
         WRITE(6,"(2X,A)") "================================================================================="
         WRITE(6,*)        "    "
      END IF
!DEC$ ELSE
      IF (sigma_coil_bnorm < bigno) THEN
         sigma_coil_bnorm = bigno
         IF (myid == master) THEN
            WRITE(6,*) '!!!!!!!!!!!!!!!!!!!! WARNING !!!!!!!!!!!!!!!!!!!!!!!!!'
            WRITE(6,*) '  Coil optimization with the COILOPT++'
            WRITE(6,*) '  code has been disabled.  Coil optimziation'
            WRITE(6,*) '  has been turned off.  Contact your vendor for'
            WRITE(6,*) '  further information.'
         END IF
      END IF
!DEC$ ENDIF
!DEC$ IF DEFINED (REGCOIL)
      IF (myid == master .and. (ANY(sigma_regcoil_chi2_b < bigno) .or. &
                                (sigma_regcoil_current_density < bigno) )) THEN
         WRITE(6,*)        " Stellarator REGCOIL Optimization provided by: "
         WRITE(6,"(2X,A)") "================================================================================="
         WRITE(6,"(2X,A)") "=========                            REGCOIL                            ========="
         WRITE(6,"(2X,A)") "=========                        (M. Landreman)                         ========="
         WRITE(6,"(2X,A)") "=========               Matt dot Landreman at gmail dot com             ========="
         WRITE(6,"(2X,A)") "================================================================================="
         WRITE(6,*)        "    "
      END IF
!DEC$ ELSE
      IF (myid == master .and. (ANY(sigma_regcoil_chi2_b < bigno) .or. &
                                (sigma_regcoil_current_density < bigno) ) ) THEN
         sigma_regcoil_chi2_b = bigno
         sigma_regcoil_current_density = bigno
         WRITE(6,*) '!!!!!!!!!!!!!!!!!!!! WARNING !!!!!!!!!!!!!!!!!!!!!!!!!'
         WRITE(6,*) '  Coil optimization with the REGCOIL'
         WRITE(6,*) '  code has been disabled.  Coil optimziation'
         WRITE(6,*) '  has been turned off.  Contact your vendor for'
         WRITE(6,*) '  further information.'
      END IF
!DEC$ ENDIF
!DEC$ IF DEFINED (DKES_OPT)
      IF (myid == master .and. ( ANY(sigma_dkes < bigno) .or. &
                                 ANY(sigma_dkes_Erdiff < bigno) .or. &
                                 ANY(sigma_dkes_alpha < bigno) ) ) THEN
         WRITE(6,*)        " Drift-Kinetic Equation Solver (DKES) provided by: "
         WRITE(6,"(2X,A)") "================================================================================="
         WRITE(6,"(2X,A)") "=========           Drift Kinetic Equation Solver, Variational          ========="
         WRITE(6,"(2X,A)") "=========                    (S.Hirshman, W. van Rij)                   ========="
         WRITE(6,"(2X,A)") "=========                      hirshmansp@ornl.gov                      ========="
         WRITE(6,"(2X,A)") "================================================================================="
         WRITE(6,*)        "    "
      END IF
!DEC$ ELSE
      IF (ANY(sigma_dkes < bigno) .or. ANY(sigma_dkes_Erdiff < bigno) .or. ANY(sigma_dkes_alpha < bigno)) THEN
         sigma_dkes(:) = bigno
         sigma_dkes_Erdiff(:) = bigno
         sigma_dkes_alpha(:) = bigno
         IF (myid == master) THEN
            WRITE(6,*) '!!!!!!!!!!!!!!!!!!!! WARNING !!!!!!!!!!!!!!!!!!!!!!!!!'
            WRITE(6,*) '  Drift-kinetic optimization with the DKES'
            WRITE(6,*) '  code has been disabled.  Drift-kinetic optimziation'
            WRITE(6,*) '  has been turned off.  Contact your vendor for'
            WRITE(6,*) '  further information.'
         END IF
      END IF
!DEC$ ENDIF
      IF (myid == master .and. (ANY(sigma_fluxloop < bigno) .or. ANY(sigma_bprobe < bigno) .or. ANY(sigma_segrog < bigno) )) THEN
         WRITE(6,*)        " Magnetic Diagnostic calculation provided by: "
         WRITE(6,"(2X,A)") "================================================================================="
         WRITE(6,"(2X,A,F5.2,A)") "=========                    DIAGNO (v",DIAGNO_VERSION,")                             ========="
         WRITE(6,"(2X,A)") "=========            (S.Lazerson, H Gardner, J. Geiger)                 ========="
         WRITE(6,"(2X,A)") "=========                   lazerson@pppl.gov                           ========="
         WRITE(6,"(2X,A)") "=========       http://vmecwiki.pppl.wikispaces.net/DIAGNO              ========="
         WRITE(6,"(2X,A)") "================================================================================="
         WRITE(6,*)        "    "
      END IF
!DEC$ IF DEFINED (NEO_OPT)
      IF (myid == master .and. ANY(sigma_neo < bigno)) THEN
         WRITE(6,*)        " Neoclassical Transport calculation provided by: "
         WRITE(6,"(2X,A)") "================================================================================="
         WRITE(6,"(2X,A)") "=========                      NEO (v3.02)                              ========="
         WRITE(6,"(2X,A)") "=========               (W.Kernbichler and S.Kasilov)                   ========="
         WRITE(6,"(2X,A)") "=========               kernbichler@itp.tu-graz.ac.at                   ========="
         WRITE(6,"(2X,A)") "================================================================================="
         WRITE(6,*)        "    "
      END IF
!DEC$ ELSE
      IF (ANY(sigma_neo < bigno)) THEN
         sigma_neo(:) = bigno
         IF (myid == master) THEN
            WRITE(6,*) '!!!!!!!!!!!!!!!!!!!! WARNING !!!!!!!!!!!!!!!!!!!!!!!!!'
            WRITE(6,*) '  Neoclassical transport optimization with the NEO'
            WRITE(6,*) '  code has been disabled.  Neoclassical optimziation'
            WRITE(6,*) '  has been turned off.  Contact your vendor for'
            WRITE(6,*) '  further information.'
         END IF
      END IF
!DEC$ ENDIF
!DEC$ IF DEFINED (TXPORT_OPT)
      IF (myid == master .and. ANY(sigma_txport < bigno)) THEN
         WRITE(6,*)        " Geometry Interface to Turbulent Transport provided by: "
         WRITE(6,"(2X,A)") "================================================================================="
         WRITE(6,"(2X,A)")       "=========        Geometry Interface for Stellarators and Tokamaks       ========="
         WRITE(6,"(2X,A)")       "=========          (P.Xanthopoulos, W.A.Cooper, and Yu.Turkin)          ========="
         WRITE(6,"(2X,A)")       "=========          pax@ipp.mpg.de  http://www.ipp.mpg.de/~pax/          ========="
         WRITE(6,"(2X,A)") "================================================================================="
         WRITE(6,*)        "    "
         WRITE(6,"(2X,A)") "     NOTICE: New TXPORT variables now used to control execution COORDINATES,"
         WRITE(6,"(2X,A)") "             IN_OUT, and SETUP namelists are ignored.  LGLOBAL_TXPORT,"
         WRITE(6,"(2X,A)") "             NZ_TXPORT, NALPHA_TXPORT, ALPHA0_TXPORT have been added to the"
         WRITE(6,"(2X,A)") "             OPTIMUM namelist."
         WRITE(6,*)        "    "
      END IF
!DEC$ ELSE
      IF (ANY(sigma_txport < bigno)) THEN
         sigma_txport(:) = bigno
         IF (myid == master) THEN
            WRITE(6,*) '!!!!!!!!!!!!!!!!!!!! WARNING !!!!!!!!!!!!!!!!!!!!!!!!!'
            WRITE(6,*) '  Turbulent transport optimization with the GIST/TXPORT'
            WRITE(6,*) '  code has been disabled.  Turbulent optimziation'
            WRITE(6,*) '  has been turned off.  Contact your vendor for'
            WRITE(6,*) '  further information.'
         END IF
      END IF
!DEC$ ENDIF
!DEC$ IF DEFINED (AEOPT)
      CALL tolower(txport_proxy)
      IF (myid == master .and. ANY(sigma_txport < bigno) .and. (TRIM(txport_proxy(1:11)) == 'availenergy') ) THEN
         WRITE(6,*)        " Turbulent Transport calculation provided by: "
         WRITE(6,"(2X,A)") "================================================================================="
         WRITE(6,"(2X,A)") "=========            Trapped Particle Available Energy Code             ========="
         WRITE(6,"(2X,A)") "=========            (R. Mackenbach, S. Lazerson, J. Proll)             ========="
         WRITE(6,"(2X,A)") "=========         https://github.com/RalfMackenbach/STELLOPT_AE         ========="
         WRITE(6,"(2X,A)") "================================================================================="
         WRITE(6,*)        "    "
      END IF
!DEC$ ELSE
      CALL tolower(txport_proxy)
      IF (ANY(sigma_txport < bigno) .and. (TRIM(txport_proxy(1:11)) == 'availenergy')) THEN
         txport_proxy = 'prox1d'
         IF (myid == master) THEN
            WRITE(6,*) '!!!!!!!!!!!!!!!!!!!! WARNING !!!!!!!!!!!!!!!!!!!!!!!!!'
            WRITE(6,*) '  STELLOPT has not been linked to the Available Energy code.      '
            WRITE(6,*) '  Optimization with Available Energy for turblent'
            WRITE(6,*) '  transport not possible.  Defaulting to proxy function'
            WRITE(6,*) '        txport_proxy = prox1d'
            WRITE(6,*) '  Code Available at:'
            WRITE(6,*) '        https://github.com/RalfMackenbach/STELLOPT_AE'
         END IF
      END IF
!DEC$ ENDIF
!DEC$ IF DEFINED (GENE)
      CALL tolower(txport_proxy)
      IF (myid == master .and. ANY(sigma_txport < bigno) .and. (TRIM(txport_proxy(1:4)) == 'gene') ) THEN
         WRITE(6,*)        " Turbulent Transport calculation provided by: "
         WRITE(6,"(2X,A)") "================================================================================="
         WRITE(6,"(2X,A)") "=========       Gyrokinetic Electromagnetic Numerical Experiment        ========="
         WRITE(6,"(2X,A)") "=========                  (F. Jenko, P. Xanthopoulos)                  ========="
         WRITE(6,"(2X,A)") "=========              http://www2.ipp.mpg.de/~fsj/gene/                ========="
         WRITE(6,"(2X,A)") "=========              GENE11  (git-branch "//release_gene//")                  ========="
         WRITE(6,"(2X,A)") "=========                      (git-master: "//svn_gene//")     ========="
         WRITE(6,"(2X,A)") "================================================================================="
         WRITE(6,*)        "    "
      END IF
!DEC$ ELSE
      CALL tolower(txport_proxy)
      IF (ANY(sigma_txport < bigno) .and. (TRIM(txport_proxy(1:4)) == 'gene')) THEN
         txport_proxy = 'prox1d'
         IF (myid == master) THEN
            WRITE(6,*) '!!!!!!!!!!!!!!!!!!!! WARNING !!!!!!!!!!!!!!!!!!!!!!!!!'
            WRITE(6,*) '  STELLOPT has not been linked to the GENE code.      '
            WRITE(6,*) '  Optimization with linear GENE for turblent'
            WRITE(6,*) '  transport not possible.  Defaulting to proxy function'
            WRITE(6,*) '        txport_proxy = prox1d'
         END IF
      END IF
!DEC$ ENDIF
      END SUBROUTINE stellopt_write_header