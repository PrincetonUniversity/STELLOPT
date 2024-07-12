!EZcdf, Easy Interface to netCDF Routine Calls
!=============================================
 
!The module is available through the NTCC Webpage,
!w3.pppl.gov/NTTC  under "Modules Library", as compressed
!tarfile, ezcdf.tar.gz, and as zip archive, ezcdf.zip.
!Alternatively it can be obtained from ftp.pppl.gov
!in pub/NTCC/.
 
 
! AUTHORS
 
!Conceived by 7/98 by Sunitinder Sekhon
!Modified by J. Menard 12/98 to run on Cray C90
!Completely re-written by C.Ludescher 2/99
!Added complex support (64 and 128 bit) by A. Pletzer 5/01
 
 
! CONTACT
 
! C. Ludescher cludescher@pppl.gov
! A. Pletzer   pletzer@pppl.gov
 
 
! REVISION HISTORY
 
!      date         Description
 
! February  1999  -- Created
! April     2000  -- A.Pletzer: Added R4
! May 01,   2000  -- C. Ludescher: Simplified by adding module ezcdf
! May 17,   2001  -- A. Pletzer: Added C8 and C16
! Interface for cdfopn to handle optional argument ier
! 04/28/00 C.Ludescher
! + ezcdf_close for symmetry (ap)
! September 2002  -- S. Hirshman, added aliases to mimic F90 I/O routines
!                    added cdf_inquire with OPTIONAL xtype argument
!                    added 'LOG' data type (nf_int) to handle logicals
!                    (user may use 'INT' interchangeably)
!
! April 2003      -- removed ONLY and some rename clauses from USE stmts
!                    (D. McCune).  These caused problems, because an
!                    application with "use ezcdf" would indirectly see
!                    "use ezcdf_inqvar" with various ONLY clauses; it
!                    was unclear which one would rule, but, at load time
!                    it was clear that the ONLY clauses were making some
!                    legitimate ezcdf entries invisible, at least when
!                    built using Lahey-Fujitsu fortran.  Also, the
!                    indirect presence of multiple "USE ezcdf_opencls"
!                    with rename lists caused problems which were resolved
!                    by removing the rename lists ans using INTERFACE
!                    statements inside ezcdf_opncls instead.
!
 
MODULE ezcdf

  ! No aliases. this caused the Intel compiler to fail, so I had to duplicate
  ! 2 interfaces: cdfPutVar <=> cdf_write and cdfGetVar <=> cdf_read (pletzer)

  USE ezcdf_GenPut
  USE ezcdf_GenGet
  USE ezcdf_attrib
  USE ezcdf_opncls

END MODULE ezcdf
