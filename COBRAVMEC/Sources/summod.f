       MODULE summod
       USE stel_kinds
       IMPLICIT NONE

       REAL(rprec), DIMENSION(:), ALLOCATABLE :: bfields, bfieldze,
     1  bfieldth, rboo, rs, rze, rth, zs, zze, zth, ccosi, ssine, lam1,
     2  zetang, thetang, arg, lambdaze, lambdath, bsubze, bsubs, bsubth,
     3  lambdas, rboo2, jacob2, rjac2i, bfield, bsupth, bsupze, aux,
     4  gtssub, gstsub, gzssub, gszsub, gtzsub, gztsub, gttsub, gzzsub,
     5  gsssub, gttsup, gzzsup, gtssup, gstsup, gzssup, gszsup, gtzsup,
     6  gztsup, gsssup, lam2, cks_vmec, ckth_vmec, bfieldi, bfield2i,
     7  lam3

       CONTAINS

       SUBROUTINE alloc_summod(npt)
       INTEGER :: npt, istat

       ALLOCATE(bfields(npt), bfieldze(npt), bfieldth(npt), rboo(npt),
     1  rs(npt), rze(npt), rth(npt), zs(npt), zze(npt), zth(npt),
     2  ccosi(npt), ssine(npt), lam1(npt), zetang(npt), thetang(npt),
     3  arg(npt), lambdaze(npt), lambdath(npt), bsubze(npt), bsubs(npt),
     4  bsubth(npt), lambdas(npt), rboo2(npt),jacob2(npt),bfield2i(npt),
     5  rjac2i(npt), bfield(npt), bsupth(npt), bsupze(npt), gtssub(npt),
     6  gstsub(npt), gzssub(npt), gszsub(npt), gtzsub(npt), gztsub(npt),
     7  gttsub(npt), gzzsub(npt), gsssub(npt), gttsup(npt), gzzsup(npt),
     8  gtssup(npt), gstsup(npt), gzssup(npt), gszsup(npt), gtzsup(npt),
     9  gztsup(npt), gsssup(npt), lam2(npt),  bfieldi(npt), aux(npt),
     A  lam3(npt), cks_vmec(npt), ckth_vmec(npt), stat= istat)

       IF (istat .ne. 0) STOP 'Allocation error in COBRA summod...'

       END SUBROUTINE alloc_summod

       SUBROUTINE free_summod
       INTEGER :: istat

       DEALLOCATE( bfields, bfieldze, bfieldi, bfield2i, aux, lam3,
     1  bfieldth, rboo, rs, rze, rth, zs, zze, zth, ccosi, ssine, lam1,
     2  zetang, thetang, arg, lambdaze, lambdath, bsubze, bsubs, bsubth,
     3  lambdas, rboo2, jacob2, rjac2i, bfield, bsupth, bsupze, lam2,
     4  gtssub, gstsub, gzssub, gszsub, gtzsub, gztsub, gttsub, gzzsub,
     5  gsssub, gttsup, gzzsup, gtssup, gstsup, gzssup, gszsup, gtzsup,
     6  gztsup, gsssup, cks_vmec, ckth_vmec, stat=istat)

       END SUBROUTINE free_summod

       END MODULE summod
