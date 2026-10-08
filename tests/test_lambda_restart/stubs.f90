module vparams
    implicit none
    integer, parameter :: rprec = kind(1.0d0)
    real(rprec), parameter :: one = 1.0_rprec, zero = 0.0_rprec
end module

module vmec_params
    use vparams
    implicit none
    integer, parameter :: ntmax = 1
    integer, parameter :: rcc=1, rss=1, rsc=1, rcs=1
    integer, parameter :: zsc=1, zcs=1, zcc=1, zss=1
    real(rprec), parameter :: mscale(0:2)=[1.0_rprec, &
        sqrt(2.0_rprec),sqrt(2.0_rprec)], nscale(0:0)=1.0_rprec
    real(rprec) :: lamscale
end module

module vmec_dim
    implicit none
    integer, parameter :: mpol1=2
end module

module vmec_input
    implicit none
    logical, parameter :: lasym=.false.
end module

module vmec_main
    use vparams
    implicit none
    logical, parameter :: lthreed=.false.
    logical :: shaped=.false.
    real(rprec), parameter :: cp5=.5_rprec
    real(rprec) :: phipf(5), sp(5), sm(5)
    real(rprec) :: pressure(5), iota(5)
end module

module parallel_include_module
    implicit none
    integer, parameter :: rank=0, t1lglob=1, t1rglob=5
    logical, parameter :: parvmec=.false.
end module

module read_wout_mod
    use vparams
    use vmec_params, only: lamscale
    use vmec_main, only: phipf, sp, sm, p5 => cp5, shaped
    implicit none
    integer, parameter :: ns=5, ntor=0, nfp=1, mnmax=3
    real(rprec), allocatable :: rmnc(:,:), zmns(:,:), lmns(:,:)
    real(rprec), allocatable :: rmns(:,:), zmnc(:,:), lmnc(:,:), xm(:), xn(:)
contains
    subroutine read_wout_file(filename, ierr)
        character(len=*), intent(in) :: filename
        integer, intent(out) :: ierr
        real(rprec) :: s, lmns1(3)
        integer :: js
        allocate(rmnc(3,5), zmns(3,5), lmns(3,5), rmns(3,5), zmnc(3,5), &
                 lmnc(3,5), xm(3), xn(3))
        xm=[0.0_rprec,1.0_rprec,2.0_rprec]; xn=0
        rmnc=0; zmns=0; lmns=0; rmns=0; zmnc=0; lmnc=0
        do js=1,ns
            s=real(js-1,rprec)/real(ns-1,rprec)
            rmnc(1,js)=6.2_rprec
            rmnc(2,js)=.62_rprec*sqrt(s)
            zmns(2,js)=.62_rprec*sqrt(s)
            if (shaped) then
                rmnc(3,js)=.08_rprec*s
                zmns(2,js)=1.7_rprec*zmns(2,js)
            end if
            ! Known internal full coefficients, independent physical lambda:
            ! lambda=.1*sqrt(s)*sin(theta)+.2*s*sin(2*theta).
            lmns1=[0.0_rprec,.1_rprec*sqrt(s),.2_rprec*s]*phipf(js)/lamscale
            include 'native_export.inc'
        end do
        include 'native_halfmesh_export.inc'
        ierr=0
    end subroutine

    subroutine read_wout_deallocate
        deallocate(rmnc, zmns, lmns, rmns, zmnc, lmnc, xm, xn)
    end subroutine
end module
