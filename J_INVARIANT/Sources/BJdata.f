       module BJdata
       use stel_kinds
       implicit none
       integer, parameter :: max_bmns = 120
       real(rprec), dimension(max_bmns) :: xm_rdcd, xn_rdcd
       integer, dimension(max_bmns) :: m_rdcd, n_rdcd
       real(rprec), dimension(max_bmns) :: blocal_rdcd
       real(rprec) theta0, phi0, iota_local, ep_mu
       real(rprec) TWOPI,PI
       real(rprec) length_factor
       integer istop, J_star_opt
       character*3 device
       end module BJdata
