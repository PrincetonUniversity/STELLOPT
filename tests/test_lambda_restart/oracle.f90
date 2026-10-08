program lambda_restart_oracle
    use vparams
    use vmec_params, only: lamscale, mscale
    use vmec_main, only: phipf, sp, sm, pressure, iota, shaped
    implicit none
    real(rprec) :: rmn(5,0:0,0:2,1), zmn(5,0:0,0:2,1)
    real(rprec) :: lmn(5,0:0,0:2,1), s, theta, expected, observed
    real(rprec) :: lambda_error, field_error, geometry_error, scale, phi
    logical :: reset
    integer :: js, k
    character(len=64) :: argument, filename
    call get_command_argument(1, argument)
    read(argument,*) scale
    call get_command_argument(2, argument)
    read(argument,*) phi
    call get_command_argument(3, argument)
    shaped=trim(argument)=="shaped"
    lamscale=scale; phipf=phi; pressure=123.0_rprec; iota=.7_rprec
    sm=0; sp=0
    do js=2,5
        s=real(js-1,rprec)/4.0_rprec
        sm(js)=sqrt(s-.125_rprec)/sqrt(s)
        if (js<5) sp(js)=sqrt(s+.125_rprec)/sqrt(s)
    end do
    sp(1)=sm(2)
    filename='synthetic.nc'; reset=.true.
    call load_xc_from_wout(rmn,zmn,lmn,reset,0,2,5,filename)
    if (reset) stop 2
    lambda_error=0; field_error=0; geometry_error=0
    do js=2,5
        s=real(js-1,rprec)/4.0_rprec
        lambda_error=max(lambda_error,abs(lmn(js,0,1,1) &
            -phi*.1_rprec*sqrt(s)/(scale*mscale(1))))
        lambda_error=max(lambda_error,abs(lmn(js,0,2,1) &
            -phi*.2_rprec*s/(scale*mscale(2))))
        geometry_error=max(geometry_error,abs(rmn(js,0,0,1)-6.2_rprec), &
            abs(mscale(1)*rmn(js,0,1,1)-.62_rprec*sqrt(s)), &
            abs(mscale(1)*zmn(js,0,1,1) &
                -merge(1.7_rprec,1.0_rprec,shaped)*.62_rprec*sqrt(s)), &
            abs(mscale(2)*rmn(js,0,2,1)-merge(.08_rprec*s,0.0_rprec,shaped)))
        do k=0,16
            theta=2.0_rprec*acos(-1.0_rprec)*real(k,rprec)/17.0_rprec
            expected=phi*(1+.1_rprec*sqrt(s)*cos(theta) &
                +.4_rprec*s*cos(2*theta))
            observed=phipf(js)+scale*(mscale(1)*lmn(js,0,1,1)*cos(theta) &
                +2*mscale(2)*lmn(js,0,2,1)*cos(2*theta))
            field_error=max(field_error,abs(observed-expected))
        end do
    end do
    print '(a,es24.16)', 'lambda_error=',lambda_error
    print '(a,es24.16)', 'field_density_error=',field_error
    print '(a,es24.16)', 'geometry_error=',geometry_error
    if (any(pressure/=123.0_rprec)) stop 3
    if (any(iota/=.7_rprec)) stop 4
    if (geometry_error>1.e-13_rprec) stop 5
    if (lambda_error>1.e-13_rprec) stop 6
    if (field_error>1.e-13_rprec) stop 7
end program
