function e = getDilatation(r,t,p)
% This function returns the volume dilatation (the trace of the 
% deformation or strain tensor) given continuum mechanics theory. Here, 
% we have solved for this analytically.

    rtot = p("Rtot");
    r0 = p("r0");
    K = p("K");
    lambda = p("lambda");
    mu = p("mu");
    phi = p("phi");
    Q = p("Q");
    nmax = p("nmax");

    e = zeros(length(r), length(t));
    c = K * (lambda + 2*mu);
    zeta_cons = (Q * phi / (2 * pi * c * rtot * r0)) * (1 - r0/rtot);
    for ii = 1:nmax
        mn = ii * pi / rtot;
        zeta_n = zeta_cons * mn^2 * besseljs(0, mn * r0);
        T = exp(-(c * mn^2 * t));
        R = besseljs(0, mn * r);
        e = e + zeta_n*R*T;
    end
end