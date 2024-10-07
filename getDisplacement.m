function u = getDisplacement(r,t,p)
% This function returns the solid displacement u of tissue
% given continuum mechanics theory. The equation has been
% solved analytically and is given in the form of an infinite series
% here solved with nMax iterations.

    r0 = p("r0");
    rtot = p("Rtot");
    K = p("K");
    lambda = p("lambda");
    mu = p("mu");
    %c = p("c");
    phi = p("phi");
    Q = p("Q");
    nmax = p("nmax");

    u = zeros(length(r), length(t));
    c = K * (lambda + 2*mu);
    %zeta = Q * phi / (4 * pi * rtot * c);
    zeta_cons = (Q * phi / (2 * pi * c * rtot * r0)) * (1 - r0/rtot);

    for ii = 1:nmax
%         mn = ii * pi / rtot;
%         T = repmat(exp(-(c * mn^2 * t)),length(r),1);
%         R_vec = besseljs(1, mn * r)/(mn)...
%             - (besseljs(1, mn * r0)/(mn)) * ((r0^2) ./ (r.^2))...
%             + ((2*mu + lambda) / (4*mu)) * ((r0^3) ./ (r.^2)) * besseljs(0,mn*r0);
%         R = repmat(R_vec,1,length(t));
%         u = u + zeta*T.*R;
        mn = ii * pi / rtot;
        zeta_n = zeta_cons * mn^2 * besseljs(0, mn * r0);

        T = exp(-(c * mn^2 * t));
        R = besseljs(1, mn * r)/(mn)...
            + ((1./r.^2) * (2 * mu + lambda)/(4 * mu))...
            * ((r0^3) * besseljs(0,mn*r0) - (r0^2/mn) * besseljs(1,mn*r0));
        u = u + zeta_n*R*T;
    end
end