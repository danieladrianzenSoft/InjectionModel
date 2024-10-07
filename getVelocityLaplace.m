function v = getVelocityLaplace(s,r,e,p)
    % this is the result of the fluid velocity in laplace space from Netti et
    % al., 2013. We later get the solution in space and time using an
    % inverse laplace function solver.
    Rtot = p("Rtot");
    r0 = p("r0");
    K = p("K");
    lambda = p("lambda");
    mu = p("mu");
    %phi_0 = p("phi");
    %phi = (e + phi_0) / (1 + e);
    phi = p("phi");
    Q = p("Q");

    c = K*(lambda+2*mu);
    al = (s/c).^(1/2);

    D = (Q*phi)./(4*pi*r0^2*c*s);
    %D = (Q*phi)./(pi*r0^2*c*s);
    C = -(phi*r0*s)/(4*mu*K);
    F = (D*r0^2)./(((1-(r0*C)).*(sinh(al*(Rtot-r0))))...
         +al*r0.*cosh(al*(Rtot-r0)));
    A = ((r0./al).*cosh(al*(r0-Rtot))-((1./(al.^2))+((2*mu+lambda)/(4*mu))*r0^2).*sinh(al*(r0-Rtot)));

    v = F.*(c/phi-s./(al.^2)).*(sinh(al*(Rtot-r))/(r^2))+F.*((c/phi)*al-(s./al)).*(cosh(al*(Rtot-r))/(r))+(F.*s.*A)/r^2;
end