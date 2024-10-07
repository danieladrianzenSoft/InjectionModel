function e = getDilatationLaplace(s,r,p)
    % this is the result of the dilatation in laplace space from Netti et
    % al., 2013. We later get the solution in space and time using an
    % inverse laplace function solver.
    Rtot = p("Rtot");
    r0 = p("r0");
    K = p("K");
    lambda = p("lambda");
    mu = p("mu");
    phi = p("phi");
    Q = p("Q");

    c = K * (lambda + 2 * mu);
    D = (Q * phi) ./ (4 * pi * r0^2 * c * s);
    %D = (Q * phi) ./ (pi * r0^2 * c * s);
    C = -(phi * r0 * s) / (4 * mu * K);

    al = (s/c).^(1/2);

    e = ((D * r0^2) ./ r)...
         .* (sinh(al*(Rtot - r))...
         ./ (((1-(r0*C)) .* (sinh(al*(Rtot-r0))))...
         + al.*r0.*cosh(al.*(Rtot-r0))));
end

%(((((0.1/10^3)/60*0.2)/(4*3.1416*0.035^2*2.0003e-05*s))*0.035^2)/0.035)*(sinh((s/2.0003e-05)^(1/2)*(1-0.035))/((1-0.035*(-(0.2*0.035*s)/(4*6.58*1e4*7.6e-11)))*(sinh((s/2.0003e-05)^(1/2)*(1-0.035))+(s/2.0003e-05)^(1/2)*0.035*cosh((s/2.0003e-05)^(1/2)*(1-0.035)))))