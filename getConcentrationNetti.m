function [DYDT] = getConcentrationNetti(t, Y, r, v_spline, dvdr_spline, params)

    dYdt = zeros(numel(r),size(Y,2));
    dYdr = zeros(numel(r),size(Y,2));
    d2Ydr2 = zeros(numel(r),size(Y,2));
    
    numr = params("numr");
    D_S = params("D_S");
    c0 = params("c0");
    Q = params("Q");
    r0 = params("r0");
        %dCdt(i) = D_T.*(-2*C(i)+C(i-1))./(dx.^2) - k_b.*C(i); %BC zero conc CORRECT

    for i = 1 %y direction BC zero flux (dcdx=0)
        dr = r(2)-r(1);
        %dYdt(i,:) = 6*D_C.*(Y(i+1,:)-Y(i,:))./(dr.^2);
        v = Q/(4*pi*r0^2);
        dvdr = -2*Q/(4*pi*r0^3);
        %dYdr(i,:) = (Y(i+1,:)-c0)./(dr); %forward difference
        %d2Ydr2(i,:) = (Y(i+2,:)-2*(Y(i+1,:))+c0)./(dr.^2); %forward difference
        dYdr(i,:) = (Y(i+1,:)-c0)./(2*dr); % central difference
        d2Ydr2(i,:) = (Y(i+1,:)-c0)./(dr.^2); % central difference
        %dCdt(i) = D_T.*(-2*C(i)+C(i-1))./(dx.^2) - k_b.*C(i); %BC zero conc CORRECT

        %dYdt(i,:) = ((2*D_C/r(i)).*dYdr(i,:))+(D_C*d2Ydr2(i,:));
        dYdt(i,:) = ((2*D_S/r(i)).*dYdr(i,:))+(D_S*d2Ydr2(i,:))-((2/r(i)).*(v.*Y(i,:)))-(Y(i,:).*dvdr)-(v.*dYdr(i,:));
    end
    for i = 2:(numr-1) %in stroma
        dr = r(i+1)-r(i);

        %v = spline(t_cm,v_cm(i,:),t);
        %dvdr = spline(t_cm,dvdr_cm(i,:),t);
        v = ppval(v_spline{i},t);
        dvdr = ppval(dvdr_spline{i},t);

        %dYdr(i) = (Y(i)-Y(i-1))./dr; %backward difference
        dYdr(i,:) = (Y(i+1,:)-Y(i-1,:))./(2*dr); %central difference
        %dYdr(i,:) = (Y(i+1,:)-Y(i,:))./dr; %forward difference
        d2Ydr2(i,:) = (Y(i+1,:)-2.*Y(i,:)+Y(i-1,:))./(dr.^2);
        %dYdt(i,:) = ((2*D_C/r(i)).*dYdr(i,:))+(D_C*d2Ydr2(i,:));
        dYdt(i,:) = ((2*D_S/r(i)).*dYdr(i,:))+(D_S*d2Ydr2(i,:))-((2/r(i)).*(v.*Y(i,:)))-(Y(i,:).*dvdr)-(v.*dYdr(i,:));
    end
    for i = numr %end of stroma, B.C.
        dr = r(i)-r(i-1);
        %dYdt(i,:) = 6*D_S*(Y(i-1,:)-Y(i,:))./(dr.^2);
        %dYdt(i,:) = D_S*(Y(i,:)-2.*Y(i-1,:)+Y(i-2,:))./(dr.^2)-(2/a)*v(i)*Y(i);
        %dYdt(i,:) = 2*D_S*(Y(i-1,:)-Y(i,:))./(dr.^2);
        %dYdt(i,:) = 2*D_S*(Y(i-1,:)-Y(i,:))./(dr.^2);
        %dYdt(i,:) = 0;
        
        v = ppval(v_spline{i},t);
        dvdr = ppval(dvdr_spline{i},t);
        %dYdr(i,:) = -Y(i-1)./dr; % 0 CONC B.C. BACKWARD DIFF
        %d2Ydr2(i,:) = (-2.*Y(i-1,:)+Y(i-2,:))./(dr.^2); % 0 CONC B.C. BACKWARD DIFF
        dYdr(i,:) = -Y(i-1)./(2*dr); % 0 CONC B.C. CENTRAL DIFF
        d2Ydr2(i,:) = (Y(i-1,:))./(dr.^2); % 0 CONC B.C. CENTRAL DIFF
        dYdt(i,:) = ((2*D_S/r(i)).*dYdr(i,:))+(D_S*d2Ydr2(i,:))-((2/r(i)).*(v.*Y(i,:)))-(Y(i,:).*dvdr)-(v.*dYdr(i,:));
    end
    
    DYDT = dYdt;

end