function [DYDT] = getConcentration_v2(t, Y, r, v_spline, dvdr_spline, p_spline, params)
    %%
    % takes in velocity, velocity gradient and params and calculates the
    % concentration in space and time using convection diffusion equation
    % in spherical coordinates. Since this function creates its own time
    % vector adaptively, we pass the velocity and velocity gradient, which
    % were calculated in CM, as splines and evaluate the spline at the
    % determined time point here
    %%
    
    dYdt = zeros(numel(r),size(Y,2));
    dYdr = zeros(numel(r),size(Y,2));
    d2Ydr2 = zeros(numel(r),size(Y,2));
    
    %numr = params("numr");
    numr = numel(r);
    D_S = params("D_S");
    c0 = params("c0");
    Q = params("Q");
    r0 = params("r0");
    ind_r0 = params("ind_r0");
    t_inj_end = params("t_inj_end");
    %dCdt(i) = D_T.*(-2*C(i)+C(i-1))./(dx.^2) - k_b.*C(i); %BC zero conc CORRECT

    if t > t_inj_end
        for i = 1 %y direction BC zero flux (dcdx=0)
            dr = r(2)-r(1);
            dYdt(i,:) = 6*D_S.*(Y(i+1,:)-Y(i,:))./(dr.^2);
%             if (ppval(v_spline{i},t) > 0)
%                 v = 2*Q/(pi*r0^2);
%                 dvdr = 0;
%             else
%                 v = 0;
%                 dvdr = 0;
%             end
    
            % ZERO FLUX B.C.
            %dYdr(i,:) = 0; %forward difference
            %d2Ydr2(i,:) = (2*Y(i+1,:)-2*Y(i,:))./(dr.^2); % central difference
    
            %dYdr(i,:) = (Y(i+1,:)-c0)./(dr); %forward difference
            %d2Ydr2(i,:) = (Y(i+2,:)-2*(Y(i+1,:))+c0)./(dr.^2); %forward difference
            %dYdr(i,:) = (Y(i+1,:)-c0)./(2*dr); % central difference
            %d2Ydr2(i,:) = (Y(i+1,:)-c0)./(dr.^2); % central difference
            %dCdt(i) = D_T.*(-2*C(i)+C(i-1))./(dx.^2) - k_b.*C(i); %BC zero conc CORRECT
            %dYdr(i,:) = 0;
            %d2Ydr2(i,:) = 0;
            %dYdt(i,:) = ((2*D_C/r(i)).*dYdr(i,:))+(D_C*d2Ydr2(i,:));
            %dYdt(i,:) = ((2*D_S/r(i)).*dYdr(i,:))+(D_S*d2Ydr2(i,:))-((2/r(i)).*(v.*Y(i,:)))-(Y(i,:).*dvdr)-(v.*dYdr(i,:));
            %dYdt(i,:) = ((2*D_S/r(i)).*dYdr(i,:))+(D_S*d2Ydr2(i,:))-((2/r(i)).*(v.*Y(i,:)))-(Y(i,:).*dvdr)-(v.*dYdr(i,:));
            %dYdt(i,:) = 6*D_S.*(Y(i+1,:)-Y(i,:))./(dr.^2);
            %dYdt(i,:) = 0;
        end
        for i = 2:ind_r0-1
            dr_p = r(i+1)-r(i);
            dr_n = r(i)-r(i-1);
%             if (ppval(v_spline{i},t) > 0)
%                 v = Q/(4*pi*r(i)^2);
%                 dvdr = -(2*Q)/(4*pi*r(i)^3);
%             else
%                 v = 0;
%                 dvdr = 0;
%             end
            dYdr(i,:) = (Y(i+1,:)-Y(i-1,:))./(dr_p+dr_n); % central difference
            %d2Ydr2(i,:) = (Y(i+1,:)-2*Y(i,:)+Y(i-1,:))./(dr_p*dr_n); % central difference
            %d2Ydr2(i,:) = 2*(Y(i+1,:)+2*Y(i,:)+Y(i-1,:)-dYdr(i,:)*(dr_p-dr_n))./(dr_p.^2+dr_n.^2); % central difference
            d2Ydr2(i,:) = ((2*(dr_n/(dr_p+dr_n))*Y(i+1,:))-2*Y(i,:)+(2*(dr_p/(dr_p+dr_n)*Y(i-1,:))))./(dr_p*dr_n);
            dYdt(i,:) = ((2*D_S/r(i)).*dYdr(i,:))+(D_S*d2Ydr2(i,:));
        end
        for i = ind_r0
            dr_p = r(i+1)-r(i);
            dr_n = r(i)-r(i-1);
            %dr = r(i)-r(i-1);
            %v = ppval(v_spline{i},t);
            %dvdr = ppval(dvdr_spline{i},t);
            dYdr(i,:) = (Y(i+1,:)-Y(i-1,:))./(dr_p+dr_n); % central difference
            %d2Ydr2(i,:) = (Y(i+1,:)-2*Y(i,:)+Y(i-1,:))./(dr_p*dr_n); % central difference
            %d2Ydr2(i,:) = 2*(Y(i+1,:)+2*Y(i,:)+Y(i-1,:)-dYdr(i,:)*(dr_p-dr_n))./(dr_p.^2+dr_n.^2); % central difference
            d2Ydr2(i,:) = ((2*(dr_n/(dr_p+dr_n))*Y(i+1,:))-2*Y(i,:)+(2*(dr_p/(dr_p+dr_n)*Y(i-1,:))))./(dr_p*dr_n);
            dYdt(i,:) = ((2*D_S/r(i)).*dYdr(i,:))+(D_S*d2Ydr2(i,:));
        end
    end
    for i = ind_r0+1
        %dr = r(i+1)-r(i);
        dr_p = r(i+1)-r(i);
        dr_n = r(i)-r(i-1);
        v = ppval(v_spline{i},t);
        dvdr = ppval(dvdr_spline{i},t);
        if t <= t_inj_end
            dYdr(i,:) = (Y(i+1,:)-c0)./(dr_p+dr_n); % central difference
            %d2Ydr2(i,:) = 2*(Y(i+1,:)+2*Y(i,:)+c0-(dYdr(i,:)*(dr_p-dr_n)))./(dr_p.^2+dr_n.^2); % central difference
            %d2Ydr2(i,:) = (Y(i+1,:)-2*Y(i,:)+c0)./(dr_p*dr_n); % central difference0
            d2Ydr2(i,:) = ((2*(dr_n/(dr_p+dr_n))*Y(i+1,:))-2*Y(i,:)+(2*(dr_p/(dr_p+dr_n)*c0)))./(dr_p*dr_n);
        else
            dYdr(i,:) = (Y(i+1,:)-Y(i-1,:))./(dr_p+dr_n); % central difference
            %d2Ydr2(i,:) = 2*(Y(i+1,:)+2*Y(i,:)+Y(i-1,:)-(dYdr(i,:)*(dr_p-dr_n)))./(dr_p.^2+dr_n.^2); % central difference
            %d2Ydr2(i,:) = (Y(i+1,:)-2*Y(i,:)+Y(i-1,:))./(dr_p*dr_n); % central difference
            d2Ydr2(i,:) = ((2*(dr_n/(dr_p+dr_n))*Y(i+1,:))-2*Y(i,:)+(2*(dr_p/(dr_p+dr_n)*Y(i-1,:))))./(dr_p*dr_n);
        end
        dYdt(i,:) = ((2*D_S/r(i)).*dYdr(i,:))+(D_S*d2Ydr2(i,:))-((2/r(i)).*(v.*Y(i,:)))-(Y(i,:).*dvdr)-(v.*dYdr(i,:));
    end
    for i = ind_r0+2:(numr-1) %in stroma
        %dr = r(i)-r(i-1);
        dr_p = r(i+1)-r(i);
        dr_n = r(i)-r(i-1);
        %v = spline(t_cm,v_cm(i,:),t);
        %dvdr = spline(t_cm,dvdr_cm(i,:),t);
        v = ppval(v_spline{i},t);
        dvdr = ppval(dvdr_spline{i},t);

        %dYdr(i) = (Y(i)-Y(i-1))./dr; %backward difference
        dYdr(i,:) = (Y(i+1,:)-Y(i-1,:))./(dr_p+dr_n); %central difference
        %dYdr(i,:) = (Y(i+1,:)-Y(i,:))./dr; %forward difference
        %d2Ydr2(i,:) = (Y(i+1,:)-2.*Y(i,:)+Y(i-1,:))./(dr_p*dr_n);
        d2Ydr2(i,:) = ((2*(dr_n/(dr_p+dr_n))*Y(i+1,:))-2*Y(i,:)+(2*(dr_p/(dr_p+dr_n)*Y(i-1,:))))./(dr_p*dr_n);
        %d2Ydr2(i,:) = 2*(Y(i+1,:)+2*Y(i,:)+Y(i-1,:)-(dYdr(i,:)*(dr_p-dr_n)))./(dr_p.^2+dr_n.^2); % central difference
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
        dYdr(i,:) = -Y(i-1)./dr; % 0 CONC B.C. BACKWARD DIFF
        d2Ydr2(i,:) = (-2.*Y(i-1,:)+Y(i-2,:))./(dr.^2); % 0 CONC B.C. BACKWARD DIFF
        %dYdr(i,:) = -Y(i-1,:)./(2*dr); % 0 CONC B.C. CENTRAL DIFF
        %d2Ydr2(i,:) = (Y(i-1,:))./(dr.^2); % 0 CONC B.C. CENTRAL DIFF
        dYdt(i,:) = ((2*D_S/r(i)).*dYdr(i,:))+(D_S*d2Ydr2(i,:))-((2/r(i)).*(v.*Y(i,:)))-(Y(i,:).*dvdr)-(v.*dYdr(i,:));
    end
    
    DYDT = dYdt;

end