function [DYDT] = getConcentration(t, Y, r, r_cavity_spline, t_cm, r_cm, v_cm, dvdr_cm, params)

    dYdt = zeros(numel(r),size(Y,2));
    dYdr = zeros(numel(r),size(Y,2));
    d2Ydr2 = zeros(numel(r),size(Y,2));
    
    numr = params("numr");
    D_C = params("D_C");
    D_S = params("D_S");
    P = params("P");
    Q = params("Q");
    r0 = params("r0");
    %Rtot = params("Rtot");
    
%     nt = binarySearchBin(r_cavity(:,1),t);
%     if (nt == -1)
%         nt = numel(r_cavity(:,1));
%     end
%     ind_Y_r_cavity = binarySearchBin(r,r_cavity(nt,2));

    cavity_radius = ppval(r_cavity_spline,t);
            
    ind_r_cavity = binarySearchBin(r,cavity_radius);
    if (ind_r_cavity == -1)
        ind_r_cavity = numel(r);
    end
    ind_Y_r_cavity = ind_r_cavity;

%     v = interp2(t_cm,r_cm,v_cm,t*ones(1,length(t_cm)),r);
%     dvdr = interp2(r_cm,t_cm,dvdr_cm,t*ones(1,length(t_cm)),r);

%     ind_r_cavity = helper.binarySearchBin(r,cavity_radius);
%     if (ind_r_cavity == -1)
%         ind_r_cavity = numel(r);
%     end
%     ind_Y_r_cavity = ind_r_cavity;

    %Concentration at cavity/stroma interface
    %C_intf_CSa = (D_C.*Y(ind_Y_r_cavity)+(D_S.*Y(ind_Y_r_cavity+1)))./((D_S.*P)+D_C); %Right before interface
    %C_intf_CSb = (D_C.*Y(ind_Y_r_cavity)+(D_S.*Y(ind_Y_r_cavity+1)))./((D_C./P)+D_S); %Right after interface
    C_intf_CSa = (D_C.*Y(ind_Y_r_cavity)+(D_S.*Y(ind_Y_r_cavity+1)))./((D_S./P)+D_C); %Right before interface
    C_intf_CSb = (D_C.*Y(ind_Y_r_cavity)+(D_S.*Y(ind_Y_r_cavity+1)))./((D_C.*P)+D_S); %Right after interface
    %C_intf_CSa = (D_C.*Y(ind_Y_r_cavity)-(D_C.*Y(ind_Y_r_cavity-1)))./((D_S./P)-D_S); %Right before interface
    %C_intf_CSb = (D_C.*Y(ind_Y_r_cavity)-(D_S.*Y(ind_Y_r_cavity-1)))./((D_S.*P)-D_S); %Right after interface
    
    for i = 1 %y direction BC zero flux (dcdx=0)
        dr = r(2)-r(1);
        dYdt(i,:) = 6*D_C.*(Y(i+1,:)-Y(i,:))./(dr.^2);
    end
    for i = 2:ind_Y_r_cavity-1 %inside cavity
        dr = r(i+1)-r(i);
        v = Q/(4*pi*r(i)^2);
        dvdr = -Q/(2*pi*r(i)^3);
        dYdr(i,:) = (Y(i,:)-Y(i-1,:))./dr; %backward difference
        %dYdr(i,:) = (Y(i+1,:)-Y(i-1,:))./(2*dr); %central difference
        %dYdr(i,:) = (Y(i+1,:)-Y(i,:))./dr; %forward difference
        d2Ydr2(i,:) = (Y(i+1,:)-2.*Y(i,:)+Y(i-1,:))./(dr.^2);
        %dYdt(i,:) = ((2*D_C/r(i)).*dYdr(i,:))+(D_C*d2Ydr2(i,:));
        dYdt(i,:) = ((2*D_C/r(i)).*dYdr(i,:))+(D_C*d2Ydr2(i,:))-((2/r(i)).*(v.*Y(i,:)))-(Y(i,:).*dvdr)-(v.*dYdr(i,:));
    end
    for i = ind_Y_r_cavity %right before interface - cavity/stroma
        dr = r(i+1)-r(i);
        v = Q/(4*pi*r(i)^2);
        dvdr = -Q/(2*pi*r(i)^3);
        dYdr(i,:) = (Y(i,:)-Y(i-1,:))./dr; %backward difference
        %dYdr(i,:) = (C_intf_CSa-Y(i-1,:))./(2*dr); %central difference
        %dYdr(i,:) = (C_intf_CSa-Y(i,:))./dr; %forward difference
        d2Ydr2(i,:) = (Y(i-1,:)-2.*Y(i,:)+C_intf_CSa)./(dr.^2);
        %dYdt(i,:) = ((2*D_C/r(i)).*dYdr(i,:))+(D_C*d2Ydr2(i,:));
        dYdt(i,:) = ((2*D_C/r(i)).*dYdr(i,:))+(D_C*d2Ydr2(i,:))-((2/r(i)).*(v.*Y(i,:)))-(Y(i,:).*dvdr)-(v.*dYdr(i,:));
    end
    for i = (ind_Y_r_cavity+1) %right after interface - cavity/stroma
        dr = r(i+1)-r(i);
        v = spline(t_cm,v_cm(i-20,:),t);
        dvdr = spline(t_cm,dvdr_cm(i-20,:),t);
        dYdr(i,:) = (Y(i,:)-C_intf_CSb)./dr; %backward difference
        %dYdr(i,:) = (Y(i+1,:)-C_intf_CSb)./(2*dr); %central difference
        %dYdr(i,:) = (Y(i+1,:)-Y(i,:))./dr; %forward difference
        d2Ydr2(i,:) = ((C_intf_CSb)-2.*Y(i,:)+Y(i+1,:))./(dr.^2);
        %dYdt(i,:) = ((2*D_C/r(i)).*dYdr(i,:))+(D_C*d2Ydr2(i,:));
        dYdt(i,:) = ((2*D_S/r(i)).*dYdr(i,:))+(D_S*d2Ydr2(i,:))-((2/r(i)).*(v.*Y(i,:)))-(Y(i,:).*dvdr)-(v.*dYdr(i,:));
    end
    for i = ind_Y_r_cavity + 2 : (numr - 1) %in stroma
        dr = r(i+1)-r(i);
        v = spline(t_cm,v_cm(i-20,:),t);
        dvdr = spline(t_cm,dvdr_cm(i-20,:),t);
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
        dYdt(i,:) = 2*D_S*(Y(i-1,:)-Y(i,:))./(dr.^2);
        %dYdt(i,:) = 0;
    end
    
    DYDT = dYdt;
end