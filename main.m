clear
clc

%% simulation settings

runCavitySimulations = 1;
runConcentrationSimulations = 1;
writeResults = 1;
e_filename = 'dilatation.mat';
u_filename = 'displacement.mat';
rc_filename = 'cavity_radius.mat';
v_filename = 'velocity.mat';
dvdr_filename = 'velocity_gradient.mat';
p_filename = 'pressure.mat';
c_filename = 'concentration.mat';

plotPressureLine = 1;
plotPressure3D = 1;
plotDilatationLine = 1;
plotDisplacementLine = 1;
plotDisplacement3D = 1;
plotVelocityLine = 1;
plotVelocity3D = 1;
plotCavityRadius = 1;
plotConcentration = 1;

%% initialize simulation params

% variableParam controls the variable that will change in your parametric analysis. you can
% make it a vector or a single value if you want to test multile values
variableParam = "K";

%K = [3.8e-11, 7.6e-11, 3.8e-10, 7.6e-10]; % m^2 / kPa = cm^2/barye
%K = 1.5e-11;
K = 3e-10; % hydraulic conductivity in CGS (cm^2/barye)
r0 = 0.03; % radius of the needle in cm, Netti 2003,
%r0 = round(0.311/(2*10),3); %nominal inner diameter of 24 gauge needle is 0.2032mm. Here we change to cm. 
Rtot = 2; % total size of simulation domain in cm, Netti 2003
%Q = (0.1/1e3)/60; % 0.1uL/min = 0.1mm^3/min Netti 2003
Q = 1/3600; % flow rate in cm3/s
Vol = 0.50; %total injection volume in cm3 -> CHANGE TO 0.5ML - 2.5ML for humans
lambda = 13.16*1e4; % lame first parameter in Barye (CGS). Netti 2003, 13.16kPa = 13.16*1e4 Barye
%mu = [6.58*1e4,1*1e5]; % Netti 2003, 6.58kPa = 6.58*1e4 Barye
mu = 6.58*1e4; % lame second parameter in Barye (CGS). 
phi = 0.2; % Porosity, dimensionless. Netti 2003
%numr = 400;
t_end = 12*60*60; %last time point for simulation in seconds.
%r = linspace(r0,Rtot,numr)';
%r = r0,
%ind_r0 = 21;

relevant_r_bound = 0.2; %spatial point for higher resolution
ind_r0 = 31;
dr_r0 = r0/ind_r0;
numr = 300; % number spatial points
r_mass = [linspace(0,r0,ind_r0), linspace(r0+dr_r0,relevant_r_bound,numr/3), ...
    linspace(relevant_r_bound+0.005,1,2*numr/3-ind_r0)]; % spatial vector for mass transport
r = r_mass(ind_r0:end); % spatial vector for everything else, we only simulate starting at radius of needle r0
t_injection_end = Vol/Q; %time point at which injection ends
t_injection = linspace(0,t_injection_end,600); % time vector for duration of injection, used for CM, MT sets its own vector 
t_relaxation = [t_injection_end+1,t_end]; %span for relaxation vector

c0 = 0.5; % intial condition 
D_S = 5*10^(-6); %(cm^2/s) % diffusion coefficient in solution
D_C = D_S/2; % diffusion coefficient in cavity
P = 1; % partition coefficient at solution/cavity boundary

c = K*(lambda+2*mu); % this just simplifies calculations, function of lame parameters.

%% make param dictionary to be used in all functions

paramNames = ["K","r0","Rtot","Q","lambda","mu","phi","numr","D_C","D_S","P","c0","ind_r0","t_inj_end"];
paramIndices = 1:length(paramNames);
paramValues = {K, r0, Rtot, Q, lambda, mu, phi, numr, D_C, D_S, P, c0, ind_r0, t_injection_end};
params = dictionary(paramNames, paramValues);
paramLabels = dictionary(paramIndices,paramNames);

%% CONTINUUM MECHANICS
% calculate dilatation, solid displacement, cavity radius and pressure
% using laplace transforms

variableParamArray = params(variableParam);
lengthVariableParam = length(variableParamArray{1});

e = cell(lengthVariableParam,1); % dilatation
u = cell(lengthVariableParam,1); % solid displacement
rc = cell(lengthVariableParam,1); % cavity radius
v = cell(lengthVariableParam,1); % fluid velocity
dvdr = cell(lengthVariableParam,1); % fluid velocity gradient
p = cell(lengthVariableParam,1); % pressure
t_m = cell(lengthVariableParam,1); % mass transport time vector
c_m = cell(lengthVariableParam,1); % mass transport ethanol concentration

if runCavitySimulations == 1
    for iter = 1:lengthVariableParam
        ep = zeros(length(r),length(t_injection));
        up = zeros(length(r),length(t_injection));
        vp = zeros(length(r),length(t_injection));
        dvdrp = zeros(length(r),length(t_injection));
        pp = zeros(length(r),length(t_injection));
        paramsIterValues = getParamValues(params,paramLabels,variableParam,iter);
        paramsIter = dictionary(paramNames,paramsIterValues);

        fprintf('Continuum Mechanics: %d\n',iter)
        tic    
        for ri = 1:length(r)
            % get solution by taking inverse laplace
            ep(ri,2:end) = inverseLaplace(@(s) getDilatationLaplace(s,r(ri),paramsIter), t_injection(2:end));
            up(ri,2:end) = inverseLaplace(@(s) getDisplacementLaplace(s,r(ri),paramsIter), t_injection(2:end));
            vp(ri,2:end) = inverseLaplace(@(s) getVelocityLaplace(s,r(ri),ep(ri,:),paramsIter), t_injection(2:end));
            dvdrp(ri,2:end) = inverseLaplace(@(s) getVelocityGradLaplace(s,r(ri),ep(ri,:),paramsIter), t_injection(2:end));
    
            %% interpolating NaN values:
            % to avoid errors
            ep(ri,1) = 0;
            ep_nan = isnan(ep(ri,:));
            er_interp = interp1(t_injection(~ep_nan),ep(ri,~ep_nan),t_injection(ep_nan),'spline','extrap');
            ep(ri,ep_nan) = er_interp;
            % u
            up(ri,1) = 0;
            up_nan = isnan(up(ri,:));
            ur_interp = interp1(t_injection(~up_nan),up(ri,~up_nan),t_injection(up_nan),'spline','extrap');
            up(ri,up_nan) = ur_interp;
            % v
            vp(ri,1) = 0;
            vp_nan = isnan(vp(ri,:));
            vr_interp = interp1(t_injection(~vp_nan),vp(ri,~vp_nan),t_injection(vp_nan),'spline','extrap');
            vp(ri,vp_nan) = vr_interp;
            % dvdr
            dvdrp(ri,1) = 0;
            dvdrp_nan = isnan(dvdrp(ri,:));
            dvdrp_interp = interp1(t_injection(~dvdrp_nan),dvdrp(ri,~dvdrp_nan),t_injection(dvdrp_nan),'spline','extrap');
            dvdrp(ri,dvdrp_nan) = dvdrp_interp;
    
            if variableParam == "mu"
                pp(ri,:) = ep(ri,:)*(lambda+2*mu(iter));
            else
                pp(ri,:) = ep(ri,:)*(lambda+2*mu);
            end
    
        end
        toc
    
        rc{iter} = r0 + up(1,:); % cavity radius is first spatial point of solid displacement field
        e{iter} = ep;
        u{iter} = up;
        v{iter} = vp;
        p{iter} = pp;
        dvdr{iter} = dvdrp;
    end

end

if writeResults == 1 && runCavitySimulations == 1
    tic
    fprintf('Writing Results CM: \n')
    save(e_filename,'e','t_injection')
    save(u_filename,'u','t_injection')
    save(rc_filename,'rc','t_injection')
    save(v_filename,'v','t_injection')
    save(p_filename,'p','t_injection')
    save(dvdr_filename,'dvdr','t_injection')
    toc
end

%% CALCULATE CONCENTRATION (MASS TRANSPORT)
% using convection-diffusion equation and finite element

% S = getSparsity(numr);
S = getSparsity(length(r_mass));
S = sparse(S);
opts1 = odeset('Vectorized','on','JPattern',S,'RelTol',1e-3,'AbsTol',1e-4,'NonNegative',true);
%opts1 = odeset('Vectorized','on','JPattern',S,'NonNegative',true);

v_spline_inj = cell(length(r_mass),lengthVariableParam);
dvdr_spline_inj = cell(length(r_mass),lengthVariableParam);
p_spline_inj = cell(length(r_mass),lengthVariableParam);
rc_spline_inj = cell(cell(1,lengthVariableParam));
v_spline_rel = cell(length(r_mass),lengthVariableParam);
dvdr_spline_rel = cell(length(r_mass),lengthVariableParam);
p_spline_rel = cell(length(r_mass),lengthVariableParam);
rc_spline_rel = cell(1,lengthVariableParam); 

if runConcentrationSimulations == 1
    for iter = 1:lengthVariableParam
    
        paramsIterValues = getParamValues(params,paramLabels,variableParam,iter);
        paramsIter = dictionary(paramNames,paramsIterValues);
    
        if runCavitySimulations == 0
            load(e_filename);
            load(u_filename);
            load(rc_filename);
            load(v_filename);
            load(p_filename);
            load(dvdr_filename);
            ep = e{iter};
            up = u{iter};
            rcp = rc{iter};
            vp = v{iter};
            pp = p{iter};
            dvdrp = dvdr{iter};
            %v_spline = spline(t,vp);
            %dvdr_spline = spline(t,dvdrp);
            %rc_spline = spline(t,rcp);
        else
            vp = v{iter};
            dvdrp = dvdr{iter};
            pp = p{iter};
            rcp = rc{iter};
        end

        rc_rel = rcp(end)*ones(length(t_relaxation),1);

        fprintf('Spline creation: %d\n',iter)
        tic

        rc_spline_inj{iter} = spline(t_injection,rcp);
        rc_spline_rel{iter} = spline(t_relaxation,rc_rel);

        for i = 1:length(r_mass)
            % if within ind_r0, use qvg velocity in needle as v,
            % corresponding dvdr, and pressure at ind_r0 (assume pressure
            % within needle is all the same). Otherwise, make a spline
            % starting at the solution for the velocity at i-ind_r0
            if i < ind_r0 && i > 1
                v_spline_inj{i,iter} = spline(t_injection,(Q/(4*pi*r_mass(i)^2))*ones(size(vp,2),1));
                dvdr_spline_inj{i,iter} = spline(t_injection,(-Q/(2*pi*r_mass(i)^3))*ones(size(vp,2),1));
                %v_spline_inj{i,iter} = spline(t_injection,(Q/(pi*r_mass(i)^2))*ones(size(vp,2),1));
                %dvdr_spline_inj{i,iter} = spline(t_injection,(-2*Q/(pi*r_mass(i)^3))*ones(size(vp,2),1));
                %v_spline_inj{i,iter} = spline(t_injection,(Q/(pi*r0^2)).*ones(size(vp,2),1));
                %dvdr_spline_inj{i,iter} = spline(t_injection,0.*ones(size(vp,2),1));
                p_spline_inj{i,iter} = spline(t_injection,pp(ind_r0,:));
            elseif i == 1
                v_spline_inj{i,iter} = spline(t_injection,(2*Q/(pi*r0^2))*ones(size(vp,2),1));
                dvdr_spline_inj{i,iter} = spline(t_injection,0.*ones(size(vp,2),1));
                p_spline_inj{i,iter} = spline(t_injection,pp(ind_r0,:));
            else            
                v_spline_inj{i,iter} = spline(t_injection,vp(i-ind_r0+1,:));
                dvdr_spline_inj{i,iter} = spline(t_injection,dvdrp(i-ind_r0+1,:));
                p_spline_inj{i,iter} = spline(t_injection,pp(i-ind_r0+1,:));
            end
        end

        save('velocity_pressure_rc_splines_inj.mat','v_spline_inj','dvdr_spline_inj','p_spline_inj','rc_spline_inj')

        toc
        
        fprintf('Mass Transport: %d\n',iter)
        tic
    
        IC = zeros(1,length(r_mass));
        IC(1:ind_r0) = c0; % assume constant intial concentration at needle tip
        [t_mp,c_mp] = ode15s(@(t_p,conc_p) getConcentration_v2(t_p, conc_p, r_mass, v_spline_inj(:,iter), dvdr_spline_inj(:,iter), p_spline_inj(:,iter), paramsIter), [t_injection(1), t_injection(end)], IC, opts1);

        if t_end > t_injection(end)
            % if we're simulating up to a time point beyond the length of
            % the injection, calculate concentration @ relaxation phase
            IC_relax = c_mp(end,:);
            for i = 1:length(r_mass)
                v_spline_rel{i,iter} = spline(t_injection,zeros(length(t_injection),1));
                dvdr_spline_rel{i,iter} = spline(t_injection,zeros(length(t_injection),1));
                p_spline_rel{i,iter} = spline(t_injection,zeros(length(t_injection),1));
            end

            save('velocity_pressure_rc_splines_rel.mat','v_spline_rel','dvdr_spline_rel','p_spline_rel','rc_spline_rel')

%             [t_mp_relax,c_mp_relax] = ode15s(@(t_p,conc_p) getConcentrationNetti(t_p, conc_p, r, v_spline, dvdr_spline, paramsIter), t_relaxation, IC_relax, opts1);
            [t_mp_relax,c_mp_relax] = ode15s(@(t_p,conc_p) getConcentration_v2(t_p, conc_p, r_mass, v_spline_rel(:,iter), dvdr_spline_rel(:,iter), p_spline_rel(:,iter), paramsIter), t_relaxation, IC_relax, opts1);
            t_m{iter} = [t_mp;t_mp_relax];
            c_m{iter} = [c_mp;c_mp_relax];
        else
            t_m{iter} = t_mp;
            c_m{iter} = c_mp;
        end
    
        toc

    end
end

if writeResults == 1 && runConcentrationSimulations == 1
    fprintf('Writing Results MT: \n')
    tic
    save(c_filename,'t_m','c_m')
    toc
end

%% PLOT PRESSURE

if plotPressureLine == 1
    if runCavitySimulations == 0
        load(p_filename);
    end
    plotResults.plotPressureLine(t_injection,p,params,variableParam)
end

if plotPressure3D == 1
    if runCavitySimulations == 0
        load(p_filename);
    end
    plotResults.plotPressure3D(t_injection,r,p,params,variableParam)
end

%% PLOT DILATATION

if plotDilatationLine == 1
    if runCavitySimulations == 0
        load(e_filename);
    end
    plotResults.plotDilatationLine(t_injection,e,params,variableParam)
end

%% PLOT SOLID DISPLACEMENT

if plotDisplacementLine == 1
    if runCavitySimulations == 0
        load(u_filename);
    end
    plotResults.plotDisplacementLine(t_injection,u,params,variableParam)
end

if plotDisplacement3D == 1
    if runCavitySimulations == 0
        load(u_filename);
    end
    plotResults.plotDisplacement3D(t_injection,r,u,params,variableParam)
end

%% PLOT FLUID VELOCITY

if plotVelocityLine == 1
    if runCavitySimulations == 0
        load(v_filename);
        load(v_spline);
    end
    t_vec = [0*60, 2*60, 5*60, 10*60, 30*60, 60*60];
    plotResults.plotVelocityLine(t_injection,v,params,variableParam)
    %plotVelocitySpline(t,t_vec,t_injection_end,r,v_splines_inj,v_splines_rel,params,variableParam
    %plotResults.plotVelocitySpline(t_m,t_vec,t_injection_end,r_mass,v_spline_inj,v_spline_rel,params,variableParam)
    %plotResults.plotSplineVsComp(t_injection,r_mass,v,v_spline_inj,params,variableParam)
    %plotResults.plotSplineVsComp(t_injection,r_mass,dvdr,dvdr_spline_inj,params,variableParam)

end

if plotVelocity3D == 1
    if runCavitySimulations == 0
        load(v_filename);
    end
    plotResults.plotVelocity3D(t_injection,r,v,params,variableParam)
end

%% PLOT CAVITY RADIUS

if plotCavityRadius == 1
    if runCavitySimulations == 0
        load(rc_filename);
    end
    plotResults.plotCavityRadius(t_injection,rc,params,variableParam)
end

%% PLOT CONCENTRATION

if plotConcentration == 1
    if runConcentrationSimulations == 0
        load(c_filename);
    end
    t_vec = [0*60, 2*60, 5*60, 10*60, 30*60, 60*60, 120*60, 240*60];
    %t_vec = [5*60, 20*60, 30*60, 60*60, 90*60];
    %t_vec = [6*60, 30*60, 60*60, 120*60];
%     plotResults.plotConcentration(t_m,t_injection,r,c_m,rc,t_vec,params,variableParam)
    plotResults.plotConcentration(t_m,t_injection,r_mass,c_m,rc,t_vec,params,variableParam)
end






