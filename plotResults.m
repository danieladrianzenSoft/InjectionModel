classdef plotResults
    methods(Static)
        function plotPressureLine(t,p,params,variableParam)
            variableParamArray = params(variableParam);
            lengthVariableParam = length(variableParamArray{1});
            f = figure;
            ax = axes("Parent",f);
            legendLabels = cell(lengthVariableParam,1);
            for iter = 1:lengthVariableParam
                pPlot = p{iter};
                plot(ax,t/60,pPlot(1,:)/(1e4), "LineWidth",5)
                hold on
                legendLabels{iter} = sprintf('%s = %.0e',variableParam, variableParamArray{1}(iter));
            end
            set(ax,"FontSize", 28)
            xlabel(ax, "Time (mins)", "FontSize", 32)
            ylabel(ax, "Pressure (kPa)", "FontSize", 32)
            xlim(ax,[t(1)/60, t(end)/60])
            legend(ax, legendLabels)
        end
        function plotPressure3D(t,r,p,params,variableParam)
            variableParamArray = params(variableParam);
            lengthVariableParam = length(variableParamArray{1});
            for iter = 1:lengthVariableParam
                f = figure;
                ax = axes("Parent",f);
                pPlot = p{iter};
                [X,Y] = meshgrid(t(1:6:end),r(1:10:end));
                Z = pPlot(1:10:end,1:6:end)/1e4;
                s = surf(X/60,Y,Z,'FaceAlpha',0.5);
                set(ax,"FontSize", 28)
                ylabel(ax, "r (cm)", "FontSize", 32)
                xlabel(ax, "Time (mins)", "FontSize", 32)
                zlabel(ax, "Pressure (kPa)","FontSize",32)
                view(-242,22)
                titleText = sprintf('%s = %.0e m^2/kPa',variableParam, variableParamArray{1}(iter));
                title(titleText,"FontSize",34)
                %zlim([0,1.4])
            end
        end
        function plotDilatationLine(t,e,params,variableParam)
            variableParamArray = params(variableParam);
            lengthVariableParam = length(variableParamArray{1});
            f = figure;
            ax = axes("Parent",f);
            legendLabels = cell(lengthVariableParam,1);
            for iter = 1:lengthVariableParam
                ePlot = e{iter};
                plot(ax,t/60,ePlot(1,:), "LineWidth",5)
                hold on
                legendLabels{iter} = sprintf('%s = %.0e',variableParam, variableParamArray{1}(iter));
            end
            set(ax,"FontSize", 28)
            xlabel(ax, "Time (mins)", "FontSize", 32)
            ylabel(ax, "Dilatation", "FontSize", 32)
            xlim(ax,[t(1)/60, t(end)/60])
            legend(ax, legendLabels)
        end
        function plotDisplacementLine(t,u,params,variableParam)
            variableParamArray = params(variableParam);
            lengthVariableParam = length(variableParamArray{1});
            f = figure;
            ax = axes("Parent",f);
            legendLabels = cell(lengthVariableParam,1);
            for iter = 1:lengthVariableParam
                uPlot = u{iter};
                plot(ax,t/60,uPlot(1,:), "LineWidth",5)
                hold on
                legendLabels{iter} = sprintf('%s = %.0e',variableParam, variableParamArray{1}(iter));
            end
            set(ax,"FontSize", 28)
            xlabel(ax, "Time (mins)", "FontSize", 32)
            ylabel(ax, "Solid displacement (cm)", "FontSize", 32)
            xlim(ax,[t(1)/60, t(end)/60])
            legend(ax, legendLabels)
        end
        function plotDisplacement3D(t,r,u,params,variableParam)
            variableParamArray = params(variableParam);
            lengthVariableParam = length(variableParamArray{1});
            for iter = 1:lengthVariableParam
                f = figure;
                ax = axes("Parent",f);
                uPlot = u{iter};
                [X,Y] = meshgrid(t(1:6:end),r(1:10:end));
                Z = uPlot(1:10:end,1:6:end)/1e4;
                s = surf(X/60,Y,Z,'FaceAlpha',0.5);
                set(ax,"FontSize", 28)
                ylabel(ax, "r (cm)", "FontSize", 32)
                xlabel(ax, "Time (mins)", "FontSize", 32)
                zlabel(ax, "Solid Displacement (cm)","FontSize",32)
                view(-242,22)
                titleText = sprintf('%s = %.0e m^2/kPa',variableParam, variableParamArray{1}(iter));
                title(titleText,"FontSize",34)
                %zlim([0,2e-7])
            end
        end
        function plotVelocityLine(t,v,params,variableParam)
            variableParamArray = params(variableParam);
            lengthVariableParam = length(variableParamArray{1});
            f = figure;
            ax = axes("Parent",f);
            legendLabels = cell(lengthVariableParam,1);
            for iter = 1:lengthVariableParam
                vPlot = v{iter};
                plot(ax,t/60,vPlot(1,:), "LineWidth",2)
                %plot(ax,r,vPlot(:,1), "LineWidth",2)
                hold on
                legendLabels{iter} = sprintf('%s = %.0e',variableParam, variableParamArray{1}(iter));
            end
            set(ax,"FontSize", 28)
            xlabel(ax, "Time (mins)", "FontSize", 32)
            %xlabel(ax, "Radius (cm)", "FontSize", 32) 
            ylabel(ax, "Velocity (cm/s)", "FontSize", 32)
            xlim(ax,[t(1)/60, t(end)/60])
            %xlim(ax,[r(1),r(end)])
            legend(ax, legendLabels)
        end
        function plotSplineVsComp(t,r,y,spline_inj,params,variableParam,options)
            arguments
               t
               r (1,:) {mustBeNumeric}
               y
               spline_inj
               params
               variableParam
               options.xLabel = 'Time (mins)'
               options.yLabel = 'Fluid Velocity (cm/s)'
            end
            variableParamArray = params(variableParam);
            lengthVariableParam = length(variableParamArray{1});
            for p_iter = 1:lengthVariableParam
                r0_iter = params("r0");
                r0 = r0_iter{p_iter};
                rVec = [find(r >= r0,1), find(r >= 0.1,1), find(r >= 0.2,1), find(r >= 0.5,1), find(r >= 1,1)];
                f = figure();
                ax = axes("Parent",f);
                %t_p = t{p_iter};
                %tPlot = t_p(t_p<=t_injection_end);
                y_spline = spline_inj(:,p_iter);
                y_plot = y{p_iter};
                legendLabels = cell(length(rVec),1);

                for r_iter = 1:length(rVec)
                    y_spline_plot = y_spline{rVec(r_iter)};
                    plot(ax,t/60,y_plot(rVec(r_iter)-rVec(1)+1,:),'LineWidth',5)
                    hold on
                    plot(ax,t/60,ppval(t,y_spline_plot),'kx')
                    legendLabels{2*r_iter-1} = sprintf('r = %.2f cm', r(rVec(r_iter)));
                    legendLabels{2*r_iter} = sprintf('r = %.2f cm - fit', r(rVec(r_iter)));
                    %yPlot = y{r_iter};
                    %plot(ax,t/60,y_plot(rVec(r_iter),:),'LineWidth',2)
                end
                ylabel(ax, options.yLabel,'FontSize',32)
                xlabel(ax, options.xLabel,'FontSize',32)
                set(ax,'FontSize',28)
                legend(legendLabels)
                xlim(ax,[0,t(end)/60])
                title(sprintf('%s = %.1e',variableParam, variableParamArray{1}(p_iter)),"FontSize", 32)
            end
        end
        function plotVelocitySpline(t,t_vec,t_injection_end,r,v_splines_inj,v_splines_rel,params,variableParam)
            %% FINISH EDITING THIS
            variableParamArray = params(variableParam);
            lengthVariableParam = length(variableParamArray{1});
            %f = figure;
            %ax = axes("Parent",f);
            %legendLabels = cell(lengthVariableParam,1);
            %plotResults.plotVelocitySpline(t_m,t_injection,t_vec,r,v_spline_inj,t_vec,params,variableParam)

            for iter = 1:lengthVariableParam
                tPlot = t{iter};
                f = figure;
                ax = axes("Parent",f);

                legendLabels = cell(length(t_vec),1);
                for t_iter=1:length(t_vec)
                    t_ind = find(tPlot >= t_vec(t_iter), 1);
                    if tPlot(t_ind) > t_injection_end
                        v_spline = v_splines_inj(:,iter);
                    else
                        v_spline = v_splines_rel(:,iter);
                    end

                    vPlot = zeros(length(r),1);
                    for r_iter = 1:length(r)
                        vPlot(r_iter) = ppval(tPlot(t_iter),v_spline{r_iter});
                    end

                    legendLabels{t_iter} = sprintf('%d mins', t_vec(t_iter)/60);
                    plot(r,vPlot,'LineWidth',5);
                    hold on
                end
                
                ylabel(ax, 'Fluid velocity (cm/s)','FontSize',32)
                xlabel(ax, 'Radius (cm)','FontSize',32)
                set(ax,'FontSize',28)
                legend(legendLabels)
                xlim(ax,[0,r(end)])
                title(sprintf('%s = %.1e',variableParam, variableParamArray{1}(iter)),"FontSize", 32)
                %plot(ax,tPlot/60,vPlot(1,:), "LineWidth",2)
                %plot(ax,r,vPlot(:,1), "LineWidth",2)
                %hold on
                %legendLabels{iter} = sprintf('%s = %.0e',variableParam, variableParamArray{1}(iter));
            end
            set(ax,"FontSize", 28)
            %xlabel(ax, "Time (mins)", "FontSize", 32)
            xlabel(ax, "Radius (cm)", "FontSize", 32)
            ylabel(ax, "Velocity (cm/s)", "FontSize", 32)
            %xlim(ax,[t(1)/60, t(end)/60])
            %xlim(ax,[r(1),r(end)])
            legend(ax, legendLabels)

            for p_iter = 1:lengthVariableParam
                r0_iter = params("r0");
                r0 = r0_iter{p_iter};
                rVec = [find(r >= r0,1), find(r >= 0.1,1), find(r >= 0.2,1), find(r >= 0.5,1), find(r >= 1,1)];
                f = figure();
                ax = axes("Parent",f);
                tPlot = t{p_iter};
                tPlot = tPlot(tPlot<=t_injection_end);
                %v = v_spline_inj()
                v_spline = v_splines_inj(:,p_iter);
                legendLabels = cell(length(rVec),1);

                for r_iter = 1:length(rVec)
                    v = v_spline{r_iter};
                    plot(ax,tPlot/60,ppval(tPlot,v),'LineWidth',2)
                    legendLabels{r_iter} = sprintf('%d cm', r(rVec(r_iter)));
                    hold on
                end
                ylabel(ax, 'Fluid velocity (cm/s)','FontSize',32)
                xlabel(ax, 'Time (mins)','FontSize',32)
                set(ax,'FontSize',28)
                legend(legendLabels)
                xlim(ax,[0,tPlot(end)/60])
                title(sprintf('%s = %.1e',variableParam, variableParamArray{1}(p_iter)),"FontSize", 32)
            end
        end
        function plotVelocity3D(t,r,v,params,variableParam)
            variableParamArray = params(variableParam);
            lengthVariableParam = length(variableParamArray{1});
            for iter = 1:lengthVariableParam
                f = figure;
                ax = axes("Parent",f);
                vPlot = v{iter};
                [X,Y] = meshgrid(t(1:6:end),r(1:10:end));
                Z = vPlot(1:10:end,1:6:end)/1e4;
                s = surf(X/60,Y,Z,'FaceAlpha',0.5);
                set(ax,"FontSize", 28)
                ylabel(ax, "r (cm)", "FontSize", 32)
                xlabel(ax, "Time (mins)", "FontSize", 32)
                zlabel(ax, "Velocity (cm/s)","FontSize",32)
                view(-242,22)
                titleText = sprintf('%s = %.0e m^2/kPa',variableParam, variableParamArray{1}(iter));
                title(titleText,"FontSize",34)
                %zlim([0,1.4])
            end
        end
        function plotCavityRadius(t,rc,params,variableParam)
            variableParamArray = params(variableParam);
            lengthVariableParam = length(variableParamArray{1});
            f = figure;
            ax = axes("Parent",f);
            r0 = params('r0');
            legendLabels = cell(lengthVariableParam,1);
            maxRc = 0;
            for iter = 1:lengthVariableParam
                rcPlot = rc{iter};
                maxRc = max([maxRc,max(rcPlot)]);
                plot(ax,t/60,rcPlot, "LineWidth",5)
                hold on
                legendLabels{iter} = sprintf('%s = %.0e',variableParam, variableParamArray{1}(iter));
            end
            set(ax,"FontSize", 28)
            xlabel(ax, "Time (mins)", "FontSize", 32)
            ylabel(ax, "Cavity Radius (cm)", "FontSize", 32)
            
            xlim(ax,[t(1)/60, t(end)/60])
            ylim(ax,[r0{1}, maxRc])
            legend(ax, legendLabels)
        end
        function plotConcentration(t_m,t,r,c_m,rc,t_vec,params,variableParam)
            variableParamArray = params(variableParam);
            lengthVariableParam = length(variableParamArray{1});

            for iter = 1:lengthVariableParam
                rcp = rc{iter};
                %rc_spline = cell(length(r),1);
                %for j = 1:length(r)
                %    rc_spline{iter} = spline(t,rcp(j,:));
                %end
                rc_spline = spline(t,rcp);
                f = figure;
                ax = axes("Parent",f);
        
                cPlot = c_m{iter};
                tPlot = t_m{iter};
        
                %t_inds = zeros(length(t_vec),1);
                legendLabels = cell(length(t_vec),1);
                for t_iter=1:length(t_vec)
                    t_ind = find(tPlot >= t_vec(t_iter), 1);
                    rc_plot = ppval(rc_spline,t_vec(t_iter));
                    legendLabels{t_iter} = sprintf('%d mins', t_vec(t_iter)/60);
                    %plot(r-rc_plot,cPlot(t_ind,:),'LineWidth',5);
                    plot(r,cPlot(t_ind,:),'LineWidth',5);
                    hold on
                end
                %title(ax, titleLabel, 'FontSize', 32, 'FontWeight', 'Bold')
                ylabel(ax, 'Ethanol Concentration (g/ml)','FontSize',32)
                xlabel(ax, 'Radius (cm)','FontSize',32)
                set(ax,'FontSize',28)
                legend(legendLabels)
                xlim(ax,[0,r(end)])
                title(sprintf('%s = %.1e',variableParam, variableParamArray{1}(iter)),"FontSize", 32)
            end
        end
    end
end
