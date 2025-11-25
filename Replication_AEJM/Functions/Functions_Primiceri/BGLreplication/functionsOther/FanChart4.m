function FanChart4(IRFdens,IRFVAR,IRFDFM,TruSer,Name,t)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Produces a fan chart with two additional point forecasts and realised
% values
%
% INPUTS
% IRFdens - nPer x nRep matrx of densities
% IRFVAR - first point forecast
% IRFDFM - second point forecast
% TruSer - realised values
% Name - name of the series
% t - labels for the horizontal axis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Setting the colors
componentsColorMap = ones(9,3) ;
componentsColorMap(:,3)=0 ; componentsColorMap(1,:)= 1 ;
componentsColorMap(:,2)= [255 240 200 160 120  120 160 200 240]/255 ;

%% Quantiles of the density
M = size(IRFdens,2);
densGrid = round([0.05 0.2:0.1:0.8 0.95]*M);
YY= IRFdens(:,densGrid) ;
YY(:,2:end)= YY(:,2:end)-YY(:,1:end-1) ;

%% Plot the fan chart

h = area(t,YY, 'LineStyle','none');
hold on
grid on
colormap(componentsColorMap)
set(gca,'Layer','top')
set(h,'BaseValue',min(min(YY)))
grid

drawnow
h=plot(t,TruSer,'g-x',t,IRFVAR,'k',t,IRFDFM,'b--', 'MarkerSize',5);
set(h(2:3),'LineWidth',2) ;
set(h(1),'LineWidth',2) ;

hold on
% set the limits
set(gca,'XLim',[min(t) max(t)])
set(gca,'Fontsize',10)
set(gca,'Fontweight','bold')

hold on
xlabel([])
title(Name)
grid on

