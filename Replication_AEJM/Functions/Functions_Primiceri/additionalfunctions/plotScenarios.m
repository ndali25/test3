%% Scorecard

red  = [228,26,28]./255;
blue = [55,126,184]./255;

plotVars         = struct;
plotVars(1).type = 'Prices';
plotVars(1).Vars = {'Personal Consumption Expenditures: Chain Price Index',...
    'PCE less Food & Energy: Chain Price Index',...
    'Spot Oil Price: West Texas Intermediate [Prior''82=Posted Price]',...
    'S&P GSCI Non-Energy Commodities Nearby Index'};

plotVars(2).type = 'Economic Activity';
plotVars(2).Vars =     {'Real Gross Domestic Product', ...
    'Real Personal Consumption Expenditures', ...
    'Real Private Nonresidential Fixed Investment', ...
    'Real Exports of Goods & Services', ...
    'Real Imports of Goods & Services'};

plotVars(3).type = 'Asset Prices';
plotVars(3).Vars =     {'Stock Price Index: Standard & Poor''s 500 Composite', ...
    '10-Year Treasury Note Yield at Constant Maturity', ...
    'Nominal Trade-Weighted Exch Value of US$ vs Major Currencies',...
    'Moody''s Seasoned Aaa Corporate Bond Yield',...
    'Moody''s Seasoned Baa Corporate Bond Yield'};

plotVars(4).type = 'Interest rates';
plotVars(4).Vars =     {'2-Year Treasury Bill Yield at Constant Maturity', ...
    '5-Year Treasury Note Yield at Constant Maturity', ...
    '10-Year Treasury Note Yield at Constant Maturity'};

plotVars(5).type = 'Labor';
plotVars(5).Vars =  {'Business Sector: Compensation Per Hour', ...
    'All Employees: Total Nonfarm', ...
    'Civilian Unemployment Rate: 16 yr +'};

% xlimit indices
iPlotStart = find(h_fore == 1, 1, 'first') - 1;
iPlotEnd   = find(h_fore == 1, 1, 'last');

% Loop through subplot categories
for j = 1:size({plotVars.type}, 2)
    
    Vars = plotVars(j).Vars;  % Variables for particular category
    
    f = figure('Visible', vis);
    hold on
    set(f, 'Position', [800, 500, 1200, 600])
    
    
    
    for k = 1:length(Vars)  % Loop through all variables
        
        idx = find(strcmp(Spec.SeriesName, Vars{k}));
        
        
        if k <= 4
            
            subplot(2, 2, k)
            quantilePlot(DateAll, squeeze(dYQQ(:, idx, :)))  % Plot bands
            plot([DateAll(iPlotStart), DateAll(iPlotEnd)], [0, 0], ...
                'LineStyle', '-', 'Color', 'k')               % Reference line
            
            titleLab = Spec.SeriesName{idx};  % Set title
            if strcmp(Spec.Transformation{idx}, 'log')
                title([titleLab ' (% chg.)'])
            else
                title([titleLab ' (difference)'])
            end
            
        else  % if k == 5, put the last series in the last panel
            hold on
            
            quantilePlot(DateAll, squeeze(dYQQ(:, idx, :)), [228,26,28]./255)  % Plot bands
            
            if strcmp(plotVars(j).type, 'Asset Prices')
                title('AAA (blue) and BAA (red) Corporate Bond Yields (% chg.)')
            else
                title('Exports (blue) and Imports (red) (% chg.)')
            end
        end
        
        if ismember(idx, find(sum(~isnan(Shock)))) && size(Shock, 2) == n % Plot dot marking impulse
            
            yImpulse = dYQQ(end-h+1:end, idx, 1);
            
            yImpulse(isnan(Shock(:, idx))) = NaN;
            
            scatter(DateAll(end-h+1:end), yImpulse, 'filled', 'MarkerFaceColor', 'b')
        end
        
        
        
        xlim([DateAll(iPlotStart), DateAll(iPlotEnd)])
     %   if strcmp(plotVars(j).type, 'Interest rates')
     %       ylim([-7, 7])
      %  end
        
        grid on
        box on
        datetick('x', 'yyyy:qq', 'keeplimits')
        
    end
    
    saveas(f, [plotDir plotVars(j).type '.png'])
    
    
end

%% Conditional and unconditional forecast comparison

iPlotStart = find(h_fore == 1, 1, 'first') - 4*4;
iPlotEnd   = find(h_fore == 1, 1, 'last');

QQ = [.159, .25, .5 .75, .841];

for iVar = 1:n  % Loop through variables
    
    figure('Visible', vis)
    hold on
    
    
    subplot(2, 1, 1)
    yQQ     = squeeze(quantile(PredY_unc(:, iVar, :), QQ, 3));
    yQQ_con = squeeze(quantile(PredY_con(:, iVar, :), QQ, 3));
    
    
    quantilePlot(DateAll, yQQ)  % Quantiles
    plot(DateAll, yQQ_con, 'Color', red, 'LineStyle', ':')        % Quantiles
    
    plot(DateAll, yQQ(:, 3),     'Color', blue, 'LineWidth', 2)   % Median
    plot(DateAll, yQQ_con(:, 3), 'Color', red, 'LineWidth', 2)    % Median
    
    if size(dataTransformed, 1) < length(Date)
        plot(DateAll, dataTransformedAll(:, iVar), 'k -', 'LineWidth', 1.5)    % Data    
    else
        plot(Date, yQQ(1:length(Date), 3), 'k -', 'LineWidth', 1.5)    % Data
    end
    xlim([DateAll(iPlotStart), DateAll(iPlotEnd)])
    datetick('keepticks', 'keeplimits')
    xlabel('Year')
    box on
    grid on
    p1 = plot(NaN, NaN, 'Color', blue, 'LineWidth', 2);
    p2 = plot(NaN, NaN, 'Color', red, 'LineWidth', 2, 'LineStyle', ':');
    p3 = plot(NaN, NaN, 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', '-');
    
    legend([p1, p2, p3], {'Unc.', 'Cond.', 'Data'}, 'Location', 'best')
    
    
    if strcmp(Spec.Transformation{iVar}, 'log')
        title([Spec.SeriesName{iVar} ' (log-level)'])
        
    else
        title([Spec.SeriesName{iVar} ' (level)'])
    end
    
    subplot(2, 1, 2)
    
    quantilePlot(DateAll, squeeze(dYQQ(:, iVar, :)))  % Plot bands
    grid on
    box on
    xlabel('Year')
    xlim([DateAll(iPlotStart), DateAll(iPlotEnd)])
    plot([DateAll(iPlotStart), DateAll(iPlotEnd)], [0, 0], ...
        'LineStyle', '-', 'Color', 'k', 'LineWidth', 1.5)   % Reference line
    
    title('[Conditional] - [Unconditional]')
    datetick('keeplimits', 'keepticks')
    saveas(gcf, [plotDir 'subplot_' Spec.SeriesID{iVar} '.png'])
    close all
end


