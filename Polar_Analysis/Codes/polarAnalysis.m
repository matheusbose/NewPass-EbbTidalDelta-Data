function polarAnalysis(datafile,origin,outPathRoot,dTheta,dRho,rhoLimOuter,thetaLims,testMultiOrigins,interest)
% polarAnalysis_clean
% Versão simplificada para:
%   (1) plotar vAnomaly em XY e rho-theta
%   (2) plotar Hovmoller de vAnomaly
%
% Mantém apenas o necessário para essas figuras.

%% ========================= USER SETTINGS =============================
fontName = 'Arial';
fontSize = 12;

figPos = [2 2 21 17];   % [left bottom width height] em cm
figRes = 300;

% escolha os anos que quer mostrar no eixo do Hovmoller
yearsToShow = [];   % ex.: [2004 2006 2010 2015 2022]
                    % deixe [] para usar automaticamente os anos de T

% máscara vertical dos dados
zMax = 0;
zMin = -5;

%% =========================== LOAD DATA ===============================
if exist(outPathRoot,'dir') == 7
    addpath(outPathRoot);
else
    mkdir(outPathRoot)
    addpath(outPathRoot)
end

load(datafile,'grd');

T = grd.year;
if iscolumn(T); T = T'; end

X = grd.x / 1000;   % km
Y = grd.y / 1000;   % km
Z = grd.dp;         % depth negative
Z(Z > zMax | Z < zMin) = NaN;

if isempty(yearsToShow)
    yearsToShow = T;
end

o0 = origin / 1000; % km
XLims = round([min(X,[],'all') max(X,[],'all')]);

%% ======================= ORIGIN DEFINITION ===========================
if testMultiOrigins
    testDist = round(XLims(2) - o0(1)) * 0.1;

    origins = [o0(1),o0(2);...
        o0(1),o0(2)+testDist;...
        o0(1)+testDist,o0(2);...
        o0(1),o0(2)-testDist;...
        o0(1)-testDist,o0(2);...
        o0(1)+testDist,o0(2)+testDist;...
        o0(1)+testDist,o0(2)-testDist;...
        o0(1)-testDist,o0(2)+testDist;...
        o0(1)-testDist,o0(2)-testDist;...
        o0(1),o0(2)+testDist*2;...
        o0(1)+testDist*2,o0(2);...
        o0(1),o0(2)-testDist*2;...
        o0(1)-testDist*2,o0(2)];

    originName = {'O';'N1';'E1';'S1';'W1';...
        'NE1';'SE1';'NW1';'SW1';...
        'N2';'E2';'S2';'W2'};

    f = newFigure([],figPos);
    hold on
    imagescn(X,Y,Z(:,:,1));
    colormap('gray')

    for i = 1:size(origins,1)
        scatter(origins(i,1),origins(i,2),'filled');
    end

    title('Origin Test Locations')
    xlabel('X (km)')
    ylabel('Y (km)')
    grid on
    legend(originName)

    applyFigureStyle(gcf,fontName,fontSize)

    pngName = 'origin_test_locations';
    savefig(fullfile(outPathRoot,[pngName '.fig']))
    exportgraphics(gcf,fullfile(outPathRoot,[pngName '.png']),'Resolution',figRes)
else
    origins = o0;
    originName = {'O'};
end

%% ===================== LOOP OVER ORIGINS =============================
for oo = 1:size(origins,1)

    if testMultiOrigins
        disp("Now analyzing origin " + originName(oo));
    end

    x0 = origins(oo,1);
    y0 = origins(oo,2);

    outPath = fullfile(outPathRoot,originName{oo});
    if ~exist(outPath,'dir')
        mkdir(outPath)
    end
    addpath(outPath)

    %% ================= WEIGHTED MEAN SURFACE ==========================
    iT0 = 2;

    Tmidpt = NaN(1,length(T)-1);
    Tmidpt(1) = T(iT0);
    ttMid = iT0:length(T)-1;
    Tmidpt(2:end) = (T(ttMid)+T(ttMid+1))/2;
    Tmidpt(length(T)) = (T(end)+T(end)+1)/2;
    Tweight = diff(Tmidpt);

    weightedMeanSurf = zeros(size(X));
    for tt = iT0:length(T)
        weightedMeanSurf = weightedMeanSurf + squeeze(Z(:,:,tt)) .* Tweight(tt-1);
    end
    weightedMeanSurf = weightedMeanSurf ./ (T(end)-T(iT0));

    %% ======================= CELL AREA ================================
    dx = unique(diff(X,1,2))*1000;
    dy = unique(diff(Y,1,1))*1000;
    areaGridCell = floor(abs(dx(1)*dy(1)));

    %% ===================== ANOMALY VOLUME =============================
    hAnomaly = nan(size(Z));
    hAnomaly(:,:,1:length(T)) = Z(:,:,1:length(T)) - weightedMeanSurf;
    vAnomaly = hAnomaly .* areaGridCell;

    %% ====================== POLAR GRID ================================
    rhoLims = [round(rhoLimOuter*interest) rhoLimOuter];
    theta = (thetaLims(1):dTheta:thetaLims(2)).*(2*pi/360);
    rho = rhoLims(1):dRho:rhoLims(2);

    [xPolar,yPolar] = pol2cart(repmat(theta,length(rho),1), ...
                               repmat(rho,length(theta),1)');
    xPolar = xPolar./1000;
    yPolar = yPolar./1000;

    xq = reshape(X,[],1);
    yq = reshape(Y,[],1);

    zPolar        = nan(length(rho)-1,length(theta)-1,length(T));
    vAnomalyPolar = nan(length(rho)-1,length(theta)-1,length(T));
    meanSurfPolar = nan(length(rho)-1,length(theta)-1);

    %% =================== BIN DATA TO POLAR CELLS =====================
    for tt = 1:length(T)
        disp(['Processing year: ' num2str(T(tt))])

        zq = reshape(Z(:,:,tt),[],1);
        vq = reshape(vAnomaly(:,:,tt),[],1);

        for jj = 1:length(theta)-1
            for ii = 1:length(rho)-1

                xv = [xPolar(ii,jj)   xPolar(ii+1,jj)   xPolar(ii+1,jj+1) ...
                      xPolar(ii,jj+1) xPolar(ii,jj)] + x0;

                yv = [yPolar(ii,jj)   yPolar(ii+1,jj)   yPolar(ii+1,jj+1) ...
                      yPolar(ii,jj+1) yPolar(ii,jj)] + y0;

                in = inpolygon(xq,yq,xv,yv);

                zPolar(ii,jj,tt)        = mean(zq(in),'omitnan');
                vAnomalyPolar(ii,jj,tt) = sum(vq(in));
            end
        end
    end
    
        %% =================== BIN WEIGHTED MEAN SURFACE ==================
    wq = reshape(weightedMeanSurf,[],1);

    for jj = 1:length(theta)-1
        for ii = 1:length(rho)-1

            xv = [xPolar(ii,jj)   xPolar(ii+1,jj)   xPolar(ii+1,jj+1) ...
                  xPolar(ii,jj+1) xPolar(ii,jj)] + x0;

            yv = [yPolar(ii,jj)   yPolar(ii+1,jj)   yPolar(ii+1,jj+1) ...
                  yPolar(ii,jj+1) yPolar(ii,jj)] + y0;

            in = inpolygon(xq,yq,xv,yv);

            meanSurfPolar(ii,jj) = mean(wq(in),'omitnan');
        end
    end
    %% ==================== COLLAPSED ANOMALIES ========================
    vAnomalyRho   = squeeze(sum(vAnomalyPolar,2,'omitnan')); % [rho x time]
    vAnomalyTheta = squeeze(sum(vAnomalyPolar,1,'omitnan')); % [theta x time]
    

    %% ================= PLOT GRID / AXIS HELPERS ======================
    dThetaPlotGrid = 22.5;
    if dTheta < 0
        dThetaPlotGrid = -dThetaPlotGrid;
    end

    thetaPlotGrid = (thetaLims(1):dThetaPlotGrid:thetaLims(2)).*(2*pi/360);
    thetaPlotGridFine = (thetaLims(1):dThetaPlotGrid/2:thetaLims(2)).*(2*pi/360);
    rhoPlotGrid = linspace(0,rhoLims(2),4);

    [xPolarPlotGrid,yPolarPlotGrid] = pol2cart(repmat(thetaPlotGrid,length(rhoPlotGrid),1), ...
                                               repmat(rhoPlotGrid,length(thetaPlotGrid),1)');
    xPolarPlotGrid = xPolarPlotGrid./1000;
    yPolarPlotGrid = yPolarPlotGrid./1000;

    [xPolarPlotGridFine,yPolarPlotGridFine] = pol2cart(repmat(thetaPlotGridFine,length(rhoPlotGrid),1), ...
                                                       repmat(rhoPlotGrid,length(thetaPlotGridFine),1)');
    xPolarPlotGridFine = xPolarPlotGridFine./1000;
    yPolarPlotGridFine = yPolarPlotGridFine./1000;

    thetaDegNorth = theta./(2*pi/360);
    thetaDegNorth = 90 - thetaDegNorth;  % nautical convention
    [thetaGrid,rhoGrid] = meshgrid(thetaDegNorth(1:end-1),rho(1:end-1)./1000);

    %% ===================== SAVE CORE DATA ============================
    save(fullfile(outPath,[originName{oo} '_polarData.mat']), ...
        'T','X','Y','Z','weightedMeanSurf','meanSurfPolar', ...
        'theta','rho','thetaDegNorth','thetaGrid','rhoGrid', ...
        'xPolar','yPolar','zPolar','vAnomalyPolar','vAnomalyTheta','vAnomalyRho', ...
        'xPolarPlotGrid','yPolarPlotGrid','xPolarPlotGridFine','yPolarPlotGridFine', ...
        'rhoLims','thetaLims');

    %% ==================== COLORMAPS / LABELS =========================
    % -------- bathy colormap (discrete) --------
    edges = [-9.8 -8 -7 -6 -5 -4 -3 -2 -1 0 1 2 3];
    cmapBath = [
        0.50 0.00 0.00
        0.80 0.00 0.00
        1.00 0.50 0.20
        1.00 0.70 0.40
        1.00 0.88 0.55
        1.00 0.97 0.80
        0.80 0.92 1.00
        0.60 0.80 1.00
        0.40 0.60 0.85
        0.25 0.45 0.85
        0.10 0.25 0.55
        0.00 0.08 0.40
    ];
    cmapBath = flipud(cmapBath);

    tickCenters = (edges(1:end-1)+edges(2:end))/2;
    tickLabels = {'-9.8–-8','-7.9–-7','-6.9–-6','-5.9–-5','-4.9–-4','-3.9–-3', ...
                  '-2.9–-2','-1.9–-1','-0.9–0','0.1–1','1.1–2','2.1–3'};

    % -------- hovmoller colormap --------
    v = [-1 -0.9 -0.7 -0.5 -0.3 -0.15 -0.07 -0.03 0 0.03 0.07 0.15 0.3 0.5 0.7 0.9 1]';
    cols = [ ...
        0.00 0.08 0.40;
        0.10 0.25 0.55;
        0.25 0.45 0.85;
        0.40 0.60 0.85;
        0.60 0.80 1.00;
        0.80 0.92 1.00;
        0.90 0.95 1.00;
        0.97 0.98 1.00;
        1.00 1.00 1.00;
        1.00 0.98 0.97;
        1.00 0.95 0.90;
        1.00 0.88 0.55;
        1.00 0.70 0.40;
        1.00 0.50 0.20;
        0.80 0.00 0.00;
        0.50 0.00 0.00;
        0.40 0.00 0.00];
    xi = linspace(-1,1,256)';
    cmapHov = interp1(v, cols, xi, 'pchip');

    %% ================================================================
    %% FIGURE 1 - BATHYMETRY OVERVIEW (XY and rho-theta)
    %% ================================================================
    for tt = 1:length(T)

        f = newFigure(1,figPos);
        tlo = tiledlayout(1,2,'TileSpacing','tight','Padding','tight');

        % -------------------- XY --------------------
        ax1 = nexttile;
        hold(ax1,'on'); grid(ax1,'on'); box(ax1,'on')

        xp = xPolar(1:end-1,1:end-1);
        yp = yPolar(1:end-1,1:end-1);

        contourf(ax1,xp,yp,zPolar(:,:,tt),edges,'LineStyle','none');
        colormap(ax1,cmapBath)
        caxis(ax1,[edges(1) edges(end)])

        contour(ax1,xp,yp,zPolar(:,:,tt),[0 0],'k-','LineWidth',1.2)
        contour(ax1,xp,yp,zPolar(:,:,tt),[-5 0],'k--','LineWidth',1.2)
        for ii = 1:length(thetaPlotGrid)
            plot(ax1,xPolarPlotGrid(:,ii),yPolarPlotGrid(:,ii),'-k')
        end
        for ii = 1:length(rhoPlotGrid)
            plot(ax1,xPolarPlotGridFine(ii,:),yPolarPlotGridFine(ii,:),'-k')
        end

        xlabel(ax1,'X (km)')
        ylabel(ax1,'Y (km)')
        title(ax1,['X-Y Space [' num2str(T(tt)) ']'])
        axis(ax1,'square')
        set(ax1,'Layer','top','GridColor','k','GridAlpha',0.35)
        % -------------------- RHO-THETA --------------------
        ax2 = nexttile;
        hold(ax2,'on'); grid(ax2,'on'); box(ax2,'on')

        contourf(ax2,thetaGrid,rhoGrid,zPolar(:,:,tt),edges,'LineStyle','none');
        colormap(ax2,cmapBath)
        caxis(ax2,[edges(1) edges(end)])

        contour(ax2,thetaGrid,rhoGrid,zPolar(:,:,tt),[0 0],'k-','LineWidth',1.2)
        contour(ax2,thetaGrid,rhoGrid,zPolar(:,:,tt),[-5 0],'k--','LineWidth',1.2)

        xlabel(ax2,'\theta (\circ)')
        ylabel(ax2,'\rho (km)')
        title(ax2,['\rho-\theta Space [' num2str(T(tt)) ']'])

        set(ax2,'XTick',min(thetaDegNorth):30:max(thetaDegNorth))
        xlim(ax2,[min(thetaDegNorth) max(thetaDegNorth)])
        ylim(ax2,[min(rho./1000) max(rho./1000)])
        axis(ax2,'square')
        set(ax2,'Layer','top','GridColor','k','GridAlpha',0.35)
        cb = colorbar(ax2,'Location','eastoutside');
        cb.Ticks = tickCenters;
        cb.TickLabels = tickLabels;
        ylabel(cb,'Elevation (m)')

        applyFigureStyle(f,fontName,fontSize)

        pngName = [originName{oo} '_polarMap_Overview_' num2str(T(tt))];
        savefig(fullfile(outPath,[pngName '.fig']))
        exportgraphics(f,fullfile(outPath,[pngName '.png']),'Resolution',figRes)
    end

    %% ================================================================
    %% FIGURE 3 - vAnomaly in XY and rho-theta
    %% ================================================================
    for tt = 1:length(T)

        f = newFigure(3,figPos);
        tlo = tiledlayout(1,2,'TileSpacing','tight','Padding','tight');

        minC = min(vAnomalyPolar,[],'all');
        maxC = max(vAnomalyPolar,[],'all');
        %climSym = [-maxC maxC];
        climSym = [-maxC/1.5 maxC/1.5];

        % -------------------- XY --------------------
        ax1 = nexttile;
        %ax1.Color = [0.97 0.97 0.97];

        hold(ax1,'on'); grid(ax1,'on'); box(ax1,'on')
        
        xp = xPolar(1:end-1,1:end-1);
        yp = yPolar(1:end-1,1:end-1);

        pc = pcolor(ax1,xp,yp,vAnomalyPolar(:,:,tt));
        pc.FaceColor = 'flat';
        pc.EdgeColor = 'none';

        for ii = 1:length(thetaPlotGrid)
            plot(ax1,xPolarPlotGrid(:,ii),yPolarPlotGrid(:,ii),'-k')
        end
        for ii = 1:length(rhoPlotGrid)
            plot(ax1,xPolarPlotGridFine(ii,:),yPolarPlotGridFine(ii,:),'-k')
        end

        colormap(ax1,cmapHov)
        clim(ax1,climSym)

        xlabel(ax1,'X (km)')
        ylabel(ax1,'Y (km)')
        title(ax1,['X-Y Space [' num2str(T(tt)) ']'])
        xlim(ax1,[min(xPolar(:)) max(xPolar(:))])
        ylim(ax1,[min(yPolar(:)) max(yPolar(:))])

        axis(ax1,'square')
        set(ax1,'Layer','top','GridColor','k','GridAlpha',0.35)
        % -------------------- RHO-THETA --------------------
        ax2 = nexttile;
        %ax2.Color = [0.97 0.97 0.97];

        hold(ax2,'on'); grid(ax2,'on'); box(ax2,'on')

        pc = pcolor(ax2,thetaGrid,rhoGrid,vAnomalyPolar(:,:,tt));
        pc.FaceColor = 'flat';
        pc.EdgeColor = 'none';

        colormap(ax2,cmapHov)
        clim(ax2,climSym)

        xlabel(ax2,'\theta (\circ)')
        ylabel(ax2,'\rho (km)')
        title(ax2,['\rho-\theta Space [' num2str(T(tt)) ']'])

        set(ax2,'XTick',min(thetaDegNorth):30:max(thetaDegNorth))
        xlim(ax2,[min(thetaDegNorth) max(thetaDegNorth)])
        ylim(ax2,[min(rho./1000) max(rho./1000)])
        axis(ax2,'square')
        set(ax2,'Layer','top','GridColor','k','GridAlpha',0.35)
        cb = colorbar(ax2,'Location','eastoutside');
        cb.Ruler.Exponent = 4;
        ylabel(cb,'V_{Anomaly} (m^3)')

        applyFigureStyle(f,fontName,fontSize)

        pngName = [originName{oo} '_polarMap_vAnomaly_' num2str(T(tt))];
        savefig(fullfile(outPath,[pngName '.fig']))
        exportgraphics(f,fullfile(outPath,[pngName '.png']),'Resolution',figRes)
    end

%% ================================================================
%% FIGURE 4 - Hovmoller (estilo fig3: anos intermediários em branco)
%% ================================================================
Tdata = T(:)';                 % anos com dado
Tfull = Tdata(1):Tdata(end);   % anos contínuos
T_hov = [Tfull Tfull(end)+1];  % borda extra para o pcolor

% se quiser escolher manualmente os ticks do eixo:
% yearsToShow = [2004 2006 2010 2015 2022];
if isempty(yearsToShow)
    yearsToShow = Tdata;
end

% monta matrizes com NaN nos anos sem dado
vTheta_full = NaN(size(vAnomalyTheta,1), numel(Tfull));   % [theta x time]
vRho_full   = NaN(size(vAnomalyRho,1),   numel(Tfull));   % [rho   x time]

for k = 1:numel(Tdata)
    idx = find(Tfull == Tdata(k));
    if ~isempty(idx)
        vTheta_full(:,idx) = vAnomalyTheta(:,k);
        vRho_full(:,idx)   = vAnomalyRho(:,k);
    end
end

cmax = max(abs([vTheta_full(:); vRho_full(:)]));

f = newFigure(4,figPos);
tlo = tiledlayout(1,2,'TileSpacing','tight','Padding','tight');

% -------------------- THETA --------------------
ax1 = nexttile;
hold(ax1,'on'); box(ax1,'on'); grid(ax1,'on')

pp = pcolor(ax1, ...
    repmat(thetaDegNorth(1:end-1)',1,length(T_hov)), ...
    repmat(T_hov,length(thetaDegNorth(1:end-1)),1), ...
    [vTheta_full NaN(size(vTheta_full,1),1)]);
pp.EdgeColor = 'none';

colormap(ax1,cmapHov)
caxis(ax1,[-cmax/1.5 cmax/1.5])

xlabel(ax1,'\theta (\circ)')
ylabel(ax1,'Year')
title(ax1,'Changes in Volume Anomaly [\theta]')

xlim(ax1,[min(thetaDegNorth) max(thetaDegNorth)])
ylim(ax1,[T_hov(1) T_hov(end)])

set(ax1,'XTick',min(thetaDegNorth):30:max(thetaDegNorth))
set(ax1,'YTick',Tfull)   % ou yearsToShow, se quiser menos anos
set(ax1,'Layer','top','GridColor','k','GridAlpha',0.35)
axis(ax1,'square')

% -------------------- RHO --------------------
ax2 = nexttile;
hold(ax2,'on'); box(ax2,'on'); grid(ax2,'on')

pp = pcolor(ax2, ...
    repmat(T_hov,length(rho(1:end-1)),1), ...
    repmat(rho(1:end-1)./1000,length(T_hov),1)', ...
    [vRho_full NaN(size(vRho_full,1),1)]);
pp.EdgeColor = 'none';

colormap(ax2,cmapHov)
caxis(ax2,[-cmax/1.5 cmax/1.5])


xlabel(ax2,'Year')
ylabel(ax2,'\rho (km)')
title(ax2,'Changes in Volume Anomaly (\rho)')

xlim(ax2,[T_hov(1) T_hov(end)])
ylim(ax2,[rhoLims(1) rhoLims(2)]./1000)

set(ax2,'XTick',Tfull)   % ou yearsToShow
set(ax2,'Layer','top','GridColor','k','GridAlpha',0.35)
axis(ax2,'square')

cb = colorbar(ax2,'Location','eastoutside');
ylabel(cb,'V_{anomaly} (m^3)')

applyFigureStyle(f,fontName,fontSize)

pngName = [originName{oo} '_YEAR_HovmullerPolarRhoAndTheta'];
savefig(fullfile(outPath,[pngName '.fig']))
exportgraphics(f,fullfile(outPath,[pngName '.png']),'Resolution',figRes)

%% ================================================================
%% FIGURE 4B - Hovmoller in theta separated by rho regions
%% Inner region: rho <= 1 km
%% Outer region: rho > 1 km
%% ================================================================

% ------------------------------------------------
% Define rho regions
% ------------------------------------------------
rhoKm = rho(1:end-1)./1000;   % same rho used in vAnomalyPolar

idxInner = rhoKm <= 0.8;
idxOuter = rhoKm > 0.8;

% ------------------------------------------------
% Collapse vAnomalyPolar in theta for each rho region
% Result dimensions: [theta x time]
% ------------------------------------------------
vAnomalyTheta_inner = squeeze(sum(vAnomalyPolar(idxInner,:,:),1,'omitnan'));
vAnomalyTheta_outer = squeeze(sum(vAnomalyPolar(idxOuter,:,:),1,'omitnan'));

% make sure output stays 2D in case one dimension collapses oddly
vAnomalyTheta_inner = reshape(vAnomalyTheta_inner,length(thetaDegNorth)-1,length(T));
vAnomalyTheta_outer = reshape(vAnomalyTheta_outer,length(thetaDegNorth)-1,length(T));

% ------------------------------------------------
% Build continuous yearly matrices with NaN for missing years
% ------------------------------------------------
Tdata = T(:)';                 
Tfull = Tdata(1):Tdata(end);   
T_hov = [Tfull Tfull(end)+1];  

if isempty(yearsToShow)
    yearsToShow = Tdata;
end

vTheta_inner_full = NaN(size(vAnomalyTheta_inner,1), numel(Tfull));   % [theta x time]
vTheta_outer_full = NaN(size(vAnomalyTheta_outer,1), numel(Tfull));   % [theta x time]

for k = 1:numel(Tdata)
    idx = find(Tfull == Tdata(k));
    if ~isempty(idx)
        vTheta_inner_full(:,idx) = vAnomalyTheta_inner(:,k);
        vTheta_outer_full(:,idx) = vAnomalyTheta_outer(:,k);
    end
end

% ------------------------------------------------
% Common color limit
% ------------------------------------------------
cmaxRegion = max(abs([vTheta_inner_full(:); vTheta_outer_full(:)]));

f = newFigure([],figPos);
tlo = tiledlayout(1,2,'TileSpacing','tight','Padding','tight');

% ------------------------------------------------
% INNER REGION: rho <= 1 km
% ------------------------------------------------
ax1 = nexttile;
hold(ax1,'on'); box(ax1,'on'); grid(ax1,'on')

pp1 = pcolor(ax1, ...
    repmat(thetaDegNorth(1:end-1)',1,length(T_hov)), ...
    repmat(T_hov,length(thetaDegNorth(1:end-1)),1), ...
    [vTheta_inner_full NaN(size(vTheta_inner_full,1),1)]);
pp1.EdgeColor = 'none';

colormap(ax1,cmapHov)
caxis(ax1,[-cmaxRegion/1.5 cmaxRegion/1.5])

xlabel(ax1,'\theta (\circ)')
ylabel(ax1,'Year')
title(ax1,'Changes in Volume Anomaly [\theta] | \rho \leq 1 km')

xlim(ax1,[min(thetaDegNorth) max(thetaDegNorth)])
ylim(ax1,[T_hov(1) T_hov(end)])

set(ax1,'XTick',min(thetaDegNorth):30:max(thetaDegNorth))
set(ax1,'YTick',Tfull)
set(ax1,'Layer','top','GridColor','k','GridAlpha',0.35)
axis(ax1,'square')

% ------------------------------------------------
% OUTER REGION: rho > 1 km
% ------------------------------------------------
ax2 = nexttile;
hold(ax2,'on'); box(ax2,'on'); grid(ax2,'on')

pp2 = pcolor(ax2, ...
    repmat(thetaDegNorth(1:end-1)',1,length(T_hov)), ...
    repmat(T_hov,length(thetaDegNorth(1:end-1)),1), ...
    [vTheta_outer_full NaN(size(vTheta_outer_full,1),1)]);
pp2.EdgeColor = 'none';

colormap(ax2,cmapHov)
caxis(ax2,[-cmaxRegion/1.5 cmaxRegion/1.5])

xlabel(ax2,'\theta (\circ)')
ylabel(ax2,'Year')
title(ax2,'Changes in Volume Anomaly [\theta] | \rho > 1 km')

xlim(ax2,[min(thetaDegNorth) max(thetaDegNorth)])
ylim(ax2,[T_hov(1) T_hov(end)])

set(ax2,'XTick',min(thetaDegNorth):30:max(thetaDegNorth))
set(ax2,'YTick',Tfull)
set(ax2,'Layer','top','GridColor','k','GridAlpha',0.35)
axis(ax2,'square')

cb = colorbar(ax2,'Location','eastoutside');
ylabel(cb,'V_{anomaly} (m^3)')

applyFigureStyle(f,fontName,fontSize)

pngName = [originName{oo} '_YEAR_HovmollerTheta_InnerOuter'];
savefig(fullfile(outPath,[pngName '.fig']))
exportgraphics(f,fullfile(outPath,[pngName '.png']),'Resolution',figRes)

%% ================================================================
%% FIGURE 4C - Hovmoller in rho separated by rho regions
%% Inner region: rho <= 0.8 km
%% Outer region: rho > 0.8 km
%% ================================================================

% ------------------------------------------------
% Define rho regions
% ------------------------------------------------
rhoKm = rho(1:end-1)./1000;

idxInner = rhoKm <= 0.8;
idxOuter = rhoKm > 0.8;

% ------------------------------------------------
% Build continuous yearly matrix with NaN for missing years
% vAnomalyRho has dimensions [rho x time]
% ------------------------------------------------
Tdata = T(:)';
Tfull = Tdata(1):Tdata(end);
T_hov = [Tfull Tfull(end)+1];

if isempty(yearsToShow)
    yearsToShow = Tdata;
end

vRho_full = NaN(size(vAnomalyRho,1), numel(Tfull));   % [rho x time]

for k = 1:numel(Tdata)
    idx = find(Tfull == Tdata(k));
    if ~isempty(idx)
        vRho_full(:,idx) = vAnomalyRho(:,k);
    end
end

% ------------------------------------------------
% Mask inner and outer regions
% ------------------------------------------------
vRho_inner_full = vRho_full;
vRho_outer_full = vRho_full;

vRho_inner_full(~idxInner,:) = NaN;
vRho_outer_full(~idxOuter,:) = NaN;

% ------------------------------------------------
% Common color limit
% ------------------------------------------------
cmaxRegion = max(abs([vRho_inner_full(:); vRho_outer_full(:)]));

f = newFigure([],figPos);
tlo = tiledlayout(1,2,'TileSpacing','tight','Padding','tight');

% ------------------------------------------------
% INNER REGION: rho <= 0.8 km
% ------------------------------------------------
ax1 = nexttile;
hold(ax1,'on'); box(ax1,'on'); grid(ax1,'on')

pp1 = pcolor(ax1, ...
    repmat(T_hov,length(rhoKm),1), ...
    repmat(rhoKm,length(T_hov),1)', ...
    [vRho_inner_full NaN(size(vRho_inner_full,1),1)]);
pp1.EdgeColor = 'none';

colormap(ax1,cmapHov)
caxis(ax1,[-cmaxRegion/1.5 cmaxRegion/1.5])

xlabel(ax1,'Year')
ylabel(ax1,'\rho (km)')
title(ax1,'Changes in Volume Anomaly [\rho] | \rho \leq 0.8 km')

xlim(ax1,[T_hov(1) T_hov(end)])
ylim(ax1,[min(rhoKm) max(rhoKm)])

set(ax1,'XTick',Tfull)
set(ax1,'Layer','top','GridColor','k','GridAlpha',0.35)
axis(ax1,'square')

% ------------------------------------------------
% OUTER REGION: rho > 0.8 km
% ------------------------------------------------
ax2 = nexttile;
hold(ax2,'on'); box(ax2,'on'); grid(ax2,'on')

pp2 = pcolor(ax2, ...
    repmat(T_hov,length(rhoKm),1), ...
    repmat(rhoKm,length(T_hov),1)', ...
    [vRho_outer_full NaN(size(vRho_outer_full,1),1)]);
pp2.EdgeColor = 'none';

colormap(ax2,cmapHov)
caxis(ax2,[-cmaxRegion/1.5 cmaxRegion/1.5])

xlabel(ax2,'Year')
ylabel(ax2,'\rho (km)')
title(ax2,'Changes in Volume Anomaly [\rho] | \rho > 0.8 km')

xlim(ax2,[T_hov(1) T_hov(end)])
ylim(ax2,[min(rhoKm) max(rhoKm)])

set(ax2,'XTick',Tfull)
set(ax2,'Layer','top','GridColor','k','GridAlpha',0.35)
axis(ax2,'square')

cb = colorbar(ax2,'Location','eastoutside');
ylabel(cb,'V_{anomaly} (m^3)')

applyFigureStyle(f,fontName,fontSize)

pngName = [originName{oo} '_YEAR_HovmollerRho_InnerOuter'];
savefig(fullfile(outPath,[pngName '.fig']))
exportgraphics(f,fullfile(outPath,[pngName '.png']),'Resolution',figRes)


%% ================================================================
%% FIGURE 4C - Single Hovmoller in theta with two bands per year
%% First band:  rho <= 1 km
%% Second band: rho > 1 km
%% ================================================================

% ------------------------------------------------
% Define rho regions
% ------------------------------------------------
rhoKm = rho(1:end-1)./1000;

idxInner = rhoKm <= 0.8;   % first band
idxOuter = rhoKm > 0.8;    % second band

% ------------------------------------------------
% Collapse in theta for each rho region
% output: [theta x time]
% ------------------------------------------------
vAnomalyTheta_inner = squeeze(sum(vAnomalyPolar(idxInner,:,:),1,'omitnan'));
vAnomalyTheta_outer = squeeze(sum(vAnomalyPolar(idxOuter,:,:),1,'omitnan'));

vAnomalyTheta_inner = reshape(vAnomalyTheta_inner,length(thetaDegNorth)-1,length(T));
vAnomalyTheta_outer = reshape(vAnomalyTheta_outer,length(thetaDegNorth)-1,length(T));

% ------------------------------------------------
% Build full yearly arrays with NaN in missing years
% ------------------------------------------------
Tdata = T(:)';
Tfull = Tdata(1):Tdata(end);
nYears = numel(Tfull);

vTheta_inner_full = NaN(length(thetaDegNorth)-1, nYears);
vTheta_outer_full = NaN(length(thetaDegNorth)-1, nYears);

for k = 1:numel(Tdata)
    idx = find(Tfull == Tdata(k));
    if ~isempty(idx)
        vTheta_inner_full(:,idx) = vAnomalyTheta_inner(:,k);
        vTheta_outer_full(:,idx) = vAnomalyTheta_outer(:,k);
    end
end

% ------------------------------------------------
% Interleave rows by year:
% row 1 = year1 inner
% row 2 = year1 outer
% row 3 = year2 inner
% row 4 = year2 outer
% ...
% Final size: [2*nYears x nTheta]
% ------------------------------------------------
nTheta = length(thetaDegNorth)-1;
vTheta_banded = NaN(2*nYears, nTheta);

for yy = 1:nYears
    vTheta_banded(2*yy-1,:) = vTheta_inner_full(:,yy)'; % first: 0-1 km
    vTheta_banded(2*yy  ,:) = vTheta_outer_full(:,yy)'; % second: >1 km
end

% ------------------------------------------------
% Prepare pcolor inputs
% pcolor works better with one extra column in C
% ------------------------------------------------
thetaPlot = thetaDegNorth(1:end-1);        % 1 x nTheta
yPlot     = 1:(2*nYears);                  % 1 x (2*nYears)

Cplot = [vTheta_banded NaN(2*nYears,1)];   % [2*nYears x (nTheta+1)]

Xplot = repmat([thetaPlot thetaPlot(end)], 2*nYears, 1);   % same size as Cplot
Yplot = repmat(yPlot', 1, nTheta+1);                       % same size as Cplot

% ------------------------------------------------
% Color limit
% ------------------------------------------------
cmaxBand = max(abs(vTheta_banded(:)));

% ------------------------------------------------
% Plot
% ------------------------------------------------
f = newFigure([],figPos);
ax = axes('Parent',f);
hold(ax,'on'); box(ax,'on'); grid(ax,'on')

pp = pcolor(ax, Xplot, Yplot, Cplot);
pp.EdgeColor = 'none';

colormap(ax,cmapHov)
caxis(ax,[-cmaxBand/1.5 cmaxBand/1.5])

xlabel(ax,'\theta (\circ)')
ylabel(ax,'Year')
title(ax,'Changes in Volume Anomaly [\theta] | 0-1 km then >1 km')

xlim(ax,[min(thetaPlot) max(thetaPlot)])
ylim(ax,[0.5 2*nYears+0.5])

set(ax,'XTick',min(thetaPlot):30:max(thetaPlot))

% year label centered between each pair of bands
yearTickPos = 1.5:2:(2*nYears-0.5);
set(ax,'YTick',yearTickPos)
set(ax,'YTickLabel',string(Tfull))

set(ax,'Layer','top','GridColor','k','GridAlpha',0.35)

cb = colorbar(ax,'Location','eastoutside');
ylabel(cb,'V_{anomaly} (m^3)')

applyFigureStyle(f,fontName,fontSize)

pngName = [originName{oo} '_YEAR_HovmollerTheta_BandedInnerOuter'];
savefig(fullfile(outPath,[pngName '.fig']))
exportgraphics(f,fullfile(outPath,[pngName '.png']),'Resolution',figRes)

%% ================================================================
%% FIGURE 4D - Single Hovmoller in rho with two bands per year
%% First band:  theta sector 1
%% Second band: theta sector 2
%% ================================================================

% ------------------------------------------------
% Define theta sectors (edit these limits as needed)
% Example:
% sector 1 = shoal-dominated sector
% sector 2 = remaining sector
% ------------------------------------------------
thetaVals = thetaDegNorth(1:end-1);

idxSector1 = thetaVals >= -170 & thetaVals <= -135;   % first band
idxSector2 = thetaVals <  -170 | thetaVals >  -135;   % second band

% ------------------------------------------------
% Collapse in rho for each theta sector
% output: [rho x time]
% ------------------------------------------------
vAnomalyRho_sector1 = squeeze(sum(vAnomalyPolar(:,idxSector1,:),2,'omitnan'));
vAnomalyRho_sector2 = squeeze(sum(vAnomalyPolar(:,idxSector2,:),2,'omitnan'));

vAnomalyRho_sector1 = reshape(vAnomalyRho_sector1,length(rho)-1,length(T));
vAnomalyRho_sector2 = reshape(vAnomalyRho_sector2,length(rho)-1,length(T));

% ------------------------------------------------
% Build full yearly arrays with NaN in missing years
% ------------------------------------------------
Tdata = T(:)';
Tfull = Tdata(1):Tdata(end);
nYears = numel(Tfull);

vRho_sector1_full = NaN(length(rho)-1, nYears);
vRho_sector2_full = NaN(length(rho)-1, nYears);

for k = 1:numel(Tdata)
    idx = find(Tfull == Tdata(k));
    if ~isempty(idx)
        vRho_sector1_full(:,idx) = vAnomalyRho_sector1(:,k);
        vRho_sector2_full(:,idx) = vAnomalyRho_sector2(:,k);
    end
end

% ------------------------------------------------
% Interleave rows by year:
% row 1 = year1 sector1
% row 2 = year1 sector2
% row 3 = year2 sector1
% row 4 = year2 sector2
% ...
% Final size: [2*nYears x nRho]
% ------------------------------------------------
nRho = length(rho)-1;
vRho_banded = NaN(2*nYears, nRho);

for yy = 1:nYears
    vRho_banded(2*yy-1,:) = vRho_sector1_full(:,yy)'; % first band
    vRho_banded(2*yy  ,:) = vRho_sector2_full(:,yy)'; % second band
end

% ------------------------------------------------
% Prepare pcolor inputs
% pcolor works better with one extra column in C
% ------------------------------------------------
rhoPlot = rho(1:end-1)./1000;             % 1 x nRho
yPlot   = 1:(2*nYears);                   % 1 x (2*nYears)

Cplot = [vRho_banded NaN(2*nYears,1)];    % [2*nYears x (nRho+1)]

Xplot = repmat([rhoPlot rhoPlot(end)], 2*nYears, 1);   % same size as Cplot
Yplot = repmat(yPlot', 1, nRho+1);                      % same size as Cplot

% ------------------------------------------------
% Color limit
% ------------------------------------------------
cmaxBand = max(abs(vRho_banded(:)));

% ------------------------------------------------
% Plot
% ------------------------------------------------
f = newFigure([],figPos);
ax = axes('Parent',f);
hold(ax,'on'); box(ax,'on'); grid(ax,'on')

pp = pcolor(ax, Xplot, Yplot, Cplot);
pp.EdgeColor = 'none';

colormap(ax,cmapHov)
caxis(ax,[-cmaxBand/1.5 cmaxBand/1.5])

xlabel(ax,'\rho (km)')
ylabel(ax,'Year')
title(ax,'Changes in Volume Anomaly [\rho] | Sector 1 then Sector 2')

xlim(ax,[min(rhoPlot) max(rhoPlot)])
ylim(ax,[0.5 2*nYears+0.5])

% year label centered between each pair of bands
yearTickPos = 1.5:2:(2*nYears-0.5);
set(ax,'YTick',yearTickPos)
set(ax,'YTickLabel',string(Tfull))

set(ax,'Layer','top','GridColor','k','GridAlpha',0.35)

cb = colorbar(ax,'Location','eastoutside');
ylabel(cb,'V_{anomaly} (m^3)')

applyFigureStyle(f,fontName,fontSize)

pngName = [originName{oo} '_YEAR_HovmollerRho_BandedSector1Sector2'];
savefig(fullfile(outPath,[pngName '.fig']))
exportgraphics(f,fullfile(outPath,[pngName '.png']),'Resolution',figRes)

%% ================================================================
%% FIGURE 4D - rho-theta vAnomaly separated by rho regions
%% One figure for rho <= 1 km
%% One figure for rho > 1 km
%% ================================================================
rhoKm = rho(1:end-1)./1000;

idxInner = rhoKm <= 0.8;
idxOuter = rhoKm > 0.8;

maxC = max(vAnomalyPolar,[],'all');
climSym = [-maxC/1.5 maxC/1.5];

thetaL = [min(thetaDegNorth) max(thetaDegNorth)];
rhoL   = [min(rho./1000) max(rho./1000)];

for tt = 1:length(T)

    % ------------------------------------------------
    % INNER REGION: rho <= 1 km
    % ------------------------------------------------
    vPlotInner = vAnomalyPolar(:,:,tt);
    vPlotInner(~idxInner,:) = NaN;

    f1 = figure;
    set(f1,'Units','centimeters')
    set(f1,'Position',[2 2 10.4 10.4])
    set(f1,'Color','w')

    ax1 = axes('Parent',f1);
    hold(ax1,'on'); box(ax1,'on'); grid(ax1,'on')

    pc1 = pcolor(ax1,thetaGrid,rhoGrid,vPlotInner);
    pc1.FaceColor = 'flat';
    pc1.EdgeColor = 'none';

    colormap(ax1,cmapHov)
    caxis(ax1,climSym)

    xlabel(ax1,'\theta (\circ)')
    ylabel(ax1,'\rho (km)')
    title(ax1,['\rho-\theta Space [' num2str(T(tt)) '] | \rho \leq 1 km'])

    xlim(ax1,thetaL)
    ylim(ax1,rhoL)
    axis(ax1,'square')
    set(ax1,'XTick',min(thetaDegNorth):30:max(thetaDegNorth))
    set(ax1,'Color',[0.98 0.98 0.98],'Layer','top','GridAlpha',0.2,'GridColor','k')

    yline(ax1,1.0,'k--','LineWidth',1.0)

    cb1 = colorbar(ax1,'southoutside');
    cb1.Label.String = 'V anomaly (m^3)';
    cb1.Ruler.Exponent = 4;

    set(findall(f1,'Type','axes'),'FontName','Arial','FontSize',11)
    set(findall(f1,'Type','text'),'FontName','Arial','FontSize',11)
    set(findall(f1,'Type','ColorBar'),'FontName','Arial','FontSize',11)

    pngName = [originName{oo} '_RhoTheta_vAnomaly_Inner_' num2str(T(tt))];
    savefig(fullfile(outPath,[pngName '.fig']))
    exportgraphics(f1,fullfile(outPath,[pngName '.png']),'Resolution',300)
    close(f1)

    % ------------------------------------------------
    % OUTER REGION: rho > 1 km
    % ------------------------------------------------
    vPlotOuter = vAnomalyPolar(:,:,tt);
    vPlotOuter(~idxOuter,:) = NaN;

    f2 = figure;
    set(f2,'Units','centimeters')
    set(f2,'Position',[2 2 10.4 10.4])
    set(f2,'Color','w')

    ax2 = axes('Parent',f2);
    hold(ax2,'on'); box(ax2,'on'); grid(ax2,'on')

    pc2 = pcolor(ax2,thetaGrid,rhoGrid,vPlotOuter);
    pc2.FaceColor = 'flat';
    pc2.EdgeColor = 'none';

    colormap(ax2,cmapHov)
    caxis(ax2,climSym)

    xlabel(ax2,'\theta (\circ)')
    ylabel(ax2,'\rho (km)')
    title(ax2,['\rho-\theta Space [' num2str(T(tt)) '] | \rho > 1 km'])

    xlim(ax2,thetaL)
    ylim(ax2,rhoL)
    axis(ax2,'square')
    set(ax2,'XTick',min(thetaDegNorth):30:max(thetaDegNorth))
    set(ax2,'Color',[0.98 0.98 0.98],'Layer','top','GridAlpha',0.2,'GridColor','k')

    yline(ax2,1.0,'k--','LineWidth',1.0)

    cb2 = colorbar(ax2,'southoutside');
    cb2.Label.String = 'V anomaly (m^3)';
    cb2.Ruler.Exponent = 4;

    set(findall(f2,'Type','axes'),'FontName','Arial','FontSize',11)
    set(findall(f2,'Type','text'),'FontName','Arial','FontSize',11)
    set(findall(f2,'Type','ColorBar'),'FontName','Arial','FontSize',11)

    pngName = [originName{oo} '_RhoTheta_vAnomaly_Outer_' num2str(T(tt))];
    savefig(fullfile(outPath,[pngName '.fig']))
    exportgraphics(f2,fullfile(outPath,[pngName '.png']),'Resolution',300)
    close(f2)
end

%% ================================================================
%% FIGURE 5 - Clean barcode for each year (theta)
%% ================================================================
for tt = 1:length(T)

    f = figure;
    set(f,'Units','centimeters')
    set(f,'Position',[2 2 21 2.2])   % largura e altura da faixa
    set(f,'Color','w')

    ax = axes(f);
    set(ax,'Position',[0.02 0.15 0.96 0.7])   % ocupa quase toda a figura
    hold(ax,'on')

    % repete a série em poucas linhas para formar a faixa
    barcodeTheta = repmat(vAnomalyTheta(:,tt)', 2, 1);

    pcolor(ax,thetaDegNorth(1:end-1),1:2,barcodeTheta)
    shading(ax,'flat')

    colormap(ax,cmapHov)
%     clim95 = prctile(abs(vAnomalyPolar(:)),95);
%     caxis(ax,[-clim95 clim95])

    caxis(ax,[-cmax/1.5 cmax/1.5])

    % remover tudo
    axis(ax,'off')

    % salvar
    pngName = [originName{oo} '_BarcodeTheta_' num2str(T(tt))];
    savefig(fullfile(outPath,[pngName '.fig']))
    exportgraphics(f,fullfile(outPath,[pngName '.png']),'Resolution',300)

    close(f)
end
%% ================================================================
%% FIGURE 6 - Clean barcode for each year (rho)
%% ================================================================
for tt = 1:length(T)

    f = figure;
    set(f,'Units','centimeters')
    set(f,'Position',[2 2 2.2 21])   % faixa vertical
    set(f,'Color','w')

    ax = axes(f);
    set(ax,'Position',[0.15 0.02 0.7 0.96])   % ocupa quase toda a figura
    hold(ax,'on')

    % repete a série em poucas colunas para formar a faixa
    barcodeRho = repmat(vAnomalyRho(:,tt), 1, 2);

    pcolor(ax,1:2,rho(1:end-1)./1000,barcodeRho)
    shading(ax,'flat')

    colormap(ax,cmapHov)
    caxis(ax,[-cmax/1.5 cmax/1.5])

    % remover tudo
    axis(ax,'off')

    % salvar
    pngName = [originName{oo} '_BarcodeRho_' num2str(T(tt))];
    savefig(fullfile(outPath,[pngName '.fig']))
    exportgraphics(f,fullfile(outPath,[pngName '.png']),'Resolution',300)

    close(f)
end

%% ================================================================
%% FIGURE SUMMARY SEPARATED - save each panel as an independent figure
%% ================================================================
maxC = max(vAnomalyPolar,[],'all');
climSym = [-maxC/1.5 maxC/1.5];

for tt = 1:length(T)

    % -------- common limits --------
    xL     = [min(xPolar(:)) max(xPolar(:))];
    yL     = [min(yPolar(:)) max(yPolar(:))];
    thetaL = [min(thetaDegNorth) max(thetaDegNorth)];
    rhoL   = [min(rho./1000) max(rho./1000)];

    xp = xPolar(1:end-1,1:end-1);
    yp = yPolar(1:end-1,1:end-1);

    %% ------------------------------------------------
    %% FIGURE A - XY bathy
    %% ------------------------------------------------
    f1 = figure;
    set(f1,'Units','centimeters')
    set(f1,'Position',[2 2 10 10])
    set(f1,'Color','w')

    ax1 = axes('Parent',f1);
    hold(ax1,'on'); box(ax1,'on'); grid(ax1,'on')

    contourf(ax1,xp,yp,zPolar(:,:,tt),edges,'LineStyle','none');
    colormap(ax1,cmapBath)
    caxis(ax1,[edges(1) edges(end)])

    contour(ax1,xp,yp,zPolar(:,:,tt),[0 0],'k-','LineWidth',1.0)
    contour(ax1,xp,yp,zPolar(:,:,tt),[-5 -5],'k--','LineWidth',1.0)

    for ii = 1:length(thetaPlotGrid)
        plot(ax1,xPolarPlotGrid(:,ii),yPolarPlotGrid(:,ii),'-', ...
            'Color',[0.25 0.25 0.25],'LineWidth',0.6)
    end
    for ii = 1:length(rhoPlotGrid)
        plot(ax1,xPolarPlotGridFine(ii,:),yPolarPlotGridFine(ii,:),'-', ...
            'Color',[0.25 0.25 0.25],'LineWidth',0.6)
    end

    xlabel(ax1,'X (km)')
    ylabel(ax1,'Y (km)')
    title(ax1,['X-Y Space [' num2str(T(tt)) ']'])

    xlim(ax1,xL)
    ylim(ax1,yL)
    axis(ax1,'square')
    set(ax1,'Color',[0.98 0.98 0.98],'Layer','top','GridAlpha',0.2,'GridColor','k')

    cb1 = colorbar(ax1,'southoutside');
    cb1.Ticks = tickCenters;
    cb1.TickLabels = tickLabels;
    cb1.Label.String = 'Elevation (m)';

    set(findall(f1,'Type','axes'),'FontName','Arial','FontSize',11)
    set(findall(f1,'Type','text'),'FontName','Arial','FontSize',11)
    set(findall(f1,'Type','ColorBar'),'FontName','Arial','FontSize',11)

    pngName = [originName{oo} '_XY_Bathy_' num2str(T(tt))];
    savefig(fullfile(outPath,[pngName '.fig']))
    exportgraphics(f1,fullfile(outPath,[pngName '.png']),'Resolution',300)
    close(f1)

    %% ------------------------------------------------
    %% FIGURE B - rho-theta bathy
    %% ------------------------------------------------
    f2 = figure;
    set(f2,'Units','centimeters')
    set(f2,'Position',[2 2 10.4 10.4])
    set(f2,'Color','w')

    ax2 = axes('Parent',f2);
    hold(ax2,'on'); box(ax2,'on'); grid(ax2,'on')

    contourf(ax2,thetaGrid,rhoGrid,zPolar(:,:,tt),edges,'LineStyle','none');
    colormap(ax2,cmapBath)
    caxis(ax2,[edges(1) edges(end)])

    contour(ax2,thetaGrid,rhoGrid,zPolar(:,:,tt),[0 0],'k-','LineWidth',1.0)
    contour(ax2,thetaGrid,rhoGrid,zPolar(:,:,tt),[-5 -5],'k--','LineWidth',1.0)

    xlabel(ax2,'\theta (\circ)')
    ylabel(ax2,'\rho (km)')
    title(ax2,['\rho-\theta Space [' num2str(T(tt)) ']'])

    xlim(ax2,thetaL)
    ylim(ax2,rhoL)
    axis(ax2,'square')
    set(ax2,'XTick',min(thetaDegNorth):30:max(thetaDegNorth))
    set(ax2,'Color',[0.98 0.98 0.98],'Layer','top','GridAlpha',0.2,'GridColor','k')

    cb2 = colorbar(ax2,'southoutside');
    cb2.Ticks = tickCenters;
    cb2.TickLabels = tickLabels;
    cb2.Label.String = 'Elevation (m)';

    set(findall(f2,'Type','axes'),'FontName','Arial','FontSize',11)
    set(findall(f2,'Type','text'),'FontName','Arial','FontSize',11)
    set(findall(f2,'Type','ColorBar'),'FontName','Arial','FontSize',11)

    pngName = [originName{oo} '_RhoTheta_Bathy_' num2str(T(tt))];
    savefig(fullfile(outPath,[pngName '.fig']))
    exportgraphics(f2,fullfile(outPath,[pngName '.png']),'Resolution',300)
    close(f2)

    %% ------------------------------------------------
    %% FIGURE C - rho-theta vAnomaly
    %% ------------------------------------------------
    f3 = figure;
    set(f3,'Units','centimeters')
    set(f3,'Position',[2 2 10.4 10.4])
    set(f3,'Color','w')

    ax3 = axes('Parent',f3);
    hold(ax3,'on'); box(ax3,'on'); grid(ax3,'on')

    pc = pcolor(ax3,thetaGrid,rhoGrid,vAnomalyPolar(:,:,tt));
    pc.FaceColor = 'flat';
    pc.EdgeColor = 'none';

    colormap(ax3,cmapHov)
    caxis(ax3,climSym)

    xlabel(ax3,'\theta (\circ)')
    ylabel(ax3,'\rho (km)')
    title(ax3,['\rho-\theta Space [' num2str(T(tt)) ']'])

    xlim(ax3,thetaL)
    ylim(ax3,rhoL)
    axis(ax3,'square')
    set(ax3,'XTick',min(thetaDegNorth):30:max(thetaDegNorth))
    set(ax3,'Color',[0.98 0.98 0.98],'Layer','top','GridAlpha',0.2,'GridColor','k')

    cb3 = colorbar(ax3,'southoutside');
    cb3.Label.String = 'V anomaly (m^3)';
    cb3.Ruler.Exponent = 4;

    set(findall(f3,'Type','axes'),'FontName','Arial','FontSize',11)
    set(findall(f3,'Type','text'),'FontName','Arial','FontSize',11)
    set(findall(f3,'Type','ColorBar'),'FontName','Arial','FontSize',11)

    pngName = [originName{oo} '_RhoTheta_vAnomaly_' num2str(T(tt))];
    savefig(fullfile(outPath,[pngName '.fig']))
    exportgraphics(f3,fullfile(outPath,[pngName '.png']),'Resolution',300)
    close(f3)

    

        %% ------------------------------------------------
    %% FIGURE E - XY weighted mean surface
    %% ------------------------------------------------
    f4 = figure;
    set(f4,'Units','centimeters')
    set(f4,'Position',[2 2 10 10])
    set(f4,'Color','w')

    ax4 = axes('Parent',f4);
    hold(ax4,'on'); box(ax4,'on'); grid(ax4,'on')

    contourf(ax4,xp,yp,meanSurfPolar,edges,'LineStyle','none');
    colormap(ax4,cmapBath)
    caxis(ax4,[edges(1) edges(end)])

    contour(ax4,xp,yp,meanSurfPolar,[0 0],'k-','LineWidth',1.0)
    contour(ax4,xp,yp,meanSurfPolar,[-5 -5],'k--','LineWidth',1.0)

    for ii = 1:length(thetaPlotGrid)
        plot(ax4,xPolarPlotGrid(:,ii),yPolarPlotGrid(:,ii),'-', ...
            'Color',[0.25 0.25 0.25],'LineWidth',0.6)
    end
    for ii = 1:length(rhoPlotGrid)
        plot(ax4,xPolarPlotGridFine(ii,:),yPolarPlotGridFine(ii,:),'-', ...
            'Color',[0.25 0.25 0.25],'LineWidth',0.6)
    end

    xlabel(ax4,'X (km)')
    ylabel(ax4,'Y (km)')
    title(ax4,'X-Y Mean Surface')

    xlim(ax4,xL)
    ylim(ax4,yL)
    axis(ax4,'square')
    set(ax4,'Color',[0.98 0.98 0.98],'Layer','top','GridAlpha',0.2,'GridColor','k')

    cb4 = colorbar(ax4,'southoutside');
    cb4.Ticks = tickCenters;
    cb4.TickLabels = tickLabels;
    cb4.Label.String = 'Elevation (m)';

    set(findall(f4,'Type','axes'),'FontName','Arial','FontSize',11)
    set(findall(f4,'Type','text'),'FontName','Arial','FontSize',11)
    set(findall(f4,'Type','ColorBar'),'FontName','Arial','FontSize',11)

    pngName = [originName{oo} '_XY_MeanSurface_' num2str(T(tt))];
    savefig(fullfile(outPath,[pngName '.fig']))
    exportgraphics(f4,fullfile(outPath,[pngName '.png']),'Resolution',300)
    close(f4)
    
        %% ------------------------------------------------
    %% FIGURE 5 - XY vAnomaly
    %% ------------------------------------------------
    f5 = figure;
    set(f5,'Units','centimeters')
    set(f5,'Position',[2 2 10 10])
    set(f5,'Color','w')

    ax5 = axes('Parent',f5);
    hold(ax5,'on'); box(ax5,'on'); grid(ax5,'on')

    pc = pcolor(ax5,xp,yp,vAnomalyPolar(:,:,tt));
    pc.FaceColor = 'flat';
    pc.EdgeColor = 'none';

    % use cmapHovFlip se quiser positivo=azul e negativo=vermelho
    colormap(ax5,cmapHov)
    caxis(ax5,climSym)

    for ii = 1:length(thetaPlotGrid)
        plot(ax5,xPolarPlotGrid(:,ii),yPolarPlotGrid(:,ii),'-', ...
            'Color',[0.25 0.25 0.25],'LineWidth',0.6)
    end
    for ii = 1:length(rhoPlotGrid)
        plot(ax5,xPolarPlotGridFine(ii,:),yPolarPlotGridFine(ii,:),'-', ...
            'Color',[0.25 0.25 0.25],'LineWidth',0.6)
    end

    xlabel(ax5,'X (km)')
    ylabel(ax5,'Y (km)')
    title(ax5,['X-Y Space [' num2str(T(tt)) ']'])

    xlim(ax5,xL)
    ylim(ax5,yL)
    axis(ax5,'square')
    set(ax5,'Color',[0.98 0.98 0.98], ...
        'Layer','top','GridAlpha',0.2,'GridColor','k')

    cb5 = colorbar(ax5,'southoutside');
    cb5.Label.String = 'V anomaly (m^3)';
    cb5.Ruler.Exponent = 4;

    set(findall(f5,'Type','axes'),'FontName','Arial','FontSize',11)
    set(findall(f5,'Type','text'),'FontName','Arial','FontSize',11)
    set(findall(f5,'Type','ColorBar'),'FontName','Arial','FontSize',11)

    pngName = [originName{oo} '_XY_vAnomaly_' num2str(T(tt))];
    savefig(fullfile(outPath,[pngName '.fig']))
    exportgraphics(f5,fullfile(outPath,[pngName '.png']),'Resolution',300)
    close(f5)
end
end
end

%% ========================== HELPERS ==================================
function f = newFigure(figNum,figPos)
if isempty(figNum)
    f = figure;
else
    f = figure(figNum);
end
clf(f)
set(f,'WindowState','normal')
set(f,'Units','centimeters')
set(f,'Position',figPos)
end

function applyFigureStyle(f,fontName,fontSize)
set(findall(f,'Type','axes'),'FontName',fontName,'FontSize',fontSize)
set(findall(f,'Type','text'),'FontName',fontName,'FontSize',fontSize)
set(findall(f,'Type','ColorBar'),'FontName',fontName,'FontSize',fontSize)
end


%% Kathy Code_Moficado por Matheus 03/04/2026
% function polarAnalysis(datafile,origin,outPathRoot,makeMiddlePlots,dTheta,dRho,rhoLimOuter,thetaLims,testMultiOrigins,interest)
% %polarAnalysis.m 
% %Conformal mapping to polar coordinates of inlet bathymetric surface data
% %to quantify and track inlet geomorphic feature migration and rotation rates. 
% 
% %INPUTS:
%     %datafile: .mat file (ex: 'bathy_data.mat') containing a structure called 'grd' with fields: 
%         %year: year for each data surface, in 'double' format. ex: [2010, 2015, 2020,
%             %2023] %%CANNOT HAVE MULTIPLES OF SAME YEAR, TO DO fix this
%         %x: gridded x data in meters, in projected coordinates, ie
%             %Eastings, size m x n
%         %y: gridded y data in meters, in projected coordinates, ie
%             %Northings, size m x n
%         %dp; gridded z data, in meters, depth negative, size m x n x year,
%             %topo data should already be masked out
%     %origin: x,y coordinates (meters) of origin for polar grid, usually center of
%         %inlet thalweg, ex: [547014,3361720]
%     %outPath: filepath to folder where results will be output
%     %makeMiddlePlots: make or don't make plots of hAnomaly and vAnomaly for all
%         %surveys. 0 or 1.
%     %dRho: polar grid spacing in radial direction (m), ex: 40 
%     %rhoLimOuter: polar grid outer limits (meters from origin). ex: rhoLim = 1500, makes a
%         %grid from 0 (origin) to 1500 meters from origin, spaced by dRho.
%         %Choose limits which cover the ebb shoal features of
%         %interest. The inner 15% of the polar grid will be removed from analysis, as the
%         %grid cells become to small to have data points in every grid cell.
%     %thetaLims: polar grid sector limits, 180 degrees, in nautical
%         %convention (0 = North, 90 = East), ex: [0 180] for an east facing
%         %inlet
%     %dTheta: polar grid angle spacing (degrees), ex: 1 or -1
%         %TO DO: automatically correct for dTheta sign, for now, do positive
%         %dTheta if in the trigonometric positive half of the unit circle
%         %(0-180, or 90-270 in nautical convention) and do negative dTheta
%         %if in the trigonometric negative half of the unit circle (-180-0,
%         %or 270-90 in nautical convention)
%     %testMultiOrigins: 0 or 1. If 1, test 12 other origins around the chosen origin. 
%         %Compares ebb tidal delta rotation rate for different origins to chose ideal origin location. 
%         %From Pearson et al., 2022 Supplementary Matieral: Well-chosen
%         %origin should have a relatively higher delta rotation rate and
%         %motion around it will be more coherent. 
%         %Outputs delta rotation rates for each test origin and the origin site most
%         %recommended for use in analysis.
%         %WILL TAKE SIGNIFICANTLY LONGER
% 
% 
% %OUTPUTS:
%     %polarData.mat, output into outPath folder
%     %Plots of Z (m) in XY and polar space for each survey
%     %Hovmoller (2D timeseries) plots of change in volume anomaly for each
%         %survey summed along the theta and rho axes, indicating rotation
%         %and migration of features, respectively, through time. %FYSA the plots 
%         %have an extra year tacked onto the end of the
%         %timeseries so that the data calculated from the last survey will
%         %show. ex: If the last year of data in the timestack is 2020, 
%         %the Hovmoller plot time axes will end at 2021.
%     %Hovmoller (2D timeseries) plots with volume anomaly peak (shoal) tracking
%     %Text file of peak (shoal) rotation estimates (deg./yr)
%     %Plots of lagged cross correlation between consecutive volume anomalies
%         %with estimates of rotation rate (deg./yr) and peak cross correlation values
%         %(deg.)
%     %'*_Delta_Rotation_Rates.mat', containing the ETD rotation rate for
%         %each year (deg./yr) and the mean ETD rotation rate for the study period
%         %(deg./yr)
% 
%     %If makeMiddlePlots == 1
%         %Plots of hAnomaly (m) (Z deviation from Zmean) in XY and polar space for each
%             %survey
%         %Plots of vAnomaly (m^3) (volume of sediment above or below Zmean)
%             %in XY and polar space for each survey
%     
%     %If testMultiOrigins == 1
%         %Outputs all of the above for each of the 13 inlet origin test
%         %locations, in separate folders labelled by the origin name, ex:
%         %'O' for user input origin, 'N1' for 'North 1', 'E2' for 'East 2',
%         %etc
%         %Also outputs figures of the origin test locations and the rotation
%         %rates at each test location
%         %Also outputs Origin_Table.txt, with the Origin Name, origin
%         %location, and Mean ETD Rotation Rate (deg./yr) calculated from
%         %that origin
% 
% 
% %Note: Currently uses functions from the Climate Data Toolbox, including
% %imagescn.m (https://www.chadagreene.com/CDT/CDT_Contents.html).
% 
% %For original study, see: 
% %Pearson, S.G., Elias, E.P.L, van Prooijen, B.C., van der Vegt, H., van der Spek, 
% %A. & Wang, Z.B. (2022). A Novel Approach to Mapping Ebb-Tidal Delta Morphodynamics 
% %and Stratigraphy. Geomorphology. https://doi.org/10.1016/j.geomorph.2022.108185
% %and Supplementary Material
% 
% 
% %Original code from:
% %https://github.com/sgpearson17/bathy2strat
% 
% %Modified by Kaitlyn McPherran, USACE ERDC CHL, 2025
% %written in MATLAB v 2024a
% 
%     
% 
% %% Load Data
% 
% 
% 
% if exist(outPathRoot) == 7
%     addpath(outPathRoot); %add the output folder to path
% else %if folder does not exist
%     mkdir(outPathRoot) %make folder
%     addpath(outPathRoot) %and add folder to path
% end
% 
% load(datafile,'grd'); %load the mat file with the bathy time stacks
% 
% T = grd.year;
% if iscolumn(T); T = T'; end %make row wise
% Y = grd.y/1000;%km
% X = grd.x/1000;%km
% Z = grd.dp; %depth negative
% 
% 
% Z(Z> 0.5 | Z < -5) = NaN; % change elevation analyses
% 
% o0 = origin/1000; %convert to km
% 
% XLims = round([min(X,[],'all') max(X,[],'all')]); %km, find x limits 
% TLim = [1 length(T)]; 
% 
% %WIP 
% % %Convert nautical convention to trigonometric
% % thetaLims = thetaLims - 90;
% % for i = 1:length(thetaLims)   
% %     % Ensure the angle is within the range [-180, 180]
% %     if thetaLims(i) > 180
% %         thetaLims(i) = thetaLims(i) - 360;
% %     elseif thetaLims(i) < -180
% %         thetaLims(i) = thetaLims(i) + 360;
% %     end
% % end
% 
% if testMultiOrigins
%     %make each test origin 10% of total raster 'away' from user chosen origin in each
%     %direction
%     testDist = round(XLims(2) - o0(1)) * 0.1;  %Pearson et al tested ~10% in x and y distance
% 	
%     origins = [o0(1),o0(2);...
%     o0(1),o0(2)+testDist;...
%     o0(1)+testDist,o0(2);...
%     o0(1),o0(2)-testDist;...
%     o0(1)-testDist,o0(2);...
%     o0(1)+testDist,o0(2)+testDist;...
%     o0(1)+testDist,o0(2)-testDist;...
%     o0(1)-testDist,o0(2)+testDist;...
%     o0(1)-testDist,o0(2)-testDist;...
%     o0(1),o0(2)+testDist*2;...
%     o0(1)+testDist*2,o0(2);...
%     o0(1),o0(2)-testDist*2;...
%     o0(1)-testDist*2,o0(2);...
% 	]; %uncomment origins(6:13,:) to test 12 sites around the user's origin instead of 4 
%      originName = {'O';'N1';'E1';'S1';'W1';%;...
%         'NE1';'SE1';'NW1';'SW1';...
%         'N2';'E2';'S2';'W2'}; 
% 	
% 	figure 
%     clf
%     hold on
%     imagescn(X,Y,Z(:,:,1)); 
%     colormap('gray');
%     for i = 1:length(origins)
%         scatter(origins(i,1),origins(i,2),'filled');
%     end
%     colororder([
%     0.0, 0.45, 0.74;   % azul
%     0.85, 0.33, 0.10;  % laranja
%     0.93, 0.69, 0.13;  % amarelo
%     0.49, 0.18, 0.56;  % roxo
%     0.47, 0.67, 0.19;  % verde
%     0.30, 0.75, 0.93;  % ciano
%     0.64, 0.08, 0.18;  % vinho
%     0.25, 0.25, 0.25;  % cinza escuro
%     0.75, 0.75, 0.75;  % cinza claro
%     1.00, 0.60, 0.78;  % rosa
%     0.00, 0.50, 0.00;  % verde escuro
%     0.00, 0.00, 0.50   % azul escuro
% ]);
%     title('Origin Test Locations'); legend(originName); grid on; xlabel('X (km)'); ylabel('Y (km)');
%     %Export figure to png file
%     pngName = 'origin_test_locations';
%     savefig([outPathRoot filesep pngName '.fig']);
%     exportgraphics(gcf,[outPathRoot filesep pngName '.png'],"Resolution",300);
% else 
% 	origins = [o0(1),o0(2)];
% 	originName = {'O'};
% end
% 
% for oo = 1:size(origins,1)
%     if testMultiOrigins; disp("Now analyzing origin " + originName(oo)); end
% 
%     x0 = origins(oo,1); y0 = origins(oo,2); 
% 
%     outPath = [outPathRoot filesep originName{oo}];
%     if ~exist(outPath); mkdir(outPath); end
%     addpath(outPath)
% 
%     % perform weighted mean surface calculation
%     iT0=2; 
%     Tmidpt = NaN*ones(1,length(T)-1);
%     Tmidpt(1)=T(iT0);
%     tt = iT0:length(T)-1;
%     Tmidpt(2:end) = (T(tt)+T(tt+1))/2;
%     Tmidpt(length(T))=(T(end)+T(end)+1)/2; % assume current year has weight of one year, given annual surveys
%     Tweight = diff(Tmidpt); 
%     
%     %TO DO: fix any nans in any surface become nans in weightedMeanSurf
%     weightedMeanSurf=zeros(size(X)); 
%     for tt = iT0:length(T)
%         weightedMeanSurf = weightedMeanSurf + squeeze(Z(:,:,tt)).*Tweight(tt-1); 
%     end
%     weightedMeanSurf =  weightedMeanSurf./(T(end)-T(iT0));
%     
%     % calculate grid cell surface area
%     dx = unique(diff(X,1,2))*1000;
%     dy = unique(diff(Y,1,1))*1000;
%     areaGridCell = floor(abs(dx(1)*dy(1))); 
%     %areaGridCell = cellsize; 
% 
%     
%     % calculate the height above the mean surface (anomaly height)
%     hAnomaly = nan(size(Z));
%     hAnomaly(:,:,TLim(1):TLim(2)) = Z(:,:,TLim(1):TLim(2))-weightedMeanSurf;
%     hActive = nan(size(Z));
%     hActive(:,:,TLim(1):TLim(2)) = Z(:,:,TLim(1):TLim(2))-Z(:,:,TLim(1));
%     vActive = hActive .* areaGridCell; % in m^3
%     
%     % calculate volume above mean surface (Anomaly volume)
%     %vAnomaly = hAnomaly .* areaGridCell; % in m^3
%     vAnomaly = hAnomaly .* areaGridCell; % in m^3
%     %% set up grid in polar coordinates
%     
%     % distribution of angles
%     %TO DO maybe: allow 360 grid, track ETD and FTD shoals?
%     %TO DO maybe: rotate all data to the same direction and make all grids and axes uniform
%     
%     %Cutoff inner 15% of polar cells, too small to have much data in them
%     %Pearson did 1 km out of 7 ~ 14%
%     rhoLims = [round(rhoLimOuter*interest) rhoLimOuter]; 
%     theta = (thetaLims(1):dTheta:thetaLims(2)).*(2*pi/360);
%     rho = rhoLims(1):dRho:rhoLims(2);
%     [xPolar,yPolar] = pol2cart(repmat(theta,length(rho),1),repmat(rho,length(theta),1)');
%     
%     % convert polar coordinates to km
%     xPolar = xPolar./1000;
%     yPolar = yPolar./1000;
%     
%     %% Calculate bathymetric properties of polar grid
%     
%     % XYZ coordinates of bathymetry data (query points)
%     xq = reshape(X,[],1);
%     yq = reshape(Y,[],1);
%     
%     % initialize arrays 
%     zPolar = nan(length(rho)-1,length(theta)-1,length(T));
%     hAnomalyPolar = nan(length(rho)-1,length(theta)-1,length(T));
%     vAnomalyPolar = nan(length(rho)-1,length(theta)-1,length(T));
%     
%     for tt = 1:length(T)
%         disp(['tt=' num2str(T(tt))]);
%         % reshape list of Z/hAnomaly/vAnomaly for a given year
%         zq = reshape(Z(:,:,tt),[],1);
%         hq = reshape(hAnomaly(:,:,tt),[],1);
%         vq = reshape(vAnomaly(:,:,tt),[],1);
%         % Loop through and tabulate contents of each grid cell at each timestep
%         for jj = 1:length(theta)-1
%             if mod(theta(jj)./(2*pi/360),10)==0
%                 disp(['T=' num2str(T(tt)) ', theta=' num2str(theta(jj)./(2*pi/360))]);
%             end
%             for ii = 1:length(rho)-1
%                 % identify coordinates of a given polar grid cell (4 corners and repeat first point)
%                 xv = [xPolar(ii,jj) xPolar(ii+1,jj) xPolar(ii+1,jj+1) xPolar(ii,jj+1) xPolar(ii,jj)]+x0;
%                 yv = [yPolar(ii,jj) yPolar(ii+1,jj) yPolar(ii+1,jj+1) yPolar(ii,jj+1) yPolar(ii,jj)]+y0;
%                 % check which XY bathy coordinates lie inside a given grid cell
%                 [in] = inpolygon(xq,yq,xv,yv); 
%                 % multiply by grid cell area 
%                 zPolar(ii,jj,tt) = mean(zq(in),'omitnan'); 
%                 hAnomalyPolar(ii,jj,tt) = mean(hq(in),'omitnan'); 
%                 vAnomalyPolar(ii,jj,tt) = sum(vq(in)); %removed omitnan, otherwise all nans get filled with zeros in final vAnomaly surfs
%             end
%         end
%     end
%     
%     % calculate dynamic volumes
%     % collapse along Rho dimension
%     vAnomalyRho = squeeze(sum(vAnomalyPolar,2,'omitnan'));
%     % collapse along Theta dimension
%     vAnomalyTheta = squeeze(sum(vAnomalyPolar,1,'omitnan'));
%     
%     hActivePolar = nan(size(zPolar));
% %     hActivePolar(:,:,TLim(1):TLim(2)) = zPolar(:,:,TLim(1):TLim(2))-min(zPolar,[],3);
%     hActivePolar(:,:,TLim(1):TLim(2)) = zPolar(:,:,TLim(1):TLim(2)) - zPolar(:,:,TLim(1));
% 
%     vActivePolar = hActivePolar .* areaGridCell; %m^3
%       
%     % === NOVO: integrar também o volume "ativo" (mudança em relação ao 1º ano) ===
%     vActiveRho   = squeeze(sum(vActivePolar,2,'omitnan'));  % [rho x tempo]
%     vActiveTheta = squeeze(sum(vActivePolar,1,'omitnan'));  % [theta x tempo]
%  
%     %% PLOT Z in POLAR COORDINATES
%     
%     fontSize = 16;
%     
%     %set up radial grids
%     dThetaPlotGrid = 22.5; % grid angle spacing
%     if dTheta<0; dThetaPlotGrid = dThetaPlotGrid*-1; end %cw vs ccw
%     thetaPlotGrid = (thetaLims(1):dThetaPlotGrid:thetaLims(2)).*(2*pi/360); 
%     thetaPlotGridFine = (thetaLims(1):dThetaPlotGrid/2:thetaLims(2)).*(2*pi/360); % where 0 deg =  East and direction is CCW
%     rhoPlotGrid = linspace(0,rhoLims(2),4);
%     [xPolarPlotGrid,yPolarPlotGrid] = pol2cart(repmat(thetaPlotGrid,length(rhoPlotGrid),1),...
%         repmat(rhoPlotGrid,length(thetaPlotGrid),1)');
%     xPolarPlotGrid = xPolarPlotGrid./1000; %convert to km
%     yPolarPlotGrid = yPolarPlotGrid./1000;
%     [xPolarPlotGridFine,yPolarPlotGridFine] = pol2cart(repmat(thetaPlotGridFine,length(rhoPlotGrid),1),...
%         repmat(rhoPlotGrid,length(thetaPlotGridFine),1)');
%     xPolarPlotGridFine = xPolarPlotGridFine./1000; %convert to km
%     yPolarPlotGridFine = yPolarPlotGridFine./1000;
%     
%     thetaDegNorth = theta./(2*pi/360); % convert back from radians to degrees
%     thetaDegNorth = 90-thetaDegNorth; % correct to 0=north, CW pos (nautical)
%     [thetaGrid,rhoGrid] = meshgrid(thetaDegNorth(1:end-1),rho(1:end-1)./1000);
% 
%     save(fullfile(outPath,[originName{oo} '_polarData.mat']),'vActivePolar','hActivePolar','vAnomalyTheta','vAnomalyRho',...
%         'zPolar','hAnomalyPolar','vAnomalyPolar','vActive','hActive',...
%         'weightedMeanSurf','theta','rho','xPolar','yPolar','zPolar','T','rhoLims','thetaLims',...
%         'thetaDegNorth', 'thetaGrid', 'rhoGrid', 'xPolarPlotGrid', 'yPolarPlotGrid', 'xPolarPlotGridFine',...
%         'yPolarPlotGridFine','fontSize','origins','originName','outPathRoot','outPath');
%     
%     % ---------------------------------------------------------------------
%     % PLOT RESAMPLED DATA IN REAL XY SPACE
%     % ===== Discrete colormap e quebras (iguais à sua legenda) =====
%     edges = [-9.8 -8 -7 -6 -5 -4 -3 -2 -1 0 1 2 3];  % 12 classes
%     cmap = [
%         0.50 0.00 0.00  %  2.1–3   (vermelho escuro)
%         0.80 0.00 0.00  %  1.1–2   (vermelho)
%         1.00 0.50 0.20  %  0.1–1   (laranja)
%         1.00 0.70 0.40  % -0.9–0   (laranja claro)
%         1.00 0.88 0.55  % -1.9–-1  (amarelo quente)
%         1.00 0.97 0.80  % -2.9–-2  (amarelo muito claro)
%         0.80 0.92 1.00  % -3.9–-3  (azul muito claro)
%         0.60 0.80 1.00  % -4.9–-4  (azul claro)
%         0.40 0.60 0.85  % -5.9–-5  (azul)
%         0.25 0.45 0.85  % -6.9–-6  (azul médio)
%         0.10 0.25 0.55  % -7.9–-7  (azul escuro)
%         0.00 0.08 0.40  % -9.8–-8  (azul bem escuro)
%     ];
%     cmap = flipud(cmap);
%     tickCenters = (edges(1:end-1)+edges(2:end))/2;
%     tickLabels = {'-9.8–-8','-7.9–-7','-6.9–-6','-5.9–-5','-4.9–-4','-3.9–-3',...
%               '-2.9–-2','-1.9–-1','-0.9–0','0.1–1','1.1–2','2.1–3'};
%     % ===== Colormap azul → branco → vermelho (simétrico) =====
%     n = 256;
%     r = [linspace(0,1,n/2)'; ones(n/2,1)];
%     g = [linspace(0,1,n/2)'; linspace(1,0,n/2)'];
%     b = [ones(n/2,1); linspace(1,0,n/2)'];
%     cmapBlueWhiteRed = [b g r];
% 
%         
%         for tt = 1:length(T)
% %             f = figure(1); f.WindowState='maximized';
% %             clf
% %             tlo1 = tiledlayout(1,2); tlo1.Padding='tight'; tlo1.TileSpacing='tight';
%               f = figure(1); 
%               set(f,'WindowState','normal')   % importante
%               set(f,'Units','centimeters')
%               set(f,'Position',[2 2 21 17])
%               tlo = tiledlayout(1,2);
%               tlo.TileSpacing = 'tight'; 
%               tlo.Padding = 'tight';
% 
%     % -------------------- XY --------------------
%             ax(1) = nexttile; hold on; grid on; box on;
%             xp = xPolar(1:end-1,1:end-1); yp = yPolar(1:end-1,1:end-1);
% 
%     % zPolar em classes discretas
%             contourf(xp, yp, zPolar(:,:,tt), edges, 'LineStyle','none');
%             colormap(ax(1), cmap); caxis(ax(1), [edges(1) edges(end)]);
%             xlabel(ax(1),'X (km)'); ylabel(ax(1),'Y (km)');
%             title(ax(1),['X-Y Space [' num2str(T(tt)) ']']);
%             set(ax(1),'FontName','Myriad Pro','FontSize',fontSize);
%             axis square
% 
%             % === contorno de 0 m ===
%             contour(xp, yp, zPolar(:,:,tt), [0 0], 'k-', 'LineWidth', 1.2);
% 
%     % grid radial por cima
%             ii = 1:length(thetaPlotGrid);
%             plot(xPolarPlotGrid(:,ii), yPolarPlotGrid(:,ii), '-k');
%             for ii = 1:length(rhoPlotGrid)
%                 plot(xPolarPlotGridFine(ii,:), yPolarPlotGridFine(ii,:), '-k');
%             end
% 
%     % -------------------- ρ–θ --------------------
%             ax(2) = nexttile; hold on; grid on; box on;
%             contourf(thetaGrid, rhoGrid, zPolar(:,:,tt), edges, 'LineStyle','none');
%             colormap(ax(2), cmap); caxis(ax(2), [edges(1) edges(end)]);
%             contour(thetaGrid, rhoGrid, zPolar(:,:,tt), [0 0], 'k-', 'LineWidth', 1.2);
%             xlabel(ax(2),'\theta ({\circ})'); ylabel(ax(2),'\rho (km)');
%             set(ax(2),'XTick',min(thetaDegNorth):30:max(thetaDegNorth));
%             set(ax(2),'Layer','top','GridColor','k','GridAlpha',0.4,'FontName','Myriad Pro','FontSize',fontSize);
%             xlim(ax(2),[min(thetaDegNorth) max(thetaDegNorth)]);
%             ylim(ax(2),[min(rho./1000) max(rho./1000)]);
%             axis square
% 
%         % colorbar coerente nas duas views
%             c = colorbar(ax(2),'Location','eastoutside');
%             c.Ticks = tickCenters; c.TickLabels = tickLabels;
%             ylabel(c,'Elevation (m)');
% 
%     % (opcional) se quiser a ordem do raso->profundo de cima pra baixo:
%     % c.Direction = 'reverse';
% 
%             title(ax(2),['\rho-\theta Space [' num2str(T(tt)) ']']);
%             set(ax(2),'FontName','Myriad Pro','FontSize',fontSize);
% 
%     % Export
%             set(gcf,'Units','centimeters')
%             set(gcf,'Position',[2 2 21 17])
%             set(findall(gcf,'Type','axes'),'FontName','Arial','FontSize',12)
% 
%             pngName = [originName{oo} '_polarMap_Overview_' num2str(T(tt))];
%             savefig([outPath filesep pngName '.fig']);
%             exportgraphics(gcf,[outPath filesep pngName '.png'],"Resolution",300);
%         end
%         
%     
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         if makeMiddlePlots
%             % ---------------------------------------------------------------------
% 
%             % PLOT hAnomaly Data in XY and Polar Space
% 
%             for tt = 1:length(T)
%                 f = figure(2); f.WindowState='maximized';
%                 tlo=tiledlayout(1,2);tlo.TileSpacing='tight';tlo.Padding='tight';
%                 ax(1) = nexttile;
%                 hold on; grid on; box on;
%                 xp = xPolar(1:end-1,1:end-1); yp = yPolar(1:end-1,1:end-1);
%                 pc2=pcolor(xp,yp,hAnomalyPolar(:,:,tt));            
%                 pc2.FaceColor = 'flat'; pc2.EdgeColor = 'none';
%                 xlabel(ax(1),'X [km]'); ylabel(ax(1),'Y [km]');
%                 title(ax(1),['X-Y Space [' num2str(T(tt)) ']']);
%                 colormap(ax(1),cmap);
%                 axis square
%                 set(ax(1),'FontName','Myriad Pro','FontSize',fontSize);
%                 for ii = 1:length(thetaPlotGrid)
%                     plot(xPolarPlotGrid(:,ii),yPolarPlotGrid(:,ii),'-k'); % plot radial lines  
%                 end
%                 for ii = 1:length(rhoPlotGrid)
%                     plot(xPolarPlotGridFine(ii,:),yPolarPlotGridFine(ii,:),'-k'); % plot circles
%                 end
% 
% 
%                 % ---------------------------------------------------------------------
%                 % PLOT IN THETA-RHO SPACE
% 
%                 ax(2)=nexttile;
%                 hold on; grid on; box on;
%                 pc2=pcolor(thetaGrid,rhoGrid,hAnomalyPolar(:,:,tt));
%                 pc2.FaceColor = 'flat'; pc2.EdgeColor = 'none';
%                 xlabel(ax(2),'\theta [{\circ}]'); ylabel(ax(2),'\rho [km]');
%                 %if thetaLims(1) < thetaLims(2); set(ax(2),'XDir','reverse'); end
%                 set(ax(2),'XTick',[min(thetaDegNorth):30:max(thetaDegNorth)]);
%                 set(ax(2),'Layer','top','GridColor','k','GridAlpha',0.4);
%                 xlim(ax(2),[min(thetaDegNorth) max(thetaDegNorth)]);
%                 ylim(ax(2),[min(rho./1000) max(rho./1000)]);
%                 colormap(ax(2),cmap); cb = colorbar;
%                 ylabel(cb,'Deviation from Z_{mean} (m)','FontName','Myriad Pro','FontSize',fontSize);
%                 minC = min(hAnomalyPolar,[],'all'); maxC=max(hAnomalyPolar,[],'all');
%                 ax(1).CLim = [-maxC maxC]; ax(2).CLim = [-maxC maxC];
%                 title(ax(2),['\rho-\theta Space [' num2str(T(tt)) ']']);
%                 set(ax(2),'FontName','Myriad Pro','FontSize',fontSize);
%                 axis square
% 
%                 % Export figure to png file
%                 set(gcf,'Units','centimeters')
%                 set(gcf,'Position',[2 2 20 16])
%     
%                 set(findall(gcf,'Type','axes'),'FontName','Arial','FontSize',12)    
% 
% 
%                 pngName = [originName{oo} '_polarMap_hAnomaly_' num2str(T(tt))];
%                 savefig([outPath filesep pngName '.fig']);
%                 exportgraphics(gcf,[outPath filesep pngName '.png'],"Resolution",300);
%             end
% 
% 
%             
%             
%             % ---------------------------------------------------------------------
%             % PLOT  vAnomaly Data in XY and Polar Space
%             
%             for tt = 1:length(T)
%                 f = figure(3); 
%                 set(f,'WindowState','normal')   % importante
%                 set(f,'Units','centimeters')
%                 set(f,'Position',[2 2 21 17])
%                 tlo = tiledlayout(1,2);
%                 tlo.TileSpacing = 'tight'; 
%                 tlo.Padding = 'tight';
%                 %f = figure(3); f.WindowState='maximized';
%                 %tlo = tiledlayout(1,2); tlo.TileSpacing='tight';tlo.Padding='tight';
%                 ax(1) = nexttile;
%                 hold on; grid on; box on;
%                 xp = xPolar(1:end-1,1:end-1); yp = yPolar(1:end-1,1:end-1);
%                 pc2=pcolor(xp,yp,vAnomalyPolar(:,:,tt));  %TO DO make zeros nans in the no data areas
%                 pc2.FaceColor = 'flat'; pc2.EdgeColor = 'none';
%                 xlabel(ax(1),'X (km)'); ylabel(ax(1),'Y (km)');
%                 title(ax(1),['X-Y Space [' num2str(T(tt)) ']']);
%                 set(ax(1),'FontName','Myriad Pro','FontSize',fontSize);
%                 colormap(ax(1),cmap);
%                 axis square
%                 ii = 1:length(thetaPlotGrid);
%                 plot(xPolarPlotGrid(:,ii),yPolarPlotGrid(:,ii),'-k'); % plot radial lines
%                 for ii = 1:length(rhoPlotGrid)
%                     plot(xPolarPlotGridFine(ii,:),yPolarPlotGridFine(ii,:),'-k'); % plot circles
%                 end
%                 
%                 % ---------------------------------------------------------------------
%                 % PLOT IN THETA-RHO SPACE
%                           
%                 ax(2) = nexttile;
%                 hold on; grid on; box on;
%                 pc2=pcolor(thetaGrid,rhoGrid,vAnomalyPolar(:,:,tt)); 
%                 pc2.FaceColor = 'flat'; pc2.EdgeColor = 'none';
%                 xlabel(ax(2),'\theta ({\circ})'); ylabel(ax(2),'\rho (km)');
%                 %if thetaLims(1) < thetaLims(2); set(gca,'XDir','reverse'); end %make theta direction more like nautical convention
%                 set(ax(2),'XTick',[min(thetaDegNorth):30:max(thetaDegNorth)]);
%                 set(ax(2),'Layer','top','GridColor','k','GridAlpha',1);
%                 xlim(ax(2),[min(thetaDegNorth) max(thetaDegNorth)]);
%                 ylim(ax(2),[min(rho./1000) max(rho./1000)]); %km
%                 %colormap(parula);
%                 colormap(ax(2),cmap);
%                 cb=colorbar;
%                 ylabel(cb,'V_{Anomaly} (m^3)','FontName','Myriad Pro','FontSize',fontSize);
%                 minC = min(vAnomalyPolar,[],'all'); maxC=max(vAnomalyPolar,[],'all');
%                 ax(1).CLim = [-maxC maxC]; 
%                 ax(2).CLim = [-maxC maxC];
%                 title(['\rho-\theta Space [' num2str(T(tt)) ']']);
%                 set(ax(2),'FontName','Myriad Pro','FontSize',fontSize);
%                 axis square
%                 
%                 % Export figure to png file
%          
%                 set(findall(gcf,'Type','axes'),'FontName','Arial','FontSize',12)
%                 pngName = [originName{oo} '_polarMap_vAnomaly_' num2str(T(tt))];
%                 savefig([outPath filesep pngName '.fig']);
%                 exportgraphics(gcf,[outPath filesep pngName '.png'],"Resolution",300);
%             end
%                                 % ---------------------------------------------------------------------
%             % PLOT vActive Data in XY and Polar Space
%             % (mudança de volume em relação ao ano de referência)
%     
%             % define limites de cor globais para deixar simétrico em torno de 0
%             minC = min(vActivePolar,[],'all','omitnan');
%             maxC = max(vActivePolar,[],'all','omitnan');
%             maxAbs = max(abs([minC maxC]));
%     
%             for tt = 1:length(T)
%                 f = figure; 
%                 f.WindowState = 'maximized';
%                 tlo = tiledlayout(1,2);
%                 tlo.TileSpacing = 'tight'; 
%                 tlo.Padding     = 'tight';
%     
%                 % -------------------- XY --------------------
%                 ax(1) = nexttile;
%                 hold on; grid on; box on;
%     
%                 xp = xPolar(1:end-1,1:end-1); 
%                 yp = yPolar(1:end-1,1:end-1);
%     
%                 pc2 = pcolor(xp,yp,vActivePolar(:,:,tt));
%                 pc2.FaceColor = 'flat'; 
%                 pc2.EdgeColor = 'none';
%     
%                 xlabel(ax(1),'X [km]');
%                 ylabel(ax(1),'Y [km]');
%                 title(ax(1),['X-Y Space V_{active} [' num2str(T(tt)) ']']);
%     
%                 % usa a paleta divergente que vc já criou (cmapHov),
%                 % ou troque por "cmap" se preferir igual ao vAnomaly
%                 colormap(ax(1),cmap);
%                 ax(1).CLim = [-maxAbs maxAbs];
%                 axis square
%                 set(ax(1),'FontName','Myriad Pro','FontSize',fontSize);
%     
%                 % plota o grid polar por cima
%                 for ii = 1:length(thetaPlotGrid)
%                     plot(xPolarPlotGrid(:,ii),yPolarPlotGrid(:,ii),'-k'); 
%                 end
%                 for ii = 1:length(rhoPlotGrid)
%                     plot(xPolarPlotGridFine(ii,:),yPolarPlotGridFine(ii,:),'-k'); 
%                 end
%     
%                 % -------------------- ρ–θ --------------------
%                 ax(2) = nexttile;
%                 hold on; grid on; box on;
%     
%                 pc2 = pcolor(thetaGrid, rhoGrid, vActivePolar(:,:,tt));
%                 pc2.FaceColor = 'flat'; 
%                 pc2.EdgeColor = 'none';
%     
%                 xlabel(ax(2),'\theta [{\circ}]');
%                 ylabel(ax(2),'\rho [km]');
%                 title(ax(2),['\rho-\theta Space V_{active} [' num2str(T(tt)) ']']);
%     
%                 % ticks e limites
%                 set(ax(2),'XTick',min(thetaDegNorth):30:max(thetaDegNorth));
%                 xlim(ax(2),[min(thetaDegNorth) max(thetaDegNorth)]);
%                 ylim(ax(2),[min(rho./1000) max(rho./1000)]);
%                 set(ax(2),'Layer','top','GridColor','k','GridAlpha',0.4,...
%                           'FontName','Myriad Pro','FontSize',fontSize);
%                 axis square
%     
%                 % mesma escala de cor dos dois painéis
%                 colormap(ax(2),cmap);
%                 ax(2).CLim = [-maxAbs maxAbs];
%     
%                 cb = colorbar(ax(2),'Location','eastoutside');
%                 ylabel(cb,'V_{active} [m^3]',...
%                     'FontName','Myriad Pro','FontSize',fontSize);
%     
%                 % salvar figura
%                 pngName = [originName{oo} '_polarMap_vActive_' num2str(T(tt))];
%                 savefig([outPath filesep pngName '.fig']);
%                 exportgraphics(gcf,[outPath filesep pngName '.png'],"Resolution",300);
%             end
% 
%             
%         end
%     
% %% ===== Colormap divergente contínua (azul < 0, branco = 0, vermelho > 0)
% % com mais pontos de controle perto de 0 para transição rápida =====
% v = [-1 -0.9 -0.7 -0.5 -0.3 -0.15 -0.07 -0.03 0 0.03 0.07 0.15 0.3 0.5 0.7 0.9 1]';
% cols = [ ...
%     0.00 0.08 0.40;  % azul bem escuro
%     0.10 0.25 0.55;
%     0.25 0.45 0.85;
%     0.40 0.60 0.85;
%     0.60 0.80 1.00;
%     0.80 0.92 1.00;  % azul muito claro
%     0.90 0.95 1.00;  % quase branco (lado <0)
%     0.97 0.98 1.00;  % quase branco
%     1.00 1.00 1.00;  % branco (zero)
%     1.00 0.98 0.97;  % quase branco (lado >0)
%     1.00 0.95 0.90;
%     1.00 0.88 0.55;  % amarelo quente
%     1.00 0.70 0.40;  % laranja
%     1.00 0.50 0.20;  % laranja forte
%     0.80 0.00 0.00;  % vermelho
%     0.50 0.00 0.00;  % vermelho escuro
%     0.40 0.00 0.00]; % extremo
% xi = linspace(-1,1,256)';                 % amostragem da paleta
% cmapHov = interp1(v, cols, xi, 'pchip');  % interpola com mais densidade perto de 0
% 
% %% Hovmoller (2D timeseries) Diagram
% 
% 
% % Add extra year so that the last year of data is visible  
% T_hov = [T T(end)+1];
% 
% %vAnomaly is in m^3
% f = figure(4); f.WindowState='maximized';
% tlo = tiledlayout(1,2);
% tlo.TileSpacing = 'tight'; tlo.Padding = 'tight';
% cmax = max(abs([vAnomalyTheta(:); vAnomalyRho(:)]));
% 
% % Theta
% ax(1) = nexttile;
% hold on; box on; grid on;
% pp = pcolor(repmat(thetaDegNorth(1:end-1)',1,length(T_hov)), ...
%             repmat(T_hov,length(thetaDegNorth(1:end-1)),1), ...
%             [vAnomalyTheta NaN*ones(length(vAnomalyTheta),1)]);
% pp.EdgeColor = 'none';
% caxis(ax(1), [-cmax cmax]);
% colormap(ax(1), cmapHov);            % << usa a paleta contínua
% axis square;
% title(ax(1),'Changes in Volume Anomaly [\theta]')
% xlabel(ax(1),'\theta [\circ]'); ylabel(ax(1),'Time [y]');
% xlim(ax(1),[min(thetaDegNorth) max(thetaDegNorth)]); 
% ylim(ax(1),[T_hov(1) T_hov(end)]);
% set(ax(1),'XTick',min(thetaDegNorth):30:max(thetaDegNorth));
% set(ax(1),'Layer','top','GridColor','k','GridAlpha',0.4);
% set(ax(1),'FontName','Myriad Pro','FontSize',fontSize)
% 
% 
% % Rho
% ax(2) = nexttile;
% hold on; box on; grid on;
% pp = pcolor(repmat(T_hov,length(rho(1:end-1)),1), ...
%             repmat(rho(1:end-1)./1000,length(T_hov),1)', ...
%             [vAnomalyRho NaN*ones(length(vAnomalyRho),1)]);
% pp.EdgeColor = 'none';
% caxis(ax(2), [-cmax cmax]);
% colormap(ax(2), cmapHov);            % << mesma paleta
% 
% axis square
% cb = colorbar('Location','eastoutside');
% ylabel(cb,'V_{anomaly} [m^3]','FontName','Myriad Pro','FontSize',fontSize,'FontWeight','bold')
% title(ax(2),'Changes in Volume Anomaly [\rho]')
% ylabel(ax(2),'\rho [km]'); xlabel(ax(2),'Time [y]');
% set(ax(2),'Layer','top','GridColor','k','GridAlpha',0.4);
% set(ax(2),'FontName','Myriad Pro','FontSize',fontSize)
% ylim(ax(2),[rhoLims(1) rhoLims(2)]./1000);
% xlim(ax(2),[T_hov(1) T_hov(end)]);
% 
% 
% 
% %Save figure
% set(gcf,'Units','centimeters')
% set(gcf,'Position',[2 2 20 16])
% set(findall(gcf,'Type','axes'),'FontName','Arial','FontSize',12)
% pngName = [originName{oo} '_HovmullerPolarRhoAndTheta'];
% savefig([outPath filesep pngName '.fig']);
% exportgraphics(gcf,[outPath filesep pngName '.png'],"Resolution",300);
% 
% 
% 
% %% Hovmoller (2D timeseries) Diagram
% 
% % ------------------------------------------------------------
% % T = grd.year  -> anos EM QUE TEM DADO (ex.: [2004 2006 2010 2014 2016 2022])
% % ------------------------------------------------------------
% Tdata = T(:)';                 % garante vetor linha
% Tfull = Tdata(1):Tdata(end);   % anos contínuos, ex.: 2004:2022
% T_hov = [Tfull Tfull(end)+1];  % bordas para o pcolor
% 
% % vAnomalyTheta  -> [Ntheta-1  x  Ndata]   (linhas = espaço, colunas = tempo)
% % vAnomalyRho    -> [Nrho-1    x  Ndata]
% 
% %yearsToShow = [2004  2006  2010  2015 2022 ];
% 
% % Cria matrizes no MESMO formato, mas com todos os anos de Tfull
% vTheta_full = NaN(size(vAnomalyTheta,1), numel(Tfull));  % [Ntheta-1 x Nfull]
% vRho_full   = NaN(size(vAnomalyRho,1),   numel(Tfull));  % [Nrho-1   x Nfull]
% 
% % Copia as colunas correspondentes aos anos em que há dado
% for k = 1:numel(Tdata)
%     idx = find(Tfull == Tdata(k));     % posição desse ano em Tfull
%     if ~isempty(idx)
%         vTheta_full(:,idx) = vAnomalyTheta(:,k);
%         vRho_full(:,idx)   = vAnomalyRho(:,k);
%     end
% end
% 
% % Escala de cores comum
% cmax = max(abs([vTheta_full(:); vRho_full(:)]));
% 
% % Figura
% f = figure(4); 
% set(f,'WindowState','normal')   % importante
% set(f,'Units','centimeters')
% set(f,'Position',[2 2 21 17])
% 
% tlo = tiledlayout(1,2);
% tlo.TileSpacing = 'tight'; 
% tlo.Padding = 'tight';
% 
% %% ======================  THETA  ======================================
% ax(1) = nexttile;
% hold on; box on; grid on;
% 
% % X: theta, Y: tempo, Z: vTheta_full (+ coluna extra de NaN no fim)
% pp = pcolor( ...
%     repmat(thetaDegNorth(1:end-1)', 1, length(T_hov)), ...
%     repmat(T_hov, length(thetaDegNorth(1:end-1)), 1), ...
%     [vTheta_full NaN(size(vTheta_full,1),1)] );
% pp.EdgeColor = 'none';
% 
% caxis(ax(1), [-cmax cmax]);
% colormap(ax(1), cmapHov);      % paleta contínua
% axis square;
% 
% title(ax(1), 'Changes in Volume Anomaly [\theta]')
% xlabel(ax(1), '\theta (\circ)'); 
% ylabel(ax(1), 'Year');
% 
% xlim(ax(1), [min(thetaDegNorth) max(thetaDegNorth)]);
% ylim(ax(1), [T_hov(1) T_hov(end)]);
% 
% set(ax(1), 'XTick', min(thetaDegNorth):30:max(thetaDegNorth));
% set(ax(1), 'YTick', Tfull);           % ticks em cada ano
% %set(ax(1),'YTick',yearsToShow)
% set(ax(1), 'Layer','top','GridColor','k','GridAlpha',0.4);
% set(ax(1), 'FontName','Myriad Pro','FontSize',fontSize);
% 
% % Se quiser contorno da isolinha zero:
% % [Xt,Yt] = deal(repmat(thetaDegNorth(1:end-1)',1,length(T_hov)), ...
% %                repmat(T_hov,length(thetaDegNorth(1:end-1)),1));
% % Zt = [vTheta_full NaN(size(vTheta_full,1),1)];
% % contour(Xt, Yt, Zt, [0 0], 'k-', 'LineWidth', 1.0);
% 
% %% ======================  RHO  ========================================
% ax(2) = nexttile;
% hold on; box on; grid on;
% 
% % X: tempo, Y: rho (km), Z: vRho_full (+ coluna extra de NaN no fim)
% pp = pcolor( ...
%     repmat(T_hov, length(rho(1:end-1)), 1), ...
%     repmat(rho(1:end-1)./1000, length(T_hov), 1)', ...
%     [vRho_full NaN(size(vRho_full,1),1)] );
% pp.EdgeColor = 'none';
% 
% caxis(ax(2), [-cmax cmax]);
% colormap(ax(2), cmapHov);      % mesma paleta
% 
% axis square;
% cb = colorbar('Location','eastoutside');
% ylabel(cb, 'V_{anomaly} (m^3)', ...
%     'FontName','Myriad Pro','FontSize',fontSize,'FontWeight','bold');
% 
% title(ax(2), 'Changes in Volume Anomaly (\rho)')
% ylabel(ax(2), '\rho (km)'); 
% xlabel(ax(2), 'Year');
% 
% set(ax(2), 'XTick', Tfull);
% %set(ax(2),'XTick',yearsToShow)
% set(ax(2), 'Layer','top','GridColor','k','GridAlpha',0.4);
% set(ax(2), 'FontName','Myriad Pro','FontSize',fontSize);
% 
% ylim(ax(2), [rhoLims(1) rhoLims(2)]./1000);
% xlim(ax(2), [T_hov(1) T_hov(end)]);
% 
% % Se quiser contorno da isolinha zero:
% % [Xr,Yr] = deal(repmat(T_hov,length(rho(1:end-1)),1), ...
% %                repmat(rho(1:end-1)./1000,length(T_hov),1)');
% % Zr = [vRho_full NaN(size(vRho_full,1),1)];
% % contour(Xr, Yr, Zr, [0 0], 'k-', 'LineWidth', 1.0);
% 
% %% =====================  SALVAR FIGURA  ===============================
% set(gcf,'WindowState','normal')   % importante
% set(gcf,'Units','centimeters')
% set(gcf,'Position',[2 2 21 17])
% 
% set(findall(gcf,'Type','axes'),'FontName','Arial','FontSize',11)
% set(findall(gcf,'Type','ColorBar'),'FontName','Arial','FontSize',12)
% set(findall(gcf,'Type','text'),'FontName','Arial','FontSize',12)
% 
% pngName = [originName{oo} '_YEAR_HovmullerPolarRhoAndTheta'];
% savefig([outPath filesep pngName '.fig']);
% exportgraphics(gcf,[outPath filesep pngName '.png'],'Resolution',300);
% 
% 
% %%
% 
% 
% 
% %%
% % %% Hovmoller (2D timeseries) Diagram - V_ACTIVE (relativo ao 1º survey)
% % 
% % % Reutiliza Tdata, Tfull e T_hov definidos acima
% % % Tdata = T(:)';                 
% % % Tfull = Tdata(1):Tdata(end);   
% % % T_hov = [Tfull Tfull(end)+1];
% % 
% % % vActiveTheta  -> [Ntheta-1 x Ndata]
% % % vActiveRho    -> [Nrho-1   x Ndata]
% % 
% % vThetaActive_full = NaN(size(vActiveTheta,1), numel(Tfull));  % [Ntheta-1 x Nfull]
% % vRhoActive_full   = NaN(size(vActiveRho,1),   numel(Tfull));  % [Nrho-1   x Nfull]
% % 
% % % Copia as colunas correspondentes aos anos em que há dado
% % for k = 1:numel(Tdata)
% %     idx = find(Tfull == Tdata(k));
% %     if ~isempty(idx)
% %         vThetaActive_full(:,idx) = vActiveTheta(:,k);
% %         vRhoActive_full(:,idx)   = vActiveRho(:,k);
% %     end
% % end
% % 
% % % Escala de cores comum
% % cmaxA = max(abs([vThetaActive_full(:); vRhoActive_full(:)]));
% % 
% % % Figura
% % f = figure; 
% % f.WindowState = 'maximized';
% % tlo = tiledlayout(1,2);
% % tlo.TileSpacing = 'tight'; 
% % tlo.Padding = 'tight';
% % 
% % %% ========  THETA (V_ACTIVE)  =========
% % ax(1) = nexttile;
% % hold on; box on; grid on;
% % 
% % pp = pcolor( ...
% %     repmat(thetaDegNorth(1:end-1)', 1, length(T_hov)), ...
% %     repmat(T_hov, length(thetaDegNorth(1:end-1)), 1), ...
% %     [vThetaActive_full NaN(size(vThetaActive_full,1),1)] );
% % pp.EdgeColor = 'none';
% % 
% % caxis(ax(1), [-cmaxA cmaxA]);
% % colormap(ax(1), cmapHov);      % mesma paleta que vc já criou
% % axis square;
% % 
% % title(ax(1), 'Changes in ACTIVE Volume [\theta]')
% % xlabel(ax(1), '\theta [\circ]'); 
% % ylabel(ax(1), 'Time [y]');
% % 
% % xlim(ax(1), [min(thetaDegNorth) max(thetaDegNorth)]);
% % ylim(ax(1), [T_hov(1) T_hov(end)]);
% % 
% % set(ax(1), 'XTick', min(thetaDegNorth):30:max(thetaDegNorth));
% % set(ax(1), 'YTick', Tfull);
% % set(ax(1), 'Layer','top','GridColor','k','GridAlpha',0.4);
% % set(ax(1), 'FontName','Myriad Pro','FontSize',fontSize);
% % 
% % %% ========  RHO (V_ACTIVE)  ===========
% % ax(2) = nexttile;
% % hold on; box on; grid on;
% % 
% % pp = pcolor( ...
% %     repmat(T_hov, length(rho(1:end-1)), 1), ...
% %     repmat(rho(1:end-1)./1000, length(T_hov), 1)', ...
% %     [vRhoActive_full NaN(size(vRhoActive_full,1),1)] );
% % pp.EdgeColor = 'none';
% % 
% % caxis(ax(2), [-cmaxA cmaxA]);
% % colormap(ax(2), cmapHov);      
% % 
% % axis square;
% % cb = colorbar('Location','eastoutside');
% % ylabel(cb, 'V_{active} [m^3]', ...
% %     'FontName','Myriad Pro','FontSize',fontSize,'FontWeight','bold');
% % 
% % title(ax(2), 'Changes in ACTIVE Volume [\rho]')
% % ylabel(ax(2), '\rho [km]'); 
% % xlabel(ax(2), 'Time [y]');
% % 
% % set(ax(2), 'XTick', Tfull);
% % set(ax(2), 'Layer','top','GridColor','k','GridAlpha',0.4);
% % set(ax(2), 'FontName','Myriad Pro','FontSize',fontSize);
% % 
% % ylim(ax(2), [rhoLims(1) rhoLims(2)]./1000);
% % xlim(ax(2), [T_hov(1) T_hov(end)]);
% % 
% % % Salvar figura
% % pngName = [originName{oo} '_YEAR_HovmullerPolarRhoAndTheta_ACTIVE'];
% % savefig([outPath filesep pngName '.fig']);
% % exportgraphics(gcf,[outPath filesep pngName '.png'], "Resolution", 300);
% 
% 
% 
% 
% %% Peakfinding Analysis
% 
% %Uses findpeaks which requires Signal Processing Toolbox
%     scaleFactor = 2;
%     minPeakDist = 8; %degrees apart
%     numPeaks = 4;
%     
%     figure
%     clf
%     for ii=1:length(T)
%         hold on
%         p(ii) = plot(thetaDegNorth(1:end-1),vAnomalyTheta(:,ii),'DisplayName',num2str(T(ii)));
%         [peaks(ii).pks,peaks(ii).locs,peaks(ii).width,peaks(ii).prom] = findpeaks(vAnomalyTheta(:,ii),...
%             'MinPeakDistance',minPeakDist,'NPeaks',numPeaks,'SortStr','descend');
%         s(ii) = scatter(thetaDegNorth(peaks(ii).locs),peaks(ii).pks,'o','filled');
%         s(ii).MarkerFaceColor = p(ii).Color;
%         text(thetaDegNorth(peaks(ii).locs)+.02,peaks(ii).pks+T(ii),num2str((1:numel(peaks(ii).pks))')); %these are numbered descending now
%     end
%     grid on
%     legend(p);
%     title('Peaks in Volume Anomaly [\theta]')
%     xlabel('\theta [\circ]');
%     ylabel('V. Anomaly \theta (m^3)');
%     % ylim([1975 2021]);
%     %xlim([-100 79]);
%     %set(gca,'XTick',[-180:22.5:180])
%     set(gca,'Layer','top','GridColor','k','GridAlpha',0.4);
%     set(gca,'FontName','Myriad Pro','FontSize',fontSize)
%     
% 
% 
% %This does not require the Signal Processing Toolbox
% x = thetaDegNorth(1:end-1);
% 
% figure('windowstate','maximized');
% clf
% hold on
% MinSeparation = 8; 
% MaxNumExtrema = 4; 
% 
% for ii = 1:size(vAnomalyTheta,2)
%     [TFmax,Pmax(:,ii)] = islocalmax(vAnomalyTheta(:,ii),...
%         'MinSeparation',MinSeparation,...
%         'MaxNumExtrema',MaxNumExtrema);
%     y = [repmat(T(ii),1,length(x))];
%     z = vAnomalyTheta(:,ii);
%     plot3(x,y,z);
%     scatter3(x(TFmax),y(TFmax),z(TFmax),'*r'); 
%     %peakLocs(:,ii) = x(TFmax);
% end 
% view(-7,57)
% grid on
% title('Changes in Volume Anomaly [\theta]')
% xlabel('\theta [\circ]');
% ylabel('Time [y]');
% zlabel('Volume Anomaly Theta (m^3)');
% % Export figure to png file
% pngName = [originName{oo} '_Vol_Anom_Theta_TimeSlices_Peakfinder'];
% savefig([outPath filesep pngName '.fig']);
% exportgraphics(gcf,[outPath filesep pngName '.png'],"Resolution",300);
% 
% 
% 
% %% Troughfinding analysis
% % figure('windowstate','maximized');
% % clf
% % hold on
% % for ii = 1:size(vAnomalyTheta,2)
% %      [TFmin,Pmin(:,ii)]  = islocalmin(vAnomalyTheta(:,ii),'MinSeparation',MinSeparation,...
% %         'MaxNumExtrema',MaxNumExtrema);
% %     y = [repmat(T(ii),1,length(x))];
% %     z = vAnomalyTheta(:,ii);
% %     plot3(x,y,z);
% %     scatter3(x(TFmin),y(TFmin),z(TFmin),'*r');
% %     troughLocs(:,ii) = x(TFmin);
% % end
% % view(-7,57)
% % grid on
% % title('Changes in Volume Anomaly [\theta]')
% % xlabel('\theta [\circ]');
% % ylabel('Time [y]');
% % zlabel('Volume Anomaly Theta (m^3)');
% % % Export figure to png file
% % pngName = [originName{oo} '_Vol_Anom_Theta_TimeSlices_Troughfinder'];
% % savefig([outPath filesep pngName '.fig']);
% % exportgraphics(gcf,[outPath filesep pngName '.png'],"Resolution",300);
% 
% 
% 
% %% Peak Plotter 
% %From Pearson et al 2022 Supplementary Material 
% %Ridges and troughs were identified in the Volume Anomaly timestack
% %(Hovmoller diagram) using findpeaks. Peaks higher than 0.1x10^6 m^3 and
% %further than 8 deg. apart at a given timestep were selected, discounting
% %endpoints. Peaks were manually connected and a linear regression was applied to
% %estimate the rotation migration rates per shoal and channel.
% %WIP! 
% 
% 
% % 
% % %Plot Hovmoller Theta with peaks and fit lines     
% % figure('windowstate','maximized');
% % clf 
% % hold on; box on; grid on;
% % %vol anom
% % pp = pcolor(repmat(x',1,length(T)),repmat(T,length(x),1),vAnomalyTheta);
% % pp.EdgeColor='none'; colormap(parula);
% % axis square; title('Peak Tracking')
% % xlabel('\theta [\circ]'); ylabel('Time [y]');
% % xlim([min(thetaDegNorth) max(thetaDegNorth)]); ylim([T(1) T(end)]);
% % set(gca,'XTick',[min(thetaDegNorth):30:max(thetaDegNorth)]);
% % set(gca,'Layer','top','GridColor','k','GridAlpha',0.4);
% % set(gca,'FontName','Myriad Pro','FontSize',fontSize);
% % cb=colorbar('Location','eastoutside');
% % ylabel(cb,'V_{anomaly} [m^3]','FontName','Myriad Pro','FontSize',fontSize,'fontweight','bold')
% % %connect peaks with lines
% % x = thetaDegNorth(1:end-1);
% % for ii = 1:size(vAnomalyTheta,2)
% %     %find the peaks
% %     [TFmax,Pmax(:,ii)] = islocalmax(vAnomalyTheta(:,ii),'MinSeparation',MinSeparation,...
% %        'MaxNumExtrema',MaxNumExtrema); 
% %     peakLocs(:,ii) = x(TFmax); %peak location at theta (deg)
% %     scatter(x(TFmax),T(ii),'*k'); %plot the peaks
% %     if ii == 1 %label peak numbers
% %         xx = sort(x(TFmax),'ascend');
% %         yy = repmat(T(ii)+0.5,length(x(TFmax)),1);
% %         txt = 1:4;
% %         for jj = 1:MaxNumExtrema
% %             text(xx(jj)+5, yy(jj), num2str(txt(jj)),'FontSize',16); %label peaks left to right
% %         end
% %     end
% % end
% % % for ii = 1:size(peakLocs,2)-1
% % %     if abs(peakLocs() > 50 %if peak is more than 50 degrees away from previous peak
% % 
% % %TO DO if peak at T2 very different than 'same' peak at T1, flag as
% % %different peaks?
% % pPeaks = cell(size(peakLocs,1),1);
% % errorEst = cell(size(peakLocs,1),1);
% % fitPeaks = cell(size(peakLocs,1),1);
% % for ii = 1:size(peakLocs,1) 
% %     %linear fit line to each peak through time
% %     %peaks are rows cols are time
% %     [pPeaks{ii},errorEst{ii}] = polyfit(T,peakLocs(ii,:),1); %linear fit
% %     fitPeaks{ii} = polyval(pPeaks{ii},T,errorEst{ii});
% %     %plot black lines between peaks
% %     plot(peakLocs(ii,:),T,'-k')
% %     %plot red fit lines
% %     plot(fitPeaks{ii},T,'--r','LineWidth',2);
% % end
% % % % Export figure to png file
% % pngName = [originName{oo} '_HovmullerPolarTheta_Peakplotter'];
% % savefig([outPath filesep pngName '.fig']);
% % exportgraphics(gcf,[outPath filesep pngName '.png'],"Resolution",300);
% 
% 
% 
% 
% 
% 
% 
% %% Trough Plotter 
% 
% % 
% % %Plot Hovmoller Theta with peaks and fit lines     
% % MinSeparation = 8; %degrees
% % MaxNumExtrema = 4; %highest four peaks, how to choose for any site?
% % MinProminence = 0.1; %x10^6 m3
% % 
% % 
% % %TO DO if next trough too far away, ignore 
% % 
% % figure('windowstate','maximized');
% % clf 
% % hold on; box on; grid on;
% % %vol anom
% % pp = pcolor(repmat(thetaDegNorth(1:end-1)',1,length(T)),...
% %     repmat(T,length(thetaDegNorth(1:end-1)),1),...
% %     vAnomalyTheta);
% % pp.EdgeColor='none';
% % colormap(parula);
% % axis square
% % title('Trough Tracking')
% % xlabel('\theta [\circ]'); 
% % ylabel('Time [y]');
% % xlim([min(thetaDegNorth) max(thetaDegNorth)]); 
% % ylim([T(1) T(end)]);
% % set(gca,'XTick',[min(thetaDegNorth):30:max(thetaDegNorth)]);
% % set(gca,'Layer','top','GridColor','k','GridAlpha',0.4);
% % set(gca,'FontName','Myriad Pro','FontSize',fontSize);
% % cb=colorbar('Location','eastoutside');
% % ylabel(cb,'V_{anomaly} [m^3]','FontName','Myriad Pro','FontSize',fontSize,'fontweight','bold')
% % %lines
% % x = thetaDegNorth(1:end-1);
% % for ii = 1:size(vAnomalyTheta,2)
% %     [TFmin,Pmin(:,ii)]  = islocalmin(vAnomalyTheta(:,ii),'MinSeparation',MinSeparation,...
% %         'MaxNumExtrema',MaxNumExtrema); %,'MinProminence',MinProminence);
% %     troughLocs(:,ii) = x(TFmin); %trough location at theta (deg)
% %     scatter(x(TFmin),T(ii),'*k'); %plot the troughs
% % end
% % for ii = 1:size(troughLocs,1) 
% %     %TROUGHS
% %     [pTroughs{ii},errorEst{ii}] = polyfit(T,troughLocs(ii,:),1); %linear fit
% %     fitTroughs{ii} = polyval(pTroughs{ii},T,errorEst{ii});
% %     %plot black lines between peaks
% %     plot(troughLocs(ii,:),T,'-k')
% %     %plot red fit lines
% %     troughFitLine = plot(fitTroughs{ii},T,'--r','LineWidth',2);
% % end
% % % Export figure to png file
% % pngName = [originName{oo} '_HovmullerPolarTheta_Troughplotter'];
% % savefig([outPath filesep pngName '.fig']);
% % exportgraphics(gcf,[outPath filesep pngName '.png'],"Resolution",300);
% 
% 
% % 
% %% Peak Rotation Estimates
% % 
% %Output file of peak rotation estimates
% %Slopes of fit lines from polyfit
% 
% % for ii = 1:MaxNumExtrema
% %     peakFitLineSlope(ii) = round(pPeaks{ii}(1),2); %round for export
% % end
% % 
% % %TO DO MAKE SURE THESE ARE NUMBERED CORRECTLY
% % peakNums = (1:MaxNumExtrema)';
% % 
% % peakTable = table(peakNums,peakFitLineSlope');
% % peakTable.Properties.VariableNames = {'Peak Number','Slope of Fit Line (deg./yr)'};
% % 
% % fName = [originName{oo} '_Rotation_Rate_Estimates_Peaks.txt'];
% % writetable(peakTable,[outPath filesep fName])
% % % 
% % %TO DO add the R^2 to table
% % 
% % 
% % %TO DO plot on t plus one vanom theta plot and see if that screws things up
% % 
% % 
% % 
% %% Trough Rotation Estimates
% % 
% % 
% % for ii = 1:length(pTroughs)
% %     troughFitLineSlope(ii) = round(pTroughs{ii}(1),2); %round for export
% % end
% % 
% % 
% % troughNums = (1:MaxNumExtrema)';
% % 
% % troughTable = table(troughNums,troughFitLineSlope');
% % troughTable.Properties.VariableNames = {'Peak Number','Slope of Fit Line'};
% % 
% % header = {'Trough Number, Clockwise around inlet origin',...
% %     'Slope of Linear Fit Line for each Volume Anomaly Trough, in degrees per year, clockwise around inlet origin'};
% % 
% % Tc = [header; table2cell(troughTable)];
% % 
% % TablewithHeader = cell2table(Tc,'VariableNames',troughTable.Properties.VariableNames);
% % 
% % fName = [originName{oo} '_Rotation_Rate_Estimates_Troughs.txt'];
% % writetable(TablewithHeader,[outPath filesep fName])
% % 
% 
% 
% % %% Peak Finder
% % %From Pearson et al 2022 Supplementary Material 
% % %Ridges and troughs were identified in the Volume Anomaly timestack
% % %(Hovmoller diagram) using findpeaks. Peaks higher than 0.1x10^6 m^3 and
% % %further than 8 deg. apart at a given timestep were selected, discounting
% % %endpoints. Peaks were manually connected and a linear regression was applied to
% % %estimate the rotation migration rates per shoal and channel.
% % 
% % scaleFactor = 2;
% % minPeakDist = 8;
% % prom = 0.1;
% % numPeaks = 4;
% % 
% % peaks = struct();
% % peakNums = 1:numPeaks;
% % 
% % figure('windowstate','maximized');
% % clf
% % %two dimensional plot
% % % for ii=1:size(vAnomalyTheta,2)
% % %     hold on
% % %     if (thetaDegNorth(2)-thetaDegNorth(1)) < 0 %x must be increasing
% % %         x = fliplr(thetaDegNorth(1:end-1)); y = flipud(vAnomalyTheta(:,ii)).*scaleFactor+T(ii);
% % %     else 
% % %         x = thetaDegNorth(1:end-1); y = vAnomalyTheta(:,ii).*scaleFactor+T(ii);
% % %     end
% % %     [peaks(ii).pks,peaks(ii).locs,peaks(ii).width,peaks(ii).prom] = findpeaks(y,x,...
% % %         'MinPeakDistance',minPeakDist,'MinPeakProminence',prom);
% % %     pl(ii) = plot(x,y,'LineWidth',1.5,'DisplayName',num2str(T(ii)));
% % %     pt(ii) = scatter(peaks(ii).locs,peaks(ii).pks,'kv','filled');
% % %     pt(ii).Annotation.LegendInformation.IconDisplayStyle = 'off';
% % %     text(peaks(ii).locs+.02,peaks(ii).pks+T(ii),num2str((1:numel(peaks(ii).pks))'))
% % % end
% % 
% % %three dimensional plot
% % for ii=1:size(vAnomalyTheta,2)
% %     hold on
% %     if (thetaDegNorth(2)-thetaDegNorth(1)) < 0 %x must be increasing
% %         x = fliplr(thetaDegNorth(1:end-1)); z = flipud(vAnomalyTheta(:,ii)).*scaleFactor+T(ii);
% %     else 
% %         x = thetaDegNorth(1:end-1); z = vAnomalyTheta(:,ii).*scaleFactor+T(ii);
% %     end
% %     y = repmat(T(ii),size(x));
% %     pl(ii) = plot3(x,y,z,'LineWidth',1.5,'DisplayName',num2str(T(ii)));
% %     [peaks(ii).pks,peaks(ii).locs,peaks(ii).width,peaks(ii).prom] = findpeaks(z,x,...
% %         'MinPeakDistance',minPeakDist,'MinPeakProminence',prom,'SortStr','descend','NPeaks',numPeaks);
% %     pt(ii) = scatter3(peaks(ii).locs,repmat(T(ii),size(peaks(ii).locs)),peaks(ii).pks,'kv','filled');
% %     pt(ii).Annotation.LegendInformation.IconDisplayStyle = 'off';
% %     %text(peaks(ii).locs+.02,repmat(T(ii),size(peaks(ii).locs)),peaks(ii).pks+T(ii),num2str((1:numel(peaks(ii).pks))'))
% %     %number peaks in order 
% %     %TO DO fix labeling of peaks
% %     % [locsSorted(:,ii),locsOrder(:,ii)] = sort(peaks(ii).locs);
% %     % peaks(ii).peakNum = peakNums(locsOrder(:,ii));
% %     % for jj = 1:numPeaks
% %     %     text(peaks(ii).locs+.02,repmat(T(ii),size(peaks(ii).locs)),peaks(ii).pks+T(ii),num2str(peaks(ii).peakNum(jj)));
% %     % end
% % end
% % title('Changes in Volume Anomaly [\theta]')
% % xlabel('\theta [\circ]'); ylabel('Time [y]'); zlabel('Vol. Anomaly (m^3)');
% % set(gca,'Layer','top','GridColor','k','GridAlpha',0.4);
% % set(gca,'FontName','Myriad Pro','FontSize',fontSize)
% % grid on
% % legend('location','northeast')
% % view(-15,60)
% % 
% % 
% % pngName = [originName{oo} '_Peak_Finding'];
% % savefig([outPath filesep pngName '.fig']);
% % exportgraphics(gcf,[outPath filesep pngName '.png'],"Resolution",300);
% % 
% % %% Trough Finding
% % 
% % % numTroughs = numPeaks - 1;
% % % troughs = struct();
% % % 
% % % figure('windowstate','maximized');
% % % clf
% % % %two dimensional plot
% % % % for ii=1:size(vAnomalyTheta,2)
% % % %     hold on
% % % %     if (thetaDegNorth(2)-thetaDegNorth(1)) < 0 %x must be increasing
% % % %         x = fliplr(thetaDegNorth(1:end-1)); y = flipud(vAnomalyTheta(:,ii)).*scaleFactor+T(ii);
% % % %     else 
% % % %         x = thetaDegNorth(1:end-1); y = vAnomalyTheta(:,ii).*scaleFactor+T(ii);
% % % %     end
% % % %     [peaks(ii).pks,peaks(ii).locs,peaks(ii).width,peaks(ii).prom] = findpeaks(y,x,...
% % % %         'MinPeakDistance',minPeakDist,'MinPeakProminence',prom);
% % % %     pl(ii) = plot(x,y,'LineWidth',1.5,'DisplayName',num2str(T(ii)));
% % % %     pt(ii) = scatter(peaks(ii).locs,peaks(ii).pks,'kv','filled');
% % % %     pt(ii).Annotation.LegendInformation.IconDisplayStyle = 'off';
% % % %     text(peaks(ii).locs+.02,peaks(ii).pks+T(ii),num2str((1:numel(peaks(ii).pks))'))
% % % % end
% % % 
% % % %three dimensional plot
% % % for ii=1:size(vAnomalyTheta,2)
% % %     hold on
% % %     if (thetaDegNorth(2)-thetaDegNorth(1)) < 0 %x must be increasing
% % %         x = fliplr(thetaDegNorth(1:end-1)); z = flipud(-vAnomalyTheta(:,ii)).*scaleFactor+T(ii);
% % %     else 
% % %         x = thetaDegNorth(1:end-1); z = -vAnomalyTheta(:,ii).*scaleFactor+T(ii);
% % %     end
% % %     y = repmat(T(ii),size(x));
% % %     pl(ii) = plot3(x,y,z,'LineWidth',1.5,'DisplayName',num2str(T(ii)));
% % %     [troughs(ii).pks,troughs(ii).locs,troughs(ii).width,troughs(ii).prom] = findpeaks(z,x,...
% % %         'MinPeakDistance',minPeakDist,'MinPeakProminence',prom,'SortStr','descend','NPeaks',numTroughs);
% % %     pt(ii) = scatter3(troughs(ii).locs,repmat(T(ii),size(troughs(ii).locs)),troughs(ii).pks,'kv','filled');
% % %     pt(ii).Annotation.LegendInformation.IconDisplayStyle = 'off';
% % %     text(troughs(ii).locs+.02,repmat(T(ii),size(troughs(ii).locs)),troughs(ii).pks+T(ii),num2str((1:numel(troughs(ii).pks))'))
% % % end
% % % title('Changes in Volume Anomaly [\theta]')
% % % xlabel('\theta [\circ]'); ylabel('Time [y]'); zlabel('Vol. Anomaly (m^3)');
% % % set(gca,'Layer','top','GridColor','k','GridAlpha',0.4);
% % % set(gca,'FontName','Myriad Pro','FontSize',fontSize)
% % % set(gca,'ZDir','reverse'); %found troughs by flipping and finding 'peaks', flip back
% % % %TO DO fix z axis upside down
% % % grid on
% % % legend('location','northeast')
% % % view(-15,60)
% % % 
% % % pngName = [originName{oo} '_Trough_Finding'];
% % % savefig([outPath filesep pngName '.fig']);
% % % exportgraphics(gcf,[outPath filesep pngName '.png'],"Resolution",300);
% % 
% % 
% % 
% % 
% % 
% % %% Peak Rotation Estimates
% % 
% % %Output file of peak rotation estimates
% % %Slopes of fit lines from polyfit
% % % peakFitLineSlope = NaN*ones(1,length(peaks));
% % % for ii = 1:length(peaks)
% % %     peakFitLineSlope(ii) = round(peaks(ii)(1),2); %round for export
% % % end
% % % 
% % % 
% % % 
% % % peakNums = (1:MaxNumExtrema)';
% % % 
% % % peakTable = table(peakNums,peakFitLineSlope');
% % % peakTable.Properties.VariableNames = {'Peak Number','Slope of Fit Line'};
% % % 
% % % header = {'Peak Number, Clockwise around inlet origin',...
% % %     'Slope of Linear Fit Line for each Volume Anomaly Peak, in degrees per year, clockwise around inlet origin'};
% % % 
% % % Tc = [header; table2cell(peakTable)];
% % % 
% % % TablewithHeader = cell2table(Tc,'VariableNames',peakTable.Properties.VariableNames);
% % % 
% % % fName = [originName{oo} '_Rotation_Rate_Estimates_Peaks.txt'];
% % % writetable(TablewithHeader,[outPath filesep fName])
% % % 
% % % %TO DO add the R^2 to table
% % % 
% % % 
% % % %TO DO plot on t plus one vanom theta plot and see if that screws things up
% % % 
% % % 
% % 
% % %% Trough Rotation Estimates
% % 
% % % troughFitLineSlope = NaN*ones(1,length(pTroughs));
% % % for ii = 1:length(pTroughs)
% % %     troughFitLineSlope(ii) = round(pTroughs{ii}(1),2); %round for export
% % % end
% % % 
% % % 
% % % troughNums = (1:MaxNumExtrema)';
% % % 
% % % troughTable = table(troughNums,troughFitLineSlope');
% % % troughTable.Properties.VariableNames = {'Peak Number','Slope of Fit Line'};
% % % 
% % % header = {'Trough Number, Clockwise around inlet origin',...
% % %     'Slope of Linear Fit Line for each Volume Anomaly Trough, in degrees per year, clockwise around inlet origin'};
% % % 
% % % Tc = [header; table2cell(troughTable)];
% % % 
% % % TablewithHeader = cell2table(Tc,'VariableNames',troughTable.Properties.VariableNames);
% % % 
% % % fName = [originName{oo} '_Rotation_Rate_Estimates_Troughs.txt'];
% % % writetable(TablewithHeader,[outPath filesep fName])
% 
% 
% %% Cross Correlation 
% %From Pearson et al 2022 Supplementary Material 
% %Migration speed of morphological features along the coast can be estimated
% %using the peak cross-correlation lag between successive bathymetric
% %transects. 
% %Quantify rotation of entire ETD around inlet by calculating the normalized
% %cross correlation for each theta volume anomaly timeslice with its
% %predecessor with lags ranging from -25 deg. to 25 deg
% %Dividing the lag at the peak of the cross-correlation delta_theta_peak by the timestep
% %length dt (gaps between surveys were irregular) yields an estimate of the delta's
% %annual rotation rate omega_ETD_t.
% %Yields estimate of entire delta's rotation rate 
% 
% timegap = 1; 
% iCutoff = 1; %test no cutoff for now, TO DO how to add to any site
% testLags = 25; %tests from -25 to 25 degree lags
% %iCutoff = 159; %bin at which to cut off calculations, avoid updrift
% %shoreline dynamics
% 
% 
% peakXC = NaN*ones(length(T)-1,1);
% peakLag = NaN*ones(length(T)-1,1);
% migRate = NaN*ones(length(T)-1,1);
% ETDrotationRate = NaN*ones(length(T)-1,1);
% xc = NaN*ones(length(T),testLags*2+1);
% lags = NaN*ones(length(T),testLags*2+1);
% 
% for tt = 2:length(T) %tt = 3:length(T)
%     figure('windowstate','maximized')
%     clf
%     subplot(2,1,1)
%     hold on; box on; grid on;
%     plot(thetaDegNorth(iCutoff:end-1),vAnomalyTheta(iCutoff:end,tt-timegap),'color','k')
%     plot(thetaDegNorth(iCutoff:end-1),vAnomalyTheta(iCutoff:end,tt),'color','b')
%     title(['(a) Volume Anomaly (' num2str(T(tt-timegap)) ' to ' num2str(T(tt)) ') [\theta]'])
%     xlabel('\theta [\circ]'); ylabel('V_{anomaly} [m^3]');
%     set(gca,'Layer','top','GridColor','k','GridAlpha',0.4);
%     set(gca,'FontName','Myriad Pro','FontSize',fontSize)
%     legend(num2str(T(tt-timegap)), num2str(T(tt)),'location','best');
% 
%     subplot(2,1,2)
%     hold on; box on; grid on;
%     [xc(tt,:),lags(tt,:)] = xcorr(vAnomalyTheta(iCutoff:end,tt),...
%         vAnomalyTheta(iCutoff:end,tt-timegap),testLags,'normalized');
%     %[R(crosscorr),lags] = xcorr(x,y) 
%     stem(lags(tt,:),xc(tt,:));
%     peakXC(tt) = max(xc(tt,:)); %peak cross corr through time
%     peakLag(tt) = lags(tt,find(xc(tt,:)==max(xc(tt,:)),1)); %peak lag through time
%     migRate(tt) = peakLag(tt)./(T(tt)-T(tt-timegap)); %time weighted peak lag deg/yr
%     %omega ETD
%     ETDrotationRate(tt) = peakLag(tt)./(T(tt)-T(tt-timegap)); 
%     title(['(b) Cross Correlation (' num2str(T(tt-timegap)) ' to ' num2str(T(tt)) ')'])
%     xlabel('Lags \theta [\circ]'); ylabel('R_{norm} [-]');
%     ylim([-1 1]);
%     text(-19,-0.8,['Peak Cross Correlation: ' num2str(peakLag(tt)) '^o'],...
%         'FontName','Myriad Pro','FontSize',10);
%     text(-19,-0.6,['Rotation Rate: ' num2str(migRate(tt),'%.1f') '^o/year'],...
%         'FontName','Myriad Pro','FontSize',10);
%     set(gca,'Layer','top','GridColor','k','GridAlpha',0.4);
%     set(gca,'FontName','Myriad Pro','FontSize',fontSize)
% 
%     % Export figure to png file
%     pngName = [originName{oo} '_HovmullerPolarTheta_CrossCorrelation_' ...
%         num2str(T(tt-timegap)) '-' num2str(T(tt))];
%     savefig([outPath filesep pngName '.fig']);
%     exportgraphics(gcf,[outPath filesep pngName '.png'],"Resolution",300);
% end
% 
% 
% 
% %% Hovmuller diagram of cross-correlation
%     
%     meanETDRotationRate = mean(ETDrotationRate(2:end)); %This is the mean migration rate of the whole ETD 
%     %Positive is counterclockwise, negative is clockwise
%     
%     figure('windowstate','maximized');
%     clf
%     subplot(2,1,1)
%     hold on; box on; grid on;
%     pp=pcolor(repmat(lags(end,:)',1,length(T)),repmat(T,length(lags),1),xc');
%     pp.EdgeColor='none';
%     xline(0,'-k','LineWidth',2.5)
%     plot(peakLag(2:end),T(2:end),'-k')
%     cb=colorbar('Location','eastoutside');
%     title('Cross-Correlation Through Time')
%     xlabel('Lag [\circ]');
%     ylabel('Time [y]');
%     ylim([T(2) T(end)]);
%     xlim([-1*testLags testLags]);
%     clim([-1 1])
%     set(gca,'XTick',[-1*testLags:5:testLags])
%     set(gca,'Layer','top','GridColor','k','GridAlpha',0.4);
%     ylabel(cb,'R_{norm}','FontName','Myriad Pro','FontSize',fontSize,'fontweight','bold')
%     set(gca,'FontName','Myriad Pro','FontSize',fontSize)
%     
%     subplot(2,1,2)
%     hold on; box on; grid on;
%     plot(T(2:end),migRate(2:end),'*-k','LineWidth',2)
%     yline(meanETDRotationRate,':r','LineWidth',2.5)
%     title('Delta Rotation Rate \omega_{ETD,t}')
%     xlabel('Time [y]');
%     ylabel('d\theta/dt [\circ/yr]');
%     set(gca,'Layer','top','GridColor','k','GridAlpha',0.4);
%     set(gca,'FontName','Myriad Pro','FontSize',fontSize)
%     legend('ETD Rotation Rate',['Mean Rotation Rate = ' num2str(meanETDRotationRate)]);
% 
%     
%     % Export figure to png file
%     pngName = [originName{oo} '_HovmullerPolarTheta_XCorr_ETD_Rotation_Rate'];
%     savefig([outPath filesep pngName '.fig']);
%     exportgraphics(gcf,[outPath filesep pngName '.png'],"Resolution",300);
% 
%     %Output ETD Rotation Rates
%     OriginName = originName{oo};
%     OriginX = origins(oo,1);
%     OriginY = origins(oo,2);
% 
%     save(fullfile(outPath,[originName{oo} '_Delta_Rotation_Rates.mat']),...
%         'ETDrotationRate','meanETDRotationRate','OriginName','OriginX',...
%         'OriginY')
% 
% 
% 
% close all
% end %testMulitOrigins loop
% 
% 
% %% Sensitivity Testing
% %Compare delta rotation rate for different origins 
% %"A well chosen origin should yield a relatively higher delta rotation rate,
% %since motion around it will be more coherent than for an arbitrary origin"
% %Pearson et al., 2022
% 
% if testMultiOrigins
%   
%     %Load origin ETD rotation rates
%      
%     rotation_mat_files = dir(fullfile(outPathRoot,'**','*Rotation_Rates.mat'));
%     %load files
%     rotationRateFiles = cell(length(rotation_mat_files),1);
%     for ii = 1:length(rotation_mat_files)
%         rotationRateFiles{ii} = load(fullfile(rotation_mat_files(ii).folder,...
%             rotation_mat_files(ii).name));
%     end
%     
%     %load variables
%     OriginNames = cell(length(rotationRateFiles),1);
%     meanETDRotationRates = NaN*ones(length(rotationRateFiles),1);
%     for ii = 1:length(rotationRateFiles)
%         OriginNames{ii} = rotationRateFiles{ii}.OriginName;
%         meanETDRotationRates(ii) = rotationRateFiles{ii}.meanETDRotationRate;
%         OriginX(ii) = rotationRateFiles{ii}.OriginX;
%         OriginY(ii) = rotationRateFiles{ii}.OriginY;
%     end
%     
%     %make table
%     originTbl = table(OriginNames,round(meanETDRotationRates,2),OriginX',OriginY');
%     originTbl.Properties.VariableNames = {'OriginName','MeanETDRotationRate',...
%         'originX','originY'};
%     
%     %find highest rotation rate   
%     [~,maxLoc] = max(abs(originTbl.MeanETDRotationRate));
%    
%     %display best origin and rate
%     disp("The highest Origin Rotation Rate is " + num2str(originTbl.MeanETDRotationRate(maxLoc)) + " deg./yr for origin " + originTbl.OriginName(maxLoc));
% 
%     %output origin names, locations, rotation rates
%     fName = 'Origin_Table.txt';
%     writetable(originTbl,[outPathRoot filesep fName])
% 
%     %Sort origin Table by rotation rate
%     originTblSorted = sortrows(originTbl,2,'descend','ComparisonMethod','abs');
% 
%     %Plot origin rotation rates
%     figure('windowstate','maximized');
%     clf
%     hold on
%     imagescn(X,Y,Z(:,:,1)); %plot example bathy surface
%     colormap('gray');
%     scatterSize = linspace(20,500,height(originTblSorted)); %size by rotation rate
%     scatterSize = sort(scatterSize,'descend');
%     for i = 1:height(originTblSorted)
%         if i == 1
%             scatter(originTblSorted.originX(i),originTblSorted.originY(i),scatterSize(i),...
%                 'filled','MarkerEdgeColor','r','LineWidth',2); 
%         else 
%             scatter(originTblSorted.originX(i),originTblSorted.originY(i),scatterSize(i),...
%                 'filled');
%         end
%     end
%     colororder('gem12');
%     title('Origin Test Locations with Rotation Rates'); 
%     legend(originTblSorted.OriginName); grid on; xlabel('X (km)'); ylabel('Y (km)');
%     %Export figure to png file
%     pngName = 'origin_test_locations_with_rotation_rates';
%     savefig([outPathRoot filesep pngName '.fig']);
%     exportgraphics(gcf,[outPathRoot filesep pngName '.png'],"Resolution",300);
% 
% end %sensitivity testing loop
% 
% 
% 
% 
% 
% end %end function


