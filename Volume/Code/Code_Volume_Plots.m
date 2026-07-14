%% ================= Paths & params =================
csv_folder    = "C:\Users\mdeassisbose\OneDrive - University of Florida\Documents\ArcGIS\Projects\Chapter_2_Shoals_Beach\sediment_budged\New_Pass\Beck_Arnold_Volume\Paper_poly_1order_20251108";
output_folder = "C:\Users\mdeassisbose\OneDrive - University of Florida\Documents\ArcGIS\Projects\Chapter_2_Shoals_Beach\sediment_budged\New_Pass\Beck_Arnold_Volume\Paper_poly_1order_20251108\results\backup";
if ~exist(output_folder, 'dir'); mkdir(output_folder); end

threshold = 60;   % cobertura mínima (%)
max_year  = 2024; % limite do eixo X

show_fig = true;   % mostrar na tela?
save_png = true;   % salvar PNG?

%%================ Choose diff pairs ================
diff_pairs = [2004 2004;
              2004 2006;
              2004 2010;
              2004 2015;
              2004 2022];

%%================= Group dictionaries ==============
polygon_groups = struct( ...
    'Channel',          ["5","6"], ...
    'ETD_Lobe',         ["8"], ...
    'Shoal_Attachment', ["7","9"], ...
    'Ebb_Tidal_Delta',  ["2","3","4","5","6","7"], ...
    'Dredging_Peat',    ["11"] ...
);
selected_groups = "Ebb_Tidal_Delta";   % [] → all groups
selected_ids    = string.empty;        % [] → no filter

%================= 1) Read files =================
frames = table();
for k = 1:size(diff_pairs,1)
    base_year   = diff_pairs(k,1);
    target_year = diff_pairs(k,2);
    if target_year > max_year, continue; end
    %stem = sprintf("volume_Difference_2order_topbath_%d_vs_%d", target_year, base_year);
    stem = sprintf("volume_Difference_topbath_%d_vs_%d", target_year, base_year);
    f_csv = csv_folder + "\" + stem + ".csv";
    f_xls = csv_folder + "\" + stem + ".xlsx";

    if isfile(f_csv)
        T = readtable(f_csv, 'TextType','string', 'VariableNamingRule','preserve');
    elseif isfile(f_xls)
        T = readtable(f_xls, 'TextType','string', 'VariableNamingRule','preserve');
    else
        warning("File not found for (%d,%d): %s", base_year, target_year, stem);
        continue;
    end

    need = ["Polygon_ID","Volume_m3","Total_Area_m2","Coverage_%"];
    missing = need(~ismember(need, string(T.Properties.VariableNames)));
    if ~isempty(missing)
        error("%s is missing columns: %s", stem, strjoin(missing,", "));
    end

    T.Polygon_ID = string(T.Polygon_ID);
    T = T(T.("Coverage_%") >= threshold, :);
    T.Year = repmat(target_year, height(T), 1);

    T.Volume_per_Area_m = T.Volume_m3 ./ T.Total_Area_m2;  % m3/m2 = m
    frames = [frames; T(:, {'Polygon_ID','Year','Volume_m3','Total_Area_m2','Volume_per_Area_m'})]; %#ok<AGROW>
end
if isempty(frames); error("No data remaining after filtering."); end

% ================ 2) Pivot absolute volume (polygon × year) =================
pivot = unstack(frames, 'Volume_m3', 'Year', ...
                'AggregationFunction', @nansum, ...
                'GroupingVariables', 'Polygon_ID');

allVars  = string(pivot.Properties.VariableNames);
yearVars = allVars(allVars ~= "Polygon_ID");
yearsNum = str2double(regexprep(yearVars, "^x", ""));
[yearsNum, ord] = sort(yearsNum);
yearVars = yearVars(ord);
pivot    = pivot(:, ["Polygon_ID", yearVars]);


% ================ 2b) Pivot Volume / Area =================
pivot_area = unstack(frames, 'Volume_per_Area_m', 'Year', ...
                'AggregationFunction', @nanmean, ...
                'GroupingVariables', 'Polygon_ID');

allVars_area  = string(pivot_area.Properties.VariableNames);
yearVars_area = allVars_area(allVars_area ~= "Polygon_ID");
yearsNum_area = str2double(regexprep(yearVars_area, "^x", ""));
[yearsNum_area, ord_area] = sort(yearsNum_area);
yearVars_area = yearVars_area(ord_area);
pivot_area    = pivot_area(:, ["Polygon_ID", yearVars_area]);
% ======= 3) Cores fixas por Polygon_ID =======
color_fixed = containers.Map( ...
    {'2','3','4','5','6','7'}, ...
    { ...
      '#4682B4', ... % 2
      '#0B1E3D', ... % 3 - ETD North Platform
      '#1F4E8C', ... % 4 - Channel
      '#1E90FF', ... % 5 - Channel Margin Linear Bar
      '#B22222', ... % 6 - Lobe
      '#D37A22'  ... % 7 - ETD South Plataforma 1
    } ...
);

% Build the final color_map (ID -> [r g b]) using 'lines()' as fallback
all_ids    = cellstr(pivot.Polygon_ID);
unique_ids = unique(all_ids);
color_map  = containers.Map('KeyType','char','ValueType','any');

fallback = lines(max(numel(unique_ids), 7));
fidx = 1;
for ii = 1:numel(unique_ids)
    key = unique_ids{ii};
    if isKey(color_fixed, key)
        rgb = hex2rgb(color_fixed(key));
    else
        rgb = fallback(fidx,:); 
        fidx = fidx + 1; 
        if fidx > size(fallback,1), fidx = 1; end
    end
    color_map(key) = rgb;
end

% ======= Determine groups =======
keys = strings(0,1); vals = cell(0,1);
fnames = fieldnames(polygon_groups);
if isempty(selected_groups)
    for i = 1:numel(fnames)
        keys(end+1,1) = string(fnames{i});
        vals{end+1,1} = string(polygon_groups.(fnames{i}));
    end
else
    for i = 1:numel(selected_groups)
        g = string(selected_groups(i));
        if isfield(polygon_groups, g)
            keys(end+1,1) = g;  vals{end+1,1} = string(polygon_groups.(g));
        else
            warning("Group does not exist: %s", g);
        end
    end
end
if ~isempty(selected_ids)
    keys(end+1,1) = "__Explicit_IDs__";  vals{end+1,1} = string(selected_ids);
end
if isempty(keys)
    keys = "ALL";  vals = {pivot.Polygon_ID};
end

% ======= Fixed markers by Polygon_ID (with fallback) =======
marker_fixed = containers.Map( ...
    {'2','3','4','5','6','7'}, ...
    {'o','s','d','^','v','>'} ...
);
marker_cycle = {'o','s','d','^','v','>','<','p','h','x','+'}; % fallback

%% ======= Figure: Absolute Volume (default style) =======
for i = 1:numel(keys)
    gname       = keys(i);
    gname_clean = strrep(gname, "_", " ");
    id_list     = string(vals{i});
    ids_pres    = intersect(id_list, pivot.Polygon_ID, 'stable');
    if isempty(ids_pres)
        fprintf("Group '%s' contains no polygons. Skipping.\n", gname_clean);
        continue;
    end

    data = pivot(ismember(pivot.Polygon_ID, ids_pres), :);
    M    = table2array(data(:, yearVars));          % V (m³)
    aggV = nansum(M, 1);                            % Total por ano

    fig = figure('Name','ETD_Volume', 'Color','w', ...
        'Visible', tern(show_fig,'on','off'), ...
        'Units','centimeters','Position',[2 2 22 9]);
    set(fig,'PaperUnits','centimeters','PaperSize',[21 29.7], 'PaperPositionMode','auto');

    ax = axes(fig); hold(ax,'on'); box(ax,'on'); grid(ax,'on');
    set(ax,'FontName','Arial','FontSize',12,'LineWidth',1)

    % Individual curves
    for r = 1:size(M,1)
        id_poly = data.Polygon_ID(r);
        y = M(r,:)/1e6; if all(isnan(y)), continue; end
        col = color_map(char(id_poly));
        mk  = marker_for_id(id_poly, marker_fixed, marker_cycle);

        plot(ax, yearsNum, y, '-', 'LineWidth',1.6, ...
            'Color',col, 'Marker',mk, ...
            'MarkerEdgeColor',col, 'MarkerFaceColor',col, ...
            'MarkerSize',5.5, 'MarkerIndices',1:numel(yearsNum), ...
            'DisplayName', string(id_poly));
    end

    % Total
    if any(~isnan(aggV))
        plot(ax, yearsNum, aggV/1e6, ':', 'Color',[0 0 0], 'LineWidth',2, ...
             'DisplayName', "Total Volume");
    end

    yline(ax,0,'--','Color',[0 0 0],'HandleVisibility','off');
    xlim(ax,[min(yearsNum) max(yearsNum)]);
    ylabel(ax,'Volume (10^6 m^3)'); xlabel(ax,'Year');
    title(ax,'b) ETD  - Volume','FontSize',12,'FontWeight','normal');

    % Robust axis limits
    vecY = [M(:); aggV(:)]/1e6; vecY = vecY(~isnan(vecY));
    if ~isempty(vecY)
        yl = prctile(vecY,[2 98]); yl = yl + 0.12*[-range(yl) range(yl)];
        ylim(ax, yl);
    end

    if save_png && isgraphics(fig)
        exportgraphics(fig, fullfile(output_folder,'ETD_volume.png'), 'Resolution',300);
    end
end


%% ======= ΔV por intervalo (não acumulado), com 2004 = 0 no gráfico =======
base_year = 2004;
[tf_base, idx_base] = ismember(base_year, yearsNum);
if ~tf_base, error('Ano base %d não encontrado nas colunas.', base_year); end

for i = 1:numel(keys)
    gname       = keys(i);
    gname_clean = strrep(gname, "_", " ");
    id_list     = string(vals{i});
    ids_pres    = intersect(id_list, pivot.Polygon_ID, 'stable');
    if isempty(ids_pres), continue; end

    data = pivot(ismember(pivot.Polygon_ID, ids_pres), :);
    M    = table2array(data(:, string("x"+yearsNum)));   % m³

    % ΔV por intervalo e anos finais
    dt     = diff(yearsNum);                 %#ok<NASGU>
    x_end  = yearsNum(2:end);                % ano final de cada intervalo
    dV     = (M(:,2:end) - M(:,1:end-1));    % m³ per polygon (variação entre medições)
    agg_dV = nansum(dV, 1);                  % Total do grupo (por intervalo)

    % Para "começar" o gráfico em 2004 = 0 (apenas como anchor visual)
    x_plot    = [yearsNum(1), x_end];
    dV_plot   = [zeros(size(dV,1),1), dV]/1e6;     % 10^6 m³
    agg_plot  = [0, agg_dV]/1e6;                   % 10^6 m³

    % ===== figura =====
    fig = figure('Name','ETD_dV_interval','Color','w', ...
        'Visible', tern(show_fig,'on','off'), ...
        'Units','centimeters','Position',[2 2 22 9.5]);
    set(fig,'PaperUnits','centimeters','PaperSize',[21 29.7],'PaperPositionMode','auto');

    ax = axes(fig); hold(ax,'on'); box(ax,'on'); grid(ax,'on');
    set(ax,'FontName','Arial','FontSize',12,'LineWidth',1)

    for r = 1:size(dV_plot,1)
        id_poly = data.Polygon_ID(r);
        y = dV_plot(r,:); if all(isnan(y)), continue; end
        col = color_map(char(id_poly)); mk = marker_for_id(id_poly, marker_fixed, marker_cycle);
        plot(ax, x_plot, y, '-', 'LineWidth',1.6, 'Color',col, ...
            'Marker',mk,'MarkerEdgeColor',col,'MarkerFaceColor',col, ...
            'MarkerSize',5.5,'MarkerIndices',1:numel(x_plot), ...
            'DisplayName',string(id_poly));
    end

    plot(ax, x_plot, agg_plot, ':', 'Color',[0 0 0], 'LineWidth',2, ...
         'DisplayName','Total ΔV (intervalos)');

    yline(ax,0,'--','Color',[0 0 0],'HandleVisibility','off');
    xlim(ax,[min(yearsNum) max(yearsNum)]);
    xlabel(ax,'Year'); ylabel(ax,'ΔV (10^6 m^3)  [change per interval]');
    title(ax, sprintf('%s — ΔV por intervalo (anchor %d = 0)', gname_clean, base_year), ...
        'FontSize',12,'FontWeight','normal');

    vecY = [dV_plot(:); agg_plot(:)]; vecY = vecY(~isnan(vecY));
    if ~isempty(vecY)
        yl = prctile(vecY,[2 98]); yl = yl + 0.12*[-range(yl) range(yl)];
        if diff(yl)==0, yl = yl + [-1 1]*0.1; end
        ylim(ax, yl);
    end

    if save_png && isgraphics(fig)
        exportgraphics(fig, fullfile(output_folder,'DV_2004_ETD_volume.png'), 'Resolution',300);
    end
end

%% ======= RATE (ΔV/Δt) no ano final, com ponto 2004 = 0 [MÉDIA entre polígonos] =======
for i = 1:numel(keys)
    gname       = keys(i);
    gname_clean = strrep(gname, "_", " ");
    id_list     = string(vals{i});
    ids_pres    = intersect(id_list, pivot.Polygon_ID, 'stable');
    if isempty(ids_pres), continue; end

    % dados
    data = pivot(ismember(pivot.Polygon_ID, ids_pres), :);
    M    = table2array(data(:, string("x"+yearsNum)));   % m³

    % intervalos e rates per polygon
    dt     = diff(yearsNum);                 % anos (1 x nInt)
    x_end  = yearsNum(2:end);                % ano final de cada intervalo (1 x nInt)
    dV     = (M(:,2:end) - M(:,1:end-1));    % m³ (nPoly x nInt)
    rates  = dV ./ dt;                       % m³/ano (nPoly x nInt)

    % ===== média do rate entre polígonos (ignora NaN) =====
    meanRate = mean(rates, 1, 'omitnan');    % 1 x nInt, m³/ano

    % vetores para plot com anchor 2004 = 0
    x_plot        = [yearsNum(1), x_end];    % 2004, 2006, 2010, ...
    meanRate_plot = [0, meanRate]/1e6;       % 10^6 m³/ano, 2004 = 0

    % ===== figura =====
    fig = figure('Name','ETD_rate_mean_ref2004','Color','w', ...
        'Visible', tern(show_fig,'on','off'), ...
        'Units','centimeters','Position',[2 2 11 9]);
    set(fig,'PaperUnits','centimeters','PaperSize',[21 29.7],'PaperPositionMode','auto');

    ax = axes(fig); hold(ax,'on'); box(ax,'on'); grid(ax,'on');
    %set(ax,'FontName','Arial','FontSize',12,'LineWidth',1)

    % --- (opcional) desenhar também cada polígono, com 2004 = 0 ---
    for r = 1:size(rates,1)
        y_poly = [0, rates(r,:)]/1e6;                 % anchor 2004 = 0
        if all(isnan(y_poly(2:end))), continue; end
        id_poly = data.Polygon_ID(r);
        col = color_map(char(id_poly));
        mk  = marker_for_id(id_poly, marker_fixed, marker_cycle);

        plot(ax, x_plot, y_poly, '-', 'LineWidth', 1.6, ...
            'Color', col, ...
            'Marker', mk, ...
            'MarkerEdgeColor', col, ...
            'MarkerFaceColor', col, ...
            'MarkerSize', 5.0, ...
            'MarkerIndices', 1:numel(x_plot), ...     % <<< marca também 2004
            'HandleVisibility','off');                % não poluir a legenda
    end

    % --- curva da MÉDIA com marcador em TODOS os pontos (inclui 2004=0) ---
    h_mean = plot(ax, x_plot, meanRate_plot, ':', ...
        'Color',[0 0 0], 'LineWidth',2.0, ...
        'MarkerIndices', 1:numel(x_plot), ...         % <<< 2004 marcado
        'DisplayName', sprintf('Mean Rate (anchor %d = 0)', base_year));

    % zero, eixos, título
    yline(ax,0,'--','Color',[0 0 0],'HandleVisibility','off');
    xlim(ax,[min(yearsNum) max(yearsNum)]);
    xticks(min(yearsNum):2:max(yearsNum));
    xtickangle(45);
    %xlabel(ax,'Year');
    set(ax,'FontName','Arial','FontSize',12,'LineWidth',1)
    ylabel(ax,'Rate of Volume Change (10^6 m^3/yr)');
    title(ax, sprintf('c) ETD', gname_clean, base_year), ...
        'FontSize',12,'FontWeight','normal');

    % Robust axis limits (polígonos + média)
    vecY = [rates(:); meanRate(:)]/1e6;
    vecY = vecY(~isnan(vecY));
    if ~isempty(vecY)
        yl = prctile(vecY,[2 98]);
        yl = yl + 0.12*[-range(yl) range(yl)];
        if diff(yl)==0, yl = yl + [-1 1]*0.1; end
        ylim(ax, [-0.1 0.2] );
    end

    legend(ax, h_mean, {'Mean Rate'}, 'Location','northeast','Box','on','FontSize',12);

    if save_png && isgraphics(fig)
        exportgraphics(fig, fullfile(output_folder,'Rate_2004_ETD_volume.png'), 'Resolution',300);
    end
end


%%
figure('Color','w','Position',[200 200 700 350])
ax = axes; hold(ax,'on'); grid(ax,'on');

yyaxis left                       % ativa eixo Y esquerdo
plot(x_end, dV(1,:), '-o', 'LineWidth',1.8, 'Color',[0 0.3 0.8])
ylabel('\DeltaV (10^6 m^3)')      % rótulo Y1 em azul
set(ax,'YColor',[0 0.3 0.8])      % cor do eixo Y1

yyaxis right                      % ativa eixo Y direito
plot(x_end, rates(1,:), '--s', 'LineWidth',1.8, 'Color',[0.85 0.33 0.1])
ylabel('Rate (10^6 m^3/yr)')      % rótulo Y2 em laranja
set(ax,'YColor',[0.85 0.33 0.1])  % cor do eixo Y2

xlabel('Year')
title('Ebb-Tidal Delta — ΔV (left) & Rate (right)')
xticks(2004:2:2022)
xtickangle(90)
legend({'ΔV','Rate'},'Location','northwest')

%% ======= Figure: Volume Normalized by Area =======

% Área fixa per polygon
poly_area = groupsummary(frames, "Polygon_ID", "mean", "Total_Area_m2");
poly_area = poly_area(:, ["Polygon_ID", "mean_Total_Area_m2"]);
poly_area.Properties.VariableNames{'mean_Total_Area_m2'} = 'Area_m2';
poly_area = sortrows(poly_area, "Area_m2", "ascend");

for i = 1:numel(keys)

    gname       = keys(i);
    gname_clean = strrep(gname, "_", " ");
    id_list     = string(vals{i});

    ids_pres = intersect(id_list, pivot_area.Polygon_ID, 'stable');

    if isempty(ids_pres)
        fprintf("Group '%s' contains no polygons. Skipping.\n", gname_clean);
        continue;
    end

    data = pivot_area(ismember(pivot_area.Polygon_ID, ids_pres), :);

    % Ordenar os polígonos pela área: menor -> maior
    area_sub = poly_area(ismember(poly_area.Polygon_ID, data.Polygon_ID), :);
    area_sub = sortrows(area_sub, "Area_m2", "ascend");

    fig = figure('Name','ETD_Volume_per_Area', 'Color','w', ...
        'Visible', tern(show_fig,'on','off'), ...
        'Units','centimeters','Position',[2 2 22 9]);

    ax = axes(fig); hold(ax,'on'); box(ax,'on'); grid(ax,'on');
    set(ax,'FontName','Arial','FontSize',12,'LineWidth',1)

    h_plot = gobjects(height(area_sub),1);

    for rr = 1:height(area_sub)

        id_poly = area_sub.Polygon_ID(rr);

        row_idx = find(data.Polygon_ID == id_poly, 1);
        y = table2array(data(row_idx, yearVars_area));

        if all(isnan(y)), continue; end

        col = color_map(char(id_poly));
        mk  = marker_for_id(id_poly, marker_fixed, marker_cycle);

        A_m2 = area_sub.Area_m2(rr);
        A_km2 = A_m2 / 1e6;

        A_1e6 = A_m2 / 1e6;

        leg_name = sprintf('%s — %.3f ×10^6 m²', ...
            char(id_poly), A_1e6);
            

        h_plot(rr) = plot(ax, yearsNum_area, y, '-', ...
            'LineWidth',1.6, ...
            'Color',col, ...
            'Marker',mk, ...
            'MarkerEdgeColor',col, ...
            'MarkerFaceColor',col, ...
            'MarkerSize',5.5, ...
            'DisplayName',leg_name);
    end

    yline(ax,0,'--','Color',[0 0 0],'HandleVisibility','off');

    xlim(ax,[min(yearsNum_area) max(yearsNum_area)]);
    xticks(min(yearsNum_area):2:max(yearsNum_area));
    xtickangle(45);

    xlabel(ax,'Year');
    ylabel(ax,'Volume / Area (m)');

    % Legenda ordenada do menor para o maior polígono
    h_valid = h_plot(isgraphics(h_plot));
    lgd = legend(ax, h_valid, ...
    'Location','eastoutside', ...
    'Box','on', ...
    'FontSize',10);

    title(lgd,'Polygon Area');

    if save_png && isgraphics(fig)
        exportgraphics(fig, fullfile(output_folder,'ETD_volume_per_area.png'), ...
            'Resolution',300);
    end
end
%% ======= EXPORTAR TABELA ΔV & Rate =======
rows = {};
for i = 1:numel(keys)
    gname    = keys(i);
    id_list  = string(vals{i});
    ids_pres = intersect(id_list, pivot.Polygon_ID, 'stable');
    if isempty(ids_pres), continue; end

    data = pivot(ismember(pivot.Polygon_ID, ids_pres), :);
    M    = table2array(data(:, yearVars));

    dt    = diff(yearsNum);
    x_end = yearsNum(2:end);

    dV     = (M(:,2:end) - M(:,1:end-1));  % m³
    rates  = dV ./ dt;                     % m³/ano

    % per polygon
    for r = 1:size(M,1)
        pid = char(data.Polygon_ID(r));
        for j = 1:numel(x_end)
            rows(end+1,:) = {char(gname), pid, yearsNum(j), x_end(j), dV(r,j), rates(r,j)}; %#ok<AGROW>
        end
    end

    % Total do grupo
    aggV    = nansum(M,1);
    agg_dV  = (aggV(2:end) - aggV(1:end-1));
    aggRate = agg_dV ./ dt;
    for j = 1:numel(x_end)
        rows(end+1,:) = {char(gname), 'TOTAL', yearsNum(j), x_end(j), agg_dV(j), aggRate(j)}; %#ok<AGROW>
    end
end

T_out = cell2table(rows, 'VariableNames', ...
    {'Group','Polygon_ID','YearStart','YearEnd','DeltaV_m3','Rate_m3_per_yr'});

out_csv = fullfile(output_folder, '_DeltaV_Rate_byInterval.csv');
writetable(T_out, out_csv);
fprintf("File saved: %s\n", out_csv);


%%
%% ======= Figure: ETD Total Volume by Year =======
for i = 1:numel(keys)

    gname       = keys(i);
    gname_clean = strrep(gname, "_", " ");
    id_list     = string(vals{i});

    ids_pres = intersect(id_list, pivot.Polygon_ID, 'stable');
    if isempty(ids_pres)
        fprintf("Group '%s' contains no polygons. Skipping.\n", gname_clean);
        continue;
    end

    data = pivot(ismember(pivot.Polygon_ID, ids_pres), :);
    M    = table2array(data(:, yearVars));     % Volume dos polígonos
    aggV = nansum(M, 1);                       % Volume total ETD por ano

    fig = figure('Name','ETD_Total_Volume_Only', 'Color','w', ...
        'Visible', tern(show_fig,'on','off'), ...
        'Units','centimeters','Position',[2 2 22 6]);

    set(fig,'PaperUnits','centimeters','PaperSize',[21 29.7], ...
        'PaperPositionMode','auto');

    ax = axes(fig); hold(ax,'on'); box(ax,'on'); grid(ax,'on');
    set(ax,'FontName','Arial','FontSize',12,'LineWidth',1)

    plot(ax, yearsNum, aggV/1e6, '--o', ...
        'Color',[0 0 0], ...
        'LineWidth',2, ...
        'MarkerFaceColor',[0 0 0], ...
        'MarkerEdgeColor',[0 0 0], ...
        'MarkerSize',6);

    yline(ax,0,'--','Color',[0 0 0],'HandleVisibility','off');

    xlim(ax,[min(yearsNum) max(yearsNum)]);
    xticks(min(yearsNum):2:max(yearsNum));
    xtickangle(45);
    
    ylim(ax,[3.5 5]);
    xlabel(ax,'Year');
    ylabel(ax,'Volume (10^6 m^3)');
    title(ax,'b) ETD Total Volume','FontSize',12,'FontWeight','normal');

    if save_png && isgraphics(fig)
        exportgraphics(fig, fullfile(output_folder,'ETD_total_volume_only.png'), ...
            'Resolution',300);
    end
end

%% ---- função ternária ----
function o = tern(cond, a, b)
if cond, o = a; else, o = b; end
end

%% ---- util: HEX -> [r g b] ----
function rgb = hex2rgb(hex)
    hex = char(hex);
    if startsWith(hex,'#'), hex = hex(2:end); end
    if numel(hex) == 3, hex = reshape([hex;hex],1,[]); end  % 'abc' -> 'aabbcc'
    r = hex2dec(hex(1:2)); g = hex2dec(hex(3:4)); b = hex2dec(hex(5:6));
    rgb = [r g b] / 255;
end

function mk = marker_for_id(pid, marker_fixed, marker_cycle)
    pid = char(string(pid));
    if isKey(marker_fixed, pid)
        mk = marker_fixed(pid);
    else
        k  = max(1, mod(str2double(pid)-1, numel(marker_cycle))+1);
        mk = marker_cycle{k};
    end
end