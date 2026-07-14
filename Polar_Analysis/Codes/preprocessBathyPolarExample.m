%Preprocess topobathymetric lidar data for use in polarAnalysis code

% get lidar data
fpath = 'C:\Users\mdeassisbose\OneDrive - University of Florida\Documents\ArcGIS\Projects\Chapter_2_Shoals_Beach\conf_mapping\DEM_EVF\paper'; %path to folder with lidar data
%fpath = 'C:\Users\mdeassisbose\OneDrive - University of Florida\Documents\ArcGIS\Projects\OR_Siuslaw\conf_map_rasters\toptoraster';
%fpath = 'C:\Users\mdeassisbose\OneDrive - University of Florida\Documents\ArcGIS\Projects\MI_FrankFort\conf_mapping'
addpath(genpath(fpath)); 
fnames = dir([fpath '\**\*.tif']); %get all tif filenames

%put into structure
fnamest = struct2table(fnames);
%sort elevation by survey date
%New Pass
fnamest.t = [2004 2006 2010 2015 2022]'; % survey dates
% FrankFort
%fnamest.t = [2018 2019 2020 2021 2022 2023]';
fnames_sorted = sortrows(fnamest,5);

%pull X Y Z data for each survey
%these should already have grid cells snapped (ie grid cells should overlap
%already)
for i = 1:height(fnames_sorted)
    [Zfull{i},R{i}] = readgeoraster(fnames_sorted.name(i),"OutputType","double","CoordinateSystemType","planar");
    info{i} = georasterinfo(fnames_sorted.name(i));
    Zfull{i} = standardizeMissing(Zfull{i},info{i}.MissingDataIndicator);
    [Xfull{i},Yfull{i}] = worldGrid(R{i});
end

figure
clf
tlo = tiledlayout('flow');
tlo.TileSpacing = 'tight'; tlo.Padding = 'tight';
for h = 1:length(Zfull)
    ax(h) = nexttile;
    imagescn(Xfull{h},Yfull{h},Zfull{h});
    c = colorbar; %ylabel(c,clab);
    title(num2str(fnamest.t(h)));
    ax(h).DataAspectRatio = [1 1 1];
end
tlo.XLabel.String = 'X (m)'; tlo.YLabel.String = 'Y (m)';
colormap('parula')


%find all x and y bounds
for n = 1:length(Xfull)
    [minXfull(n),maxXfull(n)] = bounds(Xfull{n},"all");
    [minYfull(n),maxYfull(n)] = bounds(Yfull{n},"all");
end
minXfull = min(minXfull,[],'all'); maxXfull = max(maxXfull,[],'all');
minYfull = min(minYfull,[],'all'); maxYfull = max(maxYfull,[],'all');

%xv = linspace(minXfull,maxXfull,1000);
%yv = linspace(minYfull,maxYfull,1000);

%[X,Y] = meshgrid(xv,yv);
% Defina a resolução desejada (em metros)
dx = 3;
dy = 3;

% Corrija os limites para garantir que eles "fechem" no passo correto
minX = floor(min(minXfull) / dx) * dx;
maxX = ceil(max(maxXfull) / dx) * dx;

minY = floor(min(minYfull) / dy) * dy;
maxY = ceil(max(maxYfull) / dy) * dy;

% Crie os vetores com passo fixo de 3 metros
xv = minX : dx : maxX;
yv = minY : dy : maxY;

% Crie o grid regular
[X, Y] = meshgrid(xv, yv);
for m = 1:length(Zfull)
    Z{m} = griddata(Xfull{m},Yfull{m},Zfull{m},X,Y);
end


%save as mat file for conformal code

grd = struct();
grd.x = X; 
grd.y = Y;
for i = 1:numel(Z) 
    grd.dp(:,:,i) = Z{i}; 
end
grd.year = fnames_sorted.t;

matName = 'Test_NEWPASS_2004_2022.mat';
outPathMat  = fullfile(pwd,'Paper_NEWPASS_2004_2022.mat');
save(outPathMat,'grd')

%plot grd, make sure it all looks right
% figure
% clf
% tlo = tiledlayout('flow');
% tlo.TileSpacing = 'tight'; tlo.Padding = 'tight';
% for h = 1:size(grd.dp,3)
%     ax(h) = nexttile;
%     imagescn(grd.x,grd.y,grd.dp(:,:,h));
%     c = colorbar; ylabel(c,'Depth (m)');
%     title(num2str(grd.year(h)));
%     ax(h).DataAspectRatio = [1 1 1];
%     grid on
% end
% tlo.XLabel.String = 'X (m)'; tlo.YLabel.String = 'Y (m)';
% colormap('parula');
%% modificado pr Matheus 07/17/2025

%% modificado por Matheus 07/17/2025 - versão sem MAP_Toolbox

% Preprocess topobathymetric lidar data for use in polarAnalysis code

% % get lidar data
% fpath = 'C:\Users\mdeassisbose\OneDrive - University of Florida\Documents\ArcGIS\Projects\Chapter_2_Shoals_Beach\conf_mapping\DEM_EVF'; 
% addpath(genpath(fpath)); 
% fnames = dir([fpath '\**\*.tif']); % get all tif filenames
% 
% % put into structure
% fnamest = struct2table(fnames);
% % sort elevation by survey date
% fnamest.t = [2004 2006 2010 2015 2017 2022 2024]'; % survey dates
% fnames_sorted = sortrows(fnamest, 7);

% pull X Y Z data for each survey
for i = 1:height(fnames_sorted)
    % === Ler imagem ===
    Zraw = imread(fnames_sorted.name{i});
    Zfull{i} = double(Zraw);

    % === Ler metadados sem usar Mapping Toolbox ===
    infoRaw = imfinfo(fnames_sorted.name{i});
    
    % Extrair escala e origem dos tags TIFF
    scale = infoRaw.ModelPixelScaleTag;         % [dx, dy, dz]
    tiepoint = infoRaw.ModelTiepointTag;        % [i, j, k, x0, y0, z0]

    dx = scale(1);  % resolução em X
    dy = scale(2);  % resolução em Y
    x0 = tiepoint(4);  % origem X
    y0 = tiepoint(5);  % origem Y

    % === Criar grade X/Y ===
    [nrows, ncols] = size(Zfull{i});
    x = x0 + (0:ncols-1) * dx;
    y = y0 - (0:nrows-1) * dy;  % Y invertido porque origem é canto superior esquerdo
    [Xfull{i}, Yfull{i}] = meshgrid(x, y);

    % Substituir valores inválidos se necessário (ex: -9999)
    Zfull{i}(Zfull{i} <= -9999) = NaN;
end

% === Visualização rápida ===
figure
clf
tlo = tiledlayout('flow');
tlo.TileSpacing = 'tight'; tlo.Padding = 'tight';
for h = 1:length(Zfull)
    ax(h) = nexttile;
    imagesc(Xfull{h}(1,:), Yfull{h}(:,1), Zfull{h});
    set(gca, 'YDir', 'normal');
    c = colorbar;
    title(num2str(fnamest.t(h)));
    ax(h).DataAspectRatio = [1 1 1];
end
tlo.XLabel.String = 'X (m)'; 
tlo.YLabel.String = 'Y (m)';
colormap('parula');

% === Encontrar limites espaciais ===
for n = 1:length(Xfull)
    [minXfull(n), maxXfull(n)] = bounds(Xfull{n}, "all");
    [minYfull(n), maxYfull(n)] = bounds(Yfull{n}, "all");
end
minXfull = min(minXfull); maxXfull = max(maxXfull);
minYfull = min(minYfull); maxYfull = max(maxYfull);

% === Criar grade comum para interpolação ===
xv = linspace(minXfull, maxXfull, 1000);
yv = linspace(minYfull, maxYfull, 1000);
[X, Y] = meshgrid(xv, yv);

% === Interpolar todos os rasters para a mesma grade ===
for m = 1:length(Zfull)
    Z{m} = griddata(Xfull{m}, Yfull{m}, Zfull{m}, X, Y);
end

% === Salvar como .mat para uso posterior ===
grd = struct();
grd.x = X;
grd.y = Y;
for i = 1:numel(Z)
    grd.dp(:,:,i) = Z{i};
end
grd.year = fnames_sorted.t;

outPathMat  = fullfile(pwd, 'Paper_NEWPASS_2004_2022.mat');
save(outPathMat, 'grd');
%% 3 opsao

% Caminho para os arquivos LiDAR
fpath = 'C:\Users\mdeassisbose\OneDrive - University of Florida\Documents\ArcGIS\Projects\Chapter_2_Shoals_Beach\conf_mapping\DEM_EVF';
addpath(genpath(fpath)); 

% Obter todos os arquivos .tif
fnames = dir(fullfile(fpath, '**', '*.tif'));

% Converter em tabela e adicionar datas de levantamento
fnamest = struct2table(fnames);
fnamest.t = [2004 2006 2010 2015 2017 2022 2024]'; % ajuste conforme necessário
fnames_sorted = sortrows(fnamest, 7);  % ordena por data (coluna 7 = 't')

% Inicializar e ler dados dos rasters
for i = 1:height(fnames_sorted)
    fname = fullfile(fnames_sorted.folder{i}, fnames_sorted.name{i});
    
    [Zfull{i}, R{i}] = readgeoraster(fname, ...
        "OutputType", "double", ...
        "CoordinateSystemType", "planar");

    % Substituir valores ausentes (-9999 ou outro, ajuste se necessário)
    Zfull{i}(Zfull{i} == -9999) = NaN;

    % Obter coordenadas X e Y associadas à grade
    [Xfull{i}, Yfull{i}] = worldGrid(R{i});
end

% Plotar os rasters
figure; clf
tlo = tiledlayout('flow');
tlo.TileSpacing = 'tight'; tlo.Padding = 'tight';

for h = 1:length(Zfull)
    ax(h) = nexttile;
    imagescn(Xfull{h}, Yfull{h}, Zfull{h});
    colorbar;
    title(num2str(fnames_sorted.t(h)));
    ax(h).DataAspectRatio = [1 1 1];
end
tlo.XLabel.String = 'X (m)';
tlo.YLabel.String = 'Y (m)';
colormap('parula');

% Determinar os limites espaciais globais
for n = 1:length(Xfull)
    [minXfull(n), maxXfull(n)] = bounds(Xfull{n}, "all");
    [minYfull(n), maxYfull(n)] = bounds(Yfull{n}, "all");
end

minXfull = min(minXfull); maxXfull = max(maxXfull);
minYfull = min(minYfull); maxYfull = max(maxYfull);

% Criar grade regular de interpolação
xv = linspace(minXfull, maxXfull, 1000);
yv = linspace(minYfull, maxYfull, 1000);
[X, Y] = meshgrid(xv, yv);

% Interpolar cada DEM para a nova grade
for m = 1:length(Zfull)
    Z{m} = griddata(Xfull{m}, Yfull{m}, Zfull{m}, X, Y);
end

% Salvar arquivo .mat para uso posterior
grd = struct();
grd.x = X;
grd.y = Y;
for i = 1:numel(Z)
    grd.dp(:,:,i) = Z{i};
end
grd.year = fnames_sorted.t;

outPathMat = fullfile(pwd, 'KT_New_Pass_2004_2024.mat');
save(outPathMat, 'grd');
