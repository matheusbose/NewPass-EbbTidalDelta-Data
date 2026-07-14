%Example polar analysis code run 
%INPUTS:
    %datafile: .mat file (ex: 'bathy_data.mat') containing a structure called 'grd' with fields: 
        %year: year for each data surface, in 'double' format. ex: [2010, 2015, 2020,
            %2023] %%CANNOT HAVE MULTIPLES OF SAME YEAR, TO DO fix this
        %x: gridded x data in meters, in projected coordinates, ie
            %Eastings, size m x n
        %y: gridded y data in meters, in projected coordinates, ie
            %Northings, size m x n
        %dp; gridded z data, in meters, depth negative, size m x n x year,
            %topo data should already be masked out
    %origin: x,y coordinates (meters) of origin for polar grid, usually center of
        %inlet thalweg, ex: [547014,3361720]
    %outPath: filepath to folder where results will be output
    %makeMiddlePlots: make or don't make plots of hAnomaly and vAnomaly for all
        %surveys. 0 or 1.
    %dRho: polar grid spacing in radial direction (m), ex: 40 
    %rhoLimOuter: polar grid outer limits (meters from origin). ex: rhoLim = 1500, makes a
        %grid from 0 (origin) to 1500 meters from origin, spaced by dRho.
        %Choose limits which cover the ebb shoal features of
        %interest. The inner 15% of the polar grid will be removed from analysis, as the
        %grid cells become to small to have data points in every grid cell.
    %thetaLims: polar grid sector limits, 180 degrees, in nautical
        %convention (90 = North, 0 = East), ex: [0 180] for an east facing
        %inlet
    %dTheta: polar grid angle spacing (degrees), ex: 1 or -1
        %TO DO: automatically correct for dTheta sign, for now, do positive
        %dTheta if in the trigonometric positive half of the unit circle
        %(0-180, or 90-270 in nautical convention) and do negative dTheta
        %if in the trigonometric negative half of the unit circle (-180-0,
        %or 270-90 in nautical convention)
    %testMultiOrigins: 0 or 1. If 1, test 12 other origins around the chosen origin. 
        %Compares ebb tidal delta rotation rate for different origins to chose ideal origin location. 
        %From Pearson et al., 2022 Supplementary Matieral: Well-chosen
        %origin should have a relatively higher delta rotation rate and
        %motion around it will be more coherent. 
        %Outputs delta rotation rates for each test origin and the origin site most
        %recommended for use in analysis.
        %WILL TAKE SIGNIFICANTLY LONGER

% Moriches Inlet NY

cd('C:\Users\mdeassisbose\OneDrive - University of Florida\Documents\MATLAB\Conformal_Mapping\conf_map_V02') %where the code is located 
addpath('C:\Users\mdeassisbose\OneDrive - University of Florida\Documents\MATLAB\Conformal_Mapping\conf_map_V02') %where the data is located (if different)

datafile = "Paper_NEWPASS_2004_2022";
%New Pass
origin = [342887.66, 3023740.95] ;

%
outPathRoot = 'C:\Users\mdeassisbose\OneDrive - University of Florida\Documents\MATLAB\Conformal_Mapping\conf_map_V02\Paper_NewPass\Paper_fig';
makeMiddlePlots = 1;
dTheta = 3;
rhoLimOuter = 1800;
dRho = 20;
%New Pass
thetaLims = [130 290];
%thetaLims = [275 455];
%Frankfort
% thetaLims = [120 260]; 
testMultiOrigins = 0;
interest = 0.15;  %The inner 15% of the polar grid will be removed from analysis, as the  %grid cells become to small to have data points in every grid cell.

% Original
%polarAnalysis(datafile,origin,outPathRoot,makeMiddlePlots,dTheta,dRho,rhoLimOuter,thetaLims,testMultiOrigins,interest)
%modificado
polarAnalysis(datafile,origin,outPathRoot,dTheta,dRho,rhoLimOuter,thetaLims,testMultiOrigins,interest)
