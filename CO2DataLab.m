%% Add your names in a comment here at the beginning of the code!
% Diana Hernandez & Jannitta Yao
% Instructions: Follow through this code step by step, while also referring
% to the overall instructions and questions from the lab assignment sheet.

%% 1. Read in the monthly gridded CO2 data from the .csv file
% The data file is included in your repository as “LDEO_GriddedCO2_month_flux_2006c.csv”
% Your task is to write code to read this in to MATLAB
% Hint: you can again use the function “readtable”, and use your first data lab code as an example.

CO2data = readtable('LDEO_GriddedCO2_month_flux_2006c.csv');

%% 2a. Create new 3-dimensional arrays to hold reshaped data
%Find each unique longitude, latitude, and month value that will define
%your 3-dimensional grid
longrid = unique(CO2data.LON); %finds all unique longitude values
latgrid = unique(CO2data.LAT); %<-- following the same approach, find all unique latitude values
monthgrid = unique(CO2data.MONTH);%<-- following the same approach, find all unique months

%Create empty 3-dimensional arrays of NaN values to hold your reshaped data
    %You can make these for any variables you want to extract - for this
    %lab you will need PCO2_SW (seawater pCO2) and SST (sea surface
    %temperature)
PCO2_SW = NaN([size(latgrid, 1) size(longrid, 1) size(monthgrid, 1)]);
SST = NaN([size(latgrid, 1) size(longrid, 1) size(monthgrid, 1)]);

%% 2b. Pull out the seawater pCO2 (PCO2_SW) and sea surface temperature (SST)
%data and reshape it into your new 3-dimensional arrays

for i = 2:length(CO2data.LON)
    
    tempLat = CO2data.LAT(i);
    tempLon = CO2data.LON(i);
    tempMonth = CO2data.MONTH(i);
    tempSST = CO2data.SST(i);
    tempPCO2_SW = CO2data.PCO2_SW(i);
    
    latID = find(tempLat == latgrid); 
    lonID = find(tempLon == longrid);
    monthID = find(tempMonth == monthgrid);
    
    PCO2_SW(latID, lonID, monthID) = tempPCO2_SW;
    SST(latID, lonID, monthID) = tempSST;
    
end
%% 3a. Make a quick plot to check that your reshaped data looks reasonable
%Use the imagesc plotting function, which will show a different color for
%each grid cell in your map. Since you can't plot all months at once, you
%will have to pick one at a time to check - i.e. this example is just for
%January

imagesc(flipud(SST(:,:,1)))

%% 3b. Now pretty global maps of one month of each of SST and pCO2 data.
%I have provided example code for plotting January sea surface temperature
%(though you may need to make modifications based on differences in how you
%set up or named your variables above).

figure(1); clf
worldmap world
contourfm(latgrid, longrid, SST(:,:,1),'linecolor','none');
colorbar
geoshow('landareas.shp','FaceColor','black')
title('January Sea Surface Temperature (^oC)')

%Check that you can make a similar type of global map for another month
%and/or for pCO2 using this approach. Check the documentation and see
%whether you can modify features of this map such as the contouring
%interval, color of the contour lines, labels, etc.
%June
figure(2); clf
worldmap world
contourfm(latgrid, longrid, PCO2_SW(:,:,6),'linecolor','none');
colorbar
geoshow('landareas.shp','FaceColor','black')
title('June Sea Surface Temperature (^oC)')

%% 4. Calculate and plot a global map of annual mean pCO2

PCO2mean = nanmean(PCO2_SW,3);

figure(1); clf
worldmap world
contourfm(latgrid, longrid, PCO2mean,'linecolor','none');
colorbar
geoshow('landareas.shp','FaceColor','black')
title('Annual Mean pCO2')

%% 5. Calculate and plot a global map of the difference between the annual mean seawater and atmosphere pCO2
%Reference year: 2000
%Link for Annual Mean of 2000: https://www.esrl.noaa.gov/gmd/ccgg/trends/gl_data.html?fbclid=IwAR0XJxcZWeFv4hX9TvzgXD-ZqZHZDl5HtmjGMkMzeAyTQoKwQzUns-fYR1Y
% Mean Atmospheric pCO2 for 2000: 368.84

At_pCO2mean = 368.84
difference = PCO2mean - At_pCO2mean

figure(1); clf
worldmap world
contourfm(latgrid, longrid, difference,'linecolor','none');
colorbar
geoshow('landareas.shp','FaceColor','black')
title('Difference between the Annual Mean Seawater pCO2 and Mean Atmospheric pCO2 in 2000')

%% 6. Calculate relative roles of temperature and of biology/physics in controlling seasonal cycle

tempMean = mean(SST,3);
SST_annual_mean = repmat(tempMean, 1, 1, 12);

PCO2_annual_mean = repmat(PCO2mean, 1, 1, 12);


PCO2_T = PCO2_SW .* exp(0.0423*(SST_annual_mean - SST));

PCO2_BP = PCO2_annual_mean .* exp(0.0423*(SST - SST_annual_mean)); 

%% 7. Pull out and plot the seasonal cycle data from stations of interest
%Do for BATS, Station P, and Ross Sea (note that Ross Sea is along a
%section of 14 degrees longitude - I picked the middle point)

%Lat = N or S
%Lon = E or W
% W or S means add a negative to the degrees and add 360

%BATS = 32 degrees 50 min N and 64 degrees 10 min W
%P = = 50 degrees N and 145 degrees W
%Ross = 76.5 degrees S and 169 degrees E

%-----------------BATS Station-------------------------
B_lat = 32.5;
B_lon = -64.1 + 360;

B_lat_diff = abs(latgrid - B_lat);
B_lat_index = find(B_lat_diff == min(B_lat_diff));
latgrid(B_lat_index);

B_lon_diff = abs(longrid - B_lon);
B_lon_index = find(B_lon_diff == min(B_lon_diff));
longrid(B_lon_index);

B_SST = squeeze(SST(B_lat_index, B_lon_index, :));
B_PCO2 = squeeze(PCO2_SW(B_lat_index, B_lon_index, :));
B_PCO2_T = squeeze(PCO2_T(B_lat_index, B_lon_index, :));
B_PCO2_BP = squeeze(PCO2_BP(B_lat_index, B_lon_index, :));


figure (1); clf
subplot(2,2,1)
plot(monthgrid, B_SST)
title('Sea Surface Temperature at BATS');
xlabel('Months')
ylabel('Temperature (C)')

subplot(2,2,2)
plot( monthgrid, B_PCO2)
title('PCO2 at BATS')
xlabel('Months')
ylabel('PCO2 (ppm)')

subplot(2,2,3)
plot( monthgrid, B_PCO2_T)
title('Temperature Effect on PCO2 at BATS')
xlabel('Months')
ylabel('PCO2 (ppm)')

subplot(2,2,4)
plot( monthgrid, B_PCO2_BP)
title('Bio-physical Effect on PCO2 at BATS')
xlabel('Months')
ylabel('PCO2 (ppm)')

%-----------------Papa Station------------------------
P_lat = 50;
P_lon = -145 + 360;

P_lat_diff = abs(latgrid - P_lat);
P_lat_index = find(P_lat_diff == min(P_lat_diff));
latgrid(P_lat_index);

P_lon_diff = abs(longrid - P_lon);
P_lon_index = find(P_lon_diff == min(P_lon_diff));
longrid(P_lon_index);

P_SST = squeeze(SST(P_lat_index(1), P_lon_index, :));
P_PCO2 = squeeze(PCO2_SW(P_lat_index(1), P_lon_index, :));
P_PCO2_T = squeeze(PCO2_T(P_lat_index(1), P_lon_index, :));
P_PCO2_BP = squeeze(PCO2_BP(P_lat_index(1), P_lon_index, :));


figure (2); clf
subplot(2,2,1)
plot(monthgrid, P_SST)
title('Sea Surface Temperature at Station Papa');
xlabel('Months')
ylabel('Temperature (C)')

subplot(2,2,2)
plot( monthgrid, P_PCO2)
title('PCO2 at Station Papa')
xlabel('Months')
ylabel('PCO2 (ppm)')

subplot(2,2,3)
plot( monthgrid, P_PCO2_T)
title('Temperature Effect on PCO2 at Station Papa')
xlabel('Months')
ylabel('PCO2 (ppm)')

subplot(2,2,4)
plot( monthgrid, P_PCO2_BP)
title('Bio-physical Effect on PCO2 at Station Papa')
xlabel('Months')
ylabel('PCO2 (ppm)')
%---------------------Ross--------------------------

R_lat = -76.5;
R_lon = ((169 - 177 + 360)/2);

R_lat_diff = abs(latgrid - R_lat);
R_lat_index = find(R_lat_diff == min(R_lat_diff));
latgrid(R_lat_index);

R_lon_diff = abs(longrid - R_lon);
R_lon_index = find(R_lon_diff == min(R_lon_diff));
longrid(R_lon_index);

R_SST = squeeze(SST(R_lat_index, R_lon_index, :));
R_PCO2 = squeeze(PCO2_SW(R_lat_index, R_lon_index, :));
R_PCO2_T = squeeze(PCO2_T(R_lat_index, R_lon_index, :));
R_PCO2_BP = squeeze(PCO2_BP(R_lat_index, R_lon_index, :));


figure (3); clf
subplot(2,2,1)
plot(monthgrid, R_SST)
title('Sea Surface Temperature at the Ross Sea');
xlabel('Months')
ylabel('Temperature (C)')

subplot(2,2,2)
plot( monthgrid, R_PCO2)
title('PCO2 at the Ross Sea')
xlabel('Months')
ylabel('PCO2 (ppm)')

subplot(2,2,3)
plot( monthgrid, R_PCO2_T)
title('Temperature Effect on PCO2 at the Ross Sea')
xlabel('Months')
ylabel('PCO2 (ppm)')

subplot(2,2,4)
plot( monthgrid, R_PCO2_BP)
title('Bio-physical Effect on PCO2 at the Ross Sea')
xlabel('Months')
ylabel('PCO2 (ppm)')

%% 8. Reproduce your own versions of the maps in figures 7-9 in Takahashi et al. 2002
% But please use better colormaps!!!
% Mark on thesese maps the locations of the three stations for which you plotted the
% seasonal cycle above

%<--
