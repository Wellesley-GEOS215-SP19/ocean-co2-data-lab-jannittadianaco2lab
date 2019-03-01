%% Add your names in a comment here at the beginning of the code!
%% Diana Hernandez & Jannitta Yao
% Instructions: Follow through this code step by step, while also referring
% to the overall instructions and questions from the lab assignment sheet.

%% 1. Read in the monthly gridded CO2 data from the .csv file
% The data file is included in your repository as �LDEO_GriddedCO2_month_flux_2006c.csv�
% Your task is to write code to read this in to MATLAB
% Hint: you can again use the function �readtable�, and use your first data lab code as an example.

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
%<--

%% 7. Pull out and plot the seasonal cycle data from stations of interest
%Do for BATS, Station P, and Ross Sea (note that Ross Sea is along a
%section of 14 degrees longitude - I picked the middle point)

%<--

%% 8. Reproduce your own versions of the maps in figures 7-9 in Takahashi et al. 2002
% But please use better colormaps!!!
% Mark on thesese maps the locations of the three stations for which you plotted the
% seasonal cycle above

%<--
