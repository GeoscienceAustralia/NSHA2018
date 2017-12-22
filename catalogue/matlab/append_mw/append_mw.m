% *************************************************************************
% make_ML2MW_conv.m
% 
% Develop a conversion equation between ML to MW from Australian data for WA and EA
% 
% Program: T. Allen (2011-01-12)
% *************************************************************************

% load data
if exist('mdat','var') ~= 1
    disp('Loading mdat');
    load ..\mdat.mat;
end

%% pre-define polygons modified from shown in Allen 2010 AEES manuscript

ealon = [138.3 138.3 141.0 141.0 155.5 155.5 144.0 138.3];
ealat = [-10 -29 -29 -45 -45 -22 -10 -10];
walon = [138.3 127 111.5  111.5  135.0  135.0  138.3  138.3];
walat = [-11 -11 -18 -40   -40   -29   -29   -11];
salon = [135   135   141   141   135];
salat = [-29   -40   -40   -29   -29];

%% Get events with ML & ML in WA

% read WA events with MW
[Date,lon,lat,dep,mw,src] = textread('combined_au_mw.dat','%s%f%f%f%f%s', ...
                               'headerlines',1,'delimiter',',');
dateNum = datenum(Date);
unqdate = unique(dateNum);
clear meanMw unqLon unqLat unqDep unqDateNum;

% get average Mw
for i = 1:length(unqdate)
    ind = find(dateNum == unqdate(i));
    meanMw(i) = mean(mw(ind));
    unqLon(i) = lon(ind(1));
    unqLat(i) = lat(ind(1));
    unqDep(i) = dep(ind(1));
    unqDateNum(i) = dateNum(ind(1));
end

unqML = ones(size(meanMw)) * NaN;

% find ML's in mdat
t60 = 1 / (24 * 60);
for i = 1:length(mdat)
    if mdat(i).MDAT_dateNum > datenum(1968,10,13)
        ind = find(meanMw > mdat(i).GG_Mval - 1 & meanMw < mdat(i).GG_Mval + 1.2 ...
              & unqLon > mdat(i).MDAT_lon - 1 & unqLon < mdat(i).MDAT_lon + 1 ...
              & unqLat > mdat(i).MDAT_lat - 1 & unqLat < mdat(i).MDAT_lat + 1 ...
              & unqDateNum > mdat(i).MDAT_dateNum - t60 ...
              & unqDateNum < mdat(i).MDAT_dateNum + t60);
          if ~isempty(ind) & length(ind) == 1
              disp(['match: ',datestr(mdat(i).MDAT_dateNum, 31), ' ', num2str(i)])

              % fill mw value
              %if isnan(mdat(i).mdatMW)
                  mwind = find(dateNum == unqDateNum(ind));
                  mdat(i).altMW = mw(mwind(1));
                  tsrc = src(mwind(1));
                  mdat(i).altMWsrc = tsrc{1};
                  disp(mdat(i).altMWsrc)
              %end
          end
    end
end

save ..\mdat_mw mdat;

