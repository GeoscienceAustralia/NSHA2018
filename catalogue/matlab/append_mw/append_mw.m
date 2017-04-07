% *************************************************************************
% make_ML2MW_conv.m
% 
% Develop a conversion equation between ML to MW from Australian data for WA and EA
% 
% Program: T. Allen (2011-01-12)
% *************************************************************************

% load data
if exist('mdat_pref','var') ~= 1
    disp('Loading mdat_pref12');
    load ..\preferred\mdat_pref12.mat;
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

% find ML's in mdat_pref
t60 = 1.2 / (24 * 60);
for i = 1:length(mdat_pref)
    if mdat_pref(i).MDAT_dateNum > datenum(1968,10,13)
        ind = find(meanMw > mdat_pref(i).MDAT_prefML - 1 & meanMw < mdat_pref(i).MDAT_prefML + 1.2 ...
              & unqLon > mdat_pref(i).MDAT_lon - 1 & unqLon < mdat_pref(i).MDAT_lon + 1 ...
              & unqLat > mdat_pref(i).MDAT_lat - 1 & unqLat < mdat_pref(i).MDAT_lat + 1 ...
              & unqDateNum > mdat_pref(i).MDAT_dateNum - t60 ...
              & unqDateNum < mdat_pref(i).MDAT_dateNum + t60);
          if ~isempty(ind) & length(ind) == 1
              disp(['match: ',datestr(mdat_pref(i).MDAT_dateNum)])

              % fill mw value
              %if isnan(mdat_pref(i).MDAT_prefMW)
                  mwind = find(dateNum == unqDateNum(ind));
                  mdat_pref(i).MDAT_prefMW = mw(mwind(1));
                  tsrc = src(mwind(1));
                  mdat_pref(i).MDAT_prefMWSrc = tsrc{1};
                  disp(mdat_pref(i).MDAT_prefMWSrc)
              %end
          end
    end
end

save mdat_mw_pref12 mdat_pref;

