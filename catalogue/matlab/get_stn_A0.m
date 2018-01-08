function [rhyp,repi,MLM92,WGW94,WGW96,GS86,GG91,R35,BJ84,stns] ...
         = get_stn_A0(siteDat,lat,lon,dep,mag,evdate,zone)

% function gets A0 correction factors for recording stations or those that
% may have recorded the event if no stations exist
% 
% zone = 1 > WA
% zone = 2 > EA
% zone = 3 > SA

if zone == 1 & isnan(dep)
    dep = 5;
elseif zone == 2 & isnan(dep)
    dep = 10;
elseif zone == 3 & isnan(dep)
    dep = 10;
end

%% set distance ranges to to account for offscale measurements

offScaleDist = 0;
if mag >= 4.0 && mag < 4.5
    offScaleDist = 75;
elseif mag >= 4.5 && mag < 5.0
    offScaleDist = 150;
elseif mag >= 5.0
    offScaleDist = 250;
end

% reset if modern digital recordings
if evdate > datenum(1990,1,1)
    offScaleDist = 0;
end   
 
%% make proxy stationlist if ISC stns do not exist  
%if isempty(stns) & isempty(essrng)
[rng az] = distance([siteDat.stnlat],[siteDat.stnlon],lat,lon);

% find stations within date & distnce range
ind = find([siteDat.startdate] <= evdate & [siteDat.stopdate] >= evdate ...
      & deg2km(rng') < 1500 & deg2km(rng') >= offScaleDist);
if ~isempty(ind)
    stns = {siteDat.stncode(ind)};
    repi = deg2km(rng(ind));
    rhyp = sqrt(repi.^2 + (ones(size(repi))*dep).^2);
end

%% get distance corrections
if isempty(ind)
    MLM92 = NaN;
    WGW94 = NaN;
    WGW96 = NaN;
    GS86 = NaN;
    GG91 = NaN;
    R35 = NaN;
    BJ84 = NaN;
    repi = NaN;
    rhyp = NaN;
    stns = [];
else % get corrections
    disp(length(ind))
    for j = 1:length(ind)
        % use Cuthberson's stationlist
        % get ML A0 correction
        if ~isnan(dep)
            [MLM92(j),WGW94(j),WGW96(j),BJ84(j),HB87(j),GS86(j),GG91(j),R35(j)] ...
             = getMLestimates(rhyp(j),dep,1,1);
        end            
    end
end

