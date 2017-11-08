function [rhyp,repi,MLM92,WGW94,WGW96,GS86,GG91,R35,BJ84,A10,stns] ...
         = get_stn_A0(siteDat,lat,lon,dep,evdate,stns,zone,essrng)

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
 
%% make proxy stationlist if ISC stns do not exist  
if isempty(stns) & isempty(essrng)
    [rng az] = distance([siteDat.stnlat],[siteDat.stnlon],lat,lon);
    % find stations within date & distnce range
    ind = find([siteDat.startdate] <= evdate & [siteDat.stopdate] >= evdate & deg2km(rng) < 1500);
    if ~isempty(ind)
        stns = {siteDat(ind).stncode};
%         disp (stns)
%         disp([lat lon])
%         repi = deg2km(rng(ind));
%         rhyp = sqrt(repi.^2 + (ones(size(repi))*dep).^2);
    end
% else % test if works!
%     [rng az] = distance([siteDat.stnlat],[siteDat.stnlon],lat,lon);
%     % find stations within date & distnce range
%     ind = find([siteDat.startdate] <= evdate & [siteDat.stopdate] >= evdate & deg2km(rng) < 1500);
%     if ~isempty(ind)
%         stns1 = {siteDat(ind).stncode};
%         repi = deg2km(rng(ind));
%         rhyp = sqrt(repi.^2 + (ones(size(repi))*dep).^2);
%         disp (stns)
%         disp (stns1)
%         disp(repi)
%         disp([lat lon])
%     else
%         repi = NaN;
%         rhyp = NaN;
%     end
end

%% get distance corrections
if isempty(stns)
    MLM92 = NaN;
    WGW94 = NaN;
    WGW96 = NaN;
    GS86 = NaN;
    GG91 = NaN;
    R35 = NaN;
    BJ84 = NaN;
    A10 = NaN;
    repi = NaN;
    rhyp = NaN;
else % get corrections
    for j = 1:length(stns)
        % match site
        siteind = find(strcmp({siteDat.stncode},stns(j)));
        % get distance for ansn
        if ~isempty(siteind) & isnan(essrng)
            [rng az] = distance(siteDat(siteind).stnlat,siteDat(siteind).stnlon,lat,lon);
            repi(j) = deg2km(rng);
            rhyp(j) = sqrt(deg2km(rng)^2 + dep^2);
            % get ML A0 correction
            if ~isnan(dep)
                [MLM92(j),WGW94(j),WGW96(j),BJ84(j),HB87(j),GS86(j),GG91(j),R35(j),A10(j)] ...
                 = getMLestimates(rhyp(j),dep,1,1);
            end
        % use es&s stns
        elseif ~isnan(essrng)
            repi(j) = essrng(j);
            rhyp(j) = sqrt(essrng(j)^2 + dep^2);
            % get ML A0 correction
            if ~isnan(dep)
                [MLM92(j),WGW94(j),WGW96(j),BJ84(j),HB87(j),GS86(j),GG91(j),R35(j),A10(j)] ...
                 = getMLestimates(rhyp(j),dep,1,1);
            end            
        else
            repi(j) = NaN;
            rhyp(j) = NaN;
            GG91(j) = NaN;
            MLM92(j) = NaN;
            R35(j) = NaN;
            BJ84(j) = NaN;
            GS86(j) = NaN;
            WGW94(j) = NaN;
            WGW96(j) = NaN;
        end
    end
end

