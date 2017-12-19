% Script to parse 2017 GG Cat and convert it to mdat format to use with old codes

% read csv
ggcatFile = '..\data\GGcat-170802_Aust_Local_EQ_M2.5.csv';
%ggcatFile = '../data/GGcat-170802_Aust_Local_EQ_M2.5.csv'; % for mac/linux
[auth, type, dependence, ass, year, month, day, hour, min, sec, ...
 timezone, timecor, lon, lat, dep, zcode, mx, mval, place, mtxt] ... 
 = textread(ggcatFile,'%s%s%f%s%f%f%f%f%f%f%s%f%f%f%f%s%s%f%s%s', ...
             'headerlines',1,'delimiter',',','emptyvalue',NaN);
         
% recreate mdat
for i = 1:length(auth)
    % check nan vals
    if isnan(hour(i)); hour(i)=0; end
    if isnan(min(i)); min(i)=0; end
    if isnan(sec(i)); sec(i)=0.0; end
    
    % fix time
    tmpdatetime = datenum(year(i), month(i), day(i), hour(i), min(i), sec(i)); % local time
    deltaT = timecor(i)/24;
    mdat(i).MDAT_dateNum = tmpdatetime; % already in UTC! - deltaT; % convert to UTC
    mdat(i).MDAT_dateStr = datestr(mdat(i).MDAT_dateNum, 31);
    
    mdat(i).MDAT_lat = -1*lat(i);
    mdat(i).MDAT_lon = lon(i);
    mdat(i).MDAT_dep = dep(i);
    mdat(i).MDAT_locsrc = auth{i};
    mdat(i).GG_sourceType = type{i};
    mdat(i).GG_Mtype = mx{i};
    mdat(i).GG_Mval = mval(i);
    mdat(i).GG_place = place{i};
    %mdat(i).GG_place = mtxt(i);
    mdat(i).GG_dependence = dependence(i);
    
end

save mdat mdat;