% *************************************************************************
% Program: reviseMLs.m
% 
% Aims to make consistent magnitudes across all time periods
% 
% zone = 1 > WA
% zone = 2 > EA
% zone = 3 > SA
% zone = 4 > Other
% 
% Author: T Allen (2011-01-06) - update 2018-01-03
% *************************************************************************

% load dataload ..\ML2MW\mdat11.mat;
if exist('mdat','var') ~= 1
    disp('Loading mdat_mag_types');
    load mdat_pref_mag_types.mat;
end

%% remove unecessary events
% delind = find(strcmp({mdat.MDAT_prefmbSrc},'ISC') & strcmp({mdat.MDAT_prefMLSrc},'AUST') == 0 ...
%          & strcmp({mdat.MDAT_prefMLSrc},'MEL') == 0 & strcmp({mdat.MDAT_prefMLSrc},'GG') == 0 ...
%          & strcmp({mdat.MDAT_prefMLSrc},'ADE') == 0 & strcmp({mdat.MDAT_prefMLSrc},'MGO') == 0 ...
%          & strcmp({mdat.MDAT_prefMLSrc},'MUN') == 0 & strcmp({mdat.MDAT_prefMLSrc},'AGSO') == 0 ...
%          & [mdat.MDAT_prefmb] < 5.0);
% mdat(delind) = [];
% delind = find(strcmp({mdat.MDAT_prefmbSrc},'IDC') & strcmp({mdat.MDAT_prefMLSrc},'AUST') == 0 ...
%          & strcmp({mdat.MDAT_prefMLSrc},'MEL') == 0 & strcmp({mdat.MDAT_prefMLSrc},'GG') == 0 ...
%          & strcmp({mdat.MDAT_prefMLSrc},'ADE') == 0 & strcmp({mdat.MDAT_prefMLSrc},'MGO') == 0 ...
%          & strcmp({mdat.MDAT_prefMLSrc},'MUN') == 0 & strcmp({mdat.MDAT_prefMLSrc},'AGSO') == 0 ...
%          & [mdat.MDAT_prefmb] < 5.0);
% mdat(delind) = [];
% delind = find(strcmp({mdat.MDAT_prefmbSrc},'ISC') & strcmp({mdat.MDAT_prefMLSrc},'AUST') == 0 ...
%          & strcmp({mdat.MDAT_prefMLSrc},'MEL') == 0 & strcmp({mdat.MDAT_prefMLSrc},'GG') == 0 ...
%          & strcmp({mdat.MDAT_prefMLSrc},'ADE') == 0 & strcmp({mdat.MDAT_prefMLSrc},'MGO') == 0 ...
%          & strcmp({mdat.MDAT_prefMLSrc},'MUN') == 0 & strcmp({mdat.MDAT_prefMLSrc},'AGSO') == 0 ...
%          & [mdat.MDAT_prefmb] < 5.0);
% mdat(delind) = [];
% delind = find(strcmp({mdat.MDAT_prefmbSrc},'IDC') & strcmp({mdat.MDAT_prefMLSrc},'AUST') == 0 ...
%          & strcmp({mdat.MDAT_prefMLSrc},'MEL') == 0 & strcmp({mdat.MDAT_prefMLSrc},'GG') == 0 ...
%          & strcmp({mdat.MDAT_prefMLSrc},'ADE') == 0 & strcmp({mdat.MDAT_prefMLSrc},'MGO') == 0 ...
%          & strcmp({mdat.MDAT_prefMLSrc},'MUN') == 0 & strcmp({mdat.MDAT_prefMLSrc},'AGSO') == 0 ...
%          & [mdat.MDAT_prefmb] < 5.0);
% mdat(delind) = [];
% delind = find(strcmp({mdat.MDAT_prefmbSrc},'DJA') & strcmp({mdat.MDAT_prefMLSrc},'AUST') == 0 ...
%          & strcmp({mdat.MDAT_prefMLSrc},'MEL') == 0 & strcmp({mdat.MDAT_prefMLSrc},'GG') == 0 ...
%          & strcmp({mdat.MDAT_prefMLSrc},'ADE') == 0 & strcmp({mdat.MDAT_prefMLSrc},'MGO') == 0 ...
%          & strcmp({mdat.MDAT_prefMLSrc},'MUN') == 0 & strcmp({mdat.MDAT_prefMLSrc},'AGSO') == 0 ...
%          & [mdat.MDAT_prefmb] < 5.0);
% mdat(delind) = [];
% delind = find(strcmp({mdat.MDAT_prefMLSrc},'DJA'));
% mdat(delind) = [];
% delind = find(strcmp({mdat.MDAT_prefMLSrc},'IDC'));
% mdat(delind) = [];

%%
% load mcorr coefs for converting old events based on average correction
mcorr = textread('mcorr.dat','%f','delimiter',',');

%% read station info
% disp('Reading station info...');
% sitefile = 'stationData\AUST.SeismographStations.txt';
% fid = fopen(sitefile);
% % skip header rows
% for i = 1:7
%     fgetl(fid);
% end
% 
% readsite = 1;
% i = 0;
% while readsite == 1
%     line = fgetl(fid);
%     if strcmp(line(1:4),'STOP') == 0
%         i = i + 1;
%         siteDat(i).stncode = deblank(line(1:7));
%         siteDat(i).stnlat = str2double(line(59:67));
%         siteDat(i).stnlon = str2double(line(69:78));
%         siteDat(i).stnelev = str2double(line(83:88));
%         fgetl(fid);fgetl(fid);
%     else
%         readsite = 0;
%     end
% end
% 
% %% get approximate start/stop times for ANSN
% disp('Getting station start/stop times...')
% for i = 1:length(siteDat)
%     tmparray = [];
%     for j = 1:length(mdat)
%         if ~isempty(mdat(j).ISC_stn)
%             ind = find(strcmp([mdat(j).ISC_stn],siteDat(i).stncode));
%             if ~isempty(ind)
%                 tmparray = [tmparray j];
%             end
%         end
%     end
%     if ~isempty(tmparray)
%         siteDat(i).startdate = mdat(tmparray(1)).MDAT_dateNum;
%         siteDat(i).stopdate = mdat(tmparray(end)).MDAT_dateNum;
%     else
%         siteDat(i).startdate = NaN;
%         siteDat(i).stopdate = NaN;
%     end
%     
%     if siteDat(i).stopdate > datenum(2007,7,1)
%         siteDat(i).stopdate = datenum(2011,1,1);
%     end
%     
%     if siteDat(i).startdate < datenum(1971,1,1)
%         siteDat(i).startdate = datenum(1945,1,1);
%     end
% end

%% get start/stop times for ES&S
disp('Getting ES&S station start/stop times...')

fid = fopen('Aus.net');
%[srcstns srclon srclat srcelv srcstarty,srcstartm,srcstopy,srcstopm] ...
%    = textscan(fid,'%s%f%f%f%s%s%d%d','Delimiter','\t');

[srcstns, srclon, srclat, srcelv, srcstarty,srcstartm,srcstopy,srcstopm] ...
    = textread('stationData\Aus.net','%s%f%f%f%s%s%s%s','delimiter','\t'); % file from RC

srcstart = [];
srcstop  = [];
for i = 1:length(srcstarty)   
    srcstart = [srcstart, datenum([srcstarty{i}, sprintf('%02d',srcstartm{i}), '01'],'yyyymmdd')];
    srcstop  = [srcstop,  datenum([srcstopy{i}, sprintf('%02d',srcstopm{i}), '01'],'yyyymmdd')];
end
%     
% srcstns = txt(:,1);
% for i = 1:length(mdat)
%     if isempty(mdat(i).ISC_stn) & strcmp(mdat(i).MDAT_prefMLSrc,'MEL')
%         [rng az] = distance(srclat,srclon,lat,lon);
%         ind = find(srcstart <= evdate & srcstop >= evdate & deg2km(rng) < 1500);
%         if ~isempty(ind)
%             srcstns = srcstns(ind);
%         end
%     end
% end

%% pre-define polygons as shown in Allen 2010 AEES manuscript

walon = [129  110. 110.  135.0  135.0  138.3  138.3  129];
walat = [-10 -18.5  -45   -45   -29   -29   -10   -10];
ealon = [138.3  138.3  141.0  141.0  155.5  155.5 145.5 138.3];
ealat = [-10   -29   -29   -45   -45  -18 -10   -10];
salon = [135   135   141   141   135];
salat = [-29   -45   -45   -29   -29];

%% loop through events to see if in EA or WA polygons for post-1960 events
disp('Looping thru events...')
 % for i = 57649:length(mdat)
 for i = 1:length(mdat)
    disp(datestr(mdat(i).MDAT_dateNum, 31));
       
    lat = mdat(i).MDAT_lat;
    lon = mdat(i).MDAT_lon;
    % check if in Western Australian polygon 
    if inpolygon(lon,lat,walon,walat)
        zone = 1;
    % check if in Eastern Australian polygon
    elseif inpolygon(lon,lat,ealon,ealat)
        zone = 2;
    % check if in South Australian polygon
    elseif inpolygon(lon,lat,salon,salat)
        zone = 3;
    else
        zone = 4;
    end

    dep = mdat(i).MDAT_dep;
    evdate = mdat(i).MDAT_dateNum;
    auth = mdat(i).MDAT_prefMLSrc;
    mdat(i).MDAT_MLrev = NaN;
    mdat(i).MDAT_MLrevdist = NaN;
    mdat(i).MDAT_MLminstn = [];
    
    % fix other mags to ML if ML unknown
%     if isnan(mdat(i).MDAT_prefML) & ~isnan(mdat(i).MDAT_otherM)
%         mdat(i).MDAT_prefML = mdat(i).MDAT_otherM;
%         mdat(i).MDAT_prefMLSrc = 'From Other';
%     elseif isnan(mdat(i).MDAT_prefML) & ~isnan(mdat(i).MDAT_prefmb)
%         mdat(i).MDAT_prefML = mdat(i).MDAT_prefmb;
%         mdat(i).MDAT_prefMLSrc = 'From mb';
%     end
    
    if zone < 4
        [rng az] = distance(srclat,srclon,lat,lon);
        ind = find(srcstart <= evdate & srcstop >= evdate & deg2km(rng') < 1500);
        if ~isempty(ind)
            stns = srcstns(ind);
            srcrng = deg2km(rng(ind));
        else
            srcrng = NaN;
        end
    else
        srcrng = NaN;
    end
    
    % get magnitude corrections
    if zone ~= 4 && ~isempty(ind)
        [rhyp,repi,MLM92_A0,WGW94_A0,WGW96_A0,GS86_A0,GG91_A0,R35_A0,BJ84_A0,stns] ...
        = get_stn_A0(siteDat,lat,lon,dep,evdate,stns,zone,srcrng);
    else
        % set corrections to be zero
        MLM92_A0 = 0.0;
        WGW94_A0 = 0.0;
        WGW96_A0 = 0.0;
        GS86_A0 = 0.0;
        GG91_A0 = 0.0;
        R35_A0 = 0.0;
        BJ84_A0 = 0.0;
        stns = [];
        rhyp = [];
        repi = [];
    end
    
%% get ML corrections assuming Western Australia Zone
    if zone == 1 & ~isempty(stns) & ~isnan(mdat(i).MDAT_prefML)
        
        % correct pre-1990 events to Gaull & Gregson assuming Richter (not ADE solutions)
        if mdat(i).MDAT_dateNum < datenum(1988,1,1) ...
           & strcmp(mdat(i).MDAT_prefMLSrc,'ADE') == 0
       
            % get station A
            R35_A = mdat(i).MDAT_prefML - R35_A0;
            % get revised mag
            gt50lt180 = find(rhyp >= 50 & rhyp < 180);
            if isempty(gt50lt180) % get minimum absolute Rhyp - 180 (try and avoid near-source data where possible)
                [dmindist dminind] = min(abs(rhyp-180));
                mdat(i).MDAT_MLrev = R35_A(dminind) + 1.137 * log10(rhyp(dminind)) ...
                                          + 0.000657 * rhyp(dminind) + 0.66;
                mdat(i).MDAT_MLrevdist = rhyp(dminind);                 
                mdat(i).MDAT_MLminstn = stns(dminind);
            elseif ~isempty(gt50lt180) % get between 50-180 km
                mdat(i).MDAT_MLrev = mean(R35_A(gt50lt180) + 1.137 * log10(rhyp(gt50lt180)) ...
                                          + 0.000657 * rhyp(gt50lt180) + 0.66);
                mdat(i).MDAT_MLrevdist = rhyp(gt50lt180);                 
                mdat(i).MDAT_MLminstn = stns(gt50lt180);
            end
        end
        
        % correct pre-1968 ADE events to Gaull & Gregson assuming Richter
        if mdat(i).MDAT_dateNum < datenum(1968,1,1) ...
           & strcmp(mdat(i).MDAT_prefMLSrc,'ADE') == 1
       
            % get station A
            R35_A = mdat(i).MDAT_prefML - R35_A0;
            % get revised mag
            gt50lt180 = find(rhyp >= 50 & rhyp < 180);
            if isempty(gt50lt180) % get minimum absolute Rhyp - 180 (try and avoid near-source data where possible)
                [dmindist dminind] = min(abs(rhyp-180));
                mdat(i).MDAT_MLrev = R35_A(dminind) + 1.137 * log10(rhyp(dminind)) ...
                                          + 0.000657 * rhyp(dminind) + 0.66;
                mdat(i).MDAT_MLrevdist = rhyp(dminind);                 
                mdat(i).MDAT_MLminstn = stns(dminind);
            elseif ~isempty(gt50lt180) % get between 50-180 km
                [dmindist dminind] = min(abs(rhyp-180));
                mdat(i).MDAT_MLrev = mean(R35_A(gt50lt180) + 1.137 * log10(rhyp(gt50lt180)) ...
                                          + 0.000657 * rhyp(gt50lt180) + 0.66);
                mdat(i).MDAT_MLrevdist = rhyp(gt50lt180);                 
                mdat(i).MDAT_MLminstn = stns(gt50lt180);
            end
        end
                
%% get ML corrections assuming Eastern Australia Zone
    elseif zone == 2 & ~isempty(stns) & ~isnan(mdat(i).MDAT_prefML)
        
        % keep Allen solutions as is
        if strcmp(mdat(i).MDAT_prefMLSrc,'Allen (unpublished)')
            mdat(i).MDAT_MLrev = mdat(i).MDAT_prefML;
            mdat(i).MDAT_MLrevdist = NaN;
            mdat(i).MDAT_MLminstn = [];
        end

        % correct pre-1990 events to Michael-Leiba & Manafant assuming Richter (not MEL solutions)    
        if ~isnan(mdat(i).MDAT_prefML) & mdat(i).MDAT_dateNum < datenum(1990,1,1) ...
           & strcmp(mdat(i).MDAT_prefMLSrc,'MEL') == 0
            % get station A
            R35_A = mdat(i).MDAT_prefML - R35_A0;
            % get revised mag
            gt50lt180 = find(rhyp >= 50 & rhyp < 180);
            if isempty(gt50lt180) % get minimum absolute Rhyp - 180 (try and avoid near-source data where possible)
                [dmindist dminind] = min(abs(rhyp-180));
                mdat(i).MDAT_MLrev = R35_A(dminind) + 1.34*log10(rhyp(dminind)/100) ...
                                          + 0.00055*(rhyp(dminind)-100)+ 3.0; % changed from 3.13 as assumed used maxh                                      
                mdat(i).MDAT_MLrevdist = rhyp(dminind);                 
                mdat(i).MDAT_MLminstn = stns(dminind);
            elseif ~isempty(gt50lt180) % get between 50-180 km
                [dmindist dminind] = min(abs(rhyp-180));
                mdat(i).MDAT_MLrev = mean(R35_A(gt50lt180) + 1.34*log10(rhyp(gt50lt180)/100) ...
                                          + 0.00055*(rhyp(gt50lt180)-100)+ 3.0); % changed from 3.13 as assumed used maxh
                mdat(i).MDAT_MLrevdist = rhyp(gt50lt180);                 
                mdat(i).MDAT_MLminstn = stns(gt50lt180);
            end
        end
        
        % correct pre-1994 MEL events to Michael-Leiba & Manafant assuming Richter   
        if ~isnan(mdat(i).MDAT_prefML) & mdat(i).MDAT_dateNum < datenum(1994,1,1) ...
           & strcmp(mdat(i).MDAT_prefMLSrc,'MEL') == 1
            % get station A
            R35_A = mdat(i).MDAT_prefML - R35_A0;
            % get revised mag
            gt50lt180 = find(rhyp >= 50 & rhyp < 180);
            if isempty(gt50lt180) % get minimum absolute Rhyp - 180 (try and avoid near-source data where possible)
                [dmindist dminind] = min(abs(rhyp-180));
                mdat(i).MDAT_MLrev = R35_A(dminind) + 1.34*log10(rhyp(dminind)/100) ...
                                          + 0.00055*(rhyp(dminind)-100)+ 3.0; % changed from 3.13 as assumed used maxh                                      
                mdat(i).MDAT_MLrevdist = rhyp(dminind);                 
                mdat(i).MDAT_MLminstn = stns(dminind);
            elseif ~isempty(gt50lt180) % get between 50-180 km
                [dmindist dminind] = min(abs(rhyp-180));
                mdat(i).MDAT_MLrev = mean(R35_A(gt50lt180) + 1.34*log10(rhyp(gt50lt180)/100) ...
                                          + 0.00055*(rhyp(gt50lt180)-100)+ 3.0); % changed from 3.13 as assumed used maxh
                mdat(i).MDAT_MLrevdist = rhyp(gt50lt180);                 
                mdat(i).MDAT_MLminstn = stns(gt50lt180);
            end
        end
        
        % correct pre-1994 Cuthbertson events to Michael-Leiba & Manafant assuming Richter   
        if ~isnan(mdat(i).MDAT_prefML) & mdat(i).MDAT_dateNum < datenum(1994,1,1) ...
           & strcmp(mdat(i).MDAT_prefMLSrc,'RC') == 1
            % get station A
            R35_A = mdat(i).MDAT_prefML - R35_A0;
            % get revised mag
            gt50lt180 = find(rhyp >= 50 & rhyp < 180);
            if isempty(gt50lt180) % get minimum absolute Rhyp - 180 (try and avoid near-source data where possible)
                [dmindist dminind] = min(abs(rhyp-180));
                mdat(i).MDAT_MLrev = R35_A(dminind) + 1.34*log10(rhyp(dminind)/100) ...
                                          + 0.00055*(rhyp(dminind)-100)+ 3.0; % changed from 3.13 as assumed used maxh                                      
                mdat(i).MDAT_MLrevdist = rhyp(dminind);                 
                mdat(i).MDAT_MLminstn = stns(dminind);
            elseif ~isempty(gt50lt180) % get between 50-180 km
                [dmindist dminind] = min(abs(rhyp-180));
                mdat(i).MDAT_MLrev = mean(R35_A(gt50lt180) + 1.34*log10(rhyp(gt50lt180)/100) ...
                                          + 0.00055*(rhyp(gt50lt180)-100)+ 3.0); % changed from 3.13 as assumed used maxh
                mdat(i).MDAT_MLrevdist = rhyp(gt50lt180);                 
                mdat(i).MDAT_MLminstn = stns(gt50lt180);
            end
        end
        
        % correct pre-1994 Gibson events to Michael-Leiba & Manafant assuming Richter   
        if ~isnan(mdat(i).MDAT_prefML) & mdat(i).MDAT_dateNum < datenum(1994,1,1) ...
           & strcmp(mdat(i).MDAT_prefMLSrc,'GG') == 1
            % get station A
            R35_A = mdat(i).MDAT_prefML - R35_A0;
            % get revised mag
            gt50lt180 = find(rhyp >= 50 & rhyp < 180);
            if isempty(gt50lt180) % get minimum absolute Rhyp - 180 (try and avoid near-source data where possible)
                [dmindist dminind] = min(abs(rhyp-180));
                mdat(i).MDAT_MLrev = R35_A(dminind) + 1.34*log10(rhyp(dminind)/100) ...
                                          + 0.00055*(rhyp(dminind)-100)+ 3.0; % changed from 3.13 as assumed used maxh                                      
                mdat(i).MDAT_MLrevdist = rhyp(dminind);                 
                mdat(i).MDAT_MLminstn = stns(dminind);
            elseif ~isempty(gt50lt180) % get between 50-180 km
                [dmindist dminind] = min(abs(rhyp-180));
                mdat(i).MDAT_MLrev = mean(R35_A(gt50lt180) + 1.34*log10(rhyp(gt50lt180)/100) ...
                                          + 0.00055*(rhyp(gt50lt180)-100)+ 3.0); % changed from 3.13 as assumed used maxh
                mdat(i).MDAT_MLrevdist = rhyp(gt50lt180);                 
                mdat(i).MDAT_MLminstn = stns(gt50lt180);
            end
        end
        
        % correct post-2002 MEL events to Michael-Leiba & Manafant assuming Bakun & Joyner (1984)
        if ~isnan(mdat(i).MDAT_prefML) & mdat(i).MDAT_dateNum > datenum(2002,1,1) ...
           & strcmp(mdat(i).MDAT_prefMLSrc,'MEL') == 1
            % get station A
            BJ84_A = mdat(i).MDAT_prefML - BJ84_A0;
            % get revised mag
            gt50lt180 = find(rhyp >= 50 & rhyp < 180);
            if isempty(gt50lt180) % get minimum absolute Rhyp - 180 (try and avoid near-source data where possible)
                [dmindist dminind] = min(abs(rhyp-180));
                mdat(i).MDAT_MLrev = BJ84_A(dminind) + 1.34*log10(rhyp(dminind)/100) ...
                                          + 0.00055*(rhyp(dminind)-100)+ 3.0; % changed from 3.13 as assumed used maxh                                      
                mdat(i).MDAT_MLrevdist = rhyp(dminind);                 
                mdat(i).MDAT_MLminstn = stns(dminind);
            elseif ~isempty(gt50lt180) % get between 50-180 km
                [dmindist dminind] = min(abs(rhyp-180));
                mdat(i).MDAT_MLrev = mean(BJ84_A(gt50lt180) + 1.34*log10(rhyp(gt50lt180)/100) ...
                                          + 0.00055*(rhyp(gt50lt180)-100)+ 3.0); % changed from 3.13 as assumed used maxh
                mdat(i).MDAT_MLrevdist = rhyp(gt50lt180);                 
                mdat(i).MDAT_MLminstn = stns(gt50lt180);
            end
        end
        
        % correct post-2002 RC events to Michael-Leiba & Manafant assuming Bakun & Joyner (1984)
        if ~isnan(mdat(i).MDAT_prefML) & mdat(i).MDAT_dateNum > datenum(2002,1,1) ...
           & strcmp(mdat(i).MDAT_prefMLSrc,'RC') == 1
            % get station A
            BJ84_A = mdat(i).MDAT_prefML - BJ84_A0;
            % get revised mag
            gt50lt180 = find(rhyp >= 50 & rhyp < 180);
            if isempty(gt50lt180) % get minimum absolute Rhyp - 180 (try and avoid near-source data where possible)
                [dmindist dminind] = min(abs(rhyp-180));
                mdat(i).MDAT_MLrev = BJ84_A(dminind) + 1.34*log10(rhyp(dminind)/100) ...
                                          + 0.00055*(rhyp(dminind)-100)+ 3.0; % changed from 3.13 as assumed used maxh                                      
                mdat(i).MDAT_MLrevdist = rhyp(dminind);                 
                mdat(i).MDAT_MLminstn = stns(dminind);
            elseif ~isempty(gt50lt180) % get between 50-180 km
                [dmindist dminind] = min(abs(rhyp-180));
                mdat(i).MDAT_MLrev = mean(BJ84_A(gt50lt180) + 1.34*log10(rhyp(gt50lt180)/100) ...
                                          + 0.00055*(rhyp(gt50lt180)-100)+ 3.0); % changed from 3.13 as assumed used maxh
                mdat(i).MDAT_MLrevdist = rhyp(gt50lt180);                 
                mdat(i).MDAT_MLminstn = stns(gt50lt180);
            end
        end
        
        % correct post-2002 GG events to Michael-Leiba & Manafant assuming Bakun & Joyner (1984)
        if ~isnan(mdat(i).MDAT_prefML) & mdat(i).MDAT_dateNum > datenum(2002,1,1) ...
           & strcmp(mdat(i).MDAT_prefMLSrc,'GG') == 1
            % get station A
            BJ84_A = mdat(i).MDAT_prefML - BJ84_A0;
            % get revised mag
            gt50lt180 = find(rhyp >= 50 & rhyp < 180);
            if isempty(gt50lt180) % get minimum absolute Rhyp - 180 (try and avoid near-source data where possible)
                [dmindist dminind] = min(abs(rhyp-180));
                mdat(i).MDAT_MLrev = BJ84_A(dminind) + 1.34*log10(rhyp(dminind)/100) ...
                                          + 0.00055*(rhyp(dminind)-100)+ 3.0; % changed from 3.13 as assumed used maxh                                      
                mdat(i).MDAT_MLrevdist = rhyp(dminind);                 
                mdat(i).MDAT_MLminstn = stns(dminind);
            elseif ~isempty(gt50lt180) % get between 50-180 km
                [dmindist dminind] = min(abs(rhyp-180));
                mdat(i).MDAT_MLrev = mean(BJ84_A(gt50lt180) + 1.34*log10(rhyp(gt50lt180)/100) ...
                                          + 0.00055*(rhyp(gt50lt180)-100)+ 3.0); % changed from 3.13 as assumed used maxh
                mdat(i).MDAT_MLrevdist = rhyp(gt50lt180);                 
                mdat(i).MDAT_MLminstn = stns(gt50lt180);
            end
        end
        
        % catch remaining events post-1990 (e.g. QEDB, QDM, GSQ, CQSRG)
        if ~isnan(mdat(i).MDAT_prefML) & mdat(i).MDAT_dateNum >= datenum(1990,1,1) ...
           & strcmp(mdat(i).MDAT_prefMLSrc,'MEL') == 0 & strcmp(mdat(i).MDAT_prefMLSrc,'AUST') == 0 ...
           & strcmp(mdat(i).MDAT_prefMLSrc,'RC') == 0 & strcmp(mdat(i).MDAT_prefMLSrc,'GG') == 0 ...
           & strcmp(mdat(i).MDAT_prefMLSrc,'Allen (unpublished)') == 0
           disp(mdat(i).MDAT_prefMLSrc);
            % get station A
            R35_A = mdat(i).MDAT_prefML - R35_A0;
            % get revised mag
            gt50lt180 = find(rhyp >= 50 & rhyp < 180);
            if isempty(gt50lt180) % get minimum absolute Rhyp - 180 (try and avoid near-source data where possible)
                [dmindist dminind] = min(abs(rhyp-180));
                mdat(i).MDAT_MLrev = R35_A(dminind) + 1.34*log10(rhyp(dminind)/100) ...
                                          + 0.00055*(rhyp(dminind)-100)+ 3.0; % changed from 3.13 as assumed used maxh                                      
                mdat(i).MDAT_MLrevdist = rhyp(dminind);                 
                mdat(i).MDAT_MLminstn = stns(dminind);
            elseif ~isempty(gt50lt180) % get between 50-180 km
                [dmindist dminind] = min(abs(rhyp-180));
                mdat(i).MDAT_MLrev = mean(R35_A(gt50lt180) + 1.34*log10(rhyp(gt50lt180)/100) ...
                                          + 0.00055*(rhyp(gt50lt180)-100)+ 3.0); % changed from 3.13 as assumed used maxh
                mdat(i).MDAT_MLrevdist = rhyp(gt50lt180);                 
                mdat(i).MDAT_MLminstn = stns(gt50lt180);
            end
        end
        
%         % correct 1994-1996 MEL events to Michael-Leiba & Manafant assuming Wilkie et al 94  - use restricted range
%         if ~isnan(mdat(i).MDAT_prefML) & mdat(i).MDAT_dateNum > datenum(1994,1,1) ...
%            & mdat(i).MDAT_dateNum < datenum(1996,1,1) & strcmp(mdat(i).MDAT_prefMLSrc,'MEL') == 1
%             % get station A
%             WGW94_A = mdat(i).MDAT_prefML - WGW94_A0;
%             % get revised mag
%             gt50lt180 = find(rhyp >= 50 & rhyp < 110);
%             if isempty(gt50lt180) % get minimum absolute Rhyp - 180 (try and avoid near-source data where possible)
%                 [dmindist dminind] = min(abs(rhyp-180));
%                 mdat(i).MDAT_MLrev = WGW94_A(dminind) + 1.34*log10(rhyp(dminind)/100) ...
%                                           + 0.00055*(rhyp(dminind)-100)+ 3.0; % changed from 3.13 as assumed used maxh                                      
%                 mdat(i).MDAT_MLrevdist = rhyp(dminind);                 
%                 mdat(i).MDAT_MLminstn = stns(dminind);
%             elseif ~isempty(gt50lt180) % get between 50-180 km
%                 [dmindist dminind] = min(abs(rhyp-180));
%                 mdat(i).MDAT_MLrev = mean(WGW94_A(gt50lt180) + 1.34*log10(rhyp(gt50lt180)/100) ...
%                                           + 0.00055*(rhyp(gt50lt180)-100)+ 3.0); % changed from 3.13 as assumed used maxh
%                 mdat(i).MDAT_MLrevdist = rhyp(gt50lt180);                 
%                 mdat(i).MDAT_MLminstn = stns(gt50lt180);
%             end
%         end
        
        % correct 1994-2002 MEL events to Michael-Leiba & Manafant assuming Wilkie et al 96  - use restricted range
        if ~isnan(mdat(i).MDAT_prefML) && mdat(i).MDAT_dateNum > datenum(1994,1,1) ...
           & mdat(i).MDAT_dateNum < datenum(2002,1,1) & strcmp(mdat(i).MDAT_prefMLSrc,'MEL') == 1
            % get station A
            WGW96_A = mdat(i).MDAT_prefML - WGW96_A0;
            % get revised mag
            gt50lt180 = find(rhyp >= 50 & rhyp < 110);
            if isempty(gt50lt180) % get minimum absolute Rhyp - 180 (try and avoid near-source data where possible)
                [dmindist, dminind] = min(abs(rhyp-180));
                mdat(i).MDAT_MLrev = WGW96_A(dminind) + 1.34*log10(rhyp(dminind)/100) ...
                                          + 0.00055*(rhyp(dminind)-100)+ 3.0; % changed from 3.13 as assumed used maxh                                      
                mdat(i).MDAT_MLrevdist = rhyp(dminind);                 
                mdat(i).MDAT_MLminstn = stns(dminind);
            elseif ~isempty(gt50lt180) % get between 50-180 km
                [dmindist, dminind] = min(abs(rhyp-180));
                mdat(i).MDAT_MLrev = mean(WGW96_A(gt50lt180) + 1.34*log10(rhyp(gt50lt180)/100) ...
                                          + 0.00055*(rhyp(gt50lt180)-100)+ 3.0); % changed from 3.13 as assumed used maxh
                mdat(i).MDAT_MLrevdist = rhyp(gt50lt180);                 
                mdat(i).MDAT_MLminstn = stns(gt50lt180);
            end
        end

%% get ML corrections assuming South Australia Zone
    elseif zone == 3 & ~isempty(stns) & ~isnan(mdat(i).MDAT_prefML)
        
        % correct pre-1968 ADE events to Greenhalgh 1986 assuming Richter
        if ~isnan(mdat(i).MDAT_prefML) & mdat(i).MDAT_dateNum < datenum(1968,1,1) ...
           & strcmp(mdat(i).MDAT_prefMLSrc,'ADE') == 1
            % get station A
            R35_A = mdat(i).MDAT_prefML - R35_A0;
            % get revised mag
            gt50lt180 = find(rhyp >= 50 & rhyp < 180);
            if isempty(gt50lt180) % get minimum absolute Rhyp - 180 (try and avoid near-source data where possible)
                [dmindist dminind] = min(abs(rhyp-180));
                mdat(i).MDAT_MLrev = R35_A(dminind) + 1.10 * log10(repi(dminind)/100) ...
                                          + 0.0013 * (repi(dminind) - 100) + 3.03;
                mdat(i).MDAT_MLrevdist = rhyp(dminind);                 
                mdat(i).MDAT_MLminstn = stns(dminind);
            elseif ~isempty(gt50lt180) % get between 50-180 km
                [dmindist dminind] = min(abs(rhyp-180));
                mdat(i).MDAT_MLrev = mean(R35_A(gt50lt180) + 1.10 * log10(repi(gt50lt180)/100) ...
                                          + 0.0013 * (repi(gt50lt180) - 100) + 3.03);
                mdat(i).MDAT_MLrevdist = rhyp(gt50lt180);                 
                mdat(i).MDAT_MLminstn = stns(gt50lt180);
            end
        end
        
        % correct pre-1986 non-ADE events to Greenhalgh 1986 assuming Richter
        if ~isnan(mdat(i).MDAT_prefML) & mdat(i).MDAT_dateNum < datenum(1968,1,1) ...
           & strcmp(mdat(i).MDAT_prefMLSrc,'ADE') == 0
            % get station A
            R35_A = mdat(i).MDAT_prefML - R35_A0;
            % get revised mag
            gt50lt180 = find(rhyp >= 50 & rhyp < 180);
            if isempty(gt50lt180) % get minimum absolute Rhyp - 180 (try and avoid near-source data where possible)
                [dmindist dminind] = min(abs(rhyp-180));
                mdat(i).MDAT_MLrev = R35_A(dminind) + 1.10 * log10(repi(dminind)/100) ...
                                          + 0.0013 * (repi(dminind) - 100) + 3.03;
                mdat(i).MDAT_MLrevdist = rhyp(dminind);                 
                mdat(i).MDAT_MLminstn = stns(dminind);
            elseif ~isempty(gt50lt180) % get between 50-180 km
                [dmindist dminind] = min(abs(rhyp-180));
                mdat(i).MDAT_MLrev = mean(R35_A(gt50lt180) + 1.10 * log10(repi(gt50lt180)/100) ...
                                          + 0.0013 * (repi(gt50lt180) - 100) + 3.03);
                mdat(i).MDAT_MLrevdist = rhyp(gt50lt180);                 
                mdat(i).MDAT_MLminstn = stns(gt50lt180);
            end
        end
        
        % correct post-2002 ADE events to Greenhalgh 1986 assuming Bakun & Joyner (1984)
        if ~isnan(mdat(i).MDAT_prefML) & mdat(i).MDAT_dateNum > datenum(2002,1,1) ...
           & strcmp(mdat(i).MDAT_prefMLSrc,'ADE') == 1
            % get station A
            BJ84_A = mdat(i).MDAT_prefML - BJ84_A0;
            % get revised mag
            gt50lt180 = find(rhyp >= 50 & rhyp < 180);
            if isempty(gt50lt180) % get minimum absolute Rhyp - 180 (try and avoid near-source data where possible)
                [dmindist dminind] = min(abs(rhyp-180));
                mdat(i).MDAT_MLrev = BJ84_A(dminind) + 1.10 * log10(repi(dminind)/100) ...
                                          + 0.0013 * (repi(dminind) - 100) + 3.03;                                      
                mdat(i).MDAT_MLrevdist = rhyp(dminind);                 
                mdat(i).MDAT_MLminstn = stns(dminind);
            elseif ~isempty(gt50lt180) % get between 50-180 km
                [dmindist dminind] = min(abs(rhyp-180));
                mdat(i).MDAT_MLrev = mean(BJ84_A(gt50lt180) + 1.10 * log10(repi(gt50lt180)/100) ...
                                          + 0.0013 * (repi(gt50lt180) - 100) + 3.03); % changed from 3.13 as assumed used maxh
                mdat(i).MDAT_MLrevdist = rhyp(gt50lt180);                 
                mdat(i).MDAT_MLminstn = stns(gt50lt180);
            end
        end
        
    %% get ML corrections for unknown magnitude types assuming Central & Western Australia Zone
%     elseif zone == 1 & ~isnan(mdat(i).MDAT_otherM)
%         
%         % correct WA events to Gaull & Gregson assuming Richter
%         % get station A
%         R35_A = mdat(i).MDAT_otherM - R35_A0;
%         % get revised mag
%         if ~isempty(stns)
%             gt50lt180 = find(rhyp >= 50 & rhyp < 180);
%             if isempty(gt50lt180) % get minimum absolute Rhyp - 180 (try and avoid near-source data where possible)
%                 [dmindist dminind] = min(abs(rhyp-180));
%                 mdat(i).MDAT_MLrev = R35_A(dminind) + 1.137 * log10(rhyp(dminind)) ...
%                                           + 0.000657 * rhyp(dminind) + 0.66;
%                 mdat(i).MDAT_MLrevdist = rhyp(dminind);                 
%                 mdat(i).MDAT_MLminstn = stns(dminind);
%             elseif ~isempty(gt50lt180) % get between 50-180 km
%                 mdat(i).MDAT_MLrev = mean(R35_A(gt50lt180) + 1.137 * log10(rhyp(gt50lt180)) ...
%                                           + 0.000657 * rhyp(gt50lt180) + 0.66);
%                 mdat(i).MDAT_MLrevdist = rhyp(gt50lt180);                 
%                 mdat(i).MDAT_MLminstn = stns(gt50lt180);
%             end
%         else
%             mdat(i).MDAT_MLrev = mdat(i).MDAT_otherM;
%             mdat(i).MDAT_MLrevdist = NaN;
%             mdat(i).MDAT_MLminstn = [];
%         end
%         
%     %% get ML corrections for unknown magnitude types assuming Eastern Australia Zone
%     elseif zone == 2 & ~isnan(mdat(i).MDAT_otherM)
%         
%         % correct EA events to MLM92 assuming Richter
%         % get station A
%         R35_A = mdat(i).MDAT_otherM - R35_A0;
%         % get revised mag
%         if ~isempty(stns)
%             gt50lt180 = find(rhyp >= 50 & rhyp < 180);
%             if isempty(gt50lt180) % get minimum absolute Rhyp - 180 (try and avoid near-source data where possible)
%                 [dmindist dminind] = min(abs(rhyp-180));
%                 mdat(i).MDAT_MLrev = R35_A(dminind) + 1.34*log10(rhyp(dminind)/100) ...
%                                           + 0.00055*(rhyp(dminind)-100)+ 3.0; % changed from 3.13 as assumed used maxh                                      
%                 mdat(i).MDAT_MLrevdist = rhyp(dminind);                 
%                 mdat(i).MDAT_MLminstn = stns(dminind);
%             elseif ~isempty(gt50lt180) % get between 50-180 km
%                 [dmindist dminind] = min(abs(rhyp-180));
%                 mdat(i).MDAT_MLrev = mean(R35_A(gt50lt180) + 1.34*log10(rhyp(gt50lt180)/100) ...
%                                           + 0.00055*(rhyp(gt50lt180)-100)+ 3.0); % changed from 3.13 as assumed used maxh
%                 mdat(i).MDAT_MLrevdist = rhyp(gt50lt180);                 
%                 mdat(i).MDAT_MLminstn = stns(gt50lt180);
%             end
%         else
%             mdat(i).MDAT_MLrev = mdat(i).MDAT_otherM;
%             mdat(i).MDAT_MLrevdist = NaN;
%             mdat(i).MDAT_MLminstn = [];
%         end
%         
%     %% get ML corrections for unknown magnitude types assuming South Australia Zone
%     elseif zone == 3 & ~isnan(mdat(i).MDAT_otherM)
%         
%         % correct SA events to Greenhalgh 1986 assuming Richter
%         % get station A
%         R35_A = mdat(i).MDAT_otherM - R35_A0;
%         % get revised mag
%         if ~isempty(stns)
%             gt50lt180 = find(rhyp >= 50 & rhyp < 180);
%             if isempty(gt50lt180) % get minimum absolute Rhyp - 180 (try and avoid near-source data where possible)
%                 [dmindist dminind] = min(abs(rhyp-180));
%                 mdat(i).MDAT_MLrev = R35_A(dminind) + 1.10 * log10(repi(dminind)/100) ...
%                                           + 0.0013 * (repi(dminind) - 100) + 3.03;
%                 mdat(i).MDAT_MLrevdist = rhyp(dminind);                 
%                 mdat(i).MDAT_MLminstn = stns(dminind);
%             elseif ~isempty(gt50lt180) % get between 50-180 km
%                 [dmindist dminind] = min(abs(rhyp-180));
%                 mdat(i).MDAT_MLrev = mean(R35_A(gt50lt180) + 1.10 * log10(repi(gt50lt180)/100) ...
%                                           + 0.0013 * (repi(gt50lt180) - 100) + 3.03);
%                 mdat(i).MDAT_MLrevdist = rhyp(gt50lt180);                 
%                 mdat(i).MDAT_MLminstn = stns(gt50lt180);
%             end
%         else
%             mdat(i).MDAT_MLrev = mdat(i).MDAT_otherM;
%             mdat(i).MDAT_MLrevdist = NaN;
%             mdat(i).MDAT_MLminstn = [];
%         end
        
%% Not in WA, SA or EA        
    else zone > 3 & ~isnan(mdat(i).MDAT_prefML)
        mdat(i).MDAT_MLrev = mdat(i).MDAT_prefML;
        mdat(i).MDAT_MLrevdist = NaN;
        mdat(i).MDAT_MLminstn = [];
    end
    mdat(i).zone = zone;
    
%% Use magnitude correction for those events that have no sites
    if isempty(stns) & zone ~= 4
        if isnan(mdat(i).MDAT_MLrevdist) & ~isnan(mdat(i).MDAT_prefML)
            mdat(i).MDAT_MLrev = mcorr(1) * mdat(i).MDAT_prefML + mcorr(2);
        end
        mdat(i).MDAT_MLrevdist = NaN;
        mdat(i).MDAT_MLminstn = [];
    end
end

% to fix empty structs
for i=1:length(mdat)
    if isnan(mdat(i).MDAT_MLrev) | mdat(i).MDAT_MLrev > -99
        % do nothing
    else
        mdat(i).MDAT_MLrev = NaN;
    end
end

save mdat_ml_rev mdat;
