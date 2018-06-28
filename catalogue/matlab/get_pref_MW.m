% *************************************************************************
% Program: get_pref_MW.m
% 
% Coverts all mag types to MW and selects preferred MW
% 
% Corrections based on earthquake location, as defined by:
%	https://github.com/GeoscienceAustralia/NSHA2018/blob/master/catalogue/magnitude/ml/australia_ml_regions.txt
%
% zone = 1 > WA
% zone = 2 > EA
% zone = 3 > SA
% zone = 4 > outside Australia
%
% Author: T. Allen (2011-01-11) - updated February 2018
% 
% 2018-05-07: V0.2 - Adding column for fixed hunge ML2MW conversion
%
% *************************************************************************
outfile = fullfile('..','data','NSHA18CAT.MW.V0.2.csv');

% load data
if exist('mdat','var') ~= 1
    disp('Loading mdat');
    load mdat_ml_rev.mat;
end

MS2MW = ones(size(mdat)) * NaN;
mb2MW = ones(size(mdat)) * NaN;
%ML2MWA = ones(size(mdat)) * NaN;
ML2MWG = ones(size(mdat)) * NaN;
ML2MW_BLE = ones(size(mdat)) * NaN;
ML2MW_QDE = ones(size(mdat)) * NaN;

prefFinalMW = ones(size(mdat)) * NaN;

% *************************************************************************
%% Convert MS to MW using Di Giacomo et al (2015) for zone 4

% conversion for shallow earthquakes (h < 70 km)
disp('Converting MS to MW...');
% for 3.0 <= MS <= 6.1
ind1 = find([mdat.MDAT_dep] <= 70 & [mdat.MDAT_prefMS] >= 3.0 ...
           & [mdat.MDAT_prefMS] <= 6.47 & [mdat.zone] == 4);
% include those events with NaN depth       
ind2 = find(isnan([mdat.MDAT_dep]) & [mdat.MDAT_prefMS] >= 3.0 ...
           & [mdat.MDAT_prefMS] <= 6.47 & [mdat.zone] == 4);
% include those events h > 70 and no mb   
ind3 = find([mdat.MDAT_dep] > 70  & [mdat.MDAT_prefMS] >= 3.0 ...
           & [mdat.MDAT_prefMS] <= 6.47 & isnan([mdat.MDAT_prefmb]) ...
           & [mdat.zone] == 4);

ind = [ind1, ind2, ind3];       
MS2MW(ind) = 0.67 * [mdat(ind).MDAT_prefMS] + 2.13;

% for 6.47 < MS <= 8.2! (Ignore upper lim for now)
ind1 = find([mdat.MDAT_dep] <= 70 & [mdat.MDAT_prefMS] > 6.47 ...
             & [mdat.zone] == 4);
% include those events with NaN depth       
ind2 = find(isnan([mdat.MDAT_dep]) & [mdat.MDAT_prefMS] > 6.47 ...
             & [mdat.zone] == 4);
% include those events h > 70 and no mb   
ind3 = find([mdat.MDAT_dep] > 70  & [mdat.MDAT_prefMS] > 6.47 ...
            & isnan([mdat.MDAT_prefmb]) & [mdat.zone] == 4);

ind = [ind1, ind2, ind3];       
MS2MW(ind) = 1.10 * [mdat(ind).MDAT_prefMS] - 0.67;

% *************************************************************************
%% Convert MS to MW using Allen et al (in prep) for Aust events
% 2018 ODR coeff
c1 = 0.0755045389514;
c2 = 3.33414891107;
ind = find(~isnan([mdat.MDAT_prefMS]) & [mdat.zone] ~= 4);
MS2MW(ind) = c1 * [mdat(ind).MDAT_prefMS].^2 + c2;

% *************************************************************************
%% Convert mb to MW using Di Giacomo et al (2015) for zone 4
disp('Converting mb to MW...');
% for 3.5 <= mb <= 6.2
ind = find([mdat.MDAT_prefmb] >= 3.5 & [mdat.MDAT_prefmb] <= 6.2 ...
           & [mdat.zone] == 4);
mb2MW(ind) = 1.38 * [mdat(ind).MDAT_prefmb] - 1.79;

% *************************************************************************
%% Convert mb to MW using Allen (2012) for events below latitude -13 deg - out-dated!
% c1 = 0.7362;
% c2 = 0.7374;
% c3 = 0.9707;
% mx = 5.0905;

% ind = find([mdat.MDAT_prefmb] >= 3.0 & [mdat.MDAT_prefmb] < mx ...
%        & [mdat.zone] == 4);
% mb2MW(ind) = c1 * [mdat(ind).MDAT_prefmb] + c3;
% ind = find([mdat.MDAT_prefmb] >= mx & [mdat.MDAT_prefmb] <= 6.5 ...
%             & [mdat.zone] == 4);
% mb2MW(ind) = c1 * [mdat(ind).MDAT_prefmb] ...
%              + c2 * ([mdat(ind).MDAT_prefmb] - mx) + c3;
         
%% Convert mb to MW using Ghasemi (2017) for Aust events - outdated based on 2018 catalogue
% c1 = 1.1438907424442797;
% c2 = -0.87192285009579173;

% Allen et al (in prep) ODR coeffs
c1 = 1.20025959882;
c2 = -1.1760438127;

% for 3.0 <= MS <= 7.0
ind = find(~isnan([mdat.MDAT_prefmb]) & [mdat.zone] ~= 4);
mb2MW(ind) = c1 * [mdat(ind).MDAT_prefmb] + c2;

% *************************************************************************
%% Convert ML to MW onvert using Ghasemi (2017/18)

disp('Converting ML to MW using Ghasemi...');

% use ODR polynomial of simulated data from Ghasemi & Allen (2017)
a = 0.04160769;
b = 0.48058286;
c = 1.39485216;

% add ODR bi-linear coefs from empirical data - not used in NSHA18
a_bl = 0.66053496;
b_bl = 1.20883045;
c_bl = 0.98659071;
hx_bl = 4.25;
hy_bl =  a_bl * hx_bl + b_bl;

% add ODR quadratic coefs from empirical data - not used in NSHA18
a_qd = 0.08165208;
b_qd = 0.11965149;
c_qd = 2.08418142;


for i = 1:length(mdat)
    % calculate MW for non-revised ML
    if isnan(mdat(i).MDAT_MLrev)
        % get quadratic ML2MW from simulation
        ML2MWG(i) = a*mdat(i).MDAT_prefML^2 + b*mdat(i).MDAT_prefML + c;
        
        % get quadratic ML2MW from empirical
        ML2MW_QDE(i) = a_qd*mdat(i).MDAT_prefML^2 + b_qd*mdat(i).MDAT_prefML + c_qd;
        
        % calculate bi-linear empirical
        if mdat(i).MDAT_prefML <= hx_bl
            ML2MW_BLE(i) = a_bl*mdat(i).MDAT_prefML + b_bl;
        else
            ML2MW_BLE(i) = c_bl*(mdat(i).MDAT_prefML - hx_bl) + hy_bl;
        end

    % calculate MW for revised ML
    else
        % get quadratic ML2MW from simulation
        ML2MWG(i) = a*mdat(i).MDAT_MLrev^2 + b*mdat(i).MDAT_MLrev + c;
        
        % get quadratic ML2MW from empirical
        ML2MW_QDE(i) = a_qd*mdat(i).MDAT_MLrev^2 + b_qd*mdat(i).MDAT_MLrev + c_qd;
         
         % calculate bi-linear empirical
        if mdat(i).MDAT_MLrev <= hx_bl
            ML2MW_BLE(i) = a_bl*mdat(i).MDAT_MLrev + b_bl;
        else
            ML2MW_BLE(i) = c_bl*(mdat(i).MDAT_MLrev - hx_bl) + hy_bl;
        end
    end
end

%% set fields
for i = 1:length(mdat)
    mdat(i).MS2MW = MS2MW(i);
    mdat(i).mb2MW = mb2MW(i);
    %mdat(i).ML2MWA = ML2MWA(i);
    mdat(i).ML2MWG = ML2MWG(i);
    mdat(i).ML2MW_BLE = ML2MW_BLE(i);
    mdat(i).ML2MW_QDE = ML2MW_QDE(i);
    mdat(i).prefFinalMW = prefFinalMW(i);
end

% *************************************************************************
%% Set preferred MW
% *************************************************************************

% conserve actual Mw measurements first
for i = 1:length(mdat)
    if ~isnan(mdat(i).MDAT_prefMW)
        mdat(i).prefFinalMW = mdat(i).MDAT_prefMW;
        mdat(i).prefFinalMWSrc = mdat(i).MDAT_prefMWSrc;
        mdat(i).MDAT_origMagType = 'MW';

% take MS >= 5.75
    elseif mdat(i).MS2MW > 5.75
        mdat(i).prefFinalMW = mdat(i).MS2MW;
        mdat(i).prefFinalMWSrc = 'MS2MW';
        mdat(i).MDAT_origMagType = 'MS';
        
% take ML-MW
    elseif ~isnan(mdat(i).ML2MWG)
        mdat(i).prefFinalMW = mdat(i).ML2MWG;
        mdat(i).prefFinalMWSrc = 'ML2MW';
        mdat(i).MDAT_origMagType = mdat(i).MDAT_origMLType;

% take larger of MS/mb < 5.75        
    elseif ~isnan(mdat(i).MS2MW) | ~isnan(mdat(i).mb2MW)
        maxM = max([mdat(i).MDAT_prefMS mdat(i).MDAT_prefmb]); % deliberately use orig mag here
        if mdat(i).MDAT_prefMS == maxM
            mdat(i).prefFinalMW = mdat(i).MS2MW;
            mdat(i).prefFinalMWSrc = 'MS2MW';
            mdat(i).MDAT_origMagType = 'MS';
        elseif mdat(i).MDAT_prefmb == maxM
            mdat(i).prefFinalMW = mdat(i).mb2MW;
            mdat(i).prefFinalMWSrc = 'mb2MW';
            mdat(i).MDAT_origMagType = 'mb';
        end
% else, set as NaN    
    else
        mdat(i).prefFinalMW = NaN;
        mdat(i).prefFinalMWSrc = '';
        mdat(i).MDAT_origMagType = '';
    end 
    
% ADD GG DEPENDENCE
    if isempty(mdat(i).GG_sourceType)
        mdat(i).GG_sourceType = '';
    end
    if isempty(mdat(i).GG_dependence)
        mdat(i).GG_dependence = '';
    end
%     mdat(i).GG_dependence = mdat(i).GG_dependence;
end

% *************************************************************************
%% Set preferred non-MW

% conserve actual Mw measurements first
for i = 1:length(mdat)
    if ~isnan(mdat(i).MDAT_prefMW)
        mdat(i).Mx_OrigML = mdat(i).MDAT_prefMW;
        mdat(i).Mx_RevML = mdat(i).MDAT_prefMW;
        mdat(i).Mx_RevMLSrc = mdat(i).MDAT_prefMWSrc;
        mdat(i).Mx_RevMLtype = 'MW';
         
% take MS >= 5.75
    elseif mdat(i).MDAT_prefMS >= 5.75
        mdat(i).Mx_OrigML = mdat(i).MDAT_prefMS;
        mdat(i).Mx_RevML = mdat(i).MDAT_prefMS;
        mdat(i).Mx_RevMLSrc = mdat(i).MDAT_prefMSSrc;
        mdat(i).Mx_RevMLtype = 'MS';

% take Revised ML   
    elseif ~isnan(mdat(i).MDAT_MLrev)
        mdat(i).Mx_OrigML = mdat(i).MDAT_prefML;
        mdat(i).Mx_RevML = mdat(i).MDAT_MLrev;
        mdat(i).Mx_RevMLSrc = 'REV_ML';
        mdat(i).Mx_RevMLtype = 'REV_ML';

% take Original ML   
    elseif ~isnan(mdat(i).MDAT_prefML)
        mdat(i).Mx_OrigML = mdat(i).MDAT_prefML;
        mdat(i).Mx_RevML = NaN;
        mdat(i).Mx_RevMLSrc = mdat(i).MDAT_prefMLSrc;
        mdat(i).Mx_RevMLtype = 'ML';

% take larger of MS/mb < 5.75        
    elseif ~isnan(mdat(i).MDAT_prefMS) | ~isnan(mdat(i).MDAT_prefmb)
        maxM = max([mdat(i).MDAT_prefMS mdat(i).MDAT_prefmb]);
        if mdat(i).MDAT_prefMS == maxM
            mdat(i).Mx_OrigML = mdat(i).MDAT_prefMS;
            mdat(i).Mx_RevML = mdat(i).MDAT_prefMS;
            mdat(i).Mx_RevMLSrc = mdat(i).MDAT_prefMSSrc;
            mdat(i).Mx_RevMLtype = 'MS';
        elseif mdat(i).MDAT_prefmb == maxM
            mdat(i).Mx_OrigML = mdat(i).MDAT_prefmb;
            mdat(i).Mx_RevML = mdat(i).MDAT_prefmb;
            mdat(i).Mx_RevMLSrc = mdat(i).MDAT_prefmbSrc;
            mdat(i).Mx_RevMLtype = 'mb';
        end

	else
        mdat(i).Mx_OrigML = NaN;
        mdat(i).Mx_RevML = NaN;
        mdat(i).Mx_RevMLSrc = '';
        mdat(i).Mx_RevMLtype = '';
    end
end

% *************************************************************************
% clear mdat;
disp('Saving mdat');
save mdat_mw_pref mdat;

% *************************************************************************
%% write to file

header = 'DATESTR,DATENUM,TYPE,DEPENDENCE,LON,LAT,DEP,LOCSRC,PREFMW,PREFMWSRC,PREFMS,PREFMSSRC,PREFmb,PREFmbSRC,PREFML,PREFMLSRC,MLREGION,REVML,MX_ORIGML,MX_TYPE,MX_REVML,MX_REVMLTYPE,MX_REVMLSRC,MS2MW,mb2MW,ML2MW,ML2MW_BLE,ML2MW_QDE,PREFMW,PREFMWSRC,COMM';
disp('writing to file')
dlmwrite(outfile,header,'delimiter','');
txt = [];
for i = 1:length(mdat)
    comms = ['"',mdat(i).GG_place,'"'];
    
    % set mlregion
    if mdat(i).zone == 1
        mlreg = 'WA';
    elseif mdat(i).zone == 2
        mlreg = 'SEA';
    elseif mdat(i).zone == 3
        mlreg = 'SA';
    elseif mdat(i).zone == 4
        mlreg = 'Other';
    end
    
    % get value to fill Mx Rev ML column
    if isnan(mdat(i).Mx_RevML)
        Mx_RevML = mdat(i).Mx_OrigML;
    else
        Mx_RevML = mdat(i).Mx_RevML;
    end
    
    line = [datestr(mdat(i).MDAT_dateNum,31),',',num2str(mdat(i).MDAT_dateNum),',', ...
            mdat(i).GG_sourceType,',',num2str(mdat(i).GG_dependence),',', ...
            num2str(mdat(i).MDAT_lon),',', ...
            num2str(mdat(i).MDAT_lat),',',num2str(mdat(i).MDAT_dep),',', ...
            mdat(i).MDAT_locsrc,',', ...
            num2str(mdat(i).MDAT_prefMW),',',mdat(i).MDAT_prefMWSrc,',', ...
            num2str(mdat(i).MDAT_prefMS),',',mdat(i).MDAT_prefMSSrc,',', ...
            num2str(mdat(i).MDAT_prefmb),',',mdat(i).MDAT_prefmbSrc,',', ...
            num2str(mdat(i).MDAT_prefML),',',mdat(i).MDAT_prefMLSrc,',', mlreg,',', ...
            num2str(mdat(i).MDAT_MLrev,'%0.2f'),',',num2str(mdat(i).Mx_OrigML,'%0.2f'),',', ...
            mdat(i).MDAT_origMagType,',', ...
            num2str(Mx_RevML,'%0.2f'),',',mdat(i).Mx_RevMLtype,',', ...
            mdat(i).Mx_RevMLSrc,',', ...
            num2str(mdat(i).MS2MW,'%0.2f'),',',num2str(mdat(i).mb2MW,'%0.2f'),',', ...
            num2str(mdat(i).ML2MWG,'%0.2f'),',',num2str(mdat(i).ML2MW_BLE,'%0.2f'),',', ...
            num2str(mdat(i).ML2MW_QDE,'%0.2f'),',',num2str(mdat(i).prefFinalMW,'%0.2f'),',', ...
            mdat(i).prefFinalMWSrc,',',comms,char(10)];
    txt = [txt line];
end

% remove last line
txt = txt(1:end-1);
dlmwrite(outfile,txt,'delimiter','','-append');

% *************************************************************************
%% Make GMT mag diff file for pre-1990 events
clear txt;
ind = find([mdat.zone] ~= 4 & ~isnan([mdat.MDAT_prefML]) ...
      & ~isnan([mdat.MDAT_MLrev]));
mdiff = [mdat(ind).MDAT_MLrev] - [mdat(ind).MDAT_prefML];
dat = [[mdat(ind).MDAT_lon]' [mdat(ind).MDAT_lat]' mdiff' ...
      [mdat(ind).MDAT_MLrev]'/15];

dlmwrite('ML_diff.dat',dat,'delimiter','\t','precision','%0.3f');

% *************************************************************************
%% plot histograms of ML difference

figure(1);
ind = find([mdat.MDAT_dateNum] > datenum(1940,1,1) ...
           & [mdat.MDAT_dateNum] < datenum(1990,1,1));
mldiff = [mdat(ind).MDAT_prefML] - [mdat(ind).MDAT_MLrev];
subplot(1,2,1),hist(mldiff,[-1.6:0.1:1.6]);
nn = find(~isnan(mldiff));
medmldiff = median(mldiff(nn));
meanmldiff = mean(mldiff(nn));
xlim([-1.2 1.2]);
xlabel('Old ML (Richter) - Revised ML');
ylabel('Number of earthquakes');
title('ML residual (all magnitudes), 1940-1990');
text(-1.15,8400,['Events = ',num2str(length(mldiff(nn))),char(10) ...
                'Median ML Residual = ',num2str(medmldiff,'%0.3f'),char(10), ...
                'Mean ML Residual = ',num2str(meanmldiff,'%0.3f')],'Fontsize',8);
            
ind = find([mdat.MDAT_dateNum] > datenum(1940,1,1) ...
           & [mdat.MDAT_dateNum] < datenum(1990,1,1) ...
           & [mdat.MDAT_prefML] >= 4.0);
mldiff = [mdat(ind).MDAT_prefML] - [mdat(ind).MDAT_MLrev];
subplot(1,2,2),hist(mldiff,[-1.6:0.1:1.6]);
nn = find(~isnan(mldiff));
medmldiff = median(mldiff(nn));
meanmldiff = mean(mldiff(nn));
xlim([-1.2 1.2]);
xlabel('Old ML (Richter) - Revised ML');
ylabel('Number of earthquakes');
title('ML residual [ML (Richter) >= 4.0], 1940-1990');
text(-1.15,130,['Events = ',num2str(length(mldiff(nn))),char(10) ...
                'Median ML Residual = ',num2str(medmldiff,'%0.3f'),char(10), ...
                'Mean ML Residual = ',num2str(meanmldiff,'%0.3f')],'Fontsize',8);
