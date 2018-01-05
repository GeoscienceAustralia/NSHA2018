% *************************************************************************
% Program: get_pref_MW.m
% 
% Coverts all mag types to MW and selects preferred MW
% 
% zone = 1 > WA
% zone = 2 > EA
% zone = 3 > SA
%
% Author: T. Allen (2011-01-11)
% *************************************************************************
outfile = '..\data\NSHA18CAT.MW.V0.1.csv';

% load data

if exist('mdat','var') ~= 1
    disp('Loading mdat');
    load mdat_ml_rev.mat;
end

% if exist('mdat','var') ~= 1
%     disp('Loading mdat');
%     load '..\Merge MC and QUAKES\mdat.mat';
% end

MS2MW = ones(size(mdat)) * NaN;
mb2MW = ones(size(mdat)) * NaN;
ML2MWA = ones(size(mdat)) * NaN;
ML2MWG = ones(size(mdat)) * NaN;
prefFinalMW = ones(size(mdat)) * NaN;

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

%% Convert MS to MW using Ghasemi (2017) for Aust events
c1 = 0.84896727404297323;
c2 = 1.0509630268292971;

% for 3.0 <= MS <= 7.0
ind = find(~isnan([mdat.MDAT_prefMS]) & [mdat.zone] ~= 4);
MS2MW(ind) = c1 * [mdat(ind).MDAT_prefMS] + c2;

%% Convert mb to MW using Di Giacomo et al (2015) for zone 4

disp('Converting mb to MW...');
% for 3.5 <= mb <= 6.2
ind = find([mdat.MDAT_prefmb] >= 3.5 & [mdat.MDAT_prefmb] <= 6.2 ...
           & [mdat.zone] == 4);
mb2MW(ind) = 1.38 * [mdat(ind).MDAT_prefmb] - 1.79;

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
         
%% Convert mb to MW using Ghasemi (2017) for Aust events
c1 = 1.1438907424442797;
c2 = -0.87192285009579173;

% for 3.0 <= MS <= 7.0
ind = find(~isnan([mdat.MDAT_prefmb]) & [mdat.zone] ~= 4);
mb2MW(ind) = c1 * [mdat(ind).MDAT_prefmb] + c2;

%% Convert ML to MW using Allen conversions - out-dated, but preserve in catalogue
% mx = 4.2;

disp('Converting ML to MW in CWA...');
[a1,a2,a3,mx] = textread('F:\Catalogues\ML2MW\WA.ML-MW.coef.txt','%f%f%f%f','delimiter',',');
% for ML rev
ind = find([mdat.MDAT_MLrev] <= mx & [mdat.zone] == 1 & ~isnan([mdat.MDAT_MLrev]));
ML2MWA(ind) = a1 * [mdat(ind).MDAT_MLrev] + a3;
ind = find([mdat.MDAT_MLrev] > mx & [mdat.zone] == 1 & ~isnan([mdat.MDAT_MLrev]));
ML2MWA(ind) = a1 * [mdat(ind).MDAT_MLrev] + a2 * ([mdat(ind).MDAT_MLrev] - mx) + a3;

% for pref ML
ind = find([mdat.MDAT_prefML] <= mx & [mdat.zone] == 1 & isnan([mdat.MDAT_MLrev]));
ML2MWA(ind) = a1 * [mdat(ind).MDAT_prefML] + a3;
ind = find([mdat.MDAT_prefML] > mx & [mdat.zone] == 1 & isnan([mdat.MDAT_MLrev]));
ML2MWA(ind) = a1 * [mdat(ind).MDAT_prefML] + a2 * ([mdat(ind).MDAT_prefML] - mx) + a3;

% note, changed max zone number to use SEA conversion for offshore events
disp('Converting ML to MW in eastern & south Australia...');
[a1,a2,a3,mx] = textread('F:\Catalogues\ML2MW\EA.ML-MW.coef.txt','%f%f%f%f','delimiter',',');
% for ML rev
ind = find([mdat.MDAT_MLrev] <= mx & [mdat.zone] >= 2 & [mdat.zone] <= 5 & ~isnan([mdat.MDAT_MLrev]));
ML2MWA(ind) = a1 * [mdat(ind).MDAT_MLrev] + a3;
ind = find([mdat.MDAT_MLrev] > mx & [mdat.zone] >= 2 & [mdat.zone] <= 5 & ~isnan([mdat.MDAT_MLrev]));
ML2MWA(ind) = a1 * [mdat(ind).MDAT_MLrev] + a2 * ([mdat(ind).MDAT_MLrev] - mx) + a3;

% for pref ML
ind = find([mdat.MDAT_prefML] <= mx & [mdat.zone] >= 2 & [mdat.zone] <= 5 & isnan([mdat.MDAT_MLrev]));
ML2MWA(ind) = a1 * [mdat(ind).MDAT_prefML] + a3;
ind = find([mdat.MDAT_prefML] > mx & [mdat.zone] >= 2 & [mdat.zone] <= 5 & isnan([mdat.MDAT_MLrev]));
ML2MWA(ind) = a1 * [mdat(ind).MDAT_prefML] + a2 * ([mdat(ind).MDAT_prefML] - mx) + a3;


%% Convert using Grunthal - out-dated

% disp('Converting ML to MW using Grunthal...');
% % ind = find([mdat.zone] == 4 & ~isnan([mdat.MDAT_MLrev]));
% ind = find(~isnan([mdat.MDAT_MLrev]));
% ML2MWG(ind) = 0.0376*[mdat(ind).MDAT_MLrev].^2 + 0.646*[mdat(ind).MDAT_MLrev] + 0.53;
% 
% % ind = find([mdat.zone] == 4 & isnan([mdat.MDAT_MLrev]));
% ind = find(isnan([mdat.MDAT_MLrev]));
% ML2MWG(ind) = 0.0376*[mdat(ind).MDAT_prefML].^2 + 0.646*[mdat(ind).MDAT_prefML] + 0.53;

%% Convert using Ghasemi (2017)

disp('Converting ML to MW using Ghasemi...');

% % set HG fixed mx reg coefs
% a1 = 0.66199378;
% a2 = 1.2156352;
% a3 = 1.07488336; % fixed coeff 
% mx = 4.5;
% my = a1 * mx + a2;
% 
% % py implementation
% % ans1 = (c[0] * x + c[1]) <= mx
% % yarea = c[0] * hx + c[1]
% % ans2 = (c[2] * (x-hx) + yarea) > mx
% 
% % for ML rev
% ind = find([mdat.MDAT_MLrev] <= mx & ~isnan([mdat.MDAT_MLrev]));
% ML2MWG(ind) = a1 * [mdat(ind).MDAT_MLrev] + a2;
% ind = find([mdat.MDAT_MLrev] > mx  & ~isnan([mdat.MDAT_MLrev]));
% ML2MWG(ind) = a3 * ([mdat(ind).MDAT_MLrev] - mx) + my;
% 
% % for pref ML
% ind = find([mdat.MDAT_prefML] <= mx & isnan([mdat.MDAT_MLrev]));
% ML2MWG(ind) = a1 * [mdat(ind).MDAT_prefML] + a2;
% ind = find([mdat.MDAT_prefML] > mx & isnan([mdat.MDAT_MLrev]));
% ML2MWG(ind) =  a3 * ([mdat(ind).MDAT_prefML] - mx) + my;

% use polynomial of simulated data from Ghasemi & Allen (2017)
a = 0.04160769;
b = 0.48058286;
c = 1.39485216;

for i = 1:length(mdat)
    % calculate MW for non-revised ML
	if isnan(mdat(i).MDAT_MLrev)
        ML2MWG(i) = a*mdat(i).MDAT_prefML^2 + b*mdat(i).MDAT_prefML + c;

    % calculate MW for non-revised ML
    else
         ML2MWG(i) = a*mdat(i).MDAT_MLrev^2 + b*mdat(i).MDAT_MLrev + c;
    end
end




%% set fields
for i = 1:length(mdat)
    mdat(i).MS2MW = MS2MW(i);
    mdat(i).mb2MW = mb2MW(i);
    mdat(i).ML2MWA = ML2MWA(i);
    mdat(i).ML2MWG = ML2MWG(i);
    mdat(i).prefFinalMW = prefFinalMW(i);
end

%% Set preferred MW

% conserve actual Mw measurements first
for i = 1:length(mdat)
    if ~isnan(mdat(i).MDAT_prefMW)
        mdat(i).prefFinalMW = mdat(i).MDAT_prefMW;
        mdat(i).prefFinalMWSrc = mdat(i).MDAT_prefMWSrc;
        mdat(i).MDAT_origMagType = 'MW';

% take larger of MS/mb >= 5.75
    elseif mdat(i).MS2MW > 5.75 | mdat(i).mb2MW > 5.75
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

%% Set preferred non-MW

% conserve actual Mw measurements first
for i = 1:length(mdat)
    if ~isnan(mdat(i).MDAT_prefMW)
        mdat(i).Mx_OrigML = mdat(i).MDAT_prefMW;
        mdat(i).Mx_RevML = mdat(i).MDAT_prefMW;
        mdat(i).Mx_RevMLSrc = mdat(i).MDAT_prefMWSrc;
        mdat(i).Mx_RevMLtype = 'MW';
         
% take larger of MS/mb >= 5.75
    elseif mdat(i).MDAT_prefMS >= 5.75 | mdat(i).MDAT_prefmb >= 5.75
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
% take Revised ML   
    elseif ~isnan(mdat(i).MDAT_MLrev)
        mdat(i).Mx_OrigML = mdat(i).MDAT_prefML;
        mdat(i).Mx_RevML = mdat(i).MDAT_MLrev;
        mdat(i).Mx_RevMLSrc = 'REV_ML';
        mdat(i).Mx_RevMLtype = 'REV_ML';
%         if ~isnan(mdat(i).MDAT_otherM)
%             mdat(i).Mx_OrigML = mdat(i).MDAT_otherM;
% %             mdat(i).Mx_RevMLSrc = mdat(i).MDAT_prefMLSrc;
%         end
% take Original ML   
    elseif ~isnan(mdat(i).MDAT_prefML)
        mdat(i).Mx_OrigML = mdat(i).MDAT_prefML;
        mdat(i).Mx_RevML = NaN;
        mdat(i).Mx_RevMLSrc = mdat(i).MDAT_prefMLSrc;
        mdat(i).Mx_RevMLtype = 'ML';
% take larger of MS/mb < 6.0        
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
% take Other Mag   
%     elseif ~isnan(mdat(i).MDAT_otherM)
%         mdat(i).Mx_OrigML = mdat(i).MDAT_otherM;
%         mdat(i).Mx_RevMLSrc = mdat(i).MDAT_prefMLSrc; 
%         mdat(i).Mx_RevMLtype = mdat(i).MDAT_otherMType;
    else
        mdat(i).Mx_OrigML = NaN;
        mdat(i).Mx_RevML = NaN;
        mdat(i).Mx_RevMLSrc = '';
        mdat(i).Mx_RevMLtype = '';
    end
end


% clear mdat;
disp('Saving mdat');
save mdat_mw_pref mdat;

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

% delind = find([mdat.MDAT_lon] > 160 & [mdat.MDAT_lon] > -4);
% mdat(delind) = [];
% delind = find([mdat.MDAT_lon] < 108 & [mdat.MDAT_lon] < -50);
% mdat(delind) = [];

%% write to file

header = 'DATESTR,DATENUM,TYPE,DEPENDENCE,LON,LAT,DEP,LOCSRC,PREFMW,PREFMWSRC,PREFMS,PREFMSSRC,PREFmb,PREFmbSRC,PREFML,PREFMLSRC,MLREG,REVML,MX_ORIGML,MX_TYPE,MX_REVML,MX_REVMLTYPE,MX_REVMLSRC,MS2MW,mb2MW,ML2MW,PREFMW,PREFMWSRC,COMM';
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
            num2str(mdat(i).Mx_RevML,'%0.2f'),',',mdat(i).Mx_RevMLtype,',', ...
            mdat(i).Mx_RevMLSrc,',', ...
            num2str(mdat(i).MS2MW,'%0.2f'),',',num2str(mdat(i).mb2MW,'%0.2f'),',', ...
            num2str(mdat(i).ML2MWG,'%0.2f'),',', ...
            num2str(mdat(i).prefFinalMW,'%0.2f'),',', ...
            mdat(i).prefFinalMWSrc,',',comms,char(10)];
    txt = [txt line];
end

dlmwrite(outfile,txt,'delimiter','','-append');

%% Make GMT mag diff file for pre-1990 events
clear txt;
ind = find([mdat.zone] ~= 4 & ~isnan([mdat.MDAT_prefML]) ...
      & ~isnan([mdat.MDAT_MLrev]));
mdiff = [mdat(ind).MDAT_MLrev] - [mdat(ind).MDAT_prefML];
dat = [[mdat(ind).MDAT_lon]' [mdat(ind).MDAT_lat]' mdiff' ...
      [mdat(ind).MDAT_MLrev]'/15];

dlmwrite('ML_diff.dat',dat,'delimiter','\t','precision','%0.3f');

%% plot TA vs HG mags

taml = [mdat.ML2MWA];
hgml = [mdat.ML2MWG];

figure(10);
plot(taml, hgml, 'b+')
hold on;
plot([1, 7],[1, 7],'k--')
xlabel('TA MW Conversion');
ylabel('HG (Fixed hinge) MW Conversion');

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
            




































