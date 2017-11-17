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
outfile = '..\..\data\AUSTCAT.MP.V0.12.csv';

% load data

if exist('mdat_pref','var') ~= 1
    disp('Loading mdat_mw_pref 12');
    load ..\append_mw\mdat_no_mw_pref12.mat;
end

% if exist('mdat','var') ~= 1
%     disp('Loading mdat');
%     load '..\Merge MC and QUAKES\mdat.mat';
% end

MS2MW = ones(size(mdat_pref)) * NaN;
mb2MW = ones(size(mdat_pref)) * NaN;
ML2MWA = ones(size(mdat_pref)) * NaN;
ML2MWG = ones(size(mdat_pref)) * NaN;
prefFinalMW = ones(size(mdat_pref)) * NaN;

%% Convert MS to MW using Di Giacomo et al (2015) for zone 4

% conversion for shallow earthquakes (h < 70 km)
disp('Converting MS to MW...');
% for 3.0 <= MS <= 6.1
ind1 = find([mdat_pref.MDAT_dep] <= 70 & [mdat_pref.MDAT_prefMS] >= 3.0 ...
           & [mdat_pref.MDAT_prefMS] <= 6.47 & [mdat_pref.zone] == 4);
% include those events with NaN depth       
ind2 = find(isnan([mdat_pref.MDAT_dep]) & [mdat_pref.MDAT_prefMS] >= 3.0 ...
           & [mdat_pref.MDAT_prefMS] <= 6.47 & [mdat_pref.zone] == 4);
% include those events h > 70 and no mb   
ind3 = find([mdat_pref.MDAT_dep] > 70  & [mdat_pref.MDAT_prefMS] >= 3.0 ...
           & [mdat_pref.MDAT_prefMS] <= 6.47 & isnan([mdat_pref.MDAT_prefmb]) ...
           & [mdat_pref.zone] == 4);

ind = [ind1, ind2, ind3];       
%MS2MW(ind) = 0.67 * [mdat_pref(ind).MDAT_prefMS] + 2.13;
MS2MW(ind) = mdat_pref(ind).MDAT_prefMS;

% for 6.47 < MS <= 8.2! (Ignore upper lim for now)
ind1 = find([mdat_pref.MDAT_dep] <= 70 & [mdat_pref.MDAT_prefMS] > 6.47 ...
             & [mdat_pref.zone] == 4);
% include those events with NaN depth       
ind2 = find(isnan([mdat_pref.MDAT_dep]) & [mdat_pref.MDAT_prefMS] > 6.47 ...
             & [mdat_pref.zone] == 4);
% include those events h > 70 and no mb   
ind3 = find([mdat_pref.MDAT_dep] > 70  & [mdat_pref.MDAT_prefMS] > 6.47 ...
            & isnan([mdat_pref.MDAT_prefmb]) & [mdat_pref.zone] == 4);

ind = [ind1, ind2, ind3];       
%MS2MW(ind) = 1.10 * [mdat_pref(ind).MDAT_prefMS] - 0.67;
MS2MW(ind) = mdat_pref(ind).MDAT_prefMS;

%% Convert MS to MW using Ghasemi (2017) for Aust events
c1 = 0.84896727404297323;
c2 = 1.0509630268292971;

% for 3.0 <= MS <= 7.0
ind = find(~isnan([mdat_pref.MDAT_prefMS]) & [mdat_pref.zone] ~= 4);
%MS2MW(ind) = c1 * [mdat_pref(ind).MDAT_prefMS] + c2;
MS2MW(ind) = mdat_pref(ind).MDAT_prefMS;

%% Convert mb to MW using Di Giacomo et al (2015) for zone 4

disp('Converting mb to MW...');
% for 3.5 <= mb <= 6.2
ind = find([mdat_pref.MDAT_prefmb] >= 4 & [mdat_pref.MDAT_prefmb] <= 6.2 ...
           & [mdat_pref.zone] == 4);
%mb2MW(ind) = 1.38 * [mdat_pref(ind).MDAT_prefmb] - 1.79;
mb2MW(ind) = mdat_pref(ind).MDAT_prefmb;

%% Convert mb to MW using Allen (2012) for events below latitude -13 deg - out-dated!
% c1 = 0.7362;
% c2 = 0.7374;
% c3 = 0.9707;
% mx = 5.0905;

% ind = find([mdat_pref.MDAT_prefmb] >= 3.0 & [mdat_pref.MDAT_prefmb] < mx ...
%        & [mdat_pref.zone] == 4);
% mb2MW(ind) = c1 * [mdat_pref(ind).MDAT_prefmb] + c3;
% ind = find([mdat_pref.MDAT_prefmb] >= mx & [mdat_pref.MDAT_prefmb] <= 6.5 ...
%             & [mdat_pref.zone] == 4);
% mb2MW(ind) = c1 * [mdat_pref(ind).MDAT_prefmb] ...
%              + c2 * ([mdat_pref(ind).MDAT_prefmb] - mx) + c3;
         
%% Convert mb to MW using Ghasemi (2017) for Aust events
c1 = 1.1438907424442797;
c2 = -0.87192285009579173;

% for 3.0 <= MS <= 7.0
ind = find(~isnan([mdat_pref.MDAT_prefmb]) & [mdat_pref.zone] ~= 4);
%mb2MW(ind) = c1 * [mdat_pref(ind).MDAT_prefmb] + c2;
mb2MW(ind) = mdat_pref(ind).MDAT_prefmb;

%% Convert ML to MW using Allen conversions - out-dated, but preserve in catalogue
% mx = 4.2;

disp('Converting ML to MW in CWA...');
[a1,a2,a3,mx] = textread('F:\Catalogues\ML2MW\WA.ML-MW.coef.txt','%f%f%f%f','delimiter',',');
% for ML rev
ind = find([mdat_pref.MDAT_MLrev] <= mx & [mdat_pref.zone] == 1 & ~isnan([mdat_pref.MDAT_MLrev]));
ML2MWA(ind) = a1 * [mdat_pref(ind).MDAT_MLrev] + a3;
ind = find([mdat_pref.MDAT_MLrev] > mx & [mdat_pref.zone] == 1 & ~isnan([mdat_pref.MDAT_MLrev]));
ML2MWA(ind) = a1 * [mdat_pref(ind).MDAT_MLrev] + a2 * ([mdat_pref(ind).MDAT_MLrev] - mx) + a3;

% for pref ML
ind = find([mdat_pref.MDAT_prefML] <= mx & [mdat_pref.zone] == 1 & isnan([mdat_pref.MDAT_MLrev]));
ML2MWA(ind) = a1 * [mdat_pref(ind).MDAT_prefML] + a3;
ind = find([mdat_pref.MDAT_prefML] > mx & [mdat_pref.zone] == 1 & isnan([mdat_pref.MDAT_MLrev]));
ML2MWA(ind) = a1 * [mdat_pref(ind).MDAT_prefML] + a2 * ([mdat_pref(ind).MDAT_prefML] - mx) + a3;

% note, changed max zone number to use SEA conversion for offshore events
disp('Converting ML to MW in eastern & south Australia...');
[a1,a2,a3,mx] = textread('F:\Catalogues\ML2MW\EA.ML-MW.coef.txt','%f%f%f%f','delimiter',',');
% for ML rev
ind = find([mdat_pref.MDAT_MLrev] <= mx & [mdat_pref.zone] >= 2 & [mdat_pref.zone] <= 5 & ~isnan([mdat_pref.MDAT_MLrev]));
ML2MWA(ind) = a1 * [mdat_pref(ind).MDAT_MLrev] + a3;
ind = find([mdat_pref.MDAT_MLrev] > mx & [mdat_pref.zone] >= 2 & [mdat_pref.zone] <= 5 & ~isnan([mdat_pref.MDAT_MLrev]));
ML2MWA(ind) = a1 * [mdat_pref(ind).MDAT_MLrev] + a2 * ([mdat_pref(ind).MDAT_MLrev] - mx) + a3;

% for pref ML
ind = find([mdat_pref.MDAT_prefML] <= mx & [mdat_pref.zone] >= 2 & [mdat_pref.zone] <= 5 & isnan([mdat_pref.MDAT_MLrev]));
ML2MWA(ind) = a1 * [mdat_pref(ind).MDAT_prefML] + a3;
ind = find([mdat_pref.MDAT_prefML] > mx & [mdat_pref.zone] >= 2 & [mdat_pref.zone] <= 5 & isnan([mdat_pref.MDAT_MLrev]));
ML2MWA(ind) = a1 * [mdat_pref(ind).MDAT_prefML] + a2 * ([mdat_pref(ind).MDAT_prefML] - mx) + a3;


%% Convert using Grunthal - out-dated

% disp('Converting ML to MW using Grunthal...');
% % ind = find([mdat_pref.zone] == 4 & ~isnan([mdat_pref.MDAT_MLrev]));
% ind = find(~isnan([mdat_pref.MDAT_MLrev]));
% ML2MWG(ind) = 0.0376*[mdat_pref(ind).MDAT_MLrev].^2 + 0.646*[mdat_pref(ind).MDAT_MLrev] + 0.53;
% 
% % ind = find([mdat_pref.zone] == 4 & isnan([mdat_pref.MDAT_MLrev]));
% ind = find(isnan([mdat_pref.MDAT_MLrev]));
% ML2MWG(ind) = 0.0376*[mdat_pref(ind).MDAT_prefML].^2 + 0.646*[mdat_pref(ind).MDAT_prefML] + 0.53;

%% Convert using Ghasemi (2017)

disp('Converting ML to MW using Ghasemi...');

% set HG fixed mx reg coefs
a1 = 0.66199378;
a2 = 1.2156352;
a3 = 1.07488336; % fixed coeff 
mx = 4.5;
my = a1 * mx + a2;

% py implementation
% ans1 = (c[0] * x + c[1]) <= mx
% yarea = c[0] * hx + c[1]
% ans2 = (c[2] * (x-hx) + yarea) > mx

% for ML rev
ind = find([mdat_pref.MDAT_MLrev] <= mx & ~isnan([mdat_pref.MDAT_MLrev]));
%ML2MWG(ind) = a1 * [mdat_pref(ind).MDAT_MLrev] + a2;
ML2MWG(ind) =  mdat_pref(ind).MDAT_prefML;

ind = find([mdat_pref.MDAT_MLrev] > mx  & ~isnan([mdat_pref.MDAT_MLrev]));
%ML2MWG(ind) = a3 * ([mdat_pref(ind).MDAT_MLrev] - mx) + my;
ML2MWG(ind) =  mdat_pref(ind).MDAT_prefML;

% for pref ML
ind = find([mdat_pref.MDAT_prefML] <= mx & isnan([mdat_pref.MDAT_MLrev]));
%ML2MWG(ind) = a1 * [mdat_pref(ind).MDAT_prefML] + a2;
ML2MWG(ind) =  mdat_pref(ind).MDAT_prefML;

ind = find([mdat_pref.MDAT_prefML] > mx & isnan([mdat_pref.MDAT_MLrev]));
%ML2MWG(ind) =  a3 * ([mdat_pref(ind).MDAT_prefML] - mx) + my;
ML2MWG(ind) =  mdat_pref(ind).MDAT_prefML;

%% set fields
for i = 1:length(mdat_pref)
    mdat_pref(i).MS2MW = MS2MW(i);
    mdat_pref(i).mb2MW = mb2MW(i);
    mdat_pref(i).ML2MWA = ML2MWA(i);
    mdat_pref(i).ML2MWG = ML2MWG(i);
    mdat_pref(i).prefFinalMW = prefFinalMW(i);
end

%% Set preferred MW

% conserve actual Mw measurements first
for i = 1:length(mdat_pref)
    if ~isnan(mdat_pref(i).MDAT_prefMW)
        mdat_pref(i).prefFinalMW = mdat_pref(i).MDAT_prefMW;
        mdat_pref(i).prefFinalMWSrc = mdat_pref(i).MDAT_prefMWSrc;

% take larger of MS/mb >= 5.75
    elseif mdat_pref(i).MS2MW > 5.75 | mdat_pref(i).mb2MW > 5.75
        maxM = max([mdat_pref(i).MDAT_prefMS mdat_pref(i).MDAT_prefmb]); % deliberately use orig mag here
        if mdat_pref(i).MDAT_prefMS == maxM
            mdat_pref(i).prefFinalMW = mdat_pref(i).MS2MW;
            mdat_pref(i).prefFinalMWSrc = 'MS2MW';
        elseif mdat_pref(i).MDAT_prefmb == maxM
            mdat_pref(i).prefFinalMW = mdat_pref(i).mb2MW;
            mdat_pref(i).prefFinalMWSrc = 'mb2MW';
        end
        
% take ML-MW
    elseif ~isnan(mdat_pref(i).ML2MWG)
        mdat_pref(i).prefFinalMW = mdat_pref(i).ML2MWG;
        mdat_pref(i).prefFinalMWSrc = 'ML2MWG';

% take larger of MS/mb < 6.0        
    elseif ~isnan(mdat_pref(i).MS2MW) | ~isnan(mdat_pref(i).mb2MW)
        maxM = max([mdat_pref(i).MDAT_prefMS mdat_pref(i).MDAT_prefmb]); % deliberately use orig mag here
        if mdat_pref(i).MDAT_prefMS == maxM
            mdat_pref(i).prefFinalMW = mdat_pref(i).MS2MW;
            mdat_pref(i).prefFinalMWSrc = 'MS2MW';
        elseif mdat_pref(i).MDAT_prefmb == maxM
            mdat_pref(i).prefFinalMW = mdat_pref(i).mb2MW;
            mdat_pref(i).prefFinalMWSrc = 'mb2MW';
        else
            mdat_pref(i).prefFinalMW = NaN;
            mdat_pref(i).prefFinalMWSrc = '';
        end
    end
    
% ADD GG DEPENDENCE
    if isempty(mdat_pref(i).GG_sourceType)
        mdat_pref(i).GG_sourceType = '';
    end
    if isempty(mdat_pref(i).GG_dependence)
        mdat_pref(i).GG_dependence = '';
    end
%     mdat_pref(i).GG_dependence = mdat(i).GG_dependence;
end

%% Set preferred non-MW

% conserve actual Mw measurements first
for i = 1:length(mdat_pref)
    if ~isnan(mdat_pref(i).MDAT_prefMW)
        mdat_pref(i).Mx_OrigML = mdat_pref(i).MDAT_prefMW;
        mdat_pref(i).Mx_RevML = mdat_pref(i).MDAT_prefMW;
        mdat_pref(i).Mx_RevMLSrc = mdat_pref(i).MDAT_prefMWSrc;
        mdat_pref(i).Mx_RevMLtype = 'MW';
         
% take larger of MS/mb >= 6.0
    elseif mdat_pref(i).MDAT_prefMS >= 6.0 | mdat_pref(i).MDAT_prefmb >= 6.0
        maxM = max([mdat_pref(i).MDAT_prefMS mdat_pref(i).MDAT_prefmb]);
        if mdat_pref(i).MDAT_prefMS == maxM
            mdat_pref(i).Mx_OrigML = mdat_pref(i).MDAT_prefMS;
            mdat_pref(i).Mx_RevML = mdat_pref(i).MDAT_prefMS;
            mdat_pref(i).Mx_RevMLSrc = mdat_pref(i).MDAT_prefMSSrc;
            mdat_pref(i).Mx_RevMLtype = 'MS';
        elseif mdat_pref(i).MDAT_prefmb == maxM
            mdat_pref(i).Mx_OrigML = mdat_pref(i).MDAT_prefmb;
            mdat_pref(i).Mx_RevML = mdat_pref(i).MDAT_prefmb;
            mdat_pref(i).Mx_RevMLSrc = mdat_pref(i).MDAT_prefmbSrc;
            mdat_pref(i).Mx_RevMLtype = 'mb';
        end
% take Revised ML   
    elseif ~isnan(mdat_pref(i).MDAT_MLrev)
        mdat_pref(i).Mx_OrigML = mdat_pref(i).MDAT_prefML;
        mdat_pref(i).Mx_RevML = mdat_pref(i).MDAT_MLrev;
        mdat_pref(i).Mx_RevMLSrc = 'REV_ML';
        mdat_pref(i).Mx_RevMLtype = 'REV_ML';
        if ~isnan(mdat_pref(i).MDAT_otherM)
            mdat_pref(i).Mx_OrigML = mdat_pref(i).MDAT_otherM;
%             mdat_pref(i).Mx_RevMLSrc = mdat_pref(i).MDAT_prefMLSrc;
        end
% take Original ML   
    elseif ~isnan(mdat_pref(i).MDAT_prefML)
        mdat_pref(i).Mx_OrigML = mdat_pref(i).MDAT_prefML;
        mdat_pref(i).Mx_RevML = mdat_pref(i).MDAT_prefML;
        mdat_pref(i).Mx_RevMLSrc = mdat_pref(i).MDAT_prefMLSrc;
        mdat_pref(i).Mx_RevMLtype = 'ML';
% take larger of MS/mb < 6.0        
    elseif ~isnan(mdat_pref(i).MDAT_prefMS) | ~isnan(mdat_pref(i).MDAT_prefmb)
        maxM = max([mdat_pref(i).MDAT_prefMS mdat_pref(i).MDAT_prefmb]);
        if mdat_pref(i).MDAT_prefMS == maxM
            mdat_pref(i).Mx_OrigML = mdat_pref(i).MDAT_prefMS;
            mdat_pref(i).Mx_RevML = mdat_pref(i).MDAT_prefMS;
            mdat_pref(i).Mx_RevMLSrc = mdat_pref(i).MDAT_prefMSSrc;
            mdat_pref(i).Mx_RevMLtype = 'MS';
        elseif mdat_pref(i).MDAT_prefmb == maxM
            mdat_pref(i).Mx_OrigML = mdat_pref(i).MDAT_prefmb;
            mdat_pref(i).Mx_RevML = mdat_pref(i).MDAT_prefmb;
            mdat_pref(i).Mx_RevMLSrc = mdat_pref(i).MDAT_prefmbSrc;
            mdat_pref(i).Mx_RevMLtype = 'mb';
        end
% take Other Mag   
    elseif ~isnan(mdat_pref(i).MDAT_otherM)
        mdat_pref(i).Mx_OrigML = mdat_pref(i).MDAT_otherM;
        mdat_pref(i).Mx_RevMLSrc = mdat_pref(i).MDAT_prefMLSrc; 
        mdat_pref(i).Mx_RevMLtype = mdat_pref(i).MDAT_otherMType;
    else
        mdat_pref(i).Mx_OrigML = NaN;
        mdat_pref(i).Mx_RevML = NaN;
        mdat_pref(i).Mx_RevMLSrc = '';
        mdat_pref(i).Mx_RevMLtype = '';
    end
end


% clear mdat;
disp('Saving mdat_pref');
save mdat_mw_pref12_final mdat_pref;

%% remove unecessary events
% delind = find(strcmp({mdat_pref.MDAT_prefmbSrc},'ISC') & strcmp({mdat_pref.MDAT_prefMLSrc},'AUST') == 0 ...
%          & strcmp({mdat_pref.MDAT_prefMLSrc},'MEL') == 0 & strcmp({mdat_pref.MDAT_prefMLSrc},'GG') == 0 ...
%          & strcmp({mdat_pref.MDAT_prefMLSrc},'ADE') == 0 & strcmp({mdat_pref.MDAT_prefMLSrc},'MGO') == 0 ...
%          & strcmp({mdat_pref.MDAT_prefMLSrc},'MUN') == 0 & strcmp({mdat_pref.MDAT_prefMLSrc},'AGSO') == 0 ...
%          & [mdat_pref.MDAT_prefmb] < 5.0);
% mdat_pref(delind) = [];
% delind = find(strcmp({mdat_pref.MDAT_prefmbSrc},'IDC') & strcmp({mdat_pref.MDAT_prefMLSrc},'AUST') == 0 ...
%          & strcmp({mdat_pref.MDAT_prefMLSrc},'MEL') == 0 & strcmp({mdat_pref.MDAT_prefMLSrc},'GG') == 0 ...
%          & strcmp({mdat_pref.MDAT_prefMLSrc},'ADE') == 0 & strcmp({mdat_pref.MDAT_prefMLSrc},'MGO') == 0 ...
%          & strcmp({mdat_pref.MDAT_prefMLSrc},'MUN') == 0 & strcmp({mdat_pref.MDAT_prefMLSrc},'AGSO') == 0 ...
%          & [mdat_pref.MDAT_prefmb] < 5.0);
% mdat_pref(delind) = [];
% delind = find(strcmp({mdat_pref.MDAT_prefmbSrc},'ISC') & strcmp({mdat_pref.MDAT_prefMLSrc},'AUST') == 0 ...
%          & strcmp({mdat_pref.MDAT_prefMLSrc},'MEL') == 0 & strcmp({mdat_pref.MDAT_prefMLSrc},'GG') == 0 ...
%          & strcmp({mdat_pref.MDAT_prefMLSrc},'ADE') == 0 & strcmp({mdat_pref.MDAT_prefMLSrc},'MGO') == 0 ...
%          & strcmp({mdat_pref.MDAT_prefMLSrc},'MUN') == 0 & strcmp({mdat_pref.MDAT_prefMLSrc},'AGSO') == 0 ...
%          & [mdat_pref.MDAT_prefmb] < 5.0);
% mdat_pref(delind) = [];
% delind = find(strcmp({mdat_pref.MDAT_prefmbSrc},'IDC') & strcmp({mdat_pref.MDAT_prefMLSrc},'AUST') == 0 ...
%          & strcmp({mdat_pref.MDAT_prefMLSrc},'MEL') == 0 & strcmp({mdat_pref.MDAT_prefMLSrc},'GG') == 0 ...
%          & strcmp({mdat_pref.MDAT_prefMLSrc},'ADE') == 0 & strcmp({mdat_pref.MDAT_prefMLSrc},'MGO') == 0 ...
%          & strcmp({mdat_pref.MDAT_prefMLSrc},'MUN') == 0 & strcmp({mdat_pref.MDAT_prefMLSrc},'AGSO') == 0 ...
%          & [mdat_pref.MDAT_prefmb] < 5.0);
% mdat_pref(delind) = [];
% delind = find(strcmp({mdat_pref.MDAT_prefmbSrc},'DJA') & strcmp({mdat_pref.MDAT_prefMLSrc},'AUST') == 0 ...
%          & strcmp({mdat_pref.MDAT_prefMLSrc},'MEL') == 0 & strcmp({mdat_pref.MDAT_prefMLSrc},'GG') == 0 ...
%          & strcmp({mdat_pref.MDAT_prefMLSrc},'ADE') == 0 & strcmp({mdat_pref.MDAT_prefMLSrc},'MGO') == 0 ...
%          & strcmp({mdat_pref.MDAT_prefMLSrc},'MUN') == 0 & strcmp({mdat_pref.MDAT_prefMLSrc},'AGSO') == 0 ...
%          & [mdat_pref.MDAT_prefmb] < 5.0);
% mdat_pref(delind) = [];
% delind = find(strcmp({mdat_pref.MDAT_prefMLSrc},'DJA'));
% mdat_pref(delind) = [];

delind = find([mdat_pref.MDAT_lon] > 160 & [mdat_pref.MDAT_lon] > -4);
mdat_pref(delind) = [];
delind = find([mdat_pref.MDAT_lon] < 108 & [mdat_pref.MDAT_lon] < -50);
mdat_pref(delind) = [];

%% write to file

header = 'DATESTR,DATENUM,TYPE,DEPENDENCE,LON,LAT,DEP,LOCSRC,PREFMW,PREFMWSRC,PREFMS,PREFMSSRC,PREFmb,PREFmbSRC,PREFML,PREFMLSRC,REVML,OTHERM,OTHERMTYPE,OTHERMSRC,MX_ORIGML,MX_REVML,MX_REVMLTYPE,MX_REVMLSRC,MS2MW,mb2MW,ML2MWA,ML2MWG,PREFMW,PREFMWSRC,COMM';
disp('writing to file')
dlmwrite(outfile,header,'delimiter','');
txt = [];
for i = 1:length(mdat_pref)
    if isempty(mdat_pref(i).GG_comm)
        comms = ['"',mdat_pref(i).ISC_evname,'"'];
    else
        comms = ['"',mdat_pref(i).GG_comm,'"'];
    end
    
    line = [datestr(mdat_pref(i).MDAT_dateNum,31),',',num2str(mdat_pref(i).MDAT_dateNum),',', ...
            mdat_pref(i).GG_sourceType,',',mdat_pref(i).GG_dependence,',', ...
            num2str(mdat_pref(i).MDAT_lon),',', ...
            num2str(mdat_pref(i).MDAT_lat),',',num2str(mdat_pref(i).MDAT_dep),',', ...
            mdat_pref(i).MDAT_locsrc,',', ...
            num2str(mdat_pref(i).MDAT_prefMW),',',mdat_pref(i).MDAT_prefMWSrc,',', ...
            num2str(mdat_pref(i).MDAT_prefMS),',',mdat_pref(i).MDAT_prefMSSrc,',', ...
            num2str(mdat_pref(i).MDAT_prefmb),',',mdat_pref(i).MDAT_prefmbSrc,',', ...
            num2str(mdat_pref(i).MDAT_prefML),',',mdat_pref(i).MDAT_prefMLSrc,',', ...
            num2str(mdat_pref(i).MDAT_MLrev,'%0.1f'),',', ...
            num2str(mdat_pref(i).MDAT_otherM),',',mdat_pref(i).MDAT_otherMType,',', ...
            mdat_pref(i).MDAT_otherMSrc,',',num2str(mdat_pref(i).Mx_OrigML,'%0.2f'),',', ...
            num2str(mdat_pref(i).Mx_RevML,'%0.2f'),',',mdat_pref(i).Mx_RevMLtype,',', ...
            mdat_pref(i).Mx_RevMLSrc,',', ...
            num2str(mdat_pref(i).MS2MW,'%0.2f'),',',num2str(mdat_pref(i).mb2MW,'%0.2f'),',', ...
            num2str(mdat_pref(i).ML2MWA,'%0.2f'),',',num2str(mdat_pref(i).ML2MWG,'%0.2f'),',', ...
            num2str(mdat_pref(i).prefFinalMW,'%0.2f'),',', ...
            mdat_pref(i).prefFinalMWSrc,',',comms,char(10)];
    txt = [txt line];
end

dlmwrite(outfile,txt,'delimiter','','-append');

%% Make GMT mag diff file for pre-1990 events
clear txt;
ind = find([mdat_pref.zone] ~= 4 & ~isnan([mdat_pref.MDAT_prefML]) ...
      & ~isnan([mdat_pref.MDAT_MLrev]));
mdiff = [mdat_pref(ind).MDAT_MLrev] - [mdat_pref(ind).MDAT_prefML];
dat = [[mdat_pref(ind).MDAT_lon]' [mdat_pref(ind).MDAT_lat]' mdiff' ...
      [mdat_pref(ind).MDAT_MLrev]'/15];

dlmwrite('ML_diff.dat',dat,'delimiter','\t','precision','%0.3f');

%% plot TA vs HG mags

taml = [mdat_pref.ML2MWA];
hgml = [mdat_pref.ML2MWG];

figure(10);
plot(taml, hgml, 'b+')
hold on;
plot([1, 7],[1, 7],'k--')
xlabel('TA MW Conversion');
ylabel('HG (Fixed hinge) MW Conversion');

%% plot histograms of ML difference

figure(1);
ind = find([mdat_pref.MDAT_dateNum] > datenum(1940,1,1) ...
           & [mdat_pref.MDAT_dateNum] < datenum(1990,1,1));
mldiff = [mdat_pref(ind).MDAT_prefML] - [mdat_pref(ind).MDAT_MLrev];
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
            
ind = find([mdat_pref.MDAT_dateNum] > datenum(1940,1,1) ...
           & [mdat_pref.MDAT_dateNum] < datenum(1990,1,1) ...
           & [mdat_pref.MDAT_prefML] >= 4.0);
mldiff = [mdat_pref(ind).MDAT_prefML] - [mdat_pref(ind).MDAT_MLrev];
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
            




































