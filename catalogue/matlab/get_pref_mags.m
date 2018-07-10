% script to find preferred value of each magnitude type

if exist('mdat','var') ~= 1
    disp('Loading mdat');
    load mdat_mw.mat;
end

tabtxt = '';

%% loop thru all events
for i = 1:length(mdat)
    % find preferred ML
    if ~isnan(mdat(i).Allen_ML)
        mdat(i).MDAT_prefML = mdat(i).Allen_ML;
        mdat(i).MDAT_prefMLSrc = 'Allen (unpublished)';
        mdat(i).MDAT_origMLType = 'ML';
    elseif ~isnan(mdat(i).ANSN_ml)
        mdat(i).MDAT_prefML = mdat(i).ANSN_ml;
        mdat(i).MDAT_prefMLSrc = 'AUST';
        mdat(i).MDAT_origMLType = 'ML';
    elseif strcmp(deblank(mdat(i).GG_Mtype), 'ML')
        mdat(i).MDAT_prefML = mdat(i).GG_Mval;
        mdat(i).MDAT_prefMLSrc = mdat(i).MDAT_locsrc;
       mdat(i).MDAT_origMLType = 'ML';
    elseif strcmp(deblank(mdat(i).GG_Mtype), 'MLdz') || strcmp(deblank(mdat(i).GG_Mtype), 'MLzd')
        mdat(i).MDAT_prefML = mdat(i).GG_Mval;
        mdat(i).MDAT_prefMLSrc = mdat(i).MDAT_locsrc;
        mdat(i).MDAT_origMLType = 'ML';
    elseif strcmp(deblank(mdat(i).GG_Mtype), 'mL')
        mdat(i).MDAT_prefML = mdat(i).GG_Mval;
        mdat(i).MDAT_prefMLSrc = mdat(i).MDAT_locsrc;
        mdat(i).MDAT_origMLType = 'ML';
    elseif ~isnan(mdat(i).ISC_ml)
        mdat(i).MDAT_prefML = mdat(i).ISC_ml;
        mdat(i).MDAT_prefMLSrc = mdat(i).ISC_mlSRC;
    elseif strcmp(deblank(mdat(i).GG_Mtype), 'MP') || strcmp(deblank(mdat(i).GG_Mtype), 'MD') ...
           || strcmp(deblank(mdat(i).GG_Mtype), 'M?') 
        % assume equivalence with ML
        mdat(i).MDAT_prefML = mdat(i).GG_Mval;
        mdat(i).MDAT_prefMLSrc = mdat(i).MDAT_locsrc;
        mdat(i).MDAT_origMLType = mdat(i).GG_Mtype;
    else
        mdat(i).MDAT_prefML = NaN;
        mdat(i).MDAT_prefMLSrc = '';
    end
    
    % find preferred MW
    if ~isnan(mdat(i).altMW)
        mdat(i).MDAT_prefMW = mdat(i).altMW;
        mdat(i).MDAT_prefMWSrc = mdat(i).altMWsrc;
%     elseif ~isnan(mdat(i).ANSN_mw)
%         mdat(i).MDAT_prefMW = mdat(i).ANSN_mw;
%         mdat(i).MDAT_prefMWSrc = 'AUST';
    elseif strcmp(deblank(mdat(i).GG_Mtype), 'MW')
        mdat(i).MDAT_prefMW = mdat(i).GG_Mval;
        mdat(i).MDAT_prefMWSrc = mdat(i).MDAT_locsrc;
    elseif strcmp(deblank(mdat(i).GG_Mtype), 'Mw')
        mdat(i).MDAT_prefMW = mdat(i).GG_Mval;
        mdat(i).MDAT_prefMWSrc = mdat(i).MDAT_locsrc;
    elseif ~isnan(mdat(i).ISC_mw)
        mdat(i).MDAT_prefMW = mdat(i).ISC_mw;
        mdat(i).MDAT_prefMWSrc = mdat(i).ISC_mwSRC;
%     elseif ~isnan(mdat(i).ANSN_mwp)
%         mdat(i).MDAT_prefMW = mdat(i).ANSN_mwp;
%         mdat(i).MDAT_prefMWSrc = 'AUST';
    else
        mdat(i).MDAT_prefMW = NaN;
        mdat(i).MDAT_prefMWSrc = '';
    end

    % find preferred MS
    if strcmp(deblank(mdat(i).GG_Mtype), 'MS') || strcmp(deblank(mdat(i).GG_Mtype), 'Ms')
        mdat(i).MDAT_prefMS = mdat(i).GG_Mval;
        mdat(i).MDAT_prefMSSrc = mdat(i).MDAT_locsrc;
    elseif ~isnan(mdat(i).ISC_ms)
        mdat(i).MDAT_prefMS = mdat(i).ISC_ms;
        mdat(i).MDAT_prefMSSrc = mdat(i).ISC_msSRC;
    elseif ~isnan(mdat(i).ANSN_ms)
        mdat(i).MDAT_prefMS = mdat(i).ANSN_ms;
        mdat(i).MDAT_prefMSSrc = 'AUST';
    else
        mdat(i).MDAT_prefMS = NaN;
        mdat(i).MDAT_prefMSSrc = '';
    end
    
    % find preferred mb
    if strcmp(deblank(mdat(i).GG_Mtype), 'mb')
        mdat(i).MDAT_prefmb = mdat(i).GG_Mval;
        mdat(i).MDAT_prefmbSrc = mdat(i).MDAT_locsrc;
    elseif ~isnan(mdat(i).ANSN_mb)
        mdat(i).MDAT_prefmb = mdat(i).ANSN_mb;
        mdat(i).MDAT_prefmbSrc = 'AUST';
    elseif ~isnan(mdat(i).ISC_mb)
        mdat(i).MDAT_prefmb = mdat(i).ISC_mb;
        mdat(i).MDAT_prefmbSrc = mdat(i).ISC_mbSRC;
    else
        mdat(i).MDAT_prefmb = NaN;
        mdat(i).MDAT_prefmbSrc = '';
    end
    
    % make text for Appendix of NSHA18-Cat
    if ~isnan(mdat(i).Allen_ML)
        line = [datestr(mdat(i).MDAT_dateNum, 31), ',', num2str(mdat(i).MDAT_lon), ',', num2str(mdat(i).MDAT_lat), ',', ...
                num2str(mdat(i).GG_Mval), ',', num2str(mdat(i).GG_Mtype), ',', num2str(mdat(i).ANSN_ml), ',', ...
                num2str(mdat(i).Allen_ML), ',', num2str(mdat(i).MDAT_prefMW)];
                
        tabtxt = [tabtxt line char(10)];
        
    end
end

% write alt ML file
header = 'DATESTR,LON,LAT,GGML,GGTYPE,AUSTML,ALLENML,PREFMW';
disp('writing to file')
dlmwrite('ml_comparisons.csv',header,'delimiter','');

dlmwrite('ml_comparisons.csv',tabtxt,'delimiter','','-append');

% save mat file
save mdat_pref_mag_types mdat