% script to find preferred value of each magnitude type

if exist('mdat','var') ~= 1
    disp('Loading mdat');
    load mdat_mw.mat;
end

%% loop thru all events
for i = 1:length(mdat)
    % find preferred ML
    if ~isnan(mdat(i).Allen_ML)
        mdat(i).MDAT_prefML = mdat(i).Allen_ML;
        mdat(i).MDAT_prefMLSrc = 'Allen (unpublished)';
    elseif ~isnan(mdat(i).ANSN_ml)
        mdat(i).MDAT_prefML = mdat(i).ANSN_ml;
        mdat(i).MDAT_prefMLSrc = 'AUST';
    elseif strcmp(mdat(i).GG_Mtype, 'ML')
        mdat(i).MDAT_prefML = mdat(i).GG_Mval;
        mdat(i).MDAT_prefMLSrc = mdat(i).MDAT_locsrc;
    elseif strcmp(mdat(i).GG_Mtype, 'MLdz')
        mdat(i).MDAT_prefML = mdat(i).GG_Mval;
        mdat(i).MDAT_prefMLSrc = mdat(i).MDAT_locsrc;
    elseif strcmp(mdat(i).GG_Mtype, 'MP') || strcmp(mdat(i).GG_Mtype, 'MD') ...
           || strcmp(mdat(i).GG_Mtype, 'M?') 
        % assume equivalence with ML
        mdat(i).MDAT_prefML = mdat(i).GG_Mval;
        mdat(i).MDAT_prefMLSrc = mdat(i).MDAT_locsrc;
    else
        mdat(i).MDAT_prefML = NaN;
        mdat(i).MDAT_prefMLSrc = '';
    end
    
    % find preferred MW
    if ~isnan(mdat(i).altMW)
        mdat(i).MDAT_prefMW = mdat(i).altMW;
        mdat(i).MDAT_prefMWSrc = mdat(i).altMWsrc;
    elseif ~isnan(mdat(i).ANSN_mw)
        mdat(i).MDAT_prefMW = mdat(i).ANSN_mw;
        mdat(i).MDAT_prefMWSrc = 'AUST';
    elseif strcmp(mdat(i).GG_Mtype, 'MW')
        mdat(i).MDAT_prefMW = mdat(i).GG_Mval;
        mdat(i).MDAT_prefMWSrc = mdat(i).MDAT_locsrc;
    elseif ~isnan(mdat(i).ANSN_mwp)
        mdat(i).MDAT_prefMW = mdat(i).ANSN_mwp;
        mdat(i).MDAT_prefMWSrc = 'AUST';
    else
        mdat(i).MDAT_prefMW = NaN;
        mdat(i).MDAT_prefMWSrc = '';
    end

    % find preferred MS
    if strcmp(mdat(i).GG_Mtype, 'MS')
        mdat(i).MDAT_prefMS = mdat(i).GG_Mval;
        mdat(i).MDAT_prefMSSrc = mdat(i).MDAT_locsrc;
    elseif ~isnan(mdat(i).ANSN_ms)
        mdat(i).MDAT_prefMS = mdat(i).ANSN_ms;
        mdat(i).MDAT_prefMSSrc = 'AUST';
     else
        mdat(i).MDAT_prefMS = NaN;
        mdat(i).MDAT_prefMSSrc = '';
    end
    
    % find preferred mb
    if strcmp(mdat(i).GG_Mtype, 'mb')
        mdat(i).MDAT_prefmb = mdat(i).GG_Mval;
        mdat(i).MDAT_prefmbSrc = mdat(i).MDAT_locsrc;
    elseif ~isnan(mdat(i).ANSN_mb)
        mdat(i).MDAT_prefmb = mdat(i).ANSN_mb;
        mdat(i).MDAT_prefmbSrc = 'AUST';
     else
        mdat(i).MDAT_prefmb = NaN;
        mdat(i).MDAT_prefmbSrc = '';
    end
    
end

% save mat file
save mdat_pref_mag_types mdat