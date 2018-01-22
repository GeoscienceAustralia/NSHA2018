% script to merge ANSN's earthquake catalogue to mdat

%% parse ANSN_ catalogue
disp('Parsing ANSN Catalogue...');
ANSN_catfile = fullfile('..','data','ANSN_cat_m_ge_32.csv');
[prefmag, utcdate, utctime, localdate, localtime, lat, lon, magtype, place, ...
 dep, soln, mb, ml, ms, mwp, mw, evid, oid] = ...
 textread(ANSN_catfile, '%f%s%s%s%s %f%f%s%s %f%s%f%f%f%f%f%f%f', ...
          'headerlines',1,'delimiter',',','emptyvalue',NaN);
      
for i = 1:length(prefmag)
    % get datetime
    dateSplit = str2double(strsplit(utcdate{i},'/'));
    timeSplit = str2double(strsplit(utctime{i},':'));
    %dateSplit = str2double(strsplit('/',utcdate{i})); % for mac
    %timeSplit = str2double(strsplit(':',utctime{i})); % for mac
    ANSN_dat(i).dateNum = datenum(dateSplit(3), dateSplit(2), dateSplit(1), ...
                                  timeSplit(1), timeSplit(2), timeSplit(3));
    ANSN_dat(i).lat = lat(i);
    ANSN_dat(i).lon = lon(i);
    ANSN_dat(i).dep = dep(i);
    ANSN_dat(i).prefmag = prefmag(i);
    ANSN_dat(i).evid = evid(i);
    ANSN_dat(i).oid = oid(i);
    ANSN_dat(i).mag = prefmag(i);
    ANSN_dat(i).magtype = magtype{i};
    ANSN_dat(i).mb = mb(i);
    ANSN_dat(i).ml = ml(i);
    ANSN_dat(i).ms = ms(i);
    ANSN_dat(i).mwp = mwp(i);
    ANSN_dat(i).mw = mw(i);

end

%% load mat file
if exist('mdat','var') ~= 1
    disp('Loading mdat');
    load mdat.mat;
end

%% now merge with GG Cat
disp('Merging ANSN Catalogue...');

t20 = 1/(60*24); % 1 minute threshold
j = 0;
mmin = 3.25;
for i = 1:length(mdat)
    % only consider larger events
    if mdat(i).GG_Mval >= mmin
         ind = find([ANSN_dat.dateNum] > mdat(i).MDAT_dateNum - t20 ...
                   & [ANSN_dat.dateNum] < mdat(i).MDAT_dateNum + t20 ...
                   & [ANSN_dat.lat] > mdat(i).MDAT_lat - 1 ...
                   & [ANSN_dat.lat] < mdat(i).MDAT_lat + 1 ...
                   & [ANSN_dat.lon] > mdat(i).MDAT_lon - 1 ...
                   & [ANSN_dat.lon] < mdat(i).MDAT_lon + 1 ...
                   & [ANSN_dat.prefmag] > mdat(i).GG_Mval - 0.5 ...
                   & [ANSN_dat.prefmag] < mdat(i).GG_Mval + 0.5);
               
        if length(ind) == 1
            disp(['Merging event ',datestr(ANSN_dat(ind).dateNum, 31)]);
            mdat(i).ANSN_evid = ANSN_dat(ind).evid;
            mdat(i).ANSN_oid = ANSN_dat(ind).oid;
            mdat(i).ANSN_datenum = ANSN_dat(ind).dateNum;
            mdat(i).ANSN_lat = ANSN_dat(ind).lat;
            mdat(i).ANSN_lon = ANSN_dat(ind).lon;
            mdat(i).ANSN_dep = ANSN_dat(ind).dep;
            mdat(i).ANSN_prefmag = ANSN_dat(ind).prefmag;
            mdat(i).ANSN_magType = ANSN_dat(ind).magtype;
            mdat(i).ANSN_ml = ANSN_dat(ind).ml;
            mdat(i).ANSN_mw = ANSN_dat(ind).mw;
            mdat(i).ANSN_mb = ANSN_dat(ind).mb;
            mdat(i).ANSN_ms = ANSN_dat(ind).ms;
            mdat(i).ANSN_mwp = ANSN_dat(ind).mwp;
            
        % manually select
        elseif length(ind) > 1
            disp(ind)
            txt = ['Multiple events found',char(10), ...
                   'Prev Event: ',datestr(mdat(i-1).MDAT_dateNum), ...
                   ' M ',num2str(mdat(i-1).GG_Mval),char(10), ...
                   'This Event: ',datestr(mdat(i).MDAT_dateNum), ...
                   ' M ',num2str(mdat(i).GG_Mval),char(10), ...
                   'Next Event: ',datestr(mdat(i+1).MDAT_dateNum), ...
                   ' M ',num2str(mdat(i+1).GG_Mval),char(10),char(10)];
            
            % cycle through found events
            for j = 1:length(ind)
                txt = [txt [num2str(j),' ',datestr(ANSN_dat(ind(j)).dateNum, 31), ...
                       ' M ',num2str(ANSN_dat(ind(j)).mag(1)),char(10)]];
            end
            txt = [txt [num2str(j+1),' None',char(10)]];
            txt = [txt ['Select event:',char(10)]];
            k = input(txt);
            if k > 0 & k <= j
                mdat(i).ANSN_evid = ANSN_dat(ind(k)).evid;
                mdat(i).ANSN_oid = ANSN_dat(ind(k)).oid;
                mdat(i).ANSN_datenum = ANSN_dat(ind(k)).dateNum;
                mdat(i).ANSN_lat = ANSN_dat(ind(k)).lat;
                mdat(i).ANSN_lon = ANSN_dat(ind(k)).lon;
                mdat(i).ANSN_dep = ANSN_dat(ind(k)).dep;
                mdat(i).ANSN_prefmag = ANSN_dat(ind(k)).prefmag;
                mdat(i).ANSN_magType = ANSN_dat(ind(k)).magtype;
                mdat(i).ANSN_ml = ANSN_dat(ind(k)).ml;
                mdat(i).ANSN_mw = ANSN_dat(ind(k)).mw;
                mdat(i).ANSN_mb = ANSN_dat(ind(k)).mb;
                mdat(i).ANSN_ms = ANSN_dat(ind(k)).ms;
                mdat(i).ANSN_mwp = ANSN_dat(ind(k)).mwp;

            else
                mdat(i).ANSN_evid = NaN;
                mdat(i).ANSN_oid = NaN;
                mdat(i).ANSN_datenum = NaN;
                mdat(i).ANSN_lat = NaN;
                mdat(i).ANSN_lon = NaN;
                mdat(i).ANSN_dep = NaN;
                mdat(i).ANSN_prefmag = NaN;
                mdat(i).ANSN_magType = '';
                mdat(i).ANSN_ml = NaN;
                mdat(i).ANSN_mw = NaN;
                mdat(i).ANSN_mb = NaN;
                mdat(i).ANSN_ms = NaN;
                mdat(i).ANSN_mwp = NaN;
                
            end
        else % no events
            mdat(i).ANSN_evid = NaN;
            mdat(i).ANSN_oid = NaN;
            mdat(i).ANSN_datenum = NaN;
            mdat(i).ANSN_lat = NaN;
            mdat(i).ANSN_lon = NaN;
            mdat(i).ANSN_dep = NaN;
            mdat(i).ANSN_prefmag = NaN;
            mdat(i).ANSN_magType = '';
            mdat(i).ANSN_ml = NaN;
            mdat(i).ANSN_mw = NaN;
            mdat(i).ANSN_mb = NaN;
            mdat(i).ANSN_ms = NaN;
            mdat(i).ANSN_mwp = NaN;    
        end
    else % events < mmin
        mdat(i).ANSN_evid = NaN;
        mdat(i).ANSN_oid = NaN;
        mdat(i).ANSN_datenum = NaN;
        mdat(i).ANSN_lat = NaN;
        mdat(i).ANSN_lon = NaN;
        mdat(i).ANSN_dep = NaN;
        mdat(i).ANSN_prefmag = NaN;
        mdat(i).ANSN_magType = '';
        mdat(i).ANSN_ml = NaN;
        mdat(i).ANSN_mw = NaN;
        mdat(i).ANSN_mb = NaN;
        mdat(i).ANSN_ms = NaN;
        mdat(i).ANSN_mwp = NaN;
    end
end

save mdat mdat;
