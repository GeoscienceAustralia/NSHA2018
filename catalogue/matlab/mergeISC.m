% script to merge ISC's earthquake catalogue to mdat

%% parse ISC catalogue
disp('Parsing ISC Catalogue...');
ISC_catfile = fullfile('..','data','isc_events_GT_32.txt');

fid = fopen(ISC_catfile);

% skip headers
skipline = 1;
while skipline == 1
    line = fgetl(fid);
    if strfind(line, 'ISC Bulletin') == 1
        skipline = 0;
    end
end

% Now get data
skipline = 1;
i = 0;
line = fgetl(fid);
while skipline == 1
   i = i + 1; % inc event
   nextEvent = 1;
   
   % get datenum
   ISC_dat(i).dateNum = datenum(line(1:19), 'yyyy/mm/dd HH:MM:SS');
   
   % get event loc
   ISC_dat(i).lat = str2double(line(37:44));
   ISC_dat(i).lon = str2double(line(47:54));
   
   % now loop through mag lines
   ISC_dat(i).mb = NaN;
   ISC_dat(i).ms = NaN;
   ISC_dat(i).ml = NaN;
   ISC_dat(i).mw = NaN;
   ISC_dat(i).mbSRC = '';
   ISC_dat(i).msSRC = '';
   ISC_dat(i).mlSRC = '';
   ISC_dat(i).mwSRC = '';
   while nextEvent == 1
   line = fgetl(fid);
       if strfind(line, 'mb') == 1
           ISC_dat(i).mb = str2double(line(8:10));
           ISC_dat(i).mbSRC = deblank(line(21:26));
       elseif strfind(line, 'mL') == 1
           ISC_dat(i).ml = str2double(line(8:10));
           ISC_dat(i).mlSRC = deblank(line(21:26));
       elseif strfind(line, 'ML') == 1
           ISC_dat(i).ml = str2double(line(8:10));
           ISC_dat(i).mlSRC = deblank(line(21:26));
       elseif strfind(line, 'MS') == 1
           ISC_dat(i).ms = str2double(line(8:10));
           ISC_dat(i).msSRC = deblank(line(21:26));
       elseif strfind(line, 'MW') == 1
           ISC_dat(i).mw = str2double(line(8:10));
           ISC_dat(i).mwSRC = deblank(line(21:26));
       end
   
       if strcmp(line,'')
           nextEvent = 0;
       end
   end
   
   line = fgetl(fid);
   
   if strfind(line, 'STOP') == 1
        skipline = 0;
   end
end
%% load mat file
if exist('mdat','var') ~= 1
    disp('Loading mdat');
    load mdat.mat;
end

%% now merge with ISC Cat
disp('Merging ISC Catalogue...');

t20 = 10/(60*60*24); % 20 sec threshold
j = 0;
mmin = 3.25;
for i = 1:length(mdat)
    % only consider larger events
    if mdat(i).GG_Mval >= mmin
         ind = find([ISC_dat.dateNum] > mdat(i).MDAT_dateNum - t20 ...
                   & [ISC_dat.dateNum] < mdat(i).MDAT_dateNum + t20 ...
                   & [ISC_dat.lat] > mdat(i).MDAT_lat - 1 ...
                   & [ISC_dat.lat] < mdat(i).MDAT_lat + 1 ...
                   & [ISC_dat.lon] > mdat(i).MDAT_lon - 1 ...
                   & [ISC_dat.lon] < mdat(i).MDAT_lon + 1);
               
        if length(ind) == 1
            disp(['Merging event ',datestr(ISC_dat(ind).dateNum, 31)]);
            mdat(i).ISC_lat = ISC_dat(ind).lat;
            mdat(i).ISC_lon = ISC_dat(ind).lon;
            mdat(i).ISC_ml = ISC_dat(ind).ml;
            mdat(i).ISC_mw = ISC_dat(ind).mw;
            mdat(i).ISC_mb = ISC_dat(ind).mb;
            mdat(i).ISC_ms = ISC_dat(ind).ms;
            mdat(i).ISC_mlSRC = ISC_dat(ind).mlSRC;
            mdat(i).ISC_mwSRC = ISC_dat(ind).mwSRC;
            mdat(i).ISC_mbSRC = ISC_dat(ind).mbSRC;
            mdat(i).ISC_msSRC = ISC_dat(ind).msSRC;
            
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
                txt = [txt [num2str(j),' ',datestr(ISC_dat(ind(j)).dateNum, 31), ...
                       ' M ',num2str(ISC_dat(ind(j)).ml(1)),char(10)]];
            end
            txt = [txt [num2str(j+1),' None',char(10)]];
            txt = [txt ['Select event:',char(10)]];
            k = input(txt);
            if k > 0 & k <= j
                mdat(i).ISC_lat = ISC_dat(ind(k)).lat;
                mdat(i).ISC_lon = ISC_dat(ind(k)).lon;
                mdat(i).ISC_ml = ISC_dat(ind(k)).ml;
                mdat(i).ISC_mw = ISC_dat(ind(k)).mw;
                mdat(i).ISC_mb = ISC_dat(ind(k)).mb;
                mdat(i).ISC_ms = ISC_dat(ind(k)).ms;
                mdat(i).ISC_mlSRC = ISC_dat(ind).mlSRC;
                mdat(i).ISC_mwSRC = ISC_dat(ind).mwSRC;
                mdat(i).ISC_mbSRC = ISC_dat(ind).mbSRC;
                mdat(i).ISC_msSRC = ISC_dat(ind).msSRC;

            else
                mdat(i).ISC_lat = NaN;
                mdat(i).ISC_lon = NaN;
                mdat(i).ISC_ml = NaN;
                mdat(i).ISC_mw = NaN;
                mdat(i).ISC_mb = NaN;
                mdat(i).ISC_ms = NaN;
                mdat(i).ISC_mlSRC = '';
                mdat(i).ISC_mwSRC = '';
                mdat(i).ISC_mbSRC = '';
                mdat(i).ISC_msSRC = '';
                
            end
        else % no events
            mdat(i).ISC_lat = NaN;
            mdat(i).ISC_lon = NaN;
            mdat(i).ISC_ml = NaN;
            mdat(i).ISC_mw = NaN;
            mdat(i).ISC_mb = NaN;
            mdat(i).ISC_ms = NaN; 
            mdat(i).ISC_mlSRC = '';
            mdat(i).ISC_mwSRC = '';
            mdat(i).ISC_mbSRC = '';
            mdat(i).ISC_msSRC = '';
        end
    else % events < mmin
        mdat(i).ISC_lat = NaN;
        mdat(i).ISC_lon = NaN;
        mdat(i).ISC_ml = NaN;
        mdat(i).ISC_mw = NaN;
        mdat(i).ISC_mb = NaN;
        mdat(i).ISC_ms = NaN;
        mdat(i).ISC_mlSRC = '';
        mdat(i).ISC_mwSRC = '';
        mdat(i).ISC_mbSRC = '';
        mdat(i).ISC_msSRC = '';
    end
end

%% Append gaps in ISC with PDE

disp('Appending PDE...');

pdefile = fullfile('..','data','usgs_20150201_20180122.mb.csv');
[time, latitude, longitude, depth, mag] = ...
    textread(pdefile, '%s%f%f%f%f%*[^\n]', 'headerlines',1,'delimiter',',','emptyvalue',NaN);

usgsDateNum = datenum(time, 'yyyy-mm-dd HH:MM:SS');

for i = 1:length(mdat)
    % only consider larger events
    if mdat(i).GG_Mval >= mmin && mdat(i).MDAT_dateNum >= datenum(2015,2,1)
         ind = find(usgsDateNum > mdat(i).MDAT_dateNum - t20 ...
                   & usgsDateNum < mdat(i).MDAT_dateNum + t20 ...
                   & latitude > mdat(i).MDAT_lat - 1 ...
                   & latitude < mdat(i).MDAT_lat + 1 ...
                   & longitude > mdat(i).MDAT_lon - 1 ...
                   & longitude < mdat(i).MDAT_lon + 1);
               
        if length(ind) == 1
            disp(['Merging PDE event ',datestr(usgsDateNum(ind), 31)]);
            mdat(i).ISC_lat = latitude(ind);
            mdat(i).ISC_lon = longitude(ind);
            mdat(i).ISC_mb = mag(ind);
            mdat(i).ISC_mbSRC = 'PDE';
            
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
                txt = [txt [num2str(j),' ',datestr(usgsDateNum(ind(j)), 31), ...
                       ' M ',num2str(mag(ind(j))),char(10)]];
            end
            txt = [txt [num2str(j+1),' None',char(10)]];
            txt = [txt ['Select event:',char(10)]];
            k = input(txt);
            if k > 0 & k <= j
                mdat(i).ISC_lat = latitude(ind(k));
                mdat(i).ISC_lon = longitude(ind(k));
                mdat(i).ISC_mb = mag(ind(k));
                mdat(i).ISC_mbSRC = 'PDE';
            end
        
        end
    
    end
end

%% save mat file

disp('Saving mat file...');
save mdat mdat;
