% script to merge magnitudes from alternative studies to main catalogue

%% parse Allen SEA catalogue
disp('Parsing Allen ML Catalogue...');
alt_catfile = fullfile('..','data','2017_sea_updated_ev_ml.csv');
      
[DATETIME, R35, R35STD, BJ84, BJ84STD, MLM92, MLM92STD, WGW96, WGW96STD, A16, A16STD] = ...
 textread(alt_catfile, '%s%f%f%f%f%f%f%f%f%f%f', ...
          'headerlines',1,'delimiter',',','emptyvalue',NaN);

for i = 1:length(MLM92)
    % get datetime
    alt_dat(i).dateNum = datenum(DATETIME(i), 'yyyymmddHHMM');
    alt_dat(i).ml = MLM92(i);

end

%% load mat file
if exist('mdat','var') ~= 1
    disp('Loading mdat');
    load mdat.mat;
end

%% now merge with alt mags
disp('Merging alternate ML...');

t20 = 1/(60*24); % 1 minute threshold
j = 0;
mmin = 0;
for i = 1:length(mdat)
    % only consider larger events
    if mdat(i).GG_Mval > mmin
        ind = find([alt_dat.dateNum] > mdat(i).MDAT_dateNum - t20 ...
                   & [alt_dat.dateNum] < mdat(i).MDAT_dateNum + t20);
               
        if length(ind) == 1
            disp(['Merging event ',datestr(alt_dat(ind).dateNum, 31)]);
            mdat(i).Allen_ML = alt_dat(ind).ml;
            
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
                txt = [txt [num2str(j),' ',datestr(alt_dat(ind(j)).dateNum, 31), ...
                       ' M ',num2str(alt_dat(ind(j)).ml(1)),char(10)]];
            end
            txt = [txt [num2str(j+1),' None',char(10)]];
            txt = [txt ['Select event:',char(10)]];
            k = input(txt);
            if k > 0 & k <= j
                mdat(i).Allen_ML = alt_dat(ind(k)).ml;

            else
                mdat(i).Allen_ML = NaN;
                
            end
        else % no events
            mdat(i).Allen_ML = NaN;
        end
    else % events < mmin
        mdat(i).Allen_ML = NaN;
    end
end

save mdat mdat;
