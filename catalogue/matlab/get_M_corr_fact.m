% get_M_corr_fact.m
% 
% Develop generic relations to correct magnitudes for unknown types that 
% do not have a RevML field prior to 1986
% *************************************************************************

if exist('mdat','var') ~= 1
    disp('Loading mdat');
    load mdat_ml_rev.mat;
end

ind = [];
for i = 1:length(mdat)
    if ~isnan(mdat(i).MDAT_MLrevdist(1)) ...
       & mdat(i).MDAT_dateNum > datenum(1986,1,1) ...
       & mdat(i).MDAT_lat < -10
        ind = [ind i];
    end
end
%   
plot([1, 7], [1, 7], 'k--');
hold on;
plot(([mdat(ind).MDAT_prefML]),([mdat(ind).MDAT_MLrev]),'b+');
xlabel('prefML');
ylabel('MLrev');

mrange = (1:0.2:6);
for i = 1:length(mrange)
    mind = find(([mdat(ind).MDAT_prefML]) >= mrange(i)-0.1 ...
           & ([mdat(ind).MDAT_prefML]) < mrange(i)+0.1);
    medrev(i) = median(([mdat(ind(mind)).MDAT_MLrev]));
    stdrev(i) = std(([mdat(ind(mind)).MDAT_MLrev]));
end
delnan = find(isnan(medrev));
mrange(delnan) = [];
medrev(delnan) = [];
stdrev(delnan) = [];

hold on;
errorbar((mrange),(medrev),(stdrev),'rs');
mcorr = polyfit(mrange,medrev,1);
dlmwrite('mcorr.dat',mcorr);
