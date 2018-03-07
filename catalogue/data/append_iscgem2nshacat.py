from hmtk.parsers.catalogue.csv_catalogue_parser import CsvCatalogueParser
from tools.nsha_tools import toYearFraction, get_shapely_centroid
from tools.mfd_tools import parse_hmtk_cat
from catalogue.writers import ggcat2hmtk_csv
from os import path
from datetime import timedelta

###############################################################################
# parse NSHA-Cat and ISC-GEM catalogues
###############################################################################

# parse NSHA-Cat catalogue
hmtk_csv = path.join('NSHA18CAT_V0.1_hmtk_mx_orig.csv')


nshaCat, full_neq = parse_hmtk_cat(hmtk_csv)
nshaMaxYear = toYearFraction(nshaCat[-1]['datetime'])

# parse ISC-GEM catalogue
hmtk_csv = path.join('ISC-GEM_V4_hmtk_full.csv')
iscCat, crust_neq = parse_hmtk_cat(hmtk_csv)
iscMaxYear = toYearFraction(iscCat[-1]['datetime'])

allCat = nshaCat
for ig in iscCat:
    if ig['lon'] > 100. and ig['lon'] < 160. \
       and ig['lat'] > -50. and ig['lat'] < 0.:
       
       addEvent = True
       
       for nc in nshaCat:
           if ig['datetime'] >= nc['datetime'] - timedelta(seconds=10) \
              and ig['datetime'] <= nc['datetime'] + timedelta(seconds=10) \
              and ig['prefmag'] >= nc['prefmag']-0.75 \
              and ig['prefmag'] <= nc['prefmag']+0.75:
              
              addEvent = False
              
       if addEvent == True:
           allCat.append(ig)
           print 'Adding:', ig['datetime']

# write to HMTK
mergedCatCSV = 'merged_NSHA18-ISCGEM_hmtk.csv'
ggcat2hmtk_csv(allCat, mergedCatCSV, 'mw')