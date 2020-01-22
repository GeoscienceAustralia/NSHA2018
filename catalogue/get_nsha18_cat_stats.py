from catalogue.parsers import parse_NSHA2018_catalogue
from misc_tools import dictlist2array
from numpy import where, unique

catfile = 'data//NSHA18CAT.MW.V0.1.csv'

# parse catalogue
cat = parse_NSHA2018_catalogue(catfile)

# get data arrays
mw_pref = dictlist2array(cat, 'prefmag')
mw_src = dictlist2array(cat, 'mw_src')
ml_region = dictlist2array(cat, 'ml_region')
auth = dictlist2array(cat, 'auth')

# oge3 
oge3 = len(where((mw_pref >= 3.0) & (ml_region!='Other') & (mw_src!='nan'))[0])
print 'tot', oge3

# get ML2MW
mlge3 = len(where((mw_pref >= 3.0) & (mw_src=='ML2MW') & (ml_region!='Other') & (mw_src!='nan'))[0])
print 'mlge3', mlge3, 100*mlge3 / float(oge3)

# get MS2MW
msge3 = len(where((mw_pref >= 3.0) & (mw_src=='MS2MW') & (ml_region!='Other') & (mw_src!='nan'))[0])
print 'msge3', msge3, 100*msge3 / float(oge3)

# get ML2MW
mbge3 = len(where((mw_pref >= 3.0) & (mw_src=='mb2MW') & (ml_region!='Other') & (mw_src!='nan'))[0])
print 'mbge3', mbge3, 100*mbge3 / float(oge3)

# get MW native
mwge3 = len(where((mw_pref >= 3.0) & (mw_src!='ML2MW') & (mw_src!='MS2MW') \
                & (mw_src!='mb2MW') & (ml_region!='Other') & (mw_src!='nan'))[0])
print 'mwge3', mwge3, 100*mwge3 / float(oge3), '\n'

##################

# oge4 
oge4 = len(where((mw_pref >= 4.0) & (ml_region!='Other') & (mw_src!='nan'))[0])
print 'tot', oge4

# get ML2MW
mlge4 = len(where((mw_pref >= 4.0) & (mw_src=='ML2MW') & (ml_region!='Other') & (mw_src!='nan'))[0])
print 'mlge4', mlge4, 100*mlge4 / float(oge4)

# get MS2MW
msge4 = len(where((mw_pref >= 4.0) & (mw_src=='MS2MW') & (ml_region!='Other') & (mw_src!='nan'))[0])
print 'msge4', msge4, 100*msge4 / float(oge4)

# get ML2MW
mbge4 = len(where((mw_pref >= 4.0) & (mw_src=='mb2MW') & (ml_region!='Other') & (mw_src!='nan'))[0])
print 'mbge4', mbge4, 100*mbge4 / float(oge4)

# get MW native
mwge4 = len(where((mw_pref >= 4.0) & (mw_src!='ML2MW') & (mw_src!='MS2MW') \
                & (mw_src!='mb2MW') & (ml_region!='Other') & (mw_src!='nan'))[0])
print 'mwge4', mwge4, 100*mwge4 / float(oge4), '\n'

##################

# get unique auth

uauth = unique(auth)
authtxt = ''
for a in uauth:
    authtxt += a+ '\n'

csvfile = 'authlist.csv'
f = open(csvfile, 'wb')
f.write(authtxt)
f.close()
