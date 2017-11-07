import matplotlib as mpl
mpl.use('Agg')
from matplotlib import pyplot
import numpy as np

bestb_r = [4.92642681305e-06,3.73879828925e-06,2.83747494445e-06,2.15343632834e-06,1.63430095806e-06,1.24031511235e-06,9.41308619039e-07,7.14384520074e-07,5.42165695925e-07,4.11464181513e-07,3.12271274152e-07,2.36991099206e-07,1.79858942374e-07,1.36499806365e-07,1.03593387639e-07,7.86198182134e-08,5.96667022554e-08,4.52826709466e-08,3.43662413129e-08,2.60814681928e-08,1.9793930238e-08,1.50221479624e-08,1.14007135869e-08,8.65230928473e-09,6.56647107111e-09,4.983472147e-09,3.78209153303e-09,2.87033135578e-09,2.1783719458e-09,1.65322527125e-09,1.08027714192e-09,5.92273661584e-10,1.85000115955e-10,4.82630437346e-11]

upperb_r=[3.34803245492e-06,2.63323524455e-06,2.07104559065e-06,1.62888213176e-06,1.28111955195e-06,1.00760348118e-06,7.92482460933e-07,6.23289282559e-07,4.902184577e-07,3.8555794716e-07,3.03242214329e-07,2.38500700682e-07,1.8758135094e-07,1.47533164975e-07,1.16035174385e-07,9.12619321679e-08,7.17777200507e-08,5.64533423004e-08,4.44006838701e-08,3.49212402278e-08,2.74656359486e-08,2.16017859945e-08,1.69898544867e-08,1.33625597232e-08,1.05096840291e-08,8.26589071854e-09,6.5011421068e-09,5.11316325512e-09,4.02151468833e-09,3.16293057381e-09,2.14186853413e-09,1.21697048076e-09,3.93939586384e-10,1.0650557256e-10]

lowerb_r=[1.74154752419e-06,1.27566115239e-06,9.34405379768e-07,6.84439917373e-07,5.01343432558e-07,3.67227613395e-07,2.68989501571e-07,1.97031348723e-07,1.44322927671e-07,1.05714687468e-07,7.74346483037e-08,5.67198835049e-08,4.15465848336e-08,3.04323387968e-08,2.22912965856e-08,1.63280879193e-08,1.19601143019e-08,8.76062983138e-09,6.41704862553e-09,4.70040554788e-09,3.44298658213e-09,2.52194336935e-09,1.84729106736e-09,1.35311693713e-09,9.91140745437e-10,7.25997842693e-10,5.31784078115e-10,3.89524994575e-10,2.85322046378e-10,2.08994727637e-10,1.31807001908e-10,6.97470293369e-11,2.10269022378e-11,5.29441526934e-12]

bad_mergedb_r=[2.84931951177e-05,2.21912494373e-05,1.72874807237e-05,1.34706705801e-05,1.04991092599e-05,8.18501702445e-06,6.38246958726e-06,4.97803659896e-06,3.88352136425e-06,3.03032845556e-06,2.36509290537e-06,1.84628739122e-06,1.44158780558e-06,1.12582709427e-06,8.79405544059e-07,6.87055343858e-07,5.36880197412e-07,4.1960854323e-07,3.28012711962e-07,2.56457026567e-07,2.00546126778e-07,1.56851217454e-07,1.22696919129e-07,9.59952623733e-08,7.51163662494e-08,5.87876694263e-08,4.60153907325e-08,3.60233007262e-08,2.82049775096e-08,2.20865685848e-08,1.48933544258e-08,8.42748412854e-09,2.71719760116e-09,7.3179868136e-10]

mags=np.arange(4.5,7.9,0.1)
print len(mags)
print len(bestb_r)
sum_rates=np.add(bestb_r,upperb_r)
sum_rates=np.add(sum_rates,lowerb_r)
pyplot.plot(mags,bestb_r,label='Best b w0.5')
pyplot.plot(mags,upperb_r,label='Upper b w0.3')
pyplot.plot(mags,lowerb_r,label='Lower b w0.2')
pyplot.plot(mags,sum_rates,label='Weighted sum')
pyplot.plot(mags,bad_mergedb_r,label='Old buggy weighted sum')
pyplot.xlabel('Magnitude')
pyplot.ylabel('Incremental rate')
pyplot.legend()
pyplot.savefig('test_rates.png')
