function [MLM92,WGW94,WGW96,BJ84,HB87,GS86,GG91,R35] = getMLestimates(rhyp,dep,maxh,maxv)

% Calculates local (Richter) magnitude from hypocentral distance, earthquake 
% depth, and maximum horizontal and vertical Wood-Anderson amplitudes using 
% the following relations:
%     - Michael-Leiba & Malafant (1992)
%     - Wilkie et al. (1994)
%     - Wilkie et al. (1996)
%     - Hutton & Boore (1987)
%     - Greenhalgh and Singh (1986)
%     - Gaull & Gregson (1991)
%     - Richter (1935)

repi = sqrt(rhyp^2 - dep^2);
if repi <= 0
    repi = 1;
end

% White (1968)
% W68 = log10(maxv) + ;

% Michael-Leiba & Malafant (1992)
MLM92 = log10(maxv) + 1.34*log10(rhyp/100) + 0.00055*(rhyp-100)+ 3.13;

% Wilkie et al. (1994)
WGW94 = log10(maxv) + 0.9 + log10(rhyp) + 0.0056*rhyp*exp(-0.0013*rhyp);

% Wilkie et al. (1996)
WGW96 = log10(maxv) + 0.75 + log10(rhyp) + 0.0056*rhyp*exp(-0.0013*rhyp);

% Bakun & Joyner (1984)
BJ84 = log10(maxh) + (log10(rhyp) + 0.00301 * rhyp + 0.7);

% Hutton & Boore (1987)
logA0 = 1.110 * log10(rhyp/100) + 0.00189*(rhyp - 100) + 3.0;
HB87 = log10(maxh) + logA0;

% Greenhalgh and Singh (SA,1986)
logA0 = 1.10 * log10(repi/100) + 0.0013 * (repi - 100) + 3.03;
GS86 = log10(maxh) + logA0;

% Gaull & Gregson (1991)
GG91 = log10(maxh) + 1.137 * log10(rhyp) + 0.000657 * rhyp + 0.66;


% % Allen (2010)
% c = textread('A10.GR_coeff.dat','%f');
% c1 = c(1);
% c2 = c(2);
% c3 = c(3);
% r1 = c(4);
% r2 = c(5);
% if rhyp <= r1
%     logASEA = c1*log10(rhyp/30);
% elseif rhyp > r1 & rhyp <= r2
%     logASEA = c1*log10(r1/30) + c2*log10(rhyp/r1);
% elseif rhyp > r2
%     logASEA = c1*log10(r1/30) + c2*log10(r2/r1) + c3*log10(rhyp/r2);
% end
% 
% A10 = log10(maxv) - logASEA + 2.29 + 0.18;

% Richter (1935)
logA0 = NaN;
if (repi <= 7.5)
  logA0 = 1.4;
elseif (repi > 7.5 && repi <= 12.5)
  logA0 = 1.5;
elseif (repi > 12.5 && repi <= 17.5)
  logA0 = 1.6;
elseif (repi > 17.5 && repi <= 22.5)
  logA0 = 1.7;
elseif (repi > 22.5 && repi <= 27.5)
  logA0 = 1.9;
elseif (repi > 27.5 && repi <= 32.5)
  logA0 = 2.1;
elseif (repi > 32.5 && repi <= 37.5)
  logA0 = 2.3;
elseif (repi > 37.5 && repi <= 42.5)
  logA0 = 2.4;
elseif (repi > 42.5 && repi <= 47.5)
  logA0 = 2.5;
elseif (repi > 47.5 && repi <= 52.5)
  logA0 = 2.6;
elseif (repi > 52.5 && repi <= 57.5)
  logA0 = 2.7;
elseif (repi > 57.5 && repi <= 72.5)
  logA0 = 2.8;
elseif (repi > 72.5 && repi <= 87.5)
  logA0 = 2.9;
elseif (repi > 87.5 && repi <= 105)
  logA0 = 3.0;
elseif (repi > 105 && repi <= 125)
  logA0 = 3.1;
elseif (repi > 125 && repi <= 145)
  logA0 = 3.2;
elseif (repi > 145 && repi <= 165)
  logA0 = 3.3;
elseif (repi > 165 && repi <= 185)
  logA0 = 3.4;
elseif (repi > 185 && repi <= 205)
  logA0 = 3.5;
elseif (repi > 205 && repi <= 215)
  logA0 = 3.6;
elseif (repi > 215 && repi <= 225)
  logA0 = 3.65;
elseif (repi > 225 && repi <= 245)
  logA0 = 3.7;
elseif (repi > 245 && repi <= 265)
  logA0 = 3.8;
elseif (repi > 265 && repi <= 285)
  logA0 = 3.9;
elseif (repi > 285 && repi <= 305)
  logA0 = 4.0;
elseif (repi > 305 && repi <= 325)
  logA0 = 4.1;
elseif (repi > 325 && repi <= 345)
  logA0 = 4.2;
elseif (repi > 345 && repi <= 375)
  logA0 = 4.3;
elseif (repi > 375 && repi <= 395)
  logA0 = 4.4;
elseif (repi > 395 && repi <= 425)
  logA0 = 4.5;
elseif (repi > 425 && repi <= 465)
  logA0 = 4.6;
elseif (repi > 465 && repi <= 505)
  logA0 = 4.7;
elseif (repi > 505 && repi <= 555)
  logA0 = 4.8;
elseif (repi > 555 && repi <= 605)
  logA0 = 4.9;
% Extend beyond intended distance (Eiby & Muir, 1968)
elseif (repi > 605 && repi <= 625)
  logA0 = 5.0;
elseif (repi > 625 && repi <= 675)
  logA0 = 5.1;  
elseif (repi > 675 && repi <= 725)
  logA0 = 5.2;
elseif (repi > 725 && repi <= 775)
  logA0 = 5.3;
elseif (repi > 775 && repi <= 825)
  logA0 = 5.4;
elseif (repi > 825 && repi <= 875)
  logA0 = 5.5;
elseif (repi > 875 && repi <= 925)
  logA0 = 5.55;
elseif (repi > 925 && repi <= 975)
  logA0 = 5.6;
elseif (repi > 975 && repi <= 1050)
  logA0 = 5.7;
elseif (repi > 1050 && repi <= 1150)
  logA0 = 5.8;
elseif (repi > 1150 && repi <= 1250)
  logA0 = 5.9;
elseif (repi > 1250 && repi <= 1350)
  logA0 = 6.0;
elseif (repi > 1350 && repi <= 1450)
  logA0 = 6.05;
elseif (repi > 1450 && repi <= 1550)
  logA0 = 6.15;
elseif (repi > 1550 && repi <= 1650)
  logA0 = 6.2;  
end
R35 = log10(maxh) + logA0;
