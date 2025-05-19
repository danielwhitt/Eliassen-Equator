function ox = ox_units_inv(o,salt,temp);
% OX_UNITS  converts oxygen from micromol/kg to ml/l
%           using the potential density of water sample
%
%  Usage: ox = ox_units_inv(o,salt,potentialtemp);
%
% Irene Garcia Berdeal 2001
rho = sw_dens(salt,temp,0);
ox = o*1e-3.*rho*.022403;
