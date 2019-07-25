function T = getRr(mua, musp, rho, n)

% function T = getRr_r(mua, musp, r, n)
%  EITHER mua and musp are vectors, and r is scalar
%  OR r is vector and mua,musp are scalars
%  n is always a scaler.
%  Uses the diffusion math from  Farrell et al. (eq 15)
%n=refraction index ratio of the interface
%r=source detector distance 
% calculates relative diffuse reflectance at distance rho of the source

%slightly modified from https://omlc.org/news/apr08/skinspectra/index.html
% Steven L. Jacques
% Biomedical Engineering and Dermatology
% Oregon Health and Science University (OHSU)
% Portland, Oregon, USA

ri = 0.0636*n + 0.668 +  0.710/n - 1.440/n^2;
A = (1 + ri)/(1 - ri);

zo = 1./(mua + musp);

r1 = sqrt(zo.^2 + rho.^2);
r2 = sqrt((zo*(1+4*A/3)).^2 + rho.^2);
mueff = sqrt(3*mua./zo);

alb=musp./(mua+musp); %albedo

c = zo.*(mueff + 1./r1).*exp(-r1.*mueff)./(r1.^2);
d = (zo*(1+4*A/3)).*(mueff + 1./r2).*exp(-r2.*mueff)./(r2.^2);
T = alb.*( c + d )/(4*pi);


