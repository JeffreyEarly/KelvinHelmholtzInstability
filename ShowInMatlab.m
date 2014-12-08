file = '/Users/jearly/Desktop/KelvinHelmholtzInstability.nc';
x = ncread(file, 'x');
y = ncread(file, 'y');
t = ncread(file, 't');

iTime = length(t)-1;

zeta = double(ncread(file, 'zeta'));
psi = double(ncread(file, 'psi'));
u = double(ncread(file, 'u'));
v = double(ncread(file, 'v'));

xunits = ncreadatt(file, 'x', 'units');
yunits = ncreadatt(file, 'y', 'units');
tunits = ncreadatt(file, 't', 'units');

U0 = max(max(u(:,:,1)));

figure

subplot(2,1,1)
pcolor(x, y, zeta(:,:,1)), axis equal tight, shading interp
caxis([min(min(zeta(:,:,1))) max(max(zeta(:,:,1)))])
title(sprintf('vorticity, t=0 %s', tunits))
xlabel(xunits), ylabel(xunits)

subplot(2,1,2)
pcolor(x, y, zeta(:,:,iTime)), axis equal tight, shading interp
caxis([min(min(zeta(:,:,1))) max(max(zeta(:,:,1)))])
title(sprintf('vorticity, t=%d %s', round(t(iTime)), tunits ))
xlabel(xunits), ylabel(xunits)


figure

subplot(2,1,1)
pcolor(x, y, psi(:,:,1)), axis equal tight, shading interp
caxis([min(min(psi(:,:,1))) max(max(psi(:,:,1)))])
title(sprintf('stream function, t=0 %s', tunits))
xlabel(xunits), ylabel(xunits)

subplot(2,1,2)
pcolor(x, y, psi(:,:,iTime)), axis equal tight, shading interp
caxis([min(min(psi(:,:,1))) max(max(psi(:,:,1)))])
title(sprintf('stream function, t=%d %s', round(t(iTime)), tunits ))
xlabel(xunits), ylabel(xunits)

return;

figure
subplot(2,1,1)
pcolor(x, y, u(:,:,1)), axis equal tight, shading interp
caxis([-U0 U0])
title(sprintf('u, t=0 %s', tunits))
xlabel(xunits), ylabel(xunits)

subplot(2,1,2)
pcolor(x, y, u(:,:,iTime)), axis equal tight, shading interp
caxis([-U0 U0])
title(sprintf('u, t=%d %s', round(t(iTime)), tunits ))
xlabel(xunits), ylabel(xunits)

figure
subplot(2,1,1)
pcolor(x, y, v(:,:,1)), axis equal tight, shading interp
caxis([-U0 U0])
title(sprintf('v, t=0 %s', tunits))
xlabel(xunits), ylabel(xunits)

subplot(2,1,2)
pcolor(x, y, v(:,:,iTime)), axis equal tight, shading interp
caxis([-U0 U0])
title(sprintf('v, t=%d %s', round(t(iTime)), tunits ))
xlabel(xunits), ylabel(xunits)