figure(102)
subplot(2,4,1); 
imagesc(imread('roi.tif'));colormap(flipud(bone)); title('moire image')%freezeColors;
colormap default
ax_ft=subplot(2,4,2); plot(log(fftshift(imread('fft.tif'))));set(ax_ft, 'YTick','');title('1D FFT')
%subplot_tight(2,4,3,.06); imagesc(amp);
cifimg_re = imread('ifftre.tif');
cifimg_im = imread('ifftim.tif');
subplot(2,4,3); imagesc(abs(cifimg_re + 1j*cifimg_im));title('1st order FFT filter')
cwrap = imread('dphi_w.tif');
subplot(2,4,4); imagesc(cwrap);title('wrapped fringe phase');
%subplot(2,4,5); imagesc(unwrap);
cwrap = imread('alpha2D.tif');
subplot(2,4,5); imagesc(cwrap);title('fringe deflection \alpha');
cpolys = imread('polys.tif').';
%ax1= subplot(2,4,6); plot(cpolys(1,:));title('\delta\alpha/\delta x')%xlim(ax1,[0 172])
cwpa = imread('wpa.tif');
ax2= subplot(2,4,6); plot(cwpa);title('1D \alpha'), ylabel('radians')%xlim(ax2,[0 0.0005])
%subplot(2,4,7); imagesc(out.alpha2D);
%ax2= subplot(2,4,8); plot(out.hp*1E9);title('height profile');ylabel('nanometers')%xlim(ax2,[0 0.0005])
cdphi = imread('hp.tif');
subplot(2,4,7); plot(cdphi*1E9);title('hp');
cphiy = imread('phi_y.tif');
ax2 = subplot(2,4,8); plot(cphiy*1E9);title('spherical wavefront');ylabel('nanometers')
