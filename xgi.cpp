#include "xgi.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#ifdef HAVE_TIFF
#include <tiffio.h>
#endif
#ifdef HAVE_PNG
#include <png.h>
#endif
#ifdef HAVE_HDF5
#include <hdf5.h>
#endif
#include <fftw3.h>
#include <string.h>

using namespace std;


Xgi::Xgi()
{
};

Xgi::~Xgi()
{
}

int Xgi::compute(Matrix *put_in, Xgi_ip *ip, Xgi_out *put_out)
{
    FLOATING wavelength = 12.3984191e-7 / (ip->energy);
    FLOATING magnification = 2* (ip->period_g2) / (ip->period_g1);
    FLOATING beta = - ((ip->beta_rot) - (ip->beta_fix))/2;
    Matrix hann,tmp_fit_x,tmp_fit_y;
    int rv,r,c,width,height,width_fft,height_fft;
    FLOATING val_M,val_m,val_fwhmd,val,sum;
    int pos_ind_M, pos_ind_fwhm, pos_ind_p,freq,filter_width,i_x;
    FLOATING G_width,visib;
    FLOATING Sx,Sxx,Sxy,Sy,Sxxx,Sxxxx,Sxxy,fit_a,fit_a_tmp,fit_b,fit_c,Axx,Ax,Ay;
    FLOATING alpha,alpha_mul,alpha_add;
    FLOATING polys_avg,polys_std;
    FFTC *fft_in,*fft_out,*fft_out_amp;
    FFTPLAN fft_plan_fw,fft_plan_bw;
    
    put_in->crop(&(put_out->roi),ip->roi,ip->cam_orientation);
    width = put_out->roi.width;
    height = put_out->roi.height;
    /*
    width_fft = width;
    height_fft = height * (ip->n_zeropadding);
    */
    width_fft = width * (ip->n_zeropadding);;
    height_fft = height;
    
    if (ip->verbose)
        printf("[I]Work %d %d| FFT %d %d\n",width,height,width_fft,height);
    /**HANN 
    nn=1:1:cols;
    HANN = .54 - 0.46*cos(2*pi.*nn/(cols+1));
    HANN=HANN';
    hanning_window = repmat(HANN,1,rows);
    */
    hann.resize(height,1);
    for (c=0;c<height;++c)
        hann.data[c] = 0.54 - 0.46*COS(2*C_PI*(c+1)/(height+1));
    if (ip->verbose > 9000)
        hann.save_tiff("hann.tif");    
    /**FFT 
    fimg = fft((img).*hanning_window,ip.n_zeropadding*size(img,2),dim_ana);
    */
    fft_in = (FFTC*) FFTMALLOC(sizeof(FFTC) * width_fft * height_fft);
    fft_out = (FFTC*) FFTMALLOC(sizeof(FFTC) * width_fft * height_fft);
    fft_out_amp = (FFTC*) FFTMALLOC(sizeof(FFTC) * width_fft * height_fft);
    
    if (fft_in == NULL || fft_out == NULL || fft_out_amp == NULL)
    {
        printf("[E]FFT buffer creation failed\n");
        return -1;
    }


    int rank = 1;
    int n[] = {width_fft}; 
    int howmany = (height);
    int idist=width_fft,odist=width_fft;
    int istride=1, ostride = 1;
    int *inembed = n, *onembed = n;
    /*
    int rank = 1;
    int n[] = {height_fft}; 
    int howmany = (width_fft);
    int idist=1,odist=1;
    int istride=width_fft, ostride = width_fft;
    int *inembed = n, *onembed = n;
    */
    fft_plan_fw = FFTPLANMANY(rank,n,howmany,fft_in, inembed,istride,idist,
                                          fft_out,onembed,ostride,odist,
                            FFTW_FORWARD, FFTW_MEASURE);
    fft_plan_bw = FFTPLANMANY(rank,n,howmany,fft_out, inembed,istride,idist,
                                          fft_in,onembed,ostride,odist,
                            FFTW_BACKWARD, FFTW_MEASURE);
    if (fft_plan_fw == NULL || fft_plan_bw == NULL)
    {
        printf("[E]Plan creation failed\n");
        return -1;
    }
    memset(fft_in,0,sizeof(FFTC) * width_fft * height_fft);
    for (r=0;r<height;++r)
        for (c=0;c<width;++c)
            //fft_in[r*width_fft+c][0] = (1-2*(c&1)) * hann.data[c] * work_in.data[r*width + c]; //makes fourier shifted
            fft_in[r*width_fft+c][0] = hann.data[r] * put_out->roi.data[r*width + c];


    
    FFTEXEC(fft_plan_fw);
    
    put_out->fft_re.zero(width_fft,height_fft);
    put_out->fft_im.zero(width_fft,height_fft);
    for (c=0;c<width_fft*height_fft;++c)
    {
        //put_out->fft_re.data[c] = SQRT(fft_out[c][0]*fft_out[c][0]+fft_out[c][1]*fft_out[c][1]);
        put_out->fft_re.data[c] = fft_out[c][0];
        put_out->fft_im.data[c] = fft_out[c][1];
    }
    
        
    /**AVG 
    mean_proj = mean(abs(fimg),dim_ana_t);
    out.fft = mean_proj;
    % get the position of the maximum value -> DC component position
    [val_M, pos_ind_M] = max(mean_proj);
    % get the position where the intensity is half the maximum value in order to have an estimate for the FWHM (Gaussian approximation of the DC component)
    [~,pos_ind_fwhm] = min(abs(mean_proj-val_M/2));
    % calculate from the FWHM the width of the Gaussian
    G_width = 2*abs(pos_ind_M - pos_ind_fwhm) / 2.35;
    %ensure max G_width
    fprintf(2, 'pos_ind_M=%f pos_ind_fwhm=%f  G_width=%f ', pos_ind_M, pos_ind_fwhm, G_width);
    if G_width>10
        G_width=10;    
    end
    */
    put_out->fft.resize(width_fft,1);
    val_M = 0;
    pos_ind_M = -1;
    for (c=0;c<width_fft;++c)
    {
        val = 0;
        for (r=0;r<height_fft;++r)
            val += SQRT(fft_out[r*width_fft+c][0]*fft_out[r*width_fft+c][0]+fft_out[r*width_fft+c][1]*fft_out[r*width_fft+c][1]);
        val /= (C_ONE*height_fft);
        put_out->fft.data[c] = val;
        if(val_M < val)
        {
            val_M = val;
            pos_ind_M = c;
        }
    }
    if (pos_ind_M ==-1)
    {
        printf("[E]Maximum fft abs not found!!");
        return -1;
    }
    val_fwhmd = INFINITY;
    pos_ind_fwhm = -1;
    for (c=0;c<width_fft/2;++c)
    {
        val = ABS(put_out->fft.data[c] - val_M*0.5);
        if (val_fwhmd > val)
        {
            val_fwhmd = val;
            pos_ind_fwhm = c;
        }
    }
    G_width = 2*ABS(pos_ind_M - pos_ind_fwhm) / 2.35;
    if (ip->verbose>1337){
        printf("[I]Max: %f @ %d\n",val_M,pos_ind_M);
        printf("[I]FWHM @ %d\n",pos_ind_fwhm);
        printf("[I]G width: %f\n",G_width);
    }
    G_width = FMIN(10,G_width);
    
    /** Visibility
     
    [val_m,pos_ind_p] = max(mean_proj(pos_ind_M+round(ip.excl_para*G_width):length(mean_proj)/2));
    pos_ind_p = pos_ind_p+pos_ind_M+round(ip.excl_para*G_width)-1;
    %shift
    freq = pos_ind_p-pos_ind_M;

    % finge visibility
    visib = val_m/val_M;
    fprintf(2,  'freq=%f  fringe visib=%.3f \n',freq, visib);
    */
    val_m = -1;
    pos_ind_p = -1;
    for (c=ROUND(pos_ind_M + (ip->excl_para) * G_width);c<width_fft/2;++c)
    {
        if (val_m < put_out->fft.data[c])
        {
            val_m = put_out->fft.data[c];
            pos_ind_p = c;
        }
    }
    //pos_ind_p already real index, not shifted by pos_ind_M+round()
    freq = pos_ind_p - pos_ind_M;
    visib = val_m / val_M;
    if (ip->verbose>1337){
        printf("[I]Freq: %d\n",freq);
        printf("[I]Fringe visibility: %f\n",visib);
    }


    /** line_mask = zeros(size(fimg,1), size(fimg,2));
    filter_width = round(ip.filter_para*G_width);
    % 0-order and first order component processing

    line_mask(:, 1:filter_width+1) = 1;
    line_mask(:, end-filter_width:end) = 1;
    fimg_amp = fimg;
    fimg_amp(~line_mask) = 0;

    % 0-rder -> intensity on camera
    amp = abs(ifft(fimg_amp,[],dim_ana));

    amp = amp(1:size(img,1), 1:size(img,2));
    amp = amp ./ (hanning_window);
    out.amp = amp;
    */
    filter_width = ROUND(ip->filter_para * G_width);
    if (ip->verbose>1337){
        printf("[I]Filter width: %d\n",filter_width);
    }
    for (r=0;r<height_fft;++r)
        for (c=0;c<width_fft;++c)
        {
            fft_out_amp[c+r*width_fft][0] = fft_out[c+r*width_fft][0] * ((c<filter_width+2) || (c>(width_fft - filter_width-2)));
            fft_out_amp[c+r*width_fft][1] = fft_out[c+r*width_fft][1] * ((c<filter_width+2) || (c>(width_fft - filter_width-2)));
        }
    FFTEXECA(fft_plan_bw,fft_out_amp,fft_in);
    put_out->amp.zero(width,height);
    //put_out->amp.zero(width_fft,height);
    for (r=0;r<height;++r)
        for (c=0;c<width;++c)
            put_out->amp.data[c+r*width] = SQRT(fft_in[r*width_fft+c][0]*fft_in[r*width_fft+c][0]+fft_in[r*width_fft+c][1]*fft_in[r*width_fft+c][1]) 
                                            / hann.data[r];
    /* ifimg
    fimg = circshift(fimg, [0 -round(freq)]);
    fimg(~line_mask) = 0;
    CHANGED
    +fimg = circshift(fimg, [0 -round(freq)]);
    ...
    ifimg = ifft(fimg,[],dim_ana);
    % cropping in the inverse Fourier transformed image
    ifimg = ifimg(1:size(img,1), 1:size(img,2));

    out.ifimg = ifimg;
    */

    for (r=0;r<height_fft;++r)
        for (c=0;c<width_fft;++c)
        {
            //fft_out_amp[c+r*width_fft][0] = fft_out[(c+    +width_fft)%width_fft+((r+freq+height_fft)%height_fft)*width_fft][0] * (r<=filter_width || r>=height_fft - filter_width);
            //fft_out_amp[c+r*width_fft][1] = fft_out[(c+    +width_fft)%width_fft+((r+freq+height_fft)%height_fft)*width_fft][1] * (r<=filter_width || r>=height_fft - filter_width);
            //fft_out_amp[c+r*width_fft][0] = fft_out[(c+freq+width_fft)%width_fft+((r+    +height_fft)%height_fft)*width_fft][0] * (c<filter_width+1 || c>(width_fft - filter_width-2));
            //fft_out_amp[c+r*width_fft][1] = fft_out[(c+freq+width_fft)%width_fft+((r+    +height_fft)%height_fft)*width_fft][1] * (c<filter_width+1 || c>(width_fft - filter_width-2));
            fft_out_amp[c+r*width_fft][0] = fft_out[c+r*width_fft][0] * (c+(filter_width+2)>freq && c-(filter_width+1)<freq);
            fft_out_amp[c+r*width_fft][1] = fft_out[c+r*width_fft][1] * (c+(filter_width+2)>freq && c-(filter_width+1)<freq);
        }
    FFTEXECA(fft_plan_bw,fft_out_amp,fft_in);
    put_out->ifimg_re.resize(width,height);
    put_out->ifimg_im.resize(width,height);
    //put_out->ifimg_re.resize(width_fft,height);
    //put_out->ifimg_im.resize(width_fft,height);
    for (r=0;r<height;++r)
        for (c=0;c<width;++c)
        {
            put_out->ifimg_re.data[c+r*width] = fft_in[r*width_fft+c][0]/width_fft;// / hann.data[r];
            put_out->ifimg_im.data[c+r*width] = fft_in[r*width_fft+c][1]/width_fft;// / hann.data[r];;
            //put_out->ifimg_re.data[c+r*width_fft] = fft_out_amp[r*width_fft+c][0];// / hann.data[r];
            //put_out->ifimg_im.data[c+r*width_fft] = fft_out_amp[r*width_fft+c][1];// / hann.data[r];;
        }
    
    /* phase
    dphi_wrapped = angle(ifimg); %this is atan -> values netween -pi,pi -> need to unwrap
    out.dphi_wrap = dphi_wrapped;

    dphi_unwrap = unwrap(dphi_wrapped,pi,2);
    dphi_unwrap = unwrap(dphi_unwrap,pi,1);
    %dphi_unwrap = double(Miguel_2D_unwrapper(single(dphi_wrapped)));
    out.dphi_unwrap = dphi_unwrap;
    */
    put_out->dphi_unwrap.resize(width,height);
    put_out->dphi_wrap.resize(width,height);
    for (r=0;r<height;++r)
        for (c=0;c<width;++c)
        {
            put_out->dphi_wrap.data[c+r*width] = ATAN2(fft_in[r*width_fft+c][1],fft_in[r*width_fft+c][0]);
        }
    put_out->dphi_wrap.unwrap(&put_out->dphi_unwrap);
    /*
    [ Y, X ] = ndgrid(size(img,1):-1:1, 1:size(img,2));

    % fit a first order polynom through the unwrapped phase 
    % which is due to the fringe fundamental frequency, 
    %doing this is equvialent to mTakeda1984's shifting in fourier space

    polys = zeros(2,size(img,1));
    for ii = 1:size(img,1)
        poly_tmp = polyfit(X(ii,:),dphi_unwrap(ii,:),1);
        polys(1,ii) = poly_tmp(1);
        polys(2,ii) = poly_tmp(2);
    end
    CHANGED
    -dphi = dphi_unwrap - mean(polys(1,:))*X-mean(polys(2,:));
    +%dphi = dphi_unwrap - mean(polys(1,:))*X-mean(polys(2,:));
    %do nothing
    */
    // remove average horizontal slope
    put_out->dphi.resize(width,height); 
    /* do nothing
    Sx = (width * (width+1))/2;
    Sxx = (width * (width+1) * (2*width + 1))/6;
    
    fit_a = fit_b = 0;
    for (r=0;r<height;++r)
    {
        Sy = 0;
        Sxy = 0;
        for(c=0;c<width;++c)
        {
            Sy += put_out->dphi_unwrap.data[r*width+c];
            Sxy += (c+1)*put_out->dphi_unwrap.data[r*width+c];
        }
        fit_a_tmp = (width * Sxy - Sx * Sy ) / (width * Sxx - Sx*Sx );
        fit_a += fit_a_tmp;
        fit_b += (Sy - fit_a_tmp * Sx) / width;
    }
    fit_a = fit_a / height;
    fit_b = fit_b / height;

    for (r=0;r<height;++r)
        for(c=0;c<width;++c)
            put_out->dphi.data[r*width+c] = put_out->dphi_unwrap.data[r*width+c] - fit_a*(c+1) - fit_b;
    */
    for (r=0;r<height;++r)
        for(c=0;c<width;++c)
            put_out->dphi.data[r*width+c] = put_out->dphi_unwrap.data[r*width+c];
    
    /*
    CHANGED
    -polys = zeros(1,size(img,2)*ip.ncrop);
    -for ii = ip.ncrop:size(img,2)-ip.ncrop
    -    poly_tmp = polyfit(Y(:,ii),dphi(:,ii),1);
    -    polys(1,ii-ip.ncrop+1) = poly_tmp(1);
    +polys = zeros(1,size(alpha,1)-2*ip.ncrop+2);
    +for ii = ip.ncrop:size(alpha,1)-ip.ncrop+1
    +   poly_tmp = polyfit(X(ii,:),alpha(ii,:),1);
    +   polys(1,ii-ip.ncrop+1) = poly_tmp(1);
    end
    */
    //vertical slope on central data
    //no effect
    /*
    Sx = (height * (height-1))/2;
    Sxx = (height * (height-1) * (2*height - 1))/6;
    tmp_fit_y.resize(width - 2*ip->ncrop,height);
    put_out->polys.resize(2,width - 2*ip->ncrop);
    for(c=0;c<tmp_fit_y.width;++c)
    {
        Sy = 0;
        Sxy = 0;
        for(r=0;r<height;++r)
        {
            Sy += tmp_fit_x.data[r*width+c+(ip->ncrop)];
            Sxy += (height-1-r)*tmp_fit_x.data[r*width+c+(ip->ncrop)];
        }
        
        put_out->polys.data[2*c  ] = (height * Sxy - Sx * Sy ) / (height * Sxx - Sx*Sx );
        put_out->polys.data[2*c+1] = (Sy - put_out->polys.data[2*c]  * Sx) / height;
        
        for(r=0;r<height;++r)
            tmp_fit_y.data[r*(width - 2*ip->ncrop)+c] = tmp_fit_x.data[r*width+c]
                     - (height-1-r) * (put_out->polys.data[2*(c-ip->ncrop)  ])
                     - put_out->polys.data[2*(c-ip->ncrop)+1];
        
    }
    if (ip->verbose > 9000)
        tmp_fit_y.save_tiff("fit_y.tif");
    */
    /*
    % alpha wavefront propagation angle; Eq. 2.74 Rutishauser thesis(calibration params)
    alpha = (ip.M0*cos(ip.beta_fix)/cos(ip.beta_rot) - 1)*Y*ip.pixel_size/ip.d + ip.M0*ip.p2*cos(ip.beta_fix)/(2*pi*ip.d*(cos(ip.beta_rot))^2) * dphi;
    % first order polynomial fit of the propagation angle of the wavefront through the image
    CHANGED
    -polys = zeros(1,size(alpha,2)-2*ip.ncrop+2);
    -for ii = ip.ncrop:size(alpha,2)-ip.ncrop+1
    -    poly_tmp = polyfit(Y(:,ii),alpha(:,ii),1);
    -    polys(1,ii-ip.ncrop+1) = poly_tmp(1);
    +polys = zeros(1,size(alpha,1)-2*ip.ncrop+2);
    +for ii = ip.ncrop:size(alpha,1)-ip.ncrop+1
    +   poly_tmp = polyfit(X(ii,:),alpha(ii,:),1);
    +   polys(1,ii-ip.ncrop+1) = poly_tmp(1);
    end  
    out.polys = polys;
    out.alpha2D = alpha;
    */
    put_out->polys.resize(2,height +2- 2*ip->ncrop);
    put_out->alpha2D.resize(width,height);
    alpha_mul =(magnification * COS(ip->beta_fix) / COS(ip->beta_rot) -1 ) * (ip->pixel_size)/(ip->distance);
    alpha_add = magnification * (ip->period_g2) * COS(ip->beta_fix) / (2*C_PI * (ip->distance) * COS(ip->beta_rot)* COS(ip->beta_rot));
    //alpha[c,r] = alpha_mul * r + alpha_add * dphi[c,r]
    Sx = (width * (width+1))/2;
    Sxx = (width * (width+1) * (2*width + 1))/6;
    for(c=0;c<width;++c)
        for(r=0;r<height;++r)
            put_out->alpha2D.data[c+r*width] = alpha_mul * (height-r) + alpha_add * put_out->dphi.data[r*width+c];
    polys_avg = 0;
    for(r=0;r<put_out->polys.height;++r)
    {
        Sy = 0;
        Sxy = 0;
        for(c=0;c<width;++c)
        {
            alpha = put_out->alpha2D.data[(r+(ip->ncrop-1))*width+c+0];
            Sy += alpha;
            Sxy += (c+1)*alpha;
        }
        
        put_out->polys.data[2*r  ] = (width * Sxy - Sx * Sy ) / (width * Sxx - Sx*Sx );
        put_out->polys.data[2*r+1] = (Sy - put_out->polys.data[2*c]  * Sx) / width;
        polys_avg += put_out->polys.data[2*r  ];
    }
    polys_avg /= put_out->polys.height;
    polys_std = 0;
    for(r=0;r<put_out->polys.height;++r)
    {
        val = put_out->polys.data[2*r  ] - polys_avg;
        polys_std += val*val;
    }
    polys_std = SQRT(polys_std / (put_out->polys.height - 1));

    /* ROC

    %% radius of curvature ROC (connected through a geometric consideration to the propagation angle of the wavefront) ; Eq. 2.26
    ROC = ip.pixel_size/abs(mean(polys));
    dROC = ip.pixel_size/(abs(mean(polys))^2)*std(polys);
    %
    fprintf(2, 'ROC %.3f +/- %.3f m \n', ROC, dROC);
    out.ROC=ROC;
    out.dROC=dROC;
    */
    put_out->roc = ip->pixel_size / polys_avg;
    put_out->droc =  ip->pixel_size * polys_std / (polys_avg *polys_avg);
    if (ip->verbose > 1337)
        printf("[I]ROC %e +- %e\n",put_out->roc,put_out->droc);
    
    /*center of mass the images
    CHANGED
    +X_hist=sum(img,1); X=1:size(img,2); 
    +X_hist=sum(img,2); X=1:size(img,1); 
    i_x=round(sum(X.*X_hist)/sum(X_hist));

    laxis = (1:size(img,1))*ip.pixel_size;
    out.laxis=laxis;
    CHANGED
    -slope_data = alpha(:, i_x)';   
    +slope_data = alpha(i_x,:)';   

    slope_data = slope_data-mean(slope_data);
    out.wpa=slope_data;

    residual_slope_data=detrend(slope_data,'linear');
    wf=residual_slope_data;
    wf_rms = rms(residual_slope_data*1e6);
    out.wpa = wf;
    */
    sum = 0;
    val = 0;
    put_out->laxis.resize(width,1);
    for(c=0;c<width;++c)
    {
        put_out->laxis.data[c] = (c+1) * ip->pixel_size;
        for(r=0;r<height;++r)
        {
            val += r * put_out->roi.data[r*width+c];
            sum +=     put_out->roi.data[r*width+c];
        }
    }
    i_x = ROUND(val / sum);
    //slope data
    put_out->alpha.resize(1,height);
    put_out->wpa.resize(1,width);
    sum = 0;
    for(c=0;c<width;++c)
        sum += put_out->alpha2D.data[i_x*(put_out->alpha2D.width)+c];
    sum /= width;
    for(c=0;c<width;++c)
        put_out->wpa.data[c] = put_out->alpha2D.data[i_x*(put_out->alpha2D.width)+c] - sum;
    /*%% height profile
    dphi_quant = ip.p2 / (2*pi * ip.d);
    % scale for the image [m]
    xaxis = (1:size(img,2))*ip.pixel_size;
    % scale for the image [m]
    yaxis = (1:size(img,1))*ip.pixel_size;

    slope_data_y=slope_data';
    % creating a rectangular grid with same dimensions as the data for the propagation angles
    [ Yoffset_y, Xoffset_y ] = ndgrid(1:size(slope_data_y,1), 1:size(slope_data_y,2));
    nbin=1;
    [ Yoffset_y_all, Xoffset_y_all ] = ndgrid(1:size(slope_data_y,1), 1:size(slope_data_y,2));
    phi_y_quant = ip.p2 / (ip.lambda * ip.d) * ip.pixel_size;
    phi_y_surface_profile_quant = phi_y_quant * ip.lambda/(2*pi); %lambda is in meters
    phi_y = cumsum( slope_data_y / dphi_quant, 1);
    out.phi_y=phi_y * phi_y_surface_profile_quant;
    %out.phi_y=cumsum( slope_data_y, 1);
    % second order polynomial fit (spherical component)
    py_phi = polyfit(xaxis', phi_y, 2);
   
    % aspherical
    %phi_y_aspherical = phi_y - repmat_size(polyval(py_phi, yaxis'), size(phi_y)); 
    phi_y_aspherical = phi_y - polyval(py_phi, xaxis');
    % metric dimension on the image
    eff_pix_size = ip.pixel_size;

    out.yaxis=yaxis;
    % height profile hp
    hp=phi_y_aspherical*phi_y_surface_profile_quant;
    hp = hp - mean(hp);
    out.hp = hp
    */
    Matrix phi_y;
    put_out->yaxis.resize(1,height);
    for(r=0;r<height;++r)
           put_out->yaxis.data[r] =  (r+1)* (ip->pixel_size);
    put_out->phi_y.zero(1,width);
    phi_y.zero(1,width);
    //phi_y = cumsum( slope_data_y / dphi_quant, 1); ... ip.p2 / (2*pi * ip.d);
    phi_y.data[0] = put_out->wpa.data[0]* 2*C_PI*(ip->distance) / (ip->period_g2);
    for(c=1;c<width;++c)
        phi_y.data[c] = phi_y.data[c-1] + put_out->wpa.data[c] * 2*C_PI*(ip->distance) / (ip->period_g2);
    Sx = Sxx = Sxxx = Sxxxx = 0;
    Sy = Sxy = Sxxy = 0;
    for(c=0;c<width;++c)
    {
        //out.phi_y=phi_y * phi_y_surface_profile_quant; .. phi_y_quant * ip.lambda/(2*pi); ... 
        put_out->phi_y.data[c]=phi_y.data[c] * (ip->period_g2) / (2*C_PI*(ip->distance)) * (ip->pixel_size);
        //x
        val = (c+1) * (ip->pixel_size);
        //y
        sum = phi_y.data[c];
        // /dphi_quant = ip.p2 / (2*pi * ip.d);
        // *phi_y_surface_profile_quant=;
        Sx   += val;
        Sy   += sum;
        Sxx  += val*val;
        Sxy  += val*sum;
        Sxxx += val*val*val;
        Sxxxx+= val*val*val*val;
        Sxxy += val*val*sum;
    }
    Sx /= width;
    Ax = Sx;
    Sxx /= width;
    Axx = Sxx;
    Sxxx /= width;
    Sy /= width;
    Ay = Sy;
    Sxy /= width;
    Sxxy /= width;
    Sxxxx /= width;
    Sxx -= Ax*Ax;
    Sxy -= Ax*Ay;
    Sxxx -= Ax*Axx;
    Sxxxx -= Axx*Axx;
    Sxxy -= Axx*Ay;
    fit_a = (Sxxy*Sxx - Sxy*Sxxx)/(Sxx*Sxxxx - Sxxx*Sxxx);
    fit_b = (Sxy*Sxxxx-Sxxy*Sxxx)/(Sxx*Sxxxx - Sxxx*Sxxx);
    fit_c = Ay - fit_b * Ax - fit_a*Axx;
    //phi_y_aspherical = phi_y - polyval(py_phi, xaxis');
    //phi_y_quant = ip.p2 / (ip.lambda * ip.d) * ip.pixel_size;
    //phi_y_surface_profile_quant = phi_y_quant * ip.lambda/(2*pi);
    //hp=phi_y_aspherical*phi_y_surface_profile_quant;
    //hp = hp - mean(hp);
    //out.hp = hp
    put_out->hp.resize(1,width);
    sum=0;
    for(c=0;c<width;++c)
    {
        //x
        val = (c+1) * (ip->pixel_size);
        //y = phi_y
        put_out->hp.data[c] = (phi_y.data[c] - fit_a*val*val -fit_b*val - fit_c)
                               * ((ip->period_g2)*(ip->pixel_size)) / (2*C_PI*(ip->distance));
        sum += put_out->hp.data[c];
    }
    sum /= width;
    for(c=0;c<width;++c)
        put_out->hp.data[c] -= sum;
    FFTPLAND(fft_plan_fw);
    FFTPLAND(fft_plan_bw);
    FFTFREE(fft_in);
    FFTFREE(fft_out);
    FFTFREE(fft_out_amp);
}

int Matrix::unwrap(Matrix *tgt)
{
    int adj,r,c;
    if (tgt==NULL)
    {
        printf("[E]Unwrap target is null");
        return -1;
    }
    if (data == NULL || width * height == 0)
    {
        printf("[E]Image to unwrap is empty\n");
        return -1;
    } 
    if (tgt->resize(width,height))
        return -1;
    Matrix mid;
    if (mid.resize(width,height))
        return -1;
    
    //STEP 1 - "unwrap" rows (only multiple wavelength compensation)
    for (r=0;r<height;++r)
    {
        adj = 0;
        mid.data[r*width] = 0;
        for (c=1;c<width;++c)
        {
            adj += (data[r*width+c  ] + C_DISCONT < data[r*width+c-1]);
            adj -= (data[r*width+c-1] + C_DISCONT < data[r*width+c]);
            mid.data[r*width+c] = adj;
        }    
    }
    /*
    //STEP 2 - "unwrap" central column
    adj = 0;
    for (r=1;r<height;++r)
    {
        adj += (data[(r  )*width+width/2] + C_DISCONT < data[(r-1)*width+width/2]);
        adj -= (data[(r-1)*width+width/2] + C_DISCONT < data[(r  )*width+width/2]);
        for (c=0;c<width;++c)
            tgt->data[r*width+c] += adj;
    }
    */
   //STEP 2 - "unwrap" columns
    for (c=0;c<width;++c)
    {
        adj = 0;
        tgt->data[c] = adj;
        for (r=1;r<height;++r)
        {
            adj += (C_DISCONT + data[(r  )*width+c] + 2*C_PI*mid.data[(r  )*width+c] 
                              < data[(r-1)*width+c] + 2*C_PI*mid.data[(r-1)*width+c]);
            adj -= (C_DISCONT + data[(r-1)*width+c] + 2*C_PI*mid.data[(r-1)*width+c] 
                              < data[(r  )*width+c] + 2*C_PI*mid.data[(r  )*width+c]);
            tgt->data[r*width+c] = adj;
        }    
    }
    //STEP 3 - make center 0
    
    adj = tgt->data[(height/2)*width+width/2] + mid.data[(height/2)*width+width/2];
    //don't make center 0
    adj = 0;
    for (c=0;c<width*height;++c)
        tgt->data[c] = data[c] + C_PI * 2 * (tgt->data[c] + mid.data[c] - adj);
    
    /*
    for (r=0;r<height;++r)
    {
        adj = tgt->data[(r)*width+width/2];
        for (c=0;c<width;++c)
        {
            tgt->data[r*width+c] = data[r*width+c] + C_PI * 2 * (tgt->data[r*width+c] - adj);
            //tgt->data[r*width+c] = tgt->data[r*width+c] - adj;
        }
    }
    */
    return 0;
}

Matrix::Matrix()
    : width(0)
    , height(0)
    , data(NULL)
{
}

Matrix::Matrix(int w,int h)
    : width(w)
    , height(h)
    , data((FLOATING*)malloc(SOF*w*h))
{
}

Matrix::~Matrix(){
    free();
}
void Matrix::free()
{
    if (data !=NULL)
    {
        ::free(data);
        data = NULL;
        width = 0;
        height = 0;
    }
}
    
int Matrix::resize(int w,int h){
    if (w==this->width && h==this->height)
        return 0;
    free();
    data = (FLOATING*)malloc(SOF*w*h);
    width = w;
    height = h;
    return 0 - (data==NULL);
}

int Matrix::zero(int w,int h){
    if (w==this->width && h==this->height)
    {
        memset(data,0,SOF*w*h);
        return 0;
    }
    free();
    data = (FLOATING*)malloc(SOF*w*h);
    width = w;
    height = h;
    if (data==NULL)
        return -1;
    memset(data,0,sizeof(data));
    return 0;
}

int Matrix::copy(Matrix *tgt)
{
    if (tgt==NULL)
    {
        printf("[E]Matrix copy target is null");
        return -1;
    }
    if (width == 0 || height ==0 || data == NULL)
    {
        printf("[E]Copy source is empty");
        return -1;
    }
    if (tgt->resize(width,height) != 0 )
        return -1;
    for (int i=0;i<(width)*(height);++i)
        tgt->data[i] = data[i];
    return 0;
}

int Matrix::add(Matrix *other)
{
    if (other==NULL)
    {
        printf("[E]Matrix add parameter is null");
        return -1;
    }
    if (width != other->width || height != other->height)
    {
        printf("[E]Matrix add size mismatch");
        return -1;
    }
    for (int i=0;i<width*height;++i)
        data[i] += other->data[i];
    return 0;
}

int Matrix::sub(Matrix *other)
{
    if (other==NULL)
    {
        printf("[E]Matrix sub parameter is null");
        return -1;
    }
    if (width != other->width || height != other->height)
    {
        printf("[E]Matrix sub size mismatch");
        return -1;
    }
    for (int i=0;i<width*height;++i)
        data[i] -= other->data[i];
    return 0;
}
int Matrix::mul(Matrix *other)
{
    if (other==NULL)
    {
        printf("[E]Matrix mul parameter is null");
        return -1;
    }
    if (width != other->width || height != other->height)
    {
        printf("[E]Matrix mul size mismatch");
        return -1;
    }
    for (int i=0;i<width*height;++i)
        data[i] *= other->data[i];
    return 0;
}

int Matrix::mul(FLOATING factor)
{
    for (int i=0;i<width*height;++i)
        data[i] *= factor;
    return 0;
}

int Matrix::add(Matrix *a, Matrix *b, Matrix *c)
{
    if (a == NULL || b == NULL || c==NULL)
    {
        printf("[E]Matrix add parameter is null");
        return -1;
    }
    if (a->width != b->width || a->height != b->height)
    {
        printf("[E]Matrix add size mismatch");
        return -1;
    }
    if (c->resize(a->width,a->height) != 0 )
        return -1;
    for (int i=0;i<(a->width)*(a->height);++i)
        c->data[i] = a->data[i] + b->data[i];
    return 0;
}

int Matrix::sub(Matrix *a, Matrix *b, Matrix *c)
{
    if (a == NULL || b == NULL || c==NULL)
    {
        printf("[E]Matrix sub parameter is null");
        return -1;
    }
    if (a->width != b->width || a->height != b->height)
    {
        printf("[E]Matrix sub size mismatch");
        return -1;
    }
    if (c->resize(a->width,a->height) != 0 )
        return -1;
    for (int i=0;i<(a->width)*(a->height);++i)
        c->data[i] = a->data[i] - b->data[i];
    return 0;
}

int Matrix::mul(Matrix *a, Matrix *b, Matrix *c)
{
    if (a == NULL || b == NULL || c==NULL)
    {
        printf("[E]Matrix mul parameter is null");
        return -1;
    }
    if (a->width != b->width || a->height != b->height)
    {
        printf("[E]Matrix mul size mismatch");
        return -1;
    }
    if (c->resize(a->width,a->height) != 0 )
        return -1;
    for (int i=0;i<(a->width)*(a->height);++i)
        c->data[i] = a->data[i] * b->data[i];
    return 0;
}

int Matrix::mul(Matrix *a, FLOATING b, Matrix *c)
{
    if (a == NULL || c==NULL)
    {
        printf("[E]Matrix mul parameter is null");
        return -1;
    }
    if (c->resize(a->width,a->height) != 0 )
        return -1;
    for (int i=0;i<(a->width)*(a->height);++i)
        c->data[i] = a->data[i] * b;
    return 0;
}

int Matrix::crop(Matrix *tgt, Rect roi, Orientation orientation)
{
    int rv=0;
    int r,c,tr,tc;
    if (tgt==NULL)
    {
        printf("[E]Crop target is null");
        rv = -1;
    }
    if (data == NULL || width * height == 0)
    {
        printf("[E]Image to crop is empty\n");
        rv = -1;
    }
    if (roi.x_lo >= roi.x_hi || roi.y_lo >= roi.y_hi)
    {
        printf("[E]Invalid ROI\n");
        rv = -1;
    }
    if (rv)
        return rv;
    if (orientation.tpose)
        rv = tgt->resize(roi.y_hi - roi.y_lo,roi.x_hi - roi.x_lo);
    else
        rv = tgt->resize(roi.x_hi - roi.x_lo,roi.y_hi - roi.y_lo);
    if (rv){
        printf("[E]Crop target resize error\n");
        return rv;
    }
    /*
    for (r = roi.y_lo;r<roi.y_hi;++r)
        for (c = roi.x_lo;c<roi.x_hi;++c)
        {
            //translate tgt coordinates [c,r] into my coordinates [tc,tr]
            if (orientation.tpose)
            {
                if (orientation.fliplr)
                    tc = (tgt->width) - 1 - (r - roi.y_lo);
                else
                    tc = (r - roi.y_lo);
                if (orientation.flipud)
                    tr = (tgt->height) - 1 - (c - roi.x_lo);
                else
                    tr = (c - roi.x_lo);
                    
            }else{
                if (orientation.fliplr)
                    tc = (tgt->width) - 1 - (c - roi.x_lo);
                else
                    tc = (c - roi.x_lo);
                if (orientation.flipud)
                    tr = (tgt->height) - 1 - (r - roi.y_lo); 
                else
                    tr = (r - roi.y_lo);
            }
            tgt->data[tc+tr*(tgt->width)] = data[c+r*width];
        }
    */
    
    for(tr = 0;tr<tgt->height;++tr)
        for(tc = 0;tc<tgt->width;++tc)
        {
            //translate tgt coordinates [tc,tr] into my coordinates [c,r]
            if (orientation.tpose)
            {
                if (orientation.fliplr)
                    r = height - 1 - roi.x_lo - tc;
                else
                    r = tc + roi.x_lo;
                if (orientation.flipud)
                    c = width - 1 - roi.y_lo - tr;
                else
                    c = tr + roi.y_lo;
            }else{
                if (orientation.fliplr)
                    c = width - 1 - roi.x_lo - tc;
                else
                    c = tc + roi.x_lo;
                if (orientation.flipud)
                    r = height - 1 - roi.y_lo - tr; 
                else
                    r = tr + roi.y_lo;
            }
            if (r<0 || c<0 || r>= height || c>= width)
                tgt->data[tc+tr*(tgt->width)] = 0;
            else
                tgt->data[tc+tr*(tgt->width)] = data[c+r*width];
        }
    
    return rv;
}

#ifdef HAVE_TIFF
int Matrix::save_tiff(const char *fn)
{
    if (data == NULL || width * height == 0)
    {
        printf("[E]Image to save is empty\n");
        return -1;
    }
    TIFF* tif = TIFFOpen(fn, "w");
    if (tif==NULL)
    {
        printf("[E]Could not write file %s\n",fn);
        return -1;
    }
    uint w,h;
    tdata_t buf;
    w = width;
    h = height;
    TIFFSetField(tif, TIFFTAG_IMAGEWIDTH,w);
    TIFFSetField(tif, TIFFTAG_IMAGELENGTH,h);
    TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE,32);
    TIFFSetField(tif, TIFFTAG_SAMPLESPERPIXEL,1);
    TIFFSetField(tif, TIFFTAG_SAMPLEFORMAT,SAMPLEFORMAT_IEEEFP);
    buf = _TIFFmalloc(TIFFScanlineSize(tif));
    float *f_line = (float*)buf;
    for (int row = 0; row < h; ++row)
    {
        for (int i=0;i<w;++i)
            f_line[i] = data[i+row*w];
        TIFFWriteScanline(tif, buf, row);
    }
    _TIFFfree(buf);
    TIFFClose(tif);
}
#else
int Matrix::save_tiff(const char *fn)
{
    printf("[E]No TIFF support\n",fn);
    return -1;
}
#endif

#ifdef HAVE_PNG
int Matrix::save_png(const char *fn,FLOATING data_low, FLOATING data_high)
{
    //const unsigned char col[6] = {0,0,0,255,255,255}; //black to white
    const unsigned char col[6] = {255,0,0,0,0,255}; //red to blue
    FLOATING val;
    if (data == NULL || width * height == 0)
    {
        printf("[E]Image to save is empty\n");
        return -1;
    }
    int code = 0;
    FILE *fp = NULL;
    png_structp png_ptr = NULL;
    png_infop info_ptr = NULL;
    png_bytep row = NULL;
    // Open file for writing (binary mode)
    fp = fopen(fn, "wb");
    if (fp == NULL) {
        printf("[E]Could not open file %s for writing\n", fn);
        code = 1;
        goto finalise;
    }
    // Initialize write structure
    png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    if (png_ptr == NULL) {
        printf("[E]Could not allocate write struct\n");
        code = 1;
        goto finalise;
    }

    // Initialize info structure
    info_ptr = png_create_info_struct(png_ptr);
    if (info_ptr == NULL) {
        printf("[E]Could not allocate info struct\n");
        code = 1;
        goto finalise;
    }
    // Setup Exception handling
    if (setjmp(png_jmpbuf(png_ptr))) {
        printf("[E]Png creation\n");
        code = 1;
        goto finalise;
    }
    png_init_io(png_ptr, fp);

    // Write header (8 bit colour depth)
    png_set_IHDR(png_ptr, info_ptr, width, height,
            8, PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE,
            PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);
    /*
    // Set title
    if (title != NULL) {
        png_text title_text;
        title_text.compression = PNG_TEXT_COMPRESSION_NONE;
        title_text.key = "Title";
        title_text.text = title;
        png_set_text(png_ptr, info_ptr, &title_text, 1);
    }
    */

    png_write_info(png_ptr, info_ptr);

    // Allocate memory for one row (3 bytes per pixel - RGB)
    row = (png_bytep) malloc(3 * width * sizeof(png_byte));
    if (data_low == INFINITY)
    {
        data_low = data[0];
        for (uint i=1;i<width*height;++i)
            data_low = FMIN(data_low,data[i]);
    }
    if (data_high == -INFINITY)
    {
        data_high = data[0];
        for (uint i=1;i<width*height;++i)
            data_high = FMAX(data_high,data[i]);
    }

    // Write image data
    int x, y;
    for (y=0 ; y<height ; y++)
    {
        for (x=0 ; x<width ; x++)
        {
            val = (data[x+y*width] - data_low)/(data_high-data_low+C_EPSILON);
            val = FMAX(0,val);
            val = FMIN(1,val);
            row[x*3  ] = (unsigned char)(val * col[3] + (1-val)*col[0]);
            row[x*3+1] = (unsigned char)(val * col[4] + (1-val)*col[1]);
            row[x*3+2] = (unsigned char)(val * col[5] + (1-val)*col[2]);
        }
        png_write_row(png_ptr, row);
    }

    // End write
    png_write_end(png_ptr, NULL);
    finalise:
    if (fp != NULL) fclose(fp);
    if (info_ptr != NULL) png_free_data(png_ptr, info_ptr, PNG_FREE_ALL, -1);
    if (png_ptr != NULL) png_destroy_write_struct(&png_ptr, (png_infopp)NULL);
    if (row != NULL) ::free(row);

    return code;
}
#else
int Matrix::save_png(const char *fn,FLOATING data_low, FLOATING data_high)
{
    printf("[E]No PNG support\n");
}
#endif

int Matrix::save_txt(const char *fn)
{
    FILE *f = NULL;
    if (width ==0 || height==0 || data == NULL)
    {
        printf("[E]Data to save empty\n");
        return -1;
    }
    f = fopen(fn,"w");
    if (f == NULL)
    {
        printf("[E]File open error %s\n",fn);
        return -1;
    }
    fprintf(f,"%d %d\n",width,height);
    for (int c=0;c<width*height;++c)
        fprintf(f,"%02.9e%s",data[c],(c+width-1)%width?" ":"\n");
    fclose(f);
    return 0;
}

int Matrix::load_txt(const char *fn)
{
    int w,h,rv;
    FILE *f = NULL;
    f = fopen(fn,"r");
    if (f == NULL)
    {
        printf("[E]File open error %s\n",fn);
        return -1;
    }
    if (2!=fscanf(f,"%d %d\n",&w,&h))
    {
        printf("[E]File read error %s\n",fn);
        fclose(f);
        return 0;
    }
    if (resize(w,h))
    {
        printf("[E]Resize error %s\n",fn);
        fclose(f);
        return 0;
    }
    rv = 0;
    for (int c=0;c<width*height;++c)
        rv += fscanf(f,FLOATING_SF,&(data[c]));
    fclose(f);
    return - (rv != w*h); 
}
void Matrix::print()
{
    printf("%d %d\n",width,height);
    for (int c=0;c<width*height;++c)
        printf("%02.9e%c",data[c],(c+width-1)%width?' ':'\n');
}

#ifdef HAVE_HDF5
int Matrix::load_h5_ushort(const char *fn,const char *ds, int slice)
{
    hid_t       file_id, dataset_id, space_id;  /* identifiers */
    herr_t      status;
    int         w,h, ndims;
    hsize_t     dim[4];
    ushort *data_h5;

    /* Open an existing file. */
    file_id = H5Fopen(fn, H5F_ACC_RDONLY, H5P_DEFAULT);

    /* Open an existing dataset. */
    dataset_id = H5Dopen(file_id, ds, H5P_DEFAULT);
    space_id = H5Dget_space(dataset_id);
    ndims = H5Sget_simple_extent_ndims(space_id);
    if (ndims<1)
    {
        printf("[E]H5 data dims error\n");
        status = H5Dclose(dataset_id);
        status = H5Fclose(file_id);
        return -1;
    }
    H5Sget_simple_extent_dims(space_id, dim, NULL);
    if (ndims==1){
        w = dim[0];
        h = 1;
    }else{
        w = dim[ndims-1];
        h = dim[ndims-2];
    }
    data_h5 = (ushort*)malloc(2*w*h);
    
    status = H5Dread (dataset_id, H5T_NATIVE_USHORT, H5S_ALL, H5S_ALL,
                      H5P_DEFAULT, data_h5);
    status = H5Dclose(dataset_id);
    resize(w,h);
    for (int r=0;r<h;++r)
        for (int c=0;c<w;++c)
        {
            data[r*w+c]=data_h5[r*w+c];
        }
        
    ::free(data_h5);
    /* Close the file. */
    status = H5Fclose(file_id);
    return 0;
}
#else
int Matrix::load_h5_ushort(const char *fn,const char *ds, int slice)
{
    printf("[E]No HDF5 support\n");
}
#endif
