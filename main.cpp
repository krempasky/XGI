#include <cstdio>
#include "xgi.h"
#define DPATH "."
//#define DPATH "orig"
using namespace std;
int main(int argc, char **args)
{
    Xgi x = Xgi();
    Xgi_ip ip;
    Matrix put_in;
    Xgi_out put_out;
    /* load data
    img = read_img_ProSilica(dir_path,87,'sfb_');
    dark = read_img_ProSilica(dir_path,61:62,'sfb_');
    img=img-mean(dark,3);
    */
    Matrix dark,misc;
    char fn[512];
    sprintf(fn,"%s/sfb_%04d.h5",DPATH,87);
    put_in.load_h5_ushort(fn,"/entry/instrument/detector/data",0);
    dark.zero(put_in.width,put_in.height);
    int dark_low = 61;
    int dark_high = 62;
    for (int i=dark_low;i<=dark_high;++i)
    {   
        sprintf(fn,"%s/sfb_%04d.h5",DPATH,i);
        misc.load_h5_ushort(fn,"/entry/instrument/detector/data",0);
        dark.add(&misc); 
    }
    if (dark_low <= dark_high)
        dark.mul(C_ONE/(dark_high +1 - dark_low));
    put_in.save_tiff("input.tif");
    put_in.sub(&dark);
    //ip.verbose = 9001;
    ip.n_zeropadding = 16;
    //ip.ncrop = 128;
    //put_in.data.save_tiff("testi.tif");
    //put_in.data.save_png("testi.png",INFINITY,300);
    printf("Compute ... \n");
    x.compute(&put_in,&ip,&put_out);
    printf("done\n");
    //put_out.fft_re.save_tiff("test.tif");
    put_out.roi.save_tiff("roi.tif");
    put_out.fft_re.save_tiff("fftre.tif");
    put_out.fft_im.save_tiff("fftim.tif");
    put_out.fft_re.save_txt("fftre.txt");
    put_out.fft_im.save_txt("fftim.txt");
    put_out.fft.save_tiff("fft.tif");
    put_out.amp.save_tiff("amp.tif");
    put_out.ifimg_re.save_tiff("ifftre.tif");
    put_out.ifimg_im.save_tiff("ifftim.tif");
    put_out.dphi_wrap.save_tiff("dphi_w.tif");
    put_out.dphi_unwrap.save_tiff("dphi_uw.tif");
    put_out.dphi.save_tiff("dphi.tif");
    //put_out.polys.print();
    put_out.polys.save_tiff("polys.tif");
    put_out.alpha.save_tiff("alpha.tif");
    put_out.alpha2D.save_tiff("alpha2D.tif");
    put_out.wpa.save_tiff("wpa.tif");
    put_out.phi_y.save_tiff("phi_y.tif");
    put_out.yaxis.save_tiff("xaxis.tif");
    put_out.hp.save_tiff("hp.tif");

    
}

/*
UBUNTU install libs:
libtiff5-dev,libhdf5-dev,libpng-dev,libfftw3-dev

*/
