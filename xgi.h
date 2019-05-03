#ifndef XGI_H
#define XGI_H

#ifdef USE_DOUBLE_PRECISION
#define FLOATING double
#define C_ONE 1.0
#define C_EPSILON 1e-11
#define C_PI 3.14159265358979323846
#define FMAX fmax
#define FMIN fmin
#define ROUND round
#define ABS fabs
#define COS cos
#define ATAN2 atan2
#define SQRT sqrt
#define FFTPLAN fftw_plan
#define FFTPLAN1D fftw_plan_dft_1d
#define FFTPLANMANY fftw_plan_many_dft
#define FFTPLAND fftw_destroy_plan
#define FFTC fftw_complex
#define FFTMALLOC fftw_malloc
#define FFTFREE fftw_free
#define FFTEXEC fftw_execute
#define FFTEXECA fftw_execute_dft
#define FFTTINIT fftw_init_threads
#define FFTTPLAN fftw_plan_with_nthreads
#define FFTTCLEANUP fftw_cleanup_threads
#define FLOATING_SF "%lf"
#else
#define FLOATING float
#define C_ONE 1.0f
#define C_EPSILON 1e-5
#define C_PI 3.1415926f
#define FMAX fmaxf
#define FMIN fminf
#define ROUND roundf
#define ABS fabsf
#define COS cosf
#define ATAN2 atan2f
#define SQRT sqrtf
#define FFTPLAN fftwf_plan
#define FFTPLAN1D fftwf_plan_dft_1d
#define FFTPLANMANY fftwf_plan_many_dft
#define FFTPLAND fftwf_destroy_plan
#define FFTC fftwf_complex
#define FFTMALLOC fftwf_malloc
#define FFTFREE fftwf_free
#define FFTEXEC fftwf_execute
#define FFTEXECA fftwf_execute_dft
#define FFTTINIT fftwf_init_threads
#define FFTTPLAN fftwf_plan_with_nthreads
#define FFTTCLEANUP fftwf_cleanup_threads
#define FLOATING_SF "%f"
#endif
#define SOF sizeof(FLOATING)
#define C_DISCONT (C_ONE*C_PI)

#include <math.h>

struct Rect
{
    int x_lo;
    int x_hi;
    int y_lo;
    int y_hi;
    Rect(int _x_lo,int _x_hi,int _y_lo, int _y_hi):x_lo(_x_lo),x_hi(_x_hi),y_lo(_y_lo),y_hi(_y_hi){};
    Rect():x_lo(0),x_hi(0),y_lo(0),y_hi(0){};
    
};

struct Orientation
{
    char tpose;
    char fliplr;
    char flipud;
    Orientation(char _tpose,char _fliplr,char _flipud):tpose(_tpose),fliplr(_fliplr),flipud(_flipud){};
};

struct Matrix
{
    int width;
    int height;
    FLOATING *data;

    Matrix();
    Matrix(int w,int h);
    ~Matrix();
    void free();
    int resize(int w,int h);
    int copy(Matrix *tgt);
    int zero(int w,int h);
    int add(Matrix *other);
    int sub(Matrix *other);
    int mul(Matrix *other);
    int mul(FLOATING factor);  
    int crop(Matrix *tgt,Rect roi, Orientation orientation = Orientation(0,0,0));
    int unwrap(Matrix *tgt);
    int save_tiff(const char *fn);
    int save_png(const char *fn,FLOATING data_low = INFINITY, FLOATING data_high = - INFINITY);
    int save_txt(const char *fn);
    int load_txt(const char *fn);
    void print();
    int load_h5_ushort(const char *fn,const char *ds, int slice=0);
    static int add(Matrix *a, Matrix *b, Matrix *c);
    static int sub(Matrix *a, Matrix *b, Matrix *c);
    static int mul(Matrix *a, Matrix *b, Matrix *c);
    static int mul(Matrix *a, FLOATING b, Matrix *c);
};

struct Xgi_ip
{
    //roi pixel positions: <lo,hi)
    Rect roi;
    Orientation cam_orientation;
    FLOATING beta_fix;
    FLOATING beta_rot;
    //period of grating 1 [m]
    FLOATING period_g1;
    //period of grating 2 [m]
    FLOATING period_g2;
    FLOATING distance;
    int talbot_order;
    FLOATING energy;
    FLOATING pixel_size;
    //zero padding, power of 2
    int n_zeropadding;
    //search region in fourier space
    int excl_para;
    //filter region in fourier space
    int filter_para;
    //crop after fft
    int ncrop;
    //motor position not used
    //FLOATING motor_pos;
    //threads for fftw
    int fftw_thread;
    //verbose option
    int verbose;
    Xgi_ip()
        :roi(505,827,190,525)
        ,cam_orientation(0,0,1)
        ,beta_fix(0.0271)
        ,beta_rot(-0.0185)
        ,period_g1(3.75e-6)
        ,period_g2(2.000e-6)
        ,distance(144e-3)
        ,talbot_order(11)
        ,energy(9000)
        ,pixel_size(2.857e-6)
        ,n_zeropadding(16)
        ,excl_para(13)
        ,filter_para(3)
        ,ncrop(50)
        //,motor_pos(18000)
        ,fftw_thread(1)
        ,verbose(0)
        {};
};


struct Xgi_out
{
    Matrix roi;
    Matrix fft;
    Matrix amp;
    Matrix fft_re;
    Matrix fft_im;
    Matrix ifimg_re;
    Matrix ifimg_im;
    //wrapped phase
    Matrix dphi_wrap;
    //unwrapped phase
    Matrix dphi_unwrap;
    //phase with removed average horizontal slope
    Matrix dphi;
    Matrix polys;
    Matrix alpha2D;
    Matrix alpha;
    Matrix laxis;
    Matrix wpa;
    Matrix phi_y;
    Matrix yaxis;
    Matrix hp;
    FLOATING roc;
    FLOATING droc;
    Xgi_out(){};
    ~Xgi_out(){};
    
};

class Xgi
{
    public:
        Xgi();
        ~Xgi();
    
        int compute(Matrix *put_in, Xgi_ip *ip, Xgi_out *put_out);
};

#endif //XGI_H
