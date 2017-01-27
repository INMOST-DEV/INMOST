// developed by Igor Konshin

#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>

#define NMAX 1024 // maximal pixel field dimention

static int f[NMAX][NMAX]; // pixel field
static int nm; // actual matrix dimention
static int nf; // actual field dimention
static double sc; // scaling indices from matrix to pixel field

// Initialize print of matrix portrait into PS file
void ps_init (int n)
{
    nm = n;
    nf = (n <= NMAX) ? n : NMAX;
    sc = 1;
    if (nm > 1 && nm > nf) sc = (double)(nf-1)/(double)(nm-1);
    std::cout<<"::: nm="<<nm<<" nf="<<nf<<" sc="<<sc<<std::endl;//db!
    for (int i=0; i<nf; i++)
        for (int j=0; j<nf; j++)
            f[i][j] = 0;
}

// Save 1 nonzero into some pixel (j and j are started from 0 like in C)
inline void ps_ij (int i, int j)
{
    if (i >= 0 && i < nm && j >= 0 && j < nm)
        f[(int)(j*sc)][nf-1-(int)(i*sc)] = 1;
}

// Write the collected matrix portrait into PS file
void ps_file (std::string file, int nbl = 1, int * ibl = NULL)
{
    double s1 = 50.0;
    double s = 500.0 / nf;

    std::fstream output(file.c_str(),std::ios::out);
    output << "%%BoundingBox: 0 0 600 600" << std::endl;
    output << "/m {moveto} def % x y" << std::endl;
    output << "/l {lineto} def % x y" << std::endl;
    output << "/s {stroke} def % x y" << std::endl;
    output << "/n {newpath} def % x y" << std::endl;
    output << "/c {closepath} def % x y" << std::endl;
    output << 0.1 << " setlinewidth" << std::endl;
    output << "n " << s1+s*.5 << " " << s1+s*.5 << " m "
           << s1+s*(nf+.5) << " " << s1+s*.5 << " l "
           << s1+s*(nf+.5) << " " << s1+s*(nf+.5) << " l "
           << s1+s*.5 << " " << s1+s*(nf+.5) << " l "
           << s1+s*.5 << " " << s1+s*.5 << " l s c" << std::endl;
    if (nbl > 1) {
        for (int i=1; i<nbl; i++) {
            double x = sc*ibl[i] + .5;
            double y = nf - (sc*ibl[i] + .5) + 1.0;
            output << "n " << s1+s*.5 << " " << s1+s*y << " m "
                   << s1+s*(nf+.5) << " " << s1+s*y << " l s c" << std::endl;
            output << "n " << s1+s*x << " " << s1+s*.5 << " m "
                   << s1+s*x << " " << s1+s*(nf+.5) << " l s c" << std::endl;
        }
    }
    double r = (log10(1.0*nf) + 1.0) / 4.0;
    r = s * ((r < 1.0) ? r : 1.0) / 2.0;
    r = (r > 1.0) ? r : 1.0;
    //r = .2;
    output << 2*r << " setlinewidth" << std::endl;
    output << "/d {moveto currentpoint" << r << " 0 360 arc fill} def % x y" << std::endl;
    output << "n" << std::endl;

    std::stringstream ps(std::ios::in | std::ios::out);
    //ps << std::scientific;
    //ps.precision(2);

    for (int i=0; i<nf; i++)
        for (int j=0; j<nf; j++)
            if (f[i][j] != 0)
                ps << s1+s*(i+1)-r << " " << s1+s*(j+1) << " m "
                   << s1+s*(i+1)+r << " " << s1+s*(j+1) << " l" << std::endl;

    output << ps.rdbuf();
    output << "s c" << std::endl;
    output << "showpage" << std::endl;
}

// Print of matrix portrait into PS file
void ps_crs (int n, int nbl, int * ibl, int * ia, int *ja, const std::string file)
{
    int ia0 = ia[0];
    int ja0 = n;
    for (int k=0; k<ia[nbl]-ia0; k++) ja0 = (ja0 < ja[k]) ? ja0 : ja[k];
    ps_init (n);
    for (int i=0; i<n; i++)
        for (int k=ia[i]-ia0; k<ia[i+1]-ia0; k++)
            ps_ij (i, ja[k]-ja0);
    ps_file (file, nbl, ibl);
}

#define USE_INMOST
#if defined(USE_INMOST)
// Print of INMOST matrix portrait into PS file
void ps_inmost (Sparse::Matrix & A, const std::string file)
{
    int n = A.Size();
    int row = A.GetFirstIndex();
    ps_init (n);
    for (Sparse::Matrix::iterator it = A.Begin(); it != A.End(); ++it) {
        for (Solver::Row::iterator jt = it->Begin(); jt != it->End(); ++jt)
            ps_ij (row,jt->first);
        row++;
    }
    ps_file (file, NULL, NULL);
}
#endif

//#define PS_MATR_SELFTEST_CRS
#if defined (PS_MATR_SELFTEST_CRS)
int main ()
{
    //       [  3 -1 -1  0 ]
    //   A = [ -1  3  0 -1 ]
    //       [ -1  0  3 -1 ]
    //       [  0 -1 -1  3 ]
    static int n=4; // nza=12;
    int   ia[5]  = { 0, 3, 6, 9, 12 }; //ia[n+1]
    int   ja[12] = { 1,2,3, 1,2,4, 1,3,4, 2,3,4 }; //ja[nza]
    int nbl=2, ibl[3]={0,2,4};
    ps_crs (n,nbl,ibl,ia,ja,"z.ps");
}
#endif

//#define PS_MATR_SELFTEST_0
#if defined (PS_MATR_SELFTEST_0)
int main ()
{
    int nbl=3, ibl[4]={0,30,60,100};
    ps_init (100);
    ps_ij (0,0);
    ps_ij (1,1);
    ps_ij (2,2);
    ps_ij (12,17);
    ps_ij (79,89);
    ps_ij (99,99);
    ps_file ("z.ps",nbl,ibl);
}
#endif
