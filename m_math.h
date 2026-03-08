#ifndef M_MATH
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#define M_MATH
double EPSILON = 0.00001;
double PI = 3.14159265359;
double LN2 = 0.69314718056;
double LN10 = 2.30258509299;

//---------------------Bits---------------------
//reinventing the wheel 101

typedef struct {
    u_int64_t b;
    int size;
}m_bits;

m_bits m_bits_from(char* cbits) {
    u_int64_t bits = 0;
    char c = cbits[0];
    int size;
    for(size=0;c!='\0';size++){
        switch(c) {
            case '1': bits |= 1ULL<<size; break;
            case '0': break;
            default: exit(1);
        }
        c = cbits[size];
    }
    return (m_bits){.b=bits,.size=size};
}

char* m_bits_str(m_bits bits) {
    char* cbits = malloc((bits.size+1)*sizeof(char));
    for(int i=0;i<bits.size;i++) {
        int b = bits.b << (bits.size-i-1) >> (bits.size-1);
        cbits[i] = b?'1':'0';
    }
    cbits[bits.size] = '\0';
    return cbits;
}

//---------------------Double---------------------

//double-floating point number bitwise magic

int m_double_s(double x) {
    u_int64_t m = *(u_int64_t*)&x;
    return m>>63<<63?1:-1;
}

long long m_double_mraw(double x) {
    u_int64_t m = *(u_int64_t*)&x;
    return ((m<<12)>>12);
}

double m_ipow(double,int);

double m_double_m(double x) {
    long long m = m_double_mraw(x);
    return ((double)(m | (1ULL<<52)))/m_ipow(2,52);
}

int m_double_e(double x) {
    u_int64_t e = *(u_int64_t*)&x;
    return (((e<<1) >> 53))-1023;
}

double* m_double_parts(double x) {
    int i = (int)x;
    double* parts = malloc(2*sizeof(double));
    parts[0] = (double)i;
    parts[1] = x-(double)i;
    return parts;
}

double m_ln(double);
double m_pow(double,double);

double* m_double_logRep(double x) {
    double m = m_double_m(x);
    double e = m_double_e(x);
    double* parts = m_double_parts(e*LN2/LN10);
    double n = m*m_pow(10,parts[1]);
    return (double[]){n,parts[0]};
}

//---------------------General---------------------

void m_desmos_vecPrint(double x,double y) {
    printf("(%lf,%lf),",x,y);
}

void m_desmos_print(double* a[2],int size) {
    for(int i=0;i<size;i++) {
        printf("(%lf,%lf)",a[0][i],a[1][i]);
        if (i!=size-1) printf(",");
    }
    printf("\n");
}

void m_desmos_fprint(FILE *fptr,double* a[2],int size) {
    for(int i=0;i<size;i++){
        fprintf(fptr,"(%lf,%lf)",a[0][i],a[1][i]);
        if (i!=size-1) fprintf(fptr,",");
    }
}

void m_desmos_print_fn(double fn(double),double xmin,double xmax,double xi) {
    for(double x=xmin;x<=xmax;x+=xi){
        printf("(%lf,%lf)",x,fn(x));
        if (x<=xmax-xi) printf(",");
    }
    printf("\n");
}

double m_max2(double x,double y) {
    return x>y?x:y;
}
double m_min2(double x,double y) {
    return y>x?x:y;
}

double m_abs(double);

void m_desmos_print_fn_const(double fn1(double),double fn2(double),double xmin,double xmax,double xi) {
    double x=xmin;
    double t=0;
    while(x<=xmax){
        printf("(%lf,%lf)",x,fn1(x));
        t=m_min2(m_abs(xi/fn2(x)),xi);
        x+=t;
        if (x<=xmax) printf(",");
    }
    printf("\n");
}

void m_desmos_fprint_fn(FILE *fptr,double fn(double),double xmin,double xmax,double xi) {
    for(double x=xmin;x<=xmax;x+=xi){
        fprintf(fptr,"(%lf,%lf)",x,fn(x));
        if (x<=xmax-xi) fprintf(fptr,",");
    }
}

double m_abs(double x) {
    u_int64_t y = *(u_int64_t*)&x;
    y &= ~(1ULL<<63);
    return *(double*)&y;
}

int m_sign(double x) {
    return (x<0)*-1+(x>0)+(x==0);
}

double m_ifact(int n) {
    double y = 1;
    for (double i=1;i<=n;i++) {
        y*=i;
    }
    return y;
}

double m_ipow(double x,int xn) {
    double y=1;
    int n = m_abs(xn);
    for (int i=0;i<n;i++) {
        y*=x;
    }
    return xn<0 ? 1/y : y;
}
 
double m_exp(double x) {
    if (x==0) return 1;
    double y = 0;
    double sum = 1;
    double ax = m_abs(x);
    for (int n=0;sum>EPSILON;n++) {
        sum = m_ipow(ax,n)/(m_ifact(n));
        y+=sum;
    }
    return x<0 ? 1/y : y;
}

double m_lnslow(double x) {
    assert(x>0);
    double y = 0;
    double sum = 1;
    for (int n=0;m_abs(sum)>EPSILON;n++) {
        sum = m_ipow(((x-1)/(x+1)),2*n+1)/(2*n+1);
        y += sum;
    }
    return 2*y;
}

double m_ln(double x) {
    assert(x>0);
    double y = 0;
    return m_lnslow(m_double_m(x))+LN2*m_double_e(x);
}

double m_pow(double x,double y) {
    if (m_abs(y)<EPSILON)return 1.0;
    double ax = m_abs(x);
    double z = m_exp(y*m_ln(ax));
    return x<0 ? 1/z : z;
}

double m_sqrt(double x,int n) {
    if (x==0) return 0;
    assert(x>0);
    return m_exp(m_ln(x)/(double)n);
}

double* m_quadratic_solve(double a,double b,double c) {
    double delta = b*b-a*c;
    assert(delta>0);
    return (double[2]){(-b+m_sqrt(delta,2))/(2*a),(-b-m_sqrt(delta,2))/(2*a)};
}

double m_weirdMod(double x,double ay) {
    double y = m_abs(ay);
    double z = x;
    while (z > y)
        z -=y;
    while (z < -y)
        z +=y;
    if (z==y)
        z=0;
    return z;
}

double m_cos(double ax) {
    double y = 0;
    double sum = 1;
    double x = m_weirdMod(ax,2*PI);
    for (int n=0;m_abs(sum)>EPSILON;n++) {
        sum = (m_ipow(x,2*n)*m_ipow(-1,n))/m_ifact(2*n);
        y+=sum;
    }
    return y;
}

double m_sin(double x) {
    return m_cos(x-(PI/2.0));
}

double m_tan(double x) {
    return m_sin(x)/m_cos(x);
}

double m_atan(double x) {
    double y = 0;
    double sum = 1;
    if (-1 <= x && x <= 1) {
        for (double n=0;m_abs(sum)>EPSILON;n++) {
            sum=(
                m_ipow(-1,n)*m_ipow(x,(2*n)+1)
                )/(
                (2*n)+1
                );
            y+=sum;
        }
    }else {
        for (double n=0;m_abs(sum)>EPSILON;n++) {
            sum=(
                m_ipow(-1,n)
                )/(
                ((2*n)+1)*m_ipow(x,(2*n)+1)
                );
            y+=sum;
        }
        if (x>1) y = (PI/2)-y;
        else y = (-PI/2)-y;
    }
    return y;
}

double m_asin(double x) {
    double y = 0;
    double sum = 1;
    for (int n=0;m_abs(sum)>EPSILON;n++) {
        sum = m_ifact(2*n)/((m_ipow(2,2*n)*m_ipow(m_ifact(n),2)))*((m_ipow(x,(2*n)+1))/(((2*n)+1)));
        y+=sum;
    }
    return y;
}


double m_acos(double x) {
    return (double)(PI/2.0)-m_asin(x);
}

double m_newtonMethod(double x0,double fn(double),double df(double)) {
    double x = x0;
    while (m_abs(fn(x))>EPSILON) {
        printf("(%lf,%lf),",x,fn(x));
        x -= fn(x)/(df(x)+EPSILON);
    }
    return x;
}

double m_max(double* xn,int size) {
    double x = xn[0];
    for(int i=1;i<size;i++) {
        x = (xn[i]>x)?xn[i]:x;
    }
    return x;
}

double m_min(double* xn,int size) {
    double x = xn[0];
    for(int i=1;i<size;i++) {
        x = (xn[i]<x)?xn[i]:x;
    }
    return x;
}

//finds a simple approximation of f(x) knowing 2 points and their image with x being in between these 2 values
double m_interpolate(double x,double a,double fa,double b,double fb) {
    return (m_abs(x-b)*fa+m_abs(x-a)*fb)/m_abs(a-b);
}

typedef struct {
    double x;
    double y;
}m_point;

double m_approxfn(double x,m_point* points,int size) {
    int i;
    for(i=0;x>points[i].x&&i<size;i++);
    return m_interpolate(x,points[i].x,points[i].y,points[i-1].x,points[i-1].y);
}

double m_lerp(double a,double b,double t) {
    return (1-t)*a+t*b;
}

double* m_spline(double t,double *p1,double* p2,double* joint1,double* joint2 ) {
    double *p = malloc(2*sizeof(double));
    for(int i=0;i<2;i++)
        p[i] = m_lerp(m_lerp(m_lerp(p1[i],joint1[i],t),m_lerp(joint1[i],joint2[i],t),t),m_lerp(joint2[i],p2[i],t),t);
    return p;
}

double* m_vspline(double t,double *p1,double* p2,double* v1,double* v2) {
    double *p = malloc(2*sizeof(double));
    for(int i=0;i<2;i++) {
        double d = p1[i];
        double c = v1[i];
        double b = 3*(p2[i]-p1[i]) - 2*v1[i] - v2[i];
        double a = 2*(p1[i]-p2[i]) + v1[i] + v2[i];
        p[i] = a*t*t*t + b*t*t + c*t + d;
    }
    return p;
}

//---------------------Matrix---------------------

//the information of the matrix is stored in an array of doubles
//to get data at any line and colon i,j do this:
//  <m_mat>.data[j+i*<m_mat>.cols]
// or use the function m_mat_get(<m_mat>,i,j)
typedef struct {
    double* data;
    unsigned int rows;
    unsigned int cols;
}m_mat;


m_mat m_mat_new(unsigned int rows,unsigned int cols) {
    m_mat rtn;
    rtn.data = (double*) malloc(rows*cols*sizeof(double));
    rtn.rows = rows;
    rtn.cols = cols;
    return rtn;
}

m_mat m_mat_zero(unsigned int rows,unsigned int cols) {
    m_mat rtn = m_mat_new(rows,cols);
    for(int i=0;i<rows*cols;i++) {
        rtn.data[i]=0;
    }
    return rtn;
}

double m_mat_get(m_mat A,unsigned int i,unsigned int j) {
    return A.data[j+i*A.cols];
}

void m_mat_print(m_mat A) {
    for(int i=0;i<A.rows;i++) {
        printf("[");
        for(int j=0;j<A.cols;j++) {
            printf("%lf,",A.data[j+i*A.cols]);
        }
        printf("]\n");
    }
    printf("\n");
}

m_mat m_mat_from(double *data,unsigned int rows,unsigned int cols) {
    m_mat rtn = m_mat_new(rows,cols);
    for (int i=0;i<rows*cols;i++) {
        rtn.data[i] = data[i];
    }
    return rtn;
}

m_mat m_mat_add(m_mat A,m_mat B) {
    assert(A.rows==B.rows && A.cols==B.cols);
    m_mat rtn = m_mat_new(A.rows,A.cols);
    for (int i=0;i<rtn.cols*rtn.rows;i++) {
        rtn.data[i]=A.data[i]+B.data[i];
    }
    return rtn;
}

m_mat m_mat_sub(m_mat A, m_mat B) {
    assert(A.rows==B.rows && A.cols==B.cols);
    m_mat rtn = m_mat_new(A.rows,A.cols);
    for (int i=0;i<rtn.rows*rtn.cols;i++) {
        rtn.data[i]=A.data[i]-B.data[i];
    }
    return rtn;
}

m_mat m_mat_nmul(m_mat A,double n) {
    m_mat B = m_mat_new(A.rows,A.cols);
    for (int i=0;i<A.cols*A.rows;i++) {
        B.data[i] = A.data[i]*n;
    }
    return B;
}

m_mat m_mat_ndiv(m_mat A,double n) {
    m_mat B = m_mat_new(A.rows,A.cols);
    for (int i=0;i<A.cols*A.rows;i++) {
        B.data[i] = A.data[i]/n;
    }
    return B;
}

m_mat m_mat_map(m_mat A, double fn(double)) {
    m_mat rtn = m_mat_new(A.rows,A.cols);
    for (int i=0;i<rtn.rows*rtn.cols;i++) {
        rtn.data[i] = fn(rtn.data[i]);
    }
    return rtn;
}

m_mat m_mat_mul(m_mat A,m_mat B) {
    assert(A.cols == B.rows);
    m_mat rtn = m_mat_zero(A.rows,B.cols);
    for(int i=0;i<A.rows;i++) {
        for(int j=0;j<B.cols;j++) {
            for(int k=0;k<A.cols;k++) {
                rtn.data[j+i*rtn.cols] += A.data[k+i*A.cols]*B.data[j+k*B.cols];
            }
        }
    }
    return rtn;
}

void m_mat_lineAdd(m_mat *A,unsigned int l1,unsigned int l2,double k) {
    for(int i=0;i<A->cols;i++) {
        A->data[i+l1*A->cols] += A->data[i+l2*A->cols]*k;
    }
}

void m_mat_lineMul(m_mat *A,unsigned int l,double k) {
    for(int i=0;i<A->cols;i++) {
        A->data[i+l*A->cols] = A->data[i+l*A->cols]*k;
    }
}

void m_mat_swap(m_mat *A,unsigned int l1,unsigned int l2) {
    for(int i=0;i<A->cols;i++) {
        double a = A->data[i+l1*A->cols];
        A->data[i+l1*A->cols] = A->data[i+l2*A->cols];
        A->data[i+l2*A->cols] = a;
    }
}

m_mat m_mat_I(unsigned int n) {
    m_mat rtn = m_mat_zero(n,n);
    for (int i=0;i<n;i++) {
        rtn.data[i+i*n]=1;
    }
    return rtn;
}

m_mat m_mat_cpy(m_mat A) {
    return m_mat_from(A.data,A.rows,A.cols);
}

double m_mat_Tr(m_mat A) {
    assert(A.cols == A.rows);
    unsigned int n = A.cols;
    double trace = 1;
    for (int i=0;i<n;i++) {
        trace+=A.data[i+i*n];
    }
    return trace;
}

//calculates the upper diagonal matrix of A
m_mat m_mat_diagonalize(m_mat A) {
    m_mat Ac = m_mat_cpy(A);
    assert(Ac.cols==Ac.rows);
    unsigned int n = Ac.cols;
    for(int k=0;k<n;k++) {
        int i_max = k;
        int v_max = Ac.data[k+i_max*Ac.cols];
        for (int i = k+1; i < n; i++)
            if (m_abs(Ac.data[k+i*Ac.cols]) > v_max)
                v_max = Ac.data[k+i*Ac.cols], i_max = i;
        if (i_max != k) {
            m_mat_swap(&Ac,k,i_max);
        }
    }
    for(int i=0;i<n+1;i++) {
        for(int j=1;j<n;j++) {
            if (i<j) {
                double k = -Ac.data[i+j*Ac.cols]/Ac.data[i+i*Ac.cols];
                m_mat_lineAdd(&Ac,j,i,k);
                Ac.data[i+j*Ac.cols]=0;
            }
        }
    }
    return Ac;
}

//calculates the determinant of A by diagonalizing A and calulating its trace
double m_mat_det(m_mat A) {
    assert(A.cols == A.rows);
    m_mat B = m_mat_diagonalize(A);
    unsigned int n = B.cols;
    double det = 1;
    for (int i=0;i<n;i++) {
        det*=B.data[i+i*n];
    }
    return det;
}

//performs gaussian elimination to find the solution X of:
//  AX = B
m_mat m_mat_gauss(m_mat A,m_mat B) {
    m_mat Ac = m_mat_cpy(A);
    m_mat Bc = m_mat_cpy(B);
    assert(Ac.rows==Bc.rows);
    unsigned int n = Ac.cols;
    for(int k=0;k<n;k++) {
        int i_max = k;
        int v_max = Ac.data[k+i_max*Ac.cols];
        for (int i = k+1; i < n; i++)
            if (m_abs(Ac.data[k+i*Ac.cols]) > v_max)
                v_max = Ac.data[k+i*Ac.cols], i_max = i;
        if (i_max != k) {
            m_mat_swap(&Ac,k,i_max);
            m_mat_swap(&Bc,k,i_max);
        }
    }
    for(int i=0;i<n+1;i++) {
        for(int j=1;j<n;j++) {
            if (i<j) {
                double k = -Ac.data[i+j*Ac.cols]/Ac.data[i+i*Ac.cols];
                m_mat_lineAdd(&Ac,j,i,k);
                m_mat_lineAdd(&Bc,j,i,k);
                Ac.data[i+j*Ac.cols]=0;
            }
        }
    }
    double det=m_mat_Tr(Ac);
    assert(m_abs(det)>0.0001);
    for(int i=0;i<n;i++) {
        double k = 1.0/Ac.data[i+i*n];
        m_mat_lineMul(&Ac,i,k);
        m_mat_lineMul(&Bc,i,k);
        Ac.data[i+i*n] = 1;
    }
    for(int i=0;i<n;i++) {
        for(int j=0;j<i;j++) {
            double k = -Ac.data[i+j*Ac.cols];
            m_mat_lineAdd(&Ac,j,i,k);
            m_mat_lineAdd(&Bc,j,i,k);
            Ac.data[i+j*Ac.cols]=0;
        }
    }
    return Bc;
}

//calculates the inverse of A using gaussian elimination
m_mat m_mat_invert(m_mat A) {
    assert(A.cols == A.rows);
    return m_mat_gauss(A,m_mat_I(A.cols));
}

m_mat m_mat_ipow(m_mat A,int xn) {
    assert(A.cols == A.rows);
    m_mat B = m_mat_I(A.rows);
    int n = m_abs(xn);
    for (int i=0;i<n;i++) {
        B=m_mat_mul(B,A);
    }
    return xn<0 ? m_mat_invert(B) : B;
}

//uses the taylor series of exp(x) to find exp(A)
m_mat m_mat_exp(m_mat A) {
    assert(A.rows == A.cols);
    m_mat B = m_mat_zero(A.rows,A.cols);
    m_mat sum = m_mat_I(A.cols);
    for (int n=0;m_mat_Tr(sum)>EPSILON;n++) {
        sum = m_mat_ndiv(m_mat_ipow(A,n),(m_ifact(n)));
        B = m_mat_add(B,sum);
    }
    return B;
}

//---------------------Complex Numbers---------------------

typedef struct {
    double x;
    double y;
    double r;
    double a;
}m_cmp;

void m_cmp_printAlgeb(m_cmp z) {
    printf("%lf+i%lf\n",z.x,z.y);
}

void m_cmp_printEuler(m_cmp z) {
    printf("%lfe^i%lf\n",z.r,z.a);
}

m_cmp m_cmp_zero() {
    return (m_cmp){0,0,0,0};
}

//calculates |z| , z=x+iy
double m_cmp_normxy(double x,double y) {
    return m_sqrt(x*x+y*y,2);
}

double m_cmp_norm(m_cmp z) {return m_cmp_normxy(z.x,z.y);}

//calculates arg(z) , z=x+iy
double m_cmp_argxy(double x,double y) {
    return m_abs(x) < EPSILON ? m_sign(y)*PI/2 : m_abs(y) < EPSILON ? (x<0)*PI : m_atan(y/x);
}

double m_cmp_arg(m_cmp z) {return m_cmp_argxy(z.x,z.y);}

//returns a complex number from its algebric notation z=x+iy
//this function will be significantly slower then m_cmp_euler because we need to calculate the arg(z)
//if this function is too slow make EPSILON greater (less accuracy)
m_cmp m_cmp_algeb(double x,double y) {
    return (m_cmp){x,y,m_cmp_normxy(x,y),m_cmp_argxy(x,y)};
}

//returns a complex number from its euler notation re^ia
m_cmp m_cmp_euler(double r,double a) {
    return (m_cmp){r*m_cos(a),r*m_sin(a),r,m_weirdMod(a,2*PI)};
}

m_cmp m_cmp_add(m_cmp z1,m_cmp z2) {
    return m_cmp_algeb(z1.x+z2.x,z1.y+z2.y);
}

m_cmp m_cmp_sub(m_cmp z1,m_cmp z2) {
    return m_cmp_algeb(z1.x-z2.x,z1.y-z2.y);
}

m_cmp m_cmp_mul(m_cmp z1,m_cmp z2) {
    return m_cmp_euler(z1.r*z2.r,z1.a+z2.a);
}

m_cmp m_cmp_div(m_cmp z1,m_cmp z2) {
    return m_cmp_euler(z1.r/z2.r,z1.a-z2.a);
}

m_cmp m_cmp_addn(m_cmp z,double n) {
    return m_cmp_algeb(z.x+n,z.y);
}

m_cmp m_cmp_subn(m_cmp z,double n) {
    return m_cmp_algeb(z.x-n,z.y);
}

m_cmp m_cmp_muln(m_cmp z,double n) {
    return m_cmp_algeb(z.x*n,z.y*n);
}

m_cmp m_cmp_divn(m_cmp z,double n) {
    return m_cmp_euler(z.x/n,z.y/n);
}

m_cmp m_cmp_exp(m_cmp z) {
    return m_cmp_euler(m_exp(z.x),z.y);
}

m_cmp m_cmp_ln(m_cmp z) {
    return m_cmp_algeb(m_ln(z.r),z.a);
}

//returns the conjugate of z
m_cmp m_cmp_conj(m_cmp z) {
    return m_cmp_euler(z.r,-z.a);
}

m_cmp m_cmp_ipow(m_cmp z,int n) {
    return m_cmp_euler(m_ipow(z.r,n),z.a*n);
}

m_cmp m_cmp_cpy(m_cmp z) {
    return (m_cmp){z.x,z.y,z.r,z.a};
}

int m_cmp_ne(m_cmp z1,m_cmp z2) {
    return z1.x != z2.x || z1.y != z2.y;
}

int m_cmp_eq(m_cmp z1,m_cmp z2) {
    return z1.x == z2.x && z1.y == z2.y;
}

int m_cmp_bt(m_cmp z1,m_cmp z2) {
    return z1.x > z2.x && z1.y > z2.y;
}

int m_cmp_be(m_cmp z1,m_cmp z2) {
    return z1.x >= z2.x && z1.y >= z2.y;
}

int m_cmp_lt(m_cmp z1,m_cmp z2) {
    return z1.x < z2.x && z1.y < z2.y;
}

int m_cmp_le(m_cmp z1,m_cmp z2) {
    return z1.x <= z2.x && z1.y <= z2.y;
}

m_cmp* m_cmp_circle(int n,double r,double c) {
    double a = 2.0*PI/(double)n;
    m_cmp* X = (m_cmp*) malloc(n*sizeof(m_cmp));
    for (int k=0;k<n;k++) {
        X[k] = m_cmp_euler(r,a*k+c);
    }
    return X;
}

//---------------------Polynomials---------------------

typedef struct {
    double* a;
    int size;
}m_poly;

m_poly m_poly_new(double* a,int size) {
    return (m_poly){a,size};
}

void m_poly_df(m_poly p,m_poly dp) {
    for (int n = 1;n<p.size;n++) {
        dp.a[n-1] = p.a[n]*n;
    }
}

//f(a,x) = a[0] + a[1]*x + a[2]*x^2 + ...
double m_poly_eval(m_poly p,double x) {
    double y = 0;
    double* a = p.a;
    for (int n=0;n<p.size;n++) {
        y+=a[n]*m_ipow(x,n);
    }
    return y;
}


m_cmp m_poly_evalCmp(m_poly p,m_cmp x) {
    m_cmp y = m_cmp_zero();
    double* a = p.a;
    for (int n=0;n<p.size;n++) {
        y=m_cmp_add(y,m_cmp_muln(m_cmp_ipow(x,n),a[n]));
    }
    return y;
}

m_poly m_poly_from(m_mat points[2]) {
    assert(points[0].rows==points[1].rows);
    unsigned int n = points[0].rows;
    m_mat A = m_mat_zero(n,n);
    for(int i=0;i<n;i++) {
        for(int j=0;j<n;j++) {
            A.data[j+i*n] = m_ipow(points[0].data[i],j);
        }
    }
    return (m_poly) {m_mat_gauss(A,points[1]).data,n};
}

/*
m_cmp* m_poly_durandKerner(m_poly p,int maxI) {
    int n = p.size-1;
    double *absa = malloc(p.size*sizeof(double));
    for(int i=0;i<p.size;i++) {
        absa[i]=m_abs(p.a[i]);
        printf("%lf",absa[i]);
    }
    double r = m_max(absa,p.size)/10;
    m_cmp* X = m_cmp_circle(n,r,0.1);
    for(int k=0;k<maxI;k++) {
        printf("::\n");
        m_cmp *delta = malloc(n*sizeof(m_cmp));
        for(int i=0;i<n;i++) {
            m_desmos_vecPrint(X[i].x,X[i].y);
            m_cmp Q = m_cmp_algeb(1,0);
            for(int j=0;j<n;j++) {
                if (i!=j) Q = m_cmp_mul(Q,m_cmp_sub(X[i],X[j]));
            }
            delta[i] = m_cmp_div(m_poly_evalCmp(p,X[i]),Q);
            X[i] = m_cmp_sub(X[i],delta[i]);
        }
    }
    return X;
}
*/

void m_poly_print(m_poly p) {
    for (int n=0;n<p.size;n++) {
        printf("%lfx^%d",p.a[n],n);
        if (n<p.size-1) printf("+");
    }
    printf("\n");
}
#endif