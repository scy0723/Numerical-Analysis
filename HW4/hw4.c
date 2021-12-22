#define _USE_MATH_DEFINES
#include <math.h>
#include <stdio.h>
#include <time.h>

int iteration = 0;

void zbrak(float (*fx)(float), float x1, float x2, int n, float xb1[],
    float xb2[], int* nb)
{
    int nbb, i;
    float x, fp, fc, dx;

    nbb = 0;
    dx = (x2 - x1) / n;
    fp = (*fx)(x = x1);
    for (i = 1; i <= n; i++) {
        fc = (*fx)(x += dx);
        if (fc * fp <= 0.0) {
            xb1[++nbb] = x - dx;
            xb2[nbb] = x;
            if (*nb == nbb) return;

        }
        fp = fc;
    }
    *nb = nbb;
}


#define JMAX 40
float rtbis(float (*func)(float), float x1, float x2, float xacc)
{
    int j;
    float dx, f, fmid, xmid, rtb;

    f = (*func)(x1);
    fmid = (*func)(x2);
    rtb = f < 0.0 ? (dx = x2 - x1, x1) : (dx = x1 - x2, x2);
    for (j = 1; j <= JMAX; j++) {
        fmid = (*func)(xmid = rtb + (dx *= 0.5));
        if (fmid <= 0.0) rtb = xmid;
        if (fabs(dx) < xacc || fmid == 0.0) { 
            iteration = j - 1;
            return rtb; 
        }
    }
    
    return 0.0;
}
#undef JMAX

#define MAXIT 30
float rtflsp(float (*func)(float), float x1, float x2, float xacc)
{
    int j;
    float fl, fh, xl, xh, swap, dx, del, f, rtf;

    fl = (*func)(x1);
    fh = (*func)(x2);
    if (fl < 0.0) {
        xl = x1;
        xh = x2;
    }
    else {
        xl = x2;
        xh = x1;
        swap = fl;
        fl = fh;
        fh = swap;
    }
    dx = xh - xl;
    for (j = 1; j <= MAXIT; j++) {
        rtf = xl + dx * fl / (fl - fh);
        f = (*func)(rtf);
        if (f < 0.0) {
            del = xl - rtf;
            xl = rtf;
            fl = f;
        }
        else {
            del = xh - rtf;
            xh = rtf;
            fh = f;
        }
        dx = xh - xl;
        if (fabs(del) < xacc || f == 0.0) { iteration = j - 1;  return rtf; }
    }
    return 0.0;
}

#undef MAXIT

#define MAXIT 30
float rtsec(float (*func)(float), float x1, float x2, float xacc)
{
    int j;
    float fl, f, dx, swap, xl, rts;

    fl = (*func)(x1);
    f = (*func)(x2);
    if (fabs(fl) < fabs(f)) {
        rts = x1;
        xl = x2;
        swap = fl;
        fl = f;
        f = swap;
    }
    else {
        xl = x1;
        rts = x2;
    }
    for (j = 1; j <= MAXIT; j++) {
        dx = (xl - rts) * f / (f - fl);
        xl = rts;
        fl = f;
        rts += dx;
        f = (*func)(rts);
        if (fabs(dx) < xacc || f == 0.0) {
            iteration = j - 1; return rts;
        }
    }
    return 0.0;
}
#undef MAXIT

#define JMAX 20
float rtnewt(void (*funcd)(float, float*, float*), float x1, float x2,
    float xacc)
{
    int j;
    float df, dx, f, rtn;

    rtn = 0.5 * (x1 + x2);
    for (j = 1; j <= JMAX; j++) {
        (*funcd)(rtn, &f, &df);
        dx = f / df;
        rtn -= dx;
        if (fabs(dx) < xacc) {
            iteration = j - 1; return rtn;
        }
    }
    return 0.0;
}
#undef JMAX 

#define MAXIT 100
float rtsafe(void (*funcd)(float, float*, float*), float x1, float x2,
    float xacc)
{
    int j;
    float df, dx, dxold, f, fh, fl;
    float temp, xh, xl, rts;

    (*funcd)(x1, &fl, &df);
    (*funcd)(x2, &fh, &df);

    if (fl == 0.0) return x1;
    if (fh == 0.0) return x2;
    if (fl < 0.0) {
        xl = x1;
        xh = x2;
    }
    else {
        xh = x1;
        xl = x2;
    }
    rts = 0.5 * (x1 + x2);
    dxold = fabs(x2 - x1);
    dx = dxold;
    (*funcd)(rts, &f, &df);
    for (j = 1; j <= MAXIT; j++) {
        if ((((rts - xh) * df - f) * ((rts - xl) * df - f) > 0.0)
            || (fabs(2.0 * f) > fabs(dxold * df))) {
            dxold = dx;
            dx = 0.5 * (xh - xl);
            rts = xl + dx;
            if (xl == rts) {
                iteration = j - 1; return rts;
            }
        }
        else {
            dxold = dx;
            dx = f / df;
            temp = rts;
            rts -= dx;
            if (temp == rts) {
                iteration = j - 1; return rts;
            }
        }
        if (fabs(dx) < xacc) {
            iteration = j - 1; return rts;
        }
        (*funcd)(rts, &f, &df);
        if (f < 0.0)
            xl = rts;
        else
            xh = rts;
    }

    return 0.0;
}
#undef MAXIT 

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//f(R)
float f1(float R) {
    return exp(-0.005 * R) * cos(sqrt(2000 - 0.01 * R * R) * 0.05) - 0.01;
}
void f1dx(float R, float *y, float *dy) {
    *y = f1(R);
    float root = sqrt(2000 - 0.01 * R * R) * 0.05;
    *dy = exp(-0.005 * R) * (sin(root) / root - cos(root)) * 0.005;
}

//8.32
float f2(float x) {
    return 100 * x / (8.85 * M_PI * pow(x * x + .9 * .9, 1.5)) - 1;
}
void f2dx(float x, float* y, float* dy) {
    *y = f2(x);
    *dy = -(4000 * x * x - 1620) / (177 * M_PI * pow(x * x + .9 * .9, 2.5));
}

//8.36
float f3(float R) {
    return 0.99403 + 1.671e-4 * R + 9.7215e-8 * pow(R, 2) - 9.5838e-11 * pow(R, 3) + 1.9520e-14 * pow(R, 4) -1.2;  //1.2 ÀÌÇ×
}
void f3dx(float R, float* y, float* dy) {
    *y = f3(R);
    *dy = 1.671e-4 + 9.7215e-8 * 2 * R - 9.5838e-11 * 3 * pow(R, 2) + 1.9520e-14 * 4 * pow(R, 3);
}

int main() {
    float a = 0.0;
    float b = 400;
    float xacc = 1e-4;
    int num_roots = 1;
    float root = 0;

    printf("==========R.E.=10^(-4)==========\n");

    printf("Bisection :");
    root = (*rtbis)(f1, 0, 400, xacc);
    printf("%.10f\n", root);
    printf("iteration:%d\n", iteration);
    iteration = 0;

    printf("Linear interpolation :");
    root = (*rtflsp)(f1, 0, 400, xacc);
    printf("%.10f\n", root);
    printf("iteration:%d\n", iteration);
    iteration = 0;

    printf("Secant :");
    root = (*rtsec)(f1, 0, 400, xacc);
    printf("%.10f\n", root);
    printf("iteration:%d\n", iteration); 
    iteration = 0;

    printf("Newton raphson :");
    root = (*rtnewt)(f1dx, 0, 400, xacc);
    printf("%.10f\n", root);
    printf("iteration:%d\n", iteration);
    iteration = 0;

    printf("Safe :");
    root = (*rtsafe)(f1dx, 0, 400, xacc);
    printf("%.10f\n", root);
    printf("iteration:%d\n", iteration);
    iteration = 0;

    printf("==========R.E.=10^(-6)==========\n");
    xacc = 1e-6;

    printf("Bisection :");
    root = (*rtbis)(f1, 0, 400, xacc);
    printf("%.10f\n", root);
    printf("iteration:%d\n", iteration);
    iteration = 0;

    printf("Linear interpolation :");
    root = (*rtflsp)(f1, 0, 400, xacc);
    printf("%.10f\n", root);
    printf("iteration:%d\n", iteration);
    iteration = 0;

    printf("Secant :");
    root = (*rtsec)(f1, 0, 400, xacc);
    printf("%.10f\n", root);
    printf("iteration:%d\n", iteration);
    iteration = 0;

    printf("Newton raphson :");
    root = (*rtnewt)(f1dx, 0, 400, xacc);
    printf("%.10f\n", root);
    printf("iteration:%d\n", iteration);
    iteration = 0;

    printf("Safe :");
    root = (*rtsafe)(f1dx, 0, 400, xacc);
    printf("%.10f\n", root);
    printf("iteration:%d\n", iteration);
    iteration = 0;

    printf("================================\n");
    
    printf("#8.32\n");
    root = (*rtsafe)(f2dx, 0, 1, xacc);
    printf("%.10f\n", root);
    root = (*rtsafe)(f2dx, 1, 400, xacc);
    printf("%.10f\n", root);

    printf("#8.36\n");
    root = (*rtsafe)(f3dx, 0, 10000, xacc);
    printf("%.10f\n", root);

    printf("================================\n");
    printf("Explain the concept of pointer to function and describe how you used it in your homework #3: \n");
    printf("\n  We can have pointers to functions like normal data pointers.(int)\n");
    printf("  Example: \n");
    printf("  void hello()\n  {\n      printf(\"hello\"); \n  }\n\n");
    printf("  int main()\n  {\n     void(*p)();\n      p=hello;\n      p();        ->invoking hello() \n  }");
    printf("\n\n  In HW#3:\n  void dxmethod (float (*method)(void (*func)(float, float*, float*), .....)\n  {...}");
    printf("\n  dxmethod(rtnewt, besseldx, bracket1,bracket2, num_roots, xacc);");
    printf("\n  -> method: rtnewt(...), func: bessedx(...)\n  ");
    printf("A corresponding function may be called according to the root finding method to be used.\n");

}
