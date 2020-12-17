/**************/
/* tools.h  */
/**************/

#define TINY 1.0e-50
#define PI 3.14159265359

/*definitions for ran1()*/
#define M1 259200
#define IA1 7141
#define IC1 54733
#define RM1 (1.0 / M1)
#define M2 134456
#define IA2 8121
#define IC2 28411
#define RM2 (1.0 / M2)
#define M3 243000
#define IA3 4561
#define IC3 51349

/********************/
/* matrix utilities */
/********************/
#include <iostream>

/*********************************************************************/
/* The following routines are from nrutil.c (Numerical Recipes in C )*/
/*********************************************************************/
void nerror(char *error_text)
{
    fprintf(stderr, "Run time error...\n");
    fprintf(stderr, "%s\n", error_text);
    fprintf(stderr, "Exiting to system...\n\n");
    exit(1);
}

/****************************/
int *ivector(int nl, int nh)
{
    int *v;

    v = (int *)malloc((unsigned)(nh - nl + 1) * sizeof(int));
    if (!v)
        nerror("allocation failure in ivector()");
    return v - nl;
}
/****************************/
void free_ivector(int *v, int nl)
{
    free((char *)(v + nl));
}
/****************************/
double *dvector(int nl, int nh)
{
    double *v;

    v = (double *)malloc((unsigned)(nh - nl + 1) * sizeof(double));
    if (!v)
        nerror("allocation failure in dvector()");
    return v - nl;
}
/******************************/
void free_dvector(double *v, int nl, int nh)
{
    free((char *)(v + nl));
}

/****************************/
double **dmatrix(int nrl, int nrh, int ncl, int nch)

{
    int i;
    double **m;

    m = (double **)malloc((unsigned)(nrh - nrl + 1) * sizeof(double *));
    if (!m)
        nerror("Allocation failure 1 in dmatrix()");
    m -= nrl;

    for (i = nrl; i <= nrh; i++)
    {
        m[i] = (double *)malloc((unsigned)(nch - ncl + 1) * sizeof(double));
        if (!m[i])
            nerror("Allocation failure 2 in dmatrix()");
        m[i] -= ncl;
    }
    return m;
}

/******************************/
void free_dmatrix(double **m, int nrl, int nrh, int ncl, int nch)
{
    int i;
    for (i = nrh; i >= nrl; i--)
        free((char *)(m[i] + ncl));
    free((char *)(m + nrl));
}

/****************************/
float **fmatrix(int nrl, int nrh, int ncl, int nch)
{
    int i;
    float **m;

    m = (float **)malloc((unsigned)(nrh - nrl + 1) * sizeof(float *));
    if (!m)
        nerror("Allocation failure 1 in fmatrix()");
    m -= nrl;

    for (i = nrl; i <= nrh; i++)
    {
        m[i] = (float *)malloc((unsigned)(nch - ncl + 1) * sizeof(float));
        if (!m[i])
            nerror("Allocation failure 2 in fmatrix()");
        m[i] -= ncl;
    }
    return m;
}

/******************************/
void free_fmatrix(float **m, int nrl, int nrh, int ncl, int nch)
{
    int i;
    for (i = nrh; i >= nrl; i--)
        free((char *)(m[i] + ncl));
    free((char *)(m + nrl));
}

/******************************/
int **imatrix(int nrl, int nrh, int ncl, int nch)
{
    int i;
    int **m;

    m = (int **)malloc((unsigned)(nrh - nrl + 1) * sizeof(int *));
    if (!m)
        nerror("Allocation failure 1 in imatrix()");
    m -= nrl;

    for (i = nrl; i <= nrh; i++)
    {
        m[i] = (int *)malloc((unsigned)(nch - ncl + 1) * sizeof(int));
        if (!m[i])
            nerror("Allocation failure 2 in imatrix()");
        m[i] -= ncl;
    }
    return m;
}

/******************************/
void free_imatrix(int **m, int nrl, int nrh, int ncl, int nch)
{
    int i;
    for (i = nrh; i >= nrl; i--)
        free((char *)(m[i] + ncl));
    free((char *)(m + nrl));
}

/* solvers */
/********************************/
/* LU Decomposition of a matrix */
/* From Numerical Recipes in C  */
/********************************/
void ludcmp(double **a, long n, int *indx, double *d)
{
    int i, imax, j, k;
    double big, dum, sum, temp;
    double *vv;

    vv = dvector(1, n);
    *d = 1.0;
    for (i = 1; i <= n; i++)
    {
        big = 0.0;
        for (j = 1; j <= n; j++)
            if ((temp = fabs(a[i - 1][j - 1])) > big)
                big = temp;
        if (big == 0.0)
            nerror("Singular matrix in routine LUDCMP");
        vv[i] = 1.0 / big;
    }
    for (j = 1; j <= n; j++)
    {
        for (i = 1; i < j; i++)
        {
            sum = a[i - 1][j - 1];
            for (k = 1; k < i; k++)
                sum -= a[i - 1][k - 1] * a[k - 1][j - 1];
            a[i - 1][j - 1] = sum;
        }
        big = 0.0;
        for (i = j; i <= n; i++)
        {
            sum = a[i - 1][j - 1];
            for (k = 1; k < j; k++)
                sum -= a[i - 1][k - 1] * a[k - 1][j - 1];
            a[i - 1][j - 1] = sum;
            if ((dum = vv[i] * fabs(sum)) >= big)
            {
                big = dum;
                imax = i;
            }
        }
        if (j != imax)
        {
            for (k = 1; k <= n; k++)
            {
                dum = a[imax - 1][k - 1];
                a[imax - 1][k - 1] = a[j - 1][k - 1];
                a[j - 1][k - 1] = dum;
            }
            *d = -(*d);
            vv[imax - 1] = vv[j - 1];
        }
        indx[j - 1] = imax;
        if (a[j - 1][j - 1] == 0.0)
            a[j - 1][j - 1] = TINY;
        if (j != n)
        {
            dum = 1.0 / (a[j - 1][j - 1]);
            for (i = j + 1; i <= n; i++)
                a[i - 1][j - 1] *= dum;
        }
    }
    free_dvector(vv, 1, n);
}

/**************************************/
void lubksb(double **a, long n, int *indx, double *b)
{
    int i, ii = 0, ip, j;
    double sum;

    for (i = 1; i <= n; i++)
    {
        ip = indx[i - 1];
        sum = b[ip - 1];
        b[ip - 1] = b[i - 1];
        if (ii)
            for (j = ii; j <= i - 1; j++)
                sum -= a[i - 1][j - 1] * b[j - 1];
        else if (sum)
            ii = i;
        b[i - 1] = sum;
    }
    for (i = n; i >= 1; i--)
    {
        sum = b[i - 1];
        for (j = i + 1; j <= n; j++)
            sum -= a[i - 1][j - 1] * b[j - 1];
        b[i - 1] = sum / (a[i - 1][i - 1]);
    }
}

/*****************************************/
void mprove(double **a, double **alud, int n, int *indx, double *b, double *x)
{
    int i, j;
    double sdp;
    double *r;
    r = dvector(1, n);
    for (i = 1; i <= n; i++)
    {
        sdp = -b[i];
        for (j = 1; j <= n; j++)
            sdp += a[i][j] * x[j];
        r[i] = sdp;
    }
    lubksb(alud, n, indx, r);
    for (i = 1; i <= n; i++)
        x[i] -= r[i];
    free_dvector(r, 1, n);
}

/*****************************/
void invert_matrix(double **g, double **h, int n, int *indx, double *col, double *d)
{
    int i, j;
    /* Invert "g" matrix return result in "h" matrix*/
    ludcmp(g, n, indx, d);
    for (j = 1; j <= n; j++)
    {
        for (i = 1; i <= n; i++)
            col[i - 1] = 0.0;
        col[j - 1] = 1.0;
        lubksb(g, n, indx, col);
        for (i = 1; i <= n; i++)
            h[i - 1][j - 1] = col[i - 1];
    }
}

/*****************************************/
void band_solver(double **a, double *c, int neq, int nbw)
{
    int i, j, k, l, n;
    int neqp, neqm, lim;
    double dum;
    /* Treat the case of one or more independent equations */
    if (nbw == 1)
    {
        for (n = 1; n <= neq; n++)
            c[n] /= a[n][1];
        return;
    }
    neqp = neq + 1;
    neqm = neq - 1;
    /* Forward reduction of the coefficient matrix [A] */
    for (n = 1; n <= neqm; n++)
    {
        if (nbw > (neqp - n))
            lim = (neqp - n);
        else
            lim = nbw;
        for (l = 2; l <= lim; l++)
        {
            dum = a[n][l] / a[n][1];
            if (dum != 0.0)
            { /* Modified to eliminate 0 mult.*/
                i = n + l - 1;
                j = 0;
                for (k = l; k <= lim; k++)
                {
                    j++;
                    a[i][j] -= (dum * a[n][k]);
                }
                a[n][l] = dum;
            }
        }
    }
    /* Forward reduction of the constant vector {c} */
    for (n = 1; n <= neqm; n++)
    {
        if (nbw > (neqp - n))
            lim = (neqp - n);
        else
            lim = nbw;
        for (l = 2; l <= lim; l++)
        {
            i = n + l - 1;
            c[i] -= (a[n][l] * c[n]);
        }
        c[n] /= a[n][1];
    }
    c[neq] /= a[neq][1];
    /* Back substitution. Former unknowns {x} overwrite {c} */
    for (n = neqm; n >= 1; n--)
    {
        if (nbw > (neqp - n))
            lim = (neqp - n);
        else
            lim = nbw;
        for (l = 2; l <= lim; l++)
        {
            k = n + l - 1;
            c[n] -= (a[n][l] * c[k]);
        }
    }
}

/*******************************/
void print_matrix(double **a, int m, int n)
{
    int i, j;
    for (i = 1; i <= m; i++)
    {
        for (j = 1; j <= n; j++)
        {
            printf("%12.4f", a[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

/*****************************/
void print_vector(double *v, int m)
{
    int i;
    for (i = 1; i <= m; i++)
    {
        printf("%12.4f", v[i]);
    }
    printf("\n");
}

void matrix_multiply(double **a, double **b, double **c, int l, int m, int n)
{
    int i, j, k;
    double temp;
    for (i = 1; i <= l; i++)
    {
        for (j = 1; j <= n; j++)
        {
            temp = 0.0;
            for (k = 1; k <= m; k++)
            {
                temp += a[i - 1][k - 1] * b[k - 1][j - 1];
            }
            c[i - 1][j - 1] = temp;
        }
    }
}

/*******************************/
void matrix_transpose(double **a, double **b, int m, int n)
{
    int i, j;
    for (i = 1; i <= m; i++)
    {
        for (j = 1; j <= n; j++)
        {
            a[i - 1][j - 1] = b[j - 1][i - 1];
        }
    }
}

/*******************************/
void matrix_add(double **a, double **b, int m, int n)
{
    int i, j;
    for (i = 1; i <= m; i++)
    {
        for (j = 1; j <= n; j++)
        {
            a[i][j] += b[i][j];
        }
    }
}

/*******************************/
void matrix_copy(double **a, double **b, int m, int n)
{
    int i, j;
    for (i = 1; i <= m; i++)
    {
        for (j = 1; j <= n; j++)
        {
            a[i][j] = b[i][j];
        }
    }
}

/*******************************/
void matrix_scalar_multiply(double **a, double **c, int m, int n)
{
    int i, j;
    for (i = 1; i <= m; i++)
        for (j = 1; j <= n; j++)
            a[i][j] = c[0][0]; //A Changer
}

/*******************************/
void matrix_vector_multiply(double **a, double *v, int k, int l, double *y)
{
    int i, j;
    double temp;
    for (i = 1; i <= k; i++)
    {
        temp = 0.0;
        for (j = 1; j <= l; j++)
        {
            temp += a[i][j] * v[j];
        }
        y[j] = temp;
    }
}

/*******************************/
void matrix_null(double **a, int m, int n)
{
    int i, j;
    for (i = 1; i <= m; i++)
    {
        for (j = 1; j <= n; j++)
        {
            a[i][j] = 0.0;
        }
    }
}

/*******************************/
void fmatrix_null(float **a, int k, int l, int m, int n)
{
    int i, j;
    for (i = k; i <= m; i++)
    {
        for (j = l; j <= n; j++)
        {
            a[i][j] = (float)0.0;
        }
    }
}

/*******************************/
void imatrix_null(int **a, int k, int l, int m, int n)
{
    int i, j;
    for (i = k; i <= m; i++)
    {
        for (j = l; j <= n; j++)
        {
            a[i][j] = 0;
        }
    }
}

/*******************************/
void dvector_copy(double *v, double *w, int n)
{
    int i;
    for (i = 1; i <= n; i++)
    {
        v[i] = w[i];
    }
}

/*******************************/
void vector_null(double *v, int n)
{
    int i;
    for (i = 1; i <= n; i++)
    {
        v[i] = 0.0;
    }
}

/*******************************/
void ivector_null(double *v, int n)
{
    int i;
    for (i = 1; i <= n; i++)
    {
        v[i] = 0;
    }
}

/*********************************/
/* convert double to nearest int */
/*********************************/
int d2i(double x)
{
    double intpart;
    double fraction;

    fraction = fabs(modf(x, &intpart));
    if (x >= 0.0)
    {
        if (fraction > 0.5)
            return ((int)(intpart + 1));
        return ((int)(intpart));
    }
    else
    {
        if (fraction > 0.5)
            return ((int)(intpart - 1));
        return ((int)(intpart));
    }
}
