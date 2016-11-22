void bicubic_coef(double *y, double *y1, double *y2, double *y12, double d1, double d2, double **c)
{
    const int wt[16][16] =
    {{1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0},
        { -3, 0, 0, 3, 0, 0, 0, 0, -2, 0, 0, -1, 0, 0, 0, 0},
        {2, 0, 0, -2, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0},
        {0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0},
        {0, 0, 0, 0, -3, 0, 0, 3, 0, 0, 0, 0, -2, 0, 0, -1},
        {0, 0, 0, 0, 2, 0, 0, -2, 0, 0, 0, 0, 1, 0, 0, 1},
        { -3, 3, 0, 0, -2, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, -3, 3, 0, 0, -2, -1, 0, 0},
        {9, -9, 9, -9, 6, 3, -3, -6, 6, -6, -3, 3, 4, 2, 1, 2},
        { -6, 6, -6, 6, -4, -2, 2, 4, -3, 3, 3, -3, -2, -1, -1, -2},
        {2, -2, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 2, -2, 0, 0, 1, 1, 0, 0},
        { -6, 6, -6, 6, -3, -3, 3, 3, -4, 4, 2, -2, -2, -2, -1, -1},
        {4, -4, 4, -4, 2, 2, -2, -2, 2, -2, -2, 2, 1, 1, 1, 1}
    };
    int l, k, j, i;
    double xx, d1d2, cl[16], x[16];

    d1d2 = d1 * d2;
    for (i = 0; i < 4; i++)
    {
        x[i] = y[i];
        x[i+4] = y1[i] * d1;
        x[i+8] = y2[i] * d2;
        x[i+12] = y12[i] * d1d2;
    }
    for (i = 0; i < 16; i++)
    {
        xx = 0.0;
        for (k = 0; k <= 15; k++) xx += wt[i][k] * x[k];
        cl[i] = xx;
    }
    l = 0;
    for (i = 0; i < 4; i++)
        for (j = 0; j < 4; j++) c[i][j] = cl[l++];
}