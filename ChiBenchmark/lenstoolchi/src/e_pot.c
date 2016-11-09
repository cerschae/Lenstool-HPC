#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>

/****************************************************************/
/*      nom:        e_pot               */
/*      auteur:     Jean-Paul Kneib         */
/*      date:       10/02/92            */
/*      place:      Toulouse            */
/****************************************************************/

double e_pot(struct point pi, double dlsds)
{
    //const extern  struct  g_mode          M;
    const extern struct pot      lens[];
    const extern struct g_grille G;
    double R;
    int  i;
    double z; //u,
    double za, zb;//,zas,zbc;
    struct point P, Q;
    struct polar QP;
    double tpot;

    double  pis, pa, ps, t05, t15;

    /* ici l'ellipticite c'est celle du potentielle */

    tpot = 0.;

    for (i = 0; i < G.nlens; i++)
    {
        /*positionnement au centre du potentiel*/
        P.x = pi.x - lens[i].C.x;
        P.y = pi.y - lens[i].C.y;
        /*rotation dans les axes du potentiel*/
        Q = rotation(P, lens[i].theta);

        switch (lens[i].type)
        {
            case(1):
                tpot += lens[i].b0 * dlsds *
                        sqrt((1 + lens[i].epot) * Q.x * Q.x + (1 - lens[i].epot) * Q.y * Q.y);
                break;
            case(-1):
                QP = polxy(Q);
                z = sqrt(2 * lens[i].epot * 3.);
                za = z / sqrt(1 + lens[i].epot * 3.);
                zb = z / sqrt(1 - lens[i].epot * 3.);
                tpot += lens[i].b0 * dlsds / z *
                        (sin(QP.theta) * asin(za * sin(QP.theta))
                         + cos(QP.theta) * asinh(zb * cos(QP.theta)) );
                break;
            case(-2):
                pa = psiemd(Q.x, Q.y, lens[i].epot, lens[i].b0);
                ps = pi05(Q.x, Q.y, lens[i].epot, lens[i].rcut, lens[i].b0);
                tpot += dlsds * (pa - ps);
                break;
            case(3):
                z = lens[i].rc * lens[i].rc;
                tpot += lens[i].b0 * dlsds / lens[i].rc *
                        pow(1. + (1 + lens[i].epot) * Q.x * Q.x / z + (1 - lens[i].epot) * Q.y * Q.y / z,
                            lens[i].alpha);
                break;
            case(4):
                /* ici l'ellipticite c'est celle du potentielle */
                QP = polxy(Q);
                z = lens[i].rc * lens[i].rc;
                R = QP.r * QP.r / z;
                tpot += dlsds * lens[i].b0 * lens[i].rc * (sqrt(1. + R) +
                        lens[i].epot / 2.*R / sqrt(1. + R) * cos(2.*QP.theta));
                break;
            case(41):
                /* ici l'ellipticite c'est celle du potentielle */
                /* same as 4 but with a cut
                */
                QP = polxy(Q);
                z = lens[i].rc * lens[i].rc;
                R = QP.r * QP.r / z;
                if (sqrt(R) < lens[i].rcut)
                {
                    tpot += dlsds * lens[i].b0 * lens[i].rc * (sqrt(1. + R) +
                            lens[i].epot / 2.*R / sqrt(1. + R) * cos(2.*QP.theta));
                    tpot -= dlsds * lens[i].psicut * (Q.x * Q.x + Q.y * Q.y);
                }
                else
                {
                    tpot += dlsds *
                            (lens[i].psimcut * log(Q.x * Q.x + Q.y * Q.y) / 2. + lens[i].psiccut);
                }
                break;
            case(5):
                QP = polxy(Q);
                z = lens[i].rc * lens[i].rc;
                R = QP.r * QP.r / z;
                /* NON CALCULER POUR L'INSTANT >>>> 5=4 */
                tpot += dlsds * lens[i].b0 * lens[i].rc * (sqrt(1. + R) +
                        lens[i].epot * R / sqrt(1. + R) * cos(2.*QP.theta));
                break;
            case(6):
                QP = polxy(Q);
                z = lens[i].rc * lens[i].rc;
                R = QP.r * QP.r / z;
                tpot += dlsds * lens[i].b0 * lens[i].rc * (pow(1. + R, lens[i].alpha) +
                        lens[i].epot * R / pow(1. + R, lens[i].beta) * cos(2.*QP.theta));
                break;
            case(7):
                R = Q.x * Q.x + Q.y * Q.y;
                z = (G.dx * G.dx + G.dy * G.dy) / 4.;
                /* if(R>z) */
                if ((fabs(Q.x) > G.dx / 2.) || (fabs(Q.y) > G.dy / 2.))
                    tpot += dlsds * lens[i].b0 * log(R) / 2.;
                else
                    tpot += dlsds * lens[i].b0 * (log(z) / 2. + R / z / 2. - .5);
                break;
            case(71):
                R = Q.x * Q.x + Q.y * Q.y;
                z = (G.dx * G.dx + G.dy * G.dy) / 4.;
                if ((Q.x > G.dx / 2.) || (Q.y > G.dy / 2.))
                    tpot += dlsds * lens[i].b0 * log(R) / 2.;
                else
                    tpot += dlsds * lens[i].b0 * (log(z) / 2. + R / z / 2. - .5);
                break;
            case(8):
                tpot += dlsds * pi05(Q.x, Q.y, lens[i].epot, lens[i].rc, lens[i].b0);
                break;
            case(81): /* truncated PIEMD Kovner  NOT IMPLEMENTED */
                t05 = lens[i].rcut / (lens[i].rcut - lens[i].rc);
                pa = pi05(Q.x, Q.y, lens[i].epot, lens[i].rc, lens[i].b0);
                ps = pi05(Q.x, Q.y, lens[i].epot, lens[i].rcut, lens[i].b0);
                tpot += dlsds * t05 * (pa - ps);
                break;

            case(82): /* PIEMD Kovner with a shallower central slope  NOT IMPLEMENTED */
                t05 = (lens[i].rcut + lens[i].rc) / lens[i].rcut;
                pa = pi05(Q.x, Q.y, lens[i].epot, lens[i].rc, lens[i].b0);
                ps = pi05(Q.x, Q.y, lens[i].epot, lens[i].rcut, lens[i].b0);
                tpot += dlsds * t05 * (pa + ps);
                break;

            case(83): /* EMD kovner 3/2: I1.5c*/
                tpot += dlsds * pi15(Q.x, Q.y, lens[i].epot, lens[i].rc, lens[i].b0);
                break;

            case(84): /* EMD kovner isotherme I0.5c-I0.5cut + I1.5c */
                t05 = lens[i].rcut / (lens[i].rcut - lens[i].rc);
                t15 = lens[i].rcut / lens[i].rc;
                pa = pi05(Q.x, Q.y, lens[i].epot, lens[i].rc, lens[i].b0);
                ps = pi05(Q.x, Q.y, lens[i].epot, lens[i].rcut, lens[i].b0);
                pis = pi15(Q.x, Q.y, lens[i].epot, lens[i].rc, lens[i].b0);
                tpot += dlsds * (lens[i].alpha * t15 * PI + (1 - lens[i].alpha) * t05 * (pa - ps));
                break;

            case(9):
                tpot += dlsds * lens[i].b0 * (Q.x * Q.x + Q.y * Q.y) / 2.;
                break;
            case(10):
                tpot += dlsds * sp_pot(pi);
                break;
            case(11):
                z = lens[i].rc * lens[i].rc;
                tpot += dlsds * lens[i].b0 * lens[i].rc * log(1. + (Q.x * Q.x + Q.y * Q.y) / z);
                break;
            case(12):
                break;
            case(13):
                break;
            case(14):
                break;
            default:
                fprintf(stderr, "ERROR: profil not known (e_pot) %d\n", lens[i].type);
                exit(-1);
                break;
        };


    };

    return (tpot);
}
