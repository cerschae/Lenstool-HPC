#include <iostream>
#include <string.h>
//#include <cuda_runtime.h>
#include <math.h>
#include <sys/time.h>
#include <fstream>
#include <immintrin.h>

//#include "simd_math.h"

#define __WITH_LENSTOOL
//#include "structure.h"
#include "structure_hpc.h"
//#ifdef __WITH_LENSTOOL
//#endif
//#include "setup.hpp"

extern int num_lenses = 95;
//
extern double x_c = -44.094528948812552;
extern double y_c = -36.928385887965803;
//
//extern int lense_type[]  = {81,81,81,81,81,81,81,81,81,81,81,81,81,81,81,81,81,81,81,81,81,81,81,81,81,81,81,81,81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81};
//extern int lense_type[]  = {8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8};
extern int lense_type[]  = {5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5};
//
extern double pos_x[] = {0.000000000000000, 35.000000000000000, 65.000000000000000, 130.000000000000000, -19.399999999999999, 11.025931099982012, -2.671581400009952, 8.706480300003618, -46.071898399997828, -10.335049000010734, -25.638072299983587, -21.186195800000924, 47.770210000017265, 20.889859399994865, -9.947619900012000, 28.559020299985104, -11.757235600015839, 44.290179900007189, 49.190403899989448, 36.260348400000979, -47.958373599980469, 42.325706400014852, 37.455656899984845, 35.637785899984095, 48.099282799982589, 35.578290900001157, -3.050470599984273, 3.330580999983544, 43.051886999993762, 12.070367100002324, -0.012525300009348, 0.389991100006194, 18.864752299989611, 3.278487299995037, 10.076003800004774, 11.586436500008050, 37.425482399992134, 0.486777199988778, 25.632663699981624, 64.980204499991984, 18.813228000016544, 3.462950200001255, 28.896917700002792, -21.734176000018063, -0.928577399996435, -0.731873800008366, 78.370563400014419, 39.869047099991327, 10.636794000018757, 81.825827699993738, 11.692331899990036, 65.849002200017992, 12.142102700007115, 7.184091700006435, 92.001179599988561, 9.130061900012469, 1.379771499990539, 1.090836400019363, 63.261966300009220, -26.291947899982713, 64.747063800006757, 65.182885999997012, 69.098453600009847, 63.422517399987171, 83.228941699997549, 62.814757499999601, 58.726967500001948, 109.445168699990987, 65.110011700002232, 17.873548699988238, 108.654938499989214, 4.474934500009342, 3.973924000006243, 88.196061999997639, 61.916070000008830, 87.056832499986754, 108.613092799990099, 61.803057900019461, 79.079379299996546, 75.673077399983285, 129.107267799989302, 132.913239399982729, 124.503095399984701, 143.243164500019446, 128.260958600014533, 138.446558499998758, 130.452594800015163, 134.181991500016807, 131.942247000012628, 130.527746399989439, 132.704864599988241, 118.334119300008808, 147.068493300001194, 150.936237199997493, 142.130408800003835};

extern double pos_y[] = {0.000000000000000, -10.000000000000000, 40.000000000000000, 80.000000000000000, -21.699999999999999, -59.119920000003390, -58.477319999988708, -50.961599999999407, -45.297360000006393, -42.065999999999804, -46.933559999993690, -46.580759999989141, -37.085399999989477, -38.415960000011751, -38.951640000004772, -33.509519999998361, -35.354879999999866, -22.304880000001504, -22.674599999999145, -9.754919999988942, -22.953600000010965, -21.260159999994244, -20.612880000001610, -18.947520000006079, -7.068240000003811, -1.161000000010404, -1.597679999997581, -1.477079999995112, 5.219999999994229, -3.917879999997353, 0.001799999995455, -2.880720000001702, -0.453960000010056, 2.791800000011335, 1.789919999993117, 0.972720000009986, 2.251080000010575, 20.615760000009686, 20.495160000007218, 19.365120000000502, 16.301160000006121, 17.766359999990300, 18.265679999993267, 21.366720000011696, 32.603399999999283, 21.648600000006013, 23.630759999994666, 22.818239999989487, 21.551399999995624, 22.820040000010522, 25.769519999997215, 33.740639999987820, 25.422120000004611, 23.979960000008305, 33.112080000009314, 33.039720000007833, 32.423399999993308, 33.522119999992128, 33.943680000007248, 39.039480000008098, 33.088679999991655, 36.841680000006249, 40.110120000011307, 36.187200000006214, 43.166519999994080, 37.808640000000082, 41.374079999997093, 42.308279999991782, 39.722400000007951, 43.650719999999410, 45.423719999999435, 47.048760000009793, 48.891600000007429, 45.961559999992119, 53.340119999990065, 47.024999999987926, 46.107720000011909, 52.768080000006989, 58.089600000010932, 62.472239999991075, 77.337719999999877, 72.065520000009542, 74.942999999998960, 76.529879999998229, 77.375160000002552, 85.891679999997450, 83.807999999999083, 87.973920000004568, 87.793919999998593, 95.139359999993189, 95.914439999995693, 101.069640000000049, 101.447279999987927, 107.222760000001927, 106.741440000004673};

extern double theta[] = {1.047197551196598, 1.047197551196598, 1.047197551196598, 1.047197551196598, -0.698131700797732, -1.408480706359424, 0.359537825910832, 0.027925268031909, -0.991347015132779, -0.055850536063819, -0.991347015132779, 0.036651914291881, 0.687659725285766, 0.561996019142174, 0.907571211037051, 0.054105206811824, 0.202458193231342, -1.343903524035634, 0.420624349730633, 0.246091424531200, 0.150098315671512, -0.521853446346304, -0.040142572795870, 0.980875039620813, -0.884881930761125, -1.509709802975095, 0.666715774261834, -0.514872129338327, 1.174606586592184, -0.282743338823081, 1.457349925415265, 0.260054058547155, -1.164134611080218, -0.654498469497874, 1.038470904936626, -0.453785605518526, -0.886627260013119, -0.003490658503989, 1.549852375770965, 0.134390352403563, 0.293215314335047, -1.453859266911276, 0.109955742875643, 0.823795406941324, -0.349065850398866, 1.096066770252439, -0.527089434102287, -0.766199541625511, -1.186823891356144, 1.368338133563554, 0.748746249105567, 0.066322511575785, -0.651007810993885, -1.052433538952581, 0.256563400043166, 1.321214243759707, 1.110029404268394, 0.794124809657420, 0.134390352403563, -1.132718684544320, -0.961676417848876, -1.078613477732496, -1.530653753999027, -0.191986217719376, -0.127409035395586, -0.366519142918809, -0.815068760681352, -1.469567230179226, 0.026179938779915, 1.448623279155294, -0.146607657167524, 1.118756050528365, 0.837758040957278, -0.422369678982628, 1.024508270920671, -0.464257581030492, 0.226892802759263, 1.364847475059566, 0.249582083035189, -1.115265392024376, 1.418952681871390, -0.031415926535898, -1.127482696788337, 1.083849465488479, -0.047123889803847, -0.855211333477221, -0.867428638241182, 0.841248699461267, 0.289724655831059, 1.075122819228507, 0.099483767363677, -1.197295866868110, -0.534070751110265, 0.481710873550435, -0.174532925199433};
//
extern double  epot[] =  {0.075426688904937, 0.075426688904937, 0.075426688904937, 0.075426688904937, 0.116562483442831, 0.344322344322344,0.075409836065574, 0.034285714285715, 0.158045977011495, 0.064297800338409, 0.180914512922465, 0.085616438356165, 0.449781659388646, 0.198529411764706, 0.179611650485437, 0.160975609756097, 0.492772667542707, 0.216589861751152, 0.031496062992126, 0.043407043407043, 0.376311844077961, 0.157622739018088, 0.396428571428571, 0.373650107991361, 0.241081081081081, 0.364293085655315, 0.528824833702883, 0.317901234567901, 0.338345864661654, 0.129963898916967, 0.175792507204611, 0.151898734177215, 0.059829059829059, 0.208333333333333, 0.340163934426229, 0.590493601462523, 0.039062500000001, 0.015317286652079, 0.126005361930295, 0.123473541383989, 0.308584686774942, 0.200913242009132, 0.035587188612100, 0.096280087527352, 0.167197452229299, 0.109909909909910, 0.104109589041096, 0.051409618573797, 0.087431693989071, 0.176470588235294, 0.097178683385580, 0.045769764216366, 0.074235807860262, 0.131944444444444, 0.101449275362319, 0.285714285714286, 0.367619047619048, 0.239263803680982, 0.132132132132132, 0.318295739348371, 0.078817733990148, 0.108695652173913, 0.045143638850889, 0.381165919282511, 0.099236641221374,0.240437158469945, 0.334513274336283, 0.103594080338266, 0.186246418338109, 0.328918322295806, 0.098844672657253, 0.024911032028470, 0.102719033232628, 0.462857142857143, 0.039920159680638, 0.403565640194490, 0.023391812865496, 0.315315315315315, 0.274418604651163, 0.161987041036717, 0.137254901960784, 0.044673539518900, 0.057902973395931, 0.160283687943262, 0.122340425531915, 0.520395550061805, 0.021645021645023, 0.094555873925501, 0.142045454545455, 0.145907473309608, 0.194915254237288, 0.282904689863843, 0.450331125827815, 0.135245901639344, 0.340425531914894};
//
extern double rcore[] = {10.000000000000000, 5.000000000000000, 10.000000000000000, 10.000000000000000, 1.000000000000000, 0.034052339117993, 0.030349179184560, 0.046574365790645, 0.038383710030628, 0.068256961228370, 0.020804005013518, 0.028193367328818, 0.055481341027049, 0.039097316478764, 0.031633478431847, 0.039641214329569, 0.048322227956496, 0.041128886304460, 0.028585576136945, 0.080565099864725, 0.039097316478764, 0.027425019798032, 0.034367421532497, 0.025594517675956, 0.049905332073874, 0.073815395284787, 0.052984269343359, 0.048322227956496,0.067320440611066, 0.044683475951999, 0.067320440611066, 0.022811137889776, 0.036320092822069, 0.039277781492078, 0.037857065785625, 0.021288592541707, 0.024668738981743, 0.069206510145232, 0.053474526751449, 0.048322227956496, 0.025359865727459, 0.029117021675657, 0.031199450836488, 0.051067776588639, 0.099116726380952, 0.029386438160114, 0.050833141468429, 0.034367421532497, 0.021784467563449, 0.027934888983489, 0.044273815382088, 0.052498506642202, 0.035821762466628, 0.013432193817836, 0.052257297914513, 0.037337647370553, 0.023342477548558, 0.043266020305102, 0.039641214329569, 0.061965070730691, 0.021093417552129, 0.028717521161380, 0.051067776588639, 0.032669835654398, 0.052498506642202, 0.017707080013958, 0.037509990120305, 0.030349179184560, 0.013936282772780, 0.028454237344829, 0.057829168758607, 0.036825355643309, 0.040192678553178, 0.029522079795690, 0.056512816167846, 0.034685419368686, 0.023886193697910, 0.035821762466628, 0.038738870126742, 0.034685419368686, 0.083588580432868, 0.032221588704875, 0.042672388251826, 0.052017197440259, 0.033124318340140, 0.048100207486392, 0.026677611870563, 0.023996447358569, 0.016000993535772, 0.057563467863253, 0.034526054342700, 0.045935341535597, 0.039459079494046, 0.031488135800734, 0.039277781492078};
//
extern double rcut[] = {500.000000000000000, 500.000000000000000, 500.000000000000000, 500.000000000000000, 150.000000000000000, 5.107850867699009, 4.552376877684008, 6.986154868596703, 5.757556504594200, 10.238544184255536, 3.120600752027715, 4.229005099322641, 8.322201154057367, 5.864597471814614, 4.745021764777110, 5.946182149435385, 7.248334193474441, 6.169332945669035, 4.287836420541705, 12.084764979708812, 5.864597471814614, 4.113752969704730, 5.155113229874521, 3.839177651393335, 7.485799811081068, 11.072309292717980, 7.947640401503818, 7.248334193474441, 10.098066091659968, 6.702521392799788, 10.098066091659968, 3.421670683466420, 5.448013923310371, 5.891667223811720, 5.678559867843796, 3.193288881256013, 3.700310847261404, 10.380976521784758, 8.021179012717386, 7.248334193474441, 3.803979859118900, 4.367553251348570, 4.679917625473273, 7.660166488295854, 14.867508957142872, 4.407965724017149, 7.624971220264309, 5.155113229874521, 3.267670134517328, 4.190233347523285, 6.641072307313197, 7.874775996330366, 5.373264369994223, 2.014829072675421, 7.838594687177010, 5.600647105583010, 3.501371632283681, 6.489903045765340, 5.946182149435385, 9.294760609603639, 3.164012632819346, 4.307628174206992, 7.660166488295854, 4.900475348159637, 7.874775996330366, 2.656062002093733, 5.626498518045686, 4.552376877684008, 2.090442415917007, 4.268135601724389, 8.674375313791064, 5.523803346496345, 6.028901782976731, 4.428311969353570, 8.476922425176864, 5.202812905302965, 3.582929054686518, 5.373264369994223, 5.810830519011283, 5.202812905302965, 12.538287064930271, 4.833238305731319, 6.400858237773873, 7.802579616038900, 4.968647751020976, 7.215031122958762, 4.001641780584524, 3.599467103785322, 2.400149030365801, 8.634520179487971, 5.178908151405005, 6.890301230339541, 5.918861924106924, 4.723220370110043, 5.891667223811720};
//
double b0[] = {27.686284800000006, 15.573535200000000, 21.197311800000001, 43.259820000000005, 1.401618168000000, 0.701887044377384, 0.625557486765266, 0.959991142907221, 0.791165290944736, 1.406912946824478, 0.428812292146906, 0.581121878203467, 1.143581776765294, 0.805874151883638, 0.652029471542907, 0.817084978065710, 0.996018089699853, 0.847748882881212, 0.589206088811179, 1.660608383702515, 0.805874151883638, 0.565284693698478, 0.708381525237709, 0.527554371568815, 1.028649042482203, 1.521483427216619, 1.092112118320076, 0.996018089699853, 1.387609377521935, 0.921016108754891, 1.387609377521935, 0.470183328577265, 0.748629999074688, 0.809593897959512, 0.780310096202241, 0.438800613557777, 0.508472214857728, 1.426485055525682, 1.102217307333858, 0.996018089699853, 0.522717724018210, 0.600160247852116, 0.643083291809730, 1.052609356688346, 2.042994635018428, 0.605713462253246, 1.047773056002485, 0.708381525237709, 0.449021592862168, 0.575794120800079, 0.912572182315177, 1.082099574236042, 0.738358410415477, 0.276864469886694, 1.077127778308336, 0.769603840315186, 0.481135305220354, 0.891799503367263, 0.817084978065710, 1.277224457300750, 0.434777665351592, 0.591925740547885, 1.052609356688346, 0.673390873628514, 1.082099574236042, 0.364978453094982, 0.773156170239529, 0.625557486765266, 0.287254754837836, 0.586498932739152, 1.191975217859029, 0.759044479765867, 0.828451762375228, 0.608509308563808, 1.164842549348754, 0.714936098789562, 0.492342386170850, 0.738358410415477, 0.798485853249674, 0.714936098789562, 1.722928385636979, 0.664151604471083, 0.879563555467234, 1.072178825707933, 0.682758673823120, 0.991441802267353, 0.549879116438807, 0.494614935370580, 0.329812587059260, 1.186498589205643, 0.711651265795296, 0.946819570637053, 0.813330813603161, 0.649033668246664, 0.80959389795951};
//
extern double   grad_x = -77.124413458712169;
extern double   grad_y = -46.092652444451353; 

//
//
//

/**** usefull functions for PIEMD profile : see old lenstool ****/

/**  I*w,v=0.5 Kassiola & Kovner, 1993 PIEMD, paragraph 4.1
 *
 * Global variables used :
 * - none
 */

/** @brief This module function calculates profile depended information like the impactparameter b0 and the potential ellipticity epot
 * 
 * @param lens: mass distribution for which to calculate parameters
 */
void module_readParameters_calculatePotentialparameter(Potential *lens)
{
        //std::cout << "module_readParameters_calculatePotentialparameter..." << std::endl;
        switch (lens->type)
        {

                case(5): /*Elliptical Isothermal Sphere*/
                        //impact parameter b0
                        lens->b0 = 4* pi_c2 * lens->vdisp * lens->vdisp ;
                        //ellipticity_potential 
                        lens->ellipticity_potential = lens->ellipticity/3 ;
                        break;

                case(8): /* PIEMD */
                        //impact parameter b0
                        lens->b0 = 6.*pi_c2 * lens->vdisp * lens->vdisp;
                        //ellipticity_parameter
                        if ( lens->ellipticity == 0. && lens->ellipticity_potential != 0. ){
                                // emass is (a2-b2)/(a2+b2)
                                lens->ellipticity = 2.*lens->ellipticity_potential / (1. + lens->ellipticity_potential * lens->ellipticity_potential);
                                //printf("1 : %f %f \n",lens->ellipticity,lens->ellipticity_potential);
                        }
                        else if ( lens->ellipticity == 0. && lens->ellipticity_potential == 0. ){
                                lens->ellipticity_potential = 0.00001;
                                //printf("2 : %f %f \n",lens->ellipticity,lens->ellipticity_potential);
                        }
                        else{
                                // epot is (a-b)/(a+b)
                                lens->ellipticity_potential = (1. - sqrt(1 - lens->ellipticity * lens->ellipticity)) / lens->ellipticity;
                                //printf("3 : %f %f \n",lens->ellipticity,lens->ellipticity_potential);
                        }
                        break;

                default:
                        std::cout << "ERROR: LENSPARA profil type of clump "<< lens->name << " unknown : "<< lens->type << std::endl;
                        //printf( "ERROR: LENSPARA profil type of clump %s unknown : %d\n",lens->name, lens->type);
                        break;
        };
}
//
//
//
void module_readParameters_calculatePotentialparameter_SOA(Potential_SOA *lens, int ii)
{
        //std::cout << "module_readParameters_calculatePotentialparameter..." << std::endl;
        switch (lens->type[ii])
        {

                case(5): /*Elliptical Isothermal Sphere*/
                        //impact parameter b0
                        lens->b0[ii] = 4* pi_c2 * lens->vdisp[ii] * lens->vdisp[ii] ;
                        //ellipticity_potential 
                        lens->ellipticity_potential[ii] = lens->ellipticity[ii]/3 ;
                        break;

                case(8): /* PIEMD */
                        //impact parameter b0
                        lens->b0[ii] = 6.*pi_c2 * lens->vdisp[ii] * lens->vdisp[ii];
                        //ellipticity_parameter
                        if ( lens->ellipticity[ii] == 0. && lens->ellipticity_potential[ii] != 0. ){
                                // emass is (a2-b2)/(a2+b2)
                                lens->ellipticity[ii] = 2.*lens->ellipticity_potential[ii] / (1. + lens->ellipticity_potential[ii] * lens->ellipticity_potential[ii]);
                                //printf("1 : %f %f \n",lens->ellipticity[ii],lens->ellipticity_potential[ii]);
                        }
                        else if ( lens->ellipticity[ii] == 0. && lens->ellipticity_potential[ii] == 0. ){
                                lens->ellipticity_potential[ii] = 0.00001;
                                //printf("2 : %f %f \n",lens->ellipticity[ii],lens->ellipticity_potential[ii]);
                        }
                        else{
                                // epot is (a-b)/(a+b)
                                lens->ellipticity_potential[ii] = (1. - sqrt(1 - lens->ellipticity[ii] * lens->ellipticity[ii])) / lens->ellipticity[ii];
                                //printf("3 : %f %f \n",lens->ellipticity[ii],lens->ellipticity_potential[ii]);
                        }
                        break;

                default:
                        std::cout << "ERROR: LENSPARA profil type of clump "<< lens->name << " unknown : "<< lens->type << std::endl;
                        //printf( "ERROR: LENSPARA profil type of clump %s unknown : %d\n",lens->name, lens->type);
                        break;
        };

}
//
//
//
void 
setup_jauzac_SOA(Potential_SOA *lens_soa, int* nlenses, double* x, double* y, double* sol_grad_x, double* sol_grad_y)
{
	*nlenses    = num_lenses;
	*sol_grad_x = grad_x; 
	*sol_grad_y = grad_y; 
	*x	    = x_c;
	*y	    = y_c;
	//
        //Potential_SOA lens_soa
        //
        for (int i = 0; i < *nlenses; ++i)
        {
                lens_soa->type       = (int*) malloc(*nlenses*sizeof(int));
                lens_soa->vdisp      = (double*) malloc(*nlenses*sizeof(double));
                //
                lens_soa->position_x = (double*) malloc(*nlenses*sizeof(double));
                lens_soa->position_y = (double*) malloc(*nlenses*sizeof(double));
                lens_soa->weight = (double*) malloc(*nlenses*sizeof(double));
                lens_soa->b0 = (double*) malloc(*nlenses*sizeof(double));
                lens_soa->vdisp = (double*) malloc(*nlenses*sizeof(double));
                lens_soa->ellipticity_angle = (double*) malloc(*nlenses*sizeof(double));
                lens_soa->anglecos = (double*) malloc(*nlenses*sizeof(double));
                lens_soa->anglesin = (double*) malloc(*nlenses*sizeof(double));
                lens_soa->ellipticity = (double*) malloc(*nlenses*sizeof(double));
                lens_soa->ellipticity_potential = (double*) malloc(*nlenses*sizeof(double));
                lens_soa->rcore = (double*) malloc(*nlenses*sizeof(double));
                lens_soa->rcut = (double*) malloc(*nlenses*sizeof(double));
                lens_soa->rscale = (double*) malloc(*nlenses*sizeof(double));
                lens_soa->exponent = (double*) malloc(*nlenses*sizeof(double));
                lens_soa->alpha = (double*) malloc(*nlenses*sizeof(double));
                lens_soa->einasto_kappacritic = (double*) malloc(*nlenses*sizeof(double));
                lens_soa->z = (double*) malloc(*nlenses*sizeof(double));
                lens_soa->mag = (double*) malloc(*nlenses*sizeof(double));
                lens_soa->lum = (double*) malloc(*nlenses*sizeof(double));
                lens_soa->theta = (double*) malloc(*nlenses*sizeof(double));
                lens_soa->sigma = (double*) malloc(*nlenses*sizeof(double));
        }
        //
        for (int i = 0; i < *nlenses; ++i)
        {
                //ilens = &lens[i];
                //
                lens_soa->vdisp[i]                 = 1.;
                lens_soa->position_x[i]            = pos_x[i];
                lens_soa->position_y[i]            = pos_y[i];
                lens_soa->type[i]                  = lense_type[i];
                lens_soa->ellipticity[i]           = 0.11;
                lens_soa->ellipticity_potential[i] = epot[i];
                lens_soa->ellipticity_angle[i]     = theta[i];
		lens_soa->anglecos[i]		   = cos(theta[i]);
		lens_soa->anglesin[i]		   = sin(theta[i]);
                lens_soa->rcut[i]                  = rcut[i];
                lens_soa->rcore[i]                 = rcore[i];
                lens_soa->b0[i]                    = b0[i];
                lens_soa->weight[i]                = 0;
                lens_soa->rscale[i]                = 0;
                lens_soa->exponent[i]              = 0;
                lens_soa->alpha[i]                 = 0.;
                lens_soa->einasto_kappacritic[i]   = 0;
                lens_soa->z[i]                     = 0.4;
        }
}

void
setup_jauzac(Potential** lens, int* nlenses, double* x, double* y, double* sol_grad_x, double* sol_grad_y)
{
        *nlenses    = num_lenses;
        *sol_grad_x = grad_x;
        *sol_grad_y = grad_y;
        *x          = x_c;
        *y          = y_c;
        //
       	*lens = (Potential*) malloc(sizeof(Potential)*(*nlenses)); 
        //
        for (int i = 0; i < *nlenses; ++i)
        {
                //ilens = &lens[i];
                //
		(*lens)[i].position.x            = pos_x[i];
                (*lens)[i].position.y            = pos_y[i];
		//
                (*lens)[i].vdisp                 = 1.;
                (*lens)[i].type                  = lense_type[i];
                (*lens)[i].ellipticity           = 0.11;
                (*lens)[i].ellipticity_potential = epot[i];
                (*lens)[i].ellipticity_angle     = theta[i];
                (*lens)[i].rcut                  = rcut[i];
                (*lens)[i].rcore                 = rcore[i];
                (*lens)[i].b0                    = b0[i];
                (*lens)[i].weight                = 0;
                (*lens)[i].rscale                = 0;
                (*lens)[i].exponent              = 0;
                (*lens)[i].alpha                 = 0.;
                (*lens)[i].einasto_kappacritic   = 0;
                (*lens)[i].z                     = 0.4;
        }
	std::cout << "Setup done..." << std::endl;
}

#ifdef __WITH_LENSTOOL
extern struct pot lens[NLMAX];
void
//setup_jauzac_LT(struct pot** lens, int* nlenses, double* x, double* y, double* sol_grad_x, double* sol_grad_y)
setup_jauzac_LT( int* nlenses, double* x, double* y, double* sol_grad_x, double* sol_grad_y)
{
        *nlenses    = num_lenses;
        *sol_grad_x = grad_x;
        *sol_grad_y = grad_y;
        *x          = x_c;
        *y          = y_c;
        //
        //*lens = (struct pot*) malloc(sizeof(struct pot)*(*nlenses));
        //
        for (int i = 0; i < *nlenses; ++i)
        {
                //
                lens[i].C.x            = pos_x[i];
                lens[i].C.y            = pos_y[i];
                //
                lens[i].sigma                 = 1.; 			// disp
		if (lense_type[i] == 5)
			lens[i].type                  = 1;
		else
			lens[i].type                  = lense_type[i];
                lens[i].emass                 = 0.11;
                lens[i].epot                  = epot[i];
                lens[i].theta                 = theta[i];
                lens[i].rcut                  = rcut[i];
                lens[i].rc                    = rcore[i];
                lens[i].b0                    = b0[i];
                lens[i].masse                 = 0;			// weight
                //lens[i].rc  	                 = 0;			// rscale
//               (&lens)[i].exponent              = 0;
                lens[i].alpha                 = 0.;
//                (&lens)[i].einasto_kappacritic   = 0;
                lens[i].z                     = 0.4;
        }
        std::cout << "Setup done..." << std::endl;
}
#endif








