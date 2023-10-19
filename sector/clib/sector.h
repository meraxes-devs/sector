#ifndef SECTOR_H
#define SECTOR_H

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<time.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                             *
 * SFHs related                                                                *
 *                                                                             *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
typedef long long llong_t;

typedef struct ssp_t {
    short index;
    float metals;
    float sfr;
} ssp_t;

typedef struct csp_t {
    ssp_t *bursts;
    int nBurst;
} csp_t;

typedef struct gal_params_t {
    double z;
    int nAgeStep;
    double *ageStep;
    int nGal;
    int *indices;
    csp_t *histories;
    llong_t *ids;
    int population;
} gal_params_t;


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                             *
 * SEDs and dust related                                                       *
 *                                                                             *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
typedef struct sed_params_t {
    // Raw templates
    int minZ;
    int maxZ;
    int nMaxZ;
    int nZ;
    double *Z;
    int nWaves;
    double *waves;
    int nAge;
    double *age;
    double *raw;
    // Redshift
    double z;
    // Filters
    int nFlux;
    int nObs;
    int *nFilterWaves;
    double *filterWaves;
    double *filters;
    double *centreWaves;
    double *logWaves;
    // IGM absorption
    int igm;
    // Working templates
    int nAgeStep;
    double *ageStep;
    double *integrated;
    double *ready;
    double *working;
    double *inBC;
    double *outBC;
} sed_params_t;

typedef struct dust_params_t {
    double tauUV_ISM;
    double nISM;
    double tauUV_BC;
    double nBC;
    double tBC;
} dust_params_t;

#ifndef _SECTOR_
void init_templates_raw(sed_params_t *spectra, char *fName);
void init_templates_rawIII(sed_params_t *spectra, char *fName);
void init_filters(sed_params_t *spectra,
                  double *betaBands, int nBeta, double *restBands, int nRest,
                  double *obsTrans, double *obsWaves, int *nObsWaves, int nObs, double z);
void shrink_templates_raw(sed_params_t *spectra, double maxAge);
void init_templates_integrated(sed_params_t *spectra);
void init_templates_interpolate(sed_params_t *spectra);
void init_templates_working(sed_params_t *spectra, csp_t *pHistories,
                            dust_params_t *dustParams, int iG);
void init_templates_workingIII(sed_params_t *spectra);
void fit_UV_slope(double *pTarget, double *pFit, int nGal, int nFlux,
                  double *logWaves, int nFit, int nR);
double *composite_spectra_cext(sed_params_t *spectra,
                               gal_params_t *galParams, dust_params_t *dustParams,
                               short outType, short approx, short nThread);
#endif

#ifndef _SPECTRA_
void init_templates_raw(sed_params_t *spectra, char *fName);
void shrink_templates_raw(sed_params_t *spectra, double maxAge);
void init_filters(sed_params_t *spectra,
                  double *betaBands, int nBeta, double *restBands, int nRest,
                  double *obsTrans, double *obsWaves, int *nObsWaves, int nObs, double z);
void init_templates_integrated(sed_params_t *spectra);
void init_templates_interpolate(sed_params_t *spectra);
void init_templates_working(sed_params_t *spectra, csp_t *pHistories,
                                   dust_params_t *dustParams, int iG);
void init_templates_workingIII(sed_params_t *spectra);
void free_templates_raw(sed_params_t *spectra);
void free_filters(sed_params_t *spectra);
#endif


#ifndef _DUST_
int birth_cloud_interval(double tBC, double *ageStep, int nAgeStep);
void init_templates_special(sed_params_t *spectra, double tBC, int approx);
void init_templates_specialIII(sed_params_t *spectra, int approx);
void dust_absorption(sed_params_t *spectra, dust_params_t *dustParams);
void dust_absorption_approx(double *inBCFlux, double *outBCFlux, double *centreWaves, int nFlux,
                            dust_params_t *dustParams);
#endif

#ifndef _IGM_
void add_Lyman_absorption(double *target, double *waves, int nWaves, double z);
void add_IGM_absorption_filters(sed_params_t *spectra);
void add_IGM_absorption_spectra(sed_params_t *spectra, double *pData, int nGal);
#endif

#endif
