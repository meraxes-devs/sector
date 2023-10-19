#define _SECTOR_
#include"tools.h"
#include"sector.h"

//#define SURFACE_AREA 1.1965e40 // 4*pi*(10 pc)**2 unit cm^2
//#define JANSKY(x) (3.34e4*(x)*(x))
#define M_AB(x) (-2.5*log10(x) + 8.9) // Convert Jansky to AB magnitude
#define TOL 1e-30 // Minimum Flux


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                             *
 * Profiling functions                                                         *
 *                                                                             *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#ifdef TIMING
    #define MAX_BLOCK 100
    #define INTEGRATION 0
    #define DUST 1
    #define WORKING1 2
    #define WORKING2 3
    #define SUM 4
    #define FIT 5

    static double timer[MAX_BLOCK];
    static double counter[MAX_BLOCK];
    static char blockNames[MAX_BLOCK][128];
    static struct timespec g_sTime;
    static struct timespec g_eTime;
    static struct timespec g_sTime2;
    static struct timespec g_eTime2;


    void timing_start(char* text) {
        clock_gettime(CLOCK_REALTIME, &g_sTime);
        printf("#***********************************************************\n");
        printf("# %s\n", text);
    }


    void timing_end(void) {
        clock_gettime(CLOCK_REALTIME, &g_eTime);
        double elapsedTime = g_eTime.tv_sec - g_sTime.tv_sec \
                             + (g_eTime.tv_nsec - g_sTime.tv_nsec)/1e9;
        int minute = (int)elapsedTime/60;
        printf("# 100.0%% complete!\n");
        printf("# Done!\n");
        printf("# Elapsed time: %d min %.6f sec\n", minute, elapsedTime - minute*60);
        printf("#***********************************************************\n\n");
    }


    void init_profiler(void) {
        int iB;
        for(iB = 0; iB < MAX_BLOCK; ++iB) {
            timer[iB] = 0.;
            counter[iB] = 0;
        }
    }


    void profiler_start(char* name, int blockIdx) {
        clock_gettime(CLOCK_REALTIME, &g_sTime2);
        strcpy(blockNames[blockIdx], name);
    }


    void profiler_end(int blockIdx) {
        clock_gettime(CLOCK_REALTIME, &g_eTime2);
        timer[blockIdx] += g_eTime2.tv_sec - g_sTime2.tv_sec \
                           + (g_eTime2.tv_nsec - g_sTime2.tv_nsec)/1e9;
        counter[blockIdx] += 1;
    }


    void profiler_summary(void) {
        int iB, ncall;
        printf("#***********************************************************\n");
        for(iB = 0; iB < MAX_BLOCK; ++iB) {
            ncall = counter[iB];
            if (ncall == 0)
                continue;
            printf("# %s\n", blockNames[iB]);
            printf("#  call: %6d  total: %2.3f sec  mean: %2.3f ms\n",
                   ncall, timer[iB], timer[iB]/ncall*1e3);
        }
        printf("#***********************************************************\n\n");
    }


    inline void report(int i, int tot) {
        int n = tot > 10 ? tot/10 : 1;
        if (i%n == 0) {
            printf("# %5.1f%% complete!\r", 100.*(i + 1)/tot);
            fflush(stdout);
        }
    }
#endif


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                             *
 * Galaxy properites related                                                   *
 *                                                                             *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
void trim_gal_params(gal_params_t *galParams, int minZ, int maxZ) {
    /* Set the metallicity of each SSP to the given range */
    int iB, iG;

    float f_minZ = (minZ + 1.)/1000.;
    float f_maxZ = (maxZ + 1.)/1000.;
    int metals;
    int nGal = galParams->nGal;
    csp_t *pHistories = galParams->histories;
    int nBurst;
    ssp_t *pBursts;

    for(iG = 0; iG < nGal; ++iG) {
        nBurst = pHistories->nBurst;
        pBursts = pHistories->bursts;
        for(iB = 0; iB < nBurst; ++iB) {
            metals = (int)(pBursts->metals*1000 - .5);
            if (metals < minZ)
                pBursts->metals = f_minZ;
            else if (metals > maxZ)
                pBursts->metals = f_maxZ;
            ++pBursts;
        }
        ++pHistories;
    }
}


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                             *
 * Primary Functions                                                           *
 *                                                                             *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
void fit_UV_slope(double *pTarget, double *pFit, int nGal, int nFlux,
                  double *logWaves, int nFit, int nR) {
    #ifdef TIMING
        profiler_start("Slope fit", FIT);
    #endif
    int iF, iG;
    struct linResult result;
    double *logf = malloc(nFit*sizeof(double));

    for(iG = 0; iG < nGal; ++iG) {
        for(iF = 0; iF < nFit; ++iF)
            logf[iF] = log(pFit[iF]);
        pFit += nFlux;
        result = linregress(logWaves, logf, nFit);
        *pTarget++ = (double)result.slope;
        if (nR > 1)
            *pTarget++ = (double)result.intercept;
        if (nR > 2)
            *pTarget++ = (double)result.R;
    }
    #ifdef TIMING
        profiler_end(FIT);
    #endif
}


void compute_spectra(double *target, sed_params_t *spectra,
                     gal_params_t *galParams, dust_params_t *dustParams,
                     int approx, short nThread, int population) {
    // NOTE: if dealing with multiple galaxies, they have to be the same population
    // Initialise SED templates
    spectra->ready = NULL;
    spectra->working = NULL;

    // Integrate SED templates over given time steps
    spectra->nAgeStep = galParams->nAgeStep;
    spectra->ageStep = galParams->ageStep;
    if (population == 2)
        init_templates_integrated(spectra);
    else
        init_templates_interpolate(spectra);

    // Initialise templates for birth cloud if necessary
    if (dustParams == NULL) {
        spectra->inBC = NULL;
        spectra->outBC = NULL;
    }

    #ifdef TIMING
        profiler_start("Summation over progenitors", SUM);
    #endif
    #pragma omp parallel \
    default(none) \
    firstprivate(spectra, galParams, dustParams, target, approx, population) \
    num_threads(nThread)
    {
        int iF, iG;
        int offset;
        double *pTarget;

        int iP, nProg;
        int nAgeStep = galParams->nAgeStep;
        int nGal = galParams->nGal;
        csp_t *histories = galParams->histories;
        csp_t *pHistories;
        ssp_t *pBursts;
        double sfr;
        int metals;

        int nFlux = spectra->nFlux;
        sed_params_t omp_spectra;
        memcpy(&omp_spectra, spectra, sizeof(sed_params_t));
        double *readyData = malloc(spectra->nZ*nAgeStep*spectra->nWaves*sizeof(double));
        double *workingData = malloc(spectra->nMaxZ*nAgeStep*nFlux*sizeof(double));
        omp_spectra.ready = readyData;
        omp_spectra.working = workingData;


        if (approx) {
            if (population == 2){
               init_templates_working(&omp_spectra, histories, NULL, -1);
               init_templates_special(&omp_spectra, dustParams->tBC, approx);
            }
            else{
               init_templates_workingIII(&omp_spectra);
               init_templates_specialIII(&omp_spectra, approx);
            }

            int iAge, iAgeBC;
            // Assume that tBC is the same for every galaxy
            if (population == 2)
                iAgeBC = birth_cloud_interval(dustParams->tBC, galParams->ageStep, nAgeStep);
            else
                iAgeBC = 0;

            double *inBC = omp_spectra.inBC;
            double *outBC = omp_spectra.outBC;
            double *inBCFlux = malloc(nFlux*sizeof(double));
            double *outBCFlux = malloc(nFlux*sizeof(double));

            #pragma omp for schedule(static, 1)
            for(iG = 0; iG < nGal; ++iG) {
                pHistories = histories + iG;

                // Sum contributions from all progenitors
                nProg = pHistories->nBurst;
                pTarget = target + iG*nFlux;
                for(iF = 0; iF < nFlux; ++iF) {
                    inBCFlux[iF] = 0.;
                    outBCFlux[iF] = 0.;
                }
                for(iP = 0; iP < nProg; ++iP) {
                    pBursts = pHistories->bursts + iP;
                    iAge = pBursts->index;
                    sfr = pBursts->sfr;
                    if (population == 2)
                        metals = (int)(pBursts->metals*1000 - .5);
                    else
                        metals = 0;
                    if (iAge > iAgeBC) {
                        offset = (metals*nAgeStep + iAge)*nFlux;
                        for(iF = 0; iF < nFlux; ++iF)
                            outBCFlux[iF] += sfr*workingData[offset + iF];
                    }
                    else if (iAge == iAgeBC) {
                        offset = metals*nFlux;
                        for(iF = 0; iF < nFlux; ++iF) {
                            inBCFlux[iF] += sfr*inBC[offset + iF];
                            outBCFlux[iF] += sfr*outBC[offset + iF];
                        }
                    }
                    else {
                        offset = (metals*nAgeStep + iAge)*nFlux;
                        for(iF = 0; iF < nFlux; ++iF)
                            inBCFlux[iF] += sfr*workingData[offset +iF];
                    }
                }

                // Apply dust absorption
                if (population == 2)
                    dust_absorption_approx(inBCFlux, outBCFlux, omp_spectra.centreWaves, nFlux,
                                       dustParams + iG);

                for(iF = 0; iF < nFlux; ++iF)
                    pTarget[iF] += inBCFlux[iF] + outBCFlux[iF];
                #ifdef TIMING
                    report(iG, nGal);
                #endif
            }
            free(inBC);
            free(outBC);
            free(inBCFlux);
            free(outBCFlux);
        }
        else {
            if (population == 2){
                if (dustParams == NULL)
                    init_templates_working(&omp_spectra, histories, NULL, -1);
                else
                    //  -Assume that tBC is the same for every galaxy
                    init_templates_special(&omp_spectra, dustParams->tBC, approx);
            }
            else{
                if (dustParams == NULL)
                    init_templates_workingIII(&omp_spectra);
                else
                    init_templates_specialIII(&omp_spectra, approx);
            }

            #pragma omp for schedule(static, 1)
            for(iG = 0; iG < nGal; ++iG) {
                pHistories = histories + iG;
                // Add dust absorption to SED templates
                if (population == 2)
                    init_templates_working(&omp_spectra, pHistories, dustParams, iG);
                else
                    init_templates_workingIII(&omp_spectra);
                // Sum contributions from all progenitors
                nProg = pHistories->nBurst;
                pTarget = target + iG*nFlux;
                for(iP = 0; iP < nProg; ++iP) {
                    pBursts = pHistories->bursts + iP;
                    sfr = pBursts->sfr;
                    if (population == 2)
                        metals = (int)(pBursts->metals*1000 - .5);
                    else
                        metals = 0;
                    offset = (metals*nAgeStep + pBursts->index)*nFlux;
                    for(iF = 0 ; iF < nFlux; ++iF)
                        pTarget[iF] += sfr*workingData[offset + iF];
                }
                #ifdef TIMING
                    report(iG, nGal);
                #endif
            }
        }
        free(readyData);
        free(workingData);
    }
    free(spectra->integrated);

    // Add Lyman absorption to full spectra if applicable
    add_IGM_absorption_spectra(spectra, target, galParams->nGal);

    #ifdef TIMING
        profiler_end(SUM);
    #endif
}


double *composite_spectra_cext(sed_params_t *spectra,
                               gal_params_t *galParams, dust_params_t *dustParams,
                               short outType, short approx, short nThread, int population) {
    // Trim the metallicity of each SSP such that it is within the range of
    // input SED templates
    if (population == 2)
        trim_gal_params(galParams, spectra->minZ, spectra->maxZ);

    #ifdef TIMING
        init_profiler();
        timing_start("Compute magnitudes");
    #endif

    // Initialise outputs
    int iGF;
    int nGal = galParams->nGal;
    int nFlux = spectra->nFlux;
    int nGF = nGal*nFlux;
    double *output = malloc(nGF*sizeof(double));
    double *pOutput = output;
    for(iGF = 0; iGF < nGF; ++iGF)
        *pOutput++ = TOL;

    // Compute spectra for given inputs
    //   -Only use approximation when considering dust
    compute_spectra(output, spectra, galParams, dustParams,
                    (short)(approx && dustParams != NULL), nThread, population);

    if (outType == 0) {
        // Convert to AB magnitude
        pOutput = output;
        for(iGF = 0; iGF < nGF; ++iGF) {
            *pOutput = M_AB(*pOutput);
            ++pOutput;
        }
        #ifdef TIMING
            timing_end();
            profiler_summary();
        #endif
        return output;
    }
    else if (outType == 1) {
        #ifdef TIMING
            timing_end();
            profiler_summary();
        #endif
        return output;
    }

    // Fit UV slopes
    int iG;
    int nR = 3;
    output = (double*)realloc(output, (nFlux + nR)*nGal*sizeof(double));
    fit_UV_slope(output + nFlux*nGal, output, nGal, nFlux, spectra->logWaves, nFlux - 1, nR);

    // Convert to AB magnitude
    pOutput = output + nFlux - 1;
    for(iG = 0; iG < nGal; ++iG) {
        *pOutput = M_AB(*pOutput);
        pOutput += nFlux;
    }
    #ifdef TIMING
        timing_end();
        profiler_summary();
    #endif
    return output;
}
