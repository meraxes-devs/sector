#define _DUST_
#include"tools.h"
#include"sector.h"

int birth_cloud_interval(double tBC, double *ageStep, int nAgeStep) {
    if (tBC >= ageStep[nAgeStep - 1])
        return nAgeStep;
    else if (tBC < ageStep[0])
        return 0;
    else
        return bisection_search(tBC, ageStep, nAgeStep) + 1;
}


void init_templates_special(sed_params_t *spectra, double tBC, int approx) {
    /* Special templates are for birth cloud 
     * Dimension: metallicity × wavelength
     */

    // Find the time inverval containning the birth cloud
    int nZ = spectra->nZ;
    int nWaves = spectra->nWaves;
    int nZW = nZ*nWaves;
    double *age = spectra->age;
    int nAgeStep = spectra->nAgeStep;
    double *ageStep = spectra->ageStep;
    int iAgeBC = birth_cloud_interval(tBC, ageStep, nAgeStep);
    double t0, t1;

    if (iAgeBC == nAgeStep) {
        spectra->inBC = (double*)calloc(nZW, sizeof(double));
        spectra->outBC = (double*)calloc(nZW, sizeof(double));
        return;
    }
    else if (iAgeBC == 0) {
        t0 = age[0];
        t1 = ageStep[0];
        if (tBC < t0)
            tBC = t0;
    }
    else {
        t0 = ageStep[iAgeBC - 1];
        t1 = ageStep[iAgeBC];
    }

    // Integrate special templates over the time step of the birth cloud
    int iZW;
    int nAge = spectra->nAge;
    double *rawData = spectra->raw;
    double *inBC = malloc(nZW*sizeof(double));
    double *outBC = malloc(nZW*sizeof(double));

    for(iZW = 0; iZW < nZW; ++iZW) {
        inBC[iZW] = trapz_table(rawData + iZW*nAge, age, nAge, t0, tBC);
        outBC[iZW] = trapz_table(rawData + iZW*nAge, age, nAge, tBC, t1);
    }
    if (!approx) {
        spectra->inBC = inBC;
        spectra->outBC = outBC;
        return;
    }

    // Intgrate special templates over filters
    int iF, iZ;
    int nFlux = spectra->nFlux;
    double *refInBC = malloc(nFlux*nZ*sizeof(double));
    double *refOutBC = malloc(nFlux*nZ*sizeof(double));
    double *pInBC = refInBC;
    double *pOutBC = refOutBC;

    int nFW;
    int *nFilterWaves = spectra->nFilterWaves;
    double *pFilterWaves = spectra->filterWaves;
    double *pFilters = spectra->filters;
    double *waves = spectra->waves;

    for(iF = 0; iF < nFlux; ++iF) {
        nFW = nFilterWaves[iF];
        for(iZ = 0; iZ < nZ; ++iZ) {
            pInBC[iZ] = trapz_filter(pFilters, pFilterWaves, nFW,
                                               inBC + iZ*nWaves, waves, nWaves);
            pOutBC[iZ] = trapz_filter(pFilters, pFilterWaves, nFW,
                                                outBC + iZ*nWaves, waves, nWaves);
        }
        pFilterWaves += nFW;
        pFilters += nFW;
        pInBC += nZ;
        pOutBC += nZ;
    }
    free(inBC);
    free(outBC);

    // Interploate special templates along metallicities
    int minZ = spectra->minZ;
    int nMaxZ = spectra->nMaxZ;
    double *Z = spectra->Z;
    double interpZ;

    inBC = (double*)malloc(nMaxZ*nFlux*sizeof(double));
    outBC = (double*)malloc(nMaxZ*nFlux*sizeof(double));
    pInBC = inBC;
    pOutBC = outBC;
    for(iZ = 0; iZ < nMaxZ; ++iZ) {
        interpZ = (minZ + iZ + 1.)/2000.;
        for(iF = 0; iF < nFlux; ++iF) {
            pInBC[iF] = interp(interpZ, Z, refInBC + iF*nZ, nZ);
            pOutBC[iF] = interp(interpZ, Z, refOutBC + iF*nZ, nZ);
        }
        pInBC += nFlux;
        pOutBC += nFlux;
    }
    spectra->inBC = inBC;
    spectra->outBC = outBC;
    free(refInBC);
    free(refOutBC);
}


void dust_absorption(sed_params_t *spectra, dust_params_t *dustParams) {
    /* tBC:   life time of the birth clound
     * nu:    fraction of ISM dust absorption
     * tauUV: UV-band absorption optical depth
     * nBC:   power law index of tauBC
     * nISM:  power law index of tauISM
     *
     * Reference: da Cunha et al. 2008
     */
    int iA, iW, iZ;
    int nZ = spectra->nZ;
    int nWaves = spectra->nWaves;
    double *waves = spectra->waves;
    int nAgeStep = spectra->nAgeStep;
    double *pData = spectra->ready;
    double *pInBC = spectra->inBC;
    double *pOutBC = spectra->outBC;

    double tauUV_ISM = dustParams->tauUV_ISM;
    double nISM = dustParams->nISM;
    double tauUV_BC = dustParams->tauUV_BC;
    double nBC = dustParams->nBC;
    int iAgeBC = birth_cloud_interval(dustParams->tBC, spectra->ageStep, nAgeStep);
    double *transISM = malloc(nWaves*sizeof(double));
    double *transBC = malloc(nWaves*sizeof(double));
    double ratio;

    // Compute the optical depth of both birth cloud and ISM
    for(iW = 0; iW < nWaves; ++iW) {
        ratio = waves[iW]/1600.;
        transISM[iW] = exp(-tauUV_ISM*pow(ratio, nISM));
        transBC[iW] = exp(-tauUV_BC*pow(ratio, nBC));
    }

    for(iZ = 0; iZ < nZ; ++iZ) {
        for(iA = 0; iA < nAgeStep; ++iA) {
            // Apply optical depth of birth cloud
            if (iA < iAgeBC)
                for(iW = 0; iW < nWaves; ++iW)
                    pData[iW] *= transBC[iW];
            else if (iA == iAgeBC)
                for(iW = 0; iW < nWaves; ++iW)
                    pData[iW] = transBC[iW]*pInBC[iW] + pOutBC[iW];
            // Apply optical depth of ISM
            for(iW = 0; iW < nWaves; ++iW)
                pData[iW] *= transISM[iW];
            pData += nWaves;
        }
        pInBC += nWaves;
        pOutBC += nWaves;
    }

    free(transISM);
    free(transBC);
}


void dust_absorption_approx(double *inBCFlux, double *outBCFlux, double *centreWaves, int nFlux,
                            dust_params_t *dustParams) {
    // Apply reddening directly on the flux in fitlers as an approximation
    int iF;
    double tauUV_ISM = dustParams->tauUV_ISM;
    double nISM = dustParams->nISM;
    double tauUV_BC = dustParams->tauUV_BC;
    double nBC = dustParams->nBC;
    double *transISM = malloc(nFlux*sizeof(double));
    double *transBC = malloc(nFlux*sizeof(double));
    double ratio;

    // Compute the optical depth of both birth cloud and ISM
    for(iF = 0; iF < nFlux; ++iF) {
        ratio = centreWaves[iF]/1600.;
        transISM[iF] = exp(-tauUV_ISM*pow(ratio, nISM));
        transBC[iF] = exp(-tauUV_BC*pow(ratio, nBC));
    }

    // Apply optical depth of birth cloud
    for(iF = 0; iF < nFlux; ++iF)
        inBCFlux[iF] *= transBC[iF];

    // Apply optical depth of ISM
    for(iF = 0; iF < nFlux; ++iF) {
        inBCFlux[iF] *= transISM[iF];
        outBCFlux[iF] *= transISM[iF];
    }

    free(transISM);
    free(transBC);
}
