#define _SPECTRA_
#include"tools.h"
#include"sector.h"
#include"hdf5.h"
#include"hdf5_hl.h"


int *age_flag(csp_t *histories, int nAgeStep) {
    int iA, iB;

    int *ageFlag = malloc(nAgeStep*sizeof(int));
    int nB = histories->nBurst;
    ssp_t *bursts = histories->bursts;

    for(iA = 0; iA < nAgeStep; ++iA)
        ageFlag[iA] = 1;
    for(iB = 0; iB < nB; ++iB)
        ageFlag[bursts[iB].index] = 0;

    return ageFlag;
}


int *Z_flag(csp_t *histories, int nMaxZ) {
    int iZ, iB;

    int *ZFlag = malloc(nMaxZ*sizeof(int));
    int nB = histories->nBurst;
    ssp_t *bursts = histories->bursts;

    for(iZ = 0; iZ < nMaxZ; ++iZ)
        ZFlag[iZ] = 1;
    for(iB = 0; iB < nB; ++iB)
        ZFlag[(int)(2000*bursts[iB].metals - 0.5)] = 0;

    return ZFlag;
}


void init_templates_raw(sed_params_t *spectra, char *fName) {
    /* File name should be "sed_library.hdf5" and contain:
     * "/metals" Metallicity grid (1-D dataset)
     * "/waves" Wavelength grid [AA] (1-D dataset)
     * "/age" Stellar age grid [yr] (1-D dataset)
     * "/flux" Flux density [erg/s/AA/cm^2] (3-D dataset)
     *
     * Flux densities should be normlised by the surface area of 10 pc
     * sphere. The first, second, and third dimensions should be metallicty,
     * wavelength and stellar age respectively.
     */
    hid_t file_id;
    hsize_t dims[3];

    printf("#***********************************************************\n");
    printf("# Read SED templates\n");
    file_id = H5Fopen(fName, H5F_ACC_RDONLY, H5P_DEFAULT);
    // Read dimensions
    H5LTget_dataset_info(file_id, "/flux", dims, NULL, NULL);
    int nZ = dims[0];
    int nWaves = dims[1];
    int nAge = dims[2];
    // Read metallicity range
    spectra->nZ = nZ;
    spectra->Z = (double*)malloc(nZ*sizeof(double));
    H5LTread_dataset_double(file_id, "/metals", spectra->Z);
    spectra->minZ = (int)(spectra->Z[0]*2000 - .5);
    spectra->maxZ = (int)(spectra->Z[nZ - 1]*2000 - .5);
    spectra->nMaxZ = spectra->maxZ - spectra->minZ + 1;
    printf("# Metallicity range:\n#\t%.4f to %.4f\n", spectra->Z[0], spectra->Z[nZ - 1]);
    // Read wavelength
    spectra->nWaves = nWaves;
    spectra->waves = (double*)malloc(nWaves*sizeof(double));
    H5LTread_dataset_double(file_id, "/waves", spectra->waves);
    printf("# Wavelength range:\n#\t%.1f AA to %.1f AA\n",
           spectra->waves[0], spectra->waves[nWaves - 1]);
    // Read stellar age
    spectra->nAge = nAge;
    spectra->age = (double*)malloc(nAge*sizeof(double));
    H5LTread_dataset_double(file_id, "/age", spectra->age);
    printf("# Stellar age range:\n#\t%.2f Myr to %.2f Myr\n",
           spectra->age[0]*1e-6, spectra->age[nAge - 1]*1e-6);
    // Read flux
    spectra->raw = (double*)malloc(nZ*nWaves*nAge*sizeof(double));
    H5LTread_dataset_double(file_id, "/flux", spectra->raw);
    H5Fclose(file_id);
    printf("#***********************************************************\n\n");
}


void shrink_templates_raw(sed_params_t *spectra, double maxAge) {
    if (spectra->filters == NULL)
        return;

    int iA, iF, iW, iZ;

    int nZ = spectra->nZ;
    int nWaves = spectra->nWaves;
    double *waves = spectra->waves;
    int nAge = spectra->nAge;
    double *age = spectra->age;
    double *raw = spectra->raw;
    int nFlux = spectra->nFlux;
    int nFW;
    int *nFilterWaves = spectra->nFilterWaves;
    double *filterWaves = spectra->filterWaves;
    double *pFilterWaves;

    int inFlag = 0;
    int outFlag;
    int nNewWaves = 0;
    int *wavesIndices = malloc(nWaves*sizeof(int));
    double w;
    int nNewAge = 0;
    double *newWaves;
    double *newRaw;
    double *pNewRaw;

    // Find shrinked wavelength ranges
    printf("#***********************************************************\n");
    printf("# Shrinked wavelength ranges:\n");
    for(iW = 1; iW < nWaves; ++iW) {
        w = waves[iW];
        pFilterWaves = filterWaves;
        outFlag = 1;
        for(iF = 0; iF < nFlux; ++iF) {
            nFW = nFilterWaves[iF];
            if (w > pFilterWaves[0] && w < pFilterWaves[nFW - 1]) {
                outFlag = 0;
                if (!inFlag) {
                    printf("#\t%.1f AA to ", waves[iW - 1]);
                    wavesIndices[nNewWaves++] = iW - 1;
                    wavesIndices[nNewWaves++] = iW;
                    inFlag = 1;
                }
                else {
                    wavesIndices[nNewWaves++] = iW;
                }
                break;
            }
            pFilterWaves += nFW;
        }
        if (inFlag && outFlag) {
            printf("%.1f AA\n", waves[iW]);
            wavesIndices[nNewWaves++] = iW;
            inFlag = 0;
        }
    }
    printf("# Original nWaves: %d\n", nWaves);
    printf("# Shrinked nWaves: %d\n", nNewWaves);
    // Find nNewAge
    printf("#\n# Shrinked age range:\n");
    for(iA = 0; age[iA] < maxAge; ++iA);
    nNewAge = iA + 1;
    spectra->nAge = nNewAge;
    printf("#\t%.2f Myr to %.2f Myr\n", age[0]*1e-6, age[iA]*1e-6);
    printf("# Original nAge: %d\n", nAge);
    printf("# Shrinked nAge: %d\n", nNewAge);
    printf("#***********************************************************\n\n");
    // Construct new wavelengths
    newWaves = (double*)malloc(nNewWaves*sizeof(double));
    for(iW = 0; iW < nNewWaves; ++iW)
        newWaves[iW] = waves[wavesIndices[iW]];
    spectra->nWaves = nNewWaves;
    spectra->waves = newWaves;
    free(waves);
    // Construct new raw templates
    newRaw = (double*)malloc(nZ*nNewWaves*nNewAge*sizeof(double));
    pNewRaw = newRaw;
    for(iZ = 0; iZ < nZ; ++iZ)
        for(iW = 0; iW < nNewWaves; ++iW)
            for(iA = 0; iA < nNewAge; ++iA)
                *pNewRaw++ = raw[(iZ*nWaves + wavesIndices[iW])*nAge + iA];
    spectra->raw = newRaw;
    free(raw);

    free(wavesIndices);
}


void init_filters(sed_params_t *spectra,
                  double *betaBands, int nBeta, double *restBands, int nRest,
                  double *obsTrans, double *obsWaves, int *nObsWaves, int nObs, double z) {
    // Initialise filter parameters
    //   -This function should be called after ``init_templates_raw``
    int iF, iW;
    double *pData = betaBands;

    int nFlux = nBeta + nRest + nObs;
    int nWaves = spectra->nWaves;
    double *waves = spectra->waves;
    int nNewWaves;
    double *newWaves = malloc(nWaves*sizeof(double));
    double w;
    double lower, upper;
    double norm;

    int *nFilterWaves = malloc(nFlux*sizeof(int));
    double *filterWaves = NULL;
    double *filters = NULL;
    int offset = 0;

    double *centreWaves = malloc(nFlux*sizeof(double));
    double *logWaves = malloc(nFlux*sizeof(double));

    // Initialise rest-frame filters and those that are related to the UV slope
    nRest += nBeta;
    for(iF = 0; iF < nRest; ++iF) {
        if (iF == nBeta)
            pData = restBands;
        //   -Construct wavelengths for a new filter
        lower = *pData++;
        upper = *pData++;
        if (lower < waves[0] || upper > waves[nWaves - 1]) {
            printf("Filter wavelength range are beyond SED templates\n");
            exit(0);
        }
        nNewWaves = 1;
        newWaves[0] = lower;
        for(iW = 0; iW < nWaves; ++iW) {
            w = waves[iW];
            if (w > lower && w < upper)
                newWaves[nNewWaves++] = w;
        }
        newWaves[nNewWaves++] = upper;
        nFilterWaves[iF] = nNewWaves;
        filterWaves = realloc(filterWaves, (offset + nNewWaves)*sizeof(double));
        memcpy(filterWaves + offset, newWaves, nNewWaves*sizeof(double));
        //  -Construct transmissions
        filters = realloc(filters, (offset + nNewWaves)*sizeof(double));
        if (iF < nBeta) {
            norm = 1./(upper - lower);
            for(iW = 0; iW < nNewWaves; ++iW)
                *(filters + offset + iW) = norm;
        }
        else {
            norm = 1./log(upper/lower);
            for(iW = 0; iW < nNewWaves; ++iW)
                *(filters + offset + iW) = 3.34e4*newWaves[iW]*norm;
        }
        offset += nNewWaves;
        //
        w = (lower + upper)/2.;
        centreWaves[iF] = w;
        logWaves[iF] = log(w);
    }
    free(newWaves);

    // Initialise observer-frame filters
    if (nObs > 0) {
        memcpy(nFilterWaves + nRest, nObsWaves, nObs*sizeof(int));
        nNewWaves = 0;
        for(iF = 0; iF < nObs; ++iF)
            nNewWaves += nObsWaves[iF];
        filterWaves = realloc(filterWaves, (offset + nNewWaves)*sizeof(double));
        memcpy(filterWaves + offset, obsWaves, nNewWaves*sizeof(double));
        filters = realloc(filters, (offset + nNewWaves)*sizeof(double));
        memcpy(filters + offset, obsTrans, nNewWaves*sizeof(double));
        for(iF = 0; iF < nObs; ++iF) {
            nNewWaves = nObsWaves[iF];
            lower = filterWaves[offset];
            upper = filterWaves[offset + nNewWaves - 1];
            if (lower < waves[0]*(1. + z) || upper > waves[nWaves - 1]*(1. + z)) {
                printf("Filter wavelength range are beyond SED templates\n");
                exit(0);
            }
            //   -Divide the transmission by the wavelength to calculate the normalisation
            for(iW = 0; iW < nNewWaves; ++iW)
                filters[offset + iW] /= filterWaves[offset + iW];
            norm = 1./trapz_table(filters + offset, filterWaves + offset, nNewWaves, lower, upper);
            //   -Compute the pivot wavelength and transform it to the rest-frame
            for(iW = 0; iW < nNewWaves; ++iW)
                filters[offset + iW] *= filterWaves[offset + iW]*filterWaves[offset + iW];
            w = sqrt(norm*trapz_table(filters + offset, filterWaves + offset,
                                      nNewWaves, lower, upper))/(1. + z);
            centreWaves[nRest + iF] = w;
            logWaves[nRest + iF] = log(w);
            //  -Transform the filter such that it can give AB magnitude
            for(iW = 0; iW < nNewWaves; ++iW)
                filters[offset + iW] *= 3.34e4*norm;
            //  -Transform the filter to the rest-frame
            for(iW = 0; iW < nNewWaves; ++iW)
                filterWaves[offset + iW] /= 1. + z;
            offset += nNewWaves;
        }
    }
    else
        //   -Disable IGM absorption if there is no observer-frame filter
        spectra->igm = 0;
    // Set properties
    spectra->nFlux = nFlux;
    spectra->nObs = nObs;
    spectra->nFilterWaves = nFilterWaves;
    spectra->filterWaves = filterWaves;
    spectra->filters = filters;
    spectra->centreWaves = centreWaves;
    spectra->logWaves = logWaves;
    // Add Lyman absorption
    add_IGM_absorption_filters(spectra);
}


void init_templates_integrated(sed_params_t *spectra) {
    int iA, iW, iZ;
    double *pData;

    int nAgeStep = spectra->nAgeStep;
    double *ageStep = spectra->ageStep;
    int nAge = spectra->nAge;
    double *age = spectra->age;
    int nWaves = spectra->nWaves;
    int nZ = spectra->nZ;
    double *data = spectra->raw;
    // Spectra after integration over time
    // The first dimension refers to metallicites and ages
    // The last dimension refers to wavelengths
    double *intData = malloc(nZ*nAgeStep*nWaves*sizeof(double));

    for(iZ = 0; iZ < nZ; ++iZ)
        for(iA = 0; iA < nAgeStep; ++iA) {
            pData = intData + (iZ*nAgeStep + iA)*nWaves;
            for(iW = 0; iW < nWaves; ++iW) {
                if (iA == 0)
                    // The first time step of SED templates is typicall not zero
                    // Here assumes that the templates is zero beween zero
                    // and the first time step
                    pData[iW] = trapz_table(data + (iZ*nWaves + iW)*nAge, age, nAge,
                                            age[0], ageStep[iA]);
                else
                    pData[iW] = trapz_table(data + (iZ*nWaves + iW)*nAge, age, nAge,
                                            ageStep[iA - 1], ageStep[iA]);
            }
        }
    spectra->integrated = intData;
}


void init_templates_working(sed_params_t *spectra, csp_t *pHistories,
                            dust_params_t *dustParams, int iG) {
    int *ageFlag = NULL;
    int *ZFlag = NULL;

    int minZ = spectra->minZ;
    int nMaxZ = spectra->nMaxZ;
    int nZ = spectra->nZ;
    int nWaves = spectra->nWaves;
    int nAgeStep = spectra->nAgeStep;
    int nFlux = spectra->nFlux;
    double *readyData = spectra->ready;

    if (dustParams != NULL) {
        memcpy(readyData, spectra->integrated, nZ*nAgeStep*nWaves*sizeof(double));
        ageFlag = age_flag(pHistories, nAgeStep);
        ZFlag = Z_flag(pHistories, nMaxZ);
        dust_absorption(spectra, dustParams + iG);
    }
    else if (iG == -1) {
        memcpy(readyData, spectra->integrated, nZ*nAgeStep*nWaves*sizeof(double));
        ageFlag = calloc(nAgeStep, sizeof(int));
        ZFlag = calloc(nMaxZ, sizeof(int));
    }
    else
        return;

    int iA, iF, iW, iZ, i, n;
    double *pData;
    // Spectra to be interploated along metallicities
    // Dimension: filter/wavelength × age × metallicity
    double *refSpectra = NULL;

    if (spectra->filters == NULL) {
        double z = spectra->z;

        n = nZ*nAgeStep;
        refSpectra = malloc(nWaves*n*sizeof(double));
        if (spectra->nObs > 0) {
            // Transform everything to observer frame
            // Note the fluxes in this case is a function of wavelength
            // Therefore the fluxes has a factor of 1/(1 + z)
            for(i = 0; i < n; ++i) {
                pData = readyData + i*nWaves;
                for(iW = 0; iW < nWaves; ++iW)
                   pData[iW] /= 1. + z;
            }
        }
        // Tranpose the templates such that the last dimension is the metallicity
        for(i = 0; i < n; ++i) {
            pData = readyData + i*nWaves;
            for(iW = 0; iW < nWaves; ++iW)
                refSpectra[(iW*nAgeStep + i%nAgeStep)*nZ + i/nAgeStep] = pData[iW];
        }

    }
    else {
        // Intgrate SED templates over filters
        int nFW;
        int *nFilterWaves = spectra->nFilterWaves;
        double *pFilterWaves = spectra->filterWaves;
        double *pFilters = spectra->filters;
        double *waves = spectra->waves;

        n = nZ*nAgeStep;
        refSpectra = malloc(nFlux*n*sizeof(double));
        for(iF = 0; iF < nFlux; ++iF) {
            nFW = nFilterWaves[iF];
            for(i = 0; i < n; ++i) {
                iA = i/nZ;
                if (ageFlag[iA])
                    continue;
                refSpectra[iF*n + i] = trapz_filter(pFilters, pFilterWaves, nFW,
                                                    readyData + (i%nZ*nAgeStep + iA)*nWaves,
                                                    waves, nWaves);
            }
            pFilterWaves += nFW;
            pFilters += nFW;
        }
    }
    // Interploate SED templates along metallicities
    double *Z = spectra->Z;
    double *workingData = spectra->working;
    double interpZ;

    n = nMaxZ*nAgeStep;
    for(i = 0; i < n; ++i) {
        iZ = i/nAgeStep;
        if (ZFlag[iZ])
            continue;
        iA = i%nAgeStep;
        if (ageFlag[iA])
            continue;
        interpZ = (minZ + iZ + 1.)/2000.;
        pData = workingData + i*nFlux;
        for(iF = 0; iF < nFlux; ++iF)
            pData[iF] = interp(interpZ, Z, refSpectra + (iF*nAgeStep+ iA)*nZ, nZ);
    }

    free(ageFlag);
    free(ZFlag);
    free(refSpectra);
}


void free_templates_raw(sed_params_t *spectra) {
    free(spectra->Z);
    free(spectra->age);
    free(spectra->waves);
    free(spectra->raw);
}


void free_filters(sed_params_t *spectra) {
    free(spectra->nFilterWaves);
    free(spectra->filterWaves);
    free(spectra->filters);
    free(spectra->centreWaves);
    free(spectra->logWaves);
}
