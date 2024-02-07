# cython: c_string_encoding = ascii

import numpy as np
from time import time
from pandas import DataFrame
from dragons import meraxes


global sTime

def timing_start(text):
    global sTime
    sTime = time()
    print "#***********************************************************"
    print text


def timing_end():
    global sTime
    elapsedTime = time() - sTime
    minute = int(elapsedTime)/60
    print "# Done!"
    print "# Elapsed time: %i min %.6f sec"%(minute, elapsedTime - minute*60)
    print "#***********************************************************\n"


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#                                                                               #
# Primary functions                                                             #
#                                                                               #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
DEF MAX_NODE = 100000

cdef class galaxy_tree_meraxes:
    def __cinit__(self, fname, int snapMax, double h):
        #=====================================================================
        # Load model output
        #=====================================================================
        self.fname = fname
        self.h = h
        cdef:
            int snapNum = snapMax+ 1
            int snapMin = snapMax
            int snap, N
            int[:] intMemview1, intMemview2
            float[:] floatMemview1, floatMemview2
        timing_start("# Read meraxes output")
        self.firstProgenitor = <int**>malloc(snapNum*sizeof(int*))
        self.nextProgenitor = <int**>malloc(snapMax*sizeof(int*))
        # Unit: 1e10 M_sun (New metallicity tracer)
        self.metals = <float**>malloc(snapNum*sizeof(float*))
        # Unit: M_sun/yr
        self.sfr = <float**>malloc(snapNum*sizeof(float*))
        meraxes.set_little_h(h = h)
        for snap in xrange(snapMax, -1, -1):
            try:
                # Copy metallicity and star formation rate to the pointers
                gals = meraxes.io.read_gals(fname, snap,
                                            props = ["ColdGas", "MetalsColdGas", "Sfr"])
                print ''
                metals = gals["MetalsColdGas"]/gals["ColdGas"]
                metals[np.isnan(metals)] = 0.001
                self.metals[snap] = init_1d_float(metals)
                self.sfr[snap] = init_1d_float(gals["Sfr"])
                snapMin = snap
                gals = None
            except IndexError:
                print "# No galaxies in snapshot %d"%snap
                break;
        print "# snapMin = %d"%snapMin
        for snap in xrange(snapMin, snapNum):
            # Copy first progenitor indices to the pointer
            self.firstProgenitor[snap] = \
            init_1d_int(meraxes.io.read_firstprogenitor_indices(fname, snap))
            # Copy next progenitor indices to the pointer
            if snap < snapMax:
                self.nextProgenitor[snap] = \
                init_1d_int(meraxes.io.read_nextprogenitor_indices(fname, snap))
        self.snapMin = snapMin
        self.snapMax = snapMax
        timing_end()
        # This varible is used to trace progenitors
        self.bursts = <ssp_t*>malloc(MAX_NODE*sizeof(ssp_t))


    def __dealloc__(self):
        cdef int iS
        # Free nextProgenitor. There is no indices in nextProgenitor[snapMax]
        for iS in xrange(self.snapMin, self.snapMax):
            free(self.nextProgenitor[iS])
        # Free other pointers
        for iS in xrange(self.snapMin, self.snapMax + 1):
            free(self.firstProgenitor[iS])
            free(self.metals[iS])
            free(self.sfr[iS])
        free(self.firstProgenitor)
        free(self.nextProgenitor)
        free(self.metals)
        free(self.sfr)
        # This varible is used to trace progenitors
        free(self.bursts)


    cdef void trace_progenitors(self, int snap, int galIdx):
        cdef:
            float sfr
            ssp_t *pBursts
            int nProg
        if galIdx >= 0:
            sfr = self.sfr[snap][galIdx]
            if sfr > 0.:
                self.nBurst += 1
                nProg = self.nBurst
                if (nProg >= MAX_NODE):
                    raise MemoryError("Number of progenitors exceeds MAX_NODE")
                pBursts = self.bursts + nProg
                pBursts.index = self.tSnap - snap
                pBursts.metals = self.metals[snap][galIdx]
                pBursts.sfr = sfr
            self.trace_progenitors(snap - 1, self.firstProgenitor[snap][galIdx])
            self.trace_progenitors(snap, self.nextProgenitor[snap][galIdx])


    cdef csp_t *trace_properties(self, int tSnap, int[:] indices):
        cdef:
            int iG
            int nGal = indices.shape[0]
            csp_t *histories = <csp_t*>malloc(nGal*sizeof(csp_t))
            csp_t *pHistories
            ssp_t *bursts = self.bursts
            int galIdx
            int nProg
            float sfr
            size_t memSize
            size_t totalMemSize = 0
        timing_start("# Read galaxies properties")
        self.tSnap = tSnap
        for iG in xrange(nGal):
            galIdx = indices[iG]
            nProg = -1
            sfr = self.sfr[tSnap][galIdx]
            if sfr > 0.:
                nProg += 1
                bursts.index = 0
                bursts.metals = self.metals[tSnap][galIdx]
                bursts.sfr = sfr
            self.nBurst = nProg
            self.trace_progenitors(tSnap - 1, self.firstProgenitor[tSnap][galIdx])
            nProg = self.nBurst + 1
            pHistories = histories + iG
            pHistories.nBurst = nProg
            if nProg == 0:
                pHistories.bursts = NULL
                #print "Warning: snapshot %d, index %d"%(tSnap, galIdx)
                #print "         the star formation rate is zero throughout the histroy"
            else:
                memSize = nProg*sizeof(ssp_t)
                pHistories.bursts = <ssp_t*>malloc(memSize)
                memcpy(pHistories.bursts, bursts, memSize)
                totalMemSize += memSize
        print "# %.1f MB memory has been allocted"%(totalMemSize/1024./1024.)
        timing_end()
        return histories


    def get_galaxy_ID(self, int tSnap, int[:] indices):
        return meraxes.io.read_gals(
            self.fname, tSnap, props = ["ID"], indices = indices, quiet = True
        )["ID"]


cdef void free_csp(csp_t *histories, int nGal):
    cdef int iG
    for iG in xrange(nGal):
        free(histories[iG].bursts)


cdef void copy_csp(csp_t *newH, csp_t *gpH, int nGal):
    cdef:
        int iG
        int nB
        ssp_t *bursts = NULL

    for iG in xrange(nGal):
        nB = gpH[iG].nBurst
        newH[iG].nBurst = nB

        bursts = <ssp_t*>malloc(nB*sizeof(ssp_t))
        memcpy(bursts, gpH[iG].bursts, nB*sizeof(ssp_t))
        newH[iG].bursts = bursts


cdef void save_gal_params(gal_params_t *galParams, char *fname):
    cdef:
        int iA, iG
        FILE *fp

        double z = galParams.z
        int nAgeStep = galParams.nAgeStep
        double *ageStep = galParams.ageStep
        int nGal = galParams.nGal
        int *indices = galParams.indices
        csp_t *histories = galParams.histories
        llong_t *ids = galParams.ids

        int nBurst

    fp = fopen(fname, 'wb')
    # Write redshift
    fwrite(&z, sizeof(double), 1, fp)
    # Write ageStep
    fwrite(&nAgeStep, sizeof(int), 1, fp)
    fwrite(ageStep, sizeof(double), nAgeStep, fp)
    # Write indices
    fwrite(&nGal, sizeof(int), 1, fp)
    fwrite(indices, sizeof(int), nGal, fp)
    # Write histories
    for iG in xrange(nGal):
        nBurst = histories[iG].nBurst
        fwrite(&nBurst, sizeof(int), 1, fp)
        fwrite(histories[iG].bursts, sizeof(ssp_t), nBurst, fp)
    # Write ids
    fwrite(ids, sizeof(llong_t), nGal, fp)
    fclose(fp)


cdef void read_gal_params(gal_params_t *galParams, char *fname):
    cdef:
        int iG
        FILE *fp

        double z
        int nAgeStep
        double *ageStep
        int nGal
        int *indices
        csp_t *histories
        llong_t *ids

        int nBurst

    timing_start("# Read galaxy properties")
    fp = fopen(fname, 'rb')
    if fp == NULL:
        raise IOError("Fail to open the input file")
    # Read redshift
    fread(&z, sizeof(double), 1, fp)
    # Read ageStep
    fread(&nAgeStep, sizeof(int), 1, fp)
    ageStep = <double*>malloc(nAgeStep*sizeof(double))
    fread(ageStep, sizeof(double), nAgeStep, fp)
    # Read indices
    fread(&nGal, sizeof(int), 1, fp)
    indices = <int*>malloc(nGal*sizeof(int))
    fread(indices, sizeof(int), nGal, fp)
    # Read histories
    histories = <csp_t*>malloc(nGal*sizeof(csp_t))
    pHistories = histories
    for iG in xrange(nGal):
        fread(&nBurst, sizeof(int), 1, fp)
        histories[iG].nBurst = nBurst
        histories[iG].bursts = <ssp_t*>malloc(nBurst*sizeof(ssp_t))
        fread(histories[iG].bursts, sizeof(ssp_t), nBurst, fp)
        pHistories += 1
    # Read ids
    ids = <llong_t*>malloc(nGal*sizeof(llong_t))
    fread(ids, sizeof(llong_t), nGal, fp)
    fclose(fp)

    galParams.z = z
    galParams.nAgeStep = nAgeStep
    galParams.ageStep = ageStep
    galParams.nGal = nGal
    galParams.indices = indices
    galParams.histories = histories
    galParams.ids = ids

    timing_end()


cdef void free_gal_params(gal_params_t *galParams):
    free_csp(galParams.histories, galParams.nGal)
    free(galParams.histories)
    free(galParams.ageStep)
    free(galParams.indices)
    free(galParams.ids)


cdef class stellar_population:
    property ID:
        def __get__(self):
            return np.array(<llong_t[:self.gp.nGal]>self.gp.ids)


    property indices:
        def __get__(self):
            return np.array(<int[:self.gp.nGal]>self.gp.indices)


    property timeStep:
        def __get__(self):
            return np.array(<double[:self.gp.nAgeStep]>self.gp.ageStep)


    property z:
        def __get__(self):
            return self.gp.z


    cdef void _update_age_step(self, double[:] newStep):
        cdef:
            int nAgeStep = newStep.shape[0]
            double *ageStep = <double*>malloc(nAgeStep*sizeof(double))

        memcpy(ageStep, &newStep[0], nAgeStep*sizeof(double))
        free(self.gp.ageStep)
        self.gp.ageStep = ageStep
        self.gp.nAgeStep = nAgeStep


    cdef _new_age_step(self, int nAvg):
        cdef:
            int iA
            int iNS = 0
            int nDfStep = self.nDfStep
            double *dfStep = self.dfStep
            int nA = nDfStep/nAvg if nDfStep%nAvg == 0 else nDfStep/nAvg + 1
            double[:] newStep = np.zeros(nA)

        for iA in xrange(nDfStep):
            if (iA + 1)%nAvg == 0 or iA == nDfStep - 1:
                newStep[iNS] = dfStep[iA]
                iNS += 1
        return np.asarray(newStep)


    cdef void _average_csp(self, csp_t *newH, csp_t *gpH, int nMax, int nAvg, double[:] newStep):
        cdef:
            int iA
            int nDfStep = self.nDfStep
            double *dfStep = self.dfStep
            double[:] dfInterval = np.zeros(nDfStep)

        dfInterval[0] = dfStep[0]
        for iA in xrange(1, nDfStep):
            dfInterval[iA] = dfStep[iA] - dfStep[iA - 1]

        cdef double[:] timeInterval = np.zeros(nMax)

        timeInterval[0] = newStep[0]
        for iA in xrange(1, nMax):
            timeInterval[iA] = newStep[iA] - newStep[iA - 1]

        cdef:
            int iB, iNB
            int iLow = 0
            int iHigh = nAvg
            int nB = gpH.nBurst
            int nNB = 0
            ssp_t *bursts = gpH.bursts
            ssp_t *tmpB = <ssp_t*>malloc(nMax*sizeof(ssp_t))
            ssp_t *newB = NULL
            int index
            double sfr, metals, dm, dt

        for iNB in xrange(nMax):
            sfr, metals, dm = 0., 0., 0.
            dt = timeInterval[iNB]
            for iB in xrange(nB):
                index = bursts[iB].index
                if index >= iLow and index < iHigh:
                    dm = bursts[iB].sfr*dfInterval[index]
                    sfr += dm
                    metals += bursts[iB].metals*dm
            if sfr != 0.:
                tmpB[nNB].index = iNB
                tmpB[nNB].metals = metals/sfr
                tmpB[nNB].sfr = sfr/dt
                nNB += 1

            iLow += nAvg
            iHigh += nAvg
            if iLow >= nDfStep:
                break
            if iHigh > nDfStep:
                iHigh = nDfStep

        if nNB > 0:
            newB = <ssp_t*>malloc(nNB*sizeof(ssp_t))
            memcpy(newB, tmpB, nNB*sizeof(ssp_t))
            newH.nBurst = nNB
            newH.bursts = newB
        else:
            newH.nBurst = 0
            newH.bursts = NULL
        free(tmpB)


    cdef _moment(self, double maxAge, int order):
        """Compute the n-th moment of the star formation history.

           Parameters
           ----------
           maxAge: double
               Upper limit of the integral. It should be in a unit of yr.
           order: int
               Order of the moment.

           Returns
           -------
           res: ndarray
               A 1-D array which stores the results.
        """
        cdef:
            int iG, iB

            int nAgeStep = self.gp.nAgeStep
            double *ageStep = self.gp.ageStep
            int nGal = self.gp.nGal
            csp_t *gpH = self.gp.histories
            int nBurst
            ssp_t *bursts
            int index

            int slope = order + 1
            double tLower
            double tUpper

            double[:] res = np.zeros(nGal)

        if (maxAge > ageStep[nAgeStep - 1]):
            raise ValueError("Maximum timestep is exceeded!")

        for iG in xrange(nGal):
            nBurst = gpH[iG].nBurst
            bursts = gpH[iG].bursts
            I = 0.
            for iB in xrange(nBurst):
                index = bursts[iB].index
                if index == 0:
                    tLower = 0.
                else:
                    tLower = ageStep[index - 1]
                tUpper = ageStep[index]
                if maxAge < tLower:
                    continue
                elif maxAge < tUpper:
                    tUpper = maxAge
                I += bursts[iB].sfr*(tUpper**slope - tLower**slope)/slope
            res[iG] = I

        return np.asarray(res)


    cdef void _build_data(self):
        cdef:
            int iG, iB
            int nGal = self.gp.nGal
            int nB

            csp_t *pH = self.gp.histories
            ssp_t *pB = NULL

        data = np.empty(nGal, dtype = object)
        for iG in xrange(nGal):
            nB = pH.nBurst
            pB = pH.bursts
            arr = np.zeros(nB, dtype = [('index', 'i4'), ('metallicity', 'f8'), ('sfr', 'f8')])
            for iB in xrange(nB):
                arr[iB] = pB.index, pB.metals, pB.sfr
                pB += 1
            data[iG] = arr[np.argsort(arr["index"])]
            pH += 1
        self.data = data

    cdef void _reset_gp(self):
        cdef gal_params_t *gp = &self.gp
        free_csp(gp.histories, gp.nGal)
        free(gp.ageStep)
        gp.nAgeStep = self.nDfStep
        gp.ageStep = <double*>malloc(gp.nAgeStep*sizeof(double))
        memcpy(gp.ageStep, self.dfStep, gp.nAgeStep*sizeof(double))
        self.data = None


    cdef gal_params_t *pointer(self):
        return &self.gp


    def save(self, name):
        save_gal_params(self.pointer(), name)


    def reconstruct(self, timeGrid = 1):
        if timeGrid >= 0 and timeGrid < self.nDfStep:
            timeGrid = int(timeGrid)
            self._reset_gp()
        else:
            raise ValueError("timeGrid should be between 0 and %d!"%self.nDfStep)

        cdef:
            int iG
            int nGal = self.gp.nGal
            int nAgeStep = 0
            csp_t *newH = self.gp.histories
            csp_t *dfH = self.dfH
            double[:] newStep

        if timeGrid == 0:
            copy_csp(newH, dfH, nGal)
            return
        if timeGrid == 1:
            nAgeStep = self.nDfStep
            newStep = <double[:nAgeStep]>self.dfStep
        elif timeGrid > 1:
            newStep = self._new_age_step(timeGrid)
            nAgeStep = newStep.shape[0]
            self._update_age_step(newStep)
        for iG in xrange(nGal):
            self._average_csp(newH + iG, dfH + iG, nAgeStep, timeGrid, newStep)


    def mean_SFR(self, meanAge = 100.):
        meanAge *= 1e6 # Convert Myr to yr
        return self._moment(meanAge, 0)/meanAge


    def __getitem__(self, idx):
        if self.data is None:
            self._build_data()
        return self.data[idx]


    def __cinit__(self, gals, snapshot = None, indices = None):
        cdef:
            gal_params_t *gp = &self.gp
            galaxy_tree_meraxes galData = None

        if isinstance(gals, str):
            # Read SFHs from files
            read_gal_params(gp, gals)
        else:
            # Read SFHs from meraxes outputs
            galData = <galaxy_tree_meraxes>gals
            # Read redshift
            gp.z = meraxes.io.grab_redshift(galData.fname, snapshot)
            # Read lookback time
            gp.nAgeStep = snapshot
            timeStep = meraxes.io.read_snaplist(galData.fname, galData.h)[2]*1e6 # Convert Myr to yr
            ageStep = np.zeros(snapshot, dtype = 'f8')
            for iA in xrange(snapshot):
                ageStep[iA] = timeStep[snapshot - iA - 1] - timeStep[snapshot]
            gp.ageStep = init_1d_double(ageStep)
            # Store galaxy indices
            indices = np.asarray(indices, dtype = 'i4')
            gp.nGal = len(indices)
            gp.indices = init_1d_int(indices)
            # Read SFHs
            gp.histories = galData.trace_properties(snapshot, indices)
            # Read galaxy IDs
            gp.ids = init_1d_llong(galData.get_galaxy_ID(snapshot, indices))
        #
        self.dfH = <csp_t*>malloc(gp.nGal*sizeof(csp_t))
        copy_csp(self.dfH, gp.histories, gp.nGal)
        #
        self.nDfStep = gp.nAgeStep
        self.dfStep = <double*>malloc(gp.nAgeStep*sizeof(double))
        memcpy(self.dfStep, gp.ageStep, gp.nAgeStep*sizeof(double))
        #
        self.data = None


    def __dealloc__(self):
        free_gal_params(&self.gp)
        free_csp(self.dfH, self.gp.nGal)
        free(self.dfH)
        free(self.dfStep)
