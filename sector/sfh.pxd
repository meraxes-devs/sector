cdef extern from "clib/sector.h":
    ctypedef struct ssp_t:
        short index
        float metals
        float sfr

    ctypedef struct csp_t:
        ssp_t *bursts
        int nBurst

    ctypedef struct gal_params_t:
        double z
        int nAgeStep
        double *ageStep
        int nGal
        int *indices
        csp_t *histories

cdef int *init_1d_int(int[:] memview)
cdef float *init_1d_float(float[:] memview)
cdef double *init_1d_double(double[:] memview)
cdef void read_gal_params(gal_params_t *galParams, char *fname)
cdef void free_gal_params(gal_params_t *galParams)

cdef class galaxy_tree_meraxes:
    cdef:
        object fname
        double h
        # Poniters to store data
        int **firstProgenitor
        int **nextProgenitor
        float **metals
        float **sfr
        # Trace parameters
        int tSnap
        ssp_t *bursts
        int nBurst
        #
        int snapMin
        int snapMax

    cdef void trace_progenitors(self, int snap, int galIdx)
    cdef csp_t *trace_properties(self, int tSnap, int[:] indices)

cdef class stellar_population:
    cdef:
        gal_params_t gp
        int nDfStep
        double *dfStep
        csp_t *dfH
        object data

    cdef void _update_age_step(self, double[:] newStep)
    cdef _new_age_step(self, int nAvg)
    cdef void _average_csp(self, csp_t *newH, csp_t *gpH, int nMax, int nAvg, double[:] newStep)
    cdef build_data(self)
    cdef gal_params_t *pointer(self)
    cdef void _reset_gp(self)
