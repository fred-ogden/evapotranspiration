#ifndef NOAHMP_STOMATA_H
#define NOAHMP_STOMATA_H
/* Standalone C99 translation of Noah-MP STOMATA; not yet PET/BMI integrated. */
typedef struct {
    double bp, folnmx, qe25, kc25, akc, ko25, ako;
    double vcmx25, avcmx, c3psn, mp;
} NoahmpStomataParameters;
typedef struct {
    double mpe, apar, foln, tv, ei, ea, sfctmp, sfcpres;
    double o2, co2, igs, rb;
} NoahmpStomataInputs;
typedef struct {
    double rs;  /* leaf stomatal resistance, s/m */
    double psn; /* photosynthesis, umol CO2 / m2 / s */
} NoahmpStomataOutputs;
/* Returns zero on success; -1 on invalid input/nonfinite results. BTRAN=1. */
int noahmp_stomata_potential(const NoahmpStomataParameters *p,
                            const NoahmpStomataInputs *in,
                            NoahmpStomataOutputs *out);
#endif
