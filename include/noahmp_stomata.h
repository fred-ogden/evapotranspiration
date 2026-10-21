#ifndef NOAHMP_STOMATA_H
#define NOAHMP_STOMATA_H
/* Standalone C99 translation of Noah-MP STOMATA (Ball-Berry stomatal resistance and
 * Collatz photosynthesis). FLO 10/21: BTRAN is now an input; it scales VCMX exactly as
 * in Noah-MP. In PET method 6 it is the BMI input root_stress_factor_0_to_1. */
typedef struct {
    double bp, folnmx, qe25, kc25, akc, ko25, ako;
    double vcmx25, avcmx, c3psn, mp;
} NoahmpStomataParameters;
typedef struct {
    double mpe, apar, foln, tv, ei, ea, sfctmp, sfcpres;
    double o2, co2, igs, rb;
    double btran;   /* soil moisture (root) stress factor [0,1]; FLO 10/21 */
} NoahmpStomataInputs;
typedef struct {
    double rs;  /* leaf stomatal resistance, s/m */
    double psn; /* photosynthesis, umol CO2 / m2 / s */
} NoahmpStomataOutputs;
/* Returns zero on success; -1 on invalid input/nonfinite results. */
int noahmp_stomata(const NoahmpStomataParameters *p,
                            const NoahmpStomataInputs *in,
                            NoahmpStomataOutputs *out);
#endif
