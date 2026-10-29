#include "noahmp_stomata.h"
#include "ball_berry_canopy_physics.h"
#include <math.h>
#include <stddef.h>

/* Direct translation of STOMATA in supplied module_sf_noahmplsm(3).F.
 * Retain Noah-MP's three CI iterations and algebraic quadratic root.
 * FLO 10/21: BTRAN (in->btran, root stress factor 0 to 1) scales VCMX as in
 * Noah-MP: VCMX = VCMX25/F2(TC)*FNF*BTRAN*F1(AVCMX,TC). BTRAN = 1 gives the
 * unstressed (potential) case; BTRAN < 1 is NOT a linear reduction of
 * transpiration, because it acts through photosynthesis and the Ball-Berry
 * conductance, and the leaf energy balance responds.
 * This isolated port is NOT a full canopy or PET calculation.
 * FLO 10/28: in the form Tapenade differentiates (dcfe-noahmp-stack): one return (a status variable),
 * bb_is_finite() for isfinite(), no ternaries, the nested fmin and fmax written as two statements, and
 * every chain of two or more multiplications and divisions written one operation per statement, in C's
 * left-to-right order (Tapenade regroups chains, and its primal would differ in the last bits).
 * Byte-identical.
 */
int noahmp_stomata(const NoahmpStomataParameters *p,
                            const NoahmpStomataInputs *in,
                            NoahmpStomataOutputs *out)
{
    double cf, fnf, tc, ppf, j, kc, ko, awc, cp, vcmx;
    double ci, rlb, cea, wj, wc, we, psn, cs, a, b, c, q, r1, r2, rs;
    double nitrogen_ratio, smaller_vapor_pressure, smaller_rate;
    double product, second_product, denominator, exponent_argument, discriminant;
    int iter, status;
    status = 0;
    if (!p || !in || !out) status = -1;
    if (status == 0) {
        if (in->sfctmp <= 0.0 || in->sfcpres <= 0.0 ||
            in->mpe <= 0.0 || p->bp <= 0.0 || p->folnmx < 0.0 ||
            p->ko25 <= 0.0 || p->kc25 <= 0.0 || p->akc <= 0.0 ||
            p->ako <= 0.0 || p->avcmx <= 0.0 || p->c3psn < 0.0 ||
            p->c3psn > 1.0 || in->ei <= 0.0 || in->rb < 0.0 ||
            !(in->btran >= 0.0 && in->btran <= 1.0))
            status = -1;
    }
    if (status == 0) {
        denominator = 8.314 * in->sfctmp;
        cf = in->sfcpres / denominator;
        cf = cf * 1.0e6;
        rs = cf / p->bp;
        psn = 0.0;
        if (in->apar <= 0.0) {
            out->rs = rs;
            out->psn = psn;
            if (!bb_is_finite(rs)) status = -1;
        } else {
            nitrogen_ratio = in->foln / fmax(in->mpe, p->folnmx);
            fnf = fmin(nitrogen_ratio, 1.0);
            tc = in->tv - 273.16;
            ppf = 4.6 * in->apar;
            j = ppf * p->qe25;
            kc = p->kc25 * pow(p->akc, (tc - 25.0) / 10.0);
            ko = p->ko25 * pow(p->ako, (tc - 25.0) / 10.0);
            awc = kc * (1.0 + in->o2 / ko);
            cp = 0.5 * kc;
            cp = cp / ko;
            cp = cp * in->o2;
            cp = cp * 0.21;
            denominator = 8.314 * (tc + 273.16);
            exponent_argument = (-2.2e5 + 710.0 * (tc + 273.16)) / denominator;
            vcmx = p->vcmx25 / (1.0 + exp(exponent_argument));
            vcmx = vcmx * fnf;
            vcmx = vcmx * in->btran;
            vcmx = vcmx * pow(p->avcmx, (tc - 25.0) / 10.0);
            product = 0.7 * in->co2;
            product = product * p->c3psn;
            second_product = 0.4 * in->co2;
            second_product = second_product * (1.0 - p->c3psn);
            ci = product + second_product;
            rlb = in->rb / cf;
            smaller_vapor_pressure = fmin(in->ea, in->ei);
            product = 0.25 * in->ei;
            product = product * p->c3psn;
            second_product = 0.40 * in->ei;
            second_product = second_product * (1.0 - p->c3psn);
            cea = fmax(product + second_product, smaller_vapor_pressure);
            for (iter = 0; iter < 3; ++iter) {
                product = fmax(ci - cp, 0.0) * j;
                product = product / (ci + 2.0 * cp);
                product = product * p->c3psn;
                wj = product + j * (1.0 - p->c3psn);
                product = fmax(ci - cp, 0.0) * vcmx;
                product = product / (ci + awc);
                product = product * p->c3psn;
                wc = product + vcmx * (1.0 - p->c3psn);
                product = 0.5 * vcmx;
                product = product * p->c3psn;
                second_product = 4000.0 * vcmx;
                second_product = second_product * ci;
                second_product = second_product / in->sfcpres;
                second_product = second_product * (1.0 - p->c3psn);
                we = product + second_product;
                smaller_rate = fmin(wc, we);
                psn = fmin(wj, smaller_rate) * in->igs;
                product = 1.37 * rlb;
                product = product * in->sfcpres;
                product = product * psn;
                cs = fmax(in->co2 - product, in->mpe);
                product = p->mp * psn;
                product = product * in->sfcpres;
                product = product * cea;
                denominator = cs * in->ei;
                product = product / denominator;
                a = product + p->bp;
                product = p->mp * psn;
                product = product * in->sfcpres;
                product = product / cs;
                b = (product + p->bp) * rlb - 1.0;
                c = -rlb;
                product = 4.0 * a;
                product = product * c;
                discriminant = b*b - product;
                if (b >= 0.0)
                    q = -0.5 * (b + sqrt(discriminant));
                else
                    q = -0.5 * (b - sqrt(discriminant));
                r1 = q / a;
                r2 = c / q;
                rs = fmax(r1, r2);
                product = psn * in->sfcpres;
                product = product * 1.65;
                product = product * rs;
                ci = fmax(cs - product, 0.0);
            }
            out->rs = rs * cf;
            out->psn = psn;
            if (!(bb_is_finite(out->rs) && bb_is_finite(out->psn))) status = -1;
        }
    }
    return status;
}
