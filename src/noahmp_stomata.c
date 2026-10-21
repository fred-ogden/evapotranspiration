#include "noahmp_stomata.h"
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
 */
int noahmp_stomata(const NoahmpStomataParameters *p,
                            const NoahmpStomataInputs *in,
                            NoahmpStomataOutputs *out)
{
    double cf, fnf, tc, ppf, j, kc, ko, awc, cp, vcmx;
    double ci, rlb, cea, wj, wc, we, psn, cs, a, b, c, q, r1, r2, rs;
    int iter;
    if (!p || !in || !out || in->sfctmp <= 0.0 || in->sfcpres <= 0.0 ||
        in->mpe <= 0.0 || p->bp <= 0.0 || p->folnmx < 0.0 ||
        p->ko25 <= 0.0 || p->kc25 <= 0.0 || p->akc <= 0.0 ||
        p->ako <= 0.0 || p->avcmx <= 0.0 || p->c3psn < 0.0 ||
        p->c3psn > 1.0 || in->ei <= 0.0 || in->rb < 0.0 ||
        !(in->btran >= 0.0 && in->btran <= 1.0))
        return -1;
    cf = in->sfcpres / (8.314 * in->sfctmp) * 1.0e6;
    rs = cf / p->bp;
    psn = 0.0;
    if (in->apar <= 0.0) {
        out->rs = rs;
        out->psn = psn;
        return isfinite(rs) ? 0 : -1;
    }
    fnf = fmin(in->foln / fmax(in->mpe, p->folnmx), 1.0);
    tc = in->tv - 273.16;
    ppf = 4.6 * in->apar;
    j = ppf * p->qe25;
    kc = p->kc25 * pow(p->akc, (tc - 25.0) / 10.0);
    ko = p->ko25 * pow(p->ako, (tc - 25.0) / 10.0);
    awc = kc * (1.0 + in->o2 / ko);
    cp = 0.5 * kc / ko * in->o2 * 0.21;
    vcmx = p->vcmx25 /
           (1.0 + exp((-2.2e5 + 710.0 * (tc + 273.16)) /
                      (8.314 * (tc + 273.16)))) * fnf * in->btran *
           pow(p->avcmx, (tc - 25.0) / 10.0);
    ci = 0.7 * in->co2 * p->c3psn + 0.4 * in->co2 * (1.0 - p->c3psn);
    rlb = in->rb / cf;
    cea = fmax(0.25 * in->ei * p->c3psn +
               0.40 * in->ei * (1.0 - p->c3psn), fmin(in->ea, in->ei));
    for (iter = 0; iter < 3; ++iter) {
        wj = fmax(ci - cp, 0.0) * j / (ci + 2.0 * cp) * p->c3psn +
             j * (1.0 - p->c3psn);
        wc = fmax(ci - cp, 0.0) * vcmx / (ci + awc) * p->c3psn +
             vcmx * (1.0 - p->c3psn);
        we = 0.5 * vcmx * p->c3psn +
             4000.0 * vcmx * ci / in->sfcpres * (1.0 - p->c3psn);
        psn = fmin(wj, fmin(wc, we)) * in->igs;
        cs = fmax(in->co2 - 1.37 * rlb * in->sfcpres * psn, in->mpe);
        a = p->mp * psn * in->sfcpres * cea / (cs * in->ei) + p->bp;
        b = (p->mp * psn * in->sfcpres / cs + p->bp) * rlb - 1.0;
        c = -rlb;
        if (b >= 0.0)
            q = -0.5 * (b + sqrt(b*b - 4.0*a*c));
        else
            q = -0.5 * (b - sqrt(b*b - 4.0*a*c));
        r1 = q / a;
        r2 = c / q;
        rs = fmax(r1, r2);
        ci = fmax(cs - psn * in->sfcpres * 1.65 * rs, 0.0);
    }
    out->rs = rs * cf;
    out->psn = psn;
    return (isfinite(out->rs) && isfinite(out->psn)) ? 0 : -1;
}
