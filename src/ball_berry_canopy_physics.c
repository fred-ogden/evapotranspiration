/*
 * ball_berry_canopy_physics.c - the physics of one PET method 6 (Ball-Berry canopy energy balance)
 * time step on resolved inputs, in the form Tapenade differentiates (see ball_berry_canopy_physics.h).
 *
 * Author: Fred L. Ogden, Ph.D., P.E.
 *         NOAA/National Weather Service
 *
 * This software was developed by a U.S. Government employee as part of
 * official duties and is not subject to U.S. copyright protection under
 * 17 U.S.C. Section 105.
 *
 * Licensed under the Apache License, Version 2.0.
 * SPDX-License-Identifier: Apache-2.0
 */

/*
 * FLO 10/28: moved here from ball_berry_helpers.c (the functions of the VEGE_FLUX loop1 replica:
 * ALBEDO/TWOSTREAM/SURRAD, ESAT, SFCDIF1, RAGRB, the linearized foliage step and the exact solver),
 * from ball_berry_method6_inputs.c (the science of bb6_build_canopy_input()) and from pet.c (the
 * coupling outputs of ball_berry_method6_core_step()).  The arithmetic is unchanged.  The changes are
 * the Tapenade-ready form:
 *   - one return per function (a status variable; the checks that returned early now skip the rest);
 *   - bb_is_finite() for isfinite(); no printf (the exact solver records its trace rows);
 *   - no whole-structure assignment or memset (bb_copy_<type>, bb_zero_<type>, generated);
 *   - no ternaries, no break (the loop1 exit is a flag), no NAN;
 *   - fmin(fmax()) written as two statements;
 *   - copysign(1.0e-6, sigma) in TWOSTREAM written with a comparison (the same except at sigma = -0.0);
 *   - every chain of two or more multiplications and divisions (and every division by a product, and every
 *     product with a quotient in parentheses) written one operation per statement, in C's left-to-right
 *     order: the tangent and reverse code then recompute the primal bitwise (Tapenade regroups chains).
 * Byte-identical to the code before 10/28.  ASCII only.
 */
#include <math.h>
#include <stddef.h>
#include "ball_berry_canopy_physics.h"
#include "ball_berry_struct_copies.h"

/* Constants from module_sf_noahmplsm.F (SFCDIF1). */
#define VKC 0.40
#define GRAV 9.80616
#define CPAIR 1004.64

/* Noah-MP constants of the canopy input (as in ball_berry_method6_inputs.c before 10/28) */
#define BB6_SPECIFIC_HEAT_AIR_J_PER_KG_K 1004.64        /* CPAIR */
#define BB6_GAS_CONSTANT_DRY_AIR_J_PER_KG_K 287.04      /* RAIR */
#define BB6_LATENT_HEAT_VAPORIZATION_J_PER_KG 2.5104e6  /* HVAP */
#define BB6_LATENT_HEAT_SUBLIMATION_J_PER_KG 2.8440e6   /* HSUB */
#define BB6_STEFAN_BOLTZMANN_W_PER_M2_K4 5.67e-8        /* SB */
#define BB6_FREEZING_POINT_K 273.16                     /* TFRZ */
#define BB6_MINIMUM_POSITIVE_VALUE 1.0e-6               /* MPE */
#define BB6_BARE_SOIL_ROUGHNESS_M 0.002                 /* Z0 in ENERGY */
#define BB6_SNOW_SURFACE_ROUGHNESS_M 0.002              /* Z0SNO */
#define BB6_MINIMUM_WIND_SPEED_M_PER_S 1.0              /* UR >= 1 in ENERGY */
#define BB6_MAXIMUM_EFFECTIVE_AREA_INDEX 6.0            /* VAIE = MIN(6, VAI) */
#define BB6_SOLAR_CONSTANT_W_PER_M2 1360.0

int bb_is_finite(double x)
{
    int finite;
    finite = 0;
    if (x - x == 0.0) finite = 1;   /* false for an infinity and for NaN */
    return finite;
}

/* ===================================================================== SURRAD (OPT_RAD = 3) */
int noahmp_surrad(const NoahmpSurradInput *in, NoahmpSurradOutput *out)
{
    double cad[2], cai[2], trd, tri, absorbed, laifra, product;
    int ib, status;
    status = 0;
    if (in == 0 || out == 0) status = 1;
    if (status == 0) {
        if (!bb_is_finite(in->mpe) || in->mpe <= 0.0) status = 1;
    }
    if (status == 0) {
        bb_zero_NoahmpSurradOutput(out);
        for (ib = 0; ib < 2; ++ib) {
            cad[ib] = in->solad[ib] * in->fabd[ib];
            cai[ib] = in->solai[ib] * in->fabi[ib];
            out->sav = out->sav + cad[ib] + cai[ib];
            out->fsa = out->fsa + cad[ib] + cai[ib];
            trd = in->solad[ib] * in->ftdd[ib];
            tri = in->solad[ib] * in->ftid[ib] + in->solai[ib] * in->ftii[ib];
            absorbed = trd * (1.0 - in->albgrd[ib]) + tri * (1.0 - in->albgri[ib]);
            out->sag = out->sag + absorbed;
            out->fsa = out->fsa + absorbed;
        }
        laifra = in->elai / fmax(in->vai, in->mpe);
        if (in->fsun > 0.0) {
            product = (cad[0] + in->fsun * cai[0]) * laifra;
            out->parsun = product / fmax(in->laisun, in->mpe);
            product = in->fsha * cai[0];
            product = product * laifra;
            out->parsha = product / fmax(in->laisha, in->mpe);
        } else {
            out->parsun = 0.0;
            product = (cad[0] + cai[0]) * laifra;
            out->parsha = product / fmax(in->laisha, in->mpe);
        }
        out->fsr = in->albd[0] * in->solad[0] + in->albi[0] * in->solai[0]
                 + in->albd[1] * in->solad[1] + in->albi[1] * in->solai[1];
        for (ib = 0; ib < 2; ++ib) {
            out->fsrv = out->fsrv + (in->frevd[ib] * in->solad[ib] + in->frevi[ib] * in->solai[ib]);
            out->fsrg = out->fsrg + (in->fregd[ib] * in->solad[ib] + in->fregi[ib] * in->solai[ib]);
        }
    }
    return status;
}

/* ===================================================================== TWOSTREAM (OPT_RAD = 3)
 * Algebra and evaluation sequence follow TWOSTREAM in module_sf_noahmplsm(3).F, lines 3095-3341.
 * Noah-MP REAL constants are represented here in double precision.  Writes the entries of band in->ib
 * of out (and gdir); the other band's entries are left as they are. */
int noahmp_twostream(const NoahmpTwostreamInput *in, NoahmpTwostreamOutput *out)
{
    int ib, status;
    double gap, kopen, coszi, chil, phi1, phi2, gdir, ext, avmu;
    double omegal, asu, betadl, betail, omega, betad, betai;
    double tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9;
    double b, c, d, f, h, sigma, p1, p2, p3, p4, s1, s2, u1, u2, u3;
    double d1, d2, h1, h2, h3, h4, h5, h6, h7, h8, h9, h10;
    double ftds, ftis, fres, freveg, frebar, alb, fab;
    double product, quotient, denominator, first_term, second_term, parenthesis;
    status = 0;
    ib = 0;
    if (in == NULL || out == NULL) status = 1;
    if (status == 0) {
        if (in->ib < 0 || in->ib > 1 ||
            (in->ic != 0 && in->ic != 1) || !bb_is_finite(in->cosz) ||
            in->cosz <= 0.0 || in->cosz > 1.0 || !bb_is_finite(in->vai) ||
            in->vai <= 0.0 || !bb_is_finite(in->fveg) || in->fveg < 0.0 ||
            in->fveg > 1.0 || !bb_is_finite(in->temperature_k) ||
            !bb_is_finite(in->fwet) || in->fwet < 0.0 || in->fwet > 1.0 ||
            !bb_is_finite(in->xl)) status = 1;
    }
    if (status == 0) {
        ib = in->ib;
        if (!bb_is_finite(in->rho[ib]) || !bb_is_finite(in->tau[ib]) ||
            !bb_is_finite(in->albgrd[ib]) || !bb_is_finite(in->albgri[ib]) ||
            in->albgrd[ib] <= 0.0 || in->albgri[ib] <= 0.0 ||
            in->albgrd[ib] >= 1.0 || in->albgri[ib] >= 1.0) status = 1;
    }
    if (status == 0) {
        gap = 1.0 - in->fveg;
        kopen = 1.0 - in->fveg;
        coszi = fmax(0.001, in->cosz);
        chil = fmax(in->xl, -0.4);
        chil = fmin(chil, 0.6);
        if (fabs(chil) <= 0.01) chil = 0.01;
        product = 0.330 * chil;
        product = product * chil;
        phi1 = 0.5 - 0.633 * chil - product;
        phi2 = 0.877 * (1.0 - 2.0 * phi1);
        gdir = phi1 + phi2 * coszi;
        ext = gdir / coszi;
        quotient = phi1 / phi2;
        quotient = quotient * log((phi1 + phi2) / phi1);
        avmu = (1.0 - quotient) / phi2;
        omegal = in->rho[ib] + in->tau[ib];
        if (omegal <= 0.0 || !bb_is_finite(avmu)) status = 1;
    }
    if (status == 0) {
        tmp0 = gdir + phi2 * coszi;
        tmp1 = phi1 * coszi;
        quotient = tmp1 / tmp0;
        quotient = quotient * log((tmp1 + tmp0) / tmp1);
        parenthesis = 1.0 - quotient;
        product = 0.5 * omegal;
        product = product * gdir;
        product = product / tmp0;
        asu = product * parenthesis;
        denominator = omegal * avmu;
        denominator = denominator * ext;
        quotient = (1.0 + avmu * ext) / denominator;
        betadl = quotient * asu;
        product = 0.5 * (in->rho[ib] + in->tau[ib] +
                  (in->rho[ib] - in->tau[ib]) * pow((1.0 + chil) / 2.0, 2.0));
        betail = product / omegal;
        if (in->temperature_k > 273.16) {
            tmp0 = omegal;
            tmp1 = betadl;
            tmp2 = betail;
        } else {
            tmp0 = (1.0 - in->fwet) * omegal + in->fwet * in->omegas[ib];
            if (tmp0 == 0.0) {
                status = 1;
            } else {
                first_term = (1.0 - in->fwet) * omegal;
                first_term = first_term * betadl;
                second_term = in->fwet * in->omegas[ib];
                second_term = second_term * in->betads;
                tmp1 = (first_term + second_term) / tmp0;
                first_term = (1.0 - in->fwet) * omegal;
                first_term = first_term * betail;
                second_term = in->fwet * in->omegas[ib];
                second_term = second_term * in->betais;
                tmp2 = (first_term + second_term) / tmp0;
            }
        }
    }
    if (status == 0) {
        omega = tmp0; betad = tmp1; betai = tmp2;
        b = 1.0 - omega + omega * betai;
        c = omega * betai;
        tmp0 = avmu * ext;
        d = tmp0 * omega;
        d = d * betad;
        f = tmp0 * omega;
        f = f * (1.0 - betad);
        tmp1 = b * b - c * c;
        if (tmp1 < 0.0) status = 1;
    }
    if (status == 0) {
        h = sqrt(tmp1) / avmu;
        sigma = tmp0 * tmp0 - tmp1;
        if (fabs(sigma) < 1.0e-6) {
            if (sigma < 0.0) sigma = -1.0e-6;
            else sigma = 1.0e-6;
        }
        p1 = b + avmu * h;
        p2 = b - avmu * h;
        p3 = b + tmp0;
        p4 = b - tmp0;
        s1 = exp(-h * in->vai);
        s2 = exp(-ext * in->vai);
        if (s1 == 0.0) status = 1;
    }
    if (status == 0) {
        if (in->ic == 0) alb = in->albgrd[ib];
        else alb = in->albgri[ib];
        u1 = b - c / alb;
        u2 = b - c * alb;
        u3 = f + c * alb;
        tmp2 = u1 - avmu * h;
        tmp3 = u1 + avmu * h;
        first_term = p1 * tmp2;
        first_term = first_term / s1;
        second_term = p2 * tmp3;
        second_term = second_term * s1;
        d1 = first_term - second_term;
        tmp4 = u2 + avmu * h;
        tmp5 = u2 - avmu * h;
        d2 = tmp4 / s1 - tmp5 * s1;
        if (d1 == 0.0 || d2 == 0.0) status = 1;
    }
    if (status == 0) {
        h1 = -d * p4 - c * f;
        quotient = h1 * p3;
        quotient = quotient / sigma;
        tmp6 = d - quotient;
        quotient = h1 / sigma;
        quotient = quotient * (u1 + tmp0);
        tmp7 = (d - c - quotient) * s2;
        quotient = tmp6 * tmp2;
        quotient = quotient / s1;
        h2 = (quotient - p2 * tmp7) / d1;
        product = tmp6 * tmp3;
        product = product * s1;
        h3 = -(product - p1 * tmp7) / d1;
        h4 = -f * p3 - c * d;
        tmp8 = h4 / sigma;
        tmp9 = (u3 - tmp8 * (u2 - tmp0)) * s2;
        quotient = tmp8 * tmp4;
        quotient = quotient / s1;
        h5 = -(quotient + tmp9) / d2;
        product = tmp8 * tmp5;
        product = product * s1;
        h6 = (product + tmp9) / d2;
        product = c * tmp2;
        denominator = d1 * s1;
        h7 = product / denominator;
        product = -c * tmp3;
        product = product * s1;
        h8 = product / d1;
        denominator = d2 * s1;
        h9 = tmp4 / denominator;
        product = -tmp5 * s1;
        h10 = product / d2;
        if (in->ic == 0) {
            ftds = s2 * (1.0 - gap) + gap;
            quotient = h4 * s2;
            quotient = quotient / sigma;
            ftis = (quotient + h5 * s1 + h6 / s1) * (1.0 - gap);
            fres = (h1 / sigma + h2 + h3) * (1.0 - gap) + in->albgrd[ib] * gap;
            freveg = (h1 / sigma + h2 + h3) * (1.0 - gap);
            frebar = in->albgrd[ib] * gap;
        } else {
            ftds = 0.0;
            ftis = (h9 * s1 + h10 / s1) * (1.0 - kopen) + kopen;
            fres = (h7 + h8) * (1.0 - kopen) + in->albgri[ib] * kopen;
            freveg = fres;
            frebar = 0.0;
        }
        fab = 1.0 - fres - (1.0 - in->albgrd[ib]) * ftds -
              (1.0 - in->albgri[ib]) * ftis;
        if (!bb_is_finite(fab) || !bb_is_finite(fres) ||
            !bb_is_finite(ftds) || !bb_is_finite(ftis)) status = 1;
    }
    if (status == 0) {
        out->ftd[ib] = ftds;
        out->fti[ib] = ftis;
        out->fre[ib] = fres;
        out->frev[ib] = freveg;
        out->freg[ib] = frebar;
        out->fab[ib] = fab;
        out->gdir = gdir;
    }
    return status;
}

/* ===================================================================== ALBEDO + RADIATION (OPT_RAD = 3) */
int noahmp_albedo_radiation(const NoahmpAlbedoInput *in, NoahmpAlbedoOutput *out)
{
    const double mpe = 1.0e-6;
    NoahmpTwostreamInput tw_in;
    NoahmpTwostreamOutput tw_out;
    NoahmpSurradInput *rad;
    double vai, wl, ws, rho[2], tau[2], ext, gdir, quotient;
    int ib, ic, status;
    status = 0;
    vai = 0.0;
    if (in == NULL || out == NULL) status = 1;
    if (status == 0) {
        if (!bb_is_finite(in->cosz) || !bb_is_finite(in->elai) || !bb_is_finite(in->esai) ||
            !bb_is_finite(in->fveg) || !bb_is_finite(in->fwet) ||
            in->elai < 0.0 || in->esai < 0.0 || in->fveg < 0.0 ||
            in->fveg > 1.0 || in->fwet < 0.0 || in->fwet > 1.0) status = 1;
    }
    if (status == 0) {
        vai = in->elai + in->esai;
        /* The TWOSTREAM port requires positive VAI (no exposed canopy is handled by the caller). */
        if (vai <= 0.0 || in->cosz > 1.0) status = 1;
    }
    if (status == 0) {
        bb_zero_NoahmpAlbedoOutput(out);
        bb_zero_NoahmpTwostreamInput(&tw_in);
        bb_zero_NoahmpTwostreamOutput(&tw_out);
        rad = &out->radiation_inputs;
        /* FLO 10/26: sun down: ALBEDO returns its zero initialization (IF(COSZ <= 0) GOTO 100), FSUN = 0,
         * so LAISUN = 0, LAISHA = ELAI, and SURRAD absorbs nothing */
        if (in->cosz <= 0.0) {
            for (ib = 0; ib < 2; ib++) {
                rad->solad[ib] = in->solad[ib];
                rad->solai[ib] = in->solai[ib];
            }
            out->fsun = 0.0;
            rad->fsun = 0.0;
            rad->fsha = 1.0;
            rad->elai = in->elai;
            rad->vai = vai;
            rad->laisun = 0.0;
            rad->laisha = in->elai;
            rad->mpe = mpe;
            if (noahmp_surrad(rad, &out->radiation) != 0) status = 1;
        } else {
            wl = in->elai / fmax(vai, mpe);
            ws = in->esai / fmax(vai, mpe);
            tw_in.cosz = in->cosz;
            tw_in.vai = vai;
            tw_in.fveg = in->fveg;
            tw_in.fwet = in->fwet;
            tw_in.temperature_k = in->tv_k;
            tw_in.xl = in->xl;
            tw_in.betads = in->betads;
            tw_in.betais = in->betais;
            for (ib = 0; ib < 2; ib++) {
                rho[ib] = fmax(in->rhol[ib] * wl + in->rhos[ib] * ws, mpe);
                tau[ib] = fmax(in->taul[ib] * wl + in->taus[ib] * ws, mpe);
                tw_in.rho[ib] = rho[ib];
                tw_in.tau[ib] = tau[ib];
                tw_in.albgrd[ib] = in->albgrd[ib];
                tw_in.albgri[ib] = in->albgri[ib];
                tw_in.omegas[ib] = in->omegas[ib];
                rad->albgrd[ib] = in->albgrd[ib];
                rad->albgri[ib] = in->albgri[ib];
                rad->solad[ib] = in->solad[ib];
                rad->solai[ib] = in->solai[ib];
                if (!bb_is_finite(rho[ib]) || !bb_is_finite(tau[ib]) ||
                    !bb_is_finite(in->solad[ib]) || !bb_is_finite(in->solai[ib])) status = 1;
            }
            for (ib = 0; ib < 2; ib++) {
                for (ic = 0; ic < 2; ic++) {
                    if (status == 0) {
                        tw_in.ib = ib;
                        tw_in.ic = ic;
                        if (noahmp_twostream(&tw_in, &tw_out) != 0) status = 1;
                    }
                    if (status == 0) {
                        if (ic == 0) {
                            rad->fabd[ib] = tw_out.fab[ib];
                            rad->albd[ib] = tw_out.fre[ib];
                            rad->ftdd[ib] = tw_out.ftd[ib];
                            rad->ftid[ib] = tw_out.fti[ib];
                            rad->frevd[ib] = tw_out.frev[ib];
                            rad->fregd[ib] = tw_out.freg[ib];
                        } else {
                            rad->fabi[ib] = tw_out.fab[ib];
                            rad->albi[ib] = tw_out.fre[ib];
                            rad->ftii[ib] = tw_out.fti[ib];
                            rad->frevi[ib] = tw_out.frev[ib];
                            rad->fregi[ib] = tw_out.freg[ib];
                        }
                    }
                }
            }
            if (status == 0) {
                gdir = tw_out.gdir;
                quotient = gdir / in->cosz;
                ext = quotient * sqrt(1.0 - rho[0] - tau[0]);
                if (!bb_is_finite(ext)) status = 1;
            }
            if (status == 0) {
                out->fsun = (1.0 - exp(-ext * vai)) / fmax(ext * vai, mpe);
                if (out->fsun < 0.01) out->fsun = 0.0;
                rad->fsun = out->fsun;
                rad->fsha = 1.0 - out->fsun;
                rad->elai = in->elai;
                rad->vai = vai;
                rad->laisun = in->elai * out->fsun;
                rad->laisha = in->elai * (1.0 - out->fsun);
                rad->mpe = mpe;
                if (noahmp_surrad(rad, &out->radiation) != 0) status = 1;
            }
        }
    }
    return status;
}

/* ===================================================================== ESAT
 * ESAT polynomial from module_sf_noahmplsm.F, with water/ice switch at 0 C.
 * FLO 10/20: argument formed exactly as VEGE_FLUX does before calling ESAT:
 * TDC(T) = MIN(50, MAX(-50, T - TFRZ)) with TFRZ = 273.16 K (was T - 273.15,
 * unclamped; that gave a 0.066 percent high EI at 285 K). */
int noahmp_esat(double temperature_k, double *esat_pa, double *slope_pa_per_k)
{
    const double noahmp_freezing_point_k = 273.16;   /* TFRZ */
    double t, es, slope;
    int status;
    status = 0;
    if (esat_pa == 0 || slope_pa_per_k == 0) status = -1;
    if (status == 0) {
        if (!bb_is_finite(temperature_k)) status = -1;
    }
    if (status == 0) {
        t = fmax(-50.0, temperature_k - noahmp_freezing_point_k);
        t = fmin(50.0, t);
        if (t > 0.0) {
            es = 100.0 * (6.107799961 + t * (4.436518521e-1 + t *
                 (1.428945805e-2 + t * (2.650648471e-4 + t *
                 (3.031240396e-6 + t * (2.034080948e-8 + t * 6.136820929e-11))))));
            slope = 100.0 * (4.438099984e-1 + t * (2.857002636e-2 + t *
                 (7.938054040e-4 + t * (1.215215065e-5 + t *
                 (1.036561403e-7 + t * (3.532421810e-10 - t * 7.090244804e-13))))));
        } else {
            es = 100.0 * (6.109177956 + t * (5.034698970e-1 + t *
                 (1.886013408e-2 + t * (4.176223716e-4 + t *
                 (5.824720280e-6 + t * (4.838803174e-8 + t * 1.838826904e-10))))));
            slope = 100.0 * (5.030305237e-1 + t * (3.773255020e-2 + t *
                 (1.267995369e-3 + t * (2.477563108e-5 + t *
                 (3.005693132e-7 + t * (2.158542548e-9 + t * 7.131097725e-12))))));
        }
        if (!bb_is_finite(es) || !bb_is_finite(slope) || es <= 0.0 || slope <= 0.0) status = -1;
    }
    if (status == 0) {
        *esat_pa = es;
        *slope_pa_per_k = slope;
    }
    return status;
}

/* ===================================================================== RAGRB
 * Port of VEGE_FLUX UC assignment and RAGRB, OPT_SFC=1 path.
 * Noah-MP constants VKC = 0.40, GRAV = 9.80616, CPAIR = 1004.64 (module_sf_noahmplsm.F).
 * FLO 10/21: GRAV and CPAIR were 9.81 and 1004.6, which made RAHG differ from Noah-MP by
 * about 6e-6 relative from the second canopy iteration on (through the under-canopy
 * Monin-Obukhov length); with the Noah-MP values the exact solver matches the r8 reference.
 * Invalid or singular domains are rejected instead of generating NaNs. */
int noahmp_ragrb(const NoahmpRagrbInputs *in,
                 NoahmpRagrbState *state, NoahmpRagrbOutputs *out)
{
    const double vkc = 0.4;
    const double grav = 9.80616;     /* GRAV  (m/s2) */
    const double cpair = 1004.64;    /* CPAIR (J/kg/K) */
    double uc, tmp1, molg, mozg, fhgnew, fhg, cwpc, tmp2;
    double tmprah2, kh, rb, rahg, log_denom;
    double product, quotient, denominator;
    int status;
    status = 0;
    if (in == NULL || state == NULL || out == NULL) status = -1;
    if (status == 0) {
        if (in->iteration < 1 || !bb_is_finite(in->reference_wind_speed_m_per_s) ||
            !bb_is_finite(in->reference_height_m) || !bb_is_finite(in->momentum_roughness_m) ||
            !bb_is_finite(in->displacement_height_m) || !bb_is_finite(in->canopy_height_m) ||
            !bb_is_finite(in->vegetation_area_index) || !bb_is_finite(in->air_density_kg_per_m3) ||
            !bb_is_finite(in->ground_sensible_heat_w_per_m2) ||
            !bb_is_finite(in->canopy_air_temperature_k) ||
            !bb_is_finite(in->ground_momentum_roughness_m) ||
            !bb_is_finite(in->ground_heat_roughness_m) ||
            !bb_is_finite(in->canopy_heat_roughness_m) ||
            !bb_is_finite(in->friction_velocity_m_per_s) || !bb_is_finite(in->cwp) ||
            !bb_is_finite(in->leaf_dimension_m) || !bb_is_finite(in->mpe) ||
            in->reference_wind_speed_m_per_s <= 0.0 ||
            in->reference_height_m <= in->momentum_roughness_m ||
            in->momentum_roughness_m <= 0.0 || in->canopy_height_m <= in->displacement_height_m ||
            in->ground_momentum_roughness_m <= 0.0 || in->ground_heat_roughness_m < 0.0 ||
            in->canopy_heat_roughness_m < 0.0 || in->vegetation_area_index <= 0.0 ||
            in->air_density_kg_per_m3 <= 0.0 || in->canopy_air_temperature_k <= 0.0 ||
            in->friction_velocity_m_per_s <= 0.0 || in->cwp <= 0.0 ||
            in->leaf_dimension_m <= 0.0 || in->mpe <= 0.0) status = -1;
    }
    if (status == 0) {
        log_denom = log(in->reference_height_m / in->momentum_roughness_m);
        if (log_denom == 0.0) status = -1;
    }
    if (status == 0) {
        /* The second (overwriting) UC expression in VEGE_FLUX is authoritative. */
        product = in->reference_wind_speed_m_per_s *
                  log((in->canopy_height_m - in->displacement_height_m +
                       in->momentum_roughness_m) / in->momentum_roughness_m);
        uc = product / log_denom;
        if (!bb_is_finite(uc) || uc <= 0.0) status = -1;
    }
    if (status == 0) {
        mozg = 0.0;
        if (in->iteration > 1) {
            quotient = grav / in->canopy_air_temperature_k;
            product = vkc * quotient;
            product = product * in->ground_sensible_heat_w_per_m2;
            denominator = in->air_density_kg_per_m3 * cpair;
            tmp1 = product / denominator;
            if (fabs(tmp1) <= in->mpe) tmp1 = in->mpe;
            molg = -pow(in->friction_velocity_m_per_s, 3.0) / tmp1;
            mozg = fmin((in->displacement_height_m -
                         in->ground_momentum_roughness_m) / molg, 1.0);
        }
        if (mozg < 0.0) fhgnew = pow(1.0 - 15.0 * mozg, -0.25);
        else fhgnew = 1.0 + 4.7 * mozg;
        if (in->iteration == 1) fhg = fhgnew;
        else fhg = 0.5 * (state->fhg + fhgnew);
        product = in->cwp * in->vegetation_area_index;
        product = product * in->canopy_height_m;
        product = product * fhg;
        cwpc = sqrt(product);
        if (!bb_is_finite(cwpc) || cwpc <= 0.0) status = -1;
    }
    if (status == 0) {
        product = -cwpc * in->ground_heat_roughness_m;
        product = product / in->canopy_height_m;
        tmp1 = exp(product);
        product = -cwpc * (in->canopy_heat_roughness_m + in->displacement_height_m);
        product = product / in->canopy_height_m;
        tmp2 = exp(product);
        product = in->canopy_height_m * exp(cwpc);
        product = product / cwpc;
        tmprah2 = product * (tmp1 - tmp2);
        product = vkc * in->friction_velocity_m_per_s;
        product = product * (in->canopy_height_m - in->displacement_height_m);
        kh = fmax(product, in->mpe);
        rahg = tmprah2 / kh;
        product = cwpc * 50.0;
        rb = product / (1.0 - exp(-cwpc / 2.0));
        rb = fmax(rb * sqrt(in->leaf_dimension_m / uc), 20.0);
        if (!bb_is_finite(rahg) || !bb_is_finite(rb)) status = -1;
    }
    if (status == 0) {
        state->mozg = mozg;
        state->fhg = fhg;
        out->canopy_wind_speed_m_per_s = uc;
        out->ground_momentum_resistance_s_per_m = 0.0;
        out->ground_heat_resistance_s_per_m = rahg;
        out->ground_vapor_resistance_s_per_m = rahg;
        out->leaf_boundary_resistance_s_per_m = rb;
    }
    return status;
}

/* ===================================================================== SFCDIF1 */
static void stability_functions(double moz, double *fm, double *fh)
{
    double x, a, b;
    if (moz < 0.0) {
        x = pow(1.0 - 16.0 * moz, 0.25);
        a = log((1.0 + x * x) / 2.0);
        b = log((1.0 + x) / 2.0);
        *fm = 2.0 * b + a - 2.0 * atan(x) + 1.5707963;
        *fh = 2.0 * a;
    } else {
        *fm = -5.0 * moz;
        *fh = *fm;
    }
}

int noahmp_sfcdif1(const NoahmpSfcdif1Inputs *in,
                   NoahmpSfcdif1State *state, NoahmpSfcdif1Outputs *out)
{
    double old_moz, moz2, mol, tmp1, tvir;
    double tmpcm, tmpch, tmpcm2, tmpch2;
    double fmnew, fhnew, fm2new, fh2new;
    double cmfm, chfh, cm2fm2, ch2fh2;
    double product, quotient, denominator, numerator;
    int status;
    status = 0;
    moz2 = 0.0;
    if (in == NULL || state == NULL || out == NULL) status = -1;
    if (status == 0) {
        if (in->iteration < 1 || in->reference_height_m <= in->displacement_height_m ||
            in->momentum_roughness_m <= 0.0 || in->heat_roughness_m <= 0.0 ||
            in->air_density_kg_per_m3 <= 0.0 || in->air_temperature_k <= 0.0 ||
            in->wind_speed_m_per_s < 0.0 || in->mpe <= 0.0) status = -1;
    }
    if (status == 0) {
        old_moz = state->moz;
        tmpcm = log((in->reference_height_m-in->displacement_height_m)/in->momentum_roughness_m);
        tmpch = log((in->reference_height_m-in->displacement_height_m)/in->heat_roughness_m);
        tmpcm2 = log((2.0+in->momentum_roughness_m)/in->momentum_roughness_m);
        tmpch2 = log((2.0+in->heat_roughness_m)/in->heat_roughness_m);
        if (in->iteration == 1) {
            state->friction_velocity_m_per_s = 0.0;
            state->moz = 0.0;
            moz2 = 0.0;
        } else {
            tvir = (1.0+0.61*in->specific_humidity_kg_per_kg)*in->air_temperature_k;
            quotient = GRAV/tvir;
            product = VKC*quotient;
            product = product*in->sensible_heat_w_per_m2;
            denominator = in->air_density_kg_per_m3*CPAIR;
            tmp1 = product/denominator;
            if (fabs(tmp1) <= in->mpe) tmp1 = in->mpe;
            mol = -pow(state->friction_velocity_m_per_s,3.0)/tmp1;
            if (mol == 0.0) {
                status = -1;
            } else {
                state->moz = fmin((in->reference_height_m-in->displacement_height_m)/mol,1.0);
                moz2 = fmin((2.0+in->heat_roughness_m)/mol,1.0);
            }
        }
    }
    if (status == 0) {
        if (old_moz*state->moz < 0.0) state->moz_sign_changes++;
        if (state->moz_sign_changes >= 2) {
            state->moz=0.0; state->fm=0.0; state->fh=0.0;
            moz2=0.0; state->fm2=0.0; state->fh2=0.0;
        }
        stability_functions(state->moz,&fmnew,&fhnew);
        stability_functions(moz2,&fm2new,&fh2new);
        if (in->iteration == 1) {
            state->fm=fmnew; state->fh=fhnew;
            state->fm2=fm2new; state->fh2=fh2new;
        } else {
            state->fm=0.5*(state->fm+fmnew);
            state->fh=0.5*(state->fh+fhnew);
            state->fm2=0.5*(state->fm2+fm2new);
            state->fh2=0.5*(state->fh2+fh2new);
        }
        state->fh=fmin(state->fh,0.9*tmpch);
        state->fm=fmin(state->fm,0.9*tmpcm);
        state->fh2=fmin(state->fh2,0.9*tmpch2);
        state->fm2=fmin(state->fm2,0.9*tmpcm2);
        cmfm=tmpcm-state->fm; chfh=tmpch-state->fh;
        cm2fm2=tmpcm2-state->fm2; ch2fh2=tmpch2-state->fh2;
        if (fabs(cmfm)<=in->mpe) cmfm=in->mpe;
        if (fabs(chfh)<=in->mpe) chfh=in->mpe;
        if (fabs(cm2fm2)<=in->mpe) cm2fm2=in->mpe;
        if (fabs(ch2fh2)<=in->mpe) ch2fh2=in->mpe;
        numerator=VKC*VKC;
        denominator=cmfm*cmfm;
        out->cm=numerator/denominator;
        denominator=cmfm*chfh;
        out->ch=numerator/denominator;
        denominator=cm2fm2*ch2fh2;
        out->ch2=numerator/denominator;
        state->friction_velocity_m_per_s=in->wind_speed_m_per_s*sqrt(out->cm);
        out->friction_velocity_m_per_s=state->friction_velocity_m_per_s;
        product=VKC*state->friction_velocity_m_per_s;
        out->ch2=product/ch2fh2;
        if (!bb_is_finite(out->cm) || !bb_is_finite(out->ch) || !bb_is_finite(out->ch2)) status = -1;
    }
    return status;
}

/* ===================================================================== VEGE_FLUX foliage step */
int noahmp_vege_flux_foliage_step(const NoahmpVegeFluxInputs *in,
                                 NoahmpVegeFluxOutputs *out)
{
    double cah, cvh, cgh, cond, ata, bta, csh;
    double caw, cew, ctw, cgw, aea, bea, cev, ctr;
    double air_term, cir, residual, jacobian, dtv, tv, estv, tah, eah;
    double irc, shc, evc, tr, latent_limit;
    double transpiration_ramp_factor, heat_storage, fusion_heat, liquid_end_mm, freezing_rate;
    double tanh_argument, tanh_value, ramp_center_k, ramp_scale_k;
    double transpiration_w_per_m2;
    double product, second_product;
    const double freezing_point_k = 273.16;            /* TFRZ */
    const double latent_heat_fusion_j_per_kg = 0.3336e6;   /* HFUS */
    int status;
    status = 0;
    tah = 0.0;      /* (set before use; initialized for the compiler's flow analysis) */
    eah = 0.0;
    if (in == 0 || out == 0) status = -1;
    if (status == 0) {
        if (!(in->fveg > 0.0 && in->fveg <= 1.0 && in->vai >= 0.0 &&
              in->lai_sun >= 0.0 && in->lai_shade >= 0.0 &&
              in->wet_fraction >= 0.0 && in->wet_fraction <= 1.0 &&
              in->rb_s_per_m > 0.0 && in->rs_sun_s_per_m >= 0.0 &&
              in->rs_shade_s_per_m >= 0.0 && in->rahc_s_per_m > 0.0 &&
              in->rahg_s_per_m > 0.0 && in->rawc_s_per_m > 0.0 &&
              in->rawg_s_per_m > 0.0 && in->rsurf_s_per_m >= 0.0 &&
              in->air_density_kg_per_m3 > 0.0 && in->cp_air_j_per_kg_k > 0.0 &&
              in->psychrometric_constant_pa_per_k > 0.0 &&
              in->timestep_s > 0.0 && in->latent_heat_j_per_kg > 0.0 &&
              in->vegetation_emissivity >= 0.0 && in->vegetation_emissivity <= 1.0 &&
              in->ground_emissivity >= 0.0 && in->ground_emissivity <= 1.0 &&
              in->stefan_boltzmann_w_per_m2_k4 > 0.0)) status = -1;
    }
    if (status == 0) {
        if (!(in->canopy_heat_capacity_j_per_m2_k >= 0.0) || !(in->phase_change_window_k >= 0.0) ||
            !(in->frozen_transpiration_ramp_k >= 0.0) || !(in->canopy_water_total_mm >= 0.0) ||
            !(in->canopy_liquid_start_mm >= 0.0)) status = -1;
    }
    if (status == 0) {
        tv = in->foliage_temperature_k;
        estv = in->foliage_saturation_vapor_pressure_pa;
        cah = 1.0 / in->rahc_s_per_m;
        product = 2.0 * fmin(6.0, in->vai);
        cvh = product / in->rb_s_per_m;
        cgh = 1.0 / in->rahg_s_per_m;
        cond = cah + cvh + cgh;
        ata = (in->air_temperature_k * cah + in->ground_temperature_k * cgh) / cond;
        bta = cvh / cond;
        product = (1.0 - bta) * in->air_density_kg_per_m3;
        product = product * in->cp_air_j_per_kg_k;
        csh = product * cvh;
        caw = 1.0 / in->rawc_s_per_m;
        product = in->wet_fraction * fmin(6.0, in->vai);
        cew = product / in->rb_s_per_m;
        ctw = (1.0 - in->wet_fraction) *
              (fmin(6.0, in->lai_sun) / (in->rb_s_per_m + in->rs_sun_s_per_m) +
               fmin(6.0, in->lai_shade) / (in->rb_s_per_m + in->rs_shade_s_per_m));
        /* FLO 10/24: frozen canopy, stomata closed: CTW ramps smoothly to 0 below TFRZ (logistic centered
         * half way down the span, 98 percent open at TFRZ, 2 percent at TFRZ - span). */
        transpiration_ramp_factor = 1.0;
        if (in->frozen_transpiration_ramp_k > 0.0) {
            ramp_center_k = freezing_point_k - 0.5 * in->frozen_transpiration_ramp_k;
            ramp_scale_k = in->frozen_transpiration_ramp_k / 8.0;
            transpiration_ramp_factor = 1.0 / (1.0 + exp(-(tv - ramp_center_k) / ramp_scale_k));
            ctw = ctw * transpiration_ramp_factor;
        }
        cgw = 1.0 / (in->rawg_s_per_m + in->rsurf_s_per_m);
        cond = caw + cew + ctw + cgw;
        aea = (in->air_vapor_pressure_pa * caw +
               in->ground_saturation_vapor_pressure_pa * cgw) / cond;
        bea = (cew + ctw) / cond;
        product = (1.0 - bea) * cew;
        product = product * in->air_density_kg_per_m3;
        product = product * in->cp_air_j_per_kg_k;
        cev = product / in->psychrometric_constant_pa_per_k;
        product = (1.0 - bea) * ctw;
        product = product * in->air_density_kg_per_m3;
        product = product * in->cp_air_j_per_kg_k;
        ctr = product / in->psychrometric_constant_pa_per_k;
        tah = ata + bta * tv;
        eah = aea + bea * estv;
        product = -in->vegetation_emissivity *
                  (1.0 + (1.0 - in->vegetation_emissivity) * (1.0 - in->ground_emissivity));
        product = product * in->incoming_longwave_w_per_m2;
        second_product = in->vegetation_emissivity * in->ground_emissivity;
        second_product = second_product * in->stefan_boltzmann_w_per_m2_k4;
        second_product = second_product * pow(in->ground_temperature_k, 4.0);
        air_term = product - second_product;
        product = (2.0 - in->vegetation_emissivity * (1.0 - in->ground_emissivity)) *
                  in->vegetation_emissivity;
        cir = product * in->stefan_boltzmann_w_per_m2_k4;
        irc = in->fveg * (air_term + cir * pow(tv, 4.0));
        product = in->fveg * in->air_density_kg_per_m3;
        product = product * in->cp_air_j_per_kg_k;
        product = product * cvh;
        shc = product * (tv - tah);
        product = in->fveg * in->air_density_kg_per_m3;
        product = product * in->cp_air_j_per_kg_k;
        product = product * cew;
        product = product * (estv - eah);
        evc = product / in->psychrometric_constant_pa_per_k;
        product = in->fveg * in->air_density_kg_per_m3;
        product = product * in->cp_air_j_per_kg_k;
        product = product * ctw;
        product = product * (estv - eah);
        tr = product / in->psychrometric_constant_pa_per_k;
        if (in->evaporation_phase_from_latent_heat) {
            /* FLO 10/24: the phase that evaporates is the one the latent heat was chosen for (HSUB: ice),
             * fixed for the step, so the limit does not switch as the iterate TV crosses TFRZ */
            if (in->latent_heat_j_per_kg > 0.5 * (2.5104e6 + 2.8440e6)) {
                product = in->canopy_ice_mm * in->latent_heat_j_per_kg;
            } else {
                product = in->canopy_liquid_mm * in->latent_heat_j_per_kg;
            }
        } else if (tv > 273.16) {
            product = in->canopy_liquid_mm * in->latent_heat_j_per_kg;
        } else {
            product = in->canopy_ice_mm * in->latent_heat_j_per_kg;
        }
        latent_limit = product / in->timestep_s;
        /* EVC = MIN(latent limit, EVC), written as a comparison: fmin() of +0 and -0 may return either zero,
         * and the compiler may swap its arguments (it did when this function changed, 10/28), which changed
         * the sign of a zero EVC; the comparison keeps EVC in a tie, as the code before 10/28 did */
        if (latent_limit < evc) evc = latent_limit;
        residual = in->absorbed_shortwave_w_per_m2 - irc - shc - evc - tr + in->additional_heat_w_per_m2;
        product = 4.0 * cir;
        product = product * pow(tv, 3.0);
        jacobian = in->fveg * (product + csh +
                   (cev + ctr) * in->foliage_saturation_slope_pa_per_k);
        /* FLO 10/24: canopy heat storage, C_v*(TV - TV0)/DT, removed from the net heating */
        heat_storage = 0.0;
        if (in->canopy_heat_capacity_j_per_m2_k > 0.0) {
            product = in->canopy_heat_capacity_j_per_m2_k * (tv - in->previous_foliage_temperature_k);
            heat_storage = product / in->timestep_s;
            residual = residual - heat_storage;
            jacobian = jacobian + in->canopy_heat_capacity_j_per_m2_k / in->timestep_s;
        }
        /* FLO 10/24: phase change of intercepted water (enthalpy form). At the end of the step the liquid
         * is (CANLIQ0+CANICE0)*0.5*(1+tanh((TV-TFRZ)/w)); the fusion heat HFUS*(CANLIQ0 - CANLIQ_end)/DT
         * (+ released by freezing) enters the canopy energy balance. It is self-limiting: freezing warms the
         * canopy toward TFRZ, melting cools it. Smooth in TV, and exactly energy conserving. */
        fusion_heat = 0.0;
        freezing_rate = 0.0;
        liquid_end_mm = in->canopy_liquid_start_mm;
        if (in->phase_change_window_k > 0.0 && in->canopy_water_total_mm > 0.0) {
            tanh_argument = (tv - freezing_point_k) / in->phase_change_window_k;
            tanh_value = tanh(tanh_argument);
            product = in->canopy_water_total_mm * 0.5;
            liquid_end_mm = product * (1.0 + tanh_value);
            freezing_rate = (in->canopy_liquid_start_mm - liquid_end_mm) / in->timestep_s;
            fusion_heat = latent_heat_fusion_j_per_kg * freezing_rate;      /* mm/s = kg m-2 s-1 */
            residual = residual + fusion_heat;
            product = latent_heat_fusion_j_per_kg * in->canopy_water_total_mm;
            product = product * 0.5;
            product = product * (1.0 - tanh_value * tanh_value);
            product = product / in->phase_change_window_k;
            product = product / in->timestep_s;
            jacobian = jacobian + product;
        }
        if (!bb_is_finite(jacobian) || fabs(jacobian) < 1e-15) status = -1;
    }
    if (status == 0) {
        dtv = residual / jacobian;
        product = in->fveg * ctr;
        product = product * in->foliage_saturation_slope_pa_per_k;
        product = product * dtv;
        transpiration_w_per_m2 = tr + product;
        if (!bb_is_finite(transpiration_w_per_m2) || !bb_is_finite(dtv)) status = -1;
    }
    if (status == 0) {
        out->energy_balance_residual_w_per_m2 = residual;
        out->unlinearized_transpiration_w_per_m2 = tr;
        product = in->fveg * 4.0;
        product = product * cir;
        product = product * pow(tv, 3.0);
        product = product * dtv;
        out->net_longwave_w_per_m2 = irc + product;
        product = in->fveg * csh;
        product = product * dtv;
        out->sensible_heat_w_per_m2 = shc + product;
        product = in->fveg * cev;
        product = product * in->foliage_saturation_slope_pa_per_k;
        product = product * dtv;
        out->interception_evaporation_w_per_m2 = evc + product;
        out->transpiration_w_per_m2 = transpiration_w_per_m2;
        out->canopy_air_temperature_k = tah;
        out->canopy_air_vapor_pressure_pa = eah;
        out->foliage_temperature_increment_k = dtv;
        out->updated_foliage_temperature_k = tv + dtv;
        out->transpiration_conductance_m_per_s = ctw;
        out->heat_storage_w_per_m2 = heat_storage;
        out->fusion_heat_w_per_m2 = fusion_heat;
        out->canopy_liquid_end_mm = liquid_end_mm;
        out->freezing_rate_mm_per_s = freezing_rate;
        out->transpiration_ramp_factor = transpiration_ramp_factor;
    }
    return status;
}

/* ===================================================================== VEGE_FLUX loop1 replica
 * FLO 10/20: Noah-MP exact iteration mode (see the description in ball_berry_helpers.h). Variable names
 * in comments are the Noah-MP names in module_sf_noahmplsm.F.  FLO 10/28: the trace rows are recorded in
 * trace (noahmp_vege_flux_exact_solve(), ball_berry_helpers.c, writes them). */
int noahmp_vege_flux_exact_solve_traced(const NoahmpStabilityInput *input, NoahmpStabilityOutput *output,
                                        NoahmpVegeFluxExactTrace *trace)
{
    const int maximum_canopy_iterations = NOAHMP_EXACT_MAXIMUM_CANOPY_ITERATIONS;   /* NITERC */
    const int minimum_iterations_before_exit = 5;          /* ITER >= 5 */
    const double exit_temperature_change_k = 0.01;         /* ABS(DTV) <= 0.01 */
    NoahmpSfcdif1Inputs surface_layer_inputs;
    NoahmpSfcdif1State surface_layer_state;
    NoahmpRagrbInputs ground_canopy_inputs;
    NoahmpRagrbState ground_canopy_state;
    NoahmpAlbedoInput radiation_inputs;
    NoahmpAlbedoOutput radiation_outputs;
    NoahmpStomataInputs leaf_inputs;
    NoahmpStomataOutputs sunlit_stomata, shaded_stomata;
    NoahmpVegeFluxInputs foliage_inputs;
    NoahmpVegeFluxOutputs foliage_step;
    NoahmpVegeFluxExactTraceRecord *row;
    double vegetation_temperature_k;                   /* TV  */
    double canopy_air_temperature_k;                   /* TAH */
    double canopy_air_vapor_pressure_pa;               /* EAH */
    double ground_saturation_vapor_pressure_pa;        /* ESTG */
    double ground_saturation_slope_pa_per_k;
    double leaf_saturation_vapor_pressure_pa;          /* ESTV */
    double leaf_saturation_slope_pa_per_k;             /* DESTV */
    double canopy_to_atmosphere_sensible_heat_w_per_m2;  /* H  */
    double ground_to_canopy_sensible_heat_w_per_m2;      /* HG */
    double above_canopy_heat_resistance_s_per_m;         /* RAHC */
    double air_density_kg_per_m3, specific_heat_air_j_per_kg_k, wind_speed_m_per_s;
    double photosynthesis_sunlit, photosynthesis_shaded, product, denominator;
    int iteration, last_iteration_flag, loop_done, status;

    status = 0;
    canopy_to_atmosphere_sensible_heat_w_per_m2 = 0.0;
    ground_to_canopy_sensible_heat_w_per_m2 = 0.0;
    above_canopy_heat_resistance_s_per_m = 0.0;
    last_iteration_flag = 0;
    iteration = 1;
    if (input == NULL || output == NULL || trace == NULL) status = -1;
    if (status == 0) {
        trace->started = 0;
        trace->count = 0;
        bb_zero_NoahmpStabilityOutput(output);
        bb_zero_NoahmpSfcdif1State(&surface_layer_state);
        bb_zero_NoahmpRagrbState(&ground_canopy_state);
        bb_zero_NoahmpStomataOutputs(&sunlit_stomata);
        bb_zero_NoahmpStomataOutputs(&shaded_stomata);
        bb_copy_NoahmpSfcdif1Inputs(&input->surface_layer, &surface_layer_inputs);
        bb_copy_NoahmpRagrbInputs(&input->ground_canopy, &ground_canopy_inputs);
        bb_copy_NoahmpVegeFluxInputs(&input->canopy.foliage, &foliage_inputs);
        bb_copy_NoahmpAlbedoInput(&input->canopy.radiation_stomata.radiation, &radiation_inputs);
        bb_copy_NoahmpStomataInputs(&input->canopy.radiation_stomata.canopy_air, &leaf_inputs);

        air_density_kg_per_m3 = foliage_inputs.air_density_kg_per_m3;
        specific_heat_air_j_per_kg_k = foliage_inputs.cp_air_j_per_kg_k;
        wind_speed_m_per_s = surface_layer_inputs.wind_speed_m_per_s;          /* UR */
        if (!(wind_speed_m_per_s > 0.0) || !(air_density_kg_per_m3 > 0.0) ||
            !(specific_heat_air_j_per_kg_k > 0.0)) status = -1;
    }
    if (status == 0) {
        vegetation_temperature_k = foliage_inputs.foliage_temperature_k;
        canopy_air_temperature_k = ground_canopy_inputs.canopy_air_temperature_k;
        if (!(canopy_air_temperature_k > 0.0))
            canopy_air_temperature_k = foliage_inputs.air_temperature_k;
        canopy_air_vapor_pressure_pa = leaf_inputs.ea;

        /* Radiation: RADIATION is called by ENERGY before VEGE_FLUX, once, at the incoming TV. */
        radiation_inputs.tv_k = vegetation_temperature_k;
        if (noahmp_albedo_radiation(&radiation_inputs, &radiation_outputs) != 0) status = -1;
    }
    if (status == 0) {
        foliage_inputs.lai_sun = radiation_outputs.radiation_inputs.laisun;
        foliage_inputs.lai_shade = radiation_outputs.radiation_inputs.laisha;
        foliage_inputs.absorbed_shortwave_w_per_m2 = radiation_outputs.radiation.sav;

        /* ESTG at the (prescribed) ground temperature, before loop1. */
        if (noahmp_esat(foliage_inputs.ground_temperature_k, &ground_saturation_vapor_pressure_pa,
                        &ground_saturation_slope_pa_per_k) != 0) status = -1;
    }
    if (status == 0) {
        foliage_inputs.ground_saturation_vapor_pressure_pa = ground_saturation_vapor_pressure_pa;
        trace->started = 1;
    }

    loop_done = 0;
    if (status != 0) loop_done = 1;
    while (loop_done == 0 && iteration <= maximum_canopy_iterations) {
        /* SFCDIF1: aerodynamic exchange between ZLVL and ZPD+Z0M, using H from the last iteration. */
        surface_layer_inputs.iteration = iteration;
        surface_layer_inputs.sensible_heat_w_per_m2 = canopy_to_atmosphere_sensible_heat_w_per_m2;
        if (noahmp_sfcdif1(&surface_layer_inputs, &surface_layer_state,
                           &output->surface_layer) != 0) status = -1;
        if (status == 0) {
            denominator = output->surface_layer.ch * wind_speed_m_per_s;
            above_canopy_heat_resistance_s_per_m = fmax(1.0, 1.0 / denominator);   /* RAHC */

            /* RAGRB: under-canopy resistance and leaf boundary resistance, using HG and TAH from the last iteration. */
            ground_canopy_inputs.iteration = iteration;
            ground_canopy_inputs.friction_velocity_m_per_s = output->surface_layer.friction_velocity_m_per_s;
            ground_canopy_inputs.ground_sensible_heat_w_per_m2 = ground_to_canopy_sensible_heat_w_per_m2;
            ground_canopy_inputs.canopy_air_temperature_k = canopy_air_temperature_k;
            if (noahmp_ragrb(&ground_canopy_inputs, &ground_canopy_state,
                             &output->ground_canopy) != 0) status = -1;
        }
        if (status == 0) {
            /* ESTV and DESTV at the current TV. */
            if (noahmp_esat(vegetation_temperature_k, &leaf_saturation_vapor_pressure_pa,
                            &leaf_saturation_slope_pa_per_k) != 0) status = -1;
        }
        /* STOMATA only on the first iteration, with the incoming TV, ESTV, EAH and this RB. */
        if (status == 0 && iteration == 1) {
            leaf_inputs.tv = vegetation_temperature_k;
            leaf_inputs.ei = leaf_saturation_vapor_pressure_pa;
            leaf_inputs.ea = canopy_air_vapor_pressure_pa;
            leaf_inputs.rb = output->ground_canopy.leaf_boundary_resistance_s_per_m;
            leaf_inputs.apar = radiation_outputs.radiation.parsun;
            if (noahmp_stomata(&input->canopy.radiation_stomata.physiology,
                               &leaf_inputs, &sunlit_stomata) != 0) status = -1;
            if (status == 0) {
                leaf_inputs.apar = radiation_outputs.radiation.parsha;
                if (noahmp_stomata(&input->canopy.radiation_stomata.physiology,
                                   &leaf_inputs, &shaded_stomata) != 0) status = -1;
            }
            if (status == 0) {
                foliage_inputs.rs_sun_s_per_m = sunlit_stomata.rs;
                foliage_inputs.rs_shade_s_per_m = shaded_stomata.rs;
            }
        }
        if (status == 0) {
            /* One linearized foliage energy balance step (TAH, EAH, fluxes, DTV). */
            foliage_inputs.foliage_temperature_k = vegetation_temperature_k;
            foliage_inputs.foliage_saturation_vapor_pressure_pa = leaf_saturation_vapor_pressure_pa;
            foliage_inputs.foliage_saturation_slope_pa_per_k = leaf_saturation_slope_pa_per_k;
            foliage_inputs.rb_s_per_m = output->ground_canopy.leaf_boundary_resistance_s_per_m;
            foliage_inputs.rahc_s_per_m = above_canopy_heat_resistance_s_per_m;
            foliage_inputs.rawc_s_per_m = above_canopy_heat_resistance_s_per_m;        /* RAWC = RAHC */
            foliage_inputs.rahg_s_per_m = output->ground_canopy.ground_heat_resistance_s_per_m;
            foliage_inputs.rawg_s_per_m = output->ground_canopy.ground_vapor_resistance_s_per_m;
            if (noahmp_vege_flux_foliage_step(&foliage_inputs, &foliage_step) != 0) status = -1;
        }
        if (status == 0) {
            vegetation_temperature_k = foliage_step.updated_foliage_temperature_k;     /* TV = TV + DTV */
            canopy_air_temperature_k = foliage_step.canopy_air_temperature_k;          /* TAH (pre-update TV) */
            canopy_air_vapor_pressure_pa = foliage_step.canopy_air_vapor_pressure_pa;  /* EAH */

            /* H and HG for the next iteration's Monin-Obukhov lengths. */
            product = air_density_kg_per_m3 * specific_heat_air_j_per_kg_k;
            product = product * (canopy_air_temperature_k - foliage_inputs.air_temperature_k);
            canopy_to_atmosphere_sensible_heat_w_per_m2 = product / above_canopy_heat_resistance_s_per_m;
            product = air_density_kg_per_m3 * specific_heat_air_j_per_kg_k;
            product = product * (foliage_inputs.ground_temperature_k - canopy_air_temperature_k);
            ground_to_canopy_sensible_heat_w_per_m2 = product / output->ground_canopy.ground_heat_resistance_s_per_m;

            /* the trace row of this iteration */
            row = &trace->record[trace->count];
            row->iteration = iteration;
            row->vegetation_temperature_k = vegetation_temperature_k;
            row->vegetation_temperature_change_k = foliage_step.foliage_temperature_increment_k;
            row->canopy_air_temperature_k = canopy_air_temperature_k;
            row->canopy_air_vapor_pressure_pa = canopy_air_vapor_pressure_pa;
            row->sunlit_stomatal_resistance_s_per_m = foliage_inputs.rs_sun_s_per_m;
            row->shaded_stomatal_resistance_s_per_m = foliage_inputs.rs_shade_s_per_m;
            row->leaf_boundary_layer_resistance_s_per_m = foliage_inputs.rb_s_per_m;
            row->above_canopy_heat_resistance_s_per_m = above_canopy_heat_resistance_s_per_m;
            row->under_canopy_heat_resistance_s_per_m = foliage_inputs.rahg_s_per_m;
            row->momentum_drag_coefficient = output->surface_layer.cm;
            row->heat_exchange_coefficient = output->surface_layer.ch;
            row->friction_velocity_m_per_s = output->surface_layer.friction_velocity_m_per_s;
            row->monin_obukhov_stability_parameter = surface_layer_state.moz;
            row->canopy_net_longwave_w_per_m2 = foliage_step.net_longwave_w_per_m2;
            row->canopy_sensible_heat_w_per_m2 = foliage_step.sensible_heat_w_per_m2;
            row->transpiration_heat_w_per_m2 = foliage_step.transpiration_w_per_m2;
            trace->count = trace->count + 1;

            /* the loop1 exit: one more iteration once ITER >= 5 and ABS(DTV) <= 0.01 */
            if (last_iteration_flag == 1) {
                loop_done = 1;
            } else {
                if (iteration >= minimum_iterations_before_exit &&
                    fabs(foliage_step.foliage_temperature_increment_k) <= exit_temperature_change_k &&
                    last_iteration_flag == 0) last_iteration_flag = 1;
                iteration = iteration + 1;
            }
        } else {
            loop_done = 1;
        }
    }
    if (status == 0) {
        if (iteration > maximum_canopy_iterations) iteration = maximum_canopy_iterations;

        /* Outputs, in the same structure as the nested solver. */
        bb_copy_NoahmpAlbedoOutput(&radiation_outputs, &output->canopy.radiation_stomata.radiation);
        bb_copy_NoahmpStomataOutputs(&sunlit_stomata, &output->canopy.radiation_stomata.sunlit);
        bb_copy_NoahmpStomataOutputs(&shaded_stomata, &output->canopy.radiation_stomata.shaded);
        photosynthesis_sunlit = sunlit_stomata.psn * radiation_outputs.radiation_inputs.laisun;
        photosynthesis_shaded = shaded_stomata.psn * radiation_outputs.radiation_inputs.laisha;
        output->canopy.radiation_stomata.photosynthesis_umol_m2_s = photosynthesis_sunlit + photosynthesis_shaded;
        bb_copy_NoahmpVegeFluxOutputs(&foliage_step, &output->canopy.foliage.foliage);
        output->canopy.foliage.foliage.canopy_air_temperature_k = canopy_air_temperature_k;
        output->canopy.foliage.foliage.canopy_air_vapor_pressure_pa = canopy_air_vapor_pressure_pa;
        output->canopy.foliage.foliage_temperature_k = vegetation_temperature_k;
        output->canopy.foliage.ground_skin_temperature_k = foliage_inputs.ground_temperature_k;
        product = foliage_step.transpiration_w_per_m2 / foliage_inputs.latent_heat_j_per_kg;
        output->canopy.foliage.transpiration_m_per_s = product / 1000.0;
        output->canopy.iterations = iteration;
        output->canopy.frozen_sunlit_resistance_s_per_m = sunlit_stomata.rs;
        output->canopy.frozen_shaded_resistance_s_per_m = shaded_stomata.rs;
        output->canopy.final_temperature_increment_k = foliage_step.foliage_temperature_increment_k;
        output->canopy.energy_balance_residual_w_per_m2 = foliage_step.energy_balance_residual_w_per_m2;
        output->stability_iterations = iteration;
        output->relative_resistance_change = 0.0;
        output->reference_canopy_resistance_s_per_m = above_canopy_heat_resistance_s_per_m;
        if (!bb_is_finite(output->canopy.foliage.transpiration_m_per_s)) status = -1;
    }
    return status;
}

/* ===================================================================== the canopy input (ATM, PHENOLOGY,
 * GROUNDALB, ENERGY): the science of bb6_build_canopy_input() before 10/28 */
int bb6_canopy_input_physics(const BallBerryMethod6Parameters *parameters, const BallBerryStepInputs *step,
                             NoahmpStabilityInput *canopy_input, BallBerryStepResults *results)
{
    const BallBerryCouplingInputs *coupling;
    NoahmpAlbedoInput *radiation;
    NoahmpStomataInputs *canopy_air;
    NoahmpVegeFluxInputs *foliage;
    NoahmpSfcdif1Inputs *surface_layer;
    NoahmpRagrbInputs *ground_canopy;
    double air_temperature_k, surface_pressure_pa, specific_humidity_kg_per_kg, wind_speed_m_per_s;
    double downward_shortwave_w_per_m2, downward_longwave_w_per_m2;
    double air_vapor_pressure_pa, air_density_kg_per_m3;
    double leaf_area_index, stem_area_index, vegetation_area_index, green_vegetation_fraction;
    double zero_plane_displacement_m, reference_height_m;
    double initial_vegetation_temperature_k, initial_canopy_air_temperature_k;
    double initial_canopy_air_vapor_pressure_pa, latent_heat_vegetation_j_per_kg;
    double mixing_ratio_kg_per_kg, qv, ground_shortwave, ground_albedo_direct, ground_albedo_diffuse;
    double solad[2], solai[2], INC, wind_speed_hypot_m_per_s, db_lower;
    double numerator, denominator, product, ratio, snow_term, soil_term;
    /* Noah-MP names used in the science code below */
    double BTRAN, TG, RSURF, FSNO, SNOWH, FWET, CANLIQ, CANICE;
    double ALBSND[2], ALBSNI[2], ALBSOD[2], ALBSOI[2];
    double HVT, HVB, DB, FB, SNOWHC, ELAI, ESAI, ZPDG, Z0MG, EG, COSZ, COSZ_HORIZ;
    int band;

    coupling = &step->coupling;
    bb_zero_NoahmpStabilityInput(canopy_input);

    /* ---- the coupling inputs (long names) under the Noah-MP names */
    BTRAN  = coupling->root_stress_factor_0_to_1;   /* soil moisture factor for stomatal resistance */
    RSURF  = coupling->surface_resistance_s_per_m;  /* ground surface resistance to evaporation (s/m) */
    FSNO   = coupling->snow_cover_fraction;         /* fraction of the ground covered by snow */
    SNOWH  = coupling->snow_depth_m;                /* snow depth (m) */
    ALBSND[0] = coupling->snow_albedo_direct_visible;   /* snow albedo, direct, visible */
    ALBSND[1] = coupling->snow_albedo_direct_nir;       /* snow albedo, direct, near infrared */
    ALBSNI[0] = coupling->snow_albedo_diffuse_visible;  /* snow albedo, diffuse, visible */
    ALBSNI[1] = coupling->snow_albedo_diffuse_nir;      /* snow albedo, diffuse, near infrared */
    FWET   = coupling->canopy_wetted_fraction;      /* wetted or snowed fraction of the canopy */
    CANLIQ = coupling->canopy_liquid_water_mm;      /* intercepted liquid water (mm) */
    CANICE = coupling->canopy_ice_mm;               /* intercepted ice (mm) */

    /* ---- the forcing, prepared as Noah-MP ATM and ENERGY do */
    air_temperature_k = step->air_temperature_K;                       /* SFCTMP */
    surface_pressure_pa = step->surface_pressure_Pa;                   /* SFCPRS */
    specific_humidity_kg_per_kg = step->specific_humidity_kg_per_kg;   /* QAIR */
    downward_shortwave_w_per_m2 = step->shortwave_down_W_per_m2;       /* SOLDN */
    downward_longwave_w_per_m2 = step->longwave_down_W_per_m2;         /* LWDN */
    wind_speed_hypot_m_per_s = hypot(step->wind_x_m_per_s, step->wind_y_m_per_s);
    wind_speed_m_per_s = fmax(wind_speed_hypot_m_per_s, BB6_MINIMUM_WIND_SPEED_M_PER_S);   /* UR */
    product = specific_humidity_kg_per_kg * surface_pressure_pa;
    air_vapor_pressure_pa = product / (0.622 + 0.378 * specific_humidity_kg_per_kg);     /* EAIR */
    denominator = BB6_GAS_CONSTANT_DRY_AIR_J_PER_KG_K * air_temperature_k;
    air_density_kg_per_m3 = (surface_pressure_pa - 0.378 * air_vapor_pressure_pa) / denominator;   /* RHOAIR */

    /* ---- FLO 10/26: snow-free soil albedo, configured or GROUNDALB from the top-disc soil moisture */
    for (band = 0; band < 2; ++band) {
        if (parameters->soil_albedo_option == 1) {
            INC = fmax(0.11 - 0.40 * coupling->soil_top_disc_moisture_content, 0.0);
            ALBSOD[band] = fmin(parameters->soil_albedo_saturated[band] + INC, parameters->soil_albedo_dry[band]);
            ALBSOI[band] = ALBSOD[band];
        } else {
            ALBSOD[band] = parameters->soil_albedo_direct[band];
            ALBSOI[band] = parameters->soil_albedo_diffuse[band];
        }
    }

    /* ---- FLO 10/26: solar geometry (COSZ from the preparation). Option 1: the shortwave is corrected to
     *      the sloped terrain as NOM ATM does (SWDOWN = SOLDN*COSZ/COSZ_HORIZ, at most the solar constant
     *      times COSZ; 0 when the sun is down); with zero slope COSZ = COSZ_HORIZ and SWDOWN = SOLDN. */
    COSZ = step->cosine_solar_zenith_angle;
    if (parameters->solar_geometry_option == 1) {
        COSZ_HORIZ = step->cosine_solar_zenith_angle_horizontal;
        if (COSZ <= 0.0 || COSZ_HORIZ <= 0.0) {
            downward_shortwave_w_per_m2 = 0.0;
        } else if (downward_shortwave_w_per_m2 / COSZ_HORIZ > BB6_SOLAR_CONSTANT_W_PER_M2) {
            downward_shortwave_w_per_m2 = BB6_SOLAR_CONSTANT_W_PER_M2 * COSZ;
        } else {
            ratio = COSZ / COSZ_HORIZ;
            downward_shortwave_w_per_m2 = downward_shortwave_w_per_m2 * ratio;
        }
    }
    results->coupling_outputs.cosine_solar_zenith_angle = COSZ;

    /* ---- leaf and stem area index (table or BMI values, from the preparation): PHENOLOGY small-value
     *      checks, applied to the working copies */
    leaf_area_index = step->leaf_area_index;
    stem_area_index = step->stem_area_index;
    if (stem_area_index < 0.05) stem_area_index = 0.0;
    if (leaf_area_index < 0.05 || stem_area_index == 0.0) leaf_area_index = 0.0;

    /* ---- burial of the canopy by snow, as in PHENOLOGY (ELAI, ESAI = exposed indices) */
    HVT = parameters->canopy_top_height_m;
    HVB = parameters->canopy_bottom_height_m;
    db_lower = fmax(SNOWH - HVB, 0.0);
    DB = fmin(db_lower, HVT - HVB);
    FB = DB / fmax(1.0e-06, HVT - HVB);
    if (HVT > 0.0 && HVT <= 1.0) {                  /* short vegetation (MB: 1.0 and 0.2 in Noah-MP) */
        SNOWHC = HVT * exp(-SNOWH / 0.2);
        if (SNOWHC <= SNOWH) FB = 1.0;
        else FB = SNOWH / SNOWHC;
    }
    ELAI = leaf_area_index * (1.0 - FB);
    ESAI = stem_area_index * (1.0 - FB);
    if (ESAI < 0.05) ESAI = 0.0;
    if (ELAI < 0.05 || ESAI == 0.0) ELAI = 0.0;
    leaf_area_index = ELAI;                         /* from here on, the exposed indices */
    stem_area_index = ESAI;
    vegetation_area_index = leaf_area_index + stem_area_index;
    green_vegetation_fraction = step->green_vegetation_fraction;   /* FVEG */
    /* FLO 10/26: PHENOLOGY DVEG = 1 floor, and FVEG = 0 without an exposed canopy (NOAHMP_SFLX) */
    if (parameters->dynamic_vegetation_option == 1 && green_vegetation_fraction <= 0.05) green_vegetation_fraction = 0.05;
    if (vegetation_area_index == 0.0) green_vegetation_fraction = 0.0;
    results->coupling_outputs.vegetated_fraction = green_vegetation_fraction;

    /* FLO 10/26: surface resistance with snow, OPT_RSF = 4 (EnergyMain) */
    if (parameters->surface_resistance_snow_blend_option == 1) {
        ratio = 1.0 / parameters->surface_resistance_snow_s_per_m;
        snow_term = FSNO * ratio;
        ratio = 1.0 / fmax(RSURF, 0.001);
        soil_term = (1.0 - FSNO) * ratio;
        RSURF = 1.0 / (snow_term + soil_term);
    }
    results->coupling_outputs.surface_resistance_with_snow_s_per_m = RSURF;

    /* FLO 10/26: no exposed canopy (buried by snow, or bare): Noah-MP does not call VEGE_FLUX. The
     * ground absorbs the shortwave without a canopy; TAH and EAH stay as ATM set them; TV is not updated. */
    results->no_canopy_this_step = 0;
    if (!(vegetation_area_index > 0.0)) {
        qv = specific_humidity_kg_per_kg / (1.0 - specific_humidity_kg_per_kg);
        ground_shortwave = 0.0;
        if (COSZ > 0.0) {
            product = downward_shortwave_w_per_m2 * parameters->direct_shortwave_fraction;
            solad[0] = product * parameters->visible_shortwave_fraction;
            product = downward_shortwave_w_per_m2 * parameters->direct_shortwave_fraction;
            solad[1] = product * (1.0 - parameters->visible_shortwave_fraction);
            product = downward_shortwave_w_per_m2 * (1.0 - parameters->direct_shortwave_fraction);
            solai[0] = product * parameters->visible_shortwave_fraction;
            product = downward_shortwave_w_per_m2 * (1.0 - parameters->direct_shortwave_fraction);
            solai[1] = product * (1.0 - parameters->visible_shortwave_fraction);
            for (band = 0; band < 2; ++band) {
                ground_albedo_direct = ALBSOD[band] * (1.0 - FSNO) + ALBSND[band] * FSNO;
                ground_albedo_diffuse = ALBSOI[band] * (1.0 - FSNO) + ALBSNI[band] * FSNO;
                ground_shortwave += solad[band] * (1.0 - ground_albedo_direct) + solai[band] * (1.0 - ground_albedo_diffuse);
            }
        }
        results->no_canopy_this_step = 1;
        results->no_canopy_ground_shortwave_W_per_m2 = ground_shortwave;
        results->no_canopy_ground_emissivity = parameters->soil_emissivity * (1.0 - FSNO) + FSNO;
        results->no_canopy_air_temperature_K = air_temperature_k;
        product = surface_pressure_pa * qv;
        results->no_canopy_air_vapor_pressure_Pa = product / (0.622 + qv);
        if (step->vegetation_temperature_known == 1)
            results->no_canopy_vegetation_temperature_K = coupling->vegetation_temperature_start_of_step_K;
        else
            results->no_canopy_vegetation_temperature_K = air_temperature_k;
    } else {
        /* ---- ENERGY geometry: ZPD = 0.65*HVT, raised to SNOWH by deep snow; ZLVL = MAX(ZPD,HVT) + ZREF,
         *      raised if the snow surface (ZPDG = SNOWH) reaches it; Z0MG blends bare soil and snow. */
        zero_plane_displacement_m = 0.65 * HVT;
        if (SNOWH > zero_plane_displacement_m) zero_plane_displacement_m = SNOWH;
        reference_height_m = fmax(zero_plane_displacement_m, HVT) + parameters->forcing_height_above_canopy_top_m;
        ZPDG = SNOWH;
        if (ZPDG >= reference_height_m) reference_height_m = ZPDG + parameters->forcing_height_above_canopy_top_m;
        Z0MG = BB6_BARE_SOIL_ROUGHNESS_M * (1.0 - FSNO) + FSNO * BB6_SNOW_SURFACE_ROUGHNESS_M;

        /* ---- initial canopy state: cold start (forcing or configured values) ... */
        if (parameters->initialize_canopy_state_from_forcing == 1) {
            initial_vegetation_temperature_k = air_temperature_k;
            initial_canopy_air_temperature_k = air_temperature_k;
            initial_canopy_air_vapor_pressure_pa = air_vapor_pressure_pa;
        } else if (parameters->initialize_canopy_state_from_forcing == 2) {
            /* FLO 10/27: NOM cold start: TV configured (NOM 298 K), TAH = SFCTMP and EAH from the
             * mixing ratio, as NOM sets them before its first step */
            mixing_ratio_kg_per_kg = specific_humidity_kg_per_kg / (1.0 - specific_humidity_kg_per_kg);
            initial_vegetation_temperature_k = parameters->initial_vegetation_temperature_k;
            initial_canopy_air_temperature_k = air_temperature_k;
            product = surface_pressure_pa * mixing_ratio_kg_per_kg;
            initial_canopy_air_vapor_pressure_pa = product / (0.622 + mixing_ratio_kg_per_kg);
        } else {
            initial_vegetation_temperature_k = parameters->initial_vegetation_temperature_k;
            initial_canopy_air_temperature_k = parameters->initial_canopy_air_temperature_k;
            initial_canopy_air_vapor_pressure_pa = parameters->initial_canopy_air_vapor_pressure_pa;
        }
        /* FLO 10/24: ... unless the start-of-step state is known (carried from the previous step or set) */
        if (step->vegetation_temperature_known == 1)
            initial_vegetation_temperature_k = coupling->vegetation_temperature_start_of_step_K;
        if (step->canopy_air_temperature_known == 1)
            initial_canopy_air_temperature_k = coupling->canopy_air_temperature_start_of_step_K;
        if (step->canopy_air_vapor_pressure_known == 1)
            initial_canopy_air_vapor_pressure_pa = coupling->canopy_air_vapor_pressure_start_of_step_Pa;
        if (initial_vegetation_temperature_k > BB6_FREEZING_POINT_K)
            latent_heat_vegetation_j_per_kg = BB6_LATENT_HEAT_VAPORIZATION_J_PER_KG;   /* LATHEAV */
        else
            latent_heat_vegetation_j_per_kg = BB6_LATENT_HEAT_SUBLIMATION_J_PER_KG;

        /* ---- radiation (ALBEDO, TWOSTREAM, SURRAD) */
        radiation = &canopy_input->canopy.radiation_stomata.radiation;
        radiation->cosz = COSZ;
        radiation->elai = leaf_area_index;
        radiation->esai = stem_area_index;
        radiation->fveg = green_vegetation_fraction;
        radiation->fwet = FWET;
        radiation->tv_k = initial_vegetation_temperature_k;
        radiation->xl = parameters->leaf_orientation_index;
        radiation->omegas[0] = 0.8;     /* OMEGAS, MPTABLE rad_parameters (weighted by FWET below TFRZ) */
        radiation->omegas[1] = 0.4;
        radiation->betads = 0.5;        /* BETADS */
        radiation->betais = 0.5;        /* BETAIS */
        for (band = 0; band < 2; ++band) {       /* 0 = visible, 1 = near infrared */
            radiation->rhol[band] = parameters->leaf_reflectance[band];
            radiation->rhos[band] = parameters->stem_reflectance[band];
            radiation->taul[band] = parameters->leaf_transmittance[band];
            radiation->taus[band] = parameters->stem_transmittance[band];
            /* GROUNDALB: ALBGRD = ALBSOD*(1-FSNO) + ALBSND*FSNO, ALBGRI = ALBSOI*(1-FSNO) + ALBSNI*FSNO */
            radiation->albgrd[band] = ALBSOD[band] * (1.0 - FSNO) + ALBSND[band] * FSNO;
            radiation->albgri[band] = ALBSOI[band] * (1.0 - FSNO) + ALBSNI[band] * FSNO;
        }
        product = downward_shortwave_w_per_m2 * parameters->direct_shortwave_fraction;
        radiation->solad[0] = product * parameters->visible_shortwave_fraction;
        product = downward_shortwave_w_per_m2 * parameters->direct_shortwave_fraction;
        radiation->solad[1] = product * (1.0 - parameters->visible_shortwave_fraction);
        product = downward_shortwave_w_per_m2 * (1.0 - parameters->direct_shortwave_fraction);
        radiation->solai[0] = product * parameters->visible_shortwave_fraction;
        product = downward_shortwave_w_per_m2 * (1.0 - parameters->direct_shortwave_fraction);
        radiation->solai[1] = product * (1.0 - parameters->visible_shortwave_fraction);

        /* ---- physiology (STOMATA): the BMI values, including MP and VCMX25 */
        bb_copy_NoahmpStomataParameters(&step->stomata, &canopy_input->canopy.radiation_stomata.physiology);

        canopy_air = &canopy_input->canopy.radiation_stomata.canopy_air;
        canopy_air->mpe = BB6_MINIMUM_POSITIVE_VALUE;
        canopy_air->foln = parameters->foliage_nitrogen_percent;
        if (initial_vegetation_temperature_k > parameters->minimum_temperature_for_photosynthesis_k)
            canopy_air->igs = 1.0;
        else
            canopy_air->igs = 0.0;
        canopy_air->ea = initial_canopy_air_vapor_pressure_pa;   /* EAH seen by STOMATA on iteration 1 */
        canopy_air->sfctmp = air_temperature_k;
        canopy_air->sfcpres = surface_pressure_pa;
        canopy_air->co2 = parameters->co2_mole_fraction * surface_pressure_pa;   /* CO2AIR */
        canopy_air->o2 = parameters->o2_mole_fraction * surface_pressure_pa;     /* O2AIR */
        canopy_air->btran = BTRAN;                                               /* root stress factor */

        /* ---- foliage energy balance */
        foliage = &canopy_input->canopy.foliage;
        foliage->fveg = green_vegetation_fraction;
        foliage->vai = vegetation_area_index;
        foliage->wet_fraction = FWET;
        foliage->rsurf_s_per_m = RSURF;
        foliage->air_temperature_k = air_temperature_k;
        /* TG from soil freeze-thaw once supplied through BMI; before that, air temperature + offset */
        if (step->ground_temperature_known == 0)
            TG = air_temperature_k + parameters->ground_temperature_offset_from_air_k;
        else
            TG = coupling->ground_skin_temperature_K;
        foliage->ground_temperature_k = TG;
        foliage->foliage_temperature_k = initial_vegetation_temperature_k;
        foliage->air_vapor_pressure_pa = air_vapor_pressure_pa;
        foliage->air_density_kg_per_m3 = air_density_kg_per_m3;
        foliage->cp_air_j_per_kg_k = BB6_SPECIFIC_HEAT_AIR_J_PER_KG_K;
        numerator = BB6_SPECIFIC_HEAT_AIR_J_PER_KG_K * surface_pressure_pa;
        denominator = 0.622 * latent_heat_vegetation_j_per_kg;
        foliage->psychrometric_constant_pa_per_k = numerator / denominator;      /* GAMMAV */
        foliage->incoming_longwave_w_per_m2 = downward_longwave_w_per_m2;
        foliage->vegetation_emissivity = 1.0 - exp(-(leaf_area_index + stem_area_index));        /* EMV */
        EG = parameters->soil_emissivity;                                                         /* EG(1) */
        foliage->ground_emissivity = EG * (1.0 - FSNO) + 1.0 * FSNO;                              /* EMG */
        foliage->stefan_boltzmann_w_per_m2_k4 = BB6_STEFAN_BOLTZMANN_W_PER_M2_K4;
        foliage->canopy_liquid_mm = CANLIQ;
        foliage->canopy_ice_mm = CANICE;
        foliage->latent_heat_j_per_kg = latent_heat_vegetation_j_per_kg;
        foliage->timestep_s = step->time_step_s;                                                  /* DT */
        foliage->additional_heat_w_per_m2 = coupling->precipitation_advected_heat_to_canopy_W_per_m2;   /* PAHV */

        /* ---- FLO 10/24: revised canopy energy physics (option 2). With option 1 these stay 0, which is
         *      exactly the Noah-MP calculation. */
        if (parameters->canopy_energy_physics_option == 2) {
            /* C_v = FVEG * c * (ELAI+ESAI) + heat capacity of the intercepted water (CANLIQ, CANICE per ground area) */
            product = green_vegetation_fraction * parameters->canopy_heat_capacity_per_unit_area_index_j_per_m2_k;
            foliage->canopy_heat_capacity_j_per_m2_k =
                product * vegetation_area_index + (4.188e6 * CANLIQ + 2.094e6 * CANICE) / 1000.0;   /* CWAT, CICE */
            foliage->previous_foliage_temperature_k = initial_vegetation_temperature_k;
            foliage->phase_change_window_k = parameters->canopy_phase_change_window_k;
            foliage->canopy_water_total_mm = CANLIQ + CANICE;
            foliage->canopy_liquid_start_mm = CANLIQ;
            foliage->frozen_transpiration_ramp_k = parameters->frozen_transpiration_ramp_k;
            foliage->evaporation_phase_from_latent_heat = 1;
            canopy_input->canopy.hold_radiation_at_start_temperature = 1;
            canopy_input->noahmp_outer_iteration_rule = 1;     /* never fails in stable, calm air (see ball_berry_helpers.h) */
        }

        /* ---- SFCDIF1 surface layer */
        surface_layer = &canopy_input->surface_layer;
        surface_layer->air_temperature_k = air_temperature_k;
        surface_layer->air_density_kg_per_m3 = air_density_kg_per_m3;
        surface_layer->specific_humidity_kg_per_kg = specific_humidity_kg_per_kg;
        surface_layer->reference_height_m = reference_height_m;                    /* ZLVL */
        surface_layer->displacement_height_m = zero_plane_displacement_m;          /* ZPD */
        surface_layer->momentum_roughness_m = parameters->momentum_roughness_length_m;   /* Z0M = Z0MVT */
        surface_layer->heat_roughness_m = parameters->momentum_roughness_length_m;       /* Z0H = Z0M */
        surface_layer->wind_speed_m_per_s = wind_speed_m_per_s;
        surface_layer->mpe = BB6_MINIMUM_POSITIVE_VALUE;

        /* ---- RAGRB under canopy and leaf boundary layer */
        ground_canopy = &canopy_input->ground_canopy;
        ground_canopy->reference_wind_speed_m_per_s = wind_speed_m_per_s;
        ground_canopy->reference_height_m = reference_height_m;
        ground_canopy->momentum_roughness_m = parameters->momentum_roughness_length_m;
        ground_canopy->displacement_height_m = zero_plane_displacement_m;
        ground_canopy->canopy_height_m = HVT;                                         /* HCAN = HVT */
        ground_canopy->vegetation_area_index = fmin(BB6_MAXIMUM_EFFECTIVE_AREA_INDEX, vegetation_area_index);  /* VAIE */
        ground_canopy->air_density_kg_per_m3 = air_density_kg_per_m3;
        ground_canopy->canopy_air_temperature_k = initial_canopy_air_temperature_k;
        ground_canopy->ground_momentum_roughness_m = Z0MG;
        ground_canopy->ground_heat_roughness_m = Z0MG;                                /* Z0HG = Z0MG */
        ground_canopy->canopy_heat_roughness_m = parameters->momentum_roughness_length_m;   /* Z0H = Z0M */
        ground_canopy->leaf_dimension_m = parameters->characteristic_leaf_dimension_m;     /* DLEAF */
        ground_canopy->cwp = step->canopy_wind_parameter;                              /* BMI variable CWP */
        ground_canopy->mpe = BB6_MINIMUM_POSITIVE_VALUE;
    }
    return 0;
}

/* ===================================================================== one step with the exact solver:
 * the method 6 branch of ball_berry_method6_core_step() (pet.c) before 10/28 */
int ball_berry_method6_physics_step(const BallBerryMethod6Parameters *parameters, const BallBerryStepInputs *step,
                                    NoahmpStabilityInput *canopy_input, NoahmpStabilityOutput *solution,
                                    NoahmpVegeFluxExactTrace *trace, BallBerryStepResults *results)
{
    const NoahmpVegeFluxInputs *canopy_inputs;
    const NoahmpVegeFluxOutputs *canopy_fluxes;
    BallBerryCouplingOutputs *coupling_outputs;
    double TV, EMV, EMG, LWDN, stefan_boltzmann_w_per_m2_k4, longwave_part, longwave_emitted;
    int status;

    status = bb6_canopy_input_physics(parameters, step, canopy_input, results);
    coupling_outputs = &results->coupling_outputs;
    trace->started = 0;
    trace->count = 0;
    results->freezing_rate_defined = 1;
    results->transpiration_m_per_s = 0.0;
    results->pet_m_per_s = 0.0;
    if (status == 0 && results->no_canopy_this_step == 1) {
        /* FLO 10/26: no exposed canopy (ELAI + ESAI = 0, buried by snow or bare): Noah-MP skips VEGE_FLUX.
         * No transpiration or canopy evaporation; the ground absorbs SAG without a canopy and sees the
         * downward longwave directly (EMV = 0); TAH, EAH as ATM sets them; TV is not updated. RAHG is not
         * defined: 1 s/m is written (the SFT column uses it only under a canopy). */
        TV = results->no_canopy_vegetation_temperature_K;
        coupling_outputs->canopy_transpiration_flux_m_per_s = 0.0;
        coupling_outputs->canopy_evaporation_heat_flux_W_per_m2 = 0.0;
        coupling_outputs->ground_shortwave_absorbed_below_canopy_W_per_m2 = results->no_canopy_ground_shortwave_W_per_m2;
        coupling_outputs->ground_longwave_absorbed_below_canopy_W_per_m2 =
            results->no_canopy_ground_emissivity * step->longwave_down_W_per_m2;
        coupling_outputs->canopy_air_temperature_K = results->no_canopy_air_temperature_K;
        coupling_outputs->canopy_air_vapor_pressure_Pa = results->no_canopy_air_vapor_pressure_Pa;
        coupling_outputs->below_canopy_aerodynamic_resistance_s_per_m = 1.0;
        coupling_outputs->vegetation_temperature_K = TV;
        if (TV > 273.16) coupling_outputs->canopy_latent_heat_J_per_kg = 2.5104e6;
        else coupling_outputs->canopy_latent_heat_J_per_kg = 2.8440e6;
        coupling_outputs->canopy_freezing_rate_mm_per_s = 0.0;
        if (parameters->canopy_energy_physics_option != 2) results->freezing_rate_defined = 0;
        coupling_outputs->canopy_fusion_heat_flux_W_per_m2 = 0.0;
        coupling_outputs->canopy_heat_storage_flux_W_per_m2 = 0.0;
    } else if (status == 0) {
        status = noahmp_vege_flux_exact_solve_traced(canopy_input, solution, trace);   /* VEGE_FLUX loop1 replica */
        if (status == 0) {
            /* FLO 10/21: coupling outputs. Noah-MP names on the right; see ball_berry_state.h for the area basis of
             * each. The longwave absorbed by the under-canopy ground is the absorbed part of Noah-MP's IRG
             * (VEGE_FLUX loop2), evaluated at the final vegetation temperature:
             *   IRG = EMG*SB*TG**4 - EMG*(1-EMV)*LWDN - EMG*EMV*SB*TV**4   (+ = to atmosphere) */
            canopy_inputs = &canopy_input->canopy.foliage;
            canopy_fluxes = &solution->canopy.foliage.foliage;
            TV   = solution->canopy.foliage.foliage_temperature_k;                    /* vegetation temperature */
            EMV  = canopy_inputs->vegetation_emissivity;                              /* vegetation emissivity */
            EMG  = canopy_inputs->ground_emissivity;                                  /* ground emissivity */
            LWDN = canopy_inputs->incoming_longwave_w_per_m2;                         /* downward longwave */
            stefan_boltzmann_w_per_m2_k4 = canopy_inputs->stefan_boltzmann_w_per_m2_k4;   /* Noah-MP SB */
            coupling_outputs->canopy_transpiration_flux_m_per_s = solution->canopy.foliage.transpiration_m_per_s;  /* ETRAN */
            coupling_outputs->canopy_evaporation_heat_flux_W_per_m2 = canopy_fluxes->interception_evaporation_w_per_m2; /* EVC */
            coupling_outputs->ground_shortwave_absorbed_below_canopy_W_per_m2 =
                solution->canopy.radiation_stomata.radiation.radiation.sag;                                        /* SAG */
            longwave_part = EMG * (1.0 - EMV);
            longwave_part = longwave_part * LWDN;
            longwave_emitted = EMG * EMV;
            longwave_emitted = longwave_emitted * stefan_boltzmann_w_per_m2_k4;
            longwave_emitted = longwave_emitted * pow(TV, 4.0);
            coupling_outputs->ground_longwave_absorbed_below_canopy_W_per_m2 = longwave_part + longwave_emitted;
            coupling_outputs->canopy_air_temperature_K = canopy_fluxes->canopy_air_temperature_k;                 /* TAH */
            coupling_outputs->canopy_air_vapor_pressure_Pa = canopy_fluxes->canopy_air_vapor_pressure_pa;          /* EAH */
            coupling_outputs->below_canopy_aerodynamic_resistance_s_per_m =
                solution->ground_canopy.ground_heat_resistance_s_per_m;                                            /* RAHG */
            coupling_outputs->vegetation_temperature_K = TV;                                                       /* TV */
            coupling_outputs->canopy_latent_heat_J_per_kg = canopy_inputs->latent_heat_j_per_kg;                   /* LATHEAV */
            /* FLO 10/24: phase change and heat storage of the revised canopy energy physics, at the converged TV */
            if (parameters->canopy_energy_physics_option == 2) {
                coupling_outputs->canopy_freezing_rate_mm_per_s = canopy_fluxes->freezing_rate_mm_per_s;
                coupling_outputs->canopy_fusion_heat_flux_W_per_m2 = canopy_fluxes->fusion_heat_w_per_m2;
                coupling_outputs->canopy_heat_storage_flux_W_per_m2 = canopy_fluxes->heat_storage_w_per_m2;
            } else {
                coupling_outputs->canopy_freezing_rate_mm_per_s = 0.0;     /* NaN in the BMI: CIESM uses its own rate */
                results->freezing_rate_defined = 0;
                coupling_outputs->canopy_fusion_heat_flux_W_per_m2 = 0.0;
                coupling_outputs->canopy_heat_storage_flux_W_per_m2 = 0.0;
            }
            results->transpiration_m_per_s = solution->canopy.foliage.transpiration_m_per_s;
            results->pet_m_per_s = results->transpiration_m_per_s;
            if (results->pet_m_per_s < 0.0) results->pet_m_per_s = 0.0;
        }
    }
    return status;
}
