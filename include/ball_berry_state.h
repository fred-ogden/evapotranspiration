#ifndef BALL_BERRY_STATE_H
#define BALL_BERRY_STATE_H

/* Independent, externally supplied canopy state for PET method 6.
 * Values are set by the driver/BMI, never inferred from each other.
 * LAI and SAI: m2 plant area per m2 ground; SHDFAC: fraction [0,1].
 */
typedef struct {
    double lai;
    double sai;
    double shdfac;
} BallBerryCanopyState;

typedef enum {
    BB_STATE_LAI = 0,
    BB_STATE_SAI = 1,
    BB_STATE_SHDFAC = 2
} BallBerryStateVariable;

/* Returns 0 on success, -1 for invalid pointer, value, or variable.
 * On failure, state is unchanged.
 */
int ball_berry_state_set(BallBerryCanopyState *state,
                         BallBerryStateVariable variable, double value);
int ball_berry_state_validate(const BallBerryCanopyState *state);

/* FLO 10/21: BMI inputs that PET method 6 receives from the modules it is coupled to.
 * Inside the Noah-MP science code these carry the Noah-MP names given in brackets.
 * Initial values come from the configuration file at BMI Initialize; after that the
 * current BMI values are used at every update. */
typedef struct {
    double root_stress_factor_0_to_1;        /* [BTRAN]  from the hydrologic model, previous step; 1 = unstressed */
    double ground_skin_temperature_K;        /* [TG]     from soil freeze-thaw, previous step; NaN until set, then
                                                         TG = air temperature + ground_temperature_offset_from_air_k */
    double surface_resistance_s_per_m;       /* [RSURF]  from bare soil evaporation, previous step; with snow the
                                                         provider supplies the snow-weighted value (Noah-MP OPT_RSF = 4) */
    double snow_cover_fraction;              /* [FSNO]   from the snow module */
    double snow_depth_m;                     /* [SNOWH]  from the snow module */
    double snow_albedo_direct_visible;       /* [ALBSND(1)] from the snow module */
    double snow_albedo_direct_nir;           /* [ALBSND(2)] */
    double snow_albedo_diffuse_visible;      /* [ALBSNI(1)] */
    double snow_albedo_diffuse_nir;          /* [ALBSNI(2)] */
    double canopy_wetted_fraction;           /* [FWET]   from the canopy interception module (CIESM) */
    double canopy_liquid_water_mm;           /* [CANLIQ] from CIESM */
    double canopy_ice_mm;                    /* [CANICE] from CIESM */
} BallBerryCouplingInputs;

/* FLO 10/21: calendar time of the current step, set by the BMI driver before each update.
 * UTC, at the START of the time step (the forcing time stamp). Doubles holding whole
 * numbers (the setter rejects anything else; time_second may be fractional). NaN until set.
 * Method 6 uses them only for the monthly LAI/SAI table (PHENOLOGY, DVEG = 4). */
typedef struct {
    double time_year;            /* e.g. 2016 */
    double time_month;           /* 1 to 12 */
    double time_day_of_month;    /* 1 to 31 (checked against the month at update) */
    double time_hour;            /* 0 to 23 */
    double time_minute;          /* 0 to 59 */
    double time_second;          /* 0 to less than 60 */
} BallBerryTimeInputs;

/* FLO 10/21: method 6 parameters and constants that are not BMI variables. Read once from the
 * configuration file in BMI Initialize; update() uses only this struct and the BMI variables.
 * Noah-MP names in brackets. */
typedef struct {
    int    use_monthly_lai_sai_table;                  /* 1: LAI, SAI from the monthly table every step */
    double latitude_degrees;                             /* only the sign is used (hemisphere) */
    double leaf_area_index_by_month[12];               /* [LAIM] January..December */
    double stem_area_index_by_month[12];               /* [SAIM] */
    double canopy_top_height_m;                        /* [HVT] */
    double canopy_bottom_height_m;                     /* [HVB] */
    double momentum_roughness_length_m;                /* [Z0MVT] */
    double forcing_height_above_canopy_top_m;          /* [ZREF] */
    double characteristic_leaf_dimension_m;            /* [DLEAF] */
    double leaf_orientation_index;                     /* [XL] */
    double cosine_solar_zenith_angle;                  /* [COSZ] fixed for now */
    double direct_shortwave_fraction;                  /* 0.7 in Noah-MP ATM */
    double visible_shortwave_fraction;                 /* 0.5 in Noah-MP ATM */
    double soil_albedo_direct[2];                      /* [ALBSOD] snow-free, 0 = visible, 1 = near infrared */
    double soil_albedo_diffuse[2];                     /* [ALBSOI] */
    double leaf_reflectance[2];                        /* [RHOL] */
    double stem_reflectance[2];                        /* [RHOS] */
    double leaf_transmittance[2];                      /* [TAUL] */
    double stem_transmittance[2];                      /* [TAUS] */
    double soil_emissivity;                            /* [EG(1)] */
    double foliage_nitrogen_percent;                   /* [FOLN] */
    double minimum_temperature_for_photosynthesis_k;   /* [TMIN] sets IGS */
    double co2_mole_fraction;                          /* [CO2] */
    double o2_mole_fraction;                           /* [O2] */
    double ground_temperature_offset_from_air_k;       /* TG = air temperature + offset until TG is set */
    int    initialize_canopy_state_from_forcing;       /* 1: TV = TAH = SFCTMP, EAH = EAIR */
    double initial_vegetation_temperature_k;           /* [TV]  used when the flag above is 0 */
    double initial_canopy_air_temperature_k;           /* [TAH] */
    double initial_canopy_air_vapor_pressure_pa;       /* [EAH] */
} BallBerryMethod6Parameters;

/* FLO 10/21: BMI outputs of PET method 6 for the modules it is coupled to. NaN until
 * the first update. Fluxes follow Noah-MP: TR and EVC are per unit grid-cell area
 * (they include the FVEG weight); SAG is the single ground-absorbed shortwave that
 * Noah-MP applies to both the under-canopy and the bare ground (OPT_RAD = 3);
 * the absorbed longwave, TAH, EAH and RAHG describe the under-canopy (vegetated tile)
 * ground. Positive transpiration and evaporation are water leaving the surface. */
typedef struct {
    double canopy_transpiration_flux_m_per_s;               /* [ETRAN = TR/(LATHEAV*rho_water)] at the supplied BTRAN */
    double canopy_evaporation_heat_flux_W_per_m2;           /* [EVC]  evaporation/sublimation of intercepted water */
    double ground_shortwave_absorbed_below_canopy_W_per_m2; /* [SAG] */
    double ground_longwave_absorbed_below_canopy_W_per_m2;  /* EMG*((1-EMV)*LWDN + EMV*SB*TV**4), from IRG in VEGE_FLUX */
    double canopy_air_temperature_K;                        /* [TAH] */
    double canopy_air_vapor_pressure_Pa;                    /* [EAH] */
    double below_canopy_aerodynamic_resistance_s_per_m;     /* [RAHG] = RAWG, ground to canopy air */
    double vegetation_temperature_K;                        /* [TV] */
} BallBerryCouplingOutputs;
#endif
