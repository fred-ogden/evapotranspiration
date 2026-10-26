# Basic Model Interface (BMI) for potential evapotranspiration (PET) functions:
1. Radiation balance (energy balance) method
2. Aerodynamic method
3. Combination method
4. Priestley Taylor method
5. Penman Monteith method
6. Ball-Berry canopy transpiration, Noah-MP physics (BMI only; see [doc/ball_berry_method6.md](doc/ball_berry_method6.md))

Methods 1-5 return potential evapotranspiration. Method 6 returns the canopy transpiration from
the Noah-MP leaf energy balance with Ball-Berry stomatal resistance, given the root stress factor
of the hydrologic model, plus the coupling outputs a soil freeze-thaw, snow, canopy interception
and bare soil evaporation module need.

# The configuration files (./configs/method_N_*.cfg)
There is one configuration file per method, all with the same layout:

| File | Method |
|---|---|
| `configs/method_1_radiation_balance.cfg` | 1 |
| `configs/method_2_aerodynamic_method.cfg` | 2 |
| `configs/method_3_combination_method.cfg` | 3 |
| `configs/method_4_priestley_taylor_method.cfg` | 4 |
| `configs/method_5_penman_monteith_method.cfg` | 5 |
| `configs/method_6_ball_berry.cfg` | 6 |

Each file has four commented sections: 1 control, 2 time and 3 site, with the same keys in all
six files, then 4 the parameters of the method. The format is `key = value`, one per line, and
`#` starts a comment. All six set `forcing_file = BMI` (forcing passed with `set_value`) and name a
`forcing_csv_file` that the test drivers read. For methods 1-5, `forcing_file` may instead name a
CSV file that PET reads itself. Each instance of a PET model should have its own configuration file.

The variables an instance exposes on its BMI depend on its method: methods 1-5 have the seven
AORC forcing inputs and one output, `water_potential_evaporation_flux`; method 6 has 32 inputs
(forcing, vegetation state, coupling inputs, calendar time) and 10 outputs.

# Compiling and testing
With CMake (shared library `petbmi` and all tests):

```
cmake -S . -B build && cmake --build build && (cd build && ctest --output-on-failure)
```

The test scripts, run from the repository root:

| Script | Test |
|---|---|
| `test/make_and_run_bmi_unit_test.sh` | BMI unit test, all six methods ([test/README.md](test/README.md)) |
| `./test_bmi.sh` | BMI regression and negative tests (method 5) |
| `./make_and_run_read_forcings.sh` | standalone, PET reads its own forcing (methods 1-5) |
| `./make_and_run_pass_forcings.sh` | standalone, the AORC BMI module passes the forcing (all six methods); run `git submodule update --init` first |
| `./noahmp_comparison/check_C_against_stored_reference.sh` | method 6 against the stored NWM Noah-MP results |
| `./noahmp_comparison/run_ball_berry_C_and_Fortran.sh` | method 6 against NWM Noah-MP built from its source (not in this repository) |

Each script prints its usage with `-h`. The read-forcing script prints the methods 1-5 unit test
results next to their expected values. An environment file for Cheyenne is included (`source env_cheyenne.sh`).

To build this code for use in the [Next Generation Water Resources Modeling Framework](https://github.com/NOAA-OWP/ngen), please follow the build instructions in [INSTALL.md](INSTALL.md).

# This rough code outline shows a basic outline of workflow. 
The `bmi_pet.c` file runs BMI functions that initialize, update and finalize an instance of a PET model. It also includes descriptive functions to interpret specifics of the model, such as variable names, units, time/timestep, etc. And it also allows a user (or framework) to get and set values in this model. The `bmi_pet.c` code interacts with the `pet.c` code, which sets up the model based on the PET method chosen. For instance, the aerodynamic method does not calculate the net radiation before calling the PET subroutine, but the other PET methods do. This pet.c file then calls one of the six PET methods. When a method is called it returns a value for PET in m/s, and that is set directly to the BMI model structure.
![code_flow](./figs/bmi_pet.png)

# Notes from author of the PET functions
Version 1.0 by Fred L. Ogden, NOAA-NWS-OWP, May, 2020 includes five different methods to calculate ET, from Chow, Maidment & Mays Textbook, and UNFAO Penman-Monteith:  
1. energy balance method
2. aerodynamic method
3. combination method, which combines 1 & 2.
4. Priestley-Taylor method, which assumes the ratio between 1 & 2, and only calculates 1.
5. Penman-Monteith method, which requires a value of canopy resistance term, and does not rely on 1 or 2.
This subroutine requires a considerable amount of meteorological data as input.
   * a) temperature and (relative-humidity or specific humidity) and the heights at which they are measured.
   * b) near surface wind speed measurement and the height at which it was measured.
   * c) the ambient atmospheric temperature lapse rate
   * d) the fraction of the sky covered by clouds
   * e) (optional) the height above ground to the cloud base. If not provided, then assumed.
   * f) the day of the year (1-366) and time of day (UTC only!)
   * g) the skin temperature of the earth's surface.
   * h) the zero-plane roughness height of the atmospheric boundary layer assuming log-law behavior (from land cover)
   * i) the average root zone soil temperature, or near-surface water temperature in the case of lake evaporation.
   * j) the incoming solar (shortwave) radiation.  If not provided it is computed from d,e,f, using an updated method similar to the one presented in Bras, R.L. Hydrology.  Requires value of the Linke atmospheric turbidity factor, which varies from 2 for clear mountain air to 5 for smoggy air.  According to Hove & Manyumbu 2012, who calculated values over Zimbabwe that varied from 2.14 to 3.71.  Other values exist in the literature.  

# Parameter estimation
Values for roughness length and zero plane displacement height can be estimated from land cover data. For example, the momentum roughness length table below was taken from an [archived website](https://web.archive.org/web/20200103171837/https://websites.pmc.ucsc.edu/~jnoble/wind/extrap/) based on the original European Wind Atlas. Users may access other empirical tables via an updated [European Wind Atlas](https://gmd.copernicus.org/articles/13/5053/2020/) or [other sources](https://collaboration.cmc.ec.gc.ca/science/rpn/gem/gem-climate/Version_3.3.0/dictionary_geophys.pdf). 

| Roughness class  | Momentum Roughness length (m)  | Landscape Type                                                                             |  
| ---------- | ----------- | ------------------------------------------------------------------------------------------- |  
| 0         | 0.0002     | Smooth water surface                                                                       |  
| 0.2       | 0.0005     | Inlet water                                                                                |  
| 0.5       | 0.0024     | Completely open terrain, smooth surface, e.g. concrete runways in airports, mowed grass, etc.|  
| 1         | 0.03       | Open agricultural area without fences and hedgerows and very scattered buildings. Only softly rounded hills |  
| 1.5       | 0.055      | Agricultural land with some houses and 8 metre tall sheltering hedgerows with a distance of approximately 1250 metres |  
| 2         | 0.1        | Agricultural land with some houses and 8 metre tall sheltering hedgerows with a distance of approximately 500 metres |  
| 2.5       | 0.2        | Agricultural land with many houses, shrubs and plants, or 8 metre tall sheltering hedgerows with a distance of approximately 250 metres |  
| 3         | 0.4        | Villages, small towns, agricultural land with many or tall sheltering hedgerows, forests and very rough and uneven terrain |  
| 3.5       | 0.8        | Larger cities with tall buildings|  
| 4         | 1.6        | Very large cities with tall buildings and skyscrapers|  


You may also estimate various parameters using the UN FAO Penman-Monteith example [here](http://www.fao.org/3/X0490E/x0490e06.htm#aerodynamic%20resistance%20). 

The above source indicates the zero plane roughness length, `d` can be approximated as 2/3 of the vegetation height (`H`): `d=2/3*H` for grassland/cropland. For taller vegetation (e.g., coniferous forests), a `2/3*H` estimate may produce a `d` that is greater than the wind speed measurement height, which will cause the PET calculations to fail (you cannot take the logarithm of a negative number). In this case, the model automatically limits `d` to be 2/3 of the wind speed measurement height if it detects that the land-cover-based estimate of the former is greater than the latter.  
The momentum roughness length `zom` can be estimated as `0.123*H` or `0.1845*d` if not taken from the above table or similar empirical sources.  
The heat transfer roughness length `zoh` can be approximated as `0.1*zom`.  

# A note on code adaptation for BMI
This code was minimally changed from the author's original version. These minor changes were made by Nextgen NWM formulation team:
* Much of this C code was moved to `*.h` files, with the intention of being more easily integrated into the Nextgen Framework. It turned out that this step was not strictly necessary, and that this standalone module could use the standard `*.c` files. See known issues below for a discussion on turning these back to `*.c` files.
* At one point the C code was slightly modified to compile as C++ code. This was intended for easier integration with the Nextgen Framework. This was reversed when developed as a standalone module.
* Functions that are necessary for PET calculations, but are not part of an individual method (e.g., `calculate_solar_radiation`), were moved to this file: `./include/pet_tools.h`
* The PET code was split up so that the functions were in standalone files, and they could be called independently. An example: `EtPenmanMonteithMethod.h`
* Variables within the PET functions (e.g., `psychrometric_constant_Pa_per_C`) are passed as part of the model structure (e.g., `model->inter_vars.psychrometric_constant_Pa_per_C`)

# Code issues as of May 5th 2021
* Much of this code was moved from a single original `*.c` file to several `*.h` files, for easier integration into the Nextgen Framework. This is likely not the best practice, and should probably be reversed, so that functional code is stored in `*.c` files and declarations are stored in `*.h` files.
