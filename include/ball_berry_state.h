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
#endif
