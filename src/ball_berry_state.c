#include "ball_berry_state.h"
#include <math.h>
#include <stddef.h>

int ball_berry_state_validate(const BallBerryCanopyState *state)
{
    if (state == NULL || !isfinite(state->lai) ||
        !isfinite(state->sai) || !isfinite(state->shdfac) ||
        state->lai < 0.0 || state->sai < 0.0 ||
        state->shdfac < 0.0 || state->shdfac > 1.0)
        return -1;
    return 0;
}

int ball_berry_state_set(BallBerryCanopyState *state,
                         BallBerryStateVariable variable, double value)
{
    BallBerryCanopyState next;
    if (state == NULL || !isfinite(value))
        return -1;
    next = *state;
    switch (variable) {
    case BB_STATE_LAI:
        next.lai = value;
        break;
    case BB_STATE_SAI:
        next.sai = value;
        break;
    case BB_STATE_SHDFAC:
        next.shdfac = value;
        break;
    default:
        return -1;
    }
    if (ball_berry_state_validate(&next) != 0)
        return -1;
    *state = next;
    return 0;
}
