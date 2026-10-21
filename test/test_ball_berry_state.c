#include "ball_berry_state.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>

int main(void)
{
    BallBerryCanopyState s = {2.0, 0.5, 0.8};
    assert(ball_berry_state_validate(&s) == 0);
    assert(ball_berry_state_set(&s, BB_STATE_LAI, 0.0) == 0);
    assert(s.lai == 0.0 && s.sai == 0.5 && s.shdfac == 0.8);
    assert(ball_berry_state_set(&s, BB_STATE_SAI, 1.1) == 0);
    assert(s.lai == 0.0 && s.sai == 1.1 && s.shdfac == 0.8);
    assert(ball_berry_state_set(&s, BB_STATE_SHDFAC, 0.3) == 0);
    assert(s.lai == 0.0 && s.sai == 1.1 && s.shdfac == 0.3);
    assert(ball_berry_state_set(&s, BB_STATE_SHDFAC, 1.1) == -1);
    assert(ball_berry_state_set(&s, BB_STATE_LAI, -0.1) == -1);
    assert(ball_berry_state_set(&s, BB_STATE_SAI, NAN) == -1);
    assert(ball_berry_state_set(&s, (BallBerryStateVariable)99, 0.2) == -1);
    assert(s.lai == 0.0 && s.sai == 1.1 && s.shdfac == 0.3);
    puts("PASS: independent LAI, SAI, SHDFAC updates and validation");
    return 0;
}
