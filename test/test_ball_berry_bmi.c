#include <stdio.h>
#include <math.h>
#include "bmi_pet.h"

static int check(Bmi *bmi, const char *name, double good, double bad)
{
    double value = -1.0;
    int index = 0;
    if (bmi->set_value(bmi, name, &good) != BMI_SUCCESS) return 1;
    if (bmi->get_value(bmi, name, &value) != BMI_SUCCESS) return 2;
    if (value != good) return 3;
    if (bmi->set_value(bmi, name, &bad) != BMI_FAILURE) return 4;
    if (bmi->get_value(bmi, name, &value) != BMI_SUCCESS || value != good) return 5;
    if (bmi->set_value_at_indices(bmi, name, &index, 1, &bad) != BMI_FAILURE) return 6;
    if (bmi->get_value(bmi, name, &value) != BMI_SUCCESS || value != good) return 7;
    return 0;
}

int main(void)
{
    Bmi bmi = {0};
    int count = 0;
    if (!register_bmi_pet(&bmi)) return 1;
    if (bmi.get_input_item_count(&bmi, &count) != BMI_SUCCESS || count != 13) return 2;
    const char *names[] = {"LAI", "SAI", "SHDFAC", "CWP", "MP", "VCMX25"};
    const double good[] = {2.3, 0.7, 0.6, 0.67, 9.0, 60.0};
    const double bad[] = {-1.0, -1.0, 1.2, -0.1, -2.0, -3.0};
    for (int i = 0; i < 6; ++i) {
        int result = check(&bmi, names[i], good[i], bad[i]);
        if (result != 0) {
            fprintf(stderr, "FAIL: %s test %d\n", names[i], result);
            return 3;
        }
    }
    double nan_value = NAN;
    if (bmi.set_value(&bmi, "LAI", &nan_value) != BMI_FAILURE) return 4;
    for (int i = 0; i < 6; ++i) {
        double value = -1.0;
        if (bmi.get_value(&bmi, names[i], &value) != BMI_SUCCESS || value != good[i]) return 5;
    }
    if (bmi.finalize(&bmi) != BMI_SUCCESS) return 6;
    puts("PASS: six independent Ball-Berry BMI inputs and validation");
    return 0;
}
