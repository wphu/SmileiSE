#ifndef USERFUNCTIONS_H
#define USERFUNCTIONS_H

class userFunctions {

public:

    static double erfinv(double x);
    static double erfinv2(double x);
    static bool is_nan(double dVal);
    static int isFiniteNumber(double d);

    //generate a gaussian random number R: mean=0, mean square=1.0
    //Then you can get another gaussian random number R2 = R*V + E: mean=E, mean square=V;
    static double gaussrand();
private:


};
#endif
