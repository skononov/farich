#ifndef Minuit2FCN_H
# define Minuit2FCN_H

#include <vector>
#include "Minuit2/FCNBase.h"

class Minuit2FCN : public ROOT::Minuit2::FCNBase
{
private:
    double errorDef;

public:
    Minuit2FCN(double e=0.1) : errorDef(e) {}
    ~Minuit2FCN() {}

    double Up() const { return errorDef; }
    void SetErrorDef(double e) { errorDef=e; }

    double operator()(const std::vector<double>& x) const;
};

#endif
