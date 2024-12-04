#include "vec.cpp"
class Ray
{
public:
    Vec o, d;
    // double t;
    Ray(Vec _o, Vec _d) : o(_o), d(_d) {}
};