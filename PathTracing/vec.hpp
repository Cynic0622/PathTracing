#include<cmath>

class Vec
{
public:
    double x, y, z;
public:
    Vec();
    Vec(double _x, double _y, double _z);
    Vec operator+(const Vec &b) const;
    Vec operator-(const Vec &b) const;
    double dot(const Vec &b) const;
    Vec operator%(const Vec &b) const;
    Vec operator*(const double &b) const;
    Vec mult(const Vec &b) const;
    Vec& norm();
};