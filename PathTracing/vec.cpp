#include "vec.hpp"

Vec::Vec() : x(0), y(0), z(0) {
    // x = y = z = 0;
}
Vec::Vec(double _x, double _y, double _z) 
{
    x = _x; y = _y; z = _z;
}

Vec Vec::operator+(const Vec &b) const 
{
	return Vec(x + b.x, y + b.y, z + b.z);
}

Vec Vec::operator-(const Vec &b) const
{
    return Vec(x - b.x, y - b.y, z - b.z);
}

Vec Vec::operator*(const double &b) const
{
    return Vec(x * b, y * b, z * b);
}

double Vec::dot(const Vec &b) const
{
    return x * b.x + y * b.y + z * b.z;
}

Vec Vec::mult(const Vec &b) const
{
    return Vec(x * b.x, y * b.y, z * b.z);
}

Vec Vec::operator%(const Vec &b) const
{
    return Vec(y * b.z - z * b.y, z * b.x - x * b.z, x * b.y - y * b.x);
}

Vec& Vec::norm()
{
    return *this = *this * (1 / sqrt(this->dot(*this)));
} 
