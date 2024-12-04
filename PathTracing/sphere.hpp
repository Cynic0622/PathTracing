#include "ray.hpp"
#include <string>
class Sphere
{
public:
    double radius;
    Vec pos;
    Vec color;
    Vec emission;
    std::string material;
    float roughness;
    Sphere(double _radius, Vec _pos, Vec _color, Vec _emission, std::string _material, float _roughness) : radius(_radius), pos(_pos), color(_color), emission(_emission), material(_material), roughness(_roughness) {}
    double intersect(const Ray& r) const
    {
        // x = o + td, |x-p| = r, (o + td - p) * (o + td - p) - r^2 = 0 求t
        // det = b^2 - 4ac, t^2 * d^2 - 2(o - p)td + (o - p)^2 - r^2 = 0
        // det = 4d^2(o - p)^2 - 4d^2((o - p)^2 - r^2)
        // Vec op = pos - r.o;
        // double b = op.dot(r.d), a = r.d.dot(r.d), c = op.dot(op) - radius * radius, t1, t2, eps = 1e-4;
        // double det = 4 * b * b - 4 * a * c;
        // if (det < 0) return 0;
        // else det = sqrt(det);
        // t1 = (-b - det) / (2 * a);
        // t2 = (-b + det) / (2 * a);
        // if (t1 > eps && t1 < t2) return t1;
        // return t2 > eps && t2 < t1 ? t2 : 0;
        Vec op = pos - r.o; // Solve t^2*d.d + 2*t*(o-p).d + (o-p).(o-p)-R^2 = 0
        double t, eps = 1e-4, b = op.dot(r.d), det = b * b - op.dot(op) + radius * radius;
        if (det < 0)
            return 0;
        else
            det = sqrt(det);
        return (t = b - det) > eps ? t : ((t = b + det) > eps ? t : 0);
    }
};