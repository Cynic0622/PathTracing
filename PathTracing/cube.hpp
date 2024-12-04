#include "ray.hpp"
#include <string>
#include<vector>
#include<array>
class Cube {
private:
    Vec color;
    Vec emission;
    std::string material;
    std::vector<Vec> points;
    float pMin[3], pMax[3];

public:
    Cube(Vec _color, Vec _emission, std::string _material, Vec _p1, Vec _p2, Vec _p3, Vec _p4, Vec _p5, Vec _p6, Vec _p7, Vec _p8) : color(_color), emission(_emission), material(_material) 
    {
        points.push_back(_p1);
        points.push_back(_p2);
        points.push_back(_p3);
        points.push_back(_p4);
        points.push_back(_p5);
        points.push_back(_p6);
        points.push_back(_p7);
        points.push_back(_p8);
        pMin[0] = std::min(std::min(std::min(std::min(std::min(std::min(std::min(points[0].x, points[1].x), points[2].x), points[3].x), points[4].x), points[5].x), points[6].x), points[7].x);
        pMin[1] = std::min(std::min(std::min(std::min(std::min(std::min(std::min(points[0].y, points[1].y), points[2].y), points[3].y), points[4].y), points[5].y), points[6].y), points[7].y);
        pMin[2] = std::min(std::min(std::min(std::min(std::min(std::min(std::min(points[0].z, points[1].z), points[2].z), points[3].z), points[4].z), points[5].z), points[6].z), points[7].z);
        pMax[0] = std::max(std::max(std::max(std::max(std::max(std::max(std::max(points[0].x, points[1].x), points[2].x), points[3].x), points[4].x), points[5].x), points[6].x), points[7].x);
        pMax[1] = std::max(std::max(std::max(std::max(std::max(std::max(std::max(points[0].y, points[1].y), points[2].y), points[3].y), points[4].y), points[5].y), points[6].y), points[7].y);
        pMax[2] = std::max(std::max(std::max(std::max(std::max(std::max(std::max(points[0].z, points[1].z), points[2].z), points[3].z), points[4].z), points[5].z), points[6].z), points[7].z);
    }

    bool intersect(const Ray& r, float& inter) const
    {
        std::array<int, 3> sign;
        for (int i = 0; i < 3; i++)
        {
            sign[i] = (r.d.x > 0) ? 1 : 0;
        }

        float t_x_min = (pMin[0] - r.o.x) / r.d.x;
        float t_x_max = (pMax[0] - r.o.x) / r.d.x;
        float t_y_min = (pMin[1] - r.o.y) / r.d.y;
        float t_y_max = (pMax[1] - r.o.y) / r.d.y;
        float t_z_min = (pMin[2] - r.o.z) / r.d.z;
        float t_z_max = (pMax[2] - r.o.z) / r.d.z;

        if (!sign[0])
        {
            std::swap(t_x_min, t_x_max);
        }
        
        if (!sign[1])
        {
            std::swap(t_y_min, t_y_max);
        }

        if (!sign[2])
        {
            std::swap(t_z_min, t_z_max);
        }

        float t_enter = std::max(t_x_min, std::max(t_y_min, t_z_min));
        float t_exit = std::min(t_x_max, std::min(t_y_max, t_z_max));
        if (t_enter < t_exit && t_exit >= 0) 
        {
            inter = t_enter;
            return true;
        }
        return false;
    }
};