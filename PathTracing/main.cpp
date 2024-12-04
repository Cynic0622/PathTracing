#include "sphere.hpp"
#include <random>
#include <iostream>
#include <vector>
// #include "cube.hpp"

std::vector<Sphere> spheres = {
    // radius, position, color, emission, material
	Sphere(1e5, Vec(1e5 + 1, 40.8, 81.6), Vec(.75, .25, .25), Vec(0, 0, 0), "DIFF", 1),
	Sphere(1e5, Vec(-1e5 + 99, 40.8, 81.6), Vec(.25, .25, .75), Vec(0, 0, 0), "DIFF", 1),
    Sphere(1e5, Vec(50, 40.8, 1e5), Vec(.75, .75, .75), Vec(0, 0, 0), "DIFF", 1),
    Sphere(1e5, Vec(50, 40.8, -1e5 + 170), Vec(0, 0, 0), Vec(0, 0, 0), "DIFF", 1),
    // Sphere(600, Vec(50, 40.8, -511.33 - 75 - 10 - 3.5), Vec(0, 0, 0), Vec(1, 1, 0), "DIFF"),
    Sphere(1e5, Vec(50, 1e5, 81.6), Vec(.75, .75, .75), Vec(0, 0, 0), "DIFF", 1),
    Sphere(1e5, Vec(50, -1e5 + 81.6, 81.6), Vec(.75, .75, .75), Vec(0, 0, 0), "DIFF", 1),
    Sphere(16.5, Vec(27, 16.5, 47), Vec(3, 1, 1) * .999, Vec(0, 0, 0), "GlOSSY", 0.01),
    Sphere(16.5, Vec(73, 16.5, 78), Vec(1, 1, 1) * .999, Vec(0, 0, 0), "GLOSSY", 0.99),
    Sphere(600, Vec(50, 681.6 - .27, 81.6), Vec(0, 0, 0), Vec(12, 12, 12), "DIFF", 1),
    // glossy
    // Sphere(16.5, Vec(73, 50.0, 78), Vec(30, 10, 10), Vec(0, 0, 0), "GLOSSY"),
    Sphere(6, Vec(27 - 12, 50.0, 78), Vec(30, 10, 10), Vec(0, 0, 0), "GLOSSY", 0.1),
    Sphere(6, Vec(27, 50.0, 78), Vec(30, 10, 10), Vec(0, 0, 0), "GLOSSY", 0.3),
    Sphere(6, Vec(39, 50.0, 78), Vec(30, 10, 10), Vec(0, 0, 0), "GLOSSY", 0.4),
    Sphere(6, Vec(39 + 12, 50.0, 78), Vec(30, 10, 10), Vec(0, 0, 0), "GLOSSY", 0.5),
    Sphere(6, Vec(39 + 24, 50.0, 78), Vec(30, 10, 10), Vec(0, 0, 0), "GLOSSY", 0.6),
    Sphere(6, Vec(39 + 36, 50.0, 78), Vec(30, 10, 10), Vec(0, 0, 0), "GLOSSY", 0.7),
    Sphere(6, Vec(39 + 48, 50.0, 78), Vec(30, 10, 10), Vec(0, 0, 0), "GLOSSY", 0.8)
    
};
// Vec Cen(50,40.8,-860);
// std::vector<Sphere> spheres = {//Scene: radius, position, emission, color, material
//   // center 50 40.8 62
//   // floor 0
//   // back  0

//    Sphere(1600, Vec(1,0,2)*3000, Vec(1,.9,.8)*1.2e1*1.56*2,Vec(), "DIFF", 1), // sun
//    Sphere(1560, Vec(1,0,2)*3500,Vec(1,.5,.05)*4.8e1*1.56*2, Vec(),  "DIFF", 1), // horizon sun2
// //   Sphere(10000,Cen+Vec(0,0,-200), Vec(0.0627, 0.188, 0.569)*6e-2*8, Vec(.7,.7,1)*.25,  DIFF), // sky
//    Sphere(10000,Cen+Vec(0,0,-200), Vec(0.00063842, 0.02001478, 0.28923243)*6e-2*8, Vec(.7,.7,1)*.25,  "DIFF", 1), // sky

//   Sphere(100000, Vec(50, -100000, 0),  Vec(),Vec(.3,.3,.3),"DIFF", 1), // grnd
//   Sphere(110000, Vec(50, -110048.5, 0),  Vec(.9,.5,.05)*4,Vec(),"DIFF", 1),// horizon brightener
//   Sphere(4e4, Vec(50, -4e4-30, -3000),  Vec(),Vec(.2,.2,.2),"DIFF", 1),// mountains
// //  Sphere(3.99e4, Vec(50, -3.99e4+20.045, -3000),  Vec(),Vec(.7,.7,.7),DIFF),// mountains snow

//    Sphere(26.5,Vec(22,26.5,42),   Vec(),Vec(1,1,1)*.596, "SPEC", 0), // white Mirr
//    Sphere(13,Vec(75,13,82),   Vec(),Vec(.96,.96,.96)*.96, "REFR", 0),// Glas
//   Sphere(22,Vec(87,22,24),   Vec(),Vec(.6,.6,.6)*.696, "REFR", 0)    // Glas2
// };

// std::vector<Cube> cubes = {
//     // color, emission, material, p1, p2, p3, p4, p5, p6, p7, p8
//     // w 3 h 3 d 3
//     Cube(Vec(0, 0, 0), Vec(1, 1, 0), "DIFF", Vec(73 - 34, 50.0 - 2, 78 + 2), Vec(73 - 34 + 4, 50.0 - 2, 80), Vec(73 - 34, 50 + 2, 80), Vec(73 - 34 + 4, 50.0 + 2, 80),
//         Vec(73 - 34, 50.0 - 2, 78 - 2), Vec(73 - 34 + 4, 50.0 - 2, 78 - 2), Vec(73 - 34, 50 + 2, 78 - 2), Vec(73 - 34 + 4, 50 + 2, 78 - 2))
// };

inline double erand48()
{
    static std::random_device rd;
    static std::mt19937 gen(rd());
    static std::uniform_real_distribution<double> dis(0, 1);
    return dis(gen);
}
inline double clamp(double x)
{
    return x < 0 ? 0 : x > 1 ? 1 : x;
}
inline int toInt(double x)
{
    return int(pow(clamp(x), 1 / 2.2) * 255 + .5);
}
Vec importanceSampleGGX(float x, float y, float rougless, Vec N)
{
    float a = rougless * rougless;

    float phi = 2.0 * M_PI * x;
    float cosTheta = std::sqrt((1.0 - x) / (1.0 + (a * a - 1.0) * y));
    float sinTheta = std::sqrt(1.0 - cosTheta * cosTheta);

    Vec H;
    H.x = cos(phi) * sinTheta;
    H.y = sinTheta * sin(phi);
    H.z = cosTheta;

    // from tagent-space vector to world-space
    Vec up = abs(N.z) < 0.999 ? Vec(0.0, 0.0, 1.0) : Vec(1.0, 1.0, 1.0);
    Vec tangent = (N % H).norm();
    Vec bitangent = tangent % N;

    Vec sampleVec = tangent * H.x + bitangent * H.y + N * H.z; 
    return sampleVec.norm();
}

bool cubeInter = false;
bool intersect(const Ray& r, double& t, int& id)
{
    t = 1e20;
    // double eps = 1e-4;
    bool res = false;
    // int n = 9;
    for (int i = 0; i < spheres.size(); ++i)
    {
        double t0 = spheres[i].intersect(r);
        if (t0 && t0 < t)
        {
            t = t0;
            id = i;
            res = true;
        }
    }
    // float inter;
    
    // for (int i = 0; i < cubes.size(); ++i)
    // {
    //     bool flag = cubes[i].intersect(r, inter);
    //     if (flag && t > inter)
    //     {
    //         t = inter;
    //         id = i;
    //         res = true;
    //         cubeInter = true;
    //     }
    // }
    return res;
}
Vec radiance(Ray r, int depth)
{
    // std::cout << depth << std::endl;
    double t;
    int id = 0;
    if (!intersect(r, t, id)) return Vec(0, 0, 0); // missed

    const Sphere& obj = spheres[id];
    Vec x = r.o + r.d * t, n = (x - obj.pos).norm(), f = obj.color;
    Vec nl = n.dot(r.d) < 0 ? n : n * -1;
    // double p = std::fmax(f.x, std::fmax(f.y, f.z));
    double p = f.x > f.y && f.x > f.z ? f.x : f.y > f.z ? f.y : f.z;
    if (++depth > 5)
    {
        if (erand48() < p)
        {
            f = f * (1 / p);
        }
        else
        {
            return obj.emission;
        }
    }
    if (obj.material == "DIFF")
    { // nl作用， wuv坐标系
        Vec w = nl, u = ((fabs(w.x) > .1 ? Vec(0, 1, 0) : Vec(1, 0, 0)) % w).norm(), v = w % u;
        double r1 = 2 * M_PI * erand48(), r2 = erand48(), r2s = sqrt(r2);
        Vec d = (u * cos(r1) * r2s + v * sin(r1) * r2s + w * sqrt(1 - r2)).norm();
        return obj.emission + obj.color.mult(radiance(Ray(x, d), depth));
    }
    else if (obj.material == "SPEC")
    {
        // -d + refl = 2 * n * (-d.dot(n)) ==> refl = d - 2 * n * n.dot(d)
        Vec refl = r.d - n * 2 * n.dot(r.d);
        refl.norm();
        return obj.emission + obj.color.mult(radiance(Ray(x, refl),depth));
    }
    else if (obj.material == "GLOSSY")
    {
        auto sampleVec = importanceSampleGGX(erand48(), erand48(), obj.roughness, nl);
        return obj.emission + obj.color.mult(radiance(Ray(x, sampleVec), depth));
    }
    // n1 * sin(theta1) = n2 * sin(theta2)  n1 = 1.0, n2 = 1.5, cos(theta1) = -d.dot(n)
    // sin(theta1) = sqrt(1 - cos(theta1)^2) ==> cos(theta2) = sqrt(1 - (n1 / n2)^2 * (1 - cos(theta1)^2)) 
    // sin(theta1) > n2 / n1 * sin(90) ==> 1 - cos(theta1)^2 > (n2 / n1)^2 ==> 1 - cos(theta1)^2 - (n2 / n1)^2 > 0
    double n1 = 1.0, n2 = 1.5;
    // bool into = n1 < n2 ? false : true;
    bool into = n.dot(nl) > 0;
    Vec refl = r.d - n * 2 * n.dot(r.d);
    // cos(theta2) = sqrt(1 - n1 * n1 / n2 / n2 * (1 - cos(theta1)^2)) < 0
    // if (1 - (-r.d.dot(nl)) * (-r.d.dot(nl)) - (n1 / n2) * (n1 / n2) > 0)
    double nnt = into ? n1 / n2 : n2 / n1, ddn = r.d.dot(nl), cos2t;
    if ((cos2t = 1 - nnt * nnt * (1 - ddn * ddn)) < 0)
    {
        return obj.emission + obj.color.mult(radiance(Ray(x, refl), depth));
    }
    // T = B * sin(theta) + (-N) * cos(theta)
    // B = (D - (-N) * (-D * N)).norm() ==> T = (D - (-N) * (-D * N)).norm() * sin(theta) + (-N) * cos(theta)
    // sin(theta) = n1 / n2 * sqrt(1 - (-d * n) * (-d * n))
    // Vec d = ((r.d - n * (r.d.dot(n))).norm() * (n1 / n2 * sqrt(1 - (-r.d.dot(n)) * (-r.d.dot(n))))) - n * sqrt(1 - n1 * n1 / n2 * n2 * n2 * (1 - (-r.d.dot(n)) * (-r.d.dot(n))));
    // d.norm();
    // B = (D - (-N) * (-ddn)).norm() = (D - (-N) *(-ddn)) / sqrt(1 - cos^2(theta))
    // Vec d = (r.d - nl * ddn) * nnt - n * sqrt(cos2t);
    Vec d = (r.d * nnt - n * ((into ? 1 : -1) * (ddn * nnt + sqrt(cos2t))));
    d.norm();
    // 当入射角为0时，反射比例R0为 (n1 - n2)^2 / (n1 + n2)^2
    // 当入射角为θ时，反射比例Rθ为 R0 + (1 - R0) * (1 - cos(θ))^5
    double a = n1 + n2, b = n1 - n2, R0 = b * b / (a * a), c = 1 - (into ? -ddn : d.dot(n)); // bug
    double Re = R0 + (1 - R0) * pow(c, 5); 
    // Tr = 1 - R
    double Tr = 1 - Re, P = .25 + .5 * Re, RP = Re / P, TP = Tr / (1 - P);
    if (depth <= 2)
    {
        return obj.emission + obj.color.mult(radiance(Ray(x, refl), depth) * Re + radiance(Ray(x, d), depth) * Tr);
    }
    else
    {
        double p = erand48();
        return obj.emission + obj.color.mult(p < P ? radiance(Ray(x, refl), depth) * RP : radiance(Ray(x, d), depth) * TP);
    }

}   

int main(int argc, char** argv)
{
    int w = 1024, h = 768, samps = 5;
    if (argc > 1) samps = atoi(argv[1]);
    // samps = 25;
    Ray cam(Vec(50, 52, 295.6), Vec(0, -0.042612, -1).norm());
    // std::cout << cam.d.x << cam.d.y << cam.d.z << std::endl;
    Vec cx = Vec(w * 0.5135 / h, 0, 0), cy = ((cx % cam.d).norm()) * 0.5135, res;
    auto c = new Vec[1024 * 768];
    // for (int i = 0; i < 10; ++i) std::cout << c[i].x << c[i].y << c[i].z << std::endl;
#pragma omp parallel for schedule(dynamic, 1) private(res)
    for(int y = 0; y < h; ++y)
    {
        fprintf(stderr, "\rRendering (%d spp) %5.2f%%", samps * 4, 100. * y / (h - 1));
        for(int x = 0; x < w; ++x)
        {
            for(int sy = 0, i = (h - y - 1) * w + x; sy < 2; ++sy)
            {
                for(int sx = 0; sx < 2; ++sx)
                {
                    res = Vec();
                    for(int s = 0; s < samps; ++s)
                    {
                        double r1 = 2 * erand48(), dx = r1 < 1 ? sqrt(r1) - 1 : 1 - sqrt(2 - r1);
                        double r2 = 2 * erand48(), dy = r2 < 1 ? sqrt(r2) - 1 : 1 - sqrt(2 - r2);
                        // Vec d = cx * ((x + 0.5 + dx) / w - 0.5) + cy * ((y + 0.5 + dy) / h - 0.5) + cam.d;
                        Vec d = cx * (((sx + .5 + dx) / 2 + x) / w - .5) + cy * (((sy + .5 + dy) / 2 + y) / h - .5) + cam.d;
                        // std::cout << std::endl << d.x << " " << d.y << " " << d.z << std::endl << cx.x << " " << cx.y << " " << cx.z << std::endl << cy.x << " " << cy.y << " " << cy.z << std::endl << cam.d.x << " " << cam.d.y << " " << cam.d.z << std::endl;
                        res = res + radiance(Ray(cam.o + d * 140, d.norm()), 0) * (1. / samps);
                    }
                    c[i] = c[i] + Vec(clamp(res.x), clamp(res.y), clamp(res.z)) * .25;
                }
            }
            // c[y * w + x] += res * .25;
        }
    }
    FILE *f = fopen("image.ppm", "w");
    fprintf(f, "P3\n%d %d\n%d\n", w, h, 255);
    for (int i = 0; i < w * h; i++)
        fprintf(f, "%d %d %d ", toInt(c[i].x), toInt(c[i].y), toInt(c[i].z));
	return 0;
}