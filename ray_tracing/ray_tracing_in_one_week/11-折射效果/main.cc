/*
 * @Author: geekli
 * @Date: 2021-01-06 23:55:02
 * @LastEditTime: 2021-01-07 01:26:18
 * @LastEditors: your name
 * @Description: 
 * @FilePath: /ray_tracing/ray_tracing_in_one_week/11-折射效果/main.cc
 */
//inclue 结构发生变化
#include "rtweekend.hpp"
#include "color.hpp"
#include "hittable_list.hpp"
#include "sphere.hpp"
#include "camera.hpp"
#include "material.hpp"
#include <iostream>

color ray_color(const ray& r, const hittable& world, int depth) {
    hit_record rec;

    // 光线放射限制，不再收集更多光线
    if (depth <= 0)
        return color(0,0,0);

    //0.001 避免太过于接近0的情况 0.0000001 速度快了很多
    if (world.hit(r, 0.001, infinity, rec)) {
        // point3 target = rec.p + rec.normal + random_unit_vector();
        // point3 target = rec.p + random_in_hemisphere(rec.normal);

        ray scattered;
        color attenuation;
        if (rec.mat_ptr->scatter(r, rec, attenuation, scattered))
            //递归将深度减一
            return attenuation * ray_color(scattered, world, depth-1);
        return color(0,0,0);
    }
    vec3 unit_direction = unit_vector(r.direction());
    auto t = 0.5*(unit_direction.y() + 1.0);
    return (1.0-t)*color(1.0, 1.0, 1.0) + t*color(0.5, 0.7, 1.0);
}

int main() {

    // Image
    const auto aspect_ratio = 16.0 / 9.0;
    const int image_width = 400;
    const int image_height = static_cast<int>(image_width / aspect_ratio);
    const int samples_per_pixel = 100;
    //深度50
    const int max_depth = 50;

    // World
    hittable_list world;
    auto material_ground = make_shared<lambertian>(color(0.8, 0.8, 0.0));
    /*
    auto material_center = make_shared<lambertian>(color(0.7, 0.3, 0.3));
    //散射值材质加成
    auto material_left   = make_shared<metal>(color(0.8, 0.8, 0.8), 0.3);
    */
    auto material_center = make_shared<dielectric>(1.5);
    auto material_left   = make_shared<dielectric>(1.5);
    auto material_right  = make_shared<metal>(color(0.8, 0.6, 0.2), 1.0);

    world.add(make_shared<sphere>(point3( 0.0, -100.5, -1.0), 100.0, material_ground));
    world.add(make_shared<sphere>(point3( 0.0,    0.0, -1.0),   0.5, material_center));
    world.add(make_shared<sphere>(point3(-1.0,    0.0, -1.0),   0.5, material_left));
    world.add(make_shared<sphere>(point3(-1.0,    0.0, -1.0),  -0.4, material_left));
    world.add(make_shared<sphere>(point3( 1.0,    0.0, -1.0),   0.5, material_right));

    // Camera
    camera cam;

    // Render

    std::cout << "P3\n" << image_width << ' ' << image_height << "\n255\n";

    for (int j = image_height-1; j >= 0; --j) {
        std::cerr << "\rScanlines remaining: " << j << ' ' << std::flush;
        for (int i = 0; i < image_width; ++i) {
            color pixel_color(0, 0, 0);
            for (int s = 0; s < samples_per_pixel; ++s) {
                auto u = (i + random_double()) / (image_width-1);
                auto v = (j + random_double()) / (image_height-1);
                ray r = cam.get_ray(u, v);
                //参数加上深度
                pixel_color += ray_color(r, world, max_depth);
            }
            //标准化输出颜色 + 抗锯齿
            write_color(std::cout, pixel_color, samples_per_pixel);
        }
    }

    std::cerr << "\nDone.\n";
}