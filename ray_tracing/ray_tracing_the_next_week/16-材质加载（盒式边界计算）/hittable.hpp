/*
 * @Author: geekli
 * @Date: 2021-01-07 00:22:03
 * @LastEditTime: 2021-01-07 14:29:04
 * @LastEditors: your name
 * @Description: 
 * @FilePath: /ray_tracing/ray_tracing_the_next_week/16-材质加载（盒式边界计算）/hittable.hpp
 */
#ifndef HITTABLE_H
#define HITTABLE_H

#include "ray.hpp"
#include "rtweekend.hpp"
#include "aabb.hpp"
class material;

struct hit_record {
    point3 p;
    vec3 normal;
    shared_ptr<material> mat_ptr;
    double t;
    double u;
    double v;
    bool front_face;

    inline void set_face_normal(const ray& r, const vec3& outward_normal) {
        front_face = dot(r.direction(), outward_normal) < 0;
        normal = front_face ? outward_normal :-outward_normal;
    }
};

class hittable {
    public:
        virtual bool hit(const ray& r, double t_min, double t_max, hit_record& rec) const = 0;
        //边界盒
        virtual bool bounding_box(double time0, double time1, aabb& output_box) const = 0;

};

#endif