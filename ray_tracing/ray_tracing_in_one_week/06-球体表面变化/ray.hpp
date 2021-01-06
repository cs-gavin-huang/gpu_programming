/*
 * @Author: geekli
 * @Date: 2021-01-07 00:02:33
 * @LastEditTime: 2021-01-07 00:02:34
 * @LastEditors: your name
 * @Description: 
 * @FilePath: /ray_tracing/ray_tracing_in_one_week/04-添加ray类/ray.hpp
 */
#ifndef RAY_H
#define RAY_H

#include "vec3.hpp"

class ray {
    public:
        ray() {}
        ray(const point3& origin, const vec3& direction)
            : orig(origin), dir(direction)
        {}

        point3 origin() const  { return orig; }
        vec3 direction() const { return dir; }

        point3 at(double t) const {
            return orig + t*dir;
        }

    public:
        point3 orig;
        vec3 dir;
};

#endif