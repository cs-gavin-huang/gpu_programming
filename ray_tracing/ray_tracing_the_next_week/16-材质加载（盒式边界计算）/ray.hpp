/*
 * @Author: geekli
 * @Date: 2021-01-07 00:02:33
 * @LastEditTime: 2021-01-07 01:49:51
 * @LastEditors: your name
 * @Description: 
 * @FilePath: /ray_tracing/ray_tracing_the_next_week/15-动态模糊/ray.hpp
 */
#ifndef RAY_H
#define RAY_H

#include "vec3.hpp"

class ray {
    public:
        ray() {}
        ray(const point3& origin, const vec3& direction, double time = 0.0)
            : orig(origin), dir(direction), tm(time)
        {}

        point3 origin() const  { return orig; }
        vec3 direction() const { return dir; }
        double time() const    { return tm; }
        
        point3 at(double t) const {
            return orig + t*dir;
        }

    public:
        point3 orig;
        vec3 dir;
        double tm;
};

#endif