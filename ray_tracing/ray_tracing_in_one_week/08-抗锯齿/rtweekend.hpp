/*
 * @Author: geekli
 * @Date: 2021-01-07 00:26:38
 * @LastEditTime: 2021-01-07 00:26:49
 * @LastEditors: your name
 * @Description: 
 * @FilePath: /ray_tracing/ray_tracing_in_one_week/07-多个球体功能（ground）/rtweekend.hpp
 */
#ifndef RTWEEKEND_H
#define RTWEEKEND_H

#include <cmath>
#include <limits>
#include <memory>


// Usings

using std::shared_ptr;
using std::make_shared;
using std::sqrt;

// Constants

const double infinity = std::numeric_limits<double>::infinity();
const double pi = 3.1415926535897932385;

// Utility Functions

inline double degrees_to_radians(double degrees) {
    return degrees * pi / 180.0;
}

// Common Headers

#include "ray.hpp"
#include "vec3.hpp"

#endif