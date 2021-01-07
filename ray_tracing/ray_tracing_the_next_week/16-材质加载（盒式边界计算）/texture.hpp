/*
 * @Author: geekli
 * @Date: 2021-01-07 14:28:34
 * @LastEditTime: 2021-01-07 14:28:35
 * @LastEditors: your name
 * @Description: 
 * @FilePath: /ray_tracing/ray_tracing_the_next_week/16-材质加载（盒式边界计算）/texture.hpp
 */
#ifndef TEXTURE_H
#define TEXTURE_H

#include "rtweekend.hpp"

class texture {
    public:
        virtual color value(double u, double v, const point3& p) const = 0;
};

class solid_color : public texture {
    public:
        solid_color() {}
        solid_color(color c) : color_value(c) {}

        solid_color(double red, double green, double blue)
          : solid_color(color(red,green,blue)) {}

        virtual color value(double u, double v, const vec3& p) const override {
            return color_value;
        }

    private:
        color color_value;
};

#endif