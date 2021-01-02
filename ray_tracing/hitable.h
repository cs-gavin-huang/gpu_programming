/*
 * @Author: geekli
 * @Date: 2021-01-02 08:29:54
 * @LastEditTime: 2021-01-02 08:30:13
 * @LastEditors: your name
 * @Description: 
 * @FilePath: /ray_tracing/hitable.h
 */
#ifndef HITABLE_H
#define HITABLE_H

#include "ray.h"

class material;

struct hit_record
{
    float t;
    vec3 p;
    vec3 normal;
    material *mat_ptr;
};

class hitable
{
public:
    virtual bool hit(const ray& r,float t_min,float t_max,hit_record & rec)const =0;
};


#endif 