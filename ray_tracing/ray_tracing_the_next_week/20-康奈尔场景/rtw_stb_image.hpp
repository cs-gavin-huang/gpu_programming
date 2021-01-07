/*
 * @Author: geekli
 * @Date: 2021-01-07 15:24:25
 * @LastEditTime: 2021-01-07 15:25:43
 * @LastEditors: your name
 * @Description: 
 * @FilePath: /ray_tracing/ray_tracing_the_next_week/18-图片纹理加载/rtw_stb_image.hpp
 */
#ifndef RTWEEKEND_STB_IMAGE_H
#define RTWEEKEND_STB_IMAGE_H

// Disable pedantic warnings for this external library.
#ifdef _MSC_VER
    // Microsoft Visual C++ Compiler
    #pragma warning (push, 0)
#endif

#define STB_IMAGE_IMPLEMENTATION
//外部的库
#include "external/stb_image.h"

// Restore warning levels.
#ifdef _MSC_VER
    // Microsoft Visual C++ Compiler
    #pragma warning (pop)
#endif

#endif