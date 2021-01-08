## 输出图像

![img](https://pic2.zhimg.com/80/v2-d23bb56acdadb81f85b89d8fc74f48e1_720w.jpg)



输出这该格式的C++代码:

```cpp
#include <iostream>

int main() {
    const int image_width = 200;
    const int image_height = 100;

    std::cout << "P3\n" << image_width << ' ' << image_height << "\n255\n";

    for (int j = image_height-1; j >= 0; --j) {
        for (int i = 0; i < image_width; ++i) {
            auto r = double(i) / image_width;
            auto g = double(j) / image_height;
            auto b = 0.2;
            int ir = static_cast<int>(255.999 * r);
            int ig = static_cast<int>(255.999 * g);
            int ib = static_cast<int>(255.999 * b);
            std::cout << ir << ' ' << ig << ' ' << ib << '\n';
        }
    }
}
```

代码里有一些要注意的事情:

1. 对于像素来说, 每一行是从左往右写入的。
2. 行从上开始往下写入的。
3. 通常把RGB通道的值限定在0.0到1.0。之后计算颜色值的时候将使用一个动态的范围, 这个范围并不是0到1。但是在使用这段代码输出图像之前, 将把颜色映射到0到1
4. 下方的红色从左到右由黑边红, 左侧的绿色从上到下由黑到绿。红+绿变黄, 所以右上角应该是黄的。



打开输出的文件(Mac系统, 用Viewer打开的）win默认的看图软件不支持ppm格式, 只要查一下"蜂蜜浏览器"装个就行。打开后的结果如下:

![img](https://pic3.zhimg.com/80/v2-fb2906fa70e6a50c6814f52351d1d8d2_720w.jpg)

```text
P3
200 100
255
0 253 51
1 253 51
2 253 51
3 253 51
5 253 51
6 253 51
7 253 51
8 253 51
```



## vec3向量类

​	几乎所有的图形程序都使用类似的类来储存几何向量和颜色。这些向量是四维的(对于位置或者几何向量来说是三维的齐次拓展, 对于颜色来说是RGB加透明通道)。

用一个`vec3`类来储存所有的颜色, 位置, 方向, 位置偏移。

`vec3`的头文件:

```cpp
#include <iostream>

class vec3 {
    public:
        vec3() : e{0,0,0} {}
        vec3(double e0, double e1, double e2) : e{e0, e1, e2} {}

        double x() const { return e[0]; }
        double y() const { return e[1]; }
        double z() const { return e[2]; }

        vec3 operator-() const { return vec3(-e[0], -e[1], -e[2]); }
        double operator[](int i) const { return e[i]; }
        double& operator[](int i) { return e[i]; }

        vec3& operator+=(const vec3 &v) {
            e[0] += v.e[0];
            e[1] += v.e[1];
            e[2] += v.e[2];
            return *this;
        }

        vec3& operator*=(const double t) {
            e[0] *= t;
            e[1] *= t;
            e[2] *= t;
            return *this;
        }

        vec3& operator/=(const double t) {
            return *this *= 1/t;
        }

        double length() const {
            return sqrt(length_squared());
        }

        double length_squared() const {
            return e[0]*e[0] + e[1]*e[1] + e[2]*e[2];
        }

        void write_color(std::ostream &out) {
            // Write the translated [0,255] value of each color component.
            out << static_cast<int>(255.999 * e[0]) << ' '
                << static_cast<int>(255.999 * e[1]) << ' '
                << static_cast<int>(255.999 * e[2]) << '\n';
        }

    public:
        double e[3];
};
```

使用双精度浮点double, 有些光线追踪器使用单精度浮点float



头文件的第二部分包括一些向量操作工具函数

```cpp
inline std::ostream& operator<<(std::ostream &out, const vec3 &v) {
    return out << v.e[0] << ' ' << v.e[1] << ' ' << v.e[2];
}

inline vec3 operator+(const vec3 &u, const vec3 &v) {
    return vec3(u.e[0] + v.e[0], u.e[1] + v.e[1], u.e[2] + v.e[2]);
}

inline vec3 operator-(const vec3 &u, const vec3 &v) {
    return vec3(u.e[0] - v.e[0], u.e[1] - v.e[1], u.e[2] - v.e[2]);
}

inline vec3 operator*(const vec3 &u, const vec3 &v) {
    return vec3(u.e[0] * v.e[0], u.e[1] * v.e[1], u.e[2] * v.e[2]);
}

inline vec3 operator*(double t, const vec3 &v) {
    return vec3(t*v.e[0], t*v.e[1], t*v.e[2]);
}

inline vec3 operator*(const vec3 &v, double t) {
    return t * v;
}

inline vec3 operator/(vec3 v, double t) {
    return (1/t) * v;
}

inline double dot(const vec3 &u, const vec3 &v) {
    return u.e[0] * v.e[0]
         + u.e[1] * v.e[1]
         + u.e[2] * v.e[2];
}

inline vec3 cross(const vec3 &u, const vec3 &v) {
    return vec3(u.e[1] * v.e[2] - u.e[2] * v.e[1],
                u.e[2] * v.e[0] - u.e[0] * v.e[2],
                u.e[0] * v.e[1] - u.e[1] * v.e[0]);
}

inline vec3 unit_vector(vec3 v) {
    return v / v.length();
}
```



main函数改成这样

```cpp
#include "vec3.hpp"

#include <iostream>

int main() {
    const int image_width = 200;
    const int image_height = 100;

    std::cout << "P3\n" << image_width << ' ' << image_height << "\n255\n";

    for (int j = image_height-1; j >= 0; --j) {
        std::cerr << "\rScanlines remaining: " << j << ' ' << std::flush;
        for (int i = 0; i < image_width; ++i) {
            vec3 color(double(i)/image_width, double(j)/image_height, 0.2);
            color.write_color(std::cout);
        }
    }

    std::cerr << "\nDone.\n";
}
```



## 光线, 简单摄像机, 以及背景

所有的光线追踪器都有个一个ray类, 假定光线的公式为 ![[公式]](https://www.zhihu.com/equation?tex=%5Cmathbf%7Bp%7D%28t%29+%3D+%5Cmathbf%7Ba%7D+%2B+t+%5Cvec%7B%5Cmathbf%7Bb%7D%7D) 。



 ![[公式]](https://www.zhihu.com/equation?tex=%5Cmathbf%7Bp%7D) 是三维射线上的一个点。![[公式]](https://www.zhihu.com/equation?tex=%5Cmathbf%7Ba%7D) 是射线的原点, ![[公式]](https://www.zhihu.com/equation?tex=%5Cvec%7B%5Cmathbf%7Bb%7D%7D) 是射线的方向。



类中的变量 ![[公式]](https://www.zhihu.com/equation?tex=t) 是一个实数(代码中为double类型)。 ![[公式]](https://www.zhihu.com/equation?tex=p%28t%29) 接受任意的 ![[公式]](https://www.zhihu.com/equation?tex=t) 做为变量, 返回射线上的对应点。如果允许 ![[公式]](https://www.zhihu.com/equation?tex=t) 取负值可以得到整条直线。



对于一个正数 ![[公式]](https://www.zhihu.com/equation?tex=t) , 只能得到原点前部分 ![[公式]](https://www.zhihu.com/equation?tex=%5Cmathbf%7Ba%7D) , 这常常被称为半条直线, 或者说射线。

![img](https://pic3.zhimg.com/80/v2-0a13240fad9e360645339165054c4e62_720w.jpg)





在代码中使用复杂命名, 将函数p(t)扩写为ray::at(t)

```cpp
#ifndef RAY_H
#define RAY_H

#include "vec3.hpp"

class ray {
    public:
        ray() {}
        ray(const vec3& origin, const vec3& direction)
            : orig(origin), dir(direction)
        {}

        vec3 origin() const    { return orig; }
        vec3 direction() const { return dir; }

        vec3 at(double t) const {
            return orig + t*dir;
        }

    public:
        vec3 orig;
        vec3 dir;
};
#endif
```



光线追踪器的核心是从像素发射射线, 并计算这些射线得到的颜色。这包括如下的步骤:

1. 将射线从视点转化为像素坐标
2. 计算光线是否与场景中的物体相交
3. 如果有, 计算交点的颜色。在做光线追踪器的初期, 弄个简单摄像机让代码能跑起来，编写一个简单的color(ray)函数来返回背景颜色值(一个简单的渐变色)



射线r现在只是近似的从各个像素的中心射出

```cpp
#include "ray.hpp"

#include <iostream>

vec3 ray_color(const ray& r) {
    vec3 unit_direction = unit_vector(r.direction());
    auto t = 0.5*(unit_direction.y() + 1.0);
    return (1.0-t)*vec3(1.0, 1.0, 1.0) + t*vec3(0.5, 0.7, 1.0);
}

int main() {
    const int image_width = 200;
    const int image_height = 100;

    std::cout << "P3\n" << image_width << " " << image_height << "\n255\n";
    vec3 lower_left_corner(-2.0, -1.0, -1.0);
    vec3 horizontal(4.0, 0.0, 0.0);
    vec3 vertical(0.0, 2.0, 0.0);
    vec3 origin(0.0, 0.0, 0.0);
    for (int j = image_height-1; j >= 0; --j) {
        std::cerr << "\rScanlines remaining: " << j << ' ' << std::flush;
        for (int i = 0; i < image_width; ++i) {
            auto u = double(i) / image_width;
            auto v = double(j) / image_height;
            ray r(origin, lower_left_corner + u*horizontal + v*vertical);
            vec3 color = ray_color(r);
            color.write_color(std::cout);
        }
    }

    std::cerr << "\nDone.\n";
}
```

ray_color(ray)函数根据y值将蓝白做了个线性插值的混合, 这里把射线做了个单位化, 以保证y的取值范围(-1.0<y<1.0)。因为使用y轴做渐变, 所以可以看到这个蓝白渐变也是竖直的。

接下来使用了一个标准的小技巧将y的范围从-1.0 ≤ y ≤ 1.0映射到了0 ≤ y ≤ 1.0。

t=1.0时就是蓝色, 而t=0.0时就是白色。

在蓝白之间想要一个混合效果(blend)。现在采用的是线性混合(linear blend)或者说线性插值(liner interpolation)。或者简称其为lerp。

一个lerp一般来说会是下面的形式:

![[公式]](https://www.zhihu.com/equation?tex=%5Ctext%7BblendedValue%7D+%3D+%281-t%29%5Ccdot%5Ctext%7BstartValue%7D+%2B+t%5Ccdot%5Ctext%7BendValue%7D)

当t从0到1, 会渲染出这样的图像



![sKXB9S.png](https://s3.ax1x.com/2021/01/09/sKXB9S.png)



## 加入球体



在图形学中,  对于到球面上的点 ![[公式]](https://www.zhihu.com/equation?tex=%5Cmathbf%7BP%7D+%3D+%28x%2Cy%2Cz%29) 到球心 ![[公式]](https://www.zhihu.com/equation?tex=%5Cmathbf%7Bc%7D+%3D+%28%5Cmathbf%7Bc%7D_x%2C%5Cmathbf%7Bc%7D_y%2C%5Cmathbf%7Bc%7D_z%29) 的距离可以使用向量表示为 ![[公式]](https://www.zhihu.com/equation?tex=%28%5Cmathbf%7Bp%7D+-+%5Cmathbf%7Bc%7D%29) , 于是就有

![[公式]](https://www.zhihu.com/equation?tex=%28%5Cmathbf%7Bp%7D+-+%5Cmathbf%7Bc%7D%29+%5Ccdot+%28%5Cmathbf%7Bp%7D+-+%5Cmathbf%7Bc%7D%29+%3D+%28x-%5Cmathbf%7Bc%7D_x%29%5E2+%2B+%28y-%5Cmathbf%7Bc%7D_y%29%5E2+%2B+%28z-%5Cmathbf%7Bc%7D_z%29%5E2)

得到球面方程的向量形式:

![[公式]](https://www.zhihu.com/equation?tex=+%28%5Cmathbf%7Bp%7D+-+%5Cmathbf%7Bc%7D%29+%5Ccdot+%28%5Cmathbf%7Bp%7D+-+%5Cmathbf%7Bc%7D%29+%3D+R%5E2+)

可以将其解读为"满足方程上述方程的任意一点p一定位于球面上"。

要知道射线 ![[公式]](https://www.zhihu.com/equation?tex=p%28t%29+%3D+%5Cmathbf%7Ba%7D+%2B+t%5Cvec%7B%5Cmathbf%7Bb%7D%7D) 是否与球体相交。

​	如果说它相交了, 那么肯定有一个t使直线上的点p(t)满足球面方程。

先来计算满足条件的任意t值:

![[公式]](https://www.zhihu.com/equation?tex=%28p%28t%29+-+%5Cmathbf%7Bc%7D%29%5Ccdot%28p%28t%29+-+%5Cmathbf%7Bc%7D%29+%3D+R%5E2)

或者将p(t)展开为射线方程:

![[公式]](https://www.zhihu.com/equation?tex=%28%5Cmathbf%7Ba%7D+%2B+t+%5Cvec%7B%5Cmathbf%7Bb%7D%7D+-+%5Cmathbf%7Bc%7D%29+%5Ccdot+%28%5Cmathbf%7Ba%7D+%2B+t+%5Cvec%7B%5Cmathbf%7Bb%7D%7D+-+%5Cmathbf%7Bc%7D%29+%3D+R%5E2)

展开表达式并移项, 得:

![[公式]](https://www.zhihu.com/equation?tex=t%5E2+%5Cvec%7B%5Cmathbf%7Bb%7D%7D%5Ccdot%5Cvec%7B%5Cmathbf%7Bb%7D%7D+%2B+2t+%5Cvec%7B%5Cmathbf%7Bb%7D%7D+%5Ccdot+%5Cvec%7B%28%5Cmathbf%7Ba%7D-%5Cmathbf%7Bc%7D%29%7D+%2B+%5Cvec%7B%28%5Cmathbf%7Ba%7D-%5Cmathbf%7Bc%7D%29%7D+%5Ccdot+%5Cvec%7B%28%5Cmathbf%7Ba%7D-%5Cmathbf%7Bc%7D%29%7D+-+R%5E2+%3D+0)

方程中的向量和半径R都是已知的常量, 唯一的未知数就是t, 并且这个等式是关于t的一个一元二次方程, 可以用求根公式来判别交点个数, 为正则2个交点, 为负则无交点, 为0则一个交点。



使用代码来求解, 并使用红色来表示射线击中放在(0,0,-1)的小球:

```cpp
bool hit_sphere(const vec3& center, double radius, const ray& r) {
    vec3 oc = r.origin() - center;
    auto a = dot(r.direction(), r.direction());
    auto b = 2.0 * dot(oc, r.direction());
    auto c = dot(oc, oc) - radius*radius;
    auto discriminant = b*b - 4*a*c;
    return (discriminant > 0);
}

vec3 ray_color(const ray& r) {
    if (hit_sphere(vec3(0,0,-1), 0.5, r))
        return vec3(1, 0, 0);
    vec3 unit_direction = unit_vector(r.direction());
    auto t = 0.5*(unit_direction.y() + 1.0);
    return (1.0-t)*vec3(1.0, 1.0, 1.0) + t*vec3(0.5, 0.7, 1.0);
}
```



![img](https://pic4.zhimg.com/80/v2-73a372b979105344547ea89345970093_720w.jpg)

得到一个简单的红色球体



## 面法相与复数物体

为了来给球体着色, 首先来定义一下面法相。面法相应该是一种垂直于交点所在平面的三维向量。关于面法相存在两个设计抉择。首先是是否将其设计为单位向量, 这样对于着色器来说, 所以会说"yes!"但是并没有在代码里这么做, 这部分差异可能会导致一些潜在的bug。所以记住, 这个是个人喜好, 大多数的人喜好使用单位法相。对于球体来说, 朝外的法相是直线与球的交点减去球心:

![img](https://pic1.zhimg.com/80/v2-009e0ab00e758c137a5a9a6dbc639a8c_720w.jpg)

球体的表面法相

说到底, 其实就是从球心到交点再向外延伸的那个方向。让把这部分转变成代码并开始着色。暂时还没有光源这样的东西, 所以让直接将法相值作为颜色输出吧。对于法相可视化来说, 常常将xyz分量的值先映射到0到1的范围(假定 ![[公式]](https://www.zhihu.com/equation?tex=vecN) 是一个单位向量, 它的取值范围是-1到1的),再把它赋值给rgb。对于法相来说, 光能判断射线是否与球体相交是不够的, 还需求出交点的坐标。在有两个交点的情况下, 选取最近的交点smallest(t)。计算与可视化球的法向量的代码如下:

```cpp
//main.cc 球体表面法相
double hit_sphere(const vec3& center, double radius, const ray& r) {
    vec3 oc = r.origin() - center;
    auto a = dot(r.direction(), r.direction());
    auto b = 2.0 * dot(oc, r.direction());
    auto c = dot(oc, oc) - radius*radius;
    auto discriminant = b*b - 4*a*c;
    if (discriminant < 0) {
        return -1.0;
    } else {
        return (-b - sqrt(discriminant) ) / (2.0*a);
    }
}

vec3 ray_color(const ray& r) {
    auto t = hit_sphere(vec3(0,0,-1), 0.5, r);
    if (t > 0.0) {
        vec3 N = unit_vector(r.at(t) - vec3(0,0,-1));
        return 0.5*vec3(N.x()+1, N.y()+1, N.z()+1);
    }
    vec3 unit_direction = unit_vector(r.direction());
    t = 0.5*(unit_direction.y() + 1.0);
    return (1.0-t)*vec3(1.0, 1.0, 1.0) + t*vec3(0.5, 0.7, 1.0);
}
```

这会得到下面的结果:

![img](https://s3.ax1x.com/2021/01/09/sKXWNV.png)



使用法相作为颜色值输出

```cpp
//main.cc
vec3 oc = r.origin() - center;
auto a = dot(r.direction(), r.direction());
auto b = 2.0 * dot(oc, r.direction());
auto c = dot(oc, oc) - radius*radius;
auto discriminant = b*b - 4*a*c;
```

首先, 回想一下一个向量与自己的点积就是它的长度的平方(都是 ![[公式]](https://www.zhihu.com/equation?tex=x%5E2%2By%5E2%2Bz%5E2) )

其次, 注意其实的`b`有一个系数2, 设b=2h，有:

![[公式]](https://www.zhihu.com/equation?tex=%5Cfrac%7B-b+%5Cpm+%5Csqrt%7Bb%5E2+-+4ac%7D%7D%7B2a%7D)

![[公式]](https://www.zhihu.com/equation?tex=%3D+%5Cfrac%7B-2h+%5Cpm+%5Csqrt%7B%282h%29%5E2+-+4ac%7D%7D%7B2a%7D+)

![[公式]](https://www.zhihu.com/equation?tex=%3D+%5Cfrac%7B-2h+%5Cpm+2%5Csqrt%7Bh%5E2+-+ac%7D%7D%7B2a%7D)

![[公式]](https://www.zhihu.com/equation?tex=%3D+%5Cfrac%7B-h+%5Cpm+%5Csqrt%7Bh%5E2+-+ac%7D%7D%7Ba%7D)

射线与球体求交的代码简化:

```cpp
//main.cc
vec3 oc = r.origin() - center;
auto a = r.direction().length_squared();
auto half_b = dot(oc, r.direction());
auto c = oc.length_squared() - radius*radius;
auto discriminant = half_b*half_b - a*c;

if (discriminant < 0) {
    return -1.0;
} else {
    return (-half_b - sqrt(discriminant) ) / a;
}
```



多个物体（hittable）

hittable类理应有个接受射线为参数的函数, 许多光线追踪器为了便利, 加入了一个区间 ![[公式]](https://www.zhihu.com/equation?tex=t_%7Bmin%7D%3Ct%3Ct_%7Bmax%7D) 来判断相交是否有效。

对于一开始的光线来说, 这个t值总是正的。现在有个设计上的问题: 是否在每次计算求交的时候都要去计算法相?  但其实只需要计算离射线原点最近的那个交点的法相就行了, 后面的东西会被遮挡。

hittable抽象类(计算的结果存在的结构体)

```cpp
//hittable.hpp
#ifndef HITTABLE_H
#define HITTABLE_H

#include "ray.hpp"

struct hit_record {
    vec3 p;
    vec3 normal;
    double t;
};

class hittable {
    public:
        virtual bool hit(const ray& r, double t_min, double t_max, hit_record& rec) const = 0;
};

#endif
```

继承自它的sphere球体类:

```cpp
//sphere.hpp
#ifndef SPHERE_H
#define SPHERE_H

#include "hittable.hpp"
#include "vec3.hpp"

class sphere: public hittable {
    public:
        sphere() {}
        sphere(vec3 cen, double r) : center(cen), radius(r) {};

        virtual bool hit(const ray& r, double tmin, double tmax, hit_record& rec) const;

    public:
        vec3 center;
        double radius;
};

bool sphere::hit(const ray& r, double t_min, double t_max, hit_record& rec) const {
    vec3 oc = r.origin() - center;
    auto a = r.direction().length_squared();
    auto half_b = dot(oc, r.direction());
    auto c = oc.length_squared() - radius*radius;
    auto discriminant = half_b*half_b - a*c;

    if (discriminant > 0) {
        auto root = sqrt(discriminant);
        auto temp = (-half_b - root)/a;
        if (temp < t_max && temp > t_min) {
            rec.t = temp;
            rec.p = r.at(rec.t);
            rec.normal = (rec.p - center) / radius;
            return true;
        }
        temp = (-half_b + root) / a;
        if (temp < t_max && temp > t_min) {
            rec.t = temp;
            rec.p = r.at(rec.t);
            rec.normal = (rec.p - center) / radius;
            return true;
        }
    }
    return false;
}
#endif
```



面法相的朝向问题。

​	对于现在来说, 如果光线从球体外部击中球体, 那么法相也是朝外的, 与射线的方向相反(不是数学意义上的严格相反, 只是大致逆着)。如果光线从内部射向球面时, 此时的面法相依然朝外, 与射线方向相同。相对的, 也可以总是让法相向量与射线方向相反, 即射线从外部射向球面, 法向量朝外, 射线从内部射向球面, 法向量向着球心。

​	在着色前, 需要仔细考虑一下采用哪种方式, 

	1. 如果决定让法相永远朝外, 那在就得在射入的时候判断是从表面的哪一侧射入的, 可以简单的将光线与法相做点乘来判断。如果法相与光线方向相同, 那就是从内部击中内表面, 如果相反则是从外部击中外表面。
	2. 如果永远让法相与入射方向相反, 不用去用点乘来判断射入面是内侧还是外侧了, 这时需要用一个变量储存射入面的信息:

```cpp
bool front_face;
if (dot(ray_direction, outward_normal) > 0.0) {
    // ray is inside the sphere
    normal = -outward_normal;
    front_face = false;
}
else {
    // ray is outside the sphere
    normal = outward_normal;
    front_face = true;
}
```

在结构体hit_record中加入front_face变量

```cpp
//hittable.hpp 加入时间与面朝向
ifndef HITTABLE_H
#define HITTABLE_H

#include "ray.hpp"

struct hit_record {
    vec3 p;
    vec3 normal;
    double t;
    bool front_face;

    inline void set_face_normal(const ray& r, const vec3& outward_normal) {
        front_face = dot(r.direction(), outward_normal) < 0;
        normal = front_face ? outward_normal :-outward_normal;
    }
};

class hittable {
    public:
        virtual bool hit(const ray& r, double t_min, double t_max, hit_record& rec) const = 0;
};
#endif
```

求交时加入射入面的判别:

```cpp
//sphere.hpp 加入射入面判别
bool sphere::hit(const ray& r, double t_min, double t_max, hit_record& rec) const {
    vec3 oc = r.origin() - center;
    auto a = r.direction().length_squared();
    auto half_b = dot(oc, r.direction());
    auto c = oc.length_squared() - radius*radius;
    auto discriminant = half_b*half_b - a*c;

    if (discriminant > 0) {
        auto root = sqrt(discriminant);
        auto temp = (-half_b - root)/a;
        if (temp < t_max && temp > t_min) {
            rec.t = temp;
            rec.p = r.at(rec.t);
            vec3 outward_normal = (rec.p - center) / radius;
            rec.set_face_normal(r, outward_normal);
            return true;
        }
        temp = (-half_b + root) / a;
        if (temp < t_max && temp > t_min) {
            rec.t = temp;
            rec.p = r.at(rec.t);
            vec3 outward_normal = (rec.p - center) / radius;
            rec.set_face_normal(r, outward_normal);
            return true;
        }
    }
    return false;
}
```

加入存放物体的列表

```cpp
//hittable_list.hpp
#ifndef HITTABLE_LIST_H
#define HITTABLE_LIST_H

#include "hittable.hpp"
#include <memory>
#include <vector>

using std::shared_ptr;
using std::make_shared;

class hittable_list: public hittable {
    public:
        hittable_list() {}
        hittable_list(shared_ptr<hittable> object) { add(object); }

        void clear() { objects.clear(); }
        void add(shared_ptr<hittable> object) { objects.push_back(object); }

        virtual bool hit(const ray& r, double tmin, double tmax, hit_record& rec) const;

    public:
        std::vector<shared_ptr<hittable>> objects;
};

bool hittable_list::hit(const ray& r, double t_min, double t_max, hit_record& rec) const {
    hit_record temp_rec;
    bool hit_anything = false;
    auto closest_so_far = t_max;

    for (const auto& object : objects) {
        if (object->hit(r, t_min, closest_so_far, temp_rec)) {
            hit_anything = true;
            closest_so_far = temp_rec.t;
            rec = temp_rec;
        }
    }

    return hit_anything;
}

#endif
```



## 常用的常数与工具

需要在头文件中定义一些常用的常数。目前为止只需要定义无穷。但是先把pi在这里定义好, 之后要用的。对于pi来说并没有什么跨平台的标准定义,  所以自己来定义一下。



在rtweekend.h中给出了一些未来常用的常数和函数:

```cpp
//rtweekend.hpp
#ifndef RTWEEKEND_H
#define RTWEEKEND_H

#include <cmath>
#include <cstdlib>
#include <limits>
#include <memory>


// Usings

using std::shared_ptr;
using std::make_shared;

// Constants

const double infinity = std::numeric_limits<double>::infinity();
const double pi = 3.1415926535897932385;

// Utility Functions

inline double degrees_to_radians(double degrees) {
    return degrees * pi / 180;
}

inline double ffmin(double a, double b) { return a <= b ? a : b; }
inline double ffmax(double a, double b) { return a >= b ? a : b; }

// Common Headers

#include "ray.hpp"
#include "vec3.hpp"

#endif
```

更新main函数:

```cpp
//main.cc
#include "rtweekend.hpp"

#include "hittable_list.hpp"
#include "sphere.hpp"

#include <iostream>
vec3 ray_color(const ray& r, const hittable& world) {
    hit_record rec;
    if (world.hit(r, 0, infinity, rec)) {
        return 0.5 * (rec.normal + vec3(1,1,1));
    }
    vec3 unit_direction = unit_vector(r.direction());
    auto t = 0.5*(unit_direction.y() + 1.0);
    return (1.0-t)*vec3(1.0, 1.0, 1.0) + t*vec3(0.5, 0.7, 1.0);
}

int main() {
    const int image_width = 200;
    const int image_height = 100;

    std::cout << "P3\n" << image_width << ' ' << image_height << "\n255\n";

    vec3 lower_left_corner(-2.0, -1.0, -1.0);
    vec3 horizontal(4.0, 0.0, 0.0);
    vec3 vertical(0.0, 2.0, 0.0);
    vec3 origin(0.0, 0.0, 0.0);

    hittable_list world;
    world.add(make_shared<sphere>(vec3(0,0,-1), 0.5));
    world.add(make_shared<sphere>(vec3(0,-100.5,-1), 100));

    for (int j = image_height-1; j >= 0; --j) {
        std::cerr << "\rScanlines remaining: " << j << ' ' << std::flush;
        for (int i = 0; i < image_width; ++i) {
            auto u = double(i) / image_width;
            auto v = double(j) / image_height;
            ray r(origin, lower_left_corner + u*horizontal + v*vertical);

            vec3 color = ray_color(r, world);

            color.write_color(std::cout);
        }
    }

    std::cerr << "\nDone.\n";
}
```



得到一张使用法向作为球体颜色值的图片。



![sKXE6J.png](https://s3.ax1x.com/2021/01/09/sKXE6J.png)



## 抗锯齿

真实世界中的摄像机拍摄出来的照片是没有像素状的锯齿的。

​	因为边缘像素是由背景和前景混合而成的，可以在程序中简单的对每个边缘像素多次采样取平均达到类似的效果。

需要一个能够返回真随机数的一个随机数生成器。

​	默认来说这个函数应该返回0≤r<1的随机数



使用<cstdlib>中的rand()函数，返回0到RAND_MAX中的一个任意整数。

将下面的一小段代码加到rtweekend.h中, 得到想要的随机函数:

```cpp
//rtweekend.hpp
#include <cstdlib>
...

inline double random_double() {
    // Returns a random real in [0,1).
    return rand() / (RAND_MAX + 1.0);
}

inline double random_double(double min, double max) {
    // Returns a random real in [min,max).
    return min + (max-min)*random_double();
}
```

对于给定的像素, 发射多条射线进行多次采样，然后对颜色结果求一个平均值:

![img](https://pic2.zhimg.com/80/v2-e54c77c97131c90f04a27e5a5de7ec55_720w.jpg)

 对一个像素进行多次采样

简单的轴对齐摄像机类进行了一次封装:

```cpp
//camera.hpp
#ifndef CAMERA_H
#define CAMERA_H

#include "rtweekend.hpp"

class camera {
    public:
        camera() {
            lower_left_corner = vec3(-2.0, -1.0, -1.0);
            horizontal = vec3(4.0, 0.0, 0.0);
            vertical = vec3(0.0, 2.0, 0.0);
            origin = vec3(0.0, 0.0, 0.0);
        }

        ray get_ray(double u, double v) {
            return ray(origin, lower_left_corner + u*horizontal + v*vertical - origin);
        }

    public:
        vec3 origin;
        vec3 lower_left_corner;
        vec3 horizontal;
        vec3 vertical;
};
#endif
```

为了对多重采样的颜色值进行计算, 升级了vec3::write_color()函数。不会在每次发出射线采样时都计算一个0-1之间的颜色值, 而是一次性把所有的颜色都加在一起, 然后最后只需要简单的除以采样点个数。

```cpp
//vec3.hpp
...
#include "rtweekend.hpp"
...
void write_color(std::ostream &out, int samples_per_pixel) {
    // Divide the color total by the number of samples.
    auto scale = 1.0 / samples_per_pixel;
    auto r = scale * e[0];
    auto g = scale * e[1];
    auto b = scale * e[2];

    // Write the translated [0,255] value of each color component.
    out << static_cast<int>(256 * clamp(r, 0.0, 0.999)) << ' '
        << static_cast<int>(256 * clamp(g, 0.0, 0.999)) << ' '
        << static_cast<int>(256 * clamp(b, 0.0, 0.999)) << '\n';
}
```

头文件rtweekend.h加入了一个新函数clamp(x,min,max), 用来将x限制在[min,max]区间之中:

```cpp
//rtweekend.hpp
inline double clamp(double x, double min, double max) {
    if (x < min) return min;
    if (x > max) return max;
    return x;
}
```

main函数变化:

```cpp
//main.cc
int main() {
    const int image_width = 200;
    const int image_height = 100;
    const int samples_per_pixel = 100;

    std::cout << "P3\n" << image_width << " " << image_height << "\n255\n";

    hittable_list world;
    world.add(make_shared<sphere>(vec3(0,0,-1), 0.5));
    world.add(make_shared<sphere>(vec3(0,-100.5,-1), 100));
    camera cam;
    for (int j = image_height-1; j >= 0; --j) {
        std::cerr << "\rScanlines remaining: " << j << ' ' << std::flush;
        for (int i = 0; i < image_width; ++i) {
            vec3 color(0, 0, 0);
            for (int s = 0; s < samples_per_pixel; ++s) {
                auto u = (i + random_double()) / image_width;
                auto v = (j + random_double()) / image_height;
                ray r = cam.get_ray(u, v);
                color += ray_color(r, world);
            }
            color.write_color(std::cout, samples_per_pixel);
        }
    }

    std::cerr << "\nDone.\n";
}
```



![img](https://pic2.zhimg.com/80/v2-807248338e5835236f2c884fa9ed4105_720w.jpg)



加入抗锯齿效果后前背景颜色混合的像素





## 漫反射材质

漫反射材质不仅仅接受其周围环境的光线, 还会在散射时使光线变成自己本身的颜色。

光线射入漫反射材质后, 其反射方向是随机的。

如果为下面这两个漫发射的球射入三条光线, 光线都会有不同的反射角度:

![img](https://pic2.zhimg.com/80/v2-f61ad9d024071668f27db1fc248abcf5_720w.jpg)



并且大部分的光线都会被吸收, 而不是被反射。表面越暗, 吸收就越有可能发生。

使用任意的算法生成随机的反射方向, 就能让其看上去像一个粗糙不平的漫反射材质。

需要一个算法来生成球体内的随机点。会采用最简单的做法:否定法(rejection method)。

首先, 在一个xyz取值范围为-1到+1的单位立方体中选取一个随机点, 如果这个点在球外就重新生成直到该点在球内:

```cpp
//vec3.hpp
class vec3 {
  public:
    ...
    inline static vec3 random() {
        return vec3(random_double(), random_double(), random_double());
    }

    inline static vec3 random(double min, double max) {
        return vec3(random_double(min,max), random_double(min,max), random_double(min,max));
    }
```

```cpp
//vec3.hpp
vec3 random_in_unit_sphere() {
    while (true) {
        auto p = vec3::random(-1,1);
        if (p.length_squared() >= 1) continue;
        return p;
    }
}
```

更新一下的ray_color()函数（使用新的生成随机反射方向的函数来）:

```cpp
//main.cc
vec3 ray_color(const ray& r, const hittable& world) {
    hit_record rec;

    if (world.hit(r, 0, infinity, rec)) {
        vec3 target = rec.p + rec.normal + random_in_unit_sphere();
        return 0.5 * ray_color(ray(rec.p, target - rec.p), world);
    }

    vec3 unit_direction = unit_vector(r.direction());
    auto t = 0.5*(unit_direction.y() + 1.0);
    return (1.0-t)*vec3(1.0, 1.0, 1.0) + t*vec3(0.5, 0.7, 1.0);
}
```



注意ray_color函数是一个递归函数。

​	那么递归终止的条件是什么呢 ：没有击中任何东西。

在某些条件下, 达到这个终止条件的时间会非常长, 长到足够大于函数栈，为了避免这种情况的发生, 使用一个变量depth限制递归层数。

当递归层数达到限制值时终止递归, 返回黑色:

```cpp
//main.cc
vec3 ray_color(const ray& r, const hittable& world, int depth) {
    hit_record rec;

    if (depth <= 0)
        return vec3(0,0,0);

    if (world.hit(r, 0, infinity, rec)) {
        vec3 target = rec.p + rec.normal + random_in_unit_sphere();
        return 0.5 * ray_color(ray(rec.p, target - rec.p), world, depth-1);
    }

    vec3 unit_direction = unit_vector(r.direction());
    auto t = 0.5*(unit_direction.y() + 1.0);
    return (1.0-t)*vec3(1.0, 1.0, 1.0) + t*vec3(0.5, 0.7, 1.0);
}
...
int main() {
    const int image_width = 200;
    const int image_height = 100;
    const int samples_per_pixel = 100;
    const int max_depth = 50;

    ...
    for (int j = image_height-1; j >= 0; --j) {
        std::cerr << "\rScanlines remaining: " << j << ' ' << std::flush;
        for (int i = 0; i < image_width; ++i) {
            vec3 color(0, 0, 0);
            for (int s = 0; s < samples_per_pixel; ++s) {
                auto u = (i + random_double()) / image_width;
                auto v = (j + random_double()) / image_height;
                ray r = cam.get_ray(u, v);
               color += ray_color(r, world, max_depth);
            }
            color.write_color(std::cout, samples_per_pixel);
        }
    }

    std::cerr << "\nDone.\n";
}
```

得到:

![img](https://pic3.zhimg.com/80/v2-01259cbd92f891d2db907a773d51e7aa_720w.jpg)



第一次渲染出漫反射材质的球体

注意球下面是有影子的，非常的暗, 球在散射的时候只吸收了一半。

现实世界中的这个球明显是应该更加亮一些的，所有的看图软件都默认图像已经经过了伽马校正(gamma corrected)，即在图片存入字节之前, 颜色值发生了一次转化。

使用"gamma 2"空间, 就意味着最终的颜色值要加上指数 ![[公式]](https://www.zhihu.com/equation?tex=1%2Fgamma) , 在的例子里就是 ½, 即开平方根:

```cpp
//vec3.hpp
void write_color(std::ostream &out, int samples_per_pixel) {

    auto scale = 1.0 / samples_per_pixel;
    auto r = sqrt(scale * e[0]);
    auto g = sqrt(scale * e[1]);
    auto b = sqrt(scale * e[2]);

    out << static_cast<int>(256 * clamp(r, 0.0, 0.999)) << ' '
        << static_cast<int>(256 * clamp(g, 0.0, 0.999)) << ' '
        << static_cast<int>(256 * clamp(b, 0.0, 0.999)) << '\n';
}
```



![img](https://pic4.zhimg.com/80/v2-f0ab33e768ac857087b5e2547f47b8b3_720w.jpg)



伽马校正后的漫反射球体

潜在bug：有些物体反射的光线会在t=0时再次击中自己。然而由于精度问题, 这个值可能是t=-0.000001或者是t=0.0000000001或者任意接近0的浮点数。



要忽略掉0附近的一部分范围, 防止物体发出的光线再次与自己相交：

```cpp
//main.cc
if (world.hit(r, 0.001, infinity, rec)) {
```

避免阴影痤疮(shadow ance)的产生。

拒绝法生成的点是单位球体积内的的随机点, 这样生成的向量大概率上会和法线方向相近, 并且极小概率会沿着入射方向反射回去。这个分布律的表达式有一个 ![[公式]](https://www.zhihu.com/equation?tex=%5Ccos%5E3+%28%5Cphi%29) 的系数, 其中 ![[公式]](https://www.zhihu.com/equation?tex=%5Cphi) 是反射光线距离法向量的夹角。这样当光线从一个离表面很小的角度射入时, 也会散射到一片很大的区域, 对最终颜色值的影响也会更低。

事实上的lambertian的分布律的系数是 ![[公式]](https://www.zhihu.com/equation?tex=%5Ccos+%28%5Cphi%29) ，lambertian散射后的光线距离法相比较近的概率会更高, 分布律会更加均衡。

因为选取的是单位球面上的点。可以通过在单位球内选取一个随机点, 然后将其单位化来获得该点。

```cpp
//vec3.hpp
vec3 random_unit_vector() {
    auto a = random_double(0, 2*pi);
    auto z = random_double(-1, 1);
    auto r = sqrt(1 - z*z);
    return vec3(r*cos(a), r*sin(a), z);
}
```

![img](https://pic4.zhimg.com/80/v2-2a37ce5c0bb8e7e3f51cd3f947ce4f8f_720w.jpg)



在球面上生成一个随机的向量

使用新函数random_unit_vector()替换现存的random_unit_sphere():

```cpp
//main.cc
vec3 ray_color(const ray& r, const hittable& world, int depth) {
    hit_record rec;

    if (depth <= 0)
        return vec3(0,0,0);

    if (world.hit(r, 0.001, infinity, rec)) {
        vec3 target = rec.p + rec.normal + random_unit_vector();
        return 0.5 * ray_color(ray(rec.p, target - rec.p), world, depth-1);
    }

    vec3 unit_direction = unit_vector(r.direction());
    auto t = 0.5*(unit_direction.y() + 1.0);
    return (1.0-t)*vec3(1.0, 1.0, 1.0) + t*vec3(0.5, 0.7, 1.0);
}
```

会得到这样的图片, 和之前很相像:

![img](https://pic3.zhimg.com/80/v2-cc4cdd3d687349fa0d472a701f541cfe_720w.jpg)



正确的lambertian球体

由于场景太简单, 区分这两种方法是比较难的，但应该能够注意到视觉上的一些差异:

1.阴影部分少了
2.大球和小球都变亮了

这些变化都是由散射光线的单位规整化引起的, 现在更少的光线会朝着发现方向散射。对于漫发射的物体来说, 他们会变得更亮。因为更多光线朝着摄像机反射。对于阴影部分来说, 更少的光线朝上反射, 所以小球下方的大球区域会变得更加明亮。



简单来说两种方法都选取了一个随机方向的向量

​	一种是从单位球体内取的, 其长度是随机的

​	另一种是从单位球面上取的, 长度固定为单位向量长度



另一种具有启发性的方法是, 直接从入射点开始选取一个随机的方向, 然后再判断是否在法向量所在的那个半球。

在使用lambertian漫发射模型前, 早期的光线追踪论文中大部分使用的都是这个方法:

```cpp
//vec3.hpp
vec3 random_in_hemisphere(const vec3& normal) {
    vec3 in_unit_sphere = random_in_unit_sphere();
    if (dot(in_unit_sphere, normal) > 0.0) // In the same hemisphere as the normal
        return in_unit_sphere;
    else
        return -in_unit_sphere;
}
```

将的新函数套入ray_color()函数:

```cpp
//vec3.hpp
vec3 ray_color(const ray& r, const hittable& world, int depth) {
    hit_record rec;

    if (depth <= 0)
        return vec3(0,0,0);

    if (world.hit(r, 0.001, infinity, rec)) {
        vec3 target = rec.p + random_in_hemisphere(rec.normal);
        return 0.5 * ray_color(ray(rec.p, target - rec.p), world, depth-1);
    }

    vec3 unit_direction = unit_vector(r.direction());
    auto t = 0.5*(unit_direction.y() + 1.0);
    return (1.0-t)*vec3(1.0, 1.0, 1.0) + t*vec3(0.5, 0.7, 1.0);
}
```

会得到如下的图片:

![img](https://pic4.zhimg.com/80/v2-733c6ea8353952f6f50e33a1990c8a67_720w.jpg)

使用半球面向量渲染漫反射球体



## 材质

如果想让不同的物体能拥有不同的材质, 又面临着一个设计上的抉择。可以设计一个宇宙无敌大材质, 这个材质里面有数不胜数的参数和材质类型可供选择。这样其实也不错, 但还可以设计并封装一个抽象的材质类。反正喜欢后面一种, 对于的程序来说, 一个材质类应该封装两个功能进去:

1.生成散射后的光线(或者说它吸收了入射光线)
2.如果发生散射, 决定光线会变暗多少(attenuate)

下面来看一下这个抽象类:

```cpp
//material.hpp
class material {
    public:
        virtual bool scatter(
            const ray& r_in, const hit_record& rec, vec3& attenuation, ray& scattered
        ) const = 0;
};
```

在函数中使用hit_record作为传入参数, 就可以不用传入一大堆变量了。当然如果想传一堆变量进去的话也行。

物体和材质还要能够联系在一起。在C++中只要告诉编译器, 在`hit_record`里面存了个材质的指针。

```cpp
//hittable.hpp
#ifndef HITTABLE_H
#define HITTABLE_H

#include "rtweekend.hpp"

class material;

struct hit_record {
    vec3 p;
    vec3 normal;
    shared_ptr<material> mat_ptr;
    double t;
    bool front_face;


    inline void set_face_normal(const ray& r, const vec3& outward_normal) {
        front_face = dot(r.direction(), outward_normal) < 0;
        normal = front_face ? outward_normal :-outward_normal;
    }
};

class hittable {
    public:
        virtual bool hit(const ray& r, double t_min, double t_max, hit_record& rec) const = 0;
};

#endif
```

光线会如何与表面交互是由具体的材质所决定的。

hit_record在设计上就是为了把一堆要传的参数给打包在了一起。当光线射入一个表面(比如一个球体), hit_record中的材质指针会被球体的材质指针所赋值, 而球体的材质指针是在main()函数中构造时传入的。当color()函数获取到hit_record时, 可以找到这个材质的指针, 然后由材质的函数来决定光线是否发生散射, 怎么散射。

必须在球体的构造函数和变量区域中加入材质指针, 以便之后传给`hit_record`。

```cpp
class sphere: public hittable {
    public:
        sphere() {}
+        sphere(vec3 cen, double r, shared_ptr<material> m)
+            : center(cen), radius(r), mat_ptr(m) {};

        virtual bool hit(const ray& r, double tmin, double tmax, hit_record& rec) const;

    public:
        vec3 center;
        double radius;
+        shared_ptr<material> mat_ptr;
};

bool sphere::hit(const ray& r, double t_min, double t_max, hit_record& rec) const {
    vec3 oc = r.origin() - center;
    auto a = r.direction().length_squared();
    auto half_b = dot(oc, r.direction());
    auto c = oc.length_squared() - radius*radius;
    auto discriminant = half_b*half_b - a*c;

    if (discriminant > 0) {
        auto root = sqrt(discriminant);
        auto temp = (-half_b - root)/a;
        if (temp < t_max && temp > t_min) {
            rec.t = temp;
            rec.p = r.at(rec.t);
            vec3 outward_normal = (rec.p - center) / radius;
            rec.set_face_normal(r, outward_normal);
            rec.mat_ptr = mat_ptr;
            return true;
        }
        temp = (-half_b + root) / a;
        if (temp < t_max && temp > t_min) {
            rec.t = temp;
            rec.p = r.at(rec.t);
            vec3 outward_normal = (rec.p - center) / radius;                
            rec.set_face_normal(r, outward_normal);
+            rec.mat_ptr = mat_ptr;
            return true;
        }
    }
    return false;
}
```

对于之前写过的Lambertian(漫反射)材质来说, 这里有3种理解方法, 

 1. 光线永远发生散射, 每次散射衰减至R

 2.  光线并不衰减, 转而物体吸收(1-R)的光线

 3. 两种的结合

    

可以写出Lambertian的材质类:

```cpp
//material.hpp
class lambertian : public material {
    public:
        lambertian(const vec3& a) : albedo(a) {}

        virtual bool scatter(
            const ray& r_in, const hit_record& rec, vec3& attenuation, ray& scattered
        ) const {
            vec3 scatter_direction = rec.normal + random_unit_vector();
            scattered = ray(rec.p, scatter_direction);
            attenuation = albedo;
            return true;
        }

    public:
        vec3 albedo;
};
```

也可以让光线根据一定的概率p发生散射并使光线的衰减率(代码中的attenuation)为 ![[公式]](https://www.zhihu.com/equation?tex=albedo%2Fp) 

对于光滑的金属材质来说, 光线是不会像漫反射那样随机散射的, 而是产生反射



![img](https://pic3.zhimg.com/80/v2-8029e9a049925eeec08d6884d44de012_720w.jpg)



反射方向的向量如图所示为 ![[公式]](https://www.zhihu.com/equation?tex=%5Cvec%7BV%7D%2B2%5Cvec%7BB%7D) , 其中规定向量 ![[公式]](https://www.zhihu.com/equation?tex=%5Cvec%7BN%7D) 是单位向量, 但 ![[公式]](https://www.zhihu.com/equation?tex=%5Cvec%7BV%7D) 不一定是。向量B的长度应为 ![[公式]](https://www.zhihu.com/equation?tex=%5Cvec%7BV%7D%5Ccdot%5Cvec%7BN%7D) , 因为向量 ![[公式]](https://www.zhihu.com/equation?tex=%5Cvec%7BV%7D) 与向量 ![[公式]](https://www.zhihu.com/equation?tex=%5Cvec%7BN%7D) 的方向相反, 这里需要再加上一个负号, 于是有:

```cpp
//vec3.hpp
vec3 reflect(const vec3& v, const vec3& n) {
    return v - 2*dot(v,n)*n;
}
```

金属材质使用上面的公式来计算反射方向:

```cpp
//material.h
class metal : public material {
    public:
        metal(const vec3& a) : albedo(a) {}

        virtual bool scatter(
            const ray& r_in, const hit_record& rec, vec3& attenuation, ray& scattered
        ) const {
            vec3 reflected = reflect(unit_vector(r_in.direction()), rec.normal);
            scattered = ray(rec.p, reflected);
            attenuation = albedo;
            return (dot(scattered.direction(), rec.normal) > 0);
        }

    public:
        vec3 albedo;
};
```

修改一下color:

```cpp
//main.cc
vec3 ray_color(const ray& r, const hittable& world, int depth) {
    hit_record rec;

    // If we've exceeded the ray bounce limit, no more light is gathered.
    if (depth <= 0)
        return vec3(0,0,0);

    if (world.hit(r, 0.001, infinity, rec)) {
        ray scattered;
        vec3 attenuation;
        if (rec.mat_ptr->scatter(r, rec, attenuation, scattered))
            return attenuation * ray_color(scattered, world, depth-1);
        return vec3(0,0,0);
    }

    vec3 unit_direction = unit_vector(r.direction());
    auto t = 0.5*(unit_direction.y() + 1.0);
    return (1.0-t)*vec3(1.0, 1.0, 1.0) + t*vec3(0.5, 0.7, 1.0);
}
```

场景加入一些金属球:

```cpp
//main.cc
int main() {
    const int image_width = 200;
    const int image_height = 100;
    const int samples_per_pixel = 100;
    const int max_depth = 50;

    std::cout << "P3\n" << image_width << " " << image_height << "\n255\n";

    hittable_list world;

    world.add(make_shared<sphere>(
        vec3(0,0,-1), 0.5, make_shared<lambertian>(vec3(0.7, 0.3, 0.3))));

    world.add(make_shared<sphere>(
        vec3(0,-100.5,-1), 100, make_shared<lambertian>(vec3(0.8, 0.8, 0.0))));

    world.add(make_shared<sphere>(vec3(1,0,-1), 0.5, make_shared<metal>(vec3(0.8, 0.6, 0.2))));
    world.add(make_shared<sphere>(vec3(-1,0,-1), 0.5, make_shared<metal>(vec3(0.8, 0.8, 0.8))));

    camera cam;
    for (int j = image_height-1; j >= 0; --j) {
        std::cerr << "\rScanlines remaining: " << j << ' ' << std::flush;
        for (int i = 0; i < image_width; ++i) {
            vec3 color(0, 0, 0);
            for (int s = 0; s < samples_per_pixel; ++s) {
                auto u = (i + random_double()) / image_width;
                auto v = (j + random_double()) / image_height;
                ray r = cam.get_ray(u, v);
                color += ray_color(r, world, max_depth);
            }
            color.write_color(std::cout, samples_per_pixel);
        }
    }

    std::cerr << "\nDone.\n";
}
```

得到:

![img](https://pic1.zhimg.com/80/v2-02f0718b49ea2f298ca637e70f1e2ae4_720w.jpg)

金属球



给反射方向加入一点点随机性, 只要在算出反射向量后, 在其终点为球心的球内随机选取一个点作为最终的终点:

当然这个球越大, 金属看上去就更加模糊(fuzzy, 或者说粗糙)。所以这里引入一个变量来表示模糊的程度(fuzziness)(所以当fuzz=0时不会产生模糊)。如果fuzz, 也就是随机球的半径很大, 光线可能会散射到物体内部去。这时候可以认为物体吸收了光线。

```cpp
//material.hpp
class metal : public material {
    public:
        metal(const vec3& a, double f) : albedo(a), fuzz(f < 1 ? f : 1) {}

        virtual bool scatter(
            const ray& r_in, const hit_record& rec, vec3& attenuation, ray& scattered
        ) const {
            vec3 reflected = reflect(unit_vector(r_in.direction()), rec.normal);
            scattered = ray(rec.p, reflected + fuzz*random_in_unit_sphere());
            attenuation = albedo;
            return (dot(scattered.direction(), rec.normal) > 0);//dot<0认为吸收
        }

    public:
        vec3 albedo;
        double fuzz;
};
```

可以将模糊值设置为0.3和1.0, 图片会变成这样:

![img](https://pic4.zhimg.com/80/v2-e362c4148c81abe70ebefa2d089c3eef_720w.jpg)

模糊的金属



## 折射

透明的材料, 例如水, 玻璃, 和钻石都是绝缘体。当光线击中这类材料时, 一条光线会分成两条, 一条发生反射, 一条发生折射。会采取这样的策略: 每次光线与物体相交时, 要么反射要么折射, 一次只发生一种情况,随机选取。反正最后采样次数多, 会给这些结果取个平均值。

折射部分是最难去debug的部分。常常一开始让所有的光线只发生折射来调试。在这个项目中, 加入了两个这样的玻璃球, 并且得到下图(还没教怎么弄出这样的玻璃球, 先往下读, 一会儿就知道了):

![img](https://pic3.zhimg.com/80/v2-e61f5ce6422cdb6fd13190900fc30896_720w.jpg)

这图看上去是对的么? 玻璃球在现实世界中看上去和这差不多。但是, 其实这图不对。玻璃球应该会翻转上下, 也不会有这种奇怪的黑圈。输出了图片中心的一条光线来debug, 发现它完全错了, 调试的时候也可以这样来。

折射法则是由Snell法则定义的:

![[公式]](https://www.zhihu.com/equation?tex=+%5Ceta+%5Ccdot+%5Csin%5Ctheta+%3D+%5Ceta%27+%5Ccdot+%5Csin%5Ctheta%27)

![[公式]](https://www.zhihu.com/equation?tex=%5Ctheta) 与 ![[公式]](https://www.zhihu.com/equation?tex=%5Ctheta%27) 是入射光线与折射光线距离法相的夹角, ![[公式]](https://www.zhihu.com/equation?tex=%5Ceta) 与 ![[公式]](https://www.zhihu.com/equation?tex=%5Ceta%27) (读作eta和eta prime)是介质的折射率(规定空气为1.0, 玻璃为1.3-1.7,钻石为2.4), 如图:

![img](https://pic1.zhimg.com/80/v2-4f08d1085c6bc318a23e135512b6def0_720w.jpg)

折射光线示意图

为了解出折射光线的方向, 需要解出 ![[公式]](https://www.zhihu.com/equation?tex=%5Csin%5Ctheta) :

![[公式]](https://www.zhihu.com/equation?tex=+%5Csin%5Ctheta%27+%3D+%5Cfrac%7B%5Ceta%7D%7B%5Ceta%27%7D+%5Ccdot+%5Csin%5Ctheta+)

在折射介质部分有射线光线 ![[公式]](https://www.zhihu.com/equation?tex=%5Cmathbf%7BR%27%7D) 与法向量 ![[公式]](https://www.zhihu.com/equation?tex=%5Cmathbf%7BN%27%7D) , 它们的夹角为 ![[公式]](https://www.zhihu.com/equation?tex=%5Ctheta%27) 。可以把光线 ![[公式]](https://www.zhihu.com/equation?tex=%5Cmathbf%7BR%27%7D) 分解成垂直和水平与法向量 ![[公式]](https://www.zhihu.com/equation?tex=%5Cmathbf%7BN%27%7D) 的两个向量:

![[公式]](https://www.zhihu.com/equation?tex=%5Cmathbf%7BR%27%7D+%3D+%5Cmathbf%7BR%27%7D_%7B%5Cparallel%7D+%2B+%5Cmathbf%7BR%27%7D_%7B%5Cbot%7D+)

有:

![[公式]](https://www.zhihu.com/equation?tex=%5Cmathbf%7BR%27%7D_%7B%5Cparallel%7D+%3D+%5Cfrac%7B%5Ceta%7D%7B%5Ceta%27%7D+%28%5Cmathbf%7BR%7D+%2B+%5Ccos%5Ctheta+%5Cmathbf%7BN%7D%29)

![[公式]](https://www.zhihu.com/equation?tex=%5Cmathbf%7BR%27%7D_%7B%5Cbot%7D+%3D+-%5Csqrt%7B1+-+%7C%5Cmathbf%7BR%27%7D_%7B%5Cparallel%7D%7C%5E2%7D+%5Cmathbf%7BN%7D+)

再解 ![[公式]](https://www.zhihu.com/equation?tex=%5Ccos%5Ctheta) , 

![[公式]](https://www.zhihu.com/equation?tex=%5Cmathbf%7BA%7D+%5Ccdot+%5Cmathbf%7BB%7D+%3D+%7C%5Cmathbf%7BA%7D%7C+%7C%5Cmathbf%7BB%7D%7C+%5Ccos%5Ctheta)

将 ![[公式]](https://www.zhihu.com/equation?tex=%5Cmathbf%7BA%7D) 与 ![[公式]](https://www.zhihu.com/equation?tex=%5Cmathbf%7BB%7D) 归一化为单位向量:

![[公式]](https://www.zhihu.com/equation?tex=%5Cmathbf%7BA%7D+%5Ccdot+%5Cmathbf%7BB%7D+%3D+%5Ccos%5Ctheta+)

表达垂直的那个向量:

![[公式]](https://www.zhihu.com/equation?tex=%5Cmathbf%7BR%27%7D_%7B%5Cparallel%7D+%3D+%5Cfrac%7B%5Ceta%7D%7B%5Ceta%27%7D+%28%5Cmathbf%7BR%7D+%2B+%28%5Cmathbf%7B-R%7D+%5Ccdot+%5Cmathbf%7BN%7D%29+%5Cmathbf%7BN%7D%29)

根据上述公式, 就能写出计算折射光线 ![[公式]](https://www.zhihu.com/equation?tex=%5Cmathbf%7BR%27%7D) 的函数:

```cpp
//vec3.hpp
vec3 refract(const vec3& uv, const vec3& n, double etai_over_etat) {
    auto cos_theta = dot(-uv, n);
    vec3 r_out_parallel =  etai_over_etat * (uv + cos_theta*n);
    vec3 r_out_perp = -sqrt(1.0 - r_out_parallel.length_squared()) * n;
    return r_out_parallel + r_out_perp;
}
```

发生折射的绝缘体材质为:

```cpp
//material.hpp
class dielectric : public material {
    public:
        dielectric(double ri) : ref_idx(ri) {}

        virtual bool scatter(
            const ray& r_in, const hit_record& rec, vec3& attenuation, ray& scattered
        ) const {
            attenuation = vec3(1.0, 1.0, 1.0);
            double etai_over_etat;
            if (rec.front_face) {
                etai_over_etat = 1.0 / ref_idx;
            } else {
                etai_over_etat = ref_idx;
            }

            vec3 unit_direction = unit_vector(r_in.direction());
            vec3 refracted = refract(unit_direction, rec.normal, etai_over_etat);
            scattered = ray(rec.p, refracted);
            return true;
        }

        double ref_idx;
};
```

![img](https://pic4.zhimg.com/80/v2-368925cc729af1d7d580d5b92fc8f2e3_720w.jpg)

只发生折射的玻璃材质



当光线从高折射律介质射入低折射率介质时, 对于上述的Snell方程可能没有实解【 ![[公式]](https://www.zhihu.com/equation?tex=%5Csin%5Ctheta%3E1) 】。这时候就不会发生折射, 所以就会出现许多小黑点。回头看一下snell法则的式子:

![[公式]](https://www.zhihu.com/equation?tex=%5Csin%5Ctheta%27+%3D+%5Cfrac%7B%5Ceta%7D%7B%5Ceta%27%7D+%5Ccdot+%5Csin%5Ctheta)

如果光线从玻璃( ![[公式]](https://www.zhihu.com/equation?tex=%5Ceta+%3D+1.5) )射入空气( ![[公式]](https://www.zhihu.com/equation?tex=%5Ceta+%3D+1.0) )

![[公式]](https://www.zhihu.com/equation?tex=%5Csin%5Ctheta%27+%3D+%5Cfrac%7B1.5%7D%7B1.0%7D+%5Ccdot+%5Csin%5Ctheta)

又因为 ![[公式]](https://www.zhihu.com/equation?tex=%5Csin%5Ctheta%27) 是不可能比1大的,一旦这种情况发生:

![[公式]](https://www.zhihu.com/equation?tex=+%5Cfrac%7B1.5%7D%7B1.0%7D+%5Ccdot+%5Csin%5Ctheta+%3E+1.0+)

方程无解，认为光线无法发生折射的时候, 发生了反射:

```cpp
//material.hpp
if(etai_over_etat * sin_theta > 1.0) {
    // Must Reflect
    ...
}
else {
    // Can Refract
    ...
}
```



常常在实心物体的内部发生, 所以称这种情况被称为"全内反射"。这也当浸入水中时, 发现水与空气的交界处看上去像一面镜子的原因。

可以用三角函数解出sin_theta

![[公式]](https://www.zhihu.com/equation?tex=%5Csin%5Ctheta+%3D+%5Csqrt%7B1+-+%5Ccos%5E2%5Ctheta%7D)

其中的cos_theta为

![[公式]](https://www.zhihu.com/equation?tex=%5Ccos%5Ctheta+%3D+%5Cmathbf%7BR%7D+%5Ccdot+%5Cmathbf%7BN%7D)

```cpp
//material.hpp
double cos_theta = ffmin(dot(-unit_direction, rec.normal), 1.0);
double sin_theta = sqrt(1.0 - cos_theta*cos_theta);
if(etai_over_etat * sin_theta > 1.0) {
    // Must Reflect
    ...
}
else {
    // Can Refract
    ...
}
```

一个在可以偏折的情况下总是偏折, 其余情况发生反射的绝缘体材质为:

```cpp
//material.hpp
class dielectric : public material {
    public:
        dielectric(double ri) : ref_idx(ri) {}

        virtual bool scatter(
            const ray& r_in, const hit_record& rec, vec3& attenuation, ray& scattered
        ) const {
            attenuation = vec3(1.0, 1.0, 1.0);
            double etai_over_etat = (rec.front_face) ? (1.0 / ref_idx) : (ref_idx);

            vec3 unit_direction = unit_vector(r_in.direction());
            double cos_theta = ffmin(dot(-unit_direction, rec.normal), 1.0);
            double sin_theta = sqrt(1.0 - cos_theta*cos_theta);
            if (etai_over_etat * sin_theta > 1.0 ) {
                vec3 reflected = reflect(unit_direction, rec.normal);
                scattered = ray(rec.p, reflected);
                return true;
            }

            vec3 refracted = refract(unit_direction, rec.normal, etai_over_etat);
            scattered = ray(rec.p, refracted);
            return true;
        }

    public:
        double ref_idx;
};
```

这里的光线衰减率为1——就是不衰减, 玻璃表面不吸收光的能量。使用下面的参数:

```cpp
//main.cc
world.add(make_shared<sphere>(
    vec3(0,0,-1), 0.5, make_shared<lambertian>(vec3(0.1, 0.2, 0.5))));

world.add(make_shared<sphere>(
    vec3(0,-100.5,-1), 100, make_shared<lambertian>(vec3(0.8, 0.8, 0.0))));

world.add(make_shared<sphere>(vec3(1,0,-1), 0.5, make_shared<metal>(vec3(0.8, 0.6, 0.2), 0.0)));
world.add(make_shared<sphere>(vec3(-1,0,-1), 0.5, make_shared<dielectric>(1.5)));
```

得到:

![img](https://pic2.zhimg.com/80/v2-7484d6b5ddb4744f741747ac725e7749_720w.jpg)

现实世界中的玻璃, 发生折射的概率会随着入射角而改变——从一个很狭窄的角度去看玻璃窗, 它会变成一面镜子。

用由Christophe Schlick提出的式子化简:

```cpp
double schlick(double cosine, double ref_idx) {
    auto r0 = (1-ref_idx) / (1+ref_idx);
    r0 = r0*r0;
    return r0 + (1-r0)*pow((1 - cosine),5);
}
```

完整的材质代码:

```cpp
//material.hpp
class dielectric : public material {
    public:
        dielectric(double ri) : ref_idx(ri) {}

        virtual bool scatter(
            const ray& r_in, const hit_record& rec, vec3& attenuation, ray& scattered
        ) const {
            attenuation = vec3(1.0, 1.0, 1.0);
            double etai_over_etat = (rec.front_face) ? (1.0 / ref_idx) : (ref_idx);

            vec3 unit_direction = unit_vector(r_in.direction());
            double cos_theta = ffmin(dot(-unit_direction, rec.normal), 1.0);
            double sin_theta = sqrt(1.0 - cos_theta*cos_theta);
            if (etai_over_etat * sin_theta > 1.0 ) {
                vec3 reflected = reflect(unit_direction, rec.normal);
                scattered = ray(rec.p, reflected);
                return true;
            }
            double reflect_prob = schlick(cos_theta, etai_over_etat);
            if (random_double() < reflect_prob)
            {
                vec3 reflected = reflect(unit_direction, rec.normal);
                scattered = ray(rec.p, reflected);
                return true;
            }
            vec3 refracted = refract(unit_direction, rec.normal, etai_over_etat);
            scattered = ray(rec.p, refracted);
            return true;
        }

    public:
        double ref_idx;
};
```



将球的半径设为负值, 形状看上去并没变化, 法相翻转到内部，用这个特性来做出一个通透的玻璃球:

```cpp
world.add(make_shared<sphere>(vec3(0,0,-1), 0.5, make_shared<lambertian>(vec3(0.1, 0.2, 0.5))));
world.add(make_shared<sphere>(
    vec3(0,-100.5,-1), 100, make_shared<lambertian>(vec3(0.8, 0.8, 0.0))));
world.add(make_shared<sphere>(vec3(1,0,-1), 0.5, make_shared<metal>(vec3(0.8, 0.6, 0.2), 0.3)));
world.add(make_shared<sphere>(vec3(-1,0,-1), 0.5, make_shared<dielectric>(1.5)));
world.add(make_shared<sphere>(vec3(-1,0,-1), -0.45, make_shared<dielectric>(1.5)));
```



![img](https://pic1.zhimg.com/80/v2-1c249bd477ba8b2aeb5f49d845486fbc_720w.jpg)

一个通透的玻璃球



## 可自定义位置的摄像机

摄像机总是和绝缘体一样难以debug，所以要总是一步步搭建的摄像机类。

首先, 使摄像机能调整其视野范围(field of view, fov)。

​	fov是的视角。

图片不是方的, 所以垂直和水平的fov值是不同的。总是使用垂直方向的fov，并且总是使用角度制来传参, 在构造函数中再将其转化为弧度



![img](https://pic1.zhimg.com/80/v2-d8f93c297ec72b9e072535b4f6dcb9c8_720w.jpg)

摄像机示意图

![[公式]](https://www.zhihu.com/equation?tex=h+%3D+%5Ctan%28%5Cfrac%7B%5Ctheta%7D%7B2%7D%29) 的摄像机类:

```cpp
//camera.hpp
class camera {
    public:
        camera(
            double vfov, // top to bottom, in degrees
            double aspect
        ) {
            origin = vec3(0.0, 0.0, 0.0);

            auto theta = degrees_to_radians(vfov);
            auto half_height = tan(theta/2);
            auto half_width = aspect * half_height;

            lower_left_corner = vec3(-half_width, -half_height, -1.0);

            horizontal = vec3(2*half_width, 0.0, 0.0);
            vertical = vec3(0.0, 2*half_height, 0.0);
        }

        ray get_ray(double u, double v) {
            return ray(origin, lower_left_corner + u*horizontal + v*vertical - origin);
        }

    public:
        vec3 origin;
        vec3 lower_left_corner;
        vec3 horizontal;
        vec3 vertical;
};
```

使用一个cam(90, double(image_width)/image_height)的摄像机 拍下面的球:

```cpp
auto R = cos(pi/4);
hittable_list world;
world.add(make_shared<sphere>(vec3(-R,0,-1), R, make_shared<lambertian>(vec3(0, 0, 1))));
world.add(make_shared<sphere>(vec3( R,0,-1), R, make_shared<lambertian>(vec3(1, 0, 0))));
```

得到:

![img](https://pic4.zhimg.com/80/v2-a5ee444b1b746b41d932c4fc403d6ee3_720w.jpg)

一个90°的广角镜头





摄像机所在的位置叫 lookfrom, 看向的点叫lookat，需要一个变量去描述摄像机的倾斜程度, 或者说摄像机绕着轴lookfrom - lookat旋转的角度

可以使用任意的方向向量, 将其投影到上图的平面中来获得摄像机的up vector（vup向量）。

经过一系列的点乘操作, 会有完整的u,v,w三个向量来描述摄像机的旋向vup。vup, v, w处于同一平面内。和先前的摄像机面对着-Z方向一样, 修改后的任意视角摄像机面对着-w方向。



使用世界坐标系的上方向向量(0,1,0)指定vup

```cpp
class camera {
    public:
        camera(
            vec3 lookfrom, vec3 lookat, vec3 vup,
            double vfov, // top to bottom, in degrees
            double aspect
        ) {
            origin = lookfrom;
            vec3 u, v, w;

            auto theta = degrees_to_radians(vfov);
            auto half_height = tan(theta/2);
            auto half_width = aspect * half_height;
            w = unit_vector(lookfrom - lookat);
            u = unit_vector(cross(vup, w));
            v = cross(w, u);

            lower_left_corner = origin - half_width*u - half_height*v - w;

            horizontal = 2*half_width*u;
            vertical = 2*half_height*v;
        }

        ray get_ray(double s, double t) {
            return ray(origin, lower_left_corner + s*horizontal + t*vertical - origin);
        }

    public:
        vec3 origin;
        vec3 lower_left_corner;
        vec3 horizontal;
        vec3 vertical;
};
```

可以改变的视角:

```cpp
//main.cc
const auto aspect_ratio = double(image_width) / image_height;
...
camera cam(vec3(-2,2,1), vec3(0,0,-1), vup, 90, aspect_ratio);
```

会得到:



![img](https://pic4.zhimg.com/80/v2-3a76513e86c48fc03887ee41abb03b53_720w.jpg)

从远处看

改变一下fov:

![img](https://pic1.zhimg.com/80/v2-3b1f9d7c41cfdf7c32c4452f82e80d48_720w.jpg)

放大看

## 散焦模糊

现实世界中的摄像机产生对焦模糊的原因是因为他们需要一个很大的孔, 而不是一个针眼大小的小孔来聚集光线。这会导致所有的东西都被散焦了。但如果在孔内加入一块透镜, 在一段距离内的所有物体都会被对焦。可以这样来想象透镜:所有的光线从同一点分散射出, 击中透镜后又聚焦在图像传感器上的一个点上。

​	现实世界中的摄像机的透镜组是很复杂的。但对于写代码来说, 只需要模拟上述的顺序: 图像传感器, 透镜, 快门, 然后射出光线, 最后翻转图片(进过透镜成像会被上下翻转)。只要从一个虚拟的透镜范围中发射光线到的摄像机平面就能模拟了,这个透镜与平面的距离成为焦(focus_dist)

之前所有的光线都是从lookfrom发出的, 但现在加入了散焦模糊, 所有光线都从内部的一个虚拟透镜发出, 经过lookfrom点, 这个透镜的半径越大, 图像就越模糊。可以认为之前的摄像机, 这个半径为0。

```cpp
//vec3.hpp 从一个单位小圆盘射出光线
vec3 random_in_unit_disk() {
    while (true) {
        auto p = vec3(random_double(-1,1), random_double(-1,1), 0);
        if (p.length_squared() >= 1) continue;
        return p;
    }
}
```

完整的camera类

```cpp
class camera {
    public:
        camera(
            vec3 lookfrom, vec3 lookat, vec3 vup,
            double vfov, // top to bottom, in degrees
            double aspect, double aperture, double focus_dist
        ) {
            origin = lookfrom;
            lens_radius = aperture / 2;

            auto theta = degrees_to_radians(vfov);
            auto half_height = tan(theta/2);
            auto half_width = aspect * half_height;

            w = unit_vector(lookfrom - lookat);
            u = unit_vector(cross(vup, w));
            v = cross(w, u);
            lower_left_corner = origin
                              - half_width * focus_dist * u
                              - half_height * focus_dist * v
                              - focus_dist * w;

            horizontal = 2*half_width*focus_dist*u;
            vertical = 2*half_height*focus_dist*v;
        }

        ray get_ray(double s, double t) {
            vec3 rd = lens_radius * random_in_unit_disk();
            vec3 offset = u * rd.x() + v * rd.y();

            return ray(
                origin + offset,
                lower_left_corner + s*horizontal + t*vertical - origin - offset
           );
        }

    public:
        vec3 origin;
        vec3 lower_left_corner;
        vec3 horizontal;
        vec3 vertical;
        vec3 u, v, w;
        double lens_radius;
};
```

使用一个快门光圈:

```cpp
//main.cc
const auto aspect_ratio = double(image_width) / image_height;
...
vec3 lookfrom(3,3,2);
vec3 lookat(0,0,-1);
vec3 vup(0,1,0);
auto dist_to_focus = (lookfrom-lookat).length();
auto aperture = 2.0;

camera cam(lookfrom, lookat, vup, 20, aspect_ratio, aperture, dist_to_focus);
```

就有:

![img](https://pic3.zhimg.com/80/v2-8d5daba641a0c74f8af90cd9b47f9e42_720w.jpg)

加入景深效果

## 

总结：

```cpp
//main.cc
hittable_list random_scene() {
    hittable_list world;

    world.add(make_shared<sphere>(
        vec3(0,-1000,0), 1000, make_shared<lambertian>(vec3(0.5, 0.5, 0.5))));

    int i = 1;
    for (int a = -11; a < 11; a++) {
        for (int b = -11; b < 11; b++) {
            auto choose_mat = random_double();
            vec3 center(a + 0.9*random_double(), 0.2, b + 0.9*random_double());
            if ((center - vec3(4, 0.2, 0)).length() > 0.9) {
                if (choose_mat < 0.8) {
                    // diffuse
                    auto albedo = vec3::random() * vec3::random();
                    world.add(
                        make_shared<sphere>(center, 0.2, make_shared<lambertian>(albedo)));
                } else if (choose_mat < 0.95) {
                    // metal
                    auto albedo = vec3::random(.5, 1);
                    auto fuzz = random_double(0, .5);
                    world.add(
                        make_shared<sphere>(center, 0.2, make_shared<metal>(albedo, fuzz)));
                } else {
                    // glass
                    world.add(make_shared<sphere>(center, 0.2, make_shared<dielectric>(1.5)));
                }
            }
        }
    }

    world.add(make_shared<sphere>(vec3(0, 1, 0), 1.0, make_shared<dielectric>(1.5)));

    world.add(
        make_shared<sphere>(vec3(-4, 1, 0), 1.0, make_shared<lambertian>(vec3(0.4, 0.2, 0.1))));

    world.add(
        make_shared<sphere>(vec3(4, 1, 0), 1.0, make_shared<metal>(vec3(0.7, 0.6, 0.5), 0.0)));

    return world;
}

int main() {
    ...
    auto world = random_scene();

    vec3 lookfrom(13,2,3);
    vec3 lookat(0,0,0);
    vec3 vup(0,1,0);
    auto dist_to_focus = 10.0;
    auto aperture = 0.1;

    camera cam(lookfrom, lookat, vup, 20, aspect_ratio, aperture, dist_to_focus);
    ...
}
```

得到：

![img](https://pic2.zhimg.com/80/v2-7483e528431ca10622ddd31ce8ebbba9_720w.jpg)

最终场景

可能会发现玻璃球没有阴影, 使得他们看上去像漂浮在空中。

这不是bug(在现实世界中很少有机会见到真正的玻璃球, 它们看起来的确就是这样的)。



# 下阶段目标：

1.光照。可以使用阴影光线来显式实现这部分, 也可以使用产生光线的材质来隐式实现。

2.偏移散射光线, 然后降低这些光线的权重来消除偏移。

3.加入三角形。大部分模型都是三角网格。

4.表面纹理。贴墙纸一样把图片贴到物体上去。

5.固体纹理。

6.体积体(volumes 即雾等）与其他介质。

7.并行优化。