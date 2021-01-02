/*
 * @Author: geekli
 * @Date: 2021-01-02 07:42:25
 * @LastEditTime: 2021-01-02 07:50:25
 * @LastEditors: your name
 * @Description: 
 * @FilePath: /ray_tracing/ray_tracing_in_one_week/ppm_test.cc
 */
#include <iostream>

using namespace std;


int main()
{
    int nx =200;
    int ny=100;
    cout<<"P3\n"<<nx<<" "<<ny<<"\n255\n";
    for(int j=ny-1;j>=0;j--)
    {
        for(int i=0;i<nx;i++)
        {
            float r=float(i)/float(nx);
            float g=float(j)/float(ny);
            float b=0.2;

            int ir=int(255.99*r);
            int ig=int(255.99*g);
            int ib=int(255.99*b);
            cout<<ir<<" "<<ig<<" "<<ib<<"\n";
        }
    }
}