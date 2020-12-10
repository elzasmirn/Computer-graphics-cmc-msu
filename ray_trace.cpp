#include <iostream>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <vector>
#include "geo_figs.h"

using namespace std;

#define INF 1000000

bool refract(vec3& ray_dir, vec3 nHit, double ior) 
{
    double eta = 1.0/ior; // eta = in_IOR/out_IOR,  IOR-index of refraction
    double cosi = -nHit*ray_dir;
    if(cosi < 0)
    {
        cosi *= -1;
        nHit = -nHit;
        eta = 1.0/eta;
    }
    double k = 1.0 - eta*eta*(1.0-cosi*cosi);
    if(k >= 0)
        ray_dir = ( ray_dir * eta +  nHit * (cosi*eta - sqrt(k)) ).normalize();
    return (k > 0);
}

#define MAX_DEPTH 3
#define FIGNUM_MAX 8
#define FIGNUM 8
//#define FIGNUM 3
#define LIGHTNUM 3

Color ray_trace (Ray& ray, Figure** figures, Light** lights, const int depth)
{
    double t0, t1, tnear = INF;
    bool isInter = false; //intersection flag
    bool isShadow = false;
    Ray sRay; //shadow ray
    Color resultCol;
    int nearInd = -1;
    for (int i = 0; i < FIGNUM; ++i) 
    {
        t0 = t1 = INF;
        if ( figures[i]->intersect(ray, t0 , t1) )
        {
            if (t0 < 0) t0 = t1;
            if (t0 < tnear) 
            {
                isInter = true;
                tnear = t0;
                nearInd = i;
            }
        }
    }
    
    if ( !isInter ) return resultCol; //no intersection

        double eps = 1e-4;
        vec3 pHit = ray.orig + ray.dir * tnear; //intersect point 
        vec3 nHit; //normal at pHit
        figures[nearInd]->getNorm (pHit,nHit);
        if ( ray.dir*nHit > 0)
            nHit = -nHit;
        //reflection or refraction
        if ( (figures[nearInd]->reflect > 0 || figures[nearInd]->transparent > 0) && depth < MAX_DEPTH)
        {
            Color reflCol, refrCol;
        //fresnel effect
            double facingratio = -ray.dir*nHit;
            double fresnel = mix(pow(1 - facingratio, 3), 1, 0.1);
        //reflection
            vec3 refldir = ray.dir - nHit * 2. * (ray.dir*nHit);
            refldir.normalize();
            Ray reflRay = Ray(pHit + nHit * eps, refldir);
            reflCol = ray_trace(reflRay, figures,lights, depth+1) * figures[nearInd]->reflect; 
        //refraction   
            if (figures[nearInd]->transparent)
                if ( refract(ray.dir,nHit,2) )
//                    refrCol = ray_trace(ray, figures, depth + 1) * figures[nearInd]->transparent;
                    refrCol = ray_trace(ray, figures,lights, depth + 1) * figures[nearInd]->transparent;
            resultCol = reflCol*fresnel + refrCol*(1 - fresnel);
        }
        //shadow rays
        for(int k = 0; k < LIGHTNUM; ++k)
        {
            vec3 lightdir = lights[k]->PT - pHit;
            lightdir.normalize();
            sRay = Ray (pHit + nHit*eps, lightdir); 
            isShadow = false;
            for(int i = 0; i < FIGNUM; ++i)
                if ( figures[i]->intersect(sRay, t0, t1)){
                     isShadow = true; 
                     break;
                }
            if( !isShadow ){
                double dif = nHit*( (lights[k]->PT - pHit).normalize() ); 
                resultCol = resultCol + figures[nearInd]->color * max(0, dif)*lights[k]->intens;
            }
        }
    return resultCol;
}

          

