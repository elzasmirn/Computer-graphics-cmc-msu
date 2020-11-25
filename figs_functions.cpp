#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include "geo_figs.h"

using namespace std;

double max(const double a, const double b)
{
    return a > b ? a : b;
}

double min(const double a, const double b)
{
    return a < b ? a : b;
}

double mix(const double &a, const double &b, const double &mix)
{
        return b * mix + a * (1 - mix);
}

bool Box:: intersect (const Ray& r, double& t0, double& t1) const 
    {
        vec3 inv_dir = 1/r.dir;
        double lo = inv_dir.x_*(boxMin.x_ - r.orig.x_);
        double hi = inv_dir.x_*(boxMax.x_ - r.orig.x_);
        t0  = min(lo, hi);
        t1 = max(lo, hi);
 
        double lo1 = inv_dir.y_*(boxMin.y_ - r.orig.y_);
        double hi1 = inv_dir.y_*(boxMax.y_ - r.orig.y_);
        t0 = max(t0, min(lo1, hi1));
        t1 = min(t1, max(lo1, hi1));
 
        double lo2 = inv_dir.z_*(boxMin.z_ - r.orig.z_);
        double hi2 = inv_dir.z_*(boxMax.z_ - r.orig.z_);
        t0  = max(t0, min(lo2, hi2));
        t1 = min(t1, max(lo2, hi2));
 
        return (t0 <= t1) && (t1 > 0);
    }

void Box::getNorm (const vec3& pHit, vec3& nHit) const
    {
        nHit = vec3(0);
        double eps = 1e-7;
        if (      fabs( pHit.x_ - boxMax.x_ ) < eps ) //right
            nHit.x_ = 1;
        else if ( fabs( pHit.x_ - boxMin.x_ ) < eps ) //left
            nHit.x_ = -1;
        else if ( fabs( pHit.y_ - boxMax.y_ ) < eps ) //down
            nHit.y_ = 1;
        else if ( fabs( pHit.y_ - boxMin.y_ ) < eps ) //up
            nHit.y_ = -1;
        else if ( fabs( pHit.z_ - boxMax.z_ ) < eps ) //back
            nHit.z_ = -1;
        else if ( fabs( pHit.z_ - boxMin.z_ ) < eps ) //front
            nHit.z_ = 1;
		else { cerr << "Box::getNorm: not at surface?" <<pHit<< endl; } 
    }