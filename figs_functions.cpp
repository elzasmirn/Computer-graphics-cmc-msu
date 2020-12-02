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

bool QuadEq(const double &a, const double &b, const double &c, double &x0, double &x1) 
{ 
    if (a == 0) cerr << "QuadEq :: not equation" << endl;
    double discr = b * b - 4 * a * c; 
    if (discr < 0) return false; 
    else if (discr == 0) { 
        x0 = x1 = - 0.5 * b / a; 
    } 
    else { 
        x0 = 0.5 * (-b - sqrt(discr)) / a; 
        x1 = 0.5 * (-b + sqrt(discr)) / a; 
    } 
    return true; 
} 

bool Box::intersect (const Ray& r, double& t0, double& t1) const 
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

bool Plane::intersect (const Ray& r, double& t0, double& t1) const 
    { //Ax+By+Cz+D=0
        double fr =  N.x_*r.orig.x_ + N.y_*r.orig.y_ + N.z_*r.orig.z_ + D;
        double t = N.x_*r.dir.x_ + N.y_*r.dir.y_ + N.z_*r.dir.z_ ;
		if (t==0) {
			cerr << "Plane::intersect: ray parallel?" << endl;
			return false;
		}
        t0 =  -fr / t;
        if (t0 < 0) return false;
        else return true;
    }       
void Plane::getNorm (const vec3& pHit, vec3& nHit) const
    {
        nHit = N;
    }

bool Cylindr::intersect (const Ray& r, double& t0, double& t1) const 
    { //(x-a)^ + (z-b)^ = r^2
        //top and bot intesection
        double t;
        if ( r.dir.y_>0. && top.intersect(r,t,t) ||  r.dir.y_<0. && bot.intersect(r,t,t) )
        {
            vec3 hit = r.orig + r.dir * t - center;
            if ( hit.x_*hit.x_ + hit.z_*hit.z_   < radius*radius )
            {
                t0 = t;
                return true;
            }
        }
        
        //lateral intersection
		vec3 rc_proj=r.orig-center;
		rc_proj.y_=0.;
		vec3 rd_proj=r.dir;
		rd_proj.y_=0.;
		double a=rd_proj.length2();
		double b=2.*(rc_proj*rd_proj);
		double c=rc_proj.length2()-radius * radius;
		if(b>=0. || c<0.) return false;
        if ( !QuadEq(a, b, c, t0, t1) ) return false;
		vec3 hit = r.orig + r.dir * t0;
        if ( fabs (hit.y_ - center.y_) > height/2  ) return false;
		return true;
    }
void Cylindr::getNorm (const vec3& pHit, vec3& nHit) const
    {
        nHit = pHit - center;
        if ( nHit.y_ >= height/2 )
			nHit = vec3(0,1,0);
		else if (nHit.y_ <= -height/2)
			nHit = vec3(0,-1,0);
        else{
            nHit.y_  = 0;
            nHit.normalize();
        }
    }

bool Cone::intersect (const Ray& r, double& t0, double& t1) const 
    { //x^2 + z^2 - y^2 = 0;
        double t;
        if ( bot.intersect(r,t,t) )
        {
            vec3 hit = r.orig + r.dir * t - center;
            if ( hit.x_*hit.x_ + hit.z_*hit.z_   < radius*radius )
            {
                t0 = t;
                return true;
            }
        }
        double rx = r.orig.x_ - center.x_;
        double rz = r.orig.z_ - center.z_;
        double ry = r.orig.y_ - center.y_;
        double a = r.dir.x_*r.dir.x_ + r.dir.z_*r.dir.z_ - r.dir.y_*r.dir.y_;
        double b = 2*( rx*r.dir.x_ + rz*r.dir.z_ - ry*r.dir.y_);
        double c = rx*rx + rz*rz - ry*ry;
        if ( !QuadEq(a, b, c, t0, t1) ) return false;
        if (t0 < 0) t0 = t1;
        if (t0 > 0)
        {
            vec3 hit;
            if (t0 < t1) hit = r.orig + r.dir * t0;
            else hit = r.orig + r.dir * t1;
            if ( fabs (hit.y_ - center.y_) > height || hit.y_ < center.y_  ) return false;
            return true;
        }
        return false;
    }
void Cone::getNorm (const vec3& pHit, vec3& nHit) const 
    {
        nHit = pHit - center;
        if ( fabs (nHit.y_) >= height )
            nHit = vec3(0,-1,0);
        else{
            nHit.y_  = -(nHit.x_*nHit.x_+nHit.z_*nHit.z_)/nHit.y_;
            nHit.normalize();
        }
    }

bool Sphere::intersect (const Ray& r, double& t0, double& t1) const 
    {   
        vec3 v = center - r.orig;
        double tca = v*r.dir;
        if (tca < 0) return false;
        double d2 = v.length2() - tca*tca;
        if(d2 > radius2) return false;
        double thc = sqrt ( radius2-d2 );
        t0 = tca - thc;
        t1 = tca + thc;

        return true;
    }
void Sphere::getNorm(const vec3& pHit, vec3& nHit) const
    {
        nHit = pHit - center;
        nHit.normalize();
    }

bool Bump_Sphere::intersect (const Ray& r, double& t0, double& t1) const 
    {   
        vec3 v = center - r.orig;
        double tca = v*r.dir;
        if (tca < 0) return false;
        double d2 = v.length2() - tca*tca;
        if(d2 > radius2) return false;
        double thc = sqrt ( radius2-d2 );
        t0 = tca - thc;
        t1 = tca + thc;
	
        return true;
    }
