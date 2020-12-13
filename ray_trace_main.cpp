#include <iostream>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <vector>
#include "geo_figs.h"

using namespace std;

#define INF 1000000

Ray build_ray(double x, double y, double w, double h)
{
    double fov = 3.141592654/(2.0);
    vec3 ray_dir;
    ray_dir.x_ = x+0.5 - (w/2.0);
    ray_dir.y_ = y+0.5 - (h/2.0);
    ray_dir.z_ = -(w)/tan(fov/2.0);
    return Ray( vec3(0) , ray_dir.normalize() );
}

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
    double t0,t1,tnear = INF;
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

int main(int argc, char* argv[])
{
    int width,height;
	if (argc > 1)
    {
        if ( *argv[1] == 'c')
            width = 500, height = 500;
        else if ( *argv[1] == 'f')
            width = 1000, height = 1000;
        else if ( *argv[1] == 'g')
            width = 2000, height = 2000;
        else if ( *argv[1] == 'b')
            width = 3000, height = 3000;
        else width = 1000, height = 1000;
    }
    else width = 1000, height = 1000;

    Color pixcolor;
    Light *lights[ LIGHTNUM ]; //array of lights
    Light L1 = Light( vec3 ( 250, -500, -30), 0.4);
    Light L2 = Light( vec3 (-300, -360, -40), 0.6);
    Light L3 = Light( vec3 (30, -200, -5), 0.2);
    lights[0] = &L1;
    lights[1] = &L2;
    lights[2] = &L3;
    Figure *figures[ FIGNUM_MAX ]; //array of objects
 //   Sphere sp1 = Sphere ( vec3(-150, -45, -400), 45, Color(0,120,150), 0.9, 0 ) ;
    Sphere sp2 = Sphere ( vec3(-7, 25, -125),  15,  Color(255,51,153), 0, 0.4 ) ;
 //   Sphere sp3 = Sphere ( vec3(26, -90, -300),  30,  Color(250,250,250), 0.3, 0 ) ;
 //   Cylindr cl1 = Cylindr ( vec3(-25,18,-145), 20, 12, Color(151,48,48), 0.5 , 0 ) ;
    Cylindr cl1 = Cylindr ( vec3(35,32,-160), 10, 22, Color(153,204,255), 0.5 , 0 ) ;
 //   Cylindr cl2 = Cylindr ( vec3(-20,33,-150), 10, 30, Color(135,58,203), 0.3 , 0 ) ;
    Cone co1 = Cone ( vec3(-20,0,-150), 20, 12, Color(48,151,48), 0 , 0 ) ;
    Plane  pl1 = Plane ( vec3(0, -1, 0), 38, Color(47, 47, 47), 0.4, 0) ;
 //   Plane  pl2 = Plane ( vec3(0, 0, -1), -400, Color(30,30,27), 2, 0) ;
    Plane  pl2 = Plane ( vec3(0, 0, -1), -400, Color(47, 47, 47), 5, 0) ;
    Box bx1 = Box ( vec3(-22,38,-140), vec3(-57,3,-190),  Color(153,102,255), 0.1, 0 ) ;

	// points for Rot_B_spline
	int n_pn_spl=5;
	vector <double> y_ps(n_pn_spl), r_ps(n_pn_spl);
	//y_ps[0]=8.-60; r_ps[0]=5.;
	//y_ps[1]=18.-60; r_ps[1]=15.;
	//y_ps[2]=28.-60; r_ps[2]=10.;
	//y_ps[3]=38.-60; r_ps[3]=1.;
	y_ps[0]=8.-20.; r_ps[0]=0.;
	y_ps[1]=18.-20.; r_ps[1]=15.;
	y_ps[2]=23.-20.; r_ps[2]=5.;
	y_ps[3]=28.-20.; r_ps[3]=10.;
	y_ps[4]=38.-20.; r_ps[4]=0.;
	Rot_B_spline rbs(vec3(35,-5,-160), y_ps, r_ps, Color(204, 255, 102), 0.5 , 0);
//	Rot_B_spline rbs(vec3(-25,18,-145), y_ps,r_ps, Color(151,48,48), 0.5 , 0);
	Box bx2( vec3(-18,-4,-130), vec3(-22,-12,-140),  Color(0,204,153), 0.1, 0 ) ;
	Sphere sp4 = Sphere ( vec3(-22, -10, -145),  15,  Color(0, 204, 153), 0.5, 0 ) ;
	Sphere sp5 = Sphere ( vec3(-22, 0, -145),  15,  Color(0, 204, 153), 0.5, 0 ) ;
	//Box bx2( vec3(-18,-34,-135), vec3(22,-42,-145),  Color(0,255,0), 0.1, 0 ) ;
	//Sphere sp4 = Sphere ( vec3(-22, -40, -150),  15,  Color(0,255,0), 0.5, 0 ) ;
	//Sphere sp5 = Sphere ( vec3(-22, -30, -150),  15,  Color(0,255,0), 0.5, 0 ) ;
	CSG3 csgg(sp4, bx2, sp5, Color(0, 204, 153), 0.5, 0);

    Bump_Sphere sp_bump = Bump_Sphere ( vec3(15, -50, -200),  30, 12,8,0.1, Color(250,250,250), 0.3, 0 ) ;

   figures[0] = &rbs;
//	figures[0]=&csgg;
//	figures[0]=&sp3;
 //	figures[0]=&sp_bump;
    figures[1] = &pl1;
    figures[2] = &pl2;
    figures[3] = &csgg;
    figures[4] = &sp_bump;
    figures[5] = &sp2;
    figures[6] = &cl1;
    figures[7] = &bx1;
	
	//creating image;
    ofstream image("picture.ppm");
	image << "P3" << endl;
	image << width << " " << height << endl;
	image << "255" << endl;

	cout<<"y_height="<<height<<"\nDone: ";
	for(int y = 0; y < height; y++){
		for(int x = 0 ; x < width; x++)
        {            
			//antialiasing: locally in central part of scene: nsub*nsub subpixels 
			if(fabs(x-width/2.) < width/4. && fabs(y-height/2.) < height/4.){
				Color aapixcolor(0,0,0);
				int n_sub=2;
				for(int iax=0; iax<n_sub; ++iax)
					for(int iay=0; iay<n_sub; ++iay){
						Ray aaRay = build_ray(n_sub*x+iax,n_sub*y+iay,n_sub*width,n_sub*height);
						Color subpixcolor = ray_trace (aaRay, figures,lights, 0);
						aapixcolor = aapixcolor+subpixcolor;
					}
				pixcolor=aapixcolor*(1./(n_sub*n_sub));
			}
			else{
				Ray pRay = build_ray(x,y,width,height);
				pixcolor = ray_trace (pRay, figures,lights, 0);
			}
            image << (int) min(255, pixcolor.r) << " " << 
                     (int) min(255, pixcolor.g) << " " << 
                     (int) min(255, pixcolor.b) << endl;
        }
		if ((y+1)%100 == 0) cout<<y+1<<flush;
		else if ((y+1)%10 == 0)cout<<'*'<<flush;
	}
	cout<<'\n';
        
    image.close();

	return 0;
}
