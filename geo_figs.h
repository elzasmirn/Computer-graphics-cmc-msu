#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <vector>

using namespace std;

double max(const double, const double);
double min(const double, const double);
double mix(const double&, const double&, const double&);

class vec3
{
public:
	double x_, y_, z_;
	vec3() : x_(0), y_(0), z_(0) {}
	vec3(double x) : x_(x), y_(x), z_(x) {}
	vec3(double x, double y, double z) : x_(x), y_(y), z_(z) {}

	vec3& normalize()
	{
		double nor2 = length2();
		if (nor2 > 0)
		{
			double invNor = 1 / sqrt(nor2);
			x_ *= invNor, y_ *= invNor, z_ *= invNor;
		}
		return *this;
	}
	vec3 operator * (const double &a) const { return vec3(x_ * a, y_ * a, z_ * a); }
	double operator * (const vec3 &v) const { return x_ * v.x_ + y_ * v.y_ + z_ * v.z_; }
	vec3 operator & (const vec3 &v) const { return vec3(y_*v.z_ - z_ * v.y_, -x_ * v.z_ + z_ * v.x_, x_*v.y_ - y_ * v.x_); }
	vec3 operator + (const vec3 &v) const { return vec3(x_ + v.x_, y_ + v.y_, z_ + v.z_); }
	vec3 operator - (const vec3 &v) const { return vec3(x_ - v.x_, y_ - v.y_, z_ - v.z_); }
	vec3 operator - () const { return vec3(-x_, -y_, -z_); }
	double length() const { return  sqrt(length2()); }
	double length2() const { return  x_ * x_ + y_ * y_ + z_ * z_; }
	friend vec3 operator / (const double& r, const vec3& v) { return vec3(r / v.x_, r / v.y_, r / v.z_); }
	friend ostream& operator << (ostream &s, const vec3& v)
	{
		s << "(" << v.x_ << " " << v.y_ << " " << v.z_ << ");";
		return s;
	}
};

struct Color
{
	double r, g, b;
	Color(double x = 0, double y = 0, double z = 0) : r(x), g(y), b(z) {}

	Color operator * (const double &a) const { return Color(r * a, g * a, b * a); }
	Color operator * (const Color &v) const { return Color(r * v.r, g * v.g, b * v.b); }
	Color operator + (const Color &v) const { return Color(r + v.r, g + v.g, b + v.b); }
};

struct Ray
{
	vec3 orig;
	vec3 dir;

	Ray(const vec3& origin, const vec3& direction) : orig(origin), dir(direction)
	{}
	Ray()
	{
		dir = orig = vec3();
	}
};

struct Figure //abstract 
{
	Color color;
	double reflect, transparent;
	virtual bool intersect(const Ray&, double&, double&) const = 0;
	virtual void getNorm(const vec3&, vec3&) const = 0;
};

struct Box : Figure
{
	vec3 boxMin, boxMax;
	Box()
	{
		boxMin = vec3(0, 0, -150);
		boxMax = vec3(50, -50, -200);
	}
	Box(const vec3& vmin, const vec3& vmax,
		const Color col, const double ref, const double trans)
	{
		boxMin = vmin;
		boxMax = vmax;
		color = col;
		reflect = ref;
		transparent = trans;
	}

	bool intersect(const Ray& r, double& t0, double& t1) const;

	void getNorm(const vec3& pHit, vec3& nHit) const;
};

struct Plane : Figure
{
	vec3 N;
	double D;
	Plane()
	{
		N = vec3(0, 1, 0);
		D = 20;
	}
	Plane(vec3 normal, double d,
		const Color col, const double ref, const double trans)
	{
		N = normal;
		D = d;
		color = col;
		reflect = ref;
		transparent = trans;
	}

	bool intersect(const Ray& r, double& t0, double& t1) const;
	void getNorm(const vec3& pHit, vec3& nHit) const;
};

struct Cylindr : Figure
{
	vec3 center;
	double radius, height;
	Plane top, bot;
	Cylindr(vec3 c, const double h, const double r,
		const Color col, const double ref, const double trans)
	{
		center = c;
		height = h;
		radius = r;
		color = col;
		reflect = ref;
		transparent = trans;
		double eps = 1e-4;
		top = Plane(vec3(0, -1, 0), c.y_ - h / 2 - eps, col, ref, trans);
		bot = Plane(vec3(0, -1, 0), c.y_ + h / 2 + eps, col, ref, trans);
	}

	bool intersect(const Ray& r, double& t0, double& t1) const;
	void getNorm(const vec3& pHit, vec3& nHit) const;
};

struct Cone : Figure
{
	vec3 center;
	double radius, height;
	Plane bot;
	Cone(vec3 c, double h, double r,
		const Color col, const double ref, const double trans)
	{
		center = c;
		radius = r;
		height = h;
		color = col;
		reflect = ref;
		transparent = trans;
		bot = Plane(vec3(0, 1, 0), c.y_ - h, col, ref, trans);
		//?        bot = Plane(vec3(0,1,0), c.y_+h, col, ref, trans);

	}
	bool intersect(const Ray& r, double& t0, double& t1) const;
	void getNorm(const vec3& pHit, vec3& nHit) const;
};

struct Sphere : Figure
{
	vec3 center;
	double radius, radius2;

	Sphere()
	{
		transparent = reflect = 0;
		center = vec3();
		radius = 100;
		radius2 = 100 * 100;
	}
	Sphere(const vec3 cen, const double r,
		const Color col, const double ref, const double trans)
	{
		center = cen;
		radius = r;
		radius2 = r * r;
		color = col,
			reflect = ref;
		transparent = trans;
	}

	bool intersect(const Ray& r, double& t0, double& t1) const;
	void getNorm(const vec3& pHit, vec3& nHit) const;
};

struct Bump_Sphere : Figure
{
	vec3 center;
	double radius, radius2;
	int nphi, ntheta;
	double dr_bump;

	Bump_Sphere()
	{
		transparent = reflect = 0;
		center = vec3();
		radius = 100;
		radius2 = 100 * 100;
		ntheta = 8;
		nphi = 4 * ntheta;
		dr_bump = 0.05;
	}
	Bump_Sphere(const vec3 cen, const double r,
		const int np, const int nt, const double dr,
		const Color col, const double ref, const double trans)
	{
		center = cen;
		radius = r;
		radius2 = r * r;
		nphi = np; ntheta = nt; dr_bump = dr;
		color = col,
			reflect = ref;
		transparent = trans;
	}

	bool intersect(const Ray& r, double& t0, double& t1) const;
	void getNorm(const vec3& pHit, vec3& nHit) const;
};