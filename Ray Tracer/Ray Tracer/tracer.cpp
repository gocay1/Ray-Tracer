#include <Windows.h>
#include <GL/glut.h>
#include <math.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>

using namespace std;

#define ImageW 600
#define ImageH 600

float framebuffer[ImageH][ImageW][3];

//RANDOM STRUCTS
class Pt {
public:
	float x, y, z;
	Pt()
	{
		x = 0;
		y = 0;
		z = 0;
	}
	Pt (float x1, float y1, float z1)
	{
		x = x1;
		y = y1;
		z = z1;
	}
};
class Vector {
public:	
	float x, y, z;
	Vector() {
		x = 0;
		y = 0;
		z = 0;
	}
	Vector(float x1, float y1, float z1) {
		x = x1;
		y = y1;
		z = z1;
	}
};

//VECTOR MATH

Vector normalize(Vector n) {
	Vector ret;

	float div = sqrt(n.x*n.x + n.y*n.y + n.z*n.z);
	ret.x = n.x / div;
	ret.y = n.y / div;
	ret.z = n.z / div;

	return ret;
}
Vector cross_product(Vector vect1, Vector vect2) {
	Vector ret;

	ret.x = vect1.y*vect2.z - vect1.z*vect2.y;
	ret.y = vect1.z*vect2.x - vect1.x*vect2.z;
	ret.z = vect1.x*vect2.y - vect1.y*vect2.x;

	return ret;
}
float dot_product(Vector vect1, Vector vect2) {
	float ret;

	ret = vect1.x*vect2.x + vect1.y*vect2.y + vect1.z*vect2.z;
	return ret;
}

class sphere {
public:
	Pt origin;
	float radius;
	Vector ambient, diffuse, specular;
	float reflectionC;

	sphere() {
		origin = Pt(0,0,0);
		radius = 1;
		ambient = Vector(0, 0, 0.1);
		diffuse = Vector(0, 0, 0.7);
		specular = Vector(0.5, 0.5, 0.5);
		reflectionC = 0;
	}
	sphere(Pt or, float r) {
		origin = or;
		radius = r;
		ambient = Vector(0, 0, 0.1);
		diffuse = Vector(0, 0, 0.7);
		specular = Vector(0.5, 0.5, 0.5);
		reflectionC = 0;
	}
};
class plane {
public:
	Pt planePoint;
	Vector normal;
	Vector ambient, diffuse, specular;
	float reflectionC;

	plane() {
		planePoint = Pt(0,0,10);
		normal = Vector(0,0,-1);
		ambient = Vector(0, 0.1, 0.1);
		diffuse = Vector(0.7, 0.9, 0.9);
		ambient = Vector(0.6, 0.6, 0.6);
		reflectionC = 0;
	}
	plane(Pt origin, Vector norm) {
		planePoint = origin;
		normal = norm;
		ambient = Vector(0, 0.1, 0.1);
		diffuse = Vector(0.7, 0.9, 0.9);
		ambient = Vector(0.6, 0.6, 0.6);
		reflectionC = 0;
	}
};
class cylinder{
public:
	Vector a;
	Pt c;
	float r;
	Vector ambient, diffuse, specular;
	float reflectionC;

	cylinder() {
		a = Vector(0,1,0);
		c = Pt(0,0,0);
		r = 1;
		ambient = Vector(0.1, 0.1, 0);
		diffuse = Vector(0.7, 0.7, 0);
		specular = Vector(0.5, 0.5, 0.5);
		reflectionC = 0;
	}
	cylinder(Vector dir, Pt center, float rad) {
		a = dir;
		c = center;
		r = rad;
		ambient = Vector(0.1, 0.1, 0);
		diffuse = Vector(0.7, 0.7, 0);
		specular = Vector(0.5, 0.5, 0.5);
		reflectionC = 0;
	}
};
class matrix {
public:
	float mat[4][4];
	matrix() {
		mat[0][0] = 0;
		mat[0][1] = 0;
		mat[0][2] = 0;
		mat[0][3] = 0;
		mat[1][0] = 0;
		mat[1][1] = 0;
		mat[1][2] = 0;
		mat[1][3] = 0;
		mat[2][0] = 0;
		mat[2][1] = 0;
		mat[2][2] = 0;
		mat[2][3] = 0;
		mat[3][0] = 0;
		mat[3][1] = 0;
		mat[3][2] = 0;
		mat[3][3] = 0;
	}
};
class elipse {
public:
	Pt origin;
	float radius;
	Vector deformed;
	Vector ambient, diffuse, specular;
	matrix def;
	float reflectionC;

	elipse() {
		origin = Pt(0, 0, 0);
		deformed = Vector(2,1,1);
		radius = 1;
		ambient = Vector(0.6, 0, 0.6);
		diffuse = Vector(0.7, 0, 0.7);
		specular = Vector(0.5, 0.5, 0.5);
		calcMat();
		reflectionC = 0;
	}
	elipse(Pt or , float r, Vector deform) {
		origin = or ;
		radius = r;
		deformed = deform;
		ambient = Vector(0, 0, 0.1);
		diffuse = Vector(0, 0, 0.7);
		specular = Vector(0.5, 0.5, 0.5);
		calcMat();
		reflectionC = 0;
	}
	void calcMat() {

		def.mat[0][0] = 1.0 / deformed.x;
		def.mat[1][1] = 1.0 / deformed.y;
		def.mat[2][2] = 1.0 / deformed.z;

		float m;
		if (deformed.x != 1) {
			m = deformed.x;
		}
		if (deformed.y != 1) {
			m = deformed.y;
		}
		if (deformed.z != 1) {
			m = deformed.z;
		}

		def.mat[0][3] = (1.0 - (1 / m)) - origin.x;
		def.mat[1][3] = (1.0 - (1 / m)) - origin.y;
		def.mat[2][3] = (1.0 - (1 / m)) - origin.z;

	}
};
class screen {
public:
	Pt tl, br;

	screen() {
		tl = Pt(-1,1,-2);
		br = Pt(1,-1,-2);
	}
};
class ray {
public:
	Pt origin, p2;
	Vector direction;

	ray(){}
	ray(Pt pp1, Pt pp2) {
		origin = pp1;
		p2 = pp2;
		direction.x = p2.x - origin.x;
		direction.y = p2.y - origin.y;
		direction.z = p2.z - origin.z;
		direction = normalize(direction);

	}
};
class rayTraceOut {
public:
	float t;
	Vector ambient;
	Vector diffuse;
	Vector specular;
	Vector normal;
	Pt calc;
	float reflectionC;

	bool operator < (const rayTraceOut &e)
	{
			return t < e.t;
	}
};
class light {
public:
	Pt source;
	Vector intensity;
	light() {
		source = Pt(10,10,-10);
		intensity = Vector(1, 1, 1);
	}
	light(Pt sour, Vector inte) {
		source = sour;
		intensity = inte;
	}

};

//Global Definitions
screen scr;
Pt viewer;
vector<sphere> spheres;
vector<plane> planes;
vector<cylinder> cylinders;
vector<elipse> elipses;
Vector ambient;
vector<light> lightSources;

//BABY DRAW ME SOME BUFFER LOVE!!

void setFramebuffer(int x, int y, float z, float R, float G, float B)
{
	//if (zBuffer[y][x][0] > z) {
	//	zBuffer[y][x][0] = z;
		// changes the origin from the lower-left corner to the upper-left corner
		y = ImageH - 1 - y;
		if (R <= 1.0)
			if (R >= 0.0)
				framebuffer[y][x][0] = R;
			else
				framebuffer[y][x][0] = 0.0;
		else
			framebuffer[y][x][0] = 1.0;
		if (G <= 1.0)
			if (G >= 0.0)
				framebuffer[y][x][1] = G;
			else
				framebuffer[y][x][1] = 0.0;
		else
			framebuffer[y][x][1] = 1.0;
		if (B <= 1.0)
			if (B >= 0.0)
				framebuffer[y][x][2] = B;
			else
				framebuffer[y][x][2] = 0.0;
		else
			framebuffer[y][x][2] = 1.0;
	//}
}
// Draws the scene
void drawit(void)
{
	glDrawPixels(ImageW, ImageH, GL_RGB, GL_FLOAT, framebuffer);
	glFlush();
}
// Clears framebuffer to black
void clearFramebuffer()
{
	int i, j;

	for (i = 0; i<ImageH; i++) {
		for (j = 0; j<ImageW; j++) {
			framebuffer[i][j][0] = 0.0;
			framebuffer[i][j][1] = 0.0;
			framebuffer[i][j][2] = 0.0;
		}
	}
}

//LIGHTING NONSENSE

float angle_vect(Vector vect1, Vector vect2) {
	float dot = dot_product(vect1, vect2);
	float nrm = sqrt(vect1.x*vect1.x + vect1.y*vect1.y + vect1.z*vect1.z)*
		sqrt(vect2.x*vect2.x + vect2.y*vect2.y + vect2.z*vect2.z);
	float angle = acos(dot / nrm);

	return angle;
}

bool lightRayScan(Pt light, Pt space) {
	ray fire = ray(space, light);
	vector<rayTraceOut> traceResults;
	for (sphere s : spheres) {

		float a, b, c;

		a = dot_product(fire.direction, fire.direction);

		Vector v = fire.direction;
		v.x = v.x * 2.0;
		v.y = v.y * 2.0;
		v.z = v.z * 2.0;

		Vector l;
		l.x = fire.origin.x - s.origin.x;
		l.y = fire.origin.y - s.origin.y;
		l.z = fire.origin.z - s.origin.z;

		b = dot_product(v, l);
		c = dot_product(l, l) - (s.radius*s.radius);

		//cout << a << "     "  << b << "     " << c << endl;
		if (((b*b) - (4.0*a*c)) > 0) {
			rayTraceOut r1;
			r1.t = ((-b) + sqrt((b*b) - (4.0 * a*c))) / (2.0*a);

			Pt calc;
			calc.x = (fire.origin.x + r1.t*fire.direction.x);
			calc.y = (fire.origin.y + r1.t*fire.direction.y);
			calc.z = (fire.origin.z + r1.t*fire.direction.z);

			if (r1.t > 0.01) {
				traceResults.push_back(r1);
			}
			r1.t = ((-b) - sqrt((b*b) - (4.0 * a*c))) / (2.0*a);

			calc.x = (fire.origin.x + r1.t*fire.direction.x);
			calc.y = (fire.origin.y + r1.t*fire.direction.y);
			calc.z = (fire.origin.z + r1.t*fire.direction.z);
			
			if (r1.t > 0.01) {
				traceResults.push_back(r1);
			}
		}
		else if (((b*b) - (4.0 * a*c)) == 0) {
			rayTraceOut r1;
			r1.t = (-b) / (2.0*a);
			r1.ambient = s.ambient;
			r1.diffuse = s.diffuse;
			r1.diffuse = s.specular;

			Pt calc;
			calc.x = (fire.origin.x + r1.t*fire.direction.x);
			calc.y = (fire.origin.y + r1.t*fire.direction.y);
			calc.z = (fire.origin.z + r1.t*fire.direction.z);
			
			if (r1.t > 0.01) {
				traceResults.push_back(r1);
			}
		}
	}
	for (plane p : planes) {
		Pt top;
		top.x = p.planePoint.x - fire.origin.x;
		top.y = p.planePoint.y - fire.origin.y;
		top.z = p.planePoint.z - fire.origin.z;

		float t = p.normal.x * top.x + p.normal.y * top.y + p.normal.z * top.z;
		t = t / dot_product(p.normal, fire.direction);
		if (t > 0.01) {
			rayTraceOut out;
			out.t = t;
			out.normal = p.normal;
			traceResults.push_back(out);
		}
	}
	for (cylinder cy : cylinders) {
		float a, b, c;

		cy.a = normalize(cy.a);

		//THE A

		a = dot_product(fire.direction, cy.a);

		Vector a1 = cy.a;

		a1.x = a1.x * a;
		a1.y = a1.y * a;
		a1.z = a1.z * a;

		Vector a2 = a1;

		a2.x = fire.direction.x - a2.x;
		a2.y = fire.direction.y - a2.y;
		a2.z = fire.direction.z - a2.z;

		a = dot_product(a2, a2);

		//THE B

		Vector b1;
		Vector deltaP;
		deltaP.x = fire.origin.x - cy.c.x;
		deltaP.y = fire.origin.y - cy.c.y;
		deltaP.z = fire.origin.z - cy.c.z;
		//deltaP = normalize(deltaP);

		b = dot_product(deltaP, cy.a);
		b1.x = deltaP.x - (b*cy.a.x);
		b1.y = deltaP.y - (b*cy.a.y);
		b1.z = deltaP.z - (b*cy.a.z);

		b = 2 * dot_product(a2, b1);

		//THE C

		c = dot_product(b1, b1) - (cy.r*cy.r);

		if (((b*b) - (4.0*a*c)) > 0) {
			rayTraceOut r1;
			r1.t = ((-b) + sqrt((b*b) - (4.0 * a*c))) / (2.0*a);
			r1.ambient = cy.ambient;
			r1.diffuse = cy.diffuse;
			r1.specular = cy.specular;

			Pt calc;
			calc.x = (fire.origin.x + r1.t*fire.direction.x);
			calc.y = (fire.origin.y + r1.t*fire.direction.y);
			calc.z = (fire.origin.z + r1.t*fire.direction.z);


			if (r1.t > 0.01) {
				traceResults.push_back(r1);
			}
			r1.t = ((-b) - sqrt((b*b) - (4.0 * a*c))) / (2.0*a);

			calc.x = (fire.origin.x + r1.t*fire.direction.x);
			calc.y = (fire.origin.y + r1.t*fire.direction.y);
			calc.z = (fire.origin.z + r1.t*fire.direction.z);

			if (r1.t > 0.01) {
				traceResults.push_back(r1);
			}
		}
		else if (((b*b) - (4.0 * a*c)) == 0) {
			rayTraceOut r1;
			r1.t = (-b) / (2.0*a);
			r1.ambient = cy.ambient;
			r1.diffuse = cy.diffuse;
			r1.diffuse = cy.specular;

			Pt calc;
			calc.x = (fire.origin.x + r1.t*fire.direction.x);
			calc.y = (fire.origin.y + r1.t*fire.direction.y);
			calc.z = (fire.origin.z + r1.t*fire.direction.z);

			if (r1.t > 0.01) {
				traceResults.push_back(r1);
			}
		}

	}
	for (elipse s : elipses) {

		ray alt_fire = fire;
		alt_fire.origin.x = alt_fire.origin.x*s.def.mat[0][0] + s.def.mat[0][3];
		alt_fire.origin.y = alt_fire.origin.y*s.def.mat[1][1] + s.def.mat[1][3];
		alt_fire.origin.z = alt_fire.origin.z*s.def.mat[2][2] + s.def.mat[2][3];

		alt_fire.direction.x = alt_fire.direction.x*s.def.mat[0][0];
		alt_fire.direction.y = alt_fire.direction.y*s.def.mat[1][1];
		alt_fire.direction.z = alt_fire.direction.z*s.def.mat[2][2];
		alt_fire.direction = normalize(alt_fire.direction);

		float a, b, c;

		a = dot_product(alt_fire.direction, alt_fire.direction);

		Vector v = alt_fire.direction;
		v.x = v.x * 2.0;
		v.y = v.y * 2.0;
		v.z = v.z * 2.0;

		Vector l;
		l.x = alt_fire.origin.x - s.origin.x;
		l.y = alt_fire.origin.y - s.origin.y;
		l.z = alt_fire.origin.z - s.origin.z;

		b = dot_product(v, l);
		c = dot_product(l, l) - (s.radius*s.radius);

		//cout << a << "     "  << b << "     " << c << endl;
		if (((b*b) - (4.0*a*c)) > 0) {
			rayTraceOut r1;
			r1.t = ((-b) + sqrt((b*b) - (4.0 * a*c))) / (2.0*a);
			r1.ambient = s.ambient;
			r1.diffuse = s.diffuse;
			r1.specular = s.specular;

			Pt calc;
			calc.x = (alt_fire.origin.x + r1.t*alt_fire.direction.x);
			calc.y = (alt_fire.origin.y + r1.t*alt_fire.direction.y);
			calc.z = (alt_fire.origin.z + r1.t*alt_fire.direction.z);

			if (r1.t > 0.01) {
				traceResults.push_back(r1);
			}
			r1.t = ((-b) - sqrt((b*b) - (4.0 * a*c))) / (2.0*a);

			calc.x = (alt_fire.origin.x + r1.t*alt_fire.direction.x);
			calc.y = (alt_fire.origin.y + r1.t*alt_fire.direction.y);
			calc.z = (alt_fire.origin.z + r1.t*alt_fire.direction.z);

			if (r1.t > 0.01) {
				traceResults.push_back(r1);
			}
		}
		else if (((b*b) - (4.0 * a*c)) == 0) {
			rayTraceOut r1;
			r1.t = (-b) / (2.0*a);
			r1.ambient = s.ambient;
			r1.diffuse = s.diffuse;
			r1.diffuse = s.specular;

			Pt calc;
			calc.x = (alt_fire.origin.x + r1.t*alt_fire.direction.x);
			calc.y = (alt_fire.origin.y + r1.t*alt_fire.direction.y);
			calc.z = (alt_fire.origin.z + r1.t*alt_fire.direction.z);

			if (r1.t > 0.01) {
				traceResults.push_back(r1);
			}
		}
	}
	sort(traceResults.begin(), traceResults.end());
	if (traceResults.size() > 1) {
		
		return true;
	}
	else {
		return false;
	}
}

Vector find_color(Vector normal, Vector eye, Vector ambientc, Vector diffusionc, Vector specularc, Pt space) {

	Vector ret;

	ret.x = ret.x + ambientc.x*ambient.x;
	ret.y = ret.y + ambientc.y*ambient.y;
	ret.z = ret.z + ambientc.z*ambient.z;

	for (light lightS : lightSources) {
		
		Vector light;
		light.x = lightS.source.x - space.x;
		light.y = lightS.source.y - space.y;
		light.z = lightS.source.z - space.z;

		//light = Vector(1, 1, -1);
		light = normalize(light);

		bool inShadow = lightRayScan(lightS.source, space);

			float multiplier = 2.0*dot_product(light, normal);
			Vector r = normal;
			r.x = r.x*multiplier - light.x;
			r.y = r.y*multiplier - light.y;
			r.z = r.z*multiplier - light.z;
			r = normalize(r);

			//ANGULAR WAY

			float alpha = angle_vect(r, eye);
			float angle = angle_vect(normal, light);


			if (inShadow && (diffusionc.x*cos(angle) + specularc.x*pow(cos(alpha), 5)) < 0) {
				ret.x += lightS.intensity.x*(diffusionc.x*cos(angle) + specularc.x*pow(cos(alpha), 5));
			}
			if (inShadow && (diffusionc.y*cos(angle) + specularc.y*pow(cos(alpha), 5)) < 0) {
				ret.y += lightS.intensity.y*(diffusionc.y*cos(angle) + specularc.y*pow(cos(alpha), 5));
			}
			if (inShadow && (diffusionc.z*cos(angle) + specularc.z*pow(cos(alpha), 5)) < 0) {
				ret.z += lightS.intensity.z*(diffusionc.z*cos(angle) + specularc.z*pow(cos(alpha), 5));
			}
			if (!inShadow) {
				ret.x += lightS.intensity.x*(diffusionc.x*cos(angle) + specularc.x*pow(cos(alpha), 5));
				ret.y += lightS.intensity.y*(diffusionc.y*cos(angle) + specularc.y*pow(cos(alpha), 5));
				ret.z += lightS.intensity.z*(diffusionc.z*cos(angle) + specularc.z*pow(cos(alpha), 5));
			}

			//if (ret.x > ret.z && ret.y > ret.z) {
			//	cout << alpha << endl;
			//	cout << ret.x << endl;
			//	cout << ret.y << endl;
			//	cout << ret.z << endl;
			//	cout << "---------" << endl;
			//}
		
	}
	return ret;
}

//DRAWING ALGORITHMS

Pt window_coor(float x, float y) {
	Pt ret;
	ret.x = (((float)ImageW - 1.0) / 2.0)*(1.0 + x);
	ret.y = (((float)ImageH - 1.0) / 2.0)*(1.0 - y);
	return ret;
}

Vector subtractVectors(Pt a, Pt b) {
	Vector ret;

	ret.x = a.x - b.x;
	ret.y = a.y - b.y;
	ret.z = a.z - b.z;

	return ret;
}
Vector subtractVectors(Vector a, Pt b) {
	Vector ret;

	ret.x = a.x - b.x;
	ret.y = a.y - b.y;
	ret.z = a.z - b.z;

	return ret;
}
Vector subtractVectors(Vector a, Vector b) {
	Vector ret;

	ret.x = a.x - b.x;
	ret.y = a.y - b.y;
	ret.z = a.z - b.z;

	return ret;
}

Vector scalarMult(Vector v, float x) {
	Vector ret;
	ret.x = v.x*x;
	ret.y = v.y*x;
	ret.z = v.z*x;
	return ret;
}

Vector cylinderNormal(Pt point, Pt origin, Vector a, float r) {
	Vector ret;
	Vector sub;

	sub = subtractVectors(point, origin);
	ret = subtractVectors(sub, scalarMult(a,dot_product(sub,a)));
	ret = scalarMult(ret, 1/r);
	
	return ret;
}

Vector sphereNormal(Pt origin, Pt other) {
	Vector normal;
	normal.x = other.x - origin.x;
	normal.y = other.y - origin.y;
	normal.z = other.z - origin.z;
	normal = normalize(normal);
	return normal;
};


Vector reflectionColor(Pt picked, Vector eye, Vector normal, int depth) {
	Vector reflected = scalarMult(normal, (2 * dot_product(eye, normal)));
	reflected = subtractVectors(reflected, eye);
	normalize(reflected);

	ray fire;
	fire.origin = picked;
	fire.direction = reflected;
	vector<rayTraceOut> traceResults;
	for (sphere s : spheres) {

		float a, b, c;

		a = dot_product(fire.direction, fire.direction);

		Vector v = fire.direction;
		v.x = v.x * 2.0;
		v.y = v.y * 2.0;
		v.z = v.z * 2.0;

		Vector l;
		l.x = fire.origin.x - s.origin.x;
		l.y = fire.origin.y - s.origin.y;
		l.z = fire.origin.z - s.origin.z;

		b = dot_product(v, l);
		c = dot_product(l, l) - (s.radius*s.radius);

		//cout << a << "     "  << b << "     " << c << endl;
		if (((b*b) - (4.0*a*c)) > 0) {
			rayTraceOut r1;
			r1.t = ((-b) + sqrt((b*b) - (4.0 * a*c))) / (2.0*a);
			r1.ambient = s.ambient;
			r1.diffuse = s.diffuse;
			r1.specular = s.specular;

			Pt calc;
			calc.x = (fire.origin.x + r1.t*fire.direction.x);
			calc.y = (fire.origin.y + r1.t*fire.direction.y);
			calc.z = (fire.origin.z + r1.t*fire.direction.z);


			r1.normal = sphereNormal(s.origin, calc);
			r1.calc = calc;
			r1.reflectionC = s.reflectionC;

			if (r1.t > 0.01) {
				traceResults.push_back(r1);
			}
			r1.t = ((-b) - sqrt((b*b) - (4.0 * a*c))) / (2.0*a);

			calc.x = (fire.origin.x + r1.t*fire.direction.x);
			calc.y = (fire.origin.y + r1.t*fire.direction.y);
			calc.z = (fire.origin.z + r1.t*fire.direction.z);

			r1.normal = sphereNormal(s.origin, calc);
			r1.calc = calc;
			r1.reflectionC = s.reflectionC;

			if (r1.t > 0.01) {
				traceResults.push_back(r1);
			}
		}
		else if (((b*b) - (4.0 * a*c)) == 0) {
			rayTraceOut r1;
			r1.t = (-b) / (2.0*a);
			r1.ambient = s.ambient;
			r1.diffuse = s.diffuse;
			r1.diffuse = s.specular;

			Pt calc;
			calc.x = (fire.origin.x + r1.t*fire.direction.x);
			calc.y = (fire.origin.y + r1.t*fire.direction.y);
			calc.z = (fire.origin.z + r1.t*fire.direction.z);

			r1.normal = sphereNormal(s.origin, calc);
			r1.calc = calc;
			r1.reflectionC = s.reflectionC;
			if (r1.t > 0.01) {
				traceResults.push_back(r1);
			}
		}
	}
	for (plane p : planes) {
		Pt top;
		top.x = p.planePoint.x - fire.origin.x;
		top.y = p.planePoint.y - fire.origin.y;
		top.z = p.planePoint.z - fire.origin.z;

		float t = p.normal.x * top.x + p.normal.y * top.y + p.normal.z * top.z;
		if (fabs(dot_product(p.normal, fire.direction)) > 0.00001) {
			t = t / dot_product(p.normal, fire.direction);
		}
		else {
			t = -1;
		}
		if (t > 0.01) {
			rayTraceOut out;
			out.t = t;
			out.ambient = p.ambient;
			out.diffuse = p.diffuse;
			out.specular = p.specular;
			out.normal = p.normal;

			Pt calc;
			calc.x = (fire.origin.x + out.t*fire.direction.x);
			calc.y = (fire.origin.y + out.t*fire.direction.y);
			calc.z = (fire.origin.z + out.t*fire.direction.z);

			out.calc = calc;
			out.reflectionC = p.reflectionC;

			traceResults.push_back(out);
		}
	}
	for (cylinder cy : cylinders) {
		float a, b, c;

		cy.a = normalize(cy.a);

		//THE A

		a = dot_product(fire.direction, cy.a);

		Vector a1 = cy.a;

		a1.x = a1.x * a;
		a1.y = a1.y * a;
		a1.z = a1.z * a;

		Vector a2 = a1;

		a2.x = fire.direction.x - a2.x;
		a2.y = fire.direction.y - a2.y;
		a2.z = fire.direction.z - a2.z;

		a = dot_product(a2, a2);

		//THE B

		Vector b1;
		Vector deltaP;
		deltaP.x = fire.origin.x - cy.c.x;
		deltaP.y = fire.origin.y - cy.c.y;
		deltaP.z = fire.origin.z - cy.c.z;
		//deltaP = normalize(deltaP);

		b = dot_product(deltaP, cy.a);
		b1.x = deltaP.x - (b*cy.a.x);
		b1.y = deltaP.y - (b*cy.a.y);
		b1.z = deltaP.z - (b*cy.a.z);

		b = 2 * dot_product(a2, b1);

		//THE C

		c = dot_product(b1, b1) - (cy.r*cy.r);

		if (((b*b) - (4.0*a*c)) > 0) {
			rayTraceOut r1;
			r1.t = ((-b) + sqrt((b*b) - (4.0 * a*c))) / (2.0*a);
			r1.ambient = cy.ambient;
			r1.diffuse = cy.diffuse;
			r1.specular = cy.specular;

			Pt calc;
			calc.x = (fire.origin.x + r1.t*fire.direction.x);
			calc.y = (fire.origin.y + r1.t*fire.direction.y);
			calc.z = (fire.origin.z + r1.t*fire.direction.z);


			r1.normal = cylinderNormal(calc, cy.c, cy.a, cy.r);
			r1.calc = calc;
			r1.reflectionC = cy.reflectionC;

			if (r1.t > 0.01) {
				traceResults.push_back(r1);
			}
			r1.t = ((-b) - sqrt((b*b) - (4.0 * a*c))) / (2.0*a);

			calc.x = (fire.origin.x + r1.t*fire.direction.x);
			calc.y = (fire.origin.y + r1.t*fire.direction.y);
			calc.z = (fire.origin.z + r1.t*fire.direction.z);

			r1.normal = cylinderNormal(calc, cy.c, cy.a, cy.r);
			r1.calc = calc;
			r1.reflectionC = cy.reflectionC;

			if (r1.t > 0.01) {
				traceResults.push_back(r1);
			}
		}
		else if (((b*b) - (4.0 * a*c)) == 0) {
			rayTraceOut r1;
			r1.t = (-b) / (2.0*a);
			r1.ambient = cy.ambient;
			r1.diffuse = cy.diffuse;
			r1.diffuse = cy.specular;

			Pt calc;
			calc.x = (fire.origin.x + r1.t*fire.direction.x);
			calc.y = (fire.origin.y + r1.t*fire.direction.y);
			calc.z = (fire.origin.z + r1.t*fire.direction.z);

			r1.normal = cylinderNormal(calc, cy.c, cy.a, cy.r);
			r1.calc = calc;
			r1.reflectionC = cy.reflectionC;

			if (r1.t > 0.01) {
				traceResults.push_back(r1);
			}
		}

	}
	for (elipse s : elipses) {

		ray alt_fire = fire;
		alt_fire.origin.x = alt_fire.origin.x*s.def.mat[0][0] + s.def.mat[0][3];
		alt_fire.origin.y = alt_fire.origin.y*s.def.mat[1][1] + s.def.mat[1][3];
		alt_fire.origin.z = alt_fire.origin.z*s.def.mat[2][2] + s.def.mat[2][3];

		alt_fire.direction.x = alt_fire.direction.x*s.def.mat[0][0];
		alt_fire.direction.y = alt_fire.direction.y*s.def.mat[1][1];
		alt_fire.direction.z = alt_fire.direction.z*s.def.mat[2][2];
		alt_fire.direction = normalize(alt_fire.direction);

		float a, b, c;

		a = dot_product(alt_fire.direction, alt_fire.direction);

		Vector v = alt_fire.direction;
		v.x = v.x * 2.0;
		v.y = v.y * 2.0;
		v.z = v.z * 2.0;

		Vector l;
		l.x = alt_fire.origin.x - s.origin.x;
		l.y = alt_fire.origin.y - s.origin.y;
		l.z = alt_fire.origin.z - s.origin.z;

		b = dot_product(v, l);
		c = dot_product(l, l) - (s.radius*s.radius);

		//cout << a << "     "  << b << "     " << c << endl;
		if (((b*b) - (4.0*a*c)) > 0) {
			rayTraceOut r1;
			r1.t = ((-b) + sqrt((b*b) - (4.0 * a*c))) / (2.0*a);
			r1.ambient = s.ambient;
			r1.diffuse = s.diffuse;
			r1.specular = s.specular;

			Pt calc;
			calc.x = (alt_fire.origin.x + r1.t*alt_fire.direction.x);
			calc.y = (alt_fire.origin.y + r1.t*alt_fire.direction.y);
			calc.z = (alt_fire.origin.z + r1.t*alt_fire.direction.z);


			r1.normal = sphereNormal(s.origin, calc);

			r1.normal.x = r1.normal.x*s.def.mat[0][0];
			r1.normal.y = r1.normal.y*s.def.mat[1][1];
			r1.normal.z = r1.normal.z*s.def.mat[2][2];

			r1.normal = normalize(r1.normal);

			calc.x = (fire.origin.x + r1.t*fire.direction.x);
			calc.y = (fire.origin.y + r1.t*fire.direction.y);
			calc.z = (fire.origin.z + r1.t*fire.direction.z);

			r1.calc = calc;
			r1.reflectionC = s.reflectionC;
			if (r1.t > 0.01) {
				traceResults.push_back(r1);
			}
			r1.t = ((-b) - sqrt((b*b) - (4.0 * a*c))) / (2.0*a);

			calc.x = (alt_fire.origin.x + r1.t*alt_fire.direction.x);
			calc.y = (alt_fire.origin.y + r1.t*alt_fire.direction.y);
			calc.z = (alt_fire.origin.z + r1.t*alt_fire.direction.z);

			r1.normal = sphereNormal(s.origin, calc);
			r1.normal.x = r1.normal.x*s.def.mat[0][0];
			r1.normal.y = r1.normal.y*s.def.mat[1][1];
			r1.normal.z = r1.normal.z*s.def.mat[2][2];
			r1.normal = normalize(r1.normal);

			calc.x = (fire.origin.x + r1.t*fire.direction.x);
			calc.y = (fire.origin.y + r1.t*fire.direction.y);
			calc.z = (fire.origin.z + r1.t*fire.direction.z);

			r1.calc = calc;
			r1.reflectionC = s.reflectionC;
			if (r1.t > 0.01) {
				traceResults.push_back(r1);
			}
		}
		else if (((b*b) - (4.0 * a*c)) == 0) {
			rayTraceOut r1;
			r1.t = (-b) / (2.0*a);
			r1.ambient = s.ambient;
			r1.diffuse = s.diffuse;
			r1.diffuse = s.specular;

			Pt calc;
			calc.x = (alt_fire.origin.x + r1.t*alt_fire.direction.x);
			calc.y = (alt_fire.origin.y + r1.t*alt_fire.direction.y);
			calc.z = (alt_fire.origin.z + r1.t*alt_fire.direction.z);

			r1.normal = sphereNormal(s.origin, calc);

			r1.normal.x = r1.normal.x*s.def.mat[0][0];
			r1.normal.y = r1.normal.y*s.def.mat[1][1];
			r1.normal.z = r1.normal.z*s.def.mat[2][2];

			r1.normal = normalize(r1.normal);

			calc.x = (fire.origin.x + r1.t*fire.direction.x);
			calc.y = (fire.origin.y + r1.t*fire.direction.y);
			calc.z = (fire.origin.z + r1.t*fire.direction.z);

			r1.calc = calc;
			r1.reflectionC = s.reflectionC;

			if (r1.t > 0.01) {
				traceResults.push_back(r1);
			}
		}
	}
	sort(traceResults.begin(), traceResults.end());
	if (traceResults.size() > 0) {
		float t1 = traceResults[0].t;

		Pt calc = traceResults[0].calc;

		Vector normal = traceResults[0].normal;

		Vector color = find_color(normal, eye, traceResults[0].ambient,
			traceResults[0].diffuse,
			traceResults[0].specular, calc);
		if (color.x < 0) {
			color.x = 0;
		}
		if (color.y < 0) {
			color.y = 0;
		}
		if (color.z < 0) {
			color.z = 0;
		}
		depth++;
		if (traceResults[0].reflectionC > 0 && depth < 6) {
			Vector reflectedColor = reflectionColor(traceResults[0].calc, eye, traceResults[0].normal, depth);
			color.x += reflectedColor.x*traceResults[0].reflectionC;
			color.y += reflectedColor.y*traceResults[0].reflectionC;
			color.z += reflectedColor.z*traceResults[0].reflectionC;
		}

		return color;
	}
	return Vector(0,0,0);


}


void rayTrace() {
	for (float x = -1; x < 1; x += 0.0033333333) {
		for (float y = -1; y < 1; y += 0.0033333333) {
			Pt window = window_coor(x, y);
			int windowx = window.x;
			int windowy = window.y;
			//setFramebuffer(windowx, windowy, 0, 0.5, 0);
			ray fire = ray(viewer, Pt(x, y, scr.tl.z));
			vector<rayTraceOut> traceResults;
			for (sphere s : spheres) {

				float a, b, c;

				a = dot_product(fire.direction, fire.direction);

				Vector v = fire.direction;
				v.x = v.x * 2.0;
				v.y = v.y * 2.0;
				v.z = v.z * 2.0;

				Vector l;
				l.x = fire.origin.x - s.origin.x;
				l.y = fire.origin.y - s.origin.y;
				l.z = fire.origin.z - s.origin.z;

				b = dot_product(v, l);
				c = dot_product(l, l) - (s.radius*s.radius);

				//cout << a << "     "  << b << "     " << c << endl;
				if (((b*b) - (4.0*a*c)) > 0) {
					rayTraceOut r1;
					r1.t = ((-b) + sqrt((b*b) - (4.0 * a*c))) / (2.0*a);
					r1.ambient = s.ambient;
					r1.diffuse = s.diffuse;
					r1.specular = s.specular;

					Pt calc;
					calc.x = (fire.origin.x + r1.t*fire.direction.x);
					calc.y = (fire.origin.y + r1.t*fire.direction.y);
					calc.z = (fire.origin.z + r1.t*fire.direction.z);

					
					r1.normal = sphereNormal(s.origin, calc);
					r1.calc = calc;
					r1.reflectionC = s.reflectionC;

					if (r1.t > 0) {
						traceResults.push_back(r1);
					}
					r1.t = ((-b) - sqrt((b*b) - (4.0 * a*c))) / (2.0*a);

					calc.x = (fire.origin.x + r1.t*fire.direction.x);
					calc.y = (fire.origin.y + r1.t*fire.direction.y);
					calc.z = (fire.origin.z + r1.t*fire.direction.z);

					r1.normal = sphereNormal(s.origin, calc);
					r1.calc = calc;
					r1.reflectionC = s.reflectionC;
					
					if (r1.t > 0) {
						traceResults.push_back(r1);
					}
				}
				else if (((b*b) - (4.0 * a*c)) == 0) {
					rayTraceOut r1;
					r1.t = (-b) / (2.0*a);
					r1.ambient = s.ambient;
					r1.diffuse = s.diffuse;
					r1.diffuse = s.specular;

					Pt calc;
					calc.x = (fire.origin.x + r1.t*fire.direction.x);
					calc.y = (fire.origin.y + r1.t*fire.direction.y);
					calc.z = (fire.origin.z + r1.t*fire.direction.z);

					r1.normal = sphereNormal(s.origin, calc);
					r1.calc = calc;
					r1.reflectionC = s.reflectionC;
					if (r1.t) {
						traceResults.push_back(r1);
					}
				}
			}
			for (plane p : planes) {
				Pt top;
				top.x = p.planePoint.x - fire.origin.x;
				top.y = p.planePoint.y - fire.origin.y;
				top.z = p.planePoint.z - fire.origin.z;

				float t = p.normal.x * top.x + p.normal.y * top.y + p.normal.z * top.z;
				if (fabs(dot_product(p.normal, fire.direction)) > 0.00001) {
					t = t / dot_product(p.normal, fire.direction);
				}
				else {
					t = -1;
				}
				if (t > 0) {
					rayTraceOut out;
					out.t = t;
					out.ambient = p.ambient;
					out.diffuse = p.diffuse;
					out.specular = p.specular;
					out.normal = p.normal;

					Pt calc;
					calc.x = (fire.origin.x + out.t*fire.direction.x);
					calc.y = (fire.origin.y + out.t*fire.direction.y);
					calc.z = (fire.origin.z + out.t*fire.direction.z);

					out.calc = calc;
					out.reflectionC = p.reflectionC;

					traceResults.push_back(out);
				}
			}
			for (cylinder cy : cylinders) {
				float a, b, c;

				cy.a = normalize(cy.a);

				//THE A

				a = dot_product(fire.direction, cy.a);

				Vector a1 = cy.a;

				a1.x = a1.x * a;
				a1.y = a1.y * a;
				a1.z = a1.z * a;

				Vector a2 = a1;

				a2.x = fire.direction.x - a2.x;
				a2.y = fire.direction.y - a2.y;
				a2.z = fire.direction.z - a2.z;

				a = dot_product(a2, a2);

				//THE B

				Vector b1;
				Vector deltaP;
				deltaP.x = fire.origin.x - cy.c.x;
				deltaP.y = fire.origin.y - cy.c.y;
				deltaP.z = fire.origin.z - cy.c.z;
				//deltaP = normalize(deltaP);

				b = dot_product(deltaP, cy.a);
				b1.x = deltaP.x - (b*cy.a.x);
				b1.y = deltaP.y - (b*cy.a.y);
				b1.z = deltaP.z - (b*cy.a.z);

				b = 2 * dot_product(a2, b1);

				//THE C

				c = dot_product(b1,b1) - (cy.r*cy.r);

				if (((b*b) - (4.0*a*c)) > 0) {
					rayTraceOut r1;
					r1.t = ((-b) + sqrt((b*b) - (4.0 * a*c))) / (2.0*a);
					r1.ambient = cy.ambient;
					r1.diffuse = cy.diffuse;
					r1.specular = cy.specular;

					Pt calc;
					calc.x = (fire.origin.x + r1.t*fire.direction.x);
					calc.y = (fire.origin.y + r1.t*fire.direction.y);
					calc.z = (fire.origin.z + r1.t*fire.direction.z);


					r1.normal = cylinderNormal(calc,cy.c,cy.a,cy.r);
					r1.calc = calc;
					r1.reflectionC = cy.reflectionC;
					
					if (r1.t > 0) {
						traceResults.push_back(r1);
					}
					r1.t = ((-b) - sqrt((b*b) - (4.0 * a*c))) / (2.0*a);

					calc.x = (fire.origin.x + r1.t*fire.direction.x);
					calc.y = (fire.origin.y + r1.t*fire.direction.y);
					calc.z = (fire.origin.z + r1.t*fire.direction.z);

					r1.normal = cylinderNormal(calc, cy.c, cy.a, cy.r);
					r1.calc = calc;
					r1.reflectionC = cy.reflectionC;

					if (r1.t > 0) {
						traceResults.push_back(r1);
					}
				}
				else if (((b*b) - (4.0 * a*c)) == 0) {
					rayTraceOut r1;
					r1.t = (-b) / (2.0*a);
					r1.ambient = cy.ambient;
					r1.diffuse = cy.diffuse;
					r1.diffuse = cy.specular;

					Pt calc;
					calc.x = (fire.origin.x + r1.t*fire.direction.x);
					calc.y = (fire.origin.y + r1.t*fire.direction.y);
					calc.z = (fire.origin.z + r1.t*fire.direction.z);

					r1.normal = cylinderNormal(calc, cy.c, cy.a, cy.r);
					r1.calc = calc;
					r1.reflectionC = cy.reflectionC;

					if (r1.t > 0) {
						traceResults.push_back(r1);
					}
				}

			}
			for (elipse s : elipses) {

				ray alt_fire = fire;
				alt_fire.origin.x = alt_fire.origin.x*s.def.mat[0][0] + s.def.mat[0][3];
				alt_fire.origin.y = alt_fire.origin.y*s.def.mat[1][1] + s.def.mat[1][3];
				alt_fire.origin.z = alt_fire.origin.z*s.def.mat[2][2] + s.def.mat[2][3];

				alt_fire.direction.x = alt_fire.direction.x*s.def.mat[0][0];
				alt_fire.direction.y = alt_fire.direction.y*s.def.mat[1][1];
				alt_fire.direction.z = alt_fire.direction.z*s.def.mat[2][2];
				alt_fire.direction = normalize(alt_fire.direction);

				float a, b, c;

				a = dot_product(alt_fire.direction, alt_fire.direction);

				Vector v = alt_fire.direction;
				v.x = v.x * 2.0;
				v.y = v.y * 2.0;
				v.z = v.z * 2.0;

				Vector l;
				l.x = alt_fire.origin.x - s.origin.x;
				l.y = alt_fire.origin.y - s.origin.y;
				l.z = alt_fire.origin.z - s.origin.z;

				b = dot_product(v, l);
				c = dot_product(l, l) - (s.radius*s.radius);

				//cout << a << "     "  << b << "     " << c << endl;
				if (((b*b) - (4.0*a*c)) > 0) {
					rayTraceOut r1;
					r1.t = ((-b) + sqrt((b*b) - (4.0 * a*c))) / (2.0*a);
					r1.ambient = s.ambient;
					r1.diffuse = s.diffuse;
					r1.specular = s.specular;

					Pt calc;
					calc.x = (alt_fire.origin.x + r1.t*alt_fire.direction.x);
					calc.y = (alt_fire.origin.y + r1.t*alt_fire.direction.y);
					calc.z = (alt_fire.origin.z + r1.t*alt_fire.direction.z);


					r1.normal = sphereNormal(s.origin, calc);
					
					r1.normal.x = r1.normal.x*s.def.mat[0][0];
					r1.normal.y = r1.normal.y*s.def.mat[1][1];
					r1.normal.z = r1.normal.z*s.def.mat[2][2];

					r1.normal = normalize(r1.normal);
					
					calc.x = (fire.origin.x + r1.t*fire.direction.x);
					calc.y = (fire.origin.y + r1.t*fire.direction.y);
					calc.z = (fire.origin.z + r1.t*fire.direction.z);

					r1.calc = calc;
					r1.reflectionC = s.reflectionC;
					if (r1.t > 0) {
						traceResults.push_back(r1);
					}
					r1.t = ((-b) - sqrt((b*b) - (4.0 * a*c))) / (2.0*a);

					calc.x = (alt_fire.origin.x + r1.t*alt_fire.direction.x);
					calc.y = (alt_fire.origin.y + r1.t*alt_fire.direction.y);
					calc.z = (alt_fire.origin.z + r1.t*alt_fire.direction.z);

					r1.normal = sphereNormal(s.origin, calc);
					r1.normal.x = r1.normal.x*s.def.mat[0][0];
					r1.normal.y = r1.normal.y*s.def.mat[1][1];
					r1.normal.z = r1.normal.z*s.def.mat[2][2];
					r1.normal = normalize(r1.normal);
					
					calc.x = (fire.origin.x + r1.t*fire.direction.x);
					calc.y = (fire.origin.y + r1.t*fire.direction.y);
					calc.z = (fire.origin.z + r1.t*fire.direction.z);
					
					r1.calc = calc;
					r1.reflectionC = s.reflectionC;
					if (r1.t > 0) {
						traceResults.push_back(r1);
					}
				}
				else if (((b*b) - (4.0 * a*c)) == 0) {
					rayTraceOut r1;
					r1.t = (-b) / (2.0*a);
					r1.ambient = s.ambient;
					r1.diffuse = s.diffuse;
					r1.diffuse = s.specular;

					Pt calc;
					calc.x = (alt_fire.origin.x + r1.t*alt_fire.direction.x);
					calc.y = (alt_fire.origin.y + r1.t*alt_fire.direction.y);
					calc.z = (alt_fire.origin.z + r1.t*alt_fire.direction.z);

					r1.normal = sphereNormal(s.origin, calc);

					r1.normal.x = r1.normal.x*s.def.mat[0][0];
					r1.normal.y = r1.normal.y*s.def.mat[1][1];
					r1.normal.z = r1.normal.z*s.def.mat[2][2];

					r1.normal = normalize(r1.normal);
	
					calc.x = (fire.origin.x + r1.t*fire.direction.x);
					calc.y = (fire.origin.y + r1.t*fire.direction.y);
					calc.z = (fire.origin.z + r1.t*fire.direction.z);

					r1.calc = calc;
					r1.reflectionC = s.reflectionC;

					if (r1.t > 0) {
						traceResults.push_back(r1);
					}
				}
			}
			sort(traceResults.begin(), traceResults.end());
			if (traceResults.size() > 0) {
				//cout << "-------" << endl;
				//cout << t.size() << endl;
				float t1 = traceResults[0].t;

				Pt calc = traceResults[0].calc;

				Vector normal = traceResults[0].normal;
				
				Vector eye;
				eye.x = viewer.x - calc.x;
				eye.y = viewer.y - calc.y;
				eye.z = viewer.z - calc.z;
				eye = normalize(eye);

				Vector color = find_color(normal, eye, traceResults[0].ambient,
					traceResults[0].diffuse,
					traceResults[0].specular, calc);
				//color = Vector(1,1,1);

				if (traceResults[0].reflectionC > 0) {
					Vector reflectedColor = reflectionColor(traceResults[0].calc, eye, traceResults[0].normal, 0);
					color.x += reflectedColor.x*traceResults[0].reflectionC;
					color.y += reflectedColor.y*traceResults[0].reflectionC;
					color.z += reflectedColor.z*traceResults[0].reflectionC;
				}

				setFramebuffer(windowx, windowy, calc.z, color.x, color.y, color.z);
			}
		}
	}
}

void display(void)
{
	rayTrace();
	drawit();
}

void keyboard(unsigned char key, int x, int y)
{
	spheres.clear();
	cylinders.clear();
	elipses.clear();
	lightSources.clear();
	planes.clear();
	clearFramebuffer();

	if (key == '1') {
		sphere stan = sphere(Pt(1,1,2), 2);
		plane steve = plane(Pt(0, -7, 3), Vector(0, 1, 0));
		elipse schafer = elipse(Pt(-1, -1, 3), 2, Vector(2,1,1));
		schafer.ambient = Vector(0.3,0,0.3);
		schafer.diffuse = Vector(0.7,0,0.7);
		steve.ambient = Vector(0,0.4,0);
		steve.diffuse = Vector(0, 0.7, 0);

		lightSources.push_back(light(Pt(10,10,-10), Vector(0.5, 0.5, 0.5)));
		lightSources.push_back(light(Pt(-10,10,-10), Vector(0.5,0.5,0.5)));
		spheres.push_back(stan);
		planes.push_back(steve);
		elipses.push_back(schafer);
	}
	if (key == '2') {

		lightSources.push_back(light(Pt(5, 5, -10), Vector(0.5, 0.5, 0.5)));
		lightSources.push_back(light(Pt(-10, 0, -10), Vector(0.5, 0.5, 0.5)));
		sphere stan;
		//stan.origin = Pt(0.4, 0.7, 2);
		//stan.radius = 3;
		//spheres.push_back(stan);

		stan.origin = Pt(-5.5, 1, 5);
		stan.radius = 3;
		stan.diffuse = Vector(0, 0.7, 0);
		stan.ambient = Vector(0, 0.6, 0);
		//spheres.push_back(stan);
		
		stan.origin = Pt(0.7, 0, -1);
		stan.radius = 0.5;
		stan.ambient = Vector(0.6, 0, 0);
		stan.diffuse = Vector(0.7, 0, 0);
		spheres.push_back(stan);

		plane steve;
		steve.planePoint = Pt(0,0,10);
		planes.push_back(steve);

		cylinder sebastian;
		sebastian.c = Pt(0, 0, 2);
		sebastian.a = Vector(0.3,1,0);
		cylinders.push_back(sebastian);

	}
	if (key == '3') {
		sphere stan = sphere(Pt(4,0,4), 3);
		stan.ambient = Vector(0.05, 0.05, 0.05);
		stan.diffuse = Vector(0.1, 0.1, 0.1);
		stan.reflectionC = 1;
		spheres.push_back(stan);
		spheres.push_back(sphere(Pt(-4,0,4), 3));

		plane steve = plane(Pt(0,-10,3),Vector(0,1,0.5));
		steve.ambient = Vector(0.3,0,0);
		steve.diffuse = Vector(0.7, 0, 0);
		
		lightSources.push_back(light(Pt(0, 10, -10), Vector(1, 1, 1)));
		planes.push_back(steve);
	}
	if (key == '5') {
		
		plane steve = plane(Pt(0, -4, 10), Vector(0, 0, -1));
		//planes.push_back(steve);
		steve = plane(Pt(0,-5,5), Vector(0,1,0));
		steve.ambient = Vector(0,0.3,0);
		steve.diffuse = Vector(0,0.6,0);
		//planes.push_back(steve);

		cylinder sebastian = cylinder(Vector(0, 1, 0), Pt(11, 0, 5), 1);
		sebastian.ambient = Vector(0.3, 0, 0.3);
		sebastian.diffuse = Vector(0.7, 0, 0.7);
		cylinders.push_back(sebastian);

		sebastian.c = Pt(-11,0,5);
		sebastian.ambient = Vector(0.3, 0, 0);
		sebastian.diffuse = Vector(0.7, 0, 0);
		cylinders.push_back(sebastian);

		sebastian.c = Pt(0, 0, 5);
		sebastian.ambient = Vector(0, 0.3, 0);
		sebastian.diffuse = Vector(0, 0.7, 0);
		cylinders.push_back(sebastian);

		sebastian.c = Pt(3.6666666, 0, 5);
		sebastian.ambient = Vector(0, 0, 0.3);
		sebastian.diffuse = Vector(0, 0, 0.7);
		cylinders.push_back(sebastian);

		sebastian.c = Pt(7.3333333, 0, 5);
		sebastian.ambient = Vector(0, 0.3, 0.3);
		sebastian.diffuse = Vector(0, 0.7, 0.7);
		cylinders.push_back(sebastian);

		sebastian.c = Pt(-3.6666666, 0, 5);
		sebastian.ambient = Vector(0.3, 0.3, 0);
		sebastian.diffuse = Vector(0.7, 0.7, 0);
		cylinders.push_back(sebastian);

		sebastian.c = Pt(-7.3333333, 0, 5);
		sebastian.ambient = Vector(0.3, 0.15, 0);
		sebastian.diffuse = Vector(0.7, 0.35, 0);
		cylinders.push_back(sebastian);

		lightSources.push_back(light(Pt(10, 10, -10), Vector(0.6, 0.6, 0.6)));
		lightSources.push_back(light(Pt(-10, 10, -10), Vector(0.6, 0.6, 0.6)));
	}
	if (key == '4') {
		plane steve = plane(Pt(0, -7, 3), Vector(0, 1, 0));
		steve.ambient = Vector(0, 0.4, 0);
		steve.diffuse = Vector(0, 0.7, 0);
		planes.push_back(steve);

		steve = plane(Pt(0, 0, 7), Vector(0, 0, -1));
		steve.ambient = Vector(0, 0.3, 0.4);
		steve.diffuse = Vector(0, 0.5, 0.7);
		
		sphere stan = sphere(Pt(0,2,7), 4);
		stan.ambient = Vector(0.1,0.1,0.1);
		stan.diffuse = Vector(0.2, 0.2, 0.2);
		stan.reflectionC = 1;
		spheres.push_back(stan);

		cylinder sebastian = cylinder(Vector(0,1,0),Pt(-4,0,0),1);
		sebastian.ambient = Vector(0.2, 0.1, 0);
		sebastian.diffuse = Vector(0.5, 0.25,0);
		cylinders.push_back(sebastian);

		sebastian = cylinder(Vector(0, 1, 0), Pt(4, 0, 0), 1);
		sebastian.ambient = Vector(0.2, 0.1, 0);
		sebastian.diffuse = Vector(0.5, 0.25, 0);
		cylinders.push_back(sebastian);

		lightSources.push_back(light(Pt(0,0,-10), Vector(1, 1, 1)));

		planes.push_back(steve);
	}
	glutPostRedisplay();
}

void init(void)
{
	clearFramebuffer();
	viewer = Pt(0,0,-2.5);
	ambient = Vector(0.5, 0.5, 0.5);

}

int main(int argc, char** argv)
{
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
	glutInitWindowSize(ImageW, ImageH);
	glutInitWindowPosition(100, 100);
	glutCreateWindow("Tunc Gocay - Homework 5");
	init();

	glutKeyboardFunc(keyboard);
	glutDisplayFunc(display);

	glutMainLoop();
	return 0;
}
