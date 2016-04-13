#include "sphere.h"

#include <cmath>

#include  "../bsdf.h"
#include "../misc/sphere_drawing.h"

namespace CGL { namespace StaticScene {

bool Sphere::test(const Ray& r, double& t1, double& t2) const {

  // TODO Part 1, task 4:
  // Implement ray - sphere intersection test.
  // Return true if there are intersections and writing the
  // smaller of the two intersection times in t1 and the larger in t2.
 

  //return false;

}

bool Sphere::intersect(const Ray& r) const {

  // TODO Part 1, task 4:
  // Implement ray - sphere intersection.
  // Note that you might want to use the the Sphere::test helper here.

  
   Vector3D oc = r.o - o;
  double a = r.d[0]*r.d[0]+r.d[1]*r.d[1]+r.d[2]*r.d[2];
  double b = 2*(oc[0]*r.d[0]+oc[1]*r.d[1]+oc[2]*r.d[2]);
  double c = oc[0]*oc[0] + oc[1]*oc[1] + oc[2]*oc[2] - r2;
  double b24ac = b*b - 4*a*c;
  if(b24ac<0){
   // cout<<"222"<<endl;
    return false;
  }else{
   // cout<<"111"<<endl;
    double tt1 = (-b+sqrt(b24ac))/(2*a);
    double tt2 = (-b-sqrt(b24ac))/(2*a);
    r.min_t = min(tt1,tt2);
    r.max_t = max(tt1,tt2);
    return true;
  }
  
  //return false;

}

bool Sphere::intersect(const Ray& r, Intersection *i) const {

  // TODO Part 1m task 4:
  // Implement ray - sphere intersection.
  // Note again that you might want to use the the Sphere::test helper here.
  // When an intersection takes place, the Intersection data should be updated
  // correspondingly.


/////////////////////////////////
  double a=dot(r.d,r.d);
  double b=2*dot(r.o-this->o,r.d);
  double c=dot(r.o-this->o,r.o-this->o)-this->r2;
  double delta=b*b-4*a*c;

  if(delta>=0)
  {
    if(-b-sqrt(delta)>=0)
    {

      double t=(-b-sqrt(delta))/2/a;
      if ((t>=r.min_t)&&(t<=r.max_t))
      {
        r.max_t=t;
        i->t=t;
        i->n=this->normal(r.o+t*r.d);
        i->primitive=this;
        i->bsdf=i->primitive->get_bsdf();
        return true;
      }
    }
    else
    {

      double t=(-b+sqrt(delta))/2/a;
      if ((t>=r.min_t)&&(t<=r.max_t))
      {

        r.max_t=t;
        i->t=t;
        i->n=this->normal(r.o+t*r.d);
        i->primitive=this;
        i->bsdf=i->primitive->get_bsdf();
        return true;
      }
    }
  }

  return false;

}

void Sphere::draw(const Color& c) const {
  Misc::draw_sphere_opengl(o, r, c);
}

void Sphere::drawOutline(const Color& c) const {
    //Misc::draw_sphere_opengl(o, r, c);
}


} // namespace StaticScene
} // namespace CGL
