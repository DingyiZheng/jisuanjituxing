#include "triangle.h"

#include "CGL/CGL.h"
#include "GL/glew.h"

namespace CGL { namespace StaticScene {

Triangle::Triangle(const Mesh* mesh, size_t v1, size_t v2, size_t v3) :
    mesh(mesh), v1(v1), v2(v2), v3(v3) { }

BBox Triangle::get_bbox() const {

  Vector3D p1(mesh->positions[v1]), p2(mesh->positions[v2]), p3(mesh->positions[v3]);
  BBox bb(p1);
  bb.expand(p2); 
  bb.expand(p3);
  return bb;

}

bool Triangle::intersect(const Ray& r) const {
  
  // TODO Part 1, task 3: implement ray-triangle intersection
  Vector3D p1(mesh->positions[v1]), p2(mesh->positions[v2]), p3(mesh->positions[v3]);
  double b1;
  double b2;
  double t;
  Vector3D E1 = p2 - p1;
  Vector3D E2 = p3 - p1;
  Vector3D S = r.o - p1;
  Vector3D S1 = cross(r.d,E2);
  Vector3D S2 = cross(S,E1);
 
  if(!dot(S1,E1)){
  double k = 1.0/dot(S1,E1);
  t = k*dot(S2,E2);
  b1 = k*dot(S1,S);
  b2 = k*dot(S2,r.d);
  if(b1>=.0&&b1<=1.0&&b2>=.0&&b2<=1.0&&(b1+b2)<=1.0&&t>=.0&&t<=r.max_t){
    r.max_t = t;
    return true;

  }else{
    return false;
  }
  }
  return false;
  //return false;
}

bool Triangle::intersect(const Ray& r, Intersection *isect) const {  
  // TODO Part 1, task 3:   
  // implement ray-triangle intersection. When an intersection takes
  // place, the Intersection data should be updated accordingly
	  	  //	cout<<"!!!!!!!!!!!!!!!!"<<endl;
  Vector3D p1(mesh->positions[v1]), p2(mesh->positions[v2]), p3(mesh->positions[v3]);
  Vector3D n1(mesh->normals[v1]), n2(mesh->normals[v2]), n3(mesh->normals[v3]);
  Vector3D E1 = p2 - p1;
  Vector3D E2 = p3 - p1;
  Vector3D S = r.o - p1;
  Vector3D S1 = cross(r.d,E2);
  Vector3D S2 = cross(S,E1);
  double k = dot(S1,E1);


  if(k!=.0){
      double t = dot(S2,E2)/k;
  double b1 = dot(S1,S)/k;
  double b2 = dot(S2,r.d)/k;
  Vector3D inte=Vector3D(dot(S2,E2),dot(S1,S),dot(S2,r.d))/k;
    	  	//cout<<b1<<endl;
     if((t>=r.min_t)&&(inte[0]<=r.max_t)&&(b1>0)&&(b2>0)&&((b1+b2)<1)){
  	  	  
  	r.max_t = t;
  	isect->t = t;
  	Vector3D n = (1.0-b1-b2)*n1+b1*n2+b2*n3;
  	isect->n = n;
  	Vector3D hitpoint = r.d*t;
 	isect->primitive = this;
 	//isect->primitive = new vector<Primitive *>(prims);
  	isect->bsdf = get_bsdf();
  	  	//cout<<"!!!!!!!!!!!!!!!!"<<endl;
  	return true;


    }
  }else{
    return false;
  }  
  //return false;
}

void Triangle::draw(const Color& c) const {
  glColor4f(c.r, c.g, c.b, c.a);
  glBegin(GL_TRIANGLES);
  glVertex3d(mesh->positions[v1].x,
             mesh->positions[v1].y,
             mesh->positions[v1].z);
  glVertex3d(mesh->positions[v2].x,
             mesh->positions[v2].y,
             mesh->positions[v2].z);
  glVertex3d(mesh->positions[v3].x,
             mesh->positions[v3].y,
             mesh->positions[v3].z);
  glEnd();
}

void Triangle::drawOutline(const Color& c) const {
  glColor4f(c.r, c.g, c.b, c.a);
  glBegin(GL_LINE_LOOP);
  glVertex3d(mesh->positions[v1].x,
             mesh->positions[v1].y,
             mesh->positions[v1].z);
  glVertex3d(mesh->positions[v2].x,
             mesh->positions[v2].y,
             mesh->positions[v2].z);
  glVertex3d(mesh->positions[v3].x,
             mesh->positions[v3].y,
             mesh->positions[v3].z);
  glEnd();
}



} // namespace StaticScene
} // namespace CGL
