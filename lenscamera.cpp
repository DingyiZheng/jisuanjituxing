#include "lenscamera.h"

#include "image.h"

using namespace std;

namespace CGL {


/****** Helpers ******/
  

// Extract the R, G, or B channel out of an RGBA color stored in a single 32bit integer
static uint32_t red_channel(uint32_t color) {
    return (255 & (color >> (0)));
}

static uint32_t green_channel(uint32_t color) {
    return (255 & (color >> (8)));
}

static uint32_t blue_channel(uint32_t color) {
    return (255 & (color >> (16)));
}

// Convert from millimeters to meters
static const double scale = .001;







/****** LensElement functions ******/


bool LensElement::pass_through(Ray &r, double &prev_ior) const {
  // Part 1 Task 1: Implement this. It takes r and passes it through this lens element.
  //cout<<"before"<<r.o<<endl;
  Vector3D hit_p;
/*  if(radius == .0){
    return true;
  }*/
  //return true;
    //cout<<"r.d="<<r.d<<endl;
  if(!intersect(r,&hit_p)){
    return false;
  }else{
    if(!refract(r,hit_p,prev_ior)){
      return false;
    }else{
      if(radius!=0.)
    prev_ior=ior;
      //cout<<"after"<<r.o<<endl;
      return true;
    }
  }

}
bool LensElement::intersect(const Ray &r, Vector3D *hit_p) const {
  // Part 1 Task 1: Implement this. It intersects r with this spherical lens elemnent 
  // (or aperture diaphragm). You'll want to reuse some sphere intersect code.
  if(radius == 0){
    double t = (center - r.o.z)/r.d.z;
    Vector3D hp = r.o + t*r.d;
    double hp2z = sqrt(hp.x*hp.x + hp.y*hp.y);
    if(hp2z<0.5*aperture){
      return true;
    }
    return false;
  }
  

  Vector3D Center = Vector3D(.0,.0,center);
  Vector3D oc = r.o - Center;
  double a = r.d[0]*r.d[0]+r.d[1]*r.d[1]+r.d[2]*r.d[2];
  double b = 2*(oc[0]*r.d[0]+oc[1]*r.d[1]+oc[2]*r.d[2]);
  double c = oc[0]*oc[0] + oc[1]*oc[1] + oc[2]*oc[2] - radius*radius;
  double b24ac = b*b - 4*a*c;
  if(b24ac<=0){
   // cout<<"222"<<endl;
    return false;
  }else{
   // cout<<"111"<<endl;
    double tt1 = (-b+sqrt(b24ac))/(2*a);
    double tt2 = (-b-sqrt(b24ac))/(2*a);
   // r.min_t = min(tt1,tt2);
    //r.max_t = max(tt1,tt2);
    Vector3D hit_t1 = r.o + tt1*r.d;
    Vector3D hit_t2 = r.o + tt2*r.d;
   // double hit1_2_z = sqrt(hit_t1.x*hit_t1.x + hit_t1.y*hit_t1.y);
  //  double hit2_2_z = sqrt(hit_t2.x*hit_t2.x + hit_t2.y*hit_t2.y);
    // if(min(hit1_2_z,hit2_2_z)>=0.5*aperture){
    //   return false;
    // }


    // if(hit1_2_z<hit2_2_z){
    //    *hit_p = hit_t1;

    //    return true;
    // }else{
    //    *hit_p = hit_t2;
    //    return true;
    // }
    if(helperhit(radius,aperture,center,hit_t1)&&helperhit(radius,aperture,center,hit_t2)){
      *hit_p = helpermin(r,tt1,tt2);
     // cout<<"123"<<endl;
      return true;
    }else{
      if(helperhit(radius,aperture,center,hit_t1)){
        *hit_p = hit_t1;
        return true;
      }
      if(helperhit(radius,aperture,center,hit_t2)){
        *hit_p = hit_t2;
        return true;
      }
    }

    return false;

  }


 // return true;
  
}

bool LensElement::helperhit(double radius, double aperture,double center, Vector3D hit1) const{
  if(radius<0){
    if(hit1.z>center&&sqrt(hit1.x*hit1.x+hit1.y*hit1.y)<0.5*aperture){
      return true;

    }
    return false;
  }

  if(radius >0 ){
    if(hit1.z<center&&sqrt(hit1.x*hit1.x+hit1.y*hit1.y)<0.5*aperture){
      return true;

    }
    return false;
  }

} 

Vector3D LensElement::helpermin(Ray r, double t1,double t2) const{
  // if(v1.z<v2.z){
  //   return v1;
  // }
  return r.o + min(t1,t2)*r.d;

}
bool LensElement::refract(Ray& r, const Vector3D& hit_p, const double& prev_ior) const {
  // Part 1 Task 1: Implement this. It refracts the Ray r with this lens element or 
  // does nothing at the aperture element.
  // You'll want to consult your refract function from the previous assignment.

  if(radius == 0){
   return true;

  }


  Vector3D normal = Vector3D(hit_p.x, hit_p.y, hit_p.z - center);
  normal = normal.unit();
  //Vector3D rd = r.d.unit();
   double judge = dot(normal,-r.d);
 // double judge = normal.x*rd.x + normal.y*rd.y+ normal.z*rd.z; 
  if(judge <0 ){
    normal = -normal;
  }

  // Matrix3x3 w2o;
  // make_coord_space(w2o, normal);
  // Matrix3x3 o2w = w2o.T();

    Matrix3x3 o2w;
  make_coord_space(o2w, normal);
  Matrix3x3 w2o = o2w.T();

  //Vector3D hit_p = r.o + r.d * isect.t;
  Vector3D w_in = w2o * (-r.d);
  //Vector3D w_out;

  
  double eta = prev_ior/ior;
  if(r.d.z>0){
    eta = ior/prev_ior;
  }
  //Flip this ratio if the ray is traveling backwards.

  double sin2theta =  eta*eta*fmax(0.0, 1.0 - w_in.z * w_in.z);
  if(sin2theta > 1.0){
    return false;
  }

  double cosine = sqrt(1.0 - sin2theta);
  Vector3D newdirection = Vector3D(-eta*w_in.x, -eta*w_in.y, -cosine);

  // if(dot(newdirection,normal)>=0){
  //   return false;
  // }
  newdirection = o2w*newdirection;
  newdirection = newdirection.unit();
  


  r = Ray(hit_p,newdirection);

  return true;


//  return true;





}






/****** Lens functions ******/



void Lens::parse_lens_file(std::string filename) {

  ifstream infile(filename);
  string line;
  double z_coord = 0;
  double z_ap;
  vector<LensElement> backwards;
  elts.clear();
  bool first = true;
  while (getline(infile, line)) {
    if (first) {
      cout << "[Lens] Loading lens file " << line << endl;
      first = false;
    }
    if (line[0] == '#')
      continue;
    stringstream ss(line);
    LensElement lens;
    double offset;
    ss >> lens.radius >> offset >> lens.ior >> lens.aperture;
    lens.center = z_coord;
    if (!lens.radius) {
      z_ap = z_coord;
    }
    z_coord += offset;
    backwards.push_back(lens);
  }
  for (int i = backwards.size() - 1; i >= 0; --i) {
    LensElement l = backwards[i];
    l.center = (l.center - z_ap) + l.radius;
    if (i) l.ior = backwards[i-1].ior;
    else l.ior = 1;
    if (!l.ior) l.ior = 1;
    elts.push_back(l);
    if (!l.radius)
      ap_i = elts.size()-1;
    // cout << "Lens element edge first " << (l.center - l.radius) << " " 
    //   << l.radius << " " << l.center << " " << l.ior << " " << l.aperture << endl;
  }
  double c = elts.front().center, r = elts.front().radius, a = elts.front().aperture * .5;
  back_elt = c - (r>0?1:-1) * sqrt(r*r-a*a);
  ap_radius = ap_original = elts[ap_i].aperture;

  // Get infinity and close focus depths, also get focal length.
  set_focus_params();
  // Focus at infinity to start.
  //cout<<infinity_focus<<endl;
  sensor_depth = infinity_focus;
       
}


void Lens::set_focus_params() {

  // Part 1 Task 2: Implement this. 
  // After this function is called, the three variables
  // infinity_focus, near_focus, and focal_length
  // should be set correctly.
  Ray r1 = Ray(Vector3D(1.0,.0,-100.0),Vector3D(.0,.0,1.0));
  //Ray r2 = Ray(Vector3D(.0,.0,-100.0),Vector3D(.0,.0,1.0));
 // Ray r1 = Ray(Vector3D(1.0,1.0,-10.0),Vector3D(.0,.0,1.0));
  vector<Vector3D> trace1;
  //trace1.push_back(r1.o);
  //vector<Vector3D> trace2;
  trace_backwards(r1, &trace1);
  //cout<<"r1"<<r1.o<<r1.d<<endl;
 // trace_backwards(Ray &r2, &trace2);
  //double t = (r2.o.z - r1.o.z)/(r1.d.z - r2.d.z);
  double t;
  if(r1.d.x != .0){
    t = -(r1.o.x)/(r1.d.x);
  }else{
   t = -(r1.o.y)/(r1.d.y);
  }
  Vector3D focus_p = r1.o + t*r1.d;
  //cout<<"focus_p"<<focus_p<<endl;
  infinity_focus = focus_p.z;

  focal_length = abs(1.0/r1.d.x*r1.d.z);
 
   vector<Vector3D> trace2;
  Ray r2 = Ray(Vector3D(.0,.0,-5*focal_length),Vector3D(.01,.0,1.0).unit());
  trace_backwards(r2,&trace2);
  // if(r2.d.x != .0){
  //   t = -(r2.o.x)/(r2.d.x);
  // }else{
  //  t = -(r2.o.y)/(r2.d.y);
  // }
  // Vector3D focus_p2 = r2.o + t*r2.d;
  // //cout<<"focus_p"<<focus_p<<endl;
  // near_focus = focus_p2.z;
  near_focus = -r2.o.x/r2.d.x*r2.d.z+r2.o.z;


  /////////////////
  // double z = elts.back().center - elts.back().radius - (1 + log(focal_length)) * focal_length;
  //  Vector3D near_orig(0, 0, z);
  //  Vector3D near_dir = (Vector3D(elts.back().aperture / 100.0, 0, elts.back().center - elts.back().radius) -
  //                    near_orig).unit();
  //  Ray r2 = Ray(near_orig, near_dir);
  //  trace_backwards(r2,&trace2);
  //  near_focus = -r2.o.x/r2.d.x*r2.d.z+r2.o.z;




  cout << "[Lens] Infinity focus depth is " << infinity_focus << endl;
  cout << "[Lens] Close focus depth is " << near_focus << endl;
  cout << "[Lens] True focal length is " << focal_length << endl;
}




bool Lens::trace(Ray &r, std::vector<Vector3D> *trace) const {
  // Part 1 Task 1: Implement this. It traces a ray from the sensor out into the world.
  double prev_ior = 1;
  for (int i = 0; i <= elts.size()-1; i++){
   // prev_ior = i > 0 ? elts[i-1].ior : 1;
   if (!elts[i].pass_through(r,prev_ior))
       return false;
   // elts[i].pass_through(r,prev_ior);
    trace->push_back(r.o);
  
  } 


  return true;
}

bool Lens::trace_backwards(Ray &r, std::vector<Vector3D> *trace) const {
  // Part 1 Task 1: Implement this. It traces a ray from the world backwards through 
  // the lens towards the sensor.
    double prev_ior;
  for (int i = elts.size()-1; i >= 0; --i){
    prev_ior = i > 0 ? elts[i-1].ior : 1;
   if (!elts[i].pass_through(r,prev_ior))
       return false;
   // elts[i].pass_through(r,prev_ior);
    trace->push_back(r.o);
  
  }


  return true;
}

float Lens::focus_depth(float d) const {

  // Part 1 Task 2: Implement this. Should find the conjugate of a ray
  // starting from the sensor at depth d.
   Ray r = Ray(Vector3D(0,0,d),Vector3D(0.01,0,-1).unit());
   vector<Vector3D> trace3;
   trace(r,&trace3);
   Vector3D focus_depth = r.o+(-r.o.x/r.d.x)*r.d;

   return focus_depth.z;

  //return 0;
}

Vector3D Lens::back_lens_sample() const {

  // Part 1 Task 2: Implement this. Should return a point randomly sampled
  // on the back element of the lens (the element closest to the sensor)
  double theta = 4*M_PI_2*random_uniform();
  double r = elts[0].aperture*0.5*sqrt(random_uniform());

  return Vector3D(r*cos(theta),r* sin(theta),elts[0].center - elts[0].radius);

}



/****** LensCamera functions ******/


LensCamera::LensCamera(): pt(NULL) {
  string path = string(__FILE__).substr(0,string(__FILE__).find_last_of('/')+1) + "../lenses/";
  static const vector<string> lens_files = {"dgauss.50mm.dat", "wide.22mm.dat", "telephoto.250mm.dat", "fisheye.10mm.dat"};
  for (string lens_file : lens_files)
    lenses.emplace_back(path + lens_file);

  mount_lens(0);
}


Ray LensCamera::generate_ray(double x, double y) const {
  
  Ray r = Ray(Vector3D(),Vector3D() );
  if (lens_ind >= 0) {
//cout<<"generate_ray23333"<<endl;
    // Part 1 Task 2: Implement this. It generates a ray from sensor pixel (x,y)
    // pointing toward the back element of the lens (use back_lens_sample) and traces
    // it through the Lens (using your "trace" function)
    double film_d = sqrt(24*24+36*36);
    double screen_d = sqrt(screenW*screenW + screenH*screenH);
    double film_w = film_d * screenW / screen_d;
    double film_h = film_d * screenH / screen_d;
    Vector3D sensor_point(-(x-0.5)*film_w, -(y-0.5)*film_h, lenses[lens_ind].sensor_depth);

    Vector3D backsample;
    backsample = lenses[lens_ind].back_lens_sample();
    r.o = sensor_point;
    r.d = Vector3D(backsample - sensor_point);
    r.d = r.d.unit();


        std::vector<Vector3D> tra;
    if(!lenses[lens_ind].trace(r,&tra))
    {
      r.d=c2w*Vector3D(0,0,1.0);
      //lenses[lens_ind].trace_backwards(r,&tra);
      r.d=r.d.unit();
      return r;
    }


    /***** end of your code ******/


    // This code converts the ray you traced through the lens into world coordinates.
    r.o = pos + c2w * r.o * scale;
    r.d = (c2w * r.d).unit();

  } else {

    // Generate ray for a pinhole camera. Same as in the previous assignment.
    x = 2*(x-.5); y = 2*(y-.5);
    r = Ray(pos,(c2w*Vector3D(x*tan(radians(hFov)*.5),y*tan(radians(vFov)*.5),-1)).unit());

  }

  r.min_t = nClip; r.max_t = fClip;
  return r;
}



void LensCamera::move_sensor(float delta) {
  if (lens_ind < 0) return;
  curr_lens().sensor_depth += delta;
  cout << "[LensCamera] Sensor plane moved to " << curr_lens().sensor_depth
       << ", focus now at " << lenses[lens_ind].focus_depth(lenses[lens_ind].sensor_depth) << endl;
}

void LensCamera::stop_down(float ratio) {
  float ap = curr_lens().ap_radius * ratio;
  if (ap > curr_lens().ap_original) ap = curr_lens().ap_original;
  curr_lens().ap_radius = ap;
  cout << "[LensCamera] Aperture is now " << curr_lens().ap_radius << "mm" << endl;
}

void LensCamera::mount_lens(int i) {
  lens_ind = i;
  if (i >= 0) {
    cout << "[LensCamera] Switched to lens #" << (i+1) 
         << " with focal length " << curr_lens().focal_length << "mm" << endl;
  } else {
    cout << "[LensCamera] Switched to pinhole camera" << endl;
  }
}



// A dummy function to demonstrate how to work with the image buffer.
// Calculates the average value of the green color channel in the image.
// You'll have to remember your 2D array indexing in order to take differences
// of neighboring pixels in a more sophisticated metric function.
static double mean_green(const ImageBuffer& ib) {
  double sum = 0;
  for (int i = 0; i < ib.w * ib.h; ++i) {
      sum += green_channel(ib.data[i]);
  }
  double mean = sum / (ib.w * ib.h);
  
  return mean;
}

double LensCamera::focus_metric(const ImageBuffer& ib) const {

  // Part 2 Task 1: Implement this. Design a metric to judge how "in-focus"
  // the image patch stored in the provided ImageBuffer is.

  return mean_green(ib); //  A meaningless standin
}


void LensCamera::autofocus() {


  // Part 2 Task 2: Implement this. Design a global search using your 
  // focus metric to set the sensor to be at the depth where the 
  // render cell is most "in focus". Provided code shows how to 
  // move the sensor, request a render of the cell, and evaluate the focus metric.

  // This call ensures that your pathtracer is rendering at high enough quality.
  // Increase samples per pixel to 16 and samples per light to 16.
  pt->bump_settings();

  // Example code. Nothing to do with your actual implementation except to 
  // demonstrate functionality.
  ImageBuffer ib;
  curr_lens().sensor_depth += 1;
  pt->raytrace_cell(ib);
  cout << "[LensCamera] The mean green is " << focus_metric(ib) << endl;


  
}





void LensCamera::dump_settings(string filename) {
  ofstream file(filename);
  file << hFov << " " << vFov << " " << ar << " " << nClip << " " << fClip << endl;
  for (int i = 0; i < 3; ++i)
    file << pos[i] << " ";
  for (int i = 0; i < 3; ++i)
    file << targetPos[i] << " ";
  file << endl;
  file << phi << " " << theta << " " << r << " " << minR << " " << maxR << endl;
  for (int i = 0; i < 9; ++i)
    file << c2w(i/3, i%3) << " ";
  file << endl;
  file << screenW << " " << screenH << " " << screenDist << endl;

  file << lens_ind << endl;
  for (Lens &lens : lenses) {
    file << lens.sensor_depth << " ";
  }
  file << endl;

  cout << "[LensCamera] Dumped settings to " << filename << endl;
}

void LensCamera::load_settings(string filename) {
  ifstream file(filename);

  file >> hFov >> vFov >> ar >> nClip >> fClip;
  for (int i = 0; i < 3; ++i)
    file >> pos[i];
  for (int i = 0; i < 3; ++i)
    file >> targetPos[i];
  file >> phi >> theta >> r >> minR >> maxR;
  for (int i = 0; i < 9; ++i)
    file >> c2w(i/3, i%3);
  file >> screenW >> screenH >> screenDist;

  file >> lens_ind;
  for (Lens &lens : lenses) {
    file >> lens.sensor_depth;
  }

  cout << "[LensCamera] Loaded settings from " << filename << endl;
}


} // namespace CGL

