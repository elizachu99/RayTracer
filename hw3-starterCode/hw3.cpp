/* **************************
 * CSCI 420
 * Assignment 3 Raytracer
 * Name: Elizabeth Chu
 * *************************
*/

#ifdef WIN32
  #include <windows.h>
#endif

#if defined(WIN32) || defined(linux)
  #include <GL/gl.h>
  #include <GL/glut.h>
#elif defined(__APPLE__)
  #include <OpenGL/gl.h>
  #include <GLUT/glut.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef WIN32
  #define strcasecmp _stricmp
#endif

#include <imageIO.h>
#include <math.h>
#include <glm/glm.hpp>
#include <algorithm>
#include <limits>

#define MAX_TRIANGLES 20000
#define MAX_SPHERES 100
#define MAX_LIGHTS 100

char * filename = NULL;

//different display modes
#define MODE_DISPLAY 1
#define MODE_JPEG 2

int mode = MODE_DISPLAY;

//you may want to make these smaller for debugging purposes
#define WIDTH 640
#define HEIGHT 480

//the field of view of the camera
#define fov 60.0

// other global parameters
#define PI 3.14159265
#define EPSILON 0.001
/********* NOTE: comment out line below on SHADOW_OFFSET to disable soft shadows */
#define SHADOW_OFFSET 0.2  // recommended: 0.02 for close lights, 0.2 for far lights
typedef enum { SPHERE, TRIANGLE } OBJECT_TYPE;

unsigned char buffer[HEIGHT][WIDTH][3];

struct Vertex
{
  double position[3];
  double color_diffuse[3];
  double color_specular[3];
  double normal[3];
  double shininess;
};

struct Triangle
{
  Vertex v[3];
  glm::vec3 v1, v2, v3;
  glm::vec3 normal;
  double alpha, beta, gamma;
  // alpha=v1, beta=v2, gamma=v3
};

struct Sphere
{
  double position[3];
  double color_diffuse[3];
  double color_specular[3];
  double shininess;
  double radius;
};

struct Light
{
  double position[3];
  double color[3];
};

Triangle triangles[MAX_TRIANGLES];
Sphere spheres[MAX_SPHERES];
Light lights[MAX_LIGHTS];
double ambient_light[3];

int num_triangles = 0;
int num_spheres = 0;
int num_lights = 0;

void plot_pixel_display(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plot_pixel_jpeg(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plot_pixel(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void triangleCalculations();
bool intersectSpheres(glm::vec3 &origin, glm::vec3 &dir, double &t, int &hit, bool quickReturn);
bool intersectTriangles(glm::vec3 &origin, glm::vec3 &dir, double &t, int &hit, bool quickReturn);
bool solveQuadratic(double &a, double &b, double &c, double &t0, double &t1);
bool barycentric(glm::vec3 &p, Triangle &tri);
void traceShadow(glm::vec3 &point, glm::vec3 &N, glm::vec3 &V,
		 int hit, OBJECT_TYPE obj, glm::vec3 &pcolor);
void computeLight(glm::vec3 &light, glm::vec3 &point, glm::vec3 &N,
		  glm::vec3 &V, glm::vec3 &kd, glm::vec3 &ks,
		  double &shiny, glm::vec3 &lcolor, glm::vec3 &pcolor);

/* Given light position and a point on an
   object, determine the color of the point:
   if in shadow, don't update. Else, use 
   Phong shading */
void computeLight(glm::vec3 &light, glm::vec3 &point, glm::vec3 &N,
		  glm::vec3 &V, glm::vec3 &kd, glm::vec3 &ks,
		  double &shiny, glm::vec3 &lcolor, glm::vec3 &pcolor) {
  double distance, dotLN, dotRVtoS;
  int NOT_USED;
  glm::vec3 L, R, RGB;
  
  // compute unit vector L from POINT to light source
  L = light - point;
  distance = glm::length(L);
  L = glm::normalize(L);
    
  if (!intersectSpheres(point, L, distance, NOT_USED, true) &&
      !intersectTriangles(point, L, distance, NOT_USED, true)) {
    // reflected vector R = 2(L.N)N - L
    dotLN = glm::dot(L, N);
    R = glm::vec3(2 * dotLN) * N - L;
    // clamp the L.N and R.V to 0 if negative
    if (dotLN < 0) dotLN = 0;
    if (glm::dot(R, V) < 0) dotRVtoS = 0;
    else dotRVtoS = pow(glm::dot(R, V), shiny);
    // Phong illumination
    RGB.r = lcolor.r * (kd.r * dotLN + ks.r * dotRVtoS);
    RGB.g = lcolor.g * (kd.g * dotLN + ks.g * dotRVtoS);
    RGB.b = lcolor.b * (kd.b * dotLN + ks.b * dotRVtoS);
    pcolor = pcolor + RGB;
  }
}

/* Shoot shadow rays originating from POINT to
   each light source. If unobstructed, update
   COLOR to add the Phong component of the
   light source. If SHADOW_OFFSET is defined, 
   change each light source to area light  */
void traceShadow(glm::vec3 &point, glm::vec3 &N, glm::vec3 &V,
		 int hit, OBJECT_TYPE obj, glm::vec3 &pcolor) {
  glm::vec3 kd, ks, light, lcolor, R;
  double shiny, offset;
  int NOT_USED;

  if (obj == SPHERE) {
    kd.r = spheres[hit].color_diffuse[0];
    kd.g = spheres[hit].color_diffuse[1];
    kd.b = spheres[hit].color_diffuse[2];
    ks.r = spheres[hit].color_specular[0];
    ks.g = spheres[hit].color_specular[1];
    ks.b = spheres[hit].color_specular[2];
    shiny = spheres[hit].shininess;
  } else if (obj == TRIANGLE) {
    Triangle *tri = &triangles[hit];
    kd.r = (tri->alpha * tri->v[0].color_diffuse[0] +
	    tri->beta  * tri->v[1].color_diffuse[0] +
	    tri->gamma * tri->v[2].color_diffuse[0]);
    kd.g = (tri->alpha * tri->v[0].color_diffuse[1] +
	    tri->beta  * tri->v[1].color_diffuse[1] +
	    tri->gamma * tri->v[2].color_diffuse[1]);
    kd.b = (tri->alpha * tri->v[0].color_diffuse[2] +
	    tri->beta  * tri->v[1].color_diffuse[2] +
	    tri->gamma * tri->v[2].color_diffuse[2]);
    ks.r = (tri->alpha * tri->v[0].color_specular[0] +
	    tri->beta  * tri->v[1].color_specular[0] +
	    tri->gamma * tri->v[2].color_specular[0]);
    ks.g = (tri->alpha * tri->v[0].color_specular[1] +
	    tri->beta  * tri->v[1].color_specular[1] +
	    tri->gamma * tri->v[2].color_specular[1]);
    ks.b = (tri->alpha * tri->v[0].color_specular[2] +
	    tri->beta  * tri->v[1].color_specular[2] +
	    tri->gamma * tri->v[2].color_specular[2]);
    shiny = (tri->alpha * tri->v[0].shininess +
	     tri->beta  * tri->v[1].shininess +
	     tri->gamma * tri->v[2].shininess);	     
  }

  for (int i=0; i<num_lights; i++) {
    // soft shadow: add 24 lights around the original light
    light.x = lights[i].position[0];
    light.y = lights[i].position[1];
    light.z = lights[i].position[2];
    #ifdef SHADOW_OFFSET
    lcolor.r = lights[i].color[0]/25.0;
    lcolor.g = lights[i].color[1]/25.0;
    lcolor.b = lights[i].color[2]/25.0;
    #else
    lcolor.r = lights[i].color[0];
    lcolor.g = lights[i].color[1];
    lcolor.b = lights[i].color[2];
    #endif
    computeLight(light, point, N, V, kd, ks, shiny, lcolor, pcolor);
    #ifdef SHADOW_OFFSET
    offset = SHADOW_OFFSET;
    for (int j=0; j<3; j++) {
      light.x = lights[i].position[0] + offset;
      computeLight(light, point, N, V, kd, ks, shiny, lcolor, pcolor);
      light.y = lights[i].position[1] + offset;
      computeLight(light, point, N, V, kd, ks, shiny, lcolor, pcolor);
      light.x = lights[i].position[0];
      computeLight(light, point, N, V, kd, ks, shiny, lcolor, pcolor);
      light.x = lights[i].position[0] - offset;
      computeLight(light, point, N, V, kd, ks, shiny, lcolor, pcolor);
      light.y = lights[i].position[1];
      computeLight(light, point, N, V, kd, ks, shiny, lcolor, pcolor);
      light.y = lights[i].position[1] - offset;
      computeLight(light, point, N, V, kd, ks, shiny, lcolor, pcolor);
      light.x = lights[i].position[0];
      computeLight(light, point, N, V, kd, ks, shiny, lcolor, pcolor);
      light.x = lights[i].position[0] + offset;
      computeLight(light, point, N, V, kd, ks, shiny, lcolor, pcolor);
      offset += SHADOW_OFFSET;
    }
    #endif
  }
  pcolor = glm::clamp(pcolor, glm::vec3(0.0), glm::vec3(1.0));
}

/* Given intersection point and triangle, find the barycentric
   coordinates using Christer Ericson's method in his book 
   Real-Time Collision Detection */
bool barycentric(glm::vec3 &p, Triangle &tri) {
  // compute vectors
  glm::vec3 v12 = tri.v3 - tri.v1; // v0
  glm::vec3 v13 = tri.v2 - tri.v1; // v1
  glm::vec3 v1p = p - tri.v1;      // v2

  // dot products
  double dot00 = glm::dot(v12, v12);
  double dot01 = glm::dot(v12, v13);
  double dot02 = glm::dot(v12, v1p);
  double dot11 = glm::dot(v13, v13);
  double dot12 = glm::dot(v13, v1p);

  // barycentric coordinates
  double invDenom = 1 / (dot00 * dot11 - dot01 * dot01);
  double u = (dot11 * dot02 - dot01 * dot12) * invDenom;
  double v = (dot00 * dot12 - dot01 * dot02) * invDenom;

  // check if it is inside triangle
  if (u >= 0 && v >= 0 && u+v < 1) {
    tri.gamma = u;
    tri.beta = v;
    tri.alpha = 1 - u - v;
    return true;
  }
  return false;
}

/* Given a ray that goes from ORIGIN and a normalized
   direction vector DIR, compute if this ray intersects
   with any triangle. If so, find nearest intersection
   point, if the nearest is closer than nearest sphere. */
bool intersectTriangles(glm::vec3 &origin, glm::vec3 &dir,
			double &t, int &hit, bool isShadow) {
  double NdotDir, t1;
  double nearestT = std::numeric_limits<double>::max();
  glm::vec3 point;
  
  for (int i=0; i<num_triangles; i++) {
    // if ray parallel to triangle, no intersection
    NdotDir = glm::dot(triangles[i].normal, dir);
    if (std::abs(NdotDir) < EPSILON) continue;
    
    // compute t 
    t1 = -1 * glm::dot(triangles[i].normal, origin - triangles[i].v1) / NdotDir;
    if (t1 < EPSILON) continue; // triangle behind ray or is origin
    
    // intersection point
    point = origin + glm::vec3(t1) * dir;
    
    // Barycentric coordinates
    if (t1 < nearestT && barycentric(point, triangles[i])) {
      // intersection point is inside triangle
      if (isShadow && t1 < t - EPSILON) return true;
      nearestT = t1;
      hit = i;
    }
  } // end for each triangle
  if (nearestT < std::numeric_limits<double>::max() && nearestT < t) {
    t = nearestT;
    return true;
  }
  return false;
}

/* Given a, b, c, solve quadratic equation and store
   results in t0, t1. If no real solutions exist, 
   return false. */
bool solveQuadratic(double &a, double &b, double &c,
		    double &x0, double &x1) {
  double b24ac = b*b - 4*a*c;
  if (b24ac < 0) return false;
  else if (b24ac == 0) {
    x0 = x1 = (-b) / (2*a);
    return true;
  } else {
    x0 = (-b + sqrt(b24ac)) / (2*a);
    x1 = (-b - sqrt(b24ac)) / (2*a);
    return true;
  }
}

/* Given a ray that goes from ORIGIN and a normalized 
   direction vector DIR, compute if this ray intersects 
   with any sphere; i.e., if there is a 't' that 
   satisfies equation p0 + dt = point on sphere, and 
   t > 0 and t > epsilon. If so, find the nearest t 
   and return true. If isShadow is true, return 
   true immediately once any intersection that is
   between the origin and light source is found. */
bool intersectSpheres(glm::vec3 &origin, glm::vec3 &dir,
		      double &t, int &hit, bool isShadow) {
  double nearestT = std::numeric_limits<double>::max();
  for (int i = 0; i<num_spheres; i++) {
    glm::vec3 center(spheres[i].position[0],
		     spheres[i].position[1],
		     spheres[i].position[2]);
    double t0, t1;
    double a = 1.0;
    double b = 2 * glm::dot(dir, (origin-center));
    double c = (pow(origin.x-center.x, 2) +
		pow(origin.y-center.y, 2) +
		pow(origin.z-center.z, 2) -
		pow(spheres[i].radius, 2));
    if (solveQuadratic(a, b, c, t0, t1)) {
      if (t0 > t1) std::swap(t0, t1);
      if (t0 < EPSILON) {
	if (t1 >= EPSILON && t1 < nearestT) {
	  // intersected object should be between light & origin
	  if (isShadow && t1 < t - EPSILON) return true; 
	  nearestT = t1;
	  hit = i;
	}
      } else if (t0 < nearestT) {
	if (isShadow && t0 < t - EPSILON) return true;
	nearestT = t0;
	hit = i;
      }
    }
  }
  if (nearestT < std::numeric_limits<double>::max()) {
    t = nearestT;
    return true;
  }
  return false;
}

//MODIFY THIS FUNCTION
void draw_scene()
{
  double aspectRatio = (double)WIDTH / (double)HEIGHT;
  double scale = tan(fov/2 * PI / 180);
  glm::vec3 origin(0.0f);
  glm::vec3 ambient(ambient_light[0], ambient_light[1], ambient_light[2]);
  
  // ray tracing implementation
  for(unsigned int x=0; x<WIDTH; x++) {
    glPointSize(2.0);  
    glBegin(GL_POINTS);
    for(unsigned int y=0; y<HEIGHT; y++) {
      double camX, camY, rad, t = std::numeric_limits<double>::max();
      int hitIndex;
      glm::vec3 pixelColor, dir, normal, viewvec, pt;
      // convert raster coordinates to camera coordinates
      camX = (2*(x+0.5)/(double)WIDTH - 1) * aspectRatio * scale;
      camY = -(1 - 2*(y+0.5)/(double)HEIGHT) * scale;
      // trace ray
      dir = glm::normalize(glm::vec3(camX, camY, -1));
      pixelColor = glm::vec3(1.0f);
      if (intersectSpheres(origin, dir, t, hitIndex, false)) {
	// get intersect pt & normal, then shoot shadow rays
	rad = spheres[hitIndex].radius;
	pt = glm::vec3(t) * dir;
	normal.x = (pt.x - spheres[hitIndex].position[0]) / rad;
	normal.y = (pt.y - spheres[hitIndex].position[1]) / rad;
	normal.z = (pt.z - spheres[hitIndex].position[2]) / rad;
	viewvec = glm::normalize(origin - pt);
	pixelColor = ambient;
	traceShadow(pt, normal, viewvec, hitIndex, SPHERE, pixelColor);
      }
      if (intersectTriangles(origin, dir, t, hitIndex, false)) {
	// get intersect pt & normal, then shoot shadow rays
	pt = glm::vec3(t) * dir;
	viewvec = glm::normalize(origin - pt);
	pixelColor = ambient;
	// interpolate normals using barycentric coordinates
	Triangle *tri = &triangles[hitIndex];
	normal.x = (tri->alpha * tri->v[0].normal[0] +
		    tri->beta * tri->v[1].normal[0] +
		    tri->gamma * tri->v[2].normal[0]);
	normal.y = (tri->alpha * tri->v[0].normal[1] +
		    tri->beta * tri->v[1].normal[1] +
		    tri->gamma * tri->v[2].normal[1]);
	normal.z = (tri->alpha * tri->v[0].normal[2] +
		    tri->beta * tri->v[1].normal[2] +
		    tri->gamma * tri->v[2].normal[2]);
	normal = glm::normalize(normal);
	traceShadow(pt, normal, viewvec, hitIndex, TRIANGLE, pixelColor);
      }

      // plot pixel color
      plot_pixel(x, y, (int)(255*pixelColor.r),
		 (int)(255*pixelColor.g), (int)(255*pixelColor.b));
    }
    glEnd();
    glFlush();
  }
  printf("Done!\n"); fflush(stdout);
}

void plot_pixel_display(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
  glColor3f(((float)r) / 255.0f, ((float)g) / 255.0f, ((float)b) / 255.0f);
  glVertex2i(x,y);
}

void plot_pixel_jpeg(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
  buffer[y][x][0] = r;
  buffer[y][x][1] = g;
  buffer[y][x][2] = b;
}

void plot_pixel(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
  plot_pixel_display(x,y,r,g,b);
  if(mode == MODE_JPEG)
    plot_pixel_jpeg(x,y,r,g,b);
}

void save_jpg()
{
  printf("Saving JPEG file: %s\n", filename);

  ImageIO img(WIDTH, HEIGHT, 3, &buffer[0][0][0]);
  if (img.save(filename, ImageIO::FORMAT_JPEG) != ImageIO::OK)
    printf("Error in Saving\n");
  else 
    printf("File saved Successfully\n");
}

/* For each triangle in scene, calculate their normals 
   and v1, v2, v3 as glm::vec3's. */
void triangleCalculations() {
  glm::vec3 v12, v13;
  for (int i=0; i<num_triangles; i++) {
    triangles[i].v1 = glm::vec3(triangles[i].v[0].position[0],
		      triangles[i].v[0].position[1],triangles[i].v[0].position[2]);
    triangles[i].v2 = glm::vec3(triangles[i].v[1].position[0],
		      triangles[i].v[1].position[1],triangles[i].v[1].position[2]);
    triangles[i].v3 = glm::vec3(triangles[i].v[2].position[0],
		      triangles[i].v[2].position[1],triangles[i].v[2].position[2]);
    v12 = triangles[i].v2 - triangles[i].v1;
    v13 = triangles[i].v3 - triangles[i].v1;
    triangles[i].normal = glm::normalize(glm::cross(v12, v13));
  }
}

void parse_check(const char *expected, char *found)
{
  if(strcasecmp(expected,found))
  {
    printf("Expected '%s ' found '%s '\n", expected, found);
    printf("Parse error, abnormal abortion\n");
    exit(0);
  }
}

void parse_doubles(FILE* file, const char *check, double p[3])
{
  char str[100];
  fscanf(file,"%s",str);
  parse_check(check,str);
  fscanf(file,"%lf %lf %lf",&p[0],&p[1],&p[2]);
  printf("%s %lf %lf %lf\n",check,p[0],p[1],p[2]);
}

void parse_rad(FILE *file, double *r)
{
  char str[100];
  fscanf(file,"%s",str);
  parse_check("rad:",str);
  fscanf(file,"%lf",r);
  printf("rad: %f\n",*r);
}

void parse_shi(FILE *file, double *shi)
{
  char s[100];
  fscanf(file,"%s",s);
  parse_check("shi:",s);
  fscanf(file,"%lf",shi);
  printf("shi: %f\n",*shi);
}

int loadScene(char *argv)
{
  FILE * file = fopen(argv,"r");
  int number_of_objects;
  char type[50];
  Triangle t;
  Sphere s;
  Light l;
  fscanf(file,"%i", &number_of_objects);

  printf("number of objects: %i\n",number_of_objects);

  parse_doubles(file,"amb:",ambient_light);

  for(int i=0; i<number_of_objects; i++)
  {
    fscanf(file,"%s\n",type);
    printf("%s\n",type);
    if(strcasecmp(type,"triangle")==0)
    {
      printf("found triangle\n");
      for(int j=0;j < 3;j++)
      {
        parse_doubles(file,"pos:",t.v[j].position);
        parse_doubles(file,"nor:",t.v[j].normal);
        parse_doubles(file,"dif:",t.v[j].color_diffuse);
        parse_doubles(file,"spe:",t.v[j].color_specular);
        parse_shi(file,&t.v[j].shininess);
      }

      if(num_triangles == MAX_TRIANGLES)
      {
        printf("too many triangles, you should increase MAX_TRIANGLES!\n");
        exit(0);
      }
      triangles[num_triangles++] = t;
    }
    else if(strcasecmp(type,"sphere")==0)
    {
      printf("found sphere\n");

      parse_doubles(file,"pos:",s.position);
      parse_rad(file,&s.radius);
      parse_doubles(file,"dif:",s.color_diffuse);
      parse_doubles(file,"spe:",s.color_specular);
      parse_shi(file,&s.shininess);

      if(num_spheres == MAX_SPHERES)
      {
        printf("too many spheres, you should increase MAX_SPHERES!\n");
        exit(0);
      }
      spheres[num_spheres++] = s;
    }
    else if(strcasecmp(type,"light")==0)
    {
      printf("found light\n");
      parse_doubles(file,"pos:",l.position);
      parse_doubles(file,"col:",l.color);

      if(num_lights == MAX_LIGHTS)
      {
        printf("too many lights, you should increase MAX_LIGHTS!\n");
        exit(0);
      }
      lights[num_lights++] = l;
    }
    else
    {
      printf("unknown type in scene description:\n%s\n",type);
      exit(0);
    }
  }
  return 0;
}

void display()
{
}

void init()
{
  glMatrixMode(GL_PROJECTION);
  //glOrtho(left, right, bottom, top, near, far)
  glOrtho(0,WIDTH,0,HEIGHT,1,-1);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
  glClear(GL_COLOR_BUFFER_BIT);
  triangleCalculations();
}

void idle()
{
  //hack to make it only draw once
  static int once=0;
  if(!once)
  {
    draw_scene();
    if(mode == MODE_JPEG)
      save_jpg();
  }
  once=1;
}

int main(int argc, char ** argv)
{
  if ((argc < 2) || (argc > 3))
  {  
    printf ("Usage: %s <input scenefile> [output jpegname]\n", argv[0]);
    exit(0);
  }
  if(argc == 3)
  {
    mode = MODE_JPEG;
    filename = argv[2];
  }
  else if(argc == 2)
    mode = MODE_DISPLAY;

  glutInit(&argc,argv);
  loadScene(argv[1]);

  glutInitDisplayMode(GLUT_RGBA | GLUT_SINGLE);
  glutInitWindowPosition(0,0);
  glutInitWindowSize(WIDTH,HEIGHT);
  int window = glutCreateWindow("Ray Tracer");
  glutDisplayFunc(display);
  glutIdleFunc(idle);
  init();
  glutMainLoop();
}

