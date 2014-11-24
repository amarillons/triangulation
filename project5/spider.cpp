#include <cfloat>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <vector>
#include <map>

class Vertex{
public:
  double x;
  double y;
  double z;
  int index;
  const void printCoordinate(){
    printf("v: %.2lf %.2lf %.2lf\n",x,y,z);
  }
};

class Triangle {
public:
  int ind1, ind2, ind3;      // vertex indices
  int removed;
};

class Circle{
public:
  int ind1,ind2,ind3;
  double x,y,z;
  double r;
};

std::vector<Vertex> Vertices;
std::vector<Triangle> Triangles;
std::vector<Circle> Circles;

float sign(Vertex p1, Vertex p2, Vertex p3){
  return (p1.x - p3.x) * (p2.y - p3.y) - (p2.x - p3.x) * (p1.y - p3.y);
}

bool divideOnlyOnceforEachVertex;

// bool containsTheVertex(Triangle _t, Vertex _pt){

//   int b1, b2, b3;
//   Vertex v1 = _t.p1;
//   Vertex v2 = _t.p2;
//   Vertex v3 = _t.p3;  

//   b1 = sign(_pt, v1, v2) < 0.0f;
//   b2 = sign(_pt, v2, v3) < 0.0f;
//   b3 = sign(_pt, v3, v1) < 0.0f;

//   return ((b1 == b2) && (b2 == b3));
// }

void showNumberofTriangles(){

  int num = 0;
  std::vector<class Triangle>::iterator it;
  it = Triangles.begin();
  while(it!=Triangles.end()){
    num++;
    it++;
  }
  printf("# number of Triangles:%d\n",num);
}

void assignIndex(){

  std::vector<class Vertex>::iterator it;
  it = Vertices.begin();
  int num = 0;
  while(it!=Vertices.end()){
    it->index = num;
    num++;
    it++;
  }
}

void placeVertex2(){

  Vertex v1;  Triangle tr1;

  v1.x = -12.;  v1.y = 12.; v1.z = 0.; Vertices.push_back(v1);   
  v1.x = +12.;  v1.y = 12.; v1.z = 0.; Vertices.push_back(v1);   
  v1.x =   0.;  v1.y =-12.; v1.z = 0.; Vertices.push_back(v1);  
  tr1.ind1 = 0;   tr1.ind2 = 1;   tr1.ind3 = 2; tr1.removed = 0;
  Triangles.push_back(tr1);

  v1.x = 0.; v1.y = 2; v1.z = 0.; Vertices.push_back(v1);
  v1.x = 4.; v1.y = 0.1; v1.z = 0.; Vertices.push_back(v1);
  v1.x = -4.; v1.y = -1.; v1.z = 0.; Vertices.push_back(v1);
  v1.x = 0.; v1.y = 5.; v1.z = 0.; Vertices.push_back(v1);
  v1.x = 4.; v1.y = 9.; v1.z = 0.; Vertices.push_back(v1);
  v1.x = -1.; v1.y = 11.; v1.z = 0.; Vertices.push_back(v1);

  v1.x = -5; v1.y = 5.; v1.z = 0.; Vertices.push_back(v1);
  v1.x = 0.; v1.y = -4.; v1.z = 0.; Vertices.push_back(v1);

  double mag = 0.03;
  for(int i=0;i<Vertices.size();i++){
    Vertices[i].z = -1.0*mag*(Vertices[i].x*Vertices[i].x+
			      Vertices[i].y*Vertices[i].y);
  }


}

void placeVertex1(){

  Vertex v1;
  Triangle tr1;
  v1.x = -19.;  v1.y = 19.; v1.z = -1.; Vertices.push_back(v1);   
  v1.x = +19.;  v1.y = 19.; v1.z = +1.; Vertices.push_back(v1);   
  v1.x =   0.;  v1.y =-19.; v1.z = -2.; Vertices.push_back(v1);  
  tr1.ind1 = 0;   tr1.ind2 = 1;   tr1.ind3 = 2; tr1.removed = 0;
  Triangles.push_back(tr1);

  // v1.x = 0.; v1.y = -0.2; v1.z = 0.; Vertices.push_back(v1);
  // v1.x = 4.; v1.y = 0.1; v1.z = 0.1; Vertices.push_back(v1);
  // v1.x = -5.; v1.y = 3.2; v1.z = -0.1; Vertices.push_back(v1);
  // v1.x = 2.; v1.y = 3.; v1.z = 0.02; Vertices.push_back(v1);
  // v1.x = 4.; v1.y = 3.; v1.z = 0.08; Vertices.push_back(v1);
  // v1.x = 4.; v1.y = 0.; v1.z = -0.1; Vertices.push_back(v1);

  // v1.x = 3.; v1.y = 2.; v1.z = 2.; Vertices.push_back(v1);
  // v1.x = -1.; v1.y = -1.; v1.z = -1.; Vertices.push_back(v1);
  // v1.x = 1.; v1.y = 4.; Vertices.push_back(v1);
  // v1.x = -5.; v1.y = 5.; Vertices.push_back(v1);
  // v1.x = 1.; v1.y = -5.; Vertices.push_back(v1);
  // v1.x = 1.; v1.y = 7.; Vertices.push_back(v1);
  // v1.x = 5.; v1.y = 7.; Vertices.push_back(v1);

  double rad = 8.;
  double x,y,z;
  for(int i=0;i<20;i++){
    double r = rad+1;
    while(r>rad){
      x = 2.*drand48()-1.;
      y = 2.*drand48()-1.;
      z = 2.*drand48()-1.;
      r = sqrt(x*x+y*y+z*z);
    }
    v1.x = x/r*rad;
    v1.y = y/r*rad;
    v1.z = z/r*rad;

    printf("%lf %lf %lf\n",x/r*rad,y/r*rad,z/r*rad);
    Vertices.push_back(v1);
  }
}

class Delaunay2d{

public:
  static void getDelaunay2d();
  static void outputGnuplot();
  static int containsTheVertex(class Triangle _t, class Vertex _newV);
  static bool containsTheVertex2d(class Triangle _t, class Vertex _newV);
  static class Circle getCircumScribedCircle(class Triangle _t);
  static double angleAtVertex(class Triangle _t, int _index);
  static int containsVertexCombination(class Triangle _t, int a, int b, int c);
};

double Delaunay2d::angleAtVertex(class Triangle _t, int _index){

  // measure the angle at the point. 
  int A_index, O_index, B_index;
  if(_t.ind1==_index){
    O_index=_index; A_index=_t.ind2; B_index=_t.ind3;
  }else if(_t.ind2==_index){
    O_index=_index; A_index=_t.ind1; B_index=_t.ind3;
  }else if(_t.ind3==_index){
    O_index=_index; A_index=_t.ind1; B_index=_t.ind2;
  }else{
    printf("index %d is not included. in triangle (%d %d %d)\n",
	   _index, _t.ind1, _t.ind2, _t.ind3); exit(1);
  }

  double ax = Vertices[A_index].x;
  double ay = Vertices[A_index].y;
  double az = Vertices[A_index].z;
  double ox = Vertices[O_index].x;
  double oy = Vertices[O_index].y;
  double oz = Vertices[O_index].z;
  double bx = Vertices[B_index].x;
  double by = Vertices[B_index].y;
  double bz = Vertices[B_index].z;
  
  double oax = ax - ox; 
  double oay = ay - oy; 
  double oaz = az - oz; 
  double oan = sqrt(oax*oax+oay*oay+oaz*oaz);
  double obx = bx - ox; 
  double oby = by - oy; 
  double obz = bz - oz; 
  double obn = sqrt(obx*obx+oby*oby+obz*obz);
  double ipr = oax*obx+oay*oby+oaz*obz;
  double cosaob = ipr/(oan*obn);

   // printf("computing angles=%.2lf: %.1lf %.1lf %.1lf %.1lf %.1lf %.1lf %.1lf %.1lf %.1lf\n",cosaob,ax,ay,az,ox,oy,oz,bx,by,bz);
  
  return cosaob;
}

int Delaunay2d::containsVertexCombination(class Triangle _t, int a, int b, int c){

  int result = 0;

  if((_t.ind1==a && _t.ind2==b && _t.ind3==c) ||
     (_t.ind1==a && _t.ind2==c && _t.ind3==b) ||
     (_t.ind1==b && _t.ind2==a && _t.ind3==c) ||
     (_t.ind1==b && _t.ind2==c && _t.ind3==a) ||
     (_t.ind1==c && _t.ind2==a && _t.ind3==b) ||
     (_t.ind1==c && _t.ind2==b && _t.ind3==a)
     ){
    result = 1;
  }

  return result;
}


class Circle Delaunay2d::getCircumScribedCircle(class Triangle _t){

  class Vertex v1 = Vertices[_t.ind1];
  class Vertex v2 = Vertices[_t.ind2];
  class Vertex v3 = Vertices[_t.ind3];

  double v12midx = 0.5*(v1.x+v2.x);
  double v13midx = 0.5*(v1.x+v3.x);
  double v12midy = 0.5*(v1.y+v2.y);
  double v13midy = 0.5*(v1.y+v3.y);
  double v12midz = 0.5*(v1.z+v2.z);
  // double v13midz = 0.5*(v1.z+v3.z);

  double a12x = v2.x - v1.x;
  double a12y = v2.y - v1.y;
  double a12z = v2.z - v1.z;
  double a13x = v3.x - v1.x;
  double a13y = v3.y - v1.y;
  double a13z = v3.z - v1.z;
  // double a23x = v3.x - v2.x;
  // double a23y = v3.y - v2.y;
  // double a23z = v3.z - v2.z;

  double a12 = sqrt(a12x*a12x+a12y*a12y+a12z*a12z);
  double a13 = sqrt(a13x*a13x+a13y*a13y+a13z*a13z);
  // double a23 = sqrt(a23x*a23x+a23y*a23y+a23z*a23z);

  double a12nx = a12x/a12;
  double a12ny = a12y/a12;
  double a12nz = a12z/a12;

  double a13nx = a13x/a13;
  double a13ny = a13y/a13;
  double a13nz = a13z/a13;

  // double a23nx = a23x/a23;
  // double a23ny = a23y/a23;
  // double a23nz = a23z/a23;

  double normal3on12x = a13x - (a13x*a12nx+a13y*a12ny+a13z*a12nz)*a12nx;
  double normal3on12y = a13y - (a13x*a12nx+a13y*a12ny+a13z*a12nz)*a12ny;
  double normal3on12z = a13z - (a13x*a12nx+a13y*a12ny+a13z*a12nz)*a12nz;
  double norm_normal3on12 = sqrt(normal3on12x*normal3on12x+normal3on12y*normal3on12y+normal3on12z*normal3on12z);
  double norm3on12x = normal3on12x/norm_normal3on12;
  double norm3on12y = normal3on12y/norm_normal3on12;
  double norm3on12z = normal3on12z/norm_normal3on12;

  double normal2on13x = a12x - (a12x*a13nx+a12y*a13ny+a12z*a13nz)*a13nx;
  double normal2on13y = a12y - (a12x*a13nx+a12y*a13ny+a12z*a13nz)*a13ny;
  double normal2on13z = a12z - (a12x*a13nx+a12y*a13ny+a12z*a13nz)*a13nz;
  double norm_normal2on13 = sqrt(normal2on13x*normal2on13x+normal2on13y*normal2on13y+normal2on13z*normal2on13z);
  double norm2on13x = normal2on13x/norm_normal2on13;
  double norm2on13y = normal2on13y/norm_normal2on13;
  // double norm2on13z = normal2on13z/norm_normal2on13;

  double ax = v12midx; 
  double ay = v12midy; 
  // double az = v12midz; 
  double bx = norm3on12x;
  double by = norm3on12y;
  // double bz = norm3on12z;
  double cx = v13midx;
  double cy = v13midy;
  // double cz = v13midz;
  double dx = norm2on13x;
  double dy = norm2on13y;
  // double dz = norm2on13z;
  double ex = ax - cx;
  double ey = ay - cy;
  // double ez = az - cz;
  // double s = (-by*ex+bx*ey)/(-dx*by+dy*bx);
  double t = (-dy*ex+dx*ey)/(-dx*by+dy*bx);

  // printf("s,t=%lf,%lf\n",s,t);
  // printf("set arrow from %lf,%lf,%lf to %lf,%lf,%lf\n",
  // 	 v12midx,v12midy,v12midz,
  // 	 v12midx+norm3on12x*t,v12midy+norm3on12y*t,v12midz+norm3on12z*t);

  // printf("set arrow from %lf,%lf,%lf to %lf,%lf,%lf\n",
  // 	 v13midx,v13midy,v13midz,
  // 	 v13midx+norm2on13x*s,v13midy+norm2on13y*s,v13midz+norm2on13z*s);

  double ccx = v12midx+norm3on12x*t;
  double ccy = v12midy+norm3on12y*t;
  double ccz = v12midz+norm3on12z*t;

  // printf("# norm3on12z:%lf t:%lf s:%ld v123.z = %lf %lf %lf ccz:%lf\n",norm3on12z,t,s,v1.z,v2.z,v3.z,ccz);

  double di1 = sqrt((ccx-v1.x)*(ccx-v1.x)+(ccy-v1.y)*(ccy-v1.y)+(ccz-v1.z)*(ccz-v1.z));
  // double di2 = sqrt((ccx-v2.x)*(ccx-v2.x)+(ccy-v2.y)*(ccy-v2.y)+(ccz-v2.z)*(ccz-v2.z)		    );
  // double di3 = sqrt((ccx-v3.x)*(ccx-v3.x)+(ccy-v3.y)*(ccy-v3.y)+(ccz-v3.z)*(ccz-v3.z));

  class Circle c;
  c.x = ccx;
  c.y = ccy;
  c.z = ccz;
  c.r = di1;
  c.ind1 = _t.ind1;
  c.ind2 = _t.ind2;
  c.ind3 = _t.ind3;
  Circles.push_back(c);

  return c;
}

bool Delaunay2d::containsTheVertex2d(class Triangle _t, Vertex _pt){

  int b1, b2, b3;
  Vertex v1 = Vertices[_t.ind1];
  Vertex v2 = Vertices[_t.ind2];
  Vertex v3 = Vertices[_t.ind3];  

  b1 = sign(_pt, v1, v2) < 0.0f;
  b2 = sign(_pt, v2, v3) < 0.0f;
  b3 = sign(_pt, v3, v1) < 0.0f;

  return ((b1 == b2) && (b2 == b3));

}

int Delaunay2d::containsTheVertex(class Triangle _t, class Vertex _newV){

  // printf("parentT: %d %d %d\n",_t.ind1,_t.ind2,_t.ind3);
  // projection. if contained return 1. otherwise return 0.
  Vertex v1 = Vertices[_t.ind1];
  Vertex v2 = Vertices[_t.ind2];
  Vertex v3 = Vertices[_t.ind3];

  // double intx = (v1.x+v2.x+v3.x)/3.;
  double inty = (v1.y+v2.y+v3.y)/3.;
  double intz = (v1.z+v2.z+v3.z)/3.;

  double a12x = v2.x - v1.x;
  double a12y = v2.y - v1.y;
  double a12z = v2.z - v1.z;
  double a13x = v3.x - v1.x;
  double a13y = v3.y - v1.y;
  double a13z = v3.z - v1.z;
  // double a23x = v3.x - v2.x;
  double a23y = v3.y - v2.y;
  double a23z = v3.z - v2.z;

  double nx = a12y*a13z - a12z*a13y;
  double ny = a12z*a13x - a12x*a13z;
  double nz = a12x*a13y - a12y*a13x;

  double norm = sqrt(nx*nx+ny*ny+nz*nz);
  nx = nx/norm;
  ny = ny/norm;
  nz = nz/norm;

  double projection = (_newV.x-v1.x)*nx + (_newV.y-v1.y)*ny + (_newV.z-v1.z)*nz;
  // double qprojx = _newV.x - projection*nx;
  double qprojy = _newV.y - projection*ny;
  double qprojz = _newV.z - projection*nz;

  // int sign12int = (gaiseki12intx)/(gaiseki12qprojx);
  double dsign12 = (a12y*(intz-v1.z) - a12z*(inty-v1.y))/(a12y*(qprojz-v1.z) - a12z*(qprojy-v1.y));
  double dsign23 = (a23y*(intz-v2.z) - a23z*(inty-v2.y))/(a23y*(qprojz-v2.z) - a23z*(qprojy-v2.y));
  double dsign13 = (a13y*(intz-v1.z) - a13z*(inty-v1.y))/(a13y*(qprojz-v1.z) - a13z*(qprojy-v1.y));
  int sign12 = 1; if(dsign12<0.) sign12 = -1;
  int sign23 = 1; if(dsign23<0.) sign23 = -1;
  int sign13 = 1; if(dsign13<0.) sign13 = -1;

  printf("# signs: %d %d %d\n",sign12,sign13,sign23);

  // printf("set arrow from %.1lf,%.1lf,%.1lf to %.1lf,%.1lf,%.1lf\n",_newV.x,_newV.y,_newV.z,qprojx,qprojy,qprojz);
  // printf("set arrow from %.1lf,%.1lf,%.1lf to %.1lf,%.1lf,%.1lf\n",_newV.x,_newV.y,_newV.z,intx,inty,intz);

  if(sign12==1 && sign23==1 && sign13==1){
    return 1; 
  }else return 0;

}

void Delaunay2d::getDelaunay2d(){

  // add and remove triangles.
  std::vector<class Vertex>::iterator it_v;
  it_v = Vertices.begin()+3;
  while(it_v!=Vertices.end()){

    printf("# unset obj;unset arrow\n");
    Vertex newV = *it_v;
    std::vector<class Triangle>::iterator it_t;
    it_t = Triangles.begin();
    bool divideOnlyOnceforEachVertex = true;

    while(it_t!=Triangles.end() && divideOnlyOnceforEachVertex){

      Triangle parentT = *it_t;
      printf("# newly added vertex: %lf %lf %lf\n",it_v->x,it_v->y,it_v->z);

      if(Delaunay2d::containsTheVertex2d(parentT,newV) && parentT.removed==0){
	
	printf("vertex %d in (%d %d %d)\n",newV.index,
	       parentT.ind1,parentT.ind2,parentT.ind3);
	// it_t->removed = 1;
	class Triangle t1,t2,t3;
	t1.ind1 = newV.index; t1.ind2 = parentT.ind1; t1.ind3 = parentT.ind2;
	t2.ind1 = newV.index; t2.ind2 = parentT.ind1; t2.ind3 = parentT.ind3;
	t3.ind1 = newV.index; t3.ind2 = parentT.ind2; t3.ind3 = parentT.ind3;
	t1.removed = 0; t2.removed = 0; t3.removed = 0;
	it_t->removed=1;

	// Triangles.erase(it_t);
	it_t = Triangles.end();
	it_t = Triangles.insert(it_t,t1);
	it_t = Triangles.insert(it_t,t2);
	it_t = Triangles.insert(it_t,t3);

	divideOnlyOnceforEachVertex = false;
      }
      it_t++;
    }

    // printf("adding now the vertex %d\n",it_v->index);

    int flipped = 0;
    // circum circle. 
    it_t = Triangles.begin(); it_t++;

    while(it_t!=Triangles.end()){

      if(it_t->removed!=0){it_t++;continue;}

      flipped = 0;

      class Circle c;
      class Triangle t; t = *it_t;
      c = Delaunay2d::getCircumScribedCircle(t);

      std::vector<class Vertex>::iterator it_checkv;
      it_checkv = Vertices.begin();

      while((it_checkv->index)<=(it_v->index) && flipped==0){
	printf("checkv:%d it_v:%d\n",it_checkv->index,it_v->index);
	
	Vertex checkv = *it_checkv;
	// check if this vertex is inside each circles. 
	double distance = 
			   sqrt((checkv.x-c.x)*(checkv.x-c.x)+
				(checkv.y-c.y)*(checkv.y-c.y)+
				(checkv.z-c.z)*(checkv.z-c.z));

	printf("Tsize:%ld it_v:%d checkv:%d checking distance (%d %d %d),%d...",
	       Triangles.size(),it_v->index,it_checkv->index,t.ind1,t.ind2,t.ind3,checkv.index);	

	if(distance<c.r*0.99 && flipped==0){

	  printf("in. it_v:%d checkv:%d checking distance (%d %d %d),%d...",
		 it_v->index,it_checkv->index,t.ind1,t.ind2,t.ind3,checkv.index);	
	  std::vector<class Triangle>::iterator it_checkt;
	  it_checkt = Triangles.begin(); it_checkt++;
	  while(it_checkt!=Triangles.end() && flipped==0){

	    if(it_checkt==it_t || it_checkt->removed!=0 || t.removed!=0){
	      it_checkt++;continue;
	    }

	    Triangle flipT = *it_checkt;
	    // printf("iiiit_checkt:(%d %d %d)\n",flipT.ind1,flipT.ind2,flipT.ind3);

	    if(flipped==0 && 
	       containsVertexCombination(flipT, checkv.index, t.ind1, t.ind2)){
	      printf("flip1.");
	      printf("because %d (%d %d %d) is inside cc (%d %d %d), "
		     "sharing (%d %d)\n",
		     checkv.index,flipT.ind1,flipT.ind2,flipT.ind3,
		     t.ind1,t.ind2,t.ind3,t.ind1,t.ind2
		     );
	      printf("flip to (%d %d %d)(%d %d %d)\n",
		     t.ind1,t.ind3,checkv.index,
		     t.ind2,t.ind3,checkv.index);
	      Triangle insertT, insertS;
	      insertT.ind1 = t.ind1;
	      insertT.ind2 = t.ind3;
	      insertT.ind3 = checkv.index;
	      insertS.ind1 = t.ind2;
	      insertS.ind2 = t.ind3;
	      insertS.ind3 = checkv.index;
	      insertT.removed = 0; insertS.removed = 0;
	      it_t->removed = 1;
	      it_checkt->removed = 1;
	      Triangles.push_back(insertT);
	      Triangles.push_back(insertS);
	      flipped++;
	    }

	    if(flipped==0 && 
	       containsVertexCombination(flipT, checkv.index, t.ind1, t.ind3)){
	      printf("flip2.");
	      printf("because %d (%d %d %d) is inside cc (%d %d %d), "
		     "sharing (%d %d)\n",
		     checkv.index,flipT.ind1,flipT.ind2,flipT.ind3,
		     t.ind1,t.ind2,t.ind3,t.ind1,t.ind3
		     );
	      printf("flip to (%d %d %d)(%d %d %d)\n",
		     t.ind1,t.ind2,checkv.index,
		     t.ind2,t.ind3,checkv.index);
	      Triangle insertT, insertS;
	      insertT.ind1 = t.ind1;
	      insertT.ind2 = t.ind2;
	      insertT.ind3 = checkv.index;
	      insertS.ind1 = t.ind2;
	      insertS.ind2 = t.ind3;
	      insertS.ind3 = checkv.index;
	      insertT.removed = 0; insertS.removed = 0;
	      it_t->removed = 1;
	      it_checkt->removed = 1;
	      Triangles.push_back(insertT);
	      Triangles.push_back(insertS);
	      flipped++;
	    }

	    if(flipped==0 && 
	       containsVertexCombination(flipT, checkv.index, t.ind2, t.ind3)){
	      printf("flip3.");
	      printf("because %d (%d %d %d) is inside cc (%d %d %d), "
		     "sharing (%d %d)\n",
		     checkv.index,flipT.ind1,flipT.ind2,flipT.ind3,
		     t.ind1,t.ind2,t.ind3,t.ind2,t.ind3
		     );
	      printf("flip to (%d %d %d)(%d %d %d)\n",
		     t.ind1,t.ind2,checkv.index,
		     t.ind1,t.ind3,checkv.index);
	      Triangle insertT, insertS;
	      insertT.ind1 = t.ind1;
	      insertT.ind2 = t.ind2;
	      insertT.ind3 = checkv.index;
	      insertS.ind1 = t.ind1;
	      insertS.ind2 = t.ind3;
	      insertS.ind3 = checkv.index;
	      insertT.removed = 0; insertS.removed = 0;
	      it_t->removed = 1;
	      it_checkt->removed = 1;
	      Triangles.push_back(insertT);
	      Triangles.push_back(insertS);
	      flipped++;
	    }
	    it_checkt++;
	  }
	}
	it_checkv++;
      }
      if(flipped){
	printf(".");
	it_t=Triangles.begin();
      }
      // if(Triangles.size()>100.) break;

      it_t++;
    }

    it_v++;
  }

}

void Delaunay2d::outputGnuplot(){
  
  FILE *fp;
  if((fp=fopen("out.gp","w"))==NULL){
    printf("cannot open file.\n"); exit(1);
  }

  fprintf(fp,"set view 30,26\n");
  fprintf(fp,"set isosamples 30;\n");  
  fprintf(fp,"unset arrow;unset obj;unset label\n");
  fprintf(fp,"set xrange [-15:15];set yrange [-15:15];set zrange [-8:8]\n");
  fprintf(fp,"set parametric;set urange[0:2*pi];set vrange[-pi/2:pi/2]\n");
  fprintf(fp,"r=0.5;set lt 1 lc rgb \'blue\' lt 1 lw 2;");
  fprintf(fp,"fx(v,u) = r*cos(v)*cos(u);fy(v,u) = r*cos(v)*sin(u);fz(v) = r*sin(v);");

  int show_vertices = 1;
  int show_triangles = 1;

  if(show_triangles){
    std::vector<class Triangle>::iterator it_t;
    it_t = Triangles.begin();
    while(it_t!=Triangles.end()){
      class Vertex v1 = Vertices[it_t->ind1];
      class Vertex v2 = Vertices[it_t->ind2];
      class Vertex v3 = Vertices[it_t->ind3];
      if(it_t->removed==0){
	printf("drawing (%d %d %d)\n",it_t->ind1,it_t->ind2,it_t->ind3);
	fprintf(fp,"# drawing (%d %d %d)\n",it_t->ind1,it_t->ind2,it_t->ind3);
	fprintf(fp,"set arrow from %.1lf, %.1lf, %.1lf to %.1lf, %.1lf, %.1lf nohead;",v1.x,v1.y,v1.z,v3.x,v3.y,v3.z);
	fprintf(fp,"set arrow from %.1lf, %.1lf, %.1lf to %.1lf, %.1lf, %.1lf nohead;",v1.x,v1.y,v1.z,v2.x,v2.y,v2.z);
	fprintf(fp,"set arrow from %.1lf, %.1lf, %.1lf to %.1lf, %.1lf, %.1lf nohead;",v2.x,v2.y,v2.z,v3.x,v3.y,v3.z);
      }
      it_t++;
    }
  }

  if(show_vertices){
    FILE *vfp;
    if((vfp=fopen("vertex.dat","w"))==NULL){
      printf("cannot open file.\n"); exit(1);
    }
    std::vector<class Vertex>::iterator it_v;
    it_v = Vertices.begin();
    fprintf(vfp,"# x y z index\n");
    while(it_v!=Vertices.end()){
      fprintf(vfp,"%.1lf %.1lf %.1lf %d\n",it_v->x,it_v->y,it_v->z,it_v->index);
      fprintf(fp,"set label \"%d\" at %lf,%lf,%lf\n",it_v->index,it_v->x+0.2,it_v->y+0.5,it_v->z+0.5);
      it_v++;
    }
    fclose(vfp);
  }

  std::vector<class Circle>::iterator it_c;
  it_c = Circles.begin();
  fprintf(fp,"set parametric;set urange [0:2*pi];set vrange [-pi:pi];");
  fprintf(fp,"fx(v,u) = cos(v)*cos(u);fy(v,u) = cos(v)*sin(u);fz(v)   = sin(v)\n");
  fprintf(fp,"splot \"vertex.dat\" pt 7 ps 2 lc 9, ");
  while(it_c!=Circles.end()){
    if(it_c->ind1==1 || it_c->ind2==1 || it_c->ind3==1){
      // fprintf(fp,"%lf*fx(v,u)+%lf,%lf*fy(v,u)+%lf,%lf*fz(v)+%lf ls 2 lw 1, ",it_c->r,it_c->x,it_c->r,it_c->y,it_c->r,it_c->z);
    }
    it_c++;
  }
  // fprintf(fp,"%lf*fx(v,u)+%lf,%lf*fy(v,u)+%lf,%lf*fz(v)+%lf ls 2 lw 1",it_c->r,it_c->x,it_c->r,it_c->y,it_c->r,it_c->z);

  fclose(fp);

}

int main(){

  srand48((unsigned)time(NULL));

  placeVertex2();

  assignIndex();

  Delaunay2d::getDelaunay2d();

  Delaunay2d::outputGnuplot();

  return 0;
}
