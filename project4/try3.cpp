#include <cfloat>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <vector>

#include <map>

class Circle;
class Vertex{
public:
  double x;
  double y;
  double z;
  const void printCoordinate(){
    printf("v: %.2lf %.2lf %.2lf\n",x,y,z);
  }
};

class Circle 
{
public:
  Vertex center;  // center coordinate
  double radius;  // radius
};

class Triangle {
public:
  Vertex p1, p2, p3;  // 頂点座標
  const void printCoordinate(){
    printf("# printing coordinates:\n");
    p1.printCoordinate();
    p2.printCoordinate();
    p3.printCoordinate();
  }
};

std::vector<Vertex> Vertices;
std::vector<Triangle> Triangles;

float sign(Vertex p1, Vertex p2, Vertex p3){
  return (p1.x - p3.x) * (p2.y - p3.y) - (p2.x - p3.x) * (p1.y - p3.y);
}

bool containsTheVertex(Triangle _t, Vertex _pt){
  int b1, b2, b3;
  Vertex v1 = _t.p1;
  Vertex v2 = _t.p2;
  Vertex v3 = _t.p3;  

  b1 = sign(_pt, v1, v2) < 0.0f;
  b2 = sign(_pt, v2, v3) < 0.0f;
  b3 = sign(_pt, v3, v1) < 0.0f;

  return ((b1 == b2) && (b2 == b3));
}


int main(){

  srand48((unsigned)time(NULL));

  Vertex v1;
  Triangle tr1;
  for(int i=0;i<2;i++){
    double s = drand48();     double t = drand48();
    v1.x = -4. + s*8.; v1.y = 4.-8.*t; Vertices.push_back(v1);
    // v1.x = -5.*10.*drand48() - ;  v1.y = 5.*(drand48()-1.); Vertices.push_back(v1);
  }

  v1.x = -9.;  v1.y = 9.; Vertices.push_back(v1);   tr1.p1.x = v1.x;  tr1.p1.y = v1.y; 
  v1.x = +9.;  v1.y = 9.; Vertices.push_back(v1);   tr1.p2.x = v1.x;  tr1.p2.y = v1.y; 
  v1.x =   0.;  v1.y =-9.; Vertices.push_back(v1);  tr1.p3.x = v1.x;  tr1.p3.y = v1.y; 
  Triangles.push_back(tr1);

  printf("unset arrow;unset obj;");
  std::vector<class Vertex>::iterator it_v;
  it_v = Vertices.begin();
  while(it_v!=Vertices.end()){
    printf("set obj circle at %.2lf, %.2lf size %.2lf fc rgb \"blue\" fs transparent solid 0.5\n",
	   it_v->x,it_v->y,0.2);
    // printf("vertex: %.2lf %.2lf\n",it->x,it->y);
    it_v++;
  }

  // add and remove triangles.
  Triangle t1,t2,t3;
  std::vector<class Vertex>::iterator it_v2;
  it_v2 = Vertices.begin();

  Vertex newV = *it_v2;
  std::vector<class Triangle>::iterator it_t2;
  it_t2 = Triangles.begin();

  Triangle parentT = *it_t2;

  Triangles.erase(it_t2);

  t1.p1.x = newV.x; t1.p2.x = parentT.p1.x; t1.p3.x = parentT.p2.x;
  t1.p1.y = newV.y; t1.p2.y = parentT.p1.y; t1.p3.y = parentT.p2.y;

  t2.p1.x = newV.x; t2.p2.x = parentT.p1.x; t2.p3.x = parentT.p3.x;
  t2.p1.y = newV.y; t2.p2.y = parentT.p1.y; t2.p3.y = parentT.p3.y;

  t3.p1.x = newV.x; t3.p2.x = parentT.p2.x; t3.p3.x = parentT.p3.x;
  t3.p1.y = newV.y; t3.p2.y = parentT.p2.y; t3.p3.y = parentT.p3.y;

  it_t2 = Triangles.insert(it_t2,t1);
  it_t2 = Triangles.insert(it_t2,t2);
  it_t2 = Triangles.insert(it_t2,t3);


  // second vertex.
  it_v2++;

  it_t2 = Triangles.begin()+2;
  parentT = *it_t2;
  newV = *it_v2;
  Triangles.erase(it_t2);

  t1.p1.x = newV.x; t1.p2.x = parentT.p1.x; t1.p3.x = parentT.p2.x;
  t1.p1.y = newV.y; t1.p2.y = parentT.p1.y; t1.p3.y = parentT.p2.y;

  t2.p1.x = newV.x; t2.p2.x = parentT.p1.x; t2.p3.x = parentT.p3.x;
  t2.p1.y = newV.y; t2.p2.y = parentT.p1.y; t2.p3.y = parentT.p3.y;

  t3.p1.x = newV.x; t3.p2.x = parentT.p2.x; t3.p3.x = parentT.p3.x;
  t3.p1.y = newV.y; t3.p2.y = parentT.p2.y; t3.p3.y = parentT.p3.y;

  it_t2 = Triangles.insert(it_t2,t1);
  it_t2 = Triangles.insert(it_t2,t2);
  it_t2 = Triangles.insert(it_t2,t3);

  std::vector<class Triangle>::iterator it_t;
  it_t = Triangles.begin();
  while(it_t!=Triangles.end()){
    Triangle t = *it_t;  // 三角形取得
    printf("set arrow from %.2lf, %.2lf to %.2lf, %.2lf nohead;",t.p1.x,t.p1.y,t.p2.x,t.p2.y);
    printf("set arrow from %.2lf, %.2lf to %.2lf, %.2lf nohead;",t.p1.x,t.p1.y,t.p3.x,t.p3.y);
    printf("set arrow from %.2lf, %.2lf to %.2lf, %.2lf nohead;",t.p2.x,t.p2.y,t.p3.x,t.p3.y);
    printf("\n");
    it_t++;
  }
  printf("set xrange [-10:10];set yrange [-10:10]\n");
  printf("plot x-1e10\n");



  return 0;
}
