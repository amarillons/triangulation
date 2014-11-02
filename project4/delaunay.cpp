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
    // printf("v: %.2lf %.2lf %.2lf\n",x,y,z);
  }
  // equivalence 
  bool operator==(Vertex& _v){
    return (fabs(x-_v.x)<fabs(x)*0.01 && fabs(y-_v.y)<fabs(y)*0.01);
  }
  bool operator<(Vertex& _v){
    return x != _v.x ? x < _v.x : y < _v.y;
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
  // http://stackoverflow.com/questions/2049582/how-to-determine-a-point-in-a-triangle
  float sign(class Vertex p1, class Vertex p2, class Vertex p3){
    return (p1.x - p3.x) * (p2.y - p3.y) - (p2.x - p3.x) * (p1.y - p3.y);
  }
  // decide whether this is identical with triangle _t.
  bool operator==(Triangle& _t){
    return(((p1 == _t.p1) && (p2 == _t.p2) && (p3 == _t.p3)) ||
	   ((p1 == _t.p2) && (p2 == _t.p3) && (p3 == _t.p1)) ||
	   ((p1 == _t.p3) && (p2 == _t.p1) && (p3 == _t.p2)) ||
	   ((p1 == _t.p3) && (p2 == _t.p2) && (p3 == _t.p1)) ||
	   ((p1 == _t.p2) && (p2 == _t.p1) && (p3 == _t.p3)) ||
	   ((p1 == _t.p1) && (p2 == _t.p3) && (p3 == _t.p2)));

  }
  int containsTheVertex(class Vertex pt);
  Triangle();
};

Triangle::Triangle(){


}

int Triangle::containsTheVertex(class Vertex pt){
  int b1, b2, b3;
  class Vertex v1 = this->p1;
  class Vertex v2 = this->p2;
  class Vertex v3 = this->p3;  

  b1 = sign(pt, v1, v2) < 0.0f;
  b2 = sign(pt, v2, v3) < 0.0f;
  b3 = sign(pt, v3, v1) < 0.0f;

  return ((b1 == b2) && (b2 == b3));
}

std::vector<Vertex> Vertices;
std::vector<Triangle> Triangles;

class Delaunay2d{

public:
  static void sayHello(){
    printf("this is sayhello.\n");
  }
  // typedef const std::vector<Vertex> VerticesType;
  // typedef const std::vector<Triangle> TrianglesType;

  static void testForDelaunay(){

    // std::vector<class Vertex>::iterator addVerticesOnebyOneIterator;
    // addVerticesOnebyOneIterator = Vertices.begin();

    // Vertex newVertex = *addVerticesOnebyOneIterator;
    Vertex newVertex = Vertices[0];

    std::vector<class Triangle>::iterator findParentTriangleAndDivideIterator;
    findParentTriangleAndDivideIterator = Triangles.begin();
    while(findParentTriangleAndDivideIterator!=Triangles.end()){
      // search for the triangle that newVertex is in and divide it. store it. remove tIter.

      if(findParentTriangleAndDivideIterator->containsTheVertex(newVertex)){
	printf("# This triangle contains the vertex.\n");
      }
      Triangle parentTriangle = *findParentTriangleAndDivideIterator;
      // divide it. 
      // parentTriangle.printCoordinate();

      Triangle t1,t2,t3;

      // koko ni mondai ga
      // t1.p1.x = 0.1; 	t1.p1.y = 0.1; 	t1.p1.z = 0.1;
      // t1.p2.x = 0.1; 	t1.p2.y = 0.1; 	t1.p2.z = 0.1;
      // t1.p3.x = 0.1; 	t1.p3.y = 0.1; 	t1.p3.z = 0.1;

      t1.p1 = newVertex; t1.p2 = parentTriangle.p1; t1.p3 = parentTriangle.p2;
      t2.p1 = newVertex; t2.p2 = parentTriangle.p2; t2.p3 = parentTriangle.p3;
      t3.p1 = newVertex; t3.p2 = parentTriangle.p1; t3.p3 = parentTriangle.p3;

      Triangles.pop_back();

      Triangles.push_back(t1);
      Triangles.push_back(t2);
      Triangles.push_back(t3);

      // Triangles.erase(findParentTriangleAndDivideIterator+1);
      // findParentTriangleAndDivideIterator++;
      return;
    }
  }

  static void getDelaunayTriangles(){
  
    std::vector<class Vertex>::iterator addVerticesOnebyOneIterator;
    addVerticesOnebyOneIterator = Vertices.begin();
    while(addVerticesOnebyOneIterator!=Vertices.end()){

      Vertex newVertex = *addVerticesOnebyOneIterator;
      Triangle t1,t2,t3;

      std::vector<class Triangle>::iterator findParentTriangleAndDivideIterator;
      findParentTriangleAndDivideIterator = Triangles.begin();
      while(findParentTriangleAndDivideIterator!=Triangles.end()){

      	// search for the triangle that newVertex is in and divide it. store it. remove tIter.
      	if(findParentTriangleAndDivideIterator->containsTheVertex(newVertex)){
      	  printf("# This triangle contains the vertex.\n");

      	  Triangle parentTriangle = *findParentTriangleAndDivideIterator;
      	  // divide it. 
      	  parentTriangle.printCoordinate(); 

      	  t1.p1 = newVertex; t1.p2 = parentTriangle.p1; t1.p3 = parentTriangle.p2;
      	  t2.p1 = newVertex; t2.p2 = parentTriangle.p2; t2.p3 = parentTriangle.p3;
      	  t3.p1 = newVertex; t3.p2 = parentTriangle.p1; t3.p3 = parentTriangle.p3;

      	  Triangles.erase(findParentTriangleAndDivideIterator);

	  findParentTriangleAndDivideIterator = Triangles.insert(t1);
	  findParentTriangleAndDivideIterator = Triangles.insert(t2);
	  findParentTriangleAndDivideIterator = Triangles.insert(t3);
	}
      	findParentTriangleAndDivideIterator++;
      }
      addVerticesOnebyOneIterator++;
    }
  }

};

int main(){

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

  // Delaunay2d::testForDelaunay();
  Delaunay2d::getDelaunayTriangles();

  printf("# Triangle size:%ld\n",Triangles.size());

  // exit(0);
  std::vector<class Triangle>::iterator findParentTriangleAndDivideIterator;
  findParentTriangleAndDivideIterator = Triangles.begin();
  while(findParentTriangleAndDivideIterator!=Triangles.end()){

    findParentTriangleAndDivideIterator++;
  }

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
