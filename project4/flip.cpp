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
  Vertex p1, p2, p3;  // vertices
  int ind1, ind2, ind3;      // vertex indices
  // double cx, cy, cz;  // circum circle
  // double radius;      // radius of circum circle
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

void drawt(){

  std::vector<class Triangle>::iterator it_t;
  it_t = Triangles.begin();
  while(it_t!=Triangles.end()){
    Triangle t = *it_t;  // 三角形取得
    printf("set arrow from %.1lf, %.1lf to %.1lf, %.1lf nohead;",t.p1.x,t.p1.y,t.p2.x,t.p2.y);
    printf("set arrow from %.1lf, %.1lf to %.1lf, %.1lf nohead;",t.p1.x,t.p1.y,t.p3.x,t.p3.y);
    printf("set arrow from %.1lf, %.1lf to %.1lf, %.1lf nohead;",t.p2.x,t.p2.y,t.p3.x,t.p3.y);
    printf("\n");
    it_t++;
  }
  // printf("set xrange [-30:30];set yrange [-30:30]\n");
  printf("set xrange [-10:10];set yrange [-10:10]\n");
  printf("plot x-1e10\n");
}


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
  // printf("index end.\n");
}

void placeVertex1(){

  Vertex v1;
  v1.x = 2.; v1.y = 3.; Vertices.push_back(v1);
  v1.x = -1.; v1.y = -1.; Vertices.push_back(v1);
  v1.x = 1.; v1.y = 4.; Vertices.push_back(v1);
  v1.x = -5.; v1.y = 5.; Vertices.push_back(v1);
  v1.x = 1.; v1.y = -5.; Vertices.push_back(v1);
  v1.x = 1.; v1.y = 7.; Vertices.push_back(v1);
  v1.x = 5.; v1.y = 7.; Vertices.push_back(v1);

}

void placeVertex2(){

  Vertex v1;

  double theta = 0.;

  v1.x = 0.; v1.y = 0.; Vertices.push_back(v1);
  double r = 2.;
  for(theta=0.;theta<2.;theta+=0.3){
    v1.x = r*cos(theta); v1.y = r*sin(theta); Vertices.push_back(v1);
  }

}

void placeVertex3(){

  Vertex v1;

  // int i;
  // for(i=0;i<10;i++){
  //   v1.x = (double)(i%10);
  //   v1.y = i - (double)(i%10);
  //   v1.x *= 0.5;
  //   v1.y *= 0.4;
  //   Vertices.push_back(v1);
  //   printf("# v1:%lf %lf\n",v1.x,v1.y);
  // }
  v1.x = -3.; v1.y = 0.; Vertices.push_back(v1);
  v1.x = -2.; v1.y = 0.; Vertices.push_back(v1);
  v1.x = -1.; v1.y = 0.5; Vertices.push_back(v1);
  v1.x = 0.; v1.y = 1.0; Vertices.push_back(v1);
  v1.x = 1.; v1.y = -1.0; Vertices.push_back(v1);
  v1.x = 2.; v1.y = 0.0; Vertices.push_back(v1);
  v1.x = 3.; v1.y = -1.0; Vertices.push_back(v1);
  v1.x = 0.; v1.y = 3.0; Vertices.push_back(v1);

  v1.x = -5.; v1.y = 5.; Vertices.push_back(v1);
  // exit(1);
}

void placeVertex4(){

  Vertex v1;

  int i;
  for(i=0;i<5;i++){
    // v1.x = (double)(i%10);
    // v1.y = i - (double)(i%10);
    v1.x = 5.*drand48();
    v1.y = 5.*drand48();
    // v1.x *= 0.5;
    // v1.y *= 0.4;
    Vertices.push_back(v1);
    // printf("# v1:%lf %lf\n",v1.x,v1.y);
  }
}

int main(){

  srand48((unsigned)time(NULL));

  Vertex v1;
  Triangle tr1;

  v1.x = -9.;  v1.y = 9.; Vertices.push_back(v1);   tr1.p1.x = v1.x;  tr1.p1.y = v1.y; 
  v1.x = +9.;  v1.y = 9.; Vertices.push_back(v1);   tr1.p2.x = v1.x;  tr1.p2.y = v1.y; 
  v1.x =   0.;  v1.y =-9.; Vertices.push_back(v1);  tr1.p3.x = v1.x;  tr1.p3.y = v1.y; 
  tr1.ind1 = 0;   tr1.ind2 = 1;   tr1.ind3 = 2;
  Triangles.push_back(tr1);

  // printf("tr1 index:%d %d %d\n",tr1.ind1,tr1.ind2,tr1.ind3);
  // for(int i=0;i<8;i++){
  //   // double s = drand48(); double t = drand48();
  //   v1.x = 10.*(2.*drand48()-1.);  v1.y = 10.*(2.*drand48()-1.); Vertices.push_back(v1);
  // }
  
  placeVertex4();

  assignIndex();

  printf("unset arrow;unset obj;unset label\n");
  std::vector<class Vertex>::iterator it_v;
  it_v = Vertices.begin();
  while(it_v!=Vertices.end()){
    printf("set obj circle at %.1lf, %.1lf size %.1lf fc rgb \"blue\" fs transparent solid 0.5\n",
  	   it_v->x,it_v->y,0.2);
    printf("set label \"%d\" at %lf,%lf\n",it_v->index,it_v->x,it_v->y+0.5);
    // printf("vertex: %.1lf %.1lf\n",it->x,it->y);
    it_v++;
  }

  // add and remove triangles.
  Triangle t1,t2,t3;
  std::vector<class Vertex>::iterator it_v2;
  it_v2 = Vertices.begin()+3;
  while(it_v2!=Vertices.end()){
    
    printf("unset obj;unset arrow\n");
    Vertex newV = *it_v2;
    std::vector<class Triangle>::iterator it_t2;
    it_t2 = Triangles.begin();
    bool eachvertex = true;
    while(it_t2!=Triangles.end() && eachvertex){

      Triangle parentT = *it_t2;
      // printf("parentT index:%d %d %d\n",parentT.ind1,parentT.ind2,parentT.ind3);

      if(containsTheVertex(parentT,newV)){

	printf("# contains.\n");

	t1.p1.x = newV.x; t1.p2.x = parentT.p1.x; t1.p3.x = parentT.p2.x;
	t1.p1.y = newV.y; t1.p2.y = parentT.p1.y; t1.p3.y = parentT.p2.y;
	t1.ind1 = newV.index; t1.ind2 = parentT.ind1; t1.ind3 = parentT.ind2;

	t2.p1.x = newV.x; t2.p2.x = parentT.p1.x; t2.p3.x = parentT.p3.x;
	t2.p1.y = newV.y; t2.p2.y = parentT.p1.y; t2.p3.y = parentT.p3.y;
	t2.ind1 = newV.index; t2.ind2 = parentT.ind1; t2.ind3 = parentT.ind3;

	t3.p1.x = newV.x; t3.p2.x = parentT.p2.x; t3.p3.x = parentT.p3.x;
	t3.p1.y = newV.y; t3.p2.y = parentT.p2.y; t3.p3.y = parentT.p3.y;
	t3.ind1 = newV.index; t3.ind2 = parentT.ind2; t3.ind3 = parentT.ind3;

	Triangles.erase(it_t2);
	printf("# object index: %d %d %d, %d %d %d, %d %d %d\n",
	       t1.ind1,t1.ind2,t1.ind3,
	       t2.ind1,t2.ind2,t2.ind3,
	       t3.ind1,t3.ind2,t3.ind3
	       );
	printf("set obj circle at %.1lf, %.1lf size %.1lf fc rgb \"red\" fs transparent solid 0.5\n",
	       newV.x,newV.y,0.2);

	it_t2 = Triangles.end();
	it_t2 = Triangles.insert(it_t2,t1);
	it_t2 = Triangles.insert(it_t2,t2);
	it_t2 = Triangles.insert(it_t2,t3);
	eachvertex = false;
      }
      it_t2++;
    }

    // circum circle. 
    std::vector<class Triangle>::iterator it_t3;
    it_t3 = Triangles.begin();
    while(it_t3!=Triangles.end()){

      Triangle t = *it_t3;
      double x1 = t.p1.x;  double y1 = t.p1.y;
      double x2 = t.p2.x;  double y2 = t.p2.y;
      double x3 = t.p3.x;  double y3 = t.p3.y;

      double m = 2.0*((x2-x1)*(y3-y1)-(y2-y1)*(x3-x1));  
      double x = ((y3-y1)*(x2*x2-x1*x1+y2*y2-y1*y1)
    		  +(y1-y2)*(x3*x3-x1*x1+y3*y3-y1*y1))/m;
      double y = ((x1-x3)*(x2*x2-x1*x1+y2*y2-y1*y1)
    		  +(x2-x1)*(x3*x3-x1*x1+y3*y3-y1*y1))/m;

      double cx = x; 		// circum circle coordinates
      double cy = y;
      
      double dx   = t.p1.x - x;
      double dy   = t.p1.y - y;
      double radius2 = dx * dx + dy * dy;

      // printf("m:%.1lf x,y:%.1lf %.1lf\n",m,x,y);
      printf("set obj circle at %.1lf, %.1lf size %.1lf fc rgb \"blue\"\n",
      	     cx,cy,sqrt(radius2));
      printf("set label \"(%d %d %d)\" at %lf,%lf\n",it_t3->ind1,it_t3->ind2,it_t3->ind3,cx,cy);

      std::vector<class Vertex>::iterator it_v3;
      it_v3 = Vertices.begin()+3;
      int number = 0;
      while(it_v3!=it_v2+1){

      	Vertex checkv = *it_v3;
      	double dist2 = (checkv.x-cx)*(checkv.x-cx) + (checkv.y-cy)*(checkv.y-cy);
      	// printf("dist:%.1lf radius:%.1lf",sqrt(dist2),sqrt(radius2));
      	// printf("# point %.1lf %.1lf inside.",checkv.x,checkv.y);
	// printf("# dist2:(%.1lf %.1lf) %.1lf"
	//        " radius2:(the circle (%.1lf,%.1lf),(%.1lf,%.1lf),(%.1lf,%.1lf)) %.1lf "
	//        "# %.1lf %.1lf.\n",
	//        checkv.x,checkv.y,dist2,t.p1.x,t.p1.y,t.p2.x,t.p2.y,t.p3.x,t.p3.y,
	//        radius2,checkv.x,checkv.y);

      	if(dist2<radius2*0.99){
	  printf("# radius2: (%d %d %d) > dist2 vertex %d\n",t.ind1,t.ind2,t.ind3,checkv.index);
	  std::vector<class Triangle>::iterator it_t4;
	  it_t4 = Triangles.begin();
	  while(it_t4<Triangles.end()){
	    // find triangle checkv + (two of vertices of it_t3);
	    // then flip.
	    Triangle flipT = *it_t4;
	    
	    if(
	       (flipT.ind1==checkv.index && flipT.ind2==it_t3->ind1 && flipT.ind3==it_t3->ind2) ||
	       (flipT.ind1==checkv.index && flipT.ind2==it_t3->ind2 && flipT.ind3==it_t3->ind1) ||
	       (flipT.ind2==checkv.index && flipT.ind1==it_t3->ind1 && flipT.ind3==it_t3->ind2) ||
	       (flipT.ind2==checkv.index && flipT.ind1==it_t3->ind2 && flipT.ind3==it_t3->ind1) ||
	       (flipT.ind3==checkv.index && flipT.ind1==it_t3->ind1 && flipT.ind3==it_t3->ind2) ||
	       (flipT.ind3==checkv.index && flipT.ind1==it_t3->ind2 && flipT.ind3==it_t3->ind1)
	       ){
	      printf("# overlap: (%d %d %d) & %d\n",flipT.ind1,flipT.ind2,flipT.ind3,checkv.index);
	      printf("# 3need to erase (%d %d %d),(%d %d %d)\n",
		     it_t3->ind1,it_t3->ind2,it_t3->ind3,it_t4->ind1,it_t4->ind2,it_t4->ind3);
	      printf("# 3need to insert (%d %d %d), (%d %d %d)\n",
		     it_t3->ind1,it_t3->ind3,checkv.index,it_t3->ind2,it_t3->ind3,checkv.index);

	      Triangle insertT,insertS;
	      insertT.ind1 = it_t3->ind1; 
	      insertT.ind2 = it_t3->ind3; 
	      insertT.ind3 = checkv.index;
	      insertT.p1.x = it_t3->p1.x; insertT.p1.y = it_t3->p1.y; insertT.p1.z = it_t3->p1.z;
	      insertT.p2.x = it_t3->p3.x; insertT.p2.y = it_t3->p3.y; insertT.p2.z = it_t3->p3.z;
	      insertT.p3.x = checkv.x; insertT.p3.y = checkv.y; insertT.p3.z = checkv.z;

	      insertS.ind1 = it_t3->ind2; 
	      insertS.ind2 = it_t3->ind3; 
	      insertS.ind3 = checkv.index;
	      insertS.p1.x = it_t3->p2.x; insertS.p1.y = it_t3->p2.y; insertS.p1.z = it_t3->p2.z;
	      insertS.p2.x = it_t3->p3.x; insertS.p2.y = it_t3->p3.y; insertS.p2.z = it_t3->p3.z;
	      insertS.p3.x = checkv.x; insertS.p3.y = checkv.y; insertS.p3.z = checkv.z;

	      Triangles.erase(it_t4);
	      Triangles.erase(it_t3);

	      it_t4 = Triangles.end();
	      it_t4 = Triangles.insert(it_t4,insertT);
	      it_t4 = Triangles.insert(it_t4,insertS);

	    }else if(
		     (flipT.ind1==checkv.index && flipT.ind2==it_t3->ind1 && flipT.ind3==it_t3->ind3) ||
		     (flipT.ind1==checkv.index && flipT.ind2==it_t3->ind3 && flipT.ind3==it_t3->ind1) ||
		     (flipT.ind2==checkv.index && flipT.ind1==it_t3->ind1 && flipT.ind3==it_t3->ind3) ||
		     (flipT.ind2==checkv.index && flipT.ind1==it_t3->ind3 && flipT.ind3==it_t3->ind1) ||
		     (flipT.ind3==checkv.index && flipT.ind1==it_t3->ind1 && flipT.ind3==it_t3->ind3) ||
		     (flipT.ind3==checkv.index && flipT.ind1==it_t3->ind3 && flipT.ind3==it_t3->ind1)
		     ){
	      printf("# overlap: (%d %d %d) & %d\n",flipT.ind1,flipT.ind2,flipT.ind3,checkv.index);
	      printf("# 3need to erase (%d %d %d),(%d %d %d)\n",
		     it_t3->ind1,it_t3->ind2,it_t3->ind3,it_t4->ind1,it_t4->ind2,it_t4->ind3);
	      printf("# 3need to insert (%d %d %d), (%d %d %d)\n",
		     it_t3->ind1,it_t3->ind2,checkv.index,it_t3->ind2,it_t3->ind3,checkv.index);


	      Triangle insertT,insertS;
	      insertT.ind1 = it_t3->ind1; 
	      insertT.ind2 = it_t3->ind2; 
	      insertT.ind3 = checkv.index;
	      insertT.p1.x = it_t3->p1.x; insertT.p1.y = it_t3->p1.y; insertT.p1.z = it_t3->p1.z;
	      insertT.p2.x = it_t3->p2.x; insertT.p2.y = it_t3->p2.y; insertT.p2.z = it_t3->p2.z;
	      insertT.p3.x = checkv.x; insertT.p3.y = checkv.y; insertT.p3.z = checkv.z;

	      insertS.ind1 = it_t3->ind2; 
	      insertS.ind2 = it_t3->ind3; 
	      insertS.ind3 = checkv.index;
	      insertS.p1.x = it_t3->p2.x; insertS.p1.y = it_t3->p2.y; insertS.p1.z = it_t3->p2.z;
	      insertS.p2.x = it_t3->p3.x; insertS.p2.y = it_t3->p3.y; insertS.p2.z = it_t3->p3.z;
	      insertS.p3.x = checkv.x; insertS.p3.y = checkv.y; insertS.p3.z = checkv.z;

	      Triangles.erase(it_t4);
	      Triangles.erase(it_t3);

	      it_t4 = Triangles.end();
	      it_t4 = Triangles.insert(it_t4,insertT);
	      it_t4 = Triangles.insert(it_t4,insertS);

	    }else if(
		     (flipT.ind1==checkv.index && flipT.ind2==it_t3->ind2 && flipT.ind3==it_t3->ind3) ||
		     (flipT.ind1==checkv.index && flipT.ind2==it_t3->ind3 && flipT.ind3==it_t3->ind2) ||
		     (flipT.ind2==checkv.index && flipT.ind1==it_t3->ind2 && flipT.ind3==it_t3->ind3) ||
		     (flipT.ind2==checkv.index && flipT.ind1==it_t3->ind3 && flipT.ind3==it_t3->ind2) ||
		     (flipT.ind3==checkv.index && flipT.ind1==it_t3->ind2 && flipT.ind3==it_t3->ind3) ||
		     (flipT.ind3==checkv.index && flipT.ind1==it_t3->ind3 && flipT.ind3==it_t3->ind2)
		     ){
	      printf("# flipT: (%d %d %d), it_3:(%d %d %d), checkv:%d\n",flipT.ind1,flipT.ind2,flipT.ind3,
		     it_t3->ind1,it_t3->ind2,it_t3->ind3,checkv.index);
	      printf("# overlap: (%d %d %d) & %d\n",flipT.ind1,flipT.ind2,flipT.ind3,checkv.index);
	      printf("# 3need to erase (%d %d %d),(%d %d %d)\n",
		     it_t3->ind1,it_t3->ind2,it_t3->ind3,it_t4->ind1,it_t4->ind2,it_t4->ind3);
	      printf("# 3need to insert (%d %d %d), (%d %d %d)\n",
		     it_t3->ind1,it_t3->ind2,checkv.index,it_t3->ind1,it_t3->ind3,checkv.index);

	      Triangle insertT,insertS;
	      insertT.ind1 = it_t3->ind1; 
	      insertT.ind2 = it_t3->ind2; 
	      insertT.ind3 = checkv.index;
	      insertT.p1.x = it_t3->p1.x; insertT.p1.y = it_t3->p1.y; insertT.p1.z = it_t3->p1.z;
	      insertT.p2.x = it_t3->p2.x; insertT.p2.y = it_t3->p2.y; insertT.p2.z = it_t3->p2.z;
	      insertT.p3.x = checkv.x; insertT.p3.y = checkv.y; insertT.p3.z = checkv.z;

	      insertS.ind1 = it_t3->ind1; 
	      insertS.ind2 = it_t3->ind3; 
	      insertS.ind3 = checkv.index;
	      insertS.p1.x = it_t3->p1.x; insertS.p1.y = it_t3->p1.y; insertS.p1.z = it_t3->p1.z;
	      insertS.p2.x = it_t3->p3.x; insertS.p2.y = it_t3->p3.y; insertS.p2.z = it_t3->p3.z;
	      insertS.p3.x = checkv.x; insertS.p3.y = checkv.y; insertS.p3.z = checkv.z;

	      Triangles.erase(it_t4);
	      Triangles.erase(it_t3);

	      it_t4 = Triangles.end();
	      it_t4 = Triangles.insert(it_t4,insertT);
	      it_t4 = Triangles.insert(it_t4,insertS);
	    }
	    it_t4++;
	  }

      	  // printf("dist2:%.1lf radius2:%.1lf # %.1lf %.1lf inside.\n",
	  // 	 dist2,radius2,checkv.x,checkv.y);
      	  // remove it_t3 from Triangles.
      	  // add flipped. 
      	}
      	number++;
      	it_v3++;
      }
      printf("# %d particles has been added.\n",number);
      showNumberofTriangles();
      it_t3 ++;
    }
    it_v2 ++;
  }
  drawt();

  // printf("vertex size is: %ld\n",Vertices.size());

  return 0;
}
