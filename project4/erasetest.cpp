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
};

class myClass{
 public:
  Vertex v1;
  Vertex v2;
  Vertex v3;
};

int main(){


  std::vector<myClass>myClasses;
  myClass mine1;

  Vertex myVertex;
  myVertex.x = 1.0; myVertex.y = 1.1; myVertex.z = 1.2;
  mine1.v1 = myVertex;
  myVertex.x = 2.0; myVertex.y = 2.1; myVertex.z = 2.2;
  mine1.v1 = myVertex;
  myVertex.x = 3.0; myVertex.y = 3.1; myVertex.z = 3.2;
  mine1.v1 = myVertex;

  myClasses.push_back(mine1);

  myVertex.x = 11.0; myVertex.y = 11.1; myVertex.z = 11.2;
  mine1.v1 = myVertex;
  myVertex.x = 12.0; myVertex.y = 12.1; myVertex.z = 12.2;
  mine1.v1 = myVertex;
  myVertex.x = 13.0; myVertex.y = 13.1; myVertex.z = 13.2;
  mine1.v1 = myVertex;

  myClasses.push_back(mine1);

  myVertex.x = 21.0; myVertex.y = 21.1; myVertex.z = 21.2;
  mine1.v1 = myVertex;
  myVertex.x = 22.0; myVertex.y = 22.1; myVertex.z = 22.2;
  mine1.v1 = myVertex;
  myVertex.x = 23.0; myVertex.y = 23.1; myVertex.z = 23.2;
  mine1.v1 = myVertex;

  myClasses.push_back(mine1);

  std::vector<class myClass>::iterator it;
  it = myClasses.begin();
  myClasses.erase(it);

  std::vector<class myClass>::iterator it2;
  it2 = myClasses.begin();
  while(it2!=myClasses.end()){

    myClass disp = *it2;

    printf("%lf %lf %lf\n",disp.v1.x,disp.v1.y,disp.v1.z);
    it2++;

  }

  // myClass mine1;
  // mine1.x = 1.0; mine1.y = 1.1; mine1.z = 1.2;
  // myClasses.push_back(mine1);
  // mine1.x = 2.0; mine1.y = 2.1; mine1.z = 2.2;
  // myClasses.push_back(mine1);
  // mine1.x = 3.0; mine1.y = 3.1; mine1.z = 3.2;
  // myClasses.push_back(mine1);

  // std::vector<class myClass>::iterator it;
  // it = myClasses.begin();
  // mine1.x = 4.0; mine1.y = 4.1; mine1.z = 4.2;
  // myClasses.push_back(mine1);
  // myClasses.erase(it+1);

  // std::vector<class myClass>::iterator it2;
  // it2 = myClasses.begin();
  // while(it2!=myClasses.end()){

  //   myClass disp = *it2;

  //   printf("%lf %lf %lf\n",disp.x,disp.y,disp.z);
  //   it2++;
  // }




  return 0;
}
