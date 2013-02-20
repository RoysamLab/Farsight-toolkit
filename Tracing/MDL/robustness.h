
#include <iostream>
#include <algorithm>
#include <vector>
using namespace std;

bool myfunction (float i,float j) { return (i<j); }

struct myclass {
  bool operator() (float i,float j) { return (i<j);}
} myobject;


int test (void);
double TrimmedMean(float a[], int elements, double r);
int test (void) {
  float myints[] = {32,71,12,45,26,80,53,33};
 /* vector<int> myvector (myints, myints+8);               // 32 71 12 45 26 80 53 33
  vector<int>::iterator it;

  // using default comparison (operator <):
  sort (myvector.begin(), myvector.begin()+4);           //(12 32 45 71)26 80 53 33

  // using function as comp
  sort (myvector.begin()+4, myvector.end(), myfunction); // 12 32 45 71(26 33 53 80)

  // using object as comp
  sort (myvector.begin(), myvector.end(), myobject);     //(12 26 32 33 45 53 71 80)

  // print out content:
  cout << "myvector contains:";
  for (it=myvector.begin(); it!=myvector.end(); ++it)
    cout << " " << *it;

  cout << endl;
*/
  double Tmean;
  Tmean=TrimmedMean(myints,7,0.1);
  printf("Tmean = %f  ", Tmean);
  return 0;
}

double TrimmedMean(float a[], int elements, double r)
{
   //printf("hello!");
   int i =0;
   printf("%d",elements);
   float *b = new float[elements];

   for (int i = 0; i < elements; ++i) 
   { //printf("%f ", a[i]);
     b[i] = a[i];
   }

  vector<float> myvector (b, b+elements);
  vector<float>::iterator it;
  int semi = int ((elements)/2);
  sort (myvector.begin(), myvector.begin()+semi);           //(12 32 45 71)26 80 53 33
  // using function as comp
  sort (myvector.begin()+semi, myvector.end(), myfunction); // 12 32 45 71(26 33 53 80)
  // using object as comp
  sort (myvector.begin(), myvector.end(), myobject);     //(12 26 32 33 45 53 71 80)
  i=0;
  for (it=myvector.begin(); it!=myvector.end(); ++it)
  {
   b[i]=*it;
   i++;
   //printf("%f ",b[i]);
   
  }
  // printf("%f  \n",b[0]);
  // printf("%f  \n",b[elements-1]);
   int g = 0;
   if(r>=0 && r<0.5)
     g = int (r*elements);
   else g = 0;
   
   double Tmean =0;
   for (i=g;i<elements-g;i++)
	   Tmean +=(double)b[i];
   Tmean /=(double)(elements-2*g);  
   return Tmean;
  
}


double WinsorizedMean(float a[], int elements, double r)
{
   //printf("hello!");
   int i =0;
   printf("%d",elements);
   float *b = new float[elements];

   for (int i = 0; i < elements; ++i) 
   { //printf("%f ", a[i]);
     b[i] = a[i];
   }

  vector<float> myvector (b, b+elements);
  vector<float>::iterator it;
  int semi = int ((elements)/2);
  sort (myvector.begin(), myvector.begin()+semi);          
  sort (myvector.begin()+semi, myvector.end(), myfunction); 
  sort (myvector.begin(), myvector.end(), myobject);     
  i=0;
  for (it=myvector.begin(); it!=myvector.end(); ++it)
  {
   b[i]=*it;
   i++;
  }
   int g = 0;
   if(r>=0 && r<0.5)
     g = int (r*elements);
   else g = 0;
   
   double Wmean =0;
   for (i=0;i<elements;i++)
   {   if(i<g)
	     Wmean +=(double)b[g-1];
       if (i>elements-g)
         Wmean +=(double)b[elements-g];
	   else 
		 Wmean +=b[i];
   }
   Wmean /=(double)(elements);  
   return Wmean;
  
}

