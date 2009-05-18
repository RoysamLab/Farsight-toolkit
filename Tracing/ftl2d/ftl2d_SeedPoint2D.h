#ifndef SEEDPOINT2D_H
#define SEEDPOINT2D_H

class SeedPoint2D	{
	public :
	unsigned int getx() {return(x);}
	unsigned int gety() {return(y);}
	float getScale() {return(Scale);}
	float getIntensity() {return(Intensity);}
	
	void setx(unsigned int i) {x = i;}
	void sety(unsigned int i) {y = i;}
	void setScale(float i) {Scale = i;}
	void setIntensity(float i) {Intensity = i;}
	void PrintSelf();
	
	private:
	unsigned int x;
    unsigned int y;
    float Intensity;
    float Scale;
 };
 
 #endif
