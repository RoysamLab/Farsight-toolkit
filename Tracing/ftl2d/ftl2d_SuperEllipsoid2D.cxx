/*=========================================================================
Copyright 2009 Rensselaer Polytechnic Institute
Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License. 
=========================================================================*/

#include "ftl2d_SuperEllipsoid2D.h"

#define DEBUG 0

SuperEllipsoid2D::SuperEllipsoid2D(double d)	{
	MaxAspectRatio = d;
}

SuperEllipsoid2D::~SuperEllipsoid2D()	{
	/*delte these variables if they are not automatically deleted in ITK
	delete s;
	delete S;
	delete U;
	delete F;
	delete A;
	*/
}

void SuperEllipsoid2D::generateConvexHull(double& a1, double& a2, double& e1, Vect2& mu, double& q)	{

	Mat22 R;
	R(0,0) = cos(q); 	R(0,1) = sin(q); 
	R(1,0) = -sin(q);	R(1,1) = cos(q); 

	double T, rho;
	double x,y,x1,y1, rx, ry;
	double mx, my, mx_last, my_last, rx_last, ry_last;
	double tempx, tempy;
	double mx0, my0, rx0, ry0;

	step = 2*vnl_math::pi/static_cast<double>(NN+1);
//	std::cout<<"Convex Hull Testing "<<std::endl;
	for (unsigned int i=0;i<=NN;i++)	{
		T = i*step;
	
		x1 = cos(T);		this->U(i,0) = x1;
		y1 = sin(T);		this->U(i,1) = y1;

		x = x1*a1;		y = y1*a2;

		rho = pow(fabs(x1),2/e1) + pow(fabs(y1),2/e1);
		rho = pow(rho,-0.5*e1);

		rx = rho*x;	ry = rho*y;

		mx = rx*R(0,0) + ry*R(0,1) + mu(0);
		my = rx*R(1,0) + ry*R(1,1) + mu(1);
		
		//std::cout<<i<<": "<<rho <<"\t"<< rx <<" " << ry<<"\t" << mx <<" " << my <<std::endl;
		//compute the areas (A) and centers (S)
		if (i>0)	{
			tempx = mx_last - mx;
			tempy = my_last - my;
			this->A[i-1] = vnl_math_hypot(tempx,tempy);
			//std::cout<<i<<": "<<mx <<"\t"<< mx_last <<"\t" << my<<"\t" << my_last <<"\t"<<vnl_math_hypot(tempx,tempy)<<std::endl;
			
			this->s(i-1,0) = 0.5*(mx_last + mx);
			this->s(i-1,1) = 0.5*(my_last + my);

			this->S(i-1,0) = 0.5*(rx_last + rx);
			this->S(i-1,1) = 0.5*(ry_last + ry);

		}

		else	{
			rx0 = rx;
			ry0 = ry;
			mx0 = mx;
			my0 = my;
		}

		//for the last connection to complete the ellipse
		if (i==NN)	{
			tempx = mx-mx0;
			tempy = my-my0;
			this->A[i] = vnl_math_hypot(tempx,tempy);

			this->s(i,0) = 0.5*(mx + mx0);
			this->s(i,1) = 0.5*(my + my0);

			this->S(i,0) = 0.5*(rx + rx0);
			this->S(i,1) = 0.5*(ry + ry0);
		}

		rx_last = rx;
		ry_last = ry;
		mx_last = mx;
		my_last = my;
	
	}
	
//	for(int i=0;i<=NN;i++)
//		std::cout<<i<<"- U: "<<U(i,0)<<" "<<U(i,1)<<"  "<<"\ts: "<<s(i,0)<<" "<<s(i,1)<<"  " <</*"\tS: "<<S(i,0)<<" "<<S(i,1)<<"  "<<*/"\tA: "<<A(i)<<std::endl;
//	std::cout<<std::endl;
}

bool SuperEllipsoid2D::Interpolation(ImageType::Pointer im)	{

	float x,y, fx, fy, t1,t2;
	unsigned int px, py;
	ImageType::IndexType ndx11, ndx12, ndx21, ndx22;
	ImageType::RegionType::SizeType sz = im->GetBufferedRegion().GetSize();
	
	//std::cout<<"Interpolation Testing "<<std::endl;
	for (int i = 0;i<this->s.rows();i++)	{
		if (this->s(i,0) < 0) {return 0;}
		if (this->s(i,1) < 0) {return 0;}
		if (this->s(i,0) > sz[0]-1) {return 0;}
		if (this->s(i,1) > sz[1]-1) {return 0;}
		
		px = int(this->s(i,0));			py = int(this->s(i,1));
		fx = this->s(i,0) - float(px);	fy = this->s(i,1) - float(py);
		
		ndx11[0] = px;	ndx12[0] = px;		ndx21[0] = px+1;	ndx22[0] = px+1;
		ndx11[1] = py;	ndx12[1] = py+1;	ndx21[1] = py;		ndx22[1] = py+1;
		
		t1 = (1-fy)*im->GetPixel(ndx11) + fy*im->GetPixel(ndx12);
		t2 = (1-fy)*im->GetPixel(ndx21) + fy*im->GetPixel(ndx22);
		
		this->F[i] = static_cast<double>((1-fx)*t1 + (fx)*t2);
		
		//std::cout << i << "- " <<this->s(i,0)<<" "<<this->s(i,1) <<"\t"<< im->GetPixel(ndx22) <<std::endl;
		//std::cout << i << "- " <<this->s(i,0)<<" "<<this->s(i,1) <<"\t"<< this->F[i] <<std::endl;
	}
	//std::cout <<std::endl;
	return 1;
}

void SuperEllipsoid2D::UpdateLikelihood(float& Fg, float& Bg)	{

	double f_a = 0.5;
	
	for (unsigned int i = 0; i<=NN;  i++)	
		this->L[i] = ((1-f_a)*fabs((this->F[i]) - double(Fg))) - (f_a * fabs((this->F[i]) - double(Bg)));
	
	this->L /= (fabs(Fg-Bg));
	this->L = element_product(this->L, this->A);
	
//	std::cout << "LIKELIHOOD TEST " <<(fabs(Fg-Bg))<< std::endl;
//	for (unsigned int i = 0; i<=NN;  i++)
//		std::cout << this->L[i] << std::endl;
}

void SuperEllipsoid2D::gradient_superquad2d(double& a1, double& a2, double& e1, Vect2& mu, double& q)	{

	double x, y, x1, y1;
	double cosq, sinq;
	double dT1x, dT1y, dT2x, dT2y, dT1, dT2;  
	double Gx, Gy, MagG;
	cosq = cos(q);
	sinq = sin(q);
	
	//std::cout << "NORMAL TEST" <<std::endl;
	
	for (int i = 0; i <=NN; i++)	{
		x1 = (s(i,0) - mu[0]);
		y1 = (s(i,1) - mu[1]);
		x = (cosq*x1 - sinq*y1)/a1;
		y = (sinq*x1 + cosq*y1)/a2;

		dT1x = cosq/a1;
		dT1y = -sinq/a1;
		dT2x = sinq/a2;
		dT2y = cosq/a2;

		// 2.*(abs(x)).^((2./e1)-1)./e1.*sign(x)
		dT1 = 2*pow(fabs(x), (2/e1)-1)*(1/e1)* vnl_math_sgn(x);
		dT2 = 2*pow(fabs(y), (2/e1)-1)*(1/e1)* vnl_math_sgn(y);
		
		//std::cout << i << "- "<< pow(fabs(x), (2/e1)-1)*(1/e1)*(x/fabs(x)) << " " << pow(fabs(y), (2/e1)-1)*(1/e1)*(y/fabs(y)) << std::endl;
		//std::cout << i << "- "<< (1/e1)*(x/fabs(x)) << " " << (1/e1)*(y/fabs(y)) << std::endl;

		Gx = dT1*dT1x + dT2*dT2x + 0.00001;
		Gy = dT1*dT1y + dT2*dT2y + 0.00001;
		MagG = vnl_math_hypot(Gx,Gy);
		this->N(i,0) = Gx/MagG;
		this->N(i,1) = Gy/MagG;
		//std::cout<<	i <<" "<<MagG <<" " << pow(fabs(x),2/e1) + pow(fabs(y),2/e1) <<std::endl;
	}
	
	//std::cout << "NORMAL TEST" << std::endl;
	//for (unsigned int i = 0; i<=NN;  i++)
	//	std::cout << i << "- "<< this->N(i,0) << " " << this->N(i,1)<< std::endl;
}

void SuperEllipsoid2D::UpdateOrientation(double &q, unsigned int _debug)	{

	const double dt2 = 3;

	//vnl_matrix_fixed <double,NN+1,2> dq;
	vnl_vector_fixed <double,NN+1> dq1;

	Mat22 dR;
	dR(0,0) = -1*sin(q); 
	dR(0,1) = cos(q); 
	dR(1,0) = -1*cos(q); 
	dR(1,1) = -1*sin(q); 
	
	for (unsigned int i = 0; i<=NN;  i++)
		dq1[i] = (dR(0,0)*S(i,0) + dR(0,1)*S(i,1))*N(i,0) + (dR(1,0)*S(i,0) + dR(1,1)*S(i,1))*N(i,1);
		
	//dq = (this->S) * dR.transpose();
	//dq = element_product(dq, this->N);
	//dq1 = dq.get_column(0) + dq.get_column(1);

	q += -dt2*dot_product(dq1.normalize(), this-> L) / (this->A.sum()); //Suspicious factor 3
	
	if (_debug == 1)	{
		std::cout << "Orientation update " << -dt2*dot_product(dq1, this-> L)/(this->A.sum()) << "  " << (this->A.sum()) << std::endl;
		//std::cout << "dR= [" << dR <<"]"<<std::endl;
		//for (unsigned int i = 0; i<=NN;  i++)
		//	std::cout << i << "- "<< this->L[i]<<"\t" << dq1[i] <<"\t" << this->L[i]*dq1[i]  << std::endl;
		//std::cout << i << "- "<< this->S(i,0)<<"\t" << S(i,1) <<std::endl;
	}
		
}

void SuperEllipsoid2D::UpdateShapeParameter(double &e1, double &a1, double &a2, double &q)	{

	double x,y,a,b, de2, de1, sx, sy, ssx, ssy;
	double ssx_0, ssx_last, ssy_0, ssy_last;
	double eps = 0.0001;
	double ee = e1;
	double cosq, sinq;
	double delta_e1;
	vnl_matrix_fixed <double,NN+1,2> de;
	
	double dt2 = 3;
	double dt = 0.2;

	cosq = cos(q);
	sinq = sin(q);

	for (int i=0; i<=NN; i++)	{
		x = this->U(i,0);
		y = this->U(i,1);

		a = pow(fabs(x),2/ee);
		b = pow(fabs(y),2/ee);
		de2 = ((-2*a)/(ee*ee))*log(fabs(x+eps)) + ((-2*b)/(ee*ee))*log(fabs(y+eps)) ;
		de1 = pow((a+b),-0.5*ee) * (-0.5*log(a+b+eps) -0.5*ee*de2/(a+b));
		//std::cout << i << "- "<< de1 << std::endl;
		
		sx = a1*x*de1;
		sy = a2*y*de1;
		ssx =  cosq*sx + sinq*sy;
		ssy = -sinq*sx + cosq*sy;

		if (i>0)	{
			de(i-1,0) = 0.5*(ssx_last + ssx);
			de(i-1,1) = 0.5*(ssy_last + ssy);
		}

		else	{
			ssx_0 = ssx;
			ssy_0 = ssy;
		}

		if (i==NN)	{
			de(i,0) = 0.5*(ssx_0 + ssx);
			de(i,1) = 0.5*(ssy_0 + ssy);
		}

		ssx_last = ssx;
		ssy_last = ssy;
	}

	de = element_product(de, this->N);
	delta_e1 = -dt*dt * dot_product((de.get_column(0)+de.get_column(1)), this->L);
	
	//std::cout << "delta e1 "<< delta_e1 << std::endl;
	
	if ((e1 + delta_e1) > 1)
		delta_e1 = 0;

	e1 = e1 + delta_e1;
	e1 = vnl_math_max(e1,0.25);

}


void SuperEllipsoid2D::UpdateScale(double &q, double &a1, double &a2)	{
	
	double dt = 0.2;
	double delta_a1, delta_a2;
	double r_a1, r_a2, sum_a;
	
	double R11 = cos(q);
	double R21 = -sin(q);
	double R12 = sin(q);
	double R22 = cos(q);
	vnl_matrix_fixed <double,NN+1,2> da1, da2;

	for(int i=0; i<=NN; i++)	{
		da1(i,0) = R11*S(i,0)/a1;
		da1(i,1) = R21*S(i,0)/a1;
		da2(i,0) = R12*S(i,1)/a2;
		da2(i,1) = R22*S(i,1)/a2;
		
	}

	da1 = element_product(da1, this->N);
	da2 = element_product(da2, this->N);

	delta_a1 = -1*(dt*dt) * dot_product((da1.get_column(0)+da1.get_column(1)), this->L);
	delta_a2 = -1*(dt*dt) * dot_product((da2.get_column(0)+da2.get_column(1)), this->L);
	
	sum_a = a1+a2;
	r_a1 = a1/sum_a;
	r_a2 = a2/sum_a;
	
	r_a1 = 0.1*sum_a*(1/(1+this->MaxAspectRatio) - r_a1);
	r_a2 = 0.1*sum_a*(this->MaxAspectRatio/(1+this->MaxAspectRatio) - r_a2);
	
	a1 += (delta_a1 + r_a1);
	a2 += (delta_a2 + r_a2);
	//if (((a2+delta_a2 < this->MaxAspectRatio*a1) && (delta_a2 > 0)) || (delta_a2 < 0)) 
	//	a2 += delta_a2;
	
	a1 = fabs(a1);
	a2 = fabs(a2);
	a1 = vnl_math_max(1.5, a1);
	a2 = vnl_math_max(1.5, a2);
	//a2 = vnl_math_min(a2,MaxAspectRatio*a1);
	
	//std::cout << "SCALE UPDATE TESTING" << std::endl;
	//for(int i=0; i<=NN; i++)	
	//	std::cout << i << "- "<< S(i,0)<< "\t" << S(i,1) << std::endl;
	//std::cout << "delta a1 "<< delta_a1 << std::endl;
	//std::cout << "delta a2 "<< delta_a2 << std::endl;
		
	//this is the place where numerical deviation from the Matlab Implementation was first observed

}

void SuperEllipsoid2D::UpdatePosition(Vect2& mu, unsigned char IsInit, double &q)	{
	
	double dt = 0.2;
	Vect2 dmu;
	dmu[0] = -1*(dt*dt)*dot_product(this->N.get_column(0),this->L);
	dmu[1] = -1*(dt*dt)*dot_product(this->N.get_column(1),this->L);
	
	if (!IsInit)	{
		double sinq = sin(q);
		double cosq = cos(q);
		dmu[0] = dmu[0]*(1 - sinq*sinq + cosq*sinq);
		dmu[1] = dmu[1]*(1 - cosq*cosq + cosq*sinq);
	}
	
	mu += dmu;
	//std::cout << "delta MU "<< dmu << std::endl;
}

void SuperEllipsoid2D::UpdateAxis(double& a1, double& a2, double& q, Vect2& last_dir)	{
	double temp;
	if (a1 > 1.1*a2)	{
		temp = a1;
		a1 = a2;
		a2 = temp;
		q += vnl_math::pi/2;
	}

	Vect2 dir;
	dir[0] = sin(q);
	dir[1] = cos(q);
	if (dot_product(last_dir,dir) < 0)	{
		q = q + vnl_math::pi;
	}
}

void SuperEllipsoid2D::PrintSelf()	{
	std::cout<<std::endl<<"Printing Superellipsoid Parameters" <<std::endl;
	std::cout << "s\t\tS\t\tF\t\tN\t\tL\t\tA"<<std::endl;
	for (unsigned int i=0; i<=NN; i++)
		std::cout <<s[i]<<"\t"<<S[i]<<"\t"<<F[i]<<"\t"<<N[i]<<"\t"<<L[i]<<"\t"<<A[i]<<std::endl;
}
