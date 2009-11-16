/* 
 * Copyright 2009 Rensselaer Polytechnic Institute
 * This program is free software; you can redistribute it and/or modify 
 * it under the terms of the GNU General Public License as published by 
 * the Free Software Foundation; either version 2 of the License, or 
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but 
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY 
 * or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License 
 * for more details.
 * 
 * You should have received a copy of the GNU General Public License along 
 * with this program; if not, write to the Free Software Foundation, Inc., 
 * 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 */

#include "EM_Project_3_cpp_3D_comp.h"
#include <stdlib.h> /* malloc, free */



//This function is used to initialized some parameter.. It is called from the main function
void Initialize_Parameters(std::vector<std::vector<double> > *X,int** SEEDS,double** U,double*** Segma,double* PI,double** Z,int num_points,int num_components)
{
    bool converged = false;
    double thresh = 1.0;
    double d;
    double tmp;
    int ind_min;
    double u1, u2,u3;
    double ZZ_SUM;
	ind_min = 0;
    
    //Initialize vars and set the initial means to the seeds
    double* z_sum = (double *) malloc(num_components*sizeof(double));
    double** z_x_sum = (double **) malloc(num_components*sizeof(double*));    
    for(int i=0; i<num_components; i++)
    {
        z_sum[i] = 0.0;
        z_x_sum[i] = (double *) malloc(3*sizeof(double));
        for(int j=0; j<3; j++)
        {
            z_x_sum[i][j] = 0.0;
            U[i][j] = SEEDS[i][j];
        }
    }

    //The main convergence loop	
    while(!converged)
    {			
        //reset Z, z_sum, and z_x_sum
        for(int j=0; j<num_components; j++)
        {            
            z_sum[j] = 0.0;
            z_x_sum[j][0] = 0.0;
            z_x_sum[j][1] = 0.0;
            z_x_sum[j][2] = 0.0;
            for(int i=0; i<num_points; i++)
                Z[i][j] = 0.0;
                
        }
        //Assign each point to one component of U
        for(int i=0; i<num_points; i++)
        {
            d=1000.0;
            for(int j=0; j<num_components; j++)
            {				
                tmp = sqrt(pow(X->at(i).at(0)-U[j][0],2)+pow(X->at(i).at(1)-U[j][1],2)+pow(X->at(i).at(2)-U[j][2],2));
                if(tmp<d)
                {
                    d = tmp;
                    ind_min = j;
                }
            }           
            Z[i][ind_min] = X->at(i).at(3);
            z_sum[ind_min] = z_sum[ind_min] + X->at(i).at(3);   
            z_x_sum[ind_min][0] +=(Z[i][ind_min]*X->at(i).at(0));
            z_x_sum[ind_min][1] +=(Z[i][ind_min]*X->at(i).at(1));
            z_x_sum[ind_min][2] +=(Z[i][ind_min]*X->at(i).at(2));
        }
        
        //Recalculate U
        d = 0.0;
        tmp = 0.0;
        for(int i=0; i<num_components; i++)
        {
            u1 = z_x_sum[i][0]/z_sum[i];
            u2 = z_x_sum[i][1]/z_sum[i];
            u3 = z_x_sum[i][2]/z_sum[i];
            tmp = sqrt(pow(u1-U[i][0],2)+pow(u2-U[i][1],2)+pow(u3-U[i][2],2));
            U[i][0]=u1;
            U[i][1]=u2;
            U[i][2]=u3;
            if(tmp>d)
                d = tmp;
        }
        //check for convergence
        if(d<thresh)
            converged = true;
    }
    //once done with U, compute the initial values of the other parameters
    ZZ_SUM = 0.0; //added
    for(int i=0; i<num_components; i++) //added
        ZZ_SUM = ZZ_SUM + z_sum[i]; //added
    for(int i=0; i<num_components; i++)
    {
        PI[i] = (double) z_sum[i] / ZZ_SUM;//was / num_points
        for(int j=0; j<num_points; j++)
        {
            Segma[i][0][0]+= (Z[j][i]*pow(X->at(j).at(0)-U[i][0],2));
            Segma[i][1][1]+= (Z[j][i]*pow(X->at(j).at(1)-U[i][1],2));
            Segma[i][2][2]+= (Z[j][i]*pow(X->at(j).at(2)-U[i][2],2));
            Segma[i][0][1]+= (Z[j][i]*(X->at(j).at(0)-U[i][0])*(X->at(j).at(1)-U[i][1]));            
            Segma[i][0][2]+= (Z[j][i]*(X->at(j).at(0)-U[i][0])*(X->at(j).at(2)-U[i][2]));            
            Segma[i][1][2]+= (Z[j][i]*(X->at(j).at(1)-U[i][1])*(X->at(j).at(2)-U[i][2]));            
        }        
        
        Segma[i][0][0] = Segma[i][0][0]/z_sum[i];
        Segma[i][1][1] = Segma[i][1][1]/z_sum[i];
        Segma[i][2][2] = Segma[i][2][2]/z_sum[i];
        Segma[i][0][1] = Segma[i][0][1]/z_sum[i];
        Segma[i][0][2] = Segma[i][0][2]/z_sum[i];
        Segma[i][1][2] = Segma[i][1][2]/z_sum[i];
        Segma[i][1][0] = Segma[i][0][1];
        Segma[i][2][0] = Segma[i][0][2];
        Segma[i][2][1] = Segma[i][1][2]; //here is the first error: switch 2 and 1...FIXED        
    }

	//free memory	
	free(z_sum);
	for(int i=0; i<num_components; i++)
		free(z_x_sum[i]);
	free(z_x_sum);
}

//A 3-D multivariate gaussian
double Multivar_Norm3D(double X, double Y, double Z, double Ux, double Uy, double Uz, double S00, double S01, double S11, double S02, double S12, double S22)
{
    double det_segma = S00*(S11*S22-S12*S12) - S01*(S01*S22-S12*S02) + S02*(S01*S12-S11*S02);
    double Sinv00 = (S22*S11-S12*S12)/det_segma;
    double Sinv01 = -(S22*S01-S12*S02)/det_segma;
	double Sinv02 = (S12*S01-S11*S02)/det_segma;
    double Sinv11 = (S22*S00-S02*S02)/det_segma;
	double Sinv12 = -(S12*S00-S01*S02)/det_segma;
	double Sinv22 = (S11*S00-S01*S01)/det_segma;
	X = X-Ux;
	Y = Y-Uy;
	Z = Z-Uz;
    double Mah_Dist = X*(X*Sinv00+Y*Sinv01+Z*Sinv02) + Y*(X*Sinv01+Y*Sinv11+Z*Sinv12) + Z*(X*Sinv02+Y*Sinv12+Z*Sinv22);
    double Ex = exp(-.5*Mah_Dist);
    double A = 1/(pow((2*3.1416),1.5)*sqrt(det_segma));
    
    return (A*Ex);
    
}
//The first step of the EM: E-step
void Post_Expect(double **Gamma, std::vector<std::vector<double> > *X,double *PI,double **U,double ***Segma,int num_points,int num_components)
{
    double A_sum;
    double *A = (double *) malloc(num_components*sizeof(double));
    for(int j=0; j<num_components; j++)
        A[j]=0;
    
    for(int i=0; i<num_points; i++)
    {
        A_sum=0.0;
        for(int j=0; j<num_components; j++)
        {
			//Here.. we have to see if the point is too far from closest seeds, then set the prob A[j] to epsilon
			//double D = sqrt(pow(X->at(i).at(0)-U[j][0],2) + pow(X->at(i).at(1)-U[j][1],2) + pow(X->at(i).at(2)-U[j][2],2));
			//if(D>MAX_SZ)
			//	A[j] = numeric_limits<long double>::epsilon();
			//else
			//{
				A[j] = PI[j]*Multivar_Norm3D(X->at(i).at(0),X->at(i).at(1),X->at(i).at(2),U[j][0],U[j][1],U[j][2],Segma[j][0][0],Segma[j][1][0], Segma[j][1][1],Segma[j][0][2],Segma[j][1][2],Segma[j][2][2]);
			//	if(A[j] < numeric_limits<long double>::epsilon())
			//		A[j] = numeric_limits<long double>::epsilon();
			//}

            A_sum += A[j];
        }
        for(int k=0; k<num_components; k++)
        {
            Gamma[i][k] = X->at(i).at(3)*(A[k]/A_sum);  //here another chage.. multiply by the intensity X[i][3]
        }
    }
    free(A);
}

//The second step of the EM: M-step
void Param_Maximization(std::vector<std::vector<double> > *X,double **Gamma,double** U,double*** Segma,double* PI, int num_points, int num_components)
{
    double numpnts; //
	double gam_sum; //..
	double gam_x_sum[3];//..
    for(int k=0; k<num_components; k++)
    {
		gam_sum = 0.0;//..
		gam_x_sum[0] = gam_x_sum[1] = gam_x_sum[2] = 0.0;//..
		Segma[k][0][0] = 0;//..
        Segma[k][0][1] = 0;//..
		Segma[k][0][2] = 0;//..
        Segma[k][1][0] = 0;//..
        Segma[k][1][1] = 0;//..
		Segma[k][1][2] = 0;//..
		Segma[k][2][0] = 0;//..
		Segma[k][2][1] = 0;//..
        Segma[k][2][2] = 0;//..
        numpnts = 0.0;
        for(int i=0; i<num_points; i++)
        {
            numpnts = numpnts + X->at(i).at(3);//            
			gam_sum += Gamma[i][k];
            gam_x_sum[0] += (Gamma[i][k]*X->at(i).at(0));//..
            gam_x_sum[1] += (Gamma[i][k]*X->at(i).at(1));//..
            gam_x_sum[2] += (Gamma[i][k]*X->at(i).at(2));//..
			Segma[k][0][0]+= (Gamma[i][k]*pow(X->at(i).at(0)-U[k][0],2));//..
            Segma[k][1][1]+= (Gamma[i][k]*pow(X->at(i).at(1)-U[k][1],2));//..
            Segma[k][2][2]+= (Gamma[i][k]*pow(X->at(i).at(2)-U[k][2],2));//..
            Segma[k][0][1]+= (Gamma[i][k]*(X->at(i).at(0)-U[k][0])*(X->at(i).at(1)-U[k][1]));//..                                                
            Segma[k][0][2]+= (Gamma[i][k]*(X->at(i).at(0)-U[k][0])*(X->at(i).at(2)-U[k][2]));//..           
            Segma[k][1][2]+= (Gamma[i][k]*(X->at(i).at(1)-U[k][1])*(X->at(i).at(2)-U[k][2]));//..            
        }
        
		U[k][0] = gam_x_sum[0]/gam_sum;//..
        U[k][1] = gam_x_sum[1]/gam_sum;//..
        U[k][2] = gam_x_sum[2]/gam_sum;//..
		Segma[k][0][0] = Segma[k][0][0]/gam_sum;//..
        Segma[k][0][1] = Segma[k][0][1]/gam_sum;//..
        Segma[k][1][0] = Segma[k][0][1];//..
        Segma[k][1][1] = Segma[k][1][1]/gam_sum; //..       
        Segma[k][2][2] = Segma[k][2][2]/gam_sum; //..       
        Segma[k][0][2] = Segma[k][0][2]/gam_sum;//..
		Segma[k][2][0] = Segma[k][0][2];//..
        Segma[k][1][2] = Segma[k][1][2]/gam_sum; //..       
        Segma[k][2][1] = Segma[k][2][1]; //..
        
        PI[k] = gam_sum / numpnts;//..
    }   
}

void Compute_Probabilities(std::vector<std::vector<double> > *X, int num_points, double** U,double*** Segma,double* PI)
{
    //long int index = 0;
    for(int i=0; i<num_points; i++)
	{		
		X->at(i).at(4) = PI[0]*Multivar_Norm3D(X->at(i).at(0),X->at(i).at(1),X->at(i).at(2),U[0][0],U[0][1],U[0][2],Segma[0][0][0],Segma[0][1][0], Segma[0][1][1],Segma[0][0][2],Segma[0][1][2],Segma[0][2][2]);	
        X->at(i).at(5) = PI[1]*Multivar_Norm3D(X->at(i).at(0),X->at(i).at(1),X->at(i).at(2),U[0][0],U[0][1],U[0][2],Segma[1][0][0],Segma[1][1][0], Segma[1][1][1],Segma[1][0][2],Segma[1][1][2],Segma[1][2][2]);                    		
    }
}


void EM_Gmm(std::vector<std::vector<double> > *X,int** SEEDS,int num_points, int num_components)
{   
    //declare the different variables
    double** U;
    double** U_old;
    double*** Segma;
    double*** Segma_old;
    double*  PI;
    double*  PI_old;
    double Du, Ds;
    double** Gamma;                
    
    //Initialize
    U = (double **) malloc(num_components*sizeof(double*));   
    U_old = (double **) malloc(num_components*sizeof(double*));   
    for(int i=0; i<num_components; i++)
    {
        U[i] = (double *) malloc(3*sizeof(double));
        U_old[i] = (double *) malloc(3*sizeof(double));
        for(int j=0; j<3; j++)
        {
            U[i][j] = 0;
            U_old[i][j] = 0;
        }
    }
    Segma = (double ***) malloc(num_components*sizeof(double**));   
    Segma_old = (double ***) malloc(num_components*sizeof(double**));   
    for(int i=0; i<num_components; i++)
    {
        Segma[i] = (double **) malloc(3*sizeof(double*));
        Segma_old[i] = (double **) malloc(3*sizeof(double*));
        for(int j=0; j<3; j++)
        {
            Segma[i][j] = (double *) malloc(3*sizeof(double));
            Segma_old[i][j] = (double *) malloc(3*sizeof(double));
            for(int k=0; k<3; k++)
            {
                Segma[i][j][k] = 0;            
                Segma_old[i][j][k] = 0;            
            }
        }
    }
    PI = (double *) malloc(num_components*sizeof(double));   
    PI_old = (double *) malloc(num_components*sizeof(double));   
    for(int i=0; i<num_components; i++)
    {
        PI[i] = 0;   
        PI_old[i] = 0;   
    }
           
    Gamma = (double **) malloc(num_points*sizeof(double*));
    for(int i=0; i<num_points; i++)
    {        
        Gamma[i] = (double *) malloc(num_components*sizeof(double));
        for(int j=0; j<num_components; j++)
        {            
            Gamma[i][j]=0;
        }
    }
    //Call the function that uses K-means algorithm to compute initial estimates
    //of the parameters U[],Segma[][]and PI[] as well as the assigment matrix Z
    Initialize_Parameters_V2(X,SEEDS,U,Segma,PI,Gamma,num_points,num_components); //Change by Yousef on 04-03-2008: use the same matrix for Gamma and Z
	
    //The EM steps
    bool Converged = false;     
    int iter = 0;
    while(!Converged)
    {  
        iter++;
        if(iter>15)
            break;
        
        //copy the current parameters
        for(int i=0; i<num_components; i++)
        {
            PI_old[i] = PI[i];
            for(int j=0; j<3;j++)
            {
                U_old[i][j] = U[i][j];
                for(int k=0; k<3; k++)
                    Segma_old[i][j][k] = Segma[i][j][k];
            }
        }
        //1-E step
        //using the estimated parameter values, compute the initial values of the
        //posterior probabilities for Zn,k=1: Gamma[][]
        //Gamma = Post_Expect(X,PI,U,Segma,num_points, num_components);
		Post_Expect(Gamma,X,PI,U,Segma,num_points, num_components);
        //2-M step
        //using the computed Gamma, recomute the estimated values of the
        //parameters  
        Param_Maximization_V2(X,Gamma,U,Segma,PI, num_points, num_components);    
        //3- Check for convergence
        Du = 0;
        Ds = 0;
        double tmpU, tmpS;
        for(int i=0; i<num_components; i++)
        {            
            for(int j=1; j<3;j++)
            {
                tmpU = abs((int)U_old[i][j]-(int)U[i][j]);
                if(tmpU>Du)
                    Du = tmpU;
                for(int k=0; k<3; k++)
                {
                    tmpS = abs((int)Segma_old[i][j][k]-(int)Segma[i][j][k]);
                    if(tmpS>Ds)
                        Ds = tmpS;
                }
            }
        }  
        if(Du<2 && Ds<6)
            Converged = true; 
					
    }
 
    //Compute_Probabilities(X, num_points, U, Segma, PI);
	AssignToComponent(X, num_points, num_components, U,Segma,PI);

	//free memory	
	for(int i=0; i<num_components; i++)
    {
		free(U[i]);	
		free(U_old[i]);
		free(Gamma[i]);
        for(int j=0; j<3; j++)        
		{
            free(Segma[i][j]);
			free(Segma_old[i][j]);
		}
	}
	free(PI);
	free(PI_old);
	free(U);
	free(U_old);
	free(Segma);
	free(Segma_old);
}

//The code after this line was added by Yousef on 9/6/2009
//These are different versions of the functions above, but with the means (U) fixed at the seed points
//while the covariance matrices and P_I are estimated iteratively
//This function is used to initialized the parameter
void Initialize_Parameters_V2(std::vector<std::vector<double> > *X,int** SEEDS,double** U,double*** Segma,double* PI,double** Z,int num_points,int num_components)
{
    //bool converged = false;
    //double thresh = 1.0;
    double d;
    double tmp;
    int ind_min;
    //double u1, u2,u3;
    double ZZ_SUM;
	ind_min = 0;
    
    //Initialize vars and set the initial means to the seeds
    double* z_sum = (double *) malloc(num_components*sizeof(double));
    double** z_x_sum = (double **) malloc(num_components*sizeof(double*));    
    for(int i=0; i<num_components; i++)
    {
        z_sum[i] = 0.0;
        z_x_sum[i] = (double *) malloc(3*sizeof(double));
        for(int j=0; j<3; j++)
        {
            z_x_sum[i][j] = 0.0;
            U[i][j] = SEEDS[i][j];
        }
    }

    //The main convergence loop	
    //while(!converged)
    //{			
        //reset Z, z_sum, and z_x_sum
        for(int j=0; j<num_components; j++)
        {            
            z_sum[j] = 0.0;
            z_x_sum[j][0] = 0.0;
            z_x_sum[j][1] = 0.0;
            z_x_sum[j][2] = 0.0;
            for(int i=0; i<num_points; i++)
                Z[i][j] = 0.0;
                
        }
        //Assign each point to one component of U (1-NN classifier)
        for(int i=0; i<num_points; i++)
        {
            d=1000.0;
            for(int j=0; j<num_components; j++)
            {				
                tmp = sqrt(pow(X->at(i).at(0)-U[j][0],2)+pow(X->at(i).at(1)-U[j][1],2)+pow(X->at(i).at(2)-U[j][2],2));
                if(tmp<d)
                {
                    d = tmp;
                    ind_min = j;
                }
            }           
            Z[i][ind_min] = X->at(i).at(3);
            z_sum[ind_min] = z_sum[ind_min] + X->at(i).at(3);   
            z_x_sum[ind_min][0] +=(Z[i][ind_min]*X->at(i).at(0));
            z_x_sum[ind_min][1] +=(Z[i][ind_min]*X->at(i).at(1));
            z_x_sum[ind_min][2] +=(Z[i][ind_min]*X->at(i).at(2));
        }
        
    //    //Recalculate U
    //    d = 0.0;
    //    tmp = 0.0;
    //    for(int i=0; i<num_components; i++)
    //    {
    //        u1 = z_x_sum[i][0]/z_sum[i];
    //        u2 = z_x_sum[i][1]/z_sum[i];
    //        u3 = z_x_sum[i][2]/z_sum[i];
    //        tmp = sqrt(pow(u1-U[i][0],2)+pow(u2-U[i][1],2)+pow(u3-U[i][2],2));
    //        U[i][0]=u1;
    //        U[i][1]=u2;
    //        U[i][2]=u3;
    //        if(tmp>d)
    //            d = tmp;
    //    }
    //    //check for convergence
    //    if(d<thresh)
    //        converged = true;
    //}
    //Compute the initial values of the other parameters (Segma and P_I)
    ZZ_SUM = 0.0; //added
    for(int i=0; i<num_components; i++) //added
        ZZ_SUM = ZZ_SUM + z_sum[i]; //added
    for(int i=0; i<num_components; i++)
    {
        PI[i] = (double) z_sum[i] / ZZ_SUM;//was / num_points
        for(int j=0; j<num_points; j++)
        {
            Segma[i][0][0]+= (Z[j][i]*pow(X->at(j).at(0)-U[i][0],2));
            Segma[i][1][1]+= (Z[j][i]*pow(X->at(j).at(1)-U[i][1],2));
            Segma[i][2][2]+= (Z[j][i]*pow(X->at(j).at(2)-U[i][2],2));
            Segma[i][0][1]+= (Z[j][i]*(X->at(j).at(0)-U[i][0])*(X->at(j).at(1)-U[i][1]));            
            Segma[i][0][2]+= (Z[j][i]*(X->at(j).at(0)-U[i][0])*(X->at(j).at(2)-U[i][2]));            
            Segma[i][1][2]+= (Z[j][i]*(X->at(j).at(1)-U[i][1])*(X->at(j).at(2)-U[i][2]));            
        }        
        
        Segma[i][0][0] = Segma[i][0][0]/z_sum[i];
        Segma[i][1][1] = Segma[i][1][1]/z_sum[i];
        Segma[i][2][2] = Segma[i][2][2]/z_sum[i];
        Segma[i][0][1] = Segma[i][0][1]/z_sum[i];
        Segma[i][0][2] = Segma[i][0][2]/z_sum[i];
        Segma[i][1][2] = Segma[i][1][2]/z_sum[i];
        Segma[i][1][0] = Segma[i][0][1];
        Segma[i][2][0] = Segma[i][0][2];
        Segma[i][2][1] = Segma[i][1][2]; //here is the first error: switch 2 and 1...FIXED        
    }

	//free memory	
	free(z_sum);
	for(int i=0; i<num_components; i++)
		free(z_x_sum[i]);
	free(z_x_sum);
}

//The second step of the EM: M-step
void Param_Maximization_V2(std::vector<std::vector<double> > *X,double **Gamma,double** U,double*** Segma,double* PI, int num_points, int num_components)
{
    double numpnts; //
	double gam_sum; //..
	double gam_x_sum[3];//..
    for(int k=0; k<num_components; k++)
    {
		gam_sum = 0.0;//..
		gam_x_sum[0] = gam_x_sum[1] = gam_x_sum[2] = 0.0;//..
		Segma[k][0][0] = 0;//..
        Segma[k][0][1] = 0;//..
		Segma[k][0][2] = 0;//..
        Segma[k][1][0] = 0;//..
        Segma[k][1][1] = 0;//..
		Segma[k][1][2] = 0;//..
		Segma[k][2][0] = 0;//..
		Segma[k][2][1] = 0;//..
        Segma[k][2][2] = 0;//..
        numpnts = 0.0;
        for(int i=0; i<num_points; i++)
        {
            numpnts = numpnts + X->at(i).at(3);//            
			gam_sum += Gamma[i][k];
            gam_x_sum[0] += (Gamma[i][k]*X->at(i).at(0));//..
            gam_x_sum[1] += (Gamma[i][k]*X->at(i).at(1));//..
            gam_x_sum[2] += (Gamma[i][k]*X->at(i).at(2));//..
			Segma[k][0][0]+= (Gamma[i][k]*pow(X->at(i).at(0)-U[k][0],2));//..
            Segma[k][1][1]+= (Gamma[i][k]*pow(X->at(i).at(1)-U[k][1],2));//..
            Segma[k][2][2]+= (Gamma[i][k]*pow(X->at(i).at(2)-U[k][2],2));//..
            Segma[k][0][1]+= (Gamma[i][k]*(X->at(i).at(0)-U[k][0])*(X->at(i).at(1)-U[k][1]));//..                                                
            Segma[k][0][2]+= (Gamma[i][k]*(X->at(i).at(0)-U[k][0])*(X->at(i).at(2)-U[k][2]));//..           
            Segma[k][1][2]+= (Gamma[i][k]*(X->at(i).at(1)-U[k][1])*(X->at(i).at(2)-U[k][2]));//..            
        }
        
		//U[k][0] = gam_x_sum[0]/gam_sum;//..
  //      U[k][1] = gam_x_sum[1]/gam_sum;//..
  //      U[k][2] = gam_x_sum[2]/gam_sum;//..
		Segma[k][0][0] = Segma[k][0][0]/gam_sum;//..
        Segma[k][0][1] = Segma[k][0][1]/gam_sum;//..
        Segma[k][1][0] = Segma[k][0][1];//..
        Segma[k][1][1] = Segma[k][1][1]/gam_sum; //..       
        Segma[k][2][2] = Segma[k][2][2]/gam_sum; //..       
        Segma[k][0][2] = Segma[k][0][2]/gam_sum;//..
		Segma[k][2][0] = Segma[k][0][2];//..
        Segma[k][1][2] = Segma[k][1][2]/gam_sum; //..       
        Segma[k][2][1] = Segma[k][2][1]; //..
        
        PI[k] = gam_sum / numpnts;//..
    }   
}

void AssignToComponent(std::vector<std::vector<double> > *X, int num_points, int num_components, double** U,double*** Segma,double* PI)
{
    //compute the likelihood that a point belong to each gaussian component and assign to the max
	//the assignement value overwrites the LoG response of the point just to save memory	
    for(int i=0; i<num_points; i++)
	{		
		double max_prob = -1;
		for(int j=0; j<num_components; j++)
		{
			double prob = PI[j]*Multivar_Norm3D(X->at(i).at(0),X->at(i).at(1),X->at(i).at(2),U[j][0],U[j][1],U[j][2],Segma[j][0][0],Segma[j][1][0], Segma[j][1][1],Segma[j][0][2],Segma[j][1][2],Segma[j][2][2]);	
			if(prob>max_prob)
			{
				max_prob = prob;
				X->at(i).at(3) = j+1;
			}
		}
    }
}
