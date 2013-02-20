#include "ClusClus.h" //Dendrogram

ClusClus::ClusClus(vtkSmartPointer<vtkTable> TableReceivedFromSampleEditor)
{
	this->seed = 2;
	this->num_trials = 20;
	this->num_gaps = 20;
	this->naive_plot_flag = 1;
	table = TableReceivedFromSampleEditor;
	ClusterProgressTable = vtkSmartPointer<vtkTable>::New();
	this->Cluster_Hierarchical4_Master_("new_Test.txt");
	
	this->Mergers_To_Progress_((this->num_data),this->num_feats,"mergers.txt");
	this->OptimalLeafOrder1D_("mergers.txt");
	this->Prepare_MATLAB_Trees_("mergers.txt");
	/*std::cout<< "The ClusterProgressTable is created:"<<std::endl; 
	ClusterProgressTable->Dump(3);
	*/
	

}


void ClusClus::Cluster_Hierarchical4_Master_(char *root)
{	
	int i=0;
	Determine_File_Chars__(root,&num_data,&num_feats);
	/* If this class is being called from a class passing a VTKTbale then this function need not be invoked*/
	//ConvertTxtVtk((std::string)root,table,num_data,num_feats+2);
	
	//std::cout << "I am going to dump the file"<<std::endl;
	//this->table->Dump(3);

	std::cout<<"The file has" <<"\t"<<num_data <<"\t"<<"lines and" <<"\t"<<num_feats<<"\t"<< "features ."<<std::endl;
	std::cout<<"There are "<<num_data<<"data and"<< num_feats<<"features in"<< "Test_data.txt"<<std::endl;
	std::cout<<"********************************************"<<std::endl;
	std::cout<<"*   Enter the link mode:                   *"<<std::endl;
	std::cout<<"*      1 -- average link                   *"<<std::endl;
	std::cout<<"*      2 -- single link                    *"<<std::endl;
	std::cout<<"*      3 -- complete link                  *"<<std::endl;
	std::cout<<"********************************************"<<std::endl;
	std::cout<<"        Enter link mode:   ";
	std::cin >> link_mode;
	std::cout<<std::endl;
	std::cout<<"      Do you want a GAP statistic  (y/n)?   ";
	this->gap_flag = My_Yes_No_Querry_();
	if (this->gap_flag == 1) 
		{
			std::cout<<"      Enter number of trials  (e.g., 10):   ";
			std::cin>>num_trials;
		}
	else 
		{
			std::cout<<"      Do you want to estimate # clusters?   ";
			this->naive_plot_flag = My_Yes_No_Querry_();
		}
	std::cout<<std::endl;
	if (num_data-1 < num_gaps) num_gaps = num_data-1;
	MatrixAllocate_(&my_data, num_data, num_feats+3); // feature are 14 here
	MatrixAllocate_(&thy_data, num_data, num_feats+3); // feature are 14 here
	MatrixAllocate_(&metrix, num_data, 2);
	MatrixAllocate_(&GAP, num_data, 8);
	MatrixAllocate_(&HASTER, num_data, num_trials);
	MatrixAllocate_(&mergers, num_data, 5);//last column is distance
	VectorAllocate_(&disposition, num_trials);
	Zero_Matrix_(num_data, num_trials, HASTER);
	Zero_Matrix_(num_data, 2, metrix);
	Zero_Matrix_(num_data, 8, GAP);
	Zero_Matrix_(num_data, 5, mergers);
	/*std::cout<<"The value of the num_feats+3 is"<<num_feats+3<<std::endl;*/
	///////////////// Function to create the Data_Matrix ///////////////////
	CreateDataMatrix(table,my_data);
	/*table->Dump(3);*/

	/*for(int i=0;i<num_data-1;i++)
	{
	  printf("\n");
	  for(int j=0;j<num_feats+2;j++)
	  {
		 printf("%lf",my_data[i][j]);
		 printf("\t");
	  }
	}
	*/
	Cluster_Hierarchical4_Slave_(num_data, num_feats, my_data, metrix, mergers, link_mode, 1);
	mergers[num_data-1][4] = mergers[num_data-2][4];
	Print_Matrixx_(num_data, 5, mergers, "mergers.txt");
	/*std::cout<<std::endl<<"The value stored in merges.txt is :"<<std::endl;*/
	/*for(int i=0;i<num_data-1;i++)
	{
		std::cout<<std::endl;
		for(int j=0;j<5;j++)
		{
			std::cout << mergers[i][j]<<"\t";	
		}
	}*/

	Bolivar_Sort_(num_data, 2, 0, metrix, 1);
	Print_Matrixx_(num_data, 2, metrix, "metrix.txt");
	ratio = 2.0*metrix[0][1]/mergers[num_data-1][4]/3.0;
	
	for (k = 0; k < num_data-1; k++) 
	{
    GAP[k][0] = metrix[k][0];
    GAP[k][1] = metrix[k][1];
    GAP[k][7] = ratio*mergers[num_data-k-1][4];
	}
	
	for (k = 0; k < num_data-1; k++) 
	{
    GAP[k][2] = GAP[k+1][1]-GAP[k][1];
    GAP[k][5] = GAP[k][7]-GAP[k+1][7]; //dummy entry for closest cluster
	}
	
	piv_index = -999;
	mem = -999;
	/*std::cout <<"The value of num_gaps is :"<<num_gaps<<std::endl;*/
	for (k = 0; k < num_gaps; k++) 
	{
		if (GAP[k][5] > mem) 
			{
				mem = GAP[k][5];
				piv_index = k+1;
			}
	}
	std::cout << "See whether i reached here or not"<<std::endl;
	mem = 0;
	std::cout<<std::endl<<std::endl<<"***************************************************"<<std::endl;
	std::cout<<"      *  There are"<<"\t"<<piv_index <<"\t"<< "clusters (cluster merger index)   *"<<std::endl;
	std::cout<<"***************************************************"<<std::endl;
	
	if (gap_flag == 1) 
	{
		std::cout<<std::endl;
		for (l = 0; l < num_trials; l++)
			{
			Generate_Synthetic_Data_(num_data, num_feats, seed+l-1, my_data, thy_data);
			dispersion = Cluster_Hierarchical4_Slave_(num_data, num_feats, thy_data, metrix, mergers, link_mode, -1);
			std::cout<<"      Working on syntetic clustering iteration"<<l+1<<"\t"<<"dispersion =" <<dispersion;
			Bolivar_Sort_(num_data, 2, 0, metrix, 1);
			for (i = 0; i < num_data; i++) 
			HASTER[i][l] = metrix[i][1];
			}
		for (i = 0; i < num_gaps; i++) 
		{
			for (l = 0; l < num_trials; l++) 
				disposition[l] = HASTER[i][l];
			Stand_Devv_(num_trials, &avg, &std, disposition);
			GAP[i][3] = avg;
		}
		off_set = GAP[0][3]-GAP[0][1];
		for (i = 0; i < num_gaps; i++) GAP[i][3] = GAP[i][3] - off_set;       //align plots
		for (i = 0; i < num_gaps; i++) GAP[i][4] = GAP[i][3] - GAP[i][1];     //get GAP
		for (i = 0; i < num_gaps-1; i++) GAP[i][5] = GAP[i+1][4] - GAP[i][4]; //gap difference
		for (i = 0; i < num_gaps-2; i++) GAP[i][6] = GAP[i+1][5] - GAP[i][5]; //gap difference
		for (k = 0; k < num_gaps-2; k++) 
		{
			if (fabs(GAP[k][6]) > mem) 
			{
				if (GAP[k][6] < 0) 
				{
				mem = fabs(GAP[k][6]);
				piv_index = (int) GAP[k][0]+1;
				}
			}
		}
		std::cout<<std::endl<<std::endl<<"*****************************************"<<std::endl;
		std::cout<<"      *  There are"<<"\t"<<piv_index<< "clusters" <<mem<<std::endl;
		std::cout<<"      *****************************************"<<std::endl;
  }
	std::cout <<std::endl<<std::endl<<"   Wrote mergers.txt, GAP.txt, and GAP_plot.plt"<<std::endl;
	Print_Matrixx_(num_gaps, 8, GAP, "GAP.txt");
	Bolivar_Sort_(num_data, 8, 6, GAP, 1);
	/*for (i = 0; i < num_gaps; i++) 
	{
		if (GAP[i][6] < 0)
			std::cout<<"         k = "<<(int) GAP[i][0] <<  "  delta_dat_delta = "<<GAP[i][6]<<std::endl;
	}*/
  //toc();
	/*this->Mergers_To_Progress_(num_data-1,num_feats,"mergers.txt")*/; ///Converts the mergers.txt to progress.txt
	if (naive_plot_flag == 1) 
	{
    GAP1_plot_(gap_flag);
    SCommand_("gnuplot GAP_plot.plt", "", "", "", "");
	}
	free(disposition);
	MatrixFree_(mergers, num_data);
	MatrixFree_(HASTER, num_data);
	MatrixFree_(GAP, num_data);
	MatrixFree_(metrix, num_data);
	MatrixFree_(thy_data, num_data);
	MatrixFree_(my_data, num_data);


}
/*********************************************************/
/* Module #802: Optimal leaf order for 1D dendorgram     */
/* 	          Created by: Chris Gatti (gattic@rpi.edu)	 */
/*		          Last revised: 20.MAY.2011				 */
/*              use as <doe matlab_chris.txt 802>        */
/*********************************************************/
void ClusClus::OptimalLeafOrder1D_(char *root)
{

	int      i, num_data, num_feats, k;
	VECTOR   Tnums, pickedup, kids;
	MATRIX   T;

	std::cout<<"   Determine optimal leaf ordering for 1D dendrogram"<<std::endl;
	Determine_File_Chars__(root, &num_data, &num_feats); //get num_data and num_feats from mergers
	printf("   There are %3d data and %3d features in %s\n", num_data, num_feats, root);
	num_data = num_data + 1;                           //num_data is actually number of leaves
	VectorAllocate_(&Tnums, num_data-1);
	VectorAllocate_(&pickedup, num_data);
	VectorAllocate_(&kids, num_data);
	MatrixAllocate_(&T, num_data-1, 4);
  //Read data files
	Read_Meta_Data_(num_data-1, 4, T, root, 0);         //read tree.txt
  //Initialize variables
	k = 2*(num_data-1);                                //maximum cluster number
	for (i = 0; i < num_data-1; i++) Tnums[i] = num_data+i;
	for (i = 0; i < num_data; i++) pickedup[i] = -1;
	for (i = 0; i < num_data; i++) kids[i] = -1;
  //Get kids
	GetKids_(num_data, k, Tnums, pickedup, kids, T);
	Print_IVector_(num_data, kids, "optimal_leaf_order_1D.txt");
	std::cout<<std::endl<<"    Wrote optimal_leaf_order_1D.txt (a vector)"<<std::endl;
	MatrixFree_(T, num_data-1);
	free(kids);
	free(pickedup);
	free(Tnums);
}

void ClusClus::Print_IVector_(int num_data, VECTOR y, char *root) 
{
	int    i;
	FILE   *fpotv;

	fpotv = FDeclare_(root, "", 'w');
	for (i = 0; i < num_data; i++) fprintf(fpotv, "%7d\n", (int) y[i]);
	fclose(fpotv);
}

//Accessory function to OptimalLeafOrder1D: Recursively works through tree to pick up leaf nodes
void ClusClus::GetKids_(int num_data, int k, VECTOR Tnums, VECTOR pickedup, VECTOR kids, MATRIX T) 
{
	int      i = 0, j = 0, k1, k2;

	if (k < num_data) 
	{					   //if node number is a leaf
		while (j == 0) 
		{       //find first open position in pickedup
			if (kids[i] == -1) 
			{
			j = 1;
			kids[i] = k;
			pickedup[i] = k;   //may not need this
			}
		i++;
		}
	}
  else 
  {                   //if node number is not a leaf
	for (j = 0; j < num_data-1; j++) 
	{
      if (Tnums[j] == k) i = j;
	}                     //find index of k in Tnums
	k1 = (int)T[i][0];
    k2 = (int)T[i][1];
	 GetKids_(num_data, k1, Tnums, pickedup, kids, T); //pick up kids from node k1
	 GetKids_(num_data, k2, Tnums, kids, kids, T);     //pick up kids from node k2
  }
}

/***********************************************************/
/*		Optimal leaf order for 1D dendrogram {PMT}		   */
/*													       */
/*											               */
/***********************************************************/
void ClusClus::Prepare_MATLAB_Trees_(char *root) {
   Mergers_2_Chris_(root);
  
   OptimalLeafOrder1D_("matlab_chris.txt");
}

/*********************************************************************/
/*      Convert mergers.txt to MATLAB format for dendrogram			 */
/* 	          Created by Chris Gatti (gattic@rpi.edu)	    		 */
/*              use as <doe mergers.txt 801>                         */
/*		          Last revised: 20.MAY.2011						     */
/*********************************************************************/
void ClusClus::Mergers_2_Chris_(char *root) 
{
	int      i, ii, j, jj, jjlast, indj, k, kk, kklast, indk;
	int 	   found, found1, found2, num_data, num_feats;
	VECTOR   distances;
	MATRIX   mergers, matlab_chris, progress;

	std::cout<<"   Convert mergers.txt to MATLAB format for dedrogram"<<std::endl;
	std::cout<<"   Requires progress.txt (created using operator 388)"<<std::endl;
	Determine_File_Chars__(root, &num_data, &num_feats);
	std::cout<<"   There are"<<"\t"<<num_data<<"\t"<<"data and"<<num_feats<< "features in "<<root<<std::endl;
  //Allocate memory
	MatrixAllocate_(&mergers, num_data, 5);           //mergers data
	MatrixAllocate_(&matlab_chris, num_data-1, 4);    //matrix for new tree
	MatrixAllocate_(&progress, num_data, num_data+1); //matrix for progress
	VectorAllocate_(&distances, num_data);            //cluster distances
  //Read mergers.txt and progress.txt
	Read_Meta_Data_(num_data, 5, mergers, root, 0);                     //read mergers.txt
	Read_Meta_Data_(num_data, num_data+1, progress, "progress.txt", 0); //read progress.txt
  //Initialize
	for (i = 0; i < num_data; i++) distances[i] = mergers[i][4];       //get distances
	for (i = 0; i < num_data-1; i++) 
	{
    matlab_chris[i][0] = -1; //placeholder for cluster number, all real cluster numbers are non-neg
    matlab_chris[i][1] = -1; //placeholder for cluster number, all real cluster numbers are non-neg
	}
  //main loop: work through columns of progress
	for (i = 0; i < (num_data-1); i++) 
	{
		j = (int)mergers[i][1];                        //low cluster number
		k = (int)mergers[i][2];                        //high cluster number
	 //Find parent for cluster number j
		for (ii = 0; ii < num_data; ii++) 
			if ((int)progress[ii][i+1] == j) 
				indj = ii;
	 //find row in matlab_chris that includes indj
		jj = -1;
		for (ii = 0; ii < num_data-1; ii++) 
			{
			if (matlab_chris[ii][0] == indj) jj = ii;
			if (matlab_chris[ii][1] == indj) jj = ii;
			}
	 //find ultimate parent --------------------
		if (jj == -1) 
			jjlast = indj;                  //if not in tree, add new cluster number
		else 
		{ 								                 //if currently in tree, walk through tree
		found = -1; //found1 = -1; found2 = -1;
			while (found == -1) 
			{
			jj = num_data + jj;
			jjlast = jj;                              //save jj in jjlast incase jj is not currently in tree
			//look in tree, see if jj is in either column
			found1 = FindRowIndex_(num_data, 0, jj, matlab_chris); //find jj in matlab_chris col 0
			found2 = FindRowIndex_(num_data, 1, jj, matlab_chris); //find jj in matlab_chris col 1
			if (found1 != -1) jj = found1;                        //if found, save index
			if (found2 != -1) jj = found2;                        //if found, save index
			if ((found1 == -1) && (found2 == -1)) found = -2;     //if not found in either, exit loop
		}
    }
    //Find parent for cluster number k --------
	 for (ii = 0; ii < num_data; ii++) 
		 if ((int)progress[ii][i+1] == k) 
			 indk = ii;
    //find row in matlab_chris that includes indk
	 kk = -1;
	 for (ii = 0; ii < num_data-1; ii++) 
	 {
	   if (matlab_chris[ii][0] == indk) 
		   kk = ii;
	   if (matlab_chris[ii][1] == indk) 
		   kk = ii;
	 }
	 //find ultimate parent ----------------------
    if (kk == -1) kklast = indk; 				                //if not in tree, add new cluster number
    else { 									                         //if currently in tree, walk through tree
	   found = -1; //found1 = -1; found2 = -1;
	   while (found == -1) {
		  kk = num_data + kk;
        kklast = kk;                            //save kk in kklast incase kk is not currently in tree
		  //look in tree, see if kk is in either column
		  found1 = FindRowIndex_(num_data, 0, kk, matlab_chris);//find kk in matlab_chris col 0
		  found2 = FindRowIndex_(num_data, 1, kk, matlab_chris);//find kk in matlab_chris col 1
		  if (found1 != -1) kk = found1;                       //if found, save index
        if (found2 != -1) kk = found2;                       //if found, save index
		  if ((found1 == -1) && (found2 == -1)) found = -2;    //if not found in either, exit loop
	   }
    }
    //add cluster numbers to matlab_chris ------
	 matlab_chris[i][0] = Minimilian_(jjlast, kklast);
	 matlab_chris[i][1] = Maximilian_(jjlast, kklast);
  }
  for (ii = 0; ii < num_data-1; ii++) 
  {
	 matlab_chris[ii][2] = distances[ii];                    //add distance column
	 matlab_chris[ii][3] = num_data+ii;                      //add cluster number
  }
  Print_IMatrixx_(num_data-1, 4, matlab_chris, "matlab_chris.txt");
  
  std::cout<<std::endl<<"      Wrote matlab_chris.txt"<<std::endl<<std::endl;
  free(distances);
  MatrixFree_(progress, num_data);
  MatrixFree_(matlab_chris, num_data-1);
  MatrixFree_(mergers, num_data);
}

//returns the maximum of two numbers
double ClusClus::Maximilian_(double a1, double a2) 
{
	double diogenes;

	diogenes = a2;
	if (a1 > a2) diogenes = a1;
	return diogenes;
}

//returns the minimum of two numbers
double ClusClus::Minimilian_(double a1, double a2) 
{
	double diogenes;

	diogenes = a2;
	if (a1 < a2) diogenes = a1;
	return diogenes;
}

// Find row index f which matches c in column col of matrix A
// Only for unique value of c within column
// Returns index f (non-negative) if c is found, or -1 if c is not found
int ClusClus::FindRowIndex_(int num_data, int col, int c, MATRIX A)
{
	int     i, ii;

	ii = -1;
	for (i = 0; i < num_data-1; i++) 
	{
		if (A[i][col] == c) ii = i;
	}
	return ii;
}

vtkSmartPointer<vtkTable> ClusClus::GetProgressTable()
{
return this->ClusterProgressTable;
}
ClusClus::~ClusClus()
{
}

void ClusClus::Stand_Devv_(int num_data, double *avg, double *std, VECTOR my_vector) 
{
	int      i;
	double   tempd, sum;

	sum = 0.0;
	for (i = 0; i < num_data; i++) sum += my_vector[i];
	*avg = sum/num_data;
	sum = 0.0;
	for (i = 0; i < num_data; i++) 
	{
    tempd = my_vector[i]-*avg;
    sum = sum + tempd*tempd;
	}
	*std = sqrt(sum/(num_data-1));
}
void ClusClus::Mergers_To_Progress_(int num_data, int num_feats, char *root) 
{
	int      i, k, pivot1, pivot2;
	MATRIX   mergers, progress;
	
	/*std::cout<<"The number of features are num_feats :"<<num_feats<<std::endl;
	std::cout<<"The number of features are num_data :"<<num_data<<std::endl;*/
	std::cout<<"The file to be opened is :"<<root<<std::endl;
	std::cout<<"   Convert mergers.txt to progress.txt"<<num_data<< "clusters"<<std::endl;
	if (strcmp(root, "mergers.txt")!= 0) 
	{
    printf("!!!Has to operate on a mergers.txt file\n");
	}
	MatrixAllocate_(&mergers, num_data, 5);
	MatrixAllocate_(&progress, num_data, num_data+1);
	Zero_Matrix_(num_data, num_data+1, progress);
	Read_Meta_Data_(num_data, 5, mergers, "mergers.txt", 0);  
	/*std::cout<< "The data stored inside the mergers is :"<<std::endl;
	for(int i=0;i<num_data;i++)
	{
		std::cout<<std::endl;
		for(int j=0;j<5;j++)
		{
			std::cout << mergers[i][j]<<"\t";	
		}
	}
	*/
	
	////////////Directly reading the mergers.txt and converting it to the ClusterProgressTable//////////////

	for (i = 0; i < num_data; i++) 
	{
    progress[i][0] = i;
    progress[i][1] = i;
	}
	for (k = 0; k < num_data-1; k++) 
	{
    pivot1 = (int) mergers[k][1];
    pivot2 = (int) mergers[k][2];
    //printf("k = %2d\tpivot1 = %d\tpivot2=%d\n", k, pivot1, pivot2);getche();
		for (i = 0; i < num_data; i++) 
			progress[i][k+2] = progress[i][k+1];
		for (i = 0; i < num_data; i++) 
			{
				if (progress[i][k+2] == pivot2) 
					progress[i][k+2] = pivot1;
				if (progress[i][k+2] > pivot2)
					progress[i][k+2] -= 1;
			}
    //for (i = 0; i < num_data-1; i++) printf("%3d\t", (int) progress[i][k+2]);
    //printf("%3d\n", (int) progress[num_data-1][k+2]);
	}
	this->ClusterProgressTable->Initialize();

	std::cout <<"The Elements inside the progress.txt is "<<std::endl;
	for(int i=0;i<num_data;i++)
	{
		std::cout<<std::endl;
		for(int j=0;j<num_data+1;j++)
		{
			std::cout << progress[i][j]<<"\t";	
		}
	}
	
	for ( unsigned int i = 0; i <num_data+1; i++ )
	{
		vtkSmartPointer<vtkVariantArray> col = 
		vtkSmartPointer<vtkVariantArray>::New();
			for ( unsigned int j=0;j<num_data;j++)
				{
				col->InsertNextValue ( vtkVariant ( progress[j][i] ) );
				}
		this->ClusterProgressTable->AddColumn ( col );
    }
	/*std::cout <<"The progress vtkTable so created is :"<<std::endl;
ClusterProgressTable->Dump(3);*/

	Print_IMatrixx_(num_data, num_data+1, progress, "progress.txt");
	std::cout<<"   Wrote progress.txt"<<std::endl;
	MatrixFree_(progress, num_data);
	MatrixFree_(mergers, num_data);
}
void ClusClus::Read_Meta_Data_(int num_data, int num_feats, MATRIX my_data, char *root,int mode) 
{
	int      i, j;
	double    templf;
	FILE     *fpin;

	fpin = FDeclare_(root, "", 'r');
	if (mode == 1) fscanf(fpin, "%lf", &templf);// not yet understood
	for (i = 0; i < num_data; i++) 
	{
		for (j = 0; j < num_feats; j++) 
			{
			fscanf(fpin, "%lf", &templf);
			my_data[i][j] = templf;
			}
	}
  fclose(fpin);
}


void ClusClus::Print_IMatrixx_(int num_data, int num_feats, MATRIX gerissa, char *root) 
{
	int      i, j;
	FILE     *fpotl;

	fpotl = FDeclare_(root, "", 'w');
	for (i = 0; i < num_data; i++) 
	{
		for (j = 0; j < num_feats-1; j++) 
			fprintf(fpotl, "%3d\t", (int) gerissa[i][j]);
		fprintf(fpotl, "%3d\n", (int) gerissa[i][num_feats-1]);
	}
	fclose(fpotl);
}
void ClusClus::GAP1_plot_(int gap_flag)
{
	FILE     *fpot;

	fpot = FDeclare_("GAP_plot", "plt", 'w');
	fprintf(fpot, "set title \"GAP INDEX PLOT\"\n");
	fprintf(fpot, "set terminal windows \"Helvetica Bold\" 18\n");
	fprintf(fpot, "set xlabel \"Number of CLUSTERS\"\n");
	fprintf(fpot, "set ylabel \"GAP\"\n");
	fprintf(fpot, "#set time\n");
	fprintf(fpot, "set autoscale\n");
	fprintf(fpot, "set parametric\n");
	fprintf(fpot, "#set grid\n");
	fprintf(fpot, "set key graph 0.18, 0.95\n");
	if (gap_flag == 1) 
	{
    fprintf(fpot, "plot 'GAP.txt' using 1:2 ti \"dispersion\" with lp 1,\\\n");
    fprintf(fpot, "'GAP.txt' using 1:4 ti \"synthetic\" with lp 2,\\\n");
    fprintf(fpot, "'GAP.txt' using 1:5 ti \"GAP\" with lp 3,\\\n");
    fprintf(fpot, "'GAP.txt' using 1:6 ti \"GAP_diff\" with lp 4\n");
	}
	else 
	{
    fprintf(fpot, "plot 'GAP.txt' using 1:2 ti \"dispersion\" with lp 1,\\\n");
    fprintf(fpot, "'GAP.txt' using 1:8 ti \"distance\" with lp 2\n");
	}
	fprintf(fpot, "replot\n");
	fprintf(fpot, "pause -1\n");
	fclose(fpot);
}
	
void ClusClus::SCommand_(char *string1, char *string2, char *string3, char *string4,char *string5) 
{
	char     command[80];

	strcpy(command, string1);
	strcat(command, string2);
	strcat(command, string3);
	strcat(command, string4);
	strcat(command, string5);
	system(command);
}
	
void ClusClus::Generate_Synthetic_Data_(int num_data, int num_feats, int seed,MATRIX my_data, MATRIX thy_data) 
{
	int      i, j;
	double   minn, maxx;

	srand(seed);
	for (j = 0; j < num_feats; j++) 
	{
    minn = 999;
    maxx = -999;
		for (i = 0; i < num_data; i++) 
		{
			if (my_data[i][j] < minn) minn = my_data[i][j];
			if (my_data[i][j] > maxx) maxx = my_data[i][j];
		}
    //for (i = 0; i < num_data; i++) thy_data[i][j] = minn+unifrand()*(maxx-minn);
		for (i = 0; i < num_data; i++) 
			thy_data[i][j] = minn+randomm_(32000)/32000.0*(maxx-minn);
	}
	for (i = 0; i < num_data; i++) 
	{
    thy_data[i][num_feats] = my_data[i][num_feats];
    thy_data[i][num_feats+1] = my_data[i][num_feats+1];
	}
}
int ClusClus::randomm_(int parameter) {
  return(rand()%parameter);
}


void ClusClus::MatrixFree_(MATRIX matrix, int nRows)
{
	int i;
	for (i = 0; i < nRows; i++)
    free(matrix[i]);
	free(matrix);
}
void ClusClus::Bolivar_Sort_(int numRows, int numColumns, int sort_column, MATRIX big, int mode) 
{
	int      i, j, index;
	MATRIX   thy_data;
	VECTOR   zz, yy, response;

	MatrixAllocate_(&thy_data, numRows, numColumns);
	VectorAllocate_(&yy, numRows);
	VectorAllocate_(&zz, numRows);
	VectorAllocate_(&response, numRows);
	for (i = 0; i < numRows; i++) 
	{
		for (j = 0; j < numColumns; j++) 
			thy_data[i][j] = big[i][j];
	}
	for (i = 0; i < numRows; i++) 
	{
    zz[i] = big[i][sort_column];
    response[i] = zz[i];
	}
	//Arrays.sort(zz);
	qsort(zz, numRows, sizeof(double), &compare_); //Quick Sort
	for (i = 0; i < numRows; i++) 
	{
		for (j = 0; j < numRows; j++) 
			{
			if (zz[i] == response[j]) 
				{
				yy[i] = (double) j;
				response[j]=3.14159;
				break;
				}
			}
	}
	if (mode <= 0) 
	{
		for (i = 0; i < numRows; i++) 
			zz[i] = yy[i];
		for (i = 0; i < numRows; i++) 
			yy[i] = zz[numRows-1-i];
	}
	for (i = 0; i < numRows; i++) 
	{
		index = (int) yy[i];
		for (j = 0; j < numColumns; j++) 
			big[i][j] = thy_data[index][j];
	}
  free(response);
  free(zz);
  free(yy);
  MatrixFree_(thy_data, numRows);
}

int ClusClus::compare_(const void *a, const void *b) 
{
   //return( strcmp((double *) a, (double *) b));
   if (*(double *) a  < *(double *) b) return -1;
   if (*(double *) a == *(double *) b) return 0;
   if (*(double *) a  > *(double *) b) return 1;
   return 1;
}

void ClusClus::Print_Matrixx_(int num_data, int num_feats, MATRIX gerissa, char *root) 
{
	int      i, j;
	FILE     *fpotl;

	fpotl = FDeclare_(root, "", 'w');
	for (i = 0; i < num_data; i++) 
	{
		for (j = 0; j < num_feats-1; j++) 
			fprintf(fpotl, "%14.7e\t", gerissa[i][j]);
     fprintf(fpotl, "%14.7e\n", gerissa[i][num_feats-1]);
	}
	fclose(fpotl);
	
}


double ClusClus::Cluster_Hierarchical4_Slave_(int num_data, int num_feats, MATRIX my_data,MATRIX metrix, MATRIX mergers, int link_mode, int print_mode) 
{
	/*std::cout << "I have entered in the function Cluster_Hierarchical_Slave_"<<std::endl;
	std::cout << "The value of num_data is :"<<num_data<<std::endl;
	std::cout << "The value of num_feats is :"<<num_feats<<std::endl;
	std::cout << "The value of link_mode is :"<<link_mode<<std::endl;

	std::cout << "The value of my_data is :"<<std::endl;*/
	/*for(int i=0;i<num_data;i++)
	{
		std::cout <<std::endl;
		for(int j=0;j<num_feats;j++)
			{
				std::cout<<my_data[i][j];
				std::cout<<"\t";
			}
	}
	*/
	 /*std::cout << "I have reached here"<<std::endl;*/
	int        i, it, num_protos, pivot1 = 0, pivot2 = 0;

	double     dispersion;

	VECTOR     num_cluster_data, data_distancer;

	VECTOR     proto_distancer;

	if(print_mode==1)
		{
			for(int i=0;i<num_data;i++)
				{					
					printf("\n");
					for(int j=0;j<num_feats+2;j++)
						{
						printf("%lf\t",my_data[i][j]);
						}
				}
		}
  VectorAllocate_(&proto_distancer, num_data*(num_data+1)/2);

  VectorAllocate_(&data_distancer, num_data*(num_data+1)/2);

  Zero_Vector_(num_data*(num_data+1)/2, data_distancer);

  num_protos = num_data;

  //initialize clusters, initially one at each datapoint

  for (i = 0; i < num_data; i++) {

    my_data[i][num_feats+2] = i;

    metrix[i][0] = num_data-i;

    mergers[i][0] = i;

  }

  Data_Distances_V_(num_data, num_feats, my_data, data_distancer);// It calculates the dissimilarity matrix based on Euclidean distance\
  
  

   if(print_mode==1)
	  {

	  printf("The Euclidean distance matrix is:\n\n\n\t");
  for(int i=0;i<num_data*(num_data+1)/2;i++)
	  printf("%lf\t",data_distancer[i]);
	  }

  for (it = 0; it < num_data; it++) {

    VectorAllocate_(&num_cluster_data, num_protos);

    if (it == 0) {


      dispersion = Calculate_Proto_Distances_V_(num_protos, num_data, num_feats,

        num_cluster_data, my_data, proto_distancer, data_distancer, link_mode);
	 
	  //printf("The dispersion value inside slave is :%lf \n",dispersion);
    }
    else {
      dispersion = Update_Proto_Distances_V_(num_protos, num_data, num_feats,
        pivot1, pivot2, num_cluster_data, my_data, proto_distancer,
        data_distancer, link_mode);
    }
    if (print_mode == 1)

      printf("      It number = %2d   num_clusters = %2d   dispersion = %8.3lf\r\n", it, num_protos, dispersion);

    //calculate distances between prototypes, merge clusters in dataspace & reduce # clusters in protos

    mergers[it][4] = Merge_ZE_Clusters_V_(num_data, &num_protos, num_feats, my_data,

      proto_distancer, &pivot1, &pivot2);

    mergers[it][1] = pivot1;

    mergers[it][2] = pivot2;

    mergers[it][3] = dispersion;

    metrix[it][1] = dispersion;

    free(num_cluster_data);

  }

  free(data_distancer);

  free(proto_distancer);
 
  return dispersion;

}

	   
double ClusClus::Merge_ZE_Clusters_V_(int num_data, int *num_clusters, int num_feats,MATRIX my_data, VECTOR cluster_distancer, int *pivot1, int *pivot2) 
{
	int      i, k, kk, pirot, num_protos;
	double   minn;

	num_protos = *num_clusters;
	minn = 1.0e10;
  /*
  for (k = 0; k < num_protos; k++) {
    for (kk = k+1; kk < num_protos; kk++) { //assumed lower triangle form and swapped indices
      if (cluster_distancer[kk*(kk+1)/2+k] < minn) {
        minn = cluster_distancer[kk*(kk+1)/2+k];
        *pivot2 = kk;
        *pivot1 = k;
      }
    }
  }
  */
  for (k = 0; k < num_protos; k++) {
    for (kk = 0; kk < k; kk++) { //assumed lower triangle form and swapped indices
      if (cluster_distancer[k*(k+1)/2+kk] < minn) {
        minn = cluster_distancer[k*(k+1)/2+kk];
        *pivot2 = k;
        *pivot1 = kk;
      }
    }
  }
  for (i = 0; i < num_data; i++) {
    pirot = (int) my_data[i][num_feats+2];
    if (pirot == *pivot2) my_data[i][num_feats+2] = *pivot1;
    else if (pirot > *pivot2) my_data[i][num_feats+2] = pirot-1;
  }
  *num_clusters = num_protos-1;
  return minn;
}
double ClusClus::Update_Proto_Distances_V_(int num_clusters, int num_data, int num_feats,int pirot1, int pirot2, VECTOR num_cluster_data, MATRIX my_data,VECTOR cluster_distancer, VECTOR data_distancer, int link_mode) 
{
int        i, ik, il, k, l, num_data_k, num_data_l, counter, pivot1, pivot2;
double     temp_distance, dispersion, minn, maxx, raincheck;
VECTOR     my_data_k_ID, my_data_l_ID;

  //updata by reusing
  //both conditions in the inner loop ar not obvious and very tricky
  for (k = 0; k < num_clusters+1; k++) {
    for (l = pirot2; l < k; l++) { //can by speeded up by l < k
      cluster_distancer[k*(k+1)/2+l] = cluster_distancer[k*(k+1)/2+l+1];
    }
  }
  for (k = pirot2; k < num_clusters; k++) {
    for (l = 0; l < k+1; l++) { //can be speeded up by l < k+1
      cluster_distancer[k*(k+1)/2+l] = cluster_distancer[(k+1)*(k+2)/2+l];
    }
  }
  //make new prototypes quick by updating prototypes
  //calculate number of data in clusters quickly
  for (i = 0; i < num_data ; i++) {
    pivot1 = (int) my_data[i][num_feats+2];
    num_cluster_data[pivot1] += 1;

  }
  //update by recalculating cl_dis[pirot1][pirots] and cl_dis[pirot1][pirot]
  //calculate distances between clusters
  num_data_k = (int) num_cluster_data[pirot1];
  VectorAllocate_(&my_data_k_ID, num_data_k);
  Zero_Vector_(num_data_k, my_data_k_ID);
  counter = 0;
  for (i = 0; i < num_data; i++) {
    if ((int) my_data[i][num_feats+2] == pirot1) {
      my_data_k_ID[counter] = i;
      counter += 1;
    }
  }
  for (l = 0; l < num_clusters; l++) {
    minn =  1.0e10;
    maxx = -1.0e10;
    num_data_l = (int) num_cluster_data[l];
    VectorAllocate_(&my_data_l_ID, num_data_l);
    Zero_Vector_(num_data_l, my_data_l_ID);
    counter = 0;
    for (i = 0; i < num_data; i++) {
      if ((int) my_data[i][num_feats+2] == l) {
        my_data_l_ID[counter] = i;
        counter += 1;
      }
    }
    temp_distance = 0.0;
    for (ik = 0; ik < num_data_k; ik++) {
      pivot1 = (int) my_data_k_ID[ik];
      for (il = 0; il < num_data_l; il++) {
        pivot2 = (int) my_data_l_ID[il];
        if (pivot1 > pivot2) raincheck = data_distancer[pivot1*(pivot1+1)/2+pivot2];
        else raincheck = data_distancer[pivot2*(pivot2+1)/2+pivot1];
        temp_distance += raincheck;
        if (raincheck < minn) minn = raincheck; //if (l != pirot1) {
        if (raincheck > maxx) maxx = raincheck;
      }
    }
    if (link_mode == 1) raincheck = sqrt(temp_distance*temp_distance/num_data_l/num_data_k);
    else if (link_mode == 2) raincheck = minn;
    else if (link_mode == 3) raincheck = maxx;
    if (l == pirot1) cluster_distancer[pirot1*(pirot1+1)/2+l] = temp_distance*temp_distance/num_data_l/num_data_l;
    else if (l > pirot1) { //changed this for lower triangle mode (swapped l and pirot1 indices)
      cluster_distancer[l*(l+1)/2+pirot1] = raincheck;
    }
    else if (l < pirot1) {
      cluster_distancer[pirot1*(pirot1+1)/2+l] = raincheck;
    }
    free(my_data_l_ID);
  }
  free(my_data_k_ID);
  //for (l = 0; l < num_clusters; l++) //not necessary for lower triangle mode
  //  cluster_distancer[l][pirot1] = cluster_distancer[pirot1][l];
  dispersion = 0.0;
  for (k = 0; k < num_clusters; k++) dispersion += cluster_distancer[k*(k+1)/2+k]*num_cluster_data[k];
 /* printf("The dispersion inside update is:%lf\n",dispersion);*/
  return sqrt(dispersion/num_data);
}

void ClusClus::Zero_Vector_(int N, VECTOR A)
{
	int      i;
	for (i = 0; i < N; i++) 
		A[i] = 0.0;
}

	   
double ClusClus::Calculate_Proto_Distances_V_(int num_clusters, int num_data, int num_feats,VECTOR num_cluster_data, MATRIX my_data, VECTOR cluster_distancer,VECTOR data_distancer, int link_mode) 
{
int        i, ik, il, k, l, num_data_k, num_data_l, counter, pivot1, pivot2;
double     temp_distance, dispersion, minn, maxx, raincheck;
VECTOR     my_data_k_ID, my_data_l_ID;

  //Zero_Matrix(num_clusters, num_clusters, cluster_distancer);
  //calculate number of data in each cluster faster
  Zero_Vector_(num_clusters, num_cluster_data);
  for (i = 0; i < num_data ; i++) {
	 /* printf("The value of my_data[%d][%d] is:%lf \n",i,num_feats+2,my_data[i][num_feats+2]);*/
    pivot1 = (int) my_data[i][num_feats+2];
//printf("The value of pivot1 is :%d \n",pivot1);
    num_cluster_data[pivot1] += 1;
	
	/*printf("The value of num_cluster_data[%d] is :%lf \n",i,num_cluster_data[i]);*/
  }
  //calculate distances between clusters
  for (k = 0; k < num_clusters; k++) {
    num_data_k = (int) num_cluster_data[k];
	VectorAllocate_(&my_data_k_ID, num_data_k);
    Zero_Vector_(num_data_k, my_data_k_ID);
    counter = 0;
    for (i = 0; i < num_data; i++) {
      if ((int) my_data[i][num_feats+2] == k) {
        my_data_k_ID[counter] = i;
        counter += 1;
      }
    }
    for (l = 0; l < k+1; l++) {
      maxx = -1.0e10;
      minn = 1.0e10;
      num_data_l = (int) num_cluster_data[l];
      VectorAllocate_(&my_data_l_ID, num_data_l);
      Zero_Vector_(num_data_l, my_data_l_ID);
      counter = 0;
      for (i = 0; i < num_data; i++) {
        if ((int) my_data[i][num_feats+2] == l) {
          my_data_l_ID[counter] = i;
          counter += 1;
        }
      }
      temp_distance = 0.0;
      for (ik = 0; ik < num_data_k; ik++) {
        pivot1 = (int) my_data_k_ID[ik];
        for (il = 0; il < num_data_l; il++) {
          pivot2 = (int) my_data_l_ID[il];
          if (pivot1 > pivot2) raincheck = data_distancer[pivot1*(pivot1+1)/2+pivot2];
          else raincheck = data_distancer[pivot2*(pivot2+1)/2+pivot1];
          //data_distancer[(pivot1-1)*num_data+pivot2];
          temp_distance += raincheck;
          if (raincheck < minn) minn = raincheck; //if (l != k) {
          if (raincheck > maxx) maxx = raincheck;
        }
      }
      if (l == k) cluster_distancer[k*(k+1)/2+l] = temp_distance*temp_distance/num_data_k/num_data_l;
      else if (l > k) { //swaped k,l indices for lower triangle mode
        if (link_mode == 1) cluster_distancer[l*(l+1)/2+k] = sqrt(temp_distance*temp_distance/num_data_l/num_data_k);
        else if (link_mode == 2) cluster_distancer[l*(l+1)/2+k] = minn;
        else if (link_mode == 3) cluster_distancer[l*(l+1)/2+k] = maxx;
      }
      else if (l < k) {
        if (link_mode == 1) cluster_distancer[k*(k+1)/2+l] = sqrt(temp_distance*temp_distance/num_data_l/num_data_k);
        else if (link_mode == 2) cluster_distancer[k*(k+1)/2+l] = minn;
        else if (link_mode == 3) cluster_distancer[k*(k+1)/2+l] = maxx;
      }
      free(my_data_l_ID);
    }
    free(my_data_k_ID);
  }
  dispersion = 0.0;
  for (k = 0; k < num_clusters; k++) dispersion += cluster_distancer[k*(k+1)/2+k]*num_cluster_data[k];
 /* printf("The dispersion is :%lf \n",dispersion);*/
  return sqrt(dispersion/num_data);
}

void ClusClus::Data_Distances_V_(int num_data, int num_feats, MATRIX my_data, VECTOR data_distancer) 
{
	int        i, j, ii;
	VECTOR     AA, BB;

	VectorAllocate_(&AA, num_feats);
	VectorAllocate_(&BB, num_feats);
	for (i = 0; i < num_data; i++) 
	{
		for (j = 0; j < num_feats; j++) AA[j] = my_data[i][j];
			for (ii = 0; ii < i+1; ii++) 
			{
    //for (ii = 0; ii < num_data; ii++) {
			for (j = 0; j < num_feats; j++) BB[j] = my_data[ii][j];
			 data_distancer[i*(i+1)/2+ii] = Distance_(num_feats, AA, BB); 
			}
	}
  free(BB);
  free(AA);
}

double ClusClus::Distance_(int num_feats, VECTOR A, VECTOR B) 
{
int      j;
double   sum, distance;

  sum = 0.0;
  for (j = 0; j < num_feats; j++) sum += (A[j]-B[j])*(A[j]-B[j]);
  distance = sqrt(sum);
  return distance;
}

void ClusClus::CreateDataMatrix(vtkSmartPointer<vtkTable> table,MATRIX my_data)
{

	std::cout << "The number of rows are :"<<table->GetNumberOfRows()<<std::endl;
	std::cout << "The number of columns are :"<<table->GetNumberOfColumns()<<std::endl;

for(vtkIdType r = 0; r < table->GetNumberOfRows(); r++ )
  {
    for(vtkIdType c = 1; c < table->GetNumberOfColumns(); c++ )
    {
    vtkVariant v = table->GetValue( r,c);
	my_data[r][c-1]=v.ToDouble();
    }
  }
for(vtkIdType r = 0; r < table->GetNumberOfRows(); r++ )
	{
    vtkVariant v = table->GetValue( r,0);
	my_data[r][table->GetNumberOfColumns()-1]=v.ToDouble();
	}
std::cout <<"The data contained in the my_data matrix is now as follows"<<std::endl;
for(int i=0;i<num_data;i++)
  {
	  std::cout<<std::endl;
	  for(int j=0;j<num_feats+2;j++)
	  {
		  std::cout<<my_data[i][j];
		  std::cout<<"\t";
	  }
  }
std::cout<<std::endl;
  std::cout<<"I completed the printing of the my_data matrix"<<std::endl;
std::cout << "I came out of the CreateDataMatrix()"<<std::endl;
return;
}
void ClusClus::AllocateCols_(PFLOAT matrix[], int nRows, int nCols) 
{
	int i;
	for (i = 0; i < nRows; i++)
    VectorAllocate_(&matrix[i], nCols);
}

void ClusClus::Zero_Matrix_(int N, int M, MATRIX A) {
int i, j;
  for (i = 0; i < N; i++) {
    for (j = 0; j < M; j++) A[i][j] = 0.0;
  }
}

void ClusClus::MatrixAllocate_(MATRIX *pmatrix, int nRows, int nCols) {
  if ((*pmatrix=(MATRIX) calloc(nRows, sizeof(PFLOAT))) == NULL) {
    printf("Not enough memmry\n");
    exit(1);
   }
   AllocateCols_(*pmatrix, nRows, nCols);
}

int ClusClus::My_Yes_No_Querry_() 
{
	int      querry_flag;
	char     thresholdstring[80], nostring[80];
	std::cin >> thresholdstring;
	/*scanf("%s", &thresholdstring)*/;
	strcpy(nostring, "n");
	querry_flag = strcmp(thresholdstring, nostring);
	if (querry_flag != 0) querry_flag = 1;
  //printf("  querry_flag = %10d\n", querry_flag);
	return querry_flag;
}

void ClusClus::ConvertTxtVtk(std::string hname,vtkSmartPointer<vtkTable> table,int num_row,int num_col)
{	
	/*std::cout<<"The value of num_row is "<<num_row<<std::endl;
	std::cout<< "The file name is "<<hname<<std::endl;*/

	const int MAXLINESIZE = 102400;
	char line[MAXLINESIZE];
	
	
	ifstream read_file;
	this->table->Initialize();
	/*std::vector<std::vector<std::string> > file;
	file.resize(num_row);*/
	/*for(int i=0;i<num_row;i++)
		file[i].resize(num_col);*/


	read_file.open(hname.c_str());
	

	if ( !read_file.is_open() )
		return ;

	vtkSmartPointer<vtkDoubleArray> column = vtkSmartPointer<vtkDoubleArray>::New();
	read_file.getline(line, MAXLINESIZE);
	//while ( !read_file.eof() ) //Get all values
	//{
		std::string h;
		char * pch = strtok (line," \t");
		/*std::cout<<"The value of pch is "<<pch<<std::endl;*/
		while (pch != NULL)
		{
			h = pch;
			column = vtkSmartPointer<vtkDoubleArray>::New();
			column->SetName( h.c_str() );
			this->table->AddColumn(column);

			pch = strtok (NULL, " \t");
		}
	/*	read_file.getline(line, MAXLINESIZE);
	}
	*/
		/*table->Dump(3);*/
	read_file.getline(line,MAXLINESIZE);
	/*std::cout <<line<<std::endl;*/
	/*if(read_file.is_open())
	{*/
	//std::cout << line<<std::endl;
		/*std::cout << "The file exists"<<std::endl;
		std::cout << hname<<std::endl;*/
	int k=0;
	while ( !read_file.eof() ) //Get all values
	{
		vtkSmartPointer<vtkVariantArray> row = vtkSmartPointer<vtkVariantArray>::New();
		
		char * pch = strtok (line," \t");
		
		while (pch != NULL)
		{
			
			/*std::cout << "I am printing the value of pch"<<pch<<std::endl;*/
			row->InsertNextValue( vtkVariant( atof(pch) ) );
		
			pch = strtok (NULL, " \t");
		}
		/*std::cout << " I came out the while loop"<<std::endl;*/
			if(k<=num_row-2)
			this->table->InsertNextRow(row);
			/*table->Dump(3);*/
			k++;
			
		read_file.getline(line, MAXLINESIZE);
	}
	//}
	/*std::cout << "I am going to close the file pointer"<<std::endl; */
	read_file.close();
	

	/*if(read_file.is_open())
		{
		
			for ( unsigned int i = 0; i <num_col; i++ )
				{
				vtkSmartPointer<vtkVariantArray> col = 
				vtkSmartPointer<vtkVariantArray>::New();
				for ( unsigned int j=0;j<num_row;j++)
					{
						std::cout << getline(read_file,file[i][j])<<std::endl;
					col->InsertNextValue ( vtkVariant ( getline(read_file,file[i][j],' ')) );
					}
				std::cout << "The program crashes here"<<std::endl;
				table->AddColumn ( col );
				}



		}*/
	/*std::cout<<" The table inside the convertTxtVtk is :"<<std::endl;*/
	//this->table->Dump(3);

}
int ClusClus::Determine_File_Chars__(char *root, int *num_data, int *numfeats) {
int      i, num_features = -2, num_lines = 0, num_chars = 0;
char     ch = 0;
VECTOR   cum_sum;
FILE     *fpin;
printf("I have entered in the function Determine_File_Chars");
  fpin = FDeclare_(root, "", 'r');
  while (1) {
    ch = 1;
    fscanf(fpin, "%c", &ch);
    if (feof(fpin)) break;
    if (ch == '\n') num_lines += 1;
  }
  rewind(fpin);
  printf("!!!There are %2d lines\n", num_lines);
  ch = 1;
  while (ch!='\n') {
    fscanf(fpin, "%c", &ch);
    num_chars += 1;
  }
  printf("the number of characters are:%d\n",num_chars);
  rewind(fpin);
  VectorAllocate_(&cum_sum, num_chars);
  cum_sum[0] = 0;
  for (i = 1; i < num_chars; i++) {
    fscanf(fpin, "%c", &ch);
    if ((ch == ' ')||(ch == '\t')) cum_sum[i] = 0;
    else cum_sum[i] = 1;
  }
  fclose(fpin);
  for (i = 1; i < num_chars; i++) {
    if ((cum_sum[i-1]==0)&&(cum_sum[i]==1)) num_features +=1;
  }
  *numfeats = num_features;
  *num_data = num_lines;
  free(cum_sum);
  return num_features;
}

FILE * ClusClus::FDeclare_(char *root, char *extension, char key) {
char string1[80];
FILE *fpot;

  strcpy(string1, root);
  if (strcmp(extension,"")!=0) strcat(string1, ".");
  strcat(string1, extension);
  if ((key != 'w')&&(key!='r')) {
    printf("either read or write, wrong key %c \n", key);
    exit(1);
  }
  if (key == 'r') {
    //if (( (FILE *) fpot = fopen(string1, "r")) == NULL) {
    //changed for DOS?UNIX compatibility
    if ((fpot = fopen(string1, "r")) == NULL) {
      printf("\n");
      printf("input file %s does not exist\n", string1);
      exit(1); 
    }
  }
  if (key == 'w') {
    if ((fpot = fopen(string1, "w")) == NULL) {
      printf("\n");
      printf("output file %s does not exist\n", string1);
      exit(1);
    }
  }
  return fpot;
}
void ClusClus::VectorAllocate_(VECTOR *vector, int nCols) {
  if ((*vector=(VECTOR) calloc(nCols, sizeof(double))) == NULL) {
    printf("Not enough memmry\n");
    exit(1);
  }
}