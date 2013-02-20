
////////////////Code to create a dendrogram by reading agglomerative clustering information//////////////////////////////

#include "Dendrogram.h"


Dendrogram::Dendrogram(QWidget * parent, Qt::WindowFlags flags)
: QMainWindow(parent,flags)
{
	
	Level_Detected=0;
	num_row=0;
	colour_child=0;
	count=0;
	Selection= new ObjectSelection();//
	/*QVTKWidget widget;
	widget.resize(256,256);*/
	/*clock_t t1=clock();
	Determine_File_Chars("progress.txt",&num_row,&num_col);
	std::cout << "Total time taken for determining the num_row and num_col is: " << ( (clock() - t1)/(float) CLOCKS_PER_SEC) << endl;*/

	/*colour_child=(int *) calloc(num_row,sizeof(int));*/

	
	graph_Layout = vtkSmartPointer<vtkMutableUndirectedGraph>::New();
	points = vtkSmartPointer<vtkPoints>::New();
	table = vtkSmartPointer<vtkTable>::New();
	lookupTable = vtkSmartPointer<vtkLookupTable>::New();
	vertexColors = vtkSmartPointer<vtkIntArray>::New();
	theme = vtkSmartPointer<vtkViewTheme>::New();
	/*createDataForDendogram();
	update();*/

	//MatrixAllocate(&my_data, num_row, num_col,1);
	//MatrixAllocate(&connect_Data_Tree,num_row-1,4,1);
	//
	//MatrixAllocate(&Processed_Coordinate_Data_Tree,(2*num_row)-1, 5,1);
	//MatrixAllocate(&Optimal_Leaf_Nodes,num_row,1,1);

	////////////////Progress.txt read/////////////////////
	//clock_t t2=clock();
	//Read_Meta_Data(num_row,num_col,my_data,"progress.txt");
	//std::cout << "Total time taken for reading progress.txt is: " << ( (clock() - t2)/(float) CLOCKS_PER_SEC) << endl;
	//
	//////////////Optimal_Leaf_Nodes read///////////////////////
	//Read_Meta_Data(num_row,1,Optimal_Leaf_Nodes,"optimal_leaf_order_1D.txt");

	//// Set up sizes. (HEIGHT x WIDTH x DEPTH )
}//needs comments
void Dendrogram::createDataForDendogram()
{
	//std::cout <<"I am entering the CreateDataForDendrogram()"<<std::endl;

    std::cout<< "Enter the number of clusters you want"<<std::endl;
    std::cin>> num_cluster;
   
   
    clock_t t1=clock();
    Determine_File_Chars("progress.txt",&(num_row),&(num_col)); //function returns the number of columns and number of rows in the data
   
    std::cout << "Total time taken for determining the num_row and num_col is: " << ( (clock() - t1)/(float) CLOCKS_PER_SEC) << endl;

    v = new vtkIdType [(num_row+1)*num_row/2]; // an array to store the vertices id- num_row leaf nodes+ 1 node per each level
   
    colour_child=(int *) calloc(2*(this->num_row)-1,sizeof(int));

    MatrixAllocate((&my_data),(num_row),(num_col),1);
    MatrixAllocate(&connect_Data_Tree,num_row-1,4,1);
   
    MatrixAllocate(&Processed_Coordinate_Data_Tree,(2*(this->num_row))-1, 5,1);
    MatrixAllocate(&Optimal_Leaf_Nodes,this->num_row,1,1);
    MatrixAllocate(&distance_data, num_row, 5,1);

    //////////////read Mergers.txt to get the merging distance of the clusters//////////////////////
    clock_t tm=clock();
    Read_Meta_Data(num_row,5,distance_data,"mergers.txt");
    std::cout << "Total time taken for reading mergers.txt is: " << ( (clock() - tm)/(float) CLOCKS_PER_SEC) << endl;
    double y_distance=0;

    ///////////////finding the maximum distance in the mergers.txt/////////////////
    double max=0;
    for(int i=0; i<num_row;i++)
    max= (max>distance_data[i][4])?max:distance_data[i][4];
    std::cout<<"the max is"<<max<<std::endl;

    ////////////////normalizing the data(distance) in mergers.txt//////////////////
    for(int i=0;i<num_row;i++)
        distance_data[i][4]=(distance_data[i][4]/max)*(num_row);


    ////////////reading Progress.txt to get the cluster information /////////////////////
    clock_t t2=clock();
    Read_Meta_Data(num_row,this->num_col,this->my_data,"progress.txt");
    std::cout << "Total time taken for reading progress.txt is: " << ( (clock() - t2)/(float) CLOCKS_PER_SEC) << endl;
   
    ////////////read Optimal_Leaf_Nodes to get the order of the leaf nodes///////////////////////
    Read_Meta_Data(this->num_row,1,this->Optimal_Leaf_Nodes,"optimal_leaf_order_1D.txt");
   
   
   
    ////////////////////// allocating space and sizing the Tree3D matrix ///////////////////////////////
    ////////// Tree3D matrix is an upper diagonal matrix,size= num_row*num_row*6 //////////////////////
    //////// Tree3D[][][0]= id of the node, Tree3D[][][1]=x-coordinate, Tree3D[][][2]=level //////////
    /////// Tree3D[][][3]= z-cooridnate, Tree3D[][][4]= level, Tree3D[][][5]=y-coordinate ///////////
   
    int temp_width=this->num_row;
    this->Tree3D.resize(this->num_row);
    for (int i = 0; i < this->num_row; i++)
        {
            this->Tree3D[i].resize(temp_width);
            for (int j = 0; j < temp_width; j++)
                this->Tree3D[i][j].resize(6);
                temp_width--;
        }
   
    ///////////////////////storing the leaf nodes information in Tree3D matrix////////////////////////////////
    //////the first leaf node gets x-coordinate=0, each node is spaced at 1 unit distance horizontally////////

    for(int i=0;i<this->num_row;i++)
    {

        this->Tree3D[0][i][0]=Optimal_Leaf_Nodes[i][0];
        this->Tree3D[0][i][1]=i;
        this->Tree3D[0][i][2]=0;
        this->Tree3D[0][i][3]=0;
        this->Tree3D[0][i][4]=0;
        this->Tree3D[0][i][5]=0;

    }

    ///////////creating the vtkTable from the data in my_data//////////////
   
 
         
    for ( unsigned int i = 0; i < this->num_col; i++ )
    {
        vtkSmartPointer<vtkVariantArray> col =
        vtkSmartPointer<vtkVariantArray>::New();
            for ( unsigned int j=0;j<this->num_row;j++)
                {
                col->InsertNextValue ( vtkVariant ( this->my_data[j][i] ) );
                }
        this->table->AddColumn ( col );
    }

    //////////////Creating Tree3D Matrix using the table//////////////////
    double record[2];
    clock_t t3=clock();
    int k; /// to keep track of the next row being created in the Tree3D Matrix

    for(vtkIdType c = 2; c < this->num_col; c++)//iterating the columns in my_data
    {       
        int i;
       
            for( i=0;i<this->num_row-c+2;i++)// scanning the row already created in Tree3D matrix
                {
                k=0;               
                int temp=this->Tree3D[c-2][i][0];
                vtkVariant v1=this->table->GetValue(temp,c); //get the ID of the element which is present in the tree matrix
                double temp1= v1.ToDouble();
       
                count=0; // Temporary count of the the Original Object IDs while scanning a particular column of vtkTable
     
                int j=0; // While loop control
                    while(j<this->num_row-c+2)
                        {
                        vtkVariant v = this->table->GetValue( this->Tree3D[c-2][j][0],c);
                        double chk=v.ToDouble();                       
                            if(temp1==chk)
                            {   
                            record[count]=this->Tree3D[c-2][j][0];// storing the two id's that are merging in record
                            count++;
                       
                            if(count==2)// if count==2 implies both have the same id's implies they are merging together
                            break;
                            }
               
                        j++;
                        }
                if(count==2)
                        break;// since we break twice the value of i where it merges is retained
             }
    for(int l=0;l<this->num_row-c+1;l++)
            {
            if(l<i)//if l<i implies its not merging so push the same information into the next level
                {
                    this->Tree3D[c-1][l][0]=Tree3D[c-2][l][0];
                    this->Tree3D[c-1][l][1]=Tree3D[c-2][l][1];
                    this->Tree3D[c-1][l][2]=Tree3D[c-2][l][2];
                    this->Tree3D[c-1][l][3]=Tree3D[c-2][l][3];
                    this->Tree3D[c-1][l][4]=c-1;// incrementing the level by one
                    this->Tree3D[c-1][l][5]=Tree3D[c-2][l][5];


                }
            else if(l==i)
                {
                y_distance=distance_data[c-2][4];
                this->Tree3D[c-1][l][0]=(record[0]>record[1])?record[0]:record[1];// assigning the maximum of the two ID's to the ID of the element thats forming after merging 
                x1=this->Tree3D[c-2][l][1];
                x2=this->Tree3D[c-2][l][2];
                x3=this->Tree3D[c-2][l][3];

                y1=this->Tree3D[c-2][l+1][1];
                y2=this->Tree3D[c-2][l+1][2];
                y3=this->Tree3D[c-2][l+1][3];
               
                ////////////////assigning the coordinates to the node(cluster)/////////////////////////////
                ////////////////thats getting created after merging two nodes(clusters)///////////////////

                this->Tree3D[c-1][l][1]=(x1+y1)/2; ///midpoint of the x-coordinates
                this->Tree3D[c-1][l][2]=c-1; //level information
                this->Tree3D[c-1][l][3]=(x3+y3)/2; //midpoint of the y-coordinates
                this->Tree3D[c-1][l][4]=c-1; //level information
                this->Tree3D[c-1][l][5]=y_distance; // y-coordinate, obtained by getting the merging distance

                /////////////connect_Data_Tree is a matrix that stores the ids of the nodes merging/////////


                this->connect_Data_Tree[c-2][0]=record[0];
                this->connect_Data_Tree[c-2][1]=record[1];
                this->connect_Data_Tree[c-2][2]=x2;// storing the level of the element which is merging
                this->connect_Data_Tree[c-2][3]=y2;

               

                }
            else
                {
                    this->Tree3D[c-1][l][0]=Tree3D[c-2][l+1][0];
                    this->Tree3D[c-1][l][1]=Tree3D[c-2][l+1][1];
                    this->Tree3D[c-1][l][2]=Tree3D[c-2][l+1][2];
                    this->Tree3D[c-1][l][3]=Tree3D[c-2][l+1][3];
                    this->Tree3D[c-1][l][4]=c-1;
                    this->Tree3D[c-1][l][5]=Tree3D[c-2][l+1][5];
                }   

            }

    }
std::cout << "Total time taken for data generation is: " << ( (clock() - t3)/(float) CLOCKS_PER_SEC) << endl;

 

std::cout << "I have reached my check point"<<std::endl;
  int temp=num_row; // for temporary storage of the num_row
 

  ///////////////////////////////process the points////////////////////////
  //////////Processed_Coordinate_Tree is matrix that reads all the unique points in Tree3D//////
  //////////Processed_Coordinate_Tree size = ((2*num_row)-1)*5/////////////
  //////////Processed_Coordinate_Tree[][0]=id of the node/point///////////
  //////////Processed_Coordinate_Tree[][1]=x-coordinate////////////
  //////////Processed_Coordinate_Tree[][2]=y-coordiante///////////
  //////////Processed_Coordinate_Tree[][3]=z-coordiante//////////
  //////////Processed_Coordinate_Tree[][4]=level////////////////


    double level=0.0;  // To keep track of the tree level in the dendrogram
    int temp1=this->num_row; // temporary storage of the num_row
   

   
    clock_t t4=clock();
    for(int i=0;i<((2*this->num_row)-1);i++)
    {

    ///////////////////pushing the leaf nodes information into Processed_Coordinate_Tree////////
        if(i<this->num_row)
        {
            this->Processed_Coordinate_Data_Tree[i][0]= Tree3D[0][i][0];
            this->Processed_Coordinate_Data_Tree[i][1]=Tree3D[0][i][1];
            this->Processed_Coordinate_Data_Tree[i][2]=Tree3D[0][i][5];
            this->Processed_Coordinate_Data_Tree[i][3]=Tree3D[0][i][3];
            this->Processed_Coordinate_Data_Tree[i][4]=Tree3D[0][i][4];

        }
       
    /////////////checking for all unique nodes/points in Tree3D matrix/////////////////////
    /////////////and pushing them into Processed_Coordinate_Tree//////////////////////////
        else
        {

            level++;
           
                for(int k=0;k<temp1;k++)
                {
                    if(this->Tree3D[level][k][2] == level)//there is one unique node per level
                        {
                       
                        this->Processed_Coordinate_Data_Tree[i][0]=this->Tree3D[level][k][0];
                        this->Processed_Coordinate_Data_Tree[i][1]=this->Tree3D[level][k][1];
                        this->Processed_Coordinate_Data_Tree[i][2]=this->Tree3D[level][k][5];
                        this->Processed_Coordinate_Data_Tree[i][3]=this->Tree3D[level][k][3];
                        this->Processed_Coordinate_Data_Tree[i][4]=this->Tree3D[level][k][4];
                        break;
                        }
                   
                   
                }
               
                temp1=temp1-1;
               
           
        }   
    }
   
   

std::cout << "Total time taken for processing the points is: " << ( (clock() - t4)/(float) CLOCKS_PER_SEC) << endl;

//////////////////////////printing Processed_Coordinate_Tree matrix////////////////////////////////
//for(int i=0; i<((2*this->num_row)-1); i++)
//std:: cout << "the processed coordinates are " << Processed_Coordinate_Data_Tree[i][0] <<"    "<< Processed_Coordinate_Data_Tree[i][1] <<"    "<< Processed_Coordinate_Data_Tree[i][2] <<"    "<<Processed_Coordinate_Data_Tree[i][3] <<"     "<< Processed_Coordinate_Data_Tree[i][4] <<std::endl;

//std::cout << "I am comming out of the createDataForDendrogram"<<std::endl;
}



void Dendrogram::update()
{
    //num_cluster=3;
    clust_level= num_row - num_cluster;

	///////////////////////Adding vertices to the graph//////////////////

    for(int i=0; i<((2*this->num_row)-1);i++)
    {
     
       
      v[i]=graph_Layout->AddVertex();
       
        this->points->InsertNextPoint(this->Processed_Coordinate_Data_Tree[i][1],this->Processed_Coordinate_Data_Tree[i][2],this->Processed_Coordinate_Data_Tree[i][3]);// assigning the coordintaes to the vertex after reading it from Processed_Coordinate_Data_Tree
       // std:: cout<<" the point inserted is " << Processed_Coordinate_Data_Tree[i][1]<< "  "<<Processed_Coordinate_Data_Tree[i][2]<< "    " << Processed_Coordinate_Data_Tree[i][3]<< "    " << Processed_Coordinate_Data_Tree[i][4] <<"    " << Processed_Coordinate_Data_Tree[i][5]<< std:: endl;
    }
   
    this->graph_Layout->SetPoints(this->points);
   
   
    ///////////////coloring/////////////////////
    /*vertexColors = vtkSmartPointer<vtkIntArray>::New();*/
    vertexColors->SetNumberOfComponents(1);
    vertexColors->SetName("Color");
 
    /*lookupTable = vtkSmartPointer<vtkLookupTable>::New();*/
    lookupTable->SetNumberOfTableValues(2*(this->num_row)-1);
    for(int i=0; i<2*(this->num_row)-1;i++)

    {
    lookupTable->SetTableValue(i, 0, 0, 1.0); // color the vertices- blue
    }
   
    lookupTable->Build();
   
   
   vtkSmartPointer<vtkFloatArray> scales = vtkSmartPointer<vtkFloatArray>::New(); /// scales for vertex size
    scales->SetNumberOfComponents(1);
     scales->SetName("Scales");

    for(int j=0;j<2*(this->num_row)-1;j++)
    {
   
    vertexColors->InsertNextValue(j);
    scales->InsertNextValue(1.3);

    }
   

//////////////////Generating original ID's//////////////////////////

    this->CharLabel.resize(this->num_row);
    for(int i=0;i<this->num_row;i++)
    {
        this->CharLabel[i].resize(this->num_row+1);
    }
    for(int i=0;i<this->num_row;i++)
    {
        for( int j=0;j<this->num_row+1;j++)
        {
            this->CharLabel[i][j]=-1;
        }
    }

    for(int i=0;i<this->num_row;i++)
        this->CharLabel[0][i]=this->Optimal_Leaf_Nodes[i][0];

    //////////////////For reading the labels////////////

    clock_t t6=clock();
   
    for(vtkIdType c = 2; c < this->num_col; c++)
    {
        double temp=this->connect_Data_Tree[c-2][0];
        vtkVariant v1=this->table->GetValue(temp,c); //get the ID of the element which is present in the tree matrix
        double temp1= v1.ToDouble();
        int m=0;
        for(int l=0;l<(this->table->GetNumberOfRows());l++)
            {
                vtkVariant v2=this->table->GetValue(l,c);
                if( v2==temp1)
                    {
                        this->CharLabel[c-1][m]=(this->table->GetValue(l,0)).ToDouble();
                        m++;
                    }
               
               
            }
    }

std::cout << "Total time taken for processing the points is: " << ( (clock() - t6)/(float) CLOCKS_PER_SEC) << endl;


///////////////////////////setting a default colour///////////////////
   
    //for(int i=0; i<num_row;i++)
    //{
    //    double temp= CharLabel[0][i];
    //    vtkVariant v1=table->GetValue(temp,num_row+1-num_cluster); //get the ID of the element which is present in the tree matrix
    //    double id= v1.ToDouble();
    //    lookupTable->SetTableValue(i, 0, (id/(num_cluster-1)), 0);
    //}

////////////////////////////////////////////////////////////////////////////
     
       this->graph_Layout->GetVertexData()->AddArray(vertexColors);
    this->graphLayoutView = vtkSmartPointer<vtkGraphLayoutView>::New();
    this->graphLayoutView->AddRepresentationFromInput(graph_Layout);
    this->graphLayoutView->SetLayoutStrategy("Pass Through");
    
    //this->graph_Layout->GetVertexData()->AddArray(scales);
    this->graph_Layout->GetVertexData()->AddArray(scales);
    this->graphLayoutView->ScaledGlyphsOn();
    this->graphLayoutView->SetScalingArrayName("Scales");
     vtkRenderedGraphRepresentation::SafeDownCast(this->graphLayoutView->GetRepresentation()) ->SetGlyphType(vtkGraphToGlyphs::CIRCLE);


    this->graphLayoutView->ResetCamera();
    this->graphLayoutView->SetVertexColorArrayName("Color");
    this->graphLayoutView->ColorVerticesOn();

    ///////////////////////setting the shape of vertices//////////////////////////////////

    theme->SetPointLookupTable(lookupTable);
    theme->SetBackgroundColor(1,1,1);   /// to set the background colour of the dendrogram window
    this->graphLayoutView->ApplyViewTheme(theme);
   
    std::cout<<"the no of vertices"<<graph_Layout->GetNumberOfVertices()<<std::endl;



   

    balloonRep1 = vtkSmartPointer<vtkBalloonRepresentation>::New();

    balloonWidget1 = vtkSmartPointer<vtkBalloonWidget>::New();
    balloonRep = vtkSmartPointer<vtkBalloonRepresentation>::New();
    balloonWidget = vtkSmartPointer<vtkBalloonWidget>::New();

    ////////////////////////////Labelling the Leaf Nodes////////////////////////////////////////


    //vtkSmartPointer<vtkTextActor> textActor =
 //   vtkSmartPointer<vtkTextActor>::New();
    //textActor->GetTextProperty()->SetFontSize ( 24 );
    //
    //textActor->SetDisplayPosition(a[0], a[1]);
 //   //textActor->SetPosition2 ( 1000, 500 );
    //this->graphLayoutView->GetRenderer()->AddActor2D ( textActor );
 //   textActor->SetInput ( "1" );
 //   textActor->GetTextProperty()->SetColor ( 1.0,0.0,0.0 );
//
//vtkSmartPointer<vtkAxisActor2D> axisActor2D1 =
//    vtkSmartPointer<vtkAxisActor2D>::New();
//axisActor2D1->SetPoint1(0,50);
//axisActor2D1->SetPoint2(50,50);
//axisActor2D1->AxisVisibilityOn();
//this->graphLayoutView->GetRenderer()->AddActor2D(axisActor2D1);


///////////////////////////Displaying CharLabel Matrix///////////////////////////////
//for(int i=0;i<num_row;i++)
//    {
//        std::cout << std::endl;
//        for(int j=0;j<num_row+1;j++)
//            std::cout <<CharLabel[i][j]<<"\t";
//    }
/////////////////connecting points/////////////////////////////////////////////
    double p1[3];
    double p2[3];
    double p3[3];
    double p4[3];
   
    clock_t t5=clock();

   

    for(int i=0;i<(num_row-1);i++)
    {
        double temp1=this->connect_Data_Tree[i][0];
        double temp2=this->connect_Data_Tree[i][1];
        double temp3=this->connect_Data_Tree[i][2];
        double temp4=this->connect_Data_Tree[i][3];

        for(int j=0;j<(2*(this->num_row))-1;j++)
        {
            if(this->Processed_Coordinate_Data_Tree[j][0]==temp1 && this->Processed_Coordinate_Data_Tree[j][4]==temp3)
                {
                p1[0]=this->Processed_Coordinate_Data_Tree[j][1];
                p1[1]=this->Processed_Coordinate_Data_Tree[j][2];
                p1[2]=this->Processed_Coordinate_Data_Tree[j][3];
                }   
            if(this->Processed_Coordinate_Data_Tree[j][0]==temp2 && this->Processed_Coordinate_Data_Tree[j][4]==temp4)
                {
                p2[0]=this->Processed_Coordinate_Data_Tree[j][1];
                p2[1]=this->Processed_Coordinate_Data_Tree[j][2];
                p2[2]=this->Processed_Coordinate_Data_Tree[j][3];
                }   
                           
        }
        p3[0]=p1[0];
        p3[1]=Processed_Coordinate_Data_Tree[(num_row+i)][2];
        p3[2]=p1[2];

        p4[0]=p2[0];
        p4[1]=Processed_Coordinate_Data_Tree[(num_row+i)][2];
        p4[2]=p2[2];
    vtkSmartPointer<vtkLineSource> lineSource1 = vtkSmartPointer<vtkLineSource>::New();
    lineSource1->SetPoint1(p1);
    lineSource1->SetPoint2(p3);
   
    vtkSmartPointer<vtkPolyDataMapper> mapper1 =
    vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper1->SetInputConnection(lineSource1->GetOutputPort());

    vtkSmartPointer<vtkActor> actor1 =
    vtkSmartPointer<vtkActor>::New();
    actor1->GetProperty()->SetColor(0.5,0.7,0); //set the colours of the line of the dendrogram
    actor1->GetProperty()->SetLineWidth(1.5);
    actor1->SetMapper(mapper1);
    this->graphLayoutView->GetRenderer()->AddActor(actor1);


    vtkSmartPointer<vtkLineSource> lineSource2 = vtkSmartPointer<vtkLineSource>::New();
    lineSource2->SetPoint1(p2);
    lineSource2->SetPoint2(p4);
   
    vtkSmartPointer<vtkPolyDataMapper> mapper2 =
    vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper2->SetInputConnection(lineSource2->GetOutputPort());

    vtkSmartPointer<vtkActor> actor2 =
    vtkSmartPointer<vtkActor>::New();
    actor2->SetMapper(mapper2);
    actor2->GetProperty()->SetColor(0.5,0.7,0);
    actor2->GetProperty()->SetLineWidth(1.5);
    this->graphLayoutView->GetRenderer()->AddActor(actor2);

    vtkSmartPointer<vtkLineSource> lineSource3 = vtkSmartPointer<vtkLineSource>::New();
    lineSource3->SetPoint1(p3);
    lineSource3->SetPoint2(p4);
   
    vtkSmartPointer<vtkPolyDataMapper> mapper3 =
    vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper3->SetInputConnection(lineSource3->GetOutputPort());

    vtkSmartPointer<vtkActor> actor3 =
    vtkSmartPointer<vtkActor>::New();
    actor3->GetProperty()->SetColor(0.5,0.7,0);
    actor3->GetProperty()->SetLineWidth(1.5);
    actor3->SetMapper(mapper3);
     ///////////////////////////////////

     balloonRep1->SetBalloonLayoutToImageRight();
 
   
    balloonWidget1->SetInteractor(graphLayoutView->GetInteractor());
    balloonWidget1->SetRepresentation(balloonRep1);
    char str[30];
    int buffer;
    buffer=i+1;
    int k=0;
    int num_leaf=0;
    while(CharLabel[buffer][k]!= -1)
    {
        num_leaf=num_leaf+1;
        k=k+1;
    }

    char *str1;
    str1="Merge Level : ";
    //itoa (buffer,str,10);
    std::string dummy;
    dummy = ftk::NumToString(buffer);
    strcpy (str,dummy.c_str());
    char str2[80];
    strcpy (str2,str1);
    strcat (str2,str);
    char *str3;
    char str4[30];
    str3="\nNo of Leaf Nodes : ";
    dummy.clear();
    dummy = ftk::NumToString(num_leaf);
    strcpy (str4,dummy.c_str());
    //itoa (num_leaf,str4,10);
    char str5[80];
    strcpy (str5,str3);
    strcat (str5,str4);
    strcat (str2,str5);
    puts (str2);
    balloonWidget1->SetTimerDuration(0.001);
    balloonWidget1->AddBalloon(actor3,str2, NULL);
   
    this->graphLayoutView->GetRenderer()->AddActor(actor3);
    }
    std::cout << "Total time taken for placing the points is: " << ( (clock() - t5)/(float) CLOCKS_PER_SEC) << endl;

    ///////////////////////cutting the dendrogram by drawing a line//////////////////////
    if(num_cluster >1&& num_cluster < this->num_row)
    {
    double emp1=Processed_Coordinate_Data_Tree[clust_level+num_row][2];
    /*std::cout << "Now i am going to draw the line for cutting the dendrogrm"<<std::endl;*/
    double emp2=Processed_Coordinate_Data_Tree[clust_level+num_row-1][2];
    double emp3=(emp1+emp2)/2;
    double p5[3];
    double p6[3];
    p5[0]=-1;
    p5[1]=emp3;
    p5[2]=0;
    p6[0]=num_row+1;
    p6[1]=emp3;
    p6[2]=0;
       
    vtkSmartPointer<vtkLineSource> lineSource4 = vtkSmartPointer<vtkLineSource>::New();
    lineSource4->SetPoint1(p5);
    lineSource4->SetPoint2(p6);
   
    vtkSmartPointer<vtkPolyDataMapper> mapper4 =
    vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper4->SetInputConnection(lineSource4->GetOutputPort());

    vtkSmartPointer<vtkActor> actor4 =
    vtkSmartPointer<vtkActor>::New();
    actor4->GetProperty()->SetColor(1,0,0);
    actor4->GetProperty()->SetLineWidth(2);
    actor4->SetMapper(mapper4);
   
   
    balloonRep->SetBalloonLayoutToImageRight();

   
 
    balloonWidget->SetInteractor(graphLayoutView->GetInteractor());
    balloonWidget->SetRepresentation(balloonRep);
    balloonWidget->SetTimerDuration(0.001);
    balloonWidget->AddBalloon(actor4, "This is a cutting line",NULL);
   
 //   //renderer->AddActor(sphereActor);
    this->graphLayoutView->GetRenderer()->AddActor(actor4);
    //graphLayoutView->Render();
    graphLayoutView->GetRenderer()->GradientBackgroundOff();

    this->graphLayoutView->Render();
    balloonWidget1->EnabledOn();
    balloonWidget->EnabledOn();
   
    }
    else
        exit(0);
    /*std::cout << "Now i have drawn the line for cutting the dendrogrm"<<std::endl;*/


    ////////////////Interactor ////////////////////
    this->selectionCallback = vtkSmartPointer<vtkCallbackCommand>::New();
    this->selectionCallback->SetClientData(this);
    this->selectionCallback->SetCallback (SelectionCallbackFunction);
    this->graphLayoutView->GetRepresentation()->GetAnnotationLink()->AddObserver("AnnotationChangedEvent", this->selectionCallback);
    /////////////////////////////////////////////////////////
    connect(this,SIGNAL(selection_Changed(int)),this,SLOT(GetSelectedLevel(int)));
   
    /*widget.SetRenderWindow(graphLayoutView->GetRenderWindow());
    widget.show();*/
    this->graphLayoutView->GetInteractor()->Start();
    // app.exec();
 
}
Dendrogram::~Dendrogram()
{
}
int Dendrogram::Determine_File_Chars(char *root, int *num_data, int *numfeats) {
int      i, num_features = 0, num_lines = 0, num_chars = 0;
char     ch = 0;
VECTOR   cum_sum;
FILE     *fpin;
//printf("I have entered in the function Determine_File_Chars");
  fpin = FDeclare(root, "", 'r');
  while (1) {
    ch = 1;
    fscanf(fpin, "%c", &ch);
    if (feof(fpin)) break;
    if (ch == '\n') num_lines += 1;
  }
  rewind(fpin);
 /* printf("!!!There are %2d lines\n", num_lines);*/
  ch = 1;
  while (ch!='\n') {
    fscanf(fpin, "%c", &ch);
    num_chars += 1;
  }
 /* printf("the number of characters are:%d\n",num_chars);*/
  rewind(fpin);
  VectorAllocate(&cum_sum, num_chars);
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
FILE* Dendrogram::FDeclare(char *root, char *extension, char key) {
char string1[80];
FILE *fpot;

  strcpy(string1, root);
  if (strcmp(extension,"")!=0) strcat(string1, ".");
  strcat(string1, extension);
  if ((key != 'w')&&(key!='r')) {
    printf("either read or write, wrong key %c \n", key);//getche();
    exit(1);
  }
  if (key == 'r') {
    //if (( (FILE *) fpot = fopen(string1, "r")) == NULL) {
    //changed for DOS?UNIX compatibility
    if ((fpot = fopen(string1, "r")) == NULL) {
		std::cout << std::endl;
		std::cout << "input file "<<string1<< "does not exist" <<std::endl;
      exit(1); //getche();
    }
  }
  if (key == 'w') {
    if ((fpot = fopen(string1, "w")) == NULL) {
      printf("\n");
      printf("output file %s does not exist\n", string1);
      exit(1);//getche();
    }
  }
  return fpot;
}

void Dendrogram::VectorAllocate(VECTOR *vector, int nCols) {
  if ((*vector=(VECTOR) calloc(nCols, sizeof(double))) == NULL) {
    printf("Not enough memmry\n");
    exit(1);
  }
}
void Dendrogram::Read_Meta_Data(int num_data, int num_feats, MATRIX my_data, char *root) {
int      i, j;
double  templf;
FILE     *fpin;

std::cout << "I am here";
  fpin = FDeclare(root,".",'r');
    for (i = 0; i < num_data; i++) {
    for (j = 0; j < num_feats; j++) {
      fscanf(fpin, "%lf", &templf);
      my_data[i][j] = templf;
    }
  }
  fclose(fpin);
}
void Dendrogram::MatrixAllocate(MATRIX *pmatrix, int nRows, int nCols,int mode) {
	if ((*pmatrix=(MATRIX) calloc(nRows, sizeof(std::string))) == NULL) {
    printf("Not enough memmry\n");
    exit(1);
	}
	if(mode==1)
   AllocateCols(*pmatrix, nRows, nCols);
}
void Dendrogram::AllocateCols(PFLOAT matrix[], int nRows, int nCols) {
int i;
  for (i = 0; i < nRows; i++)
   VectorAllocate(&matrix[i], nCols);
}

void Dendrogram::SelectionCallbackFunction(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData )
{
 
    Dendrogram* Dendro = (Dendrogram*)clientData;
    vtkAnnotationLink* annotationLink =
    static_cast<vtkAnnotationLink*>(caller);
    vtkSmartPointer<vtkDataObject> view =
    vtkSmartPointer<vtkDataObject>::New();
    vtkSelection* selection = annotationLink->GetCurrentSelection();
 
    vtkSelectionNode* vertices;
    vtkSelectionNode* edges;
    if(selection->GetNode(0)->GetFieldType() == vtkSelectionNode::VERTEX)
        {
        vertices = selection->GetNode(0);
        }
    else if(selection->GetNode(0)->GetFieldType() == vtkSelectionNode::EDGE)
        {
        edges = selection->GetNode(0);
        }
 
    if(selection->GetNode(1)->GetFieldType() == vtkSelectionNode::VERTEX)
        {
        vertices = selection->GetNode(1);
        }
    else if(selection->GetNode(1)->GetFieldType() == vtkSelectionNode::EDGE)
        {
        edges = selection->GetNode(1);
        }
 
    vtkIdTypeArray* vertexList = vtkIdTypeArray::SafeDownCast(vertices->GetSelectionList());
    //std::cout << "There are " << vertexList->GetNumberOfTuples() << " vertices selected." << std::endl;
 
    if(vertexList->GetNumberOfTuples() > 0)
        {
        //std::cout << "Vertex Ids: ";
        }

    if(vertexList->GetNumberOfTuples()==0)
        {
            for(int n=0;n<2*(Dendro->num_row)-1;n++)
            {
            Dendro->lookupTable->SetTableValue(n, 0.0, 0, 1.0);
            }
            emit Dendro->selection_Changed(0);
        }
    else
        {   
                int p;
            for(vtkIdType i = 0; i < vertexList->GetNumberOfTuples(); i++)
            {
            //std::cout << vertexList->GetValue(i) << " ";
            vtkIdType v=vertexList->GetValue(i);
            p=v;
                if(v>(Dendro->num_row-1))
                    {
     
                       
                           
                            int j=0;
                            double x1=Dendro->num_row;
                            double x2=0;
                            double temp1=0;
                            //double temp2=num_row;
                            double y1=Dendro->Processed_Coordinate_Data_Tree[v][2];
                            for(int i=0;i<2*(Dendro->num_row)-1;i++)
                                Dendro->colour_child[i]=0;

                            //int count=0;
                            //int flag=0;
                            int k;
                            while((Dendro->CharLabel[p-(Dendro->num_row-1)][j])!=-1)
                                {
                               
                                //std::cout << "The children are "<<(Dendro->CharLabel[p-(Dendro->num_row-1)][j])<<"\t";
                               
                                for(k=0;k<Dendro->num_row;k++)
                                    {
                                   
                                    double rep= (Dendro->Tree3D[0][k][0]);
                                    if((rep)==(Dendro->CharLabel[p-(Dendro->num_row-1)][j]))
                                        {
                                        Dendro->colour_child[k]=1;
                                        temp1=Dendro->Tree3D[0][k][1];
                                        if(x1>=temp1)
                                            x1=temp1;
                                        if(x2<=temp1)
                                            x2=temp1;
                                        //count=count+1;
                                        //flag=1;
                                       
                                        }
                                    /*else
                                        flag=0;*/
                                   
                                   
                                    }
                               
                                       
                                    /*
                                    if(flag==1)
                                        x2=Dendro->Tree3D[1][k][0];*/
                               
           
                                j++;

                               }
                            //x2=Dendro->Tree3D[0][k][1];
                            //std::cout<<"x1 "<<x1<<"x2 "<<x2<<" y1 "<<y1<<std::endl;
                            //std::cout<<"the slected ids are" << std:: endl;
                            for(int m=0;m<2*(Dendro->num_row)-1;m++)
                            {
                                if((Dendro->Processed_Coordinate_Data_Tree[m][1])>=x1 && (Dendro->Processed_Coordinate_Data_Tree[m][1])<=x2 && (Dendro->Processed_Coordinate_Data_Tree[m][2])<=y1)
                                {
                                    //std::cout<< m << '\t';
                                    Dendro->colour_child[m]=1;
                                }
                            }
                            /*for(int  m=0;m<2*(Dendro->num_row)-1;m++)
                                std::cout<< Dendro->colour_child[m] << '\t';*/

                            for(int n=0;n<2*(Dendro->num_row)-1;n++)
                            {
                                if(Dendro->colour_child[n]==1)
                                Dendro->lookupTable->SetTableValue(n, 1.0, 0.0, 1.0);
                                else
                                Dendro->lookupTable->SetTableValue(n, 0.0, 0, 1.0);
                            }
                            emit Dendro->selection_Changed(p);

                           
                       
                        }
                 else
                 {
                    for(int n=0;n<2*(Dendro->num_row)-1;n++)
                    Dendro->lookupTable->SetTableValue(n, 0.0, 0, 1.0);
                    emit Dendro->selection_Changed(p);

                }

            }
        }
 
  vtkIdTypeArray* edgeList = vtkIdTypeArray::SafeDownCast(edges->GetSelectionList());
 
}
vtkSmartPointer<vtkGraphLayoutView> Dendrogram::GetGraphLayoutView()
{

	return graphLayoutView;
}



void Dendrogram::GetSelectedLevel(int level)
{

	Level_Detected=&level;
	int p=*Level_Detected;
	std::cout << "The value of the selected level:"<<p<<std::endl;
	//get set of leaf node ids 
	std::set<long int> IDs;
	if(p!=0)
	{
	
		if(p < this->num_row)
		{	
			for(int i=0;i<1;i++)
			{
			IDs.insert(p);
			}
		}
		else
		{
		/*for(int i=0;i<this->num_row;i++)*/
			int i=0;
			while(this->CharLabel[p-(this->num_row-1)][i]!=-1)
			{
			IDs.insert(this->CharLabel[p-(this->num_row-1)][i]+1);
			i++;
			}
		}
	}
	else
		IDs.insert(0);
	
	this->setSelectedIds(IDs);  

}



void Dendrogram::setSelectedIds(std::set<long int> IDs)
{
	
	if(IDs.size()>1)
	this->Selection->select(IDs);
	else
		this->Selection->clear();
}
void Dendrogram::GetSelecectedIDs()
{
	std::set<long int> selection = this->Selection->getSelections();
	//highlight selections in dendrogram view here
	std::set<long int>::iterator it;

	int i=0;
	for(it=selection.begin();it!=selection.end();it++)
	{
		//lookupTable->SetTableValue((*it)-1,1,0,1);
		std::cout <<"The selected ID is :"<<(*it)-1<<std::endl;
		
		
	}
	
	
	this->graphLayoutView->ColorVerticesOn(); 

	vtkSmartPointer<vtkViewTheme> theme = vtkSmartPointer<vtkViewTheme>::New() ;
	theme->SetPointLookupTable(lookupTable);
	theme->SetBackgroundColor(1,1,1);   /// to set the background colour of the dendrogram window
	this->graphLayoutView->ApplyViewTheme(theme);
	
	
	lookupTable->Build();
	//std::cout << *it<<std::endl;


}

void Dendrogram::setModels(vtkSmartPointer<vtkTable> tbl, ObjectSelection * sels)
{

	std::cout << "I am here"<<std::endl;
	tbl->Dump(3);
	Cluster = new ClusClus(tbl); // Clustering is performed by this class on the feature table
	this->table=Cluster->GetProgressTable();
	//this->table = tbl;
	std::cout <<"The table containing the progress.txt info is"<<std::endl;
	table->Dump(3);					// The table contains progress.txt in vtkTable format
	std::cout << "I am inside setModels"<<std::endl;
		this->createDataForDendogram();
	if(!sels)
		this->Selection = new ObjectSelection();
	else
		this->Selection = sels;
	connect(Selection, SIGNAL(changed()), this, SLOT(GetSelecectedIDs()));

	this->update();
}





