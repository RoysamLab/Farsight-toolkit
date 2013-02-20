#if defined(_MSC_VER)
#pragma warning(disable : 4996)
#define strcasecmp _stricmp
#endif

#include "NESqliteFactory.h"
#ifdef _MSC_VER
    #include <direct.h>
    #define GetCurDir _getcwd
#else
    #include <unistd.h>
    #define GetCurDir getcwd
#endif

namespace ftk
{

//void sqliteExecuteQuery(sqlite3 *db)
//{ 
//	sqlite3_stmt *ppStmt;
//	string queryStr,queryStr1,queryStr2;
//	char  *name, *img_name;
//	int   exeStatus,img_id;
//	const char *tail;
//	char  *zErrMsg ;
//	name  = "PleaseChalja8";
//
//	int  nrow;       /* Number of result rows written here */
//	int  ncol;       /* Number of result columns written here */
//	char **result;
//	
//	queryStr = "insert into IMAGE (IMG_NAME, IMG_LOCATION)"
//		       "values (:IMG_NAME, :IMG_LOCATION)";
//	
//	sqlite3_prepare(db, queryStr.c_str(), queryStr.size(), &ppStmt, &tail);
//	sqlite3_bind_text( ppStmt,sqlite3_bind_parameter_index(ppStmt,":IMG_NAME"),name,strlen(name),SQLITE_TRANSIENT );
//	sqlite3_bind_text( ppStmt,sqlite3_bind_parameter_index(ppStmt,":IMG_LOCATION"),name,strlen(name),SQLITE_TRANSIENT );
//	exeStatus = sqlite3_step(ppStmt);
//	sqlite3_finalize(ppStmt);
//
//	/*2)Select IMG_ID from the IMAGE */
//	queryStr = "select IMG_ID,IMG_NAME from IMAGE where IMG_NAME = :IMG_NAME";         
//	exeStatus = sqlite3_prepare(db, queryStr.c_str(), queryStr.size(), &ppStmt, &tail);
//	if( exeStatus != SQLITE_OK )
//		fprintf(stderr, "Sqlite3_prepare error: %s\n", sqlite3_errmsg(db));
//	else
//	{	sqlite3_bind_text(ppStmt,sqlite3_bind_parameter_index(ppStmt,":IMG_NAME"),
//		name,strlen(name),SQLITE_TRANSIENT );
//		exeStatus = sqlite3_step(ppStmt);
//		while(exeStatus == SQLITE_ROW) 
//		{ img_id   = sqlite3_column_int(ppStmt, 0);
//		  img_name = (char*)sqlite3_column_text(ppStmt,1);
//		  exeStatus = sqlite3_step(ppStmt);
//		  fprintf(stderr,"img_id=%i \n", img_id);
//		}
//		sqlite3_finalize(ppStmt);
//	}
//
//	sqlite3_get_table(db, "select * from IMAGE_TEST", &result, &nrow, &ncol, &zErrMsg);
//	fprintf(stderr, "Number of row:%i  and Columns:%i\n \n",nrow,ncol);
//	/*
//	int sqlite3_get_table(
//		sqlite3 *db,           An open database 
//		const char *sql,      SQL to be evaluated 
//		char ***result,     Results of the query 
//		int *Row,            Number of result rows written here 
//		int *Column,         Number of result columns written here 
//		char **Errmsg        Error msg written here 
//		)
//	*/  
//
//	/*3)Inserting data in IMAGE_TEST */	
//	/*queryStr = "insert into IMAGE_TEST" 
//		"(IMG_ID, CELL_ID, VOLUME,INTEGRATED_INTENSITY,ECCENTRICITY,ELONGATION,"
//		"ORIENTATION,BOUNDING_BOX_VOLUME,SUM,MEAN,MEDIAN,MINIMUM,MAXIMUM,"
//		"SIGMA,VARIANCE,SURFACE_GRADIENT,INTERIOR_GRADIENT,SURFACE_INTENSITY,INTERIOR_INTENSITY,"
//		"INTENSITY_RATIO,RADIUS_VARIATION,SURFACE_AREA,SHAPE,SHARED_BOUNDARY,SKEW,ENERGY,ENTROPY,T_ENERGY,"
//		"T_ENTROPY,INVERSE_DIFF_MOMENT,INERTIA,CLUSTER_SHADE,CLUSTER_PROMINENCE)"
//		"values"
//		"(:IMG_ID,:CELL_ID,:VOLUME,:INTEGRATED_INTENSITY,:ECCENTRICITY,:ELONGATION,"
//		":ORIENTATION,:BOUNDING_BOX_VOLUME,:SUM,:MEAN,:MEDIAN,:MINIMUM,:MAXIMUM,"
//		":SIGMA,:VARIANCE,:SURFACE_GRADIENT,:INTERIOR_GRADIENT,:SURFACE_INTENSITY,:INTERIOR_INTENSITY,"
//		":INTENSITY_RATIO,:RADIUS_VARIATION,:SURFACE_AREA,:SHAPE,:SHARED_BOUNDARY,:SKEW,:ENERGY,:ENTROPY,:T_ENERGY,"
//		":T_ENTROPY,:INVERSE_DIFF_MOMENT,:INERTIA,:CLUSTER_SHADE,:CLUSTER_PROMINENCE)";*/
//
//	queryStr1  = "insert into IMAGE_TEST values" 
//	    "(:1,:2,:3,:4,:5,:6,"
//		":7,:8,:9,:10,:11,:12,:13,"
//		":14,:15,:16,:17,:18,:19,"
//		":20,:21,:22,:23,:24,:25,:26,:27,:28,"
//		":29,:30,:31,:32,:33";
//	queryStr = queryStr1 +")";
//
//	exeStatus = sqlite3_prepare(db, queryStr.c_str(), queryStr.size(), &ppStmt, &tail);
//	if( exeStatus != SQLITE_OK )
//		fprintf(stderr, "IMAGE_TEST:Error during sqlite3_prepare: %s\n", sqlite3_errmsg(db));
//	else
//	{
//		int cell_id = 1;
//	//	for (int i = 1 ; i<=ncol ; i++)
//		//{  // string intTOchar;
//		    char numstr[21];
//		   // intTOchar = itoa(i,numstr,10);
//			//std::string counter_string;
//			//counter_string = ":" + intTOchar; counter_string.c_str()
//
//			int i = 5;
//            const char *counter_string = itoa(i,numstr,10);
//			sqlite3_bind_double( ppStmt,sqlite3_bind_parameter_index(ppStmt, counter_string), 9.22);
//			/*sqlite3_bind_double( ppStmt,sqlite3_bind_parameter_index(ppStmt, ":1"), img_id);
//			sqlite3_bind_double( ppStmt,sqlite3_bind_parameter_index(ppStmt, ":2"), cell_id );
//			sqlite3_bind_double( ppStmt,sqlite3_bind_parameter_index(ppStmt, ":3"), 10.22);
//			sqlite3_bind_double( ppStmt,sqlite3_bind_parameter_index(ppStmt, ":4"), 10);
//			sqlite3_bind_double( ppStmt,sqlite3_bind_parameter_index(ppStmt, ":5"), 10.22);
//			sqlite3_bind_double( ppStmt,sqlite3_bind_parameter_index(ppStmt, ":6"), 10);
//			sqlite3_bind_double( ppStmt,sqlite3_bind_parameter_index(ppStmt, ":7"), 10.22);
//			sqlite3_bind_double( ppStmt,sqlite3_bind_parameter_index(ppStmt, ":8"), 10);
//			sqlite3_bind_double( ppStmt,sqlite3_bind_parameter_index(ppStmt, ":9"), 10.22);
//			sqlite3_bind_double( ppStmt,sqlite3_bind_parameter_index(ppStmt, ":10"), 10);
//			sqlite3_bind_double( ppStmt,sqlite3_bind_parameter_index(ppStmt, ":11"), 10.22);
//			sqlite3_bind_double( ppStmt,sqlite3_bind_parameter_index(ppStmt, ":12"), 10);
//			sqlite3_bind_double( ppStmt,sqlite3_bind_parameter_index(ppStmt, ":13"), 10.22);
//			sqlite3_bind_double( ppStmt,sqlite3_bind_parameter_index(ppStmt, ":14"), 10);
//			sqlite3_bind_double( ppStmt,sqlite3_bind_parameter_index(ppStmt, ":15"), 10.22);
//			sqlite3_bind_double( ppStmt,sqlite3_bind_parameter_index(ppStmt, ":16"), 10);
//			sqlite3_bind_double( ppStmt,sqlite3_bind_parameter_index(ppStmt, ":17"), 10.22);
//			sqlite3_bind_double( ppStmt,sqlite3_bind_parameter_index(ppStmt, ":18"), 10);
//			sqlite3_bind_double( ppStmt,sqlite3_bind_parameter_index(ppStmt, ":19"), 10.22);
//			sqlite3_bind_double( ppStmt,sqlite3_bind_parameter_index(ppStmt, ":20"), 10);
//			sqlite3_bind_double( ppStmt,sqlite3_bind_parameter_index(ppStmt, ":21"), 10.22);
//			sqlite3_bind_double( ppStmt,sqlite3_bind_parameter_index(ppStmt, ":22"), 10);
//			sqlite3_bind_double( ppStmt,sqlite3_bind_parameter_index(ppStmt, ":23"), 10.22);
//			sqlite3_bind_double( ppStmt,sqlite3_bind_parameter_index(ppStmt, ":24"), 10);
//			sqlite3_bind_double( ppStmt,sqlite3_bind_parameter_index(ppStmt, ":25"), 10.22);
//			sqlite3_bind_double( ppStmt,sqlite3_bind_parameter_index(ppStmt, ":26"), 10);
//			sqlite3_bind_double( ppStmt,sqlite3_bind_parameter_index(ppStmt, ":27"), 10.22);
//			sqlite3_bind_double( ppStmt,sqlite3_bind_parameter_index(ppStmt, ":28"), 10);
//			sqlite3_bind_double( ppStmt,sqlite3_bind_parameter_index(ppStmt, ":29"), 10.22);
//			sqlite3_bind_double( ppStmt,sqlite3_bind_parameter_index(ppStmt, ":30"), 10.22);
//			sqlite3_bind_double( ppStmt,sqlite3_bind_parameter_index(ppStmt, ":31"), 10.22);
//			sqlite3_bind_double( ppStmt,sqlite3_bind_parameter_index(ppStmt, ":32"), 10.22);
//			sqlite3_bind_double( ppStmt,sqlite3_bind_parameter_index(ppStmt, ":33"), 10.22);*/ 
//			
//		//}
//		exeStatus =sqlite3_step(ppStmt);
//		if( exeStatus != SQLITE_OK )
//		fprintf(stderr, "IMAGE_TEST:Error during sqlite3_step(ppStmt: %s\n", sqlite3_errmsg(db));
//         
//		sqlite3_reset(ppStmt);
//		//cell_id = cell_id + 1;
//		fprintf(stderr,"Inserted data for img id=%i ,cellid =%i\n", img_id,cell_id);
//		sqlite3_finalize(ppStmt);
//	}
//}


void GenericSearch(sqlite3 *db,const char * sql)
{ 
	int  nrow;       /* Number of result rows written here */
	int  ncol;       /* Number of result columns written here */
	char *zErr;      /* Error msg written here */
    int  exeStatus;
	char **result;
   
    exeStatus = sqlite3_get_table(db, sql, &result, &nrow, &ncol, &zErr);
	
	if( exeStatus != SQLITE_OK )
		fprintf(stderr, "Farsight Sql Error :Incorrect Sql Query: %s\n", sqlite3_errmsg(db));
	else
	{ fprintf(stderr, "Number of row:%i  and Columns:%i\n \n",nrow,ncol);
	 for(int i=0; i <= nrow; i++) 
	 {
	  for(int j=0; j < ncol; j++) 
	  {fprintf(stdout, "%s ", result[(i)*ncol + j]);
	   fprintf(stdout, ",");
	  }
	  fprintf(stdout, "\n \n");
     }
	}

    
   /* Code to write the query data in CSV file :: STARTS*/
     FILE * fp = fopen("mydata.csv", "w");  // opens file to write ("w"), using the FILE pointer "fp"
     for(int i=0; i <= nrow; i++) 
	 {
	  for(int j=0; j < ncol; j++) 
	  {fprintf(fp, "%s ", result[(i)*ncol + j]);
	   fprintf(fp, ",");
	  }
	  fprintf(fp, "\n");
	 }
	
    fclose(fp);  // close the file after you're done
   /* Code to write the query data in CSV file :: ENDS*/
	/* 
	  A result table should be deallocated using sqlite3_free_table()
	*/ 
	sqlite3_free_table(result);
}

 sqlite3 * sqliteGetConnection(const char *dbName)
{   
	int  exeStatus ;
	sqlite3 *dbConn; 

	exeStatus = sqlite3_open(dbName,&dbConn);
	
	if( exeStatus!=SQLITE_OK )
	  fprintf(stderr, "Can't open database: %s\n", sqlite3_errmsg(dbConn));
	else if(exeStatus == SQLITE_OK )
	  fprintf(stderr, "Open database %s successfully\n",dbName);

	return dbConn;
}

 void AlterTable(sqlite3 *dbConn, std::string table , std::string column)
{  	char *zErr; 
	int  exeStatus;
	std::string pSQL;
    
	//strcat( pSQL, "Alter table" );
	pSQL = "ALTER TABLE " + table + " ADD COLUMN " + column + " FLOAT DEFAULT '0.0'";

	char *pSQL_cstr = new char [pSQL.size()+1];
	strcpy(pSQL_cstr, pSQL.c_str());
	
    exeStatus = sqlite3_exec(dbConn, pSQL_cstr, 0, 0, &zErr);
    if( exeStatus != SQLITE_OK && exeStatus != SQLITE_DONE ){
		std::cerr<<"SQL error: "<< zErr << std::endl;
		std::cout<<"The statement is: "<< pSQL << std::endl;
		sqlite3_free(zErr);
    }
 
	//sqlite3_close(dbConn);
}

 void checkForUpdate(sqlite3 *dbConn, std::vector<std::string> column_names)
{   
	int  nrow;       /* Number of result rows written here */
	int  ncol;       /* Number of result columns written here */
	char *zErr;      /* Error msg written here */
	std::string MasterTableName1;
	char *MasterTableName;
    int  exeStatus;
	char **result;
	
	/*strcat( MasterTableName, project_name );
	strcat( MasterTableName, "_MASTER");*/
	MasterTableName1="IMAGE_TEST";
	MasterTableName = new char[MasterTableName1.length() + 1];
	strcpy( MasterTableName, MasterTableName1.c_str() );

	char *CheckforUpdateQuery;
	std::string CheckforUpdateQuery1;
    /*
	 1) Master table for each project must have format "Projectname_MASTER"
	 2) They must have the column names IMG_ID corresponding to IMAGE
	 3) We will insert a dummy record for 1st row with IMG_ID = 1 when we ship the database
	*/

	// make sure the table exists
	std::string createTableQuery = "CREATE TABLE IF NOT EXISTS ";
	createTableQuery += MasterTableName;
	createTableQuery += " (IMG_ID INTEGER, CELL_ID INTEGER,";
	createTableQuery += " IMG_NAME TEXT, IMG_LOCATION TEXT)";
	exeStatus = sqlite3_get_table(dbConn, createTableQuery.c_str(), &result, &nrow, &ncol, &zErr);
	if( exeStatus != SQLITE_OK )
	{
		std::cerr<<"Problem creating SQLite table: "<<sqlite3_errmsg(dbConn)<<std::endl;
	}

	//strcat( CheckforUpdateQuery, "Select * from" );
	CheckforUpdateQuery1 = "PRAGMA table_info(" + MasterTableName1 + ");";
	CheckforUpdateQuery = new char[CheckforUpdateQuery1.length() + 1];
	strcpy( CheckforUpdateQuery, CheckforUpdateQuery1.c_str() );
       
    exeStatus = sqlite3_get_table(dbConn, CheckforUpdateQuery, &result, &nrow, &ncol, &zErr);

	if( exeStatus != SQLITE_OK )
		std::cerr<<"Farsight Sql Error :Incorrect Sql Query: "<<sqlite3_errmsg(dbConn)<<std::endl;

	int ind = 7;
	std::cout<<"The column headers in the selected table on the database are:\n";
	for( int i=0; i<nrow; ++i ){
		std::string colname;
		colname = result[ind+i*6];
		std::cout<< colname << "\t";
	}

	for(std::vector<std::string>::iterator theit=column_names.begin(); theit!=column_names.end(); ++theit ){
		bool header_found = false;
		for( int i=0; i<nrow; ++i ){
			if( strcasecmp( (*theit).c_str() , result[ind+i*6] ) == 0 ){
				header_found = true;
				break;
			}
		}
		if( !header_found ){
			//char *cstr = new char [(*theit).size()+1];
			//strcpy (cstr,(ftk::GetStringInCaps(*theit)).c_str());
			//cstr = toupper(cstr);
			AlterTable( dbConn, MasterTableName1 , ftk::GetStringInCaps(*theit) );
		}
	}

   /*  
     1) The first 'ncol' values in 'result' are column names 
	 2) We can compare them to 'columnArray'array passed in this function.
	 3) Then we call function 'AlterTable(bConn,table ,column).
	    Note that we call this function each time for one column found
    */

	std::cerr << "Headers checked\n";

}


 int GenericInsert( sqlite3 *db, char * Img_name, const char * project_name,
	                char *path, std::vector<double> ProcessedImgRowArray,
					int ncol, int nrow, std::vector<std::string> colNames )
{ 
	sqlite3_stmt *ppStmt;
	std::string queryStr,queryStrPart1 ,queryStrPart2 = "?",table_name;
	//char  *img_name;
	int   exeStatus,img_id;
	const char *tail;
	int  qnrow=0;       // Number of result rows written here
	int  qncol=0;       // Number of result columns written here 
	char *zErr;      // Error msg written here 
	char **result;
//	char  *zErrMsg,*name;
    //---------------------------------------------------------------------------------------------------------
    /*PART1)Every Project will have this generic table called IMAGE 
	        which hold image name and autogenerted IMG_ID*/
	//Check if image already exists in the list
	std::string check_str;
	check_str = "SELECT * FROM ";
	check_str += project_name;
	check_str += " WHERE IMG_LOCATION='";
	char *check_str_cstr = new char [check_str.size()+strlen(path)+1];
	strcpy (check_str_cstr, check_str.c_str());
	strcat( check_str_cstr, path );
	strcat( check_str_cstr, "'" );
	exeStatus = sqlite3_get_table(db, check_str_cstr, &result, &qnrow, &qncol, &zErr);
	if( qnrow ){
		img_id = atoi( result[qnrow*qncol] );
		check_str.clear();
		check_str = "DELETE FROM ";
		check_str += project_name;
		check_str += " WHERE IMG_ID="+ftk::NumToString(img_id);
		sqlite3_prepare(db, check_str.c_str(), check_str.size(), &ppStmt, &tail);
		exeStatus = sqlite3_step(ppStmt);
		if( exeStatus != SQLITE_OK && exeStatus != SQLITE_DONE ){
			std::cerr<<"Sqlite3_prepare error:"<<sqlite3_errmsg(db)<<std::endl;
			return -1;
		}
	}
	else
	{
		std::string img_id_query;
		img_id_query = "SELECT MAX(IMG_ID) FROM ";
		img_id_query += project_name;
		exeStatus = sqlite3_get_table(db, img_id_query.c_str(), &result, &qnrow, &qncol, &zErr);
		if( qnrow > 1)
		{
			img_id = atoi( result[qnrow*qncol] );
			++img_id;
		}
		else
		{
			img_id = 1;
		}
		queryStr = "insert into ";
		queryStr += project_name;
		queryStr += " (IMG_ID, IMG_NAME, IMG_LOCATION) values (:IMG_ID, :IMG_NAME, :IMG_LOCATION)";
		//Location variable not decided
		sqlite3_prepare(db, queryStr.c_str(), queryStr.size(), &ppStmt, &tail);
		sqlite3_bind_double( ppStmt,sqlite3_bind_parameter_index(ppStmt,":IMG_ID"), img_id);
		sqlite3_bind_text( ppStmt,sqlite3_bind_parameter_index(ppStmt,":IMG_NAME"),Img_name,strlen(Img_name),SQLITE_TRANSIENT );
		sqlite3_bind_text( ppStmt,sqlite3_bind_parameter_index(ppStmt,":IMG_LOCATION"),path,strlen(path),SQLITE_TRANSIENT );
		exeStatus = sqlite3_step(ppStmt);
		sqlite3_finalize(ppStmt);
		if( exeStatus != SQLITE_OK && exeStatus != SQLITE_DONE ){
			std::cerr<<"Sqlite3_prepare error:"<<sqlite3_errmsg(db)<<std::endl;
			return -1;
		}
		exeStatus = sqlite3_get_table(db, check_str_cstr, &result, &qnrow, &qncol, &zErr);
		if( exeStatus != SQLITE_OK && exeStatus != SQLITE_DONE ){
			std::cerr<<"Sqlite3_prepare error:"<<sqlite3_errmsg(db)<<std::endl;
			return -1;
		}
	}


    //---------------------------------------------------------------------------------------------------------
	/*PART2) Then we retrive the autogenerted IMG_ID and insert the data corresponding to that
	        id in master data table in Part3*/
	
	//queryStr = "select IMG_ID from IMAGE where IMG_LOCATION = :IMG_LOCATION";         
	//exeStatus = sqlite3_prepare(db, queryStr.c_str(), queryStr.size(), &ppStmt, &tail);
	//if( exeStatus != SQLITE_OK ){
	//	std::cerr<<"Sqlite3_prepare error:"<<sqlite3_errmsg(db)<<std::endl;
	//	return -1;
	//}
	//else{
	//	sqlite3_bind_text(ppStmt,sqlite3_bind_parameter_index(ppStmt,":IMG_NAME"),path,strlen(path),SQLITE_TRANSIENT );
	//	exeStatus = sqlite3_step(ppStmt);
	//	while(exeStatus == SQLITE_ROW){
	//		img_id   = sqlite3_column_int(ppStmt, 0);
	//		img_name = (char*)sqlite3_column_text(ppStmt,1);
	//		exeStatus = sqlite3_step(ppStmt);
	//		std::cerr<<"img_id="<<img_id<<std::endl;
	//	}
	//	sqlite3_finalize(ppStmt);
	//}

    //---------------------------------------------------------------------------------------------------------
	/*PART3) 1)Now we will insert the data in master data table
	         2) For this we generate parameterised queryStr  
			 3) Refer page 224 of book i have sent in email*/

	table_name.clear();

	std::string proj_nm = project_name;
	table_name = "IMAGE_TEST";
	
	queryStrPart1 = "insert into "+table_name+"(IMG_ID, CELL_ID";//+" values";
	for( int i=0; i<(int)colNames.size(); ++i )
		queryStrPart1 = queryStrPart1 + ", " + ftk::GetStringInCaps( colNames.at(i) );
	queryStrPart1 = queryStrPart1+")values(:IMG_ID, :CELL_ID";
	for( int i=0; i<(int)colNames.size(); ++i )
		queryStrPart1 = queryStrPart1 + ", :" + ftk::GetStringInCaps( colNames.at(i) );
	queryStrPart1 = queryStrPart1 + ")";
	//std::cout<<queryStrPart1<<std::endl;


    //Save the ncols for the table in varibale when calling 'checkForUpdate'
    //for (int i=1;i<ncol; ++ncol){	
	//	queryStrPart2 = queryStrPart2 +",?" ; // 
    //}

	queryStr = queryStrPart1;//+"("+queryStrPart2+")";
	
	exeStatus = sqlite3_prepare(db, queryStr.c_str(), queryStr.size(), &ppStmt, &tail);
	if( exeStatus != SQLITE_OK && exeStatus != SQLITE_DONE )
		std::cerr<<"IMAGE_TEST:Error during sqlite3_prepare: "<<sqlite3_errmsg(db)<<std::endl;
	else
	{
		//int cell_id = 1;img_id
		for(int i=0; i<nrow; ++i){
			sqlite3_bind_double( ppStmt,sqlite3_bind_parameter_index(ppStmt,":IMG_ID"), img_id);
			sqlite3_bind_double( ppStmt,sqlite3_bind_parameter_index(ppStmt, ":CELL_ID"), ProcessedImgRowArray[i*ncol] );
			for (int j=1; j<ncol; ++j){
				std::string col_nm;
				col_nm = ":" + ftk::GetStringInCaps( colNames.at(j-1) );
				sqlite3_bind_double( ppStmt,sqlite3_bind_parameter_index(ppStmt, col_nm.c_str() ), ProcessedImgRowArray[i*ncol+j] );
			}
			exeStatus =sqlite3_step(ppStmt);
			if( exeStatus != SQLITE_OK && exeStatus != SQLITE_DONE )
				std::cerr<<"IMAGE_TEST:Error during sqlite3_step(ppStmt: "<<sqlite3_errmsg(db)<<std::endl;
			sqlite3_reset(ppStmt);
		}
		sqlite3_finalize(ppStmt);
	}
	return img_id;
}


void sqliteCloseConnection(sqlite3 *dbConn)
{   /*
	  sqlite3_close() routine is the destructor for the sqlite3 object. 
	  Calls to sqlite3_close() return SQLITE_OK if the sqlite3 object is successfullly destroyed 
	  and all associated resources are deallocated.
	*/
	sqlite3_close(dbConn);
	fprintf(stderr,"Execution complete!\n");

}


//void Testing(sqlite3 *db,const char * sql,char *project_name)
//{ 
//	int  nrow;       /* Number of result rows written here */
//	int  ncol;       /* Number of result columns written here */
//	char *zErr;      /* Error msg written here */
//    int  exeStatus;
//	char **result;
//
//	sqlite3_stmt *ppStmt;
//	string queryStr,queryStrPart1 ,queryStrPart2 = "?",table_name;
//	char  *name, *img_name;
//	const char *tail;
//	char  *zErrMsg ;
//   
//    exeStatus = sqlite3_get_table(db, "select * from IMAGE_TEST", &result, &nrow, &ncol, &zErr);
//	
//	if( exeStatus != SQLITE_OK )
//		fprintf(stderr, "Farsight Sql Error :Incorrect Sql Query: %s\n", sqlite3_errmsg(db));
//	
//	sqlite3_free_table(result);
//
//
//	
//	//table_name = project_name + "_MASTER"
//	std::string proj_nm = project_name;
//	queryStrPart1 = "insert into"+proj_nm+"values";
//
//    //Save the ncols for the table in varibale when calling 'checkForUpdate'
//    for (int i=1;i<ncols;++ncols){	
//		queryStrPart2 = queryStrPart2 +",?" ; // 
//    }
//    
//	queryStr = queryStrPart1+"("+queryStrPart2+")";
//	
//	exeStatus = sqlite3_prepare(db, queryStr.c_str(), queryStr.size(), &ppStmt, &tail);
//	if( exeStatus != SQLITE_OK )
//		fprintf(stderr, "IMAGE_TEST:Error during sqlite3_prepare: %s\n", sqlite3_errmsg(db));
//	else
//	{
//		int cell_id = 1;
//		for (int i = 1 ; i<=ncol ; i++)
//		{  sqlite3_bind_double( ppStmt,i, 1.2);
//		}
//
//		exeStatus =sqlite3_step(ppStmt);
//		if( exeStatus != SQLITE_OK )
//		   fprintf(stderr, "IMAGE_TEST:Error during sqlite3_step(ppStmt: %s\n", sqlite3_errmsg(db));
//        //sqlite3_reset(ppStmt);
//		sqlite3_finalize(ppStmt);
//		cell_id = cell_id + 1;//Use this in ur C++ code when ur invoking another row of data
//	}
//
//		sqlite3_close(db);
//}

sqlite3 * sqliteOpenConnection()
{   
	int  exeStatus ;
	sqlite3 *dbConn; 

	char cCurrentPath[FILENAME_MAX];
	GetCurDir(cCurrentPath, sizeof(cCurrentPath));
	cCurrentPath[sizeof(cCurrentPath) - 1] = '/';

	std::cout<<cCurrentPath<<"/NE.s3db"<<std::endl;
	exeStatus = sqlite3_open("./NE.s3db",&dbConn);

 	// int sqlite3_open(
    //    const char *filename,  Database filename (UTF-8) 
    //     sqlite3 **ppDb         OUT: SQLite db handle 
    // )
	
	if( exeStatus!=SQLITE_OK )
		std::cerr << "Can't open database: " << sqlite3_errmsg(dbConn) << std::endl;
	else if(exeStatus == SQLITE_OK )
		std::cerr << "Open database successfully\n";
	return dbConn;
}


void sqliteExecuteQuery2(sqlite3 *db,const char * sql, std::string path, int qnum)
{ 
	int  nrow;       // Number of result rows written here
	int  ncol;       // Number of result columns written here 
	char *zErr;      // Error msg written here 
    int  exeStatus;
	char **result;
   
    exeStatus = sqlite3_get_table(db, sql, &result, &nrow, &ncol, &zErr);
	
	//int sqlite3_get_table(
	//	sqlite3 *db,           An open database 
	//	const char *sql,      SQL to be evaluated 
	//	char ***result,     Results of the query 
	//	int *Row,            Number of result rows written here 
	//	int *Column,         Number of result columns written here 
	//	char **Errmsg        Error msg written here 
	//	)
	  

	if( exeStatus != SQLITE_OK )
		std::cerr << "Farsight Sql Error :Incorrect Sql Query: " << sqlite3_errmsg(db) << std::endl;
	else{
		std::cerr << "Number of row: "<< nrow <<" and Columns: " << ncol << std::endl;
	 //for(int i=0; i <= nrow; i++) 
	 //{
	 //for(int j=0; j < ncol; j++) 
	 //{fprintf(stdout, "%s ", result[(i)*ncol + j]);
	 //fprintf(stdout, ",");
	 //}
	 //fprintf(stdout, "\n \n");
     //}
	}
	path = path + "FarsightData" + ftk::NumToString(qnum).c_str() + ".csv";
	FILE * fp = fopen(path.c_str(),"w");

     for(int i=0; i <= nrow; i++) 
	 {
	  for(int j=0; j < ncol; j++) 
	  {fprintf(fp, "%s ", result[(i)*ncol + j]);
	   fprintf(fp, ",");
	  }
	  fprintf(fp, "\n");
     }
	 fclose(fp);

	 
	//  A result table should be deallocated using sqlite3_free_table()
	 
	sqlite3_free_table(result);
}

} //End namespace ftk

