#ifndef _NESqliteFactory_
#define _NESqliteFactory_
#include <iostream>
#include <string>
#include <vector>
#include "SQLite/sqlite3.h"
#include "ftkCommon/ftkUtils.h"

namespace ftk
{
void GenericSearch(sqlite3 *db,const char * sql);
sqlite3 * sqliteGetConnection(const char *dbName);
sqlite3 * sqliteOpenConnection();
void AlterTable(sqlite3 *dbConn,char * table ,char *column);
void checkForUpdate(sqlite3 *dbConn, std::vector<std::string> column_names);
int GenericInsert( sqlite3 *db, char * Img_name, const char * project_name, char *path, std::vector<double> ProcessedImgRowArray, int ncol, int nrow, std::vector<std::string> colNames );
void sqliteCloseConnection(sqlite3 *dbConn);
void sqliteExecuteQuery2(sqlite3 *db,const char * sql, std::string path, int qnum);
}

#endif