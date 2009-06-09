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

#include <iostream>
#include "GraphCut.h"
#include "mex.h"
#include "matrix.h"

/* memory management */
void* operator new(size_t size)
{
    void *ptr = NULL;
//    mexWarnMsgTxt("Overloaded new operator");
    ptr = malloc(size);//ptr = mxMalloc(size);
    //mxAssert(ptr!=NULL,"memory allocation error");
    //mexMakeMemoryPersistent(ptr);
    return ptr;
}
void* operator new[](size_t size)
{
    void *ptr = NULL;
//    mexWarnMsgTxt("Overloaded new[] operator");
    ptr = malloc(size);//ptr = mxMalloc(size);
    //mxAssert(ptr!=NULL,"memory allocation error");
    //mexMakeMemoryPersistent(ptr);
    return ptr;
}
void operator delete(void* ptr)
{
//    mexWarnMsgTxt("Overloaded delete operator");
    free(ptr);//mxFree(ptr);
}
void operator delete[](void* ptr)
{
//    mexWarnMsgTxt("Overloaded delete[] operator");
    free(ptr);//mxFree(ptr);
}
