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

///////////////////////////////////////////////////////////////////////
// File: Queue.h
//
// This file contains the declaration for queue and point classes. 
// These classes are used to keep track of edge and centerline points.
//
#ifndef queue_h
#define queue_h

#include <stdlib.h>
#include <iostream>


///////////////////////////////////////////////////
// a class used to encapsulate the data of a queue. 
// This class is for the sole use of CQueue and should not be used
// by anybody else. I can enforce that by declaring all its member 
// functions (including CTOR) as private, but no need.
template <class DataType>
class CQNode
{
public:
	// CTOR. take the object pointed to by the argument as
	// mine. 
	CQNode(DataType *aData) : data(aData), next(0) {}
	// copy CTOR
	CQNode(CQNode &aNode)
		: data(0), next(0)
	{
		if(aNode.data)
			data = new DataType(*(aNode.data));
	}

	~CQNode() 
	{
		if(data)
			delete data;
		if(next)
			delete next;
	}

	// data
	DataType *data;
	CQNode *next;
};

template <class DataType>
class CQueue
{
public:
	// default CTOR, empty queue
	CQueue(): head(0), length(0) {}

	// copy CTOR
	CQueue(CQueue &aQueue)
		: head(0), length(0)
	{
		CQNode<DataType> *NodePtr = aQueue.head;
		while(NodePtr)
		{
			Push(*(NodePtr->data));
			NodePtr = NodePtr->next;
		}
	}
	// delete the entire queue
	~CQueue()
	{
		if(head)
			delete head;
	}

	// Add a point to the top of the queue
	inline void Push(DataType &point)
	{
		CPoint *temp1 = new CPoint(point);
		CQNode<DataType> *temp2 = new CQNode<DataType>(temp1);
		temp2->next = head;
		head = temp2;
		++length;
	}
	inline DataType *Pop()
	{
		DataType *result = 0;
		if(length > 0)
		{
			length--;
			result = head->data;
			head = head->next;
		}
		return result;
	}
	// retrun the length of the queue
	inline int GetLength() { return length; }

	// print the queue
	void Print(ostream &out=cout)
	{
		CQNode<DataType> *temp = head;
		while(temp)
		{
			temp->data->Print(out);
			temp = temp->next;
		}
	}

	// empty the queue
	void Flush()
	{
		if(length > 0)
		{
			delete head;
		}
		head = NULL;
		length = 0;
	}

	// data
	CQNode<DataType> *head;
	int length;
	
};

#endif
