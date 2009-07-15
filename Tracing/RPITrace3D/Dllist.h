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

//////////////////////////////////////////////////////////////////////////
// File: Dllist.h
// 
// This file contains the declaration and definintion of a double linked 
// list template class.
//
#ifndef dllist_h
#define dllist_h

/////////////////////////////////////////////////////////////////////
// Class: CLNode
// 
// a template class that represents a node in a double linked list. 
// a node contains two pointers pointing to the objects to its left
// and right. It also contain a pointer to the data
template <class DataType>
class CLNode
{
public:
	// Default CTOR
	CLNode() : before(0), after(0), data(0), OwnMyData(1)
	{
	}
	// CTOR
	CLNode(DataType* aData) : before(0), after(0), OwnMyData(1)
	{
		data = new DataType(*aData);
	}
	// copy CTOR
	CLNode(CLNode& aNode) : before(aNode.before), after(aNode.after),
		OwnMyData(1)
	{
		data = new DataType(*aNode.data);
	}

	// DTOR
	~CLNode()
	{
		if (OwnMyData)
		{
			if (data)
				delete data;

			// adjust the pointers of the nodes to my left and right
			if (before)
				before->after = after;
			if (after)
				after->before = before;
		}
	}
	// set my data to the given pointer
	// valid only if I do not own my data
	inline void SetData(DataType* aData)
	{
		if (OwnMyData)
			cout << "ERROR: the object own its data. You can not "
				<< "use this method" << endl;
		else
			data = aData;
	}

	// pointers to other nodes
	CLNode* before;
	CLNode* after;
	DataType* data;

	// Set to 1 if I own my data. This is the default behavior
	int OwnMyData;
};


//////////////////////////////////////////////////////////////////////
// Class: Double Linked List
//
// a template class representing a doubly linked list
template <class DataType>
class CDLList
{
public:

	// default CTOR
	CDLList() : head(0), tail(0), current(0), length(0), OwnMyData(1),
		m_lUserFlag(0)
	{
	}
	// CTOR with a single element
	CDLList(DataType* aData)
	{
		head = new CLNode<DataType>(aData);
		tail = head;
		current = head;
		length = 1;
		OwnMyData = 1;
		m_lUserFlag = 0;
	}

	// Copy CTOR
	CDLList(CDLList<DataType>& rhs) : head(0), tail(0), length(0),
		OwnMyData(1), m_lUserFlag(0)
	{
		SetData(rhs);
	}

	// DTOR
	~CDLList()
	{
		ClearData();
	}

	// copy operator
	CDLList<DataType>& operator =(CDLList<DataType>& rhs)
	{
		// clear my data
		ClearData();
		SetData(rhs);

		return *this;
	}

	void ClearData()
	{
		CLNode<DataType>* temp;
		while (head)
		{
			temp = head;
			head = head->after;
			if (head)
				head->before = 0;
			delete temp;
		}
		head = 0;
		tail = 0;
		length = 0;
		OwnMyData = 1;
		m_lUserFlag = 0;
	}

	// copy the data of rhs to this
	void SetData(CDLList<DataType>& rhs)
	{
		CLNode<DataType>* temp = rhs.tail;
		while (temp)
		{
			AddElementOnTop(temp->data);
			temp = temp->before;
		}
		length = rhs.length;
		m_lUserFlag = rhs.m_lUserFlag;
	}

	// Add element to the list. (on top)
	inline void AddSortedDontCopy(DataType* anElem)
	{
		if (OwnMyData)
		{
			cout << "Error: This object has its own copy of all its Data"
				<< " members. Try using AddElementOnTop or AddElementOnEnd"
				<< "instead. Nothing is added. " << endl;
		}

		CLNode<DataType>* temp = new CLNode<DataType>;
		// the node does not own the data
		temp->OwnMyData = 0;
		temp->SetData(anElem);

		if (head)
		{
			// find the position were to add
			CLNode<DataType>* temp2 = head;
			while ((temp2 && (*temp2->data) > (*temp->data)))
				temp2 = temp2->after;

			// the element belong at the end
			if (temp2 == NULL)
			{
				tail->after = temp;
				temp->before = tail;
				tail = temp;
			}
			else
			{
				// temp2 points at where the point should be added
				temp->after = temp2;

				// if tempe is not the first point
				if (temp2->before)
				{
					temp2->before->after = temp;
					temp->before = temp2->before;
				}
				else
				{
					head = temp;
				}
				temp2->before = temp;
			}
		}
		else
		{
			head = temp;
			tail = temp;
		}
		length++;
		current = head;
	}
	// Add an element to the list (always add on on the head)
	// always make your own copy
	void Add(DataType* anElem)
	{
		CLNode<DataType>* temp = new CLNode<DataType>(anElem);
		if (head)
		{
			temp->after = head;
			head->before = temp;
			head = temp;
		}
		else
		{
			head = temp;
			tail = temp;
		}
		length++;
		current = head;
	}
	// Add an element to the list (always add on on the head)
	// always make your own copy
	void AddElementOnTop(DataType* anElem)
	{
		CLNode<DataType>* temp = new CLNode<DataType>(anElem);
		if (head)
		{
			temp->after = head;
			head->before = temp;
			head = temp;
		}
		else
		{
			head = temp;
			tail = temp;
		}
		length++;
		current = head;
	}
	// Add an element to the list (at the end)
	// always make your own copy
	void AddElementOnEnd(DataType* anElem)
	{
		CLNode<DataType>* temp = new CLNode<DataType>(anElem);
		if (tail)
		{
			temp->before = tail;
			tail->after = temp;
			tail = temp;
		}
		else
		{
			head = temp;
			tail = temp;
		}
		length++;
		current = head;
	}
	void AddNodeOnTop(CLNode<DataType>* aNode)
	{
		if (head)
		{
			aNode->after = head;
			head->before = aNode;
			head = aNode;
			head->before = NULL;
		}
		else
		{
			head = aNode;
			tail = aNode;
			head->before = NULL;
			tail->after = NULL;
		}
		length++;
		current = head;
	}

	void ExtendListOnTopFromTop(CDLList<DataType>& aList)
	{
		CLNode<DataType>* temp = aList.head;
		CLNode<DataType>* temp2 = temp->after;
		while (temp)
		{
			AddNodeOnTop(temp);
			temp = temp2;
			if (temp2)
				temp2 = temp2->after;
		}
	}
	void ExtendListOnTopFromEnd(CDLList<DataType>& aList)
	{
		CLNode<DataType>* temp = aList.tail;
		CLNode<DataType>* temp2 = temp->before;
		while (temp)
		{
			AddNodeOnTop(temp);
			temp = temp2;
			if (temp2)
				temp2 = temp2->before;
		}
	}

	void ExtendListOnEndFromTop(CDLList<DataType>& aList)
	{
		CLNode<DataType>* temp = aList.head;
		CLNode<DataType>* temp2 = temp->after;
		while (temp)
		{
			AddNodeOnEnd(temp);
			temp = temp2;
			if (temp2)
				temp2 = temp2->after;
		}
	}
	void ExtendListOnEndFromEnd(CDLList<DataType>& aList)
	{
		CLNode<DataType>* temp = aList.tail;
		CLNode<DataType>* temp2 = temp->before;
		while (temp)
		{
			AddNodeOnEnd(temp);
			temp = temp2;
			if (temp2)
				temp2 = temp2->before;
		}
	}

	void AddNodeOnEnd(CLNode<DataType>* aNode)
	{
		if (tail)
		{
			aNode->before = tail;
			tail->after = aNode;
			tail = aNode;
			tail->after = NULL;
		}
		else
		{
			head = aNode;
			tail = aNode;
			head->after = NULL;
			head->before = NULL;
		}
		length++;
		current = head;
	}
	// Print the nodes in the list
	void Print(ostream& out = cout)
	{
		CLNode<DataType>* temp = head;
		while (temp)
		{
			temp->data->Print(out);
			temp = temp->after;
		}
	}
	void inline Reset()
	{
		current = head;
	}

	DataType* GetNext()
	{
		DataType* temp = NULL;
		if (current)
		{
			temp = current->data;
			current = current->after;
		}
		return temp;
	}
	//
	// get the point in the middle of the linked list
	DataType* GetMiddlePoint()
	{
		int count = 0;
		int pointIndex = length / 2.0;
		CLNode<DataType>* temp = head;
		while (temp)
		{
			if (count == pointIndex)
				break;

			count++;
			temp = temp->after;
		}
		if (temp && count == pointIndex)
			return temp->data;
		else
			return NULL;
	}

	int IsMember(DataType& rhs)
	{
		CLNode<DataType>* temp = head;
		while (temp)
		{
			if (*temp->data == rhs)
				return 1;

			temp = temp->after;
		}
		return 0;
	}

	// break my list into two at the given point. return a pointer to the
	// head of the new list if succesfull, NULL otherwise
	CDLList<DataType>* BreakIntoTwo(CLNode<DataType>* aNode)
	{
		// the list is unbreakable if any of the conditions is satisfied
		if (! aNode || aNode == head || aNode == tail)
			return NULL;

		CDLList<DataType>* result = NULL;
		CLNode<DataType>* temp = head;
		int segmentLength = 0;
		while (temp)
		{
			if (temp == aNode)
			{
				result = new CDLList<DataType>;
				result->head = temp;
				result->tail = tail;

				CLNode<DataType>* ExtraNode = new CLNode<DataType>(temp->data);
				tail = temp->before;
				AddNodeOnEnd(ExtraNode);
				result->head->before = NULL;
				result->length = length - segmentLength;

				tail->after = NULL;
				length = segmentLength;

				break;
			}

			temp = temp->after;
			segmentLength++;
		}

		return result;
	}

	CDLList<DataType>& operator=(CDLList<DataType>* rhs)
	{
		if (rhs)
		{
			head = rhs->head;
			tail = rhs->tail;
			length = rhs->length;
			current = head;
		}
		return *this;
	}

	// data
	CLNode<DataType>* head;
	CLNode<DataType>* tail;
	CLNode<DataType>* current;
	int length;

	// this flag is very important. Since this class provides two
	// methods to add to the list, one adds the given pointer as is
	// and the others makes their own copy, it is very important
	// that we do not use both methods for a given class. This flag
	// serves this purpose, and is set if the instant owns its 
	// data members (default behavior)
	int OwnMyData;

	// this data member is initialized to zero and can be used
	// for whatever purpose the user required
	long m_lUserFlag;
};

#endif
