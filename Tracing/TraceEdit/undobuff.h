
#ifndef UNDOBUFFER_H
#define UNDOBUFFER_H

#include <iostream>

template <class T>
class undoBuffer
{
	private:
		unsigned int count;
		unsigned int size;
		T *buff;
	public:
		static const int UNDO = 0;
		static const int REDO = 1;		
		undoBuffer();
		undoBuffer(int sz);
		~undoBuffer();
		void Add(T newObject);
		bool UndoOrRedo(int flag);
		void print();
		unsigned int getCurrPos();
		unsigned int getNextPos();
		unsigned int getLastPos();
		T getState();
		unsigned int upperBound;
		unsigned int lowerBound;
		bool filledFlag;
};

#endif

template <class T>
undoBuffer<T>::undoBuffer()
{
	count = 0;
	size = 10;
	buff = new T[size];
	upperBound = lowerBound = 0;
	filledFlag = 0;
}

template <class T>
undoBuffer<T>::undoBuffer(int sz)
{
	count = 0;
	size = sz;
	buff = new T[size];
	upperBound = lowerBound = 0;
	filledFlag = 0;
}


template <class T>
undoBuffer<T>::~undoBuffer()
{
	delete [] buff;
}

template <class T>
void undoBuffer<T>::Add(T newObject)
{
	unsigned int pos = getNextPos();
	std::cout << "Adding object to position " << pos << std::endl;
	buff[pos] = newObject;
	count++;
	
	//count = ++count % 10;
	
	if(count == 10)
	{
		filledFlag = 1;
	}
	if(!filledFlag)
	{
		lowerBound = 0;
		upperBound = getNextPos();
	}
	else
	{
	upperBound = getNextPos();
	lowerBound = getCurrPos();
	}
	
}
template <class T>
bool undoBuffer<T>::UndoOrRedo(int flag)
{
	T* state;
	unsigned int currpos = getCurrPos();
	unsigned int nextpos = getNextPos();
	unsigned int lastpos = getLastPos();
	
	if(flag == UNDO)
	{
		//std::cout << "Undoing to position " << lastpos << " from position: " << currpos << std::endl;
		state = (buff + lastpos);
		if(state != NULL)
		{
			if(lastpos != lowerBound)
			{
				std::cout << "Undoing to position " << lastpos << " from position: " << currpos << std::endl;
				if(count == 0)
				{
					count = 9;
				}
				else
				{
					count --;
				}
				return 1;
			}
		}
		else return 0;
			
	}
	else if (flag == REDO)
	{
		//std::cout << "Redoing to position " << nextpos << " from position: " << currpos << std::endl;
		state = (buff + nextpos);
		if(state != NULL)
		{
			if(nextpos != upperBound)
			{
				std::cout << "Redoing to position " << nextpos << " from position: " << currpos << std::endl;
				count++;
				return 1;
			}
		}
		else return 0;
	}
	return 0;
}

template <class T>
void undoBuffer<T>::print()
{
	std::cout << "[";
	for(unsigned int i = 0; i < size; i++)
	{
		if((buff + i) != NULL)
		{
			std::cout << buff[i] << " ";
		}
		else
		{
			std::cout << "0" << " ";
		}
	}
	std::cout << "]" << " ";
	std::cout << upperBound << " " << lowerBound << std::endl;
}

template <class T>
unsigned int undoBuffer<T>::getCurrPos()
{
	if (count == 0)
	{
		count = 9;
	}
	return ((count - 1) % size);	
}

template <class T>
unsigned int undoBuffer<T>::getNextPos()
{
	return (count % size);
}

template <class T>
unsigned int undoBuffer<T>::getLastPos()
{
	unsigned int find;
	if(count == 1)
	{
		find = 9;
		return find;
	}
	else if(count == 0)
	{
		count = 9;	
	}
	return((count - 2) % size);
}

template <class T>
T undoBuffer<T>::getState()
{
	unsigned int pos = getCurrPos();
	return buff[pos];
}
	