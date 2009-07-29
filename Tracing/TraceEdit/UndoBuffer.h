#include <stdio.h>

template<class T>
class UndoBuffer {

public:
	static const int UNDO = 0;
	static const int REDO = 1;
	UndoBuffer(int Size){
		BufferSize = Size;
		Buffer = new T[BufferSize];
		RedoCounter = UndoCounter = RedoPosition = UndoPosition = 0; 
		emptyflag = true;
	}
//If the buffer is to be initialized with saved states. Not sure if such a functionality is need
/*	UndoBuffer<T>(T* Initializer, int Size){
		int temp = 0;
		this.BufferSize = Size;
		Buffer = new T[BufferSize];
		CurrentOpenPosition = temp =  sizeof(Initializer);
		UndoPosition = temp-1;
		for( int i = 0; i < temp; i++) Buffer[i] = Initializer[i];
	}

	UndoBuffer<T>(UndoBuffer<T>* Initializer){
		this.Buffer = Initializer.Buffer;
		this.BufferSize = Initializer.BufferSize;
		this.CurrentOpenPosition = Initializer.CurrentOpenPosition;
		this.UndoPosition = Initializer.UndoPosition;
		this.RedoPosition = Initializer.RedoPosition;
	}*/
//*******************************************************************************************************
//This function should be called everytime a change is made to the trace object by the user 
//*******************************************************************************************************
template<typename S> S UndoOrRedoandGetState(UndoBuffer<S>* buff,int Flag);
int	Add(T Object){
		T* NewObject = new T(Object); 
		Buffer[CurrentOpenPosition] = *NewObject;
		UndoPosition = CurrentOpenPosition;
		CurrentOpenPosition = MyMod((CurrentOpenPosition +1),BufferSize);
		if(UndoCounter != (BufferSize-1) && emptyflag == false) UndoCounter++;
		if(RedoCounter != 0) RedoCounter = 0;
		RedoPosition = UndoPosition;
		emptyflag = false;
		return 0;
	}
int	GetCurrentOpenPositionslot(){return CurrentOpenPosition;}
T	GetCurrentState(){return Buffer[MyMod(CurrentOpenPosition-1,BufferSize)];}
T	GetUndoPosition(){return Buffer[UndoPosition];}
T 	GetRedoPosition(){return Buffer[RedoPosition];}
int	Undo(){return DeleteAnObject();}
int	Redo(){return VisitOldFriend();}
void	print(){
		printf("[");	
		for(int i = 0; i < BufferSize; i++){
			printf("%i,",Buffer[i]);
		}
		printf("]  %i   %i   %i   %i  %i\n",CurrentOpenPosition,UndoPosition,UndoCounter ,RedoPosition,RedoCounter);

	}
private:
	T* Buffer;
	int BufferSize;
	int CurrentOpenPosition;
	int UndoPosition; 
	int RedoPosition;
	int UndoCounter;
	int RedoCounter;
	bool emptyflag;
	int	DeleteAnObject(){
		if(UndoCounter > 0){
			UndoPosition = MyMod((UndoPosition -1), BufferSize);
			UndoCounter--;
			CurrentOpenPosition = MyMod((CurrentOpenPosition -1), BufferSize);
			RedoPosition = UndoPosition; 
			if(RedoCounter != BufferSize-1) RedoCounter++;
			return 1;
		}
		else return -1;
	}
	int     VisitOldFriend(){
			if(RedoCounter > 0){
				RedoPosition = MyMod((RedoPosition +1), BufferSize);
				CurrentOpenPosition = MyMod((CurrentOpenPosition +1), BufferSize);
				UndoPosition =  MyMod((UndoPosition +1), BufferSize);
				if(UndoCounter != BufferSize-1) UndoCounter++;
				RedoCounter--;
				return 1;
			}
			else return -1;
		}
	int    MyMod(int a, int b){
			int temp = a;
			if (a >= b)
				for(;;) if(temp >= b) temp -= b; else break; 
			else if (a < 0)
				for(;;) if(temp < 0) temp += b; else break; 
			return temp;
		}
};