#include"UndoBuffer.h"

//Flags arew defined in the UndoBuffer.h file
template<typename S>
S UndoOrRedoandGetState(UndoBuffer<S>* buff,int Flag){
    	int ValCheck = 0;
    	S State = NULL;
    	if (Flag == buff->UNDO){
   		if(ValCheck = (buff->Undo() == 1)) State = buff->GetUndoPosition();
   		else if(ValCheck == -1) return NULL;
	}
    	else if (Flag == buff->REDO){
		if(ValCheck = (buff->Redo() == 1)) State = buff->GetUndoPosition();
   		else if(ValCheck == -1) return NULL;
	}
  	else return NULL;
	
   	S ImageState(State); 	//ImageState
   	return State;
}



//This is test code need to mark out the line ImageState if not working with objects. 
/*int main(int argc, char *argv[]){ 
	UndoBuffer<int>* buff = new UndoBuffer<int>(10);
	printf("add 15 elements\n");
	for(int i =0; i < 15; i++){
		buff->Add(i+1);
		buff->print();
	}
	printf("undo 8 times\n");
	for(int i =0; i < 8; i++){		
		int tempo = UndoOrRedoandGetState(buff,buff->UNDO);
		printf("%i\n",tempo);
		buff->print();
	}
	printf("add 5 more elements\n");
	for(int i =0; i < 5; i++){
		buff->Add((i+1)+16);
		buff->print();
	}
	printf("now undo 6 times\n");
	for(int i =0; i < 8; i++){
		int tempo = UndoOrRedoandGetState(buff,buff->UNDO);
		printf("%i\n",tempo);
		buff->print();
	}
	printf("now Redo 6 times\n");
	for(int i =0; i < 9; i++){
		int tempo = UndoOrRedoandGetState(buff,buff->REDO);
		printf("%i\n",tempo);
		buff->print();
	}
	return 0; 
}*/
