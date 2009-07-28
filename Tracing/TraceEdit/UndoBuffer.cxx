#include"UndoBuffer.h"
template<typename S>
S UndoandGetState(UndoBuffer<S>* buff){
    if(buff->Undo() == 1) return buff->GetUndoPosition();
    else return NULL;
}
template<typename S>
S RedoandGetState(UndoBuffer<S>* buff){
    if(buff->Redo() == 1) return buff->GetRedoPosition();
    else return NULL;
}

/*
//This is test code
int main(int argc, char *argv[]){ 
	UndoBuffer<int>* buff = new UndoBuffer<int>(10);
	printf("add 15 elements\n");
	for(int i =0; i < 15; i++){
		buff->Add(i+1);
		buff->print();
	}
	printf("undo 8 times\n");
	for(int i =0; i < 8; i++){		
		int tempo = UndoandGetState(buff);
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
		int tempo = UndoandGetState(buff);
		printf("%i\n",tempo);
		buff->print();
	}
	printf("now Redo 6 times\n");
	for(int i =0; i < 9; i++){
		int tempo = RedoandGetState(buff);
		printf("%i\n",tempo);
		buff->print();
	}
	return 0; 
}*/
