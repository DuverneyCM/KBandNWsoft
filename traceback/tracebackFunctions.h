#include <stdio.h>
#include <stdlib.h>

//***	Declaration of Defines
#define dimPackSeq		10000000
#define DIAG			3
#define LEFT			2
#define UP				1
#define BOTH			0
#define NOTARROW		0

//***   Declaration of Functions
int getFileSize(FILE *inputFile);
int isFasta(FILE *inputFile);
int getPosFirstEnter(FILE *inputFile, int dimFile, int dimSearchWindow);
int getLineSize(FILE *inputFile, int dimFile, int posFirstEnter, int dimSearchWindow);
int getLastLineSize(int dimFile, int posFirstEnter, int dimLine);
int getDimSeq(int dimFile, int posFirstEnter, int dimLine);
int getNoArrowRows(int dimFile, int noPEs);
int getArrowRowFromFILE(FILE *inputFile, int indexRow, int noPEs, int *arrowRow);
int getLastDirection(int noPEs, int *arrowRow);
int getPosLastArrow(int noPEs, int *arrowRow);


int getArrow(int noPEs, int *arrowRow, int *posArrow, int *dirArrow, int *currentArrowRow, int posSeqA, int posSeqB, char tempSymA, char tempSymB, int lastARROW);
	int filterLastArrow(int posSeqA, int posSeqB, int ARROW, char lastSymA, char lastSymB, int lastArrow);
	void decoArrow(int* posArrow, int* arrowRow, int* dirHV, int ARROW);



//***   Definition of Functions
int maxAB(int A, int B) {
	if (A > B) return A; else return B;
}
int minAB(int A, int B) {
	if (A < B) return A; else return B;
}

int getFileSize(FILE *inputFile) {
    int beginFile, endFile, dimFile;
    fseek(inputFile, 0, SEEK_SET);
        beginFile = ftell(inputFile);
    fseek(inputFile, 0, SEEK_END);
        endFile = ftell(inputFile);
	dimFile = endFile - beginFile;
    return dimFile;
}

int isFasta(FILE *inputFile) {
	char firstChar = 0;
	fseek(inputFile, 0, SEEK_SET);
	fread(&firstChar,1,1,inputFile);
	if (firstChar == '>') return 1;
	else return 0;
}

int getPosFirstEnter(FILE *inputFile, int dimFile, int dimSearchWindow) {
	int dimSearchPacket;
	char SearchPacket[dimSearchWindow];
	if (dimFile > dimSearchWindow)
			dimSearchPacket = dimSearchWindow;
	else	dimSearchPacket = dimFile;
	fseek(inputFile, 0, SEEK_SET);
	fread(&SearchPacket,1,dimSearchPacket,inputFile);
	for (int i=0; i<dimSearchPacket; i++){
		if (SearchPacket[i] == '\n') return i;
	}
	return 0;
}

int getLineSize(FILE *inputFile, int dimFile, int posFirstEnter, int dimSearchWindow) {
	int dimSearchPacket;
	char SearchPacket[dimSearchWindow];
	if (dimFile-posFirstEnter > dimSearchWindow)
			dimSearchPacket = dimSearchWindow;
	else	dimSearchPacket = dimFile-posFirstEnter;
	fseek(inputFile, posFirstEnter+1, SEEK_SET);
	fread(&SearchPacket,1,dimSearchPacket,inputFile);
	for (int i=0; i<dimSearchPacket; i++){
		if (SearchPacket[i] == '\n') return i;
	}
	return dimFile-posFirstEnter;
}

int getLastLineSize(int dimFile, int posFirstEnter, int dimLine) {
	int seqSegment = dimFile - (posFirstEnter + 1);
	int noLines = seqSegment / (dimLine+1);
	int dimLastLine = seqSegment % (dimLine+1);
	return dimLastLine;
}

int getDimSeq(int dimFile, int posFirstEnter, int dimLine) {
	int seqSegment = dimFile - (posFirstEnter + 1);
	int noLines = seqSegment / (dimLine+1);
	int dimLastLine = seqSegment % (dimLine+1);
	int dimSeq = dimLine*noLines + dimLastLine;
	return dimSeq;
}

int getNoArrowRows(int dimFile, int noPEs) {
	return dimFile / (2*noPEs/8);
}

int getArrowRowFromFILE(FILE *inputFile, int indexRow, int noPEs, int *arrowRow){
	int isHV = 0, dirHV = 0;
	int noBytesByRow = 2*noPEs/8, dimUInt=sizeof(unsigned int);
	int noIntRegsByRow = noBytesByRow/dimUInt;
	int indexRowInBytes = indexRow*noBytesByRow;
	if (indexRowInBytes % (2*noPEs/8) == 0){
		fseek(inputFile, indexRow*noBytesByRow, SEEK_SET);
		fread(arrowRow,4,noIntRegsByRow,inputFile);
		return indexRow;
	}
	return 0;
}

int getLastDirection(int noPEs, int *arrowRow) {
	const int dirHcode = 1431655765, dirVcode = 2863311530;
	int sizeValid = 0, isRead = 0, isHV = 0, dirHV = 0;
	int noBytesByRow = 2*noPEs/8, dimUInt=sizeof(unsigned int);
	int noIntRegsByRow = noBytesByRow/dimUInt;

	for (int i=0;i<noIntRegsByRow;i++){
		if (arrowRow[i] == dirHcode) isHV++;
		if (arrowRow[i] == dirVcode) isHV--;
	}
	if (isHV == noIntRegsByRow) dirHV = 1;		//H
	if (isHV == -noIntRegsByRow) dirHV = -1;	//V
	
	if (dirHV == 0) 
		return 0;
	else
		return dirHV;
}

int getPosLastArrow(int noPEs, int *arrowRow) {
	int i=0, j=0, dimUInt=sizeof(unsigned int), noBitsUInt=8*dimUInt;
	int noBytesByRow = 2*noPEs/8, noUIntRegsByRow = noBytesByRow/dimUInt;
	int noIntRegsByRow = noBytesByRow/dimUInt;
	int dupla=0, noArrows=0, posLastArrow=-1;
	//ARROW = 0 @ ERROR, 1 @ UP, 2 @ LEFT, 3 @ DIAGONAL
	for (i=0; i<noUIntRegsByRow; i++){
		for (j=0; j<noBitsUInt/2; j++){
			dupla = (arrowRow[i]>>2*j) & 3;
			if (dupla != 0) {
				noArrows++;
				posLastArrow = i*8*dimUInt + 2*j;
			}
		}
	}
	return posLastArrow;
}

void swap(int *x, int *y) {
   int temp;
   temp = *x;    /* save the value at address x */
   *x = *y;      /* put y into x */
   *y = temp;    /* put temp into y */
   printf("temp = %d\n",temp);
     return;
}

int getPacketSeq(FILE *inputFile, char *packetSeq, int *lastUnreadSym, int posFirstEnterFile){
	int sizePacket = 2*dimPackSeq;
	static char bufferPacketSeq[2*dimPackSeq];
	//printf("sizePacket = %d\n",sizePacket);
	int indexPacketSeq=0;
	int indexStartReadPacket = *lastUnreadSym - sizePacket;
	if (indexStartReadPacket < posFirstEnterFile){
		indexStartReadPacket = posFirstEnterFile;
		sizePacket = *lastUnreadSym - posFirstEnterFile;
	}
	fseek(inputFile, indexStartReadPacket, SEEK_SET);
	fread(bufferPacketSeq,1,sizePacket,inputFile);
	int indexPacketFile=sizePacket-1;
	while (indexPacketSeq<dimPackSeq && indexPacketFile>0){
		if (bufferPacketSeq[indexPacketFile] != '\n'){
			packetSeq[indexPacketSeq] = bufferPacketSeq[indexPacketFile];
			indexPacketSeq++;
		}
		indexPacketFile--;
	}
	*lastUnreadSym = *lastUnreadSym - sizePacket + indexPacketFile;
	return indexPacketSeq;
}


//WHILE LOOP
int getArrow(int noPEs, int *arrowRow, int *posArrow, int *dirArrow, int *currentArrowRow, int posSeqA, int posSeqB, char tempSymA, char tempSymB, int lastARROW){
	int dimUInt=sizeof(unsigned int);
	int posInRow = *posArrow / (8*dimUInt);
	int posInReg = (*posArrow % (8*dimUInt));
	int ARROW = (arrowRow[posInRow]>>posInReg) & 3;
		//if (posSeqA<300) printf("%d",ARROW);
	int filterARROW = filterLastArrow(posSeqA, posSeqB, ARROW, tempSymA, tempSymB, lastARROW);
	decoArrow(posArrow, currentArrowRow, dirArrow, filterARROW);
	return filterARROW;
}

int filterLastArrow(int posSeqA, int posSeqB, int ARROW, char lastSymA, char lastSymB, int lastArrow) {
	if 		(posSeqA <= 0)	lastArrow = LEFT;
	else if (posSeqB <= 0)	lastArrow = UP;
	else if (lastSymA == lastSymB && lastArrow == DIAG)
		lastArrow = DIAG;
	else if (ARROW == 0)		{
		if (lastSymA == lastSymB)
			lastArrow = DIAG;
		else
			lastArrow = lastArrow;
	}
	else	lastArrow = ARROW;
	return lastArrow;
}

void decoArrow(int* posArrow, int* posRow, int* dirHV, int ARROW) {
	if (ARROW == 0) {	//Error
		*posRow = *posRow;//*offsetFila = 0;
		*posArrow = *posArrow;//*offsetPos = 0;
		printf("error!!! ARROW = 0\n");
	}
	else if (*dirHV == 1){ //V
		if (ARROW == UP) {	//Up
			*posRow = *posRow - 1;//*offsetFila = 1;
			*posArrow = *posArrow;//*offsetPos = 0;
			*dirHV = -1;
		}
		else if (ARROW == LEFT) {	//Left
			*posRow = *posRow - 1;//*offsetFila = 1;
			*posArrow = *posArrow +2;//*offsetPos = 2;
			*dirHV = -1;
		}
		else if (ARROW == DIAG) {	//Diagonal
			*posRow = *posRow - 2;//*offsetFila = 2;
			*posArrow = *posArrow;//*offsetPos = 0;
			*dirHV = *dirHV;
		}
	}
	else if (*dirHV == -1){ //H
		if (ARROW == UP) {	//Up
			*posRow = *posRow - 1;
			*posArrow = *posArrow - 2;//*offsetPos = -2;
			*dirHV = 1;
		}
		else if (ARROW == LEFT) {	//Left
			*posRow = *posRow - 1;//*offsetFila = 1;
			*posArrow = *posArrow;//*offsetPos = 0;
			*dirHV = 1;
		}
		else if (ARROW == DIAG) {	//Diagonal
			*posRow = *posRow - 2;//*offsetFila = 2;
			*posArrow = *posArrow;//*offsetPos = 0;
			*dirHV = *dirHV;
		}
	}
}

void GetAlignSymbols(int ARROW, char *symA, char *symB, char tempSymA, char tempSymB, int *posSeqA, int *posSeqB, int *posPacketA, int *posPacketB) {
	if ( ARROW == DIAG || ARROW == UP ) { 
		*symA = tempSymA;
		*posSeqA = *posSeqA - 1; 
		*posPacketA = *posPacketA - 1;
	}
	else *symA = '_';
	if ( ARROW == DIAG || ARROW == LEFT ) { 
		*symB = tempSymB;
		*posSeqB = *posSeqB - 1; 
		*posPacketB = *posPacketB - 1;
	}
	else *symB = '_';
}

void getSimilarityAndDistance(char symA, char symB, int *similarity, int *distance) {
	if (symA == symB && symA != '_') {
		*similarity = *similarity + 1;
	}else {
		*similarity = *similarity - 1;
		*distance = *distance + 1;
	}
}


/*
void updateSeqFile(int *cntSeq, int *cnt){
	if (*cntSeq == 0 && *cnt != 0){
		//printf("posSeqA**** = %d\n",posSeqA);
		for (i=0; i<dimPackSeq; i++){
			fwrite(&VfnewseqA[i], 1, 1, fnewseqA);
		}
	}else if (posSeqA == 0 && posSeqB == 0){
		//////Antes o despues del siguiente bloque
		if (gapA == 0) VfnewseqA[cntLocal] = VseqA[(posSeqA)%dimPackSeq];
			else VfnewseqA[cntLocal] = '_';
		for (i=0; i<=cntLocal; i++){
			fwrite(&VfnewseqA[i], 1, 1, fnewseqA);
		}
	}
}
*/