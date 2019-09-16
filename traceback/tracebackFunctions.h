#include <stdio.h>
#include <stdlib.h>

//***   Declaration's Functions
int getFileSize(FILE *inputFile);
int isFasta(FILE *inputFile);
int getPosFirstEnter(FILE *inputFile, int dimFile, int dimSearchWindow);
int getLineSize(FILE *inputFile, int dimFile, int posFirstEnter, int dimSearchWindow);
int getLastLineSize(int dimFile, int posFirstEnter, int dimLine);
int getDimSeq(int dimFile, int posFirstEnter, int dimLine);
int getNoArrowRows(int dimFile, int noPEs);
int getArrowRow(FILE *inputFile, int indexRow, int noPEs, int *arrowRow);
int getLastDirection(int noPEs, int *arrowRow);
int getPosLastArrow(int noPEs, int *arrowRow);
int getArrow(int noPEs, int *arrowRow, int posArrow, int dirArrow);

void decoArrow(int* offsetPos, int* offsetFila, int* dirHV, int ARROW, int* posSeqA, int* posSeqB, int* gapA, int* gapB);



//***   Definition's Functions
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

int getArrowRow(FILE *inputFile, int indexRow, int noPEs, int *arrowRow){
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

int getArrow(int noPEs, int *arrowRow, int posArrow, int dirArrow){
	int dimUInt=sizeof(unsigned int);
	int posInRow = posArrow / (8*dimUInt);
	int posInReg = (posArrow % (8*dimUInt));
	int ARROW = (arrowRow[posInRow]>>posInReg) & 3;
	return ARROW;
}



void decoArrow(int* offsetPos, int* offsetFila, int* dirHV, int ARROW, int* posSeqA, int* posSeqB, int* gapA, int* gapB) {
	if (*dirHV == 1){ //V
		if (ARROW == 0) {	//Error
			*offsetFila = 0;
			*offsetPos = 0;
			printf("posSeqA** = %d\n",*posSeqA);
			printf("posSeqB** = %d\n",*posSeqB);
			if (*posSeqA != 0) {
				*posSeqA = *posSeqA - 1;
				*gapA = 0;
			} else *gapA = 1;
			if (*posSeqB != 0) {
				*posSeqB = *posSeqB - 1;
				*gapB = 0;
			} else *gapB = 1;
		}
		else if (ARROW == 1) {	//Up
			*offsetFila = 1;
			*offsetPos = 0;
			*dirHV = -1;
			*posSeqA = *posSeqA - 1;
			*gapA = 0;
			*gapB = 1;
		}
		else if (ARROW == 2) {	//Left
			*offsetFila = 1;
			*offsetPos = 2;
			*dirHV = -1;
			*posSeqB = *posSeqB - 1;
			*gapB = 0;
			*gapA = 1;
		}
		else if (ARROW == 3) {	//Diagonal
			*offsetFila = 2;
			*offsetPos = 0;
			*posSeqA = *posSeqA - 1;
			*gapA = 0;
			*posSeqB = *posSeqB - 1;
			*gapB = 0;
		}
	}
	else if (*dirHV == -1){ //H
		if (ARROW == 0) {	//Error
			*offsetFila = 0;
			*offsetPos = 0;
			printf("posSeqA** = %d\n",*posSeqA);
			printf("posSeqB** = %d\n",*posSeqB);
			if (*posSeqA != 0) {
				*posSeqA = *posSeqA - 1;
				*gapA = 0;
			} else *gapA = 1;
			if (*posSeqB != 0) {
				*posSeqB = *posSeqB - 1;
				*gapB = 0;
			} else *gapB = 1;
		}
		if (ARROW == 1) {	//Up
			*offsetFila = 1;
			*offsetPos = -2;
			*dirHV = 1;
			*posSeqA = *posSeqA - 1;
			*gapA = 0;
			*gapB = 1;
		}
		if (ARROW == 2) {	//Left
			*offsetFila = 1;
			*offsetPos = 0;
			*dirHV = 1;
			*posSeqB = *posSeqB - 1;
			*gapB = 0;
			*gapA = 1;
		}
		if (ARROW == 3) {	//Diagonal
			*offsetFila = 2;
			*offsetPos = 0;
			*posSeqA = *posSeqA - 1;
			*gapA = 0;
			*posSeqB = *posSeqB - 1;
			*gapB = 0;
		}
	}
}