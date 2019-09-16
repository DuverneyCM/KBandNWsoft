
#include <stdio.h>
#include <stdlib.h>
#include "tracebackFunctions.h"

#define dimPacket		512
#define NoPacketFile	200
#define NooRegs			64
#define dimPackSeq		10000000
#define dirHcode		1431655765 //"010101..."
#define dirVcode		2863311530 //"101010..."
#define CLOCKS_BY_SEC	800000

int main(int argc, char *argv[]) {
	FILE *fseqA, *fseqB, *farrows;  //read, input files
	FILE *fnewseqA, *fnewseqB;      //write, output files
	int noPEs, noFile, noRegs, dimUInt;
	int similitud = 0, distancia = 0, cntLocal = 0;
	int i = 0, j = 0; 				//Indice de los registros empezando por el último
	char arrowFILEname[32];
	int searchWindow = 200;


	//	Read Parameters. argv[0]	name of the executable file
	noPEs = atoi(argv[1]);				//Number of PEs. Usually 1024 PEs
	fseqA =	fopen(argv[2],"r");			//Name and path of the Sequence A
	fseqB =	fopen(argv[3],"r");			//Name and path of the Sequence B
										//argv[4] is used to make the name of the arrow file
	noFile =	atoi(argv[5]);			//Number of arrows files generated in FWP
	fnewseqA	= fopen(argv[6],"w+");	//Name and path of the new Sequence A (with gaps)
	fnewseqB	= fopen(argv[7],"w+");	//Name and path of the new Sequence B (with gaps)
	
	dimUInt = sizeof(unsigned int);
	noRegs = 2*noPEs/(8*dimUInt);

	//	Get features of the SeqA
	int dimFileSeqA = getFileSize(fseqA);
	int	validFileA = isFasta(fseqA);
	int posFirstEnterFileA = getPosFirstEnter(fseqA, dimFileSeqA, searchWindow);
	int dimLineSeqA = getLineSize(fseqA, dimFileSeqA, posFirstEnterFileA, searchWindow);
	int dimLastLineSeqA = getLastLineSize(dimFileSeqA, posFirstEnterFileA, dimLineSeqA);
	int dimSeqA = getDimSeq(dimFileSeqA, posFirstEnterFileA, dimLineSeqA);
		printf("dimFileSeqA = %d\n",dimFileSeqA);
		printf("validFileA = %d\n",validFileA);
		printf("posFirstEnterFileA = %d\n",posFirstEnterFileA);
		printf("dimLineSeqA = %d\n",dimLineSeqA);
		printf("dimLastLineSeqA = %d\n",dimLastLineSeqA);
		printf("dimSeqA = %d\n",dimSeqA);
	//	Get features of the  SeqB
	int dimFileSeqB = getFileSize(fseqB);
	int	validFileB = isFasta(fseqB);
	int posFirstEnterFileB = getPosFirstEnter(fseqB, dimFileSeqB, searchWindow);
	int dimLineSeqB = getLineSize(fseqB, dimFileSeqB, posFirstEnterFileB, searchWindow);
	int dimLastLineSeqB = getLastLineSize(dimFileSeqB, posFirstEnterFileB, dimLineSeqB);
	int dimSeqB = getDimSeq(dimFileSeqB, posFirstEnterFileB, dimLineSeqB);
		printf("dimFileSeqB = %d\n",dimFileSeqB);
		printf("validFileB = %d\n",validFileB);
		printf("posFirstEnterFileB = %d\n",posFirstEnterFileB);
		printf("dimLineSeqB = %d\n",dimLineSeqB);
		printf("dimLastLineSeqB = %d\n",dimLastLineSeqB);
		printf("dimSeqB = %d\n",dimSeqB);
	//	Get features of the Arrow files
	sprintf(arrowFILEname, "%sp%d.bin", argv[4], noFile);
	farrows =	fopen(arrowFILEname,"rb");
		printf("File = %s\n",arrowFILEname);
	int dimFileArrows = getFileSize(farrows);
	int	lastDir = getLastDirection(farrows, dimFileArrows, noPEs);
	int NoArrowRows = getNoArrowRows(dimFileArrows, noPEs);
		printf("dimFileArrows = %d\n",dimFileArrows);
		printf("lastDir = %d\n",lastDir);
		printf("NoArrowRows = %d\n",NoArrowRows);
	

/*
	char dataIn1[dimPacket + 1];
		unsigned int fila[NoRegs];
		static unsigned int bufferArrows[NooRegs*NoPacketFile*(2*dimPacket)];
		static char VseqA[dimPackSeq+1], VseqB[dimPackSeq+1];
		static char VfnewseqA[dimPackSeq+1], VfnewseqB[dimPackSeq+1];

    int currentFILA; //lastPosMemory, firstPosMemory,



	//////VERIFICAR QUE LA ÚLTIMA FILA CONTIENE LA DIRECCIÓN H/V DEL ÚLTIMO CICLO DEL KBAND
	sprintf(arrowFILEname, "%sp%d.bin", argv[4], noFile);
	farrows =	fopen(arrowFILEname,"rb");
		printf("File = %s\n",arrowFILEname);
		
	fseek(farrows, 0, SEEK_END);
	currentFILA = ftell(farrows)/(NoRegs*dimUInt); //pos next arrow
	fseek(farrows, 0, SEEK_SET);
	fread(&bufferArrows,dimUInt,currentFILA*NoRegs,farrows);
	currentFILA = currentFILA - 1; // next row
	
	//fseek(farrows, -1*NoRegs*dimUInt, SEEK_CUR);
	

	int isHV = 0;
	int dirHV = 0;
	//fread(&fila,dimUInt,NoRegs,farrows); //Lee la fila desde el archivo
	for (i=0; i<NoRegs; i++)
		fila[i] = bufferArrows[i + currentFILA*NoRegs];
	for (i=0;i<NoRegs;i++){
		//printf("%X\n",fila[i]);
		if (fila[i] == (unsigned int)dirH) isHV++;
		if (fila[i] == (unsigned int)dirV) isHV--;
	}
		//si dirHV sigue siendo cero, error en el archivo
	if (isHV == NoRegs) dirHV = 1;		// V
	if (isHV == -NoRegs) dirHV = -1;	// H
		printf("isHV = %d\n",dirHV);
		printf("dirHV = %d\n",dirHV);
*/
}