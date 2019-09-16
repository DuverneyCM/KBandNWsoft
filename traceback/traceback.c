
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
	int noBytesByRow = 2*noPEs/8, arrowRow[noBytesByRow];
	sprintf(arrowFILEname, "%sp%d.bin", argv[4], noFile);
	farrows =	fopen(arrowFILEname,"rb");
		printf("File = %s\n",arrowFILEname);
	int dimFileArrows = getFileSize(farrows);
	int NoArrowRows = getNoArrowRows(dimFileArrows, noPEs);
	int	dirArrow = getLastDirection(farrows, dimFileArrows, noPEs, arrowRow);
	int posArrow = getPosLastArrow(farrows, dimFileArrows, noPEs, arrowRow);
	int ARROW = getArrow(noPEs, arrowRow, posArrow, dirArrow);
		printf("dimFileArrows = %d\n",dimFileArrows);
		printf("NoArrowRows = %d\n",NoArrowRows);
		printf("dirArrow = %d\n",dirArrow);
		printf("posArrow = %d\n",posArrow);
		printf("ARROW = %d\n",ARROW);
/*	
	//////BUSCAR LA POSICIÓN DE LA FLECHA EN LA PRIMERA FILA DE FLECHAS (next row). SOLO DEBE HABER UNA FLECHA
	currentFILA = currentFILA - 1; // next row
	for (i=0; i<NoRegs; i++)
		fila[i] = bufferArrows[i + currentFILA*NoRegs];
	//fseek(farrows, -2*NoRegs*dimUInt, SEEK_CUR);
	//fread(&fila,dimUInt,NoRegs,farrows); //Lee la fila desde el archivo
	
	int dupla = 0;
	int noArrows = 0;
	int posArrow = -1;
	int ARROW = 0;
	//int prueba = (dirH>>6) & 3;
	//	printf("%8X\n",prueba);
	for (i=0;i<NoRegs;i++){
		for (j=0;j<dimUInt*8/2;j++){
			dupla = (fila[i]>>2*j) & 3;
			if (dupla != 0) {
				noArrows++;
				posArrow = i*8*dimUInt + 2*j;
				ARROW = dupla; //ARROW = 0 @ ERROR, 1 @ UP, 2 @ LEFT, 3 @ DIAGONAL
			}
		}
	}
		//Si noArrows es cero, no hay flechas en la fila.
		//Si es mayor a 1, hay varias flechas en la fila, lo cual no corresponde con la última fila de flechas
		printf("ARROW = %d\n",ARROW);
		printf("noArrows = %d\n",noArrows);
		printf("posArrow = %d\n",posArrow);
		
	//definir variables de destino
		int posSeqA = dimSeqA;
		int posSeqB = dimSeqB;
		int gapA = 0;
		int gapB = 0;
		int offsetPos;
		int offsetFila;
	
	int packA, packB;
		if (dimSeqA % dimPackSeq == 0) packA = dimSeqA / dimPackSeq;
		else packA = dimSeqA / dimPackSeq + 1;
		if (dimSeqB % dimPackSeq == 0) packB = dimSeqB / dimPackSeq;
		else packB = dimSeqB / dimPackSeq + 1;
		
		//LEER SEQA Y SEQB POR PAQUETES DE TAMAÑO dimPackSeq. Primer paquete (este) es de residuo
			fseek(fseqA, posFirstEnterA + 1 + dimPackSeq*(packA-1), SEEK_SET);
			fread(&VseqA,dimSeqA%dimPackSeq,1,fseqA);
			printf("packA = %d\n",packA);
			printf("range packA = %d %d\n", dimPackSeq*(packA-1), (int)ftell(fseqA) - (posFirstEnterA + 1) );
			packA--;
			
			fseek(fseqB, posFirstEnterB + 1 + dimPackSeq*(packB-1), SEEK_SET);
			fread(&VseqB,dimSeqB%dimPackSeq,1,fseqB);
			printf("packB = %d\n",packB);
			printf("range packB = %d %d\n",dimPackSeq*(packB-1), (int)ftell(fseqB) - (posFirstEnterB + 1) );
			packB--;
		
	int cnt = 0;
	int flag = 0;


	char dataIn1[dimPacket + 1];
		unsigned int fila[NoRegs];
		static unsigned int bufferArrows[NooRegs*NoPacketFile*(2*dimPacket)];
		static char VseqA[dimPackSeq+1], VseqB[dimPackSeq+1];
		static char VfnewseqA[dimPackSeq+1], VfnewseqB[dimPackSeq+1];
    int currentFILA; //lastPosMemory, firstPosMemory,
*/
}