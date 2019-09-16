
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
//	Get features of the SeqB
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
	int noArrowRows = getNoArrowRows(dimFileArrows, noPEs);
	int currentArrowRow = noArrowRows;
	currentArrowRow = getArrowRow(farrows, currentArrowRow-1, noPEs, arrowRow);
	int	dirArrow = getLastDirection(noPEs, arrowRow);
	currentArrowRow = getArrowRow(farrows, currentArrowRow-1, noPEs, arrowRow);
	int posArrow = getPosLastArrow(noPEs, arrowRow);
	int ARROW = getArrow(noPEs, arrowRow, posArrow, dirArrow);
		printf("dimFileArrows = %d\n",dimFileArrows);
		printf("noArrowRows = %d\n",noArrowRows);
		printf("dirArrow = %d\n",dirArrow);
		printf("posArrow = %d\n",posArrow);
		printf("ARROW = %d\n",ARROW);
//

	//////BUSCAR LA POSICIÓN DE LA FLECHA EN LA PRIMERA FILA DE FLECHAS (next row). SOLO DEBE HABER UNA FLECHA
	/*
		
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
	



	int cnt = 0, flag = 0;
	while ( posSeqA != 0 || posSeqB != 0 ) {
		cntLocal = cnt%dimPackSeq;
				
		//////DECODIFICAR LA FLECHA
		if 		(posSeqA == 0)	decoArrow(&offsetPos, &offsetFila, &dirHV, 2, &posSeqA, &posSeqB, &gapA, &gapB);
		else if (posSeqB == 0)	decoArrow(&offsetPos, &offsetFila, &dirHV, 1, &posSeqA, &posSeqB, &gapA, &gapB);
		else	decoArrow(&offsetPos, &offsetFila, &dirHV, ARROW, &posSeqA, &posSeqB, &gapA, &gapB);
		
		
				
		//////GUARDAR LAS NUEVAS SECUENCIAS EN LOS ARCHIVOS
		if (cntLocal == 0 && cnt != 0){
			printf("posSeqA**** = %d\n",posSeqA);
			printf("posSeqB**** = %d\n",posSeqB);
			for (i=0; i<dimPackSeq; i++){
				fwrite(&VfnewseqA[i], 1, 1, fnewseqA);
				fwrite(&VfnewseqB[i], 1, 1, fnewseqB);
			}
		}else if (posSeqA == 0 && posSeqB == 0){
			//////Antes o despues del siguiente bloque
			if (gapA == 0) VfnewseqA[cntLocal] = VseqA[(posSeqA)%dimPackSeq];
				else VfnewseqA[cntLocal] = '_';
			if (gapB == 0) VfnewseqB[cntLocal] = VseqB[(posSeqB)%dimPackSeq];
				else VfnewseqB[cntLocal] = '_';

			for (i=0; i<=cntLocal; i++){
				fwrite(&VfnewseqA[i], 1, 1, fnewseqA);
				fwrite(&VfnewseqB[i], 1, 1, fnewseqB);
			}
		}
		

		//LEER SEQA Y SEQB POR PAQUETES DE TAMAÑO dimPackSeq
		if ((posSeqA) % dimPackSeq == 0) {
		//if ((dimSeqA - posSeqA) % dimPackSeq == 0) {
			if (packA != 0){
				fseek(fseqA, posFirstEnterA + 1 + dimPackSeq*(packA-1), SEEK_SET);
				fread(&VseqA,dimPackSeq,1,fseqA);
				printf("packA = %d\n",packA);
				//printf("range packA = %d %d\n",posFirstEnterA + 1 + dimPackSeq*(packA-1), (int)ftell(fseqA) );
				printf("range packA = %d %d\n", dimPackSeq*(packA-1), (int)ftell(fseqA) - (posFirstEnterA + 1) );
				packA--;
			}
		}

		if ((posSeqB) % dimPackSeq == 0) {
		//if ((dimSeqB - posSeqB) % dimPackSeq == 0) {
			if (packB != 0){
				fseek(fseqB, posFirstEnterB + 1 + dimPackSeq*(packB-1), SEEK_SET);
				fread(&VseqB,dimPackSeq,1,fseqB);
				printf("packB = %d\n",packB);
				//printf("range packB = %d %d\n",posFirstEnterB + 1 + dimPackSeq*(packB-1), (int)ftell(fseqB) );
				printf("range packB = %d %d\n",dimPackSeq*(packB-1), (int)ftell(fseqB) - (posFirstEnterB + 1) );
				packB--;
			}
		}
		
		
		//////Antes o despues del siguiente bloque
		//Agregar a la secuencia alineada el nuevo elemento: nucleotido o gap
		if (gapA == 0) VfnewseqA[cntLocal] = VseqA[(posSeqA)%dimPackSeq];
			else VfnewseqA[cntLocal] = '_';
		if (gapB == 0) VfnewseqB[cntLocal] = VseqB[(posSeqB)%dimPackSeq];
			else VfnewseqB[cntLocal] = '_';
		//calcular distancia/similitud segun alineamiento óptimo
		if (VfnewseqA[cntLocal] == VfnewseqB[cntLocal]) {
			if (VfnewseqA[cntLocal] != '_') similitud++;
		}
		else {
			similitud--;
			distancia++;
		}

		
		//////OBTENER PRÓXIMA FLECHA
		//currentFILA = ftell(farrows)/(NoRegs*dimUInt); //pos next arrow
		if (currentFILA < offsetFila){ //&& noFile != 1
			printf("cnt = %d\n",cnt);
			printf("posSeqA = %d\n",posSeqA);
			printf("posSeqB = %d\n",posSeqB);
			printf("ARROW = %d\n",ARROW);
			printf("posArrow = %d\n",posArrow);
			
			flag = 1;
			//////GUARDAR FILA A INICIAR EN EL PRÓXIMO ARCHIVO
				offsetFila = (offsetFila - currentFILA); //más constante
				printf("offsetFila = %d\n",offsetFila);
			//////CERRAR EL ARCHIVO ACTUAL Y ABRIR EL SIGUIENTE ARCHIVO
			fclose(farrows);
			noFile--;
			sprintf(arrowFILEname, "%sp%d.bin", argv[4], noFile);
			farrows =	fopen(arrowFILEname,"rb");
			printf("File = %s\n",arrowFILEname);

			//////MOVER EL CURSOR A LA POSICIÓN DE LA NUEVA FLECHA
				fseek(farrows, 0, SEEK_END);
				//currentFILA = ftell(farrows)/(NoRegs*dimUInt); //pos next arrow
				currentFILA = ftell(farrows)/(NoRegs*dimUInt); //pos next arrow
				fseek(farrows, 0, SEEK_SET);
				fread(&bufferArrows,dimUInt,currentFILA*NoRegs,farrows);
		}
		
		
		//fseek(farrows, -(1+offsetFila)*NoRegs*dimUInt, SEEK_CUR);
		
		currentFILA = currentFILA - offsetFila; //pos next arrow
		fseek(farrows, -(1+offsetFila)*NoRegs*dimUInt, SEEK_CUR);
		//fread(&fila,dimUInt,NoRegs,farrows); //Lee la fila desde el archivo
		for (i=0; i<NoRegs; i++)
			fila[i] = bufferArrows[i + currentFILA*NoRegs];
		
		posArrow = posArrow + offsetPos;
			i = (posArrow) / (8*dimUInt); //cociente, #reg de la fila
			j = (posArrow) % (8*dimUInt); //residuo, #pos en el reg
		ARROW = (fila[i]>>j) & 3; //ARROW = 0 @ ERROR, 1 @ UP, 2 @ LEFT, 3 @ DIAGONAL
			//printf("ARROW = %d, posArrow = %d, dirHV = %d\n",ARROW,posArrow,dirHV);
		
		
		cnt++;
	}

	printf("distancia = %d, similitud = %d\n", distancia,similitud);

	
	//cerrar archivos 
	fclose(fseqA);
	fclose(fseqB);
	fclose(farrows);
	fclose(fnewseqA);
	fclose(fnewseqB);
	fclose(fnewIndexA);
	fclose(fnewIndexB);
	
	clock_t endT = clock();
	double timeSpent = (double)(beginT - endT)/ CLOCKS_BY_SEC;
	printf("timeSpent = %f \n", timeSpent );
	
	return 0;
}








	char dataIn1[dimPacket + 1];
		unsigned int fila[NoRegs];
		static unsigned int bufferArrows[NooRegs*NoPacketFile*(2*dimPacket)];
		static char VseqA[dimPackSeq+1], VseqB[dimPackSeq+1];
		static char VfnewseqA[dimPackSeq+1], VfnewseqB[dimPackSeq+1];
    int currentFILA; //lastPosMemory, firstPosMemory,
*/
}