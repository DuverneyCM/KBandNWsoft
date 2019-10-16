
#include <stdio.h>
#include <stdlib.h>
#include "tracebackFunctions.h"

int main(int argc, char *argv[]) {
//	Declaration of variables
	FILE *fseqA, *fseqB, *farrows;  //read, input files
	FILE *fnewseqA, *fnewseqB;      //write, output files
	int noPEs, noFile, noRegs, dimUInt;
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
//	Get features of the SeqA and SeqB
	int dimFileSeqA = getFileSize(fseqA);
	int dimFileSeqB = getFileSize(fseqB);
		//printf("dimFileSeq{A,B} = {%d,%d}\n",dimFileSeqA,dimFileSeqB);
	int	validFileA = isFasta(fseqA);
	int	validFileB = isFasta(fseqB);
		//printf("validFile{A,B} = {%d,%d}\n",validFileA,validFileB);
	int posFirstEnterFileA = getPosFirstEnter(fseqA, dimFileSeqA, searchWindow);
	int posFirstEnterFileB = getPosFirstEnter(fseqB, dimFileSeqB, searchWindow);
		//printf("posFirstEnterFile{A,B} = {%d,%d}\n",posFirstEnterFileA,posFirstEnterFileB);
	int dimLineSeqA = getLineSize(fseqA, dimFileSeqA, posFirstEnterFileA, searchWindow);
	int dimLineSeqB = getLineSize(fseqB, dimFileSeqB, posFirstEnterFileB, searchWindow);
		//printf("dimLineSeq{A,B} = {%d,%d}\n",dimLineSeqA,dimLineSeqB);
	int dimLastLineSeqA = getLastLineSize(dimFileSeqA, posFirstEnterFileA, dimLineSeqA);
	int dimLastLineSeqB = getLastLineSize(dimFileSeqB, posFirstEnterFileB, dimLineSeqB);
		//printf("dimLastLineSeq{A,B} = {%d,%d}\n",dimLastLineSeqA,dimLastLineSeqB);
	int dimSeqA = getDimSeq(dimFileSeqA, posFirstEnterFileA, dimLineSeqA);
	int dimSeqB = getDimSeq(dimFileSeqB, posFirstEnterFileB, dimLineSeqB);
		//printf("dimSeq{A,B} = {%d,%d}\n",dimSeqA,dimSeqB);
//	Read SeqA and SeqB by packets
	static char vSeqA[dimPackSeq+1], vNewSeqA[dimPackSeq+1];
	static char vSeqB[dimPackSeq+1], vNewSeqB[dimPackSeq+1];
	int lastUnreadSymFileSeqA = dimFileSeqA;
	int lastUnreadSymFileSeqB = dimFileSeqB;
	int dimPacketSeqA = getPacketSeq(fseqA, vSeqA, &lastUnreadSymFileSeqA, posFirstEnterFileA);
	int dimPacketSeqB = getPacketSeq(fseqB, vSeqB, &lastUnreadSymFileSeqB, posFirstEnterFileB);
	//swap(&lastUnreadSymFileSeqA, &lastUnreadSymFileSeqB);
	//printf("lastUnreadSymFileSeq{A,B} = {%d,%d}\n",lastUnreadSymFileSeqA, lastUnreadSymFileSeqB);
	//printf("dimPacketSeq{A,B} = {%d,%d}\n",dimPacketSeqA,dimPacketSeqB);
//	Get features of the Arrow files
	int noBytesByRow = 2*noPEs/8, arrowRow[noBytesByRow];
	sprintf(arrowFILEname, "%sp%d.bin", argv[4], noFile);
	farrows =	fopen(arrowFILEname,"rb");
		printf("File = %s\n",arrowFILEname);
	int dimFileArrows = getFileSize(farrows);
	int noArrowRows = getNoArrowRows(dimFileArrows, noPEs);
	int currentArrowRow = noArrowRows;
	//getArrowRowFromRAM();
	currentArrowRow = getArrowRowFromFILE(farrows, currentArrowRow-1, noPEs, arrowRow);
	int	dirArrow = getLastDirection(noPEs, arrowRow);
	currentArrowRow = getArrowRowFromFILE(farrows, currentArrowRow-1, noPEs, arrowRow);
	int posArrow = getPosLastArrow(noPEs, arrowRow);
		//printf("dimFileArrows = %d\n",dimFileArrows);
		//printf("noArrowRows = %d\n",noArrowRows);
		//printf("dirArrow = %d\n",dirArrow);
		//printf("posArrow = %d\n",posArrow);
//	While Loop: Traceback Process
	int similarity = 0, distance = 0, cntLocal = 0, cnt = 0;
	int posSeqA = dimSeqA, posSeqB = dimSeqB;
	int posPacketA = dimPacketSeqA, posPacketB = dimPacketSeqB;
	int gapA = 0, gapB = 0;
	int offsetPos, offsetRow;
	int ARROW = 0;
	char symA, symB, tempSymA, tempSymB;
	static char VfnewseqA[dimPackSeq+1], VfnewseqB[dimPackSeq+1];
	 
	//dimPacketSeqA, vSeqA, 
	//dirHV = dirArrow
	//offsetPos >> posArrow
	//offsetFila >> currentArrowRow
	while ( posSeqA > 0 || posSeqB > 0 ) {
		cntLocal = cnt % dimPackSeq;
		//load new packet if current packet of Seq or Arrows is already processed
		if (posPacketA == 0 && posSeqA>0){
			dimPacketSeqA = getPacketSeq(fseqA, vSeqA, &lastUnreadSymFileSeqA, posFirstEnterFileA);
			posPacketA = dimPacketSeqA;
		}
		if (posPacketB == 0 && posSeqB>0) {
			dimPacketSeqB = getPacketSeq(fseqB, vSeqB, &lastUnreadSymFileSeqB, posFirstEnterFileB);
			posPacketB = dimPacketSeqB;
		}
		if (currentArrowRow < 0){
			// Close current arrow file
			offsetRow = currentArrowRow;
			fclose(farrows);	noFile--;
			//	Open next arrow file
			sprintf(arrowFILEname, "%sp%d.bin", argv[4], noFile);
			farrows =	fopen(arrowFILEname,"rb");
				printf("File = %s\n",arrowFILEname);
			dimFileArrows = getFileSize(farrows);
			noArrowRows = getNoArrowRows(dimFileArrows, noPEs);
			//getArrowRowFromRAM();
			currentArrowRow = currentArrowRow + offsetRow;
		}
		
		tempSymA = vSeqA[dimPacketSeqA-posPacketA];
		tempSymB = vSeqB[dimPacketSeqB-posPacketB];
		//**	get arrow	**
		currentArrowRow = getArrowRowFromFILE(farrows, currentArrowRow, noPEs, arrowRow);
			ARROW = getArrow(noPEs, arrowRow, &posArrow, &dirArrow, &currentArrowRow, posSeqA, posSeqB, tempSymA, tempSymB, ARROW);
		//**	GetAlignSymbols	**
		GetAlignSymbols(ARROW, &symA, &symB, tempSymA, tempSymB, &posSeqA, &posSeqB, &posPacketA, &posPacketB);
		VfnewseqA[cntLocal] = symA;
		VfnewseqB[cntLocal] = symB;
		//**	GetSimilarityAndDistance	**
		getSimilarityAndDistance(symA, symB, &similarity, &distance);
		//**	Save new sequences in files	**
		if ( (cntLocal == 0 && cnt != 0) || (posSeqA == -1 && posSeqB == -1) ){
			for (i=0; i<dimPackSeq; i++){
				fwrite(&VfnewseqA[i], 1, 1, fnewseqA);
				fwrite(&VfnewseqB[i], 1, 1, fnewseqB);
			}
			printf("posSeqA**** = %d\n",posSeqA);
			printf("posSeqB**** = %d\n",posSeqB);
		}

		cnt++;
	}
	printf("length of new sequences = %d\n",cnt);
	printf("distance = %d, similarity = %d\n", distance, similarity);

	//close files 
	fclose(fseqA);	fclose(fseqB);	fclose(farrows);
	fclose(fnewseqA);	fclose(fnewseqB);
	return 0;
}