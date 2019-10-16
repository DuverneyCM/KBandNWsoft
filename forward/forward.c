
#include <stdio.h>
#include <stdlib.h>
#include "../traceback/tracebackFunctions.h"
#include "forwardFunctions.h"

int main(int argc, char *argv[]) {
//	Declaration of variables
	FILE *fseqA, *fseqB, *farrows;  //read, input files
	int noPEs, offsetBand, noFile=1, noRegs32ByRow, dimUInt;
	int i = 0, j = 0; 				//Indice de los registros empezando por el último
	char arrowFILEname[32];
	int searchWindow = 200;
	dimUInt = sizeof(unsigned int);
	noRegs32ByRow = 2*noPEs/(8*dimUInt);
//	Read Parameters. argv[0]	name of the executable file
	noPEs = atoi(argv[1]);				//Number of PEs. Usually 1024 PEs
	fseqA =	fopen(argv[2],"r");			//Name and path of the Sequence A
	fseqB =	fopen(argv[3],"r");			//Name and path of the Sequence B
	sprintf(arrowFILEname, "%sp%d.bin", argv[4], noFile);   //Name of the arrow file
	farrows =	fopen(arrowFILEname,"wb");  
		printf("File = %s\n",arrowFILEname);
		//Initialize first arrow file with TWO zero rows
		int dataArrows = 0;
		for (i=0; i<(2)*noRegs32ByRow; i++) {	fwrite(&dataArrows, dimUInt, 1, farrows);	}	
	offsetBand =	atoi(argv[5]);			//Number of arrows files generated in FWP
	
//	Get features of the SeqA and SeqB
	int dimFileSeqA = getFileSize(fseqA);
	int dimFileSeqB = getFileSize(fseqB);
		printf("dimFileSeq{A,B} = {%d,%d}\n",dimFileSeqA,dimFileSeqB);
	int	validFileA = isFasta(fseqA);
	int	validFileB = isFasta(fseqB);
		printf("validFile{A,B} = {%d,%d}\n",validFileA,validFileB);
	int posFirstEnterFileA = getPosFirstEnter(fseqA, dimFileSeqA, searchWindow);
	int posFirstEnterFileB = getPosFirstEnter(fseqB, dimFileSeqB, searchWindow);
		printf("posFirstEnterFile{A,B} = {%d,%d}\n",posFirstEnterFileA,posFirstEnterFileB);
	int dimLineSeqA = getLineSize(fseqA, dimFileSeqA, posFirstEnterFileA, searchWindow);
	int dimLineSeqB = getLineSize(fseqB, dimFileSeqB, posFirstEnterFileB, searchWindow);
		printf("dimLineSeq{A,B} = {%d,%d}\n",dimLineSeqA,dimLineSeqB);
	int dimLastLineSeqA = getLastLineSize(dimFileSeqA, posFirstEnterFileA, dimLineSeqA);
	int dimLastLineSeqB = getLastLineSize(dimFileSeqB, posFirstEnterFileB, dimLineSeqB);
		printf("dimLastLineSeq{A,B} = {%d,%d}\n",dimLastLineSeqA,dimLastLineSeqB);
	int dimSeqA = getDimSeq(dimFileSeqA, posFirstEnterFileA, dimLineSeqA);
	int dimSeqB = getDimSeq(dimFileSeqB, posFirstEnterFileB, dimLineSeqB);
		printf("dimSeq{A,B} = {%d,%d}\n",dimSeqA,dimSeqB);
	int diffSeqAB = dimSeqB - dimSeqA;
		if (diffSeqAB >= noPEs) printf("ERROR: diffSeq{A,B} is greater than noPEs");
		offsetBand = offsetBand + diffSeqAB/2;
//	Calculate Number of SoC Packets (innecesary)
	int dimFirstPacketA, dimFirstPacketB, dimLastPacketA, dimLastPacketB;
	int noSoCPacketsA = getNoSoCPackets(dimSeqA,&dimFirstPacketA,&dimLastPacketA,dimPackSeqSocMax);
	int noSoCPacketsB = getNoSoCPackets(dimSeqA,&dimFirstPacketB,&dimLastPacketB,dimPackSeqSocMax);
	int noSoCPacketsMax = maxAB(noSoCPacketsA,noSoCPacketsB);
	int noSoCPacketsBySeqPackets = dimPackSeq/dimPackSeqSocMax;
//	Read SeqA and SeqB by megapackets
	static char vSeqA[dimPackSeq+1], vNewSeqA[dimPackSeq+1];
	static char vSeqB[dimPackSeq+1], vNewSeqB[dimPackSeq+1];
	int lastUnreadSymFileSeqA = dimFileSeqA;
	int lastUnreadSymFileSeqB = dimFileSeqB;
	int dimPacketSeqA, dimPacketSeqB;
	//First Packet
	if (dimSeqA >= dimSeqB){
		dimPacketSeqA = getPacketSeqSoC(fseqA, vSeqA, &lastUnreadSymFileSeqA, posFirstEnterFileA, 0);
		dimPacketSeqB = getPacketSeqSoC(fseqB, vSeqB, &lastUnreadSymFileSeqB, posFirstEnterFileB, offsetBand);
	}	
	else{
		dimPacketSeqA = getPacketSeqSoC(fseqA, vSeqA, &lastUnreadSymFileSeqA, posFirstEnterFileA, offsetBand);
		dimPacketSeqB = getPacketSeqSoC(fseqB, vSeqB, &lastUnreadSymFileSeqB, posFirstEnterFileB, 0);
	}
		printf("lastUnreadSymFileSeq{A,B} = {%d,%d}\n",lastUnreadSymFileSeqA, lastUnreadSymFileSeqB);
		printf("dimPacketSeq{A,B} = {%d,%d}\n",dimPacketSeqA,dimPacketSeqB);

//	While Loop: Forward Process
	int dirHV=0, indexDiagRow=0, pack=0;
	int shRegA[noPEs], shRegB[noPEs], H1[noPEs], H2[noPEs], ARROWrow[noPEs];
	while(lastUnreadSymFileSeqA>posFirstEnterFileA || lastUnreadSymFileSeqB>posFirstEnterFileB){
	//	load new packet if current packet of Seq is already processed	
		dimPacketSeqA = 0; dimPacketSeqB = 0;
		if (lastUnreadSymFileSeqA > posFirstEnterFileA)
			dimPacketSeqA = getPacketSeqSoC(fseqA, vSeqA, &lastUnreadSymFileSeqA, posFirstEnterFileA, 0);
		if (lastUnreadSymFileSeqB > posFirstEnterFileB)
			dimPacketSeqB = getPacketSeqSoC(fseqB, vSeqB, &lastUnreadSymFileSeqB, posFirstEnterFileB, 0);
		if (dimPacketSeqA < dimPacketSeqB)
			for (i=dimPacketSeqA; i<dimPacketSeqB; i++)		vSeqA[i] = 0;
		if (dimPacketSeqB < dimPacketSeqA)
			for (i=dimPacketSeqB; i<dimPacketSeqA; i++)		vSeqB[i] = 0;
	//	Save arrows in a new file if the number of SoC packs reaches 200
		int maxNoPacketsByFile = 400*dimPackSeqSocMax/noPEs; //200
		if (pack % maxNoPacketsByFile == 0){
			fclose(farrows);
			noFile++;
			sprintf(arrowFILEname, "%sp%d.bin", argv[4], noFile);
			printf("%s\n",arrowFILEname);
			farrows =	fopen(arrowFILEname,"wb");
		}
	//	Process Packet
		for (i=0; i<dimPacketSeqA+dimPacketSeqB-1; i++){
			runOneRowNWALinear(&dirHV, noPEs, vSeqA, vSeqB, shRegA, shRegB, H1, H2, ARROWrow);
			fwrite(&ARROWrow, sizeof(int), noRegs32ByRow, farrows);
			indexDiagRow++;
		}
		pack++;
	}
	fclose(fseqA);	fclose(fseqB);	fclose(farrows);
	//Memory Free
	return 0;
}