
#include <stdio.h>
#include <stdlib.h>
#include "../traceback/tracebackFunctions.h"
#include "forwardFunctions.h"

int main(int argc, char *argv[]) {
//	Declaration of variables
	FILE *fseqA, *fseqB, *farrows;  //read, input files
	int noPEs, offsetBand, noFile=1, noRegs, dimUInt;
	int i = 0, j = 0; 				//Indice de los registros empezando por el último
	char arrowFILEname[32];
	int searchWindow = 200;
	dimUInt = sizeof(unsigned int);
	noRegs = 2*noPEs/(8*dimUInt);
//	Read Parameters. argv[0]	name of the executable file
	noPEs = atoi(argv[1]);				//Number of PEs. Usually 1024 PEs
	fseqA =	fopen(argv[2],"r");			//Name and path of the Sequence A
	fseqB =	fopen(argv[3],"r");			//Name and path of the Sequence B
	sprintf(arrowFILEname, "%sp%d.bin", argv[4], noFile);   //Name of the arrow file
	farrows =	fopen(arrowFILEname,"wb");  
		printf("File = %s\n",arrowFILEname);
		//Initialize first arrow file with TWO zero rows
		int dataArrows = 0;
		for (i=0; i<(2)*noRegs; i++) {	fwrite(&dataArrows, dimUInt, 1, farrows);	}	
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
//	Calculate Number of SoC Packets
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
	int dimPacketSeqA = getPacketSeq(fseqA, vSeqA, &lastUnreadSymFileSeqA, posFirstEnterFileA);
	int dimPacketSeqB = getPacketSeq(fseqB, vSeqB, &lastUnreadSymFileSeqB, posFirstEnterFileB);
		printf("lastUnreadSymFileSeq{A,B} = {%d,%d}\n",lastUnreadSymFileSeqA, lastUnreadSymFileSeqB);
		printf("dimPacketSeq{A,B} = {%d,%d}\n",dimPacketSeqA,dimPacketSeqB);

//	While Loop: Forward Process
	for (int pack=1; pack<=noSoCPacketsMax; pack++) {
	//load new packet if current packet of Seq is already processed
		if (pack % noSoCPacketsBySeqPackets == 0 && pack < noSoCPacketsA){
			dimPacketSeqA = getPacketSeq(fseqA, vSeqA, &lastUnreadSymFileSeqA, posFirstEnterFileA);
		}
		if (pack % noSoCPacketsBySeqPackets == 0 && pack < noSoCPacketsB) {
			dimPacketSeqB = getPacketSeq(fseqB, vSeqB, &lastUnreadSymFileSeqB, posFirstEnterFileB);
		}
	//Save arrows in a new file if the number of SoC packs reaches 200
		if (pack % 200 == 0){
			fclose(farrows);
			noFile++;
			sprintf(arrowFILEname, "%sp%d.bin", argv[4], noFile);
			printf("%s\n",arrowFILEname);
			farrows =	fopen(arrowFILEname,"wb");
		}
	//	Arrange SoC Packets
	//	Process Packet

	}
	fclose(fseqA);	fclose(fseqB);	fclose(farrows);
	//Memory Free
	return 0;
}