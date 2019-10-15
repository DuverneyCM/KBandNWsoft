#define NoPacketFile		200
#define dimPackSeqSocMax	512
#define MATCH				3
#define MISMATCH			-1
#define OGAP				-2
#define EGAP				-1
#define VALID				1
#define INVALID				0


getNoSoCPackets(int dimSeq, int dimFirstPacket, int dimLastPacket, int dimPacket){
	int noPackets=0;
	return noPackets;
}
void shiftRegisterRL(int dirHV, int noPEs, int *shRegA, int *shRegB, int codeSymA, int codeSymB){
	int i;
	if (dirHV == 0) {
		for (i=1; i<noPEs; i++){
			shRegA[i] = shRegA[i-1];
			shRegA[0] = codeSymA;
		}
	}
	else{
		for (i=0; i<noPEs-1; i++){
			shRegA[i] = shRegA[i+1];
			shRegA[noPEs-1] = codeSymB;
		}
	}
}
void getScoreAB(int noPEs, int *shRegA, int *shRegB, int *scoreAB, int *validAB){
	int i;
	for (i=0; i<noPEs; i++){
		if (shRegA[i]<1 || shRegA[i]>5 || shRegB[i]<1 || shRegB[i]>5){
			scoreAB[i] = 0;
			validAB[i] = INVALID;
		}
		else {
			validAB[i] = VALID;
			if (shRegA[i] == shRegB[i])
				scoreAB[i] = MATCH;
			else
				scoreAB[i] = MISMATCH;
		}
	}
}
void getScoreH(int dirHV, int noPEs, int *shRegA, int *shRegB, int *scoreAB, int *validAB, int *H1, int *H2, int *ARROW){
	int i;
	int hD, hU, hL, hInit, tempH, tempA;
	int hOut[noPEs];
	for (i=0; i<noPEs; i++){
	//Allias
		hD = H2[i] + scoreAB[i];
		hInit = H1[i];
		if (dirHV == 0){
			hL = H1[i];
			if (i>0) hU = H1[i-1];
			else hU = H1[i] - 1;
		}
		else{
			hU = H1[i];
			if (i<noPEs-1) hL = H1[i+1];
			else hL = H1[i] - 1;
		}
	//perform NWA
		if (hL-hU < 0){
			if (hD-hU < 0){
				tempH = hU;
				tempA = UP;
			}
			else if (hD == hU){
				tempH = hU;
				tempA = BOTH;
			}
			else{
				tempH = hD;
				tempA = DIAG;
			}
		}
		else{
			if (hD-hL < 0){
				tempH = hL;
				tempA = LEFT;
			}
			else if (hD == hL){
				tempH = hL;
				tempA = BOTH;
			}
			else{
				tempH = hD;
				tempA = DIAG;
			}
		}
	//Initialization Function
		if (validAB[i] == 1){
			hOut[i] = tempH + OGAP;
			ARROW[i] = tempA;
		}
		else{
			hOut[i] = hInit + OGAP;
			ARROW[i] = NOTARROW;
		}
	}
	//Update H values
	for (i=0; i<noPEs; i++){
		H2[i] = H1[i];
		H1[i] = hOut[i];
	}
}
void compressARROW(int noPEs, int *tempARROW, int *ARROW){
	int i;
	const int noArrowsByRow = 16;
	for (i=0; i<noPEs; i++){
		ARROW[i/noArrowsByRow] = (ARROW[i/noArrowsByRow]<<2) | tempARROW[i];
	}
}

runOneRowNWALinear(int *dirHV, int noPEs){
	int codeSymA, codeSymB;
	int shRegA[noPEs], shRegB[noPEs], scoreAB[noPEs], validAB[noPEs];
	int H1[noPEs], H2[noPEs], tempARROW[noPEs];
	int noRegs = 2*noPEs/(8*sizeof(unsigned int));
	int ARROW[noRegs];
	
	shiftRegisterRL(dirHV, noPEs, shRegA, shRegB, codeSymA, codeSymB);
	getScoreAB(noPEs, shRegA, shRegB, scoreAB, validAB);
	getScoreH(dirHV, noPEs, shRegA, shRegB, scoreAB, validAB, H1, H2, tempARROW);
	compressARROW (noPEs, tempARROW, ARROW);
	*dirHV = (*dirHV + 1) % 2;
}