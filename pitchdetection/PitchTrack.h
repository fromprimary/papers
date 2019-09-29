#ifndef _PITCH_TRACK_H
#define _PITCH_TRACK_H

#include "..\\Common_DSP\\vectoroperation.h"

class CPitchTrack  
{
public:
	int m_bValid;
	int m_nPfStartBin;
private:
	int m_nHalfFFTLen;
    float m_nfPitchStart;
	float m_nfPitchEnd;
	float m_nfPitchStep;
	float m_nfTrackEnd;
	//int m_nfEnhEndBin; //max num of inter loop
	int m_nfPitchCnt;
	float m_nInvBinFFT;
	float* m_nfWeightArray;
	int m_nMaxHarmonicsCnt;

	int *m_pHarmonicsNum;		//m_nfPitchCnt
	float *m_pfPitch;			//m_nfPitchCnt
	float *m_pHarmonicsBuf;		//m_nfPitchCnt
	float *m_pHarmonicsFlagBuf;	//m_nHalfFFTLen
    int m_nEnhanceLowBin;

	float m_nfPitchFre; //out pitch
    int m_nfMaxEnhBin; // out max num of enhanced bin
	
public:	
	CPitchTrack(int m_nFFTLen=512, int m_nHalfFFTLen=257,int m_nProcessFs=16000);
	~CPitchTrack();
	void GetPitch(float *fAmp);
	void GetPitchGain(float *fAmp,float *fModGain, float *fPwrRatio, float *fStateNoies);
	float *GetHarmonicsBin(){return m_pHarmonicsFlagBuf;};
};

#endif
