// PitchTrack.cpp: implementation of the CPitchTrack class.
//
//////////////////////////////////////////////////////////////////////

#include "PitchTrack.h"
#include <stdio.h>
#include <string.h>

CPitchTrack::CPitchTrack(int nFFTLen, int nHalfFFTLen,int fs)
				:m_nHalfFFTLen(nHalfFFTLen)
{   
	m_pfPitch = NULL;	m_pHarmonicsNum = NULL;
	m_nPfStartBin = 3;
	m_bValid =1;
	m_nfPitchStart = 100.0f;
	m_nfPitchEnd   = 400.0f;
	m_nfPitchStep  = 5.0f;
	m_nfPitchCnt   = int((m_nfPitchEnd-m_nfPitchStart)/m_nfPitchStep+1.0f);
	m_nMaxHarmonicsCnt = 15;
	m_nfTrackEnd   = 1250.0f;
	m_nInvBinFFT = 1.0f*nFFTLen/fs;

	int i;
	m_pfPitch			= new float[m_nfPitchCnt*2+m_nHalfFFTLen+m_nMaxHarmonicsCnt];
	if(!m_pfPitch) {m_bValid =0;	return;}
	m_pHarmonicsBuf		= m_pfPitch				+ m_nfPitchCnt;
	m_pHarmonicsFlagBuf = m_pHarmonicsBuf		+ m_nfPitchCnt;
	m_nfWeightArray		= m_pHarmonicsFlagBuf	+ m_nHalfFFTLen;
	m_nfWeightArray[0]= 1.f;
    m_nfWeightArray[1]= 0.915f;
    for(i=2;i<m_nMaxHarmonicsCnt;i++){
        m_nfWeightArray[i]=m_nfWeightArray[i-1]*m_nfWeightArray[1];
	}
    
	m_pHarmonicsNum = new int[m_nfPitchCnt];
	if(!m_pHarmonicsNum) {m_bValid =0;	return;}
	for( i=0;i<m_nfPitchCnt;i++){
		m_pfPitch[i]  = m_nfPitchStart+ i*m_nfPitchStep;
		m_pHarmonicsNum[i] = int(m_nfTrackEnd/(m_pfPitch[i]));	// get floor int
	}
}

CPitchTrack::~CPitchTrack()
{
	SAFE_DELETE_POINTER(m_pfPitch);
	SAFE_DELETE_POINTER(m_pHarmonicsNum);
}

void CPitchTrack::GetPitch(float *fAmp)
{
	int index,j,i;
	memset(m_pHarmonicsBuf,0,sizeof(float)*m_nfPitchCnt);
	for(i=0;i<m_nfPitchCnt;i++){
		for(j=1;j<=m_pHarmonicsNum[i];j++){
			index = int(m_pfPitch[i]*j*m_nInvBinFFT+0.5f);
			m_pHarmonicsBuf[i] += m_nfWeightArray[j-1]*fAmp[index];
		}
	}
	int max_bin;
	float maxHarmonics2 = V_GetMax(m_pHarmonicsBuf,max_bin,m_nfPitchCnt);
	m_nfPitchFre = m_pfPitch[max_bin];

	//m_nfMaxEnhBin =int(3.2*(m_pHarmonicsNum[max_bin]))+3;//When 16k,DSP process
	m_nfMaxEnhBin =int(3.2*(m_pHarmonicsNum[max_bin])+0.5f);//3.2*m_nfTrackEnd=4000Hz;
}

void CPitchTrack::GetPitchGain(float *fAmp,float *fModGain, float *fPwrRatio, float *fStateNoies)
{
	int j,i,index;
	GetPitch(fAmp);
	m_nEnhanceLowBin= int(600.0/m_nfPitchFre+0.5f);
    int BinUp = int(m_nfPitchFre*m_nfMaxEnhBin*m_nInvBinFFT+0.5f)+1;

	//m_pHarmonicsFlagBuf as tmp-buf
	V_Mult(fPwrRatio+m_nPfStartBin,fModGain+m_nPfStartBin,m_pHarmonicsFlagBuf+m_nPfStartBin,(BinUp-m_nPfStartBin));//
 	float nMinGain = V_GetMin(m_pHarmonicsFlagBuf+m_nPfStartBin,(BinUp-m_nPfStartBin));
	
	memset(m_pHarmonicsFlagBuf,0,sizeof(float)*m_nHalfFFTLen);
	for(j=1;j<=m_nEnhanceLowBin;j++){
		index = int(m_nfPitchFre*j*m_nInvBinFFT+0.5f);
		if((fStateNoies[index]<0.3f))//Do H
			fModGain[index] = 1;
		m_pHarmonicsFlagBuf[index-1] =  1;
		m_pHarmonicsFlagBuf[index]   =  1;
		m_pHarmonicsFlagBuf[index+1] =  1;
	}
	
	for(j=m_nEnhanceLowBin+1;j<=m_nfMaxEnhBin;j++)
	{
		index = int(m_nfPitchFre*j*m_nInvBinFFT+0.5f);
		if((fStateNoies[index]<0.3f)){//Do H
			fModGain[index] = fPwrRatio[index];
		}	
		m_pHarmonicsFlagBuf[index-1] =  1;
		m_pHarmonicsFlagBuf[index]   =  1;
		m_pHarmonicsFlagBuf[index+1] =  1;
	}
	
	for (i=m_nPfStartBin;i<BinUp;i++){
		fModGain[i] = ((1-m_pHarmonicsFlagBuf[i])* nMinGain + m_pHarmonicsFlagBuf[i]*fModGain[i]);
	}
}

//void CPitchTrack::GetPitchGain(float *fAmp,float *fModGain, float *fPwrRatio, float *fStateNoies)
//{
//	int j,i,index;
//	GetPitch(fAmp);
//	m_nEnhanceLowBin= int(600.0/m_nfPitchFre+0.5f);
//    int BinUp = int(m_nfPitchFre*m_nfMaxEnhBin*m_nInvBinFFT+0.5f)+1;
//
//	//m_pHarmonicsFlagBuf as tmp-buf
//	V_Mult(fPwrRatio+m_nPfStartBin,fModGain+m_nPfStartBin,m_pHarmonicsFlagBuf+m_nPfStartBin,(BinUp-m_nPfStartBin));//
// 	float nMinGain = V_GetMin(m_pHarmonicsFlagBuf+m_nPfStartBin,(BinUp-m_nPfStartBin));
//	
//	memset(m_pHarmonicsFlagBuf,0,sizeof(float)*m_nHalfFFTLen);
//	for(j=1;j<=m_nEnhanceLowBin;j++){
//		index = int(m_nfPitchFre*j*m_nInvBinFFT+0.5f);
//		if((fStateNoies[index]<0.3f))//Do H
//			fModGain[index] = 1;
//		m_pHarmonicsFlagBuf[index-1] =  1;
//		m_pHarmonicsFlagBuf[index]   =  1;
//		m_pHarmonicsFlagBuf[index+1] =  1;
//	}
//	
//	for(j=m_nEnhanceLowBin+1;j<=m_nfMaxEnhBin;j++)
//	{
//		index = int(m_nfPitchFre*j*m_nInvBinFFT+0.5f);
//		if((fStateNoies[index]<0.3f)){//Do H
//			fModGain[index] = fPwrRatio[index];
//		}	
//		m_pHarmonicsFlagBuf[index-1] =  1;
//		m_pHarmonicsFlagBuf[index]   =  1;
//		m_pHarmonicsFlagBuf[index+1] =  1;
//	}
//	
//	for (i=m_nPfStartBin;i<BinUp;i++){
//		fModGain[i] = ((1-m_pHarmonicsFlagBuf[i])* nMinGain + m_pHarmonicsFlagBuf[i]*fModGain[i]);
//	}
// }
