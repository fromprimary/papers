function [Fund_f pitch_index]=hamonic_iteration_find(y_segAmp,FrmNum,Fs,FFTLen,HalfFrmLen,VadFlag)
%%%%%%%%%%%%%%%%%%%%%%%%%
Vad_count=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%
y_segAmp1 = y_segAmp;
pitch_index=zeros(HalfFrmLen,FrmNum);
fstart=100;detf=5;fend=400;
f=fstart:detf:fend;  h=0.915;   f_Enh = Fs/2;   Nnh_Track=round(1250./f);  Nnh_Enh =floor(((Fs-200)/2)./f);
enhance_up_fre=4000;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%tracking
Tracking=1;
f_bin=round(f/31.25+1);
BandUp=max( f_bin);
Score=zeros(BandUp,FrmNum);
Fund_f=zeros(1,FrmNum);
argmaxBand=Score;
maxScore=zeros(BandUp,2);
unVoice=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ÇåÒôÅÐ¶Ï
Surd_process=0;
if Surd_process==1
    temp=sum(y_segAmp1(1:48,:).^2)./sum(y_segAmp1(49:end,:).^2);
    Surd=temp>0.95;
    clear temp
else %if Surd_process==1
    Surd=ones(FrmNum,1);
end %if Surd_process==1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for indx_frame=1:FrmNum
    if (VadFlag(indx_frame)==1)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        if Surd(indx_frame)==0&&Tracking==1
            unVoice=unVoice+1;
        else% if Surd(indx_frame)==0&&Tracking==1
            start=unVoice<6;
            unVoice=0;
        end% if Surd(indx_frame)==0&&Tracking==1
        %%%%%%%%%%%%%%%%%%%%%%%%%
        maxScore=zeros(BandUp,2);
        for ii=1:length(f)
            i_f=(round(f(ii)*[1:Nnh_Track(ii)]/Fs*FFTLen)+1);
            H(ii)=sum(h.^([0:Nnh_Track(ii)-1]').*y_segAmp1(i_f,indx_frame));
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% to find max-value in fre
            if H(ii)>maxScore(f_bin(ii),1)&&Tracking==1
                maxScore(f_bin(ii),1)=H(ii);
                maxScore(f_bin(ii),2)=f(ii);
            end% if H(ii)>maxScore(f_bin(ii),1)&&Tracking==1
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end% for ii=1:length(f)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if indx_frame==1&&Tracking==1&&unVoice==0
            Score(:,indx_frame)=(maxScore(:,1));
        elseif Tracking==1&&unVoice==0 %if indx_frame==1&&Tracking==1&&unVoice==0
            [tempV,tempB]=max(Score([1 2],indx_frame-1));%find max bin in previous frame
            Score(1,indx_frame)=maxScore(1,1)+start*tempV;%iteration
            argmaxBand(1,indx_frame)=tempB;%for back tracking
            ii=2;
            while ii<BandUp
                [tempV,tempB]=max(Score([ii-1:ii+1],indx_frame-1));%find max bin in previous frame
                Score(ii,indx_frame)=maxScore(ii,1)+start*tempV;%iteration
                argmaxBand(ii,indx_frame)=ii+tempB-2;
                ii=ii+1;
            end%if indx_frame==1&&Tracking==1&&unVoice==0
            [tempV,tempB]=max(Score([BandUp-1:BandUp],indx_frame-1));%find max bin in previous frame
            Score(BandUp,indx_frame)=maxScore(BandUp,1)+start*tempV;%iteration
            argmaxBand(BandUp,indx_frame)=BandUp+tempB-2;
        end%if indx_frame==1&&Tracking==1&&unVoice==0
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%to est pitch
        if Tracking==1
            if unVoice==0
                [tempV roughF]=max(Score(:,indx_frame));
                f_f=maxScore(roughF,2);%to est pitch
                if f_f>fstart&&f_f<fend
                    f_m=(f_f-fstart)/detf+1;
                    flag=1;
                else
                    flag=0;
                end%if f_f>fstart&&f_f<fend
                Fund_f(indx_frame)=f_f;
            else%if unVoice==0
                flag=0;
            end%if unVoice==0
        else %    if Tracking==1
            %%%%%%%%%%%%%%%%%%% for no tracking
            [amp_m,f_m]=max(H);
            f_f=f(f_m);
            Fund_f(indx_frame)=f_f;
            flag=1;
        end %    if Tracking==1
        if flag
            i_f_f=round(f_f*[1:Nnh_Enh(f_m)]/Fs*FFTLen)+1;
            pitch_index(i_f_f,indx_frame)=1;
            pitch_index(i_f_f+1,indx_frame)=0;
            pitch_index(i_f_f-1,indx_frame)=0;
        end % if flag
    else %if ((TdNlmsVadFlag_P(indx_frame)==0)&&(tbrr_Mean_buf(indx_frame)>tbrr_Mean_TH)&&(mean_gg_Flag(indx_frame)==0));
        unVoice=unVoice+1;%20091026
    end %if ((TdNlmsVadFlag_P(indx_frame)==0)&&(tbrr_Mean_buf(indx_frame)>tbrr_Mean_TH)&&(mean_gg_Flag(indx_frame)==0));
end %for indx_frame=1:FrmNum



