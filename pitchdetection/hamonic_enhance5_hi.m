function [GainFunc Fund_f]=hamonic_enhance5_hi(GainFunc,y_segAmp,ratio,FrmNum,Fs,FFTLen,GminBufFunc,HalfFrmLen,TdNlmsVadFlag_P,Stat_q_Est_V_fram,mean_gg_Flag)
%在hamonic_enhance3上增加Gain shaping,使y_segAmp乘Gain为一指定单频的谱。 no good
%     可能会受到时域平滑的影响，试验一下：平滑和平滑的区别
%加入清音判别 20090618
%20090623 增加pitch tracking
%20090630 对4500hz以上的bin不做置最小处理
%replace NoiseShap by VadR and gain contral method is modified 20091015


RemoveNoise=0;  NR_more=0;LowfrequencyTone=0;




%%%%%%%%%%%%%%%%%%%%%%%%%%

Vad_count=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%


Gain_Enh = GainFunc;
y_segAmp1 = y_segAmp.*Gain_Enh;


fstart=100;detf=5;fend=400;
f=fstart:detf:fend;  h=0.915;   f_Enh = Fs/2;   Nnh_Track=round(1250./f);  Nnh_Enh =floor(((Fs-200)/2)./f);
enhance_up_fre=4000;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%tracking
Tracking=1;
%            f_bin(1,:)=length(f);
f_bin=round(f/31.25+1);
BandUp=max( f_bin);
Score=zeros(BandUp,FrmNum);
Fund_f=zeros(1,FrmNum);
argmaxBand=Score;
maxScore=zeros(BandUp,2);%(:,1)=max(H); (:,2)=agrmax(H)
unVoice=0;
Amp=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 清音判断
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
    if indx_frame==2339
        indx_frame;
    end
    if (TdNlmsVadFlag_P(indx_frame)==0)&&mean_gg_Flag(indx_frame)==0% modified for td+pf 20091013

        %         Vad_count=min(4,Vad_count+1);%20091026

        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        if Surd(indx_frame)==0&&Tracking==1
            unVoice=unVoice+1;
        else% if Surd(indx_frame)==0&&Tracking==1
            start=unVoice<6;
            unVoice=0;
        end% if Surd(indx_frame)==0&&Tracking==1
        %%%%%%%%%%%%%%%%%%%%%%%%%
        %        if ((TdNlmsVadFlag_P(indx_frame)==0)&&(tbrr_Mean_buf(indx_frame)>tbrr_Mean_TH))&&(~index1(indx_frame));
        maxScore=zeros(BandUp,2);
        for ii=1:length(f)
            i_f=(round(f(ii)*[1:Nnh_Track(ii)]/Fs*FFTLen)+1);
            H(ii)=sum(h.^([0:Nnh_Track(ii)-1]').*y_segAmp1(i_f,indx_frame));
            %             H2(ii)=sum(h.^([0:Nnh_Track(ii)-1]').*y_segAmp2(i_f,indx_frame));
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%if
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
            [tempV,tempB]=max(Score([1 2],indx_frame-1));
            Score(1,indx_frame)=maxScore(1,1)+start*tempV;
            argmaxBand(1,indx_frame)=tempB;
            ii=2;
            while ii<BandUp
                [tempV,tempB]=max(Score([ii-1:ii+1],indx_frame-1));
                Score(ii,indx_frame)=maxScore(ii,1)+start*tempV;
                argmaxBand(ii,indx_frame)=ii+tempB-2;
                ii=ii+1;
            end%if indx_frame==1&&Tracking==1&&unVoice==0
            [tempV,tempB]=max(Score([BandUp-1:BandUp],indx_frame-1));
            Score(BandUp,indx_frame)=maxScore(BandUp,1)+start*tempV;
            argmaxBand(BandUp,indx_frame)=BandUp+tempB-2;
        end%if indx_frame==1&&Tracking==1&&unVoice==0
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if Tracking==1
            if unVoice==0
                [tempV roughF]=max(Score(:,indx_frame));
                f_f=maxScore(roughF,2);
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
            [amp_m,f_m]=max(H);
            f_f=f(f_m);
            Fund_f(indx_frame)=f_f;
            flag=1;
        end %    if Tracking==1
        if flag

            i_f_f=round(f_f*[1:Nnh_Enh(f_m)]/Fs*FFTLen)+1;

        end % if flag
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        enhance_up_m= round(enhance_up_fre./f_f);
        enhance_low_m=round(600./f_f);
        if (Surd(indx_frame)==1)&&flag;

            Stat_q_temp=(Stat_q_Est_V_fram(:,indx_frame)>0.3);

            Gain_Enh(i_f_f(1:enhance_low_m),indx_frame)=(1-Stat_q_temp(i_f_f(1:enhance_low_m)));

            Gain_Enh(i_f_f(enhance_low_m+1:enhance_up_m),indx_frame)=(1-Stat_q_temp(i_f_f(enhance_low_m+1:enhance_up_m))).*ratio(i_f_f(enhance_low_m+1:enhance_up_m),indx_frame);

    %%%%%%%%%%%%%%%%%%%%old test2
            ind=ones(HalfFrmLen,1);
            ind(i_f_f)=0; ind(i_f_f+1)=0;ind(i_f_f-1)=0;  index_up=round(4000/Fs*FFTLen);
            ind(index_up:end)=0;
             ind=logical(ind);low=min(Gain_Enh(4:enhance_up_m,indx_frame).*ratio(4:enhance_up_m,indx_frame));
             Gain_Enh(ind,indx_frame) =max(low, GminBufFunc(ind,indx_frame)) ;
%               Gain_Enh(ind,indx_frame) =max(min(Gain_Enh(ind,indx_frame).*ratio(ind,indx_frame)), GminBufFunc(ind,indx_frame)) ;
              %%%%%%%%%%oldtest2
              %
            

          


        end %    if (Surd(indx_frame)==1)&&flag;

    else %if ((TdNlmsVadFlag_P(indx_frame)==0)&&(tbrr_Mean_buf(indx_frame)>tbrr_Mean_TH)&&(mean_gg_Flag(indx_frame)==0));
        unVoice=unVoice+1;%20091026
        Amp=0;
        %         Vad_count=max(Vad_count-2,0);
        if indx_frame>1   %%%%20090720 避免历史记录不连续
            Score(:,indx_frame)=Score(:,indx_frame-1);
        end

        Gain_Enh(22:end,indx_frame) =max(Gain_Enh(22:end,indx_frame).*ratio(22:end,indx_frame),GminBufFunc(22:end,indx_frame));%test2


    end %if ((TdNlmsVadFlag_P(indx_frame)==0)&&(tbrr_Mean_buf(indx_frame)>tbrr_Mean_TH)&&(mean_gg_Flag(indx_frame)==0));


end %for indx_frame=1:FrmNum

GainFunc = Gain_Enh;

