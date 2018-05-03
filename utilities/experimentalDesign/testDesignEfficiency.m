function e = testDesignEfficiency(source,TR,nEvents,glmType)
optimiseRepeats = 1000;
binSize = 1;
glmType = 'hrfModel';
 
for i = 1:optimiseRepeats
            sequence_Opti = zeros(nEvents,length(source));
            for n=1:nEvents
                seq = find(source==n);
                sequence_Opti(n,seq) = 1;
            end
       %     desMat_RandOp =sequence_RandOp;
    
    loopLength = nEvents/binSize;
    c = 1;
    if binSize ==1
        desMatBin_Opti = sequence_Opti;
    else
        for n = 1:loopLength
            desMatBin_Opti(n,:) = sum(sequence_Opti(c:c+binSize-1,:),1);
            c = c + binSize;
        end
        desMatBin_Opti(desMatBin_Opti>1) = 1;
    end
    
    switch(glmType)
        case('hrfModel')
            
            for n = 1:nEvents/binSize
                desMat_Opti(n,:) = convolveModelResponseWithHRF(desMatBin_Opti(n,:),hrf);
            end
        case('revCorr')
            desMat_Opti = desMatBin_Opti;
            for j=1:Nh-1
                desMat_Opti = [desMat_Opti;circshift(desMatBin_Opti,[0,j])];
            end
    end
    desMat_Opti_e = desMat_Opti-repmat(mean(desMat_Opti,2),1,length(desMat_Opti));
    e_RandOp(i) = 1/trace(inv(desMat_Opti_e * desMat_Opti_e'));
    
    SeqMeanRand = sequence_Opti-repmat(mean(sequence_Opti,2),1,length(sequence_Opti));
    eSeqRand(i) = 1/trace(inv(SeqMeanRand * SeqMeanRand'));
    %     imagesc(desMat_Opti')
    %     colorbar
    %     drawnow
    %     pause(0.5)
    
    %find shift value for max efficiency
    if e_RandOp(i)>eMax_Opti
        eMax_Opti = e_RandOp(i);
        eMax_opti_seq= eSeqRand(i);
        repeat = i;
        sequenceMax_Opti = sequence_Opti;
        DesignMatMax_Opti = desMat_Opti;
    end         

end
%% Voxel duty cycle
voxelDC = voxelDutyCycle(sequenceMax_Opti,nEvents,sigma);
% duty cycle of each condition weighted by a gaussian function with TW == pTW
% sigma = sig;
% mu = 16/2;
% x = 1:nEvents;
% gaus = 1 * exp(-(x - mu).^2/2/sigma^2);
% seqAv = mean(sequenceMax_Opti,2);
% voxelDC = mean(seqAv'.*gaus);
%% Sequence duty cycle
DutyCycle_Opti = mean(mean(sequenceMax_Opti));
figure('name',sprintf('Optimised Random sequence %s , nEvents = %i, Duty Cycle = %.4f',sequenceType,nEvents,DutyCycle_Opti));
subplot(2,2,1:2)
hist(e_RandOp);
title(sprintf('Efficiency: Design = %.4f,Sequence  = %.4f',eMax_Opti,eMax_opti_seq));
corrSeqMaxu = corr(sequenceMax_Opti');
corrDmaxu = corr(DesignMatMax_Opti');
subplot(2,2,3);
imagesc(corrSeqMaxu);
title(sprintf('N = %i',length(sequence_Opti)));
caxis([-1 1]);
colorbar;

