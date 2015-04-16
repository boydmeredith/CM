function [res] = CM_assessMotion(subs)

for snum = subs
    subjID = sprintf('CM%03d', snum);
    cd(sprintf('/Users/Jesse/fMRI/COUNTERMEASURES/Data/Functional/CM%03d/Results_model2/Test',snum));
    load(sprintf('CM%03d_TEST_MvmtParamsConcat_AND_sessFX.mat',snum));
    
%     res.x1Max(snum) = max(R(:,8));
%     res.x2Max(snum) = max(R(:,9));
%     res.x3Max(snum) = max(R(:,10));
    res.x1Max(snum) = max(R(:,11));
    res.x2Max(snum) = max(R(:,12));
    res.x3Max(snum) = max(R(:,13));
    
    res.x1Min(snum) = min(R(:,11));
    res.x2Min(snum) = min(R(:,12));
    res.x3Min(snum) = min(R(:,13));
    
    res.x1Range(snum) = res.x1Max(snum) -  res.x1Min(snum);
    res.x2Range(snum) = res.x2Max(snum) -  res.x2Min(snum);
    res.x3Range(snum) = res.x3Max(snum) -  res.x3Min(snum);
    
    thisMaxSessJump = zeros(3,1);
    for rnum = 1:8
        idx = 1+256*(rnum-1):256*rnum;
        thisX1 = R(idx,11);
        thisX2 = R(idx,12);
        thisX3 = R(idx,13);
        diffX1 = thisX1 - [0; thisX1(1:end-1)];
        diffX2 = thisX2 - [0; thisX2(1:end-1)];
        diffX3 = thisX3 - [0; thisX3(1:end-1)];
        maxDiffX1 = max(diffX1);
        maxDiffX2 = max(diffX2);
        maxDiffX3 = max(diffX3);
        if maxDiffX1 > thisMaxSessJump(1)
            thisMaxSessJump(1) = maxDiffX1;
        end
        if maxDiffX2 > thisMaxSessJump(2)
            thisMaxSessJump(2) = maxDiffX2;
        end
        if maxDiffX3 > thisMaxSessJump(3)
            thisMaxSessJump(3) = maxDiffX3;
        end
        medianX3JumpByRun(rnum, snum) = median(diffX3);
        maxX3JumpByRun(rnum,snum) = maxDiffX3;
        
        medianX1JumpByRun(rnum, snum) = median(diffX1);
        medianX2JumpByRun(rnum, snum) = median(diffX2);

            
    end
    maxSessJump(:,snum) = thisMaxSessJump;
    
    
    
end
    close all;

    figure(1)
    maxes = [res.x1Max; res.x2Max; res.x3Max];%; res.x4Max; res.x5Max; res.x6Max];
    xticklabels = 1:26;
    xticks = linspace(1,size(maxes,2),numel(xticklabels));
    imagesc(maxes);
    set(gca, 'XTick', xticks, 'XTickLabel', xticklabels);
    title('Maxes');
    colorbar;
    
    figure(2)
    diffMaxes = maxSessJump;
    imagesc(diffMaxes);
    set(gca, 'XTick', xticks, 'XTickLabel', xticklabels);
    title('diff Maxes');
    colorbar;
    
    figure(3)
    imagesc(medianX3JumpByRun);
    set(gca, 'XTick', xticks, 'XTickLabel', xticklabels);
    title('median x3 Jump By Run');
    colorbar;
    
        
    figure(4)
    imagesc(medianX2JumpByRun);
    set(gca, 'XTick', xticks, 'XTickLabel', xticklabels);
    title('median x2 Jump By Run');
    colorbar;
        
    figure(5)
    imagesc(medianX1JumpByRun);
    set(gca, 'XTick', xticks, 'XTickLabel', xticklabels);
    title('median x1 Jump By Run');
    colorbar;
    
    
    
     figure(6)
    imagesc(maxX3JumpByRun);
    set(gca, 'XTick', xticks, 'XTickLabel', xticklabels);
    title('max X3 Jump By Run');
    colorbar;
    
    
%     figure(2)
%     ranges = [res.x1Range; res.x2Range; res.x3Range];%; res.x4Max; res.x5Max; res.x6Max];
%     imagesc(ranges);
%     set(gca, 'XTick', xticks, 'XTickLabel', xticklabels);
%     title('Ranges');
%     colorbar;
    