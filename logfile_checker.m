function logfile_checker(s_arr)
exptdir = '~/cm/Data/Functional/'
if isempty(s_arr)
    s_arr = [1 3:10 12:26]
end

for s = s_arr
    
    display(sprintf('Checking subject %02d...',s));
    [~, ~, idx] = CM_fMRIBehAnalysis(s);
    %check that we have the right number of trials
    
    checkBehAnalysisOutput(idx, s);
    
    
    mvpadir = fullfile(exptdir, sprintf('CM%03d',s),'mvpa');
    
    checkMvpaOnsets(mvpadir,idx,s);
    
    pause();
    
end
    
%     
%     sidx = find(s_arr == s);
%     time0 = 
%     assert(onsets(find(strcmp(names,'time0'))) == 400); %make sure there are the right number of trials
%     for r = 1:8
%         runix{r} = 
end

function checkBehAnalysisOutput(idx, s)
    assert(length(idx.time0) == 400);
    
    leftrun = (idx.hand == 'L');
    lefthand = idx.leftHand;
    righthand = idx.rightHand;
    
 
    
    imagesc([idx.ex idx.cm leftrun idx.wrongHand  idx.old, idx.indexFinger, idx.s_old, idx.hit, idx.new,idx.middleFinger, idx.s_new,(idx.indexFinger+idx.middleFinger)==0]);
    title(sprintf('s%02d Hands. [ex, cm, lr assigned, wrong hand used, OLD, index, old, hit, NEW, middle, new, nr]',s));
    
end

function checkMvpaOnsets(mvpadir,idx,s)
    auto_ons = load(fullfile(mvpadir,'mvpa_ons'));
    man_ons = load(fullfile(mvpadir, 'onsets_mvpa'));
    
    getons = @(x,y) x.onsets(find(strcmp(x.names,y)));
    
    assert(length(man_ons.onsets) == 17);
    
    
    auto_aamisscount = getons(auto_ons,'aaMiss');
    auto_aamisscount = length(auto_aamisscount{1});
    manu_aamisscount = length(man_ons.onsets{2}) + length(man_ons.onsets{10});
    display(sprintf('s%02d AA miss trials\nauto:%d\nmanu:%d \n',s,auto_aamisscount,manu_aamisscount));
    
    auto_eamisscount = getons(auto_ons,'eaMiss');
    auto_eamisscount = length(auto_eamisscount{1});
    manu_eamisscount = length(man_ons.onsets{4}) + length(man_ons.onsets{12});
    display(sprintf('s%02d EA miss trials\nauto:%d\nmanu:%d \n',s,auto_aamisscount,manu_aamisscount));
    
    auto_aafacount = getons(auto_ons,'aaFa');
    auto_aafacount = length(auto_aafacount{1});
    manu_aafacount = length(man_ons.onsets{6}) + length(man_ons.onsets{14});
    display(sprintf('s%02d AA fa trials\nauto:%d\nmanu:%d \n',s,auto_aafacount,manu_aafacount));
    
    auto_eafacount = getons(auto_ons,'eaFa');
    auto_eafacount = length(auto_eafacount{1});
    manu_eafacount = length(man_ons.onsets{8}) + length(man_ons.onsets{16});
    display(sprintf('s%02d EA fa trials\nauto:%d\nmanu:%d \n',s,auto_aafacount,manu_aafacount));
    
    
    
end
    
    