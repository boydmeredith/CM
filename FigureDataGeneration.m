
%generate data for across subject hit v cr classification masked by ex hit
%> cr and cm hit > cr effects inclusively masked at the 2nd level across
%the held-out 24 subjects.

for i =1:6
    display(['Now running across subject classification for TR ' num2str(i) '... ']);
    trW = [0 0 0 0 0 0];
    trW(i) = 1;
    CM_wrap_mvpa_across_subjects_masks('EXminus',[],[], 'trWeights',trW);
end
%results will likely require that I sort out where the condensed runs live
%in the resulting res results struct