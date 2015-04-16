function CM_run_and_report_mvpa(subj_arr,runs_to_report,xvalIterToReport, varargin)
res = CM_run_mvpa_v4_tbm(subj_arr, '', varargin{:});
cd ('/Users/Jesse/fMRI/COUNTERMEASURES/Data/Functional/mvpa_results');
CM_mvpaPostProc(res, xvalIterToReport, runs_to_report,'ret', 1, pwd, subj_arr);
end
