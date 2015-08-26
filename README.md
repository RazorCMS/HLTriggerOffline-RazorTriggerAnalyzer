HLTriggerOffline-RazorTriggerAnalyzer
=====================================
   * analyzer for trigger output
   * keep a version of HLT_ZeroBias (and ensure that it accepts *all* events) in order to compute trigger rates correctly

    cmsrel CMSSW_7_2_1_patch3
    cd CMSSW_7_2_1_patch3/src
    git cms-addpkg HLTriggerOffline/SUSYBSM
    git clone git@github.com:RazorCMS/HLTriggerOffline-RazorTriggerAnalyzer HLTriggerOffline/RazorTriggerAnalyzer
    scram b
   
