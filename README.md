# summer_student_project

## Recipe
cmsrel CMSSW_12_3_0

cd CMSSW_12_3_0/src/

cmsenv

git clone https://ahussein1194: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@github.com/ahussein1194/summerStudentProject_Ahmed.git

cd summerStudentProject_Ahmed/RPC2TMAna

scram b -j 8

cmsRun python/ConfFile_cfg.py
