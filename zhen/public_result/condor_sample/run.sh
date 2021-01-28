source $HOME/.bashrc

cd QCDAnalysis/AnalysisStep
gridenv
cmsRun test/gg_H2ZZTo4L.py
cmsRun test/tqH_H2ZZTo4L.py

cp gg_HToZZTo4L_M125_13TeV.root tqH_HToZZTo4L_M125.root $HOME/public/
