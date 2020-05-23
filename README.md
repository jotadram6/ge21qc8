# gemqc8
If you want to analyse QC8 real data:
```bash
scram p -n QC8Test CMSSW CMSSW_10_6_0
cd QC8Test/src
cmsenv
git clone https://github.com/giovanni-mocellin/gemqc8.git
mv gemqc8/* .
rm -rf gemqc8
scram b -j 4
cd Analysis/GEMQC8/test
python launcher_hot_dead_strips.py #run_number#
python launcher_certify_events.py #run_number#
python launcher_fast_efficiency.py #run_number#
python launcher_alignment.py #run_number#
python launcher_validation.py #run_number# noAlignment yesMasks
python launcher_validation.py #run_number# noAlignment noMasks
python launcher_validation.py #run_number# yesAlignment yesMasks
python launcher_validation.py #run_number# yesAlignment noMasks
```

If you want to run simulations, download the package following these instructions:
```bash
scram p -n QC8Test CMSSW CMSSW_10_6_0
cd QC8Test/src
cmsenv
git clone https://github.com/giovanni-mocellin/gemqc8.git
mv gemqc8/* .
rm -rf gemqc8
rm -rf EventFilter
scram b -j 4
cd Analysis/GEMQC8/test
python launcher_sim_fast_efficiency.py #run_number#
python launcher_sim.py #run_number#
```
