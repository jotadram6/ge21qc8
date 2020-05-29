iSpy CMSSW Analyzers for QC8

The Part of CMS iSpy code is modified and implemented in here to make QC8 event display. Please be aware this is originally created and developed by cms-outreach team. Original package can be found in http://cms-outreach.github.io/ispy/

To use:

* Clone the necessary source code (under QC8Test/src):

```
    git clone https://github.com/cms-outreach/ispy-services.git ISpy/Services
```

* Compile the code with the command:

```
    scram b
```

* Once compiled, change to ISpy/Analyzers:

```
    cd ISpy/Analyzers/python
```

* Run configuration file to make Geometry display:

```
    cmsRun ispy_qc8_geometry.py
```

* Run configuration file to make Event display:

```
    cmsRun ispy_qc8.py
```

Dataformat for ispy is '.ig' file. Each run will return below file,  

qc8_$num.ig: Event display with Track, Hit, Detector with hit.  

qc8-geometry.ig: QC8 Geometry without any event information.


* Load the ouput file in http://cern.ch/ispy-webgl  

But GEM dectector geometry is not currently readable in official page. To read QC8 Geometry, please download https://github.com/tt8888tt/ispy-webgl in local and open 'index.html' and load the file.

