//#include "EventFilter/GEMRawToDigi/plugins/GEMLocalModeDataSource.h"
#include "EventFilter/GEMRawToDigi/interface/AMC13Event.h"
#include "EventFilter/GEMRawToDigi/interface/VFATdata.h"
#include "Analysis/Unpacker/plugins/GEMLocalModeDataSource.h"
#include "DataFormats/FEDRawData/interface/FEDRawDataCollection.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"

//#include "EventFilter/GEMRawToDigi/interface/AMC13Event.h"
//#include "EventFilter/GEMRawToDigi/interface/VFATdata.h"

#include <fstream>
#include <algorithm>

GEMLocalModeDataSource::GEMLocalModeDataSource(const edm::ParameterSet & pset, edm::InputSourceDescription const &desc) :
edm::ProducerSourceFromFiles(pset,desc,true), // true - RealData
m_hasFerolHeader( pset.getUntrackedParameter<bool>("hasFerolHeader",false)),
m_fedid( pset.getUntrackedParameter<int>("fedId", 12345)),
m_filenames( pset.getUntrackedParameter<std::vector<std::string> >("fileNames" )),
m_fileindex(0),
m_runnumber( pset.getUntrackedParameter<int>("runNumber",-1)),
m_currenteventnumber(1),
m_processEvents(),
m_nProcessedEvents(0),
m_nGoodEvents(0),
m_goodEvents()
{
	produces<FEDRawDataCollection>("gemLocalModeDataSource");
	
	if (pset.exists("processEvents")) {
		m_processEvents= pset.getUntrackedParameter<std::vector<unsigned int> >("processEvents",std::vector<unsigned int>());
	}
	
	int containsListFile=0;
	std::cout << "there are " << m_filenames.size() << " files: ";
	for (unsigned int i=0; i<m_filenames.size(); i++) {
		std::cout << " " << m_filenames.at(i);
		if (m_filenames.at(i).find(".lst")!=std::string::npos)
		containsListFile=1;
	}
	std::cout << "\n";
	
	if (containsListFile) {
		std::cout << "input file definitions contains list file\n";
		std::vector<std::string> inpfnames (m_filenames);
		m_filenames.clear();
		for (unsigned int i=0; i<inpfnames.size(); i++) {
			std::string fn= inpfnames[i];
			if (fn.find(".lst")==std::string::npos) {
				m_filenames.push_back(fn);
				std::cout << " -- keeping " << fn << "\n";
			}
			else {
				std::cout << "loading list from file <" << fn << ">\n";
				std::ifstream fin(fn.c_str());
				std::string line;
				while (!fin.eof() && getline(fin,line)) {
					m_filenames.push_back(line);
					std::cout << " -- adding " << line << "\n";
				}
				fin.close();
			}
		}
		std::cout << std::endl;
	}
	
	IOOffset size = -1;
	StorageFactory::getToModify()->enableAccounting(true);
	
	for (unsigned int i=0; i<m_filenames.size(); i++) {
		std::string fname = m_filenames.at(i);
		bool exists = StorageFactory::get()->check(fname, &size);
		std::cout << "file " << fname << " size " << size << std::endl;
		if (!exists) {
			std::stringstream ss;
			ss << "file not found (" << fname << ")";
			throw cms::Exception(ss.str());
		}
	}
	
	std::string currentfilename = m_filenames[m_fileindex];
	std::cout << "examining " << currentfilename << std::endl;
	m_fileindex++;
	
	// open file stream
	storage = StorageFactory::get()->open(currentfilename);
	// (throw if storage is 0)
	if (!storage) {
		std::stringstream ss;
		ss << "GEMLocalModeDataSource: failed to open the file";
		throw cms::Exception(ss.str());
	}
	
	m_runnumber = m_fileindex; // dummy run number
}


GEMLocalModeDataSource::~GEMLocalModeDataSource()
{
	std::cout << "GEMLocalModeDataSource::~GEMLocalModeDataSource nGoodEvents=" << m_nGoodEvents << ", nProcessedEvents=" << m_nProcessedEvents << std::endl;
	std::cout << " their numbers (may be limited to 100)\n";
	for (unsigned int i=0; i<m_goodEvents.size(); i++) {
		std::cout << " , " << m_goodEvents[i];
		if (i>99) break;
	}
	std::cout << "\n";
}


void GEMLocalModeDataSource::debug_studyFile()
{
	std::cout << "\ndebug_studyFile" << std::endl;
	std::cout << "m_fileindex=" << m_fileindex << ", m_filenames.size=" << m_filenames.size() << std::endl;
	std::string currentfilename= m_filenames[m_fileindex-1];
	std::cout << "checkpoint A" << std::endl;
	std::string fname= currentfilename.substr(5,currentfilename.size());
	std::cout << "trying to open the file " << fname << std::endl;
	std::ifstream fin(fname.c_str(),std::ios::binary);
	if (!fin.is_open()) {
		std::cout << "failed" << std::endl;
		return;
	}
	
	const int bufsize=150;
	uint64_t buf[bufsize];
	if (m_hasFerolHeader && 1) {
		fin.read((char*)&buf,3*sizeof(uint64_t));
		
		for (uint64_t i=0; i<3*sizeof(uint64_t); i++) {
			char c= ((char*)&buf)[i];
			std::cout << "i=" << i << ", code=0x" << std::hex << buf[i] << std::dec << ", char=" << c << std::endl;
		}
	}
	fin.read((char*)&buf,bufsize*sizeof(uint64_t));
	
	if (0) {
		uint64_t n= fin.gcount();
		std::cout << "gcount=" << n << "\n";
		if (n>205) n=205;
		for (uint64_t i=0; i<n; i++) {
			// for some reason this if does not work
			if (i>200) { std::cout << "limiting to " << i << std::endl; break; }
			char c= ((char*)&buf)[i];
			std::cout << "i=" << i << ", code=0x" << std::hex << buf[i] << std::dec << ", char=" << c << std::endl;
		}
	}
	
	std::cout << "got words\n";
	gem::AMC13Event amc13Event;
	for (int i=0; i<bufsize-1; i++) {
		amc13Event.setCDFHeader(buf[i]);
		amc13Event.setAMC13Header(buf[i+1]);
		
		int ok=0;
		//if (buf[i] > ((uint64_t)(1)<<60)) std::cout << "b ";
		//else { std::cout << "g "; ok=1; }
		//if (amc13Event.get_cb5()==5) ok=1;
		std::cout << ((ok) ? "g " : "b ");
		std::cout << "0x" << std::hex << buf[i] << std::dec << "  ";
		/*
		 if (ok) {
		 std::cout << "cb5=" << amc13Event.get_cb5() << " cb0=" << amc13Event.get_cb0()
		 << " nAMC=" << amc13Event.get_nAMC();
		 }
		 */
		std::cout << std::endl;
	}
	
	fin.close();
	std::cout << "debug_studyFile done" << std::endl;
}

void GEMLocalModeDataSource::fillDescriptions(edm::ConfigurationDescriptions & descriptions)
{
	edm::ParameterSetDescription desc;
	desc.setComment("Reads GEM data saved in 'local mode' to produce FEDRawDataCollection.");
	desc.addUntracked<std::vector<std::string> >("fileNames")
	->setComment("Names of files to be processed.");
	desc.addUntracked<unsigned int>("skipEvents")
	->setComment("Number of events to skip before next processing.");
	desc.addUntracked<bool>("hasFerolHeader", false)
	->setComment("Whether additional header is present.");
	desc.addUntracked<int>("fedId", 12345)
	->setComment("FedID value to embed into events.");
	desc.addUntracked<int>("runNumber", -1)
	->setComment("Which runNumber to embed:\n -1 - get from filename,\n other - use this value.");
	desc.addUntracked<std::vector<unsigned int> >("processEvents", std::vector<unsigned int>());
	desc.addUntracked<std::vector<edm::LuminosityBlockID>>("firstLuminosityBlockForEachRun", {});
	//ProductSelectorRules::fillDescription(desc, "inputCommands");
	
	descriptions.addDefault(desc);
}


bool GEMLocalModeDataSource::setRunAndEventInfo(edm::EventID &id, edm::TimeValue_t &time, edm::EventAuxiliary::ExperimentType &)
{
	//std::cout << "\nsetRunAndEventInfo m_fileindex=" << m_fileindex << std::endl;
	if (0 && (m_currenteventnumber==1)) {
		debug_studyFile();
		return false;
	}
	
	if (storage->eof()) {
		storage->close();
		
		if (m_fileindex >= m_filenames.size()) {
			std::cout << "end of last file" << std::endl;
			return false;
		}
		std::string currentfilename = m_filenames[m_fileindex];
		std::cout << "processing " << currentfilename << std::endl;
		m_fileindex++;
		storage = StorageFactory::get()->open(currentfilename);
		if (!storage) {
			std::cout << "failed to open next file" << std::endl;
			return false;
		}
	}
	
	Storage & inpFile = *storage;
	
	// create product (raw data)
	buffers.reset( new FEDRawDataCollection );
	
	// assume 1 record is 1 event
	std::vector<uint64_t> buf;
	const int tmpBufSize=40;
	uint64_t tmpBuf[tmpBufSize];
	
	do {
		buf.clear();
		//std::cout << "GEMLocalModeDataSource::setRunAndEventInfo m_currenteventnumber=" << m_currenteventnumber << std::endl;
		
		int prn=0;
		
		if (1 && m_processEvents.size() &&
				(std::find(m_processEvents.begin(),m_processEvents.end(),m_currenteventnumber-1)!=m_processEvents.end())) {
			prn=1;
		}
		if (prn ==1)std::cout << "dd"<<std::endl;
		//std::cout << "the number= " << sizeof(uint64_t) << std::endl;
		
		if (m_hasFerolHeader) {
			if (inpFile.read((char*)tmpBuf,3*sizeof(uint64_t))!=3*sizeof(uint64_t)) {
				std::cout << "failed to read next FerolHeader" << std::endl;
				storage->close();
				return false;
			}
			if (0) {
				std::cout << "headerSkipped" << std::endl;
				for (int i=0; i<3; i++) std::cout << "h" << (i+1) << " 0x" << std::hex << tmpBuf[i] << std::dec << "\n";
			}
		}
		
		// get CDFHeader and AMC13Header
		int n=inpFile.read((char*)tmpBuf,2*sizeof(uint64_t)); // number of bytes to read
		if (n!=2*sizeof(uint64_t)) {
			std::cout << "failed to read next 2 words" << std::endl;
			storage->close();
			return false;
		}
		m_nProcessedEvents++;
		
		buf.push_back(tmpBuf[0]);
		buf.push_back(tmpBuf[1]);
		//std::cout << "0x" << std::hex << tmpBuf[0] << " 0x" << tmpBuf[1] << std::dec << std::endl;
		
		gem::AMC13Event amc13Event;
		amc13Event.setCDFHeader(tmpBuf[0]);
		amc13Event.setAMC13Header(tmpBuf[1]);
		
		/*if (uint8_t(tmpBufSize) < amc13Event.nAMC()) {
		 std::cout << "code correction is needed: nAMC=" << amc13Event.get_nAMC() << std::endl;
		 return false;
		 }*/
		
		//std::cout << "GEMLocalModeDataSource: cb5=" << amc13Event.get_cb5() << std::endl;
		//std::cout << "GEMLocalModeDataSource: cb0=" << amc13Event.get_cb0() << std::endl;
		
		/*if ((amc13Event.get_cb5()!=5) || (amc13Event.get_cb0()!=0)) {
		 std::cout << "data format error (cb5,cb0)" << std::endl;
		 return false;
		 }*/
		
		
		/*if (prn) {
		 std::cout << "cdfHeader sourceId=" << amc13Event.get_sourceId() << " bxId=" << amc13Event.get_bxId() << " lv1Id=" << amc13Event.get_lv1Id() << " eventType=" << amc13Event.get_eventType() << " cb5=" << amc13Event.get_cb5() << "\n";
		 std::cout << "amc13EventHeader cb0=" << amc13Event.get_cb0() << " nAMC=" << amc13Event.get_nAMC() << std::endl;
		 }*/
		
		// read AMC headers
		//std::cout << "GEMLocalModeDataSource: nAMC=" << amc13Event.get_nAMC() << std::endl;
		n = inpFile.read((char*)tmpBuf, amc13Event.nAMC()*sizeof(uint64_t));
		/*if ((uint32_t)(n)!=amc13Event.get_nAMC()*8) {
		 std::cout << " ERROR got " << n << " chk=" << amc13Event.get_nAMC()*8 << std::endl;
		 }*/
		for (uint8_t ii=0; ii<amc13Event.nAMC(); ii++) {
			amc13Event.addAMCheader(tmpBuf[ii]);
			buf.push_back(tmpBuf[ii]);
			//std::cout << "ii=" << ii << " cb0=" << amc13Event.get_cb0(ii) << " amcNr=" << amc13Event.get_amcNr(ii) << " dataSize=" << amc13Event.get_dataSize(ii) << std::endl;
		}
		
		// read AMC payloads
		for (uint8_t iamc = 0; iamc<amc13Event.nAMC(); ++iamc) {
			gem::AMCdata amcData;
			n = inpFile.read((char*)tmpBuf, 3*sizeof(uint64_t));
			for (int ii=0; ii<3; ii++) buf.push_back(tmpBuf[ii]);
			amcData.setAMCheader1(tmpBuf[0]);
			amcData.setAMCheader2(tmpBuf[1]);
			amcData.setGEMeventHeader(tmpBuf[2]);
			/*if (prn) {
			 std::cout << "iamc=" << iamc << " amcData.amcNr=" << (uint32_t)(amcData.amcNum()) << " amcData.davCnt=" << (uint32_t)(amcData.davCnt()) << std::endl;
			 std::cout << " --  " << amcData.getAMCheader1_str() << "\n -- " << amcData.getAMCheader2_str() << std::endl;
			 std::cout << " --  " << amcData.getGEMeventHeader_str() << std::endl;
			 }*/
			
			// read GEB
			for (uint8_t igeb=0; igeb<amcData.davCnt(); igeb++) {
				gem::GEBdata gebData;
				n = inpFile.read((char*)tmpBuf, sizeof(uint64_t));
				buf.push_back(tmpBuf[0]);
				gebData.setChamberHeader(tmpBuf[0]);
				
				/*if (prn) {
				 std::cout << "igeb=" << uint32_t(igeb) << " vfatWordCount=" << gebData.get_vfatWordCnt() << std::endl;
				 std::cout << " -- " << gebData.getChamberHeader_str() << std::endl;
				 }*/
				
				if (gebData.vfatWordCnt()%3!=0) {
					throw cms::Exception("gebData.vfatWordCnt()%3!=0");
				}
				
				// check if one buffer accommodates information
				if (gebData.vfatWordCnt() <= tmpBufSize) {
					// one buffer is enough
					n = inpFile.read((char*)tmpBuf, gebData.vfatWordCnt() * sizeof(uint64_t));
					for (int ii=0; ii < gebData.vfatWordCnt(); ii++) {
						buf.push_back(tmpBuf[ii]);
					}
					
					/*if (prn) {
					 std::cout << "gebData.inputID= " << gebData.get_inputID() << "\n";
					 if (gebData.get_inputID()!=0) std::cout << " not ZERO!!\n\n";
					 for (int ii=0; ii<gebData.vfatWordCnt(); ii+=3) {
					 gem::VFATdata vfd;
					 vfd.read_fw(tmpBuf[ii]);
					 vfd.read_sw(tmpBuf[ii+1]);
					 vfd.read_tw(tmpBuf[ii+2]);
					 std::cout << " ii/3=" << ii/3 << " " << " vfatPos=" << vfd.get_pos() << "\n";
					 }
					 }*/
				}
				else {
					// use buffer several times
					const int allowNumBufs= 32000/tmpBufSize; // 32k is firmware limit
					if (allowNumBufs * tmpBufSize<gebData.vfatWordCnt()) {
						std::cout << "update code: tmpBufSize=" << tmpBufSize << ", allowNumBufs=" << allowNumBufs << ", tmpBufSize*allowNumBufs=" << (tmpBufSize*allowNumBufs) << " (firmware limit), gebData.vfatWordCnt=" << gebData.vfatWordCnt() << std::endl;
						std::cout << "current file name is " << m_filenames[m_fileindex-1] << std::endl;
						return false;
					}
					
					int neededBufs= gebData.vfatWordCnt()/tmpBufSize;
					if (gebData.vfatWordCnt()%tmpBufSize>0) neededBufs++;
					if (neededBufs>allowNumBufs) {
						std::cout << "code error (neededBufs>allowNumBufs)\n";
						return false;
					}
					for (int iUse=0; iUse<neededBufs; iUse++) {
						uint32_t chunkSize= tmpBufSize;
						if (iUse*tmpBufSize+chunkSize > gebData.vfatWordCnt()) {
							chunkSize= gebData.vfatWordCnt() - iUse*tmpBufSize;
						}
						//std::cout << "iUse=" << iUse << ", chunkSize=" << chunkSize << "\n";
						n = inpFile.read((char*)tmpBuf, chunkSize * sizeof(uint64_t));
						for (unsigned int ii=0; ii < chunkSize; ii++) {
							buf.push_back(tmpBuf[ii]);
						}
						
						/* if (prn) {
						 std::cout << "gebData.inputID= " << gebData.get_inputID() << "\n";
						 if (gebData.get_inputID()!=0) std::cout << " not ZERO!!\n\n";
						 for (unsigned int ii=0; ii<chunkSize; ii+=3) {
						 gem::VFATdata vfd;
						 vfd.read_fw(tmpBuf[ii]);
						 vfd.read_sw(tmpBuf[ii+1]);
						 vfd.read_tw(tmpBuf[ii+2]);
						 std::cout << " ii/3=" << ii/3 << " " << " vfatPos=" << vfd.get_pos() << "\n";
						 }
						 }*/
					}
				}
				
				// read gebData trailer
				n = inpFile.read((char*)tmpBuf, sizeof(uint64_t));
				buf.push_back(tmpBuf[0]);
				
				// check
				gebData.setChamberTrailer(tmpBuf[0]);
				//if (prn) std::cout << " -- " << gebData.getChamberTrailer_str() << std::endl;
				if (gebData.vfatWordCntT() != gebData.vfatWordCnt()) {
					std::cout << "corrupt data? gebData vfatWordCnt does not match" << std::endl;
				}
				
			} // end of geb loop
			
			// read GEMeventTrailer and AMCTrailer
			n = inpFile.read((char*)tmpBuf, 2*sizeof(uint64_t));
			buf.push_back(tmpBuf[0]);
			buf.push_back(tmpBuf[1]);
			/*if (prn) { // check
			 amcData.setGEMeventTrailer(tmpBuf[0]);
			 amcData.setAMCTrailer(tmpBuf[1]);
			 std::cout << " -- " << amcData.getGEMeventTrailer_str() << "\n -- " << amcData.getAMCtrailer_str() << std::endl;
			 }*/
		} // end of amc loop
		
		// read AMC13trailer and CDFTrailer
		n = inpFile.read((char*)tmpBuf, 2*sizeof(uint64_t));
		buf.push_back(tmpBuf[0]);
		buf.push_back(tmpBuf[1]);
		/*if (prn) {
		 amc13Event.setAMC13Trailer(tmpBuf[0]);
		 amc13Event.setCDFTrailer(tmpBuf[1]);
		 std::cout << " -- " << amc13Event.getAMC13Trailer_str() << "\n"
		 << " -- " << amc13Event.getCDFTrailer_str() << std::endl;
		 }*/
		// end of amc13Event
		
		//std::cout << "GEMLocalModeDataSource got " << buf.size() << " words\n";
		if (buf.size()>12) {
			m_nGoodEvents++;
			if (m_goodEvents.size()<100)
			m_goodEvents.push_back(m_currenteventnumber-1);
		}
		
		m_currenteventnumber++;
		
		if (m_processEvents.size() &&
				(std::find(m_processEvents.begin(),m_processEvents.end(),m_currenteventnumber-1)!=m_processEvents.end())) {
			std::cout << "got it" << std::endl;
			break;
		}
	}
	while (m_processEvents.size());
	
	//
	// create FEDRawData
	
	auto rawData = std::make_unique<FEDRawData>(sizeof(uint64_t)*buf.size());
	unsigned char *dataptr = rawData->data();
	
	for (uint16_t i=0; i<buf.size(); i++) {
		((uint64_t*)dataptr)[i] = buf[i];
	}
	
	FEDRawData & fedRawData = buffers->FEDData( m_fedid );
	fedRawData = *rawData;
	
	// get real event number
	uint32_t realeventno = synchronizeEvents();
	int set_runnumber= (m_runnumber!=0) ? m_runnumber : id.run();
	id = edm::EventID(set_runnumber, id.luminosityBlock(), realeventno);
	return true;
}


void GEMLocalModeDataSource::produce(edm::Event &event) {
	//std::cout << "GEMLocalModeDataSource::produce" << std::endl;
	event.put(std::move(buffers), "gemLocalModeDataSource");
	buffers.reset();
}

uint32_t GEMLocalModeDataSource::synchronizeEvents() {
	//std::cout << "GEMLocalModeDataSource::synchronizeEvents" << std::endl;
	int32_t result= m_currenteventnumber -1;
	return(uint32_t) result;
}


#include "FWCore/Framework/interface/InputSourceMacros.h"
DEFINE_FWK_INPUT_SOURCE(GEMLocalModeDataSource);

