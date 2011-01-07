#include "DQMOffline/PFTau/plugins/BenchmarkClient.h"

#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"


#include "DQMServices/Core/interface/DQMStore.h"
#include "DQMServices/Core/interface/MonitorElement.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

//
// -- Constructor
//
BenchmarkClient::BenchmarkClient(const edm::ParameterSet& parameterSet)  
{
  folderNames_ = parameterSet.getParameter< std::vector<std::string> >( "FolderNames" );
  histogramNames_  = parameterSet.getParameter< std::vector<std::string> >( "HistogramNames" );
}
//
// -- BeginJob
//
void BenchmarkClient::beginJob() {

  dqmStore_ = edm::Service<DQMStore>().operator->();
  for (std::vector<std::string>::const_iterator it = folderNames_.begin();
       it != folderNames_.end(); it++) {
    std::cout << (*it) << std::endl;
  } 
}
//
// -- EndJobBegin Run
// 
void BenchmarkClient::endRun(edm::Run const& run, edm::EventSetup const& eSetup) {
  doSummaries();
}
//
// -- EndJob
// 
void BenchmarkClient::endJob() {

}
//
// -- Create Summaries
//
void BenchmarkClient::doSummaries() {

  for (std::vector<std::string>::const_iterator ifolder = folderNames_.begin();
                                    ifolder != folderNames_.end(); ifolder++) {
    std::string path  = "ParticleFlow/"+(*ifolder);

    for (std::vector<std::string>::const_iterator ihist = histogramNames_.begin();
       ihist != histogramNames_.end(); ihist++) {
      std::string hname = (*ihist); 
      createResolutionPlots(path, hname);
    }
  }
}
//
// -- Create Resolution Plots
//
void BenchmarkClient::createResolutionPlots(std::string& folder, std::string& name) {     
  MonitorElement* me = dqmStore_->get(folder+"/"+name);
  if (!me) return;
  MonitorElement* me_avarage;
  MonitorElement* me_rms;
  MonitorElement* me_mean;
  MonitorElement* me_sigma;
  if ( (me->kind() == MonitorElement::DQM_KIND_TH2F) ||
       (me->kind() == MonitorElement::DQM_KIND_TH2S) ||
       (me->kind() == MonitorElement::DQM_KIND_TH2D) ) {
    TH2* th = me->getTH2F();
    size_t nbinx = me->getNbinsX();
    size_t nbiny = me->getNbinsY();
    
    float ymin = th->GetYaxis()->GetXmin();
    float ymax = th->GetYaxis()->GetXmax();
    std::string xtit = th->GetXaxis()->GetTitle();
    std::string ytit = th->GetYaxis()->GetTitle();
    float* xbins = new float[nbinx+1];
    for (size_t ix = 1; ix < nbinx+1; ++ix) {
       xbins[ix-1] = th->GetBinLowEdge(ix);
       if (ix == nbinx) xbins[ix] = th->GetXaxis()->GetBinUpEdge(ix);
    }    

    std::string tit_new;
    dqmStore_->setCurrentFolder(folder);
    MonitorElement* me_slice = dqmStore_->book1D("PFlowSlice","PFlowSlice",nbiny,ymin,ymax); 
    
    tit_new = ";"+xtit+";Avarage_"+ytit; 
    me_avarage = dqmStore_->book1D("avarage_"+name,tit_new, nbinx, xbins); 
    tit_new = ";"+xtit+";RMS_"+ytit; 
    me_rms     = dqmStore_->book1D("rms_"+name,tit_new, nbinx, xbins); 
    tit_new = ";"+xtit+";Mean_"+ytit; 
    me_mean    = dqmStore_->book1D("mean_"+name,tit_new, nbinx, xbins); 
    tit_new = ";"+xtit+";Sigma_"+ytit; 				 
    me_sigma   = dqmStore_->book1D("sigma_"+name,tit_new, nbinx, xbins); 
				 
    double  avarage, rms, mean, sigma;
    for (size_t ix = 1; ix < nbinx+1; ++ix) {
      me_slice->Reset();
      for (size_t iy = 1; iy < nbiny+1; ++iy) {
	me_slice->setBinContent(iy,th->GetBinContent(ix,iy)); 
      }
      getHistogramParameters(me_slice, avarage, rms, mean, sigma);
      //      std::cout << "pt " << th->GetBinCenter(ix) << "  : " 
      //                << avarage << " " << rms << " " << mean << " " << sigma << std::endl;
      me_avarage->setBinContent(ix,avarage);
      me_rms->setBinContent(ix,rms);
      me_mean->setBinContent(ix,mean);
      me_sigma->setBinContent(ix,sigma);
    }
    if (me_slice) dqmStore_->removeElement(me_slice->getName());
    delete [] xbins;
  }
}
//
// -- Get Histogram Parameters
//
void BenchmarkClient::getHistogramParameters(MonitorElement* me_slice, double& avarage, 
					     double& rms,double& mean, double& sigma) {
  avarage = 0.0;
  rms     = 0.0;
  mean    = 0.0;
  sigma   = 0.0;

  if (!me_slice) return;
  if  (me_slice->kind() == MonitorElement::DQM_KIND_TH1F) {
    avarage = me_slice->getMean();
    rms     = me_slice->getRMS();  
    TH1F* th_slice = me_slice->getTH1F();
    if (th_slice && th_slice->GetEntries() > 0) {
      th_slice->Fit( "gaus","Q0");
      TF1* gaus = th_slice->GetFunction( "gaus" );
      if (gaus) {
	sigma = gaus->GetParameter(2);
        mean  = gaus->GetParameter(1);
      }
    }
  }
}
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE (BenchmarkClient) ;
