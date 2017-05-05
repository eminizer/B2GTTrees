#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Common/interface/Handle.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "TH1.h"
#include "TH2.h"
#include "TLorentzVector.h"
#include <iostream>

class  EventCounter : public edm::EDAnalyzer {
public:
  EventCounter(const edm::ParameterSet& iConfig)
  {
    edm::Service<TFileService> fs;
    h_pileup_ = fs->make<TH1D>("pileup",";Pile-up", 100,0,100);
    h_weightsign_ = fs->make<TH1D>("weightsign",";Weight sign", 2, -1,1);
    h_totweight_ = fs->make<TH1D>("totweight",";Sum of Weights", 1,0,1);

    edm::EDGetTokenT<int>(mayConsume<int>(edm::InputTag("eventUserData", "puNtrueInt")));
    edm::EDGetTokenT<float>(mayConsume<float>(edm::InputTag("extraVar", "evtGenWeight")));

    h_cstar_ = fs->make<TH2D>("cstar",";Event type (0=qqbar, 1=qg, 2=gg);cstar", 3,0.,3.,100,-1.0,1.0);
    h_x_F_ = fs->make<TH2D>("x_F",";Event type (0=qqbar, 1=qg, 2=gg);x_F", 3,0.,3.,100,-1.0,1.0);
    h_M_ = fs->make<TH2D>("M",";Event type (0=qqbar, 1=qg, 2=gg);M", 3,0.,3.,500,0.0,5000.0);
    h_initial_parton_IDs_ = fs->make<TH2D>("h_initial_parton_IDs",";Parton 1 ID;Parton 2 ID", 44,-22,22,44,-22,22);

    edm::EDGetTokenT<float>(mayConsume<float>(edm::InputTag("extraVar", "MCcstar")));
    edm::EDGetTokenT<float>(mayConsume<float>(edm::InputTag("extraVar", "MCxF")));
    edm::EDGetTokenT<float>(mayConsume<float>(edm::InputTag("extraVar", "MCMtt")));
    edm::EDGetTokenT<float>(mayConsume<float>(edm::InputTag("extraVar", "MCpart1ID")));
    edm::EDGetTokenT<float>(mayConsume<float>(edm::InputTag("extraVar", "MCpart2ID")));

    h_cstar_vs_beta_qqbar_ = fs->make<TH2D>("cstar_vs_beta_qqbar","c* vs #beta (q#bar{q} events); #beta; c*", 10,0.,1.,20,-1.,1.);
    h_cstar_vs_beta_gg_    = fs->make<TH2D>("cstar_vs_beta_gg","c* vs #beta (qg/gg events); #beta; c*", 10,0.,1.,20,-1.,1.);

    edm::EDGetTokenT<float>(mayConsume<float>(edm::InputTag("extraVar", "MCtpt")));
    edm::EDGetTokenT<float>(mayConsume<float>(edm::InputTag("extraVar", "MCteta")));
    edm::EDGetTokenT<float>(mayConsume<float>(edm::InputTag("extraVar", "MCtphi")));
    edm::EDGetTokenT<float>(mayConsume<float>(edm::InputTag("extraVar", "MCtE")));
    edm::EDGetTokenT<float>(mayConsume<float>(edm::InputTag("extraVar", "MCtbarpt")));
    edm::EDGetTokenT<float>(mayConsume<float>(edm::InputTag("extraVar", "MCtbareta")));
    edm::EDGetTokenT<float>(mayConsume<float>(edm::InputTag("extraVar", "MCtbarphi")));
    edm::EDGetTokenT<float>(mayConsume<float>(edm::InputTag("extraVar", "MCtbarE")));

    h_mu_R_sf_up_   = fs->make<TH1D>("mu_R_sf_up",  ";mu_R_sf_up",   100,-1.,3.);
    h_mu_R_sf_down_ = fs->make<TH1D>("mu_R_sf_down",";mu_R_sf_down", 100,-1.,3.);
    h_mu_F_sf_up_   = fs->make<TH1D>("mu_F_sf_up",  ";mu_F_sf_up",   100,-1.,3.);
    h_mu_F_sf_down_ = fs->make<TH1D>("mu_F_sf_down",";mu_F_sf_down", 100,-1.,3.);
    h_scale_comb_sf_up_   = fs->make<TH1D>("scale_comb_sf_up",  ";scale_comb_sf_up",   100,-1.,3.);
    h_scale_comb_sf_down_ = fs->make<TH1D>("scale_comb_sf_down",";scale_comb_sf_down", 100,-1.,3.);
    h_pdf_alphas_sf_      = fs->make<TH1D>("pdf_alphas_sf",     ";pdf_alphas_sf",      100,-1.,3.);
    h_pdf_alphas_sf_up_   = fs->make<TH1D>("pdf_alphas_sf_up",  ";pdf_alphas_sf_up",   100,-1.,3.);
    h_pdf_alphas_sf_down_ = fs->make<TH1D>("pdf_alphas_sf_down",";pdf_alphas_sf_down", 100,-1.,3.);

    edm::EDGetTokenT<std::vector<float>>(mayConsume<std::vector<float>>(edm::InputTag("extraVar", "scaleWeights")));
    edm::EDGetTokenT<std::vector<float>>(mayConsume<std::vector<float>>(edm::InputTag("extraVar", "pdfWeights")));
    edm::EDGetTokenT<std::vector<float>>(mayConsume<std::vector<float>>(edm::InputTag("extraVar", "alphasWeights")));
  }

private:
  TH1D* h_pileup_;
  TH1D* h_weightsign_;
  TH1D* h_totweight_;
  TH2D* h_cstar_;
  TH2D* h_x_F_;
  TH2D* h_M_;
  TH2D* h_initial_parton_IDs_;
  TH2D* h_cstar_vs_beta_qqbar_;
  TH2D* h_cstar_vs_beta_gg_;
  TH1D* h_mu_R_sf_up_;
  TH1D* h_mu_R_sf_down_;
  TH1D* h_mu_F_sf_up_;
  TH1D* h_mu_F_sf_down_;
  TH1D* h_scale_comb_sf_up_;
  TH1D* h_scale_comb_sf_down_;
  TH1D* h_pdf_alphas_sf_;
  TH1D* h_pdf_alphas_sf_up_;
  TH1D* h_pdf_alphas_sf_down_;
  void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
    if ( iEvent.eventAuxiliary().isRealData() ) {
      h_weightsign_->Fill(0);
      h_totweight_->Fill(0);
    } else {

      //Get the generated weight to apply to filling the histograms
      edm::Handle<float> evt_Gen_Weight;
      iEvent.getByLabel(edm::InputTag("extraVar", "evtGenWeight"), evt_Gen_Weight);

      // Save Pileup distribution for MC
      edm::Handle<int>  NTrueInt;
      iEvent.getByLabel(edm::InputTag("eventUserData", "puNtrueInt"), NTrueInt);
      h_pileup_->Fill(*NTrueInt,*evt_Gen_Weight);
      
      // Save sign of weights and total weight
      h_weightsign_->Fill((*evt_Gen_Weight>=0)-1);
      h_totweight_->Fill(double(0), *evt_Gen_Weight);

      //Save the 2D histogram of initial parton IDs
      edm::Handle<float>  MCp1ID;
      edm::Handle<float>  MCp2ID;
      iEvent.getByLabel(edm::InputTag("extraVar", "MCpart1ID"), MCp1ID);
      iEvent.getByLabel(edm::InputTag("extraVar", "MCpart2ID"), MCp2ID);
      h_initial_parton_IDs_->Fill(*MCp1ID, *MCp2ID, *evt_Gen_Weight);

      //Save MC truth cstar, x_F, and M distributions
      edm::Handle<float>  MCtruthcstar;
      iEvent.getByLabel(edm::InputTag("extraVar", "MCcstar"), MCtruthcstar);
      edm::Handle<float>  MCtruthxF;
      iEvent.getByLabel(edm::InputTag("extraVar", "MCxF"), MCtruthxF);
      edm::Handle<float>  MCtruthM;
      iEvent.getByLabel(edm::InputTag("extraVar", "MCMtt"), MCtruthM);
      //And do the event type split
      if ((*MCp1ID)+(*MCp2ID)==0) { //qqbar
        h_cstar_->Fill(0.,*MCtruthcstar,*evt_Gen_Weight);
        h_x_F_->Fill(0.,*MCtruthxF,*evt_Gen_Weight);
        h_M_->Fill(0.,*MCtruthM,*evt_Gen_Weight);
      }
      else if ((*MCp1ID)==21 && (*MCp2ID)==21) { //gg
        h_cstar_->Fill(2.,*MCtruthcstar,*evt_Gen_Weight);
        h_x_F_->Fill(2.,*MCtruthxF,*evt_Gen_Weight);
        h_M_->Fill(2.,*MCtruthM,*evt_Gen_Weight);
      }
      else if ((*MCp1ID)==21 || (*MCp2ID)==21) { //qg
        h_cstar_->Fill(1.,*MCtruthcstar,*evt_Gen_Weight);
        h_x_F_->Fill(1.,*MCtruthxF,*evt_Gen_Weight);
        h_M_->Fill(1.,*MCtruthM,*evt_Gen_Weight);
      }

      //Save the MC truth cstar vs. beta distributions
      edm::Handle<float>  MCtpt;
      edm::Handle<float>  MCteta;
      edm::Handle<float>  MCtphi;
      edm::Handle<float>  MCtE;
      edm::Handle<float>  MCtbarpt;
      edm::Handle<float>  MCtbareta;
      edm::Handle<float>  MCtbarphi;
      edm::Handle<float>  MCtbarE;
      iEvent.getByLabel(edm::InputTag("extraVar", "MCtpt"), MCtpt);
      iEvent.getByLabel(edm::InputTag("extraVar", "MCteta"), MCteta);
      iEvent.getByLabel(edm::InputTag("extraVar", "MCtphi"), MCtphi);
      iEvent.getByLabel(edm::InputTag("extraVar", "MCtE"), MCtE);
      iEvent.getByLabel(edm::InputTag("extraVar", "MCtbarpt"), MCtbarpt);
      iEvent.getByLabel(edm::InputTag("extraVar", "MCtbareta"), MCtbareta);
      iEvent.getByLabel(edm::InputTag("extraVar", "MCtbarphi"), MCtbarphi);
      iEvent.getByLabel(edm::InputTag("extraVar", "MCtbarE"), MCtbarE);
      TLorentzVector MCt, MCtbar;
      MCt.SetPtEtaPhiE(*MCtpt,*MCteta,*MCtphi,*MCtE);
      MCtbar.SetPtEtaPhiE(*MCtbarpt,*MCtbareta,*MCtbarphi,*MCtbarE);
      double M2_1 = MCt.Mag2();
      double M2_2 = MCtbar.Mag2();
      double beta=0.0;
      //std::cout<<"MC truth M = "<<*MCtruthM<<"\n"; //DEBUG
      if (*MCtruthM!=-9999 && (1. - 2.*(M2_1+M2_2)/((*MCtruthM)*(*MCtruthM)) + (M2_1-M2_2)*(M2_1-M2_2)/((*MCtruthM)*(*MCtruthM)*(*MCtruthM)*(*MCtruthM)))>0.) 
        beta = sqrt(1. - 2.*(M2_1+M2_2)/((*MCtruthM)*(*MCtruthM)) + (M2_1-M2_2)*(M2_1-M2_2)/((*MCtruthM)*(*MCtruthM)*(*MCtruthM)*(*MCtruthM)));
      //depending on production mechanism
      if ((*MCp1ID)+(*MCp2ID)==0) { //qqbar
        h_cstar_vs_beta_qqbar_->Fill(beta,*MCtruthcstar,0.5*(*evt_Gen_Weight));
        h_cstar_vs_beta_qqbar_->Fill(beta,-1.*(*MCtruthcstar),0.5*(*evt_Gen_Weight));
      }
      else if ((*MCp1ID)==21 && (*MCp2ID)==21) { //gg
        h_cstar_vs_beta_gg_->Fill(beta,*MCtruthcstar,0.5*(*evt_Gen_Weight));
        h_cstar_vs_beta_gg_->Fill(beta,-1.*(*MCtruthcstar),0.5*(*evt_Gen_Weight));
      }
      else if ((*MCp1ID)==21 || (*MCp2ID)==21) { //qg
        h_cstar_vs_beta_gg_->Fill(beta,*MCtruthcstar,*evt_Gen_Weight);
      }

      //Save the scale/PDF/alpha_s weights
      edm::Handle<std::vector<float>>  scaleWeights_hand;
      edm::Handle<std::vector<float>>  pdfWeights_hand;
      edm::Handle<std::vector<float>>  alphasWeights_hand;
      iEvent.getByLabel(edm::InputTag("extraVar", "scaleWeights"), scaleWeights_hand);
      iEvent.getByLabel(edm::InputTag("extraVar", "pdfWeights"), pdfWeights_hand);
      iEvent.getByLabel(edm::InputTag("extraVar", "alphasWeights"), alphasWeights_hand);
      //scale weights are pretty simple because they're right in the structure
      const std::vector<float> scaleWeights = * (scaleWeights_hand.product());
      if (scaleWeights.size()>=6) {
        h_mu_R_sf_up_->Fill(scaleWeights[2],*evt_Gen_Weight);
        h_mu_R_sf_down_->Fill(scaleWeights[4],*evt_Gen_Weight);
        h_mu_F_sf_up_->Fill(scaleWeights[0],*evt_Gen_Weight);
        h_mu_F_sf_down_->Fill(scaleWeights[1],*evt_Gen_Weight);
        h_scale_comb_sf_up_->Fill(scaleWeights[3],*evt_Gen_Weight);
        h_scale_comb_sf_down_->Fill(scaleWeights[5],*evt_Gen_Weight);
      }
      //pdf/alpha_s weights are more complicated
      const std::vector<float> pdfWeights = * (pdfWeights_hand.product());
      const std::vector<float> alphasWeights = * (alphasWeights_hand.product());
      size_t pdflen = pdfWeights.size();
      float sumweights  = 0.;
      float sumweights2 = 0.;
      for (size_t i=0; i<pdflen; ++i) {
        float val = pdfWeights[i];
        sumweights+=val;
        sumweights2+=val*val;
      }
      if (pdflen>0) {
        float pdfmean = sumweights/pdflen;
        float pdfunc = sqrt(abs(sumweights2)/pdflen-pdfmean*pdfmean);
        h_pdf_alphas_sf_->Fill(pdfmean,*evt_Gen_Weight);
        if (alphasWeights.size()>1) {
          float alphas_unc_up   = abs(alphasWeights[1]*0.75-1.);
          float alphas_unc_down = abs(alphasWeights[0]*0.75-1.);
          h_pdf_alphas_sf_up_->Fill(pdfmean+sqrt(pdfunc*pdfunc+alphas_unc_up*alphas_unc_up),*evt_Gen_Weight);
          h_pdf_alphas_sf_down_->Fill(pdfmean-sqrt(pdfunc*pdfunc+alphas_unc_down*alphas_unc_down),*evt_Gen_Weight);
        }
        else {
          h_pdf_alphas_sf_up_->Fill(pdfmean+pdfunc,*evt_Gen_Weight);
          h_pdf_alphas_sf_down_->Fill(pdfmean-pdfunc,*evt_Gen_Weight);
        }
      }
      else {
        h_pdf_alphas_sf_->Fill(1,*evt_Gen_Weight);
        h_pdf_alphas_sf_up_->Fill(1,*evt_Gen_Weight);
        h_pdf_alphas_sf_down_->Fill(1,*evt_Gen_Weight);
      }
    }
  }
};

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(EventCounter);
