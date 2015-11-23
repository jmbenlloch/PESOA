
#include<MCDisplay.h>

#include <iomanip>

ClassImp(MCDisplay)

//==========================================================================
MCDisplay::MCDisplay(gate::VLEVEL vl, std::string label) : 
IAlgo(vl,"MCDisplay",0,label){
//==========================================================================


}

//==========================================================================
MCDisplay::MCDisplay(const gate::ParamStore& gs, 
			   gate::VLEVEL vl, std::string label) :
  IAlgo(gs,vl,"MCDisplay",0,label){
//==========================================================================


}

//==========================================================================
bool MCDisplay::initialize(){
//==========================================================================

  _m.message("Intializing algorithm",this->getAlgoLabel(),gate::NORMAL);
  
  gate::Centella::instance()
    ->hman()->h1(this->alabel("EvtID"),"EvtID",10,0,100);

  return true;

}

//==========================================================================
bool MCDisplay::execute(gate::Event& evt){
//==========================================================================

  _m.message("Executing algorithm",this->getAlgoLabel(),gate::VERBOSE);
  
  _m.message("Event number:",evt.GetEventID(),gate::VERBOSE);
  
  gate::Centella::instance()
    ->hman()->fill(this->alabel("EvtID"),evt.GetEventID());

  EventInfo(evt, std::cout);
  
  return true;

}

//==========================================================================
bool MCDisplay::finalize(){
//==========================================================================

  _m.message("Finalising algorithm",this->getAlgoLabel(),gate::NORMAL);
  
  return true;

}

void MCDisplay::EventInfo(const gate::Event& event, ostream& s)
{
  s << std::endl; 
  s << "event number = " << event.GetID() << std::endl;

  // std::vector <gate::Track*> tracks = event.GetTracks();
  // int tothits=0;
  // for (unsigned int itrack=0; itrack<tracks.size(); ++itrack) {
  //   gate::Track* mytrack = tracks.at(itrack);     
  //   tothits = tothits +  mytrack->GetMCHits().size();
  // }
  s << "Total deposited energy of the event = " << event.GetMCEnergy() << " MeV "<< std::endl;
  s << "event has " << event.GetMCHits().size() << " true hits" << std::endl;
  s << "event has " <<  event.GetMCSensHits().size() << " sensor hits" << std::endl;
  s << "event has " <<  event.GetMCParticles().size() << " particles" << std::endl;
  s << std::endl;  

  s << " List of sensor hits in the event"
    << "------------------------------------" << std::endl;
  s << std::endl;  
  const std::vector<gate::Hit*> sensorhits =  event.GetMCSensHits();
  for (unsigned int ihit=0; ihit<sensorhits.size(); ++ihit) {
    const gate::Hit* myhit = sensorhits.at(ihit);
    SensorHitInfo(*myhit, s);
    s << std::endl;
  }

    s << " List of true hits in the event"
      << "------------------------------------" << std::endl;
    s << std::endl;   
    const std::vector<gate::MCHit*> mchits =  event.GetMCHits();
    for (unsigned int ihit=0; ihit<mchits.size(); ++ihit) {
      const gate::MCHit* myhit = mchits.at(ihit);
      TrueHitInfo(*myhit, s);
      s << std::endl;
    }

    s << " List of particles in the event"
      << "------------------------------------" << std::endl;
    s << std::endl;
    const std::vector<gate::MCParticle*> particles =  event.GetMCParticles();
    for (unsigned int ipart=0; ipart<particles.size(); ++ipart) {
      const gate::MCParticle* mypart=particles.at(ipart);
      ParticleInfo(*mypart, s);
    }
  
    s << std::endl;
}

void MCDisplay::SensorHitInfo(const gate::Hit& hit, ostream& s)
{
  s << hit.GetLabel() << " hit, ID = " << hit.GetSensorID()  
    << std::endl;
  s  << " x (mm)    y (mm)    z (mm)    " << std::endl;
  s << std::setw(5) << hit.GetPosition().x() <<"     " 
    << std::setw(5) << hit.GetPosition().y() <<"     " 
    << std::setw(5) << hit.GetPosition().z() << std::endl;
  s << "total charge = "<< hit.GetAmplitude() << " pes"<< std::endl;
  s << "Waveform" << std::endl;
  std::vector< std::pair<unsigned short, unsigned short> > wvf = hit.GetWaveform().GetData();
  for (unsigned int smp=0; smp<wvf.size(); ++smp) {   
    s <<  wvf[smp].first << ", "  << wvf[smp].second << std::endl;
  }

}

void MCDisplay::TrueHitInfo(const gate::MCHit& hit, ostream& s)
{
  s  << " x (mm)    y (mm)    z (mm)    " << std::endl;
  s << std::setw(4) << hit.GetPosition().x() <<"     " 
    << std::setw(4) << hit.GetPosition().y() <<"     " 
    << std::setw(4) << hit.GetPosition().z() << std::endl;
  s << "deposited energy = "<< hit.GetAmplitude()  << " MeV"<< std::endl;
  s << "time = "<< hit.GetTime() << std::endl;
  s << "--------------------------------------" << std::endl;
}

void MCDisplay::ParticleInfo(const gate::MCParticle& particle, ostream& s)
{
  s << std::endl;

  s << "Particle name = " << GetParticleName(particle.GetPDG()) 
    << ", PDG code = " << particle.GetPDG() 
    <<  ", mass (MeV) = " << GetParticleMass(particle.GetPDG()) 
    << ", charge = "<< GetParticleCharge(particle.GetPDG()) << std::endl;
    s << "++++at production vertex ++++" << std::endl;
    s << "particle 3 momentum (MeV) =" << std::endl;
    s << "(" << particle.GetInitialMom().x() << "," << particle.GetInitialMom().y() 
      << "," << particle.GetInitialMom().z() << ")" << std::endl;
    s << " momentum (MeV) = " 
      << std::sqrt(std::pow( particle.GetInitialMom().x(),2) + std::pow( particle.GetInitialMom().y(),2) + std::pow( particle.GetInitialMom().z(),2)) << std::endl;
    s << " energy (MeV)= " << particle.GetInitialMom().GetE() << std::endl;
    s << " vertex (mm)= " <<  std::endl;
    s << "(" << particle.GetInitialVtx().x() << "," << particle.GetInitialVtx().y() <<
      "," << particle.GetInitialVtx().z() << ")" << std::endl;
    s << "++++at decay vertex ++++" << std::endl;
     s << "particle 3 momentum (MeV) =" << std::endl;
    s << "(" << particle.GetFinalMom().x() << "," << particle.GetFinalMom().y() 
      << "," << particle.GetFinalMom().z() << ")" << std::endl;
    s << " momentum (MeV) = " 
      << std::sqrt(std::pow( particle.GetFinalMom().x(),2) + std::pow( particle.GetFinalMom().y(),2) + std::pow( particle.GetFinalMom().z(),2)) << std::endl;
    s << " energy (MeV)= " << particle.GetFinalMom().GetE() << std::endl;
    s << " vertex (mm)= " <<  std::endl;
    s << "(" << particle.GetFinalVtx().x() << "," << particle.GetFinalVtx().y() <<
      "," << particle.GetFinalVtx().z() << ")" << std::endl;
    s << "creator process = " << particle.GetCreatorProc() << std::endl;
    s << "origin volume = " << particle.GetInitialVol() << std::endl;
    s << "decay volume = " << particle.GetFinalVol() << std::endl;
    s << "Particle ID = " << particle.GetID() << std::endl;

    if (particle.IsPrimary()) {
      s << "particle is primary " << std::endl;
    } else {
      s << "particle is secondary" << std::endl;
     
      s << "mother of particle is " << GetParticleName(particle.GetMother().GetPDG() )<< std::endl;
      s << "with 3 momentum (MeV) ="  << std::endl;
      s << "(" << particle.GetMother().GetInitialMom().x() << "," << 
        particle.GetMother().GetInitialMom().y() << "," << 
        particle.GetMother().GetInitialMom().z() << ")" << std::endl;
      s << "and energy (MeV) = " << particle.GetMother().GetInitialMom().GetE() << std::endl;
      
    }

    s << " List of secondary particles "
      << "-----------------------------" << std::endl;
 
    for (unsigned int i=0; i<particle.GetDaughters().size(); ++i){
      const gate::MCParticle* p = particle.GetDaughters().at(i);
      s << "daughter name = " << GetParticleName(p->GetPDG())
        << ", daughter mass (MeV) = " << GetParticleMass(p->GetPDG())
        << ", daughter charge = " << GetParticleCharge(p->GetPDG())
        << std::endl;     
     
      s << "particle 3 momentum (MeV) = " << std::endl;
      s << "(" << p->GetInitialMom().x() << "," << 
        p->GetInitialMom().y() << "," << 
        p->GetInitialMom().z() << ")" << std::endl;
      s << "particle momentum (MeV) =" 
	<<  std::sqrt(std::pow(p->GetInitialMom().x(),2) + std::pow(p->GetInitialMom().y(),2) + std::pow(p->GetInitialMom().z(),2)) << std::endl;      
      s << "particle energy (MeV) = " << p->GetInitialMom().GetE()  << std::endl;      
    }

    s << std::endl;
    s << "List of true hits of the particle"
      << " ----------------------" << std::endl;

    const std::vector<const gate::MCTrack*> tracks 
      = (const std::vector<const gate::MCTrack*>&)particle.GetTracks();
    
    for (unsigned int i=0; i<tracks.size(); ++i) {
      const gate::MCTrack* t = tracks.at(i);
      s << std::endl;
      s << "Detector " << t->GetLabel() << std::endl;
      s << std::endl;
      std::vector<gate::MCHit*> hits = (std::vector<gate::MCHit*> &)t->GetHits() ;
      for (unsigned int ihit=0; ihit<hits.size(); ++ihit) {
	const gate::MCHit* myhit = hits.at(ihit);
	TrueHitInfo(*myhit, s);
      }
      s << "Total energy of the track associated with this particle (MeV) = " << t->GetEnergy() << std::endl;
    }
    s << "Length of the total associated track (mm) = " << particle.GetPathLength() << std::endl;
}

std::string MCDisplay::GetParticleName(int pdg)
{
  std::string name;
  if (pdg == 11) name = "e-";
  else if (pdg == -11) name = "e+";
  else if (pdg == 12) name = "nu_e";
  else if (pdg == -12) name = "anti_nu_e";
  else if (pdg == 22) name = "gamma";

  return name;
}

double MCDisplay::GetParticleMass(int pdg)
{
  double mass;
  if (pdg == 11) mass = 0.510998902;
  else if (pdg == -11) mass = 0.510998902;
  else if (pdg == 12) mass = 0;
  else if (pdg == -12) mass = 0.;
  else if (pdg == 22) mass = 0.;

  return mass;
}

double MCDisplay::GetParticleCharge(int pdg)
{
  double charge;
  if (pdg == 11) charge = -1.;
  else if (pdg == -11) charge = +1.;
  else if (pdg == 12) charge = 0;
  else if (pdg == -12) charge = 0.;
  else if (pdg == 22) charge = 0.;

  return charge;
}


