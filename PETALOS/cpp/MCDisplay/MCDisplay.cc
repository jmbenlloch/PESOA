
#include<MCDisplay.h>

#include <iomanip>

ClassImp(MCDisplay)

//==========================================================================
MCDisplay::MCDisplay(gate::VLEVEL vl, std::string label) : 
IAlgo(vl,"MCDisplay",0,label)
//==========================================================================
{

}

//==========================================================================
MCDisplay::MCDisplay(const gate::ParamStore& gs, 
			   gate::VLEVEL vl, std::string label) :
  IAlgo(gs,vl,"MCDisplay",0,label)
//==========================================================================
{
    LoadGeometry(gs,vl,label);
}

//==========================================================================
bool MCDisplay::initialize()
//==========================================================================
{
  _m.message("Intializing algorithm",this->getAlgoLabel(),gate::NORMAL);

  gate::Centella::instance()
    ->hman()->h1(this->alabel("EvtID"),"EvtID",10,0,100);

  gate::Centella::instance()
    ->hman()->h1(this->alabel("Edep"),"Total deposited energy of the event (keV)",100,0,600);

  gate::Centella::instance()
    ->hman()->h1(this->alabel("TrueHits"),"True hits",100,0,100);

  gate::Centella::instance()
    ->hman()->h1(this->alabel("SensorHits"),"Sensor hits",100,0,100);

  gate::Centella::instance()
    ->hman()->h1(this->alabel("MCParticles"),"Monte Carlo Particles",25,0,25);

  _m.message("Defining TrueHitsXY, xmin, xmax,ymin,ymax",this->fetch_dstore("xmin"),
    this->fetch_dstore("xmax"),this->fetch_dstore("ymin"), this->fetch_dstore("ymax"),
    gate::NORMAL);

  gate::Centella::instance()
    ->hman()->h2(this->alabel("TrueHitsXY"),"x-y coordinates, true hits (mm)",
      25, this->fetch_dstore("xmin"), this->fetch_dstore("xmax"),
      25,this->fetch_dstore("ymin"), this->fetch_dstore("ymax"));

  gate::Centella::instance()
    ->hman()->h1(this->alabel("TrueHitsZ"),"z coordinates true hits",25,this->fetch_dstore("zmin1"),
      this->fetch_dstore("zmax2"));

  // gate::Centella::instance()
  //   ->hman()->h1(this->alabel("TrueHitsZBox1"),"TrueHitsZBox1",25,this->fetch_dstore("zmin1"),
  //     this->fetch_dstore("zmax1"));

  // gate::Centella::instance()
  //   ->hman()->h1(this->alabel("TrueHitsZBox2"),"TrueHitsZBox2",25,this->fetch_dstore("zmin2"),
  //     this->fetch_dstore("zmax2"));

  gate::Centella::instance()
    ->hman()->h1(this->alabel("TrueHitsTime"),"TrueHitstime",25,0,
      1e+6);

  return true;

}

//==========================================================================
bool MCDisplay::execute(gate::Event& evt)
//==========================================================================
{
  _m.message("Executing algorithm",this->getAlgoLabel(),gate::VERBOSE);
  
  
  EventInfo(evt);
  
  return true;

}

//==========================================================================
bool MCDisplay::finalize()
//==========================================================================
{ 
  _m.message("Finalizing algorithm",this->getAlgoLabel(),gate::NORMAL);
  return true;
}

//==========================================================================
  void MCDisplay::LoadGeometry(const gate::ParamStore& gs,
         gate::VLEVEL,std::string label)
//==========================================================================
/* geometry of boxes

  z1 = -130 ---- z2 = -100  ---- z3 = 100 ---- z4 = 130 (all in mm)
  (x1,y1) = (-12.8, -12.8), (x2,y2) = (12.8,12.8)
 
  -12.8 -12.8 -130 (x1,y1,z1)
  -12.8 12.8 -130  (x1,y2,z1)
  12.8 -12.8 -130  (x2,y1,z1)
  12.8 12.8 -130   (x2,y1,z1)
  -12.8 -12.8 -100 (x1,y1,z2)
  -12.8 12.8 -100  (x1,y2,z2)
  12.8 -12.8 -100  (x2,y1,z2)
  12.8 12.8 -100   (x2,y2,z2)

  -12.8 -12.8 100 (x1,y1,z3)
  -12.8 12.8 100  (x1,y2,z3)
  12.8 -12.8 100  (x2,y1,z3)
  12.8 12.8 100   (x2,y1,z3)
  -12.8 -12.8 130 (x1,y1,z4)
  -12.8 12.8 130  (x1,y2,z4)
  12.8 -12.8 130  (x2,y1,z4)
  12.8 12.8 130   (x2,y2,z4)

  box labels: Box1V1 to Box1V8 and Box2V1 to Box2V8
   
  */
{
  std::vector<double> coord;

  {
    
    for(int i=1; i<9; i++)
    {
      std::string blb="Box1V";
      gate::to_string(blb, i);
      _m.message("loading now parameter",  blb, gate::VERBOSE);
      coord  = gs.fetch_dvstore(blb);
      fCoordBox1.push_back(coord);
    }
  }

  {
    
    for(int i=1; i<9; i++)
    {
      std::string blb="Box2V";
      gate::to_string(blb, i);
      _m.message("loading now parameter",  blb, gate::VERBOSE);
      coord  = gs.fetch_dvstore(blb);
      fCoordBox2.push_back(coord);
    }
  }
  
  for(int i=0; i<fCoordBox1.size(); i++)
    _m.message("Coordinates of Box 1",  gate::print_vector(fCoordBox1.at(i)), gate::NORMAL);

  for(int i=0; i<fCoordBox2.size(); i++)
    _m .message("Coordinates of Box 2", gate::print_vector(fCoordBox2.at(i)), gate::NORMAL);

       
  std::vector<double> b1v1 = fCoordBox1.at(0); 
  _m.message("Coordinates of box1 vertex 1 ",gate::print_vector(b1v1),gate::VERBOSE);
  double xmin = b1v1.at(0);
  double ymin = b1v1.at(1);
  double zmin1 = b1v1.at(2);

  std::vector<double> b1v8 = fCoordBox1.at(7); 
  _m.message("Coordinates of box 1, vertex 8",gate::print_vector(b1v8),gate::VERBOSE);
  double xmax = b1v8.at(0);
  double ymax = b1v8.at(1);
  double zmax1 = b1v8.at(2);

  std::vector<double> b2v1 = fCoordBox2.at(0); 
  _m.message("Coordinates of box2 vertex 1 ",gate::print_vector(b2v1),gate::VERBOSE);
  double zmin2 = b2v1.at(2);

  std::vector<double> b2v8 = fCoordBox2.at(7); 
  _m.message("Coordinates of box 2, vertex 8",gate::print_vector(b2v8),gate::VERBOSE);
  double zmax2 = b2v8.at(2);


  _m.message(" XY Coordinates of box1-box2: xmin, xmax, ymin, ymax ",xmin,xmax,ymin,ymax,
    gate::VERBOSE);

   _m.message(" z Coordinates of box1-box2: zmin1, zmax1, zmin2, zmax2 ",zmin1,zmax1,zmin2,zmax2,
    gate::VERBOSE);
   this->store("xmin",xmin);
   this->store("xmax",xmax);
   this->store("ymin",ymin);
   this->store("ymax",ymax);
   this->store("zmin1",zmin1);
   this->store("zmax1",zmax1);
   this->store("zmin2",zmin2);
   this->store("zmax2",zmax2);
}
//==========================================================================
void MCDisplay::EventInfo(const gate::Event& event)
//==========================================================================
{

  _m.message("event number = ",event.GetID(),gate::NORMAL);

  gate::Centella::instance()
    ->hman()->fill(this->alabel("EvtID"),event.GetEventID());


  _m.message("Total deposited energy of the event (keV) = ",
    event.GetMCEnergy()/gate::keV,gate::VERBOSE);

  gate::Centella::instance()
    ->hman()->fill(this->alabel("Edep"),event.GetMCEnergy()/gate::keV);

  _m.message("Number of true hits = ",
    event.GetMCHits().size(),gate::VERBOSE);

  gate::Centella::instance()
    ->hman()->fill(this->alabel("TrueHits"),event.GetMCHits().size());


  _m.message("Number of sensor hits = ",
    event.GetMCSensHits().size(),gate::VERBOSE);

  gate::Centella::instance()
    ->hman()->fill(this->alabel("SensorHits"),event.GetMCSensHits().size());


  _m.message("Number of MC particles = ",
    event.GetMCParticles().size(),gate::VERBOSE);
  gate::Centella::instance()
    ->hman()->fill(this->alabel("MCParticles"),event.GetMCParticles().size());

  _m.message("------------------------------------",gate::VERBOSE);

  Wait();
  _m.message("List of true hits in the event",gate::VVERBOSE);
  _m.message("------------------------------------",gate::VVERBOSE);
  
  const std::vector<gate::MCHit*> mchits =  event.GetMCHits();

  for (unsigned int ihit=0; ihit<mchits.size(); ++ihit) 
  {
    _m.message("Hit number ",ihit, gate::VVERBOSE);
    const gate::MCHit* myhit = mchits.at(ihit);
    std::string thits = TrueHitInfo(*myhit);
    _m.message("true hits-->",thits, gate::VVERBOSE);
  }

  Wait();

  _m.message("List of sensor hits in the event",gate::VVERBOSE);
  _m.message("------------------------------------",gate::VVERBOSE);
    
  const std::vector<gate::Hit*> sensorhits =  event.GetMCSensHits();
  for (unsigned int ihit=0; ihit<sensorhits.size(); ++ihit) 
  {
    _m.message("Hit number ",ihit, gate::VVERBOSE);
    const gate::Hit* myhit = sensorhits.at(ihit);
    std::string thits = SensorHitInfo(*myhit);
    _m.message("sensor hits-->",thits, gate::VVERBOSE);   
  }

  Wait();
  _m.message("List of particles in the event",gate::VERBOSE);
  _m.message("------------------------------------",gate::VERBOSE);

  const std::vector<gate::MCParticle*> particles =  event.GetMCParticles();
  for (unsigned int ipart=0; ipart<particles.size(); ++ipart) 
  {
    const gate::MCParticle* mypart=particles.at(ipart);
    std::string part = ParticleInfo(*mypart);
    _m.message("particles--->",part, gate::VERBOSE); 
  }
  Wait();
}
//==========================================================================
std::string MCDisplay::TrueHitInfo(const gate::MCHit& hit)
//==========================================================================
{
  
  _m.message("Entering TrueHitInfo",gate::VVERBOSE);
  _m.message("x,y,z (mm) = ",hit.GetPosition().x(), hit.GetPosition().y(), hit.GetPosition().z(),
    gate::VERBOSE);

  _m.message(" energy (MeV), time (?) = ",hit.GetAmplitude(),hit.GetTime(),
    gate::VERBOSE);

  _m.message("Filling histogram Z",gate::VVERBOSE);
    gate::Centella::instance()
    ->hman()->fill(this->alabel("TrueHitsZ"),hit.GetPosition().z()/gate::mm);

  _m.message("Filling histogram Time",gate::VVERBOSE);
  gate::Centella::instance()
    ->hman()->fill(this->alabel("TrueHitsTime"),hit.GetTime());

  _m.message("Filling histogram XY",gate::VVERBOSE);

  gate::Centella::instance()
    ->hman()->fill2d(this->alabel("TrueHitsXY"),hit.GetPosition().x()/gate::mm,
      hit.GetPosition().y()/gate::mm);


  _m.message("serialize info",gate::VVERBOSE);

  std::ostringstream s;
  s << std::endl;

  s  << " x (mm)    y (mm)    z (mm)    " << std::endl;
  s << std::setw(4) << hit.GetPosition().x() <<"     " 
    << std::setw(4) << hit.GetPosition().y() <<"     " 
    << std::setw(4) << hit.GetPosition().z() << std::endl;
  s << "deposited energy = "<< hit.GetAmplitude()  << " MeV"<< std::endl;
  s << "time = "<< hit.GetTime() << std::endl;
  s << "--------------------------------------" << std::endl;


  return s.str();
}
//==========================================================================
std::string MCDisplay::SensorHitInfo(const gate::Hit& hit)
//==========================================================================
{
  std::ostringstream s;
  s << std::endl;

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
  return s.str();

}

//==========================================================================
std::string MCDisplay::ParticleInfo(const gate::MCParticle& particle)
//==========================================================================
{
  std::ostringstream s;
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

  if (particle.IsPrimary()) 
  {
    s << "particle is primary " << std::endl;
  } 
  else 
  {
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
 
  for (unsigned int i=0; i<particle.GetDaughters().size(); ++i)
  {
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
	     <<  std::sqrt(std::pow(p->GetInitialMom().x(),2) + 
        std::pow(p->GetInitialMom().y(),2) + std::pow(p->GetInitialMom().z(),2)) 
       << std::endl;      
    s << "particle energy (MeV) = " << p->GetInitialMom().GetE()  << std::endl;      
  }

  s << std::endl;
  s << "List of true hits of the particle"
    << " ----------------------" << std::endl;

  const std::vector<const gate::MCTrack*> 
  tracks = (const std::vector<const gate::MCTrack*>&)particle.GetTracks();
    
  for (unsigned int i=0; i<tracks.size(); ++i) 
  {
    const gate::MCTrack* t = tracks.at(i);
    s << std::endl;
    s << "Detector " << t->GetLabel() << std::endl;
    s << std::endl;
    std::vector<gate::MCHit*> hits = (std::vector<gate::MCHit*> &)t->GetHits() ;
    for (unsigned int ihit=0; ihit<hits.size(); ++ihit) 
    {
	   const gate::MCHit* myhit = hits.at(ihit);
	   TrueHitInfo(*myhit);
    }
    s << "Total energy of the track associated with this particle (MeV) = " << t->GetEnergy() << std::endl;
  }
  
  s << "Length of the total associated track (mm) = " << particle.GetPathLength() << std::endl;
    return s.str();
}
//==========================================================================
std::string MCDisplay::GetParticleName(int pdg)
//==========================================================================
{
  std::string name;
  if (pdg == 11) name = "e-";
  else if (pdg == -11) name = "e+";
  else if (pdg == 12) name = "nu_e";
  else if (pdg == -12) name = "anti_nu_e";
  else if (pdg == 22) name = "gamma";

  return name;
}
//==========================================================================
double MCDisplay::GetParticleMass(int pdg)
//==========================================================================
{
  double mass;
  if (pdg == 11 or pdg == -11) 
    mass = 0.510998902;
  else  
    mass = 0;
 
  return mass;
}
//==========================================================================
double MCDisplay::GetParticleCharge(int pdg)
//==========================================================================
{
  double charge;
  if (pdg == 11 or pdg == -11) 
    charge = -1.;
  else 
    charge = 0;

  return charge;
}

//==========================================================================
void MCDisplay::Wait()
//==========================================================================
{
  std::cout << "Hit enter to continue\n";
  std::cin.ignore();
}
