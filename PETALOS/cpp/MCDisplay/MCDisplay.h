#ifndef _MCDisplay__
#define _MCDisplay__

#include <GATE/Centella.h>

class MCDisplay : public gate::IAlgo {

 public:
  
  //! default contructor
  MCDisplay(gate::VLEVEL=gate::NORMAL,
	       std::string label="MCDisplayInstance");
  
  //! constructor with store with input parameters 
  MCDisplay(const gate::ParamStore& gs,
	       gate::VLEVEL=gate::NORMAL,
	       std::string label="MCDisplayInstance");
  
  //! destructor
  virtual ~MCDisplay(){};
  
  //! initialize algorithm
  bool initialize();        
  
  //! execute algorithm: process current event
  bool execute(gate::Event& evt);  
  
  //! finalize algorithm
  bool finalize();     

  //! display event information
  void EventInfo(const gate::Event& event, ostream& s);
  void SensorHitInfo(const gate::Hit& hit, ostream& s);
  void TrueHitInfo(const gate::MCHit& hit, ostream& s);
  void ParticleInfo(const gate::MCParticle& hit, ostream& s);

  // temporary 
  std::string GetParticleName(int pdg);
  double GetParticleMass(int pdg);
  double GetParticleCharge(int pdg);
  
 private:
  
  ClassDef(MCDisplay,0)
    
};

//ostream& operator << (ostream& s, const gate::Event& ev);

#endif
