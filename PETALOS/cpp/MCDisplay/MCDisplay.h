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

  // load geometry
  void LoadGeometry(const gate::ParamStore& gs,
         gate::VLEVEL=gate::NORMAL,
         std::string label="MCDisplayInstance");

  //! display event information
  void EventInfo(const gate::Event& event);
  std::string SensorHitInfo(const gate::Hit& hit);
  std::string TrueHitInfo(const gate::MCHit& hit);
  std::string ParticleInfo(const gate::MCParticle& hit);

  // temporary 
  std::string GetParticleName(int pdg);
  double GetParticleMass(int pdg);
  double GetParticleCharge(int pdg);
  void Wait();
  
 private:
  std::vector<std::vector<double>> fCoordBox1; // coordinates of Box1
  std::vector<std::vector<double>> fCoordBox2; // coordinates of Box2
  
  ClassDef(MCDisplay,0)
    
};

//ostream& operator << (ostream& s, const gate::Event& ev);

#endif
