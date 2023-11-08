#ifndef Event_H
#define Event_H
#include <vector>
#include "Particle.h"

class Event {
 public:
  Event(Int_t nTrack, Float_t b, Int_t nPartP, Int_t nPartT, Int_t nElP, Int_t nInelP, Int_t nElT, Int_t nInelT, vector<Particle> particles)
    : _NumOfTracks(nTrack),
      _ImpactParameter(b),
      _nPartProjectile(nPartP),
      _nPartTarget(nPartT),
      _nPartProjectileElastic(nElP),
      _nPartProjectileInelastic(nInelP),
      _nPartTargetElastic(nElT),
      _nPartTargetInelastic(nInelT),
      _Particles(particles){};
    
  Event() = default;
  Particle getParticle(Int_t ith) { return _Particles.at(ith); };
  Int_t getNtrack() { return _NumOfTracks; };
  Float_t getb() { return _ImpactParameter; };
  Int_t getnPartP() { return _nPartProjectile; };
  Int_t getnPartT() { return _nPartTarget; };
  Int_t getnPart() { return _nPartProjectile+_nPartTarget; };
  Int_t getnElP() { return _nPartProjectileElastic; };
  Int_t getnInelP() { return _nPartProjectileInelastic; };
  Int_t getnElT() { return _nPartTargetElastic; };
  Int_t getnInelT() { return _nPartTargetInelastic; };

 private:
    std::vector<Particle> _Particles;
    Int_t _NumOfTracks;
    Float_t _ImpactParameter;
    Int_t _nPartProjectile;
    Int_t _nPartTarget;
    Int_t _nPartProjectileElastic;
    Int_t _nPartProjectileInelastic;
    Int_t _nPartTargetElastic;
    Int_t _nPartTargetInelastic;
};
#endif
