#ifndef Particle_H
#define Particle_H
#define _USE_MATH_DEFINES
#include <cmath>

class Particle {
public:
    Particle(Int_t pid, Int_t specFlag, Float_t px, Float_t py, Float_t pz, Float_t mass, Float_t charge, Float_t x, Float_t y, Float_t z, Float_t t)
    : _ParticleIdentification(pid),
    _SpectatorFlag(specFlag),
    _Px(px),
    _Py(py),
    _Pz(pz),
    _ParticleMass(mass),
    _ParticleCharge(charge),
    _X(x),
    _Y(y),
    _Z(z),
    _time(t){};
    
    Particle() = default;
    
    Int_t getPid() { return _ParticleIdentification; }
    Int_t getSpecflag() { return _SpectatorFlag; }
    Float_t getPx() { return _Px; }
    Float_t getPy() { return _Py; }
    Float_t getPz() { return _Pz; }
    Float_t getMass() { return _ParticleMass; }
    Float_t getCharge() { return _ParticleCharge; }
    Float_t getX() { return _X; }
    Float_t getY() { return _Y; }
    Float_t getZ() { return _X; }
    Float_t getT() { return _time; }
    Float_t getPt() { return sqrt(_Px*_Px+_Py*_Py); }
    Float_t getEnergy() { return sqrt(_Px*_Px + _Py*_Py + _Pz*_Pz + _ParticleMass); }
    Float_t getRapidity() { return -1./2.*log((sqrt(_Px*_Px + _Py*_Py + _Pz*_Pz)+_Pz)/(sqrt(_Px*_Px + _Py*_Py + _Pz*_Pz)-_Pz)); }
    Float_t getPhi() {return atan2(_Py,_Px)+M_PI; }
    
private:
    Int_t _ParticleIdentification{0};
    Int_t _SpectatorFlag{-999};
    Float_t _Px{-999.};
    Float_t _Py{-999.};
    Float_t _Pz{-999.};
    Float_t _ParticleMass{-999.};
    Float_t _ParticleCharge{-999.};
    Float_t _X{-999.};
    Float_t _Y{-999.};
    Float_t _Z{-999.};
    Float_t _time{-999.};
};
#endif
