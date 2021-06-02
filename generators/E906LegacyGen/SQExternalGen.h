#ifndef __SQExternalGen_H__
#define __SQExternalGen_H__

#include <vector>
#include <TString.h>
#include <g4main/PHG4ParticleGeneratorBase.h>

class PHCompositeNode;
class PHG4ParticleGeneratorBase;
class PHG4InEvent;
class PHG4Particle;

class TFile;
class TTree;
class TClonesArray;

class SQMCEvent;
class SQPrimaryVertexGen;

class ExtParticle
{
public:
  ExtParticle(int evtID, int pdg, const TVector3& pos, const TVector3& mom);

public:
  int _evtID;
  int _pdg;
  TVector3 _pos;
  TVector3 _mom;
};

class SQExternalGen: public PHG4ParticleGeneratorBase
{
public:
  SQExternalGen(const std::string& name = "ExtGen");
  virtual ~SQExternalGen();

  int Init(PHCompositeNode* topNode);
  int InitRun(PHCompositeNode* topNode);
  int process_event(PHCompositeNode* topNode);

  void setExtInputFile(const TString& name) { _extFileName = name; }
  void setBucketSize(double n) { _bucketSize = n; }

private:
  //! read n events and fill the particle container, return false if there is no enough events left
  bool readExtTree(int nEvents);

  //! pointer to the instance of vertex generator
  SQPrimaryVertexGen* _vertexGen;

  //! MC truth container
  PHG4InEvent* _ineve;
  SQMCEvent*   _mcevt;

  //! internal container for quick particle cache
  std::vector<ExtParticle> _extParticles;

  //! Number of protons in each RF bucket - TODO: may need to incorporate a profile instead of fixed number
  int _bucketSize;

  //! External input file and parameters
  TString _extFileName;
  TFile* _extFile;
  TTree* _extTree;
  int _extEvtID;
  int _readIdx;
  int _nExtPar;
  int _extPDG[10000];
  TClonesArray* _extPos;
  TClonesArray* _extMom;
};

#endif