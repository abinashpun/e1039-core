#ifndef _SQGFRaveVertexing_H
#define _SQGFRaveVertexing_H

#include <fun4all/SubsysReco.h>
#include <TRandom1.h>
#include <map>                   
#include <string>
#include <vector>
#include "GFFitter.h"

namespace genfit
{
class GFRaveVertex;
class GFRaveVertexFactory;
class GFTrack;
} /* namespace genfit */


namespace SQGenFit
{
class GFFitter;
}


class TVector3;
class SQTrack;
class SQTrackVector;
class SQDimuonVector;
class SRecEvent;
class SRecTrack;
class SRecDimuon;
class PHCompositeNode;
class PHG4TruthInfoContainer;

class PHField;
class TFile;
class TTree;
class TGeoManager;
class TClonesArray;


class SQGFRaveVertexing: public SubsysReco
{
public:
  SQGFRaveVertexing(const std::string& name = "SQGFRaveVertexing");
  ~SQGFRaveVertexing();

  int Init(PHCompositeNode* topNode);
  int InitRun(PHCompositeNode* topNode);
  int process_event(PHCompositeNode* topNode);
  int End(PHCompositeNode* topNode);

  const std::string& get_vertexing_method() const
  {
    return _vertexing_method;
  }

  void set_vertexing_method(const std::string& vertexingMethod)
  {
    _vertexing_method = vertexingMethod;
  }



  const std::string& get_geom_file_name() const { return _geom_file_name; }
  void set_geom_file_name(const std::string& geomFileName) { _geom_file_name = geomFileName; }

  void set_legacy_rec_container(const bool enable = true) { legacyContainer = enable; }
  void set_vtx_smearing(const double r) { vtxSmearing = r; }

private:
  int MakeNodes(PHCompositeNode* topNode);
  int GetNodes(PHCompositeNode* topNode);
  int InitField(PHCompositeNode* topNode);
  int InitGeom(PHCompositeNode* topNode);

  bool buildRecDimuon(double z_vtx, SRecTrack* posTrack, SRecTrack* negTrack, SRecDimuon* dimuon);
  double swimTrackToVertex(SRecTrack* track, double z, TVector3* pos = nullptr, TVector3* mom = nullptr);

  TRandom1 rndm;

  bool legacyContainer;
  double vtxSmearing;

  SRecEvent*      recEvent;
  SQTrackVector*  recTrackVec;
  SQTrackVector*  truthTrackVec;
  SQDimuonVector* truthDimuonVec;

  SQDimuonVector* recDimuonVec;


//Abi add
 genfit::GFRaveVertexFactory* _vertex_finder;
 //! https://rave.hepforge.org/trac/wiki/RaveMethods
 std::string _vertexing_method;

  //const std::vector<genfit::GFRaveVertex*>& rave_vertices,
  SQGenFit::GFFitter* _gfitter;

  PHField* _phfield;
  SQGenFit::GFField* _gfield;
   recoConsts* rc;


  std::string  _geom_file_name;
  TGeoManager* _t_geo_manager;

  genfit::Track* TranslateSRecTrackToGenFitTrack(SRecTrack* srec_track);


};

#endif
