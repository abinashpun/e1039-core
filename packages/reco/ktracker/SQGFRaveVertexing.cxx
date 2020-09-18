#include "SQGFRaveVertexing.h"

#include <phfield/PHFieldConfig_v3.h>
#include <phfield/PHFieldUtility.h>
#include <phgeom/PHGeomUtility.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHIODataNode.h>
#include <phool/PHRandomSeed.h>
#include <phool/getClass.h>
#include <phool/recoConsts.h>
#include <phgeom/PHGeomTGeo.h>
#include <interface_main/SQTrack_v1.h>
#include <interface_main/SQDimuon_v1.h>
#include <interface_main/SQTrackVector_v1.h>
#include <interface_main/SQDimuonVector_v1.h>


#include <GenFit/FitStatus.h>                     // for FitStatus
#include <GenFit/GFRaveTrackParameters.h>         // for GFRaveTrackParameters
#include <GenFit/GFRaveVertex.h>
#include <GenFit/GFRaveVertexFactory.h>
#include <GenFit/KalmanFittedStateOnPlane.h>      // for KalmanFittedStateOn...
#include <GenFit/KalmanFitterInfo.h>
#include <GenFit/MeasuredStateOnPlane.h>
#include <GenFit/RKTrackRep.h>
#include <GenFit/Track.h>
#include <GenFit/TrackPoint.h>                    // for TrackPoint

#include <TMatrixDSymfwd.h>                       // for TMatrixDSym
#include <TMatrixTSym.h>                          // for TMatrixTSym
#include <TMatrixTUtils.h>                        // for TMatrixTRow
#include <TVector3.h>

#include "SRecEvent.h"
#include "GFTrack.h"

#include <iostream>
#include <map>
#include <memory>
#include <utility>
#include <vector>

SQGFRaveVertexing::SQGFRaveVertexing(const std::string& name):
  SubsysReco(name),
  legacyContainer(false),
  vtxSmearing(-1.),
  recEvent(nullptr),
  recTrackVec(nullptr),
  truthTrackVec(nullptr),
  truthDimuonVec(nullptr),
  _vertex_finder(nullptr),
  //_vertexing_method("kalman-smoothing:1"),
  _vertexing_method("avr-smoothing:1-minweight:0.5-primcut:9-seccut:9"),
  recDimuonVec(nullptr)
{
 rc = recoConsts::instance();
}

SQGFRaveVertexing::~SQGFRaveVertexing()
{
 delete _vertex_finder;
}

int SQGFRaveVertexing::Init(PHCompositeNode* topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int SQGFRaveVertexing::InitRun(PHCompositeNode* topNode)
{
  int ret = GetNodes(topNode);
  if(ret != Fun4AllReturnCodes::EVENT_OK) return ret;

  ret = MakeNodes(topNode);
  if(ret != Fun4AllReturnCodes::EVENT_OK) return ret;

  ret = InitField(topNode);
  if(ret != Fun4AllReturnCodes::EVENT_OK) return ret;

  ret = InitGeom(topNode);
  if(ret != Fun4AllReturnCodes::EVENT_OK) return ret;

  if(vtxSmearing) rndm.SetSeed(PHRandomSeed());

  //create vertex factory 
  //_vertex_finder = new genfit::GFRaveVertexFactory(Verbosity());
  _vertex_finder = new genfit::GFRaveVertexFactory(1);
 //set the method to find vertex
  _vertex_finder->setMethod(_vertexing_method.data());
  //_vertex_finder->setBeamspot();

  if (!_vertex_finder)
  {
    std::cerr << PHWHERE << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }


  return Fun4AllReturnCodes::EVENT_OK;
}



int SQGFRaveVertexing::InitField(PHCompositeNode* topNode)
{
  if(Verbosity() > 1)
  {
    std::cout << "SQGFRaveVertexing::InitField" << std::endl;
 
  }

  std::unique_ptr<PHFieldConfig> default_field_cfg(new PHFieldConfig_v3(rc->get_CharFlag("fMagFile"), rc->get_CharFlag("kMagFile"), rc->get_DoubleFlag("FMAGSTR"), rc->get_DoubleFlag("KMAGSTR"), 5.));
  _phfield = PHFieldUtility::GetFieldMapNode(default_field_cfg.get(), topNode, 0);

   _gfield = new SQGenFit::GFField(_phfield);

  return Fun4AllReturnCodes::EVENT_OK;
}

int SQGFRaveVertexing::InitGeom(PHCompositeNode* topNode)
{
  if(Verbosity() > 1) 
  {
    std::cout << "SQGFRaveVertexing::InitGeom" << std::endl;
   
  }

  PHGeomTGeo* dstGeom = PHGeomUtility::GetGeomTGeoNode(topNode, true); //hacky way to bypass PHGeoUtility's lack of exception throwing
  if(!dstGeom->isValid())
  {
    if(_geom_file_name == "") return Fun4AllReturnCodes::ABORTEVENT;

    if(Verbosity() > 1) std::cout << "SQGFRaveVertexing::InitGeom - create geom from " << _geom_file_name << std::endl;
    int ret = PHGeomUtility::ImportGeomFile(topNode, _geom_file_name);
    if(ret != Fun4AllReturnCodes::EVENT_OK) return ret;
  }
  else
  {
    if(Verbosity() > 1) std::cout << "SQGFRaveVertexing::InitGeom - use geom from NodeTree." << std::endl;
  }

  _t_geo_manager = PHGeomUtility::GetTGeoManager(topNode);
  return Fun4AllReturnCodes::EVENT_OK;
}




int SQGFRaveVertexing::process_event(PHCompositeNode* topNode)
{


  std::vector<genfit::Track*> gf_tracks;
 //loop over the tracks
   
  int nRecTracks = legacyContainer ? recEvent->getNTracks() : recTrackVec->size();
   for(int i = 0; i < RecTracks; ++i)
	{
  
        SRecTrack* recTrack = legacyContainer ? &(recEvent->getTrack(recTrackIdx)):dynamic_cast<SRecTrack*>(recTrackVec->at(recTrackIdx));

       ///translating SRecTrack to GenFit Track
         auto genfit_track = TranslateSRecTrackToGenFitTrack(recTrack);
        if (!genfit_track)
          continue;
    
         gf_tracks.push_back(const_cast<genfit::Track*>(genfit_track));
 	
	}

  std::cout<<"GFRaveVertex: track size: "<<gf_tracks.size()<<std::endl;

  std::vector<genfit::GFRaveVertex*> rave_vertices;


///Now find vertex from GFRave
  if (gf_tracks.size() >= 2)
  {
    try
    {
      _vertex_finder->findVertices(&rave_vertices, gf_tracks);
    }
    catch (...)
    {
      //if (Verbosity() > 1)
       std::cout << PHWHERE << "GFRaveVertexFactory::findVertices failed!";
    }
  }

   std::cout<<"GFRaveVertex: vertex size: "<<rave_vertices.size()<<std::endl;

    for (unsigned int i=0; i<rave_vertices.size(); ++i) {  
      genfit::GFRaveVertex* vtx = static_cast<genfit::GFRaveVertex*>(rave_vertices[i]);
      std::cout<<"(X, Y, Z): ("<<vtx->getPos().X()<<", "vtx->getPos().Y()<<", "<<vtx->getPos().Z() <<std::endl;
     }

  return Fun4AllReturnCodes::EVENT_OK;
}

int SQGFRaveVertexing::End(PHCompositeNode* topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int SQGFRaveVertexing::GetNodes(PHCompositeNode* topNode)
{
  truthTrackVec = findNode::getClass<SQTrackVector>(topNode, "SQTruthTrackVector");
  if(!truthTrackVec)
  {
    std::cerr << Name() << ": failed finding truth track info, abort." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  truthDimuonVec = findNode::getClass<SQDimuonVector>(topNode, "SQTruthDimuonVector");
  if(!truthDimuonVec)
  {
    std::cout << Name() << ": failed finding truth dimuon info, rec dimuon will be for reference only. " << std::endl;
  }

  if(legacyContainer)
  {
    recEvent = findNode::getClass<SRecEvent>(topNode, "SRecEvent");
  }
  else
  {
    recTrackVec = findNode::getClass<SQTrackVector>(topNode, "SQRecTrackVector");
  }

  if((!recEvent) && (!recTrackVec))
  {
    std::cerr << Name() << ": failed finding rec track info, abort." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int SQGFRaveVertexing::MakeNodes(PHCompositeNode* topNode)
{
  if(!legacyContainer)
  {
    PHNodeIterator iter(topNode);
    PHCompositeNode* dstNode = static_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
    if(!dstNode) 
    {
      std::cerr << Name() << ": cannot locate DST node, abort." << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }

    recDimuonVec = new SQDimuonVector_v1();
    dstNode->addNode(new PHIODataNode<PHObject>(recDimuonVec, "SQRecDimuonVector", "PHObject"));
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

double SQGFRaveVertexing::swimTrackToVertex(SRecTrack* track, double z, TVector3* pos, TVector3* mom)
{
  SQGenFit::GFTrack gftrk(*track);

  TVector3 p, m;
  double chi2 = gftrk.swimToVertex(z, &p, &m);
  if(chi2 < 0.) return chi2;

  if(pos == nullptr)
  {
    track->setChisqVertex(chi2);
    track->setVertexPos(p);
    track->setVertexMom(m);
  }
  else
  {
    pos->SetXYZ(p.X(), p.Y(), p.Z());
    mom->SetXYZ(m.X(), m.Y(), m.Z());
  }

  return chi2;
}

bool SQGFRaveVertexing::buildRecDimuon(double z_vtx, SRecTrack* posTrack, SRecTrack* negTrack, SRecDimuon* dimuon)
{
  TVector3 p_mom, p_pos = posTrack->getVertexPos();
  double p_chi2;
  if(fabs(p_pos.Z() - z_vtx) > 1.)
  {
    p_chi2 = swimTrackToVertex(posTrack, z_vtx, &p_pos, &p_mom);
    if(p_chi2 < 0.) return false;
  }
  else
  {
    p_mom = posTrack->getVertexMom();
    p_chi2 = posTrack->getChisqVertex();
  }

  TVector3 m_mom, m_pos = negTrack->getVertexPos();
  double m_chi2;
  if(fabs(m_pos.Z() - z_vtx) > 1.)
  {
    m_chi2 = swimTrackToVertex(negTrack, z_vtx, &p_pos, &p_mom);
    if(m_chi2 < 0.) return false;
  }
  else
  {
    m_mom = negTrack->getVertexMom();
    m_chi2 = negTrack->getChisqVertex();
  }

  dimuon->p_pos.SetVectM(p_mom, M_MU);
  dimuon->p_neg.SetVectM(m_mom, M_MU);
  dimuon->chisq_kf = p_chi2 + m_chi2;
  dimuon->vtx.SetXYZ(0., 0., z_vtx);
  dimuon->vtx_pos = p_pos;
  dimuon->vtx_neg = m_pos;
  dimuon->proj_target_pos = posTrack->getTargetPos();
  dimuon->proj_dump_pos = posTrack->getDumpPos();
  dimuon->proj_target_neg = negTrack->getTargetPos();
  dimuon->proj_dump_neg = negTrack->getDumpPos();
  dimuon->chisq_target = posTrack->getChisqTarget() + negTrack->getChisqTarget();
  dimuon->chisq_dump = posTrack->getChisqDump() + negTrack->getChisqDump();
  dimuon->chisq_upstream = posTrack->getChisqUpstream() + negTrack->getChisqUpstream();
  dimuon->chisq_single = 0.;
  dimuon->chisq_vx = 0.;
  dimuon->calcVariables();

  return true;
}


genfit::Track* SQGFRaveVertexing::TranslateSRecTrackToGenFitTrack(SRecTrack* srec_track)
{

   try
  {

    TVector3 pos = srec_track->getVertexPos();//look at the node [0] for both pos and momentum
    TVector3 mom = srec_track->getVertexMom();
    TMatrixDSym cov(6);// = srec_track->getGFCov(0);//get it from srec_track

      double uncertainty[6] = {10., 10., 10., 3., 3., 10.};
   for(int i = 0; i < 6; i++)
   {
    for(int j = 0; j < 6; j++)
    { 
       //std::cout <<"cov....."<<srec_track->getGFCov(0)[i][j] << "  ";
      cov[i][j] = uncertainty[i]*uncertainty[j];
     }
     std::cout << std::endl;
   }

    int _pdg = srec_track->getCharge() > 0 ? -13 : 13;
    genfit::AbsTrackRep* rep = new genfit::RKTrackRep(_pdg);
    genfit::Track* genfit_track = new genfit::Track(rep, TVector3(0, 0, 0), TVector3(0, 0, 0));


    genfit::FitStatus* fs = new genfit::FitStatus();
    fs->setCharge(srec_track->get_charge());
    fs->setChi2(srec_track->getChisqVertex());
    fs->setNdf(srec_track->getNHits()-5);
    fs->setIsFitted(true);
    fs->setIsFitConvergedFully(true);

    genfit_track->setFitStatus(fs, rep);

    genfit::TrackPoint* tp = new genfit::TrackPoint(genfit_track);

    genfit::KalmanFitterInfo* fi = new genfit::KalmanFitterInfo(tp, rep);
    tp->setFitterInfo(fi);

    genfit::MeasuredStateOnPlane* ms = new genfit::MeasuredStateOnPlane(rep);
    ms->setPosMomCov(pos, mom, cov);
    
    genfit::KalmanFittedStateOnPlane* kfs = new genfit::KalmanFittedStateOnPlane(*ms, 1., 1.);
  
     //< Acording to the special order of using the stored states
    fi->setForwardUpdate(kfs);

    genfit_track->insertPoint(tp);
    return genfit_track;
  }
  catch (...)
  {
    std::cout << PHWHERE <<"TranslateSRecTrackToGenFitTrack failed!";
  }

  return nullptr;

}


