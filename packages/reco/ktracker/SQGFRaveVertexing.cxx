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

using namespace std;

namespace 
{
  //static flag to indicate the initialized has been done
  static bool inited = false;

  //static flag of kmag strength
  static double FMAGSTR;
  static double KMAGSTR;

  //Beam position and shape
  static double X_BEAM;
  static double Y_BEAM;
  static double SIGX_BEAM;
  static double SIGY_BEAM;

  //Simple swimming settings 
  static int NSTEPS_TARGET;
  static int NSTEPS_SHIELDING;
  static int NSTEPS_FMAG;

  //Geometric constants
  static double Z_TARGET;
  static double Z_DUMP;
  static double Z_UPSTREAM;
  static double Z_DOWNSTREAM;

  //initialize global variables
  void initGlobalVariables()
  {
    if(!inited) 
      {
	inited = true;

	recoConsts* rc = recoConsts::instance();
	FMAGSTR = rc->get_DoubleFlag("FMAGSTR");
	KMAGSTR = rc->get_DoubleFlag("KMAGSTR");
            
	X_BEAM = rc->get_DoubleFlag("X_BEAM");
	Y_BEAM = rc->get_DoubleFlag("Y_BEAM");
	SIGX_BEAM = rc->get_DoubleFlag("SIGX_BEAM");
	SIGY_BEAM = rc->get_DoubleFlag("SIGY_BEAM");

	NSTEPS_TARGET = rc->get_IntFlag("NSTEPS_TARGET");
	NSTEPS_SHIELDING = rc->get_IntFlag("NSTEPS_SHIELDING");
	NSTEPS_FMAG = rc->get_IntFlag("NSTEPS_FMAG");

	Z_TARGET = rc->get_DoubleFlag("Z_TARGET");
	Z_DUMP = rc->get_DoubleFlag("Z_DUMP");
	Z_UPSTREAM = rc->get_DoubleFlag("Z_UPSTREAM");
	Z_DOWNSTREAM = rc->get_DoubleFlag("Z_DOWNSTREAM");
      }
  }
}

SQGFRaveVertexing::SQGFRaveVertexing(const std::string& name):
  SubsysReco(name),
  legacyContainer(true),
  vtxSmearing(-1.),
  recEvent(nullptr),
  recTrackVec(nullptr),
  truthTrackVec(nullptr),
  truthDimuonVec(nullptr),
  _vertex_finder(nullptr),
  _vertexing_method("kalman-smoothing:1"),
  //_vertexing_method("avr-smoothing:1-minweight:0.5-primcut:9-seccut:9"),
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


  /* //get instance of fitter
     _fitter = PHGenFit::Fitter::getInstance(tgeo_manager,
     field, "DafRef",
     "RKTrackRep", false);
     _fitter->set_verbosity(Verbosity());

     if (!_fitter)
     {
     cerr << PHWHERE << endl;
     return Fun4AllReturnCodes::ABORTRUN;
     }
  */
  //create vertex factory 
  //_vertex_finder = new genfit::GFRaveVertexFactory(Verbosity());
  _vertex_finder = new genfit::GFRaveVertexFactory(1);
  //set the method to find vertex
  _vertex_finder->setMethod(_vertexing_method.data());
  
  TMatrixDSym beam_spot_cov(3);
  beam_spot_cov(0, 0) = SIGX_BEAM*SIGX_BEAM;;      // dx * dx
  beam_spot_cov(1, 1) = SIGY_BEAM*SIGY_BEAM;      // dy * dy
  beam_spot_cov(2, 2) = 100.0 * 100.0;  // dz * dz
  _vertex_finder->setBeamspot(TVector3(X_BEAM, Y_BEAM, -300.), beam_spot_cov);
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

  if(Verbosity() > Fun4AllBase::VERBOSITY_A_LOT) 
    {
      /* std::cout << "PHField check: " << "-------" << std::endl;
	 std::ofstream field_scan("field_scan.csv");
	 _phfield->identify(field_scan);
	 field_scan.close();
      */
    }

  //if(_fitter_type != SQGFRaveVertexing::LEGACY) 
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
  std::map<int, SRecTrack*> posTracks;
  std::map<int, SRecTrack*> negTracks;

  std::vector<genfit::Track*> gf_tracks;

  //loop over the tracks
  int nTracks = truthTrackVec->size();
  int nRecTracks = legacyContainer ? recEvent->getNTracks() : recTrackVec->size();
  for(int i = 0; i < nTracks; ++i)
    {
      SQTrack* trk = truthTrackVec->at(i);
      int recTrackIdx = trk->get_rec_track_id();
      if(recTrackIdx < 0 || recTrackIdx >= nRecTracks) continue;

      SRecTrack* recTrack = legacyContainer ? &(recEvent->getTrack(recTrackIdx)):dynamic_cast<SRecTrack*>(recTrackVec->at(recTrackIdx));


      if(recTrack->getCharge() > 0)
	{
	  posTracks[trk->get_track_id()] = recTrack;
	}
      else
	{
	  negTracks[trk->get_track_id()] = recTrack;
	}

       
      /*       SQGenFit::GFTrack track(*recTrack);         
	       genfit::Track* genfit_track = track.getGenFitTrack();                   
        
      */
      auto genfit_track = TranslateSRecTrackToGenFitTrack(recTrack);

      if (!genfit_track) continue;

      if (Verbosity())
	{
	  TVector3 pos, mom;
	  TMatrixDSym cov;

	  genfit_track->getFittedState().getPosMomCov(pos, mom, cov);

	  cout << "Track getCharge = " << genfit_track->getFitStatus()->getCharge() << " getChi2 = " << genfit_track->getFitStatus()->getChi2() << " getNdf = " << genfit_track->getFitStatus()->getNdf() << endl;
	  pos.Print();
	  mom.Print();
	  cov.Print();
	}
            
      gf_tracks.push_back(const_cast<genfit::Track*>(genfit_track));
 	
    }
     
  std::cout<<"GFRaveVertex: track size: "<<gf_tracks.size()<<std::endl;

  std::vector<genfit::GFRaveVertex*> rave_vertices;
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

  TVector3 *abi_vtx = new TVector3(9999.,9999.,9999.);
  std::cout<<"GFRaveVertex: vertex size: "<<rave_vertices.size()<<std::endl;
  for (unsigned int i=0; i<rave_vertices.size(); ++i) {

    genfit::GFRaveVertex* rvtx = static_cast<genfit::GFRaveVertex*>(rave_vertices[i]);
    
    std::cout<<"(X, Y, Z): ("<<rvtx->getPos().X()<<", "<<rvtx->getPos().Y()<<", "<<rvtx->getPos().Z()<<")"<<std::endl;
    abi_vtx->SetXYZ(rvtx->getPos().X(),rvtx->getPos().Y(),rvtx->getPos().Z());
   
  }


  if(posTracks.empty() || negTracks.empty()) return Fun4AllReturnCodes::EVENT_OK;
  if(truthDimuonVec)
    {
      int nTrueDimuons = truthDimuonVec->size();
      for(int i = 0; i < nTrueDimuons; ++i)
	{
	  SQDimuon* trueDimuon = truthDimuonVec->at(i);
	  int pid = trueDimuon->get_track_id_pos();
	  int mid = trueDimuon->get_track_id_neg();
	  if(posTracks.find(pid) == posTracks.end() || negTracks.find(mid) == negTracks.end()) continue;

	  SRecDimuon recDimuon;
	  //double z_vtx = trueDimuon->get_pos().Z() + (vtxSmearing>0. ? rndm.Gaus(0., vtxSmearing) : 0.);
	  //if(!buildRecDimuon(z_vtx, posTracks[pid], negTracks[mid], &recDimuon)) continue;
	  recDimuon.trackID_pos = pid;
	  recDimuon.trackID_neg = mid;

	  if(rave_vertices.size()==1)recDimuon.vtx.SetXYZ(abi_vtx->X(), abi_vtx->Y(), abi_vtx->Z());
	  else recDimuon.vtx.SetXYZ(9999., 9999., 9999.);
	  std::cout<<"--- Dimuon Z: "<<recDimuon.vtx.Z()<<std::endl;
	  trueDimuon->set_rec_dimuon_id(legacyContainer ? recEvent->getNDimuons() : recDimuonVec->size());
	  if(!legacyContainer)
	    recDimuonVec->push_back(&recDimuon);
	  else
	    recEvent->insertDimuon(recDimuon);
	}
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
  
      SQGenFit::GFTrack track(*srec_track);

      TMatrixDSym seed_cov = track.getGenFitTrack()->getCovSeed(); 
      TVectorD seed_st = track.getGenFitTrack()->getStateSeed();    

      int _pdg = srec_track->getCharge() > 0 ? -13 : 13;
      genfit::AbsTrackRep* rep = new genfit::RKTrackRep(_pdg);

      genfit::Track* genfit_track = new genfit::Track(rep, seed_st , seed_cov );

      genfit::FitStatus* fs = new genfit::FitStatus();
      fs->setCharge(srec_track->get_charge());
      fs->setChi2(srec_track->getChisq());
      fs->setNdf(srec_track->getNHits()-5);
      fs->setIsFitted(true);
      fs->setIsFitConvergedFully(true);
      genfit_track->setFitStatus(fs, rep);

      genfit::TrackPoint* tp = new genfit::TrackPoint(genfit_track);

      genfit::KalmanFitterInfo* fi = new genfit::KalmanFitterInfo(tp, rep);
      tp->setFitterInfo(fi);

      int nHits = srec_track->getNHits();
      std::vector<genfit::MeasuredStateOnPlane> _fitstates;
      for(int i = 0; i < nHits; ++i){
	genfit::SharedPlanePtr detPlane(new genfit::DetPlane(srec_track->getGFPlaneO(i), srec_track->getGFPlaneU(i), srec_track->getGFPlaneV(i)));
	_fitstates.push_back(genfit::MeasuredStateOnPlane(srec_track->getGFState(i),srec_track->getGFCov(i),detPlane,rep,srec_track->getGFAuxInfo(i)));
   
      }
 
      genfit::KalmanFittedStateOnPlane* kfs = new genfit::KalmanFittedStateOnPlane(_fitstates.front(), 1., 1.);
      
      //< Acording to the special order of using the stored states
      fi->setForwardUpdate(kfs);
      //fi->setBackwardUpdate(kfs);
      genfit_track->insertPoint(tp);

      return genfit_track;
    }
  catch (...)
    {
      std::cout << PHWHERE <<"TranslateSRecTrackToGenFitTrack failed!";
    }

  return nullptr;

}




