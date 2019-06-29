#ifndef __DECO_PARAM_H__
#define __DECO_PARAM_H__
#include <chan_map/ChanMapTaiwan.h>
#include <chan_map/ChanMapV1495.h>
#include <chan_map/ChanMapScaler.h>
//#include "assert.h"
class EventInfo;

/** This class contains all decoder parameters.
 * Two types of parameters are included.
 *  - For controlling the decoder behavior.
 *  - For monitoring the condition of the decoder.
 */
struct DecoParam {
  ///
  /// Control parameters
  ///
  std::string fn_in;
  std::string dir_param;
  int sampling;
  int verbose;
  int time_wait; //< waiting time in second to pretend the online data flow.

  ChanMapTaiwan chan_map_taiwan;
  ChanMapV1495  chan_map_v1495;
  ChanMapScaler chan_map_scaler;

  ///
  /// Monitoring parameters
  ///
  unsigned short int runID;

  int spillID;
  int spillID_cntr; // spillID from spill-counter event
  int spillID_slow; // spillID from slow-control event
  short targPos;
  short targPos_slow; // from slow-control event

  unsigned int codaID; //< current Coda event ID
  short spillType; //< current spill type
  short rocID; //< current ROC ID
  int eventIDstd; //< current event ID of standard physics events
  long int hitID; //< current hit ID, commonly used by Hit and TriggerHit.

  bool has_1st_bos; //< Set 'true' at the 1st BOS event.
  bool at_bos; //< Set 'true' at BOS, which triggers parsing the data in one spill.

  typedef enum {
    WORD_ONLY89     = 0x01,
    WORD_OVERFLOW   = 0x02,
    MULTIPLE_HEADER = 0x04,
    NO_EVT_ID       = 0x08,
    START_WO_STOP   = 0x10,
    START_NOT_RISE  = 0x20,
    DIRTY_FINISH    = 0x40
  } CodaPhysEvtStatus_t;
  int coda_phys_evt_status; //< Bits to hold the status of an Coda event.
  
  /// Max turnOnset in a spill.  Used to drop NIM3 events that came after the beam stops.  See elog 15010
  unsigned int turn_id_max;

  // Benchmarking Values
  time_t timeStart, timeEnd;
  
  DecoParam();
  ~DecoParam() {;}
  int InitMapper();
};

#endif // __DECO_PARAM_H__
