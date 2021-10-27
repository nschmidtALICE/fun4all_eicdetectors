#include "PHG4FoCalSubsystem.h"
#include "PHG4FoCalDetector.h"
#include "PHG4FoCalDisplayAction.h"
#include "PHG4FoCalSteppingAction.h"

#include <g4main/PHG4DisplayAction.h>       // for PHG4DisplayAction
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4SteppingAction.h>      // for PHG4SteppingAction
#include <g4main/PHG4Subsystem.h>           // for PHG4Subsystem
#include <g4main/PHG4Utils.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>             // for PHIODataNode
#include <phool/PHNode.h>                   // for PHNode
#include <phool/PHNodeIterator.h>           // for PHNodeIterator
#include <phool/PHObject.h>                 // for PHObject
#include <phool/getClass.h>

#include <set>                              // for set
#include <sstream>

class PHG4Detector;

using namespace std;

//_______________________________________________________________________
PHG4FoCalSubsystem::PHG4FoCalSubsystem(const std::string& name, const int lyr)
  : PHG4DetectorSubsystem(name)
  , m_Detector(nullptr)
  , m_SteppingAction(nullptr)
  , m_DisplayAction(nullptr)
  , active(1)
  , absorber_active(0)
  , blackhole(0)
  , tungstenpreshower(false)
  , maxlengthtower(false)
  , fullgeodet(false)
  , fillhollow(false)
  , detector_type(name)
  , mappingfile_("")
{
}

//_______________________________________________________________________
PHG4FoCalSubsystem::~PHG4FoCalSubsystem()
{
  delete m_DisplayAction;
}

//_______________________________________________________________________
int PHG4FoCalSubsystem::InitRunSubsystem(PHCompositeNode* topNode)
{
  PHNodeIterator iter(topNode);
  PHCompositeNode* dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));

  // create display settings before detector
  m_DisplayAction = new PHG4FoCalDisplayAction(Name());
  // create detector
  m_Detector = new PHG4FoCalDetector(this, topNode, Name());
  m_Detector->SetActive(active);
  m_Detector->SetAbsorberActive(absorber_active);
  m_Detector->BlackHole(blackhole);
  m_Detector->UseMaxLength(maxlengthtower);
  m_Detector->UseFullGeometry(fullgeodet);
  m_Detector->FillHollowSpace(fillhollow);
  m_Detector->OverlapCheck(CheckOverlap());
  m_Detector->Verbosity(Verbosity());
  m_Detector->SetTowerMappingFile(mappingfile_);
  m_Detector->SuperDetector(SuperDetector());
  m_Detector->SetUsePreshowerTungsten(tungstenpreshower);





  std::set<std::string> nodes;
  if (active)
  {
    PHNodeIterator dstIter(dstNode);
    PHCompositeNode* DetNode = dstNode;
    if (SuperDetector() != "NONE")
    {
      DetNode = dynamic_cast<PHCompositeNode*>(dstIter.findFirst("PHCompositeNode", SuperDetector()));
      if (!DetNode)
      {
        DetNode = new PHCompositeNode(SuperDetector());
        dstNode->addNode(DetNode);
      }
    }
    // create hit output node
    std::string nodename;
    if (SuperDetector() != "NONE")
    {
      nodename = "G4HIT_" + SuperDetector();
    }
    else
    {
      nodename = "G4HIT_" + Name();
    }
    nodes.insert(nodename);
    if (absorber_active)
    {
      if (SuperDetector() != "NONE")
      {
        nodename = "G4HIT_ABSORBER_" + SuperDetector();
      }
      else
      {
        nodename = "G4HIT_ABSORBER_" + Name();
      }
      nodes.insert(nodename);
    }
    for (auto thisnode : nodes)
    {
      PHG4HitContainer* g4_hits = findNode::getClass<PHG4HitContainer>(topNode, thisnode);
      if (!g4_hits)
      {
        g4_hits = new PHG4HitContainer(thisnode);
        DetNode->addNode(new PHIODataNode<PHObject>(g4_hits, thisnode, "PHObject"));
      }
    }
    // create stepping action
    m_SteppingAction = new PHG4FoCalSteppingAction(m_Detector,absorber_active);
    m_SteppingAction->Verbosity(0);
    m_Detector->SetSteppingAction(dynamic_cast<PHG4FoCalSteppingAction*>(m_SteppingAction));
  }

  return 0;
}

//_______________________________________________________________________
int PHG4FoCalSubsystem::process_event(PHCompositeNode* topNode)
{
  // pass top node to stepping action so that it gets
  // relevant nodes needed internally
  if (m_SteppingAction)
  {
    m_SteppingAction->SetInterfacePointers(topNode);
  }
  return 0;
}

//_______________________________________________________________________
PHG4Detector* PHG4FoCalSubsystem::GetDetector() const
{
  return m_Detector;
}


void PHG4FoCalSubsystem::SetDefaultParameters()
{
  set_default_double_param("place_x", 0.);
  set_default_double_param("place_y", 0.);
  set_default_double_param("place_z", 375.);
  set_default_double_param("tower_dx", 0.3);
  set_default_double_param("tower_dy", 0.3);
  set_default_double_param("tower_dz", 150.);
  set_default_double_param("dz", 150.);
  set_default_double_param("rMin1", 20.);
  set_default_double_param("rMax1", 220.);
  set_default_double_param("rMin2", 20.);
  set_default_double_param("rMax2", 220.);
  // set_default_double_param("wls_dw", 0.3);
  // set_default_double_param("support_dw", 0.2);
  set_default_double_param("rot_x", 0.);
  set_default_double_param("rot_y", 0.);
  set_default_double_param("rot_z", 0.);
  // set_default_double_param("thickness_absorber", 2.);
  // set_default_double_param("thickness_scintillator", 0.231);
  // set_default_string_param("scintillator", "G4_POLYSTYRENE");
  // set_default_string_param("absorber", "G4_Fe");
  // set_default_string_param("support", "G4_Fe");
  return;
}


void PHG4FoCalSubsystem::SetTowerMappingFile(const std::string& filename)
{
  // set_string_param("mapping_file", filename);
  // set_string_param("mapping_file_md5", PHG4Utils::md5sum(get_string_param("mapping_file")));
  mappingfile_ = filename;
}