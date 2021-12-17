#include "PHG4FoCalDetector.h"
#include "PHG4FoCalDisplayAction.h"
#include "PHG4FoCalSteppingAction.h"

#include <g4main/PHG4Detector.h>       // for PHG4Detector
#include <g4main/PHG4DisplayAction.h>  // for PHG4DisplayAction
#include <g4main/PHG4Subsystem.h>

#include <Geant4/G4Box.hh>
#include <Geant4/G4Tubs.hh>
#include <Geant4/G4SubtractionSolid.hh>
#include <Geant4/G4Cons.hh>
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4PVPlacement.hh>
#include <Geant4/G4RotationMatrix.hh>  // for G4RotationMatrix
#include <Geant4/G4String.hh>              // for G4String
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4ThreeVector.hh>      // for G4ThreeVector
#include <Geant4/G4Transform3D.hh>      // for G4Transform3D
#include <Geant4/G4Types.hh>               // for G4double, G4int
#include <Geant4/G4Torus.hh>               // for G4double, G4int
#include <Geant4/G4Para.hh>               // for G4double, G4int
#include <Geant4/G4VPhysicalVolume.hh>  // for G4VPhysicalVolume
#include <Geant4/G4PVParameterised.hh>
#include <Geant4/G4PVReplica.hh>
#include <Geant4/G4NistManager.hh>
#include <Geant4/G4OpticalSurface.hh>
#include <Geant4/G4LogicalSkinSurface.hh>
#include <Geant4/G4GeometryTolerance.hh>
#include <Geant4/G4LogicalBorderSurface.hh>

#include <TSystem.h>

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <utility>  // for pair, make_pair

class G4VSolid;
class PHCompositeNode;

using namespace std;

//_______________________________________________________________________
PHG4FoCalDetector::PHG4FoCalDetector(PHG4Subsystem* subsys, PHCompositeNode* Node, const std::string& dnam)
  : PHG4Detector(subsys, Node, dnam)
  , m_DisplayAction(dynamic_cast<PHG4FoCalDisplayAction*>(subsys->GetDisplayAction()))
  , m_SteppingAction(0)
  , _place_in_x(0.0 * mm)
  , _place_in_y(0.0 * mm)
  , _place_in_z(4000.0 * mm)
  , _center_offset_x(0.0 * mm)
  , _center_offset_y(0.0 * mm)
  , _quadratic_detector(0)
  , _rot_in_x(0.0)
  , _rot_in_y(0.0)
  , _rot_in_z(0.0)
  , _rMin1(50 * mm)
  , _rMax1(2620 * mm)
  , _rMin2(50 * mm)
  , _rMax2(3369 * mm)
  , _dZ(1000 * mm)
  , _sPhi(0)
  , _dPhi(2 * M_PIl)
  , _tower_type(0)
  , _tower_readout(0.5 * mm)
  , _tower_dx(100 * mm)
  , _tower_dy(100 * mm)
  , _tower_dz(1000.0 * mm)
  , _scintFiber_diam(1.0 * mm)
  , _cerenkovFiber_diam(1.0 * mm)
  , _cerenkovFiber_material(0)
  , _tower_makeNotched(0)
  , _absorber_Material(0)
  , _materialScintillator("G4_POLYSTYRENE")
  , _materialAbsorber("G4_Fe")
  , _active(1)
  , _absorberactive(0)
  , _preshowertungstenactive(false)
  , _usemaxlength(false)
  , _usefullgeom(false)
  , _fillhollow(false)
  , _layer(0)
  , _blackhole(0)
  , _towerlogicnameprefix("hfocalTower")
  , _superdetector("NONE")
  , _mapping_tower_file("")
  , m_doLightProp(true)
{
}
//_______________________________________________________________________
int PHG4FoCalDetector::IsInFoCal(G4VPhysicalVolume* volume) const
{
  if (volume->GetName().find(_towerlogicnameprefix) != string::npos)
  {
    if (volume->GetName().find("scintillator") != string::npos)
    {
      if (_active)
        return 1;
      else
        return 0;
    }
    else if (volume->GetName().find("cherenkov") != string::npos)
    {
      if (_active)
        return 1;
      else
        return 0;
    }
    //only record energy in actual absorber- drop energy lost in air gaps inside focal envelope
    else if (volume->GetName().find("absorber") != string::npos)
    {
      if (_absorberactive)
        return -1;
      else
        return 0;
    }
    else if (volume->GetName().find("envelope") != string::npos)
    {
      return 0;
    }
  }

  return 0;
}

//_______________________________________________________________________
void PHG4FoCalDetector::ConstructMe(G4LogicalVolume* logicWorld)
{
  // Verbosity(2);
  if (Verbosity() > 0)
  {
    cout << "PHG4FoCalDetector: Begin Construction" << endl;
  }

  //Read parameters for detector construction from file
  ParseParametersFromTable();

  //Create the cone envelope = 'world volume' for the calorimeter
  G4Material* Air = G4Material::GetMaterial("G4_AIR");

  int nPlates_W = 20;
  G4double thickness_W = 3.5*mm;
  G4double spacing_W = 2.5*mm;
  G4double AirPad = 0.6*mm;
  G4double BackplaneThickness = 1.0*mm;
  G4double AirGapFromFoCalE = 5.0*mm;
  G4double BoxThickness = 1.0*mm;
  G4double AirGapToFoCalH = 10.0*mm;
  G4double focal_E_dz = nPlates_W*thickness_W+nPlates_W*spacing_W+AirPad+BackplaneThickness+AirGapFromFoCalE+BoxThickness+AirGapToFoCalH;

  G4VSolid* focalE_envelope_solid = new G4Box("hfocalE_envelope_solid_precut",
                                        _rMax1,
                                        _rMax1,
                                        (focal_E_dz) / 2.0);

  G4LogicalVolume* focalE_envelope_log = new G4LogicalVolume(focalE_envelope_solid, Air, G4String("hfocalE_envelope"), 0, 0, 0);

  m_DisplayAction->AddVolume(focalE_envelope_log, "FfocalEnvelope");

  //Define rotation attributes for envelope cone
  G4RotationMatrix *focal_rotm = new G4RotationMatrix();
  focal_rotm->rotateX(_rot_in_x);
  focal_rotm->rotateY(_rot_in_y);
  focal_rotm->rotateZ(_rot_in_z);

  //Place envelope cone in simulation
  ostringstream name_envelopeE;
  name_envelopeE.str("");
  name_envelopeE << _towerlogicnameprefix << "_envelopeE" << endl;

  new G4PVPlacement(focal_rotm, G4ThreeVector(_place_in_x, _place_in_y, _place_in_z-_tower_dz/2),
                    focalE_envelope_log, name_envelopeE.str().c_str(), logicWorld, 0, false, OverlapCheck());

  if(_preshowertungstenactive){
    if(_quadratic_detector){
      // // focal E tungsten plates
      G4VSolid* focal_E_plate_cutout_solid = new G4Box("focal_E_plate_cutout_solid", _rMin1+0.01*cm,  _rMin1+0.01*cm, (thickness_W));
      G4VSolid* focal_E_plate_solid = new G4Box("focal_E_plate_solid_tmp", _rMax1,  _rMax1, (thickness_W) / 2.0);
      if(_rMin1>0){
        focal_E_plate_cutout_solid = new G4Box("focal_E_plate_cutout_solid", _rMin1,  _rMin1, (thickness_W));
        focal_E_plate_solid = new G4SubtractionSolid(G4String("focal_E_plate_solid"), focal_E_plate_solid, focal_E_plate_cutout_solid, 0 ,G4ThreeVector(0 ,0.,0)); // top right
      }
      G4LogicalVolume* focal_E_plate_log = new G4LogicalVolume(focal_E_plate_solid, G4Material::GetMaterial("G4_W"), G4String("focal_E_plate_log"), 0, 0, 0);
      for(int iplate=0;iplate<nPlates_W;iplate++){
        new G4PVPlacement(0, G4ThreeVector(0, 0, -(focal_E_dz)/2+thickness_W/2+iplate*thickness_W+iplate*spacing_W),
                            focal_E_plate_log,"focal_E_plate_physical_"+std::to_string(iplate), focalE_envelope_log, 0, false, OverlapCheck());
      }
      m_DisplayAction->AddVolume(focal_E_plate_log, "Tungsten");

      G4VSolid* focal_E_backplane_solid = new G4Box("focal_E_backplane_solid_tmp", _rMax1,  _rMax1, (BackplaneThickness) / 2.0);
      if(_rMin1>0){
        focal_E_backplane_solid = new G4SubtractionSolid(G4String("focal_E_plate_solid"), focal_E_backplane_solid, focal_E_plate_cutout_solid, 0 ,G4ThreeVector(0 ,0.,0)); // top right
      }
      G4LogicalVolume* focal_E_backplane_log = new G4LogicalVolume(focal_E_backplane_solid, G4Material::GetMaterial("G4_Al"), G4String("focal_E_backplane_log"), 0, 0, 0); 
      new G4PVPlacement(0, G4ThreeVector(0, 0, -(focal_E_dz)/2+thickness_W/2+20*thickness_W+20*spacing_W+AirPad+BackplaneThickness/2),
                            focal_E_backplane_log,"focal_E_backplane_physical", focalE_envelope_log, 0, false, OverlapCheck());
      new G4PVPlacement(0, G4ThreeVector(0, 0, -(focal_E_dz)/2+thickness_W/2+20*thickness_W+20*spacing_W+AirPad+BackplaneThickness+AirGapFromFoCalE+BoxThickness/2),
                            focal_E_backplane_log,"focal_E_backplane_physical", focalE_envelope_log, 0, false, OverlapCheck());
      m_DisplayAction->AddVolume(focal_E_backplane_log, "Aluminium");
    } else {
      // // focal E tungsten plates
      G4VSolid* focal_E_plate_solid = new G4Cons("focal_E_plate_solid",
                                        _rMin1, _rMax1,
                                        _rMin1, _rMax1,
                                        thickness_W / 2.0,
                                        0, 2 * M_PIl);
      G4LogicalVolume* focal_E_plate_log = new G4LogicalVolume(focal_E_plate_solid, G4Material::GetMaterial("G4_W"), G4String("focal_E_plate_log"), 0, 0, 0);
      for(int iplate=0;iplate<nPlates_W;iplate++){
        new G4PVPlacement(0, G4ThreeVector(0, 0, -(focal_E_dz)/2+thickness_W/2+iplate*thickness_W+iplate*spacing_W),
                            focal_E_plate_log,"focal_E_plate_physical_"+std::to_string(iplate), focalE_envelope_log, 0, false, OverlapCheck());
      }
      m_DisplayAction->AddVolume(focal_E_plate_log, "Tungsten");

      G4VSolid* focal_E_backplane_solid = new G4Cons("focal_E_backplane_solid",
                                        _rMin1, _rMax1,
                                        _rMin1, _rMax1,
                                        BackplaneThickness / 2.0,
                                        0, 2 * M_PIl);
      G4LogicalVolume* focal_E_backplane_log = new G4LogicalVolume(focal_E_backplane_solid, G4Material::GetMaterial("G4_Al"), G4String("focal_E_backplane_log"), 0, 0, 0); 
      new G4PVPlacement(0, G4ThreeVector(0, 0, -(focal_E_dz)/2+thickness_W/2+20*thickness_W+20*spacing_W+AirPad+BackplaneThickness/2),
                            focal_E_backplane_log,"focal_E_backplane_physical", focalE_envelope_log, 0, false, OverlapCheck());
      new G4PVPlacement(0, G4ThreeVector(0, 0, -(focal_E_dz)/2+thickness_W/2+20*thickness_W+20*spacing_W+AirPad+BackplaneThickness+AirGapFromFoCalE+BoxThickness/2),
                            focal_E_backplane_log,"focal_E_backplane_physical", focalE_envelope_log, 0, false, OverlapCheck());
      m_DisplayAction->AddVolume(focal_E_backplane_log, "Aluminium");


    }
  }
    G4VSolid* focalH_envelope_solid = new G4Box("hfocalH_envelope_solid_precut",
                                          _rMax1,
                                          _rMax1,
                                          (_tower_dz) / 2.0);

    G4LogicalVolume* focalH_envelope_log = new G4LogicalVolume(focalH_envelope_solid, Air, G4String("hfocalH_envelope"), 0, 0, 0);

    m_DisplayAction->AddVolume(focalH_envelope_log, "FfocalEnvelope");

    //Define rotation attributes for envelope cone
    G4RotationMatrix *focalH_rotm = new G4RotationMatrix();
    focalH_rotm->rotateX(_rot_in_x);
    focalH_rotm->rotateY(_rot_in_y);
    focalH_rotm->rotateZ(_rot_in_z);

    //Place envelope cone in simulation
    ostringstream name_envelopeH;
    name_envelopeH.str("");
    name_envelopeH << _towerlogicnameprefix << "_envelopeH" << endl;

    new G4PVPlacement(focalH_rotm, G4ThreeVector(_place_in_x, _place_in_y, _place_in_z+focal_E_dz/2),
                      focalH_envelope_log, name_envelopeH.str().c_str(), logicWorld, 0, false, OverlapCheck());



  if(_tower_type<10){
    G4LogicalVolume* singletower;
    if(_tower_type==6){
      singletower = ConstructTowerLayered(0);
    } else {
      singletower = ConstructTower(0);
    }

    G4Material* material_air = G4Material::GetMaterial("G4_AIR");
    // number of towers in radial direction (on y axis)
    int rowNtow = (int) ( (_rMax1-(_tower_dy/2)) / _tower_dy);
    for(int row=rowNtow;row>=-rowNtow;row--){
      // pythagoras -> get available length in circular mother volume for towers
      // divide given length by tower width -> get number of towers that can be placed
      int currRowNtow = (int) ( ( 2* sqrt(pow(_rMax1,2)-pow( (abs(row)*_tower_dy) ,2)) ) / _tower_dy );
      if(_quadratic_detector){
        currRowNtow = (int) ( (2*_rMax1) / _tower_dy);;
      }
      if(currRowNtow==0) continue;

      // we want an odd number of towers to be symmetrically centered around 0
      if ( currRowNtow % 2 == 0) currRowNtow-=1;

      if( ( (abs(row)*_tower_dy) + _tower_dy/2 ) < _rMin1){
        // pythagoras -> get available length in circular mother volume for towers
        // divide given length by tower width -> get number of towers that can be placed
        int currRowNtowInner = (int) ( ( 2* sqrt(pow(_rMin1,2)-pow( (abs(row)*_tower_dy) - (_tower_dy/2.0) ,2)) ) / _tower_dy );
        if(_quadratic_detector){
          currRowNtowInner = (int) ( (2*_rMin1) / _tower_dy);
        }
        // we want an odd number of towers to be symmetrically centered around 0
        if ( currRowNtowInner % 2 == 0) currRowNtowInner+=1;
        int currRowNtowMod = currRowNtow;
        if(_quadratic_detector){
          currRowNtowMod = 2*rowNtow + 1;
        }
        // create mother volume with space for currRowNtow towers along x-axis
        auto DRCalRowSolid    = new G4Box("DRCalRowBox" + std::to_string(row), (currRowNtowMod - currRowNtowInner) / 2 * _tower_dx / 2.0,_tower_dy / 2.0,_tower_dz / 2.0);
        auto DRCalRowLogical  = new G4LogicalVolume(DRCalRowSolid,material_air,"DRCalRowLogical" + std::to_string(row));
        // replicate singletower tower design currRowNtow times along x-axis
        new G4PVReplica("DRCalRowPhysical" + std::to_string(row),singletower,DRCalRowLogical,
                        kXAxis,(currRowNtowMod - currRowNtowInner) / 2,_tower_dx);

        ostringstream name_row_twr;
        name_row_twr.str("");
        name_row_twr << _towerlogicnameprefix << "_row_" << row << "_left" << endl;
        new G4PVPlacement(0, G4ThreeVector( - ( ( currRowNtowInner / 2.0 ) * _tower_dx ) - ( (currRowNtowMod - currRowNtowInner) / 2 * _tower_dx / 2.0 ), (row*_tower_dy), 0),
                      DRCalRowLogical, name_row_twr.str().c_str(), focalH_envelope_log, 0, false, OverlapCheck());

        ostringstream name_row_twr2;
        name_row_twr2.str("");
        name_row_twr2 << _towerlogicnameprefix << "_row_" << row << "_left" << endl;
        new G4PVPlacement(0, G4ThreeVector( ( ( currRowNtowInner / 2.0 ) * _tower_dx ) + ( (currRowNtowMod - currRowNtowInner) / 2 * _tower_dx / 2.0 ), (row*_tower_dy), 0),
                      DRCalRowLogical, name_row_twr2.str().c_str(), focalH_envelope_log, 0, false, OverlapCheck());
      } else {
        if(_quadratic_detector){
          // cout << currRowNtow << endl;
          // create mother volume with space for currRowNtow towers along x-axis
          auto DRCalRowSolid    = new G4Box("DRCalRowBox" + std::to_string(row), currRowNtow * _tower_dx / 2.0,_tower_dy / 2.0,_tower_dz / 2.0);
          auto DRCalRowLogical  = new G4LogicalVolume(DRCalRowSolid,material_air,"DRCalRowLogical");
          // replicate singletower tower design currRowNtow times along x-axis
          new G4PVReplica("DRCalRowPhysical" + std::to_string(row),singletower,DRCalRowLogical,
                          kXAxis,currRowNtow,_tower_dx);

          ostringstream name_row_twr;
          name_row_twr.str("");
          name_row_twr << _towerlogicnameprefix << "_row_" << row << endl;
          new G4PVPlacement(0, G4ThreeVector(0, (row*_tower_dy), 0),
                        DRCalRowLogical, name_row_twr.str().c_str(), focalH_envelope_log, 0, false, OverlapCheck());
        } else {
          // create mother volume with space for currRowNtow towers along x-axis
            // cout << currRowNtow << endl;
            auto DRCalRowSolid    = new G4Box("DRCalRowBox" + std::to_string(row), currRowNtow * _tower_dx / 2.0,_tower_dy / 2.0,_tower_dz / 2.0);
            auto DRCalRowLogical  = new G4LogicalVolume(DRCalRowSolid,material_air,"DRCalRowLogical");
            // replicate singletower tower design currRowNtow times along x-axis
            new G4PVReplica("DRCalRowPhysical" + std::to_string(row),singletower,DRCalRowLogical,
                            kXAxis,currRowNtow,_tower_dx);

            ostringstream name_row_twr;
            name_row_twr.str("");
            name_row_twr << _towerlogicnameprefix << "_row_" << row << endl;
            new G4PVPlacement(0, G4ThreeVector(0, (row*_tower_dy), 0),
                          DRCalRowLogical, name_row_twr.str().c_str(), focalH_envelope_log, 0, false, OverlapCheck());
        }
      }
    }
  } else if(_tower_type==69){
    ConstructCopperTiltedFibers(0,focalH_envelope_log);
  } else {
    ConstructCapillaryRowDetector(0,focalH_envelope_log);
    // place all 40 rows with 36 capillary tubes each
    // for(int irow=0;irow<40;irow++){
    //   G4ThreeVector vec_translation_xy(-_tower_dx/2 + (irow%2==0 ? - 2.5*mm / 4 : 2.5*mm / 4),-_tower_dy/2 +irow*(sqrt(3)*2.5*mm/2),_place_in_z+focal_E_dz/2);
    //   capillaryRow->MakeImprint( logicWorld, vec_translation_xy,0 );
    // }
  }
  return;
}

//_______________________________________________________________________
bool PHG4FoCalDetector::ConstructCopperTiltedFibers(int type,G4LogicalVolume* logic_envelope)
{
  G4double TowerDx = _rMax1*2;
  G4double TowerDz = _tower_dz;
  G4double fiber_spacing = 2*mm;
  G4double edge_distance = 5*mm;

  G4double fiber_thickness = 1.0*mm;
  G4double min_fiber_bending_radius = 1.25*cm;
  // provide an angle for the fiber loop tilt -> must be greater than Moliere spread
  G4double tilting_angle = 0.25 * M_PIl;
  // total height of the tilted loop for the given angle within the tower dimensions
  G4double loop_total_height = TowerDx / sin(tilting_angle);
  loop_total_height -= 2*sqrt(pow((fiber_thickness + fiber_spacing)/2.0,2)-pow(((fiber_thickness + fiber_spacing) * sin(tilting_angle))/2.0,2));
  // number of loops possible with the minimum bending radius
  int nLoopsFiber = (int) (loop_total_height/(2 * min_fiber_bending_radius));
  // odd number of loops required to have SiPMs on same side of tower
  if(nLoopsFiber%2 == 0) nLoopsFiber-=1;


  G4VSolid* solid_full_loop_para_mother = new G4Para("solid_full_loop_para_mother",
                                                (fiber_thickness + fiber_spacing) / 2.0, TowerDx / 2.0, TowerDx / 2.0,
                                                tilting_angle, 0, 0);
  G4double paraBoxWidth = (fiber_thickness + fiber_spacing) * sin(tilting_angle);

  G4VSolid* solid_full_loop_para = new G4Box("solid_full_loop_para",
                                                TowerDx / 2.0, loop_total_height / 2.0, paraBoxWidth / 2.0);
  G4LogicalVolume* logic_full_loop_para = new G4LogicalVolume(solid_full_loop_para,
                                                              G4Material::GetMaterial("G4_AIR"), "logic_full_loop_para",
                                                              0, 0, 0);
  m_DisplayAction->AddVolume(logic_full_loop_para, "Invisible");

  G4VSolid* solid_full_absorber_para_left = new G4Box("solid_full_absorber_para_left",
                                                (min_fiber_bending_radius + edge_distance) / 2.0, loop_total_height / 2.0, paraBoxWidth / 2.0);
  G4VSolid* solid_full_absorber_para_center = new G4Box("solid_full_absorber_para_center",
                                                (TowerDx - 2.0 * min_fiber_bending_radius - 2.0 * edge_distance) / 2.0, loop_total_height / 2.0, paraBoxWidth / 2.0);
  G4VSolid* solid_full_absorber_para_right = new G4Box("solid_full_absorber_para_right",
                                                (min_fiber_bending_radius + edge_distance) / 2.0, loop_total_height / 2.0, paraBoxWidth / 2.0);

  G4LogicalVolume* logic_full_loop_para_mother = new G4LogicalVolume(solid_full_loop_para_mother,
                                                              G4Material::GetMaterial("G4_AIR"), "logic_full_loop_para_mother",
                                                              0, 0, 0);
  m_DisplayAction->AddVolume(logic_full_loop_para_mother, "ParaEnvelope");


  G4VSolid* solid_left_half_loop = new G4Torus("solid_left_half_loop",
                                                0, fiber_thickness / 2.0,
                                                min_fiber_bending_radius,
                                                0.5 * M_PIl, 1.0 * M_PIl);
  // G4VSolid* solid_right_half_loop = new G4Torus("solid_right_half_loop",
  //                                               0, fiber_thickness / 2.0,
  //                                               min_fiber_bending_radius,
  //                                               -0.5 * M_PIl, 1.0 * M_PIl);
  G4VSolid* solid_long_fiber_straight = new G4Tubs("solid_long_fiber_straight",
                                                    0,
                                                    fiber_thickness / 2.0,
                                                    (TowerDx - 2.0 * min_fiber_bending_radius - 2.0*edge_distance) / 2.0,
                                                    0, 2.0 * M_PIl);
  G4VSolid* solid_xlong_fiber_straight = new G4Tubs("solid_xlong_fiber_straight",
                                                    0,
                                                    fiber_thickness / 2.0,
                                                    (TowerDx - min_fiber_bending_radius - edge_distance) / 2.0,
                                                    0, 2.0 * M_PIl);

  G4LogicalVolume* logic_left_half_loop = new G4LogicalVolume(solid_left_half_loop,
                                                              GetScintillatorMaterial(), "logic_left_half_loop",
                                                              0, 0, 0);
  // G4LogicalVolume* logic_right_half_loop = new G4LogicalVolume(solid_right_half_loop,
  //                                                             GetScintillatorMaterial(), "logic_right_half_loop",
  //                                                             0, 0, 0);
  G4LogicalVolume* logic_long_fiber_straight = new G4LogicalVolume(solid_long_fiber_straight,
                                                              GetScintillatorMaterial(), "logic_long_fiber_straight",
                                                              0, 0, 0);
  G4LogicalVolume* logic_xlong_fiber_straight = new G4LogicalVolume(solid_xlong_fiber_straight,
                                                              GetScintillatorMaterial(), "logic_xlong_fiber_straight",
                                                              0, 0, 0);

  // m_DisplayAction->AddVolume(logic_right_half_loop, "Scintillator");
  m_DisplayAction->AddVolume(logic_left_half_loop, "Scintillator");
  m_DisplayAction->AddVolume(logic_long_fiber_straight, "Scintillator");
  m_DisplayAction->AddVolume(logic_xlong_fiber_straight, "Scintillator");
  SurfaceTable(logic_left_half_loop, "logic_left_half_loop");
  // SurfaceTable(logic_right_half_loop, "logic_right_half_loop");
  SurfaceTable(logic_long_fiber_straight, "logic_long_fiber_straight");
  SurfaceTable(logic_xlong_fiber_straight, "logic_xlong_fiber_straight");

  G4VSolid* solid_cutout_left_half_loop = new G4Torus("solid_cutout_left_half_loop",
                                                0, 1.07*fiber_thickness / 2.0,
                                                min_fiber_bending_radius,
                                                0.4 * M_PIl, 1.2 * M_PIl);
  // G4VSolid* solid_cutout_right_half_loop = new G4Torus("solid_cutout_right_half_loop",
  //                                               0, 1.07*fiber_thickness / 2.0,
  //                                               min_fiber_bending_radius,
  //                                               -0.6 * M_PIl, 1.2 * M_PIl);
  G4VSolid* solid_cutout_long_fiber_straight = new G4Tubs("solid_cutout_long_fiber_straight",
                                                    0,
                                                    1.07*fiber_thickness / 2.0,
                                                    1.05*(TowerDx - 2.0 * min_fiber_bending_radius - 2*edge_distance) / 2.0,
                                                    0, 2.0 * M_PIl);
  G4VSolid* solid_cutout_xlong_fiber_straight = new G4Tubs("solid_cutout_xlong_fiber_straight",
                                                    0,
                                                    1.07*fiber_thickness / 2.0,
                                                    1.05*(TowerDx - min_fiber_bending_radius - edge_distance) / 2.0,
                                                    0, 2.0 * M_PIl);
  for(int iloop=-(nLoopsFiber-1)/2; iloop<=(nLoopsFiber-1)/2; iloop++){
    if(iloop!=0){
      G4RotationMatrix* rotm_fibr = new G4RotationMatrix();
      rotm_fibr->rotateY(90 * deg);
      G4double yshiftTmp = min_fiber_bending_radius;
      if(iloop>0) yshiftTmp = -min_fiber_bending_radius;
        solid_full_absorber_para_center = new G4SubtractionSolid("solid_full_absorber_para_mod1_" + std::to_string(iloop),
                                                solid_full_absorber_para_center, solid_cutout_long_fiber_straight,
                                                rotm_fibr, G4ThreeVector(0, iloop*(2*min_fiber_bending_radius) + yshiftTmp, 0));
      if(abs(iloop)==(nLoopsFiber-1)/2){
        solid_full_absorber_para_center = new G4SubtractionSolid("solid_full_absorber_para_mod2_" + std::to_string(iloop),
                                                solid_full_absorber_para_center, solid_cutout_xlong_fiber_straight,
                                                rotm_fibr, G4ThreeVector(0, iloop*(2*min_fiber_bending_radius) - yshiftTmp, 0));
        solid_full_absorber_para_right = new G4SubtractionSolid("solid_full_absorber_para_mod2r_" + std::to_string(iloop),
                                                solid_full_absorber_para_right, solid_cutout_xlong_fiber_straight,
                                                rotm_fibr, G4ThreeVector(0, iloop*(2*min_fiber_bending_radius) - yshiftTmp, 0));
      }
    }
    if(iloop%2 == 0){
        solid_full_absorber_para_left = new G4SubtractionSolid("solid_full_absorber_para_mod3_" + std::to_string(iloop),
                                                solid_full_absorber_para_left, solid_cutout_left_half_loop,
                                                0, G4ThreeVector(+min_fiber_bending_radius/2.0+edge_distance/2.0, iloop*(2.0*min_fiber_bending_radius), 0));
    } else {
      G4RotationMatrix* rotm_rcu = new G4RotationMatrix();
      rotm_rcu->rotateY(180 * deg);
        solid_full_absorber_para_right = new G4SubtractionSolid("solid_full_absorber_para_mod4_" + std::to_string(iloop),
                                                solid_full_absorber_para_right, solid_cutout_left_half_loop,
                                                rotm_rcu, G4ThreeVector(-min_fiber_bending_radius/2.0-edge_distance/2.0, iloop*(2.0*min_fiber_bending_radius), 0));
    }
  }

  G4LogicalVolume* logic_full_absorber_para_left = new G4LogicalVolume(solid_full_absorber_para_left,
                                                              G4Material::GetMaterial("G4_Cu"), "logic_full_absorber_para_left",
                                                              0, 0, 0);
  G4LogicalVolume* logic_full_absorber_para_center = new G4LogicalVolume(solid_full_absorber_para_center,
                                                              G4Material::GetMaterial("G4_Cu"), "logic_full_absorber_para_center",
                                                              0, 0, 0);
  G4LogicalVolume* logic_full_absorber_para_right = new G4LogicalVolume(solid_full_absorber_para_right,
                                                              G4Material::GetMaterial("G4_Cu"), "logic_full_absorber_para_right",
                                                              0, 0, 0);
  m_DisplayAction->AddVolume(logic_full_absorber_para_left, "Absorber");
  m_DisplayAction->AddVolume(logic_full_absorber_para_center, "Absorber");
  m_DisplayAction->AddVolume(logic_full_absorber_para_right, "Absorber");
  new G4PVPlacement(0, G4ThreeVector(-(TowerDx - min_fiber_bending_radius - edge_distance) / 2.0, 0, 0),
                    logic_full_absorber_para_left,
                    "physvol_full_absorber_para_left",
                    logic_full_loop_para,
                    0, 0, OverlapCheck());
  new G4PVPlacement(0, G4ThreeVector(0, 0, 0),
                    logic_full_absorber_para_center,
                    "physvol_full_absorber_para_center",
                    logic_full_loop_para,
                    0, 0, OverlapCheck());
  new G4PVPlacement(0, G4ThreeVector((TowerDx - min_fiber_bending_radius - edge_distance) / 2.0, 0, 0),
                    logic_full_absorber_para_right,
                    "physvol_full_absorber_para_right",
                    logic_full_loop_para,
                    0, 0, OverlapCheck());

  G4RotationMatrix* rotm_full_loop = new G4RotationMatrix();
  rotm_full_loop->rotateZ(tilting_angle);
  rotm_full_loop->rotateY(0.5*M_PIl);
  //G4VPhysicalVolume* physvol_full_loop_para =

  new G4PVPlacement(rotm_full_loop, G4ThreeVector(0, 0, 0),
                    logic_full_loop_para,
                    "physvol_full_loop_para",
                    logic_full_loop_para_mother,
                    0, 0, OverlapCheck());
  G4VPhysicalVolume* physvol_long_fiber_straight[nLoopsFiber+2] = {nullptr};
  G4VPhysicalVolume* physvol_xlong_fiber_straight[nLoopsFiber+2] = {nullptr};
  G4VPhysicalVolume* physvol_left_half_loop[nLoopsFiber+2] = {nullptr};
  G4VPhysicalVolume* physvol_right_half_loop[nLoopsFiber+2] = {nullptr};
  int indxarr = 0;
  for(int iloop=-(nLoopsFiber-1)/2-1; iloop<=(nLoopsFiber-1)/2+1; iloop++){
    if(iloop!=0){
      G4RotationMatrix* rotm_fibr = new G4RotationMatrix();
      rotm_fibr->rotateY(90 * deg);
      G4double yshiftTmp = min_fiber_bending_radius;
      if(iloop>0) yshiftTmp = -min_fiber_bending_radius;
      if(abs(iloop)!=(nLoopsFiber-1)/2+1){
        physvol_long_fiber_straight[indxarr] =
        new G4PVPlacement(rotm_fibr, G4ThreeVector(0, iloop*(2.0*min_fiber_bending_radius) + yshiftTmp, 0),
                                                                    logic_long_fiber_straight,
                                                                    "physvol_long_fiber_straight_" + std::to_string(iloop),
                                                                    logic_full_loop_para,
                                                                    0, 0, OverlapCheck());
      }
      if(abs(iloop)==(nLoopsFiber-1)/2){
        physvol_xlong_fiber_straight[indxarr] =
        new G4PVPlacement(rotm_fibr, G4ThreeVector((min_fiber_bending_radius + edge_distance) / 2.0, iloop*(2*min_fiber_bending_radius) - yshiftTmp, 0),
        // new G4PVPlacement(rotm_fibr, G4ThreeVector((min_fiber_bending_radius + edge_distance) / 2.0, iloop*(2*min_fiber_bending_radius) - yshiftTmp, 0),
                                                                    logic_xlong_fiber_straight,
                                                                    "physvol_xlong_fiber_straight_" + std::to_string(iloop),
                                                                    logic_full_loop_para,
                                                                    0, 0, OverlapCheck());
      }
    }
    indxarr++;
  }
  indxarr = 0;
  for(int iloop=-(nLoopsFiber-1)/2-1; iloop<=(nLoopsFiber-1)/2+1; iloop++){
    if(abs(iloop)<(nLoopsFiber-1)/2+1){
      if(iloop%2 == 0){
        physvol_left_half_loop[indxarr] =
        new G4PVPlacement(0, G4ThreeVector(-(TowerDx - min_fiber_bending_radius - edge_distance) / 2.0+min_fiber_bending_radius/2.0 + edge_distance/2.0, iloop*(2.0*min_fiber_bending_radius), 0),
                                                                    logic_left_half_loop,
                                                                    "physvol_left_half_loop_" + std::to_string(iloop),
                                                                    logic_full_loop_para,
                                                                    0, 0, OverlapCheck());
        if(physvol_long_fiber_straight[indxarr+1]){
          MakeBoundaryFibers(physvol_long_fiber_straight[indxarr+1],physvol_left_half_loop[indxarr],"long_left_+1_boundary" + std::to_string(iloop));
          MakeBoundaryFibers(physvol_left_half_loop[indxarr],physvol_long_fiber_straight[indxarr+1],"left_long_+1_boundary" + std::to_string(iloop));
        }
        if(physvol_long_fiber_straight[indxarr]){
          MakeBoundaryFibers(physvol_long_fiber_straight[indxarr],physvol_left_half_loop[indxarr],"long_left_0_boundary" + std::to_string(iloop));
          MakeBoundaryFibers(physvol_left_half_loop[indxarr],physvol_long_fiber_straight[indxarr],"left_long_0_boundary" + std::to_string(iloop));
        }
        if(indxarr>0){
          if(physvol_long_fiber_straight[indxarr-1]){
            MakeBoundaryFibers(physvol_long_fiber_straight[indxarr-1],physvol_left_half_loop[indxarr],"long_left_-1_boundary" + std::to_string(iloop));
            MakeBoundaryFibers(physvol_left_half_loop[indxarr],physvol_long_fiber_straight[indxarr-1],"left_long_-1_boundary" + std::to_string(iloop));
          }
        }
        if(physvol_xlong_fiber_straight[indxarr+1]){
          MakeBoundaryFibers(physvol_xlong_fiber_straight[indxarr+1],physvol_left_half_loop[indxarr],"xlong_left_+1_boundary" + std::to_string(iloop));
          MakeBoundaryFibers(physvol_left_half_loop[indxarr],physvol_xlong_fiber_straight[indxarr+1],"left_xlong_+1_boundary" + std::to_string(iloop));
        }
        if(physvol_xlong_fiber_straight[indxarr]){
          MakeBoundaryFibers(physvol_xlong_fiber_straight[indxarr],physvol_left_half_loop[indxarr],"xlong_left_0_boundary" + std::to_string(iloop));
          MakeBoundaryFibers(physvol_left_half_loop[indxarr],physvol_xlong_fiber_straight[indxarr],"left_xlong_0_boundary" + std::to_string(iloop));
        }
        if(indxarr>0){
          if(physvol_xlong_fiber_straight[indxarr-1]){
            MakeBoundaryFibers(physvol_xlong_fiber_straight[indxarr-1],physvol_left_half_loop[indxarr],"xlong_left_-1_boundary" + std::to_string(iloop));
            MakeBoundaryFibers(physvol_left_half_loop[indxarr],physvol_xlong_fiber_straight[indxarr-1],"left_xlong_-1_boundary" + std::to_string(iloop));
          }
        }
      } else {
      G4RotationMatrix* rotm_loop2 = new G4RotationMatrix();
      rotm_loop2->rotateY(180 * deg);
      rotm_loop2->rotateX(180 * deg);
        physvol_right_half_loop[indxarr] =
        new G4PVPlacement(rotm_loop2, G4ThreeVector((TowerDx - min_fiber_bending_radius - edge_distance) / 2.0-min_fiber_bending_radius/2.0 - edge_distance/2.0, iloop*(2.0*min_fiber_bending_radius), 0),
                                                                    logic_left_half_loop,
                                                                    "physvol_right_half_loop_" + std::to_string(iloop),
                                                                    logic_full_loop_para,
                                                                    0, 0, OverlapCheck());
        // new G4PVPlacement(0, G4ThreeVector((TowerDx - min_fiber_bending_radius - edge_distance) / 2.0-min_fiber_bending_radius/2.0 - edge_distance/2.0, iloop*(2.0*min_fiber_bending_radius), 0),
        //                                                             logic_right_half_loop,
        //                                                             "physvol_right_half_loop_" + std::to_string(iloop),
        //                                                             logic_full_loop_para,
        //                                                             0, 0, OverlapCheck());
        if(physvol_long_fiber_straight[indxarr+1]){
          MakeBoundaryFibers(physvol_long_fiber_straight[indxarr+1],physvol_right_half_loop[indxarr],"long_right_+1_boundary" + std::to_string(iloop));
          MakeBoundaryFibers(physvol_right_half_loop[indxarr],physvol_long_fiber_straight[indxarr+1],"right_long_+1_boundary" + std::to_string(iloop));
        }
        if(physvol_long_fiber_straight[indxarr]){
          MakeBoundaryFibers(physvol_long_fiber_straight[indxarr],physvol_right_half_loop[indxarr],"long_right_0_boundary" + std::to_string(iloop));
          MakeBoundaryFibers(physvol_right_half_loop[indxarr],physvol_long_fiber_straight[indxarr],"right_long_0_boundary" + std::to_string(iloop));
        }
        if(indxarr>0){
          if(physvol_long_fiber_straight[indxarr-1]){
            MakeBoundaryFibers(physvol_long_fiber_straight[indxarr-1],physvol_right_half_loop[indxarr],"long_right_-1_boundary" + std::to_string(iloop));
            MakeBoundaryFibers(physvol_right_half_loop[indxarr],physvol_long_fiber_straight[indxarr-1],"right_long_-1_boundary" + std::to_string(iloop));
          }
        }
      }
    }
    indxarr++;
  }
  G4double length_replsolid = TowerDz - sqrt(pow(loop_total_height,2)-pow(TowerDx,2))-1*cm;
  G4VSolid* solid_TF_Replica = new G4Para("DRCalRowBox",
                                                length_replsolid / 2.0, TowerDx / 2.0, TowerDx / 2.0,
                                                tilting_angle, 0, 0);
  auto logic_TF_Replica  = new G4LogicalVolume(solid_TF_Replica,G4Material::GetMaterial("G4_AIR"),"logic_TF_Replica");
  int nLayers = (int) (length_replsolid) / (fiber_thickness + fiber_spacing);
  for(int ilay = 0; ilay<nLayers;ilay++){
    int copynumber = ilay;
    new G4PVPlacement(0, G4ThreeVector(-length_replsolid/2+(fiber_thickness + fiber_spacing)/2+ilay*(fiber_thickness + fiber_spacing), 0, 0),
                            logic_full_loop_para_mother, "placed_mother_layer_" + std::to_string(ilay), logic_TF_Replica, 0, copynumber, OverlapCheck());

  }

  // G4VSolid* solid_TF_Replica = new G4Para("DRCalRowBox",
  //                                               (TowerDz) / 2.0, TowerDx / 2.0, TowerDx / 2.0,
  //                                               tilting_angle, 0, 0);
  // // auto solid_TF_Replica    = new G4Box("DRCalRowBox",1.01*TowerDz / 2.0, 1.01*TowerDx / 2.0,1.01*TowerDx / 2.0);
  // auto logic_TF_Replica  = new G4LogicalVolume(solid_TF_Replica,G4Material::GetMaterial("G4_AIR"),"logic_TF_Replica");
  m_DisplayAction->AddVolume(logic_TF_Replica, "ParaEnvelope");
  // // replicate singletower tower design currRowNtow times along x-axis
  // // int nReplicas = (int) (TowerDz - sqrt(pow(loop_total_height,2)-pow(TowerDx,2))) / paraBoxWidth;
  // int nReplicas = (int) (TowerDz) / (fiber_thickness + fiber_spacing);
  // // int nReplicas = (int) (TowerDz - sqrt(pow(loop_total_height,2)-pow(TowerDx,2))) / (fiber_thickness + fiber_spacing);
  // new G4PVReplica("DRCalRowPhysical",logic_full_loop_para_mother,logic_TF_Replica,
  //                 kXAxis,nReplicas,(fiber_thickness + fiber_spacing),0);
  G4RotationMatrix* rotm_loop_stack = new G4RotationMatrix();
  rotm_loop_stack->rotateY(0.5*M_PIl);
  new G4PVPlacement(rotm_loop_stack, G4ThreeVector(0, 0, 0),
                          logic_TF_Replica, "blubb", logic_envelope, 0, false, OverlapCheck());
  return true;
}
//_______________________________________________________________________
bool PHG4FoCalDetector::ConstructCapillaryRowDetector(int type,G4LogicalVolume* envelope)
{
  if (Verbosity() > 0)
  {
    cout << "PHG4FoCalDetector: Build assembly volume for capillary row..." << endl;
  }


  //create geometry volumes to place inside single_tower

  G4double copperTubeOuterDiam = 2.5*mm;
  G4double copperTubeInnerDiam = 1.2*mm;
  G4double diameter_fiber = 1.0*mm;

  // fibre cutout
  G4VSolid* solid_filledcap  = new G4Tubs(G4String("ttl_solid_filledcap"),
                                            0,
                                            copperTubeOuterDiam / 2.0,
                                            _tower_dz / 2.0,
                                            0.,2*M_PIl*rad);

  G4VSolid* solid_absorber  = new G4Tubs(G4String("ttl_copper_tube_solid"),
                                            copperTubeInnerDiam / 2.0,
                                            copperTubeOuterDiam / 2.0,
                                            _tower_dz / 2.0,
                                            0.,2*M_PIl*rad);

  G4VSolid* solid_scintillator  = new G4Tubs(G4String("single_scintillator_fiber"),
                                            0,
                                            diameter_fiber / 2.0,
                                            _tower_dz / 2.0,
                                            0.,2*M_PIl*rad);

  // G4VSolid* solid_airgap  = new G4Tubs(G4String("single_airgap"),
  //                                           diameter_fiber / 2.0,
  //                                           copperTubeInnerDiam / 2.0,
  //                                           _tower_dz / 2.0,
  //                                           0.,2*M_PIl*rad);


  G4Material* material_scintillator = GetScintillatorMaterial();


  G4NistManager* man = G4NistManager::Instance();
  G4Material* material_absorber;
  if(_absorber_Material==0)material_absorber = man->FindOrBuildMaterial(_materialAbsorber.c_str());
  else if(_absorber_Material==1)material_absorber = man->FindOrBuildMaterial("G4_W");
  else if(_absorber_Material==2)material_absorber = man->FindOrBuildMaterial("G4_Cu");
  else if(_absorber_Material==3)material_absorber = man->FindOrBuildMaterial("G4_Pb");
  else if(_absorber_Material==4)material_absorber = man->FindOrBuildMaterial("G4_Fe");
  else material_absorber = man->FindOrBuildMaterial(_materialAbsorber.c_str());


  G4LogicalVolume* logic_filledcap = new G4LogicalVolume(solid_filledcap,
                                                        G4Material::GetMaterial("G4_AIR"),
                                                        "logic_filledcap",
                                                        0, 0, 0);

  // G4LogicalVolume* logic_airgap = new G4LogicalVolume(solid_airgap,
  //                                                       G4Material::GetMaterial("G4_AIR"),
  //                                                       "logic_airgap",
  //                                                       0, 0, 0);

  G4LogicalVolume* logic_absorber = new G4LogicalVolume(solid_absorber,
                                                        material_absorber,
                                                        "absorber_solid_logic",
                                                        0, 0, 0);

  G4LogicalVolume* logic_scint = new G4LogicalVolume(solid_scintillator,
                                                    material_scintillator,
                                                    "hfocal_single_scintillator_fiber_logic",
                                                    0, 0, 0);
  if(m_doLightProp){
    SurfaceTable(logic_scint, "hfocal_single_scintillator_fiber_logic");
    // SurfaceTable(logic_absorber);
  }
  m_DisplayAction->AddVolume(logic_filledcap, "FfocalEnvelope");
  m_DisplayAction->AddVolume(logic_absorber, "Absorber");
  m_DisplayAction->AddVolume(logic_scint, "Scintillator");
  // m_DisplayAction->AddVolume(logic_airgap, "Invisible");


  new G4PVPlacement(0, G4ThreeVector(0, 0, 0),
                  logic_absorber,
                  _towerlogicnameprefix + "absorbersolid_copper",
                  logic_filledcap,
                  0, 0, OverlapCheck());
  //G4VPhysicalVolume* placed_logic_scint = 
  new G4PVPlacement(0, G4ThreeVector(0, 0, 0),
                  logic_scint,
                  _towerlogicnameprefix + "singlescintillatorfiber_scint",
                  logic_filledcap,
                  0, 0, OverlapCheck());
  //G4VPhysicalVolume* placed_logic_airgap = 
  // new G4PVPlacement(0, G4ThreeVector(0, 0, 0),
  //                 logic_airgap,
  //                 _towerlogicnameprefix + "single_airgap",
  //                 logic_filledcap,
  //                 0, 0, OverlapCheck());
  // if(m_doLightProp){
  //   MakeBoundary(placed_logic_scint,placed_logic_airgap);
  // }
  //place physical volumes for absorber and scintillator fiber

  /* Loop over all tower positions in vector and place tower */
  for (std::map<std::string, towerposition>::iterator iterator = m_TowerPostionMap.begin(); iterator != m_TowerPostionMap.end(); ++iterator)
  {
    if (Verbosity() > 0)
    {
      std::cout << "PHG4ForwardHcalDetector: Place tower " << iterator->first
                << " idx_j = " << iterator->second.idx_j << ", idx_k = " << iterator->second.idx_k
                << " at x = " << iterator->second.x << " , y = " << iterator->second.y << " , z = " << iterator->second.z << std::endl;
    }

    int copyno = (iterator->second.idx_j << 16) + iterator->second.idx_k;
    new G4PVPlacement(0, G4ThreeVector(iterator->second.x, iterator->second.y, iterator->second.z),
                      logic_filledcap,
                      iterator->first,
                      envelope,
                      0, copyno, OverlapCheck());
  }


  // G4double xposition,yposition;
  // int icapil = 0;
  // for(int icap_x=-18;icap_x<18;icap_x++){
  //   xposition = icap_x * copperTubeOuterDiam; //offset for odd/even rows
  //   for(int icap_y=-20;icap_y<20;icap_y++){
  //     yposition = icap_y * copperTubeOuterDiam * sqrt(3)/2;
  //     icap_y%2==0 ? xposition+=copperTubeOuterDiam/2 : xposition-=copperTubeOuterDiam/2;
  //     // new G4PVPlacement(0, G4ThreeVector(xposition, yposition, 0),
  //     //                 logic_absorber,
  //     //                 _towerlogicnameprefix + "absorbersolid_" + std::to_string(icap_x) + "_"  + std::to_string(icap_y),
  //     //                 envelope,
  //     //                 0, icapil, OverlapCheck());
  //     new G4PVPlacement(0, G4ThreeVector(xposition, yposition, 0),
  //                     logic_filledcap,
  //                     _towerlogicnameprefix + "tower_j_" + std::to_string(icap_x+18) + "_k_"  + std::to_string(icap_y+20),
  //                     envelope,
  //                     0, icapil, OverlapCheck());
  //     // new G4PVPlacement(0, G4ThreeVector(xposition, yposition, 0),
  //     //                 logic_scint,
  //     //                 _towerlogicnameprefix + "singlescintillatorfiber_j_" + std::to_string(icap_x) + "_k_"  + std::to_string(icap_y),
  //     //                 envelope,
  //     //                 0, icapil, OverlapCheck());
  //     icapil++;
  //   }
  // }
  bool addEnclosure = false;
  if(addEnclosure){
    G4double enclosure_width = 1.5*mm;
    G4double enclosure_height = 40*copperTubeOuterDiam* sqrt(3)/2;
    G4VSolid* enclosure_solid_sides = new G4Box(G4String("enclosure_solid_sides"),
                                            enclosure_width / 2.0,
                                            enclosure_height / 2.0,
                                            _tower_dz / 2.0);

    G4LogicalVolume* enclosure_logic_sides = new G4LogicalVolume(enclosure_solid_sides,
                                                              material_absorber,
                                                              "enclosure_logic_sides",
                                                              0, 0, 0);

    m_DisplayAction->AddVolume(enclosure_logic_sides, "Absorber");
    float xposition_enc = -18 * copperTubeOuterDiam - enclosure_width/2 - copperTubeOuterDiam/2;
    new G4PVPlacement(0, G4ThreeVector(xposition_enc, - copperTubeOuterDiam/2, 0),
                    enclosure_logic_sides,
                    "enclosuresolid_leftside",
                    envelope,
                    0, 0, OverlapCheck());
    xposition_enc = 17 * copperTubeOuterDiam + enclosure_width/2 + copperTubeOuterDiam;
    new G4PVPlacement(0, G4ThreeVector(xposition_enc, - copperTubeOuterDiam/2, 0),
                    enclosure_logic_sides,
                    "enclosuresolid_rightside",
                    envelope,
                    0, 0, OverlapCheck());

    enclosure_height = 1.5*mm;
    enclosure_width = 37*copperTubeOuterDiam+enclosure_height;
    G4VSolid* enclosure_solid_sides2 = new G4Box(G4String("enclosure_solid_sides2"),
                                            enclosure_width / 2.0,
                                            enclosure_height / 2.0,
                                            _tower_dz / 2.0);

    G4LogicalVolume* enclosure_logic_sides2 = new G4LogicalVolume(enclosure_solid_sides2,
                                                              material_absorber,
                                                              "enclosure_logic_sides2",
                                                              0, 0, 0);

    m_DisplayAction->AddVolume(enclosure_logic_sides2, "Absorber");
    float yposition_enc = -20 * (copperTubeOuterDiam * sqrt(3)/2) - enclosure_height/2-copperTubeOuterDiam/2;
    new G4PVPlacement(0, G4ThreeVector(-enclosure_height/2, yposition_enc, 0),
                    enclosure_logic_sides2,
                    "enclosuresolid_top",
                    envelope,
                    0, 0, OverlapCheck());
    yposition_enc = 19 * (copperTubeOuterDiam * sqrt(3)/2) + enclosure_height/2+copperTubeOuterDiam/2;
    new G4PVPlacement(0, G4ThreeVector(-enclosure_height/2, yposition_enc, 0),
                    enclosure_logic_sides2,
                    "enclosuresolid_bottom",
                    envelope,
                    0, 0, OverlapCheck());
  

    // yposition = icap_y * copperTubeOuterDiam * sqrt(3)/2;
  }


  if (Verbosity() > 0)
  {
    cout << "PHG4FoCalDetector: Building logical volume for single tower done." << endl;
  }
  return true;
}

//_______________________________________________________________________
G4LogicalVolume*
PHG4FoCalDetector::ConstructTower(int type)
{
  if (Verbosity() > 0)
  {
    cout << "PHG4FoCalDetector: Build logical volume for single tower..." << endl;
  }

  //create logical volume for single tower
  G4Material* material_air = G4Material::GetMaterial("G4_AIR");
  // 2x2 tower base element
  G4VSolid* single_tower_solid = new G4Box(G4String("single_tower_solid"),
                                          _tower_dx / 2.0,
                                          _tower_dy / 2.0,
                                          _tower_dz / 2.0);

  G4LogicalVolume* single_tower_logic = new G4LogicalVolume(single_tower_solid,
                                                            material_air,
                                                            "single_tower_logic",
                                                            0, 0, 0);

  //create geometry volumes to place inside single_tower

  G4double copperTubeDiam = _tower_dx/2;
  G4double diameter_fiber = _scintFiber_diam;
  G4double diameter_fiber_cherenkov = _cerenkovFiber_diam;

  // fibre cutout
  G4VSolid* solid_absorber  = new G4Tubs(G4String("ttl_copper_tube_solid"),
                                            0,
                                            copperTubeDiam / 2.0,
                                            _tower_dz / 2.0,
                                            0.,2*M_PIl*rad);

  G4VSolid* solid_scintillator  = new G4Tubs(G4String("single_scintillator_fiber"),
                                            0,
                                            diameter_fiber / 2.0,
                                            _tower_dz / 2.0,
                                            0.,2*M_PIl*rad);


  G4VSolid* solid_cherenkov  = new G4Tubs(G4String("single_cherenkov_fiber"),
                                            0,
                                            diameter_fiber_cherenkov / 2.0,
                                            _tower_dz / 2.0,
                                            0.,2*M_PIl*rad);

  G4VSolid* solid_fiberMargin  = new G4Tubs(G4String("solid_fiberMargin"),
                                            0,
                                            (diameter_fiber_cherenkov+0.05*mm) / 2.0,
                                            (_tower_dz*1.1) / 2.0,
                                            0.,2*M_PIl*rad);
  G4VSolid* solid_absorberMargin  = new G4Tubs(G4String("solid_absorberMargin"),
                                            0,
                                            (copperTubeDiam+0.05*mm) / 2.0,
                                            (_tower_dz*1.1) / 2.0,
                                            0.,2*M_PIl*rad);


  G4Material* material_scintillator = GetScintillatorMaterial();


  G4VSolid* solid_fill_hollowspace = new G4Box(G4String("solid_fill_hollowspace"),
                                          copperTubeDiam / 2.0,
                                          copperTubeDiam / 4.0,
                                          _tower_dz / 2.0);
  solid_fill_hollowspace = new G4SubtractionSolid(G4String("solid_fill_hollowspace_2"), solid_fill_hollowspace, solid_absorberMargin,
                                                  0 ,G4ThreeVector( copperTubeDiam / 2.0, -copperTubeDiam / 4.0 ,0.)); // top right
  solid_fill_hollowspace = new G4SubtractionSolid(G4String("solid_fill_hollowspace_3"), solid_fill_hollowspace, solid_absorberMargin,
                                                  0 ,G4ThreeVector( -copperTubeDiam / 2.0 +0.01*mm, -copperTubeDiam / 4.0 ,0.)); // top right
  solid_fill_hollowspace = new G4SubtractionSolid(G4String("solid_fill_hollowspace_1"), solid_fill_hollowspace, solid_fiberMargin,
                                                  0 ,G4ThreeVector( -0.01*mm , copperTubeDiam / 4.0 ,0.0)); // cut out fiber


  G4NistManager* man = G4NistManager::Instance();
  G4Material* material_absorber;
  if(_absorber_Material==0)material_absorber = man->FindOrBuildMaterial(_materialAbsorber.c_str());
  else if(_absorber_Material==1)material_absorber = man->FindOrBuildMaterial("G4_W");
  else if(_absorber_Material==2)material_absorber = man->FindOrBuildMaterial("G4_Cu");
  else if(_absorber_Material==3)material_absorber = man->FindOrBuildMaterial("G4_Pb");
  else if(_absorber_Material==4)material_absorber = man->FindOrBuildMaterial("G4_Fe");
  else material_absorber = man->FindOrBuildMaterial(_materialAbsorber.c_str());


  G4LogicalVolume* logic_absorber = new G4LogicalVolume(solid_absorber,
                                                        material_absorber,
                                                        "absorber_solid_logic",
                                                        0, 0, 0);

  G4LogicalVolume* logic_scint = new G4LogicalVolume(solid_scintillator,
                                                    material_scintillator,
                                                    "hfocal_single_scintillator_fiber_logic",
                                                    0, 0, 0);

  G4Material *material_cherenkov;
  // G4Material *material_cherenkov;
  if(_cerenkovFiber_material==0) material_cherenkov = GetPMMAMaterial();
  else if(_cerenkovFiber_material==1) material_cherenkov = GetQuartzMaterial();
  else material_cherenkov = GetPMMAMaterial();

  G4LogicalVolume* logic_cherenk = new G4LogicalVolume(solid_cherenkov,
                                                    material_cherenkov,
                                                    "hfocal_single_cherenkov_fiber_logic",
                                                    0, 0, 0);

  G4LogicalVolume* logic_fill_hollowspace = new G4LogicalVolume(solid_fill_hollowspace,
                                                          material_absorber,
                                                          "logic_fill_hollowspace",
                                                          0, 0, 0);
  m_DisplayAction->AddVolume(logic_absorber, "Absorber");
  m_DisplayAction->AddVolume(logic_fill_hollowspace, "Fill");
  m_DisplayAction->AddVolume(logic_scint, "Scintillator");
  m_DisplayAction->AddVolume(logic_cherenk, "Cherenkov");
  m_DisplayAction->AddVolume(single_tower_logic, "FfocalEnvelope");

  //place physical volumes for absorber and scintillator fiber

  ostringstream name_absorber;
  name_absorber.str("");
  name_absorber << _towerlogicnameprefix << "absorbersolid" << endl;

  ostringstream name_scintillator;
  name_scintillator.str("");
  name_scintillator << _towerlogicnameprefix << "singlescintillatorfiber" << endl;

  ostringstream name_cherenkov;
  name_cherenkov.str("");
  name_cherenkov << _towerlogicnameprefix << "_single_cherenkov_fiber"  << endl;

  ostringstream name_fill;
  name_fill.str("");
  name_fill << _towerlogicnameprefix << "fillsolid" << endl;

  new G4PVPlacement(0, G4ThreeVector( -_tower_dx/4,  _tower_dx/4 , 0),
                    logic_absorber,
                    name_absorber.str().c_str(),
                    single_tower_logic,
                    0, 0, OverlapCheck());
  new G4PVPlacement(0, G4ThreeVector( -_tower_dx/4,  -_tower_dx/4 , 0),
                    logic_absorber,
                    name_absorber.str().c_str(),
                    single_tower_logic,
                    0, 0, OverlapCheck());
  new G4PVPlacement(0, G4ThreeVector( _tower_dx/4,  _tower_dx/4 , 0),
                    logic_absorber,
                    name_absorber.str().c_str(),
                    single_tower_logic,
                    0, 0, OverlapCheck());
  new G4PVPlacement(0, G4ThreeVector( _tower_dx/4,  -_tower_dx/4 , 0),
                    logic_absorber,
                    name_absorber.str().c_str(),
                    single_tower_logic,
                    0, 0, OverlapCheck());
  if(_fillhollow){
    // place hollow space filling
    new G4PVPlacement(0, G4ThreeVector( 0,  -copperTubeDiam / 4.0 , 0),
                      logic_fill_hollowspace,
                      name_fill.str().c_str(),
                      single_tower_logic,
                      0, 0, OverlapCheck());
    new G4PVPlacement(0, G4ThreeVector( 0, _tower_dx/2 -copperTubeDiam / 4.0 , 0),
                      logic_fill_hollowspace,
                      name_fill.str().c_str(),
                      single_tower_logic,
                      0, 0, OverlapCheck());
    G4RotationMatrix *filling_rotup = new G4RotationMatrix();
    filling_rotup->rotateZ(M_PIl);
    new G4PVPlacement(filling_rotup, G4ThreeVector( 0,  copperTubeDiam / 4.0 , 0),
                      logic_fill_hollowspace,
                      name_fill.str().c_str(),
                      single_tower_logic,
                      0, 0, OverlapCheck());
    new G4PVPlacement(filling_rotup, G4ThreeVector( 0,  -_tower_dx/2 + copperTubeDiam / 4.0 , 0),
                      logic_fill_hollowspace,
                      name_fill.str().c_str(),
                      single_tower_logic,
                      0, 0, OverlapCheck());
    G4RotationMatrix *filling_rotleft = new G4RotationMatrix();
    filling_rotleft->rotateZ(M_PIl/2);
    new G4PVPlacement(filling_rotleft, G4ThreeVector( _tower_dx/2 - copperTubeDiam / 4.0,  0 , 0),
                      logic_fill_hollowspace,
                      name_fill.str().c_str(),
                      single_tower_logic,
                      0, 0, OverlapCheck());
    new G4PVPlacement(filling_rotleft, G4ThreeVector( _tower_dx/2 - copperTubeDiam / 4.0,  _tower_dx/2  , 0),
                      logic_fill_hollowspace,
                      name_fill.str().c_str(),
                      single_tower_logic,
                      0, 0, OverlapCheck());
    G4RotationMatrix *filling_rotright = new G4RotationMatrix();
    filling_rotright->rotateZ(-M_PIl/2);
    new G4PVPlacement(filling_rotright, G4ThreeVector( -_tower_dx/2 + copperTubeDiam / 4.0,  0 , 0),
                      logic_fill_hollowspace,
                      name_fill.str().c_str(),
                      single_tower_logic,
                      0, 0, OverlapCheck());
    new G4PVPlacement(filling_rotright, G4ThreeVector( -_tower_dx/2 + copperTubeDiam / 4.0,  -_tower_dx/2  , 0),
                      logic_fill_hollowspace,
                      name_fill.str().c_str(),
                      single_tower_logic,
                      0, 0, OverlapCheck());
  }
    // place scintillator fibers (top left, bottom right)
    new G4PVPlacement(0, G4ThreeVector( 0,  0 , 0),
                      logic_scint,
                      name_scintillator.str().c_str(),
                      single_tower_logic,
                      0, 0, OverlapCheck());
    new G4PVPlacement(0, G4ThreeVector( _tower_dx / 2, _tower_dx / 2 , 0),
                      logic_scint,
                      name_scintillator.str().c_str(),
                      single_tower_logic,
                      0, 0, OverlapCheck());

    // place cherenkov fibers (top right, bottom left)
    if(_tower_type==5){
      new G4PVPlacement(0, G4ThreeVector( _tower_dx / 2, 0 , 0),
                        logic_scint,
                        name_scintillator.str().c_str(),
                        single_tower_logic,
                        0, 0, OverlapCheck());
      new G4PVPlacement(0, G4ThreeVector( 0, _tower_dx / 2 , 0),
                        logic_scint,
                        name_scintillator.str().c_str(),
                        single_tower_logic,
                        0, 0, OverlapCheck());
    } else {
      new G4PVPlacement(0, G4ThreeVector( _tower_dx / 2, 0 , 0),
                        logic_cherenk,
                        name_cherenkov.str().c_str(),
                        single_tower_logic,
                        0, 0, OverlapCheck());
      new G4PVPlacement(0, G4ThreeVector( 0, _tower_dx / 2 , 0),
                        logic_cherenk,
                        name_cherenkov.str().c_str(),
                        single_tower_logic,
                        0, 0, OverlapCheck());
    }
  if (Verbosity() > 0)
  {
    cout << "PHG4FoCalDetector: Building logical volume for single tower done." << endl;
  }

  return single_tower_logic;
}


//_______________________________________________________________________
G4LogicalVolume*
PHG4FoCalDetector::ConstructTowerLayered(int type)
{
  if (Verbosity() > 0)
  {
    cout << "PHG4FoCalDetector: Build logical volume for single tower..." << endl;
  }

  //create logical volume for single tower
  G4Material* material_air = G4Material::GetMaterial("G4_AIR");
  // 2x2 tower base element
  G4VSolid* single_tower_solid = new G4Box(G4String("single_tower_solid"),
                                          _tower_dx / 2.0,
                                          _tower_dy / 2.0,
                                          _tower_dz / 2.0);

  G4LogicalVolume* single_tower_logic = new G4LogicalVolume(single_tower_solid,
                                                            material_air,
                                                            "single_tower_logic",
                                                            0, 0, 0);

  //create geometry volumes to place inside single_tower

  G4double thickness_absorber = 3*cm;
  G4double thickness_scintillator = 0.2*cm;

  // fibre cutout
  G4VSolid* solid_absorber  = new G4Box(G4String("single_tower_solid"),
                                          _tower_dx / 2.0,
                                          _tower_dy / 2.0,
                                          thickness_absorber / 2.0);

  G4VSolid* solid_scintillator  =  new G4Box(G4String("single_tower_solid"),
                                          _tower_dx / 2.0,
                                          _tower_dy / 2.0,
                                          thickness_scintillator / 2.0);


  G4Material* material_scintillator = GetScintillatorMaterial();


  G4NistManager* man = G4NistManager::Instance();
  G4Material* material_absorber;
  if(_absorber_Material==0)material_absorber = man->FindOrBuildMaterial(_materialAbsorber.c_str());
  else if(_absorber_Material==1)material_absorber = man->FindOrBuildMaterial("G4_W");
  else if(_absorber_Material==2)material_absorber = man->FindOrBuildMaterial("G4_Cu");
  else if(_absorber_Material==3)material_absorber = man->FindOrBuildMaterial("G4_Pb");
  else if(_absorber_Material==4)material_absorber = man->FindOrBuildMaterial("G4_Fe");
  else material_absorber = man->FindOrBuildMaterial(_materialAbsorber.c_str());


  G4LogicalVolume* logic_absorber = new G4LogicalVolume(solid_absorber,
                                                        material_absorber,
                                                        "absorber_solid_logic",
                                                        0, 0, 0);

  G4LogicalVolume* logic_scint = new G4LogicalVolume(solid_scintillator,
                                                    material_scintillator,
                                                    "hfocal_single_scintillator_fiber_logic",
                                                    0, 0, 0);

  m_DisplayAction->AddVolume(logic_absorber, "Absorber");
  m_DisplayAction->AddVolume(logic_scint, "Scintillator");
  m_DisplayAction->AddVolume(single_tower_logic, "FfocalEnvelope");

  /* place physical volumes for absorber and scintillator plates */
  G4int nlayers = _tower_dz / (thickness_absorber + thickness_scintillator);;
  G4double zpos_i = (-1 * _tower_dz / 2.0) + thickness_absorber / 2.0;

  std::string name_absorber = _towerlogicnameprefix + "_single_plate_absorber";

  std::string name_scintillator = _towerlogicnameprefix + "_single_plate_scintillator";

  for (int i = 1; i <= nlayers; i++)
  {
    new G4PVPlacement(0, G4ThreeVector(0, 0, zpos_i),
                      logic_absorber,
                      name_absorber,
                      single_tower_logic,
                      0, 0, OverlapCheck());

    zpos_i += (thickness_absorber / 2. + thickness_scintillator / 2.);

    new G4PVPlacement(0, G4ThreeVector(0, 0, zpos_i),
                      logic_scint,
                      name_scintillator,
                      single_tower_logic,
                      0, 0, OverlapCheck());

    zpos_i += (thickness_absorber / 2. + thickness_scintillator / 2.);
  }

  if (Verbosity() > 0)
  {
    cout << "PHG4FoCalDetector: Building logical volume for single tower done." << endl;
  }

  return single_tower_logic;
}

//_______________________________________________________________________
G4Material*
PHG4FoCalDetector::GetScintillatorMaterial()
{
  G4Material *material_G4_POLYSTYRENE = G4Material::GetMaterial("G4_POLYSTYRENE_FOCAL");

  if (!material_G4_POLYSTYRENE){
    if (Verbosity() > 0)
    {
      cout << "PHG4FoCalDetector: Making Scintillator material..." << endl;
    }
    G4double density;
    G4int ncomponents;
    material_G4_POLYSTYRENE = new G4Material("G4_POLYSTYRENE_FOCAL", density = 1.05 * g / cm3, ncomponents = 2);
    material_G4_POLYSTYRENE->AddElement(G4Element::GetElement("C"), 8);
    material_G4_POLYSTYRENE->AddElement(G4Element::GetElement("H"), 8);
    material_G4_POLYSTYRENE->GetIonisation()->SetBirksConstant(0.126*mm/MeV);


    G4MaterialPropertiesTable *tab = new G4MaterialPropertiesTable();

    const G4int ntab = 31;
    tab->AddConstProperty("FASTTIMECONSTANT", 2.8*ns); // was 6
    // tab->AddConstProperty("SCINTILLATIONYIELD", 13.9/keV); // was 200/MEV nominal  10
    tab->AddConstProperty("SCINTILLATIONYIELD", 2/keV);//1.39/keV); // was 200/MEV nominal  10
    // tab->AddConstProperty("SCINTILLATIONYIELD", 200/MeV); // was 200/MEV nominal, should maybe be 13.9/keV
    tab->AddConstProperty("RESOLUTIONSCALE", 1.0);

    G4double opt_en[] =
      { 1.37760*eV, 1.45864*eV, 1.54980*eV, 1.65312*eV, 1.71013*eV, 1.77120*eV, 1.83680*eV, 1.90745*eV, 1.98375*eV, 2.06640*eV,
        2.10143*eV, 2.13766*eV, 2.17516*eV, 2.21400*eV, 2.25426*eV, 2.29600*eV, 2.33932*eV, 2.38431*eV, 2.43106*eV, 2.47968*eV,
        2.53029*eV, 2.58300*eV, 2.63796*eV, 2.69531*eV, 2.75520*eV, 2.81782*eV, 2.88335*eV, 2.95200*eV, 3.09960*eV, 3.54241*eV,
        4.13281*eV }; // 350 - 800 nm
    G4double scin_fast[] =
      { 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0003, 0.0008, 0.0032,
        0.0057, 0.0084, 0.0153, 0.0234, 0.0343, 0.0604, 0.0927, 0.1398, 0.2105, 0.2903,
        0.4122, 0.5518, 0.7086, 0.8678, 1.0000, 0.8676, 0.2311, 0.0033, 0.0012, 0.0000,
        0 };
    tab->AddProperty("FASTCOMPONENT", opt_en, scin_fast, ntab);

    G4double opt_r[] =
      { 1.5749, 1.5764, 1.5782, 1.5803, 1.5815, 1.5829, 1.5845, 1.5862, 1.5882, 1.5904,
        1.5914, 1.5924, 1.5935, 1.5947, 1.5959, 1.5972, 1.5986, 1.6000, 1.6016, 1.6033,
        1.6051, 1.6070, 1.6090, 1.6112, 1.6136, 1.6161, 1.6170, 1.6230, 1.62858, 1.65191,
        1.69165 };
    tab->AddProperty("RINDEX", opt_en, opt_r, ntab);

    G4double opt_abs[] =
      { 2.714*m, 3.619*m, 5.791*m, 4.343*m, 7.896*m, 5.429*m, 36.19*m, 17.37*m, 36.19*m, 5.429*m,
        13.00*m, 14.50*m, 16.00*m, 18.00*m, 16.50*m, 17.00*m, 14.00*m, 16.00*m, 15.00*m, 14.50*m,
        13.00*m, 12.00*m, 10.00*m, 8.000*m, 7.238*m, 4.000*m, 1.200*m, 0.500*m, 0.200*m, 0.200*m,
        0.100*m };
    tab->AddProperty("ABSLENGTH", opt_en, opt_abs, ntab);

    material_G4_POLYSTYRENE->SetMaterialPropertiesTable(tab);




    // const G4int nEntries = 50;
    // G4double photonEnergy[nEntries] =
    //     {2.00 * eV, 2.03 * eV, 2.06 * eV, 2.09 * eV, 2.12 * eV,
    //      2.15 * eV, 2.18 * eV, 2.21 * eV, 2.24 * eV, 2.27 * eV,
    //      2.30 * eV, 2.33 * eV, 2.36 * eV, 2.39 * eV, 2.42 * eV,
    //      2.45 * eV, 2.48 * eV, 2.51 * eV, 2.54 * eV, 2.57 * eV,
    //      2.60 * eV, 2.63 * eV, 2.66 * eV, 2.69 * eV, 2.72 * eV,
    //      2.75 * eV, 2.78 * eV, 2.81 * eV, 2.84 * eV, 2.87 * eV,
    //      2.90 * eV, 2.93 * eV, 2.96 * eV, 2.99 * eV, 3.02 * eV,
    //      3.05 * eV, 3.08 * eV, 3.11 * eV, 3.14 * eV, 3.17 * eV,
    //      3.20 * eV, 3.23 * eV, 3.26 * eV, 3.29 * eV, 3.32 * eV,
    //      3.35 * eV, 3.38 * eV, 3.41 * eV, 3.44 * eV, 3.47 * eV};
    // G4double scintilFast[nEntries] =
    //     {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    //      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    //      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    //      1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    //      1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};

    // const G4int ntab = 4;

    // G4double wls_Energy[] = {2.00 * eV, 2.87 * eV, 2.90 * eV,
    //                          3.47 * eV};

    // G4double rIndexPstyrene[] = {1.5, 1.5, 1.5, 1.5};
    // // G4double absorption1[] = {2. * m, 2. * m, 2. * m, 2. * m};
    // G4MaterialPropertiesTable* fMPTPStyrene = new G4MaterialPropertiesTable();
    // fMPTPStyrene->AddProperty("RINDEX", wls_Energy, rIndexPstyrene, ntab);
    // fMPTPStyrene->AddProperty("ABSLENGTH", opt_en, opt_abs, ntab);
    // // fMPTPStyrene->AddProperty("ABSLENGTH", wls_Energy, absorption1, ntab);
    // fMPTPStyrene->AddProperty("SCINTILLATIONCOMPONENT1", photonEnergy, scintilFast, nEntries);
    // fMPTPStyrene->AddConstProperty("SCINTILLATIONYIELD", 10. / keV);
    // fMPTPStyrene->AddConstProperty("RESOLUTIONSCALE", 1.0);
    // fMPTPStyrene->AddConstProperty("SCINTILLATIONTIMECONSTANT", 10. * ns);

    // material_G4_POLYSTYRENE->SetMaterialPropertiesTable(fMPTPStyrene);



    if (Verbosity() > 0)
    {
      cout << "PHG4FoCalDetector:  Making Scintillator material done." << endl;
    }
  }
  return material_G4_POLYSTYRENE;
}





//_______________________________________________________________________
G4Material*
PHG4FoCalDetector::GetPMMAMaterial()
{
  G4Material *material_PMMA = G4Material::GetMaterial("PMMA_FOCAL");
  if (!material_PMMA){
    if (Verbosity() > 0)
    {
      cout << "PHG4FoCalDetector: Making PMMA material..." << endl;
    }

    G4double density;
    G4int ncomponents;

    material_PMMA = new G4Material("PMMA_FOCAL", density = 1.18 * g / cm3, ncomponents = 3);
    material_PMMA->AddElement(G4Element::GetElement("C"), 5);
    material_PMMA->AddElement(G4Element::GetElement("H"), 8);
    material_PMMA->AddElement(G4Element::GetElement("O"), 2);

    const G4int nEntries = 31;

    G4double photonEnergy[nEntries] =
      { 1.37760*eV, 1.45864*eV, 1.54980*eV, 1.65312*eV, 1.71013*eV, 1.77120*eV, 1.83680*eV, 1.90745*eV, 1.98375*eV, 2.06640*eV,
        2.10143*eV, 2.13766*eV, 2.17516*eV, 2.21400*eV, 2.25426*eV, 2.29600*eV, 2.33932*eV, 2.38431*eV, 2.43106*eV, 2.47968*eV,
        2.53029*eV, 2.58300*eV, 2.63796*eV, 2.69531*eV, 2.75520*eV, 2.81782*eV, 2.88335*eV, 2.95200*eV, 3.09960*eV, 3.54241*eV,
        4.13281*eV };
    G4double refractiveIndexWLSfiber[nEntries] =
      { 1.4852, 1.4859, 1.4867, 1.4877, 1.4882, 1.4888, 1.4895, 1.4903, 1.4911, 1.4920,
        1.4924, 1.4929, 1.4933, 1.4938, 1.4943, 1.4948, 1.4954, 1.4960, 1.4966, 1.4973,
        1.4981, 1.4989, 1.4997, 1.5006, 1.5016, 1.5026, 1.5038, 1.5050, 1.5052, 1.5152,
        1.5306 };

    G4double absWLSfiber[nEntries] =
      { 0.414*m, 0.965*m, 2.171*m, 4.343*m, 1.448*m, 4.343*m, 14.48*m, 21.71*m, 8.686*m, 39.48*m,
        48.25*m, 54.29*m, 57.91*m, 54.29*m, 33.40*m, 31.02*m, 43.43*m, 43.43*m, 41.36*m, 39.48*m,
        37.76*m, 36.19*m, 36.19*m, 33.40*m, 31.02*m, 28.95*m, 25.55*m, 24.13*m, 21.71*m, 2.171*m,
        0.434*m };


    // Add entries into properties table
    G4MaterialPropertiesTable* mptWLSfiber = new G4MaterialPropertiesTable();
    mptWLSfiber->AddProperty("RINDEX",photonEnergy,refractiveIndexWLSfiber,nEntries);
    mptWLSfiber->AddProperty("ABSLENGTH",photonEnergy,absWLSfiber,nEntries);
    material_PMMA->SetMaterialPropertiesTable(mptWLSfiber);
    if (Verbosity() > 0)
    {
      cout << "PHG4FoCalDetector:  Making PMMA material done." << endl;
    }
  }
  return material_PMMA;
}

//_______________________________________________________________________
G4Material*
PHG4FoCalDetector::GetQuartzMaterial()
{
  if (Verbosity() > 0)
  {
    cout << "PHG4FoCalDetector: Making Scintillator material..." << endl;
  }

  G4MaterialPropertiesTable* mptWLSfiber = new G4MaterialPropertiesTable();

	G4double density;
	G4int ncomponents;

  // G4Material* material_Quartz = G4Material::GetMaterial("Quartz");
  G4Material *material_Quartz = new G4Material("Quartz", density = 2.200 * g / cm3, ncomponents = 2);
  material_Quartz->AddElement(G4Element::GetElement("Si"), 1);
  material_Quartz->AddElement(G4Element::GetElement("O"), 2);

  const G4int nEntriesQuartz = 279;

  G4double photonEnergyQuartz[nEntriesQuartz] =
    { 5.9040*eV, 5.7667*eV, 5.6356*eV, 5.5104*eV, 5.3906*eV, 5.2759*eV, 5.1660*eV, 5.0606*eV, 4.9594*eV, 4.8621*eV,
      4.7686*eV, 4.6786*eV, 4.5920*eV, 4.5085*eV, 4.4280*eV, 4.3503*eV, 4.2753*eV, 4.2029*eV, 4.1328*eV, 4.0651*eV,
      3.9995*eV, 3.9360*eV, 3.8745*eV, 3.8149*eV, 3.7571*eV, 3.7010*eV, 3.6466*eV, 3.5937*eV, 3.5424*eV, 3.4925*eV,
      3.4440*eV, 3.3968*eV, 3.3509*eV, 3.3062*eV, 3.2627*eV, 3.2204*eV, 3.1791*eV, 3.1388*eV, 3.0996*eV, 3.0613*eV,
      3.0240*eV, 2.9876*eV, 2.9520*eV, 2.9173*eV, 2.8834*eV, 2.8502*eV, 2.8178*eV, 2.7862*eV, 2.7552*eV, 2.7249*eV,
      2.6953*eV, 2.6663*eV, 2.6380*eV, 2.6102*eV, 2.5830*eV, 2.5564*eV, 2.5303*eV, 2.5047*eV, 2.4797*eV, 2.4551*eV,
      2.4311*eV, 2.4075*eV, 2.3843*eV, 2.3616*eV, 2.3393*eV, 2.3175*eV, 2.2960*eV, 2.2749*eV, 2.2543*eV, 2.2339*eV,
      2.2140*eV, 2.1944*eV, 2.1752*eV, 2.1562*eV, 2.1377*eV, 2.1194*eV, 2.1014*eV, 2.0838*eV, 2.0664*eV, 2.0493*eV,
      2.0325*eV, 2.0160*eV, 1.9997*eV, 1.9837*eV, 1.9680*eV, 1.9525*eV, 1.9373*eV, 1.9222*eV, 1.9074*eV, 1.8929*eV,
      1.8785*eV, 1.8644*eV, 1.8505*eV, 1.8368*eV, 1.8233*eV, 1.8100*eV, 1.7969*eV, 1.7839*eV, 1.7712*eV, 1.7463*eV,
      1.7220*eV, 1.6984*eV, 1.6755*eV, 1.6531*eV, 1.6314*eV, 1.6102*eV, 1.5895*eV, 1.5694*eV, 1.5498*eV, 1.5307*eV,
      1.5120*eV, 1.4938*eV, 1.4760*eV, 1.4586*eV, 1.4417*eV, 1.4251*eV, 1.4089*eV, 1.3931*eV, 1.3776*eV, 1.3625*eV,
      1.3477*eV, 1.3332*eV, 1.3190*eV, 1.3051*eV, 1.2915*eV, 1.2782*eV, 1.2651*eV, 1.2524*eV, 1.2398*eV, 1.2276*eV,
      1.2155*eV, 1.2037*eV, 1.1922*eV, 1.1808*eV, 1.1697*eV, 1.1587*eV, 1.1480*eV, 1.1375*eV, 1.1271*eV, 1.1170*eV,
      1.1070*eV, 1.0972*eV, 1.0876*eV, 1.0781*eV, 1.0688*eV, 1.0597*eV, 1.0507*eV, 1.0419*eV, 1.0332*eV, 1.0247*eV,
      1.0163*eV, 1.0080*eV, 0.9999*eV, 0.9919*eV, 0.9840*eV, 0.9763*eV, 0.9686*eV, 0.9611*eV, 0.9537*eV, 0.9464*eV,
      0.9393*eV, 0.9322*eV, 0.9253*eV, 0.9184*eV, 0.9116*eV, 0.9050*eV, 0.8984*eV, 0.8920*eV, 0.8856*eV, 0.8793*eV,
      0.8731*eV, 0.8670*eV, 0.8610*eV, 0.8551*eV, 0.8492*eV, 0.8434*eV, 0.8377*eV, 0.8321*eV, 0.8266*eV, 0.8211*eV,
      0.8157*eV, 0.8104*eV, 0.8051*eV, 0.7999*eV, 0.7948*eV, 0.7897*eV, 0.7847*eV, 0.7798*eV, 0.7749*eV, 0.7701*eV,
      0.7653*eV, 0.7606*eV, 0.7560*eV, 0.7514*eV, 0.7469*eV, 0.7424*eV, 0.7380*eV, 0.7336*eV, 0.7293*eV, 0.7251*eV,
      0.7208*eV, 0.7167*eV, 0.7126*eV, 0.7085*eV, 0.7045*eV, 0.7005*eV, 0.6965*eV, 0.6926*eV, 0.6888*eV, 0.6850*eV,
      0.6812*eV, 0.6775*eV, 0.6738*eV, 0.6702*eV, 0.6666*eV, 0.6630*eV, 0.6595*eV, 0.6560*eV, 0.6525*eV, 0.6491*eV,
      0.6458*eV, 0.6424*eV, 0.6391*eV, 0.6358*eV, 0.6326*eV, 0.6294*eV, 0.6262*eV, 0.6230*eV, 0.6199*eV, 0.6168*eV,
      0.6138*eV, 0.6108*eV, 0.6078*eV, 0.6048*eV, 0.6019*eV, 0.5990*eV, 0.5961*eV, 0.5932*eV, 0.5904*eV, 0.5876*eV,
      0.5848*eV, 0.5821*eV, 0.5794*eV, 0.5767*eV, 0.5740*eV, 0.5714*eV, 0.5687*eV, 0.5661*eV, 0.5636*eV, 0.5610*eV,
      0.5585*eV, 0.5560*eV, 0.5535*eV, 0.5510*eV, 0.5486*eV, 0.5462*eV, 0.5438*eV, 0.5414*eV, 0.5391*eV, 0.5367*eV,
      0.5344*eV, 0.5321*eV, 0.5298*eV, 0.5276*eV, 0.5254*eV, 0.5231*eV, 0.5209*eV, 0.5188*eV, 0.5166*eV, 0.5145*eV,
      0.5123*eV, 0.5102*eV, 0.5081*eV, 0.5061*eV, 0.5040*eV, 0.5020*eV, 0.4999*eV, 0.4979*eV, 0.4959*eV};
  G4double refractiveIndexQuartz[nEntriesQuartz] =
    { 1.5384, 1.5332, 1.5285, 1.5242, 1.5202, 1.5166, 1.5133, 1.5103, 1.5074, 1.5048,
      1.5024, 1.5001, 1.4980, 1.4960, 1.4942, 1.4924, 1.4908, 1.4892, 1.4878, 1.4864,
      1.4851, 1.4839, 1.4827, 1.4816, 1.4806, 1.4796, 1.4787, 1.4778, 1.4769, 1.4761,
      1.4753, 1.4745, 1.4738, 1.4731, 1.4725, 1.4719, 1.4713, 1.4707, 1.4701, 1.4696,
      1.4691, 1.4686, 1.4681, 1.4676, 1.4672, 1.4668, 1.4663, 1.4660, 1.4656, 1.4652,
      1.4648, 1.4645, 1.4641, 1.4638, 1.4635, 1.4632, 1.4629, 1.4626, 1.4623, 1.4621,
      1.4618, 1.4615, 1.4613, 1.4610, 1.4608, 1.4606, 1.4603, 1.4601, 1.4599, 1.4597,
      1.4595, 1.4593, 1.4591, 1.4589, 1.4587, 1.4586, 1.4584, 1.4582, 1.4580, 1.4579,
      1.4577, 1.4576, 1.4574, 1.4572, 1.4571, 1.4570, 1.4568, 1.4567, 1.4565, 1.4564,
      1.4563, 1.4561, 1.4560, 1.4559, 1.4558, 1.4556, 1.4555, 1.4554, 1.4553, 1.4551,
      1.4549, 1.4546, 1.4544, 1.4542, 1.4540, 1.4539, 1.4537, 1.4535, 1.4533, 1.4531,
      1.4530, 1.4528, 1.4527, 1.4525, 1.4523, 1.4522, 1.4520, 1.4519, 1.4518, 1.4516,
      1.4515, 1.4513, 1.4512, 1.4511, 1.4509, 1.4508, 1.4507, 1.4505, 1.4504, 1.4503,
      1.4502, 1.4500, 1.4499, 1.4498, 1.4497, 1.4496, 1.4494, 1.4493, 1.4492, 1.4491,
      1.4490, 1.4489, 1.4487, 1.4486, 1.4485, 1.4484, 1.4483, 1.4482, 1.4481, 1.4479,
      1.4478, 1.4477, 1.4476, 1.4475, 1.4474, 1.4473, 1.4471, 1.4470, 1.4469, 1.4468,
      1.4467, 1.4466, 1.4465, 1.4464, 1.4462, 1.4461, 1.4460, 1.4459, 1.4458, 1.4457,
      1.4455, 1.4454, 1.4453, 1.4452, 1.4451, 1.4450, 1.4449, 1.4447, 1.4446, 1.4445,
      1.4444, 1.4443, 1.4441, 1.4440, 1.4439, 1.4438, 1.4437, 1.4435, 1.4434, 1.4433,
      1.4432, 1.4431, 1.4429, 1.4428, 1.4427, 1.4426, 1.4424, 1.4423, 1.4422, 1.4420,
      1.4419, 1.4418, 1.4417, 1.4415, 1.4414, 1.4413, 1.4411, 1.4410, 1.4409, 1.4407,
      1.4406, 1.4405, 1.4403, 1.4402, 1.4401, 1.4399, 1.4398, 1.4397, 1.4395, 1.4394,
      1.4392, 1.4391, 1.4389, 1.4388, 1.4387, 1.4385, 1.4384, 1.4382, 1.4381, 1.4379,
      1.4378, 1.4376, 1.4375, 1.4373, 1.4372, 1.4370, 1.4369, 1.4367, 1.4366, 1.4364,
      1.4363, 1.4361, 1.4360, 1.4358, 1.4357, 1.4355, 1.4353, 1.4352, 1.4350, 1.4349,
      1.4347, 1.4345, 1.4344, 1.4342, 1.4340, 1.4339, 1.4337, 1.4335, 1.4334, 1.4332,
      1.4330, 1.4328, 1.4327, 1.4325, 1.4323, 1.4322, 1.4320, 1.4318, 1.4316, 1.4314,
      1.4313, 1.4311, 1.4309, 1.4307, 1.4305, 1.4304, 1.4302, 1.4300, 1.4298 };
  mptWLSfiber->AddProperty("RINDEX",photonEnergyQuartz,refractiveIndexQuartz,nEntriesQuartz);


  const G4int nEntries_Quartz = 14;
  G4double PhotonEnergy_Quartz[nEntries_Quartz] =
    { 2.21*eV, 2.30*eV, 2.38*eV, 2.48*eV, 2.58*eV, 2.70*eV, 2.82*eV, 2.95*eV, 3.10*eV, 3.26*eV,
      3.44*eV, 3.65*eV, 3.88*eV, 4.13*eV };
  G4double Quartz_Abs[nEntries_Quartz] =
    { 550.7*mm, 530.7*mm, 590.1*mm, 490.7*mm, 470.7*mm, 520.3*mm, 500.0*mm, 470.7*mm, 450.5*mm, 270.5*mm,
      190.1*mm,  60.9*mm,  10.6*mm,   4.0*mm};
  mptWLSfiber->AddProperty("ABSLENGTH",  PhotonEnergy_Quartz, Quartz_Abs,  nEntries_Quartz);

  material_Quartz->SetMaterialPropertiesTable(mptWLSfiber);

  if (Verbosity() > 0)
  {
    cout << "PHG4FoCalDetector:  Making Scintillator material done." << endl;
  }

  return material_Quartz;
}

//_____________________________________________________________________________
void PHG4FoCalDetector::SurfaceTable(G4LogicalVolume *vol,  G4String  name) {

  G4OpticalSurface *surface = new G4OpticalSurface(name + "_OPSurf");

  new G4LogicalSkinSurface( name + "_LogSkinSurf", vol, surface);

  // surface->SetType(dielectric_dielectric);
  surface->SetType(dielectric_metal);
  // surface->SetType(dielectric_LUT);
  surface->SetFinish(polished);
  // surface->SetModel(LUT);
  surface->SetModel(glisur);

  //crystal optical surface

  //surface material
  const G4int ntab = 2;
  G4double opt_en[] = {1.451*eV, 3.645*eV}; // 350 - 800 nm
  G4double reflectivity[] = {0.99, 0.99};
  G4double efficiency[] = {0.99, 0.99};
  G4MaterialPropertiesTable *surfmat = new G4MaterialPropertiesTable();
  surfmat->AddProperty("REFLECTIVITY", opt_en, reflectivity, ntab);
  surfmat->AddProperty("EFFICIENCY", opt_en, efficiency, ntab);
  surface->SetMaterialPropertiesTable(surfmat);
  //csurf->DumpInfo();

}//SurfaceTable


//_____________________________________________________________________________
void PHG4FoCalDetector::MakeBoundaryFibers(G4VPhysicalVolume *fromVol, G4VPhysicalVolume *toVol, G4String name) {

  //optical boundary between the fromVol and optical photons detector

  G4OpticalSurface *surf = new G4OpticalSurface(name +"_OpDetS1");
  surf->SetType(dielectric_dielectric); // photons go to the detector, must have rindex defined
  // surf->SetType(dielectric_metal); // photon is absorbed when reaching the detector, no material rindex required
  //surf->SetFinish(ground);
  surf->SetFinish(polished);
  //surf->SetModel(unified);
  surf->SetModel(glisur);

  new G4LogicalBorderSurface(name +"OpDetBS", fromVol, toVol, surf);

  const G4int ntab = 2;
  G4double opt_en[] = {1.051*eV, 4.545*eV}; // 350 - 800 nm
  //G4double reflectivity[] = {0., 0.};
  G4double reflectivity[] = {0.0, 0.0};
  //G4double reflectivity[] = {1., 1.};
  // G4double efficiency[] = {1., 1.};
  G4double transmittance[] = { 1.0, 1.0 };
  G4double RefractiveIndexBoundary[] = { 1.5, 1.5};

  G4MaterialPropertiesTable *surfmat = new G4MaterialPropertiesTable();
  // surfmat->AddProperty("EFFICIENCY", opt_en, efficiency, ntab);
  surfmat->AddProperty("RINDEX", opt_en, RefractiveIndexBoundary, ntab);
  surfmat->AddProperty("REFLECTIVITY", opt_en, reflectivity, ntab);
  surfmat->AddProperty("TRANSMITTANCE", opt_en, transmittance, ntab);
  surf->SetMaterialPropertiesTable(surfmat);



}//MakeBoundary



int PHG4FoCalDetector::ParseParametersFromTable()
{
  //Open the datafile, if it won't open return an error
  ifstream istream_mapping;
  istream_mapping.open(_mapping_tower_file);
  if (!istream_mapping.is_open())
  {
    std::cout << "ERROR in PHG4FoCalDetector: Failed to open mapping file " << _mapping_tower_file << std::endl;
    gSystem->Exit(1);
  }
  /* loop over lines in file */
  std::string line_mapping;
  while (getline(istream_mapping, line_mapping))
  {
    /* Skip lines starting with / including a '#' */
    if (line_mapping.find("#") != std::string::npos)
    {
      if (Verbosity() > 0)
      {
        std::cout << "PHG4FoCalDetector: SKIPPING line in mapping file: " << line_mapping << std::endl;
      }
      continue;
    }

    std::istringstream iss(line_mapping);

    /* If line starts with keyword Tower, add to tower positions */
    if (line_mapping.find("Tower ") != std::string::npos)
    {
      unsigned idx_j, idx_k, idx_l;
      G4double pos_x, pos_y, pos_z;
      G4double size_x, size_y, size_z;
      G4double rot_x, rot_y, rot_z;
      G4double dummy;
      std::string dummys;

      /* read string- break if error */
      if (!(iss >> dummys >> dummy >> idx_j >> idx_k >> idx_l >> pos_x >> pos_y >> pos_z >> size_x >> size_y >> size_z >> rot_x >> rot_y >> rot_z))
      {
        std::cout << "ERROR in PHG4FoCalDetector: Failed to read line in mapping file " << _mapping_tower_file << std::endl;
        gSystem->Exit(1);
      }

      /* Construct unique name for tower */
      /* Mapping file uses cm, this class uses mm for length */
      std::ostringstream towername;
      towername.str("");
      towername << _towerlogicnameprefix << "_j_" << idx_j << "_k_" << idx_k;
      // cout << _towerlogicnameprefix << "_j_" << idx_j << "_k_" << idx_k << endl;
      /* Add Geant4 units */
      pos_x = pos_x * cm;
      pos_y = pos_y * cm;
      pos_z = pos_z * cm;

      /* insert tower into tower map */
      towerposition tower_new;
      tower_new.x = pos_x;
      tower_new.y = pos_y;
      tower_new.z = pos_z;
      tower_new.idx_j = idx_j;
      tower_new.idx_k = idx_k;
      m_TowerPostionMap.insert(make_pair(towername.str(), tower_new));
    }
    else
    {
      /* If this line is not a comment and not a tower, save parameter as string / value. */
      std::string parname;
      G4double parval;

      /* read string- break if error */
      if (!(iss >> parname >> parval))
      {
        std::cout << "ERROR in PHG4ForwardHcalDetector: Failed to read line in mapping file " << _mapping_tower_file << std::endl;
        gSystem->Exit(1);
      }

      m_GlobalParameterMap.insert(make_pair(parname, parval));
    }
  }


  //Update member variables for global parameters based on parsed parameter file
  std::map<string, G4double>::iterator parit;


  parit = m_GlobalParameterMap.find("Gtype");
  if (parit != m_GlobalParameterMap.end())
  {
    _tower_type = parit->second;
  }

  cout << __LINE__ << endl;
  parit = m_GlobalParameterMap.find("Gtower_dx");
  if (parit != m_GlobalParameterMap.end())
  {
    _tower_dx = parit->second * cm;
    m_SteppingAction->SetTowerSize(_tower_dx);
  }
  cout << __LINE__ << endl;

  parit = m_GlobalParameterMap.find("Gtower_readout");
  if (parit != m_GlobalParameterMap.end())
  {
    _tower_readout = parit->second * cm;
    m_SteppingAction->SetTowerReadout(_tower_readout);
  }

  parit = m_GlobalParameterMap.find("Gtower_dy");
  if (parit != m_GlobalParameterMap.end())
  {
    _tower_dy = parit->second * cm;
  }

  parit = m_GlobalParameterMap.find("Gtower_dz");
  if (parit != m_GlobalParameterMap.end())
  {
    _tower_dz = parit->second * cm;
  }

  // new start
  parit = m_GlobalParameterMap.find("Scint_Diam");
  if (parit != m_GlobalParameterMap.end())
  {
    _scintFiber_diam = parit->second * cm;
  }

  parit = m_GlobalParameterMap.find("Cerenkov_Diam");
  if (parit != m_GlobalParameterMap.end())
  {
    _cerenkovFiber_diam = parit->second * cm;
  }

  parit = m_GlobalParameterMap.find("Cerenkov_Material");
  if (parit != m_GlobalParameterMap.end())
  {
    _cerenkovFiber_material = parit->second;
  }

  parit = m_GlobalParameterMap.find("NotchCutout");
  if (parit != m_GlobalParameterMap.end())
  {
    _tower_makeNotched = parit->second;
  }

  parit = m_GlobalParameterMap.find("Absorber_Material");
  if (parit != m_GlobalParameterMap.end())
  {
    _absorber_Material = parit->second;
  }
  // new end

  parit = m_GlobalParameterMap.find("Gr1_inner");
  if (parit != m_GlobalParameterMap.end())
  {
    _rMin1 = parit->second * cm;
  }

  parit = m_GlobalParameterMap.find("Gr1_outer");
  if (parit != m_GlobalParameterMap.end())
  {
    _rMax1 = parit->second * cm;
    m_SteppingAction->SetDetectorSize(_rMax1);
  }

  parit = m_GlobalParameterMap.find("Gr2_inner");
  if (parit != m_GlobalParameterMap.end())
  {
    _rMin2 = parit->second * cm;
  }

  parit = m_GlobalParameterMap.find("Gr2_outer");
  if (parit != m_GlobalParameterMap.end())
  {
    _rMax2 = parit->second * cm;
  }

  parit = m_GlobalParameterMap.find("Gdz");
  if (parit != m_GlobalParameterMap.end())
  {
    _dZ = parit->second * cm;
  }

  parit = m_GlobalParameterMap.find("Gx0");
  if (parit != m_GlobalParameterMap.end())
  {
    _place_in_x = parit->second * cm;
  }

  parit = m_GlobalParameterMap.find("Gy0");
  if (parit != m_GlobalParameterMap.end())
  {
    _place_in_y = parit->second * cm;
  }

  parit = m_GlobalParameterMap.find("Gz0");
  if (parit != m_GlobalParameterMap.end())
  {
    _place_in_z = parit->second * cm;
  }

  parit = m_GlobalParameterMap.find("Center_Offset_x");
  if (parit != m_GlobalParameterMap.end())
  {
    _center_offset_x = parit->second * cm;
  }

  parit = m_GlobalParameterMap.find("Center_Offset_y");
  if (parit != m_GlobalParameterMap.end())
  {
    _center_offset_y = parit->second * cm;
  }
  parit = m_GlobalParameterMap.find("Quadratic_Detector");
  if (parit != m_GlobalParameterMap.end())
  {
    _quadratic_detector = parit->second;
  }

  parit = m_GlobalParameterMap.find("Grot_x");
  if (parit != m_GlobalParameterMap.end())
  {
    _rot_in_x = parit->second;
  }

  parit = m_GlobalParameterMap.find("Grot_y");
  if (parit != m_GlobalParameterMap.end())
  {
    _rot_in_y = parit->second;
  }

  parit = m_GlobalParameterMap.find("Grot_z");
  if (parit != m_GlobalParameterMap.end())
  {
    _rot_in_z = parit->second;
  }

  return 0;
}
