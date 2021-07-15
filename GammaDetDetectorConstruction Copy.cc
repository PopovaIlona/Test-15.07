
#include "GammaDetDetectorConstruction.hh"

#include "G4Material.hh"
#include "G4Element.hh"
#include "G4SystemOfUnits.hh"
#include "G4NistManager.hh"
#include "G4VisAttributes.hh"
#include "G4PhysicalConstants.hh"

#include "G4Box.hh"
#include "G4Trap.hh"
#include "G4Tubs.hh"
#include "G4GenericTrap.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"

#include "G4Region.hh"
#include "G4RegionStore.hh"
#include "G4UserLimits.hh"

#include <G4SDManager.hh>
#include <G4PSEnergyDeposit.hh>
#include "EnergyTimeSD.hh"
#include <G4GlobalMagFieldMessenger.hh>
 Внесение каких-то изменений тест 
 
 Внесение еще каких-то изменений тест
 
 Внесение изменений в ветке 1
 И еще немного ветка 1
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

GammaDetDetectorConstruction::GammaDetDetectorConstruction()
 : G4VUserDetectorConstruction()
{
    ncr_x=3;
    ncr_y=3;
   
    Sizecr_x=2.5*cm;
    Sizecr_y=2.5*cm;
    Sizecr_z=14*cm;
   
    SizeSi_x=5*mm;
    SizeSi_y=5*mm;
    SizeSi_z=1*mm;
   
    
    SizePl_x=ncr_x*Sizecr_x;
    SizePl_y=ncr_y*Sizecr_y;
    SizePl_z=1*cm;
   
    SizeAl_x=ncr_x*Sizecr_x;
    SizeAl_y=ncr_y*Sizecr_y;
    SizeAl_z=3*mm;
   
    SizeVac_x=ncr_x*Sizecr_x;
    SizeVac_y=ncr_y*Sizecr_y;
    SizeVac_z=3*cm;
   
    PaperThick=0.2*mm;
    Detgap=1*mm;
   
      
  fExpHall_x = ncr_x*Sizecr_x+4*cm;
  fExpHall_y = ncr_y*Sizecr_y+4*cm;
  fExpHall_z = Sizecr_z+SizeSi_z+SizePl_z+SizeAl_z+SizeVac_z+Detgap+4*cm;
  

  
   
   
    PosVac_z=-fExpHall_z/2+2*cm+SizeVac_z/2;
	PosAl_z=PosVac_z+SizeAl_z/2+SizeVac_z/2;
   

    PosPl_z=PosAl_z+Detgap+SizePl_z/2+SizeAl_z/2;
   
    Poscr_z=PosPl_z+SizePl_z/2+Sizecr_z/2;
 
    PosSi_z=Poscr_z+Sizecr_z/2+SizeSi_z/2;
  
 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

GammaDetDetectorConstruction::~GammaDetDetectorConstruction(){;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* GammaDetDetectorConstruction::Construct()
{
 G4NistManager* nist = G4NistManager::Instance();
 // ------------- Materials -------------
 
 G4double a, z, density, temperature, pressure;
 G4int nelements;

 G4Element* N = new G4Element("Nitrogen", "N", z = 7, a = 14.01*g/mole);
 G4Element* O = new G4Element("Oxygen", "O", z = 8, a = 16.00*g / mole);
 G4Element* Si = new G4Element("Silicon", "Si" , z= 14., a= 28.09*g/mole);
 G4Element* H = new G4Element("Hydrogen", "H", z = 1, a = 1.01*g/mole);
 G4Element* Na = new G4Element( "Na", "Na" , z = 11,  a = 23.0*g/mole);
 G4Element* Ba = new G4Element("Barium","Ba",z=56, a= 137.327*g/mole);
 G4Element* K = new G4Element("Potassium","K",z=19, a=39.098*g/mole);
 G4Element* B = new G4Element("Boron", "B",z=5, a=10.81*g/mole);
 G4Element* As = new G4Element("Arsen",  "As", z=33, a= 74.92160*g/mole);
 G4Element* C = new G4Element("Carbon", "C", z = 6, a = 12.0116 *g/mole);
 G4Element* Pb = new G4Element("Lead", "Pb", z = 82, a = 207.19*g / mole);
 G4Element* F = new G4Element("Fluorine", "F", z = 9, a = 19.00*g / mole);
 G4Element* Li = new G4Element( "Lithium"   , "Li" ,  z = 3,  a =  7.0 *g/mole );
 // Air
 //
 
 G4Material* air = new G4Material("Air", density = 1.29*mg / cm3, nelements = 2);
 air->AddElement(N, 70.*perCent);
 air->AddElement(O, 30.*perCent);
 
 // Aluminium
 //  
 G4Material* Al = new G4Material("Aluminum", z = 13., a = 26.98*g / mole, density = 2.7*g / cm3);
 
 // Silicon
 //  
 G4Material* Sil = new G4Material("Silicon", z = 14., a = 28.09*g / mole, density = 2.330*g / cm3);
 
 
 //
 // vacuum
 G4Material*Galactic = new G4Material("Galactic", z = 1., a = 1.01*g / mole, density = 1.e-25*g / cm3, kStateGas,
  temperature = 0.1*kelvin, pressure = 1.e-19*pascal);

 
 // other way CH2O
 G4Material* CH2O = new G4Material("CH2O", density = 1.41*g / cm3, nelements = 3);
 CH2O->AddElement(C, 1);
 CH2O->AddElement(H, 2);
 CH2O->AddElement(O, 1);
 
  //PbF2

 G4Material * PbF2 = new G4Material("PbF2", density = 7.77*g / cm3, nelements = 2);
 PbF2->AddElement(Pb, 1);
 PbF2->AddElement(F, 2);

//
// ------------ Generate & Add Material Properties Table ------------
//
 G4double photonEnergy[] =
            { 2.034*eV, 2.068*eV, 2.103*eV, 2.139*eV,
              2.177*eV, 2.216*eV, 2.256*eV, 2.298*eV,
              2.341*eV, 2.386*eV, 2.433*eV, 2.481*eV,
              2.532*eV, 2.585*eV, 2.640*eV, 2.697*eV,
              2.757*eV, 2.820*eV, 2.885*eV, 2.954*eV,
              3.026*eV, 3.102*eV, 3.181*eV, 3.265*eV,
              3.353*eV, 3.446*eV, 3.545*eV, 3.649*eV,
              3.760*eV, 3.877*eV, 4.002*eV, 4.136*eV };

  const G4int nEntries = sizeof(photonEnergy)/sizeof(G4double);

//
// PbF2
//
  G4double refractiveIndex1[] =
            { 1.7637, 1.7650,  1.7663, 1.7676,  1.7691,
              1.7707,  1.7723, 1.7740,  1.7759, 1.7779,
              1.7800, 1.7823, 1.7847,   1.7874, 1.7902,
              1.7932, 1.7965, 1.8002, 1.8040, 1.8083,
              1.8130, 1.8182,  1.8239, 1.8303,  1.8373,
              1.8453, 1.8545,  1.8649, 1.8772,  1.8920,
              1.9107,   1.9366};

  assert(sizeof(refractiveIndex1) == sizeof(photonEnergy));
    G4MaterialPropertiesTable* myMPT1 = new G4MaterialPropertiesTable();

  myMPT1->AddProperty("RINDEX",       photonEnergy, refractiveIndex1,nEntries)
        ->SetSpline(true);
		 G4cout << "PbF2 G4MaterialPropertiesTable" << G4endl;
  myMPT1->DumpTable();

  PbF2->SetMaterialPropertiesTable(myMPT1);
  
  //
// Air
//
  G4double refractiveIndex2[] =
            { 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
              1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
              1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
              1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
              1.00, 1.00, 1.00, 1.00 };

  G4MaterialPropertiesTable* myMPT2 = new G4MaterialPropertiesTable();
  myMPT2->AddProperty("RINDEX", photonEnergy, refractiveIndex2, nEntries);

  G4cout << "Air G4MaterialPropertiesTable" << G4endl;
  myMPT2->DumpTable();

  air->SetMaterialPropertiesTable(myMPT2);
  
   G4double transmissionsS[nEntries] =

  { 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
   0, 0, 0 };

  G4double refractiveIndexS[nEntries] =

  { 1.988,     2.452,     3.12,     4.087,     4.888,     5.02,     5.01,     
  5.016,     5.065,     5.156,     5.296,     5.61,     6.522,     6.709,    
  6.062,     5.57,     5.222,     4.961,     4.753,     4.583,     4.442,   
  4.32,     4.215,     4.123,     4.042,     3.969,     3.906,     3.847,   
  3.796,     3.752,     3.714,     3.673
  };


 G4double absorption[nEntries] =
           {3.448*nm,  4.082*nm,  6.329*nm,  9.174*nm, 12.346*nm, 13.889*nm,
           15.152*nm, 17.241*nm, 18.868*nm, 20.000*nm, 26.316*nm, 35.714*nm,
           45.455*nm, 47.619*nm, 52.632*nm, 52.632*nm, 55.556*nm, 52.632*nm,
           52.632*nm, 47.619*nm, 45.455*nm, 41.667*nm, 37.037*nm, 33.333*nm,
           30.000*nm, 28.500*nm, 27.000*nm, 24.500*nm, 22.000*nm, 19.500*nm,
           17.500*nm, 14.500*nm };

 G4MaterialPropertiesTable* tabS = new G4MaterialPropertiesTable();

  tabS->AddProperty("RINDEX", photonEnergy, refractiveIndexS, nEntries)
     ->SetSpline(true);

  tabS->AddProperty("TRANSMISSION", photonEnergy, transmissionsS, nEntries)
     ->SetSpline(true);
  tabS->AddProperty("ABSLENGTH",    photonEnergy, absorption,     nEntries)
         ->SetSpline(true);

  Sil->SetMaterialPropertiesTable(tabS);
 
 // ------------- Volumes --------------

     // The experimental Hall
 //
 G4Box* expHall_box = new G4Box("World", fExpHall_x/2, fExpHall_y/2, fExpHall_z/2);

 G4LogicalVolume* expHall_log
  = new G4LogicalVolume(expHall_box, air, "World", 0, 0, 0);

 G4VPhysicalVolume* expHall_phys
  = new G4PVPlacement(0, G4ThreeVector(), expHall_log, "World", 0, false, 0);

 // invisible world
 //expHall_log->SetVisAttributes(G4VisAttributes::Invisible);

  
      // The Vacuum

  G4Box* vacuum_box = new G4Box("vacuum",SizeVac_x/2,SizeVac_y/2,SizeVac_z/2);

  G4LogicalVolume* vacuum_log
    = new G4LogicalVolume(vacuum_box, Galactic,"vacuum",0,0,0);

G4VPhysicalVolume* vacuum_phys =
     new G4PVPlacement(0,G4ThreeVector(0,0,PosVac_z),vacuum_log,"vacuum",

                       expHall_log,false,0);
					   
// The AL

  G4Box* al_box = new G4Box("al",SizeAl_x/2,SizeAl_y/2,SizeAl_z/2);

  G4LogicalVolume* al_log
    = new G4LogicalVolume(al_box, Al,"al",0,0,0);

G4VPhysicalVolume* al_phys =
     new G4PVPlacement(0,G4ThreeVector(0,0,PosAl_z),al_log,"al",

                       expHall_log,false,0);					   
					   
// The Pl

  G4Box* pl_box = new G4Box("pl",SizePl_x/2,SizePl_y/2,SizePl_z/2);

  G4LogicalVolume* pl_log
    = new G4LogicalVolume(pl_box, CH2O,"pl",0,0,0);

G4VPhysicalVolume* pl_phys =
     new G4PVPlacement(0,G4ThreeVector(0,0,PosPl_z),pl_log,"pl",

                       expHall_log,false,0);	
					   
// The cr

  G4Box* cr_box = new G4Box("cr",Sizecr_x/2-0.01*mm,Sizecr_y/2-0.01*mm,Sizecr_z/2);

  G4LogicalVolume* cr_log
    = new G4LogicalVolume(cr_box, PbF2,"cr",0,0,0);
	
// The Sil

  G4Box* Sil_box = new G4Box("Sil",SizeSi_x/2,SizeSi_y/2,SizeSi_z/2);

  G4LogicalVolume* Sil_log
    = new G4LogicalVolume(Sil_box, Sil,"Sil",0,0,0);	
	
	for (G4int i=1; i<=ncr_x; i++)
	{
	for (G4int j=1; j<=ncr_y; j++)
	{ G4double pos_x=-SizeAl_x/2-(Sizecr_x/2)+Sizecr_x*i;
	G4double pos_y=-SizeAl_y/2-(Sizecr_y/2)+Sizecr_y*j;
	G4VPhysicalVolume* cr_phys =
     new G4PVPlacement(0,G4ThreeVector(pos_x,pos_y,Poscr_z),cr_log,"cr",

                       expHall_log,false,0);
G4VPhysicalVolume* Sil_phys =
     new G4PVPlacement(0,G4ThreeVector(pos_x,pos_y,PosSi_z),Sil_log,"Sil",

                       expHall_log,false,0);						   
	}
	}

 // G4Colour color3(0, 0, 1);
 // G4VisAttributes* vacvistat = new  G4VisAttributes(color3);
 // vacuum_log->SetVisAttributes(vacvistat);

    EnergyTimeSD* absorberET = new EnergyTimeSD("absorberET");
	G4SDManager* SDman = G4SDManager::GetSDMpointer();
    SDman->AddNewDetector(absorberET);
    Sil_log->SetSensitiveDetector(absorberET);
 

 
 
 
  return expHall_phys;
}


void GammaDetDetectorConstruction::ConstructSDandField()
{

}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
