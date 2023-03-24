This is the note for studying the RE01 example.

# Detector construction
In this part, it has many geometric attributes. Its declaration and initialization are included from other files. Its initializer only initializes these geometric attributes. The geometry is defined as follows: for each object, define a solid, logical volume, and physical volume. Then instantiate a `G4VisAttributes` defining color, and assign it to the logical volume. Then instantiate a `G4Region` object, assign attributes to it through `defaultRInfo`. It looks like this:
* `G4PhysicalVolume`
    * `G4LogicalVolume`
        * `G4Solid`
        * `G4VisAttributes`
* `G4Region`
    * `RE01RegionInformation`
In the specification, 
* The world volume `experimentalHall_phys` is a box, of size `fExpHall_x`, `fExpHall_y`, `fExpHall_z` equaling to 600cm, 600cm, 600cm, respectively.
* The tracer `tracker_phys` is a tub, of size `fTrkTubs_rmin`, `fTrkTubs_rmax`, `fTrkTubs_dz` equaling to 20cm, 50cm, 100cm, respectively. The $\phi$ angle spans 360 degree. The material is `arGas`. It is placed at the center of the world volume.
* The `trackerLayer_tubs` is a tub of the same size as the tub above. But the material is `silicon`. It involves special classes `G4VParameterisation` and `G4VParameterised`.
* The `calorimeter_rubs` is a tub of inner radius, outer radius, length of 50cm, 300cm, and 200cm, respectively. It spans all 360 degree. It is of material `scinti`, which is composed of carbon and hydrogen. It is also placed at the center of the world wolume.
* The `caloLayer_tubs` is also a tub of the same shape as above. It is made of lead. It is parameterized to the `fCalorimeter_log`.
Now, I'de like to redraw the relationships between these volumes.
* `experimentalHall_phys` / `experimentalHall_log` / `experimentalHall_box`. `experimentalHallVisAtt` is associated to the logical volume. and `defaultRegion` is defined standalone.
    * `tracker_phys` (not embodied) / `tracker_log` / `tracker_tubs`. `tracker_logVisAtt` is assigned to the logical volume.  `trackerRegion` is defined standalone, and `trackerInfo` and `tracker_log` are associated to it.
    * ~ / `fTrackerLayer_log` / `trackerLayer_tubs`. It is noted that the logical volume is not placed in the world volume as a physical volume. `trackerLayer_logVisAtt` is associated to the logical volume. The logical volume is somehow associated to the logical volume of tracker `tracker_log` through the function call: `new G4PVParameterised("trackerLayer_phys", fTrackerLayer_log, tracker_log, kXAxis, fNotrkLayers, trackerParam);`, where `trackerParam` is an object of `G4PVParameterisation`.
    * `calorimeter_phys` (not embodied) / `calorimeter_log` / `calorimeter_tubs`. `calorimeter_logVisATT` is assigned to the logical volume. `calorimeterRegion` is defined standalone, and `calorimeterInfo` and `fCalorimeter_log` are associated to it.
    * ~ / `caloLayer_log` / `caloLayer_tubs`. There is no physical volume. `caloLayer_logVisAtt` is associated to the logical volume. The logical volume is somehow asoscated to the logical volume of the calorimeter `calorimeter_phys` through the function call: `new G4PVParameterised("caloLayer_phys", caloLayer_log, fCalorimeter_log, kXAxis, fNocaloLayers, calorimeterParam)`. `calorimeterParam` is an object of `G4VParameterisation`.
The two classes `RE01TrackerParametrisation` and `RE01CaloriParameterisation` does two things: 1) translate to origin. 2) Set the calorimeterLayer with new size.

Updates on Monday: I figured out that the geometry is a layer-by-layer one. An experimental hall contains two groups of detectors: tracker and calorimeter. 1) Tracker. As specified by `RE01TrackerParametrisation`, The tracker layers are of sizes (rmin, rmax, length) of (25, 25.5, 25), (30, 30.5, 30), (35, 35.5, 35), (40, 40.5, 40), (45, 45.5, 45). The tracker is made of `arGas`, and the tracker layers are made of `silicon`. 2) Calorimeter. As specified by `RE01CalorimeterParametrisation`, the number of calorimeter layers are (300 - 50) / (3 + 2) = 50. The sizes (rmin, rmax, length) are (50, 53, 200), (55, 58, 200), (60, 63, 200) ...

In sensitive detector construction `void RE01DetectorConstruction::ConstructSDandField`, it initializes `RE01TrackerSD * trackerSD = new RE01TrackerSD(trackerSDname);` and assign it to `fTrackerLayer_log`. Also, it initializes a field `RE01Filed`. It is mentioned that the Calorimeter SD is defined in the parallel world.

The initializer of `RE01TrackerSD` does one thing: to add the name `"trackerCollection"` to the member `collectionName`. Its member function `void RE01TrackerSD::Initialize(G4HCofThisEvent* HCE)` initializes a `fTrackerCollection` and assign it to the `HCE` of type `G4HCofThisEvent*`. The `Initialize` member function is a virtual function of the base class, which is invoked at the beginning of each event. In this case, it creates a hits collection and set it to the `G4HCofThisEvent` obejct.

`fTrackerCollection` is of type `RE01TrackerHitsCollection`, which is defined as `G4THitsCollection<RE01TrackerHit>`. `RE01TrackerHit` is derived from the base class `G4VHit`. It has three member values: `G4double fEdep;`, `G4ThreeVector fPos;`, `G4int fTrackID;`. Its initializer sets these values. The member function `void RE01TrackerHit::Draw()` method creates a circle, and draw it in the `G4VVisManager`. The member function `const std::map<G4String,G4AttDef>* RE01TrackerHit::GetAttDefs() const` returns a description of the quantities "HitType", "TrackID", "Energy", "ETrack", and "Pos". The member function `std::vector<G4AttValue>* RE01TrackerHit::CreateAttValues() const` returns a vector of values of these quantities. The member function `void RE01TrackerHit::Print()` just prints out the member values. The member function `G4bool RE01TrackerSD::ProcessHits(G4Step* aStep, G4TouchableHistory*)` extracts the `edep`, position, user information, trackID from the `G4Step* aStep`, and insert the `newHit` into the `fTrackerCollection`. The `EndOfEvent` method does nothing.

The class `RE01Filed` definition is simple. The initializer sets the value of `fBz`, `fRmax_sq`, and `fZmax`. The member function `void RE01Field::GetFieldValue(const double point[3],double *bfield) const` Sets the field inside a cylinder to be `fBz`, and 0 elsewhere.

# Calorimeter construction
The initializer initializes the base class `G4VUserParallelWorld` with parameter `G4String& parallelWorldName`. And initializes the geometry parameters by including `"RE01DetectorParameterDef.icc"`. Its member function `void RE01CalorimeterROGeometry::Construct()` creates a readout geometry of calorimeter readout geometry. Specifically, it divides the tub into replicas along `kPhi` and `kZAxis`.

In the member function `void RE01CalorimeterROGeometry::ConstructSD()`, it creates a class `RE01CalorimeterSD`, which is derived from the base class `G4VSensitiveDetector`. It has private member values: `RE01CalorimeterHitsCollection *fCalCollection`, `const int fNumberOfCellsInZ;`, and `const int fNumberOfCellsInPhi;`. The class `RE01CalorimeterHitsCollection` is defined as `G4THitsCollection<RE01CalorimeterHit>`. `RE01CalorimeterHit` is derived from `G4VHit`. It has several member values `fZCellID`, `fPhiCellID`, `fEdep`, `fPos`, `fRot`, `const G4LogicalVolume* fPLogV`, `fEdepByATrack`, and `TrackInfo`. 

`RE01CalorimeterSD` overrides the base virtual methods. `void RE01CalorimeterSD::Initialize(G4HCofThisEvent* HCE)` initializes all values of `fCellID` to -1, initializes `fCalCollection = new RE01CalorimeterHitsCollection(SensitiveDetectorName, collectionName[0]);`, and assign it to `HCE`.

`G4bool RE01CalorimeterSD::ProcessHits(G4Step*aStep, G4TouchableHistory*)` Does something more complex. Whenever a step occurs, extracts its energy deposition `G4double edep`. If it equals to 0, then return. Then it gets the z and y replica id, `copyIDinZ` and `copyIDinPhi`, respectively. if `fCellID[copyIDinZ][copyIDinPhi] == -1`, then we insert a hit to `fCalCollection`, and register the index of that hit in `fCalCollection`. Otherwise, we get the id in `fCalCollection`, get the hit, and add the energy deposition in the hit, and `SetTrackInformation(aStep->GetTrack());`

In summary, the geometry part defines a solid world and a parallel world. For each world, it defines a sensitive detector. In solid world, the sensitive detector is of type `RE01TrackerSD`. the `initialize` method creates a `HitsCollection` for every event. And the `ProcessHits` method extracts `edep`, `fPos`, `fTrackID`, and assign them to `newHit`. In parallel world, the sensitive detector is of type `RE01CalorimeterSD`. Its `Initialize` method initializes `fCalCollection` and `fCellID`. Its `ProcessHits` method extracts the information `edep`, translation, rotation, track information etc.

# Action construction
## `RE01PrimaryGeneratorAction`
The initializer initializes its member values, including two `G4VPrimaryGenerator`, one `fMessenger`, and `fUseHEPEvt`. The `GeneratePrimaries` method conditions on the member variable `fUseHEPEvt`. If it is true, then `fHEPEvt->GeneratePrimaryVertex(anEvent);`, else `fParticleGun->GeneratePrimaryVertex(anEvent);`

Its `RunAction` does nearly nothing.

Its `RE01EventAction` has several member values, `G4int fTrackerCollID` and `G4int fCalorimeterCollID`. Its initializer sets the two values to -1. Its `BeginOfEventAction` method gets the collection ID with the collection name, and assign them to the two member values. Its `EndOfEventAction` prints a lot of information, which might be helpful in further understanding.

Something confusing is about the use of the muted. For example, it uses `G4AutoLock lock(&RE01PrimGenDestrMutex);` and `G4AutoLock lock(&RE01PrimGenMutex);`. For now, We don't know what they mean.

The `G4UserStackingAction` class has two member values: `G4int fStage` adn `G4int fCalorimeterHitsColID`. The initializer sets the two values to 0 and -1, respectively. Its method `ClassifyNewTrack(const G4Track * aTrack)` determines the urgent level. The default value is `fUrgent`. If `fStage != 0`, then just return the default value. If `fStage==0`, then if `aTrack->GetTrackStatus()==fSuspend`, then `classification = fWaiting`, and set `trackInfo`. Else, if `aTrack` is the primary particle, then only set the `trackInfo`. The `NewStage()` member function gets the calorimeter hits collection `RE01CalorimeterHitsCollection* CHC`, and prints out the energy depositions of every cells, and the accumulated energy. Personally, I think `RE01StackingAction` is hard to understand.

`RE01TrackingAction` class has two member methods: `PreUserTrackingAction` and `PostUserTrackingAction`. Thr former creates a `RE01Trajectory` for each `G4Track* aTrack` under some conditions. The latter creates a `new RE01TrackInformation(info)` and assign it to every secondaries.

The class `RE01SteppingAction` has primarily one member function `UserSteppingAction(const G4Step* theStep`. It extracts its `thePreRInfo` and `thePostRInfo`. Depending on their values, do something.

# Detailed understanding
## Sensitive detector
The initializer has one parameter: `name`. In the solid class of the initializer, the user has to pass the `name` to the initializer of the base class. Also, `collectionName` must be filled with the name(s) of hits collections created by this particular sensitive detector.

The `Initialize` and `EndOfEvent` methods. In this example, the hits collection class `RE01TrackerHitsCollection * fTrackerCollection` is initialized. The initializer takes two arguments: `SensitiveDetectorName` and `collectionName[0]`, where the former is initialized in the sensitive detector initializer. The initialized `fTrackerCollection` is then set to `G4HCofThisEvent` object, which is the argument for both `Initialize` and `EndOfEvent` methods.`G4HCofThisEvent` object is essentially a vector of `G4VHitsCollection`.

In `ProcessHits` method, a `newHit` object is initialized and added to `fTrackerCollection` by `fTrackerCollection->insert(newHits)`; In this example, `RE01TrackerHit` is derived from `G4Hit`. It has several member values: `fEdep`, `fPos`, `fTrackID`. its base class `G4VHit` is trivial. The user has the option define two member functions `GetAttDefs` and `CreateAttValues`.

## Run action
In this example, the solid class `RE01RunAction` is trivial. However, in the comments of `G4Run`, The user can override the methods `RecordEvent`, `Merge`, which are responsible for recording events and merging local run into global run, respectively.

## Event action
The initializer is trivial. In `BeginOfEventAction` method gets the collection ID of the two collections and assign to its member values `fTrackerCollID` and `fCalorimeterCollID`.

`G4Event` has several member values: `G4HCofThisEvent* HC = nullptr;`, `G4DCofThisEvent* D = nullpr;`, `G4TrajectoryContainer* trajectoryContainer = nullptr`. The user can initialize these containers.

The `EndOfEventAction` prints out all hits, trajectory, and primary vertex.