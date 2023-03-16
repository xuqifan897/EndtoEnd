This is the note studying the RE02 example, which does dose calculation.

* Detector construction.
    The detector construction is invoked in the `main` function: `RE02DetectorConstruction* detector = new RE02DetectorConstruction;`. The class `RE02DetectorConstruction` is derived from `G4VUserDetectorConstruction`. Besides its constructor and destructor, it has the following member functions:
    ```
    virtual G4VPhysicalVolume* Construct();  // override the base method
    virtual ConstructSDandField();  // override the base method
    G4ThreeVector fPhantomSize;  // Size of Water Phantom
    G4int fNx, fNy, fNz  // Number of segmentation of water phantom
    G4bool fInsertLead;  // Flag for inserting lead plate in water phantom
    G4LogicalVolume* fLVPhantomSens;
    ```
    The constructor initializes the `fPhantomSize` to be (200mm, 200mm, 400mm), and fNx = fNy = fNz = 100. Then in the main function, the `fPhantomSize` is reset to (200mm, 200mm, 400mm), and fNx, fNy, fNz to 100, 100, 200, respectively. `fInsertLead` is set to `true`. Subsequently, the `detector` is registered to `runManager`: `runManager->SetUserInitialization(detector);`. Then a physics list is registered to `runManager`: `runManager->SetUserInitialization(new QGS_BIG());`. Then the visual manager is initialized. Then action initialization is initialized and registered to runManager: `runManager->SetUserInitialization(new RE02ActionInitialization);`.

    * `userActionInitialization->Build();`
        * `SetUserAction(new RE02PrimaryGeneratorAction)`
        This function initializes several member functions. It sets the momentum direction, particle position, and energy of the particle gun `fParticleGun`. It also sets the member value `fSigmaPosition`. The member function `virtual void GeneratePrimaries(G4Event*)` defines the behavior of primary particle generation.
        * `SetUserAction(new RE02RunAction)`
        This function sets the members `fNx`, `fNy`, `fNz` to 0, and push the name of `MultiFunctionalDetector` to the `fSDName`. Besides the initializer and destructor, it has several member functions. `virtual G4Run* GenerateRun()` returns a `RE02Run` object, which is derived from `G4Run` base object.

        Besides its constructor and destructor, it has several other member variables: `std::vector<G4String> fCollName;`, `std::vector<G4int> fCollID;`, `std::vector<G4HitsMap<G4double>*> fRunMap`. In its initializer, it iterates over all the sensitive detectors and subsequently scorers of the detector, and pushes back the `fullCollectionName`, `collectionID`, `new G4HitsMap<G4double>(detName, collectionName)` to the members `fCollName`, `fCollID`, and `fRunMap`, respectively.

        `virtual void RE02Run::RecordEvent(const G4Event*);` overrides the base class. Essentially, it sums up HitsMap of this event into HitsMap of this run.

        `virtual void Merge(const G4Run*);` merges the local `G4HitsMap` to the global `G4HitsMap`.

        `void DumpAllScorers();` prints the number of entries for each sensitive detector and hits map.

        Besides, the `RE02RunAction` has other member functions. `void RE02RunAction::BeginOfRunAction(const G4Run* aRun)` just prints that a certain run starts. `void RE02RunAction::EndOfRunAction(const G4Run* aRun)` is only implemented by the master thread. The it revises the member values `fNx`, `fNy`, `fNz` according to the detector configurations. And then outputs the value and write to file. At this point, we know that the `HitsMap` has a `[]` operator, which just acts like a n-dimensional matrix.

        * `SetUserAction(new RE02EventAction)`
        Trivial, just print.

    After the above steps, `runManager->Initialize();` is invoked, which performs detector construction, creates physics process, and setup the run.

* `runManager->Initialize()`
    Firstly, it checks the `currentState`, and assures that it is either `G4State_PreInit` or `G4State_Idle`. Then it consecutively initializes geometry and physics.

    * `InitializeGeometry()`
        In the geometry initialization, `G4RunManager::InitializeGeometry()`, firstly, it ensures that the detector is registered. Then the subsequent functions are invoked:
        ```
        userDetector->Construct();
        userDetector->ConstructSDandField();
        userDetector->ConstructParallelGeometris();
        userDetector->ConstructParallelSD();
        ```
        In the `userDetector->Construct();` function, the world physical volume and water phantom physical volume were created. The water physical volume is registered using the logical world volume. Then phantom replication is done by firstly setting a y direction replica `logYRep` object then x direction replica `logXRep`. Finally the logical volume `fLVPhantomSens` was registered to `logXRep`. Visualization attributes were set and associated to different structures. The physical world volume `physiWorld` is returned, and the pointer to logical volume named "phantomSens" is assigned to the member `fLVPhantomSens`.

        In the `userDetector->ConstructSDandField();` function, the multifunctional detector `mFDet` is defined and assigned to the sensitive detector manager `pSDman` and the logical volume `fLVPhantomSens`:
        ```
        G4MultiFunctionalDetector* mFDet = new G4MultiFunctionalDetector(phantomSDname);
        pSDman->AddNewDetector(mFDet);
        fLVPhantomSens->SetSensitiveDetector(mFDet);
        ```
        Proton filter `protonFilter` and electron filter `electronFilter` are defined, as well as a charged particle filter `chargedFilter`. Subsequently, many scores are defined. `scorer0` and `scorer1` are of type `G4PSEnergyDeposit3D`. `protonFilter` is registered to `scorer1`. `scorer2` is a step scorer of type `G4PSNofStep3D`, to which `protonNStep` is registered. `scorer3`, `scorer4`, and `scorer5` are passage cell flux, cell flux, flat surface flux respectively. All scorers are registered to `mFDet`. Then it creates and registers other primitive scorers with different energy bins and particle type gamma.

        `ConstructParallelGeometries()` and `ConstructParallelSD()` are not defined in the solid class.

    * `InitializePhysics()`
        As above, it firstly checks the current state. If the `currentState` is either `G4State_PreInit` or `G4State_Idle`, it is set to the new state `G4State_Init`. For now, we do not plan to dive too deep into the physics part.

