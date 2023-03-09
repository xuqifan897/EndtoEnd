This note is to help me organize my understanding of the code.

* User interface initialization.
    ```
    ui = new G4UIExecutive(argc, argv);
    ```
    In the header of `G4UIExecutive`, there is an enumeration type
    ```
    enum SessionType { kNone, kQt, kXm, kWin32, kTcsh, kCsh };
    ```
    and a member variable `selected`. At the beginning, the program outputs the available UIs. By default, `selected` is initialized to `kNone`. If the parameter `type` is specified, then the `selected` is initialized according to `type`. Else, if `selected` is still `kNone`, then it is initialized according to the environmental variable. Else if `selected` is still `kNone`, then get the application name, and determine `selected` from the application name as well as a file `${HOME}/.g4session`.

* create run mannager
    ```
    auto* runManager = G4RunManagerFactory::CreateRunManager(G4RunManagerType::Default);
    ```
    There is a enum class `enum class G4RunManagerType: G4int` taking values `{Serial, SerialOnly, MT, MTOnly, Tasking, TaskingOnly, TBB, TBBOnly, Default}`. If the `_type` parameter is `...Only`, then `tail_if_unavail = true;`.

    Else, to determin `rm_type`:
    ```
    rm_type = G4GetEnv<std::string>("G4RUN_MANAGER_TYPE", GetName(_type), "Overriding G4RunManager type...");
    ```
    If the environment variable `"G4RUN_MANAGER_TYPE"` is defined and has value, then it is returned. Otherwise the default value `GetName(_type)` is returned.
    ```
    auto force_rm = G4GetEnv<std::string>("G4FORCE_RUN_MANAGER_TYPE", "", "Forcing G4RunManager type...");
    ```
    get another run manager type `force_rm` with another environment variable `"G4FORCE_RUN_MANAGER_TYPE"`. If the `force_rm` has nontrivial value, then `rm_type` is set to `force_rm`, and `trail_if_unavail` is set to true. Then if both `force_rm` and `rm_type` are trivial, then `rm_type` is set to default, `Tasking`.

    Then, at this point, `rm_type` should have a defined value. Here we get a variable `opts = GetOptions();`. It is a lambda function, which is used to check wheter the `rm_type` is supported.

    `_type = GetYpe(rm_type);` converts the `string` variable back to `enum` and assign to `_type`. The pointer to the run manager is defined: `G4RunManager* rm = nullptr;`. Then the value is initialized according to the value of `_type`. If at this point `rm` is still `nullptr`, then report an error.

    Finally, it tries to cast the `rm` pointer to `G4MTRunManager*` type and assign to mtrm. Static variables are `master_run_manager`, `mt_master_run_manager`, `master_run_manager_kernel` are set according to `rm` and `mtrm`.

    During debugging, we know that `_type == G4RunManagerType::Tasking;`

* Detector construction
    ```
    runManager->SetUserInitialization(new B2a::DetectorConstruction());
    ```
    We firstly look into the function `B2a::DetectorConstruction();` It is rather simple, three lines in total:
    ```
    fMessenger = new DetectorMessenger(this);
    fNbOfChambers = 5;
    fLogicChamber = new G4LogicalVolume*[fNbOfChambers];
    ```
    Then we take a look at the first line. The `DetectorMessenger` class is a user defined class derived from `G4UImessenger`. Its initializer takes a pointer `DetectorConstruction* det`, and assign it to its member value `DetectorConstruction* fDetectorConstruction`.

    `DetectorMessenger` has several member variables
    ```
    DetectorConstruction* fDetectorConstruction
    G4UIdirectory* fDirectory = nullptr;
    G4UIdirectory* fDetDirectory = nullptr;
    G4UIcmdWithAString* fTargMatCmd = nullptr;
    G4UIcmdWithAString* fChamMatCmd = nullptr;
    G4UIcmdWithADoubleAndUnit* fStepMaxCmd = nullptr;
    ```
    In the `DetectorMessenger` initialization, `fDirectory` and `fDetDirectory` are initialized with similar ways, 1) set the `CommandPath` and 2) set the `commandsToBeBroadcasted`. `fTargMatCmd`, `fChamMatCmd`, and `fStepMaxCmd` seems to be specific commands. They are firstly initialized by setting the `commandPath` and the pointer to the messenger. Then set guidance, `fTargMatCmd->SetGuidance("...")`. Then SetParameterName `fTargMatCmd->SetParameterName("choice", false);`. Then specify the states on which the command can be executed `fTargmatCmd->AvailableForStates(G4State_PreInit, G4State_Idle)`.

    The `SetUserInitialization` function is simple, it assigns `userInit` to `userDetector`.

* Physics list initialization
    It only entails three lines of code:
    ```
    G4VModularPhysicsList* physicsList = new FTFP_BERT;
    physicsList->RegisterPhyscs(new G4StepLimiterPhysics());
    runManager->SetUserInitialization(physicsList);
    ```
    In the first line, `FTFP_BERT` is a physics list provided by Geant4. The second line registers a physics to the physics list. The third line associates the physics list to the run manager, which includes two lines:
    ```
    physicsList = userInit;
    kernel->SetPhysics(userInit);
    ```

* Action initialization
    It involves only one line:
    ```
    runManager->SetUserInitialization(new B2::ActionInitialization());
    ```
    The ActionInitialization class is derived from `G4VUserActionInitialization`. The base class has three virtual member functions to be override by the derived class: 
    ```
    virtual void Build() const = 0;
    virtual void BuildForMaster() const;
    virtual G4VSteppingVerbose* InitializeSteppingVerbose() const;
    ```
    The first function is to instantiate user action class. The second function is to instantiate user run action class to be used by G4MTRunManager. The user should not use this method to instantiate user action class except for user run action. I cannot understand the third function at this time.

    In the user defined implementation, only the first two lines are override. the `Build()` function has three lines:
    ```
    SetUserAction(new PrimaryGeneratorAction);
    SetUserAction(new RunAction);
    SetUserAction(new EventAction);
    ```
    `PrimaryGeneratorAction` is derived from `G4VUserPrimaryGeneratorAction`. Besides the initializer and destructor, it has three member functions:
    ```
    void GeneratePrimaries(G4Event*) override;
    G4ParticleGun* GetParticleGun() {return fParticleGun;}
    void SetRandomFlag(G4bool);  # which is not defined

    G4ParticleGun* fParticleGun = nullptr;
    ```
    The `PrimaryGeneratorAction` initializer does the following things: 1) Initialize the particle gun: `fParticleGun = new G4Particle(nofParticles);`, where `nofParticles` defines the number of particles per each invocation. 2) Define the particle, and associate the particle to the fParticleGun. 3) Set the direction and the momentum. As for the position of the particle gun, maybe we'll see it later.

    the `GeneratePrimaries` function, as described in the annotation, is called at the beginning of the event. It sets the position of the particle gun, and invoke the function `fParticleGun->generatePrimaryVertex(anEvent);`

    The `RunAction` class is derived from `G4UserRunAction`. Besides its initializer and destructor, it has two other member functions: `BeginOfRunAction` and `EndOfRunAction`. The initializer does one thing: set printing event number per each 100 events (which I don't totally understand), `G4RunManager::GetRunManager()->SetPrintProgress(1000);`. In this example, the `BeginOfRunAction` does one thing: to inform the runManager to save random number seed: `G4RunManager:GetRunManager()->SetRandomNumberStore(false);`. The function `EndOfRunAction` does nothing. Perhaps other implementation has different customizations.

    The base class of `RunAction`, `G4UserRunAction` has the following methods and members besides what are mentioned above:
    ```
    virtual G4Run* GenerateRun();
    inline virtual void SetMaster(G4bool val = true) { isMaster = val; }
    G4bool isMaster = true;
    ```
    The default return value of `GenerateRun` is a `nullptr`

    The `EventAction` is derived from `G4UserEventAction`. Besides its initializer and destructor, it has two other member functions: `BeginOfEventAction` and `EndOfEventAction`. In this example, `BeginOfEventAction` is empty. `EndOfEventAction` does the following: print the number of trajectories and the number of hits for certain eventID.

    The base class of `EventAction`, `G4UserEventAction` has the following methods and members besides what are mentioned above: `virtual void SetEventManager(G4EventManager* value);` and `G4EventManager* fpEventManager = nullptr;`

* Initialize visualization
    ```
    G4VisManager* visManager = new G4VisExecutive;
    ```
    `G4VisExecutive` is derived from `G4VisManager`. It has several member functions: 1) initializer `G4visExecutive (const G4String& verbosityString = "warnings")`. 2) `void RegisterGraphicsSystems();`, and 3) `void RegisterModelFactories();`. Its initializer only initializes its base class. The two register functions does a lot of registrations, which might not be of interest.