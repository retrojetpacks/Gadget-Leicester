## /rds/user/dc-naya1/rds-dirac-dp005/
#######################################################################
#  Look at end of file for a brief guide to the compile-time options. #
#######################################################################

#--------------------------------------- Basic operation mode of code
#OPT	+=  -DPERIODIC
OPT	+=  -DCOOLING
OPT	+=  -DSFR
OPT	+=  -DUNEQUALSOFTENINGS

# ------ Virtual particles choices ----
OPT += -DVIRTUAL           #Always needed for virtual particles = photons or wind particles. 

# ---- momentum (wind) feedback options
#OPT += -DFB_MOMENTUM		#this models momentum transfer only by wind particles. Can be turned off or on.

# To run full rad transfer, enable VIRTUAL_HEATING, VIRTUAL_RE_EMISSION, LIMITED_ACCRETION_AND_FEEDBACK

#-- radiation transfer options
#OPT += -DPHOTONS 		#turns on RADIATION TRANSFER VIA Monte Carlo packets
#OPT += -DVIRTUAL_RE_EMISSION	#turns on re-emission of photons, otherwise they are just absorbed
#OPT += -DOTHIN_ACCELERATOR #TURNS OFF EMISSION OF PHOTONS, ASSUMING OPTICALLY THIN CASE
#OPT += -DPHOTON_ENERGY_LEDD   #CHOOSE ENERGY OF ONE PHOTON
#OPT += -DBACKGROUND_ILLUMINATION_TEMP #a constant background radiation field of a given T_blackbody



# ------ cooling options
#OPT += -DISOTHERM               #isothermal equation of state (ignored if VIRTUAL_HEATING enabled)
#OPT += -DPOLYTROP               #polytropic equation of state(ignored if VIRTUAL_HEATING enabled)
#OPT += -DADIABATIC         #adiabatic
#OPT += -DCONSTANT_TCOOL  #only works with ADIABATIC and does what it says.
#OPT += -DPOLYTROP_REDUCE               #change K const with time
#OPT += -DCOMMON_COOLING               #COOLING process by Galvagni et al (2012), SHC
#OPT += -DOPACITY_MODIFCATION	     #Modifies dust opacity as a function of rho/rho_tidal. Needs COMMON_COOLING
#OPT += -DH2_DISSOCIATIVE_COLLAPSE     #set maximum T to 2,000 K to mimick Hs dissociative collapse
#OPT += -DSECOND_CORE     #mimicks finite H2 dissociative energy
#OPT += -DREDUCE_PLANET_OPACITY	#reduces optical depth at very high densities if common_cooling used
#OPT += -DVIRTUAL_HEATING 	#heating of SPH particles is calculated via fb_particles.c
#OPT += -DBETA_COOLING   #Beta-cooling model like Rice et al 2005.
#OPT += -DEVAPORATION #Heats gas below Evap_dens parameter, evaporating gas clouds
OPT += -DEVAPORATION_RADIAL #Radial dependence for evaporation
#OPT += -DBETA_VARY      #Beta varies with time -- turn off BETA_DIP_AND_GROW!
#OPT += -DBETA_DIP_AND_GROW      #Beta first decreases and then increases with time -- turn off BETA_VARY!
#OPT += -DBETA_COOLING_TAPPER_OFF #turn off cooling at very high densities
#OPT += -DDO_NOT_PROTECT    #decrease protection level from a sudden entropy change to conserve energy better
#OPT += -DSMART_HEATING    #temporary to test heating/cooling within smart_diffusion
#OPT += -DSUPPRESS_COOLING_TIDAL  #increase tcool if rho>rho_tid--requires CUSP_POTENTIAL,BETA_COOLING
#OPT += -DQUASAR_HEATING		  #Sazonov et al Quasar heating rate
#OPT += -DDUST_OPACITY_FIT        #Uses k(T) = k0 (T/T0) opacity law
OPT += -DFIND_SMBH 	  #Find location of the black hole particle, may be needed for central irradiation, etc.
#OPT += -DSMBH_IRRADIATION #FindS minimum temperature of the gas as a function of distance from SMBH
#OPT += -DPLANET_IRRADIATION #FindS minimum temperature of the gas as a function of distance from a planet

#--- feedback options
#OPT += -DFRACTION_OF_LEDD_FEEDBACK #emits feedback as VirtualFeedBack x L_Edd for the given P.Mass
#OPT += -DSTARTUP_LUMINOSITY   #neglect emission of photons if new_virt << 1 (useful if too many low L sinks appear)
#OPT += -DFRACTION_OF_LSOLAR_FB

#---------- increase SPH particle mass with time for external disc feeding tests
#OPT += -DINCREASE_SPH_MASS #needs parameter DeltaTimeGrow in the input file
#OPT += -DINCREASE_SPH_MASS_LIMIT #only allow SPH mass to increase by 20% and then keep it constant


#-- gas accretion options. SWALLOW_GAS needs to be enabled for accretion
OPT += -DACCRETION_RADIUS       #accretion within the accretion radius
OPT += -DACCRETION_DENSITY       #REQUIRES GAS TO BE DENSER THAN A MINIMUM FOR ACCRETION
OPT += -DBH_MASS_REAL           #uses real BH mass in calculating accretion rate (P.Mass not P.BH_Mass)
OPT += -DLIMITED_ACCRETION_AND_FEEDBACK #a reservoir of mass stores super-edd accreted material to smooth out feedback
OPT += -DTMP_FEEDBACK          #SN: thermal feedback BlackHoleFeedbackFactor * Mdot, where BHFF is erg/g released
OPT += -DACCRETION_OF_DUST_ONLY #SN: only disk particles are allowed to be accreted
#OPT += -DPLANET_ACCRETION_FEEDBACK #SN: thermal feedback G M_p Mdot/R_p, where R_p = 2R_Jupiter
#OPT += -DPLANET_ACCRETION_REPORTING #SN: prints in the log file only the planet's mass and accretion rate
#OPT += -DPLANET_INITIAL_MASS	#SN: redefines the smaller BH/planet initial mass using the input file
#OPT += -DSTAR_FROMATION_FEEDBACK #uses stellar accr luminosity for mdot_edd
#OPT += -DNO_ACCRETION           #turn off accretion in blackhole.c
#OPT += -DACCRETION_AT_EDDINGTON_RATE #fixes mdot at mdot_eddington

#--external grav potential: NOTE: only one or none to be used !@*
#OPT += -DSGRA_POTENTIAL         #ADD STELLAR CUSP AS EXTERNAL POTETNIAL 
#OPT += -DCUSP_POTENTIAL         #ADD STELLAR CUSP AS EXTERNAL POTETNIAL 
#OPT += -DSIS_POTENTIAL		 #Singular Isothermal Sphere pot.
#OPT += -DNFW_POTENTIAL          # Add background NFW potential

#--dynamics options
#OPT += -DBH_FIXED               #zeros all kicks that the SMBH receives to have it fixed in space in 0,0,0
#OPT += -DBH_KICK		 #gives a kick to BH initially
#OPT += -DGAS_FIXED              #Fix the gas particles' position for tests

#--star formation options, mergers 
OPT += -DBH_FORM    	        # Turns sph particles into BHs, disables default Gadget star formation
#OPT += -DSMBH_MASS_INIT         #sets P[i].BH_Mass P[i].Mass at t=0. requires BLACK_HOLES
#OPT += -DJEANS_MASS_SF	        #introduce sink particles in accord with Jeans criterium
#OPT += -DNO_BH_MERGERS 	#turns off mergers of BHs in blackhole.c
OPT += -DBH_MERGERS_WITHIN_H	#enforces mergers when BHs are separated by less than softening H
#OPT += -DTIDAL_JEANS_SF	#for SF near SMBH; tidal density for sf, cooling -- requires CUSP_POTENTIAL


#--photon/virt particle propagation options 
OPT +=-DVARIABLE_PHOTON_SEARCH_RADIUS # CBP (06/11/2008)
OPT += -DIGNORE_NON_UNIQUE 	#ignore non-unique IDs in the initial conditions file
#OPT += -DEXCLUDE_VIRTUAL_FROM_TREE
#OPT += -DNO_VIRTUAL_GRAVITY     #turns off gravity acting on virtual particles


OPT += -DVIRTUAL_SET_MASS	#sets OriginalGasMass to an sph particle mass, assumed all equal.
OPT += -DVIRTUAL_FLY_THORUGH_EMPTY_SPACE #allow for large delta t or delta r of photons in empty space.
#OPT += -DNO_BH_ACCRETION 	#TURNS OFF BLACKHOLE ACCRETION; USEFUL FOR SOME TESTING

#OPT += -DRAD_ACCEL		# Radiation acceleration calculation switch
#OPT += -DVIRTUAL_OTHIN 	# calculates the momentum passed to SPH partcles, but doesn't reduce photon momentum
#OPT += -DISOTROPIC #if set, then emit particles isotropically in fb_particles.c

#OPT += -DIGNORE_NFU 		#temporary to understand crashes in predict.c
#OPT += -DDUMP_BY_NSTEP
#OPT += -DVIRTUAL_LARGE_DTIME # rough timestep for photons in timestep.c
OPT += -DCONSTANT_MEAN_MOLECULAR_WEIGHT # using a constat value as the mu
#OPT += -DCHECK_ENERGY_CONSERVATION # Checking the energy conservation

#OPT += -DNO_DAUGHTER_PHOTONS	#traditional monte-carlo where photons are absorbed completely and then re-emitted
#OPT += -DSPH_EMIT		#allow SPH particles to emit photons
#OPT += -DREAL_EOS              #For more realistic EOS by SHC
OPT += -DDUST                   #For work with dust grain dynamics by SHC and SN
OPT += -DDUST_TIMESTEP
#OPT += -DDUST_GROWTH           #allows grain growth
#OPT += -DDUST_GROWTH_FIXED_SIZE           #treats grain growth via change in pebble mass rather than size (overtakes DUST_GROWTH)
#OPT += -DDUST_TWO_POPULATIONS   #dust is divided on microscopic grains (moves with gas) + pebbles (independent)
#OPT += -DDUST_PEBBLES_BORN   #pebbles are born on the fly. Requires dust_two_populations
#OPT += -DDUST_EPSTEIN		#friction in the Epstein regime only, for tests
#OPT += -DDUST_VAPORIZE         #allows grain vaporization
#OPT += -DDUST_FE_AND_ICE_GRAINS      #even ID grains are FE, uneven water ICE. Requires dust_vaporize!
#OPT += -DDUST_ENERGY_CONSERVATION #enforces en conservation by injecting "missing" energy in gas in sfr_eff.c. 
#OPT += -DDUST_2ND_POPULATION    #pop. of small particles tightly coupled to gas, requires DUST_GROWTH
#OPT += -DDUST_REAL_PEBBLE_COLLISIONS #calculate velocity dispersion of pebbles
#OPT += -DDUST_SINK_ON_FLY #introduce a dust sink particle by turning a dust particle into the sink
#OPT += -DDUST_NO_FRICTION_HEATING #neglect heating of gas due to friction -- useful for tests
#OPT += -DDUST_MDUST_GROW  #dust particle mass grows with time; uses VirtualStart in the input file
OPT += -DDUST_POWERLAW  #RJH set initial dust population size distribution after injection


#OPT += -DDECREASE_POISSON_NOISE  #averages rad energy dens in time to decrease fluctuations
#OPT += -DRADIATIVE_EQUILIBRIUM   #assumes that du/dt term is small in energy equilibrium

#OPT += -DNON_GRAY_RAD_TRANSFER  #every photon carries its own opacity with it


#--------------------------------------- TreePM Options
#OPT	+=  -DPMGRID=256
#OPT	+=  -DGRIDBOOST=2
#OPT	+=  -DASMTH=1.25
#OPT	+=  -DRCUT=4.5
#OPT	+=  -DPLACEHIGHRESREGION=3
#OPT	+=  -DENLARGEREGION=1.2


#--------------------------------------- Multi-Domain and Top-Level Tree options
OPT	+=  -DMULTIPLEDOMAINS=16
#OPT	+=  -DTOPNODEFACTOR=100.0


#--------------------------------------- Things that are always recommended
OPT	+=  -DPEANOHILBERT
OPT	+=  -DWALLCLOCK
OPT	+=  -DMYSORT
OPT	+=  -DPEDANTIC_MEMORY_HANDLER   # this enables an internal memory handler (most recently allocated block needs to be freed first)
#OPT	+=  -DCPUSPEEDADJUSTMENT


#---------------------------------------- Single/Double Precision
OPT	+=  -DDOUBLEPRECISION
#OPT	+=  -DDOUBLEPRECISION_FFTW
#OPT	+=  -DOUTPUT_IN_DOUBLEPRECISION # snapshot files will be written in double precision
#OPT	+=  -DINPUT_IN_DOUBLEPRECISION  # initial conditions are in double precision

#OPT	+=  -DFLTROUNDOFFREDUCTION      # enables (expensive!) `double-double' round-off reduction in particle sums
#OPT	+=  -DSOFTDOUBLEDOUBLE          # needs to be set if a C++ software implementation of 128bit double-double precision should be used


#---------------------------------------- On the fly FOF groupfinder 
#OPT	+=  -DFOF                                # enable FoF output
#OPT	+=  -DFOF_PRIMARY_LINK_TYPES=2           # 2^type for the primary dark matter type
#OPT	+=  -DFOF_SECONDARY_LINK_TYPES=(1+16+32) # 2^type for the types linked to nearest primaries
#OPT	+=  -DSUBFIND                            # enables substructure finder
#OPT	+=  -DORDER_SNAPSHOTS_BY_ID

#--------------------------------------- SFR/feedback model
#OPT	+=  -DSOFTEREQS
OPT	+=  -DMOREPARAMS
#OPT	+=  -DMETALS
OPT	+=  -DSTELLARAGE            # SN: I use it for virtual particle destruction only
#OPT	+=  -DWINDS
#OPT	+=  -DQUICK_LYALPHA
#OPT	+=  -DISOTROPICWINDS
#OPT	+=  -DMHM


#-------------------------------------- AGN stuff
OPT	+=  -DBLACK_HOLES             # enables Black-Holes (master switch)
#OPT	+=  -DBONDI                   # Bondi-Hoyle style accretion model
#OPT	+=  -DENFORCE_EDDINGTON_LIMIT # put a hard limit on the maximum accretion rate
OPT	+=  -DBH_THERMALFEEDBACK      # couple a fraction of the BH luminosity into surrounding gas
#OPT	+=  -DBH_DRAG                 # Drag on black-holes due to accretion
OPT	+=  -DSWALLOWGAS              # Enables stochastic accretion of gas particles consistent with growth rate of hole
#OPT	+=  -DEVALPOTENTIAL           # computes gravitational potential
#OPT	+=  -DREPOSITION_ON_POTMIN    # repositions hole on potential minimum (requires EVALPOTENTIAL)


#-------------------------------------- AGN-Bubble feedback
#OPT	+=  -DBUBBLES                 # generation of hot bubbles in an isolated halo or the the biggest halo in the run
#OPT	+=  -DMULTI_BUBBLES 	      # hot bubbles in all haloes above certain mass threshold (works only without FOF and BUBBLES)
#OPT	+=  -DEBUB_PROPTO_BHAR        # Energy content of the bubbles with cosmic time evolves as an integrated BHAR(z) over a Salpeter time (Di Matteo 2003 eq. [11])


#-------------------------------------- Viscous gas treatment 
#OPT	+=  -DNAVIERSTOKES            # Braginskii-Spitzer parametrization of the shear viscosity: mu = f x T^{5/2}
#OPT	+=  -DNAVIERSTOKES_CONSTANT   # Shear viscosity set constant for all gas particles
#OPT	+=  -DNAVIERSTOKES_BULK       # Bulk viscosity set constant for all gas particles. To run with bulk visocity only one has to set shear viscosity to zero in the parameterfile.
#OPT	+=  -DVISCOSITY_SATURATION    # Both shear and bulk viscosities are saturated, so that unphysical accelerations and entropy increases are avoided. Relevant for the cosmological simulations.
#OPT	+=  -DNS_TIMESTEP             # Enables timestep criterion based on entropy increase due to internal friction forces
#OPT	+=  -DOUTPUTSTRESS            # Outputs diagonal and offdiagonal components of viscous shear stress tensor
#OPT	+=  -DOUTPUTBULKSTRESS        # Outputs viscous bulk stress tensor
#OPT	+=  -DOUTPUTSHEARCOEFF        # Outputs variable shear viscosity coefficient in internal code units


#-------------------------------------------- Things for special behaviour
#OPT	+=  -DPEDANTIC_MEMORY_CEILING=400   # this should be set to define the memory ceiling in MByte
#OPT	+=  -DDOMAIN_MEMORY_CEILING=400
#OPT    +=  -DNO_ISEND_IRECV_IN_DOMAIN
OPT	+=  -DFIX_PATHSCALE_MPI_STATUS_IGNORE_BUG
#OPT	+=  -DNOGRAVITY
#OPT	+=  -DNOACCEL
#OPT	+=  -DNOISMPRESSURE
#OPT	+=  -DNOVISCOSITYLIMITER
#OPT	+=  -DNOTREERND
#OPT	+=  -DNOSTOP_WHEN_BELOW_MINTIMESTEP
#OPT	+=  -DNOPMSTEPADJUSTMENT
#OPT	+=  -DNOTYPEPREFIX_FFTW
OPT	+=  -DNO_TREEDATA_IN_RESTART
#OPT	+=  -DNOWINDTIMESTEPPING            # Disable wind reducing timestep (not recommended)
#OPT	+=  -DISOTHERM=200                  # adds potential of an isothermal sphere
OPT	+=  -DCOMPUTE_POTENTIAL_ENERGY
#OPT	+=  -DALLOWEXTRAPARAMS
#OPT	+=  -DLONGIDS
OPT	+=  -DINHOMOG_GASDISTR_HINT         # if the gas is distributed very different from collisionless particles, this can helps to avoid problems in the domain decomposition
#OPT	+=  -DLONG_X=100
#OPT	+=  -DLONG_Y=5
#OPT	+=  -DLONG_Z=0.2
#OPT	+=  -DTWODIMS
#OPT	+=  -DSPH_BND_PARTICLES
#OPT	+=  -DNEW_RATES                     # switches in updated cooling rates from Naoki
#OPT	+=  -DRADIATIVE_RATES               # used in non-equilibrium chemistry model
#OPT	+=  -DREAD_HSML                     # reads hsml from IC file
OPT   	+=  -DADAPTIVE_GRAVSOFT_FORGAS      # allows variable softening length for gas particles (requires UNEQUALSOFTENINGLENGTH)
OPT	+=  -DADAPTIVE_GRAVSOFT_FORGAS_HSML # this sets the gravitational softening for SPH particles equal to the SPH smoothing (requires ADAPTIVE_GRAVSOFT_FORGAS)
#OPT	+=  -DGENERATE_GAS_IN_ICS
#OPT    +=  -DNEUTRINOS                     # Option for special integration of light neutrino species 


#--------------------------------------- Time integration options
#OPT	+=  -DALTERNATIVE_VISCOUS_TIMESTEP


#--------------------------------------- Output/Input options
#OPT	+=  -DOUTPUTPOTENTIAL
#OPT	+=  -DRECOMPUTE_POTENTIAL_ON_OUTPUT # update potential every output even it EVALPOTENTIAL is set
#OPT	+=  -DOUTPUTACCELERATION
#OPT	+=  -DOUTPUTCHANGEOFENTROPY
#OPT	+=  -DOUTPUTTIMESTEP
#OPT	+=  -DOUTPUTCOOLRATE                # outputs cooling rate, and conduction rate if enabled
OPT	+=  -DHAVE_HDF5                     # needed when HDF5 I/O support is desired
#OPT	+=  -DOUTPUTBSMOOTH
#OPT	+=  -DOUTPUTDENSNORM
#OPT	+=  -DXXLINFO                       # Enables additional output for viscosityand bfield
#OPT	+=  -DOUTPUTLINEOFSIGHT             # enables on-the-fly output of Ly-alpha absorption spectra
#OPT	+=  -DOUTPUTLINEOFSIGHT_SPECTRUM
#OPT	+=  -DOUTPUTLINEOFSIGHT_PARTICLES
#OPT	+=  -DOUTPUT_TIDALTENSOR            #  outputs tidal tensor (=matrix of second derivatives of grav. potential)


#--------------------------------------- Testing and Debugging options
#OPT	+=  -DFORCETEST=0.1
#OPT	+=  -DDEBUG                     # enables core-dumps and FPU exceptions
#OPT	+=  -DPARTICLE_DEBUG            # auxiliary communication of IDs
#OPT	+=  -DVERBOSE


#--------------------------------------- Static NFW Potential
#OPT	+=  -DSTATICNFW
#OPT	+=  -DNFW_C=12
#OPT	+=  -DNFW_M200=100.0
#OPT	+=  -DNFW_Eps=0.01
#OPT	+=  -DNFW_DARKFRACTION=0.87


#--------------------------------------- Static Hernquist Potential
#OPT	+=  -DSTATICHQ
#OPT	+=  -DHQ_M200=1.0
#OPT	+=  -DHQ_C=10
#OPT	+=  -DHQ_DARKFRACTION=0.9


#--------------------------------------- Thermal conduction
#OPT	+=  -DCONDUCTION
#OPT	+=  -DCONDUCTION_CONSTANT
#OPT	+=  -DCONDUCTION_SATURATION


#--------------------------------------- Dark energy
#OPT	+=  -DDARKENERGY # Enables Dark Energy
#OPT	+=  -DTIMEDEPDE  # read w(z) from a DE file
#OPT	+=  -DRESCALEVINI # rescale v_ini in read_ic / read_ic_cluster
#OPT    +=  -DEXTERNALHUBBLE # reads the hubble function from the DE file
#OPT    +=  -DTIMEDEPGRAV # resacles H and G according to DE model
#OPT    +=  -DDARKENERGY_DEBUG # enable writing of drift/kick table


#--------------------------------------- Long-range scalar field
#OPT	+=  -DSCALARFIELD


#--------------------------------------- SPH viscosity options
#OPT	+=  -DCONVENTIONAL_VISCOSITY     # enables the old viscosity
#OPT	+=  -DTIME_DEP_ART_VISC          # Enables time dependend viscosity
#OPT	+=  -DNO_SHEAR_VISCOSITY_LIMITER # Turns of the shear viscosity supression
#OPT	+=  -DHIGH_ART_VISC_START        # Start with high rather than low viscosity
#OPT	+=  -DALTVISCOSITY               # enables alternative viscosity based on div(v)


#--------------------------------------- Magnetic Field options
#OPT	+=  -DMAGNETIC
#OPT	+=  -DMAGFORCE
#OPT	+=  -DBRIOWU
#OPT	+=  -DARTBPRES
#OPT	+=  -DDIVBFORCE
#OPT	+=  -DTRACEDIVB
#OPT	+=  -DDBOUTPUT
#OPT	+=  -DCORRECTDB
#OPT	+=  -DCORRECTBFRC
#OPT	+=  -DBINISET
#OPT	+=  -DBSMOOTH
#OPT	+=  -DBFROMROTA
#OPT	+=  -DIGNORE_PERIODIC_IN_ROTA
#OPT	+=  -DMAGNETIC_SIGNALVEL
#OPT	+=  -DMAGNETIC_DISSIPATION
#OPT	+=  -DMAGDISSIPATION_PERPEN
#OPT	+=  -DTIME_DEP_MAGN_DISP
#OPT	+=  -DHIGH_MAGN_DISP_START
#OPT	+=  -DDIVBCLEANING_DEDNER
#OPT	+=  -DSMOOTH_PHI


#--------------------------------------- Glass making
#OPT	+=  -DMAKEGLASS


#--------------------------------------- Raditive Transfer
#OPT	+=  -DRADTRANSFER


#--------------------------------------- Distortion Tensor
#OPT	+=  -DDISTORTIONTENSOR


#---------------------------------------- nonequilibrium proimodal chemisitry
#OPT	+=  -DNONEQUILIBRIUM
#OPT	+=  -DCHEMISTRY
#OPT	+=  -DCMB
#OPT	+=  -DRADIATION


#---------------------------------------- Cosmic Rays (Martin)
#OPT	+=  -DCOSMIC_RAYS               # Cosmic Rays Master Switch
#OPT	+=  -DCR_IC                     # IC files contain CR information
#OPT	+=  -DCR_IC_PHYSICAL
#OPT	+=  -DCR_DISSIPATION            # Catastrophic losses
#OPT	+=  -DCR_THERMALIZATION         # Coulomb cooling
#OPT	+=  -DCR_SHOCK=2                # Shock energy is directed into CR
			                # 2 = Mach-Number dependent shocks, Mach-number derived for thermal gas/CR composite
			                # 3 = Mach-Number dependent shocks, Mach-number derived for thermal gas
#OPT	+=  -DCR_DIFFUSION              # Cosmic Ray diffusion
#OPT	+=  -DCR_DIFFUSION_GREEN        # alternative diffusion model
#OPT	+=  -DUPDATE_PARANOIA=1         # 1 = Update on every predict, 2 = Update on every energy injection and on every predict
#OPT	+=  -DCR_OUTPUT_INJECTION       # Output energy injection rate in snapshots
#OPT	+=  -DCR_NO_CHANGE              # Compute changes to CR, but do not execute them, useful for estimating the size of effects
#OPT	+=  -DCOSMIC_RAY_TEST           # starts a test routine instead of the simulation
#OPT	+=  -DCR_NOPRESSURE             # computes CRs as usual, but ignores the pressure in the dynamics


#---------------------------------------- Mach number finder (Christoph)
#OPT	+=  -DMACHNUM                   # Mach number Master Switch
#OPT	+=  -DMACHSTATISTIC             # Dissipated thermal energy at shocks
#OPT	+=  -DCR_OUTPUT_JUMP_CONDITIONS # CR: density and thermal energy jump at shocks




#--------------------------------------- Select target Computer

#SYSTYPE="McKenzie"
#SYSTYPE="CITA"
#SYSTYPE="Stella"
#SYSTYPE="Sauron"
#SYSTYPE="Sauron-gcc"
#SYSTYPE="Mako"
#SYSTYPE="MPA"
#SYSTYPE="Regatta"
#SYSTYPE="Solaris"
#SYSTYPE="RZG_LinuxCluster"
#SYSTYPE="RZG_LinuxCluster-gcc"
#SYSTYPE="CINECA32"
#SYSTYPE="CINECA64"
#SYSTYPE="DEI32"
#SYSTYPE="OpteronMPA"
#SYSTYPE="Ingeld_LinuxCluster"
#SYSTYPE="OPA-Cluster32"
#SYSTYPE="OPA-Cluster64"
#SYSTYPE="Titania"
#SYSTYPE="Quintor"
#SYSTYPE="Zijl"
#SYSTYPE="hpcf"
#SYSTYPE="LeicesterOpteron"






SYSTYPE="ALICE"
SYSTYPE = "DARWIN_NEW_IMPI"


CC       = mpicc        # sets the C-compiler (default)
OPTIMIZE = -Wall  -g   # optimization and warning flags (default)

MPICHLIB = -lmpich

ifeq (SOFTDOUBLEDOUBLE,$(findstring SOFTDOUBLEDOUBLE,$(OPT)))
CC       =   mpiCC     # default C++ compiler
OPTIMIZE =   -g 
OPT     +=  -DX86FIX   # only needed for 32-bit intel/amd systems
endif

ifeq ($(SYSTYPE),"MyLaptop")
CC       =  mpicc

OPTIMIZE =   -O3 -Wall
GSL_INCL =  -I/usr/include
GSL_LIBS =  -L/usr/lib 
FFTW_INCL=  
FFTW_LIBS=
MPICHLIB =
HDF5INCL = -I/opt/hdf5/include
HDF5LIB  = -L/opt/hdf5/lib -static -lhdf5 -lz
endif

ifeq ($(SYSTYPE),"DARWIN")
CC       =  mpicc
OPTIMIZE =   -O3 -Wall -g
GSL_INCL =  -I/usr/local/Cluster-Apps/gsl/1.9/include
GSL_LIBS =  -L/usr/local/Cluster-Apps/gsl/1.9/lib
FFTW_INCL=  -I/usr/local/Cluster-Apps/fftw/intel/64/2.1.5/float/include
FFTW_LIBS=  -L/usr/local/Cluster-Apps/fftw/intel/64/2.1.5/float/lib
MPICHLIB =
#HDF5INCL =  -I/usr/local/Cluster-Apps/hdf5/impi/1.8.0/include
#HDF5LIB  =  -L/usr/local/Cluster-Apps/hdf5/impi/1.8.0/lib -lhdf5
HDF5INCL = -I/sw/lib     -I/opt/local/include -DH5_USE_16_API #-DUSE_SSE   
HDF5LIB  = -L/sw/lib -L/opt/local/lib  -lhdf5 -lz

endif

ifeq ($(SYSTYPE),"DARWIN_NEW_IMPI")
CC       = mpicc 
##OPTIMIZE = -O1 -xHOST -g 
OPTIMIZE = -O3 -Wall -xHOST -g 
HDF5LIB  = -lhdf5 -lz 
HDF5INCL = -DH5_USE_16_API 
MPICHLIB = -lmpi
##GSL_INCL =  -I/usr/local/Cluster-Apps/gsl/1.9/include
##GSL_LIBS =  -L/usr/local/Cluster-Apps/gsl/1.9/lib
EXTRA_LDFLAGS = -shared-intel $(OMPI_LDFLAGS) # only for intel compiler
endif


ifeq ($(SYSTYPE),"ALICE")
CC       =  mpicc
OPTIMIZE =   -O3 -Wall -m64
#GSL_INCL =  -I/home/c/cbp1/include/
#GSL_LIBS =  -L/home/c/cbp1/lib
#GSL_INCL =  -I/home/s/sn85/gsl/include
#GSL_LIBS =  -L/home/s/sn85/gsl/lib
GSL_INCL =  -I/scratch/astro1/shared/gsl/include
GSL_LIBS =  -L/scratch/astro1/shared/gsl/lib
#FFTW_INCL=  -I/home/s/sn85/fftw/include/
#FFTW_LIBS=  -L/home/s/sn85/fftw/lib
FFTW_INCL=
FFTW_LIBS=
MPICHLIB =
#HDF5INCL = -I/opt/hdf5/include
#HDF5LIB  = -L/opt/hdf5/lib -static -lhdf5 -lzendif
HDF5INCL = -I/cm/shared/apps/hdf5/current/include  -DH5_USE_16_API
HDF5LIB = -L/cm/shared/apps/hdf5/current/lib -lhdf5 -lz

endif


ifeq ($(SYSTYPE),"LeicesterOpteron")
#CC       =    /usr/local/mpich-mx_GNU/bin/mpicc   # sets the C-compiler
#CC       =  /opt/openmpi/gnu/bin/mpicc # sets the C-compiler

OPTIMIZE =   -O3 -Wall -m64
GSL_INCL =  -L/usr/local/include
GSL_LIBS =  -L/usr/local/lib
FFTW_INCL=
FFTW_LIBS=
MPICHLIB =
HDF5INCL = -I/opt/hdf5/include
HDF5LIB  = -L/opt/hdf5/lib -static -lhdf5 -lz
endif












ifeq ($(SYSTYPE),"McKenzie")
CC       =  mpicc   # sets the C-compiler
OPTIMIZE =  -g -O3 #-Wall
GSL_INCL = -I/usr/include
GSL_LIBS = -L/usr/lib
FFTW_INCL= -I/opt/fftw/intel_8.1/2.1.5/include
FFTW_LIBS= -L/opt/fftw/intel_8.1/2.1.5/lib #-ldrfftw_mpi
MPICHLIB = -L/opt/lam-7.1.2b24-g77/lib -lmpi
HDF5INCL = -I/opt/hdf5-oscar-1.6.4/include
HDF5LIB  = -L/opt/hdf5-oscar-1.6.4/lib -lhdf5 -lz
endif

ifeq ($(SYSTYPE),"CITA")
CC       =  mpicc
OPTIMIZE =  -O3 -Wall
GSL_INCL =  -I/usr/include/gsl
GSL_LIBS =  -L/usr/lib/libgsl
FFTW_INCL=  -I/opt/fftw-2.1.5/include
FFTW_LIBS=  -L/opt/fftw-2.1.5/lib
MPICHLIB =  -L/usr/lib/libmpi
HDF5INCL =  -I/usr/include
HDF5LIB  =  -L/usr/lib/libhdf5 -static -lhdf5 -lz
endif 


ifeq ($(SYSTYPE),"Titania")
CC       = mpcc
#PTIMIZE = -fast -xtarget=ultra3 -xarch=v9b -xcache=64/32/4:8192/512/1
OPTIMIZE = -g -xtarget=ultra3 -xarch=v9b
GSL_INCL = -I/opt/local/gsl/64/current/include
GSL_LIBS = -L/opt/local/gsl/64/current/lib
FFTW_INCL= -I/opt/local/fftw/64/current/include
FFTW_LIBS= -L/opt/local/fftw/64/current/lib
MPICHLIB = -lmpi
HDF5INCL = -I/opt/local/hdf5/64/current/include
HDF5LIB  = -L/opt/local/hdf5/64/current/lib -lhdf5
RLIBS    = -R/opt/local/gsl/64/current/lib:/opt/local/fftw/64/current/lib:/opt/local/hdf5/64/current/lib
endif


ifeq ($(SYSTYPE),"Quintor")
CC       = mpcc
OPTIMIZE = -fast -xtarget=ultra3i -xarch=v9b -xcache=64/32/4:1024/64/4
#OPTIMIZE = -g -xtarget=ultra3i -xarch=v9b
GSL_INCL = -I/opt/local/gsl/64/current/include
GSL_LIBS = -L/opt/local/gsl/64/current/lib
FFTW_INCL= -I/opt/local/fftw/64/current/include
FFTW_LIBS= -L/opt/local/fftw/64/current/lib
MPICHLIB = -lmpi
HDF5INCL = -I/opt/local/hdf5/64/current/include
HDF5LIB  = -L/opt/local/hdf5/64/current/lib -lhdf5
RLIBS    = -R/opt/local/gsl/64/current/lib:/opt/local/fftw/64/current/lib:/opt/local/hdf5/64/current/lib
endif


ifeq ($(SYSTYPE),"Stella")
CC       =  mpicc
OPTIMIZE =  -O3 -Wall
GSL_INCL =  -I/home/schaye/libs/include
GSL_LIBS =  -L/home/schaye/libs/lib -static
FFTW_INCL=  -I/home/schaye/libs/include
FFTW_LIBS=  -L/home/schaye/libs/lib
MPICHLIB =
HDF5INCL =
HDF5LIB  =
OPT      +=  -DNOCALLSOFSYSTEM
endif


ifeq ($(SYSTYPE),"OPA-Cluster32")
CC       =   mpiccg
ifeq (SOFTDOUBLEDOUBLE,$(findstring SOFTDOUBLEDOUBLE,$(OPT)))
CC       =   mpiCCg
OPT     +=  -DX86FIX
endif
OPTIMIZE =  -O3 -Wall
GSL_INCL =  -I/afs/rzg/bc-b/vrs/opteron32/include
GSL_LIBS =  -L/afs/rzg/bc-b/vrs/opteron32/lib -Xlinker -R -Xlinker /afs/rzg/bc-b/vrs/opteron32/lib
FFTW_INCL=  -I/afs/rzg/bc-b/vrs/opteron32/include
FFTW_LIBS=  -L/afs/rzg/bc-b/vrs/opteron32/lib
MPICHLIB =
HDF5INCL =
HDF5LIB  =
OPT      +=  -DNOCALLSOFSYSTEM
endif

ifeq ($(SYSTYPE),"OPA-Cluster64")
CC       =   mpiccg 
ifeq (SOFTDOUBLEDOUBLE,$(findstring SOFTDOUBLEDOUBLE,$(OPT))) 
CC       =   mpiCCg
endif
OPTIMIZE =  -O3 -g -Wall -m64
GSL_INCL =  -I/afs/rzg/home/v/vrs/Libs/opteron64/mvapich2-0.9.3/include
GSL_LIBS =  -L/afs/rzg/home/v/vrs/Libs/opteron64/mvapich2-0.9.3/lib  -Wl,"-R /afs/rzg/home/v/vrs/Libs/opteron64/mvapich2-0.9.3/lib"
FFTW_INCL=  -I/afs/rzg/home/v/vrs/Libs/opteron64/mvapich2-0.9.3/include
FFTW_LIBS=  -L/afs/rzg/home/v/vrs/Libs/opteron64/mvapich2-0.9.3/lib
MPICHLIB =
HDF5INCL =  -I/afs/rzg/home/v/vrs/Libs/opteron64/mvapich2-0.9.3/include
HDF5LIB  =  -L/afs/rzg/home/v/vrs/Libs/opteron64/mvapich2-0.9.3/lib  -lhdf5 -lz 
OPT      +=  -DNOCALLSOFSYSTEM
endif


ifeq ($(SYSTYPE),"Sauron-gcc")
CC       =   mpicc.gcc   # sets the C-compiler
OPTIMIZE =   -O3 -funroll-loops -march=k8 -msse2 -static
GSL_INCL =   -I/usr/local/gsl.gcc/include
GSL_LIBS =   -L/usr/local/gsl.gcc/lib -static -lgsl -lgslcblas
FFTW_INCL=   -I/usr/local/fftw.gcc/include
FFTW_LIBS=   -L/usr/local/fftw.gcc/lib -static -lrfftw_mpi -lfftw_mpi -lrfftw -lfftw
MPICHLIB =
endif


ifeq ($(SYSTYPE),"Sauron")
CC       =  mpicc  -m64 # sets the C-compiler
ifeq (SOFTDOUBLEDOUBLE,$(findstring SOFTDOUBLEDOUBLE,$(OPT)))
CC       =  mpiCC  -m64
endif
OPTIMIZE =   -g
GSL_INCL =
GSL_LIBS =
FFTW_INCL=
FFTW_LIBS=
MPICHLIB =
endif


ifeq ($(SYSTYPE),"Solaris")
CC       =   mpcc   # sets the C-compiler
#OPTIMIZE =   -fast -xO3 -xtarget=ultra3i -xcache=64/32/4:1024/64/4 -xarch=v9b
OPTIMIZE =  -g -xarch=v9b
R_PATH= -R/opt/local/gsl/64/current/lib:/opt/local/fftw/64/current/lib
#GSL_INCL = $(GSL64INCL)
#GSL_LIBS = $(R_PATH) $(GSL64LIB)
#FFTW_INCL= $(FFTW64INCL)
#FFTW_LIBS= $(FFTW64LIB)
GSL_INCL = -I/opt/local/gsl/64/current/include
GSL_LIBS = $(R_PATH) -L/opt/local/gsl/64/current/lib
FFTW_INCL= -I/opt/local/fftw/64/current/include
FFTW_LIBS= -L/opt/local/fftw/64/current/lib
MPICHLIB = -lmpi
HDF5INCL = -I/opt/local/hdf5/include/sparcv9
HDF5LIB  = -L/opt/local/hdf5/lib/sparcv9 -static -lhdf5 -lz
endif


ifeq ($(SYSTYPE),"OpteronMPA")
CC       =  mpicc  -m64 # sets the C-compiler
ifeq (SOFTDOUBLEDOUBLE,$(findstring SOFTDOUBLEDOUBLE,$(OPT)))
CC       =  mpiCC  -m64
endif
OPTIMIZE =   -O3 -Wall
GSL_INCL =  -L/usr/local/include
GSL_LIBS =  -L/usr/local/lib -static
FFTW_INCL=
FFTW_LIBS=
MPICHLIB =
HDF5INCL = -I/opt/hdf5/include
HDF5LIB  = -L/opt/hdf5/lib -static -lhdf5 -lz
endif


ifeq ($(SYSTYPE),"MPA")
CC       =  mpicc   # sets the C-compiler
ifeq (SOFTDOUBLEDOUBLE,$(findstring SOFTDOUBLEDOUBLE,$(OPT)))
CC       =  mpiCC
OPT     +=  -DX86FIX
endif
OPTIMIZE =   -g -Wall 
GSL_INCL = -I/usr/common/pdsoft/include
GSL_LIBS = -L/usr/common/pdsoft/lib
FFTW_INCL=
FFTW_LIBS=
MPICHLIB =
HDF5INCL = -I/usr/common/pdsoft/appl/hdf/include
HDF5LIB  = -L/usr/common/pdsoft/appl/hdf/lib -L/usr/common/pdsoft/appl/zlib/lib -static -lhdf5 -lz
endif


ifeq ($(SYSTYPE),"Mako")
CC       =  mpicc   # sets the C-compiler
ifeq (SOFTDOUBLEDOUBLE,$(findstring SOFTDOUBLEDOUBLE,$(OPT)))
CC       =  mpiCC
OPT     +=  -DX86FIX
endif
OPTIMIZE =   -O3 -march=athlon-mp  -mfpmath=sse
GSL_INCL =
GSL_LIBS =
FFTW_INCL=
FFTW_LIBS=
MPICHLIB =
endif


ifeq ($(SYSTYPE),"Regatta")
CC       =  mpcc_r -g # -qflttrap=enable:zerodivide:nanq # sets the C-compiler
ifeq (FLTROUNDOFFREDUCTION,$(findstring FLTROUNDOFFREDUCTION,$(OPT)))
CC       =  mpcc_r  -qldbl128 -lC128_r  # sets the C-compiler
#            (this compiler has native support for 128bit long double, SOFTDOUBLEDOUBLE not needed)
endif
OPTIMIZE =  -qstrict -q64 -qcpluscmt  -O3 -qipa
GSL_INCL = -I/afs/rzg/u/vrs/gsl_psi64/include
GSL_LIBS = -L/afs/rzg/u/vrs/gsl_psi64/lib
FFTW_INCL= -I/afs/rzg/u/vrs/fftw_psi64/include
FFTW_LIBS= -L/afs/rzg/u/vrs/fftw_psi64/lib  -q64 # -qipa
MPICHLIB =
HDF5INCL = -I/afs/rzg/u/vrs/hdf5_psi64/include
HDF5LIB  = -L/afs/rzg/u/vrs/hdf5_psi64/lib  -lhdf5 -lz
endif


ifeq ($(SYSTYPE),"CINECA32")
CC       =  mpcc_r   # sets the C-compiler
OPTIMIZE =  -O5 -qstrict -qipa -bmaxdata:500000000
GSL_INCL = -I/sfs/sanfs/home/userinaf/inapd006/include
GSL_LIBS = -L/sfs/sanfs/home/userinaf/inapd006/lib32
FFTW_INCL= -I/sfs/sanfs/home/userinaf/inapd006/include
FFTW_LIBS= -L/sfs/sanfs/home/userinaf/inapd006/lib32 -bmaxdata:500000000 -qipa
MPICHLIB =
endif

ifeq ($(SYSTYPE),"CINECA64")
CC       =   mpcc_r   # sets the C-compiler
OPTIMIZE =   -O5 -qstrict -qipa -q64
GSL_INCL = -I/sfs/sanfs/home/userinaf/inats004/include
GSL_LIBS = -L/sfs/sanfs/home/userinaf/inats004/lib
FFTW_INCL= -I/sfs/sanfs/home/userinaf/inats004/include
FFTW_LIBS= -L/sfs/sanfs/home/userinaf/inats004/lib -q64 -qipa
MPICHLIB =
endif


ifeq ($(SYSTYPE),"DEI32")
CC       =  mpcc   # sets the C-compiler
OPTIMIZE =  -O3 -qarch=pwr3 -qtune=pwr3 -qstrict -bmaxdata:1000000000
GSL_INCL = -I/home/kdolag/include
GSL_LIBS = -L/home/kdolag/lib
FFTW_INCL= -I/home/kdolag/include
FFTW_LIBS= -L/home/kdolag/lib -bmaxdata:1000000000
MPICHLIB =
endif


ifeq ($(SYSTYPE),"RZG_LinuxCluster")
CC       =   mpicci   # sets the C-compiler
ifeq (SOFTDOUBLEDOUBLE,$(findstring SOFTDOUBLEDOUBLE,$(OPT)))
CC       =   mpiCCi
OPT     +=  -DX86FIX
endif
OPTIMIZE =   -O3 -ip #-Qlongdouble
# Don't use the "-rcd" optimization of Intel's compiler! Makes the code crash at times!
GSL_INCL = -I/afs/rzg/u/vrs/gsl_linux/include
GSL_LIBS = -L/afs/rzg/u/vrs/gsl_linux/lib -Wl,"-R /afs/rzg/u/vrs/gsl_linux/lib"
FFTW_INCL= -I/afs/rzg/u/vrs/fftw_linux/include
FFTW_LIBS= -L/afs/rzg/u/vrs/fftw_linux/lib
HDF5INCL = -I/afs/rzg/u/vrs/hdf5_linux/include
HDF5LIB  = -L/afs/rzg/u/vrs/hdf5_linux/lib -lhdf5 -lz -Wl,"-R /afs/rzg/u/vrs/hdf5_linux/lib"
OPT      += -DX86FIX
endif


ifeq ($(SYSTYPE),"RZG_LinuxCluster-gcc")
CC       =   mpiccg -g  # sets the gcc C-compiler
ifeq (SOFTDOUBLEDOUBLE,$(findstring SOFTDOUBLEDOUBLE,$(OPT)))
CC       =   mpiCCg
OPT     +=  -DX86FIX
endif
OPTIMIZE =  -Wall    #-march=pentium4
GSL_INCL = -I/afs/rzg/u/vrs/gsl_linux_gcc3.2/include
GSL_LIBS = -L/afs/rzg/u/vrs/gsl_linux_gcc3.2/lib -Wl,"-R /afs/rzg/u/vrs/gsl_linux_gcc3.2/lib"
FFTW_INCL= -I/afs/rzg/u/vrs/fftw_linux_gcc3.2/include
FFTW_LIBS= -L/afs/rzg/u/vrs/fftw_linux_gcc3.2/lib
HDF5INCL = -I/afs/rzg/u/vrs/hdf5_linux/include
HDF5LIB  = -L/afs/rzg/u/vrs/hdf5_linux/lib -lhdf5 -lz -Wl,"-R /afs/rzg/u/vrs/hdf5_linux/lib"
OPT      +=  -DX86FIX
endif

ifeq ($(SYSTYPE),"Ingeld_LinuxCluster")
CC       =   mpicc  # sets the C-compiler
OPTIMIZE =   -O3 -Wall

GSL_INCL = -I/home/patricia/my-gsl/include
GSL_LIBS = -L/home/patricia/my-gsl/lib -static
FFTW_INCL= -I/home/patricia/my-fftw/include
FFTW_LIBS= -L/home/patricia/my-fftw/lib
endif

ifeq ($(SYSTYPE),"hpcf")
CC       =   mpicc  # sets the C-compiler
OPT     +=  -DFIX_PATHSCALE_MPI_STATUS_IGNORE_BUG
OPTIMIZE =  -O3
GSL_INCL = -I/home/gadget/Libs/include
GSL_LIBS = -L/home/gadget/Libs/lib
FFTW_INCL= -I/home/gadget/Libs/include
FFTW_LIBS= -L/home/gadget/Libs/lib
endif


ifneq (HAVE_HDF5,$(findstring HAVE_HDF5,$(OPT)))
HDF5INCL =
HDF5LIB  =
endif


OPTIONS = $(OPTIMIZE) $(OPT)

#EXEC = main_evapR_dust
#EXEC = main_evapR_dust_sinks
EXEC = main_evapR_dust_sink
#EXEC = main_evapR
#EXEC = main_poly

OBJS   = fof.o subfind.o subfind_vars.o subfind_collective.o subfind_serial.o subfind_so.o subfind_cont.o \
	 subfind_distribute.o subfind_findlinkngb.o subfind_nearesttwo.o subfind_loctree.o \
	 domain.o allvars.o main.o greenf_diffusion.o  \
	 subfind_potential.o subfind_density.o lineofsight.o kinfb_mhm.o sfr_mhm.o \
	 run.o predict.o begrun.o endrun.o global.o chemistry_noneq.o blackhole.o \
	 timestep.o init.o restart.o io.o sfr_eff.o fb_particles.o momentum_feedback.o \
	 accel.o read_ic.o cooling.o ngb.o parallel_sort.o read_ic_cluster.o \
	 system.o allocate.o density.o bsmooth.o bubbles.o  \
	 gravtree.o hydra.o eddington.o radtransfer.o get_Je.o driftfac.o darkenergy.o \
	 potential.o  forcetree.o  peano.o gravtree_forcetest.o \
	 pm_periodic.o pm_nonperiodic.o longrange.o mymalloc.o dust_density_on_gas.o \
	 cosmic_rays.o machfinder.o b_from_rot_a.o smooth_simple.o othin_accelerator.o ene_per_gram.o dust.o
###virtual_density.o virtual_spread_momentum.o \

INCL   = allvars.h proto.h forcetree.h cooling.h domain.h  cosmic_rays.h chemistry.h \
	 machfinder.h subfind.h dd.h fof.h Makefile

CFLAGS = $(OPTIONS) $(GSL_INCL) $(FFTW_INCL) $(HDF5INCL)

#ifeq (NOTYPEPREFIX_FFTW,$(findstring NOTYPEPREFIX_FFTW,$(OPT)))    # fftw installed with type prefix?
#  FFTW_LIB = $(FFTW_LIBS) -lrfftw_mpi -lfftw_mpi -lrfftw -lfftw
#else
#ifeq (DOUBLEPRECISION_FFTW,$(findstring DOUBLEPRECISION_FFTW,$(OPT)))
#  FFTW_LIB = $(FFTW_LIBS) -ldrfftw_mpi -ldfftw_mpi -ldrfftw -ldfftw
#else
#  FFTW_LIB = $(FFTW_LIBS) -lsrfftw_mpi -lsfftw_mpi -lsrfftw -lsfftw
#endif
#endif


LIBS   = -lm $(HDF5LIB) -g $(MPICHLIB) $(GSL_LIBS) -lgsl -lgslcblas $(FFTW_LIB)
#LIBS   = -lm $(HDF5LIB) -g $(MPICHLIB) $(GSL_LIBS) -lgsl -lgslcblas 


$(EXEC): $(OBJS)
	$(CC) $(OPTIONS) $(OBJS) $(LIBS) $(RLIBS) -o $(EXEC)

$(OBJS): $(INCL)

clean:
	rm -f $(OBJS) $(EXEC)




###############################################################################
#
# at compile-time. From the list below, please activate/deactivate the
# options that apply to your run. If you modify any of these options,
# make sure that you recompile the whole code by typing "make clean;
# make".
#
# Main code options:
#
#     These affect the physical model that is simulated.
#
#     - PERIODIC:   Set this if you want to have periodic boundary conditions.
#     - COOLING:    This enables radiative cooling and heating. It also enables
#                   an external UV background which is read from a file.
#     - SFR:        This enables star formation using an effective multiphase
#                   models. This option requires cooling.
#     - METALS:     This model activates the tracking of enrichment in gas and
#                   stars. Note that metal-line cooling is not included yet.
#     - STELLARAGE: This stores the formation redshift of each star particle.
#     - WINDS:      This activates galactic winds. Requires star formation.
#     - ISOTROPICWINDS: This makes the wind isotropic. If not set the wind is
#                       spawned in an axial way. Requires winds to be activated.
#     - NOGRAVITY:  This switches off gravity. Makes only sense for pure
#                   SPH simulations in non-expanding space.
#
# Options for SPH:
#
#     - NOFIXEDMASSINKERNEL:  If set, the number of SPH particles in the kernel
#                             is kept constant instead of the mass.
#     - NOGRADHSML:           If actived, an equation of motion without grad(h)
#                             terms is used.
#            Note: To have the default "entropy"-formulation of SPH (Springel &
#                  Hernquist), the switches NOFIXEDMASSINKERNEL and NOGRADHSML
#                  should *not* be set.
#     - NOVISCOSITYLIMITER:   If this is set, there is no explicit upper limit
#                             on the viscosity that tries to prevent particle
#                             'reflection' in case of poor timestepping.
#
# Numerical options:
#
#     - PMGRID:     This enables the TreePM method, i.e. the long-range force
#                   is computed with a PM-algoritthm, and the short range force
#                   with the tree. The parameter has to be set to the size of the
#                   mesh that should be used, (e.g. 64, 96, 128, etc). The mesh
#                   dimensions need not necessarily be a power of two.
#                   Note: If the simulation is not in a periodic box, then a FFT
#                   method for vacuum boundaries is employed, using a mesh with
#                   dimension twice that specified by PMGRID.
#     - PLACEHIGHRESREGION: If this option is set (will only work together
#                   with PMGRID), then the long range force is computed in two
#                   stages: One Fourier-grid is used to cover the whole simulation
#                   volume, allowing the computation of the large-scale force.
#                   A second Fourier mesh is placed on the region occupied by
#                   "high-resolution" particles, allowing the computation of an
#                   intermediate scale force. Finally, the force on very small
#                   scales is supplemented by the tree. This procedure can be useful
#                   for "zoom-simulations", where the majority of particles (the
#                   high-res particles) are occupying only a small fraction of the
#                   volume. To activate this option, the parameter needs to be set
#                   to an integer that encodes the particle types that represent the
#                   high-res particles in the form of a bit mask. For example, if
#                   types 0, 1, and 4 form the high-res particles, set the parameter
#                   to PLACEHIGHRESREGION=1+2+16. The spatial region covered by the
#                   high-res grid is determined automatically from the initial
#                   conditions. Note: If a periodic box is used, the high-res zone
#                   may not intersect the box boundaries.
#     - ENLARGEREGION: The spatial region covered by the high-res zone has a fixed
#                   size during the simulation, which initially is set to the
#                   smallest region that encompasses all high-res particles. Normally, the
#                   simulation will be interrupted, if high-res particles leave this
#                   region in the course of the run. However, by setting this parameter
#                   to a value larger than one, the high-res region can be expanded.
#                   For example, setting it to 1.4 will enlarge its side-length by
#                   40% (it remains centered on the high-res particles). Hence, with
#                   such a setting, the high-res region may expand or move by a
#                   limited amount. If in addition SYNCHRONIZATION is activated, then
#                   the code will be able to continue even if high-res particles
#                   leave the initial high-res grid. In this case, the code will
#                   update the size and position of the grid that is placed onto
#                   the high-resolution region automatically. To prevent that this
#                   potentially happens every single PM step, one should nevertheless
#                   assign a value slightly larger than 1 to ENLARGEREGION.
#     - DOUBLEPRECISION: This makes the code store and compute internal
#                        particle data in double precision. Note that output
#                        files are nevertheless written by converting to single
#                        precision.
#     - NOTREERND:       If this is not set, the tree construction will succeed
#                        even when there are a few particles at identical
#                        locations. This is done by `rerouting' particles once
#                        the node-size has fallen below 1.0e-3 of the softening
#                        length. When this option is activated, this will be
#                        surpressed and the tree construction will always fail
#                        if there are particles at extremely close coordinates.
#     - NOSTOP_WHEN_BELOW_MINTIMESTEP: If this is activated, the code will not
#                        terminate when the timestep falls below the value of
#                        MinSizeTimestep specified in the parameterfile. This
#                        is useful for runs where one wants to enforce a
#                        constant timestep for all particles. This can be done
#                        by activating this option, and by setting Min- and
#                        MaxSizeTimestep to an equal value.
#     - PSEUDOSYMMETRIC: When this option is set, the code will try to "anticipate"
#                        timestep changes by extrapolating the change of the
#                        acceleration into the future. This in general improves the
#                        long-term integration behaviour of periodic orbits.
#     - SYNCHRONIZATION: When this is set, particles may only increase their
#                        timestep if the new timestep will put them into
#                        synchronization with the higher time level. This typically
#                        means that only on half of the timesteps of a particle
#                        an increase of the step may occur.
#     - NOPMSTEPADJUSTMENT: When this is set, the long-range timestep for the
#                        PM force computation is always determined by MaxSizeTimeStep.
#                        Otherwise, it is set to the minimum of MaxSizeTimeStep and
#                        the timestep obtained for the maximum long-range force with
#                        an effective softening scale equal to the PM smoothing-scale.
# - LONG_X/Y/Z:
#     These options can be used together with PERIODIC and NOGRAVITY only.
#     When set, the options define numerical factors that can be used to
#     distorts the periodic simulation cube into a parallelepiped of
#     arbitrary aspect ratio. This can be useful for idealized SPH tests.
#
# - TWODIMS:
#     This effectively switches of one dimension in SPH, i.e. the code
#     follows only 2d hydrodynamics in the xy-, yz-, or xz-plane. This
#     only works with NOGRAVITY, and if all coordinates of the third
#     axis are exactly equal. Can be useful for idealized SPH tests.
#
# - SPH_BND_PARTICLES:
#     If this is set, particles with a particle-ID equal to zero do not
#     receive any SPH acceleration. This can be useful for idealized
#     SPH tests, where these particles represent fixed "walls".
#
#
# Architecture options:
#
#     - T3E:       The code assumes that sizeof(int)=4 holds. A few machines
#                  (like Cray T3E) have sizeof(int)=8. In this case, set the
#                  T3E flag.
#     - NOTYPEPREFIX_FFTW: If this is set, the fftw-header/libraries are accessed
#                  without type prefix (adopting whatever was chosen as default at compile
#                  of fftw). Otherwise, the type prefix 'd' for double is used.
#
# Input options:
#
#     - MOREPARAMS:  Activate this to allow a set of additional parameters in
#                    the parameterfile which control the star formation and
#                    feedback sector. This option must be activated when star
#                    formation is switched on.
#
# Output options:
#
#     - OUTPUTPOTENTIAL: This will force the code to compute gravitational
#                        potentials for all particles each time a snapshot file
#                        is generated. This values are then included in the
#                        snapshot file. Note that the computation of the
#                        values of the potential costs additional time.
#     - OUTPUTACCELERATION: This will include the physical acceleration of
#                        each particle in snapshot files.
#     - OUTPUTCHANGEOFENTROPY: This will include the rate of change of entropy
#                        of gas particles in snapshot files.
#     - OUTPUTTIMESTEP:  This will include an output of the timesteps actually
#                        taken by each particle.
#     - OUTPUT_D2POT: This will output the second derivatives of the potential
#                     for each particle in the snapshot file.
#
# Miscellaneous options:
#
#     - PEANOHILBERT:    This is a tuning option. When set, the code will bring
#                        the particles after each domain decomposition into
#                        Peano-Hilbert order. This improves cache utilization
#                        and performance.
#     - WALLCLOCK:       If set, a wallclock timer is used by the code to
#                        measure internal time consumption (see cpu-log file).
#                        Otherwise a timer that measures consumed processor
#                        ticks is used.
#
# Debugging/testing options:
#
#     - FORCETEST:       This can be set to check the force accuracy of the
#                        code. The option needs to be set to a number between
#                        0 and 1 (e.g. 0.01), which is taken to specify a
#                        random fraction of particles for which at each
#                        timestep forces by direct summation are computed. The
#                        normal tree-forces and the "correct" direct summation
#                        forces are collected in a file. Note that the
#                        simulation itself is unaffected by this option, but it
#                        will of course run much(!) slower
#                        if FORCETEST*NumPart*NumPart >> NumPart. Note: Particle
#                        IDs must be set to numbers >=1 for this to work.
#
###############################################################################




# - FLTROUNDOFFREDUCTION     enables round off reduction in particle sums
#		             if DOUBLEPRECISION is set, these sums are done in 'long double'
#                            if single precision is used, they are done in 'double'
#                            This should in principle allow to make computations
#                            *exactly* invariant to different numbers of CPUs.
#
# - SOFTDOUBLEDOUBLE         when this is set, a software implementation of
#                            128bit double-double addition is used, implemented as a c++ class.
#                            Hence, this option requires compilation with a c++ compiler
#


#     - QUICK_LYALPHA:   This only works for cosmological simulations in periodic boxes
#                        with COOLING & SFR. (WINDS, METALS should be deselected).
#                        It will simply convert all gas particles above overdensity
#                        CritPhysOverdensity and with Temperature below 10^5 K to stars.
#                        This should still leave the Ly-Alpha forest largely unaffected,
#                        but should be faster. It is recommended to set GENERATIONS equal
#                        to 1 for maximum speed-up.
#
#     - TIDALTENSOR:     Calculates the tidal tensor for each particle. 
