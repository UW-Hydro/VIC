headers = '''
FILE *LOG_DEST;
void finalize_logging(void);
void get_current_datetime(char *cdt);
void get_logname(const char *path, int id, char *filename);
void initialize_log(void);
void setup_logging(int id);
extern size_t NR;
extern size_t NF;
enum
{
    ASCII,
    BINARY
};
enum
{
    LITTLE,
    BIG
};
enum
{
    DENS_BRAS,
    DENS_SNTHRM
};
enum
{
    ARNO,
    NIJSSEN2001
};
enum
{
    AR_406,
    AR_406_LS,
    AR_406_FULL,
    AR_410
};
enum
{
    GF_406,
    GF_410
};
enum
{
    FROM_DEFAULT,
    FROM_VEGLIB,
    FROM_VEGPARAM,
    FROM_VEGHIST
};
enum
{
    RC_JARVIS,
    RC_PHOTO
};
enum
{
    PS_FARQUHAR,
    PS_MONTEITH
};
enum
{
    PHOTO_C3,
    PHOTO_C4
};
enum
{
    AIR_TEMP,
    ALBEDO,
    CATM,
    CHANNEL_IN,
    DENSITY,
    FCANOPY,
    FDIR,
    LAI_IN,
    LWDOWN,
    PAR,
    PREC,
    PRESSURE,
    QAIR,
    REL_HUMID,
    SWDOWN,
    VP,
    WIND,
    SKIP,
    N_FORCING_TYPES
};
enum
{
    OUT_ASAT,
    OUT_LAKE_AREA_FRAC,
    OUT_LAKE_DEPTH,
    OUT_LAKE_ICE,
    OUT_LAKE_ICE_FRACT,
    OUT_LAKE_ICE_HEIGHT,
    OUT_LAKE_MOIST,
    OUT_LAKE_SURF_AREA,
    OUT_LAKE_SWE,
    OUT_LAKE_SWE_V,
    OUT_LAKE_VOLUME,
    OUT_ROOTMOIST,
    OUT_SMFROZFRAC,
    OUT_SMLIQFRAC,
    OUT_SNOW_CANOPY,
    OUT_SNOW_COVER,
    OUT_SNOW_DEPTH,
    OUT_SOIL_ICE,
    OUT_SOIL_LIQ,
    OUT_SOIL_MOIST,
    OUT_SOIL_WET,
    OUT_SURFSTOR,
    OUT_SURF_FROST_FRAC,
    OUT_SWE,
    OUT_WDEW,
    OUT_ZWT,
    OUT_ZWT_LUMPED,
    OUT_BASEFLOW,
    OUT_DELINTERCEPT,
    OUT_DELSOILMOIST,
    OUT_DELSURFSTOR,
    OUT_DELSWE,
    OUT_EVAP,
    OUT_EVAP_BARE,
    OUT_EVAP_CANOP,
    OUT_INFLOW,
    OUT_LAKE_BF_IN,
    OUT_LAKE_BF_IN_V,
    OUT_LAKE_BF_OUT,
    OUT_LAKE_BF_OUT_V,
    OUT_LAKE_CHAN_IN,
    OUT_LAKE_CHAN_IN_V,
    OUT_LAKE_CHAN_OUT,
    OUT_LAKE_CHAN_OUT_V,
    OUT_LAKE_DSTOR,
    OUT_LAKE_DSTOR_V,
    OUT_LAKE_DSWE,
    OUT_LAKE_DSWE_V,
    OUT_LAKE_EVAP,
    OUT_LAKE_EVAP_V,
    OUT_LAKE_PREC_V,
    OUT_LAKE_RCHRG,
    OUT_LAKE_RCHRG_V,
    OUT_LAKE_RO_IN,
    OUT_LAKE_RO_IN_V,
    OUT_LAKE_VAPFLX,
    OUT_LAKE_VAPFLX_V,
    OUT_PET,
    OUT_PREC,
    OUT_RAINF,
    OUT_REFREEZE,
    OUT_RUNOFF,
    OUT_SNOW_MELT,
    OUT_SNOWF,
    OUT_SUB_BLOWING,
    OUT_SUB_CANOP,
    OUT_SUB_SNOW,
    OUT_SUB_SURFACE,
    OUT_TRANSP_VEG,
    OUT_WATER_ERROR,
    OUT_ALBEDO,
    OUT_BARESOILT,
    OUT_FDEPTH,
    OUT_LAKE_ICE_TEMP,
    OUT_LAKE_SURF_TEMP,
    OUT_RAD_TEMP,
    OUT_SALBEDO,
    OUT_SNOW_PACK_TEMP,
    OUT_SNOW_SURF_TEMP,
    OUT_SNOWT_FBFLAG,
    OUT_SOIL_TEMP,
    OUT_SOIL_TNODE,
    OUT_SOIL_TNODE_WL,
    OUT_SOILT_FBFLAG,
    OUT_SURF_TEMP,
    OUT_SURFT_FBFLAG,
    OUT_TCAN_FBFLAG,
    OUT_TDEPTH,
    OUT_TFOL_FBFLAG,
    OUT_VEGT,
    OUT_ADV_SENS,
    OUT_ADVECTION,
    OUT_DELTACC,
    OUT_DELTAH,
    OUT_ENERGY_ERROR,
    OUT_FUSION,
    OUT_GRND_FLUX,
    OUT_IN_LONG,
    OUT_LATENT,
    OUT_LATENT_SUB,
    OUT_MELT_ENERGY,
    OUT_LWNET,
    OUT_SWNET,
    OUT_R_NET,
    OUT_RFRZ_ENERGY,
    OUT_SENSIBLE,
    OUT_SNOW_FLUX,
    OUT_AERO_COND,
    OUT_AERO_COND1,
    OUT_AERO_COND2,
    OUT_AERO_RESIST,
    OUT_AERO_RESIST1,
    OUT_AERO_RESIST2,
    OUT_AIR_TEMP,
    OUT_CATM,
    OUT_COSZEN,
    OUT_DENSITY,
    OUT_FCANOPY,
    OUT_FDIR,
    OUT_LAI,
    OUT_LWDOWN,
    OUT_PAR,
    OUT_PRESSURE,
    OUT_QAIR,
    OUT_REL_HUMID,
    OUT_SWDOWN,
    OUT_SURF_COND,
    OUT_TSKC,
    OUT_VP,
    OUT_VPD,
    OUT_WIND,
    OUT_ADV_SENS_BAND,
    OUT_ADVECTION_BAND,
    OUT_ALBEDO_BAND,
    OUT_DELTACC_BAND,
    OUT_GRND_FLUX_BAND,
    OUT_IN_LONG_BAND,
    OUT_LATENT_BAND,
    OUT_LATENT_SUB_BAND,
    OUT_MELT_ENERGY_BAND,
    OUT_LWNET_BAND,
    OUT_SWNET_BAND,
    OUT_RFRZ_ENERGY_BAND,
    OUT_SENSIBLE_BAND,
    OUT_SNOW_CANOPY_BAND,
    OUT_SNOW_COVER_BAND,
    OUT_SNOW_DEPTH_BAND,
    OUT_SNOW_FLUX_BAND,
    OUT_SNOW_MELT_BAND,
    OUT_SNOW_PACKT_BAND,
    OUT_SNOW_SURFT_BAND,
    OUT_SWE_BAND,
    OUT_APAR,
    OUT_GPP,
    OUT_RAUT,
    OUT_NPP,
    OUT_LITTERFALL,
    OUT_RHET,
    OUT_NEE,
    OUT_CLITTER,
    OUT_CINTER,
    OUT_CSLOW,
    N_OUTVAR_TYPES
};
enum
{
    OUT_TYPE_DEFAULT,
    OUT_TYPE_CHAR,
    OUT_TYPE_SINT,
    OUT_TYPE_USINT,
    OUT_TYPE_INT,
    OUT_TYPE_FLOAT,
    OUT_TYPE_DOUBLE
};
enum
{
    AGG_TYPE_AVG,
    AGG_TYPE_BEG,
    AGG_TYPE_END,
    AGG_TYPE_MAX,
    AGG_TYPE_MIN,
    AGG_TYPE_SUM
};
enum
{
    DISP_VERSION,
    DISP_COMPILE_TIME,
    DISP_ALL
};
enum calendars
{
    CALENDAR_STANDARD,
    CALENDAR_GREGORIAN,
    CALENDAR_PROLEPTIC_GREGORIAN,
    CALENDAR_NOLEAP,
    CALENDAR_365_DAY,
    CALENDAR_360_DAY,
    CALENDAR_JULIAN,
    CALENDAR_ALL_LEAP,
    CALENDAR_366_DAY
};
enum time_units
{
    TIME_UNITS_SECONDS,
    TIME_UNITS_MINUTES,
    TIME_UNITS_HOURS,
    TIME_UNITS_DAYS
};
typedef struct {
    FILE *forcing[2];
    FILE *globalparam;
    FILE *constants;
    FILE *domain;
    FILE *init_state;
    FILE *lakeparam;
    FILE *snowband;
    FILE *soilparam;
    FILE *statefile;
    FILE *veglib;
    FILE *vegparam;
    FILE *logfile;
} filep_struct;
typedef struct {
    char forcing[2][2048];
    char f_path_pfx[2][2048];
    char global[2048];
    char domain[2048];
    char constants[2048];
    char init_state[2048];
    char lakeparam[2048];
    char result_dir[2048];
    char snowband[2048];
    char soil[2048];
    char statefile[2048];
    char veg[2048];
    char veglib[2048];
    char log_path[2048];
} filenames_struct;
typedef struct {
    short AboveTreelineVeg;
    unsigned short int AERO_RESIST_CANSNOW;
    _Bool BLOWING;
    _Bool BLOWING_VAR_THRESHOLD;
    _Bool BLOWING_CALC_PROB;
    _Bool BLOWING_SIMPLE;
    _Bool BLOWING_FETCH;
    _Bool BLOWING_SPATIAL_WIND;
    _Bool CARBON;
    _Bool CLOSE_ENERGY;
    _Bool COMPUTE_TREELINE;
    _Bool CONTINUEONERROR;
    _Bool CORRPREC;
    _Bool EQUAL_AREA;
    _Bool EXP_TRANS;
    _Bool FROZEN_SOIL;
    _Bool FULL_ENERGY;
    unsigned short int GRND_FLUX_TYPE;
    _Bool IMPLICIT;
    _Bool JULY_TAVG_SUPPLIED;
    _Bool LAKES;
    size_t Ncanopy;
    size_t Nfrost;
    size_t Nlakenode;
    size_t Nlayer;
    size_t Nnode;
    _Bool NOFLUX;
    size_t NVEGTYPES;
    unsigned short int RC_MODE;
    size_t ROOT_ZONES;
    _Bool QUICK_FLUX;
    _Bool QUICK_SOLVE;
    _Bool SHARE_LAYER_MOIST;
    unsigned short int SNOW_DENSITY;
    size_t SNOW_BAND;
    _Bool SPATIAL_FROST;
    _Bool SPATIAL_SNOW;
    _Bool TFALLBACK;
    _Bool BASEFLOW;
    unsigned short int GRID_DECIMAL;
    _Bool VEGLIB_FCAN;
    _Bool VEGLIB_PHOTO;
    _Bool VEGPARAM_ALB;
    _Bool VEGPARAM_FCAN;
    _Bool VEGPARAM_LAI;
    unsigned short int ALB_SRC;
    unsigned short int FCAN_SRC;
    unsigned short int LAI_SRC;
    _Bool LAKE_PROFILE;
    _Bool ORGANIC_FRACT;
    _Bool BINARY_STATE_FILE;
    _Bool INIT_STATE;
    _Bool SAVE_STATE;
    _Bool ALMA_OUTPUT;
    _Bool BINARY_OUTPUT;
    _Bool COMPRESS;
    _Bool MOISTFRACT;
    size_t Noutfiles;
    _Bool PRT_HEADER;
    _Bool PRT_SNOW_BAND;
} option_struct;
typedef struct {
    double wind_h;
    double resolution;
    double dt;
    double snow_dt;
    double runoff_dt;
    double atmos_dt;
    double out_dt;
    size_t model_steps_per_day;
    size_t snow_steps_per_day;
    size_t runoff_steps_per_day;
    size_t atmos_steps_per_day;
    size_t output_steps_per_day;
    unsigned short int endday;
    unsigned short int endmonth;
    unsigned short int endyear;
    unsigned short int forceday[2];
    unsigned int forcesec[2];
    unsigned short int forcemonth[2];
    unsigned short int forceoffset[2];
    unsigned int forceskip[2];
    unsigned short int forceyear[2];
    size_t nrecs;
    unsigned short int skipyear;
    unsigned short int startday;
    unsigned int startsec;
    unsigned short int startmonth;
    unsigned short int startyear;
    unsigned short int stateday;
    unsigned short int statemonth;
    unsigned short int stateyear;
    unsigned short int calendar;
    unsigned short int time_units;
    double time_origin_num;
} global_param_struct;
typedef struct {
    double LAPSE_RATE;
    double GAUGE_HEIGHT;
    double HUGE_RESIST;
    double ALBEDO_BARE_SOIL;
    double ALBEDO_H20_SURF;
    double EMISS_GRND;
    double EMISS_VEG;
    double EMISS_ICE;
    double EMISS_SNOW;
    double EMISS_H2O;
    double SOIL_RARC;
    double SOIL_RESID_MOIST;
    double SOIL_SLAB_MOIST_FRACT;
    double SOIL_WINDH;
    double VEG_LAI_SNOW_MULTIPLIER;
    double VEG_MIN_INTERCEPTION_STORAGE;
    double VEG_LAI_WATER_FACTOR;
    double VEG_RATIO_DH_HEIGHT;
    double VEG_RATIO_RL_HEIGHT;
    double CANOPY_CLOSURE;
    double CANOPY_RSMAX;
    double CANOPY_VPDMINFACTOR;
    double LAKE_TMELT;
    double LAKE_MAX_SURFACE;
    double LAKE_BETA;
    double LAKE_FRACMIN;
    double LAKE_FRACLIM;
    double LAKE_DM;
    double LAKE_SNOWCRIT;
    double LAKE_ZWATER;
    double LAKE_ZSNOW;
    double LAKE_RHOSNOW;
    double LAKE_CONDI;
    double LAKE_CONDS;
    double LAKE_LAMISW;
    double LAKE_LAMILW;
    double LAKE_LAMSSW;
    double LAKE_LAMSLW;
    double LAKE_LAMWSW;
    double LAKE_LAMWLW;
    double LAKE_A1;
    double LAKE_A2;
    double LAKE_QWTAU;
    int LAKE_MAX_ITER;
    double SVP_A;
    double SVP_B;
    double SVP_C;
    double PHOTO_OMEGA;
    double PHOTO_LAIMAX;
    double PHOTO_LAILIMIT;
    double PHOTO_LAIMIN;
    double PHOTO_EPAR;
    double PHOTO_FCMAX;
    double PHOTO_FCMIN;
    double PHOTO_ZENITHMIN;
    double PHOTO_ZENITHMINPAR;
    double PHOTO_ALBSOIPARMIN;
    double PHOTO_MINMAXETRANS;
    double PHOTO_MINSTOMCOND;
    double PHOTO_FCI1C3;
    double PHOTO_FCI1C4;
    double PHOTO_OX;
    double PHOTO_KC;
    double PHOTO_KO;
    double PHOTO_EC;
    double PHOTO_EO;
    double PHOTO_EV;
    double PHOTO_ER;
    double PHOTO_ALC3;
    double PHOTO_FRDC3;
    double PHOTO_EK;
    double PHOTO_ALC4;
    double PHOTO_FRDC4;
    double PHOTO_THETA;
    double PHOTO_FRLEAF;
    double PHOTO_FRGROWTH;
    double SRESP_E0_LT;
    double SRESP_T0_LT;
    double SRESP_WMINFM;
    double SRESP_WMAXFM;
    double SRESP_WOPTFM;
    double SRESP_RHSAT;
    double SRESP_RFACTOR;
    double SRESP_TAULITTER;
    double SRESP_TAUINTER;
    double SRESP_TAUSLOW;
    double SRESP_FAIR;
    double SRESP_FINTER;
    double SNOW_MAX_SURFACE_SWE;
    double SNOW_LIQUID_WATER_CAPACITY;
    double SNOW_NEW_SNOW_DENSITY;
    double SNOW_DENS_DMLIMIT;
    double SNOW_DENS_MAX_CHANGE;
    double SNOW_DENS_ETA0;
    double SNOW_DENS_C1;
    double SNOW_DENS_C2;
    double SNOW_DENS_C5;
    double SNOW_DENS_C6;
    double SNOW_DENS_F;
    double SNOW_MIN_SWQ_EB_THRES;
    double SNOW_A1;
    double SNOW_A2;
    double SNOW_L1;
    double SNOW_L2;
    double SNOW_NEW_SNOW_ALB;
    double SNOW_ALB_ACCUM_A;
    double SNOW_ALB_ACCUM_B;
    double SNOW_ALB_THAW_A;
    double SNOW_ALB_THAW_B;
    double SNOW_TRACESNOW;
    double SNOW_CONDUCT;
    double SNOW_MAX_SNOW_TEMP;
    double SNOW_MIN_RAIN_TEMP;
    double BLOWING_KA;
    double BLOWING_CSALT;
    double BLOWING_UTHRESH;
    double BLOWING_KIN_VIS;
    int BLOWING_MAX_ITER;
    int BLOWING_K;
    double BLOWING_SETTLING;
    int BLOWING_NUMINCS;
    double TREELINE_TEMPERATURE;
    double SNOW_DT;
    double SURF_DT;
    double SOIL_DT;
    double CANOPY_DT;
    double CANOPY_VP;
    double TOL_GRND;
    double TOL_OVER;
    int FROZEN_MAXITER;
    int NEWT_RAPH_MAXTRIAL;
    double NEWT_RAPH_TOLX;
    double NEWT_RAPH_TOLF;
    double NEWT_RAPH_R_MAX;
    double NEWT_RAPH_R_MIN;
    double NEWT_RAPH_RELAX1;
    double NEWT_RAPH_RELAX2;
    double NEWT_RAPH_RELAX3;
    double NEWT_RAPH_EPS2;
    int ROOT_BRENT_MAXTRIES;
    int ROOT_BRENT_MAXITER;
    double ROOT_BRENT_TSTEP;
    double ROOT_BRENT_T;
} parameters_struct;
typedef struct {
    _Bool FS_ACTIVE;
    double Ds;
    double Dsmax;
    double Ksat[3];
    double Wcr[3];
    double Wpwp[3];
    double Ws;
    double AlbedoPar;
    double alpha[50];
    double annual_prec;
    double avg_temp;
    double avgJulyAirTemp;
    double b_infilt;
    double beta[50];
    double bubble[3];
    double bubble_node[50];
    double bulk_density[3];
    double bulk_dens_min[3];
    double bulk_dens_org[3];
    double c;
    double depth[3];
    double dp;
    double dz_node[50];
    double Zsum_node[50];
    double expt[3];
    double expt_node[50];
    double frost_fract[10];
    double frost_slope;
    double gamma[50];
    double init_moist[3];
    double max_infil;
    double max_moist[3];
    double max_moist_node[50];
    double max_snow_distrib_slope;
    double phi_s[3];
    double porosity[3];
    double quartz[3];
    double organic[3];
    double resid_moist[3];
    double rough;
    double snow_rough;
    double soil_density[3];
    double soil_dens_min[3];
    double soil_dens_org[3];
    double *BandElev;
    double *AreaFract;
    double *Pfactor;
    double *Tfactor;
    _Bool *AboveTreeLine;
    double elevation;
    double lat;
    double lng;
    double cell_area;
    double time_zone_lng;
    unsigned int gridcel;
    double zwtvmoist_zwt[5][11];
    double zwtvmoist_moist[5][11];
    double slope;
    double aspect;
    double ehoriz;
    double whoriz;
} soil_con_struct;
typedef struct {
    double Cv;
    double Cv_sum;
    double root[3];
    double *zone_depth;
    double *zone_fract;
    int veg_class;
    size_t vegetat_type_num;
    double sigma_slope;
    double lag_one;
    double fetch;
    int LAKE;
    double *CanopLayerBnd;
} veg_con_struct;
typedef struct {
    _Bool overstory;
    double LAI[12];
    double fcanopy[12];
    double Wdmax[12];
    double albedo[12];
    double displacement[12];
    double emissivity[12];
    size_t NVegLibTypes;
    double rad_atten;
    double rarc;
    double rmin;
    double roughness[12];
    double trunk_ratio;
    double wind_atten;
    double wind_h;
    double RGL;
    unsigned short int veg_class;
    char Ctype;
    double MaxCarboxRate;
    double MaxETransport;
    double CO2Specificity;
    double LightUseEff;
    _Bool NscaleFlag;
    double Wnpp_inhib;
    double NPPfactor_sat;
} veg_lib_struct;
typedef struct {
    double *albedo;
    double *LAI;
    double *fcanopy;
} veg_hist_struct;
typedef struct {
    double *air_temp;
    double *Catm;
    double *channel_in;
    double *density;
    double *fdir;
    double *longwave;
    double out_prec;
    double out_rain;
    double out_snow;
    double *par;
    double *prec;
    double *pressure;
    double *shortwave;
    _Bool *snowflag;
    double *vp;
    double *vpd;
    double *wind;
} atmos_data_struct;
typedef struct {
    unsigned short int day;
    unsigned short int day_in_year;
    unsigned short int month;
    int year;
    unsigned int dayseconds;
} dmy_struct;
typedef struct {
    double Cs;
    double T;
    double bare_evap_frac;
    double evap;
    double ice[10];
    double kappa;
    double moist;
    double phi;
    double zwt;
} layer_data_struct;
typedef struct {
    double aero_resist[2];
    double asat;
    double baseflow;
    double CLitter;
    double CInter;
    double CSlow;
    double inflow;
    double pot_evap;
    double runoff;
    layer_data_struct layer[3];
    double RhLitter;
    double RhLitter2Atm;
    double RhInter;
    double RhSlow;
    double RhTot;
    double rootmoist;
    double wetness;
    double zwt;
    double zwt_lumped;
} cell_data_struct;
typedef struct {
    double AlbedoLake;
    double AlbedoOver;
    double AlbedoUnder;
    double Cs[2];
    double Cs_node[50];
    double fdepth[3];
    _Bool frozen;
    double ice[50];
    double kappa[2];
    double kappa_node[50];
    double moist[50];
    size_t Nfrost;
    size_t Nthaw;
    double T[50];
    _Bool T_fbflag[50];
    unsigned int T_fbcount[50];
    int T1_index;
    double Tcanopy;
    _Bool Tcanopy_fbflag;
    unsigned int Tcanopy_fbcount;
    double tdepth[3];
    double Tfoliage;
    _Bool Tfoliage_fbflag;
    unsigned int Tfoliage_fbcount;
    double Tsurf;
    _Bool Tsurf_fbflag;
    unsigned int Tsurf_fbcount;
    double unfrozen;
    double advected_sensible;
    double advection;
    double AtmosError;
    double AtmosLatent;
    double AtmosLatentSub;
    double AtmosSensible;
    double canopy_advection;
    double canopy_latent;
    double canopy_latent_sub;
    double canopy_refreeze;
    double canopy_sensible;
    double deltaCC;
    double deltaH;
    double error;
    double fusion;
    double grnd_flux;
    double latent;
    double latent_sub;
    double longwave;
    double LongOverIn;
    double LongUnderIn;
    double LongUnderOut;
    double melt_energy;
    double NetLongAtmos;
    double NetLongOver;
    double NetLongUnder;
    double NetShortAtmos;
    double NetShortGrnd;
    double NetShortOver;
    double NetShortUnder;
    double out_long_canopy;
    double out_long_surface;
    double refreeze_energy;
    double sensible;
    double shortwave;
    double ShortOverIn;
    double ShortUnderIn;
    double snow_flux;
} energy_bal_struct;
typedef struct {
    double albedo;
    double canopyevap;
    double LAI;
    double throughfall;
    double fcanopy;
    double Wdew;
    double Wdmax;
    double *NscaleFactor;
    double *aPARLayer;
    double *CiLayer;
    double *rsLayer;
    double aPAR;
    double Ci;
    double rc;
    double NPPfactor;
    double GPP;
    double Rphoto;
    double Rdark;
    double Rmaint;
    double Rgrowth;
    double Raut;
    double NPP;
    double Litterfall;
    double AnnualNPP;
    double AnnualNPPPrev;
} veg_var_struct;
typedef struct {
    double albedo;
    double canopy_albedo;
    double coldcontent;
    double coverage;
    double density;
    double depth;
    unsigned int last_snow;
    double max_snow_depth;
    _Bool MELTING;
    double pack_temp;
    double pack_water;
    _Bool snow;
    double snow_canopy;
    double store_coverage;
    _Bool store_snow;
    double store_swq;
    double surf_temp;
    unsigned int surf_temp_fbcount;
    _Bool surf_temp_fbflag;
    double surf_water;
    double swq;
    double snow_distrib_slope;
    double tmp_int_storage;
    double blowing_flux;
    double canopy_vapor_flux;
    double mass_error;
    double melt;
    double Qnet;
    double surface_flux;
    double transport;
    double vapor_flux;
} snow_data_struct;
typedef struct {
    size_t numnod;
    double z[21];
    double basin[21];
    double Cl[21];
    double b;
    double maxdepth;
    double mindepth;
    double maxvolume;
    double minvolume;
    double bpercent;
    double rpercent;
    double wfrac;
    double depth_in;
    int lake_idx;
} lake_con_struct;
typedef struct {
    unsigned short int activenod;
    double dz;
    double surfdz;
    double ldepth;
    double surface[21];
    double sarea;
    double sarea_save;
    double volume;
    double volume_save;
    double temp[20];
    double tempavg;
    double areai;
    double new_ice_area;
    double ice_water_eq;
    double hice;
    double tempi;
    double swe;
    double swe_save;
    double surf_temp;
    double pack_temp;
    double coldcontent;
    double surf_water;
    double pack_water;
    double SAlbedo;
    double sdepth;
    double aero_resist;
    double density[20];
    double baseflow_in;
    double baseflow_out;
    double channel_in;
    double evapw;
    double ice_throughfall;
    double prec;
    double recharge;
    double runoff_in;
    double runoff_out;
    double snowmlt;
    double vapor_flux;
    snow_data_struct snow;
    energy_bal_struct energy;
    cell_data_struct soil;
} lake_var_struct;
typedef struct {
    cell_data_struct **cell;
    energy_bal_struct **energy;
    lake_var_struct lake_var;
    snow_data_struct **snow;
    veg_var_struct **veg_var;
} all_vars_struct;
typedef struct {
    double total_soil_moist;
    double surfstor;
    double swe;
    double wdew;
} save_data_struct;
typedef struct {
    char varname[20];
    _Bool write;
    char format[10];
    unsigned short int type;
    double mult;
    unsigned short int aggtype;
    unsigned int nelem;
    double *data;
    double *aggdata;
} out_data_struct;
typedef struct {
    char prefix[20];
    char filename[2048];
    FILE *fh;
    size_t nvars;
    unsigned int *varid;
} out_data_file_struct;
typedef struct {
    atmos_data_struct *atmos;
    double dt;
    energy_bal_struct *energy;
    filep_struct filep;
    size_t rec;
    out_data_struct *out_data;
    out_data_file_struct *out_data_files;
    snow_data_struct *snow;
    soil_con_struct soil_con;
    veg_con_struct *veg_con;
    veg_var_struct *veg_var;
} Error_struct;
void advect_carbon_storage(double, double, lake_var_struct *,
                           cell_data_struct *);
void advect_snow_storage(double, double, double, snow_data_struct *);
void advect_soil_veg_storage(double, double, double, double *,
                             soil_con_struct *, veg_con_struct *,
                             cell_data_struct *, veg_var_struct *,
                             lake_con_struct);
double advected_sensible_heat(double, double, double, double, double);
void alblake(double, double, double *, double *, double *, double *, double,
             double, double, unsigned int *, double, _Bool *, unsigned short int,
             double);
double arno_evap(layer_data_struct *, double, double, double, double, double,
                 double, double, double, double, double, double *);
bool assert_close_double(double x, double y, double rtol, double abs_tol);
bool assert_close_float(float x, float y, float rtol, float abs_tol);
double calc_atmos_energy_bal(double, double, double, double, double, double,
                             double, double, double, double, double, double,
                             double, double *, double *, double *, double *,
                             double *, double *, _Bool *, unsigned int*);
double calc_density(double);
double calc_energy_balance_error(int, double, double, double, double, double);
double calc_latent_heat_of_sublimation(double temp);
double calc_latent_heat_of_vaporization(double temp);
int calc_layer_average_thermal_props(energy_bal_struct *, layer_data_struct *,
                                     soil_con_struct *, size_t, double *);
double calc_outgoing_longwave(double temp, double emis);
double calc_scale_height(double tair, double elevation);
double calc_sensible_heat(double atmos_density, double t1, double t0,
                          double Ra);
void calc_Nscale_factors(char, double *, double, double, double, double,
                         unsigned short int, double *);
double calc_rainonly(double, double, double, double);
double calc_rc(double, double, double, double, double, double, double, char);
void calc_rc_ps(char, double, double, double, double *, double, double,
                double *, double, double, double *, double, double, double,
                double *, double *);
double calc_snow_coverage(_Bool *, double, double, double, double, double,
                          double, double, double *, double, double *, double *,
                          double *);
int calc_soil_thermal_fluxes(int, double *, double *, char *, unsigned int *,
                             double *, double *, double *, double *, double *,
                             double *, double *, double *, double *, double *,
                             double *, int, int, int);
double calc_surf_energy_bal(double, double, double, double, double, double,
                            double, double, double, double, double, double,
                            double, double, double, double, double, double,
                            double, double, double, double, double, double *,
                            double *, double *, double *, double *, double *, double,
                            double *, double *, double, double *, double *, int,
                            int, size_t, size_t, double, size_t,
                            unsigned short int, int, int, unsigned short int,
                            double *, double *, atmos_data_struct *,
                            dmy_struct *, energy_bal_struct *,
                            layer_data_struct *, snow_data_struct *,
                            soil_con_struct *, veg_var_struct *);
double calc_veg_displacement(double);
double calc_veg_height(double);
double calc_veg_roughness(double);
double calc_water_balance_error(int, double, double, double);
unsigned short int calendar_from_chars(char *cal_chars);
int CalcAerodynamic(_Bool, double, double, double, double, double, double *,
                    double *, double *, double *, double *);
double CalcBlowingSnow(double, double, unsigned int, double, double, double,
                       double, double, double, double, double, double, double,
                       double, int, int, double, double, double, double *);
double CalcIcePackEnergyBalance(double Tsurf, ...);
double CalcSnowPackEnergyBalance(double Tsurf, ...);
double CalcSubFlux(double EactAir, double es, double Zrh, double AirDens,
                   double utshear, double ushear, double fe, double Tsnow,
                   double Tair, double U10, double Zo_salt, double F,
                   double *Transport);
void canopy_assimilation(char, double, double, double, double *, double, double,
                         double *, double, double, double *, double, char *,
                         double *, double *, double *, double *, double *,
                         double *, double *, double *, double *, double *);
double canopy_evap(layer_data_struct *, veg_var_struct *, _Bool,
                   unsigned short int, double *, double, double, double, double,
                   double, double, double, double, double *, double *, double *,
                   double *, double *, double *, double, double, double *);
void colavg(double *, double *, double *, double, double *, int, double,
            double);
void collect_eb_terms(energy_bal_struct, snow_data_struct, cell_data_struct,
                      int *, int *, int *, int *, int *, double, double, double,
                      int, int, double, int, int, double *, double,
                      out_data_struct *);
void collect_wb_terms(cell_data_struct, veg_var_struct, snow_data_struct,
                      double, double, double, int, double, int, double *,
                      double *, out_data_struct *);
double compute_coszen(double, double, double, unsigned short int, unsigned int);
void compute_pot_evap(double, double, double, double, double, double, double,
                      double, double, double, double *, char, double, double,
                      double, double *);
void compute_runoff_and_asat(soil_con_struct *, double *, double, double *,
                             double *);
void compute_soil_resp(int, double *, double, double, double *, double *,
                       double, double, double, double *, double *, double *);
void compute_soil_layer_thermal_properties(layer_data_struct *, double *,
                                           double *, double *, double *,
                                           double *, double *, double *,
                                           double *, size_t);
double compute_zwt(soil_con_struct *, int, double);
void correct_precip(double *, double, double, double, double);
double darkinhib(double);
int distribute_node_moisture_properties(double *, double *, double *, double *,
                                        double *, double *, double *, double *,
                                        double *, double *, double *, double *,
                                        double *, double *, double *, double *,
                                        double *, int, int, char);
void eddy(int, double, double *, double *, double, int, double, double);
void energycalc(double *, double *, int, double, double, double *, double *,
                double *);
int estimate_layer_ice_content(layer_data_struct *, double *, double *,
                               double *, double *, double *, double *, double *,
                               double, size_t, size_t, char);
int estimate_layer_ice_content_quick_flux(layer_data_struct *, double *, double,
                                          double, double, double, double *,
                                          double *, double *, double *, double,
                                          char);
double estimate_T1(double, double, double, double, double, double, double,
                   double, double, double);
void faparl(double *, double, double, double, double, double *, double *);
void fda_heat_eqn(double *, double *, int, int, ...);
void fdjac3(double *, double *, double *, double *, double *, void (*vecfunc)(
                double *, double *, int, int, ...), int);
void find_0_degree_fronts(energy_bal_struct *, double *, double *, int);
double (*funcd)(double z, double es, double Wind, double AirDens, double ZO,
                double EactAir, double F, double hsalt, double phi_r,
                double ushear,
                double Zrh);
int get_depth(lake_con_struct, double, double *);
double get_prob(double Tair, double Age, double SurfaceLiquidWater, double U10);
int get_sarea(lake_con_struct, double, double *);
void get_shear(double x, double *f, double *df, double Ur, double Zr);
double get_thresh(double Tair, double SurfaceLiquidWater, double Zo_salt);
int get_volume(lake_con_struct, double, double *);
double hiTinhib(double);
int initialize_lake(lake_var_struct *, lake_con_struct, soil_con_struct *,
                    cell_data_struct *, double, int);
int ice_melt(double, double, double *, double, snow_data_struct *,
             lake_var_struct *, double, double, double, double, double, double,
             double, double, double, double, double, double, double, double,
             double, double *, double *, double *, double *, double *, double *,
             double *, double *, double *);
void iceform(double *, double *, double, double, double *, int, double, double,
             double, double *, double *, double *, double *, double);
void icerad(double, double, double, double *, double *, double *);
int lakeice(double, double, double, double, double, double *, double, double *,
            double *, double, double);
void latent_heat_from_snow(double, double, double, double, double, double,
                           double, double *, double *, double *, double *,
                           double *);
void latsens(double, double, double, double, double, double, double, double,
             double *, double *, double);
double linear_interp(double, double, double, double, double);
double lkdrag(double, double, double, double, double);
void MassRelease(double *, double *, double *, double *);
double maximum_unfrozen_water(double, double, double, double);
double new_snow_density(double);
int newt_raph(void (*vecfunc)(double *, double *, int, int,
                              ...), double *, int);
double penman(double, double, double, double, double, double, double);
void photosynth(char, double, double, double, double, double, double, double,
                double, double, char *, double *, double *, double *, double *,
                double *);
void polint(double xa[], double ya[], int n, double x, double *y, double *dy);
void prepare_full_energy(int, all_vars_struct *, soil_con_struct *, double *,
                         double *);
int put_data(all_vars_struct *, atmos_data_struct *, soil_con_struct *,
             veg_con_struct *, veg_lib_struct *veg_lib, lake_con_struct *,
             out_data_struct *, save_data_struct *, int);
double qromb(
    double (*sub_with_height)(), double es, double Wind, double AirDens, double ZO, double EactAir, double F, double hsalt, double phi_r, double ushear, double Zrh, double a,
    double b);
void rescale_snow_energy_fluxes(double, double, snow_data_struct *,
                                energy_bal_struct *);
void rescale_snow_storage(double, double, snow_data_struct *);
void rescale_soil_veg_fluxes(double, double, cell_data_struct *,
                             veg_var_struct *);
void rhoinit(double *, double);
double rtnewt(double x1, double x2, double xacc, double Ur, double Zr);
int runoff(cell_data_struct *, energy_bal_struct *, soil_con_struct *, double,
           double *, int);
void set_node_parameters(double *, double *, double *, double *, double *,
                         double *, double *, double *, double *, double *,
                         double *, int, int);
void shear_stress(double U10, double ZO, double *ushear, double *Zo_salt,
                  double utshear);
double snow_albedo(double, double, double, double, double, int, char);
double snow_density(snow_data_struct *, double, double, double, double);
int snow_intercept(double, double, double, double, double, double, double,
                   double, double, double, double *, double *, double *,
                   double *, double *, double *, double *, double *, double *,
                   double *, double *, double *, double *, double *, double *,
                   double *, _Bool *, unsigned int *, double *, double *,
                   double *, double *, double *, double *, double *, int, int,
                   int, int, int, int, unsigned short int, double *, double *,
                   atmos_data_struct *, layer_data_struct *, soil_con_struct *,
                   veg_var_struct *);
int snow_melt(double, double, double, double, double *, double, double *,
              double, double, double, double, double, double, double, double,
              double, double, double, double, double, double *, double *,
              double *, double *, double *, double *, double *, double *,
              double *, double *, double *, double *, int, int, int, int,
              snow_data_struct *);
void soil_carbon_balance(soil_con_struct *, energy_bal_struct *,
                         cell_data_struct *, veg_var_struct *);
double soil_conductivity(double, double, double, double, double, double, double,
                         double);
int solve_lake(double, double, double, double, double, double, double, double,
               double, double, lake_var_struct *, soil_con_struct, double,
               double, dmy_struct, double);
double solve_snow(char, double, double, double, double, double, double, double,
                  double, double, double *, double *, double *, double *,
                  double *, double *, double *, double *, double *, double *,
                  double *, double *, double *, double *, double *, double *,
                  double *, double *, double *, double *, double *, double *,
                  double *, double *, double *, double *, double *, double *,
                  int, size_t, unsigned short int, unsigned short int, double,
                  size_t, size_t, int, int *, double *, double *, dmy_struct *,
                  atmos_data_struct *, energy_bal_struct *, layer_data_struct *,
                  snow_data_struct *, soil_con_struct *, veg_var_struct *);
double solve_surf_energy_bal(double Tsurf, ...);
int solve_T_profile(double *, double *, char *, unsigned int *, double *,
                    double *, double *, double *, double, double *, double *,
                    double *, double *, double *, double *, double *, double,
                    int, int *, int, int, int);
int solve_T_profile_implicit(double *, double *, char *, unsigned int *,
                             double *, double *, double *, double *, double,
                             double *, double *, double *, double *, double *,
                             double *, double *, double, int, int *, int, int,
                             double *, double *, double *, double *, double *,
                             double *, double *);
double specheat(double);
double StabilityCorrection(double, double, double, double, double, double);
double sub_with_height(double z, double es, double Wind, double AirDens,
                       double ZO, double EactAir, double F, double hsalt,
                       double phi_r, double ushear, double Zrh);
int surface_fluxes(_Bool, double, double, double, double, double *, double *,
                   double *, double *, double *, double *, double *, double *,
                   double *, double *, double *, double *, double *, size_t,
                   size_t, unsigned short int, double, unsigned short int,
                   size_t, unsigned short int, atmos_data_struct *,
                   dmy_struct *, energy_bal_struct *, global_param_struct *,
                   cell_data_struct *, snow_data_struct *, soil_con_struct *,
                   veg_var_struct *, double, double, double, double *);
double svp(double);
double svp_slope(double);
void temp_area(double, double, double, double *, double *, double *, double *,
               double, double *, int, double, double, double *, double *,
               double *);
void tracer_mixer(double *, int *, double *, int, double, double, double *);
void transpiration(layer_data_struct *, veg_var_struct *, unsigned short int,
                   double, double, double, double, double, double, double,
                   double, double *, double *, double *, double *, double *,
                   double *, double, double, double *);
double transport_with_height(double z, double es, double Wind, double AirDens,
                             double ZO, double EactAir, double F, double hsalt,
                             double phi_r, double ushear, double Zrh);
double trapzd(
    double (*funcd)(), double es, double Wind, double AirDens, double ZO, double EactAir, double F, double hsalt, double phi_r, double ushear, double Zrh, double a, double b,
    int n);
void tridia(int, double *, double *, double *, double *, double *);
void tridiag(double *, double *, double *, double *, unsigned int);
int vic_run(int, atmos_data_struct *, all_vars_struct *, dmy_struct *,
            global_param_struct *, lake_con_struct *, soil_con_struct *,
            veg_con_struct *, veg_lib_struct *, veg_hist_struct *veg_hist);
double volumetric_heat_capacity(double, double, double, double);
int water_balance(lake_var_struct *, lake_con_struct, double, all_vars_struct *,
                  int, int, int, double, soil_con_struct, veg_con_struct);
int water_energy_balance(int, double *, double *, double, double, double,
                         double, double, double, double, double, double, double,
                         double, double, double, double *, double *, double *,
                         double *, double *, double *, double *, double,
                         double *, double *, double *, double *, double *,
                         double);
int water_under_ice(int, double, double, double *, double *, double, int,
                    double, double, double, double *, double *, double *,
                    double *, int, double, double, double, double *);
void wrap_compute_zwt(soil_con_struct *, cell_data_struct *);
void write_layer(layer_data_struct *, int, double *);
void write_vegvar(veg_var_struct *, int);
void zero_output_list(out_data_struct *);
typedef struct {
    size_t N_ELEM;
    _Bool SIGNED;
    _Bool SUPPLIED;
    double multiplier;
    char varname[2048];
} force_type_struct;
typedef struct {
    force_type_struct TYPE[N_FORCING_TYPES];
    double FORCE_DT[2];
    size_t force_steps_per_day[2];
    unsigned short int FORCE_ENDIAN[2];
    int FORCE_FORMAT[2];
    int FORCE_INDEX[2][N_FORCING_TYPES];
    size_t N_TYPES[2];
} param_set_struct;
double all_30_day_from_dmy(dmy_struct *dmy);
double all_leap_from_dmy(dmy_struct *dmy);
void calc_root_fractions(veg_con_struct *veg_con, soil_con_struct *soil_con);
void compute_treeline(atmos_data_struct *, dmy_struct *, double, double *,
                      _Bool *);
void cmd_proc(int argc, char **argv, char *globalfilename);
void compress_files(char string[]);
out_data_struct *create_output_list(void);
double date2num(double origin, dmy_struct *date, double tzoffset,
                unsigned short int calendar, unsigned short int time_units);
void dmy_all_30_day(double julian, dmy_struct *dmy);
void dmy_all_leap(double julian, dmy_struct *dmy);
void dmy_julian_day(double julian, unsigned short int calendar,
                    dmy_struct *dmy);
void dmy_no_leap_day(double julian, dmy_struct *dmy);
void dt_seconds_to_time_units(unsigned short int time_units, double dt_seconds,
                              double *dt_time_units);
void display_current_settings(int);
double fractional_day_from_dmy(dmy_struct *dmy);
void free_all_vars(all_vars_struct *all_vars, int Nveg);
void free_dmy(dmy_struct **dmy);
void free_out_data_files(out_data_file_struct **);
void free_out_data(out_data_struct **);
void free_vegcon(veg_con_struct **veg_con);
double get_dist(double lat1, double long1, double lat2, double long2);
void get_parameters(FILE *paramfile);
void init_output_list(out_data_struct *out_data, int write, char *format,
                      int type, double mult);
void initialize_filenames(void);
void initialize_fileps(void);
void initialize_global(void);
void initialize_options(void);
void initialize_parameters(void);
void initialize_snow(snow_data_struct **snow, size_t veg_num);
void initialize_soil(cell_data_struct **cell, soil_con_struct *soil_con,
                     size_t veg_num);
void initialize_time(void);
void initialize_veg(veg_var_struct **veg_var, size_t nveg);
double julian_day_from_dmy(dmy_struct *dmy, unsigned short int calendar);
_Bool leap_year(unsigned short int year, unsigned short int calendar);
all_vars_struct make_all_vars(size_t nveg);
cell_data_struct **make_cell_data(size_t veg_type_num);
dmy_struct *make_dmy(global_param_struct *global);
energy_bal_struct **make_energy_bal(size_t nveg);
void make_lastday(unsigned short int calendar, unsigned short int year,
                  unsigned short int lastday[]);
snow_data_struct **make_snow_data(size_t nveg);
veg_var_struct **make_veg_var(size_t veg_type_num);
double no_leap_day_from_dmy(dmy_struct *dmy);
void num2date(double origin, double time_value, double tzoffset,
              unsigned short int calendar, unsigned short int time_units,
              dmy_struct *date);
FILE *open_file(char string[], char type[]);
void parse_nc_time_units(char *nc_unit_chars, unsigned short int *units,
                         dmy_struct *dmy);
void print_cell_data(cell_data_struct *cell, size_t nlayers, size_t nfrost);
void print_dmy(dmy_struct *dmy);
void print_energy_bal(energy_bal_struct *eb, size_t nnodes, size_t nfronts);
void print_filenames(filenames_struct *fnames);
void print_filep(filep_struct *fp);
void print_force_type(force_type_struct *force_type);
void print_global_param(global_param_struct *gp);
void print_lake_con(lake_con_struct *lcon, size_t nlnodes);
void print_lake_var(lake_var_struct *lvar, size_t nlnodes, size_t nfronts,
                    size_t nlayers, size_t nnodes, size_t nfrost);
void print_layer_data(layer_data_struct *ldata, size_t nfrost);
void print_license(void);
void print_option(option_struct *option);
void print_out_data(out_data_struct *out, size_t nelem);
void print_out_data_file(out_data_file_struct *outf);
void print_param_set(param_set_struct *param_set);
void print_parameters(parameters_struct *param);
void print_save_data(save_data_struct *save);
void print_snow_data(snow_data_struct *snow);
void print_soil_con(soil_con_struct *scon, size_t nlayers, size_t nnodes,
                    size_t nfrost, size_t nbands, size_t nzwt);
void print_veg_con(veg_con_struct *vcon, size_t nroots, char blowing, char lake,
                   char carbon, size_t ncanopy);
void print_veg_lib(veg_lib_struct *vlib, char carbon);
void print_veg_var(veg_var_struct *vvar, size_t ncanopy);
void print_version(char *);
void print_usage(char *);
int set_output_var(out_data_file_struct *, int, int, out_data_struct *, char *,
                   int, char *, int, double);
void soil_moisture_from_water_table(soil_con_struct *soil_con, size_t nlayers);
unsigned short int timeunits_from_chars(char *units_chars);
int update_step_vars(all_vars_struct *, veg_con_struct *, veg_hist_struct *);
int valid_date(unsigned short int calendar, dmy_struct *dmy);
void validate_parameters(void);
int flag;
size_t NR;
size_t NF;
global_param_struct global_param;
option_struct options;
parameters_struct param;
filenames_struct filenames;
filep_struct filep;
param_set_struct param_set;
'''
