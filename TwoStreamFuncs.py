
class TwoStreamFuncs:
    
    def __init__(self):


        return

    def generate_ice_physical_params(density,grain_rds,layer_type,dz,algae,solzen, incoming_i, DIRECT):

        """
        takes user-defined params and collects them into a named tuple for passing to call_snicar()
        
        """

        import collections

        params = collections.namedtuple("params","rho_layers, grain_rds, layer_type, dz, mss_cnc_glacier_algae, solzen")
        params.grain_rds = grain_rds
        params.rho_layers = density
        params.layer_type = layer_type
        params.dz = dz
        params.mss_cnc_glacier_algae = algae
        params.solzen = solzen
        params.incoming_i = incoming_i
        params.DIRECT = DIRECT

        return params


    def generate_incoming_irradiance(params):

        import xarray as xr

        incoming_i = params.incoming_i
        DIRECT = params.DIRECT
        solzen = params.solzen 

        dir_fsds = '/home/joe/Code/BioSNICAR_GO_PY/Data/Mie_files/480band/fsds/'
        zen = str('SZA'+str(solzen).rjust(2,'0'))
        
        if DIRECT:
            if incoming_i == 0:
                Incoming_file = xr.open_dataset(str(dir_fsds + "swnb_480bnd_mlw_clr_"+zen+".nc")) 
                
            elif incoming_i == 1:
                Incoming_file = xr.open_dataset(str(dir_fsds + "swnb_480bnd_mls_clr_"+zen+".nc"))
                
            elif incoming_i == 2:
                Incoming_file = xr.open_dataset(str(dir_fsds + "swnb_480bnd_saw_clr_"+zen+".nc"))
                
            elif incoming_i == 3:
                Incoming_file = xr.open_dataset(str(dir_fsds + "swnb_480bnd_sas_clr_"+zen+".nc"))
                
            elif incoming_i == 4:
                Incoming_file = xr.open_dataset(str(dir_fsds + "swnb_480bnd_smm_clr_"+zen+".nc"))
                
            elif incoming_i == 5:
                Incoming_file = xr.open_dataset(str(dir_fsds + "swnb_480bnd_hmn_clr_"+zen+".nc"))
                
            elif incoming_i == 6:
                Incoming_file = xr.open_dataset(str(dir_fsds + "swnb_480bnd_toa_clr.nc"))

            else:
                raise ValueError ("Invalid choice of atmospheric profile")


        else:

            if incoming_i == 0:
                Incoming_file = xr.open_dataset(str(dir_fsds + "swnb_480bnd_mlw_cld.nc"))
            elif incoming_i == 1:
                Incoming_file = xr.open_dataset(str(dir_fsds + "swnb_480bnd_mls_cld.nc"))
            elif incoming_i == 2:
                Incoming_file = xr.open_dataset(str(dir_fsds + "swnb_480bnd_saw_cld.nc"))
            elif incoming_i == 3:
                Incoming_file = xr.open_dataset(str(dir_fsds + "swnb_480bnd_sas_cld.nc"))
            elif incoming_i == 4:
                Incoming_file = xr.open_dataset(str(dir_fsds + "swnb_480bnd_smm_cld.nc"))
            elif incoming_i == 5:
                Incoming_file = xr.open_dataset(str(dir_fsds + "swnb_480bnd_hmn_cld.nc"))   
            elif incoming_i == 6:
                Incoming_file = xr.open_dataset(str(dir_fsds + "swnb_480bnd_toa_cld.nc"))

            else:
                raise ValueError ("Invalid choice of atmospheric profile") 


        #flx_dwn_sfc is the spectral irradiance in W m-2 and is pre-calculated (flx_frc_sfc*flx_bb_sfc in original code)
        incoming = Incoming_file['flx_dwn_sfc'].values 
        incoming[incoming<=0]=1e-30
        incoming = incoming[10:]

        return incoming


    def call_snicar(params):

        from SNICAR_feeder import snicar_feeder

    # set dir_base to the location of the BioSNICAR_GO_PY folder

        dir_base = '/home/joe/Code/BioSNICAR_GO_PY/'
        savepath = dir_base # base path for saving figures
        
        TOON = False # toggle Toon et al tridiagonal matrix solver
        ADD_DOUBLE = True # toggle adding-doubling solver
        layer_type = params.layer_type
        DIRECT   = 1        # 1= Direct-beam incident flux, 0= Diffuse incident flux
        APRX_TYP = 1        # 1= Eddington, 2= Quadrature, 3= Hemispheric Mean
        DELTA    = 1        # 1= Apply Delta approximation, 0= No delta
        solzen   = params.solzen      # if DIRECT give solar zenith angle (degrees from 0 = nadir, 90 = horizon)
        rf_ice = 2        # define source of ice refractive index data. 0 = Warren 1984, 1 = Warren 2008, 2 = Picard 2016
        incoming_i = 4
        nbr_lyr = len(params.dz)  # number of snow layers
        R_sfc = 0.15 # reflectance of underlying surface - set across all wavelengths
        rwater = [0]*len(params.dz) # if  using Mie calculations, add radius of optional liquid water coating
        grain_shp =[0]*len(params.dz) # grain shape(He et al. 2016, 2017)
        shp_fctr = [0]*len(params.dz) # shape factor (ratio of aspherical grain radii to that of equal-volume sphere)
        grain_ar = [0]*len(params.dz) # aspect ratio (ratio of width to length)
        side_length = 0
        depth=0
        grain_rds = params.grain_rds
        rho_layers = params.rho_layers
        dz = params.dz
        
        mss_cnc_soot1 = [0]*len(params.dz)    # uncoated black carbon (Bohren and Huffman, 1983)
        mss_cnc_soot2 = [0]*len(params.dz)    # coated black carbon (Bohren and Huffman, 1983)
        mss_cnc_brwnC1 = [0]*len(params.dz)   # uncoated brown carbon (Kirchstetter et al. (2004).)
        mss_cnc_brwnC2 = [0]*len(params.dz)   # sulfate-coated brown carbon (Kirchstetter et al. (2004).)
        mss_cnc_dust1 = [0]*len(params.dz)    # dust size 1 (r=0.05-0.5um) (Balkanski et al 2007)
        mss_cnc_dust2 = [0]*len(params.dz)    # dust size 2 (r=0.5-1.25um) (Balkanski et al 2007)
        mss_cnc_dust3 = [0]*len(params.dz)    # dust size 3 (r=1.25-2.5um) (Balkanski et al 2007)
        mss_cnc_dust4 = [0]*len(params.dz)    # dust size 4 (r=2.5-5.0um)  (Balkanski et al 2007)
        mss_cnc_dust5 = [0]*len(params.dz)    # dust size 5 (r=5.0-50um)  (Balkanski et al 2007)
        mss_cnc_ash1 = [0]*len(params.dz)    # volcanic ash size 1 (r=0.05-0.5um) (Flanner et al 2014)
        mss_cnc_ash2 = [0]*len(params.dz)    # volcanic ash size 2 (r=0.5-1.25um) (Flanner et al 2014)
        mss_cnc_ash3 = [0]*len(params.dz)    # volcanic ash size 3 (r=1.25-2.5um) (Flanner et al 2014)
        mss_cnc_ash4 = [0]*len(params.dz)    # volcanic ash size 4 (r=2.5-5.0um) (Flanner et al 2014)
        mss_cnc_ash5 = [0]*len(params.dz)    # volcanic ash size 5 (r=5.0-50um) (Flanner et al 2014)
        mss_cnc_ash_st_helens = [0]*len(params.dz)   # ash from Mount Saint Helen's
        mss_cnc_Skiles_dust1 = [0]*len(params.dz)    # Colorado dust size 1 (Skiles et al 2017)
        mss_cnc_Skiles_dust2 = [0]*len(params.dz)    # Colorado dust size 2 (Skiles et al 2017)
        mss_cnc_Skiles_dust3 = [0]*len(params.dz)    # Colorado dust size 3 (Skiles et al 2017)
        mss_cnc_Skiles_dust4 = [0]*len(params.dz)  # Colorado dust size 4 (Skiles et al 2017)
        mss_cnc_Skiles_dust5 = [0]*len(params.dz)  # Colorado dust size 5 (Skiles et al 2017)
        mss_cnc_GreenlandCentral1 = [0]*len(params.dz) # Greenland Central dust size 1 (Polashenski et al 2015)
        mss_cnc_GreenlandCentral2 = [0]*len(params.dz) # Greenland Central dust size 2 (Polashenski et al 2015)
        mss_cnc_GreenlandCentral3 = [0]*len(params.dz) # Greenland Central dust size 3 (Polashenski et al 2015)
        mss_cnc_GreenlandCentral4 = [0]*len(params.dz) # Greenland Central dust size 4 (Polashenski et al 2015)
        mss_cnc_GreenlandCentral5 = [0]*len(params.dz) # Greenland Central dust size 5 (Polashenski et al 2015)
        mss_cnc_Cook_Greenland_dust_L = [0]*len(params.dz)
        mss_cnc_Cook_Greenland_dust_C = [0]*len(params.dz)
        mss_cnc_Cook_Greenland_dust_H = [0]*len(params.dz)
        mss_cnc_snw_alg = [0]*len(params.dz)    # Snow Algae (spherical, C nivalis) (Cook et al. 2017)
        mss_cnc_glacier_algae = params.mss_cnc_glacier_algae   # glacier algae type1 (Cook et al. 2020)

        nbr_aer = 30

        # Set names of files containing the optical properties of these LAPs:
        FILE_soot1  = 'mie_sot_ChC90_dns_1317.nc'
        FILE_soot2  = 'miecot_slfsot_ChC90_dns_1317.nc'
        FILE_brwnC1 = 'brC_Kirch_BCsd.nc'
        FILE_brwnC2 = 'brC_Kirch_BCsd_slfcot.nc'
        FILE_dust1  = 'dust_balkanski_central_size1.nc'
        FILE_dust2  = 'dust_balkanski_central_size2.nc'
        FILE_dust3  = 'dust_balkanski_central_size3.nc'
        FILE_dust4  = 'dust_balkanski_central_size4.nc'
        FILE_dust5 = 'dust_balkanski_central_size5.nc'
        FILE_ash1  = 'volc_ash_eyja_central_size1.nc'
        FILE_ash2 = 'volc_ash_eyja_central_size2.nc'
        FILE_ash3 = 'volc_ash_eyja_central_size3.nc'
        FILE_ash4 = 'volc_ash_eyja_central_size4.nc'
        FILE_ash5 = 'volc_ash_eyja_central_size5.nc'
        FILE_ash_st_helens = 'volc_ash_mtsthelens_20081011.nc'
        FILE_Skiles_dust1 = 'dust_skiles_size1.nc'
        FILE_Skiles_dust2 = 'dust_skiles_size2.nc'
        FILE_Skiles_dust3 = 'dust_skiles_size3.nc'
        FILE_Skiles_dust4 = 'dust_skiles_size4.nc'
        FILE_Skiles_dust5 = 'dust_skiles_size5.nc'
        FILE_GreenlandCentral1 = 'dust_greenland_central_size1.nc'
        FILE_GreenlandCentral2 = 'dust_greenland_central_size2.nc'
        FILE_GreenlandCentral3 = 'dust_greenland_central_size3.nc'
        FILE_GreenlandCentral4 = 'dust_greenland_central_size4.nc'
        FILE_GreenlandCentral5  = 'dust_greenland_central_size5.nc'
        FILE_Cook_Greenland_dust_L = 'dust_greenland_Cook_LOW_20190911.nc'
        FILE_Cook_Greenland_dust_C = 'dust_greenland_Cook_CENTRAL_20190911.nc'
        FILE_Cook_Greenland_dust_H = 'dust_greenland_Cook_HIGH_20190911.nc'
        FILE_snw_alg  = 'snw_alg_r025um_chla020_chlb025_cara150_carb140.nc'
        FILE_glacier_algae = 'Glacier_Algae_480.nc'


        #######################################
        # IF NO INPUT ERRORS --> FUNCTION CALLS
        #######################################

            
        [wvl, albedo, BBA, BBAVIS, BBANIR, abs_slr, heat_rt, F_btm_net] =\
        snicar_feeder(dir_base,\
        rf_ice, incoming_i, DIRECT, layer_type,\
        APRX_TYP, DELTA, solzen, TOON, ADD_DOUBLE, R_sfc, dz, rho_layers, grain_rds,\
        side_length, depth, rwater, nbr_lyr, nbr_aer, grain_shp, shp_fctr, grain_ar,\
        mss_cnc_soot1, mss_cnc_soot2, mss_cnc_brwnC1, mss_cnc_brwnC2, mss_cnc_dust1,\
        mss_cnc_dust2, mss_cnc_dust3, mss_cnc_dust4, mss_cnc_dust5, mss_cnc_ash1, mss_cnc_ash2,\
        mss_cnc_ash3, mss_cnc_ash4, mss_cnc_ash5, mss_cnc_ash_st_helens, mss_cnc_Skiles_dust1, mss_cnc_Skiles_dust2,\
        mss_cnc_Skiles_dust3, mss_cnc_Skiles_dust4, mss_cnc_Skiles_dust5, mss_cnc_GreenlandCentral1,\
        mss_cnc_GreenlandCentral2, mss_cnc_GreenlandCentral3, mss_cnc_GreenlandCentral4,\
        mss_cnc_GreenlandCentral5, mss_cnc_Cook_Greenland_dust_L, mss_cnc_Cook_Greenland_dust_C,\
        mss_cnc_Cook_Greenland_dust_H, mss_cnc_snw_alg, mss_cnc_glacier_algae, FILE_soot1,\
        FILE_soot2, FILE_brwnC1, FILE_brwnC2, FILE_dust1, FILE_dust2, FILE_dust3, FILE_dust4, FILE_dust5,\
        FILE_ash1, FILE_ash2, FILE_ash3, FILE_ash4, FILE_ash5, FILE_ash_st_helens, FILE_Skiles_dust1, FILE_Skiles_dust2,\
        FILE_Skiles_dust3, FILE_Skiles_dust4, FILE_Skiles_dust5, FILE_GreenlandCentral1,\
        FILE_GreenlandCentral2, FILE_GreenlandCentral3, FILE_GreenlandCentral4, FILE_GreenlandCentral5,\
        FILE_Cook_Greenland_dust_L, FILE_Cook_Greenland_dust_C, FILE_Cook_Greenland_dust_H, FILE_snw_alg, FILE_glacier_algae)

        return albedo, BBA, F_btm_net
