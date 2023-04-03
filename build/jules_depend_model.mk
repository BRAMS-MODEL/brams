cellarea_calc_mod.o : $(JULES_01)/cellarea_calc_mod.F90 conversions_mod_jls.o parkind1.o \
	planet_constants_mod_jls.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

clim_calc.o : $(JULES_01)/clim_calc.F90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

delta_temp.o : $(JULES_01)/delta_temp.F90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

day_calc.o : $(JULES_01)/day_calc.F90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

diff_atmos_ch4.o : $(JULES_01)/diff_atmos_ch4.F90 imogen_constants.o imogen_progs.o imogen_run.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

diffcarb_land_ch4.o : $(JULES_01)/diffcarb_land_ch4.F90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

diffcarb_land_co2.o : $(JULES_01)/diffcarb_land_co2.F90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

drdat.o : $(JULES_01)/drdat.F90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

gcm_anlg.o : $(JULES_01)/gcm_anlg.F90 imogen_map.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

imogen_check.o : $(JULES_01)/imogen_check.F90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

imogen_update_clim.o : $(JULES_01)/imogen_update_clim.F90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

imogen_update_carb.o : $(JULES_01)/imogen_update_carb.F90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

invert.o : $(JULES_01)/invert.F90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

ocean_co2.o : $(JULES_01)/ocean_co2.F90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

radf_ch4.o : $(JULES_01)/radf_ch4.F90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

radf_co2.o : $(JULES_01)/radf_co2.F90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

radf_non_co2.o : $(JULES_01)/radf_non_co2.F90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

response.o : $(JULES_01)/response.F90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

redis.o : $(JULES_01)/redis.F90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

solang.o : $(JULES_01)/solang.F90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

rndm.o : $(JULES_01)/rndm.F90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

sunny.o : $(JULES_01)/sunny.F90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

solpos.o : $(JULES_01)/solpos.F90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

imogen_clim.o : $(JULES_01)/var/imogen_clim.F90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

imogen_drive_vars.o : $(JULES_01)/var/imogen_drive_vars.F90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

imogen_io_vars.o : $(JULES_01)/var/imogen_io_vars.F90 imogen_constants.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

imogen_map.o : $(JULES_01)/var/imogen_map.F90 ancil_info.o imogen_constants.o io_constants.o \
	jules_fields_mod.o model_grid_mod.o theta_field_sizes_mod.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

imogen_progs.o : $(JULES_01)/var/imogen_progs.F90 imogen_constants.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

CN_utils_mod.o : $(JULES_02)/CN_utils_mod.F90 crop_utils_mod.o ereport_mod.o jules_surface_types_mod.o \
	jules_vegetation_mod.o pftparm_mod.o trif_mod.o um_types.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

aero.o : $(JULES_02)/aero.F90 parkind1.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

ancil_info.o : $(JULES_02)/ancil_info.F90 parkind1.o um_types.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

bvoc_vars.o : $(JULES_02)/bvoc_vars.F90 parkind1.o um_types.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

cable_prognostic_info_mod.o : $(JULES_02)/cable_prognostic_info_mod.F90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

c_elevate.o : $(JULES_02)/c_elevate.F90 jules_print_mgr.o max_dimensions.o missing_data_mod.o um_types.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

cable_surface_types_mod.o : $(JULES_02)/cable_surface_types_mod.F90 ereport_mod.o jules_print_mgr.o \
	 max_dimensions.o missing_data_mod.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

calc_c_comps_triffid_mod.o : $(JULES_02)/calc_c_comps_triffid_mod.F90 jules_surface_mod.o \
	jules_vegetation_mod.o pftparm_mod.o um_types.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

calc_litter_flux_mod.o : $(JULES_02)/calc_litter_flux_mod.F90 um_types.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

coastal.o : $(JULES_02)/coastal.F90 parkind1.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

crop_utils_mod.o : $(JULES_02)/crop_utils_mod.F90 cropparm_mod.o um_types.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

crop_vars_mod.o : $(JULES_02)/crop_vars_mod.F90 parkind1.o um_types.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

datetime_utils_mod.o : $(JULES_02)/datetime_utils_mod.F90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

fire_vars_mod.o : $(JULES_02)/fire_vars_mod.F90 parkind1.o um_types.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

fluxes.o : $(JULES_02)/fluxes.F90 jules_vegetation_mod.o parkind1.o um_types.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

jules_deposition_mod.o : $(JULES_02)/jules_deposition_mod.F90 ereport_mod.o jules_print_mgr.o \
	jules_surface_mod.o max_dimensions.o missing_data_mod.o parkind1.o um_types.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

jules_hydrology_mod.o : $(JULES_02)/jules_hydrology_mod.F90 ereport_mod.o jules_print_mgr.o um_types.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

jules_irrig_mod.o : $(JULES_02)/jules_irrig_mod.F90 ereport_mod.o errormessagelength_mod.o io_constants.o \
	jules_hydrology_mod.o jules_print_mgr.o jules_surface_types_mod.o jules_vegetation_mod.o logging_mod.o \
	max_dimensions.o missing_data_mod.o parkind1.o string_utils_mod.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

jules_mod.o : $(JULES_02)/jules_mod.F90 atm_fields_bounds_mod.o parkind1.o um_types.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

jules_model_environment_mod.o : $(JULES_02)/jules_model_environment_mod.F90 ereport_mod.o \
	errormessagelength_mod.o io_constants.o jules_print_mgr.o logging_mod.o missing_data_mod.o \
	string_utils_mod.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

jules_plant_n_uptake_mod.o : $(JULES_02)/jules_plant_n_uptake_mod.F90 ereport_mod.o jules_surface_mod.o \
	jules_vegetation_mod.o string_utils_mod.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

jules_print_mgr.o : $(JULES_02)/jules_print_mgr.F90 logging_mod.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

jules_radiation_mod.o : $(JULES_02)/jules_radiation_mod.F90 ancil_info.o ereport_mod.o jules_print_mgr.o \
	missing_data_mod.o um_types.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

jules_rivers_mod.o : $(JULES_02)/jules_rivers_mod.F90 ereport_mod.o jules_irrig_mod.o jules_print_mgr.o \
	missing_data_mod.o parkind1.o um_types.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

jules_science_fixes_mod.o : $(JULES_02)/jules_science_fixes_mod.F90 ereport_mod.o errormessagelength_mod.o \
	io_constants.o jules_print_mgr.o logging_mod.o string_utils_mod.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

jules_snow_mod.o : $(JULES_02)/jules_snow_mod.F90 ereport_mod.o jules_print_mgr.o jules_surface_mod.o \
	jules_surface_types_mod.o jules_vegetation_mod.o max_dimensions.o missing_data_mod.o um_types.o \
	water_constants_mod_jls.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

jules_sea_seaice_mod.o : $(JULES_02)/jules_sea_seaice_mod.F90 c_kappai_mod.o ereport_mod.o jules_print_mgr.o \
	missing_data_mod.o um_types.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

jules_soil_biogeochem_mod.o : $(JULES_02)/jules_soil_biogeochem_mod.F90 ereport_mod.o jules_print_mgr.o \
	jules_surface_mod.o jules_vegetation_mod.o um_types.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

jules_soil_ecosse_mod.o : $(JULES_02)/jules_soil_ecosse_mod.F90 ancil_info.o conversions_mod_jls.o \
	ecosse_param_mod.o ereport_mod.o jules_plant_n_uptake_mod.o jules_soil_mod.o jules_vegetation_mod.o \
	max_dimensions.o timestep_mod.o um_types.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

jules_soil_mod.o : $(JULES_02)/jules_soil_mod.F90 ereport_mod.o jules_irrig_mod.o jules_print_mgr.o \
	jules_surface_mod.o max_dimensions.o um_types.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

jules_surface_mod.o : $(JULES_02)/jules_surface_mod.F90 ereport_mod.o jules_print_mgr.o \
	jules_surface_types_mod.o missing_data_mod.o planet_constants_mod_jls.o um_types.o \
	water_constants_mod_jls.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

jules_surface_types_mod.o : $(JULES_02)/jules_surface_types_mod.F90 ereport_mod.o \
	errormessagelength_mod.o io_constants.o jules_print_mgr.o logging_mod.o max_dimensions.o \
	missing_data_mod.o string_utils_mod.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

jules_vegetation_mod.o : $(JULES_02)/jules_vegetation_mod.F90 conversions_mod_jls.o ereport_mod.o \
	jules_hydrology_mod.o jules_print_mgr.o jules_surface_mod.o jules_surface_types_mod.o max_dimensions.o \
	missing_data_mod.o timestep_mod.o um_types.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

jules_water_resources_mod.o : $(JULES_02)/jules_water_resources_mod.F90 ereport_mod.o jules_irrig_mod.o \
	jules_rivers_mod.o jules_soil_mod.o missing_data_mod.o um_types.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

lake_mod.o : $(JULES_02)/lake_mod.F90 atm_fields_bounds_mod.o parkind1.o um_types.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

land_tile_ids.o : $(JULES_02)/land_tile_ids.F90 ereport_mod.o errormessagelength_mod.o jules_print_mgr.o \
	jules_snow_mod.o jules_surface_mod.o jules_surface_types_mod.o max_dimensions.o missing_data_mod.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

next_gen_biogeochem_mod.o : $(JULES_02)/next_gen_biogeochem_mod.F90 conversions_mod_jls.o veg3_field_mod.o \
	veg3_litter_mod.o veg3_param_mod.o veg3_red_dynamic_mod.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

max_dimensions.o : $(JULES_02)/max_dimensions.F90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

ozone_vars.o : $(JULES_02)/ozone_vars.F90 parkind1.o um_types.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

overbank_inundation_mod.o : $(JULES_02)/overbank_inundation_mod.F90 ereport_mod.o jules_print_mgr.o \
	jules_rivers_mod.o missing_data_mod.o um_types.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

p_s_parms.o : $(JULES_02)/p_s_parms.F90 parkind1.o um_types.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

prognostics.o : $(JULES_02)/prognostics.F90 parkind1.o um_types.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

soil_biogeochem_control_mod.o : $(JULES_02)/soil_biogeochem_control_mod.F90 ancil_info.o ecosse_control_mod.o \
	ereport_mod.o jules_hydrology_mod.o jules_soil_biogeochem_mod.o jules_soil_ecosse_mod.o jules_soil_mod.o \
		jules_surface_types_mod.o model_time_mod.o parkind1.o soil_ecosse_vars_mod.o timestep_mod.o um_types.o \
		yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

river_control.o : $(JULES_02)/river_control.F90 ancil_info.o atm_fields_bounds_mod.o diag_swchs.o ereport_mod.o \
	jules_irrig_mod.o jules_rivers_mod.o jules_soil_mod.o jules_surface_types_mod.o model_grid_mod.o \
	overbank_inundation_mod.o overbank_update_mod.o parallel_mod.o parkind1.o rivers_route_mod.o timestep_mod.o \
	um_types.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

surf_couple_explicit_mod.o : $(JULES_02)/surf_couple_explicit_mod.F90 aero.o ancil_info.o atm_fields_bounds_mod.o \
bvoc_vars.o coastal.o crop_vars_mod.o dust_param_mod.o ereport_mod.o fluxes.o jules_griddiag_sf_explicit_jls.o \
jules_gridinit_sf_explicit_jls.o jules_land_sf_explicit_jls.o jules_mod.o jules_model_environment_mod.o \
jules_print_mgr.o jules_sea_seaice_mod.o jules_soil_mod.o jules_ssi_sf_explicit_jls.o jules_surface_types_mod.o \
jules_vegetation_mod.o lake_mod.o metstats_mod.o ozone_vars.o p_s_parms.o parkind1.o prognostics.o sf_diags_mod.o \
switches.o trif_vars_mod.o trifctl.o um_types.o urban_param_mod.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

jules_griddiag_sf_explicit_jls.o : $(JULES_26)/jules_griddiag_sf_explicit_jls.F90 sf_aero.o switches.o \
	sf_orog_gb_jls.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

jules_gridinit_sf_explicit_jls.o : $(JULES_26)/jules_gridinit_sf_explicit_jls.F90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

soil_ecosse_vars_mod.o : $(JULES_02)/soil_ecosse_vars_mod.F90 parkind1.o um_types.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

surf_couple_extra_mod.o : $(JULES_02)/surf_couple_extra_mod.F90 ancil_info.o atm_fields_bounds_mod.o \
	conversions_mod_jls.o crop_mod.o crop_vars_mod.o ereport_mod.o fao_evapotranspiration_mod.o fire_mod.o \
	fire_timestep_mod.o fire_vars_mod.o fluxes.o forcing.o gridbox_mean_mod.o hydrol_jls_mod.o inferno_io_mod.o \
	inferno_mod.o irrigation_mod.o jules_hydrology_mod.o jules_irrig_mod.o jules_mod.o jules_model_environment_mod.o \
	jules_print_mgr.o jules_rivers_mod.o jules_snow_mod.o jules_soil_biogeochem_mod.o jules_soil_mod.o \
	jules_surface_mod.o jules_surface_types_mod.o jules_vegetation_mod.o jules_water_resources_mod.o \
	lake_mod.o metstats_mod.o metstats_timestep.o model_grid_mod.o model_time_mod.o next_gen_biogeochem_mod.o \
	p_s_parms.o parkind1.o prognostics.o river_control.o sf_diags_mod.o snow_mod.o soil_ecosse_vars_mod.o \
	theta_field_sizes_mod.o top_pdm.o trif_vars_mod.o trifctl.o um_types.o urban_param_mod.o veg3_field_mod.o \
	veg3_param_mod.o veg_control.o water_constants_mod_jls.o water_resources_control_mod.o yomhook.o zenith_mod.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

surf_couple_implicit_mod.o : $(JULES_02)/surf_couple_implicit_mod.F90 aero.o ancil_info.o atm_fields_bounds_mod.o \
	coastal.o crop_vars_mod.o ereport_mod.o fluxes.o jules_grid_update_implicit_jls.o \
	jules_griddiag_sf_implicit_jls.o jules_land_sf_implicit.jls.o jules_mod.o jules_model_environment_mod.o \
	jules_print_mgr.o jules_sea_seaice_mod.o jules_soil_mod.o jules_ssi_sf_implicit.jls.o jules_surface_types_mod.o \
	lake_mod.o parkind1.o prognostics.o sf_diags_mod.o um_types.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

jules_griddiag_sf_implicit_jls.o : $(JULES_26)/jules_griddiag_sf_implicit_jls.F90
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

jules_grid_update_implicit_jls.o : $(JULES_26)/jules_grid_update_implicit_jls.F90 sf_melt_jls.o screen_tq_jls.o \
	im_sf_pt2_jls.o sice_htf_jls.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

surf_couple_radiation_mod.o : $(JULES_02)/surf_couple_radiation_mod.F90 ancil_info.o coastal.o ereport_mod.o \
	jules_land_albedo_jls_mod.o jules_mod.o jules_model_environment_mod.o jules_print_mgr.o jules_sea_seaice_mod.o \
	jules_ssi_albedo_jls_mod.o lake_mod.o p_s_parms.o parkind1.o prognostics.o theta_field_sizes_mod.o um_types.o \
	urban_param_mod.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

switches.o : $(JULES_02)/switches.F90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

switches_urban.o : $(JULES_02)/switches_urban.F90 ereport_mod.o errormessagelength_mod.o io_constants.o \
	jules_print_mgr.o jules_radiation_mod.o jules_surface_mod.o logging_mod.o string_utils_mod.o urban_param_mod.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

tilepts_jls.o : $(JULES_02)/tilepts_jls.F90 ancil_info.o jules_surface_mod.o jules_surface_types_mod.o parkind1.o \
	um_types.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

time_info_mod.o : $(JULES_02)/time_info_mod.F90 datetime_mod.o model_time_mod.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

top_pdm.o : $(JULES_02)/top_pdm.F90 parkind1.o um_types.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

trifctl.o : $(JULES_02)/trifctl.F90 parkind1.o um_types.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

trif_vars_mod.o : $(JULES_02)/trif_vars_mod.F90 parkind1.o um_types.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

veg3_field_mod.o : $(JULES_02)/veg3_field_mod.F90 ancil_info.o calc_c_comps_triffid_mod.o conversions_mod_jls.o \
	gridbox_mean_mod.o jules_surface_types_mod.o jules_vegetation_mod.o model_time_mod.o pftparm_mod.o \
	prognostics.o timestep_mod.o trif_mod.o um_types.o veg3_param_mod.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

veg_control.o : $(JULES_02)/veg_control.F90 ancil_info.o conversions_mod_jls.o ereport_mod.o \
	errormessagelength_mod.o jules_soil_biogeochem_mod.o jules_soil_mod.o jules_surface_types_mod.o \
	jules_vegetation_mod.o model_time_mod.o parkind1.o soil_biogeochem_control_mod.o soil_ecosse_vars_mod.o \
	timestep_mod.o trif_vars_mod.o trifctl.o um_types.o urban_param_mod.o veg-veg1a_jls_mod.o veg-veg2a_jls_mod.o \
	veg_soil_index_mod.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

water_resources_control_mod.o : $(JULES_02)/water_resources_control_mod.F90 ancil_info.o atm_fields_bounds_mod.o \
	crop_date_mod.o crop_vars_mod.o ereport_mod.o jules_irrig_mod.o jules_soil_mod.o jules_surface_types_mod.o \
	jules_water_resources_mod.o model_time_mod.o parkind1.o theta_field_sizes_mod.o timestep_mod.o um_types.o \
	water_resources_drive.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

fao_evapotranspiration_mod.o : $(JULES_27)/fao_evapotranspiration_mod.F90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

jules.o : $(JULES_03)/jules.F90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

control.o : $(JULES_03)/control.F90 surf_couple_radiation_mod.o surf_couple_explicit_mod.o \
	surf_couple_implicit_mod.o surf_couple_extra_mod.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

jules_fields_mod.o : $(JULES_03)/jules_fields_mod.F90 aero.o ancil_info.o coastal.o crop_vars_mod.o fire_vars_mod.o \
	jules_mod.o p_s_parms.o prognostics.o soil_ecosse_vars_mod.o top_pdm.o trif_vars_mod.o trifctl.o urban_param_mod.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

jules_final_mod.o : $(JULES_03)/jules_final_mod.F90 jules_deposition_mod.o jules_soil_biogeochem_mod.o \
	jules_water_resources_mod.o logging_mod.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

next_time.o : $(JULES_03)/next_time.F90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

jules_vectlib_mod.o : $(JULES_03)/jules_vectlib_mod.F90 parkind1.o um_types.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

atm_fields_bounds_mod.o : $(JULES_03)/var/atm_fields_bounds_mod.F90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

disaggregated_precip.o : $(JULES_03)/var/disaggregated_precip.F90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

ereport_mod.o : $(JULES_03)/var/ereport_mod.F90 logging_mod.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

errormessagelength_mod.o : $(JULES_03)/var/errormessagelength_mod.F90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

forcing.o : $(JULES_03)/var/forcing.F90 parkind1.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

gridmean_fluxes.o : $(JULES_03)/var/gridmean_fluxes.F90 parkind1.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

model_time_mod.o : $(JULES_03)/var/model_time_mod.F90 datetime_mod.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

model_grid_mod.o : $(JULES_03)/var/model_grid_mod.F90 grid_utils_mod.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

timestep_mod.o : $(JULES_03)/var/timestep_mod.F90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

solinc_data.o : $(JULES_03)/var/solinc_data.F90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

zenith_mod.o : $(JULES_03)/zenith_mod.F90 conversions_mod_jls.o datetime_utils_mod.o model_grid_mod.o \
	theta_field_sizes_mod.o time_info_mod.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

allocate_cable_arrays_mod.o : $(JULES_04)/allocate_cable_arrays_mod.F90 ancil_info.o cable_other_constants_mod.o \
	cable_prognostic_info_mod.o cable_surface_types_mod.o cable_types_mod.o ereport_mod.o jules_print_mgr.o \
	jules_snow_mod.o jules_soil_mod.o parkind1.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

allocate_jules_arrays.o : $(JULES_04)/allocate_jules_arrays.F90 aero.o ancil_info.o atm_fields_bounds_mod.o \
	bvoc_vars.o c_z0h_z0m_mod.o coastal.o crop_vars_mod.o cropparm_mod.o deposition_species_mod.o \
	dust_parameters_mod_jls.o fire_vars_mod.o fluxes.o forcing.o gridmean_fluxes.o jules_deposition_mod.o \
	jules_irrig_mod.o jules_mod.o jules_radiation_mod.o jules_rivers_mod.o jules_sea_seaice_mod.o jules_snow_mod.o \
	jules_soil_biogeochem_mod.o jules_soil_mod.o jules_surface_mod.o jules_surface_types_mod.o jules_vegetation_mod.o \
	jules_water_resources_mod.o lake_mod.o metstats_mod.o nlsizes_namelist_mod.o nvegparm_mod.o ozone_vars.o \
	p_s_parms.o parkind1.o pftparm_mod.o prognostics.o soil_ecosse_vars_mod.o switches_urban.o theta_field_sizes_mod.o \
	top_pdm.o trif_mod.o trif_vars_mod.o trifctl.o urban_param_mod.o veg3_field_mod.o veg3_param_mod.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

freeze_soil.o : $(JULES_04)/freeze_soil.F90 conversions_mod_jls.o ereport_mod.o parkind1.o um_types.o \
	water_constants_mod_jls.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

init_flake_ancils_mod.o : $(JULES_05)/init_flake_ancils_mod.F90 dump_mod.o errormessagelength_mod.o input_mod.o \
	io_constants.o jules_surface_mod.o logging_mod.o missing_data_mod.o model_interface_mod.o string_utils_mod.o \
	templating_mod.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

init_ancillaries_mod.o : $(JULES_05)/init_ancillaries_mod.F90 input_mod.o logging_mod.o init_flake_ancils_mod.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

init_grid_mod.o : $(JULES_06)/init_grid_mod.F90 logging_mod.o time_varying_input_mod.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

init.o : $(JULES_39)/init.F90 aero.o allocate_cable_arrays_mod.o ancil_info.o coastal.o crop_vars_mod.o dump_mod.o \
	fire_vars_mod.o init_ancillaries_mod.o init_cable_mod.o init_check_compatibility.o init_deposition_mod.o \
	init_grid_mod.o init_hydrology.o init_jules_sf_diags_mod.o init_model_environment_mod.o init_output_mod.o \
	init_params_mod.o init_parms.o init_plant_n_uptake_mod.o init_radiation.o init_rivers.o init_snow.o init_soil.o \
	init_soil_biogeochem_mod.o init_soil_ecosse_mod.o init_surface.o init_surface_types.o init_vegetation.o \
	init_water_resources_mod.o initial_conditions_mod.o jules_mod.o jules_model_environment_mod.o \
	jules_science_fixes_mod.o jules_surface_mod.o jules_surface_types_mod.o land_tile_ids.o logging_mod.o \
	metstats_init.o model_interface_mod.o model_time_mod.o p_s_parms.o prognostics.o soil_ecosse_vars_mod.o \
	spinup_mod.o time_varying_input_mod.o top_pdm.o trif_vars_mod.o trifctl.o urban_param_mod.o veg3_field_mod.o \
	veg3_param_mod.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

init_params_mod.o : $(JULES_40)/init_params_mod.F90 pftparm_io_mod.o nvegparm_io_mod.o trif_io_mod.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

init_cable_mod.o : $(JULES_39)/init_cable_mod.F90 ancil_info.o cable_types_mod.o jules_fields_mod.o \
	jules_surface_types_mod.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

init_check_compatibility.o : $(JULES_39)/init_check_compatibility.F90 check_unavailable_options_mod.o \
	errormessagelength_mod.o jules_irrig_mod.o jules_radiation_mod.o jules_soil_mod.o jules_surface_types_mod.o \
	jules_vegetation_mod.o logging_mod.o nvegparm_mod.o switches_urban.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

check_unavailable_options_mod.o : $(JULES_03)/check_unavailable_options_mod.F90
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

init_deposition_mod.o : $(JULES_39)/init_deposition_mod.F90 errormessagelength_mod.o io_constants.o \
	jules_deposition_mod.o logging_mod.o string_utils_mod.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

init_drive.o : $(JULES_39)/init_drive.F90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

init_fire.o : $(JULES_39)/init_fire.F90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

init_hydrology.o : $(JULES_39)/init_hydrology.F90 errormessagelength_mod.o io_constants.o jules_hydrology_mod.o \
	logging_mod.o string_utils_mod.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

init_imogen.o : $(JULES_39)/init_imogen.F90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

init_irrigation.o : $(JULES_39)/init_irrigation.F90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

init_jules_sf_diags_mod.o : $(JULES_39)/init_jules_sf_diags_mod.F90 ancil_info.o ereport_mod.o jules_radiation_mod.o \
	jules_sea_seaice_mod.o jules_soil_mod.o jules_vegetation_mod.o parkind1.o sf_diags_mod.o theta_field_sizes_mod.o \
	yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

init_model_environment_mod.o : $(JULES_39)/init_model_environment_mod.F90 jules_model_environment_mod.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

init_output_mod.o : $(JULES_39)/init_output_mod.F90 ancil_info.o datetime_mod.o errormessagelength_mod.o io_constants.o \
	jules_snow_mod.o jules_soil_biogeochem_mod.o jules_soil_ecosse_mod.o jules_vegetation_mod.o \
	jules_water_resources_mod.o logging_mod.o mem_brams_jules.o model_interface_mod.o model_time_mod.o output_mod.o \
	overbank_inundation_mod.o sf_diags_mod.o string_utils_mod.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

init_parms.o : $(JULES_39)/init_parms.F90 ancil_info.o coastal.o fluxes.o infiltration_rate_mod.o jules_mod.o \
	jules_sea_seaice_mod.o jules_surface_mod.o p_s_parms.o prognostics.o sparm_jls_mod.o theta_field_sizes_mod.o \
	urban_param_mod.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

init_plant_n_uptake_mod.o : $(JULES_39)/init_plant_n_uptake_mod.F90 jules_plant_n_uptake_mod.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

init_prescribed_data.o : $(JULES_39)/init_prescribed_data.F90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

init_radiation.o : $(JULES_39)/init_radiation.F90 errormessagelength_mod.o io_constants.o jules_radiation_mod.o \
	logging_mod.o string_utils_mod.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

init_rivers.o : $(JULES_39)/init_rivers.F90 io_constants.o jules_rivers_mod.o logging_mod.o overbank_inundation_mod.o \
	string_utils_mod.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

init_snow.o : $(JULES_39)/init_snow.F90 errormessagelength_mod.o io_constants.o jules_snow_mod.o logging_mod.o \
	string_utils_mod.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

init_soil.o : $(JULES_39)/init_soil.F90 errormessagelength_mod.o io_constants.o jules_soil_mod.o jules_surface_mod.o \
	logging_mod.o mem_brams_jules.o string_utils_mod.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

init_soil_biogeochem_mod.o : $(JULES_39)/init_soil_biogeochem_mod.F90 errormessagelength_mod.o io_constants.o \
	jules_soil_biogeochem_mod.o jules_soil_mod.o logging_mod.o string_utils_mod.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

init_soil_ecosse_mod.o : $(JULES_39)/init_soil_ecosse_mod.F90 ancil_info.o errormessagelength_mod.o io_constants.o \
	jules_soil_biogeochem_mod.o jules_soil_ecosse_mod.o logging_mod.o string_utils_mod.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

init_surface.o : $(JULES_39)/init_surface.F90 errormessagelength_mod.o io_constants.o jules_surface_mod.o \
	logging_mod.o string_utils_mod.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

init_surface_types.o : $(JULES_39)/init_surface_types.F90 cable_surface_types_mod.o errormessagelength_mod.o \
	io_constants.o jules_surface_types_mod.o logging_mod.o string_utils_mod.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

init_urban.o : $(JULES_39)/init_urban.F90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

init_time.o : $(JULES_39)/init_time.F90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

init_vegetation.o : $(JULES_39)/init_vegetation.F90 errormessagelength_mod.o io_constants.o jules_vegetation_mod.o \
	logging_mod.o metstats_mod.o string_utils_mod.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

init_vars_tmp.o : $(JULES_39)/init_vars_tmp.F90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

init_water_resources_mod.o : $(JULES_39)/init_water_resources_mod.F90 errormessagelength_mod.o io_constants.o \
	jules_water_resources_mod.o logging_mod.o string_utils_mod.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

driver_ascii_mod.o : $(JULES_10)/driver_ascii_mod.F90 io_constants.o logging_mod.o string_utils_mod.o \
	file_ascii_generic_sync_mod.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

file_ascii_generic_sync_mod.o : $(JULES_10)/file_ascii_generic_sync_mod.F90 logging_mod.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

file_gridded_mod.o : $(JULES_11)/file_gridded_mod.F90 file_mod.o grid_utils_mod.o io_constants.o logging_mod.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

file_mod.o : $(JULES_12)/file_mod.F90 driver_ascii_mod.o driver_ncdf_mod.o logging_mod.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

input_mod.o : $(JULES_13)/input_mod.F90 grid_utils_mod.o logging_mod.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

file_ts_mod.o : $(JULES_14)/file_ts_mod.F90 datetime_mod.o dictionary_mod.o file_gridded_mod.o grid_utils_mod.o \
	io_constants.o logging_mod.o string_utils_mod.o templating_mod.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

interpolation_mod.o : $(JULES_15)/interpolation_mod.F90 logging_mod.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

imogen_anlg_vals.o : $(JULES_16)/imogen_anlg_vals.F90 missing_data_mod.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

imogen_constants.o : $(JULES_16)/imogen_constants.F90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

imogen_run.o : $(JULES_16)/imogen_run.F90 io_constants.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

imogen_time.o : $(JULES_16)/imogen_time.F90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

cable_maths_constants_mod.o : $(JULES_17)/cable_maths_constants_mod.F90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

cable_other_constants_mod.o : $(JULES_17)/cable_other_constants_mod.F90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

cable_phys_constants_mod.o : $(JULES_17)/cable_phys_constants_mod.F90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

cable_types_mod.o : $(JULES_17)/cable_types_mod.F90 cable_other_constants_mod.o max_dimensions.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

bl_option_mod.o : $(JULES_18)/bl_option_mod.F90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

chemistry_constants_mod.o : $(JULES_18)/chemistry_constants_mod.F90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

conversions_mod_jls.o : $(JULES_18)/conversions_mod_jls.F90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

diag_swchs.o : $(JULES_18)/diag_swchs.F90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

dust_parameters_mod_jls.o : $(JULES_18)/dust_parameters_mod_jls.F90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

io_constants.o : $(JULES_18)/io_constants.F90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

missing_data_mod.o : $(JULES_18)/missing_data_mod.F90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

nlsizes_namelist_mod.o : $(JULES_18)/nlsizes_namelist_mod.F90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

planet_constants_mod_jls.o : $(JULES_18)/planet_constants_mod_jls.F90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

precision_mod.o : $(JULES_18)/precision_mod.F90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

um_types.o : $(JULES_18)/um_types.F90 precision_mod.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

stochastic_physics_run_mod.o : $(JULES_18)/stochastic_physics_run_mod.F90 max_dimensions.o missing_data_mod.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

water_constants_mod_jls.o : $(JULES_18)/water_constants_mod_jls.F90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

canadian_mod.o : $(JULES_19)/canadian_mod.F90 um_types.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

fire_allocate.o : $(JULES_19)/fire_allocate.F90 ancil_info.o fire_mod.o model_grid_mod.o parkind1.o \
	theta_field_sizes_mod.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

fire_calc_daily.o : $(JULES_19)/fire_calc_daily.F90 canadian_mod.o conversions_mod_jls.o fire_mod.o mcarthur_mod.o \
	metstats_mod.o nesterov_mod.o parkind1.o um_types.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

fire_init.o : $(JULES_19)/fire_init.F90 fire_mod.o jules_soil_mod.o metstats_mod.o parkind1.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

fire_mod.o : $(JULES_19)/fire_mod.F90 um_types.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

inferno_io_mod.o : $(JULES_19)/inferno_io_mod.F90 ancil_info.o calc_c_comps_triffid_mod.o fire_vars_mod.o \
	inferno_mod.o jules_surface_types_mod.o jules_vegetation_mod.o parkind1.o pftparm_mod.o qsat_mod.o \
	timestep_mod.o um_types.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

fire_timestep_mod.o : $(JULES_19)/fire_timestep_mod.F90 fire_calc_daily.o fire_mod.o metstats_mod.o parkind1.o \
	um_types.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

mcarthur_mod.o : $(JULES_19)/mcarthur_mod.F90 um_types.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

inferno_mod.o : $(JULES_19)/inferno_mod.F90 ancil_info.o conversions_mod_jls.o jules_soil_biogeochem_mod.o \
	parkind1.o um_types.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

nesterov_mod.o : $(JULES_19)/nesterov_mod.F90 um_types.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

data_parameters_mod.o : $(JULES_20)/data_parameters_mod.F90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

flake_albedo_ref_mod.o : $(JULES_20)/flake_albedo_ref_mod.F90 data_parameters_mod.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

flake_configure_mod.o : $(JULES_20)/flake_configure_mod.F90 data_parameters_mod.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

flake_derivedtypes_mod.o : $(JULES_20)/flake_derivedtypes_mod.F90 data_parameters_mod.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

flake_driver_mod.o : $(JULES_20)/flake_driver_mod.F90 data_parameters_mod.o flake_configure_mod.o flake_mod.o \
	flake_parameters_mod.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

flake_mod.o : $(JULES_20)/flake_mod.F90 data_parameters_mod.o flake_parameters_mod.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

flake_parameters_mod.o : $(JULES_20)/flake_parameters_mod.F90 data_parameters_mod.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

flake_paramoptic_ref_mod.o : $(JULES_20)/flake_paramoptic_ref_mod.F90 data_parameters_mod.o flake_derivedtypes_mod.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

flake_radflux.o : $(JULES_20)/flake_radflux.F90 data_parameters_mod.o flake_derivedtypes_mod.o flake_mod.o \
	flake_parameters_mod.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

jules_flake_interface_1D.o : $(JULES_20)/jules_flake_interface_1D.F90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

blend_h_mod.o : $(JULES_21)/blend_h_mod.F90 um_types.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

c_bvoc_mod.o : $(JULES_21)/c_bvoc_mod.F90 um_types.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

c_kappai_mod.o : $(JULES_21)/c_kappai_mod.F90 um_types.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

c_rmol_mod.o : $(JULES_21)/c_rmol_mod.F90 um_types.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

c_sicehc_mod.o : $(JULES_21)/c_sicehc_mod.F90 um_types.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

c_topog_mod.o : $(JULES_21)/c_topog_mod.F90 um_types.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

c_surf_mod.o : $(JULES_21)/c_surf_mod.F90 um_types.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

c_z0h_z0m_mod.o : $(JULES_21)/c_z0h_z0m_mod.F90 parkind1.o um_types.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

ccarbon_mod.o : $(JULES_21)/ccarbon_mod.F90 um_types.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

cropparm_io_mod.o : $(JULES_21)/cropparm_io_mod.F90 max_dimensions.o um_types.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

cropparm_mod.o : $(JULES_21)/cropparm_mod.F90 parkind1.o um_types.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

csigma_mod.o : $(JULES_21)/csigma_mod.F90 um_types.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

deposition_species_io_mod.o : $(JULES_21)/deposition_species_io_mod.F90 deposition_species_mod.o jules_print_mgr.o \
	max_dimensions.o um_types.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

deposition_species_mod.o : $(JULES_21)/deposition_species_mod.F90 ereport_mod.o jules_deposition_mod.o \
	missing_data_mod.o parkind1.o um_types.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

descent_mod.o : $(JULES_21)/descent_mod.F90 um_types.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

dust_param_mod.o : $(JULES_21)/dust_param_mod.F90 dust_parameters_mod_jls.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

ecosse_param_mod.o : $(JULES_21)/ecosse_param_mod.F90 um_types.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

nvegparm_io_mod.o : $(JULES_21)/nvegparm_io_mod.F90 jules_print_mgr.o max_dimensions.o missing_data_mod.o um_types.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

nvegparm_mod.o : $(JULES_21)/nvegparm_mod.F90 parkind1.o um_types.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

pftparm_io_mod.o : $(JULES_21)/pftparm_io_mod.F90 jules_print_mgr.o max_dimensions.o missing_data_mod.o um_types.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

pftparm_mod.o : $(JULES_21)/pftparm_mod.F90 parkind1.o um_types.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

theta_field_sizes_mod.o : $(JULES_21)/theta_field_sizes_mod.F90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

trif_io_mod.o : $(JULES_21)/trif_io_mod.F90 jules_print_mgr.o max_dimensions.o missing_data_mod.o um_types.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

trif_mod.o : $(JULES_21)/trif_mod.F90 parkind1.o um_types.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

urban_param_mod.o : $(JULES_21)/urban_param_mod.F90 jules_print_mgr.o parkind1.o um_types.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

veg3_param_mod.o : $(JULES_21)/veg3_param_mod.F90 conversions_mod_jls.o jules_surface_types_mod.o jules_vegetation_mod.o \
	model_time_mod.o pftparm_mod.o timestep_mod.o trif_mod.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

veg_param_mod.o : $(JULES_21)/veg_param_mod.F90 um_types.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

Jin11_osa_mod.o : $(JULES_22)/Jin11_osa_mod.F90 ereport_mod.o errormessagelength_mod.o jules_science_fixes_mod.o \
	parkind1.o um_types.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

albpft_jls_mod.o : $(JULES_22)/albpft_jls_mod.F90 ancil_info.o ereport_mod.o jules_radiation_mod.o \
	jules_surface_types_mod.o jules_vegetation_mod.o parkind1.o pftparm_mod.o um_types.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

albsnow_jls_mod.o : $(JULES_22)/albsnow_jls_mod.F90 jules_snow_mod.o jules_surface_mod.o jules_surface_types_mod.o \
	parkind1.o um_types.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

albsnow_ts_jls_mod.o : $(JULES_22)/albsnow_ts_jls_mod.F90 jules_science_fixes_mod.o jules_snow_mod.o parkind1.o \
	um_types.o water_constants_mod_jls.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

calc_direct_albsoil_mod.o : $(JULES_22)/calc_direct_albsoil_mod.F90 parkind1.o um_types.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

canyonalb_mod.o : $(JULES_22)/canyonalb_mod.F90 jules_print_mgr.o matinv_mod.o um_types.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

jules_land_albedo_jls_mod.o : $(JULES_22)/jules_land_albedo_jls_mod.F90 albpft_jls_mod.o albsnow_jls_mod.o \
	albsnow_ts_jls_mod.o ancil_info.o calc_direct_albsoil_mod.o canyonalb_mod.o ereport_mod.o \
	errormessagelength_mod.o jules_print_mgr.o jules_radiation_mod.o jules_snow_mod.o jules_surface_mod.o \
	jules_surface_types_mod.o jules_vegetation_mod.o lake_mod.o missing_data_mod.o nvegparm_mod.o parkind1.o \
	pftparm_mod.o set_soil_alb_components.o switches_urban.o um_types.o water_constants_mod_jls.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

jules_ssi_albedo_jls_mod.o : $(JULES_22)/jules_ssi_albedo_jls_mod.F90 Jin11_osa_mod.o ancil_info.o c_kappai_mod.o \
	ereport_mod.o errormessagelength_mod.o fluxes.o jules_radiation_mod.o jules_science_fixes_mod.o \
	jules_sea_seaice_mod.o parkind1.o theta_field_sizes_mod.o um_types.o water_constants_mod_jls.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

matinv_mod.o : $(JULES_22)/matinv_mod.F90 ereport_mod.o um_types.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

urbanemis_mod.o : $(JULES_22)/urbanemis_mod.F90 jules_print_mgr.o matinv_mod.o parkind1.o um_types.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

set_soil_alb_components.o : $(JULES_22)/set_soil_alb_components.F90 calc_direct_albsoil_mod.o jules_radiation_mod.o \
	parkind1.o um_types.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

areaver_mod.o : $(JULES_23)/areaver_mod.F90 conversions_mod_jls.o parkind1.o um_types.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

overbank_update_mod.o : $(JULES_23)/overbank_update_mod.F90 conversions_mod_jls.o jules_rivers_mod.o \
	overbank_inundation_mod.o planet_constants_mod_jls.o timestep_mod.o um_types.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

rivers_route_mod.o : $(JULES_23)/rivers_route_mod.F90 ancil_info.o ereport_mod.o jules_print_mgr.o jules_rivers_mod.o \
	model_grid_mod.o parallel_mod.o rivers_regrid_mod.o rivers_route_rfm_mod.o rivers_route_trip_mod.o um_types.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

rivers_regrid_mod.o : $(JULES_23)/rivers_regrid_mod.F90 areaver_mod.o jules_rivers_mod.o model_grid_mod.o \
	parallel_mod.o um_types.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

rivers_route_rfm_mod.o : $(JULES_23)/rivers_route_rfm_mod.F90 conversions_mod_jls.o jules_print_mgr.o \
	jules_rivers_mod.o planet_constants_mod_jls.o timestep_mod.o um_types.o water_constants_mod_jls.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

rivers_route_trip_mod.o : $(JULES_23)/rivers_route_trip_mod.F90 jules_rivers_mod.o rivers_route_utils_mod.o \
	timestep_mod.o um_types.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

canopysnow_mod.o : $(JULES_24)/canopysnow_mod.F90 jules_snow_mod.o parkind1.o um_types.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

rivers_route_utils_mod.o : $(JULES_23)/rivers_route_utils_mod.F90 conversions_mod_jls.o jules_rivers_mod.o \
	parkind1.o um_types.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

layersnow_mod.o : $(JULES_24)/layersnow_mod.F90 jules_snow_mod.o parkind1.o um_types.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

compactsnow_mod.o : $(JULES_24)/compactsnow_mod.F90 jules_snow_mod.o parkind1.o planet_constants_mod_jls.o \
	um_types.o water_constants_mod_jls.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

relayersnow_mod.o : $(JULES_24)/relayersnow_mod.F90 ereport_mod.o jules_radiation_mod.o jules_snow_mod.o \
	layersnow_mod.o parkind1.o um_types.o water_constants_mod_jls.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

snow_mod.o : $(JULES_24)/snow_mod.F90 ancil_info.o canopysnow_mod.o compactsnow_mod.o jules_radiation_mod.o \
	jules_snow_mod.o jules_surface_types_mod.o layersnow_mod.o parkind1.o relayersnow_mod.o sf_diags_mod.o \
	snowgrain_mod.o snowpack_mod.o snowtherm_mod.o um_types.o water_constants_mod_jls.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

snowpack_mod.o : $(JULES_24)/snowpack_mod.F90 ancil_info.o jules_snow_mod.o jules_soil_mod.o jules_surface_mod.o \
	jules_surface_types_mod.o lake_mod.o parkind1.o sf_diags_mod.o tridag_mod.o um_types.o water_constants_mod_jls.o \
	yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

snowgrain_mod.o : $(JULES_24)/snowgrain_mod.F90 c_rmol_mod.o conversions_mod_jls.o jules_snow_mod.o parkind1.o \
	um_types.o water_constants_mod_jls.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

snowtherm_mod.o : $(JULES_24)/snowtherm_mod.F90 jules_snow_mod.o parkind1.o um_types.o water_constants_mod_jls.o \
	yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

tridag_mod.o : $(JULES_24)/tridag_mod.F90 parkind1.o um_types.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

calc_baseflow_jls_mod.o : $(JULES_25)/calc_baseflow_jls_mod.F90 jules_hydrology_mod.o jules_print_mgr.o \
	parkind1.o um_types.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

bedrock_jls_mod.o : $(JULES_25)/bedrock_jls_mod.F90 conversions_mod_jls.o jules_soil_mod.o parkind1.o um_types.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

calc_baseflow_jules_mod.o : $(JULES_25)/calc_baseflow_jules_mod.F90 jules_hydrology_mod.o jules_print_mgr.o \
	jules_soil_mod.o parkind1.o um_types.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

calc_fsat_mod.o : $(JULES_25)/calc_fsat_mod.F90 c_topog_mod.o jules_hydrology_mod.o parkind1.o um_types.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

calc_zw_inund_jls_mod.o : $(JULES_25)/calc_zw_inund_jls_mod.F90 jules_hydrology_mod.o parkind1.o um_types.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

calc_zw_jls_mod.o : $(JULES_25)/calc_zw_jls_mod.F90 jules_hydrology_mod.o parkind1.o um_types.o \
	water_constants_mod_jls.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

ch4_tdep_jls_mod.o : $(JULES_25)/ch4_tdep_jls_mod.F90 jules_soil_biogeochem_mod.o parkind1.o um_types.o \
	water_constants_mod_jls.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

ch4_microbe_mod.o : $(JULES_25)/ch4_microbe_mod.F90 jules_soil_biogeochem_mod.o jules_soil_mod.o parkind1.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

ch4_tdep_layers_jls_mod.o : $(JULES_25)/ch4_tdep_layers_jls_mod.F90 jules_soil_biogeochem_mod.o parkind1.o \
	um_types.o water_constants_mod_jls.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

ch4_wetl_jls_mod.o : $(JULES_25)/ch4_wetl_jls_mod.F90 ancil_info.o ch4_microbe_mod.o ch4_tdep_jls_mod.o \
	ch4_tdep_layers_jls_mod.o jules_soil_biogeochem_mod.o jules_soil_mod.o parkind1.o um_types.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

darcy_ch_mod.o : $(JULES_25)/darcy_ch_mod.F90 hyd_con_ch_mod.o jules_soil_mod.o parkind1.o um_types.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

darcy_ic_mod.o : $(JULES_25)/darcy_ic_mod.F90 darcy_ch_mod.o darcy_vg_jls_mod.o jules_soil_mod.o parkind1.o \
	um_types.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

darcy_vg_jls_mod.o : $(JULES_25)/darcy_vg_jls_mod.F90 hyd_con_vg_jls_mod.o jules_soil_mod.o parkind1.o \
	um_types.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

elev_htc_jls_mod.o : $(JULES_25)/elev_htc_jls_mod.F90 ancil_info.o jules_snow_mod.o jules_soil_mod.o \
	parkind1.o um_types.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

frunoff_jls_mod.o : $(JULES_25)/frunoff_jls_mod.F90 parkind1.o um_types.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

gauss_jls_mod.o : $(JULES_25)/gauss_jls_mod.F90 parkind1.o um_types.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

heat_con_jls_mod.o : $(JULES_25)/heat_con_jls_mod.F90 jules_snow_mod.o jules_soil_mod.o parkind1.o \
	um_types.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

hyd_con_ch_mod.o : $(JULES_25)/hyd_con_ch_mod.F90 parkind1.o um_types.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

hyd_con_ic_mod.o : $(JULES_25)/hyd_con_ic_mod.F90 hyd_con_ch_mod.o hyd_con_vg_jls_mod.o jules_soil_mod.o \
	parkind1.o um_types.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

hyd_con_vg_jls_mod.o : $(JULES_25)/hyd_con_vg_jls_mod.F90 parkind1.o um_types.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

hyd_psi_mod.o : $(JULES_25)/hyd_psi_mod.F90 jules_soil_mod.o parkind1.o planet_constants_mod_jls.o um_types.o \
	water_constants_mod_jls.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

hydrol_jls_mod.o : $(JULES_25)/hydrol_jls_mod.F90 ancil_info.o calc_baseflow_jules_mod.o calc_zw_inund_jls_mod.o \
	ch4_wetl_jls_mod.o elev_htc_jls_mod.o ereport_mod.o ice_htc_jls_mod.o jules_hydrology_mod.o jules_irrig_mod.o \
	jules_soil_biogeochem_mod.o jules_soil_mod.o jules_surface_mod.o jules_vegetation_mod.o model_time_mod.o \
	n_leach_mod.o parkind1.o soil_htc_jls_mod.o soil_hyd_jls_mod.o soil_hyd_update_mod.o soil_hyd_wt_mod.o \
	soilmc_jls_mod.o soilt_jls_mod.o surf_hyd_jls_mod.o um_types.o update_mod.o water_constants_mod_jls.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

ice_htc_jls_mod.o : $(JULES_25)/ice_htc_jls_mod.F90 gauss_jls_mod.o jules_snow_mod.o jules_soil_mod.o \
	jules_surface_mod.o parkind1.o um_types.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

infiltration_rate_mod.o : $(JULES_25)/infiltration_rate_mod.F90 ancil_info.o jules_surface_mod.o \
	jules_surface_types_mod.o nvegparm_mod.o parkind1.o pftparm_mod.o um_types.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

n_leach_mod.o : $(JULES_25)/n_leach_mod.F90 ancil_info.o jules_hydrology_mod.o jules_soil_biogeochem_mod.o \
	jules_soil_mod.o jules_surface_types_mod.o parkind1.o um_types.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

pdm_jls_mod.o : $(JULES_25)/pdm_jls_mod.F90 jules_hydrology_mod.o jules_soil_mod.o parkind1.o um_types.o \
	water_constants_mod_jls.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

sieve_jls_mod.o : $(JULES_25)/sieve_jls_mod.F90 parkind1.o um_types.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

soil_htc_jls_mod.o : $(JULES_25)/soil_htc_jls_mod.F90 bedrock_jls_mod.o conversions_mod_jls.o gauss_jls_mod.o \
	heat_con_jls_mod.o jules_irrig_mod.o jules_snow_mod.o jules_soil_mod.o parkind1.o um_types.o \
	water_constants_mod_jls.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

soil_hyd_jls_mod.o : $(JULES_25)/soil_hyd_jls_mod.F90 darcy_ic_mod.o gauss_jls_mod.o hyd_con_ic_mod.o \
	jules_hydrology_mod.o jules_soil_mod.o parkind1.o um_types.o water_constants_mod_jls.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

soil_hyd_update_mod.o : $(JULES_25)/soil_hyd_update_mod.F90 jules_hydrology_mod.o parkind1.o um_types.o \
	water_constants_mod_jls.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

soil_hyd_wt_mod.o : $(JULES_25)/soil_hyd_wt_mod.F90 calc_zw_jls_mod.o parkind1.o um_types.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

soilmc_jls_mod.o : $(JULES_25)/soilmc_jls_mod.F90 jules_soil_mod.o parkind1.o um_types.o water_constants_mod_jls.o \
	yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

soilt_jls_mod.o : $(JULES_25)/soilt_jls_mod.F90 jules_soil_mod.o parkind1.o um_types.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

surf_hyd_jls_mod.o : $(JULES_25)/surf_hyd_jls_mod.F90 ancil_info.o frunoff_jls_mod.o jules_surface_mod.o \
	parkind1.o pdm_jls_mod.o sieve_jls_mod.o um_types.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

vmc_from_head_mod.o : $(JULES_25)/vmc_from_head_mod.F90 jules_soil_mod.o planet_constants_mod_jls.o um_types.o \
	water_constants_mod_jls.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

ecosse_control_mod.o : $(JULES_25)_biogeochem/ecosse/ecosse_control_mod.F90 ancil_info.o ecosse_decomposition_mod.o \
	ecosse_denitrification_mod.o ecosse_leaching_mod.o ecosse_nitrification_mod.o ecosse_prepare_mod.o \
	ecosse_source_sink_mod.o ecosse_utils_mod.o ereport_mod.o jules_soil_ecosse_mod.o jules_surface_types_mod.o \
	parkind1.o soil_ecosse_vars_mod.o um_types.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

ecosse_decomposition_mod.o : $(JULES_25)_biogeochem/ecosse/ecosse_decomposition_mod.F90 ancil_info.o \
	ecosse_param_mod.o ecosse_rate_modifier_mod.o ecosse_utils_mod.o ereport_mod.o jules_soil_ecosse_mod.o parkind1.o \
	string_utils_mod.o um_types.o veg_param_mod.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

ecosse_denitrification_mod.o : $(JULES_25)_biogeochem/ecosse/ecosse_denitrification_mod.F90 ancil_info.o \
	ecosse_param_mod.o ecosse_utils_mod.o jules_soil_ecosse_mod.o parkind1.o um_types.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

ecosse_leaching_mod.o : $(JULES_25)_biogeochem/ecosse/ecosse_leaching_mod.F90 ancil_info.o ecosse_param_mod.o \
	ecosse_utils_mod.o jules_hydrology_mod.o jules_soil_ecosse_mod.o um_types.o water_constants_mod_jls.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

ecosse_nitrification_mod.o : $(JULES_25)_biogeochem/ecosse/ecosse_nitrification_mod.F90 ancil_info.o \
	conversions_mod_jls.o ecosse_param_mod.o ecosse_rate_modifier_mod.o ecosse_utils_mod.o jules_soil_ecosse_mod.o \
	parkind1.o um_types.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

ecosse_prepare_mod.o : $(JULES_25)_biogeochem/ecosse/ecosse_prepare_mod.F90 ancil_info.o conversions_mod_jls.o \
	ecosse_param_mod.o ecosse_utils_mod.o jules_soil_biogeochem_mod.o jules_soil_ecosse_mod.o jules_soil_mod.o \
	jules_surface_types_mod.o jules_vegetation_mod.o parkind1.o pftparm_mod.o root_frac_jls_mod.o \
	soil_ecosse_vars_mod.o trif_mod.o um_types.o veg_param_mod.o vmc_from_head_mod.o water_constants_mod_jls.o \
	yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

ecosse_rate_modifier_mod.o : $(JULES_25)_biogeochem/ecosse/ecosse_rate_modifier_mod.F90 ecosse_param_mod.o \
	jules_soil_biogeochem_mod.o um_types.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

ecosse_source_sink_mod.o : $(JULES_25)_biogeochem/ecosse/ecosse_source_sink_mod.F90 ancil_info.o ecosse_utils_mod.o \
	ereport_mod.o jules_plant_n_uptake_mod.o jules_soil_ecosse_mod.o jules_surface_types_mod.o parkind1.o \
	soil_ecosse_vars_mod.o um_types.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

ecosse_utils_mod.o : $(JULES_25)_biogeochem/ecosse/ecosse_utils_mod.F90 ancil_info.o ecosse_param_mod.o \
	jules_soil_ecosse_mod.o jules_soil_mod.o parkind1.o um_types.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

soil_biogeochem_utils_mod.o : $(JULES_25)_biogeochem/soil_biogeochem_utils_mod.F90 ancil_info.o conversions_mod_jls.o \
	ecosse_utils_mod.o ereport_mod.o jules_soil_biogeochem_mod.o jules_soil_ecosse_mod.o jules_soil_mod.o \
	jules_surface_types_mod.o pftparm_mod.o root_frac_jls_mod.o soil_ecosse_vars_mod.o um_types.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

soil_inorg_n_mod.o : $(JULES_25)_biogeochem/soil_inorg_n_mod.F90 ancil_info.o ereport_mod.o jules_soil_biogeochem_mod.o \
	jules_soil_ecosse_mod.o jules_surface_types_mod.o jules_vegetation_mod.o um_types.o veg_param_mod.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

bvoc_emissions_mod.o : $(JULES_26)/bvoc_emissions_mod.F90 c_bvoc_mod.o conversions_mod_jls.o crop_utils_mod.o \
	jules_surface_mod.o jules_surface_types_mod.o jules_vegetation_mod.o parkind1.o pftparm_mod.o um_types.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

calc_air_dens_jls_mod.o : $(JULES_26)/calc_air_dens_jls_mod.F90 atm_fields_bounds_mod.o parkind1.o \
	planet_constants_mod_jls.o um_types.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

can_drag_mod.o : $(JULES_26)/can_drag_mod.F90 atm_fields_bounds_mod.o conversions_mod_jls.o jules_print_mgr.o \
	jules_surface_mod.o jules_vegetation_mod.o parkind1.o planet_constants_mod_jls.o theta_field_sizes_mod.o um_types.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

cancap_jls.o : $(JULES_26)/cancap_jls.F90 crop_utils_mod.o jules_surface_mod.o jules_surface_types_mod.o \
	jules_vegetation_mod.o parkind1.o pftparm_mod.o um_types.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

dewpnt_jls.o : $(JULES_26)/dewpnt_jls.F90 parkind1.o planet_constants_mod_jls.o qsat_mod.o um_types.o \
	water_constants_mod_jls.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

dustresb_jls.o : $(JULES_26)/dustresb_jls.F90 atm_fields_bounds_mod.o chemistry_constants_mod.o conversions_mod_jls.o \
	dust_parameters_mod_jls.o jules_science_fixes_mod.o parkind1.o planet_constants_mod_jls.o um_types.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

elevate.o : $(JULES_26)/elevate.F90 atm_fields_bounds_mod.o dewpnt_jls.o parkind1.o planet_constants_mod_jls.o \
	qsat_mod.o theta_field_sizes_mod.o um_types.o water_constants_mod_jls.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

fcdch.o : $(JULES_26)/fcdch.F90 atm_fields_bounds_mod.o bl_option_mod.o can_drag_mod.o jules_sea_seaice_mod.o \
	jules_surface_mod.o parkind1.o phi_m_h.o phi_m_h_vol.o planet_constants_mod_jls.o sea_rough_int.o \
	theta_field_sizes_mod.o um_types.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

ice_formdrag_lupkes_mod.o : $(JULES_26)/ice_formdrag_lupkes_mod.F90 atm_fields_bounds_mod.o bl_option_mod.o \
	jules_sea_seaice_mod.o jules_surface_mod.o parkind1.o planet_constants_mod_jls.o um_types.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

get_us.o : $(JULES_26)/get_us.F90 um_types.o urban_param_mod.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

im_sf_pt2_jls.o : $(JULES_26)/im_sf_pt2_jls.F90 atm_fields_bounds_mod.o jules_science_fixes_mod.o \
	jules_sea_seaice_mod.o jules_surface_mod.o parkind1.o planet_constants_mod_jls.o theta_field_sizes_mod.o \
	um_types.o water_constants_mod_jls.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

jules_land_sf_explicit_jls.o : $(JULES_26)/jules_land_sf_explicit_jls.F90 ancil_info.o atm_fields_bounds_mod.o \
	bl_option_mod.o c_elevate.o c_z0h_z0m_mod.o calc_air_dens_jls_mod.o can_drag_mod.o csigma_mod.o dust_param_mod.o \
	elevate.o ereport_mod.o errormessagelength_mod.o fcdch.o generate_anthrop_heat_jls_mod.o heat_con_jls_mod.o \
	jules_irrig_mod.o jules_science_fixes_mod.o jules_sea_seaice_mod.o jules_snow_mod.o jules_soil_biogeochem_mod.o \
	jules_soil_mod.o jules_surface_mod.o jules_surface_types_mod.o jules_vegetation_mod.o metstats_mod.o ozone_vars.o \
	parkind1.o physiol_jls_mod.o planet_constants_mod_jls.o qsat_mod.o sf_diags_mod.o sf_flux_mod.o sf_orog_jls.o \
	sf_resist_jls.o sf_rib.o sfl_int_mod.o snowtherm_mod.o solinc_data.o stdev1.o stochastic_physics_run_mod.o \
	switches_urban.o theta_field_sizes_mod.o tilepts_jls.o timestep_mod.o um_types.o urban_param_mod.o urbanz0.o \
	veg_param_mod.o water_constants_mod_jls.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

jules_land_sf_implicit.jls.o : $(JULES_26)/jules_land_sf_implicit.jls.F90 ancil_info.o atm_fields_bounds_mod.o \
	crop_vars_mod.o csigma_mod.o im_sf_pt2_jls.o jules_print_mgr.o jules_snow_mod.o jules_surface_mod.o \
	jules_surface_types_mod.o parkind1.o planet_constants_mod_jls.o screen_tq_jls.o sf_diags_mod.o sf_evap_jls.o \
	sf_melt_jls.o sice_htf_jls.o solinc_data.o theta_field_sizes_mod.o timestep_mod.o um_types.o \
	water_constants_mod_jls.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

generate_anthrop_heat_jls_mod.o : $(JULES_26)/generate_anthrop_heat_jls_mod.F90 tilepts_jls.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

jules_ssi_sf_explicit_jls.o : $(JULES_26)/jules_ssi_sf_explicit_jls.F90 ancil_info.o atm_fields_bounds_mod.o \
	bl_option_mod.o c_kappai_mod.o calc_air_dens_jls_mod.o csigma_mod.o fcdch.o ice_formdrag_lupkes_mod.o \
	jules_science_fixes_mod.o jules_sea_seaice_mod.o jules_surface_mod.o parkind1.o planet_constants_mod_jls.o \
	qsat_mod.o sf_diags_mod.o sf_flux_mod.o sf_rib.o sfl_int_mod.o stdev1.o theta_field_sizes_mod.o timestep_mod.o \
	um_types.o water_constants_mod_jls.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

jules_ssi_sf_implicit.jls.o : $(JULES_26)/jules_ssi_sf_implicit.jls.F90 ancil_info.o atm_fields_bounds_mod.o \
	csigma_mod.o fluxes.o im_sf_pt2_jls.o jules_sea_seaice_mod.o jules_snow_mod.o jules_surface_mod.o parkind1.o \
	planet_constants_mod_jls.o screen_tq_jls.o sf_diags_mod.o sf_melt_jls.o sice_htf_jls.o theta_field_sizes_mod.o \
	timestep_mod.o um_types.o water_constants_mod_jls.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

leaf_jls_mod.o : $(JULES_26)/leaf_jls_mod.F90 c_rmol_mod.o ereport_mod.o jules_surface_mod.o jules_vegetation_mod.o \
	parkind1.o pftparm_mod.o um_types.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

leaf_limits_mod.o : $(JULES_26)/leaf_limits_mod.F90 ereport_mod.o jules_surface_mod.o jules_vegetation_mod.o \
	parkind1.o pftparm_mod.o planet_constants_mod_jls.o um_types.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

leaf_lit_jls.o : $(JULES_26)/leaf_lit_jls.F90 parkind1.o pftparm_mod.o um_types.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

microbe_jls.o : $(JULES_26)/microbe_jls.F90 ancil_info.o jules_soil_biogeochem_mod.o jules_soil_mod.o parkind1.o \
	sf_diags_mod.o um_types.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

phi_m_h.o : $(JULES_26)/phi_m_h.F90 atm_fields_bounds_mod.o jules_surface_mod.o parkind1.o theta_field_sizes_mod.o \
	um_types.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

phi_m_h_vol.o : $(JULES_26)/phi_m_h_vol.F90 atm_fields_bounds_mod.o jules_surface_mod.o parkind1.o \
	theta_field_sizes_mod.o um_types.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

physiol_jls_mod.o : $(JULES_26)/physiol_jls_mod.F90 albpft_jls_mod.o ancil_info.o atm_fields_bounds_mod.o \
	calc_direct_albsoil_mod.o cancap_jls.o conversions_mod_jls.o cropparm_mod.o ereport_mod.o jules_hydrology_mod.o \
	jules_irrig_mod.o jules_print_mgr.o jules_radiation_mod.o jules_soil_biogeochem_mod.o jules_soil_mod.o \
	jules_surface_mod.o jules_surface_types_mod.o jules_vegetation_mod.o leaf_lit_jls.o metstats_mod.o microbe_jls.o \
	missing_data_mod.o nvegparm_mod.o parkind1.o pftparm_mod.o raero_jls.o root_frac_jls_mod.o set_soil_alb_components.o \
	sf_diags_mod.o sf_stom_jls_mod.o smc_ext_jls_mod.o soil_evap_jls.o stochastic_physics_run_mod.o switches_urban.o \
	theta_field_sizes_mod.o timestep_mod.o um_types.o urban_param_mod.o urbanemis_mod.o water_constants_mod_jls.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

qsat_data_mod.o : $(JULES_26)/qsat_data_mod.F90 um_types.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

qsat_mod.o : $(JULES_26)/qsat_mod.F90 planet_constants_mod_jls.o qsat_data_mod.o um_types.o water_constants_mod_jls.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

raero_jls.o : $(JULES_26)/raero_jls.F90 atm_fields_bounds_mod.o jules_surface_mod.o parkind1.o planet_constants_mod_jls.o \
	theta_field_sizes_mod.o um_types.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

root_frac_jls_mod.o : $(JULES_26)/root_frac_jls_mod.F90 parkind1.o pftparm_mod.o um_types.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

screen_tq_jls.o : $(JULES_26)/screen_tq_jls.F90 atm_fields_bounds_mod.o csigma_mod.o jules_surface_mod.o missing_data_mod.o \
	parkind1.o qsat_mod.o sf_diags_mod.o theta_field_sizes_mod.o um_types.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

sea_rough_int.o : $(JULES_26)/sea_rough_int.F90 atm_fields_bounds_mod.o jules_sea_seaice_mod.o parkind1.o phi_m_h.o \
	planet_constants_mod_jls.o theta_field_sizes_mod.o um_types.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

sf_aero.o : $(JULES_26)/sf_aero.F90 atm_fields_bounds_mod.o dust_param_mod.o dustresb_jls.o jules_science_fixes_mod.o \
	parkind1.o theta_field_sizes_mod.o um_types.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

sf_diags_mod.o : $(JULES_26)/sf_diags_mod.F90 um_types.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

sf_evap_jls.o : $(JULES_26)/sf_evap_jls.F90 ancil_info.o atm_fields_bounds_mod.o jules_irrig_mod.o parkind1.o \
	planet_constants_mod_jls.o sf_diags_mod.o theta_field_sizes_mod.o um_types.o water_constants_mod_jls.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

sf_flux_mod.o : $(JULES_26)/sf_flux_mod.F90 atm_fields_bounds_mod.o csigma_mod.o jules_science_fixes_mod.o \
	jules_sea_seaice_mod.o jules_snow_mod.o jules_surface_mod.o jules_surface_types_mod.o jules_vegetation_mod.o \
	parkind1.o planet_constants_mod_jls.o sf_diags_mod.o switches_urban.o theta_field_sizes_mod.o um_types.o \
	water_constants_mod_jls.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

sf_melt_jls.o : $(JULES_26)/sf_melt_jls.F90 atm_fields_bounds_mod.o jules_snow_mod.o parkind1.o planet_constants_mod_jls.o \
	theta_field_sizes_mod.o um_types.o water_constants_mod_jls.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

sf_orog_gb_jls.o : $(JULES_26)/sf_orog_gb_jls.F90 atm_fields_bounds_mod.o bl_option_mod.o c_surf_mod.o jules_surface_mod.o \
	missing_data_mod.o parkind1.o planet_constants_mod_jls.o sf_diags_mod.o theta_field_sizes_mod.o um_types.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

sf_orog_jls.o : $(JULES_26)/sf_orog_jls.F90 atm_fields_bounds_mod.o bl_option_mod.o c_surf_mod.o jules_surface_mod.o \
	missing_data_mod.o parkind1.o planet_constants_mod_jls.o theta_field_sizes_mod.o um_types.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

sf_resist_jls.o : $(JULES_26)/sf_resist_jls.F90 atm_fields_bounds_mod.o jules_irrig_mod.o jules_snow_mod.o parkind1.o \
	theta_field_sizes_mod.o um_types.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

sf_rib.o : $(JULES_26)/sf_rib.F90 atm_fields_bounds_mod.o parkind1.o planet_constants_mod_jls.o theta_field_sizes_mod.o \
	um_types.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

sf_stom_jls_mod.o : $(JULES_26)/sf_stom_jls_mod.F90 CN_utils_mod.o atm_fields_bounds_mod.o bvoc_emissions_mod.o \
	c_rmol_mod.o ccarbon_mod.o conversions_mod_jls.o crop_utils_mod.o cropparm_mod.o ereport_mod.o jules_surface_mod.o \
		jules_surface_types_mod.o jules_vegetation_mod.o leaf_jls_mod.o leaf_limits_mod.o parkind1.o pftparm_mod.o qsat_mod.o \
		theta_field_sizes_mod.o um_types.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

sfl_int_mod.o : $(JULES_26)/sfl_int_mod.F90 atm_fields_bounds_mod.o ereport_mod.o jules_surface_mod.o parkind1.o phi_m_h.o \
	planet_constants_mod_jls.o sf_diags_mod.o theta_field_sizes_mod.o um_types.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

sice_htf_jls.o : $(JULES_26)/sice_htf_jls.F90 atm_fields_bounds_mod.o c_kappai_mod.o c_sicehc_mod.o jules_print_mgr.o \
	jules_sea_seaice_mod.o missing_data_mod.o parkind1.o sf_diags_mod.o timestep_mod.o um_types.o water_constants_mod_jls.o \
	yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

smc_ext_jls_mod.o : $(JULES_26)/smc_ext_jls_mod.F90 hyd_psi_mod.o jules_vegetation_mod.o parkind1.o pftparm_mod.o um_types.o \
	yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

soil_evap_jls.o : $(JULES_26)/soil_evap_jls.F90 jules_irrig_mod.o parkind1.o um_types.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

stdev1.o : $(JULES_26)/stdev1.F90 atm_fields_bounds_mod.o parkind1.o planet_constants_mod_jls.o theta_field_sizes_mod.o \
	um_types.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

urbanz0.o : $(JULES_26)/urbanz0.F90 get_us.o jules_print_mgr.o jules_surface_types_mod.o um_types.o urban_param_mod.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

vgrav_jls.o : $(JULES_26)/vgrav_jls.F90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

compete_jls.o : $(JULES_27)/compete_jls.F90 descent_mod.o jules_surface_types_mod.o jules_vegetation_mod.o parkind1.o \
	trif_mod.o um_types.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

crop_date_mod.o : $(JULES_27)/crop_date_mod.F90 conversions_mod_jls.o crop_vars_mod.o csigma_mod.o datetime_utils_mod.o \
	jules_surface_types_mod.o time_info_mod.o timestep_mod.o um_types.o water_constants_mod_jls.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

crop_mod.o : $(JULES_27)/crop_mod.F90 ancil_info.o crop_utils_mod.o cropparm_mod.o datetime_utils_mod.o jules_surface_types_mod.o \
	jules_vegetation_mod.o pft_sparm_jls_mod.o tilepts_jls.o time_info_mod.o um_types.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

decay_jls.o : $(JULES_27)/decay_jls.F90 descent_mod.o jules_soil_mod.o parkind1.o um_types.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

develop.o : $(JULES_27)/develop.F90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

dpm_rpm_jls.o : $(JULES_27)/dpm_rpm_jls.F90 jules_surface_types_mod.o parkind1.o trif_mod.o um_types.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

emerge.o : $(JULES_27)/emerge.F90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

jules_subroutine.o : $(JULES_27)/jules_subroutine.F90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

irrigation_mod.o : $(JULES_27)/irrigation_mod.F90 ancil_info.o atm_fields_bounds_mod.o conversions_mod_jls.o crop_date_mod.o \
	crop_vars_mod.o datetime_utils_mod.o ereport_mod.o jules_hydrology_mod.o jules_irrig_mod.o jules_soil_mod.o \
	jules_surface_types_mod.o rivers_route_mod.o theta_field_sizes_mod.o time_info_mod.o timestep_mod.o um_types.o water_constants_mod_jls.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

lotka_eq_jls.o : $(JULES_27)/lotka_eq_jls.F90 jules_surface_types_mod.o jules_vegetation_mod.o parkind1.o trif_mod.o um_types.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

lotka_jls.o : $(JULES_27)/lotka_jls.F90 compete_jls.o jules_surface_types_mod.o jules_vegetation_mod.o parkind1.o pftparm_mod.o \
	trif_mod.o um_types.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

lotka_noeq_jls.o : $(JULES_27)/lotka_noeq_jls.F90 jules_surface_types_mod.o jules_vegetation_mod.o parkind1.o trif_mod.o \
	um_types.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

lotka_noeq_subset_jls.o : $(JULES_27)/lotka_noeq_subset_jls.F90 jules_surface_types_mod.o jules_vegetation_mod.o parkind1.o \
	trif_mod.o um_types.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

partition.o : $(JULES_27)/partition.F90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

mem_brams_jules.o : $(JULES_27)/mem_brams_jules.F90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

pft_sparm_jls_mod.o : $(JULES_27)/pft_sparm_jls_mod.F90 can_drag_mod.o jules_vegetation_mod.o parkind1.o pftparm_mod.o \
	stochastic_physics_run_mod.o um_types.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

phenol_jls.o : $(JULES_27)/phenol_jls.F90 jules_vegetation_mod.o parkind1.o pftparm_mod.o trif_mod.o um_types.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

plant_growth_n.o : $(JULES_27)/plant_growth_n.F90 CN_utils_mod.o descent_mod.o jules_surface_mod.o jules_vegetation_mod.o \
	parkind1.o pftparm_mod.o trif_mod.o um_types.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

soilcarb_jls.o : $(JULES_27)/soilcarb_jls.F90 ancil_info.o decay_jls.o dpm_rpm_jls.o jules_soil_biogeochem_mod.o jules_soil_mod.o \
	jules_surface_types_mod.o jules_vegetation_mod.o parkind1.o um_types.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

soilcarb_mix_jls_mod.o : $(JULES_27)/soilcarb_mix_jls_mod.F90 ancil_info.o conversions_mod_jls.o jules_soil_mod.o parkind1.o \
	um_types.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

soilcarb_layers_jls_mod.o : $(JULES_27)/soilcarb_layers_jls_mod.F90 ancil_info.o decay_jls.o dpm_rpm_jls.o jules_soil_biogeochem_mod.o \
	jules_soil_mod.o jules_surface_types_mod.o jules_vegetation_mod.o parkind1.o pftparm_mod.o root_frac_jls_mod.o soilcarb_mix_jls_mod.o \
	um_types.o veg_param_mod.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

sow.o : $(JULES_27)/sow.F90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

sparm_jls_mod.o : $(JULES_27)/sparm_jls_mod.F90 blend_h_mod.o c_z0h_z0m_mod.o dust_param_mod.o jules_snow_mod.o jules_surface_mod.o \
	jules_surface_types_mod.o jules_vegetation_mod.o nvegparm_mod.o parkind1.o pft_sparm_jls_mod.o stochastic_physics_run_mod.o \
	switches_urban.o um_types.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

veg-veg1a_jls_mod.o : $(JULES_27)/veg-veg1a_jls_mod.F90 ancil_info.o conversions_mod_jls.o infiltration_rate_mod.o \
	jules_surface_types_mod.o jules_vegetation_mod.o parkind1.o phenol_jls.o sparm_jls_mod.o tilepts_jls.o um_types.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

triffid_jls.o : $(JULES_27)/triffid_jls.F90 CN_utils_mod.o ancil_info.o calc_c_comps_triffid_mod.o calc_litter_flux_mod.o \
	jules_soil_biogeochem_mod.o jules_soil_mod.o jules_surface_mod.o jules_surface_types_mod.o jules_vegetation_mod.o lotka_eq_jls.o \
	lotka_jls.o lotka_noeq_jls.o lotka_noeq_subset_jls.o parkind1.o pftparm_mod.o soil_biogeochem_utils_mod.o soil_inorg_n_mod.o \
	soilcarb_jls.o soilcarb_layers_jls_mod.o trif_mod.o um_types.o vegcarb_jls.o woodprod.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

veg-veg2a_jls_mod.o : $(JULES_27)/veg-veg2a_jls_mod.F90 ancil_info.o calc_c_comps_triffid_mod.o conversions_mod_jls.o descent_mod.o \
	infiltration_rate_mod.o jules_soil_biogeochem_mod.o jules_soil_mod.o jules_surface_mod.o jules_surface_types_mod.o \
	jules_vegetation_mod.o parkind1.o phenol_jls.o soilcarb_mix_jls_mod.o sparm_jls_mod.o tilepts_jls.o trif_vars_mod.o triffid_jls.o \
	um_types.o veg_param_mod.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

veg3_litter_mod.o : $(JULES_27)/veg3_litter_mod.F90 veg3_field_mod.o veg3_param_mod.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

veg3_red_dynamic_mod.o : $(JULES_27)/veg3_red_dynamic_mod.F90 veg3_field_mod.o veg3_param_mod.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

veg_soil_index_mod.o : $(JULES_27)/veg_soil_index_mod.F90 jules_surface_types_mod.o jules_vegetation_mod.o um_types.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

vegcarb_jls.o : $(JULES_27)/vegcarb_jls.F90 CN_utils_mod.o jules_surface_mod.o jules_surface_types_mod.o jules_vegetation_mod.o \
	parkind1.o pftparm_mod.o plant_growth_n.o trif_mod.o um_types.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

woodprod.o : $(JULES_27)/woodprod.F90 jules_surface_types_mod.o trif_mod.o um_types.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

irrigation_water_mod.o : $(JULES_28)/irrigation_water_mod.F90 ancil_info.o irrigation_mod.o jules_soil_mod.o jules_surface_types_mod.o \
	jules_water_resources_mod.o timestep_mod.o um_types.o water_constants_mod_jls.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

water_resources_drive.o : $(JULES_28)/water_resources_drive.F90 ancil_info.o atm_fields_bounds_mod.o irrigation_water_mod.o \
	jules_soil_mod.o jules_surface_types_mod.o jules_water_resources_mod.o parkind1.o theta_field_sizes_mod.o um_types.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

#cable_pack_mod.o : $(JULES_29)/cable_pack_mod.F90 ancil_info.o cable_types_mod.o jules_fields_mod.o jules_surface_types_mod.o
#	@cp -f $< $(<F:.f90=.f90)
#	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
#	@mv -f $(<F:.f90=.f90) ../doc/src

datetime_mod.o : $(JULES_29)/datetime_mod.F90 conversions_mod_jls.o datetime_utils_mod.o logging_mod.o precision_mod.o string_utils_mod.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

dictionary_mod.o : $(JULES_29)/dictionary_mod.F90 logging_mod.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

logging_mod.o : $(JULES_29)/logging_mod.F90 errormessagelength_mod.o io_constants.o mem_brams_jules.o string_utils_mod.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

gridbox_mean_mod.o : $(JULES_30)/gridbox_mean_mod.F90 ancil_info.o jules_surface_mod.o jules_surface_types_mod.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

metstats_init.o : $(JULES_30)/metstats/metstats_init.F90 ancil_info.o jules_surface_mod.o jules_surface_types_mod.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

metstats_mod.o : $(JULES_30)/metstats/metstats_mod.F90 conversions_mod_jls.o parkind1.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

metstats_timestep.o : $(JULES_30)/metstats/metstats_timestep.F90 conversions_mod_jls.o dewpnt_jls.o jules_vegetation_mod.o \
	metstats_mod.o parkind1.o qsat_mod.o yomhook.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

string_utils_mod.o : $(JULES_29)/string_utils_mod.F90 io_constants.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

parkind1.o : $(JULES_31)/parkind1.F90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

templating_mod.o : $(JULES_29)/templating_mod.F90 datetime_mod.o io_constants.o string_utils_mod.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

yomhook.o : $(JULES_31)/yomhook.F90 parkind1.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

mpi_mod.o : $(JULES_32)/mpi_mod.F90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

mpi_routines.o : $(JULES_32)/mpi_routines.F90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

driver_ncdf_mod.o : $(JULES_33)/driver_ncdf_mod.F90 io_constants.o logging_mod.o string_utils_mod.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

parallel_mod.o : $(JULES_03)/parallel/parallel_mod.F90 logging_mod.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

output_mod.o : $(JULES_34)/output_mod.F90 data_cube_mod.o datetime_mod.o file_ts_mod.o grid_utils_mod.o io_constants.o \
	logging_mod.o string_utils_mod.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

data_cube_mod.o : $(JULES_07)/data_cube_mod.F90
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

spinup_mod.o : $(JULES_03)/spinup/spinup_mod.F90 logging_mod.o model_interface_mod.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

grid_utils_mod.o : $(JULES_35)/grid_utils_mod.F90 io_constants.o logging_mod.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

initial_conditions_mod.o : $(JULES_08)/initial_conditions_mod.F90 logging_mod.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

time_varying_input_mod.o : $(JULES_36)/time_varying_input_mod.F90 data_cube_mod.o datetime_mod.o file_ts_mod.o \
	input_mod.o io_constants.o logging_mod.o interpolation_mod.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

update_mod.o : $(JULES_03)/update/update_mod.F90 missing_data_mod.o sparm_jls_mod.o imogen_drive_vars.o \
	freeze_soil.o disaggregated_precip.o model_grid_mod.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

dump_mod.o : $(JULES_37)/dump_mod.F90 io_constants.o logging_mod.o file_mod.o model_interface_mod.o \
	output_mod.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

model_interface_mod.o : $(JULES_38)/model_interface_mod.F90 io_constants.o logging_mod.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

ecosse_init_mod.o : $(JULES_08)/ecosse_init_mod.F90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

flake_init_mod.o : $(JULES_08)/flake_init_mod.F90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src
