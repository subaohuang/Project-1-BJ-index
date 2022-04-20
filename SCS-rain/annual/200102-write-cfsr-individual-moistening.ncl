;;;************************************************************
;; read CFSR latent heating and convection heating
;; convert to nc data
;;avg region over scs and philipian sea
;; linshuheng 2019.12.31
;;;*********************************************************
function get_data(var1,month)
begin
    if month.ne."1" then
    var1_annual = month_to_season(var1, month)
    end if
    
    var1_annual  = month_to_annual(var1, 1)

    var1_annual!0="time" 
    dimz = dimsizes(var1_annual)
    yr_length = dimz(0)

    var1_annual_rc = regCoef_n(ispan(1,yr_length,1), var1_annual,0,0)
    copy_VarCoords(var1_annual(1,:,:,:), var1_annual_rc)
    var1_annual_rc = var1_annual_rc*yr_length

    var1_clim = dim_avg_n_Wrap(var1_annual,0)
    ;;;区域平均

    g_latS=0.
    g_latN=15.
    g_lonW=110.
    g_lonE=160.

    var1_clim_area = wgt_areaave_Wrap(var1_clim(:,{g_latS:g_latN},{g_lonW:g_lonE}),1.0,1.0,0)
    var1_annual_rc_area = wgt_areaave_Wrap(var1_annual_rc(:,{g_latS:g_latN},{g_lonW:g_lonE}),1.0,1.0,0)

    y = var1_annual(:1,:,0,0)
    y(0,:)= (/var1_clim_area/)
    y(1,:)= y(0,:)+(/var1_annual_rc_area/)
  return y
end 
; ==============================================================
; Mission 1. read data
; ==============================================================

 begin
  m=""+ispan(1979,2008,1)+""

  fils=new((/30,12/),string)
  do i=0,29
  fils(i,:)= systemfunc("ls /home/ys17-19/data-observation/CFSR/DIAB_HEAT/diabf01.gdas."+m(i)+"*"+".grb2")
  end do
  

  fils_new=ndtooned(fils)
  fileh=new(360,string)
  fileh=fils_new
  f1=addfiles(fileh, "join")
  printVarSummary(f1)
  
  ;;;bulid a nc file 
    g_data_out_path = "/home/ys17-19/lsh/Project/SCS-rain/individual_moistening_CFSR_ANN_clm_change_profile.nc"
    system("rm "+g_data_out_path)
    ncdf = addfile(g_data_out_path ,"c")  ; open output netCDF file
  ;;;

  ;;read diabiatic heating
     ListSetType (f1, "join")             ; concatenate or "merge" (default)
    temp  = f1[:]->LRGMR_P8_L100_GLL0     ; Large scale condensate moistening rate
    var1  = temp
     delete(temp)
    var1!0="time"
    var1!1="lev"
    var1!2="lat"
    var1!3="lon"
    var1@long_name="Large scale moistening rate"
    
    ; output variables directly
    ncdf->lscmr_q = get_data(var1,"1")

    delete(var1)
   
    temp  = f1[:]->CNVMR_P8_L100_GLL0     ; Deep convective moistening rate 
    var1  = temp

    delete(temp)
    var1!0="time"
    var1!1="lev"
    var1!2="lat"
    var1!3="lon"
    var1@long_name="Deep convective moistening rate"
      ; output variables directly
    ncdf->dcmr_q = get_data(var1,"1")

    delete(var1)

   
    temp  = f1[:]->SHAMR_P8_L100_GLL0     ; Shallow convective moistening rate
    var1  = temp
    delete(temp)
    var1!0="time"
    var1!1="lev"
    var1!2="lat"
    var1!3="lon"
    var1@long_name="Shallow convective moistening rate"
      ; output variables directly
    ncdf->scmr_q = get_data(var1,"1")

    delete(var1)


    temp  = f1[:]->VDFMR_P8_L100_GLL0   ; Vertical diffusion moistening rate
    var1  = temp
    delete(temp)

    var1!0="time"
    var1!1="lev"
    var1!2="lat"
    var1!3="lon"
    var1@long_name="Vertical diffusion moistening rate"
    ; output variables directly
    ncdf->vdmr_q = get_data(var1,"1")

    delete(var1)



   print("finish read data")


  

end 
