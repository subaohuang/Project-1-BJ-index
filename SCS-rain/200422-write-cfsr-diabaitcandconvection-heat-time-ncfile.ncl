;;;************************************************************
;; read CFSR latent heating and convection heating
;; convert to nc data
;;avg region over scs and philipian sea
;; linshuheng 2019.12.31
;;;*********************************************************
; ==============================================================
; Mission 1. read data
; ==============================================================
 begin
  diri="/mnt/f/share/MERRA/MERRA/"
   ;;;读入MERRA的数据
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
  ;;read diabiatic heating
     ListSetType (f1, "join")             ; concatenate or "merge" (default)
    temp  = f1[:]->LRGHR_P8_L100_GLL0     ; Large scale condensate heating rate
    var1  = temp
    delete(temp)
    print("var1")
    temp  = f1[:]->CNVHR_P8_L100_GLL0     ; Deep convective heating rate 
    var1  = var1+temp
    delete(temp)
    temp  = f1[:]->SHAHR_P8_L100_GLL0     ; Shallow convective heating rate
    var1  = var1+temp
    delete(temp)
    temp  = f1[:]->VDFHR_P8_L100_GLL0     ; Vertical diffusion heating rate
    var1  = var1+temp
    delete(temp)
    temp  = f1[:]->SWHR_P8_L100_GLL0     ; Solar radiative heating rate
    var1  = var1+temp
    delete(temp)
    temp  = f1[:]->LWHR_P8_L100_GLL0     ; Long-wave radiative heating rate
    var1  = var1+temp
    delete(temp)

    var1!0="time"
    var1!1="lev"
    var1!2="lat"
    var1!3="lon"
    var1@long_name="Diabatic Process Heating Rate"
print("finish read data")
; ==============================================================
; Mission 2. post-post process data to deserved results
; ==============================================================

;-----------Follow is the vertical integral-----------
;    pver=var1&lev
;    dp=dpres_plevel(pver,100000,100,0)              ; caculate dp
;    dp_cfm=conform_dims(dimsizes(var1),dp,1)        ; confirm dp to caculate sum(var*dp)
;    wgt_var1=var1(:,0,:,:)
;    wgt_var1=dim_sum_n(dp_cfm*var1,1)               ; sum(var*dp)
;    wgt_var1=86400*wgt_var1/sum(dp)                 ; Integral=sum(var*dp)/sum(dp)
;    printVarSummary(wgt_var1)
;-----------Follow to caculate the trend--------------
    var1_annual = month_to_annual(var1,1)
    ; var1_annual = month_to_season(var1, "MAM")

    var1_annual!0="time" 
    dimz = dimsizes(var1_annual)
    yr_length = dimz(0)

    var1_annual_rc = regCoef_n(ispan(1,yr_length,1), var1_annual,0,0)
    copy_VarCoords(var1_annual(1,:,:,:), var1_annual_rc)
    var1_annual_rc = var1_annual_rc*yr_length

    ; var1_clim = dim_avg_n_Wrap(var1_annual,0)
    ;;;区域平均

    g_latS=0.
    g_latN=15.
    g_lonW=110.
    g_lonE=160.
    
    var1_annual_area = var1_annual(:,:,1,1)
    var1_annual_area = wgt_areaave_Wrap(var1_annual(:,:,{g_latS:g_latN},{g_lonW:g_lonE}),1.0,1.0,0)

    
    g_data_out_path = "/home/ys17-19/lsh/Project/SCS-rain/diabatic_and_deepcon_heating_CFSR_ANN_SCS.nc"

    system("rm "+g_data_out_path)
    ncdf = addfile(g_data_out_path ,"c")  ; open output netCDF file
    ; output variables directly
    ncdf->diab_q_time = var1_annual_area

;;;读入 deep convection 的数据
    print("***********deep convection数据*******************")

    deepcon_q = f1[:]->CNVHR_P8_L100_GLL0     ; Deep convective heating rate 
    var2 = deepcon_q

    var2_annual = month_to_annual(var2,1)
    var2_annual!0="time" 
    dimz = dimsizes(var2_annual)
    yr_length = dimz(0)
    var2_annual_area = var2_annual(:,:,1,1)
    var2_annual_area = wgt_areaave_Wrap(var2_annual(:,:,{g_latS:g_latN},{g_lonW:g_lonE}),1.0,1.0,0)

    ncdf->deepcon_q_time = var2_annual_area
;;;;
end 
