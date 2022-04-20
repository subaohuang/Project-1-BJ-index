 ;;;************************************************************
;; read CFSR latent heating and convection heating
;; convert to nc data
;;avg region over scs and philipian sea
;; linshuheng 2019.12.31
;;;*********************************************************
undef("remove_enso_signal")
function remove_enso_signal(data)
  begin  
    ;;;读入HadiSST海温数据
    diri="/home/ys17-19/lsh/data/sst/"
    fils1=systemfunc("ls "+diri+"HadISST_sst.nc")
    f1=addfile(fils1, "r")
    date:=cd_calendar(f1->time, 1)

    timeind:=ind(date.le.201012.and.date.ge.197901)
    sst =  lonFlip(f1->sst(timeind,:,:))
    ssta=(/rmMonAnnCycTLL(sst)/)
    copy_VarCoords(sst, ssta)
    ssta_de = ssta
    ssta_de=(/dtrend_msg_n(ispan(1,dimsizes(ssta&time),1),ssta,False,True,0)/)
    ssta_index = dim_avg_n_Wrap(ssta_de(:,{5:-5},{190:240}), (/1,2/))
    printVarSummary(ssta_index)
    nino3_4 = month_to_season(ssta_index, "DJF")   ;;;第一年是JF平均，第二年是DJF平均
    printVarSummary(nino3_4)
    ; nino3_4_new = nino3_4                ;;;;前期ENSO
    ; nino3_4_new(0:33) = (/nino3_4(1:)/)
    ; nino3_4_new(34) =ssta_index(35*12-1) 
    nino3_4 = (/dim_standardize_n(nino3_4, 1, 0)/)
;;;;;
    dim = dimsizes(data)
    data_de=(/dtrend_msg_n(ispan(1,dim(0),1),data,False,True,0)/) 
    data_rc = regCoef_n(nino3_4,data_de, 0, 0)
    copy_VarCoords(data(1,:,:,:), data_rc)
    data_rc_3d  = conform_dims(dimsizes(data), data_rc, (/1,2,3/))
    m = data_rc_3d 
    m = (/data_rc@yintercept(0)/)

    data_relate_enso  = data_de*data_rc_3d + m 
    data_new = data 
    data_new = (/data - data_relate_enso/)
  return data_new
end 
; ==============================================================
; Mission 1. read data
; ==============================================================
 begin

 ;;;decide season 

  season = "DJF"
  
  diri="/mnt/f/share/MERRA/MERRA/"
   ;;;读入MERRA的数据
  m=""+ispan(1979,2010,1)+""

  fils=new((/32,12/),string)
  do i=0,31
  fils(i,:)= systemfunc("ls /home/ys17-19/data-observation/CFSR/DIAB_HEAT/diabf01.gdas."+m(i)+"*"+".grb2")
  end do
  

  fils_new=ndtooned(fils)
  fileh=new(32*12,string)
  fileh=fils_new
  f1=addfiles(fileh, "join")
  printVarSummary(f1)
  ;;read diabiatic heating
     ListSetType (f1, "join")             ; concatenate or "merge" (default)
    temp  = f1[:]->LRGHR_P8_L100_GLL0     ; Large scale condensate heating rate
    var1  = temp
    delete(temp)
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

    ; var1_annual = month_to_annual(var1,1)

    var1_annual = month_to_season(var1, season)
    var1_annual!0="time" 
    var1_clim = dim_avg_n_Wrap(var1_annual,0)


    dimz = dimsizes(var1_annual)
    yr_length = dimz(0)

    var1_annual = (/remove_enso_signal(var1_annual)/)

    var1_annual_rc = regCoef_n(ispan(1,yr_length,1), var1_annual,0,0)
    copy_VarCoords(var1_annual(1,:,:,:), var1_annual_rc)
    var1_annual_rc = var1_annual_rc*yr_length

  
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

    g_data_out_path = "/home/ys17-19/lsh/Project/SCS-rain/season/diabatic_heating_CFSR_"+season+"_clm_change_profile_2010_remove_enso.nc"

    system("rm "+g_data_out_path)
    ncdf = addfile(g_data_out_path ,"c")  ; open output netCDF file
    ; output variables directly
    ncdf->diab_q = y

;;;读入 deep convection 的数据
    print("***********deep convection数据*******************")

    deepcon_q = f1[:]->CNVHR_P8_L100_GLL0     ; Deep convective heating rate 
    var2 = deepcon_q

    ; var2_annual = month_to_annual(var2,1)
    var2_annual = month_to_season(var2, season)
    var2_annual!0="time" 
    var2_clim = dim_avg_n_Wrap(var2_annual,0)


    dimz = dimsizes(var2_annual)
    yr_length = dimz(0)

    var2_annual = (/remove_enso_signal(var2_annual)/)

    var2_annual_rc = regCoef_n(ispan(1,yr_length,1), var2_annual,0,0)
    copy_VarCoords(var2_annual(1,:,:,:), var2_annual_rc)
    var2_annual_rc = var2_annual_rc*yr_length

    
    var2_clim_area = wgt_areaave_Wrap(var2_clim(:,{g_latS:g_latN},{g_lonW:g_lonE}),1.0,1.0,0)

    var2_annual_rc_area = wgt_areaave_Wrap(var2_annual_rc(:,{g_latS:g_latN},{g_lonW:g_lonE}),1.0,1.0,0)

    y2 = var2_annual(:1,:,0,0)
    y2(0,:)= (/var2_clim_area/)
    y2(1,:)= y2(0,:)+(/var2_annual_rc_area/)

    g_data_out_path = "/home/ys17-19/lsh/Project/SCS-rain/season/deepcon_heating_CFSR_"+season+"_clm_change_profile_2010_remove_enso.nc"

    system("rm "+g_data_out_path)
    ncdf = addfile(g_data_out_path ,"c")  ; open output netCDF file
    ; output variables directly
    ncdf->deepcon_q = y2
;;;;
end 
