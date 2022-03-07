load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
;==========================================================
begin
;========================================================== Basics settings
  ; model_name = (/"B2000_alt_north_year_CTRL"/)
  ; diri = "/home/yangsong3/data-model/lsh/CESM/B/"

  ; model_name = (/"lsh_B2000_alt_north_db_year_80_20"/)
  diri = "/home/ys17-23/Extension/SODA_region/"


  ; dim_models = dimsizes(model_name)

  opt = False

  his_start = 187101
  his_end   = 201012

  lat1 = -30.0
  lat2 = 30.0
  lon1 = 135.0
  lon2  = 280.0

;========================================================== Read data
    path = diri
    ; var_name = (/"pottmp","ucur","vcur","dzdt","thflx"/)

    var_name = (/"TEMP","U","V","W"/)
    ; var_name = (/"thflx"/)

    ; var_name2 = (/"TEMP","UVEL","VVEL","WVEL","SHF"/)
    var_name2 = (/"TEMP","UVEL","VVEL","WVEL"/)

    var_name_out1 = var_name2 +"_mixed_avg"
    var_name_out2 = var_name2 +"_mixed_bot"

    filein = path +var_name
;;;;

 ;;lev: 0:4  5.0 15.07 25.28 35.76 46.61 57.98
;;;
;========================================================== thetao
 do i = 0, dimsizes(var_name)-1

    f1 := diri +"soda_"+var_name(i)+"_1871-2010.nc"

    fin_thetao_his := addfile(f1, "r") 

    time_thetao   = fin_thetao_his->TIME
    time_thetao = (/time_thetao - 1/)
   printVarSummary(time_thetao)

    date :=cd_calendar(time_thetao, 1)
    timeind :=ind(date.ge.his_start.and.date.le.his_end)
   
    if var_name(i).ne."thflx"
    thetao_his_org := fin_thetao_his->$var_name(i)$(timeind,0:5,{lat1:lat2},{lon1:lon2})  
    else
    thetao_his_org := fin_thetao_his->$var_name(i)$(timeind,{lat1:lat2},{lon1:lon2})
    end if 
    
    delete(thetao_his_org&LEV1_6)
    delete(thetao_his_org&LAT92_211)
    delete(thetao_his_org&LON271_570)
    delete(thetao_his_org&TIME)

    thetao_his_org!0 = "time" 
    thetao_his_org!1 = "depth" 
    thetao_his_org!2 = "lat" 
    thetao_his_org!3 = "lon" 
    
    thetao_his_org&time = (/tofloat(fin_thetao_his->TIME)/)
    thetao_his_org&time@units = time_thetao@units
  
   
    thetao_his_org&depth = (/tofloat(fin_thetao_his->LEV1_6)/)
    thetao_his_org&lat = (/fspan(-29.75,29.75,120)/)
    thetao_his_org&lon = (/fspan(135.25, 279.75, 290)/)

    printVarSummary(thetao_his_org)

    printMinMax(thetao_his_org, opt)

    if var_name(i).ne."thflx" then
    thetao_his := thetao_his_org(timeind,:,:,:)
    ; thetao_his = where(thetao_his.ge.10000.,thetao_his@_FillValue,thetao_his)
    dim_thetao := dimsizes(thetao_his)
    thetao_his_mixed = dim_avg_n_Wrap(thetao_his(:,0:dim_thetao(1)-2,:,:), 1) ;;;5-45 mean
    printVarSummary(thetao_his_mixed)
    thetao_his_bot = thetao_his(:,dim_thetao(1)-1,:,:)
    printVarSummary(thetao_his_bot)
    else 
    end if

    system("rm ./"+var_name2(i)+"_mixed_soda.nc")
    fout_thetao_his = addfile("./"+var_name2(i)+"_mixed_soda.nc", "c")

    if var_name(i).ne."thflx" then 
    fout_thetao_his->$var_name_out1(i)$ = thetao_his_mixed
    fout_thetao_his->$var_name_out2(i)$ = thetao_his_bot
    else
    fout_thetao_his->$var_name_out1(i)$ = thetao_his_mixed
    end if 
end do
end
