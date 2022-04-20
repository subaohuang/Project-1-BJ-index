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
  diri = "/home/ys17-23/Extension/20thCRV2/"


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
    file_name = (/"shtfl","lhtfl","uswrf.sfc","ulwrf.sfc","dswrf.sfc","dlwrf.sfc"/)
    var_name2 = (/"shtfl","lhtfl","uswrf","ulwrf","dswrf","dlwrf"/)

    var_name_out1 = var_name2

    filein = path +file_name
;;;;


;========================================================== 
 do i = 0, dimsizes(var_name2)-1

    f1 := diri +file_name(i)+".mon.mean.nc"

    fin_thetao_his := addfile(f1, "r") 

    time_thetao   = fin_thetao_his->time 

    date :=cd_calendar(time_thetao, 1)
    timeind :=ind(date.ge.his_start.and.date.le.his_end)
    
    thetao_his_org := fin_thetao_his->$var_name2(i)$(timeind,{lat1:lat2},{lon1:lon2})
    printMinMax(thetao_his_org, opt)

    ;;;interpolate 
    xo=fspan(135.25, 279.75, 290)
    yo=fspan(-29.75,29.75,120)
    
    
    print("*****************开始水平插值************************")  
    thetao_his := linint2_Wrap(thetao_his_org&lon, thetao_his_org&lat, thetao_his_org, \
                      True, xo, yo, 0)

    thetao_his!1="lat"
    thetao_his!2="lon"

    thetao_his&lat@units="degrees_north"
    thetao_his&lon@units="degrees_east"


    system("rm ./"+var_name2(i)+"_chazhi.nc")
    fout_thetao_his = addfile("./"+var_name2(i)+"_chazhi.nc", "c")

    
    fout_thetao_his->$var_name_out1(i)$ = thetao_his
end do
end
