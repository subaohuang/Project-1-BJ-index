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
  diri = "/home/ys17-23/CMIP6/2.5x2.5/CESM2/piControl/Omon/"


  ; dim_models = dimsizes(model_name)

  opt = False

  ; his_start = 101
  ; his_end   = 99912


  his_start = 20101
  his_end   = 120012

  lat1 = -30.0
  lat2 = 30.0
  lon1 = 120.0
  lon2  = 280.0

;========================================================== Read data
    path = diri
    ; var_name = (/"pottmp","ucur","vcur","dzdt","hfds"/)

    var_name = (/"thetao","uo","vo","wo","hfds"/)
    ; var_name = (/"hfds"/)

    ; var_name2 = (/"TEMP","UVEL","VVEL","WVEL","SHF"/)
    var_name2 = (/"TEMP","UVEL","VVEL","WVEL","SHF"/)

    var_name_out1 = var_name2 +"_mixed_avg"
    var_name_out2 = var_name2 +"_mixed_bot"

    filein = path +var_name
;;;;

 ;;lev: 500, 1500, 2500, 3500, 4500, 5500, 6500 
 ;;wo: 0, 1000, 2000, 3000, 4000, 5000, 6000 
;;;
;========================================================== thetao
 do i = 0, dimsizes(var_name)-1
    
    ; if var_name(i).ne."thetao" then 
    ; f1 := systemfunc("ls " + diri +var_name(i)+"_Omon_CESM2_piControl_r1i1p1f1_gn_0[0-9]"+"*")
    ; else 
    ; f1 := systemfunc("ls " + diri +var_name(i)+"_Omon_CESM2_piControl_r1i1p1f1_gr_0[0-9]"+"*")
    ; end if 

    if var_name(i).ne."thetao" then 
    f1 := systemfunc("ls " + diri +var_name(i)+"_Omon_CESM2_piControl_r1i1p1f1_gn_"+"*")
    else 
    f1 := systemfunc("ls " + diri +var_name(i)+"_Omon_CESM2_piControl_r1i1p1f1_gr_"+"*")
    end if 

    fin_thetao_his := addfiles(f1, "r") 
    time_thetao   = fin_thetao_his[:]->time
    printVarSummary(time_thetao)

    date :=cd_calendar(time_thetao, 1)

    timeind :=ind(date.ge.his_start.and.date.le.his_end)
   
    if var_name(i).ne."hfds"
    thetao_his_org := lonFlip(fin_thetao_his[:]->$var_name(i)$(timeind,0:5,{lat1:lat2},:))
    else
    thetao_his_org := lonFlip(fin_thetao_his[:]->$var_name(i)$(timeind,{lat1:lat2},:))
    end if 
    
    printVarSummary(thetao_his_org)
    printMinMax(thetao_his_org, opt)

    if var_name(i).ne."hfds" then
      print(i)
    thetao_his := thetao_his_org(:,:,:,{lon1:lon2})
    ; thetao_his = where(thetao_his.ge.10000.,thetao_his@_FillValue,thetao_his)
    dim_thetao := dimsizes(thetao_his)
    thetao_his_mixed = dim_avg_n_Wrap(thetao_his(:,0:dim_thetao(1)-2,:,:), 1) ;;;5-45 mean
    printVarSummary(thetao_his_mixed)
    thetao_his_bot = thetao_his(:,dim_thetao(1)-1,:,:)
    printVarSummary(thetao_his_bot)
    else 
    thetao_his := thetao_his_org(:,:,{lon1:lon2})
    thetao_his_mixed = thetao_his
    end if

    system("rm ./"+var_name2(i)+"_mixed_picontrol_201-1200.nc")
    fout_thetao_his = addfile("./"+var_name2(i)+"_mixed_picontrol_201-1200.nc", "c")

    if var_name(i).ne."hfds" then 
    fout_thetao_his->$var_name_out1(i)$ = thetao_his_mixed
    fout_thetao_his->$var_name_out2(i)$ = thetao_his_bot
    else
    fout_thetao_his->$var_name_out1(i)$ = thetao_his_mixed
    end if 
end do
end
