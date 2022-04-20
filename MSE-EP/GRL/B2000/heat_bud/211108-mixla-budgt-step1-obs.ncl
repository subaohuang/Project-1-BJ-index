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
  diri = "/home/yangsong3/data-observation/GODAS/monthly/float/"


  ; dim_models = dimsizes(model_name)

  opt = False

  his_start = 198001
  his_end   = 202012

;========================================================== Read data
    path = diri
    var_name = (/"pottmp","ucur","vcur","dzdt","thflx"/)
    ; var_name = (/"thflx"/)

    var_name2 = (/"TEMP","UVEL","VVEL","WVEL","SHF"/)

    var_name_out1 = var_name2 +"_mixed_avg"
    var_name_out2 = var_name2 +"_mixed_bot"

    filein = path +var_name
;;;;

;========================================================== thetao
 do i = 0, dimsizes(var_name)-1

    f1 := systemfunc("ls "+ filein(i)+"."+"{1980..2020}"+"*")

    fin_thetao_his := addfiles(f1, "r") 

    time_thetao   = fin_thetao_his[:]->time 

    date :=cd_calendar(time_thetao, 1)
    timeind :=ind(date.ge.his_start.and.date.le.his_end)


    if var_name(i).ne."thflx"
    thetao_his_org := fin_thetao_his[:]->$var_name(i)$(:,0:5,{-30:30},:)
    else
    thetao_his_org := fin_thetao_his[:]->$var_name(i)$(:,{-30:30},:)
    end if 

    printMinMax(thetao_his_org, opt)

    if var_name(i).ne."thflx".and.var_name(i).ne."dzdt" then
    thetao_his := thetao_his_org(timeind,:,:,:)
    ; thetao_his = where(thetao_his.ge.10000.,thetao_his@_FillValue,thetao_his)
    dim_thetao := dimsizes(thetao_his)
    thetao_his_mixed = dim_avg_n_Wrap(thetao_his(:,0:dim_thetao(1)-2,:,:), 1)
    printVarSummary(thetao_his_mixed)
    thetao_his_bot = thetao_his(:,dim_thetao(1)-1,:,:)
    printVarSummary(thetao_his_bot)
    else 
    if var_name(i).eq."dzdt" then 
    thetao_his := thetao_his_org(timeind,0:4,:,:)
    ; thetao_his = where(thetao_his.ge.10000.,thetao_his@_FillValue,thetao_his)

    dim_thetao := dimsizes(thetao_his)
   
    thetao_his_mixed = dim_avg_n_Wrap(thetao_his(:,0:dim_thetao(1)-1,:,:), 1)
    printVarSummary(thetao_his_mixed)
    thetao_his_bot = thetao_his(:,dim_thetao(1)-1,:,:)
    printVarSummary(thetao_his_bot)
    else 
    thetao_his := thetao_his_org(timeind,:,:)
    ; thetao_his = where(thetao_his.ge.10000.,thetao_his@_FillValue,thetao_his)

    dim_thetao := dimsizes(thetao_his)
   
    thetao_his_mixed :=  thetao_his
    end if
    end if

    system("rm ./obs/"+var_name2(i)+"_mixed_godas.nc")
    fout_thetao_his = addfile("./obs/"+var_name2(i)+"_mixed_godas.nc", "c")

    if var_name(i).ne."thflx" then 
    fout_thetao_his->$var_name_out1(i)$ = thetao_his_mixed
    fout_thetao_his->$var_name_out2(i)$ = thetao_his_bot
    else
    fout_thetao_his->$var_name_out1(i)$ = thetao_his_mixed
    end if 
end do
end
