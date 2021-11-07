load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
;==========================================================
begin
;========================================================== Basics settings
  ; model_name = (/"B2000_alt_north_year_CTRL"/)
  ; diri = "/home/yangsong3/data-model/lsh/CESM/B/"

  model_name = (/"lsh_B2000_alt_north_db_year_80_20"/)
  diri = "/home/yangsong3/CMIP6/linshh/CESM-data/B/"


  dim_models = dimsizes(model_name)

  opt = False

  his_start = 25601
  his_end   = 28012

;========================================================== Read data
    path = diri+model_name+"/"
    var_name = (/"TEMP","UVEL","VVEL","WVEL","SHF"/)

    var_name_out1 = (/"TEMP","UVEL","VVEL","WVEL","SHF"/) +"_mixed_avg"
    var_name_out2 = (/"TEMP","UVEL","VVEL","WVEL","SHF"/) +"_mixed_bot"

    filein = path + model_name+"."+var_name+"."+"fvgrid.0251-0280.nc"


;========================================================== thetao
 do i = 0, dimsizes(var_name)-1

    fin_thetao_his := addfile(filein(i), "r")

    time_thetao   = fin_thetao_his->time 
    time_thetao = (/time_thetao-1/)

    date :=cd_calendar(time_thetao, 1)
    timeind :=ind(date.ge.his_start.and.date.le.his_end)

    thetao_his_org := fin_thetao_his->$var_name(i)$

    printVarSummary(thetao_his_org)
    printMinMax(thetao_his_org, opt)

    if var_name(i).ne."SHF".and.var_name(i).ne."WVEL" then
    thetao_his := thetao_his_org(timeind,:,{-30:30},:)
    ; thetao_his = where(thetao_his.ge.10000.,thetao_his@_FillValue,thetao_his)
    dim_thetao := dimsizes(thetao_his)
    thetao_his_mixed = dim_avg_n_Wrap(thetao_his(:,0:dim_thetao(1)-2,:,:), 1)
    printVarSummary(thetao_his_mixed)
    thetao_his_bot = thetao_his(:,dim_thetao(1)-1,:,:)
    printVarSummary(thetao_his_bot)

    else 
    if var_name(i).eq."WVEL" then 
    thetao_his := thetao_his_org(timeind,:,{-30:30},:)
    ; thetao_his = where(thetao_his.ge.10000.,thetao_his@_FillValue,thetao_his)

    dim_thetao := dimsizes(thetao_his)
   
    thetao_his_mixed = dim_avg_n_Wrap(thetao_his(:,0:dim_thetao(1)-1,:,:), 1)
    printVarSummary(thetao_his_mixed)
    thetao_his_bot = thetao_his(:,dim_thetao(1)-1,:,:)
    printVarSummary(thetao_his_bot)
    else 
    thetao_his := thetao_his_org(timeind,{-30:30},:)
    ; thetao_his = where(thetao_his.ge.10000.,thetao_his@_FillValue,thetao_his)

    dim_thetao := dimsizes(thetao_his)
   
    thetao_his_mixed :=  thetao_his
    end if
    end if

    system("rm ./cesmmodel/"+model_name+"/"+var_name(i)+"_mixed_"+model_name+".nc")
    fout_thetao_his = addfile("./cesmmodel/"+model_name+"/"+var_name(i)+"_mixed_"+model_name+".nc", "c")

    if var_name(i).ne."SHF" then 
    fout_thetao_his->$var_name_out1(i)$ = thetao_his_mixed
    fout_thetao_his->$var_name_out2(i)$ = thetao_his_bot
    else
    fout_thetao_his->$var_name_out1(i)$ = thetao_his_mixed
    end if 
end do
end
