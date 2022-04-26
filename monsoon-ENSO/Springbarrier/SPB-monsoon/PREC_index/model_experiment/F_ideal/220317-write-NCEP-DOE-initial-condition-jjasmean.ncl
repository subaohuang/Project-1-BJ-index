
;;*************************************************************
;Aim: Output global forcing in certain period
;Mission:
;
;   1. Read data from datasets
;   2. Post-post process to deserved results
;   3. Output the data
;                     
;*************************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/popRemap.ncl" 
function split_12(data)
begin
    dim = dimsizes(data) 
    yearnum = dim(0)/12

    data_12 = new((/12,yearnum,dim(1)/),float)
    do i=0,11   
    data_12(i,:,:) = data(i::12,:)
    end do 
return data_12
end 

function get_clm_data(diri,var,timestart,timeend)
begin

  fils = diri + var+ ".mon.mean.nc"

  f1 = addfile(fils, "r")
  time = f1->time 
  date = cd_calendar(time, 1)
  timeind = ind(date.ge.timestart.and.date.le.timeend)

  data = f1->$var$(timeind,::-1,{-90:90},{0:360})
  data_clm = clmMonTLLL(data)

  data_clm_jjas = dim_avg_n_Wrap(data_clm(5:8,:,:,:), 0)

return data_clm_jjas
end 

begin
 
; ==============================================================
; Mission 1. read data
; ==============================================================
;;;read UVT monthly data from linshuheng

  diri = "/home/ys17-23/Extension/NCEP-NCAR/monthly/"
  var = (/"uwnd","vwnd","air"/)
  timestart = 197901
  timeend = 202012
  u_clm_jjas = get_clm_data(diri,var(0),timestart,timeend)
  v_clm_jjas = get_clm_data(diri,var(1),timestart,timeend)  
  t_clm_jjas = get_clm_data(diri,var(2),timestart,timeend)
  t_clm_jjas = (/t_clm_jjas + 273.15/)
  printVarSummary(u_clm_jjas)
;;;;


;;; read data from model intial file 

  fils = systemfunc("ls /home/ys17-23/lsh/Project/monsoon-ENSO/Springbarrier/SPB-monsoon/PREC_index/model_experiment/F_ideal/cami_0000-01-01_1.9x2.5_L30_c070703.nc") ;;;daily数据计算的q1

  f3 = addfile(fils, "r")

  US = f3->US
  VS = f3->VS
  T = f3->T

  hyam = f3->hyam
  hybm = f3->hybm
  P0 = f3->P0
  PSI = f3->PS
  print(P0)
;;;;

;;;;read notopo data 

  fils = systemfunc("ls /home/ys17-23/lsh/Project/monsoon-ENSO/Springbarrier/SPB-monsoon/PREC_index/model_experiment/F_ideal/F_dynamic_f19_ctrl_20years.cam.h0.0002-01.nc") ;;;daily数据计算的q1

  f4 = addfile(fils, "r")

  PS = f4->PS(0,:,:)

  printVarSummary(P0)
  printVarSummary(PS)
  printVarSummary(PSI)


  printVarSummary(US)
  printVarSummary(VS)

;;;;

;;;=========================================================
;;;intepolate t 
;;;=========================================================

  ;;; Horizontal
  lati = t_clm_jjas&lat
  lato = T&lat
  
  loni = t_clm_jjas&lon
  lono = T&lon

  to_hr = linint2_Wrap(loni,lati,t_clm_jjas,True, lono,lato, 0)
  printVarSummary(to_hr)
  ;;;

  ;;;Vertical pres2hybrid
  psfc = PS
  pi = 100.0*t_clm_jjas&level 
  printVarSummary(pi)

  to_final   = pres2hybrid(pi,psfc,P0,to_hr,hyam,hybm,4)  ;;;ti four demensions
  printVarSummary(to_final)
;;;;

;;;=========================================================
;;;intepolate US,VS
;;;=========================================================
  ;;; Horizontal
  lati := u_clm_jjas&lat
  lato := US&slat
  
  loni = u_clm_jjas&lon
  lono = US&lon

  USo_hr = linint2_Wrap(loni,lati,u_clm_jjas, True, lono,lato, 0)
  printVarSummary(USo_hr)

  ;;;for v
  lati := v_clm_jjas&lat
  lato := VS&lat  

  loni := v_clm_jjas&lon
  lono := VS&slon

  VSo_hr = linint2_Wrap(loni,lati,v_clm_jjas, True, lono,lato, 0)
  printVarSummary(VSo_hr)
  

  ;;;for PS 
  lati := PS&lat
  lato := US&slat

  loni := US&lon
  lono := US&lon
  
  PS_hr1 = linint2_Wrap(loni,lati,PS, True, lono,lato, 0)
  printVarSummary(PS_hr1)

  

  lati := PS&lat
  lato := VS&lat

  loni := PS&lon
  lono := VS&slon
  
  PS_hr2 = linint2_Wrap(loni,lati,PS, True, lono,lato, 0)
  printVarSummary(PS_hr2)
  
;;;;

  ;;;Vertical
  pi := 100.0*u_clm_jjas&level 
   
  USo_final   = pres2hybrid_Wrap(pi,PS_hr1,P0,USo_hr,hyam,hybm,4)
  printVarSummary(USo_final)

  VSo_final   = pres2hybrid_Wrap(pi,PS_hr2,P0,VSo_hr,hyam,hybm,4)
  printVarSummary(VSo_final)
  ;;;
    ;;;
  delete(u_clm_jjas)
  delete(v_clm_jjas)
  delete(t_clm_jjas)
;;;;

  T_diff = T
  T_diff(0,:,:,:) = T(0,:,:,:) -to_final

  ; printMinMax(T_diff, 1)
  ; exit()

  T(0,:,:,:) = (/to_final/)
  US(0,:,:,:) = (/USo_final/)
  VS(0,:,:,:) = (/VSo_final/)

;=====================================================================

;;==================================================================
;;;WRITE DATA
;;================================================================

 fils4 = systemfunc("ls //home/ys17-23/lsh/Project/monsoon-ENSO/Springbarrier/SPB-monsoon/PREC_index/model_experiment/F_ideal/cami_0000-01-01_1.9x2.5_L30_c070703_JJASMEAN_noto.nc") 
 f4 = addfile(fils4, "w")

 


;overwrite a netCDF file

  f4->T = (/T/)
  f4->US = (/US/)
  f4->VS = (/VS/)
  f4->PS(0,:,:) = (/PS/)
  print("Done!")



end 
   


