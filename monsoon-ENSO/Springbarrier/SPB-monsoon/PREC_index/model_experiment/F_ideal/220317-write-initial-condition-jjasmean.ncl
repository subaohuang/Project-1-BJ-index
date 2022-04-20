
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

  fils = diri + var+ "_mon_1x1_1979_2020.nc"

  f1 = addfile(fils, "r")
  time = f1->time 
  date = cd_calendar(time, 1)
  timeind = ind(date.ge.timestart.and.date.le.timeend)

  data = short2flt(f1->$var$(timeind,:,{-90:90},{0:360}))
  data_clm = clmMonTLLL(data)

  ; data_clm_jjas = dim_avg_n_Wrap(data_clm(5:8,:,:,:), 0)
  data_clm_jjas = data_clm(0,:,:,:)

return data_clm_jjas
end 

begin
 
; ==============================================================
; Mission 1. read data
; ==============================================================
;;;read UVT monthly data from linshuheng

  diri = "/home/ys17-23/Extension/ERA5/ERA5-monthly/pressure/"
  var = (/"u","v","t"/)
  timestart = 197901
  timeend = 202012
  u_clm_jjas = get_clm_data(diri,var(0),timestart,timeend)
  v_clm_jjas = get_clm_data(diri,var(1),timestart,timeend)  
  t_clm_jjas = get_clm_data(diri,var(2),timestart,timeend)
  printVarSummary(u_clm_jjas)
;;;;


;;; read data from model intial file 

  fils = systemfunc("ls /home/ys17-23/lsh/Project/monsoon-ENSO/Springbarrier/SPB-monsoon/PREC_index/model_experiment/F_ideal/cami_0000-01-01_1.9x2.5_L30_c070703.nc") 

  f3 = addfile(fils, "r")

  US = f3->US
  VS = f3->VS
  T = f3->T

  hyam = f3->hyam
  hybm = f3->hybm
  P0 = f3->P0
  PS = f3->PS

  printVarSummary(P0)
  printVarSummary(PS)

  printVarSummary(US)
  printVarSummary(VS)

;;;;

;;;=========================================================
;;;intepolate t 
;;;=========================================================
  ;;;Vertical pres2hybrid
  psfc = PS

  pi = 100.0*t_clm_jjas&level 
  printVarSummary(pi)
           
  ti = t_clm_jjas
  
  
  to   = pres2hybrid_Wrap(pi,psfc,p0,ti,hyam,hybm,0)
  printVarSummary(to)
  ;;;;

  ;;; Horizontal
  lati = t_clm_jjas&latitude
  lato = T&lat
  
  loni = t_clm_jjas&longitude
  lono = T&lon

  to_final = linint2_Wrap(loni,lati,to, True, lono,lato, 0)
  printVarSummary(to_final)
  ;;;
;;;;

;;;=========================================================
;;;intepolate US,VS
;;;=========================================================
  ;;;Vertical
  pi := 100.0*u_clm_jjas&level 
           
  Ui := u_clm_jjas

  USo   = pres2hybrid_Wrap(pi,psfc,p0,Ui,hyam,hybm,0)
  printVarSummary(USo)
  
  Vi := v_clm_jjas

  VSo   = pres2hybrid_Wrap(pi,psfc,p0,Vi,hyam,hybm,0)
  printVarSummary(VSo)
  
  ;;;;

  ;;; Horizontal
  lati := u_clm_jjas&latitude
  lato := US&slat
  
  loni = u_clm_jjas&longitude
  lono = US&lon

  USo_final = linint2_Wrap(loni,lati,USo, True, lono,lato, 0)
  printVarSummary(USo_final)

  ;;;for v
  lati := v_clm_jjas&latitude
  lato := VS&lat  
  VSo_final = linint2_Wrap(loni,lati,VSo, True, lono,lato, 0)
  printVarSummary(VSo_final)
  ;;;
  delete(u_clm_jjas)
  delete(v_clm_jjas)
  delete(t_clm_jjas)


  T_diff = T
  T_diff(0,:,:,:) = T(0,:,:,:) - flt2dble(to_final)

  printMinMax(T_diff, 1)
  ; exit()

  T(0,:,:,:) = (/flt2dble(to_final)/)
  US(0,:,:,:) = (/flt2dble(USo_final)/)
  VS(0,:,:,:) = (/flt2dble(VSo_final)/)

;=====================================================================

;;==================================================================
;;;WRITE DATA
;;================================================================

 ; fils4 = systemfunc("ls //home/ys17-23/lsh/Project/monsoon-ENSO/Springbarrier/SPB-monsoon/PREC_index/model_experiment/F_ideal/cami_0000-01-01_1.9x2.5_L30_c070703_JJASMEAN.nc") 
 fils4 = systemfunc("ls //home/ys17-23/lsh/Project/monsoon-ENSO/Springbarrier/SPB-monsoon/PREC_index/model_experiment/F_ideal/cami_0000-01-01_1.9x2.5_L30_c070703_era5.nc") 
 f4 = addfile(fils4, "w")



;overwrite a netCDF file

  f4->T = (/T/)
  f4->US = (/US/)
  f4->VS = (/VS/)
  print("Done!")



end 
   


