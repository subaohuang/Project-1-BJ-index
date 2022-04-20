;*************************************************
; spec_3.ncl
;
; Concepts illustrated:
;   - Calculating confidence intervals
;************************************************
;
; These files are loaded by default in NCL V6.2.0 and newer
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl" 
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;
; This file still has to be loaded manually
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
;************************************************
begin
;************************************************
; variable and file handling
;************************************************
  yearstart = 195001
  yearend = 201912
  
  runmean_time = 11

  f1 = addfile("/home/ys17-19/lsh/Project/MSE-EP/MSE_necp1.nc", "r")

  date:=cd_calendar(f1->time, 1)

  timeind:=ind(date.le.yearend.and.date.ge.yearstart)

  MSE = f1->MSE(timeind,:,:,:)
                                       
  ;;;MSE 
    plev = MSE&level
    ptop = min(plev)
    psfc = 1013.0
    dz = dpres_plevel(plev, psfc, ptop, 0)

    MSE_col =  dim_sum_wgt_n_Wrap(MSE, dz, 1, 1)

    MSE_col = (/rmMonAnnCycTLL(MSE_col)/)

    MSE_col = (/dtrend_msg_n(ispan(1,dimsizes(MSE_col&time),1),MSE_col,False,True,0)/)

    ; MSE_col = (/runave_n_Wrap(MSE_col, 9*12, 0, 0)/)

    ; MSE_col = dim_avg_n_Wrap(MSE_col, (/1,2/))
    MSE_col = (/MSE_col/10^6/)
    MSE_col_annual = month_to_annual(MSE_col, 1)
    MSE_col_annual!0 = "time"

    MSE_col_MAM = month_to_season(MSE_col, "MAM")
    MSE_col_JJA = month_to_season(MSE_col, "JJA")
    MSE_col_SON = month_to_season(MSE_col, "SON")
    MSE_col_DJF = month_to_season(MSE_col, "DJF")
    
    ; MSE_col_MAM = (/dim_standardize_n_Wrap(MSE_col_MAM, 0, 0)/)

    ; MSE_col_MAM = (/runave_n(MSE_col_MAM, runmean_time, 0, 0)/)

    MSE_nep_mam = dim_avg_n_Wrap(MSE_col_MAM(:,{0:20},{240:280}), (/1,2/))

    MSE_nep_mam =(/runave_n(MSE_nep_mam, runmean_time, 0, 0)/)

  ;;;


;************************************************
; set function arguments
;************************************************
  d   = 0    ; detrending opt: 0=>remove mean 1=>remove mean + detrend
  sm  = 3   ; smooth: should be at least 3 and odd
  pct = 0.10 ; percent taper: (0.0 <= pct <= 1.0) 0.10 common. 
;************************************************
; calculate spectrum
;************************************************
  sdof = specx_anal(MSE_nep_mam(5:dimsizes(MSE_nep_mam)-6),d,sm,pct)
;************************************************
; plotting parameters
;************************************************
   wks  = gsn_open_wks("pdf","11year-rm-spec_MSE_MAM")             ; send graphics to PNG file
   plot = new(1,graphic)                         ; create graphic array 

   r               = True                        ; plot mods desired
   r@gsnDraw       = False                       ; do not draw
   r@gsnFrame      = False                       ; do not advance frame
   r@tiMainString  = "MSE MAM"                       ; title
   r@tiXAxisString = "Period (year)"  ; xaxis
   r@tiYAxisString = "Variance"                  ; yaxis
; ;************************************************
; ; first plot
; ;************************************************
;    plot(0)=gsn_csm_xy(wks, sdof@frq, sdof@spcx, r)  ; create plot 
; ;***********************************************
; ; second plot: Generate quick-n-dirty plot showing 
; ;              5 and 95% "red noise" confidence 
; ;        bounds
; ;***********************************************
   splt = specx_ci(sdof, 0.90, 0.95)             ; calc confidence interval
   ; plot(1) = gsn_csm_xy(wks,sdof@frq, splt,r)
;***********************************************
; Third Plot: Generate fancier plot showing 
;             "red noise" confidence bounds
;             (a) solid for spectrum and Markov, 
;             (b) dash for bounds 
;***********************************************
   r@xyLineThicknesses   = (/2.,1.,1.,1./)       ; Define line thicknesses 
   r@xyDashPatterns      = (/0,0,1,1/)           ; Dash patterns 
   ; r@xyLineColors        = (/"foreground","green","blue","red"/)
   r@xyLineColors        = (/"foreground","green","white","red"/)
   plot = gsn_csm_xy(wks,1/sdof@frq, splt,r)        ;;;
;    splt(0,:) - input spectrum
; splt(1,:) - Markov "Red Noise" spectrum
; splt(2,:) - lower confidence bound for Markov
; splt(3,:) - upper confidence bound for Markov
; ;***********************************************
; ; Fouth Plot: Generate color plot showing 
; ;             "red noise" confidence bounds
; ;***********************************************
; 
;    plot(3) = gsn_csm_xy(wks,sdof@frq, splt,r)
; ;***********************************************
;    res_P             = True                        ; panel mods desired
;    res_P@gsnMaximize = True                        ; blow up plot
;    gsn_panel(wks,plot,(/2,2/),res_P)               ; create panel plots
  draw(plot)
  frame(wks)
end

