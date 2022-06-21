;;*************************************************
;; draw circulation pattern during monsoon year 
;;*************************************************************
function split_12(data)
begin
    dim = dimsizes(data) 
    yearnum = dim(0)/12

    data_12 = new((/12,yearnum/),float)
    do i=0,11   
    data_12(i,:) = data(i::12)
    end do 
return data_12
end 

undef("remove_enso_signal")
function remove_enso_signal(nino3_4,  data)
  begin  
   ;;;;;
      ; data_de=(/dtrend_msg_n(ispan(1,dimsizes(data&time),1),data,False,True,0)/) 
      data_de  = data
      data_rc = regCoef_n(nino3_4,data_de, 0, 0)
      copy_VarCoords(data(1,:,:), data_rc)
      data_rc_3d  = conform_dims(dimsizes(data), data_rc, (/1,2/))
      m = reshape(data_rc@yintercept,dimsizes(data_rc_3d)) 
      nino3_4_3d = conform(data_rc_3d,nino3_4,(/0/))
  
      data_relate_enso  = nino3_4_3d*data_rc_3d + m 
      data_new = data 
      data_new = (/data - data_relate_enso/)
  return data_new
end 
undef("enso_eof")
function enso_eof(data)
  begin  
     ;;;;EOF 
    optEOF      = True
    optEOF@jopt = 0
    neval  = 3                   ; calculate 2 EOFs out of 7 ;specifies the number of eigenvalues and eigenvectors to be returned
    ev_sst = eofunc_n_Wrap(data,neval,optEOF,0)    ;;;得到的标准化后的场，要使得不标准化
    ev_sst&lon@units = "degrees_east"
    ev_sst&lat@units = "degrees_north"
    ev_timef = eofunc_ts_n_Wrap(data, ev_sst, False,0)

    do ne=0,neval-1
       ev_sst(ne,:,:) = ev_sst(ne,:,:)*sqrt(ev_sst@eval(ne)) ;;;使得量级能够相同
    end do

    eval=ev_sst@eval
    ts_model = ev_timef(0,:)/eval(0)
    ; ts_model = ev_timef(0,:)
    ; ts_model=ev_timef(1,:)/eval(1)
    ; tos_enino = ev_sst(0,:,:)
    ts_model_std = dim_standardize(ts_model, 0)
    tos_enino  =  regCoef_n(ts_model_std, data, 0, 0)
    tos_enino@var = ev_sst@pcvar
    copy_VarCoords(data(1,:,:), tos_enino)

    method = ev_sst@method
   
  return [/tos_enino,ts_model_std/]
end 

;;*************************************************
;; main program 
;;*************************************************************
begin
    startmon = 101
    endmon = 99912
    startyear = 1
    endyear = 999

;;;nino3
    lat1 = -5
    lat2 = 5
    lon1 = 190
    lon2 = 240
;;;;

; ;;;nino4
;     lat1 = -5
;     lat2 = 5
;     lon1 = 150
;     lon2 = 200
; ;;;;

; ;;;读入海温数据
    ncout = addfile("ssta_picontrol_highpass_bw.nc", "r")
    ssta = ncout->ssta 
;;


;;;;;
    nino34 = dim_avg_n_Wrap(ssta(:,{lat1:lat2},{lon1:lon2}), (/1,2/))
    nino34_12_r =  split_12(nino34)
    nino34_djf = month_to_season(nino34, "DJF")

    nino34_djf(0:dimsizes(nino34_djf)-2) = nino34_djf(1:) 
    nino34_djf(dimsizes(nino34_djf)-1) = nino34(dimsizes(nino34)-1) 

    nino34_djf = (/dim_avg_n_Wrap(nino34_12_r(5:8,:), 0)/)  ;;;; JJAS 

    std_nino = dim_stddev(nino34_djf)   
    sst_total_year = ispan(startyear, endyear, 1) ;;;0-199
;;;;

;;;;;
    nino3 = dim_avg_n_Wrap(ssta(:,{-5:5},{210:270}), (/1,2/))
    nino3_12_r =  split_12(nino3)
    nino3_djf = month_to_season(nino3, "DJF")

    nino3_djf(0:dimsizes(nino3_djf)-2) = nino3_djf(1:) 
    nino3_djf(dimsizes(nino3_djf)-1) = nino3(dimsizes(nino3)-1) 

    ; nino3_djf = (/dim_avg_n_Wrap(nino3_12_r(5:8,:), 0)/)  ;;;; JJAS 
;;;;

;;;;;
    nino4 = dim_avg_n_Wrap(ssta(:,{-5:5},{160:210}), (/1,2/))
    nino4_12_r =  split_12(nino4)
    nino4_djf = month_to_season(nino4, "DJF")

    nino4_djf(0:dimsizes(nino4_djf)-2) = nino4_djf(1:) 
    nino4_djf(dimsizes(nino4_djf)-1) = nino4(dimsizes(nino4)-1) 

    nino4_djf = (/dim_avg_n_Wrap(nino4_12_r(5:8,:), 0)/)  ;;;; JJAS 
;;;;

;;;;;
    nino1_2 = dim_avg_n_Wrap(ssta(:,{-10:0},{270:280}), (/1,2/))
    nino1_2_12_r =  split_12(nino1_2)
    nino1_2_djf = month_to_season(nino1_2, "DJF")

    nino1_2_djf(0:dimsizes(nino1_2_djf)-2) = nino1_2_djf(1:) 
    nino1_2_djf(dimsizes(nino1_2_djf)-1) = nino1_2(dimsizes(nino1_2)-1) 

    nino1_2_djf = (/dim_avg_n_Wrap(nino1_2_12_r(5:8,:), 0)/)  ;;;; JJAS 
;;;;

; ;;; Kao and Yu 2009

;    ;;;去掉Nino4去信号
    
;     ; ssta_rm_nino1_2 = ssta(:,{-20:20},{120:280})
;     ; ssta_rm_nino1_2 = (/remove_enso_signal(nino1_2,ssta(:,{-20:20},{120:280}))/)

;     ssta_rm_nino4 = ssta(:,{-20:20},{120:280})
;     ssta_rm_nino4 = (/remove_enso_signal(nino4,ssta(:,{-20:20},{120:280}))/)
;   ;;;
;    print("******start EOF*************")
;      ;;;;EOF 
;     ; eof_cp = enso_eof(ssta_rm_nino1_2)
;     ; eof_cp_pattern = eof_cp[0]
;     ; ; eof_cp_pattern = (/-1*eof_cp_pattern/)
;     ; pc_cp = eof_cp[1]
;     ; var1 = eof_cp_pattern@var

;     eof_ep = enso_eof(ssta_rm_nino4)
;     eof_ep_pattern = eof_ep[0]
;     pc_ep = eof_ep[1]
;     var2 = eof_ep_pattern@var

;     copy_VarCoords(nino1_2, pc_ep)
;     ; copy_VarCoords(nino1_2, pc_cp)
; ;;;;
;   ;;输出PC的数据
;     ; opt="CP_PC1_picontrol_bw.ascii"    
;     ; asciiwrite(opt, sprintf ("%6.2f,",pc_cp))  

;     opt="EP_PC1_picontrol_bw.ascii"    
;     asciiwrite(opt, sprintf ("%6.2f,",pc_ep))  

;  print("******end EOF*************")

;  exit()
;;;    
;;;;;;

;;;;
  ;;输入PC的数据
    opt="CP_PC1_picontrol_bw.ascii"    
    pc_cp =  asciiread(opt, -1,"float")  

    opt="EP_PC1_picontrol_bw.ascii"    
    pc_ep =  asciiread(opt, -1, "float")  

    tos_enino_ep  =  regCoef_n(pc_ep, ssta, 0, 0)
    copy_VarCoords(ssta(1,:,:), tos_enino_ep)

    tos_enino_cp  =  regCoef_n(pc_cp, ssta, 0, 0)
    copy_VarCoords(ssta(1,:,:), tos_enino_cp)

;;;;;

; ;;;pattern 
;     pc_ep_jjas  =  month_to_season(pc_ep, "DJF")
;     pc_ep_jjas  =  (/dim_avg_n_Wrap(pc_ep(5:8,:), 0)/)  ;;;; JJAS 
    
;     pc_cp_jjas  =  month_to_season(pc_cp, "DJF")
;     pc_cp_jjas  =  (/dim_avg_n_Wrap(pc_cp(5:8,:), 0)/)  ;;;; JJAS 




   ; wks = gsn_open_wks("pdf", "ATL warming-pattern EOF-CP&EPelnino")
   wks = gsn_open_wks("pdf", "Picontrol-EOF-CP&EPelnino")

    plot = new(2,graphic)

;;;;;;
    res=True
    res@gsnDraw=False
    res@gsnFrame=False
    res@gsnAddCyclic=False


    res@cnLineColors="Transparent"  ;
    res@lbLabelBarOn=False
    res@cnFillOn = True 
    res@cnLinesOn = False

    res@tmXBMode="Explicit"
    res@tmYLMode="Explicit"

    res@tmXBLabelFontHeightF      = 0.025    ;font height of tick labels
    res@tmYLLabelFontHeightF      = 0.025

    res@gsnContourNegLineDashPattern = 1
    res@gsnContourPosLineDashPattern = 0
    res@gsnContourZeroLineThicknessF = 2.0



    res@cnLineColor = "black"
    res@cnLineLabelInterval = 2
    res@cnLineLabelsOn=False

    ; res@cnFillPalette = "precip_diff_12lev"
    res@cnFillPalette = "BlueWhiteOrangeRed"
    ; res@cnFillPalette = "BlueDarkRed18"

    res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour level
    res@cnMinLevelValF       =   -0.6
    res@cnMaxLevelValF       =   0.6
    res@cnLevelSpacingF      =  0.1
    ; res@cnFillPalette = "NCV_jaisnd"
    ; res@cnFillPalette = "sunshine_9lev"
    

    res@tmXTOn=False
    res@tmYROn=False

    res@mpLimitMode           = "LatLon"
    res@mpMinLatF             = -30.0
    res@mpMaxLatF             = 30.0
    res@mpMinLonF             = 125.0
    res@mpMaxLonF             = 280.0
   
    res@mpCenterLonF          = (/res@mpMinLonF + res@mpMaxLonF/)/2


    res@mpFillOn = True
    res@mpGeophysicalLineColor      = "Grey18"
    res@mpGeophysicalLineThicknessF = 1.5
    

    ; res@gsnRightString = var1(0) +""
    res@gsnLeftString = "CTRL CP El Nino"
    res@gsnLeftString = "ATL warming CP El Nino"

    res@gsnLeftStringFontHeightF = 0.03

    ; res@vpHeightF = 0.3
    ; res@vpWidthF = 0.7
    res@mpShapeMode = "FreeAspect"

    plot(0) = gsn_csm_contour_map(wks,tos_enino_ep,res)
    ; res@gsnLeftString = "CTRL EP El Nino"
    res@gsnLeftString = "ATL warming EP El Nino"
    ; res@gsnRightString =var2(0) +""
    plot(1) = gsn_csm_contour_map(wks,tos_enino_cp,res)


    pres = True
    pres@gsnMaximize = True
    pres@gsnPanelLabelBar = True
    pres@pmLabelBarHeightF = 0.05
    pres@pmLabelBarWidthF  = 0.6
    pres@pmLabelBarOrthogonalPosF = -0.02
    ; pres@lbOrientation = "Vertical"
    pres@lbLabelFontHeightF = 0.01
    gsn_panel(wks, plot, (/2,1/), pres)


end 
