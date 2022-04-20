load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/popRemap.ncl" 
function  rc_test(ssta_rc)
begin
  tval = onedtond(ssta_rc@tval , dimsizes(ssta_rc))
   df   = onedtond(ssta_rc@nptxy, dimsizes(ssta_rc)) - 2
   b = tval    ; b must be same size as tval (and df)
   b = 0.5
   prob = betainc(df/(df+tval^2),df/2.0,b)    ; prob(nlat,nlon)
   copy_VarCoords(ssta_rc, prob)
return prob
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
    neval  = 10                   ; calculate 2 EOFs out of 7 ;specifies the number of eigenvalues and eigenvectors to be returned
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
    
  return [/tos_enino,ts_model_std/]
end 


begin
;读入数据
 ;;;time span
    startmon =25601
    endmon = 28012

    startyear_r = str_split_by_length(startmon, 3)  ;; 将时间string分解成年和月
    endyear_r = str_split_by_length(endmon, 3 )

    startyear = stringtoint(startyear_r(0))
    endyear = stringtoint(endyear_r(0))

    yearnum = endyear - startyear +1
  ;;
;;;region
   lat1 = -90
   lat2 = 90
   lon1 = 0
   lon2 = 360
;;;;

;;;;读入控制实验结果
    ; fils1 = "/home/yangsong3/data-observation/linshh/CESM-data/B/B2000_alt_north_year_CTRL/B2000_alt_north_year_CTRL.TEMP_taux_tauy.h1.0251-0280.nc"\
    fils1 = "/home/yangsong3/CMIP6/linshh/CESM-data/B/lsh_B2000_alt_north_db_year_80_20/lsh_B2000_alt_north_db_year_80_20.TEMP.h1.0251-0280.nc"
    ; fils1 = "/home/yangsong3/data-model/CESM_CTRL/B2000_F19G16_CAM4_CTRL/pre/B2000_f19g16_CP_CTRL.TEMP.fvgrid.0251-0280.nc"
    f1=addfile(fils1, "r")
    time = f1->time 
    time = (/time-1/)
    date=cd_calendar(time, 1)

    timeind=ind(date.le.endmon.and.date.ge.startmon)

    sst_r = f1->TEMP(timeind,0,:,:)  
    ; sst = f1->TEMP(timeind,0,:,:)  
   ;;差值TEMP 
    sst = PopLatLon(sst_r,"gx1v6","fv1.9x2.5","bilin","da","090206")
    printVarSummary(sst)
    
    sst_ct = sst(:,{lat1:lat2},{lon1:lon2})

    sst_ct = (/rmMonAnnCycTLL(sst_ct)/)


    sst_ct_djf = month_to_season(sst_ct, "DJF")
    sst_ct_djf(0:dimsizes(sst_ct_djf&time)-2,:,:) = sst_ct_djf(1:dimsizes(sst_ct_djf&time)-1,:,:)
    sst_ct_djf(dimsizes(sst_ct_djf&time)-1,:,:) = sst_ct(dimsizes(sst_ct&time)-1,:,:)


    tos_djf = sst_ct_djf
;;;;


    ;;;去掉Nino4去信号
    nino1_2 = dim_avg_n_Wrap(tos_djf(:,{0:-10},{270:280}), (/1,2/))
    nino4 = dim_avg_n_Wrap(tos_djf(:,{-5:5},{160:210}), (/1,2/))
    
    tos_djf_rm_nino12 = tos_djf
    tos_djf_rm_nino12 = (/remove_enso_signal(nino1_2,tos_djf)/)

    tos_djf_rm_nino4 = tos_djf
    tos_djf_rm_nino4 = (/remove_enso_signal(nino4,tos_djf)/)


    ;;;;EOF 
    eof_cp = enso_eof(tos_djf)
    eof_cp_pattern = eof_cp[0]
    ; eof_cp_pattern = (/eof_cp_pattern/)

    pc_cp = eof_cp[1]

    var1 = eof_cp_pattern@var

    eof_ep = enso_eof(tos_djf)
    eof_ep_pattern = eof_ep[0]
    pc_ep = eof_ep[1]

    var2 = eof_ep_pattern@var
    


    ; ;;输出PC的数据
    ; opt="CP_PC1_B2000_alt_north_year_ATL.ascii"    
    ; asciiwrite(opt, sprintf ("%6.2f,",pc_cp))  

    ; opt="EP_PC1_B2000_alt_north_year_ATL.ascii"    
    ; asciiwrite(opt, sprintf ("%6.2f,",pc_ep))  

    ; print(ind(pc_cp.ge.1.0))
    ; print(ind(pc_ep.ge.1.0))

;;;;


; ;;;输出数据

;     diri = "/home/ys17-19/lsh/Project/ITCZ_bias/data/CMIP5/EOF/CP/"


;     system("rm "+diri+"CP-EOF-ssta-djf-"+model+".nc")
;     out = addfile(diri+"CP-EOF-ssta-djf-"+model+".nc", "c")
;     out->ts_model = ts_model
;     out->pattern_model = tos_enino
;     out->pcvar = ev_sst@pcvar(0)
; ;;;
   ; wks = gsn_open_wks("pdf", "ATL warming-pattern EOF-CP&EPelnino")
   ; wks = gsn_open_wks("pdf", "controlrun-pattern-EOF1-elnino")
   wks = gsn_open_wks("pdf", "ATLwarmingrun-pattern-EOF1-elnino")

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
    ; res@cnMinLevelValF       =   -1
    ; res@cnMaxLevelValF       =   1
    ; res@cnLevelSpacingF      =  0.2
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
    

    res@gsnRightString = var1(0) +""
    res@gsnLeftString = "(a) Control run EOF1 pattern of SST"
    ; res@gsnLeftString = "ATL warming CP El Nino"

    res@gsnLeftStringFontHeightF = 0.02
    res@gsnRightStringFontHeightF = 0.02

    ; res@vpHeightF = 0.3
    ; res@vpWidthF = 0.7
    res@mpShapeMode = "FreeAspect"

    plot(0) = gsn_csm_contour_map(wks,eof_cp_pattern,res)
    ; res@gsnLeftString = "CTRL EP El Nino"
    res@gsnLeftString = "ATL warming EP El Nino"
    res@gsnRightString =var2(0) +""
    plot(1) = gsn_csm_contour_map(wks,eof_ep_pattern,res)


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
