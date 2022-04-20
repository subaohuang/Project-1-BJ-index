;;********************************************************
;;画GPCP的降水趋势
;;
;;lnshuheng 2019.12.24
;**********************************************************
function split_12(data)
begin
    dim = dimsizes(data) 
    yearnum = dim(0)/12

    data_12 = new((/12,yearnum,dim(1),dim(2)/),float)
    do i=0,11   
    data_12(i,:,:,:) = data(i::12,:,:)
    end do 
return data_12
end 
; function cal_season_trend(data,season)
;   begin
;     p_djf_r=month_to_season(data, season)

;     if season.eq."DJF"
;       p_djf_r(0:dimsizes(p_djf_r&time)-2,:,:) = p_djf_r(1:,:,:)
;       p_djf_r(dimsizes(p_djf_r&time)-1,:,:) = data(dimsizes(data&time)-1,:,:)
;     end if
;     dim = dimsizes(p_djf_r)
;     p_mk=trend_manken(p_djf_r, False, 0)
;     p_rc = p_mk(1,:,:)
;     copy_VarCoords(p_djf_r(1,:,:), p_rc)
;     p_rc=(/p_rc*dim(0)/)
;     p_rc@prob = p_mk(0,:,:)
;   return p_rc
; end 

begin
  
   
    f2=addfile("/home/yangsong3/data-observation/GPCP/precip.mon.mean.nc", "r")
    date=cd_calendar(f2->time, 1)
    timeind=ind(date.ge.197901.and.date.le.201612)
    p=f2->precip(timeind,:,:)
   
    p = (/rmMonAnnCycTLL(p)/)

    p = (/dtrend_msg_n(ispan(1,dimsizes(p&time),1),p,False,True,0)/)
    ; p_ind = dim_avg_n_Wrap(p(:,{5:25},{60:100}), (/1,2/))
    p_ind_12 = split_12(p)
    p_jjas = dim_avg_n_Wrap(p_ind_12(:,6:7,:,:), 1)

    std  = dim_stddev_n_Wrap(p_jjas, 0)
    ; std = dim_avg_n_Wrap(dim_stddev_n_Wrap(p_ind_12(:,7:8,:,:), 0),0)
    ; 



    
wks   = gsn_open_wks("pdf", "GPCP-PRECIP-STD")
;;**************************************************************
;画图
;;***************************************************************
      
    res=True
    res@gsnFrame=False
    res@gsnDraw=False
    res@gsnLeftStringFontHeightF=20.0
    res@cnFillOn=True
    res@cnLinesOn=False

    res@cnLineLabelsOn=False
    res@lbLabelBarOn=True
    res@lbLabelStride=2
    res@mpLimitMode           = "LatLon"
    res@mpMinLatF             = -30
    res@mpMaxLatF             = 30
    res@mpMinLonF             = 0
    res@mpMaxLonF             = 360
    res@mpCenterLonF          = (res@mpMinLonF+res@mpMaxLonF)/2.0

    res@tmXTOn=False
    res@tmYROn=False
    res@tmXBOn=True
    res@tmYROn=False

    res@tmXBMode="Explicit"
    res@tmYLMode="Explicit"

    ; res@tmXBLabelFontHeightF      = 0.025    ;font height of tick labels
    ; res@tmYLLabelFontHeightF      = 0.025
    ; res@gsnLeftStringFontHeightF = 0.03
    res@tmYLMinorOn = True

    res@lbBoxLinesOn              = True
    res@lbBoxLineThicknessF       = 0.


    ; res@tmXBValues = (/100,120,140,160/) 
    ; res@tmXBLabels = (/100,120,140,160/) +""+"E" 


    ; res@tmXBLabelFontHeightF=0.015
    ; res@tmYLLabelFontHeightF=0.015
    ; res@tmXBMajorThicknessF=1.0
    ; res@tmYLMajorThicknessF=1.0
    ; res@tmBorderThicknessF=1.0

    ; res@mpShapeMode="FreeAspect"
    ; res@vpHeightF=0.3
    ; res@vpWidthF=0.7

    res@cnFillPalette = "precip_diff_12lev"
    res@cnLevelSelectionMode  = "ExplicitLevels" ; set manual contour levels
    ; res@cnLevels = (/-1.5,-1.0,-0.5,-0.25,0.25,0.5,1.0,1.5,2.0,2.5/)
    ; res@cnFillColors= (/4,5,6,7,8,9,10,11,12,13,14/) - 2

    res@cnLevels = ispan(-20,20,4)*0.1
    ; res@cnFillColors= (/5,6,7,8,10,11,12,13,14/) - 2





  res@gsnLeftString="(a) MAM Precip trend " 
  plot = new(1,graphic)
  plot(0) = gsn_csm_contour_map(wks, std, res)


   str1 = unique_string("dum1")  ; "unique_string" will return a unique    
 

 draw(plot)
 frame(wks)
end