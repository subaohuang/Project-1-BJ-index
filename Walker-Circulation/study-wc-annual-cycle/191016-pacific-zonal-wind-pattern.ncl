;;;************************************************************
;; draw the annual cycle of u10 
;；探究不同经度的10-m
;; linshuheng 2019.10.16
;;;*********************************************************
function annual_cycle(u10_inten)
local u10_inten_month, nmos,nmo
begin
    u10_inten_month = new(12,typeof(u10_inten))
    
    nmos = 12
    ntim = 456
    do nmo=0,nmos-1
    u10_inten_month(nmo) = dim_avg( u10_inten(nmo:ntim-1:nmos) )
    end do
return u10_inten_month
end
;;---------------------------------------
;;read  wind and slp data
;----------------------------------------
begin
    fils = systemfunc("ls /home/ys17-19/lsh/data/era/sur_1.5/era_interim_sur1.5_"+"*"+".nc")
    f3 = addfiles(fils, "r")

    time=f3[:]->time
    date=cd_calendar(time, 1)
    tind=ind(date.le.201612.and.date.ge.197901)
    u10 = short2flt(f3[:]->u10(tind,:,:))
    slp = short2flt(f3[:]->msl(tind,:,:))
    printVarSummary(u10)     

    u10 = (/dtrend_msg_n(ispan(1,dimsizes(u10&time),1),u10,False,True,0)/)
    slp = (/dtrend_msg_n(ispan(1,dimsizes(slp&time),1),slp,False,True,0)/)
    
    u10_clm = clmMonTLL(u10)
    u10_inten = dim_avg_n_Wrap(u10_clm(:,{-5:5},{80:280}),1)
    ; u10_inten = dim_avg_n_Wrap(u10(:,{-5:5},{120:170}),(/1,2/))


    slp_clm = clmMonTLL(slp)
    slp_inten = dim_avg_n_Wrap(slp_clm(:,{-5:5},{80:280}),1)
    copy_VarCoords(slp(:,1,:), slp_inten)
     
    ;;;remove annual mean 
    slp_inten = (/dim_rmvmean_n(slp_inten, 0)/)
    u10_inten = (/dim_rmvmean_n(u10_inten, 0)/)    
    
    ;;;


    printVarSummary(u10_inten)
    printVarSummary(slp_inten)

;;;===================画图=================================
wks=gsn_open_wks("pdf", "pacific_zonal_wind_seanson_panel_pattern")
;;;
    res            = True                         ; plot mods desired
    res@gsnDraw    = False                        ; don't draw
    res@gsnFrame   = False                        ; don't advance frame
    
    res@cnLineLabelsOn =False
    res@cnFillDrawOrder="PreDraw"
    res@cnLinesOn=False
    res@cnFillOn=True
    res@lbLabelBarOn=False     
    res@cnInfoLabelOn=False
    res@cnFillPalette="GMT_polar"
  

    ; res@cnLevelSelectionMode = "ManualLevels"     ; manually select contours
    ; res@cnMinLevelValF       = -8.0              ; min level
    ; res@cnMaxLevelValF       =  8.0               ; max leve
    ; res@cnLevelSpacingF      =1.0


    res@tmYLMode      = "Explicit"                ; Define own tick mark labels.
    res@tmYLLabelsOn=True
    res@tmYUseLeft=True
    res@tmYLValues    = ispan(0,11,1)
    res@tmYLLabels    = (/"J","F","M","A","M","J","J","A","S","O","N","D"/)
    res@tmXTOn=False
    res@tmYROn=False

    res@tmXBValues=ispan(80,280,20)
    res@tmXBLabels= ""+ ispan(80,280,20)

    res@tmXBMinorOn=True
    res@tmXBMinorValues=ispan(80,280,5)

    res@vpHeightF=0.2
    res@vpWidthF=0.4
  
    ;;;LineArrow
    res@gsnScalarContour        = True
    
    res@vcGlyphStyle             = "LineArrow"      ; turn on curly vectors
    ; res@vcGlyphStyle="CurlyVector"
    ; res@vcMinDistanceF           = 0.015           ; thin vectors
    ; res@vcMinDistanceF           = 0.02           ; thin vectors
    res@vcMinMagnitudeF          =0.1
    res@vcRefAnnoOrthogonalPosF  = -1.0              ; move ref vector down -1.0在图的右上
    res@vcLineArrowThicknessF=1.0
    ; res@vcLineArrowHeadMaxSizeF=0.01
    ; res@vcLineArrowHeadMinSizeF=0.01
    res@vcVectorDrawOrder="PostDraw"

    ;;;CurlyVector
    res@vcLineArrowThicknessF=1.5
    res@vcMinFracLengthF=0.3
    res@vcMinDistanceF           = 0.015           ; thin vectors
    ; res@vcLineArrowThicknessF=1.0
    ; res@vcMinDistanceF           = 0.02          ; thin vectors
    ; res@vcMinMagnitudeF          =0.05
    ; res@vcRefAnnoOrthogonalPosF  = -1.0              ; move ref vector down -1.0在图的右上
    res@vcRefMagnitudeF          = 0.8
    res@vcRefLengthF             = 0.04     ; define length of vec ref

    v10_inten = u10_inten
    v10_inten = 0
    plot = gsn_csm_vector_scalar(wks, u10_inten,v10_inten,slp_inten, res)
 
 draw(plot)
 frame(wks)
 end

    

