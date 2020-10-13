function one_test(pv_ano,east_time)
begin
  
    pv_east=dim_avg_n_Wrap( pv_ano(east_time,:,:),0)
    n=dimsizes(east_time)
    var_pv_ano_east    =  dim_variance_n(pv_ano(east_time,:,:), 0)    ; 13.477
    var_pv_ano_east  =  sqrt(var_pv_ano_east/n) 
    var_pv_ano_east  =where(var_pv_ano_east.eq.0, var_pv_ano_east@_FillValue, var_pv_ano_east)

    t_pv_ano_east = (pv_east)/var_pv_ano_east

    confi_pv_east=student_t(t_pv_ano_east , n-1)

    copy_VarCoords(pv_east, confi_pv_east)
 return confi_pv_east
 end

;***********************************************************
function get_an_ao(sst)
begin
    ssta=(/rmMonAnnCycTLL(sst)/)
    copy_VarCoords(sst, ssta)
    ssta_de=dtrend_msg_n(ispan(1,dimsizes(ssta&time),1),ssta,False,True,0)
    copy_VarCoords(ssta, ssta_de)
    
    ssta_annual = month_to_annual(ssta_de, 1)
    ssta_annual!0 ="time"
    ssta_annual&time = ispan(1, 38, 1)
return ssta_annual
end 



;*********************************************************
;; Main code
;********************************************************

begin
;************************************************************
;;;read location dataset
;************************************************************
    ;;;read
    locationxi = asciiread("/mnt/f/share/rea_msf/era_1.5wlocationxi.ascii", 456, "float")
    locationxi = (/dtrend_msg_n(ispan(1,456,1),locationxi,False,True,0)/)
    locationxi_anu=month_to_annual(locationxi, 1)
    locationxi_anu_std=dim_standardize(locationxi_anu, 1)
    locationxi_anu_std!0="time"
    printVarSummary(locationxi_anu_std)
    locationxi_anu_std&time=ispan(1,38,1)
    ;;;提取数据
   
    east_time = ind(locationxi_anu_std.lt.-1)
    west_time=ind(locationxi_anu_std.gt.1)
    climate_mean_state=ind(locationxi_anu_std.ge.-0.5.and.locationxi_anu_std.le.0.5)

;************************************************************
;;;read SLP  U V 10-M and SST
;************************************************************
   ;;read sur data
    fils = systemfunc("ls /mnt/f/share/era/data/data_series/sur_1.5/era_interim_sur1.5_"+"*"+".nc")
    f1 = addfiles(fils, "r")
    time = f1[:]->time
    date = cd_calendar(time, 1)
    timeind=ind(date.ge.197901.and.date.le.201612)
  
    u10 = short2flt(f1[:]->u10(timeind,{30:-30},:))
    slp = short2flt(f1[:]->msl(timeind,{30:-30},:))
    slp = (/slp/100/) ;;change to hpa
    v10 = short2flt(f1[:]->v10(timeind,{30:-30},:))
    

    delete(fils)
    delete(f1)
    delete(date)
    delete(timeind)
    delete(time)
    ;;;

  ;;;read sst
    f2=addfile("/mnt/f/share/ersst/HadISST_sst.nc", "r")
    time = f2->time
    date = cd_calendar(time, 1)
    timeind=ind(date.ge.197901.and.date.le.201612)
    sst=lonFlip(f2->sst(timeind,{30:-30},:))
  ;;;
  
;************************************************************
;;;get anomaly value and phse 
;************************************************************
  
   u10_a = get_an_ao(u10)

   v10_a = get_an_ao(v10)

   slp_a = get_an_ao(slp)

   sst_a = get_an_ao(sst)
    
    ;;; get phase value
    u10_east = u10_a(east_time,:,:)
    u10_west = u10_a(west_time,:,:)
    
    v10_east = v10_a(east_time,:,:)
    v10_west = v10_a(west_time,:,:)

    slp_east = slp_a(east_time,:,:)
    slp_west = slp_a(west_time,:,:)


    sst_east = sst_a(east_time,:,:)
    sst_west = sst_a(west_time,:,:)
    
    

    v10_east_avg = dim_avg_n_Wrap(v10_east,0)
    u10_east_avg = dim_avg_n_Wrap(u10_east,0)

    v10_west_avg = dim_avg_n_Wrap(v10_west, 0)
    u10_west_avg = dim_avg_n_Wrap(u10_west, 0)

    slp_east_avg = dim_avg_n_Wrap(slp_east, 0)
    slp_west_avg = dim_avg_n_Wrap(slp_west, 0)

    sst_east_avg = dim_avg_n_Wrap(sst_east, 0)
    sst_west_avg = dim_avg_n_Wrap(sst_west, 0)


    u10_clm = dim_avg_n_Wrap(u10,0)

    v10_clm = dim_avg_n_Wrap(v10,0)
    sst_clm = dim_avg_n_Wrap(sst, 0)
    slp_clm = dim_avg_n_Wrap(slp, 0)


;;;w的t检验
    ; confi_w_east = one_test(w_latavg,east_time,w_clm_avg_ne)
    ; confi_w_west = one_test(w_latavg,west_time,w_clm_avg_ne)

    ; printVarSummary(confi_w_east)
    ; ; confi_w_east_se=confi_w_east(::2,::12)


    ; prob1d=ndtooned(confi_w_east)
    ; i1d=ind(prob1d.lt.0.05)
    ; i=ind_resolve(i1d, dimsizes(confi_w_east))

    ; ilat=confi_w_east&level(i(:,0))
    ; ilon=confi_w_east&longitude(i(:,1))
;;;

;;;u的t检验
  confi_u10_east = one_test(u10,east_time)
  confi_u10_west = one_test(u10,west_time)

  u10_east_avg = where(confi_u10_east.gt.0.05, u10_east_avg@_FillValue, u10_east_avg)
  u10_west_avg = where(confi_u10_west.gt.0.05, u10_west_avg@_FillValue, u10_west_avg)
;;;

;;; 
  wks=gsn_open_wks("pdf", "/mnt/f/share/era/data/WL指数定义的sur图")

;;--------------------------------------------------------------------------------
;;;画w和流场图
   gsn_define_colormap(wks,"hotcolr_19lev")
 ; gsn_define_colormap(wks,"WhiteBlueGreenYellowRed")
    cmap=read_colormap_file("MPL_PuOr")
    cmap=(/cmap(::-1,0:3)/)

    
    res=True
    res@gsnDraw=False
    res@gsnFrame=False
    res@gsnAddCyclic=False
 
    res@lbLabelBarOn=False
    res@cnFillOn = True 
    res@cnLinesOn =False
    res@cnLineLabelsOn = False
    res@cnInfoLabelOn =False

    res@cnFillPalette="GMT_polar"
    res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour levels
    ; res@cnFillDrawOrder="PreDraw"
    res@cnMinLevelValF       =   -1.0
    res@cnMaxLevelValF       =   1.0
    res@cnLevelSpacingF      =  0.1
    
    res@tmXTOn = False
    res@tmYROn = False
    res@tmXBLabelFontHeightF=0.012
    res@tmYLLabelFontHeightF=0.012

    res@tmXBMajorOutwardLengthF=0
    res@tmXBMajorLengthF=0.007


    res@tmYLMajorOutwardLengthF=0
    res@tmYLMajorLengthF=0.01

    res@tmYLMajorThicknessF=2.0
    res@tmXBMajorThicknessF=2.0
    res@tmBorderThicknessF=2.5

    res@tmXTMode      = "Explicit"                ; Define own tick mark labels.
    res@tmXTValues = ispan(0,360,60)
    res@tmXTLabelsOn = True
    res@tmXTLabels = (/"0E","60~S~o~N~E","120~S~o~N~E","180~S~o~N~E","120~S~o~N~W","60~S~o~N~W","0E"/)
    res@tmXTLabelDeltaF=-0.6
    

    res@tmYLMode      = "Explicit"                ; Define own tick mark labels.
    res@tmYLValues=ispan(30,-30,10)
    res@tmYLLabels=(/"30~S~o~N~N","20~S~o","10~S~o","0","10~S~o","20~S~o","30~S~o~N~S"/)
    res@tmYLLabelDeltaF=-0.6

    res@mpFillDrawOrder="PostDraw"
    res@mpFillOn=False

    res@mpLimitMode           = "LatLon"
    res@mpMinLatF             = -30
    res@mpMaxLatF             = 30
    res@mpMinLonF             = 0
    res@mpMaxLonF             = 360
    res@mpCenterLonF          = (res@mpMinLonF+res@mpMaxLonF)/2.0

    res@mpShapeMode="FreeAspect"
    res@mpFillDrawOrder  = "PostDraw"



    ;;;LineArrow
    res@gsnScalarContour        = True
    res@vcGlyphStyle             = "LineArrow"      ; turn on curly vectors   
    ; res@vcGlyphStyle = "FillArrow"
    res@vcRefMagnitudeF          = 3
    res@vcRefLengthF             = 0.05     ; define length of vec ref
    res@vcMinDistanceF           = 0.02           ; thin vectors
    res@vcMinMagnitudeF = 0.5


    res@vcRefAnnoOrthogonalPosF  = -1.0              ; move ref vector down -1.0在图的右上
    res@vcLineArrowThicknessF = 1.0
    res@vcLineArrowHeadMaxSizeF=0.01
    res@vcLineArrowHeadMinSizeF=0.01
    res@vcVectorDrawOrder="PostDraw"

    res@vcRefAnnoString1="3m/s"
    res@vcRefAnnoOn = False
    res@vcRefAnnoFontHeightF=0.01

    res@vpHeightF = 0.4
    res@vpWidthF = 0.6
   
    ;---draw plot 
        
  
    plot=new(3,graphic)
    res@gsnLeftString="(a)West phase"
    plot(0) = gsn_csm_vector_scalar_map(wks,u10_east_avg,v10_east_avg,sst_east_avg,res)
    res@gsnLeftString="(b)West phase"
    plot(1) = gsn_csm_vector_scalar_map(wks,u10_west_avg,v10_west_avg,sst_west_avg,res)
    res@gsnLeftString="(c)climate_mean_state"
    
    res@cnMinLevelValF = 20.0
    res@cnMaxLevelValF = 29.0
    res@cnLevelSpacingF = 1.0
 
    ; res@vcRefMagnitudeF          = 10
    delete(res@vcRefMagnitudeF)
    res@vcMinDistanceF           = 0.04           ; thin vectors
    res@vcRefLengthF             = 0.04     ; define length of vec ref
    plot(2)  = gsn_csm_vector_scalar_map(wks,u10_clm,v10_clm,sst_clm,res)
    ;;;;

    ;;;;; 画slp 等值线

        res_slp=True
        res_slp@gsnDraw=False
        res_slp@gsnFrame=False
        res_slp@gsnLeftString = ""
        res_slp@gsnRightString = ""
   
        res_slp@cnLinesOn=True
        res_slp@cnLevelSelectionMode = "ManualLevels" ; set manual contour levels
        res_slp@cnMinLevelValF      = -1.0        ; set min contour level
        res_slp@cnMaxLevelValF      = 1.0      ; set max contour level
        res_slp@cnLevelSpacingF     = 0.2      ; set contour spacin

        res_slp@cnLineLabelsOn     = True
        res_slp@cnLineLabelBackgroundColor="Transparent"

        
        colors=new(11,string)
        ; colors(0:3)="darkgoldenrod"
        ; colors(4)="black"
        ; colors(5:)="royalblue2"
         
        colors(0:4)="blueviolet"
        ; colors(0:4)="springgreen"
        ; colors(0:4)="darkgreen"
        
        ; colors(5)="sandybrown"
        colors(5)="black"
        ; colors(5)="white"
        ; colors(6:)="seagreen4"
        ; colors(6:)="springgreen"
        colors(6:)="orange"
        ; colors(6:)="deepskyblue"

        res_slp@cnLineColors=colors
        res_slp@cnMonoLineColor=False
        res_slp@gsnContourNegLineDashPattern=12
        res_slp@cnLineLabelInterval        = 1             ; default = 2
     
        res_slp@cnLineLabelFormat= "0@;*.2f" 
        res_slp@cnInfoLabelOn=False
        res_slp@cnLineLabelFont=3.0
        res_slp@cnLineLabelFontHeightF=0.008
        res_slp@cnInfoLabelOn=False
        res_slp@cnLineThicknessF=2.5
        res_slp@gsnContourZeroLineThicknessF=2.0

        plot_slp=new(3,graphic)
  
       
        plot_slp(0)=gsn_csm_contour(wks,slp_east_avg ,res_slp)
        
        plot_slp(1)=gsn_csm_contour(wks,slp_west_avg,res_slp)  
         
        
     ; set contour spacin
        delete(res_slp@cnLineLabelFormat)  
        res_slp@cnMonoLineColor=True
        res_slp@cnLineColor = "black"
        res_slp@cnLevelSelectionMode ="ManualLevels"
        res_slp@cnMinLevelValF = 1008.0
        res_slp@cnMaxLevelValF = 1021.0
        res_slp@cnLevelSpacingF = 2.0
        plot_slp(2)=gsn_csm_contour(wks,slp_clm,res_slp)
        
        do i=0,2
        overlay(plot(i), plot_slp(i))
        end do


;;;panel
    pres=True
    pres@gsnMaximize = True
    pres@gsnPanelLabelBar=True
    pres@lbOrientation ="Vertical"
    ; pres@pmLabelBarWidthF=0.5
    gsn_panel(wks, plot, (/1,3/), pres)

  end
