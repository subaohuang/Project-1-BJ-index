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

function  sst_cold_annual(sst,timeind)
begin 
    year1=1979
    year2=2005
    ;;;先得到秋季和冬季的SONDJF海温
    sst_son=month_to_season(sst(time|:,lat|:,lon|:), "SON")
    sst_djf=month_to_season(sst(time|:,lat|:,lon|:), "DJF")
    printVarSummary(sst_son)

    sst_son&time=ispan(year1,year2,1)
    sst_djf&time=ispan(year1,year2,1)  ;;;因为这里第一个是1958 1月和2月的平均，第二个是195812,195901,195902的平均
    ;;注意点
    sst_djf(0:dimsizes(sst_djf&time)-2,:,:)=(/sst_djf(1:,:,:)/)
    sst_djf(dimsizes(sst_djf&time)-1,:,:)=(/sst({200512},:,:)/)
    
    ;;;
    sst_mean_annual=new((/dimsizes(timeind)/12,dimsizes(sst&lat),dimsizes(sst&lon)/),typeof(sst))
    printVarSummary(sst_mean_annual)
    ; sst_mean_annual=month_to_annual(sst_mean, 1)
     ; sst_mean_annual=(/sst_son(:,:,:)+sst_djf(:,:,:)/)*0.5
     sst_mean_annual=(/sst_djf(:,:,:)/)*1.0
    delete(sst_son)
    delete(sst_djf)

    sst_mean_annual!0="time"
    sst_mean_annual!1="lat"
    sst_mean_annual!2="lon"
    sst_mean_annual&time=ispan(0,dimsizes(timeind)/12-1,1)
    sst_mean_annual&lon=sst&lon
    sst_mean_annual&lat=sst&lat
    printVarSummary(sst_mean_annual)
return sst_mean_annual
end

function sst_cha_trend(sst,year1,year2)
begin
 sst_cha= dim_avg_n_Wrap(sst(:,{-5:5},{80:160}),(/1,2/)) - dim_avg_n_Wrap(sst(:,{-5:5},{200:280}),(/1,2/))
    sst_cha!0="time"
    sst_cha&time=ispan(year1,year2,1)
    printVarSummary(sst_cha)

    rc_sst=regCoef_n(ispan(0,dimsizes(sst_cha&time)-1,1),sst_cha , 0, 0)
    rc_sst=(/rc_sst*10/)
return rc_sst
end 

 begin
;;读入数据
    print("*****************HadISST*********************")
    yearstart=197901
    yearend=200512
    
    year1=1979
    year2=2005
    
    latmin=-30
    latmax=30
    lonmin=90
    lonmax=300

    fils1=systemfunc("ls /mnt/f/share/ersst/HadISST_sst.nc")
    f1=addfile(fils1, "r")
    

    vNames= getfilevarnames(f1) ; get variable names of all groups on file
    print(vNames)
    
    date=cd_calendar(f1->time, 1)
    print("******************日期*********************")
    timeind=ind(date.le.yearend.and.date.ge.yearstart)

    sst_r_had=lonFlip(f1->sst(timeind,:,:))
    dNames=getVarDimNames(sst_r_had)
  
    if any(dNames.eq."latitude") then 
    sst_r_had!0="time"
    sst_r_had!1="lat"
    sst_r_had!2="lon"
    end if 
    

    printVarSummary(sst_r_had)
    
    delete(sst_r_had&time)
    sst_r_had&time=date(timeind)
    sst_had=sst_cold_annual(sst_r_had(:,{latmax:latmin},{lonmin:lonmax}),timeind)

    delete(sst_r_had)
    delete(date)
    delete(f1)
    delete(fils1)
    
 print("*****************ERSST*********************")
    ;;;读入ERSST海温数据
    m=""+ispan(1979,2005,1)+""
    diri="/mnt/f/share/ersst/"

    fils=new((/27,12/),string)
    do i=0,26
    fils1=systemfunc("ls "+diri+"ersst."+m(i)+"*"+".nc")
    fils(i,:)= fils1
    end do
    fils_new=onedtond(fils, 27*12)
    fileh=new(324,string)
    fileh(0:323)=fils_new

    f1=addfiles(fileh, "join")
    printVarSummary(f1)
    date=cd_calendar(f1[:]->time, 1)
    print("******************日期*********************")
    timeind=ind(date.le.yearend.and.date.ge.yearstart)

    sst_r_ersst=f1[:]->sst(timeind,0,:,:)
    ; sst_r=(/sst_r-273.15/)
    printVarSummary(sst_r_ersst)

    delete(sst_r_ersst&time)
    sst_r_ersst&time=date(timeind)
    sst_ersst=sst_cold_annual(sst_r_ersst(:,{latmax:latmin},{lonmin:lonmax}),timeind)

    delete(sst_r_ersst)
    delete(date)
    delete(f1)
    delete(fils1)

 print("*****************COBE*********************")

    fils1=systemfunc("ls /mnt/f/share/ersst/COBE_SST.nc")
    f1=addfile(fils1, "r")

    date=cd_calendar(f1->time, 1)
    print("******************日期*********************")
    timeind=ind(date.le.yearend.and.date.ge.yearstart)

    sst_r_cobe=f1->sst(timeind,:,:)
    ; sst_r=(/sst_r-273.15/)
    printVarSummary(sst_r_cobe)
    
    delete(sst_r_cobe&time)
    sst_r_cobe&time=date(timeind)
    sst_cobe=sst_cold_annual(sst_r_cobe(:,{latmax:latmin},{lonmin:lonmax}),timeind)

    delete(sst_r_cobe)
    delete(date)
    delete(f1)
    delete(fils1)
   
rc_hadisst=sst_cha_trend(sst_had,year1,year2)
rc_ersst=sst_cha_trend(sst_ersst,year1,year2)
rc_cobe=sst_cha_trend(sst_cobe,year1,year2)
rc_mean=(/rc_hadisst+rc_ersst+rc_cobe/)/3.0


;;;;
   print("*****************读入MRE的风场和海表气压场*********************")
   
    fils1=systemfunc("ls /mnt/f/share/MRE/mre_u_v_chazhi_1000.nc")
    f1=addfile(fils1, "r")
    date=cd_calendar(f1->time, 1)
    print("******************日期*********************")
    timeind=ind(date.le.yearend.and.date.ge.yearstart)
    
    uwnd_mre_r=f1->uwnd(timeind,:,:)
    vwnd_mre_r=f1->vwnd(timeind,:,:)

    delete(uwnd_mre_r&time)
    delete(vwnd_mre_r&time)
    uwnd_mre_r&time=date(timeind)
    vwnd_mre_r&time=date(timeind)
    uwnd_mre=sst_cold_annual(uwnd_mre_r(:,{latmax:latmin},{lonmin:lonmax}),timeind)
    vwnd_mre=sst_cold_annual(vwnd_mre_r(:,{latmax:latmin},{lonmin:lonmax}),timeind)
     
    uwnd_mre&lat@units="degree_north"
    uwnd_mre&lon@units="degree_east"
    vwnd_mre&lat@units="degree_north"
    vwnd_mre&lon@units="degree_east"

    delete(uwnd_mre_r)
    delete(vwnd_mre_r)
    delete(date)
    delete(f1)
    delete(fils1)

  
    ;;;气压场
    fils1=systemfunc("ls /mnt/f/share/MRE/mre_slp_chazhi.nc")
    f1=addfile(fils1, "r")
    date=cd_calendar(f1->time, 1)
    print("******************日期*********************")
    timeind=ind(date.le.yearend.and.date.ge.yearstart)
    
    slp_r=f1->slp(timeind,:,:)

    delete(slp_r&time)
    slp_r&time=date(timeind)
    slp=sst_cold_annual(slp_r(:,{latmax:latmin},{lonmin:lonmax}),timeind)
      
    slp&lat@units="degree_north"
    slp&lon@units="degree_east"
    delete(slp_r)
    delete(date)
    delete(f1)
    delete(fils1)
;;;


;;;
    print("*****************计算趋势************************")  

    sst_rc1=regCoef_n(ispan(1,dimsizes(sst_had&time),1),sst_had(:,:,:),0,0)  
    copy_VarCoords(sst_had(0,:,:), sst_rc1)
    printVarSummary(sst_rc1)
    
    prob1=rc_test(sst_rc1)
    copy_VarCoords(sst_rc1, prob1)
    printVarSummary(prob1)
    
    sst_rc1=(/sst_rc1/)*10
    
  

    sst_rc2=regCoef_n(ispan(1,dimsizes(sst_ersst&time),1),sst_ersst(:,:,:),0,0)  
    copy_VarCoords(sst_ersst(0,:,:), sst_rc2)
    printVarSummary(sst_rc2)
  
    prob2=rc_test(sst_rc2)
    copy_VarCoords(sst_rc2, prob2)
    printVarSummary(prob2)
   
    sst_rc2=(/sst_rc2/)*10


    sst_rc3=regCoef_n(ispan(1,dimsizes(sst_cobe&time),1),sst_cobe(:,:,:),0,0)  
    copy_VarCoords(sst_cobe(0,:,:), sst_rc3)
    printVarSummary(sst_rc3)

    prob3=rc_test(sst_rc3)
    copy_VarCoords(sst_rc3, prob3)
    printVarSummary(prob3)

    sst_rc3=(/sst_rc3/)*10


    uwnd_rc1=regCoef_n(ispan(1,dimsizes(uwnd_mre&time),1),uwnd_mre(:,:,:),0,0)  
    copy_VarCoords(uwnd_mre(0,:,:), uwnd_rc1)
    printVarSummary(uwnd_rc1)
  

    prob4=rc_test(uwnd_rc1)
    copy_VarCoords(uwnd_rc1, prob4)
    printVarSummary(prob4)
    

    uwnd_rc1=(/uwnd_rc1/)*10
    
    ;;;************************************
      
    vwnd_rc1=regCoef_n(ispan(1,dimsizes(vwnd_mre&time),1),vwnd_mre(:,:,:),0,0)  
    copy_VarCoords(vwnd_mre(0,:,:), vwnd_rc1)
    printVarSummary(vwnd_rc1)
    
    prob5=rc_test(vwnd_rc1)
    copy_VarCoords(vwnd_rc1, prob5)
    printVarSummary(prob5)
    vwnd_rc1=(/vwnd_rc1/)*10

    uwnd_rc1=where(prob4.le.0.1.or.prob5.le.0.1, uwnd_rc1, uwnd_rc1@_FillValue)


    ;;;***********************************
    slp_rc1=regCoef_n(ispan(1,dimsizes(slp&time),1),slp(:,:,:),0,0)  
    copy_VarCoords(slp(0,:,:), slp_rc1)
    printVarSummary(slp_rc1)
    slp_rc1=(/slp_rc1/)*10

  

    
    print("**********************画图*****************************")

    wks=gsn_open_wks("pdf", "/mnt/f/share/CMIP5/tos_group/观测的sst趋势变化")
    ; colormap=read_colormap_file("BlWhRe")
    ; printVarSummary(colormap)
    ; cmap=new((/18,4/),typeof(colormap))
    ; printVarSummary(colormap)
    ; cmap(0:4,:)=colormap(0:40:10,:)
    ; cmap(5:17,:)=colormap(52:100:4,:)
    colormap=read_colormap_file("GMT_jet")
    printVarSummary(colormap)
    cmap=new((/20,4/),typeof(colormap))
    printVarSummary(colormap)
    cmap(0:6,:)=colormap(2:92:15,:)
    cmap(7:19,:)=colormap(106:250:12,:)


    res=True
    res@gsnDraw=False
    res@gsnFrame=False
    res@gsnAddCyclic=False
 
   
    res@gsnCenterString ="HadISST"
    res@gsnCenterStringFontHeightF=0.02
    res@cnLineColors="Transparent"  ;
    res@lbLabelBarOn=False
    res@cnFillOn = True 
    res@cnLinesOn =False
    res@cnLineLabelsOn=False

    res@lbLabelBarOn=False
  
    res@cnFillPalette=cmap
    res@cnLevelSelectionMode = "ManualLevels"      
    res@cnMinLevelValF       =  -0.3
    res@cnMaxLevelValF       =   0.6
    res@cnLevelSpacingF      =  0.05

    res@tmXTOn=False
    res@tmYROn=False
  
    res@mpFillDrawOrder="PostDraw"
    res@mpLimitMode           = "LatLon"
    res@mpMinLatF             = latmin
    res@mpMaxLatF             = latmax
    res@mpMinLonF             = lonmin
    res@mpMaxLonF             = lonmax
    res@mpCenterLonF          = (res@mpMinLonF+res@mpMaxLonF)/2.0

    res@gsnScalarContour        = True
    res@vcRefMagnitudeF          = 0.5 
    res@vcRefLengthF             = 0.03     ; define length of vec ref
    res@vcGlyphStyle             = "LineArrow"      ; turn on curly vectors
    res@vcMinDistanceF           = 0.015           ; thin vectors
    res@vcRefAnnoOrthogonalPosF  = -1.0              ; move ref vector down -1.0在图的右上
   
    res@vcLineArrowThicknessF=2.0
    res@vcLineArrowHeadMaxSizeF=0.01
    res@vcLineArrowHeadMinSizeF=0.01
    res@mpShapeMode="FreeAspect"
    res@vpWidthF = 0.7
    res@vpHeightF = 0.3
    plot=new(3,graphic)

    plot(0)= gsn_csm_vector_scalar_map(wks,uwnd_rc1,vwnd_rc1,sst_rc1,res)
 
    res@gsnCenterString ="ERSST"
    plot(1)= gsn_csm_vector_scalar_map(wks,uwnd_rc1,vwnd_rc1,sst_rc2,res)
    res@gsnCenterString ="COBE_SST"
    plot(2)= gsn_csm_vector_scalar_map(wks,uwnd_rc1,vwnd_rc1,sst_rc3,res)


         ;;;打点
        res_pro=True
        res_pro@gsnDraw=False
        res_pro@gsnFrame=False
        res_pro@cnLevelSelectionMode = "ExplicitLevels"
        res_pro@cnLinesOn=False
        res_pro@cnLevelSelectionMode = "ManualLevels" ; set manual contour levels
        res_pro@cnMinLevelValF      = 0.00        ; set min contour level
        res_pro@cnMaxLevelValF      = 1.00        ; set max contour level
        res_pro@cnLevelSpacingF     = 0.05       ; set contour spacing

        res_pro@cnLineLabelsOn=False
        res_pro@tmYRMode = "Automatic"
        res_pro@cnInfoLabelOn=False

        iplot=new(3,graphic)

        iplot(0)=gsn_csm_contour(wks,prob1,res_pro)
        iplot(1)=gsn_csm_contour(wks,prob2,res_pro)
        iplot(2)=gsn_csm_contour(wks,prob3,res_pro)
     
        opt                      = True
        opt@gsnShadeFillType     = "pattern"
        opt@gsnShadeMid     = 17         ; stipple
        ; opt@gsnShadeHigh   = "green"
        opt@gsnShadeFillScaleF   = 0.2         ; add extra density
        opt@gsnShadeFillDotSizeF = 0.0      ; make dots larger
        do i=0,2
        iplot(i) = gsn_contour_shade(iplot(i), 0,0.05, opt)
    end do
        ;;;
        do i=0,2  
        overlay(plot(i), iplot(i))
    end do

        ;;;画slp的等值线
        res_slp=True
        res_slp@gsnDraw=False
        res_slp@gsnFrame=False
   
        res_slp@cnLinesOn=True
        res_slp@cnLevelSelectionMode = "ManualLevels" ; set manual contour levels
        res_slp@cnMinLevelValF      = -0.8        ; set min contour level
        res_slp@cnMaxLevelValF      = 0.8        ; set max contour level
        res_slp@cnLevelSpacingF     = 0.2       ; set contour spacing
        res_slp@cnLineLabelsOn     =True
        res_slp@cnLineLabelBackgroundColor="Transparent"
        res_slp@gsnContourNegLineDashPattern=2
        res_slp@cnLineLabelInterval        = 1             ; default = 2
        res_slp@cnLineLabelPlacementMode   = "constant"    ; choose constant label method
        res_slp@cnLineLabelFormat= "0@;*.1f" 
        res_slp@cnLineLabelFont=3.0
        res_slp@cnLineLabelFontHeightF=0.008
        res_slp@cnInfoLabelOn=False
        res_slp@cnLineColor="black"
        res_slp@cnLineThicknessF=1.0
    plot_slp=new(3,graphic)
    do i=0,2
     plot_slp(i)=gsn_csm_contour(wks,slp_rc1,res_slp)   
 end do
        do i=0,2
        overlay(plot(i), plot_slp(i))
        end do
       pres=True
       pres@gsnPanelLabelBar=True

       gsn_panel(wks, plot, (/3,1/), pres)
    end


  
