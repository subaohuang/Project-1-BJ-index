function level_mean(msf_time)
begin
    plev=msf_time&level
    psfc=conform_dims((/dimsizes(msf_time&time),dimsizes(msf_time&lon)/), 1013.0, 0)
    ptop=min(plev)
    dp   = dpres_plevel_Wrap(plev, psfc, ptop, 0)  ; dp(30)
    msfdp=msf_time(level|:,time|:,lon|:)*dp
    copy_VarCoords(msf_time(level|:,time|:,lon|:), msfdp)

    msf_mean=dim_sum_n_Wrap(msfdp, 0)/dim_sum_n_Wrap(dp, 0)
    copy_VarCoords(msf_time(time|:,lon|:,level|1), msf_mean)
return(msf_mean)
end 

function  search_location(msf_mean_annual)  ;;;msf_mean_annual(time,lon)
begin
    locationxi=new(dimsizes(msf_mean_annual&time),typeof(msf_mean_annual&lon))
    locationxi!0="time"
    locationxi&time=ispan(0,dimsizes(msf_mean_annual&time)-1,1)
    locationdong=locationxi
    msf_mean_annual_xi=msf_mean_annual(:,{120:220})
    msf_mean_annual_dong=msf_mean_annual(:,{230:280})
;;;judge locationxi
    do t=0,dimsizes(msf_mean_annual&time)-1
        locationxi(t)=msf_mean_annual_xi&lon(maxind(msf_mean_annual_xi(t,:)))
        locationdong(t)=msf_mean_annual_dong&lon(minind(msf_mean_annual_dong(t,:)))
     end do
     location = [/ locationxi, locationdong /]   ; return two variables as a type 'list'
return location
end 

function draw_msf_hov(fils1,wks,LeftString,season)
begin
    diri="/mnt/f/share/"
    ; diri_pi="/mnt/f/share/era/data/"
    ; fils1="/mnt/f/share/ncep/ncep1/ncep1_m_zonal_mer.nc"
   
    f1=addfile(fils1, "r")
  
    date=cd_calendar(f1->time, 1)
    timeind=ind(date.ge.195801.and.date.le.201212)
    ; timeind=ind(date.ge.197901.and.date.le.201212)
    m_zonal_r=f1->m_zonal(timeind,:,{5:-5},:)
    ; m_zonal_r=f1->omega(timeind,:,{-5:5},:)
    printVarSummary(m_zonal_r)
    
    m_zonal=dim_avg_n_Wrap(m_zonal_r, 2)
    ; m_zonal=(/-1*m_zonal/)
    msf_mean=m_zonal(:,{500},:)

    msf_mean_annual=month_to_annual(msf_mean, 1)
    printVarSummary(msf_mean_annual)
    msf_mean_annual!0="time"
    msf_mean_annual&time=ispan(0,dimsizes(timeind)/12-1,1)


    location=search_location(msf_mean_annual)
    locationxi=location[0]
    locationdong=location[1]

    rc_xi=regCoef_n(ispan(0,dimsizes(locationxi)-1,1),locationxi , 0, 0)
    rc_xi_line=rc_xi*ispan(0,dimsizes(locationxi)-1,1) + rc_xi@yintercept

    rc_dong=regCoef_n(ispan(0,dimsizes(locationdong)-1,1),locationdong, 0, 0)
    rc_dong_line=rc_dong*ispan(0,dimsizes(locationdong)-1,1) + rc_dong@yintercept

    width=locationdong - locationxi

    width=(/runave(width, 3, -1)/)
    width=(/dim_standardize_n(width, 1, 0)/)
    time=locationxi&time

    rc_width=regCoef_n(ispan(0,dimsizes(locationxi)-1,1),width , 0, 0)
    rc_width_line=rc_width*ispan(0,dimsizes(locationxi)-1,1) + rc_width@yintercept

;***********************************************
;draw hov
;***********************************************
   
    ; cmap=read_colormap_file("MPL_PuOr")
    ; cmap=(/cmap(::-1,0:3)/)
  
	res            = True                         ; plot mods desired
	res@gsnDraw    = False                        ; don't draw
	res@gsnFrame   = False                        ; don't advance frame

	; res@gsnContourZeroLineThicknessF=4.0
 ;    res@gsnContourNegLineDashPattern=1
 ;    res@gsnContourPosLineDashPattern=0

    res@gsnLeftString="Ncep1"
	res@cnLineLabelsOn =False
	res@cnFillDrawOrder="PreDraw"
	res@cnLinesOn=False
    ; res@cnFillPalette="BlueDarkRed18"
    res@cnFillPalette="testcmap"
	res@cnFillOn=True
    res@lbLabelBarOn=True
    res@cnInfoLabelOn=False
	; res@cnLevelSelectionMode = "ManualLevels"     ; manually select contours
	; res@cnMinLevelValF       = -0.0028               ; min level
	; res@cnMaxLevelValF       =  0.0028               ; max leve
	; res@cnLevelSpacingF      =0.0004
	res@lbOrientation         ="Vertical"
    res@lbLabelBarOn=False

	res@tmYLMode      = "Explicit"                ; Define own tick mark labels.
	res@tmYLLabelsOn=True
	res@tmYUseLeft=True

	res@tmYLValues    = ispan(2,dimsizes(timeind)/12-1,5)    
	res@tmYLLabels    = ispan(1960,2010,5)
    ; res@tmYLValues    = ispan(1,dimsizes(timeind)/12-1,5)
    ; res@tmYLLabels=ispan(1980, 2010, 5)

	res@tmXTOn=False
	res@tmYROn=False

	res@tmXBValues=(/0,90,180,270,357.5/)
    res@tmXBLabels=(/"0","90~S~o~N~E","180","90~S~o~N~W","0"/)
    res@tmXBMinorOn=True
    res@tmXBMinorValues=ispan(0,360,10)
	res@tmXBMinorOn=True

	res@vpHeightF=0.25
	res@vpWidthF=0.4
    plot = gsn_csm_contour(wks, msf_mean_annual, res)   ; create plot
    

    res@cnLevelSelectionMode = "ExplicitLevels"
    res@cnLinesOn=True
    res@cnLevels = 0
    res@cnLineThicknessF = 4.
    res@cnFillOn        = False

    plot1 = gsn_csm_contour(wks, msf_mean_annual, res)
    
    overlay(plot, plot1)

    ;;;;;add trend
    lineres=True
    lineres@gsLineColor="red"
    lineres@gsLineDashPattern=1
    lineres@gsLineThicknessF=3.0
   
    str1 = unique_string("dum1")  ; "unique_string" will return a unique
                                     ; string every time it is called from
                                     ;  within a single NCL session.
    str2 = unique_string("dum2")
    str3 = unique_string("dum3")
    str4 = unique_string("dum4")
    plot@$str1$=gsn_add_polyline(wks, plot, (/rc_xi_line/), ispan(0,54,1), lineres)
    plot@$str2$=gsn_add_polyline(wks, plot, (/rc_dong_line/), ispan(0,54,1), lineres)
    
    plot@$str3$=gsn_add_polyline(wks, plot, (/locationxi/), ispan(0,54,1), lineres)
    plot@$str4$=gsn_add_polyline(wks, plot, (/locationdong/), ispan(0,54,1), lineres)
    ;;;;draw_attch_plot***************************************************************************
    resxy=True
    resxy@gsnDraw=False
    resxy@gsnFrame=False
    resxy@tmXBMode="Explicit"
    resxy@tmXTOn=False
    resxy@tmYLOn=False
    resxy@tmYROn=False
    resxy@trYMaxF=dimsizes(timeind)/12-1
    resxy@trXMinF=-2.0
    resxy@trXMaxF=3.0
    resxy@tmXBValues=ispan(-1,3,1)
    resxy@tmXBLabels=ispan(-1,3,1)
    resxy@tmXBMinorValues=fspan(-2,3,26 )
    resxy@tmXBMinorOn=True

    resxy@tmBorderThicknessF=3.0     ;调整边框粗细
    resxy@tmYUseLeft=False
 
    resxy@xyLineColors      = (/"blue","red"/)          ; change line color
    resxy@xyMarkLineModes   = (/"lines","lines"/)
    resxy@xyDashPatterns    =(/0,0/)
    resxy@xyLineThicknesses=(/2.0,2.0/)
    
    resxy@vpHeightF=0.25
    resxy@vpWidthF=0.1
    plot_line=gsn_csm_xy(wks, (/width,rc_width_line/),locationxi&time, resxy)

    pres1=True
    pres2=True
    plot3= gsn_attach_plots(plot,plot_line,pres1,pres2)
 
    return(plot)
    end  
;;**************************************************************
;; main program
;;;***************************************************************
 begin 
    fils1="ncep/20CR/20cr_m_zonal_meri.nc"
    fils2="ncep/ncep1/ncep1_m_zonal_mer.nc"
    ; fils3="era/era40/msf_time_era40_in_global.nc"
    ; fils4="JRA-55/msf_time_jra55_global.nc"

    ; wks=gsn_open_wks("pdf", "/mnt/f/share/era/data/4个再分析资料(1979-2012)计算的msf_hov")
    wks=gsn_open_wks("pdf", "/mnt/f/share/era/data/500hpa&4个再分析资料计算的local_wo_hov")
    plot=new(2,graphic)
    season=0
    plot(0)=draw_msf_hov(fils1,wks,"20CR",season)   ;1 means winter 0 means aunnal 2 means summer
    plot(1)=draw_msf_hov(fils2,wks,"NCPE1",season)
    ; plot(2)=draw_msf_hov(fils3,wks,"ERA",season)
    ; plot(3)=draw_msf_hov(fils4,wks,"JRA-55",season)  

    pres=True
    pres@gsnMaximize=True
    pres@gsnPanelLabelBar=True
    pres@pmLabelBarWidthF=0.8
    pres@pmLabelBarHeightF=0.05
    pres@pmLabelBarOrthogonalPosF =0.01
    pres@gsnPanelDebug=True
    gsn_panel(wks, plot, (/2,1/), pres)
 
    end



