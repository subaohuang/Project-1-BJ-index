function level_mean(msf_time)
begin

plev=msf_time&level
psfc=conform_dims((/dimsizes(msf_time&time),dimsizes(msf_time&lon)/), 1013.0, 0)
ptop=min(plev)
dp   = dpres_plevel_Wrap(plev, psfc, ptop, 0)  ; dp(30)
msfdp=msf_time(level|:,time|:,lon|:)*dp
copy_VarCoords(msf_time(level|:,time|:,lon|:), msfdp)
print("dp="+dp(:,1,:))

msf_mean=dim_sum_n_Wrap(msfdp, 0)/dim_sum_n_Wrap(dp, 0)
copy_VarCoords(msf_time(time|:,lon|:,level|1), msf_mean)
return(msf_mean)
end 

function  search_location(msf_mean_annual)
begin

locationxi=new(dimsizes(msf_mean_annual&time),typeof(msf_mean_annual&lon))
locationxi!0="time"
locationxi&time=ispan(0,dimsizes(msf_mean_annual&time)-1,1)
locationdong=locationxi
msf_mean_annual_xi=msf_mean_annual(:,{120:220})
msf_mean_annual_dong=msf_mean_annual(:,{230:300})
;;;judge locationxi
do t=0,dimsizes(msf_mean_annual&time)-1
    do i=0,dimsizes(msf_mean_annual_xi&lon)-13
        if  (msf_mean_annual_xi(t,i+1).ge.0).and.((msf_mean_annual_xi(t,i).le.0).and.(msf_mean_annual_xi(t,i+12).gt.0))  then 
        locationxi(t)=(msf_mean_annual_xi&lon(i)+msf_mean_annual_xi&lon(i+1))/2
        end if    
        if ismissing(locationxi(t)).ne.True  then 
        break
        end if
     end do
end do
;;;judge locationdong
do t=0,dimsizes(msf_mean_annual&time)-1
    do i=0,dimsizes(msf_mean_annual_dong&lon)-2
        if  (msf_mean_annual_dong(t,i+1).le.0).and.((msf_mean_annual_dong(t,i).ge.0).and.(msf_mean_annual_dong(t,i+3).lt.0))  then 
         ; if  (msf_mean_annual_dong(t,i+1).le.0).and.((msf_mean_annual_dong(t,i).ge.0))  then 
        locationdong(t)=(msf_mean_annual_dong&lon(i)+msf_mean_annual_dong&lon(i+1))/2
        end if    
        if ismissing(locationdong(t)).ne.True  then 
        break
        end if
     end do
end do
; locationxi=(/linmsg(locationxi, -1)/)
; locationdong=(/linmsg(locationdong, -1)/)

 location = [/ locationxi, locationdong /]   ; return two variables as a type 'list'
return location
end 


  undef("fit") ;;;得到回归的值
  function fit(slp_20cr_mean,ipindex1)
  begin
  rc_slp_20cr=regCoef_n(ipindex1, slp_20cr_mean ,0, 0)
  slp_20cr_re=new(dimsizes(slp_20cr_mean),typeof(slp_20cr_mean))
  do i=0,dimsizes(slp_20cr_mean&lon)-1,1
  slp_20cr_re(:,i)=rc_slp_20cr(i)*ipindex1+rc_slp_20cr@yintercept(i)
  end do 
  copy_VarCoords(slp_20cr_mean, slp_20cr_re)
  return(slp_20cr_re)
  end

;********************************************************
;draw_msf_hov
;*******************************************************
function draw_msf_hov(fils1,wks,LeftString,mode)
begin
diri="/mnt/f/share/"
; diri_pi="/mnt/f/share/era/data/"

f1=addfile(fils1, "r")
date=cd_calendar(f1->time, 1)
timeind=ind(date.ge.195801.and.date.le.201212)

; msf_time=f1->msf_time({400:200},:,timeind)

msf_time=f1->msf_time({700},:,timeind)   ;;;测试cmip5数据时，取了700
msf_time=(/msf_time/10^9/)
; msf_mean=level_mean(msf_time)
msf_mean=msf_time(time|:,lon|:)
; msf_mean=dim_avg_n_Wrap(msf_time(time|:,level|:,lon|:), 1)
msf_mean_annual=month_to_annual(msf_mean, 1)
printVarSummary(msf_mean_annual)
msf_mean_annual!0="time"
msf_mean_annual&time=ispan(0,54,1)
; msf_mean_annual=(/runave_n(msf_mean_annual, 3,-1 , 0)/)
;;;
ipindex1=asciiread("/mnt/f/share/era/data/GHCN_HAD_p_ud_timef1_1958_2012.ascii", 55, "float")
ipindex2=asciiread("/mnt/f/share/era/data/GHCN_HAD_p_ud_timef2_1958_2012.ascii", 55, "float")
;;;
  ;;;;回归
    if mode.eq.0 then 
    msf_mean_annual=(/fit(msf_mean_annual,ipindex1)/)
    end if
    if mode.eq.1 then
    msf_mean_annual=(/fit(msf_mean_annual,ipindex2)/)
    end if 
    if mode.eq.3
    msf_mean_annual=((/fit(msf_mean_annual,ipindex1)/)+(/fit(msf_mean_annual,ipindex2)/))/2.0   ;;;这里有点问题
    end if
location=search_location(msf_mean_annual)
locationxi=location[0]
locationdong=location[1]

rc_xi=regCoef_n(ispan(0,dimsizes(locationxi)-1,1),locationxi , 0, 0)
rc_xi_line=rc_xi*ispan(0,dimsizes(locationxi)-1,1) + rc_xi@yintercept
rc_xi=(/rc_xi*10/)
rc_dong=regCoef_n(ispan(0,dimsizes(locationdong)-1,1),locationdong, 0, 0)
rc_dong_line=rc_dong*ispan(0,dimsizes(locationdong)-1,1) + rc_dong@yintercept
rc_dong=(/rc_dong*10/)
width=locationdong - locationxi

  
    width=(/dim_rmvmean_n(width,0)/)
    ; width=(/dim_standardize_n(width, 1, 0)/)
    width=(/runave(width, 3, -1)/)
    time=locationxi&time

    rc_width=regCoef_n(ispan(0,dimsizes(locationxi)-1,1),width , 0, 0)
    rc_width_line=rc_width*ispan(0,dimsizes(locationxi)-1,1) + rc_width@yintercept
    rc_width=(/rc_width*10/)
  
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

    res@gsnLeftString=LeftString
	res@cnLineLabelsOn =False
	res@cnFillDrawOrder="PreDraw"
	res@cnLinesOn=False
    res@gsnRightString="rc_xi:"+sprintf("%2.2f",rc_xi)+"~S~o~N~/dacade"+"&rc_dong:"+sprintf("%2.2f",rc_dong)+"~S~o~N~/dacade"
    res@cnFillPalette="testcmap"
	res@cnFillOn=True
    res@lbLabelBarOn=False     
    res@cnInfoLabelOn=False
	res@cnLevelSelectionMode = "ManualLevels"     ; manually select contours
	res@cnMinLevelValF       = -7.0               ; min level
	res@cnMaxLevelValF       =  7.0               ; max leve
	res@cnLevelSpacingF      =0.5
	res@lbOrientation         ="Vertical"
   
	res@tmYLMode      = "Explicit"                ; Define own tick mark labels.
	res@tmYLLabelsOn=True
	res@tmYUseLeft=True
	res@tmYLValues    = ispan(2,52,5)
	res@tmYLLabels    = ispan(1960,2010,5)
	res@tmXTOn=False
	res@tmYROn=False

	; res@tmXBValues=(/0,90,180,270,357.5/)
 ;    res@tmXBLabels=(/"0","90~S~o~N~E","180","90~S~o~N~W","0"/)
    res@tmXBValues=(/0,90,180,270/)
    res@tmXBLabels=(/"0","90~S~o~N~E","180","90~S~o~N~W"/)
    
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

    plot@$str1$=gsn_add_polyline(wks, plot, rc_xi_line, ispan(0,54,1), lineres)
    plot@$str2$=gsn_add_polyline(wks, plot, rc_dong_line, ispan(0,54,1), lineres)

    ;;;;draw_attch_plot***************************************************************************
    resxy=True
    resxy@gsnDraw=False
    resxy@gsnFrame=False
    resxy@tmXBMode="Explicit"
    resxy@tmXTOn=False
    resxy@tmYLOn=False
    resxy@tmYROn=False
    resxy@trYMaxF=54
    ; resxy@trXMinF=-2.0
    ; resxy@trXMaxF=3.0
    ; resxy@tmXBValues=ispan(-1,3,1)
    ; resxy@tmXBLabels=ispan(-1,3,1)
    ; resxy@tmXBMinorValues=fspan(-2,3,26 )
    ; resxy@tmXBMinorOn=True

    ; if fils1.ne."era/era40/msf_time_era40_in_global.nc" then 
     if mode.eq.1 then
    resxy@trXMinF=-12.0
    resxy@trXMaxF=12.0
    resxy@tmXBValues=ispan(-12,12,4)
    resxy@tmXBLabels=ispan(-12,12,4)
    resxy@tmXBMinorValues=fspan(-12,12,25)
    end if

    if mode.eq.0 then
    resxy@trXMinF=-6.0
    resxy@trXMaxF=6.0
    resxy@tmXBValues=ispan(-6,6,2)
    resxy@tmXBLabels=ispan(-6,6,2)
    resxy@tmXBMinorValues=fspan(-6,6,13)
    end if
    ; else 
    ; resxy@trXMinF=-20.0
    ; resxy@trXMaxF=20.0
    ; resxy@tmXBValues=ispan(-20,20,10)
    ; resxy@tmXBLabels=ispan(-20,20,10)
    ; resxy@tmXBMinorValues=fspan(-20,20,21)
    ; end if
    resxy@gsnRightString=""+sprintf("%2.2f",rc_width)+"~S~o~N~/dacade"

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

    printVarSummary(wks)

 ; return[/plot,rc_xi_line,rc_dong_line/]
 return(plot)
 end 

 ;***************************************************************
 ;main program
 ;****************************************************************
  begin 
    ; fils1="ncep/20CR/msf_time_20cr_global_6--6.nc"
    fils1="CMIP5/msf_history_rcp45/msf_time_GFDL-CM3_r1i1p1_195001-209912.nc"
    fils2="ncep/ncep1/msf_time_ncep1_global.nc"
    fils3="era/era40/msf_time_era40_in_global.nc"
    fils4="JRA-55/msf_time_jra55_global.nc"

    ; wks=gsn_open_wks("pdf", "/mnt/f/share/era/data/500hpa4个再分析资料计算的msf_hov回归GW")
    wks=gsn_open_wks("pdf", "/mnt/f/share/era/data/700hpa_cmip5资料计算的msf_hov回归GW")

    plot=new(4,graphic)
    mode=0  ;;;mode=o express regress to GW. mode=1 means regress to IPO
    plot(0)=draw_msf_hov(fils1,wks,"20CR",mode)
    plot(1)=draw_msf_hov(fils2,wks,"NCPE1",mode)
    plot(2)=draw_msf_hov(fils3,wks,"ERA",mode)
    plot(3)=draw_msf_hov(fils4,wks,"JRA-55",mode)  
    
    pres=True
    pres@gsnMaximize=True
    pres@gsnPanelLabelBar=True
    pres@pmLabelBarWidthF=0.8
    pres@pmLabelBarHeightF=0.05
    pres@pmLabelBarOrthogonalPosF =0.01
    pres@gsnPanelDebug=True
    gsn_panel(wks, plot, (/2,2/), pres)
  ; draw(plot(0))
  ; frame(wks)
    end