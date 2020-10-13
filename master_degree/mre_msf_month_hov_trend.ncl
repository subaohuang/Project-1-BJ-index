
function level_mean(data)
local data,msf_mean
begin
plev=data&level
psfc=conform_dims((/dimsizes(data&time),dimsizes(data&lon)/), 1013.0, 0)
ptop=min(plev)
dp   = dpres_plevel_Wrap(plev, psfc, ptop, 0)  ; dp(30)
msfdp=data(level|:,time|:,lon|:)*dp
copy_VarCoords(data(level|:,time|:,lon|:), msfdp)

msf_mean=dim_sum_n_Wrap(msfdp, 0)/dim_sum_n_Wrap(dp, 0)
copy_VarCoords(data(time|:,lon|:,level|1), msf_mean)
return(msf_mean)
end 

function  rc_test(ssta_rc)
local ssta_rc
begin
  tval = onedtond(ssta_rc@tval , dimsizes(ssta_rc))
   df   = onedtond(ssta_rc@nptxy, dimsizes(ssta_rc)) - 2
   b = tval    ; b must be same size as tval (and df)
   b = 0.5
   prob = betainc(df/(df+tval^2),df/2.0,b)    ; prob(nlat,nlon)
   copy_VarCoords(ssta_rc, prob)
return prob
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

    f1=addfile(diri+fils1, "r")
    date=cd_calendar(f1->time, 1)
    timeind=ind(date.ge.195801.and.date.le.201212)

    msf_time=f1->msf_time({500},:,timeind)
    msf_time=(/msf_time/10^9/)
    msf_mean=msf_time(time|:,lon|:)
    ; msf_mean=level_mean(msf_time)      ;;msf_mean_month(time,lon)  

    msf_month_trend=new((/12,dimsizes(msf_mean&lon)/),double)
    prob=msf_month_trend

   
    do i=0,11
    msf_month_trend(i,:)=regCoef_n(ispan(0,54,1), msf_mean(i:659:12,:), 0, 0)
    prob(i,:)=rc_test(msf_month_trend(i,:))
    end do


    msf_time_month=clmMonLLT(f1->msf_time)
    msf_time_month=(/msf_time_month/10^9/)

    printVarSummary(msf_time_month)
    msf_time_month!2="time"

    msf_mean_month=level_mean(msf_time_month({600:400},:,:))    ;;msf_mean_month(time,lon)
    ; msf_mean_month_500=msf_time_month({500},:,:)
    ; msf_mean_month=msf_mean_month_500(time|:,lon|:)

    copy_VarCoords(msf_mean_month, msf_month_trend)
    copy_VarCoords(msf_mean_month, prob)


    ;;;回归GW和IPO信号
    ;;;
    ipindex1=asciiread("/mnt/f/share/era/data/GHCN_HAD_p_ud_timef1_1958_2012.ascii", 55, "float")
    ipindex2=asciiread("/mnt/f/share/era/data/GHCN_HAD_p_ud_timef2_1958_2012.ascii", 55, "float")
    ;;;    
           if mode.ne.0 then 
             ; msf_re1=msf_mean
             ; msf_re1_trend=msf_month_trend
             ; prob_re1=prob
             ; msf_re1_month=msf_mean_month
             msf_re1=new(dimsizes(msf_mean),typeof(msf_mean))
             msf_re1_trend=new(dimsizes(msf_month_trend),typeof(msf_month_trend))
             prob_re1=new(dimsizes(prob),typeof(prob))
             msf_re1_month=new(dimsizes(msf_mean_month),typeof(msf_mean_month))
             ;;;GW
             if mode.eq.1 then 
            do i=0,11
            msf_re1(i:659:12,:)=fit(msf_mean(i:659:12,:),ipindex1)
            end do
            do i=0,11  
            msf_re1_month(i,:)=dim_avg_n_Wrap(msf_re1(i:659:12,:), 0)
            end do
            do i=0,11
            msf_re1_trend(i,:)=regCoef_n(ispan(0,54,1), msf_re1(i:659:12,:), 0, 0)
            prob_re1(i,:)=rc_test(msf_re1_trend(i,:))
            end do
            msf_mean_month=(/msf_re1_month/)
            msf_month_trend=(/msf_re1_trend/)
            end if
            ;;;IPO
            if mode.eq.2 then 
            do i=0,11
            msf_re1(i:659:12,:)=fit(msf_mean(i:659:12,:),ipindex2)
            end do
            do i=0,11  
            msf_re1_month(i,:)=dim_avg_n_Wrap(msf_re1(i:659:12,:), 0)
            end do
            do i=0,11
            msf_re1_trend(i,:)=regCoef_n(ispan(0,54,1), msf_re1(i:659:12,:), 0, 0)
            prob_re1(i,:)=rc_test(msf_re1_trend(i,:))
            end do
            msf_mean_month=(/msf_re1_month/)
            msf_month_trend=(/msf_re1_trend/)
            end if
        end if
    ;;;
 
 ;;;

            

;***********************************************
;draw hov
;***********************************************

    ; cmap=read_colormap_file("MPL_PuOr")
    ; cmap=(/cmap(::-1,0:3)/)
  
	res            = True                         ; plot mods desired
	res@gsnDraw    = False                        ; don't draw
	res@gsnFrame   = False                        ; don't advance frame
    res@gsnLeftString=LeftString
    if mode.eq.1 then 
    res@gsnRightString="Regress GW"
    end if 
      if mode.eq.2 then 
    res@gsnRightString="Regress IPO"
    end if 
	res@cnLineLabelsOn =False
	res@cnFillDrawOrder="PreDraw"
	res@cnLinesOn=False
    res@cnFillOn=True
    res@lbLabelBarOn=False     
    res@cnInfoLabelOn=False
    res@cnFillPalette="GMT_polar"
    if mode.eq.2 then 
          
    res@cnLevelSelectionMode = "ManualLevels"     ; manually select contours
    res@cnMinLevelValF       = -0.018               ; min level
    res@cnMaxLevelValF       =  0.018               ; max leve
    res@cnLevelSpacingF      =0.002
    else    

	res@cnLevelSelectionMode = "ManualLevels"     ; manually select contours
	res@cnMinLevelValF       = -0.05              ; min level
	res@cnMaxLevelValF       =  0.05               ; max leve
	res@cnLevelSpacingF      =0.01
	end if



   
	res@tmYLMode      = "Explicit"                ; Define own tick mark labels.
	res@tmYLLabelsOn=True
	res@tmYUseLeft=True
	res@tmYLValues    = ispan(0,11,1)
	res@tmYLLabels    = (/"J","F","M","A","M","J","J","A","S","O","N","D"/)
	res@tmXTOn=False
	res@tmYROn=False

	res@tmXBValues=(/0,90,180,270,357.5/)
    res@tmXBLabels=(/"0","90~S~o~N~E","180","90~S~o~N~W","0"/)
    res@tmXBMinorOn=True
    res@tmXBMinorValues=ispan(0,360,10)

	res@vpHeightF=0.2
	res@vpWidthF=0.4
    plot = gsn_csm_contour(wks, msf_month_trend, res)   ; create plot



    ;;;;打点上去
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
   
    iplot=gsn_csm_contour(wks,prob,res_pro)

    opt                      = True
    opt@gsnShadeFillType     = "pattern"
    opt@gsnShadeMid     = 17         ; stipple
    ; opt@gsnShadeHigh   = "green"
    opt@gsnShadeFillScaleF   = 0.2         ; add extra density
    opt@gsnShadeFillDotSizeF = 0.0      ; make dots larger
     
    iplot= gsn_contour_shade(iplot, 0,0.05, opt)

     overlay(plot, iplot)
    ; res@cnLevelSelectionMode = "ManualLevels"     ; manually select contours
    ; res@cnMinLevelValF       = -4.0               ; min level
    ; res@cnMaxLevelValF       =  4.0               ; max leve
    ; res@cnLevelSpacingF      =1.0

    res@cnLinesOn=True
    res@cnFillOn        = False
    res@cnLineLabelsOn=True
    res@gsnContourZeroLineThicknessF=4.0
    res@gsnContourNegLineDashPattern=1
    res@gsnContourPosLineDashPattern=0
     
    res@cnLevelSelectionMode = "ManualLevels"     ; manually select contours
    res@cnMinLevelValF       = -6.0               ; min level
    res@cnMaxLevelValF       =  6.0               ; max leve
    res@cnLevelSpacingF      =1.0


    plot1 = gsn_csm_contour(wks, msf_mean_month, res)
    
    overlay(plot, plot1)

 return(plot)
 end 

 ;***************************************************************
 ;main program
 ;****************************************************************
  begin 
    fils1="ncep/20CR/msf_time_20cr_global_6--6.nc"
    fils2="ncep/ncep1/msf_time_ncep1_global.nc"
    fils3="era/era40/msf_time_era40_in_global.nc"
    fils4="JRA-55/msf_time_jra55_global.nc"

    wks=gsn_open_wks("pdf", "/mnt/f/share/era/data/500hpa_4个再分析资料计算的msf_month_hov_IPO")

    plot=new(4,graphic)
    mode=0
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
     gsn_panel(wks, plot, (/2,2/), pres)

  
    end