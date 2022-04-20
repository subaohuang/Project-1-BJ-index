
;;***********************************************************************
;;draw 暖池SST的时间经度图
;;;**********************************************************************
;;function
;;; 
function closest_val_AnyOrder(xVal[1]:numeric, x:numeric)
local xAbsDif, xMinVal, iClose
begin
  xAbsDif = abs(xVal-x)
  xMinVal = min(xAbsDif)
  iClose  = (ind(xAbsDif.eq.xMinVal))
  return(iClose(0))   ; original doc says "first occurrence"
end

undef("search_location")
function  search_location(msf_mean_annual)
local msf_mean_annual_xi 
begin
locationxi=new(dimsizes(msf_mean_annual&time),typeof(msf_mean_annual&lon))
locationxi!0="time"
locationxi&time=ispan(0,dimsizes(msf_mean_annual&time)-1,1)
msf_mean_annual_xi=msf_mean_annual(:,{170:280})
; co=new(dimsizes(msf_mean_annual&time),integer)
;;;judge locationxi
do t=0,dimsizes(msf_mean_annual&time)-1
     co=closest_val_AnyOrder(28.0,msf_mean_annual_xi(t,:))  
     locationxi(t)=msf_mean_annual_xi&lon(co)
    ; delete(co)
end do
return locationxi
end 

;;set
  begin
    diri="/mnt/f/share/ersst/"
	fils1=systemfunc("ls "+diri+"sst.mnmean.v4.nc")
	f1=addfile(fils1, "r")
	date=cd_calendar(f1->time, 1)
	timeind=ind(date.le.201212.and.date.ge.197901)
	sst=f1->sst(timeind,{10:-10},{70:280})
	sst!0="time"
	sst!1="lat"
	sst!2="lon"
    sst&time=ispan(1,dimsizes(timeind),1)
	sst_mean=dim_avg_n_Wrap(sst, 1)
    sst_lo=search_location(sst_mean)

 

    lo_trend=regline(sst_mean&time, sst_lo)
    y = lo_trend*sst_mean&time + lo_trend@yintercept   ; NCL array notation 

    


wks=gsn_open_wks("pdf", "/mnt/f/share/era/data/暖池SST的时间经度图")

;;;;;draw
	res            = True                         ; plot mods desired
	res@gsnDraw    = False                        ; don't draw
	res@gsnFrame   = False                        ; don't advance frame
    res@cnLineLabelsOn =False
    res@cnFillDrawOrder="PreDraw"
    res@cnLinesOn=False

    res@cnFillPalette  = "BlueDarkRed18"
 
    res@cnFillOn=True
    res@lbLabelBarOn=False     
    res@cnInfoLabelOn=False
    res@cnLevelSelectionMode = "ManualLevels"     ; manually select contours
    res@cnMinLevelValF       = 20.0               ; min level
    res@cnMaxLevelValF       = 32.0               ; max leve
    res@cnLevelSpacingF      =1.0
    res@lbOrientation         ="Vertical"

    res@tmYLMode      = "Explicit"                ; Define own tick mark labels.
    res@tmYLLabelsOn=True
    res@tmYUseLeft=True
    res@tmYLValues    = ispan(1,dimsizes(timeind),12*5)
    res@tmYLLabels=ispan(1979, 2012,5)

    res@tmXTOn=False
    res@tmYROn=False
      ; res@tmXBMode       ="Explicit"
    ; res@tmXBValues=ispan(100, 280, 60)
    ; res@tmXBValues=(/30,90,150,210,270,330,387.5/)
    ; res@tmXBLabels=(/"30~S~o~N~E","90~S~o~N~E","150~S~o~N~E","150~S~o~N~W"\
    ;                 ,"90~S~o~N~W","30~S~o~N~W","30~S~o~N~E"/)
   
    ; res@tmXBMinorOn=True
    ; res@tmXBMinorValues=ispan(30,390,10)
    ; res@tmYLMinorOn=True
    ; res@tmYLMinorValues=ispan(1,34,1)
    res@tmXBLabelFontHeightF=0.02
    res@tmXBLabelFontThicknessF=5.0
    res@tmYLLabelFontThicknessF=5.0
    res@tmYLLabelFontHeightF=0.02
    res@tmYLMajorThicknessF=4.0
    res@tmXBMajorThicknessF=4.0
    res@tmBorderThicknessF=4.0
    res@vpHeightF=0.6
    res@vpWidthF=0.3
    plot = gsn_csm_contour(wks, sst_mean(:,:), res)   ; create plot\

    res1=res
    res1@tmXBOn=False
    res1@cnLevelSelectionMode = "ExplicitLevels"     ; manually select contours
    res1@cnLinesOn=True
    res1@cnLevels=(/28.75/)
    res1@cnFillOn=False
    res1@cnLineThicknesses =30.0
    
    plot1 = gsn_csm_contour(wks, sst_mean(:,{70:280}), res1)   ; create plot\    
    
    overlay(plot, plot1)

    res_line=True
    res_line@gsLineColor="red"
    res_line@gsLineThicknessF=2.0
    dum=gsn_add_polyline(wks, plot,sst_lo,sst_mean&time , res_line)
    res_line@gsLineColor="black"
    dum1=gsn_add_polyline(wks,plot,y,sst_mean&time , res_line)

     
    draw(plot)
    frame(wks)

    end 

    

