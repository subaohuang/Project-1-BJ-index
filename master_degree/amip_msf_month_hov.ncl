
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


;********************************************************
;draw_msf_hov
;*******************************************************
function draw_msf_hov(wks,LeftString)
begin 
print("************************"+LeftString+"**************************")
    diri="/mnt/f/share/CMIP5/msf_AMIP/"
    fils1=diri+"msf_time_"+LeftString +"_amip_r1i1p1_197901-200812.nc"
    f1=addfile(fils1, "r")
    date=cd_calendar(f1->time, 1)
    timeind=ind(date.ge.197901.and.date.le.200512)

  
    msf_time_month=clmMonLLT(f1->msf_time(:,:,timeind))
    msf_time_month=(/msf_time_month/10^9/)

    printVarSummary(msf_time_month)
    msf_time_month!2="time"

    ; msf_mean_month=level_mean(msf_time_month({850:700},:,:))    ;;msf_mean_month(time,lon)
    levind=ind(msf_time_month&level.eq.500)
    msf_mean_month=msf_time_month(time|:,level|levind,lon|:) 
    printVarSummary(msf_mean_month)


;***********************************************
;draw hov
;***********************************************

    ; cmap=read_colormap_file("MPL_PuOr")
    ; cmap=(/cmap(::-1,0:3)/)
  
	res            = True                         ; plot mods desired
	res@gsnDraw    = False                        ; don't draw
	res@gsnFrame   = False                        ; don't advance frame
    res@gsnLeftString=LeftString
   	res@cnLineLabelsOn =False
	res@cnFillDrawOrder="PreDraw"
	res@cnLinesOn=False
    res@cnFillOn=True
    res@lbLabelBarOn=False     
    res@cnInfoLabelOn=False
    res@cnFillPalette="GMT_polar"
  

	res@cnLevelSelectionMode = "ManualLevels"     ; manually select contours
	res@cnMinLevelValF       = -8.0              ; min level
	res@cnMaxLevelValF       =  8.0               ; max leve
	res@cnLevelSpacingF      =1.0


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
    plot = gsn_csm_contour(wks, msf_mean_month, res)   ; create plot


    res@cnLinesOn=True
    res@cnFillOn        = False
    res@cnLevelSelectionMode="ExplicitLevels"
    res@cnLevels=0.0
    res@cnLineThicknessF=2.0
 
    plot1 = gsn_csm_contour(wks, msf_mean_month, res)
    
    overlay(plot, plot1)

 return(plot)
 end 

 ;***************************************************************
 ;main program
 ;****************************************************************
  begin 
   rea=(/"20cr","era","jra55","MERRA","ncep1","ncep2"/)

   model=(/"ACCESS1-0","BCC-CSM1-1","BNU-ESM","CanAM4","CCSM4","CMCC-CM","CNRM-CM5","CSIRO-Mk3-6-0","EC-EARTH","FGOALS-s2","FGOALS-g2","GFDL-CM3","GFDL-HIRAM-C180","HadGEM2-A","INMCM4","IPSL-CM5A-LR",\
    "IPSL-CM5A-MR","IPSL-CM5B-LR","MIROC5","MPI-ESM-LR","MPI-ESM-MR","MRI-AGCM3-2H","MRI-CGCM3","NorESM1-M"/)



    wks1=gsn_open_wks("pdf", "/mnt/f/share/era/data/amip模式再分析资料沃克环流月气候态")
 

    plot1=new(dimsizes(model),graphic)
  
    do i=0,dimsizes(model)-1
    plot1(i)=draw_msf_hov(wks1,model(i))
    end do
    
    pres=True
    pres@gsnMaximize=True
    pres@gsnPanelLabelBar=True
    pres@pmLabelBarWidthF=0.8
    pres@pmLabelBarHeightF=0.05
    pres@pmLabelBarOrthogonalPosF =0.01
    gsn_panel(wks1, plot1, (/6,4/), pres)
 

  
    end