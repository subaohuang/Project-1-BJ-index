
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
        ; if  (msf_mean_annual_dong(t,i+1).le.0).and.((msf_mean_annual_dong(t,i).ge.0).and.(msf_mean_annual_dong(t,i+3).lt.0))  then 
         if  (msf_mean_annual_dong(t,i+1).le.0).and.((msf_mean_annual_dong(t,i).ge.0))  then 
        locationdong(t)=(msf_mean_annual_dong&lon(i)+msf_mean_annual_dong&lon(i+1))/2
        end if    
        if ismissing(locationdong(t)).ne.True  then 
        break
        end if
     end do
end do
locationxi=(/linmsg(locationxi, -1)/)
locationdong=(/linmsg(locationdong, -1)/)

 location = [/ locationxi, locationdong /]   ; return two variables as a type 'list'
return location
end 


function  search_location_1d(msf_mean_annual)
begin

msf_mean_annual_xi=msf_mean_annual({120:250})
msf_mean_annual_dong=msf_mean_annual({230:300})
locationxi=new(1,typeof(msf_mean_annual&lon))
locationdong=new(1,typeof(msf_mean_annual&lon))
;;;judge locationxi
    do i=0,dimsizes(msf_mean_annual_xi&lon)-9
        if  (msf_mean_annual_xi(i+1).ge.0).and.((msf_mean_annual_xi(i).le.0).and.(msf_mean_annual_xi(i+8).gt.0))  then 
        locationxi=(msf_mean_annual_xi&lon(i)+msf_mean_annual_xi&lon(i+1))/2
        end if    
        if ismissing(locationxi).ne.True  then 
        break
        end if
end do
;;;judge locationdong
    do i=0,dimsizes(msf_mean_annual_dong&lon)-2
        ; if  (msf_mean_annual_dong(t,i+1).le.0).and.((msf_mean_annual_dong(t,i).ge.0).and.(msf_mean_annual_dong(t,i+3).lt.0))  then 
         if  (msf_mean_annual_dong(i+1).le.0).and.((msf_mean_annual_dong(i).ge.0))  then 
        locationdong=(msf_mean_annual_dong&lon(i)+msf_mean_annual_dong&lon(i+1))/2
        end if    
        if ismissing(locationdong).ne.True  then 
        break
        end if
end do
; locationxi=(/linmsg(locationxi, -1)/)
; locationdong=(/linmsg(locationdong, -1)/)
 location = [/ locationxi, locationdong /]   ; return two variables as a type 'list'
return location
end 












function select_win(in_20cr)
begin
in_20cr_win=new(3,typeof(in_20cr))
   in_20cr_win(0:1)=in_20cr(9:10)     ;;;;1,2,12月
   in_20cr_win(2)=in_20cr(11)

return(in_20cr_win)
end  
;********************************************************
;draw_msf_hov
;*******************************************************
function draw_msf_hov(wks,LeftString)
begin 
print("************************"+LeftString+"**************************")
    diri="/mnt/f/share/rea_msf/"
    fils1=diri+"msf_time_"+LeftString+"_global.nc"
    f1=addfile(fils1, "r")
    date=cd_calendar(f1->time, 1)
    timeind=ind(date.ge.197901.and.date.le.201212)

    ; levind=ind(f1->level.eq.500)
    level=f1->level
    
    levind=ind(level.le.600.and.level.ge.400)

    
    msf_time=f1->msf_time(:,:,timeind)
    delete(f1)
    delete(fils1)
    msf_time=(/msf_time/10^9/)
 
    delete(msf_time&time)
    msf_time&time=date(timeind)
    msf_time=(/dtrend_msg_n(ispan(1,dimsizes(msf_time&time),1),msf_time,False,True,2)/)
   
   msf_time_mean=level_mean(msf_time(time|:,level|levind,lon|:))


    msf_annual=month_to_annual(msf_time(time|:,level|:,lon|:), 1)
    msf_annual!0="time"
    msf_annual&time=ispan(1979,2012,1)
    msf_djf=month_to_season(msf_time(time|:,level|:,lon|:), "DJF")
    msf_djf(0:dimsizes(msf_djf&time)-2,:,:)=(/msf_djf(1:,:,:)/)
    msf_djf(dimsizes(msf_djf&time)-1,:,:)=(/msf_time(:,:,{201212})/)
    msf_djf&time=ispan(1979,2012,1)  ;;;因为这里第一个是1958 1月和2月的平均，第二个是195812,195901,195902的平均
    ; location=search_location(msf_djf(:,{500},:))
    msf_djf_mean=level_mean(msf_djf(time|:,level|levind,lon|:))
    location=search_location(msf_time_mean(time|:,lon|:))
   
    locationxi=location[0]
    locationdong=location[1]
    width=locationdong -locationxi

    delete(locationxi&time)
    locationxi&time=date(timeind)
    
    delete(locationdong&time)
    locationdong&time=date(timeind)
    copy_VarCoords(locationxi, width)
     
    delete([/date,timeind/])

    locationxi_annual=month_to_annual(locationxi, 1)
    locationdong_annual=month_to_annual(locationdong, 1)
    width_annual=month_to_annual(width, 1)



    locationxi_djf=month_to_season(locationxi, "DJF")
    locationxi_djf(0:dimsizes(locationxi_djf&time)-2)=(/locationxi_djf(1:)/)
    locationxi_djf(dimsizes(locationxi_djf&time)-1)=(/locationxi({201212})/)
    locationxi_djf&time=ispan(1979,2012,1)  ;;;因为这里第一个是1958 1月和2月的平均，第二个是195812,195901,195902的平均
    
    locationdong_djf=month_to_season(locationdong, "DJF")
    locationdong_djf(0:dimsizes(locationdong_djf&time)-2)=(/locationdong_djf(1:)/)
    locationdong_djf(dimsizes(locationdong_djf&time)-1)=(/locationdong({201212})/)
    locationdong_djf&time=ispan(1979,2012,1)  ;;;因为这里第一个是1958 1月和2月的平均，第二个是195812,195901,195902的平均

    width_djf=month_to_season(width, "DJF")
    width_djf(0:dimsizes(width_djf&time)-2)=(/width_djf(1:)/)
    width_djf(dimsizes(width_djf&time)-1)=(/width({201212})/)
    width_djf&time=ispan(1979,2012,1)  ;;;因为这里第一个是1958 1月和2月的平均，第二个是195812,195901,195902的平均



    printVarSummary(locationxi_djf)
 
    indxi=ind_nearest_coord(locationxi_djf, msf_djf&lon, 0)
    inddong=ind_nearest_coord(locationdong_djf, msf_djf&lon, 0)

msf_inten=new(dimsizes(msf_djf&time),typeof(msf_djf))



do i=0,dimsizes(msf_djf&time)-1
msf_inten(i)=(/dim_avg_n_Wrap(msf_djf_mean(i,{150:210}),0)/)
end do 



    ; rc1=regCoef(locationxi_annual,msf_inten)
    ; line = rc1*locationxi_annual+rc1@yintercept
     rc1=regCoef(width_djf,msf_inten)
     line = rc1*width_djf+rc1@yintercept

   ;;;corr
    r=escorc(width_djf, msf_inten)

  ;;;read sst data
    diri="/mnt/f/share/ersst/"
    fils1=systemfunc("ls "+diri+"sst.mnmean.v4.nc")
    f1=addfile(fils1, "r")
    date=cd_calendar(f1->time, 1)
    timeind=ind(date.le.201312.and.date.ge.197901)

    sst=f1->sst(timeind,{30:-30},:)
    sst!0="time"
    sst!1="lat"
    sst!2="lon"
    delete(date)
    delete(f1)
  
    ssta=(/rmMonAnnCycTLL(sst)/)
    copy_VarCoords(sst, ssta)
    ssta_de=dtrend_msg_n(ispan(1,dimsizes(ssta&time),1),ssta,False,True,0)
    copy_VarCoords(ssta, ssta_de)
    ssta_djf=month_to_season(ssta_de, "DJF")
    ssta_son=month_to_season(ssta_de, "SON")

    ssta_s=(/ssta_djf(1:,:,:)/)*1.0
    ssta_q=(/ssta_son(0:dimsizes(ssta_djf&time)-2,:,:)/)*1.0

    copy_VarCoords(ssta_djf(0:dimsizes(ssta_djf&time)-2,:,:), ssta_s)
    copy_VarCoords(ssta_djf(0:dimsizes(ssta_djf&time)-2,:,:), ssta_q)
  
    emi=dim_avg_n(ssta_q(:,{10:-10},{165:220}),(/1,2/))-0.5*dim_avg_n(ssta_q(:,{20:-10},{125:145}), (/1,2/))-0.5*dim_avg_n(ssta_q(:,{5:-15},{250:290}),(/1,2/))
    emi_std=dim_standardize_n(emi, 1, 0)
 
    nino3=dim_avg_n(ssta_s(:,{5:-5},{210:270}), (/1,2/))
    nino3_std=dim_standardize_n(nino3, 1, 0)
   
    modokiind=ind(emi_std.ge.1)
    caind=ind(nino3_std.ge.1)
    laind=ind(nino3_std.le.(/-1/))
    noind=ind(nino3_std.ge.(/-0.5/).and.nino3_std.le.(/0.5/))
    
width_djf_cano=width_djf(caind)
width_djf_lani=width_djf(laind)
width_djf_modo=width_djf(modokiind)
width_djf_nor=width_djf(noind)


msf_inten_djf_cano=msf_inten(caind)
msf_inten_djf_lani=msf_inten(laind)
msf_inten_djf_modo=msf_inten(modokiind)
msf_inten_djf_nor=msf_inten(noind)
printVarSummary(msf_djf)
msf_cano=msf_djf({1997},:,:)

location_djf=search_location(msf_djf(:,{500},:))
location_djf_xi=location_djf[0]
location_djf_dong=location_djf[1]
delete(location_djf_xi&time)
delete(location_djf_dong&time)
copy_VarCoords(msf_djf(:,1,1), location_djf_xi)
copy_VarCoords(msf_djf(:,1,1), location_djf_dong)


location1d=search_location_1d(msf_djf({1997},{500},:))
print("locationxi1d="+location1d[0])
print("locatiodong1d="+location1d[1])


print("lo_xi="+location_djf_xi({1997}))
print("lo_xi="+location_djf_dong({1997}))
print("width_djf="+width_djf({1997}))




;;;画散点图


  
    str1 = unique_string("dum1")  ; "unique_string" will return a unique                                     ; string every time it is called fro                                 ;  within a single NCL session.
    str2 = unique_string("dum2")
    str3 = unique_string("dum3")
    str4 = unique_string("dum4")

       res = True
    res@gsnDraw=False
    res@gsnFrame=False
    res@gsnLeftString=LeftString
    res@gsnRightString="10~S~9~N~Kg s~S~-1"
    res@gsnMaximize=False
    ; res@lbLabelStrings=ispan(-6,6,1)
    res@cnFillPalette  = "BlueDarkRed18"
    res@cnLinesOn=True
    res@cnLevelSelectionMode ="ManualLevels"
    res@cnMinLevelValF = -6
    res@cnMaxLevelValF = 6
    res@cnLevelSpacingF = 1.0
    res@cnFillOn = True 
    ; res@gsnYRefLine=500
    ; res@gsnYRefLineColor="black"

    res@gsnContourNegLineDashPattern=1
    res@gsnContourZeroLineThicknessF=2.5

    res@cnLineLabelsOn=False
    res@cnMonoFillScale=True
    res@cnInfoLabelOn=False

    ; res@cnExplicitLineLabelsOn=True
    ; res@cnLineLabelStrings=""+ispan(-2,6,1)+""
    ; res@cnLabelScalingMode="ScaleFactor"
    ; res@cnLabelScaleValueF=10^9

    res@lbLabelBarOn=False
    ; res@cnLineColors=(/0,0,0,0,0,0,18,18,18,18,18,18/)
    res@tmYRMode="Automatic"
    res@tmXTOn=False
    res@tmYROn=False
    res@tmXBMode="Explicit"
    res@tmXBValues=(/0,90,180,270,357.5/)
    res@tmXBLabels=(/"0","90~S~o~N~E","180","90~S~o~N~W","0"/)
    res@tmXBMinorOn=True
    res@tmXBMinorValues=ispan(0,360,10)
     
    res@tiYAxisString=""

    res@vpHeightF=0.3
    res@vpWidthF=0.7

   plot=gsn_csm_pres_hgt(wks,msf_cano,res)

return (plot)
end
;;;主程序

begin
rea=(/"20cr","era","jra55","MERRA","ncep1","ncep2"/)

plot=new(6,graphic)

wks=gsn_open_wks("pdf", "/mnt/f/share/era/data/6个再分析资料djf 1997年质量流函数")

do i=0,5
plot(i)=draw_msf_hov(wks,rea(i))
end do


  xleg = (/0.2,0.4,0.6,0.8/)   ; Location of

  xtxt = (/0.25,0.46,0.65,0.85/)   ; legend markers

  yleg = (/0.1,0.1,0.1,0.1/)   ; and text

  ytxt = (/0.1,0.1,0.1,0.1/)   ; strings.
 
  colors=(/"red","blue","pink","black"/)
  labels=(/" EP Enino","EP Lanina","CP Enino","Neutral"/)

  txres               = True
  txres@txFontHeightF = 0.015
  
   gsres               = True
  gsres@gsMarkerIndex = 16          ; Use filled dots for markers.

    do i = 0, 3
      gsres@gsMarkerColor      = colors(i)
      ; gsres@gsMarkerThicknessF = 0.7
      if i.eq.3 then
         gsres@gsMarkerSizeF=0.005
     else
     gsres@gsMarkerSizeF=10.0
 end if
      gsn_polymarker_ndc(wks, xleg(i),yleg(i),gsres)
      gsn_text_ndc      (wks,labels(i),xtxt(i),ytxt(i),txres)
    end do


pres=True
pres@gsnMaximize=True
pres@gsnPanelDebug=True
gsn_panel(wks, plot, (/3,2/), pres)
end 