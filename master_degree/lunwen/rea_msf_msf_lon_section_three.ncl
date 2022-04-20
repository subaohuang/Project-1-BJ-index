;************************************************
; 画20cr,ncep,era40及jRA-55的msf的变化趋势。
;************************************************
;;;
undef("level_mean")
function level_mean(msf_time)   ;;;;select 850hpa-200hpa 
begin
; msf_time=(/msf_time/10^9/)
plev=msf_time&level
psfc=conform_dims((/dimsizes(msf_time&time),dimsizes(msf_time&lon)/), 1013.0, 0)
ptop=min(plev)
dp   = dpres_plevel_Wrap(plev, psfc, ptop, 0)  ; dp(30)
; dp=(/50/)
msfdp=msf_time(level|:,time|:,lon|:)*dp
copy_VarCoords(msf_time(level|:,time|:,lon|:), msfdp)

msf_mean=dim_sum_n_Wrap(msfdp, 0)/dim_sum_n_Wrap(dp, 0)
copy_VarCoords(msf_time(time|:,lon|:,level|1), msf_mean)
return(msf_mean)
end
undef("search_location_2d")
function  search_location_2d(msf_mean_annual)
begin
msf_mean_annual_xi=msf_mean_annual({120:250})
msf_mean_annual_dong=msf_mean_annual({230:300})
locationxi=msf_mean_annual@_FillValue
locationdong=locationxi
;;;judge locationxi
    do i=0,dimsizes(msf_mean_annual_xi&lon)-9
        if  (msf_mean_annual_xi(i+1).ge.0).and.((msf_mean_annual_xi(i).le.0).and.(msf_mean_annual_xi(i+8).gt.0))  then 
         a=(msf_mean_annual_xi(i+1)-msf_mean_annual_xi(i))/(msf_mean_annual_xi&lon(i+1)-msf_mean_annual_xi&lon(i))
         locationxi=tofloat(msf_mean_annual_xi&lon(i)- msf_mean_annual_xi(i)/a)
        ; locationxi=(msf_mean_annual_xi&lon(i)+msf_mean_annual_xi&lon(i+1))/2
        end if    
        if  (msf_mean_annual_xi(i).eq.0)  then 
        locationxi=msf_mean_annual_dong&xi(i)
        end if  
        if ismissing(locationxi.ne.True)  then 
        break
        end if
     end do
;;;judge locationdong
    do i=0,dimsizes(msf_mean_annual_dong&lon)-2
        ; if  (msf_mean_annual_dong(t,i+1).le.0).and.((msf_mean_annual_dong(t,i).ge.0).and.(msf_mean_annual_dong(t,i+3).lt.0))  then 
         if  (msf_mean_annual_dong(i+1).le.0).and.((msf_mean_annual_dong(i).ge.0))  then 
         a=(msf_mean_annual_dong(i+1)-msf_mean_annual_dong(i))/(msf_mean_annual_dong&lon(i+1)-msf_mean_annual_dong&lon(i))
         locationdong=tofloat(msf_mean_annual_dong&lon(i)- msf_mean_annual_dong(i)/a)
        ; locationdong=(msf_mean_annual_dong&lon(i)+msf_mean_annual_dong&lon(i+1))/2
        end if    
         if  (msf_mean_annual_dong(i).eq.0)  then 
        locationdong=msf_mean_annual_dong&lon(i)
        end if    
        if ismissing(locationdong.ne.True)  then 
        break
      end if
        if ismissing(locationdong).eq.True  then 
        co_in=closest_val(0,msf_mean_annual_dong)
        locationdong=msf_mean_annual_dong&lon(co_in)
        delete(co_in)
        end if
     end do
  location = [/ locationxi, locationdong /]   ; return two variables as a type 'list'
return location
end 

undef("polygon_xy")
function polygon_xy(msf_rea1)
begin

 nlon=dimsizes(msf_rea1&lon)
 xp    = new( (/2*nlon/), float )
 yp    = new( (/2*nlon/), float )

 do k=0,nlon-1
   yp(k) = max(msf_rea1(:,k))
   xp(k) = msf_rea1&lon(k)
   xp(2*nlon-1-k) = msf_rea1&lon(k)
   yp(2*nlon-1-k) = min(msf_rea1(:,k))
 end do
return [/xp,yp/]
end
;;;;;;;
undef("search_max_min")
function search_max_min(w_rea1_mean)
begin
w_rea1_xi=w_rea1_mean&lon(maxind(w_rea1_mean))
w_rea1_dong=w_rea1_mean&lon(minind(w_rea1_mean))

return [/w_rea1_xi,w_rea1_dong/]
end
;;;;
undef("draw_msf_trend")
function draw_msf_trend(LeftString,wks)
begin
    diri="/mnt/f/share/rea_msf/"
    f1=addfile(diri+"msf_time_"+LeftString+"_global.nc", "r")
    date=cd_calendar(f1->time, 1)
    timeind=ind(date.ge.197901.and.date.le.201212)
    level=f1->level
    levind=ind(f1->level.eq.500)
    levind2=ind(level.ge.400.and.level.le.600)
    

    msf_time_20cr=f1->msf_time(:,:,timeind)
    msf_time_20cr=(/msf_time_20cr/10^9/)
    delete(msf_time_20cr&time)
    msf_time_20cr&time=date(timeind)

    w_time=f1->w_time(:,:,timeind)
    delete(w_time&time)
    w_time&time=date(timeind)

    delete(f1) 
    delete(diri)
    delete(date)
    delete(timeind)
;;;去趋势
    ; msf_time_20cr_de=dtrend_msg_n(ispan(1,dimsizes(msf_time_20cr&time),1),msf_time_20cr,False,True,2)
    msf_time_20cr_de=msf_time_20cr
    copy_VarCoords(msf_time_20cr, msf_time_20cr_de)

    ;;;先得到秋季和冬季eninio年的沃克环流情况SONDJF
    msf_son=month_to_season(msf_time_20cr_de(time|:,level|:,lon|:), "SON")
    msf_djf=month_to_season(msf_time_20cr_de(time|:,level|:,lon|:), "DJF")

    msf_son&time=ispan(1979,2012,1)
    msf_djf&time=ispan(1979,2012,1)  ;;;因为这里第一个是1958 1月和2月的平均，第二个是195812,195901,195902的平均
    ;;注意点
    msf_djf(0:dimsizes(msf_djf&time)-2,:,:)=(/msf_djf(1:,:,:)/)
    msf_djf(dimsizes(msf_djf&time)-1,:,:)=(/msf_time_20cr_de(:,:,{201212})/)
    ;;;
    
    ; msf_djf=(/month_to_annual(msf_time_20cr_de(time|:,level|:,lon|:), 1)/)
    w_djf=month_to_season(w_time(time|:,level|:,lon|:), "DJF")
    w_djf&time=ispan(1979,2012,1)  ;;;因为这里第一个是1958 1月和2月的平均，第二个是195812,195901,195902的平均

    ;;注意点
    w_djf(0:dimsizes(w_djf&time)-2,:,:)=(/w_djf(1:,:,:)/)
    w_djf(dimsizes(w_djf&time)-1,:,:)=(/w_time(:,:,{201212})/)
    printVarSummary(msf_djf)
    ;;;


    msf_djf_500=level_mean(msf_djf(time|:,level|levind2,lon|:))
   

    w_djf_500=w_djf(time|:,level|levind,lon|:)

;;;;水平插值
   xo=fspan(0,357.5,144)
    printVarSummary(msf_djf)
    printVarSummary(msf_djf&lon)
    print("*****************开始水平插值************************")  
    msf_djf_c=tofloat(linint1_n_Wrap(msf_djf_500&lon, msf_djf_500, \
            True, xo, 0,1))
    msf_djf_c!1="lon"

    w_djf_c=tofloat(linint1_n_Wrap(w_djf_500&lon, w_djf_500, \
            True, xo, 0,1))
    w_djf_c!1="lon"

;;;画图
if LeftString.ne."mre" then
return [/msf_djf_c,w_djf_c/]
else 
return [/msf_djf_c,w_djf_c,msf_djf/]
end if
end

;***********************************************************
;主程序
;***********************************************************         
begin


wks=gsn_open_wks("pdf", "/mnt/f/share/era/data/rea_msf与lon的section图_3子图")


;;;画图
lo1=draw_msf_trend("20cr",wks)
lo2=draw_msf_trend("era",wks)
lo3=draw_msf_trend("jra55",wks)
lo4=draw_msf_trend("MERRA",wks)
lo5=draw_msf_trend("ncep1",wks)
lo6=draw_msf_trend("ncep2",wks)

lo_mre=draw_msf_trend("mre",wks)

msf_djf_mre=lo_mre[2]
printVarSummary(msf_djf_mre)

msf_djf_mre1=dim_avg_n_Wrap(msf_djf_mre({1979:1994},:,:),0)
msf_djf_mre2=dim_avg_n_Wrap(msf_djf_mre({1995:},:,:),0)



; lo5=draw_msf_trend(fils5,wks,"mre",pa)


 lo7 = NewList("lifo")
;;;MRE还是取平均比较稳
do i=0,1
 ListPush(lo7, (lo1[i]+lo2[i]+lo3[i]+lo4[i]+lo5[i]+lo6[i])/6.0)
end do

msf_rea=new((/6,34,144/),float)
w_rea=msf_rea
msf_rea(0,:,:)=lo1[0]
msf_rea(1,:,:)=lo2[0]
msf_rea(2,:,:)=lo3[0]
msf_rea(3,:,:)=lo4[0]
msf_rea(4,:,:)=lo5[0]
msf_rea(5,:,:)=lo6[0]

w_rea(0,:,:)=lo1[1]
w_rea(1,:,:)=lo2[1]
w_rea(2,:,:)=lo3[1]
w_rea(3,:,:)=lo4[1]
w_rea(4,:,:)=lo5[1]
w_rea(5,:,:)=lo6[1]

msf_rea!0="rea"
msf_rea!1="time"
msf_rea!2="lon"
msf_rea&rea=ispan(1,6,1)
msf_rea&time=ispan(1979,2012,1)
msf_rea&lon=fspan(0,357.5,144)

msf_rea1=msf_rea(:,{1979:1994},{120:280})
msf_rea2= msf_rea(:,{1995:},{120:280})



msf_rea1_mean=dim_avg_n_Wrap(msf_rea1, (/0,1/))
msf_rea2_mean=dim_avg_n_Wrap(msf_rea2, (/0,1/))


msf_rea1_lo=search_location_2d(msf_rea1_mean)
msf_rea2_lo=search_location_2d(msf_rea2_mean)

msf_rea1_xi=conform(ispan(-6,6,1), msf_rea1_lo[0],0)
msf_rea1_dong=conform(ispan(-6,6,1), msf_rea1_lo[1],0)

msf_rea2_xi=conform(ispan(-6,6,1), msf_rea2_lo[0],0)
msf_rea2_dong=conform(ispan(-6,6,1), msf_rea2_lo[1],0)

print("msf_rea1_xi="+msf_rea1_lo[0])
print("msf_rea1_dong="+msf_rea1_lo[1])

print("msf_rea2_xi="+msf_rea2_lo[0])
print("msf_rea2_dong="+msf_rea2_lo[1])


w_rea!0="rea"
w_rea!1="time"
w_rea!2="lon"
w_rea&rea=ispan(1,6,1)
w_rea&time=ispan(1979,2012,1)
w_rea&lon=fspan(0,357.5,144)

w_rea1=w_rea(:,{1979:1994},{100:280})
w_rea2= w_rea(:,{1995:},{100:280})

w_rea1_mean=dim_avg_n_Wrap(w_rea1, (/0,1/))
w_rea2_mean=dim_avg_n_Wrap(w_rea2, (/0,1/))


w_rea1_lo=search_max_min(w_rea1_mean)
w_rea2_lo=search_max_min(w_rea2_mean)



w_rea1_xi=conform(ispan(-2,4,1), w_rea1_lo[0],0)
w_rea1_dong=conform(ispan(-2,4,1), w_rea1_lo[1],0)

w_rea2_xi=conform(ispan(-2,4,1), w_rea2_lo[0],0)
w_rea2_dong=conform(ispan(-2,4,1), w_rea2_lo[1],0)

print("w_rea1_xi="+w_rea1_lo[0])
print("w_rea1_dong="+w_rea1_lo[1])

print("w_rea2_xi="+w_rea2_lo[0])
print("w_rea2_dong="+w_rea2_lo[1])


printVarSummary(w_rea1_xi)
;;;polygon_xy

  po_xy1=polygon_xy(dim_avg_n_Wrap(msf_rea1,1))
  po_xy2=polygon_xy(dim_avg_n_Wrap(msf_rea2,1))

  po_xy3=polygon_xy(dim_avg_n_Wrap(w_rea1,1))
  po_xy4=polygon_xy(dim_avg_n_Wrap(w_rea2,1))

;;;

res=True
res@gsnDraw=False
res@gsnFrame=False
res@gsnYRefLineDashPattern=1
res@tmXTOn=False
res@tmYROn=False
res@tmXBMode="Explicit"
res@trXMaxF=280
res@trXMinF=120

res@tmXBValues=ispan(120,270,30)
res@tmXBLabels=(/"120~S~o~N~E","150~S~o~N~E","180","150~S~o~N~W","120~S~o~N~W","90~S~o~N~W"/)
res@tmXBMinorOn=True
res@tmXBMinorValues=ispan(120,280,10)


res@xyDashPatterns=(/0,0/)
res@xyLineThicknesses=(/3,3/)
res@xyLineColors=(/"blue","deeppink"/)

res@vpHeightF=0.16
res@vpWidthF=0.55
res@vpXF=0.27
res@vpYF=0.65


plot=new(3,graphic)
plot(1)=gsn_csm_xy(wks, msf_rea1_mean&lon,(/msf_rea1_mean,msf_rea2_mean/), res)


gsres                   = True                        ; poly res
gsres@gsFillColor       = "blue"                 ; color chosen
gsres@gsFillOpacityF    =0.1
gsres@tfPolyDrawOrder="PreDraw"
dummy = gsn_add_polygon (wks,plot(1),po_xy1[0],po_xy1[1],gsres)
gsres@gsFillOpacityF    =0.05
gsres@gsFillColor       = "deeppink"   
dummy1 = gsn_add_polygon (wks,plot(1),po_xy2[0],po_xy2[1],gsres)
gsres@gsFillColor       = "blue"      
gsres@gsFillOpacityF    =0.1                      ; color chosen

polyline=True
polyline@gsLineThicknessF=2.0
polyline@tfPolyDrawOrder="PreDraw"
dummy4 = gsn_add_polyline(wks,plot(1),msf_rea1_mean&lon,conform(msf_rea1_mean&lon, 0, 0),polyline)

;;;添加位置
polyline@gsLineThicknessF=1.0
polyline@gsLineDashPattern=1
polyline@gsLineColor="blue"

res1=True
res1@gsnDraw=False
res1@gsnFrame=False
res1@xyDashPatterns=(/1,1/)
res1@xyLineColors=(/"blue","blue"/)


dummy6 = gsn_csm_xy(wks,(/msf_rea1_xi,msf_rea1_dong/), ispan(-6,6,1),res1)
res1@xyLineColors=(/"deeppink","deeppink"/)
dummy7 = gsn_csm_xy(wks,(/msf_rea2_xi,msf_rea2_dong/), ispan(-6,6,1),res1)
res1@xyLineColors=(/"blue","blue"/)

; overlay(plot(1), dummy6)
; overlay(plot(1), dummy7)

;;;draw plot(0)

    res2 = True
    res2@gsnDraw=False
    res2@gsnFrame=False
    res2@gsnMaximize=False
    res2@cnLinesOn=True
    res2@cnInfoLabelOn=False
    res2@lbLabelBarOn=False
    res2@cnLineLabelsOn=False
    res2@tmYRMode="Automatic"
    res2@tmXTOn=False
    res2@tmXBOn=False
    res2@tmYROn=False
    res2@tiYAxisString=""
    res2@cnLevelSelectionMode ="ExplicitLevels"
    res2@cnLevels=(/0,1,2,3,4,5,6/)
    res2@cnFillOn=False
    res2@cnLineColor="blue"
    res2@cnMonoLineColor=True
    res2@cnLineThicknessF =2.0
    ; res2@gsnContourZeroLineThicknessF=4.0

    res2@cnLineDashPattern=2.0
    res2@vpHeightF=0.31
    res2@vpWidthF=0.55
    res2@vpXF=0.27
    res2@vpYF=0.99

    plot(0)=gsn_csm_pres_hgt(wks, msf_djf_mre1(:,{120:280}), res2)
    res2@cnLineColor="red"
    plot1=gsn_csm_pres_hgt(wks, msf_djf_mre2(:,{120:280}), res2)
    overlay(plot(0), plot1)

    res2@cnLevelSelectionMode ="ExplicitLevels"
    delete(res2@cnLevels)
    res2@cnLevels=(/0/)
    res2@cnLineThicknessF =4.0
    res2@cnLineDashPattern=0
    res2@cnLineColor="blue"
    plot3=gsn_csm_pres_hgt(wks, msf_djf_mre1(:,{120:280}), res2)
    res2@cnLineColor="red"
    plot4=gsn_csm_pres_hgt(wks, msf_djf_mre2(:,{120:280}), res2)
    overlay(plot(0), plot3)
    overlay(plot(0), plot4)
 ;***********************************************
; legend resources
;***********************************************
  lgres                    = True
  lgres@lgLineColors     = (/"blue","red"/)
  lgres@lgLineThicknessF   = res2@cnLineThicknessF  ; legend line thickness = contour line thickness
  lgres@lgMonoDashIndex=True
  lgres@lgDashIndex=0
  lgres@lgLabelFontHeightF = .08            ; set the legend label font thickness
  lgres@vpWidthF           = 0.15           ; width of legend (NDC)
  lgres@vpHeightF          = 0.1            ; height of legend (NDC)
  lgres@lgPerimColor       = "orange"       ; draw the box perimeter in orange
  lgres@lgPerimOn=False
  lgres@lgPerimThicknessF  = 5.0            ; thicken the box perimeter
  lbid = gsn_create_legend(wks,2,(/" 1979-1994","1995-2012"/),lgres)         ; create legend

  amres = True
  amres@amParallelPosF   = -0.35         ; move legend to the right
  amres@amOrthogonalPosF = -0.38             ; move the legend down
  annoid1 = gsn_add_annotation(plot(0),lbid,amres)   ; add legend to plot

;;;添加直线
   x_red1=(/0.37,0.37/)
   y_red1=(/0.49,0.72/)
   
   x_red2=(/0.75,0.75/)
   y_red2=(/0.49,0.75/)

   x_blue1=(/0.425,0.425/)
   y_blue1=(/0.49,0.73/)

   x_blue2=(/0.77,0.77/)
   y_blue2=(/0.49,0.76/)

    polyline@gsLineThicknessF=1.0
    polyline@gsLineDashPattern=2
    polyline@gsLineColor="red"
gsn_polyline_ndc(wks, x_red1, y_red1, polyline)
gsn_polyline_ndc(wks, x_red2, y_red2, polyline)
  polyline@gsLineColor="blue"
gsn_polyline_ndc(wks, x_blue1, y_blue1, polyline)
gsn_polyline_ndc(wks, x_blue2, y_blue2, polyline)

;;;添加地形
    ployres                  = True          ; polyline mods desiblue
    ployres@gsLineThicknessF = 16.0            ; line thickness   
    ployres@gsLineColor      = "black"         ; color of lines
    ployres@gsLineDashPattern=1
   ;---Attach the polyline
    ploylevel=1000
    ploylon2=fspan(100, 149, 4)
    ploylon3=fspan(280, 315, 25)
    ploylevel2=conform_dims(dimsizes(ploylon2), ploylevel, 0)
    ploylevel3=conform_dims(dimsizes(ploylon3), ploylevel, 0)  
     
  
    dum1=new(2,graphic)
    dum2=dum1
    dum3=dum1

    i=0
    ployres@gsLineDashPattern=16
    dum2(i)= gsn_add_polyline(wks,plot(i),ploylon2,ploylevel2,  ployres)
    ployres@gsLineDashPattern=0
    dum3(i)= gsn_add_polyline(wks,plot(i),ploylon3,ploylevel3,  ployres)

    ploylevel=-6
    ploylon2=fspan(100, 149, 4)
    ploylon3=fspan(280, 315, 25)
    ploylevel2=conform_dims(dimsizes(ploylon2), ploylevel, 0)
    ploylevel3=conform_dims(dimsizes(ploylon3), ploylevel, 0)  
    i=1
    ployres@gsLineDashPattern=16
    dum2(i)= gsn_add_polyline(wks,plot(i),ploylon2,ploylevel2,  ployres)
    ployres@gsLineDashPattern=0
    dum3(i)= gsn_add_polyline(wks,plot(i),ploylon3,ploylevel3,  ployres)
delete(res)
delete(res2)
delete(ploylon3)
delete(ploylon2)

delete(ploylevel2)
delete(ploylevel3)


;;*************************************************************************************************
;;;画图3
  diri="/mnt/f/share/rea_msf/"
    LeftString="mre"
    fils1=systemfunc("ls "+diri+"msf_time_"+LeftString+"_global.nc")
    ; fils1=systemfunc("ls "+diri+"/MERRA/MERRA/msf_time_MERRA_global.nc")
    ; fils1=systemfunc("ls "+diri+"MRE/msf_time_mre_global.nc")
    f1=addfile(fils1, "r")
    date=cd_calendar(f1->time, 1)
    timeind=ind(date.ge.197901.and.date.le.201212)
    levind1=ind(f1->level.eq.400)
    levind2=ind(f1->level.eq.600)
    msf_time_20cr=f1->msf_time(:,:,timeind)
    ud_time_20cr=f1->ud_time(:,:,timeind)
    w_time_20cr=f1->w_time(:,:,timeind)
    msf_time_20cr=(/msf_time_20cr/10^9/)
    
    delete(msf_time_20cr&time)
    msf_time_20cr&time=date(timeind)

    printVarSummary(ud_time_20cr)
    printVarSummary(msf_time_20cr)
   
    delete(fils1)
    delete(f1)
    delete(date)
    delete(timeind)
;;;
      ;;;先得到秋季和冬季eninio年的沃克环流情况SONDJF
      msf_son=month_to_season(msf_time_20cr(time|:,level|:,lon|:), "SON")
      msf_jja=month_to_season(msf_time_20cr(time|:,level|:,lon|:), "JJA")
      msf_djf=month_to_season(msf_time_20cr(time|:,level|:,lon|:), "DJF")
      msf_year=month_to_annual(msf_time_20cr(time|:,level|:,lon|:), 1)
      
      msf_son&time=ispan(1979,2012,1)
      msf_year&year=ispan(1979,2012,1)
      msf_jja&time=ispan(1979,2012,1)
      msf_djf&time=ispan(1979,2012,1)  ;;;因为这里第一个是1958 1月和2月的平均，第二个是195812,195901,195902的平均
      ;;注意点
      msf_djf(0:dimsizes(msf_djf&time)-2,:,:)=(/msf_djf(1:,:,:)/)
      msf_djf(dimsizes(msf_djf&time)-1,:,:)=(/msf_time_20cr(:,:,{201212})/)
      ; msf_djf=(/msf_jja/)
      ;;;
      printVarSummary(msf_djf)
      msf_clm_20cr=dim_avg_n_Wrap(month_to_season(msf_time_20cr(time|:,level|:,lon|:),"DJF"), 0)
      msf_djf_1=dim_avg_n_Wrap(msf_djf({1979:1994},:,:), 0)
      msf_djf_2=dim_avg_n_Wrap(msf_djf({1995:2012},:,:),0)

      msf_jja_1=dim_avg_n_Wrap(msf_jja({1979:1994},:,:), 0)
      msf_jja_2=dim_avg_n_Wrap(msf_jja({1995:2012},:,:),0)

      msf_year_1=dim_avg_n_Wrap(msf_year({1979:1994},:,:), 0)
      msf_year_2=dim_avg_n_Wrap(msf_year({1995:2012},:,:),0)


    ud_clm_20cr=dim_avg_n_Wrap(month_to_season(ud_time_20cr(time|:,level|:,lon|:),"DJF"), 0)
    w_clm_20cr=dim_avg_n_Wrap(month_to_season(w_time_20cr(time|:,level|:,lon|:),"DJF"), 0)

    
;************************************************
; 画图气候态
;************************************************
 
    res = True
    res@gsnDraw=False
    res@gsnFrame=False
    
   
    res@gsnMaximize=False

    res@cnLinesOn=True
    res@cnInfoLabelOn=False
    res@lbLabelBarOn=False
    res@cnLineLabelsOn=False
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

    res@cnLevelSelectionMode ="ExplicitLevels"
    res@cnLevels=(/-1.0/)

    res@cnFillOn=False
    res@cnLineColor="blue"
    res@cnMonoLineColor=True
    res@cnLineThicknessF =4.0


    res@vpHeightF=0.24
    res@vpWidthF=0.55
    res@vpXF=0.27
    res@vpYF=0.43

    res@cnLinesOn=False
     plot(2)=gsn_csm_pres_hgt(wks, msf_djf_1(:,:), res)
     res@cnLinesOn=True
     plot_blue1=gsn_csm_pres_hgt(wks, msf_djf_1(:,{0:180}), res)
     res@cnLevels=(/0/)
     plot_blue2=gsn_csm_pres_hgt(wks, msf_djf_1(:,{120:270}), res)
     delete(res@cnLevels)
     res@cnLevels=(/0/)
     plot_blue3=gsn_csm_pres_hgt(wks, msf_djf_1(:,{270:}), res)
     delete(res@cnLevels)
     res@cnLevels=(/-1/)
     plot_blue4=gsn_csm_pres_hgt(wks, msf_djf_1(:,{180:}), res)

     overlay(plot(2), plot_blue1)
     overlay(plot(2), plot_blue2)
     overlay(plot(2), plot_blue3)
     overlay(plot(2), plot_blue4)

    res@cnLineColor="red"
    res@cnLinesOn=False
    ; plot1=gsn_csm_pres_hgt(wks, msf_djf_2(:,:), res)
     res@cnLinesOn=True
     delete(res@cnLevels)
     res@cnLevels=(/-1.0/)
     plot_red1=gsn_csm_pres_hgt(wks, msf_djf_2(:,{0:180}), res)
     delete(res@cnLevels)
     res@cnLevels=(/0/)
     plot_red2=gsn_csm_pres_hgt(wks, msf_djf_2(:,{120:270}), res)
     delete(res@cnLevels)
     res@cnLevels=(/0/)
     plot_red3=gsn_csm_pres_hgt(wks, msf_djf_2(:,{270:}), res)
     delete(res@cnLevels)
     res@cnLevels=(/-1/)
     plot_red4=gsn_csm_pres_hgt(wks, msf_djf_2(:,{180:}), res)
     
     overlay(plot(2), plot_red1)
     overlay(plot(2), plot_red2)
     overlay(plot(2), plot_red3)
     overlay(plot(2), plot_red4)





 ;;---Resources for the polyline  添加下垫面的陆地
      ployres                  = True          ; polyline mods desired
      ployres@gsLineThicknessF = 15.0            ; line thickness   
      ployres@gsLineColor      = "black"         ; color of lines
      ployres@gsLineDashPattern=15
 ;---Attach the polyline
      ploylevel=1000
      ploylon1=fspan(9, 39, 25)
      ploylon2=fspan(100, 149, 25)
      ploylon3=fspan(280, 315, 25)


      ploylevel1=conform_dims(dimsizes(ploylon1), ploylevel, 0)
      ploylevel2=conform_dims(dimsizes(ploylon2), ploylevel, 0)
      ploylevel3=conform_dims(dimsizes(ploylon3), ploylevel, 0)
      
    str1 = unique_string("dum1")  ; "unique_string" will return a unique                                     ; string every time it is called fro                                 ;  within a single NCL session.
    str2 = unique_string("dum2")
    str3 = unique_string("dum3")
  
    ployres@gsLineDashPattern=0
     
    plot@$str1$= gsn_add_polyline(wks,plot(2),  ploylon1,ploylevel1,  ployres)

    ployres@gsLineDashPattern=16

    plot@$str2$= gsn_add_polyline(wks,plot(2),  ploylon2,ploylevel2,  ployres)

    ployres@gsLineDashPattern=0
    plot@$str3$= gsn_add_polyline(wks,plot(2),  ploylon3,ploylevel3,  ployres)
    
    text=(/"IWC","PWC","AWC"/)
    x=(/90,210,330/)
    y=(/400,400,400/)
   text1=new(3,graphic)
    textres=True
    textres@txFontThicknessF=0.1
    textres@txFontHeightF=0.02
    do i=0,2
    text1(i)=gsn_add_text(wks, plot(2), text(i), x(i), y(i), textres)
    end do

 ;------------------------------

;add Text
;------------------------------
txid_tr = new(3, graphic)
amid_tr = new(3, graphic)
txres                       = True
txres@txPerimOn             = True
txres@txFontHeightF         = 0.02

amres_tr                  = True
amres_tr@amParallelPosF   =  -0.46    ; This is the right edge of the plot.
amres_tr@amOrthogonalPosF =  0.4    ; This is the top edge of the plot.

letter = (/"a","b","c"/)
do i = 0, 2
;---Create text strings
    tr_label =  letter(i) 

    txres@txBackgroundFillColor = "white"
    txres@txPerimOn = True
    if i.eq.1 then
    amres_tr@amOrthogonalPosF =  0.32    ; This is the top edge of the plot.
  else
  amres_tr@amOrthogonalPosF =  0.4    ; This is the top edge of the plot.
end if
    txid_tr(i) = gsn_create_text(wks, tr_label, txres)

;---Attach text strings to plot
    amid_tr(i) = gsn_add_annotation(plot(i), txid_tr(i), amres_tr)

end do




     draw(plot)
     frame(wks)







end

