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
return [/msf_djf_c,w_djf_c/]
end

;***********************************************************
;主程序
;***********************************************************         
begin


wks=gsn_open_wks("pdf", "/mnt/f/share/era/data/rea_lon的section图")


;;;画图
lo1=draw_msf_trend("20cr",wks)
lo2=draw_msf_trend("era",wks)
lo3=draw_msf_trend("jra55",wks)
lo4=draw_msf_trend("MERRA",wks)
lo5=draw_msf_trend("ncep1",wks)
lo6=draw_msf_trend("ncep2",wks)
; 
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

msf_rea1=msf_rea(:,{1979:1994},{100:280})
msf_rea2= msf_rea(:,{1995:},{100:280})



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

; w_rea1=w_rea(:,{1979:1994},{100:280})
; w_rea2= w_rea(:,{1995:},{100:280})

w_rea1=w_rea(:,{1979:1994},:)
w_rea2= w_rea(:,{1995:},:)


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
res@gsnLeftString="(a)400-600hpa average of Zonal Streamfunction"
res@tmXTOn=False
res@tmYROn=False
res@tmXBMode="Explicit"
; res@trXMaxF=280
; res@trXMinF=100
; res@trYMaxF=30
; res@trYMinF=-30
; res@tmXBValues=ispan(1,34,10)
; res@tmXBLabels=ispan(1979,2012,10)
; res@tmXBMinorValues=ispan(1,34,1)
res@vpHeightF=0.2
res@vpWidthF=0.7
res@xyDashPatterns=(/0,0/)
res@xyLineThicknesses=(/3,3/)
res@xyLineColors=(/"blue","deeppink"/)


; res@tfPolyDrawOrder="PostDraw"
;;;
;;;画图

plot=new(2,graphic)
plot(0)=gsn_csm_xy(wks, msf_rea1_mean&lon,(/msf_rea1_mean,msf_rea2_mean/), res)
res@gsnLeftString="(b)Vertical Wind at 500 hPa"
plot(1)=gsn_csm_xy(wks, w_rea1_mean&lon,(/w_rea1_mean,w_rea2_mean/), res)

gsres                   = True                        ; poly res
gsres@gsFillColor       = "blue"                 ; color chosen
gsres@gsFillOpacityF    =0.1
gsres@tfPolyDrawOrder="PreDraw"
dummy = gsn_add_polygon (wks,plot(0),po_xy1[0],po_xy1[1],gsres)
gsres@gsFillOpacityF    =0.05
gsres@gsFillColor       = "deeppink"   
dummy1 = gsn_add_polygon (wks,plot(0),po_xy2[0],po_xy2[1],gsres)
 gsres@gsFillColor       = "blue"      
 gsres@gsFillOpacityF    =0.1                      ; color chosen
dummy2 = gsn_add_polygon (wks,plot(1),po_xy3[0],po_xy3[1],gsres)
gsres@gsFillColor       = "deeppink"
gsres@gsFillOpacityF    =0.05
dummy3 = gsn_add_polygon (wks,plot(1),po_xy4[0],po_xy4[1],gsres)

polyline=True
polyline@gsLineThicknessF=2.0
polyline@tfPolyDrawOrder="PreDraw"
dummy4 = gsn_add_polyline(wks,plot(0),msf_rea1_mean&lon,conform(msf_rea1_mean&lon, 0, 0),polyline)
dummy5 = gsn_add_polyline(wks,plot(1),w_rea1_mean&lon,conform(w_rea1_mean&lon, 0, 0),polyline)

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
dummy8 = gsn_csm_xy(wks,(/w_rea1_xi,w_rea1_dong/), ispan(-2,4,1),res1)
res1@xyLineColors=(/"deeppink","deeppink"/)
dummy9 = gsn_csm_xy(wks,(/w_rea2_xi,w_rea2_dong/), ispan(-2,4,1),res1)
overlay(plot(0), dummy6)
overlay(plot(0), dummy7)
overlay(plot(1), dummy8)
overlay(plot(1), dummy9)

;;;

;;;添加地形
    ployres                  = True          ; polyline mods desired
    ployres@gsLineThicknessF = 16.0            ; line thickness   
    ployres@gsLineColor      = "black"         ; color of lines
    ployres@gsLineDashPattern=1
   ;---Attach the polyline
    ploylevel=-6
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
    
    ploylevel=-2
    ploylon2=fspan(100, 149, 4)
    ploylon3=fspan(280, 315, 25)
    ploylevel2=conform_dims(dimsizes(ploylon2), ploylevel, 0)
    ploylevel3=conform_dims(dimsizes(ploylon3), ploylevel, 0)  
    
    i=1
    ployres@gsLineDashPattern=16
    dum2(i)= gsn_add_polyline(wks,plot(i),ploylon2,ploylevel2,  ployres)
    ployres@gsLineDashPattern=0
    dum3(i)= gsn_add_polyline(wks,plot(i),ploylon3,ploylevel3,  ployres)
      

pres=True
gsn_panel(wks, plot, (/2,1/), pres)

end

