load "/mnt/d/ncl/lib/ncarg/nclscripts/csm/arrow.ncl"

;************************************************
; 画20cr,ncep,era的msf气候态
;************************************************
;;;
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
undef("draw_msf_trend")
function draw_msf_trend(LeftString,wks)
begin
    diri="/mnt/f/share/rea_msf/"
    f1=addfile(diri+"msf_time_"+LeftString+"_global.nc", "r")
    date=cd_calendar(f1->time, 1)
    timeind=ind(date.ge.197901.and.date.le.201212)
    level=f1->level
    levind=ind(f1->level.eq.500)
    
    w_time=f1->w_time(:,:,timeind)
    delete(w_time&time)
    w_time&time=date(timeind)

    delete(f1) 
    delete(diri)
    delete(date)
    delete(timeind)
;
    w_djf=month_to_season(w_time(time|:,level|:,lon|:), "DJF")
    w_djf&time=ispan(1979,2012,1)  ;;;因为这里第一个是1958 1月和2月的平均，第二个是195812,195901,195902的平均

    ;;注意点
    w_djf(0:dimsizes(w_djf&time)-2,:,:)=(/w_djf(1:,:,:)/)
    w_djf(dimsizes(w_djf&time)-1,:,:)=(/w_time(:,:,{201212})/)
  
    ;;;

    w_djf_500=w_djf(time|:,level|levind,lon|:)


;;;;水平插值
   xo=fspan(0,357.5,144)
   
    print("*****************开始水平插值************************")  

    w_djf_c=tofloat(linint1_n_Wrap(w_djf_500&lon, w_djf_500, \
            True, xo, 0,1))
    w_djf_c!1="lon"
    w_djf_c&lon=fspan(0,357.5,144)
   printVarSummary(w_djf_c)
    ; w_djf_c_pi=lonPivot(w_djf_c, 30)
    ; printVarSummary(w_djf_c_pi)
;;;画图
return [/w_djf_c/]
end

;***********************************************************
;主程序
;***********************************************************         
begin

wks=gsn_open_wks("pdf", "/mnt/f/share/era/data/rea_msf_compare_w_lon的图")


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
 i=0
 ListPush(lo7, (lo1[i]+lo2[i]+lo3[i]+lo4[i]+lo5[i]+lo6[i])/6.0)


w_rea=new((/6,34,144/),float)

w_rea(0,:,:)=lo1[0]
w_rea(1,:,:)=lo2[0]
w_rea(2,:,:)=lo3[0]
w_rea(3,:,:)=lo4[0]
w_rea(4,:,:)=lo5[0]
w_rea(5,:,:)=lo6[0]



w_rea!0="rea"
w_rea!1="time"
w_rea!2="lon"
w_rea&rea=ispan(1,6,1)
w_rea&time=ispan(1979,2012,1)
w_rea&lon=fspan(0,357.5,144)

w_rea1=w_rea(:,{1979:1994},{120:290})
w_rea2= w_rea(:,{1995:},{120:290})

w_rea1_mean=dim_avg_n_Wrap(w_rea1, (/0,1/))
w_rea2_mean=dim_avg_n_Wrap(w_rea2, (/0,1/))

;;;polygon_xy
  po_xy3=polygon_xy(dim_avg_n_Wrap(w_rea1,1))
  po_xy4=polygon_xy(dim_avg_n_Wrap(w_rea2,1))

;;;

res=True
res@gsnDraw=False
res@gsnFrame=False
res@gsnYRefLineDashPattern=1
res@gsnLeftString="(b)Vertical Wind at 500 hPa"
res@tmXTOn=False
res@tmYROn=False
res@tmXBMode="Explicit"
res@trXMaxF=290
res@trXMinF=120
res@tmXBMode="Explicit"
res@tmXBValues=ispan(120, 290, 30)
res@tmXBLabels=(/"120~S~o~N~E","150~S~o~N~E","180","150~S~o~N~W","120~S~o~N~W"\
                ,"90~S~o~N~W"/)
res@tmXBMinorOn=True
res@tmXBMinorValues=ispan(120,290,10)
res@vpHeightF=0.2
res@vpWidthF=0.7
res@xyDashPatterns=(/0,0/)
res@xyLineThicknesses=(/3,3/)
res@xyLineColors=(/"blue","deeppink"/)
;;;
;;;画图

plot=new(2,graphic)
plot(1)=gsn_csm_xy(wks, w_rea1_mean&lon,(/w_rea1_mean,w_rea2_mean/), res)

gsres                   = True                        ; poly res
gsres@gsFillColor       = "blue"                 ; color chosen
gsres@gsFillOpacityF    =0.1
gsres@tfPolyDrawOrder="PreDraw"
gsres@gsFillColor       = "blue"      
gsres@gsFillOpacityF    =0.1                      ; color chosen
dummy2 = gsn_add_polygon (wks,plot(1),po_xy3[0],po_xy3[1],gsres)
gsres@gsFillColor       = "deeppink"
gsres@gsFillOpacityF    =0.05
dummy3 = gsn_add_polygon (wks,plot(1),po_xy4[0],po_xy4[1],gsres)

polyline=True
polyline@gsLineThicknessF=2.0
polyline@tfPolyDrawOrder="PreDraw"
dummy5 = gsn_add_polyline(wks,plot(1),w_rea1_mean&lon,conform(w_rea1_mean&lon, 0, 0),polyline)

;;;

;;;添加地形
    ployres                  = True          ; polyline mods desired
    ployres@gsLineThicknessF = 16.0            ; line thickness   
    ployres@gsLineColor      = "black"         ; color of lines
    ployres@gsLineDashPattern=1
   ;---Attach the polyline
    ploylevel=-2
    ploylon1=fspan(9, 30, 25)
    ploylon2=fspan(100, 149, 25)
    ploylon3=fspan(280, 315, 25)
    ploylon4=fspan(30, 39, 25)
  

    ploylevel1=conform_dims(dimsizes(ploylon1), ploylevel, 0)
    ploylevel2=conform_dims(dimsizes(ploylon2), ploylevel, 0)
    ploylevel3=conform_dims(dimsizes(ploylon3), ploylevel, 0)
    ploylevel4=conform_dims(dimsizes(ploylon4), ploylevel, 0)
      

    str1 = unique_string("dum1")  ; "unique_string" will return a unique                                     ; string every time it is called fro                                 ;  within a single NCL session.
    str2 = unique_string("dum2")
    str3 = unique_string("dum3")
    str4 = unique_string("dum4")
    ployres@gsLineDashPattern=0
    plot@$str1$= gsn_add_polyline(wks,plot(1),  ploylon1,ploylevel1,  ployres)
    ployres@gsLineDashPattern=16
    plot@$str2$= gsn_add_polyline(wks,plot(1),  ploylon2,ploylevel2,  ployres)
    ployres@gsLineDashPattern=0
    plot@$str3$= gsn_add_polyline(wks,plot(1),  ploylon3,ploylevel3,  ployres)
    ployres@gsLineDashPattern=0
    plot@$str4$= gsn_add_polyline(wks,plot(1),  ploylon4,ploylevel4,  ployres)
    delete(str1)
    delete(str2)
    delete(str3)
    delete(str4)

;;;读入MSF mre的数据

    diri="/mnt/f/share/rea_msf/"
    LeftString="mre"
    fils1=systemfunc("ls "+diri+"msf_time_"+LeftString+"_global.nc")
    f1=addfile(fils1, "r")
    date=cd_calendar(f1->time, 1)
    timeind=ind(date.ge.197901.and.date.le.201212)
    levind1=ind(f1->level.eq.400)
    levind2=ind(f1->level.eq.600)
    msf_time_20cr_r=f1->msf_time(:,:,timeind)


    msf_time_20cr_pi=lonPivot(msf_time_20cr_r(time|:,level|:,lon|:), 30)
    msf_time_20cr=msf_time_20cr_pi(level|:,lon|:,time|:)
 
    msf_time_20cr=(/msf_time_20cr/10^9/)
    delete(msf_time_20cr&time)
    msf_time_20cr&time=date(timeind)
    printVarSummary(msf_time_20cr)
   
    delete(fils1)
    delete(f1)
    delete(date)
    delete(timeind)


;;; 
      ;;;先得到秋季和冬季eninio年的沃克环流情况SONDJF

      msf_djf=month_to_season(msf_time_20cr(time|:,level|:,lon|:), "DJF")
      msf_year=month_to_annual(msf_time_20cr(time|:,level|:,lon|:), 1)
      msf_djf&time=ispan(1979,2012,1)  ;;;因为这里第一个是1958 1月和2月的平均，第二个是195812,195901,195902的平均
      ;;注意点
      msf_djf(0:dimsizes(msf_djf&time)-2,:,:)=(/msf_djf(1:,:,:)/)
      msf_djf(dimsizes(msf_djf&time)-1,:,:)=(/msf_time_20cr(:,:,{201212})/)
      ;;;
      printVarSummary(msf_djf)
      msf_clm_20cr=dim_avg_n_Wrap(month_to_season(msf_time_20cr(time|:,level|:,lon|:),"DJF"), 0)
      msf_djf_1=dim_avg_n_Wrap(msf_djf({1979:1994},:,:), 0)
      msf_djf_2=dim_avg_n_Wrap(msf_djf({1995:2012},:,:),0)
    
;************************************************
; 画图气候态
;************************************************
  





    res2 = True
    res2@gsnDraw=False
    res2@gsnFrame=False
    res2@gsnMaximize=False
    res2@cnLinesOn=True
    res2@cnInfoLabelOn=False
    res2@lbLabelBarOn=False
    res2@cnLineLabelsOn=False
    res2@tmYRMode="Automatic"
    res2@gsnLeftString="(a)Mass streamfunction"
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
    res2@vpHeightF=0.4
    res2@vpWidthF=0.7

    plot(0)=gsn_csm_pres_hgt(wks, msf_djf_1(:,{120:290}), res2)
    res2@cnLineColor="red"
    plot1=gsn_csm_pres_hgt(wks, msf_djf_2(:,{120:290}), res2)
    overlay(plot(0), plot1)

    res2@cnLevelSelectionMode ="ExplicitLevels"
    delete(res2@cnLevels)
    res2@cnLevels=(/0/)
    res2@cnLineThicknessF =4.0
    res2@cnLineDashPattern=0
    res2@cnLineColor="blue"
    plot3=gsn_csm_pres_hgt(wks, msf_djf_1(:,{120:290}), res2)
    res2@cnLineColor="red"
    plot4=gsn_csm_pres_hgt(wks, msf_djf_2(:,{120:290}), res2)
    overlay(plot(0), plot3)
    overlay(plot(0), plot4)


  


 ;;---Resources for the polyline  添加下垫面的陆地
      ployres                  = True          ; polyline mods desired
      ployres@gsLineThicknessF = 15.0            ; line thickness   
      ployres@gsLineColor      = "black"         ; color of lines
      ployres@gsLineDashPattern=15
 ;---Attach the polyline
      ploylevel=1000

      ploylon1=fspan(9, 30, 25)
      ploylon2=fspan(100, 149, 25)
      ploylon3=fspan(280, 315, 25)
      ploylon4=fspan(30, 39, 25)

      ploylevel1=conform_dims(dimsizes(ploylon1), ploylevel, 0)
      ploylevel2=conform_dims(dimsizes(ploylon2), ploylevel, 0)
      ploylevel3=conform_dims(dimsizes(ploylon3), ploylevel, 0)
      ploylevel4=conform_dims(dimsizes(ploylon4), ploylevel, 0)
      




      str1 = unique_string("dum1")  ; "unique_string" will return a unique                                     ; string every time it is called fro                                 ;  within a single NCL session.
      str2 = unique_string("dum2")
      str3 = unique_string("dum3")
      str4 = unique_string("dum4")
      ployres@gsLineDashPattern=0

      plot@$str1$= gsn_add_polyline(wks,plot(0),  ploylon1,ploylevel1,  ployres)

      ployres@gsLineDashPattern=16

      plot@$str2$= gsn_add_polyline(wks,plot(0),  ploylon2,ploylevel2,  ployres)

      ployres@gsLineDashPattern=0
      plot@$str3$= gsn_add_polyline(wks,plot(0),  ploylon3,ploylevel3,  ployres)

      ployres@gsLineDashPattern=0
      plot@$str4$= gsn_add_polyline(wks,plot(0),  ploylon4,ploylevel4,  ployres)

      text=(/"IWC","PWC","AWC"/)
      x=(/90,210,330/)
      y=(/400,400,400/)
      text1=new(3,graphic)
      textres=True
      textres@txFontThicknessF=0.1
      textres@txFontHeightF=0.02
      ; do i=0,2
      ; text1(i)=gsn_add_text(wks, plot(0), text(i), x(i), y(i), textres)
      ; end do

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
      amres@amParallelPosF   = -0.38         ; move legend to the right
      amres@amOrthogonalPosF = -0.38             ; move the legend down
      annoid1 = gsn_add_annotation(plot(0),lbid,amres)   ; add legend to plot

pres=True
gsn_panel(wks, plot, (/2,1/), pres)
    
 end