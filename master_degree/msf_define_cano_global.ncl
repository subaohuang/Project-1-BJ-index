
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
;;; 
undef("search_location")
function  search_location(msf_mean_annual)
local msf_mean_annual_xi , msf_mean_annual_dong
begin
locationxi=new(dimsizes(msf_mean_annual&time),typeof(msf_mean_annual&lon))
locationxi!0="time"
locationxi&time=ispan(0,dimsizes(msf_mean_annual&time)-1,1)
locationdong=locationxi
msf_mean_annual_xi=msf_mean_annual(:,{120:250})
msf_mean_annual_dong=msf_mean_annual(:,{230:280})
;;;judge locationxi
do t=0,dimsizes(msf_mean_annual&time)-1
    do i=0,dimsizes(msf_mean_annual_xi&lon)-9
        if  (msf_mean_annual_xi(t,i+1).ge.0).and.((msf_mean_annual_xi(t,i).le.0).and.(msf_mean_annual_xi(t,i+8).gt.0))  then 
         a=(msf_mean_annual_xi(t,i+1)-msf_mean_annual_xi(t,i))/(msf_mean_annual_xi&lon(i+1)-msf_mean_annual_xi&lon(i))
         locationxi(t)=tofloat(msf_mean_annual_xi&lon(i)- msf_mean_annual_xi(t,i)/a)
        ; locationxi(t)=(msf_mean_annual_xi&lon(i)+msf_mean_annual_xi&lon(i+1))/2
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
         a=(msf_mean_annual_dong(t,i+1)-msf_mean_annual_dong(t,i))/(msf_mean_annual_dong&lon(i+1)-msf_mean_annual_dong&lon(i))
         locationdong(t)=tofloat(msf_mean_annual_dong&lon(i)- msf_mean_annual_dong(t,i)/a)
        ; locationdong(t)=(msf_mean_annual_dong&lon(i)+msf_mean_annual_dong&lon(i+1))/2
        end if    
        if ismissing(locationdong(t)).ne.True  then 
        break
        end if
     end do
   ; if ismissing(locationdong(t)).eq.True  then 
   ; co_in=closest_val(0,msf_mean_annual_dong(t,:))      ；；；插值的问题比较明显
   ; locationdong(t)=msf_mean_annual_dong&lon(co_in)
   ; delete(co_in)
   ; end if

end do
locationxi=(/linmsg(locationxi, -1)/)
locationdong=(/linmsg(locationdong, -1)/)
 location = [/ locationxi, locationdong /]   ; return two variables as a type 'list'
return location
end 
;;;2d search
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
         locationxi= tofloat(msf_mean_annual_xi&lon(i)- msf_mean_annual_xi(i)/a)
        ; locationxi=(msf_mean_annual_xi&lon(i)+msf_mean_annual_xi&lon(i+1))/2
        end if    
        if  (msf_mean_annual_xi(i).eq.0)  then 
        locationxi=tofloat(msf_mean_annual_dong&xi(i))
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
        locationdong=tofloat(msf_mean_annual_dong&lon(i))
        end if    
        if ismissing(locationdong.ne.True)  then 
        break
        end if
    if ismissing(locationdong.eq.True)  then 
   co_in=closest_val(0,msf_mean_annual_dong)
   locationdong=tofloat(msf_mean_annual_dong&lon(co_in))
   delete(co_in)
   end if
     end do
locationxi=(/linmsg(locationxi, -1)/)
; locationdong=(/linmsg(locationdong, -1)/)
 location = [/ locationxi, locationdong /]   ; return two variables as a type 'list'
return location
end 
;************************************************
; 画20cr,ncep,era40及jRA-55的msf的变化趋势。
;************************************************
begin
fils1="rea_msf/msf_time_era_global.nc"
    diri="/mnt/f/share/"
    f1=addfile(diri+fils1, "r")
    date=cd_calendar(f1->time, 1)
    timeind=ind(date.ge.197901.and.date.le.201212)
    
    ; msf_time_20cr=f1->msf_time(:,{100:280},timeind)
    msf_time_20cr=f1->msf_time(:,:,timeind)
    msf_time_20cr=(/msf_time_20cr/10^9/)
    delete(msf_time_20cr&time)
    msf_time_20cr&time=date(timeind)
    delete(f1) 
    delete(date)
    delete(timeind)
    
    ;读入ERSST海温数据
    m=""+ispan(1979,2013,1)+""
    diri="/mnt/f/share/ersst/"
    fils=new((/dimsizes(m),12/),string)
    do i=0,dimsizes(m)-1
    fils2=systemfunc("ls "+diri+"ersst."+m(i)+"*"+".nc")
    fils(i,:)= fils2
    end do
    fils_new=onedtond(fils, dimsizes(m)*12)
    fileh=new(dimsizes(m)*12,string)
    fileh(0:dimsizes(m)*12-1)=fils_new
    f2=addfiles(fileh, "join")
    printVarSummary(f2)

    date=cd_calendar(f2[:]->time, 1)
    timeind=ind(date.le.201312.and.date.ge.197901)
    ; m=""+ispan(1979,2013,1)+""
    ; diri="/mnt/f/share/ersst/"
    ; fils1=systemfunc("ls "+diri+"sst.mnmean.v4.nc")
    ; ; fils1=systemfunc("ls "+diri+"HadISST_sst.nc")
    ; f2=addfile(fils1, "r")
    ; date=cd_calendar(f2->time, 1)
    ; timeind=ind(date.le.201312.and.date.ge.197901)
    ; date=cd_calendar(f2->time, 1)

    timeind=ind(date.le.201312.and.date.ge.197901)
    
    
    ; sst=f2->sst(timeind,{20:-20},:)
    sst=f2[:]->sst(timeind,0,{20:-20},:)

    ssta=(/rmMonAnnCycTLL(sst)/)
    copy_VarCoords(sst, ssta)
    ssta_de=dtrend_msg_n(ispan(1,dimsizes(ssta&time),1),ssta,False,True,0)
    copy_VarCoords(ssta, ssta_de)
    ssta_djf=month_to_season(ssta_de, "DJF")
    
    ssta_s=(/ssta_djf(1:,:,:)/)*1.0
    copy_VarCoords(ssta_djf(0:dimsizes(m)-2,:,:), ssta_s)

    emi=dim_avg_n(ssta_s(:,{10:-10},{165:220}),(/1,2/))-0.5*dim_avg_n(ssta_s(:,{20:-10},{125:145}), (/1,2/))-0.5*dim_avg_n(ssta_s(:,{5:-15},{250:290}),(/1,2/))
   
    emi_std=dim_standardize_n(emi, 1, 0)
    std_emi=dim_stddev_n(emi, 0)
   
    nino3=dim_avg_n(ssta_s(:,{5:-5},{210:270}), (/1,2/))
    nino3_std=dim_standardize_n(nino3, 1, 0)
    std_nino3=dim_stddev_n(nino3, 0)
    delete(diri)
    diri="/mnt/f/share/rea_msf/"

;;;

  
;;;去趋势
   ; msf_time_20cr_de=dtrend_msg_n(ispan(1,dimsizes(msf_time_20cr&time),1),msf_time_20cr,False,True,2)
   msf_time_20cr_de=msf_time_20cr
   copy_VarCoords(msf_time_20cr, msf_time_20cr_de)
 
;;;   
  
;;;先得到秋季和冬季eninio年的沃克环流情况SONDJF
 msf_son=month_to_season(msf_time_20cr_de(time|:,level|:,lon|:), "SON")
 msf_djf=month_to_season(msf_time_20cr_de(time|:,level|:,lon|:), "DJF")

 msf_son&time=ispan(1979,2012,1)
 msf_djf&time=ispan(1979,2012,1)  ;;;因为这里第一个是1958 1月和2月的平均，第二个是195812,195901,195902的平均
 ;;注意点
    msf_djf(0:dimsizes(msf_djf&time)-2,:,:)=(/msf_djf(1:,:,:)/)
    msf_djf(dimsizes(msf_djf&time)-1,:,:)=(/msf_time_20cr_de(:,:,{201212})/)
 ;;;

 ; msf_clm=dim_avg_n((/msf_son(:,:,:)+msf_djf(:,:,:)/)*0.5,0)
 msf_clm=dim_avg_n((/msf_djf(:,:,:)/)*1.0,0)
 copy_VarCoords(msf_son(1,:,:), msf_clm)
 msf_cold=msf_son
;;;
 msf_cold=(/msf_djf(:,:,:)*1.0/) 
 
 ; msf_cold=(/month_to_annual(msf_time_20cr_de(time|:,level|:,lon|:), 1)/)
 ; msf_djf=msf_cold
 msf_cold_mean=level_mean(msf_cold(:,{600:400},:))
 ; levind=ind(msf_cold&level.eq.500)
 ; msf_cold_mean=msf_cold(time|:,level|levind,lon|:)
 location=search_location(msf_cold_mean)
;;
;;;挑选年份 modoki year
 modokiind=ind(emi_std.ge.0.7)
 msf_modoki=dim_avg_n((/msf_djf(modokiind,:,:)+msf_djf(modokiind,:,:)/)*0.5,0)
 ; msf_modoki=dim_avg_n((/msf_djf(modokiind+1,:,:)/),0)   ;;;;djf取的要加1
 copy_VarCoords(msf_djf(1,:,:), msf_modoki)
;;;

;;;挑选年份 canino  enino year
 caind=ind(nino3_std.ge.1)
 msf_canonical=dim_avg_n((/msf_djf(caind,:,:)+msf_djf(caind,:,:)/)*0.5,0)
 ; msf_canonical=dim_avg_n((/msf_djf(caind+1,:,:)/),0)
 copy_VarCoords(msf_djf(1,:,:), msf_canonical)
;;;

;;;挑选年份 canino  lanina year
 laind=ind(nino3_std.le.(/-1/))
 msf_lanina=dim_avg_n((/msf_djf(laind,:,:)+msf_djf(laind,:,:)/)*0.5,0)
 ; msf_canonical=dim_avg_n((/msf_djf(caind+1,:,:)/),0)
 copy_VarCoords(msf_djf(1,:,:), msf_lanina)
;;; normal time
 noind=ind(nino3_std.ge.(/-0.5/).and.nino3_std.le.(/0.5/))
 msf_noind=dim_avg_n((/msf_djf(noind,:,:)+msf_djf(noind,:,:)/)*0.5,0)
 ; msf_canonical=dim_avg_n((/msf_djf(caind+1,:,:)/),0)
 copy_VarCoords(msf_djf(1,:,:), msf_noind)


;;;挑选年份 lanina modoki year
 lamodokiind=ind(emi_std.le.(/-1/))
 msf_lamodoki=dim_avg_n((/msf_djf(lamodokiind,:,:)+msf_djf(lamodokiind,:,:)/)*0.5,0)
 ; msf_canonical=dim_avg_n((/msf_djf(caind+1,:,:)/),0)
 copy_VarCoords(msf_djf(1,:,:), msf_lamodoki)
;;
 

;;;计算响应的结果，即用响应的nino指数做标准化
   nino3_can=avg(nino3(caind))
   nino3_lan=avg(nino3(laind))
   emi_modoki=avg(emi(modokiind))
;;;

lo_clm=search_location_2d(msf_clm({500},:))
width_clm=lo_clm[1]-lo_clm[0]


;;;画图
wks=gsn_open_wks("pdf", "/mnt/f/share/era/data/ERA_冬季各类事件的情况")

    res2 = True
    res2@gsnDraw=False
    res2@gsnFrame=False
    res2@gsnLeftString="Canonical Enino"

    res2@cnFillPalette  = "BlueDarkRed18"
    res2@cnLinesOn=True
    res2@cnLevelSelectionMode ="ManualLevels"
    res2@cnMinLevelValF = -6
    res2@cnMaxLevelValF = 6
    res2@cnLevelSpacingF = 1.0
    res2@cnFillOn=True


    res2@cnLineLabelsOn=False
    res2@cnInfoLabelOn=False
    

    res2@tiYAxisString=""
    res2@tmXTOn=False
    res2@tmYROn=False
    res2@tmXBMode="Explicit"
    res2@tmXBValues=(/0,90,180,270,357.5/)
    res2@tmXBLabels=(/"0","90~S~o~N~E","180","90~S~o~N~W","0"/)
    res2@tmXBMinorOn=True
    res2@tmXBMinorValues=ispan(0,360,10)

    ; res2@tmXBValues=ispan(100,310,40)
    ; res2@tmXBLabels=(/"100~S~o~N~E","140~S~o~N~E","180~S~o~N~E","140~S~o~N~W","100~S~o~N~W","60~S~o~N~W"/)
    ; res2@tmXBMinorOn=True
    ; res2@tmXBMinorValues=ispan(100,310,10)
    res2@tmYRMode="Automatic"
     
    res2@vpHeightF=0.35
    res2@vpWidthF=0.65
   
    plot=new(4,graphic)
    
    plot(0)=gsn_csm_pres_hgt(wks,msf_canonical,res2)
      res2@gsnLeftString="Enino modoki"
    ; res2@cnLineColor="burlywood4"
    plot(1)=gsn_csm_pres_hgt(wks,msf_modoki,res2)
      res2@gsnLeftString="Canonical Lanina"
    ; res2@cnLineColor="black"
    ; plot2=gsn_csm_pres_hgt(wks,msf_clm,res2)    
    ; res2@cnLineColor="blue"
    plot(2)=gsn_csm_pres_hgt(wks,msf_lanina,res2)  
      res2@gsnLeftString="Neutral"
    ; res2@cnLineColor="green"
    plot(3)=gsn_csm_pres_hgt(wks,msf_noind,res2)  

;;;添加地形
    ployres                  = True          ; polyline mods desired
    ployres@gsLineThicknessF = 8.0            ; line thickness   
    ployres@gsLineColor      = "black"         ; color of lines
    ployres@gsLineDashPattern=1
   ;---Attach the polyline
    ploylevel=1000
    ploylon1=fspan(9, 39, 25)
    ploylon2=fspan(100, 149, 10)
    ploylon3=fspan(280, 315, 25)
    ploylevel1=conform_dims(dimsizes(ploylon1), ploylevel, 0)
    ploylevel2=conform_dims(dimsizes(ploylon2), ploylevel, 0)
    ploylevel3=conform_dims(dimsizes(ploylon3), ploylevel, 0)  
 
    str1 = unique_string("dum1")  ; "unique_string" will return a unique                                     ; string every time it is called fro                                 ;  within a single NCL session.
    str2 = unique_string("dum2")
    str3 = unique_string("dum3")

    ployres@gsLineDashPattern=0
    do i=0,3
    plot@$str1$= gsn_add_polyline(wks,plot(i),ploylon1,ploylevel1,  ployres)
    ployres@gsLineDashPattern=16
    plot@$str2$= gsn_add_polyline(wks,plot(i),ploylon2,ploylevel2,  ployres)
    ployres@gsLineDashPattern=0
    plot@$str3$= gsn_add_polyline(wks,plot(i),ploylon3,ploylevel3,  ployres)
    str1 = unique_string("dum1")  ; "unique_string" will return a unique                                     ; string every time it is called fro                                 ;  within a single NCL session.
    str2 = unique_string("dum2")
    str3 = unique_string("dum3")

  end do


  pres=True
  pres@gsnPanelLabelBar=False
  pres@gsnPanelDebug=True
  gsn_panel(wks, plot, (/2,2/), pres)



end



