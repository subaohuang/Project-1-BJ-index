
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
msf_mean_annual_xi=msf_mean_annual(:,{120:220})
msf_mean_annual_dong=msf_mean_annual(:,{230:280})
;;;judge locationxi
do t=0,dimsizes(msf_mean_annual&time)-1
    do i=0,dimsizes(msf_mean_annual_xi&lon)-13
        if  (msf_mean_annual_xi(t,i+1).ge.0).and.((msf_mean_annual_xi(t,i).le.0).and.(msf_mean_annual_xi(t,i+12).gt.0))  then 
         a=(msf_mean_annual_xi(t,i+1)-msf_mean_annual_xi(t,i))/(msf_mean_annual_xi&lon(i+1)-msf_mean_annual_xi&lon(i))
         locationxi(t)=dble2flt((/msf_mean_annual_xi&lon(i)- msf_mean_annual_xi(t,i)/a/))
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
         locationdong(t)=dble2flt((/msf_mean_annual_dong&lon(i)- msf_mean_annual_dong(t,i)/a/))
        ; locationdong(t)=(msf_mean_annual_dong&lon(i)+msf_mean_annual_dong&lon(i+1))/2
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
;;;;二维的位置查找
undef("search_location_2d")
function  search_location_2d(msf_mean_annual)
begin
msf_mean_annual_xi=msf_mean_annual({120:220})
msf_mean_annual_dong=msf_mean_annual({230:280})
locationxi=msf_mean_annual@_FillValue
locationdong=locationxi
;;;judge locationxi
    do i=0,dimsizes(msf_mean_annual_xi&lon)-13
        if  (msf_mean_annual_xi(i+1).ge.0).and.((msf_mean_annual_xi(i).le.0).and.(msf_mean_annual_xi(i+12).gt.0))  then 
         a=(msf_mean_annual_xi(i+1)-msf_mean_annual_xi(i))/(msf_mean_annual_xi&lon(i+1)-msf_mean_annual_xi&lon(i))
         locationxi=msf_mean_annual_xi&lon(i)- msf_mean_annual_xi(i)/a
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
         locationdong=msf_mean_annual_dong&lon(i)- msf_mean_annual_dong(i)/a    
        ; locationdong=(msf_mean_annual_dong&lon(i)+msf_mean_annual_dong&lon(i+1))/2
        end if    
         if  (msf_mean_annual_dong(i).eq.0)  then 
        locationdong=msf_mean_annual_dong&lon(i)
        end if    
        if ismissing(locationdong.ne.True)  then 
        break
        end if
     end do
 location = [/ locationxi, locationdong /]   ; return two variables as a type 'list'
return location
end 


;************************************************
; 画20cr,ncep,era40及jRA-55的msf的变化趋势。
;************************************************
undef("draw_msf_trend")
function draw_msf_trend(fils1,wks,LeftString)
begin
    diri="/mnt/f/share/"
    f1=addfile(fils1, "r")
    date=cd_calendar(f1->time, 1)
    timeind=ind(date.ge.195801.and.date.le.201212)
    
    msf_time_20cr=f1->msf_time(:,{100:310},timeind)
    msf_time_20cr=(/msf_time_20cr/10^9/)
    delete(msf_time_20cr&time)
    msf_time_20cr&time=date(timeind)
    delete(f1) 
    delete(date)
    delete(timeind)
    
    f2=addfile("/mnt/f/share/ersst/HadISST_sst.nc", "r")
    date=cd_calendar(f2->time, 1)
    timeind=ind(date.ge.195801.and.date.le.201312)
    sst=lonFlip(f2->sst(timeind,{20:-20},:))
    ssta=(/rmMonAnnCycTLL(sst)/)
    copy_VarCoords(sst, ssta)
    ssta_de=dtrend_msg_n(ispan(1,dimsizes(ssta&time),1),ssta,False,True,0)
    copy_VarCoords(ssta, ssta_de)
    ssta_son=month_to_season(ssta_de, "SON")
    ssta_djf=month_to_season(ssta_de, "DJF")
    
    
    ; ssta_s=(/ssta_son(0:54,:,:)+ssta_djf(1:,:,:)/)*0.5
    ssta_s=(/ssta_djf(1:,:,:)/)*0.5
    copy_VarCoords(ssta_son(0:54,:,:), ssta_s)

    emi=dim_avg_n(ssta_s(:,{10:-10},{165:220}),(/1,2/))-0.5*dim_avg_n(ssta_s(:,{20:-10},{125:145}), (/1,2/))-0.5*dim_avg_n(ssta_s(:,{5:-15},{250:290}),(/1,2/))
   
    emi_std=dim_standardize_n(emi, 1, 0)
    std_emi=dim_stddev_n(emi, 0)

    nino3=dim_avg_n(ssta_s(:,{5:-5},{210:270}), (/1,2/))
    nino3_std=dim_standardize_n(nino3, 1, 0)
    std_nino3=dim_stddev_n(nino3, 0)
 
;;;

  
;;;去趋势
   ; msf_time_20cr_de=dtrend_msg_n(ispan(1,dimsizes(msf_time_20cr&time),1),msf_time_20cr,False,True,2)
   msf_time_20cr_de=msf_time_20cr
   copy_VarCoords(msf_time_20cr, msf_time_20cr_de)
 
;;;   
  
 
;;;先得到秋季和冬季eninio年的沃克环流情况SONDJF
 msf_son=month_to_season(msf_time_20cr_de(time|:,level|:,lon|:), "SON")
 msf_djf=month_to_season(msf_time_20cr_de(time|:,level|:,lon|:), "DJF")

 msf_son&time=ispan(1958,2012,1)
 msf_djf&time=ispan(1958,2012,1)  ;;;因为这里第一个是1958 1月和2月的平均，第二个是195812,195901,195902的平均
 ;;注意点
    msf_djf(0:dimsizes(msf_djf&time)-2,:,:)=(/msf_djf(1:,:,:)/)
    msf_djf(dimsizes(msf_djf&time)-1,:,:)=(/msf_time_20cr_de(:,:,{201212})/)
 ;;;

 ; msf_clm=dim_avg_n((/msf_son(:,:,:)+msf_djf(:,:,:)/)*0.5,0)
 msf_clm=dim_avg_n((/msf_djf(:,:,:)/)*1.0,0)
 copy_VarCoords(msf_son(1,:,:), msf_clm)
 msf_cold=msf_son
 ; msf_cold=(/(msf_son(:,:,:)+msf_djf(:,:,:))*0.5/)
;;;
msf_cold=(/msf_djf(:,:,:))*1.0) 
 ; msf_cold_mean=level_mean(msf_cold(:,{600:400},:))

 levind=ind(msf_cold&level.eq.500)
 msf_cold_mean=msf_cold(time|:,level|levind,lon|:)
 location=search_location(msf_cold_mean)
;;

;;输出每个数据location的位置
opt=diri+LeftString+"_msf_cold_locationxi_node_index.ascii"    
asciiwrite(opt, sprintf ("%6.2f,",location[0]))  

opt1=diri+LeftString+"_msf_cold_locationdong_node_index.ascii"    
asciiwrite(opt1, sprintf ("%6.2f,",location[1]))  



; ;;;输出SONDJF每个月的位置情况
;     msf_de_sondjf=new((/dimsizes(msf_time_20cr_de&level),dimsizes(msf_time_20cr_de&lon),328/),typeof(msf_time_20cr_de))
;     do i=0,53
;     msf_de_sondjf(:,:,i*6:5+6*i)=msf_time_20cr_de(:,:,8+12*i:13+12*i)
;     end do
;     msf_de_sondjf(:,:,324:327)=msf_time_20cr_de(:,:,656:659)
 

;     msf_de_sondjf!0="level"
;     msf_de_sondjf!1="lon"
;     msf_de_sondjf!2="time"
;     msf_de_sondjf&time=ispan(1,328,1)
;     msf_de_sondjf&level=msf_time_20cr_de&level
;     msf_de_sondjf&lon=msf_time_20cr_de&lon
   
;     printVarSummary(msf_de_sondjf)
;     msf_cold_month=level_mean(msf_de_sondjf({600:400},:,:))
;     location_month=search_location(msf_cold_month)
  
;     opt2=diri+LeftString+"_msf_cold_locationxi_month_index.ascii"   
;     asciiwrite(opt2, sprintf ("%6.2f,",location_month[0]))  
;     opt3=diri+LeftString+"_msf_cold_locationdong_month_index.ascii"  
;     asciiwrite(opt3, sprintf ("%6.2f,",location_month[1]))  
; ;;;


;;;挑选年份 modoki year
 modokiind=ind(emi_std.ge.std_emi)
 msf_modoki=dim_avg_n((/msf_son(modokiind,:,:)+msf_djf(modokiind,:,:)/)*0.5,0)
 ; msf_modoki=dim_avg_n((/msf_djf(modokiind+1,:,:)/),0)   ;;;;djf取的要加1
 copy_VarCoords(msf_son(1,:,:), msf_modoki)
;;;

;;;挑选年份 canino  enino year
 caind=ind(nino3_std.ge.std_nino3)
 msf_canonical=dim_avg_n((/msf_son(caind,:,:)+msf_djf(caind,:,:)/)*0.5,0)
 ; msf_canonical=dim_avg_n((/msf_djf(caind+1,:,:)/),0)
 copy_VarCoords(msf_son(1,:,:), msf_canonical)
;;;

;;;挑选年份 canino  lanina year
 laind=ind(nino3_std.le.(/-1*std_nino3/))
 msf_lanina=dim_avg_n((/msf_son(laind,:,:)+msf_djf(laind,:,:)/)*0.5,0)
 ; msf_canonical=dim_avg_n((/msf_djf(caind+1,:,:)/),0)
 copy_VarCoords(msf_son(1,:,:), msf_lanina)
;;; normal time
 noind=ind(nino3_std.ge.(/-0.5*std_nino3/).and.nino3_std.le.(/0.5*std_nino3/))
 msf_noind=dim_avg_n((/msf_son(noind,:,:)+msf_djf(noind,:,:)/)*0.5,0)
 ; msf_canonical=dim_avg_n((/msf_djf(caind+1,:,:)/),0)
 copy_VarCoords(msf_son(1,:,:), msf_noind)

;;;挑选年份 lanina modoki year
 lamodokiind=ind(emi_std.le.(/-1*std_emi/))
 msf_lamodoki=dim_avg_n((/msf_son(lamodokiind,:,:)+msf_djf(lamodokiind,:,:)/)*0.5,0)
 ; msf_canonical=dim_avg_n((/msf_djf(caind+1,:,:)/),0)
 copy_VarCoords(msf_son(1,:,:), msf_lamodoki)
;;
 

;;;计算响应的结果，即用响应的nino指数做标准化
   nino3_can=avg(nino3(caind))
   nino3_lan=avg(nino3(laind))
   emi_modoki=avg(emi(modokiind))
;;;

lo_clm=search_location_2d(msf_clm({500},:))
width_clm=lo_clm[1]-lo_clm[0]

if  fils1.eq."MRE/msf_time_mre_global.nc" then 
print("******cano***********")
 lo=search_location_2d(msf_canonical({500},:))
 rsp1=(lo[0]-lo_clm[0])/nino3_can
 rsp2=(lo[1]-lo_clm[1])/nino3_can
 rsp3=(lo[1]-lo[0]-width_clm)/nino3_can

 print("lo_xi="+lo[0])
 print("lo_dong="+lo[1])
 print("rsp1="+rsp1)
 print("rsp2="+rsp2)
 print("rsp3="+rsp3)
 delete(lo)
 delete(rsp1)
 delete(rsp2)
 delete(rsp3)
print("******lanina***********")
lo=search_location_2d(msf_lanina({500},:))
 rsp1=(lo[0]-lo_clm[0])/nino3_lan
 rsp2=(lo[1]-lo_clm[1])/nino3_lan
 rsp3=(lo[1]-lo[0]-width_clm)/nino3_lan
 print("lo_xi="+lo[0])
 print("lo_dong="+lo[1])
 print("rsp1="+rsp1)
 print("rsp2="+rsp2)
print("rsp3="+rsp3)
delete(lo)
delete(rsp1)
delete(rsp2)
delete(rsp3)
print("******modoki***********")
lo=search_location_2d(msf_modoki({500},:))
 rsp1=(lo[0]-lo_clm[0])/emi_modoki
 rsp2=(lo[1]-lo_clm[1])/emi_modoki
 rsp3=(lo[1]-lo[0]-width_clm)/emi_modoki
 print("rsp1="+rsp1)
 print("rsp2="+rsp2)
print("rsp3="+rsp3)
 print("lo_xi="+lo[0])
 print("lo_dong="+lo[1])

 delete(lo)
  delete(rsp1)
 delete(rsp2)
 delete(rsp3)

 print("******clm***********")
lo=search_location_2d(msf_clm({500},:))
 print("lo_xi="+lo[0])
 print("lo_dong="+lo[1])
 delete(lo)
  
end if


;;;画图

    res2 = True
    res2@gsnDraw=False
    res2@gsnFrame=False
    res2@gsnLeftString=LeftString
    res2@cnLinesOn=True
    res2@cnLevelSelectionMode ="ExplicitLevels"
    res2@cnLevels=0.0
    res2@cnFillOn=False
    res2@cnLineColor="red"
    res2@cnMonoLineColor=True
    res2@cnLineThicknesses=4.0
    res2@cnLineLabelsOn=False
    res2@cnInfoLabelOn=False
    

    res2@tiYAxisString=""
    res2@tmXTOn=False
    res2@tmYROn=False
    res2@tmXBMode="Explicit"
  
    ; res2@tmXBValues=(/0,90,180,270,357.5/)
    ; res2@tmXBLabels=(/"0","90~S~o~N~E","180","90~S~o~N~W","0"/)
    ; res2@tmXBMinorOn=True
    ; res2@tmXBMinorValues=ispan(0,360,10)

    res2@tmXBValues=ispan(100,310,40)
    res2@tmXBLabels=(/"100~S~o~N~E","140~S~o~N~E","180~S~o~N~E","140~S~o~N~W","100~S~o~N~W","60~S~o~N~W"/)
    res2@tmXBMinorOn=True
    res2@tmXBMinorValues=ispan(100,310,10)
    res2@tmYRMode="Automatic"
     
    res2@vpHeightF=0.3
    res2@vpWidthF=0.7
   
    plot=new(1,graphic)
    
    plot(0)=gsn_csm_pres_hgt(wks,msf_canonical,res2)
    
    res2@cnLineColor="green"

    plot1=gsn_csm_pres_hgt(wks,msf_modoki,res2)
    ; res2@cnLineColor="black"
    ; plot2=gsn_csm_pres_hgt(wks,msf_clm,res2)    
    res2@cnLineColor="blue"
    plot2=gsn_csm_pres_hgt(wks,msf_lanina,res2)  
    res2@cnLineColor="yellow"
    plot3=gsn_csm_pres_hgt(wks,msf_noind,res2)  
    overlay(plot,plot1)
    overlay(plot,plot2)
    overlay(plot,plot3)
    
    
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
    plot@$str1$= gsn_add_polyline(wks,plot,ploylon1,ploylevel1,  ployres)
    ployres@gsLineDashPattern=16
    plot@$str2$= gsn_add_polyline(wks,plot,ploylon2,ploylevel2,  ployres)
    ployres@gsLineDashPattern=0
    plot@$str3$= gsn_add_polyline(wks,plot,ploylon3,ploylevel3,  ployres)

return plot
end

;***********************************************************
;主程序
;***********************************************************         
begin

; fils1="ncep/20CR/msf_time_20cr_global_6--6.nc"
; fils2="ncep/ncep1/msf_time_ncep1_global.nc"
; ; fils3="era/era40/msf_time_era40_in_global.nc"
; fils3="era/era40/msf_time_in_era40_global.nc"
; fils4="JRA-55/msf_time_jra55_global.nc"
; fils5="MRE/msf_time_mre_global.nc"


fils1="rea_msf/msf_time_20cr_global.nc"
fils2="rea_msf/msf_time_ncep1_global.nc"
; fils3="era/era40/msf_time_era40_in_global.nc"
fils3="rea_msf/msf_time_era_global.nc"
fils4="rea_msf/msf_time_jra55_global.nc"
fils5="rea_msf/msf_time_ncep2_global.nc"
fils6="rea_msf/msf_time_merra_global.nc"
exit()


wks=gsn_open_wks("pdf", "/mnt/f/share/era/data/4个再分析资料_自定义的modoki位置秋冬")

plot=new(5,graphic)
plot(0)=draw_msf_trend(fils1,wks,"20CR")
plot(1)=draw_msf_trend(fils2,wks,"NCPE1")
plot(2)=draw_msf_trend(fils3,wks,"ERA")
plot(3)=draw_msf_trend(fils4,wks,"JRA-55")
plot(4)=draw_msf_trend(fils5,wks,"MRE")



    pres=True
    pres@gsnPanelLabelBar=False
    gsn_panel(wks, plot(0:3), (/2,2/), pres)
end



