
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
undef("search_location_alt") ;;;altlantic
function  search_location_alt(msf_mean_annual)
local msf_mean_annual_xi , msf_mean_annual_dong
begin
locationxi=new(dimsizes(msf_mean_annual&time),typeof(msf_mean_annual&lon))
locationxi!0="time"
locationxi&time=ispan(0,dimsizes(msf_mean_annual&time)-1,1)
locationdong=locationxi
msf_mean_annual_xi=msf_mean_annual(:,{-90:-30})
msf_mean_annual_dong=msf_mean_annual(:,{-40:10})
;;;judge locationxi
do t=0,dimsizes(msf_mean_annual&time)-1
    do i=0,dimsizes(msf_mean_annual_xi&lon)-9
        if  (msf_mean_annual_xi(t,i+1).ge.0).and.((msf_mean_annual_xi(t,i).le.0).and.(msf_mean_annual_xi(t,i+2).gt.0))  then 
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
   if ismissing(locationxi(t)).eq.True  then 
   co_in=closest_val(max(msf_mean_annual(t,:))/10,msf_mean_annual_dong(t,:))
   test=msf_mean_annual_xi(t,{-70:-30})
   co_in=closest_val(0,msf_mean_annual_xi(t,{-70:-30}))
   locationxi(t)=test&lon(co_in)
   delete(co_in)
  end if

end do
; locationxi=(/linmsg(locationxi, -1)/)
locationdong=(/linmsg(locationdong, -1)/)
 location = [/ locationxi, locationdong /]   ; return two variables as a type 'list'
return location
end 

;;; 
undef("search_location_id") ;;;altlantic
function  search_location_id(msf_mean_annual)
local msf_mean_annual_xi , msf_mean_annual_dong
begin
locationxi=new(dimsizes(msf_mean_annual&time),typeof(msf_mean_annual&lon))
locationxi!0="time"
locationxi&time=ispan(0,dimsizes(msf_mean_annual&time)-1,1)
locationdong=locationxi
msf_mean_annual_xi=msf_mean_annual(:,{35:55})
msf_mean_annual_dong=msf_mean_annual(:,{100:220})
msf_mean_annual_xi=(/msf_mean_annual_xi+0.5/)
msf_mean_annual_dong=(/msf_mean_annual_dong+0.5/)


;;;judge locationxi
do t=0,dimsizes(msf_mean_annual&time)-1
   co_in=closest_val(0,msf_mean_annual_xi(t,:))
   locationxi(t)=msf_mean_annual_xi&lon(co_in)
   delete(co_in)
end do

;;;judge locationdong
do t=0,dimsizes(msf_mean_annual&time)-1
    do i=0,dimsizes(msf_mean_annual_dong&lon)-2
        ; if  (msf_mean_annual_dong(t,i+1).le.0).and.((msf_mean_annual_dong(t,i).ge.0).and.(msf_mean_annual_dong(t,i+3).lt.0))  then 
         if  (msf_mean_annual_dong(t,i+1).ge.0).and.((msf_mean_annual_dong(t,i).le.0))  then 
         a=(msf_mean_annual_dong(t,i+1)-msf_mean_annual_dong(t,i))/(msf_mean_annual_dong&lon(i+1)-msf_mean_annual_dong&lon(i))
         locationdong(t)=tofloat(msf_mean_annual_dong&lon(i)- msf_mean_annual_dong(t,i)/a)
        ; locationdong(t)=(msf_mean_annual_dong&lon(i)+msf_mean_annual_dong&lon(i+1))/2
        end if    
        if ismissing(locationdong(t)).ne.True  then 
        break
        end if
     end do
  ;  if ismissing(locationxi(t)).eq.True  then 
  ;  co_in=closest_val(max(msf_mean_annual(t,:))/10,msf_mean_annual_dong(t,:))
  ;  test=msf_mean_annual_xi(t,{-70:-30})
  ;  co_in=closest_val(0,msf_mean_annual_xi(t,{-70:-30}))
  ;  locationxi(t)=test&lon(co_in)
  ;  delete(co_in)
  ; end if
end do
; locationxi=(/linmsg(locationxi, -1)/)
; locationdong=(/linmsg(locationdong, -1)/)
 location = [/ locationxi, locationdong /]   ; return two variables as a type 'list'
return location
end 

;************************************************
; 画20cr,ncep,era40及jRA-55的msf的变化趋势。
;************************************************
undef("draw_msf_trend")
function draw_msf_trend(LeftString,wks)
begin
    diri="/mnt/f/share/rea_msf/"
    f1=addfile(diri+"msf_time_"+LeftString+"_global.nc", "r")
    date=cd_calendar(f1->time, 1)
    timeind=ind(date.ge.197901.and.date.le.201212)
    levind=ind(f1->level.eq.500)
    
    msf_time_r=f1->msf_time(:,:,timeind)
    models=LeftString
    ; print("*************"+models(i)+"*************")

    msf_time_f=lonPivot((msf_time_r(time|:,level|:,lon|:)),180)
    delete(msf_time_r)
    if models.eq."MERRA" then
    msf_time_f&lon=fspan(-180,177.5,143)
    end if
    msf_time_20cr=msf_time_f(level|:,lon|:,time|:)
    delete(msf_time_f)
    printVarSummary(msf_time_20cr)

    msf_time_20cr=(/msf_time_20cr/10^9/)
    delete(msf_time_20cr&time)
    msf_time_20cr&time=date(timeind)
    delete(f1) 
    delete(date)
    delete(timeind) 
;;;

;;;去趋势
   ; msf_time_20cr_de=dtrend_msg_n(ispan(1,dimsizes(msf_time_20cr&time),1),msf_time_20cr,False,True,2)
   msf_time_20cr_de=msf_time_20cr
   copy_VarCoords(msf_time_20cr, msf_time_20cr_de)

;;;   
  ;  location_ceshi=search_location(msf_time_20cr(time|:,level|levind,lon|:))

  ;  location_xi=location_ceshi[0]
  ; delete(location_xi&time)
  ;   location_xi!0="time"
  ;   location_xi&time=date(timeind) ;
  ;  delete(date)
  ;  delete(timeind)
  ;   loxi_djf=month_to_season(location_xi, "DJF")
  ;   loxi_djf!0="time"
  ;   loxi_djf&time=ispan(1979,2012,1)  ;

  ;   ;;注意点
  ;   loxi_djf(0:dimsizes(loxi_djf&time)-2)=(/loxi_djf(1:)/)
  ;   loxi_djf(dimsizes(loxi_djf&time)-1)=(/location_xi({201212})/)
  ;   ;;;
  ;   loxi_djf=(/dim_rmvmean(loxi_djf)/)

;;;先得到秋季和冬季eninio年的沃克环流情况SONDJF
 msf_son=month_to_season(msf_time_20cr_de(time|:,level|:,lon|:), "SON")
 msf_djf=month_to_season(msf_time_20cr_de(time|:,level|:,lon|:), "DJF")

 msf_son&time=ispan(1979,2012,1)
 msf_djf&time=ispan(1979,2012,1)  ;;;因为这里第一个是1958 1月和2月的平均，第二个是195812,195901,195902的平均
 ;;注意点
    msf_djf(0:dimsizes(msf_djf&time)-2,:,:)=(/msf_djf(1:,:,:)/)
    msf_djf(dimsizes(msf_djf&time)-1,:,:)=(/msf_time_20cr_de(:,:,{201212})/)
 ;;;

 msf_clm=dim_avg_n((/msf_djf(:,:,:)/)*1.0,0)
 copy_VarCoords(msf_son(1,:,:), msf_clm)
 msf_cold=msf_djf
;;;
 msf_cold=(/msf_djf(:,:,:)/) 
 msf_cold_mean=level_mean(msf_cold(:,{600:400},:))

 ; levind=ind(msf_cold&level.eq.500)
 ; msf_cold_mean=msf_cold(time|:,level|levind,lon|:)
 
 location=search_location_id(msf_cold_mean)
 print(location[0])
;;
 
    locationxi=dim_rmvmean(location[0])
    locationdong=dim_rmvmean(location[1])
    width=location[1] - location[0]
    width=(/dim_rmvmean(width)/)
    ;;;读入ENSO海温
    diri="/mnt/f/share/ersst/"
    fils2=systemfunc("ls "+diri+"sst.mnmean.v4.nc")
    ; fils1=systemfunc("ls "+diri+"HadISST_sst.nc")
    f1=addfile(fils2, "r")
    date=cd_calendar(f1->time, 1)
    timeind=ind(date.le.201312.and.date.ge.197901)
    sst=f1->sst(timeind,{20:-20},:)

    delete(fils2)
    ssta=(/rmMonAnnCycTLL(sst)/)
    copy_VarCoords(sst, ssta)
    ssta_de=dtrend_msg_n(ispan(1,dimsizes(ssta&time),1),ssta,False,True,0)
    copy_VarCoords(ssta, ssta_de)
    ssta_djf=month_to_season(ssta_de, "DJF")
     m=""+ispan(1979,2013,1)+""
    ssta_s=(/ssta_djf(1:,:,:)/)*1.0
    copy_VarCoords(ssta_djf(0:dimsizes(m)-2,:,:), ssta_s)

    nino3_4=dim_avg_n(ssta_s(:,{5:-5},{190:240}), (/1,2/))  
   
    ;;;remove enso cycle

     ; a1=regCoef_n(nino3_4, locationxi, 0, 0)
     ; locationxi_enso=a1*nino3_4+a1@yintercept
     ; locationxi=(/locationxi- locationxi_enso/)
    
     ; a2=regCoef_n(nino3_4, locationdong, 0, 0)
     ; locationdong_enso=a2*nino3_4+a2@yintercept
     ; locationdong=(/locationdong- locationdong_enso/)
     
     ; a3=regCoef_n(nino3_4, width, 0, 0)
     ; width_enso=a3*nino3_4+a3@yintercept
     ; width=(/width- width_enso/)
    ;;;
    r=escorc(nino3_4, width)
    print("************"+LeftString+"*************")
    print("r="+r)

     



    ;;;


;;输出每个数据location的位置
; opt=diri+LeftString+"_msf_cold_locationxi_node_index_1979.ascii"    
; asciiwrite(opt, sprintf ("%6.2f,",location[0]))  

; opt1=diri+LeftString+"_msf_cold_locationdong_node_index_1979.ascii"    
; asciiwrite(opt1, sprintf ("%6.2f,",location[1]))  

locationxi=(/runave(locationxi, 3, 1)/)
locationdong=(/runave(locationdong, 3, 1)/)
width=(/runave(width, 3, 1)/)


;;;画图
return [/locationxi ,locationdong,width/]
end

;***********************************************************
;主程序
;***********************************************************         
begin


wks=gsn_open_wks("pdf", "/mnt/f/share/era/data/大西洋all_msf_位置的趋势")


;;;画图
lo1=draw_msf_trend("mre",wks)
lo2=draw_msf_trend("era",wks)
lo3=draw_msf_trend("jra55",wks)
lo4=draw_msf_trend("MERRA",wks)
lo5=draw_msf_trend("ncep1",wks)
lo6=draw_msf_trend("ncep2",wks)
; 
; lo5=draw_msf_trend(fils5,wks,"mre",pa)


 lo7 = NewList("lifo")
;;;MRE还是取平均比较稳
do i=0,2
 ListPush(lo7, (lo1[i]+lo2[i]+lo3[i]+lo4[i]+lo5[i]+lo6[i])/6.0)
end do
time=ispan(1,34,1)

rc1=regCoef_n(time, lo7[2],0,0)
rc_line1=rc1*time+rc1@yintercept

rc2=regCoef_n(time, lo7[1],0,0)
rc_line2=rc2*time+rc2@yintercept

rc3=regCoef_n(time, lo7[0],0,0)
rc_line3=rc3*time+rc3@yintercept

;;;
h1=(/lo1[0],lo2[0],lo3[0],lo4[0],lo5[0],lo6[0],lo7[2] ,rc_line1/)

h2=(/lo1[1],lo2[1],lo3[1],lo4[1],lo5[1],lo6[1] ,lo7[1], rc_line2/)

h3=(/lo1[2],lo2[2],lo3[2],lo4[2],lo5[2],lo6[2] ,lo7[0], rc_line3/)
;;;
res=True
res@gsnDraw=False
res@gsnFrame=False
res@gsnYRefLine=0.0
; res@gsnYRefLineDashPattern=1
res@gsnLeftString="(a)western edge"
res@tmXTOn=False
res@tmYROn=False
res@tmXBMode="Explicit"

res@trXMaxF=34
res@trXMinF=1

; res@trYMaxF=30
; res@trYMinF=-30

res@tmXBValues=ispan(1,34,10)
res@tmXBLabels=ispan(1979,2012,10)
res@tmXBMinorValues=ispan(1,34,1)


res@vpHeightF=0.2
res@vpWidthF=0.7


res@xyLineColors=(/"blue","green","orange","burlywood4","red" ,"blueviolet" ,"black","azure4"/)
res@xyDashPattern=(/0,0,0,0,0,0,0,1/)
res@xyLineThicknesses=(/2,2,2,2,2,2,2,2/)

;;;
;;;画图
plot=new(3,graphic)
plot(0)=gsn_csm_xy(wks, time, h1, res)
res@gsnLeftString="(b)eastern edge"
; res@trYMaxF=12
; res@trYMinF=-12

plot(1)=gsn_csm_xy(wks, time, h2 , res)
res@gsnLeftString="(c)width"
; res@trYMaxF=30
; res@trYMinF=-30
plot(2)=gsn_csm_xy(wks, time, h3 , res)

text=(/"20CR","ERA","JRA-55","MERRA","NCEP1","NCEP2","ENS"/)
x=(/6,10,15,20,24,28,32/)
y=(/27,27,27,27,27,27,27/)
color=(/"blue","green","orange","burlywood4","red" ,"blueviolet" ,"black","azure4"/)
m=new(7,graphic)
textres=True
textres@txFontHeightF =12.0
do i=0,6
textres@txFontColor=color(i)
m(i)=gsn_add_text(wks,plot(0),text(i),x(i),y(i),textres)
end do

pres=True
gsn_panel(wks, plot, (/3,1/), pres)

end






