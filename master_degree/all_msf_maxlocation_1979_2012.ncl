
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

msf_mean_annual_xi=msf_mean_annual(:,{120:250})
;;;judge locationxi
do t=0,dimsizes(msf_mean_annual&time)-1
  locationxi(t)=msf_mean_annual_xi&lon(maxind(msf_mean_annual_xi(t,:)))
end do
 location = [/ locationxi /]   ; return two variables as a type 'list'
return location
end 
;;; 
;;; 
undef("search_location_alt")
function  search_location_alt(msf_mean_annual)
local msf_mean_annual_xi , msf_mean_annual_dong
begin
locationxi=new(dimsizes(msf_mean_annual&time),typeof(msf_mean_annual&lon))
locationxi!0="time"
locationxi&time=ispan(0,dimsizes(msf_mean_annual&time)-1,1)

msf_mean_annual_xi=msf_mean_annual(:,{300:350})
;;;judge locationxi
do t=0,dimsizes(msf_mean_annual&time)-1
  locationxi(t)=msf_mean_annual_xi&lon(maxind(msf_mean_annual_xi(t,:)))
end do
 location = [/ locationxi /]   ; return two variables as a type 'list'
return location
end 
;;; 
undef("search_location_id")
function  search_location_id(msf_mean_annual)
local msf_mean_annual_xi , msf_mean_annual_dong
begin
locationxi=new(dimsizes(msf_mean_annual&time),typeof(msf_mean_annual&lon))
locationxi!0="time"
locationxi&time=ispan(0,dimsizes(msf_mean_annual&time)-1,1)

msf_mean_annual_xi=msf_mean_annual(:,{40:150})
;;;judge locationxi
do t=0,dimsizes(msf_mean_annual&time)-1
  locationxi(t)=msf_mean_annual_xi&lon(minind(msf_mean_annual_xi(t,:)))
end do
 location = [/ locationxi /]   ; return two variables as a type 'list'
return location
end 
;;; 
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

    msf_time_20cr=f1->msf_time(:,:,timeind)
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

    ;;;先得到秋季和冬季eninio年的沃克环流情况SONDJF
    msf_son=month_to_season(msf_time_20cr_de(time|:,level|:,lon|:), "SON")
    msf_djf=month_to_season(msf_time_20cr_de(time|:,level|:,lon|:), "DJF")

    msf_son&time=ispan(1979,2012,1)
    msf_djf&time=ispan(1979,2012,1)  ;;;因为这里第一个是1958 1月和2月的平均，第二个是195812,195901,195902的平均
    ; ;;注意点
    ; msf_djf(0:dimsizes(msf_djf&time)-2,:,:)=(/msf_djf(1:,:,:)/)
    ; msf_djf(dimsizes(msf_djf&time)-1,:,:)=(/msf_time_20cr_de(:,:,{201212})/)
    ; ;;;
    msf_cold=msf_djf
    ;;;
    msf_cold=(/msf_djf(:,:,:)/) 
    ; msf_cold_mean=level_mean(msf_cold(:,:,:))
    msf_cold_mean=msf_cold(:,{500},:)
 

    location=search_location(msf_cold_mean)
    location_id=search_location_id(msf_cold_mean)
 
;;
    width_pa= location[0]
    width_pa=(/dim_rmvmean(width_pa)/)
 
    width_id= location_id[0]
    width_id=(/dim_rmvmean(width_id)/)
    
    printVarSummary(msf_cold_mean)
    msf_cold_mean_pi=lonPivot(msf_cold_mean, 30)
    printVarSummary(msf_cold_mean_pi)

    location_al=search_location_alt(msf_cold_mean_pi)

    width_alt=location_al[0]
    width_alt=(/dim_rmvmean(width_alt)/)
    

  

;;;读入海温数据
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
;;;处理成异常的冬季平均场
  ssta=(/rmMonAnnCycTLL(sst(:,:,:))/)
  copy_VarCoords(sst(:,:,:),ssta)
  ssta_de=ssta
  ssta_djf=month_to_season(ssta_de, "DJF")
  ssta_s=(/ssta_djf(1:,:,:)/)*1.0
  copy_VarCoords(ssta_djf(0:dimsizes(ssta_djf&time)-2,:,:), ssta_s)
  nino3_4=dim_avg_n(ssta_s(:,{5:-5},{190:240}), (/1,2/))  

r_id_nino3_4=escorc_n(width_id, nino3_4, 0, 0)
r_pa_nino3_4=escorc_n(width_pa, nino3_4, 0, 0)

prob_id_nino3_4=rtest(r_id_nino3_4, 34, 0)
prob_pa_nino3_4=rtest(r_pa_nino3_4, 34, 0)
print("r_id_nino3_4="+r_id_nino3_4)
print("prob_id_nino3_4="+prob_id_nino3_4)
print("r_pa_nino3_4"+r_pa_nino3_4)
print("prob_pa_nino3_4"+prob_pa_nino3_4)
;;;;

print("******"+"width"+"****")

r_width_pa_id=escorc_n(width_id, width_pa, 0, 0)
print("r_width_pa_id="+r_width_pa_id)



    width_pa=(/runave(width_pa, 3, 0)/)
    width_id=(/runave(width_id, 3, 0)/)
    width_alt=(/runave(width_alt, 3, 0)/)
;;;画图
return [/width_id ,width_pa,width_alt/]
end

;***********************************************************
;主程序
;***********************************************************         
begin


wks=gsn_open_wks("pdf", "/mnt/f/share/era/data/all_msf_3个环流的最大量级位置变化")


;;;画图
lo1=draw_msf_trend("20cr",wks)
lo2=draw_msf_trend("era",wks)
lo3=draw_msf_trend("jra55",wks)
lo4=draw_msf_trend("MERRA",wks)
lo5=draw_msf_trend("ncep1",wks)
lo6=draw_msf_trend("ncep2",wks)
; 
; lo5=draw_msf_trend(fils5,wks,"mre",pa)
width_id_ncep2=lo6[0]
width_pa_ncep2=lo6[1]
width_id_jra55=lo3[0]
width_pa_jra55=lo3[1]


width_id_ncep2(32)=width_id_jra55(32)
width_pa_ncep2(32)=width_pa_jra55(32)
lo6[0]=width_id_ncep2
lo6[1]=width_pa_ncep2


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
res@gsnLeftString="(a)IWCW"
res@tmXTOn=False
res@tmYROn=False
res@tmXBMode="Explicit"

res@trXMaxF=34
res@trXMinF=1

res@trYMaxF=30
res@trYMinF=-30

res@tmXBValues=ispan(1,34,10)
res@tmXBLabels=ispan(1979,2012,10)
res@tmXBMinorValues=ispan(1,34,1)


res@vpHeightF=0.2
res@vpWidthF=0.7


; res@xyLineColors=(/"blue","green","orange","burlywood4","red" ,"blueviolet" ,"black","azure4"/)
; res@xyDashPattern=(/0,0,0,0,0,0,0,1/)
; res@xyLineThicknesses=(/2,2,2,2,2,2,2,2/)
res@xyLineColors=(/"blue","green","orange","burlywood4","red" ,"blueviolet" ,"black"/)
res@xyDashPattern=(/0,0,0,0,0,0,0/)
res@xyLineThicknesses=(/2,2,2,2,2,2,2/)

;;;
;;;画图
plot=new(3,graphic)
plot(0)=gsn_csm_xy(wks, time, h1(0:6,:), res)
res@gsnLeftString="(b)PWCW"
res@trYMaxF=30
res@trYMinF=-30

plot(1)=gsn_csm_xy(wks, time, h2(0:6,:) , res)
res@gsnLeftString="(c)AWCW"
res@trYMaxF=16
res@trYMinF=-16
plot(2)=gsn_csm_xy(wks, time, h3(0:6,:) , res)

text=(/"20CR","ERA","JRA-55","MERRA","NCEP1","NCEP2","ENS"/)
x=(/6,10,15,20,24,28,32/)
y=(/27,27,27,27,27,27,27/)
; color=(/"blue","green","orange","burlywood4","red" ,"blueviolet" ,"black","azure4"/)
color=(/"blue","green","orange","burlywood4","red" ,"blueviolet" ,"black"/)
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






