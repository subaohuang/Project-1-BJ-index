
function slid_cor(m,x,y,dimnum);滑动相关,一维的情况
begin  
if dimnum.eq.1 then
length = dimsizes(x)
R=new(length-m+1,typeof(x))
do t=m,length
  R(t-m)=escorc(x(t-m:t-1),y(t-m:t-1))   
end do
end if

if dimnum.eq.2 then 
length = dimsizes(x(0,:))
R=new((/length-m+1,dimsizes(x(0,:))/),typeof(x))
do t=m,length
  R(t-m,:)=escorc_n(x(t-m:t-1,:),y(t-m:t-1,:),0,0)   
end do
end if
return(R)
end

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


;;;

undef("draw_msf_ip_pa")
function draw_msf_ip_pa(fils1,wks,LeftString)
begin
;************************************************
; 读入数据
;************************************************
diri="/mnt/f/share/"

f1=addfile(fils1, "r")
date=cd_calendar(f1->time, 1)
; timeind=ind(date.ge.195801.and.date.le.201212)
timeind=ind(date.ge.197901.and.date.le.201212)


msf_time=f1->msf_time(:,:,timeind)
msf_time=(/msf_time/10^9/)
msf_mean=level_mean(msf_time({600:400},:,:))
;;;

; ;;;选择位置，用位置定义一个精确的msf
; location=search_location(msf_mean)
; locationxi=location[0]
; locationdong=location[1]

; indxi=ind_nearest_coord(locationxi, msf_mean&lon, 0)
; inddong=ind_nearest_coord(locationdong, msf_mean&lon, 0)
; ;;;

msf_mean_iod=dim_avg_n_Wrap(msf_mean(:,{50:120}), 1)
msf_mean_aoc=dim_avg_n_Wrap(msf_mean(:,{290:350}), 1)   ;;70w-10w   
; msf_mean_pac=dim_avg_n_Wrap(msf_mean(:,{180:250}), 1)
msf_mean_pac=dim_avg_n_Wrap(msf_mean(:,{150:210}), 1)

; msf_mean_pac=msf_mean_iod

; do i=0,dimsizes(msf_mean&time)-1
; msf_mean_pac(i)=(/dim_avg_n_Wrap(msf_mean(i,indxi(i):inddong(i)), 0)/)
; end do 

; msf_haiyang=dim_avg_n_Wrap(msf_mean, dims)

opt=diri+LeftString+"_pa_msf_index_1958_2012.ascii"    ;;;输出强度指数
asciiwrite(opt, sprintf ("%6.2f,",msf_mean_pac))  ;;;输出强度指数

opt1=diri+LeftString+"_id_msf_index_1958_2012.ascii"
asciiwrite(opt1, sprintf ("%6.2f,",msf_mean_iod))

opt2=diri+LeftString+"_id_msf_index_1958_2012.ascii"
asciiwrite(opt2, sprintf ("%6.2f,",msf_mean_aoc))

msf_mean_pac_anu=month_to_annual(msf_mean_pac, 1)
msf_mean_iod_anu=month_to_annual(msf_mean_iod, 1)
msf_mean_aoc_anu=month_to_annual(msf_mean_aoc, 1)
msf_mean_pac_anu&year=ispan(1, dimsizes(msf_mean_pac_anu), 1)
msf_mean_iod_anu&year=ispan(1, dimsizes(msf_mean_iod_anu), 1)
msf_mean_aoc_anu&year=ispan(1, dimsizes(msf_mean_aoc_anu), 1)
  
;;;计算趋势

; rc1=regCoef_n(msf_mean_pac&time, msf_mean_pac, 0, 0)

; msf_pac_line=rc1*msf_mean_pac&time + rc1@yintercept

; rc2=regCoef_n(msf_mean_iod&time, msf_mean_iod, 0, 0)

; msf_iod_line=rc2*msf_mean_iod&time + rc2@yintercept
;;;


;;;
rc1=regCoef_n(msf_mean_pac_anu&year, msf_mean_pac_anu, 0, 0)

msf_pac_line=rc1*msf_mean_pac_anu&year + rc1@yintercept

rc2=regCoef_n(msf_mean_iod_anu&year, msf_mean_iod_anu, 0, 0)

msf_iod_line=rc2*msf_mean_iod_anu&year + rc2@yintercept

rc3=regCoef_n(msf_mean_aoc_anu&year, msf_mean_aoc_anu, 0, 0)

msf_aoc_line=rc3*msf_mean_aoc_anu&year + rc3@yintercept


; ;;;太平洋分段趋势
; printVarSummary(msf_mean_pac_anu)
; rc4=regCoef_n(ispan(1,21,1), msf_mean_pac_anu(0:20), 0, 0)
; msf_pac_line1=rc4*ispan(1,21,1) + rc4@yintercept

; rc5=regCoef_n(ispan(22,55,1), msf_mean_pac_anu(21:54), 0, 0)
; msf_pac_line2=rc5*ispan(22,55,1) + rc5@yintercept
; ;;;indian
; rc6=regCoef_n(ispan(1,21,1), msf_mean_iod_anu(0:20), 0, 0)
; msf_iod_line1=rc6*ispan(1,21,1) + rc6@yintercept

; rc7=regCoef_n(ispan(22,55,1), msf_mean_iod_anu(21:54), 0, 0)
; msf_iod_line2=rc7*ispan(22,55,1) + rc7@yintercept
; ;;;altantic 
; rc8=regCoef_n(ispan(1,21,1), msf_mean_aoc_anu(0:20), 0, 0)
; msf_aoc_line1=rc8*ispan(1,21,1) + rc8@yintercept

; rc9=regCoef_n(ispan(22,55,1), msf_mean_aoc_anu(21:54), 0, 0)
; msf_aoc_line2=rc9*ispan(22,55,1) + rc9@yintercept
; ;;;

; ;************************************************
; ;  计算滑动相关情况
; ;************************************************

;     nyear=15
;     corr=slid_cor(nyear,(/-1*msf_mean_iod_anu/),msf_mean_pac_anu,1)
;     printVarSummary(corr)
;     t=2.101
;     r_test=sqrt(t^2/(nyear-2+t^2))                     ;0.05显著性检验的t值为2.101

;     r_test_new=conform_dims(dimsizes(corr), r_test,0)
; ;;;
;    printVarSummary(msf_mean_pac_anu)
;********************************************************************
;draw
;******************************************************************** 
    res=True
    res@gsnDraw=False
    res@gsnFrame=False
    res@gsnLeftString=LeftString
    res@gsnYRefLine=0.0
    res@gsnYRefLineDashPattern=1
    res@tmXTOn=False
    res@tmYROn=False
    res@tmXBMode="Explicit"
    
    ; res@trXMaxF=55
    ; res@trXMinF=1

    res@trYMaxF=6.0
    res@trYMinF=-6.0
    res@tmXBMinorValues=ispan(1,dimsizes(msf_mean_pac_anu&year),1)
    res@tmXBValues=ispan(3,dimsizes(msf_mean_pac_anu&year),10)

    ; res@tmXBLabels=ispan(1960,2010,10)
    res@tmXBLabels=ispan(1980,2010,5)
     

    res@vpHeightF=0.4
    res@vpWidthF=0.6
     

    res@xyLineColors=(/"blue","red","green","blue","red","green"/)
    res@xyDashPattern=(/0,0,0,1,1,1/)
    res@xyLineThicknesses=(/4,4,4,1,1,1/)
  
 
    
    h=(/msf_mean_pac_anu,msf_mean_iod_anu,msf_mean_aoc_anu,msf_pac_line,msf_iod_line,msf_aoc_line/)

    plot=gsn_csm_xy(wks, msf_mean_pac_anu&year, h , res)

    ; resline=True
    ; resline@gsLineThicknessF=2.0
    ; resline@gsLineColor="black"
    ; resline@gsLineDashPattern=1

    ; Polyline1=gsn_add_polyline(wks, plot, ispan(1, 21, 1), msf_pac_line1, resline)
    ; res@gsLineColor="green"
    ; Polyline2=gsn_add_polyline(wks, plot, ispan(22, 55, 1), msf_pac_line2, resline)

    ; Polyline3=gsn_add_polyline(wks, plot, ispan(1, 21, 1), msf_iod_line1, resline)
    ; res@gsLineColor="green"
    ; Polyline4=gsn_add_polyline(wks, plot, ispan(22, 55, 1), msf_iod_line2, resline)

    ; Polyline5=gsn_add_polyline(wks, plot, ispan(1, 21, 1), msf_aoc_line1, resline)
    ; res@gsLineColor="green"
    ; Polyline6=gsn_add_polyline(wks, plot, ispan(22, 55, 1), msf_aoc_line2, resline)

    
return(plot)
end 
;**********************************************************
; main program
;******************
begin 
fils1="ncep/20CR/msf_time_20cr_global_6--6.nc"
fils2="ncep/ncep1/msf_time_ncep1_global.nc"
fils3="era/era40/msf_time_era40_in_global.nc"
fils4="JRA-55/msf_time_jra55_global.nc"

wks=gsn_open_wks("pdf", "/mnt/f/share/era/data/all_msf定义的pa和id指数_1979-2012")

plot=new(4,graphic)
plot(0)=draw_msf_ip_pa(fils1,wks,"20CR")
plot(1)=draw_msf_ip_pa(fils2,wks,"NCPE1")
plot(2)=draw_msf_ip_pa(fils3,wks,"ERA")
plot(3)=draw_msf_ip_pa(fils4,wks,"JRA-55")

pres=True
gsn_panel(wks, plot, (/2,2/), pres)

end








;;;