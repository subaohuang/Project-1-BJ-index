
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

;;;
begin
;************************************************
; 读入数据
;************************************************
diri="/mnt/f/share/ncep/20CR/"
diri_pi="/mnt/f/share/era/data/"
fils1=systemfunc("ls "+diri+"msf_time_20cr_global.nc")
f1=addfile(fils1, "r")
; msf_time=f1->msf_time(:,:,348:731)
; msf_time=f1->msf_time(:,:,0:347)
msf_time=f1->msf_time(:,:,0:731)
plev=msf_time&level
psfc=conform_dims((/dimsizes(msf_time&time),dimsizes(msf_time&lon)/), 1013.0, 0)
ptop=min(plev)
dp   = dpres_plevel_Wrap(plev, psfc, ptop, 0)  ; dp(30)
printVarSummary(dp)
msfdp=msf_time(level|:,time|:,lon|:)*dp
copy_VarCoords(msf_time(level|:,time|:,lon|:), msfdp)
printVarSummary(msfdp)
msf_mean=dim_sum_n_Wrap(msfdp, 0)/dim_sum_n_Wrap(dp, 0)
copy_VarCoords(msf_time(time|:,lon|:,level|1), msf_mean)
printVarSummary(msf_mean)
;;;

msf_mean_iod=dim_avg_n_Wrap(msf_mean(:,{40:110}), 1)
msf_mean_pac=dim_avg_n_Wrap(msf_mean(:,{150:260}), 1)


; msf_haiyang=dim_avg_n_Wrap(msf_mean, dims)

opt=diri_pi+"pa_msf_index_1950_2010.ascii"    ;;;输出强度指数
asciiwrite(opt, sprintf ("%6.2f,",msf_mean_pac))  ;;;输出强度指数

opt1=diri_pi+"id_msf_index_1950_2010.ascii"
asciiwrite(opt1, sprintf ("%6.2f,",msf_mean_iod))

msf_mean_pac_anu=month_to_annual(msf_mean_pac, 1)
msf_mean_iod_anu=month_to_annual(msf_mean_iod, 1)
msf_mean_pac_anu&year=ispan(1, dimsizes(msf_mean_pac_anu), 1)
msf_mean_iod_anu&year=ispan(1, dimsizes(msf_mean_iod_anu), 1)
  
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

;;;
printVarSummary(msf_mean_pac_anu)
rc3=regCoef_n(ispan(1,30,1), msf_mean_pac_anu(0:29), 0, 0)
msf_pac_line1=rc3*ispan(1,30,1) + rc3@yintercept

rc4=regCoef_n(ispan(1,31,1), msf_mean_pac_anu(30:), 0, 0)
msf_pac_line2=rc4*ispan(1,31,1) + rc4@yintercept


rc5=regCoef_n(ispan(1,29,1), msf_mean_iod_anu(0:28), 0, 0)
msf_iod_line1=rc5*ispan(1,30,1) + rc5@yintercept

rc6=regCoef_n(ispan(1,32,1), msf_mean_iod_anu(29:), 0, 0)
msf_iod_line2=rc6*ispan(1,32,1) + rc6@yintercept
print("r6="+rc6)
;;;

;************************************************
;  计算滑动相关情况
;************************************************

    nyear=15
    corr=slid_cor(nyear,(/-1*msf_mean_iod_anu/),msf_mean_pac_anu,1)
    printVarSummary(corr)
    t=2.101
    r_test=sqrt(t^2/(nyear-2+t^2))                     ;0.05显著性检验的t值为2.101

    r_test_new=conform_dims(dimsizes(corr), r_test,0)
;;;
   printVarSummary(msf_mean_pac_anu)

    h1= escorc(msf_mean_pac_anu(0:29),msf_mean_iod_anu(0:29))
    h2= escorc(msf_mean_pac_anu(30:60),msf_mean_iod_anu(30:60))  
    h3= escorc(msf_mean_pac_anu(29:49),msf_mean_iod_anu(29:49))


    print("h1="+h1)
    print("h2="+h2)
    print("h3="+h3)
  exit()



;;;画图

    wks=gsn_open_wks("pdf", ""+diri_pi+"msf定义的pa指数")
    wks1=gsn_open_wks("pdf", ""+diri_pi+"msf定义的id指数")
    wks2=gsn_open_wks("pdf", ""+diri_pi+"msf定义指数滑动相关")

    res=True
    res@gsnDraw=False
    res@gsnFrame=False
    ; plot=gsn_csm_xy(wks, msf_mean_pac&time, (/msf_mean_pac,msf_pac_line/), res)
    ; plot1=gsn_csm_xy(wks1, msf_mean_iod&time, (/msf_mean_iod,msf_iod_line/), res)

    plot=gsn_csm_xy(wks, msf_mean_pac_anu&year, (/msf_mean_pac_anu,msf_pac_line/), res)
    plot1=gsn_csm_xy(wks1, msf_mean_iod_anu&year, (/msf_mean_iod_anu,msf_iod_line/), res)

    resline=True
    resline@gsLineColor="black"
    resline@gsLineThicknessF=2.0
    res@gsLineColor="green"

    Polyline1=gsn_add_polyline(wks, plot, ispan(1, 30, 1), msf_pac_line1, resline)
    res@gsLineColor="green"
    Polyline2=gsn_add_polyline(wks, plot, ispan(31, 61, 1), msf_pac_line2, resline)

    Polyline3=gsn_add_polyline(wks, plot1, ispan(1, 29, 1), msf_iod_line1, resline)
    res@gsLineColor="green"
    Polyline4=gsn_add_polyline(wks, plot1, ispan(30, 61, 1), msf_iod_line2, resline)



    draw(plot)
    draw(plot1)
    frame(wks)
    frame(wks1)
    
    res2=True
    res2@tmXBMode="Explicit"
    res2@tmXBMinorValues=ispan(0,dimsizes(corr)-1,1)
    res2@tmXBValues=ispan(2,dimsizes(corr)-1,10)
    res2@tmXBLabels=fspan(1960,2000,5)
    res2@xyDashPatterns=(/0,1/)
    plot2=gsn_csm_xy(wks2, ispan(0,dimsizes(corr)-1,1), (/corr,r_test_new/), res2)
    

  
    draw(plot2)
    frame(wks2)



end








;;;