
;************************************************
; 计算20cr,ncep,era40及jRA-55的msf的变化趋势。
;************************************************
undef("draw_msf_trend")
function draw_msf_trend(fils1,fils2,LeftString,pa)
begin
    diri="/mnt/f/share/"
    locationxi_r=asciiread(diri+fils1, 55, "float")
    locationdong_r=asciiread(diri+fils2, 55, "float")   
    width_r=locationdong_r-locationxi_r
     
    locationxi=dim_rmvmean(locationxi_r)
    locationdong=dim_rmvmean(locationdong_r)
    
    
    ; locationxi=locationxi_r
    ; locationdong=locationdong_r
    
    width=dim_rmvmean(width_r)

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
 
    
    ssta_s=(/ssta_son(0:54,:,:)+ssta_djf(1:,:,:)/)*0.5
    

    copy_VarCoords(ssta_son(0:54,:,:), ssta_s)
    emi=dim_avg_n(ssta_s(:,{10:-10},{165:220}),(/1,2/))-0.5*dim_avg_n(ssta_s(:,{20:-10},{125:145}), (/1,2/))-0.5*dim_avg_n(ssta_s(:,{5:-15},{250:290}),(/1,2/))
    
    emi_center=dim_avg_n(ssta_s(:,{10:-10},{165:220}),(/1,2/))

    emi_std=dim_standardize_n(emi, 1, 0)
    std_emi=dim_stddev_n(emi, 0)

    nino3=dim_avg_n(ssta_s(:,{5:-5},{210:270}), (/1,2/))
    nino3_std=dim_standardize_n(nino3, 1, 0)
    std_nino3=dim_stddev_n(nino3, 0)
      
    nino3_4=dim_avg_n(ssta_s(:,{5:-5},{190:240}), (/1,2/))    
;;;
   
    diri="/mnt/f/share/"
    locationxi_r=asciiread(diri+fils1, 55, "float")
    locationdong_r=asciiread(diri+fils2, 55, "float")   

   

;;;挑选年份 modoki year
 modokiind=ind(emi_std.ge.std_emi)
 loxi_modoki=avg(locationxi(modokiind))
 lodong_modoki=avg(locationdong(modokiind))
 ; print("loxi_modoki_mean="+loxi_modoki_mean)
 ; print("lodong_modoki_mean="+lodong_modoki_mean)
 ; print("with_modoki="+with_modoki)
;;;

;;;挑选年份 canino  enino year
 caind=ind(nino3_std.ge.std_nino3)
 loxi_cano=avg(locationxi(caind))
 lodong_cano=avg(locationdong(caind))
 
;;;
;;;挑选年份 canino  lanina year
 laind=ind(nino3_std.le.(/-1*std_nino3/))
 loxi_lano=avg(locationxi(laind))
 lodong_lano=avg(locationdong(laind))
;;;

;;; normal time
 
 loxi_no=avg(locationxi)
 lodong_no=avg(locationdong)

; ;;;挑选年份 lanina modoki year
;  lamodokiind=ind(emi_std.le.(/-1*std_emi/))
 
;  loxi_no=locationxi(noind)
;  loxi_no_mean=avg(loxi_no)

;  lodong_no=locationdong(noind)
;  lodong_no_mean=avg(lodong_no)

;  with_no=avg(width(noind))
;;

;;; 计算与nino3_4的相关性

   ; r1=escorc(locationxi, nino3_4)
   ; r2=escorc(locationdong, nino3_4)
   ; r3=escorc(width, nino3_4)

   r1=escorc(locationxi, emi_center)
   r2=escorc(locationdong, emi_center)
   r3=escorc(width, emi_center)
   ;;;
   r=[/r1,r2,r3/]
;;;;

if pa.eq."nino3" then 
    nino3=nino3
else 
    nino3=emi
end if
;;;;计算与nino3偏相关性
   r_xi_nino3=escorc(locationxi, nino3)
   r_dong_nino3=escorc(locationdong, nino3)
   r_xi_dong=escorc(locationxi, locationdong)
   r3=escorc(width, nino3)

    r_xi_nino3_pa=(/(r_xi_nino3- r_xi_dong*r_dong_nino3)/sqrt((1- r_xi_dong^2)*(1- r_dong_nino3^2))/)
    r_dong_nino3_pa=(/(r_dong_nino3- r_xi_dong*r_xi_nino3)/sqrt((1- r_xi_dong^2)*(1- r_xi_nino3^2))/)
;;;
    r_pa=[/r_xi_nino3_pa, r_dong_nino3_pa/]



locationxi=(/runave(locationxi, 3, 1)/)
locationdong=(/runave(locationdong, 3, 1)/)
width=(/runave(width, 3, 1)/)



return [/locationxi,locationdong,width,loxi_cano,lodong_cano,loxi_lano,lodong_lano,loxi_modoki,lodong_modoki,loxi_no,lodong_no,r,r_pa/]

end

;***********************************************************
;主程序
;***********************************************************         
begin

; fils1="20CR_msf_cold_locationxi_index.ascii"
; fils2="NCPE1_msf_cold_locationxi_index.ascii"
; fils3="ERA_msf_cold_locationxi_index.ascii"
; fils4="JRA-55_msf_cold_locationxi_index.ascii"
; fils5="MRE_msf_cold_locationxi_index.ascii"
; ; fils5="MRE_msf_mean_cold_locationxi_index.ascii"


; fils6="20CR_msf_cold_locationdong_index.ascii"
; fils7="NCPE1_msf_cold_locationdong_index.ascii"
; fils8="ERA_msf_cold_locationdong_index.ascii"
; fils9="JRA-55_msf_cold_locationdong_index.ascii"
; fils10="MRE_msf_cold_locationdong_index.ascii"
; ; fils10="MRE_msf_mean_cold_locationdong_index2.ascii"

;;;
fils1="20CR_msf_cold_locationxi_node_index.ascii"
fils2="NCPE1_msf_cold_locationxi_node_index.ascii"
fils3="ERA_msf_cold_locationxi_node_index.ascii"
fils4="JRA-55_msf_cold_locationxi_node_index.ascii"
fils5="MRE_msf_cold_locationxi_node_index.ascii"
; fils5="MRE_msf_mean_cold_locationxi_index.ascii"


fils6="20CR_msf_cold_locationdong_node_index.ascii"
fils7="NCPE1_msf_cold_locationdong_node_index.ascii"
fils8="ERA_msf_cold_locationdong_node_index.ascii"
fils9="JRA-55_msf_cold_locationdong_node_index.ascii"
fils10="MRE_msf_cold_locationdong_node_index.ascii"
; fils10="MRE_msf_mean_cold_locationdong_index2.ascii"






pa="emi"
lo1=draw_msf_trend(fils1,fils6,"20CR",pa)
lo2=draw_msf_trend(fils2,fils7,"NCPE1",pa)
lo3=draw_msf_trend(fils3,fils8,"ERA",pa)
lo4=draw_msf_trend(fils4,fils9,"JRA-55",pa)
lo5=draw_msf_trend(fils5,fils10,"MRE",pa)

; ;;;输出位置指数的MRE的结果
;  m=new((/3,dimsizes(lo1[0])/),typeof(lo1[0]))
;  m(0,:)=(lo1[0]+lo2[0]+lo3[0]+lo4[0])/4.0
;  m(1,:)=(lo1[1]+lo2[1]+lo3[1]+lo4[1])/4.0
;  m(2,:)=(lo1[2]+lo2[2]+lo3[2]+lo4[2])/4.0

; opt="/mnt/f/share/MRE_msf_mean_cold_locationxi_index.ascii"    
; asciiwrite(opt, sprintf ("%6.2f,",m(0,:)))  
; opt1="/mnt/f/share/MRE_msf_mean_cold_locationdong_index.ascii"     
; asciiwrite(opt1, sprintf ("%6.2f,",m(1,:)))  
; opt2="/mnt/f/share/MRE_msf_mean_cold_locationwidth_index.ascii"    
; asciiwrite(opt2, sprintf ("%6.2f,",m(2,:)))  

; exit()

wks=gsn_open_wks("pdf", "/mnt/f/share/era/data/4个再分析资料_自定义的西边界东边界位置的年际变化没去趋势")

 ;    ;计算MRE的相关性
 ;    m=new((/3,dimsizes(lo1[0])/),typeof(lo1[0]))
 ;    do i=0,2
 ;    m(i,:)=(lo1[i]+lo2[i]+lo3[i]+lo4[i])/4.0
 ;    r1=escorc(m(i,:), lo1[11])
 ;    print("r1="+r1)
 ;    end do
 ; delete(m)

    ;;计算各自的相关性
    ; do i=1,5
    ;  m=[/lo1[11],lo2[11],lo3[11],lo4[11],lo5[11]/]
    ; print("****"+i+"******")
    ; r=m[i-1]
    ; print("r1="+r[0])
    ; print("r2="+r[1])
    ; print("r3="+r[2])
    ; delete(r)
    ; end do
;    ;;输入各偏相关系数
    ; do i=1,5
    ; m=[/lo1[12],lo2[12],lo3[12],lo4[12],lo5[12]/]
    ; print("****"+i+"******")
    ; r=m[i-1]
    ; print("r1="+r[0])
    ; print("r2="+r[1])
    ; delete(r)
    ; end do
 


do i=3,10
; mean=(lo1[i]+lo2[i]+lo3[i]+lo4[i])*0.25
mean=lo4[i]
; if i.le.9 then 
; h=lo2[i+1]-lo2[i]
; print("h="+h)
; delete(h)
; end if 
print("lo"+i+"="+mean)

end do

 

 k1=lo2[4]-lo2[3]
 k2=lo2[8]-lo2[7]
 
print("k1="+k1)
print("k2="+k2)


year=dimsizes(lo1[0])
time=ispan(1,year,1)
;;;draw three picture
print("jixus")

    res=True
    res@gsnDraw=False
    res@gsnFrame=False
    res@gsnYRefLine=0.0
    ; res@gsnYRefLineDashPattern=1
    res@gsnLeftString="western edge"
    res@tmXTOn=False
    res@tmYROn=False
    res@tmXBMode="Explicit"
    
    res@trXMaxF=55
    res@trXMinF=1
    
    res@trYMaxF=20
    res@trYMinF=-20

    ; res@trYMaxF=6.0
    ; res@trYMinF=-6.0
    ; res@tmXBMinorValues=ispan(1,dimsizes(year),1)
    ; res@tmXBValues=ispan(3,dimsizes(year),5)
    ; res@tmXBLabels=ispan(1960,2010,5)
    
    res@tmXBValues=ispan(1,53,10)
    res@tmXBLabels=ispan(1958,2010,10)
    res@tmXBMinorValues=ispan(1,55,1)


    res@vpHeightF=0.2
    res@vpWidthF=0.7
     

    res@xyLineColors=(/"blue","green","orange","red","black","gray"/)
    res@xyDashPattern=(/0,0,0,0,0,1/)
    res@xyLineThicknesses=(/2,2,2,2,2,2/)
  
    

    ;;;MRE还是取平均比较稳
    do i=0,2
    lo5[i]=(lo1[i]+lo2[i]+lo3[i]+lo4[i])*0.25
   end do

   rc1=regCoef_n(time, lo5[0],0,0)
   rc_line1=rc1*time+rc1@yintercept

   rc2=regCoef_n(time, lo5[1],0,0)
   rc_line2=rc2*time+rc2@yintercept
   
   rc3=regCoef_n(time, lo5[2],0,0)
   rc_line3=rc3*time+rc3@yintercept
   
    ;;;
    h1=(/lo1[0],lo2[0],lo3[0],lo4[0],lo5[0],rc_line1/)
   
    h2=(/lo1[1],lo2[1],lo3[1],lo4[1],lo5[1],rc_line2/)
 
    h3=(/lo1[2],lo2[2],lo3[2],lo4[2],lo5[2],rc_line3/)
;;;



;;;
;;;画图
    plot=new(3,graphic)
    plot(0)=gsn_csm_xy(wks, time, h1, res)
     res@gsnLeftString="eastern edge"
    res@trYMaxF=12
    res@trYMinF=-12

    plot(1)=gsn_csm_xy(wks, time, h2 , res)
       res@gsnLeftString="width"
       res@trYMaxF=20
       res@trYMinF=-20
    plot(2)=gsn_csm_xy(wks, time, h3 , res)
    
   text=(/"20CR","NCEP1","ERA","JRA-55","MRE"/)
   x=(/10,20,30,40,50/)
   y=(/17,17,17,17,17/)
   color=(/"blue","green","orange","red","black"/)
   m=new(5,graphic)
   textres=True
   textres@txFontHeightF =15.0
   do i=0,4
   textres@txFontColor=color(i)
    m(i)=gsn_add_text(wks,plot(0),text(i),x(i),y(i),textres)
   end do
    
   
    pres=True
    gsn_panel(wks, plot, (/3,1/), pres)

end



