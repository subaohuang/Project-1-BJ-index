
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

undef("search_location_2d_level")
function  search_location_2d_level(msf_mean_annual)
local msf_mean_annual_xi , msf_mean_annual_dong
begin
locationxi=new(dimsizes(msf_mean_annual&level),typeof(msf_mean_annual&lon))
locationxi!0="level"
locationxi&level=msf_mean_annual&level
locationdong=locationxi
msf_mean_annual_xi=msf_mean_annual(:,{120:250})
msf_mean_annual_dong=msf_mean_annual(:,{230:280})
;;;judge locationxi
do t=0,dimsizes(msf_mean_annual&level)-1
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
do t=0,dimsizes(msf_mean_annual&level)-1
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
   if ismissing(locationdong(t)).eq.True  then 
   co_in=closest_val(0,msf_mean_annual_dong(t,:))      ;;插值的问题比较明显
   locationdong(t)=msf_mean_annual_dong&lon(co_in)
   delete(co_in)
   end if
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
function draw_msf_trend(fils1,wks,LeftString1,LeftString2)
begin
    diri="/mnt/f/share/"
    print("fils1="+fils1)
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
    
    ; ; ;读入ERSST海温数据
    m=""+ispan(1979,2013,1)+""
    ; diri="/mnt/f/share/ersst/"
    ; fils=new((/dimsizes(m),12/),string)
    ; do i=0,dimsizes(m)-1
    ; fils2=systemfunc("ls "+diri+"ersst."+m(i)+"*"+".nc")
    ; fils(i,:)= fils2
    ; end do

    ; fils_new=onedtond(fils, dimsizes(m)*12)
    ; fileh=new(dimsizes(m)*12,string)
    ; fileh(0:dimsizes(m)*12-1)=fils_new
    ; f2=addfiles(fileh, "join")
    ; printVarSummary(f2)

    ; date=cd_calendar(f2[:]->time, 1)
    ; timeind=ind(date.le.201312.and.date.ge.197901)
    ; sst=f2->sst(timeind,{20:-20},:)
    ; sst=f2[:]->sst(timeind,0,{20:-20},:)
    ;  f2=addfile("/mnt/f/share/ersst/HadISST_sst.nc", "r")
    ; date=cd_calendar(f2->time, 1)
    ; timeind=ind(date.ge.197901.and.date.le.201312)
    ; sst=lonFlip(f2->sst(timeind,{20:-20},:))
  
    ;;;


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
    
    ssta_s=(/ssta_djf(1:,:,:)/)*1.0
    copy_VarCoords(ssta_djf(0:dimsizes(m)-2,:,:), ssta_s)

    emi=dim_avg_n(ssta_s(:,{10:-10},{165:220}),(/1,2/))-0.5*dim_avg_n(ssta_s(:,{20:-10},{125:145}), (/1,2/))-0.5*dim_avg_n(ssta_s(:,{5:-15},{250:290}),(/1,2/))
   
    emi_std=dim_standardize_n(emi, 1, 0)
    std_emi=dim_stddev_n(emi, 0)
   
    nino3=dim_avg_n(ssta_s(:,{5:-5},{210:270}), (/1,2/))
    nino3_std=dim_standardize_n(nino3, 1, 0)
    std_nino3=dim_stddev_n(nino3, 0)
    

    nino4=dim_avg_n(ssta_s(:,{5:-5},{160:210}), (/1,2/))
    nino4_std=dim_standardize_n(nino4, 1, 0)
    std_nino4=dim_stddev_n(nino4, 0)

    delete(diri)
    diri="/mnt/f/share/rea_msf/"

;;;

  
;;;去趋势
   msf_time_20cr_de=dtrend_msg_n(ispan(1,dimsizes(msf_time_20cr&time),1),msf_time_20cr,False,True,2)
   ; msf_time_20cr_de=msf_time_20cr
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


;;输出每个数据location的位置
opt=diri+LeftString1+"_msf_djf_locationxi_node_index.ascii"    
asciiwrite(opt, sprintf ("%6.2f,",location[0]))  

opt1=diri+LeftString1+"_msf_djf_locationdong_node_index.ascii"    
asciiwrite(opt1, sprintf ("%6.2f,",location[1]))  



cano_year=(/1982,1991,1997/)
modoki_year=(/1994,2002,2004,2006,2009/)
lanina_year=(/1983,1988,1995,1998,2005,2007,2010/)
noind_year=(/1979,1980,1981,1984,1985,1989,1990,1992,1993,1996,1999,2000,2001,2003,2008,2011,2012/)

caind=get1Dindex(msf_djf&time,cano_year)
modokiind=get1Dindex(msf_djf&time,modoki_year)
laind=get1Dindex(msf_djf&time,lanina_year)
noind=get1Dindex(msf_djf&time,noind_year)


;;;挑选年份 modoki year
 ; modokiind=ind(emi_std.ge.0.7)
 msf_modoki=dim_avg_n((/msf_djf(modokiind,:,:)+msf_djf(modokiind,:,:)/)*0.5,0)
 ; msf_modoki=dim_avg_n((/msf_djf(modokiind+1,:,:)/),0)   ;;;;djf取的要加1
 copy_VarCoords(msf_djf(1,:,:), msf_modoki)
;;;


;;;挑选年份 canino  enino year
 ; caind=ind(nino3_std.ge.1)
 msf_canonical=dim_avg_n((/msf_djf(caind,:,:)+msf_djf(caind,:,:)/)*0.5,0)
 ; msf_canonical=dim_avg_n((/msf_djf(caind+1,:,:)/),0)
 copy_VarCoords(msf_djf(1,:,:), msf_canonical)
;;;

;;;挑选年份 canino  lanina year
 ; laind=ind(nino3_std.le.(/-1/))
 msf_lanina=dim_avg_n((/msf_djf(laind,:,:)+msf_djf(laind,:,:)/)*0.5,0)
 ; msf_canonical=dim_avg_n((/msf_djf(caind+1,:,:)/),0)
 copy_VarCoords(msf_djf(1,:,:), msf_lanina)
;;; normal time
 ; noind=ind(nino3_std.ge.(/-0.5/).and.nino3_std.le.(/0.5/))
 msf_noind=dim_avg_n((/msf_djf(noind,:,:)+msf_djf(noind,:,:)/)*0.5,0)
 ; msf_canonical=dim_avg_n((/msf_djf(caind+1,:,:)/),0)
 copy_VarCoords(msf_djf(1,:,:), msf_noind)

;;;挑选年份 lanina modoki year
 lamodokiind=ind(emi_std.le.(/-1/))
 msf_lamodoki=dim_avg_n((/msf_djf(lamodokiind,:,:)+msf_djf(lamodokiind,:,:)/)*0.5,0)
 ; msf_canonical=dim_avg_n((/msf_djf(caind+1,:,:)/),0)
 copy_VarCoords(msf_djf(1,:,:), msf_lamodoki)
;;
 nino_year=(/1982,1991,1994,1997,2002,2004,2006,2009/)
 caind1=ind(nino3_std.gt.nino4_std.and.nino3_std.ge.1)
 ep_time=msf_djf&time(caind1)
 modoind1=ind(nino3_std.lt.nino4_std)
 cp_time=msf_djf&time(modoind1)
 print(ep_time)
 print(cp_time)


;;;计算响应的结果，即用响应的nino指数做标准化
   nino3_can=avg(nino3(caind))
   nino3_lan=avg(nino3(laind))
   emi_modoki=avg(emi(modokiind))
;;;

lo_clm=search_location_2d(msf_clm({500},:))
width_clm=lo_clm[1]-lo_clm[0]



if  fils1.eq."rea_msf/msf_time_mre_global.nc"then 
print("******cano***********")
 lo=search_location_2d(msf_canonical({500},:))
 rsp1=(lo[0]-lo_clm[0])/nino3_can
 rsp2=(lo[1]-lo_clm[1])/nino3_can
 rsp3=(lo[1]-lo[0]-width_clm)/nino3_can
 width=lo[1]-lo[0]
 print("lo_xi="+lo[0])
 print("lo_dong="+lo[1])
 print("width="+width)
 print("rsp1="+rsp1)
 print("rsp2="+rsp2)
 print("rsp3="+rsp3)
 delete(lo)
 delete(rsp1)
 delete(rsp2)
 delete(rsp3)
 delete(width)
print("******lanina***********")
lo=search_location_2d(msf_lanina({500},:))
 width=lo[1]-lo[0]
 rsp1=(lo[0]-lo_clm[0])/nino3_lan
 rsp2=(lo[1]-lo_clm[1])/nino3_lan
 rsp3=(lo[1]-lo[0]-width_clm)/nino3_lan
 print("lo_xi="+lo[0])
 print("lo_dong="+lo[1])
  print("width="+width)
 print("rsp1="+rsp1)
 print("rsp2="+rsp2)
print("rsp3="+rsp3)
delete(lo)
delete(rsp1)
delete(rsp2)
delete(rsp3)
delete(width)
print("******modoki***********")
lo=search_location_2d(msf_modoki({500},:))
 width=lo[1]-lo[0]
 rsp1=(lo[0]-lo_clm[0])/emi_modoki
 rsp2=(lo[1]-lo_clm[1])/emi_modoki
 rsp3=(lo[1]-lo[0]-width_clm)/emi_modoki
  print("width="+width)
print("rsp1="+rsp1)
print("rsp2="+rsp2)
print("rsp3="+rsp3)
print("lo_xi="+lo[0])
print("lo_dong="+lo[1])

delete(lo)
delete(rsp1)
delete(rsp2)
delete(rsp3)
delete(width)
 print("******normal***********")
lo=search_location_2d(msf_noind({500},:))
 width=lo[1]-lo[0]
 print("width="+width)
 print("lo_xi="+lo[0])
 print("lo_dong="+lo[1])
 delete(lo)
end if


;;;画图

    res2 = True
    res2@gsnDraw=False
    res2@gsnFrame=False
    res2@gsnLeftString=LeftString2
    res2@cnLinesOn=True
    res2@cnLevelSelectionMode ="ExplicitLevels"
    res2@cnLevels=0.0
    res2@cnFillOn=False
    res2@cnLineColor="blue"
    res2@cnMonoLineColor=True
    res2@cnLineThicknessF =4.0
    
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
     
    res2@vpHeightF=0.35
    res2@vpWidthF=0.65
   
    plot=new(1,graphic)
    
    plot=gsn_csm_pres_hgt(wks,msf_lanina(:,{120:285}),res2)

    res2@cnLineColor="burlywood4"
    plot1=gsn_csm_pres_hgt(wks,msf_modoki(:,{120:285}),res2)
    ; res2@cnLineColor="black"
    ; plot2=gsn_csm_pres_hgt(wks,msf_clm,res2)    
    res2@cnLineColor="red"
    if fils1.eq."rea_msf/msf_time_jra55_global.nc"
    plot2=gsn_csm_pres_hgt(wks,msf_canonical(:,{120:260}),res2)  
    else
    plot2=gsn_csm_pres_hgt(wks,msf_canonical(:,{120:285}),res2)  
    end if
    res2@cnLineColor="green"
    plot3=gsn_csm_pres_hgt(wks,msf_noind(:,{120:285}),res2)  
    if fils1.eq."rea_msf/msf_time_ncep2_global.nc"
    overlay(plot,plot1)
    overlay(plot,plot3)   
  else 
    overlay(plot,plot1)
    overlay(plot,plot2)
    overlay(plot,plot3)   
  end if
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
    
    str4 = unique_string("dum4")
    str5 = unique_string("dum5")
    ;;;添加东边界位置
    if fils1.eq."rea_msf/msf_time_ncep2_global.nc"
 
    ployres@gsLineThicknessF = 4.0            ; line thickness   
    ployres@gsLineColor      = "red"         ; color of lines
    location_can=search_location_2d_level(msf_canonical(:,:))
    locationdong_can=location_can[1]
    locationxi_can=location_can[0]


    plot@$str4$= gsn_add_polyline(wks,plot,locationdong_can,locationdong_can&level, ployres)
    plot@$str5$= gsn_add_polyline(wks,plot,locationxi_can,locationxi_can&level, ployres)
    end if
    
      if fils1.eq."rea_msf/msf_time_jra55_global.nc"
    ployres@gsLineThicknessF = 4.0            ; line thickness   
    ployres@gsLineColor      = "red"         ; color of lines
    location_can=search_location_2d_level(msf_canonical(:,:))
    locationdong_can=location_can[1]
    locationxi_can=location_can[0]
    plot@$str4$= gsn_add_polyline(wks,plot,locationdong_can,locationdong_can&level, ployres)
    end if

return plot
end

;***********************************************************
;主程序
;***********************************************************         
begin




fils1="rea_msf/msf_time_20cr_global.nc"

; fils3="era/era40/msf_time_era40_in_global.nc"
fils2="rea_msf/msf_time_era_global.nc"
fils3="rea_msf/msf_time_jra55_global.nc"
fils4="rea_msf/msf_time_MERRA_global.nc"
fils5="rea_msf/msf_time_ncep1_global.nc"
fils6="rea_msf/msf_time_ncep2_global.nc"
fils7="rea_msf/msf_time_mre_global.nc"

wks=gsn_open_wks("pdf", "/mnt/f/share/era/data/6个再分析资料_自定义的各类事件位置冬")

plot=new(7,graphic)
plot(0)=draw_msf_trend(fils1,wks,"20CR","(a)20CR")
plot(1)=draw_msf_trend(fils2,wks,"ERA","(b)ERAIM")
plot(2)=draw_msf_trend(fils3,wks,"JRA55","(c)JRA55")
plot(3)=draw_msf_trend(fils4,wks,"MERRA","(d)MERRA")
plot(4)=draw_msf_trend(fils5,wks,"NCEP1","(e)NCEP1")
plot(5)=draw_msf_trend(fils6,wks,"NCEP2","(f)NCEP2")
plot(6)=draw_msf_trend(fils7,wks,"MRE","ENS")

  

  
  txres               = True
  txres@txFontHeightF = 0.02
  

  xtxt = (/0.2,0.45,0.7,0.9/)   ; legend markers
  ytxt = (/0.065,0.065,0.065,0.065/)   ; strings.
 

  xtxt = (/0.25,0.45,0.65,0.8/)   ; legend markers
  ytxt = (/0.38,0.38,0.38,0.38/)   ; strings.
    txres@txFontHeightF = 0.015
 
   labels=(/"EP El Nino","La Nina","CP El Nino","Neutral"/)
   color=(/"red","blue","burlywood4","green"/)
  do i=0,3
  txres@txFontColor = color(i)  
  gsn_text_ndc(wks,labels(i),xtxt(i),ytxt(i),txres)
  end do
    ; pres=True
    ;  pres@gsnPanelLabelBar=False
    ;  pres@gsnPanelDebug=True
    ; gsn_panel(wks, plot(0:5), (/3,2/), pres)

draw(plot(6))
frame(wks)

end



