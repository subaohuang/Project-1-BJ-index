
function level_mean(data)
local data,msf_mean
begin
plev=data&level
psfc=conform_dims((/dimsizes(data&time),dimsizes(data&lon)/), 1013.0, 0)
ptop=min(plev)
dp   = dpres_plevel_Wrap(plev, psfc, ptop, 0)  ; dp(30)
msfdp=data(level|:,time|:,lon|:)*dp
copy_VarCoords(data(level|:,time|:,lon|:), msfdp)

msf_mean=dim_sum_n_Wrap(msfdp, 0)/dim_sum_n_Wrap(dp, 0)
copy_VarCoords(data(time|:,lon|:,level|1), msf_mean)
return(msf_mean)
end 

undef("search_location")
function  search_location(msf_mean_annual)
local msf_mean_annual_xi , msf_mean_annual_dong
begin
locationxi=new(dimsizes(msf_mean_annual&time),typeof(msf_mean_annual&lon))
locationxi!0="time"
locationxi&time=ispan(0,dimsizes(msf_mean_annual&time)-1,1)
locationdong=locationxi
msf_mean_annual_xi=msf_mean_annual(:,{120:250})
msf_mean_annual_dong=msf_mean_annual(:,{230:300})
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
   if ismissing(locationdong(t)).eq.True  then 
   co_in=closest_val(0,msf_mean_annual_dong(t,:))
   locationdong(t)=msf_mean_annual_dong&lon(co_in)
   delete(co_in)
end if

end do
locationxi=(/linmsg(locationxi, -1)/)
; locationdong=(/linmsg(locationdong, -1)/)
location = [/ locationxi, locationdong /]   ; return two variables as a type 'list'
return location
end 

function select_win(in_20cr)
begin
  in_20cr_win=new(3,typeof(in_20cr))
  in_20cr_win(0:1)=in_20cr(9:10)     ;;;;1,2,12月
  in_20cr_win(2)=in_20cr(11)

return(in_20cr_win)
end  
;********************************************************
;draw_msf_hov
;*******************************************************
function draw_msf_hov(wks,LeftString1,LeftString2)
begin 
print("************************"+LeftString1+"**************************")
yearstart=197901
year=1979
    diri="/mnt/f/share/rea_msf/"
    fils1=diri+"msf_time_"+LeftString1+"_global.nc"
    f1=addfile(fils1, "r")
    date=cd_calendar(f1->time, 1)
    timeind=ind(date.ge.yearstart.and.date.le.201212)

    ; levind=ind(f1->level.eq.500)
    level=f1->level
    levind1=ind(level.eq.500)
    levind2=ind(level.le.600.and.level.ge.400)
    msf_time=f1->msf_time(:,:,timeind)
    delete(f1)
    delete(fils1)
    msf_time=(/msf_time/10^9/)
 
    delete(msf_time&time)
    msf_time&time=date(timeind)
    msf_time=(/dtrend_msg_n(ispan(1,dimsizes(msf_time&time),1),msf_time,False,True,2)/)
   
    msf_annual=month_to_annual(msf_time(time|:,level|:,lon|:), 1)
    msf_annual!0="time"
    msf_annual&time=ispan(year,2012,1)

    msf_djf=month_to_season(msf_time(time|:,level|:,lon|:), "DJF")
    msf_djf(0:dimsizes(msf_djf&time)-2,:,:)=(/msf_djf(1:,:,:)/)
    msf_djf(dimsizes(msf_djf&time)-1,:,:)=(/msf_time(:,:,{201212})/)
    msf_djf&time=ispan(year,2012,1)  ;;;因为这里第一个是1958 1月和2月的平均，第二个是195812,195901,195902的平均
    ; location=search_location(msf_djf(:,{500},:))
    ; msf_djf_mean=level_mean(msf_djf(time|:,level|levind,lon|:))
    ; msf_time_mean=level_mean(msf_djf(time|:,level|levind,lon|:))
    msf_mean=level_mean(msf_djf(time|:,level|:,lon|:))
    msf_mean2=level_mean(msf_djf(time|:,level|levind2,lon|:))
    ; location=search_location(msf_djf(time|:,level|levind1,lon|:))
    location=search_location(msf_mean2)
    locationxi=location[0]
    locationdong=location[1]
    width=locationdong -locationxi

    delete(locationxi&time)
    locationxi&time=ispan(year,2012,1)
    
    delete(locationdong&time)
    locationdong&time=ispan(year,2012,1)
    copy_VarCoords(locationxi, width)
     
    delete([/date,timeind/])

    ; locationxi_annual=month_to_annual(locationxi, 1)
    ; locationdong_annual=month_to_annual(locationdong, 1)
    ; width_annual=month_to_annual(width, 1)

    ; locationxi_djf=month_to_season(locationxi, "DJF")
    ; locationxi_djf(0:dimsizes(locationxi_djf&time)-2)=(/locationxi_djf(1:)/)
    ; locationxi_djf(dimsizes(locationxi_djf&time)-1)=(/locationxi({201212})/)
    ; locationxi_djf&time=ispan(year,2012,1)  ;;;因为这里第一个是1958 1月和2月的平均，第二个是195812,195901,195902的平均
    
    ; locationdong_djf=month_to_season(locationdong, "DJF")
    ; locationdong_djf(0:dimsizes(locationdong_djf&time)-2)=(/locationdong_djf(1:)/)
    ; locationdong_djf(dimsizes(locationdong_djf&time)-1)=(/locationdong({201212})/)
    ; locationdong_djf&time=ispan(year,2012,1)  ;;;因为这里第一个是1958 1月和2月的平均，第二个是195812,195901,195902的平均

    ; width_djf=month_to_season(width, "DJF")
    ; width_djf(0:dimsizes(width_djf&time)-2)=(/width_djf(1:)/)
    ; width_djf(dimsizes(width_djf&time)-1)=(/width({201212})/)
    ; width_djf&time=ispan(year,2012,1)  ;;;因为这里第一个是1958 1月和2月的平均，第二个是195812,195901,195902的平均
    ; printVarSummary(locationxi_djf)
 
    indxi=ind_nearest_coord(locationxi, msf_djf&lon, 0)
    inddong=ind_nearest_coord(locationdong, msf_djf&lon, 0)
    
   indmax=inddong
  do i=0,dimsizes(msf_djf&time)-1
  indmax(i)=maxind(msf_djf(i,{500},:))
  end do 
 

    msf_inten=new(dimsizes(msf_djf&time),typeof(msf_djf))



do i=0,dimsizes(msf_djf&time)-1
msf_inten(i)=(/dim_avg_n_Wrap(msf_mean(i,{150:210}),0)/)
; msf_inten(i)=(/dim_avg_n_Wrap(msf_annual(i,{500},{150:210}),0)/)
end do 

    ; rc1=regCoef(locationxi_annual,msf_inten)
    ; line = rc1*locationxi_annual+rc1@yintercept
     rc1=regline_stats(width,msf_inten)
     line = rc1*width+rc1@yintercept
     R=rc1
     
     r=escorc(msf_inten, width)

    ; Yest_025 = rc1@YPI025 
    ; Yest_975 = rc1@YPI975
     ; Yest_025 = rc1@y95(0) + rc1@b95(0)*width  
    ; Yest_975 = rc1@y95(1) + rc1@b95(1)*width  

    ;;;画ploygon
    ; mono = 1                            ; ascending=1 , descending=-1
    ; ii   =  dim_pqsort_n(width,mono,0)
    ; xx   = width(ii)                        ; ascending order 
    ; yy1   = Yest_025(ii)
    ; yy2   =Yest_975(ii)

    ; n=dimsizes(width)
    ; xp    = new( (/2*dimsizes(width)/), typeof(width))
    ; yp    = new( (/2*dimsizes(width)/), typeof(Yest_975) )
    ; do k=0,dimsizes(width)-1
    ; xp(k) = xx(k)
    ; yp(k)= yy2(k)
    ; xp(2*n-1-k) = xx(k)
    ; yp(2*n-1-k) = yy1(k)
    ; end do
   ;;;corr

  ;read sst data
    diri="/mnt/f/share/ersst/"
    fils1=systemfunc("ls "+diri+"sst.mnmean.v4.nc")
     ; fils1=systemfunc("ls "+diri+"HadISST_sst.nc")
    f1=addfile(fils1, "r")
    date=cd_calendar(f1->time, 1)
    timeind=ind(date.le.201312.and.date.ge.yearstart)
   
    ; print("*****************ERSST*********************")
    ; ;;;读入ERSST海温数据
    ; m=""+ispan(1979,2013,1)+""
    ; diri="/mnt/f/share/ersst/"

    ; fils=new((/35,12/),string)
    ; do i=0,34
    ; fils1=systemfunc("ls "+diri+"ersst."+m(i)+"*"+".nc")
    ; fils(i,:)= fils1
    ; end do
    ; fils_new=onedtond(fils, 35*12)
    ; fileh=new(35*12,string)
    ; fileh(0:35*12-1)=fils_new

    ; f1=addfiles(fileh, "join")
    ; printVarSummary(f1)
    ; date=cd_calendar(f1[:]->time, 1)
    ;  timeind=ind(date.le.201312.and.date.ge.yearstart)
    ; sst=f1[:]->sst(timeind,0,{30:-30},:)
    sst=f1->sst(timeind,{30:-30},:)
    ; sst=lonFlip(f1->sst(timeind,{30:-30},:))

    sst!0="time"
    sst!1="lat"
    sst!2="lon"
    delete(date)
    delete(f1)
    

    ssta=(/rmMonAnnCycTLL(sst)/)
    copy_VarCoords(sst, ssta)
    ssta_de=dtrend_msg_n(ispan(1,dimsizes(ssta&time),1),ssta,False,True,0)
    printVarSummary(ssta_de)
    copy_VarCoords(ssta, ssta_de)
    ssta_djf=month_to_season(ssta_de, "DJF")
    ssta_son=month_to_season(ssta_de, "SON")
    
    ssta_s=(/ssta_djf(1:,:,:)/)*1.0
   
    ; ssta_s=(/month_to_annual(ssta_de(0:407,:,:), 1)/)

    ssta_q=(/ssta_son(0:dimsizes(ssta_djf&time)-2,:,:)/)*1.0
    ; ssta_s=(/ssta_q/)

    copy_VarCoords(ssta_djf(0:dimsizes(ssta_djf&time)-2,:,:), ssta_s)
    copy_VarCoords(ssta_djf(0:dimsizes(ssta_djf&time)-2,:,:), ssta_q)
  
    emi=dim_avg_n(ssta_s(:,{10:-10},{165:220}),(/1,2/))-0.5*dim_avg_n(ssta_s(:,{20:-10},{125:145}), (/1,2/))-0.5*dim_avg_n(ssta_s(:,{5:-15},{250:290}),(/1,2/))
    ; emi=dim_avg_n(ssta_q(:,{10:-10},{165:220}),(/1,2/))-0.5*dim_avg_n(ssta_q(:,{20:-10},{125:145}), (/1,2/))-0.5*dim_avg_n(ssta_q(:,{5:-15},{250:290}),(/1,2/))
    emi_std=dim_standardize_n(emi, 1, 0)
 
    nino3=dim_avg_n(ssta_s(:,{5:-5},{210:270}), (/1,2/))
    nino3_std=dim_standardize_n(nino3, 1, 0)
    
    ; modokiind=ind(emi_std.ge.0.7)
    ; modokilaind=ind(emi_std.le.-1)
    ; caind=ind(nino3_std.ge.1)
    ; laind=ind(nino3_std.le.(/-1/))
    ; noind=ind(nino3_std.ge.(/-0.5/).and.nino3_std.le.(/0.5/))
    
    cano_year=(/1982,1991,1997/)
    modoki_year=(/1994,2002,2004,2006,2009/)
    lanina_year=(/1983,1988,1995,1998,2005,2007,2010/)
    noind_year=(/1979,1980,1981,1984,1985,1989,1990,1992,1993,1996,1999,2000,2001,2003,2008,2011,2012/)

    caind=get1Dindex(width&time,cano_year)
    modokiind=get1Dindex(width&time,modoki_year)
    laind=get1Dindex(width&time,lanina_year)
    noind=get1Dindex(width&time,noind_year)

    ; can_year=(/1958,1965,1972,1976,1982,1986,1997/)
    ; can_year=(/1982,1986,1997/)
    ; ; ; modo_year=(/1963,1968,1977,1987,1990,1991,1992,1993,1994,2002,2004,2009/)
    ; modo_year=(/1987,1990,1991,1992,1994,2002,2004,2009/)
    ; caind=get1Dindex(width&time,can_year)
    ; modokiind=get1Dindex(width&time,modo_year)

    width_djf_cano=width(caind)
    width_djf_lani=width(laind)
    width_djf_modo=width(modokiind)
    width_djf_nor=width(noind)

      
    locationdong_djf_cano=locationdong(caind)
    locationdong_djf_lani=locationdong(laind)
    locationdong_djf_modo=locationdong(modokiind)
    locationdong_djf_nor=locationdong(noind)


      
    locationxi_djf_cano=locationxi(caind)
    locationxi_djf_lani=locationxi(laind)
    locationxi_djf_modo=locationxi(modokiind)
    locationxi_djf_nor=locationxi(noind)



    msf_inten_djf_cano=msf_inten(caind)
    msf_inten_djf_lani=msf_inten(laind)
    msf_inten_djf_modo=msf_inten(modokiind)
    msf_inten_djf_nor=msf_inten(noind)


print( "canoind="+width&time(caind))
print("modokiind="+width&time(modokiind))
print( "laind="+width&time(laind))

; print("********"+LeftString+"_各事件的平均"+"*************")

; print("locationxi_cano="+dim_avg_n_Wrap(locationxi_djf_cano, 0))
; print("locationdong_cano="+dim_avg_n_Wrap(locationdong_djf_cano, 0))
; print("locationwidth_cano="+dim_avg_n_Wrap(width_djf_cano, 0))

; print("locationxi_lani="+dim_avg_n_Wrap(locationxi_djf_lani, 0))
; print("locationdong_lani="+dim_avg_n_Wrap(locationdong_djf_lani, 0))
; print("locationwidth_lani="+dim_avg_n_Wrap(width_djf_lani, 0))

; print("locationxi_modo="+dim_avg_n_Wrap(locationxi_djf_modo, 0))
; print("locationdong_modo="+dim_avg_n_Wrap(locationdong_djf_modo, 0))
; print("locationwidth_modo="+dim_avg_n_Wrap(width_djf_modo, 0))

;;;相关性
print("********"+LeftString1+"_各事件的相关性"+"*************")
 r=escorc(width, msf_inten)



 r1=escorc(locationxi, nino3)
 r2=escorc(locationdong, nino3)
 r3=escorc(width, nino3)

 r4=escorc(locationxi, emi)
 r5=escorc(locationdong, emi)
 r6=escorc(width, emi)

    print("r="+sprintf("%4.2f", r))
    print("r1="+sprintf("%4.2f", r1))
    print("r2="+sprintf("%4.2f", r2))
    print("r3="+sprintf("%4.2f", r3))
    print("r4="+sprintf("%4.2f", r4))
    print("r5="+sprintf("%4.2f", r5))
    print("r6="+sprintf("%4.2f", r6))
  print("************偏相关*************")
    r7=(/(r1- r*r2)/sqrt((1- r^2)*(1- r2^2))/)
    r8=(/(r2- r*r1)/sqrt((1- r^2)*(1- r1^2))/)
    
     print("r7="+sprintf("%4.2f", r7))
    print("r8="+sprintf("%4.2f", r8))
    r9=(/(r4- r*r5)/sqrt((1- r^2)*(1- r5^2))/)
    r10=(/(r5- r*r4)/sqrt((1- r^2)*(1- r4^2))/)
    print("r9="+sprintf("%4.2f", r9))
    print("r10="+sprintf("%4.2f", r10))

 
  siglvl = 0.05            ; a-priori specified sig level
  pr     = rtest((/r1,r2,r3,r4,r5,r6/), dimsizes(locationxi), 0)    ; ===> pr=0.046
  print("pr="+pr)
print("*********************************************")

;;;画散点图


    res                   = True                     ; plot mods desired
    res@gsnDraw=False
    res@gsnFrame=False
    res@xyLineThicknesses = (/ 1.0,2.0/)          ; make second line thicker
    res@xyLineColors      = (/"black"  ,"black"/)          ; change line color
    res@xyMarkerColor = "black"
    res@xyMarkLineModes   = (/"Markers","lines"/)
    res@xyDashPattern=0
    ; res@xyExplicitLegendLabels = (/"Anual change", "mon_change", "Linear trend"/) 
    res@xyMarkers         =  16                      ; choose type of marker  
    res@xyMarkerSizeF     = 0.005        ; Marker size (default 0.01)
    res@gsnLeftString=LeftString2
    res@gsnRightString="r="+sprintf("%4.2f", r)
    res@vpHeightF=0.3
    res@vpWidthF=0.7
    ; res@tfPolyDrawOrder    = "Predraw"               ; put line on top 
    ; res@trYMinF     = -3.0        ; set minimum Y-axis value
    ; res@trYMaxF     =  6.0        ; set maximum Y-axis value
    ; res@trXMinF     = 80        ; set minimum X-axis value
    ; res@trXMaxF     = 140        ; set maximum X-axis value

     plot  = gsn_csm_xy (wks,width,(/msf_inten,line/),res) ; create plot
 
    polyres=True
    polyres                   = True
    polyres@gsMarkerIndex     = 16          ; polymarker style
    polyres@gsMarkerSizeF     = 10        ; polymarker size
    polyres@gsMarkerColor     ="pink"
    ; polyres@tfPolyDrawOrder   = "Postdraw" 
  
    str1 = unique_string("dum1")  ; "unique_string" will return a unique                                     ; string every time it is called fro                                 ;  within a single NCL session.
    str2 = unique_string("dum2")
    str3 = unique_string("dum3")
    str4 = unique_string("dum4")
    
  

    plot@$str1$=gsn_add_polymarker(wks, plot, width_djf_modo, msf_inten_djf_modo, polyres)
    polyres@gsMarkerColor     ="blue"
    plot@$str2$=gsn_add_polymarker(wks, plot, width_djf_lani, msf_inten_djf_lani, polyres)
    polyres@gsMarkerColor     ="red"
    plot@$str3$=gsn_add_polymarker(wks, plot, width_djf_cano, msf_inten_djf_cano, polyres)

   


locationxi_cano=dim_avg_n_Wrap(locationxi_djf_cano, 0)
locationdong_cano=dim_avg_n_Wrap(locationdong_djf_cano, 0)
width_cano=dim_avg_n_Wrap(width_djf_cano, 0)

locationxi_lani=dim_avg_n_Wrap(locationxi_djf_lani, 0)
locationdong_lani=dim_avg_n_Wrap(locationdong_djf_lani, 0)
width_lani=dim_avg_n_Wrap(width_djf_lani, 0)

locationxi_modo=dim_avg_n_Wrap(locationxi_djf_modo, 0)
locationdong_modo=dim_avg_n_Wrap(locationdong_djf_modo, 0)
width_modo=dim_avg_n_Wrap(width_djf_modo, 0)


locationxi_nor=dim_avg_n_Wrap(locationxi_djf_nor, 0)
locationdong_nor=dim_avg_n_Wrap(locationdong_djf_nor, 0)
width_nor=dim_avg_n_Wrap(width_djf_nor, 0)

lo=(/locationxi_cano,locationxi_lani,locationxi_modo,locationxi_nor, locationdong_cano,locationdong_lani,locationdong_modo,locationdong_nor\
,width_cano,width_lani,width_modo,width_nor/)


return [/plot,lo/]
end
;;;*************************主程序**************************
;;;
;;;*********************************************************


begin
rea=(/"20cr","era","jra55","MERRA","ncep1","ncep2"/)
; rea=(/"20cr","in_era40","jra55","ncep1"/)
LeftString2=(/"(a)20CR","(b)ERAIM","(c)JRA55","(d)MERRA","(e)NCEP1","(d)NCEP2"/)
plot=new(dimsizes(rea),graphic)

wks=gsn_open_wks("pdf", "/mnt/f/share/era/data/500hpa_6个再分析资料djf强度与宽度的散点图")

m=NewList("lifo")
do i=0,dimsizes(rea)-1
ListPush(m,draw_msf_hov(wks,rea(i),LeftString2(i)))
end do
printVarSummary(m)

plot_20cr=m[5]
plot_era=m[4]
plot_jra55=m[3]
plot_merra=m[2]
plot_ncep1=m[1]
plot_ncep2=m[0]
printVarSummary(plot_era)

plot(0)=plot_20cr[0]
plot(1)=plot_era[0]
plot(2)=plot_jra55[0]
plot(3)=plot_merra[0]
plot(4)=plot_ncep1[0]
plot(5)=plot_ncep2[0]


lo_20cr=plot_20cr[1]
lo_era=plot_era[1]
lo_jra55=plot_jra55[1]
lo_merra=plot_merra[1]
lo_ncep1=plot_ncep1[1]
lo_ncep2=plot_ncep2[1]

lo_mean=new(12,typeof(lo_20cr))
lo=new((/12,6/),typeof(lo_20cr))
printVarSummary(lo_20cr)
do i=0,11
  lo_mean(i)=tofloat((/lo_20cr(i)+lo_era(i)+lo_jra55(i)+lo_merra(i)+lo_ncep1(i)+lo_ncep2(i)/)/6)
  lo(i,0)=tofloat(lo_20cr(i))
  lo(i,1)=tofloat(lo_era(i))
  lo(i,2)=tofloat(lo_jra55(i))
  lo(i,3)=tofloat(lo_merra(i))
  lo(i,4)=tofloat(lo_ncep1(i))
  lo(i,5)=tofloat(lo_ncep2(i))
end do

print("lo_mean="+lo_mean)

std=dim_stddev_n(lo, 1)
print("std="+std)


  xleg = (/0.16,0.36,0.6,0.8/)   ; Location of

  xtxt = (/0.23,0.43,0.67,0.85/)   ; legend markers

  yleg = (/0.1,0.1,0.1,0.1/)   ; and text

  ytxt = (/0.1,0.1,0.1,0.1/)   ; strings.
 
  colors=(/"red","blue","pink","black"/)
  labels=(/" EP El Nino"," La Nina","CP El Nino ","Neutral"/)

  txres               = True
  txres@txFontHeightF = 0.015
  
   gsres               = True
  gsres@gsMarkerIndex = 16          ; Use filled dots for markers.

    do i = 0, 3
      gsres@gsMarkerColor      = colors(i)
      ; gsres@gsMarkerThicknessF = 0.7
      if i.eq.3 then
         gsres@gsMarkerSizeF=0.005
     else
     gsres@gsMarkerSizeF=10.0
 end if
      gsn_polymarker_ndc(wks, xleg(i),yleg(i),gsres)
      gsn_text_ndc      (wks,labels(i),xtxt(i),ytxt(i),txres)
    end do


pres=True
pres@gsnMaximize=True
pres@gsnPanelDebug=False
gsn_panel(wks, plot, (/3,2/), pres)
end 