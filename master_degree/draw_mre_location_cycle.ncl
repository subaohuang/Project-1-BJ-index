
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
msf_mean_annual_dong=msf_mean_annual(:,{230:290})
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
   co_in=closest_val(max(msf_mean_annual(t,:))/8,msf_mean_annual_dong(t,:))
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

function select_season(locationxi)
begin
location_seanson=new(4,typeof(locationxi))
locationxi_chun=avg(locationxi(2:4))
locationxi_xia=avg(locationxi(5:7))
locationxi_qiu=avg(locationxi(8:10))
locationxi_dong=(locationxi(11)+locationxi(1)+locationxi(0))/3.0

location_seanson(0)=locationxi_chun
location_seanson(1)=locationxi_xia
location_seanson(2)=locationxi_qiu
location_seanson(3)=locationxi_dong
return(location_seanson)
end  




;********************************************************
;draw_msf_hov
;*******************************************************
function draw_msf_hov(wks,LeftString)
begin 
print("************************"+LeftString+"**************************")
    diri="/mnt/f/share/rea_msf/"
    fils1=diri+"msf_time_"+LeftString+"_global.nc"
    f1=addfile(fils1, "r")
    date=cd_calendar(f1->time, 1)
    timeind=ind(date.ge.197901.and.date.le.201212)

    level=f1->level
    levind=ind(level.eq.500)
    msf_time=f1->msf_time(:,:,timeind)
    msf_time=(/msf_time/10^9/)
     
    ; printVarSummary(msf_time)
    ; msf_time_month!2="time"
    msf_mean_time=level_mean(msf_time({600:400},:,:))    ;;msf_mean_month(time,lon)
    
    ; printVarSummary(msf_mean_month)
    ; location_r=search_location(msf_time(time|:,level|levind,lon|:))
    ; locationxi_r=location_r[0]
    ; locationdong_r=location_r[1]
    ; locationxi=new(12,typeof(location_r[0]))
    ; locationdong=new(12,typeof(location_r[0]))

     msf_mean_month=new((/12,dimsizes(msf_mean_time&lon)/),typeof(msf_mean_time))
     do i=0,11
     msf_mean_month(i,:)=dim_avg_n_Wrap(msf_mean_time(i:dimsizes(msf_mean_time&time)-1:12,:), 0)
     ; locationxi(i)=dim_avg_n(locationxi_r(i:dimsizes(locationxi_r)-1:12),0)
     ; locationdong(i)=dim_avg_n(locationdong_r(i:dimsizes(locationxi_r)-1:12),0)
 end do  
 ; location=[/locationxi,locationdong/]
 msf_mean_month!0="time"
 msf_mean_month&time=ispan(1, 12, 1)
printVarSummary(msf_mean_month)

    location_r=search_location(msf_mean_month(time|:,lon|:))
    locationxi=location_r[0]
    locationdong=location_r[1]
    
    width = locationdong - locationxi
indxi=ind_nearest_coord(locationxi, msf_mean_month&lon, 0)
inddong=ind_nearest_coord(locationdong, msf_mean_month&lon, 0)
msf_inten_month=new(12,typeof(msf_mean_month))
do i=0,11
msf_inten_month(i)=(/dim_avg_n_Wrap(msf_mean_month(i,{180:240}), 0)/)
end do 


    lo_in=[/locationxi,locationdong,msf_inten_month,width/]
 
 return(lo_in)
 end 

 ;***************************************************************
 ;main program
 ;****************************************************************
  begin 
   rea=(/"20cr","era","jra55","MERRA","ncep1","ncep2"/)
    wks=gsn_open_wks("pdf", "/mnt/f/share/era/data/再分析资料强度与宽度的年循环_测试")
    
     wks1=gsn_open_wks("pdf", "/mnt/f/share/era/data/再分析资料宽度的月变化")
     wks2=gsn_open_wks("pdf", "/mnt/f/share/era/data/再分析资料强度的月变化")
   
     my_list = NewList("lifo")

    do i=0,5
       ListPush(my_list,draw_msf_hov(wks,rea(i)))
    end do
  
   cr20=my_list[5]
   era=my_list[4]
   JRA55=my_list[3]
   merra=my_list[2]
   ncep1=my_list[1]
   ncep2=my_list[0]

   xi_20cr=cr20[0]
   xi_era=era[0]
   xi_jra55=JRA55[0]
   xi_merra=merra[0]
   xi_ncep1=ncep1[0]
   xi_ncep2=ncep2[0]

  dong_20cr=cr20[1]
  dong_era=era[1]
  dong_jra55=JRA55[1]
  dong_merra=merra[1]
  dong_ncep1=ncep1[1]
  dong_ncep2=ncep2[1]

   in_20cr=cr20[2]
   in_era=era[2]
   in_jra55=JRA55[2]
   in_merra=merra[2]
   in_ncep1=ncep1[2]
   in_ncep2=ncep2[2]
   
    width_20cr=cr20[3]
    width_era=era[3]
    width_jra55=JRA55[3]
    width_merra=merra[3]
    width_ncep1=ncep1[3]
    width_ncep2=ncep2[3]



  in_20cr_win=select_win(in_20cr)
  in_era_win=select_win(in_era)
  in_jra55_win=select_win(in_jra55)
  in_merra_win=select_win(in_merra)
  in_ncep1_win=select_win(in_ncep1)
  in_ncep2_win=select_win(in_ncep2)

  xi_20cr_win=select_win(xi_20cr)
  xi_era_win=select_win(xi_era)
  xi_jra55_win=select_win(xi_jra55)
  xi_merra_win=select_win(xi_merra)
  xi_ncep1_win=select_win(xi_ncep1)
  xi_ncep2_win=select_win(xi_ncep2)
  
  dong_20cr_win=select_win(dong_20cr)
  dong_era_win=select_win(dong_era)
  dong_jra55_win=select_win(dong_jra55)
  dong_merra_win=select_win(dong_merra)
  dong_ncep1_win=select_win(dong_ncep1)
  dong_ncep2_win=select_win(dong_ncep2)


  
  in_20cr_season=select_season(in_20cr)
  in_era_season=select_season(in_era)
  in_jra55_season=select_season(in_jra55)
  in_merra_season=select_season(in_merra)
  in_ncep1_season=select_season(in_ncep1)
  in_ncep2_season=select_season(in_ncep2)

  xi_20cr_season=select_season(xi_20cr)
  xi_era_season=select_season(xi_era)
  xi_jra55_season=select_season(xi_jra55)
  xi_merra_season=select_season(xi_merra)
  xi_ncep1_season=select_season(xi_ncep1)
  xi_ncep2_season=select_season(xi_ncep2)
  
  dong_20cr_season=select_season(dong_20cr)
  dong_era_season=select_season(dong_era)
  dong_jra55_season=select_season(dong_jra55)
  dong_merra_season=select_season(dong_merra)
  dong_ncep1_season=select_season(dong_ncep1)
  dong_ncep2_season=select_season(dong_ncep2)

 
  width_20cr_season=select_season(width_20cr)
  width_era_season=select_season(width_era)
  width_jra55_season=select_season(width_jra55)
  width_merra_season=select_season(width_merra)
  width_ncep1_season=select_season(width_ncep1)
  width_ncep2_season=select_season(width_ncep2)


 
 xi_mre=(/xi_20cr_win ,xi_era_win,xi_jra55_win,xi_merra_win,xi_ncep1_win,xi_ncep2_win/)


 dong_mre=(/dong_20cr_win ,dong_era_win,dong_jra55_win,dong_merra_win,dong_ncep1_win,dong_ncep2_win/)
 
 width_mre=(/width_20cr_season ,width_era_season,width_jra55_season,width_merra_season,width_ncep1_season,width_ncep2_season/)
 in_mre=(/in_20cr_season ,in_era_season,in_jra55_season,in_merra_season,in_ncep1_season,in_ncep2_season/)
  ; width_mre=(/width_20cr ,width_era,width_jra55,width_merra,width_ncep1,width_ncep2/)
  ; in_mre=(/in_20cr,in_era,in_jra55,in_merra,in_ncep1,in_ncep2/)
   ; data=(/"20cr","era","jra55","merra",""/)
   r1=escorc(width_20cr_season,in_20cr_season)
   r2=escorc(width_era_season,in_era_season)
   r3=escorc(width_jra55_season,in_jra55_season)
   r4=escorc(width_merra_season,in_merra_season)
   r5=escorc(width_ncep1_season,in_ncep1_season)
   r6=escorc(width_ncep2_season,in_ncep2_season)
   r=(/r1,r2,r3,r5,r6/)
   r_mean=avg(r)
   print("r_rea="+r)
   print("r="+r_mean)
   ; exit()






    res                      = True 

    res@gsnMaximize          = True    ; Maximize plot in window.

    res@gsnDraw              = False   ; Don't draw plot
    res@gsnFrame             = False   ; Don't advance frame

    ;  res@gsnPaperOrientation  = "portrait"

    res@gsnLeftString        = ""
    res@gsnRightString       = ""

    ; res@xyMarkLineModes       = (/"MarkLines","MarkLines","MarkLines","MarkLines","MarkLines","MarkLines"/) 
    ; res@xyMarkers            = (/16,16,16,16,16,16/)         ; filled dot
    ; res@xyMarkerSizeF        = 0.0135
    ; res@xyMarkerThicknessF   = 0.6
    ; res@xyMarkerColors       =(/"black","blue","orange","green","yellow","red"/)
    
  res@xyMarkLineModes       = "MarkLines"  ; markers and lines
  res@xyMonoMarkLineMode=True
  res@xyMarkers            = 16           ; filled dot
  res@xyMarkerSizeF        = 0.0135
  res@xyMarkerThicknessF   = 0.6
  res@xyMonoDashPattern =True
  res@xyDashPattern = 0
  res@xyLineColors = (/"blue","green","orange","burlywood4","red" ,"blueviolet" ,"black","azure4"/)
    ; plot = gsn_csm_xy (wks,in_era_new ,lo_era[0],res)    ; Create plot, but don't draw it yet

   plot = gsn_csm_xy(wks,dble2flt(in_mre),width_mre,res)    ; Create plot, but don't draw it yet
; Attach some text strings next to the markers.
  txres               = True
  txres@txFontHeightF = 0.0275
  txres@txJust        = "TopLeft"

  labels            = ""+ispan(1,4,1)
  ; labels(11)        = ""   ; Don't do the last one.

  txid = gsn_add_text(wks,plot,labels,dble2flt(in_mre),width_mre,txres)
  draw(plot)    ; Drawing the plot will also cause text strings to appear
  frame(wks)    ; Advance the frame.   


  ; res2=True
  ; ; res2@xyLineColors=(/"blue",""/)
  ; plot1=gsn_csm_xy(wks1,ispan(1, 12, 1),(/lo_20cr[0],lo_era[0],lo_ncep1[0],lo_ncep2[0],lo_jra55[0],lo_jra55[0]/),res2)
  ; plot2=gsn_csm_xy(wks2,ispan(1, 12, 1),(/in_20cr,in_era,in_ncep1,in_ncep2,in_jra55,in_jra55/),res2)

  end 
   