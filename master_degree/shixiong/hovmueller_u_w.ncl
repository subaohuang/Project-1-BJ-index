begin
diri="/mnt/f/share/era/data/shixiong/"
fils=systemfunc("ls "+diri+"era_iod_pacific_w.nc")
f1=addfile(fils, "r")
time=f1->time
date=cd_calendar(time, 1)
tind=ind(date.le.198312.and.date.ge.198201)
w_r=f1->w(tind,{20:-20},{30:280})

w=(w_r@scale_factor*w_r+w_r@add_offset)*100 - 2.0
copy_VarMeta(w_r, w)

delete(time)
delete(date)
delete(tind)
;;;读入纬向风数据
fils1=systemfunc("ls "+diri+"era_iod_pacific_u10.nc")
f2=addfile(fils1, "r")
time=f2->time
date=cd_calendar(time, 1)
tind=ind(date.le.198312.and.date.ge.198201)
u10_r=f2->u10(tind,{20:-20},{30:280})
u10=u10_r@scale_factor*u10_r+u10_r@add_offset
copy_VarMeta(u10_r, u10)
;;;去除掉陆地上的纬向风的影响
   fils2=systemfunc("ls "+diri+"landsea.nc")
   a=addfile(fils2, "r")
   lsdata = a->LSMASK
   lsm  = landsea_mask(lsdata,u10&latitude,u10&longitude)   ;;得到的结果，1表示land   
   u10 = mask(u10,lsm.eq.1,False)      
   
  
;;;------------------------
w_avg=dim_avg_n_Wrap(w, 1)
u10_avg=dim_avg_n_Wrap(u10, 1)

w_avg&time=ispan(1,24,1)
u10_avg&time=ispan(1,24,1)
printVarSummary(w_avg)

printMinMax(u10_avg, 1)
v10_avg=new(dimsizes(u10_avg),typeof(u10_avg))
v10_avg=0.0


;******************************
;画图
;**********************************
 wks  = gsn_open_wks ("png", ""+diri+"hov" )           ; send graphics to PNG file
 res            = True                         ; plot mods desired
 res@gsnDraw    = False                        ; don't draw
 res@gsnFrame   = False                        ; don't advance frame
 ; res@cnLineLabelsOn =False
 res@cnFillOn=False
 res@cnMonoLineColor  = False           ; Tells NCL not to draw contour lines in one color
 res@cnLineColors     = span_color_rgba ("BlueDarkRed18",11)   ; NCV_jet has 256 colors; span it to get 11 colors
 res@cnLevelSelectionMode = "ManualLevels"     ; manually select contours
 res@cnMinLevelValF       = -5.               ; min level
 res@cnMaxLevelValF       =  5.               ; max leve
 res@cnLevelSpacingF      =1.0
 res@tmYLMode      = "Explicit"                ; Define own tick mark labels.
 ; res@tmYLValues    = ispan(1,24,)
 ; res@tmYLLabels    = (/"DEC","JAN","FEB","MAR" ,"APR","MAY" /)
 
  plot = gsn_csm_hov(wks, w_avg, res)   ; create plot

 res_line =True
 res_line@gsnDraw =False
 res_line@gsnFrame =False
 res_line@vcMapDirection=False
 res_line@vcMinDistanceF=0.03
 

 plot2=gsn_csm_streamline(wks,u10_avg,v10_avg,res_line)


draw(plot)
frame(wks)
end