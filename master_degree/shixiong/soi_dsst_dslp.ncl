function   xiangguan(intension,location)
local r1,prob,yes_no,siglvl
begin 
  r1=escorc_n(intension, location, 0, 0)
  prob=rtest(r1,dimsizes(intension),0)
  siglvl=0.05
  yes_no = where(prob.lt.siglvl, 1.0, 0.0)
   ; print("500hpa以上msf_r1=" + r1 + "yes_no=" +yes_no)
   result=new(2,double)
   result(0)=r1
   result(1)=yes_no
return result
end 

begin 
;*********************************************
diri="/mnt/f/share/era/data/"
fils=systemfunc("ls "+diri+"era_u_sst_slp.nc")
f1=addfile(fils, "r")
time=f1->time
date=cd_calendar(time, 1)
tind=ind(date.le.201412.and.date.ge.197901)
sst_r=f1->sst(tind,{10:-10},:)
msl_r=f1->msl(tind,{5:-5},:)
u10_r=f1->u10(tind,{10:-10},:)

sst=sst_r@scale_factor*sst_r+sst_r@add_offset-274.15
copy_VarMeta(sst_r, sst)
 printVarSummary(msl_r)
msl=(msl_r@scale_factor*msl_r+msl_r@add_offset)/100.0
copy_VarMeta(msl_r, msl)

u10=u10_r@scale_factor*u10_r+u10_r@add_offset
copy_VarMeta(u10_r, u10)
    

;------------------------------------------------
;读入soi
;-----------------------------------------------
diri1="/mnt/f/share/era/data/shixiong/"
params=systemfunc("ls "+diri1+"SOI.txt")
ncols =13
nrows=67
nhead=2
  rsoi=readAsciiTable(params, ncols, "float", nhead)
  soi=rsoi(:,1:)
  soi@_FillValue=-999.9
  soi!0="yr"
  soi!1="month"
  soi&yr=ispan(1951, 2017, 1)
  soi_month=ndtooned(soi({1979:2014},:))
 printVarSummary(soi_month)
 delete(params)
;-----------------------------------------------------
 diri2="/mnt/f/share/era/data/shixiong/"
params=systemfunc("ls "+diri1+"nino3.4.txt")
ncols =5
nrows=812
nhead=1
rnino3_4=readAsciiTable(params, ncols, "float", nhead)


nino3_4=rnino3_4(348:779,2)
print("nino3_4="+nino3_4(0))
printVarSummary(nino3_4)
;--------------------------------------------
  ; diri2="/mnt/f/share/era/data/"
  ; nino3_4=asciiread(diri2+"era_nino3_4.ascii", 432, "float")
  ; printVarSummary(nino3_4)
  ; nino3_4!0="time"
  ; nino3_4&time=ispan(1,432,1)




;********************************************
 ;dslp:取5°S-5°N,100°-165°E以及 165E-290
;********************************************
dslp=dim_avg_n(msl(:,:,{195:290}),(/1,2/))-dim_avg_n(msl(:,:,{100:195}),(/1,2/))

;********************************************
 ;dsst:取10°S-10°N,90°-150°E以及 90°-150°W
;********************************************
dsst=dim_avg_n(sst(:,:,{210:270}),(/1,2/))-dim_avg_n(sst(:,:,{90:150}),(/1,2/))


;********************************************
 ;东风强度:取10°S-10°N，100-290°E
;********************************************
wind=-1*dim_avg_n(u10, (/1,2/))



;********************************************
 ;相关性:
;********************************************
 sst_slp=xiangguan(dsst,dslp)
 sst_soi=xiangguan(dsst,soi_month)
 sst_wind=xiangguan(dsst,wind)
 slp_soi=xiangguan(dslp,soi_month)
 slp_wind=xiangguan(dslp,wind)
 soi_wind=xiangguan(soi_month,wind)
 nino3_4_dslp=xiangguan(nino3_4,dslp)
 



print("sst_slp="+sst_slp)
print("sst_soi="+sst_soi)
print("sst_wind="+sst_wind)
print("slp_soi="+slp_soi)
print("slp_wind="+slp_wind)
print("soi_wind="+soi_wind)
print("nino3_4_dslp="+nino3_4_dslp)

;-------------------------------------------------
 ;读入位置
;--------------------------------------------
 
diri3=diri+"shixiong/locationxi.ascii"
diri4=diri+"shixiong/locationdong.ascii"

locationxi=asciiread(diri3, 432, "float") 
locationdong=asciiread(diri4, 432, "float")

 mxlag    = 12
 x_Lead_y = esccr(locationxi,nino3_4,mxlag)
 y_Lead_x = esccr(nino3_4,locationxi,mxlag)    ; switch the order of the series
 
 ccr = new ( 2*mxlag+1, double)    
 ccr(0:mxlag-1) = y_Lead_x(1:mxlag:-1)  ; "negative lag", -1 reverses order
 ccr(mxlag:)    = x_Lead_y(0:mxlag)     ; "positive lag"

 print(ccr(16))

 x=ispan(-12,12,1)     
 wks=gsn_open_wks("png", ""+diri2+"位置和强度相关性")
 res=True
 res@gsnXRefLine=(/4,0/)
 res@vpHeightF=0.3
 res@vpWidthF=0.7
 plot=gsn_csm_xy(wks,x,ccr,res)
 

 ;---------------------
 dslp!0="time"
 dslp&time=ispan(1,432,1)
 copy_VarCoords(dslp, nino3_4)
 printVarSummary(dslp)
 printVarSummary(nino3_4)


;--------------------------------------------
;画nino3.4与slp的时间序列图：
;---------------------------------------

    resxy=True
    resxy@gsnDraw=False
    resxy@gsnFrame=False
    resxy@gsnYRefLine=(/27/)
    resxy@gsnYRefLineThicknesses=1.5
    resxy@tmYLMode="Explicit"
    resxy@trYMaxF=30
    ; resxy@trYMinF=24
    resxy@tmYLValues=ispan(24,30,1)
    resxy@tmYLLabels=ispan(24,30,1)
    resxy@tmYLMinorOn=False
    
    
    resxy@tmXBMode="Explicit"
    
    resxy@trXMaxF=dimsizes(dslp&time)
    resxy@tmXBValues=ispan(13,dimsizes(dslp&time),60)
    resxy@tmXBLabels=ispan(1980,2015,5)
    resxy@tmXBMinorValues=ispan(1,432,12)
    
    resxy@vpHeightF=0.3
    resxy@vpWidthF=0.7
    resxy@vpXF=0.15
    resxy@vpYF=0.9
    ; resxy@xyLineThicknesses = (/ 2.0,1.0,2.0/)          ; make second line thicker
    resxy@xyLineColor      = (/"black"/)          ; change line color
    resxy@xyMarkLineModes   = (/"lines"/)
    ; resxy@xyExplicitLegendLabels = (/"Anual change", "mon_change", "Linear trend"/)   

    ; resxy@xyDashPatterns    =(/0,0,2/)
    
    resxy@tiYAxisString ="nino3.4" 
    resxy@tmXBLabelFontAspectF=2.5
    resxy@tmYLLabelFontAspectF=2.5
    resxy@tmXTOn=False
    resxy@tmBorderThicknessF=3.0     ;调整边框粗细
    resxy@tmXBLabelFontAspectF=1.0 
    resxy@tmXBLabelFontHeightF=0.008
    resxy@tmXBLabelFont=30     ;黑体
    resxy@tmXBLabelFontThicknessF=2.0
    resxy@tmXBMajorLengthF=0.01        ;調整刻度的长度
    resxy@tmXBMajorOutwardLengthF=0.01  ;调整刻度的长度，注意2个命令要同时使用

    resxy@tmXBMinorLengthF  =0.005        ;調整刻度的长度
    resxy@tmXBMinorOutwardLengthF =0.005  ;调整刻度的长度，注意2个命令要同时使用
    
    resxy@tmYLLabelFontHeightF=0.01
    resxy@tmYLLabelFont=30     ;黑体
    resxy@tmYLLabelFontThicknessF=2.0

    resxy@xyLineThicknesses=(/3.0/)
    resxy@pmLegendDisplayMode    = "Always"            ; turn on legend
 
  ; resxy@pmLegendParallelPosF   = 0.45                 ; move units right
  ; resxy@pmLegendOrthogonalPosF = -0.4                ; move units down

  ; resxy@pmLegendWidthF         = 0.015                ; Change width and
  ; resxy@pmLegendHeightF        = 0.018                ; height of legend.
  ; resxy@lgLabelFontHeightF     = .03                 ; change font height
  ; resxy@lgTitleOn              = True                ; turn on legend title

  ; resxy@lgTitleFontHeightF     = .025                ; font of legend title
 
  ;   ; resxy@lgPerimOn              = False               ; turn off box around
  ; resxy@lgLabelFontHeightF     = .008                 ; label font height

  resxy@pmLegendOrthogonalPosF = -1.28    ; Move into plot
  resxy@pmLegendParallelPosF   = 0.858       ; Move to right
  resxy@pmLegendWidthF         = 0.1       ; Decrease width
  resxy@pmLegendHeightF        = 0.08       ; Decrease height
  resxy@lgBoxMinorExtentF      = 0.2       ; Shorten the legend lines
  resxy@lgLabelFontHeightF     = 0.015     ; Change the font size
  resxy@lgPerimOn=False

  resxy@xyExplicitLegendLabels = (/"nino3.4"/)         ; create explicit labels
    
   
    resxy2=resxy
    resxy2@xyLineColor      = (/"red"/)  
    resxy2@xyLineThicknesses=(/3.0/)
    resxy2@tmYRMode ="Explicit"
    resxy2@trYMinF=-2
    resxy2@trYMaxF=6
    resxy2@tmYRValues=ispan(-2,6,2)
    resxy2@tmYRLabels=ispan(-2,6,2)
    resxy2@tmYRMinorOn=False
    resxy2@tmYRLabelFontHeightF=0.01
    resxy2@tmYRLabelFont=30     ;黑体
    resxy2@tmYRLabelFontThicknessF=2.0
    
    resxy2@tiYAxisString ="dslp"


    resxy2@pmLegendOrthogonalPosF = -1.1    ; Move into plot
    resxy2@pmLegendParallelPosF   = 0.835       ; Move to right
    resxy2@pmLegendWidthF         = 0.1       ; Decrease width
    resxy2@pmLegendHeightF        = 0.08       ; Decrease height
    resxy2@lgPerimOn=False
    resxy2@lgBoxMinorExtentF      = 0.2      ; Shorten the legend lines
    resxy2@lgLabelFontHeightF     = 0.015     ; Change the font size
    resxy2@xyExplicitLegendLabels = (/"dslp"/)         ; create explicit labels
   
    



   wks2=gsn_open_wks("png", ""+diri2+"nino3_4&dslp")

   plot2=gsn_csm_xy2(wks2, dslp&time, runave_n(nino3_4,5,0,0), runave_n(dslp,5,0,0), resxy, resxy2)
draw(plot2)
frame(wks2)
;-----------------------------
;画上升支位置和dslp的位置关系,由超前滞后相关系数可知，最大的是西边位置超前slp2个月，相关系数0处为-0.46、-0.477
;超前2个月的为
;--------------------------
 
  resxy3=resxy
  delete(resxy3@tmYLValues)
  delete(resxy3@tmYLLabels)
  delete(resxy3@tmXBLabels)
  delete(resxy3@tmXBValues)
  delete(resxy3@tmXBMinorValues)
  delete(resxy3@trYMaxF)
  delete(resxy3@gsnYRefLine)
  resxy3@xyLineColor=(/"blue"/)
  resxy3@trYMinF=100
  resxy3@tmYLValues=ispan(100,220,20)
  resxy3@tmYLLabels=""+ispan(100,220,20)+"~S~O~N~E"
  resxy3@trXMaxF=dimsizes(dslp&time)
  resxy3@tmXBValues=ispan(13,dimsizes(dslp&time),60)
  resxy3@tmXBLabels=ispan(1980,2015,5)
  resxy3@tmXBMinorValues=ispan(1,432,12)
  resxy3@gsnYRefLine=avg(locationxi)
   resxy3@tiYAxisString ="WL" 
  resxy3@xyExplicitLegendLabels = (/"WL_index"/)         ; create explicit labels
  resxy3@gsnLeftString="(a)"


  resxy4=resxy2
  resxy4@trYMinF=-2
  resxy4@trYReverse=True
  resxy4@pmLegendOrthogonalPosF = -1.23    ; Move into plot
  resxy4@pmLegendParallelPosF   = 0.7       ; Move to right
  wks3=gsn_open_wks("pdf", ""+diri2+"WL&dslp")


  

  plot3=gsn_csm_xy2(wks3, dslp&time, runave_n(locationxi,5,0,0), runave_n(dslp  ,5,0,0), resxy3, resxy4)
draw(plot3)
frame(wks3)

time3=ispan(1,428,1)

   resxy5=resxy3
   resxy5@tiYAxisString ="WL2(advanced)" 
   resxy5@gsnLeftString="(b)"
   resxy5@tmXBMinorValues=ispan(1,428,12)
   resxy5@xyExplicitLegendLabels = (/"WL2"/)  
   delete(resxy4@tmYRValues)
   delete(resxy4@tmYRLabels)
   delete(resxy4@gsnYRefLine)
   resxy6=resxy4
   resxy6@tmYRValues=ispan(-24,30,2)
   resxy6@tmYRLabels=ispan(-24,30,2)
   resxy6@trYMaxF=30
   resxy6@trYMinF=24

   resxy6@xyExplicitLegendLabels = (/"nino3.4"/)         ; create explicit labels
   resxy6@tiYAxisString ="nino3.4" 
   locationdong_advance=   locationdong(0:427) - 160.0

  wks4=gsn_open_wks("png", ""+diri2+"WL2&nino3_4")
  plot4=gsn_csm_xy2(wks4, time3, runave_n(locationdong_advance(:),5,0,0), runave_n(nino3_4(4:431),5,0,0), resxy5, resxy6)
  draw(plot4)
  frame(wks4)

end
