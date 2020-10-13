;draw_sst_hov
load "/mnt/d/ncl/lib/ncarg/nclscripts/csm/taylor_diagram.ncl"   
load "/mnt/d/ncl/lib/ncarg/nclscripts/csm/taylor_diagram_lin_sst _scs.ncl" 
load "/mnt/d/ncl/Dot-Matrix-master/WCHAR.ncl" ; All related function tar into this lib 
;
begin
;画时间序列图的时候，把BCC-CSM1-1给去掉

models=(/"BNU-ESM","CanESM2", "CESM1-BGC", "CESM1-CAM5","CMCC-CMS","FIO-ESM","GFDL-ESM2G","GFDL-ESM2M","GISS-E2-H",\
    "GISS-E2-H-CC","GISS-E2-R","GISS-E2-R-CC","HadGEM2-AO","HadGEM2-CC",\
    "HadGEM2-ES","MIROC-ESM-CHEM","NorESM1-ME","Hadsst","ersst"/)
;;;
;;
rea=(/"Hadsst","ersst"/)

dim=dimsizes(models)

; sst_total=new((/dim,25,144/), float)
sst_total=new((/dim,25,85/), float)

yearstart=197901
yearend=200812  
year1=1979
year2=2008
; p  = 0.975
p=0.95
year=year2-year1+1

;******************************************************************************************************
;循环读入数据
;******************************************************************************************************
do i=0,dimsizes(models)-1
 ; model="CanESM2"
 model=models(i)
 print("************************"+model+"*****************************")
 if any(rea.eq.model) then 
   diri="/mnt/f/share/ersst/"
  if model.eq."ersst"  then
 fils1=systemfunc("ls "+diri+"ersst_v4b.chazhi.nc")
 else
 fils1=systemfunc("ls "+diri+"HadISST_sst_chazhi_2012.nc")
end if
;;;********************************************* 
    f1=addfile(fils1, "r")
    delete(fils1)
    date=cd_calendar(f1->time, 1)
    timeind=ind(date.le.yearend.and.date.ge.yearstart)
    sst_time= f1->sst(timeind,{-30:30},{80:290})
    delete(f1)
    printVarSummary(sst_time)
    delete(sst_time&time)
    sst_time&time=date(timeind)
    delete(date)
    ; sst_mean=level_mean(sst_time(:,:,:))
    ; sst_mean=dim_avg_n_Wrap(sst_time(time|:,level|:,lon|:),1)
    sst_mean=sst_time
    printVarSummary(sst_mean)
    delete(sst_time)
    delete(timeind)
  ;;;******************************************************

    else

    f1=addfile("/mnt/f/share/CMIP5/tos_chazhi/tos_Omon_"+model+"_chazhi_r1i1p1_197901-201212.nc", "r")
    date=cd_calendar(f1->time, 1)
    timeind=ind(date.le.yearend.and.date.ge.yearstart)
    sst_time= f1->sst(timeind,{-30:30},{80:290})
    delete(f1)
    printVarSummary(sst_time)
    delete(sst_time&time)
    sst_time&time=date(timeind)
    delete(date)
    sst_mean=sst_time
    printVarSummary(sst_mean)
    delete(sst_time)
    delete(timeind)
    end if
    sst_clm=dim_avg_n_Wrap(sst_mean, 0)
    delete(sst_mean)
    sst_total(i,:,:)=sst_clm(:,:)
    delete(sst_clm)
end do

sst_cmip=dim_avg_n_Wrap(sst_total(0:dim-3,:,:), 0)
printVarSummary(sst_cmip)
; sst_mre=sst_total(dim-1,:,:)
sst_mre=dim_avg_n_Wrap(sst_total(dim-2:,:,:), 0)
printVarSummary(sst_mre)


pa_cor=new(dim+1,float)
ratio=new(dim+1,float)

;;权重系数

 rad    = 4.0*atan(1.0)/180.0
 clat   = cos(sst_cmip&lat*rad)
 w = clat            ; clat(nlat) ==> cos(rad*lat) 

do i=0,dim-1
  taylor_stat = taylor_stats(sst_total(i,:,:),sst_mre,w,0)
  pa_cor(i)=tofloat(taylor_stat(0))
  ratio(i)=tofloat(taylor_stat(1))
  delete(taylor_stat)
  print("*********i=************"+i)
end do

  taylor_stat = taylor_stats(sst_cmip,sst_mre,w,0)
  pa_cor(dim)=tofloat(taylor_stat(0))
  ratio(dim)=tofloat(taylor_stat(1))
  delete(taylor_stat)

;;pattern_cor
ncase=dim+1
cc=new((/ncase,1/),typeof(pa_cor))
;;ratio
ra=new((/ncase,1/),typeof(ratio))
;;
cc=pa_cor
ra=ratio
Markers=new(ncase,integer)
Markers=16
Colors=new(ncase,string)
Colors="white"
;**********************************
; create plot
;**********************************

  ty_opt   = True                           ; taylor diagram with options
  ty_opt@taylorDraw   =False
  ty_opt@taylorFrame  =False
  ty_opt@Markers       = Markers         ; make all solid fill

  ty_opt@Colors        = Colors          
  ty_opt@varLabels     = ""
  ty_opt@caseLabels    = ""
  ty_opt@markerTxYOffset = 0.02                    ; offset btwn marker & label
  ty_opt@varLabelsYloc = 1.5                ; Move location of variable labels [default 0.45]
  ty_opt@caseLabelsFontHeightF = 0.14       ; make slight larger   [default=0.12 ]
  ty_opt@varLabelsFontHeightF  = 0.011      ; make slight smaller  [default=0.013]
 
  ty_opt@tiMainString  = "SST"           ; title
  ty_opt@stnRad        = (/ 0.5, 1.5 /)     ; additional standard radii
  ty_opt@ccRays        = (/ 0.90, 0.99 /)     ; correllation rays
  ty_opt@centerDiffRMS = True               ; RMS 'circles'

  wks   = gsn_open_wks("pdf", "/mnt/f/share/era/data/scs_fig_exam.2b")
  plot  = taylor_diagram_lin(wks,ra,cc,ty_opt)


;;***********************************************************************
;;;添加注释
;;;**********************************************************************
alfa= ""+ispan(1, 17, 1) +""
labels2= alfa + ":" +(/"BNU-ESM","CanESM2", "CESM1-BGC", "CESM1-CAM5","CMCC-CMS","FIO-ESM","GFDL-ESM2G","GFDL-ESM2M","GISS-E2-H",\
    "GISS-E2-H-CC","GISS-E2-R","GISS-E2-R-CC","HadGEM2-AO","HadGEM2-CC",\
    "HadGEM2-ES","MIROC-ESM-CHEM","NorESM1-ME"/)


    txres2               = True
    txres2@txFontHeightF = 0.010
    txres2@txFontColor="black"
    txres2@txJust="TopLeft"

    txres2@txFontThicknessF=5.0
    txres2@txFontHeightF = 0.015
    txres2@txFontColor="black"

     gsn_text_ndc (wks,"CMIP5",0.72,0.94,txres2)

  dim2=dimsizes(labels2)
  ytxt2=fspan(0.3, 0.9, dim2)
  ytxt2=ytxt2(::-1)
  xtxt=0.72
  txres2@txFontHeightF = 0.015
  txres2@txFontColor="black"
  do i=0,dim2-1
  gsn_text_ndc (wks,labels2(i),xtxt,ytxt2(i),txres2)
  end do

 ;;;;再分析的注释
  index_maker=(/7,8,16/)
  labels3=(/"HadISST","ERSST","MME"/)
  xleg = 0.90
  xtxt_rea = 0.92  ; legend markers
  yleg = (/0.9,0.86,0.82/)

  ytxt_rea = yleg  ; strings.

  do i=0,2
  gsres               = True
  gsres@gsMarkerIndex = index_maker(i)         ; Use filled dots for markers.
  if i.eq.2 then 
  gsres@gsMarkerColor      = "green"
  else
  gsres@gsMarkerColor      = "red"
  end if
  gsres@gsMarkerThicknessF = 0.015
  gsres@gsMarkerSizeF=10
   
; Add marker and text for the legend.
  gsn_polymarker_ndc(wks, xleg,yleg(i),gsres)
   txres2@txJust="CenterLeft"
   if i.eq.2 then
    txres2@txFontThicknessF=5.0
    txres2@txFontHeightF = 0.015
   txres2@txFontColor="green"
   else
    txres2@txFontColor="red"
  end if
  gsn_text_ndc      (wks,labels3(i),xtxt_rea,ytxt_rea(i),txres2)
  end do

  txres2@txFontColor="black"
  txres2@txFontHeightF = 0.03
  txres2@txJust="TopLeft"
  gsn_text_ndc (wks,"b",0.1,0.95,txres2)

   ;;;;添加中文注释
    fontres=True
    fontres@FontQuality="High"
    fontres@fontAngleF=90
    fontres@fontSize=14
    fontres@fontNDC=True

    fontX(wks,plot, "标准偏差", 0.04,0.6, fontres)
    fontres@fontAngleF=-45
    fontres@fontNDC=False
    fontX(wks,plot, "相关系数", 1.30,1.30, fontres)



;;;;************************************************************************************
draw(plot)
frame(wks)

end 
  