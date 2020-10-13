;draw_sst_hov
load "/mnt/d/ncl/lib/ncarg/nclscripts/csm/taylor_diagram.ncl"   
load "/mnt/d/ncl/lib/ncarg/nclscripts/csm/taylor_diagram_lin_sst.ncl" 
;
begin
;画时间序列图的时候，把BCC-CSM1-1给去掉

models=(/"ACCESS1-0","ACCESS1-3","BCC-CSM1-1","BCC-CSM1-1-M","CCSM4","CMCC-CM","CNRM-CM5",\
        "CSIRO-Mk3-6-0","GFDL-CM3","INMCM4","IPSL-CM5A-LR",\
    "IPSL-CM5A-MR","IPSL-CM5B-LR","MIROC5","MPI-ESM-LR","MPI-ESM-MR","MRI-CGCM3","NorESM1-M","Hadsst"/)
;;;
;;
rea=(/"Hadsst"/)

dim=dimsizes(models)

sst_total=new((/dim,60,360/), float)

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
 fils1=systemfunc("ls "+diri+"HadISST_sst.nc")


;;;********************************************* 
    f1=addfile(fils1, "r")
    date=cd_calendar(f1->time, 1)
    timeind=ind(date.le.yearend.and.date.ge.yearstart)
    sst_time= lonFlip(f1->sst(timeind,{-30:30},:))
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

    fils1=systemfunc("ls /mnt/i/CMIP5/history/tos_new/"+model+"/r1i1p1/tos_Omon_"+model+"_historical_r1i1p1_"+"*"+".nc")
    fils2=systemfunc("ls /mnt/i/CMIP5/rcp45/tos_new/"+model+"/r1i1p1/tos_Omon_"+model+"_rcp45_r1i1p1_"+"*"+".nc")

    printVarSummary(fils1)
    printVarSummary(fils2)
    print(fils1)
    print(fils2)
    
    fils=new(2,string)
    fils(0)=fils1
    fils(1)=fils2

 
    f1=addfiles(fils, "r")
    delete(fils)
    vNames= getfilevarnames(f1[0]) ; get variable names of all groups on file
    print(vNames)
    
    if (any(vNames.eq."LONN179_180")) then 
    time=f1[:]->TIME
    date=cd_calendar(time, 1)

    print("******************日期*********************")

    timeind=ind(date.le.yearend.and.date.ge.yearstart)

    if model.eq."FGOALS-s2" then
    sst_time=f1[:]->TOS(timeind,{-30:30},:)
    else
    sst_time=lonFlip(f1[:]->TOS(timeind,{-30:30},:))
    end if    
    
     if model.eq."MIROC5" then
      sst_time=where(sst_time.eq.0, sst_time@_FillValue, sst_time)
    end if
    sst_time=(/sst_time-273.15/)
    printVarSummary(sst_time)
    dNames = getvardims(f1[:]->TOS)
      delete(f1)
    print(dNames)

    sst_time!0="time"
    sst_time!1="latitude"
    sst_time!2="longitude"

    sst_time&time=date(timeind)
    delete(date)
    sst_mean=sst_time
    printVarSummary(sst_mean)
    delete(sst_time)
    delete(timeind)
    delete(vNames)
    delete(dNames)
    end if
   end if

    sst_clm=dim_avg_n_Wrap(sst_mean, 0)
    sst_total(i,:,:)=sst_clm(:,:)
    delete(sst_clm)
end do

sst_cmip=dim_avg_n_Wrap(sst_total(0:dim-2,:,:), 0)
printVarSummary(sst_cmip)
; sst_mre=sst_total(dim-1,:,:)
sst_mre=sst_total(dim-1,:,:)
printVarSummary(sst_mre)


pa_cor=new(dim,float)
ratio=new(dim,float)

;;权重系数

 rad    = 4.0*atan(1.0)/180.0
 clat   = cos(sst_cmip&latitude*rad)
 w = clat            ; clat(nlat) ==> cos(rad*lat) 

do i=0,dim-2
  taylor_stat = taylor_stats(sst_total(i,:,:),sst_mre,w,0)
  pa_cor(i)=tofloat(taylor_stat(0))
  ratio(i)=tofloat(taylor_stat(1))
  delete(taylor_stat)
  print("*********i=************"+i)
end do

  taylor_stat = taylor_stats(sst_cmip,sst_mre,w,0)
  pa_cor(dim-1)=tofloat(taylor_stat(0))
  ratio(dim-1)=tofloat(taylor_stat(1))
  delete(taylor_stat)

;;pattern_cor
ncase=dim
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
  ty_opt@markerTxYOffset = 0.0                    ; offset btwn marker & label
  ty_opt@varLabelsYloc = 1.5                ; Move location of variable labels [default 0.45]
  ty_opt@caseLabelsFontHeightF = 0.14       ; make slight larger   [default=0.12 ]
  ty_opt@varLabelsFontHeightF  = 0.011      ; make slight smaller  [default=0.013]
 
  ty_opt@tiMainString  = "SST"           ; title
  ty_opt@stnRad        = (/ 0.5, 1.5 /)     ; additional standard radii
  ty_opt@ccRays        = (/ 0.85, 0.95 /)     ; correllation rays
  ty_opt@centerDiffRMS = True               ; RMS 'circles'

  wks   = gsn_open_wks("pdf", "/mnt/f/share/era/data/泰勒图_cmip模式_热带SST")
  plot  = taylor_diagram_lin(wks,ra,cc,ty_opt)


;;***********************************************************************
;;;添加注释
;;;**********************************************************************
alfa= (/"a","b","c","d","e","f","g","h","i","j","k","l","m","n","o","p","q","r"/)
labels2= alfa + ":" + (/"ACCESS1-0","ACCESS1-3","BCC-CSM1-1","BCC-CSM1-1-M","CCSM4","CMCC-CM","CNRM-CM5",\
        "CSIRO-Mk3-6-0","GFDL-CM3","INMCM4","IPSL-CM5A-LR",\
    "IPSL-CM5A-MR","IPSL-CM5B-LR","MIROC5","MPI-ESM-LR","MPI-ESM-MR","MRI-CGCM3","NorESM1-M"/)



    txres2               = True
    txres2@txFontHeightF = 0.010
    txres2@txFontColor="black"
    txres2@txJust="TopLeft"

    txres2@txFontThicknessF=5.0
    txres2@txFontHeightF = 0.015
    txres2@txFontColor="blue"

     gsn_text_ndc (wks,"CMIP",0.72,0.94,txres2)

  dim2=dimsizes(labels2)
  ytxt2=fspan(0.3, 0.9, dim2)
  ytxt2=ytxt2(::-1)
  xtxt=0.72
  txres2@txFontHeightF = 0.010
  txres2@txFontColor="blue"
  do i=0,dim2-1
  gsn_text_ndc (wks,labels2(i),xtxt,ytxt2(i),txres2)
  end do

 ;;;;再分析的注释
   index_maker=(/16/)
  labels3=(/"MME"/)
  xleg = 0.86
  xtxt_rea = 0.88  ; legend markers
  yleg = 0.90

  ytxt_rea = yleg  ; strings.

  i=0
  gsres               = True
  gsres@gsMarkerIndex = index_maker(i)         ; Use filled dots for markers.
  if i.eq.0 then 
  gsres@gsMarkerColor      = "green"
  else
  gsres@gsMarkerColor      = "red"
  end if
  gsres@gsMarkerThicknessF = 0.015
  gsres@gsMarkerSizeF=10
   
; Add marker and text for the legend.
  gsn_polymarker_ndc(wks, xleg,yleg(i),gsres)
   txres2@txJust="CenterLeft"
   if i.eq.6 then
    txres2@txFontThicknessF=5.0
    txres2@txFontHeightF = 0.015
   txres2@txFontColor="green"
   else
    txres2@txFontColor="blue"
  end if
  gsn_text_ndc(wks,labels3(i),xtxt_rea,ytxt_rea(i),txres2)

  


;;;;************************************************************************************
draw(plot)
frame(wks)

end 
  