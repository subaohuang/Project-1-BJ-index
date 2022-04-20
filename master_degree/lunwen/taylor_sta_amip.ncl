;draw_msf_hov
load "/mnt/d/ncl/lib/ncarg/nclscripts/csm/taylor_diagram.ncl"   
load "/mnt/d/ncl/lib/ncarg/nclscripts/csm/taylor_diagram_lin.ncl" 
;
begin
;画时间序列图的时候，把BCC-CSM1-1给去掉
models=(/"ACCESS1-0","ACCESS1-3","BCC-CSM1-1","BCC-CSM1-1-M","BNU-ESM","CanAM4","CCSM4","CMCC-CM","CNRM-CM5","CSIRO-Mk3-6-0","EC-EARTH","FGOALS-g2","FGOALS-s2","GFDL-CM3","GFDL-HIRAM-C180","HadGEM2-A","INMCM4","IPSL-CM5A-LR",\
    "IPSL-CM5A-MR","IPSL-CM5B-LR","MIROC5","MPI-ESM-LR","MPI-ESM-MR","MRI-AGCM3-2H","MRI-CGCM3","NorESM1-M","20cr",\
    "era","jra55","MERRA","ncep1","ncep2","mre"/)

rea=(/"20cr","era","jra55","MERRA","ncep1","ncep2","mre"/)
; models=(/"BNU-ESM"/)
dim=dimsizes(models)


msf_total=new((/dim,12,144/), float)
; msf_total=new((/dim,12,57/), float)
; msf_total=new((/dim,12,97/), float)
; msf_total=new((/dim,12,28/), float)
; msf_total=new((/dim,12,105/), float)


yearstart=197901
yearend=200812  
year1=1979
year2=2008
; p  = 0.975
p=0.95
year=year2-year1+1
locationxi_rea=new((/dim,year/), float)
;******************************************************************************************************
;循环读入数据
;******************************************************************************************************
do i=0,dimsizes(models)-1
 ; model="CanESM2"
 model=models(i)
 print("************************"+model+"*****************************")
 if any(rea.eq.model) then 
 diri="/mnt/f/share/rea_msf/"
 fils1=systemfunc("ls "+diri+"msf_time_"+model+"_global_new.nc")
 else
 diri="/mnt/f/share/CMIP5/msf_AMIP/"
 fils1=systemfunc("ls "+diri+"msf_time_"+model+"_amip_r1i1p1_197901-200812.nc")
 end if



;;;********************************************* 
    f1=addfile(fils1, "r")
    date=cd_calendar(f1->time, 1)
    timeind=ind(date.le.yearend.and.date.ge.yearstart)
    ; timeind=ind(date.ge.197901.and.date.le.201212)
    ; msf_time=f1->msf_time({850:500},:,timeind)
    z0=(/1000.0,925.0,850.0,700.0,600.0,500.0,400.0,300.0,250.0,200.0,150.0,100.0/)
    levind=get1Dindex(f1->level,z0)
   
    msf_time_r= f1->msf_time(levind,:,timeind)
    delete(f1)

    msf_time_r=(/msf_time_r/10^9/)
    msf_time=tofloat(msf_time_r)
    copy_VarMeta(msf_time_r, msf_time)
    printVarSummary(msf_time)
    delete(msf_time_r)
    delete(msf_time&time)
    msf_time&time=date(timeind)
    delete(date)
    ; msf_mean=level_mean(msf_time(:,:,:))
    ; msf_mean=dim_avg_n_Wrap(msf_time(time|:,level|:,lon|:),1)
    msf_mean=msf_time
    printVarSummary(msf_mean)
    delete(msf_time)

    ;;;先得到秋季和冬季eninio年的沃克环流情况SONDJF
    msf_son=month_to_season(msf_mean(time|:,level|:,lon|:), "SON")
    msf_djf=month_to_season(msf_mean(time|:,level|:,lon|:), "DJF")
    printVarSummary(msf_son)
    msf_son&time=ispan(year1,year2,1)
    msf_djf&time=ispan(year1,year2,1)  ;;;因为这里第一个是1958 1月和2月的平均，第二个是195812,195901,195902的平均
    ;;注意点
    msf_djf(0:dimsizes(msf_djf&time)-2,:,:)=(/msf_djf(1:,:,:)/)
    msf_djf(dimsizes(msf_djf&time)-1,:,:)=(/msf_mean(:,:,{yearend})/)
    ;;;
    ; msf_djf=msf_son
    ; delete(msf_mean)
    delete(msf_djf@_FillValue_original)
    ; msf_djf_clm = dim_avg_n_Wrap(msf_djf, 0)
    msf_djf_clm = dim_avg_n_Wrap(msf_mean(time|:,level|:,lon|:), 0)
    delete(msf_mean)


    delete(msf_djf)
    delete(msf_son)
    msf_total(i,:,:)=msf_djf_clm(:,:)
    delete(msf_djf_clm)
end do

msf_amip=dim_avg_n_Wrap(msf_total(0:dim-8,:,:), 0)
printVarSummary(msf_amip)
; msf_mre=msf_total(dim-1,:,:)
msf_mre=dim_avg_n_Wrap(msf_total(dim-7:dim-2,:,:),0)
printVarSummary(msf_mre)


pa_cor=new(dim,float)
ratio=new(dim,float)

;;权重系数
plev=msf_amip&level
psfc=1013.0
ptop=min(plev)
dp   = dpres_plevel_Wrap(plev, psfc, ptop, 0)  ; dp(30)
printVarSummary(dp)
w=dp/sum(dp)
printVarSummary(w)
;;;

do i=0,dim-2
  taylor_stat = taylor_stats(msf_total(i,:,:),msf_mre,w,0)
  pa_cor(i)=taylor_stat(0)
  ratio(i)=taylor_stat(1)
  delete(taylor_stat)
end do

  taylor_stat = taylor_stats(msf_amip,msf_mre,w,0)
  pa_cor(dim-1)=taylor_stat(0)
  ratio(dim-1)=taylor_stat(1)
  delete(taylor_stat)

;;pattern_cor
case=models(0:dim-2)
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
  ty_opt@markerTxYOffset = 0.02                    ; offset btwn marker & label
  ty_opt@varLabelsYloc = 1.5                ; Move location of variable labels [default 0.45]
  ty_opt@caseLabelsFontHeightF = 0.14       ; make slight larger   [default=0.12 ]
  ty_opt@varLabelsFontHeightF  = 0.011      ; make slight smaller  [default=0.013]
 
  ty_opt@tiMainString  = "MSF"           ; title
  ty_opt@stnRad        = (/ 0.5, 1.5 /)     ; additional standard radii
  ty_opt@ccRays        = (/ 0.85, 0.95 /)     ; correllation rays
  ty_opt@centerDiffRMS = True               ; RMS 'circles'

  wks   = gsn_open_wks("pdf", "/mnt/f/share/era/data/泰勒图_amip")
  plot  = taylor_diagram_lin(wks,ra,cc,ty_opt)


;;***********************************************************************
;;;添加注释
;;;**********************************************************************
labels2=(/"a: ACCESS1-0","b: ACCESS1-3","c: BCC-CSM1-1","d: BCC-CSM1-1-M",\
         "e: BNU-ESM","f: CanAM4","g: CCSM4","h: CMCC-CM","i: CNRM-CM5","j: CSIRO-Mk3-6-0","k: EC-EARTH","l: FGOALS-g2","m: FGOALS-s2",\
         "n: GFDL-CM3","o: GFDL-HIRAM-C180","p: HadGEM2-A","q: INMCM4","r: IPSL-CM5A-LR",\
    "s: IPSL-CM5A-MR","t: IPSL-CM5B-LR","u: MIROC5","v: MPI-ESM-LR","w: MPI-ESM-MR","x: MRI-AGCM3-2H","y: MRI-CGCM3","z: NorESM1-M"/)




    txres2               = True
    txres2@txFontHeightF = 0.010
    txres2@txFontColor="black"
    txres2@txJust="TopLeft"

    txres2@txFontThicknessF=5.0
    txres2@txFontHeightF = 0.015
    txres2@txFontColor="black"

     gsn_text_ndc (wks,"AMIP",0.72,0.94,txres2)

  dim2=dimsizes(labels2)
  ytxt2=fspan(0.3, 0.9, dim2)
  ytxt2=ytxt2(::-1)
  xtxt=0.72
  txres2@txFontHeightF = 0.010
  txres2@txFontColor="black"
  do i=0,dim2-1
  gsn_text_ndc (wks,labels2(i),xtxt,ytxt2(i),txres2)
  end do

 ;;;;再分析的注释
     index_maker=(/2,15,5,6,7,8,16/)
  labels3=(/"20CR","ERA-Interium","JRA55","MERRA","NCEP1","NCEP2","MME"/)
  xleg = 0.86
  xtxt_rea = 0.88  ; legend markers
  yleg = fspan(0.65, 0.9, 7)
  yleg=yleg(::-1)
  ytxt_rea = yleg  ; strings.

 do i=0,6
  gsres               = True
  gsres@gsMarkerIndex = index_maker(i)         ; Use filled dots for markers.
  if i.eq.6 then 
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
    txres2@txFontColor="red"
  end if
  gsn_text_ndc      (wks,labels3(i),xtxt_rea,ytxt_rea(i),txres2)
end do
  
  
     txres2@txFontColor="red"
  txres2@txFontHeightF = 0.015
  txres2@txJust="TopLeft"
  gsn_text_ndc (wks,"Reanalysis",0.86,0.94,txres2)

;;;;************************************************************************************
draw(plot)
frame(wks)

end 
  