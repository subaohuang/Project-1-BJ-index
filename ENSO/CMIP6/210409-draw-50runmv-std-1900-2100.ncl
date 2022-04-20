
;;**************************************************************
;;;main program
;;***************************************************************
function  ploy(std1)
begin

dim  = dimsizes(std1)
mnmx  = new ((/3,dim(1)/), typeof(std1))
;;;计算标准偏差 
sigma1  = dim_stddev_n(std1, 0)

mnmx(0,0:) =  dim_avg_n_Wrap(std1, 0)
mnmx(1,0:) = mnmx(0,:) + sigma1
mnmx(2,0:) = mnmx(0,:) - sigma1

return mnmx
end 


begin
;读入数据

; model_name = (/"ACCESS-CM2","ACCESS-ESM1-5","AWI-CM-1-1-MR","BCC-CSM2-MR", "BCC-ESM1",\
;               "CAMS-CSM1-0","CanESM5","CESM2","CESM2-WACCM","CNRM-CM6-1","CNRM-ESM2-1"\
;     "E3SM-1-0", "EC-Earth3","EC-Earth3-Veg","FGOALS-f3-L",  "FGOALS-g3", \
;     "GFDL-CM4", "GFDL-ESM4",    "GISS-E2-1-G",  "GISS-E2-1-G-CC",   "GISS-E2-1-H"\
;         "HadGEM3-GC31-LL","IPSL-CM6A-LR",   "MCM-UA-1-0","MIROC6", \
;         "MIROC-ES2L","MPI-ESM1-2-HR","MRI-ESM2-0","NESM3","NorCPM1",  "NorESM2-LM"/)
model_name = (/"ACCESS-CM2","ACCESS-ESM1-5", "BCC-CSM2-MR",\
              "CAMS-CSM1-0","CanESM5", "CAS-ESM2-0" ,"CESM2","CESM2-WACCM","CIESM","CMCC-CM2-SR5","CMCC-ESM2","CNRM-CM6-1"\
   , "EC-Earth3","EC-Earth3-Veg","FGOALS-f3-L",  "FGOALS-g3", \
     "GISS-E2-1-G","HadGEM3-GC31-LL","IPSL-CM6A-LR","MIROC-ES2L","MIROC6", \
    "MPI-ESM1-2-HR","MPI-ESM1-2-LR","MRI-ESM2-0","NESM3",  "NorESM2-LM"/)

f = addfile("std_diff_SSP-his-2050-2100minus1950-2000.nc", "r")

std_diff  = f ->std_diff

std_diff_nino3 = dim_avg_n_Wrap(std_diff(:,{-5:5},{210:270}), (/1,2/))

model_ind = ind(std_diff_nino3.ge.0.0)
f1 = addfile("runmvstd_1900-2099-rmanucycl-eachperiod.nc", "r")
; f1 = addfile("annual_runmvstd_1900-2099-rmanucycl.nc", "r")
; std_all = f1->std_all(model_ind,:,:,:)
std_all = f1->std_all(:,:,:,:)

std_nino3 = dim_avg_n_Wrap(std_all(:,:,{-5:5},{210:270}), (/2,3/))


; mnmx1 = ploy(std_nino3(:,0:66))
; mnmx2 = ploy(std_nino3(:,66:))

mnmx1 = ploy(std_nino3(:,0:116))
mnmx2 = ploy(std_nino3(:,116:))


x = ispan(1850, 2050, 1)
printVarSummary(x)
printVarSummary(std_nino3)

; x1 = x(0:66)
; x2 = x(66:)

x1 = x(0:116)
x2 = x(116:)
; ind1 = ind(x.eq.1999)
; ind2 = ind(x.eq.2099)
; print(std_nino3(:,ind2) - std_nino3(:,ind1))



;;;;draw 
;;*************************************************
; create plot
;**************************************************
 wks = gsn_open_wks("pdf", "50runmvstd-1900-2100-nino3-rmancycle")
  res               = True                   ; plot mods desired
  res@gsnDraw = False
  res@gsnFrame = False
  res@tiMainString  = "50-year running periods STD of Nino3"   ; title
  res@tiMainFontHeightF = 0.02
  res@tiYAxisFontHeightF = 0.02
  res@tiXAxisFontHeightF = 0.02
  res@tiYAxisString = "SST anomaly s.d."           ; y axis title
  res@tiXAxisString = "Year"
    
  ; res@tmXBLabelFontHeightF=0.02
  ; res@tmYLLabelFontHeightF=0.02
  ; res@tmXBMajorThicknessF=1.0
  ; res@tmYLMajorThicknessF=1.0
  ; res@tmBorderThicknessF=1.0

  res@xyLineColor      = (/"blue"/)  ; line colors
  res@xyLineThicknessF = (/4.0/)        ; line thicknesses
  res@xyDashPattern    = (/0.0/)        ; line patterns

  ; res@xyMarkLineModes   = (/"Markers","Lines"/) ; markers?
  ; res@xyMarkerColors    = (/"red",    "black"/)     ; color
  ; res@xyMarkers         = (/2,0/)          ; style

  res@tmXBMode="Explicit"

  res@trXMaxF = 2050
  ; res@trXMinF = 1950
  res@trXMinF = 1850
  res@trYMinF = 0.2
  res@trYMaxF = 1.8

  res@tmXTOn = False
  res@tmYROn  = False


  res@vpWidthF = 0.3
  res@vpHeightF = 0.5

  ; res@tmXBValues=ispan(1950,2100,50)
  ; res@tmXBLabels=ispan(1950,2100,50) + ""
  ; res@tmXBMinorOn=True
  ; res@tmXBMinorValues=ispan(1950,2100,25)
  ; res@tmYRMode="Automatic"

   res@tmXBValues=ispan(1850,2050,50)
  res@tmXBLabels=ispan(1850,2050,50) + ""
  res@tmXBMinorOn=True
  res@tmXBMinorValues=ispan(1850,2050,25)
  res@tmYRMode="Automatic"
 
  
 printVarSummary(x1)
 printVarSummary(mnmx1)
 printVarSummary(x2)
 printVarSummary(mnmx2)
  mid_plot1= gsn_csm_xy(wks,x1,mnmx1(0,:),res)       ; create line plot
  res@xyLineColor      = (/"red"/)  ; line colors
  mid_plot2= gsn_csm_xy(wks,(/x2/),(/mnmx2(0,:)/),res)       ; create line plot


  
  delete(res@xyLineColors)
  res@gsnXYFillColors = "blue"
  res@xyLineColor     = -1                           ; We don't want the line, so make it transparent.
  res@gsnXYFillOpacities=0.1

  plot_std1= gsn_csm_xy(wks,x1,mnmx1(1:2,:),res)       ; create line plot



  res@gsnXYFillColors = "pink"
  res@xyLineColor     = -1                           ; We don't want the line, so make it transparent.
  ; res@gsnXYFillOpacities=0.1

  plot_std2= gsn_csm_xy(wks,x2,mnmx2(1:2,:),res)       ; create line plot
  
  overlay(mid_plot1, plot_std1)
  overlay(mid_plot2, plot_std2)

  
  ; ;;添加legend 
  ; gres = True
  ; gres@YPosPercent = 95.    ; expressed as %, 0->100, sets position of top border of legend 
  ;             ;  when gres@Position is set to its default setting of "Top" (Default = 95.)
  ; gres@XPosPercent = 8    ; expressed as %, 0->100, sets position of left border of legend(Default = 5.)

  ; lineres = True
  ; ; lineres@lgLineColors = (/"black"，/) ; line colors
  ; lineres@lgLineColors       = (/"Blue","red"/) ; colors chosen
  ; lineres@lgLineThicknesses = 4.0                       ; line thicknesses
  ; lineres@LineLengthPercent = 9.                         ; expressed as %, 0->100, length of line

  ; textres = True

  ; textres@lgLabels = (/"1900-1950","1960-2014"/)
  ; textres@lgLabelFontHeightF = 0.01
  ; ; textres@lgLabelOpacity = (/0/) 


  ; mid = simple_legend(wks,mid_plot,gres,lineres,textres)   
 draw(mid_plot1)
 draw(mid_plot2)
 frame(wks)


end 


