load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/popRemap.ncl"
; ==============================================================
opt = False  

latS   =  -20
latN   =  30
lonL   =  175 ; 175E
lonR   =  265 ; 95W

yrStrt = 251
yrLast = 349
year_count  = fspan(yrStrt,yrLast+1,(yrLast-yrStrt +1)*12+1)
years = year_count(:(yrLast-yrStrt +1)*12-1)
 ;;;time span
startmon =25101
endmon = 34912

startyear_r = str_split_by_length(startmon, 3)  ;; 将时间string分解成年和月
endyear_r = str_split_by_length(endmon, 3 )

startyear = stringtoint(startyear_r(0))
endyear = stringtoint(endyear_r(0))

yearnum = endyear - startyear +1
  ;;
begin

lsm_in = addfile("$NCARG_ROOT/lib/ncarg/data/cdf/landsea.nc","r")
lsdata = lsm_in->LSMASK
printVarSummary(lsdata)
printMinMax(lsdata, opt)

; ==============================================================
; fin = addfile("/home/yangsong3/data-observation/linshh/CESM-data/B/B2000_f19g16_CTRL/B2000_f19g16_CP_CTRL_TEMP.h1.0251-0350.nc", "r")

; sst_org = fin->sst(26*12:,:,:)

;;;;读入控制实验结果
    diri="/home/yangsong3/data-observation/linshh/CESM-data/B/B2000_f19g16_CTRL/"
    fils1=systemfunc("ls "+diri+"B2000_f19g16_CP_CTRL_TEMP.h1.0251-0350.nc")
    f1=addfile(fils1, "r")
    time = f1->time
    time = (/time - 1/)
    date=cd_calendar(time, 1)

    timeind=ind(date.le.endmon.and.date.ge.startmon)
   
    sst_r = f1->TEMP(timeind,0,:,:)  
   ;;;;差值TEMP 
    sst = PopLatLon(sst_r,"gx1v6","fv1.9x2.5","bilin","da","090206")
    printVarSummary(sst)
    
    lat = sst&lat
    lon = sst&lon
    lsm = landsea_mask(lsdata,lat,lon)
    printVarSummary(lsm)
    printMinMax(lsm, opt)


    diri2="/home/yangsong3/data-observation/linshh/CESM-data/B/B2000_f19g16_CTRL/"
    fils2=systemfunc("ls "+diri+"B2000_f19g16_CP_CTRL.cam.h1.0250-0349_chazhi.nc")
    f2=addfile(fils2, "r")
    time := f2->time
    time = (/time - 1/)
    date:=cd_calendar(time, 1)

    timeind:=ind(date.le.endmon.and.date.ge.startmon)

    uwnd = f2->U(timeind,0,:,:)
 
    vwnd = f2->V(timeind,0,:,:)


; ==============================================================
    sst_clm = clmMonTLL(sst(:,{latS:latN},{lonL:lonR}))
    sst_anom = calcMonAnomTLL(sst(:,{latS:latN},{lonL:lonR}), sst_clm)
    printVarSummary(sst_anom)

    lsm_3d = conform(sst_anom, lsm({latS:latN},{lonL:lonR}), (/1,2/))
    copy_VarCoords(sst_anom, lsm_3d)
    printVarSummary(lsm_3d)

    sst_anom = where(lsm_3d.eq.0, sst_anom, sst_anom@_FillValue)
    printMinMax(sst_anom, opt)

    cti = dim_avg_n_Wrap(sst_anom(:,{-5:5},{180:270}), (/1,2/))
    cti_3d = conform(sst_anom, cti, 0)
    
    ;;;去掉ENSO 信号
    reg_cti_sst = regCoef_n(cti, sst_anom, 0, 0)
    reg_cti_sst_3d = conform(sst_anom, reg_cti_sst, (/1,2/))
    sst_anom_rmcti = sst_anom - cti_3d * reg_cti_sst_3d
    copy_VarCoords(sst_anom, sst_anom_rmcti)

    uwnd_clm = clmMonTLL(uwnd(:,{latS:latN},{lonL:lonR}))
    uwnd_anom = calcMonAnomTLL(uwnd(:,{latS:latN},{lonL:lonR}), uwnd_clm)
    printVarSummary(uwnd_anom)
    uwnd_anom = where(lsm_3d.eq.0, uwnd_anom, uwnd_anom@_FillValue)
    printMinMax(uwnd_anom, opt)

    reg_cti_uwnd = regCoef_n(cti, uwnd_anom, 0, 0)
    reg_cti_uwnd_3d = conform(uwnd_anom, reg_cti_uwnd, (/1,2/))
    uwnd_anom_rmcti = uwnd_anom - cti_3d * reg_cti_uwnd_3d
    copy_VarCoords(uwnd_anom, uwnd_anom_rmcti)

    vwnd_clm = clmMonTLL(vwnd(:,{latS:latN},{lonL:lonR}))
    vwnd_anom = calcMonAnomTLL(vwnd(:,{latS:latN},{lonL:lonR}), vwnd_clm)
    printVarSummary(vwnd_anom)
    vwnd_anom = where(lsm_3d.eq.0, vwnd_anom, vwnd_anom@_FillValue)
    printMinMax(vwnd_anom, opt)

    reg_cti_vwnd = regCoef_n(cti, vwnd_anom, 0, 0)
    reg_cti_vwnd_3d = conform(vwnd_anom, reg_cti_vwnd, (/1,2/))
    vwnd_anom_rmcti = vwnd_anom - cti_3d * reg_cti_vwnd_3d
    copy_VarCoords(vwnd_anom, vwnd_anom_rmcti)

; ==============================================================
dimsize_uv = dimsizes(uwnd_anom)
nmonths = dimsize_uv(0)

surf_wind = new((/dimsize_uv(0),dimsize_uv(1),dimsize_uv(2),2/),"float")
surf_wind(:,:,:,0) = uwnd_anom_rmcti
surf_wind(:,:,:,1) = vwnd_anom_rmcti
surf_wind!3 = "var"
printVarSummary(surf_wind)
printMinMax(surf_wind, opt)

;===========================================================================
nmca = 2 ; how many MCA we need
var1_region = sst_anom_rmcti(:,{0:latN},{lonL:lonR});把ts数据reshape成二维
var1_size = dimsizes(var1_region)
print(var1_size)
n_var1_size = var1_size(1)*var1_size(2)

var1_ano_line_org = reshape(var1_region,(/var1_size(0),n_var1_size/))
var1_ano_exclude_1d = get1Dindex_Exclude(var1_ano_line_org(0,:), var1_ano_line_org@_FillValue)
size_1d_var1 = dimsizes(var1_ano_exclude_1d)
print(size_1d_var1)
var1_ano_line = new((/var1_size(0),size_1d_var1/), "float")
do i = 0, size_1d_var1-1
  var1_ano_line(:,i) = var1_ano_line_org(:,var1_ano_exclude_1d(i))
end do
printMinMax(var1_ano_line, False)

var1_ano_line!0 = "time" 
var1_ano_line!1 = "pts"
printVarSummary(var1_ano_line)
printMinMax(var1_ano_line, False)

var2_region = surf_wind(:,{0:latN},{lonL:lonR},:);把u_jas_anom数据reshape成二维
var2_size = dimsizes(var2_region)
print(var2_size)
n_var2_size = var2_size(1)*var2_size(2)*var2_size(3)

var2_ano_line_org = reshape(var2_region,(/var2_size(0),n_var2_size/))
var2_ano_exclude_1d = get1Dindex_Exclude(var2_ano_line_org(0,:), var2_ano_line_org@_FillValue)
size_1d_var2 = dimsizes(var2_ano_exclude_1d)
print(size_1d_var2)
var2_ano_line = new((/var2_size(0),size_1d_var2/), "float")
do i = 0, size_1d_var2-1
  var2_ano_line(:,i) = var2_ano_line_org(:,var2_ano_exclude_1d(i))
end do
printMinMax(var2_ano_line, False)

var2_ano_line!0 = "time" 
var2_ano_line!1 = "pts"
printVarSummary(var2_ano_line)
printMinMax(var2_ano_line, False)

printVarSummary(var1_ano_line(pts|:,time|:))
printVarSummary(var2_ano_line(pts|:,time|:))

;===========================================================================
ntime   = nmonths                            ; # time steps
ncols   = n_var1_size                        ; # columns (stations or grid pts) for S
ncolz   = n_var2_size                        ; # columns (stations or grid pts) for Z
nsvd    = 2                                  ; # svd patterns to calculate 
                                             ; [nsvd <= min(ncols, ncolz) ]
; xmsg    = -999.9                             ; missing value

homlft = new((/nmca,size_1d_var1/),"float")
hetlft = new((/nmca,size_1d_var1/),"float")
homrgt = new((/nmca,size_1d_var2/),"float")
hetrgt = new((/nmca,size_1d_var2/),"float")

print("start svd")
x = svdstd(var1_ano_line(pts|:,time|:),var2_ano_line(pts|:,time|:),nsvd,homlft,hetlft,homrgt,hetrgt) 
print("svdcov: percent variance= " + x)
printVarSummary(x)

svLeft = new((/nmca,size_1d_var1/),"float")
svRight = new((/nmca,size_1d_var2/),"float")       

x2 = svdstd_sv(var1_ano_line(pts|:,time|:), var2_ano_line(pts|:,time|:), nsvd, svLeft, svRight)
sv = x2@sv(0:1)
print("svdcov_sv: singular values = " + sv)
printVarSummary(sv)

ak   = (/-1*onedtond(x@ak,(/nsvd,ntime/))/)
bk   = (/-1*onedtond(x@bk,(/nsvd,ntime/))/)
; ak   = onedtond(x@ak,(/nsvd,ntime/))
; bk   = onedtond(x@bk,(/nsvd,ntime/))
ak!0 := "sv"
ak!1 := "time"
bk!0 := "sv"
bk!1 := "time"

; homlft_out = new((/nmca, n_var1_size/), "float", sst@_FillValue)
; do i = 0, size_1d_var1-1
;   homlft_out(:,var1_ano_exclude_1d(i)) = homlft(:,i)
; end do
; homlft_data = onedtond(homlft_out, (/nmca,var1_size(1),var1_size(2)/))
; copy_VarCoords(var1_region(0:1,:,:), homlft_data)
; printVarSummary(homlft_data)
; printMinMax(homlft_data, False)

; homrgt_out = new((/nmca, n_var2_size/), "float", sst@_FillValue)
; do i = 0, size_1d_var2-1
;   homrgt_out(:,var2_ano_exclude_1d(i)) = homrgt(:,i)
; end do
; homrgt_data = onedtond(homrgt_out, (/nmca,var2_size(1),var2_size(2),var2_size(3)/))
; copy_VarCoords(var2_region(0:1,:,:,:), homrgt_data)
; printVarSummary(homrgt_data)
; printMinMax(homrgt_data, False)

; hetlft_out = new((/nmca, n_var1_size/), "float", sst@_FillValue)
; do i = 0, size_1d_var1-1
;   hetlft_out(:,var1_ano_exclude_1d(i)) = hetlft(:,i)
; end do
; hetlft_data = onedtond(hetlft_out, (/nmca,var1_size(1),var1_size(2)/))
; copy_VarCoords(var1_region(0:1,:,:), hetlft_data)
; printVarSummary(hetlft_data)
; printMinMax(hetlft_data, False)

; hetrgt_out = new((/nmca, n_var2_size/), "float", sst@_FillValue)
; do i = 0, size_1d_var2-1
;   hetrgt_out(:,var2_ano_exclude_1d(i)) = hetrgt(:,i)
; end do
; hetrgt_data = onedtond(hetrgt_out, (/nmca,var2_size(1),var2_size(2),var2_size(3)/))
; copy_VarCoords(var2_region(0:1,:,:,:), hetrgt_data)
; printVarSummary(hetrgt_data)
; printMinMax(hetrgt_data, False)

ccr = escorc(ak(0,:), bk(0,:))

; ak_std = dim_standardize_Wrap(ak,1)
; bk_std = dim_standardize_Wrap(bk,1)
; printVarSummary(ak_std)
; printVarSummary(bk_std)
ak_std = ak
bk_std = bk

svd_sv = sv
do ns = 0, 1
    ; ak_std(ns,:) = ak(ns,:)/sqrt(svd_sv(ns))
    ; bk_std(ns,:) = bk(ns,:)/sqrt(svd_sv(ns))
    ak_std(ns,:) = dim_standardize_n_Wrap(ak(ns,:), 1, 0)
    bk_std(ns,:) = dim_standardize_n_Wrap(bk(ns,:), 1, 0)
end do
printVarSummary(ak_std)
printVarSummary(bk_std)  

; system("rm pmm_mode_his_left.nc")
; pmm_mode_left = addfile("pmm_mode_his_left.nc", "c")
; pmm_mode_left->homleft = homlft_data
; pmm_mode_left->hetlft = hetlft_data

; system("rm pmm_mode_his_right.nc")
; pmm_mode_right = addfile("pmm_mode_his_right.nc", "c")
; pmm_mode_right->homrgt = homrgt_data
; pmm_mode_right->hetrgt = hetrgt_data

data = new((/2,nmonths/),"float")
data(0,:) = ak(0,:)
; data(1,:) = bk(0,:)
printMinMax(data, False)

data_std = new((/2,nmonths/),"float")
data_std(0,:) = ak_std(0,:)
data_std(1,:) = bk_std(0,:)

reg_uwnd = regCoef_n(ak_std(0,:), uwnd_anom_rmcti, 0, 0)
copy_VarCoords(uwnd_anom(0,:,:), reg_uwnd)
printVarSummary(reg_uwnd)
tval_uwnd = onedtond(reg_uwnd@tval, dimsizes(reg_uwnd))   ;t-statistic of reg_uwnd_m1
df_uwnd   = onedtond(reg_uwnd@nptxy, dimsizes(reg_uwnd)) - 2.  ;自由度
b_uwnd = tval_uwnd    ; b_uwnd must be same size as tval_uwnd (and df_uwnd)
b_uwnd = 0.5
prob_reg_uwnd = betainc(df_uwnd/(df_uwnd+tval_uwnd^2),df_uwnd/2.0,b_uwnd) 
copy_VarCoords(reg_uwnd, prob_reg_uwnd)

reg_vwnd = regCoef_n(ak_std(0,:), vwnd_anom_rmcti, 0, 0)
copy_VarCoords(vwnd_anom(0,:,:), reg_vwnd)
printVarSummary(reg_vwnd)
tval_vwnd = onedtond(reg_vwnd@tval, dimsizes(reg_vwnd))   ;t-statistic of reg_vwnd_m1
df_vwnd   = onedtond(reg_vwnd@nptxy, dimsizes(reg_vwnd)) - 2.  ;自由度
b_vwnd = tval_vwnd    ; b_vwnd must be same size as tval_vwnd (and df_vwnd)
b_vwnd = 0.5
prob_reg_vwnd = betainc(df_vwnd/(df_vwnd+tval_vwnd^2),df_vwnd/2.0,b_vwnd) 
copy_VarCoords(reg_vwnd, prob_reg_vwnd)

reg_sst = regCoef_n(ak_std(0,:), sst_anom_rmcti, 0, 0)
copy_VarCoords(sst_anom(0,:,:), reg_sst)
printVarSummary(reg_sst)
tval_sst = onedtond(reg_sst@tval, dimsizes(reg_sst))   ;t-statistic of reg_sst_m1
df_sst   = onedtond(reg_sst@nptxy, dimsizes(reg_sst)) - 2.  ;自由度
b_sst = tval_sst    ; b_sst must be same size as tval_sst (and df_sst)
b_sst = 0.5
prob_reg_sst = betainc(df_sst/(df_sst+tval_sst^2),df_sst/2.0,b_sst) 
copy_VarCoords(reg_sst, prob_reg_sst)

reg_uwnd_test = where(prob_reg_uwnd.lt.0.05.or.prob_reg_vwnd.lt.0.05, reg_uwnd, reg_uwnd@_FillValue)
copy_VarCoords(reg_uwnd, reg_uwnd_test)

reg_vwnd_test = where(prob_reg_uwnd.lt.0.05.or.prob_reg_vwnd.lt.0.05, reg_vwnd, reg_vwnd@_FillValue)
copy_VarCoords(reg_vwnd, reg_vwnd_test)

reg_sst_test = where(prob_reg_sst.lt.0.05, reg_sst, reg_sst@_FillValue)
copy_VarCoords(reg_sst, reg_sst_test)

; system("rm pmm_index_CTRL.nc")
; pmm_out = addfile("pmm_index_CTRL.nc", "c")

; pmm_out->x = x
; pmm_out->sv = sv
; pmm_out->sst_pmm_std = ak_std
; pmm_out->wind_pmm_std = bk_std
; pmm_out->sst_pmm = ak
; pmm_out->wind_pmm = bk

asciiwrite("pmm_index_CTRL.txt", ak_std(0,:))

; delete([/reg_sst@nptxy,reg_sst@rstd,reg_sst@yintercept,reg_sst@tval/])
; pmm_out->reg_sst = reg_sst
; pmm_out->prob_reg_sst = prob_reg_sst
; delete([/reg_uwnd@nptxy,reg_uwnd@rstd,reg_uwnd@yintercept,reg_uwnd@tval/])
; pmm_out->reg_uwnd = reg_uwnd
; pmm_out->prob_reg_uwnd = prob_reg_uwnd
; delete([/reg_vwnd@nptxy,reg_vwnd@rstd,reg_vwnd@yintercept,reg_vwnd@tval/])
; pmm_out->reg_vwnd = reg_vwnd
; pmm_out->prob_reg_vwnd = prob_reg_vwnd

system("rm pmm_reg_sst.nc")
pmm_out = addfile("pmm_reg_sst.nc", "c")
delete([/reg_sst@nptxy,reg_sst@rstd,reg_sst@yintercept,reg_sst@tval/])
pmm_out->reg_sst = reg_sst

;===========================================================================
wks  = gsn_open_wks("pdf","CTRL")
gsn_define_colormap(wks,"BlueWhiteOrangeRed")

res                             = True 
res@gsnDraw                     = False
res@gsnFrame                    = False      ;don't advance frame
res@gsnAddCyclic                = False      ;EOF data is not cyclic
res@cnFillOn                    = True         ; turn on color fill
res@cnLinesOn                   = False
res@cnLineLabelsOn              = False      ;turn off cn line labels
res@cnInfoLabelOn               = False      ;turn off contour information label
res@mpFillOn                    = True        ; turn off map fill
res@mpMinLatF                   = latS
res@mpMaxLatF                   = latN
res@mpMinLonF                   = lonL+1
res@mpMaxLonF                   = lonR-1
res@mpCenterLonF                = (lonL+lonR)/2

res@lbLabelBarOn        = True
res@pmLabelBarHeightF   = 0.08
res@pmLabelBarWidthF    = 0.55
res@lbLabelFontHeightF  = 0.018
res@pmLabelBarOrthogonalPosF    = 0.14

res@tmXBLabelFontHeightF        = 0.022    ;font height of tick labels
res@tmYLLabelFontHeightF        = 0.022
res@tmXBTickSpacingF            = 15.     ;label spacing
res@tmYLTickSpacingF            = 10.    ;label spacing

res@gsnLeftStringFontHeightF    = 0.022
res@gsnRightStringFontHeightF   = 0.022
res@gsnLeftString               = "(b) CTRL"
res@gsnRightString              = sprintf("%5.2f",x(0))+"%"
res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = -0.5
res@cnMaxLevelValF = 0.5
res@cnLevelSpacingF = 0.1

plot = gsn_csm_contour_map(wks,reg_sst,res)
print("end of plot")

; plot = gsn_csm_contour_map(wks,sst_anom(0,:,:),res)
; draw(plot)
; frame(wks)
   
vcres                 =True
vcres@gsnDraw         = False     
vcres@gsnFrame        = False
vcres@gsnAddCyclic    = False

vcres@vcRefMagnitudeF           = 0.5              ; define vector ref mag
vcres@vcRefLengthF              = 0.045          ; define length of vec ref
vcres@vcRefAnnoFontHeightF      = 0.02
vcres@vcRefAnnoOrthogonalPosF   = -1.34            ; move ref vector
vcres@vcRefAnnoArrowLineColor   = "black"         ; change ref vector color
vcres@vcRefAnnoArrowUseVecColor = False           ; don't use vec color for ref
vcres@vcRefAnnoString2     = ""

vcres@vcLineArrowThicknessF   = 2.0               ; change vector thickness
vcres@vcLineArrowColor        = "black"
      
vcres@tiMainString          = " "
vcres@gsnLeftString = ""
vcres@gsnRightString = ""

print("start of plot_vc")
plot_vc = gsn_csm_vector(wks, reg_uwnd_test(::1,::1), reg_vwnd_test(::1,::1),vcres)
print("end of plot_vc")

overlay(plot, plot_vc)
print("end of overlay")

rts           = True
rts@gsnDraw   = False       ; don't draw yet
rts@gsnFrame  = False       ; don't advance frame yet
rts@vpXF     = 0.15
rts@vpWidthF = 0.7
rts@vpHeightF= 0.35
rts@tiYAxisString = "Standardized"         ; y-axis label 
rts@tmXBLabelFontHeightF = 0.02    ;font height of tick labels
rts@tmYLLabelFontHeightF = 0.02    
rts@gsnLeftStringFontHeightF = 0.02
rts@gsnRightStringFontHeightF = 0.02

rts@gsnYRefLine        = 0.              ; reference line   
rts@xyLineColors       = (/"red","blue"/) ; colors chosen
rts@xyLineThicknesses  = (/3.0,3.0/)      ; line thicknesses
rts@xyDashPatterns     = (/0.,0./)      ; make all lines solid

rts@trXMinF          = yrStrt  ; leave a margin for legend
rts@trXMaxF          = yrLast+1 
rts@tmXBTickSpacingF = 15

rts@trYMinF      = -6.5            ; min value on x-axis
rts@trYMaxF      = 6.5            ; max value on x-axis

rts@gsnLeftString = "(a) Time Series of MCA1 (R = "+sprintf("%5.2f",ccr)+")" 
rts@gsnRightString = sprintf("%5.2f",x(0))+"%"
plot_ts = gsn_csm_xy(wks,years,data_std,rts)
print("end of plot_ts")

;---------------------------Add plot legend-----------------------------------
res_lines                   = True                  ; polyline mods desired
res_lines@gsLineDashPattern = 0.                    ; solid line
res_lines@gsLineThicknessF  = 3.                    ; line thicker

res_text                    = True                  ; text mods desired
res_text@txFontHeightF      = 0.015                 ; change text size
res_text@txJust             = "CenterLeft"          ; text justification

res_lines@gsLineColor       = "red"
yy = (/6.0,6.0/)
xx = (/278,283/)
text_pos = 285
dum11 = gsn_add_polyline(wks,plot_ts,xx,yy,res_lines)              ; add polyline
dum12 = gsn_add_text(wks,plot_ts,"PMMsst",text_pos,yy(0),res_text); add text

res_lines@gsLineColor       = "blue"
yy = yy - 0.8
dum21 = gsn_add_polyline(wks,plot_ts,xx,yy,res_lines)              ; add polyline
dum22 = gsn_add_text(wks,plot_ts,"PMMwind",text_pos,yy(0),res_text); add text

;************************************************
; panel plot only resources
;************************************************
resP                     = True         ; modify the panel plot
resP@gsnPanelLabelBar    = False         ; add common colorbar
resP@txString            = " "

gsn_panel(wks,(/plot,plot_ts/),(/2,1/),resP)    

end