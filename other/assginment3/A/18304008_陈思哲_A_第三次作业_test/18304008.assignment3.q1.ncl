begin

    ; 读取 .nc 文件
    ftp          = addfile("tp.assignment3.nc","r")
    fuwnd        = addfile("uwnd.assignment3.nc", "r")
    fvwnd        = addfile("vwnd.assignment3.nc", "r")
    
    tp           = ftp->tp
    time_tp      = ftp->time
    uwnd         = fuwnd->uwnd
    time_uwnd    = fuwnd->time
    vwnd         = fvwnd->vwnd
    time_vwnd    = fvwnd->time

    
    ; ========================================================================
    ; 提取 .nc 文件指定的时间位置
    ; (建议使用 cd_calendar 及 ind)
    
    ; 2018年6月27日
    time_tp_ARR     = cd_calendar(time_tp, -2)
    time_tp_IND1    = ind(time_tp_ARR.eq.20180627)

    time_uwnd_ARR   = cd_calendar(time_uwnd, -2)
    time_uwnd_IND1  = ind(time_uwnd_ARR.eq.20180627)

    time_vwnd_ARR   = cd_calendar(time_vwnd, -2)
    time_vwnd_IND1  = ind(time_vwnd_ARR.eq.20180627)
    
    
    ; 2018年6月15日至7月15日
    time_tp_IND2    = ind(time_tp_ARR.ge.20180615 .and. time_tp_ARR .le. 20180715)

    time_uwnd_IND2  = ind(time_uwnd_ARR.ge.20180615 .and. time_uwnd_ARR .le. 20180715)

    time_vwnd_IND2  = ind(time_vwnd_ARR.ge.20180615 .and. time_vwnd_ARR .le. 20180715)
    
    ; ========================================================================
    ; 提取 2018年6月27日的降水数据
    ; 单位需要由 m 转为 mm
    tp_wanted1_m    = tp(time_tp_IND1,:,:)
    tp_wanted1_mm   = tp_wanted1_m*1000
    copy_VarCoords(tp_wanted1_m,tp_wanted1_mm)    ;copy_VarCoords可以复制前一个文件的各个维度的单位到新的文件
    tp_wanted1_mm@unit = "mm"                    ;改变单位为毫米
    
    
    ; 提取 2018年6月15日至7月15日北纬35度，东经127.5的降水数据
    ; 单位需要由 m 转为 mm
    tp_wanted2_m    = tp(time_tp_IND2,{35},{127.5})
    tp_wanted2_mm   = tp_wanted2_m*1000
    copy_VarCoords(tp_wanted2_m,tp_wanted2_mm) ;同上
    tp_wanted2_mm@unit = "mm"
    
    
    
    ; 提取 2018年6月27日 700百帕的风场
    uwnd_wanted     = uwnd(time_uwnd_IND1,{700},:,:)
    vwnd_wanted     = vwnd(time_vwnd_IND1,{700},:,:)
    
    ; ========================================================================
    ; 打开工作空间
    wks             = gsn_open_wks("eps", "18304008.assignment3.q1")
    
    
    ; =========================================================================
    ; 定义绘图属性
    
    res                                 = True
    resMark=res
    res@gsnDraw                         = False
    res@gsnFrame                        = False
    res@gsnLeftString                   = ""
    res@gsnRightString                  = ""
    
    resxy                               = res
    resxy2                              =res   ;用于xy图的叠加，在图上绘制横的参考线
    
    res@gsnAddCyclic                    = False
    rescn                               = res
    resvc                               = res
                
    res@vpWidthF                        = 0.7
                
    res@cnFillOn                        = True      
    res@cnLinesOn                       = False


    ; !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ; res: 图 (a) 底图
    
    ; 标题
    res@gsnCenterString                 ="(a) Total Precipitation (mm) & Wind Vector (m/s) on June 27,2018"
    
    
    ; 地图范围
    ; mpMinLatF, mpMaxLatF, mpMinLonF, mpMaxLonF
    res@mpMinLatF=15
    res@mpMaxLatF=45
    res@mpMinLonF=100
    res@mpMaxLonF=140
    
    
    ; 高分辨率地图
    ; mpDataBaseVersion
    res@mpDataBaseVersion               ="HighRes"
    res@mpDataResolution                ="CoarseResolution"   ;额外定义使用地图的精度为"CoarseResolution"
    res@mpGeophysicalLineThicknessF     =0.6 ;把地理的线的粗细调低，使图像美观
    
    
    ; 色板
    ; cnFillPalette
    res@cnFillPalette                   = "WhiteBlue"
    res@lbBoxEndCapStyle                = "TriangleHighEnd" ;因为降水的最低值就是0，而最高值可能更高，所以低的一端是方形而高的一端设置为三角形
    res@pmLabelBarHeightF               = 0.1
    res@pmLabelBarOrthogonalPosF        = 0.09
    
    
    ; 等值线间距
    ; cnLevelSelectionMode 及相应属性
    ; res@cnLevelSelectionMode="ManualLevels"
    ; res@cnMaxLevelValF=100
    ; res@cnMinLevelValF=0
    ; res@cnLevelSpacingF=5
    res@cnLevelSelectionMode            = "ExplicitLevels"            
    res@cnLevels                        = (/0,2.5,5,10,25,50,100,200/)       ;这一部分等值线间距的设定参照中央气象台所使用的的降水间距
    
    
    ; (可选) 字体大小
    ; XY轴：tmXBLabelFontHeightF, tmYLLabelFontHeightF
    ; 主标题: tiMainFontHeightF
    ; 副标题: gsnCenterStringFontHeightF
    res@tmXBLabelFontHeightF            =0.02
    res@tmYLLabelFontHeightF            =0.02
    res@gsnCenterStringFontHeightF      =0.02

    
    
    ; (可选) 主标题Y方向移动位置
    ; tiMainOffsetYF
    ;res@tiMainOffsetYF=0.003

    
    
    ; ---------------------------------------------------------------------------
    ; rescn: 图 (a) 等值线：降水=30mm

    ; 等值线间距
    ; cnLevelSelectionMode 及相应属性
    rescn@cnLevelSelectionMode          = "ExplicitLevels"
    rescn@cnLevels                      = 30               ;只有30的等值线做特别设定，稍后进行叠加
    
    ; 等值线颜色及粗细度
    ; cnLineColor, cnLineThicknessF
    rescn@cnLineColor                   = "Blue"
    rescn@cnLineThicknessF              = 2.
    
    ; 不显示等值线标签及信息
    ; cnLineLabelsOn, cnInfoLabelOn
    rescn@cnLineLabelsOn                = False
    rescn@cnInfoLabelOn                 = False

    
    
    ; ---------------------------------------------------------------------------
    ; resvc: 图 (a) 矢量
    
    ; 矢量样式
    ; vcGlyphStyle
    resvc@vcGlyphStyle                  ="CurlyVector"    
    
    
    ; (可选) 矢量的粗细度及最短距离(NDC单位)
    ; vcLineArrowThicknessF, vcMinDistanceF
    resvc@vcLineArrowThicknessF         =1.5
    resvc@vcMinDistanceF                =0.04
    
    

    ; 参考矢量的长度(NDC单位)及其表示的数值
    ; vcRefLengthF, vcRefMagnitudeF
    resvc@vcRefLengthF                  =0.055
    resvc@vcRefMagnitudeF               =20.
    ;resvc@vcMonoLineArrowColor=False
    ;resvc@vcLevelPalette="wgne15"
    
    
    ; 参考矢量信息的X,Y方向的位置
    ; vcRefAnnoParallelPosF, vcRefAnnoOrthogonalPosF
    resvc@vcRefAnnoParallelPosF         =0.5             ;把参考信息设置在右上角
    resvc@vcRefAnnoOrthogonalPosF       =-0.5
    resvc@vcRefAnnoZone                 =0
    
    ; (可选) 参考矢量信息的字体大小
    ;vcRefAnnoFontHeightF
    
            
    ; ---------------------------------------------------------------------------
    ; resxy: 图 (b) XY 图
    
    resxy@vpWidthF                      = 0.7
    resxy@vpHeightF                     = 0.55          ;这里修改了老师预先设定的图片高度，为了使展板两张图保持一致
    
    resxy@tmXTOn                        = False
    resxy@tmYROn                        = False
    
    resxy@trYMinF                       = 0.
    
    ; !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ; 以柱状图表示
    ; 标题
    resxy@tiMainString                  ="(b) Total Precipitation at 35~S~o~N~N,127.5~S~o~N~E During June 15-July 15,2018"
    ;resxy@tiMainOffsetYF=0.002
    ; X,Y轴标题
    ; tiXAxisString, tiYAxisString
    resxy@tiXAxisString                 ="time (date) "
    resxy@tiYAxisString                 ="Total Precipitation (mm)"
    ;resxy@tmXBLabelAngleF               =30.
    
    
    
    ; gsnXYBarChart
    resxy@gsnXYBarChart                 =True   ;使画成柱状图
    
    
    ; X 轴
    ; tmXBMode 及相应属性
    resxy@tmXBMode                      ="Explicit"
    resxy@tmXBValues                    =(/1,6,11,16,21,26/)    ;设置主刻度的位置，实际上x轴是定义了一个0到30的序列，只是再对应的位置加上序列值所代表的日期
    resxy@tmXBLabels                    =(/"June 16","June 21","June 26","July 1","July 6","July 11"/)
    resxy@tmXBMinorValues               =ispan(0, 30, 1)       ;设置副刻度，显示所有日子（也就是所有序列值）

    
    resxy2@gsnYRefLine                  = (/20,40,60,80/)      ;在叠加的第二张图上绘制横向的参考线（可能是只有折线图可以绘制这个，柱状图不行）
    resxy2@gsnYRefLineThicknessF        = 2.
    resxy2@gsnYRefLineColor             = "Gray"
    resxy2@gsnYRefLineDashPattern       = 11
    resxy2@trYMinF                      = 0.
    tp_y2                               =fspan(0, 0, 31)        ;将叠加的第二张图上的数值全部设置为0，这样叠加折线图的时候折线与x轴重合，就不影响原来的柱状图

    

    
    ; (可选) XY线及柱状的颜色
    resxy@gsnYRefLine                   =0
    resxy@gsnXYBarChartColors           ="Blue"
    resxy@gsnXYBarChartBarWidth         =0.8    ;将宽度设置为0.8，让每一个柱之间有一定间隔，方便观看
    resxy@gsnXYBarChartOutlineThicknessF=1.2    ;将每个柱的边线加粗
    
    ; (可选) 字体大小
    ; 主标题: tiMainFontHeightF
    ; X,Y轴标题: tiXAxisFontHeightF, tiYAxisFontHeightF
    ;resxy@gsnCenterStringFontHeightF=0.015
    resxy@tiXAxisFontHeightF            =0.02
    resxy@tiYAxisFontHeightF            =0.02
    resxy@tiMainFontHeightF             =0.02
    
    
    ; ---------------------------------------------------------------------------
    ; 标记 (北纬35度，东经127.5度)
    ; gsMarkerIndex, gsMarkerColor, gsMarkerSizeF
    ; (可选) gsMarkerThicknessF
    resmark                             = True
    resMark@gsMarkerIndex               = 15
    resMark@gsMarkerThicknessF          = 2.
    resMark@gsMarkerColor               = "Red"
    resMark@gsMarkerSizeF               = 0.015

    resText                     = True
    resText@txFontHeightF       = 0.02
    resText@txFontColor         = "Red" 
    text                        =(/"35~S~o~N~N,127.5~S~o~N~E"/)   ;在标记的上方加上表示标记的位置的注释

    
    ; =========================================================================
    ; 画图
    
    ; 图 (a)
    
    ; 降水地图
    plottp=gsn_csm_contour_map_ce(wks, tp_wanted1_mm, res)
    
    ; 降水=30 mm
    plotcn=gsn_csm_contour(wks,tp_wanted1_mm,rescn)
    
    ; 700百帕风场
    plotvc=gsn_csm_vector(wks, uwnd_wanted, vwnd_wanted, resvc)
    
    ; 添加标记
    x=127.5
    y=35         ;所需添加标记的位置
    plotMarker                  = gsn_add_polymarker(wks, plottp, x, y, resMark) ;添加标记，用之前标记的设定

    plotText                    = gsn_add_text(wks,plottp,text,x,y+1.6,resText)  ;在标记的上方一点的位置添加设定好的注释




    overlay(plottp, plotcn)
    overlay(plottp, plotvc)          ;将三张图叠加为一张
    ;draw(plottp)
    ;frame(wks)
    
    ; --------------------------------------------------------------------------
    ; 图 (b)
    tp_x    =ispan(0, 30, 1)                          ;如前所述，设定一个0到30的序列当做xy图的横坐标
    plotxy  =gsn_csm_xy(wks, tp_x, tp_wanted2_mm, resxy)
    plotxy2 =gsn_csm_xy(wks, tp_x, tp_y2, resxy2)     ;这里叠加的第二张xy图，所使用的数据tp_y2全部为0，只是为了叠加横向的参考线，模仿图a叠加30mm等值线的方法
    overlay(plotxy, plotxy2)                          ;叠加两张图
    ;draw(plotxy)
    ;frame(wks)
    
    
    ; =========================================================================
    ; 绘画展板
    ; gsnPanelXF 设定图(a)及(b)在X方向的位置(NDC坐标)
    ; gsnPanelYF 设定图(a)及(b)在Y方向的位置(NDC坐标)
    resPanel                            =True                    ;设定展板的相关设置
    resPanel@gsnPanelMainString         ="18304008 Chen Sizhe"
    resPanel@gsnPanelMainFontHeightF    =0.025
    resPanel@gsnPanelXWhiteSpacePercent =5.                      ;两张图之间设定一定的间隔

    plot                                =new(2,graphic)
    plot(0)                             =plottp
    plot(1)                             =plotxy
    gsn_panel(wks,plot,(/1,2/),resPanel)
    
    
    
    ; =========================================================================
    ; 去掉空白空间
    delete (wks)            ; 必须先关闭 workstation

    pltName="18304008.assignment3.q1"
    system("convert -geometry 3000x3000 -density 600x600 -trim " + pltName + ".eps " + pltName + ".png")   
    system("rm " + pltName + ".eps -f") 

    

end