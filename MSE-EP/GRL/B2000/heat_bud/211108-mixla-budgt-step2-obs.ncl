load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
;==========================================================
begin
;========================================================== Basics settings
  diri = "/home/ys17-19/lsh/Project/MSE-EP/GRL/B2000/heat_bud/obs/"
  ; dim_models = dimsizes(model_name)

  opt = False

  his_start = 198001
  his_end   = 202012
;;;

;;;;clm 1982-2004 

;========================================================== Read data
 var_name = (/"TEMP","UVEL","VVEL","WVEL","SHF"/)

 i=0
;*******************************************
;  Read data
;*******************************************
    
    fin_thetao_his = addfile("./obs/"+var_name(0)+"_mixed_godas.nc", "r")


    thetao_his_mixed = fin_thetao_his->TEMP_mixed_avg
    printVarSummary(thetao_his_mixed)
    thetao_his_bot = fin_thetao_his->TEMP_mixed_bot
    printVarSummary(thetao_his_bot)

    thetao_his_mixed_CLM = clmMonTLL(thetao_his_mixed)
    printVarSummary(thetao_his_mixed_CLM)
    thetao_his_mixed_anom = calcMonAnomTLL(thetao_his_mixed, thetao_his_mixed_CLM)
    printVarSummary(thetao_his_mixed_anom)

    thetao_his_bot_CLM = clmMonTLL(thetao_his_bot)
    printVarSummary(thetao_his_bot_CLM)
    thetao_his_bot_anom = calcMonAnomTLL(thetao_his_bot, thetao_his_bot_CLM)
    printVarSummary(thetao_his_bot_anom)

;==========================================================
    fin_uo_his = addfile("./obs/"+var_name(1)+"_mixed_godas.nc", "r")
    uo_his_mixed = fin_uo_his->UVEL_mixed_avg
    printVarSummary(uo_his_mixed)

    uo_his_mixed_CLM = clmMonTLL(uo_his_mixed)
    printVarSummary(uo_his_mixed_CLM)
    uo_his_mixed_anom = calcMonAnomTLL(uo_his_mixed, uo_his_mixed_CLM)
    printVarSummary(uo_his_mixed_anom)

;==========================================================
    fin_vo_his = addfile("./obs/"+var_name(2)+"_mixed_godas.nc", "r")
    vo_his_mixed = fin_vo_his->VVEL_mixed_avg
    printVarSummary(vo_his_mixed)

    vo_his_mixed_CLM = clmMonTLL(vo_his_mixed)
    printVarSummary(vo_his_mixed_CLM)
    vo_his_mixed_anom = calcMonAnomTLL(vo_his_mixed, vo_his_mixed_CLM)
    printVarSummary(vo_his_mixed_anom)

;==========================================================
    fin_wo_his = addfile("./obs/"+var_name(3)+"_mixed_godas.nc", "r")
    wo_his_bot = fin_wo_his->WVEL_mixed_bot
    printVarSummary(wo_his_bot)

    wo_his_bot_CLM_org = clmMonTLL(wo_his_bot)
    wo_his_bot_CLM = where(wo_his_bot_CLM_org.le.0, 0, wo_his_bot_CLM_org)
    copy_VarCoords(wo_his_bot_CLM_org, wo_his_bot_CLM)
    printVarSummary(wo_his_bot_CLM)

    wo_his_bot_positive = where(wo_his_bot.le.0, 0, wo_his_bot)
    copy_VarCoords(wo_his_bot, wo_his_bot_positive)
    wo_his_bot_anom = calcMonAnomTLL(wo_his_bot_positive, wo_his_bot_CLM)
    printVarSummary(wo_his_bot_anom)

;==========================================================
    fin_hfds_his = addfile("./obs/"+var_name(4)+"_mixed_godas.nc", "r")
    hfds_his = fin_hfds_his->SHF_mixed_avg
    printVarSummary(hfds_his)

    hfds_his_CLM = clmMonTLL(hfds_his)
    printVarSummary(hfds_his_CLM)
    hfds_his_anom = calcMonAnomTLL(hfds_his, hfds_his_CLM)
    printVarSummary(hfds_his_anom)
;==========================================================

    dim_anom = dimsizes(hfds_his_anom)
    print(dim_anom)

;;;;calculate 
;==========================================================

    thetao_his_mixed_clm = thetao_his_mixed_anom
    do k = 0, 11
      do j = 0, dim_anom(0)/12-1
        thetao_his_mixed_clm(j*12+k,:,:) = (/thetao_his_mixed_CLM(k,:,:)/)
      end do
    end do
    printVarSummary(thetao_his_mixed_clm)

    thetao_his_bot_clm = thetao_his_bot_anom
    do k = 0, 11
      do j = 0, dim_anom(0)/12-1
        thetao_his_bot_clm(j*12+k,:,:) = (/thetao_his_bot_CLM(k,:,:)/)
      end do
    end do
    printVarSummary(thetao_his_bot_clm)

    uo_his_mixed_clm = uo_his_mixed_anom
    do k = 0, 11
      do j = 0, dim_anom(0)/12-1
        uo_his_mixed_clm(j*12+k,:,:) = (/uo_his_mixed_CLM(k,:,:)/)
      end do
    end do
    printVarSummary(uo_his_mixed_clm)

    vo_his_mixed_clm = vo_his_mixed_anom
    do k = 0, 11
      do j = 0, dim_anom(0)/12-1
        vo_his_mixed_clm(j*12+k,:,:) = (/vo_his_mixed_CLM(k,:,:)/)
      end do
    end do
    printVarSummary(vo_his_mixed_clm)

    wo_his_bot_clm = wo_his_bot_anom
    do k = 0, 11
      do j = 0, dim_anom(0)/12-1
        wo_his_bot_clm(j*12+k,:,:) = (/wo_his_bot_CLM(k,:,:)/)
      end do
    end do
    printVarSummary(wo_his_bot_clm)

;====================================================== dx & dy
    pi = get_pi("float")
    Re = 6371393

    dlon = 1./360. * 2 * pi   ; 格点间经度（横x）间距：0.5度
    dlat = 1./360. * 2 * pi   ; 格点间纬度（竖y）间距：0.5度

    lat_max = doubletofloat(max(thetao_his_mixed_anom&lat)) 

    dx = new((/dim_anom(0),dim_anom(1)-2,dim_anom(2)-2/), "float")
    do l = 1, dim_anom(1)-2
      lat_dx = (l-lat_max)
      print(lat_dx)
      dx(:,l-1,:) = dlon * Re * cos(lat_dx*pi/180)
    end do

    dy = dlat * Re

;========================================================== dTdx/y/z anom
    dTdxdy = grad_latlon_cfd(thetao_his_mixed_anom,thetao_his_mixed_anom&lat,thetao_his_mixed_anom&lon,True,False)
   
    dTdx_his_anom  = dTdxdy[1]
    dTdy_his_anom  = dTdxdy[0]

    dTdz_his_anom = (thetao_his_mixed_anom - thetao_his_bot_anom)/(60-25)
    copy_VarCoords(thetao_his_mixed_anom, dTdz_his_anom)
    printVarSummary(dTdz_his_anom)

;========================================================== dTdx/y/z clm

    dTclmdxdy = grad_latlon_cfd(thetao_his_mixed_clm,thetao_his_mixed_clm&lat,thetao_his_mixed_clm&lon,True,False)
   
    dTdx_his_clm  = dTclmdxdy[1]
    dTdy_his_clm  = dTclmdxdy[0]

    dTdz_his_clm =  (thetao_his_mixed_clm - thetao_his_bot_clm)/(60-25)
    copy_VarCoords(thetao_his_mixed_clm, dTdz_his_clm)
    printVarSummary(dTdz_his_clm)

;========================================================== dTdt
    dTdt_his = thetao_his_mixed_anom
    dTdt_his = (/center_finite_diff_n(thetao_his_mixed_anom,1,False, 1, 0)/)

;================================================================================= dynamics U
    Uanom_Txbar_his = uo_his_mixed_anom
    Uanom_Txbar_his = -uo_his_mixed_anom * dTdx_his_clm * 86400 * 30
    copy_VarCoords(uo_his_mixed_anom, Uanom_Txbar_his)
    printVarSummary(Uanom_Txbar_his)
    printMinMax(Uanom_Txbar_his, opt)

    Ubar_Txanom_his = uo_his_mixed_clm
    Ubar_Txanom_his = -uo_his_mixed_clm * dTdx_his_anom * 86400 * 30
    copy_VarCoords(uo_his_mixed_anom, Ubar_Txanom_his)
    printVarSummary(Ubar_Txanom_his)
    printMinMax(Ubar_Txanom_his, opt)

    Uanom_Txanom_his = -uo_his_mixed_anom * dTdx_his_anom  * 86400 * 30
    copy_VarCoords(uo_his_mixed_anom, Uanom_Txanom_his)
    printVarSummary(Uanom_Txanom_his)
    printMinMax(Uanom_Txanom_his, opt)

;================================================================================= dynamics V
    Vanom_Tybar_his = vo_his_mixed_anom
    Vanom_Tybar_his = -vo_his_mixed_anom * dTdy_his_clm * 86400 * 30
    copy_VarCoords(vo_his_mixed_anom, Vanom_Tybar_his)
    printVarSummary(Vanom_Tybar_his)
    printMinMax(Vanom_Tybar_his, opt)

    Vbar_Tyanom_his = vo_his_mixed_clm
    Vbar_Tyanom_his = -vo_his_mixed_clm * dTdy_his_anom * 86400 * 30
    copy_VarCoords(vo_his_mixed_anom, Vbar_Tyanom_his)
    printVarSummary(Vbar_Tyanom_his)
    printMinMax(Vbar_Tyanom_his, opt)

    Vanom_Tyanom_his = -vo_his_mixed_anom * dTdy_his_anom  * 86400 * 30
    copy_VarCoords(vo_his_mixed_anom, Vanom_Tyanom_his)
    printVarSummary(Vanom_Tyanom_his)
    printMinMax(Vanom_Tyanom_his, opt)

;================================================================================= dynamics W
    Wanom_Tzbar_his = wo_his_bot_anom
    Wanom_Tzbar_his = -wo_his_bot_anom * dTdz_his_clm * 86400 * 30
    copy_VarCoords(wo_his_bot_anom, Wanom_Tzbar_his)
    printVarSummary(Wanom_Tzbar_his)
    printMinMax(Wanom_Tzbar_his, opt)

    Wbar_Tzanom_his = wo_his_bot_clm
    Wbar_Tzanom_his = -wo_his_bot_clm * dTdz_his_anom * 86400 * 30
    copy_VarCoords(wo_his_bot_anom, Wbar_Tzanom_his)
    printVarSummary(Wbar_Tzanom_his)
    printMinMax(Wbar_Tzanom_his, opt)

    Wanom_Tzanom_his = -wo_his_bot_anom * dTdz_his_anom  * 86400 * 30
    copy_VarCoords(wo_his_bot_anom, Wanom_Tzanom_his)
    printVarSummary(Wanom_Tzanom_his)
    printMinMax(Wanom_Tzanom_his, opt)
;================================================================================= Nonlinear terms

    Uanom_Txbar_his = uo_his_mixed_anom
    Uanom_Txbar_his = -uo_his_mixed_anom * dTdx_his_clm * 86400 * 30
    copy_VarCoords(uo_his_mixed_anom, Uanom_Txbar_his)
    printVarSummary(Uanom_Txbar_his)
    printMinMax(Uanom_Txbar_his, opt)

    
    Wanom_Tzanom_his = wo_his_bot_anom
    Wanom_Tzanom_his = -wo_his_bot_anom * dTdz_his_anom * 86400 * 30
    copy_VarCoords(wo_his_bot_anom, Wanom_Tzanom_his)


   
;================================================================================= thermodynamics Q
    Hm = 50
    Cp = 3940;1996
    rho0 = 1022.4;1026
    Q_his = hfds_his_anom/(Cp*rho0*Hm)  * 86400 * 30
    copy_VarCoords(hfds_his_anom, Q_his)
    printVarSummary(Q_his)
    printMinMax(Q_his, False)

;================================================================================= output
    UATXB = Uanom_Txbar_his
    VATYB = Vanom_Tybar_his
    WATZB = Wanom_Tzbar_his

    UBTXA = Ubar_Txanom_his
    VBTYA = Vbar_Tyanom_his
    WBTZA = Wbar_Tzanom_his

    UATXA = Uanom_Txanom_his
    VATYA = Vanom_Tyanom_his
    WATZA = Wanom_Tzanom_his

    Q = Q_his


    ; UATXB = (/UATXB/100.0/)  ;;;centimeter to meter 
    ; VATYB = (/VATYB/100.0/)

    ; WATZB = (/WATZB/100.0/)

    ; UBTXA = (/UBTXA/100.0/)
    ; VBTYA = (/VBTYA/100.0/)
    ; WBTZA = (/WBTZA/100.0/)

    ; UATXA = (/UATXA/100.0/)
    ; VATYA = (/VATYA/100.0/)

    ; WATZA = (/WATZA/100.0/)

;=================================================================================
    system("rm ./obs/"+"hb_godas.nc")
    out = addfile("./obs/"+"hb_godas.nc", "c")
    out->dTdt   = dTdt_his
    out->UATXB  = UATXB
    out->VATYB  = VATYB
    out->WATZB  = WATZB

    out->UBTXA  = UBTXA
    out->VBTYA  = VBTYA
    out->WBTZA  = WBTZA

    out->UATXA  = UATXA
    out->VATYA  = VATYA
    out->WATZA  = WATZA

    out->Q = Q
  ; end do

end
