;;;************************************************************
;;
;; linshuheng 2019.10.22
;;;*********************************************************
begin
  fils1 = "/home/ys17-19/lsh/CESM-data/F/F_2000_CTRL/CTRL.cam.h0.0101-4012.nc"
  fils2 = "/home/ys17-19/lsh/CESM-data/F/F_2000_addindianwarming/F_2000_indian.cam.h0.0101-4012.nc"  
  

    f1=addfile(fils1, "r")
    f2=addfile(fils2, "r")

    date=cd_calendar(f1->time, 1)
    timeind=ind(date.ge.601.and.date.le.4012)
    

    PSL_CT  =  f1->PSL(timeind,:,:)
    PSL_CT = (/PSL_CT/100/)

    PSL_CT = (/rmMonAnnCycTLL(PSL_CT)/)

    PSL_CT_enso = dim_avg_n_Wrap(PSL_CT(:,{-5:5},{160:280}),(/1,2/)) - dim_avg_n_Wrap(PSL_CT(:,{-5:5},{80:160}),(/1,2/))

    copy_VarCoords(PSL_CT(:,1,1), PSL_CT_enso)

    PSL_CT_enso_annual = month_to_annual(PSL_CT_enso, 1)
   
    PSL_CT_enso_annual = (/dim_standardize(PSL_CT_enso_annual, 1)/)

    printVarSummary(PSL_CT_enso_annual)
    


    res = True
    wks = gsn_open_wks("pdf", "PSL_CT_enso")
    plot = gsn_csm_xy(wks, ispan(1, dimsizes(PSL_CT_enso_annual), 1), PSL_CT_enso_annual, res)



  end
  

