;;计算气候态的风场数据，作为nudging的对象
;;
;;
;;                                         lnshuheng 2020.03.10
;**********************************************************
begin 

modelname = "B2000_alt_north_year_CTRL"
    diri="/home/ys17-19/lsh/CESM-data/B/"+ modelname+"/"
    fils1=systemfunc("ls "+diri+""+modelname+".cam.h0.0251-0280_UV.nc")
    f1=addfile(fils1, "r")
    
    U = f1->U
    V = f1->V
    lat = f1->lat
    lon = f1->lon
   

    UMC = clmMonTLLL(U)
    VMC = clmMonTLLL(V)

    UMC_DAY = clmMon2clmDay(UMC, 0, 0)
    VMC_DAY = clmMon2clmDay(VMC, 0, 0)
    
    printMinMax(UMC_DAY, 1)
    printVarSummary(UMC_DAY)
    printVarSummary(VMC_DAY)

    UMC_DAY!0 = "time"
    VMC_DAY!0 = "time"

    UMC_DAY&time = ispan(1, 365, 1)
    VMC_DAY&time = ispan(1, 365, 1)
    time = ispan(1, 365, 1)
    system("rm uv_nudging_to_WNPCLM.nc")
    fout = addfile("uv_nudging_to_WNPCLM.nc", "c")


    fout->lat = lat 
    fout->lon = lon
    fout->time = time
    fout->U = UMC_DAY    
    fout->V = VMC_DAY
end 