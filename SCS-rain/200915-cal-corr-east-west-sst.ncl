;;**************************************************************
;;;main program
;;***************************************************************
begin

;;;读入HadiSST海温数据
    diri="/home/ys17-19/lsh/data/sst/"
    fils1=systemfunc("ls "+diri+"HadISST_sst.nc")
    f1=addfile(fils1, "r")
    date:=cd_calendar(f1->time, 1)

    timeind:=ind(date.le.201912.and.date.ge.197901)

    sst=lonFlip(f1->sst(timeind,:,:))
    printVarSummary(sst)

    ssta=(/rmMonAnnCycTLL(sst)/)
    copy_VarCoords(sst, ssta)

    ssta_jja = month_to_season(ssta, "JJA")
    ;;;
    delete(date)
    delete(timeind)
;;;


;;;;读入OISST 数据

fils2 = "/home/yangsong3/data-observation/xll/NOAA/sst.mnmean.nc"

f2 = addfile(fils2, "r")
date:=cd_calendar(f2->time, 1)
timeind:=ind(date.le.201912.and.date.ge.198201)

print(date(timeind))

sst_oi = short2flt(f2->sst(timeind,:,:))
printVarSummary(sst_oi)


ssta_oi = (/rmMonAnnCycTLL(sst_oi)/)
copy_VarCoords(sst_oi, ssta_oi)
ssta_oi_jja = month_to_season(ssta_oi, "JJA")
;;;


;;;区域平均
 ssta_west = dim_avg_n_Wrap(ssta_jja(:,{-8:15},{115:150}), (/1,2/))
 ssta_east = dim_avg_n_Wrap(ssta_jja(:,{-10:12},{170:250}), (/1,2/))

 ssta_oi_west = dim_avg_n_Wrap(ssta_oi_jja(:,{-8:15},{115:150}), (/1,2/))
 ssta_oi_east = dim_avg_n_Wrap(ssta_oi_jja(:,{-10:12},{170:250}), (/1,2/))
;;;

;;;相关系数

r_hadisst =  escorc_n(ssta_west, ssta_east, 0, 0)

r_oisst =  escorc_n(ssta_oi_west, ssta_oi_east, 0, 0)


print("r_hadisst=" + r_hadisst)
print("r_oisst=" + r_oisst)

end 