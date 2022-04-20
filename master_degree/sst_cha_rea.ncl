;;;;读入hadisst 的海温
begin 

    print("*****************读入海温场*********************")
   ;;;海温场
    fils1=systemfunc("ls /mnt/f/share/ersst/haderscobe_sst_chazhi.nc")
    f1=addfile(fils1, "r")
    date=cd_calendar(f1->time, 1)
    yearstart=197901
    yearend=201212
    print("******************日期*********************")
    timeind=ind(date.le.yearend.and.date.ge.yearstart) 
    sst_mre_had=f1->sst(timeind,{-50:50},:)
    delete(sst_mre_had&time)
    printVarSummary(sst_mre_had)
    sst_mre_had&time = date(timeind)
    delete(date)
    delete(f1)
    delete(fils1)
    sst_mre_had_annual=month_to_annual(sst_mre_had, 1)
    sst_mre_had_annual!0="time"
    sst_mre_had_annual&time=ispan(1979,2012,1)
    printVarSummary(sst_mre_had_annual)

    sst_cha= dim_avg_n_Wrap(sst_mre_had_annual(:,{-5:5},{80:160}),(/1,2/)) - dim_avg_n_Wrap(sst_mre_had_annual(:,{-5:5},{200:280}),(/1,2/))
    sst_cha!0="time"
    sst_cha&time=ispan(1979,2012,1)
    printVarSummary(sst_cha)

    nmodel=1
    ntime=34

;***************************************************
;写入nc文件
;*************************************************
    diro="/mnt/f/share/ersst/"
    filo1="sstca_rea_r1i1p1_djf.nc"
    system("/bin/rm -f "+diro + filo1+"")
    fout=addfile(diro + filo1,"c")
    setfileoption(fout,"DefineMode",True)
   ;======================================文件名
    fAtt               = True            ; assign file attributes
    fAtt@title         = "location,msf djf平均后计算得到的位置结果"    
    fAtt@creation_date = "2018.4.26"
    fileattdef( fout, fAtt)            ; copys file attributes    
    ;----------------------------坐标变量的定义
    dimnames=(/"model","time"/)
    dimsize=(/nmodel,ntime/)
    dimunlim=(/False,False/)
    filedimdef(fout, dimnames, dimsize, dimunlim)
 ;-----------------------------定义变量名，类型和维数
    filevardef(fout, "model","integer",(/"model"/))
    filevardef(fout, "time", typeof(sst_cha&time),(/"time"/))

    filevardef(fout, "sst_cha", typeof(sst_cha),(/"model","time"/))

    
    filevarattdef(fout, "sst_cha", sst_cha)    
    setfileoption(fout,"DefineMode",False)
    fout->model=1
    fout->time=(/sst_cha&time/)
    fout->sst_cha=(/sst_cha/)


end






