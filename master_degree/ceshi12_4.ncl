begin

f1=addfile("/mnt/f/share/CMIP5/sst_cha/sstca_cmip_r1i1p1_djf.nc", "r")
sst_cha=f1->sst_cha
model=f1->model
time=f1->time
res=True
wks=gsn_open_wks("pdf","/mnt/f/share/era/data/cmip_2012_ssta_cha")
plot=gsn_csm_xy(wks, time, sst_cha(7,:), res)
end