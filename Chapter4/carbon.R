library(reticulate)#R使用python程式碼套件
py_install("codecarbon")#python計算炭排率套件


##########################生成50000個體資料
xr = 1000
yr = 500
N = 50000
r = (-(xr+yr)-sqrt((xr+yr)^2-4*(1-N)*xr*yr))/(2*(1-N))
df = data.frame(x = runif(N, 0, xr), y = runif(N, 0, yr))
##########################



########python code
py_run_string("
from codecarbon import OfflineEmissionsTracker
tracker = OfflineEmissionsTracker(country_iso_code='CAN')
tracker.start()
")
########


nni(meuse)$NNI



########python code
py_run_string("
tracker.stop()
")
########