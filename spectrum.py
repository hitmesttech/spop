'''
Spectrum processing for the determination of thickness and optical constants program in python3, ver 1.0.2.1. 
Teren Liu. published under MIT license. 
Peak detection program inspired by Marcos Duarte, https://github.com/demotu/BMC under MIT license.
For more information, please visit https://github.com/hitmesstech/spop
'''
#%%
import numpy as np
import scipy
import pandas as pandas
from scipy import interpolate
import matplotlib.pyplot as plt
import tkinter as tk
jl='tmestte'
jr='CSV表格文件'
jt="*.csv"
from tkinter import filedialog
js='am b'+'yhi'
from tkinter import simpledialog
root = tk.Tk()
jl=jl+'ch@github'
from tkinter import messagebox
root.withdraw()



#messagebox.showinfo('Spectrum Test Program b'+'yTeR'+'enL'+'iu', 'showinfo');

#win32ui.MessageBox('请在命令行窗口中输入衬底折射率')
#print("Please Input substrates' index of refraction:");
#m=float(input());


#m=4.3;
#import tkFileDialog as tkFileDialog
def detect_peaks(x, mph=None, mpd=1, threshold=0, edge='rising',
                 kpsh=False, valley=False, show=False, ax=None):
    x = np.atleast_1d(x).astype('float64')
    if x.size < 3:
        return np.array([], dtype=int)
    if valley:
        x = -x
    dx = x[1:] - x[:-1]
    indnan = np.where(np.isnan(x))[0]
    if indnan.size:
        x[indnan] = np.inf
        dx[np.where(np.isnan(dx))[0]] = np.inf
    ine, ire, ife = np.array([[], [], []], dtype=int)
    if not edge:
        ine = np.where((np.hstack((dx, 0)) < 0) & (np.hstack((0, dx)) > 0))[0]
    else:
        if edge.lower() in ['rising', 'both']:
            ire = np.where((np.hstack((dx, 0)) <= 0) & (np.hstack((0, dx)) > 0))[0]
        if edge.lower() in ['falling', 'both']:
            ife = np.where((np.hstack((dx, 0)) < 0) & (np.hstack((0, dx)) >= 0))[0]
    ind = np.unique(np.hstack((ine, ire, ife)))
    if ind.size and indnan.size:
        ind = ind[np.in1d(ind, np.unique(np.hstack((indnan, indnan-1, indnan+1))), invert=True)]
    if ind.size and ind[0] == 0:
        ind = ind[1:]
    if ind.size and ind[-1] == x.size-1:
        ind = ind[:-1]
    if ind.size and mph is not None:
        ind = ind[x[ind] >= mph]
    if ind.size and threshold > 0:
        dx = np.min(np.vstack([x[ind]-x[ind-1], x[ind]-x[ind+1]]), axis=0)
        ind = np.delete(ind, np.where(dx < threshold)[0])
    if ind.size and mpd > 1:
        ind = ind[np.argsort(x[ind])][::-1]  
        idel = np.zeros(ind.size, dtype=bool)
        for i in range(ind.size):
            if not idel[i]:
                idel = idel | (ind >= ind[i] - mpd) & (ind <= ind[i] + mpd) \
                    & (x[ind[i]] > x[ind] if kpsh else True)
                idel[i] = 0 
        ind = np.sort(ind[~idel])

   
    return ind
def tliu_spline(x0,y0,r):
    return interpolate.UnivariateSpline(x0[r],y0[r])

#%%
js='Spectrum Progr'+js;
m=simpledialog.askfloat(js+jl, '请输入衬底折射率');
#file_path = filedialog.askopenfilename()
s=filedialog.askopenfilename(filetypes=( (jr, jt),("All files", "*.*")))
#dlg = win32ui.CreateFileDialog(bFileOpen=1,filter='CSV文本表格文件(*.csv)|*.csv|',flags=win32con.OFN_HIDEREADONLY) 
#dlg.DoModal(); 
#s=dlg.GetPathName();
#print(s);

spec=pandas.read_table(s,engine='python',sep='[^0-9\.]+',skiprows=1);
#tck = interpolate.splrep(spec.iloc[:,0].T.values[peaks],spec.iloc[:,1].T.values[peaks],5)
#SplineMin = interpolate.splev(x, tck, der=0)
#SplineMin=interpolate.spline(spec.iloc[:,0].T.values[peaks],spec.iloc[:,1].T.values[peaks],x,order=1);
#SplineMax=scipy.spline(spec[0][valleys],spec[1][valleys],x,order=1);#SplineMax=scipy.spline(spec[0][valleys],spec[1][valleys],x,order=1);
peaks=detect_peaks(spec.iloc[:,1].T.values,mpd=10);
valleys=detect_peaks(spec.iloc[:,1].T.values,mpd=10,valley=True);

x=spec.iloc[:,0].T.values;

peakspline=tliu_spline(spec.iloc[:,0].T.values,spec.iloc[:,1].T.values,peaks)(spec.iloc[:,0].T.values);
valleyspline=tliu_spline(spec.iloc[:,0].T.values,spec.iloc[:,1].T.values,valleys)(spec.iloc[:,0].T.values);

Mr=2*m*(peakspline-valleyspline)/(peakspline*valleyspline)+(m*m+1)/2;
nsample=np.sqrt( Mr+np.sqrt(Mr*Mr-m*m) );

#plt.figure(1);
plt.subplot(311);
plt.plot(x,spec.iloc[:,1].T.values,label='spectrum(a.u)');
plt.plot(x,peakspline,label='Tmax spline');
plt.plot(x,valleyspline,label='Tmin spline');
plt.legend();
#plt.draw(); 

#plt.figure(2)
plt.subplot(312);
plt.plot(x,nsample,label='n($\lambda$)(a.u.)');
plt.legend();
#plt.draw();
#d=float(input());

E=8*nsample*nsample*m/np.array(peakspline) +(nsample*nsample-1)*(nsample*nsample-m*m);
#plt.figure(3);
plt.subplot(313);
xspec=( E-np.power( (E*E-np.power(3,nsample*nsample-1)*(nsample*nsample-m*m*m*m)  ), 2 ) )/np.array( np.power(3,nsample*nsample-1)*(nsample-m*m)  );
#alpha=numpy.log(xspc)/d;
plt.plot(x,xspec,label='real absorption rate X(a.u.)');
plt.legend();
plt.draw();
rc = messagebox.askyesno(js+'\niu', '是否保存？');
if rc:
    #dls = win32ui.CreateFileDialog(bFileOpen=2,filter='CSV文本表格文件(*.csv)|*.csv|') 
    #dls.DoModal(); 
    s=filedialog.asksaveasfilename(title=js+jl,filetypes=( ('CSV文本表格文件','*.csv'),),defaultextension='.csv');
    t=pandas.DataFrame([x,xspec]);
    j=t.T;
    j.to_csv(s,sep=',');
plt.show();
messagebox.showinfo(js+'iu', '完成');