
import numpy as np

def runningmedian(x,y,xlolim=-1.e20,ylolim=-1.e20,bins=10):
        xp = x[(x>xlolim)&(y>ylolim)]
        yp = y[(x>xlolim)&(y>ylolim)]
        hist,bin_edges=np.histogram(xp,bins)
        ymed = []
        ymean = []
        ysigma = []
        for i in range(0,len(bin_edges[:-1])):
                xsub = xp[xp>bin_edges[i]]
                ysub = yp[xp>bin_edges[i]]
                ysub = ysub[xsub<bin_edges[i+1]]
                ymed.append(np.median(10**ysub))
                ymean.append(np.mean(10**ysub))
                ysigma.append(np.std(10**ysub))
        bin_cent = 0.5*(bin_edges[1:]+bin_edges[:-1])
        ymean = np.asarray(ymed)
        ysiglo = np.maximum(ymean-ysigma,ymean*0.1)
        ysiglo = np.log10(ymean)-np.log10(ysiglo)
        ysighi = np.log10(ymean+ysigma)-np.log10(ymean)
        ymean = np.log10(ymean)
        #print bin_cent,ymean,ysiglo,ysighi
        #plt.plot(bin_cent,ymed,'ro',ms=12,color='c')
        #plt.plot(bin_cent,ymean,'--',lw=3,color='m')
        #plt.errorbar(bin_cent,ymean,yerr=[ysiglo,ysighi],fmt='ro')
	return bin_cent,ymean,ysiglo,ysighi

