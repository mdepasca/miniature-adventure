description="Perform SN light curve fit"
""" #########################################################################
2013/09/13 Start
""" #########################################################################
import os,sys,shutil
import time
from numpy import *
from scipy import stats
import StringIO
import argparse
from pyraf import iraf
import pylab
import pickle
import pymysql as mdb

sys.path.append('/home/enrico/scripts/alice')
import alice
sys.path.append('/home/enrico/scripts/python_scripts')
import cosmo

xx = iraf.stsdas(_doprint=0,Stdout=1)
iraf.hst_calib(_doprint=0)
iraf.synphot(_doprint=0)

asa_dir = '/home/supern/asa/'

if __name__ == "__main__":
    start_time = time.time()
    parser = argparse.ArgumentParser(description=description,\
         formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("candidate",help='candidate id')
    parser.add_argument("-t","--template",help='template sn name',default='all')
    parser.add_argument("-g", "--guess",dest="guess",action="store_true",\
           default=False,help='guess mode (coarse steps)')
    parser.add_argument("-z","--zrange",help='redshift range',\
                            default='0.01,0.8')
    parser.add_argument("-p","--phrange",help='phase range',\
                            default='m25,25')
    parser.add_argument("-a","--avrange",help='AV range',\
                            default='0.0,2.0')
    parser.add_argument("-v", "--verbose",dest="verbose",action="store_true",\
           default=False,help='Enable task progress report')
    parser.add_argument("-d", "--debug",dest="debug",action="store_true",\
           default=False,help='Enable debug mode')
    parser.add_argument("-c", "--chi",dest="chi",action="store_true",\
           default=False,help='Show chisquare map')

    args = parser.parse_args()

zra_kk,zstep_kk = [0.0,1.0],0.01
#kcor_range = {'g': ('UB',  [1.00, 0.19]), \
#              'r': ('UBV', [1.00, 0.55, 0.24]), \
#              'i': ('UBVR',[1.00, 0.88, 0.55, 0.26])} 

uerr = 1.                           # error adopted for upper limits
kcor_range = {}  
kcor_range['BVR'] = {'g': ('B',  [1.00]), \
              'r': ('BV', [1.00, 0.24]), \
              'i': ('BVR',[1.00, 0.51, 0.24])} 
kcor_range['gri'] = {'g': ('g',  [1.00]), \
              'r': ('gr', [1.00, 0.13]), \
              'i': ('gri',[1.00, 0.38, 0.10])} 

class assign_filter:
    def __init__(self,filtro,zed,tbands):
        if all(array([x  in tbands for x in 'BVR'])): 
            self._kcor_range = kcor_range['BVR']
        elif all(array([x  in tbands for x in 'gri'])):
            self._kcor_range = kcor_range['gri']
        iz = -len([x for x in self._kcor_range[filtro][1] if x<zed])-1
        while self._kcor_range[filtro][0][iz] not in tbands: iz += 1
        self.iz = iz
        self.filtro = filtro

    def f(self):
        return self._kcor_range[self.filtro][0][self.iz]

def mdb_conn():           ####################   connect to mysql database  
    conn = mdb.connect('passto.oapd.inaf.it', 'sudare', 'pignatubo','sudare')  
    return conn

def DictCursor(cursor,righe,one):

    field_name = [i[0] for i in cursor.description]
    if one == 1: righe = [righe]

    output = []
    for r in  righe:
        rr = {}
        for i,f in enumerate(field_name):
            if isinstance(r[i],unicode): rr[f] = r[i].encode('ascii') 
            else: rr[f] = r[i]
        output.append(rr)
    if one == 1: output = output[0]        

    return output

def read_lc_mysql(trid,verbose):

    try:
        conn = mdb_conn()
        cursor = conn.cursor()     
        cursor.execute('SELECT lightcurve from master WHERE trid=%s',(trid))
        _riga = cursor.fetchone()
    except mdb.Error,e:
        return "Error %d: %s" % (e.args[0],e.args[1])

    if not _riga:
        print "!!! ERROR: candidate",trid,"not in master database !!!"
        sys.exit()

    riga = DictCursor(cursor,_riga,True)

    righe = riga['lightcurve'].split('\n')
    idx = {}
    for i in range(len(righe)):
        if len(righe[i].strip()):
            if '#' ==righe[i][0]: idx[righe[i][1]] = i

    jd,mag,merr,snr = {},{},{},{}
    for i,f in enumerate('gri'):
        jd[f],mag[f],merr[f],snr[f] = [],[],[],[]
        if f == 'g':  _righe = righe[idx['g']+1:idx['r']]
        elif f == 'r':_righe = righe[idx['r']+1:idx['i']]
        elif f == 'i':_righe = righe[idx['i']+1:]        

        for r in _righe:
            if len(r.strip()):
                _snr = float(r.split()[8])
                magt = float(r.split()[9])
                jd[f].append(float(r.split()[1]))
                snr[f].append(_snr)
                if r.split()[2]=='INDEF' or r.split()[3]=='INDEF':
                    mag[f].append(magt)
                    merr[f].append(-1.)
                elif r.split()[2]=='REF':
                    mag[f].append(magt)
                    merr[f].append(-2.)
                else:
                    if _snr<2:
                        mag[f].append(magt)
                        merr[f].append(-1.)
                    else:
                        mag[f].append(float(r.split()[2]))
                        merr[f].append(float(r.split()[3]))
 
        jd[f],mag[f],merr[f],snr[f] = array(jd[f]),array(mag[f]),\
            array(merr[f]),array(snr[f])

    if verbose:
        print ">>> read lc   jd=",min(jd['r']),"-",max(jd['r']),\
            "filtro r ",min(mag['r']),"-",max(mag['r'])
    return jd,mag,merr,snr

def read_template_lc(sndata,tsn):  #  read candidate light curves

    if all(array([x  in sndata['bands'] for x in 'BVR'])): _bands,b = 'BVR','B'
    elif all(array([x  in sndata['bands'] for x in 'gri'])):_bands,b = 'gri','g'

    jdmax = sndata['jdmax'][b]
    Rv = 3.1
    if sndata['ABi'] > 99:
        print '!!! WARNING: for SN',_sn,'ABi=',str(sndata['ABi']),\
                   ' (not available) ==> set to 0'
        sndata['ABi'] = 0.0

    ph,absmag = {},{}

    for b in _bands:
        if b not in sndata['bands']:
            print "!!! ERROR: band",b,"not available for SN",sndata['sn']
            sys.exit()
        abx = alice.AX(b,sndata['ABg']+sndata['ABi'],R_V=Rv)
        jd,mag,source = array(sndata['jd'][b]),array(sndata['mag'][b]),\
            array(sndata['source'][b])
        __ph = jd-jdmax

        ii =where(source>=0)
        _ph,_absmag = __ph[ii],mag[ii]-tsn['mu']-abx
        jj = argsort(_ph)
        ph[b],absmag[b] = _ph[jj],_absmag[jj]

    EBV = (sndata['ABg']+sndata['ABi'])/4.1

    return ph,absmag,jdmax,sndata['bands'],EBV

def read_template_spectra(tsnlist,tsn):

    tspe = {}
    for sn in tsnlist:
        tspe[sn] = {}
        for c in ['spec','ph','wmin','wmax','flag']:
            tspe[sn][c] = []

    ff = open('sn_template_spec.csv')
    righe = ff.readlines()
    for r in righe[4:]:
        sn = r.split(',')[1]
        if sn in tsnlist:
            ph = float(r.split(',')[2])-2400000-tsn[sn]['jdmax']

            tspe[sn]['spec'].append(r.split(',')[0])
            tspe[sn]['ph'].append(ph)
            for i,c in enumerate(['wmin','wmax','flag']):
                tspe[sn][c].append(float(r.split(',')[i+3].strip('\n')))

    return tspe


def kmm(i,filt_in,filt_ou,zrange,zstep,z0,wmin,wmax,verbose):
 
    err = StringIO.StringIO()
    if verbose: stderr = 0
    else: stderr = err
    pref = ''
    lamc_in,fwhm_in = alice.bandpar[filt_in][2],alice.bandpar[filt_in][3]
    lamc_ou,fwhm_ou = alice.bandpar[filt_ou][2],alice.bandpar[filt_ou][3]
    if wmin/(1+z0) < lamc_in-fwhm_in/2. and wmax/(1+z0) > lamc_in+fwhm_in/2.: 
        if filt_in in 'ugriz': pref = 'sdss,'
        _out =  iraf.calcphot(pref+filt_in,'z(/tmp/_tmp.'+str(i)+","+\
                       str(-z0)+')','vegamag',Stdout=1,Stderr=stderr)
    else: return ''

    mag_in = float(_out[6].split()[1])
    if mag_in == 100: return ''

    pref = ''
    if filt_ou in 'ugriz': pref='sdss,'
    if zstep == 0: vzero=str(zrange[0]-z0)
    else: vzero = str(zrange[0]-z0)+'-'+str(zrange[1]-z0)+'x'+str(zstep)
    _out =  iraf.calcphot(pref+filt_ou,'z(/tmp/_tmp.'+str(i)+',$0)','vegamag',\
          vzero=vzero,Stdout=1,Stderr=stderr)

    kcor = []
    for i,m in enumerate(_out[6:]):
        if  wmin*(1+(zrange[0]-z0+i*zstep)) < lamc_ou-fwhm_ou/2. and \
            wmax*(1+zrange[0]-z0+i*zstep) > lamc_ou+fwhm_ou/2. and \
            float(m.split()[1]) != 100: 
            kcor.append(mag_in-float(m.split()[1])) 
        else: kcor.append(100)
 
    return kcor

def kcor_compute(sn,zrange,zstep,tspe,tsn,verbose): #########

    lam = arange(3000,10000,1)
    ablam = array([alice.cardelli(l,3.1) for l in lam])

    if all(array([x  in tsn['bands'] for x in 'BVR'])): 
        _kcor_range = kcor_range['BVR']
    elif all(array([x  in tsn['bands'] for x in 'gri'])):
        _kcor_range = kcor_range['gri']

    kph,kz,kcor = {},{},{}
    for f in 'gri':
        for g in _kcor_range[f][0]:
            kph[g+f],kcor[g+f] =[],[]

    for i,spe in enumerate(tspe['spec']):
        spe = spe[:-5]
        if verbose: print spe,
        iraf.delete('/tmp/_tmp?'+str(i),verify=False)
        iraf.wspec(asa_dir+sn+"/"+spe+'[*,1,1]',"/tmp/_tmp."+str(i),\
                       header=False)

        ff = open("/tmp/_tmp."+str(i))
        righe = ff.readlines()
        ll,fl = [],[]
        for r in righe: 
            if tspe['wmin'][i]<float(r.split()[0])<tspe['wmax'][i]: 
                ll.append(float(r.split()[0]))
                fl.append(float(r.split()[1]))

        if tsn['EBV']:  
            ll,_fl = array(ll),array(fl)
            iablam = interp(ll,lam,ablam)
            fl = _fl*10**(0.4*iablam*tsn['EBV']*3.1)

        ff = open('/tmp/_tmp.'+str(i),'w')
        for l in range(len(ll)):
            ff.write(str(ll[l])+' '+str(fl[l])+' \n')
        ff.close()

        for f in 'gri':
            z0 = 0
            for k in range(len(_kcor_range[f][1]))[::-1]:
                g = _kcor_range[f][0][k]
                z1 = _kcor_range[f][1][k]
                zra = [z0,z1]
                if verbose: print '-',f+g,'-',
                z0 = z1
                kmmout = kmm(i,g,f,zra,zstep,\
                    0.0,tspe['wmin'][i],tspe['wmax'][i],False)
   # to be consistent with photometry no correction for SN template redshift
                if kmmout:
                    if any(array(kmmout)<100):
                        kcor[g+f].append(kmmout)
                        kph[g+f].append(tspe['ph'][i])
                        kz[g+f] = arange(zra[0],zra[1]+zstep/2,zstep)    

        if verbose: print

    for gf in kph.keys():     
        kph[gf],kcor[gf] = array(kph[gf]),array(kcor[gf])

    for gf in kph.keys():                    #     clean output
        for i in range(len(kz[gf]))[::-1]:
            if all(array(kcor[gf])[:,i]==100):
                kz[gf] = delete(kz[gf],i)
                kcor[gf] = delete(kcor[gf],i,1)
                
    return kph,kz,kcor

def kcor_interpolate(tph,z,kph,kz,kcor):

    zz = argmin(abs(kz-z))
    ii = where(kcor[:,zz]<100)
    _kcor = interp(tph,kph[ii],kcor[:,zz][ii])
    return _kcor

def mag_observer_frame(tph,tabsmag,tzed,z,kcor,tEXT):

    phobs = tph*(1-tzed+float(z))
    magobs = tabsmag-kcor+cosmo.mu(float(z))+tEXT

    return phobs,magobs

# --------------------------------------------------------------------
# Best fit
# --------------------------------------------------------------------
def best_fit(sn, t, jd, 
             mag, merr, snr, 
             tph, tabsmag, tzed, tbands, 
             kph, kz, kcor, 
             _zrange, _phrange, _avrange):
    # Guess max to start: it's the JD with the lowest 
    # magnitude in 'r' band (magnitude is an inverse scale)
    jdmax = jd['r'][argmin(mag['r'])] 

    # arange([start,]stop[,step,][,dtype])
    zrange  = arange(_zrange[0], _zrange[1]+_zrange[2]/2., _zrange[2])      # zed range 
    phrange = arange(_phrange[0], _phrange[1]+_phrange[2]/2., _phrange[2]) # phase range
    avrange = arange(_avrange[0], _avrange[1]+_avrange[2]/2., _avrange[2]) # extiction (reddening) range

    chi = zeros((len(phrange),len(zrange),len(avrange))) # chi is a 3-Dimensional array
    """
    enumerate(sequence, start=0) generates a tuple containing a count and the values 
    obtained from iterating over sequence:
    [(0,sequence[0]),(1,sequence[1]),...,(len(sequence)-1,sequence[len(sequence)-1])].
    Thus, we need two indexes in each of the 'for loop'
    """
    for i, phmax in enumerate(phrange) :
        for j, zed in enumerate(zrange):
            for k, AV in enumerate(avrange):
                _chi = 0
                np = 0.0        # number of points used in fitting and chi calculation see below (selection on error and SNR)
                # scanning each photometric band
                for f in 'gri':
                    g = assign_filter(f, zed, tbands).f() # k correction range
                    _kcor = kcor_interpolate(tph[g], zed,kph[g+f], kz[g+f], kcor[g+f])

                    avfact = alice.cardelli(alice.bandpar[f][2]/(1+zed),3.1)

                    phobs, magobs = mag_observer_frame(tph[g], tabsmag[g], tzed, zed, _kcor, avfact*AV)

                    ph = jd[f] - jdmax + phmax

                    # compute chisquare
                    ii = where((merr[f] > 0) & (snr[f] >= 3)) # SNR >= 3
                    magfit = interp(ph[ii], phobs, magobs, left=25, right=25) # Numpy function. left and right are the values to return if ph[ii] (x) is outside phobs (xp)
                    _chif = sum(((mag[f][ii]-magfit)/merr[f][ii])**2)

                    hh = where((merr[f] > 0) & (snr[f] < 3)) # SNR < 3 #snr<3 merr*(4-snr)
                    magfit = interp(ph[hh],phobs,magobs,left=25,right=25)
                    _chif += sum(((mag[f][hh]-magfit)/((4-snr[f][hh])*merr[f][hh])**2))

                    # print f,_chif,len(ii[0])

                    # Deling with UPPER LIMITS
                    jj = where(merr[f] < 0) # magnitude error < 0
                    kk = [[]]
                    if len(jj[0]):    # add upper limits
                        uph,umag = ph[jj],mag[f][jj]
                        umagfit = interp(uph,phobs,magobs,left=25,right=25)
                        kk = where(umag-umagfit>0)
                        _uchif = sum(((umag[kk]-umagfit[kk])/uerr)**2)

                        _chif+=_uchif
                    # print f,_uchif,len(kk[0])
                    _chi += _chif
                    np += len(ii[0])+len(hh[0])+len(kk[0]) # number of points used
                chi[i,j,k] = _chi/(np-3)

    imin,jmin,kmin = unravel_index(argmin(chi),chi.shape)

###  compute  confidence levels ........
    chimin = chi[imin,jmin,kmin]
    dchi90 = stats.f.ppf(0.90,2,np-3)
    dchi95 = stats.f.ppf(0.95,2,np-3)
    dchi99 = stats.f.ppf(0.99,2,np-3)
    tchi = []
    for i in range(len(phrange)):
        tchi.append(chi[i,:,:].min())
    ii = where((array(tchi)/chimin-1)*(np-3)<dchi95)
    phlo,phup= phrange[ii][0],phrange[ii][-1]

    tchi = []
    for i in range(len(zrange)):
        tchi.append(chi[:,i,:].min())
    ii = where((array(tchi)/chimin-1)*(np-3)<dchi95)
    zlo,zup= zrange[ii][0],zrange[ii][-1]

    tchi = []
    for i in range(len(avrange)):
        tchi.append(chi[:,:,i].min())
    ii = where((array(tchi)/chimin-1)*(np-3)<dchi95)
    avlo,avup= avrange[ii][0],avrange[ii][-1]

##########
    phflag,zflag,avflag = '','',''
    if len(phrange)>1:
        if imin==0: phflag='ph-'
        if imin==len(phrange)-1: phflag='ph+'
    if len(zrange)>1:
        if jmin==0: zflag='z-'
        if jmin==len(zrange)-1: zflag='z+'
    if len(avrange)>1:
        if kmin==len(avrange)-1: avflag='AV+'

    print "%5s %6s   %4.2f    %7.1f+%4.1f    %4.2f    %4.2f   %5s %5s %5s "%\
        (t,sn,zrange[jmin],jdmax,phrange[imin],avrange[kmin],chimin,phflag,\
             zflag,avflag)

    if args.chi:
        fig = pylab.figure()
        ax = pylab.axes([.15,.1,.8,.85])
        pylab.imshow(chi[:,:,kmin],vmin=0,vmax=100)

        circ = pylab.Circle((jmin,imin),radius=.3,fill=True,color='w')
        ax.add_patch(circ)
        CS  = pylab.contour(chi[:,:,kmin],[(chimin+1)*dchi99/(np+3.),\
                  (chimin+1)*dchi95/(np+3.),(chimin+1)*dchi90/(np+3)])     
        pylab.clabel(CS,inline=1,fontsize=12)
        xti = pylab.getp(ax,'xticks')
        yti = pylab.getp(ax,'yticks')
        pylab.yticks(yti[1:-1],array(yti[1:-1])*_phrange[2]+_phrange[0])
        pylab.xticks(xti[1:-1],array(xti[1:-1])*_zrange[2]+_zrange[0])
        pylab.ylabel('phase')
        pylab.xlabel('redshift')

        ax2 = pylab.axes([.75,.7,.2,.2],frameon=False)
        ax2.plot(avrange,chi[imin,jmin,:],'-w')
        ax2.spines['bottom'].set_color('w')
        ax2.set_xlabel(r'A$_V$',color='w')
        ax2.set_yticks([])
        xti = pylab.getp(ax2,'xticks')
        ax2.set_xticks(arange(xti[0],xti[-1],.5))
        ax2.tick_params(axis='x', colors='w')
        
    return [phrange[imin],phlo,phup],[zrange[jmin],zlo,zup],\
        [avrange[kmin],avlo,avup],chimin,np

def plot_lc(jd,mag,merr,snr,sn,tph,tabsmag,tsn,kph,kz,kcor,\
                phmax,z,AV,chi):

    pylab.rcParams['font.size']=16
    pylab.rcParams['legend.fontsize']=8
    pylab.figure()
    avfact = {}
    for g in 'UBVRgri':
        avfact[g] = alice.cardelli(alice.bandpar[g][2],3.1)

    g = {}
    for f in 'rgi': g[f] = ''
    
    jdmax = jd['r'][argmin(mag['r'])]

    for f in 'rgi':
        if not sn: break
        g[f] = assign_filter(f,z,tsn['bands']).f()

        _kcor = kcor_interpolate(tph[g[f]],float(z),kph[g[f]+f],\
                               kz[g[f]+f],kcor[g[f]+f])

        if   f=='r': pylab.subplot(211)
        elif f=='g': pylab.subplot(223)
        elif f=='i': pylab.subplot(224)

        phfit,magfit = mag_observer_frame(tph[g[f]],tabsmag[g[f]],\
                  tsn['zed'],z,_kcor,avfact[g[f]]*AV)
        pylab.plot(phfit,magfit,'o',\
              label="%s %s z=%4.2f max=%7.1f+%4.1f AV=%3.1f chi=%5.2f "% \
              (sn,tsn['type'],z,jdmax,phmax,AV,chi),mec='r',mfc='w',mew=1)
        if f =='r': pylab.legend(numpoints=1)

    pylab.subplots_adjust(hspace=0.15)
    pylab.subplot(211)
    pylab.title(args.candidate,fontsize=16)
    xmin,xmax = 100,-100
    ymin,ymax = 0,30
    for f in 'rgi':
        if f=='g': pylab.subplot(223)
        elif f=='i': pylab.subplot(224)
        ph = jd[f]-jdmax+phmax
        ii = (merr[f]>0)&(snr[f]>=3)
        if True in ii.tolist():
            pylab.errorbar(ph[ii],mag[f][ii],yerr=merr[f][ii],fmt='o',mfc='b',\
               mec='b',ecolor='b')
        ii = (merr[f]>0)&(snr[f]<3)
        if True in ii.tolist():
            pylab.errorbar(ph[ii],mag[f][ii],yerr=merr[f][ii],fmt='o',mfc='w',\
               mec='b',ecolor='b')
        if False in ii.tolist():
            ii = merr[f]==-1
            if True in ii.tolist():
                yerru = zeros(len(ph[ii]))
                yerrl = yerru-.5
                pylab.errorbar(ph[ii],mag[f][ii],yerr=[yerrl,yerru],\
                  ecolor='b',lolims=True,fmt='o',mfc='b',mec='b') 
            ii = merr[f]==-2
            if True in ii.tolist():
                yerru = zeros(len(ph[ii]))
                yerrl = yerru-.5
                pylab.errorbar(ph[ii],mag[f][ii],yerr=[yerrl,yerru],\
                  ecolor='g',lolims=True,fmt='o',mfc='b',mec='b') 

        if len(ph):
            xmin,xmax = min(xmin,min(ph)-20),max(xmax,max(ph)+10)
            ymin,ymax = max(mag[f])+.5,min(mag[f])-.5
            if ymin-ymax<3: ymin,ymax = ymin+1,ymax-1
            pylab.ylim(ymin,ymax)
    
        pylab.text(min(ph)-20,max(mag[f])+.3,f+' ['+g[f]+']')
        pylab.yticks(range(int(ymin),int(ymax),-1))

    for ii in [211,223,224] :
        pylab.subplot(ii)
        pylab.xlim(xmin,xmax)
        pylab.xticks(range(int((xmin-20)/50)*50,int((xmax+10)/50)*50,50))
         
    pylab.subplot(223)
    pylab.xlabel('       phase [observer frame]') 


def plot_kkcor(sn,kz,kph,kcor,tsn,tph,tabsmag):
    pylab.ion()

    avfact = {}
    pylab.ion()
    for g in 'UBVR':
        avfact[g] = alice.cardelli(alice.bandpar[g][2],3.1)

    for gf in kz.keys():
        for i,zed in enumerate(kz[gf]):
            ii = kcor[gf][:,i]<100
            pylab.plot(kph[gf][ii],kcor[gf][:,i][ii],'-b')
            pylab.plot(kph[gf][ii],kcor[gf][:,i][ii],'ob')
            pylab.text(kph[gf][ii][-1],kcor[gf][:,i][ii][-1],zed)
        pylab.title(sn+'   '+gf)
        raw_input('...')
        pylab.clf()

        for i,ph in enumerate(kph[gf]):
            ii = kcor[gf][i,:]<100
            pylab.plot(kz[gf][ii],kcor[gf][i,:][ii],'-b')
            pylab.plot(kz[gf][ii],kcor[gf][i,:][ii],'ob')
            pylab.text(kz[gf][ii][-1],kcor[gf][i,:][ii][-1],int(ph))
        pylab.title(sn+'   '+gf)
        raw_input('...')
        pylab.clf()

        for f in 'gri':     #!!!!!!!!  test kcor continuity 
            magfit = []
            for z in arange(.1,.8,.01):
                g = assign_filter(f,z,tsn['bands']).f()
                _kcor = kcor_interpolate(tph[g],z,kph[g+f],\
                              kz[g+f],kcor[g+f])
                phobs,magobs = mag_observer_frame(tph[g],\
                                 tabsmag[g],tsn['zed'],z,_kcor,0)
                magfit.append(interp(array([0]),phobs,magobs))
            pylab.plot(arange(.1,.8,.01),magfit,'ro')
            pylab.xlabel('z')
            pylab.title(f)
            raw_input('...')
            pylab.cla()

def plot_template_lc(sn,tph,tabsmag):
    pylab.ion()
    pylab.title(sn)
    off,ymin,ymax = 0.5,-100,0
    phfit = arange(-20,300,1)
    for f in 'UBVR':
        if f in tph.keys():
            pylab.plot(tph[f],tabsmag[f]+off,'o',label=f)
            magfit = interp(phfit,tph[f],tabsmag[f]+off)
            pylab.plot(phfit,magfit,'-')
            ymin = max(ymin,max(tabsmag[f])+off)
            ymax = min(ymax,min(tabsmag[f])+off)
            off += .5
    pylab.ylim(ymin+1,ymax-1)
    pylab.legend(numpoints=1)
    raw_input('.. next ..')

####################################################################
# MAIN ROUTINE
#-------------------------------------------------------------------
def fit_lc(jd,mag,merr,snr,
           tsnlist,             # list of SN templates
           zrange,phrange,avrange,verbose):
    
    ff = open('fit_template.list')          # read template list/info. File with 4 columns: SN | type | vhel | mu
    righe = ff.readlines()
    tsn = {}                    # set empty dictionary
    _alltypes = []              # set empty array
    # The template file is read into a dictionary
    for r in righe:
        if r[0] == '#': continue  # ignores commented lines, jump to the next iteration of for loop
        sn = r.split(None)[0]     # SN template ID: string
        tsn[sn] = {}              # set empty sub-dictionary
        """
        each SN sub-dictionary has this keys
        sn : SNID = {
                 'type' : sn type,
                 'zed'  : redshift,
                 'mu'   : ...
                 }
        """
        tsn[sn]['type'] = r.split()[1]
        tsn[sn]['zed'] = float(r.split()[2])
        if abs(tsn[sn]['zed']) > 1: tsn[sn]['zed'] *= 1/300000.
        tsn[sn]['mu'] = float(r.split()[3])
        _alltypes.append(tsn[sn]['type'])

    alltypes = set(_alltypes)

    if tsnlist[0] == 'all': tsnlist=tsn.keys() 
    elif tsnlist[0] in alltypes: 
        _tsnlist = []           # set empty array
        for t in tsn.keys():    # t goes on SNID values
            if tsn[t]['type'] == tsnlist[0]: _tsnlist.append(t) # extract from template dictionary all SN with specified type
        tsnlist = _tsnlist[:]                                   # overwrite tsnlist (!!! input variable !!!) now it's equal to the list of SNID
     
        
    # read template light curves
    tph,tabsmag = {},{}    # set empty dictionaries
    for sn in tsnlist:
        sndata = alice.leggifile(alice.alice_data+sn)
        tph[sn],tabsmag[sn],tsn[sn]['jdmax'],tsn[sn]['bands'],tsn[sn]['EBV'] =\
                            read_template_lc(sndata,tsn[sn])
        if verbose: 
            print "    SN=",sndata['sn'],tsn[sn]['bands'],\
                " jdmax=",tsn[sn]['jdmax'],\
                "  mu=",sndata['mu'],"  E(B-V)=",'%4.2f' % tsn[sn]['EBV']
            
    if args.debug: #!!!!!!!!
        for sn in tsnlist: plot_template_lc(sn,tph[sn],tabsmag[sn])
        pylab.close()

    tspe = read_template_spectra(tsnlist,tsn)   # read template spectra

    if args.debug: #!!!!!!!!
        for sn in tsnlist: 
            ff = open('_tmp.txt','w')
            for s in tspe[sn]['spec']: ff.write(asa_dir+sn+'/'+s+' \n')
            ff.close()
            iraf.specplot('@_tmp.txt')

    # read or compute k-correction
    kph,kz,kcor = {},{},{}      # set empty dictionaries
    if os.path.exists('kcor.save'):
        ff = open('kcor.save')
        _zra_kk   = pickle.load(ff)
        _zstep_kk = pickle.load(ff)
        if (_zra_kk == zra_kk) and (_zstep_kk == zstep_kk): # zra_kk and zstep_kk defined in the header
            print "Using kcor.save"
            kph  = pickle.load(ff)
            kz  = pickle.load(ff)
            kcor = pickle.load(ff)

    for sn in tsnlist:             
        if sn not in kph.keys(): 
            print 'compute Kcor ',sn
            kph[sn],kz[sn],kcor[sn] = kcor_compute(sn,zra_kk,zstep_kk,\
                  tspe[sn],tsn[sn],verbose) 

    ff = open('kcor.save','w')
    for var in  [zra_kk,zstep_kk,kph,kz,kcor]: pickle.dump(var,ff)

    if args.debug:   #!!!!!!!!  ## plot kcor
        pylab.ion()
        for sn in tsnlist:  
            plot_kkcor(sn,kz[sn],kph[sn],kcor[sn],tsn[sn],tph[sn],tabsmag[sn])

    phmax,zbest,AVbest,chi,np,phfit,magfit = {},{},{},{},{},{},{} # set empty dictionaries

    print " %-5s %-6s %-8s %-15s %-8s %-8s" %\
        ('type','sn','redshift','jdmax','A_V','chi')
    for t in alltypes:
        for sn in tsnlist:                       #  fit template light curves
            if t == tsn[sn]['type']:
                phmax[sn],zbest[sn],AVbest[sn],chi[sn],np[sn] = \
                   best_fit(sn,t,jd,mag,merr,snr,tph[sn],tabsmag[sn],\
                   tsn[sn]['zed'],tsn[sn]['bands'],kph[sn],kz[sn],kcor[sn],\
                   zrange,phrange,avrange)

    # find the lowest chi values from the fits
    snbest = [chi.keys()[argmin(chi.values())]][0]

    dof = 0
    if  zrange[1] >  zrange[0]: dof+1
    if phrange[1] > phrange[0]: dof+1
    if avrange[1] > avrange[0]: dof+1

    # statistical test: other SN type with similar variance .........
    if len(tsnlist)>1:
        prob = {}
        for sn in tsnlist:                      
            prob[sn] = stats.f.pdf(chi[sn]/chi[snbest],\
                                       np[sn]-dof,np[snbest]-dof)
            
        tipi = []
        for sn in tsnlist: tipi.append(tsn[sn]['type'])
        _tipi=set(tipi)
        print 80*'='
        print ">>> also fit:"
        for t in _tipi:
            ss,pp = [],[]
            for sn in tsnlist:                      
                if t == tsn[sn]['type']:
                    ss.append(sn)
                    pp.append(prob[sn])
            imin = argmin(pp)
            if pp[imin] >0.05 and sn != snbest: 
                print '%10s %5s %10s %4.2f'% ('',t,ss[imin],pp[imin])

    #****************
    print '='*80
    print "*** BEST FIT %s %s z=%4.2f [%4.2f,%4.2f]   \n             ph=%4.1f [%4.2f,%4.2f] AV=%3.1f [%3.1f,%3.1f] chi=%5.2f \n"% \
        (snbest,tsn[snbest]['type'],zbest[snbest][0],zbest[snbest][1],\
         zbest[snbest][2],\
         phmax[snbest][0],phmax[snbest][1],phmax[snbest][2],\
       AVbest[snbest][0],AVbest[snbest][1],AVbest[snbest][2],chi[snbest])

    plot_lc(jd,mag,merr,snr,snbest,tph[snbest],tabsmag[snbest],tsn[snbest],\
         kph[snbest],kz[snbest],kcor[snbest],phmax[snbest][0],\
         zbest[snbest][0],AVbest[snbest][0],chi[snbest])

    return

###########################################################################
if __name__ == "__main__":

    if 'm' in args.phrange: args.phrange = '-'+args.phrange[1:]

    zrange  = [float(x) for x in args.zrange.split(',')]
    phrange = [float(x) for x in args.phrange.split(',')]
    avrange = [float(x) for x in args.avrange.split(',')]

    if len(zrange)==1: zrange += [zrange[0]]
    if len(phrange)==1: phrange += [phrange[0]]
    if len(avrange)==1: avrange += [avrange[0]]

    if args.guess:
        zrange .append(.03)
        phrange.append(3.)
        avrange.append(.3)
    else:
        zrange .append(.01)
        phrange.append(1.)
        avrange.append(.1)  

    jd,mag,merr,snr = read_lc_mysql(args.candidate,args.verbose)

    pylab.ion()
    if args.template: 
        tsnlist = args.template.split(',')
        fit_lc(jd,mag,merr,snr,tsnlist,zrange,phrange,avrange,args.verbose)
    else:
        plot_lc(jd,mag,merr,snr,'',{},{},{},{},{},{},0.,'','','')

    print "********** Completed in ",int(time.time()-start_time),"sec"
    _answ = 'n'
    answ = raw_input('>>> save figure ? [n] ')
    if len(answ)>0: _answ = answ  
    if _answ =='y': pylab.savefig('bestfit/'+args.candidate+'.png')
