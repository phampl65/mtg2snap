import hdf5plugin
import numpy as np
from netCDF4 import Dataset, default_fillvals
import glob
import tkinter as tk
import tkinter.font as tkfont
from tkinter import filedialog as fd
import warnings
import os
import sys
import json 
import pandas as pd
import math
import datetime

AUnit = 149597870.700         # astronomical unit in km
jsonfilename = 'mtg2snap.json'

class gvariableclass:
    gui = True
    ncmask = "W*.nc"                 # basic mask to search for NC4 files with MTG chunks
    outputmask = './DT-SLOT-TM-CH-RES.nc'  # output file mask, following codes will be replaced:
                                            # DT - nominal date YYYYMMDD
                                            # TM - nominal starttime as  hhmm in 10 minute step
                                            # SLOT - slot code 001-144 with leading zeros
                                            # CH - channel or group code. If:
                                            #    - every channel in separate file - channel code
                                            #    - every type (vis, nir, ir, vw) in separate file, then typecode
                                            #    - every group (vis+nir, ir+wv) in separate file,
                                            #      then groupcode (visnir, irwv)
                                            # RES - resolution:
                                            #     - LOW = 11136 for vis+nir, 5568 for ir+wv
                                            #     - HI = 22272 for vis+nir, 11136 for ir+vw
                                            #     - reflects original res. of source (FDHSI, HRFI)
    istarttime = 7                    # word position in filename:  scan start  YYYYMMDDhhmmss
    iendtime = 8                      #                             scan end
    islot = 12                        #                             slot number 001-144
    ichunk = 13                       #                             chunk in slot 001-040
    apppath = ''
    jsonnamebase = ''         # base for json configuration filename. Either filled from ARGV or "mtg2snap"
    win = None
    lbox_date = None
    lbox_slot = None
    scrbar_date = None
    scrbar_slot = None
    wordlist = []
    slots_chunks = []
    datum = None
    slot = None
    resolution = None
    irchannels = []
    vischannels = []
    channeldims = {}
    longitudes = {}     # lon, lat and solar zenith angle cosinus arrays dict
    latitudes = {}      # the key for dict is image dimension
    coszenits = {}      #
    times = {}          # pixel time array dict for every dimension
    chandata = {}
    wholelon = None
    wholelat = None
    lbldate = None
    lblslot = None
    fchannels = None
    fselection = None
    filelist = []
    channels_vars = {}
    channels_btns = {}
    filedistribution = None
    lonlat = None 
    costheta = None
    timechan = None
    outputformat =None 
    slotselection = None
    indexmap = None
    folderdates = []
    origdir = ''
    jsondata = {}
    solar_irradiance = {}
    sundist = None
    crops = {}
    cropcode = {}
    selectedcrop = None
    selectedcropname = ''
    croprows = None
    croplines = None
    mindimension = 5568
    pdstartrow = None
    pdendrow = None
    ncattr = {}


#############################################
# #
# #  constant attributes added into every file regardless of data content
# #
staticAtt = {
    'summary': 'Flexible Combined Imager (FCI) Level 1c Rectified Radiance dataset',
    'format': 'NC4'

}

#############################################
# #
# #   global attrs transferred into every file
# #
transferredAtt = [
    'title','keywords', 'institution','location_indicator','data_designator',
    'platform','data_source','processing_level','type','subtype',
    'processor_version','algorithm_version','format_version','source','facility_or_tool',
    'references','creator_type','creator_institution','creator_name','creator_email','creator_url',
    'standard_name_vocabulary','project','cdm_datatype','processing_mode','disposition_mode',
    'environment','repeat_cycle_in_day','date_time_position','time_position','group_tag',
    'product_id','instrument_configuration_id','instrument_configuration_id_version',
    'mtg_name'
]

###############################################
# #
# #   transferred channel attributes
channelAtt = [
    'central_wavelength_specified','spectral_width_specified',
    'central_wavelength_actual', 'spectral_width_actual',
    'ssd'
]

###############################################
# #
# #  transferred channel attributes from the "Measured" group
# #
channelMeasuredAtt = [
    'radiance_unit_conversion_coefficient',
    'radiance_to_bt_conversion_coefficient_a',
    'radiance_to_bt_conversion_coefficient_b',
    'radiance_to_bt_conversion_coefficient_wavenumber',
    'radiance_to_bt_conversion_constant_c1',
    'radiance_to_bt_conversion_constant_c2',
    'channel_effective_solar_irradiance'
]
    
    
gvars = gvariableclass()

##################################################
# #
# #    this is something copied from MTG reader
# #
# #
#def _ensure_dataarray(arr):
#    if not isinstance(arr, xr.DataArray):
#        attrnames = arr.ncattrs()
#        attrs = { i:arr.getncattr(i) for i in attrnames }
#        pp = xr.DataArray(da.from_array(arr), dims=arr.dimensions, attrs=attrs, name=arr.name)
#    return pp

#####################################################
# #
# #  store different attributes from source file to global variable
# #
def storeStaticAttrs():
    global gvars, staticAtt
    for att in staticAtt:
        if att not in gvars.ncattr:
            gvars.ncattr[att]= staticAtt[att]
########################################################
# #
# #  saves all stored attributes to the Dataset as global attributes
# #
def storeGlobalAtt(ds):
    global gvars, transferredAtt
    for att in transferredAtt:
        if att not in gvars.ncattr:
            v = ds.getncattr(att)
            gvars.ncattr[att]=v
        
#######################################################
# #
# #    stores saved channel attributes to specified channel as its locals
# #
def storeChannelAtt(ds,ch):
    global gvars, channelAtt, channelMeasuredAttr
    dt = ds.groups['data']
    vs = dt.groups[ch]
    for vname in channelAtt:
        fullname = ch+'_'+vname
        if fullname not in gvars.ncattr:
            v=vs.variables[vname]
            vl = str(v[0])
            lst = v.ncattrs()
            if 'units' in lst:
                vl = vl + ' ' + v.units
            gvars.ncattr[fullname]=vl
    ms = vs.groups['measured']
    for vname in channelMeasuredAtt:
        fullname = ch+'_'+vname
        if fullname not in gvars.ncattr:
            v=ms.variables[vname]
            vl = str(v[0])
            lst = v.ncattrs()
            if 'units' in lst:
                vl = vl + ' ' + v.units
            gvars.ncattr[fullname]=vl
    
#############################################
# #
# #   save all stored attrs into destination file
# #
def saveAttrs(dsout):
    global gvars
    for att in gvars.ncattr:
        dsout.setncattr(att,gvars.ncattr[att])


########################################
#  procedure used for sorting
def sortCycles(x):
    return '_'.join(x)

##########################################
# # procedure returns text form of nominal start time for the given slot
# # ex.:  0000, 0010, 0020, ... , 2340, 2350
# # 
def getStartTimeForSlot(sl):
    isl= int(sl)
    poc_h= (isl-1) // 6
    poc_m= ((isl-1) % 6) * 10
    res = "{:02d}{:02d}".format(poc_h,poc_m)
    return res

###########################################
# #   procedure returns text form of nominal time interval for the given slot
# #   ex.:  0000 -> 0010
# #
def getTimeIntervalForSlot(sl):
    isl= int(sl)
    poc_h= (isl-1) // 6
    poc_m= ((isl-1) % 6) * 10
    kon_h = isl // 6
    kon_m = (isl % 6) *10
    res = "{:02d}:{:02d} -> {:02d}:{:02d}".format(poc_h,poc_m,kon_h,kon_m)
    return res

##################################################
# #   procedure returns number of chunks available for the given slot 
# #  in the gvars.wordlist
# #
def getNumberOfChunks(dt,slot):
    global gvars
    chunks = [x[13] for x in gvars.wordlist
              if (x[7][:8] == dt) 
              and (x[12] == slot[0]) 
              and (x[1][5]==slot[2])
              and (x[1][9] == 'BODY')
              ]
    return len(chunks)
    
########################################################
# #  the procedure displays list of available slots and number of their chunks
# #  for selected date / GUI
# #
def getCyclesForDate(dt):   # dt is datestring in the format YYYYMMDD
    global gvars
    gvars.lbox_slot.delete(0,tk.END)
    slots = list(set([(x[12], 
                       getTimeIntervalForSlot(x[12]) ,
                       x[1][5]) 
                      for x in gvars.wordlist 
                      if (x[7][:8] == dt)]))
    slots.sort(key = sortCycles)
    gvars.slots_chunks = [(x[0],x[1],x[2],getNumberOfChunks(dt,x)) for x in slots]
    for x in gvars.slots_chunks:
        #print(x)
        (a,b,c,d) = x
        gvars.lbox_slot.insert(tk.END,'{:4s} {:s} {: ^5} {:3d}'.format(a,b,c,d))
        
        
##############################################################
# #  procedure returns list of slots available for given date and resolution
# #
def getCyclesForDateAndRes(dt,res):
    global gvars
    slots = list(set([x[12] for x in gvars.wordlist
             if (x[7][:8] == dt) and (x[1][5] == res)
             ]))
    slots.sort()
    return slots
        
##############################################################
# # 
# #   when in GUI is selected different date
# #  this proc finds all slots and chunks available
# #  and displays them on GUI
# 
def datechange(event):
    global gvars
    t = gvars.lbox_date.curselection()
    if len(t)>0:
        i=t[0] 
        txt = gvars.lbox_date.get(i)
        gvars.datum = txt
        gvars.lbldate.config(text = txt)
        gvars.lblslot.config(text = '')
        print(i, txt)
        getCyclesForDate(txt)
        
#########################################################
# #
# #  read channel dimensions either from TRAIL or from BODY files
# #
def getAllDims(dt,slot, resolution):
    global gvars
    traillist = [x for x in gvars.wordlist
               if (x[7][:8] == dt)
               and (x[12]==slot)
               and (x[1][5]==resolution)
               and (x[1][9]=='TRAIL')
               ]
    if len(traillist)>0:
        print('number of files =', len(traillist))
        line = traillist[0].copy()
        line[1]='-'.join(line[1])
        name = '_'.join(line)
        print('getting info from: ',name)
        name = name + '.nc'
        getDimsFromTrail(name)
    else:
        print('TRAIL not found')
        bodylist = [x for x in gvars.wordlist
                   if (x[7][:8] == dt)
                   and (x[12]==slot)
                   and (x[1][5]==resolution)
                   and (x[1][9]=='BODY')
                   ]
        if len(bodylist)>0:
            print('number of files =', len(bodylist))
            line = bodylist[0].copy()
            line[1]='-'.join(line[1])
            name = '_'.join(line)
            print('getting info from: ', name)
            name = name + '.nc'
            getDimsFromBody(name)
        else:
            print('even BODY not found, some strange error!')
        
################################################################
# #
# #   gets list of files for given date, slot and resolution
# #  and stores it in gvars.filelist
# #
def getFileList(dt,slot, resolution):
    global gvars
    bodylist = [x for x in gvars.wordlist
               if (x[7][:8] == dt)
               and (x[12]==slot)
               and (x[1][5]==resolution)
               and (x[1][9]=='BODY')
               ]
    gvars.filelist = []
    if len(bodylist)>0:
        print('number of files =', len(bodylist))
        for r in bodylist:
            line = r.copy()
            line[1]='-'.join(r[1])
            name = '_'.join(line)
            name = name + '.nc'
            gvars.filelist.append(name)

##########################################################
# #
# #   reads dimensions of channels from the TRAIL nc4 file
# # 
def getDimsFromTrail(fname):
    global gvars
    ds = Dataset(fname,'r')
    channels = ds.getncattr('subsettable_groups_present')
    channels = channels.split()
    channels = [x for x in channels if x[:3] in ['vis','ir_','wv_','nir']]
    print(channels)
    gvars.irchannels = [i for i in channels if (i[:3]=='ir_' or i[:3]=='wv_')]
    print(gvars.irchannels)
    gvars.vischannels = [i for i in channels if (i[:4]=='nir_' or i[:4]=='vis_')]
    print(gvars.vischannels)
    dt = ds.groups['data']
    gvars.channeldims = {}
    for ch in channels:
        cg = dt.groups[ch]
        ch_xdim = cg.variables['number_of_columns'][:].data.item()
        ch_ydim = cg.variables['number_of_rows'][:].data.item()
        ch_y0 = cg.variables['start_row_number'][:].data.item()
        ch_y1 = cg.variables['end_row_number'][:].data.item()
        tpl = (ch_xdim, ch_ydim, ch_y0, ch_y1)
        #print("proc getdimsfromtrail: ",ch,' - ', tpl)
        gvars.channeldims[ch]= tpl
    
    ds.close()

##########################################################
# # 
# #   reads dimensions of channels from the BODY nc4 files
# #
def getDimsFromBody(fname):
    global gvars
    ds = Dataset(fname,'r')
    channels = ds.getncattr('subsettable_groups_present')
    channels = channels.split()
    channels = [x for x in channels if x[:3] in ['vis','ir_','wv_','nir']]
    print(channels)
    gvars.irchannels = [i for i in channels if (i[:3]=='ir_' or i[:3]=='wv_')]
    print(gvars.irchannels)
    gvars.vischannels = [i for i in channels if (i[:4]=='nir_' or i[:4]=='vis_')]
    print(gvars.vischannels)
    dt = ds.groups['data']
    gvars.channeldims = {}
    for ch in channels:
        cg = dt.groups[ch]
        ch_xdim = cg.dimensions['x'].size
        ch_ydim = ch_xdim
        tpl = (ch_xdim, ch_ydim)
        #print('proc getdimsfrombody: ', ch,' - ', tpl)
        gvars.channeldims[ch]= tpl
    
    ds.close()

####################################################################
# # 
# #  reads start/endrow for every chunk and every channel
# #
# #  this is to decide which chunks to use when making crops
# #
def getChunkRows():
    global gvars
    chlist = [str(i+1) for i in range(40)]
    gvars.pdstartrow=pd.DataFrame(index = chlist)
    gvars.pdendrow=pd.DataFrame(index = chlist)
    for fl in gvars.filelist:
        chunk = str(int(fl[-5:-3]))
        print(chunk,end=' ', flush=True)
        ds = Dataset(fl, 'r')
        channels = ds.getncattr('subsettable_groups_present')
        channels = channels.split()
        dt = ds.groups['data']
        for ch in channels:
            chgrp = dt.groups[ch]
            meas = chgrp.groups['measured']
            startrow = meas.variables['start_position_row'][:]
            endrow = meas.variables['end_position_row'][:]
            gvars.pdstartrow.at[chunk,ch] = startrow
            gvars.pdendrow.at[chunk,ch] = endrow
        ds.close()
    print()
    #print('startrow')
    #print(gvars.pdstartrow.iloc[:,:13])
    #print('endrow')
    #print(gvars.pdendrow.iloc[:,:13])
            

######################################################################
# #
# #  reads X/Y variables from NC4 chunk and computes lon/lat arrays from them
# # 
def getGeoChunk(DS,chname):
    global  gvars, curlon, curlat
    
    dt = DS.groups['data']

    vs = dt.groups[chname]
    meas = vs.groups['measured']
    startrow = meas.variables['start_position_row'][:]
    endrow = meas.variables['end_position_row'][:]
    startcol = meas.variables['start_position_column'][:]
    endcol = meas.variables['end_position_column'][:]

    # decide, if the chunk and crop collide
    yn,ys = gvars.croprows
    if endrow<ys or startrow>yn:
        return 

    # projection parameters
    # 
    projvar = dt.variables['mtg_geos_projection']
    axis_a = float(projvar.getncattr('semi_major_axis'))
    axis_b = float(projvar.getncattr('semi_minor_axis'))
    sat_height = float(projvar.getncattr('perspective_point_height'))
    axis_sat = sat_height + axis_a
    sat_lon = float(projvar.getncattr('longitude_of_projection_origin'))

    pole_lams = meas.variables['x']
    pole_fis = meas.variables['y']
    sinfis = np.sin(pole_fis)
    sinfis = np.transpose(np.atleast_2d(sinfis))
    sin2fis = np.square(sinfis)
    cosfis = np.cos(pole_fis)
    cosfis = np.transpose(np.atleast_2d(cosfis))
    cos2fis = np.square(cosfis)
    sinlams = np.atleast_2d(np.sin(pole_lams))
    coslams = np.atleast_2d(np.cos(pole_lams))
    coslamscosfis = coslams * cosfis
    s4 = (axis_a*axis_a)/(axis_b*axis_b)
    s5 = (axis_sat*axis_sat)-(axis_a*axis_a)
    sd = np.sqrt(np.square(axis_sat*coslamscosfis) -((cos2fis+(s4*sin2fis))*s5))
    sn = ((axis_sat*coslamscosfis)-sd)/(cos2fis+(s4*sin2fis))
    s3 = sn*sinfis
    s1 = axis_sat-(sn*coslamscosfis)
    s2 = 0-(sn*sinlams*cosfis)
    s3 =sn*sinfis
    sxy = np.hypot(s1,s2)
    lon = np.rad2deg(np.arctan(s2/s1))+ sat_lon
    lat = np.rad2deg(np.arctan(s4*s3/sxy))
    dim_x,dim_y = gvars.channeldims[chname][:2]
    gvars.longitudes[dim_x][dim_y-endrow:dim_y-startrow+1,startcol-1:endcol] = np.flip(lon,axis=0)
    gvars.latitudes[dim_x][dim_y-endrow:dim_y-startrow+1,startcol-1:endcol] = np.flip(lat,axis=0)


#################################################
# #
# #   reads the nc4 variable index_map, subtracts the offset to make it start from zero,
# #   and stores it in gvars.indexmap
# #
def getIndexMapForChunk(DS,chname):
    global gvars
    indexsize = DS.dimensions['index'].size
    indexoffset = DS.variables['index_offset'][:].data.item()
    dt = DS.groups['data']
    vs = dt.groups[chname]
    meas = vs.groups['measured']
    gvars.indexmap = meas.variables['index_map'][:,:]
    gvars.indexmap = gvars.indexmap - indexoffset
    #  WORKAROUND BECAUSE OF DATA ERROR / there exist indexvalues that are larger than the array dim
    gvars.indexmap = np.where(gvars.indexmap<indexsize, gvars.indexmap, 0)
    gvars.indexmap = np.where(gvars.indexmap<indexsize, gvars.indexmap, 0)
    #  END OF WORKAROUND
    #gvars.indexmap = np.flip(gvars.indexmap, axis=0)

##########################################################
# #
# #    transforms timedif from 1.1.2000 0:0:0
# #    into integer in numeric format DHHmmss
# #    The day is there so that even low datetimes have all places
# #
def visualiseTime(t):
    base = datetime.datetime(2000,1,1,0,0,0)
    y=t.shape
    t2=np.zeros(y,dtype='u4')
    for i in range(y[0]):
        delta = datetime.timedelta(seconds=t[i]-86400)  # ten posun je o jeden den vetsi nez by mel byt
        mt=base+delta
        txt=mt.strftime('%d%H%M%S')
        t2[i]=int(txt)
    return t2



###############################################################
# # 
# #  reads the pixel timestamp and transforms it into a readable numeric format
# #  it is stored in gvars.times for each dimension available
# #
def getTimeChunk(DS,chname):
    global  gvars, curlon, curlat
    
    dt = DS.groups['data']

    vs = dt.groups[chname]
    meas = vs.groups['measured']
    startrow = meas.variables['start_position_row'][:]
    endrow = meas.variables['end_position_row'][:]
    startcol = meas.variables['start_position_column'][:]
    endcol = meas.variables['end_position_column'][:]

    # decide, if the chunk and the crop collide
    yn,ys = gvars.croprows
    if endrow<ys or startrow>yn:
        return 

    getIndexMapForChunk(DS,chname)
    timearray = DS.variables['time'][:]
    timeintarray = visualiseTime(timearray.data)
    pixtimes = timeintarray[gvars.indexmap]



    dimension_x,dimension_y = gvars.channeldims[chname][:2]
    gvars.times[dimension_x][dimension_y-endrow:dimension_y-startrow+1,startcol-1:endcol] = np.flip(pixtimes,axis=0)

###########################################################
# #
# #  creates timestamp array for each pixel
# #
def getTimeForChannel(chname):
    global gvars
    doit = False
    wasmade = False
    #print('proc getTimeForChannel')
    dim_x, dim_y = gvars.channeldims[chname][:2]
    if (dim_x in gvars.times):
        wasmade = True
    else:
        celektim = np.zeros((dim_y,dim_x),np.float32)
        gvars.times[dim_x]=celektim
    if not wasmade:
        print('    ----> TIME')
        for fl in gvars.filelist:
            chunk = str(int(fl[-5:-3]))
            yn,ys = gvars.croprows
            if gvars.pdstartrow.at[chunk,chname]<=yn and gvars.pdendrow.at[chunk,chname]>=ys:
                print(int(fl[-5:-3]), end=' ', flush=True)
                ds = Dataset(fl, 'r')
                getTimeChunk(ds,chname)
                wasmade = True
                ds.close()
            
        print()
    #print('wasmade =',wasmade)
    return wasmade
   
#####################################################
# #
# #   creates LON and LAT array for each pixel
# #  
def getGeoForChannel(chname):
    global gvars
    doit = False
    wasmade = False
    #print('proc getGeoForChannel')
    dim_x, dim_y = gvars.channeldims[chname][:2]
    if (dim_x in gvars.longitudes) and (dim_x not in gvars.latitudes):
        wasmade = True
    if dim_x not in gvars.longitudes:
        wholelon = np.zeros((dim_y,dim_x),np.float32)
        gvars.longitudes[dim_x]=wholelon
    if dim_x not in gvars.latitudes:
        wholelat = np.zeros((dim_y,dim_x),np.float32)
        gvars.latitudes[dim_x]=wholelat
    if not wasmade:
        wasmade = wasmade or getTimeForChannel(chname)
        if wasmade:
            print('    ----> LONLAT')
            for fl in gvars.filelist:
                chunk = str(int(fl[-5:-3]))
                yn,ys = gvars.croprows
                if gvars.pdstartrow.at[chunk,chname]<=yn and gvars.pdendrow.at[chunk,chname]>=ys:
                    print(int(fl[-5:-3]), end=' ', flush=True)
                    ds = Dataset(fl, 'r')
                    getGeoChunk(ds,chname)
                    wasmade = True
                    ds.close()
            
        print()
    #print('wasmade =',wasmade)
    return wasmade
            
######################################################
# #  
# # compute solar zenith angle cos for pixels in chunk
# # 
def getZenitChunk(DS,chname):

    global gvars, celekcosr, curcosr

    #indexsize = DS.dimensions['index'].size

    dim_x,dim_y = gvars.channeldims[chname][:2]
    dt = DS.groups['data']

    vs = dt.groups[chname]
    meas = vs.groups['measured']
    startrow = meas.variables['start_position_row'][:]
    endrow = meas.variables['end_position_row'][:]
    startcol = meas.variables['start_position_column'][:]
    endcol = meas.variables['end_position_column'][:]

    # decide, if crop and chunk collide
    yn,ys = gvars.croprows
    if endrow<ys or startrow>yn:
        return 

    
    
    state = DS.groups['state']
    celes = state.groups['celestial']
    sun_long = celes.variables['subsolar_longitude'][:]
    sun_lat  = celes.variables['subsolar_latitude'][:]

    # read the coordinates for chunk
    curlat = gvars.latitudes[dim_x][dim_y-endrow:dim_y-startrow+1,startcol-1:endcol]
    curlon = gvars.longitudes[dim_x][dim_y-endrow:dim_y-startrow+1,startcol-1:endcol]

    # COSINUS of SOLAR ZENITH ANGLE
    # 
    getIndexMapForChunk(DS, chname)
    indexmap = np.flip(gvars.indexmap, axis=0)
    sunlong = sun_long[indexmap]
    sunlat = sun_lat[indexmap]
    #print('sunlong:',sunlong[0,2000],'     sunlat:',sunlat[0,2000])
    sunlat = np.deg2rad(sunlat)
    lat = np.deg2rad(curlat)
    diflon = np.deg2rad(curlon-sunlong)
    cosr = (np.sin(lat)*np.sin(sunlat)) + (np.cos(lat)*np.cos(sunlat)*np.cos(diflon))
    cosr = np.where(cosr>0, cosr,0)
    gvars.coszenits[dim_x][dim_y-endrow:dim_y-startrow+1,startcol-1:endcol] = cosr 
    # ^^^ there is no need for flip here, source array is correctly oriented already
    curcosr = cosr


######################################################
# #
# #   create an array of solar zenith angle cos
# #
def getZenithForChannel(chname):
    global gvars
    #print('proc getZenithForChannel')
    wasmade = False
    dim_x, dim_y = gvars.channeldims[chname][:2]
    #print('    dim X,Y:', dim_x, dim_y)
    #print('    COStheta hotove pro:', gvars.coszenits.keys())
    if dim_x not in gvars.coszenits:
        celekzen = np.zeros((dim_y,dim_x),np.float32)
        gvars.coszenits[dim_x]=celekzen
    else:
        wasmade=True
    if not wasmade:
        wasmade = wasmade or getGeoForChannel(chname)
        if wasmade:
            print('    ----> COS SOLAR ZENIT')
            #print(gvars.filelist)
            #print('............................')
            for fl in gvars.filelist:
                chunk = str(int(fl[-5:-3]))
                yn,ys = gvars.croprows
                if gvars.pdstartrow.at[chunk,chname]<=yn and gvars.pdendrow.at[chunk,chname]>=ys:
                    print(int(fl[-5:-3]), end=' ', flush=True)
                    ds = Dataset(fl, 'r')
                    getZenitChunk(ds,chname)
                    wasmade = True
                    ds.close()
            print()
    #print('wasmade =',wasmade)
    return wasmade


########################################################
# #
# #  get IR/WV channel data for one chunk
# #
def getThermalChunk(DS,chname):
    global gvars, celek, curdata
    dt = DS.groups['data']
    vs = dt.groups[chname]
    meas = vs.groups['measured']
    startrow = meas.variables['start_position_row'][:]
    endrow = meas.variables['end_position_row'][:]
    #print('rows: ',startrow, endrow)
    startcol = meas.variables['start_position_column'][:]
    endcol = meas.variables['end_position_column'][:]
    #print('cols: ',startcol, endcol)

    # there is the decision here, if the crop contains some of those lines 
    yn,ys = gvars.croprows
    if endrow<ys or startrow>yn:
        return 

    dim_x, dim_y = gvars.channeldims[chname][:2]

    state = DS.groups['state']
    celes = state.groups['celestial']
    sun_dist = celes.variables['earth_sun_distance'][0]
    sun_dist = sun_dist / AUnit
    gvars.sundist = sun_dist
    
    v = meas.variables['effective_radiance']
    koef_vc = meas.variables['radiance_to_bt_conversion_coefficient_wavenumber'][:]
    koef_a = meas.variables['radiance_to_bt_conversion_coefficient_a'][:]
    koef_b = meas.variables['radiance_to_bt_conversion_coefficient_b'][:]
    koef_c1 = meas.variables['radiance_to_bt_conversion_constant_c1'][:]
    koef_c2 = meas.variables['radiance_to_bt_conversion_constant_c2'][:]
    koef_c2vc = koef_c2 * koef_vc
    koef_c1vc3 = koef_c1 * (koef_vc**3)
    koef_ba = koef_b / koef_a
    attrnames = v.ncattrs()
    attrs = { i:v.getncattr(i) for i in attrnames }
    pole = v[:].data
    fv = default_fillvals.get(pole.dtype.str[1:], np.nan)
    fv=0
    vr = attrs.get('valid_range',[-np.inf, np.inf])
    pole = np.where(pole>= vr[0], pole, fv)
    pole = np.where(pole<= vr[1], pole, fv)
    indata = np.asarray(pole.data)
    msk = np.where(indata==fv,0,1)
    #indata = np.where(msk, indata*konv*1e6,0)  # prevod jednotek, i na metry z mikronu
    #indata = np.where(msk,(koef_k2/(koef_a*np.log(1+(koef_k1/indata))))-koef_ba,0)
    indata = np.where(msk,(koef_c2vc/(koef_a*np.log(1+(koef_c1vc3/indata))))-koef_ba,0)
    #print()
    gvars.chandata[chname][dim_y-endrow:dim_y-startrow+1,startcol-1:endcol] = np.flip(indata,axis=0)
    curdata = indata

#################################################
# #
# #    get IR/WV channel data for scene
# #
def getThermalForChannel(chname):
    global gvars
    doit = True
    wasmade = False
    attrsWereStored = False
    print('    ---->  ', chname)
    dim_x, dim_y = gvars.channeldims[chname][:2]
    print('    dim X,Y:', dim_x, dim_y)
    #print('    COStheta made pro:', gvars.coszenits.keys())
    celek = np.zeros((dim_y,dim_x),np.float32)
    gvars.chandata[chname]=celek
    #doit = doit or getGeoForChannel(chname)
    #print('doit1 =',doit)
    doit = doit and getZenithForChannel(chname)
    #print('doit2 =',doit)
    if doit:
        #print('............................')
        for fl in gvars.filelist:
            chunk = str(int(fl[-5:-3]))
            yn,ys = gvars.croprows
            if gvars.pdstartrow.at[chunk,chname]<=yn and gvars.pdendrow.at[chunk,chname]>=ys:
                print(int(fl[-5:-3]), end=' ', flush=True)
                ds = Dataset(fl, 'r')
                if not attrsWereStored:
                    storeGlobalAtt(ds)
                    storeChannelAtt(ds,chname)
                    attrsWereStored = True
                getThermalChunk(ds,chname)
                wasmade = True
                ds.close()
        print()
    # if no chunk collided with crop, the array is removed
    if not wasmade:
        gvars.chandata.pop(chname)
    return wasmade
            

##########################################################
# #
# #     get VIS/NIR channel data for chunk
# #
def getVisChunk(DS,chname):
    global gvars,celek, celekcorr

    dim_x,dim_y = gvars.channeldims[chname][:2]

    indexsize = DS.dimensions['index'].size

    dt = DS.groups['data']
    vs = dt.groups[chname]
    meas = vs.groups['measured']

    startrow = meas.variables['start_position_row'][:]
    endrow = meas.variables['end_position_row'][:]
    #print('rows: ',startrow, endrow)
    startcol = meas.variables['start_position_column'][:]
    endcol = meas.variables['end_position_column'][:]
    #print('cols: ',startcol, endcol)

    # decide if chunk and crop collide
    yn,ys = gvars.croprows
    if endrow<ys or startrow>yn:
        return 


    state = DS.groups['state']
    celes = state.groups['celestial']
    sun_dist = celes.variables['earth_sun_distance'][0]
    sun_dist = sun_dist / AUnit
    gvars.sundist = sun_dist
    
    koef_ilam = meas.variables['channel_effective_solar_irradiance'][:]
    gvars.solar_irradiance[chname]=koef_ilam
    koef_multi = np.pi * gvars.sundist * gvars.sundist / koef_ilam


    v = meas.variables['effective_radiance']
    attrnames = v.ncattrs()
    attrs = { i:v.getncattr(i) for i in attrnames }
    #pole = _ensure_dataarray(v)
    pole = v[:].data
    fv = default_fillvals.get(pole.dtype.str[1:], np.nan)
    #print('fillval = ',fv)
    fv=0
    vr = attrs.get('valid_range',[-np.inf, np.inf])
    #print('vr = ',vr)
    pole = np.where(pole>= vr[0], pole, fv)
    pole = np.where(pole<= vr[1], pole, fv)
    indata = np.asarray(pole.data)
    msk = np.where(indata==fv,0,1)
    indatacorr = np.where(msk,(koef_multi*indata),0)
    gvars.chandata[chname][dim_y-endrow:dim_y-startrow+1,startcol-1:endcol] = np.flip(indatacorr,axis=0)
    

def getVisForChannel(chname):
    global gvars
    doit = True
    wasmade = False
    attrsWereStored=False
    dim_x, dim_y = gvars.channeldims[chname][:2]
    #print('    dim X,Y:', dim_x, dim_y)
    #print('    COStheta hotove pro:', gvars.coszenits.keys())
    celek = np.zeros((dim_y,dim_x),np.float32)
    gvars.chandata[chname]=celek
    #doit = doit or getGeoForChannel(chname)
    #print('doit1 =',doit)
    doit = doit and getZenithForChannel(chname)
    #print('doit2 =',doit)
    print('    ---->  ', chname)
    if doit:
        #print('............................')
        for fl in gvars.filelist:
            chunk = str(int(fl[-5:-3]))
            yn,ys = gvars.croprows
            if gvars.pdstartrow.at[chunk,chname]<=yn and gvars.pdendrow.at[chunk,chname]>=ys:
                print(int(fl[-5:-3]),end=' ',flush=True)
                ds = Dataset(fl, 'r')
                if not attrsWereStored:
                    storeGlobalAtt(ds)
                    storeChannelAtt(ds,chname)
                    attrsWereStored = True
                getVisChunk(ds,chname)
                wasmade = True
                ds.close()
        print()
    # if there is no chunk that fits, get rid of the array.
    if not wasmade:
        gvars.chandata.pop(chname)
    return wasmade    
            





#  GUI - selects all IR and WV channels
def setIRall():
    global gvars
    for ch in gvars.channels_vars:
        if ch[:3] in ['ir_', 'wv_']:
            #print(ch)
            gvars.channels_vars[ch].set(1) 

# GUI - DEselects all IR and WV channels
def setIRnone():
    global gvars
    for ch in gvars.channels_vars:
        if ch[:3] in ['ir_', 'wv_']:
            #print(ch)
            gvars.channels_vars[ch].set(0)
            
# GUI - selects all VIS and NIR channels
def setVISall():
    global gvars
    for ch in gvars.channels_vars:
        if ch[:4] in ['vis_', 'nir_']:
            #print(ch)
            gvars.channels_vars[ch].set(1)

# GUI - DEselects all VIS and NIR channels
def setVISnone():
    global gvars
    for ch in gvars.channels_vars:
        if ch[:4] in ['vis_', 'nir_']:
            #print(ch)
            gvars.channels_vars[ch].set(0)
            
             
# GUI - displays available channels regarding to the resolution
def listAvailChannels():
    global gvars
    if gvars.fchannels is not None:
        gvars.fchannels.destroy()
    btnfont = tkfont.Font(family = 'Lucida Sans Typewriter', size = 10)
    gvars.fchannels = tk.Frame(gvars.fselection, pady =10)
    gvars.fchannels.columnconfigure(2,minsize = 30)
    gvars.channels_btns.clear()
    gvars.channels_vars.clear()
    btnIRall = tk.Button(gvars.fchannels, text='All', command = setIRall)
    btnIRnone = tk.Button(gvars.fchannels, text='None', command = setIRnone)
    btnVISall = tk.Button(gvars.fchannels, text='All', command = setVISall)
    btnVISnone = tk.Button(gvars.fchannels, text='None', command = setVISnone)
    btnIRall.grid(in_ = gvars.fchannels, column = 0, row = 0, sticky = tk.E)
    btnIRnone.grid(in_ = gvars.fchannels, column = 1, row = 0, sticky = tk.W)
    btnVISall.grid(in_ = gvars.fchannels, column = 3, row = 0, sticky = tk.E)
    btnVISnone.grid(in_ = gvars.fchannels, column = 4, row = 0, sticky = tk.W)
    row = 1
    for ch in gvars.irchannels:
        #print(ch)
        ctlvar = tk.IntVar()
        gvars.channels_vars[ch]=ctlvar
        chk = tk.Checkbutton(gvars.fchannels,variable = gvars.channels_vars[ch], 
                             onvalue = 1, offvalue = 0, text = ch, font=btnfont)
        chk.grid(in_ = gvars.fchannels, column = 0, row = row, sticky = tk.W)
        gvars.channels_btns[ch]=chk
        res = str(gvars.channeldims[ch][0])
        lbl = tk.Label(text = res)
        lbl.grid(in_ = gvars.fchannels, column=1, row=row)
        row += 1
    row = 1
    for ch in gvars.vischannels:
        #print(ch)
        ctlvar = tk.IntVar()
        gvars.channels_vars[ch]=ctlvar
        chk = tk.Checkbutton(gvars.fchannels,variable = gvars.channels_vars[ch], 
                             onvalue = 1, offvalue = 0, text = ch, font=btnfont)
        chk.grid(in_ = gvars.fchannels, column = 3, row = row, sticky = tk.W)
        gvars.channels_btns[ch]=chk
        res = str(gvars.channeldims[ch][0])
        lbl = tk.Label(text = res)
        lbl.grid(in_ = gvars.fchannels, column=4, row=row)
        row += 1
    gvars.fchannels.pack(side = tk.TOP, fill=tk.X)
    
##########################################################
# #
# #  GUI / when user selects some slot in the list, it prints its properties and
# #        displays available channels to select
# #
def changeTheSlot(event):
    global gvars
    t=gvars.lbox_slot.curselection()
    if len(t)>0:
        i=t[0]
        tup = gvars.slots_chunks[i]
        gvars.slot, casint, gvars.resolution, chunks = tup
        gvars.lblslot.config(text = 'slot '+gvars.slot+': '+casint+'  '
                             +gvars.resolution+', chunks: '+str(chunks))
        getAllDims(gvars.datum, gvars.slot, gvars.resolution)
        print('================\nIR channels')
        for i in gvars.irchannels:
            print(i,' ',gvars.channeldims[i][:2])
        print('================\nVIS channels')
        for i in gvars.vischannels:
            print(i,' ',gvars.channeldims[i][:2])
        
        listAvailChannels()

######################################################
# #
# #    proc saves channel into NC4 file
def saveChannel(chname,dsout,tplcrop):
    global gvars
    #print('gvars.chandata =',list(gvars.chandata.keys()))
    if chname in gvars.chandata:
        pole = gvars.chandata[chname]
        #typ = pole.dtype
        typ = 'u2'
        #print('saving ',chname, ' into file ',dsout.filepath())
        vout = dsout.createVariable(chname,typ,dimensions=('y','x'), compression = 'zlib', fill_value=0)
        if chname.startswith('ir_') or chname.startswith('wv_'):
            vout.setncattr('scale_factor', 0.01)  # 1/100 of degree
            vout.setncattr('add_offset',0.0)
            vout[:,:] = pole[tplcrop[2]:tplcrop[3],tplcrop[0]:tplcrop[1]]
        else:
            # VIS a NIR se ulozi nejdriv do vlastniho u2 pole a to se zkopiruje do promenne.
            # teprve potom se nastavi scale a offset
            datamax = np.max(gvars.chandata[chname])
            nonzerotpl = np.nonzero(gvars.chandata[chname])
            if len(nonzerotpl)>0:
                datamin =np.min(gvars.chandata[chname][nonzerotpl])
            else:
                datamin=0
            #print('datamin = ', datamin)
            #print('datamax = ', datamax)
            scale = (datamax-datamin)/8190
            offs = datamin - scale
            res = (pole - offs)/scale
            res = np.where(pole>0, res,0)
            intpole = res.astype(np.uint16)
            vout.set_auto_scale(False)
            vout[:,:] = intpole[tplcrop[2]:tplcrop[3],tplcrop[0]:tplcrop[1]]
            vout.setncattr('scale_factor', scale)   # reflectance precision
            vout.setncattr('add_offset',offs-scale)
            vout.set_auto_scale(True)

        if chname.startswith('ir_') or chname.startswith('wv_'):
            vout.long_name = chname + ' radiation temperature'
        if chname.startswith('nir_') or chname.startswith('vis_'):
            vout.long_name = chname + ' radiance'
        dsout.sync() 
        if chname.startswith('nir_') or chname.startswith('vis_'):
            dsout.setncattr(chname+'_solar_irradiance', gvars.solar_irradiance[chname])

###################################################
# #
# #   proc saves LON and LAT into NC4 file
# #
def saveGeo(res,dsout,tplcrop):
    global gvars
    pole = gvars.longitudes[res]
    typ = pole.dtype
    xw,xe,yn,ys=tplcrop
    vout = dsout.createVariable('longitude',typ,dimensions=('y','x'), compression = 'zlib', least_significant_digit = 4, fill_value = -2147483648)
    vout.setncattr('valid_max',180.0)
    vout.setncattr('valid_min',-180.0)
    vout[:,:] = pole[yn:ys,xw:xe]
    vout.setncattr('units','degrees_east')
    vout.set_auto_mask(True)
    pole = gvars.latitudes[res]
    vout = dsout.createVariable('latitude',typ,dimensions=('y','x'), compression = 'zlib', least_significant_digit = 4,  fill_value = -2147483648)
    vout.setncattr('valid_max',90.0)
    vout.setncattr('valid_min',-90.0)
    vout[:,:] = pole[yn:ys,xw:xe]
    vout.setncattr('units','degrees_north')
    vout.set_auto_mask(True)
    dsout.sync()
    
###################################################
# #
# #   proc saves solar zenith angle cos into NC4 file
# #
def saveCosZenits(res,dsout,tplcrop):
    global gvars
    xw,xe,yn,ys=tplcrop
    pole = gvars.coszenits[res]
    typ = pole.dtype
    vout = dsout.createVariable('cosSolZenith',typ,dimensions=('y','x'), compression = 'zlib' )
    vout[:,:] = pole[yn:ys,xw:xe]
    dsout.sync()
    
###################################################
# #
# #   proc saves pixel scanning time array into NC4 file
# #
def saveTimes(res,dsout,tplcrop):
    global gvars
    xw,xe,yn,ys=tplcrop
    pole = gvars.times[res]
    typ = pole.dtype
    vout = dsout.createVariable('time',typ,dimensions=('y','x'), compression = 'zlib' )
    vout[:,:] = pole[yn:ys,xw:xe]
    dsout.sync()
    
###################################################
# #
# #   proc saves creates an AUX file and stores required  arrays there
# #
def saveAUX(fname):
    global gvars    

    lonlat = gvars.lonlat.get() if gvars.gui else gvars.lonlat
    costheta = gvars.costheta.get() if gvars.gui else gvars.costheta
    timechan = gvars.timechan.get() if gvars.gui else gvars.timechan
    for x in gvars.longitudes:
        res = str(x)
        tplcrop = cropForDim(int(x))
        xw,xe,yn,ys = tplcrop
        xv = xe-xw
        yv = ys-yn
        curfname = fname.replace('CH','aux'+'_'+res)
        testpath = os.path.dirname(curfname)
        if not os.path.exists(testpath):
            os.makedirs(testpath)
        dsout = Dataset(curfname, "w", format="NETCDF4")
        dsout.createDimension('y', yv)
        dsout.createDimension('x', xv)
        if (lonlat == 'aux'):
            print('saving LONLAT for ', res, ' into',curfname)
            saveGeo(x,dsout,tplcrop)
            dsout.sync()
        if (timechan == 'aux'):
            print('saving TIME for ', res, ' into',curfname)
            saveTimes(x,dsout,tplcrop)
            dsout.sync()
        if (costheta == 'aux'):
            saveCosZenits(x,dsout,tplcrop)
            dsout.sync()
        dsout.close()
    gvars.longitudes.clear()
    gvars.latitudes.clear()
    gvars.coszenits.clear()

########################################################
# #
# #  gets starting and ending row for selected crop
# #  / this is needed for selecting chunks that are to be read and processed
# #
def croprowsForChannel(chname):
    global gvars
    tpl = gvars.crops[gvars.selectedcrop.get()] if gvars.gui else gvars.crops[gvars.selectedcropname]
    
    yn = tpl[2]
    ys = tpl[3]
    dm = gvars.channeldims[chname][0]
    yn = (yn*dm) // gvars.mindimension 
    ys = (ys*dm) // gvars.mindimension
    yn = dm-yn+2
    ys = dm-ys+2 
    gvars.croprows=(yn,ys)

###########################################################
# #
# #   this recalculates the crop limits regarding to the required resolution
# 
def cropForDim(dm):
    global gvars
    tpl = gvars.crops[gvars.selectedcrop.get()] if gvars.gui else gvars.crops[gvars.selectedcropname]
    xw = tpl[0]
    xe = tpl[1]
    yn = tpl[2]
    ys = tpl[3]
    xw = (xw*dm) // gvars.mindimension 
    xe = (xe*dm) // gvars.mindimension
    yn = (yn*dm) // gvars.mindimension 
    ys = (ys*dm) // gvars.mindimension
    return(xw,xe,yn,ys)

########################################################
# #
# # this is the main process to decode and store selected channels etc.
# #
def doProcess(pslot = None):
    global gvars
    ##
    # # if GUI is running, get values from GUI
    if gvars.gui:
        gvars.win.withdraw()
        channels_vars = {i:gvars.channels_vars[i].get() for i in gvars.channels_vars}
        gvars.jsondata['channels']=channels_vars
        division = gvars.filedistribution.get()
        gvars.jsondata['division']=division
        lonlat = gvars.lonlat.get()
        gvars.jsondata['lonlat']=lonlat
        costheta = gvars.costheta.get()
        gvars.jsondata['cossolar']=costheta
        timechan = gvars.timechan.get()
        gvars.jsondata['timechan']=timechan
        outputformat = gvars.outputformat.get()
        gvars.jsondata['output']=outputformat
        vyr = gvars.selectedcrop.get()
        gvars.jsondata['crop']=vyr
        datum = gvars.datum
        print('datum = ',datum)
        slot = gvars.slot
        print('slot = ',slot)
        starttime = getStartTimeForSlot(slot)
        res = gvars.resolution
        jsonname = gvars.jsonnamebase + '_'+res+'.cfg'
        jsontxt = json.dumps(gvars.jsondata, indent=4, ensure_ascii=False)
        with open(gvars.apppath+'/'+jsonname,'w',encoding='utf-8') as out:
            out.write(jsontxt)
    else:
        ##   if running as commandline, get values from JSON file
        res = gvars.resolution
        jsonname = gvars.apppath+'/'+gvars.jsonnamebase + '_'+res+'.cfg'
        print('jsonname:',jsonname)
        with open(jsonname,'r') as fl:
            x = fl.read()
            print('x:')
            print(x)
        gvars.jsondata= json.loads(x)
        channels_vars = gvars.jsondata['channels']
        division = gvars.jsondata['division']
        lonlat = gvars.jsondata['lonlat']
        gvars.lonlat = lonlat
        costheta = gvars.jsondata['cossolar']
        gvars.costheta = costheta
        timechan = gvars.jsondata['timechan']
        gvars.timechan = timechan
        outputformat = gvars.jsondata['output']
        ## crop is found either in gvars.crops or in gvars.cropcode
        cropfromjson = gvars.jsondata['crop']
        if cropfromjson in gvars.cropcode:
            cropfromjson = gvars.cropcode[cropfromjson]
        if cropfromjson in gvars.crops:
            gvars.selectedcropname = cropfromjson
        else:
            print('parameters for crop named "'+cropfromjson+'" not found')
            exit()
        datum = gvars.datum
        slot = gvars.slot
        starttime = getStartTimeForSlot(slot)

    print('res = ', res)
    resolution = 'HI' if res=='HRFI' else 'LOW'
    getFileList(datum,slot,res)  # this fills gvars.filelist with list of filenames
                                 # for given combination of date, slot and resolution
    getChunkRows() # creates database with line numbers for files
        # from gvars.filelist and all channels.
        # Result is in gvars.startrow[chunk,ch]
        #       and in gvars.endrow[chunk,ch] 
    #
    # check if some available chunks get into crop 
    shallprocess = False
    print('channels_vars:', channels_vars)
    chans = list(channels_vars.keys())
    print('channels:', chans)
    chname = chans[0]
    print('columns:', gvars.pdstartrow.columns)
    croprowsForChannel(chname)
    yn,ys = gvars.croprows
    print('crop lines:', yn,ys)
    for ind in gvars.pdstartrow.index:
        strow = gvars.pdstartrow.at[ind,chname]
        endrow = gvars.pdendrow.at[ind,chname]
        if not(math.isnan(strow) or math.isnan(endrow)):
            #print('index:',strow,endrow)
            if (endrow>=ys) and (strow<=yn):
                shallprocess = True
    if not shallprocess:
        print('=========================================================================')
        print('Slot',slot,'does not contain the required data, it will not be processed.')
        print('=========================================================================')
        return
    fname = outputformat.replace('DT',datum)
    fname=fname.replace('SLOT',slot)
    fname=fname.replace('TM',starttime)
    fname=fname.replace('RES',resolution)
    print('=========================================')
    print('channels_vars:')
    for i in channels_vars:
        print(i,'  ',channels_vars[i])
    print('-------------------------')
    print('filedistribution:',division)
    print('lonlat:', lonlat) 
    print("costheta:", costheta)
    print('outputformat:', outputformat)
    print('datum = ', datum)
    print('slot = ', slot)
    print('starttime = ', starttime)
    print('res = ', res)
    print('resolution = ', resolution)
    print('fname = ',fname)
    print('=========================================')
    
    if division == 'every':   # every channel into its own file
        tplcrop = None
        for chname in channels_vars:
            listofmade = []
            wasmade = False
            gvars.ncattr.clear()
            storeStaticAttrs()
            croprowsForChannel(chname)
            x,y = gvars.channeldims[chname][:2]
            tplcrop = cropForDim(x)
            print('tplcrop =', tplcrop)
            # crop size
            xv=tplcrop[1]-tplcrop[0]
            yv=tplcrop[3]-tplcrop[2]
            print('channel:',chname,'  x:', xv, '  y:',yv)
            if channels_vars[chname]==1:
                #print('    processing channel ', chname )
                if chname.startswith('ir_') or chname.startswith('wv_'):
                    wasmade = getThermalForChannel(chname)
                if chname.startswith('vis_') or chname.startswith('nir_'):
                    wasmade = getVisForChannel(chname)
                    
                if wasmade:
                    listofmade.append(chname)
                            
                if wasmade:  # if channel was processed, it is stored immediately
                    curfname = fname.replace('CH',chname)
                    testpath = os.path.dirname(curfname)
                    if not os.path.exists(testpath):
                        os.makedirs(testpath)
                    dsout = Dataset(curfname, "w", format="NETCDF4")
                    dsout.createDimension('y', int(yv))
                    dsout.createDimension('x', int(xv))
                    dsout.setncattr('solar_distance_AU', gvars.sundist)
    
                    print('saving channel ', chname, 'into file ',curfname)
                    print('tplcrop =', tplcrop)
                    saveChannel(chname,dsout,tplcrop)         
                    dsout.setncattr('bands_available', ' '.join(listofmade))           
                    if lonlat=='every':
                        print('saving lonlat into file ',curfname)
                        listofmade.append('longitude')
                        listofmade.append('latitude')
                        saveGeo(x,dsout,tplcrop)
                        dsout.setncattr('lonlat_available','1')
                    else:
                        dsout.setncattr('lonlat_available','0')
                    if timechan=='every':
                        print('saving time info into file ',curfname)
                        listofmade.append('timechan')
                        saveTimes(x,dsout,tplcrop)
                        dsout.setncattr('timechan_available','1')
                    else:
                        dsout.setncattr('timechan_available','0')
                    if (costheta == 'every') or (
                        (costheta == 'vis') and (
                        (chname.startswith('vis_'))or(chname.startswith('nir_')))):
                        print('saving cos solar zenit into file ',curfname)
                        saveCosZenits(x,dsout,tplcrop)
                        dsout.setncattr('cosSolarZen_available','1')
                    else:
                        dsout.setncattr('cosSolarZen_available','0')
                    print('-------')
                    dsout.sync()
                    dsout.close()
                    gvars.chandata.pop(chname)    # smaznou se data kanalu, aby zbytecne nezabirala misto v pameti            
        if wasmade:
            if (lonlat=='aux' or costheta=='aux') and (tplcrop is not None):
                saveAUX(fname, tplcrop)
                
                
    if division == 'typ':   # one file for every type IR, WV, VIS, NIR
        for chtype in ['ir_','wv_','vis_','nir_']:
            # first determine number of selected channels of this type
            somethingwasmade = False
            gvars.ncattr.clear()
            storeStaticAttrs()
            llist = [i for i in channels_vars if (i.startswith(chtype)) and (channels_vars[i]==1)]
            cnt = len(llist)
            tplcrop = None
            listofmade = []
            if cnt>0:
                print('chtype:', chtype)
                print(llist)
                for chname in llist:
                    croprowsForChannel(chname)
                    x,y = gvars.channeldims[chname][:2]
                    print('    processing channel ', chname )
                    if chname.startswith('ir_') or chname.startswith('wv_'):
                        wasmade = getThermalForChannel(chname)
                    if chname.startswith('vis_') or chname.startswith('nir_'):
                        wasmade = getVisForChannel(chname)
                    if wasmade:
                        listofmade.append(chname)
                if len(listofmade)>0:
                    curfname = fname.replace('CH',chtype[:-1])
                    tplcrop = cropForDim(x)
                    # crop size
                    xv=tplcrop[1]-tplcrop[0]
                    yv=tplcrop[3]-tplcrop[2]
                    testpath = os.path.dirname(curfname)
                    if not os.path.exists(testpath):
                        os.makedirs(testpath)
                    dsout = Dataset(curfname, "w", format="NETCDF4")
                    saveAttrs(dsout)
                    dsout.createDimension('y', int(yv))
                    dsout.createDimension('x', int(xv))
                    #print('solar distance =',gvars.sundist)
                    dsout.solar_distance_AU = gvars.sundist
                    for chname in listofmade:
                        print('saving channel', chname, 'into file',curfname)
                        saveChannel(chname,dsout,tplcrop)
                    dsout.setncattr('bands_available', ' '.join(listofmade))           
                    somethingwasmade = True

                if len(listofmade)>0:
                    if lonlat=='every':
                        print('saving lonlat into file ',curfname)
                        saveGeo(x,dsout,tplcrop)
                        dsout.setncattr('lonlat_available','1')
                    else:
                        dsout.setncattr('lonlat_available','0')
                    if timechan=='every':
                        print('saving time info into file ',curfname)
                        saveTimes(x,dsout,tplcrop)
                        dsout.setncattr('timechan_available','1')
                    else:
                        dsout.setncattr('timechan_available','0')
                    if (costheta == 'every') or ((costheta == 'vis') and ((chtype =='vis_')or(chtype =='nir_'))):
                        print('saving costheta into file',curfname)
                        saveCosZenits(x,dsout,tplcrop)
                        dsout.setncattr('cosSolarZen_available','1')
                    else:
                        dsout.setncattr('cosSolarZen_available','0')
                    dsout.close()
        if somethingwasmade:
            if (lonlat=='aux' or costheta=='aux'):
                saveAUX(fname)
                
    if division == 'group':      # one file for every group (IR,WV), (VIS,NIR)
        for group in [['ir_','wv_'],['vis_','nir_']]:
            fulllist = []
            somethingwasmade = False
            print('group:', group)
            listofmade = []
            for chtype in group:
                wasmade = False
                llist = [i for i in channels_vars if (i.startswith(chtype)) and (channels_vars[i]==1)]
                fulllist.extend(llist)
                cnt = len(llist)
                if cnt>0:
                    print('    chtype:', chtype)
                    print('    llist:', llist)
                    for chname in llist:
                        croprowsForChannel(chname)
                        x,y = gvars.channeldims[chname][:2]
                        print('        processing channel ', chname )
                        if chname.startswith('ir_') or chname.startswith('wv_'):
                            wasmade = getThermalForChannel(chname) or wasmade
                        if chname.startswith('vis_') or chname.startswith('nir_'):
                            wasmade = getVisForChannel(chname) or wasmade
                        if wasmade:
                            listofmade.append(chname)
            if len(listofmade)>0:
                curfname = fname.replace('CH',group[0][:-1]+group[1][:-1])
                tplcrop = cropForDim(x)
                # crop size
                xv=tplcrop[1]-tplcrop[0]
                yv=tplcrop[3]-tplcrop[2]
                testpath = os.path.dirname(curfname)
                if not os.path.exists(testpath):
                    os.makedirs(testpath)
                dsout = Dataset(curfname, "w", format="NETCDF4")
                dsout.createDimension('y', int(yv))
                dsout.createDimension('x', int(xv))
                dsout.setncattr('solar_distance_AU', gvars.sundist)
                for chname in listofmade:
                    print('saving channel', chname, 'into file',curfname)
                    saveChannel(chname,dsout, tplcrop)
                dsout.setncattr('bands_available', ' '.join(listofmade))           
                if lonlat=='every':
                    print('saving lonlat into file ',curfname)
                    saveGeo(x,dsout,tplcrop)
                    dsout.setncattr('lonlat_available','1')
                else:
                    dsout.setncattr('lonlat_available','0')
                if timechan=='every':
                    print('saving time info into file ',curfname)
                    saveTimes(x,dsout,tplcrop)
                    dsout.setncattr('timechan_available','1')
                else:
                    dsout.setncattr('timechan_available','0')
                if (costheta == 'every') or ((costheta == 'vis') and ('vis_' in group)):
                    print('saving cosSolarZen into file',curfname)
                    saveCosZenits(x,dsout,tplcrop)
                    dsout.setncattr('cosSolarZen_available','1')
                else:
                    dsout.setncattr('cosSolarZen_available','0')
                dsout.close()
        if lonlat=='aux' or costheta=='aux':
            saveAUX(fname)
    gvars.chandata.clear()
    gvars.longitudes.clear()
    gvars.latitudes.clear()
    gvars.coszenits.clear()
    print('##############################################')
    print('FINISHED')
    if gvars.gui:
        gvars.win.deiconify()

#############################################################
# #
# #  processes all available slots in selected date
def processDatum():
    global gvars
    
    datum = gvars.datum
    print('datum = ',datum)
    res = gvars.resolution
    print('res = ', res)
    dateslots = getCyclesForDateAndRes(datum,res)
    for s in dateslots:
        gvars.slot = s
        print('*********************************')
        print('                 SLOT   ',s)
        print('*********************************')
        doProcess()
        
#############################################################
# #
# #  processes slots written in the editbox for selected date
def processSlots():
    global gvars
    selectedslots = gvars.slotselection.get()
    txtslots = decodeslotlist(selectedslots)
    res = gvars.resolution
    #print('txtsloty =',txtsloty)
    existing = [x[0] for x in gvars.slots_chunks if x[2]==res]
    doproc = [x for x in txtslots if x in existing]
    print('slots to process =',doproc)
    for s in doproc:
        gvars.slot = s
        print('*********************************')
        print('                 SLOT   ',s)
        print('*********************************')
        doProcess()
        
        
    
##############################################################
# #
# #  proc reads the current working folder and finds all dates
# #  for which there are MTG files there
# #
def readFolder():
    global gvars

    seznam = glob.glob(gvars.ncmask)
    
    gvars.folderdates = []
    gvars.wordlist = [x[:-3].split('_') for x in seznam]  # cut off the extension ".nc"
    for x in gvars.wordlist:
        x[1]= x[1].split('-')
    gvars.folderdates = list(set([x[gvars.istarttime][:8] for x in gvars.wordlist]))
    gvars.folderdates.sort()

    # if GUI is running, the listboxes are emptied
    if gvars.gui:
        gvars.lbox_date.delete(0,tk.END)
        gvars.lbox_slot.delete(0,tk.END)
        #  now the listbox containing dates is filled
        for d in gvars.folderdates:
            gvars.lbox_date.insert(tk.END,d)
        
        
    
############################################################
# #
# #    proc opens the "select folder" dialog and tries to
# #    change the working dir into the result
# # 
def setFolder():
  global gvars
  res = fd.askdirectory(initialdir = '.', mustexist = tk.TRUE)
  if res != '':
      os.chdir(res)
      readFolder()

###########################################################
# #
# #  proc gets the slot list written into either the editbox
# #  or from commandline and converts it to classic list without intervals
def decodeslotlist(parslot):
    lslotlist = []
    try:
        parslot = parslot.replace(' ','')   # get rid of spaces
        parslotlist = parslot.split(',')
        print('parslotlist:',parslotlist)
        for prvek in parslotlist:
            border = prvek.split('-')
            print('border:',border)
            iborder = [int(x) for x in border]
            print('iborder:',iborder)
            if len(iborder)>1:     # if this is interval
                print('interval')
                for x in range(iborder[0],iborder[1]+1):
                    if (x >0) and (x<= 144):
                        lslotlist.append(x)
            else:                  # if single slot
                print('single')
                if iborder[0]>0 and iborder[0]<=144:
                    lslotlist.append(iborder[0])
            print('lslotlist:',lslotlist)
        lslotlist = list(set(lslotlist))   # only unique values
        lslotlist.sort()
        txtslotlist = ['{:04d}'.format(x) for x in lslotlist]
    except:
        print('incorrect format found')
        txtslotlist = []
    return txtslotlist



def onexit():
    global gvars
    os.chdir(gvars.origdir)
    gvars.win.destroy()                  


##########################################################################
##########################################################################
##                                                                      ##
##                                                                      ##
##                                 MAIN                                 ##
##                                                                      ##
##                                                                      ##
##########################################################################
##########################################################################



n=len(sys.argv)
print('number of arguments:',n)
print('argument 0: ', sys.argv[0])


gvars.jsonnamebase = sys.argv[1] if n>1 else 'mtg2snap'  # first argument is the base for JSON filename
gvars.gui = False if n>2 else True                       # next arguments are processing parameters
                                                         # datum, slot/s, resolution

np.seterr(invalid='ignore', divide='ignore')
warnings.filterwarnings('ignore','divide')
warnings.filterwarnings('ignore','invalid value')


gvars.origdir = os.getcwd()

appname = sys.argv[0]
head,tail = os.path.split(appname)
gvars.apppath =head
print('head =', head)
print('tail =', tail)

## load the JSON file with crop definitions
# 
jsonname = gvars.apppath+'/'+'mtg2snap_crops.json'
print('jsonname:',jsonname)
with open(jsonname,'r') as fl:
    x = fl.read()
    print('x:')
    print(x)
jsondata= json.loads(x)
gvars.crops = jsondata['crops']
gvars.cropcode = jsondata['cropcode']




######################################
## 
##  vytvoreni tkinter GUI, pokud je to pozadovano
##

if gvars.gui: 
    gvars.win = tk.Tk()
    
    listfont = tkfont.Font(family = 'Lucida Sans Typewriter', size = 12)
    btnfont = tkfont.Font(family = 'Lucida Sans Typewriter', size = 10)
    boldfont = tkfont.Font(family='Helvetica', size=10, weight='bold')
    slotcolor = '#ffc'
    
    fdatum = tk.Frame(gvars.win, bd = 5, relief = tk.RIDGE)
    btn = tk.Button(fdatum, text='Change the folder', command=setFolder)
    btn.pack(side = tk.BOTTOM, fill=tk.X)
    gvars.lbox_date = tk.Listbox(fdatum,selectmode = tk.SINGLE, 
                                  exportselection = 0,
                                  width = 10, font=listfont)
    gvars.scrbar_date = tk.Scrollbar(
            fdatum,
            orient=tk.VERTICAL,
            command=gvars.lbox_date.yview
        )
    gvars.lbox_date['yscrollcommand'] = gvars.scrbar_date.set
    gvars.lbox_date.pack(side=tk.LEFT, fill = tk.Y)
    gvars.scrbar_date.pack(side=tk.LEFT, fill=tk.Y)
    fdatum.pack(side=tk.LEFT, fill = tk.Y)
    
    fslot = tk.Frame(gvars.win, bd = 5, relief = tk.RIDGE)
    gvars.lbox_slot = tk.Listbox(fslot,selectmode = tk.SINGLE, 
                                 exportselection=0,
                                 width = 30, font=listfont)
    gvars.scrbar_slot = tk.Scrollbar(
            fslot,
            orient=tk.VERTICAL,
            command=gvars.lbox_slot.yview
        )
    gvars.lbox_slot['yscrollcommand'] = gvars.scrbar_slot.set
    gvars.lbox_slot.pack(side=tk.LEFT, fill = tk.Y)
    gvars.scrbar_slot.pack(side=tk.LEFT, fill=tk.Y)
    fslot.pack(side=tk.LEFT, fill = tk.Y)
    
    gvars.lbox_date.bind('<<ListboxSelect>>',datechange)
    gvars.lbox_slot.bind('<<ListboxSelect>>',changeTheSlot)
    
    # frame for displaying of selected date and slot, channel selection and other actions
    gvars.fselection = tk.Frame(gvars.win )
    gvars.fselection.pack(side=tk.LEFT, fill=tk.Y)
    # datums and slots
    fdatslot = tk.Frame(gvars.fselection,bd = 5, relief = tk.RIDGE, padx = 10, bg = slotcolor)
    fdatslot.pack(side = tk.TOP, fill=tk.X)
    gvars.lbldate = tk.Label(fdatslot, text = '', font = btnfont, padx = 5, bg = slotcolor)
    gvars.lblslot = tk.Label(fdatslot, text = '', font = btnfont, padx = 5, bg = slotcolor)
    gvars.lbldate.pack(side=tk.LEFT)
    gvars.lblslot.pack(side=tk.LEFT)
    fdatslot.pack(side = tk.TOP, fill=tk.X)
    
    # frame for storage ways and others
    faux = tk.Frame(gvars.fselection, bd = 5,relief=tk.RIDGE)
    faux.pack(side=tk.BOTTOM, fill=tk.X)
    
    # how everything will be stored into files
    gvars.filedistribution = tk.StringVar()
    lblfiles = tk.Label(faux, text = 'how to store:', font=boldfont)
    rbfileevery = tk.Radiobutton(faux,text = 'every channel one file', variable = gvars.filedistribution, value = 'every')
    rbfiletyp = tk.Radiobutton(faux,text = 'vis,nir,ir,wv', variable = gvars.filedistribution, value = 'typ')
    rbfilegroup = tk.Radiobutton(faux,text = '(vis+nir),(ir+wv)', variable = gvars.filedistribution, value = 'group')
    gvars.filedistribution.set('group')
    lblfiles.grid(in_ = faux, row=0, column=0, columnspan=3, sticky=tk.W)
    rbfileevery.grid(in_ = faux, row = 1, column = 0, sticky = tk.W)
    rbfiletyp.grid(in_ = faux, row = 1, column = 1, sticky = tk.W)
    rbfilegroup.grid(in_ = faux, row = 1, column = 2, sticky = tk.W)
    
    # what will happen with LONLAT
    gvars.lonlat = tk.StringVar()
    lbllonlat = tk.Label(faux, text = 'LON and LAT:', font=boldfont)
    rblonlatevery = tk.Radiobutton(faux,text = 'every file', variable = gvars.lonlat, value = 'every')
    rblonlataux = tk.Radiobutton(faux,text = 'special AUX file', variable = gvars.lonlat, value = 'aux')
    rblonlatnone = tk.Radiobutton(faux,text = 'none', variable = gvars.lonlat, value = 'none')
    gvars.lonlat.set('aux')
    lbllonlat.grid(in_ = faux, row=2, column=0, columnspan=3, sticky=tk.W)
    rblonlatevery.grid(in_ = faux, row = 3, column = 0, sticky = tk.W)
    rblonlataux.grid(in_ = faux, row = 3, column = 2, sticky = tk.W)
    rblonlatnone.grid(in_ = faux, row = 3, column = 3, sticky = tk.W)
    
    
    # what will happen with cos of solar zenith angle
    gvars.costheta = tk.StringVar()
    lbltheta = tk.Label(faux, text = 'cos of solar zenith angle:', font=boldfont)
    rbthetaevery = tk.Radiobutton(faux,text = 'every file', variable = gvars.costheta, value = 'every')
    rbthetavis = tk.Radiobutton(faux,text = 'every VIS and NIR', variable = gvars.costheta, value = 'vis')
    rbthetaaux = tk.Radiobutton(faux,text = 'special AUX file', variable = gvars.costheta, value = 'aux')
    rbthetanone = tk.Radiobutton(faux,text = 'none', variable = gvars.costheta, value = 'none')
    gvars.costheta.set('aux')
    lbltheta.grid(in_ = faux, row=4, column=0, columnspan=3, sticky=tk.W)
    rbthetaevery.grid(in_ = faux, row = 5, column = 0, sticky = tk.W)
    rbthetavis.grid(in_ = faux, row = 5, column = 1, sticky = tk.W)
    rbthetaaux.grid(in_ = faux, row = 5, column = 2, sticky = tk.W)
    rbthetanone.grid(in_ = faux, row = 5, column = 3, sticky = tk.W)
    
    # co bude s daty timechannel - pixel scanning time
    gvars.timechan = tk.StringVar()
    lblTimechan = tk.Label(faux, text = 'pixel scanning time:', font=boldfont)
    rbTimechanevery = tk.Radiobutton(faux,text = 'every file', variable = gvars.timechan, value = 'every')
    rbTimechanaux = tk.Radiobutton(faux,text = 'special AUX file', variable = gvars.timechan, value = 'aux')
    rbTimechannone = tk.Radiobutton(faux,text = 'none', variable = gvars.timechan, value = 'none')
    gvars.timechan.set('aux')
    lblTimechan.grid(in_ = faux, row=6, column=0, columnspan=3, sticky=tk.W)
    rbTimechanevery.grid(in_ = faux, row = 7, column = 0, sticky = tk.W)
    rbTimechanaux.grid(in_ = faux, row = 7, column = 2, sticky = tk.W)
    rbTimechannone.grid(in_ = faux, row = 7, column = 3, sticky = tk.W)
    
    # path and format of output fie name:
    # 
    # DT = datum YYYYMMDD
    # TM = nominal time of start HHmm  - step by 10 minutes
    # SLOT = slot number 0001-0144
    # CH = if single channel per file, then channel name
    #    = if group ir, wv, vis, nir nebo irwv, visnir per file, then this group
    #    = aux for file with additional info
    # RES = resolution LOWRES, HIRES 
    
    entrylabel = tk.Label(faux, text='output path and filename:', font=boldfont)
    entrylabel.grid(in_ = faux,row = 8, column = 0, columnspan=3, sticky=tk.W)
    gvars.outputformat = tk.StringVar(faux,'./DT-SLOT-TM-CH-RES.nc')
    entryformat = tk.Entry(faux, textvariable=gvars.outputformat, font=boldfont)
    entryformat.grid(in_ = faux,row = 9, column = 0, columnspan=3, sticky=tk.W+tk.E)
    txt = 'Formatting codes:\n\n'+ \
        'DT - datum in format YYYYMMDD\n' + \
        'TM - nominal scanning starttime in format HHmm\n' + \
        'SLOT -  slot number 0001-0144\n' + \
        'CH - channel or grop name\n' + \
        '     (vis_09, wv_78, ir, nir, irwv, visnir, ...)\n' + \
        'RES - resolution LOW, HI'
    infolabel = tk.Label(faux, text=txt, anchor=tk.NW, justify = tk.LEFT)
    infolabel.grid(in_=faux,row=10, column=0, rowspan = 6, columnspan=3, stick = tk.N+tk.S+tk.E+tk.W)
    
    #  crop selection
    lblcrop = tk.Label(gvars.win,text = 'select crop:')
    lblcrop.pack(side = tk.TOP)
    optionlist = list(gvars.crops.keys())
    print('optionlist =', optionlist)
    gvars.selectedcrop=tk.StringVar()
    gvars.selectedcrop.set(optionlist[0])
    mnucrop = tk.OptionMenu(gvars.win,gvars.selectedcrop, *optionlist)
    mnucrop.pack(side = tk.TOP, fill=tk.X)
    
    btndelejdatum = tk.Button(gvars.win, text = ' process full day ', command = processDatum)
    btndelejdatum.pack(side = tk.BOTTOM)
    btndelejsloty = tk.Button(gvars.win, text = ' process written slots ', command = processSlots)
    btndelejsloty.pack(side = tk.BOTTOM)
    gvars.slotselection = tk.StringVar()
    gvars.slotselection.set('')
    entrysloty = tk.Entry(gvars.win, textvariable=gvars.slotselection, font=boldfont)
    entrysloty.pack(side = tk.BOTTOM, fill = tk.X)
    btndelej = tk.Button(gvars.win, text = 'Process', command = doProcess)
    btndelej.pack(side = tk.BOTTOM)
    
    
    
    readFolder()
    gvars.win.protocol("WM_DELETE_WINDOW", onexit)
    gvars.win.mainloop()    

else:      ##   NO GUI   ###########################################
    
    parametry = sys.argv[2:]
    print('parameters:')
    print(parametry)
    
    readFolder()
    
    # now parameters are processed in triplets
    # 
    for pardatum,parslot,parres in zip(parametry[0::3],parametry[1::3],parametry[2::3]):
        #
        # resolution
        # 
        parres = parres.upper()
        if parres not in ['FDHSI', 'HRFI', 'ALL', 'BOTH']:
            print('resolution',parres,'is invalid, go on for next combination')
            continue
        if parres in ['ALL', 'BOTH']:
            ress = ['FDHSI', 'HRFI']
        else:
            ress = [parres]   # in "ress" there is the list of resolutions that are to be processed
        #
        # SLOT
        #       may be given as single number,
        #       as list with commas 
        #       nebo jako interval s pomlckou  or as interval a-b
        # 
        #       the list with commas is the main structure,
        #       its elements may be either single values or intervals
        #
        txtslotlist = decodeslotlist(parslot)
        gvars.datum = pardatum
        for myres in ress:
            gvars.resolution = myres
            slotsavail = getCyclesForDateAndRes(pardatum, myres)        
            for myslot in txtslotlist:
                if myslot in slotsavail:
                    gvars.slot=myslot
                    getFileList(pardatum,myslot,myres)  # this fills "gvars.filelist", which is 
                                                        # list of files for given combination
                    getAllDims(pardatum,myslot,myres)  # here dimensions of all channels are determined
                        # dimensions are stored in gvars.channeldims[ch]
                        # as a tuple (cols,rows,startrow,endrow)
                    # getChunkRows() # creates database of start+endlines for chunks
                        # from gvars.filelist and all channels
                        # Results are in the pandas dataframes gvars.startrow[chunk,ch]
                        #                                  and gvars.endrow[chunk,ch] 
                    doProcess()
                else:
                    print('--------------------------------------------------------------')
                    print('slot', myslot, 'is not available for given date and resolution')
                    print('--------------------------------------------------------------')
                
    
    