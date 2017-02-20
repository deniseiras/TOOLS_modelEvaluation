#!usr/bin/python
# -*- coding: UTF8 -*-
import numpy as np
from math import *
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import datetime
from netCDF4 import Dataset as Dataset, date2index, num2date
import argparse

model = [None] * 2
mcase = [None] * 2
scase = [None] * 2
varname = [None] * 2
level = [None] * 2
sdd = [None] * 2
smm = [None] * 2
syy = [None] * 2
shh = [None] * 2
sddForecast = [None] * 2
smmForecast = [None] * 2
syyForecast = [None] * 2
shhForecast = [None] * 2


def getBaseDirectory():
    aBaseDirectory = '/home2/denis/new_aerosols'
    return aBaseDirectory


def getDirectory(aBaseDirectory, amodel, amcase, ahh):
    aDirectory = "{}/{}/{}/r{}/".format(aBaseDirectory, amodel, amcase, ahh)
    return aDirectory


def getFileName(amodel, amcase, ascase, ayy, amm, add, ahh):
    aBaseDirectory = getBaseDirectory()
    aDirectory = getDirectory(aBaseDirectory, amodel, amcase, ahh)
    aFileName = "{}{}_{}_{}_{}{}{}0000.nc".format(aDirectory, amodel, amcase, ascase, ayy, amm, add, ahh)
    return aFileName


def openNetCdfFile(aFileName):
    try:
        print("Opening file ... " + aFileName)
        ncFile = Dataset(aFileName, 'r')
    except:
        print("Error opening file ", aFileName)
    return ncFile


def plotOrSaveFigure(isPlot, strFile, plt):
    if isPlot is False:
        print("Generating figure ... {}".format(strFile))
        plt.savefig(strFile, bbox_inches='tight')
    else:
        print("Viewing figure ... ")
        plt.show()


def drawMapLines(m):
    m.drawcoastlines()
    # m.fillcontinents()
    m.drawmapboundary()
    m.drawparallels(np.arange(-90., 120., 30.), labels=[1, 0, 0, 0])
    m.drawmeridians(np.arange(-180., 180., 60.), labels=[0, 0, 0, 1])


def generateFigure(saveDirectory, isPlot, amodel, amcase, ascase, avarname, level, add, amm, ayy, ahh, addf, ammf, ayyf, ahhf):

    dateInitial = datetime.datetime(int(ayy), int(amm), int(add), int(ahh))
    dateForecast = datetime.datetime(int(ayyf), int(ammf), int(addf), int(ahhf))

    aFileName = getFileName(amodel, amcase, ascase, ayy, amm, add, ahh)
    ncfile = openNetCdfFile(aFileName)

    aVariable = ncfile.variables[avarname]
    lat = ncfile.variables['lat'][:]
    lon = ncfile.variables['lon'][:]
    tim = ncfile.variables['time']

    idxDateIni = date2index(dateInitial, tim, calendar='standard')
    idxDateForecast = date2index(dateForecast, tim, calendar='standard')
    dates = num2date(tim[idxDateIni:idxDateForecast + 1], tim.units, calendar='standard')
    date = dates[idxDateForecast]

    if aVariable.ndim == 3:
        data = aVariable[idxDateForecast, :, :]
    elif aVariable.ndim == 4:
        data = aVariable[idxDateForecast, level, :, :]

    plt.figure(figsize=(20, 20), facecolor='white')
    m = Basemap(projection='mill', lat_ts=10, llcrnrlon=lon.min(),
                urcrnrlon=lon.max(), llcrnrlat=lat.min(), urcrnrlat=lat.max(),
                resolution='c')

    # convert the lat[i]/lon[i] values to x/y projections.
    x, y = m(*np.meshgrid(lon, lat))

    m.pcolormesh(x, y, data, shading='flat', cmap=plt.cm.jet)
    cbar = m.colorbar(location='bottom', pad="10%")

    if level is None:
        levelStr = ""
    else:
        levelStr = ", level = " + str(level)

    cbar.set_label(avarname + r"($" + aVariable.units + "$)" + levelStr)
    drawMapLines(m)

    strCaso = "Case: {}, {}, {}, on {}-{}-{} {}:00:00, forecast {} ".format(amodel, amcase, ascase, ayy, amm, add, ahh, str(object=date))
    plt.title(strCaso)
    strFile = "{}{}_{}_{}_{}_{}-{}-{}_{}:00:00_fore_{}-{}-{}_{}:00:00.png".format(saveDirectory, amodel, amcase, ascase, avarname, ayy, amm, add, ahh,
                                                                                  ayyf, ammf, addf, ahhf)
    plotOrSaveFigure(isPlot, strFile, plt)
    ncfile.close()


def generateComparison(saveDirectory, isPlot, strFunction, model, mcase, scase, varname, level, dd, mm, yy, hh, ddForecast, mmForecast, yyForecast, hhForecast):

    dateInitial = [None] * 2
    dateForecast = [None] * 2
    fileName = [None] * 2
    ncfile = [None] * 2
    variable = [None] * 2
    lat = [None] * 2
    lon = [None] * 2
    tim = [None] * 2
    idxDateIni = [None] * 2
    idxDateForecast = [None] * 2
    dates = [None] * 2
    data = [None] * 2
    date = [None] * 2

    if len(level) != 2:
        print("\n\n Error: -level and -level2 paramters are mandatory for comparison! ")
        exit()

    for i in [0, 1]:

        dateInitial[i] = datetime.datetime(yy[i], mm[i], dd[i], hh[i])
        dateForecast[i] = datetime.datetime(yyForecast[i], mmForecast[i], ddForecast[i], int(hhForecast[i]))

        fileName[i] = getFileName(model[i], mcase[i], scase[i], syy[i], smm[i], sdd[i], shh[i])
        ncfile[i] = openNetCdfFile(fileName[i])

        try:
            variable[i] = ncfile[i].variables[varname[i]]
        except Exception:
            print("\n\n Error: variable {} not found".format(varname[i]))
            exit()

        variable[i] = ncfile[i].variables[varname[i]]

        lat[i] = ncfile[i].variables['lat'][:]
        lon[i] = ncfile[i].variables['lon'][:]
        tim[i] = ncfile[i].variables['time']

        idxDateIni[i] = date2index(dateInitial[i], tim[i], calendar='standard')
        idxDateForecast[i] = date2index(dateForecast[i], tim[i], calendar='standard')
        dates[i] = num2date(tim[i][idxDateIni[i]:idxDateForecast[i] + 1], tim[i].units, calendar='standard')
        date[i] = dates[i][idxDateForecast[i]]

        # in case of one date, use the forecast date ...
        if variable[i].ndim == 3:
            data[i] = variable[i][idxDateForecast[i], :, :]
        elif variable[i].ndim == 4:
            data[i] = variable[i][idxDateForecast[i], level[i], :, :]

    i = 0
    plt.figure(figsize=(20, 20), facecolor='white')
    m = Basemap(projection='mill', lat_ts=10, llcrnrlon=lon[i].min(),
                urcrnrlon=lon[i].max(), llcrnrlat=lat[i].min(), urcrnrlat=lat[i].max(),
                resolution='c')

    x, y = m(*np.meshgrid(lon[i], lat[i]))

    datadiff = np.array(data)
    if strFunction == "bias":
        strFunctionTitle = "Bias"
        datadiff = datadiff[0] - datadiff[1]
    elif strFunction == "avg":
        strFunctionTitle = "Average"
        datadiff = np.mean(datadiff, axis=0)
    # elif strFunction == "rmse":
        # preciso de mais dados pra calcular as 2 medias
        # datadiff = sqrt(pow(media1-media2, 2))
    elif strFunction == "std":
        strFunctionTitle = "Standard deviation"
        datadiff = np.std(datadiff, axis=0)
    elif strFunction == "var":
        strFunctionTitle = "Variance"
        datadiff = np.var(datadiff, axis=0)

    m.pcolormesh(x, y, datadiff, shading='flat', cmap=plt.cm.jet)
    cbar = m.colorbar(location='bottom', pad="10%")

    varUnits = [" - "] * 2
    for uni in [0, 1]:
        if variable[uni].units != '':
            varUnits[uni] = variable[uni].units

    levelStr = [None] * 2
    if level is None:
        levelStr[i] = ""
        levelStr[i + 1] = ""
    else:
        levelStr[i] = ", level = " + str(level[i])
        levelStr[i + 1] = ", level = " + str(level[i + 1])

    strLabel = varname[i] + r"($" + varUnits[i] + r"$)" + levelStr[i] + " / " + varname[i + 1] + r"($" + varUnits[i + 1] + r"$)" + levelStr[i + 1]
    cbar.set_label(strLabel)
    drawMapLines(m)

    strCaso = strFunctionTitle + ":\nCase 1: {}, {}, {} {}-{}-{} {}:00:00, forecast {}-{}-{} {}:00:00\nCase 2: {}, {}, {} {}-{}-{} {}:00:00, forecast {}-{}-{} {}:00:00" \
        .format(model[i], mcase[i], scase[i], syy[i], smm[i], sdd[i], shh[i], syyForecast[i], smmForecast[i], sddForecast[i], shhForecast[i],
                model[i + 1], mcase[i + 1], scase[i + 1], syy[i + 1], smm[i + 1], sdd[i + 1], shh[i + 1], syyForecast[i + 1], smmForecast[i + 1], sddForecast[i + 1],
                shhForecast[i + 1])
    plt.title(strCaso)

    strFile = "{}{}_{}_{}_{}_{}_{}-{}-{}_{}:00:00_fore_{}-{}-{}_{}:00:00__x__{}_{}_{}_{}_{}-{}-{}_{}:00:00_fore_{}-{}-{}_{}:00:00.png" \
        .format(saveDirectory, strFunction, model[i], mcase[i], scase[i], varname[i], syy[i], smm[i], sdd[i], shh[i], syyForecast[i], smmForecast[i], sddForecast[i],
                shhForecast[i], model[i + 1], mcase[i + 1], scase[i + 1], varname[i + 1], syy[i + 1], smm[i + 1], sdd[i + 1], shh[i + 1], syyForecast[i + 1],
                smmForecast[i + 1], sddForecast[i + 1], shhForecast[i + 1])

    plotOrSaveFigure(isPlot, strFile, plt)
    for i in [0, 1]:
        ncfile[i].close()


# main =======================================

parser = argparse.ArgumentParser(description='Evaluate models')

parser.add_argument("-function", type=str, help='function: "bias" for Bias, "avg" for Average, "std" for Standard deviation, "var" for Variance: ')
parser.add_argument('-plot', dest='isPlot', action='store_true', help='Plot the first case figure on the screen')
parser.add_argument('-no-plot', dest='isPlot', action='store_false', help='Save the first case figure or function figure, whether specified, to a file')
parser.set_defaults(isPlot=True)

parser.add_argument("-model", type=str, required=True, help='model center file (ncep, nasa, ecmwf) of the first case')
parser.add_argument("-case", type=str, required=True, help='case (dust, smoke, pollution) of the first case')
parser.add_argument("-subcase", type=str, required=True, help='subcase (interactive, direct, indirect, noaerosols) of the first case')
parser.add_argument("-variable", type=str, required=True, help='2D or 3D variable (aod, dlwf, dswf, dustmass, temp2m, wdir, wmag) of the first case')
parser.add_argument("-level", type=str, help='level of the 3D variable (aod, dlwf, dswf, dustmass, temp2m, wdir, wmag) of the first case')
parser.add_argument("-year", type=str, required=True, help='year of the first case')
parser.add_argument("-month", type=str, required=True, help='month of the first case')
parser.add_argument("-day", type=str, required=True, help='day of the first case')
parser.add_argument("-hour", type=str, required=True, help='hour of the first case')
parser.add_argument("-foreyear", type=str, required=True, help=' forecast year of the first case')
parser.add_argument("-foremonth", type=str, required=True, help='forecast month of the first case')
parser.add_argument("-foreday", type=str, required=True, help='forecast day of the first case')
parser.add_argument("-forehour", type=str, required=True, help='forecast hour of the first case')

parser.add_argument("-model2", type=str, default='ncep', help='2 model center file (ncep, nasa, ecmwf) of the second case')
parser.add_argument("-case2", type=str, help='2 case (dust, smoke, pollution) of the second case')
parser.add_argument("-subcase2", type=str, help='2 subcase (interactive, direct, indirect, noaerosols) of the second case')
parser.add_argument("-variable2", type=str, help='2 2D variable (aod, dlwf, dswf, dustmass, temp2m, wdir, wmag) of the second case')
parser.add_argument("-level2", type=str, help='level of the 3D variable (aod, dlwf, dswf, dustmass, temp2m, wdir, wmag) of the first case')
parser.add_argument("-year2", type=str, help='year of the second case')
parser.add_argument("-month2", type=str, help='month of the second case')
parser.add_argument("-day2", type=str, help='day of the second case')
parser.add_argument("-hour2", type=str, help='hour of the second case')
parser.add_argument("-foreyear2", type=str, help='forecast year of the second case')
parser.add_argument("-foremonth2", type=str, help='forecast month of the second case')
parser.add_argument("-foreday2", type=str, help='forecast day of the second case')
parser.add_argument("-forehour2", type=str, help='forecast hour of the second case')

args = parser.parse_args()

saveDirectory = "./"
func = args.function
isPlot = args.isPlot

model[0] = args.model
mcase[0] = args.case
scase[0] = args.subcase
varname[0] = args.variable
level[0] = args.level
syy[0] = args.year
smm[0] = args.month
sdd[0] = args.day
shh[0] = args.hour
syyForecast[0] = args.foreyear
smmForecast[0] = args.foremonth
sddForecast[0] = args.foreday
shhForecast[0] = args. forehour

model[1] = args.model2
mcase[1] = args.case2
scase[1] = args.subcase2
varname[1] = args.variable2
level[1] = args.level2
syy[1] = args.year2
smm[1] = args.month2
sdd[1] = args.day2
shh[1] = args.hour2
syyForecast[1] = args.foreyear2
smmForecast[1] = args.foremonth2
sddForecast[1] = args.foreday2
shhForecast[1] = args.forehour2

yy = [int(i) for i in syy if i is not None]
mm = [int(i) for i in smm if i is not None]
dd = [int(i) for i in sdd if i is not None]
hh = [int(i) for i in shh if i is not None]
yyForecast = [int(i) for i in syyForecast if i is not None]
mmForecast = [int(i) for i in smmForecast if i is not None]
ddForecast = [int(i) for i in sddForecast if i is not None]
hhForecast = [int(i) for i in shhForecast if i is not None]

if func is None:
    generateFigure(saveDirectory, isPlot, model[0], mcase[0], scase[0], varname[0], level[0], sdd[0], smm[0], syy[0], shh[0], sddForecast[0], smmForecast[0], syyForecast[0], shhForecast[0])
    # generateFigure(isPlot, model[1], mcase[1], scase[1], varname[1], sdd[1], smm[1], syy[1], shh[1], sddForecast[1], smmForecast[1], syyForecast[1], shhForecast[1])
else:
    generateComparison(saveDirectory, isPlot, func, model, mcase, scase, varname, level, dd, mm, yy, hh, ddForecast, mmForecast, yyForecast, hhForecast)


# This is a simple example which reads a small dummy array, from a
# netCDF data file created by the companion program simple_xy_wr.py.

# This example demonstrates the netCDF Python API.
# It will work either with the Scientific Python NetCDF version 3 interface
# (http://dirac.cnrs-orleans.fr/ScientificPython/)
# (http://netcdf4-python.googlecode.com/svn/trunk/docs/netCDF4_classic-module.html)
# To switch from one to another, just comment/uncomment the appropriate
# import statements at the beginning of this file.

# Jeff Whitaker <jeffrey.s.whitaker@noaa.gov> 20070201

#  the Scientific Python netCDF 3 interface
#  http://dirac.cnrs-orleans.fr/ScientificPython/


# from Scientific.IO.NetCDF import NetCDFFile as Dataset


# the 'classic' version of the netCDF4 python interface
# http://code.google.com/p/netcdf4-python/
# from numpy import arange as # array module from http://numpy.scipy.org
# from numpy.testing import assert_array_equal, assert_array_almost_equal
