#!/usr/bin/env python

'''
EnvironmentalLoggerAnalyser.py

----------------------------------------------------------------------------------------
This module will read data generated by Environmental Sensor and convert to netCDF file
----------------------------------------------------------------------------------------
Prerequisite:
1. Python (2.7+ recommended)
2. netCDF4 module for Python (and many other supplements such as numpy, scipy and HDF5 if needed)
----------------------------------------------------------------------------------------

Usage (in bash):

python ${HOME}/terraref/computing-pipeline/scripts/hyperspectral/EnvironmentalLoggerAnalyser.py ${DATA}/terraref/input ${DATA}/terraref/output

Please note EnvironmentalLoggerAnalyser.py will take the second parameter as the input folder (containing JSON files,
but it can also be one single file) and the third parameter as the output folder (will dump netCDF files here).
If the output folder does not exist, EnvironmentalLoggerAnalyser.py will create one for you.

----------------------------------------------------------------------------------------
Update 4.29

The output JSON file is now completely composed by variables
2D spectrometer variables (wavelength and spectrum) are available in the exported file
----------------------------------------------------------------------------------------
Thanks for the advice from Professor Zender and testing data from Mr. Maloney
----------------------------------------------------------------------------------------
'''

import json
import time
import sys
import os
from netCDF4 import Dataset

_UNIT_DICTIONARY = {u'm': 'meter', u"hPa": "hecto-Pascal", u"DegCelsius": "Celsius",
                    u's': 'second', u'm/s': 'meter second-1', u"mm/h": 'milimeters hour-1',
                    u"relHumPerCent": "PerCent", u"?mol/(m^2*s)": "micromole meters-2 second-1",
                    u'kilo Lux': 'kilo Lux', u'degrees': 'degrees', '': ''}
_NAMES = {'sensor par': 'Sensor Photosynthetical Active Radiation'}


def formattingTheJSONFileAndReturnWavelengthAndSpectrum(fileLocation):
    '''
    This function will format the source JSON file including multiple JSON objects
    into a file of JSON array
    '''
    with open(fileLocation, 'r+') as fileHandler:
        tempList, linePosition, wavelengthList, spectrumList, j, k =\
            fileHandler.read().split('\n'), list(), [[]], [[]], 0, 0
        for i in range(len(tempList)):
            if "environment_sensor_set_reading" in tempList[i] and i > 2:
                linePosition.append(i - 1)
            if "wavelength" in tempList[i]:
                wavelengthList[j].append(
                    float(tempList[i][tempList[i].find(':') + 1: -2]))
            if "wavelength" not in tempList[i] and "wavelength" in tempList[i - 4]\
                    and "band" not in tempList[i] and "," not in tempList[i]:
                wavelengthList.append([])
                j += 1
                spectrumList.append([])
                k += 1
            if "spectrum" in tempList[i]:
                spectrumList[k].append(
                    float(tempList[i][tempList[i].find(':') + 1: -2]))
          
        wavelengthList.remove([])
        spectrumList.remove([])

        for line in linePosition:
            del tempList[line]
            tempList.insert(line, "},{")

        fileHandler.seek(0)
        fileHandler.truncate()

        if '[' not in tempList[0] and ']' not in tempList[-1]:
            tempList.insert(0, '[')
            tempList.append(']')
            fileHandler.write('\n'.join(tempList))
        else:
            fileHandler.write('\n'.join(tempList))
    return wavelengthList, spectrumList


def JSONHandler(fileLocation):
    '''
    Main JSON handler, write JSON file to a Python list with standard JSON module
    '''
    wavelength, spectrum = formattingTheJSONFileAndReturnWavelengthAndSpectrum(fileLocation)
    with open(fileLocation, 'r') as fileHandler:
        return json.loads(fileHandler.read()), wavelength, spectrum


def renameTheValue(name):
    '''
    Rename the value so it becomes legal in netCDF
    '''
    if type(name) is unicode:
        name = name.encode('ascii', 'ignore')
    if name in _UNIT_DICTIONARY:
        name = _UNIT_DICTIONARY[name]
    elif name in _NAMES:
        name = _NAMES[name]

    returningString = str()
    for letters in name:
        if letters == ' ':
            returningString += '_'
        else:
            returningString += letters
    return returningString


def getSpectrometerInformation(arrayOfJSON):
    '''
    Collect information from spectrometer with special care
    '''
    maxFixedIntensity = [int(intensityMembers["spectrometer"]["maxFixedIntensity"]) for intensityMembers in
                         arrayOfJSON]
    integrationTime = [int(integrateMembers["spectrometer"]["integration time in ?s"]) for integrateMembers in
                       arrayOfJSON]

    # TODO Reformat the spectrometer[bands]; illegal JSON formatting

    return maxFixedIntensity, integrationTime


def getListOfValue(arrayOfJSON, dataName):
    '''
    Collect data from JSON objects which have "value" member
    '''
    return [float(valueMembers[dataName]['value'].encode('ascii', 'ignore')) for valueMembers in arrayOfJSON]


def getListOfRawValue(arrayOfJSON, dataName):
    '''
    Collect data from JSON objects which have "rawValue" member
    '''
    return [float(valueMembers[dataName]['rawValue'].encode('ascii', 'ignore')) for valueMembers in arrayOfJSON]


def main(JSONArray, outputFileName, wavelength=None, spectrum=None):
    '''
    Main netCDF handler, write data to the netCDF file indicated.
    '''
    netCDFHandler = Dataset(outputFileName, 'w', format='NETCDF4')
    dataMemberList = [JSONMembers[u"environment_sensor_set_reading"]
                      for JSONMembers in JSONArray]
    timeStampList = [JSONMembers[u'timestamp']
                     for JSONMembers in dataMemberList]
    timeDimension = netCDFHandler.createDimension("time", None)
    tempTimeVariable = netCDFHandler.createVariable('time', str, ('time',))
    for i in range(len(timeStampList)):  # Assign Times
        tempTimeVariable[i] = timeStampList[i]

    for data in dataMemberList[0]:
        if data != 'spectrometer' and type(dataMemberList[0][data]) not in (str, unicode):
            tempVariable = netCDFHandler.createVariable(
                renameTheValue(data), 'f4', ('time',))
            tempVariable[:] = getListOfValue(
                dataMemberList, data)  # Assign "values"
            if 'unit' in dataMemberList[0][data]:  # Assign Units
                setattr(tempVariable, 'units', _UNIT_DICTIONARY[
                        dataMemberList[0][data]['unit']])
            if 'rawValue' in dataMemberList[0][data]:  # Assign "rawValues"
                netCDFHandler.createVariable(renameTheValue(data) + '_rawValue', 'f4', ('time',))[:] =\
                    getListOfRawValue(dataMemberList, data)
        elif type(dataMemberList[0][data]) in (str, unicode):
            netCDFHandler.createVariable(renameTheValue(data), str)[
                0] = dataMemberList[0][data]

        if data == 'spectrometer':  # Special care for spectrometers :)
            netCDFHandler.createVariable('Spectrometer_maxFixedIntensity', 'f4', ('time',))[:] =\
                getSpectrometerInformation(dataMemberList)[0]
            netCDFHandler.createVariable('Spectrometer_Integration_Time_In_Microseconds', 'f4', ('time',))[:] =\
                getSpectrometerInformation(dataMemberList)[1]

    if wavelength and spectrum:
    	netCDFHandler.createDimension("wavelength", len(wavelength[0]))
    	netCDFHandler.createVariable("wavelength", 'f4', ('time', 'wavelength'))[:,:] = wavelength
    	netCDFHandler.createVariable("spectrum", 'f4', ('time', 'wavelength'))[:,:] = spectrum
    	print "assigned"

    netCDFHandler.close()


if __name__ == '__main__':
    fileInputLocation, fileOutputLocation = sys.argv[1], sys.argv[2]
    if not os.path.exists(fileOutputLocation):
        os.mkdir(fileOutputLocation)  # Create folder

    if not os.path.isdir(fileInputLocation):
        tempJSONMasterList = JSONHandler(fileInputLocation)
        main(tempJSONMasterList, fileInputLocation.strip('.json') + '.nc')
    else:  # Read and Export netCDF to folder
        for filePath, fileDirectory, fileName in os.walk(fileInputLocation):
            for members in fileName:
                if os.path.join(filePath, members).endswith('.json'):
                    outputFileName = members.strip('.json') + '.nc'
                    tempJSONMasterList, wavelength, spectrum = JSONHandler(
                        os.path.join(filePath, members))
                    main(tempJSONMasterList, os.path.join(
                        fileOutputLocation, outputFileName), wavelength, spectrum)
