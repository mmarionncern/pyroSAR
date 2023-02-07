###############################################################################
# Convenience functions for SAR image batch processing with ESA SNAP

# Copyright (c) 2016-2022, the pyroSAR Developers.

# This file is part of the pyroSAR Project. It is subject to the
# license terms in the LICENSE.txt file found in the top-level
# directory of this distribution and at
# https://github.com/johntruckenbrodt/pyroSAR/blob/master/LICENSE.txt.
# No part of the pyroSAR project, including this file, may be
# copied, modified, propagated, or distributed except according
# to the terms contained in the LICENSE.txt file.
###############################################################################
import os
import re
import datetime as dt
from subprocess import Popen
import shutil
from ..drivers import identify, identify_many, ID
from .auxil import parse_recipe, parse_node, gpt, groupbyWorkers, writer, \
    windows_fileprefix, orb_parametrize, geo_parametrize, sub_parametrize, \
    mli_parametrize, dem_parametrize

from spatialist.ancillary import dissolve

import logging

log = logging.getLogger(__name__)



def geocode(infile, outdir, t_srs=4326, spacing=20, polarizations='all', shapefile=None, scaling='dB',
            geocoding_type='Range-Doppler', removeS1BorderNoise=True, removeS1BorderNoiseMethod='pyroSAR',
            removeS1ThermalNoise=True, offset=None, allow_RES_OSV=False, demName='SRTM 1Sec HGT',
            externalDEMFile=None, externalDEMNoDataValue=None, externalDEMApplyEGM=True, terrainFlattening=True,
            basename_extensions=None, test=False, export_extra=None, groupsize=1, cleanup=True, tmpdir=None,
            gpt_exceptions=None, gpt_args=None, returnWF=False, nodataValueAtSea=True,
            demResamplingMethod='BILINEAR_INTERPOLATION', imgResamplingMethod='BILINEAR_INTERPOLATION',
            alignToStandardGrid=False, standardGridOriginX=0, standardGridOriginY=0,
            speckleFilter=False, refarea='gamma0', clean_edges=False, clean_edges_npixels=1,
            rlks=None, azlks=None, decomposition_modes=None, dynamic_cleaning=False):
    """
    general function for geocoding of SAR backscatter images with SNAP.
    
    This function performs the following steps:
    
    - (if necessary) identify the SAR scene(s) passed via argument `infile` (:func:`pyroSAR.drivers.identify`)
    - (if necessary) create the directories defined via `outdir` and `tmpdir`
    - (if necessary) download Sentinel-1 OSV files
    - parse a SNAP workflow (:class:`pyroSAR.snap.auxil.Workflow`)
    - write the workflow to an XML file in `outdir`
    - execute the workflow (:func:`pyroSAR.snap.auxil.gpt`)

    Note
    ----
    The function may create workflows with multiple `Write` nodes. All nodes are parametrized to write data in ENVI format,
    in which case the node parameter `file` is going to be a directory. All nodes will use the same temporary directory,
    which will be created in `tmpdir`.
    Its name is created from the basename of the `infile` (:meth:`pyroSAR.drivers.ID.outname_base`)
    and a suffix identifying each processing node of the workflow (:meth:`pyroSAR.snap.auxil.Workflow.suffix`).
    
    For example: `S1A__IW___A_20180101T170648_NR_Orb_Cal_ML_TF_TC`.
    
    Parameters
    ----------
    infile: str or ~pyroSAR.drivers.ID or list
        The SAR scene(s) to be processed; multiple scenes are treated as consecutive acquisitions, which will be
        mosaicked with SNAP's SliceAssembly operator.
    outdir: str
        The directory to write the final files to.
    t_srs: int or str or osgeo.osr.SpatialReference
        A target geographic reference system in WKT, EPSG, PROJ4 or OPENGIS format.
        See function :func:`spatialist.auxil.crsConvert()` for details.
        Default: `4326 <https://spatialreference.org/ref/epsg/4326/>`_.
    spacing: int or float, optional
        The target pixel spacing in meters. Default is 20
    polarizations: list[str] or str
        The polarizations to be processed; can be a string for a single polarization, e.g. 'VV', or a list of several
        polarizations, e.g. ['VV', 'VH']. With the special value 'all' (default) all available polarizations are
        processed.
    shapefile: str or :py:class:`~spatialist.vector.Vector` or dict, optional
        A vector geometry for subsetting the SAR scene to a test site. Default is None.
    scaling: {'dB', 'db', 'linear'}, optional
        Should the output be in linear or decibel scaling? Default is 'dB'.
    geocoding_type: {'Range-Doppler', 'SAR simulation cross correlation'}, optional
        The type of geocoding applied; can be either 'Range-Doppler' (default) or 'SAR simulation cross correlation'
    removeS1BorderNoise: bool, optional
        Enables removal of S1 GRD border noise (default). Will be ignored if SLC scenes are processed.
    removeS1BorderNoiseMethod: str, optional
        The border noise removal method to be applied if `removeS1BorderNoise` is True.
        See :func:`pyroSAR.S1.removeGRDBorderNoise` for details. One of the following:
        
         - 'ESA': the pure implementation as described by ESA
         - 'pyroSAR': the ESA method plus the custom pyroSAR refinement (default)
    removeS1ThermalNoise: bool, optional
        Enables removal of S1 thermal noise (default).
    offset: tuple, optional
        A tuple defining offsets for left, right, top and bottom in pixels, e.g. (100, 100, 0, 0); this variable is
        overridden if a shapefile is defined. Default is None.
    allow_RES_OSV: bool
        (only applies to Sentinel-1) Also allow the less accurate RES orbit files to be used?
        The function first tries to download a POE file for the scene.
        If this fails and RES files are allowed, it will download the RES file.
        The selected OSV type is written to the workflow XML file.
        Processing is aborted if the correction fails (Apply-Orbit-File parameter continueOnFail set to false).
    demName: str
        The name of the auto-download DEM. Default is 'SRTM 1Sec HGT'. Is ignored when `externalDEMFile` is not None.
        Supported options:
        
         - ACE2_5Min
         - ACE30
         - ASTER 1sec GDEM
         - CDEM
         - Copernicus 30m Global DEM
         - Copernicus 90m Global DEM
         - GETASSE30
         - SRTM 1Sec Grid
         - SRTM 1Sec HGT
         - SRTM 3Sec
    externalDEMFile: str or None, optional
        The absolute path to an external DEM file. Default is None. Overrides `demName`.
    externalDEMNoDataValue: int, float or None, optional
        The no data value of the external DEM. If not specified (default) the function will try to read it from the
        specified external DEM.
    externalDEMApplyEGM: bool, optional
        Apply Earth Gravitational Model to external DEM? Default is True.
    terrainFlattening: bool
        Apply topographic normalization on the data?
    basename_extensions: list of str or None
        Names of additional parameters to append to the basename, e.g. ['orbitNumber_rel'].
    test: bool, optional
        If set to True the workflow xml file is only written and not executed. Default is False.
    export_extra: list or None
        A list of image file IDs to be exported to outdir. The following IDs are currently supported:
        
         - incidenceAngleFromEllipsoid
         - localIncidenceAngle
         - projectedLocalIncidenceAngle
         - DEM
         - layoverShadowMask
         - scatteringArea (requires ``terrainFlattening=True``)
         - gammaSigmaRatio (requires ``terrainFlattening=True`` and ``refarea=['sigma0', 'gamma0']``)
    groupsize: int
        The number of workers executed together in one gpt call.
    cleanup: bool
        Should all files written to the temporary directory during function execution be deleted after processing?
        Default is True.
    tmpdir: str or None
        Path of custom temporary directory, useful to separate output folder and temp folder. If `None`, the `outdir`
        location will be used. The created subdirectory will be deleted after processing if ``cleanup=True``.
    gpt_exceptions: dict or None
        A dictionary to override the configured GPT executable for certain operators;
        each (sub-)workflow containing this operator will be executed with the define executable;
        
         - e.g. ``{'Terrain-Flattening': '/home/user/snap/bin/gpt'}``
    gpt_args: list or None
        A list of additional arguments to be passed to the gpt call.
        
        - e.g. ``['-x', '-c', '2048M']`` for increased tile cache size and intermediate clearing
    returnWF: bool
        Return the full name of the written workflow XML file?
    nodataValueAtSea: bool
        Mask pixels acquired over sea? The sea mask depends on the selected DEM.
    demResamplingMethod: str
        One of the following:
        
         - 'NEAREST_NEIGHBOUR'
         - 'BILINEAR_INTERPOLATION'
         - 'CUBIC_CONVOLUTION'
         - 'BISINC_5_POINT_INTERPOLATION'
         - 'BISINC_11_POINT_INTERPOLATION'
         - 'BISINC_21_POINT_INTERPOLATION'
         - 'BICUBIC_INTERPOLATION'
    imgResamplingMethod: str
        The resampling method for geocoding the SAR image; the options are identical to demResamplingMethod.
    alignToStandardGrid: bool
        Align all processed images to a common grid?
    standardGridOriginX: int or float
        The x origin value for grid alignment
    standardGridOriginY: int or float
        The y origin value for grid alignment
    speckleFilter: str
        One of the following:
        
         - 'Boxcar'
         - 'Median'
         - 'Frost'
         - 'Gamma Map'
         - 'Refined Lee'
         - 'Lee'
         - 'Lee Sigma'
    refarea: str or list
        'sigma0', 'gamma0' or a list of both
    clean_edges: bool
        erode noisy image edges? See :func:`pyroSAR.snap.auxil.erode_edges`.
        Does not apply to layover-shadow mask.
    clean_edges_npixels: int
        the number of pixels to erode.
    rlks: int or None
        the number of range looks. If not None, overrides the computation done by function
        :func:`pyroSAR.ancillary.multilook_factors` based on the image pixel spacing and the target spacing.
    azlks: int or None
        the number of azimuth looks. Like `rlks`.
    decomposition_modes: str, list or None
        'H-alpha', 'c2' or list of both

    Returns
    -------
    str or None
        Either the name of the workflow file if ``returnWF == True`` or None otherwise
    
    
    .. figure:: figures/snap_geocode.svg
        :align: center
        
        Function geocode workflow diagram for processing Sentinel-1 scenes.
        Dashed lines depict optional steps. The output is sigma or gamma nought
        backscatter with ellipsoid or radiometric terrain correction (suffix elp/rtc)
        as well as several optional ancillary datasets (controlled via argument `export_extra`).

    Examples
    --------
    geocode a Sentinel-1 scene and export the local incidence angle map with it

    >>> from pyroSAR.snap import geocode
    >>> filename = 'S1A_IW_GRDH_1SDV_20180829T170656_20180829T170721_023464_028DE0_F7BD.zip'
    >>> geocode(infile=filename, outdir='outdir', tr=20, scaling='dB',
    >>>         export_extra=['DEM', 'localIncidenceAngle'], t_srs=4326)

    See Also
    --------
    :class:`pyroSAR.drivers.ID`,
    :class:`spatialist.vector.Vector`,
    :func:`spatialist.auxil.crsConvert()`
    """
    if clean_edges:
        try:
            import scipy
        except ImportError:
            raise RuntimeError('please install scipy to clean edges')
    
    if isinstance(infile, ID):
        id = infile
        ids = [id]
    elif isinstance(infile, str):
        id = identify(infile)
        ids = [id]
    elif isinstance(infile, list):
        ids = identify_many(infile, sortkey='start')
        id = ids[0]
    else:
        raise TypeError("'infile' must be of type str, list or pyroSAR.ID")

    if id.is_processed(outdir):
        log.info('scene {} already processed'.format(id.outname_base()))
        return
    
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    ############################################
    # general setup
    process_S1_SLC = False
    
    if id.sensor in ['ASAR', 'ERS1', 'ERS2']:
        formatName = 'ENVISAT'
    elif id.sensor in ['S1A', 'S1B']:
        if id.product == 'SLC':
            removeS1BorderNoise = False
            process_S1_SLC = True
        formatName = 'SENTINEL-1'
    else:
        raise RuntimeError('sensor not supported (yet)')

    if decomposition_modes is not None:
        if not process_S1_SLC:
            raise RuntimeError("Signal decomposition is only valid for SLC data")
        
        valid_decompositions = ["C2","H-alpha"]
        if set(valid_decompositions+decomposition_modes)!=set(valid_decompositions):
            raise RuntimeError(f"Decomposition modes {decomposition_modes} are not valid, authorized modes are {valid_decompositions}")
        
        #thermal nois removal not suited for signal decomposition
        removeS1ThermalNoise = False
        #unrelevant terrain flattening
        terrainFlattening = False
        #sclaing has to be linear
        scaling = 'linear'
        if isinstance(decomposition_modes,str):
            decompositions = [decomposition_modes]
        elif isinstance(decomposition_modes, list):
            decompositions = decomposition_modes
        else:
            raise RuntimeError("decomposition_modes must be of type str or list")
    else:
        decompositions=[]
        
    # several options like resampling are modified globally for the whole workflow at the end of this function
    # this list gathers IDs of nodes for which this should not be done because they are configured individually
    resampling_exceptions = []
    ######################
    if isinstance(polarizations, str):
        if polarizations == 'all':
            polarizations = id.polarizations
        else:
            if polarizations in id.polarizations:
                polarizations = [polarizations]
            else:
                raise RuntimeError('polarization {} does not exists in the source product'.format(polarizations))
    elif isinstance(polarizations, list):
        polarizations = [x for x in polarizations if x in id.polarizations]
    else:
        raise RuntimeError('polarizations must be of type str or list')
    
    swaths = None
    if process_S1_SLC:
        if id.acquisition_mode == 'IW':
            swaths = ['IW1', 'IW2', 'IW3']
        elif id.acquisition_mode == 'EW':
            swaths = ['EW1', 'EW2', 'EW3', 'EW4', 'EW5']
        elif re.search('S[1-6]', id.acquisition_mode):
            pass
        else:
            raise RuntimeError('acquisition mode {} not supported'.format(id.acquisition_mode))
    
    bandnames = dict()
    bandnames['beta0'] = ['Beta0_' + x for x in polarizations]
    bandnames['gamma0'] = ['Gamma0_' + x for x in polarizations]
    bandnames['sigma0'] = ['Sigma0_' + x for x in polarizations]
    bandnames['int'] = ['Intensity_' + x for x in polarizations]
    ############################################
    ############################################
    # parse base workflow
    workflow = parse_recipe('blank')
    ############################################
    if not isinstance(infile, list):
        infile = [infile]
    
    last = None
    collect = []
    for i in range(0, len(infile)):
        ############################################
        # Read node configuration
        read = parse_node('Read')
        workflow.insert_node(read)
        read.parameters['file'] = ids[i].scene
        read.parameters['formatName'] = formatName
        last = read
        ############################################
        # Remove-GRD-Border-Noise node configuration
        if id.sensor in ['S1A', 'S1B'] and id.product == 'GRD' and removeS1BorderNoise:
            bn = parse_node('Remove-GRD-Border-Noise')
            workflow.insert_node(bn, before=last.id)
            bn.parameters['selectedPolarisations'] = polarizations
            last = bn
        ############################################
        # Calibration node configuration
        cal = parse_node('Calibration')
        workflow.insert_node(cal, before=last.id)
        cal.parameters['selectedPolarisations'] = polarizations
        if isinstance(refarea, str):
            refarea = [refarea]
        for item in refarea:
            if item not in ['sigma0', 'gamma0']:
                raise ValueError('unsupported value for refarea: {}'.format(item))
        if terrainFlattening:
            cal.parameters['outputBetaBand'] = True
            cal.parameters['outputSigmaBand'] = False
        else:
            for opt in refarea:
                cal.parameters['output{}Band'.format(opt[:-1].capitalize())] = True
        if id.sensor in ['ERS1', 'ERS2', 'ASAR']:
            cal.parameters['createBetaBand'] = True
        if len(decompositions)!=0:
            cal.parameters['outputImageInComplex'] = True
        last = cal
        ############################################
        # ThermalNoiseRemoval node configuration
        if id.sensor in ['S1A', 'S1B'] and removeS1ThermalNoise:
            tn = parse_node('ThermalNoiseRemoval')
            workflow.insert_node(tn, before=last.id)
            tn.parameters['selectedPolarisations'] = polarizations
            last = tn
        collect.append(last.id)
    ############################################
    # SliceAssembly node configuration
    if len(collect) > 1:
        sliceAssembly = parse_node('SliceAssembly')
        sliceAssembly.parameters['selectedPolarisations'] = polarizations
        workflow.insert_node(sliceAssembly, before=collect)
        last = sliceAssembly
    ############################################
    # TOPSAR-Deburst node configuration
    if process_S1_SLC and swaths is not None:
        deb = parse_node('TOPSAR-Deburst')
        workflow.insert_node(deb, before=last.id)
        deb.parameters['selectedPolarisations'] = polarizations
        last = deb
    ############################################
    # Apply-Orbit-File node configuration
    orb = orb_parametrize(scene=id, formatName=formatName, allow_RES_OSV=allow_RES_OSV)
    workflow.insert_node(orb, before=last.id)
    last = orb
    ############################################
    # Subset node configuration
    if shapefile is not None or offset is not None:
        sub = sub_parametrize(scene=id, geometry=shapefile, offset=offset, buffer=0.01)
        workflow.insert_node(sub, before=last.id)
        last = sub
    ############################################
    last_ids=[]
    # Matrix decomposition node configurations
    if "H-alpha" in decompositions or "C2" in decompositions:
        ##create C2 covariance matrix
        pol_m = parse_node("Polarimetric-Matrices")
        pol_m.parameters["matrix"] = "C2"
        workflow.insert_node(pol_m, before=last)
        last = pol_m
        bands = ["C11", "C12_real", "C12_imag", "C22"]
        
    ############################################
    # Multilook node configuration
    if id.sensor in ['ERS1', 'ERS2', 'ASAR']:
        bands = bandnames['beta0'] + bandnames['sigma0']
    else:
        bands = None
    ml = mli_parametrize(scene=id, spacing=spacing, rlks=rlks, azlks=azlks,
                         sourceBands=bands)
    if ml is not None:
        workflow.insert_node(ml, before=last.id)
        last = ml
        if not speckleFilter:
            last_ids.append(last.id)
    else:
        if not speckleFilter:
            last_ids.append(pol_m.id)
        
    ############################################
    # Terrain-Flattening node configuration
    tf = None
    if terrainFlattening:
        tf = parse_node('Terrain-Flattening')
        workflow.insert_node(tf, before=last.id)
        tf.parameters['sourceBands'] = bandnames['beta0']
        if 'reGridMethod' in tf.parameters.keys():
            if externalDEMFile is None:
                tf.parameters['reGridMethod'] = True
            else:
                tf.parameters['reGridMethod'] = False
        if 'sigma0' in refarea:
            try:
                tf.parameters['outputSigma0'] = True
            except KeyError:
                raise RuntimeError("The Terrain-Flattening node does not accept "
                                   "parameter 'outputSigma0'. Please update S1TBX.")
        last = tf
    ############################################
    # Speckle-Filter node configuration
    speckleFilter_options = ['Boxcar',
                             'Median',
                             'Frost',
                             'Gamma Map',
                             'Refined Lee',
                             'Lee',
                             'Lee Sigma']
    polarimetricSpeckleFilter_options = ['Box Car Filter',
                                         'IDAN Filter',
                                         'Refined Lee Filter',
                                         'Improved Lee Sigma Filter']
    
    if speckleFilter:
        message = '{0} must be one of the following:\n- {1}'
        if speckleFilter not in speckleFilter_options and len(decompositions)==0:
            raise ValueError(message.format('speckleFilter', '\n- '.join(speckleFilter_options)))
        if speckleFilter not in polarimetricSpeckleFilter_options and len(decompositions)!=0:
            raise ValueError(message.format('speckleFilter', '\n- '.join(speckleFilter_options)))

        if len(decompositions)==0:
            sf = parse_node('Speckle-Filter')
            workflow.insert_node(sf, before=last.id)
            sf.parameters['sourceBands'] = None
            sf.parameters['filter'] = speckleFilter
            last = sf
            last_ids = [last.id]
        else:
            ##polaricmetric speckle filtering
            sf = parse_node("Polarimetric-Speckle-Filter")
            workflow.insert_node(sf, before=last.id)
            sf.parameters["filter"] = speckleFilter
            last = sf
            last_ids.append(last.id)

    ############################################
    # Dual polarization H-alpha decomposition
    if "H-alpha" in decompositions:
        pol_dc = parse_node("Polarimetric-Decomposition")
        workflow.insert_node(pol_dc, before=last.id)
        pol_dc.parameters["decomposition"] = "H-Alpha Dual Pol Decomposition"
        pol_dc.parameters["windowSize"] = 5
        pol_dc.parameters["outputHAAlpha"] = True
        last = pol_dc
    
    ############################################
    # merge bands to pass them to Terrain-Correction
    bm_tc = None
    if len(decompositions)==0:
        bands = dissolve([bandnames[opt] for opt in refarea])
    else:
        bands=[]
        nodes=[]
        if "H-alpha" in decompositions:
            bands.extend(["Alpha", "Entropy", "Anisotropy"])
            nodes.append(last.id)
        if "C2" in decompositions:
            bands.extend(["C11", "C12_real", "C12_imag", "C22"])
            nodes.append(sf.id)
        bm_tc = parse_node('BandMerge')
        workflow.insert_node(bm_tc, before=nodes)
        bm_tc.parameters['sourceBands']=bands
        last = bm_tc
        
    if len(refarea) > 1 and terrainFlattening and 'scatteringArea' in export_extra:
        bm_tc = parse_node('BandMerge')
        workflow.insert_node(bm_tc, before=[last.source, last.id])
        sources = bm_tc.source
        gamma_index = sources.index('Terrain-Flattening')
        sigma_index = abs(gamma_index - 1)
        s1_id = os.path.basename(os.path.splitext(id.scene)[0])
        bands_long = []
        for band in bands:
            comp = [band + '::']
            if shapefile is not None:
                comp.append('Subset_')
            comp.append(s1_id)
            if band.startswith('Gamma'):
                comp.append('_' + workflow.suffix(stop=sources[gamma_index]))
            else:
                comp.append('_' + workflow.suffix(stop=sources[sigma_index]))
            bands_long.append(''.join(comp))
        bm_tc.parameters['sourceBands'] = bands_long
        last = bm_tc
        last_ids.append(last.id)
    ############################################
    # configuration of node sequence for specific geocoding approaches
    tc = geo_parametrize(spacing=spacing, t_srs=t_srs,
                         tc_method=geocoding_type, sourceBands=bands,
                         alignToStandardGrid=alignToStandardGrid,
                         standardGridOriginX=standardGridOriginX,
                         standardGridOriginY=standardGridOriginY)
    workflow.insert_node(tc, before=last)
    if isinstance(tc, list):
        last = tc = tc[-1]
    else:
        last = tc
    ############################################
    # (optionally) add node for conversion from linear to db scaling
    if scaling not in ['dB', 'db', 'linear']:
        raise RuntimeError('scaling must be  a string of either "dB", "db" or "linear"')
    
    if scaling in ['dB', 'db']:
        lin2db = parse_node('LinearToFromdB')
        workflow.insert_node(lin2db, before=last.id)
        lin2db.parameters['sourceBands'] = bands
        last = lin2db
    ############################################
    # parametrize write node
    # create a suffix for the output file to identify processing steps performed in the workflow
    suffix = workflow.suffix()
    if tmpdir is None:
        tmpdir = outdir
    basename = os.path.join(tmpdir, id.outname_base(basename_extensions))
    outname = basename + '_' + suffix
    
    write = parse_node('Write')
    workflow.insert_node(write, before=last.id)
    write.parameters['file'] = outname
    write.parameters['formatName'] = 'ENVI'
    ############################################
    ############################################
    if export_extra is not None:
        tc_options = ['incidenceAngleFromEllipsoid',
                      'localIncidenceAngle',
                      'projectedLocalIncidenceAngle',
                      'DEM',
                      'layoverShadowMask']
        tc_selection = []
        for item in export_extra:
            if item in tc_options:
                key = 'save{}{}'.format(item[0].upper(), item[1:])
                tc.parameters[key] = True
                if item == 'DEM':
                    tc_selection.append('elevation')
                else:
                    tc_selection.append(item)
            elif item == 'scatteringArea':
                if not terrainFlattening:
                    raise RuntimeError('scatteringArea can only be created if terrain flattening is performed')
                area_select = parse_node('BandSelect')
                workflow.insert_node(area_select, before=tf.source, resetSuccessorSource=False)
                area_select.parameters['sourceBands'] = bandnames['beta0']
                
                area_merge1 = parse_node('BandMerge')
                workflow.insert_node(area_merge1, before=[tf.id, area_select.id], resetSuccessorSource=False)
                
                math = parse_node('BandMaths')
                workflow.insert_node(math, before=area_merge1.id, resetSuccessorSource=False)
                
                pol = polarizations[0]  # the result will be the same for each polarization
                area = 'scatteringArea_{0}'.format(pol)
                expression = 'Beta0_{0} / Gamma0_{0}'.format(pol)
                
                math.parameters.clear_variables()
                exp = math.parameters['targetBands'][0]
                exp['name'] = area
                exp['type'] = 'float32'
                exp['expression'] = expression
                exp['noDataValue'] = 0.0
                
                if len(refarea) > 1:
                    bm_tc.source = bm_tc.source + [math.id]
                else:
                    bm_tc = parse_node('BandMerge')
                    workflow.insert_node(bm_tc, before=[tf.id, math.id], resetSuccessorSource=False)
                    tc.source = bm_tc.id
                
                # modify Terrain-Correction source bands
                tc_bands = tc.parameters['sourceBands'] + ',' + area
                tc.parameters['sourceBands'] = tc_bands
                
                # add scattering Area to list of band directly written from Terrain-Correction
                tc_selection.append(area)
            elif item == 'gammaSigmaRatio':
                if not terrainFlattening:
                    raise RuntimeError('gammaSigmaRatio can only be created if terrain flattening is performed')
                if sorted(refarea) != ['gamma0', 'sigma0']:
                    raise ValueError("For export_extra layer 'gammaSigmaRatio' 'refarea' "
                                     "must contain both sigma0 and gamma0")
                math = parse_node('BandMaths')
                workflow.insert_node(math, before=tf.id, resetSuccessorSource=False)
                
                pol = polarizations[0]  # the result will be the same for each polarization
                ratio = 'gammaSigmaRatio_{0}'.format(pol)
                expression = 'Sigma0_{0} / Gamma0_{0}'.format(pol)
                
                math.parameters.clear_variables()
                exp = math.parameters['targetBands'][0]
                exp['name'] = ratio
                exp['type'] = 'float32'
                exp['expression'] = expression
                exp['noDataValue'] = 0.0
                
                if len(refarea) > 1:
                    bm_tc.source = bm_tc.source + [math.id]
                else:
                    bm_tc = parse_node('BandMerge')
                    workflow.insert_node(bm_tc, before=[tf.id, math.id], resetSuccessorSource=False)
                    tc.source = bm_tc.id
                
                # modify Terrain-Correction source bands
                tc_bands = tc.parameters['sourceBands'] + ',' + ratio
                tc.parameters['sourceBands'] = tc_bands
                
                # add scattering Area to list of band directly written from Terrain-Correction
                tc_selection.append(ratio)
            else:
                raise RuntimeError("ID '{}' not valid for argument 'export_extra'".format(item))
        # directly write export_extra layers to avoid dB scaling
        if scaling in ['db', 'dB'] and len(tc_selection) > 0:
            tc_write = parse_node('Write')
            workflow.insert_node(tc_write, before=tc.id, resetSuccessorSource=False)
            tc_write.parameters['file'] = outname
            tc_write.parameters['formatName'] = 'ENVI'
            tc_select = parse_node('BandSelect')
            workflow.insert_node(tc_select, after=tc_write.id)
            tc_select.parameters['sourceBands'] = tc_selection
    ############################################
    ############################################
    # DEM handling
    dem_parametrize(workflow=workflow, demName=demName,
                    externalDEMFile=externalDEMFile,
                    externalDEMNoDataValue=externalDEMNoDataValue,
                    externalDEMApplyEGM=externalDEMApplyEGM)
    ############################################
    ############################################
    # configure the resampling methods
    
    options_img = ['NEAREST_NEIGHBOUR',
                   'BILINEAR_INTERPOLATION',
                   'CUBIC_CONVOLUTION',
                   'BISINC_5_POINT_INTERPOLATION',
                   'BISINC_11_POINT_INTERPOLATION',
                   'BISINC_21_POINT_INTERPOLATION',
                   'BICUBIC_INTERPOLATION']
    options_dem = options_img + ['DELAUNAY_INTERPOLATION']
    
    message = '{0} must be one of the following:\n- {1}'
    if demResamplingMethod not in options_dem:
        raise ValueError(message.format('demResamplingMethod', '\n- '.join(options_dem)))
    if imgResamplingMethod not in options_img:
        raise ValueError(message.format('imgResamplingMethod', '\n- '.join(options_img)))
    
    workflow.set_par('demResamplingMethod', demResamplingMethod)
    workflow.set_par('imgResamplingMethod', imgResamplingMethod,
                     exceptions=resampling_exceptions)
    ############################################
    ############################################
    # additional parameter settings applied to the whole workflow
    
    workflow.set_par('nodataValueAtSea', nodataValueAtSea)
    ############################################
    ############################################
    # write workflow to file and optionally execute it
    log.debug('writing workflow to file')
    
    wf_name = outname.replace(tmpdir, outdir) + '_proc.xml'
    workflow.write(wf_name)
    
    # execute the newly written workflow
    if not test:
        try:
            p=None
            if groupsize>1 and dynamic_cleaning:
                tmp_name +="/sub"
                p = Popen(['python', 'dynamicCOHCleaning.py', tmp_name])

            groups = groupbyWorkers(wf_name, groupsize)
            gpt(wf_name, groups=groups, cleanup=cleanup, tmpdir=outname,
                gpt_exceptions=gpt_exceptions, gpt_args=gpt_args,
                removeS1BorderNoiseMethod=removeS1BorderNoiseMethod)
            writer(xmlfile=wf_name, outdir=outdir, basename_extensions=basename_extensions,
                   clean_edges=clean_edges, clean_edges_npixels=clean_edges_npixels)
            if p:
                p.terminate()
        except Exception as e:
            log.info(str(e))
            with open(wf_name.replace('_proc.xml', '_error.log'), 'w') as logfile:
                logfile.write(str(e))
        finally:
            if cleanup and os.path.isdir(outname):
                log.info('deleting temporary files')
                shutil.rmtree(outname, onerror=windows_fileprefix)
        log.info('done')
    if returnWF:
        return wf_name


def noise_power(infile, outdir, polarizations, spacing, t_srs, refarea='sigma0', tmpdir=None, test=False, cleanup=True,
                demName='SRTM 1Sec HGT', externalDEMFile=None, externalDEMNoDataValue=None, externalDEMApplyEGM=True,
                alignToStandardGrid=False, standardGridOriginX=0, standardGridOriginY=0, groupsize=1,
                clean_edges=False, clean_edges_npixels=1, rlks=None, azlks=None):
    """
    Generate Sentinel-1 noise power images for each polarization, calibrated to either beta, sigma or gamma nought.
    The written GeoTIFF files will carry the suffix NEBZ, NESZ or NEGZ respectively.

    Parameters
    ----------
    infile: str
        The SAR scene(s) to be processed
    outdir: str
        The directory to write the final files to.
    polarizations: list[str]
        The polarizations to be processed, e.g. ['VV', 'VH'].
    spacing: int or float
        The target pixel spacing in meters.
    t_srs: int or str or osgeo.osr.SpatialReference
        A target geographic reference system in WKT, EPSG, PROJ4 or OPENGIS format.
    refarea: str
        either 'beta0', 'gamma0' or 'sigma0'.
    tmpdir: str
        Path of custom temporary directory, useful to separate output folder and temp folder. If `None`, the `outdir`
        location will be used. The created subdirectory will be deleted after processing if ``cleanup=True``.
    test: bool
        If set to True the workflow xml file is only written and not executed. Default is False.
    cleanup: bool
        Should all files written to the temporary directory during function execution be deleted after processing?
        Default is True.
    demName: str
        The name of the auto-download DEM. Default is 'SRTM 1Sec HGT'. Is ignored when `externalDEMFile` is not None.
        Supported options:
        
         - ACE2_5Min
         - ACE30
         - ASTER 1sec GDEM
         - CDEM
         - Copernicus 30m Global DEM
         - Copernicus 90m Global DEM
         - GETASSE30
         - SRTM 1Sec Grid
         - SRTM 1Sec HGT
         - SRTM 3Sec
    externalDEMFile: str or None, optional
        The absolute path to an external DEM file. Default is None. Overrides `demName`.
    externalDEMNoDataValue: int, float or None, optional
        The no data value of the external DEM. If not specified (default) the function will try to read it from the
        specified external DEM.
    externalDEMApplyEGM: bool, optional
        Apply Earth Gravitational Model to external DEM? Default is True.
    alignToStandardGrid: bool
        Align all processed images to a common grid?
    standardGridOriginX: int or float
        The x origin value for grid alignment
    standardGridOriginY: int or float
        The y origin value for grid alignment
    groupsize: int
        The number of workers executed together in one gpt call.
    clean_edges: bool
        erode noisy image edges? See :func:`pyroSAR.snap.auxil.erode_edges`.
        Does not apply to layover-shadow mask.
    clean_edges_npixels: int
        the number of pixels to erode.
    rlks: int or None
        the number of range looks. If not None, overrides the computation done by function
        :func:`pyroSAR.ancillary.multilook_factors` based on the image pixel spacing and the target spacing.
    azlks: int or None
        the number of azimuth looks. Like `rlks`.
    
    Returns
    -------

    """
    if clean_edges:
        try:
            import scipy
        except ImportError:
            raise RuntimeError('please install scipy to clean edges')
    
    if refarea not in ['beta0', 'sigma0', 'gamma0']:
        raise ValueError('refarea not supported')
    
    id = identify(infile)
    
    if id.sensor not in ['S1A', 'S1B']:
        raise RuntimeError('this function is for Sentinel-1 only')
    
    os.makedirs(outdir, exist_ok=True)
    if tmpdir is not None:
        os.makedirs(tmpdir, exist_ok=True)
    
    wf = parse_recipe('blank')
    
    read = parse_node('Read')
    read.parameters['file'] = infile
    wf.insert_node(read)
    ############################################
    orb = orb_parametrize(scene=id, workflow=wf, before=read.id,
                          formatName='SENTINEL-1', allow_RES_OSV=True)
    ############################################
    cal = parse_node('Calibration')
    wf.insert_node(cal, before=orb.id)
    cal.parameters['selectedPolarisations'] = polarizations
    cal.parameters['outputBetaBand'] = False
    cal.parameters['outputSigmaBand'] = False
    cal.parameters['outputGammaBand'] = False
    
    inband = refarea.capitalize()
    cal.parameters['output{}Band'.format(inband[:-1])] = True
    
    tnr = parse_node('ThermalNoiseRemoval')
    wf.insert_node(tnr, before=cal.id)
    if 'outputNoise' in tnr.parameters.keys():
        tnr.parameters['outputNoise'] = True
    last = tnr
    ############################################
    if id.product == 'SLC' and id.acquisition_mode in ['EW', 'IW']:
        deb = parse_node('TOPSAR-Deburst')
        wf.insert_node(deb, before=tnr.id)
        last = deb
    ############################################
    select = parse_node('BandSelect')
    wf.insert_node(select, before=last.id)
    measure = 'NE{}Z'.format(refarea.capitalize()[0])
    bands = ['{}_{}'.format(measure, pol) for pol in polarizations]
    select.parameters['sourceBands'] = bands
    last = select
    ############################################
    # Multilook node configuration
    ml = mli_parametrize(scene=id, spacing=spacing, rlks=rlks, azlks=azlks)
    if ml is not None:
        wf.insert_node(ml, before=last.id)
        last = ml
    ############################################
    tc = geo_parametrize(spacing=spacing, t_srs=t_srs, demName=demName,
                         externalDEMFile=externalDEMFile,
                         externalDEMNoDataValue=externalDEMNoDataValue,
                         externalDEMApplyEGM=externalDEMApplyEGM,
                         alignToStandardGrid=alignToStandardGrid,
                         standardGridOriginX=standardGridOriginX,
                         standardGridOriginY=standardGridOriginY)
    wf.insert_node(tc, before=last.id)
    last = tc
    ############################################
    
    suffix = wf.suffix()
    if tmpdir is None:
        tmpdir = outdir
    basename = id.outname_base() + '_' + suffix
    procdir = os.path.join(tmpdir, basename)
    outname = os.path.join(procdir, basename + '.dim')
    
    write = parse_node('Write')
    wf.insert_node(write, before=last.id)
    write.parameters['file'] = outname
    write.parameters['formatName'] = 'BEAM-DIMAP'
    
    wf_name = os.path.join(outdir, basename + '_proc.xml')
    wf.write(wf_name)
    
    if not test:
        groups = groupbyWorkers(wf_name, groupsize)
        gpt(xmlfile=wf_name, tmpdir=procdir, groups=groups, cleanup=cleanup)
        writer(xmlfile=wf_name, outdir=outdir, clean_edges=clean_edges,
               clean_edges_npixels=clean_edges_npixels)
        if cleanup:
            if os.path.isdir(procdir):
                shutil.rmtree(procdir, onerror=windows_fileprefix)


def halpha(infile, swaths=["IW1","IW2","IW3"], t_srs=4326, demName='SRTM 1Sec HGT',
           allow_RES_OSV=False, demResamplingMethodBGC="BICUBIC_INTERPOLATION",
           decomposition_modes=["H-alpha","C2"], speckleFilter="IDAN Filter",
           maskOutAreaWithoutElevation=False, externalDEMFile=None,
           externalDEMNoDataValue=None, externalDEMApplyEGM=True,
           demResamplingMethod='BILINEAR_INTERPOLATION',
           imgResamplingMethod='BILINEAR_INTERPOLATION',
           alignToStandardGrid=False, standardGridOriginX=0, standardGridOriginY=0,
           gpt_exceptions=None, gpt_args=None, rlks=None, azlks=None, spacing=20,
           cleanup=True, tmpdir=None, outName=None, test=False, enhancedSpectralDiversity=True,
           clean_edges=False, clean_edges_npixels=1, outdir=None, basename_extensions=None,
           returnWF=False,groupsize=2, removeS1BorderNoiseMethod='pyroSAR',
           dynamic_cleaning=True,geocoding_type='Range-Doppler'):

    if not isinstance(infile, str):
        raise RuntimeError("'infile' must be of type str")
                
    id = identify(infile)

    if id.is_processed(outdir):
        log.info(f'scene {id.outname_base()} already processed')
        return

    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    
    ############################################
    # general setup
    process_S1_SLC = False
    if id.sensor not in ['S1A', 'S1B', 'S1C'] or id.product!="SLC" :
        raise RuntimeError('Insar coherence only available for Sentinel mission in SLC mode')

    formatName = 'SENTINEL-1'

    if id.acquisition_mode != 'IW':
        raise RuntimeError(f"acquisition mode {id_1.acquisition_mode}/{id_2.acquisition_mode} not supported")

    ######################

    if externalDEMFile is None and externalDEMNoDataValue is None:
        externalDEMNoDataValue=0.

    if tmpdir is None:
        tmpdir = outdir
    
    for iw in swaths:

        # parse base workflow
        workflow = parse_recipe('blank')
        ############################################

        ############################################
        # Read node configuration
        read = parse_node('Read')
        workflow.insert_node(read)
        read.parameters['file'] = id.scene
        read.parameters['formatName'] = formatName
        last = read

        ############################################
        #TOP-SAR split
        ts_split = parse_node('TOPSAR-Split')
        workflow.insert_node(ts_split, before=last.id)
        ts_split.parameters["subswath"] = iw
        ts_split.parameters["selectedPolarisations"] = id.polarizations
        last = ts_split

        ############################################
        # Calibration node configuration
        cal = parse_node('Calibration')
        workflow.insert_node(cal, before=last.id)
        cal.parameters['selectedPolarisations'] = id.polarizations
        cal.parameters['outputImageInComplex'] = True
        last = cal
        ############################################
        # TOPSAR-Deburst node configuration
        deb = parse_node('TOPSAR-Deburst')
        workflow.insert_node(deb, before=last.id)
        deb.parameters['selectedPolarisations'] = id.polarizations
        last = deb

        ############################################
        # Apply-Orbit-File node configuration in coherence estimation
        orb = orb_parametrize(scene=id, formatName=formatName, allow_RES_OSV=allow_RES_OSV)
        workflow.insert_node(orb, before=last.id)
        #read_iw_orb_ids.append(orb.id)
        last = orb

        ############################################
        # Write temporary file
    
        date= id.start.split("T")[0]
        tmp_out_name= f"S1_HAlpha_{iw}_{date}_tmp"
        tmp_out= os.path.join(tmpdir, tmp_out_name)
    
        write_c2=parse_node("Write")
        workflow.insert_node(write_c2, before= last.id)
        write_c2.parameters["file"]= tmp_out
        write_c2.parameters["formatName"]= "BEAM-DIMAP"

        tmp_wf_name = f"{tmpdir}/Coh_tmp_prep_graph.xml"
        workflow.write(tmp_wf_name)
        if not test:
            gpt(tmp_wf_name, groups=None, cleanup=False, tmpdir=tmpdir,
                gpt_exceptions=gpt_exceptions, gpt_args=gpt_args,
                removeS1BorderNoiseMethod=removeS1BorderNoiseMethod)
    
    import glob
    tmp_c2_files= glob.glob(f"{tmpdir}/S1_HAlpha_*_{date}_tmp.dim")
    print(tmp_c2_files)
    
    workflow = parse_recipe('blank')
    reads=[]
    for c2f in tmp_c2_files:
        ############################################
        # Read node configuration for 2nd part
        read = parse_node('Read')
        workflow.insert_node(read, resetSuccessorSource=False)
        read.parameters['file'] = c2f
        read.parameters['formatName'] = "BEAM-DIMAP"
        reads.append(read.id)
        
    ############################################
    # merge sub-swaths node configuration
    merge=parse_node("TOPSAR-Merge")
    merge.parameters["selectedPolarisations"]=id.polarizations
    workflow.insert_node(merge, before=reads)
    last = merge

    ############################################
    ##create C2 covariance matrix
    pol_m = parse_node("Polarimetric-Matrices")
    workflow.insert_node(pol_m, before=last.id)
    pol_m.parameters["matrix"] = "C2"
    last = pol_m
    bands = ["C11", "C12_real", "C12_imag", "C22"]

    tmp_out_name= f"S1_HAlpha_c2_{date}_tmp"
    tmp_out= os.path.join(tmpdir, tmp_out_name)

    write_c2=parse_node("Write")
    workflow.insert_node(write_c2, before= last.id)
    write_c2.parameters["file"]= tmp_out
    write_c2.parameters["formatName"]= "BEAM-DIMAP"

    tmp_wf_name = f"{tmpdir}/Coh_tmp_c2_prep_graph.xml"
    workflow.write(tmp_wf_name)
    if not test:
        gpt(tmp_wf_name, groups=None, cleanup=False, tmpdir=tmpdir,
            gpt_exceptions=gpt_exceptions, gpt_args=gpt_args,
            removeS1BorderNoiseMethod=removeS1BorderNoiseMethod)

    ############################################
    # last phase        
    workflow = parse_recipe('blank')

    read = parse_node('Read')
    workflow.insert_node(read, resetSuccessorSource=False)
    read.parameters['file'] = f"{tmpdir}/S1_HAlpha_c2_{date}_tmp.dim"
    read.parameters['formatName'] = "BEAM-DIMAP"
    last = read
    
    
    ############################################
    # Multilook node configuration
    bands = None
    ml = mli_parametrize(scene=id, spacing=spacing, rlks=rlks, azlks=azlks,
                         sourceBands=bands, outputIntensity=False)
    if ml is not None:
        workflow.insert_node(ml, before=last.id)
        last = ml




        
    ############################################
    # Speckle-Filter node configuration
    polarimetricSpeckleFilter_options = ['Box Car Filter',
                                         'IDAN Filter',
                                         'Refined Lee Filter',
                                         'Improved Lee Sigma Filter']

    if speckleFilter:
        message = '{0} must be one of the following:\n- {1}'
    if speckleFilter not in polarimetricSpeckleFilter_options:
        raise ValueError(message.format('speckleFilter', '\n- '.join(speckleFilter_options)))

    ##polaricmetric speckle filtering
    sf = parse_node("Polarimetric-Speckle-Filter")
    workflow.insert_node(sf, before=last.id)
    sf.parameters["filter"] = speckleFilter
    last = sf
    #last_ids.append(last.id)
    
    ############################################
    # Dual polarization H-alpha decomposition
    if "H-alpha" in decomposition_modes:
        pol_dc = parse_node("Polarimetric-Decomposition")
        workflow.insert_node(pol_dc, before=last.id)
        pol_dc.parameters["decomposition"] = "H-Alpha Dual Pol Decomposition"
        pol_dc.parameters["windowSize"] = 5
        pol_dc.parameters["outputHAAlpha"] = True
        last = pol_dc

    ############################################
    # merge bands to pass them to Terrain-Correction
    if len(decomposition_modes)>1:
        bm_tc = None
        bands=[]
        nodes=[]
        if "H-alpha" in decomposition_modes:
            bands.extend(["Alpha", "Entropy", "Anisotropy"])
            nodes.append(last.id)
        if "C2" in decomposition_modes:
            bands.extend(["C11", "C12_real", "C12_imag", "C22"])
            nodes.append(sf.id)
        bm_tc = parse_node('BandMerge')
        workflow.insert_node(bm_tc, before=nodes)
        bm_tc.parameters['sourceBands']=bands
        last = bm_tc

    ############################################
    # configuration of node sequence for specific geocoding approaches
    tc = geo_parametrize(spacing=spacing, t_srs=t_srs,
                         tc_method=geocoding_type, sourceBands=bands,
                         alignToStandardGrid=alignToStandardGrid,
                         standardGridOriginX=standardGridOriginX,
                         standardGridOriginY=standardGridOriginY)
    workflow.insert_node(tc, before=last.id)
    last = tc

    ############################################
    # parametrize write node
    # create a suffix for the output file to identify processing steps performed in the workflow
    suffix = workflow.suffix()
    if tmpdir is None:
        tmpdir = outdir
    basename = os.path.join(tmpdir, id.outname_base(basename_extensions))
    outname = basename + '_' + suffix

    write = parse_node('Write')
    workflow.insert_node(write, before=last.id)
    write.parameters['file'] = outname
    write.parameters['formatName'] = 'ENVI'

    ############################################
    # DEM handling
    dem_parametrize(workflow=workflow, demName=demName,
                    externalDEMFile=externalDEMFile,
                    externalDEMNoDataValue=externalDEMNoDataValue,
                    externalDEMApplyEGM=externalDEMApplyEGM)

    ############################################
    # write workflow to file and optionally execute it
    log.debug('writing workflow to file')

    tmp_name = outname.replace(tmpdir, outdir)
    wf_name = tmp_name + '_proc.xml'
    workflow.write(wf_name)

    # execute the newly written workflow
    if not test:
        #try:
        if True:
            groups=None
            if groupsize is not None:
                groups = groupbyWorkers(wf_name, groupsize)
            gpt(wf_name, groups=groups, cleanup=cleanup, tmpdir=outname,
                gpt_exceptions=gpt_exceptions, gpt_args=gpt_args,
                removeS1BorderNoiseMethod=removeS1BorderNoiseMethod)
            writer(xmlfile=wf_name, outdir=outdir, basename_extensions=basename_extensions,
                   clean_edges=clean_edges, clean_edges_npixels=clean_edges_npixels)
        #except Exception as e:
        #    log.info(str(e))
        #    with open(wf_name.replace('_proc.xml', '_error.log'), 'w') as logfile:
        #        logfile.write(str(e))
        #finally:
        #    if cleanup and os.path.isdir(outname):
        #        log.info('deleting temporary files')
        #        shutil.rmtree(outname, onerror=windows_fileprefix)
        #log.info('done')


            
    
def insar_coherence(infiles, swaths=["IW1","IW2","IW3"], polarizations='all', t_srs=4326, demName='SRTM 1Sec HGT',
                    allow_RES_OSV=False, 
                    demResamplingMethodBGC="BICUBIC_INTERPOLATION", maskOutAreaWithoutElevation=False,
                    externalDEMFile=None, externalDEMNoDataValue=None, externalDEMApplyEGM=True,
                    demResamplingMethod='BILINEAR_INTERPOLATION', imgResamplingMethod='BILINEAR_INTERPOLATION',
                    alignToStandardGrid=False, standardGridOriginX=0, standardGridOriginY=0,
                    gpt_exceptions=None, gpt_args=None, rlks=None, azlks=None, cohWinRg=5, cohWinAz=5, spacing=20,
                    cleanup=True, tmpdir=None, outName=None, test=False, enhancedSpectralDiversity=True,
                    clean_edges=False, clean_edges_npixels=1, outdir=None, basename_extensions=None,
                    returnWF=False,groupsize=2, removeS1BorderNoiseMethod='pyroSAR',
                    dynamic_cleaning=True):

    if not isinstance(infiles, list):
        raise RuntimeError("'infiles' must be of type list")
    if len(infiles)!=2:
        raise RuntimeError("Two files are mandatory to compute Insar coherence, less or more were provided")
    ids = identify_many(infiles, sortkey='start')
    id_1 = ids[0]
    id_2 = ids[1]

    if id_1.is_processed(outdir):
        log.info(f'scene {id_1.outname_base()} already processed')
        return

    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    
    ############################################
    # general setup
    process_S1_SLC = False
    if id_1.sensor not in ['S1A', 'S1B', 'S1C'] or id_1.product!="SLC" or\
       id_2.sensor not in ['S1A', 'S1B', 'S1C'] or id_2.product!="SLC" :
        raise RuntimeError('Insar coherence only available for Sentinel mission in SLC mode')

    formatName = 'SENTINEL-1'

    if id_1.acquisition_mode != 'IW' or id_2.acquisition_mode != 'IW':
        raise RuntimeError(f"acquisition mode {id_1.acquisition_mode}/{id_2.acquisition_mode} not supported")

    #check s1 frame compatibility
    # second image has to be the first image used for coherence estimation
    delta_t = (dt.datetime.strptime(id_2.start,"%Y%m%dT%H%M%S")-dt.datetime.strptime(id_1.start,"%Y%m%dT%H%M%S")).total_seconds()-12*86400
    #delta_day = (dt.datetime.strptime(id_1.start,"%Y%m%dT%H%M%S")-dt.datetime.strptime(id_2.start,"%Y%m%dT%H%M%S")).days
    #delta_sec_p = abs((dt.datetime.strptime(id_1.start,"%Y%m%dT%H%M%S")-dt.datetime.strptime(id_2.start,"%Y%m%dT%H%M%S")).seconds)
    #delta_sec_m = abs((dt.datetime.strptime(id_2.start,"%Y%m%dT%H%M%S")-dt.datetime.strptime(id_1.start,"%Y%m%dT%H%M%S")).seconds)
    if id_1.orbitNumber_rel!=id_2.orbitNumber_rel or\
       id_1.orbit!=id_2.orbit or\
       id_1.sensor!=id_2.sensor or \
       abs(delta_t)>10:
       #(delta_day!=12 and not (delta_day==11 and delta_sec_p>86345)) or
       #min(delta_sec_p, delta_sec_m)>10:
       
       import logging
       logger = logging.getLogger("my_logger")
       logger.error(delta_t)
       raise RuntimeError(f"Invalid S1 frame comparison, wrong geographical matching")
    ######################
    if isinstance(polarizations, str):
        if polarizations == 'all':
            polarizations = id_1.polarizations
        else:
            if polarizations in id_1.polarizations:
                polarizations = [polarizations]
            else:
                raise RuntimeError('polarization {} does not exists in the source product'.format(polarizations))
    elif isinstance(polarizations, list):
        polarizations = [x for x in polarizations if x in id.polarizations]
    else:
        raise RuntimeError('polarizations must be of type str or list')


    if externalDEMFile is None and externalDEMNoDataValue is None:
        externalDEMNoDataValue=0.

    if tmpdir is None:
        tmpdir = outdir
        
    ############################################

    for pol in polarizations:
        for iw in swaths:

            # parse base workflow
            workflow = parse_recipe('blank')
            ############################################
            
            read_iw_orb_ids=[]
            for i in range(0, len(infiles)):
                ############################################
                # Read node configuration
                read = parse_node('Read')
                workflow.insert_node(read)
                read.parameters['file'] = ids[i].scene
                read.parameters['formatName'] = formatName

                ############################################
                # TOP-SAR split
                ts_split = parse_node('TOPSAR-Split')
                workflow.insert_node(ts_split, before=read.id)
                ts_split.parameters["subswath"] = iw
                ts_split.parameters["selectedPolarisations"] = [pol]

                ############################################
                # Apply-Orbit-File node configuration in coherence estimation
                orb = orb_parametrize(scene=ids[i], formatName=formatName, allow_RES_OSV=allow_RES_OSV)
                workflow.insert_node(orb, before=ts_split.id)
                read_iw_orb_ids.append(orb.id)
                
            ############################################
            # Back geocoding node configuration for coherence estimation
            bgc = parse_node("Back-Geocoding")
            workflow.insert_node(bgc, before=read_iw_orb_ids)
            bgc.parameters["demName"]= demName
            bgc.parameters["demResamplingMethod"]= demResamplingMethodBGC
            bgc.parameters["externalDEMFile"]= externalDEMFile
            bgc.parameters["externalDEMNoDataValue"]= externalDEMNoDataValue
            bgc.parameters["resamplingType"]= "BISINC_5_POINT_INTERPOLATION"
            bgc.parameters["maskOutAreaWithoutElevation"]=maskOutAreaWithoutElevation
            last = bgc
            ############################################
            # Enhanced spectral diversity node configuration
            if enhancedSpectralDiversity:
                esd = parse_node("Enhanced-Spectral-Diversity")
                workflow.insert_node(esd, before=last.id)
                last = esd
        
            ############################################
            # Coherence estimation node configuration
            coh = parse_node("Coherence")
            workflow.insert_node(coh, before=last.id)
            coh.parameters["subtractFlatEarthPhase"]= True
            coh.parameters["singleMaster"]= True
            coh.parameters["cohWinRg"]= cohWinRg
            coh.parameters["cohWinAz"]= cohWinAz
            coh.parameters["demName"]= demName
            coh.parameters["subtractTopographicPhase"]= True
            coh.parameters["externalDEMFile"]= externalDEMFile
            coh.parameters["externalDEMNoDataValue"]= externalDEMNoDataValue
            coh.parameters["externalDEMApplyEGM"]= True
            last = coh
        
            ############################################
            deb = parse_node('TOPSAR-Deburst')
            workflow.insert_node(deb, before=last.id)
            deb.parameters["selectedPolarisations"]=[pol]

            ############################################
            ##create out_name and write

            ##extract dates as str from filename for the day and the full datetime
            date1= ids[0].start.split("T")[0]
            date2= ids[1].start.split("T")[0]
            
            tmp_out_name= f"S1_{iw}_COH_{pol}_{date1}_{date2}_tmp"
            tmp_out= os.path.join(tmpdir, tmp_out_name)

            write_coh=parse_node("Write")
            workflow.insert_node(write_coh, before= deb.id)
            write_coh.parameters["file"]= tmp_out
            write_coh.parameters["formatName"]= "BEAM-DIMAP"

            tmp_wf_name = f"{tmpdir}/Coh_tmp_prep_graph_{pol}_{iw}.xml"
            workflow.write(tmp_wf_name)
            if not test:
                gpt(tmp_wf_name, groups=None, cleanup=False, tmpdir=tmpdir,
                    gpt_exceptions=gpt_exceptions, gpt_args=gpt_args,
                    removeS1BorderNoiseMethod=removeS1BorderNoiseMethod)

        #sys.exit(0)
        #back to all iws, single polarization
        import glob
        tmp_coh_files= glob.glob(f"{tmpdir}/S1_*_COH_{pol}_{date1}_{date2}_tmp.dim")
        print(tmp_coh_files)

        workflow = parse_recipe('blank')
        reads=[]
        for cohf in tmp_coh_files:
            ############################################
            # Read node configuration
            read = parse_node('Read')
            workflow.insert_node(read, resetSuccessorSource=False)
            read.parameters['file'] = cohf
            read.parameters['formatName'] = "BEAM-DIMAP"
            reads.append(read.id)

        ############################################
        # merge sub-swaths node configuration
        merge=parse_node("TOPSAR-Merge")
        merge.parameters["selectedPolarisations"]=pol
        workflow.insert_node(merge, before=reads)
        last = merge

        ############################################
        # Multilook node configuration
        ml = mli_parametrize(scene=id_1, spacing=spacing, rlks=rlks, azlks=azlks)
        if ml is not None:
            workflow.insert_node(ml, before=last.id)
            last = ml
        ############################################
        tc = geo_parametrize(spacing=spacing, t_srs=t_srs, demName=demName,
                             externalDEMFile=externalDEMFile,
                         externalDEMNoDataValue=externalDEMNoDataValue,
                         externalDEMApplyEGM=externalDEMApplyEGM,
                         alignToStandardGrid=alignToStandardGrid,
                         standardGridOriginX=standardGridOriginX,
                         standardGridOriginY=standardGridOriginY)
        workflow.insert_node(tc, before=last.id)
        last = tc
        ############################################

        ############################################
        # parametrize write node
        # create a suffix for the output file to identify processing steps performed in the workflow
    
        suffix = workflow.suffix()
        if tmpdir is None:
            tmpdir = outdir
        basename = os.path.join(tmpdir, id_1.outname_base(basename_extensions))
        outname = f"{basename}_{pol}_{suffix}"
    
        write = parse_node('Write')
        workflow.insert_node(write, before=last.id)
        write.parameters['file'] = outname
        write.parameters['formatName'] = 'ENVI'
        ############################################
    
        ############################################
        # write workflow to file and optionally execute it
        log.debug('writing workflow to file')

        tmp_name = outname.replace(tmpdir, outdir)
        wf_name = tmp_name + '_proc.xml'
        workflow.write(wf_name)
    
        # execute the newly written workflow
        if not test:
            try:
                gpt(wf_name, groups=None, cleanup=cleanup, tmpdir=outname,
                    gpt_exceptions=gpt_exceptions, gpt_args=gpt_args,
                    removeS1BorderNoiseMethod=removeS1BorderNoiseMethod)
                writer(xmlfile=wf_name, outdir=outdir, basename_extensions=basename_extensions,
                       clean_edges=clean_edges, clean_edges_npixels=clean_edges_npixels)
            except Exception as e:
                log.info(str(e))
                with open(wf_name.replace('_proc.xml', '_error.log'), 'w') as logfile:
                    logfile.write(str(e))
            finally:
                if cleanup and os.path.isdir(outname):
                    log.info('deleting temporary files')
                    shutil.rmtree(outname, onerror=windows_fileprefix)
            log.info('done')


        

"""
            
    coh_iws=[]
    for iw in swaths:

        read_iw_orb_ids=[]
        #reads=[]
        for i in range(0, len(infiles)):
            ############################################
            # Read node configuration
            read = parse_node('Read')
            workflow.insert_node(read)
            read.parameters['file'] = ids[i].scene
            read.parameters['formatName'] = formatName
            #reads.append(read)

        #for read in reads: #FIXME MM?
            ############################################
            # TOP-SAR split
            ts_split = parse_node('TOPSAR-Split')
            workflow.insert_node(ts_split, before=read.id)#, resetSuccessorSource=False)
            ts_split.parameters["subswath"] = iw
            ts_split.parameters["selectedPolarisations"] = polarizations

            ############################################
            # Apply-Orbit-File node configuration in coherence estimation
            orb = orb_parametrize(scene=ids[i], formatName=formatName, allow_RES_OSV=allow_RES_OSV)
            workflow.insert_node(orb, before=ts_split.id)
            
            read_iw_orb_ids.append(orb.id)

        ############################################
        # Back geocoding node configuration for coherence estimation
        bgc = parse_node("Back-Geocoding")
        workflow.insert_node(bgc, before=read_iw_orb_ids)
        bgc.parameters["demName"]= demName
        bgc.parameters["demResamplingMethod"]= demResamplingMethodBGC
        bgc.parameters["externalDEMFile"]= externalDEMFile
        bgc.parameters["externalDEMNoDataValue"]= externalDEMNoDataValue
        bgc.parameters["resamplingType"]= "BISINC_5_POINT_INTERPOLATION"
        bgc.parameters["maskOutAreaWithoutElevation"]=maskOutAreaWithoutElevation
        last = bgc
        ############################################
        # Enhanced spectral diversity node configuration
        if enhancedSpectralDiversity:
            esd = parse_node("Enhanced-Spectral-Diversity")
            workflow.insert_node(esd, before=last.id)
            last = esd
        
        ############################################
        # Coherence estimation node configuration
        coh = parse_node("Coherence")
        workflow.insert_node(coh, before=last.id)
        coh.parameters["subtractFlatEarthPhase"]= True
        coh.parameters["singleMaster"]= True
        coh.parameters["cohWinRg"]= cohWinRg
        coh.parameters["cohWinAz"]= cohWinAz
        coh.parameters["demName"]= demName
        coh.parameters["subtractTopographicPhase"]= True
        coh.parameters["externalDEMFile"]= externalDEMFile
        coh.parameters["externalDEMNoDataValue"]= externalDEMNoDataValue
        coh.parameters["externalDEMApplyEGM"]= True
        last = coh
        
        ############################################
        deb = parse_node('TOPSAR-Deburst')
        workflow.insert_node(deb, before=last.id)
        deb.parameters["selectedPolarisations"]=polarizations
        coh_iws.append(deb.id)

    ############################################
    # merge sub-swaths node configuration
    merge=parse_node("TOPSAR-Merge")
    merge.parameters["selectedPolarisations"]=polarizations
    workflow.insert_node(merge, before=coh_iws)
    last = merge

    ############################################
    # Multilook node configuration
    ml = mli_parametrize(scene=id_1, spacing=spacing, rlks=rlks, azlks=azlks)
    if ml is not None:
        workflow.insert_node(ml, before=last.id)
        last = ml
    ############################################
    tc = geo_parametrize(spacing=spacing, t_srs=t_srs, demName=demName,
                         externalDEMFile=externalDEMFile,
                         externalDEMNoDataValue=externalDEMNoDataValue,
                         externalDEMApplyEGM=externalDEMApplyEGM,
                         alignToStandardGrid=alignToStandardGrid,
                         standardGridOriginX=standardGridOriginX,
                         standardGridOriginY=standardGridOriginY)
    workflow.insert_node(tc, before=last.id)
    last = tc
    ############################################

    ############################################
    # parametrize write node
    # create a suffix for the output file to identify processing steps performed in the workflow
    
    suffix = workflow.suffix()
    if tmpdir is None:
        tmpdir = outdir
    basename = os.path.join(tmpdir, id_1.outname_base(basename_extensions))
    outname = basename + '_' + suffix

    print(basename, "--", suffix, "--", outname)
    
    write = parse_node('Write')
    workflow.insert_node(write, before=last.id)
    write.parameters['file'] = outname
    write.parameters['formatName'] = 'ENVI'
    ############################################
    
    ############################################
    # write workflow to file and optionally execute it
    log.debug('writing workflow to file')

    tmp_name = outname.replace(tmpdir, outdir)
    wf_name = tmp_name + '_proc.xml'
    workflow.write(wf_name)
    print(wf_name)
    print(groupbyWorkers(wf_name, groupsize))
    
    # execute the newly written workflow
    if not test:
        try:
            p=None
            if groupsize>1 and dynamic_cleaning:
                tmp_name +="/sub"
                p = Popen(['python', 'dynamicCOHCleaning.py', tmp_name])

            groups = groupbyWorkers(wf_name, groupsize)
            gpt(wf_name, groups=groups, cleanup=cleanup, tmpdir=outname,
                gpt_exceptions=gpt_exceptions, gpt_args=gpt_args,
                removeS1BorderNoiseMethod=removeS1BorderNoiseMethod)
            writer(xmlfile=wf_name, outdir=outdir, basename_extensions=basename_extensions,
                   clean_edges=clean_edges, clean_edges_npixels=clean_edges_npixels)
            if p:
                p.terminate()
        except Exception as e:
            log.info(str(e))
            with open(wf_name.replace('_proc.xml', '_error.log'), 'w') as logfile:
                logfile.write(str(e))
        finally:
            if cleanup and os.path.isdir(outname):
                log.info('deleting temporary files')
                shutil.rmtree(outname, onerror=windows_fileprefix)
        log.info('done')
    


    if returnWF:
        return wf_name
"""
