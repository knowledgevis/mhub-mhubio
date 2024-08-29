"""
-------------------------------------------------
MHub - WSI DICOM Importer Module

Description:  Import a selected resolution from a multi-resolution DICOM WSI image. 
This is based off of the DicomImporter module, but with the focus on DICOM-WSI only. 
-------------------------------------------------

-------------------------------------------------
Author: Leonard Nürnberg, Curtis Lisle
Email:  leonard.nuernberg@maastrichtuniversity.nl
Email:  curtislisle@knowledgevis.com
-------------------------------------------------
"""

import os
import shutil
from glob import glob

from enum import Enum
from typing import List
from mhubio.core import Module, Instance, InstanceData, DataType, Meta, FileType, DirectoryChain, IO, SEG
import pydicom

class InputDirStructure(Enum):
    """
    Input directory structure options.
    """
    FLAT = 1
    SERIES = 2
    UNKNOWN = 3

@IO.Config('source_dir', str, 'input_data', the="source input directory containing the (unsorted) dicom data")
@IO.Config('import_dir', str, 'sorted_data', the="output directory where the imported (sorted / organized) dicom data will be placed")
@IO.Config('import_data', bool, True, the="flag to enable the import process")
@IO.Config('structure', str, "%SeriesInstanceUID/dicom/%SOPInstanceUID.dcm", the="schema used for sorting the dicom data")
@IO.Config('magnification', int, 10, the="the desired magnification level of the output WSI image")
@IO.Config('magTolerance', float, 2.0, the="the limit of the mag difference allows in the input WSI image")
@IO.Config('meta', dict, {'mod': '%Modality'}, the="meta data used for every imported instance")
class WsiDicomImporter(Module):
    """
    Import dicom data into instances.  Allow two standard organizations (flat or series of instances in separate directories)
    For now, the static schema is: %SeriesInstanceUID/dicom/%SOPInstanceUID.dcm

    To override the metadata via a custom config, pass a json string dictionary:
     ... '--config:modules.DicomImporter.meta={"mod": "MR"}'
    """

    source_dir: str
    import_dir: str
    sort_data: bool
    import_data: bool
    structure: str
    magnification: int 
    magTolerance: float
    meta: dict

    def updateMeta(self, dicom_data: InstanceData) -> None:

        if not any(v.startswith('%') for v in dicom_data.type.meta.values()):
            return

        # pick first file
        dicom_file = os.listdir(dicom_data.abspath)[0]
        dicom_file_path = os.path.join(dicom_data.abspath, dicom_file)

        # load dicom meta data
        ds = pydicom.read_file(dicom_file_path)

        # update meta (lookup dicom placeholders starting with %)
        meta_update = {}
        for k, v in dicom_data.type.meta.items():
            if v.startswith('%'):
                dicom_field = v[1:]
                if hasattr(ds, dicom_field):
                    meta_update[k] = getattr(ds, dicom_field)
                else:
                    self.v(f">> dicom field not found: {dicom_field}")
                    meta_update[k] = "" # empty string as placeholder

        # update meta of the dicom data
        dicom_data.type.meta += meta_update

    

    def scanSourceDir(self, input_dir: str) -> InputDirStructure:

        # check for only dicom files / only folders
        hasOnlyDicomFiles: bool = True
        hasOnlyFolders: bool = True

        # scan first level of input dir
        for d in os.listdir(input_dir):
            if not os.path.isfile(os.path.join(input_dir, d)) or not d.endswith(".dcm"):
                hasOnlyDicomFiles = False
            if not os.path.isdir(os.path.join(input_dir, d)):
                hasOnlyFolders = False

        # return
        if hasOnlyDicomFiles:
            return InputDirStructure.FLAT
        elif hasOnlyFolders: 
            return InputDirStructure.SERIES
        else:
            self.v("unknown directory structure: OnlyDicomFiles: ", hasOnlyDicomFiles, " OnlyFolders: ", hasOnlyFolders)
            return InputDirStructure.UNKNOWN


    # ---- begin of support dicom WSI directory analysis routines ----

    # extract header information from a dicom whole slide image file
    def extract_image_info_dicom(self,dicom_file):
        correctExtraction = True
        try:
            ds = pydicom.dcmread(dicom_file,stop_before_pixels=True)
            ObjectiveLensPower = float(ds[0x0048,0x0105][0][0x0048,0x0112].value)
            TotalColumns = ds.TotalPixelMatrixColumns
            TotalRows = ds.TotalPixelMatrixRows
            ImageSize = TotalRows * TotalColumns
            ImageColumns = ds.Columns
            ImageRows = ds.Rows
            return correctExtraction, ImageSize,ObjectiveLensPower,TotalColumns,TotalRows,ImageColumns,ImageRows
        except:
            print('error extracting image metadata from:',dicom_file)
            correctExtraction = False
            return correctExtraction, None,None,None,None,None,None

    # this routine looks through all dicom images in a directory and extracts some basic information
    # about the images.  It returns a list of dictionaries with the filename, in-file magnification 
    # (which always lists the highest magnification in the stack), columns and rows

    def pydicom_extract_information(self,basedir):
        os.chdir(basedir)
        imageList = []
        goodImageCount = 0
        for filename in glob('*dcm'):
            print(filename)
            #print_image_info_dicom(filename)
            success, size,mag,cols,rows,icols,irows = self.extract_image_info_dicom(filename)
            if success:
                rec = {'filename': filename,'size':size, 'mag': mag, 'cols': cols, 'rows': rows, 'icols': icols, 'irows': irows}
                imageList.append(rec)
                goodImageCount += 1
        print('good image count:',goodImageCount)
        return imageList

    # post process the image list to determine the image with the highest magnification. Return only the
    # multi-tiled images that are part of the multi-tiled image stack and calculate their true magnification.
    # NOTE: this assumes square tiles

    def analyze_image_list(self,imageList):
        # first sort the list by total number of pixels, so the largest file is first
        imageList.sort(key=lambda x: x['size'],reverse=True)
        # now that the images are sorted by size, we look for ones with the same image rows&cols size
        # because these are multi-tiled images
        tiledImages = []
        for i in range(len(imageList)):
            # if the image tile is square, it is a multi-tiled image
            if (imageList[i]['icols'] == imageList[i]['irows'] ):
                #print('multi-tiled image:',imageList[i]['filename'])
                #print('total pixels:',imageList[i]['size'])
                magnification = imageList[0]['mag']*(imageList[i]['cols']/imageList[0]['cols'])
                #print('magnification:',magnification)
                rec = {'filename': imageList[i]['filename'],'size':imageList[i]['size'], 'mag': magnification, 'cols': imageList[i]['cols'], 'rows': imageList[i]['rows'], 'icols': imageList[i]['icols'], 'irows': imageList[i]['irows']}
                tiledImages.append(rec)
        return tiledImages

    # if there is an image with a magnification within default tolerance = 2.0 
    # of the requested magnification, return it
    def return_magnitude_match(self,imageList,targetMag,magTolerance=2.0):
        for i in range(len(imageList)):     
            if (abs(float(imageList[i]['mag']) - float(targetMag)) < magTolerance):
                self.v('found image with magnification near target:',targetMag)
                self.v('filename:',imageList[i]['filename'])
                #print('total pixels:',imageList[i]['size'])
                self.v('actual magnification:',imageList[i]['mag'])
                return imageList[i]
        return None

# ---- end of support dicom WSI directory analysis routines ----

    def importSingleInstance(self, input_dir: str, sorted_dir: str) -> None:
    
        # verbose
        self.v("> importing single instance")

        # create new instance
        #  as we just have one instance, we use the base folder
        instance = Instance(sorted_dir)
        instance.attr['sid'] = 'inst0' 

        # append instance to data handler to resolve dc chains corectly / automatically
        self.config.data.addInstance(instance)

        # decide which image to import, there are images of different resolution in this directory
        imageInfo = self.pydicom_extract_information(input_dir)
        tiledImages = self.analyze_image_list(imageInfo)
        self.v('found ',len(tiledImages),'tiled images in this directory')
        matchingImage = self.return_magnitude_match(tiledImages,self.magnification,self.magTolerance)
        if matchingImage is None:
            self.v('no image matching the target magnification found')
        else:
            self.v(matchingImage['filename'])

            # create dicom data
            dicom_data_meta = self.meta.copy()
            dicom_data_type = DataType(FileType.DICOM, dicom_data_meta)
            dicom_data = InstanceData('dicom', dicom_data_type, instance)

            # copy the dicom data
            # For the radiology data, we copy the entire instance root directory. However, with
            # pathology data (dicom-wsi), the files can be extremely large. In this case, we only 
            # copy the input image that matches the desired magnification level.

            #shutil.copytree(input_dir, dicom_data.abspath)
            if not os.path.exists(dicom_data.abspath):
                os.makedirs(dicom_data.abspath)
            shutil.copyfile(os.path.join(input_dir,matchingImage['filename']),os.path.join(dicom_data.abspath,matchingImage['filename']))

            # update meta with dicom placeholders
            self.updateMeta(dicom_data)

            # confirm data is where we expect it to be
            if os.path.isdir(dicom_data.abspath):
                dicom_data.confirm()

    def importMultipleInstances(self, input_dir: str, sorted_dir: str) -> None:
        
        # iterate sorted data and generate a new instance with the dicom data 
        instances_n = len(os.listdir(input_dir))
        for i, sid in enumerate(os.listdir(input_dir)):
            self.v(f"> importing instance ({i+1}/{instances_n}): ", sid)
            
            # create new instance
            instance = Instance(os.path.join(sorted_dir, sid))
            instance.attr['sid'] = sid
            
            # append instance to data handler to resolve dc chains corectly / automatically
            self.config.data.addInstance(instance)

            # decide which image to import, there are images of different resolution in this directory
            instance_dir = os.path.join(input_dir, sid)
            imageInfo = self.pydicom_extract_information(instance_dir)
            tiledImages = self.analyze_image_list(imageInfo)
            self.v('found ',len(tiledImages),'tiled images in this directory')
            matchingImage = self.return_magnitude_match(tiledImages,self.magnification,self.magTolerance)
            if matchingImage is None:
                self.v('no image matching the desired target magnification was found in directory',input_dir)
            else:
                self.v('found matching image to import:',matchingImage['filename'])

            # create dicom data
            dicom_data_meta = self.meta.copy()
            dicom_data_type = DataType(FileType.DICOM, dicom_data_meta)
            dicom_data = InstanceData('dicom', dicom_data_type, instance)

            # copy the dicom data
            # For the radiology data, we copy the entire instance root directory. However, with
            # pathology data (dicom-wsi), the files can be extremely large. In this case, we only 
            # copy the input image that matches the desired magnification level.

            if not os.path.exists(dicom_data.abspath):
                os.makedirs(dicom_data.abspath)
            shutil.copyfile(os.path.join(input_dir,os.path.join(instance_dir,matchingImage['filename'])),os.path.join(dicom_data.abspath,matchingImage['filename']))

            # update meta with dicom placeholders
            self.updateMeta(dicom_data)

            # confirm 
            if os.path.isdir(dicom_data.abspath):
                dicom_data.confirm()

    
    def task(self) -> None:

        self.v('> WsiDicomImporter task')
        # resolve the input directory
        source_dc = DirectoryChain(path=self.source_dir, parent=self.config.data.dc)
        import_dc = DirectoryChain(path=self.import_dir, parent=self.config.data.dc)

        self.v("> source input dir: ", self.source_dir, " --> ", source_dc.abspath)
        self.v("> import sort  dir: ", self.import_dir, " --> ", import_dc.abspath)

        # import the data
        self.v()
        self.v("> importing dicom whole slide imaging data")

        # scan input structure
        source_dir_structure = self.scanSourceDir(source_dc.abspath)
        print("> source dir structure: ", source_dir_structure)

        #self.importSingleInstance(source_dc.abspath, import_dc.abspath) 

        # import data depending on input structure
        if source_dir_structure == InputDirStructure.FLAT:
            self.importSingleInstance(source_dc.abspath, import_dc.abspath) # Note: import_dc.abspath overrides instance with an abspath // self.import_dir instead?
        elif source_dir_structure == InputDirStructure.SERIES:
            self.importMultipleInstances(source_dc.abspath, import_dc.abspath)
        else:
            raise ValueError("Error: input directory structure is unknown. Cannot determine if this is a single series or multiple series.")
        
