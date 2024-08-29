"""
-------------------------------------------------
MHub - Select Whole Slide Image Magnitude  Filter Module.
-------------------------------------------------

-------------------------------------------------
Author: Curtis Lisle
Email:  curtislisle@knowledgevis.com
-------------------------------------------------
"""

import os, uuid, json
from mhubio.core import Module, Instance, InstanceDataCollection, InstanceData, IO, DataType, FileType, Meta, InstanceDataBundle, SEG
from typing import List, Union
import shutil
from glob import glob
import pydicom


@IO.Config("target_magnification", float, 10.0, the="image magnitude to filter on")
@IO.Config("magnification_tolerance", float, 1.0, the="tolerance for image magnitude to allow")
@IO.ConfigInput('target_dicom', 'dicom:mod=sm', the='Dicom WSI image file.')
@IO.ConfigInput('in_datas', 'dicom:mod=sm', the='Dicom segmentation file.')
@IO.Config('bundle', str, 'extracted', the='Bundle name converted data will be added to')
class WsiMagnificationExtractor(Module):

    target_magnification: float
    magnification_tolerance: float
    bundle: str

    @IO.Instance()
    @IO.Input("target_dicom", the="input dicom data to extract files from")
    @IO.Inputs("in_datas", the="wsi files related to the target dicom that will be examined and filtered")
    @IO.Outputs("out_datas", "extracted", "dicom:mod=sm")
    def task(self, instance: Instance, target_dicom: InstanceData, in_datas: InstanceDataCollection, out_datas: InstanceDataCollection) -> None:
        
        # create bundle for output data
        bundle = instance.getDataBundle(self.bundle)
        bundle.dc.makedirs(is_file=False)
        
        # check the input data directory
        print('input directory:',target_dicom.abspath)
     
        # iterate through the datasets in this instance
        for in_data in in_datas:
            print('found input data:',in_data.abspath)
            self.process_data(in_data, bundle, out_datas)
        


            
    def process_data(self, in_data: InstanceData, bundle: InstanceDataBundle, out_datas: InstanceDataCollection) -> None:
            
        # log
        self.log("in_data", in_data)
        self.log("bundle", bundle)
        self.log("bundle path: ", bundle.abspath, os.path.exists(bundle.abspath))
        
        # random 10 digit run id
        run_id = uuid.uuid4().hex[:10]
        
        imageInfo = self.pydicom_extract_information(in_data.abspath)
        tiledImages = self.analyze_image_list(imageInfo)
        self.v('found ',len(tiledImages),'tiled images in this directory')
        print('tiled images:',tiledImages)
        matchingImage = self.return_magnitude_match(tiledImages,self.target_magnification,self.magnification_tolerance)
        # if we found a matching image, we will include it in the output data
        if matchingImage: 
            wsiFile = matchingImage['filename']
            print('found matching image:',wsiFile)
            # we are generating output metadata for the WSI file that will pass through
            out_data_type = DataType(FileType.DICOM, Meta(origin="sm"))

            # specify the output data file in the output instance
            out_data = InstanceData(
                path=os.path.join(bundle.abspath, wsiFile),
                type=out_data_type,
                bundle=bundle,
                auto_increment=True
            )
            print('out_data path:',out_data.abspath)

            # copy the input data file to the output data file
            outbasename = os.path.basename(out_data.abspath)
            outdirname = os.path.dirname(out_data.abspath)
            if not os.path.exists(outdirname):
                os.makedirs(outdirname)
            shutil.copyfile(wsiFile,os.path.join(outdirname,outbasename))

            # add output data file to collection
            out_datas.add(out_data)





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