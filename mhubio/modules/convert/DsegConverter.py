"""
-------------------------------------------------
MHub - DicomSeg Conversion Module
-------------------------------------------------

-------------------------------------------------
Author: Leonard Nürnberg
Email:  leonard.nuernberg@maastrichtuniversity.nl
-------------------------------------------------
"""

from typing import List
from mhubio.core import Module, Instance, InstanceData, InstanceDataCollection, DataType, IO
import os, subprocess

from segdb.tools import DcmqiDsegConfigGenerator

@IO.Config('source_segs', List[DataType], ['nifti:mod=seg:roi=*', 'nrrd:mod=seg:roi=*'], factory=IO.F.list(DataType.fromString), the='target segmentation files to convert to dicomseg')
@IO.Config('target_dicom', DataType, 'dicom:mod=ct', factory=DataType.fromString, the='dicom data all segmentations align to')
@IO.Config('skip_empty_slices', bool, True, the='flag to skip empty slices')
@IO.Config('converted_file_name', str, 'seg.dcm', the='name of the converted file')
@IO.Config('bundle_name', str, None, the="bundle name converted data will be added to")
@IO.Config('model_name', str, 'MHub-Model', the="model name populated in the dicom seg SeriesDescription attribute")
@IO.Config('json_config_path', str, None, the='path to the dicomseg json config file')    
@IO.Config('segment_id_meta_key', str, 'roi', the='meta key used to identify the roi in the dicomseg json config file')
class DsegConverter(Module):

    source_segs: List[DataType]
    target_dicom: DataType
    skip_empty_slices: bool
    converted_file_name: str
    bundle_name: str
    model_name: str
    json_config_path: str
    segment_id_meta_key: str

    @IO.Instance()
    @IO.Inputs("in_segs", IO.C("source_segs"), the="input data to convert to dicomseg")
    @IO.Input("in_dicom", IO.C("target_dicom"), the="input dicom data to convert to dicomseg")
    @IO.Output('out_data', path=IO.C('converted_file_name'), dtype='dicomseg:mod=seg', data='in_dicom', bundle=IO.C('bundle_name'), auto_increment=True, the="converted data")
    def task(self, instance: Instance, in_segs: InstanceDataCollection, in_dicom: InstanceData, out_data: InstanceData) -> None:
        
        # either use a custom json config or generate based on meta label (default key: roi)
        if self.json_config_path is not None:

            # sort files alphabetically
            file_list = sorted([in_seg.abspath for in_seg in in_segs])   
            json_config_path = self.json_config_path

        else:

            # generate json meta generator instance
            generator = DcmqiDsegConfigGenerator(
                model_name = self.model_name
            )

            # extract and populate data 
            for in_seg in in_segs:
                generator.addItem(
                    file = in_seg.abspath, 
                    segment_ids = in_seg.type.meta[self.segment_id_meta_key].split(","),
                    model_name = in_seg.type.meta.getValue('model'),
                )

            # store json in temp dir
            tmp_dir = self.config.data.requestTempDir('dseg_converter')
            json_config_path = os.path.join(tmp_dir, "temp-meta.json")
            generator.save(config_file=json_config_path, overwrite=True)

            # get file list (comma separated string)
            file_list = ",".join(generator.getFileList())

            # create outdir if required
            # TODO: can we handle this during bundle creation or in IO decorator?
            if not os.path.isdir(os.path.dirname(out_data.abspath)):
                os.makedirs(os.path.dirname(out_data.abspath))

        # build command
        bash_command  = ["itkimage2segimage"]
        bash_command += ["--inputImageList", file_list]
        bash_command += ["--inputDICOMDirectory", in_dicom.abspath]
        bash_command += ["--outputDICOM", out_data.abspath]
        bash_command += ["--inputMetadata", json_config_path]

        # add skip empty slices flag
        if self.c["skip_empty_slices"] == True:
            bash_command += ["--skip"]

        # run command
        try: 
            self.v(">> run: ", " ".join(bash_command))
            _ = subprocess.run(bash_command, check = True, text = True)
        except Exception as e:
            self.v("Error while running dicomseg-conversion for instance " + str(instance) + ": " + str(e))