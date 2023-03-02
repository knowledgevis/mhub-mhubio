"""
-------------------------------------------------
MHub - File type class for mhubio instance data
-------------------------------------------------

-------------------------------------------------
Author: Leonard Nürnberg (27.02.2023)
Email:  leonard.nuernberg@maastrichtuniversity.nl
-------------------------------------------------
"""

from enum import Enum

class FileType(Enum):
    NONE        = None
    NRRD        = "nrrd"
    NIFTI       = "nifti"
    DICOM       = "dicom"
    DICOMSEG    = "dicomseg"
    RTSTRUCT    = "rtstruct"
    CSV         = "csv"
    TXT         = "txt"
    JSON        = "json"

    def __str__(self) -> str:
        return self.name