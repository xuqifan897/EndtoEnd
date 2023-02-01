#include <stdio.h>
#include <iostream>
#include "./rtimages.h"
#include "./rtstruct.h"

void print_usage(char* argv[]) {
    printf("Anonymize and save new Dicom series\n\n"
           "  Usage\n"
           "    %s <indir> <outdir>\n\n"
           "    indir:  path to dir containing dicom image files\n"
           "    outdir: path to dir where new dicom image files will be saved"
           , argv[0]);
}

int main(int argc, char* argv[]) {
    if (argc != 3) {
        print_usage(argv);
    }
    char* indir = argv[1];
    char* outdir = argv[2];
    RTImage *Image = new RTImage;

    Image->setDicomDirectory(indir);
    if (!Image->loadDicomInfo()) {
        RTStruct rtstruct = RTStruct();
        rtstruct.setDicomDirectory(indir);
        rtstruct.loadDicomInfo();
        rtstruct.loadRTStructInfo();
        rtstruct.saveRTStructData(outdir, true);
    } else {
        Image->loadRTImageData();
        Image->saveRTImageData(outdir, true);
    }

    return 0;
}