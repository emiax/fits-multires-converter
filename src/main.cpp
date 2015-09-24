#include <CCfits>
#include <iostream>
#include <fitsconverter.h>
#include <lz4compressor.h>

int main() {

    // Create lz4 compressed output
    {
        FitsConverter fc;
        fc.setInFolder("input/");
        fc.setOutFile("two_compressed.data");

        Compressor* compressor = new Lz4Compressor();

        fc.setCompressor(compressor);
        fc.setBrickSize(glm::ivec2(64, 64));
        fc.setPadding(glm::ivec2(1, 1));
        fc.setInputRectangle(glm::ivec2(0, 0), glm::ivec2(4096, 4096));
        fc.convertFolder();
        delete compressor;
    }

    // Create uncompressed output
    {
        FitsConverter fc;
        fc.setInFolder("input/");
        fc.setOutFile("two_raw.data");

        fc.setCompressor(nullptr);
        fc.setBrickSize(glm::ivec2(64, 64));
        fc.setPadding(glm::ivec2(1, 1));
        fc.setInputRectangle(glm::ivec2(0, 0), glm::ivec2(4096, 4096));
        fc.convertFolder();
    }

}
