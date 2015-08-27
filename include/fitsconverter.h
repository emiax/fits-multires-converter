#ifndef __FITSCONVERTER_H__
#define __FITSCONVERTER_H__

#include <boost/filesystem.hpp>
#include <fstream>
#include <glm/glm.hpp>

namespace fs = boost::filesystem;

class FitsConverter {
 public:
    FitsConverter();
    void setInFolder(const std::string& inFolder);
    void setOutFile(const std::string& outFile);
    void setBrickSize(const glm::ivec2& brickSize);
    void setPadding(const glm::ivec2& padding);
    void setInputRectangle(const glm::ivec2& topLeft, const glm::ivec2& size);
    bool convertFolder();

    struct CommonMetaData {
        glm::ivec2 originalSize;
        glm::ivec2 croppedSize;
    };

    struct ImageMetaData {
        std::string filename;
        unsigned int timestamp;
    };
 private:
    bool validateInput();
    bool convertFile(std::string filename, const CommonMetaData& common, std::fstream& file);
    bool readCommonMetaData(std::string filename, CommonMetaData& metaData);
    bool readImageMetaData(std::string filename, ImageMetaData& imageMetaData, const CommonMetaData& common);
    bool createHeader(const CommonMetaData& common, const std::vector<ImageMetaData>& imageMetaData, std::fstream& out);

    std::string _inFolderName;
    std::string _outFileName;
    glm::ivec2 _brickSize;
    glm::ivec2 _padding;
    glm::ivec2 _inputRectangleTopLeft;
    glm::ivec2 _inputRectangleSize;

};


#endif // __FITSCONVERTER_H__
