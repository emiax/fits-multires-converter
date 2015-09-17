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
        unsigned int nTimesteps;
        glm::ivec2 originalSize;
        glm::ivec2 croppedSize;
        bool hasMinMaxEnergy = false;
        bool hasMinMaxExpTime = false;
        double minFlux;
        double maxFlux;
        int16_t maxEnergy;
        int16_t minEnergy;
        double minExpTime;
        double maxExpTime;
        unsigned int startTime;
        unsigned int endTime;
    };

    struct ImageMetaData {
        std::string filename;
        unsigned int timestamp;
        double expTime;
    };
 private:
    bool validateInput();
    bool convertFile(std::string filename, ImageMetaData& metaData, CommonMetaData& common, std::fstream& file);
    bool readCommonMetaData(std::string filename, CommonMetaData& metaData);
    bool readImageMetaData(std::string filename, ImageMetaData& imageMetaData, CommonMetaData& common);
    bool createHeader(const CommonMetaData& common, const std::vector<ImageMetaData>& imageMetaData, std::fstream& out, bool log = true);
    bool createHeaderPlaceholder(const CommonMetaData& common, std::fstream& out);

    std::string _inFolderName;
    std::string _outFileName;
    glm::ivec2 _brickSize;
    glm::ivec2 _padding;
    glm::ivec2 _inputRectangleTopLeft;
    glm::ivec2 _inputRectangleSize;

};


#endif // __FITSCONVERTER_H__
