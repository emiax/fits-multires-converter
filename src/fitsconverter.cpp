#include <fitsconverter.h>
#include <set>
#include <iostream>
#include <CCfits>

namespace {
    int cartesianToLinear(glm::ivec2 in, glm::ivec2 size) {
        return in.x + in.y * size.x;
    }

    int clampCartesianToLinear(glm::ivec2 in, glm::ivec2 size) {
        glm::ivec2 clamped = glm::clamp(in, glm::ivec2(0, 0), size - glm::ivec2(1, 1));
        return clamped.x + clamped.y * size.x;
    }

    glm::ivec2 linearToCartesian(int in, glm::ivec2 size) {
        return glm::ivec2(in % size.x, in / size.x);
    }

    bool isPowerOfTwo(unsigned int in) {
        return (in & (in - 1)) == 0;
    }

    glm::ivec2 zIndexToCartesian(uint32_t in) {
        glm::ivec2 out(0);
        for (int i = 0; i < 16; i++) {
            out.x += (i+1) * (in & (i << (i*2)));
            out.y += (i+1) * (in & (i << (i*2 + 1)));
        }
        return out;
    }
}

FitsConverter::FitsConverter()
    : _inFolderName("")
    , _outFileName("")
    , _brickSize(0)
    , _padding(0)
    , _inputRectangleTopLeft(0)
    , _inputRectangleSize(-1) {}

void FitsConverter::setInFolder(const std::string& inFolder) {
    _inFolderName = inFolder;
}

void FitsConverter::setOutFile(const std::string& outFile) {
    _outFileName = outFile;
}

void FitsConverter::setBrickSize(const glm::ivec2& brickSize) {
    _brickSize = brickSize;
}

void FitsConverter::setPadding(const glm::ivec2& padding) {
    _padding = padding;
}

void FitsConverter::setInputRectangle(const glm::ivec2& topLeft, const glm::ivec2& size) {
    _inputRectangleTopLeft = topLeft;
    _inputRectangleSize = size;
}
bool FitsConverter::convertFolder() {
    if (!validateInput()) {
        return false;
    }

    std::set<fs::path> filenames;
    for (fs::directory_iterator it(_inFolderName); it!=fs::directory_iterator(); ++it) {
        filenames.insert(it->path());
    }

    unsigned int nTimesteps;
    nTimesteps = filenames.size();

    if (nTimesteps == 0) {
        std::cerr << "No input files found in directory" << std::endl;
        return false;
    } else if (nTimesteps == 1) {
        std::cout << "Found 1 timestep" << std::endl;
    } else {
        std::cout << "Found " << nTimesteps << " timesteps" << std::endl;
    }

    CommonMetaData common;
    std::vector<ImageMetaData> imageMetaData;
    imageMetaData.reserve(nTimesteps);

    if (!readCommonMetaData(filenames.begin()->string(), common)) {
        std::cerr << "Could not load common meta data from " << *filenames.begin() << std::endl;
        return false;
    }

    for (const fs::path& filename : filenames) {
        ImageMetaData metaData;
        if (!readImageMetaData(filename.string(), metaData, common)) {
            std::cerr << "Failed to load meta data from " << filename << std::endl;
            return false;
        }
        imageMetaData.push_back(metaData);
    }

    std::fstream out(_outFileName, std::fstream::out | std::fstream::binary);
    createHeader(common, imageMetaData, out);
    int i = 0;
    for (const fs::path& filename : filenames) {
        if (convertFile(filename.string(), common, out)) {
            i++;
            std::cout << "Converted timestep " << i << "/" << nTimesteps << " from " << filename << std::endl;
        } else {
            std::cerr << "Failed to convert " << filename << std::endl;
            out.close();
            return false;
        }
    }

    if (nTimesteps == 1) {
        std::cout << "Successfully converted 1 timestep from \""
                  << _inFolderName << "\" to \"" << _outFileName << "\"" << std::endl;
    } else {
        std::cout << "Successfully converted " << nTimesteps << " timesteps from "
                  << _inFolderName << "\" to \"" << _outFileName << "\"" << std::endl;
    }
    out.close();
}


bool FitsConverter::validateInput() {
  if (!fs::exists(_inFolderName)) {
    std::cerr << "The input folder path \"" << _inFolderName << "\" does not exist" << std::endl;
    return false;
  }

  if (!fs::is_directory(_inFolderName)) {
    std::cerr << "The path "<< _inFolderName << " is not a directory" << std::endl;
    return false;
  }

  if (_brickSize.x <= 0 || _brickSize.y <= 0) {
    std::cerr << "The brick size must be a vector of two numbers larger than zero" << std::endl;
    return false;
  }

  if (_padding.x < 0 || _padding.y < 0) {
    std::cerr << "The brick size must be a vector of two non-negative numbers" << std::endl;
    return false;
  }

  return true;
}


bool FitsConverter::readCommonMetaData(std::string filename, CommonMetaData& metaData) {
    CCfits::FITS file(filename, CCfits::Read, true);
    CCfits::PHDU& image = file.pHDU();
    int nAxes = image.axes();
    if (nAxes != 2) {
        std::cerr << "Fits file does not have 2 axes. It has " << nAxes << std::endl;
        return false;
    }
    int width = image.axis(0);
    int height = image.axis(1);

    metaData.originalSize = glm::ivec2(width, height);

    if (_inputRectangleSize.x == -1 && _inputRectangleSize.y == -1) {
        metaData.croppedSize = metaData.originalSize - _inputRectangleTopLeft;
    } else {
        metaData.croppedSize = _inputRectangleSize;
    }

    if (metaData.croppedSize.x > metaData.originalSize.x ||
        metaData.croppedSize.y > metaData.originalSize.y) {
        std::cerr << "Orginal image does not fit entire input rectangle" << std::endl;
        return false;
    }

    if (!(metaData.croppedSize % _brickSize == glm::ivec2(0, 0))) {
        std::cerr << "Input rectangle dimensions (" << metaData.croppedSize.x << ", " << metaData.croppedSize.y << ") is not a multiple of the brick dimenstions  (" << _brickSize.x << ", " << _brickSize.y << ")" << std::endl;
        return false;
    }

    glm::ivec2 bricksPerDim = metaData.croppedSize / _brickSize;
    if (bricksPerDim.x != bricksPerDim.y) {
        std::cerr << "Number of bricks per dimension is not equal" << std::endl;
        return false;
    }

    if (!isPowerOfTwo(bricksPerDim.x)) {
        std::cerr << "Number of bricks per dimension (" << bricksPerDim.x << ") is not a power of two. " << std::endl;
        return false;
    }

    std::cout << "Number of bricks per dimenssion: " << bricksPerDim.x << std::endl;
    glm::ivec2 paddedBrickSize = _brickSize + _padding*2;
    std::cout << "Padded brick size: (" << paddedBrickSize.x << ", " << paddedBrickSize.y << ")" << std::endl;


    return true;
}


bool FitsConverter::readImageMetaData(std::string filename, ImageMetaData& imageMetaData, const CommonMetaData& common) {
    CCfits::FITS file(filename, CCfits::Read, true);
    CCfits::PHDU& image = file.pHDU();
    int nAxes = image.axes();

    if (nAxes != 2) {
        std::cerr << "Fits file does not have 2 axes. It has " << nAxes << std::endl;
        return false;
    }

    int width = image.axis(0);
    int height = image.axis(1);
    glm::ivec2 size(width, height);

    if (size != common.originalSize) {
        std::cerr << "Image size of " << filename << "(" << width << ", " << height << " is not same as first input file (" << common.originalSize.x << ", " << common.originalSize.y << ")." << std::endl;
        return false;
    }

    imageMetaData.filename = filename;
    imageMetaData.timestamp = 0; // todo: actually extract this from file.
    return true;
}


bool FitsConverter::createHeader(const CommonMetaData& common, const std::vector<ImageMetaData>& imageMetaData, std::fstream& out) {

    int32_t nTimesteps = imageMetaData.size();
    int32_t nBricksPerDim = common.croppedSize.x / _brickSize.x;
    int32_t brickWidth = _brickSize.x;
    int32_t brickHeight = _brickSize.y;
    int32_t paddingX = _padding.x;
    int32_t paddingY = _padding.y;

    out.write(reinterpret_cast<char*>(&nTimesteps), sizeof(int32_t));
    out.write(reinterpret_cast<char*>(&nBricksPerDim), sizeof(int32_t));
    out.write(reinterpret_cast<char*>(&brickWidth), sizeof(int32_t));
    out.write(reinterpret_cast<char*>(&brickHeight), sizeof(int32_t));
    out.write(reinterpret_cast<char*>(&paddingX), sizeof(int32_t));
    out.write(reinterpret_cast<char*>(&paddingY), sizeof(int32_t));

    // todo: output some meta data for each image.
}


bool FitsConverter::convertFile(std::string inFilename, const CommonMetaData& common, std::fstream& out) {

    CCfits::FITS file(inFilename, CCfits::Read, true);
    CCfits::PHDU& image = file.pHDU();

    glm::ivec2 originalSize = common.originalSize;
    glm::ivec2 croppedSize = common.croppedSize;

    std::valarray<int16_t> contents;
    image.readAllKeys();
    image.read(contents);

    unsigned int bricksPerDim = croppedSize.x / _brickSize.x;
    unsigned int nLevels = log2(bricksPerDim) + 1;


    std::vector<std::vector<int16_t> > levels(nLevels);

    // Read from contents and crop data if necessary.
    std::vector<int16_t>& leafLevel = levels[0];
    leafLevel.resize(croppedSize.x * croppedSize.y);
    for (unsigned int y = 0; y < croppedSize.y; y++) {
        for (unsigned int x = 0; x < croppedSize.x; x++) {
            leafLevel[cartesianToLinear(glm::ivec2(x, y), croppedSize)] =
                contents[cartesianToLinear(glm::ivec2(x, y) + _inputRectangleTopLeft, originalSize)];
        }
    }

    // Downsample data for all levels
    for (unsigned int levelIndex = 1; levelIndex < nLevels; levelIndex++) {
        std::vector<int16_t>& level = levels[levelIndex];
        std::vector<int16_t>& childLevel = levels[levelIndex - 1];
        glm::ivec2 levelSize = croppedSize / static_cast<int>(pow(2, levelIndex));
        glm::ivec2 childLevelSize = croppedSize / static_cast<int>(pow(2, levelIndex - 1));

        level.resize(levelSize.x * levelSize.y);
        for (unsigned int y = 0; y < levelSize.y; y++) {
            for (unsigned int x = 0; x < levelSize.x; x++) {
                level[cartesianToLinear(glm::ivec2(x, y), levelSize)] =
                    (childLevel[cartesianToLinear(glm::ivec2(x*2, y*2), childLevelSize)] +
                     childLevel[cartesianToLinear(glm::ivec2(x*2 + 1, y*2), childLevelSize)] +
                     childLevel[cartesianToLinear(glm::ivec2(x*2, y*2 + 1), childLevelSize)] +
                     childLevel[cartesianToLinear(glm::ivec2(x*2 + 1, y*2 + 1), childLevelSize)]) * 0.25;
            }
        }
    }

    // Write quad tree top-down, in z-order.
    glm::ivec2 padding = _padding;
    glm::ivec2 paddedBrickSize = _brickSize + 2*padding;
    unsigned int nBricksInQuadTree = (pow(4, nLevels) - 1) / 3;

    unsigned int nBrickVals = paddedBrickSize.x * paddedBrickSize.y;

    std::vector<int16_t> brickData(nBricksInQuadTree * nBrickVals);
    for (int levelIndex = nLevels - 1; levelIndex >= 0; levelIndex--) {
        std::vector<int16_t>& level = levels[levelIndex];
        unsigned int bricksInLevel = pow(4, nLevels - levelIndex - 1);

        glm::ivec2 levelSize = croppedSize / static_cast<int>(pow(2, levelIndex));
        unsigned int levelOffset = (pow(4, nLevels - levelIndex - 1) - 1) / 3;

        for (int zIndex = 0; zIndex < bricksInLevel; zIndex++) {
            glm::ivec2 brickCoord = zIndexToCartesian(zIndex);
            glm::ivec2 brickMin = brickCoord * _brickSize;
            glm::ivec2 brickMax = brickCoord * (_brickSize + glm::ivec2(1));

            for (int brickY = 0; brickY < paddedBrickSize.y; brickY++) {
                for (int brickX = 0; brickX < paddedBrickSize.x; brickX++) {
                    int globalY = brickMin.y - padding.y + brickY;
                    int globalX = brickMin.x - padding.x + brickX;
                    unsigned int brickOffset = (levelOffset + zIndex) * nBrickVals;
                    unsigned int voxelOffset = cartesianToLinear(glm::ivec2(brickX, brickY), paddedBrickSize);
                    unsigned int inputOffset = clampCartesianToLinear(glm::ivec2(globalX, globalY), levelSize);
                    brickData[brickOffset + voxelOffset] = level[inputOffset];
                }
            }
        }
    }

    out.write(reinterpret_cast<char*>(&brickData[0]), sizeof(int16_t) * brickData.size());
    return true;
}

