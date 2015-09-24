#include <fitsconverter.h>
#include <set>
#include <iostream>
#include <CCfits>
#include <lz4/lz4.h>
#include <cassert>
#include <sstream>
#include <string>
#include <stdint.h>
#include <vector>
#include <compressor.h>

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
            out.x |= (in >> (i*2) & 1) << i;
            out.y |= (in >> (i*2 + 1) & 1) << i;
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
    , _inputRectangleSize(-1)
    , _compressor(nullptr) {}

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

void FitsConverter::setCompressor(Compressor* compressor) {
    _compressor = compressor;
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
    if (!readCommonMetaData(filenames.begin()->string(), common)) {
        std::cerr << "Could not load common meta data from " << *filenames.begin() << std::endl;
        return false;
    }

    glm::ivec2 croppedSize = common.croppedSize;
    unsigned int bricksPerDim = croppedSize.x / _brickSize.x;
    unsigned int nLevels = log2(bricksPerDim) + 1;
    unsigned int nBricksInQuadTree = (pow(4, nLevels) - 1) / 3;

    common.compressionType = _compressor ? static_cast<unsigned int>(_compressor->compressionType()) : 0;
    common.nTimesteps = nTimesteps;
    common.nBricks = common.nTimesteps * nBricksInQuadTree;

    std::vector<ImageMetaData> imageMetaData;
    std::vector<BrickMetaData> brickMetaData;

    imageMetaData.resize(common.nTimesteps);
    brickMetaData.resize(common.nBricks);

    int i = 0;
    for (const fs::path& filename : filenames) {
        ImageMetaData metaData;
        if (!readImageMetaData(filename.string(), metaData, common)) {
            std::cerr << "Failed to load meta data from " << filename << std::endl;
            return false;
        }
        imageMetaData[i] = metaData;
        i++;
    }


    std::fstream out(_outFileName, std::fstream::out | std::fstream::binary);
    createHeaderPlaceholder(common, out);

    i = 0;
    for (const fs::path& filename : filenames) {
        std::cout << "\rConverting timestep " << (i+1) << "/" << nTimesteps << " from " << filename << "..." << std::flush;
        if (convertFile(filename.string(), common, imageMetaData[i], brickMetaData, out)) {
            i++;
        } else {
            std::cerr << "Failed to convert " << filename << std::endl;
            out.close();
            return false;
        }
    }

    std::cout << "\r" << std::endl;

    if (nTimesteps == 1) {
        std::cout << "Successfully converted 1 timestep from \""
                  << _inFolderName << "\" to \"" << _outFileName << "\"" << std::endl;
    } else {
        std::cout << "Successfully converted " << nTimesteps << " timesteps from "
                  << _inFolderName << "\" to \"" << _outFileName << "\"" << std::endl;
    }
    createHeader(common, imageMetaData, brickMetaData, out);
    out.close();
    return true;
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


bool FitsConverter::readImageMetaData(std::string filename, ImageMetaData& imageMetaData, CommonMetaData& common) {
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

    double expTime;
    image.readKey("EXPTIME", expTime);
    imageMetaData.expTime = expTime;

    if (common.hasMinMaxExpTime) {
        common.minExpTime = std::min(common.minExpTime, expTime);
        common.maxExpTime = std::max(common.maxExpTime, expTime);
    } else {
        common.minExpTime = common.maxExpTime = expTime;
        common.hasMinMaxExpTime = true;
    }

    return true;
}

bool FitsConverter::createHeaderPlaceholder(const CommonMetaData& common, std::fstream& out) {
    int size = 68 + common.nTimesteps*sizeof(double) + common.nBricks*sizeof(unsigned int);
    char* zeros = new char[size];
    out.write(zeros, size);
    return true;
}


bool FitsConverter::createHeader(const CommonMetaData& common, const std::vector<ImageMetaData>& imageMetaData, const std::vector<BrickMetaData>& brickMetaData, std::fstream& out, bool log) {
    out.seekg(0);
    uint32_t nTimesteps = common.nTimesteps;
    uint32_t nBricks = common.nBricks;
    uint32_t compressionType = common.compressionType;
    uint32_t nBricksPerDim = common.croppedSize.x / _brickSize.x;
    uint32_t brickWidth = _brickSize.x;
    uint32_t brickHeight = _brickSize.y;
    uint32_t paddingX = _padding.x;
    uint32_t paddingY = _padding.y;
    int16_t minEnergy = common.minEnergy;
    int16_t maxEnergy = common.maxEnergy;
    double minFlux = common.minFlux;
    double maxFlux = common.maxFlux;
    double minExpTime = common.minExpTime;
    double maxExpTime = common.maxExpTime;

    out.write(reinterpret_cast<char*>(&nTimesteps), sizeof(uint32_t));
    out.write(reinterpret_cast<char*>(&nBricks), sizeof(uint32_t));
    out.write(reinterpret_cast<char*>(&compressionType), sizeof(uint32_t));
    out.write(reinterpret_cast<char*>(&nBricksPerDim), sizeof(uint32_t));
    out.write(reinterpret_cast<char*>(&brickWidth), sizeof(uint32_t));
    out.write(reinterpret_cast<char*>(&brickHeight), sizeof(uint32_t));
    out.write(reinterpret_cast<char*>(&paddingX), sizeof(uint32_t));
    out.write(reinterpret_cast<char*>(&paddingY), sizeof(uint32_t));
    out.write(reinterpret_cast<char*>(&minEnergy), sizeof(int16_t));
    out.write(reinterpret_cast<char*>(&maxEnergy), sizeof(int16_t));
    out.write(reinterpret_cast<char*>(&minFlux), sizeof(double));
    out.write(reinterpret_cast<char*>(&maxFlux), sizeof(double));
    out.write(reinterpret_cast<char*>(&minExpTime), sizeof(double));
    out.write(reinterpret_cast<char*>(&maxExpTime), sizeof(double));

    int col1 = 30;
    int col2 = 0;
    if (log) {
        std::cout << std::setw(col1) << "Timesteps: " << std::setw(col2) << nTimesteps << std::endl;
        std::cout << std::setw(col1) << "Bricks per dimension: " << std::setw(col2) << nBricksPerDim << std::endl;
        std::cout << std::setw(col1) << "Brick dimensions: " << std::setw(col2) << brickWidth << "x" << brickHeight << std::endl;
        std::cout << std::setw(col1) << "Padding: " << std::setw(col2) << paddingX << "x" << paddingY << " (on each side)"<< std::endl;
        std::cout << std::setw(col1) << "Min energy: " << std::setw(col2) << minEnergy << std::endl;
        std::cout << std::setw(col1) << "Max energy: " << std::setw(col2) << maxEnergy << std::endl;
        std::cout << std::setw(col1) << "Min flux: " << std::setw(col2) << minFlux << std::endl;
        std::cout << std::setw(col1) << "Max flux: " << std::setw(col2) << maxFlux << std::endl;
        std::cout << std::setw(col1) << "Min exposure time: " << std::setw(col2) << minExpTime << std::endl;
        std::cout << std::setw(col1) << "Max exposure time: " << std::setw(col2) << maxExpTime << std::endl;
    }

    for (unsigned int i = 0; i < imageMetaData.size(); i++) {
        const ImageMetaData& data = imageMetaData[i];
        double expTime = data.expTime;
        out.write(reinterpret_cast<char*>(&expTime), sizeof(double));
    }

    for (unsigned int i = 0; i < brickMetaData.size(); i++) {
        const BrickMetaData& data = brickMetaData[i];
        unsigned int dataPosition = data.dataPosition;
        out.write(reinterpret_cast<char*>(&dataPosition), sizeof(unsigned int));
    }


    return true;
}


bool FitsConverter::convertFile(std::string inFilename, CommonMetaData& common, ImageMetaData& metaData, std::vector<BrickMetaData>& brickMetaData, std::fstream& out) {
    CCfits::FITS file(inFilename, CCfits::Read, true);
    CCfits::PHDU& image = file.pHDU();

    glm::ivec2 originalSize = common.originalSize;
    glm::ivec2 croppedSize = common.croppedSize;

    std::valarray<int16_t> contents;
    image.readAllKeys();
    image.read(contents);

    unsigned int bricksPerDim = croppedSize.x / _brickSize.x;
    unsigned int nLevels = log2(bricksPerDim) + 1;

    std::vector<std::vector<double> > levels(nLevels);

    // Read from contents and crop data if necessary.
    std::vector<double>& leafLevel = levels[0];
    leafLevel.resize(croppedSize.x * croppedSize.y);
    for (int y = 0; y < croppedSize.y; y++) {
        for (int x = 0; x < croppedSize.x; x++) {

            int16_t low = 0;
            int16_t energy = glm::max(contents[cartesianToLinear(glm::ivec2(x, y) + _inputRectangleTopLeft, originalSize)], low);

            double flux = static_cast<double>(energy)/static_cast<double>(metaData.expTime);
            if (common.hasMinMaxEnergy) {
                common.minEnergy = std::min(common.minEnergy, energy);
                common.maxEnergy = std::max(common.maxEnergy, energy);
                common.minFlux = std::min(common.minFlux, flux);
                common.maxFlux = std::max(common.maxFlux, flux);
            } else {
                common.minEnergy = common.maxEnergy = energy;
                common.minFlux = common.maxFlux = flux;
                common.hasMinMaxEnergy = true;
            }
            leafLevel[cartesianToLinear(glm::ivec2(x, y), croppedSize)] = energy;
        }
    }




    // Downsample data for all levels
    for (unsigned int levelIndex = 1; levelIndex < nLevels; levelIndex++) {
        std::vector<double>& level = levels[levelIndex];
        std::vector<double>& childLevel = levels[levelIndex - 1];
        glm::ivec2 levelSize = croppedSize / static_cast<int>(pow(2, levelIndex));
        glm::ivec2 childLevelSize = croppedSize / static_cast<int>(pow(2, levelIndex - 1));

        level.resize(levelSize.x * levelSize.y);
        for (int y = 0; y < levelSize.y; y++) {
            for (int x = 0; x < levelSize.x; x++) {
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

    unsigned int nBrickVals = paddedBrickSize.x * paddedBrickSize.y;

    std::vector<int16_t> brickData(nBrickVals);

    for (int levelIndex = nLevels - 1; levelIndex >= 0; levelIndex--) {
        std::vector<double>& level = levels[levelIndex];
        unsigned int bricksInLevel = pow(4, nLevels - levelIndex - 1);

        glm::ivec2 levelSize = croppedSize / static_cast<int>(pow(2, levelIndex));
        for (unsigned int zIndex = 0; zIndex < bricksInLevel; zIndex++) {
            glm::ivec2 brickCoord = zIndexToCartesian(zIndex);

            glm::ivec2 brickMin = brickCoord * _brickSize;

            for (int brickY = 0; brickY < paddedBrickSize.y; brickY++) {
                for (int brickX = 0; brickX < paddedBrickSize.x; brickX++) {
                    int globalY = brickMin.y - padding.y + brickY;
                    int globalX = brickMin.x - padding.x + brickX;
                    unsigned int pixelOffset = cartesianToLinear(glm::ivec2(brickX, brickY), paddedBrickSize);
                    unsigned int inputOffset = clampCartesianToLinear(glm::ivec2(globalX, globalY), levelSize);
                    brickData[pixelOffset] = glm::clamp(static_cast<int>(std::round(level[inputOffset])), 0, UINT16_MAX);
                }
            }

            char* output = nullptr;
            size_t outputSize = 0;
            if (_compressor) {
                // _compressor.compress allocates a vector of the required size
                output = _compressor->compress(brickData.data(), paddedBrickSize, outputSize);
            } else {
                outputSize = nBrickVals * sizeof(int16_t);
                output = new char[outputSize];
                memcpy(output, reinterpret_cast<char*>(brickData.data()), outputSize);
            }

            brickMetaData[common.currentBrickIndex].dataPosition = out.tellg();
            out.write(output, outputSize);

            delete[] output;
            common.currentBrickIndex++;
        }
    }
    return true;
}

