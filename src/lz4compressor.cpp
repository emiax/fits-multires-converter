#include <lz4compressor.h>

#include <lz4/lz4.h>
#include <vector>
#include <cstring>

namespace {
    int cartesianToLinear(glm::ivec2 in, glm::ivec2 size) {
        return in.x + in.y * size.x;
    }
}

Lz4Compressor::Lz4Compressor() {}

Lz4Compressor::~Lz4Compressor() {}

CompressionType Lz4Compressor::compressionType() {
    return CompressionType::LZ4_2D_PRED;
}

char* Lz4Compressor::compress(int16_t* brickData, glm::ivec2 inSize, size_t &outSize) {
    std::vector<int16_t> predErrorData(inSize.x * inSize.y);
    for (int y = 0; y < inSize.y; y++) {
        for (int x = 0; x < inSize.x; x++) {
            unsigned int here = cartesianToLinear(glm::ivec2(x, y), inSize);
            unsigned int left = cartesianToLinear(glm::ivec2(x - 1, y), inSize);
            unsigned int top = cartesianToLinear(glm::ivec2(x, y - 1), inSize);
            unsigned int topLeft = cartesianToLinear(glm::ivec2(x - 1, y - 1), inSize);
            int16_t prediction = 0;
            if (x > 0) {
                if (y > 0) {
                    prediction = brickData[left] + brickData[top] - brickData[topLeft];
                } else {
                    prediction = brickData[left];
                }
            } else if (y > 0) {
                prediction = brickData[top];
            }
            predErrorData[here] = brickData[here] - prediction;
        }
    }

    int predErrorDataSize = sizeof(int16_t) * predErrorData.size();
    int maxCompressionSize = LZ4_compressBound(predErrorDataSize);
    
    char* compressed = new char[maxCompressionSize];
    outSize = LZ4_compress(reinterpret_cast<char*>(predErrorData.data()), compressed, predErrorDataSize);
    
    char* outData = new char[outSize];
    memcpy(outData, compressed, outSize);

    delete compressed;
    return outData;
}


