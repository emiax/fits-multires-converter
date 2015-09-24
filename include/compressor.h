#ifndef __COMPRESSOR_H__
#define __COMPRESSOR_H__

#include <glm/glm.hpp>

enum class CompressionType: unsigned int {
    NONE = 0,
    LZ4_2D_PRED = 1
};

class Compressor {
 public:
    Compressor() {};
    virtual ~Compressor() {};
    virtual CompressionType compressionType() = 0;
    virtual char* compress(int16_t* brickData, glm::ivec2 inSize, size_t &outSize) = 0;
};


#endif // __COMPRESSOR_H__
