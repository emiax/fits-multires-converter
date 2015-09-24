#ifndef __LZ4COMPRESSOR_H__
#define __LZ4COMPRESSOR_H__

#include <compressor.h>

class Lz4Compressor : public Compressor {
 public:
    Lz4Compressor();
    virtual ~Lz4Compressor();
    virtual CompressionType compressionType();
    virtual char* compress(int16_t* brickData, glm::ivec2 inSize, size_t &outSize);
};

#endif //  __LZ4COMPRESSOR_H__

