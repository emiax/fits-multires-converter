#include <CCfits>
#include <iostream>
#include <fitsconverter.h>

int main() {
    FitsConverter fc;
    fc.setInFolder("input");
    fc.setOutFile("out.data");

    fc.setBrickSize(glm::ivec2(510, 510));
    fc.setPadding(glm::ivec2(1, 1));
    fc.setInputRectangle(glm::ivec2(8, 8), glm::ivec2(4080, 4080));
    fc.convertFolder();

    return 0;
}
