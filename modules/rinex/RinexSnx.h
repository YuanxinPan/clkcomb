#ifndef RINEXSNX_H
#define RINEXSNX_H

#include "../coord/CartCoor.h"

#include <map>
#include <string>

class RinexSnx {
public:
    RinexSnx() {}
    ~RinexSnx(){}

    bool read(const std::string &path);
    bool find_pos(const std::string &site, CartCoor &pos) const;

private:
    std::map<std::string, CartCoor> poss_;
};

#endif // RINEXSNX_H
