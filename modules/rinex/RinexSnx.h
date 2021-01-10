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
    bool find_dome(const std::string &site, std::string &dome) const;

private:
    std::map<std::string, std::string> domes_;
    std::map<std::string, CartCoor> poss_;
};

#endif // RINEXSNX_H
