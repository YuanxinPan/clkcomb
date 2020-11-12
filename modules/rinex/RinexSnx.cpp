#include "RinexSnx.h"
#include "../io/io.h"

#include <stdio.h>
#include <string.h>
#include <algorithm>
#include <cctype>

bool RinexSnx::read(const std::string &path)
{
    FILE *fp = fopen(path.c_str(), "r");
    if (fp == nullptr) {
        fprintf(stderr, MSG_ERR "no such file: %s\n", path.c_str());
        return false;
    }

    char buf[BUFSIZ];
    while (fgets(buf, sizeof(buf), fp)) {
        if (strncmp(buf, "+SOLUTION/ESTIMATE", 18) == 0) {
            fgets(buf, sizeof(buf), fp);
            break;
        }
    }

    std::string site;
    CartCoor pos;
    while (fgets(buf, sizeof(buf), fp)) {
        if (strncmp(buf, "-SOLUTION/ESTIMATE", 18) == 0)
            break;
        else if (strncmp(buf+7, "STA", 3) != 0)
            continue;

        switch (buf[10]) {
            case 'X':
                pos.x = atof(buf+47);
                break;
            case 'Y':
                pos.y = atof(buf+47);
                break;
            case 'Z':
                pos.z = atof(buf+47);
                site.assign(buf+14, 4);
                std::transform(site.begin(), site.end(), site.begin(),
                            [](unsigned char c) { return std::toupper(c); });
                poss_[site] = pos;
                break;
            default:
                break;
        }
    }

    // printf("%s\n", path.c_str());
    // for (auto it=poss_.begin(); it!=poss_.end(); ++it) {
    //     printf("%4s %14.3f %14.3f %14.3f\n", it->first.c_str(), it->second.x, it->second.y, it->second.z);
    // }
    // printf("\n");

    return true;
}

bool RinexSnx::find_pos(const std::string &site, CartCoor &pos) const
{
    if (poss_.find(site) == poss_.end())
        return false;
    else {
        auto it = poss_.find(site);
        pos = it->second;
        return true;
    }
}
