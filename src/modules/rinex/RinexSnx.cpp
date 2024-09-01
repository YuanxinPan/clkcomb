// clkcomb - Clock and phase bias products Combination
// Copyright (C) 2021 Yuanxin Pan
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice,
//    this list of conditions and the following disclaimer.
// 2. Redistributions in binary form must reproduce the above copyright notice,
//    this list of conditions and the following disclaimer in the documentation
//    and/or other materials provided with the distribution.
// 3. Neither the name of the copyright holder nor the names of its contributors
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

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

    // DOME
    char buf[BUFSIZ];
    while (fgets(buf, sizeof(buf), fp)) {
        if (strncmp(buf, "+SITE/ID", 8) == 0) {
            // fgets(buf, sizeof(buf), fp);
            break;
        }
    }

    std::string site;
    std::string dome;
    while (fgets(buf, sizeof(buf), fp)) {
        if (strncmp(buf, "-SITE/ID", 8) == 0)
            break;
        else if (strncmp(buf, "*CODE", 5) == 0)
            continue;

        site.assign(buf+1, 4);
        dome.assign(buf+9, 9);
        std::transform(site.begin(), site.end(), site.begin(),
                    [](unsigned char c) { return std::toupper(c); });
        domes_[site] = dome;
    }

    // ESTIMATE
    while (fgets(buf, sizeof(buf), fp)) {
        if (strncmp(buf, "+SOLUTION/ESTIMATE", 18) == 0) {
            // fgets(buf, sizeof(buf), fp);
            break;
        }
    }

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

bool RinexSnx::find_dome(const std::string &site, std::string &dome) const
{
    if (domes_.find(site) == domes_.end())
        return false;
    else {
        auto it = domes_.find(site);
        dome = it->second;
        return true;
    }
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
