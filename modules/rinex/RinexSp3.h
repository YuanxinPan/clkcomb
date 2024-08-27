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

#ifndef RINEXSP3_H
#define RINEXSP3_H

#include "../coord/coord.h"
#include "../chrono/chrono.h"

#include <vector>
#include <string>

class RinexSp3
{
public:
    struct Sp3_t {
        MJD t;
        CartCoor pos;
    };

public:
    RinexSp3(){}
    ~RinexSp3(){}

    //static std::vector<std::string> satPRN(const std::string &path);
    bool read(const std::string &path);
    bool read(const std::vector<std::string> &paths);
    bool satPos(MJD t, const std::string &prn, double *pos)const;
    bool satPos(MJD t, const std::string &prn, CartCoor &pos)const;
    bool satPosVel(MJD t, const std::string &prn, double *pos, double *vel)const;
    bool satPosVel(MJD t, const std::string &prn, CartCoor &pos, CartCoor &vel)const;

    const std::vector<std::string> &sat_prns() { return prns_; }

private:
    std::vector<std::string> prns_;
    std::vector<std::vector<Sp3_t>> satposs_;
};

bool operator<(const RinexSp3::Sp3_t &lhs, const RinexSp3::Sp3_t &rhs);

#endif
