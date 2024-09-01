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

#include "io.h"

#include <string>
#include <string.h>

void skip_nline(FILE *fp, int n)
{
    static char buf[BUFSIZ];
    while (n-- > 0)
        fgets(buf, sizeof(buf), fp);
}

void skip_header(FILE *fp, int shift)
{
    char buf[128];
    do {
        fgets(buf, sizeof(buf), fp);
    } while (strncmp(buf+shift, "END OF HEADER", 13));
}

//void skip_nline(std::ifstream &in, int n)
//{
//    static std::string buf;
//    while (n-- > 0)
//        std::getline(in, buf);
//}

void backspace(FILE *fp)
{
    int  num = 0;
    long pos = ftell(fp);
    char tmp;
    while (true)
    {
        --pos;
        fseek(fp, pos, SEEK_SET);
        tmp = fgetc(fp);
        if (pos == -1) break;
        if (tmp == '\n') {
            ++num;
            if (num == 2) break;
            while (tmp == '\n') {
                pos--;
                fseek(fp, pos, SEEK_SET);
                tmp = fgetc(fp);
            }
        }
    }
    fseek(fp, ++pos, SEEK_SET);
}

