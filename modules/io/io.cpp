#include "io.h"

#include <string>
#include <string.h>

void skip_nline(FILE *fp, int n)
{
    static char buf[BUFSIZ];
    while (n-- > 0)
        fgets(buf, sizeof(buf), fp);
}

void skip_header(FILE *fp)
{
    char buf[128];
    do {
        fgets(buf, sizeof(buf), fp);
    } while (strncmp(buf+60, "END OF HEADER", 13));
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

