#include "RinexAtx.h"
#include "../io/io.h"
#include "../coord/coord.h"
#include "../const.h"

#include <string.h>
#include <stdlib.h>
#include <algorithm>
#include <iterator>

bool RinexAtx::atx_t::pco(const std::string &f, double *pco)const
{
    auto it = std::find(freqs.begin(), freqs.end(), f);
    int i = it - freqs.begin();
    if (it == freqs.end())
        return false;

    memcpy(pco, pcos[i].data(), 3*sizeof(double));
    // pco[0] = pcos[i][0];
    // pco[1] = pcos[i][1];
    // pco[2] = pcos[i][2];
    return true;
}

void RinexAtx::atx_t::pco(double *pco, const std::string &f1, const std::string &f2)const
{
    double f[2] = { 0, 0 };
    // if (name.size() == 20u) {
    //     f[0] = GPS_f1; f[1] = GPS_f2;
    // } else if (name[0] == 'E') {
    //     f[0] = GAL_f1; f[1] = GAL_f5;
    // } else if (name[0] == 'R') {
    //     f[0] = GLS_f1; f[1] = GLS_f2;
    // } else if (name[0] == 'C') {
    //     f[0] = BDS_f1; f[1] = BDS_f2;
    // } else { // G J
    //     f[0] = GPS_f1; f[1] = GPS_f2;
    // }
    switch (f1[0])
    {
    case 'E':
        f[0] = GAL_f1; f[1] = GAL_f5;
        break;
    case 'R':
        f[0] = GLS_f1; f[1] = GLS_f2;
        break;
    case 'C':
        f[0] = BDS_f1; f[1] = BDS_f2;
        break;
    case 'G':
    case 'J':
    default:
        f[0] = GPS_f1; f[1] = GPS_f2;
        break;
    }

    double tmp = f[0]*f[0] - f[1]*f[1];
    double coef[2] = { f[0]*f[0]/tmp, -f[1]*f[1]/tmp };

    int i=0, j=0;
    auto it = std::find(freqs.begin(), freqs.end(), f1);
    i = it - freqs.begin();
    it = std::find(freqs.begin(), freqs.end(), f2);
    j = it - freqs.begin();

    // fprintf(stderr, "%3s %3s:%2d %3s:%2d\n", name.c_str(), f1.c_str(), i, f2.c_str(), j);

    pco[0] = coef[0]*pcos[i][0] + coef[1]*pcos[j][0];
    pco[1] = coef[0]*pcos[i][1] + coef[1]*pcos[j][1];
    pco[2] = coef[0]*pcos[i][2] + coef[1]*pcos[j][2];
}

double RinexAtx::atx_t::pcv(double zen, double azi)const
{
    zen *= R2D;
    azi *= R2D;
    if (zen > zen2)
        zen = zen2;
    else if (zen < zen1)
        zen = zen1;
    if (azi < 0.0)
        azi += 360.0;

    int irow = dazi==0.0 ? 0 : static_cast<int>(azi/dazi)+1;
    int icol = static_cast<int>((zen - zen1)/dzen);

    double pcv[2];
    for (int i=0; i<2; ++i)
    {
        double v[2];
        double ratio;
        if (irow == 0) {
            v[0] = pcvs[i][irow][icol];
            v[1] = pcvs[i][irow][icol+1];
        }
        else {
            double v1 = pcvs[i][irow  ][icol];
            double v2 = pcvs[i][irow+1][icol];
            ratio = (azi - dazi*irow)/dazi;
            v[0] = v1 + ratio*(v2 - v1);

            v1 = pcvs[i][irow  ][icol+1];
            v2 = pcvs[i][irow+1][icol+1];
            ratio = (azi - dazi*irow)/dazi;
            v[1] = v1 + ratio*(v2 - v1);
        }
        ratio = (zen - zen1 - icol*dzen)/dzen;
        pcv[i] = v[0] + ratio*(v[1] - v[0]);
    }
    //printf("zen:%8.3f\n", zen);
    //printf("val:%8.3f%8.3f\n", 1000*v[0], 1000*v[1]);
    //printf("ratio:%8.3f\n", ratio);
    return IF_coef1*pcv[0] + IF_coef2*pcv[1];
    //return v[0] + ratio*(v[1] - v[0]);
}

void RinexAtx::close()
{
    if (atxFile_ != nullptr) {
        fclose(atxFile_);
        atxFile_ = nullptr;
    }
}

bool RinexAtx::open(const std::string &path)
{
    // in case of redundant open
    if (atxFile_ != nullptr) {
        fclose(atxFile_);
        atxFile_ = nullptr;
    }
    atxFile_ = fopen(path.c_str(), "r");
    if (atxFile_ == nullptr) {
        fprintf(stderr, MSG_ERR "RinexAtx::read: no such file: %s\n", path.c_str());
        return false;
    }

    // check whether absolute atx file
    char buf[128];
    fgets(buf, sizeof(buf), atxFile_);
    fgets(buf, sizeof(buf), atxFile_);
    if (buf[0] != 'A') {
        fprintf(stderr, MSG_ERR "RinexAtx::read: not an absolute atx: %s\n", path.c_str());
        return false;
    }
    atxs_.reserve(150);
    return true;
}

const RinexAtx::atx_t *RinexAtx::atx(MJD t, const std::string &ant)const
{
    // check cache
    auto it = std::find(atxs_.begin(), atxs_.end(), ant);
    if (it != atxs_.end())
        return &*it;

    atx_t atx;
    if (find_atx(t, ant, atx)) {
        atxs_.push_back(atx);
        return &atxs_.back();
    }
    // receiver antenna without dome
    else if (ant.size() == 20u && ant.substr(16, 4)!="NONE") {
        std::string name = ant;
        name.replace(16, 4, "NONE");
        fprintf(stdout, ANSI_BOLD_BLUE ":: " ANSI_RESET
                "RinexAtx::atx:: try %s\n", name.c_str());
        if (find_atx(t, name, atx)) {
            atxs_.push_back(atx);
            return &atxs_.back();
        }
    }
    fprintf(stdout, MSG_WAR "RinexAtx::atx:: no %s\n", ant.c_str());
    return nullptr;
}

bool RinexAtx::find_atx(MJD t, const std::string &ant, atx_t &atx)const
{
    rewind(atxFile_);
    skip_header(atxFile_);

    bool issatatx = ant.size()==3u;
    char buf[128];
    MJD beg, end;
    CalendTime calend;
    while (fgets(buf, sizeof(buf), atxFile_))
    {
        if (strncmp(buf+60, "START OF ANTENNA", 16) == 0) {
            beg = end = MJD();
        } else if (strncmp(buf+60, "TYPE / SERIAL NO", 16) == 0) {
            if (issatatx && buf[20]==' ')
                break;
            else if (!issatatx && buf[20]!=' ')
                skip_atx(atxFile_);
            else {
                if (issatatx) atx.name.assign(buf+20, 3);
                else          atx.name.assign(buf, 20);
            }
            // check whether required
            if (atx.name != ant)
                skip_atx(atxFile_);
            else if (issatatx) {
                atx.svn_.assign(buf+40, 4);
                atx.blk_.assign(buf, 12);
            }
        } else if (strncmp(buf+60, "DAZI", 4) == 0) {
            sscanf(buf, "%lf", &atx.dazi);
        } else if (strncmp(buf+60, "ZEN1 / ZEN2 / DZEN", 18) == 0) {
            sscanf(buf, "%lf%lf%lf", &atx.zen1, &atx.zen2, &atx.dzen);
        } else if (strncmp(buf+60, "# OF FREQUENCIES", 16) == 0) {
            sscanf(buf, "%d", &atx.nfreq);
        } else if (strncmp(buf+60, "VALID FROM", 10) == 0) {
            calend.read(buf);
            beg = calend;
            if (beg > t)
                skip_atx(atxFile_);
        } else if (strncmp(buf+60, "VALID UNTIL", 11) == 0) {
            calend.read(buf);
            end = calend;
            if (end < t)
                skip_atx(atxFile_);
        } else if (strncmp(buf+60, "SINEX CODE", 10) == 0) {
            read_atx(atxFile_, atx);
            //int year, doy;
            //mjd2doy(beg.d, &year, &doy);
            //printf("%5d%4d\n", year, doy);
            return true;
        }
    }
    return false;
}

void RinexAtx::read_atx(FILE *fp, atx_t &atx)const
{
    int nzen = static_cast<int>((atx.zen2-atx.zen1)/atx.dzen) + 1;
    int nazi = atx.dazi==0.0 ? 1 : static_cast<int>(360.0/atx.dazi) + 2;  // NOAZI

    atx.freqs.reserve(atx.nfreq);
    atx.pcos.reserve(atx.nfreq);
    atx.pcvs.reserve(atx.nfreq);

    std::array<double, 3> pco;
    std::vector<double> row;
    std::vector<std::vector<double>> table;
    row.reserve(nzen);
    table.reserve(nazi);

    char buf[BUFSIZ];
    while (fgets(buf, sizeof(buf), fp))
    {
        if (strncmp(buf+60, "START OF FREQUENCY", 18) == 0) {
            atx.freqs.emplace_back(buf+3, 3);
            table.clear();
        } else if (strncmp(buf+60, "NORTH / EAST / UP", 17) == 0) {
            // read PCO
            sscanf(buf, "%lf%lf%lf", &pco[0], &pco[1], &pco[2]);
            pco[0]/=1000; pco[1]/=1000; pco[2]/=1000; // unit: mm => m
            atx.pcos.push_back(pco);
            // read PCV
            double pcv;
            for (int i=0; i<nazi; ++i) {
                row.clear();
                fgets(buf, sizeof(buf), atxFile_);
                for (int j=0; j<nzen; ++j) {
                    pcv = atof(buf+8+8*j)/1000.0;
                    row.push_back(pcv);
                }
                table.push_back(row);
            }
            atx.pcvs.push_back(table);
        }
        //else if (strncmp(buf+60, "END OF FREQUENCY", 16) == 0) {
        //    atx.pcvs.push_back(table);
        //}
        else if (strncmp(buf+60, "END OF ANTENNA", 14) == 0) {
            break;
        }
    }
}

void RinexAtx::skip_atx(FILE *fp)const
{
    char buf[BUFSIZ];
    while (fgets(buf, sizeof(buf), fp) &&
           strncmp(buf+60, "END OF ANTENNA", 14) != 0);
}
