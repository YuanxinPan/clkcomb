#include "RinexAtx.h"
#include "../io/io.h"
#include "../coord/coord.h"

#include <string.h>
#include <stdlib.h>
#include <algorithm>
#include <iterator>
#include <Eigen/Dense>
#include <pppx/const.h>

bool RinexAtx::atx_t::pco(const std::string &f, double *pco)const
{
    auto it = std::find(freq.begin(), freq.end(), f);
    int i = it - freq.begin();
    if (it == freq.end())
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
    auto it = std::find(freq.begin(), freq.end(), f1);
    i = it - freq.begin();
    it = std::find(freq.begin(), freq.end(), f2);
    j = it - freq.begin();

    // fprintf(stderr, "%3s %3s:%2d %3s:%2d\n", name.c_str(), f1.c_str(), i, f2.c_str(), j);

    pco[0] = coef[0]*pcos[i][0] + coef[1]*pcos[j][0];
    pco[1] = coef[0]*pcos[i][1] + coef[1]*pcos[j][1];
    pco[2] = coef[0]*pcos[i][2] + coef[1]*pcos[j][2];
}

void RinexAtx::atx_t::pco(const double *xsat, const double *vsat,
                          const double *xsun, double *pco, double *R)const
{
    Eigen::Matrix3d _R;
    if (R != nullptr) {
        // copy rotation matrix
        memcpy(_R.data(), R, sizeof(double)*9); // dst src byte
        memcpy(R_, R, sizeof(double)*9);

        Eigen::Vector3d antpco;
        antpco[0] = IF_coef1*pcos[0][0] + IF_coef2*pcos[1][0];
        antpco[1] = IF_coef1*pcos[0][1] + IF_coef2*pcos[1][1];
        antpco[2] = IF_coef1*pcos[0][2] + IF_coef2*pcos[1][2];

        Eigen::Vector3d val = _R*antpco;
        pco[0] = val[0];
        pco[1] = val[1];
        pco[2] = val[2];
        return;
    }

    double *xscf = _R.data();
    double *yscf = xscf+3;
    double *zscf = xscf+6;

    //double vscf[3], u_orbn[3], u_e2sun[3];
    //sc-fixed z-axis, from sc to earth
    unit_vector(xsat, zscf);
    zscf[0] = -zscf[0];
    zscf[1] = -zscf[1];
    zscf[2] = -zscf[2];

    //cos of the angle between the sv radius and sun radius vectors
    //unit_vector(xsun, u_e2sun);
    //double svbcos = dot(zscf, u_e2sun);
    //Sun angle w.r.t orbital plane
    //unit_vector(vsat, vscf);
    //cross(zscf, vscf, u_orbn);
    //unit_vector(u_orbn, u_orbn);

    //unint vector from sc to sun
    double u_sc2sun[3];
    u_sc2sun[0] = xsun[0] - xsat[0];
    u_sc2sun[1] = xsun[1] - xsat[1];
    u_sc2sun[2] = xsun[2] - xsat[2];
    unit_vector(u_sc2sun, u_sc2sun);

    //the sc-fixed y-axis
    cross(zscf, u_sc2sun, yscf);
    unit_vector(yscf, yscf);
    //the sc-fixed x-axis
    cross(yscf, zscf, xscf);

    // copy rotation matrix
    memcpy(R_, _R.data(), sizeof(double)*9);

    Eigen::Vector3d antpco;
    antpco[0] = IF_coef1*pcos[0][0] + IF_coef2*pcos[1][0];
    antpco[1] = IF_coef1*pcos[0][1] + IF_coef2*pcos[1][1];
    antpco[2] = IF_coef1*pcos[0][2] + IF_coef2*pcos[1][2];

    Eigen::Vector3d val = _R*antpco;
    pco[0] = val[0];
    pco[1] = val[1];
    pco[2] = val[2];
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
        fprintf(stderr, ANSI_BOLD_RED "error: " ANSI_RESET
                "RinexAtx::read: no such file: %s\n", path.c_str());
        return false;
    }

    // check whether absolute atx file
    char buf[128];
    fgets(buf, sizeof(buf), atxFile_);
    fgets(buf, sizeof(buf), atxFile_);
    if (buf[0] != 'A') {
        fprintf(stderr, ANSI_BOLD_RED "error: " ANSI_RESET
                "RinexAtx::read: not an absolute atx: %s\n", path.c_str());
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
    fprintf(stdout, ANSI_BOLD_YELLOW "warning: " ANSI_RESET
            "RinexAtx::atx:: no %s\n", ant.c_str());
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
        if (strncmp(buf+60, "START OF ANTENNA", 16) == 0)
            beg = end = MJD();
        else if (strncmp(buf+60, "TYPE / SERIAL NO", 16) == 0) {
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
        }
        else if (strncmp(buf+60, "DAZI", 4) == 0)
            sscanf(buf, "%lf", &atx.dazi);
        else if (strncmp(buf+60, "ZEN1 / ZEN2 / DZEN", 18) == 0)
            sscanf(buf, "%lf%lf%lf", &atx.zen1, &atx.zen2, &atx.dzen);
        else if (strncmp(buf+60, "# OF FREQUENCIES", 16) == 0)
            sscanf(buf, "%d", &atx.nfreq);
        else if (strncmp(buf+60, "VALID FROM", 10) == 0) {
            calend.read(buf);
            beg = calend;
            if (beg > t)
                skip_atx(atxFile_);
        }
        else if (strncmp(buf+60, "VALID UNTIL", 11) == 0) {
            calend.read(buf);
            end = calend;
            if (end < t)
                skip_atx(atxFile_);
        }
        else if (strncmp(buf+60, "SINEX CODE", 10) == 0) {
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

    atx.freq.reserve(atx.nfreq);
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
            atx.freq.emplace_back(buf+3, 3);
            table.clear();
        }
        else if (strncmp(buf+60, "NORTH / EAST / UP", 17) == 0) {
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
    static char buf[BUFSIZ];
    while (fgets(buf, sizeof(buf), fp) &&
           strncmp(buf+60, "END OF ANTENNA", 14) != 0);
}

#if 0
bool RinexAtx::read(const std::string &path, const std::vector<std::string> &prns, MJD t)
{
    atxFile_ = fopen(path.c_str(), "r");
    if (atxFile_ == nullptr) {
        fprintf(stderr, ANSI_BOLD_RED "error: " ANSI_RESET
                "RinexAtx::read: no such file: %s\n", path.c_str());
        return false;
    }

    char buf[256];
    fgets(buf, sizeof(buf), atxFile_);
    fgets(buf, sizeof(buf), atxFile_);
    if (buf[0] != 'A') {
        fprintf(stderr, ANSI_BOLD_RED "error: " ANSI_RESET
                "RinexAtx::read: not an absolute atx: %s\n", path.c_str());
        return false;
    }
    while (strncmp(buf+60, "END OF HEADER", 13)) {
        fgets(buf, sizeof(buf), atxFile_);
    }

    return read_satatx(prns, t); // && read_rcvatx(prns);
}

bool RinexAtx::read_satatx(const std::vector<std::string> &prns, MJD t)
{
    char buf[256];
    SatAtx satAtx;
    CalendTime calend;
    MJD beg, end;
    int ifreq = 0;
    satAtxs_.reserve(prns.size());
    while (fgets(buf, sizeof(buf), atxFile_))
    {
        if (strncmp(buf+60, "START OF ANTENNA", 16) == 0) {
            ifreq = 0;
            beg = end = MJD();
        }
        else if (strncmp(buf+60, "TYPE / SERIAL NO", 16) == 0) {
            if (buf[20] == ' ') {  // receiver antenna
                backspace(atxFile_);
                break;
            }
            satAtx.name.assign(buf+20, 3);
            // check whether required
            if (!std::binary_search(prns.begin(), prns.end(), satAtx.name)) {
                while (fgets(buf, sizeof(buf), atxFile_) &&
                       strncmp(buf+60, "END OF ANTENNA", 14) != 0);
            }
        }
        else if (strncmp(buf+60, "ZEN1 / ZEN2 / DZEN", 18) == 0) {
            sscanf(buf, "%lf%lf%lf", &satAtx.zen1, &satAtx.zen2, &satAtx.dzen);
        }
        else if (strncmp(buf+60, "VALID FROM", 10) == 0) {
            calend.read(buf);
            beg = calend;
        }
        else if (strncmp(buf+60, "VALID UNTIL", 11) == 0) {
            calend.read(buf);
            end = calend;
        }
        else if (strncmp(buf+60, "START OF FREQUENCY", 18) == 0) {
            continue;
        }
        else if (strncmp(buf+60, "NORTH / EAST / UP", 17) == 0) {
            sscanf(buf, "%lf%lf%lf", satAtx.pco[ifreq], satAtx.pco[ifreq]+1, satAtx.pco[ifreq]+2);
            for (int i=0; i<3; ++i)
                satAtx.pco[ifreq][i] /= 1000; // unit: mm => m
            fgets(buf, sizeof(buf), atxFile_);
            int n = static_cast<int>((satAtx.zen2 - satAtx.zen1)/satAtx.dzen) + 1;
            if (n > MaxSatZen) {
                fprintf(stderr, ANSI_BOLD_RED "error: " ANSI_RESET
                        "RinexAtx::read: %2d > MaxSatZen(%2d)\n", n, MaxSatZen);
                return false;
            }
            for (int i=0; i<n; ++i) {
                sscanf(buf+8+8*i, "%lf", satAtx.pcv[ifreq]+i);
                satAtx.pcv[ifreq][i] /= 1000; // unit: mm => m
            }
        }
        else if (strncmp(buf+60, "END OF FREQUENCY", 16) == 0) {
            ++ifreq;
        }
        else if (strncmp(buf+60, "END OF ANTENNA", 14) == 0) {
            if (ifreq > MaxSatFreq) {
                fprintf(stderr, ANSI_BOLD_RED "error: " ANSI_RESET
                        "RinexAtx::read: %2d > MaxSatFreq(%2d)\n", ifreq, MaxSatFreq);
                return false;
            }
            satAtx.nfreq = ifreq;
            if ((t>=beg && t<=end) || (t>=beg && end.d==0))
                satAtxs_.push_back(satAtx);
        }
    }

    std::sort(satAtxs_.begin(), satAtxs_.end());
    //for (auto atx : satAtxs_) {
    //    printf("%3s %8.3f %8.3f %8.3f\n", atx.name.c_str(), 1000*atx.pco[0][0], 1000*atx.pco[0][1], 1000*atx.pco[0][2]);
    //    int n = static_cast<int>((atx.zen2-atx.zen1)/atx.dzen) + 1;
    //    for (int i=0; i<n; ++i)
    //        printf("%6.1f ", 1000*atx.pcv[0][i]);
    //    printf("\n\n");
    //}
    return check(prns);
}

bool RinexAtx::read_rcvatx(const std::vector<std::string> &rcvs)
{
    char buf[256];
    RcvAtx rcvAtx;
    int ifreq = 0;
    rcvAtxs_.reserve(rcvs.size());
    while (fgets(buf, sizeof(buf), atxFile_))
    {
        if (strncmp(buf+60, "START OF ANTENNA", 16) == 0) {
            ifreq = 0;
        }
        else if (strncmp(buf+60, "TYPE / SERIAL NO", 16) == 0) {
            rcvAtx.name.assign(buf, 20);
            //if (rcvAtx.name[16] == ' ')
            //    rcvAtx.name.replace(16, 4, "NONE");
            // check whether required
            if (!std::binary_search(rcvs.begin(), rcvs.end(), rcvAtx.name)) {
                while (fgets(buf, sizeof(buf), atxFile_) &&
                       strncmp(buf+60, "END OF ANTENNA", 14) != 0);
            }
        }
        else if (strncmp(buf+60, "ZEN1 / ZEN2 / DZEN", 18) == 0) {
            sscanf(buf, "%lf%lf%lf", &rcvAtx.zen1, &rcvAtx.zen2, &rcvAtx.dzen);
        }
        else if (strncmp(buf+60, "START OF FREQUENCY", 18) == 0) {
            continue;
        }
        else if (strncmp(buf+60, "NORTH / EAST / UP", 17) == 0) {
            sscanf(buf, "%lf%lf%lf", rcvAtx.pco[ifreq], rcvAtx.pco[ifreq]+1, rcvAtx.pco[ifreq]+2);
            for (int i=0; i<3; ++i)
                rcvAtx.pco[ifreq][i] /= 1000; // unit: mm => m
            fgets(buf, sizeof(buf), atxFile_);
            int n = static_cast<int>((rcvAtx.zen2 - rcvAtx.zen1)/rcvAtx.dzen) + 1;
            if (n > MaxRcvZen) {
                fprintf(stderr, ANSI_BOLD_RED "error: " ANSI_RESET
                        "RinexAtx::read: %2d > MaxRcvZen(%2d)\n", n, MaxRcvZen);
                return false;
            }
            for (int i=0; i<n; ++i) {
                //sscanf(buf+8+8*i, "%lf", rcvAtx.pcv[ifreq]+i);
                //rcvAtx.pcv[ifreq][i] /= 1000; // unit: mm => m
            }
        }
        else if (strncmp(buf+60, "END OF FREQUENCY", 16) == 0) {
            ++ifreq;
        }
        else if (strncmp(buf+60, "END OF ANTENNA", 14) == 0) {
            if (ifreq > MaxSatFreq) {
                fprintf(stderr, ANSI_BOLD_RED "error: " ANSI_RESET
                        "RinexAtx::read: %2d > MaxSatFreq(%2d)\n", ifreq, MaxSatFreq);
                return false;
            }
            rcvAtx.nfreq = ifreq;
            rcvAtxs_.push_back(rcvAtx);
        }
    }
    std::sort(rcvAtxs_.begin(), rcvAtxs_.end());
    return true;
}

bool RinexAtx::check(const std::vector<std::string> &prns)
{
    if (satAtxs_.size() != prns.size()) {
        std::vector<std::string> prns_;
        std::vector<std::string> diff;
        prns_.reserve(prns.size());
        diff.reserve(prns.size()/4);
        for (auto it=satAtxs_.begin(); it!=satAtxs_.end(); ++it) {
            prns_.push_back(it->name);
        }
        std::set_difference(prns.begin(), prns.end(), prns_.begin(), prns_.end(), std::back_inserter(diff));

        fprintf(stderr, ANSI_BOLD_RED "error: " ANSI_RESET
                "RinexAtx::read: no atx for:");
        for (auto it=diff.begin(); it!=diff.end(); ++it)
            fprintf(stderr, " %3s", it->c_str());
        fprintf(stderr, "\n");
        return false;
    }
    return true;
}
#endif
