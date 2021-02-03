#include <stdio.h>
#include <string.h>
#include <map>
#include <vector>
#include <string>

class Satellite
{
public:
    Satellite() {}
    explicit Satellite(const std::string &_prn) : prn(_prn) {}

    std::string prn;
    std::string svn;
    double f[2];
    double waveLen[2];
};

const double LightSpeed = 299792458;

// GPS frequency
const double GPS_f1 = 1.57542E+09;
const double GPS_f2 = 1.22760E+09;
const double GPS_f5 = 1.17645E+09;
// GLS frequency
const double GLS_f1  = 1.6020E+09;
const double GLS_f2  = 1.2460E+09;
const double GLS_df1 = 0.5625E+06;
const double GLS_df2 = 0.4375E+06;
// BDS frequency
const double BDS_f1 = 1.561098E+09;
const double BDS_f2 = 1.207140E+09;
const double BDS_f3 = 1.268520E+09;
// GAL frequency
const double GAL_f1 = 1.575420E+09;
const double GAL_f5 = 1.176450E+09;
const double GAL_f6 = 1.278750E+09;
const double GAL_f7 = 1.207140E+09;
const double GAL_f8 = 1.191795E+09;

bool sat_freq(const std::string &prn, double f[])
{
    switch (prn[0]) {
        case 'G':
        case 'J':
            f[0] = GPS_f1;
            f[1] = GPS_f2;
            break;
        case 'E':
            f[0] = GAL_f1;
            f[1] = GAL_f5; // ?
            break;
        case 'C':
            f[0] = BDS_f1;
            f[1] = BDS_f3;
            break;
        case 'R':
            f[0] = GLS_f1;
            f[1] = GLS_f2;
            break;
        default:
            fprintf(stderr, "sat_freq: unsupported system %c\n", prn[0]);
            return false;
    }
    return true;
}

bool init_sats(std::vector<Satellite> &sats)
{
    for (auto it=sats.begin(); it!=sats.end(); ++it) {
        if (!sat_freq(it->prn, it->f))
            return false;
        it->waveLen[0] = LightSpeed/it->f[0];
        it->waveLen[1] = LightSpeed/it->f[1];

        it->svn = "SVNX";
        // atx = rnxatx.atx(t, sat.prn);
        // if (nullptr != atx) {
        //     sat.svn = atx->SVN();
        // } else {
        //     fprintf(stderr, "init_sats: no SVN for %3s\n", sat.prn.c_str());
        // }
    }
    return true;
}

static void convert2raw_bias(double wl, double nl, const double f[], double raw[])
{
    // double nl  = 0.0;
    double dcb = 0.0; // ns
    double g = f[0]/f[1];
    double dw = 1/(f[0] - f[1]);
    double dn = 1/(f[0] + f[1]);
    double alpha = g*g/(g*g - 1);
    double beta = 1 - alpha;

    raw[0] = beta*dcb; // C1W
    raw[1] = -alpha*dcb; // C2W
    raw[2] = (-wl/g*dw + (1+1/g)*nl*dn)*1E9 - beta*dcb; // L1W
    raw[3] = (-wl*g*dw + (1 + g)*nl*dn)*1E9 + alpha*dcb ; // L2W
}

bool write_bias(const std::string &path, int y, int doy,
                const std::vector<Satellite> &sats,
                const std::vector<double> &wl_bias,
                const std::vector<double> &nl_bias)
{
    FILE *fp = fopen(path.c_str(), "w");
    if (fp == nullptr) {
        fprintf(stderr, "failed to create file: %s\n", path.c_str());
        return false;
    }

    // header
    fprintf(fp, "+BIAS/SOLUTION\n");
    fprintf(fp, "*BIAS SVN_ PRN STATION__ OBS1 OBS2 BIAS_START____ BIAS_END______ UNIT __ESTIMATED_VALUE____ _STD_DEV___ __ESTIMATED_SLOPE____ _STD_DEV___\n");

    const size_t nsat = sats.size();
    // int y, doy;
    double raw[4];
    const char *types[4] = { "C1W", "C2W", "L1C", "L2W" };
    const char *etypes[4] = { "C1X", "C5X", "L1X", "L5X" };
    const char *ctypes[4] = { "C2I", "C6I", "L2I", "L6I" };
    // mjd2doy(t.d, &y, &doy);
    for (size_t i=0; i!=nsat; ++i) {
        convert2raw_bias(wl_bias[i], nl_bias[i], sats[i].f, raw);
        for (int j=0; j!=4; ++j)
            if (sats[i].prn[0] == 'G')
            fprintf(fp, " OSB  %4s %3s %9s %-4s %-4s %4d:%03d:%05d %4d:%03d:%05d %-4s %21.3f %11.3f\n",
                    sats[i].svn.c_str(), sats[i].prn.c_str(), "", types[j], "", y, doy, 0, y, doy+1, 0, "ns", raw[j], 0.0);
            else if (sats[i].prn[0] == 'E')
            fprintf(fp, " OSB  %4s %3s %9s %-4s %-4s %4d:%03d:%05d %4d:%03d:%05d %-4s %21.3f %11.3f\n",
                    sats[i].svn.c_str(), sats[i].prn.c_str(), "", etypes[j], "", y, doy, 0, y, doy+1, 0, "ns", raw[j], 0.0);
            else if (sats[i].prn[0] == 'C')
            fprintf(fp, " OSB  %4s %3s %9s %-4s %-4s %4d:%03d:%05d %4d:%03d:%05d %-4s %21.3f %11.3f\n",
                    sats[i].svn.c_str(), sats[i].prn.c_str(), "", ctypes[j], "", y, doy, 0, y, doy+1, 0, "ns", raw[j], 0.0);
    }

    fprintf(fp, "-BIAS/SOLUTION\n");
    fprintf(fp, "%%=ENDBIA\n");

    fclose(fp);
    return true;
}

int main(int argc, char *argv[])
{
    if (argc != 5) {
        fprintf(stderr, "usage: zqy2snx upd_wl upd_nl year doy\n");
        return 1;
    }

    FILE *fp = fopen(argv[1], "r");
    if (fp == nullptr) {
        fprintf(stderr, "no such file: %s\n", argv[1]);
        return 1;
    }
    const int year = atoi(argv[3]);
    const int doy  = atoi(argv[4]);

    char buf[BUFSIZ];
    int index = 0;
    double val;
    std::string prn;
    std::vector<Satellite> sats;
    std::vector<double> wl_bias;
    std::map<std::string, int> map;
    while (fgets(buf, sizeof(buf), fp))
    {
        if (buf[0]==' ' && buf[4]==' ') {
            prn.assign(buf+1, 3);
            sats.emplace_back(prn);
            map[prn] = index++;

            val = atof(buf+5);
            wl_bias.push_back(val);
            // fprintf(stdout, "WL  %3s %8.3f\n", prn.c_str(), wl_bias.back());
        }
    }
    fclose(fp);
    init_sats(sats);

    // open upd_nl
    fp = fopen(argv[2], "r");
    if (fp == nullptr) {
        fprintf(stderr, "no such file: %s\n", argv[2]);
        return 1;
    }
    fgets(buf, sizeof(buf), fp);

    // Output
    std::string path = "bias";
    FILE *fout = fopen(path.c_str(), "w");
    if (fout == nullptr) {
        fprintf(stderr, "failed to create file: %s\n", path.c_str());
        return false;
    }
    fprintf(fout, "+BIAS/SOLUTION\n");
    fprintf(fout, "*BIAS SVN_ PRN STATION__ OBS1 OBS2 BIAS_START____ BIAS_END______ UNIT __ESTIMATED_VALUE____ _STD_DEV___ __ESTIMATED_SLOPE____ _STD_DEV___\n");

    const char *types[4] = { "C1W", "C2W", "L1C", "L2W" };
    const char *etypes[4] = { "C1X", "C5X", "L1X", "L5X" };
    const char *ctypes[4] = { "C2I", "C6I", "L2I", "L6I" };

    double raw[4], nl_bias;
    int sod=0;
    while (fgets(buf, sizeof(buf), fp))
    {
        if (strncmp(buf, " EPOCH-TIME", 11) == 0) {
            // sscanf(buf+13, "%d %d %d", &h, &m, &sec);
            // sod = 3600*h + 60*m + sec;
            sod = atoi(buf+20);
            continue;
        } else /*if (buf[0] == 'P')*/ {
            if (sod == 0 || buf[0]!=' ')
                continue;
            prn.assign(buf+1, 3);
            nl_bias = -atof(buf+5);
        }

        index = map[prn];
        // fprintf(stdout, "NL  %3s %6d %8.3f %8.3f\n", prn.c_str(), sod, wl_bias[index], nl_bias);
        convert2raw_bias(wl_bias[index], nl_bias, sats[index].f, raw);
        for (int j=0; j!=4; ++j)
            if (prn[0] == 'G')
            fprintf(fout, " OSB  %4s %3s %9s %-4s %-4s %4d:%03d:%05d %4d:%03d:%05d %-4s %21.3f %11.3f\n",
                    sats[index].svn.c_str(), prn.c_str(), "",  types[j], "", year, doy, sod-900, year, doy, sod, "ns", raw[j], 0.0);
            else if (prn[0] == 'E')
            fprintf(fout, " OSB  %4s %3s %9s %-4s %-4s %4d:%03d:%05d %4d:%03d:%05d %-4s %21.3f %11.3f\n",
                    sats[index].svn.c_str(), prn.c_str(), "", etypes[j], "", year, doy, sod-900, year, doy, sod, "ns", raw[j], 0.0);
            else if (prn[0] == 'C')
            fprintf(fout, " OSB  %4s %3s %9s %-4s %-4s %4d:%03d:%05d %4d:%03d:%05d %-4s %21.3f %11.3f\n",
                    sats[index].svn.c_str(), prn.c_str(), "", ctypes[j], "", year, doy, sod-900, year, doy, sod, "ns", raw[j], 0.0);
    }

    fprintf(fout, "-BIAS/SOLUTION\n");
    fprintf(fout, "%%=ENDBIA\n");
    fclose(fp);
    fclose(fout);
    return 0;
}