#include "CartCoor.h"

#if 0
std::ostream &operator<<(std::ostream &out, const CartCoor c)
{
    using std::setw;
    std::streamsize strsize = out.precision();
    out << std::setprecision(4) << std::fixed;
    out << setw(15) << c.x
        << setw(15) << c.y
        << setw(15) << c.z;

    out.precision(strsize);
    return out;
}
#endif
