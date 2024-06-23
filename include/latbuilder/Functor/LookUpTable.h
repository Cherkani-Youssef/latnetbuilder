#ifndef LATBUILDER__FUNCTOR__LOOKUPTABLE_H
#define LATBUILDER__FUNCTOR__LOOKUPTABLE_H
#include <vector>
#include <cmath>
#include <stdexcept>
/**
 * This class computes the Walsh Figure of Merit (WAFOM) for a Digital Net in
 * Base 2 using a lookup table to accelerate the computation.
 *
 * For Matsumoto's definition
 *
 * WAFOM(P) = -1 + \frac{1}{|P|} \left( \sum_{\mathbf{x}_i \in P}
 * \prod_{j=1}^{s} \prod_{l=1}^{w} \left( 1 + \eta(x_{i,j,l}) 2^{-l} \right)
 * \right)
 *
 * The same formula applies to Yoshiki's definition; simply replace \( l \) with
 * \( l + 1 \).
 *
 * WAFOM(P) = -1 + \frac{1}{|P|} \left( \sum_{\mathbf{x}_i \in P}
 * \prod_{j=1}^{s} \prod_{l=1}^{w} \left( 1 + \eta(x_{i,j,l}) 2^{-(l+1)} \right)
 * \right) .
 *
 * For Goda's definition:
 *
 * \[ \mathcal{W}(P, \mu) = \sqrt{-1 + \frac{1}{|P|} \left( \sum_{\mathbf{x}_i
 * \in P} \prod_{j=1}^{s} \prod_{l=1}^{w} \left( 1 + \eta(x_{i,j,l}) 2^{-2l}
 * \right) \right)} \] The same formula applies to Yoshiki's definition; simply
 * replace \( l \) with \( l + 1 \).
 *
 * The same formula applies to Yoshiki's definition; simply replace \( l \) with
 * \( l + 1 \).
 *
 * \[ \mathcal{W}(P, \mu + h) = \sqrt{-1 + \frac{1}{|P|} \left(
 * \sum_{\mathbf{x}_i \in P} \prod_{j=1}^{s} \prod_{l=1}^{w} \left( 1 +
 * \eta(x_{i,j,l}) 2^{-2(l+1)} \right) \right)} \]
 *
 * The method, proposed by Shin Harase in "A search for extensible low-WAFOM
 * point sets," Monte Carlo Methods and Applications 22.4 (2016), pp. 349–357,
 * involves dividing a coordinate of a point with 'w' precision into 'q'
 * subparts, each containing 'l' elements if q is a divisor of the number of
 * 'w'. When 'q' is not a divisor of the number of output digits then will
 * another be a segment of lehgth t where \(t = w \mod q\) .
 *
 * A coordinate of a point is expressed as \(X^i = d_1^i, \ldots, d_q^i,
 * d_{q+1}^i\), where for \(1 \leq c \leq q\), the length of each segment is \(l
 * = w \div q\), and the remainder \(d_{q+1}\) is of length \(t = w \mod q\).
 * Each segment \(d_c^i\) is represented as \(x_{i, (c-1)*l+1}, \ldots, x_{i,
 * c*l}\).
 *
 * Instead of computing \(\prod_{j=1}^{w} [1 + (-1)^{(-1)^{x_{i,j}}} 2^{-j}] -
 * 1\), we use \(\prod_{1 \leq c \leq q} \text{table}_c[d_c^i]\), where
 * \(\text{table}_c[d_c^i]\) are precomputed values: \(\text{table}_c[d_c^i] =
 * \prod_{1 \leq j \leq l} (1 + (-1)^{x_{i, (c-1)*l + j}} \cdot 2^{-((c-1)*l + j
 * + 1)})\). This significantly reduces computation time.
 *
 * We use Kahan summation algorithm improves numerical accuracy by compensating
 * for floating-point errors, reducing round-off error accumulation.
 **/
class LookUpTable
{

public:
    /**
     * Constructor for LookUpTableWafom.
     *
     * @param n      the number of bits
     * @param q      a positive integer that divides n into l segements
     * @param l      length of the segment
     * @param h      for Matsumoto's definition set to 0, and for Yoshiki's set to 1
     * @param factor set to 1 for wafom, set to 2 to get Wafom for RMSE
     * */
    LookUpTable(int n, int q, int l, double h, double factor)
        : n(n), q(q), l(l), h(h), factor(factor)
    {
        if (h != 0 && h != 1)
        {
            throw std::invalid_argument("h must be either 0 or 1: for Matsumoto's definition set to 0, and for Yoshiki's set to 1");
        }
        if (factor != 1 && factor != 2)
        {
            throw std::invalid_argument("factor must be set to 1 for wafom, and set to 2 for Wafom for RMSE");
        }
        if (n < q)
        {
            throw std::invalid_argument("n must be larger than or equal to q");
        }

        long tableSize = q * (1L << l) + (n % q != 0 ? (1L << (n % q)) : 0);
        lookupTable.resize(tableSize);

        generate();
    }

    double get(int c, int e) const
    {
        if (c == 0 || c > q + (n % q != 0 ? 1 : 0))
        {
            throw std::out_of_range("c is out of range");
        }
        if (e >= (1 << l) || (c == q + 1 && e >= (1 << (n % q))))
        {
            throw std::out_of_range("e is out of range");
        }
        return lookupTable[index(c, e)];
    }

private:
    int q;
    int l;
    std::vector<double> lookupTable;
    double h;
    double factor;
    int n;

    int index(int m, int e) const
    {
        if (m <= q)
        {
            return (m - 1) * (1 << l) + e;
        }
        else
        {
            return q * (1 << l) + e;
        }
    }

    void generate()
    {
        int length = 1 << l;
        for (int c = 1; c <= q; ++c)
        {
            for (int e = 0; e < length; ++e)
            {
                lookupTable[index(c, e)] = computeProduct(c, e, l);
            }
        }
        if (n % q != 0)
        {
            int remainderLength = n % q;
            for (int e = 0; e < (1 << remainderLength); ++e)
            {
                lookupTable[index(q + 1, e)] = computeProduct(q + 1, e, remainderLength);
            }
        }
    }

    double computeProduct(int c, int e, int length) const
    {
        double product = 1.0;
        for (int j = 1; j <= length; ++j)
        {
            int e_j = (e >> (length - j)) & 1;
            double two_exponent = std::pow(2.0, -factor * ((c - 1) * l + (j + h)));
            product *= (1 + (1 - 2 * e_j) * two_exponent);
        }
        return product;
    }
};
#endif
