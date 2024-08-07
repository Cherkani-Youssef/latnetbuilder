#ifndef GET_COLS_REVERSE_CACHE_H
#define GET_COLS_REVERSE_CACHE_H

#include <unordered_map>
#include <vector>
#include <stdexcept>
#include <cmath>
#include <boost/dynamic_bitset.hpp>
#include "netbuilder/Types.h"
#include "netbuilder/DigitalNet.h"
#include "latbuilder/LFSR258.h"

namespace NetBuilder
{

    class GetColsReverseCache
    {

    public:
        GetColsReverseCache(const AbstractDigitalNet &net) : m_net(net)
        {
            int dim = net.dimension();
            columns.resize(dim); // Resize the columns vector
            for (int i = 0; i < dim; i++)
            {
                auto matrix = m_net.generatingMatrix(i);
                columns[i] = matrix.getColsReverse();
            }
        }



        /***
         * @param i: index of the point in the pointset
         * @param dim: dimension of the pointset
         * @param cachedCurPoint array to store the point
         */
        void getPoint(const int &i, const int &dim, uint64_t cachedCurPoint[])
        {
            // Iterate over each dimension
            for (int j = 0; j < dim; j++)
            {
                uint64_t res = 0;
                // const std::vector<uInteger> &getColsReverse = get(j);

                for (size_t c = 0; c <  columns[j].size(); c++)
                    res ^= ((i >> c) & 1) * columns[j][c];
                // res ^= ((i >> c) & 1) * getColsReverse[c];

                cachedCurPoint[j] = res;
            }
        }

        /***
         * @param i: index of the point in the pointset
         * @param projection: subset of the dimension of the pointset
         * @param cachedCurPoint array to store the point
         *
         */
        void getPointForProj(const int &i, const LatticeTester::Coordinates &projection, uint64_t cachedCurPoint[])
        {
            uint64_t res;
            int j = 0;
            // Iterate over each dimension
            for (auto dim : projection)
            {

                res = 0;
                // const std::vector<uInteger> &getColsReverse = get(dim);

                for (size_t c = 0; c < columns[j].size(); c++)
                {
                    // res ^= ((i >> c) & 1) * getColsReverse[c];
                     res ^= ((i >> c) & 1) * columns[dim][c];
                }
                cachedCurPoint[j] = res;
                j++;
            }
        }

        /**
         * @param uint64_t x_i integer representation of the point,
         * @param c is the  index of the segment
         * @param w is the total number of bits
         * @param q number of segments
         */
        uInteger compute_d_c(uint64_t x_i, uInteger c, uInteger w, uInteger q) const
        {
            uInteger l = w / q;                              // Length of each full segment
            uInteger remainder = w % q;                      // Length of the remainder segment, if any
            uInteger segments = q + (remainder > 0 ? 1 : 0); // Total number of segments including remainder if exists

            // Boundary check for c
            if (c < 1 || c > segments)
            {
                throw std::invalid_argument("Segment index c is out of bounds.");
            }

            if (c <= q)
            {
                uInteger startBitFromMsb = (c * l) - 1;
                uInteger shiftAmount = w - startBitFromMsb - 1;
                return (x_i >> shiftAmount) & ((1 << l) - 1);
            }
            else
            {
                // Handle the remainder segment correctly
                return x_i & ((1 << remainder) - 1); // Extract the remainder segment value
            }
        }

    private:
        std::unordered_map<int, std::vector<uInteger>> cache;
        const AbstractDigitalNet &m_net;

        std::vector<std::vector<uInteger>> columns;

        /**
         * @param GeneratingMatrix matrix
         * @return  std::vector<uInteger>  the integer representations of the a column of a matrix
         */
        std::vector<uInteger> computeGetColsReverse(const GeneratingMatrix &matrix) const
        {
            const unsigned int nCols = matrix.nCols();
            const unsigned int nRows = matrix.nRows();
            std::vector<uInteger> colsReverse(nCols, 0);

            for (unsigned int j = 0; j < nCols; ++j)
            {
                uInteger s = 0;
                for (unsigned int i = 0; i < nRows; ++i)
                {
                    s += matrix(i, j) << (nRows - i - 1);
                }
                colsReverse[j] = s;
            }
            return colsReverse;
        }
    };

} // namespace NetBuilder

#endif // GET_COLS_REVERSE_CACHE_H
