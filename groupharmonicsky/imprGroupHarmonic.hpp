/*
 * GroupHarmonicCloseness.hpp
 *
 * Created on: 15.12.2020
 *     Author: Eugenio Angriman <angrimae@hu-berlin.de>
 */

// networkit-format

#ifndef IMPRGROUPHARMONIC_HPP
#define IMPRGROUPHARMONIC_HPP

#include <memory>
#include <vector>

#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>

using namespace std;
/**
 * @ingroup centrality
 */
class imprGroupHarmonic {

public:
    /**
     * Approximation algorithm for the group-harmonic maximization problem. The computed solutions
     * have a guaranteed $\\lambda(1 - \\frac{1}{2e})$ (directed graphs) and
     * $\\lambda(1 - \\frac{1}/{e})/2$ (undirected graphs) approximation ratio,
     * where $\\lambda$ is the ratio
     * between the minimal and the maximal edge weight. The algorithm is the one proposed in
     * Angriman et al., ALENEX 2021. The worst-case running time of this approach is quadratic, but
     * usually much faster in practice.
     *
     * @param G The input graph.
     * @param k Size of the group of nodes.
     */
    imprGroupHarmonic(const NetworKit::Graph &G, NetworKit::count k = 1);

    ~imprGroupHarmonic()=default;

    /**
     * Runs the algorithm.
     */
    void run(bool usesky, vector<NetworKit::node> &skyvec);

    /**
     * Returns the computed group.
     *
     * @return The computed group.
     */
    const vector<NetworKit::node> &groupMaxHarmonicCloseness() const;

    /**
     * Computes the group-harmonic score of the group of nodes in the given range.
     *
     * @param graph The input Graph.
     * @param first,last The range that contains the vertices in the group.
     *
     * @returns The score of the group of nodes in the given range.
     */
    template <class InputIt>
    static double scoreOfGroup(const NetworKit::Graph &graph, InputIt first, InputIt last);

    class GroupHarmonicClosenessInterface{
    public:
        std::vector<NetworKit::node> group;
        virtual void inrun(bool usesky, vector<NetworKit::node> &skyvec) = 0;
    };


private:
    const bool weighted;
    // Always one between GroupHarmonicClosenessUnweighted and GroupHarmonicClosenessWeighted, see
    // implementation file.
    unique_ptr<GroupHarmonicClosenessInterface> impl;

    static double scoreOfGroup(const NetworKit::Graph &graph, const vector<NetworKit::node> &inputGroup);

};

template <class InputIt>
double imprGroupHarmonic::scoreOfGroup(const NetworKit::Graph &graph, InputIt first, InputIt last) {
    return scoreOfGroup(graph, {first, last});
}

#endif // NETWORKIT_CENTRALITY_GROUP_HARMONIC_CLOSENESS_HPP_