#include <algorithm>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <unordered_map>
#include <utility>
#include <vector>
#include <map>
#include <stdexcept>
#include <string>

// Computes the 13-prime seed histogram for the C45 Romanoff certificate from
// first principles. The computation does not hard-code the histogram.
//
// For a base modulus M0 = product(q), q in BASE_PRIMES, and T0 = lcm ord_q(2),
// each residue class a mod M0 determines
//   nu(a) = #{r mod T0 : gcd(a - 2^r, M0) = 1}.
// This program computes the histogram #{a mod M0 : nu(a)=nu} exactly.
//
// The dynamic program processes primes one at a time. At each stage it tracks
// counts of still-admissible exponent classes modulo the lcm of the remaining
// orders. It quotients states by cyclic translation of the exponent variable,
// which is valid because all transition choices are translation-equivariant.

struct Key {
    std::vector<uint16_t> v;
};

struct Hash {
    std::size_t operator()(Key const& k) const noexcept {
        uint64_t h = 1469598103934665603ull;
        for (uint16_t x : k.v) {
            h ^= x;
            h *= 1099511628211ull;
        }
        return static_cast<std::size_t>(h);
    }
};

struct Eq {
    bool operator()(Key const& a, Key const& b) const noexcept {
        return a.v == b.v;
    }
};

using Count = unsigned long long;
using StateMap = std::unordered_map<Key, Count, Hash, Eq>;

struct HistogramResult {
    std::map<int, Count> histogram;
    Count modulus = 1;
    int period = 1;
};

static int gcd_int(int a, int b) {
    while (b) {
        int t = a % b;
        a = b;
        b = t;
    }
    return a;
}

static int lcm_int(int a, int b) {
    return a / gcd_int(a, b) * b;
}

static int ord_mod_2(int q) {
    int x = 2 % q;
    int d = 1;
    while (x != 1) {
        x = (x * 2) % q;
        ++d;
    }
    return d;
}

static Key canonical_rotation(std::vector<uint16_t> v) {
    const int L = static_cast<int>(v.size());
    if (L <= 1) return Key{std::move(v)};

    int best = 0;
    for (int shift = 1; shift < L; ++shift) {
        bool less = false;
        for (int i = 0; i < L; ++i) {
            const uint16_t a = v[(i + shift) % L];
            const uint16_t b = v[(i + best) % L];
            if (a < b) {
                less = true;
                break;
            }
            if (a > b) break;
        }
        if (less) best = shift;
    }

    if (best == 0) return Key{std::move(v)};
    std::vector<uint16_t> w(L);
    for (int i = 0; i < L; ++i) w[i] = v[(i + best) % L];
    return Key{std::move(w)};
}

static std::vector<int> remaining_lcms(std::vector<std::pair<int, int>> const& order) {
    std::vector<int> Ls(order.size() + 1);
    for (std::size_t i = 0; i <= order.size(); ++i) {
        int L = 1;
        for (std::size_t j = i; j < order.size(); ++j) L = lcm_int(L, order[j].second);
        Ls[i] = L;
    }
    return Ls;
}

static std::vector<std::pair<int, int>> base_order() {
    // Processing order chosen to keep the exact state space small. The second
    // component is checked against ord_q(2), not trusted.
    return {
        {17, 8}, {29, 28}, {61, 60}, {37, 36}, {19, 18}, {73, 9},
        {3, 2}, {7, 3}, {5, 4}, {31, 5}, {11, 10}, {13, 12}, {41, 20},
    };
}

static StateMap step(StateMap const& states, int q, int d, int L, int next_L) {
    const int non_power_weight = q - d;
    StateMap next;
    next.reserve(states.size() * static_cast<std::size_t>(std::min(d + 1, 16)));

    for (auto const& kv : states) {
        std::vector<uint16_t> const& vec = kv.first.v;
        const Count weight = kv.second;

        // Residues modulo q that are not powers of 2 remove no exponent class.
        if (non_power_weight) {
            std::vector<uint16_t> out;
            if (next_L == L) {
                out = vec;
            } else {
                out.assign(next_L, 0);
                for (int x = 0; x < L; ++x) out[x % next_L] += vec[x];
            }
            next[canonical_rotation(std::move(out))] += weight * static_cast<Count>(non_power_weight);
        }

        // A power residue 2^s mod q removes the class r == s mod d.
        for (int s = 0; s < d; ++s) {
            std::vector<uint16_t> out;
            if (next_L == L) {
                out = vec;
                for (int x = s; x < L; x += d) out[x] = 0;
            } else if (next_L == 1) {
                unsigned int total = 0;
                for (int x = 0; x < L; ++x) {
                    if (x % d != s) total += vec[x];
                }
                out.assign(1, static_cast<uint16_t>(total));
            } else {
                out.assign(next_L, 0);
                for (int x = 0; x < L; ++x) {
                    if (x % d != s) out[x % next_L] += vec[x];
                }
            }
            next[canonical_rotation(std::move(out))] += weight;
        }
    }
    return next;
}

static HistogramResult compute_histogram(
    std::vector<std::pair<int, int>> const& order,
    bool trace
) {
    Count modulus = 1;
    int period = 1;
    for (auto [q, d] : order) {
        const int actual = ord_mod_2(q);
        if (actual != d) {
            throw std::runtime_error(
                "order mismatch for q=" + std::to_string(q)
                + ": got " + std::to_string(actual)
                + ", expected " + std::to_string(d)
            );
        }
        modulus *= static_cast<Count>(q);
        period = lcm_int(period, d);
    }

    const std::vector<int> Ls = remaining_lcms(order);
    StateMap states;
    states.reserve(1024);
    std::vector<uint16_t> initial(Ls[0], static_cast<uint16_t>(period / Ls[0]));
    states[canonical_rotation(std::move(initial))] = 1;

    if (trace) {
        std::cerr << "13-prime seed histogram generator\n";
        std::cerr << "M0 = " << modulus << "\n";
        std::cerr << "T0 = " << period << "\n";
    }

    for (std::size_t i = 0; i < order.size(); ++i) {
        const int q = order[i].first;
        const int d = order[i].second;
        states = step(states, q, d, Ls[i], Ls[i + 1]);
        if (trace) {
            std::cerr << "step " << (i + 1)
                      << ": q=" << q
                      << ", ord=" << d
                      << ", L=" << Ls[i] << "->" << Ls[i + 1]
                      << ", states=" << states.size() << "\n";
        }
    }

    std::map<int, Count> histogram;
    Count total = 0;
    for (auto const& kv : states) {
        if (kv.first.v.size() != 1) {
            throw std::runtime_error("internal error: final state length is not 1");
        }
        histogram[kv.first.v[0]] += kv.second;
        total += kv.second;
    }

    if (total != modulus) {
        throw std::runtime_error(
            "histogram total mismatch: " + std::to_string(total)
            + " != " + std::to_string(modulus)
        );
    }

    if (trace) {
        std::cerr << "histogram entries = " << histogram.size() << "\n";
        std::cerr << "histogram total = " << total << "\n";
    }

    return HistogramResult{histogram, modulus, period};
}

static std::map<int, Count> direct_histogram(std::vector<std::pair<int, int>> const& order) {
    Count modulus = 1;
    int period = 1;
    std::vector<std::vector<int>> powers;
    powers.reserve(order.size());

    for (auto [q, d] : order) {
        modulus *= static_cast<Count>(q);
        period = lcm_int(period, d);
        std::vector<int> residues(d);
        int x = 1 % q;
        for (int r = 0; r < d; ++r) {
            residues[r] = x;
            x = (x * 2) % q;
        }
        powers.push_back(std::move(residues));
    }

    std::map<int, Count> histogram;
    for (Count a = 0; a < modulus; ++a) {
        int nu = 0;
        for (int r = 0; r < period; ++r) {
            bool admissible = true;
            for (std::size_t i = 0; i < order.size(); ++i) {
                const int q = order[i].first;
                const int d = order[i].second;
                if (static_cast<int>(a % q) == powers[i][r % d]) {
                    admissible = false;
                    break;
                }
            }
            if (admissible) ++nu;
        }
        histogram[nu] += 1;
    }
    return histogram;
}

static void run_self_test() {
    const std::vector<std::pair<int, int>> tiny_order = {
        {3, 2}, {5, 4}, {7, 3}, {13, 12},
    };
    const HistogramResult dp = compute_histogram(tiny_order, false);
    const std::map<int, Count> direct = direct_histogram(tiny_order);
    if (dp.histogram != direct) {
        throw std::runtime_error("self-test failed: DP histogram differs from direct enumeration");
    }
    if (dp.modulus != 1365 || dp.period != 12) {
        throw std::runtime_error("self-test failed: unexpected tiny modulus or period");
    }
    std::cerr << "seed histogram self-test passed"
              << " (modulus=" << dp.modulus
              << ", period=" << dp.period
              << ", entries=" << dp.histogram.size() << ")\n";
}

int main(int argc, char** argv) {
    bool self_test_only = false;
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--self-test") {
            self_test_only = true;
        } else if (arg == "--help" || arg == "-h") {
            std::cerr << "usage: " << argv[0] << " [--self-test]\n";
            return 0;
        } else {
            std::cerr << "unknown argument: " << arg << "\n";
            return 2;
        }
    }

    try {
        if (self_test_only) {
            run_self_test();
            return 0;
        }

        const HistogramResult result = compute_histogram(base_order(), true);

        // Data format for the Python verifier: one "nu count" pair per line.
        for (auto const& [nu, count] : result.histogram) {
            std::cout << nu << " " << count << "\n";
        }
    } catch (std::exception const& exc) {
        std::cerr << exc.what() << "\n";
        return 1;
    }
    return 0;
}
